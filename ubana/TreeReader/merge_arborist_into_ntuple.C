// merge_arborist_into_ntuple.C
//
// Usage from shell:
//   root -l -b -q 'merge_arborist_into_ntuple.C("arborist_out.root","source_ntuple.root","out.root")'
//
// What it does:
// - Copies source_ntuple.root -> out.root (preserves all directories/objects)
// - Adds a new TTree to out.root:
//     TTree "spline_weights"
//   containing Arborist's mcweight content in a simple/ROOT-stable format:
//     - run, subrun, event, entry
//     - one branch per mcweight key, each a vector<double> (per event)
//
// Notes:
// - This does NOT modify existing trees/directories in the source file; it only
//   appends/overwrites the "spline_weights" tree in the output file.
// - Events are matched by Arborist's "event" number, which equals (source entry + 1).
//   (Art assigns sequential event IDs starting at 1 when reading from TreeReader.)
// - If an event exists in source but not in Arborist output (e.g., skipped due to
//   processing errors), the weights are filled with sentinel value -9898.

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
  const char* kDefaultSourceTreePath = "singlephotonana/eventweight_tree";
  const char* kDefaultArboristTreePath = "arborist/eventweight_tree";
  const char* kArboristBranchName = "mcweight";
  const char* kOutputTreeName = "spline_weights";
  const double kMissingSentinel = -9898.0;

  std::string SanitizeBranchName(const std::string& in) {
    // ROOT branch names can be finicky with certain characters. We'll keep
    // alnum + '_' and replace everything else with '_'. Also avoid leading digits.
    std::string out;
    out.reserve(in.size());
    for (char c : in) {
      const bool ok =
          (c >= 'a' && c <= 'z') ||
          (c >= 'A' && c <= 'Z') ||
          (c >= '0' && c <= '9') ||
          (c == '_');
      out.push_back(ok ? c : '_');
    }
    if (out.empty()) out = "w";
    if (out[0] >= '0' && out[0] <= '9') out = std::string("w_") + out;
    return out;
  }

}

int merge_arborist_into_ntuple(const char* arborist_file,
                              const char* source_file,
                              const char* output_file,
                              const char* source_tree_path = kDefaultSourceTreePath,
                              const char* arborist_tree_path = kDefaultArboristTreePath) {
  if (!arborist_file || !source_file || !output_file) {
    std::cerr << "ERROR: missing args\n";
    std::cerr << "Usage:\n"
              << "  root -l -b -q 'merge_arborist_into_ntuple.C(\"arb.root\",\"src.root\",\"out.root\")'\n";
    return 1;
  }

  const std::string src(source_file);
  const std::string out(output_file);
  if (src == out) {
    std::cerr << "ERROR: output_file must be different from source_file (would overwrite in-place)\n";
    return 2;
  }

  // Copy the source file to preserve *all* directories/objects.
  std::cout << "Copying '" << src << "' -> '" << out << "' ...\n";
  const int cp_rc = gSystem->CopyFile(src.c_str(), out.c_str(), /*overwrite*/ kTRUE);
  if (cp_rc != 0) {
    std::cerr << "ERROR: CopyFile failed with code " << cp_rc << "\n";
    return 3;
  }

  TFile fa(arborist_file, "READ");
  if (fa.IsZombie()) {
    std::cerr << "ERROR: failed to open arborist_file: " << arborist_file << "\n";
    return 4;
  }

  TFile fout(output_file, "UPDATE");
  if (fout.IsZombie()) {
    std::cerr << "ERROR: failed to open output_file for UPDATE: " << output_file << "\n";
    return 5;
  }

  auto* ta = dynamic_cast<TTree*>(fa.Get(arborist_tree_path));
  if (!ta) {
    std::cerr << "ERROR: could not find arborist tree at '" << arborist_tree_path << "' in " << arborist_file << "\n";
    fa.ls();
    return 6;
  }

  auto* tr = dynamic_cast<TTree*>(fout.Get(source_tree_path));
  if (!tr) {
    std::cerr << "ERROR: could not find source tree at '" << source_tree_path << "' in " << output_file << "\n";
    fout.ls();
    return 7;
  }

  const Long64_t na = ta->GetEntries();
  const Long64_t nr = tr->GetEntries();
  std::cout << "Arborist entries: " << na << "\n";
  std::cout << "Source   entries: " << nr << "\n";
  
  // Allow different entry counts - we'll match by (run, subrun, event)
  if (na != nr) {
    std::cout << "INFO: Entry counts differ. Will match events by (run, subrun, event).\n";
    std::cout << "      Missing events will be filled with sentinel value " << kMissingSentinel << "\n";
  }

  // Hook up the arborist branch we want to copy.
  std::map<std::string, std::vector<double>>* mcweight = nullptr;
  std::cout << "Setting arborist branch address '" << kArboristBranchName << "'...\n";
  if (ta->SetBranchAddress(kArboristBranchName, &mcweight) < 0) {
    std::cerr << "ERROR: could not SetBranchAddress('" << kArboristBranchName << "') on arborist tree\n";
    ta->Print();
    return 9;
  }
  std::cout << "OK\n";

  // Hook up event branch on arborist tree for building lookup map
  // Art assigns event IDs starting at 1, so arborist_event = source_entry + 1
  Int_t arb_event = 0;
  if (ta->SetBranchAddress("event", &arb_event) < 0) {
    std::cerr << "ERROR: arborist tree is missing branch 'event'\n";
    ta->Print();
    return 15;
  }

  // Build lookup map: arborist_event -> arborist tree index
  // arborist_event = source_entry + 1 (art's sequential event numbering)
  std::cout << "Building event lookup map from Arborist output...\n";
  std::map<Int_t, Long64_t> arboristEventMap;
  for (Long64_t i = 0; i < na; ++i) {
    ta->GetEntry(i);
    arboristEventMap[arb_event] = i;
  }
  std::cout << "Indexed " << arboristEventMap.size() << " events from Arborist.\n";

  // Peek at the first arborist entry to discover which keys exist.
  // We'll create one output branch per key.
  if (na > 0) {
    ta->GetEntry(0);
  }
  if (!mcweight) {
    std::cerr << "ERROR: Arborist mcweight pointer is null after reading entry 0.\n";
    return 13;
  }
  if (mcweight->empty()) {
    std::cerr << "ERROR: Arborist mcweight map is empty in entry 0; nothing to write.\n";
    return 14;
  }

  // Hook up source identifiers (used only for convenience in the output tree).
  Int_t run = 0;
  Int_t subrun = 0;
  Int_t event = 0;
  if (tr->SetBranchAddress("run", &run) < 0) {
    std::cerr << "ERROR: source tree is missing branch 'run'\n";
    tr->Print();
    return 10;
  }
  if (tr->SetBranchAddress("subrun", &subrun) < 0) {
    std::cerr << "ERROR: source tree is missing branch 'subrun'\n";
    tr->Print();
    return 11;
  }
  if (tr->SetBranchAddress("event", &event) < 0) {
    std::cerr << "ERROR: source tree is missing branch 'event'\n";
    tr->Print();
    return 12;
  }

  // (Re)create the output tree at the ROOT file top-level.
  fout.cd();
  fout.Delete(std::string(std::string(kOutputTreeName) + ";*").c_str());

  TTree ts(kOutputTreeName, "Spline weights from Arborist mcweight (matched by entry number)");

  Long64_t entry = 0;

  ts.Branch("run", &run, "run/I");
  ts.Branch("subrun", &subrun, "subrun/I");
  ts.Branch("event", &event, "event/I");
  ts.Branch("entry", &entry, "entry/L");

  struct WeightBranch {
    std::string key;        // original map key
    std::string branchName; // sanitized/unique branch name
    std::vector<double>* buf = nullptr;
  };
  std::vector<WeightBranch> branches;
  branches.reserve(mcweight->size());

  std::set<std::string> usedBranchNames;
  for (const auto& kv : *mcweight) {
    const std::string& key = kv.first;
    std::string bname = SanitizeBranchName(key);

    // Ensure uniqueness after sanitization.
    if (usedBranchNames.count(bname)) {
      int idx = 2;
      std::string tryname;
      do {
        tryname = bname + "_" + std::to_string(idx++);
      } while (usedBranchNames.count(tryname));
      bname = tryname;
    }
    usedBranchNames.insert(bname);

    auto* v = new std::vector<double>();
    branches.push_back(WeightBranch{key, bname, v});
    ts.Branch(bname.c_str(), v);
  }

  // Determine how many weight values each branch should have (from first valid entry)
  std::map<std::string, size_t> expectedSizes;
  ta->GetEntry(0);
  for (const auto& kv : *mcweight) {
    expectedSizes[kv.first] = kv.second.size();
  }

  // Fill output tree entry-by-entry, matching by arborist event number
  // Art assigns event = entry + 1 when reading from TreeReader
  Long64_t missingCount = 0;
  for (Long64_t i = 0; i < nr; ++i) {
    entry = i;
    tr->GetEntry(i); // loads run/subrun/event from source

    // Reset buffers
    for (auto& br : branches) br.buf->clear();

    // Look up this entry in the arborist map (arborist event = source entry + 1)
    Int_t expectedArbEvent = static_cast<Int_t>(i + 1);
    auto it = arboristEventMap.find(expectedArbEvent);

    if (it != arboristEventMap.end()) {
      // Found in arborist - load the weights
      ta->GetEntry(it->second);
      
      if (mcweight) {
        for (auto& br : branches) {
          auto wit = mcweight->find(br.key);
          if (wit != mcweight->end()) {
            *(br.buf) = wit->second;
          } else {
            br.buf->clear();
          }
        }
      }
    } else {
      // Not found in arborist - fill with sentinel values
      missingCount++;
      for (auto& br : branches) {
        size_t sz = expectedSizes.count(br.key) ? expectedSizes[br.key] : 1;
        br.buf->assign(sz, kMissingSentinel);
      }
    }

    ts.Fill();
  }

  fout.cd();
  ts.Write(kOutputTreeName, TObject::kOverwrite);

  // Cleanup heap buffers
  for (auto& br : branches) {
    delete br.buf;
    br.buf = nullptr;
  }

  fout.Close();
  fa.Close();

  std::cout << "Wrote merged output: " << output_file << "\n";
  std::cout << "Added/updated TTree: " << kOutputTreeName << "\n";
  std::cout << "Events matched: " << (nr - missingCount) << " / " << nr << "\n";
  if (missingCount > 0) {
    std::cout << "Events with sentinel weights (" << kMissingSentinel << "): " << missingCount << "\n";
  }
  return 0;
}


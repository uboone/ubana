#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  ./merge_arborist_into_ntuple.sh <arborist_out.root> <source_ntuple.root> <output.root>

Example:
  ./merge_arborist_into_ntuple.sh arborist_20251219T195436_out.root random_ncpi0_ntuple_file.root random_ncpi0_with_arborist.root

What it does:
  - Copies <source_ntuple.root> -> <output.root> (preserves all original directories/objects)
  - Adds/overwrites a top-level TTree named:
      spline_weights
    with branches:
      run (int), subrun (int), event (int), entry (long64)
      one branch per mcweight key (each is vector<double>)
    filled from Arborist's arborist/eventweight_tree branch: mcweight.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -ne 3 ]]; then
  usage
  exit 2
fi

ARB_FILE="$1"
SRC_FILE="$2"
OUT_FILE="$3"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Look for the macro in the script dir first (dev environment),
# then fall back to the installed source location (via UBANA_DIR).
if [[ -f "${SCRIPT_DIR}/merge_arborist_into_ntuple.C" ]]; then
  MACRO="${SCRIPT_DIR}/merge_arborist_into_ntuple.C"
elif [[ -n "${UBANA_DIR:-}" && -f "${UBANA_DIR}/source/ubana/TreeReader/merge_arborist_into_ntuple.C" ]]; then
  MACRO="${UBANA_DIR}/source/ubana/TreeReader/merge_arborist_into_ntuple.C"
else
  echo "ERROR: Cannot find merge_arborist_into_ntuple.C" >&2
  echo "Looked in: ${SCRIPT_DIR}/" >&2
  [[ -n "${UBANA_DIR:-}" ]] && echo "       and: ${UBANA_DIR}/source/ubana/TreeReader/" >&2
  exit 3
fi

ROOT_EXIT_CODE=0
root -l -b <<EOF || ROOT_EXIT_CODE=$?
.L ${MACRO}
int rc = merge_arborist_into_ntuple("${ARB_FILE}","${SRC_FILE}","${OUT_FILE}");
gSystem->Exit(rc);
EOF

if [[ "${ROOT_EXIT_CODE}" -ne 0 ]]; then
  echo "ERROR: ROOT macro failed with exit code ${ROOT_EXIT_CODE}" >&2
  exit 4
fi

# Verify the spline_weights tree exists in the output
echo "Verifying spline_weights tree in output..."
VERIFY_CODE=0
root -l -b -q <<VERIFY_EOF || VERIFY_CODE=$?
{
  TFile f("${OUT_FILE}", "READ");
  if (f.IsZombie()) {
    std::cerr << "ERROR: Cannot open output file: ${OUT_FILE}" << std::endl;
    gSystem->Exit(10);
  }
  TTree* t = (TTree*)f.Get("spline_weights");
  if (!t) {
    std::cerr << "ERROR: spline_weights tree NOT FOUND in output file!" << std::endl;
    f.ls();
    gSystem->Exit(11);
  }
  Long64_t nentries = t->GetEntries();
  if (nentries == 0) {
    std::cerr << "ERROR: spline_weights tree has 0 entries!" << std::endl;
    gSystem->Exit(12);
  }
  std::cout << "OK: spline_weights tree found with " << nentries << " entries" << std::endl;
  gSystem->Exit(0);
}
VERIFY_EOF

if [[ "${VERIFY_CODE}" -ne 0 ]]; then
  echo "ERROR: Verification failed - spline_weights tree missing or empty (code ${VERIFY_CODE})" >&2
  exit 5
fi

echo "Merge completed successfully."



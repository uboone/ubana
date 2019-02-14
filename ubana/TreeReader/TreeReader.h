#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"

#include <fstream>
#include <vector>
#include <map>

#include "TreeInterface.h"

class TH1D;
class TTree;
class TFile;

namespace uboone {

class TreeReader {
public:
  /** Constructor for an InputSource. */
  TreeReader(fhicl::ParameterSet const &pset,
             art::ProductRegistryHelper &helper,
             art::SourceHelper const &pm);
  
  /** Required by FileReaderSource. */
  void closeCurrentFile();

  void readFile(std::string const &name,
      	        art::FileBlock* &fb);

  bool readNext(art::RunPrincipal* const &inR,
      	        art::SubRunPrincipal* const &inSR,
      	        art::RunPrincipal* &outR,
      	        art::SubRunPrincipal* &outSR,
      	        art::EventPrincipal* &outE);

  void endJob();

private:
  art::SourceHelper   fSourceHelper;
  art::SubRunID       fCurrentSubRunID;

  bool                fVerbose;
  bool                fDLMode;
  uint32_t            fEventCounter; 
  uint32_t            fEntry;
  int                 fMaxEvents;  //fhicl parameter.  Maximum number of events.
  uint32_t            fSkipEvents; // fhicl parameter.  Number of events to skip.
  std::string         fInputType; //fhicl parameter.  Maximum number of events.
  float               fPOT;

  std::string         fTreeName;
  std::string         fDataProductName;
  fhicl::ParameterSet fBranchDef;

  TreeInterface*      fInterface;
  TFile*              fInputFile; 
};

}  // namespace uboone


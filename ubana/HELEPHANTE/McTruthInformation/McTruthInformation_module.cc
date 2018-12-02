#ifndef MCTRUTHINFORMATION_module
#define MCTRUTHINFORMATION_module

#include "McTruthInformation.h"

// Analyzer class
class McTruthInformation : public art::EDAnalyzer
{
public:
  explicit McTruthInformation(fhicl::ParameterSet const & pset);
  virtual ~McTruthInformation();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void GetTruthInformation(art::Event const & evt);
private:
  // Declare fhiclcpp variables
  std::string fNuLabel;
  std::string fCosmicLabel;
  bool fRecordCosmics;
  bool fVerbose;

  // Declare trees
  TTree *tDataTree;
  int run, subrun, event;
  int int_CCNC, int_mode, int_interactionType, int_target, int_hitNuc, int_hitQuark, nu_pdgCode, lepton_pdgCode;
  std::vector<int> cosmic_pdgCode;
  double int_w, int_x, int_y, int_qSqr, int_pT, int_theta;
  double nu_vX, nu_vY, nu_vZ, nu_T, nu_endX, nu_endY, nu_endZ, nu_endT, nu_pX, nu_pY, nu_pZ, nu_e, nu_p, nu_pT, nu_mass, nu_endPx, nu_endPy, nu_endPz, nu_endE;
  double lepton_vX, lepton_vY, lepton_vZ, lepton_T, lepton_endX, lepton_endY, lepton_endZ, lepton_endT, lepton_pX, lepton_pY, lepton_pZ, lepton_e, lepton_p, lepton_pT, lepton_mass, lepton_endPx, lepton_endPy, lepton_endPz, lepton_endE;
  std::vector<double> cosmic_vX, cosmic_vY, cosmic_vZ, cosmic_T, cosmic_endX, cosmic_endY, cosmic_endZ, cosmic_endT, cosmic_pX, cosmic_pY, cosmic_pZ, cosmic_e, cosmic_p, cosmic_pT, cosmic_mass, cosmic_endPx, cosmic_endPy, cosmic_endPz, cosmic_endE;

  void ClearData();
}; // End class McTruthInformation

McTruthInformation::McTruthInformation(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fNuLabel(pset.get<std::string>("nuLabel")),
    fCosmicLabel(pset.get<std::string>("cosmicLabel")),
    fRecordCosmics(pset.get<bool>("recordCosmics")),
    fVerbose(pset.get<bool>("verbose"))      
{} // END constructor McTruthInformation

McTruthInformation::~McTruthInformation()
{} // END destructor McTruthInformation

void McTruthInformation::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run,"run/I");
  tDataTree->Branch("subrun",&subrun,"subrun/I");
  tDataTree->Branch("event",&event,"event/I");
  tDataTree->Branch("int_CCNC",&int_CCNC);
  tDataTree->Branch("int_mode",&int_mode);
  tDataTree->Branch("int_interactionType",&int_interactionType);
  tDataTree->Branch("int_target",&int_target);
  tDataTree->Branch("int_hitNuc",&int_hitNuc);
  tDataTree->Branch("int_hitQuark",&int_hitQuark);
  tDataTree->Branch("int_w",&int_w);
  tDataTree->Branch("int_x",&int_x);
  tDataTree->Branch("int_y",&int_y);
  tDataTree->Branch("int_qSqr",&int_qSqr);
  tDataTree->Branch("int_pT",&int_pT);
  tDataTree->Branch("int_theta",&int_theta);

  tDataTree->Branch("nu_pdgCode",&nu_pdgCode);
  tDataTree->Branch("nu_vX",&nu_vX);
  tDataTree->Branch("nu_vY",&nu_vY);
  tDataTree->Branch("nu_vZ",&nu_vZ);
  tDataTree->Branch("nu_T",&nu_T);
  tDataTree->Branch("nu_endX",&nu_endX);
  tDataTree->Branch("nu_endY",&nu_endY);
  tDataTree->Branch("nu_endZ",&nu_endZ);
  tDataTree->Branch("nu_endT",&nu_endT);
  tDataTree->Branch("nu_pX",&nu_pX);
  tDataTree->Branch("nu_pY",&nu_pY);
  tDataTree->Branch("nu_pZ",&nu_pZ);
  tDataTree->Branch("nu_e",&nu_e);
  tDataTree->Branch("nu_p",&nu_p);
  tDataTree->Branch("nu_pT",&nu_pT);
  tDataTree->Branch("nu_mass",&nu_mass);
  tDataTree->Branch("nu_endPx",&nu_endPx);
  tDataTree->Branch("nu_endPy",&nu_endPy);
  tDataTree->Branch("nu_endPz",&nu_endPz);
  tDataTree->Branch("nu_endE",&nu_endE);

  tDataTree->Branch("lepton_pdgCode",&lepton_pdgCode);
  tDataTree->Branch("lepton_vX",&lepton_vX);
  tDataTree->Branch("lepton_vY",&lepton_vY);
  tDataTree->Branch("lepton_vZ",&lepton_vZ);
  tDataTree->Branch("lepton_T",&lepton_T);
  tDataTree->Branch("lepton_endX",&lepton_endX);
  tDataTree->Branch("lepton_endY",&lepton_endY);
  tDataTree->Branch("lepton_endZ",&lepton_endZ);
  tDataTree->Branch("lepton_endT",&lepton_endT);
  tDataTree->Branch("lepton_pX",&lepton_pX);
  tDataTree->Branch("lepton_pY",&lepton_pY);
  tDataTree->Branch("lepton_pZ",&lepton_pZ);
  tDataTree->Branch("lepton_e",&lepton_e);
  tDataTree->Branch("lepton_p",&lepton_p);
  tDataTree->Branch("lepton_pT",&lepton_pT);
  tDataTree->Branch("lepton_mass",&lepton_mass);
  tDataTree->Branch("lepton_endPx",&lepton_endPx);
  tDataTree->Branch("lepton_endPy",&lepton_endPy);
  tDataTree->Branch("lepton_endPz",&lepton_endPz);
  tDataTree->Branch("lepton_endE",&lepton_endE);

  tDataTree->Branch("cosmic_pdgCode",&cosmic_pdgCode);
  tDataTree->Branch("cosmic_vX",&cosmic_vX);
  tDataTree->Branch("cosmic_vY",&cosmic_vY);
  tDataTree->Branch("cosmic_vZ",&cosmic_vZ);
  tDataTree->Branch("cosmic_T",&cosmic_T);
  tDataTree->Branch("cosmic_endX",&cosmic_endX);
  tDataTree->Branch("cosmic_endY",&cosmic_endY);
  tDataTree->Branch("cosmic_endZ",&cosmic_endZ);
  tDataTree->Branch("cosmic_endT",&cosmic_endT);
  tDataTree->Branch("cosmic_pX",&cosmic_pX);
  tDataTree->Branch("cosmic_pY",&cosmic_pY);
  tDataTree->Branch("cosmic_pZ",&cosmic_pZ);
  tDataTree->Branch("cosmic_e",&cosmic_e);
  tDataTree->Branch("cosmic_p",&cosmic_p);
  tDataTree->Branch("cosmic_pT",&cosmic_pT);
  tDataTree->Branch("cosmic_mass",&cosmic_mass);
  tDataTree->Branch("cosmic_endPx",&cosmic_endPx);
  tDataTree->Branch("cosmic_endPy",&cosmic_endPy);
  tDataTree->Branch("cosmic_endPz",&cosmic_endPz);
  tDataTree->Branch("cosmic_endE",&cosmic_endE);
} // END function beginJob

void McTruthInformation::endJob()
{
} // END function endJob

void McTruthInformation::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  cosmic_pdgCode.clear();
  cosmic_vX.clear();
  cosmic_vY.clear();
  cosmic_vZ.clear();
  cosmic_T.clear();
  cosmic_endX.clear();
  cosmic_endY.clear();
  cosmic_endZ.clear();
  cosmic_endT.clear();
  cosmic_pX.clear();
  cosmic_pY.clear();
  cosmic_pZ.clear();
  cosmic_e.clear();
  cosmic_p.clear();
  cosmic_pT.clear();
  cosmic_mass.clear();
  cosmic_endPx.clear();
  cosmic_endPy.clear();
  cosmic_endPz.clear();
  cosmic_endE.clear();
} // END function ClearData

void McTruthInformation::GetTruthInformation(art::Event const & evt)
{
  // Prepare handle labels
  art::InputTag nuTag {fNuLabel};
  art::InputTag cosmicTag {fCosmicLabel};

  // Find nu mcTruth
  const auto& mctNuHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(nuTag);
  if ((*mctNuHandle).size()>1) printf("WAIT!!! There are %i MCTruth in this event!", (int) (*mctNuHandle).size());
  for (auto const& mct : (*mctNuHandle))
  {
    const simb::MCNeutrino & mcn = mct.GetNeutrino();
    const simb::MCParticle & nu = mcn.Nu();
    const simb::MCParticle & lepton = mcn.Lepton();

    // Interaction info
    int_CCNC = mcn.CCNC();
    int_mode = mcn.Mode();
    int_interactionType = mcn.InteractionType();
    int_target = mcn.Target();
    int_hitNuc = mcn.HitNuc();
    int_hitQuark = mcn.HitQuark();
    int_w = mcn.W();
    int_x = mcn.X();
    int_y = mcn.Y();
    int_qSqr = mcn.QSqr();
    int_pT = mcn.Pt();
    int_theta = mcn.Theta();

    // Incoming neutrino info
    nu_pdgCode = nu.PdgCode();
    nu_vX = nu.Vx();
    nu_vY = nu.Vy();
    nu_vZ = nu.Vz();
    nu_T = nu.T();
    nu_endX = nu.EndX();
    nu_endY = nu.EndY();
    nu_endZ = nu.EndZ();
    nu_endT = nu.EndT();
    nu_pX = nu.Px();
    nu_pY = nu.Py();
    nu_pZ = nu.Pz();
    nu_e = nu.E();
    nu_p = nu.P();
    nu_pT = nu.Pt();
    nu_mass = nu.Mass();
    nu_endPx = nu.EndPx();
    nu_endPy = nu.EndPy();
    nu_endPz = nu.EndPz();
    nu_endE = nu.EndE();

    // Outgoing lepton info
    lepton_pdgCode = lepton.PdgCode();
    lepton_vX = lepton.Vx();
    lepton_vY = lepton.Vy();
    lepton_vZ = lepton.Vz();
    lepton_T = lepton.T();
    lepton_endX = lepton.EndX();
    lepton_endY = lepton.EndY();
    lepton_endZ = lepton.EndZ();
    lepton_endT = lepton.EndT();
    lepton_pX = lepton.Px();
    lepton_pY = lepton.Py();
    lepton_pZ = lepton.Pz();
    lepton_e = lepton.E();
    lepton_p = lepton.P();
    lepton_pT = lepton.Pt();
    lepton_mass = lepton.Mass();
    lepton_endPx = lepton.EndPx();
    lepton_endPy = lepton.EndPy();
    lepton_endPz = lepton.EndPz();
    lepton_endE = lepton.EndE();

  } // End of mctruth loop

  // Find cosmic mcTruth
  if ( fRecordCosmics )
  {
    const auto& mctCosmicHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(cosmicTag);
    if ((*mctCosmicHandle).size()>1) printf("WAIT!!! There are %i cosmic MCTruth in this event!", (int) (*mctCosmicHandle).size());
    for (auto const& mct : (*mctCosmicHandle))
    {
      int nCosmics = mct.NParticles();
      for (int j=0; j<nCosmics; j++)
      {
        const simb::MCParticle& cosmic = mct.GetParticle(j);
        cosmic_pdgCode.push_back(cosmic.PdgCode());
        cosmic_vX.push_back(cosmic.Vx());
        cosmic_vY.push_back(cosmic.Vy());
        cosmic_vZ.push_back(cosmic.Vz());
        cosmic_T.push_back(cosmic.T());
        cosmic_endX.push_back(cosmic.EndX());
        cosmic_endY.push_back(cosmic.EndY());
        cosmic_endZ.push_back(cosmic.EndZ());
        cosmic_endT.push_back(cosmic.EndT());
        cosmic_pX.push_back(cosmic.Px());
        cosmic_pY.push_back(cosmic.Py());
        cosmic_pZ.push_back(cosmic.Pz());
        cosmic_e.push_back(cosmic.E());
        cosmic_p.push_back(cosmic.P());
        cosmic_pT.push_back(cosmic.Pt());
        cosmic_mass.push_back(cosmic.Mass());
        cosmic_endPx.push_back(cosmic.EndPx());
        cosmic_endPy.push_back(cosmic.EndPy());
        cosmic_endPz.push_back(cosmic.EndPz());
        cosmic_endE.push_back(cosmic.EndE());
      } // End of mcpart in mctruth loop
    } // End of mctruth loop
  } // End of if recording cosmics

  // Diagnostic information
  std::map<int, std::string> CCNCTable;
  CCNCTable[0] = std::string("CC");
  CCNCTable[1] = std::string("NC");
  if (fVerbose)
  {
    printf("Found a %s %i interaction with a %i incoming neutrino and %i outgoing lepton.\n", CCNCTable[int_CCNC].c_str(), int_interactionType, nu_pdgCode, lepton_pdgCode);
    if (fRecordCosmics ) printf("Event contains %i cosmics.\n", (int) cosmic_pdgCode.size());
    else printf("Not analyzing cosmics.\n");
  }
  return;
} // END function GetTruthParticles


void McTruthInformation::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  if (fVerbose) printf("\n|-----------------------------------------------------|");
  if (fVerbose) printf("\n|   MCTRUTHINFORMATION MODULE                         |");
  if (fVerbose) printf("\n|-----------------------------------------------------|\n\n");  
  
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  if (fVerbose) printf("||  INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]  ||\n", event, run, subrun);


  // Assign mcTruth values to global variables
  GetTruthInformation(evt);

  // Fill tree and finish event loop
  tDataTree->Fill();
  if (fVerbose) printf("\n|-----------------------------------------------------|\n\n");  
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(McTruthInformation)

#endif // END def McTruthInformation_module
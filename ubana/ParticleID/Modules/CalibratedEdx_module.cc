////////////////////////////////////////////////////////////////////////
// Class:       UBPID::CalibratedEdx
// Plugin Type: producer (art v2_05_01)
// File:        UBPID::CalibratedEdx_module.cc
//
// Generated at Thu Jun  7 10:24:02 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 *
 * \class ParticleId
 *
 * \brief Calibration producer module
 *
 * \author Kirsty Duffy (kduffy@fnal.gov), Adam Lister (a.lister1@lancaster.ac.uk)
 *
 * \date 2018/04/18
 *
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

// local includes
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"

// Other stuff to select neutrino slices
#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "Pandora/PdgTable.h"

#include "TRandom3.h"

#include <memory>

namespace UBPID{
  class CalibratedEdx;
}

class UBPID::CalibratedEdx : public art::EDProducer {
  public:
    explicit CalibratedEdx(fhicl::ParameterSet const & p);

    // Plugins should not be copied or assigned.
    CalibratedEdx(UBPID::CalibratedEdx const &) = delete;
    CalibratedEdx(UBPID::CalibratedEdx &&) = delete;
    CalibratedEdx & operator = (UBPID::CalibratedEdx const &) = delete;
    CalibratedEdx & operator = (UBPID::CalibratedEdx &&) = delete;
    
    typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

    // Required functions.
    void produce(art::Event & e) override;

  private:

    // Declare member data here.
    std::string fCaloLabel;
    std::string fTrackLabel;
    std::string fPFParticleLabel;

    bool fIsSimSmear;
    bool fIsSetSeed;
    int fRandomSeed;
    std::vector<double> fSimGausSmearWidth;

    bool fIsDataNewRecomb;
    double fBoxrecomb_betap;
    double fBoxrecomb_alpha;
    double fBoxrecomb_rho;
    double fBoxrecomb_wion;
    double fBoxrecomb_efield;
    double fBoxrecomb_ADCtoe;

};


UBPID::CalibratedEdx::CalibratedEdx(fhicl::ParameterSet const & p)
  // :
  // Initialize member data here.
{

  fhicl::ParameterSet p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");

  fCaloLabel = p_labels.get< std::string > ("CalorimetryLabel");
  fTrackLabel = p_labels.get< std::string > ("TrackLabel");
  fPFParticleLabel = p_labels.get<std::string>("PFParticleLabel");
  fIsSimSmear = p.get< bool > ("IsSimulationSmearing");
  fSimGausSmearWidth = p.get< std::vector<double> > ("SimulationGausSmearWidth");
  fIsSetSeed = p.get< bool > ("IsSetSeed");
  fRandomSeed = p.get< int > ("RandomSeed");
  fIsDataNewRecomb = p.get<bool> ("IsDataNewRecombination");
  fBoxrecomb_betap = p.get<double> ("NewRecombinationBoxBeta",0.183592);
  fBoxrecomb_alpha = p.get<double>("NewRecombinationBoxAlpha",0.921969);
  fBoxrecomb_rho = p.get<double>("NewRecombinationBoxRho",1.383);
  fBoxrecomb_wion = p.get<double>("NewRecombinationBoxWion",23.6e-6);
  fBoxrecomb_efield = p.get<double>("NewRecombinationBoxEfield",0.273);
  fBoxrecomb_ADCtoe = p.get<double>("ConversionFactordQdxADCtoe",243.);

  produces< std::vector<anab::Calorimetry> >();
  produces< art::Assns< recob::Track, anab::Calorimetry> >();

  std::cout << "[CalibratedEdx] >> Track Label: " << fTrackLabel << std::endl;
  std::cout << "[CalibratedEdx] >> Calo Label: " << fCaloLabel << std::endl;
  std::cout << "[CalibratedEdx] The following is only useful if running simulated data" << std::endl;
  std::cout << "[CalibratedEdx] >> Do simulation smearing? " << fIsSimSmear << std::endl;
  std::cout << "[CalibratedEdx] >> Smearing simulated data by ["
    << fSimGausSmearWidth.at(0) << ", " << fSimGausSmearWidth.at(1)
    << ", " << fSimGausSmearWidth.at(2) << "]" << std::endl;
  std::cout << "[CalibratedEdx] The following is only useful if running real data" << std::endl;
  std::cout << "[CalibratedEdx] >> Apply new recombination model to data? " << fIsDataNewRecomb << std::endl;
  std::cout << "[CalibratedEdx] >> Using Box model parameters " << std::endl
    << "\t betap = " << fBoxrecomb_betap << std::endl
    << "\t alpha = " << fBoxrecomb_alpha << std::endl
    << "\t rho = " << fBoxrecomb_rho << std::endl
    << "\t Wion = " << fBoxrecomb_wion << std::endl
    << "\t Efied = " << fBoxrecomb_efield << std::endl;
  std::cout << "[CalibratedEdx] >> Using dQdx ADC->e correction factor " << fBoxrecomb_ADCtoe << std::endl;

}

void UBPID::CalibratedEdx::produce(art::Event & e)
{

  bool isData = e.isRealData();
  std::cout << "In producer loop" << std::endl;

  // get seed
  int run = e.run();
  int sub_run = e.subRun();
  int event = e.event();

  // setup random
  TRandom3 r;

  if (fIsSetSeed){
    r.SetSeed(fRandomSeed);
  }
  else{
    std::string randomSeedString =  Form("%i%i%i", event, sub_run, run);
    r.SetSeed(std::atoi(randomSeedString.c_str()));
  }

  std::cout << "[CalibratedEdx] Run " << run << " SubRun " << sub_run << " Event " << event << std::endl;
  std::cout << "[CalibratedEdx] The random number seed is " << r.GetSeed() << std::endl;

  std::unique_ptr< std::vector<anab::Calorimetry> > calorimetryCollection( new std::vector<anab::Calorimetry> );
  std::unique_ptr< art::Assns <recob::Track, anab::Calorimetry> > trackCalorimetryAssn( new art::Assns<recob::Track, anab::Calorimetry> );

  // tracks...
  art::Handle < std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr<recob::Track> > trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  // calorimetry object...
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackHandle, e, fCaloLabel);

  ////// Determine if track is in neutrino slice /////////
  //
  // First, get PFParticle handle
  art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
  std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
  art::fill_ptr_vector(pfParticleVector, pfParticleHandle);

  if (!pfParticleHandle.isValid() ) {
      std::cout << "Bad PFParticle handle" << std::endl;
      return;
  }

  // Create PFParicle map
  PFParticleIdMap pfParticleMap;

  for (unsigned int i = 0; i < pfParticleHandle->size(); ++i) {
    const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
    if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second) {
        std::cout << "Coudln't get PFParticleIdMap" << std::endl;
        return;
    }   
  }

  // Separate PFParticles into cosmic and neutrino-induced
  std::vector< art::Ptr<recob::PFParticle> > crParticles;
  std::vector< art::Ptr<recob::PFParticle> > nuParticles;

  for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it) {
        const art::Ptr<recob::PFParticle> pParticle(it->second);

        // Only look for primary particles
        if (!pParticle->IsPrimary()) continue;

        // Check if this particle is identified as the neutrino
        const int pdg(pParticle->PdgCode());
        const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);


        // If it is, lets get the vertex position
        /*
        if(isNeutrino){
            this->GetVertex(pfParticlesToVerticesMap, pParticle );

        }
        */

        // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
        if (!isNeutrino)
        {
            crParticles.push_back(pParticle);
            continue;
        }

        // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
        //       If this is not the case please handle accordingly
        if (!nuParticles.empty())
        {
            //throw cet::exception("SinglePhoton") << "  This event contains multiple reconstructed neutrinos!";
            std::cerr << "Exception" << std::endl;
        }

        // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
        for (const size_t daughterId : pParticle->Daughters())
        {
            if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                //throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";
                std::cerr << "Exception" << std::endl;

            nuParticles.push_back(pfParticleMap.at(daughterId));
        }
    }

    std::cout << "NuParticles: " << nuParticles.size() << std::endl;
    std::cout << "CrParticles: " << crParticles.size() << std::endl;

    // Get PFP <-> Track associations
    art::FindOneP< recob::PFParticle > pfParticlePerTrack(trackHandle, e, fPFParticleLabel);
    art::FindOneP< recob::Track > trackPerPFParticle(pfParticleHandle, e, fTrackLabel);
    std::map< art::Ptr<recob::Track>, art::Ptr<recob::PFParticle> > trackToPFParticleMap;
    std::map< art::Ptr<recob::PFParticle>, art::Ptr<recob::Track> > pfParticleToTrackMap;

    // Fill track -> pfparticle map
    for (size_t i = 0; i < trackCollection.size(); i++) {
        art::Ptr<recob::Track> thisTrack = trackCollection.at(i);
        trackToPFParticleMap[thisTrack] = pfParticlePerTrack.at(thisTrack.key());
    }

    // Fill pfparticle -> track map
    for (size_t i = 0; i < pfParticleVector.size(); i++) {
        art::Ptr<recob::PFParticle> thisPFP = pfParticleVector.at(i);
        pfParticleToTrackMap[thisPFP] = trackPerPFParticle.at(thisPFP.key());
    }

    // Find neutrino tracks in neutrino slice
    std::vector< art::Ptr<recob::Track> > nuTrackCollection;
    for (art::Ptr<recob::Track> track : trackCollection) {
        art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];
        for (art::Ptr<recob::PFParticle> &p : nuParticles ) {
            if (pfp == p) {
                std::cout << "Found pfp in nu slice" << std::endl;
                nuTrackCollection.push_back(pfParticleToTrackMap[p]);
                //isNu = true;
            }
        }
    }
  
  //bool isNu = false;
  //for (auto& track : trackCollection){
  for (art::Ptr<recob::Track> track : nuTrackCollection){

    /*
    art::Ptr<recob::PFParticle> pfp = trackToPFParticleMap[track];
    for (art::Ptr<recob::PFParticle> &p : nuParticles ) {
        if (pfp == p) {
            std::cout << "Found pfp in nu slice" << std::endl;
            //isNu = true;
        }
    }

    if (!isNu) {
        std::cout << "Track is not in neutrino slice. Skipping" << std::endl;
        continue;
    }
    */

    if (!caloFromTracks.isValid()){
      std::cout << "[CalibratedEdx] Calorimetry<->Track associations are not valid for this track. Skipping." << std::endl;
      continue;
    }

    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());

    art::Ptr< anab:: Calorimetry > calo;
    int planenum;
    for (auto c : caloFromTrack){
      planenum = c->PlaneID().Plane;
      calo = c;

      std::vector<double> dEdx = calo->dEdx();
      std::vector<double> dQdx = calo->dQdx();

      if (!calo || planenum < 0 || planenum > 2){
        std::cout << "[CalibratedEdx] Calorimetry on plane " << planenum << " is unavailable. Not smearing or applying new recombination." << std::endl;
        //continue;
        // Don't continue, just do nothing and return original calo object, with no smearing/recombination corrections.
      }
      else{
        /**
         * If this is simulated data, then we want to smear the dE/dx values
         * by a gaussian with a width defined in fcl
         */

        if (isData == false && fIsSimSmear == true){

          for (size_t i = 0; i < dEdx.size(); i++){

            double simulationSmear = r.Gaus(1., fSimGausSmearWidth.at(planenum));
            dEdx.at(i) = dEdx.at(i) * simulationSmear;

          }

        }

        /**
         * If this is real data, then we want to generate new dEdx values
         * based on dQdx and the Box model with parameters defined in fcl
         */

        if (isData==true && fIsDataNewRecomb==true){

          for (size_t i = 0; i < dQdx.size(); i++){

            double dqdx_e = dQdx.at(i) * fBoxrecomb_ADCtoe;

            double dEdx_newrecomb = (exp(dqdx_e*(fBoxrecomb_betap/(fBoxrecomb_rho*fBoxrecomb_efield))*fBoxrecomb_wion)-fBoxrecomb_alpha)/(fBoxrecomb_betap/(fBoxrecomb_rho*fBoxrecomb_efield));

            dEdx.at(i) = dEdx_newrecomb;

          }

        }

      }

      /** build calorimetry object */
      anab::Calorimetry caloObj(
          calo->KineticEnergy(),
          dEdx,
          calo->dQdx(),
          calo->ResidualRange(),
          calo->DeadWireResRC(),
          calo->Range(),
          calo->TrkPitchVec(),
          calo->PlaneID()
          );

      calorimetryCollection->push_back(caloObj);

      util::CreateAssn(*this, e, *calorimetryCollection, track, *trackCalorimetryAssn);


    }


  }

  e.put(std::move(calorimetryCollection));
  e.put(std::move(trackCalorimetryAssn));
}

DEFINE_ART_MODULE(UBPID::CalibratedEdx)

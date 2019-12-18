#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "Objects/CartesianVector.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "ubreco/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "helpers/PandoraInterfaceHelper.h"
//#include "helpers/FlashHelper.h"

#include "TTree.h"

#include <memory>

class SingleTrackFilter;


class SingleTrackFilter : public art::EDFilter {
public:
  explicit SingleTrackFilter(fhicl::ParameterSet const & p);

  SingleTrackFilter(SingleTrackFilter const &) = delete;
  SingleTrackFilter(SingleTrackFilter &&) = delete;
  SingleTrackFilter & operator = (SingleTrackFilter const &) = delete;
  SingleTrackFilter & operator = (SingleTrackFilter &&) = delete;

    //for charge stuff copied from flash id tool
  class Deposition
  {
  public:
    /**
     *  @brief  Default constructor
     *
     *  @param  x the x-component of the charge position
     *  @param  y the z-component of the charge position
     *  @param  z the z-component of the charge position
     *  @param  charge the charge deposited
     *  @param  nPhotons the estimated numer of photons produced
     */
      Deposition(const float x, const float y, const float z, const float charge, const float nPhotons);

    float m_x;        ///< The x-component of the charge position
    float m_y;        ///< The z-component of the charge position
    float m_z;        ///< The z-component of the charge position
    float m_charge;   ///< The charge deposited
    float m_nPhotons; ///< The estimated numer of photons produced
    };


  typedef std::vector<Deposition> DepositionVector;
  bool filter(art::Event & e) override;
  void endJob() override;
  
  bool GetReconstructed(art::Event const &e);
  bool GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                       const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                       const art::FindManyP<anab::ParticleID> &trackPIDAssn);
  bool GetNeutrinoSliceMuonCandidate(int nuPDG,float vtxx,float vtxy,float vtxz , int n_trk,int n_shwr );
  bool GetNeutrinoSliceProtonCandidate(int nuPDG, int n_trk,int n_shwr, float p_startx, float p_starty, float p_startz, float p_endx, float p_endy, float p_endz, float p_length, float pid);
  bool GetNonNeutrinoSliceProtonCandidate(art::Event & e);
  DepositionVector GetDepositionVector(const lar_pandora::PFParticleMap &particleMap, const lar_pandora::PFParticlesToSpacePoints &particlesToSpacePoints,const lar_pandora::SpacePointsToHits &spacePointsToHits, const art::Ptr<recob::PFParticle> &particle);
  pandora::CartesianVector GetChargeWeightedCenter(const DepositionVector &depositionVector);
  flashana::QCluster_t GetLightCluster(const DepositionVector &depositionVector);
  float GetTotalCharge(const DepositionVector &depositionVector);
  bool FlashMatched(float deltay, float deltaz, float deltaysigma, float deltazsigma, float xclvariable, float flashscore,float p_startx, float p_starty, float p_startz, float p_endx, float p_endy, float p_endz, float p_length, float pid);
  void clearEvent();

private:
  
  PandoraInterfaceHelper pandoraInterfaceHelper;
  // TrackHelper trackHelper;
  
  // LAr Pandora Helper fields
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::TrackVector pftracks, trackVector2;
  lar_pandora::ShowerVector pfshowers;
  lar_pandora::VertexVector pfvertices;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::PFParticlesToTracks particlesToTracks;
  lar_pandora::TracksToHits tracksToHits;
  lar_pandora::PFParticlesToClusters particlesToClusters;
  lar_pandora::PFParticlesToShowers particlesToShowers;
  lar_pandora::PFParticlesToHits pfParticlesToHits;
  lar_pandora::HitsToPFParticles hitsToPfParticles;
  lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints;
  lar_pandora::ClustersToHits clustersToHits;
  lar_pandora::SpacePointVector spacePoints;
  lar_pandora::HitsToSpacePoints hitsToSpacePoints;
  lar_pandora::SpacePointsToHits spacePointsToHits;

  // producer datalabels
  std::string m_pfp_producer;
  std::string m_hitfinder_producer;
  std::string m_flash_producer;
  bool m_muon;
  bool m_proton;
  bool _debug;
  double m_beamwindow_low;
  double m_beamwindow_high;
  // variables for cut
 
  int  _nu_n_pfp,_n_track,_n_shower, _nu_PDG;
  float _nu_vtxx,_nu_vtxy,_nu_vtxz; 

  float _p_startx, _p_starty, _p_startz;
  float _p_endx, _p_endy, _p_endz;
  float _p_length,_p_chi2_p_2;
  float _flash_score;
  float _deltaY, _deltaZ, _deltaYSigma, _deltaZSigma, _xclVariable;

  // cut values for both muon and proton
  int m_nu_n_pfp;
  int m_n_track;
  int m_n_shower;
  int m_nu_PDG;
 
  //cuts for muon vertex
  float m_start_x; //vertex start
  float m_start_y;
  float m_start_z;
  float m_end_x; //vertex end
  float m_end_y;
  float m_end_z;

  //cut for protons only
  float m_p_start_x;
  float m_p_start_y;
  float m_p_start_z;
  float m_p_end_x;
  float m_p_end_y;
  float m_p_end_z;
  float m_p_length_min,m_p_length_max;
  float m_p_chi2_p_2;
  float m_flash_score;
  float m_deltaY, m_deltaZ, m_deltaYSigma, m_deltaZSigma,m_xclVariable_low, m_xclVariable_high;
  
};


SingleTrackFilter::SingleTrackFilter(fhicl::ParameterSet const & pset) :
    EDFilter(pset)
{
    //producer datalabels
    m_pfp_producer = pset.get<std::string>("pfp_producer", "pandora");
    m_hitfinder_producer = pset.get<std::string>("hitfinder_producer", "gaushit");
    m_flash_producer    = pset.get<std::string>("flash_producer", "simpleFlashBeam");
    m_muon         = pset.get<bool>("DoMuon", true);
    m_proton         = pset.get<bool>("DoProton", true);
    _debug              = pset.get<bool>("DebugMode","false");
    //cut variables
    m_nu_n_pfp = pset.get<int>("NPFP", 1);
    m_n_track = pset.get<int>("NTrack", 1);
    m_n_shower = pset.get<int>("NShower", 0);
    m_nu_PDG = pset.get<int>("nuPDG", 14);
    m_start_x = pset.get<double>("StartX", -9999);
    m_start_y = pset.get<double>("StartY", -9999);
    m_start_z = pset.get<double>("StartZ", -9999);
    m_end_x = pset.get<double>("Endx", 9999);
    m_end_y = pset.get<double>("Endy", 9999);
    m_end_z = pset.get<double>("Endz", 9999);

    //for protons
    m_p_start_x = pset.get<double>("ProtonStartX", -9999);
    m_p_start_y = pset.get<double>("ProtonStartY", -9999);
    m_p_start_z = pset.get<double>("ProtonStartZ", -9999);
    m_p_end_x = pset.get<double>("ProtonEndx", 9999);
    m_p_end_y = pset.get<double>("ProtonEndy", 9999);
    m_p_end_z = pset.get<double>("ProtonEndz", 9999);
    m_p_length_min = pset.get<double>("ProtonLengthMin", 9999);
    m_p_length_max = pset.get<double>("ProtonLengthMax", 9999);
    m_p_chi2_p_2 = pset.get<double>("ProtonChi2", -9999);
    m_deltaY    = pset.get<double>("DeltaY", 95);
    m_deltaZ    = pset.get<double>("DeltaZ", 95);
    m_deltaYSigma    = pset.get<double>("DeltaYSigma", 95);
    m_deltaZSigma    = pset.get<double>("DeltaZSigma", 95);
    m_xclVariable_low    = pset.get<double>("xclVariableLow", 95);
    m_xclVariable_high    = pset.get<double>("xclVariableHigh", 95);
    m_flash_score = pset.get<double>("FlashScore", 10);
    m_beamwindow_low    = pset.get<double>("Beam_low", 3.57);
  m_beamwindow_high   = pset.get<double>("Beam_high", 5.25);
    //   verbose_ = pset.get<int>("verbose");
}

bool SingleTrackFilter::filter(art::Event & e)
{
  //if(_debug) 
std::cout<<"start filtering...."<<std::endl;
  clearEvent();
  larpandora.CollectPFParticleMetadata(e, m_pfp_producer, pfparticles, particlesToMetadata);
   larpandora.BuildPFParticleMap(pfparticles, particleMap);
  larpandora.CollectTracks(e,m_pfp_producer, trackVector2, tracksToHits);
  larpandora.CollectTracks(e, m_pfp_producer, pftracks, particlesToTracks);
  larpandora.CollectShowers(e, m_pfp_producer, pfshowers, particlesToShowers);
  larpandora.CollectVertices(e, m_pfp_producer, pfvertices, particlesToVertices);
  larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
  //  larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, "pandoracalipidSCE");
  art::Handle<std::vector<recob::PFParticle>> pfparticles_handle;
  e.getByLabel(m_pfp_producer, pfparticles_handle);
  art::FindManyP<recob::Track> PFPTrackAsso(pfparticles_handle, e, "pandora");
  art::FindManyP<recob::Shower> PFPShowerAsso(pfparticles_handle, e, "pandora");
 
  if(_debug) std::cout<<"=======run/subrun/event: "<<e.id().run()<<" "<<e.id().subRun()<<" "<<e.id().event()<<"========"<<std::endl;
   

  
  if (pfparticles.size() == 0)
    {
    return false;
    }
  else
    {
      std::cout<<"pfneutrinos.size()"<<pfneutrinos.size()<<std::endl;
      if (pfneutrinos.size() == 1)
	{
	  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
	  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
	  const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
	  _nu_vtxx = neutrino_vtx.X();
	  _nu_vtxy = neutrino_vtx.Y();
	  _nu_vtxz = neutrino_vtx.Z();
	  _nu_PDG = pfnu->PdgCode();
	  const std::vector<size_t> &daughterIDs = pfnu->Daughters();
          _nu_n_pfp=daughterIDs.size();
	  _n_track = 0;
	  _n_shower = 0;
	  // count number of tracks and showers as neutrino daughter in neutrino slice
	  for(int j = 0; j< _nu_n_pfp; j++)
	    {
	   
	      auto Iterator = particleMap.find(pfnu->Daughters().at(j));
	      auto this_pfp = Iterator->second;
	      auto assoTrack = PFPTrackAsso.at(this_pfp.key());
	      auto assoShower = PFPShowerAsso.at(this_pfp.key());
	      if(this_pfp->IsPrimary()) continue;
	      if(assoTrack.size()==1)
		{
		  _n_track++;
		}
	      if(assoShower.size()==1)
		{
		  _n_shower++;
		}
	    }
	  std::cout<<"pfp/track/shower "<<_nu_n_pfp<<" "<<_n_track<<" "<<_n_shower<<std::endl;
	  // get track information if there's only 1 track in the neutrino slice
	  // pid, length,start, end
	  if(_n_track==1)
	    {
	      for (unsigned int n = 0; n < pfparticles.size(); ++n)
                {
		  const art::Ptr<recob::PFParticle> particle = pfparticles.at(n);
		  if(!particle->IsPrimary())
		    {
		      const auto parentIterator = particleMap.find(particle->Parent());
		      const int parentPDG = std::abs(parentIterator->second->PdgCode());
		      if(abs(parentPDG) != 14) continue;
		      
		      lar_pandora::PFParticlesToTracks::const_iterator trkIter = particlesToTracks.find(particle);
		      if (particlesToTracks.end() != trkIter)
			{
			  const lar_pandora::TrackVector &pftracks = trkIter->second;
			  if (!pftracks.empty())
			    {
			      if (pftracks.size() !=1)
				{
				  std::cout << " Warning: Found particle with more than one associated track "<<pftracks.size() << std::endl;
				  continue;
				}
			      const art::Ptr<recob::Track> track = *(pftracks.begin());
			      const auto &trackVtxPosition = track->Vertex();
			      const auto &trackEndPosition = track->End();
			      _p_startx=trackVtxPosition.x();
			      _p_starty=trackVtxPosition.y();
			      _p_startz=trackVtxPosition.z();
			      _p_endx=trackEndPosition.x();
			      _p_endy=trackEndPosition.y();
			      _p_endz=trackEndPosition.z();
			      _p_length=track->Length();
			      std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track.key());
			      if (trackPID.size() != 0)
				{
				  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
				  for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
				    {
				      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
				      int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
				      if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward)
					{
					  if (AlgScore.fAlgName == "Chi2"&&TMath::Abs(AlgScore.fAssumedPdg) == 2212)
					    {
					      if(planenum==2)    _p_chi2_p_2=AlgScore.fValue;
					    }
					}
				    }
				}//end of pid
			    }
			}
		    }
		}
	    
	    }//end of single track info
	  
	  if(m_muon)
            { 
              if(GetNeutrinoSliceMuonCandidate(_nu_PDG,_nu_vtxx,_nu_vtxy,_nu_vtxz, _n_track, _n_shower))
                { return true;}
              else
                { return false;}
            }
	  if(m_proton)
	    {
	      if(_debug) 
		{
		  std::cout<<"-------------neutrino----------------"<<std::endl;
		  std::cout<<" parent PDG: "<<_nu_PDG<<std::endl;
		  std::cout<<" # of track: "<<_n_track<<std::endl;
		  std::cout<<" # of shower: "<<_n_shower<<std::endl;
		  std::cout<<" start x: "<<_p_startx<<std::endl;
		  std::cout<<" start y: "<<_p_starty<<std::endl;
		  std::cout<<" start z: "<<_p_startz<<std::endl;
		  std::cout<<" end x: "<<_p_endx<<std::endl;
		  std::cout<<" end y: "<<_p_endy<<std::endl;
		  std::cout<<" end z: "<<_p_endz<<std::endl;
		  std::cout<<" length: "<<_p_length<<std::endl;
		  std::cout<<" pid: "<<_p_chi2_p_2<<std::endl;
		}
	      if(GetNeutrinoSliceProtonCandidate(_nu_PDG, _n_track, _n_shower, _p_startx,_p_starty,_p_startz,_p_endx,_p_endy,_p_endz,_p_length,_p_chi2_p_2)) return true;
//	   }//end of neutrino slice
              else //for non-neutrino slice protons
	        {
	          std::cout<<"non-neutrino"<<std::endl;
	          if(m_proton && GetNonNeutrinoSliceProtonCandidate(e)) return true;
	        }// end of non-neutrino slice protons
    
            } // end of proton session
       
        } // end of pfneutrino == 1
     
    }//end of pfparticle.size>0
  
  return false;

}


bool  SingleTrackFilter::GetNeutrinoSliceMuonCandidate(int nuPDG,float vtxx,float vtxy,float vtxz , int n_trk,int n_shwr )
 {
   if(abs(nuPDG) != m_nu_PDG) return false;
   if (vtxx < m_start_x) return false;
   if (vtxx > m_end_x) return false;
   if (vtxy < m_start_y) return false;
   if (vtxy > m_end_y) return false;
   if (vtxz < m_start_z) return false;
   if (vtxz > m_end_z) return false;
   if (n_trk !=m_n_track) return false;
   if (n_shwr !=m_n_shower) return false;
   else return true;
 }

bool  SingleTrackFilter::GetNeutrinoSliceProtonCandidate(int nuPDG, int n_trk,int n_shwr, float p_startx, float p_starty, float p_startz, float p_endx, float p_endy, float p_endz, float p_length, float pid)
{
   if(abs(nuPDG) != m_nu_PDG) return false;
   if (p_startx < m_p_start_x) return false;
   if (p_startx > m_p_end_x) return false;
   if (p_starty < m_p_start_y) return false;
   if (p_starty > m_p_end_y) return false;
   if (p_startz < m_p_start_z) return false;
   if (p_startz > m_p_end_z) return false;
   if (p_endx < m_p_start_x) return false;
   if (p_endx > m_p_end_x) return false;
   if (p_endy < m_p_start_y) return false;
   if (p_endy > m_p_end_y) return false;
   if (p_endz < m_p_start_z) return false;
   if (p_endz > m_p_end_z) return false;
   if (p_length < m_p_length_min) return false;
  if (p_length > m_p_length_max) return false;
   if (pid > m_p_chi2_p_2) return false;
   if (n_trk !=m_n_track) return false;
   if (n_shwr !=m_n_shower) return false;
   else 
     {
       if(_debug) std::cout<<"got a proton candidate from neutrino slice"<<std::endl;
       return true;
     }
}
bool SingleTrackFilter::FlashMatched(float deltay, float deltaz, float deltaysigma, float deltazsigma, float xclvariable, float flashscore,float p_startx, float p_starty, float p_startz, float p_endx, float p_endy, float p_endz, float p_length, float pid)
{
  if (p_startx < m_start_x) return false;
  if (p_startx > m_end_x) return false;
  if (p_starty < m_start_y) return false;
  if (p_starty > m_end_y) return false;
  if (p_startz < m_start_z) return false;
  if (p_startz > m_end_z) return false;
  if (p_endx < m_start_x) return false;
  if (p_endx > m_end_x) return false;
  if (p_endy < m_start_y) return false;
  if (p_endy > m_end_y) return false;
  if (p_endz < m_start_z) return false;
  if (p_endz > m_end_z) return false;
  if (p_length < m_p_length_min) return false;
  if (p_length > m_p_length_max) return false;
  if (pid > m_p_chi2_p_2) return false;
  if (abs(deltay)>m_deltaY) return false;
  if (abs(deltaz)>m_deltaZ) return false;
  if (abs(deltaysigma)>m_deltaYSigma) return false;
  if (abs(deltazsigma)>m_deltaZSigma) return false;
  if (xclvariable<m_xclVariable_low) return false;
  if (xclvariable>m_xclVariable_high) return false;
  if (flashscore>m_flash_score) return false;
  else return true;
}

bool SingleTrackFilter::GetNonNeutrinoSliceProtonCandidate(art::Event & e)
{
  //flash info
  clearEvent();
  art::Handle<std::vector<recob::OpFlash>> flash_handle_beam;
  std::vector<art::Ptr<recob::OpFlash> > flash_beam_vec;
  if ( e.getByLabel(m_flash_producer,flash_handle_beam)) { art::fill_ptr_vector(flash_beam_vec,flash_handle_beam); }
  
  float maxTotalPE=0;
  float flashbeam_Ywidth=-999;
  float flashbeam_Zwidth=-999;
  float flashbeam_Ycenter=-999;
  float flashbeam_Zcenter=-999;
  float flashbeam_Time=-999;
  float flashbeam_TotalPE=-999;
  float flash_brightest_Ywidth=-999;
  float flash_brightest_Zwidth=-999;
  float flash_brightest_Ycenter=-999;
  float flash_brightest_Zcenter=-999;
  float flash_brightest_TotalPE=-999;

  for (int i_fb = 0; i_fb < int(flash_beam_vec.size()); ++i_fb)
    {
      art::Ptr<recob::OpFlash> CurrentBeamFlash = flash_beam_vec.at(i_fb);
      flashbeam_Ywidth=CurrentBeamFlash->YWidth();
      flashbeam_Zwidth=CurrentBeamFlash->ZWidth();
      flashbeam_Ycenter=CurrentBeamFlash->YCenter();
      flashbeam_Zcenter=CurrentBeamFlash->ZCenter();
      flashbeam_Time=CurrentBeamFlash->Time();
      flashbeam_TotalPE=CurrentBeamFlash->TotalPE();
      if(flashbeam_Time>m_beamwindow_low&&flashbeam_Time<m_beamwindow_high)
	{
	  if(flashbeam_TotalPE<maxTotalPE)
	    {
	      if(_debug) std::cout<<"[flash!] "<<CurrentBeamFlash->TotalPE()<<" is not as bright as "<<maxTotalPE<<std::endl;
	    }
	  else
	    {
	      maxTotalPE=flashbeam_TotalPE;
	      flash_brightest_TotalPE=flashbeam_TotalPE;
	      flash_brightest_Zcenter=flashbeam_Zcenter;
	      flash_brightest_Ycenter=flashbeam_Ycenter;
	      flash_brightest_Zwidth= flashbeam_Zwidth;
	      flash_brightest_Ywidth= flashbeam_Ywidth;
	    }
	}
    }
  
  ///
larpandora.CollectPFParticleMetadata(e, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  larpandora.CollectTracks(e, m_pfp_producer, pftracks, particlesToTracks);
  larpandora.CollectPFParticles(e, m_pfp_producer, pfparticles, particlesToSpacePoints);
  larpandora.CollectSpacePoints(e, m_pfp_producer, spacePoints, spacePointsToHits);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, "pandoracalipidSCE");
  art::Handle<std::vector<recob::PFParticle>> pfparticles_handle;
  e.getByLabel(m_pfp_producer, pfparticles_handle);

  art::FindManyP<anab::T0> nuFlashScoreAsso(pfparticles_handle, e, "flashmatch");
    
  int np = 0;
  std::cout<<pfparticles.size()<< "pfps in total"<<std::endl;
  for (unsigned int n = 0; n < pfparticles.size(); ++n)
    {
      _deltaY =-9999;
      _deltaZ =-9999;
      _deltaYSigma =-9999;
      _deltaZSigma =-9999;
      _xclVariable =-9999;
      _flash_score = 9999;
      _p_startx =-9999;
      _p_starty =-9999;
      _p_startz =-9999;
      _p_endx =-9999;
      _p_endy =-9999;
      _p_endz =-9999;
      _p_length =-9999;
      _p_chi2_p_2=9999;
  
      const art::Ptr<recob::PFParticle> particle = pfparticles.at(n);
      if(particle->IsPrimary())
	{
	  lar_pandora::PFParticlesToTracks::const_iterator trkIter = particlesToTracks.find(particle);
	  const std::vector<art::Ptr<anab::T0>> T0_flashchi_v2 = nuFlashScoreAsso.at(particle.key());
	  if (T0_flashchi_v2.size() == 1)
	    {
	      _flash_score = T0_flashchi_v2.at(0)->TriggerConfidence();
	    }
	  const auto chargeDeposition=GetDepositionVector(particleMap, particlesToSpacePoints, spacePointsToHits, particle);
	  float m_totalCharge = GetTotalCharge(chargeDeposition);
	  if( m_totalCharge>0)
	    {
	      const auto chargeCenter(this->GetChargeWeightedCenter(chargeDeposition));
	      float  m_centerX = chargeCenter.GetX();
	      float  m_centerY = chargeCenter.GetY();
	      float  m_centerZ = chargeCenter.GetZ();
	      
	      _deltaY=m_centerY-flash_brightest_Ycenter;
	      _deltaZ=m_centerZ-flash_brightest_Zcenter;
	      if(flash_brightest_Ywidth>0)
		{
		  _deltaYSigma=_deltaY/flash_brightest_Ywidth;
		}
	      if(flash_brightest_Zwidth>0)
		{
		  _deltaZSigma=_deltaZ/flash_brightest_Ywidth;
		}
	      
	      float _chargeToLightRatio=m_totalCharge/flash_brightest_TotalPE;
	      _xclVariable=270*log10(_chargeToLightRatio) - m_centerX;//CoefXCL= 270 from flash_neutrino_id.fcl
	    }
	  if (particlesToTracks.end() != trkIter)
	    {
	      const lar_pandora::TrackVector &pftracks = trkIter->second;
	      if (!pftracks.empty())
		{
		  if (pftracks.size() !=1)
		    {
		      std::cout << " Warning: Found particle with more than one associated track "<<pftracks.size() << std::endl;
		      continue;
		    }
		  const art::Ptr<recob::Track> track = *(pftracks.begin());
		  const auto &trackVtxPosition = track->Vertex();
		  const auto &trackEndPosition = track->End();
		  _p_startx=trackVtxPosition.x();
		  _p_starty=trackVtxPosition.y();
		  _p_startz=trackVtxPosition.z();
		  _p_endx=trackEndPosition.x();
		  _p_endy=trackEndPosition.y();
		  _p_endz=trackEndPosition.z();
		  _p_length=track->Length();
		  std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track.key());
		  if (trackPID.size() != 0)
		    {
		      std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
		      for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
			{
			  anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
			  int planenum = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);
			  if (anab::kVariableType(AlgScore.fVariableType) == anab::kGOF && anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward)
			    {
			      if (AlgScore.fAlgName == "Chi2"&&TMath::Abs(AlgScore.fAssumedPdg) == 2212)
				{
				  if(planenum==2)    _p_chi2_p_2=AlgScore.fValue;
				}
			    }
			}
		    }//end of pid
		}//!pftracks.empty()
	    }//particlesToTracks.end() != trkIter
	}//particle->IsPrimary()
     
      if(FlashMatched(_deltaY, _deltaZ, _deltaYSigma, _deltaZSigma, _xclVariable, _flash_score,_p_startx, _p_starty, _p_startz, _p_endx, _p_endy, _p_endz, _p_length, _p_chi2_p_2))
	{
	  np++;
	 if(_debug) 
	{
	  std::cout<<"-------------non-neutrino----------------"<<std::endl;
	  std::cout<<n<<"th PFParticle"<<std::endl;
	  std::cout<<" start x: "<<_p_startx<<std::endl;
	  std::cout<<" start y: "<<_p_starty<<std::endl;
	  std::cout<<" start z: "<<_p_startz<<std::endl;
	  std::cout<<" end x: "<<_p_endx<<std::endl;
	  std::cout<<" end y: "<<_p_endy<<std::endl;
	  std::cout<<" end z: "<<_p_endz<<std::endl;
	  std::cout<<" length: "<<_p_length<<std::endl;
	  std::cout<<" pid: "<<_p_chi2_p_2<<std::endl;
	  std::cout<<" deltaY : "<<_deltaY<<std::endl;
   	  std::cout<<" deltaZ: "<<_deltaZ<<std::endl;
   	  std::cout<<" deltaYsigma: "<<_deltaYSigma<<std::endl;
   	  std::cout<<" deltaZSigma: "<<_deltaZSigma<<std::endl;
   	  std::cout<<" xclvariable: "<<_xclVariable<<std::endl;
   	  std::cout<<" flash score: "<<_flash_score<<std::endl;
   
	}
	}
    }//end of looping over pfps
  if(np==0) return false;
  else 
    {
      std::cout<<"non neutrino slice true!"<<std::endl;      
      return true;
    }
}

void SingleTrackFilter::clearEvent()
{ 
  pfparticles.clear();
  pfneutrinos.clear();
  pfdaughters.clear();
  pftracks.clear();
  particleMap.clear();
  particlesToMetadata.clear();
  particlesToVertices.clear();
  particlesToTracks.clear();
  pfparticles.clear();
  pfshowers.clear();
  particlesToClusters.clear();
  particlesToSpacePoints.clear();
  particlesToShowers.clear();
  particlesToTracks.clear();
  clustersToHits.clear();
  hitsToSpacePoints.clear();
  spacePointsToHits.clear();
  tracksToHits.clear();
 _deltaY =-9999;
      _deltaZ =-9999;
      _deltaYSigma =-9999;
      _deltaZSigma =-9999;
      _xclVariable =-9999;
      _flash_score = 9999;
      _p_startx =-9999;
      _p_starty =-9999;
      _p_startz =-9999;
      _p_endx =-9999;
      _p_endy =-9999;
      _p_endz =-9999;
      _p_length =-9999;
      _p_chi2_p_2=9999;
}

void SingleTrackFilter::endJob(){
  // Implementation of optional member function here.
}
SingleTrackFilter::DepositionVector SingleTrackFilter::GetDepositionVector(const lar_pandora::PFParticleMap &particleMap, const lar_pandora::PFParticlesToSpacePoints &particlesToSpacePoints,const lar_pandora::SpacePointsToHits &spacePointsToHits, const art::Ptr<recob::PFParticle> &particle)
{
  DepositionVector depositionVector;
  // Get the associated spacepoints
  const auto &partToSpacePointIter(particlesToSpacePoints.find(particle));
  if (partToSpacePointIter == particlesToSpacePoints.end())
    std::cout<<"[Single Proton] particle not found"<<std::endl;
  else
    {
      for (const auto &spacePoint : partToSpacePointIter->second)
        {
          // Get the associated hit
          const auto &spacePointToHitIter(spacePointsToHits.find(spacePoint));
          if (spacePointToHitIter == spacePointsToHits.end())
            continue;

          // Only use hits from the collection plane
          const auto &hit(spacePointToHitIter->second);
          if (hit->View() != geo::kZ)
            continue;

          // Add the charged point to the vector
          const auto &position(spacePoint->XYZ());
          const auto charge(hit->Integral());
          const float nphotons = charge * 164.;//164 comes from flash_neutrino_id.fcl 
          depositionVector.emplace_back(position[0], position[1], position[2], charge, nphotons);
          //      if(_debug) std::cout<<"[depositionvector]"<<position[0]<<" "<<charge<<" "<<nphotons<<std::endl;
        }
    }
  return depositionVector;
}
SingleTrackFilter::Deposition::Deposition(const float x, const float y, const float z, const float charge, const float nPhotons) : m_x(x),
                                                                                                                                                 m_y(y),
                                                                                                                                                 m_z(z),
                                                                                                                                                 m_charge(charge),
                                                                                                                                                 m_nPhotons(nPhotons)
{
}
pandora::CartesianVector SingleTrackFilter::GetChargeWeightedCenter(const DepositionVector &depositionVector)
{
    pandora::CartesianVector center(0.f, 0.f, 0.f);
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
    {
        center += pandora::CartesianVector(chargePoint.m_x, chargePoint.m_y, chargePoint.m_z) * chargePoint.m_charge;
        totalCharge += chargePoint.m_charge;
    }

    if (totalCharge <= std::numeric_limits<float>::epsilon())
        throw cet::exception("FlashNeutrinoId") << "Can't find charge weighted center of slice with zero total charge" << std::endl;

    center *= (1.f / totalCharge);

    return center;
}
flashana::QCluster_t SingleTrackFilter::GetLightCluster(const DepositionVector &depositionVector)
{
    flashana::QCluster_t lightCluster;

    for (const auto &point : depositionVector)
        lightCluster.emplace_back(point.m_x, point.m_y, point.m_z, point.m_nPhotons);

    return lightCluster;
}
float SingleTrackFilter::GetTotalCharge(const DepositionVector &depositionVector)
{
    float totalCharge(0.f);

    for (const auto &chargePoint : depositionVector)
        totalCharge += chargePoint.m_charge;

    return totalCharge;
}



DEFINE_ART_MODULE(SingleTrackFilter)

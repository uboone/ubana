

#ifndef FILLTREEVARIABLES_H
#define FILLTREEVARIABLES_H

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "TTree.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "RecoMCMatching.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "ParticleAssociations.h"
#include "EnergyHelper.h"


class FillTreeVariables {

  std::string fmcordata;  
  std::string ftrack_producer;
  std::string fswtrigger_product;
  std::string fshower_producer;
  std::string fhit_producer;
  std::string frmcmassociation_producer;
  std::string fopflash_producer;
  RecoMCMatching const * frmcm;

  geoalgo::AABox ftpc_volume;
  double foffset;
  geoalgo::AABox ffiducial_volume;
  geoalgo::GeoAlgo const falgo;

  int const fwire_plane;

  bool fverbose;

  lee::EnergyHelper energyHelper;

  TTree * fevent_tree;
  TTree * fvertex_tree;

  int run_number;
  int subrun_number;
  int event_number;
  int nu_pdg;
  double nu_energy;
  int lep_pdg;
  double lep_energy;
  int ccnc;
  int mode;
  int interaction_type;

  int exiting_photon_number;
  int exiting_proton_number;
  int exiting_neutron_number;
  int exiting_electron_number;
  int exiting_antielectron_number;
  int exiting_muon_number;
  int exiting_antimuon_number;
  int exiting_piplus_number;
  int exiting_piminus_number;
  int exiting_pi0_number;
  int total_exiting_particles;
  std::vector<int> exiting_particle_vector;

  int is_single_photon;
  int is_delta_rad;
  int delta_true_pdg;
  double delta_true_energy;
  double delta_photon_energy;
  double delta_proton_energy;
  int delta_mcshower_true_pdg;
  double delta_mcshower_true_energy;
  double delta_mcshower_detprofile_energy;
  int delta_mctrack_true_pdg;
  double delta_mctrack_true_energy;
  double delta_mctrack_true_length;
  
  int mctracknumber;
  int mcshowernumber;
  int tracknumber;
  int showernumber;

  int external_singlephoton_protons_neutrons_only;

  int passed_swtrigger;
  double totalpe_sum;  
  double totalpe_ibg_sum;
  double totalpe_bbg_sum;
  double shower_dist_to_flashzcenter;
  double shower_dist_to_flashzcenter_splitshowers;  

  double most_energetic_reco_shower;
  double second_most_energetic_reco_shower;

  double true_nuvertx;
  double true_nuverty;
  double true_nuvertz;

  double true_nu_E;

  int true_nu_vtx_tpc_contained;
  int true_nu_vtx_fid_contained;

  int selected;
  int number_of_selected_vertices;
  int number_of_selected_vertices_cut;

  double reco_nuvertx;
  double reco_nuverty;
  double reco_nuvertz;

  double reco_true_nuvert_dist;

  double reco_nu_vtx_dist_to_closest_tpc_wall;
  int reco_nu_vtx_fid_contained;

  int reco_asso_tracks;
  int reco_asso_showers;

  double longest_asso_track_matching_ratio;
  int longest_asso_track_matched_to_mcshower;
  int longest_asso_track_matched_to_mctrack;
  int longest_asso_track_matched_to_mcparticle;
  int longest_asso_track_true_pdg;
  int longest_asso_track_true_parent_pdg;
  int longest_asso_track_true_ancestor_pdg;
  int longest_asso_track_true_origin;
  double longest_asso_track_true_startx; 
  double longest_asso_track_true_starty; 
  double longest_asso_track_true_startz;    
  double longest_asso_track_true_endx; 
  double longest_asso_track_true_endy; 
  double longest_asso_track_true_endz;    
  int longest_asso_track_from_ncdeltarad_interaction;
  int longest_asso_track_matched_to_ncdeltarad_photon;
  int longest_asso_track_matched_to_ncdeltarad_proton;
  double longest_asso_track_true_thetayx;
  double longest_asso_track_true_thetaxz;
  double longest_asso_track_true_thetayz;
  double longest_asso_track_true_energy;

  double longest_asso_track_length;
  double longest_asso_track_reco_dirx;
  double longest_asso_track_reco_diry;
  double longest_asso_track_reco_dirz;
  double longest_asso_track_thetayx;
  double longest_asso_track_thetaxz;
  double longest_asso_track_thetayz;
  double longest_asso_track_phi;
  double longest_asso_track_theta;

  std::vector<double> longest_asso_track_calo_dEdx;
  std::vector<double> longest_asso_track_calo_resrange;

  double shower_matching_ratio;
  int shower_matched_to_mcshower;
  int shower_matched_to_mctrack;
  int shower_matched_to_mcparticle;
  int shower_true_pdg;
  int shower_true_parent_pdg;
  int shower_true_ancestor_pdg;
  int shower_true_origin;
  double shower_true_startx; 
  double shower_true_starty; 
  double shower_true_startz;    
  double shower_true_endx; 
  double shower_true_endy; 
  double shower_true_endz;
  int shower_from_ncdeltarad_interaction;
  int shower_matched_to_ncdeltarad_photon;
  int shower_matched_to_ncdeltarad_proton;
  double shower_true_dist; 
  double shower_true_distx;
  double shower_true_disty;
  double shower_true_distz;
  double shower_true_thetayx;
  double shower_true_thetaxz;
  double shower_true_thetayz;
  double shower_true_energy;
  double shower_detprofile_thetayx;
  double shower_detprofile_thetaxz;
  double shower_detprofile_thetayz;
  double shower_detprofile_energy;
 
  double closest_asso_shower_dist_to_flashzcenter;

  double most_energetic_shower_reco_startx;
  double most_energetic_shower_reco_starty;
  double most_energetic_shower_reco_startz;
  double most_energetic_shower_reco_dist;
  double most_energetic_shower_reco_distx;
  double most_energetic_shower_reco_disty;
  double most_energetic_shower_reco_distz;
  double most_energetic_shower_reco_thetayx;
  double most_energetic_shower_reco_thetaxz;
  double most_energetic_shower_reco_thetayz;
  double most_energetic_shower_reco_width0;
  double most_energetic_shower_reco_width1;
  double most_energetic_shower_reco_opening_angle;
  double most_energetic_shower_reco_length;
  double most_energetic_shower_reco_dirx;
  double most_energetic_shower_reco_diry;
  double most_energetic_shower_reco_dirz;
  double most_energetic_shower_reco_energy;
  double most_energetic_shower_helper_energy;
  std::vector<double> reco_shower_energy_vector;
  double most_energetic_shower_bp_dist_to_tpc;
  std::vector<double> reco_shower_dedx_vector;
  double reco_shower_dedx_plane0;
  double reco_shower_dedx_plane1;
  double reco_shower_dedx_plane2;
  double reco_shower_dedx_best_plane;
  double closest_shower_dedx_plane2;
  double closest_shower_dedx_best_plane;

  double shortest_asso_shower_to_vert_dist;
  double summed_associated_reco_shower_energy;
  double summed_associated_helper_shower_energy;
  double summed_associated_helper_track_energy;

  int pass_all_cuts;

  std::map<std::string, std::vector<double *>> weight_branch_map; 
  //EventWeight
  std::map<std::string, std::vector<double>> fmcweight;


  double fweight_genie_ncelaxial_p1sigma;
  double fweight_genie_ncelaxial_m1sigma;
  double fweight_genie_nceleta_p1sigma;
  double fweight_genie_nceleta_m1sigma;
  double fweight_genie_qema_p1sigma;
  double fweight_genie_qema_m1sigma;
  double fweight_genie_qevec_p1sigma;
  double fweight_genie_qevec_m1sigma;
  double fweight_genie_ccresaxial_p1sigma;
  double fweight_genie_ccresaxial_m1sigma;
  double fweight_genie_ccresvector_p1sigma;
  double fweight_genie_ccresvector_m1sigma;
  double fweight_genie_resganged_p1sigma;
  double fweight_genie_resganged_m1sigma;
  double fweight_genie_ncresaxial_p1sigma;
  double fweight_genie_ncresaxial_m1sigma;
  double fweight_genie_ncresvector_p1sigma;
  double fweight_genie_ncresvector_m1sigma;
  double fweight_genie_cohma_p1sigma;
  double fweight_genie_cohma_m1sigma;
  double fweight_genie_cohr0_p1sigma;
  double fweight_genie_cohr0_m1sigma;
  double fweight_genie_nonresrvp1pi_p1sigma;
  double fweight_genie_nonresrvp1pi_m1sigma;
  double fweight_genie_nonresrvbarp1pi_p1sigma;
  double fweight_genie_nonresrvbarp1pi_m1sigma;
  double fweight_genie_nonresrvp2pi_p1sigma;
  double fweight_genie_nonresrvp2pi_m1sigma;
  double fweight_genie_nonresrvbarp2pi_p1sigma;
  double fweight_genie_nonresrvbarp2pi_m1sigma;
  double fweight_genie_resdecaygamma_p1sigma;
  double fweight_genie_resdecaygamma_m1sigma;
  double fweight_genie_resdecayeta_p1sigma;
  double fweight_genie_resdecayeta_m1sigma;
  double fweight_genie_resdecaytheta_p1sigma;
  double fweight_genie_resdecaytheta_m1sigma;
  double fweight_genie_nc_p1sigma;
  double fweight_genie_nc_m1sigma;
  double fweight_genie_disath_p1sigma;
  double fweight_genie_disath_m1sigma;
  double fweight_genie_disbth_p1sigma;
  double fweight_genie_disbth_m1sigma;
  double fweight_genie_discv1u_p1sigma;
  double fweight_genie_discv1u_m1sigma;
  double fweight_genie_discv2u_p1sigma;
  double fweight_genie_discv2u_m1sigma;
  double fweight_genie_disnucl_p1sigma;
  double fweight_genie_disnucl_m1sigma;
  double fweight_genie_agkyxf_p1sigma;
  double fweight_genie_agkyxf_m1sigma;
  double fweight_genie_agkypt_p1sigma;
  double fweight_genie_agkypt_m1sigma;
  double fweight_genie_formzone_p1sigma;
  double fweight_genie_formzone_m1sigma;
  double fweight_genie_fermigasmodelkf_p1sigma;
  double fweight_genie_fermigasmodelkf_m1sigma;
  double fweight_genie_fermigasmodelsf_p1sigma;
  double fweight_genie_fermigasmodelsf_m1sigma;
  double fweight_genie_intranukenmfp_p1sigma;
  double fweight_genie_intranukenmfp_m1sigma;
  double fweight_genie_intranukencex_p1sigma;
  double fweight_genie_intranukencex_m1sigma;
  double fweight_genie_intranukenel_p1sigma;
  double fweight_genie_intranukenel_m1sigma;
  double fweight_genie_intranukeninel_p1sigma;
  double fweight_genie_intranukeninel_m1sigma;
  double fweight_genie_intranukenabs_p1sigma;
  double fweight_genie_intranukenabs_m1sigma;
  double fweight_genie_intranukenpi_p1sigma;
  double fweight_genie_intranukenpi_m1sigma;
  double fweight_genie_intranukepimfp_p1sigma;
  double fweight_genie_intranukepimfp_m1sigma;
  double fweight_genie_intranukepicex_p1sigma;
  double fweight_genie_intranukepicex_m1sigma;
  double fweight_genie_intranukepiel_p1sigma;
  double fweight_genie_intranukepiel_m1sigma;
  double fweight_genie_intranukepiinel_p1sigma;
  double fweight_genie_intranukepiinel_m1sigma;
  double fweight_genie_intranukepiabs_p1sigma;
  double fweight_genie_intranukepiabs_m1sigma;
  double fweight_genie_intranukepipi_p1sigma;
  double fweight_genie_intranukepipi_m1sigma;

  void ResetEvent();
  void ResetVertex();
  void FillWeights(art::Event const & e);
  void FillEventWeights(art::Event const & e);

public:

  FillTreeVariables();

  void SetVerbose(bool const verbose = true) {fverbose = verbose;}
  void SetProducers(std::string const & mcordata,
		    std::string const & track_producer,
		    std::string const & shower_producer,
		    std::string const & hit_producer,
		    std::string const & opflash_producer,
		    std::string const & trigger_product,
		    std::string const & rmcmassociations_producer,
		    RecoMCMatching const * rmcm);
  void SetUpTreeBranches();
  bool SinglePhotonFilter(art::Event const & e,
			  size_t & mctruth_index);
  bool DeltaRadFilter(art::Event const & e,
		      size_t & mctruth_index);
  void FillTruth(art::Event const & e,
		 size_t & delta_rad_mct_index);
  bool PassedSWTrigger(art::Event const & e) const;
  void GetExitingDeltaPhotonAndNucleon(art::Event const & e,
				       size_t const delta_rad_mct_index,
				       size_t & delta_photon_index,
				       size_t & delta_proton_index);
  void GetDeltaMCShowerMCTrackIndices(art::Event const & e,
				      size_t const delta_rad_mct_index,
				      size_t & delta_photon_index,
				      size_t & delta_mcshower_index,
				      size_t & delta_proton_index,
				      size_t & delta_mctrack_index);
  double DistToClosestTPCWall(geoalgo::Point_t const & pos);
  double GetShowerHelperEnergy(art::Event const & e,
			       size_t const shower_index);
  double GetTrackHelperEnergy(art::Event const & e,
			      size_t const track_index);

  std::pair<std::vector<double>,std::vector<double>> GetTrackCaloInfo(art::Event const & e,    size_t const track_index);


  int GetBestShowerPlane(art::Event const & e,
			 size_t const shower_index);
  void FilldEdx(art::Event const & e,
		size_t const most_energetic_associated_shower_index,
		size_t const closest_associated_shower_index);
  double ShowerZDistToClosestFlash(art::Event const & e,
				   int const most_energetic_shower_index);
  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const & e, int const geant_track_id);
  void FillShowerRecoMCMatching(art::Event const & e,
				size_t const most_energetic_associated_shower_index,
				size_t const delta_rad_mct_index,
				size_t const delta_mcshower_index,
				size_t const delta_mctrack_index);
  void FillTrackRecoMCMatching(art::Event const & e,
			       size_t const longest_asso_track_index,
			       size_t const delta_rad_mct_index,
			       size_t const delta_mcshower_index,
			       size_t const delta_mctrack_index);
  void FillVertexTree(art::Event const & e,
		      ParticleAssociations const & pas,
		      size_t const pn,
		      size_t const delta_rad_mct_index,
		      size_t const delta_mcshower_index,
		      size_t const delta_mctrack_index);
  void Fill(art::Event const & e,
	    ParticleAssociations const & pas);

};


#endif

#ifndef SINGLE_PHOTON_ANALYSIS
#define SINGLE_PHOTON_ANALYSIS

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/simb.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Helper function for PID stuff
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"

#include "Pandora/PdgTable.h"
#include <chrono>

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>

#include "bad_channel_matching.h"
//------------------------------------------------------------------------------------------------------------------------------------------

namespace single_photon
{

    template <typename T>
        std::vector<size_t> sort_indexes(const std::vector<T> &v) {

            std::vector<size_t> idx(v.size());
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in v
            std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

            return idx;
        }


    /**
     *  @brief  SinglePhoton class
     */
    class SinglePhoton : public art::EDAnalyzer
    {
        public:
            typedef art::ValidHandle< std::vector<recob::PFParticle> > PFParticleHandle;
            typedef std::vector< art::Ptr<recob::PFParticle> > PFParticleVector;
            typedef std::vector< art::Ptr<recob::Track> > TrackVector;
            typedef std::vector< art::Ptr<recob::Shower> > ShowerVector;
            typedef std::map< size_t, art::Ptr<recob::PFParticle>> PFParticleIdMap;

            /**
             *  @brief  Constructor
             *
             *  @param  pset the set of input fhicl parameters
             */
            SinglePhoton(fhicl::ParameterSet const &pset);

            /**
             *  @brief  Configure memeber variables using FHiCL parameters
             *
             *  @param  pset the set of input fhicl parameters
             */
            void reconfigure(fhicl::ParameterSet const &pset);

            /**
             *  @brief  Analyze an event!
             *
             *  @param  evt the art event to analyze
             */
            void analyze(const art::Event &evt);

            /**
             *  @brief  Begin the job, setting up !
             *
             */
            void beginJob();

            /**
             *  @brief  End the job, setting down !
             *
             */
            void endJob();


            void beginSubRun(art::SubRun const & sr);

        private:
            void ClearVertex();
            /**
             *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
             *
             *  @param  pfParticleHandle the handle for the PFParticle collection
             *  @param  pfParticleMap the mapping from ID to PFParticle
             */
            void GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap);

            /**
             * @brief Print out scores in PFParticleMetadata
             *
             * @param evt the art event to analyze
             * @param pfParticleHandle the handle for the PFParticle collection
             */
            void PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const;

            /**
             *  @brief  Produce a mapping from PFParticle ID to the art ptr to the PFParticle itself for fast navigation
             *
             *  @param  pfParticleMap the mapping from ID to PFParticle
             *  @param  crParticles a vector to hold the top-level PFParticles reconstructed under the cosmic hypothesis
             *  @param  nuParticles a vector to hold the final-states of the reconstruced neutrino
             */
            void GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap,const lar_pandora::PFParticlesToVertices &particlesToVertices, PFParticleVector &crParticles, PFParticleVector &nuParticles);

            /**
             *  @brief  Collect associated tracks and showers to particles in an input particle vector
             *
             *  @param  particles a vector holding PFParticles from which to find the associated tracks and showers
             *  @param  pfParticleHandle the handle for the PFParticle collection
             *  @param  evt the art event to analyze
             *  @param  tracks a vector to hold the associated tracks
             *  @param  showers a vector to hold the associated showers
             */
            void CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap );

            void GetVertex(const lar_pandora::PFParticlesToVertices & particlestoVertices, const art::Ptr<recob::PFParticle> & particle );

            void CollectCalo(const art::Event &evt,const art::Ptr<recob::Shower> &shower);


            /*
             *@brief Calculated the shower energy by looping over all the hits and summing the charge
             *@param hits -  an art pointer of all the hits in a shower
             *
             *
             * */
            double CalcEShower(std::vector<art::Ptr<recob::Hit>> hits);

            /**
             *@brief Takes a hit and multiplies the charge by the gain
             *@param thishitptr art pointer to a hit
             *@param plane the plane the hit is on
             **/
            double GetQHit(art::Ptr<recob::Hit> thishitptr, int plane);

            /**
             * @brief Calculate the E value in MeV for a given hit
             * @param thishit - an individual hit 
             * 
             *
             * */
            double QtoEConversionHit(art::Ptr<recob::Hit> thishitptr, int plane);

            /**
             * @brief Calculate the E value in MeV from a given Q value
             * @param q - the charge value
             * 
             * */
            double QtoEConversion(double q);


            /**
             *@brief Takes a vector of dQ/dx values and converts to dE/dx
             *@param dqdx - vector of dqdx points
             *
             * */
            std::vector<double> CalcdEdxFromdQdx(std::vector<double> dqdx);

            /**
             *
             *@brief For a single shower, calculates the dQdx for each hit in the clusters in the shower for a single plane
             *@param shower - a Pandora shower
             *@param clusters - all of the clusters in the shower
             *@param clusterToHitMap - a map between each cluster and all of the hits in the cluster
             *@param plane - a single plane
             * * */

            std::vector<double> CalcdQdxShower(
                    const art::Ptr<recob::Shower>& shower,
                    const std::vector<art::Ptr<recob::Cluster>> & clusters, 
                    std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane);
            /**
             *@brief Gets the pitch between the 3D reconstructed shower direction and the wires for a given plane (the dx in dQdx)
             *@param shower_dir - the 3D shower direction
             *@param plane - a single plane
             * */
            double getPitch(TVector3 shower_dir, int plane);

            /**
             *@brief Calculates the four corners of a box of given length and width around a cluster given the start point and axis direction
             *@param cluster_start - the start position of a cluster in CM
             *@param cluster_axis - calculated from the cluster end minus the cluster start
             *@param width - typically ~1cm
             *@param length - typically a few cm
             *
             * */
            std::vector<std::vector<double>> buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length);

            /**
             *@brief For a 2d point on a plane in cm and a rectangle, returns true if the point is inside of the rectangle
             *@param thishit_pos - 2d location of a hit in cm
             *@param rectangle - vector of the positions of the four corners of the rectangle
             *
             * */
            bool insideBox(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);

            /**
             *
             *@brief For a 2d point on a plane in cm and a rectangle, returns true if ponint is inside or on the boundary
             *uses triangle area check
             *
             * */
            bool isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);

            double areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3);

            /***
             *@brief returns the value at the median position in a vector of doubles, returns nan for vector of size <= 0
             *@param thisvector - vector of doubles
             *
             * */
            double getMedian(std::vector<double> thisvector);


            //----------------  Templatees ----------------------------
            void AnalyzeTemplates();
            void ClearTemplates();
            void ResizeTemplates(size_t);
            void CreateTemplateBranches();



            //----------------  Flashes ----------------------------
            void AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes);
            void ClearFlashes();
            void ResizeFlashes(size_t);
            void CreateFlashBranches();

            //----------------  Tracks ----------------------------
            void AnalyzeTracks(const std::vector<art::Ptr<recob::Track>>& tracks, std::map<art::Ptr<recob::Track>, art::Ptr<recob::PFParticle>> & tracktopfparticlemap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfparticletospacepointmap , std::map<int, art::Ptr<simb::MCParticle> > &  MCParticleToTrackIdMap);
            void ClearTracks();
            void ResizeTracks(size_t);
            void CreateTrackBranches();
            void AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::map<art::Ptr<recob::Track>,art::Ptr<anab::Calorimetry>> &trackToCaloMap);
            void RecoMCTracks(const std::vector<art::Ptr<recob::Track>>& tracks,  std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle>> & trackToPFParticleMap, std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > & trackToMCParticleMap,  std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector, std::map< int, art::Ptr<simb::MCParticle> > &      MCParticleToTrackIdMap);

            void CollectPID(std::vector<art::Ptr<recob::Track>> & tracks,std::map< art::Ptr<recob::Track>, art::Ptr<anab::ParticleID>> & trackToPIDMap);
            TGraph proton_length2energy_tgraph;

            //----------------  Showers ----------------------------

            void AnalyzeShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>>> & pfParticleToHitMap,std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > & pfParticleToClusterMap, std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > & clusterToHitMap );
            void ClearShowers();
            void ResizeShowers(size_t);
            void CreateShowerBranches();

            void RecoMCShowers(const std::vector<art::Ptr<recob::Shower>>& showers,  std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>> & showerToPFParticleMap, std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap,  std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>> & MCParticleToMCTruthMap,
                    std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector);

            std::vector<double> showerRecoMCmatching(std::vector<art::Ptr<recob::Shower>>& objectVector,
            std::map<art::Ptr<recob::Shower>,art::Ptr<simb::MCParticle>>& objectToMCParticleMap,
            std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle>>& objectToPFParticleMap,
            std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> >& pfParticleToHitsMap,
            art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData>& mcparticles_per_hit,
            std::vector<art::Ptr<simb::MCParticle>>& mcParticleVector,
            std::map< size_t, art::Ptr<recob::PFParticle>> & pfParticleIdMap,
            std::map< int ,art::Ptr<simb::MCParticle> >  &  MCParticleToTrackIdMap );




            //---------------- MCTruths ----------------------------

            void AnalyzeMCTruths(std::vector<art::Ptr<simb::MCTruth>> & mcTruthVector,  std::vector<art::Ptr<simb::MCParticle>> & mcParticleVector );
            void ClearMCTruths();
            void ResizeMCTruths(size_t);
            void CreateMCTruthBranches();

            std::map<int,std::string> is_delta_map;


            //These three are shameless steals from LArPandorHelper But overlays dont work so this is a direct clone. We will filter out later.
            void CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector);
            void CollectMCParticles(const art::Event &evt, const std::string &label, std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>              &particlesToTruth, std::map< int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIdMap);
            void BuildMCParticleHitMaps(const art::Event &evt, const std::string &label, const std::vector<art::Ptr<recob::Hit>> &hitVector,   std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap);


            //-------------- Slices/Pandora Metadata ---------------//
            void  ClearSlices();
            void  ResizeSlices(size_t size); 
            void CreateSliceBranches();
            void AnalyzeSlices( std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > & pfParticleToMetadataMap,  PFParticleIdMap &pfParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec);
            int GetShowerSlice(art::Ptr<recob::Shower>& this_shower, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>>& showerToPFParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec);
            int GetTrackSlice(art::Ptr<recob::Track>& this_track, std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>& trackToPFParticleMap, std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec);
            //can also look at things like shower energy, conversion length, etc.
            void FindSignalSlice(std::string signal_def, std::map<int, art::Ptr<simb::MCParticle>> & MCParticleToTrackIDMap,std::map<art::Ptr<recob::Shower>,art::Ptr<recob::PFParticle> > & showerToPFParticleMap,  std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> & allPFPSliceIdVec, std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > & showerToMCParticleMap, std::map<art::Ptr<recob::Track>,art::Ptr<recob::PFParticle> > & trackToNuPFParticleMap, std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > &trackToMCParticleMap);
            
            int  m_reco_slice_num; //total number of slices in the event
            std::vector<double> m_reco_slice_nuscore; //vector of the neutrino score for each slice in an event
            int m_sim_shower_num_matched_signal; //the number of sim showers matched an MCP in the signal def
            int m_sim_track_num_matched_signal; //the number of sim showers matched an MCP in the signal def
            //std::map<art::Ptr<recob::PFParticle>, double > & pfParticleToNuScoreMap;//is filled during analyze slices

            //------------------ Delaunay triangle tools -----------//

            double triangle_area(double a1, double a2, double b1, double b2, double c1, double c2);
            int quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area);
            int delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area);


            int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected);
            int spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected);
            int spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input);

            //databased http://dbdata0vm.fnal.gov:8186/uboonecon_prod/app/data?f=channelstatus_data&t=357812824
            std::vector<std::pair<int,int>> bad_channel_list_fixed_mcc9;
            std::map<int,bool> bad_channel_map_fixed_mcc9;

            std::string m_pandoraLabel;         ///< The label for the pandora producer
            std::string m_trackLabel;           ///< The label for the track producer from PFParticles
            std::string m_showerLabel;          ///< The label for the shower producer from PFParticles
            std::string m_caloLabel;            ///< The label for calorimetry associations producer
            std::string m_flashLabel;
            std::string m_geantModuleLabel;
            std::string m_backtrackerLabel;
            std::string m_hitfinderLabel;
            std::string m_hitMCParticleAssnsLabel;
            std::string m_potLabel;
            std::string m_generatorLabel;
            std::string m_badChannelLabel;
            std::string m_badChannelProducer;
            std::string m_mcTrackLabel;
            std::string m_mcShowerLabel;
            std::string m_pidLabel;            ///< For PID stuff
            bool m_use_PID_algorithms;
            bool m_use_delaunay;
            bool m_is_verbose;
            bool m_is_data;
            bool m_is_overlayed;

            double m_track_calo_min_dEdx;
            double m_track_calo_max_dEdx;
            double m_track_calo_min_dEdx_hits;
            double m_track_calo_trunc_fraction;

            detinfo::DetectorProperties const * theDetector ;// = lar::providerFrom<detinfo::DetectorPropertiesService>();
            detinfo::DetectorClocks    const *  detClocks   ;//= lar::providerFrom<detinfo::DetectorClocksService>();
            spacecharge::SpaceCharge const * SCE;
            geo::GeometryCore const * geom;
            double m_work_function;
            double m_recombination_factor;
            //double m_gain;
            std::vector<double> m_gain_mc; 
            std::vector<double> m_gain_data; 
            double m_wire_spacing;

            int m_Cryostat;
            int m_TPC;

            double m_width_dqdx_box;
            double m_length_dqdx_box;

            TTree* pot_tree;
            TTree* vertex_tree;

            //------------ POT related variables --------------
            int m_number_of_events;
            double m_pot_count;
            int m_number_of_vertices;

            //------------ Event Related Variables -------------
            int m_run_number;
            int m_subrun_number;
            int m_event_number;

            //------------ Vertex Related variables -------------
            int m_reco_vertex_size;
            double m_vertex_pos_x;
            double m_vertex_pos_y;
            double m_vertex_pos_z;
            int m_reco_asso_showers;

            double m_reco_vertex_to_nearest_dead_wire_plane0;
            double m_reco_vertex_to_nearest_dead_wire_plane1;
            double m_reco_vertex_to_nearest_dead_wire_plane2;


            //-------------- Flash related variables -------------
            int m_reco_num_templates;
            std::vector<double> m_reco_template;


            //-------------- Flash related variables -------------
            std::vector<double> m_reco_flash_total_pe;
            std::vector<double> m_reco_flash_time;
            std::vector<double> m_reco_flash_time_width;
            std::vector<double> m_reco_flash_abs_time;
            std::vector<int>    m_reco_flash_frame;
            std::vector<double> m_reco_flash_ycenter;
            std::vector<double> m_reco_flash_ywidth;
            std::vector<double> m_reco_flash_zcenter;
            std::vector<double> m_reco_flash_zwidth;
            std::vector<double> m_reco_flash_total_pe_in_beamgate;
            std::vector<double> m_reco_flash_time_in_beamgate;
            std::vector<double> m_reco_flash_ycenter_in_beamgate;
            std::vector<double> m_reco_flash_zcenter_in_beamgate;

            int m_reco_num_flashes;
            int m_reco_num_flashes_in_beamgate;

            double m_beamgate_flash_start;
            double m_beamgate_flash_end;

            //------------ Track related Variables -------------
            int m_reco_asso_tracks;
            std::vector<double> m_reco_track_length;
            std::vector<double> m_reco_track_dirx;
            std::vector<double> m_reco_track_diry;
            std::vector<double> m_reco_track_dirz;
            std::vector<double> m_reco_track_startx;
            std::vector<double> m_reco_track_starty;
            std::vector<double> m_reco_track_startz;
            std::vector<double> m_reco_track_endx;
            std::vector<double> m_reco_track_endy;
            std::vector<double> m_reco_track_endz;
            std::vector<double>   m_reco_track_theta_yz;
            std::vector<double>   m_reco_track_phi_yx;


            std::vector<int> m_reco_track_num_trajpoints;
            std::vector<int> m_reco_track_num_spacepoints;
            std::vector<double> m_reco_track_proton_kinetic_energy;
            std::vector<size_t>  m_reco_track_ordered_energy_index;
            std::vector<double> m_reco_track_spacepoint_principal0;
            std::vector<double> m_reco_track_spacepoint_principal1;
            std::vector<double> m_reco_track_spacepoint_principal2;

            std::vector<double> m_reco_track_spacepoint_chi;
            std::vector<double> m_reco_track_spacepoint_max_dist;

            std::vector<double> m_reco_track_mean_dEdx;
            std::vector<double> m_reco_track_mean_dEdx_start_half;
            std::vector<double> m_reco_track_mean_dEdx_end_half;
            std::vector<int> m_reco_track_good_calo;
            std::vector<double> m_reco_track_mean_trunc_dEdx;
            std::vector<double> m_reco_track_mean_trunc_dEdx_start_half;
            std::vector<double> m_reco_track_mean_trunc_dEdx_end_half;
            std::vector<double> m_reco_track_trunc_PIDA;
            std::vector<std::vector<double>> m_reco_track_resrange;
            std::vector<std::vector<double>> m_reco_track_dEdx;


            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane0;
            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane1;
            std::vector<double> m_reco_track_end_to_nearest_dead_wire_plane2;

            std::vector<int> m_sim_track_matched;
            std::vector<double> m_sim_track_energy;
            std::vector<double> m_sim_track_mass;
            std::vector<double> m_sim_track_kinetic_energy;
            std::vector<int> m_sim_track_pdg;
            std::vector<int> m_sim_track_parent_pdg;
            std::vector<int> m_sim_track_origin;
            std::vector<std::string> m_sim_track_process;
            std::vector<double> m_sim_track_startx;
            std::vector<double> m_sim_track_starty;
            std::vector<double> m_sim_track_startz;
            std::vector<int> m_sim_track_trackID;


            //------------ Shower related Variables  -------------

            std::vector<double>   m_reco_shower_startx;
            std::vector<double>   m_reco_shower_starty;
            std::vector<double>   m_reco_shower_startz;
            std::vector<double>   m_reco_shower_dirx;
            std::vector<double>   m_reco_shower_diry;
            std::vector<double>   m_reco_shower_dirz;
            std::vector<double>   m_reco_shower_theta_yz;
            std::vector<double>   m_reco_shower_phi_yx;

            std::vector<double> m_reco_shower_openingangle;
            std::vector<double> m_reco_shower_length;
            std::vector<double> m_reco_shower_conversion_distance;

            std::vector<int> m_reco_shower_delaunay_num_triangles_plane0;
            std::vector<int> m_reco_shower_delaunay_num_triangles_plane1;
            std::vector<int> m_reco_shower_delaunay_num_triangles_plane2;

            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane0;
            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane1;
            std::vector<double> m_reco_shower_start_to_nearest_dead_wire_plane2;


            //shower flash matching

            std::vector<double> m_reco_shower_flash_shortest_distz;
            std::vector<double> m_reco_shower_flash_shortest_disty;
            std::vector<double> m_reco_shower_flash_shortest_distyz;

            std::vector<int> m_reco_shower_flash_shortest_index_z;
            std::vector<int> m_reco_shower_flash_shortest_index_y;
            std::vector<int> m_reco_shower_flash_shortest_index_yz;

            //end flash matching
            std::vector<int> m_reco_shower_num_hits_plane0;
            std::vector<int> m_reco_shower_num_hits_plane1;
            std::vector<int> m_reco_shower_num_hits_plane2;


            std::vector<double> m_reco_shower_delaunay_area_plane0;
            std::vector<double> m_reco_shower_delaunay_area_plane1;
            std::vector<double> m_reco_shower_delaunay_area_plane2;

            std::vector<int> m_sim_shower_matched;
            std::vector<double> m_sim_shower_energy;
            std::vector<double> m_sim_shower_kinetic_energy;
            std::vector<double> m_sim_shower_mass;
            std::vector<int> m_sim_shower_pdg;
            std::vector<int> m_sim_shower_trackID;
            std::vector<int> m_sim_shower_parent_pdg;
            std::vector<int> m_sim_shower_parent_trackID;
            std::vector<int> m_sim_shower_origin;
            std::vector<std::string> m_sim_shower_process;
            std::vector<std::string> m_sim_shower_end_process;
            std::vector<double> m_sim_shower_start_x;
            std::vector<double> m_sim_shower_start_y;
            std::vector<double> m_sim_shower_start_z;
            std::vector<double> m_sim_shower_vertex_x;
            std::vector<double> m_sim_shower_vertex_y;
            std::vector<double> m_sim_shower_vertex_z;

            std::vector<double> m_sim_shower_px;
            std::vector<double> m_sim_shower_py;
            std::vector<double> m_sim_shower_pz;


            std::vector<int> m_sim_shower_is_true_shower;
            std::vector<int> m_sim_shower_best_matched_plane;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane0;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane1;
            std::vector<double> m_sim_shower_matched_energy_fraction_plane2;
            std::vector<double> m_sim_shower_overlay_fraction;


            //------------ MCTruth related Variables  -------------
            int m_mctruth_num;
            int m_mctruth_origin;
            double m_mctruth_nu_E;
            double m_mctruth_nu_vertex_x;
            double m_mctruth_nu_vertex_y;
            double m_mctruth_nu_vertex_z;
            double m_mctruth_lepton_E;
            int m_mctruth_nu_pdg;
            int m_mctruth_lepton_pdg;
            int m_mctruth_mode ;
            int m_mctruth_interaction_type ;
            int m_mctruth_ccnc;
            double m_mctruth_qsqr;

            int m_mctruth_num_daughter_particles;
            std::vector<int> m_mctruth_daughters_pdg;
            std::vector<double> m_mctruth_daughters_E;

            int     m_mctruth_num_exiting_photons ;
            int      m_mctruth_num_exiting_protons ;
            int    m_mctruth_num_exiting_pi0 ;
            int   m_mctruth_num_exiting_pipm ;
            int   m_mctruth_num_exiting_neutrons; 
            int   m_mctruth_num_exiting_delta0; 
            int   m_mctruth_num_exiting_deltapm; 
            int   m_mctruth_num_exiting_deltapp; 

            int m_mctruth_is_delta_radiative;
            int m_mctruth_delta_radiative_1g1p_or_1g1n;
            double m_mctruth_delta_photon_energy;
            double m_mctruth_delta_proton_energy;
            double m_mctruth_delta_neutron_energy;
            std::vector<int> m_mctruth_exiting_delta0_num_daughters;

            std::vector<int> m_mctruth_exiting_photon_trackID;
            std::vector<int> m_mctruth_exiting_photon_mother_trackID;
            std::vector<int> m_mctruth_exiting_photon_from_delta_decay;
            std::vector<double> m_mctruth_exiting_photon_energy;
            std::vector<int> m_mctruth_exiting_proton_trackID;
            std::vector<int> m_mctruth_exiting_proton_mother_trackID;
            std::vector<int> m_mctruth_exiting_proton_from_delta_decay;
            std::vector<double> m_mctruth_exiting_proton_energy;


            std::vector<double>        m_mctruth_exiting_pi0_E;
            std::vector<double>        m_mctruth_exiting_pi0_px;
            std::vector<double>        m_mctruth_exiting_pi0_py;
            std::vector<double>        m_mctruth_exiting_pi0_pz;

            std::string  m_truthmatching_signaldef;

            //the calo calculated quantities 
            std::vector<double> m_reco_shower_energy; //for each hit in a shower, converts Q->E, and sums
            std::vector<size_t>  m_reco_shower_ordered_energy_index;
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane0; //for each shower, looks at the hits for all clusters in the plane, stores the dQ/dx for each hit 
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane1;
            std::vector<std::vector<double>> m_reco_shower_dQdx_plane2;
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane0; //dE/dx from the calculated dQ/dx for each hit in shower on plane 	
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane1;
            std::vector<std::vector<double>> m_reco_shower_dEdx_plane2;

            std::vector<double> m_reco_shower_dEdx_plane0_mean;
            std::vector<double> m_reco_shower_dEdx_plane1_mean;
            std::vector<double> m_reco_shower_dEdx_plane2_mean;
            std::vector<double> m_reco_shower_dEdx_plane0_max;
            std::vector<double> m_reco_shower_dEdx_plane1_max;
            std::vector<double> m_reco_shower_dEdx_plane2_max;
            std::vector<double> m_reco_shower_dEdx_plane0_min;
            std::vector<double> m_reco_shower_dEdx_plane1_min;
            std::vector<double> m_reco_shower_dEdx_plane2_min;
            std::vector<double> m_reco_shower_dEdx_plane0_median;
            std::vector<double> m_reco_shower_dEdx_plane1_median;
            std::vector<double> m_reco_shower_dEdx_plane2_median;
            std::vector<double> m_reco_shower_dEdx_plane0_nhits;
            std::vector<double> m_reco_shower_dEdx_plane1_nhits;
            std::vector<double> m_reco_shower_dEdx_plane2_nhits;

            double _time2cm;//value modeled from David's shower code

            // PID-related variables
            std::vector<double> m_reco_track_pid_bragg_likelihood_plane2;
            std::vector<double> m_reco_track_pid_pida_plane2;
            std::vector<double> m_reco_track_pid_chi_plane2;



    };

    DEFINE_ART_MODULE(SinglePhoton)

} // namespace lar_pandora
#endif
//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

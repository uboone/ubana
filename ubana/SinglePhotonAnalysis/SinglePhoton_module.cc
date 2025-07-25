#include "SinglePhoton_module.h"
#include "analyze_OpFlashes.h"
#include "analyze_Tracks.h"
#include "analyze_Showers.h"
#include "analyze_Template.h"
#include "analyze_MCTruth.h"
#include "analyze_EventWeight.h"
#include "analyze_Slice.h"
#include "analyze_Geant4.h"
#include "fiducial_volume.h"
#include "second_shower_search.h"
#include "isolation.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()

namespace single_photon
{

    SinglePhoton::SinglePhoton(fhicl::ParameterSet const &pset) : 
      art::EDFilter(pset),
      detClocks(art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob()),
      theDetector(art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(detClocks))
    {
        this->reconfigure(pset);
        SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        geom = lar::providerFrom<geo::Geometry>();
        m_channelMap = &art::ServiceHandle<geo::WireReadout const>()->Get();
    }

    //Reconfigure the internal class parameters from .fcl parameters
    void SinglePhoton::reconfigure(fhicl::ParameterSet const &pset)
    {
        //input parameters for what file/mode were running in
        m_print_out_event = pset.get<bool>("PrintOut", false);
        m_is_verbose = pset.get<bool>("Verbose",false);
        m_is_data = pset.get<bool>("isData",false);
        m_is_overlayed = pset.get<bool>("isOverlayed",false);
        m_is_textgen = pset.get<bool>("isTextGen",false);

        //some specific additonal info, default not include
        m_use_PID_algorithms = pset.get<bool>("usePID",false);
        m_use_delaunay = pset.get<bool>("useDelaunay",false);
        m_delaunay_max_hits = pset.get<int>("maxDelaunayHits",1000);

        //When we moved to filter instead of analyzer, kept ability to run 2g filter. A bit depreciated (not in code, but in our use case)
        m_fill_trees          = pset.get<bool>("FillTrees",true);
        m_run_pi0_filter      = pset.get<bool>("RunPi0Filter",false);
        m_run_pi0_filter_2g1p = pset.get<bool>("FilterMode2g1p",false);
        m_run_pi0_filter_2g0p = pset.get<bool>("FilterMode2g0p",false);

        if(m_run_pi0_filter) m_is_data = true;// If running in filter mode, treat all as data

        //Some output for logging
        std::cout<<"SinglePhoton::reconfigure || whats configured? "<<std::endl;
        std::cout<<"SinglePhoton::reconfigure || m_is_data: "<<m_is_data<<std::endl;
        std::cout<<"SinglePhoton::reconfigure || m_is_overlayed: "<<m_is_overlayed<<std::endl;
        std::cout<<"SinglePhoton::reconfigure || m_is_textgen: "<<m_is_textgen<<std::endl;

        //Save all spacepoints?
        m_bool_save_sp = pset.get<bool>("SaveSpacepoints",true);

        //Instead of running over ALL evnts in a file, can pass in a file to run over just the run:subrun:event listed
        m_runSelectedEvent    = pset.get<bool>("SelectEvent", false);
        m_selected_event_list = pset.get<std::string>("SelectEventList", "");

        //Studies for photo nuclear EventWeights
        m_runPhotoNuTruth = pset.get<bool>("RunPhotoNu",false); 

        //Ability to save some FULL eventweight components, rather than run later. Useful for systematic studies. Harcoded to two currently (TODO)
        m_runTrueEventweight = pset.get<bool>("RunTrueEventWeight",false); 
        m_true_eventweight_label = pset.get<std::string>("true_eventweight_label","eventweight");
        m_Spline_CV_label = pset.get<std::string>("SplineCVLabel", "eventweight");//4to4aFix");


        //Input ArtRoot data products 
        m_pandoraLabel = pset.get<std::string>("PandoraLabel");
        m_trackLabel = pset.get<std::string>("TrackLabel");
        m_sliceLabel = pset.get<std::string>("SliceLabel","pandora");
        m_showerLabel = pset.get<std::string>("ShowerLabel");
        m_caloLabel = pset.get<std::string>("CaloLabel");
        m_flashLabel = pset.get<std::string>("FlashLabel");
        m_potLabel = pset.get<std::string>("POTLabel");
        m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
        m_badChannelLabel = pset.get<std::string>("BadChannelLabel","badmasks");
        m_showerKalmanLabel = pset.get<std::string>("ShowerTrackFitter","pandoraKalmanShower");
        m_showerKalmanCaloLabel =  pset.get<std::string>("ShowerTrackFitterCalo","pandoraKalmanShowercali");
        m_generatorLabel = pset.get<std::string>("GeneratorLabel","generator");
        m_mcTrackLabel = pset.get<std::string>("MCTrackLabel","mcreco");
        m_mcShowerLabel = pset.get<std::string>("MCShowerLabel","mcreco");
        m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
        m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
        m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");
        m_pidLabel = pset.get<std::string>("ParticleIDLabel","pandoracalipidSCE");
        m_shower3dLabel = pset.get<std::string>("Shower3DLabel","shrreco3d");

        //Flash related variables. 
        m_beamgate_flash_start = pset.get<double>("beamgateStartTime",3.2); //Defaults to MC for now. Probably should change
        m_beamgate_flash_end = pset.get<double>("beamgateEndTime",4.8);

        // m_badChannelProducer = pset.get<std::string>("BadChannelProducer","nfspl1");
        m_badChannelProducer = pset.get<std::string>("BadChannelProducer","simnfspl1");

        //CRT related variables, should run only for RUN3+ enabled
        m_runCRT = pset.get<bool>("runCRT",false);
        m_CRTTzeroLabel = pset.get<std::string>("CRTTzeroLabel","crttzero");
        m_CRTVetoLabel = pset.get<std::string>("CRTVetoLabel","crtveto");
        m_CRTHitProducer = pset.get<std::string>("CRTHitProducer", "crthitcorr");
        m_DTOffset = pset.get<double>("DTOffset" , 68600.); //us, taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter.fcl
        m_Resolution = pset.get<double>("Resolution" ,  1.0); //us, taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter.fcl
        m_DAQHeaderProducer = pset.get<std::string>("DAQHeaderProducer" ,  "daq");

        //Some track calorimetry parameters
        m_track_calo_min_dEdx = pset.get<double>("Min_dEdx",0.005);
        m_track_calo_max_dEdx = pset.get<double>("Max_dEdx", 30);
        m_track_calo_min_dEdx_hits = pset.get<double>("Min_dEdx_hits",5); //might be good?
        m_track_calo_trunc_fraction = pset.get<double>("TruncMeanFraction",20.0);

        //Some shower calorimetry parameters
        m_work_function = pset.get<double>("work_function");
        m_recombination_factor =pset.get<double>("recombination_factor");
        m_gain_mc =pset.get<std::vector<double>>("gain_mc");
        m_gain_data =pset.get<std::vector<double>>("gain_data");
        m_wire_spacing = pset.get<double>("wire_spacing");
        m_width_dqdx_box = pset.get<double>("width_box");
        m_length_dqdx_box = pset.get<double>("length_box");
        m_truthmatching_signaldef = pset.get<std::string>("truthmatching_signaldef");

        //A seperate mode to run over AllPFPs and not just slice particles 
        m_run_all_pfps = pset.get<bool>("runAllPFPs",false);


        //Some paramaters for counting protons & photons
        m_exiting_photon_energy_threshold = pset.get<double>("exiting_photon_energy");
        m_exiting_proton_energy_threshold = pset.get<double>("exiting_proton_energy");
        m_mass_pi0_mev =  139.57;

        //SEAviwer Settings for shower clustering and proton stub finding
        //Have two sets:
        //Base SEAview is for Second Shower Veto
        m_runSEAview = pset.get<bool>("runSEAviewShower", false);
        m_SEAviewHitThreshold = pset.get<double>("SEAviewShowerHitThreshold",25);
        m_SEAviewPlotDistance = pset.get<double>("SEAviewShowerPlotDistance",80);
        m_SEAviewDbscanMinPts = pset.get<double>("SEAviewShowerDBSCANMinPts",8);
        m_SEAviewDbscanEps = pset.get<double>("SEAviewShowerDBSCANEps",4);
        m_SEAviewMaxPtsLinFit = pset.get<double>("SEAviewShowerMaxHitsLinFit",20.0);
        m_SEAviewMakePDF = pset.get<bool>("SEAviewShowerMakePDF",false);
        m_SEAviewNumRecoShower = pset.get<int>("SEAviewShowerNumRecoShower", -1);
        m_SEAviewNumRecoTrack = pset.get<int>("SEAviewShowerNumRecoTrack", -1);

        // Second set is for Proton Stub finding
        m_runSEAviewStub = pset.get<bool>("runSEAviewStub", false);
        m_SEAviewStubHitThreshold = pset.get<double>("SEAviewStubHitThreshold",25);
        m_SEAviewStubPlotDistance = pset.get<double>("SEAviewStubPlotDistance",80);
        m_SEAviewStubDbscanMinPts = pset.get<double>("SEAviewStubDBSCANMinPts",1);
        m_SEAviewStubDbscanEps = pset.get<double>("SEAviewStubDBSCANEps",1);
        m_SEAviewStubMakePDF = pset.get<bool>("SEAviewStubMakePDF",false);
        m_SEAviewStubNumRecoShower = pset.get<int>("SEAviewStubNumRecoShower", -1);
        m_SEAviewStubNumRecoTrack = pset.get<int>("SEAviewStubNumRecoTrack", -1);

        bool_make_sss_plots = true;

        //Misc setup 
        this->setTPCGeom(); 
        rangen = new TRandom3(22);

        //Whats a Delta?
        std::vector<std::string> delta_names = {"Delta++","Delta+","Delta-","Delta0"};
        std::vector<int> delta_pdg_list = {2224,2214,1114,2114};
        for(size_t i=0; i< delta_pdg_list.size(); ++i){
            is_delta_map[delta_pdg_list[i]] = delta_names[i];
            is_delta_map[-delta_pdg_list[i]] ="Anti-"+delta_names[i];
        }


        //Text print event? Depreciated at the moment.
        if (m_print_out_event ){
            out_stream.open("v12_ncdelta_missing_trackshower_events.txt");
            if (!out_stream.is_open()){
                std::cout<<"ERROR output file not open"<<std::endl;
                exit(0);
            }
        }

        // Depreciated SSV run inline
        //std::vector<std::string> inputVars = { "sss_candidate_num_hits", "sss_candidate_num_wires", "sss_candidate_num_ticks", "sss_candidate_PCA", "log10(sss_candidate_impact_parameter)", "log10(sss_candidate_min_dist)", "sss_candidate_impact_parameter/sss_candidate_min_dist", "sss_candidate_energy*0.001", "cos(sss_candidate_angle_to_shower)", "sss_candidate_fit_slope", "sss_candidate_fit_constant", "sss_candidate_plane", "sss_reco_shower_energy*0.001", "2*0.001*0.001*sss_reco_shower_energy*sss_candidate_energy*(1-cos(sss_candidate_angle_to_shower))", "log10(2*0.001*0.001*sss_reco_shower_energy*sss_candidate_energy*(1-cos(sss_candidate_angle_to_shower)))", "sss_candidate_energy*0.001/(sss_reco_shower_energy*0.001)", "sss_candidate_closest_neighbour" };
        //sssVetov1 = new ReadBDT(inputVars);

    }

    //------------------------------------------------------------------------------------------------------------------------------------------



    //--------------------------------------- Primary Filter------------------------------------------------------------------------------------
    // Runs over every artroot event
    bool SinglePhoton::filter(art::Event &evt)
    {

        std::cout<<"---------------------------------------------------------------------------------"<<std::endl;
        std::cout<<"SinglePhoton::analyze()\t||\t On entry: "<<m_number_of_events<<std::endl;

        //Clear all output branches 
        this->ClearVertex();

        //Some event based properties
        m_subrun_counts++;
        m_number_of_events++;
        m_run_number = evt.run();
        m_subrun_number = evt.subRun();
        m_event_number = evt.id().event();


        //if module is run in selected-event mode, and current event is not in the list, skip it
        if(m_runSelectedEvent && !IsEventInList(m_run_number, m_subrun_number, m_event_number)){
            std::cout << "SinglePhoton::analyze()\t||\t event " << m_run_number << "/" << m_subrun_number << "/" << m_event_number << " is not in the list, skip it" << std::endl;
            return true;
        }

        //Timing and TPC info
        auto const ID = *geom->begin<geo::TPCID>();
        m_Cryostat = ID.Cryostat;
        m_TPC = ID.TPC;

        _time2cm = sampling_rate(detClocks) / 1000.0 * theDetector.DriftVelocity( theDetector.Efield(), theDetector.Temperature() );//found in ProtoShowerPandora_tool.cc


        //******************************Setup*****************Setup**************************************/
        //******************************Setup*****************Setup**************************************/
        // OK in this section we will get a bunch of stuff we need later, general associations and maps. These are either from pandora helper stuff or from direct associations. 
        // Make sure under the hood you understand this!
        // ------------------------
        // The basic idea is that here we get every possible vector of data products we might need, along with all maps. e.g 
        // tracks->pfparticles->hits
        // tracks->pfparticles->spacepoints ..etc..
        // And then later we just write pseudo-independant code assuming you have every objecets you want (see analyze_Tracks.h) and assume you have access to everything. 
        //TODO: Think about making these class members, we can access them in the pseudo-indepenant code without passing messy maps.


        // Collect all the hits. We will need these. Lets grab both the handle as well as a vector of art::Ptr as I like both. 
        art::ValidHandle<std::vector<recob::Hit>> const & hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(m_hitfinderLabel); 
        std::vector<art::Ptr<recob::Hit>> hitVector;
        art::fill_ptr_vector(hitVector,hitHandle);

        //Lets do "THE EXACT SAME STUFF" for Optical Flashes
        art::ValidHandle<std::vector<recob::OpFlash>> const & flashHandle  = evt.getValidHandle<std::vector<recob::OpFlash>>(m_flashLabel);
        std::vector<art::Ptr<recob::OpFlash>> flashVector;
        art::fill_ptr_vector(flashVector,flashHandle);

        //tracks
        art::ValidHandle<std::vector<recob::Track>> const & trackHandle  = evt.getValidHandle<std::vector<recob::Track>>(m_trackLabel);
        std::vector<art::Ptr<recob::Track>> trackVector;
        art::fill_ptr_vector(trackVector,trackHandle);

        //BadChannels// Fill later
        art::Handle<std::vector<int> > badChannelHandle;
        std::vector<int> badChannelVector;
        if(evt.getByLabel(m_badChannelProducer, m_badChannelLabel, badChannelHandle)){
            badChannelVector            = *(badChannelHandle);
        }

        //Collect the PFParticles from the event. This is the core!
        art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
        std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
        art::fill_ptr_vector(pfParticleVector,pfParticleHandle);
        //So a cross check
        if (!pfParticleHandle.isValid())
        { mf::LogDebug("SinglePhoton") << "  Failed to find the PFParticles.\n";
            return (m_run_pi0_filter ? false : true) ;
        }

        //get the cluster handle for the dQ/dx calc
        art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pandoraLabel);
        std::vector< art::Ptr<recob::Cluster> > clusterVector;
        art::fill_ptr_vector(clusterVector,clusterHandle);

        // This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;
        // typedef std::map< size_t, art::Ptr<recob::PFParticle>>
        // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
        PFParticleIdMap pfParticleMap;
        this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);

        //Slices
        art::ValidHandle<std::vector<recob::Slice>> const & sliceHandle  = evt.getValidHandle<std::vector<recob::Slice>>(m_pandoraLabel);
        std::vector<art::Ptr<recob::Slice>> sliceVector;
        art::fill_ptr_vector(sliceVector,sliceHandle);

        //And some associations
        art::FindManyP<recob::PFParticle> pfparticles_per_slice(sliceHandle, evt, m_pandoraLabel);
        art::FindManyP<recob::Hit> hits_per_slice(sliceHandle, evt, m_pandoraLabel);

        //Slice to PFParticle
        std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::PFParticle>> > sliceToPFParticlesMap;
        std::map<int, std::vector<art::Ptr<recob::PFParticle>> > sliceIDToPFParticlesMap;
        for(size_t i=0; i< sliceVector.size(); ++i){
            auto slice = sliceVector[i];
            sliceToPFParticlesMap[slice] =pfparticles_per_slice.at(slice.key());
            sliceIDToPFParticlesMap[slice->ID()] = pfparticles_per_slice.at(slice.key());
        }

        //Slice to Hits
        std::map< art::Ptr<recob::Slice>, std::vector<art::Ptr<recob::Hit>> > sliceToHitsMap;
        std::map<int, std::vector<art::Ptr<recob::Hit>> > sliceIDToHitsMap;
        for(size_t i=0; i< sliceVector.size(); ++i){
            auto slice = sliceVector[i];
            sliceToHitsMap[slice] =hits_per_slice.at(slice.key());
            sliceIDToHitsMap[slice->ID()] = hits_per_slice.at(slice.key());
        }

        //And some verticies.        
        art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(m_pandoraLabel);
        std::vector<art::Ptr<recob::Vertex>> vertexVector;
        art::fill_ptr_vector(vertexVector,vertexHandle);
        if(vertexVector.size()>0)  m_number_of_vertices++;

        //PFParticle to Vertices
        art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
        std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            auto pfp = pfParticleVector[i];
            pfParticlesToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
        }

        //------- 3D showers
        art::FindOneP<recob::Shower> showerreco3D_per_pfparticle(pfParticleHandle, evt, m_shower3dLabel);
        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> pfParticlesToShowerReco3DMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            auto pfp = pfParticleVector[i];
            if(!showerreco3D_per_pfparticle.at(pfp.key()).isNull()){
                pfParticlesToShowerReco3DMap[pfp] = showerreco3D_per_pfparticle.at(pfp.key());
            }

        }
        //---------Kalman Track Showers
        art::FindOneP<recob::Track> showerKalman_per_pfparticle(pfParticleHandle, evt, m_showerKalmanLabel);
        std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Track>> pfParticlesToShowerKalmanMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            auto pfp = pfParticleVector[i];
            if(!showerKalman_per_pfparticle.at(pfp.key()).isNull()){ 
                pfParticlesToShowerKalmanMap[pfp] =showerKalman_per_pfparticle.at(pfp.key());
            }
        }

        //----- kalmon Cali
        art::ValidHandle<std::vector<recob::Track>> const & kalmanTrackHandle  = evt.getValidHandle<std::vector<recob::Track>>(m_showerKalmanLabel);
        std::vector<art::Ptr<recob::Track>> kalmanTrackVector;
        art::fill_ptr_vector(kalmanTrackVector,kalmanTrackHandle);

        art::FindManyP<anab::Calorimetry> cali_per_kalmantrack(kalmanTrackHandle, evt, m_showerKalmanCaloLabel);
        std::map<art::Ptr<recob::Track>,std::vector<art::Ptr<anab::Calorimetry>>> kalmanTrackToCaloMap;
        for(size_t i=0; i< kalmanTrackVector.size(); ++i){
            auto trk = kalmanTrackVector[i];
            if(cali_per_kalmantrack.at(trk.key()).size()!=0){
                kalmanTrackToCaloMap[trk] =cali_per_kalmantrack.at(trk.key());
            }
        }

        // Once we have actual verticies, lets concentrate on JUST the neutrino PFParticles for now:
        //--------------------------------
        // Produce two PFParticle vectors containing final-state particles:
        // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
        // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
        std::vector< art::Ptr<recob::PFParticle> > crParticles;
        std::vector< art::Ptr<recob::PFParticle> > nuParticles;
        this->GetFinalStatePFParticleVectors(pfParticleMap, pfParticlesToVerticesMap, crParticles, nuParticles);

        //Fill some more vertex info
        m_vertex_pos_wire_p0 =calcWire(m_vertex_pos_y, m_vertex_pos_z, 0, m_TPC, m_Cryostat, *m_channelMap);
        m_vertex_pos_wire_p1 =calcWire(m_vertex_pos_y, m_vertex_pos_z, 1, m_TPC, m_Cryostat, *m_channelMap);
        m_vertex_pos_wire_p2 =calcWire(m_vertex_pos_y, m_vertex_pos_z, 2, m_TPC, m_Cryostat, *m_channelMap);
        m_vertex_pos_tick = calcTime(m_vertex_pos_x, 2, m_TPC,m_Cryostat, theDetector);
        //std::cout<<calcTime(m_vertex_pos_x, 2, m_TPC,m_Cryostat, *theDetector)<<" "<<calcTime(m_vertex_pos_x, 0, m_TPC,m_Cryostat, *theDetector)<<" "<<calcTime(m_vertex_pos_x, 1, m_TPC,m_Cryostat, *theDetector)<<std::endl;






        //if not running over neutrino slice only, use all pfp's in event
        if (m_run_all_pfps ==true){
            nuParticles = pfParticleVector;
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Spacepoints"<<std::endl;
        //Spacepoint associaitions
        art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
        for(size_t i=0; i< nuParticles.size(); ++i){
            const art::Ptr<recob::PFParticle> pfp = nuParticles[i];
            pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get PandoraMetadata"<<std::endl;
        //add the associaton between PFP and metadata, this is important to look at the slices and scores
        art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt,  m_pandoraLabel);
        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> > pfParticleToMetadataMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            const art::Ptr<recob::PFParticle> pfp = pfParticleVector[i];
            pfParticleToMetadataMap[pfp] =  pfPartToMetadataAssoc.at(pfp.key());
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Clusters"<<std::endl;
        //Get a map between the PFP's and the clusters  they're imporant for the shower dQ/dx
        //Also need a map between clusters and hits
        art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
        art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, m_pandoraLabel);
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
        std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;
        //fill map PFP to Clusters
        for(size_t i=0; i< nuParticles.size(); ++i){
            auto pfp = nuParticles[i];
            pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
        }
        //fill map Cluster to Hits
        for(size_t i=0; i< clusterVector.size(); ++i){
            auto cluster = clusterVector[i];
            clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
        }
        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Build hits to PFP Maps"<<std::endl;


        //OK Here we build two IMPORTANT maps for the analysis, (a) given a PFParticle get a vector of hits..
        //and (b) given a single hit, get the PFParticle it is in (MARK: is it only one? always? RE-MARK: Yes)
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;


        //use pfp->cluster and cluster->hit to build pfp->hit map
        //for each PFP
        for(size_t i=0; i<nuParticles.size(); ++i){
            auto pfp = nuParticles[i];

            //get the associated clusters
            std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp] ;

            //make empty vector to store hits
            std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};


            //for each cluster, get the associated hits
            for (art::Ptr<recob::Cluster> cluster: clusters_vec){
                std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

                //insert hits into vector
                hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
            }

            //fill the map
            pfParticleToHitsMap[pfp] = hits_for_pfp;

        }//for each pfp



        /**************************************************************************
         * For SEAview: grab cosmic-related PFPaticles and recob::Hits
         *
         **************************************************************************/
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > cr_pfParticleToClustersMap;
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > cr_pfParticleToHitsMap;

        //first, collect all daughters of primary cosmic
        int num_primary_cosmic_particle = crParticles.size();
        for(int i =0; i!=num_primary_cosmic_particle; ++i){
            auto& pParticle = crParticles[i];
            for(const size_t daughterId : pParticle->Daughters())
            {
                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";

                crParticles.push_back(pfParticleMap.at(daughterId));
            }
        }

        //second, build PFP to hits map for cosmic-related PFParticles
        for(size_t i=0; i< crParticles.size(); ++i){
            auto pfp = crParticles[i];
            cr_pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
        }
        for(size_t i=0; i< crParticles.size(); ++i){
            auto pfp = crParticles[i];

            // std::cout<<"starting to match to hits for pfp "<<pfp->Self()<<std::endl;
            //get the associated clusters
            std::vector<art::Ptr<recob::Cluster>> clusters_vec  = cr_pfParticleToClustersMap[pfp] ;

            //make empty vector to store hits
            std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};

            // std::cout<<"-- there are "<<clusters_vec.size()<<" associated clusters"<<std::endl;

            //for each cluster, get the associated hits
            for (art::Ptr<recob::Cluster> cluster: clusters_vec){
                std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];

                //   std::cout<<"looking at cluster in pfp "<<pfp->Self()<<" with "<<hits_vec.size() <<" hits"<<std::endl;
                //insert hits into vector
                hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
            }

            //fill the map
            cr_pfParticleToHitsMap[pfp] = hits_for_pfp;
            //std::cout<<"saving a total of "<<hits_for_pfp.size()<<" hits for pfp "<<pfp->Self()<<std::endl;

        }//for each pfp
        std::cout << "Guanqun SEAview test: initial crParticle size: " << num_primary_cosmic_particle << ", current size: " << crParticles.size() << std::endl;
        /********************************************************************
         * End of SEAview: grab cosmic-related PFPaticles and recob::Hits
         *
         **************************************************************************/


        //  Test ground for some slice stuff, dont have it run automatically, ignore for now, mark

        if(false){
            std::cout<<"SliceTest: there are "<<sliceVector.size()<<" slices in this event"<<std::endl;
            for(size_t s =0; s<sliceVector.size(); s++){
                auto slice = sliceVector[s];
                auto pfps = sliceToPFParticlesMap[slice]; 

                std::cout<<"SliceTest: On Slice "<<s<<" it has "<<pfps.size()<<" pfparticles"<<std::endl;
                std::vector<float> nu_scores;
                bool isSelectedSlice = false;
                int primaries = 0;
                int primary_pdg = 0;

                for(auto &pfp: pfps){
                    std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadatas = pfParticleToMetadataMap[pfp];
                    for(auto &meta: metadatas){
                        std::map<std::string, float> propertiesmap  = meta->GetPropertiesMap();
                        //for each of the things in the list
                        if(propertiesmap.count("NuScore")==1){
                            nu_scores.push_back(propertiesmap["NuScore"]);
                        }
                        if(propertiesmap.count("IsNeutrino")==1){
                            isSelectedSlice = true; 
                        }
                    }

                    if (pfp->IsPrimary()) {
                        primaries++;
                        primary_pdg = (pfp->PdgCode());    
                    }
                    /*if (!pfp->IsPrimary()) continue;
                    // Check if this particle is identified as the neutrino
                    const int pdg(pfp->PdgCode());
                    const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
                    if(isNeutrino){
                    isSelectedSlice = true; 
                    }*/
                }

                if(nu_scores.size()>0){
                    double mean  = std::accumulate(nu_scores.begin(), nu_scores.end(), 0.0)/(double)nu_scores.size();
                    if(mean!=nu_scores.front()){
                        std::cout<<"ERROR! Somehow the pfp's in this slice have different nu-scores? IMpossible."<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                    std::cout<<"SliceTest: -- and has a nu_score of "<<nu_scores.front()<<std::endl;
                    std::cout<<"SliceTest: -- with "<<primaries<<" primaries: pdg last: "<<primary_pdg<<std::endl;
                }else{
                    std::cout<<"SliceTest: -- and does not have a nu_score of. "<<std::endl;
                }
                if(isSelectedSlice) std::cout<<"SliceTest: -- -- And is the Selected Neutrino Slice"<<std::endl;

            }
        }// End test of slice metadata





        //******************************  Collecting Info and Analysis  **************************************/
        //******************************  Collecting Info and Analysis  **************************************/

        // These are the vectors to hold the tracks and showers for the final-states of the reconstructed neutrino
        //At this point, nuParticles is a std::vector< art::Ptr<recon::PFParticle>> of the PFParticles that we are interested in.
        //tracks is a vector of recob::Tracks and same for showers.
        //Implicitly, tracks.size() + showers.size() =  nuParticles.size(); At this point I would like two things.
        std::vector< art::Ptr<recob::Track> > tracks;
        std::vector< art::Ptr<recob::Shower> > showers;
        std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap; 
        std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap;

        if(m_is_verbose) std::cout<<"SinglePhoton::analyze() \t||\t Get Tracks and Showers"<<std::endl;

        //Helper function (can be found below) to collect tracks and showers in neutrino slice
        this->CollectTracksAndShowers(nuParticles, pfParticleMap,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);

        bool local_save_sp = true; 
        if(tracks.size()+showers.size()<3 && m_bool_save_sp ){local_save_sp = true;}else{local_save_sp = false;}

        //Track Calorimetry. Bit odd here but bear with me, good to match and fill here
        art::FindManyP<anab::Calorimetry> calo_per_track(trackHandle, evt, m_caloLabel);
        std::map<art::Ptr<recob::Track>, std::vector<art::Ptr<anab::Calorimetry>> > trackToCalorimetryMap;
        //So a cross check
        if (!calo_per_track.isValid())
        {
            mf::LogDebug("SinglePhoton") << "  Failed to get Assns between recob::Track and anab::Calorimetry.\n";
            return (m_run_pi0_filter ? false : true);
        }

        //Loop over all tracks we have to fill calorimetry map
        for(size_t i=0; i< tracks.size(); ++i){
            if(calo_per_track.at(tracks[i].key()).size() ==0){
                std::cerr<<"Track Calorimetry Breaking! the vector of calo_per_track is of length 0 at this track."<<std::endl;
            }
            //size_t calo_size = calo_per_track.at(tracks[i].key()).size();
            //std::cout<<"Track Calo from producer: "<<m_caloLabel<<" has "<<calo_size<<" anab::Calorimetry objects associaed."<<std::endl;
            trackToCalorimetryMap[tracks[i]] = calo_per_track.at(tracks[i].key());
        }

        art::FindOneP<anab::ParticleID> pid_per_track(trackHandle, evt, m_pidLabel);
        std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID> > trackToPIDMap;

        // If we want PID algorithms to run. do so here
        // Build a map to get PID from PFParticles, then call PID collection function
        if(m_use_PID_algorithms){
            for(size_t i=0; i< tracks.size(); ++i){
                trackToPIDMap[tracks[i]] = pid_per_track.at(tracks[i].key());
            }
        }



        //**********************************************************************************************/
        //**********************************************************************************************/
        //---------------------------------- MC TRUTH, MC Only---------------------------
        //**********************************************************************************************/
        //**********************************************************************************************/

        //Get the MCtruth handles and vectors
        std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
        std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;

        //Then build a map from MCparticles to Hits and vice versa
        std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  mcParticleToHitsMap;
        std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  hitToMCParticleMap;

        //Apparrently a MCParticle doesn't know its origin (thanks Andy!)
        //I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa
        //Note which map is which!       //First  is one-to-many.         //Second is one-to-one
        std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;
        std::map<int, art::Ptr<simb::MCParticle> >                                     MCParticleToTrackIdMap;

        std::vector<art::Ptr<sim::MCTrack>> mcTrackVector;
        std::vector<art::Ptr<sim::MCShower>> mcShowerVector;

        std::vector<art::Ptr<simb::MCParticle>> matchedMCParticleVector;
        std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > trackToMCParticleMap;
        std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > showerToMCParticleMap;

        //Given a simb::MCParticle we would like a map to either a sim::MCTrack or sim::MCShower
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCTrack> > MCParticleToMCTrackMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCShower> > MCParticleToMCShowerMap;

        if(m_is_verbose){
            std::cout << "SinglePhoton::analyze()\t||\t Consolidated event summary:" << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t - Number of neutrino final-state PFParticles : " << nuParticles.size() << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t    ... of which are track-like   : " << tracks.size() << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t    ... of which are showers-like : " << showers.size() << "\n";
        }

        //**********************************************************************************************/
        //**********************************************************************************************/

        //and now get the simb::MCparticle to both MCtrack and MCshower maps (just for the MCparticles matched ok).
        if(!m_run_pi0_filter) badChannelMatching<art::Ptr<recob::Track>>(badChannelVector, tracks, trackToNuPFParticleMap, pfParticleToHitsMap,bad_channel_list_fixed_mcc9);
        badChannelMatching<art::Ptr<recob::Track>>(badChannelVector, tracks, trackToNuPFParticleMap, pfParticleToHitsMap,bad_channel_list_fixed_mcc9);


        //*******************************Slices***************************************************************/
        //*******************************Slices***************************************************************/

        //these are all filled in analyze slice
        std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
        std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> primaryPFPSliceIdVec; //stores a pair of only the primary PFP's in the event and the slice ind
        std::map<int, double> sliceIdToNuScoreMap; //map between a slice Id and neutrino score
        std::map<art::Ptr<recob::PFParticle>, bool> PFPToClearCosmicMap; //returns true for clear cosmic, false otherwise
        std::map<art::Ptr<recob::PFParticle>, int> PFPToSliceIdMap; //returns the slice id for all PFP's
        std::map<art::Ptr<recob::PFParticle>,bool> PFPToNuSliceMap;
        std::map<art::Ptr<recob::PFParticle>,double> PFPToTrackScoreMap;
        std::map<int, int> sliceIdToNumPFPsMap;
        std::cout<<"SinglePhoton::analyze::AnalyzeSlice()\t||\t Starting"<<std::endl;

        //Slice helper found in analyze_Slice.h
        this->AnalyzeSlices(pfParticleToMetadataMap, pfParticleMap,  primaryPFPSliceIdVec, sliceIdToNuScoreMap, PFPToClearCosmicMap, PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap);

        if (PFPToSliceIdMap.size() < 1)  std::cout<<"ERROR, not storing PFP's in PFPToSliceIdMap"<<std::endl;
        if(m_is_verbose){

            std::cout<<"SinglePhoton::analyze\t||\tthe number of PPF's with stored clear cosmic info is "<<PFPToClearCosmicMap.size()<<std::endl;
            std::cout<<"SinglePhoton::analyze\t||\tthe number of PFP's stored in the PFPToSliceIdMap is "<<PFPToSliceIdMap.size()<<std::endl;

            for (auto pair:PFPToNuSliceMap){
                auto pfp = pair.first;
                auto is_nuslice = pair.second;
                if (is_nuslice){
                    std::cout<<"pfp in nuslice "<<pfp->Self()<<std::endl;
                }

            }


            for (auto pair:sliceIDToPFParticlesMap){ 
                std::vector<art::Ptr<recob::PFParticle>> pfp_vec = pair.second;
                int slice_id = pair.first;
                //if (slice_vec[0]->Slice() != PFPToSliceIdMap[pfp] )
                for(auto pfp: pfp_vec){
                    if (slice_id != PFPToSliceIdMap[pfp] && PFPToSliceIdMap[pfp]>=0){
                        std::cout<<"sliceIDToPFParticlesMap[slice->ID()] for pfp "<<pfp->Self()<<" is slice "<< slice_id<< "but PFPToSliceIdMap[pfp] = "<<PFPToSliceIdMap[pfp]<<std::endl;
                    }
                }

            }

        }//end verbose



        //******************************* CRT CRT***************************************************************/
        //******************************* CRT CRT***************************************************************/

        std::map<art::Ptr<recob::OpFlash>, std::vector< art::Ptr<crt::CRTHit>>> crtvetoToFlashMap;

        if(m_runCRT){
            art::FindManyP<crt::CRTHit> crtveto_per_flash(flashHandle, evt,   m_CRTVetoLabel);
            for(size_t i=0; i< flashVector.size(); ++i){
                crtvetoToFlashMap[flashVector[i]] = crtveto_per_flash.at(flashVector[i].key());
            }
        }

        art::Handle<std::vector<crt::CRTHit>> crthit_h; //only filled when there are hits, otherwise empty
        art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
        double evt_timeGPS_nsec = -999 ;
        if(m_runCRT){
            evt.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);
            evt.getByLabel(m_CRTHitProducer, crthit_h);
            raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
            art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
            evt_timeGPS_nsec = evtTimeGPS.timeLow(); 
            std::cout<<"SinglePhoton::analyze \t||\t Got CRT hits"<<std::endl;
        }

        //Analize the CRT flashes. Found in analyze_OpFlashes.h
        this->AnalyzeFlashes(flashVector, crthit_h, evt_timeGPS_nsec, crtvetoToFlashMap);



        //******************************* Common Optical Filter **************************************************************/
        //******************************* Common Optical Filter **************************************************************/
        //Raw Optical fltr
        art::Handle<uboone::UbooneOpticalFilter> uBooNE_common_optFltr;
      if(evt.getByLabel("opfiltercommon", uBooNE_common_optFltr)){
            m_flash_optfltr_pe_beam     = uBooNE_common_optFltr->PE_Beam();
            m_flash_optfltr_pe_beam_tot = uBooNE_common_optFltr->PE_Beam_Total();
            m_flash_optfltr_pe_veto     = uBooNE_common_optFltr->PE_Veto();
            m_flash_optfltr_pe_veto_tot = uBooNE_common_optFltr->PE_Veto_Total();
        }else{
            m_flash_optfltr_pe_beam     = -999;
            m_flash_optfltr_pe_beam_tot = -999;
            m_flash_optfltr_pe_veto     = -999;
            m_flash_optfltr_pe_veto_tot = -999;
            std::cout<<"No opfiltercommon product:"<<std::endl;
        }



        //*******************************  Tracks **************************************************************/
        //*******************************  Tracks **************************************************************/
        std::cout<<"SinglePhoton::analyze \t||\t Start on Track Analysis "<<std::endl;

        //found in analyze_Tracks.h
        this->AnalyzeTracks(tracks, trackToNuPFParticleMap, pfParticleToHitsMap,  pfParticleToSpacePointsMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap,  PFPToTrackScoreMap, PFPToNuSliceMap,pfParticleMap);
        this->AnalyzeTrackCalo(tracks,   trackToCalorimetryMap);

        for(size_t i_trk = 0; i_trk<tracks.size();i_trk++){
            const art::Ptr<recob::Track> t = tracks[i_trk];
            const art::Ptr<recob::PFParticle> pfp = trackToNuPFParticleMap[t];
            const std::vector< art::Ptr<recob::SpacePoint> > trk_spacepoints = pfParticleToSpacePointsMap[pfp];

            std::vector<double> tmp_sp_x;
            std::vector<double> tmp_sp_y;
            std::vector<double> tmp_sp_z;

            if(local_save_sp){
            for(auto &sp: trk_spacepoints){
                    tmp_sp_x.push_back(sp->XYZ()[0]);
                    tmp_sp_y.push_back(sp->XYZ()[1]);
                    tmp_sp_z.push_back(sp->XYZ()[2]);
                }

                    m_reco_track_spacepoint_x.push_back(tmp_sp_x);
                    m_reco_track_spacepoint_y.push_back(tmp_sp_y);
                    m_reco_track_spacepoint_z.push_back(tmp_sp_z);

            }
        }


        //Run over PID?
        if(m_use_PID_algorithms)  this->CollectPID(tracks, trackToPIDMap);


        //*******************************  Showers **************************************************************/
        //*******************************  Showers **************************************************************/
        std::cout<<"SinglePhoton::analyze \t||\t Start on Shower Analysis "<<std::endl;

        //found in analyze_Showers.h
        this->AnalyzeShowers(showers,showerToNuPFParticleMap, pfParticleToHitsMap, pfParticleToClustersMap, clusterToHitsMap,sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap, PFPToNuSliceMap, PFPToTrackScoreMap,pfParticleMap,pfParticlesToShowerReco3DMap); 
        this->AnalyzeKalmanShowers(showers,showerToNuPFParticleMap,pfParticlesToShowerKalmanMap, kalmanTrackToCaloMap, pfParticleToHitsMap);


        //Some misc things thrown in here rather than in a proper helper function. TODO. fix
        //Calc a fake shower "end" distance. How to define an end distance? good question
        for(size_t i_shr = 0; i_shr<showers.size();i_shr++){
            const art::Ptr<recob::Shower> s = showers[i_shr];
            const art::Ptr<recob::PFParticle> pfp = showerToNuPFParticleMap[s];
            const std::vector< art::Ptr<recob::SpacePoint> > shr_spacepoints = pfParticleToSpacePointsMap[pfp];

            m_reco_shower_end_dist_to_active_TPC[i_shr] = 99999;
            m_reco_shower_end_dist_to_SCB[i_shr] = 99999;

            std::vector<double> tmp_sp_x;
            std::vector<double> tmp_sp_y;
            std::vector<double> tmp_sp_z;

            for(auto &sp: shr_spacepoints){
                std::vector<double> tmp_spt = {sp->XYZ()[0],sp->XYZ()[1] , sp->XYZ()[2]};
                m_reco_shower_end_dist_to_active_TPC[i_shr] = std::min(m_reco_shower_end_dist_to_active_TPC[i_shr], distToTPCActive(tmp_spt));
                double tmo;
                this->distToSCB(tmo,tmp_spt);
                m_reco_shower_end_dist_to_SCB[i_shr] = std::min(m_reco_shower_end_dist_to_SCB[i_shr],tmo);

                //This section runs for only 1 shower events for purpose of testing delta specifics 

                if(local_save_sp){
                    tmp_sp_x.push_back(sp->XYZ()[0]);
                    tmp_sp_y.push_back(sp->XYZ()[1]);
                    tmp_sp_z.push_back(sp->XYZ()[2]);
                }

            }
                    m_reco_shower_spacepoint_x.push_back(tmp_sp_x);
                    m_reco_shower_spacepoint_y.push_back(tmp_sp_y);
                    m_reco_shower_spacepoint_z.push_back(tmp_sp_z);

        }


        //*******************************   MCTruth  **************************************************************/
        //*******************************   MCTruth  **************************************************************/

        //Grab the backtracker info for MCTruth Matching
        art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);

        // MCTruth, MCParticle, MCNeutrino information all comes directly from GENIE.
        // MCShower and MCTrack come from energy depositions in GEANT4

        //Only run if its not data :)
        if(!m_is_data){

            std::vector<art::Ptr<simb::GTruth>> gTruthVector;
            if(!m_is_textgen){

                // if Text is in the generator label, skip it. TODO this is a bit simple but works, maybe add a boolean
                art::ValidHandle<std::vector<simb::GTruth>> const & gTruthHandle= evt.getValidHandle<std::vector<simb::GTruth>>(m_generatorLabel);
                art::fill_ptr_vector(gTruthVector,gTruthHandle);
                if(m_is_verbose){
                    for(size_t p=0; p< gTruthVector.size();p++) std::cout<<gTruthVector[p]<<" "<<*gTruthVector[p]<<std::endl;
                }
            }

            //get MCTruth (GENIE)
            art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle= evt.getValidHandle<std::vector<simb::MCTruth>>(m_generatorLabel);
            art::fill_ptr_vector(mcTruthVector,mcTruthHandle);

            //get MCPartilces (GEANT4)
            art::ValidHandle<std::vector<simb::MCParticle>> const & mcParticleHandle= evt.getValidHandle<std::vector<simb::MCParticle>>(m_geantModuleLabel);
            art::fill_ptr_vector(mcParticleVector,mcParticleHandle);

            //Found inanalyze_Geant4.h
            //Currently just saves first 2 particles. TODO have a input list of things to save. Dont want to save everything!!
            this->AnalyzeGeant4(mcParticleVector);


            //Get the MCParticles (move to do this ourselves later)
            this->CollectMCParticles(evt, m_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap, MCParticleToTrackIdMap);


            //mcc9 march miniretreat fix
            std::vector<art::Ptr<simb::MCParticle>> particle_vec; //vector of all MCParticles associated with a given hit in the reco PFP
            std::vector<anab::BackTrackerHitMatchingData const *> match_vec; //vector of some backtracker thing
            m_test_matched_hits = 0;

            for(size_t j=0; j<hitVector.size();j++){
                const art::Ptr<recob::Hit> hit = hitVector[j];

                particle_vec.clear(); match_vec.clear(); //only store per hit
                mcparticles_per_hit.get(hit.key(), particle_vec, match_vec);
                if(particle_vec.size() > 0){
                    m_test_matched_hits++;
                }
            }


            //Important map, given a MCparticle, whats the "hits" associated
            this->BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector,  mcParticleToHitsMap, hitToMCParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters,  MCParticleToTrackIdMap);


            //The recoMC was originally templated for any track shower, but sufficient differences in showers emerged to have stand alone sadly
            std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::track"<<std::endl;
            std::vector<double> trk_overlay_vec = recoMCmatching<art::Ptr<recob::Track>>( tracks, trackToMCParticleMap, trackToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector);
            std::cout<<"SinglePhoton\t||\t Starting backtracker on recob::shower"<<std::endl;
            this->showerRecoMCmatching(showers, showerToMCParticleMap, showerToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, matchedMCParticleVector, pfParticleMap,  MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap, PFPToNuSliceMap);

            //photoNuclearTesting(matchedMCParticleVector);

            std::cout<<"Starting outside RecoMCTracks "<<std::endl;
            this->RecoMCTracks(tracks, trackToNuPFParticleMap, trackToMCParticleMap, MCParticleToMCTruthMap,mcParticleVector, MCParticleToTrackIdMap, sliceIdToNuScoreMap, PFPToClearCosmicMap,  PFPToSliceIdMap,trk_overlay_vec);






            std::cout<<"Starting outside AnalyzeMCTruths "<<std::endl;
            this->AnalyzeMCTruths(mcTruthVector, mcParticleVector);


            std::cout<<"Starting AnalyzeEventWeight"<<std::endl;
            if(!m_is_textgen){// if Text is in the generator label, skip it. wont contain any

                //This is important for rebuilding MCFlux and GTruth to do eventweight 
                this->AnalyzeEventWeight(evt);


                //added since last time?
                std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind

                /*   std::vector<std::pair<art::Ptr<recob::PFParticle>,int>> allPFPSliceIdVec; //stores a pair of all PFP's in the event and the slice ind
                     std::map<int, std::vector<art::Ptr<recob::PFParticle>>> sliceIdToPFPMap; //this is an alternative, stores all the PFP's but organized by slice ID
                     std::cout<<"SinglePhoton::AnalyzeSlice()\t||\t Starting"<<std::endl;
                     this->AnalyzeSlices( pfParticleToMetadataMap, pfParticleMap, allPFPSliceIdVec, sliceIdToPFPMap);
                     std::cout<<"There are "<< allPFPSliceIdVec.size()<<" pfp-slice id matches stored in the vector"<<std::endl;
                     if (showers.size()>0){
                     std::cout<<"the shower at 0 is in slice "<<this->GetShowerSlice(showers[0], showerToNuPFParticleMap, allPFPSliceIdVec)<<std::endl;
                     }
                     */

                //this one was for testing, leaving out for now
                // this->FindSignalSlice( m_truthmatching_signaldef, MCParticleToTrackIdMap, showerToNuPFParticleMap , allPFPSliceIdVec, showerToMCParticleMap, trackToNuPFParticleMap, trackToMCParticleMap);

                //if(!m_run_pi0_filter){
                //    this->SecondShowerSearch(tracks,  trackToNuPFParticleMap, showers, showerToNuPFParticleMap, pfParticleToHitsMap, PFPToSliceIdMap, sliceIDToHitsMap,mcparticles_per_hit, matchedMCParticleVector, pfParticleMap,  MCParticleToTrackIdMap);
                //}

                std::cout<<"filling info in ncdelta slice tree"<<std::endl;
                this->AnalyzeRecoMCSlices( m_truthmatching_signaldef, MCParticleToTrackIdMap, showerToNuPFParticleMap , allPFPSliceIdVec, showerToMCParticleMap, trackToNuPFParticleMap, trackToMCParticleMap,  PFPToSliceIdMap);

                if (m_print_out_event){
                    if (m_matched_signal_shower_num != 1 || m_matched_signal_track_num != 1){
                        out_stream <<"run subrunevent "<<m_run_number<<" "<<m_subrun_number<<" "<<m_event_number<<"\n";
                    }
                }

                //This section gets some important splines and weights 
                std::cout<<"Going to grab eventweightSplines for CCQE genie fix, won't be necessary long term"<<std::endl;
                //Mark: It was long term....
                art::Handle<std::vector<evwgh::MCEventWeight>>  ev_evw ;
                if(    evt.getByLabel(m_Spline_CV_label,ev_evw)){

                    std::map<std::string, std::vector<double>> const & weight_map = ev_evw->front().fWeight;
                    if(ev_evw->size() > 1) std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"<< "WARNING: eventweight slice genie fix has more than one entry\n";

                    for (auto const& x : weight_map){
                        std::cout << x.first  // string (key)
                            << ':' 
                            << x.second.size() << std::endl ;
                        if(x.second.size()==1 && x.first == "splines_general_Spline"){
                            m_genie_spline_weight = x.second.front();
                            std::cout<<"Its a spline fix, value: "<<m_genie_spline_weight<<std::endl;
                        }
                        if(x.second.size()==1 && ( x.first == "TunedCentralValue_Genie" || x.first == "TunedCentralValue_UBGenie")){
                            m_genie_CV_tune_weight = x.second.front();
                            std::cout<<"Its a CV fix, value:  "<<m_genie_CV_tune_weight<<std::endl;
                        }
                    }

                }else{
                    std::cout<<"WARNING!!  No data product called eventweightSplines. Setting them to 1. Make sure this is what you want ok."<<std::endl;
                    m_genie_spline_weight =1.0;
                    m_genie_CV_tune_weight =1.0;
                }



                //*******************************   PhotoNuclear Absorption  **************************************************************/
                //*******************************   PhotoNuclear Absorption  **************************************************************/
                m_photonu_weight_low = -999;
                m_photonu_weight_high = -999;
                if(m_runPhotoNuTruth){
                    art::Handle<std::vector<evwgh::MCEventWeight>>  ev_evw_ph ;
                    if(    evt.getByLabel("eventweight",ev_evw_ph)){
                        std::map<std::string, std::vector<double>> const & weight_map = ev_evw_ph->front().fWeight;
                        for (auto const& x : weight_map){
                            std::cout << x.first  // string (key)
                                << ':' 
                                << x.second.size() << std::endl ;
                            if(x.first == "photonuclear_photon_PhotoNuclear"){
                                auto vec  = x.second;
                                double ph_low = vec[1];
                                double ph_high = vec[0];
                                std::cout<<"PhotoNuBit: "<<ph_low<<" "<<ph_high<<std::endl;
                                m_photonu_weight_low = ph_low;
                                m_photonu_weight_high = ph_high;
                            }
                        }

                    }
                }


                //*******************************   True EventWeight  **************************************************************/
                //*******************************   True EventWeight  **************************************************************/

                if(m_runTrueEventweight){

                    art::ValidHandle<std::vector<evwgh::MCEventWeight>> const & ev_evw_true =  evt.getValidHandle<std::vector<evwgh::MCEventWeight>>(m_true_eventweight_label);
                    std::map<std::string, std::vector<double>> const & weight_map = ev_evw_true->front().fWeight;
                    if(ev_evw_true->size() > 1) {
                        std::cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\n"
                            << "WARNING: eventweight has more than one entry\n";
                    }
                    fmcweight=weight_map;
                }

            }//end NOT textgen 

            //if it IS textgen, do we have additonal info saved?
            if(m_is_textgen){

                art::Handle<std::vector<double>> text_info_handle;
                if(evt.getByLabel(m_generatorLabel, text_info_handle)){
                std::cout<<"Textgen Additional Info Found: Size: "<<text_info_handle->size()<<"\n";
                for(size_t k=0; k<text_info_handle->size(); k++){std::cout<<text_info_handle->at(k)<<" ";
                    m_textgen_info.push_back(text_info_handle->at(k));
                }
                std::cout<<"\n";
                }

            }



            std::cout<<"SinglePhoton::analyze\t||\t finnished loop for this event"<<std::endl;
        }//end data loop



        //*******************************   Isolation (SSV precursor)  **************************************************************/
        //*******************************   Isolation (SSV precursor)  **************************************************************/
        if(!m_run_all_pfps && ! m_run_pi0_filter) this-> IsolationStudy(tracks,  trackToNuPFParticleMap, showers, showerToNuPFParticleMap, pfParticleToHitsMap, PFPToSliceIdMap, sliceIDToHitsMap);




        // ################################################### SEAview SEAview #########################################################
        // ################################################### SEAview SEAview #########################################################


        // ################################################### Proton Stub ###########################################
        // ------------- stub clustering ---------------------------
        std::cout << "----------------- Stub clustering --------------------------- " << std::endl;
        std::cout << "SEAview Stub formation: " << (m_runSEAviewStub ? "true" : "false" ) << " nshower requirement: " << m_SEAviewStubNumRecoShower << ", actual num shower: " << showers.size() << " | ntrack requirement: " << m_SEAviewStubNumRecoTrack << ", actual num track: " << tracks.size() << std::endl;

        if(!m_run_pi0_filter && m_runSEAviewStub && (m_SEAviewStubNumRecoShower == -1 || (int)showers.size()== m_SEAviewStubNumRecoShower) && (m_SEAviewStubNumRecoTrack == -1 || (int)tracks.size() == m_SEAviewStubNumRecoTrack)){    

            // grab all hits in the slice of the reco shower
            art::Ptr<recob::Shower> p_shr = showers.front();
            art::Ptr<recob::PFParticle> p_pfp = showerToNuPFParticleMap[p_shr];
            std::vector<art::Ptr<recob::Hit>> p_hits = pfParticleToHitsMap[p_pfp];

            int p_sliceid = PFPToSliceIdMap[p_pfp];
            auto p_slice_hits =    sliceIDToHitsMap[p_sliceid];

            std::string uniq_tag = "HitThres_"+ std::to_string(static_cast<int>(m_SEAviewStubHitThreshold)) + "_" + std::to_string(m_run_number)+"_"+std::to_string(m_subrun_number)+"_"+std::to_string(m_event_number);

            //Setup seaviewr object
            seaview::SEAviewer sevd("Stub_"+uniq_tag, geom, m_channelMap, theDetector );
            //Pass in any bad channels you like
            sevd.setBadChannelList(bad_channel_list_fixed_mcc9);
            //Give it a vertex to center around
            sevd.loadVertex(m_vertex_pos_x,m_vertex_pos_y, m_vertex_pos_z);

            //Add the hits from just this slice, as well as hits within 150cm of the vertex 
            sevd.addHitsToConsider(hitVector);   // std::vector<art::Ptr<recob::Hit>> 
            sevd.filterConsideredHits(150); //remve hits that're not within 150cm of the vertex on 2D view
            sevd.addHitsToConsider(p_slice_hits);

            sevd.addAllHits(hitVector); // std::vector<art::Ptr<recob::Hit>> 
            sevd.setHitThreshold(m_SEAviewStubHitThreshold); 


            //Add all the "nice " PFParticle Hits, as well as what to label
            //sevd.addPFParticleHits(p_hits, "Shower");  //std::vector<art::Ptr<recob::Hit>> and std::string
            sevd.addPFParticleHits(p_hits, "Shower", m_reco_shower_energy_max[0], m_reco_shower_conversion_distance[0], m_reco_shower_impact_parameter[0]);  //std::vector<art::Ptr<recob::Hit>> and std::string

            //and add the SingleShower we like
            sevd.addShower(p_shr); // art::Ptr<recob::Shower>

            //Add all track PFP
            int i_trk = 0;
            for(auto &trk: tracks){
                art::Ptr<recob::PFParticle> p_pfp_trk = trackToNuPFParticleMap[trk];
                std::vector<art::Ptr<recob::Hit>> p_hits_trk = pfParticleToHitsMap[p_pfp_trk];
                //sevd.addPFParticleHits(p_hits_trk,"track");
                sevd.addPFParticleHits(p_hits_trk,"track", m_reco_track_length[i_trk], m_reco_track_spacepoint_principal0[i_trk]);
                sevd.addTrack(trk);
                ++i_trk;
            }

            //Add all cosmic-relatd PFP
            for(auto &cr: crParticles){
                std::vector<art::Ptr<recob::Hit>> p_hits_cr = cr_pfParticleToHitsMap[cr];
                sevd.addPFParticleHits(p_hits_cr,"cosmic");
            }


            //We then calculate Unassociated hits, i.e the hits not associated to the "Shower" or tracksyou passed in. 
            auto vnh= sevd.calcUnassociatedHits();
            m_trackstub_num_unassociated_hits =vnh[1]+vnh[2];
            m_trackstub_unassociated_hits_below_threshold = vnh[2];
            m_trackstub_associated_hits = vnh[0]-vnh[1]-vnh[2];

            //Recluster, group unassociated hits into different clusters
            sevd.runseaDBSCAN(m_SEAviewStubDbscanMinPts, m_SEAviewStubDbscanEps);

            //And some plotting
            // If we want to plot pdfs again later, then we can't plot here
            //if(m_SEAviewStubMakePDF) sevd.Print(m_SEAviewStubPlotDistance);

            //Analyze formed clusters and save info
            std::vector<seaview::cluster> vec_SEAclusters ;
            sevd.analyzeTrackLikeClusters(m_SEAviewStubDbscanEps, showerToNuPFParticleMap, pfParticleToHitsMap, vec_SEAclusters);


            //And save to file.
            std::cout<<"After SEAview we have "<<vec_SEAclusters.size()<<" Stub clusters to chat about"<<std::endl;

            m_trackstub_num_candidates = 0;
            for(size_t c=0; c< vec_SEAclusters.size(); c++){
                auto& clu = vec_SEAclusters.at(c); //type: seaview::cluster
                int pl = clu.getPlane();
                auto hitz = clu.getHits();
                double Ep = this->CalcEShowerPlane(hitz,pl); 
                int remerge = clu.getShowerRemerge();
                seaview::cluster_score * ssscorz = clu.getScore();

                std::cout<<c<<" "<<pl<<" "<<Ep<<" "<<clu.getMinHitImpactParam()<<" "<<clu.getMinHitConvDist()<<" "<<clu.getMinHitIOC()<<" "<<clu.getMeanADC()<<" "<<clu.getADCrms()<<" "<< clu.getLinearChi() << " " << remerge<<std::endl;

                //if the cluster is too close to the recob::shower, then do not include it
                if(remerge>=0 && remerge< (int)m_reco_shower_reclustered_energy_plane2.size()){
                    //decide not to add energy of the cluster to reco shower if it's matched
                    //
                    //if(pl==0)m_reco_shower_reclustered_energy_plane0[remerge]+=Ep;
                    //if(pl==1)m_reco_shower_reclustered_energy_plane1[remerge]+=Ep;
                    //if(pl==2)m_reco_shower_reclustered_energy_plane2[remerge]+=Ep;

                    continue;// Dont include this as a viable cluster!
                }

                ++m_trackstub_num_candidates;
                //determine if this cluster is in neutrino slice
                m_trackstub_candidate_in_nu_slice.push_back(clu.InNuSlice(sliceIDToHitsMap, p_sliceid));

                //Fill All the bits
                m_trackstub_candidate_num_hits.push_back((int)hitz.size());
                m_trackstub_candidate_num_wires.push_back((int)ssscorz->n_wires);
                m_trackstub_candidate_num_ticks.push_back((int)ssscorz->n_ticks);
                m_trackstub_candidate_plane.push_back(pl);
                m_trackstub_candidate_PCA.push_back(ssscorz->pca_0);
                m_trackstub_candidate_mean_tick.push_back(ssscorz->mean_tick);
                m_trackstub_candidate_max_tick.push_back(ssscorz->max_tick);
                m_trackstub_candidate_min_tick.push_back(ssscorz->min_tick);
                m_trackstub_candidate_min_wire.push_back(ssscorz->min_wire);
                m_trackstub_candidate_max_wire.push_back(ssscorz->max_wire);
                m_trackstub_candidate_mean_wire.push_back(ssscorz->mean_wire);
                m_trackstub_candidate_min_dist.push_back(ssscorz->min_dist);
                m_trackstub_candidate_min_impact_parameter_to_shower.push_back(clu.getMinHitImpactParam());
                m_trackstub_candidate_min_conversion_dist_to_shower_start.push_back(clu.getMinHitConvDist());
                m_trackstub_candidate_min_ioc_to_shower_start.push_back(clu.getMinHitIOC());
                m_trackstub_candidate_ioc_based_length.push_back(clu.getIOCbasedLength());
                m_trackstub_candidate_wire_tick_based_length.push_back(clu.getWireTickBasedLength());
                m_trackstub_candidate_mean_ADC_first_half.push_back(clu.getMeanADCFirstHalf());
                m_trackstub_candidate_mean_ADC_second_half.push_back(clu.getMeanADCSecondHalf());
                m_trackstub_candidate_mean_ADC_first_to_second_ratio.push_back(clu.getMeanADCRatio());
                m_trackstub_candidate_track_angle_wrt_shower_direction.push_back(clu.getTrackAngleToShowerDirection());
                m_trackstub_candidate_linear_fit_chi2.push_back(clu.getLinearChi());
                m_trackstub_candidate_mean_ADC.push_back(clu.getMeanADC());
                m_trackstub_candidate_ADC_RMS.push_back(clu.getADCrms());
                m_trackstub_candidate_energy.push_back(Ep);
                m_trackstub_candidate_remerge.push_back(remerge);


                //MCTruth matching for pi0's
                if(m_is_data){
                    m_trackstub_candidate_matched.push_back(-1);
                    m_trackstub_candidate_pdg.push_back(-1);
                    m_trackstub_candidate_parent_pdg.push_back(-1);
                    m_trackstub_candidate_trackid.push_back(-1);
                    m_trackstub_candidate_true_energy.push_back(-1);
                    m_trackstub_candidate_overlay_fraction.push_back(-1);
                    m_trackstub_candidate_matched_energy_fraction_best_plane.push_back(-1);
                }else{

                    auto ssmatched = this->SecondShowerMatching(hitz, mcparticles_per_hit, mcParticleVector, pfParticleMap,  MCParticleToTrackIdMap);
                    m_trackstub_candidate_matched.push_back(ssmatched[0]);
                    m_trackstub_candidate_pdg.push_back(ssmatched[1]);
                    m_trackstub_candidate_parent_pdg.push_back(ssmatched[2]);
                    m_trackstub_candidate_trackid.push_back(ssmatched[3]);
                    m_trackstub_candidate_true_energy.push_back(ssmatched[4]);
                    m_trackstub_candidate_overlay_fraction.push_back(ssmatched[5]);
                    m_trackstub_candidate_matched_energy_fraction_best_plane.push_back(ssmatched[6]);

                    //Guanqun: print out (best-matched) truth information of the cluster
                    
                    if(m_is_verbose){

                    std::cout << "Cluster: " << m_trackstub_num_candidates-1  << " plane: " << m_trackstub_candidate_plane.back() << ", energy: " << m_trackstub_candidate_energy.back() << ", min IOC of hit(wrt shower): " << m_trackstub_candidate_min_ioc_to_shower_start.back() << "\n";
                    std::cout << "Cluster is matched: " << m_trackstub_candidate_matched.back() << ", matched PDG: " << m_trackstub_candidate_pdg.back() << " track ID: " << m_trackstub_candidate_trackid.back() << " overlay fraction: " << m_trackstub_candidate_overlay_fraction.back() << std::endl; 
                    std::cout << "===============================================================" << std::endl;
                    }
                }

                sevd.SetClusterLegend(c, m_trackstub_candidate_energy.back(),  m_trackstub_candidate_matched.back(), m_trackstub_candidate_pdg.back() , m_trackstub_candidate_overlay_fraction.back() );


            } //end of cluster loop

            // Plot the event
            if(m_SEAviewStubMakePDF){
                sevd.Print(m_SEAviewStubPlotDistance);
            }

            //group clusters HERE
            std::pair<int, std::pair<std::vector<std::vector<double>>, std::vector<double>> > group_result = GroupClusterCandidate(m_trackstub_num_candidates,  m_trackstub_candidate_plane, m_trackstub_candidate_max_tick, m_trackstub_candidate_min_tick);
            m_trackstub_num_candidate_groups = group_result.first;
            m_grouped_trackstub_candidate_indices = group_result.second.first;
            m_trackstub_candidate_group_timeoverlap_fraction = group_result.second.second;
        }


        // --------------- shower clustering --------------------------
        std::cout << "------------- Shower clustering --------------------" << std::endl;
        std::cout << "SEAview Shower cluster formation: " << (m_runSEAview ? "true" : "false" ) << " nshower requirement: " << m_SEAviewNumRecoShower << ", actual num shower: " << showers.size() << " | ntrack requirement: " << m_SEAviewNumRecoTrack << ", actual num track: " << tracks.size() << std::endl;

        if(!m_run_pi0_filter &&  m_runSEAview && (m_SEAviewNumRecoShower == -1 || (int)showers.size()== m_SEAviewNumRecoShower) && (m_SEAviewNumRecoTrack == -1 || (int)tracks.size() == m_SEAviewNumRecoTrack) ){    

            art::Ptr<recob::Shower> p_shr = showers.front();
            art::Ptr<recob::PFParticle> p_pfp = showerToNuPFParticleMap[p_shr];
            std::vector<art::Ptr<recob::Hit>> p_hits = pfParticleToHitsMap[p_pfp];


            int p_sliceid = PFPToSliceIdMap[p_pfp];
            auto p_slice_hits =    sliceIDToHitsMap[p_sliceid];

            std::string uniq_tag = "HitThres_"+ std::to_string(static_cast<int>(m_SEAviewHitThreshold)) + "_" + std::to_string(m_run_number)+"_"+std::to_string(m_subrun_number)+"_"+std::to_string(m_event_number);

            //Setup seaviewr object
            seaview::SEAviewer sevd("Shower_"+uniq_tag, geom, m_channelMap, theDetector );
            //Pass in any bad channels you like
            sevd.setBadChannelList(bad_channel_list_fixed_mcc9);
            //Give it a vertex to center around
            sevd.loadVertex(m_vertex_pos_x,m_vertex_pos_y, m_vertex_pos_z);

            //Add hits to consider for clustering 
            //sevd.addHitsToConsider(hitVector);// DONT do this yet, need something smarter for SSV
            //sevd.filterConsideredHits(150); 
            sevd.addHitsToConsider(p_slice_hits);

            //Add all hits in the events
            sevd.addAllHits(hitVector); // std::vector<art::Ptr<recob::Hit>> 
            sevd.setHitThreshold(m_SEAviewHitThreshold); 

            //Add all the "nice " PFParticle Hits, as well as what to label
            //sevd.addPFParticleHits(p_hits, "Shower");  //std::vector<art::Ptr<recob::Hit>> and std::string
            sevd.addPFParticleHits(p_hits, "Shower", m_reco_shower_energy_max[0], m_reco_shower_conversion_distance[0], m_reco_shower_impact_parameter[0]);  //std::vector<art::Ptr<recob::Hit>> and std::string

            //and add the SingleShower we like
            sevd.addShower(p_shr); // art::Ptr<recob::Shower>

            //Add all track PFP
            int i_trk = 0;
            for(auto &trk: tracks){
                art::Ptr<recob::PFParticle> p_pfp_trk = trackToNuPFParticleMap[trk];
                std::vector<art::Ptr<recob::Hit>> p_hits_trk = pfParticleToHitsMap[p_pfp_trk];
                //sevd.addPFParticleHits(p_hits_trk,"track");
                sevd.addPFParticleHits(p_hits_trk,"track", m_reco_track_length[i_trk], m_reco_track_spacepoint_principal0[i_trk]);
                sevd.addTrack(trk);
                ++i_trk;
            }

            //Add all cosmic-relatd PFP // DONT do this yet, see line 1206
            /*for(auto &cr: crParticles){
                std::vector<art::Ptr<recob::Hit>> p_hits_cr = cr_pfParticleToHitsMap[cr];
                sevd.addPFParticleHits(p_hits_cr,"cosmic");
            }
            */

            //We then calculate Unassociated hits, i.e the hits not associated to the "Shower" or tracksyou passed in. 
            auto vnh= sevd.calcUnassociatedHits();
            m_sss_num_unassociated_hits =vnh[1]+vnh[2];
            m_sss_num_unassociated_hits_below_threshold = vnh[2];
            m_sss_num_associated_hits = vnh[0]-vnh[1]-vnh[2];

            //Recluster, group unassociated hits into different clusters
            sevd.runseaDBSCAN(m_SEAviewDbscanMinPts, m_SEAviewDbscanEps);


            //This is the place I will put the new Second Shower Search
            std::vector<seaview::cluster> vec_SEAclusters ;
            sevd.analyzeShowerLikeClusters(m_SEAviewDbscanEps, showerToNuPFParticleMap, pfParticleToHitsMap, vec_SEAclusters);

            //And save to file.
            std::cout<<"After SEAview we have "<<vec_SEAclusters.size()<<" Shower clusters to chat about"<<std::endl;

            m_sss_num_candidates = 0;
            for(size_t c=0; c< vec_SEAclusters.size(); c++){
                auto clu = vec_SEAclusters.at(c); //type: seaview::cluster
                int pl = clu.getPlane();
                auto hitz = clu.getHits();
                double Ep = this->CalcEShowerPlane(hitz,pl); 
                int remerge = clu.getShowerRemerge();
                seaview::cluster_score * ssscorz = clu.getScore();

                if(m_is_verbose) std::cout<<c<<" "<<pl<<" "<<Ep<<" "<<clu.getImpactParam()<<" "<<clu.getFitSlope()<<" "<<clu.getFitCons()<<" "<<clu.getMeanADC() << " " << clu.getADCrms() << " "<<clu.getAngleWRTShower()<<" "<<remerge<<std::endl;

                //if the cluster is too close to the recob::shower, then do not include it
                if(remerge>=0 && remerge< (int)m_reco_shower_reclustered_energy_plane2.size()){
                    if(pl==0)m_reco_shower_reclustered_energy_plane0[remerge]+=Ep;
                    if(pl==1)m_reco_shower_reclustered_energy_plane1[remerge]+=Ep;
                    if(pl==2)m_reco_shower_reclustered_energy_plane2[remerge]+=Ep;

                    continue;// Dont include this as a viable cluster!
                }

                ++m_sss_num_candidates;

                //determine if this cluster is in neutrino slice
                m_sss_candidate_in_nu_slice.push_back(clu.InNuSlice(sliceIDToHitsMap, p_sliceid));

                //Fill All the bits
                m_sss_candidate_num_hits.push_back((int)hitz.size());
                m_sss_candidate_num_wires.push_back((int)ssscorz->n_wires);
                m_sss_candidate_num_ticks.push_back((int)ssscorz->n_ticks);
                m_sss_candidate_plane.push_back(pl);
                m_sss_candidate_PCA.push_back(ssscorz->pca_0);
                m_sss_candidate_impact_parameter.push_back(clu.getImpactParam());
                m_sss_candidate_fit_slope.push_back(clu.getFitSlope());
                m_sss_candidate_fit_constant.push_back(clu.getFitCons());
                m_sss_candidate_mean_tick.push_back(ssscorz->mean_tick);
                m_sss_candidate_max_tick.push_back(ssscorz->max_tick);
                m_sss_candidate_min_tick.push_back(ssscorz->min_tick);
                m_sss_candidate_min_wire.push_back(ssscorz->min_wire);
                m_sss_candidate_max_wire.push_back(ssscorz->max_wire);
                m_sss_candidate_mean_wire.push_back(ssscorz->mean_wire);
                m_sss_candidate_min_dist.push_back(ssscorz->min_dist);
                m_sss_candidate_wire_tick_based_length.push_back(clu.getWireTickBasedLength());
                m_sss_candidate_mean_ADC.push_back(clu.getMeanADC());
                m_sss_candidate_ADC_RMS.push_back(clu.getADCrms());
                m_sss_candidate_energy.push_back(Ep);
                m_sss_candidate_angle_to_shower.push_back(clu.getAngleWRTShower());
                m_sss_candidate_remerge.push_back(remerge);


                //MCTruth matching for pi0's
                if(m_is_data){
                    m_sss_candidate_matched.push_back(-1);
                    m_sss_candidate_pdg.push_back(-1);
                    m_sss_candidate_parent_pdg.push_back(-1);
                    m_sss_candidate_trackid.push_back(-1);
                    m_sss_candidate_true_energy.push_back(-1);
                    m_sss_candidate_overlay_fraction.push_back(-1);
                    m_sss_candidate_matched_energy_fraction_best_plane.push_back(-1);
                }else{

                    auto ssmatched = this->SecondShowerMatching(hitz, mcparticles_per_hit, mcParticleVector, pfParticleMap,  MCParticleToTrackIdMap);
                    m_sss_candidate_matched.push_back(ssmatched[0]);
                    m_sss_candidate_pdg.push_back(ssmatched[1]);
                    m_sss_candidate_parent_pdg.push_back(ssmatched[2]);
                    m_sss_candidate_trackid.push_back(ssmatched[3]);
                    m_sss_candidate_true_energy.push_back(ssmatched[4]);
                    m_sss_candidate_overlay_fraction.push_back(ssmatched[5]);
                    m_sss_candidate_matched_energy_fraction_best_plane.push_back(ssmatched[6]);

                    //Guanqun: print out (best-matched) truth information of the cluster
                    if(m_is_verbose){
                    std::cout << "Cluster: " << m_sss_num_candidates-1  << " plane: " << m_sss_candidate_plane.back() << ", energy: " << m_sss_candidate_energy.back() << "\n";
                    std::cout << "Cluster is matched: " << m_sss_candidate_matched.back() << ", matched PDG: " << m_sss_candidate_pdg.back() << " track ID: " << m_sss_candidate_trackid.back() << " overlay fraction: " << m_sss_candidate_overlay_fraction.back() << std::endl; 
                    std::cout << "===============================================================" << std::endl;
                    }
                }


                sevd.SetClusterLegend(c, m_sss_candidate_energy.back(),  m_sss_candidate_matched.back(), m_sss_candidate_pdg.back() , m_sss_candidate_overlay_fraction.back() );

            } //end of cluster loop

            // Plot the event
            if(m_SEAviewMakePDF){
                sevd.Print(m_SEAviewPlotDistance);
            }

        }

        for(int i =0; i<(int)showers.size(); i++){
            m_reco_shower_reclustered_energy_max[i] = std::max(m_reco_shower_reclustered_energy_plane1[i],std::max(m_reco_shower_reclustered_energy_plane0[i],m_reco_shower_reclustered_energy_plane2[i]));
        }

        // ################################################### END SEAview END SEAview #########################################################
        // #####################################################################################################################################



        // PandoraAllOutComes
        // I.e This runs over all 3D reco showers in the whole event and find second shower candidates
        if(!m_run_pi0_filter){
            std::cout<<"------------ Shower3D --------------"<<std::endl;
            /*for(auto &s : showers){
              std::cout<<"shower pfp key : "<<showerToNuPFParticleMap[s].key()<<" self: "<<showerToNuPFParticleMap[s]->Self()<<std::endl;
              }
              for(auto &s : tracks){
              std::cout<<"track pfp key : "<<trackToNuPFParticleMap[s].key()<<" self: "<<trackToNuPFParticleMap[s]->Self()<<std::endl;
              }*/

            this->SecondShowerSearch3D(showers, showerToNuPFParticleMap, tracks,trackToNuPFParticleMap,evt);


            //And cluster the 2d and 3d second showers. Very simple TODO
            this->SimpleSecondShowerCluster();

        }        


        // ################################################### Some info for Slice outside of Neutrino Slices #########################################################

        size_t n_neutrino_slice=0;
        size_t n_neutrino_candidate_pfp_id=0;

        for(size_t s=0; s< sliceVector.size(); s++){
            auto slice = sliceVector[s];
            std::vector<art::Ptr<recob::PFParticle>> pfps = sliceToPFParticlesMap[slice]; 

            int primaries=0;
            int n_dau=0;
            int found = 0;
            //std::cout<<"Starting a loop over "<<pfps.size()<<" pfparticles"<<std::endl;
            for(auto &pfp: pfps){
                //std::cout<<pfp->Self()<<" Primary: "<<pfp->IsPrimary()<<" PDG "<<pfp->PdgCode()<<" NDau: "<<pfp->NumDaughters()<<" Parent: "<<pfp->Parent()<<std::endl;

                if (!pfp->IsPrimary()) continue;
                // Check if this particle is identified as the neutrino
                const int pdg(pfp->PdgCode());
                const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
                primaries++;
                // If it is, lets get the vertex position
                if(isNeutrino){
                    found++;
                    //Ok this is neutrino candidate. 

                    std::cout<<"Found Neutrinoi Slice "<<s<<std::endl;
                    for(auto &pfp: pfps){
                        std::cout<<pfp->Self()<<" Primary: "<<pfp->IsPrimary()<<" PDG "<<pfp->PdgCode()<<" NDau: "<<pfp->NumDaughters()<<" Parent: "<<pfp->Parent()<<std::endl;
                    }
                    std::cout<<"************   Printing hierarcy "<<m_run_number<<" "<<m_subrun_number<<" "<<m_event_number<<" **************"<<std::endl;
                    n_neutrino_candidate_pfp_id = pfp->Self();
                    for (const size_t daughterId : pfp->Daughters()){
                        n_dau++;
                        auto dau = pfParticleMap[daughterId];
                        std::cout<<"---> gen1 --->"<<daughterId<<" trkScore: "<<PFPToTrackScoreMap[dau]<<" PDG: "<<dau->PdgCode()<<" NumDau: "<<dau->NumDaughters()<<std::endl;
                        auto tmp = dau;
                        int n_gen = 2;
                        for (const size_t granDaughterId : tmp->Daughters()){
                            while(tmp->NumDaughters()>0 && n_gen < 4){
                                for(int k=0; k< n_gen; k++){
                                    std::cout<<"---> ";
                                }
                                auto grandau = pfParticleMap[granDaughterId];
                                std::cout<<"gen"<<n_gen<<"  --->"<<granDaughterId<<" trkScore: "<<PFPToTrackScoreMap[grandau]<<" PDG: "<<grandau->PdgCode()<<" NumDau: "<<grandau->NumDaughters()<<std::endl;
                                tmp = grandau;    
                                n_gen++;
                            }
                            if(n_gen >=4) break;
                        }

                    }
                    std::cout<<"************   Finished hierarcy **************"<<std::endl;

                }
            }

            if(found==1){
                n_neutrino_slice = s;
                std::cout<<"Found a neutrino slice @ slice "<<n_neutrino_slice<<" ID "<<slice->ID()<<" key "<<slice.key()<<" pdfID "<<n_neutrino_candidate_pfp_id<<std::endl;
                std::cout<<"And there is "<<pfps.size()<<" PFParticles of which "<<primaries<<" are primary and "<<n_dau<<" are daughters of the Neutrino."<<std::endl;
                if((int)pfps.size() > n_dau+1){
                    std::cout<<"We're Missing Something!."<<std::endl;
                }
                m_reco_slice_objects = (int)pfps.size();
            }else if(found >1){
                throw cet::exception("DetachedVertexFinder") << "  This event contains multiple reconstructed neutrinos! Size: "<<found<<std::endl;
            }else if(found ==0){

            }
        }




        //---------------------- END OF LOOP, fill vertex ---------------------



        bool filter_pass_2g1p = Pi0PreselectionFilter();
        bool filter_pass_2g0p = Pi0PreselectionFilter2g0p();

        if (m_fill_trees &&  (  (filter_pass_2g1p && m_run_pi0_filter_2g1p) || (filter_pass_2g0p && m_run_pi0_filter_2g0p) || !m_run_pi0_filter ) ) {
            vertex_tree->Fill();
            //ncdelta_slice_tree->Fill();
            eventweight_tree->Fill();
            true_eventweight_tree->Fill();
            geant4_tree->Fill();
        }


        //Rest the vertex after filling
        this->ClearVertex();

        if(m_run_pi0_filter_2g1p)  return filter_pass_2g1p;
        else if(m_run_pi0_filter_2g0p)  return filter_pass_2g0p;

        //if not in filter mode pass all
        return true;

    }// end filter_module main class



    //-------------------------------------------------------------------------------------------
    void SinglePhoton::endJob()
    {
        if (m_print_out_event){
            out_stream.close();
        }
        pot_tree->Fill();
    }

    //-------------------------------------------------------------------------------------------

    //This runs ONCE at the start of the job and sets up all the necessary services and TTrees
    void SinglePhoton::beginJob()
    {
        mf::LogDebug("SinglePhoton") << " *** beginJob() *** " << "\n";
        art::ServiceHandle<art::TFileService> tfs;

        vertex_tree = tfs->make<TTree>("vertex_tree", "vertex_tree");
        pot_tree = tfs->make<TTree>("pot_tree", "pot_tree");
        eventweight_tree = tfs->make<TTree>("eventweight_tree", "eventweight_tree");
        ncdelta_slice_tree = tfs->make<TTree>("ncdelta_slice_tree", "ncdelta_slice_tree");
        run_subrun_tree = tfs->make<TTree>("run_subrun_tree","run_subrun_tree");
        geant4_tree = tfs->make<TTree>("geant4_tree","geant4_tree");

        //run_subrun_tree, reset some POT
        m_run = 0;
        m_subrun = 0;
        m_subrun_pot = 0;
        run_subrun_tree->Branch("run",&m_run,"run/I");
        run_subrun_tree->Branch("subrun",&m_subrun,"subrun/I");
        run_subrun_tree->Branch("subrun_pot",&m_subrun_pot,"subrun_pot/D");
        run_subrun_tree->Branch("subrun_counts",&m_subrun_counts,"subrun_counts/I");

        true_eventweight_tree = tfs->make<TTree>("true_eventweight_tree", "true_eventweight_tree");
        true_eventweight_tree->Branch("mcweight", "std::map<std::string, std::vector<double>>",&fmcweight);

        // --------------------- POT Releated variables -----------------
        m_number_of_events = 0;
        m_number_of_vertices = 0;
        m_pot_count=0;
        m_pot_per_event = 0;
        m_pot_per_subrun = 0;
        m_number_of_events_in_subrun=0;


        pot_tree->Branch("number_of_events",&m_number_of_events,"number_of_events/I");
        pot_tree->Branch("number_of_vertices",&m_number_of_vertices,"number_of_vertices/I");
        pot_tree->Branch("POT",&m_pot_count,"POT/D");

        // --------------------- Event Related variables ------------
        vertex_tree->Branch("run_number", &m_run_number, "run_number/I");
        vertex_tree->Branch("subrun_number", &m_subrun_number, "subrun_number/I");
        vertex_tree->Branch("event_number", &m_event_number, "event_number/I");

        vertex_tree->Branch("pot_per_event",&m_pot_per_event,"pot_per_event/D");
        vertex_tree->Branch("pot_per_subrun",&m_pot_per_subrun,"pot_per_subrun/D");
        vertex_tree->Branch("number_of_events_in_subrun",&m_number_of_events_in_subrun,"number_of_events_in_subrun/D");


        vertex_tree->Branch("genie_spline_weight", &m_genie_spline_weight, "genie_spline_weight/D");
        vertex_tree->Branch("genie_CV_tune_weight", &m_genie_CV_tune_weight, "genie_CV_tune_weight/D");

        vertex_tree->Branch("photonu_weight_low", &m_photonu_weight_low, "photonu_weight_low/D");
        vertex_tree->Branch("photonu_weight_high", &m_photonu_weight_high, "photonu_weight_high/D");

        vertex_tree->Branch("test_matched_hits", &m_test_matched_hits, "test_matched_hits/I");

        // --------------------- Vertex Related variables ------------
        vertex_tree->Branch("reco_vertex_size", &m_reco_vertex_size);
        vertex_tree->Branch("reco_vertex_x", &m_vertex_pos_x);
        vertex_tree->Branch("reco_vertex_y", &m_vertex_pos_y);
        vertex_tree->Branch("reco_vertex_z", &m_vertex_pos_z);
        vertex_tree->Branch("reco_vertex_wire_p0", &m_vertex_pos_wire_p0);
        vertex_tree->Branch("reco_vertex_wire_p1", &m_vertex_pos_wire_p1);
        vertex_tree->Branch("reco_vertex_wire_p2", &m_vertex_pos_wire_p2);
        vertex_tree->Branch("reco_vertex_tick", &m_vertex_pos_tick);


        vertex_tree->Branch("reco_vertex_in_SCB", &m_reco_vertex_in_SCB);
        vertex_tree->Branch("reco_vertex_dist_to_SCB",&m_reco_vertex_dist_to_SCB);
        vertex_tree->Branch("reco_vertex_dist_to_active_TPC",&m_reco_vertex_dist_to_active_TPC);
        vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane0",&m_reco_vertex_to_nearest_dead_wire_plane0);
        vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane1",&m_reco_vertex_to_nearest_dead_wire_plane1);
        vertex_tree->Branch("reco_vertex_to_nearest_dead_wire_plane2",&m_reco_vertex_to_nearest_dead_wire_plane2);

        vertex_tree->Branch("reco_slice_objects", &m_reco_slice_objects, "reco_slice_objects/I");

        vertex_tree->Branch("m_flash_optfltr_pe_beam",&m_flash_optfltr_pe_beam);
        vertex_tree->Branch("m_flash_optfltr_pe_veto",&m_flash_optfltr_pe_veto);
        vertex_tree->Branch("m_flash_optfltr_pe_veto_tot",&m_flash_optfltr_pe_veto_tot);
        vertex_tree->Branch("m_flash_optfltr_pe_beam_tot",&m_flash_optfltr_pe_beam_tot);

    
        vertex_tree->Branch("textgen_info",&m_textgen_info);
            
        //create branches as found in individual analyze_XXX.h
        this->CreateIsolationBranches();
        this->CreateSecondShowerBranches();
        this->CreateSecondShowerBranches3D();
        this->CreateStubBranches();
        this->CreateFlashBranches();
        this->CreateShowerBranches();
        this->CreateSliceBranches();
        this->CreateMCTruthBranches();
        this->CreateEventWeightBranches();
        this->CreateGeant4Branches();
        this->CreateTrackBranches();

        //hardcode some info (TODO change)
        //
        cet::search_path sp("FW_SEARCH_PATH");

        //Get the info for length->energy conversion from PSTAR database.
        TFile *fileconv;
        std::string proton_filename;
        sp.find_file("SinglePhotonAnalysis/proton_conversion.root",proton_filename);

        //some useful input data
        if(!m_run_pi0_filter){
            /*if(stat("proton_conversion.root", &buffer) == 0){
                fileconv = new TFile("proton_conversion.root", "read");
            }else{
                fileconv = new TFile((gpvm_location+"proton_conversion.root").c_str(), "read");
            }*/
            std::cout << "SinglePhoton \t||\t Loading proton length-> energy file from " <<proton_filename<< std::endl;

            fileconv = new TFile(proton_filename.c_str(),"read");

            proton_length2energy_tgraph = *(TGraph*)fileconv->Get("Graph");
            proton_length2energy_tgraph.GetMean();
            fileconv->Close();
        }

        //bad channels 
        std::string bad_channel_filename;// = "MCC9_channel_list.txt";
        sp.find_file("SinglePhotonAnalysis/MCC9_channel_list.txt",bad_channel_filename);

        if(!m_run_pi0_filter){
            /*if(stat(bad_channel_file.c_str(), &buffer) != 0){
                bad_channel_file = gpvm_location+bad_channel_file;
            }*/


            std::cout << "SinglePhoton \t||\t Loading channel id filename from "<<bad_channel_filename<< std::endl;
            std::ifstream bc_file(bad_channel_filename);

            if (bc_file.is_open())
            {
                std::string line;
                while ( getline (bc_file,line) )
                {
                    std::vector<int> res;
                    std::istringstream iss(line);
                    for(std::string s; iss >> s; )
                        res.push_back( std::stof(s));

                    std::pair<int,int> t(res[0],res[1]);
                    bad_channel_list_fixed_mcc9.push_back(t);
                }
                bc_file.close();
            }
        }


        //------------------- List of Selected Events to run --------
        if(m_runSelectedEvent){
            std::cout << "SinglePhoton \t||\t Running in selected-event only mode " << std::endl;

            std::ifstream infile(m_selected_event_list);
            if(!infile){
                std::cerr << "Fail to open file: " << m_selected_event_list << std::endl;
                return;
            }

            //read from file, run number, subrun number ,event number that should be run
            m_selected_set.clear();
            std::string line;
            while(std::getline(infile, line)){
                std::istringstream ss(line);

                std::vector<int> event_info; 
                for(int i; ss >> i; ) event_info.push_back(i);

                m_selected_set.insert(event_info);
            }

            infile.close();

            if(m_is_verbose){
                std::cout << "Selected Events: " << std::endl;
                std::cout << "Run \t SubRun \t Event" << std::endl;
                for(auto & v: m_selected_set){
                    std::for_each(v.begin(), v.end(), [](int n){std::cout << n<<" \t "; });
                    std::cout << std::endl;
                }
            }
        }

        std::cout<<"SinglePhoton \t||\t beginJob() is complete"<<std::endl;

    }




    //-------------------------------------------------------------------------------------------
    void SinglePhoton::ClearVertex(){
        //Clears and resets all output verticies in TTree

        //------------ Event related Variables -------------
        m_event_number = -99;
        m_subrun_number = -99;
        m_run_number = -99;
        m_test_matched_hits = 0;

        m_pot_per_event = 0;
        m_pot_per_subrun = m_subrun_pot;
        m_number_of_events_in_subrun = 0;

        m_genie_spline_weight = 1.0;

        m_textgen_info.clear();

        //------------ Vertex related Variables -------------
        m_reco_vertex_size = 0;
        m_vertex_pos_x=-99999;
        m_vertex_pos_y=-99999;
        m_vertex_pos_z=-99999;
        m_vertex_pos_tick=-9999;
        m_vertex_pos_wire_p0=-9999;
        m_vertex_pos_wire_p1=-9999;
        m_vertex_pos_wire_p2=-9999;
        m_reco_vertex_in_SCB = -9999;
        m_reco_vertex_dist_to_SCB = -9999;
        m_reco_vertex_dist_to_active_TPC= -9999;

        m_reco_vertex_to_nearest_dead_wire_plane0=-99999;
        m_reco_vertex_to_nearest_dead_wire_plane1=-99999;
        m_reco_vertex_to_nearest_dead_wire_plane2=-99999;

        m_reco_slice_objects = 0;

        this->ClearIsolation();
        this->ClearSecondShowers();
        this->ClearSecondShowers3D();
        this->ClearStubs();
        this->ClearFlashes();
        this->ClearTracks();
        this->ClearShowers();
        this->ClearMCTruths();
        this->ClearEventWeightBranches();
        fmcweight.clear();
        this->ClearGeant4Branches();
        this->ClearSlices();


    }


    bool SinglePhoton::beginSubRun(art::SubRun& sr) {


        m_run = sr.run();
        m_subrun = sr.subRun();

        double this_pot = 0;

        //reset subrun count 
        m_subrun_counts = 0;


        if(m_potLabel != ""){
            if(m_potLabel == "generator"){

                art::Handle<sumdata::POTSummary> gen_pot_hand;
                if(sr.getByLabel(m_potLabel,gen_pot_hand)){
                    this_pot =  gen_pot_hand->totgoodpot;
                    m_pot_count += this_pot;
                    std::cout<<"SinglePhoton::beginSubRun()\t||\t SubRun POT: "<<this_pot<<" . Current total POT this file: "<<m_pot_count<<" (label) "<<m_potLabel<<std::endl;
                }
            }else{

                art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
                if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
                    this_pot =potSummaryHandlebnbETOR875->totpot; 
                    m_pot_count += this_pot;
                    std::cout<<"SinglePhoton::beginSubRun()\t||\t SubRun POT: "<<potSummaryHandlebnbETOR875->totpot<<" . Current total POT this file: "<<m_pot_count<<" (label) "<<m_potLabel<<std::endl;
                }
            }
        }

        m_subrun_pot = this_pot; 

        return true;

    }


    bool SinglePhoton::endSubRun(art::SubRun&sr){


        run_subrun_tree->Fill();
        return true;
    }






    //-----------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::GetVertex(const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, const art::Ptr<recob::PFParticle> & particle ){

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Starting to analyze recob::Vertex\n";
        int n_vert =0;

        //std::cout<<"There are "<<pfParticlesToVerticesMap.count(particle)<<" verticies associated with this particle"<<std::endl;

        lar_pandora::PFParticlesToVertices::const_iterator vIter = pfParticlesToVerticesMap.find(particle);
        if (pfParticlesToVerticesMap.end() != vIter)
        {
            const lar_pandora::VertexVector &vertexVector = vIter->second;
            if (!vertexVector.empty())
            {
                if (vertexVector.size() !=1)
                    std::cout << " Warning: Found particle with more than one associated vertex " << "\n";

                const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                double xyz[3] = {0.0, 0.0, 0.0} ;
                vertex->XYZ(xyz);

                n_vert++;
                //std::cout<<"Vertex!"<<"\t "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n";

                m_vertex_pos_x = xyz[0];
                m_vertex_pos_y = xyz[1];
                m_vertex_pos_z = xyz[2];
                std::vector<double> tmp = {xyz[0],xyz[1],xyz[2]};
                m_reco_vertex_in_SCB = this->distToSCB(m_reco_vertex_dist_to_SCB,tmp);
                m_reco_vertex_dist_to_active_TPC = this->distToTPCActive(tmp);

                if(!m_run_pi0_filter){
                    m_reco_vertex_to_nearest_dead_wire_plane0 = distanceToNearestDeadWire(0, m_vertex_pos_y, m_vertex_pos_z,m_channelMap, bad_channel_list_fixed_mcc9);
                    m_reco_vertex_to_nearest_dead_wire_plane1 = distanceToNearestDeadWire(1, m_vertex_pos_y, m_vertex_pos_z,m_channelMap, bad_channel_list_fixed_mcc9);
                    m_reco_vertex_to_nearest_dead_wire_plane2 = distanceToNearestDeadWire(2, m_vertex_pos_y, m_vertex_pos_z,m_channelMap, bad_channel_list_fixed_mcc9);
                }

            }else{
                std::cout << " Error: vertexVector associated with this particle is empty " << "\n";
                std::cerr << " Error: vertexVector associated with this particle is empty " << "\n";
                //exit(0);

            }
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Finished. Found "<<n_vert<<" vertices.\n";
    }


    void SinglePhoton::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap)
    {
        //std::cout<<"Filling pfParticleMap with from the handle with total number "<<pfParticleHandle->size()<<std::endl;
        for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
        {
            const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
            // std::cout<<"Adding PFP to pfParticleMap with pfp id  "<<pParticle->Self()<<std::endl;
            if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second)
            {
                throw cet::exception("SinglePhoton") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
            }
        }
    }


    void SinglePhoton::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, PFParticleVector &crParticles, PFParticleVector &nuParticles )
    {

        int found = 0;
        int primaries = 0;
        int full = 0;
        for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
        {
            const art::Ptr<recob::PFParticle> pParticle(it->second);

            full++;
            // Only look for primary particles
            if (!pParticle->IsPrimary()) continue;

            // Check if this particle is identified as the neutrino
            const int pdg(pParticle->PdgCode());
            const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);


            primaries++;
            // If it is, lets get the vertex position
            if(isNeutrino){
                found++;
                this->GetVertex(pfParticlesToVerticesMap, pParticle );
                


            }

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
                throw cet::exception("SinglePhoton") << "  This event contains multiple reconstructed neutrinos!";
            }

            // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
            for (const size_t daughterId : pParticle->Daughters())
            {
                if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                    throw cet::exception("SinglePhoton") << "  Invalid PFParticle collection!";

                nuParticles.push_back(pfParticleMap.at(daughterId));
            }
        }
        std::cout<<"SinglePhoton::GetFinalStatePFParticleVectors()\t||\t Found "<<primaries<<" primary PFParticles (out of "<<full<<") of which: "<<found<<" were neutrinos."<<std::endl;
        m_reco_vertex_size = found;




    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::CollectTracksAndShowers(const PFParticleVector &particles,const PFParticleIdMap pfParticleMap, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
    {


        // Get the associations between PFParticles and tracks/showers from the event
        art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
        art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);

        //if running over the neutrino slice only 
        if (m_run_all_pfps == false){ 
            for (const art::Ptr<recob::PFParticle> &pParticle : particles) {
                const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
                const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));

                FillTracksAndShowers(associatedTracks, associatedShowers, pParticle,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);
            }
        } else{ //if running over all slices
            std::cout<<"SinglePhoton\t||\tThe total number of PFP's in the map is "<<pfParticleMap.size()<<std::endl;
            //            std::cout<<"The total number of PFP's in the vector is "<< particles.size()<<std::endl;
            for (auto pair : pfParticleMap){
                const art::Ptr<recob::PFParticle> &pParticle = pair.second;

                const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
                const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));

                FillTracksAndShowers(associatedTracks, associatedShowers, pParticle,  pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);

            }
        }


    }

    void SinglePhoton::FillTracksAndShowers( const std::vector< art::Ptr<recob::Track> > & associatedTracks, const std::vector< art::Ptr<recob::Shower> > & associatedShowers, const art::Ptr<recob::PFParticle> &pParticle , const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap)
    {

        const unsigned int nTracks(associatedTracks.size());
        const unsigned int nShowers(associatedShowers.size());


        // Check if the PFParticle has no associated tracks or showers
        if (nTracks == 0 && nShowers == 0)
        {
            //  std::cout<<"ERROR No tracks or showers were associated to PFParticle " << pParticle->Self()<<" with pdg "<<pParticle->PdgCode() <<std::endl;
            //std::cout<<"-- isPrimary = "<<pParticle->IsPrimary()<<std::endl;
            mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
            return;
        }

        // Check if there is an associated track
        if (nTracks == 1 && nShowers == 0)
        {

            tracks.push_back(associatedTracks.front());
            trackToNuPFParticleMap[tracks.back()]= pParticle;
            //std::cout<<"adding to trackToNuPFParticleMap this track with id "<<  associatedTracks.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

            return;
        }

        // Check if there is an associated shower
        if (nTracks == 0 && nShowers == 1)
        {
            showers.push_back(associatedShowers.front());
            showerToNuPFParticleMap[showers.back()] = pParticle;
            // std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  associatedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

            return;
        }

        // mcc9.10 quick fix for new pandora
        if (nTracks >0 && nShowers == 1)
        {
            showers.push_back(associatedShowers.front());
            showerToNuPFParticleMap[showers.back()] = pParticle;
            // std::cout<<"adding to showerToNuPFParticleMap this shower with id "<<  associatedShowers.front()->ID() << " and PFP "<< pParticle->Self()<<std::endl;

            return;
        }
        throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();

    }




    double SinglePhoton::triangle_area(double a1, double a2, double b1, double b2, double c1, double c2){
        double m1 = 0.3;
        double m2 = 1.0/25.0;

        return fabs((a1*m1*(b2*m2-c2*m2)+b1*m1*(c2*m2-a2*m2)+c1*m1*(a2*m2-b2*m2))/2.0);
    }

    int SinglePhoton::quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area){

        std::vector<double> z(n,0.0);

        TGraph2D *g = new TGraph2D(n,X,Y,&z[0]);
        TGraphDelaunay delan(g);
        delan.SetMarginBinsContent(0);
        delan.ComputeZ(0,0);
        delan.FindAllTriangles();
        (*num_triangles)=delan.GetNdt(); // number of Delaunay triangles found

        //Grab the locations of all the trianges. These will be intergers referencing to position in X,Y arrays
        Int_t *MT = delan.GetMTried();
        Int_t *NT = delan.GetNTried();
        Int_t *PT = delan.GetPTried();

        (*area)=0.0;
        for(int i = 0; i<delan.GetNdt(); i++){
            (*area)+=triangle_area(X[MT[i]-1],Y[MT[i]-1],X[NT[i]-1],Y[NT[i]-1],X[PT[i]-1],Y[PT[i]-1]);
        }

        delete g;
        return 0;
    }

    int SinglePhoton::delaunay_hit_wrapper(const std::vector<art::Ptr<recob::Hit>>& hits, std::vector<int> & num_hits, std::vector<int>& num_triangles, std::vector<double> & area){

        int n = hits.size();
        std::vector<double> C0,T0;
        std::vector<double> C1,T1;
        std::vector<double> C2,T2;
        size_t n_0=0;
        size_t n_1=0;
        size_t n_2=0;

        for(int i=0;i<n; i++){
            const art::Ptr<recob::Hit> hit = hits[i];
            switch(hit->View()){
                case 0:
                    C0.push_back((double)hit->WireID().Wire);         
                    T0.push_back(hit->PeakTime());         
                    n_0++;
                    break;
                case 1:
                    C1.push_back((double)hit->WireID().Wire);         
                    T1.push_back(hit->PeakTime());         
                    n_1++;
                    break;
                case 2:
                    C2.push_back((double)hit->WireID().Wire);         
                    T2.push_back(hit->PeakTime());         
                    n_2++;
                    break;
                default:
                    break;
            }
        }
        if(m_use_delaunay){
            if(n_0>0 && (int)n_0 < m_delaunay_max_hits) this->quick_delaunay_fit(n_0, &C0[0]  , &T0[0]  , &num_triangles[0],&area[0]);
            if(n_1>0 && (int)n_1 < m_delaunay_max_hits) this->quick_delaunay_fit(n_1, &C1[0]  , &T1[0]  , &num_triangles[1],&area[1]);
            if(n_2>0 && (int)n_2 < m_delaunay_max_hits) this->quick_delaunay_fit(n_2, &C2[0]  , &T2[0]  , &num_triangles[2],&area[2]);
        }
        num_hits[0] = n_0;
        num_hits[1] = n_1;
        num_hits[2] = n_2;

        //std::cout<<"Plane 0: "<<n_0<<" hits with "<<num_triangles[0]<<" triangles of area: "<< area[0]<<std::endl;
        //std::cout<<"Plane 1: "<<n_1<<" hits with "<<num_triangles[1]<<" triangles of area: "<< area[1]<<std::endl;
        //std::cout<<"Plane 2: "<<n_2<<" hits with "<<num_triangles[2]<<" triangles of area: "<< area[2]<<std::endl;

        return 0;
    }

    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected, std::vector<double> & input){
        corrected.resize(3);

        double kx = input[0];
        double ky = input[1];
        double kz = input[2];

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks.TPCG4Time2Tick(mcparticle->T())+theDetector.GetXTicksOffset(0,0,0)-trigger_offset(detClocks);

        double xtimeoffset = theDetector.ConvertTicksToX(g4Ticks,0,0,0);

        //        double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector.ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks.TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector.GetXTicksOffset(0,0,0)<<" "<<trigger_offset(detClocks)<<std::endl;
        return 0;
    }





    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);

        double kx = mcparticle->Vx();
        double ky = mcparticle->Vy();
        double kz = mcparticle->Vz();

        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks.TPCG4Time2Tick(mcparticle->T())+theDetector.GetXTicksOffset(0,0,0)-trigger_offset(detClocks);

        double xtimeoffset = theDetector.ConvertTicksToX(g4Ticks,0,0,0);

        //double xOffset = -scecorr.X() +xtimeoffset+0.6;
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=kx - scecorr.X() + xtimeoffset + 0.6; //due to sim/wirecell differences  Seev https://cdcvs.fnal.gov/redmine/projects/uboone-physics-analysis/wiki/MCC9_Tutorials 
        corrected[1]=ky+yOffset;
        corrected[2]=kz+zOffset;

        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<kx<<" "<<xOffset<<" "<<theDetector.ConvertTicksToX(g4Ticks, 0, 0, 0)<<" "<<scecorr.X()<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: "<<xOffset<<" "<<yOffset<<" "<<zOffset<<std::endl;
        //std::cout<<"SinglePhoton\t||\tTRIGGER_OFF: mcp->T(): "<<mcparticle->T()<<" TPCG4Time2Tick(): "<<detClocks.TPCG4Time2Tick(mcparticle->T())<<". "<<theDetector.GetXTicksOffset(0,0,0)<<" "<<trigger_offset(detClocks)<<std::endl;
        return 0;
    }





    int SinglePhoton::spacecharge_correction(const simb::MCParticle & mcparticle, std::vector<double> & corrected){
        corrected.resize(3);
        //Space Charge Effect! functionize this soon.
        double kx = mcparticle.Vx();
        double ky = mcparticle.Vy();
        double kz = mcparticle.Vz();
        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks.TPCG4Time2Tick(mcparticle.T())+theDetector.GetXTicksOffset(0,0,0)-trigger_offset(detClocks);

        double xtimeoffset = theDetector.ConvertTicksToX(g4Ticks,0,0,0);

        corrected[0]=kx - scecorr.X() +xtimeoffset+0.6;
        corrected[1]=ky + scecorr.Y();
        corrected[2]=kz + scecorr.Z();
        return 0;
    }

    void SinglePhoton::CollectMCParticles(const art::Event &evt, const std::string &label, std::map< art::Ptr<simb::MCTruth>, std::vector<art::Ptr<simb::MCParticle>>> &truthToParticles,        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>              &particlesToTruth, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap)
    {

        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

        art::Handle< std::vector< simb::MCParticle>  > theParticles;
        evt.getByLabel(label, theParticles);

        if (!theParticles.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find MC particles... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theParticles->size() << " MC particles " << std::endl;
        }

        art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);

        for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
        {
            const art::Ptr<simb::MCParticle> particle(theParticles, i);
            const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
            truthToParticles[truth].push_back(particle);
            particlesToTruth[particle] = truth;
            MCParticleToTrackIdMap[particle->TrackId()] = particle;
        }

        std::cout<<"SinglePhoton::CollectMCParticles() \t||\t the number of MCParticles in the event is "<<theParticles->size()<<std::endl;
    }

    void SinglePhoton::CollectSimChannels(const art::Event &evt, const std::string &label,  std::vector< art::Ptr<sim::SimChannel> >  &simChannelVector)
    {
        //    if (evt.isRealData())
        //      throw cet::exception("LArPandora") << " PandoraCollector::CollectSimChannels --- Trying to access MC truth from real data ";

        art::Handle< std::vector<sim::SimChannel> > theSimChannels;
        evt.getByLabel(label, theSimChannels);

        if (!theSimChannels.isValid())
        {
            mf::LogDebug("LArPandora") << "  Failed to find sim channels... " << std::endl;
            return;
        }
        else
        {
            mf::LogDebug("LArPandora") << "  Found: " << theSimChannels->size() << " SimChannels " << std::endl;
        }

        for (unsigned int i = 0; i < theSimChannels->size(); ++i)
        {
            const art::Ptr<sim::SimChannel> channel(theSimChannels, i);
            simChannelVector.push_back(channel);
        }
    }


    void SinglePhoton::BuildMCParticleHitMaps(const art::Event &evt, const std::string &label, const std::vector<art::Ptr<recob::Hit>> &hitVector,   std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  &particlesToHits,         std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  &hitsToParticles, const lar_pandora::LArPandoraHelper::DaughterMode daughterMode, std::map< int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap)
    {
        std::vector< art::Ptr<sim::SimChannel> >   simChannelVector;
        std::map< art::Ptr<simb::MCTruth>,     std::vector<art::Ptr<simb::MCParticle>>  >    truthToParticles;
        std::map< art::Ptr<simb::MCParticle>,  art::Ptr<simb::MCTruth> > particlesToTruth;
        std::map< art::Ptr<recob::Hit>,    std::vector< sim::TrackIDE >    >               hitsToTrackIDEs;

        this->CollectSimChannels(evt, label, simChannelVector);
        this->CollectMCParticles(evt, label, truthToParticles, particlesToTruth, MCParticleToTrackIdMap);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(evt, hitVector, simChannelVector, hitsToTrackIDEs);
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(hitsToTrackIDEs, truthToParticles, particlesToHits, hitsToParticles, daughterMode);


    }

    bool SinglePhoton::Pi0PreselectionFilter()
    {

        if(m_vertex_pos_x < 5.0 || m_vertex_pos_x > 251.) return false;
        if(m_vertex_pos_y < -112. || m_vertex_pos_y > 112.) return false;
        if(m_vertex_pos_z < 5.0 || m_vertex_pos_z > 1031.8) return false;

        if(m_reco_asso_showers!=2) return false;
        if(m_reco_asso_tracks!=1) return false;
        if(m_reco_vertex_size<1) return false;

        if(m_reco_shower_conversion_distance.size()!=2) return false;
        if(m_reco_shower_conversion_distance[0]<1. || m_reco_shower_conversion_distance[1]<1.) return false;

        return true;
    }



    bool SinglePhoton::Pi0PreselectionFilter2g0p()
    {

        if(m_vertex_pos_x < 5.0 || m_vertex_pos_x > 251.) return false;
        if(m_vertex_pos_y < -112. || m_vertex_pos_y > 112.) return false;
        if(m_vertex_pos_z < 5.0 || m_vertex_pos_z > 1031.8) return false;

        if(m_reco_asso_showers!=2) return false;
        if(m_reco_asso_tracks!=0) return false;
        if(m_reco_vertex_size<1) return false;

        if(m_reco_shower_energy_max.size()!=2) return false;
        //if the maximum energy of all showers on all planes is smaller than 30
        if(m_reco_shower_energy_max[m_reco_shower_ordered_energy_index[0]]<30.) return false;

        return true;
    }

    bool SinglePhoton::IsEventInList(int run, int subrun, int event){
        if(m_selected_set.find( {run, subrun, event} ) == m_selected_set.end()){
            if(m_selected_set.find({run, subrun})  == m_selected_set.end() ){
                if(m_selected_set.find({run}) == m_selected_set.end())
                    return false;
            }
        }
        return true;
    }

} //namespace

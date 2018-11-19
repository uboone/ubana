#include "SinglePhoton_module.h"
#include "analyze_OpFlashes.h"
#include "analyze_Tracks.h"
#include "analyze_Showers.h"
#include "analyze_Template.h"
#include "analyze_MCTruth.h"


#include "reco_truth_matching.h"
#include "bad_channel_matching.h"

namespace single_photon
{

    SinglePhoton::SinglePhoton(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
    {
        this->reconfigure(pset);
        theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
        detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
        SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        geom = lar::providerFrom<geo::Geometry>();

    }

    void SinglePhoton::reconfigure(fhicl::ParameterSet const &pset)
    {
        m_is_verbose = pset.get<bool>("Verbose",true);
        m_use_PID_algorithms = pset.get<bool>("usePID",false);

        m_pandoraLabel = pset.get<std::string>("PandoraLabel");
        m_trackLabel = pset.get<std::string>("TrackLabel");
        m_showerLabel = pset.get<std::string>("ShowerLabel");
        m_caloLabel = pset.get<std::string>("CaloLabel");
        m_flashLabel = pset.get<std::string>("FlashLabel");
        m_potLabel = pset.get<std::string>("POTLabel");
        m_beamgate_flash_start = pset.get<double>("beamgateStartTime",3.2); //Defaults to MC for now. Probably should change
        m_beamgate_flash_end = pset.get<double>("beamgateEndTime",4.8);
        m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
        m_badChannelLabel = pset.get<std::string>("BadChannelLabel","badmasks");
        m_badChannelProducer = pset.get<std::string>("BadChannelProducer","simnfspl1");

        m_generatorLabel = pset.get<std::string>("GeneratorLabel","generator");
        m_mcTrackLabel = pset.get<std::string>("MCTrackLabel","mcreco");
        m_mcShowerLabel = pset.get<std::string>("MCShowerLabel","mcreco");
        m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
        m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
        m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");


        //Some track calorimetry parameters
        m_track_calo_min_dEdx = pset.get<double>("Min_dEdx",0.01);
        m_track_calo_max_dEdx = pset.get<double>("Max_dEdx",25);
        m_track_calo_min_dEdx_hits = pset.get<double>("Min_dEdx_hits",2);
        m_track_calo_trunc_fraction = pset.get<double>("TruncMeanFraction",20.0);

        //Some shower calorimetry parameters
        m_work_function = pset.get<double>("work_function");
        m_recombination_factor =pset.get<double>("recombination_factor");
        m_gain =pset.get<double>("gain");
        m_wire_spacing = pset.get<double>("wire_spacing");
        m_width_dqdx_box = pset.get<double>("width_box");
        m_length_dqdx_box = pset.get<double>("length_box");
        m_pidLabel = pset.get<std::string>("ParticleIDLabel","particleid");
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::analyze(const art::Event &evt)
    {
        std::cout<<"---------------------------------------------------------------------------------"<<std::endl;
        std::cout<<"SinglePhoton::analyze()\t||\t On entry: "<<m_number_of_events<<std::endl;

        auto const TPC = (*geom).begin_TPC();
        auto ID = TPC.ID();
        m_Cryostat = ID.Cryostat;
        m_TPC = ID.TPC;

        this->ClearVertex();

        //******************************Setup*****************Setup**************************************/
        //***********************************************************************************************/
        // OK in this section we will get a bunch of stuff we need later, general associations and maps. These are either from pandora helper stuff or from direct associations. 
        // Make sure under the hood you understand this!
        // ------------------------
        // The basic idea is that here we get every possible vector of data products we might need, along with all maps. e.g 
        // tracks->pfparticles->hits
        // tracks->pfparticles->spacepoints ..etc..
        //
        // And then later we just write pseudo-independant code assuming you have every objecets you want (see analyze_Tracks.h) and assume you have access to everything. 


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

        //BadChannels
        art::Handle<std::vector<int> > badChannelHandle;
        evt.getByLabel(m_badChannelProducer, m_badChannelLabel, badChannelHandle);
        std::vector<int> badChannelVector = *(badChannelHandle);


        // Collect the PFParticles from the event. This is the core!
        art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
        std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
        art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

        //get the cluster handle for the dQ/dx calc
        art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pandoraLabel);
        std::vector< art::Ptr<recob::Cluster> > clusterVector;
        art::fill_ptr_vector(clusterVector,clusterHandle);

        //So a cross check
        if (!pfParticleHandle.isValid())
        {
            mf::LogDebug("SinglePhoton") << "  Failed to find the PFParticles.\n";
            return;
        }


        //Get the MCtruth handles and vectors
        art::ValidHandle<std::vector<simb::MCTruth>> const & mcTruthHandle  = evt.getValidHandle<std::vector<simb::MCTruth>>(m_generatorLabel);
        std::vector<art::Ptr<simb::MCTruth>> mcTruthVector;
        art::fill_ptr_vector(mcTruthVector,mcTruthHandle);


        //This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;
        // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
        PFParticleIdMap pfParticleMap;
        this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);


        //And some verticies.        
        art::ValidHandle<std::vector<recob::Vertex>> const & vertexHandle = evt.getValidHandle<std::vector<recob::Vertex>>(m_pandoraLabel);
        std::vector<art::Ptr<recob::Vertex>> vertexVector;
        art::fill_ptr_vector(vertexVector,vertexHandle);
        art::FindManyP<recob::Vertex> vertices_per_pfparticle(pfParticleHandle, evt, m_trackLabel);
        std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap;
        for(size_t i=0; i< pfParticleVector.size(); ++i){
            auto pfp = pfParticleVector[i];
            pfParticlesToVerticesMap[pfp] =vertices_per_pfparticle.at(pfp.key());
        }


        //Once we have actual verticies, lets concentrate on JUST the neutrino PFParticles for now:
        //--------------------------------
        // Produce two PFParticle vectors containing final-state particles:
        // 1. Particles identified as cosmic-rays - recontructed under cosmic-hypothesis
        // 2. Daughters of the neutrino PFParticle - reconstructed under the neutrino hypothesis
        std::vector< art::Ptr<recob::PFParticle> > crParticles;
        std::vector< art::Ptr<recob::PFParticle> > nuParticles;
        this->GetFinalStatePFParticleVectors(pfParticleMap, pfParticlesToVerticesMap, crParticles, nuParticles);


        //Look, here is a map that I just forced myself rather than build using helpers. Not that different is it. But for somereason I only use PFParticles.. huh,
        //Spacepoint associaitions
        art::FindManyP<recob::SpacePoint> spacePoints_per_pfparticle(pfParticleHandle, evt, m_trackLabel);
        std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>> > pfParticleToSpacePointsMap;
        for(size_t i=0; i< nuParticles.size(); ++i){
            const art::Ptr<recob::PFParticle> pfp = nuParticles[i];
            pfParticleToSpacePointsMap[pfp] = spacePoints_per_pfparticle.at(pfp.key());
        }

        //Get a map between the PFP's and the clusters. Although Mark isn't a fan of clusters, they're imporant for the shower dQ/dx
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


        //OK Here we build two IMPORTANT maps for the analysis, (a) given a PFParticle get a vector of hits..
        //and (b) given a single hit, get the PFParticle it is in (MARK: is it only one? always? RE-MARK: Yes)
        std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
        std::map<art::Ptr<recob::Hit>, art::Ptr<recob::PFParticle>>                hitToPFParticleMap;
        //Using a pandora helper here, but to be honest we should probably just build using normal associations so keep independant if pssoble
        lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, m_pandoraLabel, pfParticleToHitsMap, hitToPFParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters);


        //------------------------------------ OK, RECO - TRUTH matching stuff-----------------------------------
        //-------------------------------------------------------------------------------------------------------
        //Then build a map from MCparticles to Hits and vice versa
        std::map< art::Ptr<simb::MCParticle>,  std::vector<art::Ptr<recob::Hit> >  >  mcParticleToHitsMap;
        std::map< art::Ptr<recob::Hit>, art::Ptr<simb::MCParticle> >                  hitToMCParticleMap;
        lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(evt, m_geantModuleLabel, hitVector,  mcParticleToHitsMap, hitToMCParticleMap, lar_pandora::LArPandoraHelper::kAddDaughters);
        //This is basically just doing a FindMany..
        art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);


        //Apparrently a MCParticle doesn't know its origin (thanks Andy!)
        //I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa
        //Note which map is which!       //First  is one-to-many.         //Second is one-to-one
        std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;
        lar_pandora::LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap);
        //All of this will then come together to form Two things:



        art::ValidHandle<std::vector<sim::MCTrack>> const & mcTrackHandle  = evt.getValidHandle<std::vector<sim::MCTrack>>(m_mcTrackLabel);
        art::ValidHandle<std::vector<sim::MCShower>> const & mcShowerHandle  = evt.getValidHandle<std::vector<sim::MCShower>>(m_mcShowerLabel);
        std::vector<art::Ptr<sim::MCTrack>> mcTrackVector;
        std::vector<art::Ptr<sim::MCShower>> mcShowerVector;
        art::fill_ptr_vector(mcTrackVector,mcTrackHandle);
        art::fill_ptr_vector(mcShowerVector,mcShowerHandle);


        // These are the vectors to hold the tracks and showers for the final-states of the reconstructed neutrino
        //At this point, nuParticles is a std::vector< art::Ptr<recon::PFParticle>> of the PFParticles that we are interested in.
        //tracks is a vector of recob::Tracks and same for showers.
        //Implicitly, tracks.size() + showers.size() =  nuParticles.size(); At this point I would like two things.
        std::vector< art::Ptr<recob::Track> > tracks;
        std::vector< art::Ptr<recob::Shower> > showers;
        std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap; 
        std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap; 
        this->CollectTracksAndShowers(nuParticles, pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);

        //Track Calorimetry
        art::FindManyP<anab::Calorimetry> calo_per_track(trackHandle, evt, m_caloLabel);
        std::map<art::Ptr<recob::Track>, art::Ptr<anab::Calorimetry> > trackToCalorimetryMap;
        for(size_t i=0; i< tracks.size(); ++i){
                trackToCalorimetryMap[tracks[i]] = calo_per_track.at(tracks[i].key())[0];
        }

        art::FindOneP<anab::ParticleID> pid_per_track(trackHandle, evt, m_pidLabel);
        std::map<art::Ptr<recob::Track>, art::Ptr<anab::ParticleID> > trackToPIDMap;

        if(m_use_PID_algorithms){
            // Build a map to get PID from PFParticles, then call PID collection function
                     for(size_t i=0; i< tracks.size(); ++i){
                art::Ptr<recob::Track> track = tracks[i];
                trackToPIDMap[track] = pid_per_track.at(track.key());
            }
        }

        //**********************************************************************************************/
        //**********************************************************************************************/

        //Some event based properties

        m_number_of_events++;

        m_run_number = evt.run();
        m_event_number = evt.id().event();

        if(vertexVector.size()>0){
            m_number_of_vertices++;
        }


        //tests of reco_mc
        //Create a vectors and two maps to store the information
        std::vector<art::Ptr<simb::MCParticle>> mcParticleVector;
        std::map<art::Ptr<recob::Track>, art::Ptr<simb::MCParticle> > trackToMCParticleMap;
        std::map<art::Ptr<recob::Shower>, art::Ptr<simb::MCParticle> > showerToMCParticleMap;
        //Perform the (Templated for the recob::Objects)
        recoMCmatching<art::Ptr<recob::Track>>( tracks, trackToMCParticleMap, trackToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, mcParticleVector);
        recoMCmatching<art::Ptr<recob::Shower>>( showers, showerToMCParticleMap, showerToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit, mcParticleVector );

        //and now get the simb::MCparticle to both MCtrack and MCshower maps (just for the MCparticles matched ok).
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCTrack> > MCParticleToMCTrackMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<sim::MCShower> > MCParticleToMCShowerMap;

        std::cout<<"WAAH: track"<<std::endl;
        perfectRecoMatching<art::Ptr<sim::MCTrack>>(mcParticleVector, mcTrackVector, MCParticleToMCTrackMap);
        std::cout<<"WAAH: Shower"<<std::endl;
        perfectRecoMatching<art::Ptr<sim::MCShower>>(mcParticleVector, mcShowerVector, MCParticleToMCShowerMap);


        for(auto & track: tracks){
            auto mp = trackToMCParticleMap[track];
            //            auto mct = MCParticleToMCTrackMap[mp];
            //           auto mcs = MCParticleToMCTrackMap[mp];
            std::cout<<"CHECKTRACK: count trackmap: "<<MCParticleToMCTrackMap.count(mp)<<" "<< MCParticleToMCShowerMap.count(mp)<<std::endl;
        }
        for(auto & shower: showers){
            auto mp = showerToMCParticleMap[shower];
            //            auto mct = MCParticleToMCShowerMap[mp];
            //           auto mcs = MCParticleToMCShowerMap[mp];
            std::cout<<"CHECKSHOWER: count trackmap: "<<MCParticleToMCTrackMap.count(mp)<<" "<< MCParticleToMCShowerMap.count(mp)<<std::endl;
        }


        badChannelMatching<art::Ptr<recob::Track>>(badChannelVector, tracks, trackToNuPFParticleMap, pfParticleToHitsMap,geom,bad_channel_list_fixed_mcc9);



        if(m_is_verbose){
            std::cout << "SinglePhoton::analyze()\t||\t Consolidated event summary:" << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t - Number of neutrino final-state PFParticles : " << nuParticles.size() << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t    ... of which are track-like   : " << tracks.size() << "\n";
            std::cout << "SinglePhoton::analyze()\t||\t    ... of which are showers-like : " << showers.size() << "\n";
        }

        this->AnalyzeFlashes(flashVector);

        this->AnalyzeTracks(tracks, trackToNuPFParticleMap, pfParticleToSpacePointsMap);
        this->AnalyzeTrackCalo(tracks,   trackToCalorimetryMap);
        this->RecoMCTracks(tracks, trackToNuPFParticleMap, trackToMCParticleMap, MCParticleToMCTruthMap);
        if(m_use_PID_algorithms)  this->CollectPID(tracks, trackToPIDMap);


        this->AnalyzeShowers(showers,showerToNuPFParticleMap, pfParticleToHitsMap, pfParticleToClustersMap, clusterToHitsMap); 
        this->RecoMCShowers(showers, showerToNuPFParticleMap, showerToMCParticleMap, MCParticleToMCTruthMap);

        // MCTruth, MCParticle, MCNeutrino information all comes directly from GENIE.
        // MCShower and MCTrack come from energy depositions in GEANT4
        this->AnalyzeMCTruths(mcTruthVector);


        //---------------------- END OF LOOP, fill vertex ---------------------
        vertex_tree->Fill();

        std::cout<<"---------------------------------------------------------------------------------"<<std::endl;
    }



    //-------------------------------------------------------------------------------------------

    void SinglePhoton::endJob()
    {
        pot_tree->Fill();
    }

    //-------------------------------------------------------------------------------------------


    void SinglePhoton::beginJob()
    {
        mf::LogDebug("SinglePhoton") << " *** beginJob() *** " << "\n";

        art::ServiceHandle<art::TFileService> tfs;

        vertex_tree = tfs->make<TTree>("vertex_tree", "vertex_tree");
        pot_tree = tfs->make<TTree>("pot_tree", "pot_tree");

        // --------------------- POT Releated variables -----------------
        m_number_of_events = 0;
        m_number_of_vertices = 0;
        m_pot_count=0;
        pot_tree->Branch("number_of_events",&m_number_of_events,"number_of_events/I");
        pot_tree->Branch("number_of_vertices",&m_number_of_vertices,"number_of_vertices/I");
        pot_tree->Branch("POT",&m_pot_count,"POT/D");

        // --------------------- Event Related variables ------------
        vertex_tree->Branch("run_number", &m_run_number, "run_number/I");
        vertex_tree->Branch("event_number", &m_event_number, "event_number/I");

        // --------------------- Vertex Related variables ------------
        vertex_tree->Branch("reco_vertex_size", &m_reco_vertex_size);
        vertex_tree->Branch("reco_vertex_x", &m_vertex_pos_x);
        vertex_tree->Branch("reco_vertex_y", &m_vertex_pos_y);
        vertex_tree->Branch("reco_vertex_z", &m_vertex_pos_z);

        // --------------------- Flash Related Variables ----------------------
        this->CreateFlashBranches();

        // --------------------- Track Related variables ------------
        this->CreateTrackBranches();

        //Get the info for length->energy conversion from PSTAR database.
        TFile *fileconv = new TFile("/pnfs/uboone/resilient/users/markross/tars/proton_conversion.root", "read");
        proton_length2energy_tgraph = *(TGraph*)fileconv->Get("Graph");
        proton_length2energy_tgraph.GetMean();
        fileconv->Close();

        // --------------------- Shower Related variables ------------
        this->CreateShowerBranches();


        // ---------------------- MCTruth Related Variables ----------
        this->CreateMCTruthBranches();

        std::string bad_channel_file = "/pnfs/uboone/resilient/users/markross/tars/MCC9_channel_list.txt";
        std::ifstream bc_file(bad_channel_file);

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




    //-------------------------------------------------------------------------------------------
    void SinglePhoton::ClearVertex(){

        //------------ Event related Variables -------------
        m_event_number = 0;
        m_run_number = 0;


        //------------ Vertex related Variables -------------
        m_reco_vertex_size = 0;
        m_vertex_pos_x=0;
        m_vertex_pos_y=0;
        m_vertex_pos_z=0;

        //------------- Flash related Variables ------------------
        this->ClearFlashes();

        //------------- Track Related Variables -----------------
        this->ClearTracks();

        //------------- Track Related Variables -----------------
        this->ClearShowers();
        this->ClearMCTruths();



    }


    void SinglePhoton::beginSubRun(art::SubRun const & sr) {

        if(m_potLabel != ""){
            if(m_potLabel == "generator"){
                m_pot_count += sr.getValidHandle<sumdata::POTSummary>(m_potLabel)->totgoodpot;
            }else{
                art::Handle<sumdata::POTSummary> potSummaryHandlebnbETOR875;
                if (sr.getByLabel("beamdata","bnbETOR875",potSummaryHandlebnbETOR875)){
                    m_pot_count += potSummaryHandlebnbETOR875->totpot;
                }
            }
        }

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

            }else{
                std::cout << " Error: vertexVector associated with this particle is empty " << "\n";
                exit(0);

            }
        }

        if(m_is_verbose) std::cout<<"SinglePhoton::Getvertex()\t||\t Finished. Found "<<n_vert<<" vertices.\n";
    }


    void SinglePhoton::GetPFParticleIdMap(const PFParticleHandle &pfParticleHandle, PFParticleIdMap &pfParticleMap)
    {
        for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
        {
            const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
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

    void SinglePhoton::CollectTracksAndShowers(const PFParticleVector &particles, const PFParticleHandle &pfParticleHandle, const art::Event &evt, TrackVector &tracks, ShowerVector &showers,  std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle>>  &trackToNuPFParticleMap, std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> &showerToNuPFParticleMap )
    {


        // Get the associations between PFParticles and tracks/showers from the event
        art::FindManyP< recob::Track     > pfPartToTrackAssoc(pfParticleHandle, evt, m_trackLabel);
        art::FindManyP< recob::Shower    > pfPartToShowerAssoc(pfParticleHandle, evt, m_showerLabel);

        for (const art::Ptr<recob::PFParticle> &pParticle : particles)
        {
            const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
            const std::vector< art::Ptr<recob::Shower> > associatedShowers(pfPartToShowerAssoc.at(pParticle.key()));
            const unsigned int nTracks(associatedTracks.size());
            const unsigned int nShowers(associatedShowers.size());


            // Check if the PFParticle has no associated tracks or showers
            if (nTracks == 0 && nShowers == 0)
            {
                mf::LogDebug("SinglePhoton") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << "\n";
                continue;
            }

            // Check if there is an associated track
            if (nTracks == 1 && nShowers == 0)
            {

                tracks.push_back(associatedTracks.front());
                trackToNuPFParticleMap[tracks.back()]= pParticle;

                continue;
            }

            // Check if there is an associated shower
            if (nTracks == 0 && nShowers == 1)
            {
                showers.push_back(associatedShowers.front());
                showerToNuPFParticleMap[showers.back()] = pParticle;
                continue;
            }

            throw cet::exception("SinglePhoton") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
        }
    }




    double SinglePhoton::triangle_area(double a1, double a2, double b1, double b2, double c1, double c2){
        return fabs((a1*(b2-c2)+b1*(c2-a2)+c1*(a2-b2))/2.0);
    }

    int SinglePhoton::quick_delaunay_fit(int n, double *X, double *Y, int *num_triangles, double * area){

        std::vector<double> z(n,0.0);

        TGraph2D *g = new TGraph2D(n,X,Y,&z[0]);
        TGraphDelaunay delan(g);
        delan.SetMarginBinsContent(0);
        delan.ComputeZ(0,0);
        delan.FindAllTriangles();
        (*num_triangles)=delan.GetNdt();

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
                    C0.push_back((double)hit->Channel());         
                    T0.push_back(hit->PeakTime());         
                    n_0++;
                    break;
                case 1:
                    C1.push_back((double)hit->Channel());         
                    T1.push_back(hit->PeakTime());         
                    n_1++;
                    break;
                case 2:
                    C2.push_back((double)hit->Channel());         
                    T2.push_back(hit->PeakTime());         
                    n_2++;
                    break;
                default:
                    break;
            }
        }
        if(n_0>0) this->quick_delaunay_fit(n_0, &C0[0]  , &T0[0]  , &num_triangles[0],&area[0]);
        if(n_1>0) this->quick_delaunay_fit(n_1, &C1[0]  , &T1[0]  , &num_triangles[1],&area[1]);
        if(n_2>0) this->quick_delaunay_fit(n_2, &C2[0]  , &T2[0]  , &num_triangles[2],&area[2]);

        num_hits[0] = n_0;
        num_hits[1] = n_1;
        num_hits[2] = n_2;

        //std::cout<<"Plane 0: "<<n_0<<" hits with "<<num_triangles[0]<<" triangles of area: "<< area[0]<<std::endl;
        //std::cout<<"Plane 1: "<<n_1<<" hits with "<<num_triangles[1]<<" triangles of area: "<< area[1]<<std::endl;
        //std::cout<<"Plane 2: "<<n_2<<" hits with "<<num_triangles[2]<<" triangles of area: "<< area[2]<<std::endl;

        return 0;
    }


    int SinglePhoton::spacecharge_correction(const art::Ptr<simb::MCParticle> & mcparticle, std::vector<double> & corrected){

        corrected.resize(3);
        //Space Charge Effect! functionize this soon.
        double kx = mcparticle->Position().X();
        double ky = mcparticle->Position().Y();
        double kz = mcparticle->Position().Z();
        auto scecorr = SCE->GetPosOffsets( geo::Point_t(kx,ky,kz));
        double g4Ticks = detClocks->TPCG4Time2Tick(mcparticle->T())+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
        double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
        double yOffset = scecorr.Y();
        double zOffset = scecorr.Z();

        corrected[0]=xOffset;
        corrected[1]=yOffset;
        corrected[2]=zOffset;

        return 0;
    }



} //namespace lar_pandora

#include "SinglePhoton_module.h"
#include "analyze_OpFlashes.h"
#include "analyze_Tracks.h"
#include "analyze_Showers.h"
#include "analyze_Template.h"



#include "reco_truth_matching.h"
namespace single_photon
{

    SinglePhoton::SinglePhoton(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset)
    {
        this->reconfigure(pset);
        theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
	geom = lar::providerFrom<geo::Geometry>();
    
   }
//  -------------------------------------------------------------------------------------------------------

    void SinglePhoton::reconfigure(fhicl::ParameterSet const &pset)
    {
        m_pandoraLabel = pset.get<std::string>("PandoraLabel");
        m_trackLabel = pset.get<std::string>("TrackLabel");
        m_showerLabel = pset.get<std::string>("ShowerLabel");
        m_caloLabel = pset.get<std::string>("CaloLabel");
        m_printOutScores = pset.get<bool>("PrintOutScores",true);
        m_flashLabel = pset.get<std::string>("FlashLabel");
        m_beamgate_flash_start = pset.get<double>("beamgateStartTime",3.2); //Defaults to MC for now. Probably should change
        m_beamgate_flash_end = pset.get<double>("beamgateEndTime",4.8);
        m_potLabel = pset.get<std::string>("POTLabel");
        m_geantModuleLabel = pset.get<std::string>("GeantModule","largeant");
        m_backtrackerLabel = pset.get<std::string>("BackTrackerModule","gaushitTruthMatch");
        m_hitfinderLabel = pset.get<std::string>("HitFinderModule", "gaushit");
        m_hitMCParticleAssnsLabel = pset.get<std::string>("HitMCParticleAssnLabel","gaushitTruthMatch");
        m_useModBox = pset.get<bool>("UseModBox",true);
        m_is_verbose = pset.get<bool>("Verbose",true);
	m_work_function = pset.get<double>("work_function");
	m_recombination_factor =pset.get<double>("recombination_factor");
	m_gain =pset.get<double>("gain");

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void SinglePhoton::analyze(const art::Event &evt)
    {
        std::cout<<"---------------------------------------------------------------------------------"<<std::endl;
        std::cout<<"SinglePhoton::analyze()\t||\t On entry: "<<m_number_of_events<<std::endl;

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

        //This is another pandora helper. I don't like PFParticle ID lookups but I guess lets keep for now;
        // Produce a map of the PFParticle IDs for fast navigation through the hierarchy
        PFParticleIdMap pfParticleMap;
        this->GetPFParticleIdMap(pfParticleHandle, pfParticleMap);


        //And some verticies. I am pretty bad at consistency here, lets just abandon helper soon.        
        std::vector<art::Ptr<recob::Vertex>> vertexVector;
        std::map< art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Vertex>> > pfParticlesToVerticesMap;
        lar_pandora::LArPandoraHelper::CollectVertices(evt, m_pandoraLabel, vertexVector, pfParticlesToVerticesMap);


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
            //const art::Ptr<recob::PFParticle> pfp = nuParticles[i];
            auto pfp = nuParticles[i];
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
        art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> mcparticles_per_hit(hitHandle, evt, m_hitMCParticleAssnsLabel);

        //Apparrently a MCParticle doesn't know its origin (thanks Andy!)
        //I would also like a map from MCparticle to MCtruth and then I will be done.  and Vice Versa
        //Note which map is which!       //First  is one-to-many.         //Second is one-to-one
        std::map< art::Ptr<simb::MCTruth>,    std::vector<art::Ptr<simb::MCParticle>>>  MCTruthToMCParticlesMap;
        std::map< art::Ptr<simb::MCParticle>, art::Ptr<simb::MCTruth>>                  MCParticleToMCTruthMap;
        lar_pandora::LArPandoraHelper::CollectMCParticles(evt, m_geantModuleLabel, MCTruthToMCParticlesMap, MCParticleToMCTruthMap);
        //All of this will then come together to form Two things:


        // These are the vectors to hold the tracks and showers for the final-states of the reconstructed neutrino
        //At this point, nuParticles is a std::vector< art::Ptr<recon::PFParticle>> of the PFParticles that we are interested in.
        //tracks is a vector of recob::Tracks and same for showers.
        //Implicitly, tracks.size() + showers.size() =  nuParticles.size(); At this point I would like two things.
        std::vector< art::Ptr<recob::Track> > tracks;
        std::vector< art::Ptr<recob::Shower> > showers;
        std::map< art::Ptr<recob::Track> , art::Ptr<recob::PFParticle >> trackToNuPFParticleMap; 
        std::map< art::Ptr<recob::Shower> , art::Ptr<recob::PFParticle>> showerToNuPFParticleMap; 
        this->CollectTracksAndShowers(nuParticles, pfParticleHandle, evt, tracks, showers, trackToNuPFParticleMap, showerToNuPFParticleMap);


        //**********************************************************************************************/
        //**********************************************************************************************/

        //Some event based properties

        m_number_of_events++;
        this->ClearVertex();

        m_run_number = evt.run();
        m_event_number = evt.id().event();

        if(vertexVector.size()>0){
            m_number_of_vertices++;
        }


        //test of reco_mc

        recoMCmatching<art::Ptr<recob::Track>>( tracks, trackToNuPFParticleMap, pfParticleToHitsMap, mcparticles_per_hit );


        if(m_is_verbose){
        std::cout << "SinglePhoton::analyze()\t||\t Consolidated event summary:" << "\n";
        std::cout << "SinglePhoton::analyze()\t||\t - Number of primary cosmic-ray PFParticles   : " << crParticles.size() << "\n";
        std::cout << "SinglePhoton::analyze()\t||\t - Number of neutrino final-state PFParticles : " << nuParticles.size() << "\n";
        std::cout << "SinglePhoton::analyze()\t||\t    ... of which are track-like   : " << tracks.size() << "\n";
        std::cout << "SinglePhoton::analyze()\t||\t    ... of which are showers-like : " << showers.size() << "\n";
        }
        
        this->AnalyzeFlashes(flashVector);
        this->AnalyzeTracks(tracks, trackToNuPFParticleMap, pfParticleToSpacePointsMap);
        this->AnalyzeShowers(showers,showerToNuPFParticleMap, pfParticleToHitsMap, pfParticleToClustersMap, clusterToHitsMap); 

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
    vertex_tree->Branch("reco_nuvertx", &m_vertex_pos_x, "reco_nuvertx/D");
    vertex_tree->Branch("reco_nuverty", &m_vertex_pos_y, "reco_nuverty/D");
    vertex_tree->Branch("reco_nuvertz", &m_vertex_pos_z, "reco_nuvertz/D");

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

}




//-------------------------------------------------------------------------------------------
void SinglePhoton::ClearVertex(){

    //------------ Event related Variables -------------
    m_event_number = 0;
    m_run_number = 0;


    //------------ Vertex related Variables -------------
    m_vertex_pos_x=0;
    m_vertex_pos_y=0;
    m_vertex_pos_z=0;

    //------------- Flash related Variables ------------------
    this->ClearFlashes();

    //------------- Track Related Variables -----------------
    this->ClearTracks();

    //------------- Track Related Variables -----------------
    this->ClearShowers();




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
            //std::cout<<"Vertex!"<<"\n";
            //std::cout<<"\t "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<"\n";

            m_vertex_pos_x = xyz[0];
            m_vertex_pos_y = xyz[1];
            m_vertex_pos_z = xyz[2];


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

//------------------------------------------------------------------------------------------------------------------------------------------

/*     ---------- Currently commeted out as not in base larpandora
void SinglePhoton::PrintOutScores(const art::Event &evt, const PFParticleHandle &pfParticleHandle) const
{
    // Get the associations between PFParticles and larpandoraobj::PFParticleMetadata
    art::FindManyP< larpandoraobj::PFParticleMetadata > pfPartToMetadataAssoc(pfParticleHandle, evt, m_pandoraLabel);

    for (unsigned int i = 0; i < pfParticleHandle->size(); ++i)
    {
        const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > &pfParticleMetadataList(pfPartToMetadataAssoc.at(i));
        if (!pfParticleMetadataList.empty())
        {
            const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
            for (unsigned int j=0; j<pfParticleMetadataList.size(); ++j)
            {
                const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
                const pandora::PropertiesMap &pfParticlePropertiesMap(pfParticleMetadata->GetPropertiesMap());
                if (!pfParticlePropertiesMap.empty())
                    std::cout << " Found PFParticle " << pParticle->Self() << " with: " << "\n";
                for (pandora::PropertiesMap::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
                    std::cout << "  - " << it->first << " = " << it->second << "\n";
            }
        }
    }
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

void SinglePhoton::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, const lar_pandora::PFParticlesToVertices &pfParticlesToVerticesMap, PFParticleVector &crParticles, PFParticleVector &nuParticles )
{

    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)
    {
        const art::Ptr<recob::PFParticle> pParticle(it->second);

        // Only look for primary particles
        if (!pParticle->IsPrimary()) continue;

        // Check if this particle is identified as the neutrino
        const int pdg(pParticle->PdgCode());
        const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);


        // If it is, lets get the vertex position
        if(isNeutrino){
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


} //namespace lar_pandora

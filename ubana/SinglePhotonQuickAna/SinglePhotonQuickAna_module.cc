////////////////////////////////////////////////////////////////////////
// Class:       SinglePhotonQuickAna
// Plugin Type: analyzer (art v3_01_02)
// File:        SinglePhotonQuickAna_module.cc
//
// Generated at Mon Jan 13 05:58:20 2020 by Mark Ross-Lonergan using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/Wire.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h" 


#include <vector>

class SinglePhotonQuickAna;


class SinglePhotonQuickAna : public art::EDAnalyzer {
    public:
        explicit SinglePhotonQuickAna(fhicl::ParameterSet const& p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        SinglePhotonQuickAna(SinglePhotonQuickAna const&) = delete;
        SinglePhotonQuickAna(SinglePhotonQuickAna&&) = delete;
        SinglePhotonQuickAna& operator=(SinglePhotonQuickAna const&) = delete;
        SinglePhotonQuickAna& operator=(SinglePhotonQuickAna&&) = delete;

        // Required functions.
        void analyze(art::Event const& e) override;

        // Selected optional functions.
        void beginJob() override;
        void endJob() override;

    private:
        // Declare member data here.
        // All of these below are loaded from fcl via pset
        std::string m_pandoraLabel;
        std::string m_showerLabel;
        std::string m_showerKalmanLabel;
        std::string m_showerKalmanCaloLabel;
        double m_work_function;  
        double m_recombination_factor; 
        std::vector<double> m_gain_mc; 
        std::vector<double> m_gain_data; 
        double m_wire_spacing; 
        double m_width_dqdx_box; 
        double m_length_dqdx_box;
        double m_time2cm;

        bool m_is_data;
        bool m_is_overlayed;

        detinfo::DetectorProperties const * theDetector ;
        detinfo::DetectorClocks    const *  detClocks   ;
        spacecharge::SpaceCharge const * SCE;


        /*
         *@brief Calculated the shower energy by looping over all the hits and summing the charge
         *@param hits -  an art pointer of all the hits in a shower
         * */
        double CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits); /* returns max energy on any of the planes */
        double CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int plane); /* returns energy sum of hits on a certain plane */
        int getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);
        double getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane);
        /**
         *@brief Takes a hit and multiplies the charge by the gain
         *@param thishitptr art pointer to a hit
         *@param plane the plane the hit is on
         **/
        double GetQHit(art::Ptr<recob::Hit> thishitptr, int plane);
        
        /**
         * @brief Calculate the E value in MeV for a given hit
         * @param thishit - an individual hit 
         * */
        double QtoEConversionHit(art::Ptr<recob::Hit> thishitptr, int plane);

        /**
         * @brief Calculate the E value in MeV from a given Q value
         * @param q - the charge value
         * */
        double QtoEConversion(double q);

        /**
         *@brief Takes a vector of dQ/dx values and converts to dE/dx
         *@param dqdx - vector of dqdx points
         * */
        std::vector<double> CalcdEdxFromdQdx(std::vector<double> dqdx);

        /**
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
        double getPitch(TVector3 shower_dir, int plane);  /* distance between hit in shower direction projected on plane */
        TVector3 getWireVec(int plane); /* unit vector orthogonal to the  wire direction of plane -- usually named as wire_dir */
        double getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir); /* dot product of wire_dir and shower direction vectors */

        double degToRad(double deg);
        double radToDeg(double rad);
        /**
         *@brief Calculates the four corners of a rectangle of given length and width around a cluster given the start point and axis direction
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
         *@brief For a 2d point on a plane in cm and a rectangle, returns true if ponint is inside or on the boundary
         *uses triangle area check
         * */
        bool isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle);
        double areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3);
        /***
         *@brief returns the value at the median position in a vector of doubles, returns nan for vector of size <= 0
         *@param thisvector - vector of doubles
         *
         * */
        double getMedian(std::vector<double> thisvector);


};





SinglePhotonQuickAna::SinglePhotonQuickAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}  // ,
{
    m_pandoraLabel = p.get<std::string>("PandoraLabel");
    m_showerLabel = p.get<std::string>("ShowerLabel");

    m_is_data = p.get<bool>("isData");
    m_is_overlayed = p.get<bool>("isOverlayed");

    m_showerKalmanLabel = p.get<std::string>("ShowerTrackFitter","pandoraKalmanShower");
    m_showerKalmanCaloLabel =  p.get<std::string>("ShowerTrackFitterCalo","pandoraKalmanShowercali");

    m_work_function = p.get<double>("work_function");
    m_recombination_factor =p.get<double>("recombination_factor");
    m_gain_mc =p.get<std::vector<double>>("gain_mc");
    m_gain_data =p.get<std::vector<double>>("gain_data");
    m_wire_spacing = p.get<double>("wire_spacing");
    m_width_dqdx_box = p.get<double>("width_box");
    m_length_dqdx_box = p.get<double>("length_box");

    theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
    m_time2cm = theDetector->SamplingRate() / 1000.0 * theDetector->DriftVelocity( theDetector->Efield(), theDetector->Temperature() );

    detClocks   = lar::providerFrom<detinfo::DetectorClocksService>();
    SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

}




void SinglePhotonQuickAna::analyze(art::Event const& evt)
{

    //Collect the PFParticles from the event. This is the core!
    art::ValidHandle<std::vector<recob::PFParticle>> const & pfParticleHandle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pandoraLabel);
    std::vector<art::Ptr<recob::PFParticle>> pfParticleVector;
    art::fill_ptr_vector(pfParticleVector,pfParticleHandle);

    //So a cross check
    if (!pfParticleHandle.isValid())
    {
        mf::LogDebug("SinglePhotonQuickEnergy") << "  Failed to find the PFParticles.\n";
        return;
    }

    //get the cluster handle for direct hit based dQ/dx calc
    art::ValidHandle<std::vector<recob::Cluster>> const & clusterHandle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pandoraLabel);
    std::vector< art::Ptr<recob::Cluster> > clusterVector;
    art::fill_ptr_vector(clusterVector,clusterHandle);


    //Grab all recob::Showers in the event
    art::ValidHandle<std::vector<recob::Shower>> const & ShowerHandle = evt.getValidHandle<std::vector<recob::Shower>>(m_showerLabel);
    std::vector<art::Ptr<recob::Shower>> ShowerVector;
    art::fill_ptr_vector(ShowerVector,ShowerHandle);


    //---------Build up a map of PFParticles to RecobShowers
    art::FindOneP<recob::Shower> shower_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Shower>> pfParticlesToShowerMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];
        if(!shower_per_pfparticle.at(pfp.key()).isNull()){ 
            pfParticlesToShowerMap[pfp] =shower_per_pfparticle.at(pfp.key());
        }
    }


    //Get a map between the PFP's and the clusters  they're imporant for the shower dQ/dx
    //Also need a map between clusters and hits
    art::FindManyP<recob::Cluster> clusters_per_pfparticle(pfParticleHandle, evt, m_pandoraLabel);
    art::FindManyP<recob::Hit> hits_per_cluster(clusterHandle, evt, m_pandoraLabel);
    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Cluster>> > pfParticleToClustersMap;
    std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> > clusterToHitsMap;
    //fill map PFP to Clusters
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];
        pfParticleToClustersMap[pfp] = clusters_per_pfparticle.at(pfp.key());
    }
    //fill map Cluster to Hits
    for(size_t i=0; i< clusterVector.size(); ++i){
        auto cluster = clusterVector[i];
        clusterToHitsMap[cluster] = hits_per_cluster.at(cluster.key());
    }

    std::map<art::Ptr<recob::PFParticle>,  std::vector<art::Ptr<recob::Hit>> > pfParticleToHitsMap;
    //use pfp->cluster and cluster->hit to build pfp->hit map
    //for each PFP
    for(size_t i=0; i<pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];
        //get the associated clusters
        std::vector<art::Ptr<recob::Cluster>> clusters_vec  = pfParticleToClustersMap[pfp] ;
        //make empty vector to store hits
        std::vector<art::Ptr<recob::Hit>> hits_for_pfp = {};
        //for each cluster, get the associated hits
        for (art::Ptr<recob::Cluster> cluster: clusters_vec){
            std::vector<art::Ptr<recob::Hit>> hits_vec =  clusterToHitsMap[cluster];
            //   std::cout<<"looking at cluster in pfp "<<pfp->Self()<<" with "<<hits_vec.size() <<" hits"<<std::endl;
            //insert hits into vector
            hits_for_pfp.insert( hits_for_pfp.end(), hits_vec.begin(), hits_vec.end() );
        }
        //fill the map
        pfParticleToHitsMap[pfp] = hits_for_pfp;
        //std::cout<<"saving a total of "<<hits_for_pfp.size()<<" hits for pfp "<<pfp->Self()<<std::endl;
    }//for each pfp





    //**************************** Starting here to loop over all showers ****************************//
    int i_shr = 0;

    for (std::vector< art::Ptr<recob::PFParticle> >::const_iterator iter = pfParticleVector.begin(), iterEnd = pfParticleVector.end(); iter != iterEnd; ++iter)
    {

        //Lets get all the information we need at the start. The PFparticle, the shower (if) associated, all the hits and all the clusters in the shower. 
        const art::Ptr<recob::PFParticle>           pfp = *iter;
        if( pfParticlesToShowerMap.count(pfp)==0) continue;
        const art::Ptr<recob::Shower>               shower = pfParticlesToShowerMap[pfp];
        const std::vector<art::Ptr<recob::Hit>>     hits =  pfParticleToHitsMap[pfp];
        const std::vector<art::Ptr<recob::Cluster>> clusters = pfParticleToClustersMap[pfp];


        //------------- calorimetry ------------

        //double reco_shower_energy_max = CalcEShower(hits);
        double reco_shower_energy_plane0 = CalcEShowerPlane(hits, 0);
        double reco_shower_energy_plane1 = CalcEShowerPlane(hits, 1);
        double reco_shower_energy_plane2 = CalcEShowerPlane(hits, 2);

        int reco_shower_plane0_nhits = getNHitsPlane(hits, 0);
        int reco_shower_plane1_nhits = getNHitsPlane(hits, 1);
        int reco_shower_plane2_nhits = getNHitsPlane(hits, 2);

        std::vector<double> reco_shower_dQdx_plane0 = CalcdQdxShower(shower,clusters, clusterToHitsMap, 0 ); 
        std::vector<double> reco_shower_dQdx_plane1 = CalcdQdxShower(shower,clusters, clusterToHitsMap, 1 ); 
        std::vector<double> reco_shower_dQdx_plane2 = CalcdQdxShower(shower,clusters, clusterToHitsMap, 2 ); 

        std::vector<double> reco_shower_dEdx_plane0 = CalcdEdxFromdQdx(reco_shower_dQdx_plane0);
        std::vector<double> reco_shower_dEdx_plane1 = CalcdEdxFromdQdx(reco_shower_dQdx_plane1);
        std::vector<double> reco_shower_dEdx_plane2 = CalcdEdxFromdQdx(reco_shower_dQdx_plane2);

        double reco_shower_dEdx_plane0_median = getMedian(reco_shower_dEdx_plane0);
        double reco_shower_dEdx_plane1_median = getMedian(reco_shower_dEdx_plane1);
        double reco_shower_dEdx_plane2_median = getMedian(reco_shower_dEdx_plane2);

        std::cout<<"Shower # "<<i_shr<<" Energy : ( "<<reco_shower_energy_plane0<<","<<reco_shower_energy_plane1<<","<<reco_shower_energy_plane2<<"), Nhits : ( "<<reco_shower_plane0_nhits<<","<<reco_shower_plane1_nhits<<","<<reco_shower_plane2_nhits<<"), Median dEdx : ("<<reco_shower_dEdx_plane0_median<<","<<reco_shower_dEdx_plane1_median<<","<<reco_shower_dEdx_plane2_median<<")"<<std::endl;
        i_shr++;

    }//end Shower loop





    //**************************** On to Kalman Fitted Shower Trunks ****************************//


    //---------Kalman Track Showers
    art::FindOneP<recob::Track> showerKalman_per_pfparticle(pfParticleHandle, evt, m_showerKalmanLabel);
    std::map<art::Ptr<recob::PFParticle>, art::Ptr<recob::Track>> pfParticlesToShowerKalmanMap;
    for(size_t i=0; i< pfParticleVector.size(); ++i){
        auto pfp = pfParticleVector[i];
        if(!showerKalman_per_pfparticle.at(pfp.key()).isNull()){ 
            pfParticlesToShowerKalmanMap[pfp] =showerKalman_per_pfparticle.at(pfp.key());
        }
    }


    //----- kalman Cali (calibrated calorimetry)
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


}

void SinglePhotonQuickAna::beginJob()
{
    // Implementation of optional member function here.
}

void SinglePhotonQuickAna::endJob()
{
    // Implementation of optional member function here.
}




//Below here a lot of functions for calculating things

double SinglePhotonQuickAna::CalcEShowerPlane(const std::vector<art::Ptr<recob::Hit>>& hits, int this_plane){    
    double energy = 0.;

    //for each hit in the shower
    for (auto &thishitptr : hits){
        //check the plane
        int plane= thishitptr->View();

        //skip invalid planes     	
        if (plane != this_plane )	continue;

        //calc the energy of the hit
        double E = QtoEConversionHit(thishitptr, plane);	

        //add the energy to the plane
        energy += E;
    }//for each hit

    return energy;

}

double SinglePhotonQuickAna::CalcEShower(const std::vector<art::Ptr<recob::Hit>> &hits){    
    double energy[3] = {0., 0., 0.};

    //std::cout<<"SinglePhotonQuickAna::AnalyzeShowers() \t||\t Looking at shower with "<<hits.size() <<" hits on all planes"<<std::endl;

    //for each hit in the shower
    for (auto &thishitptr : hits){
        //check the plane
        int plane= thishitptr->View();

        //skip invalid planes     	
        if (plane > 2 || plane < 0)	continue;

        //calc the energy of the hit
        double E = QtoEConversionHit(thishitptr, plane);    

        //add the energy to the plane
        energy[plane] += E;
    }//for each hiti

    //find the max energy on a single plane
    double max = energy[0];
    for (double en: energy){
        if( en > max){
            max = en;
        }
    }
    // std::cout<<"SinglePhotonQuickAna::AnalyzeShowers() \t||\t The energy on each plane for this shower is "<<energy[0]<<", "<<energy[1]<<", "<<energy[2]<<std::endl;

    //return the highest energy on any of the planes
    return max;

}

double SinglePhotonQuickAna::GetQHit(art::Ptr<recob::Hit> thishitptr, int plane){
    double gain;
    //choose gain based on whether data/mc and by plane
    if (m_is_data == false &&  m_is_overlayed == false){
        gain = m_gain_mc[plane] ;
        //if (m_is_verbose) std::cout<<"the gain for mc on plane "<<plane<<" is "<<gain<<std::endl;
    } if (m_is_data == true ||  m_is_overlayed == true){
        gain = m_gain_data[plane] ;
        //if (m_is_verbose) std::cout<<"the gain for data on plane "<<plane<<" is "<<gain<<std::endl;

    }

    double Q = thishitptr->Integral()*gain;
    return Q;
}

double SinglePhotonQuickAna::QtoEConversionHit(art::Ptr<recob::Hit> thishitptr, int plane){
    return QtoEConversion(GetQHit(thishitptr, plane));

}

double SinglePhotonQuickAna::QtoEConversion(double Q){
    //return the energy value converted to MeV (the factor of 1e-6)
    double E = Q* m_work_function *1e-6 /m_recombination_factor;
    return E;

}


std::vector<double> SinglePhotonQuickAna::CalcdEdxFromdQdx(std::vector<double> dqdx){
    int n = dqdx.size();
    std::vector<double> dedx(n,0.0);
    for (int i = 0; i < n; i++){
        //std::cout<<"The dQ/dx is "<<dqdx[i]<<std::endl;
        dedx[i] = QtoEConversion(dqdx[i]);
        //std::cout<<"The dE/dx is "<<dedx[i]<<std::endl;
    }
    return dedx;
}


std::vector<double> SinglePhotonQuickAna::CalcdQdxShower(
        const art::Ptr<recob::Shower>& shower,
        const std::vector<art::Ptr<recob::Cluster>> & clusters, 
        std::map<art::Ptr<recob::Cluster>,    std::vector<art::Ptr<recob::Hit>> > &  clusterToHitMap ,int plane){
    //if(m_is_verbose) std::cout<<"SinglePhotonQuickAna::AnalyzeShowers() \t||\t The number of clusters in this shower is "<<clusters.size()<<std::endl;
    std::vector<double> dqdx;

    //get the 3D shower direction
    //note: in previous versions of the code there was a check to fix cases where the shower direction was inverted - this hasn't been implemented
    TVector3 shower_dir(shower->Direction().X(), shower->Direction().Y(),shower->Direction().Z());

    //calculate the pitch for this plane
    double pitch = getPitch(shower_dir, plane);	
    //if(m_is_verbose) std::cout<<"SinglePhotonQuickAna::AnalyzeShowers() \t||\t The pitch between the shower and plane "<<plane<<" is "<<pitch<<std::endl;

    //for all the clusters in the shower
    for (const art::Ptr<recob::Cluster> &thiscluster: clusters){
        //keep only clusters on the plane
        if(thiscluster->View() != plane) continue;

        //calculate the cluster direction
        std::vector<double> cluster_axis = {cos(thiscluster->StartAngle()), sin(thiscluster->StartAngle())};		

        //get the cluster start and and in CM
        //std::cout<<"for plane/tpc/cryo:"<<plane<<"/"<<m_TPC<<"/"<<m_Cryostat<<", fXTicksOffset: "<<theDetector->GetXTicksOffset(plane, m_TPC, m_Cryostat)<<" fXTicksCoefficient: "<<theDetector->GetXTicksCoefficient(m_TPC, m_Cryostat)<<std::endl;

        //convert the cluster start and end positions to time and wire coordinates
        std::vector<double> cluster_start = {thiscluster->StartWire() * m_wire_spacing,(thiscluster->StartTick() - theDetector->TriggerOffset())* m_time2cm};
        std::vector<double> cluster_end = {thiscluster->EndWire() * m_wire_spacing,(thiscluster->EndTick() - theDetector->TriggerOffset())* m_time2cm };

        //check that the cluster has non-zero length
        double length = sqrt(pow(cluster_end[0] - cluster_start[0], 2) + pow(cluster_end[1] - cluster_start[1], 2));
        //if(m_is_verbose) std::cout<<"SinglePhotonQuickAna::AnalyzeShowers() \t||\t The cluster length is "<<length<<std::endl;
        if (length <= 0){ 
            std::cout<<"skipping cluster on plane "<<plane<<", length = "<<length<<std::endl;
            continue;
        }


        //draw a rectangle around the cluster axis 
        std::vector<std::vector<double>> rectangle = buildRectangle(cluster_start, cluster_axis, m_width_dqdx_box, m_length_dqdx_box);	

        //get all the hits for this cluster
        std::vector<art::Ptr<recob::Hit>> hits =  clusterToHitMap[thiscluster];

        //for each hit in the cluster
        for (art::Ptr<recob::Hit> &thishit: hits){	
            //get the hit position in cm from the wire and time
            std::vector<double> thishit_pos = {thishit->WireID().Wire * m_wire_spacing, (thishit->PeakTime() - theDetector->TriggerOffset())* m_time2cm};

            //check if inside the box
            bool v2 = isInsidev2(thishit_pos, rectangle);
            if (v2){
                double q = GetQHit(thishit, plane); 
                double this_dqdx = q/pitch; 
                dqdx.push_back(this_dqdx);
            }//if hit falls inside the box

        }//for each hit inthe cluster
    }//for each cluster
    return dqdx;
}

double SinglePhotonQuickAna::getPitch(TVector3 shower_dir, int plane){
    //get the wire direction for this plane - values are hardcoded which isn't great but the TPC geom object gave weird values
    TVector3 wire_dir = getWireVec(plane);

    //take dot product of shower and wire dir
    double cos = getCoswrtWires(shower_dir, wire_dir);

    //want only positive values so take abs, normalize by the lengths of the shower and wire
    cos = abs(cos)/(wire_dir.Mag() * shower_dir.Mag());	

    //If the cos is 0 shower is perpendicular and therefore get infinite distance 
    if (cos == 0){ return std::numeric_limits<double>::max(); }

    //output is always >= the wire spacing
    return m_wire_spacing/cos;
}

TVector3 SinglePhotonQuickAna::getWireVec(int plane){
    TVector3 wire_dir;
    if (plane == 0){
        wire_dir = {0., -sqrt(3) / 2., 1 / 2.};
    } else if (plane == 1){
        wire_dir = {0., sqrt(3) / 2., 1 / 2.};
    } else if (plane == 2) {
        wire_dir = {0., 0., 1.};
    }
    return wire_dir;

}

double SinglePhotonQuickAna::getCoswrtWires(TVector3 shower_dir, TVector3 wire_dir){
    //take the dot product between the wire direction and the shower direction
    double cos = wire_dir.Dot(shower_dir);
    return cos;
}


double SinglePhotonQuickAna::degToRad(double deg){
    return deg * M_PI/180;
}

double SinglePhotonQuickAna::radToDeg(double rad){
    return rad * 180/M_PI;
}


double SinglePhotonQuickAna::getMeanHitWidthPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
    int nhits = 0;
    double widths = 0;
    for (art::Ptr<recob::Hit> thishitptr : hits){
        //check the plane
        int plane= thishitptr->View();

        //skip invalid planes       
        if (plane != this_plane) continue;

        widths += thishitptr->RMS(); // recob::Hit->RMS() returns RMS of the hit shape in tick units
        nhits++;


    }//for each hiti
    return   widths/(double)nhits;

}



int SinglePhotonQuickAna::getNHitsPlane(std::vector<art::Ptr<recob::Hit>> hits, int this_plane){
    int nhits = 0;
    for (art::Ptr<recob::Hit> thishitptr : hits){
        //check the plane
        int plane= thishitptr->View();

        //skip invalid planes       
        if (plane != this_plane) continue;

        nhits++;

    }//for each hiti
    return nhits;

}

std::vector<std::vector<double>> SinglePhotonQuickAna::buildRectangle(std::vector<double> cluster_start, std::vector<double> cluster_axis, double width, double length){
    std::vector<std::vector<double>> corners;

    //get the axis perpedicular to the cluster axis
    double perp_axis[2] = {-cluster_axis[1], cluster_axis[0]};

    //create a vector for each corner of the rectangle on the plane
    //c1 = bottom left corner
    std::vector<double> c1 = {cluster_start[0] + perp_axis[0] * width / 2,  cluster_start[1] + perp_axis[1] * width / 2};
    //c2 = top left corner
    std::vector<double> c2 = {c1[0] + cluster_axis[0] * length, c1[1] + cluster_axis[1] * length};
    //c3 = bottom right corner
    std::vector<double> c3 = {cluster_start[0] - perp_axis[0] * width / 2, cluster_start[1] - perp_axis[1] * width / 2};
    //c4 = top right corner
    std::vector<double> c4 ={c3[0] + cluster_axis[0] * length, c3[1] + cluster_axis[1] * length}; 

    //save each of the vectors
    corners.push_back(c1);
    corners.push_back(c2);
    corners.push_back(c4);
    corners.push_back(c3);
    return corners;
}

bool SinglePhotonQuickAna::insideBox(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle){
    //for a rectangle this is a known value but this is the most general
    int n_vertices = (int)rectangle.size();
    bool inside = false;
    int i, j = 0;
    //for each pair of vertices
    for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
        //if the hit y coordinate is between the y and x coordinates of two vertices
        if ( ((rectangle[i][1]> thishit_pos[1]) != (rectangle[j][1]>thishit_pos[1])) 
                &&(thishit_pos[0] < (rectangle[j][0]-rectangle[i][0]) * (thishit_pos[1]-rectangle[i][1]) / (rectangle[j][1]-rectangle[i][1]) + rectangle[i][0]) ){   
            if (inside == false){    
                inside = true;
            } else{
                inside = false;
            }
        }
    }
    return inside;
}

//determines if a point is inside the rectangle by summing the areas of the four triangles made by 
//if the point is inside, the sum of the triangles should exactly equal the area of the rectangle
//also returns true if the point is on the boundary
bool SinglePhotonQuickAna::isInsidev2(std::vector<double> thishit_pos, std::vector<std::vector<double >> rectangle){
    int n_vertices = (int)rectangle.size();
    //bool inside = false;
    int i, j = 0;
    double areas = 0;
    //for each pair of vertices
    for (i = 0, j = n_vertices-1; i < n_vertices; j = i++) {
        //calculate the area of a triangle with the point and two vertices
        double this_area = areaTriangle(rectangle[i][0], rectangle[i][1], rectangle[j][0], rectangle[j][1], thishit_pos[0], thishit_pos[1]);
        areas += this_area;
    }        
    //calc area of the rectangle
    double area_rectangle = m_width_dqdx_box* m_length_dqdx_box;

    //check the sum of areas match
    if (abs(areas - area_rectangle) <= 0.001 ){
        return true;
    }
    return false;
}

//area of a triangle given three vertices
double SinglePhotonQuickAna::areaTriangle(double x1, double y1, double x2, double y2, double x3, double y3){
    double num = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2);
    return abs(num)/2;
}

double SinglePhotonQuickAna::getMedian(std::vector<double> thisvector){
    //So return median if odd, average of median in even, if size==0, return the point. 

    //here the size corresponds to the max index

    int size = thisvector.size() - 1;
    //if no entries, return nonsense value
    if (size < 0) return NAN;
    if (size==0) return thisvector[size];

    //find index of median location
    double median;
    if (size%2 == 0){  // if vector has odd elements
        int ind = size/2;
        median = thisvector[ind];  
    } else{   // if vector has even number of elements
        int ind1 = size/2; 
        int ind2 = size/2-1;
        median = (thisvector[ind1]+thisvector[ind2])/2.0;
    }

    //return the value at median index
    return median;		
}





DEFINE_ART_MODULE(SinglePhotonQuickAna)

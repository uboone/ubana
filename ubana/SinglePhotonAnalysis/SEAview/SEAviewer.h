#ifndef SEAVIEWER_H
#define SEAVIEWER_H

#include <iostream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include <string>

#include <memory>

#include "larevt/SpaceChargeServices/SpaceChargeService.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"

#include "canvas/Utilities/ensurePointer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindOne.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TFile.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPrincipal.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TEllipse.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <map>
#include <sys/stat.h>


#include "seaDBSCAN.h"
namespace seaview {


    template <typename T>
        std::vector<size_t> seaview_sort_indexes(const std::vector<T> &v) {

            std::vector<size_t> idx(v.size());
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in v
            std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

            return idx;
        }



    struct cluster_score{
        int plane;
        int cluster_label;
        double point_score;
        int n_hits;

        double mean_wire;
        double max_wire;
        double min_wire;
        double mean_tick;
        double max_tick;
        double min_tick;

        double close_tick;
        double close_wire;
        double angle;//w.r.t shower primary

        double impact_parameter;

        double max_dist_tick;
        double mean_dist_tick;
        double min_dist_tick;
        double max_dist_wire;
        double mean_dist_wire;
        double min_dist_wire;

        double mean_dist;
        double max_dist;
        double min_dist;

        double pca_0;
        double pca_1;
        double pca_theta;

        int n_wires;
        int n_ticks;

        bool pass;

        cluster_score(int ip, int cl): plane(ip), cluster_label(cl){};
    };




    class cluster {

        public:

        double f_ImpactParameter;
        double f_FitSlope;
        double f_FitCons;
        double f_MeanADC;
        double f_AngleWRTShower;


        cluster(int ID, int plane, std::vector<std::vector<double>> &pts, std::vector<art::Ptr<recob::Hit>> &hits) :f_ID(ID), f_plane(plane), f_pts(pts), f_hits(hits), f_score(0,0), f_shower_remerge(-1) {

            f_npts = f_pts.size();
            if(pts.size() != hits.size()){
                std::cerr<<"seaviewer::cluster, input hits and pts not alligned"<<std::endl;
            }
            std::vector<double> wires(f_npts);
            std::vector<double> ticks(f_npts);
            for(int p =0; p< f_npts; ++p){
                wires[p]=f_pts[p][0];
                ticks[p]=f_pts[p][1];
            }
            TGraph af_graph(f_npts, &wires[0], &ticks[0]);
            f_graph = af_graph;

        };
        int setScore(cluster_score &in_score){ f_score = in_score;return 0;}

        cluster_score * getScore(){return &f_score;};
        int getID() {return f_ID;}
        int getN() {return f_npts;}
        int getPlane(){ return f_plane;}
        std::vector<std::vector<double>> getPTS(){return f_pts;}
        TGraph * getGraph(){ return &f_graph;}
        std::vector<art::Ptr<recob::Hit>>  getHits(){return f_hits;}
        int getShowerRemerge(){return f_shower_remerge;}
        int setShowerRemerge(int remerge_in){
                f_shower_remerge = remerge_in;
                return f_shower_remerge;
        }

        private:
        int f_ID;
        int f_npts;
        int f_plane;
        std::vector<std::vector<double>> f_pts;
        std::vector<art::Ptr<recob::Hit>> f_hits;
        cluster_score f_score;
        int f_shower_remerge;
        TGraph f_graph;
    };


    class SEAviewer {

        public:

            /// Default constructor
      SEAviewer(std::string tag,geo::GeometryCore const * geom,
		detinfo::DetectorPropertiesData theDetector );

            void configure(const fhicl::ParameterSet& pset){};

            int loadVertex(double m_vertex_pos_x, double m_vertex_pos_y, double m_vertex_pos_z);
            int addTrueVertex(double x, double y,double z);
            int addSliceHits(std::vector<art::Ptr<recob::Hit>>& hits);
            int addAllHits(std::vector<art::Ptr<recob::Hit>>& hits);
            int addPFParticleHits(std::vector<art::Ptr<recob::Hit>>& hits, std::string leg );
            int setBadChannelList(std::vector<std::pair<int,int>> &in);
            int addShower(art::Ptr<recob::Shower>&shr);
            int addTrack(art::Ptr<recob::Track>&trk);
            std::vector<int> calcUnassociatedHits();
            int setHitThreshold(double);
            int Print(double plot_distance);
            int runseaDBSCAN(double min_pts, double eps);

            double calcWire(double Y, double Z, int plane, int fTPC, int fCryostat, geo::GeometryCore const& geo ){
                return geo.WireCoordinate(geo::Point_t{0, Y, Z}, geo::PlaneID(fCryostat, fTPC, plane));
            }

            double calcTime(double X,int plane,int fTPC,int fCryostat, detinfo::DetectorPropertiesData const& detprop){
                return detprop.ConvertXToTicks(X, plane, fTPC,fCryostat);
            }


            std::vector<std::vector<double>> to2D(std::vector<double> & threeD);

            double dist_line_point( std::vector<double>&X1, std::vector<double>& X2, std::vector<double>& point);

            std::vector<double> analyzeClusters(double eps, std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap, std::vector<seaview::cluster> &vec_c );

            int SeaviewCompareToShowers(int p ,int cl, std::vector<art::Ptr<recob::Hit>>& hitz,double vertex_wire,double vertex_tick,   const std::vector<art::Ptr<recob::Shower>>& showers, std::map<art::Ptr<recob::Shower>,  art::Ptr<recob::PFParticle>> & showerToPFParticleMap,      const   std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::Hit>> > & pfParticleToHitsMap, double eps);


            cluster_score SeaviewScoreCluster(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hits, double vertex_wire, double vertex_tick, const art::Ptr<recob::Shower> &shower);

            TGraph* SeaviewGetNearestNpts(int p, int cl, std::vector<art::Ptr<recob::Hit>> &hitz, double vertex_wire, double vertex_tick, int Npts);



        protected:
            int n_pfps;
            int n_showers;
            int n_tracks;

            std::string tag;
            double hit_threshold;
            bool has_been_clustered;  
            std::vector<std::vector<TGraph>> vec_graphs;

            std::vector<std::string> vec_pfp_legend;
            // PFP, Plane: index
            std::vector<std::vector<std::vector<double>>> vec_ticks;
            std::vector<std::vector<std::vector<double>>> vec_chans;

            geo::GeometryCore const * geom;
            detinfo::DetectorPropertiesData theDetector ;

            double tick_shift;
            double chan_shift;

            double tick_max;
            double tick_min;
            std::vector<double> chan_max;
            std::vector<double> chan_min;

            std::vector<std::pair<int,int>> m_bad_channel_list;

            //Vertex
            std::vector<double> vertex_tick; 
            std::vector<double> vertex_chan; 
            std::vector<TGraph> vertex_graph;

            bool plot_true_vertex;
            std::vector<double> true_vertex_tick; 
            std::vector<double> true_vertex_chan; 
            std::vector<TGraph> true_vertex_graph;

            std::vector<art::Ptr<recob::Hit>> slice_hits;
            std::vector<art::Ptr<recob::Hit>> all_hits;
            std::map<art::Ptr<recob::Hit>,bool> map_unassociated_hits;
            std::map<art::Ptr<recob::Hit>, bool> map_slice_hits;

            std::vector<TGraph> vec_unass_graphs;
            std::vector<std::vector<double>> vec_unass_ticks;
            std::vector<std::vector<double>> vec_unass_chans;
            std::vector<std::vector<std::vector<double>>> vec_unass_pts;
            std::vector<std::vector<art::Ptr<recob::Hit>>> vec_unass_hits;

            std::vector<TGraph> vec_all_graphs;
            std::vector<std::vector<double>> vec_all_ticks;
            std::vector<std::vector<double>> vec_all_chans;

            std::vector<int> num_clusters;
            std::vector<std::vector<int>> cluster_labels;
            TRandom3 *rangen;

            std::vector<seaview::cluster> vec_clusters;  
            std::vector<art::Ptr<recob::Shower>> vec_showers;
            std::vector<art::Ptr<recob::Track>> vec_tracks;

    };

}// namespace

#endif

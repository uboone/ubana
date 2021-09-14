#ifndef _ClusterBuilder_h_
#define _ClusterBuilder_h_

#include <iostream>
#include <fstream>
#include <algorithm>

#include "TH2D.h"
#include "TCanvas.h"
#include <TROOT.h>
#include <TStyle.h>

namespace Connectedness {

struct Cluster {

int ID;

// Lists of bins in this cluster
std::vector<int> bins_x;
std::vector<int> bins_y;


};

class ClusterBuilder {

   public: 

      ClusterBuilder();
      ClusterBuilder(bool draw,std::string displaydir="");
      ~ClusterBuilder();

      // Setup functions

      void LoadDeadWireMaps(std::string dir="../");

      void SetThreshold(double threshold=1.8);
      void SetOffsets(int x_offset,int y_offset);
      void SetSeachArea(int x_max,int y_max);

      // Event processing functions

      void ReadData(std::vector<int> channel,std::vector<int> tick,std::vector<double> signal,std::string rse="");

      // Returns the ID of the cluster (check this against list of existing cluster IDs to
      // see if some merging has taken place). Returns -1,-1 if seed landed on empty bin 
      std::pair<int,int> MakeCluster(int seed_channel,int seed_tick,int ID);

      std::vector<Cluster> GetClusters();

      // Plane quality checks

      // Check seeds are not separated by dead wires
      bool SeedDeadWireCheck(std::vector<int> seeds_channel,std::vector<int> seeds_tick, int plane);


      // Empty the list of clusters, do before running different set of clusters for same event
      void ClearClusters();

      // Empty the list of clusters and delete histograms, do before reading a new event
      void Reset();
        
      void SetDisplayDir(std::string dir);

   private:

      bool DrawEverything=false;
      int XOffset = 0;
      int YOffset = 20;

      std::vector<int> DeadChannels_Plane0;
      std::vector<int> DeadChannels_Plane1;
      std::vector<int> DeadChannels_Plane2;

      double Threshold=1.8;

      TH2D *h_Raw = nullptr;
      TH2D *h_Binary = nullptr;
      TH2D *h_Clustered = nullptr;

      TCanvas *c;

      std::vector<Cluster> Clusters;

      std::pair<int,int> FindNearestOccupiedBin(TH2D *hist,int x, int y);
      int MaxSearchX = 2;
      int MaxSearchY = 15;

       void DeadWireFill(int plane);

       void Focus();

       std::string DisplayDir="";

   public:

      // pass = -1 (selection not applicable) , pass = 0 = not selected plane , pass = 1 = selected plane
      void DrawRaw(std::string rse="",int pass=-1);
      void DrawBinary(std::string rse="",int pass=-1);
      void DrawClustered(std::string rse="",int plane=-1,int pass=-1);
     
};

}

#endif



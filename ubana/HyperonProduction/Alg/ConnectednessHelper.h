#ifndef _ConnectednessHelper_h_
#define _ConnectednessHelper_h_


#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "cetlib_except/exception.h"

#include "lardataobj/RecoBase/Wire.h"

#include "ubana/HyperonProduction/Connectedness/Alg/ClusterBuilder.h"
#include "Position_To_Wire.h"

#include "TVector3.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Wiremap {

   int View;

   std::vector<int> Channel;
   std::vector<int> Tick;
   std::vector<double> Signal;

   void LoadActivity(std::vector<art::Ptr<recob::Wire>> wires){

      Channel.clear();
      Tick.clear();
      Signal.clear();

      // Iterate through all of the wires, record signal at every tick with nonzero signal
      for(const art::Ptr<recob::Wire> &wire : wires){

         if(wire->View() == View){

            // Get regions of interest        
            unsigned int NROI = wire->SignalROI().n_ranges();
            for(size_t i_roi=0; i_roi<NROI; ++i_roi){

               // Region of tick space with nonzero activity
               recob::Wire::RegionsOfInterest_t::datarange_t const& range = wire->SignalROI().range(i_roi);

               // Iterate through the ticks in this ROI, record signal
               unsigned int thisTick = range.begin_index();

               while(thisTick < range.end_index()){

                  Channel.push_back(wire->Channel());
                  Tick.push_back(thisTick);
                  Signal.push_back(wire->Signal().at(thisTick));

                  thisTick++;

               } // while(thisTick < range.end_index

            } // loop over ROI

         } // if(view == View)

      } // loop over wires

   }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Result performing the test on a single set of three tracks

struct CTSingleOutcome {

int Plane;

std::vector<int> SeedIndexes;
std::vector<int> OutputIndexes;
std::vector<int> OutputSizes;

std::vector<int> SeedChannels;
std::vector<int> SeedTicks;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// List of results from all the combinations of three tracks in event
// for all three planes

struct CTOutcome {

      std::vector<std::vector<int>> SeedIndexes_Plane0;
      std::vector<std::vector<int>> OutputIndexes_Plane0;
      std::vector<std::vector<int>> OutputSizes_Plane0;
      std::vector<std::vector<int>> SeedChannels_Plane0;
      std::vector<std::vector<int>> SeedTicks_Plane0;

      std::vector<std::vector<int>> SeedIndexes_Plane1;
      std::vector<std::vector<int>> OutputIndexes_Plane1;
      std::vector<std::vector<int>> OutputSizes_Plane1;
      std::vector<std::vector<int>> SeedChannels_Plane1;
      std::vector<std::vector<int>> SeedTicks_Plane1;

      std::vector<std::vector<int>> SeedIndexes_Plane2;
      std::vector<std::vector<int>> OutputIndexes_Plane2;
      std::vector<std::vector<int>> OutputSizes_Plane2;
      std::vector<std::vector<int>> SeedChannels_Plane2;
      std::vector<std::vector<int>> SeedTicks_Plane2;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ConnectednessHelper {

   public:

      ConnectednessHelper(bool draw);
      void LoadWireActivity(std::vector<art::Ptr<recob::Wire>> wires);
      void AddStartPositions(std::vector<TVector3> positions);
  
      std::vector<CTSingleOutcome> RunTest();

      std::vector<CTSingleOutcome> RunClustering(std::vector<int> indexes);        

      CTOutcome PrepareAndTestEvent(art::Event const& e,std::string wirelabel,std::vector<TVector3> trackstarts);

   private: 

      

      Connectedness::ClusterBuilder C_Plane0;
      Connectedness::ClusterBuilder C_Plane1;
      Connectedness::ClusterBuilder C_Plane2;

      bool Draw=false;

      Wiremap WM_Plane0;
      Wiremap WM_Plane1;
      Wiremap WM_Plane2;


      // Starting positions of tracks in 3D (do not correct for SC)
      std::vector<TVector3> Positions_3D;   
   
      // Starting positions of tracks in channel/tick space  
      std::vector<std::pair<int,int>> Positions_Plane0;
      std::vector<std::pair<int,int>> Positions_Plane1;
      std::vector<std::pair<int,int>> Positions_Plane2;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

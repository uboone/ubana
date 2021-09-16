#ifndef _ConnectednessHelper_cxx_
#define _ConnectednessHelper_cxx_

#include "ConnectednessHelper.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

ConnectednessHelper::ConnectednessHelper(bool draw) :
   C_Plane0(draw),
   C_Plane1(draw),
   C_Plane2(draw)
{
   Draw = draw;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConnectednessHelper::LoadWireActivity(std::vector<art::Ptr<recob::Wire>> wires){

   WM_Plane0.View = 0;
   WM_Plane0.LoadActivity(wires);
   C_Plane0.Reset();
   C_Plane0.ReadData(WM_Plane0.Channel,WM_Plane0.Tick,WM_Plane0.Signal);

   WM_Plane1.View = 1;
   WM_Plane1.LoadActivity(wires);
   C_Plane1.Reset();
   C_Plane1.ReadData(WM_Plane1.Channel,WM_Plane1.Tick,WM_Plane1.Signal);

   WM_Plane2.View = 2;
   WM_Plane2.LoadActivity(wires);
   C_Plane2.Reset();
   C_Plane2.ReadData(WM_Plane2.Channel,WM_Plane2.Tick,WM_Plane2.Signal);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ConnectednessHelper::AddStartPositions(std::vector<TVector3> positions){

   Positions_3D.clear();
   Positions_Plane0.clear();
   Positions_Plane1.clear();
   Positions_Plane2.clear();

   Positions_3D = positions;

   for(size_t i=0;i<Positions_3D.size();i++){

      Positions_Plane0.push_back(std::make_pair(U_wire(Positions_3D.at(i)),tick(Positions_3D.at(i))));
      Positions_Plane1.push_back(std::make_pair(V_wire(Positions_3D.at(i)),tick(Positions_3D.at(i))));
      Positions_Plane2.push_back(std::make_pair(Y_wire(Positions_3D.at(i)),tick(Positions_3D.at(i))));

   }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<CTSingleOutcome> ConnectednessHelper::RunClustering(std::vector<int> indexes){


   // Setup the output structs   

   CTSingleOutcome Outcome_Plane0;
   Outcome_Plane0.Plane = 0;
   Outcome_Plane0.SeedIndexes = indexes;  

   CTSingleOutcome Outcome_Plane1;
   Outcome_Plane1.Plane = 1;
   Outcome_Plane1.SeedIndexes = indexes;  

   CTSingleOutcome Outcome_Plane2;
   Outcome_Plane2.Plane = 2;
   Outcome_Plane2.SeedIndexes = indexes;  

   for(size_t i=0;i<indexes.size();i++){

      int index = indexes.at(i);

      Outcome_Plane0.SeedChannels.push_back(Positions_Plane0.at(index).first);
      Outcome_Plane0.SeedTicks.push_back(Positions_Plane0.at(index).second);

      Outcome_Plane1.SeedChannels.push_back(Positions_Plane1.at(index).first);
      Outcome_Plane1.SeedTicks.push_back(Positions_Plane1.at(index).second);

      Outcome_Plane2.SeedChannels.push_back(Positions_Plane2.at(index).first);
      Outcome_Plane2.SeedTicks.push_back(Positions_Plane2.at(index).second);

   }

   // Try generating clusters

   // Plane0 //

   for(size_t i=0;i<indexes.size();i++){
      int index = indexes.at(i);

      std::pair<int,int> id_and_size =  C_Plane0.MakeCluster(Positions_Plane0.at(index).first,Positions_Plane0.at(index).second,index);

      Outcome_Plane0.OutputIndexes.push_back(id_and_size.first);
      Outcome_Plane0.OutputSizes.push_back(id_and_size.second);

   }

   C_Plane0.ClearClusters();

   // Plane1 //

   for(size_t i=0;i<indexes.size();i++){
      int index = indexes.at(i);

      std::pair<int,int> id_and_size =  C_Plane1.MakeCluster(Positions_Plane1.at(index).first,Positions_Plane1.at(index).second,index);

      Outcome_Plane1.OutputIndexes.push_back(id_and_size.first);
      Outcome_Plane1.OutputSizes.push_back(id_and_size.second);

   }

   C_Plane1.ClearClusters();

   // Plane2 //

   for(size_t i=0;i<indexes.size();i++){
      int index = indexes.at(i);

      std::pair<int,int> id_and_size =  C_Plane2.MakeCluster(Positions_Plane2.at(index).first,Positions_Plane2.at(index).second,index);

      Outcome_Plane2.OutputIndexes.push_back(id_and_size.first);
      Outcome_Plane2.OutputSizes.push_back(id_and_size.second);

   }

   C_Plane2.ClearClusters();

   return {Outcome_Plane0,Outcome_Plane1,Outcome_Plane2};

}        

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<CTSingleOutcome> ConnectednessHelper::RunTest(){

   std::vector<CTSingleOutcome> AllOutcomes;

   // Need at least three tracks for the test
   if(Positions_3D.size() < 3) return AllOutcomes;

   // Try all 3 track combinations possible and record outcomes

   int nseeds = static_cast<int>(Positions_3D.size());

   for(int i=0;i<nseeds;i++){
      for(int j=i+1;j<nseeds;j++){
         for(int k=j+1;k<nseeds;k++){

            std::vector<CTSingleOutcome> Outcomes = RunClustering({i,j,k});

            AllOutcomes.push_back(Outcomes.at(0));
            AllOutcomes.push_back(Outcomes.at(1));
            AllOutcomes.push_back(Outcomes.at(2));

         }
      }
   }

   return AllOutcomes;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

CTOutcome ConnectednessHelper::PrepareAndTestEvent(art::Event const& e,std::string wirelabel,std::vector<TVector3> trackstarts){

   CTOutcome theOutcome;

   if(trackstarts.size() < 3) return theOutcome;

   art::Handle<std::vector<recob::Wire>> Handle_Wire;
   std::vector<art::Ptr<recob::Wire>> Vect_Wire;

   if(!e.getByLabel(wirelabel,Handle_Wire)) 
      throw cet::exception("ConnectednessHelper") << "Wire data product not found!" << std::endl;

   art::fill_ptr_vector(Vect_Wire,Handle_Wire);

   LoadWireActivity(Vect_Wire);
   AddStartPositions(trackstarts);

   // Vector containing the result for each combintation of three tracks
   std::vector<CTSingleOutcome> Outcomes = RunTest(); 


   // iterate over each combination of three tracks, store the result
   for(size_t i=0;i<Outcomes.size();i++){

      CTSingleOutcome this_Outcome = Outcomes.at(i);

      std::vector<int> this_SeedIndexes = this_Outcome.SeedIndexes;
      std::vector<int> this_OutputIndexes = this_Outcome.OutputIndexes;
      std::vector<int> this_OutputSizes = this_Outcome.OutputSizes;
      std::vector<int> this_SeedChannels = this_Outcome.SeedChannels;
      std::vector<int> this_SeedTicks = this_Outcome.SeedTicks;

      if(this_Outcome.Plane == 0){
         theOutcome.SeedIndexes_Plane0.push_back(this_SeedIndexes);
         theOutcome.OutputIndexes_Plane0.push_back(this_OutputIndexes);
         theOutcome.OutputSizes_Plane0.push_back(this_OutputSizes);
         theOutcome.SeedChannels_Plane0.push_back(this_SeedChannels);
         theOutcome.SeedTicks_Plane0.push_back(this_SeedTicks);
      } 
      else if(this_Outcome.Plane == 1){
         theOutcome.SeedIndexes_Plane1.push_back(this_SeedIndexes);
         theOutcome.OutputIndexes_Plane1.push_back(this_OutputIndexes);
         theOutcome.OutputSizes_Plane1.push_back(this_OutputSizes);
         theOutcome.SeedChannels_Plane1.push_back(this_SeedChannels);
         theOutcome.SeedTicks_Plane1.push_back(this_SeedTicks);
      } 
      else if(this_Outcome.Plane == 2){
         theOutcome.SeedIndexes_Plane2.push_back(this_SeedIndexes);
         theOutcome.OutputIndexes_Plane2.push_back(this_OutputIndexes);
         theOutcome.OutputSizes_Plane2.push_back(this_OutputSizes);
         theOutcome.SeedChannels_Plane2.push_back(this_SeedChannels);
         theOutcome.SeedTicks_Plane2.push_back(this_SeedTicks);
      } 
   }

   return theOutcome;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

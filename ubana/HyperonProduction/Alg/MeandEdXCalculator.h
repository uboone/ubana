#ifndef _MeandEdXCalculator_h_
#define _MeandEdXCalculator_h_

#include <vector>

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"

#include "TVector3.h"

//Define wire plane angles
const double plane0_wireangle = 30*6.28/360.0;
const double plane1_wireangle = -30*6.28/360.0;
const double plane2_wireangle = 90*6.28/360.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct dEdXStore {

   double Weight_Plane0 = 0;
   double Weight_Plane1 = 0;
   double Weight_Plane2 = 0;

   double Plane0 = -1;
   double Plane1 = -1;
   double Plane2 = -1;

   double ThreePlaneAverage = -1;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MeandEdXCalculator {

   public:

      MeandEdXCalculator();

      void SetTophatThresh(double thresh);

      double GetMeandEdX(art::Ptr<anab::Calorimetry> calo);
      double PlaneWeight(art::Ptr<recob::Track> track,int i_pl);
      dEdXStore ThreePlaneMeandEdX(art::Ptr<recob::Track> track,std::vector<art::Ptr<anab::Calorimetry>> calo_v);

   private:

      // Miniumum value of sin2(angle between track and wires)
      double TophatThresh = 0.175;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

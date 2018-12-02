/******************************************************************************
 * @file CalorimetryRadiusAlg.cxx
 * @brief  Inactive algorithm (it won't work if uncommented, it needs to be corrected since there's been few changes to the rest of the code). This algorithm looks at the ratio between hits from the two tracks and any other hits in a predefined circle of radius r. This is obviously more useful for background from neutrino interactions rather than cosmics.
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  CalorimetryRadiusAlg.h
 * ****************************************************************************/
#include "CalorimetryRadiusAlg.h"

namespace CalorimetryRadius
{
  // Constructor/destructor
  CalorimetryRadiusAlg::CalorimetryRadiusAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  CalorimetryRadiusAlg::~CalorimetryRadiusAlg()
  {}

  void CalorimetryRadiusAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fPfpLabel = pset.get<std::string>("PfpLabel");
    fHitLabel = pset.get<std::string>("HitLabel");
    fVerbose = pset.get<bool>("VerboseMode");
    fRadiusProfileLimits = pset.get<std::vector<double>>("RadiusProfileLimits");
    fRadiusProfileBins = pset.get<int>("RadiusProfileBins");
    fChannelNorm = pset.get<double>("ChannelNorm");
    fTickNorm = pset.get<double>("TickNorm");

    // Determine profile ticks
    double profileStep = (fRadiusProfileLimits[1] - fRadiusProfileLimits[0]) / float(fRadiusProfileBins);
    double currTick = fRadiusProfileLimits[0];
    for (int i=0; i<fRadiusProfileBins; i++)
    {
      currTick += profileStep;
      profileTicks.push_back(currTick);
    }
  }

  
  // Perform calorimetry analysis. At this stage we finally calculate all the charge deposited by the hits of track1 and track2 (or shower) within a radius from the assumed HSN decay vertex (cleanVertices[i]), for each decay vertex.
  // The second step looks at all the charge deposited by any hit in radius (which may come from hadronic interaction, in the case of background). And we finally calculate the ratio between the two (caloRatio). We would expect this ratio to be closer to 1 for signal, since HSN decaying in the detector don't interact with any particle, and we'd expect charge deposited by the two decay products to be the only charge within a certain radius from the decay point.
  // Now, we actually repeat this step for different radia in order to build up a profile. The width of the profile is given by fRadiusProfileLimits and the number of bins by fRadiusProfileBin.
  void CalorimetryRadiusAlg::PerformCalorimetry(
          art::Event const & evt,
          AuxEvent::EventTreeFiller & evd,
          std::vector<AuxVertex::DecayVertex>& decayVertices)
  {

    // Clean vectors that will be returned by function
    // evd.calo_prong1ChargeInRadius.clear();
    // evd.calo_prong2ChargeInRadius.clear();
    // evd.calo_totChargeInRadius.clear();
    // evd.calo_caloRatio.clear();

    // Prepare total hits vector (this is absolutely inefficient and potentially wrong, must be fixed)
    // art::InputTag hitTag {fHitLabel};
    // const auto& hitHandle = evt.getValidHandle< std::vector<recob::Hit> >(hitTag);
    // std::vector<art::Ptr<recob::Hit>> totHits;
    // for(std::vector<int>::size_type i=0; i!=(*hitHandle).size(); i++)
    // {
    //   art::Ptr<recob::Hit> hit(hitHandle,i);
    //   // totHits.push_back(hit);
    // }

    // // Loop through each decay vertex
    // for (std::vector<int>::size_type i=0; i!=decayVertices.size(); i++)
    // {
    //   AuxVertex::DecayVertex currentVertex = decayVertices[i];
    //   // Get useful quantities about the decay vertex currently being analyzed, like coordinates and parent indices
    //   int channel0[3] = {currentVertex.GetChannelLoc(0),currentVertex.GetChannelLoc(1),currentVertex.GetChannelLoc(2)};
    //   float tick0[3] = {currentVertex.GetTickLoc(0),currentVertex.GetTickLoc(1),currentVertex.GetTickLoc(2)};
      
    //   // Initialize the vector of hits that will be pushed back to the vector of vector of hits and returned by the function
    //   std::vector<recob::Hit const*> totHitsInMaxRadius;
    //   // Get the vector of hits for the prongs
    //   std::vector<recob::Hit const*> trackHits1 = currentVertex.GetProngHits(0);
    //   std::vector<recob::Hit const*> trackHits2 = currentVertex.GetProngHits(0);

    //   // Calculate calorimetry for track 1 within radius
    //   // parCharge1 is a vector, which contains all the charge due to particle1 in a circle of radius r around the vertex.
    //   // Each element of the vector is that integrated charge in increasing value of r
    //   // It starts out as a vector of size equal to the number of bins in radius profile, each element is equal to 0.
    //   // A loop goes then through each hit and for each radius size asks whether the hit is in it. If it is, the charge gets added to the total.
    //   // Now declare parCharge1 and fill it with zeros
    //   std::vector<float> prongCharge1;
    //   for (int j=0; j<fRadiusProfileBins; j++) prongCharge1.push_back(0.);
    //   for (recob::Hit const * hit : trackHits1)
    //   {
    //     int hitChannel = hit->Channel();
    //     double hitTick = (hit->EndTick() + hit->StartTick())/2.;
    //     int hitPlane = hit->View();
    //     for (int j=0; j<fRadiusProfileBins; j++)
    //     {
    //       double caloCut = profileTicks[j];
    //       bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(caloCut,2.));
    //       if (isInsideRadius)
    //       {
    //         double hitCharge = hit->Integral();
    //         prongCharge1[j] += hitCharge;
    //       } // END if hit is in current radius
    //     } // END loop for each radius
    //   } // END loop for each hit

    //   // Calculate calorimetry for track 2 within radius
    //   std::vector<float> prongCharge2;
    //   for (int j=0; j<fRadiusProfileBins; j++) prongCharge2.push_back(0.);
    //   for (auto hit : trackHits2)
    //   {
    //     int hitChannel = hit->Channel();
    //     double hitTick = (hit->EndTick() + hit->StartTick())/2.;
    //     int hitPlane = hit->View();
    //     for (int j=0; j<fRadiusProfileBins; j++)
    //     {
    //       double caloCut = profileTicks[j];
    //       bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(caloCut,2.));
    //       if (isInsideRadius)
    //       {
    //         double hitCharge = hit->Integral();
    //         prongCharge2[j] += hitCharge;
    //       } // END if hit is in current radius
    //     } // END loop for each radius
    //   } // END loop for each hit

    //   // Calculate total calorimetry within radius
    //   std::vector<float> totCharge;
    //   for (int j=0; j<fRadiusProfileBins; j++) totCharge.push_back(0.);
    //   for (auto hit : totHits)
    //   {
    //     int hitChannel = hit->Channel();
    //     double hitTick = (hit->EndTick() + hit->StartTick())/2.;
    //     int hitPlane = hit->View();
    //     for (int j=0; j<fRadiusProfileBins; j++)
    //     {
    //       double caloCut = profileTicks[j];
    //       bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(caloCut,2.));
    //       if (isInsideRadius)
    //       {
    //         double hitCharge = hit->Integral();
    //         totCharge[j] += hitCharge;

    //         // totHitsInMaxRadius are used to draw the evd, you need to do that only for the largest radius
    //         if (j == fRadiusProfileBins-1) totHitsInMaxRadius.push_back(hit);
    //       } // END if hit is in current radius
    //     } // END loop for each radius
    //   } // END loop for each hit

    //   // Assign the total hits in maximum radius vertex of pointers to the vertex
    //   currentVertex.SetTotHits(totHitsInMaxRadius);

    //   // Calculate the calorimetry ratio
    //   std::vector<float> thisCaloRatio;
    //   for (int j=0; j<fRadiusProfileBins; j++) thisCaloRatio.push_back((prongCharge1[j]+prongCharge2[j])/float(totCharge[j]));

    //   evd.calo_prong1ChargeInRadius.push_back(prongCharge1);
    //   evd.calo_prong2ChargeInRadius.push_back(prongCharge2);
    //   evd.calo_totChargeInRadius.push_back(totCharge);
    //   evd.calo_caloRatio.push_back(thisCaloRatio);
    // } // END loop for each decay vertex
    return;
  } // END function PerformCalorimetry
} // END namespace CalorimetryRadius
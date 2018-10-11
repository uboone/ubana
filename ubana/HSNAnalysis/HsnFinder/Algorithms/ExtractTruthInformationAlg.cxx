#include "ExtractTruthInformationAlg.h"

namespace ExtractTruthInformation
{
  // Constructor/destructor
  ExtractTruthInformationAlg::ExtractTruthInformationAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  ExtractTruthInformationAlg::~ExtractTruthInformationAlg()
  {}
  void ExtractTruthInformationAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fMinTpcBound = pset.get<std::vector<double>>("MinTpcBound");
    fMaxTpcBound = pset.get<std::vector<double>>("MaxTpcBound");
    fMcTrackLabel = pset.get<std::string>("McTrackLabel");
    fIsHSN = pset.get<bool>("IsHSN");
    fVerbose = pset.get<bool>("VerboseMode");
  }

  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void ExtractTruthInformationAlg::FillEventTreeWithTruth(
            art::Event const & evt,
            AuxEvent::EventTreeFiller & etf,
            std::vector<AuxVertex::DecayVertex> & ana_decayVertices)
  {
    // Prepare handle labels
    std::string mcTruthLabel = "generator";
    art::InputTag mcTruthTag {mcTruthLabel};
    const auto& mcTruthHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTruthTag);
    art::Ptr<simb::MCTruth> mct(mcTruthHandle,0);
    if (fIsHSN)
    {
      // Convention valid only for HSN! (first mcParticle in first mcTruth is either pi or mu from decay).
      const simb::MCParticle & mcPart = mct->GetParticle(0);
      etf.truth_vx = mcPart.Vx();
      etf.truth_vy = mcPart.Vy();
      etf.truth_vz = mcPart.Vz();
    }
    else
    {
      const simb::MCNeutrino & mcn = mct->GetNeutrino();
      const simb::MCParticle & nu = mcn.Nu();
      etf.truth_vx = nu.Vx();
      etf.truth_vy = nu.Vy();
      etf.truth_vz = nu.Vz();
    }

    // Loop through each decay vertex and find reco-truth distance
    float minDist = 1e10;
    int minDistInd = -1;
    for (std::vector<int>::size_type i=0; i!=ana_decayVertices.size(); i++)
    {
      AuxVertex::DecayVertex dv = ana_decayVertices[i];
      float distance = sqrt( pow((etf.truth_vx - dv.fX),2.) + pow((etf.truth_vy - dv.fY),2.) + pow((etf.truth_vz - dv.fZ),2.) );
      if ( distance < minDist )
      {
        minDist = distance;
        minDistInd = i;
      }
      etf.recoTruthDistances.push_back(distance);
      etf.isClosestToTruth.push_back(0);
    }
    etf.isClosestToTruth[minDistInd] = 1;
  } // END function DetermineExtractTruthInformation


  void ExtractTruthInformationAlg::FillDrawTreeWithTruth(
            art::Event const & evt,
            AuxEvent::DrawTreeFiller & dtf)
  {

    // Clear vectors
    dtf.truth_primaryTracks_start_p0_wireCoordinates.clear();
    dtf.truth_primaryTracks_start_p0_tickCoordinates.clear();
    dtf.truth_primaryTracks_start_p1_wireCoordinates.clear();
    dtf.truth_primaryTracks_start_p1_tickCoordinates.clear();
    dtf.truth_primaryTracks_start_p2_wireCoordinates.clear();
    dtf.truth_primaryTracks_start_p2_tickCoordinates.clear();
    dtf.truth_primaryTracks_end_p0_wireCoordinates.clear();
    dtf.truth_primaryTracks_end_p0_tickCoordinates.clear();
    dtf.truth_primaryTracks_end_p1_wireCoordinates.clear();
    dtf.truth_primaryTracks_end_p1_tickCoordinates.clear();
    dtf.truth_primaryTracks_end_p2_wireCoordinates.clear();
    dtf.truth_primaryTracks_end_p2_tickCoordinates.clear();
    dtf.truth_secondaryTracks_start_p0_wireCoordinates.clear();
    dtf.truth_secondaryTracks_start_p0_tickCoordinates.clear();
    dtf.truth_secondaryTracks_start_p1_wireCoordinates.clear();
    dtf.truth_secondaryTracks_start_p1_tickCoordinates.clear();
    dtf.truth_secondaryTracks_start_p2_wireCoordinates.clear();
    dtf.truth_secondaryTracks_start_p2_tickCoordinates.clear();
    dtf.truth_secondaryTracks_end_p0_wireCoordinates.clear();
    dtf.truth_secondaryTracks_end_p0_tickCoordinates.clear();
    dtf.truth_secondaryTracks_end_p1_wireCoordinates.clear();
    dtf.truth_secondaryTracks_end_p1_tickCoordinates.clear();
    dtf.truth_secondaryTracks_end_p2_wireCoordinates.clear();
    dtf.truth_secondaryTracks_end_p2_tickCoordinates.clear();
    dtf.truth_primaryShowers_start_p0_wireCoordinates.clear();
    dtf.truth_primaryShowers_start_p0_tickCoordinates.clear();
    dtf.truth_primaryShowers_start_p1_wireCoordinates.clear();
    dtf.truth_primaryShowers_start_p1_tickCoordinates.clear();
    dtf.truth_primaryShowers_start_p2_wireCoordinates.clear();
    dtf.truth_primaryShowers_start_p2_tickCoordinates.clear();
    dtf.truth_primaryShowers_end_p0_wireCoordinates.clear();
    dtf.truth_primaryShowers_end_p0_tickCoordinates.clear();
    dtf.truth_primaryShowers_end_p1_wireCoordinates.clear();
    dtf.truth_primaryShowers_end_p1_tickCoordinates.clear();
    dtf.truth_primaryShowers_end_p2_wireCoordinates.clear();
    dtf.truth_primaryShowers_end_p2_tickCoordinates.clear();
    dtf.truth_secondaryShowers_start_p0_wireCoordinates.clear();
    dtf.truth_secondaryShowers_start_p0_tickCoordinates.clear();
    dtf.truth_secondaryShowers_start_p1_wireCoordinates.clear();
    dtf.truth_secondaryShowers_start_p1_tickCoordinates.clear();
    dtf.truth_secondaryShowers_start_p2_wireCoordinates.clear();
    dtf.truth_secondaryShowers_start_p2_tickCoordinates.clear();
    dtf.truth_secondaryShowers_end_p0_wireCoordinates.clear();
    dtf.truth_secondaryShowers_end_p0_tickCoordinates.clear();
    dtf.truth_secondaryShowers_end_p1_wireCoordinates.clear();
    dtf.truth_secondaryShowers_end_p1_tickCoordinates.clear();
    dtf.truth_secondaryShowers_end_p2_wireCoordinates.clear();
    dtf.truth_secondaryShowers_end_p2_tickCoordinates.clear();

    // Prepare handle labels
    std::string mcTruthLabel = "generator";
    art::InputTag mcTruthTag {mcTruthLabel};
    const auto& mcTruthHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTruthTag);
    art::Ptr<simb::MCTruth> mct(mcTruthHandle,0);
    // And for tracks and showers
    art::InputTag mcTrackTag {fMcTrackLabel};
    const auto& mcTrackHandle = evt.getValidHandle< std::vector<sim::MCTrack> >(mcTrackTag);
    art::InputTag mcShowerTag {fMcTrackLabel};
    const auto& mcShowerHandle = evt.getValidHandle< std::vector<sim::MCShower> >(mcShowerTag);

    // Extract vertex location
    float nu_xyz[3];
    if (fIsHSN)
    {
      // Convention valid only for HSN! (first mcParticle in first mcTruth is either pi or mu from decay).
      const simb::MCParticle & mcPart = mct->GetParticle(0);
      nu_xyz[0] = mcPart.Vx();
      nu_xyz[1] = mcPart.Vy();
      nu_xyz[2] = mcPart.Vz();
    }
    else
    {
      const simb::MCNeutrino & mcn = mct->GetNeutrino();
      const simb::MCParticle & nu = mcn.Nu();
      nu_xyz[0] = nu.Vx();
      nu_xyz[1] = nu.Vy();
      nu_xyz[2] = nu.Vz();
    }

    // Placeholders
    std::vector<int> channelLoc;
    std::vector<float> tickLoc;
    // Now convert to wire-channel if inside TPC
    XYZtoWireTick(nu_xyz, channelLoc, tickLoc);
    dtf.truth_dv_p0_wireCoordinates = channelLoc[0];
    dtf.truth_dv_p0_tickCoordinates = tickLoc[0];
    dtf.truth_dv_p1_wireCoordinates = channelLoc[1];
    dtf.truth_dv_p1_tickCoordinates = tickLoc[1];
    dtf.truth_dv_p2_wireCoordinates = channelLoc[2];
    dtf.truth_dv_p2_tickCoordinates = tickLoc[2];

    // Finally, new vertex position can enlarge drawing window
    if (dtf.truth_dv_p0_wireCoordinates<dtf.p0_minWire) dtf.p0_minWire = dtf.truth_dv_p0_wireCoordinates;
    if (dtf.truth_dv_p0_wireCoordinates>dtf.p0_maxWire) dtf.p0_maxWire = dtf.truth_dv_p0_wireCoordinates;
    if (dtf.truth_dv_p1_wireCoordinates<dtf.p1_minWire) dtf.p1_minWire = dtf.truth_dv_p1_wireCoordinates;
    if (dtf.truth_dv_p1_wireCoordinates>dtf.p1_maxWire) dtf.p1_maxWire = dtf.truth_dv_p1_wireCoordinates;
    if (dtf.truth_dv_p2_wireCoordinates<dtf.p2_minWire) dtf.p2_minWire = dtf.truth_dv_p2_wireCoordinates;
    if (dtf.truth_dv_p2_wireCoordinates>dtf.p2_maxWire) dtf.p2_maxWire = dtf.truth_dv_p2_wireCoordinates;
    if (dtf.truth_dv_p0_tickCoordinates<dtf.p0_minTick) dtf.p0_minTick = dtf.truth_dv_p0_tickCoordinates;
    if (dtf.truth_dv_p0_tickCoordinates>dtf.p0_maxTick) dtf.p0_maxTick = dtf.truth_dv_p0_tickCoordinates;
    if (dtf.truth_dv_p1_tickCoordinates<dtf.p1_minTick) dtf.p1_minTick = dtf.truth_dv_p1_tickCoordinates;
    if (dtf.truth_dv_p1_tickCoordinates>dtf.p1_maxTick) dtf.p1_maxTick = dtf.truth_dv_p1_tickCoordinates;
    if (dtf.truth_dv_p2_tickCoordinates<dtf.p2_minTick) dtf.p2_minTick = dtf.truth_dv_p2_tickCoordinates;
    if (dtf.truth_dv_p2_tickCoordinates>dtf.p2_maxTick) dtf.p2_maxTick = dtf.truth_dv_p2_tickCoordinates;

    // Now do it for tracks
    int nPrimaryTracks = 0;
    int nSecondaryTracks = 0;
    for(std::vector<int>::size_type k=0; k!=(*mcTrackHandle).size(); k++)
    {
      art::Ptr<sim::MCTrack> mctrack(mcTrackHandle,k);
      bool mctrackIsPrimary = (mctrack->Process()=="primary");
      if (mctrackIsPrimary)
      {
        if ((*mctrack).size()>1)
        {
          float start[3] = {(float) mctrack->at(0).X(),(float) mctrack->at(0).Y(),(float) mctrack->at(0).Z()};
          if (IsInsideTpc(start))
          {
            nPrimaryTracks += 1;

            float end[3] = {0.,0.,0.};
            for(std::vector<int>::size_type j=0; j!=(*mctrack).size(); j++)
            {
              float xyz[3] = {(float) mctrack->at(j).X(),(float) mctrack->at(j).Y(),(float) mctrack->at(j).Z()};
              if (IsInsideTpc(xyz))
              {
                end[0] = xyz[0];
                end[1] = xyz[1];
                end[2] = xyz[2];
              }
            }

            // Placeholders
            std::vector<int> channelLoc;
            std::vector<float> tickLoc;
            // Now convert start to wire-channel
            XYZtoWireTick(start, channelLoc, tickLoc);
            dtf.truth_primaryTracks_start_p0_wireCoordinates.push_back(channelLoc[0]);
            dtf.truth_primaryTracks_start_p0_tickCoordinates.push_back(tickLoc[0]);
            dtf.truth_primaryTracks_start_p1_wireCoordinates.push_back(channelLoc[1]);
            dtf.truth_primaryTracks_start_p1_tickCoordinates.push_back(tickLoc[1]);
            dtf.truth_primaryTracks_start_p2_wireCoordinates.push_back(channelLoc[2]);
            dtf.truth_primaryTracks_start_p2_tickCoordinates.push_back(tickLoc[2]);
            // Now convert end to wire-channel 
            XYZtoWireTick(end, channelLoc, tickLoc);
            dtf.truth_primaryTracks_end_p0_wireCoordinates.push_back(channelLoc[0]);
            dtf.truth_primaryTracks_end_p0_tickCoordinates.push_back(tickLoc[0]);
            dtf.truth_primaryTracks_end_p1_wireCoordinates.push_back(channelLoc[1]);
            dtf.truth_primaryTracks_end_p1_tickCoordinates.push_back(tickLoc[1]);
            dtf.truth_primaryTracks_end_p2_wireCoordinates.push_back(channelLoc[2]);
            dtf.truth_primaryTracks_end_p2_tickCoordinates.push_back(tickLoc[2]);
          } // END IsInsideTpc condition
        } // END condition on mcpart length
      } // END particle is primary condition
      else
      {
        if ((*mctrack).size()>1)
        {
          float start[3] = {(float) mctrack->at(0).X(),(float) mctrack->at(0).Y(),(float) mctrack->at(0).Z()};
          if (IsInsideTpc(start))
          {
            nSecondaryTracks += 1;

            float end[3] = {(float) mctrack->at(0).X(),(float) mctrack->at(0).Y(),(float) mctrack->at(0).Z()};
            for(std::vector<int>::size_type j=0; j!=(*mctrack).size(); j++)
            {
              float xyz[3] = {(float) mctrack->at(j).X(),(float) mctrack->at(j).Y(),(float) mctrack->at(j).Z()};
              if (IsInsideTpc(xyz))
              {
                end[0] = xyz[0];
                end[1] = xyz[1];
                end[2] = xyz[2];
              }
            }

            // Placeholders
            std::vector<int> channelLoc;
            std::vector<float> tickLoc;
            // Now convert start to wire-channel
            XYZtoWireTick(start, channelLoc, tickLoc);
            dtf.truth_secondaryTracks_start_p0_wireCoordinates.push_back(channelLoc[0]);
            dtf.truth_secondaryTracks_start_p0_tickCoordinates.push_back(tickLoc[0]);
            dtf.truth_secondaryTracks_start_p1_wireCoordinates.push_back(channelLoc[1]);
            dtf.truth_secondaryTracks_start_p1_tickCoordinates.push_back(tickLoc[1]);
            dtf.truth_secondaryTracks_start_p2_wireCoordinates.push_back(channelLoc[2]);
            dtf.truth_secondaryTracks_start_p2_tickCoordinates.push_back(tickLoc[2]);
            // Now convert end to wire-channel 
            XYZtoWireTick(end, channelLoc, tickLoc);
            dtf.truth_secondaryTracks_end_p0_wireCoordinates.push_back(channelLoc[0]);
            dtf.truth_secondaryTracks_end_p0_tickCoordinates.push_back(tickLoc[0]);
            dtf.truth_secondaryTracks_end_p1_wireCoordinates.push_back(channelLoc[1]);
            dtf.truth_secondaryTracks_end_p1_tickCoordinates.push_back(tickLoc[1]);
            dtf.truth_secondaryTracks_end_p2_wireCoordinates.push_back(channelLoc[2]);
            dtf.truth_secondaryTracks_end_p2_tickCoordinates.push_back(tickLoc[2]);
          } // END IsInsideTpc condition
        } // END condition on mcpart length
      } // END particle is NOT primary condition
    } // END loop through mctracks
    dtf.truth_nPrimaryTracks = nPrimaryTracks;
    dtf.truth_nSecondaryTracks = nSecondaryTracks;

    // Now do it for showers
    int nPrimaryShowers = 0;
    int nSecondaryShowers = 0;
    for(std::vector<int>::size_type k=0; k!=(*mcShowerHandle).size(); k++)
    {
      art::Ptr<sim::MCShower> mcshower(mcShowerHandle,k);
      bool mcshowerIsPrimary = (mcshower->Process()=="primary");
      if (mcshowerIsPrimary)
      {
        float start[3] = {(float) mcshower->Start().X(),(float) mcshower->Start().Y(),(float) mcshower->Start().Z()};
        float end[3] = {(float) mcshower->End().X(),(float) mcshower->End().Y(),(float) mcshower->End().Z()};
        if (IsInsideTpc(start) && IsInsideTpc(end))
        {
          nPrimaryShowers += 1;
          // Placeholders
          std::vector<int> channelLoc;
          std::vector<float> tickLoc;
          // Now convert start to wire-channel
          XYZtoWireTick(start, channelLoc, tickLoc);
          dtf.truth_primaryShowers_start_p0_wireCoordinates.push_back(channelLoc[0]);
          dtf.truth_primaryShowers_start_p0_tickCoordinates.push_back(tickLoc[0]);
          dtf.truth_primaryShowers_start_p1_wireCoordinates.push_back(channelLoc[1]);
          dtf.truth_primaryShowers_start_p1_tickCoordinates.push_back(tickLoc[1]);
          dtf.truth_primaryShowers_start_p2_wireCoordinates.push_back(channelLoc[2]);
          dtf.truth_primaryShowers_start_p2_tickCoordinates.push_back(tickLoc[2]);
          // Now convert end to wire-channel 
          XYZtoWireTick(end, channelLoc, tickLoc);
          dtf.truth_primaryShowers_end_p0_wireCoordinates.push_back(channelLoc[0]);
          dtf.truth_primaryShowers_end_p0_tickCoordinates.push_back(tickLoc[0]);
          dtf.truth_primaryShowers_end_p1_wireCoordinates.push_back(channelLoc[1]);
          dtf.truth_primaryShowers_end_p1_tickCoordinates.push_back(tickLoc[1]);
          dtf.truth_primaryShowers_end_p2_wireCoordinates.push_back(channelLoc[2]);
          dtf.truth_primaryShowers_end_p2_tickCoordinates.push_back(tickLoc[2]);
        } // END IsInsideTpc condition
      } // END particle is primary condition
      else
      {
        float start[3] = {(float) mcshower->Start().X(),(float) mcshower->Start().Y(),(float) mcshower->Start().Z()};
        float end[3] = {(float) mcshower->End().X(),(float) mcshower->End().Y(),(float) mcshower->End().Z()};
        if (IsInsideTpc(start) && IsInsideTpc(end))
        {
          nSecondaryShowers += 1;
          // Placeholders
          std::vector<int> channelLoc;
          std::vector<float> tickLoc;
          // Now convert start to wire-channel
          XYZtoWireTick(start, channelLoc, tickLoc);
          dtf.truth_secondaryShowers_start_p0_wireCoordinates.push_back(channelLoc[0]);
          dtf.truth_secondaryShowers_start_p0_tickCoordinates.push_back(tickLoc[0]);
          dtf.truth_secondaryShowers_start_p1_wireCoordinates.push_back(channelLoc[1]);
          dtf.truth_secondaryShowers_start_p1_tickCoordinates.push_back(tickLoc[1]);
          dtf.truth_secondaryShowers_start_p2_wireCoordinates.push_back(channelLoc[2]);
          dtf.truth_secondaryShowers_start_p2_tickCoordinates.push_back(tickLoc[2]);
          // Now convert end to wire-channel 
          XYZtoWireTick(end, channelLoc, tickLoc);
          dtf.truth_secondaryShowers_end_p0_wireCoordinates.push_back(channelLoc[0]);
          dtf.truth_secondaryShowers_end_p0_tickCoordinates.push_back(tickLoc[0]);
          dtf.truth_secondaryShowers_end_p1_wireCoordinates.push_back(channelLoc[1]);
          dtf.truth_secondaryShowers_end_p1_tickCoordinates.push_back(tickLoc[1]);
          dtf.truth_secondaryShowers_end_p2_wireCoordinates.push_back(channelLoc[2]);
          dtf.truth_secondaryShowers_end_p2_tickCoordinates.push_back(tickLoc[2]);
        } // END IsInsideTpc condition
      } // END particle is NOT primary condition
    } // END loop through mcShowers
    dtf.truth_nPrimaryShowers = nPrimaryShowers;
    dtf.truth_nSecondaryShowers = nSecondaryShowers;

  } // END function FillDrawTree

  // Determine if coordinates are inside TPC
  bool ExtractTruthInformationAlg::IsInsideTpc(const float* xyz)
  {
    double extraEdge = 0;
    bool isInsideX = (xyz[0]>fMinTpcBound[0]+extraEdge &&
      xyz[0]<fMaxTpcBound[0]-extraEdge);
    bool isInsideY = (xyz[1]>fMinTpcBound[1]+extraEdge &&
      xyz[1]<fMaxTpcBound[1]-extraEdge);
    bool isInsideZ = (xyz[2]>fMinTpcBound[2]+extraEdge &&
      xyz[2]<fMaxTpcBound[2]-extraEdge);
    bool isInsideTpc = (isInsideX && isInsideY && isInsideZ);
    return isInsideTpc;
  } // END function IsInsideTpc

  // Convert XYZ coordinates to wire-tick coordinates
  void ExtractTruthInformationAlg::XYZtoWireTick(const float* xyz, std::vector<int>& channelLoc, std::vector<float>& tickLoc)
  {
    if (IsInsideTpc(xyz))
    {
      raw::ChannelID_t channel0 = fGeometry->NearestChannel(xyz,0);
      raw::ChannelID_t channel1 = fGeometry->NearestChannel(xyz,1);
      raw::ChannelID_t channel2 = fGeometry->NearestChannel(xyz,2);
      double tick0 = fDetectorProperties->ConvertXToTicks(xyz[0], 0, 0, 0);
      double tick1 = fDetectorProperties->ConvertXToTicks(xyz[0], 1, 0, 0);
      double tick2 = fDetectorProperties->ConvertXToTicks(xyz[0], 2, 0, 0);
      channelLoc = {(int) channel0,(int) channel1,(int) channel2};
      tickLoc = { (float) tick0, (float) tick1, (float) tick2};
    }
    else
    {
      channelLoc = {-999,-999,-999};
      tickLoc = {-999.0,-999.0,-999.0};
    }
  } // END function XYZtoWireTick
} // END namespace ExtractTruthInformation

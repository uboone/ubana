/******************************************************************************
 * @file DecayVertex.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.h
 * ****************************************************************************/

// Decay vertex header
#include "DecayVertex.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace AuxVertex
{
  DecayVertex::DecayVertex()
  {}
  DecayVertex::~DecayVertex()
  {}

  DecayVertex::DecayVertex(
            const art::Ptr<recob::Vertex> &nuVertex,
            const art::Ptr<recob::Vertex> &t1Vertex,
            const art::Ptr<recob::Vertex> &t2Vertex,
            const art::Ptr<recob::Track> &t1Track,
            const art::Ptr<recob::Track> &t2Track,
            const std::vector<art::Ptr<recob::Hit>> &t1Hits,
            const std::vector<art::Ptr<recob::Hit>> &t2Hits,
            const art::Ptr<recob::MCSFitResult> &t1Mcs,
            const art::Ptr<recob::MCSFitResult> &t2Mcs)
  {
    /*
    This function creates a HSN decay vertex candidate object.
    It is assigned when a neutrino with two (and only two) track daughters has been found.
    This object is build using the vertex pointers to the neutrino and the two start points of the tracks, the actual track objects and all the hits associated with the two track objects.
    */

    // Set default mock attributes for the Decay Vertex candidate
    // The real values for these attributes will be assigned later
    fChannelLoc = {-1,-1,-1};
    fTickLoc = {-1.,-1.,-1.};
    fProngChannelLoc = {{-1,-1,-1}, {-1,-1,-1}};
    fProngTickLoc = {{-1.,-1.,-1.}, {-1.,-1.,-1.}};
    fIsDetLocAssigned = false;
    fIsInsideTPC = false;

    // Store pointers to reconstructed objects provided as input in internal attributes
    fNuVertex = nuVertex;
    fProngVertex = {t1Vertex, t2Vertex};
    fProngTrack = {t1Track, t2Track};
    fProngHits = {t1Hits, t2Hits};
    fProngMcs = {t1Mcs, t2Mcs};

    // Use pointers to reconstructed objects to obtain start coordinates for vertices.
    double nuVertexPosition[3], t1VertexPosition[3], t2VertexPosition[3];
    nuVertex->XYZ(nuVertexPosition);
    t1Vertex->XYZ(t1VertexPosition);
    t2Vertex->XYZ(t2VertexPosition);
    fX = nuVertexPosition[0];
    fY = nuVertexPosition[1];
    fZ = nuVertexPosition[2];
    fProngX = {(float) t1VertexPosition[0], (float) t2VertexPosition[0]};
    fProngY = {(float) t1VertexPosition[1], (float) t2VertexPosition[1]};
    fProngZ = {(float) t1VertexPosition[2], (float) t2VertexPosition[2]};
    // Do the same for tracks.
    fProngStartX = {(float) fProngTrack[0]->Start().X(), (float) fProngTrack[1]->Start().X()};
    fProngStartY = {(float) fProngTrack[0]->Start().Y(), (float) fProngTrack[1]->Start().Y()};
    fProngStartZ = {(float) fProngTrack[0]->Start().Z(), (float) fProngTrack[1]->Start().Z()};
    fProngEndX = {(float) fProngTrack[0]->End()[0], (float) fProngTrack[1]->End()[0]};
    fProngEndY = {(float) fProngTrack[0]->End()[1], (float) fProngTrack[1]->End()[1]};
    fProngEndZ = {(float) fProngTrack[0]->End()[2], (float) fProngTrack[1]->End()[2]};

    // Calculate length, theta, phi and number of hits
    fProngLength = {(float) fProngTrack[0]->Length(), (float) fProngTrack[1]->Length()};
    fProngTheta = {(float) fProngTrack[0]->Theta(), (float) fProngTrack[1]->Theta()};
    fProngPhi = {(float) fProngTrack[0]->Phi(), (float) fProngTrack[1]->Phi()};
    fProngNumHits = {(int) t1Hits.size(), (int) t2Hits.size()};

    // Calculate track directions
    fProngDirX = {(float) fProngTrack[0]->VertexMomentumVector().X(), (float) fProngTrack[1]->VertexMomentumVector().X()};
    fProngDirY = {(float) fProngTrack[0]->VertexMomentumVector().Y(), (float) fProngTrack[1]->VertexMomentumVector().Y()};
    fProngDirZ = {(float) fProngTrack[0]->VertexMomentumVector().Z(), (float) fProngTrack[1]->VertexMomentumVector().Z()};

    // Calculate momenta and quantities related to them, using range and MCS method
    SetHypothesisLabels();
    SetMomentumQuantities_ByRange();
    SetMomentumQuantities_ByMCS();
  } //  END constructor DecayVertex

  // Getters
  // Pointers
  art::Ptr<recob::Vertex> DecayVertex::GetNuVertex() const {return fNuVertex;}
  art::Ptr<recob::Vertex> DecayVertex::GetProngVertex(int prong) const {return fProngVertex[prong];}
  art::Ptr<recob::Track> DecayVertex::GetProngTrack(int prong) const {return fProngTrack[prong];}
  std::vector<art::Ptr<recob::Hit>> DecayVertex::GetProngHits(int prong) const {return fProngHits[prong];}
  std::vector<art::Ptr<recob::Hit>> DecayVertex::GetTotHits() const {return fTotHitsInMaxRadius;}

  // Setters
  void DecayVertex::SetChannelLoc(int channel0, int channel1, int channel2) {fChannelLoc = {channel0,channel1,channel2}; return;}
  void DecayVertex::SetTickLoc(float tick0, float tick1, float tick2) {fTickLoc = {tick0, tick1, tick2}; return;}
  void DecayVertex::SetProngChannelLoc(int prong, int channel0, int channel1, int channel2) {fProngChannelLoc[prong] =  {channel0,channel1,channel2}; return;}
  void DecayVertex::SetProngTickLoc(int prong, float tick0, float tick1, float tick2) {fProngTickLoc[prong] = {tick0, tick1, tick2}; return;}
  void DecayVertex::SetProngXYZ(int prong, float x, float y, float z) {fProngX[prong] = x; fProngY[prong] = y; fProngZ[prong] = z; return;}
  void DecayVertex::SetIsInsideTPC(bool val) {fIsInsideTPC = val; return;}
  void DecayVertex::SetIsDetLocAssigned(bool val) {fIsDetLocAssigned = val; return;}
  void DecayVertex::SetTotHits(std::vector<art::Ptr<recob::Hit>> totHitsInMaxRadius) {fTotHitsInMaxRadius = totHitsInMaxRadius; return;}

  void DecayVertex::SetDetectorCoordinates(
    const std::vector<double>& minTpcBound,
    const std::vector<double>& maxTpcBound,
    geo::GeometryCore const* geometry,
    detinfo::DetectorProperties const* detectorProperties)
  {
    /* This functions call special geometry and detectorProperties objects
    in order to translate x,y,z coordinates in wire,tick coordinates.
    This allows the determination of the vertices location in the event display
    */

    // Get spatial coordinates and mark vertex as assigned
    float xyz[3] = {fX,fY,fZ};
    float prong1_xyz[3] = {fProngX[0],fProngY[0],fProngZ[0]};
    float prong2_xyz[3] = {fProngX[1],fProngY[1],fProngZ[1]};

    fIsDetLocAssigned = true;

    // Check whether coordinates are inside TPC
    double extraEdge = 0;

    bool nuIsInsideX = (xyz[0]>minTpcBound[0]+extraEdge &&
      xyz[0]<maxTpcBound[0]-extraEdge);
    bool nuIsInsideY = (xyz[1]>minTpcBound[1]+extraEdge &&
      xyz[1]<maxTpcBound[1]-extraEdge);
    bool nuIsInsideZ = (xyz[2]>minTpcBound[2]+extraEdge &&
      xyz[2]<maxTpcBound[2]-extraEdge);

    bool p1IsInsideX = (prong1_xyz[0]>minTpcBound[0]+extraEdge &&
      prong1_xyz[0]<maxTpcBound[0]-extraEdge);
    bool p1IsInsideY = (prong1_xyz[1]>minTpcBound[1]+extraEdge &&
      prong1_xyz[1]<maxTpcBound[1]-extraEdge);
    bool p1IsInsideZ = (prong1_xyz[2]>minTpcBound[2]+extraEdge &&
      prong1_xyz[2]<maxTpcBound[2]-extraEdge);

    bool p2IsInsideX = (prong2_xyz[0]>minTpcBound[0]+extraEdge &&
      prong2_xyz[0]<maxTpcBound[0]-extraEdge);
    bool p2IsInsideY = (prong2_xyz[1]>minTpcBound[1]+extraEdge &&
      prong2_xyz[1]<maxTpcBound[1]-extraEdge);
    bool p2IsInsideZ = (prong2_xyz[2]>minTpcBound[2]+extraEdge &&
      prong2_xyz[2]<maxTpcBound[2]-extraEdge);

    bool nuIsInside = (nuIsInsideX && nuIsInsideY && nuIsInsideZ);
    bool p1IsInside = (p1IsInsideX && p1IsInsideY && p1IsInsideZ);
    bool p2IsInside = (p2IsInsideX && p2IsInsideY && p2IsInsideZ);


    // If vertex is inside TPC, determine channel/tick coordinates and assign them
    if (nuIsInside && p1IsInside && p2IsInside)
    {
      fIsInsideTPC = true;
      raw::ChannelID_t channel0 = geometry->NearestChannel(xyz,0);
      raw::ChannelID_t channel1 = geometry->NearestChannel(xyz,1);
      raw::ChannelID_t channel2 = geometry->NearestChannel(xyz,2);
      double tick0 = detectorProperties->ConvertXToTicks(xyz[0], 0, 0, 0);
      double tick1 = detectorProperties->ConvertXToTicks(xyz[0], 1, 0, 0);
      double tick2 = detectorProperties->ConvertXToTicks(xyz[0], 2, 0, 0);
      fChannelLoc = {(int) channel0,(int) channel1,(int) channel2};
      fTickLoc = { (float) tick0, (float) tick1, (float) tick2};

      raw::ChannelID_t prong1_channel0 = geometry->NearestChannel(prong1_xyz,0);
      raw::ChannelID_t prong1_channel1 = geometry->NearestChannel(prong1_xyz,1);
      raw::ChannelID_t prong1_channel2 = geometry->NearestChannel(prong1_xyz,2);
      double prong1_tick0 = detectorProperties->ConvertXToTicks(prong1_xyz[0], 0, 0, 0);
      double prong1_tick1 = detectorProperties->ConvertXToTicks(prong1_xyz[0], 1, 0, 0);
      double prong1_tick2 = detectorProperties->ConvertXToTicks(prong1_xyz[0], 2, 0, 0);
      raw::ChannelID_t prong2_channel0 = geometry->NearestChannel(prong2_xyz,0);
      raw::ChannelID_t prong2_channel1 = geometry->NearestChannel(prong2_xyz,1);
      raw::ChannelID_t prong2_channel2 = geometry->NearestChannel(prong2_xyz,2);
      double prong2_tick0 = detectorProperties->ConvertXToTicks(prong2_xyz[0], 0, 0, 0);
      double prong2_tick1 = detectorProperties->ConvertXToTicks(prong2_xyz[0], 1, 0, 0);
      double prong2_tick2 = detectorProperties->ConvertXToTicks(prong2_xyz[0], 2, 0, 0);

      fProngChannelLoc = {{ (int) prong1_channel0, (int) prong1_channel1, (int) prong1_channel2}, { (int) prong2_channel0, (int) prong2_channel1, (int) prong2_channel2}};
      fProngTickLoc = {{ (float) prong1_tick0, (float) prong1_tick1, (float) prong1_tick2}, { (float) prong2_tick0, (float) prong2_tick1, (float) prong2_tick2}};
      return;
    }

    // Else flag it as outside the TPC and exit function
    else
      {
        fIsInsideTPC = false;
        return;
      }
  } // END function SetDetectorCoordinates


  void DecayVertex::SetHypothesisLabels()
  {
    /* Determine PDG code and mass for the two particles in two different hypotheses.
    h1: Longest track is muon (13), shortest is pion (211)
    h2: Longest track is pion (211), shortest is muon (13)
    */
    // Possible masses
    float muonMass = 0.10566;
    float pionMass = 0.13957;

    // First determine which of the two tracks is the longest.
    int longTrackInd = -1, shortTrackInd = -1;
    if (fProngLength[0]>=fProngLength[1])
      {
        longTrackInd = 0;
        shortTrackInd = 1;
      }
    else
      {
        longTrackInd = 1;
        shortTrackInd = 0;
      }

    // Assign pdg code and masses based on hypothesis
    fProngPdgCode_h1 = {-1,-1};
    fProngPdgCode_h2 = {-1,-1};
    fProngMass_h1 = {-1.,-1.};
    fProngMass_h2 = {-1.,-1.};
    fProngPdgCode_h1[longTrackInd] = 13;
    fProngPdgCode_h1[shortTrackInd] = 211;
    fProngPdgCode_h2[longTrackInd] = 211;
    fProngPdgCode_h2[shortTrackInd] = 13;
    fProngMass_h1[longTrackInd] = muonMass;
    fProngMass_h1[shortTrackInd] = pionMass;
    fProngMass_h2[longTrackInd] = pionMass;
    fProngMass_h2[shortTrackInd] = muonMass;
  } // END function SetHypothesisLabels


  // Internal setters
  void DecayVertex::SetMomentumQuantities_ByRange()
  {
    /* Use internal attributes (like prong lengths and position) for the decay vertex to determine prong momenta and quantities determined from them (total momentum, energy, invariant mass, etc.).
    Momenta in this case are determined by range. Two hypothesis are assumed.
    h1: Longest track is muon (13), shortest is pion (211)
    h2: Longest track is pion (211), shortest is muon (13)

    The algorithms starts by calculating muon mass for both, then it assumes one of them is pion and scales momentum by mass ratio (m_pi/m_mu). It does that for shortest track in h1 and for longest track in h2.
    */

    // Assign momentum using TrackMomentumCalculator for the two hypotheses
    trkf::TrackMomentumCalculator tmc;
    // First assume both are muons, for both hypotheses (13 and 2212 are the only arguments accepted by GetTrackMomentum method)
    fProngMomMag_ByRange_h1 = {(float) tmc.GetTrackMomentum(fProngLength[0],13), (float) tmc.GetTrackMomentum(fProngLength[1],13)};
    fProngMomMag_ByRange_h2 = {(float) tmc.GetTrackMomentum(fProngLength[0],13), (float) tmc.GetTrackMomentum(fProngLength[1],13)};
    // Then scale the pion candidate momentum by mass ratio.
    float e1, e2, muonMass = 0.10566;
    for(int i=0; i<2; i++)
    {
      fProngMomMag_ByRange_h1[i] = fProngMomMag_ByRange_h1[i]*fProngMass_h1[i]/muonMass;
      fProngMomMag_ByRange_h2[i] = fProngMomMag_ByRange_h2[i]*fProngMass_h2[i]/muonMass;
    }

    // Prong momentum components
    fProngMom_ByRange_h1_X = {(float) fProngDirX[0]*fProngMomMag_ByRange_h1[0], (float) fProngDirX[1]*fProngMomMag_ByRange_h1[1]};
    fProngMom_ByRange_h1_Y = {(float) fProngDirY[0]*fProngMomMag_ByRange_h1[0], (float) fProngDirY[1]*fProngMomMag_ByRange_h1[1]};
    fProngMom_ByRange_h1_Z = {(float) fProngDirZ[0]*fProngMomMag_ByRange_h1[0], (float) fProngDirZ[1]*fProngMomMag_ByRange_h1[1]};
    fProngMom_ByRange_h2_X = {(float) fProngDirX[0]*fProngMomMag_ByRange_h2[0], (float) fProngDirX[1]*fProngMomMag_ByRange_h2[1]};
    fProngMom_ByRange_h2_Y = {(float) fProngDirY[0]*fProngMomMag_ByRange_h2[0], (float) fProngDirY[1]*fProngMomMag_ByRange_h2[1]};
    fProngMom_ByRange_h2_Z = {(float) fProngDirZ[0]*fProngMomMag_ByRange_h2[0], (float) fProngDirZ[1]*fProngMomMag_ByRange_h2[1]};

    // Prong energy
    e1 = sqrt(pow(fProngMass_h1[0],2.) + pow(fProngMomMag_ByRange_h1[0],2.));
    e2 = sqrt(pow(fProngMass_h1[1],2.) + pow(fProngMomMag_ByRange_h1[1],2.));
    fProngEnergy_ByRange_h1 = {e1, e2};
    e1 = sqrt(pow(fProngMass_h2[0],2.) + pow(fProngMomMag_ByRange_h2[0],2.));
    e2 = sqrt(pow(fProngMass_h2[1],2.) + pow(fProngMomMag_ByRange_h2[1],2.));
    fProngEnergy_ByRange_h2 = {e1, e2};
    /**/
    // Calculate kinematic quantities for parent neutrino
    // Total momentum components
    fTotMom_ByRange_h1_X = fProngMom_ByRange_h1_X[0] + fProngMom_ByRange_h1_X[1];
    fTotMom_ByRange_h1_Y = fProngMom_ByRange_h1_Y[0] + fProngMom_ByRange_h1_Y[1];
    fTotMom_ByRange_h1_Z = fProngMom_ByRange_h1_Z[0] + fProngMom_ByRange_h1_Z[1];
    fTotMom_ByRange_h2_X = fProngMom_ByRange_h2_X[0] + fProngMom_ByRange_h2_X[1];
    fTotMom_ByRange_h2_Y = fProngMom_ByRange_h2_Y[0] + fProngMom_ByRange_h2_Y[1];
    fTotMom_ByRange_h2_Z = fProngMom_ByRange_h2_Z[0] + fProngMom_ByRange_h2_Z[1];
    // Total momentum magnitude
    fTotMomMag_ByRange_h1 = sqrt(pow(fTotMom_ByRange_h1_X,2.) + pow(fTotMom_ByRange_h1_Y,2.) + pow(fTotMom_ByRange_h1_Z,2.));
    fTotMomMag_ByRange_h2 = sqrt(pow(fTotMom_ByRange_h2_X,2.) + pow(fTotMom_ByRange_h2_Y,2.) + pow(fTotMom_ByRange_h2_Z,2.));
    // Total direction components
    fTotDir_ByRange_h1_X = fTotMom_ByRange_h1_X/fTotMomMag_ByRange_h1;
    fTotDir_ByRange_h1_Y = fTotMom_ByRange_h1_Y/fTotMomMag_ByRange_h1;
    fTotDir_ByRange_h1_Z = fTotMom_ByRange_h1_Z/fTotMomMag_ByRange_h1;
    fTotDir_ByRange_h2_X = fTotMom_ByRange_h2_X/fTotMomMag_ByRange_h2;
    fTotDir_ByRange_h2_Y = fTotMom_ByRange_h2_Y/fTotMomMag_ByRange_h2;
    fTotDir_ByRange_h2_Z = fTotMom_ByRange_h2_Z/fTotMomMag_ByRange_h2;
    // Total direction angles
    fTotTheta_ByRange_h1 = acos(fTotDir_ByRange_h1_Z);
    fTotPhi_ByRange_h1 = atan2(fTotDir_ByRange_h1_Y,fTotDir_ByRange_h1_X);
    fTotTheta_ByRange_h2 = acos(fTotDir_ByRange_h2_Z);
    fTotPhi_ByRange_h2 = atan2(fTotDir_ByRange_h2_Y,fTotDir_ByRange_h2_X);
    // Total energy
    fTotEnergy_ByRange_h1 = fProngEnergy_ByRange_h1[0] + fProngEnergy_ByRange_h1[1];
    fTotEnergy_ByRange_h2 = fProngEnergy_ByRange_h2[0] + fProngEnergy_ByRange_h2[1];
    // Invariant mass
    fInvMass_ByRange_h1 = sqrt(pow(fTotEnergy_ByRange_h1,2.) - pow(fTotMomMag_ByRange_h1,2.));
    fInvMass_ByRange_h2 = sqrt(pow(fTotEnergy_ByRange_h2,2.) - pow(fTotMomMag_ByRange_h2,2.));

    // Calculate opening angle
    std::vector<double> startDirection1 = {
      fProngTrack[0]->StartDirection().X(),
      fProngTrack[0]->StartDirection().Y(),
      fProngTrack[0]->StartDirection().Z()
    };
    std::vector<double> startDirection2 = {
      fProngTrack[1]->StartDirection().X(),
      fProngTrack[1]->StartDirection().Y(),
      fProngTrack[1]->StartDirection().Z()
    };
    float magnitude1 = sqrt(startDirection1[0]*startDirection1[0] + startDirection1[1]*startDirection1[1] + startDirection1[2]*startDirection1[2]);
    float magnitude2 = sqrt(startDirection2[0]*startDirection2[0] + startDirection2[1]*startDirection2[1] + startDirection2[2]*startDirection2[2]);
    float dotProduct = startDirection1[0]*startDirection2[0] + startDirection1[1]*startDirection2[1] + startDirection1[2]*startDirection2[2];
    fOpeningAngle = acos(dotProduct / (magnitude1*magnitude2));
    // Calculate start point to neutrino vertex distance
    float prong1_distance = sqrt(pow(fX - fProngX[0],2.) + pow(fY - fProngY[0],2.) + pow(fZ - fProngZ[0],2.));
    float prong2_distance = sqrt(pow(fX - fProngX[1],2.) + pow(fY - fProngY[1],2.) + pow(fZ - fProngZ[1],2.));
    fProngStartToNeutrinoDistance = {prong1_distance,prong2_distance};
    return;
  } // END function SetMomentumQuatities_ByRange


  void DecayVertex::SetMomentumQuantities_ByMCS()
  {
    /* Use internal attributes (like prong lengths and position) for the decay vertex to determine prong momenta and quantities determined from them (total momentum, energy, invariant mass, etc.).
    Momenta in this case are determined by Multiple Coulomb Scattering. Two hypothesis are assumed.
    fwd: Forward fit is used for momentum
    best: Best fit between forward and backward is used for particle momentum

    The algorithms starts by calculating muon mass for both, then it assumes one of them is pion and scales momentum by mass ratio (m_pi/m_mu). It does that for shortest track in h1 and for longest track in h2.
    */

    float e1, e2;//, muonMass = 0.10566;
    fProngPdgCodeHypothesis_ByMcs = {fProngMcs[0]->particleIdHyp(),fProngMcs[1]->particleIdHyp()};
    fProngIsBestFwd_ByMcs = {fProngMcs[0]->isBestFwd(),fProngMcs[1]->isBestFwd()};
    // Prong momentum magnitude
    fProngMomMag_ByMcs_fwd_h1 = {(float) fProngMcs[0]->fwdMomentum(), (float) fProngMcs[1]->fwdMomentum()};
    fProngMomMag_ByMcs_best_h1 = {(float) fProngMcs[0]->bestMomentum(), (float) fProngMcs[1]->bestMomentum()};
    fProngMomMag_ByMcs_fwd_h2 = {(float) fProngMcs[0]->fwdMomentum(), (float) fProngMcs[1]->fwdMomentum()};
    fProngMomMag_ByMcs_best_h2 = {(float) fProngMcs[0]->bestMomentum(), (float) fProngMcs[1]->bestMomentum()};
    // Prong momentum components
    fProngMom_ByMcs_fwd_h1_X = {fProngDirX[0]*fProngMomMag_ByMcs_fwd_h1[0],fProngDirX[1]*fProngMomMag_ByMcs_fwd_h1[1]};
    fProngMom_ByMcs_fwd_h1_Y = {fProngDirY[0]*fProngMomMag_ByMcs_fwd_h1[0],fProngDirY[1]*fProngMomMag_ByMcs_fwd_h1[1]};
    fProngMom_ByMcs_fwd_h1_Z = {fProngDirZ[0]*fProngMomMag_ByMcs_fwd_h1[0],fProngDirZ[1]*fProngMomMag_ByMcs_fwd_h1[1]};
    /**/
    fProngMom_ByMcs_best_h1_X = {fProngDirX[0]*fProngMomMag_ByMcs_best_h1[0],fProngDirX[1]*fProngMomMag_ByMcs_best_h1[1]};
    fProngMom_ByMcs_best_h1_Y = {fProngDirY[0]*fProngMomMag_ByMcs_best_h1[0],fProngDirY[1]*fProngMomMag_ByMcs_best_h1[1]};
    fProngMom_ByMcs_best_h1_Z = {fProngDirZ[0]*fProngMomMag_ByMcs_best_h1[0],fProngDirZ[1]*fProngMomMag_ByMcs_best_h1[1]};
    /**/
    fProngMom_ByMcs_fwd_h2_X = {fProngDirX[0]*fProngMomMag_ByMcs_fwd_h2[0],fProngDirX[1]*fProngMomMag_ByMcs_fwd_h2[1]};
    fProngMom_ByMcs_fwd_h2_Y = {fProngDirY[0]*fProngMomMag_ByMcs_fwd_h2[0],fProngDirY[1]*fProngMomMag_ByMcs_fwd_h2[1]};
    fProngMom_ByMcs_fwd_h2_Z = {fProngDirZ[0]*fProngMomMag_ByMcs_fwd_h2[0],fProngDirZ[1]*fProngMomMag_ByMcs_fwd_h2[1]};
    /**/
    fProngMom_ByMcs_best_h2_X = {fProngDirX[0]*fProngMomMag_ByMcs_best_h2[0],fProngDirX[1]*fProngMomMag_ByMcs_best_h2[1]};
    fProngMom_ByMcs_best_h2_Y = {fProngDirY[0]*fProngMomMag_ByMcs_best_h2[0],fProngDirY[1]*fProngMomMag_ByMcs_best_h2[1]};
    fProngMom_ByMcs_best_h2_Z = {fProngDirZ[0]*fProngMomMag_ByMcs_best_h2[0],fProngDirZ[1]*fProngMomMag_ByMcs_best_h2[1]};
    /**/
    // Prong energy
    e1 = sqrt(pow(fProngMass_h1[0],2.) + pow(fProngMomMag_ByMcs_fwd_h1[0],2.));
    e2 = sqrt(pow(fProngMass_h1[1],2.) + pow(fProngMomMag_ByMcs_fwd_h1[1],2.));
    fProngEnergy_ByMcs_fwd_h1 = {e1, e2};
    e1 = sqrt(pow(fProngMass_h1[0],2.) + pow(fProngMomMag_ByMcs_best_h1[0],2.));
    e2 = sqrt(pow(fProngMass_h1[1],2.) + pow(fProngMomMag_ByMcs_best_h1[1],2.));
    fProngEnergy_ByMcs_best_h1 = {e1, e2};
    e1 = sqrt(pow(fProngMass_h2[0],2.) + pow(fProngMomMag_ByMcs_fwd_h2[0],2.));
    e2 = sqrt(pow(fProngMass_h2[1],2.) + pow(fProngMomMag_ByMcs_fwd_h2[1],2.));
    fProngEnergy_ByMcs_fwd_h2 = {e1, e2};
    e1 = sqrt(pow(fProngMass_h2[0],2.) + pow(fProngMomMag_ByMcs_best_h2[0],2.));
    e2 = sqrt(pow(fProngMass_h2[1],2.) + pow(fProngMomMag_ByMcs_best_h2[1],2.));
    fProngEnergy_ByMcs_best_h2 = {e1, e2};
    /**/
    // Calculate kinematic quantities for parent neutrino
    // Total momentum components
    fTotMom_ByMcs_fwd_h1_X = fProngMom_ByMcs_fwd_h1_X[0] + fProngMom_ByMcs_fwd_h1_X[1];
    fTotMom_ByMcs_fwd_h1_Y = fProngMom_ByMcs_fwd_h1_Y[0] + fProngMom_ByMcs_fwd_h1_Y[1];
    fTotMom_ByMcs_fwd_h1_Z = fProngMom_ByMcs_fwd_h1_Z[0] + fProngMom_ByMcs_fwd_h1_Z[1];
    /**/
    fTotMom_ByMcs_best_h1_X = fProngMom_ByMcs_best_h1_X[0] + fProngMom_ByMcs_best_h1_X[1];
    fTotMom_ByMcs_best_h1_Y = fProngMom_ByMcs_best_h1_Y[0] + fProngMom_ByMcs_best_h1_Y[1];
    fTotMom_ByMcs_best_h1_Z = fProngMom_ByMcs_best_h1_Z[0] + fProngMom_ByMcs_best_h1_Z[1];
    /**/
    fTotMom_ByMcs_fwd_h2_X = fProngMom_ByMcs_fwd_h2_X[0] + fProngMom_ByMcs_fwd_h2_X[1];
    fTotMom_ByMcs_fwd_h2_Y = fProngMom_ByMcs_fwd_h2_Y[0] + fProngMom_ByMcs_fwd_h2_Y[1];
    fTotMom_ByMcs_fwd_h2_Z = fProngMom_ByMcs_fwd_h2_Z[0] + fProngMom_ByMcs_fwd_h2_Z[1];
    /**/
    fTotMom_ByMcs_best_h2_X = fProngMom_ByMcs_best_h2_X[0] + fProngMom_ByMcs_best_h2_X[1];
    fTotMom_ByMcs_best_h2_Y = fProngMom_ByMcs_best_h2_Y[0] + fProngMom_ByMcs_best_h2_Y[1];
    fTotMom_ByMcs_best_h2_Z = fProngMom_ByMcs_best_h2_Z[0] + fProngMom_ByMcs_best_h2_Z[1];
    // Total momentum magnitude
    fTotMomMag_ByMcs_fwd_h1 = sqrt(pow(fTotMom_ByMcs_fwd_h1_X,2.) + pow(fTotMom_ByMcs_fwd_h1_Y,2.) + pow(fTotMom_ByMcs_fwd_h1_Z,2.));
    fTotMomMag_ByMcs_best_h1 = sqrt(pow(fTotMom_ByMcs_best_h1_X,2.) + pow(fTotMom_ByMcs_best_h1_Y,2.) + pow(fTotMom_ByMcs_best_h1_Z,2.));
    fTotMomMag_ByMcs_fwd_h2 = sqrt(pow(fTotMom_ByMcs_fwd_h2_X,2.) + pow(fTotMom_ByMcs_fwd_h2_Y,2.) + pow(fTotMom_ByMcs_fwd_h2_Z,2.));
    fTotMomMag_ByMcs_best_h2 = sqrt(pow(fTotMom_ByMcs_best_h2_X,2.) + pow(fTotMom_ByMcs_best_h2_Y,2.) + pow(fTotMom_ByMcs_best_h2_Z,2.));
    // Total direction components
    fTotDir_ByMcs_fwd_h1_X = fTotMom_ByMcs_fwd_h1_X/fTotMomMag_ByMcs_fwd_h1;
    fTotDir_ByMcs_fwd_h1_Y = fTotMom_ByMcs_fwd_h1_Y/fTotMomMag_ByMcs_fwd_h1;
    fTotDir_ByMcs_fwd_h1_Z = fTotMom_ByMcs_fwd_h1_Z/fTotMomMag_ByMcs_fwd_h1;
    /**/
    fTotDir_ByMcs_best_h1_X = fTotMom_ByMcs_best_h1_X/fTotMomMag_ByMcs_best_h1;
    fTotDir_ByMcs_best_h1_Y = fTotMom_ByMcs_best_h1_Y/fTotMomMag_ByMcs_best_h1;
    fTotDir_ByMcs_best_h1_Z = fTotMom_ByMcs_best_h1_Z/fTotMomMag_ByMcs_best_h1;
    /**/
    fTotDir_ByMcs_fwd_h2_X = fTotMom_ByMcs_fwd_h2_X/fTotMomMag_ByMcs_fwd_h2;
    fTotDir_ByMcs_fwd_h2_Y = fTotMom_ByMcs_fwd_h2_Y/fTotMomMag_ByMcs_fwd_h2;
    fTotDir_ByMcs_fwd_h2_Z = fTotMom_ByMcs_fwd_h2_Z/fTotMomMag_ByMcs_fwd_h2;
    /**/
    fTotDir_ByMcs_best_h2_X = fTotMom_ByMcs_best_h2_X/fTotMomMag_ByMcs_best_h2;
    fTotDir_ByMcs_best_h2_Y = fTotMom_ByMcs_best_h2_Y/fTotMomMag_ByMcs_best_h2;
    fTotDir_ByMcs_best_h2_Z = fTotMom_ByMcs_best_h2_Z/fTotMomMag_ByMcs_best_h2;
    // Total direction angles
    fTotTheta_ByMcs_fwd_h1 = acos(fTotDir_ByMcs_fwd_h1_Z);
    fTotPhi_ByMcs_fwd_h1 = atan2(fTotDir_ByMcs_fwd_h1_Y,fTotDir_ByMcs_fwd_h1_X);
    /**/
    fTotTheta_ByMcs_best_h1 = acos(fTotDir_ByMcs_best_h1_Z);
    fTotPhi_ByMcs_best_h1 = atan2(fTotDir_ByMcs_best_h1_Y,fTotDir_ByMcs_best_h1_X);
    /**/
    fTotTheta_ByMcs_fwd_h2 = acos(fTotDir_ByMcs_fwd_h2_Z);
    fTotPhi_ByMcs_fwd_h2 = atan2(fTotDir_ByMcs_fwd_h2_Y,fTotDir_ByMcs_fwd_h2_X);
    /**/
    fTotTheta_ByMcs_best_h2 = acos(fTotDir_ByMcs_best_h2_Z);
    fTotPhi_ByMcs_best_h2 = atan2(fTotDir_ByMcs_best_h2_Y,fTotDir_ByMcs_best_h2_X);
    // Total energy
    fTotEnergy_ByMcs_fwd_h1 = fProngEnergy_ByMcs_fwd_h1[0] + fProngEnergy_ByMcs_fwd_h1[1];
    fTotEnergy_ByMcs_best_h1 = fProngEnergy_ByMcs_best_h1[0] + fProngEnergy_ByMcs_best_h1[1];
    fTotEnergy_ByMcs_fwd_h2 = fProngEnergy_ByMcs_fwd_h2[0] + fProngEnergy_ByMcs_fwd_h2[1];
    fTotEnergy_ByMcs_best_h2 = fProngEnergy_ByMcs_best_h2[0] + fProngEnergy_ByMcs_best_h2[1];
    // Invariant mass
    fInvMass_ByMcs_fwd_h1 = sqrt(pow(fTotEnergy_ByMcs_fwd_h1,2.) - pow(fTotMomMag_ByMcs_fwd_h1,2.));
    fInvMass_ByMcs_best_h1 = sqrt(pow(fTotEnergy_ByMcs_best_h1,2.) - pow(fTotMomMag_ByMcs_best_h1,2.));
    fInvMass_ByMcs_fwd_h2 = sqrt(pow(fTotEnergy_ByMcs_fwd_h2,2.) - pow(fTotMomMag_ByMcs_fwd_h2,2.));
    fInvMass_ByMcs_best_h2 = sqrt(pow(fTotEnergy_ByMcs_best_h2,2.) - pow(fTotMomMag_ByMcs_best_h2,2.));

    // Calculate opening angle
    std::vector<double> startDirection1 = {
      fProngTrack[0]->StartDirection().X(),
      fProngTrack[0]->StartDirection().Y(),
      fProngTrack[0]->StartDirection().Z()
    };
    std::vector<double> startDirection2 = {
      fProngTrack[1]->StartDirection().X(),
      fProngTrack[1]->StartDirection().Y(),
      fProngTrack[1]->StartDirection().Z()
    };
    float magnitude1 = sqrt(startDirection1[0]*startDirection1[0] + startDirection1[1]*startDirection1[1] + startDirection1[2]*startDirection1[2]);
    float magnitude2 = sqrt(startDirection2[0]*startDirection2[0] + startDirection2[1]*startDirection2[1] + startDirection2[2]*startDirection2[2]);
    float dotProduct = startDirection1[0]*startDirection2[0] + startDirection1[1]*startDirection2[1] + startDirection1[2]*startDirection2[2];
    fOpeningAngle = acos(dotProduct / (magnitude1*magnitude2));
    // Calculate start point to neutrino vertex distance
    float prong1_distance = sqrt(pow(fX - fProngX[0],2.) + pow(fY - fProngY[0],2.) + pow(fZ - fProngZ[0],2.));
    float prong2_distance = sqrt(pow(fX - fProngX[1],2.) + pow(fY - fProngY[1],2.) + pow(fZ - fProngZ[1],2.));
    fProngStartToNeutrinoDistance = {prong1_distance,prong2_distance};



    /**/
    return;
  } // END function SetMomentumQuantities_ByMCS




  // Printers
  void DecayVertex::PrintInformation() const
  {
    int fStartWire[3] = {0,2399,4798};
    printf("\n-|Vertex information|\n");
    printf("|_Vertex inside TPC: %d\n", fIsInsideTPC);
    printf("|_Spatial Coordinates: [%.1f, %.1f, %.1f]\n",fX,fY,fZ);
    printf("|_Prongs lengths: [%.1f,%.1f]\n", fProngLength[0], fProngLength[1]);
    printf("|_Prongs theta: [%.1f,%.1f]\n", fProngTheta[0], fProngTheta[1]);
    printf("|_Prongs phi: [%.1f,%.1f]\n", fProngPhi[0], fProngPhi[1]);
    printf("|_Prongs hit number: [%i,%i]\n", fProngNumHits[0], fProngNumHits[1]);
    printf("|_Detector location assigned: %d\n", fIsDetLocAssigned);
    if (fIsDetLocAssigned)
    {
      printf("|_Channel coordinates: [%i, %i, %i]\n",fChannelLoc[0]-fStartWire[0],fChannelLoc[1]-fStartWire[1],fChannelLoc[2]-fStartWire[2]);
      printf("|_Ticks coordinates: [%.1f, %.1f, %.1f]\n",fTickLoc[0],fTickLoc[1],fTickLoc[2]);
    }
    return;
  }

} // END namespace AuxVertex

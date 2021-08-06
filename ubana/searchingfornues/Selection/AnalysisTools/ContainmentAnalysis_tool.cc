#ifndef SELECTION_CONTAINMENTSELECTION_CXX
#define SELECTION_CONTAINMENTSELECTION_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// backtracking tools
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "larcore/Geometry/Geometry.h"

namespace analysis
{
////////////////////////////////////////////////////////////////////////
//
// Class:       ContainmentAnalysis
// File:        ContainmentAnalysis.cc
//
//              A basic selection example
//
// Configuration parameters:
//
// TBD
//
// Created by David Caratelli (davidc@fnal.gov) on 01/30/2019
//
////////////////////////////////////////////////////////////////////////

class ContainmentAnalysis : public AnalysisToolBase
{

public:
  /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
  ContainmentAnalysis(const fhicl::ParameterSet &pset);

  /**
     *  @brief  Destructor
     */
  ~ContainmentAnalysis(){};

  // provide for initialization
  void configure(fhicl::ParameterSet const &pset);

  /**
     * @brief Analysis function
     */
  void analyzeEvent(art::Event const &e, bool fData) override;

  /**
     * @brief Selection function
     */
  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  /**
     * @brief set branches for TTree
     */
  void setBranches(TTree *_tree) override;

  /**
     * @brief reset ttree branches
     */
  void resetTTree(TTree *_tree) override;

private:
  float DistFiducial(float x, float y, float z);
  void DistFiducialBoundaries(float x, float y, float z, std::vector<std::vector<double>> &dboundaries);
  bool isFiducial(const float x[3]) const;

  art::InputTag fMCTproducer;

  float _FV; // FV boundary to apply

  // TTree variables
  float _dvtx; // smallest distance between vertex and any boundary
  float _dtrk; // smallest distance between any track start/end point and any boundary

  std::vector<std::vector<double>> _dmc_boundary;  // smallest distance between any MCParticle and each boundary
  std::vector<std::vector<double>> _dtrk_boundary; // smallest distance between any track start/end point and each boundary
  std::vector<double> _dtrk_x_boundary;            // smallest distance between any track start/end point and x boundaries
  std::vector<double> _dtrk_y_boundary;            // smallest distance between any track start/end point and y boundaries
  std::vector<double> _dtrk_z_boundary;            // smallest distance between any track start/end point and z boundaries

  std::vector<double> _dshr_x_boundary; // smallest distance between any shower start point and x boundaries
  std::vector<double> _dshr_y_boundary; // smallest distance between any shower start point and y boundaries
  std::vector<double> _dshr_z_boundary; // smallest distance between any shower start point and z boundaries

  std::vector<double> _dvtx_x_boundary; // smallest distance between the neutrino vertex and x boundaries
  std::vector<double> _dvtx_y_boundary; // smallest distance between the neutrino vertex and y boundaries
  std::vector<double> _dvtx_z_boundary; // smallest distance between the neutrino vertex and z boundaries

  std::vector<std::vector<double>> _dvtx_boundary; // smallest distance between the neutrino vertex and each boundary
  std::vector<std::vector<double>> _dshr_boundary; // smallest distance between any shower start point and each boundary

  float contained_sps_ratio;
};

//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
ContainmentAnalysis::ContainmentAnalysis(const fhicl::ParameterSet &pset)
{
  configure(pset);
}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void ContainmentAnalysis::configure(fhicl::ParameterSet const &pset)
{
  fMCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
  _FV = pset.get<float>("FV");
}

void ContainmentAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  if (!fData)
  {
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    auto mct = mct_h->at(0);
    size_t npart = mct.NParticles();
    _dmc_boundary.clear();
    _dmc_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));
    for (size_t i = 0; i < npart; i++)
    {
      auto const &part = mct.GetParticle(i);
      if (part.StatusCode() != 1)
      {
        continue;
      }
      DistFiducialBoundaries(part.Vx(), part.Vy(), part.Vz(), _dmc_boundary);
      DistFiducialBoundaries(part.EndX(), part.EndY(), part.EndZ(), _dmc_boundary);
    }
  }
}

//----------------------------------------------------------------------------
/// selectEvent
///
/// Arguments:
///
/// art::Event
/// slice track pointer vector
/// slice shower pointer vector
///
void ContainmentAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{

  TVector3 nuvtx;
  Double_t xyz[3] = {};

  _dvtx = 1e3;
  _dtrk = 1e3;

  _dtrk_x_boundary.clear();
  _dtrk_y_boundary.clear();
  _dtrk_z_boundary.clear();
  _dshr_x_boundary.clear();
  _dshr_y_boundary.clear();
  _dshr_z_boundary.clear();
  _dvtx_x_boundary.clear();
  _dvtx_y_boundary.clear();
  _dvtx_z_boundary.clear();

  _dtrk_boundary.clear();
  _dvtx_boundary.clear();
  _dshr_boundary.clear();

  _dtrk_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));
  _dvtx_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));
  _dshr_boundary.resize(3, std::vector<double>(2, std::numeric_limits<double>::max()));

  size_t sps_fv = 0, sps_all = 0;

  for (const auto &pfp_pxy : slice_pfp_v)
  {

    auto PDG = fabs(pfp_pxy->PdgCode());

    if ((PDG == 12) || (PDG == 14))
    {

      // grab vertex
      auto vtx = pfp_pxy.get<recob::Vertex>();
      if (vtx.size() != 1)
      {
        std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
        return;
      }

      // save vertex to array
      vtx.at(0)->XYZ(xyz);
      nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

      _dvtx = DistFiducial(nuvtx.X(), nuvtx.Y(), nuvtx.Z());
      DistFiducialBoundaries(xyz[0], xyz[1], xyz[2], _dvtx_boundary);
      _dvtx_x_boundary = _dvtx_boundary[0];
      _dvtx_y_boundary = _dvtx_boundary[1];
      _dvtx_z_boundary = _dvtx_boundary[2];
    } // if neutrino PFP

    else
    { // if not the neutrino PFP

      auto spcpnts = pfp_pxy.get<recob::SpacePoint>();
      sps_all += spcpnts.size();
      for (auto &sp : spcpnts)
      {
        float _reco_nu_vtx_sce[3];

        _reco_nu_vtx_sce[0] = sp->XYZ()[0];
        _reco_nu_vtx_sce[1] = sp->XYZ()[1];
        _reco_nu_vtx_sce[2] = sp->XYZ()[2];
        searchingfornues::ApplySCECorrectionXYZ(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2], _reco_nu_vtx_sce);
        if (isFiducial(_reco_nu_vtx_sce))
          sps_fv++;
      }

      auto ntrk = pfp_pxy.get<recob::Track>().size();

      if (ntrk == 1)
      {

        auto trk = pfp_pxy.get<recob::Track>().at(0);

        auto trkstart = trk->Vertex();
        auto trkend = trk->End();

        float dstrt = DistFiducial(trkstart.X(), trkstart.Y(), trkstart.Z());
        if (dstrt < _dtrk)
          _dtrk = dstrt;

        float dend = DistFiducial(trkend.X(), trkend.Y(), trkend.Z());
        if (dend < _dtrk)
          _dtrk = dend;

        DistFiducialBoundaries(trkend.X(), trkend.Y(), trkend.Z(), _dtrk_boundary);
        DistFiducialBoundaries(trkstart.X(), trkstart.Y(), trkstart.Z(), _dtrk_boundary);
        _dtrk_x_boundary = _dtrk_boundary[0];
        _dtrk_y_boundary = _dtrk_boundary[1];
        _dtrk_z_boundary = _dtrk_boundary[2];
      }

      auto nshr = pfp_pxy.get<recob::Shower>().size();

      if (nshr == 1)
      {
        auto shr = pfp_pxy.get<recob::Shower>().at(0);
        auto shrstart = shr->ShowerStart();
        DistFiducialBoundaries(shrstart.X(), shrstart.Y(), shrstart.Z(), _dshr_boundary);
        _dshr_x_boundary = _dshr_boundary[0];
        _dshr_y_boundary = _dshr_boundary[1];
        _dshr_z_boundary = _dshr_boundary[2];
      }

    } // if not the neutrino PFP

  } // for all PFP
  contained_sps_ratio = sps_fv / ((float) sps_all);
}

void ContainmentAnalysis::setBranches(TTree *_tree)
{

  _tree->Branch("dvtx", &_dvtx, "dvtx/F");
  _tree->Branch("dtrk", &_dtrk, "dtrk/F");
  _tree->Branch("contained_sps_ratio", &contained_sps_ratio, "contained_sps_ratio/F");

  _tree->Branch("dtrk_x_boundary", "std::vector < double >", &_dtrk_x_boundary);
  _tree->Branch("dtrk_y_boundary", "std::vector < double >", &_dtrk_y_boundary);
  _tree->Branch("dtrk_z_boundary", "std::vector < double >", &_dtrk_z_boundary);
  _tree->Branch("dshr_x_boundary", "std::vector < double >", &_dshr_x_boundary);
  _tree->Branch("dshr_y_boundary", "std::vector < double >", &_dshr_y_boundary);
  _tree->Branch("dshr_z_boundary", "std::vector < double >", &_dshr_z_boundary);
  _tree->Branch("dvtx_x_boundary", "std::vector < double >", &_dvtx_x_boundary);
  _tree->Branch("dvtx_y_boundary", "std::vector < double >", &_dvtx_y_boundary);
  _tree->Branch("dvtx_z_boundary", "std::vector < double >", &_dvtx_z_boundary);

  _tree->Branch("dtrk_boundary", "std::vector < std::vector < double > >", &_dtrk_boundary);
  _tree->Branch("dvtx_boundary", "std::vector < std::vector < double > >", &_dvtx_boundary);
  _tree->Branch("dshr_boundary", "std::vector < std::vector < double > >", &_dshr_boundary);
  _tree->Branch("dmc_boundary", "std::vector < std::vector < double > >", &_dmc_boundary);

  return;
}

void ContainmentAnalysis::resetTTree(TTree *_tree)
{

  _dvtx = std::numeric_limits<int>::max();
  _dtrk = std::numeric_limits<int>::max();

  return;
}

bool ContainmentAnalysis::isFiducial(const float x[3]) const
{
  float border = 10.;
  art::ServiceHandle<geo::Geometry> geo;
  geo::TPCGeo const &thisTPC = geo->TPC();
  geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
  std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
  bool is_x =
      x[0] > (bnd[0] + border) && x[0] < (bnd[1] - border);
  bool is_y =
      x[1] > (bnd[2] + border) && x[1] < (bnd[3] - border);
  bool is_z =
      x[2] > (bnd[4] + border) && x[2] < (bnd[5] - border);

  return is_x && is_y && is_z;
}

void ContainmentAnalysis::DistFiducialBoundaries(float x, float y, float z, std::vector<std::vector<double>> &dist_boundaries)
{
  art::ServiceHandle<geo::Geometry> geo;
  geo::TPCGeo const &thisTPC = geo->TPC();
  geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
  std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  if (SCE->EnableCalSpatialSCE() == true)
  {

    auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
    x += offset.X();
    y -= offset.Y();
    z -= offset.Z();

  } // if spatial offset calibrations are enabled

  if (x - bnd[0] < dist_boundaries[0][0])
  {
    dist_boundaries[0][0] = x - bnd[0];
  }

  if (bnd[1] - x < dist_boundaries[0][1])
  {
    dist_boundaries[0][1] = bnd[1] - x;
  }

  if (y - bnd[2] < dist_boundaries[1][0])
  {
    dist_boundaries[1][0] = y - bnd[2];
  }

  if (bnd[3] - y < dist_boundaries[1][1])
  {

    dist_boundaries[1][1] = bnd[3] - y;
  }

  if (z - bnd[4] < dist_boundaries[2][0])
  {
    dist_boundaries[2][0] = z - bnd[4];
  }

  if (bnd[5] - z < dist_boundaries[2][1])
  {
    dist_boundaries[2][1] = bnd[5] - z;
  }
}

float ContainmentAnalysis::DistFiducial(float x, float y, float z)
{

  float dmin = 1e3;

  auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

  if (SCE->EnableCalSpatialSCE() == true)
  {

    auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
    x += offset.X();
    y -= offset.Y();
    z -= offset.Z();

  } // if spatial offset calibrations are enabled

  double dxl = x - 0.;
  if (dxl < dmin)
    dmin = dxl;

  double dxh = 256. - x;
  if (dxh < dmin)
    dmin = dxh;

  double dyl = y - -116.;
  if (dyl < dmin)
    dmin = dyl;

  double dyh = 116. - y;
  if (dyh < dmin)
    dmin = dyh;

  double dzl = z - 0.;
  if (dzl < dmin)
    dmin = dzl;

  double dzh = 1036. - z;
  if (dzh < dmin)
    dmin = dzh;

  return dmin;
}

DEFINE_ART_CLASS_TOOL(ContainmentAnalysis)
} // namespace analysis

#endif

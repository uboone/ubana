/**
 *  @file   NuMuCCSelectionIIAlg.h
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices
 * 
 *  @authors xiao.luo@yale.edu, tjyang@fnal.gov
 */
#ifndef NuMuCCSelectionIIAlg_h
#define NuMuCCSelectionIIAlg_h

#include "ubana/TPCNeutrinoIDFilter/Algorithms/NeutrinoIDAlgBase.h"

// LArSoft includes
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/Track.h"

// Root includes
#include "TH1D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{

/**
 *  @brief  NuMuCCSelectionIIAlg class
 */
class NuMuCCSelectionIIAlg : virtual public NeutrinoIDAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    NuMuCCSelectionIIAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~NuMuCCSelectionIIAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const&);
    
    /**
     *  @brief Set up for "beginJob" phase if requested
     */
    virtual void beginJob(art::ServiceHandle<art::TFileService>&);
    
    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::ProducesCollector&);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    virtual bool findNeutrinoCandidates(art::Event&) const;

private:
    
    bool   inFV(double x, double y, double z) const;
    
    double scaledEdx(double x, int plane, bool isdata) const;

    //double MaxDeflection(recob::Track track) ;

    double TrunMean(std::vector<double> dqdx_v) const;

    double Median(std::vector<double> input) const;
    
    /**
     *  @ brief FHICL parameters.
     */
    std::string                fTrackModuleLabel;        ///< Producer of input tracks
    std::string                fVertexModuleLabel;       ///< Producer of input vertices
    std::string                fOpFlashModuleLabel;      ///< Producer of flashes
    std::string                fCalorimetryModuleLabel;  ///< Producer of calorimetry module

    double                     fDistToEdgeX;             ///< fiducial volume - x
    double                     fDistToEdgeY;             ///< fiducial volume - y
    double                     fDistToEdgeZ;             ///< fiducial volume - z
    
    double                     fBeamMin;                 ///< Cut on min beam time
    double                     fBeamMax;                 ///< Cut on max beam time
    double                     fPEThresh;                ///< Cut on PE threshold
    double                     fTrk2FlashDist;           ///< Cut on track to flash distance
    double                     fMinTrk2VtxDist;          ///< Minimum track to vertex distance
    double                     fMinTrackLen;             ///< Minimum track length
    double                     fMaxCosineAngle;          ///< Cut on cosine angle of two longest track
    double                     fMaxCosy1stTrk;           ///< Maximum cosy of longtest track when mult>1
    double                     fMinTrackLen2ndTrk;       ///< Minimum track length of the second longest track when mult>1
    double                     fMaxCosySingle;           ///< Maximum cosy of the single track
    double                     fMinTrackLenSingle;       ///< Minimum track length of single track
    double                     fMindEdxRatioSingle;      ///< Minimum dEdx ratio of single track
    double                     fMaxTrkLengthySingle;     ///< Maximum track length in y projection
    double                     fMinStartdEdx1stTrk;      ///< Minimum dEdx of track start
    double                     fMaxEnddEdx1stTrk;        ///< Maximum dEdx of track end
   
    bool                       fIncludeMIPCuts;          ///< MIP Consistency
    bool                       fUseBNB;                  ///< MIP Consistency
    double                     fDeflection;              ///< MIP Consistency
    double                     fMIPLength;               ///< MIP Consistency
    double                     fMIPdQdx ;                ///< MIP Consistency
    std::vector<double>        fGainCorrections;         ///< Inherited Calibrated Gain Corrections

    bool                       fDoHists;                 ///< Fill histograms
    int                        fDebug;                   ///< Print out debug information
    //TH1D*                      fNFlashPerEvent;          ///< number of flashes per event
    //TH1D*                      fFlashPE;                 ///< flash photoelectrons
    //TH1D*                      fFlashTime;               ///< flash timing
    
    /// @{
    /**
     *  @brief Standard useful properties
     */
    geo::GeometryCore const*            fGeometry;           ///< pointer to the Geometry service
    /// @}
};

} // namespace lar_cluster3d
#endif

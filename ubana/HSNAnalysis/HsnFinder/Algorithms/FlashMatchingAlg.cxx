#include "FlashMatchingAlg.h"

namespace FlashMatching
{
  // Constructor/destructor
  FlashMatchingAlg::FlashMatchingAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
  }
  FlashMatchingAlg::~FlashMatchingAlg()
  {}
  void FlashMatchingAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fFlashLabel = pset.get<std::string>("FlashLabel");
  }

  // Find 2d-distance between vertex and largest flash
  void FlashMatchingAlg::AddFlashMatching(
            art::Event const & evt,
            AuxEvent::CandidateTreeFiller & ctf,
            AuxVertex::DecayVertex const & dv)
  {
    // Get current vertex projection on pmt plane
    float vY = dv.fY;
    float vZ = dv.fZ;

    // Prepare flash information
    art::InputTag flashTag {fFlashLabel};
    const auto& flashHandle = evt.getValidHandle< std::vector<recob::OpFlash> >(flashTag);
    if ((*flashHandle).size() > 0)
    {
      // Find largest flash
      float maxPE = 0;
      float maxPE_ind = 0;
      for(std::vector<int>::size_type i=0; i!=(*flashHandle).size(); i++)
      {
        art::Ptr<recob::OpFlash> flash(flashHandle,i);
        float flashPE = flash->TotalPE();
        if (flashPE>maxPE)
        {
          maxPE = flashPE;
          maxPE_ind = i;
        }
      }
      // Retrieve largest flash coordinates
      art::Ptr<recob::OpFlash> largestFlash(flashHandle,maxPE_ind);
      float flashY = largestFlash->YCenter();
      float flashZ = largestFlash->ZCenter();

      // Determine distance
      float flashDistance = Distance_2d(vY,flashY,vZ,flashZ);
      ctf.flash_flashDistance = flashDistance;
    }
    else
    {
      ctf.flash_flashDistance = -999.;
    }
  } //  END function AddFlashMatching

  float FlashMatchingAlg::Distance_2d(float x1,float x2,float y1,float y2)
  {
    return sqrt(pow(x1-x2,2.) + pow(y1-y2,2.));
  }
} // END namespace FlashMatching

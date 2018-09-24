#ifndef LREVTRECONSTRUCTION_CXX
#define LREVTRECONSTRUCTION_CXX

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TMath.h"

#include <iostream>

#include "uboone/NCElastic/Algorithms/LREvtReconstruction.h"

namespace qeselection {

  LREvtReconstruction::LREvtReconstruction()
  {
    _nendhits      = 6;

    _beamwinstartt = 3.2;
    _beamwinendt   = 4.8;
    _minflashpe    = 6.5;

    _minplength = 2.5;
    _minpscore  = 0.5;

    _fid_xmin   = 3.45;
    _fid_xmax   = 249.8;
    _fid_ymin   = -100.53;
    _fid_ymax   = 102.47;
    _fid_zmin   = 5.1;
    _fid_zmax   = 1031.9;

    _minmuscore = 0.5;
    _minpiscore = 0.5;

  }

  void LREvtReconstruction::Configure(fhicl::ParameterSet const& p)
  {
    _nendhits      = p.get<unsigned int>("NEndHits");

    _beamwinstartt = p.get<float>("BeamWindowStartT");
    _beamwinendt   = p.get<float>("BeamWindowEndT");
    _minflashpe    = p.get<float>("MinFlashPE");

    _minplength = p.get<float>("MinPLength");
    _minpscore  = p.get<float>("MinPScore");

    _fid_xmin   = p.get<float>("FiducialXMin");
    _fid_xmax   = p.get<float>("FiducialXMax");
    _fid_ymin   = p.get<float>("FiducialYMin");
    _fid_ymax   = p.get<float>("FiducialYMax");
    _fid_zmin   = p.get<float>("FiducialZMin");
    _fid_zmax   = p.get<float>("FiducialZMax");

    _minmuscore = p.get<float>("MinMuScore");
    _minpiscore = p.get<float>("MinPiScore");
                  
    _trackmodulelabel = p.get<std::string>("TrackModuleLabel");
    _flashmodulelabel = p.get<std::string>("FlashModuleLabel");
    _dtassoclabel     = p.get<std::string>("DTAssocLabel");
    _caloassoclabel   = p.get<std::string>("CaloAssocLabel");

    _lrcoefficients.clear();
    _lrcoefficients   = p.get< std::vector<float> >("LRCoefficients");

  }

  std::vector<float> LREvtReconstruction::ReconstructEvent(art::Event const& e)
  {

    float A = 17.;
    float bp1 = 0.58;
    float Mp = 0.938272;

    float Q2 = -9999.;
    float Tp = -9999.;
    float Rp = -9999.;
    float LR = -9999.;
    int TID = -1;
    float TH = -9999.;
    float PH = -9999.;
    float nux = -9999.;
    float nuy = -9999.;
    float nuz = -9999.;

    float bestLR = -9999.;
    float theta  = -9999.;
    float phi    = -9999.;

    // recover track handle
    art::Handle< std::vector<recob::Track> > trackVecHandle;
    e.getByLabel(_trackmodulelabel, trackVecHandle);
    // recover flash handle
    art::Handle< std::vector<recob::OpFlash> > flashVecHandle;
    e.getByLabel(_flashmodulelabel, flashVecHandle);

    // check for flash in beam window
    bool beamflash = false;
    if(flashVecHandle.isValid())
    {
      for(size_t flIdx = 0; flIdx < flashVecHandle->size(); flIdx++)
      {
        art::Ptr<recob::OpFlash> flash(flashVecHandle,flIdx);
        if(flash->Time() >= _beamwinstartt && flash->Time() <= _beamwinendt && flash->TotalPE() >= _minflashpe)
        {
          beamflash = true;
          break;
        }
      }
    }

    if(trackVecHandle.isValid() && beamflash)
    {
      // recover associations
      art::FindManyP<anab::Calorimetry> caloAssns(trackVecHandle, e, _caloassoclabel);
      art::FindManyP<anab::CosmicTag> dtAssns(trackVecHandle, e, _dtassoclabel); 

      // check for pi/mu bg
      float piflashdist = 999.;
      float muflashdist = 999.;
      for(size_t ctrkIdx = 0; ctrkIdx < trackVecHandle->size(); ctrkIdx++)
      {
        art::Ptr<recob::Track> ctrk(trackVecHandle,ctrkIdx);
        std::vector< art::Ptr<anab::CosmicTag> > cdtVec = dtAssns.at(ctrk.key());
        for(auto const& dttag : cdtVec)
        {
          if(dttag->CosmicType() == TAGID_PI && dttag->CosmicScore() > _minpiscore)
          {
            // found pi, is it near beam flash?
            if(MinFlashDistance(flashVecHandle,ctrk)[0] <= piflashdist)
            {
              // found pi bg
              piflashdist = MinFlashDistance(flashVecHandle,ctrk)[0];
            }
          }
          else if(dttag->CosmicType() == TAGID_MU && dttag->CosmicScore() > _minmuscore)
          {
            // found mu, is it near beam flash?
            if(MinFlashDistance(flashVecHandle,ctrk)[0] <= muflashdist)
            {
              // found mu bg
              muflashdist = MinFlashDistance(flashVecHandle,ctrk)[0];
            }
          }
        }
      }
      
      
      // Loop over input tracks
      for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
      {
        art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);
      
        // skip if too short or outside of fiducial
        if(trk->Length() < _minplength) continue;
        if(trk->Vertex().X() < _fid_xmin || trk->Vertex().X() > _fid_xmax
           || trk->Vertex().Y() < _fid_ymin || trk->Vertex().Y() > _fid_ymax
           || trk->Vertex().Z() < _fid_zmin || trk->Vertex().Z() > _fid_zmax) continue;
        if(trk->End().X() < _fid_xmin || trk->End().X() > _fid_xmax
            || trk->End().Y() < _fid_ymin || trk->End().Y() > _fid_ymax
            || trk->End().Z() < _fid_zmin || trk->End().Z() > _fid_zmax) continue;
      
        float predict_p   = -9999.;
        if(dtAssns.isValid())
        {
          std::vector< art::Ptr<anab::CosmicTag> > dtVec = dtAssns.at(trk.key());
          for(auto const& dttag : dtVec)
          {
            if(dttag->CosmicType() == TAGID_P)  predict_p   = dttag->CosmicScore();
          }
        }
        if(predict_p <= _minpscore) continue;
        // found proton candidate
        
        // check direction
        // get associated dedx
        bool tflip = false;
        int forward = 1;
        float tmpstartdedx;
        float tmpenddedx;
        std::vector<art::Ptr<anab::Calorimetry> > caloVec  = caloAssns.at(trk.key());
        if(!caloVec.empty())
        {
          // get collection plane
          std::vector<double> dedx;
          for(size_t iplane=0; iplane < caloVec.size();iplane++)
          {
            if(caloVec.at(iplane)->PlaneID().Plane == 2)
            {
              dedx = caloVec.at(iplane)->dEdx();
            }
          }
          // check if ends are nonzero for collection plane
          size_t endhits;
          if(dedx.size() >= _nendhits)
          {
            endhits = _nendhits;
          } else {
            endhits = TMath::Floor((float)dedx.size()/2.);
          }
          tmpenddedx   = 0.;
          tmpstartdedx = 0.;
          size_t nbgn = 0;
          size_t nend = 0;
          size_t idedx;
          idedx = 0;
          while(nbgn < endhits && idedx < dedx.size())
          {
            if(dedx.at(idedx) > 0.)
            {
              tmpstartdedx += dedx.at(idedx);
              nbgn += 1;
            }
            idedx += 1;
          }
          idedx = 0;
          while(nend < endhits && idedx < dedx.size())
          {
            if(dedx.at(idedx) > 0.)
            {
              tmpenddedx += dedx.at(dedx.size() - (1+idedx));
              nend += 1;
            }
            idedx += 1;
          }
          if(tmpenddedx < tmpstartdedx) tflip = true;
        }
        if(!tflip)
        {
          if(trk->Vertex().Z() > trk->End().Z()) forward = 0;
          theta  = trk->VertexDirection().Theta();
          phi    = trk->VertexDirection().Phi();
          nux    = trk->Vertex().X();
          nuy    = trk->Vertex().Y();
          nuz    = trk->Vertex().Z();
        }
        else
        {
          if(trk->End().Z() > trk->Vertex().Z()) forward = 0;
          theta  = TMath::Pi() - trk->EndDirection().Theta();
          float tmpphi = trk->EndDirection().Phi();
          phi    = (tmpphi > 0 ? tmpphi - TMath::Pi() :  tmpphi + TMath::Pi());
          nux    = trk->End().X();
          nuy    = trk->End().Y();
          nuz    = trk->End().Z();
        }
      
        // check if isolated
        float vtxdist = 999.;
        for(size_t vtrkIdx = 0; vtrkIdx < trackVecHandle->size(); vtrkIdx++)
        {
          if(vtrkIdx == trkIdx) continue;
          art::Ptr<recob::Track> vtrk(trackVecHandle,vtrkIdx);
      
          float ssdist = TMath::Sqrt(TMath::Power(trk->Vertex().X() - vtrk->Vertex().X(),2)
                                   + TMath::Power(trk->Vertex().Y() - vtrk->Vertex().Y(),2)
                                   + TMath::Power(trk->Vertex().Z() - vtrk->Vertex().Z(),2));
          float sedist = TMath::Sqrt(TMath::Power(trk->Vertex().X() - vtrk->End().X(),2)
                                   + TMath::Power(trk->Vertex().Y() - vtrk->End().Y(),2)
                                   + TMath::Power(trk->Vertex().Z() - vtrk->End().Z(),2));
          float sdist = TMath::Min(ssdist,sedist);
          float esdist = TMath::Sqrt(TMath::Power(trk->End().X() - vtrk->Vertex().X(),2)
                                   + TMath::Power(trk->End().Y() - vtrk->Vertex().Y(),2)
                                   + TMath::Power(trk->End().Z() - vtrk->Vertex().Z(),2));
          float eedist = TMath::Sqrt(TMath::Power(trk->End().X() - vtrk->End().X(),2)
                                   + TMath::Power(trk->End().Y() - vtrk->End().Y(),2)
                                   + TMath::Power(trk->End().Z() - vtrk->End().Z(),2));
          float edist = TMath::Min(esdist,eedist);
          float mdist = TMath::Min(edist,sdist);
      
          if(mdist < vtxdist)
          {
            vtxdist = mdist;
          }
        }

        std::vector<float> flashvec = MinFlashDistance(flashVecHandle,trk);
        
        float LRo = _lrcoefficients[0] + _lrcoefficients[1]*predict_p 
                  + _lrcoefficients[2]*vtxdist + _lrcoefficients[3]*flashvec[0] 
                  + _lrcoefficients[4]*flashvec[1] + _lrcoefficients[5]*forward 
                  + _lrcoefficients[6]*muflashdist + _lrcoefficients[7]*piflashdist;
        LR = TMath::Exp(LRo)/(1. + TMath::Exp(LRo));

        std::cout << "run: " << e.id().run() << ", event: " << e.id().event() << std::endl;
        std::cout << "predict_p: " << predict_p << ", vtxdist: " << vtxdist << ", zflash: " << flashvec[0] << ", yflash: " << flashvec[1] << std::endl;
        std::cout << "forward: " << forward << ", muflashdist: " << muflashdist << ", piflashdist: " << piflashdist << std::endl;
        std::cout << "LRscore: " << LR << std::endl;
        
        // found candidate
        if(LR > bestLR)
        {
          bestLR = LR;
          Rp = trk->Length();
          Tp = A/bp1* TMath::Power(Rp,bp1) /1000.;
          Q2 = 2.*Tp*Mp;
          TID = trk->ID();
          TH = theta;
          PH = phi;
        }
        
      }
      
    }

    std::vector<float> recovec = {Q2,Tp,Rp,bestLR,(float)TID,TH,PH,nux,nuy,nuz};
    return recovec;
  }

  std::vector<float> LREvtReconstruction::MinFlashDistance(
            art::Handle< std::vector<recob::OpFlash> > flashVecHandle,art::Ptr<recob::Track> trk)
  {
    
    float minz = 9999.;
    float miny = 9999.;
    for(size_t flIdx = 0; flIdx < flashVecHandle->size(); flIdx++)
    {
      art::Ptr<recob::OpFlash> flash(flashVecHandle,flIdx);
      if(flash->Time() >= _beamwinstartt && flash->Time() <= _beamwinendt 
          && flash->TotalPE() >= _minflashpe)
      {
        float trkzcenter = (trk->Vertex().Z() + trk->End().Z())/2.;
        if(TMath::Abs(trkzcenter - flash->ZCenter()) < minz)
        {
          minz = TMath::Abs(trkzcenter - flash->ZCenter());
          float trkycenter = (trk->Vertex().Y() + trk->End().Y())/2.;
          miny = TMath::Abs(trkycenter - flash->YCenter());
        }
      }
    }
    std::vector<float> flsvec = {minz,miny};
    return flsvec;
  }
}

#endif

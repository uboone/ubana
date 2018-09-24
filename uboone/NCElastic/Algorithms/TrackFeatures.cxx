#ifndef TRACKFEATURES_CXX
#define TRACKFEATURES_CXX

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
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TMath.h"

#include <iostream>
#include <map>
#include <cmath>

#include "uboone/NCElastic/Algorithms/TrackFeatures.h"

namespace dtfeatures {

  TrackFeatures::TrackFeatures()
  {
    _nendhits      = 6;
    _beamwinstartt = 3.2;
    _beamwinendt   = 4.8;

  }

  void TrackFeatures::Configure(fhicl::ParameterSet const& p)
  {
    _nendhits      = p.get<unsigned int>("NEndHits");
    _beamwinstartt = p.get<float>("BeamWindowStartT");
    _beamwinendt   = p.get<float>("BeamWindowEndT");

    _pidassoclabel        = p.get<std::string>("PIDAssocLabel");
    _caloassoclabel       = p.get<std::string>("CaloAssocLabel");
    _cosmictagassoclabel  = p.get<std::string>("CosmicTagAssocLabel");
    _containtagassoclabel = p.get<std::string>("ContainTagAssocLabel");
  }

  std::vector<std::vector<float>> TrackFeatures::CreateFeatures(art::Event const& e, 
                                       art::Handle<std::vector<recob::Track>>& trackVecHandle )
  {
    // recover associations
    art::FindManyP<anab::ParticleID> pidAssns(trackVecHandle, e, _pidassoclabel);

    art::FindManyP<anab::Calorimetry> caloAssns(trackVecHandle, e, _caloassoclabel);
    art::FindManyP<anab::CosmicTag> cosmictagAssns(trackVecHandle, e, _cosmictagassoclabel); 
    art::FindManyP<anab::CosmicTag> containtagAssns(trackVecHandle, e, _containtagassoclabel);

    std::vector< std::vector<float> > evtdata;
    if(trackVecHandle.isValid())
    {
      // Loop over tracks in event and get features
      for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
      {
        art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);

        std::vector<float> trkdata;
        
        std::vector<art::Ptr<anab::ParticleID> > pidVec    = pidAssns.at(trk.key());
        std::vector<art::Ptr<anab::Calorimetry> > caloVec  = caloAssns.at(trk.key());
        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec  = cosmictagAssns.at(trk.key());
        std::vector<art::Ptr<anab::CosmicTag> > containVec = containtagAssns.at(trk.key());
        
        int   nhits         = 0;
        float length        = -9999.;
        float starty        = -9999.;
        float startz        = -9999.;
        float endy          = -9999.;
        float endz          = -9999.;
        float theta         = -9999.;
        float phi           = -9999.;
        float distlenratio  = -9999.;
        float startdedx     = -9999.;
        float enddedx       = -9999.;
        float dedxratio     = -9999.;
        float totaldedx     = -9999.;
        float trtotaldedx     = -9999.;
        float averagededx   = -9999.;
        float traveragededx = -9999.;
        float cosmicscore   = -9999.;
        float coscontscore  = -9999.;
       
        float tmpphi;
        float tmpstartdedx;
        float tmpenddedx;
        float startx;
        float endx;
        float dist;
        
        bool tflip          = false;

        size_t bstplane = 2;
        if(!caloVec.empty())
        {
          
          std::vector<double> dedx;
          totaldedx = 0.;
          trtotaldedx = 0.;
          for(size_t iplane = 0; iplane < caloVec.size(); iplane++)
          {
            if(caloVec.at(iplane)->PlaneID().Plane == bstplane)
            {
              dedx = caloVec.at(iplane)->dEdx();
            }
          }
          nhits = dedx.size();

          if(dedx.size() > 0)
          {
            double totaldedx2 = 0.;
            for( auto& de : dedx)
            {
              totaldedx += de;
              totaldedx2 += TMath::Power(de,2);
            }
            averagededx = (float)totaldedx / (float)dedx.size();
            std::vector<double> dedx_sort = dedx;
            std::sort(dedx_sort.begin(),dedx_sort.end());

            double dedxmed = dedx_sort[std::floor(dedx.size()/2.)];
            double dedxstd = TMath::Sqrt((float)totaldedx2 / (float)dedx.size() - TMath::Power(averagededx,2));
            int ntrdedx = 0;
            for( auto& de : dedx)
            {
              if(TMath::Abs(de - dedxmed) <= dedxstd)
              {
                trtotaldedx += de;
                ntrdedx += 1;
              }
            }
            traveragededx = trtotaldedx/(float)ntrdedx;

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

            if(tmpenddedx < tmpstartdedx)
            {
              tflip     = true;
              enddedx   = tmpstartdedx;
              startdedx = tmpenddedx;
            } else {
              tflip     = false;
              enddedx   = tmpenddedx;
              startdedx = tmpstartdedx;
            }
            dedxratio = enddedx/startdedx;
          }
          else
          {
            averagededx = 0.;
            traveragededx = 0.;
            startdedx = 0.;
            enddedx = 0.;
            dedxratio = 0.;
          }
        }
        
        length = trk->Length();
        if(!tflip)
        {
          startx = trk->Vertex().X();
          starty = trk->Vertex().Y();
          startz = trk->Vertex().Z();
          endx   = trk->End().X();
          endy   = trk->End().Y();
          endz   = trk->End().Z();
          theta  = trk->VertexDirection().Theta();
          phi    = trk->VertexDirection().Phi();
        } else {
          startx = trk->End().X();
          starty = trk->End().Y();
          startz = trk->End().Z();
          endx   = trk->Start().X();
          endy   = trk->Start().Y();
          endz   = trk->Start().Z();
          theta  = TMath::Pi() - trk->EndDirection().Theta();
          tmpphi = trk->EndDirection().Phi();
          phi    = (tmpphi > 0 ? tmpphi - TMath::Pi() :  tmpphi + TMath::Pi());
        }
        
        dist = TMath::Sqrt( TMath::Power(endx-startx,2) 
                          + TMath::Power(endy-starty,2) 
                          + TMath::Power(endz-startz,2));
        distlenratio = dist/length;
        
        if(!cosmicVec.empty()) cosmicscore = cosmicVec.at(0)->CosmicScore();
        if(!containVec.empty()) coscontscore = containVec.at(0)->CosmicScore();
        
        trkdata.clear();
        
        trkdata.emplace_back(nhits);
        trkdata.emplace_back(length);
        trkdata.emplace_back(starty);
        trkdata.emplace_back(startz);
        trkdata.emplace_back(endy);
        trkdata.emplace_back(endz);
        trkdata.emplace_back(theta);
        trkdata.emplace_back(phi);
        trkdata.emplace_back(distlenratio);
        trkdata.emplace_back(startdedx);
        trkdata.emplace_back(dedxratio);
        trkdata.emplace_back(trtotaldedx);
        trkdata.emplace_back(traveragededx);
        trkdata.emplace_back(cosmicscore);
        trkdata.emplace_back(coscontscore);

        evtdata.emplace_back(trkdata);
      }
    }
    
    return evtdata; 
  }

  std::vector<anab::CosmicTagID_t> TrackFeatures::CosmicTagLabel()
  {
    std::vector<anab::CosmicTagID_t> ctlabelvec;

    const anab::CosmicTagID_t TAGID_P  = anab::CosmicTagID_t::kGeometry_YY;
    const anab::CosmicTagID_t TAGID_MU = anab::CosmicTagID_t::kGeometry_YZ;
    const anab::CosmicTagID_t TAGID_PI = anab::CosmicTagID_t::kGeometry_ZZ;
    const anab::CosmicTagID_t TAGID_EM = anab::CosmicTagID_t::kGeometry_XX;
    const anab::CosmicTagID_t TAGID_CS = anab::CosmicTagID_t::kGeometry_XY;

    ctlabelvec.clear();
    ctlabelvec.emplace_back(TAGID_P);    
    ctlabelvec.emplace_back(TAGID_MU);    
    ctlabelvec.emplace_back(TAGID_PI);    
    ctlabelvec.emplace_back(TAGID_EM);    
    ctlabelvec.emplace_back(TAGID_CS);    

    return ctlabelvec;
  }


}

#endif

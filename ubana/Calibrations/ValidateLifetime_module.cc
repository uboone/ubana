////////////////////////////////////////////////////////////////////////
// Class:       ValidateLifetime
// Module Type: analyzer
// File:        ValidateLifetime_module.cc
//
// Generated at Sun Oct  9 22:21:11 2016 by Tingjun Yang using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TH2D.h"

namespace ub {
  class ValidateLifetime;
}

class ub::ValidateLifetime : public art::EDAnalyzer {
public:
  explicit ValidateLifetime(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ValidateLifetime(ValidateLifetime const &) = delete;
  ValidateLifetime(ValidateLifetime &&) = delete;
  ValidateLifetime & operator = (ValidateLifetime const &) = delete;
  ValidateLifetime & operator = (ValidateLifetime &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  calo::CalorimetryAlg caloAlg;

  TH2D *dqdxtime[3];
  TH2D *dqdxtimecor[3];
  TH1D *hdqdx[3];
  TH1D *hdqdxcor[3];
  TH1D *hdedx[3];

  TH1D *hitdis;
  TH1D *hz;

  art::ServiceHandle<geo::Geometry> geom;
  detinfo::DetectorProperties const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double XDriftVelocity = detprop->DriftVelocity()*1e-3; //cm/ns

};


ub::ValidateLifetime::ValidateLifetime(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  ,
  fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
  fCalorimetryModuleLabel(p.get<std::string>("CalorimetryModuleLabel")),
  caloAlg(p.get< fhicl::ParameterSet >("CaloAlg"))
{}

void ub::ValidateLifetime::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track>> tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);

  for (size_t i = 0; i<tracklist.size(); ++i){
    const auto& trkstart = tracklist[i]->Vertex();
    const auto& trkend = tracklist[i]->End();
    double X = std::abs(trkstart.X()-trkend.X());
    if (X>250 && X<270){
      if (fmcal.isValid()){
        auto &calos = fmcal.at(i);
        for (size_t j = 0; j<calos.size(); ++j){
          if (!calos[j]) continue;
          if (!calos[j]->PlaneID().isValid) continue;
          int planenum = calos[j]->PlaneID().Plane;
          if (planenum<0||planenum>2) continue;
          const auto& dir_start = tracklist[i]->VertexDirection();
          const geo::WireGeo& wire = geom->TPC().Plane(planenum).MiddleWire();
          double wirestart[3], wireend[3];
          wire.GetStart(wirestart);
          wire.GetEnd(wireend);
          double cosangle = (dir_start.Y()*(wirestart[1]-wireend[1])+dir_start.Z()*(wirestart[2]-wireend[2]))/sqrt((pow(dir_start.Y(),2)+pow(dir_start.Z(),2))*(pow(wirestart[1]-wireend[1],2)+pow(wirestart[2]-wireend[2],2)));
          if (std::abs(cosangle)<0.5){
            double minx = 1e10;
            const size_t NHits = calos[j] -> dEdx().size();
            for(size_t iHit = 0; iHit < NHits; ++iHit) {
              if ((calos[j]->TrkPitchVec())[iHit]>1) continue;
              const auto& TrkPos = (calos[j] -> XYZ())[iHit];
              if (TrkPos.X()<minx)
                minx = TrkPos.X();
            }// loop NHits
            for(size_t iHit = 0; iHit < NHits; ++iHit) {
              if ((calos[j]->TrkPitchVec())[iHit]>1) continue;
              const auto& TrkPos1 = (calos[j] -> XYZ())[iHit];
              double x = TrkPos1.X()-minx; //subtract the minx to get correct t0
              double t = x/(XDriftVelocity*1000); //change the velocity units to cm/ns to cm/us
              dqdxtime[planenum]->Fill(t, (calos[j] -> dQdx())[iHit]);
              double dqdxcor = (calos[j] -> dQdx())[iHit]*caloAlg.LifetimeCorrection(t/(detprop->SamplingRate()*1.e-3)+detprop->TriggerOffset());
              dqdxtimecor[planenum]->Fill(t, dqdxcor);
              hdqdx[planenum]->Fill((calos[j] -> dQdx())[iHit]);
              hdqdxcor[planenum]->Fill(dqdxcor);
              double dedx = caloAlg.dEdx_AREA((calos[j] -> dQdx())[iHit], t/(detprop->SamplingRate()*1.e-3)+detprop->TriggerOffset(), planenum);
              hdedx[planenum]->Fill(dedx); 
              //if (planenum==2&&(calos[j] -> dQdx())[iHit]<50){
              //if (planenum==2&&dqdxcor<140){
              if (planenum==2&&dedx<1.3){
                std::cout<<"Low dQ/dx "<<evt.run()<<" "<<evt.subRun()<<" "<<evt.id().event()<<" "<<i<<" "<<TrkPos1.X()<<" "<<TrkPos1.Y()<<" "<<TrkPos1.Z()<<" "<<(calos[j] -> dQdx())[iHit]<<" "<<dqdxcor<<std::endl;
                hz->Fill(TrkPos1.Z());
              }
//              if (planenum==2){
//                std::cout<<(calos[j] -> dQdx())[iHit]<<" "<<(calos[j]->TrkPitchVec())[iHit]<<std::endl;
//              }
              if (iHit<NHits-1){
                hitdis->Fill((calos[j]->XYZ()[iHit]-calos[j]->XYZ()[iHit+1]).R());
              }
            }
          }//track not overlaping with wire
        }//all calorimetry objects
      }//found calorimetry
    }//crossing track
  }//all tracks
}
            
void ub::ValidateLifetime::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  for (int i = 0; i<3; ++i){
    dqdxtime[i] = tfs->make<TH2D>(Form("dqdxtime_%d",i), Form("Plane %d;Drift time (#mus);dQ/dx (ADC/cm)",i),100,0,2200,100,0,600);
    dqdxtimecor[i] = tfs->make<TH2D>(Form("dqdxtimecor_%d",i), Form("Plane %d, lifetime corrected;Drift time (#mus);dQ/dx (ADC/cm)",i),100,0,2200,100,0,600);
    hdqdx[i] = tfs->make<TH1D>(Form("hdqdx_%d",i), Form("Plane %d; dQ/dx (ADC/cm); Events",i), 100,0,600);
    hdqdx[i]->Sumw2();
    hdqdxcor[i] = tfs->make<TH1D>(Form("hdqdxcor_%d",i), Form("Plane %d, lifetime corrected; dQ/dx (ADC/cm); Events",i), 100,0,600);
    hdqdxcor[i]->Sumw2();
    hdedx[i] = tfs->make<TH1D>(Form("hdedx_%d",i), Form("Plane %d; dE/dx (MeV/cm); Events",i), 100,0,10);
    hdedx[i]->Sumw2();
  }
  hitdis = tfs->make<TH1D>("hitdis","hit distance",100,0,2);
  hitdis->Sumw2();
  hz = tfs->make<TH1D>("hz","hz",2000,0,1000);
}

DEFINE_ART_MODULE(ub::ValidateLifetime)

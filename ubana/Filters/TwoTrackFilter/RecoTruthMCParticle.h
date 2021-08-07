#ifndef RECOTRUTHMCPARTICLE_H
#define RECOTRUTHMCPARTICLE_H

#include <iostream>
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

class RecoTruth{
  public:
    //RecoTruth();

    void TrackToMCParticle(art::Event const & evt, const art::Handle<std::vector<recob::Hit> >& hit_handle, std::string track_producer);

}; 

#endif

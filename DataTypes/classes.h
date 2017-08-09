
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

#include "uboone/UBXSec/DataTypes/FlashMatch.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include <vector>

template class art::Assns<anab::FlashMatch,recob::PFParticle>;
template class art::Assns<recob::PFParticle,anab::FlashMatch>;

template class art::Wrapper<art::Assns<anab::FlashMatch,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, anab::FlashMatch >>;


template class std::vector<ubana::FlashMatch>;

template class art::Assns<ubana::FlashMatch,recob::PFParticle>;
template class art::Assns<recob::PFParticle,ubana::FlashMatch>;
template class art::Assns<ubana::FlashMatch,recob::Track,void>;
template class art::Assns<recob::Track,ubana::FlashMatch,void>;

template class art::Wrapper<art::Assns<ubana::FlashMatch,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, ubana::FlashMatch >>;
template class art::Wrapper<art::Assns<ubana::FlashMatch,recob::Track,void> >;
template class art::Wrapper<art::Assns<recob::Track,ubana::FlashMatch,void> >;




template class std::vector<ubana::TPCObject>;

template class art::Assns<ubana::TPCObject,recob::PFParticle>;
template class art::Assns<recob::PFParticle,ubana::TPCObject>;
template class art::Assns<ubana::TPCObject,recob::Track,void>;
template class art::Assns<recob::Track,ubana::TPCObject,void>;
template class art::Assns<ubana::TPCObject,ubana::FlashMatch,void>;
template class art::Assns<ubana::FlashMatch,ubana::TPCObject,void>;

template class art::Wrapper<art::Assns<ubana::TPCObject,  recob::PFParticle>>;
template class art::Wrapper<art::Assns<recob::PFParticle, ubana::TPCObject >>;
template class art::Wrapper<art::Assns<recob::Track,ubana::TPCObject,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,recob::Track,void> >;
template class art::Wrapper<art::Assns<ubana::TPCObject,ubana::FlashMatch,void> >;
template class art::Wrapper<art::Assns<ubana::FlashMatch,ubana::TPCObject,void> >;


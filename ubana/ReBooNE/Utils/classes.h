//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by wes
//

#include "canvas/Persistency/Common/Wrapper.h"

// data-products
// lardataobj
//#include "lardata/Utilities/AssociationUtil.h"
#include "canvas/Persistency/Common/Assns.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "ubana/ReBooNE/Utils/BooNEInfo.h"

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class std::vector<BooNEInfo>;
template class art::Wrapper<std::vector<BooNEInfo> >;

template class art::Assns<BooNEInfo, simb::MCTruth>;
template class art::Wrapper<art::Assns<BooNEInfo, simb::MCTruth> >;

template class art::Assns<simb::MCTruth, BooNEInfo>;
template class art::Wrapper<art::Assns<simb::MCTruth, BooNEInfo> >;



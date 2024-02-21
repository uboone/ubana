#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_ONLYBLIPS_H
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_ONLYBLIPS_H

// blip reco includes
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "variables.h"


namespace BVA_ana
{
	void AnalyzeBlips(blip::BlipRecoAlg* fBlipAlg, var_all& vars);


}









#endif

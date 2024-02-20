#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_INIT_BRANCHES_H
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_INIT_BRANCHES_H

#include "variables.h"

namespace BVA_ana
{
	void ClearVars(var_all& vars);
	void CreateVarsBranches(var_all& vars);

	void ClearVarsEvt(var_evt& vars);
	void CreateVarsEvtBranches(var_evt& vars);
}


#endif

#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_H
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_H

#include "TVectorD.h"
#include "variables.h"

namespace BVA_ana
{
	//Look at blips near the vertex;
	void AnalyzeBlipsNearVertex(var_all& vars);

	//subfunctions
	void SummarizeDist(var_all& vars, unordered_map<int, vector<int>>groups_dist);
	void SummarizeAngInRad(var_all& vars, TVector3 vertex, unordered_map<int, vector<int>>groups_dist);
}


#endif

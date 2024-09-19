#ifndef UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_H
#define UBANA_BLIPVERTEXANALYSIS_LIBRARIES_BLIPSNEARVERTEX_H

#include "TVectorD.h"
#include "variables.h"
#include "unordered_map"

namespace BVA_ana
{
	//Look at blips near the vertex;
	void AnalyzeBlipsNearVertex(var_all& vars);

	//subfunctions
	void SummarizeDist(var_all& vars, std::unordered_map<int, std::vector<int>>groups_dist);
	double GetThetaYZDegree(TVector3 vertex, TVector3 blip);
	double GetPhiXZDegree(TVector3 vertex, TVector3 blip);
///	double GetPhiXZDegree(TVector3 vertex, TVector3 blip, TVector3 beam);
//	void SummarizeAngInRad(var_all& vars, TVector3 vertex, std::vector<double>groups_dist);
}


#endif

#ifndef _FittedVertex_h_
#define _FittedVertex_h_

#include "TVector3.h"

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class FittedVertex : public TObject{
#else
class FittedVertex {
#endif

public:

FittedVertex() {}
~FittedVertex() {}

//indices of two tracks
int Index_1 = -1;
int Index_2 = -1;

double X;
double Y;
double Z;

double Chi2_1 = -1;
double Chi2_2 = -1;

#ifdef __MAKE_ROOT_DICT__
ClassDef(FittedVertex,1);
#endif


};

#endif

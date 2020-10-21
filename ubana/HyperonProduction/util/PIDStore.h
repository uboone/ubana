#ifndef _PIDStore_h_
#define _PIDStore_h_

#include <iostream>
//#include "TString.h"

#ifdef __MAKE_ROOT_DICT__
#include "TObject.h"
#endif

#ifdef __MAKE_ROOT_DICT__
class PIDStore : public TObject{
#else
class PIDStore {
#endif

public:

PIDStore() {}
~PIDStore() {}

//int planeID=-1; //plane no, -1 if not applicabale
//TString name;
int AssumedPDG;

//int true_pdg;
//double true_length;



#ifdef __MAKE_ROOT_DICT__
ClassDef(PIDStore,1);
#endif


};



#endif

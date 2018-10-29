// Header file including the function Varuna has been using to convert dQdx to dEdx
// Should become obsolete when we have the new dEdx calibration

#ifndef VARUNA_GETDEDX_H
#define VARUNA_GETDEDX_H

double _Rho = 1.383;
double _betap = 0.212;
double _alpha = 0.93;
double _Wion = 23.6e-6;
double _Efield = 0.273;

////////////////////////////////// Function definition ////////////////////////////////

inline double VarunaGetdEdx(double dqdx){
 return (exp(dqdx*(_betap/(_Rho*_Efield)*_Wion))-_alpha)/(_betap/(_Rho*_Efield));
}

//////////////////////////////////////////////////////////////////////////////////////



#endif

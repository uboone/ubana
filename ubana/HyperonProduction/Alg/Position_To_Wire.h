#ifndef _Position_To_Wire_h_
#define _Position_To_Wire_h_

#include "TVector3.h"

// Convert 3D positions to channels/ticks
// See docdb 25505-v2 for explanation

const double A_w = 3.33328;
const double C_U = 338.140;
const double C_V = 2732.53;
const double C_Y = 4799.19;
const double A_t = 18.2148;
const double C_t = 818.351;

// cos(60) and sin(60)
const double cos60 = 0.5;
const double sin60 = sqrt(3)/2.0;

int U_wire(TVector3 pos) { return A_w*(-sin60*pos.Y()+cos60*pos.Z())+C_U; }
int V_wire(TVector3 pos) { return A_w*(sin60*pos.Y()+cos60*pos.Z())+C_V; }
int Y_wire(TVector3 pos) { return A_w*pos.Z() + C_Y; }
int tick(TVector3 pos) { return A_t*pos.X() + C_t; }

#endif

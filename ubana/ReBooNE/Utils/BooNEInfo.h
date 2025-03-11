#ifndef BOONEINFO_H
#define BOONEINFO_H

//boone ntuple structure
class BooNEInfo {  
  
  public:

    BooNEInfo(){};
    ~BooNEInfo(){};

    //Neutrino information
    int pdg;

    double nu_vtx_x; // neutrino x position at flux window
    double nu_vtx_y; // neutrino y position at flux window
    double nu_vtx_z; // neutrino z position at flux window (z position of the flux window)
    double nu_vtx_t; // neutrino time at flux window

    double nu_startt; // neutrino start time 
    double nu_startx; // neutrino start time
    double nu_starty; // neutrino start time
    double nu_startz; // neutrino start time

    //Beam information
    float tank_pos_beam[3]; /// 3-position of the flux window;
    float targ_pos_beam[3]; /// Frame conversion from beam to flux frame
    float windowbase[3];

};

#endif

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

class Topology {

	public:

		// Default constructor
		Topology(){}

		// Default destructor
		~Topology(){}

                int TopologyLabel(int Nmuons, int Nelectrons, int Npiplus_abTH, int Npiplus_blTH, int Npiminus_abTH, int Npiminus_blTH, int Npi0, int Nprotons_abTH, int Nprotons_blTH, int PDG, int CCNC, bool ifbeam, bool ifFV);
                //int TopologyLabel(int Nmuons, int Nelectrons, int Npiplus_abTH, int Npiplus_blTH, int Npiminus_abTH, int Npiminus_blTH, int Npi0, int Nprotons_abTH, int Nprotons_blTH, int PDG = -99, int CCNC = -99, bool ifbeam, bool ifFV);


	private:


};

#endif

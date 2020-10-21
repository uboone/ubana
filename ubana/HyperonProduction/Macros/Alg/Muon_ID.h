#ifndef _Muon_ID_h_
#define _Muon_ID_h_

//returns index of muon candidate
//reuturn -1 if no muon found

int Muon_ID(std::vector<double> lengths , std::vector<double> PIDs , double pid_cut = 0.6){


int i_longest=-1;
int length=-1;

//find longest mip like track

for(size_t i=0;i<lengths.size();i++){


//if pid says track is too proton like, skip
if( PIDs.at(i) > pid_cut ) continue;

if(lengths.at(i) > length){

i_longest = i;
length = lengths.at(i);

}

}


return i_longest;

}



#endif


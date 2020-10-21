#ifndef _Track_Lengths_h_
#define _Track_Lengths_h_

bool Track_Length_Cut(std::vector<double> lengths, int i_muon , double cut_secondary=65 , double cut_tertiary=35){

//take muon candidate out of length vector
lengths.erase(lengths.begin()+i_muon);

//place remaining tracks in order of length
std::sort(lengths.begin(),lengths.end());

if(lengths.at(lengths.size()-1) > cut_secondary) return false;
if(lengths.at(lengths.size()-2) > cut_tertiary) return false;

return true;

}


#endif

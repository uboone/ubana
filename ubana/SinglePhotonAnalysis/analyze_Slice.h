#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearSlices(){
        m_reco_num_slices = 0;
        m_reco_slice.clear();
    }

    void SinglePhoton::ResizeSlices(size_t size){
        m_reco_slice.resize(size);
    }


    void SinglePhoton::CreateSliceBranches(){
        vertex_tree->Branch("reco_slice",&m_reco_slice);
    }

    void SinglePhoton::AnalyzeSlices(){
        m_reco_num_slices = 1;
        this->ResizeSlices(m_reco_num_slices);

        m_reco_slice[0]=99;

    }
}

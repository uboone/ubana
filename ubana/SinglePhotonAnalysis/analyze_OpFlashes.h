#include "SinglePhoton_module.h"

namespace single_photon
{
    void SinglePhoton::ClearFlashes(){
        m_reco_num_flashes =0;
        m_reco_num_flashes_in_beamgate =0;
        m_reco_flash_total_pe.clear();
        m_reco_flash_time.clear();
        m_reco_flash_time_width.clear();
        m_reco_flash_abs_time.clear();
        m_reco_flash_frame.clear();
        m_reco_flash_ycenter.clear();
        m_reco_flash_ywidth.clear();
        m_reco_flash_zcenter.clear();
        m_reco_flash_zwidth.clear();
        m_reco_flash_total_pe_in_beamgate.clear();
        m_reco_flash_time_in_beamgate.clear();
        m_reco_flash_ycenter_in_beamgate.clear();
        m_reco_flash_zcenter_in_beamgate.clear();
        m_CRT_min_hit_time = -999;
        m_CRT_min_hit_PE = -999;
        m_CRT_min_hit_x = -999;
        m_CRT_min_hit_y = -999;
        m_CRT_min_hit_z = -999;
        m_CRT_hits_time.clear();
        m_CRT_hits_PE.clear();
        m_CRT_hits_x.clear(); 
        m_CRT_hits_y.clear();
        m_CRT_hits_z.clear();
        m_CRT_dt = -999;

    }

    void SinglePhoton::ResizeFlashes(size_t size){
        m_reco_flash_total_pe.resize(size);
        m_reco_flash_time.resize(size);
        m_reco_flash_time_width.resize(size);
        m_reco_flash_abs_time.resize(size);
        m_reco_flash_frame.resize(size);
        m_reco_flash_ycenter.resize(size);
        m_reco_flash_ywidth.resize(size);
        m_reco_flash_zcenter.resize(size);
        m_reco_flash_zwidth.resize(size);
        m_reco_flash_total_pe_in_beamgate.resize(size);
        m_reco_flash_time_in_beamgate.resize(size);
        m_reco_flash_ycenter_in_beamgate.resize(size);
        m_reco_flash_zcenter_in_beamgate.resize(size);
        m_CRT_hits_time.resize(size);
        m_CRT_hits_PE.resize(size);
        m_CRT_hits_x.resize(size); 
        m_CRT_hits_y.resize(size);
        m_CRT_hits_z.resize(size);

    }


    void SinglePhoton::CreateFlashBranches(){
        vertex_tree->Branch("beamgate_flash_start",&m_beamgate_flash_start,"beamgate_flash_start/D");
        vertex_tree->Branch("beamgate_flash_end",&m_beamgate_flash_end,"beamgate_flash_end/D");
        vertex_tree->Branch("reco_num_flashes",&m_reco_num_flashes,"reco_num_flashes/I");
        vertex_tree->Branch("reco_num_flashes_in_beamgate",&m_reco_num_flashes_in_beamgate,"reco_num_flashes_in_beamgate/I");
        vertex_tree->Branch("reco_flash_total_pe", &m_reco_flash_total_pe);
        vertex_tree->Branch("reco_flash_time", &m_reco_flash_time);
        vertex_tree->Branch("reco_flash_time_width",&m_reco_flash_time_width);
        vertex_tree->Branch("reco_flash_abs_time",&m_reco_flash_abs_time);
        vertex_tree->Branch("reco_flash_frame",&m_reco_flash_frame);
        vertex_tree->Branch("reco_flash_ycenter",&m_reco_flash_ycenter);
        vertex_tree->Branch("reco_flash_ywidth",&m_reco_flash_ywidth);
        vertex_tree->Branch("reco_flash_zcenter",&m_reco_flash_zcenter);
        vertex_tree->Branch("reco_flash_zwidth",&m_reco_flash_zwidth);
        vertex_tree->Branch("reco_flash_total_pe_in_beamgate", &m_reco_flash_total_pe_in_beamgate);
        vertex_tree->Branch("reco_flash_time_in_beamgate", &m_reco_flash_time_in_beamgate);
        vertex_tree->Branch("reco_flash_ycenter_in_beamgate",&m_reco_flash_ycenter_in_beamgate);
        vertex_tree->Branch("reco_flash_zcenter_in_beamgate",&m_reco_flash_zcenter_in_beamgate);

        vertex_tree->Branch("CRT_min_hit_time",&m_CRT_min_hit_time,"CRT_min_hit_time/D");
        vertex_tree->Branch("CRT_min_hit_PE",&m_CRT_min_hit_PE,"CRT_min_hit_PE/D");
        vertex_tree->Branch("CRT_min_hit_x",&m_CRT_min_hit_x,"CRT_min_hit_x/D");
        vertex_tree->Branch("CRT_min_hit_y",&m_CRT_min_hit_y,"CRT_min_hit_y/D");
        vertex_tree->Branch("CRT_min_hit_z",&m_CRT_min_hit_z,"CRT_min_hit_z/D");

        vertex_tree->Branch("CRT_hits_time",&m_CRT_hits_time);
        vertex_tree->Branch("CRT_hits_PE",&m_CRT_hits_PE);
        vertex_tree->Branch("CRT_hits_x",&m_CRT_hits_x);
        vertex_tree->Branch("CRT_hits_y",&m_CRT_hits_y);
        vertex_tree->Branch("CRT_hits_z",&m_CRT_hits_z);

    }


    void SinglePhoton::AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes){

        size_t flash_size = flashes.size();
        m_reco_num_flashes = flash_size;

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeFlashes()\t||\t Beginning analysis of recob::OpFlash\n";

        this->ResizeFlashes(flash_size);

        for(size_t i = 0; i < flash_size; ++i) {


            art::Ptr<recob::OpFlash> const & flash = flashes[i];

            m_reco_flash_total_pe[i]=(flash->TotalPE());
            m_reco_flash_time[i]=(flash->Time());
            m_reco_flash_time_width[i]=flash->TimeWidth();
            m_reco_flash_abs_time[i]=flash->AbsTime();
            m_reco_flash_frame[i]=flash->Frame();
            m_reco_flash_ycenter[i]=flash->YCenter();
            m_reco_flash_ywidth[i]=flash->YWidth();
            m_reco_flash_zcenter[i]=flash->ZCenter();
            m_reco_flash_zwidth[i]=flash->ZWidth();

            if(m_reco_flash_time[i] <= m_beamgate_flash_end && m_reco_flash_time[i] >= m_beamgate_flash_start){
                m_reco_num_flashes_in_beamgate++;
                m_reco_flash_total_pe_in_beamgate[i]=(flash->TotalPE());
                m_reco_flash_time_in_beamgate[i]=(flash->Time());
                m_reco_flash_ycenter_in_beamgate[i] = flash->YCenter();
                m_reco_flash_zcenter_in_beamgate[i] = flash->ZCenter();
            }



        }

        if(m_is_verbose) std::cout<<"SinglePhoton::AnalyzeFlashes()\t||\t Finished. There was "<<flash_size<<" flashes with: "<<m_reco_num_flashes_in_beamgate<<" in the beamgate defined by: "<<m_beamgate_flash_start<<" <-> "<<m_beamgate_flash_end<<std::endl;

    }

    //fill these values only for events that have CRT information - run3 G and later
    //code taken from ubcrt/UBCRTCosmicFilter/UBCRTCosmicFilter_module.cc 


}

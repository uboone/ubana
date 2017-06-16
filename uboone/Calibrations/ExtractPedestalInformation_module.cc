///////////////////////////////////////////////////////////////////////
// Class:       ExtractPedestalInformation
// Module Type: analyzer
// File:        ExtractPedestalInformation_module.cc
//
// Author:      Adam Lister
// Email:       a.lister1@lancaster.ac.uk
//
// Pulls out pedestal information on a channel-by-channel basis in order
// to fill databases
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes 
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

// C++ includes
#include <iostream>
#include <fstream>
#include <vector>

//class ExtractPedestalInformation;

class ExtractPedestalInformation : public art::EDAnalyzer {
    public:
        explicit ExtractPedestalInformation(fhicl::ParameterSet const& p);
        ExtractPedestalInformation(ExtractPedestalInformation const &) = delete;
        ExtractPedestalInformation(ExtractPedestalInformation &&) = delete;
        ExtractPedestalInformation & operator = (ExtractPedestalInformation const &) = delete;
        ExtractPedestalInformation & operator = (ExtractPedestalInformation &&) = delete;

        // Required functions.
        virtual void beginJob() override;
        virtual void reconfigure(fhicl::ParameterSet const& parameterSet) override;
        virtual void analyze(const art::Event& e) override;
        virtual void endJob() override;

    private:

        // Declare member data here.

        // parameters
        std::string fRawDLabel;
        bool fDebug;
        int fNSigma;

        // var
        int fEvent;
        int fRun;
        int fSubRun;

        std::vector< TH1D* > fGetPedMeans;
        std::vector< TH1D* > fGetPedRmss;
        std::vector< TH1D* > fMedianPedMeans;
        std::vector< TH1D* > fMedianPedRmss;
        std::vector< TH1D* > fGausPedMeans;
        std::vector< TH1D* > fGausPedRmss;
        std::vector< TH1D* > fModePedMeans;
        //std::vector< TH1D* > fRetPedMeans;
        //std::vector< TH1D* > fRetPedRmss;

        TH2D* fMedianvGetPed;
        TH1D* fGetMinusMedianPed;
        TH2D* fGausvGetPed;
        TH1D* fGetMinusGausPed;
        TH2D* fModevGetPed;
        TH1D* fGetMinusModePed;
        //TH2D* fGetvRetPed;
        //TH1D* fRetMinusGetPed;

        int channel;

        // TTree
        TTree *pedestalInfo;

};


ExtractPedestalInformation::ExtractPedestalInformation(fhicl::ParameterSet const& p)
    :
        EDAnalyzer(p)  // ,
        // More initializers here.
{

    reconfigure(p);

}

//
// beginJob
//

void ExtractPedestalInformation::beginJob()
{

    art::ServiceHandle< art::TFileService > tfs;

    pedestalInfo = tfs->make<TTree>("pedestalInfo", "pedestalInfo");

    // TFileDirectories
    art::TFileDirectory modePedMeanDir   = tfs->mkdir("modePedMeanDir");
    art::TFileDirectory gausPedMeanDir   = tfs->mkdir("gaussPedMeanDir");
    art::TFileDirectory gausPedRmsDir    = tfs->mkdir("gausPedRmsDir");
    art::TFileDirectory medianPedMeanDir = tfs->mkdir("medianPedMeanDir");
    art::TFileDirectory medianPedRmsDir  = tfs->mkdir("medianPedRmsDir");
    art::TFileDirectory getPedMeanDir    = tfs->mkdir("getPedMeanDir");
    art::TFileDirectory getPedRmsDir     = tfs->mkdir("getPedRmsDir");

    // define branches
    pedestalInfo->Branch("Event", &fEvent, "Event/I");
    pedestalInfo->Branch("Run", &fRun, "Run/I");
    pedestalInfo->Branch("SubRun", &fSubRun, "SubRun/I");

    int NBINS = 1000;
    int LOWPED = 400;
    int HIGHPED = 500;

    // histograms
    fModevGetPed       = tfs->make<TH2D>("fModevGetPed", ";Pedestal from Get; Calculated pedestal (mode)", NBINS, LOWPED, HIGHPED, NBINS, LOWPED, HIGHPED);
    fGetMinusModePed   = tfs->make<TH1D>("fGetMinusModePed", ";|Pedestal from Get - Calculated pedestal (mode)|;", NBINS, 0, 10);
    fMedianvGetPed     = tfs->make<TH2D>("fMedianvGetPed", ";Pedestal from Get; Calculated pedestal (median);", NBINS, LOWPED, HIGHPED, NBINS, LOWPED, HIGHPED);
    fGetMinusMedianPed = tfs->make<TH1D>("fGetMinusMedianPed", ";|Pedestal from Get - Calculated pedestal (median)|;", NBINS, 0, 10);
    fGausvGetPed       = tfs->make<TH2D>("fGausvGetPed", ";Pedestal from Get; Calculated pedestal (Gaussian fit);", NBINS, LOWPED, HIGHPED, NBINS, LOWPED, HIGHPED);
    fGetMinusGausPed   = tfs->make<TH1D>("fGetMinusGausPed", ";|Pedestal from Get - Calculated pedestal (median)|;", NBINS, 0, 10);
    //fGetvRetPed        = tfs->make<TH2D>("fGetvRetPed", ";Pedestal from //DB; Pedestal from rawd.Get()", 1000, LOWPED, HIGHPED, 1000, LOWPED, HIGHPED);
    //fRetMinusGetPed    = tfs->make<TH1D>("fRetMinusGetPed", ";|Pedestal from DB - pedestal from rawd.Get()|;", 1000, 0, 10);

    // produce one histogram per channel in a TDirectoryFile

    for (int i = 0; i <8256; i++){

        TString modePedMeanName = Form("ModePedMean_Channel%i", i);
        fModePedMeans.push_back(modePedMeanDir.make<TH1D>(modePedMeanName, ";Pedestal mean (using mode);", NBINS, LOWPED, HIGHPED));

        TString gausPedMeanName = Form("GausPedMean_Channel%i", i);
        fGausPedMeans.push_back(gausPedMeanDir.make<TH1D>(gausPedMeanName, ";Pedestal mean (from Gaussian);", NBINS, LOWPED, HIGHPED));

        TString gausPedRmsName = Form("GausPedRmss_Channel%i", i);
        fGausPedRmss.push_back(gausPedRmsDir.make<TH1D>(gausPedRmsName, ";Pedestal RMS (from Gaussian);", 200, 0, 5));

        TString medianPedMeanName = Form("MedianPedMean_Channel%i", i);
        fMedianPedMeans.push_back(medianPedMeanDir.make<TH1D>(medianPedMeanName, ";Estimated pedestal mean (using median);", NBINS, LOWPED, HIGHPED));

        TString medianPedRmssName = Form("MedianPedRmss_Channel%i", i);
        fMedianPedRmss.push_back(medianPedRmsDir.make<TH1D>(medianPedRmssName, ";Estimated pedestal RMS (using median);", 200, 0, 5));

        TString getPedMeanName = Form("GetPedMean_Channel%i", i);
        fGetPedMeans.push_back(getPedMeanDir.make<TH1D>(getPedMeanName, ";Pedestal Mean (from rawd.GetPedestal());", NBINS, LOWPED, HIGHPED));

        TString getPedRmssName = Form("GetPedRmss_Channel%i", i);
        fGetPedRmss.push_back(getPedRmsDir.make<TH1D>(getPedRmssName, ";Pedestal RMS (from rawd.GetPedestal());", 200, 0, 5));

    }

}


//
// reconfigure
//

void ExtractPedestalInformation::reconfigure(fhicl::ParameterSet const& pset)
{

    fDebug     = pset.get< bool > ("Debug");
    fRawDLabel = pset.get< std::string > ("RawDLabel");
    fNSigma    = pset.get< int > ("NSigmaSignalRejection");
}

//
// analyze
//

void ExtractPedestalInformation::analyze(art::Event const & event)
{

    std::cout << "beginning" << std::endl;

    fEvent = event.id().event();
    fRun = event.run();
    fSubRun = event.subRun();

    //get pedestal conditions
    //const lariov::DetPedestalProvider& pedestalRetrievalAlg 
    //    = art::ServiceHandle<lariov::DetPedestalService>()->GetPedestalProvider();

    //float pedMeanFromDb = 0;
    //float pedRmsFromDb  = 0;

    // get handle to rawdigit
    art::Handle< std::vector<raw::RawDigit> > rawDHandle;
    event.getByLabel(fRawDLabel, rawDHandle);

    for ( auto const& rawd : (*rawDHandle) ){

        if (!rawd.Channel()) continue;

        channel = (int)rawd.Channel();
        int nTicks = rawd.NADC();

        // extract pedestal values from database
        //pedMeanFromDb = pedestalRetrievalAlg.PedMean(channel);
        //pedRmsFromDb = pedestalRetrievalAlg.PedRms(channel);

        // get adc vals and fill a histogram in order to calculate median and
        // signal removed RMS of waveform
        const std::vector<short int> adcv = rawd.ADCs();

        TH1F *h_adc = new TH1F("h_adc", "", 2600, 400, 3000);

        for (int i=0; i < nTicks; i++)
        {

            h_adc->Fill(adcv.at(i));

        }

        //
        // find the median value of the ADCs
        //

        //Manual calculation of RMS on h_adc
        double pars[3];
        double calculatedRms = 0;
        double calculatedMedian = 0;

        if (h_adc->GetSum()>0){

            double xc = 0.5-0.34;
            h_adc->GetQuantiles(1, &pars[0], &xc);

            xc = 0.5;
            h_adc->GetQuantiles(1, &pars[1], &xc);

            xc = 0.5+0.34;
            h_adc->GetQuantiles(1, &pars[2], &xc);

            calculatedRms = sqrt((pow(pars[1]-pars[0],2) + pow(pars[2]-pars[1],2))/2.);

            calculatedMedian = pars[1];


        }

        //
        // Remove ADC values that are more than 4 sigma away from mean value
        //

        double adcSigma = calculatedRms;
        if (adcSigma == 0) adcSigma = 1;

        double lowerBound = calculatedMedian - ( adcSigma * fNSigma );
        double upperBound = calculatedMedian + ( adcSigma * fNSigma );

        for (int i = 0; i < h_adc->GetNbinsX(); i++){

            if (h_adc->GetBinContent(i) > 0  &&  
                    ((h_adc->GetBinLowEdge(i) < lowerBound ) || 
                     ((h_adc->GetBinLowEdge(i) + h_adc->GetBinWidth(i)) > upperBound ))){
                h_adc->SetBinContent(i, 0);
            }

        }

        //
        // calculate mode
        //

        double modeVal = h_adc->GetBinLowEdge(h_adc->GetMaximumBin());

        //
        // fit Gaussian to remaining ADC values
        //

        double gausMean = 0;
        double gausRms = 0;

        //TF1 *gaussian = new TF1("gaussian", "[0]*exp(-0.5*((x-[1])/[2])**2)", 400, 3000);
        //gaussian->SetParameters(10, 470, 0.5);

        if (h_adc->GetMaximum() > 0){
        
            h_adc->Fit("gaus", "Q");
            gausMean = h_adc->GetFunction("gaus")->GetParameter(1);
            gausRms = h_adc->GetFunction("gaus")->GetParameter(2);
        
        }

        float truePedMean = rawd.GetPedestal();
        float truePedRmss = rawd.GetSigma();

        fModePedMeans.at(channel)->Fill(modeVal);
        fModevGetPed->Fill(truePedMean, modeVal);
        fGetMinusModePed->Fill(std::abs(truePedMean - modeVal));
        
        fGausPedMeans.at(channel)->Fill(gausMean - 0.5);
        fGausPedRmss.at(channel)->Fill(gausRms);
        fGausvGetPed->Fill(truePedMean, gausMean-0.5);
        fGetMinusGausPed->Fill(std::abs(truePedMean - gausMean+0.5));
        
        fGetPedMeans.at(channel)->Fill(truePedMean);
        fGetPedRmss.at(channel)->Fill(truePedRmss);
       
        fMedianPedMeans.at(channel)->Fill(calculatedMedian);
        fMedianvGetPed->Fill(truePedMean, calculatedMedian-0.5);
        fGetMinusMedianPed->Fill(std::abs(truePedMean - calculatedMedian+0.5));       
        fMedianPedRmss.at(channel)->Fill(calculatedRms);
        
        if (fDebug){

             std::cout << "\nRemoving ADC values above " << upperBound 
                << " and below " << lowerBound  
                << "\nSum of histo is " << h_adc->GetSum()
                << "\nTrue pedestal is " << truePedMean
                //<< "\nFrom DB; " << pedMeanFromDb 
                << "\nMedian value is " << calculatedMedian 
                << "\nModal value is " << modeVal
                << "\nGaussian Mean value is " << gausMean
                << "\nadcSigma:" << adcSigma
                << "\nfNSigma: " << fNSigma << std::endl;
        
        }

        delete h_adc;
        pedestalInfo->Fill();
    
    } // rawdigits


}

void ExtractPedestalInformation::endJob()
{

}


DEFINE_ART_MODULE(ExtractPedestalInformation)

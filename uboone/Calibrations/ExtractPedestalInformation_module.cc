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
        std::vector< TH1D* > fGausPedMeans;
        std::vector< TH1D* > fGausPedRmss;

        TH1D* fGetMinusGausPed;

        int channel;
        float truePedMean[8256];
        float truePedRmss[8256];


        // TTree
        TTree *pedestalInfo;

        long long int timeFromEpoch = 9223372036854775807;
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
    art::TFileDirectory gausPedMeanDir   = tfs->mkdir("gausPedMeanDir");
    art::TFileDirectory gausPedRmsDir    = tfs->mkdir("gausPedRmsDir");
    art::TFileDirectory getPedMeanDir    = tfs->mkdir("getPedMeanDir");
    art::TFileDirectory getPedRmsDir     = tfs->mkdir("getPedRmsDir");

    // define branches
    pedestalInfo->Branch("Event", &fEvent, "Event/I");
    pedestalInfo->Branch("Run", &fRun, "Run/I");
    pedestalInfo->Branch("SubRun", &fSubRun, "SubRun/I");

    int NBINS        = 10000;
    int LOWPED_COLL  = 400;
    int HIGHPED_COLL = 500;
    int NCHANS       = 8256;
    int LOWPED_IND   = 2000;
    int HIGHPED_IND  = 2100;

    // histograms
    fGetMinusGausPed   = tfs->make<TH1D>("fGetMinusGausPed", ";|Pedestal from Get - Calculated pedestal (median)|;", NBINS, 0, 10);

    // produce one histogram per channel in a TDirectoryFile

    for (int i = 1; i <= NCHANS; i++){

        TString gausPedRmsName = Form("GausPedRmss_Channel%i", i);
        fGausPedRmss.push_back(gausPedRmsDir.make<TH1D>(gausPedRmsName, ";Pedestal RMS (from Gaussian);", 200, 0, 5));

        TString getPedRmssName = Form("GetPedRmss_Channel%i", i);
        fGetPedRmss.push_back(getPedRmsDir.make<TH1D>(getPedRmssName, ";Pedestal RMS (from rawd.GetPedestal());", 200, 0, 5));


        if (i > 4800){

            TString gausPedMeanName = Form("GausPedMean_Channel%i", i);
            fGausPedMeans.push_back(gausPedMeanDir.make<TH1D>(gausPedMeanName, ";Pedestal mean (from Gaussian);", NBINS, LOWPED_COLL, HIGHPED_COLL));

            TString getPedMeanName = Form("GetPedMean_Channel%i", i);
            fGetPedMeans.push_back(getPedMeanDir.make<TH1D>(getPedMeanName, ";Pedestal Mean (from rawd.GetPedestal());", NBINS, LOWPED_COLL, HIGHPED_COLL));


        }
        else{

            TString gausPedMeanName = Form("GausPedMean_Channel%i", i);
            fGausPedMeans.push_back(gausPedMeanDir.make<TH1D>(gausPedMeanName, ";Pedestal mean (from Gaussian);", NBINS, LOWPED_IND, HIGHPED_IND));

            TString getPedMeanName = Form("GetPedMean_Channel%i", i);
            fGetPedMeans.push_back(getPedMeanDir.make<TH1D>(getPedMeanName, ";Pedestal Mean (from rawd.GetPedestal());", NBINS, LOWPED_IND, HIGHPED_IND));

        }

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

    if ((long long int)event.time().value() < timeFromEpoch){
        timeFromEpoch = (long long int)event.time().value();
    }
    fEvent = event.id().event();
    fRun = event.run();
    fSubRun = event.subRun();

    // get handle to rawdigit
    art::Handle< std::vector<raw::RawDigit> > rawDHandle;
    event.getByLabel(fRawDLabel, rawDHandle);

    int chan_counter = 1;
    for ( auto const& rawd : (*rawDHandle) ){

        if (!rawd.Channel()) continue;

        channel = (int)rawd.Channel();
        int nTicks = rawd.NADC();

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
        // fit Gaussian to remaining ADC values
        //

        double gausMean = 0;
        double gausRms = 0;

        if (h_adc->GetMaximum() > 0){

            h_adc->Fit("gaus", "Q");
            gausMean = h_adc->GetFunction("gaus")->GetParameter(1);
            gausRms = h_adc->GetFunction("gaus")->GetParameter(2);

        }

        truePedMean[chan_counter-1] = rawd.GetPedestal();
        truePedRmss[chan_counter-1] = rawd.GetSigma();

        fGausPedMeans.at(chan_counter-1)->Fill(gausMean - 0.5);
        fGausPedRmss.at(chan_counter-1)->Fill(gausRms);

        fGetPedMeans.at(channel-1)->Fill(truePedMean[chan_counter-1]);
        fGetPedRmss.at(channel -1)->Fill(truePedRmss[chan_counter-1]);

        if (fDebug){

            std::cout << "\nRemoving ADC values above " << upperBound 
                << " and below " << lowerBound  
                << "\nSum of histo is " << h_adc->GetSum()
                << "\nTrue pedestal is " << truePedMean[chan_counter-1]
                << "\nGaussian Mean value is " << gausMean
                << "\nadcSigma:" << adcSigma
                << "\nfNSigma: " << fNSigma << std::endl;

        }

        delete h_adc;
        pedestalInfo->Fill();

        chan_counter++;
    } // rawdigits


}

void ExtractPedestalInformation::endJob()
{

    std::ofstream pedValues;
    pedValues.open("pedestalValues.dat");
    pedValues << "# time " << timeFromEpoch << std::endl;
    pedValues << "# " << "channel mean mean_err rms rms_err" << std::endl;

    float gaussianMean = -1;
    float gaussianMeanEr = -1;
    float stddev = 0;
    float stddevEr = 0;
    for (int chan_it = 1; chan_it <= 8256; chan_it++){

        TString chanNo = Form("Channel%i", chan_it);

        TH1D *hGaus = (TH1D*)((TFile*)pedestalInfo->GetCurrentFile())->Get("getpedinfo/gausPedMeanDir/GausPedMean_"+chanNo);

        if (hGaus->Integral() > 0){

            hGaus->Fit("gaus", "Q");
            gaussianMean = hGaus->GetFunction("gaus")->GetParameter(1);
            gaussianMeanEr = hGaus->GetFunction("gaus")->GetParError(1);
            stddev = hGaus->GetFunction("gaus")->GetParameter(2);
            stddevEr = hGaus->GetFunction("gaus")->GetParError(2);

        }

        pedValues << std::fixed;
        pedValues << chan_it << " " << std::setprecision(1) << gaussianMean << " " << std::setprecision(1) << gaussianMeanEr << " " << std::setprecision(2) << stddev << " " <<  std::setprecision(2) << stddevEr << std::endl;

        fGetMinusGausPed->Fill(std::abs(truePedMean[chan_it-1] - gaussianMean));

    }



}


DEFINE_ART_MODULE(ExtractPedestalInformation)

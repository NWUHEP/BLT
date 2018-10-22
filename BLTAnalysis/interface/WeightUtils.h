/*
   Utilities for retrieving weights for PU,etc.
   */

#ifndef _WeightUtils_H
#define _WeightUtils_H

// c++ libraries
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

using namespace std;

class EfficiencyContainer: public TObject{
    public:
        EfficiencyContainer();
        EfficiencyContainer(float, float, float, float);
        virtual ~EfficiencyContainer() {};
        void SetData(float, float, float, float);
        pair<double, double> GetEff() {return make_pair(_dataEff, _mcEff);};
        pair<double, double> GetErr() {return make_pair(_dataErr, _mcErr);};
        
    private:
        float _dataEff, _mcEff;
        float _dataErr, _mcErr;
};

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string dataPeriod, string selection, bool isRealData);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSelection(string);

        float   GetPUWeight(float);
        EfficiencyContainer GetTriggerEffWeight(string, TLorentzVector&) const;
        EfficiencyContainer GetMuonRecoEff(TLorentzVector&) const; 
        EfficiencyContainer GetElectronRecoEff(TLorentzVector&) const;

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        // pileup
        TGraph  *_puReweight;

        // muon trigger and ID/ISO scale factors
        TGraphAsymmErrors *_muSF_IsoMu24_DATA_BCDEF[4], *_muSF_IsoMu24_MC_BCDEF[4];
        TGraphAsymmErrors *_muSF_ID_DATA_BCDEF[4], *_muSF_ID_MC_BCDEF[4];
        TGraphAsymmErrors *_muSF_ISO_DATA_BCDEF[4], *_muSF_ISO_MC_BCDEF[4]; 

        TGraphAsymmErrors *_muSF_IsoMu24_DATA_GH[4], *_muSF_IsoMu24_MC_GH[4];  
        TGraphAsymmErrors *_muSF_ID_DATA_GH[4], *_muSF_ID_MC_GH[4]; 
        TGraphAsymmErrors *_muSF_ISO_DATA_GH[4], *_muSF_ISO_MC_GH[4]; 

        // electron RECO/ID scale factors (what about ISO?)
        TGraphErrors *_eleSF_RECO[5], *_eleSF_ID[5];

        // electron trigger efficiencies (the bins for 2.1 < |eta| < 2.4 are copies of the 1.6 to 2.1 bins 
        TH2D *_elSF_Trigger_BCDEF, *_elSF_Trigger_GH;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif

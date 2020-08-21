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
        float   GetWWPtWeight(float, float&, float&, float&, float&);
        EfficiencyContainer GetTriggerEffWeight(string, TLorentzVector&) const;
        EfficiencyContainer GetMuonIDEff(TLorentzVector&) const; 
        EfficiencyContainer GetMuonISOEff(TLorentzVector&) const; 
        EfficiencyContainer GetElectronRecoEff(TLorentzVector&) const;
        EfficiencyContainer GetElectronIDEff(TLorentzVector&) const;
        EfficiencyContainer GetElectronPrefiringWeight(vector<TLorentzVector>, vector<TLorentzVector>) const;
        double  GetPrefiringRate(double, double, TH2D*, int) const;
        float GetEleTriggerSyst(string, TLorentzVector&) const;

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        // pileup
        TGraph  *_puReweight;

        // ww pt weights
        TH1D *_wwPtWeight, *_wwPtScaleErrUp, *_wwPtScaleErrDown, *_wwPtResumErrUp, *_wwPtResumErrDown;

        // muon trigger and ID/ISO scale factors
        TGraphAsymmErrors *_muSF_IsoMu24_DATA_BCDEF[4], *_muSF_IsoMu24_MC_BCDEF[4];
        TGraphAsymmErrors *_muSF_ID_DATA_BCDEF[4], *_muSF_ID_MC_BCDEF[4];
        TGraphAsymmErrors *_muSF_ISO_DATA_BCDEF[4], *_muSF_ISO_MC_BCDEF[4]; 

        TGraphAsymmErrors *_muSF_IsoMu24_DATA_GH[4], *_muSF_IsoMu24_MC_GH[4];  
        TGraphAsymmErrors *_muSF_ID_DATA_GH[4], *_muSF_ID_MC_GH[4]; 
        TGraphAsymmErrors *_muSF_ISO_DATA_GH[4], *_muSF_ISO_MC_GH[4]; 

        // electron RECO/ID scale factors (id includes isolation)
        TGraphErrors *_eleSF_RECO, *_eleSF_ID[5];

        // electron trigger efficiencies from Andrey Popov (https://indico.cern.ch/event/604912/)
        TH2D *_elSF_Trigger_BCDEF, *_elSF_Trigger_GH;
        TH2D *_elSF_Trigger_BCDEF_tag, *_elSF_Trigger_GH_tag;
        TH2D *_elSF_Trigger_BCDEF_probe, *_elSF_Trigger_GH_probe;

         // electron prefiring weights
        TH2D *_elSF_Prefiring_photon, *_elSF_Prefiring_jet;
        double prefiringRateSystUnc_ = 0.2; // official value https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1ECALPrefiringWeightProducer.cc#L221
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif

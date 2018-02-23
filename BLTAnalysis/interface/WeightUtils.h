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
        float   GetMuonIDEff(TLorentzVector&) const; 
        float   GetMuonISOEff(TLorentzVector&) const; 
        float   GetElectronRecoIdEff(TLorentzVector&) const;

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
        TGraphErrors *_eleSF_RECO, *_eleSF_ID[5];

        // electron trigger efficiencies (the bins for 2.1 < |eta| < 2.4 are copies of the 1.6 to 2.1 bins 
        float _elePtBins[8] = {30, 32, 35, 40, 50, 60, 120, 9999};
        float _eleEtaBins[13] = {-2.4, -2.1, -1.6, -1.4, -0.8, -0.4, 0., 0.4, 0.8, 1.4, 1.6, 2.1, 2.4};
        float _ele_trigEff_data[12][7] = {
            {0.686, 0.747, 0.798, 0.838, 0.855, 0.866, 0.884}, // not real
            {0.686, 0.747, 0.798, 0.838, 0.855, 0.866, 0.884},
            {0.714, 0.772, 0.822, 0.865, 0.882, 0.888, 0.922},
            {0.809, 0.859, 0.893, 0.924, 0.939, 0.954, 0.965},
            {0.789, 0.851, 0.885, 0.922, 0.938, 0.956, 0.970},
            {0.710, 0.782, 0.829, 0.882, 0.909, 0.940, 0.965},
            {0.711, 0.779, 0.827, 0.876, 0.906, 0.931, 0.959},
            {0.796, 0.853, 0.884, 0.918, 0.935, 0.950, 0.965},
            {0.805, 0.860, 0.896, 0.927, 0.942, 0.957, 0.971},
            {0.683, 0.763, 0.814, 0.864, 0.880, 0.880, 0.896},
            {0.676, 0.743, 0.797, 0.841, 0.860, 0.877, 0.896},
            {0.676, 0.743, 0.797, 0.841, 0.860, 0.877, 0.896} // not real
        };
        float _ele_trigEff_mc[12][7] = {
            {0.796, 0.833, 0.879, 0.909, 0.920, 0.940, 0.963}, // not real
            {0.796, 0.833, 0.879, 0.909, 0.920, 0.940, 0.963},
            {0.792, 0.836, 0.867, 0.904, 0.919, 0.933, 0.954},
            {0.845, 0.886, 0.904, 0.928, 0.947, 0.974, 0.980},
            {0.831, 0.864, 0.882, 0.912, 0.941, 0.964, 0.979},
            {0.789, 0.821, 0.848, 0.882, 0.915, 0.951, 0.985},
            {0.794, 0.815, 0.846, 0.884, 0.917, 0.949, 0.979},
            {0.833, 0.861, 0.881, 0.913, 0.942, 0.964, 0.988},
            {0.853, 0.886, 0.906, 0.930, 0.953, 0.973, 0.990},
            {0.783, 0.821, 0.867, 0.903, 0.918, 0.936, 0.940},
            {0.802, 0.848, 0.889, 0.914, 0.929, 0.944, 0.975},
            {0.802, 0.848, 0.889, 0.914, 0.929, 0.944, 0.975} // not real
        };

};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif

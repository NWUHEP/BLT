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
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

// boost libraries
//#include <boost/array.hpp>

// custom libraries
#include "../interface/TCPhysObject.h"
#include "../interface/TCPhoton.h"
#include "../interface/TCJet.h"
#include "../interface/TCGenJet.h"

using namespace std;

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string sampleName, string dataPeriod, string selection, bool isRealData);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSampleName(string);
        void    SetSelection(string);
        void    SetPassTrigger(string); 
        void    SetObjects(vector<TCPhysObject>&, vector<TCJet>&, float, string);

        float   GetPUWeight(unsigned);
        float   RecoWeight();
        float   GetElectronEff(TLorentzVector&) const;
        float   GetMuEff(TLorentzVector&) const; 
        float   GetMuTriggerEff(TLorentzVector&, TLorentzVector&) const;
        float   GetEleTriggerEff(TLorentzVector&, TLorentzVector&) const;

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        //sources
        TH1D*  puReweight;

        TGraphErrors *_muSF2012_ID[4], *_muSF2012_ISO[4];
        TGraph *_muSF2012_ID_err[4], *_muSF2012_ISO_err[4];

        TH2D    *h2_MuTriggerSFs[2]; // Good for Mu17_Mu8 or Mu17_TkMu8
        TH2D    *h2_EleMVASF;
        TH2D    *h2_DielectronMisQ;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif

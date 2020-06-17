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

// BaconAna class definitions
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"

using namespace std;
using namespace baconhep;

class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string dataPeriod, string selection);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSelection(string);

        std::map<std::string, float>  GetPUWeight(float);
        pair<float, float>   GetTriggerEffWeight(string, TLorentzVector&) const;
        //pair<float, float>   GetDoubleEGTriggerEffWeight(string, TElectron&) const;
        //pair<float, float>   GetDoubleMuonTriggerEffWeight(string, TMuon&) const;
        //float   GetDoubleEGTriggerEffWeight(string, TElectron&) const;
        pair<float, float> GetDoubleEGTriggerEffWeight(string, TElectron&) const;
        pair<float, float> GetSingleElectronTriggerEffWeight(TElectron&) const;
        float   GetDoubleMuonTriggerEffWeight(string, TMuon&) const;
        float   GetMuonIDEff(TLorentzVector&) const; 
        float   GetMuonISOEff(TLorentzVector&) const; 
        float   GetLooseMuonIDEff(TLorentzVector&) const;
        float   GetLooseMuonISOEff(TLorentzVector&) const;
        float   GetHZZMuonIDEff(TMuon&) const;
        float   GetElectronRecoIdEff(TLorentzVector&) const;
        float   GetHZZElectronRecoIdEff(TElectron&) const;
        pair<float, float> GetElectronMVARecoEff(TElectron&) const;
        pair<float, float> GetElectronMVAIdEff(TElectron&, string) const;
        //float   GetElectronMVARecoIdEff(TElectron&) const;
        pair<float, float> GetPhotonMVAIdEff(TPhoton&) const;
        pair<float, float> GetPhotonMVACSEVEff(TPhoton &) const;
        float   GetCorrectedPhotonR9(TPhoton&) const;

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;

        // rng
        TRandom3 *rng;

        TH1D  *_puReweightNom, *_puReweightUp, *_puReweightDown;
        TGraphAsymmErrors *_eff_IsoMu24_DATA[4]; 
        TGraphAsymmErrors *_muSF_ID_DATA_BCDEF[4], *_muSF_ID_MC_BCDEF[4]; 
        TGraphAsymmErrors *_muSF_ISO_DATA_BCDEF[4], *_muSF_ISO_MC_BCDEF[4]; 
        TGraphAsymmErrors *_muSF_ID_DATA_GH[4], *_muSF_ID_MC_GH[4]; 
        TGraphAsymmErrors *_muSF_ISO_DATA_GH[4], *_muSF_ISO_MC_GH[4]; 
        TGraphAsymmErrors *_muSF_Loose_ID_DATA_BCDEF[4], *_muSF_Loose_ID_MC_BCDEF[4]; 
        TGraphAsymmErrors *_muSF_Loose_ISO_DATA_BCDEF[4], *_muSF_Loose_ISO_MC_BCDEF[4]; 
        TGraphAsymmErrors *_muSF_Loose_ID_DATA_GH[4], *_muSF_Loose_ID_MC_GH[4]; 
        TGraphAsymmErrors *_muSF_Loose_ISO_DATA_GH[4], *_muSF_Loose_ISO_MC_GH[4]; 
        TGraphAsymmErrors *_eff_doubleMu_leg1_DATA[4], *_eff_doubleMu_leg1_MC[4];
        TGraphAsymmErrors *_eff_doubleMu_leg2_DATA[4], *_eff_doubleMu_leg2_MC[4];
        TH1F *_sf_doubleMu_leg1[4];
        TH1F *_sf_doubleMu_leg2[4];

        TGraphErrors *_eleSF_RECO, *_eleSF_ID[5], *_hzz_eleSF_ID[13];
        TH2F *_eleSF_RECO_2D, *_hzz_eleSF_ID_2D;
        TH2F *_eleSF_MVA_RECO_2D, *_eleSF_MVA_LOW_RECO_2D;
        TH2F *_eleSF_MVA_ID_2D_WP80, *_eleSF_MVA_ID_2D_WP90, *_eleSF_MVA_ID_2D_WPLoose;

        TGraphErrors *_mva_gammaSF_ID[5];
        TH2F *_mva_gammaSF;

        TGraph *_photon_r9_barrel, *_photon_r9_endcap;

        TH2F *_eff_doubleg_leg1_DATA, *_eff_doubleg_leg1_MC;
        TH2F *_eff_doubleg_leg2_DATA, *_eff_doubleg_leg2_MC;
        TH2F *_sf_doubleg_leg1, *_sf_doubleg_leg2;
        TH2F *_sf_single_electron_trig;

        TH2F *_hzz_muIdSF;

        //TH2D    *h2_MuTriggerSFs[2]; // Good for Mu17_Mu8 or Mu17_TkMu8
        //TH2D    *h2_EleMVASF;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif

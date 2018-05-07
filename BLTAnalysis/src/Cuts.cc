#include "BLT/BLTAnalysis/interface/Cuts.hh"

#include <iostream>
#include <stdexcept>


Cuts::Cuts() {
    leadJetPt   = 30;
    trailJetPt  = 30;
    leadMuPt    = 20;
    trailMuPt   = 10;
    leadElPt    = 20;
    trailElPt   = 10;
    gPtOverMass = 15./100.;
    gPt         = 15.;
    zMassLow    = 50;
    zMassHigh   = 999999;
    metLow      = -999999;
    metHigh     = 999999;
    zgMassLow   = 100;
    zgMassHigh  = 190;

    //combined rel ISO, 2012 Data, 0.5 GeV
    EAMu[0] =   0.674; //         eta < 1.0
    EAMu[1] =   0.565; // 1.0   < eta < 1.5
    EAMu[2] =   0.442; // 1.5   < eta < 2.0
    EAMu[3] =   0.515; // 2.0   < eta < 2.2
    EAMu[4] =   0.821; // 2.2   < eta < 2.3
    EAMu[5] =   0.660; // 2.3   < eta < 2.4

    EAEl[0] =   0.13; //         eta < 1.0
    EAEl[1] =   0.14; // 1.0   < eta < 1.5
    EAEl[2] =   0.07; // 1.5   < eta < 2.0
    EAEl[3] =   0.09; // 2.0   < eta < 2.2
    EAEl[4] =   0.11; // 2.2   < eta < 2.3
    EAEl[5] =   0.11; // 2.3   < eta < 2.4
    EAEl[6] =   0.14; // 2.4   < eta

    //   ch      nh       ph
    float EAPhoTemp[7][3] = {
        {0.012,  0.030,   0.148}, //         eta < 1.0
        {0.010,  0.057,   0.130}, // 1.0   < eta < 1.479
        {0.014,  0.039,   0.112}, // 1.479 < eta < 2.0
        {0.012,  0.015,   0.216}, // 2.0   < eta < 2.2
        {0.016,  0.024,   0.262}, // 2.2   < eta < 2.3
        {0.020,  0.039,   0.260}, // 2.3   < eta < 2.4
        {0.012,  0.072,   0.266}  // 2.4   < eta
    };

    for (unsigned int i =0; i<7; i++)
        for (unsigned int j =0; j<3; j++)
            EAPho[i][j] = EAPhoTemp[i][j];


    vetoElID.cutName                      = "vetoElID";
    vetoElID.dEtaIn[0]                    = 0.007;
    vetoElID.dPhiIn[0]                    = 0.8;
    vetoElID.sigmaIetaIeta[0]             = 0.01;
    vetoElID.HadOverEm[0]                 = 0.15;
    vetoElID.dxy[0]                       = 0.04;
    vetoElID.dz[0]                        = 0.2;
    vetoElID.fabsEPDiff[0]                = 99999;
    vetoElID.ConversionMissHits[0]        = 99999;

    vetoElID.dEtaIn[1]                    = 0.01;
    vetoElID.dPhiIn[1]                    = 0.7;
    vetoElID.sigmaIetaIeta[1]             = 0.03;
    vetoElID.HadOverEm[1]                 = 999999;
    vetoElID.dxy[1]                       = 0.04;
    vetoElID.dz[1]                        = 0.2;
    vetoElID.fabsEPDiff[1]                = 999999;
    vetoElID.ConversionMissHits[1]        = 999999;

    // tight electorn ID
    tightElID.cutName                    = "tightElID";
    tightElID.dEtaIn[0]                  = 0.004;
    tightElID.dPhiIn[0]                  = 0.03;
    tightElID.sigmaIetaIeta[0]           = 0.01;
    tightElID.HadOverEm[0]               = 0.12;
    tightElID.dxy[0]                     = 0.02;
    tightElID.dz[0]                      = 0.1;
    tightElID.fabsEPDiff[0]              = 0.05;
    tightElID.ConversionMissHits[0]      = 0;

    tightElID.dEtaIn[1]                  = 0.007;
    tightElID.dPhiIn[1]                  = 0.03;
    tightElID.sigmaIetaIeta[1]           = 0.03;
    tightElID.HadOverEm[1]               = 0.10;
    tightElID.dxy[1]                     = 0.02;
    tightElID.dz[1]                      = 0.1;
    tightElID.fabsEPDiff[1]              = 0.05;
    tightElID.ConversionMissHits[1]      = 1;

    // medium electron ID
    mediumElID.cutName                    = "mediumElID";
    mediumElID.dEtaIn[0]                  = 0.004;
    mediumElID.dPhiIn[0]                  = 0.06;
    mediumElID.sigmaIetaIeta[0]           = 0.01;
    mediumElID.HadOverEm[0]               = 0.12;
    mediumElID.dxy[0]                     = 0.02;
    mediumElID.dz[0]                      = 0.1;
    mediumElID.fabsEPDiff[0]              = 0.05;
    mediumElID.ConversionMissHits[0]      = 1;

    mediumElID.dEtaIn[1]                  = 0.007;
    mediumElID.dPhiIn[1]                  = 0.03;
    mediumElID.sigmaIetaIeta[1]           = 0.03;
    mediumElID.HadOverEm[1]               = 0.10;
    mediumElID.dxy[1]                     = 0.02;
    mediumElID.dz[1]                      = 0.1;
    mediumElID.fabsEPDiff[1]              = 0.05;
    mediumElID.ConversionMissHits[1]      = 1;

    looseElID.cutName                     = "looseElID";
    looseElID.dEtaIn[0]                   = 0.007;
    looseElID.dPhiIn[0]                   = 0.15;
    looseElID.sigmaIetaIeta[0]            = 0.01;
    looseElID.HadOverEm[0]                = 0.12;
    looseElID.dxy[0]                      = 0.02;
    looseElID.dz[0]                       = 0.2;
    looseElID.fabsEPDiff[0]               = 0.05;
    looseElID.ConversionMissHits[0]       = 1;

    looseElID.dEtaIn[1]                   = 0.009;
    looseElID.dPhiIn[1]                   = 0.10;
    looseElID.sigmaIetaIeta[1]            = 0.03;
    looseElID.HadOverEm[1]                = 0.10;
    looseElID.dxy[1]                      = 0.02;
    looseElID.dz[1]                       = 0.2;
    looseElID.fabsEPDiff[1]               = 0.05;
    looseElID.ConversionMissHits[1]       = 1;

    mvaPreElID.cutName                    = "mvaPreElID";
    mvaPreElID.dEtaIn[0]                  = 99999;
    mvaPreElID.dPhiIn[0]                  = 9999;
    mvaPreElID.sigmaIetaIeta[0]           = 0.014;
    mvaPreElID.HadOverEm[0]               = 0.15;
    mvaPreElID.dxy[0]                     = 99999;
    mvaPreElID.dz[0]                      = 99999;
    mvaPreElID.fabsEPDiff[0]              = 99999;
    mvaPreElID.ConversionMissHits[0]      = 99999;
    mvaPreElID.dr03TkSumPt[0]             = 0.2;
    mvaPreElID.dr03EcalRecHitSumEt[0]     = 0.2;
    mvaPreElID.dr03HcalTowerSumEt[0]      = 0.2;
    mvaPreElID.numberOfLostHits[0]        = 0;

    mvaPreElID.dEtaIn[1]                  = 99999;
    mvaPreElID.dPhiIn[1]                  = 99999;
    mvaPreElID.sigmaIetaIeta[1]           = 0.035;
    mvaPreElID.HadOverEm[1]               = 0.10;
    mvaPreElID.dxy[1]                     = 99999;
    mvaPreElID.dz[1]                      = 99999;
    mvaPreElID.fabsEPDiff[1]              = 99999;
    mvaPreElID.ConversionMissHits[1]      = 99999;
    mvaPreElID.dr03TkSumPt[1]             = 0.2;
    mvaPreElID.dr03EcalRecHitSumEt[1]     = 0.2;
    mvaPreElID.dr03HcalTowerSumEt[1]      = 0.2;
    mvaPreElID.numberOfLostHits[1]        = 0;
    
    looseMVAElID.cutName                  = "looseMVAElID";
    tightMVAElID.cutName                  = "tightMVAElID";

    looseElIso.cutName                    = "looseElIso";
    looseElIso.chIso04                    = 99999;
    looseElIso.nhIso04                    = 99999;
    looseElIso.phIso04                    = 99999;
    looseElIso.relCombIso04               = 0.4;
    looseElIso.chIso                    = 99999;
    looseElIso.nhIso                    = 99999;
    looseElIso.phIso                    = 99999;
    looseElIso.relCombIso               = 0.15;

    mediumElIso.cutName                   = "mediumElIso";
    mediumElIso.chIso04                   = 99999;
    mediumElIso.nhIso04                   = 99999;
    mediumElIso.phIso04                   = 99999;
    mediumElIso.relCombIso04              = 0.15;
    mediumElIso.chIso                   = 99999;
    mediumElIso.nhIso                   = 99999;
    mediumElIso.phIso                   = 99999;
    mediumElIso.relCombIso              = 0.15;

    tightElIso.cutName                   = "tightElIso";
    tightElIso.chIso04                   = 99999;
    tightElIso.nhIso04                   = 99999;
    tightElIso.phIso04                   = 99999;
    tightElIso.relCombIso04              = 0.1;
    tightElIso.chIso                   = 99999;
    tightElIso.nhIso                   = 99999;
    tightElIso.phIso                   = 99999;
    tightElIso.relCombIso              = 0.1;

    /* Muon ID */
    tightMuID.cutName                     = "tightMuID";
    tightMuID.IsPF                        = 1;
    tightMuID.IsGLB                       = 1;
    tightMuID.IsTRK                       = 1;
    tightMuID.NormalizedChi2              = 10;
    tightMuID.NumberOfValidMuonHits       = 0;
    tightMuID.NumberOfMatchedStations     = 1;
    tightMuID.NumberOfValidPixelHits      = 0;
    tightMuID.TrackLayersWithMeasurement  = 5;
    tightMuID.dxy                         = 0.2;
    tightMuID.dz                          = 0.5;

    /* Muon ISO */
    amumuMuDetIso.cutName                    = "amumuMuDetIso";
    amumuMuDetIso.hcalIso                  = 99999;
    amumuMuDetIso.ecalIso                  = 99999;
    amumuMuDetIso.trkIso                   = 0.1;
    amumuMuDetIso.relCombIso               = 0.2;

    looseMuDetIso.cutName                    = "looseMuDetIso";
    looseMuDetIso.hcalIso                  = 99999;
    looseMuDetIso.ecalIso                  = 99999;
    looseMuDetIso.trkIso                   = 99999;
    looseMuDetIso.relCombIso               = 0.2;

    tightMuDetIso.cutName                    = "tightMuDetIso";
    tightMuDetIso.hcalIso                  = 99999;
    tightMuDetIso.ecalIso                  = 99999;
    tightMuDetIso.trkIso                   = 99999;
    tightMuDetIso.relCombIso               = 0.2;

    looseMuIso.cutName                    = "looseMuIso";
    looseMuIso.chIso04                    = 99999;
    looseMuIso.nhIso04                    = 99999;
    looseMuIso.phIso04                    = 99999;
    looseMuIso.relCombIso04               = 0.2;

    tightMuIso.cutName                    = "tightMuIso";
    tightMuIso.chIso04                    = 99999;
    tightMuIso.nhIso04                    = 99999;
    tightMuIso.phIso04                    = 99999;
    tightMuIso.relCombIso04               = 0.12;

    loosePhID.cutName                     = "loosePhID";
    loosePhID.PassedEleSafeVeto[0]        = 1;
    loosePhID.HadOverEm[0]                = 0.05;
    loosePhID.sigmaIetaIeta[0]            = 0.012;

    loosePhID.PassedEleSafeVeto[1]        = 1;
    loosePhID.HadOverEm[1]                = 0.05;
    loosePhID.sigmaIetaIeta[1]            = 0.034;

    loosePhIso.cutName                    = "loosePhIso";
    loosePhIso.chIso[0]                 = 2.6;
    loosePhIso.nhIso[0]                 = 3.5;
    loosePhIso.phIso[0]                 = 1.3;

    loosePhIso.chIso[1]                 = 2.3;
    loosePhIso.nhIso[1]                 = 2.9;
    loosePhIso.phIso[1]                 = 99999;

    mediumPhID.cutName                    = "mediumPhID";
    mediumPhID.PassedEleSafeVeto[0]       = 1;
    mediumPhID.HadOverEm[0]               = 0.05;
    mediumPhID.sigmaIetaIeta[0]           = 0.011;

    mediumPhID.PassedEleSafeVeto[1]       = 1;
    mediumPhID.HadOverEm[1]               = 0.05;
    mediumPhID.sigmaIetaIeta[1]           = 0.033;

    mediumPhIso.cutName                   = "mediumPhIso";
    mediumPhIso.chIso[0]                = 1.5;
    mediumPhIso.nhIso[0]                = 1.0;
    mediumPhIso.phIso[0]                = 0.7;

    mediumPhIso.chIso[1]                = 1.2;
    mediumPhIso.nhIso[1]                = 1.5;
    mediumPhIso.phIso[1]                = 1.0;

    preSelPhID.cutName                    = "preSelPhID";
    preSelPhID.PassedEleSafeVeto[0]       = 1;
    preSelPhID.HadOverEm[0]               = 0.082;
    preSelPhID.sigmaIetaIeta[0]           = 0.014;
    preSelPhID.HcalIso[0]                 = 50;
    preSelPhID.TrkIso[0]                  = 50;
    preSelPhID.ChPfIso[0]                 = 4;

    preSelPhID.PassedEleSafeVeto[1]       = 1;
    preSelPhID.HadOverEm[1]               = 0.075;
    preSelPhID.sigmaIetaIeta[1]           = 0.034;
    preSelPhID.HcalIso[1]                 = 4;
    preSelPhID.TrkIso[1]                  = 4;
    preSelPhID.ChPfIso[1]                 = 4;

    catPhMVAID.cutName                    = "catPhMVAID";
    catPhMVAID.mvaValCat1                 = 0.126;
    catPhMVAID.mvaValCat2                 = 0.107;
    catPhMVAID.mvaValCat3                 = 0.126;
    catPhMVAID.mvaValCat4                 = 0.135;

    looseMVAPhID.cutName                  = "looseMVAPhID";
    tightMVAPhID.cutName                  = "tightMVAPhID";

    hzgMVAID.cutName                        = "hzgMVAID";
    hzgMVAID.mvaVal[0]                      = -0.9;
    hzgMVAID.mvaVal[1]                      = -0.5;
    hzgMVAID.mvaVal[2]                      = -99;
    hzgMVAID.mvaVal[3]                      = -99;
    hzgMVAID.mvaVal[4]                      = -99;
    hzgMVAID.mvaVal[5]                      = -99;
    hzgMVAID.pt[0]                          =  10;
    hzgMVAID.pt[1]                          =  20;
    hzgMVAID.eta[0]                         = -99;
    hzgMVAID.eta[1]                         = -99;
    hzgMVAID.eta[2]                         = -99;
    hzgMVAID.missHits[0]                    =  0;
    hzgMVAID.missHits[1]                    =  0;
    hzgMVAID.missHits[2]                    =  0;
    hzgMVAID.missHits[3]                    =  0;
    hzgMVAID.missHits[4]                    =  0;
    hzgMVAID.missHits[5]                    =  0;
    hzgMVAID.conversionVeto[0]              =  0;
    hzgMVAID.conversionVeto[1]              =  0;
    hzgMVAID.conversionVeto[2]              =  0;
    hzgMVAID.conversionVeto[3]              =  0;
    hzgMVAID.conversionVeto[4]              =  0;
    hzgMVAID.conversionVeto[5]              =  0;
    hzgMVAID.sip[0]                         =  0;
    hzgMVAID.sip[1]                         =  0;
    hzgMVAID.sip[2]                         =  0;
    hzgMVAID.sip[3]                         =  0;
    hzgMVAID.sip[4]                         =  0;
    hzgMVAID.sip[5]                         =  0;

    hzzMVAID.cutName                        = "hzzMVAID";
    hzzMVAID.mvaVal[0]                      =  -0.211;
    hzzMVAID.mvaVal[1]                      =  -0.396;
    hzzMVAID.mvaVal[2]                      =  -0.215;
    hzzMVAID.mvaVal[3]                      = -0.870;
    hzzMVAID.mvaVal[4]                      = -0.838;
    hzzMVAID.mvaVal[5]                      =  -0.763;
    hzzMVAID.pt[0]                          =  5.0;
    hzzMVAID.pt[1]                          =  10.0;
    hzzMVAID.eta[0]                         =  0.8;
    hzzMVAID.eta[1]                         =  1.479;
    hzzMVAID.eta[2]                         =  2.5;

    hwwMVAID.cutName                        = "hwwMVAID";
    hwwMVAID.mvaVal[0]                      =  0.00;
    hwwMVAID.mvaVal[1]                      =  0.10;
    hwwMVAID.mvaVal[2]                      =  0.62;
    hwwMVAID.mvaVal[3]                      =  0.94;
    hwwMVAID.mvaVal[4]                      =  0.85;
    hwwMVAID.mvaVal[5]                      =  0.92;
    hwwMVAID.pt[0]                          =  10;
    hwwMVAID.pt[1]                          =  20;
    hwwMVAID.eta[0]                         =  0.8;
    hwwMVAID.eta[1]                         =  1.479;
    hwwMVAID.eta[2]                         =  2.5;
    hwwMVAID.missHits[0]                    =  0;
    hwwMVAID.missHits[1]                    =  0;
    hwwMVAID.missHits[2]                    =  0;
    hwwMVAID.missHits[3]                    =  0;
    hwwMVAID.missHits[4]                    =  0;
    hwwMVAID.missHits[5]                    =  0;
    hwwMVAID.conversionVeto[0]              =  1;
    hwwMVAID.conversionVeto[1]              =  1;
    hwwMVAID.conversionVeto[2]              =  1;
    hwwMVAID.conversionVeto[3]              =  1;
    hwwMVAID.conversionVeto[4]              =  1;
    hwwMVAID.conversionVeto[5]              =  1;
    hwwMVAID.sip[0]                         = -99;
    hwwMVAID.sip[1]                         = -99;
    hwwMVAID.sip[2]                         = -99;
    hwwMVAID.sip[3]                         = -99;
    hwwMVAID.sip[4]                         = -99;
    hwwMVAID.sip[5]                         = -99;

    vbfJetID.cutName                        = "vbfJetID";
    vbfJetID.betaStarC[0]                   = 0.2;
    vbfJetID.dR2Mean[0]                     = 0.06;
    vbfJetID.betaStarC[1]                   = 0.3;
    vbfJetID.dR2Mean[1]                     = 0.05;
    vbfJetID.dR2Mean[2]                     = 0.05;
    vbfJetID.dR2Mean[3]                     = 0.055;

    looseJetID.cutName                      = "looseJetID";
    looseJetID.NHF                          = 0.99;
    looseJetID.NEMF                         = 0.99;
    looseJetID.NumConst                     = 1;
    looseJetID.MUF                          = 0.8;
    looseJetID.CHF                          = 0;
    looseJetID.CHM                          = 0;
    looseJetID.CEMF                         = 0.99;
    looseJetID.CSV                          = -99;

    bJetID.cutName                          = "bJetID";
    bJetID.NHF                              = 0.99;
    bJetID.NEMF                             = 0.99;
    bJetID.NumConst                         = 1;
    bJetID.MUF                              = 0.8;
    bJetID.CHF                              = 0;
    bJetID.CHM                              = 0;
    bJetID.CEMF                             = 0.99;
    bJetID.CSV                              = 0.898;  // medium WP
}

#include "BLT/BLTAnalysis/interface/WeightUtils.h"

WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    const std::string cmssw_base = getenv("CMSSW_BASE");

    // PU weights
    std::string puFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup_sf_2016_BCD.root";
    TFile* puFile = new TFile(puFileName.c_str(), "OPEN");
    _puReweight = (TGraph*)puFile->Get("pileup");

    // IsoMu22 efficiency ratios
    std::string triggerFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root";
    std::string filePath = "IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/";
    TFile* triggerFile = new TFile(triggerFileName.c_str(), "OPEN");
    _eff_IsoMu22_DATA[0] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu22_pass_DATA").c_str());
    _eff_IsoMu22_DATA[1] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu22_pass_DATA").c_str());
    _eff_IsoMu22_DATA[2] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu22_pass_DATA").c_str());
    _eff_IsoMu22_DATA[3] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu22_pass_DATA").c_str());

    // tight muon ID sf
    TFile* f_muRecoSF2012_ID = new TFile("../data/MuonID_Z_RunBCD_prompt80X_7p65.root", "OPEN"); 
    filePath = "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/efficienciesDATA/";
    _muSF2012_ID_DATA[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin0_&_tag_IsoMu22_pass_DATA").c_str());
    _muSF2012_ID_DATA[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin1_&_tag_IsoMu22_pass_DATA").c_str());
    _muSF2012_ID_DATA[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin2_&_tag_IsoMu22_pass_DATA").c_str());
    _muSF2012_ID_DATA[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin3_&_tag_IsoMu22_pass_DATA").c_str());

    filePath = "MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/efficienciesMC/";
    _muSF2012_ID_MC[0] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin0_&_tag_IsoMu22_pass_MC").c_str());
    _muSF2012_ID_MC[1] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin1_&_tag_IsoMu22_pass_MC").c_str());
    _muSF2012_ID_MC[2] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin2_&_tag_IsoMu22_pass_MC").c_str());
    _muSF2012_ID_MC[3] = (TGraphAsymmErrors*)f_muRecoSF2012_ID->Get((filePath + "pt_PLOT_abseta_bin3_&_tag_IsoMu22_pass_MC").c_str());
}

void WeightUtils::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}

void WeightUtils::SetDataPeriod(string dataPeriod)
{
    _dataPeriod = dataPeriod;
}

void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
}

float WeightUtils::GetPUWeight(float nPU)
{
    return _puReweight->Eval(nPU); 
}

std::pair<float,float> WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector &lepton) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(lepton.Eta()) > binningEta[i] && fabs(lepton.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }
    
    float effMC   = 1;
    float effData = 1;
    if (triggerName == "HLT_IsoMu22_v*") {
        if (lepton.Pt() < 500.) {
            effData = _eff_IsoMu22_DATA[etaBin]->Eval(lepton.Pt());
        }
    }

    return std::make_pair(effData, effMC);
}

float WeightUtils::GetMuonRecoEff(TLorentzVector& muon) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    float weight = 1;
    if (muon.Pt() < 200.) {
        weight   *= _muSF2012_ID_DATA[etaBin]->Eval(muon.Pt())/_muSF2012_ID_MC[etaBin]->Eval(muon.Pt());
        //isoWeight = 1.;//_muSF2012_ISO[etaBin]->Eval(muonP4.Pt());
    }
    
    return weight;
}

//float WeightUtils::GetElectronRecoEff(TLorentzVector& electron) const
//{
//    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
//    int etaBin = 0;
//    for (int i = 0; i < 4; ++i) {
//        if (fabs(electron.Eta()) > binningEta[i] && fabs(electron.Eta()) <= binningEta[i+1]) {
//            etaBin = i;
//            break;
//        }
//    }
//
//    float weight = 1;
//    if (electron.Pt() < 200.) {
//        weight   *= _muSF2012_ID_DATA[etaBin]->Eval(electron.Pt())/_muSF2012_ID_MC[etaBin]->Eval(electron.Pt());
//        //isoWeight = 1.;//_muSF2012_ISO[etaBin]->Eval(electronP4.Pt());
//    }
//    
//    return weight;
//}

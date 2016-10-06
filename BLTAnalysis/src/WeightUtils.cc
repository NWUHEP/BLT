#include "BLT/BLTAnalysis/interface/WeightUtils.h"

WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;

    const std::string cmssw_base = getenv("CMSSW_BASE");

    // PU weights
    std::string puFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/PUWeights_2012.root";
    TFile* puFile = new TFile(puFileName.c_str(), "OPEN");
    puReweight = (TH1D*)puFile->Get("pileup");

    // IsoMu24_eta2p1 efficiency ratios
    std::string triggerFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root";
    TFile* triggerFile = new TFile(triggerFileName.c_str(), "OPEN");
    _sf_IsoMu24_Eta2p1_data[0] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_data[1] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_data[2] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");

    _sf_IsoMu24_Eta2p1_mc[0] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_mc[1] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1_mc[2] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");

    // tight muon ID sf
    TFile* f_muRecoSF2012_ID = new TFile("../data/MuonEfficiencies_Run2012ReReco_53X.root", "OPEN"); 
    _muSF2012_ID[0] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta<0.9");
    _muSF2012_ID[1] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2");
    _muSF2012_ID[2] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1");
    _muSF2012_ID[3] = (TGraphErrors*)f_muRecoSF2012_ID->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4");

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
    return puReweight->GetBinContent(puReweight->FindBin(nPU)); 
}

std::pair<float,float> WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector &lepton) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 3; ++i) {
        if (fabs(lepton.Eta()) > binningEta[i] && fabs(lepton.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }
    
    float effMC   = 1;
    float effData = 1;
    if (triggerName == "HLT_IsoMu24_eta2p1_v*") {
        if (lepton.Pt() < 300.) {
            effMC   = _sf_IsoMu24_Eta2p1_data[etaBin]->Eval(lepton.Pt());
            effData = _sf_IsoMu24_Eta2p1_mc[etaBin]->Eval(lepton.Pt());
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
    if (muon.Pt() < 300.) {
        weight   *= _muSF2012_ID[etaBin]->Eval(muon.Pt());
        //isoWeight = 1.;//_muSF2012_ISO[etaBin]->Eval(muonP4.Pt());
    }
    
    return weight;
}

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
    _sf_IsoMu24_Eta2p1[0] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1[1] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD");
    _sf_IsoMu24_Eta2p1[2] = (TGraphErrors*)triggerFile->Get("IsoMu24_eta2p1_DATA_over_MC_TightID_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD");

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

float WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector leptonP4) const
{
    float weight = 1;
    if (triggerName == "HLT_IsoMu24_eta2p1_v*") {
        if (fabs(leptonP4.Eta()) < 0.9) {
            weight *= _sf_IsoMu24_Eta2p1[0]->Eval(leptonP4.Pt());
        } else if (fabs(leptonP4.Eta()) > 0.9 && fabs(leptonP4.Eta()) < 1.2) {
            weight *= _sf_IsoMu24_Eta2p1[1]->Eval(leptonP4.Pt());
        } else if (fabs(leptonP4.Eta()) > 1.2 && fabs(leptonP4.Eta()) < 2.4) {
            weight *= _sf_IsoMu24_Eta2p1[2]->Eval(leptonP4.Pt());
        } 
    }
    return weight;
}


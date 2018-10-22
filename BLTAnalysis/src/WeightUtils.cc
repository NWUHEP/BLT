#include "BLT/BLTAnalysis/interface/WeightUtils.h"

// helper functions
template <class Graph>
int GetBinNumber(Graph graph, float x0)
{
    Double_t x,y;
    float xmin = 1e9;
    for (int i = 0; i < graph->GetN(); i++) {
        graph->GetPoint(i, x, y);
        float diff = fabs(x - x0);
        if (diff > xmin) {
            return i - 1;
        } else {
            xmin = diff;
        }
    }
    return graph->GetN() - 1; 
}


// definitions for WeightUtils

WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;
    std::string fileName;

    const std::string cmssw_base = getenv("CMSSW_BASE");
    // PU weights
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup_sf_2016_full.root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TGraph*)puFile->Get("pileup_sf");

    // muon trigger sf (BCDEF) 
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_RunBtoF.root";
    TFile* muTriggerFile_BCDEF = new TFile(fileName.c_str(), "OPEN");

    std::string filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
    _muSF_IsoMu24_DATA_BCDEF[0] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_BCDEF[1] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_BCDEF[2] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_BCDEF[3] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/";
    _muSF_IsoMu24_MC_BCDEF[0] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_BCDEF[1] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_BCDEF[2] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_BCDEF[3] = (TGraphAsymmErrors*)muTriggerFile_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());

    // muon trigger sf (GH) 
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_Period4.root";
    TFile* muTriggerFile_GH = new TFile(fileName.c_str(), "OPEN");

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
    _muSF_IsoMu24_DATA_GH[0] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[1] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[2] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _muSF_IsoMu24_DATA_GH[3] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());

    filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesMC/";
    _muSF_IsoMu24_MC_GH[0] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[1] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[2] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());
    _muSF_IsoMu24_MC_GH[3] = (TGraphAsymmErrors*)muTriggerFile_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_MC").c_str());

    // muon tight ID sf (BCDEF)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF_ID_BCDEF = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // muon tight ID sf (GH)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_ID_GH = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // tight muon ISO sf (BCDEF)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF_ISO_BCDEF = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
    _muSF_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

    filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
    _muSF_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());

    // tight muon ISO sf (GH)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_ISO_GH = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "TightISO_TightID_pt_eta/efficienciesDATA/";
    _muSF_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_DATA").c_str());
    _muSF_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_DATA").c_str());

    filePath = "TightISO_TightID_pt_eta/efficienciesMC/";
    _muSF_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_MC").c_str());
    _muSF_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_MC").c_str());

    // electron trigger efficiencies (BCDEF and GH)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/2017Feb-A-Popov_TriggerSF_Run2016BCDEF_v1.root";
    TFile* elTriggerFile_BCDEF = new TFile(fileName.c_str(), "OPEN");
    _elSF_Trigger_BCDEF = (TH2D*)elTriggerFile_BCDEF->Get("Ele27_WPTight_Gsf");

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/2017Feb-A-Popov_TriggerSF_Run2016GH_v1.root";
    TFile* elTriggerFile_GH = new TFile(fileName.c_str(), "OPEN");
    _elSF_Trigger_GH = (TH2D*)elTriggerFile_GH->Get("Ele27_WPTight_Gsf");

    // electron reco efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/EGM2D_BtoH_low_RecoSF_Legacy2016.root";
    TFile* f_eleRecoSF_1 = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO[0] = (TGraphErrors*)f_eleRecoSF_1->Get("grSF1D_0");

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root";
    TFile* f_eleRecoSF_2 = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO[1] = (TGraphErrors*)f_eleRecoSF_2->Get("grSF1D_0");
    _eleSF_RECO[2] = (TGraphErrors*)f_eleRecoSF_2->Get("grSF1D_1");
    _eleSF_RECO[3] = (TGraphErrors*)f_eleRecoSF_2->Get("grSF1D_2");
    _eleSF_RECO[4] = (TGraphErrors*)f_eleRecoSF_2->Get("grSF1D_3");

    // electron id efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_ID_2016.root";
    TFile* f_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
    _eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
    _eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
    _eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
    _eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");

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

EfficiencyContainer WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector &lepton) const
{
    float effData = 1;
    float errData = 1;
    float effMC   = 1;
    float errMC   = 1;
    if (triggerName == "HLT_IsoMu24_v*") {
        float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
        int etaBin = 0;
        for (int i = 0; i < 4; ++i) {
            if (fabs(lepton.Eta()) > binningEta[i] && fabs(lepton.Eta()) <= binningEta[i+1]) {
                etaBin = i;
                break;
            }
        }
        
        if (_dataPeriod == "2016BtoF") {
            effData = _muSF_IsoMu24_DATA_BCDEF[etaBin]->Eval(lepton.Pt());
            effMC   = _muSF_IsoMu24_MC_BCDEF[etaBin]->Eval(lepton.Pt());

            int ptBin = GetBinNumber(_muSF_IsoMu24_DATA_BCDEF[etaBin], lepton.Pt()); 
            errData = _muSF_IsoMu24_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_IsoMu24_MC_BCDEF[etaBin]->GetErrorY(ptBin);
        } else if (_dataPeriod == "2016GH") {
            effData = _muSF_IsoMu24_DATA_GH[etaBin]->Eval(lepton.Pt());
            effMC   = _muSF_IsoMu24_MC_GH[etaBin]->Eval(lepton.Pt());

            int ptBin = GetBinNumber(_muSF_IsoMu24_DATA_GH[etaBin], lepton.Pt()); 
            errData = _muSF_IsoMu24_DATA_GH[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_IsoMu24_MC_GH[etaBin]->GetErrorY(ptBin);
        }
    } else if (triggerName == "HLT_Ele27_WPTight_Gsf_v*") {
        if (_dataPeriod == "2016BtoF") {
            int bin = _elSF_Trigger_BCDEF->FindBin(lepton.Pt(), lepton.Eta());
            effData = _elSF_Trigger_BCDEF->GetBinContent(bin);
            effMC = 1.;
            errData = _elSF_Trigger_BCDEF->GetBinError(bin);
            errMC   = 0.;
        } else if (_dataPeriod == "2016GH") {
            int bin = _elSF_Trigger_GH->FindBin(lepton.Pt(), lepton.Eta());
            effData = _elSF_Trigger_GH->GetBinContent(bin);
            effMC = 1.;
            errData = _elSF_Trigger_GH->GetBinError(bin);
            errMC   = 0.;
        }
    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetMuonRecoEff(TLorentzVector& muon) const
{
    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.Eta()) > binningEta[i] && fabs(muon.Eta()) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    float effData = 1;
    float errData = 1;
    float effMC   = 1;
    float errMC   = 1;
    if (_dataPeriod == "2016BtoF") {
        float effDataID  = _muSF_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        float effDataISO = _muSF_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effData    = effDataID*effDataISO;

        float effMCID    = _muSF_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
        float effMCISO   = _muSF_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
        effMC      = effMCID*effMCISO;

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_BCDEF[etaBin], muon.Pt()); 
        float errDataID  = _muSF_ID_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        float errDataISO = _muSF_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        float errMCID    = _muSF_ID_MC_BCDEF[etaBin]->GetErrorY(ptBin);
        float errMCISO   = _muSF_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);

        errData = effData*(pow(errDataID/effDataID, 2) + pow(errDataISO/effDataISO, 2));
        errMC   = effMC*(pow(errMCID/effMCID, 2) + pow(errMCISO/effMCISO, 2));

    } else if (_dataPeriod == "2016GH") {
        float effDataID  = _muSF_ID_DATA_GH[etaBin]->Eval(muon.Pt());
        float effDataISO = _muSF_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effData    = effDataID*effDataISO;

        float effMCID    = _muSF_ID_MC_GH[etaBin]->Eval(muon.Pt());
        float effMCISO   = _muSF_ISO_MC_GH[etaBin]->Eval(muon.Pt());
        effMC      = effMCID*effMCISO;

        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_GH[etaBin], muon.Pt()); 
        float errDataID  = _muSF_ID_DATA_GH[etaBin]->GetErrorY(ptBin);
        float errDataISO = _muSF_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        float errMCID    = _muSF_ID_MC_GH[etaBin]->GetErrorY(ptBin);
        float errMCISO   = _muSF_ISO_MC_GH[etaBin]->GetErrorY(ptBin);

        errData = effData*(pow(errDataID/effDataID, 2) + pow(errDataISO/effDataISO, 2));
        errMC   = effMC*(pow(errMCID/effMCID, 2) + pow(errMCISO/effMCISO, 2));
    }
    
    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetElectronRecoEff(TLorentzVector& electron) const
{
    float binningPtID[] = {10., 20., 35., 50., 90., 999.};
    float binningPtRECO[] = {10., 20., 45., 75., 100., 999.};
    int ptBinID = 0;
    int ptBinRECO = 0;
    for (int i = 0; i < 5; ++i) {
        if (fabs(electron.Pt()) > binningPtID[i] && fabs(electron.Pt()) <= binningPtID[i+1]) {
            ptBinID = i;
            break;
        }
        if (fabs(electron.Pt()) > binningPtRECO[i] && fabs(electron.Pt()) <= binningPtID[i+1]) {
            ptBinRECO = i;
            break;
        }
    }
    float sfId   = _eleSF_ID[ptBinID]->Eval(electron.Eta());
    float sfReco = _eleSF_RECO[ptBinRECO]->Eval(electron.Eta());

    int etaBin;
    etaBin = GetBinNumber<TGraphErrors*>(_eleSF_RECO[ptBinRECO], electron.Eta()); 
    float errReco = _eleSF_RECO[ptBinRECO]->GetErrorY(etaBin);

    etaBin = GetBinNumber<TGraphErrors*>(_eleSF_ID[ptBinID], electron.Eta()); 
    float errId  = _eleSF_ID[ptBinID]->GetErrorY(etaBin);
    float sfErr = sfReco*sfId*(pow(errReco/sfReco, 2) + pow(errId/sfId, 2));
    
    EfficiencyContainer effCont(sfReco*sfId, 1., sfErr, 0.);
    return effCont;
}

//
// Definitions for EfficiencyContainer
//

EfficiencyContainer::EfficiencyContainer()
{
    _dataEff = 0.;
    _mcEff   = 0.;
    _dataErr = 0.;
    _mcErr   = 0.;
};

EfficiencyContainer::EfficiencyContainer(float dataEff, float mcEff, float dataErr, float mcErr)
{
    _dataEff = dataEff;
    _mcEff   = mcEff;
    _dataErr = dataErr;
    _mcErr   = mcErr;
};

void EfficiencyContainer::SetData(float dataEff, float mcEff, float dataErr, float mcErr)
{
    _dataEff = dataEff;
    _mcEff   = mcEff;
    _dataErr = dataErr;
    _mcErr   = mcErr;
};


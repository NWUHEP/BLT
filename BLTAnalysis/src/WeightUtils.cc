#include "BLT/BLTAnalysis/interface/WeightUtils.h"

// helper functions
template <class Graph>
int GetBinNumber(Graph graph, float x0)
{
    Double_t x,y;
    for (int i = 0; i < graph->GetN(); i++) {
        graph->GetPoint(i, x, y);
        float diff = fabs(x - x0);
        float ex = graph->GetErrorX(i);
        //cout << i << ", " << x << ", " << y << ", " << ex << ", " << diff << endl;
        if (diff < ex) {
            return i;
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
    if (_dataPeriod == "2016"){
      fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup/pileup_sf_2016_full.root";
      TFile* puFile = new TFile(fileName.c_str(), "OPEN");
      _puReweight2016 = (TGraph*)puFile->Get("pileup_sf");
    }
    
    // PU weights
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup/pu_weights_" + _dataPeriod + ".root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TH1D*)puFile->Get("mcwei_run000001");


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

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/single_electron_sf_BCDEF.root";
    TFile* elTriggerFile_BCDEF = new TFile(fileName.c_str(), "OPEN");
    _elSF_Trigger_BCDEF        = (TH2D*)elTriggerFile_BCDEF->Get("h2_sf_BCDEF");
    _elSF_Trigger_BCDEF_tag    = (TH2D*)elTriggerFile_BCDEF->Get("h2_err_tag_BCDEF");
    _elSF_Trigger_BCDEF_probe  = (TH2D*)elTriggerFile_BCDEF->Get("h2_err_prob_BCDEF");

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/single_electron_sf_GH.root";
    TFile* elTriggerFile_GH = new TFile(fileName.c_str(), "OPEN");
    _elSF_Trigger_GH        = (TH2D*)elTriggerFile_GH->Get("h2_sf_GH");
    _elSF_Trigger_GH_tag    = (TH2D*)elTriggerFile_GH->Get("h2_err_tag_GH");
    _elSF_Trigger_GH_probe  = (TH2D*)elTriggerFile_GH->Get("h2_err_prob_GH");

    // fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/TriggerSF_Run2016BCDEF_v1.root";
    // TFile* elTriggerFile_BCDEF = new TFile(fileName.c_str(), "OPEN");
    // _elSF_Trigger_BCDEF = (TH2D*)elTriggerFile_BCDEF->Get("Ele27_WPTight_Gsf");

    // fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/TriggerSF_Run2016GH_v1.root";
    // TFile* elTriggerFile_GH = new TFile(fileName.c_str(), "OPEN");
    // _elSF_Trigger_GH = (TH2D*)elTriggerFile_GH->Get("Ele27_WPTight_Gsf");


    // electron reco efficiencies 
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_reco/egamma_eff_reco_2016.root";
    TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");
    _eleSF_RECO->Sort();
    

    // electron id efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/egamma_eff_ID_2016.root";
    TFile* f_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
    _eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
    _eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
    _eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
    _eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");
    // sort the graph points for easier access
    _eleSF_ID[0]->Sort();
    _eleSF_ID[1]->Sort();
    _eleSF_ID[2]->Sort();
    _eleSF_ID[3]->Sort();
    _eleSF_ID[4]->Sort();


    // electron l1 prefiring
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_prefiring/L1prefiring_jetpt_2016BtoH.root";
    TFile* f_elePrefiringSF_jet = new TFile(fileName.c_str(), "OPEN");
    _elSF_Prefiring_jet = (TH2D*)f_elePrefiringSF_jet->Get("L1prefiring_jetpt_2016BtoH");

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_prefiring/L1prefiring_photonpt_2016BtoH.root";
    TFile* f_elePrefiringSF_photon = new TFile(fileName.c_str(), "OPEN");
    _elSF_Prefiring_photon = (TH2D*)f_elePrefiringSF_photon->Get("L1prefiring_photonpt_2016BtoH");

    // photon mva id (90%) efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_id/photon_mva_id_2016.root";
    TFile* f_mva_gammaIdSF = new TFile(fileName.c_str(), "OPEN");
    _mva_gammaSF_ID[0] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_0");
    _mva_gammaSF_ID[1] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_1");
    _mva_gammaSF_ID[2] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_2");
    _mva_gammaSF_ID[3] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_3");
    _mva_gammaSF = (TH2F *)f_mva_gammaIdSF->Get("EGamma_SF2D");

    // photon r9 reweighting
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_iso/photon_r9_reweighting_2016.root";
    TFile* f_photon_r9 = new TFile(fileName.c_str(), "OPEN"); 
    _photon_r9_barrel = (TGraph *)f_photon_r9->Get("transffull5x5R9EB");
    _photon_r9_endcap = (TGraph *)f_photon_r9->Get("transffull5x5R9EE");


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
    if(_dataPeriod=="2016") return _puReweight2016->Eval(nPU); 
    else return _puReweight->GetBinContent(_puReweight->FindBin(nPU)); 
}

EfficiencyContainer WeightUtils::GetTriggerEffWeight(string triggerName, TLorentzVector &lepton) const
{
    float effData = 1;
    float errData = 0;
    float effMC   = 1;
    float errMC   = 0;
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
        // official numbers
        if (_dataPeriod == "2016BtoF") {
            int bin = _elSF_Trigger_BCDEF->FindBin(lepton.Pt(), lepton.Eta());
            effData = _elSF_Trigger_BCDEF->GetBinContent(bin);
            errData = _elSF_Trigger_BCDEF->GetBinError(bin);
        } else if (_dataPeriod == "2016GH") {
            int bin = _elSF_Trigger_GH->FindBin(lepton.Pt(), lepton.Eta());
            effData = _elSF_Trigger_GH->GetBinContent(bin);
            errData = _elSF_Trigger_GH->GetBinError(bin);
        }
        effMC   = 1.;
        errMC   = 0.;

        // values from Kevin, et. al.
        //int etaBin = 0;
        //int ptBin = 0;
        //for (int i = 0; i < 13; ++i) {
        //    if (lepton.Eta() > _eleEtaBins[i] && lepton.Eta() <= _eleEtaBins[i+1]) {
        //        etaBin = i;
        //        break;
        //    }
        //}
        //for (int i = 0; i < 8; ++i) {
        //    if (lepton.Pt() > _elePtBins[i] && lepton.Pt() <= _elePtBins[i+1]) {
        //        ptBin = i;
        //        break;
        //    }
        //}
        //effData = _ele_trigEff_data[etaBin][ptBin];
        //effMC   = _ele_trigEff_mc[etaBin][ptBin];
        //errData = 0.005*_ele_trigEff_data[etaBin][ptBin];
        //errMC   = 0.005*_ele_trigEff_mc[etaBin][ptBin];

    }

    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

float WeightUtils::GetEleTriggerSyst(string errType, TLorentzVector &lepton) const
 {
     float err = 0.;
     if (_dataPeriod == "2016BtoF") {
         if (errType == "tag") {
             int bin = _elSF_Trigger_BCDEF_tag->FindBin(lepton.Pt(), lepton.Eta());
             err = _elSF_Trigger_BCDEF_tag->GetBinContent(bin);
         } else if (errType == "probe") {
             int bin = _elSF_Trigger_BCDEF_probe->FindBin(lepton.Pt(), lepton.Eta());
             err = _elSF_Trigger_BCDEF_probe->GetBinContent(bin);
         }
     } else if (_dataPeriod == "2016GH") {
         if (errType == "tag") {
             int bin = _elSF_Trigger_GH_tag->FindBin(lepton.Pt(), lepton.Eta());
             err = _elSF_Trigger_GH_tag->GetBinContent(bin);
         } else if (errType == "probe") {
             int bin = _elSF_Trigger_GH_probe->FindBin(lepton.Pt(), lepton.Eta());
             err = _elSF_Trigger_GH_probe->GetBinContent(bin);
         }
     }

     return err;
 }


EfficiencyContainer WeightUtils::GetMuonIDEff(TLorentzVector& muon) const
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
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_BCDEF[etaBin], muon.Pt()); 
        effData   = _muSF_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ID_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ID_MC_BCDEF[etaBin]->GetErrorY(ptBin);

    } else if (_dataPeriod == "2016GH") {
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_GH[etaBin], muon.Pt()); 
        effData   = _muSF_ID_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ID_MC_GH[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ID_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ID_MC_GH[etaBin]->GetErrorY(ptBin);
    }
    
    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetMuonISOEff(TLorentzVector& muon) const
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
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ISO_DATA_BCDEF[etaBin], muon.Pt()); 
        effData   = _muSF_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);

    } else if (_dataPeriod == "2016GH") {
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ISO_DATA_GH[etaBin], muon.Pt()); 
        effData   = _muSF_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ISO_MC_GH[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ISO_MC_GH[etaBin]->GetErrorY(ptBin);
    }
    
    EfficiencyContainer effCont(effData, effMC, errData, errMC);
    return effCont;
}

EfficiencyContainer WeightUtils::GetElectronRecoEff(TLorentzVector& electron) const
{
    int etaBin = GetBinNumber<TGraphErrors*>(_eleSF_RECO, electron.Eta());
    float sf   = _eleSF_RECO->Eval(electron.Eta());
    float err  = _eleSF_RECO->GetErrorY(etaBin);
    
    EfficiencyContainer effCont(sf, 1., err, 0.);
    return effCont;
}

EfficiencyContainer WeightUtils::GetElectronIDEff(TLorentzVector& electron) const
{
    float binningPt[] = {10., 20., 35., 50., 90., 9999.};
    int ptBin = 0;
    for (int i = 0; i < 5; ++i) {
        if (fabs(electron.Pt()) > binningPt[i] && fabs(electron.Pt()) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }

    int etaBin = GetBinNumber<TGraphErrors*>(_eleSF_ID[ptBin], electron.Eta());
    float sf   = _eleSF_ID[ptBin]->Eval(electron.Eta());
    float err  = _eleSF_ID[ptBin]->GetErrorY(etaBin);

    //cout << electron.Eta() << ", " << etaBin << ", " << sf << ", " << err << endl;
    
    EfficiencyContainer effCont(sf, 1., err, 0.);
    return effCont;
}

EfficiencyContainer WeightUtils::GetElectronPrefiringWeight(vector<TLorentzVector> prefiring_photons, vector<TLorentzVector> prefiring_jets) const
//https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1ECALPrefiringWeightProducer.cc
{

    //Probability for the event NOT to prefire, computed with the prefiring maps per object.
    //Up and down values correspond to the resulting value when shifting up/down all prefiring rates in prefiring maps.

    double nonPrefiringProba[3] = {1., 1., 1.};  
    //0: central, 1: up, 2: down
    for (int fluct = 0; fluct<3; fluct++) {

        for (const auto& photon : prefiring_photons) {
            nonPrefiringProba[fluct] *= (1. - GetPrefiringRate(photon.Eta(), photon.Pt(), _elSF_Prefiring_photon, fluct));
        }

        //Now applying the prefiring maps to jets in the affected regions.
        for (const auto& jet : prefiring_jets) {

            

            //Loop over photons to remove overlap
            double nonprefiringprobfromoverlappingphotons = 1.;
            for (const auto& photon : prefiring_photons) {
                double dR = jet.DeltaR(photon);
                if (dR > 0.4)
                    continue;
                nonprefiringprobfromoverlappingphotons *= (1. - GetPrefiringRate(photon.Eta(), photon.Pt(), _elSF_Prefiring_photon, fluct));
            }

            double nonprefiringprobfromoverlappingjet = 1. - GetPrefiringRate(jet.Eta(), jet.Pt(), _elSF_Prefiring_jet, fluct);


            // if not overlap
            if (nonprefiringprobfromoverlappingphotons == 1.) {
                nonPrefiringProba[fluct] *= nonprefiringprobfromoverlappingjet;
            }
            //If overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
            else if (nonprefiringprobfromoverlappingphotons > nonprefiringprobfromoverlappingjet) {
                if (nonprefiringprobfromoverlappingphotons != 0.) {
                    nonPrefiringProba[fluct] *= nonprefiringprobfromoverlappingjet / nonprefiringprobfromoverlappingphotons;
                } else {
                    nonPrefiringProba[fluct] = 0.;
                }
            }
            //Last case: if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight, and do nothing.
        }
    }

    float sf   = nonPrefiringProba[0];
    float err  = abs(nonPrefiringProba[1]-nonPrefiringProba[2])/2;
    EfficiencyContainer effCont(sf, 1., err, 0.);
    // cout << sf << "," << err << endl; // about ~ 0.985775,0.00459289
    return effCont;

}

double WeightUtils::GetPrefiringRate(double eta, double pt, TH2D* h_prefmap, int fluctuation) const {
    //same as https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatUtils/plugins/L1ECALPrefiringWeightProducer.cc#L185

    //Check pt is not above map overflow
    int nbinsy = h_prefmap->GetNbinsY();
    double maxy = h_prefmap->GetYaxis()->GetBinLowEdge(nbinsy + 1);
    if (pt >= maxy)
        pt = maxy - 0.01;
    int thebin = h_prefmap->FindBin(eta, pt);


    double prefrate = h_prefmap->GetBinContent(thebin);

    double statuncty = h_prefmap->GetBinError(thebin);
    double systuncty = prefiringRateSystUnc_ * prefrate;

    if (fluctuation == 1) // up
        prefrate = std::min(1., prefrate + sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
    if (fluctuation == 2) // down
        prefrate = std::max(0., prefrate - sqrt(pow(statuncty, 2) + pow(systuncty, 2)));
    return prefrate;
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

float WeightUtils::GetPhotonMVAIdEff(TPhoton& photon) const
{
    /*float binningPt[] = {20., 35., 50., 90., 150.};
    int ptBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(photon.calibPt) > binningPt[i] && fabs(photon.calibPt) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }*/
    float tmpPhotonPt = (photon.pt < 150.) ? photon.pt : 149.;
    float weight = 1.;
    //if (photon.calibPt < 150.) {
    //    weight *= _mva_gammaSF_ID[ptBin]->Eval(photon.scEta);
    //}
    weight *= _mva_gammaSF->GetBinContent(_mva_gammaSF->FindBin(photon.scEta, tmpPhotonPt));

    // electron veto scale factor
    if (fabs(photon.scEta) <= 1.49)
        weight *= 0.9938;
    else if (fabs(photon.scEta) > 1.49)
        weight *= 0.9875;

    if (weight == 0)
        weight = 1.;
    
    return weight;
}

float WeightUtils::GetCorrectedPhotonR9(TPhoton& photon) const 
{
    float r9 = photon.r9;
    if (fabs(photon.scEta) < 1.444)
        r9 = _photon_r9_barrel->Eval(photon.r9);
    else if (fabs(photon.scEta) > 1.566)
        r9 = _photon_r9_endcap->Eval(photon.r9);
    else
        std::cout << "bad value of photon scEta: returning original r9" << std::endl;
    
    return r9;
}

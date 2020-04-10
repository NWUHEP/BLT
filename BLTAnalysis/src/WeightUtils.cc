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
    //----------------------------------
    // PU weights
    //----------------------------------
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup/pu_weights_" + _dataPeriod + ".root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TH1D*)puFile->Get("mcwei_run000001");



    //----------------------------------
    // muon trigger sf 
    //----------------------------------
    // 2016 (BCDEF) 
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
    // 2016 (GH) 
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
    // 2017
    // fix me for IsoMu24
    // only IsoMu27 are provided for 2017
    // https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root";
    TFile* muTriggerFile_2017 = new TFile(fileName.c_str(), "OPEN");
    _muSF_Trigger_2017 = (TH2D*) muTriggerFile_2017->Get("IsoMu27_PtEtaBins/pt_abseta_ratio");
    // 2018
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root";
    TFile* muTriggerFile_2018 = new TFile(fileName.c_str(), "OPEN");
    _muSF_Trigger_2018 = (TH2D*) muTriggerFile_2018->Get("IsoMu24_PtEtaBins/pt_abseta_ratio");
    // 2018 after Muon HLT Upgrade
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root";
    muTriggerFile_2018 = new TFile(fileName.c_str(), "OPEN");
    _muSF_Trigger_2018_AfterMuonHLTUpdated = (TH2D*) muTriggerFile_2018->Get("IsoMu24_PtEtaBins/pt_abseta_ratio");



    
    //----------------------------------
    // muon tight ID sf 
    //----------------------------------
    // 2016(BCDEF)
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
    // 2016(GH)
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
    // 2017,2018
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/Run"+_dataPeriod+"_SF_ID.root";
    TFile* f_muRecoSF_ID = new TFile(fileName.c_str(), "OPEN"); 
    if (_dataPeriod == "2017"){
        _muSF2D_ID = (TH2D*)f_muRecoSF_ID->Get("NUM_TightID_DEN_genTracks_pt_abseta");
    } else if (_dataPeriod == "2018"){
        _muSF2D_ID = (TH2D*)f_muRecoSF_ID->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta");   
    }
    
        

    //----------------------------------
    // tight muon ISO sf 
    //----------------------------------
    // 2016(BCDEF)
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
    // 2016(GH)
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
    // 2017,2018
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/Run"+_dataPeriod+"_SF_ISO.root";
    TFile* f_muRecoSF_ISO = new TFile(fileName.c_str(), "OPEN");
    _muSF2D_ISO = (TH2D*)f_muRecoSF_ISO->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta");

    //----------------------------------
    // electron trigger efficiencies (BCDEF and GH)
    //----------------------------------
    // 2016 BCDEF
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/single_electron_sf_BCDEF.root";
    TFile* elTriggerFile_BCDEF = new TFile(fileName.c_str(), "OPEN");
    _elSF_Trigger_BCDEF        = (TH2D*)elTriggerFile_BCDEF->Get("h2_sf_BCDEF");
    _elSF_Trigger_BCDEF_tag    = (TH2D*)elTriggerFile_BCDEF->Get("h2_err_tag_BCDEF");
    _elSF_Trigger_BCDEF_probe  = (TH2D*)elTriggerFile_BCDEF->Get("h2_err_prob_BCDEF");
    // 2016 GH
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_trigger/single_electron_sf_GH.root";
    TFile* elTriggerFile_GH = new TFile(fileName.c_str(), "OPEN");
    _elSF_Trigger_GH        = (TH2D*)elTriggerFile_GH->Get("h2_sf_GH");
    _elSF_Trigger_GH_tag    = (TH2D*)elTriggerFile_GH->Get("h2_err_tag_GH");
    _elSF_Trigger_GH_probe  = (TH2D*)elTriggerFile_GH->Get("h2_err_prob_GH");

    //----------------------------------
    // electron reco efficiencies 
    //----------------------------------
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_reco/egamma_eff_reco_"+_dataPeriod+".root";
    TFile* elRecoFile = new TFile(fileName.c_str(), "OPEN"); 
    _elSF2D_RECO = (TH2D*)elRecoFile->Get("EGamma_SF2D");

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_reco/egamma_eff_reco_2017_ls20.root";
    elRecoFile = new TFile(fileName.c_str(), "OPEN"); 
    _elSF2D_RECO_2017ls20 = (TH2D*)elRecoFile->Get("EGamma_SF2D");

    

    //----------------------------------
    // electron id efficiencies
    //----------------------------------
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/egamma_eff_ID_"+_dataPeriod+".root";
    TFile* elIDFile = new TFile(fileName.c_str(), "OPEN"); 
    _elSF2D_ID = (TH2D*)elIDFile->Get("EGamma_SF2D");
     
    //----------------------------------
    // photon mva id (90%) efficiencies
    //----------------------------------
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_id/photon_mva_id_2016.root";
    TFile* f_mva_gammaIdSF = new TFile(fileName.c_str(), "OPEN");
    _mva_gammaSF_ID[0] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_0");
    _mva_gammaSF_ID[1] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_1");
    _mva_gammaSF_ID[2] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_2");
    _mva_gammaSF_ID[3] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_3");
    _mva_gammaSF = (TH2F *)f_mva_gammaIdSF->Get("EGamma_SF2D");

    //----------------------------------
    // photon r9 reweighting
    //----------------------------------
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

    EfficiencyContainer effCont(1, 1, 0, 0);

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
            effCont.SetData(effData, effMC, errData, errMC);
        } else if (_dataPeriod == "2016GH") {
            effData = _muSF_IsoMu24_DATA_GH[etaBin]->Eval(lepton.Pt());
            effMC   = _muSF_IsoMu24_MC_GH[etaBin]->Eval(lepton.Pt());

            int ptBin = GetBinNumber(_muSF_IsoMu24_DATA_GH[etaBin], lepton.Pt()); 
            errData = _muSF_IsoMu24_DATA_GH[etaBin]->GetErrorY(ptBin);
            errMC   = _muSF_IsoMu24_MC_GH[etaBin]->GetErrorY(ptBin);
            effCont.SetData(effData, effMC, errData, errMC);
        } else if(_dataPeriod == "2017") {
            int bin = _muSF_Trigger_2017->FindBin(lepton.Pt(),abs(lepton.Eta()));
            float sf = _muSF_Trigger_2017->GetBinContent(bin);
            float err = _muSF_Trigger_2017->GetBinError(bin);
            if (_muSF_Trigger_2017->IsBinUnderflow(bin) || _muSF_Trigger_2017->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
            effCont.SetData(sf, 1., err, 0.);
            effCont.SetData(1., 1., 0., 0.);
            // if (lepton.Pt()>120 || lepton.Pt()<20) cout<< _muSF_Trigger_2017->IsBinUnderflow(bin)<< _muSF_Trigger_2017->IsBinOverflow(bin) <<"Muon Trigger: pt=" << lepton.Pt() << " GeV, sf=" << sf <<endl;
        } else if(_dataPeriod == "2018") {
            int bin = _muSF_Trigger_2018->FindBin(lepton.Pt(),abs(lepton.Eta()));
            float sf = _muSF_Trigger_2018->GetBinContent(bin);
            float err = _muSF_Trigger_2018->GetBinError(bin);
            if (_muSF_Trigger_2018->IsBinUnderflow(bin) || _muSF_Trigger_2018->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
            effCont.SetData(sf, 1., err, 0.);
            // if (lepton.Pt()>120 || lepton.Pt()<20) cout<< _muSF_Trigger_2018->IsBinUnderflow(bin)<< _muSF_Trigger_2018->IsBinOverflow(bin) << "Muon Trigger: pt=" << lepton.Pt() << " GeV, sf=" << sf <<endl;
        }

    } else if (triggerName == "HLT_Ele27_WPTight_Gsf_v*") {
        // official numbers
        if (_dataPeriod == "2016BtoF") {
            int bin = _elSF_Trigger_BCDEF->FindBin(lepton.Pt(), lepton.Eta());
            effData = _elSF_Trigger_BCDEF->GetBinContent(bin);
            errData = _elSF_Trigger_BCDEF->GetBinError(bin);
            effMC   = 1.;
            errMC   = 0.;
            effCont.SetData(effData, effMC, errData, errMC);

        } else if (_dataPeriod == "2016GH") {
            int bin = _elSF_Trigger_GH->FindBin(lepton.Pt(), lepton.Eta());
            effData = _elSF_Trigger_GH->GetBinContent(bin);
            errData = _elSF_Trigger_GH->GetBinError(bin);
            effMC   = 1.;
            errMC   = 0.;
            effCont.SetData(effData, effMC, errData, errMC);

        } else if (_dataPeriod == "2017" || _dataPeriod == "2018" ){
            effData = 1;
            errData = 0;
            effMC   = 1.;
            errMC   = 0.;
            effCont.SetData(effData, effMC, errData, errMC);
        }

        
    } 

    
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
     } else if(_dataPeriod == "2017" || _dataPeriod == "2018" ){
         err = 0;
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

    EfficiencyContainer effCont(1, 1, 0, 0);

    if (_dataPeriod == "2016BtoF") {
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_BCDEF[etaBin], muon.Pt()); 
        effData   = _muSF_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ID_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ID_MC_BCDEF[etaBin]->GetErrorY(ptBin);
        effCont.SetData(effData, effMC, errData, errMC);
        
    } else if (_dataPeriod == "2016GH") {
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ID_DATA_GH[etaBin], muon.Pt()); 
        effData   = _muSF_ID_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ID_MC_GH[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ID_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ID_MC_GH[etaBin]->GetErrorY(ptBin);
        effCont.SetData(effData, effMC, errData, errMC);

    } else if (_dataPeriod == "2017" || _dataPeriod == "2018" ){
        int bin = _muSF2D_ID->FindBin(muon.Pt(),abs(muon.Eta()));
        float sf = _muSF2D_ID->GetBinContent(bin);
        float err = _muSF2D_ID->GetBinError(bin);
        if (_muSF2D_ID->IsBinUnderflow(bin) || _muSF2D_ID->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
        effCont.SetData(sf, 1., err, 0.);
        // if (muon.Pt()>120 || muon.Pt()<20) cout<< _muSF2D_ID->IsBinUnderflow(bin)<< _muSF2D_ID->IsBinOverflow(bin) <<"Muon ID: pt=" << muon.Pt() << " GeV, sf=" << sf <<endl;
    }
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

    EfficiencyContainer effCont(1, 1, 0, 0);

    if (_dataPeriod == "2016BtoF") {
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ISO_DATA_BCDEF[etaBin], muon.Pt()); 
        effData   = _muSF_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ISO_DATA_BCDEF[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ISO_MC_BCDEF[etaBin]->GetErrorY(ptBin);
        effCont.SetData(effData, effMC, errData, errMC);

    } else if (_dataPeriod == "2016GH") {
        int ptBin = GetBinNumber<TGraphAsymmErrors*>(_muSF_ISO_DATA_GH[etaBin], muon.Pt()); 
        effData   = _muSF_ISO_DATA_GH[etaBin]->Eval(muon.Pt());
        effMC     = _muSF_ISO_MC_GH[etaBin]->Eval(muon.Pt());
        errData   = _muSF_ISO_DATA_GH[etaBin]->GetErrorY(ptBin);
        errMC     = _muSF_ISO_MC_GH[etaBin]->GetErrorY(ptBin);
        effCont.SetData(effData, effMC, errData, errMC);

    } else if (_dataPeriod == "2017" || _dataPeriod == "2018" ){
        int bin = _muSF2D_ISO->FindBin(muon.Pt(),abs(muon.Eta()));
        float sf = _muSF2D_ISO->GetBinContent(bin);
        float err = _muSF2D_ISO->GetBinError(bin);
        if (_muSF2D_ISO->IsBinUnderflow(bin) || _muSF2D_ISO->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
        effCont.SetData(sf, 1., err, 0.);
        // if (muon.Pt()>120 || muon.Pt()<20) cout<< _muSF2D_ISO->IsBinUnderflow(bin)<< _muSF2D_ISO->IsBinOverflow(bin) <<"Muon Iso: pt=" << muon.Pt() << " GeV, sf=" << sf <<endl;
    }
    
    return effCont;
}

EfficiencyContainer WeightUtils::GetElectronRecoEff(TLorentzVector& electron) const
{   
    EfficiencyContainer effCont(1, 1., 0, 0.);

    if (_dataPeriod == "2017" && electron.Pt()<20){
        int bin = _elSF2D_RECO_2017ls20->FindBin(electron.Eta(),electron.Pt());
        float sf = _elSF2D_RECO_2017ls20->GetBinContent(bin);
        float err = _elSF2D_RECO_2017ls20->GetBinError(bin);
        if (_elSF2D_RECO_2017ls20->IsBinUnderflow(bin) || _elSF2D_RECO_2017ls20->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
        effCont.SetData(sf, 1., err, 0.);
        // if (electron.Pt()>120 || electron.Pt()<20) cout<< "Electron Reco: pt=" << electron.Pt() << " GeV, sf=" << sf <<endl;
        
    } else{
        int bin = _elSF2D_RECO->FindBin(electron.Eta(),electron.Pt());
        float sf = _elSF2D_RECO->GetBinContent(bin);
        float err = _elSF2D_RECO->GetBinError(bin);
        if (_elSF2D_RECO->IsBinUnderflow(bin) || _elSF2D_RECO->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
        effCont.SetData(sf, 1., err, 0.);
        // if (electron.Pt()>120 || electron.Pt()<20) cout<< "Electron Reco: pt=" << electron.Pt() << " GeV, sf=" << sf <<endl;

    }
    return effCont;


}

EfficiencyContainer WeightUtils::GetElectronIDEff(TLorentzVector& electron) const
{   
    int bin = _elSF2D_ID->FindBin(electron.Eta(),electron.Pt());
    float sf = _elSF2D_ID->GetBinContent(bin);
    float err = _elSF2D_ID->GetBinError(bin);
    if (_elSF2D_ID->IsBinUnderflow(bin) || _elSF2D_ID->IsBinOverflow(bin)) { sf = 1.0; err = 0.0; }
    EfficiencyContainer effCont(sf, 1., err, 0.);
    // if (electron.Pt()>120 || electron.Pt()<20) cout<< _elSF2D_ID->IsBinUnderflow(bin)<< _elSF2D_ID->IsBinOverflow(bin) <<"Electron Id: pt=" << electron.Pt() << " GeV, sf=" << sf <<endl;
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

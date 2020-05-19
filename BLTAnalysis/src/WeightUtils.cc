#include "BLT/BLTAnalysis/interface/WeightUtils.h"

WeightUtils::WeightUtils(string dataPeriod, string selection)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    std::string fileName;

    rng = new TRandom3();

    const std::string cmssw_base = getenv("CMSSW_BASE");

    // PU weights
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pu_weights/pu_weights_" 
                          + _dataPeriod + "_nom.root";
    std::cout << fileName << std::endl;
    TFile* puFileNom = new TFile(fileName.c_str(), "OPEN");
    _puReweightNom = (TH1D*)puFileNom->Get("mcwei_run000001");
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pu_weights/pu_weights_" 
                          + _dataPeriod + "_up.root";
    std::cout << fileName << std::endl;
    TFile* puFileUp = new TFile(fileName.c_str(), "OPEN");
    _puReweightUp = (TH1D*)puFileUp->Get("mcwei_run000001");
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pu_weights/pu_weights_" 
                          + _dataPeriod + "_down.root";
    std::cout << fileName << std::endl;
    TFile* puFileDown = new TFile(fileName.c_str(), "OPEN");
    _puReweightDown = (TH1D*)puFileDown->Get("mcwei_run000001");

    // muon trigger efficiencies 
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_trigger/EfficienciesAndSF_BCDEF.root";
    TFile* triggerFile = new TFile(fileName.c_str(), "OPEN");

    std::string filePath = "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/";
    _eff_IsoMu24_DATA[0] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin0_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _eff_IsoMu24_DATA[1] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin1_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _eff_IsoMu24_DATA[2] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin2_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());
    _eff_IsoMu24_DATA[3] = (TGraphAsymmErrors*)triggerFile->Get((filePath + "pt_PLOT_abseta_bin3_&_Tight2012_pass_&_tag_IsoMu24_pass_DATA").c_str());

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
    
    // muon loose ID sf (BCDEF)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF_Loose_ID_BCDEF = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Loose_ID_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Loose_ID_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Loose_ID_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Loose_ID_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

    // muon loose ID sf (GH)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_Loose_ID_GH = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesDATA/";
    _muSF_Loose_ID_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_DATA").c_str());
    _muSF_Loose_ID_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_DATA").c_str());

    filePath = "MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/";
    _muSF_Loose_ID_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin0_MC").c_str());
    _muSF_Loose_ID_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin1_MC").c_str());
    _muSF_Loose_ID_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin2_MC").c_str());
    _muSF_Loose_ID_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ID_GH->Get((filePath + "pt_PLOT_abseta_bin3_MC").c_str());

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

    // loose muon ISO sf (BCDEF)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_BCDEF.root";
    TFile* f_muRecoSF_Loose_ISO_BCDEF = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "LooseISO_LooseID_pt_eta/efficienciesDATA/";
    _muSF_Loose_ISO_DATA_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_PF_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_PF_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_PF_pass_DATA").c_str());

    filePath = "LooseISO_LooseID_pt_eta/efficienciesMC/";
    _muSF_Loose_ISO_MC_BCDEF[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
    _muSF_Loose_ISO_MC_BCDEF[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin1_&_PF_pass_MC").c_str());
    _muSF_Loose_ISO_MC_BCDEF[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin2_&_PF_pass_MC").c_str());
    _muSF_Loose_ISO_MC_BCDEF[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_BCDEF->Get((filePath + "pt_PLOT_abseta_bin3_&_PF_pass_MC").c_str());

    // loose muon ISO sf (GH)
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_iso/EfficienciesAndSF_GH.root";
    TFile* f_muRecoSF_Loose_ISO_GH = new TFile(fileName.c_str(), "OPEN"); 

    filePath = "LooseISO_LooseID_pt_eta/efficienciesDATA/";
    _muSF_Loose_ISO_DATA_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());
    _muSF_Loose_ISO_DATA_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_DATA").c_str());

    filePath = "LooseISO_LooseID_pt_eta/efficienciesMC/";
    _muSF_Loose_ISO_MC_GH[0] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
    _muSF_Loose_ISO_MC_GH[1] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
    _muSF_Loose_ISO_MC_GH[2] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
    _muSF_Loose_ISO_MC_GH[3] = (TGraphAsymmErrors*)f_muRecoSF_Loose_ISO_GH->Get((filePath + "pt_PLOT_abseta_bin0_&_PF_pass_MC").c_str());
   
    // hzz muon id efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/hzz_muon_id_sf_" + _dataPeriod + ".root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    _hzz_muIdSF = (TH2F*)f_hzz_muIdSF->Get("FINAL");

    // double electron trigger efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doubleg_trigger/doubleg_leg1_" + _dataPeriod + ".root";
    TFile* f_elTrigSF_leg1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleg_leg1_DATA = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffData2D");
    _eff_doubleg_leg1_MC   = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffMC2D"); 
    _sf_doubleg_leg1 = (TH2F *)f_elTrigSF_leg1->Get("EGamma_SF2D");

    //fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doubleg_trigger/SFs_Leg2_Ele12_HZZSelection_Tag35.root";
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doubleg_trigger/doubleg_leg2_" + _dataPeriod + ".root";
    TFile* f_elTrigSF_leg2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleg_leg2_DATA = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffData2D");
    _eff_doubleg_leg2_MC   = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffMC2D");
    _sf_doubleg_leg2 = (TH2F *)f_elTrigSF_leg2->Get("EGamma_SF2D");

    // double electron 2017 trigger efficiencies from Laurent
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/triggeffcymapsRA5_Run2_ALL.root";
    TFile* f_elTrigSF_Laurent;
    if (_dataPeriod == "2017") {
        _eff_doubleg_leg1_DATA = (TH2F*)f_elTrigSF_Laurent->Get("2017BtoF/ele23|leg24_diele");
        _eff_doubleg_leg1_MC = (TH2F*)f_elTrigSF_Laurent->Get("2017BtoF/ele23|leg24_diele_mc");
        _eff_doubleg_leg2_DATA = (TH2F*)f_elTrigSF_Laurent->Get("2017BtoF/ele12_diele");
        _eff_doubleg_leg2_MC = (TH2F*)f_elTrigSF_Laurent->Get("2017BtoF/ele12_diele_mc");
    }


    // double muon trigger efficiencies
    
    const char* hist_str;
    if (_dataPeriod == "2018") {
        hist_str = "scale_factor";
    }
    else {
        hist_str = "h1";
    }
   
    // leg 1
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu17Leg_Eta0to09_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg1_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_data");
    _eff_doubleMu_leg1_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_mc"); 
    _sf_doubleMu_leg1[0] = (TH1F *)f_DoubleMuTrigSF_leg1_0->Get(hist_str);

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu17Leg_Eta09to12_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg1_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_data");
    _eff_doubleMu_leg1_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_mc"); 
    _sf_doubleMu_leg1[1] = (TH1F *)f_DoubleMuTrigSF_leg1_1->Get(hist_str);
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu17Leg_Eta12to21_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg1_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_data");
    _eff_doubleMu_leg1_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_mc"); 
    _sf_doubleMu_leg1[2] = (TH1F *)f_DoubleMuTrigSF_leg1_2->Get(hist_str);
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu17Leg_Eta21to24_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg1_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_data");
    _eff_doubleMu_leg1_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_mc"); 
    _sf_doubleMu_leg1[3] = (TH1F *)f_DoubleMuTrigSF_leg1_3->Get(hist_str);

    // leg 2 
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu8Leg_Eta0to09_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg2_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_data");
    _eff_doubleMu_leg2_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_mc"); 
    _sf_doubleMu_leg2[0] = (TH1F *)f_DoubleMuTrigSF_leg2_0->Get(hist_str);

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu8Leg_Eta09to12_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg2_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_data");
    _eff_doubleMu_leg2_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_mc"); 
    _sf_doubleMu_leg2[1] = (TH1F *)f_DoubleMuTrigSF_leg2_1->Get(hist_str);
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu8Leg_Eta12to21_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg2_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_data");
    _eff_doubleMu_leg2_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_mc"); 
    _sf_doubleMu_leg2[2] = (TH1F *)f_DoubleMuTrigSF_leg2_2->Get(hist_str);
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/mutrgsf_Mu8Leg_Eta21to24_" + _dataPeriod + ".root";
    TFile* f_DoubleMuTrigSF_leg2_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_data");
    _eff_doubleMu_leg2_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_mc"); 
    _sf_doubleMu_leg2[3] = (TH1F *)f_DoubleMuTrigSF_leg2_3->Get(hist_str);

    // electron reco efficiencies
    //fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_reco_2016" + _dataPeriod + ".root";
    //TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN"); 
    //_eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");
    
    //fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/eleReco_HZZ_Moriond17_SFs.root";
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/eleReco_HZZ_SFs_" + _dataPeriod + ".root";
    TFile* f_eleRecoHZZSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO_2D = (TH2F *)f_eleRecoHZZSF->Get("EGamma_SF2D");

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/eleReco_MVA_SFs_" + _dataPeriod + ".root";
    TFile* f_eleRecoMVASF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_MVA_RECO_2D = (TH2F *)f_eleRecoMVASF->Get("EGamma_SF2D");
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/eleReco_low_MVA_SFs_" + _dataPeriod + ".root";
    TFile* f_eleRecoLowMVASF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_MVA_LOW_RECO_2D = (TH2F *)f_eleRecoLowMVASF->Get("EGamma_SF2D");

    // electron id efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_ID_2016.root";
    TFile* f_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
    _eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
    _eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
    _eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
    _eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");

    // hzz electron id efficiencies
    /*fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/egamma_eff_hzz_ID_2016.root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleSF_ID[0] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_0");
    _hzz_eleSF_ID[1] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_1");
    _hzz_eleSF_ID[2] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_2");
    _hzz_eleSF_ID[3] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_3");
    _hzz_eleSF_ID[4] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_4");
    _hzz_eleSF_ID[5] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_5");
    _hzz_eleSF_ID[6] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_6");
    _hzz_eleSF_ID[7] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_7");
    _hzz_eleSF_ID[8] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_8");
    _hzz_eleSF_ID[9] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_9");
    _hzz_eleSF_ID[10] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_10");
    _hzz_eleSF_ID[11] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_11");
    _hzz_eleSF_ID[12] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_12");*/

    //fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/eleSelectionSF_HZZ_Moriond17.root";
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/eleSelectionSF_HZZ_" + _dataPeriod + ".root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleSF_ID_2D = (TH2F *)f_hzz_eleIdSF->Get("EGamma_SF2D");
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/eleSelectionSF_MVA_" + _dataPeriod + ".root";
    TFile* f_mva_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_MVA_ID_2D = (TH2F *)f_mva_eleIdSF->Get("EGamma_SF2D");

    // photon mva id (90%) efficiencies
    //fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_id/photon_mva_id_2016.root";
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_id/photon_mva_id_" + _dataPeriod + ".root";
    TFile* f_mva_gammaIdSF = new TFile(fileName.c_str(), "OPEN");
    _mva_gammaSF_ID[0] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_0");
    _mva_gammaSF_ID[1] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_1");
    _mva_gammaSF_ID[2] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_2");
    _mva_gammaSF_ID[3] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_3");
    _mva_gammaSF = (TH2F *)f_mva_gammaIdSF->Get("EGamma_SF2D");

    // photon r9 reweighting
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_r9/photon_r9_reweighting_" + _dataPeriod + ".root";
    TFile* f_photon_r9 = new TFile(fileName.c_str(), "OPEN"); 
    _photon_r9_barrel = (TGraph *)f_photon_r9->Get("gR9_EB");
    _photon_r9_endcap = (TGraph *)f_photon_r9->Get("gR9_EE");
}

void WeightUtils::SetDataPeriod(string dataPeriod)
{
    _dataPeriod = dataPeriod;
}

void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
}

std::map<std::string, float> WeightUtils::GetPUWeight(float nPU)
{
    std::map<std::string, float> weights;
    weights["nom"] = _puReweightNom->GetBinContent(_puReweightNom->FindBin(nPU)); 
    weights["up"] = _puReweightUp->GetBinContent(_puReweightUp->FindBin(nPU)); 
    weights["down"] = _puReweightDown->GetBinContent(_puReweightDown->FindBin(nPU)); 
    return weights; 
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
    if (triggerName == "HLT_IsoMu24_v*") {
        if (lepton.Pt() < 500.) {
            effData = _eff_IsoMu24_DATA[etaBin]->Eval(lepton.Pt());
        }
    }
    
    return std::make_pair(effData, effMC);
}

/*std::pair<float,float> WeightUtils::GetDoubleEGTriggerEffWeight(string triggerName, TElectron &electron) const
{
    float effData = 1;
    float effMC = 1;

    if (electron.calibPt < 200.) {
        if (triggerName == "HLT_DoubleEG_leg1") {
            //effData = _eff_doubleg_leg1_DATA->Interpolate(lepton.Eta(), lepton.Pt());
            //effMC   = _eff_doubleg_leg1_MC->Interpolate(lepton.Eta(), lepton.Pt());
            effData = _eff_doubleg_leg1_DATA->GetBinContent(_eff_doubleg_leg1_DATA->FindBin(electron.scEta, electron.calibPt));
            effMC = _eff_doubleg_leg1_MC->GetBinContent(_eff_doubleg_leg1_MC->FindBin(electron.scEta, electron.calibPt));
        }
        else if (triggerName == "HLT_DoubleEG_leg2") {
            //effData = _eff_doubleg_leg2_DATA->Interpolate(lepton.Eta(), lepton.Pt());
            //effMC   = _eff_doubleg_leg2_MC->Interpolate(lepton.Eta(), lepton.Pt());
            effData = _eff_doubleg_leg2_DATA->GetBinContent(_eff_doubleg_leg2_DATA->FindBin(electron.scEta, electron.calibPt));
            effMC = _eff_doubleg_leg2_MC->GetBinContent(_eff_doubleg_leg2_MC->FindBin(electron.scEta, electron.calibPt));
        }
    }

    if (effMC == 0) {
        cout << "zero value for effMC" << endl;
        effMC = 1;
    }

    return std::make_pair(effData, effMC);
}*/

float WeightUtils::GetDoubleEGTriggerEffWeight(string triggerName, TElectron &electron) const
{
    float weight = 1.;
    float tmpElePt = (electron.pt < 200.) ? electron.pt : 199.;

        if (triggerName == "HLT_DoubleEG_leg1") {
            if (_dataPeriod == "2017") {
                float eff_data, eff_mc;
                eff_data = _eff_doubleg_leg1_DATA->GetBinContent(_eff_doubleg_leg1_DATA->FindBin(electron.scEta, tmpElePt));
                eff_mc = _eff_doubleg_leg1_MC->GetBinContent(_eff_doubleg_leg1_MC->FindBin(electron.scEta, tmpElePt));
                weight *= eff_data/eff_mc;
            }
            else {
                weight *= _sf_doubleg_leg1->GetBinContent(_sf_doubleg_leg1->FindBin(electron.scEta, tmpElePt));
            }
        }
        else if (triggerName == "HLT_DoubleEG_leg2") {
            if (_dataPeriod == "2017") {
                float eff_data, eff_mc;
                eff_data = _eff_doubleg_leg2_DATA->GetBinContent(_eff_doubleg_leg2_DATA->FindBin(electron.scEta, tmpElePt));
                eff_mc = _eff_doubleg_leg2_MC->GetBinContent(_eff_doubleg_leg2_MC->FindBin(electron.scEta, tmpElePt));
                weight *= eff_data/eff_mc;
            }
            else {
                weight *= _sf_doubleg_leg2->GetBinContent(_sf_doubleg_leg2->FindBin(electron.scEta, tmpElePt));
            }
        }

        return weight;
}

/*std::pair<float,float> WeightUtils::GetDoubleMuonTriggerEffWeight(string triggerName, TMuon &muon) const
{
    float tmpMuPt = (muon.pt < 200.) ? muon.pt : 199.;
    float effData = 1;
    float effMC = 1;

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.eta) > binningEta[i] && fabs(muon.eta) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    if (triggerName == "HLT_DoubleMuon_leg1") {
        effData = _eff_doubleMu_leg1_DATA[etaBin]->Eval(muon.pt);
        effMC = _eff_doubleMu_leg1_MC[etaBin]->Eval(muon.pt);
        //effData = _eff_doubleMu_leg1_DATA[etaBin]->FindBin(tmpMuPt);
        //effMC = _eff_doubleMu_leg1_MC[etaBin]->FindBin(tmpMuPt);
    }
    else if (triggerName == "HLT_DoubleMuon_leg2") {
        effData = _eff_doubleMu_leg2_DATA[etaBin]->Eval(muon.pt);
        effMC = _eff_doubleMu_leg2_MC[etaBin]->Eval(muon.pt);
        //effData = _eff_doubleMu_leg2_DATA[etaBin]->FindBin(tmpMuPt);
        //effMC = _eff_doubleMu_leg2_MC[etaBin]->FindBin(tmpMuPt);
    }

    if (effMC == 0) {
        cout << "zero value for effMC" << endl;
        effMC = 1;
    }

    return std::make_pair(effData, effMC);
}*/

float WeightUtils::GetDoubleMuonTriggerEffWeight(string triggerName, TMuon &muon) const
{
    float weight = 1.;
    float tmpMuPt = (muon.pt < 200.) ? muon.pt : 199.;

    float binningEta[] = {0., 0.9, 1.2, 2.1, 2.4};
    int etaBin = 0;
    for (int i = 0; i < 4; ++i) {
        if (fabs(muon.eta) > binningEta[i] && fabs(muon.eta) <= binningEta[i+1]) {
            etaBin = i;
            break;
        }
    }

    if (triggerName == "HLT_DoubleMuon_leg1") 
        weight *= _sf_doubleMu_leg1[etaBin]->GetBinContent(_sf_doubleMu_leg1[etaBin]->FindBin(tmpMuPt));
    else if (triggerName == "HLT_DoubleMuon_leg2") 
        weight *= _sf_doubleMu_leg2[etaBin]->GetBinContent(_sf_doubleMu_leg2[etaBin]->FindBin(tmpMuPt));

    return weight;
}

float WeightUtils::GetMuonIDEff(TLorentzVector& muon) const
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
    float random = rng->Rndm();
    if (muon.Pt() < 200.) {
        if (random > 0.468) {
            weight   *= _muSF_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
        } else {
            weight   *= _muSF_ID_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_ID_MC_GH[etaBin]->Eval(muon.Pt());
        }
    }
    
    return weight;
}

float WeightUtils::GetLooseMuonIDEff(TLorentzVector& muon) const
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
    float random = rng->Rndm();
    if (muon.Pt() < 200.) {
        if (random > 0.468) {
            weight   *= _muSF_Loose_ID_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_Loose_ID_MC_BCDEF[etaBin]->Eval(muon.Pt());
        } else {
            weight   *= _muSF_Loose_ID_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_Loose_ID_MC_GH[etaBin]->Eval(muon.Pt());
        }
    }
    
    return weight;
}

float WeightUtils::GetMuonISOEff(TLorentzVector& muon) const
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
    float random = rng->Rndm();
    if (muon.Pt() < 200.) {
        if (random > 0.468) {
            weight   *= _muSF_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
        } else {
            weight   *= _muSF_ISO_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_ISO_MC_GH[etaBin]->Eval(muon.Pt());
        }
    }
    
    return weight;
}

float WeightUtils::GetLooseMuonISOEff(TLorentzVector& muon) const
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
    float random = rng->Rndm();
    if (muon.Pt() < 200.) {
        if (random > 0.468) {
            weight   *= _muSF_Loose_ISO_DATA_BCDEF[etaBin]->Eval(muon.Pt())/_muSF_Loose_ISO_MC_BCDEF[etaBin]->Eval(muon.Pt());
        } else {
            weight   *= _muSF_Loose_ISO_DATA_GH[etaBin]->Eval(muon.Pt())/_muSF_Loose_ISO_MC_GH[etaBin]->Eval(muon.Pt());
        }
    }
    
    return weight;
}

float WeightUtils::GetHZZMuonIDEff(TMuon& muon) const
{
    float weight = 1;
    float tmpMuPt = (muon.pt < 200.) ? muon.pt : 199.;
    
    //weight *= _hzz_muIdSF->Interpolate(muon.eta, muon.pt);
    weight *= _hzz_muIdSF->GetBinContent(_hzz_muIdSF->FindBin(muon.eta, tmpMuPt));
    
    return weight;
}



float WeightUtils::GetElectronRecoIdEff(TLorentzVector& electron) const
{
    float binningPt[] = {10., 20., 35., 50., 90., 500.};
    int ptBin = 0;
    for (int i = 0; i < 5; ++i) {
        if (fabs(electron.Pt()) > binningPt[i] && fabs(electron.Pt()) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }

    float weight = 1;
    if (electron.Pt() < 500.) {
        weight *= _eleSF_RECO->Eval(electron.Eta());
        weight *= _eleSF_ID[ptBin]->Eval(electron.Eta());
    }
    
    return weight;
}

float WeightUtils::GetHZZElectronRecoIdEff(TElectron& electron) const 
{
    float tmpElePt = (electron.calibPt < 200.) ? electron.calibPt : 199.;
    float weight = 1;

    weight *= _eleSF_RECO_2D->GetBinContent(_eleSF_RECO_2D->FindBin(electron.scEta, 50.));
    weight *= _hzz_eleSF_ID_2D->GetBinContent(_hzz_eleSF_ID_2D->FindBin(fabs(electron.scEta), tmpElePt));
    
    return weight;
}

float WeightUtils::GetElectronMVARecoIdEff(TElectron& electron) const
{
    float tmpElePt = (electron.calibPt < 200.) ? electron.calibPt : 199.;
    float weight = 1;

    if (tmpElePt < 20.) {
        weight *= _eleSF_MVA_LOW_RECO_2D->GetBinContent(_eleSF_MVA_LOW_RECO_2D->FindBin(electron.scEta, 15.));
    }
    else {
        weight *= _eleSF_MVA_RECO_2D->GetBinContent(_eleSF_MVA_RECO_2D->FindBin(electron.scEta, tmpElePt));
    }

    weight *= _eleSF_MVA_ID_2D->GetBinContent(_eleSF_MVA_ID_2D->FindBin(electron.scEta, tmpElePt));
    
    return weight;
}
    

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
    if (fabs(photon.scEta) <= 1.49) {
        if      (_dataPeriod == "2016") weight *= 0.9938;
        else if (_dataPeriod == "2017") weight *= 0.967;
        else if (_dataPeriod == "2018") weight *= 0.967;
    }
    else if (fabs(photon.scEta) > 1.49) {
        if      (_dataPeriod == "2016") weight *= 0.9875;
        else if (_dataPeriod == "2017") weight *= 0.915;
        else if (_dataPeriod == "2018") weight *= 0.915;
    }

    if (weight == 0)
        weight = 1.;
    
    return weight;
}

float WeightUtils::GetCorrectedPhotonR9(TPhoton& photon) const 
{
    float r9 = photon.r9_full5x5;
    if (fabs(photon.scEta) < 1.444)
        r9 = _photon_r9_barrel->Eval(photon.r9_full5x5);
    else if (fabs(photon.scEta) > 1.566)
        r9 = _photon_r9_endcap->Eval(photon.r9_full5x5);
    else
        std::cout << "bad value of photon scEta: returning original r9" << std::endl;
    
    return r9;
}


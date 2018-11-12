#include "BLT/BLTAnalysis/interface/WeightUtils.h"

WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;
    std::string fileName;

    rng = new TRandom3();


    const std::string cmssw_base = getenv("CMSSW_BASE");
    // PU weights
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/pileup_sf_2016_full.root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TGraph*)puFile->Get("pileup_sf");

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
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/muon_id/hzz_muon_id_sf.root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    _hzz_muIdSF = (TH2F*)f_hzz_muIdSF->Get("FINAL");

    // double electron trigger efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doubleg_trigger/SFs_Leg1_Ele23_HZZSelection_Tag35.root";
    TFile* f_elTrigSF_leg1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleg_leg1_DATA = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffData2D");
    _eff_doubleg_leg1_MC   = (TH2F*)f_elTrigSF_leg1->Get("EGamma_EffMC2D"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doubleg_trigger/SFs_Leg2_Ele12_HZZSelection_Tag35.root";
    TFile* f_elTrigSF_leg2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleg_leg2_DATA = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffData2D");
    _eff_doubleg_leg2_MC   = (TH2F*)f_elTrigSF_leg2->Get("EGamma_EffMC2D");

    // double muon trigger efficiencies
   
    // leg 1
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu17Leg_Eta0to09.root";
    TFile* f_DoubleMuTrigSF_leg1_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_data");
    _eff_doubleMu_leg1_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu17Leg_Eta09to12.root";
    TFile* f_DoubleMuTrigSF_leg1_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_data");
    _eff_doubleMu_leg1_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_mc"); 
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu17Leg_Eta12to21.root";
    TFile* f_DoubleMuTrigSF_leg1_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_data");
    _eff_doubleMu_leg1_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_mc"); 
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu17Leg_Eta21to24.root";
    TFile* f_DoubleMuTrigSF_leg1_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_data");
    _eff_doubleMu_leg1_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_mc"); 

    // leg 2 
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu8Leg_Eta0to09.root";
    TFile* f_DoubleMuTrigSF_leg2_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_data");
    _eff_doubleMu_leg2_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_mc"); 

    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu8Leg_Eta09to12.root";
    TFile* f_DoubleMuTrigSF_leg2_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_data");
    _eff_doubleMu_leg2_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_mc"); 
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu8Leg_Eta12to21.root";
    TFile* f_DoubleMuTrigSF_leg2_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_data");
    _eff_doubleMu_leg2_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_mc"); 
    
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/doublemuon_trigger/sf_Mu8Leg_Eta21to24.root";
    TFile* f_DoubleMuTrigSF_leg2_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_data");
    _eff_doubleMu_leg2_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_mc"); 

    // electron reco efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_reco_2016.root";
    TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO = (TGraphErrors*)f_eleRecoSF->Get("grSF1D_0");

    // electron id efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/egamma_eff_ID_2016.root";
    TFile* f_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_ID[0] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_0");
    _eleSF_ID[1] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_1");
    _eleSF_ID[2] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_2");
    _eleSF_ID[3] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_3");
    _eleSF_ID[4] = (TGraphErrors*)f_eleIdSF->Get("grSF1D_4");

    // hzz electron id efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/electron_id/egamma_eff_hzz_ID_2016.root";
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
    _hzz_eleSF_ID[12] = (TGraphErrors*)f_hzz_eleIdSF->Get("grSF1D_12");

    // photon mva id (90%) efficiencies
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_id/photon_mva_id_2016.root";
    TFile* f_mva_gammaIdSF = new TFile(fileName.c_str(), "OPEN");
    _mva_gammaSF_ID[0] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_0");
    _mva_gammaSF_ID[1] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_1");
    _mva_gammaSF_ID[2] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_2");
    _mva_gammaSF_ID[3] = (TGraphErrors*)f_mva_gammaIdSF->Get("grSF1D_3");
    _mva_gammaSF = (TH2F *)f_mva_gammaIdSF->Get("EGamma_SF2D");

    // photon r9 reweighting
    fileName = cmssw_base + "/src/BLT/BLTAnalysis/data/photon_r9_reweighting_2016.root";
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
    if (triggerName == "HLT_IsoMu24_v*") {
        if (lepton.Pt() < 500.) {
            effData = _eff_IsoMu24_DATA[etaBin]->Eval(lepton.Pt());
        }
    }
    
    return std::make_pair(effData, effMC);
}

std::pair<float,float> WeightUtils::GetDoubleEGTriggerEffWeight(string triggerName, TElectron &electron) const
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
}

std::pair<float,float> WeightUtils::GetDoubleMuonTriggerEffWeight(string triggerName, TMuon &muon) const
{
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

    if (muon.pt < 200.) {
        if (triggerName == "HLT_DoubleMuon_leg1") {
            effData = _eff_doubleMu_leg1_DATA[etaBin]->Eval(muon.pt);
            effMC = _eff_doubleMu_leg1_MC[etaBin]->Eval(muon.pt);
        }
        else if (triggerName == "HLT_DoubleMuon_leg2") {
            effData = _eff_doubleMu_leg2_DATA[etaBin]->Eval(muon.pt);
            effMC = _eff_doubleMu_leg2_MC[etaBin]->Eval(muon.pt);
        }
    }

    if (effMC == 0) {
        cout << "zero value for effMC" << endl;
        effMC = 1;
    }

    return std::make_pair(effData, effMC);
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
    if (muon.pt < 200.) 
        weight *= _hzz_muIdSF->Interpolate(muon.eta, muon.pt);
    
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
    float binningPt[] = {7., 15., 20., 30., 40., 50., 60., 70., 80., 100., 120., 140., 160., 200.}; 
    int ptBin = 0;
    for (int i = 0; i < 13; ++i) {
        if (fabs(electron.calibPt) > binningPt[i] && fabs(electron.calibPt) <= binningPt[i+1]) {
            ptBin = i;
            break;
        }
    }

    float weight = 1;
    if (electron.calibPt < 200.) {
        weight *= _eleSF_RECO->Eval(electron.scEta);
        weight *= _hzz_eleSF_ID[ptBin]->Eval(electron.scEta);
    }
    
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

    float weight = 1;
    //if (photon.calibPt < 150.) {
    //    weight *= _mva_gammaSF_ID[ptBin]->Eval(photon.scEta);
    //}
    weight *= _mva_gammaSF->GetBinContent(_mva_gammaSF->FindBin(photon.scEta, photon.pt));

    // electron veto scale factor
    if (fabs(photon.scEta) <= 1.49)
        weight *= 0.9938;
    else if (fabs(photon.scEta) > 1.49)
        weight *= 0.9875;
    
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


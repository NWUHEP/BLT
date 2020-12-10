#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys

''' Specify parameters '''
cfg        = bm.JobConfig
selection = 'mmg'
periods = [2016]
executable = 'execBatch.sh'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

muon_data_16 = []
muon_data_16.extend([

    # Double muon data
    cfg(data_name = 'muon_2016B_part0',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016B-17Jul2018_ver2-v1/200923_195738/0000/',
        nJobs    = 25,
        suffix   = 'muon_2016B'
       ),
    cfg(data_name = 'muon_2016B_part1',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016B-17Jul2018_ver2-v1/200923_195738/0001/',
        nJobs    = 25,
        suffix   = 'muon_2016B'
       ),

    cfg(data_name = 'muon_2016C',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016C-17Jul2018-v1/200923_195917/0000/',
        nJobs    = 50,
        suffix   = 'muon_2016C'
       ),

    cfg(data_name = 'muon_2016D_part0',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016D-17Jul2018-v1/200923_210305/0000/',
        nJobs    = 25,
        suffix   = 'muon_2016D'
       ),
    cfg(data_name = 'muon_2016D_part1',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016D-17Jul2018-v1/200923_210305/0001/',
        nJobs    = 25,
        suffix   = 'muon_2016D'
       ),
       
    cfg(data_name = 'muon_2016E',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016E-17Jul2018-v1/200923_200231/0000/',
        nJobs    = 50,
        suffix   = 'muon_2016E'
       ),
    cfg(data_name = 'muon_2016F',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016F-17Jul2018-v1/200923_200407/0000/',
        nJobs    = 50,
        suffix   = 'muon_2016F'
       ),
    cfg(data_name = 'muon_2016G',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016G-17Jul2018-v1/200923_200545/0000/',
        nJobs    = 50,
        suffix   = 'muon_2016G'
       ),
    cfg(data_name = 'muon_2016H_part0',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016H-17Jul2018-v1/200923_200722/0000/',
        nJobs    = 25,
        suffix   = 'muon_2016H'
       ),
    cfg(data_name = 'muon_2016H_part1',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleMuon/2016_data_prod_TEST_V2_DoubleMuon_Run2016H-17Jul2018-v1/200923_200722/0001/',
        nJobs    = 25,
        suffix   = 'muon_2016H'
       ),
    ])

electron_data_16 = []
electron_data_16.extend([

    # Double electron data
    cfg(data_name = 'electron_2016B_part0',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleEG/2016_data_prod_DoubleEG_B_retry_DoubleEG_Run2016B-17Jul2018_ver2-v1/201001_013436/0000/',
        nJobs    = 25,
        suffix   = 'electron_2016B'
       ),
    cfg(data_name = 'electron_2016B_part1',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleEG/2016_data_prod_DoubleEG_B_retry_DoubleEG_Run2016B-17Jul2018_ver2-v1/201001_013436/0001/',
        nJobs    = 25,
        suffix   = 'electron_2016B'
       ),
    cfg(data_name = 'electron_2016B_part2',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleEG/2016_data_prod_DoubleEG_B_retry_DoubleEG_Run2016B-17Jul2018_ver2-v1/201001_013436/0002/',
        nJobs    = 25,
        suffix   = 'electron_2016B'
       ),
    cfg(data_name = 'electron_2016B_part3',
        path     = '/eos/uscms/store/group/lpcbacon/jbueghly/DoubleEG/2016_data_prod_DoubleEG_B_retry_DoubleEG_Run2016B-17Jul2018_ver2-v1/201001_013436/0003/',
        nJobs    = 25,
        suffix   = 'electron_2016B'
       ),

    cfg(data_name = 'electron_2016C',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016C-17Jul2018-v1/200814_175814/0000/',
        nJobs    = 50,
        suffix   = 'electron_2016C'
       ),

    cfg(data_name = 'electron_2016D',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016D-17Jul2018-v1/200814_180009/0000/',
        nJobs    = 50,
        suffix   = 'electron_2016D'
       ),

    cfg(data_name = 'electron_2016E_part0',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016E-17Jul2018-v1/200814_180154/0000/',
        nJobs    = 25,
        suffix   = 'electron_2016E'
       ),
    cfg(data_name = 'electron_2016E_part1',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016E-17Jul2018-v1/200814_180154/0001/',
        nJobs    = 25,
        suffix   = 'electron_2016E'
       ),
    cfg(data_name = 'electron_2016E_recovery',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_recovery_DoubleEG_Run2016E-17Jul2018-v1/200901_230507/0000/',
        nJobs    = 21,
        suffix   = 'electron_2016E'
       ),

    cfg(data_name = 'electron_2016F',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016F-17Jul2018-v1/200814_180416/0000/',
        nJobs    = 50,
        suffix   = 'electron_2016F'
       ),
    cfg(data_name = 'electron_2016G_part0',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016G-17Jul2018-v1/200814_180731/0000/',
        nJobs    = 25,
        suffix   = 'electron_2016G'
       ),
    cfg(data_name = 'electron_2016G_part1',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016G-17Jul2018-v1/200814_180731/0001/',
        nJobs    = 25,
        suffix   = 'electron_2016G'
       ),
    cfg(data_name = 'electron_2016G_part2',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016G-17Jul2018-v1/200814_180731/0002/',
        nJobs    = 25,
        suffix   = 'electron_2016G'
       ),
    cfg(data_name = 'electron_2016G_recovery',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_recovery_DoubleEG_Run2016G-17Jul2018-v1/200902_144203/0000/',
        nJobs    = 26,
        suffix   = 'electron_2016G'
       ),

    cfg(data_name = 'electron_2016H_part0',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016H-17Jul2018-v1/200814_180922/0000/',
        nJobs    = 25,
        suffix   = 'electron_2016H'
       ),
    cfg(data_name = 'electron_2016H_part1',
        path     = '/eos/uscms/store/group/lpcbacon/mmackenz/DoubleEG/2016_data_prod_TEST_V2_DoubleEG_Run2016H-17Jul2018-v1/200814_180922/0001/',
        nJobs    = 25,
        suffix   = 'electron_2016H'
       ),
    ])
 
mc_16 = []
mc_16.extend([

    cfg(data_name = 'DYJetsToLL_M50_2016', 
        path      = '/eos/uscms/store/group/lpcbacon/corderom/james/mc_2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_prod_DYJets_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRu/200922_144944/0000/',
        nJobs     = 50,
        suffix    = 'zjets_M50_2016'
        ),

    cfg(data_name = 'ZGToLLG_2016', 
        path      = '/eos/uscms/store/group/lpcbacon/zichen/ZGToLLG_01J_LoosePtlPtg_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_prod_TEST_V2_ZGToLLG_01J_LoosePtlPtg_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond/200812_184026/0000/',
        nJobs     = 50,
        suffix    = 'zg_llg_2016'
        ),
    cfg(data_name = 'LLAJJ_EWK_M50_2016', 
        path      = '/eos/uscms/store/group/lpcbacon/zichen/LLAJJ_EWK_MLL-50_MJJ-120_13TeV-madgraph-pythia8/2016_mc_prod_TEST_V2_LLAJJ_EWK_MLL-50_MJJ-120_13TeV-madgraph-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymp/200812_184341/0000/',
        nJobs     = 50,
        suffix    = 'zg_ewk_2016'
        ),

    cfg(data_name = 'TTJets_2016_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/zichen/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/2016_mc_prod_TEST_V2_TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asym/200812_184203/0000/',
        nJobs     = 25,
        suffix    = 'ttjets_2016'
        ),
    cfg(data_name = 'TTJets_2016_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/zichen/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/2016_mc_prod_TEST_V2_TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asym/200812_184203/0001/',
        nJobs     = 25,
        suffix    = 'ttjets_2016'
        ),
    cfg(data_name = 'TTJets_2016_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/zichen/TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8/2016_mc_prod_TEST_V2_TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asym/200812_184203/0002/',
        nJobs     = 25,
        suffix    = 'ttjets_2016'
        ),

    cfg(data_name = 'hzg_gluglu_M125_2016',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/GluGluHToZG_M-125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_GluGluHToZG_M-125_13TeV_powheg_pythia8_cmkuo-ggF_HZg_M125_Summer16_minAODv3_v1-bd3e7bcff6c9bcad356e/200914_142146/0000/',
        nJobs     = 25,
        suffix    = 'hzg_gluglu_M125_2016'
        ),
    cfg(data_name = 'hzg_vbf_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/VBFHToZG_M-125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_VBFHToZG_M-125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext/200811_015754/0000/',
        nJobs     = 25,
        suffix    = 'hzg_vbf_M125_2016'
        ),
    cfg(data_name = 'hzg_wplush_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_a/200811_020242/0000/',
        nJobs     = 3,
        suffix    = 'hzg_wplush_M125_2016'
        ),
    cfg(data_name = 'hzg_wminush_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_/200811_020423/0000/',
        nJobs     = 9,
        suffix    = 'hzg_wminush_M125_2016'
        ),
    cfg(data_name = 'hzg_zh_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymp/200811_020105/0000/',
        nJobs     = 7,
        suffix    = 'hzg_zh_M125_2016'
        ),
    cfg(data_name = 'hzg_tth_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ttHToZG_ZToLL_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ttHToZG_ZToLL_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3/200811_015929/0000/',
        nJobs     = 5, 
        suffix    = 'hzg_tth_M125_2016'
        ),
    
    cfg(data_name = 'hzg_gluglu_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/GluGluHToZG_M-120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_GluGluHToZG_M-120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_/200811_020604/0000/',
        nJobs     = 14,
        suffix    = 'hzg_gluglu_M120_2016'
        ),
    cfg(data_name = 'hzg_vbf_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/VBFHToZG_M-120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_VBFHToZG_M-120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/200811_020740/0000/',
        nJobs     = 16,
        suffix    = 'hzg_vbf_M120_2016'
        ),
    cfg(data_name = 'hzg_wplush_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_a/200811_021225/0000/',
        nJobs     = 9,
        suffix    = 'hzg_wplush_M120_2016'
        ),
    cfg(data_name = 'hzg_wminush_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_/200811_021400/0000/',
        nJobs     = 9,
        suffix    = 'hzg_wminush_M120_2016'
        ),
    cfg(data_name = 'hzg_zh_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymp/200811_021050/0000/',
        nJobs     = 8,
        suffix    = 'hzg_zh_M120_2016'
        ),
    cfg(data_name = 'hzg_tth_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ttHToZG_ZToLL_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ttHToZG_ZToLL_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3/200811_020914/0000/',
        nJobs     = 6, 
        suffix    = 'hzg_tth_M120_2016'
        ),
    
    cfg(data_name = 'hzg_gluglu_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/GluGluHToZG_M-130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_GluGluHToZG_M-130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_/200811_021546/0000/',
        nJobs     = 15,
        suffix    = 'hzg_gluglu_M130_2016'
        ),
    cfg(data_name = 'hzg_vbf_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/VBFHToZG_M-130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_VBFHToZG_M-130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/200811_021720/0000/',
        nJobs     = 15,
        suffix    = 'hzg_vbf_M130_2016'
        ),
    cfg(data_name = 'hzg_wplush_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_a/200811_022154/0000/',
        nJobs     = 5,
        suffix    = 'hzg_wplush_M130_2016'
        ),
    cfg(data_name = 'hzg_wminush_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_/200811_022326/0000/',
        nJobs     = 5,
        suffix    = 'hzg_wminush_M130_2016'
        ),
    cfg(data_name = 'hzg_zh_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymp/200811_022022/0000/',
        nJobs     = 4,
        suffix    = 'hzg_zh_M130_2016'
        ),
    cfg(data_name = 'hzg_tth_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ttHToZG_ZToLL_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ttHToZG_ZToLL_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3/200811_021851/0000/',
        nJobs     = 12, 
        suffix    = 'hzg_tth_M130_2016'
        ),
    
    cfg(data_name = 'hmumu_gluglu_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/2016_mc_prod_TEST_V2_GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIISummer16MiniAODv3-PUMoriond17_/200812_152118/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_gluglu_M125_2016'
        ),
    cfg(data_name = 'hmumu_vbf_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/2016_mc_prod_TEST_V2_VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcR/200812_152247/0000/',
        nJobs     = 10,
        suffix    = 'hmumu_vbf_M125_2016'
        ),
    cfg(data_name = 'hmumu_wplush_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WPlusH_HToMuMu_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WPlusH_HToMuMu_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v/200812_152715/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_wplush_M125_2016'
        ),
    cfg(data_name = 'hmumu_wminush_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WMinusH_HToMuMu_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WMinusH_HToMuMu_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_/200812_152847/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_wminush_M125_2016'
        ),
    cfg(data_name = 'hmumu_zh_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ZH_HToMuMu_M125_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ZH_HToMuMu_M125_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/200812_152547/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_zh_M125_2016'
        ),
    cfg(data_name = 'hmumu_tth_M125_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/2016_mc_prod_TEST_V2_ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2/200812_152417/0000/',
        nJobs     = 25, 
        suffix    = 'hmumu_tth_M125_2016'
        ),
    
    cfg(data_name = 'hmumu_gluglu_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/GluGluHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/2016_mc_prod_TEST_V2_GluGluHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIISummer16MiniAODv3-PUMoriond17_/200812_153016/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_gluglu_M120_2016'
        ),
    cfg(data_name = 'hmumu_vbf_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/VBFHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/2016_mc_prod_TEST_V2_VBFHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnlo_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcR/200812_153158/0000/',
        nJobs     = 9,
        suffix    = 'hmumu_vbf_M120_2016'
        ),
    cfg(data_name = 'hmumu_wplush_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WPlusH_HToMuMu_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WPlusH_HToMuMu_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v/200812_153633/0000/',
        nJobs     = 3,
        suffix    = 'hmumu_wplush_M120_2016'
        ),
    cfg(data_name = 'hmumu_wminush_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WMinusH_HToMuMu_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WMinusH_HToMuMu_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_/200812_153802/0000/',
        nJobs     = 4,
        suffix    = 'hmumu_wminush_M120_2016'
        ),
    cfg(data_name = 'hmumu_zh_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ZH_HToMuMu_M120_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ZH_HToMuMu_M120_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/200812_153502/0000/',
        nJobs     = 12,
        suffix    = 'hmumu_zh_M120_2016'
        ),
    cfg(data_name = 'hmumu_tth_M120_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ttHToMuMu_M120_TuneCP5_PSweights_13TeV-powheg-pythia8/2016_mc_prod_TEST_V2_ttHToMuMu_M120_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2/200812_153329/0000/',
        nJobs     = 25, 
        suffix    = 'hmumu_tth_M120_2016'
        ),
    
    cfg(data_name = 'hmumu_gluglu_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/GluGluHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/2016_mc_prod_TEST_V2_GluGluHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIISummer16MiniAODv3-PUMoriond17_/200812_153956/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_gluglu_M130_2016'
        ),
    cfg(data_name = 'hmumu_vbf_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/VBFHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/2016_mc_prod_TEST_V2_VBFHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnlo_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcR/200812_154136/0000/',
        nJobs     = 5,
        suffix    = 'hmumu_vbf_M130_2016'
        ),
    cfg(data_name = 'hmumu_wplush_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WPlusH_HToMuMu_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WPlusH_HToMuMu_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v/200812_154620/0000/',
        nJobs     = 5,
        suffix    = 'hmumu_wplush_M130_2016'
        ),
    cfg(data_name = 'hmumu_wminush_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/WMinusH_HToMuMu_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_WMinusH_HToMuMu_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_/200812_154808/0000/',
        nJobs     = 7,
        suffix    = 'hmumu_wminush_M130_2016'
        ),
    cfg(data_name = 'hmumu_zh_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ZH_HToMuMu_M130_13TeV_powheg_pythia8/2016_mc_prod_TEST_V2_ZH_HToMuMu_M130_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/200812_154444/0000/',
        nJobs     = 6,
        suffix    = 'hmumu_zh_M130_2016'
        ),
    cfg(data_name = 'hmumu_tth_M130_2016',
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2016/ttHToMuMu_M130_TuneCP5_PSweights_13TeV-powheg-pythia8/2016_mc_prod_TEST_V2_ttHToMuMu_M130_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2/200812_154311/0000/',
        nJobs     = 25, 
        suffix    = 'hmumu_tth_M130_2016'
        ),

    ])

batch_16_samples = []
if selection == 'mmg':
    batch_16_samples += muon_data_16
elif selection == 'eeg':
    batch_16_samples += electron_data_16
batch_16_samples += mc_16

batch_16 = bm.BatchMaster(config_list = batch_16_samples, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = 2016,
                       executable  = executable,
                       location    = 'lpc'
                     )

muon_data_17 = []
muon_data_17.extend([

    # double muon 2017
    cfg(data_name = 'muon_2017B',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017B-31Mar2018-v1/200814_225141/0000/',
        nJobs    = 50,
        suffix   = 'muon_2017B'
       ),

    cfg(data_name = 'muon_2017C_part0',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017C-31Mar2018-v1/200814_225327/0000/',
        nJobs    = 25,
        suffix   = 'muon_2017C'
       ),
    cfg(data_name = 'muon_2017C_part1',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017C-31Mar2018-v1/200814_225327/0001/',
        nJobs    = 25,
        suffix   = 'muon_2017C'
       ),

    cfg(data_name = 'muon_2017D',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017D-31Mar2018-v1/200814_225525/0000/',
        nJobs    = 50,
        suffix   = 'muon_2017D'
       ),

    cfg(data_name = 'muon_2017E_part0',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017E-31Mar2018-v1/200814_225948/0000/',
        nJobs    = 25,
        suffix   = 'muon_2017E'
       ),
    cfg(data_name = 'muon_2017E_part1',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017E-31Mar2018-v1/200814_225948/0001/',
        nJobs    = 25,
        suffix   = 'muon_2017E'
       ),
    cfg(data_name = 'muon_2017E_recovery',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_recovery_DoubleMuon_Run2017E-31Mar2018-v1_recovery/200917_214422/0000/',
        nJobs    = 15,
        suffix   = 'muon_2017E'
       ),

    cfg(data_name = 'muon_2017F_part0',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017F-31Mar2018-v1/200814_230348/0000/',
        nJobs    = 25,
        suffix   = 'muon_2017F'
       ),
    cfg(data_name = 'muon_2017F_part1',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017F-31Mar2018-v1/200814_230348/0001/',
        nJobs    = 25,
        suffix   = 'muon_2017F'
       ),
    cfg(data_name = 'muon_2017F_part2',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017F-31Mar2018-v1/200814_230348/0002/',
        nJobs    = 25,
        suffix   = 'muon_2017F'
       ),
    cfg(data_name = 'muon_2017F_part3',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_DoubleMuon_Run2017F-31Mar2018-v1/200814_230348/0003/',
        nJobs    = 25,
        suffix   = 'muon_2017F'
       ),
    cfg(data_name = 'muon_2017F_recovery',
        path     = '/eos/uscms/store/group/lpcbacon/naodell/DoubleMuon/2017_data_prod_recovery_DoubleMuon_Run2017F-31Mar2018-v1_recovery/200917_214601/0000/',
        nJobs    = 4,
        suffix   = 'muon_2017F'
       ),
    ])

electron_data_17 = []
electron_data_17.extend([
    
    # DoubleEG 2017
    cfg(data_name = 'electron_2017B_part0',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017B-31Mar2018-v1/200828_210815/0000/',
        nJobs    = 25,
        suffix   = 'electron_2017B'
       ),
    cfg(data_name = 'electron_2017B_part1',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017B-31Mar2018-v1/200828_210815/0001/',
        nJobs    = 25,
        suffix   = 'electron_2017B'
       ),
    cfg(data_name = 'electron_2017B_part2',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017B-31Mar2018-v1/200828_210815/0002/',
        nJobs    = 25,
        suffix   = 'electron_2017B'
       ),

    cfg(data_name = 'electron_2017C_part0',
        path     = '/eos/uscms/store/group/lpcbacon/corderom/james/data_2017/DoubleEG/2017_data_prod_DoubleEG_Run2017C-31Mar2018-v1/200828_211011/0000/',
        nJobs    = 25,
        suffix   = 'electron_2017C'
       ),
    cfg(data_name = 'electron_2017C_part1',
        path     = '/eos/uscms/store/group/lpcbacon/corderom/james/data_2017/DoubleEG/2017_data_prod_DoubleEG_Run2017C-31Mar2018-v1/200828_211011/0001/',
        nJobs    = 25,
        suffix   = 'electron_2017C'
       ),
    cfg(data_name = 'electron_2017C_part2',
        path     = '/eos/uscms/store/group/lpcbacon/corderom/james/data_2017/DoubleEG/2017_data_prod_DoubleEG_Run2017C-31Mar2018-v1/200828_211011/0002/',
        nJobs    = 25,
        suffix   = 'electron_2017C'
       ),

    cfg(data_name = 'electron_2017D',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017D-31Mar2018-v1/200828_211201/0000/',
        nJobs    = 50,
        suffix   = 'electron_2017D'
       ),

    cfg(data_name = 'electron_2017E_part0',
        path     = '/eos/uscms/store/group/lpcbacon/corderom/james/data_2017/DoubleEG/2017_data_prod_DoubleEG_Run2017E-31Mar2018-v1/200828_211358/0000/',
        nJobs    = 25,
        suffix   = 'electron_2017E'
       ),
    cfg(data_name = 'electron_2017E_part1',
        path     = '/eos/uscms/store/group/lpcbacon/corderom/james/data_2017/DoubleEG/2017_data_prod_DoubleEG_Run2017E-31Mar2018-v1/200828_211358/0001/',
        nJobs    = 25,
        suffix   = 'electron_2017E'
       ),
    cfg(data_name = 'electron_2017E_part2',
        path     = '/eos/uscms/store/group/lpcbacon/corderom/james/data_2017/DoubleEG/2017_data_prod_DoubleEG_Run2017E-31Mar2018-v1/200828_211358/0002/',
        nJobs    = 25,
        suffix   = 'electron_2017E'
       ),

    cfg(data_name = 'electron_2017F_part0',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017F-31Mar2018-v1/200828_211551/0000/',
        nJobs    = 25,
        suffix   = 'electron_2017F'
       ),
    cfg(data_name = 'electron_2017F_part1',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017F-31Mar2018-v1/200828_211551/0001/',
        nJobs    = 25,
        suffix   = 'electron_2017F'
       ),
    cfg(data_name = 'electron_2017F_part2',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017F-31Mar2018-v1/200828_211551/0002/',
        nJobs    = 25,
        suffix   = 'electron_2017F'
       ),
    cfg(data_name = 'electron_2017F_part3',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017F-31Mar2018-v1/200828_211551/0003/',
        nJobs    = 25,
        suffix   = 'electron_2017F'
       ),
    cfg(data_name = 'electron_2017F_part4',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017F-31Mar2018-v1/200828_211551/0004/',
        nJobs    = 25,
        suffix   = 'electron_2017F'
       ),
    cfg(data_name = 'electron_2017F_part5',
        path     = '/eos/uscms/store/group/lpchzg/james/data_2017/2017_data_prod_DoubleEG_Run2017F-31Mar2018-v1/200828_211551/0005/',
        nJobs    = 25,
        suffix   = 'electron_2017F'
       ),
   ])

mc_17 = []
mc_17.extend([
    #cfg(data_name = 'DYJetsToLL_M50_2017_part0',
    #    path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_retry_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201021_224308/0000/',
    #    nJobs     = 25,
    #    suffix    = 'zjets_M50_2017'
    #    ),
    #cfg(data_name = 'DYJetsToLL_M50_2017_part1',
    #    path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_retry_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201021_224308/0001/',
    #    nJobs     = 25,
    #    suffix    = 'zjets_M50_2017'
    #    ),
    
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part0',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0000/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part1',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0001/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part2',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0002/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part3',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0003/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part4',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0004/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part5',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0005/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part6',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0006/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part7',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0007/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part8',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0008/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2017_ext_part9',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_zjets_extended_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/201024_152614/0009/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2017'
        ),

    cfg(data_name = 'ZGToLLG_2017_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_zg_llg_retry_ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018/201021_220028/0000/',
        nJobs     = 25,
        suffix    = 'zg_llg_2017'
        ),
    cfg(data_name = 'ZGToLLG_2017_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_zg_llg_retry_ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018/201021_220028/0001/',
        nJobs     = 25,
        suffix    = 'zg_llg_2017'
        ),

    cfg(data_name = 'LLAJJ_EWK_M50_2017', 
        path      = '/eos/uscms/store/user/corderom/JAMES/mc_2017/LLAJJ_EWK_MLL-50_MJJ-120_TuneCP5_13TeV-madgraph-pythia8/2017_mc_prod_LLAJJ_EWK_MLL-50_MJJ-120_TuneCP5_13TeV-madgraph-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_m/200825_024831/0000/',
        nJobs     = 50,
        suffix    = 'zg_ewk_2017'
        ),
    
    cfg(data_name = 'TTJets_2017_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0000/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),
    cfg(data_name = 'TTJets_2017_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0001/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),
    cfg(data_name = 'TTJets_2017_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0002/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),
    cfg(data_name = 'TTJets_2017_part3', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0003/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),
    cfg(data_name = 'TTJets_2017_part4', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0004/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),
    cfg(data_name = 'TTJets_2017_part5', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0005/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),
    cfg(data_name = 'TTJets_2017_part6', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_prod_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200929_163736/0006/',
        nJobs     = 25,
        suffix    = 'ttjets_2017'
        ),

    cfg(data_name = 'hzg_gluglu_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/GluGluHToZG_ZToLL_M-120_13TeV_powheg_pythia8/2017_mc_prod_hzg_GluGluHToZG_ZToLL_M-120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_reali/200814_184900/0000/',
        nJobs     = 25,
        suffix    = 'hzg_gluglu_M120_2017'
        ),
    cfg(data_name = 'hzg_vbf_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/VBFHToZG_ZToLL_M-120_13TeV_powheg_pythia8/2017_mc_prod_hzg_VBFHToZG_ZToLL_M-120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200814_185042/0000/',
        nJobs     = 18,
        suffix    = 'hzg_vbf_M120_2017'
        ),
    cfg(data_name = 'hzg_wplush_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8/2017_mc_prod_hzg_WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc201/200814_185713/0000/',
        nJobs     = 12,
        suffix    = 'hzg_wplush_M120_2017'
        ),
    cfg(data_name = 'hzg_wminush_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8/2017_mc_prod_hzg_WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc20/200814_185905/0000/',
        nJobs     = 9,
        suffix    = 'hzg_wminush_M120_2017'
        ),
    cfg(data_name = 'hzg_zh_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_pythia8/2017_mc_prod_hzg_ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_re/200814_185425/0000/',
        nJobs     = 10,
        suffix    = 'hzg_zh_M120_2017'
        ),
    cfg(data_name = 'hzg_tth_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ttHToZG_ZToLL_M120_13TeV_powheg_pythia8/2017_mc_prod_hzg_ttHToZG_ZToLL_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_/200814_185241/0000/',
        nJobs     = 17, 
        suffix    = 'hzg_tth_M120_2017'
        ),

    cfg(data_name = 'hzg_gluglu_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/GluGluHToZG_ZToLL_M-125_13TeV_powheg_pythia8/2017_mc_prod_hzg_GluGluHToZG_ZToLL_M-125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_reali/200814_183718/0000/',
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M125_2017'
        ),
    cfg(data_name = 'hzg_vbf_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/VBFHToZG_ZToLL_M-125_13TeV_powheg_pythia8/2017_mc_prod_hzg_VBFHToZG_ZToLL_M-125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200814_183903/0000/',
        nJobs     = 42,
        suffix    = 'hzg_vbf_M125_2017'
        ),
    cfg(data_name = 'hzg_wplush_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8/2017_mc_prod_hzg_WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc201/200814_184444/0000/',
        nJobs     = 11,
        suffix    = 'hzg_wplush_M125_2017'
        ),
    cfg(data_name = 'hzg_wminush_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8/2017_mc_prod_hzg_WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc20/200814_184632/0000/',
        nJobs     = 11,
        suffix    = 'hzg_wminush_M125_2017'
        ),
    cfg(data_name = 'hzg_zh_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_pythia8/2017_mc_prod_hzg_ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_re/200814_184234/0000/',
        nJobs     = 11,
        suffix    = 'hzg_zh_M125_2017'
        ),
    cfg(data_name = 'hzg_tth_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ttHToZG_ZToLL_M125_13TeV_powheg_pythia8/2017_mc_prod_hzg_ttHToZG_ZToLL_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_/200814_184050/0000/',
        nJobs     = 9, 
        suffix    = 'hzg_tth_M125_2017'
        ),
    
    cfg(data_name = 'hzg_gluglu_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/GluGluHToZG_ZToLL_M-130_13TeV_powheg_pythia8/2017_mc_prod_hzg_GluGluHToZG_ZToLL_M-130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_reali/200814_190106/0000/',
        nJobs     = 26,
        suffix    = 'hzg_gluglu_M130_2017'
        ),
    cfg(data_name = 'hzg_vbf_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/VBFHToZG_ZToLL_M-130_13TeV_powheg_pythia8/2017_mc_prod_hzg_VBFHToZG_ZToLL_M-130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/200814_190252/0000/',
        nJobs     = 21,
        suffix    = 'hzg_vbf_M130_2017'
        ),
    cfg(data_name = 'hzg_wplush_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8/2017_mc_prod_hzg_WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc201/200814_190900/0000/',
        nJobs     = 9,
        suffix    = 'hzg_wplush_M130_2017'
        ),
    cfg(data_name = 'hzg_wminush_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8/2017_mc_prod_hzg_WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc20/200814_191114/0000/',
        nJobs     = 8,
        suffix    = 'hzg_wminush_M130_2017'
        ),
    cfg(data_name = 'hzg_zh_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_pythia8/2017_mc_prod_hzg_ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_re/200814_190714/0000/',
        nJobs     = 12,
        suffix    = 'hzg_zh_M130_2017'
        ),
    cfg(data_name = 'hzg_tth_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ttHToZG_ZToLL_M130_13TeV_powheg_pythia8/2017_mc_prod_hzg_ttHToZG_ZToLL_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_/200814_190437/0000/',
        nJobs     = 20, 
        suffix    = 'hzg_tth_M130_2017'
        ),

    cfg(data_name = 'hmumu_gluglu_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/GluGluHToMuMu_M120_13TeV_amcatnloFXFX_pythia8/2017_mc_prod_hmumu_GluGluHToMuMu_M120_13TeV_amcatnloFXFX_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_real/200814_204026/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_gluglu_M120_2017'
        ),
    cfg(data_name = 'hmumu_vbf_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/VBFHToMuMu_M120_13TeV_amcatnlo_pythia8/2017_mc_prod_hmumu_VBFHToMuMu_M120_13TeV_amcatnlo_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v/200814_204216/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_vbf_M120_2017'
        ),
    cfg(data_name = 'hmumu_wplush_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WplusH_HToMuMu_WToAll_M120_13TeV_powheg_pythia8/2017_mc_prod_hmumu_WplusH_HToMuMu_WToAll_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_re/200825_151047/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_wplush_M120_2017'
        ),
    cfg(data_name = 'hmumu_wminush_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WminusH_HToMuMu_WToAll_M120_13TeV_powheg_pythia8/2017_mc_prod_hmumu_WminusH_HToMuMu_WToAll_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_r/200814_204948/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_wminush_M120_2017'
        ),
    cfg(data_name = 'hmumu_zh_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ZH_HToMuMu_ZToAll_M120_13TeV_powheg_pythia8/2017_mc_prod_hmumu_ZH_HToMuMu_ZToAll_M120_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realis/200814_204618/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_zh_M120_2017'
        ),
    cfg(data_name = 'hmumu_tth_M120_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ttHToMuMu_M120_TuneCP5_PSweights_13TeV-powheg-pythia8/2017_mc_prod_hmumu_ttHToMuMu_M120_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2/200814_204426/0000/',
        nJobs     = 50, 
        suffix    = 'hmumu_tth_M120_2017'
        ),

    cfg(data_name = 'hmumu_gluglu_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/GluGluHToMuMu_M125_13TeV_amcatnloFXFX_pythia8/2017_mc_prod_hmumu_GluGluHToMuMu_M125_13TeV_amcatnloFXFX_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_real/200814_202713/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_gluglu_M125_2017'
        ),
    cfg(data_name = 'hmumu_vbf_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/VBFHToMuMu_M125_13TeV_amcatnlo_pythia8/2017_mc_prod_hmumu_VBFHToMuMu_M125_13TeV_amcatnlo_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v/200814_202916/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_vbf_M125_2017'
        ),
    cfg(data_name = 'hmumu_wplush_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WplusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8/2017_mc_prod_hmumu_WplusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_re/200814_203611/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_wplush_M125_2017'
        ),
    cfg(data_name = 'hmumu_wminush_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WminusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8/2017_mc_prod_hmumu_WminusH_HToMuMu_WToAll_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_r/200814_203824/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_wminush_M125_2017'
        ),
    cfg(data_name = 'hmumu_zh_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8/2017_mc_prod_hmumu_ZH_HToMuMu_ZToAll_M125_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realis/200814_203349/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_zh_M125_2017'
        ),
    cfg(data_name = 'hmumu_tth_M125_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/2017_mc_prod_hmumuttH_ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2/200827_203453/0000/',
        nJobs     = 23, 
        suffix    = 'hmumu_tth_M125_2017'
        ),
    
    cfg(data_name = 'hmumu_gluglu_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/GluGluHToMuMu_M130_13TeV_amcatnloFXFX_pythia8/2017_mc_prod_hmumu_GluGluHToMuMu_M130_13TeV_amcatnloFXFX_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_real/200814_205153/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_gluglu_M130_2017'
        ),
    cfg(data_name = 'hmumu_vbf_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/VBFHToMuMu_M130_13TeV_amcatnlo_pythia8/2017_mc_prod_hmumu_VBFHToMuMu_M130_13TeV_amcatnlo_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v/200814_205339/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_vbf_M130_2017'
        ),
    cfg(data_name = 'hmumu_wplush_M130_2017',
        path      = '/eos/uscms/store/user/mmackenz/WplusH_HToMuMu_WToAll_M130_13TeV_powheg_pythia8/2017_mc_prod_hmumu_WplusH_HToMuMu_WToAll_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_re/',
        nJobs     = 1,
        suffix    = 'hmumu_wplush_M130_2017'
        ),
    cfg(data_name = 'hmumu_wminush_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/WminusH_HToMuMu_WToAll_M130_13TeV_powheg_pythia8/2017_mc_prod_hmumu_WminusH_HToMuMu_WToAll_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_r/200814_210403/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_wminush_M130_2017'
        ),
    cfg(data_name = 'hmumu_zh_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ZH_HToMuMu_ZToAll_M130_13TeV_powheg_pythia8/2017_mc_prod_hmumu_ZH_HToMuMu_ZToAll_M130_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realis/200814_205958/0000/',
        nJobs     = 50,
        suffix    = 'hmumu_zh_M130_2017'
        ),
    cfg(data_name = 'hmumu_tth_M130_2017',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/ttHToMuMu_M130_TuneCP5_PSweights_13TeV-powheg-pythia8/2017_mc_prod_hmumu_ttHToMuMu_M130_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2/200814_205632/0000/',
        nJobs     = 50, 
        suffix    = 'hmumu_tth_M130_2017'
        ),
    ])


batch_17_samples = []
if selection == 'mmg':
    batch_17_samples += muon_data_17
elif selection == 'eeg':
    batch_17_samples += electron_data_17
batch_17_samples += mc_17

batch_17 = bm.BatchMaster(config_list = batch_17_samples, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = 2017,
                       executable  = executable,
                       location    = 'lpc'
                     )

muon_data_18 = []
muon_data_18.extend([
    
    cfg(data_name = 'muon_2018A_part0',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part0/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_part1',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part1/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_part2',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part2/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_part3',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part3/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_part4',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part4/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_part5',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part5/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_part6',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_part6/',
        nJobs     = 25,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018A_recovery',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018A-17Sep2018-v2_recovery/',
        nJobs     = 10,
        suffix    = 'muon_2018A'
       ),

    cfg(data_name = 'muon_2018B_part0',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018B-17Sep2018-v1_part0/',
        nJobs     = 25,
        suffix    = 'muon_2018B'
       ),
    cfg(data_name = 'muon_2018B_part1',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018B-17Sep2018-v1_part1/',
        nJobs     = 25,
        suffix    = 'muon_2018B'
       ),
    cfg(data_name = 'muon_2018B_part2',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018B-17Sep2018-v1_part2/',
        nJobs     = 25,
        suffix    = 'muon_2018B'
       ),
    cfg(data_name = 'muon_2018B_recovery',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018B-17Sep2018-v1_recovery/',
        nJobs     = 16,
        suffix    = 'muon_2018B'
       ),

    cfg(data_name = 'muon_2018C',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018C-17Sep2018-v1/',
        nJobs     = 50,
        suffix    = 'muon_2018C'
       ),

    cfg(data_name = 'muon_2018D_part0',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part0/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part1',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part1/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part2',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part2/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part3',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part3/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part4',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part4/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part5',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part5/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part6',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part6/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part7',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part7/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part8',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part8/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    cfg(data_name = 'muon_2018D_part9',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/DoubleMuon_Run2018D-PromptReco-v2_part9/',
        nJobs     = 25,
        suffix    = 'muon_2018D'
       ),
    ])

electron_data_18 = []
electron_data_18.extend([

    cfg(data_name = 'electron_2018A_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0000/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0001/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0002/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part3', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0003/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part4', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0004/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part5', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0005/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part6', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018A-17Sep2018-v2/200902_141330/0006/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_recovery', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_A_recovery_EGamma_Run2018A-17Sep2018-v2/200923_214917/0000/',
        nJobs     = 25,
        suffix    = 'electron_2018A'
        ),

    cfg(data_name = 'electron_2018B_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0000/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0001/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0002/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part3', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0003/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part4', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0004/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part5', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0005/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part6', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0006/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part7', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0007/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part8', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0008/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018B_part9', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_B_retry_EGamma_Run2018B-17Sep2018-v1/201001_012459/0009/',
        nJobs     = 25,
        suffix    = 'electron_2018B'
        ),

    cfg(data_name = 'electron_2018C_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018C-17Sep2018-v1/200902_141533/0000/',
        nJobs     = 25,
        suffix    = 'electron_2018C'
        ),
    cfg(data_name = 'electron_2018C_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018C-17Sep2018-v1/200902_141533/0001/',
        nJobs     = 25,
        suffix    = 'electron_2018C'
        ),
    cfg(data_name = 'electron_2018C_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/EGamma/2018_dataABC_prod_aug_20_EGamma_AC_retry_EGamma_Run2018C-17Sep2018-v1/200902_141533/0002/',
        nJobs     = 25,
        suffix    = 'electron_2018C'
        ),

    cfg(data_name = 'electron_2018D_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part0/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part1/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part2/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part3', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part3/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part4', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part4/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part5', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part5/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part6', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma_Run2018D-22Jan2019-v2_part6/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0000/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0001/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0002/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part3', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0003/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part4', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0004/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part5', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0005/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part6', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0006/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part7', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_EGamma_Run2018D-22Jan2019-v2/200902_162204/0007/',
        nJobs     = 25,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_recovery_part8', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/EGamma/2018_dataD_prod_aug_20_EGammaD_recovery_part2_EGamma_Run2018D-22Jan2019-v2/200923_212059/0000/',
        nJobs     = 5,
        suffix    = 'electron_2018D'
        ),
    ])

mc_18 = []
mc_18.extend([
    cfg(data_name = 'DYJetsToLL_M50_2018_part0',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_/200817_152906/0000/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2018'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2018_part1',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_/200817_152906/0001/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2018'
        ),
    cfg(data_name = 'DYJetsToLL_M50_2018_part2',
        path      = '/eos/uscms/store/group/lpcbacon/mmackenz/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_/200817_152906/0002/',
        nJobs     = 25,
        suffix    = 'zjets_M50_2018'
        ),

    cfg(data_name = 'ZGToLLG_2018_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/naodell/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018/200815_001048/0000/',
        nJobs     = 25,
        suffix    = 'zg_llg_2018'
        ),
    cfg(data_name = 'ZGToLLG_2018_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/naodell/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018/200815_001048/0001/',
        nJobs     = 25,
        suffix    = 'zg_llg_2018'
        ),

    cfg(data_name = 'LLAJJ_EWK_M50_2018', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/LLAJJ_EWK_MLL-50_MJJ-120_TuneCP5_13TeV-madgraph-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 50,
        suffix    = 'zg_ewk_2018'
        ),

    cfg(data_name = 'TTJets_2018_part0', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_aug_20_TTJets_retry_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-/200902_164241/0000/',
        nJobs     = 25,
        suffix    = 'ttjets_2018'
        ),
    cfg(data_name = 'TTJets_2018_part1', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_aug_20_TTJets_retry_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-/200902_164241/0001/',
        nJobs     = 25,
        suffix    = 'ttjets_2018'
        ),
    cfg(data_name = 'TTJets_2018_part2', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_aug_20_TTJets_retry_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-/200902_164241/0002/',
        nJobs     = 25,
        suffix    = 'ttjets_2018'
        ),
    cfg(data_name = 'TTJets_2018_part3', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_aug_20_TTJets_retry_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-/200902_164241/0003/',
        nJobs     = 25,
        suffix    = 'ttjets_2018'
        ),
    cfg(data_name = 'TTJets_2018_part4', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_aug_20_TTJets_retry_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-/200902_164241/0004/',
        nJobs     = 25,
        suffix    = 'ttjets_2018'
        ),
    cfg(data_name = 'TTJets_2018_part5', 
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/2018_mc_prod_aug_20_TTJets_retry_TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-/200902_164241/0005/',
        nJobs     = 25,
        suffix    = 'ttjets_2018'
        ),

    cfg(data_name = 'hzg_gluglu_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/GluGluHToZG_M-120_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 20,
        suffix    = 'hzg_gluglu_M120_2018'
        ),
    cfg(data_name = 'hzg_vbf_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/VBFHToZG_M-120_TuneCP5_13TeV-powheg-pythia8/2018_mc_prod_aug_20_hzg_signal_retry_VBFHToZG_M-120_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/200914_143831/0000/',
        nJobs     = 20,
        suffix    = 'hzg_vbf_M120_2018'
        ),
    cfg(data_name = 'hzg_wplush_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 9,
        suffix    = 'hzg_wplush_M120_2018'
        ),
    cfg(data_name = 'hzg_wminush_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WminusH_HToZG_WToAll_ZToLL_M120_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 8,
        suffix    = 'hzg_wminush_M120_2018'
        ),
    cfg(data_name = 'hzg_zh_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZH_ZToAll_HToZG_ZToLL_TuneCP5_M120_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 15,
        suffix    = 'hzg_zh_M120_2018'
        ),
    cfg(data_name = 'hzg_tth_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/ttHToZG_ZToLL_M120_13TeV_powheg_pythia8/2018_mc_prod_aug_20_ttHToZG_ZToLL_M120_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/200821_173411/0000/',
        nJobs     = 25, 
        suffix    = 'hzg_tth_M120_2018'
        ),

    cfg(data_name = 'hzg_gluglu_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/GluGluHToZG_M-125_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 10,
        suffix    = 'hzg_gluglu_M125_2018'
        ),
    cfg(data_name = 'hzg_vbf_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/VBFHToZG_M-125_TuneCP5_13TeV-powheg-pythia8/2018_mc_prod_aug_20_hzg_vbf_pfpart_fix_VBFHToZG_M-125_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/201204_231549/0000/',
        nJobs     = 15,
        suffix    = 'hzg_vbf_M125_2018'
        ),
    cfg(data_name = 'hzg_wplush_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8/2018_mc_prod_aug_20_wplush_hzg_125_retry_WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic/200903_155731/0000/',
        nJobs     = 2,
        suffix    = 'hzg_wplush_M125_2018'
        ),
    cfg(data_name = 'hzg_wminush_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8/2018_mc_prod_aug_20_WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realisti/200821_172918/0000/',
        nJobs     = 8,
        suffix    = 'hzg_wminush_M125_2018'
        ),
    cfg(data_name = 'hzg_zh_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZH_ZToAll_HToZG_ZToLL_TuneCP5_M125_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 8,
        suffix    = 'hzg_zh_M125_2018'
        ),
    cfg(data_name = 'hzg_tth_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/ttHToZG_ZToLL_M125_13TeV_powheg_pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_ttHToZG_ZToLL_M125_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/200916_152359/0000/',
        nJobs     = 5, 
        suffix    = 'hzg_tth_M125_2018'
        ),
    
    cfg(data_name = 'hzg_gluglu_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/GluGluHToZG_M-130_TuneCP5_13TeV-powheg-pythia8/2018_mc_prod_aug_20_hzg_signal_retry_GluGluHToZG_M-130_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-/200914_144005/0000/',
        nJobs     = 19,
        suffix    = 'hzg_gluglu_M130_2018'
        ),
    cfg(data_name = 'hzg_vbf_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/VBFHToZG_M-130_TuneCP5_13TeV-powheg-pythia8/2018_mc_prod_aug_20_hzg_signal_retry_VBFHToZG_M-130_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/200914_144140/0000/',
        nJobs     = 15,
        suffix    = 'hzg_vbf_M130_2018'
        ),
    cfg(data_name = 'hzg_wplush_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 8,
        suffix    = 'hzg_wplush_M130_2018'
        ),
    cfg(data_name = 'hzg_wminush_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 11,
        suffix    = 'hzg_wminush_M130_2018'
        ),
    cfg(data_name = 'hzg_zh_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZH_ZToAll_HToZG_ZToLL_M130_TuneCP5_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 12,
        suffix    = 'hzg_zh_M130_2018'
        ),
    cfg(data_name = 'hzg_tth_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/ttHToZG_ZToLL_M130_13TeV_powheg_pythia8/2018_mc_prod_aug_20_ttHToZG_ZToLL_M130_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/200821_174343/0000/',
        nJobs     = 25, 
        suffix    = 'hzg_tth_M130_2018'
        ),

    cfg(data_name = 'hmumu_gluglu_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/GluGluHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_GluGluHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIIAutumn18MiniAOD-102X_upgrade20/200914_145913/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_gluglu_M120_2018'
        ),
    cfg(data_name = 'hmumu_vbf_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/VBFHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_VBFHToMuMu_M120_TuneCP5_PSweights_13TeV_amcatnlo_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_real/200914_150055/0000/',
        nJobs     = 20,
        suffix    = 'hmumu_vbf_M120_2018'
        ),
    cfg(data_name = 'hmumu_wplush_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/WplusH_HToMuMu_WToAll_M120_TuneCP5_PSweights_13TeV_powheg_pythia8/2018_mc_prod_aug_20_WplusH_HToMuMu_WToAll_M120_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade/200821_180649/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_wplush_M120_2018'
        ),
    cfg(data_name = 'hmumu_wminush_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WminusH_HToMuMu_WToAll_M120_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_wminush_M120_2018'
        ),
    cfg(data_name = 'hmumu_zh_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZH_HToMuMu_ZToAll_M120_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_zh_M120_2018'
        ),
    cfg(data_name = 'hmumu_tth_M120_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/ttHToMuMu_M120_TuneCP5_PSweights_13TeV-powheg-pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_ttHToMuMu_M120_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realist/200914_150232/0000/',
        nJobs     = 22, 
        suffix    = 'hmumu_tth_M120_2018'
        ),

    cfg(data_name = 'hmumu_gluglu_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/2018_mc_prod_aug_20_GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIIAutumn18MiniAOD-102X_upgrade20/200821_175027/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_gluglu_M125_2018'
        ),
    cfg(data_name = 'hmumu_vbf_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_VBFHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnlo_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_real/200914_145417/0000/',
        nJobs     = 22,
        suffix    = 'hmumu_vbf_M125_2018'
        ),
    cfg(data_name = 'hmumu_wplush_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/WplusH_HToMuMu_WToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_WplusH_HToMuMu_WToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD-102X_upgrade/200914_145735/0000/',
        nJobs     = 23,
        suffix    = 'hmumu_wplush_M125_2018'
        ),
    cfg(data_name = 'hmumu_wminush_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WminusH_HToMuMu_WToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_wminush_M125_2018'
        ),
    cfg(data_name = 'hmumu_zh_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZH_HToMuMu_ZToAll_M125_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_zh_M125_2018'
        ),
    cfg(data_name = 'hmumu_tth_M125_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_ttHToMuMu_M125_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realist/200914_145558/0000/',
        nJobs     = 23, 
        suffix    = 'hmumu_tth_M125_2018'
        ),
    
    cfg(data_name = 'hmumu_gluglu_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/GluGluHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8/2018_mc_prod_aug_20_GluGluHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIIAutumn18MiniAOD-102X_upgrade20/200821_181003/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_gluglu_M130_2018'
        ),
    cfg(data_name = 'hmumu_vbf_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/VBFHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnlo_pythia8/2018_mc_prod_aug_20_hmumu_signal_retry_VBFHToMuMu_M130_TuneCP5_PSweights_13TeV_amcatnlo_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_real/200914_150408/0000/',
        nJobs     = 25,
        suffix    = 'hmumu_vbf_M130_2018'
        ),
    cfg(data_name = 'hmumu_wplush_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WplusH_HToMuMu_WToAll_M130_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_wplush_M130_2018'
        ),
    cfg(data_name = 'hmumu_wminush_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/WminusH_HToMuMu_WToAll_M130_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_wminush_M130_2018'
        ),
    cfg(data_name = 'hmumu_zh_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/aug_20/ZH_HToMuMu_ZToAll_M130_TuneCP5_PSweights_13TeV_powheg_pythia8_RunIIAutumn18MiniAOD/',
        nJobs     = 25,
        suffix    = 'hmumu_zh_M130_2018'
        ),
    cfg(data_name = 'hmumu_tth_M130_2018',
        path      = '/eos/uscms/store/group/lpcbacon/jbueghly/ttHToMuMu_M130_TuneCP5_PSweights_13TeV-powheg-pythia8/2018_mc_prod_aug_20_ttHToMuMu_M130_TuneCP5_PSweights_13TeV-powheg-pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realist/200821_181319/0000/',
        nJobs     = 25, 
        suffix    = 'hmumu_tth_M130_2018'
        ),
    ])

batch_18_samples = []
if selection == 'mmg':
    batch_18_samples += muon_data_18
elif selection == 'eeg':
    batch_18_samples += electron_data_18
batch_18_samples += mc_18

batch_18 = bm.BatchMaster(config_list = batch_18_samples, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = 2018,
                       executable  = executable,
                       location    = 'lpc'
                     )

if 2016 in periods:
    batch_16.submit_to_batch()

if 2017 in periods:
    batch_17.submit_to_batch()

if 2018 in periods:
    batch_18.submit_to_batch()

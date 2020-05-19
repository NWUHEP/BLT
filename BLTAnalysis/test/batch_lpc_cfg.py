#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys

''' Specify parameters '''
cfg        = bm.JobConfig
selection = 'eeg'
periods = [2017]
path        = '/eos/uscms/store/group/lpcbacon/jbueghly/'
executable = 'execBatch.sh'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

muon_data_16 = []
muon_data_16.extend([

    # Double muon data
    cfg(data_name = 'muon_2016B',
        path     = '{0}/DoubleMuon_Run2016B-17Jul2018_ver2-v1_tmp'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016B'
       ),
    cfg(data_name = 'muon_2016C',
        path     = '{0}/DoubleMuon_Run2016C-17Jul2018-v1_tmp'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016C'
       ),
    cfg(data_name = 'muon_2016D',
        path     = '{0}/DoubleMuon_Run2016D-17Jul2018-v1_tmp'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016D'
       ),
    cfg(data_name = 'muon_2016E',
        path     = '{0}/DoubleMuon_Run2016E-17Jul2018-v1_tmp'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016E'
       ),
    cfg(data_name = 'muon_2016F',
        path     = '{0}/DoubleMuon_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016F'
       ),
    cfg(data_name = 'muon_2016G',
        path     = '{0}/DoubleMuon_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016G'
       ),
    cfg(data_name = 'muon_2016H',
        path     = '{0}/DoubleMuon_Run2016H-17Jul2018-v1_tmp'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016H'
       ),
    ])

electron_data_16 = []
electron_data_16.extend([

    # Double electron data
    cfg(data_name = 'electron_2016B',
        path     = '{0}/DoubleEG_Run2016B-17Jul2018_ver2-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016B'
       ),
    cfg(data_name = 'electron_2016C',
        path     = '{0}/DoubleEG_Run2016C-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016C'
       ),
    cfg(data_name = 'electron_2016D',
        path     = '{0}/DoubleEG_Run2016D-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016D'
       ),
    cfg(data_name = 'electron_2016E',
        path     = '{0}/DoubleEG_Run2016E-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016E'
       ),
    cfg(data_name = 'electron_2016F',
        path     = '{0}/DoubleEG_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016F'
       ),
    cfg(data_name = 'electron_2016G',
        path     = '{0}/DoubleEG_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016G'
       ),
    cfg(data_name = 'electron_2016H',
        path     = '{0}/DoubleEG_Run2016H-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016H'
       ),
    ])
 
mc_16 = []
mc_16.extend([

    cfg(data_name = 'DYJetsToLL_M50_2016', 
        path      = '{0}/DYJetsToLL_M-50_amcatnlo_RunIISummer16MiniAODv3_tmp'.format(path),
        nJobs     = 50,
        suffix    = 'zjets_M50_2016'
        ),
    cfg(data_name = 'ZGToLLG_2016', 
        path      = '{0}/ZGToLLG_01J_LoosePtlPtg_5f_13TeV-amcatnloFXFX_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50,
        suffix    = 'zg_llg_2016'
        ),

    cfg(data_name = 'hzg_gluglu_M125_2016',
        path      = '{0}/GluGluHToZG_M-125_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M125_2016'
        ),
    cfg(data_name = 'hzg_vbf_M125_2016',
        path      = '{0}/VBFHToZG_M-125_13TeV_powheg_RunIISummer16MiniAODv3_v3_ext'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M125_2016'
        ),
    cfg(data_name = 'hzg_wplush_M125_2016',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M125_2016'
        ),
    cfg(data_name = 'hzg_wminush_M125_2016',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M125_2016'
        ),
    cfg(data_name = 'hzg_zh_M125_2016',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M125_2016'
        ),
    cfg(data_name = 'hzg_tth_M125_2016',
        path      = '{0}/ttHToZG_ZToLL_M125_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M125_2016'
        ),
    
    cfg(data_name = 'hzg_gluglu_M120_2016',
        path      = '{0}/GluGluHToZG_M-120_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 4,
        suffix    = 'hzg_gluglu_M120_2016'
        ),
    cfg(data_name = 'hzg_vbf_M120_2016',
        path      = '{0}/VBFHToZG_M-120_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M120_2016'
        ),
    cfg(data_name = 'hzg_wplush_M120_2016',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 1,
        suffix    = 'hzg_wplush_M120_2016'
        ),
    cfg(data_name = 'hzg_wminush_M120_2016',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 1,
        suffix    = 'hzg_wminush_M120_2016'
        ),
    cfg(data_name = 'hzg_zh_M120_2016',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 2,
        suffix    = 'hzg_zh_M120_2016'
        ),
    cfg(data_name = 'hzg_tth_M120_2016',
        path      = '{0}/ttHToZG_ZToLL_M120_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M120_2016'
        ),
    
    cfg(data_name = 'hzg_gluglu_M130_2016',
        path      = '{0}/GluGluHToZG_M-130_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 5,
        suffix    = 'hzg_gluglu_M130_2016'
        ),
    cfg(data_name = 'hzg_vbf_M130_2016',
        path      = '{0}/VBFHToZG_M-130_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 3,
        suffix    = 'hzg_vbf_M130_2016'
        ),
    cfg(data_name = 'hzg_wplush_M130_2016',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 2,
        suffix    = 'hzg_wplush_M130_2016'
        ),
    cfg(data_name = 'hzg_wminush_M130_2016',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 1,
        suffix    = 'hzg_wminush_M130_2016'
        ),
    cfg(data_name = 'hzg_zh_M130_2016',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 1,
        suffix    = 'hzg_zh_M130_2016'
        ),
    cfg(data_name = 'hzg_tth_M130_2016',
        path      = '{0}/ttHToZG_ZToLL_M130_13TeV_powheg_RunIISummer16MiniAODv3'.format(path),
        nJobs     = 1, 
        suffix    = 'hzg_tth_M130_2016'
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
        path     = '{0}/DoubleMuon_Run2017B-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2017B'
       ),
    cfg(data_name = 'muon_2017C',
        path     = '{0}/DoubleMuon_Run2017C-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2017C'
       ),
    cfg(data_name = 'muon_2017D',
        path     = '{0}/DoubleMuon_Run2017D-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2017D'
       ),
    cfg(data_name = 'muon_2017E',
        path     = '{0}/DoubleMuon_Run2017E-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2017E'
       ),
    cfg(data_name = 'muon_2017F',
        path     = '{0}/DoubleMuon_Run2017F-31Mar2018-v1/'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2017F'
       ),
    ])

electron_data_17 = []
electron_data_17.extend([
    
    # DoubleEG 2017
    cfg(data_name = 'electron_2017B',
        path     = '{0}/DoubleEG_Run2017B-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2017B'
       ),
    cfg(data_name = 'electron_2017C',
        path     = '{0}/DoubleEG_Run2017C-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2017C'
       ),
    cfg(data_name = 'electron_2017D',
        path     = '{0}/DoubleEG_Run2017D-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2017D'
       ),
    cfg(data_name = 'electron_2017E',
        path     = '{0}/DoubleEG_Run2017E-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2017E'
       ),
    cfg(data_name = 'electron_2017F',
        path     = '{0}/DoubleEG_Run2017F-31Mar2018-v1/'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2017F'
       ),
   ])

mc_17 = []
mc_17.extend([
    cfg(data_name = 'DYJetsToLL_M50_2017',
        path      = '{0}/DYJetsToLL_M-50_amcatnlo_Fall17_v14_ext1-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'zjets_M50_2017'
        ),
    cfg(data_name = 'ZGToLLG_2017', 
        path      = '{0}/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX_Fall17_tmp/'.format(path),
        nJobs     = 50,
        suffix    = 'zg_llg_2017'
        ),

    cfg(data_name = 'hzg_gluglu_M120_2017',
        #path      = '{0}/GluGluHToZG_ZToLL_M-120_13TeV_powheg_Fall17_v14_ext1-v1/'.format(path),
        path      = '{0}/GluGluHToZG_ZToLL_M-120_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M120_2017'
        ),
    cfg(data_name = 'hzg_vbf_M120_2017',
        #path      = '{0}/VBFHToZG_ZToLL_M-120_13TeV_powheg_Fall17_v14_ext1-v1/'.format(path),
        path      = '{0}/VBFHToZG_ZToLL_M-120_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M120_2017'
        ),
    cfg(data_name = 'hzg_wplush_M120_2017',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M120_2017'
        ),
    cfg(data_name = 'hzg_wminush_M120_2017',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M120_2017'
        ),
    cfg(data_name = 'hzg_zh_M120_2017',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M120_2017'
        ),
    cfg(data_name = 'hzg_tth_M120_2017',
        path      = '{0}/ttHToZG_ZToLL_M120_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M120_2017'
        ),

    cfg(data_name = 'hzg_gluglu_M125_2017',
        path      = '{0}/GluGluHToZG_ZToLL_M-125_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M125_2017'
        ),
    cfg(data_name = 'hzg_vbf_M125_2017',
        path      = '{0}/VBFHToZG_ZToLL_M-125_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M125_2017'
        ),
    cfg(data_name = 'hzg_wplush_M125_2017',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M125_2017'
        ),
    cfg(data_name = 'hzg_wminush_M125_2017',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M125_2017'
        ),
    cfg(data_name = 'hzg_zh_M125_2017',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M125_2017'
        ),
    cfg(data_name = 'hzg_tth_M125_2017',
        path      = '{0}/ttHToZG_ZToLL_M125_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M125_2017'
        ),
    
    cfg(data_name = 'hzg_gluglu_M130_2017',
        path      = '{0}/GluGluHToZG_ZToLL_M-130_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M130_2017'
        ),
    cfg(data_name = 'hzg_vbf_M130_2017',
        path      = '{0}/VBFHToZG_ZToLL_M-130_13TeV_powheg_Fall17_rerun/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M130_2017'
        ),
    cfg(data_name = 'hzg_wplush_M130_2017',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M130_2017'
        ),
    cfg(data_name = 'hzg_wminush_M130_2017',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M130_2017'
        ),
    cfg(data_name = 'hzg_zh_M130_2017',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M130_2017'
        ),
    cfg(data_name = 'hzg_tth_M130_2017',
        path      = '{0}/ttHToZG_ZToLL_M130_13TeV_powheg_Fall17_v14-v1/'.format(path),
        nJobs     = 5, 
        suffix    = 'hzg_tth_M130_2017'
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
    
    cfg(data_name = 'muon_2018A',
        path      = '{0}/DoubleMuon_Run2018A-17Sep2018-v2'.format(path),
        nJobs     = 50,
        suffix    = 'muon_2018A'
       ),
    cfg(data_name = 'muon_2018B',
        path      = '{0}/DoubleMuon_Run2018B-17Sep2018-v1'.format(path),
        nJobs     = 50,
        suffix    = 'muon_2018B'
       ),
    cfg(data_name = 'muon_2018C',
        path      = '{0}/DoubleMuon_Run2018C-17Sep2018-v1'.format(path),
        nJobs     = 50,
        suffix    = 'muon_2018C'
       ),
    cfg(data_name = 'muon_2018D',
        path      = '{0}/DoubleMuon_Run2018D-PromptReco-v2'.format(path),
        nJobs     = 50,
        suffix    = 'muon_2018D'
       ),
    ])

electron_data_18 = []
electron_data_18.extend([

    cfg(data_name = 'electron_2018A_part1', 
        path      = '{0}/EGamma_Run2018A-17Sep2018-v2_part1_tmp'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part2', 
        path      = '{0}/EGamma_Run2018A-17Sep2018-v2_part2'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part3', 
        path      = '{0}/EGamma_Run2018A-17Sep2018-v2_part3'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018A_part4', 
        path      = '{0}/EGamma_Run2018A-17Sep2018-v2_part4'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018A'
        ),
    cfg(data_name = 'electron_2018B', 
        path      = '{0}/EGamma/2018_data_EG_B_EGamma_Run2018B-17Sep2018-v1/190816_130108/tmp'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018B'
        ),
    cfg(data_name = 'electron_2018C', 
        path      = '{0}/EGamma_Run2018C-17Sep2018-v1'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018C'
        ),
    cfg(data_name = 'electron_2018D_part1', 
        path      = '{0}/EGamma_Run2018D-22Jan2019-v2_part1_tmp'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part2', 
        path      = '{0}/EGamma_Run2018D-22Jan2019-v2_part2_tmp'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part3', 
        path      = '{0}/EGamma_Run2018D-22Jan2019-v2_part3_tmp'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018D'
        ),
    cfg(data_name = 'electron_2018D_part4', 
        path      = '{0}/EGamma_Run2018D-22Jan2019-v2_part4'.format(path),
        nJobs     = 50,
        suffix    = 'electron_2018D'
        ),
    ])

mc_18 = []
mc_18.extend([
    cfg(data_name = 'DYJetsToLL_M50_2018',
        path      = '{0}/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX_RunIIAutumn18MiniAOD_tmp/'.format(path),
        nJobs     = 50,
        suffix    = 'zjets_M50_2018'
        ),
    cfg(data_name = 'ZGToLLG_2018', 
        path      = '{0}/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX_RunIIAutumn18MiniAOD_tmp/'.format(path),
        nJobs     = 50,
        suffix    = 'zg_llg_2018'
        ),

    cfg(data_name = 'hzg_gluglu_M120_2018',
        path      = '{0}/GluGluHToZG_M-120_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M120_2018'
        ),
    cfg(data_name = 'hzg_vbf_M120_2018',
        path      = '{0}/VBFHToZG_M-120_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M120_2018'
        ),
    cfg(data_name = 'hzg_wplush_M120_2018',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M120_2018'
        ),
    cfg(data_name = 'hzg_wminush_M120_2018',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M120_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M120_2018'
        ),
    cfg(data_name = 'hzg_zh_M120_2018',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M120_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M120_2018'
        ),
    cfg(data_name = 'hzg_tth_M120_2018',
        path      = '{0}/ttHToZG_ZToLL_M120_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M120_2018'
        ),
    
    cfg(data_name = 'hzg_gluglu_M125_2018',
        path      = '{0}/GluGluHToZG_M-125_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M125_2018'
        ),
    cfg(data_name = 'hzg_vbf_M125_2018',
        path      = '{0}/VBFHToZG_M-125_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M125_2018'
        ),
    cfg(data_name = 'hzg_wplush_M125_2018',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M125_2018'
        ),
    cfg(data_name = 'hzg_wminush_M125_2018',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M125_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M125_2018'
        ),
    cfg(data_name = 'hzg_zh_M125_2018',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M125_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M125_2018'
        ),
    cfg(data_name = 'hzg_tth_M125_2018',
        path      = '{0}/ttHToZG_ZToLL_M125_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M125_2018'
        ),
    
    cfg(data_name = 'hzg_gluglu_M130_2018',
        path      = '{0}/GluGluHToZG_M-130_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_gluglu_M130_2018'
        ),
    cfg(data_name = 'hzg_vbf_M130_2018',
        path      = '{0}/VBFHToZG_M-130_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_vbf_M130_2018'
        ),
    cfg(data_name = 'hzg_wplush_M130_2018',
        path      = '{0}/WplusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wplush_M130_2018'
        ),
    cfg(data_name = 'hzg_wminush_M130_2018',
        path      = '{0}/WminusH_HToZG_WToAll_ZToLL_M130_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_wminush_M130_2018'
        ),
    cfg(data_name = 'hzg_zh_M130_2018',
        path      = '{0}/ZH_ZToAll_HToZG_ZToLL_M130_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50,
        suffix    = 'hzg_zh_M130_2018'
        ),
    cfg(data_name = 'hzg_tth_M130_2018',
        path      = '{0}/ttHToZG_ZToLL_M130_13TeV_powheg_RunIIAutumn18MiniAOD/'.format(path),
        nJobs     = 50, 
        suffix    = 'hzg_tth_M130_2018'
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

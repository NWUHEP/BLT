#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys
#MultilepAnalyzer /eos/uscms/store/group/lpcbacon/15/TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_0.root 10000 ttbar_2l2nu ttbar_2l2nu single_lepton 2018 1
#MultilepAnalyzer /eos/uscms/store/user/zchen/EGamma_20191120/2018_EGamma_prod_EGamma_Run2018D-PromptReco-v2/Output_1.root 100000 electron_2018D electron_2018D single_lepton 2018 1
#MultilepAnalyzer /eos/uscms/store/group/lpcbacon/15/SingleMuonRun2018C_17Sep2018_v1/SingleMuonRun2018C_17Sep2018_v1_1.root 100000 muon_2018C muon_2018C single_lepton 2018 1
''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'single_lepton'
period     = '2018'
analyzer   = 'MultilepAnalyzer'

data_samples = ['single_mu', 'single_el']
mc_samples   = ['ttbar', 'wjets', 'zjets', 't', 'diboson']

''' 
    Set job configurations.  
'''
data_dict, mc_dict = {},{}

path       = '/eos/uscms/store/group/lpcbacon/15'
data_dict['single_mu'] = [
        cfg(data_name = 'muon_2018A',
            path      = '{0}/SingleMuonRun2018A_17Sep2018_v2_v2'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2018A'
            ),
        cfg(data_name = 'muon_2018B',
            path      = '{0}/SingleMuonRun2018B_17Sep2018_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2018B'
          ),
        cfg(data_name = 'muon_2018C',
            path      = '{0}/SingleMuonRun2018C_17Sep2018_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2018C'
            ),
        cfg(data_name = 'muon_2018D',
            path      = '{0}/SingleMuonRun2018D_PromptReco_v2_v3'.format(path),
            nJobs     = 100,
            suffix    = 'muon_2018D'
            )
    ]

pathzchen = '/eos/uscms/store/user/zchen/EGamma_20191120'
data_dict['single_el'] = [
        cfg(data_name = 'electron_2018A',
            nJobs    = 50,
            path     = '{0}/2018_EGamma_prod_EGamma_Run2018A-PromptReco-v3'.format(pathzchen),
            suffix   = 'electron_2018A'
           ),

        cfg(data_name = 'electron_2018B',
            path     = '{0}/2018_EGamma_prod_EGamma_Run2018B-PromptReco-v2'.format(pathzchen),
            nJobs    = 50,
            suffix   = 'electron_2018B'
           ),

        cfg(data_name = 'electron_2018C',
            path     = '{0}/2018_EGamma_prod_EGamma_Run2018C-PromptReco-v3'.format(pathzchen),
            nJobs    = 50,
            suffix   = 'electron_2018C'
           ),
           
        cfg(data_name = 'electron_2018D',
            path     = '{0}/2018_EGamma_prod_EGamma_Run2018D-PromptReco-v2'.format(pathzchen),
            nJobs    = 100,
            suffix   = 'electron_2018D'
           ),

        cfg(data_name = 'electron_2018D',
            path     = '{0}/2018_EGamma_prod_EGamma_Run2018E-PromptReco-v1'.format(pathzchen),
            nJobs    = 10,
            suffix   = 'electron_2018D'
        )

    ]



mc_dict['zjets'] = [
    cfg(data_name = 'DY0JetsToLL',
       path     = '/eos/uscms/store/user/zchen/lpcbacon/15/DYJetsToLL_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8',
       nJobs    = 20,
       suffix   = 'DY0JetsToLL'
      ),
    cfg(data_name = 'DY1JetsToLL',
       path     = '/eos/uscms/store/user/zchen/lpcbacon/15/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8',
       nJobs    = 20,
       suffix   = 'DY1JetsToLL'
      ),
    cfg(data_name = 'DY2JetsToLL',
       path     = '/eos/uscms/store/user/zchen/lpcbacon/15/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8',
       nJobs    = 20,
       suffix   = 'DY2JetsToLL'
      ),
] # need reformate folder

mc_dict['wjets'] = [
      cfg(data_name = 'WJetsToLNu_HT_100To200',
        path     = '{0}/WJetsToLNu_HT_100To200_TuneCP5_13TeV'.format(path),
        nJobs    = 50,
        suffix   = 'WJetsToLNu_HT_100To200'
        ), 
        
      cfg(data_name = 'WJetsToLNu_HT_200To400',
        path     = '{0}/WJetsToLNu_HT_200To400_TuneCP5_13TeV'.format(path),
        nJobs    = 50,
        suffix   = 'WJetsToLNu_HT_200To400'
        ), 
      
      cfg(data_name = 'WJetsToLNu_HT_400To600',
        path     = '{0}/WJetsToLNu_HT_400To600_TuneCP5_13TeV'.format(path),
        nJobs    = 50,
        suffix   = 'WJetsToLNu_HT_400To600'
        ), 
        
      cfg(data_name = 'WJetsToLNu_HT_600To800',
        path     = '{0}/WJetsToLNu_HT_600To800_TuneCP5_13TeV'.format(path),
        nJobs    = 50,
        suffix   = 'WJetsToLNu_HT_600To800'
        ), 
      cfg(data_name = 'WJetsToLNu_HT_800To1200',
        path     = '{0}/WJetsToLNu_HT_800To1200_TuneCP5_13TeV'.format(path),
        nJobs    = 50,
        suffix   = 'WJetsToLNu_HT_800To1200'
        ), 
        
      cfg(data_name = 'WJetsToLNu_HT_1200To2500',
        path     = '{0}/WJetsToLNu_HT_1200To2500_TuneCP5_13TeV'.format(path),
        nJobs    = 50,
        suffix   = 'WJetsToLNu_HT_1200To2500'
        ), 

] # need jet binned

mc_dict['diboson'] = [
      cfg(data_name = 'ww',
        path     = '{0}/WW_TuneCP5_13TeV_pythia8'.format(path),
        nJobs    = 20,
        suffix   = 'ww'
       ),    
      cfg(data_name = 'wz',
        path     = '{0}/WZ_TuneCP5_13TeV_pythia8'.format(path),
        nJobs    = 20,
        suffix   = 'wz'
       ), 
      cfg(data_name = 'zz',
        path     = '{0}/ZZ_TuneCP5_13TeV_pythia8'.format(path),
        nJobs    = 20,
        suffix   = 'zz'
       ),    
    ]


mc_dict['ttbar'] = [
    # top
    cfg(data_name = 'ttbar_2l2nu',
        path     = '{0}/TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_2l2nu'
       ),    
    cfg(data_name = 'ttbar_semilepton',
        path     = '{0}/TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_semilepton'
       ),
    cfg(data_name = 'ttbar_hadronic',
        path     = '{0}/TTToHadronic_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_hadronic'
       )
    ]


mc_dict['t'] = [
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 10,
        suffix   = 't_tw'
       ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_tw'
       ),
    ]


batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
batch_list += sum([mc_dict[n] for n in mc_samples], []) 
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       output_dir  = '/store/user/zchen/batchout',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc',
                       analyzer    = analyzer
                     )
batch.submit_to_batch()

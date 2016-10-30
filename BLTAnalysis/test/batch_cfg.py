#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
path       = '/tthome/share/bacon/production/10'
executable = 'execBatch.sh'
selection  = 'mumu'
period     = '2016'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

data_list = []

if selection == 'mumu' or selection == 'emu':
    data_list.extend([
        cfg(data_name = 'muon_2016B',
            path      = '{0}/SingleMuon_Run2016B-PromptReco-v2'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2016B'
           ),
        cfg(data_name = 'muon_2016C',
            path      = '{0}/SingleMuon_Run2016C-PromptReco-v2'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016C'
           ),
        cfg(data_name = 'muon_2016D',
            path      = '{0}/SingleMuon_Run2016D-PromptReco-v2'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016D'
           ),
        #cfg(data_name = 'muon_2016E',
        #    path      = '{0}/SingleMuon_Run2016E-PromptReco-v2'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016E'
        #   ),
        #cfg(data_name = 'muon_2016F',
        #    path      = '{0}/SingleMuon_Run2016F-PromptReco-v1'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016F'
        #   ),
        ])
elif selection == 'ee':
    data_list.extend([
        cfg(data_name = 'electron_2012A',
            path     = '{0}/SingleElectron_2012A-22Jan2013'.format(path),
            nJobs    = 1,
            suffix   = 'electron_2012A'
           ),
        ])

mc_list = []
mc_list.extend([
    cfg(data_name = 'TTJets',
        path     = '{0}/TTJets_13TeV_amcatnloFXFX_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'ttjets'
       ),
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M_50_13TeV_amcatnloFXFX_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/Spring16_DYJetsToLL_M-10to50_TuneCUETP8M1'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_m-10to50'
       ),
    #cfg(data_name = 'ttbar_leptonic',
    #    path     = '{0}/Summer12_TTJets_FullLeptMGDecays'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_lep'
    #   ),
    #cfg(data_name = 'ttbar_semileptonic',
    #    path     = '{0}/Summer12_TTJets_SemiLeptMGDecays'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_semilep'
    #   ),
    #cfg(data_name = 'T_s-channel',
    #    path     = '{0}/'.format(path),
    #    nJobs    = 10,
    #    suffix   = 't_s'
    #   ),
    #cfg(data_name = 'Tbar_s-channel',
    #    path     = '{0}/'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'tbar_s'
    #   ),
    cfg(data_name = 'T_t-channel',
        path     = '{0}/ST_t_channel_top_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1'.format(path),
        nJobs    = 10,
        suffix   = 't_t'
       ),
    cfg(data_name = 'Tbar_t-channel',
        path     = '{0}/ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_t'
       ),
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/ST_tW_top_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M1'.format(path),
        nJobs    = 10,
        suffix   = 't_tw'
       ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/ST_tW_antitop_5f_inclusiveDecays_13TeV_powheg_pythia8_TuneCUETP8M1'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_tw'
       ),
    #cfg(data_name = 'WW',
    #    path     = '{0}/Summer12_WW_TuneZ2star'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'ww'
    #   ),
    #cfg(data_name = 'WZJetsTo2L2Q',
    #    path     = '{0}/Summer12_WZJetsTo2L2Q_TuneZ2star'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'wz_2l2q'
    #   ),
    #cfg(data_name = 'WZJetsTo3LNu',
    #    path     = '{0}/Summer12_WZJetsTo3LNu_TuneZ2'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'wz_3lnu'
    #   ),
    #cfg(data_name = 'ZZJetsTo2L2Nu',
    #    path     = '{0}/Summer12_ZZJetsTo2L2Nu_TuneZ2star'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'zz_2l2nu'
    #   ),
    #cfg(data_name = 'ZZJetsTo2L2Q',
    #    path     = '{0}/Summer12_ZZJetsTo2L2Q_TuneZ2star'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'zz_2l2q'
    #   ),
    ])

'''
sigList = []
sigList.extend([
    cfg(data_name = 'Bprime2Xb_X2mumu',
        path     = '{0}/Summer12_Bprime2Xb_X2mumu'.format(path),
        nJobs    = 5,
        suffix   = 'bprime_xb'
       ),
    ])
'''

batch_list = []
batch_list += data_list
batch_list += mc_list 
#batch_list += sig_list
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'nut3'
                     )
batch.submit_to_batch()


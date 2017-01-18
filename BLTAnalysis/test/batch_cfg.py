#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
path       = '/tthome/share/bacon/production/12'
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
        cfg(data_name = 'muon_2016B_v1',
            path      = '{0}/SingleMuon_Run2016B-23Sep2016-v1'.format(path),
            nJobs     = 20,
            suffix    = 'muon_2016B'
           ),
        cfg(data_name = 'muon_2016C_v1',
            path      = '{0}/SingleMuon_Run2016C-23Sep2016-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016C'
           ),
        cfg(data_name = 'muon_2016D_v1',
            path      = '{0}/SingleMuon_Run2016D-23Sep2016-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016D'
           ),
        cfg(data_name = 'muon_2016E_v1',
            path      = '{0}/SingleMuon_Run2016E-23Sep2016-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016E'
           ),
        cfg(data_name = 'muon_2016F_v1',
            path      = '{0}/SingleMuon_Run2016F-23Sep2016-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016F'
           ),
        cfg(data_name = 'muon_2016G',
            path      = '{0}/SingleMuon_Run2016G-23Sep2016-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016G'
           ),
        cfg(data_name = 'muon_2016H_v1',
            path      = '{0}/SingleMuon_Run2016H-PromptReco-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        cfg(data_name = 'muon_2016H_v2',
            path      = '{0}/SingleMuon_Run2016H-PromptReco-v2'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        cfg(data_name = 'muon_2016H_v3',
            path      = '{0}/SingleMuon_Run2016H-PromptReco-v3'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        ])
elif selection == 'ee':
    data_list.extend([
        cfg(data_name = 'electron_2016A',
            path     = '{0}/SingleElectron_2016A-22Jan2013'.format(path),
            nJobs    = 1,
            suffix   = 'electron_2016A'
           ),
        cfg(data_name = 'electron_2016B',
            path     = '{0}/SingleElectron_2016B-22Jan2013'.format(path),
            nJobs    = 1,
            suffix   = 'electron_2016B'
           ),
        cfg(data_name = 'electron_2016C',
            path     = '{0}/SingleElectron_2016C-22Jan2013'.format(path),
            nJobs    = 1,
            suffix   = 'electron_2016C'
           ),
        ])

mc_list = []
mc_list.extend([
    # Drell-Yan
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Spring16_DYJetsToLL_M-50_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/Spring16_DYJetsToLL_M-10to50_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_m-10to50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-50',
        path     = '{0}/Spring16_DY1JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-10to50',
        path     = '{0}/Spring16_DY1JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-10to50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-50',
        path     = '{0}/Spring16_DY2JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-10to50',
        path     = '{0}/Spring16_DY2JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-10to50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-50',
        path     = '{0}/Spring16_DY3JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-10to50',
        path     = '{0}/Spring16_DY3JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-10to50'
       ),
    cfg(data_name = 'DY4JetsToLL_M-50',
        path     = '{0}/Spring16_DY4JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_m-50'
       ),
    cfg(data_name = 'DY4JetsToLL_M-10to50',
        path     = '{0}/Spring16_DY4JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_m-10to50'
       ),

    # top
    cfg(data_name = 'ttbar_leptonic',
        path     = '{0}/Spring16_TTTo2L2Nu_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_lep'
       ),
    cfg(data_name = 'ttbar_semileptonic',
        path     = '{0}/Spring16_TTToSemiLeptonic_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_semilep'
       ),
    cfg(data_name = 'T_s-channel',
        path     = '{0}/Spring16_ST_s-channel_4f_leptonDecays_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 't_s'
       ),
    #cfg(data_name = 'Tbar_s-channel',
    #    path     = '{0}/'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'tbar_s'
    #   ),
    cfg(data_name = 'T_t-channel',
        path     = '{0}/Spring16_ST_t-channel_4f_leptonDecays_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 't_t'
       ),
    #cfg(data_name = 'Tbar_t-channel',
    #    path     = '{0}/ST_t_channel_antitop_4f_inclusiveDecays_13TeV_powhegV2_madspin_pythia8_TuneCUETP8M1'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'tbar_t'
    #   ),
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/Spring16_ST_tW_top_5f_NoFullyHadronicDecays'.format(path),
        nJobs    = 10,
        suffix   = 't_tw'
       ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/Spring16_ST_tW_antitop_5f_NoFullyHadronicDecays'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_tw'
       ),

    # Diboson
    cfg(data_name = 'WW',
        path     = '{0}/Spring16_WWTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ww'
       ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/Spring16_WZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'wz_2l2q'
       ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/Spring16_WZTo3LNu_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'wz_3lnu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/Spring16_ZZTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2nu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/Spring16_ZZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2q'
       ),
    #cfg(data_name = 'ZZJetsTo4L',
    #    path     = '{0}/Spring16_ZZTo4L_powheg'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'zz_4l'
    #   ),
    ])

sig_list = []
sig_list.extend([
    cfg(data_name = 'Bprime2Xb_X2mumu',
        path     = '{0}/Spring16_Bprime2Xb_X2mumu_tChannel'.format(path),
        nJobs    = 5,
        suffix   = 'bprime_t-channel'
       ),
    cfg(data_name = 'BB2bbdimufatjet',
        path     = '{0}/Spring16_BB2bbdimufatjet'.format(path),
        nJobs    = 5,
        suffix   = 'bprime_bb_semilep'
       ),
    ])

batch_list = []
batch_list += data_list
#batch_list += mc_list 
#batch_list += sig_list
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'nut3'
                     )
batch.submit_to_batch()


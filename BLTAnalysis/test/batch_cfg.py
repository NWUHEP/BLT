#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys

''' Specify parameters '''
cfg        = bm.JobConfig
path       = '/tthome/share/bacon/production/12a'
executable = 'execBatch.sh'
selection  = 'emu'
period     = '2016'

data_samples = []
mc_samples   = ['ttbar']

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

data_list = {}
data_list['single_mu'] = [
        cfg(data_name = 'muon_2016B_v1',
            path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver1-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016B'
           ),
        cfg(data_name = 'muon_2016B_v2',
            path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016B'
           ),
        cfg(data_name = 'muon_2016C_v1',
            path      = '{0}/SingleMuon_Run2016C-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016C'
           ),
        cfg(data_name = 'muon_2016D_v1',
            path      = '{0}/SingleMuon_Run2016D-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016D'
           ),
        cfg(data_name = 'muon_2016E_v1',
            path      = '{0}/SingleMuon_Run2016E-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016E'
           ),
        cfg(data_name = 'muon_2016F_v1',
            path      = '{0}/SingleMuon_Run2016F-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016F'
           ),
        cfg(data_name = 'muon_2016G',
            path      = '{0}/SingleMuon_Run2016G-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016G'
           ),
        cfg(data_name = 'muon_2016H_v2',
            path      = '{0}/SingleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        cfg(data_name = 'muon_2016H_v3',
            path      = '{0}/SingleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        ]

data_list['single_el'] = [
        cfg(data_name = 'electron_2016B',
            nJobs    = 30,
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver1-v1'.format(path),
            suffix   = 'electron_2016B'
           ),
        cfg(data_name = 'electron_2016B',
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver2-v2'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016B'
           ),
        cfg(data_name = 'electron_2016C',
            path     = '{0}/SingleElectron_Run2016C-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016C'
           ),
        cfg(data_name = 'electron_2016D',
            path     = '{0}/SingleElectron_Run2016D-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016D'
           ),
        cfg(data_name = 'electron_2016E',
            path     = '{0}/SingleElectron_Run2016E-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016E'
           ),
        cfg(data_name = 'electron_2016F',
            path     = '{0}/SingleElectron_Run2016F-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016F'
           ),
        cfg(data_name = 'electron_2016G',
            path     = '{0}/SingleElectron_Run2016G-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016G'
           ),
        cfg(data_name = 'electron_2016H',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016H'
           ),
        cfg(data_name = 'electron_2016H',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016H'
           ),
        ]

data_list['mueg'] = [
        cfg(data_name = 'mueg_2016B',
            path     = '{0}/MuonEG_Run2016B-03Feb2017_ver1-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016B'
           ),
        cfg(data_name = 'mueg_2016B',
            path     = '{0}/MuonEG_Run2016B-03Feb2017_ver2-v2'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016B'
           ),
        cfg(data_name = 'mueg_2016C',
            path     = '{0}/MuonEG_Run2016C-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016C'
           ),
        cfg(data_name = 'mueg_2016D',
            path     = '{0}/MuonEG_Run2016D-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016D'
           ),
        cfg(data_name = 'mueg_2016E',
            path     = '{0}/MuonEG_Run2016E-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016E'
           ),
        cfg(data_name = 'mueg_2016F',
            path     = '{0}/MuonEG_Run2016F-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016F'
           ),
        cfg(data_name = 'mueg_2016G',
            path     = '{0}/MuonEG_Run2016G-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016G'
           ),
        cfg(data_name = 'mueg_2016H',
            path     = '{0}/MuonEG_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016H'
           ),
        cfg(data_name = 'mueg_2016H',
            path     = '{0}/MuonEG_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs    = 30,
            suffix   = 'mueg_2016H'
           ),
        ]

path = '/tthome/share/bacon/production/12'
mc_dict = {}
mc_dict['zjets'] = [
    # Drell-Yan
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer16_DYJetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/Summer16_DYJetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_m-10to50'
       ),
    #cfg(data_name = 'DY1JetsToLL_M-50',
    #    path     = '{0}/Summer16_DY1JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z1jets_m-50'
    #   ),
    #cfg(data_name = 'DY1JetsToLL_M-10to50',
    #    path     = '{0}/Summer16_DY1JetsToLL_M-10to50_madgraph_concat'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z1jets_m-10to50'
    #   ),
    #cfg(data_name = 'DY2JetsToLL_M-50',
    #    path     = '{0}/Summer16_DY2JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z2jets_m-50'
    #   ),
    #cfg(data_name = 'DY2JetsToLL_M-10to50',
    #    path     = '{0}/Summer16_DY2JetsToLL_M-10to50_madgraph_concat'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z2jets_m-10to50'
    #   ),
    #cfg(data_name = 'DY3JetsToLL_M-50',
    #    path     = '{0}/Summer16_DY3JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z3jets_m-50'
    #   ),
    #cfg(data_name = 'DY3JetsToLL_M-10to50',
    #    path     = '{0}/Summer16_DY3JetsToLL_M-10to50_madgraph_concat'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z3jets_m-10to50'
    #   ),
    #cfg(data_name = 'DY4JetsToLL_M-50',
    #    path     = '{0}/Summer16_DY4JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z4jets_m-50'
    #   ),
    #cfg(data_name = 'DY4JetsToLL_M-10to50',
    #    path     = '{0}/Summer16_DY4JetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z4jets_m-10to50'
    #   ),
    ]

mc_dict['wjets'] = [
    # W+jets
    cfg(data_name = 'W1JetsToLNu',
        path     = '{0}/Summer16_W1JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w1jets'
       ),
    cfg(data_name = 'W2JetsToLNu',
        path     = '{0}/Summer16_W2JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w2jets'
       ),
    cfg(data_name = 'W3JetsToLNu',
        path     = '{0}/Summer16_W3JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w3jets'
       ),
    cfg(data_name = 'W4JetsToLNu',
        path     = '{0}/Summer16_W4JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w4jets'
       ),
    ]

mc_dict['qcd'] = [
    # QCD
    cfg(data_name = 'QCD_HT50to100',
        path     = '{0}/Summer16_QCD_HT50to100'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht50to100'
       ),
    cfg(data_name = 'QCD_HT100to200',
        path     = '{0}/Summer16_QCD_HT100to200'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht100to200'
       ),
    cfg(data_name = 'QCD_HT200to300',
        path     = '{0}/Summer16_QCD_HT200to300'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht200to300'
       ),
    cfg(data_name = 'QCD_HT300to500',
        path     = '{0}/Summer16_QCD_HT300to500'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht300to500'
       ),
    cfg(data_name = 'QCD_HT500to700',
        path     = '{0}/Summer16_QCD_HT500to700'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht500to700'
       ),
    cfg(data_name = 'QCD_HT700to1000',
        path     = '{0}/Summer16_QCD_HT700to1000'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht700to1000'
       ),
    cfg(data_name = 'QCD_HT1000to1500',
        path     = '{0}/Summer16_QCD_HT1000to1500'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht1000to1500'
       ),
    cfg(data_name = 'QCD_HT1500to2000',
        path     = '{0}/Summer16_QCD_HT1500to2000'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht1500to2000'
       ),
    cfg(data_name = 'QCD_HT2000toInf',
        path     = '{0}/Summer16_QCD_HT2000toInf'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht2000'
       ),
    ]

mc_dict['ttbar'] = [
    # top
    cfg(data_name = 'ttbar_leptonic',
        path     = '{0}/Summer16_TTTo2L2Nu_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_lep'
       ),
    cfg(data_name = 'ttbar_semileptonic',
        path     = '{0}/Summer16_TTToSemilepton_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_semilep'
       ),
    #cfg(data_name = 'ttbar_leptonic',
    #    path     = '{0}/Summer16_TTJets_DiLept_madgraph'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_lep'
    #   ),
    ]

mc_dict['t'] = [
#    cfg(data_name = 'T_s-channel',
#        path     = '{0}/Summer16_ST_s-channel_4f_leptonDecays_amcatnlo'.format(path),
#        nJobs    = 10,
#        suffix   = 't_s'
#       ),
#    cfg(data_name = 'Tbar_s-channel',
#        path     = '{0}/'.format(path),
#        nJobs    = 10,
#        suffix   = 'tbar_s'
#       ),
    #cfg(data_name = 'T_t-channel',
    #    path     = '{0}/Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
    #    nJobs    = 10,
    #    suffix   = 't_t'
    #   ),
    #cfg(data_name = 'Tbar_t-channel',
    #    path     = '{0}/Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'tbar_t'
    #   ),
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/Summer16_ST_tW_top_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 't_tw'
       ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/Summer16_ST_tW_antitop_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_tw'
       ),
    ]

mc_dict['diboson'] = [
    # Diboson
    cfg(data_name = 'WW',
        path     = '{0}/Summer16_WWTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ww'
       ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/Summer16_WZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'wz_2l2q'
       ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/Summer16_WZTo3LNu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'wz_3lnu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer16_ZZTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2nu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer16_ZZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2q'
       ),
    cfg(data_name = 'ZZJetsTo4L',
        path     = '{0}/Summer16_ZZto4L_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_4l'
       ),
    ]

batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
batch_list += sum([mc_dict[n] for n in mc_samples], []) 
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'nut3'
                     )
batch.submit_to_batch()


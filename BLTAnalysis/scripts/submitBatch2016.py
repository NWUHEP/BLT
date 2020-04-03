#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys

# to run interactively
# MultilepAnalyzer /eos/uscms/store/group/lpcbacon/12a/SingleMuon_Run2016C-03Feb2017-v1/SingleMuon_Run2016C-03Feb2017-v1_bacon_00.root 100000 muon_2016C muon_2016C single_lepton 2016 1
# MultilepAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_090.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1
# MultilepAnalyzer /eos/uscms/store/group/lpcbacon/mmackenz/crab_ntupler/WJetsToLNu_amcatnlo/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16/WJetsToLNu_1.root 100000 WJetsToLNu WJetsToLNu single_lepton 2016 1
# MultilepLLTauAnalyzer /eos/uscms/store/group/lpcbacon/12a/SingleMuon_Run2016C-03Feb2017-v1/SingleMuon_Run2016C-03Feb2017-v1_bacon_00.root 100000 muon_2016C muon_2016C single_lepton 2016 1
# MultilepLLTauAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_090.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1




''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'single_lepton'
period     = '2016'
analyzer   = 'MultilepLLTauAnalyzer'

data_samples = ['single_mu', 'single_el']
mc_samples   = ['ttbar', 'wjets', 'zjets', 't', 'diboson']


# data_samples = []
# mc_samples   = ['zjets'] #['ttbar_systematics']

''' 
    Set job configurations.  
'''
data_dict, mc_dict = {},{}

path       = '/eos/uscms/store/group/lpcbacon/12a'
data_dict['single_mu'] = [
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
       cfg(data_name = 'muon_2016C',
           path      = '{0}/SingleMuon_Run2016C-03Feb2017-v1'.format(path),
           nJobs     = 30,
           suffix    = 'muon_2016C'
          ),
       cfg(data_name = 'muon_2016D',
           path      = '{0}/SingleMuon_Run2016D-03Feb2017-v1'.format(path),
           nJobs     = 30,
           suffix    = 'muon_2016D'
          ),
       cfg(data_name = 'muon_2016E',
           path      = '{0}/SingleMuon_Run2016E-03Feb2017-v1'.format(path),
           nJobs     = 30,
           suffix    = 'muon_2016E'
          ),
       cfg(data_name = 'muon_2016F',
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

data_dict['single_el'] = [
        cfg(data_name = 'electron_2016B_v1',
            nJobs    = 30,
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver1-v1'.format(path),
            suffix   = 'electron_2016B'
           ),
        cfg(data_name = 'electron_2016B_v2',
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
        cfg(data_name = 'electron_2016H_v2',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016H'
           ),
        cfg(data_name = 'electron_2016H_v3',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016H'
           ),
        ]



path = '/eos/uscms/store/group/lpcbacon/12d'
mc_dict['zjets'] = [
    # Drell-Yan
    cfg(data_name = 'DYJetsToLL_M-50_amcatnlo',
        path     = '/eos/uscms/store/group/lpcbacon/12/Summer16_DYJetsToLL_M-50_amcatnlo',
        nJobs    = 50,
        suffix   = 'zjets_m-50_amcatnlo'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50_amcatnlo',
        path     = '/eos/uscms/store/group/lpcbacon/12/Summer16_DYJetsToLL_M-10to50_amcatnlo',
        nJobs    = 50,
        suffix   = 'zjets_m-10to50_amcatnlo'
       ),

    cfg(data_name = 'DYToLL_0J_amcatnlo',
        path     = '/eos/uscms/store/user/naodell/bacontuples/Summer16_DYToLL_0J_amcatnlo',
        nJobs    = 10,
        suffix   = 'z0jets_m-50_amcatnlo'
       ),

    cfg(data_name = 'DYToLL_1J_amcatnlo',
        path     = '/eos/uscms/store/user/naodell/bacontuples/Summer16_DYToLL_1J_amcatnlo',
        nJobs    = 10,
        suffix   = 'z1jets_m-50_amcatnlo'
       ),

    cfg(data_name = 'DYToLL_2J_amcatnlo',
        path     = '/eos/uscms/store/user/naodell/bacontuples/Summer16_DYToLL_2J_amcatnlo',
        nJobs    = 10,
        suffix   = 'z2jets_m-50_amcatnlo'
       ),


    # cfg(data_name = 'DYJetsToLL_M-50',
    #    path     = '{0}/DYJetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'zjets_m-50'
    #   ),
    # cfg(data_name = 'DYJetsToLL_M-10to50',
    #    path     = '{0}/DYJetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'zjets_m-10to50'
    #   ),

    # cfg(data_name = 'DY1JetsToLL_M-50',
    #    path     = '{0}/DY1JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z1jets_m-50'
    #   ),
    # cfg(data_name = 'DY1JetsToLL_M-10to50',
    #    path     = '{0}/DY1JetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z1jets_m-10to50'
    #   ),
    # cfg(data_name = 'DY2JetsToLL_M-50',
    #    path     = '{0}/DY2JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z2jets_m-50'
    #   ),
    # cfg(data_name = 'DY2JetsToLL_M-10to50',
    #    path     = '{0}/DY2JetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z2jets_m-10to50'
    #   ),
    # cfg(data_name = 'DY3JetsToLL_M-50',
    #    path     = '{0}/DY3JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z3jets_m-50'
    #   ),
    # cfg(data_name = 'DY3JetsToLL_M-10to50',
    #    path     = '{0}/DY3JetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z3jets_m-10to50'
    #   ),
    # cfg(data_name = 'DY4JetsToLL_M-50',
    #    path     = '{0}/DY4JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z4jets_m-50'
    #   ),
    # cfg(data_name = 'DY4JetsToLL_M-10to50',
    #    path     = '{0}/DY4JetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z4jets_m-10to50'
    #   ),
    ]

path = '/eos/uscms/store/group/lpcbacon/12'
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

mc_dict['ttbar'] = [
    # top
    cfg(data_name = 'ttbar_inclusive',
        path     = '{0}/Summer16_TT_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive'
       ),    
    # leptonic top
    cfg(data_name = 'ttbar_2l2nu',
        path     = '{0}/Summer16_TTTo2L2Nu_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_2l2nu'
       ),
    # semi leptonic top
    cfg(data_name = 'ttbar_semilepton',
        path     = '{0}/Summer16_TTToSemilepton_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_semilepton'
       )
    ]

mc_dict['ttbar_systematics'] = [
    # top full gen
    cfg(data_name = 'ttbar_inclusive_tauReweight',
        path     = '/eos/uscms/store/user/zchen/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/CRAB3/190627_132039/0000',
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_tauReweight'
       ),
    # fsr
    cfg(data_name = 'ttbar_inclusive_fsrUp',
        path     = '{0}/Summer16_TT_powheg_fsrup'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_fsrUp'
       ),
    cfg(data_name = 'ttbar_inclusive_fsrDown',
        path     = '{0}/Summer16_TT_powheg_fsrdown'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_fsrDown'
       ),
    # isr
    cfg(data_name = 'ttbar_inclusive_isrUp',
        path     = '{0}/Summer16_TT_powheg_isrup'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_isrUp'
       ),
    cfg(data_name = 'ttbar_inclusive_isrDown',
        path     = '{0}/Summer16_TT_powheg_isrdown'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_isrDown'
       ),
    # underline event
    cfg(data_name = 'ttbar_inclusive_ueUp',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4up'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_ueUp'
       ),
    cfg(data_name = 'ttbar_inclusive_ueDown',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4down'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_ueDown'
       ),
    # matrix element and parton shower matching
    cfg(data_name = 'ttbar_inclusive_mepsUp',
        path     = '{0}/Summer16_TT_powheg_hdampUP'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_mepsUp'
       ),
    cfg(data_name = 'ttbar_inclusive_mepsDown',
        path     = '{0}/Summer16_TT_powheg_hdampDOWN'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_mepsDown'
       ),
    ]

mc_dict['t'] = [
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

path = '/eos/uscms/store/user/jbueghly/sync_mc'
mc_dict['hzg'] = [
    cfg(data_name = 'hzg_gluglu',
        path      = '{0}/GluGluHToZG_M-125_powheg_calib'.format(path),
        nJobs     = 9,
        suffix    = 'hzg_gluglu'
        ),
    cfg(data_name = 'hzg_vbf',
        path      = '{0}/VBFHToZG_M-125_powheg_calib'.format(path),
        nJobs     = 26,
        suffix    = 'hzg_vbf'
        ),
    cfg(data_name = 'hzg_wplus',
        path      = '{0}/WplusHtoZG_M-125_powheg'.format(path),
        nJobs     = 61,
        suffix    = 'hzg_wh'
        ),
    cfg(data_name = 'hzg_wminus',
        path      = '{0}/WminusHtoZG_M-125_powheg'.format(path),
        nJobs     = 59,
        suffix    = 'hzg_wh'
        ),
    cfg(data_name = 'hzg_zh',
        path      = '{0}/ZHtoZG_M-125_powheg'.format(path),
        nJobs     = 60,
        suffix    = 'hzg_zh'
        ),
    cfg(data_name = 'hzg_tth',
        path      = '{0}/ttHtoZG_M-125_powheg'.format(path),
        nJobs     = 40, 
        suffix    = 'hzg_tth'
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
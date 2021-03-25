#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys


''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'single_lepton'
period     = '2016'




# # to run interactively
# SinglelepAnalyzer /eos/uscms/store/group/lpcbacon/12a/SingleMuon_Run2016C-03Feb2017-v1/SingleMuon_Run2016C-03Feb2017-v1_bacon_00.root 100000 muon_2016C muon_2016C single_lepton 2016 1
# SinglelepAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_W1JetsToLNu/Summer16_W1JetsToLNu_bacon_00.root 100000 w1jet w1jet single_lepton 2016 1
# ####################################
# analyzer     = 'SinglelepAnalyzer'
# data_samples = ['single_mu', 'single_el']
# mc_samples   = ['ttbar', 'wjets', 'zjets', 't', 'diboson','gjets', 'qcd']






# # to run interactively
# # LLTauAnalyzer /eos/uscms/store/group/lpcbacon/12a/SingleMuon_Run2016C-03Feb2017-v1/SingleMuon_Run2016C-03Feb2017-v1_bacon_00.root 100000 muon_2016C muon_2016C single_lepton 2016 1
# # LLTauAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_090.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1
# ####################################
analyzer   = 'LLTauAnalyzer'
data_samples = ['single_mu', 'single_el']
mc_samples   = ['ttbar', 'wjets', 'zjets', 't', 'diboson','gjets']









# to run interactively
# TauAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_090.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1
####################################
#analyzer     = 'TauAnalyzer'
#data_samples = []
#mc_samples   = ['ttbar', 'ttbar_systematics']




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
        nJobs    = 30,
        suffix   = 'zjets_m-50_amcatnlo'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50_amcatnlo',
        path     = '/eos/uscms/store/group/lpcbacon/12/Summer16_DYJetsToLL_M-10to50_amcatnlo',
        nJobs    = 30,
        suffix   = 'zjets_m-10to50_amcatnlo'
       ),
    ]

path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['wjets'] = [
    # W+jets
    cfg(data_name = 'W1JetsToLNu',
        path     = '{0}/Summer16_W1JetsToLNu'.format(path),
        nJobs    = 30,
        suffix   = 'w1jets'
       ),
    cfg(data_name = 'W2JetsToLNu',
        path     = '{0}/Summer16_W2JetsToLNu'.format(path),
        nJobs    = 30,
        suffix   = 'w2jets'
       ),
    cfg(data_name = 'W3JetsToLNu',
        path     = '{0}/Summer16_W3JetsToLNu'.format(path),
        nJobs    = 30,
        suffix   = 'w3jets'
       ),
    cfg(data_name = 'W4JetsToLNu',
        path     = '{0}/Summer16_W4JetsToLNu'.format(path),
        nJobs    = 30,
        suffix   = 'w4jets'
       ),
    ]



mc_dict['gjets'] = [
    # g+jets
    cfg(data_name = 'gjets_ht40to100',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_201210/*/',
        nJobs    = 10,
        suffix   = 'gjets_ht40to100'
       ),
    cfg(data_name = 'gjets_ht100to200',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_200727/*/',
        nJobs    = 10,
        suffix   = 'gjets_ht100to200'
       ),
    cfg(data_name = 'gjets_ht200to400',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_200045/*/',
        nJobs    = 10,
        suffix   = 'gjets_ht200to400'
       ),
    cfg(data_name = 'gjets_ht400to600',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_195447/*/',
        nJobs    = 10,
        suffix   = 'gjets_ht400to600'
       ),
    cfg(data_name = 'gjets_ht600toinf',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_185407/*/',
        nJobs    = 10,
        suffix   = 'gjets_ht600toinf'
       ),
    ]


mc_dict['ttbar'] = [
    # top
    cfg(data_name = 'ttbar_inclusive',
        path     = '{0}/Summer16_TT_powheg'.format(path),
        nJobs    = 30,
        suffix   = 'ttbar_inclusive'
       ), 
    ]

mc_dict['ttbar_systematics'] = [
    # top full gen
    cfg(data_name = 'ttbar_inclusive_tauReweight',
        path     = '/eos/uscms/store/user/zchen/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/CRAB3/200520_221006/0000',
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_tauReweight'
       ),
    # fsr
    cfg(data_name = 'ttbar_inclusive_fsrUp',
        path     = '{0}/Summer16_TT_powheg_fsrup'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_fsrUp'
       ),
    cfg(data_name = 'ttbar_inclusive_fsrDown',
        path     = '{0}/Summer16_TT_powheg_fsrdown'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_fsrDown'
       ),
    # isr
    cfg(data_name = 'ttbar_inclusive_isrUp',
        path     = '{0}/Summer16_TT_powheg_isrup'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_isrUp'
       ),
    cfg(data_name = 'ttbar_inclusive_isrDown',
        path     = '{0}/Summer16_TT_powheg_isrdown'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_isrDown'
       ),
    # underline event
    cfg(data_name = 'ttbar_inclusive_ueUp',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4up'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_ueUp'
       ),
    cfg(data_name = 'ttbar_inclusive_ueDown',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4down'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_ueDown'
       ),
    # matrix element and parton shower matching
    cfg(data_name = 'ttbar_inclusive_mepsUp',
        path     = '{0}/Summer16_TT_powheg_hdampUP'.format(path),
        nJobs    = 20,
        suffix   = 'ttbar_inclusive_mepsUp'
       ),
    cfg(data_name = 'ttbar_inclusive_mepsDown',
        path     = '{0}/Summer16_TT_powheg_hdampDOWN'.format(path),
        nJobs    = 20,
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

    cfg(data_name = 'T_t-channel',
        path     = '{0}/Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 't_t'
       ),
    cfg(data_name = 'Tbar_t-channel',
        path     = '{0}/Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_t'
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

mc_dict['qcd'] = [
    # qcd
    cfg(data_name = 'qcd_ht50to100',
        path     = '{0}/Summer16_QCD_HT50to100'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht50to100'
       ),
    cfg(data_name = 'qcd_ht100to200',
        path     = '{0}/Summer16_QCD_HT100to200'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht100to200'
       ),
    cfg(data_name = 'qcd_ht200to300',
        path     = '{0}/Summer16_QCD_HT200to300'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht200to300'
       ),
    
    cfg(data_name = 'qcd_ht300to500',
        path     = '{0}/Summer16_QCD_HT300to500'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht300to500'
       ),
       
    cfg(data_name = 'qcd_ht500to700',
        path     = '{0}/Summer16_QCD_HT500to700'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht500to700'
       ),
    cfg(data_name = 'qcd_ht700to1000',
        path     = '{0}/Summer16_QCD_HT700to1000'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht700to1000'
       ),
    cfg(data_name = 'qcd_ht1000to1500',
        path     = '{0}/Summer16_QCD_HT1000to1500'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht1000to1500'
       ),
    
    cfg(data_name = 'qcd_ht1500to2000',
        path     = '{0}/Summer16_QCD_HT1500to2000'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht1500to2000'
       ),

    cfg(data_name = 'qcd_ht2000toInf',
        path     = '{0}/Summer16_QCD_HT2000toInf'.format(path),
        nJobs    = 20,
        suffix   = 'qcd_ht2000toInf'
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
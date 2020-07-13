#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys

''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'single_lepton'
period     = '2016'

data_samples = []#'single_mu', 'single_el']
mc_samples   = ['gjets']#'zjets', 'zjets_ext', 'ttbar', 'diboson', 't', 'wjets']#, 'wjets_alt', 'qcd']

data_dict = {}
mc_dict = {}

''' 
    Set job configurations.  
'''

path = '/eos/uscms/store/group/lpcbacon/12a'
data_dict['single_mu'] = [
#        cfg(dataname = 'muon_2016B_v1',
#            path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver1-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016B'
#           ),
#        cfg(dataname = 'muon_2016B_v2',
#            path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016B'
#           ),
#        cfg(dataname = 'muon_2016C',
#            path      = '{0}/SingleMuon_Run2016C-03Feb2017-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016C'
#           ),
#        cfg(dataname = 'muon_2016D',
#            path      = '{0}/SingleMuon_Run2016D-03Feb2017-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016D'
#           ),
#        cfg(dataname = 'muon_2016E',
#            path      = '{0}/SingleMuon_Run2016E-03Feb2017-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016E'
#           ),
#        cfg(dataname = 'muon_2016F',
#            path      = '{0}/SingleMuon_Run2016F-03Feb2017-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016F'
#           ),
        cfg(dataname = 'muon_2016G',
            path      = '{0}/SingleMuon_Run2016G-03Feb2017-v1'.format(path),
            n_jobs     = 30,
            suffix    = 'muon_2016G'
           ),
#        cfg(dataname = 'muon_2016H_v2',
#            path      = '{0}/SingleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016H'
#           ),
#        cfg(dataname = 'muon_2016H_v3',
#            path      = '{0}/SingleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
#            n_jobs     = 30,
#            suffix    = 'muon_2016H'
#           ),
        ]

data_dict['single_el'] = [
        cfg(dataname = 'electron_2016B_v1',
            n_jobs    = 30,
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver1-v1'.format(path),
            suffix   = 'electron_2016B'
           ),
        cfg(dataname = 'electron_2016B_v2',
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver2-v2'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016B'
           ),
        cfg(dataname = 'electron_2016C',
            path     = '{0}/SingleElectron_Run2016C-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016C'
           ),
        cfg(dataname = 'electron_2016D',
            path     = '{0}/SingleElectron_Run2016D-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016D'
           ),
        cfg(dataname = 'electron_2016E',
            path     = '{0}/SingleElectron_Run2016E-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016E'
           ),
        cfg(dataname = 'electron_2016F',
            path     = '{0}/SingleElectron_Run2016F-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016F'
           ),
        cfg(dataname = 'electron_2016G',
            path     = '{0}/SingleElectron_Run2016G-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016G'
           ),
        cfg(dataname = 'electron_2016H_v2',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver2-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016H'
           ),
        cfg(dataname = 'electron_2016H_v3',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver3-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'electron_2016H'
           ),
        ]

data_dict['mueg'] = [
        cfg(dataname = 'mueg_2016B_v1',
            path     = '{0}/MuonEG_Run2016B-03Feb2017_ver1-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016B'
           ),
        cfg(dataname = 'mueg_2016B_v2',
            path     = '{0}/MuonEG_Run2016B-03Feb2017_ver2-v2'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016B'
           ),
        cfg(dataname = 'mueg_2016C',
            path     = '{0}/MuonEG_Run2016C-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016C'
           ),
        cfg(dataname = 'mueg_2016D',
            path     = '{0}/MuonEG_Run2016D-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016D'
           ),
        cfg(dataname = 'mueg_2016E',
            path     = '{0}/MuonEG_Run2016E-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016E'
           ),
        cfg(dataname = 'mueg_2016F',
            path     = '{0}/MuonEG_Run2016F-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016F'
           ),
        cfg(dataname = 'mueg_2016G',
            path     = '{0}/MuonEG_Run2016G-03Feb2017-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016G'
           ),
        cfg(dataname = 'mueg_2016H_v2',
            path     = '{0}/MuonEG_Run2016H-03Feb2017_ver2-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016H'
           ),
        cfg(dataname = 'mueg_2016H_v3',
            path     = '{0}/MuonEG_Run2016H-03Feb2017_ver3-v1'.format(path),
            n_jobs    = 30,
            suffix   = 'mueg_2016H'
           ),
        ]

# Drell-Yan
## LO madgraph sample with parton binning
#path = '/eos/uscms/store/group/lpcbacon/12d'
#mc_dict['zjets'] = [
#    cfg(dataname = 'DYJetsToLL_M-50',
#        path     = '{0}/DYJetsToLL_M-50_madgraph'.format(path),
#        n_jobs    = 50,
#        suffix   = 'zjets_m-50'
#       ),
#    cfg(dataname = 'DYJetsToLL_M-10to50',
#        path     = '{0}/DYJetsToLL_M-10to50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'zjets_m-10to50'
#       ),
#    cfg(dataname = 'DY1JetsToLL_M-50',
#        path     = '{0}/DY1JetsToLL_M-50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z1jets_m-50'
#       ),
#    cfg(dataname = 'DY1JetsToLL_M-10to50',
#        path     = '{0}/DY1JetsToLL_M-10to50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z1jets_m-10to50'
#       ),
#    cfg(dataname = 'DY2JetsToLL_M-50',
#        path     = '{0}/DY2JetsToLL_M-50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z2jets_m-50'
#       ),
#    cfg(dataname = 'DY2JetsToLL_M-10to50',
#        path     = '{0}/DY2JetsToLL_M-10to50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z2jets_m-10to50'
#       ),
#    cfg(dataname = 'DY3JetsToLL_M-50',
#        path     = '{0}/DY3JetsToLL_M-50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z3jets_m-50'
#       ),
#    cfg(dataname = 'DY3JetsToLL_M-10to50',
#        path     = '{0}/DY3JetsToLL_M-10to50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z3jets_m-10to50'
#       ),
#    cfg(dataname = 'DY4JetsToLL_M-50',
#        path     = '{0}/DY4JetsToLL_M-50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z4jets_m-50'
#       ),
#    cfg(dataname = 'DY4JetsToLL_M-10to50',
#        path     = '{0}/DY4JetsToLL_M-10to50_madgraph'.format(path),
#        n_jobs    = 10,
#        suffix   = 'z4jets_m-10to50'
#       ),
#    ]

# amc@nlo drell-yan sample
path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['zjets'] = [
    cfg(dataname = 'DYJetsToLL_M-50_alt',
        path     = '{0}/Summer16_DYJetsToLL_M-50_amcatnlo'.format(path),
        n_jobs    = 50,
        suffix   = 'zjets_m-50_alt'
       ),
    cfg(dataname = 'DYJetsToLL_M-10to50_alt',
        path     = '{0}/Summer16_DYJetsToLL_M-10to50_amcatnlo'.format(path),
        n_jobs    = 20,
        suffix   = 'zjets_m-10to50_alt'
       ),
    ]

path = '/eos/uscms/store/user/naodell/bacontuples'
mc_dict['zjets_ext'] = [
    cfg(dataname = 'DY0JetsToLL_alt',
        path     = '{0}/Summer16_DYToLL_0J_amcatnlo'.format(path),
        n_jobs    = 50,
        suffix   = 'z0jets_alt'
       ),
    cfg(dataname = 'DY1JetsToLL_alt',
        path     = '{0}/Summer16_DYToLL_1J_amcatnlo'.format(path),
        n_jobs    = 50,
        suffix   = 'z1jets_alt'
       ),
    cfg(dataname = 'DY2JetsToLL_alt',
        path     = '{0}/Summer16_DYToLL_2J_amcatnlo'.format(path),
        n_jobs    = 50,
        suffix   = 'z2jets_alt'
       ),
    ]

path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['wjets'] = [
    # W+jets (MADGRAPH jet binned)
    cfg(dataname = 'W1JetsToLNu',
        path     = '{0}/Summer16_W1JetsToLNu'.format(path),
        n_jobs    = 40,
        suffix   = 'w1jets'
       ),
    cfg(dataname = 'W2JetsToLNu',
        path     = '{0}/Summer16_W2JetsToLNu'.format(path),
        n_jobs    = 40,
        suffix   = 'w2jets'
       ),
    cfg(dataname = 'W3JetsToLNu',
        path     = '{0}/Summer16_W3JetsToLNu'.format(path),
        n_jobs    = 40,
        suffix   = 'w3jets'
       ),
    cfg(dataname = 'W4JetsToLNu',
        path     = '{0}/Summer16_W4JetsToLNu'.format(path),
        n_jobs    = 40,
        suffix   = 'w4jets'
       ),
    ]

path = '/eos/uscms/store/group/lpcbacon/mmackenz/crab_ntupler/WJetsToLNu_amcatnlo'
mc_dict['wjets_alt'] = [
    # W+jets (AMC@NLO)
    cfg(dataname = 'WJetsToLNu',
        path     = '{0}/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16'.format(path),
        n_jobs    = 40,
        suffix   = 'wjets_alt'
       ),
    cfg(dataname = 'WJetsToLNu_ext1',
        path     = '{0}/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16_ext1'.format(path),
        n_jobs    = 40,
        suffix   = 'wjets_alt_ext1'
       ),
    cfg(dataname = 'WJetsToLNu_ext2',
        path     = '{0}/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16_ext2'.format(path),
        n_jobs    = 40,
        suffix   = 'wjets_alt_ext2'
       ),
    ]

path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['qcd'] = [
    # QCD
    cfg(dataname = 'QCD_HT50to100',
        path     = '{0}/Summer16_QCD_HT50to100'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht50to100'
       ),
    cfg(dataname = 'QCD_HT100to200',
        path     = '{0}/Summer16_QCD_HT100to200'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht100to200'
       ),
    cfg(dataname = 'QCD_HT200to300',
        path     = '{0}/Summer16_QCD_HT200to300'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht200to300'
       ),
    cfg(dataname = 'QCD_HT300to500',
        path     = '{0}/Summer16_QCD_HT300to500'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht300to500'
       ),
    cfg(dataname = 'QCD_HT500to700',
        path     = '{0}/Summer16_QCD_HT500to700'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht500to700'
       ),
    cfg(dataname = 'QCD_HT700to1000',
        path     = '{0}/Summer16_QCD_HT700to1000'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht700to1000'
       ),
    cfg(dataname = 'QCD_HT1000to1500',
        path     = '{0}/Summer16_QCD_HT1000to1500'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht1000to1500'
       ),
    cfg(dataname = 'QCD_HT1500to2000',
        path     = '{0}/Summer16_QCD_HT1500to2000'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht1500to2000'
       ),
    cfg(dataname = 'QCD_HT2000toInf',
        path     = '{0}/Summer16_QCD_HT2000toInf'.format(path),
        n_jobs    = 10,
        suffix   = 'qcd_ht2000'
       ),
    ]

mc_dict['ttbar'] = [
    # top
    cfg(dataname = 'ttbar_inclusive',
        path     = '{0}/Summer16_TT_powheg'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive'
       ),
    cfg(dataname = 'ttbar_leptonic',
        path     = '{0}/Summer16_TTTo2L2Nu_powheg'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_lep'
       ),
    cfg(dataname = 'ttbar_semileptonic',
        path     = '{0}/Summer16_TTToSemilepton_powheg'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_semilep'
       ),
    #cfg(dataname = 'ttbar_inclusive_herwig',
    #    path     = '{0}/Summer16_TT_powheg_herwig'.format(path),
    #    n_jobs    = 50,
    #    suffix   = 'ttbar_inclusive_herwig'
    #    ),
    #cfg(dataname = 'ttbar_leptonic',
    #    path     = '{0}/Summer16_TTJets_DiLept_madgraph'.format(path),
    #    n_jobs    = 50,
    #    suffix   = 'ttbar_lep'
    #   ),
    ]

mc_dict['ttbar_syst'] = [
    cfg(dataname = 'ttbar_inclusive_tuneup',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4up'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_tuneup'
       ),
    cfg(dataname = 'ttbar_inclusive_tunedown',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4down'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_tunedown'
       ),
    cfg(dataname = 'ttbar_inclusive_isrup',
        path     = '{0}/Summer16_TT_powheg_isrup'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_isrup'
       ),
    cfg(dataname = 'ttbar_inclusive_isrdown',
        path     = '{0}/Summer16_TT_powheg_isrdown'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_isrdown'
       ),
    cfg(dataname = 'ttbar_inclusive_fsrup',
        path     = '{0}/Summer16_TT_powheg_fsrup'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_fsrup'
       ),
    cfg(dataname = 'ttbar_inclusive_fsrdown',
        path     = '{0}/Summer16_TT_powheg_fsrdown'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_fsrdown'
       ),
    cfg(dataname = 'ttbar_inclusive_hdampup',
        path     = '{0}/Summer16_TT_powheg_hdampUP'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_hdampup'
        ),
    cfg(dataname = 'ttbar_inclusive_hdampdown',
        path     = '{0}/Summer16_TT_powheg_hdampDOWN'.format(path),
        n_jobs    = 50,
        suffix   = 'ttbar_inclusive_hdampdown'
       ),
    ]

path = '/eos/uscms/store/group/lpcbacon/12/ttbar_systematics_ext'
mc_dict['ttbar_syst_ext'] = [
    cfg(dataname = 'ttbar_inclusive_tuneup_ext1',
        path     = '{0}/Summer16_TT_powheg_tuneup'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_tuneup_ext1'
       ),
    cfg(dataname = 'ttbar_inclusive_tunedown_ext1',
        path     = '{0}/Summer16_TT_powheg_tunedown'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_tunedown_ext1'
       ),
    #cfg(dataname = 'ttbar_inclusive_isrup_ext1',
    #    path     = '{0}/Summer16_TT_powheg_isrup'.format(path),
    #    n_jobs    = 20,
    #    suffix   = 'ttbar_inclusive_isrup_ext1'
    #   ),
    cfg(dataname = 'ttbar_inclusive_isrup_ext2',
        path     = '{0}/Summer16_TT_powheg_isrup/ext2'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_isrup_ext2'
       ),
    #cfg(dataname = 'ttbar_inclusive_isrdown_ext1',
    #    path     = '{0}/Summer16_TT_powheg_isrdown'.format(path),
    #    n_jobs    = 20,
    #    suffix   = 'ttbar_inclusive_isrdown_ext1'
    #   ),
    cfg(dataname = 'ttbar_inclusive_isrdown_ext2',
        path     = '{0}/Summer16_TT_powheg_isrdown/ext2'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_isrdown_ext2'
       ),
    cfg(dataname = 'ttbar_inclusive_fsrup_ext1',
        path     = '{0}/Summer16_TT_powheg_fsrup/ext1'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_fsrup_ext1'
       ),
    cfg(dataname = 'ttbar_inclusive_fsrup_ext2',
        path     = '{0}/Summer16_TT_powheg_fsrup/ext2'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_fsrup_ext2'
       ),
    cfg(dataname = 'ttbar_inclusive_fsrdown_ext1',
        path     = '{0}/Summer16_TT_powheg_fsrdown/ext1'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_fsrdown_ext1'
       ),
    cfg(dataname = 'ttbar_inclusive_fsrdown_ext2',
        path     = '{0}/Summer16_TT_powheg_fsrdown/ext2'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_fsrdown_ext2'
       ),
    cfg(dataname = 'ttbar_inclusive_hdampup_ext1',
        path     = '{0}/Summer16_TT_powheg_hdampup'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_hdampup_ext1'
       ),
    cfg(dataname = 'ttbar_inclusive_hdampdown_ext1',
        path     = '{0}/Summer16_TT_powheg_hdampdown'.format(path),
        n_jobs    = 20,
        suffix   = 'ttbar_inclusive_hdampdown_ext1'
       ),
    ]

path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['t'] = [
#    cfg(dataname = 'T_s-channel',
#        path     = '{0}/Summer16_ST_s-channel_4f_leptonDecays_amcatnlo'.format(path),
#        n_jobs    = 10,
#        suffix   = 't_s'
#       ),
#    cfg(dataname = 'Tbar_s-channel',
#        path     = '{0}/'.format(path),
#        n_jobs    = 10,
#        suffix   = 'tbar_s'
#       ),
    #cfg(dataname = 'T_t-channel',
    #    path     = '{0}/Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
    #    n_jobs    = 10,
    #    suffix   = 't_t'
    #   ),
    #cfg(dataname = 'Tbar_t-channel',
    #    path     = '{0}/Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
    #    n_jobs    = 10,
    #    suffix   = 'tbar_t'
    #   ),
    cfg(dataname = 'T_tW-channel',
        path     = '{0}/Summer16_ST_tW_top_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        n_jobs    = 10,
        suffix   = 't_tw'
       ),
    cfg(dataname = 'Tbar_tW-channel',
        path     = '{0}/Summer16_ST_tW_antitop_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        n_jobs    = 10,
        suffix   = 'tbar_tw'
       ),
    ]

mc_dict['diboson'] = [
    # Diboson
    cfg(dataname = 'qqWW',
        path     = '{0}/Summer16_WWTo2L2Nu_powheg'.format(path),
        n_jobs    = 10,
        suffix   = 'ww_qq'
       ),
    cfg(dataname = 'ggWW',
        path     = '/eos/uscms/store/user/tgunter/makingBacon/13TeV/COMPLETE/GluGluWWTo2L2Nu_summer',
        n_jobs    = 10,
        suffix   = 'ww_gg'
       ),
    cfg(dataname = 'WZJetsTo2L2Q',
        path     = '{0}/Summer16_WZTo2L2Q_amcatnlo'.format(path),
        n_jobs    = 10,
        suffix   = 'wz_2l2q'
       ),
    cfg(dataname = 'WZJetsTo3LNu',
        path     = '{0}/Summer16_WZTo3LNu_powheg'.format(path),
        n_jobs    = 10,
        suffix   = 'wz_3lnu'
       ),
    cfg(dataname = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer16_ZZTo2L2Nu_powheg'.format(path),
        n_jobs    = 10,
        suffix   = 'zz_2l2nu'
       ),
    cfg(dataname = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer16_ZZTo2L2Q_amcatnlo'.format(path),
        n_jobs    = 10,
        suffix   = 'zz_2l2q'
       ),
    cfg(dataname = 'ZZJetsTo4L',
        path     = '{0}/Summer16_ZZto4L_amcatnlo'.format(path),
        n_jobs    = 10,
        suffix   = 'zz_4l'
       ),
    ]

mc_dict['gjets'] = [
    # g+jets
    cfg(dataname = 'gjets_ht40to100',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_201210/*',
        n_jobs    = 10,
        suffix   = 'gjets_ht40to100'
       ),
    cfg(dataname = 'gjets_ht100to200',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_200727/*',
        n_jobs    = 10,
        suffix   = 'gjets_ht100to200'
       ),
    cfg(dataname = 'gjets_ht200to400',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_200045/*',
        n_jobs    = 10,
        suffix   = 'gjets_ht200to400'
       ),
    cfg(dataname = 'gjets_ht400to600',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_195447/*',
        n_jobs    = 10,
        suffix   = 'gjets_ht400to600'
       ),
    cfg(dataname = 'gjets_ht600toinf',
        path     = '/eos/uscms/store/user/zchen/GJets_DR-0p4_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/CRAB3/200709_185407/*',
        n_jobs    = 10,
        suffix   = 'gjets_ht600toinf'
       ),
    ]

path = '/eos/uscms/store/user/mmackenz/private_mc/2016/'
mc_dict['misc'] = [
    cfg(dataname = 'ZToMuTau',
        path     = '{0}/ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        n_jobs    = 10,
        suffix   = 'zjets_mutau'
       ),
    ]


batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
batch_list += sum([mc_dict[n] for n in mc_samples], []) 
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       output_dir  = '/store/user/naodell/batch_output',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc'
                     )
batch.submit_to_batch()

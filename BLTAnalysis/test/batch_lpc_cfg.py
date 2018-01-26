#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
#selection  = '4mu'
#selection  = '2e2mu'
selection  = 'jpsi_control'
#selection = 'mumu'
#selection = 'ee'
period     = '2016'
#path       = '/eos/uscms/store/user/jbueghly/jbueghly_data_multicrab/DoubleMuon/'
#path       = '/eos/uscms/store/user/jbueghly/jbueghly_data_multicrab/DoubleEG/'
path       = '/eos/uscms/store/user/jbueghly/jbueghly_data_multicrab/Charmonium/'
#path       = '/eos/uscms/store/user/naodell/bacontuples/12_vertex/'
#path       = '/eos/uscms/store/user/jbueghly/ref_double_mu/'
executable = 'execBatch.sh'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

data_list = []
data_list.extend([

    # Double muon data
    #cfg(data_name = 'muon_2016B_v2',
    #    path     = '{0}/DoubleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016B'
    #   ),
    #cfg(data_name = 'muon_2016C_v1',
    #    path     = '{0}/DoubleMuon_Run2016C-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016C'
    #   ),
    #cfg(data_name = 'muon_2016D_v1',
    #    path     = '{0}/DoubleMuon_Run2016D-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016D'
    #   ),
    #cfg(data_name = 'muon_2016E_v1',
    #    path     = '{0}/DoubleMuon_Run2016E-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016E'
    #   ),
    #cfg(data_name = 'muon_2016F_v1',
    #    path     = '{0}/DoubleMuon_Run2016F-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016F'
    #   ),
    #cfg(data_name = 'muon_2016G_v1',
    #    path     = '{0}/DoubleMuon_Run2016G-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016G'
    #   ),
    #cfg(data_name = 'muon_2016H_v2',
    #    path     = '{0}/DoubleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016H'
    #   ),
    #cfg(data_name = 'muon_2016H_v3',
    #    path     = '{0}/DoubleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'muon_2016H'
    #   ),

    # Double electron data
    #cfg(data_name = 'electron_2016B_v2',
    #    path     = '{0}/DoubleEG_Run2016B-03Feb2017_ver2-v2'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016B'
    #   ),
    #cfg(data_name = 'electron_2016C_v1',
    #    path     = '{0}/DoubleEG_Run2016C-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016C'
    #   ),
    #cfg(data_name = 'electron_2016D_v1',
    #    path     = '{0}/DoubleEG_Run2016D-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016D'
    #   ),
    #cfg(data_name = 'electron_2016E_v1',
    #    path     = '{0}/DoubleEG_Run2016E-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016E'
    #   ),
    #cfg(data_name = 'electron_2016F_v1',
    #    path     = '{0}/DoubleEG_Run2016F-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016F'
    #   ),
    #cfg(data_name = 'electron_2016G_v1',
    #    path     = '{0}/DoubleEG_Run2016G-03Feb2017-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016G'
    #   ),
    #cfg(data_name = 'electron_2016H_v2',
    #    path     = '{0}/DoubleEG_Run2016H-03Feb2017_ver2-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016H'
    #   ),
    #cfg(data_name = 'electron_2016H_v3',
    #    path     = '{0}/DoubleEG_Run2016H-03Feb2017_ver3-v1'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'electron_2016H'
    #   ),

    # Charmonium data
    cfg(data_name = 'muon_2016B_v2',
        path     = '{0}/Charmonium_Run2016B-03Feb2017_ver2-v2'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016B'
       ),
    cfg(data_name = 'muon_2016C_v1',
        path     = '{0}/Charmonium_Run2016C-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016C'
       ),
    cfg(data_name = 'muon_2016D_v1',
        path     = '{0}/Charmonium_Run2016D-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016D'
       ),
    cfg(data_name = 'muon_2016E_v1',
        path     = '{0}/Charmonium_Run2016E-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016E'
       ),
    cfg(data_name = 'muon_2016F_v1',
        path     = '{0}/Charmonium_Run2016F-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016F'
       ),
    cfg(data_name = 'muon_2016G_v1',
        path     = '{0}/Charmonium_Run2016G-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016G'
       ),
    cfg(data_name = 'muon_2016H_v2',
        path     = '{0}/Charmonium_Run2016H-03Feb2017_ver2-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016H'
       ),
    cfg(data_name = 'muon_2016H_v3',
        path     = '{0}/Charmonium_Run2016H-03Feb2017_ver3-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016H'
       ),
    ])

mc_list = []
path       = '/eos/uscms/store/user/jbueghly/jbueghly_mc_multicrab'
mc_list.extend([
    #cfg(data_name = 'JpsiToMuMu',
    #    path     = '{0}/JpsiToMuMu_JpsiPt8'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'JpsiToMuMu'
    #   ),
    #cfg(data_name = 'BtoJpsi',
    #    path     = '{0}/BToJPsiKMu_JpsiMuMu'.format(path),
    #    nJobs    = 20,
    #    suffix   = 'BtoJpsi'
    #   ),
    # Drell-Yan
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/DYJetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_m-10to50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-50',
        path     = '{0}/DY1JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-10to50',
        path     = '{0}/DY1JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-10to50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-50',
        path     = '{0}/DY2JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-10to50',
        path     = '{0}/DY2JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-10to50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-50',
        path     = '{0}/DY3JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-10to50',
        path     = '{0}/DY3JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-10to50'
       ),
    #cfg(data_name = 'DY4JetsToLL_M-50',
    #    path     = '{0}/DY4JetsToLL_M-50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z4jets_m-50'
    #   ),
    #cfg(data_name = 'DY4JetsToLL_M-10to50',
    #    path     = '{0}/DY4JetsToLL_M-10to50_madgraph'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'z4jets_m-10to50'
    #   ),

    # top
    #cfg(data_name = 'ttbar_leptonic',
    #    path     = '{0}/TTTo2L2Nu_powheg'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_lep'
    #   ),
#    cfg(data_name = 'ttbar_semileptonic',
#        path     = '{0}/Summer16_TTToSemilepton_powheg'.format(path),
#        nJobs    = 50,
#        suffix   = 'ttbar_semilep'
#       ),
#     cfg(data_name = 'ttbar_leptonic',
#         path     = '{0}/TTJets_DiLept_madgraph/all'.format(path),
#         nJobs    = 50,
#         suffix   = 'ttbar_lep'
#        ),
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
    #cfg(data_name = 'T_tW-channel',
    #    path     = '{0}/Summer16_ST_tW_top_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
    #    nJobs    = 10,
    #    suffix   = 't_tw'
    #   ),
    #cfg(data_name = 'Tbar_tW-channel',
    #    path     = '{0}/Summer16_ST_tW_antitop_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'tbar_tw'
    #   ),

    # Diboson
    #cfg(data_name = 'WW',
    #    path     = '{0}/WWTo2L2Nu_powheg'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'ww'
    #   ),
    #cfg(data_name = 'WZJetsTo2L2Q',
    #    path     = '{0}/WZTo2L2Q_amcatnlo'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'wz_2l2q'
    #   ),
    #cfg(data_name = 'WZJetsTo3LNu',
    #    path     = '{0}/WZTo3LNu_powheg'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'wz_3lnu'
    #   ),
    #cfg(data_name = 'ZZJetsTo2L2Nu',
    #    path     = '{0}/ZZTo2L2Nu_powheg'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'zz_2l2nu'
    #   ),
    #cfg(data_name = 'ZZJetsTo2L2Q',
    #    path     = '{0}/ZZTo2L2Q_amcatnlo'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'zz_2l2q'
    #   ),
    cfg(data_name = 'ZZJetsTo4L',
        path     = '{0}/ZZTo4L_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_4l'
       ),
#    cfg(data_name = 'GluGluHToZZTo4L',
#        path     = '{0}/GluGluHToZZTo4L_powheg/all'.format(path),
#        nJobs    = 10,
#        suffix   = 'gluglu_h_zz_4l'
#       ),
#    cfg(data_name = 'VBFHToZZTo4L',
#        path     = '{0}/VBF_HToZZTo4L_powheg/all'.format(path),
#        nJobs    = 10,
#        suffix   = 'vbf_h_zz_4l'
#       ),
    ])

#path       = '/eos/uscms/store/user/jbueghly/12_vertex'
#sig_list = []
#sig_list.extend([
#    cfg(data_name = 'zjpsi_singlet',
#        path      = '{0}/zjpsi_singlet'.format(path),
#        nJobs     = 5,
#        suffix    = 'zjpsi_singlet'
#        ),
#    ])

batch_list = []
batch_list += data_list
#batch_list += mc_list
#batch_list += sig_list

#batch = bm.BatchMaster(configList = data_list, 
#                      shortQueue = False,
#                      stageDir   = 'batch',
#                      executable = executable,
#                      selection  = selection,
#                      location   = 'lpc'
#                     )
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc'
                     )
batch.submit_to_batch()


#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys

# to run interactively
# MultileptonAnalyzer /eos/uscms/store/group/lpcbacon/12a/SingleMuon_Run2016C-03Feb2017-v1/SingleMuon_Run2016C-03Feb2017-v1_bacon_00.root 100000 muon_2016C muon_2016C single_lepton 2016 1
# MultileptonAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_000.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1
# ExampleAnalyzer /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_000.root 1000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1

''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'single_lepton'
period     = '2016'

data_samples = ['single_mu']
mc_samples   = ['ttbar', 't']


''' 
    Set job configurations.  
'''
data_dict, mc_dict = {},{}



temp = []
path = '/eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017B_17Nov2017_v1/SingleMuon/CRAB3/181208_161717'
for i in range(9):
    print('{}/{:04d}'.format(path,i))
    temp.append( 
        cfg(data_name = 'muon_2016B_v1_{:04d}'.format(i),
            path      = '{}/{:04d}'.format(path,i),
            nJobs     = 10,
            suffix    = 'muon_2016B'))

path = '/eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017C_17Nov2017_v1/SingleMuon/CRAB3/181208_161820'
for i in range(9):
    print('{}/{:04d}'.format(path,i))
    temp.append( 
        cfg(data_name = 'muon_2016B_v1_{:04d}'.format(i),
            path      = '{}/{:04d}'.format(path,i),
            nJobs     = 10,
            suffix    = 'muon_2016C'))

path = '/eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017D_17Nov2017_v1/SingleMuon/CRAB3/181208_161909'
for i in range(4):
    print('{}/{:04d}'.format(path,i))
    temp.append( 
        cfg(data_name = 'muon_2016D_v1_{:04d}'.format(i),
            path      = '{}/{:04d}'.format(path,i),
            nJobs     = 10,
            suffix    = 'muon_2016D'))

path = '/eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017E_17Nov2017_v1/SingleMuon/CRAB3/181208_162023'
for i in range(9):
    print('{}/{:04d}'.format(path,i))
    temp.append( 
        cfg(data_name = 'muon_2016E_v1_{:04d}'.format(i),
            path      = '{}/{:04d}'.format(path,i),
            nJobs     = 10,
            suffix    = 'muon_2016E'))

path = '/eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017F_17Nov2017_v1/SingleMuon/CRAB3/181209_112437'
for i in range(9):
    print('{}/{:04d}'.format(path,i))
    temp.append(
        cfg(data_name = 'muon_2016F_v1_{:04d}'.format(i),
            path      = '{}/{:04d}'.format(path,i),
            nJobs     = 10,
            suffix    = 'muon_2016F'))

data_dict['single_mu'] = temp




path = '/eos/uscms/store/group/lpcbacon/15'
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
    cfg(data_name = 'ttbar_hadron',
        path     = '{0}/TTToHadronic_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_hadron'
       )
    ]


mc_dict['t'] = [
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 20,
        suffix   = 't_tw'
       ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8'.format(path),
        nJobs    = 20,
        suffix   = 'tbar_tw'
       ),
    ]

####################
# testing
####################

# path = '/eos/uscms/store/user/zchen/lpcbacon/15'
# data_dict['single_mu'] = [
#     cfg(data_name = 'muon_2017B_v1',
#         path     = '{0}/SingleMuonRun2017B_17Nov2017_v1'.format(path),
#         nJobs    = 50,
#         suffix   = 'muon_2017B'
#        ),
#     ]

# path = '/eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017B_17Nov2017_v1/SingleMuon/CRAB3/181208_161717'
# data_dict['single_mu'] = [
#     cfg(data_name = 'muon_2017B_v1',
#         path     = '{0}/0000'.format(path),
#         nJobs    = 50,
#         suffix   = 'muon_2017B'
#        ),
#     ]


batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
batch_list += sum([mc_dict[n] for n in mc_samples], []) 
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc'
                     )
batch.submit_to_batch()

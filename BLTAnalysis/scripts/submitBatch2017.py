#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys
#MultilepAnalyzer /eos/uscms/store/group/lpcbacon/15/SingleMuonRun2017B_17Nov2017_v1/SingleMuonRun2017B_17Nov2017_v1_0.root 10000 muon_2017B muon_2017B single_lepton 2017 1

''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'single_lepton'
period     = '2017'
analyzer   = 'MultilepAnalyzer'

data_samples = ['single_mu', 'single_el']
mc_samples   = ['ttbar', 'wjets', 'zjets', 't', 'diboson']

''' 
    Set job configurations.  
'''
data_dict, mc_dict = {},{}

path       = '/eos/uscms/store/group/lpcbacon/15'
data_dict['single_mu'] = [
        cfg(data_name = 'muon_2017B',
            path      = '{0}/SingleMuonRun2017B_17Nov2017_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2017B'
            ),
        cfg(data_name = 'muon_2017C',
            path      = '{0}/SingleMuonRun2017C_17Nov2017_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2017C'
          ),
        cfg(data_name = 'muon_2017D',
            path      = '{0}/SingleMuonRun2017D_17Nov2017_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2017D'
            ),
        cfg(data_name = 'muon_2017E',
            path      = '{0}/SingleMuonRun2017E_17Nov2017_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2017E'
            ),
        cfg(data_name = 'muon_2017F',
            path      = '{0}/SingleMuonRun2017F_17Nov2017_v1'.format(path),
            nJobs     = 50,
            suffix    = 'muon_2017F'
            ),
    ]

data_dict['single_el'] = [
        cfg(data_name = 'electron_2017B',
            nJobs    = 50,
            path     = '{0}/SingleElectronRun2017B_31Mar2018_v1'.format(path),
            suffix   = 'electron_2017B'
           ),

        cfg(data_name = 'electron_2017C',
            path     = '{0}/SingleElectronRun2017C_31Mar2018_v1'.format(path),
            nJobs    = 50,
            suffix   = 'electron_2017C'
           ),

        cfg(data_name = 'electron_2017D',
            path     = '{0}/SingleElectronRun2017D_31Mar2018_v1'.format(path),
            nJobs    = 50,
            suffix   = 'electron_2017D'
           ),
           
        cfg(data_name = 'electron_2017E',
            path     = '{0}/SingleElectronRun2017E_31Mar2018_v1'.format(path),
            nJobs    = 50,
            suffix   = 'electron_2017E'
           ),

        cfg(data_name = 'electron_2017F',
            path     = '{0}/SingleElectronRun2017F_31Mar2018_v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2017F'
        )

    ]



mc_dict['zjets'] = [] # need reformate folder

mc_dict['wjets'] = [] # need jet binned

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

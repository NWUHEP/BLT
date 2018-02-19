#! /usr/bin/env python


import BLT.BLTAnalysis.BatchMaster as bm
import sys

cfg = bm.JobConfig

''' Specify parameters '''
path       = '/tthome/share/bacon/production/12a'
executable = 'fakerate_execBatch.sh'
selection  = '3mu'
period     = '2016'

''' 
    Set job configurations.  The order of arguments is:  (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

data_list = []
if selection in ['3mu']:
    data_list.extend([
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
        ])



path       = '/tthome/share/bacon/production/12'
mc_list = []

mc_list.extend([
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

    cfg(data_name = 'DY1JetsToLL_M-50',
        path     = '{0}/Summer16_DY1JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'z1jets_m-50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-10to50',
        path     = '{0}/Summer16_DY1JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-10to50'
      ),
    cfg(data_name = 'DY2JetsToLL_M-50',
        path     = '{0}/Summer16_DY2JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'z2jets_m-50'
       ),
    cfg(data_name = 'DYJ2etsToLL_M-10to50',
        path     = '{0}/Summer16_DY2JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-10to50'
      ),
    cfg(data_name = 'DY3JetsToLL_M-50',
        path     = '{0}/Summer16_DY3JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'z3jets_m-50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-10to50',
        path     = '{0}/Summer16_DY3JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-10to50'
      ),
        cfg(data_name = 'DY4JetsToLL_M-50',
        path     = '{0}/Summer16_DY4JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'z4jets_m-50'
       ),
    cfg(data_name = 'DYJ4etsToLL_M-10to50',
        path     = '{0}/Summer16_DY4JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_m-10to50'
      ),
    # wjets
    cfg(data_name = 'W1JetsToLNu',
        path     = '{0}/Summer16_W1JetsToLNu'.format(path),
        nJobs    = 10,
        suffix   = 'w1jets'
        ),
    cfg(data_name = 'W2JetsToLNu',
        path     = '{0}/Summer16_W2JetsToLNu'.format(path),
        nJobs    = 10,
        suffix   = 'w2jets'
        ),
    cfg(data_name = 'W3JetsToLNu',
        path     = '{0}/Summer16_W3JetsToLNu'.format(path),
        nJobs    = 10,
        suffix   = 'w3jets'
        ),
    cfg(data_name = 'W4JetsToLNu',
        path     = '{0}/Summer16_W4JetsToLNu'.format(path),
        nJobs    = 10,
        suffix   = 'w4jets'
        )
])



batch_list = []
batch_list += data_list
batch_list += mc_list 

batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'nut3'
                     )
batch.submit_to_batch()
# end


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


    #cfg(data_name = 'ttbar_leptonic',
    #    path     = '{0}/Summer16_TTTo2L2Nu_powheg'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_lep'
    #   ),
    
    #cfg(data_name = 'ttbar_semileptonic',
    #    path     = '{0}/Summer16_TTToSemilepton_powheg'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_semilep'
    #   ),
    #cfg(data_name = 'ttbar_semileptonic_minus',
    #    path     = '{0}/Summer16_TTJets_SingleLeptFromTbar_amcatnlo'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_semilep'
    #   ),
    #cfg(data_name = 'ttbar_leptonic',
    #    path     = '{0}/Summer16_TTJets_DiLept_madgraph'.format(path),
    #    nJobs    = 50,
    #    suffix   = 'ttbar_lep'
    #   ),


        #cfg(data_name = 'T_s-channel',
    #    path     = '{0}/Summer16_ST_s-channel_4f_leptonDecays_amcatnlo'.format(path),
    #    nJobs    = 10,
    #    suffix   = 't_s'
    #   ),
    #cfg(data_name = 'Tbar_s-channel',
    #    path     = '{0}/'.format(path),
    #    nJobs    = 10,
    #    suffix   = 'tbar_s'
    #   ),

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
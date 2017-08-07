#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = '4mu'
period     = '2016'
path       = '/eos/uscms/store/user/jbueghly/zjpsi_data_multicrab_v2_concat/DoubleMuon'
executable = 'execBatch.sh'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

data_list = []
data_list.extend([
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
    cfg(data_name = 'muon_2016E_v1',
        path     = '{0}/DoubleMuon_Run2016E-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016E'
       ),
    cfg(data_name = 'muon_2016F_v1',
        path     = '{0}/DoubleMuon_Run2016F-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016F'
       ),
    cfg(data_name = 'muon_2016G_v1',
        path     = '{0}/DoubleMuon_Run2016G-03Feb2017-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016G'
       ),
    cfg(data_name = 'muon_2016H_v2',
        path     = '{0}/DoubleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016H'
       ),
    cfg(data_name = 'muon_2016H_v3',
        path     = '{0}/DoubleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016H'
       ),
    ])

batch_list = []
batch_list += data_list

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


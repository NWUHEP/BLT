#! /usr/bin/env python
import BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = 'amumu'
period     = '2016'
path       = '/eos/uscms/store/user/naodell/bacontuples'
executable = 'execBatch.sh'
location   = 'lpc'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

dataList = []
dataList.extend([
    cfg(dataName = 'muon_2016B',
        path     = '{0}/SingleMuonRun2016B_PromptReco_v2'.format(path),
        nJobs    = 50,
        args     = 'muon_2016B muon 2016',
       ),
    ])

batch = bm.BatchMaster(configList = dataList, 
                      shortQueue = False,
                      stageDir   = 'batch',
                      executable = executable,
                      selection  = '{0}_{1}'.format(selection, period),
                      location   = 'lpc'
                     )
batch.submit_to_batch()


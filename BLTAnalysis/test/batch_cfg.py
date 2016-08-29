#! /usr/bin/env python
import BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
path       = '/tthome/share/noobs/bacon/production/10'
executable = 'execBatch.sh'
selection  = 'amumu'
period     = '2016'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

dataList = []
dataList.extend([
    cfg(dataName = 'muon_2016B',
        path     = '{0}/SingleMuonRun2016B_PromptReco_v2'.format(path),
        nJobs    = 50,
        args     = 'muon_2016B muon 2016'
       ),
    cfg(dataName = 'muon_2016C',
        path     = '{0}/SingleMuonRun2016C_PromptReco_v2'.format(path),
        nJobs    = 30,
        args     = 'muon_2016C muon 2016'
       ),
    cfg(dataName = 'muon_2016D',
        path     = '{0}/SingleMuonRun2016D_PromptReco_v2'.format(path),
        nJobs    = 30,
        args     = 'muon_2016D muon 2016'
       ),
    #cfg(dataName = 'muon_2016E',
    #    path     = '{0}/SingleMuonRun2016E_PromptReco_v2'.format(path),
    #    nJobs    = 30,
    #    args     = 'muon_2016E muon 2016'
    #   ),
    #cfg(dataName = 'muon_2016F',
    #    path     = '{0}/SingleMuonRun2016F_PromptReco_v1'.format(path),
    #    nJobs    = 25,
    #    args     = 'muon_2016F muon 2016'
    #   ),
    ])

batch = bm.BatchMaster(configList = dataList, 
                      shortQueue = False,
                      stageDir   = 'batch',
                      executable = executable,
                      selection  = '{0}_{1}'.format(selection, period),
                      location   = 'nut3'
                     )
batch.submit_to_batch()


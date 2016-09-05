#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
path       = '/tthome/share/bacon/production/04'
#path       = '/eos/uscms/store/user/naodell/bacontuples'
executable = 'execBatch.sh'
selection  = 'amumu'
period     = '2012'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

dataList = []
dataList.extend([
    cfg(dataName = 'muon_2012A',
        path     = '{0}/SingleMu_2012A-22Jan2013'.format(path),
        nJobs    = 1,
        args     = 'muon_2012A muon 2012'
       ),
    cfg(dataName = 'muon_2012B',
        path     = '{0}/SingleMu_2012B-22Jan2013'.format(path),
        nJobs    = 15,
        args     = 'muon_2012B muon 2012'
       ),
    cfg(dataName = 'muon_2012C',
        path     = '{0}/SingleMu_2012C-22Jan2013'.format(path),
        nJobs    = 15,
        args     = 'muon_2012C muon 2012'
       ),
    cfg(dataName = 'muon_2012D',
        path     = '{0}/SingleMu_2012D-22Jan2013'.format(path),
        nJobs    = 15,
        args     = 'muon_2012D muon 2012'
       )
    ])

mcList = []
mcList.extend([
    cfg(dataName = 'ttbar_leptonic',
        path     = '{0}/Summer12_TTJets_FullLeptMGDecays'.format(path),
        nJobs    = 50,
        args     = 'ttbar_lep_2012 muon 2012'
       )
    ])

batchList = mcList
batch = bm.BatchMaster(configList = batchList, 
                       shortQueue = False,
                       stageDir   = 'batch',
                       executable = executable,
                       selection  = '{0}_{1}'.format(selection, period),
                       location   = 'nut3'
                     )
batch.submit_to_batch()


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
        args     = 'ttbar_lep muon 2012'
       ),
    cfg(dataName = 'DYJetsToLL_M-50',
        path     = '{0}/Summer12_DYJetsToLL_M-50_TuneZ2Star'.format(path),
        nJobs    = 50,
        args     = 'zjets_m-50 muon 2012'
       ),
    cfg(dataName = 'DYJetsToLL_M-10to50',
        path     = '{0}/Summer12_DYJetsToLL_M-10to50filter'.format(path),
        nJobs    = 10,
        args     = 'zjets_m-10to50 muon 2012'
       ),
    cfg(dataName = 'T_s-channel',
        path     = '{0}/Summer12_T_s-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 't_s muon 2012'
       ),
    cfg(dataName = 'Tbar_s-channel',
        path     = '{0}/Summer12_Tbar_s-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'tbar_s muon 2012'
       ),
    cfg(dataName = 'T_t-channel',
        path     = '{0}/Summer12_T_t-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 't_t muon 2012'
       ),
    cfg(dataName = 'Tbar_t-channel',
        path     = '{0}/Summer12_Tbar_t-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'tbar_t muon 2012'
       ),
    cfg(dataName = 'T_tW-channel',
        path     = '{0}/Summer12_T_tW-channel-DR_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 't_tw muon 2012'
       ),
    cfg(dataName = 'Tbar_tW-channel',
        path     = '{0}/Summer12_Tbar_tW-channel-DR_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'tbar_tw muon 2012'
       ),
    cfg(dataName = 'WW',
        path     = '{0}/Summer12_WW_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'ww muon 2012'
       ),
    cfg(dataName = 'WZJetsTo2L2Q',
        path     = '{0}/Summer12_WZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'wz_2l2q muon 2012'
       ),
    cfg(dataName = 'WZJetsTo3LNu',
        path     = '{0}/Summer12_WZJetsTo3LNu_TuneZ2'.format(path),
        nJobs    = 10,
        args     = 'wz_3lnu muon 2012'
       ),
    cfg(dataName = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer12_ZZJetsTo2L2Nu_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'zz_2l2nu muon 2012'
       ),
    cfg(dataName = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer12_ZZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 10,
        args     = 'zz_2l2q muon 2012'
       ),
    ])

batchList = mcList + dataList
batch = bm.BatchMaster(configList = batchList, 
                       shortQueue = False,
                       stageDir   = 'batch',
                       executable = executable,
                       selection  = '{0}_{1}'.format(selection, period),
                       location   = 'nut3'
                     )
batch.submit_to_batch()


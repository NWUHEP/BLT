#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
path       = '/tthome/share/bacon/production/04'
#path       = '/eos/uscms/store/user/naodell/bacontuples'
executable = 'execBatch.sh'
selection  = 'ee'
period     = '2012'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

dataList = []
if selection == 'mumu':
    dataList.extend([
        cfg(data_name = 'muon_2012A',
            path     = '{0}/SingleMu_2012A-22Jan2013'.format(path),
            nJobs    = 1,
            suffix   = 'muon_2012A'
           ),
        cfg(data_name = 'muon_2012B',
            path     = '{0}/SingleMu_2012B-22Jan2013'.format(path),
            nJobs    = 15,
            suffix   = 'muon_2012B'
           ),
        cfg(data_name = 'muon_2012C',
            path     = '{0}/SingleMu_2012C-22Jan2013'.format(path),
            nJobs    = 15,
            suffix   = 'muon_2012C'
           ),
        cfg(data_name = 'muon_2012D',
            path     = '{0}/SingleMu_2012D-22Jan2013'.format(path),
            nJobs    = 15,
            suffix   = 'muon_2012D'
           )
        ])
elif selection == 'ee':
    dataList.extend([
        cfg(data_name = 'electron_2012A',
            path     = '{0}/SingleElectron_2012A-22Jan2013'.format(path),
            nJobs    = 1,
            suffix   = 'electron_2012A'
           ),
        cfg(data_name = 'electron_2012B',
            path     = '{0}/SingleElectron_2012B-22Jan2013'.format(path),
            nJobs    = 15,
            suffix   = 'electron_2012B'
           ),
        cfg(data_name = 'electron_2012C',
            path     = '{0}/SingleElectron_2012C-22Jan2013'.format(path),
            nJobs    = 15,
            suffix   = 'electron_2012C'
           ),
        cfg(data_name = 'electron_2012D',
            path     = '{0}/SingleElectron_2012D-22Jan2013'.format(path),
            nJobs    = 15,
            suffix   = 'electron_2012D'
           )
        ])

mcList = []
mcList.extend([
    cfg(data_name = 'ttbar_leptonic',
        path     = '{0}/Summer12_TTJets_FullLeptMGDecays'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_lep'
       ),
    cfg(data_name = 'ttbar_semileptonic',
        path     = '{0}/Summer12_TTJets_SemiLeptMGDecays'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_semilep'
       ),
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer12_DYJetsToLL_M-50_TuneZ2Star'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/Summer12_DYJetsToLL_M-10to50filter'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_m-10to50'
       ),
    cfg(data_name = 'T_s-channel',
        path     = '{0}/Summer12_T_s-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 't_s'
       ),
    cfg(data_name = 'Tbar_s-channel',
        path     = '{0}/Summer12_Tbar_s-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_s'
       ),
    cfg(data_name = 'T_t-channel',
        path     = '{0}/Summer12_T_t-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 't_t'
       ),
    cfg(data_name = 'Tbar_t-channel',
        path     = '{0}/Summer12_Tbar_t-channel_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_t'
       ),
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/Summer12_T_tW-channel-DR_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 't_tw'
       ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/Summer12_Tbar_tW-channel-DR_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_tw'
       ),
    cfg(data_name = 'WW',
        path     = '{0}/Summer12_WW_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'ww'
       ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/Summer12_WZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'wz_2l2q'
       ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/Summer12_WZJetsTo3LNu_TuneZ2'.format(path),
        nJobs    = 10,
        suffix   = 'wz_3lnu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer12_ZZJetsTo2L2Nu_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2nu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer12_ZZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2q'
       ),
    ])

sigList = []
sigList.extend([
    cfg(data_name = 'Bprime2Xb_X2mumu',
        path     = '{0}/Summer12_Bprime2Xb_X2mumu'.format(path),
        nJobs    = 5,
        suffix   = 'bprime_xb'
       ),
    ])

batchList = []
batchList += mcList 
batchList += dataList
batchList += sigList
batch = bm.BatchMaster(config_list = batchList, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'nut3'
                     )
batch.submit_to_batch()


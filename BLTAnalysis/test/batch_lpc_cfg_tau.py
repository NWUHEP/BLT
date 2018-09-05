#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys
# FakeSelector /eos/uscms/store/group/lpcbacon/12d/DYJetsToLL_M-50_madgraph/Output_0.root 10000 DYJetsToLL_M-50 DYJetsToLL_M-50 single_lepton 2016 1

# TauSelector /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_000.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1
# FakeSelector /eos/uscms/store/group/lpcbacon/12/Summer16_TT_powheg/Summer16_TT_powheg_bacon_000.root 100000 ttbar_inclusive ttbar_inclusive single_lepton 2016 1

''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch_lpc_tau.sh'
selection  = 'single_lepton'
period     = '2016'


mc_samples   = ['ttbar','ttbar_theory']

''' 
    Set job configurations.  
'''


mc_dict = {}
path = '/eos/uscms/store/group/lpcbacon/12'


mc_dict['ttbar'] = [
    # top
    cfg(data_name = 'ttbar_inclusive',
        path     = '{0}/Summer16_TT_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive'
       )

    ]

mc_dict['ttbar_theory'] = [
    cfg(data_name = 'ttbar_inclusive_tunedown',
       path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4down'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_down'
      ),
    cfg(data_name = 'ttbar_inclusive_tuneup',
       path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4up'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_up'
      ),
    cfg(data_name = 'ttbar_inclusive_isrdown',
       path     = '{0}/Summer16_TT_powheg_isrdown'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_isrdown'
      ),
    cfg(data_name = 'ttbar_inclusive_isrup',
       path     = '{0}/Summer16_TT_powheg_isrup'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_isrup'
      ),
    cfg(data_name = 'ttbar_inclusive_fsrdown',
       path     = '{0}/Summer16_TT_powheg_fsrdown'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_fsrdown'
      ),
    cfg(data_name = 'ttbar_inclusive_fsrup',
       path     = '{0}/Summer16_TT_powheg_fsrup'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_fsrup'
      ),
    cfg(data_name = 'ttbar_inclusive_hdampdown',
       path     = '{0}/Summer16_TT_powheg_hdampDOWN'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_hdampdown'
      ),
    cfg(data_name = 'ttbar_inclusive_hdampup',
       path     = '{0}/Summer16_TT_powheg_hdampUP'.format(path),
       nJobs    = 50,
       suffix   = 'ttbar_inclusive_hdampup'
       )
    ]



batch_list = []
batch_list += sum([mc_dict[n] for n in mc_samples], []) 
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc'
                     )
batch.submit_to_batch()


#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm

import sys

''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execGenMCBatch.sh'
selection  = 'single_lepton'
period     = '2016'

''' Whether or not to actually submit '''
dryrun = False

gen_samples  = [ ]
digi_samples = [ ]
aod_samples  = [ ]
mini_samples = [ ]
bacon_samples = [ 'zmutau' ]
''' 
    Set job configurations.  
'''
gen_dict, digi_dict, aod_dict, mini_dict, bacon_dict = {},{},{},{}, {}

path = '/eos/uscms/store/user/mmackenz/private_mc/2016'
##################################################################
# Generation dictionary setup
##################################################################
gen_dict['zemu'] = [
    cfg(data_name = 'z_emu',
        path      = '',
        nJobs     = 100,
        suffix    = 'z_emu',
        arguments = 'ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016.py'
    ),
]
gen_dict['zetau'] = [
    cfg(data_name = 'z_etau',
        path      = '',
        nJobs     = 100,
        suffix    = 'z_etau',
        arguments = 'ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016.py'        
    ),
]
gen_dict['zmutau'] = [
    cfg(data_name = 'z_mutau',
        path      = '',
        nJobs     = 1,
        suffix    = 'z_mutau',
        arguments = 'ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016.py'        
    ),
]
gen_dict['hemu'] = [
    cfg(data_name = 'h_emu',
        path      = '',
        nJobs     = 100,
        suffix    = 'h_emu',
        arguments = 'HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016.py'        
    ),
]
gen_dict['hetau'] = [
    cfg(data_name = 'h_etau',
        path      = '',
        nJobs     = 100,
        suffix    = 'h_etau',
        arguments = 'HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016.py'        
    ),
]
gen_dict['hmutau'] = [
    cfg(data_name = 'h_mutau',
        path      = '',
        nJobs     = 100,
        suffix    = 'h_mutau',
        arguments = 'HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016.py'        
    ),
]

##################################################################
# Digitization dictionary setup
##################################################################
digi_dict['zemu'] = [
    cfg(data_name = 'z_emu',
        path      = '{0}/ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'z_emu',
        arguments = 'ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

digi_dict['zetau'] = [
    cfg(data_name = 'z_etau',
        path      = '{0}/ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'z_etau',
        arguments = 'ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

digi_dict['zmutau'] = [
    cfg(data_name = 'z_mutau',
        path      = '{0}/ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'z_mutau',
        arguments = 'ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

digi_dict['hemu'] = [
    cfg(data_name = 'h_emu',
        path      = '{0}/HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'h_emu',
        arguments = 'HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

digi_dict['hetau'] = [
    cfg(data_name = 'h_etau',
        path      = '{0}/HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'h_etau',
        arguments = 'HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

digi_dict['hmutau'] = [
    cfg(data_name = 'h_mutau',
        path      = '{0}/HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'h_mutau',
        arguments = 'HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

##################################################################
# AOD dictionary setup
##################################################################
aod_dict['zemu'] = [
    cfg(data_name = 'z_emu',
        path      = '{0}/ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-Digi'.format(path),
        nJobs     = 100,
        suffix    = 'z_emu',
        arguments = 'ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016_AOD.py'
    ),
]

aod_dict['zetau'] = [
    cfg(data_name = 'z_etau',
        path      = '{0}/ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-Digi'.format(path),
        nJobs     = 100,
        suffix    = 'z_etau',
        arguments = 'ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016_AOD.py'
    ),
]

aod_dict['zmutau'] = [
    cfg(data_name = 'z_mutau',
        path      = '{0}/ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-Digi'.format(path),
        nJobs     = 100,
        suffix    = 'z_mutau',
        arguments = 'ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_AOD.py'
    ),
]

aod_dict['hemu'] = [
    cfg(data_name = 'h_emu',
        path      = '{0}/HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-Digi'.format(path),
        nJobs     = 100,
        suffix    = 'h_emu',
        arguments = 'HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016_AOD.py'
    ),
]

aod_dict['hetau'] = [
    cfg(data_name = 'h_etau',
        path      = '{0}/HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-Digi'.format(path),
        nJobs     = 100,
        suffix    = 'h_etau',
        arguments = 'HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016_AOD.py'
    ),
]

aod_dict['hmutau'] = [
    cfg(data_name = 'h_mutau',
        path      = '{0}/HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-Digi'.format(path),
        nJobs     = 100,
        suffix    = 'h_mutau',
        arguments = 'HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_AOD.py'
    ),
]

aod_dict['hmutau'] = [
    cfg(data_name = 'h_mutau',
        path      = '{0}/HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-Gen'.format(path),
        nJobs     = 100,
        suffix    = 'h_mutau',
        arguments = 'HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_Digi.py'
    ),
]

##################################################################
# MINIAOD dictionary setup
##################################################################
mini_dict['zemu'] = [
    cfg(data_name = 'z_emu',
        path      = '{0}/ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-AOD'.format(path),
        nJobs     = 100,
        suffix    = 'z_emu',
        arguments = 'ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016_MINIAOD.py'
    ),
]

mini_dict['zetau'] = [
    cfg(data_name = 'z_etau',
        path      = '{0}/ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-AOD'.format(path),
        nJobs     = 100,
        suffix    = 'z_etau',
        arguments = 'ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016_MINIAOD.py'
    ),
]

mini_dict['zmutau'] = [
    cfg(data_name = 'z_mutau',
        path      = '{0}/ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-AOD'.format(path),
        nJobs     = 100,
        suffix    = 'z_mutau',
        arguments = 'ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_MINIAOD.py'
    ),
]

mini_dict['hemu'] = [
    cfg(data_name = 'h_emu',
        path      = '{0}/HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-AOD'.format(path),
        nJobs     = 100,
        suffix    = 'h_emu',
        arguments = 'HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016_MINIAOD.py'
    ),
]

mini_dict['hetau'] = [
    cfg(data_name = 'h_etau',
        path      = '{0}/HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-AOD'.format(path),
        nJobs     = 100,
        suffix    = 'h_etau',
        arguments = 'HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016_MINIAOD.py'
    ),
]

mini_dict['hmutau'] = [
    cfg(data_name = 'h_mutau',
        path      = '{0}/HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-AOD'.format(path),
        nJobs     = 100,
        suffix    = 'h_mutau',
        arguments = 'HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016_MINIAOD.py'
    ),
]

##################################################################
# Bacon dictionary setup
##################################################################
bacon_dict['zemu'] = [
    cfg(data_name = 'z_emu',
        path      = '{0}/ZToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        nJobs     = 100,
        suffix    = 'z_emu',
        arguments = 'private_mc_bacon.py'
    ),
]

bacon_dict['zetau'] = [
    cfg(data_name = 'z_etau',
        path      = '{0}/ZToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        nJobs     = 100,
        suffix    = 'z_etau',
        arguments = 'private_mc_bacon.py'
    ),
]

bacon_dict['zmutau'] = [
    cfg(data_name = 'z_mutau',
        path      = '{0}/ZToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        nJobs     = 100,
        suffix    = 'z_mutau',
        arguments = 'private_mc_bacon.py'
    ),
]

bacon_dict['hemu'] = [
    cfg(data_name = 'h_emu',
        path      = '{0}/HToEMu_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        nJobs     = 100,
        suffix    = 'h_emu',
        arguments = 'private_mc_bacon.py'
    ),
]

bacon_dict['hetau'] = [
    cfg(data_name = 'h_etau',
        path      = '{0}/HToETau_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        nJobs     = 100,
        suffix    = 'h_etau',
        arguments = 'private_mc_bacon.py'
    ),
]

bacon_dict['hmutau'] = [
    cfg(data_name = 'h_mutau',
        path      = '{0}/HToMuTau_13TeV_TuneCUETP8M1_Filter_RunII2016-MINIAOD'.format(path),
        nJobs     = 100,
        suffix    = 'h_mutau',
        arguments = 'private_mc_bacon.py'
    ),
]

gen_list   = []
gen_list  += sum([gen_dict[n]  for n in gen_samples ], []) 
digi_list  = []
digi_list += sum([digi_dict[n] for n in digi_samples], []) 
aod_list   = []
aod_list  += sum([aod_dict[n]  for n in aod_samples ], []) 
mini_list  = []
mini_list += sum([mini_dict[n] for n in mini_samples], []) 
bacon_list  = []
bacon_list += sum([bacon_dict[n] for n in bacon_samples], []) 

gen_batch = bm.BatchMaster(config_list = gen_list, 
                           stage_dir   = 'batch',
                           output_dir  = '/store/user/mmackenz/batch_output',
                           selection   = selection,
                           period      = period,
                           executable  = executable,
                           location    = 'lpc',
                           dryrun      = dryrun,
)
digi_batch = bm.BatchMaster(config_list = digi_list, 
                            stage_dir   = 'batch',
                            output_dir  = '/store/user/mmackenz/batch_output',
                            selection   = selection,
                            period      = period,
                            executable  = executable,
                            location    = 'lpc',
                            dryrun      = dryrun,
)
aod_batch = bm.BatchMaster(config_list = aod_list, 
                           stage_dir   = 'batch',
                           output_dir  = '/store/user/mmackenz/batch_output',
                           selection   = selection,
                           period      = period,
                           executable  = executable,
                           location    = 'lpc',
                           dryrun      = dryrun,
)
mini_batch = bm.BatchMaster(config_list = mini_list, 
                            stage_dir   = 'batch',
                            output_dir  = '/store/user/mmackenz/batch_output',
                            selection   = selection,
                            period      = period,
                            executable  = executable,
                            location    = 'lpc',
                            dryrun      = dryrun,
)
bacon_batch = bm.BatchMaster(config_list = bacon_list, 
                             stage_dir   = 'batch',
                             output_dir  = '/store/user/mmackenz/batch_output',
                             selection   = selection,
                             period      = period,
                             executable  = executable,
                             location    = 'lpc',
                             dryrun      = dryrun,
)

if len(gen_list) > 0:
    gen_batch.submit_to_batch()
if len(digi_list) > 0:
    digi_batch.submit_to_batch()
if len(aod_list) > 0:
    aod_batch.submit_to_batch()
if len(mini_list) > 0:
    mini_batch.submit_to_batch()
if len(bacon_list) > 0:
    bacon_batch.submit_to_batch()

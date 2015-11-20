#!/usr/bin/env python

import os
import sys
import subprocess
import logging
import string

#from CRABClient.JobType.UserTarball import UserTarball
from pack import UserTarball
from CRABClient.UserUtilities import config


# Based on CRABClient.JobType.BasicJobType and CRABClient.JobType.Analysis
class BasicJobType(object):
    '''
    BasicJobType
    '''

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.projdir = 'blt_projects'
        self.srcdir = 'sourceFiles'
        self.tarFilename = 'default.tgz'
        self.txtFilename = 'input.txt'
        self.jdlFilename = 'blt.jdl'
        self.exeFilename = 'blt.sh'
        self.jobid = 1
        self.njobs = self.config.User.njobs

        self.configurations = {
            'PROJDIR'      : self.projdir,
            'SRCDIR'       : self.srcdir,
            'TARBALL'      : self.tarFilename,
            'SOURCEFILE'   : self.txtFilename,
            'JOBAD'        : self.jdlFilename,
            'EXECUTABLE'   : self.exeFilename,
            'JOBID'        : self.jobid,
            'MACHINE'      : 'cmslpc',
            'MAXEVENTS'    : '-1',
            'JOBNAME'      : 'job',
            'LOGNAME'      : 'res',
        }

        for k, v in self.config.User.dictionary_().iteritems():
            self.configurations[string.upper(k)] = v
        self.configurations['JOBPATH'] = '{PROJDIR}/{ANALYZER}/{DATASET}'.format(**self.configurations)


    def run(self):
        print('[INFO   ] Job directory: %s' % self.configurations['JOBPATH'])

        self.make_dirs()

        self.write_tgz()

        self.write_exe()

        while self.jobid <= self.njobs:
            if self.jobid == 1:
                self.write_jdl()
            else:
                self.append_jdl()
            self.jobid += 1
            self.configurations['JOBID'] = self.jobid

        self.submit()
        return

    def make_dirs(self):
        commands = \
'''mkdir -p {JOBPATH}/
rm -rf {JOBPATH}/*
mkdir {JOBPATH}/{JOBNAME}/ {JOBPATH}/{LOGNAME}/
cp {SRCDIR}/{DATASET}.txt {SOURCEFILE}
'''.format(**self.configurations)
        commands = commands.strip().split('\n')

        #print commands
        for cmd in commands:
            subprocess.call(cmd, shell=True)
        return

    def write_tgz(self):
        cfgOutputName = None
        with UserTarball(name=self.tarFilename, logger=self.logger, config=self.config) as tb:
            inputFiles = [re.sub(r'^file:', '', file) for file in getattr(self.config.JobType, 'inputFiles', [])]
            tb.addFiles(userFiles=inputFiles, cfgOutputName=cfgOutputName)
        return

    def write_jdl(self):
        writeme = \
'''Universe                = vanilla
Notification            = Error
Executable              = {JOBNAME}/{EXECUTABLE}
Arguments               = {SOURCEFILE} {MAXEVENTS} {DATASET} {DATASETGROUP} {SELECTION} {PERIOD} {JOBID}
Transfer_Input_Files    = ../../../{TARBALL},{JOBNAME}/{EXECUTABLE},{JOBNAME}/{SOURCEFILE}
Output                  = {LOGNAME}/{JOBNAME}_$(Cluster)_$(Process).stdout
Error                   = {LOGNAME}/{JOBNAME}_$(Cluster)_$(Process).stderr
Log                     = {LOGNAME}/{JOBNAME}_$(Cluster)_$(Process).out
Requirements            = (OpSys == "LINUX") && (Arch != "DUMMY")
request_disk            = 1000000
request_memory          = 199
notify_user             = NOBODY@FNAL.GOV
x509userproxy           = $ENV(X509_USER_PROXY)
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Queue 1
'''.format(**self.configurations)

        #print writeme
        with open(self.jdlFilename, 'w') as f:
            f.write(writeme)
        return

    def append_jdl(self):
        writeme = \
'''Arguments               = {SOURCEFILE} {MAXEVENTS} {DATASET} {DATASETGROUP} {SELECTION} {PERIOD} {JOBID}
Queue 1
'''.format(**self.configurations)

        #print writeme
        with open(self.jdlFilename, 'a') as f:
            f.write(writeme)
        return

    def write_exe(self):

        writeme = \
'''#!/bin/bash

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"
echo ">>> pwd: `pwd`"

JOBID=$7

MACHINE={MACHINE}
if [ $MACHINE == "cmslpc" ]; then
    export SCRAM_ARCH=slc6_amd64_gcc491
    export CMSSW_VERSION=CMSSW_7_4_14
    source /cvmfs/cms.cern.ch/cmsset_default.sh
else
    export SCRAM_ARCH=slc6_amd64_gcc491
    export CMSSW_VERSION=CMSSW_7_4_12
    source /software/tier3/osg/cmsset_default.sh
fi

scram project CMSSW $CMSSW_VERSION
cd $CMSSW_VERSION
eval `scram runtime -sh`

tar -xzf ../default.tgz

cd -

echo ">>> ls -l:"
ls -l

cat <<EOF >> modify_source_file.py
# Modify input source
def chunk(l,n,i):
    m=(len(l)+n-1)/n
    return l[(i-1)*m:i*m]
def read_lines(f):
    return tuple(open(f))
def write_lines(f,lines):
    with open(f, 'w') as w:
        w.write(''.join(lines))
n={NJOBS}
i=$JOBID
f='{SOURCEFILE}'
write_lines(f, chunk(read_lines(f),n,i))
EOF

python modify_source_file.py

echo ">>> cat {SOURCEFILE}:"
cat {SOURCEFILE}

echo ">>> which {ANALYZER}:"
which {ANALYZER}

echo ">>> {ANALYZER} $@"
{ANALYZER} $@

echo ">>> ls -l:"
ls -l

rm input.txt modify_source_file.py
echo "Job finished on host `hostname` on `date`"
'''.format(**self.configurations)

        #print writeme
        with open(self.exeFilename, 'w') as f:
            f.write(writeme)
        os.chmod(self.exeFilename, 0744)
        return

    def submit(self):
        commands = \
'''mv {EXECUTABLE} {SOURCEFILE} {JOBAD} {JOBPATH}/{JOBNAME}/
cd {JOBPATH}/ && condor_submit {JOBNAME}/{JOBAD}
'''.format(**self.configurations)
        commands = commands.strip().split('\n')

        #print commands
        for cmd in commands:
            subprocess.call(cmd, shell=True)
        return

def main():

    if len(sys.argv) < 7:
        raise Exception('Expect 6 command line arguments, received %i' % (len(sys.argv)-1))

    print('[INFO   ] Command line arguments: %s' % (' '.join(sys.argv[1:])))

    logging.basicConfig(filename='jobify.log',level=logging.DEBUG)
    logger = logging.getLogger()

    jobconfig = config()
    jobconfig.JobType.pluginName = 'Analysis'

    jobconfig.User.analyzer     = sys.argv[1]
    jobconfig.User.dataset      = sys.argv[2]
    jobconfig.User.datasetgroup = sys.argv[3]
    jobconfig.User.selection    = sys.argv[4]
    jobconfig.User.period       = sys.argv[5]
    jobconfig.User.njobs        = int(sys.argv[6])

    job = BasicJobType(jobconfig, logger)

    if not os.path.exists('{SRCDIR}/{DATASET}.txt'.format(**job.configurations)):
        raise Exception('Cannot find source file: %s' % ('{SRCDIR}/{DATASET}.txt'.format(**job.configurations)))

    job.run()


if __name__ == '__main__':

    main()

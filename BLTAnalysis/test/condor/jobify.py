#!/usr/bin/env python

"""
A condor job submitter. This script prepares the input source file, the job
config, and the executable.
"""

import os
import sys
import subprocess

class CondorJobType(object):

    def __init__(self):

        self.projdir        = 'blt_projects'
        self.srcdir         = 'sourceFiles'
        self.tgz_name       = 'default.tgz'
        self.txt_name       = 'input.txt'
        self.jdl_name       = 'blt.jdl'
        self.exe_name       = 'blt.sh'
        self.chk_name       = '.checkfile.txt'
        self.jobid          = 1

        self.analyzer       = sys.argv[1]
        self.dataset        = sys.argv[2]
        self.datasetgroup   = sys.argv[3]
        self.selection      = sys.argv[4]
        self.period         = sys.argv[5]
        self.njobs          = int(sys.argv[6])

        self.config = {
            'PROJDIR'      : self.projdir,
            'SRCDIR'       : self.srcdir,
            'TARBALL'      : self.tgz_name,
            'SOURCEFILE'   : self.txt_name,
            'JOBAD'        : self.jdl_name,
            'EXECUTABLE'   : self.exe_name,
            'CHECKFILE'    : self.chk_name,
            'JOBID'        : self.jobid,

            'ANALYZER'     : self.analyzer,
            'DATASET'      : self.dataset,
            'DATASETGROUP' : self.datasetgroup,
            'SELECTION'    : self.selection,
            'PERIOD'       : self.period,
            'NJOBS'        : self.njobs,

            'MACHINE'      : 'cmslpc',
            'MAXEVENTS'    : '-1',
            'JOBNAME'      : 'job',
            'LOGNAME'      : 'res',
        }

        self.config['JOBPATH'] = '{PROJDIR}/{ANALYZER}/{DATASET}'.format(**self.config)

    def safety_check(self):
        if not os.path.exists('{SRCDIR}/{DATASET}.txt'.format(**self.config)):
            raise Exception('Cannot find source file: {SRCDIR}/{DATASET}.txt'.format(**job.configurations))

        if not os.path.exists('{TARBALL}'.format(**self.config)):
            raise Exception('Cannot find tarball: {TARBALL}'.format(**self.config))

    def run(self):
        self.make_dirs()

        self.write_exe()

        while self.jobid <= self.njobs:
            if self.jobid == 1:
                self.write_jdl()
            else:
                self.append_jdl()
            self.jobid += 1
            self.config['JOBID'] = self.jobid

        self.add_check()

        self.submit_jobs()
        return

    def make_dirs(self):
        commands = \
'''mkdir -p {JOBPATH}/
rm -rf {JOBPATH}/*
mkdir {JOBPATH}/{JOBNAME}/ {JOBPATH}/{LOGNAME}/
cp {SRCDIR}/{DATASET}.txt {SOURCEFILE}
'''.format(**self.config)
        commands = commands.strip().split('\n')

        #print commands
        for cmd in commands:
            subprocess.call(cmd, shell=True)
        return

    def write_jdl(self):

        if 'X509_USER_PROXY' not in os.environ:
            myproxy = '/tmp/x509up_u%s' % subprocess.check_output('id -u', shell=True)
            os.environ['X509_USER_PROXY'] = myproxy.rstrip('\n')

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
request_disk            = 2000000
request_memory          = 1024
#notify_user             = NOBODY@FNAL.GOV
x509userproxy           = $ENV(X509_USER_PROXY)
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Queue 1
'''.format(**self.config)

        #print writeme
        with open(self.jdl_name, 'w') as f:
            f.write(writeme)
        return

    def append_jdl(self):
        writeme = \
'''Arguments               = {SOURCEFILE} {MAXEVENTS} {DATASET} {DATASETGROUP} {SELECTION} {PERIOD} {JOBID}
Queue 1
'''.format(**self.config)

        #print writeme
        with open(self.jdl_name, 'a') as f:
            f.write(writeme)
        return

    def write_exe(self):

        writeme = \
'''#!/bin/bash

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

RUNTIME_AREA=`pwd`
echo ">>> RUNTIME_AREA=$RUNTIME_AREA"

JOBID=$7
echo ">>> JOBID=$JOBID"

MACHINE={MACHINE}
echo ">>> MACHINE=$MACHINE"

if [ $MACHINE == "cmslpc" ]; then
    export SCRAM_ARCH=slc6_amd64_gcc491
    export CMSSW_VERSION=CMSSW_7_4_14
    source /cvmfs/cms.cern.ch/cmsset_default.sh
else
    export SCRAM_ARCH=slc6_amd64_gcc491
    export CMSSW_VERSION=CMSSW_7_4_12
    source /software/tier3/osg/cmsset_default.sh
fi

echo ">>> X509_USER_PROXY=$X509_USER_PROXY"

# Setup CMSSW environment
scram project CMSSW $CMSSW_VERSION
cd $CMSSW_VERSION
eval `scram runtime -sh`

SOFTWARE_DIR=`pwd`
echo ">>> SOFTWARE_DIR=$SOFTWARE_DIR"

# Extract files
tar -xzf ../default.tgz

echo ">>> Listing SOFTWARE_DIR"
ls -Al $SOFTWARE_DIR

# Return to working directory
cd $RUNTIME_AREA

echo ">>> Listing RUNTIME_AREA"
ls -Al $RUNTIME_AREA

# Modify source file
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

# Run the analyzer
echo ">>> {ANALYZER} $@"
{ANALYZER} $@

# Done
echo "Exit status is $?"

echo ">>> Listing RUNTIME_AREA"
ls -Al $RUNTIME_AREA

rm input.txt modify_source_file.py
echo "Job finished on host `hostname` on `date`"
'''.format(**self.config)

        #print writeme
        with open(self.exe_name, 'w') as f:
            f.write(writeme)
        os.chmod(self.exe_name, 0744)  # make executable
        return

    def add_check(self):
        commands = \
'''echo {JOBPATH} {NJOBS} `wc -l < {SOURCEFILE}` >> {CHECKFILE}
'''.format(**self.config)
        commands = commands.strip().split('\n')

        #print commands
        for cmd in commands:
            subprocess.call(cmd, shell=True)
        return

    def submit_jobs(self):
        commands = \
'''mv {EXECUTABLE} {SOURCEFILE} {JOBAD} {JOBPATH}/{JOBNAME}/
cd {JOBPATH}/ && condor_submit {JOBNAME}/{JOBAD}
'''.format(**self.config)
        commands = commands.strip().split('\n')

        #print commands
        for cmd in commands:
            subprocess.call(cmd, shell=True)
        return


# ______________________________________________________________________________
def main():

    if len(sys.argv) < 7:
        raise Exception('Expect 6 command line arguments, received %i' % (len(sys.argv)-1))

    print('[INFO   ] Creating condor jobs ...')
    print('[INFO   ] %sPlease make sure your GRID proxy is valid!%s' % ('\033[93m', '\033[0m'))
    print('[INFO   ] Command line arguments: %s' % (' '.join(sys.argv[1:])))

    job = CondorJobType()
    job.safety_check()

    print('[INFO   ] Job directory: %s' % job.config['JOBPATH'])

    job.run()

    print('[INFO   ] %s%i jobs are submitted to condor.%s' % ('\033[92m', job.config['NJOBS'], '\033[0m'))


# ______________________________________________________________________________
if __name__ == '__main__':

    main()

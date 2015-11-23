#!/usr/bin/env python

import os
import re
import glob
import subprocess

class Retriever(object):

    def __init__(self, chk='.checkfile.txt', logname='res'):

        self.chk = chk
        self.logname = logname
        self.datasets = []
        self.target = ''

    def run(self):

        self.setup()
        self.hadd()

    def setup(self):

        if not os.path.exists(self.chk):
            raise Exception('Cannot find check file %s?!' % self.chk)

        with open(self.chk) as f:
            for line in f.readlines():
                dataset, njobs, nfiles = line.rstrip('\n').split(' ')
                self.datasets.append((dataset, int(njobs), int(nfiles)))

        # Safety check
        for dataset in self.datasets:
            logdir = os.path.join(dataset[0], self.logname)
            out = subprocess.check_output('cd %s && grep "Exit status is 0" *.stdout | wc -l' % (logdir), shell=True)
            out = int(out)
            if out != dataset[1]:
                raise Exception('Expect %i jobs be completed, found %i' % (dataset[1], out))

            out = subprocess.check_output('cd %s && grep "Processing.*\.root" *.stdout | wc -l' % (logdir), shell=True)
            out = int(out)
            if out != dataset[2]:
                raise Exception('Expect %i input files be processed, found %i' % (dataset[2], out))

    def hadd(self):

        files = []
        for dataset in self.datasets:
            files.append('%s/*.root' % dataset[0])

        if not files:
            raise Exception('No files to be added.')

        # Decide target file name using the first source file name
        files_tmp = []
        for f in files:
            files_tmp += glob.glob(f)

        if not files_tmp:
            raise Exception('No files to be added.')

        m = re.match('(.*?)_.*\.root', os.path.basename(files_tmp[0]))
        if not m:
            raise Exception('Cannot parse source file name such as this: %s' % (os.path.basename(files_tmp[0])))

        self.target = "%s.root" % (m.group(1))

        print('[INFO   ] Calling the following')
        print('hadd -f %s %s' % (self.target, ' '.join(files)))
        subprocess.call('hadd -f %s %s' % (self.target, ' '.join(files)), shell=True)


# ______________________________________________________________________________
def main():

    print('[INFO   ] Retrieving condor jobs ...')

    ret = Retriever()
    ret.run()

    print('[INFO   ] %s%s is created (%iM).%s' % ('\033[92m', ret.target, os.stat(ret.target).st_size >> 20, '\033[0m'))

# ______________________________________________________________________________
if __name__ == '__main__':

    main()

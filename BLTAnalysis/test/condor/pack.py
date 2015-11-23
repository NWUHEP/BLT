#!/usr/bin/env python

"""
Based on CERN CRAB3 UserTarball.py
"""

import os
import sys
import tarfile
import glob

class UserTarball(object):

    def __init__(self, name='default.tgz', mode='w:gz', chk='.checkfile.txt'):
        self.name = name
        self.mode = mode
        self.chk = chk
        self.directories = ['lib', 'module', 'bin']
        self.dataDirs = ['data']
        self.dataDirsRoot = ['BaconAna', 'BLT']

    def run(self):

        self.add_files()

        self.add_check()

    def add_files(self, userFiles=None):
        userFiles = userFiles or []

        with tarfile.open(self.name, self.mode) as tar:

            # Tar up whole directories
            for directory in self.directories:
                fullPath = os.path.join(os.environ['CMSSW_BASE'], directory)
                if os.path.exists(fullPath):
                    tar.add(fullPath, directory, recursive=True)

            # Search for and tar up "data" directories in src/
            srcPath = os.path.join(os.environ['CMSSW_BASE'], 'src')
            for root, dirs, files in os.walk(srcPath):
                # Prune the directories visited by os.walk
                if os.path.basename(root) == 'src':
                    dirs[:] = [d for d in dirs if d in self.dataDirsRoot]
                dirs[:] = [d for d in dirs if d != '.git']

                if os.path.basename(root) in self.dataDirs:
                    directory = root.replace(srcPath,'src')
                    tar.add(root, directory, recursive=True)

            # Tar up extra files the user needs
            for globName in userFiles:
                fileNames = glob.glob(globName)
                if not fileNames:
                    raise Exception("The input file '%s' cannot be found." % globName)
                for filename in fileNames:
                    directory = os.path.basename(filename)
                    tar.add(filename, directory, recursive=True)

    def add_check(self):
        # Make a blank file
        open(self.chk, 'w').close()


# ______________________________________________________________________________
def main():

    print('[INFO   ] Packing tarball ...')
    print('[INFO   ] Using CMSSW_BASE: %s' % (os.environ['CMSSW_BASE']))

    tb = UserTarball()
    tb.run()

    print('[INFO   ] %s%s is created (%iM).%s' % ('\033[92m', tb.name, os.stat(tb.name).st_size >> 20, '\033[0m'))


# ______________________________________________________________________________
if __name__ == '__main__':

    main()

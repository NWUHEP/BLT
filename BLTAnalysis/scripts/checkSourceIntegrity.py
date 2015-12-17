#!/usr/bin/env python

import os
import re
import glob

outLFNPrefix = '/eos/uscms/'
outLFNDirBase = '/store/group/lpchzg/makingBacon/01/%s/'
outputFile = 'Output_%d.root'


projects = [
    ('DoubleEG_Run2015D-PromptReco-v4'  , 6265),
    ('DoubleEG_Run2015D-05Oct2015-v1'   , 4242),
    ('DoubleMuon_Run2015D-PromptReco-v4', 3382),
    ('DoubleMuon_Run2015D-05Oct2015-v1' , 2114),
]

def check_project(proj, nfiles):
    flags = [0 for i in xrange(nfiles)]

    files = glob.glob((outLFNPrefix+outLFNDirBase+'*/*/*/*/*.root') % proj)
    total_size = 0

    for f in files:
        s = os.lstat(f)
        if s.st_size > 100:
            pattern = '.*' + outputFile.replace('%d', '(\d+)')
            m = re.match(pattern, f)
            if m:
               jobid = int(m.group(1)) - 1  # jobid starts from 1
               flags[jobid] = 1
        total_size += s.st_size

    nfiles_success = sum(flags)
    percent_success = float(nfiles_success) / float(nfiles) * 100.

    print("Project    : {0}".format(proj))
    print("Total size : {0:8d}G".format(total_size >> 30))
    print("Jobs status:	{0:4.2f}% ({1:4d}/{2:4d})".format(percent_success, nfiles_success, nfiles))
    for i in xrange(nfiles):
        if not flags[i]:
            print("{0:11s}   failed job {1:4d}".format("", i + 1))  # jobid starts from 1
    print("")

def check():
    for proj, nfiles in projects:
        check_project(proj, nfiles)

# ______________________________________________________________________________
if __name__ == '__main__':

    print("Checking... be patient")
    check()

    print("DONE checking.")


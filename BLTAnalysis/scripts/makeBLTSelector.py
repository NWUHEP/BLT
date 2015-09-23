#!/usr/bin/env python

import sys
from ROOT import gSystem, TFile

def main():

    # __________________________________________________________________________
    # Parse input arguments
    if len(sys.argv) < 2:
        raise Exception("Usage: %s file_name [tree_name]" % sys.argv[0])

    fname = sys.argv[1]
    tname = "Events"
    if len(sys.argv) > 2:
        tname = sys.argv[2]

    # __________________________________________________________________________
    # Load Bacon library
    libname = "libBaconAnaDataFormats.so"
    ret = gSystem.Load(libname)
    if ret != 0:
        raise Exception("ERROR: Failed to load %s." % libname)

    # __________________________________________________________________________
    # Retrieve the TTree
    infile = TFile.Open(fname)
    if infile == None:
        raise Exception("ERROR: Failed to open %s." % fname)

    tree = infile.Get(tname)
    if tree == None:
        raise Exception("ERROR: Failed to get %s." % tname)

    # Get the branches
    branches = tree.GetListOfBranches()
    nb = branches.GetEntriesFast()

    blist = []
    for i in xrange(nb):
        branch = branches.UncheckedAt(i)
        blist.append((branch.GetName(), branch.GetClassName()))

    # __________________________________________________________________________
    # Now print
    sep = "-" * 40
    fmt_branch = lambda n, c: "b_"+n+("Arr" if "Array" in c else "")
    fmt_member = lambda n, c: "f"+n[0:1].upper()+n[1:]+("Arr" if "Array" in c else "")

    print "number of branches: %i" % nb
    print sep
    print "@ Member data"
    for n, c in blist:
        print "    {0:24}*{1};".format(c, fmt_member(n, c))
    print
    for n, c in blist:
        print "    {0:24}*{1};".format("TBranch", fmt_branch(n, c))
    print sep
    print "@ Init()"
    for n, c in blist:
        print "    {0:24} = 0;".format(fmt_member(n, c))
    print
    for n, c in blist:
        print "    fChain->SetBranchAddress(\"{0}\", &{1}, &{2});".format(n, fmt_member(n, c), fmt_branch(n, c))
    print sep

if __name__=="__main__":
    main()

BLT
===

BLT is a framework based on the ROOT TSelector class specifically for working with [bacon ntuples](https://github.com/ksung25/BaconProd) developed and maintained by Kevin Sung and Phil Harris. 

Setup
=====

The master version has been tested to run with CMSSW_7_4_14 with ROOT6.  This will work cmslpc using scientific linux 6.

### Build and source environment

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
cmsenv
```

### Checkout dependencies

These forks are sync'ed to Kevin's repositories on 2015/11/17

```
git clone git@github.com:NWUHEP/BaconAna.git
```

Depending on which branch of BLT you are working with you will need to check out a specific tag of BaconAna.  These are listed in the following table:

| BLT branch | BaconAna tag/branch |
|:---:|:---:|
| master    | -- |
| amumu_2012 | 04 |
| amumu_2016 | master |

To check out a tag do the following from the top of the BaconAna repository,

```
git checkout tags/<tag>
```

### Checkout and compile BLT code

```
git clone git@github.com:NWUHEP/BLT.git
scram b -j 12
```

### (Optional) Sync forked repositories

```
cd BaconProd && git checkout master && git fetch upstream && git merge upstream/master && git push origin master && cd -
cd BaconAna && git checkout master && git fetch upstream && git merge upstream/master && git push origin master && cd -
```

## Running the demo analyzer

There are two example analyzers included in the BLT repository.  These produce a set of skimmed ntuples...

7 input arguments are mandatory: [input file] [no of events] [selection] [dataset] [datasetgroup] [period] [jobid]

```
cd BLT/BLTAnalysis/test
DemoAnalyzer Output.root 1000 DYJetsToLL_M-50 DYJetsToLL mumu 2015 0
```

or alternatively you can provide a text file with a list of input files,

```
DemoAnalyzer input.txt 1000 DYJetsToLL_M-50 DYJetsToLL mumu 2015 0
```

## Running a BLT analyzer with condor

The batch submission is handled by the BatchMaster class.  

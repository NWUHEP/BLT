BLT
===

BLT is a framework for analyzing bacon ntuples originally developed by @jiafulow.

Setup
=====

The master version has been tested to run with CMSSW_7_4_14 with ROOT6.  This will work cmslpc using scientific linux 6.

# Build and source environment

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_14
cd CMSSW_7_4_14/src
cmsenv
```

# Checkout dependencies

These forks are sync'ed to Kevin's repositories on 2015/11/17

```
git clone git@github.com:NWUHEP/BaconProd.git
cd BaconProd && git remote add upstream git@github.com:ksung25/BaconProd.git && cd -
git clone git@github.com:NWUHEP/BaconAna.git
cd BaconAna && git remote add upstream git@github.com:ksung25/BaconAna.git && cd -
```

If you plan on producing ntuples using the bacon framework, you will need to checkout some additional CMSSW dependencies,

```
source BaconProd/scripts/setup_prod.sh
# If copying via AFS is too slow, use scp instead:
#   scp -r -C <USERNAME>@lxplus.cern.ch:/afs/cern.ch/work/k/ksung/public/Development/Run2Packages/* ./
```

# Checkout and compile BLT code

```
git clone git@github.com:jiafulow/BLT.git
scram b -j 12
```

# Make Bacon ntuples

```
cd BLT/BLTAnalysis/test
cmsRun makingBacon_MC_25ns_MINIAOD.py 
```

# (Optional) Sync forked repositories

```
cd BaconProd && git checkout master && git fetch upstream && git merge upstream/master && git push origin master && cd -
cd BaconAna && git checkout master && git fetch upstream && git merge upstream/master && git push origin master && cd -
```

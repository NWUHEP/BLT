#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
DATANAME=$1
COUNT=$2

### Specify addtional arguments here ####
SUFFIX=$3
SELECTION=$4
PERIOD=$5

### Transfer files, prepare directory ###
export SCRAM_ARCH=slc6_amd64_gcc491
export CMSSW_VERSION=CMSSW_7_4_12
source /software/tier3/osg/cmsset_default.sh

# Setup CMSSW environment
scram project CMSSW $CMSSW_VERSION
cd $CMSSW_VERSION/src
eval `scram runtime -sh`

echo $PWD

cp ../../source.tar.gz .
tar -xzf source.tar.gz
ls
cd BLT/BLTAnalysis/test
cp ../../../../input_${DATANAME}_${COUNT}.txt input.txt

### Run the analyzer
#DimuonAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

### Copy output and cleanup ###
#rename .root _${DATANAME}_$COUNT.root histos/fcncHistograms*
#cp ouput_*.root ${_CONDOR_SCRATCH_DIR}

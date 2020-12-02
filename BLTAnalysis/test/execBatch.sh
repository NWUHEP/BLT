#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
DATANAME=$1
COUNT=$2
SUFFIX=$3
SELECTION=$4
PERIOD=$5

### Transfer files, prepare directory ###
TOPDIR=$PWD

# lpc
export SCRAM_ARCH=slc7_amd64_gcc820
export CMSSW_VERSION=CMSSW_11_1_2_patch3
source /cvmfs/cms.cern.ch/cmsset_default.sh

# nut3
#export SCRAM_ARCH=slc6_amd64_gcc530
#export CMSSW_VERSION=CMSSW_8_0_24_patch1
#source /software/tier3/osg/cmsset_default.sh 

# Setup CMSSW environment
eval `scramv1 project CMSSW $CMSSW_VERSION`
cd $CMSSW_VERSION

cp $TOPDIR/source.tar.gz .
tar -xzf source.tar.gz
cd src
eval `scramv1 runtime -sh`

cd BLT/BLTAnalysis/test
cp $TOPDIR/input_${DATANAME}_${COUNT}.txt input.txt

echo $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
echo $PATH
pwd
cat input.txt
### Run the analyzer
#DimuonAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#MultileptonAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#zjpsiAnalyzerV2 input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

hzgAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#hltPrintoutAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#hzgAnalyzerUCSBSync input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

#gbrTrainAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#upsilonGammaAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

### Copy output and cleanup ###
cp output_${DATANAME}_${COUNT}.root ${_CONDOR_SCRATCH_DIR}

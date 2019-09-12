#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
DATANAME=$1
COUNT=$2
SUFFIX=$3
SELECTION=$4
PERIOD=$5
OUTDIR=$6

### Transfer files, prepare directory ###
TOPDIR=$PWD

# lpc
export SCRAM_ARCH=slc6_amd64_gcc530
export CMSSW_VERSION=CMSSW_8_0_20
source /cvmfs/cms.cern.ch/cmsset_default.sh

# temporary fix
tar -xzf source.tar.gz
mv $CMSSW_VERSION tmp_source
scram project CMSSW $CMSSW_VERSION
cp -r tmp_source/src/* $CMSSW_VERSION/src
cd $CMSSW_VERSION/src
eval `scram runtime -sh`

# this used to work, now it don't
#tar -xzf source.tar.gz
#cd $CMSSW_VERSION/src/
#scramv1 b ProjectRename

cmsenv
scramv1 b -j8 #ProjectRename
cd BLT/BLTAnalysis/scripts
cp $TOPDIR/input_${DATANAME}_${COUNT}.txt input.txt

echo $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
echo $PATH
pwd
cat input.txt

### Run the analyzer
#MultilepAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#ZTauTauAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#TauAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
ZEEAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

### Copy output and cleanup ###
# cp output_${DATANAME}_${COUNT}.root ${_CONDOR_SCRATCH_DIR}



### Copy output and cleanup ###
FILE=output_${DATANAME}_${COUNT}.root
xrdcp -f ${FILE} ${OUTDIR}/${FILE} 2>&1
XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
  rm *.root
  echo "exit code $XRDEXIT, failure in xrdcp"
  exit $XRDEXIT
fi
rm ${FILE}


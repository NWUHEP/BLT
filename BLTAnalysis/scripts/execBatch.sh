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
echo '----- ls TOPDIR -----'
ls

echo '----- untar source.tar.gz -----'
tar -xzf source.tar.gz

echo '----- rename source as tmp_source -----'
mv $CMSSW_VERSION tmp_source

echo '----- ls tmp_source/src -----'
ls tmp_source/src

echo '----- scram project -----'
scram project CMSSW $CMSSW_VERSION
#cmsrel CMSSW_8_0_20
cp -r tmp_source/src/* $CMSSW_VERSION/src
cd $CMSSW_VERSION/src
eval `scram runtime -sh`

# this used to work, now it don't
#tar -xzf source.tar.gz
#cd $CMSSW_VERSION/src/
#scramv1 b ProjectRename

cmsenv
scram b -j8
cd BLT/BLTAnalysis/scripts
cp $TOPDIR/input_${DATANAME}_${COUNT}.txt input.txt

echo $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
echo $PATH
pwd
cat input.txt

### Run the analyzer
MultilepAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#DileptonSelector input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#FakeSelector input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

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


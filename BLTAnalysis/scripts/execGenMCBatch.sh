#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
NARGS=$(($#))
echo "nArgs: $NARGS"
ARGUMENTS=${@:1:$NARGS-2}
COUNT="${@:(-2):1}"
OUTDIR=${@:$NARGS:$NARGS}
SCRIPT=$1

### Transfer files, prepare directory ###
TOPDIR=$PWD

# lpc
export SCRAM_ARCH=slc6_amd64_gcc530
# # export CMSSW_VERSION=CMSSW_8_0_20
source /cvmfs/cms.cern.ch/cmsset_default.sh
tar -xzf source.tar.gz
cd CMSSW_?_?_*/src

cmsenv
scramv1 b -j8 #ProjectRename

echo "topdir: $TOPDIR"
echo "script+args: $ARGUMENTS"
echo "script: $SCRIPT"
echo "count: $COUNT"
echo "outdir: $OUTDIR"
pwd

echo "Appending process.RandomNumberGeneratorService.generator.initialSeed = ${COUNT} to $SCRIPT"
echo "process.RandomNumberGeneratorService.generator.initialSeed = ${COUNT}" >> ${SCRIPT}

### Run the script
cmsRun $ARGUMENTS

### Copy output and cleanup ###
FILE=$(ls *.root | head -n 1)
FBASE=${FILE:0:${#FILE}-5}
echo "File: $FILE"
echo "FileBase: $FBASE"
echo "xrdcp -f $FILE ${OUTDIR}/${FBASE}_${COUNT}.root 2>&1"
xrdcp -f $FILE ${OUTDIR}/${FBASE}_${COUNT}.root 2>&1
XRDEXIT=$?
if [[ $XRDEXIT -ne 0 ]]; then
  rm *.root
  echo "exit code $XRDEXIT, failure in xrdcp"
  exit $XRDEXIT
fi
rm ${FILE}


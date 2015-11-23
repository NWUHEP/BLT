#!/bin/bash

python pack.py

python jobify.py DemoAnalyzer DYJetsToLL_M-50 DYJetsToLL mumu 2015 50
python jobify.py DemoAnalyzer TTJets          TTJets     mumu 2015 50


# When the condor jobs are done, do this:
#python retrieve.py

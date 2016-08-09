#!/bin/bash

python pack.py

python jobify.py DemoAnalyzer DYJetsToLL_M-50 DYJetsToLL mumu 2015 50
#python jobify.py DemoAnalyzer TTJets          TTJets     mumu 2015 50
#
#python jobify.py DemoAnalyzer DoubleMuon_Run2015D-PromptReco-v4 DoubleMuon mumu 2015 50 
#python jobify.py DemoAnalyzer DoubleMuon_Run2015D-05Oct2015-v1  DoubleMuon mumu 2015 50 
#
#python jobify.py DemoAnalyzer DoubleEG_Run2015D-PromptReco-v4   DoubleEG   ee   2015 50 
#python jobify.py DemoAnalyzer DoubleEG_Run2015D-05Oct2015-v1    DoubleEG   ee   2015 50 

# When the condor jobs are done, do this:
#python retrieve.py

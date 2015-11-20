#!/bin/bash

python pack.py

python jobify.py DemoAnalyzer DYJetsToLL_M-50 DYJetsToLL mumu 2015 2
python jobify.py DemoAnalyzer TTJets          TTJets     mumu 2015 2

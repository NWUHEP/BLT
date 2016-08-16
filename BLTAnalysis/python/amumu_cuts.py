#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def stack_by_group(df, by, var, histtype='stepfilled', bins=30., range=None):

    g = df.groupby(by)
    data = [g.get_group(group)[var] for group in g.groups.keys()] 
    plt.hist(data, histtype=histtype, bins=bins, stacked=True)

if __name__ == '__main__':

    if len(sys.argv) > 2:
        infile  = sys.argv[1]
        cat     = sys.argv[2]
    else:
        infile = 'data/ntuple_dimuon.csv'
        cat     = '1b1f'

    data = pd.read_csv(infile)

    # Get the reference data
    if cat == '1b1f':
        data_1b1f = pd.read_csv('data/events_pf_1stbump.txt')
        en = data_1b1f.event_number.values
        data['1b1f'] = data.event_number.isin(en)
    elif cat == '1b1c':
        data_1b1f = pd.read_csv('data/events_pf_2ndbump.txt')
        en = data_1b1f.event_number.values
        data['1b1c'] = data.event_number.isin(en)

    print data.shape[0], data.event_number.isin(en).value_counts()[True]

    cut4a = 'n_bjets == 1'
    data_cut = data.query(cut4a)
    print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]

    if cat == '1b1f':
        cut4b = 'n_jets == 0'
        data_cut = data_cut.query(cut4b)
        print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]

        cut5 = 'jet_pt > 30 and abs(jet_eta) > 2.4'
        data_cut = data_cut.query(cut5)
        print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]

    elif cat == '1b1c':
        cut4b = 'n_jets == 1'
        data_cut = data_cut.query(cut4b)
        print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]

        cut5 = 'jet_pt > 30 and abs(jet_eta) <= 2.4'
        data_cut = data_cut.query(cut5)
        print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]

        cut6 = 'met_mag < 40'
        data_cut = data_cut.query(cut6)
        print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]
    
        cut7 = 'delta_phi > 2.5'
        data_cut = data_cut.query(cut7)
        print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]
    
    #cut6 = '26 < dimuon_mass < 32'
    #data_cut = data_cut.query(cut6)
    #print data_cut.shape[0], data_cut.event_number.isin(en).value_counts()[True]
    

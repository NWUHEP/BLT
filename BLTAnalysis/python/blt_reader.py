#!/usr/bin/env python

import os, sys

import numpy as np
import pandas as pd
import ROOT as r
import json
from collections import OrderedDict

from tqdm import tqdm

'''
Script for getting blt data out of ROOT files and into CSV format.
'''

def make_directory(filePath, clear=True):
    if not os.path.exists(filePath):
        os.system('mkdir -p '+filePath)

    if clear and len(os.listdir(filePath)) != 0:
        os.system('rm '+filePath+'/*')

def calculate_cos_theta(ref_p4, boost_p4, target_p4):
    '''
    Calculate cos(theta) of target in rest frame of boost relative with ref
    defining the x axis in boost's reference frame.

    Paramters:
    ==========
    ref_p4: TLorentzVector of reference particle
    boost_p4: TLorentzVector of particle defining the rest frame
    target_p4: TLorentzVector of particle whose cos(theta) value we would like to calculate
    '''
                                                                                  
    boost = -1*boost_p4.BoostVector() 
    B_clone.Boost(boost)  
    ref = ref_p4.Vect()                                        
    beam_axis = r.TLorentzVector(0, 0, 1, 1)
    beam_axis.Boost(boost)                                                        
                                                                                  
    axis_z = (-1*B_threevect).Unit() 
    axis_y = beamAxis.Vect().Cross(B_threevect).Unit()               
    axis_x = axis_y.Cross(axis_z).Unit()                       
    rotation = rotation.RotateAxes(axis_x, axis_y, axis_z).Inverse()    
                                                                                  
    target_p4.Boost(boost)                                                  
    target_p4.Transform(rotation) 
                                                                                  
    return target_p4.CosTheta() 


if __name__ == '__main__':

    ### Configuration ###
    selection    = 'mumu'
    period       = 2012
    infile       = '/tthome/naodell/batch_output/output_{0}_{1}.root'.format(selection, period)
    output_path  = '/tthome/share/noobs/flattuples/{0}_{1}'.format(selection, period)

    if period == 2016:
        dataset_list = [
                        'muon_2016B', 'muon_2016C', 'muon_2016D', #'muon_2016E', 'muon_2016F', 
                        'ttjets', 't_t', 't_tw', 'tbar_t', 'tbar_tw', 
                        'zjets_m-50', 'zjets_m-10to50',
                        ]
    elif period == 2012:
        dataset_list = [
                        'muon_2012A', 'muon_2012B', 'muon_2012C', 'muon_2012D', 
                        #'electron_2012A', 'electron_2012B', 'electron_2012C', 'electron_2012D', 
                        'ttbar_lep', 'ttbar_semilep',
                        'zjets_m-50', 'zjets_m-10to50',
                        'z1jets_m-50', 'z1jets_m-10to50',
                        'z2jets_m-50', 'z2jets_m-10to50',
                        'z3jets_m-50', 'z3jets_m-10to50',
                        'z4jets_m-50', 'z4jets_m-10to50',
                        't_s', 't_t', 't_tw', 'tbar_s', 'tbar_t', 'tbar_tw', 
                        'ww', 'wz_2l2q', 'wz_3lnu', 'zz_2l2q', 'zz_2l2nu',
                        'bprime_bb_semilep', 'bprime_t-channel', 
                        'fcnc_s-channel', 'fcnc_tt_semilep'
                       ]

    features = [
               'run_number', 'event_number', 'lumi', 'weight', 'trigger_status', 'n_pu',
               'n_jets', 'n_fwdjets', 'n_bjets', 'n_muons', 'n_electrons',

               'lepton1_pt', 'lepton1_eta', 'lepton1_phi', 
               'lepton1_iso', 'lepton1_q', 'lepton1_flavor', 'lepton1_trigger',
               'lepton2_pt', 'lepton2_eta', 'lepton2_phi',  
               'lepton2_iso', 'lepton2_q', 'lepton2_flavor', 'lepton2_trigger',
               'lepton_delta_eta', 'lepton_delta_phi', 'lepton_delta_r',
               'dilepton_mass', 'dilepton_pt', 'dilepton_eta', 'dilepton_phi', 
               'dilepton_pt_over_m',

               'met_mag', 'met_phi',
               'bjet_pt', 'bjet_eta', 'bjet_phi', 'bjet_e', 
               'bjet_d0', 'bjet_tag', 'bjet_flavor',
               'jet_pt', 'jet_eta', 'jet_phi', 'jet_e', 
               'jet_d0', 'jet_tag',  'jet_flavor',
               'jet_delta_eta', 'jet_delta_phi', 'jet_delta_r',
               'dijet_mass', 'dijet_pt', 'dijet_eta', 'dijet_phi', 
               'dijet_pt_over_m',

               'lepton1_b_mass', 'lepton1_b_pt', 
               'lepton1_b_delta_eta', 'lepton1_b_delta_phi', 'lepton1_b_delta_r',
               'lepton2_b_mass', 'lepton2_b_pt', 
               'lepton2_b_delta_eta', 'lepton2_b_delta_phi', 'lepton2_b_delta_r',

               'dilepton_j_mass', 'dilepton_j_pt', 
               'dilepton_j_delta_eta', 'dilepton_j_delta_phi', 'dilepton_j_delta_r',
               'dilepton_b_mass', 'dilepton_b_pt', 
               'dilepton_b_delta_eta', 'dilepton_b_delta_phi', 'dilepton_b_delta_r',
               'four_body_mass',
               'four_body_delta_phi', 'four_body_delta_eta', 'four_body_delta_r',

               #'t_xj', 't_xb', 't_bj',

               #'n_partons'
               #'gen_bjet_pt', 'gen_bjet_eta', 'gen_bjet_phi', 'gen_bjet_e', 
               #'gen_bjet_tag', 'gen_dilepton_b_mass',
               #'gen_jet_pt', 'gen_jet_eta', 'gen_jet_phi', 'gen_jet_e', 
               #'gen_jet_tag', 'gen_dilepton_j_mass'
              ]
    make_directory(output_path, clear=False)

    ### Get input bltuple ###
    print 'Opening file {0}'.format(infile)
    froot  = r.TFile(infile)

    event_count = {}
    for dataset in tqdm(dataset_list, 
                        desc       = 'Unpacking',
                        unit_scale = True,
                        ncols      = 75,
                        total      = len(dataset_list)
                       ):
        ecount = froot.Get('TotalEvents_{0}'.format(dataset))
        if ecount:
            event_count[dataset] = [ecount.GetBinContent(i+1) for i in range(ecount.GetNbinsX())]
        else:
            print 'Could not find dataset {0} in root file...'.format(dataset)
            continue

        tree    = froot.Get('tree_{0}'.format(dataset))
        n       = tree.GetEntriesFast()
        ntuple = OrderedDict([(f,[]) for f in features])

        for i in tqdm(xrange(n), 
                      desc       = dataset,
                      unit_scale = True,
                      ncols      = 75,
                      total      = n,
                      leave      = False
                     ):
            tree.GetEntry(i)

            # get and build physics objects
            lep1, lep2,       = tree.leptonOneP4, tree.leptonTwoP4
            bjet, jet         = tree.bjetP4, tree.jetP4
            gen_bjet, gen_jet = tree.genBJetP4, None, #tree.genJetP4
            met, met_phi      = tree.met, tree.metPhi
            dilepton          = lep1 + lep2
            dijet             = jet + bjet
            lepton1_b         = lep1 + bjet
            lepton2_b         = lep2 + bjet
            dilepton_b        = dilepton + bjet
            dilepton_j        = dilepton + jet
            fourbody          = dilepton + dijet

            ### Blind the 2016 data in the signal region (24 < M_mumu < 32)
            #if period == 2016 \
            #    and selection == 'mumu' \
            #    and abs(dilepton.M() - 28) < 4 \
            #    and dataset.split('_')[0] == 'muon':
            #    continue

            # event info
            ntuple['run_number'].append(tree.runNumber)
            ntuple['event_number'].append(tree.evtNumber)
            ntuple['lumi'].append(tree.lumiSection)
            ntuple['trigger_status'].append(tree.triggerStatus)

            if period == 2012:
                if dataset in ['zjets_m-50', 'zjets_m-10to50'] \
                    and tree.nPartons > 0 \
                    and tree.nPartons < 5:
                    ntuple['weight'].append(0.)
                else:
                    ntuple['weight'].append(tree.eventWeight)
            else:
                ntuple['weight'].append(tree.eventWeight)

            ntuple['n_pu'].append(tree.nPU)
            ntuple['n_muons'].append(tree.nMuons)
            ntuple['n_electrons'].append(tree.nElectrons)
            ntuple['n_jets'].append(tree.nJets)
            ntuple['n_fwdjets'].append(tree.nFwdJets)
            ntuple['n_bjets'].append(tree.nBJets)

            ### lepton
            ntuple['lepton1_pt'].append(lep1.Pt())
            ntuple['lepton1_eta'].append(lep1.Eta())
            ntuple['lepton1_phi'].append(lep1.Phi())
            ntuple['lepton1_q'].append(tree.leptonOneQ)
            ntuple['lepton1_iso'].append(tree.leptonOneIso)
            ntuple['lepton1_flavor'].append(tree.leptonOneFlavor)
            ntuple['lepton1_trigger'].append(tree.leptonOneTrigger)

            ntuple['lepton2_pt'].append(lep2.Pt())
            ntuple['lepton2_eta'].append(lep2.Eta())
            ntuple['lepton2_phi'].append(lep2.Phi())
            ntuple['lepton2_q'].append(tree.leptonTwoQ)
            ntuple['lepton2_iso'].append(tree.leptonTwoIso)
            ntuple['lepton2_flavor'].append(tree.leptonTwoFlavor)
            ntuple['lepton2_trigger'].append(tree.leptonTwoTrigger)

            ### dilepton 
            ntuple['lepton_delta_eta'].append(abs(lep1.Eta() - lep2.Eta()))
            ntuple['lepton_delta_phi'].append(abs(lep1.DeltaPhi(lep2)))
            ntuple['lepton_delta_r'].append(lep1.DeltaR(lep2))

            ntuple['dilepton_mass'].append(dilepton.M())
            ntuple['dilepton_pt'].append(dilepton.Pt())
            ntuple['dilepton_eta'].append(dilepton.Eta())
            ntuple['dilepton_phi'].append(dilepton.Phi())
            ntuple['dilepton_pt_over_m'].append(dilepton.Pt()/dilepton.M())

            # jets
            ntuple['bjet_pt'].append(bjet.Pt())
            ntuple['bjet_eta'].append(bjet.Eta())
            ntuple['bjet_phi'].append(bjet.Phi())
            ntuple['bjet_e'].append(bjet.E())
            ntuple['bjet_d0'].append(tree.bjetD0)
            ntuple['bjet_tag'].append(tree.bjetTag)
            ntuple['bjet_flavor'].append(tree.bjetFlavor)

            ntuple['jet_pt'].append(jet.Pt())
            ntuple['jet_eta'].append(jet.Eta())
            ntuple['jet_phi'].append(jet.Phi())
            ntuple['jet_e'].append(jet.E())
            ntuple['jet_d0'].append(tree.jetD0)
            ntuple['jet_tag'].append(tree.jetTag)
            ntuple['jet_flavor'].append(tree.jetFlavor)

            ntuple['jet_delta_eta'].append(abs(bjet.Eta() - jet.Eta()))
            ntuple['jet_delta_phi'].append(abs(bjet.DeltaPhi(jet)))
            ntuple['jet_delta_r'].append(bjet.DeltaR(jet))

            ntuple['dijet_mass'].append(dijet.M())
            ntuple['dijet_pt'].append(dijet.Pt())
            ntuple['dijet_eta'].append(dijet.Eta())
            ntuple['dijet_phi'].append(dijet.Phi())
            ntuple['dijet_pt_over_m'].append(dijet.Pt()/dijet.M() if dijet.M() > 0 else -1)

            # MET
            ntuple['met_mag'].append(met)
            ntuple['met_phi'].append(met_phi)

            # lepton + b jet variables
            ntuple['lepton1_b_mass'].append(lepton1_b.M())
            ntuple['lepton1_b_pt'].append(lepton1_b.Pt())
            ntuple['lepton1_b_delta_eta'].append(abs(lep1.Eta() - bjet.Eta()))
            ntuple['lepton1_b_delta_phi'].append(abs(lep1.DeltaPhi(bjet)))
            ntuple['lepton1_b_delta_r'].append(lep1.DeltaR(bjet))

            ntuple['lepton2_b_mass'].append(lepton2_b.M())
            ntuple['lepton2_b_pt'].append(lepton2_b.Pt())
            ntuple['lepton2_b_delta_eta'].append(abs(lep2.Eta() - bjet.Eta()))
            ntuple['lepton2_b_delta_phi'].append(abs(lep2.DeltaPhi(bjet)))
            ntuple['lepton2_b_delta_r'].append(lep2.DeltaR(bjet))

            # dilepton + b jet variables
            ntuple['dilepton_b_mass'].append(dilepton_b.M())
            ntuple['dilepton_b_pt'].append(dilepton_b.Pt())
            ntuple['dilepton_b_delta_eta'].append(abs(dilepton.Eta() - bjet.Eta()))
            ntuple['dilepton_b_delta_phi'].append(abs(dilepton.DeltaPhi(bjet)))
            ntuple['dilepton_b_delta_r'].append(dilepton.DeltaR(bjet))

            # dilepton + b jet variables
            ntuple['dilepton_j_mass'].append(dilepton_j.M())
            ntuple['dilepton_j_pt'].append(dilepton_j.Pt())
            ntuple['dilepton_j_delta_eta'].append(abs(dilepton.Eta() - jet.Eta()))
            ntuple['dilepton_j_delta_phi'].append(abs(dilepton.DeltaPhi(jet)))
            ntuple['dilepton_j_delta_r'].append(dilepton.DeltaR(jet))

            # four body variables
            ntuple['four_body_mass'].append(fourbody.M())
            ntuple['four_body_delta_eta'].append(abs(dijet.Eta() - dilepton.Eta()))
            ntuple['four_body_delta_phi'].append(abs(dijet.DeltaPhi(dilepton)))
            ntuple['four_body_delta_r'].append(dijet.DeltaR(dilepton))

            # Mendelstam variables
            #ntuple['t_xj'].append((dilepton - jet).M()**2)
            #ntuple['t_xb'].append((dilepton - bjet).M()**2)
            #ntuple['t_bj'].append((bjet - jet).M()**2)

            # generator level information (all 0 for data)
            #ntuple['n_partons'].append(tree.nPartons)
            #ntuple['gen_bjet_pt'].append(gen_bjet.Pt())
            #ntuple['gen_bjet_eta'].append(gen_bjet.Eta())
            #ntuple['gen_bjet_phi'].append(gen_bjet.Phi())
            #ntuple['gen_bjet_e'].append(gen_bjet.E())
            #ntuple['gen_bjet_tag'].append(tree.genBJetTag)
            #ntuple['gen_dilepton_b_mass'].append((dilepton + gen_bjet).M())

            #ntuple['gen_jet_pt'].append(gen_jet.Pt())
            #ntuple['gen_jet_eta'].append(gen_jet.Eta())
            #ntuple['gen_jet_phi'].append(gen_jet.Phi())
            #ntuple['gen_jet_e'].append(gen_bjet.E())
            #ntuple['gen_jet_tag'].append(tree.genJetTag)
            #ntuple['gen_dilepton_j_mass'].append((dilepton + gen_jet).M())

        df = pd.DataFrame(ntuple)
        df.to_pickle('{0}/ntuple_{1}.pkl'.format(output_path, dataset))

    fname = '{0}/event_counts.csv'.format(output_path)
    df = pd.DataFrame(event_count)
    df.to_csv(fname)

    #if not os.path.isfile(fname):
    #    df = pd.DataFrame(event_count)
    #    df.to_csv(fname)

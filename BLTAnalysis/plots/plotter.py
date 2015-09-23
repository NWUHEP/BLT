from ROOT import gROOT, gStyle, gPad, TFile, TTree, TH1F, TCanvas
import sys
import re


# ______________________________________________________________________________
donotdelete = []  # persist in memory

def make_plots(ttree, verbose):
    h_dimuon_mass = TH1F("dimuon_mass", "; M(#mu#mu) [GeV]", 200, 70, 120);
    h_dimuon_mass_gen = TH1F("dimuon_mass_gen", "; M(#mu#mu) [GeV]", 200, 70, 120);

    for ievt, evt in enumerate(ttree):
        h_dimuon_mass.Fill(evt.dimuon.M());
        h_dimuon_mass_gen.Fill(evt.genZ.M());

    h_dimuon_mass.SetLineWidth(2)
    h_dimuon_mass.SetLineColor(1)
    h_dimuon_mass_gen.SetLineWidth(2)
    h_dimuon_mass_gen.SetLineColor(2)

    h_dimuon_mass.SetMaximum(h_dimuon_mass.GetMaximum() * 1.6)
    h_dimuon_mass.Draw()
    h_dimuon_mass_gen.Draw("same")

    #print h_dimuon_mass.Integral(), h_dimuon_mass_gen.Integral()

    donotdelete.append(h_dimuon_mass)
    donotdelete.append(h_dimuon_mass_gen)
    

# ______________________________________________________________________________
if __name__=="__main__":

    # Init
    gROOT.LoadMacro("tdrstyle.C")
    gROOT.ProcessLine("setTDRStyle()")

    suffix = "DYJetsToLL"
    dataname = "DYJetsToLL_test"
    jobcount = "local"

    fname = "demoFile_%s_%s.root" % (dataname, jobcount)
    tname = "demoTree_%s" % suffix

    tfile = TFile.Open(fname)
    ttree = tfile.Get(tname)
    make_plots(ttree, verbose=1)


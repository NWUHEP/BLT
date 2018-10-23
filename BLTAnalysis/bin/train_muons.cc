#include <vector>

#include <TCut.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooConstVar.h>
#include <RooWorkspace.h>

// GBRLikelihood
#include "HiggsAnalysis/GBRLikelihood/interface/RooRevCBExp.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooGausDoubleExp.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForest.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

using namespace RooFit;

int main(int argc, char **argv) { // Trains semi-parametric MVA for muons.
    const char* infile = "output_dilepton_0.root";
    const char* outfile = "muon_training.root";
   
    RooArgList allvars; // input variables + target variable
    allvars.addOwned(*new RooRealVar("var01", "rho", 0));
    allvars.addOwned(*new RooRealVar("var02", "leptonOneEnergy", 0));
    allvars.addOwned(*new RooRealVar("var03", "leptonOneEta", 0));
    allvars.addOwned(*new RooRealVar("var04", "leptonOneTkNChi2", 0));
    allvars.addOwned(*new RooRealVar("var05", "leptonOneNTkLayers", 0));
    allvars.addOwned(*new RooRealVar("var06", "leptonOneNPixHits", 0));
    allvars.addOwned(*new RooRealVar("var07", "leptonOneNValidHits", 0));
    allvars.addOwned(*new RooRealVar("var08", "leptonOneNMatchStn", 0));
    allvars.addOwned(*new RooRealVar("var09", "leptonOnePFIsoCH", 0));
    allvars.addOwned(*new RooRealVar("var10", "leptonOnePFIsoNH", 0));
    allvars.addOwned(*new RooRealVar("var11", "leptonOnePFIsoPho", 0));
    allvars.addOwned(*new RooRealVar("var12", "leptonOnePFIsoPU", 0));
   
    RooArgList invars = allvars; // input variables only
    RooRealVar* target = new RooRealVar("target", "leptonOneERes", 1., 0.6, 1.5); // target variable
    allvars.addOwned(*target);
   
    // variables corresponding to regressed parameters
    RooRealVar mean("mean", "", 1.0);
    RooRealVar sigma("sigma", "", 0.02);
    RooRealVar alphaL("alphaL", "", 1.5);
    RooRealVar alphaR("alphaR", "", 1.5);
    mean.setConstant(false);
    sigma.setConstant(false);
    alphaL.setConstant(false);
    alphaR.setConstant(false);
   
    // non-parametric functions for each regressed parameter
    RooGBRFunctionFlex funcMean("funcMean", "");
    RooGBRFunctionFlex funcSigma("funcSigma", "");
    RooGBRFunctionFlex funcAlphaL("funcAlphaL", "");
    RooGBRFunctionFlex funcAlphaR("funcAlphaR", "");
   
    // mapping of input variables to non-parametric functions
    RooGBRTargetFlex tgtMean("tgtMean", "", funcMean, mean, invars);
    RooGBRTargetFlex tgtSigma("tgtSigma", "", funcSigma, sigma, invars);
    RooGBRTargetFlex tgtAlphaL("tgtAlphaL", "", funcAlphaL, alphaL, invars);
    RooGBRTargetFlex tgtAlphaR("tgtAlphaR", "", funcAlphaR, alphaR, invars);
   
    // parameters' bounds
    RooRealConstraint limMean("limMean", "", tgtMean, 0.6, 1.5);
    RooRealConstraint limSigma("limSigma", "", tgtSigma, 0.001, 0.4);
    RooRealConstraint limAlphaL("limAlphaL", "", tgtAlphaL, 0.2, 5.);
    RooRealConstraint limAlphaR("limAlphaR", "", tgtAlphaR, 0.2, 5.);

    // Gaussian + left exponential tail + right exponential tail
    RooAbsPdf* pdf = new RooGausDoubleExp("pdfGausDoubleExp", "", *target, limMean, limSigma, limAlphaL, limAlphaR);

    // list of mapped functions to regress
    RooArgList tgts;
    tgts.add(tgtMean);
    tgts.add(tgtSigma);
    tgts.add(tgtAlphaL);
    tgts.add(tgtAlphaR);

    // list of pdfs
    std::vector<RooAbsReal*> pdfs;
    pdfs.push_back(pdf);

    // open file and get tree with the inputs and the target
    TFile* fi = TFile::Open(infile);
    if (!fi || fi->IsZombie())
        FATAL("TFile::Open() failed");
    
    TTree* tree = dynamic_cast<TTree*>(fi->Get("mm/bltTree_dilepton"));
    //TTree* tree = dynamic_cast<TTree*>(fi->Get("mm/bltTree_zjets_m-50_amc"));
    if (!tree) FATAL("TFile::Get() failed");
   
    // create a memory-resident friend TTree with linear event numbers
    if (!gROOT->cd()) FATAL("TROOT::cd() failed");
    TTree evtree("FlatTree", "Trivial event numbers");
    evtree.SetAutoFlush(0);
    evtree.SetAutoSave(0);
    Long64_t event;
    evtree.Branch("event", &event);
    for (event = 0; event < tree->GetEntriesFast(); event++)
      evtree.Fill();
    tree->AddFriend(&evtree);
 
    //TCut cuts = "";
    TCut cuts = "event % 1 == 0";     // NOTE: take all tree entries
    // per-event weight
    RooRealVar weightvar("weightvar", "", 1.);
    weightvar.SetTitle(cuts);
   
    // list of training datasets
    RooDataSet* dataset = RooTreeConvert::CreateDataSet("data", tree, allvars, weightvar);
    std::vector<RooAbsData*> datasets;
    datasets.push_back(dataset);

    // minimum event weight per tree
    std::vector<double> minweights;
    minweights.push_back(200);
   
    // dummies
    RooConstVar etermconst("etermconst", "", 0.);
    RooRealVar r("r", "", 1.);
    r.setConstant(true);
   
    // training
    RooHybridBDTAutoPdf bdtpdfdiff("bdtpdfdiff", "", tgts, etermconst, r, datasets, pdfs);
    bdtpdfdiff.SetMinCutSignificance(5.);

    //bdtpdfdiff.SetPrescaleInit(100);
    bdtpdfdiff.SetShrinkage(0.1);
    bdtpdfdiff.SetMinWeights(minweights);
    bdtpdfdiff.SetMaxNodes(750);
    bdtpdfdiff.TrainForest(1e+6); // NOTE: valid training will stop at ~100-500 trees
   
    // save output to file
    RooWorkspace* ws = new RooWorkspace("ws_mva_muons");
    ws->import(*pdf);
    ws->writeToFile(outfile, false); // false = update output file, not recreate
   
    // NOTE: no memory cleanup for simplicity

    return EXIT_SUCCESS;
}

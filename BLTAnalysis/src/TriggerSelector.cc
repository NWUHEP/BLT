#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"

#include <cstdlib>  // for getenv

//FIXME: implement how to get trigger objects

TriggerSelector::TriggerSelector() {
    std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/";
    trigfilename += "BaconAna/DataFormats/data/HLTFile_50ns";
    _ttrigger.reset(new baconhep::TTrigger(trigfilename));
}

bool TriggerSelector::PassTrigger(const TriggerBits &iTrig) const {
    bool pass = false;
    for (const auto& triggerName: _triggerNames) {
        pass |= _ttrigger->pass(triggerName, iTrig);
    }
    return pass;
}

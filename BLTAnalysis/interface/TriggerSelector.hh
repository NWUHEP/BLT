#ifndef TRIGGERSELECTOR_HH
#define TRIGGERSELECTOR_HH


// Bacon header files
#include "BaconAna/Utils/interface/TTrigger.hh"

#include <string>
#include <vector>
#include <memory>


class TriggerSelector {
public:

    TriggerSelector();
    ~TriggerSelector() {}

    // Setters
    void SetRealData(bool isRealData)       { _isRealData = isRealData; }
    void SetTriggerNames(const std::vector<std::string>& triggerNames) { _triggerNames = triggerNames; }

    // Get trigger decisions
    // N.B. why isn't TriggerBits under the namespace baconhep?
    bool PassTrigger(const TriggerBits &iTrig) const;

private:
    std::unique_ptr<baconhep::TTrigger>  _ttrigger;
    std::vector<std::string>             _triggerNames;
    bool                                 _isRealData;
};

#endif  // TRIGGERSELECTOR_HH

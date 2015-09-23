#ifndef TRIGGERSELECTOR_HH
#define TRIGGERSELECTOR_HH

#include <string>
#include <vector>


class TriggerSelector {
public:
    typedef std::vector<std::string>  vstring;

    TriggerSelector();
    ~TriggerSelector() {}

    // Setters
    void SetTriggers(const std::vector<std::string>& triggers) { _triggers = triggers; }

    // Get trigger decisions
    bool PassTrigger() const;

private:
    vstring         _triggers;
    bool            _isRealData;
};

#endif  // TRIGGERSELECTOR_HH

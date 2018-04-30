#ifndef ELECTRONCORRECTOR_H
#define ELECTRONCORRECTOR_H

#include <map>
#include <cmath>
#include <string>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TMath.h"

#include "BaconAna/DataFormats/interface/TElectron.hh"

using namespace std;
using namespace baconhep;

/*
 * r9 bins
 * 0: r9 < 0.94
 * 1: r9 >= 0.94
 *
 * eta bins
 * 0: 0 < eta <= 1
 * 1: 1 < eta <= 1.479
 * 2: 1.479 < eta <= 2
 * 3: 2 < eta <= 2.5
 * 
 */

//============================== 

struct scaleData
{
    float scale, scaleErr;
};

struct smearData
{
    float scale, scaleErr, rho, rhoErr;
};

// [r9, eta]: smear_data
typedef map<int, map<int, smearData> > smearMap; 
// [r9, eta, run number]: scale_data
typedef map<int, map<int, map<int, scaleData> > > scaleMap; 

//============================== Main class
class EnergyScaleCorrection
{
    public:
        EnergyScaleCorrection(std::string filePath);
        EnergyScaleCorrection(){}; 
        ~EnergyScaleCorrection(void){};

        scaleData GetScaleData(TElectron* electron, int runNumber);
        smearData GetSmearData(TElectron* electron);
        float GetSmearingFactor(TElectron* electron, int nSigRho, int nSigPhi);
    private:
        vector<int> _runNumbers;
        smearMap _smearMap;
        scaleMap _scaleMap;
};

#endif

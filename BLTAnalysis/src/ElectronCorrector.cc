#include "BLT/BLTAnalysis/interface/ElectronCorrector.h"

#include <cstdlib>
#include <cfloat> 
#include <iomanip>
#include <sstream>

using namespace std;

EnergyScaleCorrection::EnergyScaleCorrection(std::string filePath)
{
    // open scale file and unpack values
    int r9Bin, etaBin, runMin, runMax;
    float scale, scaleErr, nada;

    string line;
    std::ifstream scaleFile("/tthome/naodell/work/CMSSW_7_4_12/src/BLT/BLTAnalysis/data/electron_scale.dat");
    while (getline(scaleFile, line)) {
        istringstream iss(line);
        iss >> r9Bin >> etaBin 
            >> runMin >> runMax
            >> scale >> scaleErr
            >> nada >> nada >> nada;
        runBin = 0;
        for (const auto& iRun: runRanges) {
            if (runMin == iRun) 
                break;
            ++runBin
        }

        scaleData data = {scale, scaleErr};
        this->_scaleMap[r9Bin][etaBin][runBin] = data;
    }
    scaleFile.close();

    // open smear file and unpack values
    float rho, rhoErr;
    std::ifstream smearFile("/tthome/naodell/work/CMSSW_7_4_12/src/BLT/BLTAnalysis/data/electron_smear.dat");
    while (getline(smearFile, line)) {
        istringstream iss(line);
        iss >> r9Bin >> etaBin 
            >> scale >> scaleErr
            >> rho >> rhoErr;

        smearData data = {scale, scaleErr, rho, rhoErr};
        this->_smearMap[r9Bin][etaBin] = data;
    }
    smearFile.close();
}

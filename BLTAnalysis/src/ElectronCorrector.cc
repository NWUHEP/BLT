#include "BLT/BLTAnalysis/interface/ElectronCorrector.h"

EnergyScaleCorrection::EnergyScaleCorrection(std::string filePath)
{
    // open scale file and unpack values
    int r9Bin, etaBin, runMin, runMax;
    float scale, scaleErr, nada;

    string line;
    std::ifstream scaleFile(filePath + "/electron_scale.dat");
    while (getline(scaleFile, line)) {
        istringstream iss(line);
        iss >> r9Bin >> etaBin 
            >> runMin >> runMax
            >> scale >> scaleErr
            >> nada >> nada >> nada;

        //cout << " " << r9Bin << " " << etaBin 
        //    << " " << runMin << " " << runMax
        //    << " " << scale << " " << scaleErr
        //    << endl;

        if (r9Bin == 0 and etaBin == 0) 
            _runNumbers.push_back(runMin);
        sort(_runNumbers.begin(), _runNumbers.end());

        scaleData data = {scale, scaleErr};
        this->_scaleMap[r9Bin][etaBin][runMin] = data;
    }
    scaleFile.close();

    // open smear file and unpack values
    float rho, rhoErr;
    std::ifstream smearFile(filePath + "/electron_smear.dat");
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

scaleData EnergyScaleCorrection::GetScaleData(TElectron* electron, int runNumber)
{
    unsigned r9Bin = 0;
    if (electron->r9 >= 0.94) {
        r9Bin = 1;
    }

    unsigned etaBin = 0;
    float etaBinning[] = {0, 1, 1.479, 2, 2.5};
    for (unsigned i = 0; i < 4; ++i) {
        if (fabs(electron->scEta) >= etaBinning[i] && fabs(electron->scEta) < etaBinning[i+1]) {
            etaBin = i;
            break;
        }
    }

    int runMin = 0;
    for (unsigned i = 0; i < _runNumbers.size(); ++i) {
        if (runNumber >= _runNumbers[i] && runNumber < _runNumbers[i+1]) {
            runMin = _runNumbers[i];
            break;
        }
    }

    if (runMin == 0) runMin = _runNumbers[0];

    return _scaleMap[r9Bin][etaBin][runMin];
}

smearData EnergyScaleCorrection::GetSmearData(TElectron* electron)
{
    unsigned r9Bin = 0;
    if (electron->r9 >= 0.94) {
        r9Bin = 1;
    }

    unsigned etaBin = 0;
    float etaBinning[] = {0, 1, 1.479, 2, 2.5};
    for (unsigned i = 0; i < 4; ++i) {
        if (fabs(electron->scEta) >= etaBinning[i] && fabs(electron->scEta) < etaBinning[i+1]) {
            etaBin = i;
            break;
        }
    }

    return _smearMap[r9Bin][etaBin];
}

float EnergyScaleCorrection::GetSmearingFactor(TElectron* electron, int nSigRho, int nSigPhi)
{
    smearData sdata = GetSmearData(electron);
    float rho = sdata.rho + sdata.rhoErr*nSigRho;
    float phi = TMath::Pi()/2 + nSigPhi*TMath::Pi()/2;

    float c = rho*sin(phi);
    float alpha = rho*sdata.scale*cos(phi);

    return sqrt(c*c + alpha*alpha/electron->pt);
}

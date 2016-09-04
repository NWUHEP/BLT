#include <iostream>
#include <TChain.h>
#include <TClonesArray.h>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "muresolution.h"

class rochcor2012 {
 public:
  rochcor2012();
  rochcor2012(int seed);
  ~rochcor2012();
  
  void momcor_mc(TLorentzVector&, float, int, float&);
  void momcor_data(TLorentzVector&, float, int, float&);
  
  int aetabin(double);
  int etabin(double);
  int phibin(double);
  
 private:
  
  TRandom3 eran;
  TRandom3 sran;
  
  
  //  static float netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};
  const double pi = TMath::Pi();
  static const double netabin[25];
  static const double anetabin[13];
  
  const double genm_smr = 9.10193e+01; //gen mass peak with eta dependent gaussian smearing => better match in Z mass profile vs. eta/phi
  const double genm = 91.06; //gen mass peak without smearing => Z mass profile vs. eta/phi in CMS note
  
  const double mrecm = 9.10411e+01; //rec mass peak in MC 
  const double drecm = 9.06951e+01; //rec mass peak in data 
  const double mgscl_stat = 0.0001; //stat. error of global factor for mass peak in MC 
  const double mgscl_syst = 0.0006; //syst. error of global factor for mass peak in MC  
  const double dgscl_stat = 0.0001; //stat. error of global factor for mass peak in data 
  const double dgscl_syst = 0.0008; //syst. error of global factor for mass peak in data 
        
  static const double sf[3];
  static const double sfer[3];

  //---------------------------------------------------------------------------------------------
  
  static const double dcor_bf[16][24];  
  static const double dcor_ma[16][24];
  static const double mcor_bf[16][24];
  static const double mcor_ma[16][24];
  static const double dcor_bfer[16][24];  
  static const double dcor_maer[16][24];
  static const double mcor_bfer[16][24];
  static const double mcor_maer[16][24];
  
  //=======================================================================================================
  
  static const double dmavg[16][24];  
  static const double dpavg[16][24];  
  static const double mmavg[16][24];  
  static const double mpavg[16][24];

  static const double dd[12];
  static const double dder[12];
  static const double md[12];
  static const double mder[12];  

  static const double dscl[24];
  static const double mscl[24];
  //===============================================================================================

  double mptsys_mc_dm[16][24];
  double mptsys_mc_da[16][24];
  double mptsys_da_dm[16][24];
  double mptsys_da_da[16][24];

  double gscler_mc_dev;
  double gscler_da_dev;


};
  

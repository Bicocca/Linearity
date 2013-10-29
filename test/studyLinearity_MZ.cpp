#include "setTDRStyle.h"
#include "ntpleUtils.h"
#include "geometryUtils.h"
#include "PUReweighting.h"
#include "GetScaleCorrection.h"
#include "GetExtraSmearing.h"
#include "GetCategories.h"
#include "ScaleEstimators.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "TSystem.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TArrow.h"
#include "TRandom3.h"



//-----------------
// global variables

std::vector<double>* EtBinEdges;
std::vector<double>* HtBinEdges;
double* EtBinEdgesDouble;
double* HtBinEdgesDouble;
unsigned int nEtBins;
unsigned int nHtBins;

TH1F** h_Ht_HtBin_MC;

TH1F** h_Et_EtBin_MC;
TH1F** h_Et_EtBin_DA;
TH1F** h_Et_EtBin_fit_DA;
TH1F** h_Et_EtBin_gausFit_DA;
TH1F** h_Et_EtBin_mean_DA;
TH1F** h_Et_EtBin_recursiveMean_DA;
TH1F** h_Et_EtBin_smallestInterval_DA;

TFile* outFile;
std::string extension = "pdf";

TF1* f_scaleVsEt;
TF1* f_scaleVs2Et;
TF1* f_scaleVsAvgEt;
TF1* f_invScaleVsEt;

double invScaleVsEt(double* x, double* par);
double scaleVsAvgEt(double* x, double* par);






int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if( argc != 2 )
  {
    std::cerr << ">>> studyLinearity_MZ::usage: " << argv[0] << " configFileName" << std::endl;
    return -1;
  }
  
  
  
  //----------------------
  // Parse the config file
  
  parseConfigFile(argv[1]);
  
  std::string inputFilesDA = gConfigParser -> readStringOption("Input::inputFilesDA");
  std::string inputFilesMC = gConfigParser -> readStringOption("Input::inputFilesMC");
  
  std::string year       = gConfigParser -> readStringOption("Options::year");
  std::string dataLabel  = gConfigParser -> readStringOption("Options::dataLabel");

  bool useGlobeNtuple   = gConfigParser -> readBoolOption("Options::useGlobeNtuple");
  bool useShervinNtuple = gConfigParser -> readBoolOption("Options::useShervinNtuple");
  
  int nBinsMee  = gConfigParser -> readIntOption("Options::nBinsMee");
  double meeMin = gConfigParser -> readDoubleOption("Options::meeMin");
  double meeMax = gConfigParser -> readDoubleOption("Options::meeMax");
  
  int maxEntries    = gConfigParser -> readIntOption("Options::maxEntries");
  std::string MCGen = gConfigParser -> readStringOption("Options::MCGen");
  bool runDepFlag   = gConfigParser -> readBoolOption("Options::runDepFlag");
  int runMin        = gConfigParser -> readIntOption("Options::runMin");
  int runMax        = gConfigParser -> readIntOption("Options::runMax");
  
  std::string eleIDSelection = gConfigParser -> readStringOption("Options::eleIDSelection");
  int eleIDBit = 1;
  if( eleIDSelection == "loose"  ) eleIDBit = 2;
  if( eleIDSelection == "medium" ) eleIDBit = 6;
  if( eleIDSelection == "tight"  ) eleIDBit = 14;
  if( eleIDSelection == "WP90PU" ) eleIDBit = 16;
  if( eleIDSelection == "WP80PU" ) eleIDBit = 48;
  
  bool applyPUWeight        = gConfigParser -> readBoolOption("Options::applyPUWeight");
  bool applyEnergyScaleCorr = gConfigParser -> readBoolOption("Options::applyEnergyScaleCorr");
  bool applyEnergySmearing  = gConfigParser -> readBoolOption("Options::applyEnergySmearing");
  bool applyEtaR9Reweighting = gConfigParser -> readBoolOption("Options::applyEtaR9Reweighting");
  
  std::string enCorrType          = gConfigParser -> readStringOption("Options::enCorrType");
  std::string energyScaleCorrType = gConfigParser -> readStringOption("Options::energyScaleCorrType");
  std::string energySmearingType  = gConfigParser -> readStringOption("Options::energySmearingType");
  
  std::string runRangeFile        = gConfigParser -> readStringOption("Options::runRangeFile");
  std::string ShervinScaleFile    = gConfigParser -> readStringOption("Options::ShervinScaleFile");
  std::string ShervinSmearingFile = gConfigParser -> readStringOption("Options::ShervinSmearingFile");
  std::string IJazZGlobalFolder   = gConfigParser -> readStringOption("Options::IJazZGlobalFolder");
  std::string IJazZRunDepFolder   = gConfigParser -> readStringOption("Options::IJazZRunDepFolder");
  
  std::string catType = gConfigParser -> readStringOption("Options::catType");
  int category        = gConfigParser -> readIntOption("Options::category");
  int evtsPerPoint    = gConfigParser -> readIntOption("Options::evtsPerPoint");

  std::vector<float> extHtBinEdges;
  if(category == 0) extHtBinEdges = gConfigParser -> readFloatListOption("Options::extHtBinEdges0");
  if(category == 1) extHtBinEdges = gConfigParser -> readFloatListOption("Options::extHtBinEdges1");
  if(category == 2) extHtBinEdges = gConfigParser -> readFloatListOption("Options::extHtBinEdges2");
  if(category == 3) extHtBinEdges = gConfigParser -> readFloatListOption("Options::extHtBinEdges3");
  


  bool MCClosure = gConfigParser -> readBoolOption("Options::MCClosure");
  bool MCHiggs   = gConfigParser -> readBoolOption("Options::MCHiggs");
  
  double DphiMax = gConfigParser -> readDoubleOption("Options::DphiMax");
  
  std::string outFilePath = gConfigParser -> readStringOption("Output::outFilePath");
  
  int nSteps = 1;
  float smearingSyst = 0.;
  
  
  
  
  //------------------
  // Set style options
  
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  
  
  //------------------
  // Fitting functions
  
//    f_scaleVsEt  = new TF1("f_scaleVsEt", "1.+0.000",0., 1000.);
//    f_scaleVs2Et = new TF1("f_scaleVsEt", "1.+0.000",0., 1000.);
  
  //f_scaleVsEt  = new TF1("f_scaleVsEt", "1.+0.002",0., 1000.);
  //f_scaleVs2Et = new TF1("f_scaleVsEt", "1.+0.002",0., 1000.);
 

  f_scaleVsEt = new TF1("f_scaleVsEt", "1. + [0] * (1 - exp(-[1] * (x-45.)) )",0., 1000.);
  f_scaleVsEt -> SetParameters(7.50e-03,2.00e-02);
  
  f_scaleVs2Et = new TF1("f_scaleVs2Et", "1. + [0] * (1 - exp(-[1] * (0.5*x-45.)) )",0., 1000.);
  f_scaleVs2Et -> SetParameters(7.50e-03,2.00e-02);
  
  f_scaleVsAvgEt = new TF1("f_scaleVsAvgEt",scaleVsAvgEt,0.,1000.,0);
  f_invScaleVsEt = new TF1("f_invScaleVsEt",invScaleVsEt,0.,1000.,0);
  
  ///Error funcions by Stat. toy
  TF1* graphError_SmallestInterval = new TF1("graphError_SmallestInterval","[0] + [1]/sqrt(x)", 1., 10000000);
  if(category == 0) graphError_SmallestInterval->SetParameters(0.000592, 0.0602);
  if(category == 1) graphError_SmallestInterval->SetParameters(0.000958, 0.0760);
  if(category == 2) graphError_SmallestInterval->SetParameters(0.00120, 0.0861);
  if(category == 3) graphError_SmallestInterval->SetParameters(0.00132, 0.105);

  TF1* graphError_RecursiveMean = new TF1("graphError_RecursiveMean","[0] + [1]/sqrt(x)", 1., 10000000);
  if(category == 0) graphError_RecursiveMean->SetParameters(0.0000221, 0.0413);
  if(category == 1) graphError_RecursiveMean->SetParameters(0.0000224, 0.0494);
  if(category == 2) graphError_RecursiveMean->SetParameters(0.0000657, 0.0528);
  if(category == 3) graphError_RecursiveMean->SetParameters(0.00000772, 0.0618);
  
//   ///////////////////////////////////////////////////////////////////////////////////
//   //Inject scale from stdCat in MC and look at DATA in Cic 

// TF1_pol1
// TF1_pol1_std_cat0.root
// TF1_pol1_std_cat1.root
// TF1_pol1_std_cat2.root
// TF1_pol1_std_cat3.root

  
  //----------
  // Get trees
  std::cout << std::endl;
  std::cout << ">>> Get trees" << std::endl;
  
  std::string treeNameMC;
  std::string treeNameDA;
  
  if( !useShervinNtuple && !useGlobeNtuple )
  {
    treeNameMC = "simpleNtupleEoverP/SimpleNtupleEoverP";
    treeNameDA = "simpleNtupleEoverP/SimpleNtupleEoverP";
  }
  if( useShervinNtuple )
  {
    treeNameMC = "selected";
    treeNameDA = "selected";
    if( MCClosure == 1 )
    {
      inputFilesDA = inputFilesMC;
    }
  }
  if( useGlobeNtuple )
  {
    if( MCClosure == 0 && MCHiggs == 0 )
    {
      treeNameMC = "DYJetsToLL";
      treeNameDA = "Data";
    }
    if( MCClosure == 1 && MCHiggs == 0 )
    {
      treeNameMC = "DYJetsToLL";
      treeNameDA = "DYJetsToLL";
    }
    if( MCHiggs == 1 )
    {
      treeNameMC = "ggh_m125_8TeV";
      treeNameDA = "ggh_m125_8TeV";
    }
  }
  
  TChain* ntu_MC = new TChain(treeNameMC.c_str());
  FillChain(ntu_MC,inputFilesMC);
  std::cout << ">>>   MC: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;
  
  TChain* ntu_DA = new TChain(treeNameDA.c_str());
  FillChain(ntu_DA,inputFilesDA);
  std::cout << ">>> DATA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;
  
  if( ntu_MC->GetEntries() == 0 || ntu_DA->GetEntries() == 0 )
  {
    std::cout << ">>> compareZPeaks::Error: at least one file is empty" << std::endl;
    return -1;
  }
  
  
  
  //------------------------
  // Define branch addresses
  std::cout << std::endl;
  std::cout << ">>> Define branch addresses" << std::endl;
  
  // vectors
  std::vector<double> scEta1_MC, scEta2_MC;
  std::vector<double> scEta1_DA, scEta2_DA;
  std::vector<double> weight_MC, Ht_MC;
  std::vector<double> weight_DA, Ht_DA;
  std::vector<double> Zpt_MC, mee_MC, Et1_MC, Et2_MC, R91_MC, R92_MC;
  std::vector<double> Zpt_DA, mee_DA, Et1_DA, Et2_DA, R91_DA, R92_DA;
  std::vector<double> mee_fit_DA, Et1_fit_DA, Et2_fit_DA;
  std::vector<double> mee_gausFit_DA, Et1_gausFit_DA, Et2_gausFit_DA;
  std::vector<double> mee_mean_DA, Et1_mean_DA, Et2_mean_DA;
  std::vector<double> mee_recursiveMean_DA, Et1_recursiveMean_DA, Et2_recursiveMean_DA;
  std::vector<double> mee_smallestInterval_DA, Et1_smallestInterval_DA, Et2_smallestInterval_DA;
  
  // global variables
  int runId,nPU,cat;
  float weight;
  bool HLTfire;
  
  if( useGlobeNtuple )
  {
    HLTfire = true;
    
    ntu_MC -> SetBranchStatus("*",0);
    ntu_MC -> SetBranchStatus("run",   1); ntu_MC -> SetBranchAddress("run",   &runId);
    ntu_MC -> SetBranchStatus("weight",1); ntu_MC -> SetBranchAddress("weight",&weight);
    if( std::string(catType) == "CiC") {
      ntu_MC -> SetBranchStatus("category_baseline",1); ntu_MC -> SetBranchAddress("category_baseline",&cat); }
    if( std::string(catType) == "MVA") {
      ntu_MC -> SetBranchStatus("category",1); ntu_MC -> SetBranchAddress("category",&cat); }
    
    ntu_DA -> SetBranchStatus("*",0);                         
    ntu_DA -> SetBranchStatus("run",   1); ntu_DA -> SetBranchAddress("run",   &runId);
    ntu_DA -> SetBranchStatus("weight",1); ntu_DA -> SetBranchAddress("weight",&weight);
    if( std::string(catType) == "CiC") {
      ntu_DA -> SetBranchStatus("category_baseline",1); ntu_DA -> SetBranchAddress("category_baseline",&cat); }
    if( std::string(catType) == "MVA") {
      ntu_DA -> SetBranchStatus("category",1); ntu_DA -> SetBranchAddress("category",&cat); }
  }
  if( useShervinNtuple )
  {
    ntu_MC -> SetBranchStatus("*",0);
    ntu_MC -> SetBranchStatus("HLTfire",  1);   ntu_MC -> SetBranchAddress("HLTfire",&HLTfire);
    ntu_MC -> SetBranchStatus("runNumber",1);   ntu_MC -> SetBranchAddress("runNumber",&runId);
    ntu_MC -> SetBranchStatus("nPU",      1);   ntu_MC -> SetBranchAddress("nPU",&nPU);
    
    ntu_DA -> SetBranchStatus("*",0);
    ntu_DA -> SetBranchStatus("HLTfire",  1);   ntu_DA -> SetBranchAddress("HLTfire",&HLTfire);
    ntu_DA -> SetBranchStatus("runNumber",1);   ntu_DA -> SetBranchAddress("runNumber",&runId);
  }
  if( !useGlobeNtuple && !useShervinNtuple )
  {
    HLTfire = true;
    
    ntu_MC -> SetBranchStatus("*",0);
    ntu_MC -> SetBranchStatus("runId",1); ntu_MC -> SetBranchAddress("runId",&runId);
    ntu_MC -> SetBranchStatus("PUit_TrueNumInteractions",1); ntu_MC -> SetBranchAddress("PUit_TrueNumInteractions",&nPU);
    
    ntu_DA -> SetBranchStatus("*",0);
    ntu_DA -> SetBranchStatus("runId",1); ntu_DA -> SetBranchAddress("runId",&runId);
  }
  
  // electron variables
  float scEta[2];
  float scPhi[2];
  float eta[2];
  float phi[2];
  float scE[2];
  float scERaw[2];
  float scEReg[2];
  float R9[2];
  int eleID[2];
  
  double eta1, phi1;
  double eta2, phi2;
  float scERaw1, scEReg1, scEta1, scPhi1, etaFloat1, phiFloat1, E3x31, R91;
  float scERaw2, scEReg2, scEta2, scPhi2, etaFloat2, phiFloat2, E3x32, R92;
  int eleID1, eleID2;
    
  if( useGlobeNtuple )
  {  
    ntu_DA->SetBranchStatus("pho1_energy_regr",1); ntu_DA->SetBranchAddress("pho1_energy_regr",&scEReg1);
    ntu_DA->SetBranchStatus("pho1_sceta",      1); ntu_DA->SetBranchAddress("pho1_sceta",      &scEta1);
    ntu_DA->SetBranchStatus("pho1_scphi",      1); ntu_DA->SetBranchAddress("pho1_scphi",      &scPhi1);
    ntu_DA->SetBranchStatus("pho1_eta",        1); ntu_DA->SetBranchAddress("pho1_eta",        &eta1);
    ntu_DA->SetBranchStatus("pho1_phi",        1); ntu_DA->SetBranchAddress("pho1_phi",        &phi1);
    ntu_DA->SetBranchStatus("pho1_r9",         1); ntu_DA->SetBranchAddress("pho1_r9",         &R91);
    
    ntu_DA->SetBranchStatus("pho2_energy_regr",1); ntu_DA->SetBranchAddress("pho2_energy_regr",&scEReg2);
    ntu_DA->SetBranchStatus("pho2_sceta",      1); ntu_DA->SetBranchAddress("pho2_sceta",      &scEta2);
    ntu_DA->SetBranchStatus("pho2_scphi",      1); ntu_DA->SetBranchAddress("pho2_scphi",      &scPhi2);
    ntu_DA->SetBranchStatus("pho2_eta",        1); ntu_DA->SetBranchAddress("pho2_eta",        &eta2);
    ntu_DA->SetBranchStatus("pho2_phi",        1); ntu_DA->SetBranchAddress("pho2_phi",        &phi2);
    ntu_DA->SetBranchStatus("pho2_r9",         1); ntu_DA->SetBranchAddress("pho2_r9",         &R92);
    
    ntu_MC->SetBranchStatus("pho1_energy_regr",1); ntu_MC->SetBranchAddress("pho1_energy_regr",&scEReg1);
    ntu_MC->SetBranchStatus("pho1_sceta",      1); ntu_MC->SetBranchAddress("pho1_sceta",      &scEta1);
    ntu_MC->SetBranchStatus("pho1_scphi",      1); ntu_MC->SetBranchAddress("pho1_scphi",      &scPhi1);
    ntu_MC->SetBranchStatus("pho1_eta",        1); ntu_MC->SetBranchAddress("pho1_eta",        &eta1);
    ntu_MC->SetBranchStatus("pho1_phi",        1); ntu_MC->SetBranchAddress("pho1_phi",        &phi1);
    ntu_MC->SetBranchStatus("pho1_r9",         1); ntu_MC->SetBranchAddress("pho1_r9",         &R91);
    
    ntu_MC->SetBranchStatus("pho2_energy_regr",1); ntu_MC->SetBranchAddress("pho2_energy_regr",&scEReg2);
    ntu_MC->SetBranchStatus("pho2_sceta",      1); ntu_MC->SetBranchAddress("pho2_sceta",      &scEta2);
    ntu_MC->SetBranchStatus("pho2_scphi",      1); ntu_MC->SetBranchAddress("pho2_scphi",      &scPhi2);
    ntu_MC->SetBranchStatus("pho2_eta",        1); ntu_MC->SetBranchAddress("pho2_eta",        &eta2);
    ntu_MC->SetBranchStatus("pho2_phi",        1); ntu_MC->SetBranchAddress("pho2_phi",        &phi2);
    ntu_MC->SetBranchStatus("pho2_r9",         1); ntu_MC->SetBranchAddress("pho2_r9",         &R92);
  }
  if( useShervinNtuple )
  {
    ntu_DA -> SetBranchStatus("R9Ele",         1);   ntu_DA -> SetBranchAddress("R9Ele",R9);
    ntu_DA -> SetBranchStatus("etaSCEle",      1);   ntu_DA -> SetBranchAddress("etaSCEle",scEta);
    ntu_DA -> SetBranchStatus("phiSCEle",      1);   ntu_DA -> SetBranchAddress("phiSCEle",scPhi);
    ntu_DA -> SetBranchStatus("etaEle",        1);   ntu_DA -> SetBranchAddress("etaEle",eta);
    ntu_DA -> SetBranchStatus("phiEle",        1);   ntu_DA -> SetBranchAddress("phiEle",phi);
    ntu_DA -> SetBranchStatus("rawEnergySCEle",1);   ntu_DA -> SetBranchAddress("rawEnergySCEle",scERaw);
    ntu_DA -> SetBranchStatus("energySCEle",   1);   ntu_DA -> SetBranchAddress("energySCEle",scE);
    if( enCorrType == "stdSC" )
    {
      ntu_DA -> SetBranchStatus("energySCEle",1);   ntu_DA -> SetBranchAddress("energySCEle",scEReg);
    }
    if( enCorrType == "eleTunedRegV3" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorr_ele",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorr_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV3" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorr_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorr_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV4" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV4_ele",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV4_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV4" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV4_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV4_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV5" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV5_ele",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV5_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV5" )
    {
      ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV5_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV5_pho",scEReg);
    }
    ntu_DA -> SetBranchStatus("eleID",1);   ntu_DA -> SetBranchAddress("eleID",eleID);
    
    ntu_MC -> SetBranchStatus("R9Ele",         1);   ntu_MC -> SetBranchAddress("R9Ele",R9);
    ntu_MC -> SetBranchStatus("etaSCEle",      1);   ntu_MC -> SetBranchAddress("etaSCEle",scEta);
    ntu_MC -> SetBranchStatus("phiSCEle",      1);   ntu_MC -> SetBranchAddress("phiSCEle",scPhi);
    ntu_MC -> SetBranchStatus("etaEle",        1);   ntu_MC -> SetBranchAddress("etaEle",eta);
    ntu_MC -> SetBranchStatus("phiEle",        1);   ntu_MC -> SetBranchAddress("phiEle",phi);
    ntu_MC -> SetBranchStatus("rawEnergySCEle",1);   ntu_MC -> SetBranchAddress("rawEnergySCEle",scERaw);
    ntu_MC -> SetBranchStatus("energySCEle",   1);   ntu_MC -> SetBranchAddress("energySCEle",scE);
    if( enCorrType == "stdSC" )
    {
      ntu_MC -> SetBranchStatus("energySCEle",1);   ntu_MC -> SetBranchAddress("energySCEle",scEReg);
    }
    if( enCorrType == "eleTunedRegV3" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorr_ele",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorr_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV3" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorr_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorr_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV4" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV4_ele",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV4_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV4" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV4_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV4_pho",scEReg);
    }
    if( enCorrType == "eleTunedRegV5" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV5_ele",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV5_ele",scEReg);
    }
    if( enCorrType == "phoTunedRegV5" )
    {
      ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV5_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV5_pho",scEReg);
    }
    ntu_MC -> SetBranchStatus("eleID",1);   ntu_MC -> SetBranchAddress("eleID",eleID);
  }
  if( !useGlobeNtuple && !useShervinNtuple )
  {
    if( enCorrType == "stdSC" ) {
      ntu_DA->SetBranchStatus("ele1_scE",1); ntu_DA->SetBranchAddress("ele1_scE",&scEReg1); }
    if( enCorrType == "eleTunedReg" ) {
      ntu_DA->SetBranchStatus("ele1_scE_regression",1); ntu_DA->SetBranchAddress("ele1_scE_regression",&scEReg1); }
    if( enCorrType == "phoTunedReg" ) {
      ntu_DA->SetBranchStatus("ele1_scE_regression_PhotonTuned",1); ntu_DA->SetBranchAddress("ele1_scE_regression_PhotonTuned",&scEReg1); }
    ntu_DA->SetBranchStatus("ele1_scERaw",     1); ntu_DA->SetBranchAddress("ele1_scERaw",     &scERaw1);
    ntu_DA->SetBranchStatus("ele1_scEta",      1); ntu_DA->SetBranchAddress("ele1_scEta",      &scEta1);
    ntu_DA->SetBranchStatus("ele1_scPhi",      1); ntu_DA->SetBranchAddress("ele1_scPhi",      &scPhi1);
    ntu_DA->SetBranchStatus("ele1_eta",        1); ntu_DA->SetBranchAddress("ele1_eta",        &etaFloat1);
    ntu_DA->SetBranchStatus("ele1_phi",        1); ntu_DA->SetBranchAddress("ele1_phi",        &phiFloat1);
    ntu_DA->SetBranchStatus("ele1_e3x3",       1); ntu_DA->SetBranchAddress("ele1_e3x3",       &E3x31);
    
    if( enCorrType == "stdSC" ) {
      ntu_DA->SetBranchStatus("ele2_scE",1); ntu_DA->SetBranchAddress("ele2_scE",&scEReg2); }
    if( enCorrType == "eleTunedReg" ) {
      ntu_DA->SetBranchStatus("ele2_scE_regression",1); ntu_DA->SetBranchAddress("ele2_scE_regression",&scEReg2); }
    if( enCorrType == "phoTunedReg" ) {
      ntu_DA->SetBranchStatus("ele2_scE_regression_PhotonTuned",1); ntu_DA->SetBranchAddress("ele2_scE_regression_PhotonTuned",&scEReg2); }
    ntu_DA->SetBranchStatus("ele2_scERaw",     1); ntu_DA->SetBranchAddress("ele2_scERaw",     &scERaw2);
    ntu_DA->SetBranchStatus("ele2_scEta",      1); ntu_DA->SetBranchAddress("ele2_scEta",      &scEta2);
    ntu_DA->SetBranchStatus("ele2_scPhi",      1); ntu_DA->SetBranchAddress("ele2_scPhi",      &scPhi2);
    ntu_DA->SetBranchStatus("ele2_eta",        1); ntu_DA->SetBranchAddress("ele2_eta",        &etaFloat2);
    ntu_DA->SetBranchStatus("ele2_phi",        1); ntu_DA->SetBranchAddress("ele2_phi",        &phiFloat2);
    ntu_DA->SetBranchStatus("ele2_e3x3",       1); ntu_DA->SetBranchAddress("ele2_e3x3",       &E3x32);
    
    if( enCorrType == "stdSC" ) {
      ntu_MC->SetBranchStatus("ele1_scE",1); ntu_MC->SetBranchAddress("ele1_scE",&scEReg1); }
    if( enCorrType == "eleTunedReg" ) {
      ntu_MC->SetBranchStatus("ele1_scE_regression",1); ntu_MC->SetBranchAddress("ele1_scE_regression",&scEReg1); }
    if( enCorrType == "phoTunedReg" ) {
      ntu_MC->SetBranchStatus("ele1_scE_regression_PhotonTuned",1); ntu_MC->SetBranchAddress("ele1_scE_regression_PhotonTuned",&scEReg1); }
    ntu_MC->SetBranchStatus("ele1_scERaw",     1); ntu_MC->SetBranchAddress("ele1_scERaw",     &scERaw1);
    ntu_MC->SetBranchStatus("ele1_scEta",      1); ntu_MC->SetBranchAddress("ele1_scEta",      &scEta1);
    ntu_MC->SetBranchStatus("ele1_scPhi",      1); ntu_MC->SetBranchAddress("ele1_scPhi",      &scPhi1);
    ntu_MC->SetBranchStatus("ele1_eta",        1); ntu_MC->SetBranchAddress("ele1_eta",        &etaFloat1);
    ntu_MC->SetBranchStatus("ele1_phi",        1); ntu_MC->SetBranchAddress("ele1_phi",        &phiFloat1);
    ntu_MC->SetBranchStatus("ele1_e3x3",       1); ntu_MC->SetBranchAddress("ele1_e3x3",       &E3x31);
    
    if( enCorrType == "stdSC" ) {
      ntu_MC->SetBranchStatus("ele2_scE",1); ntu_MC->SetBranchAddress("ele2_scE",&scEReg2); }
    if( enCorrType == "eleTunedReg" ) {
      ntu_MC->SetBranchStatus("ele2_scE_regression",1); ntu_MC->SetBranchAddress("ele2_scE_regression",&scEReg2); }
    if( enCorrType == "phoTunedReg" ) {
      ntu_MC->SetBranchStatus("ele2_scE_regression_PhotonTuned",1); ntu_MC->SetBranchAddress("ele2_scE_regression_PhotonTuned",&scEReg2); }
    ntu_MC->SetBranchStatus("ele2_scERaw",     1); ntu_MC->SetBranchAddress("ele2_scERaw",     &scERaw2);
    ntu_MC->SetBranchStatus("ele2_scEta",      1); ntu_MC->SetBranchAddress("ele2_scEta",      &scEta2);
    ntu_MC->SetBranchStatus("ele2_scPhi",      1); ntu_MC->SetBranchAddress("ele2_scPhi",      &scPhi2);
    ntu_MC->SetBranchStatus("ele2_eta",        1); ntu_MC->SetBranchAddress("ele2_eta",        &etaFloat2);
    ntu_MC->SetBranchStatus("ele2_phi",        1); ntu_MC->SetBranchAddress("ele2_phi",        &phiFloat2);
    ntu_MC->SetBranchStatus("ele2_e3x3",       1); ntu_MC->SetBranchAddress("ele2_e3x3",       &E3x32);
  }
  
  
  
  //-----------------------------
  // Setup data scale corrections
  std::cout << std::endl;
  std::cout << ">>> Setup data scale corrections" << std::endl;
  
  ScaleCorrector* myScaleCorrector = new ScaleCorrector(runRangeFile);
  
  if( energyScaleCorrType == "shervin" ) myScaleCorrector -> SetShervinRunDepScaleMap(ShervinScaleFile);
  if( energyScaleCorrType == "fabrice" ) myScaleCorrector -> SetIJazZGlobalScaleHisto(IJazZGlobalFolder);
  if( energyScaleCorrType == "fabrice" ) myScaleCorrector -> SetIJazZRunDepScaleHistoMap(IJazZRunDepFolder);
  
  
  
  //-----------------------
  // Setup MC extrasmearing
  std::cout << std::endl;
  std::cout << ">>> Setup MC extrasmearing" << std::endl;
  
  Smearer* mySmearer = new Smearer();
  
  if( energyScaleCorrType == "shervin" ) mySmearer -> SetShervinExtraSmearingMap(ShervinSmearingFile);
  if( energyScaleCorrType == "fabrice" ) mySmearer -> SetIJazZExtraSmearingHisto(IJazZGlobalFolder);
  
  
  
  //--------------------------
  // pileup reweighting for MC
  std::cout << std::endl;
  std::cout << ">>> Setup MC pileup reweighting" << std::endl;
  
  std::map<std::string, TH1F*>* PUWeights = ReadPUWeights(MCGen,runDepFlag,runMin,runMax);
  
  
  
  //-------------------
  // eta/R9 reweighting
  std::cout << std::endl;
  std::cout << ">>> Setup eta/R9 reweighting" << std::endl;
    
  std::string dataDir(getenv("LINEARITY"));
  TFile* etaR9reweightFile = TFile::Open((dataDir+"/data/zee_etaR9reweight.root").c_str());
  
  TH2F* etaR9reweight_lead    = (TH2F*)( etaR9reweightFile->Get("etaR9reweight_lead") );
  TH2F* etaR9reweight_sublead = (TH2F*)( etaR9reweightFile->Get("etaR9reweight_lead") );
  
  
  
  //----------------
  // define outfiles
  
  std::string plotFolderName = outFilePath + "/" + year + "/";
  gSystem->mkdir(plotFolderName.c_str());
  plotFolderName += dataLabel;
  gSystem->mkdir(plotFolderName.c_str());
  
  plotFolderName += "/MZ_";
  plotFolderName += catType; 
  plotFolderName += "_" + (useGlobeNtuple == true ? std::string("globe") : std::string("nonGlobe"));
  plotFolderName += "_" + MCGen + "-" + (runDepFlag == true ? "runDependent" : "allRange");
  plotFolderName += "_" + enCorrType;
  plotFolderName += Form("_Dphi%dp%02d",int(DphiMax),int(DphiMax*100)%100);
  if( MCClosure == true )             plotFolderName += "_MCClosure";
  if( MCHiggs   == true )             plotFolderName += "_MCHiggs";
  if( applyEtaR9Reweighting == true ) plotFolderName += "_etaR9Reweighting";
  if( smearingSyst != 0 )             plotFolderName += Form("_smearingSyst%dp%02d",int(smearingSyst),int(smearingSyst*100)%100);
  plotFolderName += "/";
  gSystem->mkdir(plotFolderName.c_str());
  
  std::string label = "cat" + std::string(Form("%d",category)) + "_" + std::string(Form("%devtsPerPoint",evtsPerPoint));
  outFile = TFile::Open((plotFolderName+"/studyLinearity_MZ_"+label+".root").c_str(),"RECREATE");
  
  
  
  
  
  
  //-----------------
  // Loop over events
  std::cout << std::endl;
  std::cout << ">>> Read data from MC sample" << std::endl;
  
  int nEntries_MC = ntu_MC -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_MC; ++ientry)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   MC entry " << ientry << " / " << nEntries_MC << "\r" << std::flush;
    ntu_MC->GetEntry(ientry);  
    
    
    // variables
    if( !useShervinNtuple )
    {
      eleID1 = 7;
      eleID2 = 7;
    }
    if( !useGlobeNtuple && !useShervinNtuple )
    {
      eta1 = (double)(etaFloat1);
      eta2 = (double)(etaFloat2);
      phi1 = (double)(phiFloat1);
      phi2 = (double)(phiFloat2);
      R91 = E3x31 / scERaw1;
      R92 = E3x32 / scERaw2;
    }
    if( useShervinNtuple )
    {
      R91 = R9[0];
      R92 = R9[1];
      scEta1 = scEta[0];
      scEta2 = scEta[1];
      scPhi1 = scPhi[0];
      scPhi2 = scPhi[1];
      eta1 = eta[0];
      eta2 = eta[1];
      phi1 = phi[0];
      phi2 = phi[1];
      scEReg1 = scEReg[0];
      scEReg2 = scEReg[1];
      eleID1 = eleID[0];
      eleID2 = eleID[1];    
    }
    if( !useGlobeNtuple )
    {
      if( year == "2011" )
      {
        if( fabs(eta1) < 1.5 ) R91 *= 1.005;
        if( fabs(eta1) > 1.5 ) R91 *= 1.004;
        if( fabs(eta2) < 1.5 ) R92 *= 1.005;
        if( fabs(eta2) > 1.5 ) R92 *= 1.004;
      }
      if( catType == "stdCat" ) cat = GetStraightCategory(scEta1,R91,scEta2,R92);
      if( catType == "CiC"    ) cat =   GetHggCiCCategory(scEta1,R91,scEta2,R92);
      if( applyPUWeight )
      {
        std::string periodLabel = getPeriodLabel(runId,runDepFlag,runMin,runMax);
        
        int ibin = (*PUWeights)[periodLabel] -> FindBin( nPU );
        if( ibin <= 1 ) ibin = 1;
        if( ibin >= (*PUWeights)[periodLabel]->GetNbinsX() ) ibin = (*PUWeights)[periodLabel]->GetNbinsX();
        weight = 1. * (*PUWeights)[periodLabel]->GetBinContent(ibin);
      }
    }
    
    float theta1 = 2*atan(exp(-eta1));
    float theta2 = 2*atan(exp(-eta2));
    float Rt1 = sin(theta1);
    float Rt2 = sin(theta2);
    
    
    // selections
    if( !HLTfire ) continue;
    if( fabs(scEta1) >= 2.5000 || fabs(scEta2) >= 2.5000  ) continue;
    if( fabs(scEta1) >  1.4442 && fabs(scEta1) <  1.5660 ) continue;
    if( fabs(scEta2) >  1.4442 && fabs(scEta2) <  1.5660 ) continue;
    if( R91 < 0.0 || R91 >= 1.0 ) continue;
    if( R92 < 0.0 || R92 >= 1.0 ) continue;
    if( scEReg1*Rt1 < 30. ) continue;
    if( scEReg2*Rt2 < 30. ) continue;
    if( ((eleID1 & eleIDBit) != eleIDBit) || ((eleID2 & eleIDBit) != eleIDBit) ) continue;
    if( cat == -1 ) continue;
    if( (category != -1) && (category != cat) ) continue;
    
    
    if( (applyEnergySmearing == true) && (MCClosure == false) )
    {
      float energySmearing1 = gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta1,R91,dataLabel,energySmearingType, 0.));
      float energySmearing2 = gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta2,R92,dataLabel,energySmearingType, 0.));

//       std::cout << " >>> err 0 = " << gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta1,R91,dataLabel,energySmearingType, 0.)) << std::endl;
//       std::cout << " >>> err 1 = " << gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta1,R91,dataLabel,energySmearingType, 1.)) << std::endl;

      scEReg1 *= energySmearing1;
      scEReg2 *= energySmearing2;
    }
    
    if( applyEtaR9Reweighting == true )
    {
      float leadEta = scEta1;
      float leadR9 = R91;
      float subleadEta = scEta2;
      float subleadR9 = R92;
      if( scEReg1*Rt1 < scEReg2*Rt2 )
      {
        leadEta = scEta2;
        leadR9 = R92;
        subleadEta = scEta1;
        subleadR9 = R91;
      }
      float etaR9Weight_lead    = etaR9reweight_lead    -> GetBinContent(etaR9reweight_lead->FindBin(leadEta,leadR9));
      float etaR9Weight_sublead = etaR9reweight_sublead -> GetBinContent(etaR9reweight_sublead->FindBin(subleadEta,subleadR9));    
      weight *= etaR9Weight_lead;
      weight *= etaR9Weight_sublead;
    }
    
    TLorentzVector p1; p1.SetPtEtaPhiE(scEReg1*Rt1,eta1,phi1,scEReg1);
    TLorentzVector p2; p2.SetPtEtaPhiE(scEReg2*Rt2,eta2,phi2,scEReg2);
    float mee = sqrt( 4. * scEReg1 * scEReg2 * pow(sin(0.5*(p1.Vect()).Angle(p2.Vect())),2) ) / 91.18;
    float Dphi = deltaPhi(scPhi1,scPhi2);
        
    
    // apply cuts
    if( MCClosure == true ) if( ientry%2 == 1 ) continue;
    if( mee < 0. || mee > 2.5 ) continue;
    if( Dphi > DphiMax ) continue;
    
    
    // fill vectors
    scEta1_MC.push_back(scEta1);
    scEta2_MC.push_back(scEta2);
    
    R91_MC.push_back(R91);
    R92_MC.push_back(R92);
    
    Zpt_MC.push_back((p1+p2).Pt());
    mee_MC.push_back(mee);
    Ht_MC.push_back(scEReg1*Rt1 + scEReg2*Rt2);
    Et1_MC.push_back(scEReg1*Rt1);
    Et2_MC.push_back(scEReg2*Rt2);
    weight_MC.push_back(weight);
  }
  std::cout << std::endl;
  
  
  
  std::cout << ">>> Read data from DATA sample" << std::endl;
  
  int nEntries_DA = ntu_DA -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_DA; ientry++)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading DATA entry " << ientry << " / " << nEntries_DA << "\r" << std::flush;
    ntu_DA->GetEntry(ientry);  
    
    
    // variables
    if( !useShervinNtuple )
    {
      eleID1 = 7;
      eleID2 = 7;
    }
    if( !useGlobeNtuple && !useShervinNtuple )
    {
      eta1 = (double)(etaFloat1);
      eta2 = (double)(etaFloat2);
      phi1 = (double)(phiFloat1);
      phi2 = (double)(phiFloat2);
      R91 = E3x31 / scERaw1;
      R92 = E3x32 / scERaw2;
    }
    if( useShervinNtuple )
    {
      R91 = R9[0];
      R92 = R9[1];
      scEta1 = scEta[0];
      scEta2 = scEta[1];
      scPhi1 = scPhi[0];
      scPhi2 = scPhi[1];
      eta1 = eta[0];
      eta2 = eta[1];
      phi1 = phi[0];
      phi2 = phi[1];
      scEReg1 = scEReg[0];
      scEReg2 = scEReg[1];
      eleID1 = eleID[0];
      eleID2 = eleID[1];    
    }
    if( !useGlobeNtuple )
    {
      if( year == "2011" )
      {
        if( fabs(eta1) < 1.5 ) R91 *= 1.005;
        if( fabs(eta1) > 1.5 ) R91 *= 1.004;
        if( fabs(eta2) < 1.5 ) R92 *= 1.005;
        if( fabs(eta2) > 1.5 ) R92 *= 1.004;
      }
      if( catType == "stdCat" ) cat = GetStraightCategory(scEta1,R91,scEta2,R92);
      if( catType == "CiC"    ) cat =   GetHggCiCCategory(scEta1,R91,scEta2,R92);
      if( applyPUWeight )
      {
        std::string periodLabel = getPeriodLabel(runId,runDepFlag,runMin,runMax);
        
        int ibin = (*PUWeights)[periodLabel] -> FindBin( nPU );
        if( ibin <= 1 ) ibin = 1;
        if( ibin >= (*PUWeights)[periodLabel]->GetNbinsX() ) ibin = (*PUWeights)[periodLabel]->GetNbinsX();
        weight = 1. * (*PUWeights)[periodLabel]->GetBinContent(ibin);
      }
    }
    
    float theta1 = 2*atan(exp(-eta1));
    float theta2 = 2*atan(exp(-eta2));
    float Rt1 = sin(theta1);
    float Rt2 = sin(theta2);
    
    
    // selections
    if( !HLTfire ) continue;
    if( fabs(scEta1) >= 2.5000 || fabs(scEta2) >= 2.5000  ) continue;
    if( fabs(scEta1) >  1.4442 && fabs(scEta1) <  1.5660 ) continue;
    if( fabs(scEta2) >  1.4442 && fabs(scEta2) <  1.5660 ) continue;
    if( R91 < 0.0 || R91 >= 1.0 ) continue;
    if( R92 < 0.0 || R92 >= 1.0 ) continue;
    if( scEReg1*Rt1 < 30. ) continue;
    if( scEReg2*Rt2 < 30. ) continue;
    if( ((eleID1 & eleIDBit) != eleIDBit) || ((eleID2 & eleIDBit) != eleIDBit) ) continue;
    if( cat == -1 ) continue;
    if( (category != -1) && (category != cat) ) continue;
    
    
    if( MCClosure == true )
    {
      scEReg1 *= ( f_scaleVsEt -> Eval(scEReg1*Rt1) );
      scEReg2 *= ( f_scaleVsEt -> Eval(scEReg2*Rt2) );
    }
    
    if( (applyEnergyScaleCorr == true) && (MCClosure == false) )
    {
      scEReg1 *= myScaleCorrector->GetScaleCorrection(scEta1,R91,runId,dataLabel,energyScaleCorrType, 0.);
      scEReg2 *= myScaleCorrector->GetScaleCorrection(scEta2,R92,runId,dataLabel,energyScaleCorrType, 0.);

//       std::cout << " err 0 = " << myScaleCorrector->GetScaleCorrection(scEta1,R91,runId,dataLabel,energyScaleCorrType, 0.) << std::endl;
//       std::cout << " err -1 = " << myScaleCorrector->GetScaleCorrection(scEta1,R91,runId,dataLabel,energyScaleCorrType, -1.) << std::endl;
    }
    
    if( applyEtaR9Reweighting == true )
    {
      float leadEta = scEta1;
      float leadR9 = R91;
      float subleadEta = scEta2;
      float subleadR9 = R92;
      if( scEReg1*Rt1 < scEReg2*Rt2 )
      {
        leadEta = scEta2;
        leadR9 = R92;
        subleadEta = scEta1;
        subleadR9 = R91;
      }
      float etaR9Weight_lead    = etaR9reweight_lead    -> GetBinContent(etaR9reweight_lead->FindBin(leadEta,leadR9));
      float etaR9Weight_sublead = etaR9reweight_sublead -> GetBinContent(etaR9reweight_sublead->FindBin(subleadEta,subleadR9));    
      weight *= etaR9Weight_lead;
      weight *= etaR9Weight_sublead;
    }
    
    TLorentzVector p1; p1.SetPtEtaPhiE(scEReg1*Rt1,eta1,phi1,scEReg1);
    TLorentzVector p2; p2.SetPtEtaPhiE(scEReg2*Rt2,eta2,phi2,scEReg2);
    float mee = sqrt( 4. * scEReg1 * scEReg2 * pow(sin(0.5*(p1.Vect()).Angle(p2.Vect())),2) ) / 91.18;
    float Dphi = deltaPhi(scPhi1,scPhi2);
        
    
    // apply cuts
    if( MCClosure == true ) if( ientry%2 == 0 ) continue;
    if( mee < 0. || mee > 2.5 ) continue;
    if( Dphi > DphiMax ) continue;
    
    
    // fill vectors
    scEta1_DA.push_back(scEta1);
    scEta2_DA.push_back(scEta2);
    
    R91_DA.push_back(R91);
    R92_DA.push_back(R92);

    Zpt_DA.push_back((p1+p2).Pt());
    mee_DA.push_back(mee);
    Ht_DA.push_back(scEReg1*Rt1 + scEReg2*Rt2);
    
    mee_fit_DA.push_back(mee);
    mee_gausFit_DA.push_back(mee);
    mee_mean_DA.push_back(mee);
    mee_recursiveMean_DA.push_back(mee);
    mee_smallestInterval_DA.push_back(mee);
    
    Et1_DA.push_back(scEReg1*Rt1);
    Et2_DA.push_back(scEReg2*Rt2);
    Et1_fit_DA.push_back(scEReg1*Rt1);
    Et2_fit_DA.push_back(scEReg2*Rt2);
    Et1_gausFit_DA.push_back(scEReg1*Rt1);
    Et2_gausFit_DA.push_back(scEReg2*Rt2);
    Et1_mean_DA.push_back(scEReg1*Rt1);
    Et2_mean_DA.push_back(scEReg2*Rt2);
    Et1_recursiveMean_DA.push_back(scEReg1*Rt1);
    Et2_recursiveMean_DA.push_back(scEReg2*Rt2);
    Et1_smallestInterval_DA.push_back(scEReg1*Rt1);
    Et2_smallestInterval_DA.push_back(scEReg2*Rt2);
    
    weight_DA.push_back(weight);
  }
  std::cout << std::endl;
  
  
  
  
  //------------
  // sort events
  std::cout << std::endl;
  std::cout << ">>> sort DATA events vs. Ht" << std::endl;
  
  int nEntries = Ht_DA.size();
  int nSavePts = 0;
  std::vector<SorterLC> sortedEntries;

  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    SorterLC dummy;
    dummy.laserCorr = Ht_DA.at(ientry);
    dummy.entry = ientry;
    sortedEntries.push_back(dummy);
    nSavePts++;   
  }
  
  std::cout << ">>>>>> Sorting variable " << "Ht" << std::endl;
  std::cout << ">>>>>> Effective entries: " << nSavePts << std::endl;
  std::cout << ">>>>>> sortedEntries.size(): " << sortedEntries.size() << std::endl;
  std::sort(sortedEntries.begin(),sortedEntries.end(),SorterLC());
  
  
  
  
  //------------
  // define bins
  std::cout << std::endl;
  std::cout << ">>> define bins" << std::endl;
   
  HtBinEdges = new std::vector<double>;


  if(evtsPerPoint > 0){
    HtBinEdges -> push_back( Ht_DA.at(sortedEntries.at(0).entry) );
  
    int nBinTempPts = 0;
    for(int iSaved = 0; iSaved < nSavePts; ++iSaved)
      {
	++nBinTempPts;
	
	if( nBinTempPts == evtsPerPoint )
	  {
	    HtBinEdges -> push_back( Ht_DA.at(sortedEntries.at(iSaved).entry) );
	    nBinTempPts = 0;
	  }
      }
    HtBinEdges -> push_back( Ht_DA.at(sortedEntries.at(nSavePts-1).entry) );
  }
  else{
    for(unsigned int posVec = 0; posVec< extHtBinEdges.size(); ++posVec)   
      HtBinEdges->push_back( extHtBinEdges.at(posVec) );
  }


  
  nHtBins = HtBinEdges->size() - 1;
  for(unsigned int i = 0; i < nHtBins; ++i)
    std::cout << ">>> Ht bin " << i << ":   [" << HtBinEdges->at(i) << "," << HtBinEdges->at(i+1) << "]" << std::endl;
  std::cout << std::endl;
  
  
  
  EtBinEdges = new std::vector<double>;
  for(unsigned int HtBinEdgeIt = 0; HtBinEdgeIt < HtBinEdges->size(); ++HtBinEdgeIt)
    EtBinEdges -> push_back( 0.5 * HtBinEdges->at(HtBinEdgeIt) );
  nEtBins = EtBinEdges->size()-1;
  
  for(unsigned int i = 0; i < nEtBins; ++i)
    std::cout << ">>> Et bin " << i << ":   [" << EtBinEdges->at(i) << "," << EtBinEdges->at(i+1) << "]" << std::endl;
  std::cout << std::endl;
  
  HtBinEdgesDouble = new double[HtBinEdges->size()];
  for(unsigned int i = 0; i < HtBinEdges->size(); ++i) HtBinEdgesDouble[i] = HtBinEdges->at(i);
  EtBinEdgesDouble = new double[EtBinEdges->size()];
  for(unsigned int i = 0; i < EtBinEdges->size(); ++i) EtBinEdgesDouble[i] = EtBinEdges->at(i);
  
  
  
  //------------------
  // define histograms
  std::cout << std::endl;
  std::cout << ">>> define histograms" << std::endl;
  
  TH1F* h_scEta_MC = new TH1F("h_scEta_MC","",500,-2.5,2.5);
  h_scEta_MC -> Sumw2();
  TH1F* h_scEta_DA = new TH1F("h_scEta_DA","",500,-2.5,2.5);
  h_scEta_DA -> Sumw2();
  TH1F* h_R9_MC = new TH1F("h_R9_MC","",500,-1.,1.);
  h_R9_MC -> Sumw2();
  TH1F* h_R9_DA = new TH1F("h_R9_DA","",500,-1.,1.);
  h_R9_DA -> Sumw2();
  TH1F* h_Ht_MC = new TH1F("h_Ht_MC","",500,0.,500.);
  h_Ht_MC -> Sumw2();
  TH1F* h_Ht_DA = new TH1F("h_Ht_DA","",500,0.,500.);
  h_Ht_DA -> Sumw2();
  TH1F* h_Zpt_MC = new TH1F("h_Zpt_MC","",300,0.,300.);
  h_Zpt_MC -> Sumw2();
  TH1F* h_Zpt_DA = new TH1F("h_Zpt_DA","",300,0.,300.);
  h_Zpt_DA -> Sumw2();
  TH1F* h_mee_MC = new TH1F("h_mee_MC","",480,60.,120.);
  h_mee_MC -> Sumw2();
  TH1F* h_mee_DA = new TH1F("h_mee_DA","",480,60.,120.);
  h_mee_DA -> Sumw2();
  
  TGraphAsymmErrors* scale_fit_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_fit_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_fit_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_gausFit_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_gausFit_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_gausFit_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_mean_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_mean_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_mean_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_recursiveMean_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_recursiveMean_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_recursiveMean_DAOverMC = new TGraphAsymmErrors();
  
  TGraphAsymmErrors* scale_smallestInterval_MC       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_smallestInterval_DA       = new TGraphAsymmErrors();
  TGraphAsymmErrors* scale_smallestInterval_DAOverMC = new TGraphAsymmErrors();
  
  std::vector<double>* mee_HtBin_MC    = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_MC = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_DA                      = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_fit_DA                  = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_gausFit_DA              = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_mean_DA                 = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_recursiveMean_DA        = new std::vector<double>[nHtBins];
  std::vector<double>* mee_HtBin_smallestInterval_DA     = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_DA                   = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_fit_DA               = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_gausFit_DA           = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_mean_DA              = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_recursiveMean_DA     = new std::vector<double>[nHtBins];
  std::vector<double>* weight_HtBin_smallestInterval_DA  = new std::vector<double>[nHtBins];
    
  TH1F** h_Zpt_HtBin_MC = new TH1F*[nHtBins];
  TH1F** h_Zpt_HtBin_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_MC = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_fit_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_gausFit_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_mean_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_recursiveMean_DA = new TH1F*[nHtBins];
  TH1F** h_mee_HtBin_smallestInterval_DA = new TH1F*[nHtBins];
  
  TF1** f_gausFit_HtBin_MC = new TF1*[nHtBins];
  TF1** f_gausFit_HtBin_DA = new TF1*[nHtBins];
  
  h_Ht_HtBin_MC = new TH1F*[nHtBins];
  
  h_Et_EtBin_MC = new TH1F*[nEtBins];
  h_Et_EtBin_DA = new TH1F*[nEtBins];
  h_Et_EtBin_fit_DA = new TH1F*[nEtBins];
  h_Et_EtBin_gausFit_DA = new TH1F*[nEtBins];
  h_Et_EtBin_mean_DA = new TH1F*[nEtBins];
  h_Et_EtBin_recursiveMean_DA = new TH1F*[nEtBins];
  h_Et_EtBin_smallestInterval_DA = new TH1F*[nEtBins];
  
  for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
  {
    char histoName[50];
    
    sprintf(histoName,"h_mee_HtBin%d_MC",HtBin);
    h_mee_HtBin_MC[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
    //h_mee_HtBin_MC[HtBin] -> SetFillColor(kGreen+2);
    //h_mee_HtBin_MC[HtBin] -> SetFillStyle(3004);
    //h_mee_HtBin_MC[HtBin] -> SetMarkerStyle(7);
    //h_mee_HtBin_MC[HtBin] -> SetMarkerColor(kGreen+2);
    //h_mee_HtBin_MC[HtBin] -> SetLineColor(kGreen+2);
    h_mee_HtBin_MC[HtBin]->Sumw2();
    
    sprintf(histoName, "h_Zpt_HtBin%d_MC",HtBin);
    h_Zpt_HtBin_MC[HtBin] = new TH1F(histoName,"",300,0.,300.);
    h_Zpt_HtBin_MC[HtBin] -> SetLineColor(kGreen+2);
    h_Zpt_HtBin_MC[HtBin] -> Sumw2();
    
    sprintf(histoName, "Ht_HtBin%d_MC",HtBin);
    h_Ht_HtBin_MC[HtBin] = new TH1F(histoName,"",5000,0.,1000.);
    h_Ht_HtBin_MC[HtBin] -> SetLineColor(kGreen+2);
    h_Ht_HtBin_MC[HtBin] -> Sumw2();
  }
  
  for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
  {
    char histoName[50];
    
    sprintf(histoName, "Et_EtBin%d_MC",EtBin);
    h_Et_EtBin_MC[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
    h_Et_EtBin_MC[EtBin] -> SetLineColor(kGreen+2);
    h_Et_EtBin_MC[EtBin] -> Sumw2();
  }
  
  
  
  
  //----------------
  // fill histograms
  std::cout << std::endl;
  std::cout << ">>> fill histograms" << std::endl;
  
  int MCEntries = mee_MC.size();
  for(int ientry = 0; ientry < MCEntries; ++ientry)
  {   
    if( (ientry%100000 == 0) ) std::cout << "reading   MC entry " << ientry << " / " << MCEntries << "\r" << std::flush;
    
    
    int HtBin = MyFindBin(Ht_MC.at(ientry),HtBinEdges);
    if( HtBin == -1 ) continue;
    
    int EtBin1 = MyFindBin(Et1_MC.at(ientry),EtBinEdges);
    int EtBin2 = MyFindBin(Et2_MC.at(ientry),EtBinEdges);
    if( EtBin1 == -1 ) continue;
    if( EtBin2 == -1 ) continue;
    
    mee_HtBin_MC[HtBin].push_back( mee_MC.at(ientry) );
    weight_HtBin_MC[HtBin].push_back( weight_MC.at(ientry) );
    
    h_Zpt_HtBin_MC[HtBin] -> Fill( Zpt_MC.at(ientry),weight_MC.at(ientry) );
    h_mee_HtBin_MC[HtBin] -> Fill( mee_MC.at(ientry),weight_MC.at(ientry) );
    h_Ht_HtBin_MC[HtBin]  -> Fill(  Ht_MC.at(ientry),weight_MC.at(ientry) );
    
    h_Et_EtBin_MC[EtBin1] -> Fill( Et1_MC.at(ientry),weight_MC.at(ientry) );
    h_Et_EtBin_MC[EtBin2] -> Fill( Et2_MC.at(ientry),weight_MC.at(ientry) );
    
    h_scEta_MC -> Fill( scEta1_MC.at(ientry),weight_MC.at(ientry) );
    h_scEta_MC -> Fill( scEta2_MC.at(ientry),weight_MC.at(ientry) );
    
    h_R9_MC -> Fill( R91_MC.at(ientry),weight_MC.at(ientry) );
    h_R9_MC -> Fill( R92_MC.at(ientry),weight_MC.at(ientry) );
    
    h_Ht_MC -> Fill( Ht_MC.at(ientry),weight_MC.at(ientry) );
    h_Ht_MC -> Fill( Ht_MC.at(ientry),weight_MC.at(ientry) );
    
    h_Zpt_MC -> Fill( Zpt_MC.at(ientry),weight_MC.at(ientry) );
    h_Zpt_MC -> Fill( Zpt_MC.at(ientry),weight_MC.at(ientry) );
    
    h_mee_MC -> Fill( mee_MC.at(ientry)*91.18,weight_MC.at(ientry) );
    h_mee_MC -> Fill( mee_MC.at(ientry)*91.18,weight_MC.at(ientry) );
  }
  std::cout << std::endl;
  
  
  
  for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
  {
    double x = h_Ht_HtBin_MC[HtBin]->GetMean();
    
    scale_fit_DA -> SetPoint(HtBin,x,1.);
    scale_fit_MC -> SetPoint(HtBin,x,1.);
    scale_fit_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_gausFit_DA -> SetPoint(HtBin,x,1.);
    scale_gausFit_MC -> SetPoint(HtBin,x,1.);
    scale_gausFit_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_mean_DA -> SetPoint(HtBin,x,1.);
    scale_mean_MC -> SetPoint(HtBin,x,1.);
    scale_mean_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_recursiveMean_DA -> SetPoint(HtBin,x,1.);
    scale_recursiveMean_MC -> SetPoint(HtBin,x,1.);
    scale_recursiveMean_DAOverMC -> SetPoint(HtBin,x,1.);
    
    scale_smallestInterval_DA -> SetPoint(HtBin,x,1.);
    scale_smallestInterval_MC -> SetPoint(HtBin,x,1.);
    scale_smallestInterval_DAOverMC -> SetPoint(HtBin,x,1.);
  }
  
  
  
  for(int step = 1; step < nSteps+1; ++step)
  {
    std::cout << std::endl;
    std::cout << "****** step " << step << " ******" << std::endl;
    
    
    TProfile* p_avgEtCorr_fit              = new TProfile("p_avgEtCorr_fit",             "",nEtBins,EtBinEdgesDouble);
    TProfile* p_avgEtCorr_gausFit          = new TProfile("p_avgEtCorr_gausFit",         "",nEtBins,EtBinEdgesDouble);
    TProfile* p_avgEtCorr_mean             = new TProfile("p_avgEtCorr_mean",            "",nEtBins,EtBinEdgesDouble);
    TProfile* p_avgEtCorr_recursiveMean    = new TProfile("p_avgEtCorr_recursiveMean",   "",nEtBins,EtBinEdgesDouble);
    TProfile* p_avgEtCorr_smallestInterval = new TProfile("p_avgEtCorr_smallestInterval","",nEtBins,EtBinEdgesDouble);
    
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
        mee_HtBin_DA[HtBin].clear();
        mee_HtBin_fit_DA[HtBin].clear();
        mee_HtBin_gausFit_DA[HtBin].clear();
        mee_HtBin_mean_DA[HtBin].clear();
        mee_HtBin_recursiveMean_DA[HtBin].clear();
        mee_HtBin_smallestInterval_DA[HtBin].clear();
        weight_HtBin_DA[HtBin].clear();
        weight_HtBin_fit_DA[HtBin].clear();
        weight_HtBin_gausFit_DA[HtBin].clear();
        weight_HtBin_mean_DA[HtBin].clear();
        weight_HtBin_recursiveMean_DA[HtBin].clear();
        weight_HtBin_smallestInterval_DA[HtBin].clear();
    }	
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
      char histoName[50];
      
      sprintf(histoName,"h_Zpt_HtBin%d_DA_step%d",HtBin,step);
      h_Zpt_HtBin_DA[HtBin] = new TH1F(histoName,"",300,0.,300.);
      h_Zpt_HtBin_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_DA_step%d",HtBin,step);
      h_mee_HtBin_DA[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
      h_mee_HtBin_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_fit_DA_step%d",HtBin,step);
      h_mee_HtBin_fit_DA[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
      h_mee_HtBin_fit_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_gausFit_DA_step%d",HtBin,step);
      h_mee_HtBin_gausFit_DA[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
      h_mee_HtBin_gausFit_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_mean_DA_step%d",HtBin,step);
      h_mee_HtBin_mean_DA[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
      h_mee_HtBin_mean_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_recursiveMean_DA_step%d",HtBin,step);
      h_mee_HtBin_recursiveMean_DA[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
      h_mee_HtBin_recursiveMean_DA[HtBin] -> Sumw2();
      
      sprintf(histoName,"h_mee_HtBin%d_smallestInterval_DA_step%d",HtBin,step);
      h_mee_HtBin_smallestInterval_DA[HtBin] = new TH1F(histoName,"",nBinsMee,meeMin,meeMax);
      h_mee_HtBin_smallestInterval_DA[HtBin] -> Sumw2();
    }
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      char histoName[50];
      
      sprintf(histoName,"h_Et_EtBin%d_DA_step%d",EtBin,step);
      h_Et_EtBin_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_fit_DA_step%d",EtBin,step);
      h_Et_EtBin_fit_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_fit_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_fit_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_gausFit_DA_step%d",EtBin,step);
      h_Et_EtBin_gausFit_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_gausFit_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_gausFit_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_mean_DA_step%d",EtBin,step);
      h_Et_EtBin_mean_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_mean_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_mean_DA[EtBin] -> Sumw2();
            
      sprintf(histoName,"h_Et_EtBin%d_recursiveMean_DA_step%d",EtBin,step);
      h_Et_EtBin_recursiveMean_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_recursiveMean_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_recursiveMean_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_Et_EtBin%d_smallestInterval_DA_step%d",EtBin,step);
      h_Et_EtBin_smallestInterval_DA[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
      h_Et_EtBin_smallestInterval_DA[EtBin] -> SetLineColor(kRed+2);
      h_Et_EtBin_smallestInterval_DA[EtBin] -> Sumw2();
    }
    
    
    int DAEntries = mee_DA.size();
    for(int ientry = 0; ientry < DAEntries; ++ientry)
    {   
      if( (ientry%100000 == 0) ) std::cout << "reading DATA entry " << ientry << " / " << DAEntries << "\r" << std::flush;
      
      
      double k1 = 1.;
      double k2 = 1.;
      double Et1 = Et1_DA.at(ientry)/k1;
      double Et2 = Et2_DA.at(ientry)/k2;
      int EtBin1 = MyFindBin(Et1,EtBinEdges);
      int EtBin2 = MyFindBin(Et2,EtBinEdges);
      double Ht = Et1 + Et2;
      int HtBin = MyFindBin(Ht,HtBinEdges);
      double mee = mee_DA.at(ientry)/sqrt(k1*k2);
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_DA[HtBin].push_back( mee );
        weight_HtBin_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        
        h_Zpt_HtBin_DA[HtBin] -> Fill( Zpt_DA.at(ientry),weight_DA.at(ientry) );
        h_mee_HtBin_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      
      k1 = MyEval(scale_fit_DAOverMC,2.*Et1_fit_DA.at(ientry));
      k2 = MyEval(scale_fit_DAOverMC,2.*Et2_fit_DA.at(ientry));
      Et1 = Et1_fit_DA.at(ientry)/k1;
      Et2 = Et2_fit_DA.at(ientry)/k2;
      Et1_fit_DA.at(ientry) = Et1;
      Et2_fit_DA.at(ientry) = Et2;
      p_avgEtCorr_fit -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_fit -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_fit_DA.at(ientry)/sqrt(k1*k2);
      mee_fit_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_fit_DA[HtBin].push_back( mee );
        weight_HtBin_fit_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_fit_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_fit_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_fit_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_mean_DAOverMC,2.*Et1_mean_DA.at(ientry));
      k2 = MyEval(scale_mean_DAOverMC,2.*Et2_mean_DA.at(ientry));
      Et1 = Et1_mean_DA.at(ientry)/k1;
      Et2 = Et2_mean_DA.at(ientry)/k2;
      Et1_mean_DA.at(ientry) = Et1;
      Et2_mean_DA.at(ientry) = Et2;
      p_avgEtCorr_mean -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_mean -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_mean_DA.at(ientry)/sqrt(k1*k2);
      mee_mean_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_mean_DA[HtBin].push_back( mee );
        weight_HtBin_mean_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_mean_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_mean_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_mean_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_gausFit_DAOverMC,2.*Et1_gausFit_DA.at(ientry));
      k2 = MyEval(scale_gausFit_DAOverMC,2.*Et2_gausFit_DA.at(ientry));
      Et1 = Et1_gausFit_DA.at(ientry)/k1;
      Et2 = Et2_gausFit_DA.at(ientry)/k2;
      Et1_gausFit_DA.at(ientry) = Et1;
      Et2_gausFit_DA.at(ientry) = Et2;
      p_avgEtCorr_gausFit -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_gausFit -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_gausFit_DA.at(ientry)/sqrt(k1*k2);
      mee_gausFit_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_gausFit_DA[HtBin].push_back( mee );
        weight_HtBin_gausFit_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_gausFit_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_gausFit_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_gausFit_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_recursiveMean_DAOverMC,2.*Et1_recursiveMean_DA.at(ientry));
      k2 = MyEval(scale_recursiveMean_DAOverMC,2.*Et2_recursiveMean_DA.at(ientry));
      Et1 = Et1_recursiveMean_DA.at(ientry)/k1;
      Et2 = Et2_recursiveMean_DA.at(ientry)/k2;
      Et1_recursiveMean_DA.at(ientry) = Et1;
      Et2_recursiveMean_DA.at(ientry) = Et2;
      p_avgEtCorr_recursiveMean -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_recursiveMean -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_recursiveMean_DA.at(ientry)/sqrt(k1*k2);
      mee_recursiveMean_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_recursiveMean_DA[HtBin].push_back( mee );
        weight_HtBin_recursiveMean_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_recursiveMean_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_recursiveMean_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_recursiveMean_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      k1 = MyEval(scale_smallestInterval_DAOverMC,2.*Et1_smallestInterval_DA.at(ientry));
      k2 = MyEval(scale_smallestInterval_DAOverMC,2.*Et2_smallestInterval_DA.at(ientry));
      Et1 = Et1_smallestInterval_DA.at(ientry)/k1;
      Et2 = Et2_smallestInterval_DA.at(ientry)/k2;
      Et1_smallestInterval_DA.at(ientry) = Et1;
      Et2_smallestInterval_DA.at(ientry) = Et2;
      p_avgEtCorr_smallestInterval -> Fill(Et1,Et1/Et1_DA.at(ientry));
      p_avgEtCorr_smallestInterval -> Fill(Et2,Et2/Et2_DA.at(ientry));
      EtBin1 = MyFindBin(Et1,EtBinEdges);
      EtBin2 = MyFindBin(Et2,EtBinEdges);
      Ht = Et1 + Et2;
      HtBin = MyFindBin(Ht,HtBinEdges);
      mee = mee_smallestInterval_DA.at(ientry)/sqrt(k1*k2);
      mee_smallestInterval_DA.at(ientry) = mee;
      if( HtBin != -1 && EtBin1 != -1 && EtBin2 != -1 )
      {
        mee_HtBin_smallestInterval_DA[HtBin].push_back( mee );
        weight_HtBin_smallestInterval_DA[HtBin].push_back( weight_DA.at(ientry) );
        
        h_mee_HtBin_smallestInterval_DA[HtBin] -> Fill( mee,weight_DA.at(ientry) );
        
        h_Et_EtBin_smallestInterval_DA[EtBin1] -> Fill( Et1,weight_DA.at(ientry) );
        h_Et_EtBin_smallestInterval_DA[EtBin2] -> Fill( Et2,weight_DA.at(ientry) );
      }
      
      
      if( step == 1 )
      {
        h_scEta_DA -> Fill( scEta1_DA.at(ientry),weight_DA.at(ientry) );
        h_scEta_DA -> Fill( scEta2_DA.at(ientry),weight_DA.at(ientry) );
        
        h_R9_DA -> Fill( R91_DA.at(ientry),weight_DA.at(ientry) );
        h_R9_DA -> Fill( R92_DA.at(ientry),weight_DA.at(ientry) );
        
        h_Ht_DA -> Fill( Ht_DA.at(ientry),weight_DA.at(ientry) );
        h_Ht_DA -> Fill( Ht_DA.at(ientry),weight_DA.at(ientry) );
        
        h_Zpt_DA -> Fill( Zpt_DA.at(ientry),weight_DA.at(ientry) );
        h_Zpt_DA -> Fill( Zpt_DA.at(ientry),weight_DA.at(ientry) );
        
        h_mee_DA -> Fill( mee_DA.at(ientry)*91.18,weight_DA.at(ientry) );
        h_mee_DA -> Fill( mee_DA.at(ientry)*91.18,weight_DA.at(ientry) );
      }
    }
    std::cout << std::endl;
    
    
    
    
    //---------------------
    // find scale estimator
    std::cout << std::endl;
    std::cout << ">>> find scale estimator" << std::endl;
    
    std::vector<double> smallestIntervalMins_MC;
    std::vector<double> smallestIntervalMaxs_MC;
    std::vector<double> smallestIntervalMins_DA;
    std::vector<double> smallestIntervalMaxs_DA;
    
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
      std::cout << ">>> HtBin::" << HtBin << std::endl;
      std::cout << ">>>>>> nEvents:                  " << mee_HtBin_DA[HtBin].size() << std::endl;
      
      double x = h_Ht_HtBin_MC[HtBin]->GetMean();
      double ex = h_Ht_HtBin_MC[HtBin]->GetMeanError();
      double exlow = ex;
      double exhig = ex;
      
      h_Zpt_HtBin_MC[HtBin] -> Scale(1./h_Zpt_HtBin_MC[HtBin]->Integral());
      h_Zpt_HtBin_DA[HtBin] -> Scale(1./h_Zpt_HtBin_DA[HtBin]->Integral());
      
      h_mee_HtBin_MC[HtBin] -> Scale(1./h_mee_HtBin_MC[HtBin]->Integral());
      h_mee_HtBin_DA[HtBin] -> Scale(1./h_mee_HtBin_DA[HtBin]->Integral());
      h_mee_HtBin_fit_DA[HtBin] -> Scale(1./h_mee_HtBin_fit_DA[HtBin]->Integral());
      h_mee_HtBin_gausFit_DA[HtBin] -> Scale(1./h_mee_HtBin_gausFit_DA[HtBin]->Integral());
      h_mee_HtBin_mean_DA[HtBin] -> Scale(1./h_mee_HtBin_mean_DA[HtBin]->Integral());
      h_mee_HtBin_recursiveMean_DA[HtBin] -> Scale(1./h_mee_HtBin_recursiveMean_DA[HtBin]->Integral());
      h_mee_HtBin_smallestInterval_DA[HtBin] -> Scale(1./h_mee_HtBin_smallestInterval_DA[HtBin]->Integral());

      
      // fit
      std::cout << ">>>>>> nEvents_fit:              " << mee_HtBin_fit_DA[HtBin].size() << std::endl;
      if( mee_HtBin_fit_DA[HtBin].size() > 3 )
      {
        double scale_MC = 1.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindTemplateFit(scale_DA,scaleErr_DA,h_mee_HtBin_MC[HtBin],h_mee_HtBin_fit_DA[HtBin]);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_fit_MC -> SetPoint(HtBin,x,scale_MC);
        scale_fit_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_fit_DA -> SetPoint(HtBin,x,scale_DA);
        scale_fit_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_fit_DAOverMC -> SetPoint(HtBin,x,y);
        scale_fit_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // gausFit
      std::cout << ">>>>>> nEvents_gausFit:          " << mee_HtBin_gausFit_DA[HtBin].size() << std::endl;
      if( mee_HtBin_gausFit_DA[HtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        char funcName_MC[50];
        sprintf(funcName_MC,"f_gausFit_HtBin%d_step%d_MC",HtBin,step);
        char funcName_DA[50];
        sprintf(funcName_DA,"f_gausFit_HtBin%d_step%d_DA",HtBin,step);
        FindGausFit(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],nBinsMee,meeMin,meeMax,&(f_gausFit_HtBin_MC[HtBin]),std::string(funcName_MC),1.);
        FindGausFit(scale_DA,scaleErr_DA,mee_HtBin_gausFit_DA[HtBin],weight_HtBin_gausFit_DA[HtBin],nBinsMee,meeMin,meeMax,&(f_gausFit_HtBin_DA[HtBin]),std::string(funcName_DA),1.);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_gausFit_MC -> SetPoint(HtBin,x,scale_MC);
        scale_gausFit_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_gausFit_DA -> SetPoint(HtBin,x,scale_DA);
        scale_gausFit_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_gausFit_DAOverMC -> SetPoint(HtBin,x,y);
        scale_gausFit_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // mean
      std::cout << ">>>>>> nEvents_mean:             " << mee_HtBin_mean_DA[HtBin].size() << std::endl;
      if( mee_HtBin_mean_DA[HtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindMean(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin]);
        FindMean(scale_DA,scaleErr_DA,mee_HtBin_mean_DA[HtBin],weight_HtBin_mean_DA[HtBin]);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_mean_MC -> SetPoint(HtBin,x,scale_MC);
        scale_mean_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_mean_DA -> SetPoint(HtBin,x,scale_DA);
        scale_mean_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_mean_DAOverMC -> SetPoint(HtBin,x,y);
        scale_mean_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // recursive mean
      std::cout << ">>>>>> nEvents_recursiveMean:    " << mee_HtBin_recursiveMean_DA[HtBin].size() << std::endl;
      if( mee_HtBin_recursiveMean_DA[HtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindRecursiveMean(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],5./91.18,0.0001);
        FindRecursiveMean(scale_DA,scaleErr_DA,mee_HtBin_recursiveMean_DA[HtBin],weight_HtBin_recursiveMean_DA[HtBin],5./91.18,0.0001);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
	std::cout << " FindRecursiveMean HtBin = " << HtBin << " x = " << x << std::endl; 
	std::cout << " MC scale = " << scale_MC << std::endl; 
	std::cout << " DA scale = " << scale_DA << std::endl; 
	std::cout << " DA/MC scale = " << y << std::endl; 

	if(scale_DA == 0 || scale_MC == 0){
	  FindRecursiveMean(scale_MC,scaleErr_MC,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],5./h_mee_HtBin_MC[HtBin]->GetMean(),0.0001);
	  FindRecursiveMean(scale_DA,scaleErr_DA,mee_HtBin_recursiveMean_DA[HtBin],weight_HtBin_recursiveMean_DA[HtBin],5./h_mee_HtBin_recursiveMean_DA[HtBin]->GetMean(),0.0001);
	  y = scale_DA / scale_MC;
	  eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
	  eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
	std::cout << " FindRecursiveMean HtBin = " << HtBin << " x = " << x << std::endl; 
	std::cout << " MC scale = " << scale_MC << std::endl; 
	std::cout << " DA scale = " << scale_DA << std::endl; 
	std::cout << " DA/MC scale = " << y << std::endl; 

	}

        scale_recursiveMean_MC -> SetPoint(HtBin,x,scale_MC);
        scale_recursiveMean_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_recursiveMean_DA -> SetPoint(HtBin,x,scale_DA);
        scale_recursiveMean_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);

	std::cout << " HtBin = " << HtBin << " x = " << x << " y = " << y << std::endl; 
	std::cout << " HtBin = " << HtBin << " exlow = " << exlow << " exhig = " << exhig << " eylow = " << eylow << " eyhig = " << eyhig << std::endl; 
        scale_recursiveMean_DAOverMC -> SetPoint(HtBin,x,y);
	//float locErrorGraph = graphError_RecursiveMean->Eval(mee_HtBin_recursiveMean_DA[HtBin].size());
	//scale_recursiveMean_DAOverMC->SetPointError(HtBin,exlow,exhig,locErrorGraph,locErrorGraph);
	scale_recursiveMean_DAOverMC->SetPointError(HtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // smallest interval
      std::cout << ">>>>>> nEvents_smallestInterval: " << mee_HtBin_smallestInterval_DA[HtBin].size() << std::endl;
      if( mee_HtBin_smallestInterval_DA[HtBin].size() > 3 )
      {
        double min_MC;
        double max_MC;
        double min_DA;
        double max_DA;
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindSmallestInterval(scale_MC,scaleErr_MC,min_MC,max_MC,0.,2.,0.0005,mee_HtBin_MC[HtBin],weight_HtBin_MC[HtBin],0.5);
        FindSmallestInterval(scale_DA,scaleErr_DA,min_DA,max_DA,0.,2.,0.0005,mee_HtBin_smallestInterval_DA[HtBin],weight_HtBin_smallestInterval_DA[HtBin],0.5);
        smallestIntervalMins_MC.push_back(min_MC);
        smallestIntervalMaxs_MC.push_back(max_MC);
        smallestIntervalMins_DA.push_back(min_DA);
        smallestIntervalMaxs_DA.push_back(max_DA);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_smallestInterval_MC -> SetPoint(HtBin,x,scale_MC);
        scale_smallestInterval_MC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_smallestInterval_DA -> SetPoint(HtBin,x,scale_DA);
        scale_smallestInterval_DA -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
        scale_smallestInterval_DAOverMC -> SetPoint(HtBin,x,y);
	//        scale_smallestInterval_DAOverMC -> SetPointError(HtBin,exlow,exhig,eylow,eyhig);
	float locErrorGraph = graphError_SmallestInterval->Eval(mee_HtBin_smallestInterval_DA[HtBin].size());
        scale_smallestInterval_DAOverMC -> SetPointError(HtBin,exlow,exhig,locErrorGraph,locErrorGraph);
      }
    }
    
    

    
    outFile -> cd();
    
    char dirName[50];
    sprintf(dirName,"step%d",step);
    TDirectory* baseDir = outFile -> mkdir(dirName);
    TDirectory* subDir = NULL;
    baseDir -> cd();
    
    scale_fit_MC           -> Write("scale_fit_MC");
    scale_fit_DA           -> Write("scale_fit_DA");
    scale_fit_DAOverMC     -> Write("scale_fit_DAOverMC");
    p_avgEtCorr_fit        -> Write();
    
    scale_gausFit_MC           -> Write("scale_gausFit_MC");
    scale_gausFit_DA           -> Write("scale_gausFit_DA");
    scale_gausFit_DAOverMC     -> Write("scale_gausFit_DAOverMC");
    p_avgEtCorr_gausFit        -> Write();
    
    scale_mean_MC           -> Write("scale_mean_MC");
    scale_mean_DA           -> Write("scale_mean_DA");
    scale_mean_DAOverMC     -> Write("scale_mean_DAOverMC");
    p_avgEtCorr_mean        -> Write();
    
    scale_recursiveMean_MC           -> Write("scale_recursiveMean_MC");
    scale_recursiveMean_DA           -> Write("scale_recursiveMean_DA");
    scale_recursiveMean_DAOverMC     -> Write("scale_recursiveMean_DAOverMC");
    p_avgEtCorr_recursiveMean        -> Write();
    
    scale_smallestInterval_MC           -> Write("scale_smallestInterval_MC");
    scale_smallestInterval_DA           -> Write("scale_smallestInterval_DA");
    scale_smallestInterval_DAOverMC     -> Write("scale_smallestInterval_DAOverMC");
    p_avgEtCorr_smallestInterval        -> Write();
    
    
    delete p_avgEtCorr_fit;
    delete p_avgEtCorr_gausFit;
    delete p_avgEtCorr_mean;
    delete p_avgEtCorr_recursiveMean;
    delete p_avgEtCorr_smallestInterval;
    
    
    baseDir -> cd();
    subDir = baseDir -> mkdir("mee_HtBin");
    subDir -> cd();
    
    std::string outputPdf_DAMC  = plotFolderName + "h_mee_HtBin_DAOverMC_" + catType+"_cat"+ Form("%d",category) + ".pdf";
    std::string outputPdf2_DAMC = plotFolderName + "h_Zpt_HtBin_DAOverMC_" + catType+"_cat"+ Form("%d",category) + ".pdf";
    
    std::string outputPdf_MC = plotFolderName + "h_mee_HtBin_MC_"   + catType+"_cat"+ Form("%d",category) + ".pdf";
    std::string outputPdf_DA = plotFolderName + "h_mee_HtBin_DATA_" + catType+"_cat"+ Form("%d",category) + ".pdf";
    
    for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
    {
      h_mee_HtBin_MC[HtBin] -> Write();
      h_mee_HtBin_DA[HtBin] -> Write();
      h_mee_HtBin_fit_DA[HtBin] -> Write();
      h_mee_HtBin_gausFit_DA[HtBin] -> Write();
      h_mee_HtBin_mean_DA[HtBin] -> Write();
      h_mee_HtBin_recursiveMean_DA[HtBin] -> Write();
      h_mee_HtBin_smallestInterval_DA[HtBin] -> Write();
      
      f_gausFit_HtBin_MC[HtBin] -> Write();
      f_gausFit_HtBin_DA[HtBin] -> Write();
      
      
      
      if( step == 1 )
      {
        TCanvas* c_DAOverMC = new TCanvas("c_Zpt");
        c_DAOverMC -> cd();
        c_DAOverMC -> SetGridx();
        c_DAOverMC -> SetGridy();
        
        char axisTitle[50];
        sprintf(axisTitle,"p_{T}(ee) [GeV]   -   H_{T} #in [%d,%d]",int(HtBinEdges->at(HtBin)),int(HtBinEdges->at(HtBin+1)));
        h_Zpt_HtBin_MC[HtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_Zpt_HtBin_MC[HtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_Zpt_HtBin_MC[HtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_Zpt_HtBin_MC[HtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_Zpt_HtBin_MC[HtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_Zpt_HtBin_MC[HtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_Zpt_HtBin_MC[HtBin] -> SetLineColor(kBlack);
        h_Zpt_HtBin_MC[HtBin] -> SetLineWidth(1);
        h_Zpt_HtBin_MC[HtBin] -> GetXaxis() -> SetRangeUser(0.,300.);
        float maximum = std::max(h_Zpt_HtBin_MC[HtBin]->GetMaximum(),h_Zpt_HtBin_DA[HtBin]->GetMaximum());
        h_Zpt_HtBin_MC[HtBin] -> SetMaximum( 1.1*maximum );
        
        h_Zpt_HtBin_MC[HtBin] -> Draw("hist");
        h_Zpt_HtBin_DA[HtBin] -> Draw("P,same");
        
        if( HtBin == 0 )         c_DAOverMC -> Print((outputPdf2_DAMC+"[").c_str());
        c_DAOverMC -> Print(outputPdf2_DAMC.c_str());
        if( HtBin == nHtBins-1 ) c_DAOverMC -> Print((outputPdf2_DAMC+"]").c_str());
        delete c_DAOverMC;
        
        
        
        c_DAOverMC = new TCanvas("c_mee");
        c_DAOverMC -> cd();
        c_DAOverMC -> SetGridx();
        c_DAOverMC -> SetGridy();
        
        sprintf(axisTitle,"m_{ee}/m_{Z}^{PDG}   -   H_{T} #in [%d,%d]",int(HtBinEdges->at(HtBin)),int(HtBinEdges->at(HtBin+1)));
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_mee_HtBin_MC[HtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_MC[HtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_MC[HtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_MC[HtBin] -> SetLineColor(kBlack);
        h_mee_HtBin_MC[HtBin] -> SetLineWidth(1);
        h_mee_HtBin_MC[HtBin] -> GetXaxis() -> SetRangeUser(0.80,1.19999);
        maximum = std::max(h_mee_HtBin_MC[HtBin]->GetMaximum(),h_mee_HtBin_DA[HtBin]->GetMaximum());
        h_mee_HtBin_MC[HtBin] -> SetMaximum( 1.1*maximum );
        
        h_mee_HtBin_MC[HtBin] -> Draw("hist");
        h_mee_HtBin_DA[HtBin] -> Draw("P,same");
        
        if( HtBin == 0 )         c_DAOverMC -> Print((outputPdf_DAMC+"[").c_str());
        c_DAOverMC -> Print(outputPdf_DAMC.c_str());
        if( HtBin == nHtBins-1 ) c_DAOverMC -> Print((outputPdf_DAMC+"]").c_str());
        delete c_DAOverMC;
        
        
        
        TCanvas* c_MC = new TCanvas("c_MC");
        c_MC -> cd();
        c_MC -> SetGridx();
        c_MC -> SetGridy();
        
        double x,y;
        
        scale_mean_MC -> GetPoint(HtBin,x,y);
        TArrow* line_mean_MC = new TArrow(y,0.,y,h_mee_HtBin_MC[HtBin]->GetBinContent(h_mee_HtBin_MC[HtBin]->FindBin(y)));
        line_mean_MC -> SetLineColor(kBlack);
        line_mean_MC -> SetLineWidth(2);
        
        scale_recursiveMean_MC -> GetPoint(HtBin,x,y);
        TArrow* line_recursiveMean_MC = new TArrow(y,0.,y,h_mee_HtBin_MC[HtBin]->GetMaximum());
        line_recursiveMean_MC -> SetLineColor(kBlue);
        line_recursiveMean_MC -> SetLineWidth(2);
        TArrow* line_recursiveMean_min_MC = new TArrow(y-5./91.18,0.,y-5./91.18,h_mee_HtBin_MC[HtBin]->GetMaximum());
        line_recursiveMean_min_MC -> SetLineColor(kBlue);
        line_recursiveMean_min_MC -> SetLineWidth(2);
        TArrow* line_recursiveMean_max_MC = new TArrow(y+5./91.18,0.,y+5./91.18,h_mee_HtBin_MC[HtBin]->GetMaximum());
        line_recursiveMean_max_MC -> SetLineColor(kBlue);
        line_recursiveMean_max_MC -> SetLineWidth(2);
        
        scale_smallestInterval_MC -> GetPoint(HtBin,x,y);
        TArrow* line_smallestInterval_MC = new TArrow(y,0.,y,h_mee_HtBin_MC[HtBin]->GetBinContent(h_mee_HtBin_MC[HtBin]->FindBin(y)));
        line_smallestInterval_MC -> SetLineColor(kRed);
        line_smallestInterval_MC -> SetLineWidth(2);
        
        h_mee_HtBin_MC[HtBin] -> Draw("HIST");
        
        TH1F* clone = (TH1F*)( h_mee_HtBin_MC[HtBin]->Clone() );
        for(int bin = 1; bin < clone->GetNbinsX(); ++bin)
          if( (clone->GetBinCenter(bin) < smallestIntervalMins_MC.at(HtBin)) ||
              (clone->GetBinCenter(bin) > smallestIntervalMaxs_MC.at(HtBin)) )
            clone -> SetBinContent(bin,0.);
        clone -> SetFillColor(kYellow);
        clone -> SetLineWidth(0);
        clone -> Draw("HIST,same");
        
        line_mean_MC              -> Draw("same");
        line_recursiveMean_MC     -> Draw("same");
        line_recursiveMean_min_MC -> Draw("same");
        line_recursiveMean_max_MC -> Draw("same");
        line_smallestInterval_MC  -> Draw("same");
        f_gausFit_HtBin_MC[HtBin] -> Draw("same");
        
        
        TCanvas* c_DA = new TCanvas("c_DA");
        c_DA -> cd();
        c_DA -> SetGridx();
        c_DA -> SetGridy();
        
        sprintf(axisTitle,"m_{ee} (GeV/c^{2})   -   H_{T} #in [%d,%d]",int(HtBinEdges->at(HtBin)),int(HtBinEdges->at(HtBin+1)));
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_mee_HtBin_DA[HtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_DA[HtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_DA[HtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_mee_HtBin_DA[HtBin] -> SetLineColor(kBlack);
        h_mee_HtBin_DA[HtBin] -> SetLineWidth(1);
        h_mee_HtBin_DA[HtBin] -> GetXaxis() -> SetRangeUser(0.80,1.19999);
        h_mee_HtBin_DA[HtBin] -> SetMaximum( 1.1*h_mee_HtBin_DA[HtBin]->GetMaximum() );
        
        
        scale_mean_DA -> GetPoint(HtBin,x,y);
        TArrow* line_mean_DA = new TArrow(y,0.,y,h_mee_HtBin_DA[HtBin]->GetBinContent(h_mee_HtBin_DA[HtBin]->FindBin(y)));
        line_mean_DA -> SetLineColor(kBlack);
        line_mean_DA -> SetLineWidth(2);
        
        scale_recursiveMean_DA -> GetPoint(HtBin,x,y);
        TArrow* line_recursiveMean_DA = new TArrow(y,0.,y,h_mee_HtBin_DA[HtBin]->GetMaximum());
        line_recursiveMean_DA -> SetLineColor(kBlue);
        line_recursiveMean_DA -> SetLineWidth(2);
        TArrow* line_recursiveMean_min_DA = new TArrow(y-5./91.18,0.,y-5./91.18,h_mee_HtBin_DA[HtBin]->GetMaximum());
        line_recursiveMean_min_DA -> SetLineColor(kBlue);
        line_recursiveMean_min_DA -> SetLineWidth(2);
        TArrow* line_recursiveMean_max_DA = new TArrow(y+5./91.18,0.,y+5./91.18,h_mee_HtBin_DA[HtBin]->GetMaximum());
        line_recursiveMean_max_DA -> SetLineColor(kBlue);
        line_recursiveMean_max_DA -> SetLineWidth(2);
        
        scale_smallestInterval_DA -> GetPoint(HtBin,x,y);
        TArrow* line_smallestInterval_DA = new TArrow(y,0.,y,h_mee_HtBin_DA[HtBin]->GetBinContent(h_mee_HtBin_DA[HtBin]->FindBin(y)));
        line_smallestInterval_DA -> SetLineColor(kRed);
        line_smallestInterval_DA -> SetLineWidth(2);
        
        h_mee_HtBin_DA[HtBin] -> Draw("HIST");
        
        clone = (TH1F*)( h_mee_HtBin_DA[HtBin]->Clone() );
        for(int bin = 1; bin < clone->GetNbinsX(); ++bin)
          if( (clone->GetBinCenter(bin) < smallestIntervalMins_DA.at(HtBin)) ||
              (clone->GetBinCenter(bin) > smallestIntervalMaxs_DA.at(HtBin)) )
            clone -> SetBinContent(bin,0.);
        clone -> SetFillColor(kYellow);
        clone -> SetLineWidth(0);
        clone -> Draw("HIST,same");
        
        line_mean_DA              -> Draw("same");
        line_recursiveMean_DA     -> Draw("same");
        line_recursiveMean_min_DA -> Draw("same");
        line_recursiveMean_max_DA -> Draw("same");
        line_smallestInterval_DA  -> Draw("same");
        f_gausFit_HtBin_DA[HtBin] -> Draw("same");
        
        
        if( HtBin == 0 )
        {
          c_MC -> Print((outputPdf_MC+"[").c_str());
          c_DA -> Print((outputPdf_DA+"[").c_str());
        }
        {
          c_MC -> RedrawAxis();
          c_MC -> Print(outputPdf_MC.c_str());
          c_DA -> RedrawAxis();
          c_DA -> Print(outputPdf_DA.c_str());
        }
        if( HtBin == nHtBins-1 )
        {
          c_MC -> RedrawAxis();
          c_MC -> Print((outputPdf_MC+"]").c_str());
          c_DA -> RedrawAxis();
          c_DA -> Print((outputPdf_DA+"]").c_str());
        }
        
        delete c_MC;
        delete c_DA;
      }
    }
    
    
    
    baseDir -> cd();
    subDir = baseDir -> mkdir("Et_EtBin");
    subDir -> cd();
    
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      h_Et_EtBin_MC[EtBin] -> Scale(1./h_Et_EtBin_MC[EtBin]->Integral());
      h_Et_EtBin_MC[EtBin] -> Write();
      
      h_Et_EtBin_DA[EtBin] -> Scale(1./h_Et_EtBin_DA[EtBin]->Integral());
      h_Et_EtBin_DA[EtBin] -> Write();
      h_Et_EtBin_fit_DA[EtBin] -> Scale(1./h_Et_EtBin_fit_DA[EtBin]->Integral());
      h_Et_EtBin_fit_DA[EtBin] -> Write();
      h_Et_EtBin_gausFit_DA[EtBin] -> Scale(1./h_Et_EtBin_gausFit_DA[EtBin]->Integral());
      h_Et_EtBin_gausFit_DA[EtBin] -> Write();
      h_Et_EtBin_mean_DA[EtBin] -> Scale(1./h_Et_EtBin_mean_DA[EtBin]->Integral());
      h_Et_EtBin_mean_DA[EtBin] -> Write();
      h_Et_EtBin_recursiveMean_DA[EtBin] -> Scale(1./h_Et_EtBin_recursiveMean_DA[EtBin]->Integral());
      h_Et_EtBin_recursiveMean_DA[EtBin] -> Write();
      h_Et_EtBin_smallestInterval_DA[EtBin] -> Scale(1./h_Et_EtBin_smallestInterval_DA[EtBin]->Integral());
      h_Et_EtBin_smallestInterval_DA[EtBin] -> Write();
    }
    
    outFile -> cd();
  }
  
  
  
  outFile -> mkdir("Ht_HtBin");
  outFile -> cd("Ht_HtBin");
  
  for(unsigned int HtBin = 0; HtBin < nHtBins; ++HtBin)
  {
    h_Ht_HtBin_MC[HtBin] -> Scale(1./h_Ht_HtBin_MC[HtBin]->Integral());
    h_Ht_HtBin_MC[HtBin] -> Write();
  }
  
  outFile -> cd();
  
  f_scaleVsEt -> Write();
  f_scaleVs2Et -> Write();
  f_scaleVsAvgEt -> Write();
  f_invScaleVsEt -> Write();
  
  h_scEta_MC -> Scale(1./h_scEta_MC->Integral());
  h_scEta_MC -> Write();
  h_scEta_DA -> Scale(1./h_scEta_DA->Integral());
  h_scEta_DA -> Write();
  
  h_R9_MC -> Scale(1./h_R9_MC->Integral());
  h_R9_MC -> Write();
  h_R9_DA -> Scale(1./h_R9_DA->Integral());
  h_R9_DA -> Write();
  
  h_Ht_MC -> Scale(1./h_Ht_MC->Integral());
  h_Ht_MC -> Write();
  h_Ht_DA -> Scale(1./h_Ht_DA->Integral());
  h_Ht_DA -> Write();
  
  h_Zpt_MC -> Scale(1./h_Zpt_MC->Integral());
  h_Zpt_MC -> Write();
  h_Zpt_DA -> Scale(1./h_Zpt_DA->Integral());
  h_Zpt_DA -> Write();
  
  h_mee_MC -> Scale(1./h_mee_MC->Integral());
  h_mee_MC -> Write();
  h_mee_DA -> Scale(1./h_mee_DA->Integral());
  h_mee_DA -> Write();
  
  outFile -> Close();
  
  
  
  return 0;
}






double invScaleVsEt(double* x, double* par)
{
  double xx = x[0];
  
  if( xx  < EtBinEdges->at(0)       ) return 1. / f_scaleVsEt -> Eval( h_Et_EtBin_DA[0]        ->GetMean() );
  if( xx >= EtBinEdges->at(nEtBins) ) return 1. / f_scaleVsEt -> Eval( h_Et_EtBin_DA[nEtBins-1]->GetMean() );
  
  int EtBin = MyFindBin(xx,EtBinEdges);
  return 1. / f_scaleVsEt -> Eval(h_Et_EtBin_DA[EtBin]->GetMean());
  
  return 0.;
}



double scaleVsAvgEt(double* x, double* par)
{
  double xx = x[0];
  
  if( xx  < EtBinEdges->at(0)       ) return f_scaleVsEt -> Eval( h_Et_EtBin_DA[0]        ->GetMean() );
  if( xx >= EtBinEdges->at(nEtBins) ) return f_scaleVsEt -> Eval( h_Et_EtBin_DA[nEtBins-1]->GetMean() );
  
  int EtBin = MyFindBin(xx,EtBinEdges);
  return f_scaleVsEt -> Eval(h_Et_EtBin_DA[EtBin]->GetMean());
  
  return 0.;
}

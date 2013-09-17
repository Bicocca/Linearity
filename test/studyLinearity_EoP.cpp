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
  if(argc != 2)
  {
    std::cerr << ">>> studyLinearity_EoP::usage: " << argv[0] << " configFileName" << std::endl;
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
  
  int nBinsEoP  = gConfigParser -> readIntOption("Options::nBinsEoP");
  double EoPMin = gConfigParser -> readDoubleOption("Options::EoPMin");
  double EoPMax = gConfigParser -> readDoubleOption("Options::EoPMax");

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
  
  //f_scaleVsEt  = new TF1("f_scaleVsEt", "1.+0.000",0., 1000.);
  //f_scaleVs2Et = new TF1("f_scaleVsEt", "1.+0.000",0., 1000.);
  
  //f_scaleVsEt  = new TF1("f_scaleVsEt", "1.+0.002",0., 1000.);
  //f_scaleVs2Et = new TF1("f_scaleVsEt", "1.+0.002",0., 1000.);
  
  f_scaleVsEt = new TF1("f_scaleVsEt", "1. + [0] * (1 - exp(-[1] * (x-45.)) )",0., 1000.);
  f_scaleVsEt -> SetParameters(7.50e-03,2.00e-02);
  
  f_scaleVs2Et = new TF1("f_scaleVs2Et", "1. + [0] * (1 - exp(-[1] * (0.5*x-45.)) )",0., 1000.);
  f_scaleVs2Et -> SetParameters(6.00e-03,2.00e-02);
  
  f_scaleVsAvgEt = new TF1("f_scaleVsAvgEt",scaleVsAvgEt,0.,1000.,0);
  f_invScaleVsEt = new TF1("f_invScaleVsEt",invScaleVsEt,0.,1000.,0);
  
  
  
  
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
  std::cout << ">>> define branch addresses" << std::endl;
  
  // vectors
  std::vector<double> weight_MC, Ht_MC;
  std::vector<double> weight_DA, Ht_DA;
  std::vector<double> Zpt_MC, EoP_MC, Et_MC, scEta_MC, R9_MC;
  std::vector<double> Zpt_DA, EoP_DA, Et_DA, scEta_DA, R9_DA;
  std::vector<double> EoP_fit_DA, Et_fit_DA;
  std::vector<double> EoP_gausFit_DA, Et_gausFit_DA;
  std::vector<double> EoP_mean_DA, Et_mean_DA;
  std::vector<double> EoP_recursiveMean_DA, Et_recursiveMean_DA;
  std::vector<double> EoP_smallestInterval_DA, Et_smallestInterval_DA;
  
  // global variables
  int runId,nPU,cat,cat1,cat2;
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
    ntu_DA -> SetBranchStatus("nvtx",  1); ntu_DA -> SetBranchAddress("nvtx",  &nPU);
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
  float tkP[2];
  
  double eta1, phi1;
  double eta2, phi2;
  float tkP1, scERaw1, scEReg1, scEta1, scPhi1, etaFloat1, phiFloat1, E3x31, R91;
  float tkP2, scERaw2, scEReg2, scEta2, scPhi2, etaFloat2, phiFloat2, E3x32, R92;
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
    ntu_DA -> SetBranchStatus("eleID",       1);   ntu_DA -> SetBranchAddress("eleID",eleID);
    ntu_DA -> SetBranchStatus("pAtVtxGsfEle",1);   ntu_DA -> SetBranchAddress("pAtVtxGsfEle",tkP);
    
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
    ntu_MC -> SetBranchStatus("pAtVtxGsfEle",1);   ntu_MC -> SetBranchAddress("pAtVtxGsfEle",tkP);
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
    ntu_DA->SetBranchStatus("ele1_tkP",        1); ntu_DA->SetBranchAddress("ele1_tkP",        &tkP1);
    
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
    ntu_DA->SetBranchStatus("ele2_tkP",        1); ntu_DA->SetBranchAddress("ele2_tkP",        &tkP2);
    
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
    ntu_MC->SetBranchStatus("ele1_tkP",        1); ntu_MC->SetBranchAddress("ele1_tkP",        &tkP1);
    
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
    ntu_MC->SetBranchStatus("ele2_tkP",        1); ntu_MC->SetBranchAddress("ele2_tkP",        &tkP2);
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
  
  plotFolderName += "/EoP_";
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
  outFile = TFile::Open((plotFolderName+"/studyLinearity_EoP_"+label+".root").c_str(),"RECREATE");
  
  
  
  
  
  
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
    
    
    // define variables
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
      tkP1 = tkP[0];
      tkP2 = tkP[1];
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
      if( catType == "stdCat" )
      {
        cat1 = GetSingleCategory(scEta1,R91);
        cat2 = GetSingleCategory(scEta2,R92);
      }
      if( catType == "CiC"    )
      {
        cat1 = GetHggCiCCategory(scEta1,R91,scEta2,R92);
        cat2 = GetHggCiCCategory(scEta1,R91,scEta2,R92);
      }
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
    if( cat1 == -1 && cat2 == -1 ) continue;
    if( (category != -1) && (category != cat1) && (category != cat2) ) continue;
    
    
    if( (applyEnergySmearing == true) && (MCClosure == false) )
    {
      float energySmearing1 = gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta1,R91,dataLabel,energySmearingType));
      float energySmearing2 = gRandom->Gaus(1.,mySmearer->GetExtraSmearing(scEta2,R92,dataLabel,energySmearingType));
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
    float EoP1 = scEReg1 / tkP1;
    float EoP2 = scEReg2 / tkP2;
    float Dphi = deltaPhi(scPhi1,scPhi2);
        
    
    // apply cuts
    if( MCClosure == true ) if( ientry%2 == 1 ) continue;
    if( EoP1 < 0. || EoP1 > 5. ) continue;
    if( EoP2 < 0. || EoP2 > 5. ) continue;
    if( Dphi > DphiMax ) continue;
    
    
    // fill vectors
    if( (category != -1) && (category == cat1) )
    {
      scEta_MC.push_back(scEta1);
      R9_MC.push_back(R91);
      Zpt_MC.push_back((p1+p2).Pt());
      EoP_MC.push_back(EoP1);
      Ht_MC.push_back(scEReg1*Rt1 + scEReg2*Rt2);
      Et_MC.push_back(scEReg1*Rt1);
      weight_MC.push_back(weight);
    }
    if( (category != -1) && (category == cat2) )
    {
      scEta_MC.push_back(scEta2);
      R9_MC.push_back(R92);
      Zpt_MC.push_back((p1+p2).Pt());
      EoP_MC.push_back(EoP2);
      Ht_MC.push_back(scEReg1*Rt1 + scEReg2*Rt2);
      Et_MC.push_back(scEReg2*Rt2);
      weight_MC.push_back(weight);
    }
  }
  std::cout << std::endl;
  
  
  
  std::cout << ">>> Read data from DATA sample" << std::endl;
  
  int nEntries_DA = ntu_DA -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_DA; ientry++)
  {
    if( maxEntries != -1 && ientry == maxEntries ) break;
    if( ientry%100000 == 0 ) std::cout << ">>>>>> reading DATA entry " << ientry << " / " << nEntries_DA << "\r" << std::flush;
    ntu_DA->GetEntry(ientry);
    
    
    // define variables
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
      tkP1 = tkP[0];
      tkP2 = tkP[1];
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
      if( catType == "stdCat" )
      {
        cat1 = GetSingleCategory(scEta1,R91);
        cat2 = GetSingleCategory(scEta2,R92);
      }
      if( catType == "CiC"    )
      {
        cat1 = GetHggCiCCategory(scEta1,R91,scEta2,R92);
        cat2 = GetHggCiCCategory(scEta1,R91,scEta2,R92);
      }
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
    if( cat1 == -1 && cat2 == -1 ) continue;
    if( (category != -1) && (category != cat1) && (category != cat2) ) continue;
    
    
    if( MCClosure == true )
    {
      scEReg1 *= ( f_scaleVsEt -> Eval(scEReg1*Rt1) );
      scEReg2 *= ( f_scaleVsEt -> Eval(scEReg2*Rt2) );
    }
    
    if( (applyEnergyScaleCorr == true) && (MCClosure == false) )
    {
      scEReg1 *= myScaleCorrector->GetScaleCorrection(scEta1,R91,runId,dataLabel,energyScaleCorrType);
      scEReg2 *= myScaleCorrector->GetScaleCorrection(scEta2,R92,runId,dataLabel,energyScaleCorrType);
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
    float EoP1 = scEReg1 / tkP1;
    float EoP2 = scEReg2 / tkP2;
    float Dphi = deltaPhi(scPhi1,scPhi2);    
    
    
    // apply cuts
    if( MCClosure == true ) if( ientry%2 == 0 ) continue;
    if( EoP1 < 0. || EoP1 > 5. ) continue;
    if( EoP2 < 0. || EoP2 > 5. ) continue;
    if( Dphi > DphiMax ) continue;
    
    
    // fill vectors
    if( (category != -1) && (category == cat1) )
    {
      scEta_DA.push_back(scEta1);
      R9_DA.push_back(R91);
      Zpt_DA.push_back((p1+p2).Pt());
      EoP_DA.push_back(EoP1);
      Ht_DA.push_back(scEReg1*Rt1 + scEReg2*Rt2);
      
      EoP_fit_DA.push_back(EoP1);
      EoP_gausFit_DA.push_back(EoP1);
      EoP_mean_DA.push_back(EoP1);
      EoP_recursiveMean_DA.push_back(EoP1);
      EoP_smallestInterval_DA.push_back(EoP1);
      
      Et_DA.push_back(scEReg1*Rt1);
      Et_fit_DA.push_back(scEReg1*Rt1);
      Et_gausFit_DA.push_back(scEReg1*Rt1);
      Et_mean_DA.push_back(scEReg1*Rt1);
      Et_recursiveMean_DA.push_back(scEReg1*Rt1);
      Et_smallestInterval_DA.push_back(scEReg1*Rt1);
      
      weight_DA.push_back(weight);
    }
    if( (category != -1) && (category == cat2) )
    {
      scEta_DA.push_back(scEta2);
      R9_DA.push_back(R92);
      Zpt_DA.push_back((p1+p2).Pt());
      EoP_DA.push_back(EoP2);
      Ht_DA.push_back(scEReg1*Rt1 + scEReg2*Rt2);
      
      EoP_fit_DA.push_back(EoP2);
      EoP_gausFit_DA.push_back(EoP2);
      EoP_mean_DA.push_back(EoP2);
      EoP_recursiveMean_DA.push_back(EoP2);
      EoP_smallestInterval_DA.push_back(EoP2);
      
      Et_DA.push_back(scEReg2*Rt2);
      Et_fit_DA.push_back(scEReg2*Rt2);
      Et_gausFit_DA.push_back(scEReg2*Rt2);
      Et_mean_DA.push_back(scEReg2*Rt2);
      Et_recursiveMean_DA.push_back(scEReg2*Rt2);
      Et_smallestInterval_DA.push_back(scEReg2*Rt2);
      
      weight_DA.push_back(weight);
    }
  }
  std::cout << std::endl;
  
  
  
  
  //------------
  // sort events
  std::cout << std::endl;
  std::cout << ">>> sort MC events vs. Et" << std::endl;
  
  int nEntries = Et_MC.size();
  int nSavePts = 0;
  std::vector<SorterLC> sortedEntries;

  for(int ientry = 0; ientry < nEntries; ++ientry)
  {
    SorterLC dummy;
    dummy.laserCorr = Et_MC.at(ientry);
    dummy.entry = ientry;
    sortedEntries.push_back(dummy);
    nSavePts++;   
  }
  
  std::cout << ">>>>>> Sorting variable " << "Et" << std::endl;
  std::cout << ">>>>>> Effective entries: " << nSavePts << std::endl;
  std::cout << ">>>>>> sortedEntries.size(): " << sortedEntries.size() << std::endl;
  std::sort(sortedEntries.begin(),sortedEntries.end(),SorterLC());
  
  
  
  
  //------------
  // define bins
  std::cout << std::endl;
  std::cout << ">>> define bins" << std::endl;
  
  EtBinEdges = new std::vector<double>;
  EtBinEdges -> push_back( Et_MC.at(sortedEntries.at(0).entry) );
  
  int nBinTempPts = 0;
  for(int iSaved = 0; iSaved < nSavePts; ++iSaved)
  {
    ++nBinTempPts;
    
    if( nBinTempPts == evtsPerPoint )
    {
      EtBinEdges -> push_back( Et_MC.at(sortedEntries.at(iSaved).entry) );
      nBinTempPts = 0;
    }
  }
  EtBinEdges -> push_back( Et_MC.at(sortedEntries.at(nSavePts-1).entry) );
  
  nEtBins = EtBinEdges->size() - 1;
  for(unsigned int i = 0; i < nEtBins; ++i)
    std::cout << ">>> Et bin " << i << ":   [" << EtBinEdges->at(i) << "," << EtBinEdges->at(i+1) << "]" << std::endl;
  std::cout << std::endl;
  
  
  
  HtBinEdges = new std::vector<double>;
  for(unsigned int EtBinEdgeIt = 0; EtBinEdgeIt < EtBinEdges->size(); ++EtBinEdgeIt)
    HtBinEdges -> push_back( 2. * EtBinEdges->at(EtBinEdgeIt) );
  nHtBins = HtBinEdges->size()-1;
  
  for(unsigned int i = 0; i < nHtBins; ++i)
    std::cout << ">>> Ht bin " << i << ":   [" << HtBinEdges->at(i) << "," << HtBinEdges->at(i+1) << "]" << std::endl;
  std::cout << std::endl;
  
  
  
  
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
  TH1F* h_EoP_MC = new TH1F("h_EoP_MC","",480,0.,5.);
  h_EoP_MC -> Sumw2();
  TH1F* h_EoP_DA = new TH1F("h_EoP_DA","",480,0.,5.);
  h_EoP_DA -> Sumw2();
  
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
  
  std::vector<double>* EoP_EtBin_MC    = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_MC = new std::vector<double>[nEtBins];
  std::vector<double>* EoP_EtBin_DA                      = new std::vector<double>[nEtBins];
  std::vector<double>* EoP_EtBin_fit_DA                  = new std::vector<double>[nEtBins];
  std::vector<double>* EoP_EtBin_gausFit_DA              = new std::vector<double>[nEtBins];
  std::vector<double>* EoP_EtBin_mean_DA                 = new std::vector<double>[nEtBins];
  std::vector<double>* EoP_EtBin_recursiveMean_DA        = new std::vector<double>[nEtBins];
  std::vector<double>* EoP_EtBin_smallestInterval_DA     = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_DA                   = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_fit_DA               = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_gausFit_DA           = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_mean_DA              = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_recursiveMean_DA     = new std::vector<double>[nEtBins];
  std::vector<double>* weight_EtBin_smallestInterval_DA  = new std::vector<double>[nEtBins];
    
  TH1F** h_Zpt_EtBin_MC = new TH1F*[nEtBins];
  TH1F** h_Zpt_EtBin_DA = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_MC = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_DA = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_fit_DA = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_gausFit_DA = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_mean_DA = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_recursiveMean_DA = new TH1F*[nEtBins];
  TH1F** h_EoP_EtBin_smallestInterval_DA = new TH1F*[nEtBins];
  
  TF1** f_gausFit_EtBin_MC = new TF1*[nEtBins];
  TF1** f_gausFit_EtBin_DA = new TF1*[nEtBins];
  
  h_Et_EtBin_MC = new TH1F*[nEtBins];
  
  h_Et_EtBin_MC = new TH1F*[nEtBins];
  h_Et_EtBin_DA = new TH1F*[nEtBins];
  h_Et_EtBin_fit_DA = new TH1F*[nEtBins];
  h_Et_EtBin_gausFit_DA = new TH1F*[nEtBins];
  h_Et_EtBin_mean_DA = new TH1F*[nEtBins];
  h_Et_EtBin_recursiveMean_DA = new TH1F*[nEtBins];
  h_Et_EtBin_smallestInterval_DA = new TH1F*[nEtBins];
  
  for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
  {
    char histoName[50];
    
    sprintf(histoName,"h_EoP_EtBin%d_MC",EtBin);
    h_EoP_EtBin_MC[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
    //h_EoP_EtBin_MC[EtBin] -> SetFillColor(kGreen+2);
    //h_EoP_EtBin_MC[EtBin] -> SetFillStyle(3004);
    //h_EoP_EtBin_MC[EtBin] -> SetMarkerStyle(7);
    //h_EoP_EtBin_MC[EtBin] -> SetMarkerColor(kGreen+2);
    //h_EoP_EtBin_MC[EtBin] -> SetLineColor(kGreen+2);
    h_EoP_EtBin_MC[EtBin]->Sumw2();
    
    sprintf(histoName, "h_Zpt_EtBin%d_MC",EtBin);
    h_Zpt_EtBin_MC[EtBin] = new TH1F(histoName,"",300,0.,300.);
    h_Zpt_EtBin_MC[EtBin] -> SetLineColor(kGreen+2);
    h_Zpt_EtBin_MC[EtBin] -> Sumw2();
    
    sprintf(histoName, "Ht_EtBin%d_MC",EtBin);
    h_Et_EtBin_MC[EtBin] = new TH1F(histoName,"",5000,0.,1000.);
    h_Et_EtBin_MC[EtBin] -> SetLineColor(kGreen+2);
    h_Et_EtBin_MC[EtBin] -> Sumw2();
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
  
  int MCEntries = EoP_MC.size();
  for(int ientry = 0; ientry < MCEntries; ++ientry)
  {   
    if( (ientry%100000 == 0) ) std::cout << "reading   MC entry " << ientry << " / " << MCEntries << "\r" << std::flush;
    
    
    int EtBin = MyFindBin(Et_MC.at(ientry),EtBinEdges);
    if( EtBin == -1 ) continue;
    
    int HtBin = MyFindBin(Ht_MC.at(ientry),HtBinEdges);
    if( HtBin == -1 ) continue;
    
    EoP_EtBin_MC[EtBin].push_back( EoP_MC.at(ientry) );
    weight_EtBin_MC[EtBin].push_back( weight_MC.at(ientry) );
    
    h_Zpt_EtBin_MC[EtBin] -> Fill( Zpt_MC.at(ientry),weight_MC.at(ientry) );
    h_EoP_EtBin_MC[EtBin] -> Fill( EoP_MC.at(ientry),weight_MC.at(ientry) );
    h_Et_EtBin_MC[EtBin] -> Fill( Et_MC.at(ientry),weight_MC.at(ientry) );
    h_scEta_MC -> Fill( scEta_MC.at(ientry),weight_MC.at(ientry) );
    h_R9_MC -> Fill( R9_MC.at(ientry),weight_MC.at(ientry) );
    h_Ht_MC -> Fill( Ht_MC.at(ientry),weight_MC.at(ientry) );
    h_Zpt_MC -> Fill( Zpt_MC.at(ientry),weight_MC.at(ientry) );
    h_EoP_MC -> Fill( EoP_MC.at(ientry),weight_MC.at(ientry) );
  }
  std::cout << std::endl;
  
  
  
  for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
  {
    double x = h_Et_EtBin_MC[EtBin]->GetMean();
    
    scale_fit_DA -> SetPoint(EtBin,x,1.);
    scale_fit_MC -> SetPoint(EtBin,x,1.);
    scale_fit_DAOverMC -> SetPoint(EtBin,x,1.);
    
    scale_gausFit_DA -> SetPoint(EtBin,x,1.);
    scale_gausFit_MC -> SetPoint(EtBin,x,1.);
    scale_gausFit_DAOverMC -> SetPoint(EtBin,x,1.);
    
    scale_mean_DA -> SetPoint(EtBin,x,1.);
    scale_mean_MC -> SetPoint(EtBin,x,1.);
    scale_mean_DAOverMC -> SetPoint(EtBin,x,1.);
    
    scale_recursiveMean_DA -> SetPoint(EtBin,x,1.);
    scale_recursiveMean_MC -> SetPoint(EtBin,x,1.);
    scale_recursiveMean_DAOverMC -> SetPoint(EtBin,x,1.);
    
    scale_smallestInterval_DA -> SetPoint(EtBin,x,1.);
    scale_smallestInterval_MC -> SetPoint(EtBin,x,1.);
    scale_smallestInterval_DAOverMC -> SetPoint(EtBin,x,1.);
  }
  
  
  
  for(int step = 1; step < nSteps+1; ++step)
  {
    std::cout << std::endl;
    std::cout << "****** step " << step << " ******" << std::endl;
    
    
    TProfile* p_avgEtCorr_fit              = new TProfile("p_avgEtCorr_fit",             "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_gausFit          = new TProfile("p_avgEtCorr_gausFit",         "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_mean             = new TProfile("p_avgEtCorr_mean",            "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_recursiveMean    = new TProfile("p_avgEtCorr_recursiveMean",   "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_smallestInterval = new TProfile("p_avgEtCorr_smallestInterval","",nEtBins,30.,1000.);
    MySetBins(p_avgEtCorr_fit,*EtBinEdges);
    MySetBins(p_avgEtCorr_gausFit,*EtBinEdges);
    MySetBins(p_avgEtCorr_mean,*EtBinEdges);
    MySetBins(p_avgEtCorr_recursiveMean,*EtBinEdges);
    MySetBins(p_avgEtCorr_smallestInterval,*EtBinEdges);
    
    TProfile* p_avgEtCorr_rebin2_fit              = new TProfile("p_avgEtCorr_rebin2_fit",             "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_gausFit          = new TProfile("p_avgEtCorr_rebin2_gausFit",         "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_mean             = new TProfile("p_avgEtCorr_rebin2_mean",            "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_recursiveMean    = new TProfile("p_avgEtCorr_rebin2_recursiveMean",   "",nEtBins,30.,1000.);
    TProfile* p_avgEtCorr_rebin2_smallestInterval = new TProfile("p_avgEtCorr_rebin2_smallestInterval","",nEtBins,30.,1000.);
    MySetBins(p_avgEtCorr_rebin2_fit,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_gausFit,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_mean,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_recursiveMean,*EtBinEdges);
    MySetBins(p_avgEtCorr_rebin2_smallestInterval,*EtBinEdges);
    p_avgEtCorr_rebin2_fit              -> Rebin(2);
    p_avgEtCorr_rebin2_gausFit          -> Rebin(2);
    p_avgEtCorr_rebin2_mean             -> Rebin(2);
    p_avgEtCorr_rebin2_recursiveMean    -> Rebin(2);
    p_avgEtCorr_rebin2_smallestInterval -> Rebin(2);
    
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
        EoP_EtBin_DA[EtBin].clear();
        EoP_EtBin_fit_DA[EtBin].clear();
        EoP_EtBin_gausFit_DA[EtBin].clear();
        EoP_EtBin_mean_DA[EtBin].clear();
        EoP_EtBin_recursiveMean_DA[EtBin].clear();
        EoP_EtBin_smallestInterval_DA[EtBin].clear();
        weight_EtBin_DA[EtBin].clear();
        weight_EtBin_fit_DA[EtBin].clear();
        weight_EtBin_gausFit_DA[EtBin].clear();
        weight_EtBin_mean_DA[EtBin].clear();
        weight_EtBin_recursiveMean_DA[EtBin].clear();
        weight_EtBin_smallestInterval_DA[EtBin].clear();
    }	
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      char histoName[50];
      
      sprintf(histoName,"h_Zpt_EtBin%d_DA_step%d",EtBin,step);
      h_Zpt_EtBin_DA[EtBin] = new TH1F(histoName,"",300,0.,300.);
      h_Zpt_EtBin_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_EoP_EtBin%d_DA_step%d",EtBin,step);
      h_EoP_EtBin_DA[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
      h_EoP_EtBin_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_EoP_EtBin%d_fit_DA_step%d",EtBin,step);
      h_EoP_EtBin_fit_DA[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
      h_EoP_EtBin_fit_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_EoP_EtBin%d_gausFit_DA_step%d",EtBin,step);
      h_EoP_EtBin_gausFit_DA[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
      h_EoP_EtBin_gausFit_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_EoP_EtBin%d_mean_DA_step%d",EtBin,step);
      h_EoP_EtBin_mean_DA[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
      h_EoP_EtBin_mean_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_EoP_EtBin%d_recursiveMean_DA_step%d",EtBin,step);
      h_EoP_EtBin_recursiveMean_DA[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
      h_EoP_EtBin_recursiveMean_DA[EtBin] -> Sumw2();
      
      sprintf(histoName,"h_EoP_EtBin%d_smallestInterval_DA_step%d",EtBin,step);
      h_EoP_EtBin_smallestInterval_DA[EtBin] = new TH1F(histoName,"",nBinsEoP,EoPMin,EoPMax);
      h_EoP_EtBin_smallestInterval_DA[EtBin] -> Sumw2();
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
    
    
    int DAEntries = EoP_DA.size();
    for(int ientry = 0; ientry < DAEntries; ++ientry)
    {   
      if( (ientry%100000 == 0) ) std::cout << "reading DATA entry " << ientry << " / " << DAEntries << "\r" << std::flush;
      
      
      double k = 1.;
      double Et = Et_DA.at(ientry)/k;
      int EtBin = MyFindBin(Et,EtBinEdges);
      double EoP = EoP_DA.at(ientry)/k;
      if( EtBin != -1 )
      {
        EoP_EtBin_DA[EtBin].push_back( EoP );
        weight_EtBin_DA[EtBin].push_back( weight_DA.at(ientry) );
        
        
        h_Zpt_EtBin_DA[EtBin] -> Fill( Zpt_DA.at(ientry),weight_DA.at(ientry) );
        h_EoP_EtBin_DA[EtBin] -> Fill( EoP,weight_DA.at(ientry) );
        
        h_Et_EtBin_DA[EtBin] -> Fill( Et,weight_DA.at(ientry) );
      }
      
      
      k = MyEval(scale_fit_DAOverMC,Et_fit_DA.at(ientry));
      Et = Et_fit_DA.at(ientry)/k;
      Et_fit_DA.at(ientry) = Et;
      p_avgEtCorr_fit -> Fill(Et,Et/Et_DA.at(ientry));
      p_avgEtCorr_rebin2_fit -> Fill(Et,Et/Et_DA.at(ientry));
      EtBin = MyFindBin(Et,EtBinEdges);
      EoP = EoP_fit_DA.at(ientry)/k;
      EoP_fit_DA.at(ientry) = EoP;
      if( EtBin != -1 )
      {
        EoP_EtBin_fit_DA[EtBin].push_back( EoP );
        weight_EtBin_fit_DA[EtBin].push_back( weight_DA.at(ientry) );
        
        h_EoP_EtBin_fit_DA[EtBin] -> Fill( EoP,weight_DA.at(ientry) );
        
        h_Et_EtBin_fit_DA[EtBin] -> Fill( Et,weight_DA.at(ientry) );
      }
      
      k = MyEval(scale_mean_DAOverMC,Et_mean_DA.at(ientry));
      Et = Et_mean_DA.at(ientry)/k;
      Et_mean_DA.at(ientry) = Et;
      p_avgEtCorr_mean -> Fill(Et,Et/Et_DA.at(ientry));
      p_avgEtCorr_rebin2_mean -> Fill(Et,Et/Et_DA.at(ientry));
      EtBin = MyFindBin(Et,EtBinEdges);
      EoP = EoP_mean_DA.at(ientry)/k;
      EoP_mean_DA.at(ientry) = EoP;
      if( EtBin != -1 )
      {
        EoP_EtBin_mean_DA[EtBin].push_back( EoP );
        weight_EtBin_mean_DA[EtBin].push_back( weight_DA.at(ientry) );
        
        h_EoP_EtBin_mean_DA[EtBin] -> Fill( EoP,weight_DA.at(ientry) );
        
        h_Et_EtBin_mean_DA[EtBin] -> Fill( Et,weight_DA.at(ientry) );
      }
      
      k = MyEval(scale_gausFit_DAOverMC,Et_gausFit_DA.at(ientry));
      Et = Et_gausFit_DA.at(ientry)/k;
      Et_gausFit_DA.at(ientry) = Et;
      p_avgEtCorr_gausFit -> Fill(Et,Et/Et_DA.at(ientry));
      p_avgEtCorr_rebin2_gausFit -> Fill(Et,Et/Et_DA.at(ientry));
      EtBin = MyFindBin(Et,EtBinEdges);
      EoP = EoP_gausFit_DA.at(ientry)/k;
      EoP_gausFit_DA.at(ientry) = EoP;
      if( EtBin != -1 )
      {
        EoP_EtBin_gausFit_DA[EtBin].push_back( EoP );
        weight_EtBin_gausFit_DA[EtBin].push_back( weight_DA.at(ientry) );
        
        h_EoP_EtBin_gausFit_DA[EtBin] -> Fill( EoP,weight_DA.at(ientry) );
        
        h_Et_EtBin_gausFit_DA[EtBin] -> Fill( Et,weight_DA.at(ientry) );
      }
      
      k = MyEval(scale_recursiveMean_DAOverMC,Et_recursiveMean_DA.at(ientry));
      Et = Et_recursiveMean_DA.at(ientry)/k;
      Et_recursiveMean_DA.at(ientry) = Et;
      p_avgEtCorr_recursiveMean -> Fill(Et,Et/Et_DA.at(ientry));
      p_avgEtCorr_rebin2_recursiveMean -> Fill(Et,Et/Et_DA.at(ientry));
      EtBin = MyFindBin(Et,EtBinEdges);
      EoP = EoP_recursiveMean_DA.at(ientry)/k;
      EoP_recursiveMean_DA.at(ientry) = EoP;
      if( EtBin != -1 )
      {
        EoP_EtBin_recursiveMean_DA[EtBin].push_back( EoP );
        weight_EtBin_recursiveMean_DA[EtBin].push_back( weight_DA.at(ientry) );
        
        h_EoP_EtBin_recursiveMean_DA[EtBin] -> Fill( EoP,weight_DA.at(ientry) );
        
        h_Et_EtBin_recursiveMean_DA[EtBin] -> Fill( Et,weight_DA.at(ientry) );
      }
      
      k = MyEval(scale_smallestInterval_DAOverMC,Et_smallestInterval_DA.at(ientry));
      Et = Et_smallestInterval_DA.at(ientry)/k;
      Et_smallestInterval_DA.at(ientry) = Et;
      p_avgEtCorr_smallestInterval -> Fill(Et,Et/Et_DA.at(ientry));
      p_avgEtCorr_rebin2_smallestInterval -> Fill(Et,Et/Et_DA.at(ientry));
      EtBin = MyFindBin(Et,EtBinEdges);
      EoP = EoP_smallestInterval_DA.at(ientry)/k;
      EoP_smallestInterval_DA.at(ientry) = EoP;
      if( EtBin != -1 )
      {
        EoP_EtBin_smallestInterval_DA[EtBin].push_back( EoP );
        weight_EtBin_smallestInterval_DA[EtBin].push_back( weight_DA.at(ientry) );
        
        h_EoP_EtBin_smallestInterval_DA[EtBin] -> Fill( EoP,weight_DA.at(ientry) );
        
        h_Et_EtBin_smallestInterval_DA[EtBin] -> Fill( Et,weight_DA.at(ientry) );
      }
      
      
      
      if( step == 1 )
      {
        h_scEta_DA -> Fill( scEta_DA.at(ientry),weight_DA.at(ientry) );
        h_R9_DA -> Fill( R9_DA.at(ientry),weight_DA.at(ientry) );
        h_Ht_DA -> Fill( Ht_DA.at(ientry),weight_DA.at(ientry) );
        h_Zpt_DA -> Fill( Zpt_DA.at(ientry),weight_DA.at(ientry) );
        h_EoP_DA -> Fill( EoP_DA.at(ientry),weight_DA.at(ientry) );
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
    
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      std::cout << ">>> EtBin::" << EtBin << std::endl;
      std::cout << ">>>>>> nEvents:                  " << EoP_EtBin_DA[EtBin].size() << std::endl;
      
      double x = h_Et_EtBin_MC[EtBin]->GetMean();
      double ex = h_Et_EtBin_MC[EtBin]->GetMeanError();
      double exlow = ex;
      double exhig = ex;
      
      h_Zpt_EtBin_MC[EtBin] -> Scale(1./h_Zpt_EtBin_MC[EtBin]->Integral());
      h_Zpt_EtBin_DA[EtBin] -> Scale(1./h_Zpt_EtBin_DA[EtBin]->Integral());
      
      h_EoP_EtBin_MC[EtBin] -> Scale(1./h_EoP_EtBin_MC[EtBin]->Integral());
      h_EoP_EtBin_DA[EtBin] -> Scale(1./h_EoP_EtBin_DA[EtBin]->Integral());
      h_EoP_EtBin_fit_DA[EtBin] -> Scale(1./h_EoP_EtBin_fit_DA[EtBin]->Integral());
      h_EoP_EtBin_gausFit_DA[EtBin] -> Scale(1./h_EoP_EtBin_gausFit_DA[EtBin]->Integral());
      h_EoP_EtBin_mean_DA[EtBin] -> Scale(1./h_EoP_EtBin_mean_DA[EtBin]->Integral());
      h_EoP_EtBin_recursiveMean_DA[EtBin] -> Scale(1./h_EoP_EtBin_recursiveMean_DA[EtBin]->Integral());
      h_EoP_EtBin_smallestInterval_DA[EtBin] -> Scale(1./h_EoP_EtBin_smallestInterval_DA[EtBin]->Integral());
      
      
      // fit
      std::cout << ">>>>>> nEvents_fit:              " << EoP_EtBin_fit_DA[EtBin].size() << std::endl;
      if( EoP_EtBin_fit_DA[EtBin].size() > 3 )
      {
        double scale_MC = 1.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindTemplateFit(scale_DA,scaleErr_DA,h_EoP_EtBin_MC[EtBin],h_EoP_EtBin_fit_DA[EtBin]);
        double y = scale_DA / scale_MC;
	double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
	double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        
        scale_fit_MC -> SetPoint(EtBin,x,scale_MC);
        scale_fit_MC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_fit_DA -> SetPoint(EtBin,x,scale_DA);
        scale_fit_DA -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_fit_DAOverMC -> SetPoint(EtBin,x,y);
        scale_fit_DAOverMC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // gausFit
      std::cout << ">>>>>> nEvents_gausFit:          " << EoP_EtBin_gausFit_DA[EtBin].size() << std::endl;
      if( EoP_EtBin_gausFit_DA[EtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        char funcName_MC[50];
        sprintf(funcName_MC,"f_gausFit_EtBin%d_step%d_MC",EtBin,step);
        char funcName_DA[50];
        sprintf(funcName_DA,"f_gausFit_EtBin%d_step%d_DA",EtBin,step);
        FindGausFit(scale_MC,scaleErr_MC,EoP_EtBin_MC[EtBin],weight_EtBin_MC[EtBin],nBinsEoP,EoPMin,EoPMax,&(f_gausFit_EtBin_MC[EtBin]),std::string(funcName_MC),1.);
        FindGausFit(scale_DA,scaleErr_DA,EoP_EtBin_gausFit_DA[EtBin],weight_EtBin_gausFit_DA[EtBin],nBinsEoP,EoPMin,EoPMax,&(f_gausFit_EtBin_DA[EtBin]),std::string(funcName_DA),1.);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_gausFit_MC -> SetPoint(EtBin,x,scale_MC);
        scale_gausFit_MC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_gausFit_DA -> SetPoint(EtBin,x,scale_DA);
        scale_gausFit_DA -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_gausFit_DAOverMC -> SetPoint(EtBin,x,y);
        scale_gausFit_DAOverMC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // mean
      std::cout << ">>>>>> nEvents_mean:             " << EoP_EtBin_mean_DA[EtBin].size() << std::endl;
      if( EoP_EtBin_mean_DA[EtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindMean(scale_MC,scaleErr_MC,EoP_EtBin_MC[EtBin],weight_EtBin_MC[EtBin]);
        FindMean(scale_DA,scaleErr_DA,EoP_EtBin_mean_DA[EtBin],weight_EtBin_mean_DA[EtBin]);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_mean_MC -> SetPoint(EtBin,x,scale_MC);
        scale_mean_MC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_mean_DA -> SetPoint(EtBin,x,scale_DA);
        scale_mean_DA -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_mean_DAOverMC -> SetPoint(EtBin,x,y);
        scale_mean_DAOverMC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // recursive mean
      std::cout << ">>>>>> nEvents_recursiveMean:    " << EoP_EtBin_recursiveMean_DA[EtBin].size() << std::endl;
      if( EoP_EtBin_recursiveMean_DA[EtBin].size() > 3 )
      {
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindRecursiveMean(scale_MC,scaleErr_MC,EoP_EtBin_MC[EtBin],weight_EtBin_MC[EtBin],0.1,0.001);
        FindRecursiveMean(scale_DA,scaleErr_DA,EoP_EtBin_recursiveMean_DA[EtBin],weight_EtBin_recursiveMean_DA[EtBin],0.1,0.001);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_recursiveMean_MC -> SetPoint(EtBin,x,scale_MC);
        scale_recursiveMean_MC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_recursiveMean_DA -> SetPoint(EtBin,x,scale_DA);
        scale_recursiveMean_DA -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_recursiveMean_DAOverMC -> SetPoint(EtBin,x,y);
        scale_recursiveMean_DAOverMC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
      }
      
      
      // smallest interval
      std::cout << ">>>>>> nEvents_smallestInterval: " << EoP_EtBin_smallestInterval_DA[EtBin].size() << std::endl;
      if( EoP_EtBin_smallestInterval_DA[EtBin].size() > 3 )
      {
        double min_MC;
        double max_MC;
        double min_DA;
        double max_DA;
        double scale_MC = 0.;
        double scale_DA = 0.;
        double scaleErr_MC = 0.;
        double scaleErr_DA = 0.;
        FindSmallestInterval(scale_MC,scaleErr_MC,min_MC,max_MC,0.,5.,0.0005,EoP_EtBin_MC[EtBin],weight_EtBin_MC[EtBin],0.5);
        FindSmallestInterval(scale_DA,scaleErr_DA,min_DA,max_DA,0.,5.,0.0005,EoP_EtBin_smallestInterval_DA[EtBin],weight_EtBin_smallestInterval_DA[EtBin],0.5);
        smallestIntervalMins_MC.push_back(min_MC);
        smallestIntervalMaxs_MC.push_back(max_MC);
        smallestIntervalMins_DA.push_back(min_DA);
        smallestIntervalMaxs_DA.push_back(max_DA);
        double y = scale_DA / scale_MC;
        double eylow = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
        double eyhig = y * sqrt( pow(scaleErr_DA/scale_DA,2) + pow(scaleErr_MC/scale_MC,2) );
                
        scale_smallestInterval_MC -> SetPoint(EtBin,x,scale_MC);
        scale_smallestInterval_MC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_smallestInterval_DA -> SetPoint(EtBin,x,scale_DA);
        scale_smallestInterval_DA -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
        scale_smallestInterval_DAOverMC -> SetPoint(EtBin,x,y);
        scale_smallestInterval_DAOverMC -> SetPointError(EtBin,exlow,exhig,eylow,eyhig);
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
    p_avgEtCorr_rebin2_fit -> Write();
    
    scale_gausFit_MC           -> Write("scale_gausFit_MC");
    scale_gausFit_DA           -> Write("scale_gausFit_DA");
    scale_gausFit_DAOverMC     -> Write("scale_gausFit_DAOverMC");
    p_avgEtCorr_gausFit        -> Write();
    p_avgEtCorr_rebin2_gausFit -> Write();
    
    scale_mean_MC           -> Write("scale_mean_MC");
    scale_mean_DA           -> Write("scale_mean_DA");
    scale_mean_DAOverMC     -> Write("scale_mean_DAOverMC");
    p_avgEtCorr_mean        -> Write();
    p_avgEtCorr_rebin2_mean -> Write();
    
    scale_recursiveMean_MC           -> Write("scale_recursiveMean_MC");
    scale_recursiveMean_DA           -> Write("scale_recursiveMean_DA");
    scale_recursiveMean_DAOverMC     -> Write("scale_recursiveMean_DAOverMC");
    p_avgEtCorr_recursiveMean        -> Write();
    p_avgEtCorr_rebin2_recursiveMean -> Write();
    
    scale_smallestInterval_MC           -> Write("scale_smallestInterval_MC");
    scale_smallestInterval_DA           -> Write("scale_smallestInterval_DA");
    scale_smallestInterval_DAOverMC     -> Write("scale_smallestInterval_DAOverMC");
    p_avgEtCorr_smallestInterval        -> Write();
    p_avgEtCorr_rebin2_smallestInterval -> Write();
    
    
    delete p_avgEtCorr_fit;
    delete p_avgEtCorr_gausFit;
    delete p_avgEtCorr_mean;
    delete p_avgEtCorr_recursiveMean;
    delete p_avgEtCorr_smallestInterval;
    delete p_avgEtCorr_rebin2_fit;
    delete p_avgEtCorr_rebin2_gausFit;
    delete p_avgEtCorr_rebin2_mean;
    delete p_avgEtCorr_rebin2_recursiveMean;
    delete p_avgEtCorr_rebin2_smallestInterval;
    
    
    baseDir -> cd();
    subDir = baseDir -> mkdir("EoP_EtBin");
    subDir -> cd();
    
    std::string outputPdf_DAMC  = plotFolderName + "h_EoP_EtBin_DAOverMC_" + catType+"_cat"+ Form("%d",category) + ".pdf";
    std::string outputPdf2_DAMC = plotFolderName + "h_Zpt_EtBin_DAOverMC_" + catType+"_cat"+ Form("%d",category) + ".pdf";
    
    std::string outputPdf_MC = plotFolderName + "h_EoP_EtBin_MC_"   + catType+"_cat"+ Form("%d",category) + ".pdf";
    std::string outputPdf_DA = plotFolderName + "h_EoP_EtBin_DATA_" + catType+"_cat"+ Form("%d",category) + ".pdf";
    
    for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
    {
      h_EoP_EtBin_MC[EtBin] -> Write();
      h_EoP_EtBin_DA[EtBin] -> Write();
      h_EoP_EtBin_fit_DA[EtBin] -> Write();
      h_EoP_EtBin_gausFit_DA[EtBin] -> Write();
      h_EoP_EtBin_mean_DA[EtBin] -> Write();
      h_EoP_EtBin_recursiveMean_DA[EtBin] -> Write();
      h_EoP_EtBin_smallestInterval_DA[EtBin] -> Write();
      
      f_gausFit_EtBin_MC[EtBin] -> Write();
      f_gausFit_EtBin_DA[EtBin] -> Write();
      
      
      
      if( step == 1 )
      {
        TCanvas* c_DAOverMC = new TCanvas("c_Zpt");
        c_DAOverMC -> cd();
        c_DAOverMC -> SetGridx();
        c_DAOverMC -> SetGridy();
        
        char axisTitle[50];
        sprintf(axisTitle,"p_{T}(ee) [GeV]   -   E_{T} #in [%d,%d]",int(EtBinEdges->at(EtBin)),int(EtBinEdges->at(EtBin+1)));
        h_Zpt_EtBin_MC[EtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_Zpt_EtBin_MC[EtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_Zpt_EtBin_MC[EtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_Zpt_EtBin_MC[EtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_Zpt_EtBin_MC[EtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_Zpt_EtBin_MC[EtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_Zpt_EtBin_MC[EtBin] -> SetLineColor(kBlack);
        h_Zpt_EtBin_MC[EtBin] -> SetLineWidth(1);
        h_Zpt_EtBin_MC[EtBin] -> GetXaxis() -> SetRangeUser(0.,300.);
        float maximum = std::max(h_Zpt_EtBin_MC[EtBin]->GetMaximum(),h_Zpt_EtBin_DA[EtBin]->GetMaximum());
        h_Zpt_EtBin_MC[EtBin] -> SetMaximum( 1.1*maximum );
        
        h_Zpt_EtBin_MC[EtBin] -> Draw("hist");
        h_Zpt_EtBin_DA[EtBin] -> Draw("P,same");
        
        if( EtBin == 0 )         c_DAOverMC -> Print((outputPdf2_DAMC+"[").c_str());
        c_DAOverMC -> Print(outputPdf2_DAMC.c_str());
        if( EtBin == nEtBins-1 ) c_DAOverMC -> Print((outputPdf2_DAMC+"]").c_str());
        
        
        
        c_DAOverMC = new TCanvas("c_EoP");
        c_DAOverMC -> cd();
        c_DAOverMC -> SetGridx();
        c_DAOverMC -> SetGridy();
        
        sprintf(axisTitle,"E/p scale   -   E_{T} #in [%d,%d]",int(EtBinEdges->at(EtBin)),int(EtBinEdges->at(EtBin+1)));
        h_EoP_EtBin_MC[EtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_EoP_EtBin_MC[EtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_EoP_EtBin_MC[EtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_EoP_EtBin_MC[EtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_EoP_EtBin_MC[EtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_EoP_EtBin_MC[EtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_EoP_EtBin_MC[EtBin] -> SetLineColor(kBlack);
        h_EoP_EtBin_MC[EtBin] -> SetLineWidth(1);
        h_EoP_EtBin_MC[EtBin] -> GetXaxis() -> SetRangeUser(0.65,1.34999);
        maximum = std::max(h_EoP_EtBin_MC[EtBin]->GetMaximum(),h_EoP_EtBin_DA[EtBin]->GetMaximum());
        h_EoP_EtBin_MC[EtBin] -> SetMaximum( 1.1*maximum );
        
        h_EoP_EtBin_MC[EtBin] -> Draw("hist");
        h_EoP_EtBin_DA[EtBin] -> Draw("P,same");
        
        if( EtBin == 0 )         c_DAOverMC -> Print((outputPdf_DAMC+"[").c_str());
        c_DAOverMC -> Print(outputPdf_DAMC.c_str());
        if( EtBin == nEtBins-1 ) c_DAOverMC -> Print((outputPdf_DAMC+"]").c_str());
        
        
        
        TCanvas* c_MC = new TCanvas();
        c_MC -> cd();
        c_MC -> SetGridx();
        c_MC -> SetGridy();
        
        double x,y;
        
        scale_mean_MC -> GetPoint(EtBin,x,y);
        TArrow* line_mean_MC = new TArrow(y,0.,y,h_EoP_EtBin_MC[EtBin]->GetBinContent(h_EoP_EtBin_MC[EtBin]->FindBin(y)));
        line_mean_MC -> SetLineColor(kBlack);
        line_mean_MC -> SetLineWidth(2);
        
        scale_recursiveMean_MC -> GetPoint(EtBin,x,y);
        TArrow* line_recursiveMean_MC = new TArrow(y,0.,y,h_EoP_EtBin_MC[EtBin]->GetMaximum());
        line_recursiveMean_MC -> SetLineColor(kBlue);
        line_recursiveMean_MC -> SetLineWidth(2);
        TArrow* line_recursiveMean_min_MC = new TArrow(y-0.05,0.,y-0.05,h_EoP_EtBin_MC[EtBin]->GetMaximum());
        line_recursiveMean_min_MC -> SetLineColor(kBlue);
        line_recursiveMean_min_MC -> SetLineWidth(2);
        TArrow* line_recursiveMean_max_MC = new TArrow(y+0.05,0.,y+0.05,h_EoP_EtBin_MC[EtBin]->GetMaximum());
        line_recursiveMean_max_MC -> SetLineColor(kBlue);
        line_recursiveMean_max_MC -> SetLineWidth(2);
        
        scale_smallestInterval_MC -> GetPoint(EtBin,x,y);
        TArrow* line_smallestInterval_MC = new TArrow(y,0.,y,h_EoP_EtBin_MC[EtBin]->GetBinContent(h_EoP_EtBin_MC[EtBin]->FindBin(y)));
        line_smallestInterval_MC -> SetLineColor(kRed);
        line_smallestInterval_MC -> SetLineWidth(2);
        
        h_EoP_EtBin_MC[EtBin] -> Draw("HIST");
        
        TH1F* clone = (TH1F*)( h_EoP_EtBin_MC[EtBin]->Clone() );
        for(int bin = 1; bin < clone->GetNbinsX(); ++bin)
          if( (clone->GetBinCenter(bin) < smallestIntervalMins_MC.at(EtBin)) ||
              (clone->GetBinCenter(bin) > smallestIntervalMaxs_MC.at(EtBin)) )
            clone -> SetBinContent(bin,0.);
        clone -> SetFillColor(kYellow);
        clone -> SetLineWidth(0);
        clone -> Draw("HIST,same");
        
        line_mean_MC              -> Draw("same");
        line_recursiveMean_MC     -> Draw("same");
        line_recursiveMean_min_MC -> Draw("same");
        line_recursiveMean_max_MC -> Draw("same");
        line_smallestInterval_MC  -> Draw("same");
        f_gausFit_EtBin_MC[EtBin] -> Draw("same");
        
        
        TCanvas* c_DA = new TCanvas();
        c_DA -> cd();
        c_DA -> SetGridx();
        c_DA -> SetGridy();
        
        sprintf(axisTitle,"E/p scale   -   E_{T} #in [%d,%d]",int(EtBinEdges->at(EtBin)),int(EtBinEdges->at(EtBin+1)));
        h_EoP_EtBin_DA[EtBin] -> GetXaxis() -> SetTitle(axisTitle);
        sprintf(axisTitle,"event fraction");
        h_EoP_EtBin_DA[EtBin] -> GetYaxis() -> SetTitle(axisTitle);
        h_EoP_EtBin_DA[EtBin] -> GetXaxis() -> SetLabelSize(0.04);
        h_EoP_EtBin_DA[EtBin] -> GetYaxis() -> SetLabelSize(0.04);
        h_EoP_EtBin_DA[EtBin] -> GetXaxis() -> SetTitleSize(0.05);
        h_EoP_EtBin_DA[EtBin] -> GetYaxis() -> SetTitleSize(0.05);
        h_EoP_EtBin_DA[EtBin] -> SetLineColor(kBlack);
        h_EoP_EtBin_DA[EtBin] -> SetLineWidth(1);
        h_EoP_EtBin_DA[EtBin] -> GetXaxis() -> SetRangeUser(0.65,1.34999);
        h_EoP_EtBin_DA[EtBin] -> SetMaximum( 1.1*h_EoP_EtBin_DA[EtBin]->GetMaximum() );
        
        
        scale_mean_DA -> GetPoint(EtBin,x,y);
        TArrow* line_mean_DA = new TArrow(y,0.,y,h_EoP_EtBin_DA[EtBin]->GetBinContent(h_EoP_EtBin_DA[EtBin]->FindBin(y)));
        line_mean_DA -> SetLineColor(kBlack);
        line_mean_DA -> SetLineWidth(2);
        
        scale_recursiveMean_DA -> GetPoint(EtBin,x,y);
        TArrow* line_recursiveMean_DA = new TArrow(y,0.,y,h_EoP_EtBin_DA[EtBin]->GetMaximum());
        line_recursiveMean_DA -> SetLineColor(kBlue);
        line_recursiveMean_DA -> SetLineWidth(2);
        TArrow* line_recursiveMean_min_DA = new TArrow(y-0.05,0.,y-0.05,h_EoP_EtBin_DA[EtBin]->GetMaximum());
        line_recursiveMean_min_DA -> SetLineColor(kBlue);
        line_recursiveMean_min_DA -> SetLineWidth(2);
        TArrow* line_recursiveMean_max_DA = new TArrow(y+0.05,0.,y+0.05,h_EoP_EtBin_DA[EtBin]->GetMaximum());
        line_recursiveMean_max_DA -> SetLineColor(kBlue);
        line_recursiveMean_max_DA -> SetLineWidth(2);
        
        scale_smallestInterval_DA -> GetPoint(EtBin,x,y);
        TArrow* line_smallestInterval_DA = new TArrow(y,0.,y,h_EoP_EtBin_DA[EtBin]->GetBinContent(h_EoP_EtBin_DA[EtBin]->FindBin(y)));
        line_smallestInterval_DA -> SetLineColor(kRed);
        line_smallestInterval_DA -> SetLineWidth(2);
        
        h_EoP_EtBin_DA[EtBin] -> Draw("HIST");
        
        clone = (TH1F*)( h_EoP_EtBin_DA[EtBin]->Clone() );
        for(int bin = 1; bin < clone->GetNbinsX(); ++bin)
          if( (clone->GetBinCenter(bin) < smallestIntervalMins_DA.at(EtBin)) ||
              (clone->GetBinCenter(bin) > smallestIntervalMaxs_DA.at(EtBin)) )
            clone -> SetBinContent(bin,0.);
        clone -> SetFillColor(kYellow);
        clone -> SetLineWidth(0);
        clone -> Draw("HIST,same");
        
        line_mean_DA              -> Draw("same");
        line_recursiveMean_DA     -> Draw("same");
        line_recursiveMean_min_DA -> Draw("same");
        line_recursiveMean_max_DA -> Draw("same");
        line_smallestInterval_DA  -> Draw("same");
        f_gausFit_EtBin_DA[EtBin] -> Draw("same");
        
        
        if( EtBin == 0 )
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
        if( EtBin == nEtBins-1 )
        {
          c_MC -> RedrawAxis();
          c_MC -> Print((outputPdf_MC+"]").c_str());
          c_DA -> RedrawAxis();
          c_DA -> Print((outputPdf_DA+"]").c_str());
        }
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
  
  
  
  outFile -> mkdir("Et_EtBin");
  outFile -> cd("Et_EtBin");
  
  for(unsigned int EtBin = 0; EtBin < nEtBins; ++EtBin)
  {
    h_Et_EtBin_MC[EtBin] -> Scale(1./h_Et_EtBin_MC[EtBin]->Integral());
    h_Et_EtBin_MC[EtBin] -> Write();
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
  
  h_EoP_MC -> Scale(1./h_EoP_MC->Integral());
  h_EoP_MC -> Write();
  h_EoP_DA -> Scale(1./h_EoP_DA->Integral());
  h_EoP_DA -> Write();
  
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

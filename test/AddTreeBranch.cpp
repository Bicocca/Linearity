//g++ -Wall -o AddTreeBranch `root-config --cflags --glibs` AddTreeBranch.cpp

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TChain.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TRandom.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>


int main()
{
  
  //------------
  // input Globe trees
  TFile* f = new TFile("../data/zeeTOT_tree.root", "read");
  TTree* treeG = (TTree*)f->Get("opttree");

  Float_t Gitype;  
  Float_t Grun;  
  Float_t Glumis;  
  Double_t Gevent;
  Float_t Gnvtx, Gpu_n;  
  Float_t Gfull_weight;
  Float_t Gfull_cat;
  Float_t Gdipho_mva_cat;
  Float_t Gmass;
  Float_t Get1, Get2, Geta1, Geta2, Gr91, Gr92;

  treeG->SetBranchStatus("*",0);
  treeG->SetBranchStatus("itype", 1);    treeG->SetBranchAddress("itype", &Gitype);
  treeG->SetBranchStatus("run", 1);    treeG->SetBranchAddress("run", &Grun);
  treeG->SetBranchStatus("lumis", 1);    treeG->SetBranchAddress("lumis", &Glumis);
  treeG->SetBranchStatus("event", 1);    treeG->SetBranchAddress("event", &Gevent);
  treeG->SetBranchStatus("nvtx", 1);    treeG->SetBranchAddress("nvtx", &Gnvtx);
  treeG->SetBranchStatus("pu_n", 1);    treeG->SetBranchAddress("pu_n", &Gpu_n);
  treeG->SetBranchStatus("full_weight", 1);    treeG->SetBranchAddress("full_weight", &Gfull_weight);
  treeG->SetBranchStatus("full_cat", 1);       treeG->SetBranchAddress("full_cat", &Gfull_cat);
  treeG->SetBranchStatus("dipho_mva_cat", 1);    treeG->SetBranchAddress("dipho_mva_cat", &Gdipho_mva_cat);
  treeG->SetBranchStatus("mass", 1);           treeG->SetBranchAddress("mass", &Gmass);
  treeG->SetBranchStatus("et1", 1);             treeG->SetBranchAddress("et1", &Get1);
  treeG->SetBranchStatus("et2", 1);             treeG->SetBranchAddress("et2", &Get2);
//   treeG->SetBranchStatus("sceta1", 1);           treeG->SetBranchAddress("sceta1", &Geta1);
//   treeG->SetBranchStatus("sceta2", 1);           treeG->SetBranchAddress("sceta2", &Geta2);
  treeG->SetBranchStatus("eta1", 1);           treeG->SetBranchAddress("eta1", &Geta1);
  treeG->SetBranchStatus("eta2", 1);           treeG->SetBranchAddress("eta2", &Geta2);
  treeG->SetBranchStatus("r91", 1);           treeG->SetBranchAddress("r91", &Gr91);
  treeG->SetBranchStatus("r92", 1);           treeG->SetBranchAddress("r92", &Gr92);


  //Shervin
  std::cout << " Shervin's ntuples " << std::endl;

  std::string inputFilesMC;
  inputFilesMC = "/gwterax2/users/martelli/CALIBRATION/MC/DYToEE_M20_powheg-Summer12-START53-ZSkim-runDependent-allRange.root";
  TChain* ntu_MC = new TChain("selected");
  ntu_MC->Add(inputFilesMC.c_str());
  std::cout << ">>>   MC: " << std::setw(8) << ntu_MC->GetEntries() << " entries" << std::endl;

  std::string inputFilesDA;
  inputFilesDA = "/gwterax2/users/martelli/CALIBRATION/data/DoubleElectron-ZSkim-RUN2012ABCD-22Jan-v1-allRange.root";
  TChain* ntu_DA = new TChain("selected");
  ntu_DA->Add(inputFilesDA.c_str());
  std::cout << ">>>   DA: " << std::setw(8) << ntu_DA->GetEntries() << " entries" << std::endl;


  // global variables 
  int runId,nPU,cat, nVtx, eventNumber, lumiBlock;
  float weight;
  bool HLTfire;

  ntu_MC -> SetBranchStatus("*",0);
  ntu_MC -> SetBranchStatus("HLTfire",  1);   ntu_MC -> SetBranchAddress("HLTfire",&HLTfire);
  ntu_MC -> SetBranchStatus("runNumber",1);   ntu_MC -> SetBranchAddress("runNumber",&runId);
  ntu_MC -> SetBranchStatus("eventNumber",1);   ntu_MC -> SetBranchAddress("eventNumber",&eventNumber);
  ntu_MC -> SetBranchStatus("lumiBlock",1);   ntu_MC -> SetBranchAddress("lumiBlock",&lumiBlock);
  ntu_MC -> SetBranchStatus("nPU",      1);   ntu_MC -> SetBranchAddress("nPU",&nPU);
  ntu_MC -> SetBranchStatus("nPV",      1);   ntu_MC -> SetBranchAddress("nPV",&nVtx);


  ntu_DA -> SetBranchStatus("*",0);
  ntu_DA -> SetBranchStatus("HLTfire",  1);   ntu_DA -> SetBranchAddress("HLTfire",&HLTfire);
  ntu_DA -> SetBranchStatus("runNumber",1);   ntu_DA -> SetBranchAddress("runNumber",&runId);
  ntu_DA -> SetBranchStatus("eventNumber",1); ntu_DA -> SetBranchAddress("eventNumber",&eventNumber);
  ntu_DA -> SetBranchStatus("lumiBlock",1);   ntu_DA -> SetBranchAddress("lumiBlock",&lumiBlock);
  ntu_DA -> SetBranchStatus("nPV",      1);   ntu_DA -> SetBranchAddress("nPV",&nVtx);

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
  int CiC_cat;
  int MVA_cat;
  ntu_DA -> SetBranchStatus("R9Ele",         1);   ntu_DA -> SetBranchAddress("R9Ele",R9);
  ntu_DA -> SetBranchStatus("etaSCEle",      1);   ntu_DA -> SetBranchAddress("etaSCEle",scEta);
  ntu_DA -> SetBranchStatus("phiSCEle",      1);   ntu_DA -> SetBranchAddress("phiSCEle",scPhi);
  ntu_DA -> SetBranchStatus("etaEle",        1);   ntu_DA -> SetBranchAddress("etaEle",eta);
  ntu_DA -> SetBranchStatus("phiEle",        1);   ntu_DA -> SetBranchAddress("phiEle",phi);
  ntu_DA -> SetBranchStatus("rawEnergySCEle",1);   ntu_DA -> SetBranchAddress("rawEnergySCEle",scERaw);
  ntu_DA -> SetBranchStatus("energySCEle",   1);   ntu_DA -> SetBranchAddress("energySCEle",scE);
  ntu_DA -> SetBranchStatus("energySCEle_regrCorrSemiParV5_pho",1);   ntu_DA -> SetBranchAddress("energySCEle_regrCorrSemiParV5_pho",scEReg);
  ntu_DA -> SetBranchStatus("eleID",1);   ntu_DA -> SetBranchAddress("eleID",eleID);
  //  ntu_DA -> SetDirectory(0);

  ntu_MC -> SetBranchStatus("R9Ele",         1);   ntu_MC -> SetBranchAddress("R9Ele",R9);
  ntu_MC -> SetBranchStatus("etaSCEle",      1);   ntu_MC -> SetBranchAddress("etaSCEle",scEta);
  ntu_MC -> SetBranchStatus("phiSCEle",      1);   ntu_MC -> SetBranchAddress("phiSCEle",scPhi);
  ntu_MC -> SetBranchStatus("etaEle",        1);   ntu_MC -> SetBranchAddress("etaEle",eta);
  ntu_MC -> SetBranchStatus("phiEle",        1);   ntu_MC -> SetBranchAddress("phiEle",phi);
  ntu_MC -> SetBranchStatus("rawEnergySCEle",1);   ntu_MC -> SetBranchAddress("rawEnergySCEle",scERaw);
  ntu_MC -> SetBranchStatus("energySCEle",   1);   ntu_MC -> SetBranchAddress("energySCEle",scE);
  ntu_MC -> SetBranchStatus("energySCEle_regrCorrSemiParV5_pho",1);   ntu_MC -> SetBranchAddress("energySCEle_regrCorrSemiParV5_pho",scEReg);
  ntu_MC -> SetBranchStatus("eleID",1);   ntu_MC -> SetBranchAddress("eleID",eleID);
  //  ntu_MC -> SetDirectory(0);


  TFile outputFile_DA("/gwterax2/users/martelli/CALIBRATION/data/DoubleElectron-ZSkim-RUN2012ABCD-22Jan-v1-allRange_Friend.root","RECREATE");
  TTree *ntu_DA_F = new TTree("selectedFriend", "" );
  ntu_DA_F->Branch("CiC_cat",&CiC_cat,"CiC_cat/I");
  ntu_DA_F->Branch("MVA_cat",&MVA_cat,"MVA_cat/I");


  //-----------------                                                                                                                 
  // Loop over events DA
  std::cout << std::endl;
  std::cout << ">>> Read data from DA sample" << std::endl;
  int nEntries_DA = ntu_DA -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_DA; ++ientry)
    {
      if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   DA entry " << ientry << " / " << nEntries_DA << "\r" << std::flush;
      ntu_DA->GetEntry(ientry);

      //      std::cout << " runId = " << runId << std::endl;
      std::cout << " S " << lumiBlock << " " << runId << " " << eventNumber << std::endl;

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

      float theta1 = 2*atan(exp(-eta1));
      float theta2 = 2*atan(exp(-eta2));
      float Rt1 = sin(theta1);
      float Rt2 = sin(theta2);
      float Et1 = scEReg1*Rt1;
      float Et2 = scEReg2*Rt2;

      CiC_cat = -1;
      MVA_cat = -1;
      //-------------------
      // loop - Globe
      for(int entry = 0; entry < treeG->GetEntries(); ++entry)
      //      for(int entry = 0; entry < 0; ++entry)
	{
	  //if( entry%10000 == 0 ) std::cout << " >>> reading entry " << entry << " / " << treeG->GetEntries() << "\r" << std::endl;
	  treeG->GetEntry(entry);
	  
	  if(Gitype != 0) continue;

	  if(Glumis != lumiBlock) continue;
	  if(Grun != runId) continue;
	  if(Gevent != eventNumber) continue;

	  std::cout << " G " << Glumis << " " << Grun << " " << Gevent << std::endl;

 	  if(fabs(1.*eta1 - 1.*Geta1) > 0.003 || fabs(1.*eta2 - 1.*Geta2) > 0.003) continue; 
 	  if(R91 < 1.*0.94 && 1.*Gr91 > 0.94) continue;
 	  if(R92 < 1.*0.94 && 1.*Gr92 > 0.94) continue;
	    std::cout << "  " << std::endl;
	    std::cout << " eta1 = " << eta1 << std::endl;
	    std::cout << " Geta1 = " << Geta1 << std::endl;
	    std::cout << " eta2 = " << eta2 << std::endl;
	    std::cout << " Geta2 = " << Geta2 << std::endl;
	    std::cout << " Gr92 = " << Gr92 << std::endl;
	    std::cout << " R92 = " << R92 << std::endl;
	    std::cout << " Gr91 = " << Gr91 << std::endl;
	    std::cout << " R91 = " << R91 << std::endl;


// 	    std::cout << " Gevent = " << Gevent << std::endl;
// 	    std::cout << " Gfull_cat = " << Gfull_cat << std::endl;
// 	    std::cout << " Gdipho_mva_cat = " << Gdipho_mva_cat << std::endl;

	    CiC_cat = int(Gfull_cat);
	    MVA_cat = int(Gdipho_mva_cat);

	    break;
	}

      ntu_DA_F->Fill();
    }

  outputFile_DA.cd();
  ntu_DA_F->Write();
  outputFile_DA.Close();


  ///////////////////////////////////////////////////////////////////

  TFile outputFile_MC("/gwterax2/users/martelli/CALIBRATION/MC/DYToEE_M20_powheg-Summer12-START53-ZSkim-allRange_Friend.root","RECREATE");
  TTree *ntu_MC_F = new TTree("selectedFriend", "");
  ntu_MC_F->Branch("CiC_cat",&CiC_cat,"CiC_cat/I");
  ntu_MC_F->Branch("MVA_cat",&MVA_cat,"MVA_cat/I");


  //-----------------                                                                                                                 
  // Loop over events MC
  std::cout << std::endl;
  std::cout << ">>> Read data from MC sample" << std::endl;
  int nEntries_MC = ntu_MC -> GetEntriesFast();
  for(int ientry = 0; ientry < nEntries_MC; ++ientry)
    {
      if( ientry%100000 == 0 ) std::cout << ">>>>>> reading   MC entry " << ientry << " / " << nEntries_MC << "\r" << std::flush;
      ntu_MC->GetEntry(ientry);

      //      std::cout << " runId = " << runId << std::endl;
      std::cout << " S " << lumiBlock << " " << runId << " " << eventNumber << std::endl;

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

      float theta1 = 2*atan(exp(-eta1));
      float theta2 = 2*atan(exp(-eta2));
      float Rt1 = sin(theta1);
      float Rt2 = sin(theta2);
      float Et1 = scEReg1*Rt1;
      float Et2 = scEReg2*Rt2;

      CiC_cat = -1;
      MVA_cat = -1;

      //-------------------
      // loop - Globe
      for(int entry = 0; entry < treeG->GetEntries(); ++entry)
      //      for(int entry = 0; entry < 0; ++entry)
	{
	  //if( entry%10000 == 0 ) std::cout << " >>> reading entry " << entry << " / " << treeG->GetEntries() << "\r" << std::endl;
	  treeG->GetEntry(entry);
	  
	  if(Gitype != -38) continue;
	  if(Grun != runId) continue;



// 	  if(nVtx != Gnvtx) continue;
// 	  if(nPU != Gpu_n) continue;

	  //	    std::cout << "  " << std::endl;
	  //	  if(eta1 == Geta1 && eta2 == Geta2 && R91 == Gr91 && R92 == Gr92){

//  	  std::cout << " Glumis = " << Glumis << std::endl;
//  	  std::cout << " Grun = " << Grun << std::endl;
//  	  std::cout << " Gevent = " << Gevent << std::endl;
//  	  std::cout << " Gfull_cat = " << Gfull_cat << std::endl;
//  	  std::cout << " Gdipho_mva_cat = " << Gdipho_mva_cat << std::endl;


//	  if(eta1 == Geta1 && eta2 == Geta2){
//	  if(fabs(eta1 - Geta1) < 0.01 && fabs(eta2 - Geta2) < 0.01 ){
 	  if(fabs(1.*eta1 - 1.*Geta1) > 0.003 || fabs(1.*eta2 - 1.*Geta2) > 0.003) continue; 
 	  if(R91 < 1.*0.94 && 1.*Gr91 > 0.94) continue;
 	  if(R92 < 1.*0.94 && 1.*Gr92 > 0.94) continue;

// 	    std::cout << "  " << std::endl;
// 	    std::cout << " eta1 = " << eta1 << std::endl;
// 	    std::cout << " Geta1 = " << Geta1 << std::endl;
// 	    std::cout << " eta2 = " << eta2 << std::endl;
// 	    std::cout << " Geta2 = " << Geta2 << std::endl;
// 	    std::cout << " Gr92 = " << Gr92 << std::endl;
// 	    std::cout << " R92 = " << R92 << std::endl;
// 	    std::cout << " Gr91 = " << Gr91 << std::endl;
// 	    std::cout << " R91 = " << R91 << std::endl;
// 	    std::cout << " G " << Glumis << " " << Grun << " " << Gevent << std::endl;

// 	    std::cout << " Gevent = " << Gevent << std::endl;
// 	    std::cout << " Gfull_cat = " << Gfull_cat << std::endl;
// 	    std::cout << " Gdipho_mva_cat = " << Gdipho_mva_cat << std::endl;

	    CiC_cat = int(Gfull_cat);
	    MVA_cat = int(Gdipho_mva_cat);

	    break;
	}

      ntu_MC_F->Fill();
    }


  
  

  outputFile_MC.cd();
  ntu_MC_F->Write();
  outputFile_MC.Close();

  return 0;
}




//g++ -Wall -o CloneMVAHggTree `root-config --cflags --glibs` CloneMVAHggTree.cpp

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
#include "TVirtualFitter.h"
#include "TRandom.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>


// Cholesky decomposition of a matrix
double* cholesky(double* A, int n);

int GetSingleCategory(const float& scEta1, const float& R91);

int main()
{
  int nSamples = 1;

  std::string type = "8TeVCiC";
  //std::string type = "8TeVMVA";
  //std::string type = "7TeVCiC";
  //std::string type = "7TeVMVA";

    std::string process = "125500";
  //  std::string process = "125400";
  //  std::string process = "125100";
  //  std::string process = "125000";

  //------------
  // input trees
  TFile* f;
  if(type == "8TeVCiC") f = new TFile("../data/Hgg_signalTrees/jackCiC_8TeV.root", "read");
  if(type == "7TeVCiC") f = new TFile("../data/Hgg_signalTrees/jackCiC_7TeV.root", "read");
  if(type == "8TeVMVA") f = new TFile("../data/Hgg_signalTrees/jackMVA_8TeV.root", "read");
  if(type == "7TeVMVA") f = new TFile("../data/Hgg_signalTrees/jackMVA_7TeV.root", "read");
  TTree** trees = new TTree*[nSamples];
  trees[0] = (TTree*)f->Get("opttree");
  
  float itype;
  float full_weight;
  float full_cat;
  float dipho_mva;
  float dipho_pt;
  float mass;
  float et1, et2, eta1, eta2, r91, r92;
  
  for(int iSample = 0; iSample < nSamples; ++iSample)
  {
    trees[iSample]->SetBranchStatus("*",0);
    trees[iSample]->SetBranchStatus("itype", 1);          trees[iSample]->SetBranchAddress("itype", &itype);
    trees[iSample]->SetBranchStatus("full_weight", 1);    trees[iSample]->SetBranchAddress("full_weight", &full_weight);
    trees[iSample]->SetBranchStatus("full_cat", 1);       trees[iSample]->SetBranchAddress("full_cat", &full_cat);
    trees[iSample]->SetBranchStatus("dipho_mva", 1);      trees[iSample]->SetBranchAddress("dipho_mva", &dipho_mva);
    trees[iSample]->SetBranchStatus("dipho_pt", 1);      trees[iSample]->SetBranchAddress("dipho_pt", &dipho_pt);
    trees[iSample]->SetBranchStatus("mass", 1);           trees[iSample]->SetBranchAddress("mass", &mass);
    trees[iSample]->SetBranchStatus("et1", 1);             trees[iSample]->SetBranchAddress("et1", &et1);
    trees[iSample]->SetBranchStatus("et2", 1);             trees[iSample]->SetBranchAddress("et2", &et2);
    trees[iSample]->SetBranchStatus("eta1", 1);           trees[iSample]->SetBranchAddress("eta1", &eta1);
    trees[iSample]->SetBranchStatus("eta2", 1);           trees[iSample]->SetBranchAddress("eta2", &eta2);
    trees[iSample]->SetBranchStatus("r91", 1);           trees[iSample]->SetBranchAddress("r91", &r91);
    trees[iSample]->SetBranchStatus("r92", 1);           trees[iSample]->SetBranchAddress("r92", &r92);
  }
  
  
  TTree* newTree;
  // = trees[0]->CloneTree(0);

  if(trees[0]->GetEntries() == 0 )
    {
      std::cout << "Error: input file is empty" << std::endl;
      return -1;
    }
  else newTree = trees[0]->CloneTree(0);
  
  //-------------------
  // loop - all samples
  int iSample = 0;
  for(int entry = 0; entry < trees[iSample]->GetEntries(); ++entry)
    {
      if( entry%10000 == 0 ) std::cout << " >>> reading entry " << entry << " / " << trees[iSample]->GetEntries() << "\r" << std::endl;
      trees[iSample]->GetEntry(entry);

      //      if(itype < -12499) std::cout << " itype = " << itype << std::endl;
      //            if(itype == -125500) newTree->Fill();
      //      if(itype == -125400) newTree->Fill();
      //      if(itype == -125100) newTree->Fill();
      //      if(itype == -125000) newTree->Fill();
//       else continue;
	// || itype == -125400 || itype == -125100 || itype == -125000) newTree->Fill();

      //      std::cout << " et1 = " << et1 << std::endl; 
      

      if(itype == (-1. * atof(process.c_str())) ) newTree->Fill();

    }


//  TFile* f = new TFile("../data/Hgg_signalTrees/Hgg_Dump_CiC_8TeV.root", "read");
//   TFile* f = new TFile("../data/Hgg_signalTrees/Hgg_Dump_CiC_7TeV.root", "read"); 
//   TFile* f = new TFile("../data/Hgg_signalTrees/Hgg_Dump_MVA_8TeV.root", "read"); 
//   TFile* f = new TFile("../data/Hgg_signalTrees/Hgg_Dump_MVA_7TeV.root", "read");   

  TFile* outputFile;
  if(type == "8TeVCiC") outputFile = new TFile(("../data/Hgg_signalTrees/Hgg_Dump_CiC_8TeV_"+process+".root").c_str(),"RECREATE");
  if(type == "7TeVCiC") outputFile = new TFile(("../data/Hgg_signalTrees/Hgg_Dump_CiC_7TeV_"+process+".root").c_str(),"RECREATE");
  if(type == "8TeVMVA") outputFile = new TFile(("../data/Hgg_signalTrees/Hgg_Dump_MVA_8TeV_"+process+".root").c_str(),"RECREATE");
  if(type == "7TeVMVA") outputFile = new TFile(("../data/Hgg_signalTrees/Hgg_Dump_MVA_7TeV_"+process+".root").c_str(),"RECREATE");
  
  outputFile->cd();
  
  newTree->SetName("opttree");
  newTree->Print();
  newTree->AutoSave();
  newTree->Write();
  outputFile->Close();
  
  return 0;
}



double* cholesky(double* A, int n)
{
  double *L = (double*)calloc(n * n, sizeof(double));
  if (L == NULL)
    exit(EXIT_FAILURE);
 
  for (int i = 0; i < n; i++)
    for (int j = 0; j < (i+1); j++) {
      double s = 0;
      for (int k = 0; k < j; k++)
	s += L[i * n + k] * L[j * n + k];
      L[i * n + j] = (i == j) ?
	sqrt(A[i * n + i] - s) :
	(1.0 / L[j * n + j] * (A[i * n + j] - s));
    }
  
  return L;
}




int GetSingleCategory(const float& scEta1, const float& R91)
{
//   std::cout << " in GetSingleCategory scEta1 = " << scEta1 << std::endl;
//   std::cout << " in GetSingleCategory R91 = " << R91 << std::endl;

  if( fabs(scEta1) < 1.4442 )
    {
      if( R91 > 0.94 ) return 0;
      if( R91 < 0.94 ) return 1;
    }

  if( fabs(scEta1) > 1.5600 )
    {
      if( R91 > 0.94 ) return 2;
      if( R91 < 0.94 ) return 3;
    }

  return -1;
}

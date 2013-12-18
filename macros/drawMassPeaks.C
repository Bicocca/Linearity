//g++ -Wall -o linearityPlots_MZ `root-config --cflags --glibs` linearityPlots_MZ.cpp

#include "/Users/Arabella/Public/setTDRStyle.h"
// #include "ntpleUtils.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TVirtualFitter.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"


int drawMassPeaks()
{
  //Check if all nedeed arguments to parse are there
  
  int nCat = 4;
  int* evtsPerPoint = new int[nCat];
  for(int iCat = 0; iCat < nCat; ++iCat)
    evtsPerPoint[iCat] = -1;


  std::vector<std::string> directories;
  directories.push_back( "MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_S0S5/");
  //  directories.push_back( "../5Dicembre/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_noEtResidual/");
  directories.push_back( "../Ottobre_2013/Linearity_Uncertainty/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_noEtdep_17Nov/");
//   directories.push_back( "MZ_stdCat_nonGlobe_powheg-runDependent_eleTunedRegV5_Dphi3p15_noEtdep/");
//   directories.push_back( "MZ_stdCat_nonGlobe_powheg-runDependent_eleTunedRegV5_Dphi3p15_Etdep/");
//   directories.push_back( "MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_AllMC_noEtdep/");
//   directories.push_back( "MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_AllMC_Etdep/");
  int nDir = directories.size();
  
  std::vector<std::string> methods;
  methods.push_back( "recursiveMean" );
  
  std::vector<std::string> nameTrial;
  nameTrial.push_back("EtDep");
  //  nameTrial.push_back("EtResidualDep");
  nameTrial.push_back("noEtDep");
//   nameTrial.push_back("eleTuned_noEtDep");
//   nameTrial.push_back("eleTuned_EtDep");
//   nameTrial.push_back("allMC_noEtDep");
//   nameTrial.push_back("allMC_EtDep");


  //--------------------
  // Set fitting options
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  
  //setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  std::vector<int> colors2;
  std::vector<int> linestyles;
  std::vector<int> markerstyles;
  colors2.push_back(kRed+1);
  colors2.push_back(kOrange+1);
  colors2.push_back(kGreen+1);
  colors2.push_back(kBlue+1);
  colors2.push_back(kCyan+1);
  colors2.push_back(kMagenta+1);
  
  std::vector<int> colors;
  std::vector<int> linestyles;
  std::vector<int> markerstyles;
  colors.push_back(kRed+1);
  colors.push_back(kOrange+1);
  colors.push_back(kGreen+1);
  colors.push_back(kBlue+1);
  linestyles.push_back(1);
  linestyles.push_back(1);
  linestyles.push_back(2);
  linestyles.push_back(2);
  markerstyles.push_back(20);
  markerstyles.push_back(26);
  markerstyles.push_back(32);
  markerstyles.push_back(32);

    
  
  //----------------
  // Define canvases
  
  TCanvas* c_all = new TCanvas("c_all","c_all");
  c_all -> Divide(2,2);
  
  TCanvas** c = new TCanvas*[nCat];
  TLegend** legend = new TLegend*[nCat]; 
  //----------------------------
  // Define infiles and canvases
  for(int iCat = 3; iCat >= 0; --iCat) {
      std::cout << "\n***************** cat: " << iCat << " *****************" << std::endl;
      
      legend[iCat] = new TLegend(0.16, 0.75, 0.30, 1.);
      legend[iCat] -> SetFillColor(kWhite);
      legend[iCat] -> SetFillStyle(1001);  
      legend[iCat] -> SetTextFont(42);  
      legend[iCat] -> SetTextSize(0.03);

      TLatex* latex = new TLatex(0.70,0.90,Form("cat%d",iCat));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextColor(colors.at(iCat));
      latex -> SetTextSize(0.04);
      
      std::cout << ">>> iCat = " << iCat << std::endl;
      std::string* inFileNames = new std::string[nDir];
      TH1F** hDA = new TH1F*[nDir];
      TH1F** hMC = new TH1F*[nDir];
      TPaveStats** pDA = new TPaveStats*[nDir];
      TPaveStats** pMC = new TPaveStats*[nDir];

      //      TCanvas* cprova = new TCanvas();
      //      TPad* pad = new TPad("pad","",0,0,1,1);

      for(unsigned int iDir = 0; iDir < nDir; ++iDir)
	{
	  std::string directory = directories.at(iDir);
    	  std::string baseFileName = "studyLinearity_MZ_";
  
	  char EvtString[50];
	  sprintf(EvtString,"cat%d_%devtsPerPoint",iCat,evtsPerPoint[iCat]);
        
	  inFileNames[iDir] = directory + "/" + baseFileName + std::string(EvtString) + ".root";
	  std::cout << ">>> inFileName: " << inFileNames[iDir] << std::endl;
	  TFile* f = TFile::Open((inFileNames[iDir]).c_str(),"READ");

	  char graphName[50];
	  sprintf(graphName,"h_mee_DA");
	  hDA[iDir] = (TH1F*)( f->Get(graphName) );
	  hDA[iDir]->SetName(("DA "+nameTrial.at(iDir)).c_str() );
	  sprintf(graphName,"h_mee_MC");
	  hMC[iDir] = (TH1F*)( f->Get(graphName) );
	  hMC[iDir]->SetName(("MC "+nameTrial.at(iDir)).c_str());

	  hDA[iDir] -> GetXaxis() -> SetTitle("m_{ee}");
	  hDA[iDir] -> GetXaxis() -> SetTitleSize(0.05);
	  hDA[iDir] -> GetYaxis() -> SetTitleSize(0.05);
	  hDA[iDir] -> GetXaxis() -> SetTitleOffset(1.10);
	  hDA[iDir] -> GetYaxis() -> SetTitleOffset(1.28);
	  hDA[iDir] -> GetXaxis() -> SetLabelSize(0.04);
	  hDA[iDir] -> GetYaxis() -> SetLabelSize(0.04);
	  hDA[iDir] -> GetYaxis() -> SetNdivisions(405);

	  hDA[iDir] -> SetLineColor(colors.at(iDir));
	  //	  hDA[iDir] -> SetLineStyle(1);
	  hDA[iDir] -> SetLineWidth(2);

	  hMC[iDir] -> SetLineColor(colors.at(iDir+2));
 	  hMC[iDir] -> SetLineWidth(2);
// 	  hMC[iDir] -> SetLineStyle(2);

	  gStyle->SetOptStat(1111);

	  if( iDir == 0 ) {
	    c[iCat] = new TCanvas();
// 	    pad->Draw();
// 	    pad->cd();
	    c[iCat] -> cd();
	    c[iCat] -> SetGridx();
	    c[iCat] -> SetGridy();
	    //c[iCat] -> SetLogx();
	  
	    //	    cprova->cd();
	    c[iCat] -> cd();
            hDA[iDir] -> Draw("hist");
	    gPad->Update();
	    pDA[iDir] = (TPaveStats*)hDA[iDir]->FindObject("stats");
//  	    pDA[iDir]->SetX1NDC(0.1);
//  	    pDA[iDir]->SetX2NDC(0.3);
//  	    pDA[iDir]->SetY1NDC(0.7);
//  	    pDA[iDir]->SetY2NDC(0.9);


	    hMC[iDir] -> Draw("hist, sames");
	    gPad->Update();
	    pMC[iDir] = (TPaveStats*)hMC[iDir]->FindObject("stats");

//  	    pMC[iDir]->SetX1NDC(0.3);
//  	    pMC[iDir]->SetX2NDC(0.5);
//  	    pMC[iDir]->SetY1NDC(0.7);
//  	    pMC[iDir]->SetY2NDC(0.9);
//  	    pMC[iDir]->SetTextColor(kRed);
//  	    pMC[iDir]->SetTextColor(kBlue);

//             c[iCat] -> cd();
//             hDA[iDir] -> Draw("hist");
//             hMC[iDir] -> Draw("hist, sames");

            c_all -> cd(iCat+1);
            hDA[iDir] -> Draw("hist"); 
	    hMC[iDir] -> Draw("hist,sames");
          }
          else
          {
            c[iCat] -> cd();
	    //	    cprova->cd();
	    hDA[iDir] -> Draw("hist, sames");
            gPad->Update();
            pDA[iDir] = (TPaveStats*)hDA[iDir]->FindObject("stats");

	    hMC[iDir] -> Draw("hist, sames");
	    gPad->Update();
	    pMC[iDir] = (TPaveStats*)hMC[iDir]->FindObject("stats");

//  	    pDA[iDir]->SetX1NDC(0.1);
//  	    pDA[iDir]->SetX2NDC(0.3);
//  	    pDA[iDir]->SetY1NDC(0.7);
//  	    pDA[iDir]->SetY2NDC(0.9);
//  	    pMC[iDir]->SetX1NDC(0.3);
// 	    pMC[iDir]->SetX2NDC(0.5);
//  	    pMC[iDir]->SetY1NDC(0.7);
//  	    pMC[iDir]->SetY2NDC(0.9);
//  	    pMC[iDir]->SetTextColor(kRed);
//  	    pMC[iDir]->SetTextColor(kBlue);


//             hDA[iDir] -> Draw("hist,sames");
// 	    hMC[iDir] -> Draw("hist,sames");
            c_all -> cd(iCat+1);
            hDA[iDir] -> Draw("hist,sames");
	    hMC[iDir] -> Draw("hist,sames");
          }

	  legend[iCat]->AddEntry(hDA[iDir], ("DA "+nameTrial.at(iDir)).c_str(), "l");
	  legend[iCat]->AddEntry(hMC[iDir], ("MC "+nameTrial.at(iDir)).c_str(), "l");

	  if(iDir == nDir-1){	 
	    c[iCat] -> cd();
	    legend[iCat] -> Draw("same");
	    latex -> Draw("same");
	    
	    char pdfName[250];
	    sprintf(pdfName,"PLOTS/MZ_%s_cat%d.pdf", (nameTrial.at(iDir)).c_str(),iCat);
	    c[iCat] -> Print(pdfName,"pdf");
	    char pngName[250];
	    sprintf(pngName,"PLOTS/MZ_%s_cat%d.png", (nameTrial.at(iDir)).c_str(),iCat);
	    c[iCat] -> Print(pngName,"png");
	  }
	}
      std::cout << ">>> inFileName: fine " << std::endl;
  }
  
}
 
 

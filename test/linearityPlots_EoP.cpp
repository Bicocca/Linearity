// g++ -Wall -o drawScaleMZPlots `root-config --cflags --glibs` setTDRStyle.cc drawScaleMZPlots.cpp

#include "setTDRStyle.h"
#include "ntpleUtils.h"

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


bool MCClosure   = false;
bool drawFitFunc = true;

std::string analysis  = "stdCat";
std::string fitMethod = "exp";
std::string directory = "data/2012/Winter2013/EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/";



int main(int argc, char** argv)
{
  //Check if all nedeed arguments to parse are there
  if( argc < 2 )
  {
    std::cerr << ">>> linearityPlots_MZ::usage: " << argv[0] << " evtsPerPoint_cat0 (evtsPerPoint_cat1 evtsPerPoint_cat2 ...)" << std::endl;
    return -1;
  }
  
  int nCat = argc - 1;
  int* evtsPerPoint = new int[nCat];
  for(int iCat = 0; iCat < nCat; ++iCat)
    evtsPerPoint[iCat] = atoi(argv[iCat+1]);
  
  
  
  //-----------------------------
  // Decide which methods to draw
  
  std::vector<std::string> methods;
  //methods.push_back( "fit" );
  //methods.push_back( "gausFit" );
  //methods.push_back( "mean" );
  //methods.push_back( "recursiveMean" );
  methods.push_back( "smallestInterval" );
  int nMethods = methods.size();
  
  
  
  //--------------------
  // Set fitting options
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.17);
  gStyle->SetLabelSize(0.04,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  std::vector<int> colors;
  std::vector<int> linestyles;
  
  for(int iCat = 0; iCat < int((nCat+1)/2); ++iCat)
  {
    colors.push_back(kBlue+iCat);
    colors.push_back(kRed+iCat);
    
    linestyles.push_back(1);
    linestyles.push_back(1);
  }  
  
  
  
  //---------------------------
  // Define infiles and canvases
  
  std::string* inFileNames = new std::string[nCat];
  
  std::string baseFileName = "studyLinearity_EoP_";
  if( MCClosure == true) baseFileName += "_MCClosure";
  
  
  TCanvas** c01 = new TCanvas*[2];
  TCanvas** c23 = new TCanvas*[2];
  
  for(int iMeth = 0; iMeth < nMethods; ++iMeth)
  {
    std::cout << ">>> method: " << methods.at(iMeth) << std::endl;
    
    TLegend* legend01 = new TLegend(0.16, 0.77, 0.30, 0.92);
    legend01 -> SetFillColor(kWhite);
    legend01 -> SetFillStyle(1001);  
    legend01 -> SetTextFont(42);  
    legend01 -> SetTextSize(0.05);
    
    TLegend* legend23 = new TLegend(0.16, 0.77, 0.30, 0.92);
    legend23 -> SetFillColor(kWhite);
    legend23 -> SetFillStyle(1001);  
    legend23 -> SetTextFont(42);  
    legend23 -> SetTextSize(0.05);    
    
    char latexText[50];
    sprintf(latexText,"analysis: %s   scale estimator: %s",analysis.c_str(),(methods.at(iMeth)).c_str());
    TLatex* latex = new TLatex(0.14,0.96,latexText);                                                                                                                              
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    
    
    
    double* scale_MZ = new double[nCat];
    //double* syst     = new double[nCat];
    //double* systErr  = new double[nCat];
    TGraphAsymmErrors** g = new TGraphAsymmErrors*[nCat];
    
    TF1** f_fit = new TF1*[nCat];
    TH1F** hint = new TH1F*[nCat];
    
    for(int iCat = 0; iCat < nCat; ++iCat)
    {
      std::cout << "\n*****************" << std::endl;
      char EvtString[50];
      sprintf(EvtString,"cat%d_%devtsPerPoint",iCat,evtsPerPoint[iCat]);
      
      inFileNames[iCat] = directory + "/" + baseFileName + std::string(EvtString) + ".root";
      TFile* f     = TFile::Open(     (inFileNames[iCat]).c_str(),"READ" );
      //std::cout << ">>> inFileName: " << inFileNames[iCat] << std::endl;
      
      char graphName[50];
      sprintf(graphName,"step1/scale_%s_DAOverMC",methods.at(iMeth).c_str());
      g[iCat] = (TGraphAsymmErrors*)( f->Get(graphName) );
      
      
      if( iCat == 0 )
      {
        c01[iMeth] = new TCanvas();
        c01[iMeth] -> cd();
        c01[iMeth] -> SetGridx();
        c01[iMeth] -> SetGridy();
        //c01[iMeth] -> SetLogx();
      }
      if( iCat == 1 )
      {
        c01[iMeth] -> cd();
      }
      if( iCat == 2 )
      {
        c23[iMeth] = new TCanvas();
        c23[iMeth] -> cd();
        c23[iMeth] -> SetGridx();
        c23[iMeth] -> SetGridy();
        //c23[iMeth] -> SetLogx();
      }
      if( iCat == 3 )
      {
        c23[iMeth] -> cd();
      }
      
      g[iCat] -> SetPoint(g[iCat]->GetN(),1500.,1.);
      g[iCat] -> GetXaxis() -> SetRangeUser(30,120.);        
      g[iCat] -> GetXaxis() -> SetMoreLogLabels();
      g[iCat] -> GetXaxis() -> SetTitle("E_{T} [GeV]");
      g[iCat] -> GetYaxis() -> SetTitle("scale_{data} / scale_{MC}");
      g[iCat] -> GetXaxis() -> SetTitleSize(0.05);
      g[iCat] -> GetYaxis() -> SetTitleSize(0.05);
      g[iCat] -> GetXaxis() -> SetTitleOffset(1.10);
      g[iCat] -> GetYaxis() -> SetTitleOffset(1.10);
      g[iCat] -> GetXaxis() -> SetLabelSize(0.04);
      g[iCat] -> GetYaxis() -> SetLabelSize(0.04);
      g[iCat] -> SetMarkerSize(0.7);
      g[iCat] -> SetMarkerStyle(20+iCat%2);
      g[iCat] -> SetMarkerColor(colors.at(iCat));
      g[iCat] -> SetLineColor(colors.at(iCat));
      g[iCat] -> SetLineStyle(linestyles.at(iCat));
      
      g[iCat] -> SetMinimum(0.980);
      g[iCat] -> SetMaximum(1.020);
      if( drawFitFunc == false )
      {
        if( iCat == 0 || iCat == 2 ) g[iCat] -> Draw("APL");
        else                       { g[iCat] -> Draw("PL,same"); g[iCat-1] -> Draw("PL,same"); }
      }
      else
      {
        if( iCat == 0 || iCat == 2 ) g[iCat] -> Draw("AP");
        else                       { g[iCat] -> Draw("P,same"); g[iCat-1] -> Draw("P,same"); }
      }
      
      char funcName[150];
      sprintf(funcName,"f_prefit_%s_%d",methods.at(iMeth).c_str(),iCat);
      TF1* f_prefit;
      if( fitMethod == "pol1")
      {
        f_prefit = new TF1(funcName,"[0]+[1]*x",65.,210.);
        f_prefit -> SetParameters(1.,1.);         
      }
      if( fitMethod == "exp")
      {
        f_prefit = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-90.)))+[2]",65.,210.);
        f_prefit -> SetParameters(0.005,0.02,0.);
        f_prefit -> SetParLimits(0,0.,1.);
        f_prefit -> SetParLimits(1,0.,0.05);
      }
      g[iCat] -> Fit(funcName,"QERNS");
      
      
      if( drawFitFunc == true )
      {
        for(int point = 0; point < g[iCat]->GetN(); ++point)
        {
          double ey = g[iCat] -> GetErrorY(point);
          g[iCat] -> SetPointEYhigh(point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
          g[iCat] -> SetPointEYlow (point,ey*sqrt(f_prefit->GetChisquare()/f_prefit->GetNDF()));
        }
      }
      
      
      sprintf(funcName,"f_fit_%s_%d",methods.at(iMeth).c_str(),iCat);
      if( fitMethod == "pol1" )
      {
        f_fit[iCat] = new TF1(funcName,"[0]+[1]*x",65.,210.);
        f_fit[iCat] -> SetParameters(1.,1.);
      }
      if( fitMethod == "exp")
      {
        f_fit[iCat] = new TF1(funcName,"1.+[0]*(1.-exp(-1.*[1]*(x-90.)))+[2]",65.,210.);
        f_fit[iCat] -> SetParameters(0.005,0.02,0.);
        f_fit[iCat] -> SetParLimits(0,-1.,1.);
        f_fit[iCat] -> SetParLimits(1,0.,0.05);
      }
      f_fit[iCat] -> SetLineColor(colors.at(iCat));
      f_fit[iCat] -> SetLineWidth(2);
      f_fit[iCat] -> SetLineStyle(linestyles.at(iCat));
      
      TFitResultPtr fitResult;
      int fitStatus = -1;
      int nTrials = 0;
      
      while( (fitStatus != 1) && (nTrials < 2) )
      {
        fitResult = g[iCat] -> Fit(funcName,"QENRS+");
        fitStatus = fitResult;
        if( fitStatus == 1 ) break;
        ++nTrials;
      }
      
      

      if( drawFitFunc == true )
      {
        f_fit[iCat] -> Draw("same");
        hint[iCat] = new TH1F("hint","",1000,65.,1065.);
        (TVirtualFitter::GetFitter()) -> GetConfidenceIntervals(hint[iCat],0.68);
        hint[iCat] -> SetMarkerSize(0);
        hint[iCat] -> SetFillColor(colors.at(iCat));
        hint[iCat] -> SetFillStyle(3001);
        hint[iCat] -> Draw("same,E4");
      }  
      
      
            
      char legendText[50];
      sprintf(legendText,"cat. %d",iCat);
      if( iCat == 0 || iCat == 1) legend01 -> AddEntry(g[iCat],legendText,"PL");
      else                      legend23 -> AddEntry(g[iCat],legendText,"PL");
      
      
      
      TH1F* h_Ht_MC = (TH1F*)( f->Get("h_Ht_MC") );
      double HT_Z = h_Ht_MC->GetMean();
      scale_MZ[iCat] = f_fit[iCat] -> Eval(HT_Z);
      
      std::cout << ">>>>>> iCat: " << iCat
                << std::fixed
                << "   rel. scale(MZ): "  << std::setprecision(4) << scale_MZ[iCat]
                << std::endl;
    }
    
    c01[iMeth] -> cd();
    legend01 -> Draw("same");
    latex -> Draw("same");
    c23[iMeth] -> cd();
    legend23 -> Draw("same");
    latex -> Draw("same");
    
    
    char pdfName[250];
    char pngName[250];
    
    if( MCClosure == false )
    {
      sprintf(pdfName,"%s/scale_%s_%s_%s_cat01_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c01[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat01_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c01[iMeth] -> Print(pngName,"png");
      
      sprintf(pdfName,"%s/scale_%s_%s_%s_cat23_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c23[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat23_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c23[iMeth] -> Print(pngName,"png");
    }
    else
    {
      sprintf(pdfName,"%s/scale_%s_%s_%s_cat01_MCClosure_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c01[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat01_MCClosure_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c01[iMeth] -> Print(pngName,"png");

      sprintf(pdfName,"%s/scale_%s_%s_%s_cat23_MCClosure_DAOverMC.pdf",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pdfName << std::endl;
      c23[iMeth] -> Print(pdfName,"pdf");
      sprintf(pngName,"%s/scale_%s_%s_%s_cat23_MCClosure_DAOverMC.png",directory.c_str(),analysis.c_str(),(methods.at(iMeth)).c_str(),fitMethod.c_str());
      std::cout << ">>> saving file " << pngName << std::endl;
      c23[iMeth] -> Print(pngName,"png");
    }
  }
}

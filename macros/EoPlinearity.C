
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TClonesArray.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TString.h"
#include "TProfile.h"
#include "TArrow.h"
#include "TText.h"
#include "TLatex.h"
#include "TPaveStats.h"

// #include "RooMath.h"
// #include "RooComplex.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

/*
TProfile *hFit = new TProfile("hFit","hFit",180,80,105,"s"); 

// Cholesky decomposition of a matrix
double *cholesky(double *A, int n) {
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

*/

void EoPlinearity(){

  // inits
  TGraphAsymmErrors *gMee[4];
  TH1F *hInt[4]; 
  char graphName[20]; 
  TF1 *user = new TF1("user","[0]+[1]*x",50.,250.); 
  user->SetLineWidth(2); 
  user->SetLineColor(kAzure+7); 

  // init graphics
  TCanvas *c = new TCanvas("c","c",50,50,700,500);
  c -> Divide(2,2);

  //read file
  TFile *_file1 = TFile::Open("MZ.root");

  // loop on categories
  for (int icat = 0; icat<4; icat++) {

    // canvasize
    c -> cd(icat+1); 
    gPad->SetGridy();
    gPad -> SetLeftMargin(.25); 
    gPad ->SetRightMargin(0.25);
    gPad -> SetTopMargin(0.2); 
    gPad ->SetBottomMargin(0.2);
    if (icat%2==0) 
      gPad ->SetRightMargin(0.002);
    else 
      gPad->SetLeftMargin(0.002); 
    
    if (icat<2) 
      gPad ->SetBottomMargin(0.005); 
    else 
      gPad->SetTopMargin(0.005);

    //    TH1F *h = (TH1F*)gPad->DrawFrame(50,-0.0051,230.,0.0051); 
    TH1F *h = (TH1F*)gPad->DrawFrame(50,1-0.0051,230.,1+0.0051); 
    h -> GetYaxis() -> SetTitleSize(0.1); 
    h -> GetYaxis() -> SetLabelSize(0.07); 
    h -> GetYaxis() -> SetTitleOffset(1.2); 
    h -> GetXaxis() -> SetTitleSize(0.1); 
    h -> GetXaxis() -> SetLabelSize(0.07); 
    //    h -> GetXaxis() -> SetTitleOffset(1.2); 
    h -> GetXaxis() -> SetTitle("H_{T} or 2E_{T} (GeV)"); 
    //    h -> GetYaxis() -> SetTitle("#Deltam_{ee}^{  Data-MC}"); 
    h -> GetYaxis() -> SetTitle("Relative Response"); 

    // get graph
    sprintf(graphName,"graph_cat%d",icat); 
    gMee[icat] = (TGraphAsymmErrors*)_file1->Get(graphName); 

    // Set plot nicer
    gMee[icat] -> SetLineWidth(1);
    gMee[icat] -> SetLineColor(1);
    gMee[icat] -> SetMarkerStyle(20); 
    gMee[icat] -> SetMarkerSize(.7); 
    if (icat==0)     gMee[icat] -> SetMarkerColor(kRed+1); 
    if (icat==2)     gMee[icat] -> SetMarkerColor(kGreen+1); 
    if (icat==1)     gMee[icat] -> SetMarkerColor(kOrange+1); 
    if (icat==3)     gMee[icat] -> SetMarkerColor(kBlue+1); 

    // fit silent, and scale errors
    gMee[icat]->Fit("user","QN"); 
    // scale; 
    double scale = sqrt(user->GetChisquare()/user->GetNDF()); 
    for (int i=0;i<gMee[icat]->GetN();i++) {
      double ey=1,ex=1;
      ey = gMee[icat]->GetErrorYhigh(i); 
      ex = gMee[icat]->GetErrorXhigh(i); 
      gMee[icat]->SetPointError(i,ex,ex,ey*scale,ey*scale);
      gMee[icat]->GetPoint(i,ex,ey);
      //      g[icat]->SetPoint(i,ex,ey-1);
    }
    // fit again to get the correct error estimate
    gMee[icat]->Fit("user","QN"); 

    // confidence intervals
    hInt[icat] = (TH1F*)h->Clone("Conf.Interval"); 
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hInt[icat],0.68);
    hInt[icat]->SetStats(kFALSE);
    hInt[icat]->SetFillStyle(3001);
    hInt[icat]->SetFillColor(kAzure+7);
    // redraw, with confidence belt

    hInt[icat]->Draw("e3 same");
    //    user->Draw("same"); 
    gMee[icat]->Draw("P");
    if (icat>1) h->GetYaxis()->SetRangeUser(0.9849,1.0159);
  }

  TLegend *tl = new TLegend(0.8,0.2,1.00,0.7); 
  tl ->SetFillColor(0);
  tl ->SetBorderSize(0); 
  tl ->SetTextFont(42); 
  //  tl->SetNColumns(2);

  tl -> AddEntry(gMee[0],"#gamma CiC0","PL"); 
  tl -> AddEntry(gMee[1],"#gamma CiC1","PL"); 
  tl -> AddEntry(gMee[2],"#gamma CiC2","PL"); 
  tl -> AddEntry(gMee[3],"#gamma CiC3","PL"); 

  c->cd(2); tl->Draw();

  TGraphAsymmErrors *gEoP[4];
  TGraphAsymmErrors *gEoP2[4];
  TFile *_file0 = TFile::Open("EoP.root");

  //  double yn[4]={1.0018,1.,1.0065,0.999};
  double yn[4]={1.0024,1.0013,1.0204,1.0140};
  for (int icat=0;icat<4;icat++) {
    sprintf(graphName,"graph_cat%d",icat); 
    gEoP[icat] = (TGraphAsymmErrors*)_file0->Get(graphName); 
    gEoP2[icat] = (TGraphAsymmErrors*)_file0->Get(graphName); 
    for (int i=0;i<gEoP[icat]->GetN();i++) {
      double y=1,x=1;
      gEoP[icat]->GetPoint(i,x,y);
      gEoP[icat]->SetPoint(i,2*x,y/yn[icat]);
      gEoP2[icat]->SetPoint(i,2*x,y);
     }

    c->cd(icat+1);
    gEoP[icat]->SetMarkerStyle(3); 
    gEoP[icat]->SetMarkerSize(1.); 
    if (icat==0)     gEoP[icat] -> SetMarkerColor(kRed+1); 
    if (icat==2)     gEoP[icat] -> SetMarkerColor(kGreen+1); 
    if (icat==1)     gEoP[icat] -> SetMarkerColor(kOrange+1); 
    if (icat==3)     gEoP[icat] -> SetMarkerColor(kBlue+1); 
    gEoP[icat]->Draw("P");
    gEoP2[icat]->Draw("P");
  }

  TLegend *tt = new TLegend(0.,0.8,0.75,.95); 
  tt ->SetFillColor(0);
  tt ->SetBorderSize(0); 
  tt ->SetTextFont(42); 
  tt->SetNColumns(2);

  tt -> AddEntry(gMee[0],"m_{ee} data/MC","PL"); 
  tt -> AddEntry(gEoP[0],"E/p data/MC","PL"); 

  c->cd(2); tt->Draw();




}




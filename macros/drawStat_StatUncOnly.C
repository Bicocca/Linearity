#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

void FindSmallestInterval(double& mean, double& meanErr, double& min, double& max,
                          TH1F* histo,
                          const double& fraction, const bool& verbosity);



void drawStat_StatUncOnly(std::string analysis, std::string HggMethod)
{


  //  gROOT->ProcessLine(".x ~/Public/rootLogon.C");
  //  gROOT->ProcessLine(".x ~/Public/setTDRStyle.C");
  //gROOT->ProcessLine(".x ~/Public/style.C");
  //  gROOT->SetStyle("Plain");

  TFile* f = TFile::Open("rescaledHgg_histos.root","READ");
  TLatex* latex;
  
  //  int rebin = 1;
  //  int rebin = 12;
  int rebin = 6;
  int nCats = 5;
  int nTrials = 1000;

  TFile** fHt = new TFile*[nCats];
  TH1F** hET = new TH1F*[nCats];
  TH1F** hHT = new TH1F*[nCats];
  TH1F** hMCorr = new TH1F*[nCats];
  TH1F** hSingPho = new TH1F*[nCats];

  std::vector<int> colors;
  for(int iCat = 0; iCat < nCats; ++iCat)
    {
      colors.push_back(kRed+1);
      colors.push_back(kOrange);
      colors.push_back(kGreen+1);
      colors.push_back(kBlue);
      colors.push_back(kBlack);
    }

  
  for(int iCat = 0; iCat < nCats; ++iCat)
  {
    fHt[iCat] = TFile::Open(("../TF1_pol1_EtScale/TF1_pol1_"+analysis+Form("_std_cat%d.root",iCat)).c_str(),"READ");
    /*
    TF1* f_std = (TF1*)(fHt[iCat]->Get(Form("HT_cat%d",iCat)) );
    f_std->SetLineColor(kBlue);
    f_std->SetLineWidth(2);

    TProfile* tf1 = (TProfile*)( f->Get(Form("envelop_cat%d",iCat)) );
    tf1 -> SetLineColor(kBlue);
    tf1 -> SetMaximum(1.002);
    tf1 -> SetMinimum(0.998);
    tf1 -> GetXaxis() -> SetRangeUser(65,220);
    tf1 -> GetXaxis() -> SetTitle("H_{T} (GeV)");
    //    tf1 -> GetYaxis() -> SetTitle("events");
    */

    hET[iCat] = (TH1F*)( f->Get(Form("ETSinglePho%d",iCat)) );
    hET[iCat]->SetLineColor(colors.at(iCat));
    hET[iCat]->SetLineWidth(3);
    //    hET -> Rebin(rebin);
    hET[iCat]->Scale(1./hET[iCat]->Integral());
    hET[iCat]->SetMaximum(1.25*hET[iCat]->GetMaximum());
    hET[iCat]->GetXaxis() -> SetRangeUser(0.,500.);
    hET[iCat]->GetXaxis() -> SetTitle("single #gamma E_{T} (GeV)");
    hET[iCat]->GetYaxis() -> SetTitle("events");

    hHT[iCat] = (TH1F*)( f->Get(Form("HTDoublePho%d",iCat)) );
    hHT[iCat]->SetLineColor(colors.at(iCat));
    hHT[iCat]->SetLineWidth(3);
    //    hHT -> Rebin(rebin);
    hET[iCat]->Scale(1./hHT[iCat]->Integral());
    hHT[iCat]->SetMaximum(1.25*hHT[iCat]->GetMaximum());
    //hHT -> GetXaxis() -> SetRangeUser(0.,500.);
    hHT[iCat]->GetXaxis() -> SetTitle("E^{#gamma 1}_{T}+ E^{#gamma 2}_{T} (GeV)");
    hHT[iCat]->GetYaxis() -> SetTitle("events");

    hMCorr[iCat] = (TH1F*)( f->Get(Form("MassCorrection%d",iCat)) );
    hMCorr[iCat] -> SetLineColor(colors.at(iCat));
    hMCorr[iCat] -> SetLineWidth(3);
    hMCorr[iCat] -> Rebin(10*rebin);
    hMCorr[iCat]->Scale(1./hMCorr[iCat]->Integral());
    hMCorr[iCat] -> SetMaximum(1.25*hMCorr[iCat]->GetMaximum());
    hMCorr[iCat] -> GetXaxis() -> SetRangeUser(0.995,1.005);
    hMCorr[iCat] -> GetXaxis() -> SetTitle("mass correction");
    hMCorr[iCat] -> GetYaxis() -> SetTitle("events");

    hSingPho[iCat] = (TH1F*)( f->Get(Form("SingleCat%d",iCat)) );
    hSingPho[iCat]->SetLineColor(colors.at(iCat));
    hSingPho[iCat]->SetLineWidth(3);
    //    hSingPho -> Rebin(rebin);
    hSingPho[iCat]->Scale(1./hSingPho[iCat]->Integral());
    hSingPho[iCat]->SetMaximum(1.25*hSingPho[iCat]->GetMaximum());
    hSingPho[iCat]->GetXaxis() -> SetRangeUser(0.,500.);
    hSingPho[iCat]->GetXaxis() -> SetTitle("single #gamma category");
    hSingPho[iCat]->GetYaxis() -> SetTitle("events");


    TH1F* h1 = (TH1F*)( f->Get(Form("Hgg_original_cat%d",iCat)) );
    h1 -> SetLineColor(kBlue);
    h1 -> SetLineWidth(3);
    h1 -> Rebin(rebin);
    h1 -> SetMaximum(1.25*h1->GetMaximum());
    h1 -> GetXaxis() -> SetRangeUser(110,140);
    h1 -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h1 -> GetYaxis() -> SetTitle("events");
    
    TH1F* h2 = (TH1F*)( f->Get(Form("Hgg_measuredDiff_cat%d",iCat)) );
    //    h2 -> SetLineColor(colors.at(iCat));
    h2 -> SetLineColor(kRed);
    h2 -> SetLineWidth(3);
    h2 -> Rebin(rebin);
    h2 -> SetMaximum(1.25*h2->GetMaximum());    
    h2 -> GetXaxis() -> SetRangeUser(110,140);
    h2 -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h2 -> GetYaxis() -> SetTitle("events");
    
    TH1F* h31 = (TH1F*)( f->Get(Form("Hgg_scaleMinus1_cat%d",iCat)) );
    h31 -> SetLineColor(kGreen+2);
    h31 -> SetLineWidth(3);
    h31 -> Rebin(rebin);
    h31 -> SetMaximum(1.25*h31->GetMaximum());
    h31 -> GetXaxis() -> SetRangeUser(110,140);
    h31 -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h31 -> GetYaxis() -> SetTitle("events");
    
    TH1F* h32 = (TH1F*)( f->Get(Form("Hgg_scalePlus1_cat%d",iCat)) );
    h32 -> SetLineColor(kBlue);
    h32 -> SetLineWidth(3);
    h32 -> Rebin(rebin);
    h32 -> SetMaximum(1.25*h32->GetMaximum());
    h32 -> GetXaxis() -> SetRangeUser(110,140);
    h32 -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h32 -> GetYaxis() -> SetTitle("events");
    
    TH1F* h41 = (TH1F*)( f->Get(Form("Hgg_smearMinus1_cat%d",iCat)) );
    h41 -> SetLineColor(kGreen+2);
    h41 -> SetLineWidth(3);
    h41 -> Rebin(rebin);
    h41 -> SetMaximum(1.25*h41->GetMaximum());
    h41 -> GetXaxis() -> SetRangeUser(110,140);
    h41 -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h41 -> GetYaxis() -> SetTitle("events");
    
    TH1F* h42 = (TH1F*)( f->Get(Form("Hgg_smearPlus1_cat%d",iCat)) );
    h42 -> SetLineColor(kBlue);
    h42 -> SetLineWidth(3);
    h42 -> Rebin(rebin);
    h42 -> SetMaximum(1.25*h42->GetMaximum());
    h42 -> GetXaxis() -> SetRangeUser(110,140);
    h42 -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h42 -> GetYaxis() -> SetTitle("events");
        
    /*
    TCanvas* c_tf1 = new TCanvas(Form("c_tf1_cat%d",iCat),Form("c_tf1_cat%d",iCat));
    tf1->SetFillColor(kAzure+1);
    tf1->SetFillStyle(3001);
    tf1->SetLineColor(kAzure+1);
    tf1->SetMarkerStyle(7);
    tf1->Draw("e3");
    f_std->SetLineWidth(2);
    f_std->Draw("same");
    c_tf1->Print(Form("c_tf1_cat%d.png",iCat),"png"); 
    */
  
  
    TCanvas* c_std = new TCanvas(Form("c_std_cat%d",iCat),Form("c_std_cat%d",iCat));
    c_std -> cd();
    c_std -> SetGridx();
    c_std -> SetGridy();
    
    h1 -> Draw("hist");
    h2 -> Draw("hist,same");    

    h2->GetXaxis()->SetRangeUser(0.,300.);
    h1->GetXaxis()->SetRangeUser(0.,300.);
    latex = new TLatex(0.15,0.90,Form("(mean_{H#gamma#gamma}^{meas. non lin.}/mean_{H#gamma#gamma}^{central sample}-1) = %1.2e",h2->GetMean()/h1->GetMean()-1));
    h2->GetXaxis()->SetRangeUser(110,140);
    h1->GetXaxis()->SetRangeUser(110,140);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(1);
    latex -> Draw("same");
    
    c_std -> Print(Form("c_std_cat%d.png",iCat),"png");
    
    std::cout << " >>>> ciao 0 " << std::endl;    
    
    TCanvas* c_scaleSyst = new TCanvas(Form("c_scaleSyst_cat%d",iCat),Form("c_scaleSyst_cat%d",iCat));
    c_scaleSyst -> cd();
    
    h2 -> Draw("hist");
    h31 -> Draw("hist,same");
    h32 -> Draw("hist,same");
    h2 -> Draw("hist,same");
    
    h31->GetXaxis()->SetRangeUser(0.,300.);
    h2->GetXaxis()->SetRangeUser(0.,300.);
    latex = new TLatex(0.15,0.90,Form("(mean_{H#gamma#gamma}^{scale syst.+1#sigma}/mean_{H#gamma#gamma}^{meas. non lin.}-1) = %1.2e",h31->GetMean()/h2->GetMean()-1));
    h31->GetXaxis()->SetRangeUser(110,140);
    h2->GetXaxis()->SetRangeUser(110,140);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(h31->GetLineColor());
    latex -> Draw("same");
    h32->GetXaxis()->SetRangeUser(0.,300.);
    h2->GetXaxis()->SetRangeUser(0.,300.);
    latex = new TLatex(0.15,0.82,Form("(mean_{H#gamma#gamma}^{scale syst.-1#sigma}/mean_{H#gamma#gamma}^{meas. non lin.}-1) = %1.2e",h32->GetMean()/h2->GetMean()-1));
    h32->GetXaxis()->SetRangeUser(110,140);
    h2->GetXaxis()->SetRangeUser(110,140);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(h32->GetLineColor());
    latex -> Draw("same");
    
    c_scaleSyst -> Print(Form("c_scaleSyst_cat%d.png",iCat),"png");
    
    std::cout << " >>>> ciao 1 " << std::endl;    
    
    TCanvas* c_smearSyst = new TCanvas(Form("c_smearSyst_cat%d",iCat),Form("c_smearSyst_cat%d",iCat));
    c_smearSyst -> cd();
    
    h2 -> Draw("hist");
    h41 -> Draw("hist,same");
    h42 -> Draw("hist,same");
    h2 -> Draw("hist,same");
    
    h41->GetXaxis()->SetRangeUser(0.,300.);
    h2->GetXaxis()->SetRangeUser(0.,300.);
    latex = new TLatex(0.15,0.90,Form("(mean_{H#gamma#gamma}^{smear syst.+1#sigma}/mean_{H#gamma#gamma}^{meas. non lin.}-1) = %1.2e",h41->GetMean()/h2->GetMean()-1));
    h41->GetXaxis()->SetRangeUser(110,140);
    h2->GetXaxis()->SetRangeUser(110,140);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(h31->GetLineColor());
    latex -> Draw("same");
    h42->GetXaxis()->SetRangeUser(0.,300.);
    h2->GetXaxis()->SetRangeUser(0.,300.);
    latex = new TLatex(0.15,0.82,Form("(mean_{H#gamma#gamma}^{smear syst.-1#sigma}/mean_{H#gamma#gamma}^{meas. non lin.}-1) = %1.2e",h42->GetMean()/h2->GetMean()-1));
    h42->GetXaxis()->SetRangeUser(110,140);
    h2->GetXaxis()->SetRangeUser(110,140);
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(h32->GetLineColor());
    latex -> Draw("same");
    
    c_smearSyst -> Print(Form("c_smearSyst_cat%d.png",iCat),"png");
    
    std::cout << " >>>> ciao 2 " << std::endl;    
    
    TCanvas* c = new TCanvas(Form("c_stat_cat%d",iCat),Form("c_stat_cat%d",iCat));
    c -> SetGridx();
    c -> SetGridy();
    c -> cd();

    std::cout << " >>>> ciao 2a " << std::endl;    
    
    //    TH1F* h_distr = new TH1F(Form("h_distr_cat%d",iCat),"",60,-0.003,0.003);
    TH1F* h_distr = new TH1F(Form("h_distr_cat%d",iCat),"",120,-0.006,0.006);
    h_distr -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    h_distr -> GetYaxis() -> SetTitle("a.u.");
    h_distr -> SetLineWidth(2);
    
    h1 -> Draw("hist");
    h2 -> Draw("hist,same");
    
    h2 -> GetXaxis() -> SetRangeUser(0.,300.);
    for(int iTrial = 0; iTrial < nTrials; ++iTrial)
    {
      TH1F* g = (TH1F*)( f->Get(Form("stat/Hgg_measuredDiff_stat_cat%d_trial%d",iCat,iTrial)) );
      g -> Rebin(rebin);
//       if(g->Integral() != h2->Integral()) {
// 	std::cout << " >>>> g->Integral() = " << g->OverFlow()Integral() << std::endl;
// 	std::cout << " >>>> h2->Integral() = " << h2->Integral() << std::endl;
//       }
//      if(g->GetMean() > h2->GetMean()) g->SetLineColor(kGreen);
//       if(g->GetMean() < h2->GetMean()) g->SetLineColor(kCyan);
      g->Draw("hist,same");
//       if(g->GetMean() == h2->GetMean()){
//       std::cout << " g->GetMean() = " << std::setprecision(6) << g->GetMean() << std::endl;
//       std::cout << " h2->GetMean() = " << std::setprecision(6) << h2->GetMean() << std::endl;
//       }
//       std::cout << " g->GetMean() = " << std::scientific << g->GetMean() << std::endl;
//       std::cout << " h2->GetMean() = " << std::scientific << h2->GetMean() << std::endl;
      h_distr -> Fill(g->GetMean()/(1.*h2->GetMean()) - 1.);
    }
    
    std::cout << " >>>> ciao 2b " << std::endl;

    h1 -> Draw("hist,same");    
    //h2 -> Draw("hist,same");

    c -> Print(Form("c_stat_cat%d.png",iCat),"png");

    
    std::cout << " >>>> ciao 3 " << std::endl;    

    TCanvas* c_distr = new TCanvas(Form("c_distr_cat%d",iCat),Form("c_distr_cat%d",iCat));
    c_distr -> SetGridx();
    c_distr -> SetGridy();
    c_distr -> cd();
    //    h_distr -> GetXaxis() -> SetRangeUser(0.975,1.025);
    h_distr -> GetXaxis() -> SetRangeUser(-0.025,0.025);
    h_distr -> Draw("");

    latex = new TLatex(0.15,0.90,Form("RMS: %1.2e",h_distr->GetRMS()));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(1);
    latex -> Draw("same");

    std::cout << " iCat = " << iCat << std::endl;
    latex = new TLatex(0.15,0.80,Form("#Delta: %1.2e",(h_distr->GetMean())) );
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(1);
    latex -> Draw("same");

    c_distr -> Print(Form("c_distr_cat%d.png",iCat),"png");

    std::cout << " >>>> ciao 4 " << std::endl;    

    
    /*
    double mean, meanErr, min, max;
    //    FindSmallestInterval(mean, meanErr,min,max,h_distr,0.68,false);
        
    TH1F* h_distr2 = (TH1F*)( h_distr->Clone("h_distr2"));
    h_distr2 -> Reset();
    for(int bin = 1; bin <= h_distr->GetNbinsX(); ++bin)
    {
      float binCenter = h_distr->GetBinCenter(bin);
      if( binCenter >= min && binCenter <= max ) h_distr2 -> SetBinContent(bin,h_distr->GetBinContent(bin));
    }
    //    h_distr2->SetFillColor(kYellow);
    //    h_distr2->Draw("same");
    */
    

  }

  /*
  if(HggMethod == "CiC"){
    hMCorr[0]->GetYaxis()->SetRangeUser(0.01, 0.1.25*hMCorr[1]->GetMaximum());
    hET[0]->GetYaxis()->SetRangeUser(0.01, 1.25*hET[1]->GetMaximum());
   hHT[0]->GetYaxis()->SetRangeUser(0.01, 1.25*hHT[1]->GetMaximum());
  }
  if(HggMethod == "MVA"){
    hMCorr[0]->GetYaxis()->SetRangeUser(0.01, 1.25*hMCorr[3]->GetMaximum());
    hET[0]->GetYaxis()->SetRangeUser(0.01, 1.25*hET[3]->GetMaximum());
    hHT[0]->GetYaxis()->SetRangeUser(0.01, 1.25*hHT[3]->GetMaximum());
  }
  */

  TCanvas* c_MCorr = new TCanvas("c_MCorr", "c_MCorr");
  hMCorr[0]->DrawNormalized("hist");
  for(int iCat = 1; iCat < nCats; ++iCat)    hMCorr[iCat]->DrawNormalized("hist, same");
  c_MCorr->Print("c_MCorr.png"); 

  TCanvas* c_SingPho = new TCanvas("c_SingPho", "c_SingPho");
  hSingPho[0]->DrawNormalized("hist");
  for(int iCat = 1; iCat < nCats; ++iCat)    hSingPho[iCat]->DrawNormalized("hist, same");
  c_SingPho->Print("c_SingPho.png"); 

  TCanvas* c_ET = new TCanvas("c_ET","c_ET");
  hET[0]->DrawNormalized("hist");
  for(int iCat = 1; iCat < nCats; ++iCat)    hET[iCat]->DrawNormalized("hist, same");
  c_ET->Print("c_ET.png"); 

  TCanvas* c_HT = new TCanvas("c_HT", "c_HT");
  hHT[0]->DrawNormalized("hist");
  for(int iCat = 1; iCat < nCats; ++iCat)    hHT[iCat]->DrawNormalized("hist, same");
  c_HT->Print("c_HT.png"); 


//   TCanvas* c_MCorr_Log = new TCanvas("c_MCorr_Log", "c_MCorr_Log");
//   gPad->SetLogy();
//   hMCorr[3]->Draw("hist");
//   for(int iCat = 0; iCat < nCats-1; ++iCat)    hMCorr[iCat]->Draw("hist, same");
//   c_MCorr_Log->Print("c_MCorr_Log.png"); 

  TCanvas* c_ET_Log = new TCanvas("c_ET_Log", "c_ET_Log");
  gPad->SetLogy();
  hET[0]->DrawNormalized("hist");
  for(int iCat = 1; iCat < nCats; ++iCat)    hET[iCat]->DrawNormalized("hist, same");
  c_ET->Print("c_ET_Log.png"); 

  TCanvas* c_HT_Log = new TCanvas("c_HT_Log", "c_HT_Log");
  gPad->SetLogy();
  hHT[0]->DrawNormalized("hist");
  for(int iCat = 1; iCat < nCats; ++iCat)    hHT[iCat]->DrawNormalized("hist, same");
  c_HT->Print("c_HT_Log.png"); 



}



void FindSmallestInterval(double& mean, double& meanErr, double& min, double& max,
                          TH1F* histo,
                          const double& fraction, const bool& verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;
  
  
  double binWidth = histo->GetBinWidth(1);
  double delta = 999999.;
  double integralMax = fraction * histo->Integral();
  
  int N = histo -> GetNbinsX();
  std::vector<float> binCenters(N);
  std::vector<float> binContents(N);
  std::vector<float> binIntegrals(N);
  for(int bin1 = 0; bin1 < N; ++bin1)
  {
    if( verbosity ) std::cout << "trying bin " << bin1 << " / " << N << "\r" << std::flush; 
    binCenters[bin1] = histo->GetBinCenter(bin1+1);
    binContents[bin1] = histo->GetBinContent(bin1+1);
    
    for(int bin2 = 0; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += binContents[bin2];
  }
  
  for(int bin1 = 0; bin1 < N; ++bin1)
  {
    for(int bin2 = bin1+1; bin2 < N; ++bin2)
    {
      if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;
      
      double tmpMin = histo -> GetBinCenter(bin1+1)-0.5*binWidth;
      double tmpMax = histo -> GetBinCenter(bin2+1)+0.5*binWidth;
      if( tmpMax-tmpMin < delta )
      {
        delta = tmpMax - tmpMin;
        min = tmpMin;
        max = tmpMax;
      }
      
      break;
    }
  }
  
  TH1F* h_temp2 = (TH1F*)( histo->Clone("h_temp2") );
  h_temp2 -> Reset();
  for(int bin = 0; bin < N; ++bin)
    if( (binCenters[bin] >= min) && (binCenters[bin] <= max) )
      h_temp2 -> SetBinContent(bin,binContents[bin]);
  
  mean = h_temp2 -> GetMean();
  meanErr = h_temp2 -> GetMeanError();  
  
  delete h_temp2;
}

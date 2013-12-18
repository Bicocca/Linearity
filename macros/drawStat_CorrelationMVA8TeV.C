#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

void FindSmallestInterval(double& mean, double& meanErr, double& min, double& max,
                          TH1F* histo,
                          const double& fraction, const bool& verbosity);



void drawStat_CorrelationMVA8TeV(std::string analysis, std::string HggMethod)
{
  // 5 x 5 categories
  TFile* f = TFile::Open("rescaledHgg_histos.root","READ");
  TLatex* latex;
  
  //  int rebin = 1;
  //  int rebin = 12;
  int rebin = 6;
  int nCats = 5;
  int nTrials = 1000;

  bool PrintHistos = false;

  TH1F** hDistr_cat0 = new TH1F*[nCats];
  TH1F** hDistr_cat1 = new TH1F*[nCats];
  TH1F** hDistr_cat2 = new TH1F*[nCats];
  TH1F** hDistr_cat3 = new TH1F*[nCats];
  TH1F** hDistr_cat4 = new TH1F*[nCats];

  TH1F** hDM_cat = new TH1F*[nCats];

  TH1F** h1 = new TH1F*[nCats];
  TH1F** h2 = new TH1F*[nCats];

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
    h1[iCat] = (TH1F*)( f->Get(Form("Hgg_original_cat%d",iCat)) );
    h1[iCat] -> SetLineColor(kBlue);
    h1[iCat] -> SetLineWidth(3);
    h1[iCat] -> Rebin(rebin);
    h1[iCat] -> SetMaximum(1.25*h1[iCat]->GetMaximum());
    //    h1[iCat] -> GetXaxis() -> SetRangeUser(110,140);
    h1[iCat] -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h1[iCat] -> GetYaxis() -> SetTitle("events");
     
    h2[iCat] = (TH1F*)( f->Get(Form("Hgg_measuredDiff_cat%d",iCat)) );
    h2[iCat] -> SetLineColor(kRed);
    h2[iCat] -> SetLineWidth(3);
    h2[iCat] -> Rebin(rebin);
    h2[iCat] -> SetMaximum(1.25*h2[iCat]->GetMaximum());    
    //    h2[iCat] -> GetXaxis() -> SetRangeUser(110,140);
    h2[iCat] -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h2[iCat] -> GetYaxis() -> SetTitle("events");
   
    hDistr_cat0[iCat] = new TH1F(Form("h_distr_cat0%d",iCat),"",120,-0.006,0.006);
    //hDistr_cat0[iCat] = new TH1F(Form("h_distr_cat0%d",iCat),"",240,-0.012,0.012);
    hDistr_cat0[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_cat0[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_cat0[iCat] -> SetLineWidth(2);

    hDistr_cat1[iCat] = new TH1F(Form("h_distr_cat1%d",iCat),"",120,-0.006,0.006);
    hDistr_cat1[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_cat1[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_cat1[iCat] -> SetLineWidth(2);

    hDistr_cat2[iCat] = new TH1F(Form("h_distr_cat2%d",iCat),"",120,-0.006,0.006);
    hDistr_cat2[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_cat2[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_cat2[iCat] -> SetLineWidth(2);

    hDistr_cat3[iCat] = new TH1F(Form("h_distr_cat3%d",iCat),"",120,-0.006,0.006);
    hDistr_cat3[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_cat3[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_cat3[iCat] -> SetLineWidth(2);

    hDistr_cat4[iCat] = new TH1F(Form("h_distr_cat4%d",iCat),"",120,-0.006,0.006);
    hDistr_cat4[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_cat4[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_cat4[iCat] -> SetLineWidth(2);

    hDM_cat[iCat] = new TH1F(Form("h_DM_cat%d",iCat),"",120,-0.006,0.006);
    hDM_cat[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDM_cat[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDM_cat[iCat] -> SetLineWidth(2);
  }

  for(int iTrial = 0; iTrial < nTrials; ++iTrial){
      TH1F** g = new TH1F*[nCats];
      for(int iCat = 0; iCat < nCats; ++iCat){      
	g[iCat] = (TH1F*)( f->Get(Form("stat/Hgg_measuredDiff_stat_cat%d_trial%d",iCat, iTrial)) );
      }
      for(int iCat = 0; iCat < nCats; ++iCat){      
      //      for(int iCat = 0; iCat < 1; ++iCat){      
 	hDistr_cat0[iCat]->Fill( (g[0]->GetMean()/(1.*h2[0]->GetMean()) - 1.) * (g[iCat]->GetMean()/(1.*h2[iCat]->GetMean()) - 1.) );
 	hDistr_cat1[iCat]->Fill( (g[1]->GetMean()/(1.*h2[1]->GetMean()) - 1.) * (g[iCat]->GetMean()/(1.*h2[iCat]->GetMean()) - 1.) );
 	hDistr_cat2[iCat]->Fill( (g[2]->GetMean()/(1.*h2[2]->GetMean()) - 1.) * (g[iCat]->GetMean()/(1.*h2[iCat]->GetMean()) - 1.) );
 	hDistr_cat3[iCat]->Fill( (g[3]->GetMean()/(1.*h2[3]->GetMean()) - 1.) * (g[iCat]->GetMean()/(1.*h2[iCat]->GetMean()) - 1.) );
 	hDistr_cat4[iCat]->Fill( (g[4]->GetMean()/(1.*h2[4]->GetMean()) - 1.) * (g[iCat]->GetMean()/(1.*h2[iCat]->GetMean()) - 1.) );
 	hDM_cat[iCat]->Fill( (g[iCat]->GetMean()/(1.*h2[iCat]->GetMean()) - 1.) );
      }
      for(int iCat = 0; iCat < nCats; ++iCat)  delete g[iCat];
    }


    float cov00 = sqrt(hDistr_cat0[0]->GetMean() - (hDM_cat[0]->GetMean() * hDM_cat[0]->GetMean()) );
    float cov01 = sqrt(hDistr_cat0[1]->GetMean() - (hDM_cat[0]->GetMean() * hDM_cat[1]->GetMean()) );
    float cov02 = sqrt(hDistr_cat0[2]->GetMean() - (hDM_cat[0]->GetMean() * hDM_cat[2]->GetMean()) );
    float cov03 = sqrt(hDistr_cat0[3]->GetMean() - (hDM_cat[0]->GetMean() * hDM_cat[3]->GetMean()) );
    float cov04 = sqrt(hDistr_cat0[4]->GetMean() - (hDM_cat[0]->GetMean() * hDM_cat[4]->GetMean()) );

    float cov10 = sqrt(hDistr_cat1[0]->GetMean() - (hDM_cat[1]->GetMean() * hDM_cat[0]->GetMean()) );
    float cov11 = sqrt(hDistr_cat1[1]->GetMean() - (hDM_cat[1]->GetMean() * hDM_cat[1]->GetMean()) );
    float cov12 = sqrt(hDistr_cat1[2]->GetMean() - (hDM_cat[1]->GetMean() * hDM_cat[2]->GetMean()) );
    float cov13 = sqrt(hDistr_cat1[3]->GetMean() - (hDM_cat[1]->GetMean() * hDM_cat[3]->GetMean()) );
    float cov14 = sqrt(hDistr_cat1[4]->GetMean() - (hDM_cat[1]->GetMean() * hDM_cat[4]->GetMean()) );

    float cov20 = sqrt(hDistr_cat2[0]->GetMean() - (hDM_cat[2]->GetMean() * hDM_cat[0]->GetMean()) );
    float cov21 = sqrt(hDistr_cat2[1]->GetMean() - (hDM_cat[2]->GetMean() * hDM_cat[1]->GetMean()) );
    float cov22 = sqrt(hDistr_cat2[2]->GetMean() - (hDM_cat[2]->GetMean() * hDM_cat[2]->GetMean()) );
    float cov23 = sqrt(hDistr_cat2[3]->GetMean() - (hDM_cat[2]->GetMean() * hDM_cat[3]->GetMean()) );
    float cov24 = sqrt(hDistr_cat2[4]->GetMean() - (hDM_cat[2]->GetMean() * hDM_cat[4]->GetMean()) );

    float cov30 = sqrt(hDistr_cat3[0]->GetMean() - (hDM_cat[3]->GetMean() * hDM_cat[0]->GetMean()) );
    float cov31 = sqrt(hDistr_cat3[1]->GetMean() - (hDM_cat[3]->GetMean() * hDM_cat[1]->GetMean()) );
    float cov32 = sqrt(hDistr_cat3[2]->GetMean() - (hDM_cat[3]->GetMean() * hDM_cat[2]->GetMean()) );
    float cov33 = sqrt(hDistr_cat3[3]->GetMean() - (hDM_cat[3]->GetMean() * hDM_cat[3]->GetMean()) );
    float cov34 = sqrt(hDistr_cat3[4]->GetMean() - (hDM_cat[3]->GetMean() * hDM_cat[4]->GetMean()) );

    float cov40 = sqrt(hDistr_cat4[0]->GetMean() - (hDM_cat[4]->GetMean() * hDM_cat[0]->GetMean()) );
    float cov41 = sqrt(hDistr_cat4[1]->GetMean() - (hDM_cat[4]->GetMean() * hDM_cat[1]->GetMean()) );
    float cov42 = sqrt(hDistr_cat4[2]->GetMean() - (hDM_cat[4]->GetMean() * hDM_cat[2]->GetMean()) );
    float cov43 = sqrt(hDistr_cat4[3]->GetMean() - (hDM_cat[4]->GetMean() * hDM_cat[3]->GetMean()) );
    float cov44 = sqrt(hDistr_cat4[4]->GetMean() - (hDM_cat[4]->GetMean() * hDM_cat[4]->GetMean()) );

    std::cout << " Cov 00 = " << cov00 << std::endl;
    std::cout << " Cor 00 = " << cov00 / sqrt(cov00*cov00) << std::endl;

    std::cout << " Cov 01 = " << cov01 << std::endl;
    std::cout << " Cor 01 = " << cov01 / sqrt(cov00*cov11) << std::endl;

    std::cout << " Cov 02 = " << cov02 << std::endl;
    std::cout << " Cor 02 = " << cov02 / sqrt(cov00*cov22) << std::endl;

    std::cout << " Cov 03 = " << cov03 << std::endl;
    std::cout << " Cor 03 = " << cov03 / sqrt(cov00*cov33) << std::endl;

    std::cout << " Cov 04 = " << cov04 << std::endl;
    std::cout << " Cor 04 = " << cov04 / sqrt(cov00*cov44) << std::endl;

    ////////////
    std::cout << " Cov 10 = " << cov10 << std::endl;
    std::cout << " Cor 10 = " << cov10 / sqrt(cov11*cov00) << std::endl;

    std::cout << " Cov 11 = " << cov11 << std::endl;
    std::cout << " Cor 11 = " << cov11 / sqrt(cov11*cov11) << std::endl;

    std::cout << " Cov 12 = " << cov12 << std::endl;
    std::cout << " Cor 12 = " << cov12 / sqrt(cov11*cov22) << std::endl;

    std::cout << " Cov 13 = " << cov13 << std::endl;
    std::cout << " Cor 13 = " << cov13 / sqrt(cov11*cov33) << std::endl;

    std::cout << " Cov 14 = " << cov14 << std::endl;
    std::cout << " Cor 14 = " << cov14 / sqrt(cov11*cov44) << std::endl;

    ////////////
    std::cout << " Cov 20 = " << cov20 << std::endl;
    std::cout << " Cor 20 = " << cov20 / sqrt(cov22*cov00) << std::endl;

    std::cout << " Cov 21 = " << cov21 << std::endl;
    std::cout << " Cor 21 = " << cov21 / sqrt(cov22*cov11) << std::endl;

    std::cout << " Cov 22 = " << cov22 << std::endl;
    std::cout << " Cor 22 = " << cov22 / sqrt(cov22*cov22) << std::endl;

    std::cout << " Cov 23 = " << cov23 << std::endl;
    std::cout << " Cor 23 = " << cov23 / sqrt(cov22*cov33) << std::endl;

    std::cout << " Cov 24 = " << cov24 << std::endl;
    std::cout << " Cor 24 = " << cov24 / sqrt(cov22*cov44) << std::endl;

    ////////////
    std::cout << " Cov 30 = " << cov30 << std::endl;
    std::cout << " Cor 30 = " << cov30 / sqrt(cov33*cov00) << std::endl;

    std::cout << " Cov 31 = " << cov31 << std::endl;
    std::cout << " Cor 31 = " << cov31 / sqrt(cov33*cov11) << std::endl;

    std::cout << " Cov 32 = " << cov32 << std::endl;
    std::cout << " Cor 32 = " << cov32 / sqrt(cov33*cov22) << std::endl;

    std::cout << " Cov 33 = " << cov33 << std::endl;
    std::cout << " Cor 33 = " << cov33 / sqrt(cov33*cov33) << std::endl;

    std::cout << " Cov 34 = " << cov34 << std::endl;
    std::cout << " Cor 34 = " << cov34 / sqrt(cov33*cov44) << std::endl;

    ////////////
    std::cout << " Cov 40 = " << cov40 << std::endl;
    std::cout << " Cor 40 = " << cov40 / sqrt(cov44*cov00) << std::endl;

    std::cout << " Cov 41 = " << cov41 << std::endl;
    std::cout << " Cor 41 = " << cov41 / sqrt(cov44*cov11) << std::endl;

    std::cout << " Cov 42 = " << cov42 << std::endl;
    std::cout << " Cor 42 = " << cov42 / sqrt(cov44*cov22) << std::endl;

    std::cout << " Cov 43 = " << cov43 << std::endl;
    std::cout << " Cor 43 = " << cov43 / sqrt(cov44*cov33) << std::endl;

    std::cout << " Cov 44 = " << cov44 << std::endl;
    std::cout << " Cor 44 = " << cov44 / sqrt(cov44*cov44) << std::endl;


  if(PrintHistos){
     TCanvas* c0 = new TCanvas(Form("c_stat_cat0%d",iCat),Form("c_stat_cat0%d",iCat));
     c0 -> SetGridx();
     c0 -> SetGridy();
     c0 -> cd();
     hDistr_cat0[iCat]-> GetXaxis() -> SetRangeUser(-0.025,0.025);
     hDistr_cat0[iCat]-> Draw("");

     latex = new TLatex(0.15,0.90,Form("RMS: %1.2e",hDistr_cat0[iCat]->GetRMS()));
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     latex = new TLatex(0.15,0.80,Form("#Delta: %1.2e",(hDistr_cat0[iCat]->GetMean())) );
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     std::cout << " RMS = " << hDistr_cat0[iCat]->GetRMS() << " #Delta = " << hDistr_cat0[iCat]->GetMean() 
	       << " Quad Sum = " << sqrt(pow(hDistr_cat0[iCat]->GetRMS(), 2) + pow(hDistr_cat0[iCat]->GetMean(), 2)) << std::endl; 
     c0->Print(Form("c_stat_cat0%d.png",iCat),"png");

     /////////
     TCanvas* c1 = new TCanvas(Form("c_stat_cat1%d",iCat),Form("c_stat_cat1%d",iCat));
     c1 -> SetGridx();
     c1 -> SetGridy();
     c1 -> cd();
     hDistr_cat1[iCat]-> GetXaxis() -> SetRangeUser(-0.025,0.025);
     hDistr_cat1[iCat]-> Draw("");

     latex = new TLatex(0.15,0.90,Form("RMS: %1.2e",hDistr_cat1[iCat]->GetRMS()));
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     latex = new TLatex(0.15,0.80,Form("#Delta: %1.2e",(hDistr_cat1[iCat]->GetMean())) );
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     std::cout << " RMS = " << hDistr_cat1[iCat]->GetRMS() << " #Delta = " << hDistr_cat1[iCat]->GetMean() 
	       << " Quad Sum = " << sqrt(pow(hDistr_cat1[iCat]->GetRMS(), 2) + pow(hDistr_cat1[iCat]->GetMean(), 2)) << std::endl; 
     c1->Print(Form("c_stat_cat1%d.png",iCat),"png");
     ////////////
     TCanvas* c2 = new TCanvas(Form("c_stat_cat2%d",iCat),Form("c_stat_cat2%d",iCat));
     c2 -> SetGridx();
     c2 -> SetGridy();
     c2 -> cd();
     hDistr_cat2[iCat]-> GetXaxis() -> SetRangeUser(-0.025,0.025);
     hDistr_cat2[iCat]-> Draw("");

     latex = new TLatex(0.15,0.90,Form("RMS: %1.2e",hDistr_cat2[iCat]->GetRMS()));
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     latex = new TLatex(0.15,0.80,Form("#Delta: %1.2e",(hDistr_cat2[iCat]->GetMean())) );
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     std::cout << " RMS = " << hDistr_cat2[iCat]->GetRMS() << " #Delta = " << hDistr_cat2[iCat]->GetMean() 
	       << " Quad Sum = " << sqrt(pow(hDistr_cat2[iCat]->GetRMS(), 2) + pow(hDistr_cat2[iCat]->GetMean(), 2)) << std::endl; 
     c2->Print(Form("c_stat_cat2%d.png",iCat),"png");
     ////////////////////
     TCanvas* c3 = new TCanvas(Form("c_stat_cat3%d",iCat),Form("c_stat_cat3%d",iCat));
     c3 -> SetGridx();
     c3 -> SetGridy();
     c3 -> cd();
     hDistr_cat3[iCat]-> GetXaxis() -> SetRangeUser(-0.025,0.025);
     hDistr_cat3[iCat]-> Draw("");

     latex = new TLatex(0.15,0.90,Form("RMS: %1.2e",hDistr_cat3[iCat]->GetRMS()));
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     latex = new TLatex(0.15,0.80,Form("#Delta: %1.2e",(hDistr_cat3[iCat]->GetMean())) );
     latex -> SetNDC();
     latex -> SetTextFont(42);
     latex -> SetTextSize(0.04);
     latex -> SetTextColor(1);
     latex -> Draw("same");

     std::cout << " RMS = " << hDistr_cat3[iCat]->GetRMS() << " #Delta = " << hDistr_cat3[iCat]->GetMean() 
	       << " Quad Sum = " << sqrt(pow(hDistr_cat3[iCat]->GetRMS(), 2) + pow(hDistr_cat3[iCat]->GetMean(), 2)) << std::endl; 
     c3->Print(Form("c_stat_cat3%d.png",iCat),"png");
  }
    
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

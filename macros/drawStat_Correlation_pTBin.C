#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

void FindSmallestInterval(double& mean, double& meanErr, double& min, double& max,
                          TH1F* histo,
                          const double& fraction, const bool& verbosity);



void drawStat_Correlation_pTBin(std::string analysis, std::string HggMethod)
{
  // 8 x 8 categories
  TFile* f = TFile::Open("rescaledHgg_histos.root","READ");
  TLatex* latex;
  
  //  int rebin = 1;
  //  int rebin = 12;
  int rebin = 6;
  int nCats = 4;
  int nTrials = 1000;

  bool PrintHistos = false;

  TH1F** hDistr_Lcat0L = new TH1F*[nCats];
  TH1F** hDistr_Lcat0H = new TH1F*[nCats];
  TH1F** hDistr_Hcat0L = new TH1F*[nCats];
  TH1F** hDistr_Hcat0H = new TH1F*[nCats];
  TH1F** hDistr_Lcat1L = new TH1F*[nCats];
  TH1F** hDistr_Lcat1H = new TH1F*[nCats];
  TH1F** hDistr_Hcat1L = new TH1F*[nCats];
  TH1F** hDistr_Hcat1H = new TH1F*[nCats];
  TH1F** hDistr_Hcat2L = new TH1F*[nCats];
  TH1F** hDistr_Hcat2H = new TH1F*[nCats];
  TH1F** hDistr_Lcat2H = new TH1F*[nCats];
  TH1F** hDistr_Lcat2L = new TH1F*[nCats];
  TH1F** hDistr_Hcat3H = new TH1F*[nCats];
  TH1F** hDistr_Hcat3L = new TH1F*[nCats];
  TH1F** hDistr_Lcat3L = new TH1F*[nCats];
  TH1F** hDistr_Lcat3H = new TH1F*[nCats];

  TH1F** hDM_Lcat = new TH1F*[nCats];
  TH1F** hDM_Hcat = new TH1F*[nCats];

  TH1F** h1L = new TH1F*[nCats];
  TH1F** h1H = new TH1F*[nCats];
  TH1F** h2L = new TH1F*[nCats];
  TH1F** h2H = new TH1F*[nCats];

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
    h1L[iCat] = (TH1F*)( f->Get(Form("Hgg_original_pTLow_cat%d",iCat)) );
    h1L[iCat] -> SetLineColor(kBlue);
    h1L[iCat] -> SetLineWidth(3);
    h1L[iCat] -> Rebin(rebin);
    h1L[iCat] -> SetMaximum(1.25*h1L[iCat]->GetMaximum());
    h1L[iCat] -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h1L[iCat] -> GetYaxis() -> SetTitle("events");

    h1H[iCat] = (TH1F*)( f->Get(Form("Hgg_original_pTHigh_cat%d",iCat)) );
    h1H[iCat] -> SetLineColor(kBlue);
    h1H[iCat] -> SetLineWidth(3);
    h1H[iCat] -> Rebin(rebin);
    h1H[iCat] -> SetMaximum(1.25*h1H[iCat]->GetMaximum());
    h1H[iCat] -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h1H[iCat] -> GetYaxis() -> SetTitle("events");
     
    h2L[iCat] = (TH1F*)( f->Get(Form("Hgg_measuredDiff_pTLow_cat%d",iCat)) );
    h2L[iCat] -> SetLineColor(kRed);
    h2L[iCat] -> SetLineWidth(3);
    h2L[iCat] -> Rebin(rebin);
    h2L[iCat] -> SetMaximum(1.25*h2L[iCat]->GetMaximum());    
    h2L[iCat] -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h2L[iCat] -> GetYaxis() -> SetTitle("events");

    h2H[iCat] = (TH1F*)( f->Get(Form("Hgg_measuredDiff_pTHigh_cat%d",iCat)) );
    h2H[iCat] -> SetLineColor(kRed);
    h2H[iCat] -> SetLineWidth(3);
    h2H[iCat] -> Rebin(rebin);
    h2H[iCat] -> SetMaximum(1.25*h2H[iCat]->GetMaximum());    
    h2H[iCat] -> GetXaxis() -> SetTitle("m_{#gamma#gamma} (GeV)");
    h2H[iCat] -> GetYaxis() -> SetTitle("events");
   
    hDistr_Lcat0L[iCat] = new TH1F(Form("h_distr_Lcat0L%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat0L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat0L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat0L[iCat] -> SetLineWidth(2);
    //
    hDistr_Lcat0H[iCat] = new TH1F(Form("h_distr_Lcat0H%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat0H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat0H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat0H[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat0L[iCat] = new TH1F(Form("h_distr_Hcat0L%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat0L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat0L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat0L[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat0H[iCat] = new TH1F(Form("h_distr_Hcat0H%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat0H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat0H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat0H[iCat] -> SetLineWidth(2);

    hDistr_Lcat1L[iCat] = new TH1F(Form("h_distr_Lcat1L%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat1L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat1L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat1L[iCat] -> SetLineWidth(2);
    //
    hDistr_Lcat1H[iCat] = new TH1F(Form("h_distr_Lcat1H%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat1H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat1H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat1H[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat1L[iCat] = new TH1F(Form("h_distr_Hcat1L%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat1L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat1L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat1L[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat1H[iCat] = new TH1F(Form("h_distr_Hcat1H%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat1H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat1H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat1H[iCat] -> SetLineWidth(2);

    hDistr_Lcat2L[iCat] = new TH1F(Form("h_distr_Lcat2L%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat2L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat2L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat2L[iCat] -> SetLineWidth(2);
    //
    hDistr_Lcat2H[iCat] = new TH1F(Form("h_distr_Lcat2H%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat2H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat2H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat2H[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat2L[iCat] = new TH1F(Form("h_distr_Hcat2L%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat2L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat2L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat2L[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat2H[iCat] = new TH1F(Form("h_distr_Hcat2H%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat2H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat2H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat2H[iCat] -> SetLineWidth(2);

    hDistr_Lcat3L[iCat] = new TH1F(Form("h_distr_LcatL3%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat3L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat3L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat3L[iCat] -> SetLineWidth(2);
    //
    hDistr_Lcat3H[iCat] = new TH1F(Form("h_distr_LcatH3%d",iCat),"",120,-0.006,0.006);
    hDistr_Lcat3H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Lcat3H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Lcat3H[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat3L[iCat] = new TH1F(Form("h_distr_HcatL3%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat3L[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat3L[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat3L[iCat] -> SetLineWidth(2);
    //
    hDistr_Hcat3H[iCat] = new TH1F(Form("h_distr_HcatH3%d",iCat),"",120,-0.006,0.006);
    hDistr_Hcat3H[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDistr_Hcat3H[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDistr_Hcat3H[iCat] -> SetLineWidth(2);

    hDM_Lcat[iCat] = new TH1F(Form("h_DM_Lcat%d",iCat),"",120,-0.006,0.006);
    hDM_Lcat[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDM_Lcat[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDM_Lcat[iCat] -> SetLineWidth(2);

    hDM_Hcat[iCat] = new TH1F(Form("h_DM_Hcat%d",iCat),"",120,-0.006,0.006);
    hDM_Hcat[iCat] -> GetXaxis() -> SetTitle("mean_{H#gamma#gamma}^{stat. err. prop.} / mean_{H#gamma#gamma}^{central sample}");
    hDM_Hcat[iCat] -> GetYaxis() -> SetTitle("a.u.");
    hDM_Hcat[iCat] -> SetLineWidth(2);
  }

  for(int iTrial = 0; iTrial < nTrials; ++iTrial){
      TH1F** gL = new TH1F*[nCats];
      TH1F** gH = new TH1F*[nCats];
      for(int iCat = 0; iCat < nCats; ++iCat){      
	gL[iCat] = (TH1F*)( f->Get(Form("stat/Hgg_measuredDiff_stat_pTLow_cat%d_trial%d",iCat, iTrial)) );
	gH[iCat] = (TH1F*)( f->Get(Form("stat/Hgg_measuredDiff_stat_pTHigh_cat%d_trial%d",iCat, iTrial)) );
      }
      for(int iCat = 0; iCat < nCats; ++iCat){      
 	hDistr_Lcat0L[iCat]->Fill( (gL[0]->GetMean()/(1.*h2L[0]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Lcat0H[iCat]->Fill( (gL[0]->GetMean()/(1.*h2L[0]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat0L[iCat]->Fill( (gH[0]->GetMean()/(1.*h2H[0]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat0H[iCat]->Fill( (gH[0]->GetMean()/(1.*h2H[0]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );

 	hDistr_Lcat1L[iCat]->Fill( (gL[1]->GetMean()/(1.*h2L[1]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Lcat1H[iCat]->Fill( (gL[1]->GetMean()/(1.*h2L[1]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat1L[iCat]->Fill( (gH[1]->GetMean()/(1.*h2H[1]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat1H[iCat]->Fill( (gH[1]->GetMean()/(1.*h2H[1]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );

 	hDistr_Lcat2L[iCat]->Fill( (gL[2]->GetMean()/(1.*h2L[2]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Lcat2H[iCat]->Fill( (gL[2]->GetMean()/(1.*h2L[2]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat2L[iCat]->Fill( (gH[2]->GetMean()/(1.*h2H[2]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat2H[iCat]->Fill( (gH[2]->GetMean()/(1.*h2H[2]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );

 	hDistr_Lcat3L[iCat]->Fill( (gL[3]->GetMean()/(1.*h2L[3]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Lcat3H[iCat]->Fill( (gL[3]->GetMean()/(1.*h2L[3]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat3L[iCat]->Fill( (gH[3]->GetMean()/(1.*h2H[3]->GetMean()) - 1.) * (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDistr_Hcat3H[iCat]->Fill( (gH[3]->GetMean()/(1.*h2H[3]->GetMean()) - 1.) * (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );

 	hDM_Lcat[iCat]->Fill( (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDM_Lcat[iCat]->Fill( (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDM_Lcat[iCat]->Fill( (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );
 	hDM_Lcat[iCat]->Fill( (gL[iCat]->GetMean()/(1.*h2L[iCat]->GetMean()) - 1.) );

 	hDM_Hcat[iCat]->Fill( (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDM_Hcat[iCat]->Fill( (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDM_Hcat[iCat]->Fill( (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
 	hDM_Hcat[iCat]->Fill( (gH[iCat]->GetMean()/(1.*h2H[iCat]->GetMean()) - 1.) );
      }
      for(int iCat = 0; iCat < nCats; ++iCat)  {delete gL[iCat]; delete gH[iCat];}
    }


    float cov0L0L = sqrt(hDistr_Lcat0L[0]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov0L0H = sqrt(hDistr_Lcat0H[0]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Hcat[0]->GetMean()) );
    float cov0H0L = sqrt(hDistr_Hcat0L[0]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov0H0H = sqrt(hDistr_Hcat0H[0]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Hcat[0]->GetMean()) );
    //
    float cov0L1L = sqrt(hDistr_Lcat0L[1]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov0L1H = sqrt(hDistr_Lcat0H[1]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Hcat[1]->GetMean()) );
    float cov0H1L = sqrt(hDistr_Hcat0L[1]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov0H1H = sqrt(hDistr_Hcat0H[1]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Hcat[1]->GetMean()) );
    //
    float cov0L2L = sqrt(hDistr_Lcat0L[2]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov0L2H = sqrt(hDistr_Lcat0H[2]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Hcat[2]->GetMean()) );
    float cov0H2L = sqrt(hDistr_Hcat0L[2]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov0H2H = sqrt(hDistr_Hcat0H[2]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Hcat[2]->GetMean()) );
    //
    float cov0L3L = sqrt(hDistr_Lcat0L[3]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov0L3H = sqrt(hDistr_Lcat0H[3]->GetMean() - (hDM_Lcat[0]->GetMean() * hDM_Hcat[3]->GetMean()) );
    float cov0H3L = sqrt(hDistr_Hcat0L[3]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov0H3L = sqrt(hDistr_Hcat0H[3]->GetMean() - (hDM_Hcat[0]->GetMean() * hDM_Hcat[3]->GetMean()) );

    float cov1L0L = sqrt(hDistr_Lcat1L[0]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov1L0H = sqrt(hDistr_Lcat1H[0]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Hcat[0]->GetMean()) );
    float cov1H0L = sqrt(hDistr_Hcat1L[0]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov1H0H = sqrt(hDistr_Hcat1H[0]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Hcat[0]->GetMean()) );
    //
    float cov1L1L = sqrt(hDistr_Lcat1L[1]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov1L1H = sqrt(hDistr_Lcat1H[1]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Hcat[1]->GetMean()) );
    float cov1H1L = sqrt(hDistr_Hcat1L[1]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov1H1H = sqrt(hDistr_Hcat1H[1]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Hcat[1]->GetMean()) );
    //
    float cov1L2L = sqrt(hDistr_Lcat1L[2]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov1L2H = sqrt(hDistr_Lcat1H[2]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Hcat[2]->GetMean()) );
    float cov1H2L = sqrt(hDistr_Hcat1L[2]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov1H2H = sqrt(hDistr_Hcat1H[2]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Hcat[2]->GetMean()) );
    //
    float cov1L3L = sqrt(hDistr_Lcat1L[3]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov1L3H = sqrt(hDistr_Lcat1H[3]->GetMean() - (hDM_Lcat[1]->GetMean() * hDM_Hcat[3]->GetMean()) );
    float cov1H3L = sqrt(hDistr_Hcat1L[3]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov1H3L = sqrt(hDistr_Hcat1H[3]->GetMean() - (hDM_Hcat[1]->GetMean() * hDM_Hcat[3]->GetMean()) );

    float cov2L0L = sqrt(hDistr_Lcat2L[0]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov2L0H = sqrt(hDistr_Lcat2H[0]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Hcat[0]->GetMean()) );
    float cov2H0L = sqrt(hDistr_Hcat2L[0]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov2H0H = sqrt(hDistr_Hcat2H[0]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Hcat[0]->GetMean()) );
    //
    float cov2L1L = sqrt(hDistr_Lcat2L[1]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov2L1H = sqrt(hDistr_Lcat2H[1]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Hcat[1]->GetMean()) );
    float cov2H1L = sqrt(hDistr_Hcat2L[1]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov2H1H = sqrt(hDistr_Hcat2H[1]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Hcat[1]->GetMean()) );
    //
    float cov2L2L = sqrt(hDistr_Lcat2L[2]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov2L2H = sqrt(hDistr_Lcat2H[2]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Hcat[2]->GetMean()) );
    float cov2H2L = sqrt(hDistr_Hcat2L[2]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov2H2H = sqrt(hDistr_Hcat2H[2]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Hcat[2]->GetMean()) );
    //
    float cov2L3L = sqrt(hDistr_Lcat2L[3]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov2L3H = sqrt(hDistr_Lcat2H[3]->GetMean() - (hDM_Lcat[2]->GetMean() * hDM_Hcat[3]->GetMean()) );
    float cov2H3L = sqrt(hDistr_Hcat2L[3]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov2H3L = sqrt(hDistr_Hcat2H[3]->GetMean() - (hDM_Hcat[2]->GetMean() * hDM_Hcat[3]->GetMean()) );

    float cov3L0L = sqrt(hDistr_Lcat3L[0]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov3L0H = sqrt(hDistr_Lcat3H[0]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Hcat[0]->GetMean()) );
    float cov3H0L = sqrt(hDistr_Hcat3L[0]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Lcat[0]->GetMean()) );
    float cov3H0H = sqrt(hDistr_Hcat3H[0]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Hcat[0]->GetMean()) );
    //
    float cov3L1L = sqrt(hDistr_Lcat3L[1]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov3L1H = sqrt(hDistr_Lcat3H[1]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Hcat[1]->GetMean()) );
    float cov3H1L = sqrt(hDistr_Hcat3L[1]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Lcat[1]->GetMean()) );
    float cov3H1H = sqrt(hDistr_Hcat3H[1]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Hcat[1]->GetMean()) );
    //
    float cov3L2L = sqrt(hDistr_Lcat3L[2]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov3L2H = sqrt(hDistr_Lcat3H[2]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Hcat[2]->GetMean()) );
    float cov3H2L = sqrt(hDistr_Hcat3L[2]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Lcat[2]->GetMean()) );
    float cov3H2H = sqrt(hDistr_Hcat3H[2]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Hcat[2]->GetMean()) );
    //
    float cov3L3L = sqrt(hDistr_Lcat3L[3]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov3L3H = sqrt(hDistr_Lcat3H[3]->GetMean() - (hDM_Lcat[3]->GetMean() * hDM_Hcat[3]->GetMean()) );
    float cov3H3L = sqrt(hDistr_Hcat3L[3]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Lcat[3]->GetMean()) );
    float cov3H3H = sqrt(hDistr_Hcat3H[3]->GetMean() - (hDM_Hcat[3]->GetMean() * hDM_Hcat[3]->GetMean()) );
//     std::cout << " hDistr_Hcat3H[3]->GetMean() = " << hDistr_Hcat3H[3]->GetMean() << std::endl;
//     std::cout << " hDM_Hcat[3]->GetMean() = " << hDM_Hcat[3]->GetMean() << std::endl;

    std::cout << " Cor 0L 0L = " << cov0L0L / sqrt(cov0L0L*cov0L0L) << std::endl;
    std::cout << " Cor 0L 0H = " << cov0L0H / sqrt(cov0L0H*cov0L0H) << std::endl;
    std::cout << " Cor 0H 0L = " << cov0H0L / sqrt(cov0H0L*cov0H0L) << std::endl;
    std::cout << " Cor 0H 0H = " << cov0H0H / sqrt(cov0H0H*cov0H0H) << std::endl;

    std::cout << " Cor 0L 1L = " << cov0L1L / sqrt(cov0L0L*cov1L1L) << std::endl;
    std::cout << " Cor 0L 1H = " << cov0L1H / sqrt(cov0L0L*cov1H1H) << std::endl;
    std::cout << " Cor 0H 1L = " << cov0H1L / sqrt(cov0H0H*cov1L1L) << std::endl;
    std::cout << " Cor 0H 1H = " << cov0H1H / sqrt(cov0H0H*cov1H1H) << std::endl;

    std::cout << " Cor 0L 2L = " << cov0L2L / sqrt(cov0L0L*cov2L2L) << std::endl;
    std::cout << " Cor 0L 2H = " << cov0L2H / sqrt(cov0L0L*cov2H2H) << std::endl;
    std::cout << " Cor 0H 2L = " << cov0H2L / sqrt(cov0H0H*cov2L2L) << std::endl;
    std::cout << " Cor 0H 2H = " << cov0H2H / sqrt(cov0H0H*cov2H2H) << std::endl;

    std::cout << " Cor 0L 3L = " << cov0L3L / sqrt(cov0L0L*cov3L3L) << std::endl;
    //    std::cout << " cov3H3H = " << cov3H3H << std::endl;
    std::cout << " Cor 0L 3H = " << cov0L3H / sqrt(cov0L0L*cov3H3H) << std::endl;
    std::cout << " Cor 0H 3L = " << cov0L3L / sqrt(cov0H0H*cov3L3L) << std::endl;
    std::cout << " Cor 0H 3H = " << cov0L3L / sqrt(cov0H0H*cov3H3H) << std::endl;

    ////////////

    std::cout << " Cor 1L 0L = " << cov1L0L / sqrt(cov1L1L*cov0L0L) << std::endl;
    std::cout << " Cor 1L 0H = " << cov1L0H / sqrt(cov1L1L*cov0H0H) << std::endl;
    std::cout << " Cor 1H 0L = " << cov1H0L / sqrt(cov1H1H*cov0L0L) << std::endl;
    std::cout << " Cor 1H 0H = " << cov1H0H / sqrt(cov1H1H*cov0H0H) << std::endl;

    std::cout << " Cor 1L 1L = " << cov1L1L / sqrt(cov1L1L*cov1L1L) << std::endl;
    std::cout << " Cor 1L 1H = " << cov1L1H / sqrt(cov1L1L*cov1H1H) << std::endl;
    std::cout << " Cor 1H 1L = " << cov1H1L / sqrt(cov1H1H*cov1L1L) << std::endl;
    std::cout << " Cor 1H 1H = " << cov1H1H / sqrt(cov1H1H*cov1H1H) << std::endl;

    std::cout << " Cor 1L 2L = " << cov1L2L / sqrt(cov1L1L*cov2L2L) << std::endl;
    std::cout << " Cor 1L 2H = " << cov1L2H / sqrt(cov1L1L*cov2H2H) << std::endl;
    std::cout << " Cor 1H 2L = " << cov1H2L / sqrt(cov1H1H*cov2L2L) << std::endl;
    std::cout << " Cor 1H 2H = " << cov1H2H / sqrt(cov1H1H*cov2H2H) << std::endl;

    std::cout << " Cor 1L 3L = " << cov1L3L / sqrt(cov1L1L*cov3L3L) << std::endl;
    std::cout << " Cor 1L 3H = " << cov1L3H / sqrt(cov1L1L*cov3H3H) << std::endl;
    std::cout << " Cor 1H 3L = " << cov1L3L / sqrt(cov1H1H*cov3L3L) << std::endl;
    std::cout << " Cor 1H 3H = " << cov1L3L / sqrt(cov1H1H*cov3H3H) << std::endl;
    ////////////
    std::cout << " Cor 2L 0L = " << cov2L0L / sqrt(cov2L2L*cov0L0L) << std::endl;
    std::cout << " Cor 2L 0H = " << cov2L0H / sqrt(cov2L2L*cov0H0H) << std::endl;
    std::cout << " Cor 2H 0L = " << cov2H0L / sqrt(cov2H2H*cov0L0L) << std::endl;
    std::cout << " Cor 2H 0H = " << cov2H0H / sqrt(cov2H2H*cov0H0H) << std::endl;

    std::cout << " Cor 2L 1L = " << cov2L1L / sqrt(cov2L2L*cov1L1L) << std::endl;
    std::cout << " Cor 2L 1H = " << cov2L1H / sqrt(cov2L2L*cov1H1H) << std::endl;
    std::cout << " Cor 2H 1L = " << cov2H1L / sqrt(cov2H2H*cov1L1L) << std::endl;
    std::cout << " Cor 2H 1H = " << cov2H1H / sqrt(cov2H2H*cov1H1H) << std::endl;

    std::cout << " Cor 2L 2L = " << cov2L2L / sqrt(cov2L2L*cov2L2L) << std::endl;
    std::cout << " Cor 2L 2H = " << cov2L2H / sqrt(cov2L2L*cov2H2H) << std::endl;
    std::cout << " Cor 2H 2L = " << cov2H2L / sqrt(cov2H2H*cov2L2L) << std::endl;
    std::cout << " Cor 2H 2H = " << cov2H2H / sqrt(cov2H2H*cov2H2H) << std::endl;

    std::cout << " Cor 2L 3L = " << cov2L3L / sqrt(cov2L2L*cov3L3L) << std::endl;
    std::cout << " Cor 2L 3H = " << cov2L3H / sqrt(cov2L2L*cov3H3H) << std::endl;
    std::cout << " Cor 2H 3L = " << cov2L3L / sqrt(cov2H2H*cov3L3L) << std::endl;
    std::cout << " Cor 2H 3H = " << cov2L3L / sqrt(cov2H2H*cov3H3H) << std::endl;

    ////////////
    std::cout << " Cor 3L 0L = " << cov3L0L / sqrt(cov3L3L*cov0L0L) << std::endl;
    std::cout << " Cor 3L 0H = " << cov3L0H / sqrt(cov3L3L*cov0H0H) << std::endl;
    std::cout << " Cor 3H 0L = " << cov3H0L / sqrt(cov3H3H*cov0L0L) << std::endl;
    std::cout << " Cor 3H 0H = " << cov3H0H / sqrt(cov3H3H*cov0H0H) << std::endl;

    std::cout << " Cor 3L 1L = " << cov3L1L / sqrt(cov3L3L*cov1L1L) << std::endl;
    std::cout << " Cor 3L 1H = " << cov3L1H / sqrt(cov3L3L*cov1H1H) << std::endl;
    std::cout << " Cor 3H 1L = " << cov3H1L / sqrt(cov3H3H*cov1L1L) << std::endl;
    std::cout << " Cor 3H 1H = " << cov3H1H / sqrt(cov3H3H*cov1H1H) << std::endl;

    std::cout << " Cor 3L 2L = " << cov3L2L / sqrt(cov3L3L*cov2L2L) << std::endl;
    std::cout << " Cor 3L 2H = " << cov3L2H / sqrt(cov3L3L*cov2H2H) << std::endl;
    std::cout << " Cor 3H 2L = " << cov3H2L / sqrt(cov3H3H*cov2L2L) << std::endl;
    std::cout << " Cor 3H 2H = " << cov3H2H / sqrt(cov3H3H*cov2H2H) << std::endl;

    std::cout << " Cor 3L 3L = " << cov3L3L / sqrt(cov3L3L*cov3L3L) << std::endl;
    std::cout << " Cor 3L 3H = " << cov3L3H / sqrt(cov3L3L*cov3H3H) << std::endl;
    std::cout << " Cor 3H 3L = " << cov3L3L / sqrt(cov3H3H*cov3L3L) << std::endl;
    std::cout << " Cor 3H 3H = " << cov3L3L / sqrt(cov3H3H*cov3H3H) << std::endl;

    ////////////


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

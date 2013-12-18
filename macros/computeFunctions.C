void computeFunctions()
{

  int nCats = 4;
  //  std::string analysis = "CiC";
  std::string analysis = "stdCat";

  TFile** fHt_std = new TFile*[nCats];
  TFile** fHt_ScaleP1 = new TFile*[nCats];
  TFile** fHt_ScaleM1 = new TFile*[nCats];

  for(int iCat = 0; iCat < nCats; ++iCat)
    {
      fHt_std[iCat] = TFile::Open(("TF1_pol1_EtScale/TF1_pol1_"+analysis+Form("_std_cat%d.root",iCat)).c_str(),"READ");
      fHt_ScaleP1[iCat] = TFile::Open(("TF1_pol1_EtScale/TF1_pol1_"+analysis+Form("_scalePlus1_cat%d.root",iCat)).c_str(),"READ");
      fHt_ScaleM1[iCat] = TFile::Open(("TF1_pol1_EtScale/TF1_pol1_"+analysis+Form("_scaleMinus1_cat%d.root",iCat)).c_str(),"READ");

      TF1* f_std = (TF1*)(fHt_std[iCat]->Get(Form("HT_cat%d",iCat)) );
      f_std->SetName("f_std");
      TF1* f_ScaleP1 = (TF1*)(fHt_ScaleP1[iCat]->Get(Form("HT_cat%d",iCat)) );
      f_ScaleP1->SetName("f_ScaleP1");
      TF1* f_ScaleM1 = (TF1*)(fHt_ScaleM1[iCat]->Get(Form("HT_cat%d",iCat)) );
      f_std->SetLineColor(kRed);
      f_ScaleM1->SetName("f_ScaleM1");
      f_std->SetLineWidth(2);

      f_ScaleP1->SetLineColor(kGreen+2);
      f_ScaleP1->SetLineWidth(2);
      f_ScaleM1->SetLineColor(kBlue);
      f_ScaleM1->SetLineWidth(2);

      TF1* f_ScaleP1_d = new TF1("f_ScaleP1_d", "f_ScaleP1 - [0] ", 20., 1000.);
      f_ScaleP1_d->FixParameter(0, f_ScaleP1->Eval(90.) - f_std->Eval(90) );
      f_ScaleP1_d->SetName(Form("HT_cat%d",iCat));
      std::cout << " cat " << iCat << std::endl;
      std::cout << " f_ScaleP1->Eval(90) " << f_ScaleP1->Eval(90) << std::endl;
      std::cout << " f_std->Eval(90) " << f_std->Eval(90) << std::endl;
      std::cout << " f_ScaleP1_d->Eval(90) " << f_ScaleP1_d->Eval(90) << std::endl;
      std::cout << " dovrebbe " << f_ScaleP1->Eval(90) - (f_ScaleP1->Eval(90) - f_std->Eval(90)) << std::endl;

      std::cout << " cat " << iCat << std::endl;
      std::cout << " cat " << iCat << std::endl;

      TF1* f_ScaleP1_d_ET = new TF1(Form("ET_cat%d", iCat), "[0]+ [1] * (x-45.)", 20., 1000.);
      f_ScaleP1_d_ET->SetParameter(0, f_ScaleP1_d->Eval(90.));
      f_ScaleP1_d_ET->SetParError(0, f_ScaleP1_d->GetParError(0));
      f_ScaleP1_d_ET->SetParameter(1, 2.* (f_ScaleP1_d->Eval(91.) - f_ScaleP1_d->Eval(90.) ) );
      f_ScaleP1_d_ET->SetParError(1, 2. * f_ScaleP1_d->GetParError(1));
      f_ScaleP1_d_ET->SetName(Form("ET_cat%d",iCat));
//       f_ScaleP1_d_ET->SetParameter(0, f_ScaleP1_d->GetParameter(0));
//       f_ScaleP1_d_ET->SetParError(0, f_ScaleP1_d->GetParError(0));
//       f_ScaleP1_d_ET->SetParameter(1, f_ScaleP1_d->GetParameter(1));
//       f_ScaleP1_d_ET->SetParError(1, f_ScaleP1_d->GetParError(1));

//       std::cout << " f_ScaleP1_d->GetParameter(0) = " << f_ScaleP1_d->GetParameter(0) << std::endl;
//       std::cout << " f_ScaleP1_d_ET->GetParameter(0) = " << f_ScaleP1_d_ET->GetParameter(0) << std::endl;
//       std::cout << " f_ScaleP1_d->GetParameter(1) = " << f_ScaleP1_d->GetParameter(1) << std::endl;
//       std::cout << " f_ScaleP1_d_ET->GetParameter(1) = " << f_ScaleP1_d_ET->GetParameter(1) << std::endl;
//       std::cout << " f_ScaleP1_d->Eval(90) " << f_ScaleP1_d->Eval(90) << std::endl;
//       std::cout << " f_ScaleP1_d_ET->Eval(90) " << f_ScaleP1_d_ET->Eval(90) << std::endl;
//       std::cout << " f_ScaleP1_d->Eval(45) " << f_ScaleP1_d->Eval(45) << std::endl;
//       std::cout << " f_ScaleP1_d_ET->Eval(45) " << f_ScaleP1_d_ET->Eval(45) << std::endl;


      TF1* f_ScaleM1_d = new TF1("f_ScaleM1_d", "f_ScaleM1 - [0]", 20., 1000.);
      f_ScaleM1_d->SetParameter(0, f_ScaleM1->Eval(90.)- f_std->Eval(90.));
      f_ScaleM1_d->SetName(Form("HT_cat%d",iCat));

      std::cout << " f_ScaleM1->Eval(90) " << f_ScaleM1->Eval(90) << std::endl;
      std::cout << " f_std->Eval(90) " << f_std->Eval(90) << std::endl;
      std::cout << " f_ScaleM1_d->Eval(90) " << f_ScaleM1_d->Eval(90) << std::endl;

      TF1* f_ScaleM1_d_ET = new TF1(Form("ET_cat%d", iCat), "[0]+ [1] * (x-45.)", 20., 1000.);
      f_ScaleM1_d_ET->SetParameter(0, f_ScaleM1_d->Eval(90.));
      f_ScaleM1_d_ET->SetParError(0, f_ScaleM1_d->GetParError(0));
      f_ScaleM1_d_ET->SetParameter(1, 2. * (f_ScaleM1_d->Eval(91.) - f_ScaleM1_d->Eval(90.) ) );
      f_ScaleM1_d_ET->SetParError(1, 2. * f_ScaleM1_d->GetParError(1));
      f_ScaleM1_d_ET->SetName(Form("ET_cat%d",iCat));

      TFile fHt_ScaleP1_n(("TF1_pol1_"+analysis+Form("_scalePlus1Diff_cat%d.root",iCat)).c_str(),"recreate");
      f_ScaleP1_d->Write();
      f_ScaleP1_d_ET->Write();
      fHt_ScaleP1_n.Close();

      TFile fHt_ScaleM1_n(("TF1_pol1_"+analysis+Form("_scaleMinus1Diff_cat%d.root",iCat)).c_str(),"recreate");
      f_ScaleM1_d->Write();
      f_ScaleM1_d_ET->Write();
      fHt_ScaleM1_n.Close();

      f_ScaleP1_d->SetLineColor(kGreen+2);
      f_ScaleP1_d->SetLineWidth(2);

      f_ScaleM1_d->SetLineColor(kBlue);
      f_ScaleM1_d->SetLineWidth(2);

      f_ScaleP1_d->SetLineStyle(2);
      f_ScaleM1_d->SetLineStyle(2);

      f_ScaleP1_d_ET->SetLineColor(kRed+2);
      f_ScaleM1_d_ET->SetLineColor(kGreen+2);

      f_ScaleP1_d_ET->SetLineWidth(2);
      f_ScaleM1_d_ET->SetLineWidth(2);
      
      TCanvas* c1 = new TCanvas(Form("c_cat%d",iCat),Form("c_cat%d",iCat));
//       f_ScaleP1_d_ET->GetXaxis()->SetRangeUser(0., 200.);         
//       f_ScaleP1_d_ET->GetYaxis()->SetRangeUser(0.995, 1.005);     
       f_std->GetYaxis()->SetRangeUser(0.998, 1.002);
       f_std->GetXaxis()->SetRangeUser(0., 200.);
       f_std->Draw();
       f_ScaleP1->Draw("same");
       f_ScaleM1->Draw("same");
       f_ScaleP1_d->Draw("same");
       f_ScaleM1_d->Draw("same");
//        f_ScaleP1_d_ET->Draw("same");
//        f_ScaleM1_d_ET->Draw("same");
       c1->Print((analysis+Form("c_cat%d.png",iCat)).c_str(), "png");


      //      std::cout << " fHt_std " << std::endl;
      delete f_std;
      delete f_ScaleP1;
      delete f_ScaleM1;
      delete f_ScaleP1_d;
      delete f_ScaleM1_d;

      delete f_ScaleP1_d_ET;
      delete f_ScaleM1_d_ET;


    }
}

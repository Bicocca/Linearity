{
  /*
  TFile* EoP_cat0 = TFile::Open("../EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/studyLinearity_EoP_cat0_-1evtsPerPoint.root");
  TFile* EoP_cat1 = TFile::Open("../EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/studyLinearity_EoP_cat1_-1evtsPerPoint.root");
  TFile* EoP_cat2 = TFile::Open("../EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/studyLinearity_EoP_cat2_-1evtsPerPoint.root");
  TFile* EoP_cat3 = TFile::Open("../EoP_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15/studyLinearity_EoP_cat3_-1evtsPerPoint.root");


  TGraphAsymmErrors** g = new TGraphAsymmErrors*[4];
  char graphName[50];
  sprintf(graphName,"step1/scale_recursiveMean_DAOverMC");
  g[0] = (TGraphAsymmErrors*)( EoP_cat0->Get(graphName) );
  g[1] = (TGraphAsymmErrors*)( EoP_cat1->Get(graphName) );
  g[2] = (TGraphAsymmErrors*)( EoP_cat2->Get(graphName) );
  g[3] = (TGraphAsymmErrors*)( EoP_cat3->Get(graphName) );

  TFile EoPCat4("EoP.root","recreate");
  for(int iCat = 0; iCat < 4; ++iCat){
    g[iCat]->Write(Form("graph_cat%d",iCat));
  }
  EoPCat4.Close();
*/

  TFile* MZ_cat0 = TFile::Open("AddSysOK/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_fewBins/studyLinearity_MZ_cat0_-1evtsPerPoint.root");
  TFile* MZ_cat1 = TFile::Open("AddSysOK/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_fewBins/studyLinearity_MZ_cat1_-1evtsPerPoint.root");
  TFile* MZ_cat2 = TFile::Open("AddSysOK/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_fewBins/studyLinearity_MZ_cat2_-1evtsPerPoint.root");
  TFile* MZ_cat3 = TFile::Open("AddSysOK/MZ_stdCat_nonGlobe_powheg-runDependent_phoTunedRegV5_Dphi3p15_fewBins/studyLinearity_MZ_cat3_-1evtsPerPoint.root");


  TGraphAsymmErrors** g2 = new TGraphAsymmErrors*[4];
  char graphName[50];
  sprintf(graphName,"step1/scale_recursiveMean_DAOverMC");
  g2[0] = (TGraphAsymmErrors*)( MZ_cat0->Get(graphName) );
  g2[1] = (TGraphAsymmErrors*)( MZ_cat1->Get(graphName) );
  g2[2] = (TGraphAsymmErrors*)( MZ_cat2->Get(graphName) );
  g2[3] = (TGraphAsymmErrors*)( MZ_cat3->Get(graphName) );

  TFile MZCat4("MZ_Etdep_AllStepinOne_noRescale.root","recreate");
  for(int iCat = 0; iCat < 4; ++iCat){
    g2[iCat]->Write(Form("graph_cat%d",iCat));
  }
  MZCat4.Close();

}

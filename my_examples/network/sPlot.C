void sPlot(){
	TCanvas *c1=new TCanvas("c1","aaa",900,700);	
	c1->SetLogy();
  //TGraph *gr=new TGraph("result/txt/r-abundance.txt","%*lg %lg %lg %*lg %*lg");
  
    ifstream fdat;
  fdat.open("result/txt/r-abundance.txt");
  Double_t nmass[141];
  Double_t stan[141];
  Double_t min[141];
  Double_t max[141];
  Double_t temp2,fmin,fmax;
  for (Int_t i=0;i<141;i++){
    fdat>>temp2>>nmass[i]>>stan[i]>>fmin>>fmax;
    min[i]=stan[i]-fmin;
    max[i]=fmax-stan[i];
  }
  TGraphAsymmErrors* grs=new TGraphAsymmErrors(141,nmass,stan,0,0,min,max);
  grs->SetMarkerStyle(20);
  grs->SetFillColor(0);
  grs->Draw("AP");
  grs->GetYaxis()->SetRangeUser(1e-13,20);
  gPad->SetLogy();
  TString gr1s="result/txt/NSM/nsm_o_y025.txt";
  TString gr2s="result/txt/NSM/nsm_o_y031.txt";
  TString gr3s="result/txt/NSM/nsm_o_y032.txt";
  TString gr4s="result/txt/NSM/nsm_o_y033.txt";
  TString gr5s="result/txt/NSM/nsm_o_y035.txt";
  TString gr6s="result/txt/NSM/nsm_o_y037.txt";
  
  TGraph *gr1=new TGraph(gr1s,"%lg %lg");
  gr1->SetLineColor(2);
  TGraph *gr2=new TGraph(gr2s,"%lg %lg");
  gr2->SetLineColor(3);
  TGraph *gr3=new TGraph(gr3s,"%lg %lg");
  gr3->SetLineColor(4);
  TGraph *gr4=new TGraph(gr4s,"%lg %lg");
  gr4->SetLineColor(5);
  TGraph *gr5=new TGraph(gr5s,"%lg %lg");
  gr5->SetLineColor(6);
  TGraph *gr6=new TGraph(gr6s,"%lg %lg");
  gr6->SetLineColor(7);
  //Draw
  //gr->GetYaxis()->SetRangeUser(1e-13,20);
  //gr->Draw("APL");  
  gr1->Draw("SAME");
  gr2->Draw("SAME");
  gr3->Draw("SAME");
  gr4->Draw("SAME");
  gr5->Draw("SAME");
  gr6->Draw("SAME"); 
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  //gStyle->SetFillColor(0);
  leg->AddEntry(gr1, gr1s);
  leg->AddEntry(gr2, gr2s);
  leg->AddEntry(gr3, gr3s);
  leg->AddEntry(gr4, gr4s);
  leg->AddEntry(gr5, gr5s);
  leg->AddEntry(gr6, gr6s);
  
  leg->Draw();
}
void sPlot2(){
	TCanvas *c2=new TCanvas("c2","aaa2",900,700);	
	c2->SetLogy();
  TGraph *gr=new TGraph("result/txt/r-abundance_z.txt","%lg %lg %*lg *lg");
  //TGraph *gr1=new TGraph("result/txt/ktuy_s50_ye30.txt","%lg %lg");
  TGraph *gr1=new TGraph("result/txt/special_o_s175_newFRDM_maxAg_z.txt","%lg %lg");
  gr1->SetLineColor(2);
  TGraph *gr2=new TGraph("result/txt/special_o_s175_newFRDM_z.txt","%lg %lg");
  gr2->SetLineColor(3);
  TGraph *gr3=new TGraph("result/txt/special_o_s175_newFRDM.txt","%lg %lg");
  gr3->SetLineColor(4);
  TGraph *gr4=new TGraph("result/txt/special_o_s195_newFRDM.txt","%lg %lg");
  gr4->SetLineColor(5);
  TGraph *gr5=new TGraph("result/txt/special_o_s115_newFRDM.txt","%lg %lg");
  gr5->SetLineColor(6);
  TGraph *gr6=new TGraph("result/txt/ktuy_s40_ye30.txt","%lg %lg");
  gr6->SetLineColor(7);
  //Draw
  gr->GetYaxis()->SetRangeUser(1e-13,20);
  gr->Draw("APL");  
  gr1->Draw("SAME");
  gr2->Draw("SAME");
  //gr3->Draw("SAME");
  //gr4->Draw("SAME");
  //gr5->Draw("SAME");
  //gr6->Draw("SAME"); 
}
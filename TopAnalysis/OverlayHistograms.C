void OverlayHistograms() {
  TFile *Top = TFile::Open("$CMSSW_BASE/src/TopLJets2015/TopAnalysis/testsel_2017.root");
  //TH1F* mlb = (TH1F*)Top->Get("inc_mlb");
  //TH1F* nvtx = (TH1F*)Top->Get("inc_nvtx");
  TH1F* nprotons_high_nvtx = (TH1F*)Top->Get("inc_nprotons_high_nvtx");
  TH1F* nprotons_low_nvtx = (TH1F*)Top->Get("inc_nprotons_low_nvtx");
  //mlb->SetLineColor(kOrange);
  //nvtx->SetLineColor(kOrange+4);
  nprotons_high_nvtx->SetLineColor(kBlue-3);
  nprotons_low_nvtx->SetLineColor(kRed+1);
  //mlb->Draw();
  //nvtx->Draw("same");
  nprotons_high_nvtx->Draw();
  nprotons_low_nvtx->Draw("same");

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetHeader("Number of protons for events of high and low nvtx","C");
  legend->AddEntry(nprotons_low_nvtx, "# of protons in events of nvtx < 27","f");
  legend->AddEntry(nprotons_high_nvtx, "# of protons in events of nvtx >= 27","f");
  legend->Draw();
  return;
}

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/TopSummer2019.h"
#include "TopLJets2015/TopAnalysis/interface/L1PrefireEfficiencyWrapper.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cmath>

#include "TMath.h"

using namespace std;


//
void RunTopSummer2019(const TString in_fname,
                      TString outname,
                      TH1F *normH, 
                      TH1F *genPU, 
                      TString era,
                      Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////
  
  bool isLowPUrun(false);
  if(in_fname.Contains("2017H")) {
    isLowPUrun=true;
    cout << "Running with low PU run settings" << endl;
  }

  //preselection cuts to apply
  float minLeptonPt( isLowPUrun ? 20. : 27);
  size_t minJetMultiplicity(4);

  //CORRECTIONS: LUMINOSITY+PILEUP
  LumiTools lumi(era,genPU);
  std::map<Int_t,Float_t> lumiPerRun=lumi.lumiPerRun();
  
  //CORRECTIONS: LEPTON EFFICIENCIES
  std::map<TString,TString> cfgMap;
  cfgMap["g_id"]="MVAwp90";
  cfgMap["m_id"]="TightID";
  cfgMap["m_iso"]="TightRelIso";
  cfgMap["m_id4iso"]="TightIDandIPCut";
  cfgMap["e_id"]="MVA90";
  EfficiencyScaleFactorsWrapper lepEffH(in_fname.Contains("Data13TeV"),era,cfgMap);

  //CORRECTIONS: L1-prefire 
  L1PrefireEfficiencyWrapper l1PrefireWR(in_fname.Contains("Data13TeV"),era);
  
  //CORRECTIONS: B-TAG CALIBRATION
  BTagSFUtil btvSF(era,BTagEntry::OperatingPoint::OP_MEDIUM,"",0);

  //PREPARE OUTPUT (BOOK SOME HISTOGRAMS)
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  HistTool ht;
  ht.setNsyst(0);
  ht.addHist("nvtx",         new TH1F("nvtx",        ";Vertex multiplicity;Events",50,0,100));
  ht.addHist("mlb",          new TH1F("mlb",         ";m(l,b) [GeV];Events",20,0,250));
  ht.addHist("mlnjets",      new TH1F("mlnjets",     ";m(l,neutrino,jets) [GeV];Events",30,0,2500)); //invariant mass of lepton, neutrino, and jets
  ht.addHist("nprotons",     new TH1F("nprotons",    ";Proton multiplicity; Events",6,0,6) );
  ht.addHist("nprotons_new", new TH1F("nprotons_new",";Proton multiplicity; Events",6,0,6) );
  ht.addHist("npdiff0",      new TH1F("npdiff0",     ";Proton multiplicity; Events",10,-5,5) );
  ht.addHist("npdiff1",      new TH1F("npdiff1",     ";Proton multiplicity; Events",10,-5,5) );
  ht.addHist("csi",          new TH1F("csi",         ";#xi = #deltap/p; Events",50,0,0.3) );
  ht.addHist("x",            new TH1F("x",           ";x  [cm]; Events",50,0,25) );
  ht.addHist("ratevsrun",    new TH1F("ratevsrun",   ";Run number; #sigma [pb]",int(lumiPerRun.size()),0,float(lumiPerRun.size())));
  //CM energy of ttbar system minus lost proton energy
  ht.addHist("Ecentral_minus_Eprotons", new TH1F("Ecentral_minus_Eprotons",";difference [TeV]; Events",50,-1,1));
  ht.addHist("Ecentral_minus_Eprotons_no_neutrino", new TH1F("Ecentral_minus_Eprotons_no_neutrino",";difference [TeV]; Events",50,-1,1));

  //Background and signal plots
  TH1F *Ecentral_minus_Eprotons_bg = new TH1F("Ecentral_minus_Eprotons_bg",";difference [TeV]; Events",50,-1,1);
  TH1F *signal_minus_bg = new TH1F("signal", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);

  //CM minus lost sliced at different CM energies
  TH1F *Ecentral_minus_Eprotons_00_02_TeV = new TH1F("Ecentral_minus_Eprotons_00_02_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_02_04_TeV = new TH1F("Ecentral_minus_Eprotons_02_04_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_04_06_TeV = new TH1F("Ecentral_minus_Eprotons_04_06_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_06_08_TeV = new TH1F("Ecentral_minus_Eprotons_06_08_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_08_10_TeV = new TH1F("Ecentral_minus_Eprotons_08_10_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_10_12_TeV = new TH1F("Ecentral_minus_Eprotons_10_12_TeV",";difference [TeV]; Events",50,-1,1);

  //Background sliced at different CM energies
  TH1F *Ecentral_minus_Eprotons_bg_00_02_TeV = new TH1F("Ecentral_minus_Eprotons_bg_00_02_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_bg_02_04_TeV = new TH1F("Ecentral_minus_Eprotons_bg_02_04_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_bg_04_06_TeV = new TH1F("Ecentral_minus_Eprotons_bg_04_06_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_bg_06_08_TeV = new TH1F("Ecentral_minus_Eprotons_bg_06_08_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_bg_08_10_TeV = new TH1F("Ecentral_minus_Eprotons_bg_08_10_TeV",";difference [TeV]; Events",50,-1,1);
  TH1F *Ecentral_minus_Eprotons_bg_10_12_TeV = new TH1F("Ecentral_minus_Eprotons_bg_10_12_TeV",";difference [TeV]; Events",50,-1,1);

  //Signal plot sliced at different CM energies
  TH1F *signal_minus_bg_00_02_TeV = new TH1F("signal_00_02_TeV", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);
  TH1F *signal_minus_bg_02_04_TeV = new TH1F("signal_02_04_TeV", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);
  TH1F *signal_minus_bg_04_06_TeV = new TH1F("signal_04_06_TeV", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);
  TH1F *signal_minus_bg_06_08_TeV = new TH1F("signal_06_08_TeV", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);
  TH1F *signal_minus_bg_08_10_TeV = new TH1F("signal_08_10_TeV", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);
  TH1F *signal_minus_bg_10_12_TeV = new TH1F("signal_10_12_TeV", ";CM energy - lost proton energy [TeV]; Events(data) - Events(bg)",50,-1,1);

  //2D histo's of CoM vs proton energy for selection and backhground
  TH2F *protons_vs_CM_energy = new TH2F("Eprotons_vs_Ecentral", "Eprotons_vs_Ecentral;CoM_energy;Proton_loss_energy", 50,0,1.2,50,0,1.2);
  protons_vs_CM_energy->SetMarkerStyle(kMultiply);
  protons_vs_CM_energy->SetMarkerColor(9);
  
  TH2F *protons_vs_CM_energy_bg = new TH2F("Eprotons_vs_Ecentral_bg", "Eprotons_vs_Ecentral_bg;CoM_energy;Proton_loss_energy", 50,0,1.2,50,0,1.2);
  protons_vs_CM_energy_bg->SetMarkerStyle(kMultiply);
  protons_vs_CM_energy_bg->SetMarkerColor(8);

  //linear line for 2D histo's
  Double_t x[100], y[100];
  Int_t n = 24;
  for (Int_t i=0;i<n;i++){
    x[i] = i*0.05;
    y[i] = i*0.05;
  }
  TGraph *linear_line = new TGraph(n,x,y);
  int i=0;
  for(auto key : lumiPerRun) {
    i++;
    ht.getPlots()["ratevsrun"]->GetXaxis()->SetBinLabel(i,Form("%d",key.first));
  }

  //INPUT
  MiniEvent_t ev;
  TFile *f = TFile::Open(in_fname);
  if(f==NULL || f->IsZombie()) {
    cout << "Corrupted or missing file " << in_fname << endl;
    return;
  }

  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = min(100000,nentries); //restrict number of entries for testing
  t->GetEntry(0);
  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //EVENT SELECTION WRAPPER (GETS LISTS OF PHYSICS OBJECTS FROM THE INPUT)
  SelectionTool selector(in_fname, false, triggerList);

  //declare variables
  //CM energy of leptons + jet system
  float mlnjets(0);
  std::vector<float> mlnjets_vect{};
  //energy lost by protons
  std::vector<float> lost_proton_energy_vect{};
  std::vector<float> rand_proton_energy_vect{};
  //# of events with CM energy in various ranges
  float nevents_00_02_TeV(0);
  float nevents_02_04_TeV(0);
  float nevents_04_06_TeV(0);
  float nevents_06_08_TeV(0);
  float nevents_08_10_TeV(0);
  float nevents_10_12_TeV(0);
  float nevents_00_02_bg_TeV(0);
  float nevents_02_04_bg_TeV(0);
  float nevents_04_06_bg_TeV(0);
  float nevents_06_08_bg_TeV(0);
  float nevents_08_10_bg_TeV(0);  
  float nevents_10_12_bg_TeV(0);

  //Random protons vector
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      //if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }

      //roman pots
      int nprotons23(0), nprotons123(0);
      int nprotons03(0), nprotons103(0);
      int ntrks( isLowPUrun ? ev.nppstrk : ev.nfwdtrk );
      for (int ift=0; ift<ntrks; ift++) {

        //single pot reconstruction
        if(!isLowPUrun && ev.fwdtrk_method[ift]!=0) continue;

        //only near (pixels) detectors
        const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123 && pot_raw_id!=03 && pot_raw_id!=103) continue;            
        //count nr of protons in each pot  
        nprotons23 += (pot_raw_id==23);
        nprotons123 += (pot_raw_id==123);
        nprotons03 += (pot_raw_id==03);
        nprotons103 += (pot_raw_id==103);
      }
 
      //proton energy loss
      float nrp23(0);
      float nrp123(0);
      //loop over forward trackers to store the number of protons in the
      //two RP:s we're considering (TODO: do this in an earlier loop)
      for (int ift=0; ift<ntrks; ift++) {
        //only near (pixels) detectors
	const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123) continue;
	
	nrp23 += (pot_raw_id==23);
	nrp123 += (pot_raw_id==123);
	
      }
      
      //select events where one proton is captured by each RP
      if (nrp23!=1) continue;
      if (nrp123!=1) continue;
      //store xi for each RP
      float xi_23(0);
      float xi_123(0);
      for (int ift=0; ift<ntrks; ift++) {
	
        //only near (pixels) detectors
        const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123) continue;
	if (pot_raw_id==23) {
	  xi_23 = (isLowPUrun ? 0. : ev.fwdtrk_xi[ift]);
	}
	if (pot_raw_id==123) {
	  xi_123 = (isLowPUrun ? 0. : ev.fwdtrk_xi[ift]);
        }
	
      }
      //calculate proton energy according to P.Meiring's Eq. (9)
      //assuming 13 TeV collisions. Unit: TeV
      float rand_proton_energy = sqrt(13*xi_23*xi_123);
      rand_proton_energy_vect.push_back(rand_proton_energy);
      
    }

  //EVENT LOOP
  //select mu+>=4 jets events triggered by a single muon trigger 
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      if(iev%1000==0) { printf("\r [%3.0f%%] done", 100.*(float)iev/(float)nentries); fflush(stdout); }
      
      //trigger
      bool hasMTrigger(false);
      if(era.Contains("2016")) hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu24_v", ev.triggerBits) );     
      if(era.Contains("2017")) {
        if(isLowPUrun) hasMTrigger=(selector.hasTriggerBit("HLT_HIMu12_v",  ev.addTriggerBits) );   
        else           hasMTrigger=(selector.hasTriggerBit("HLT_IsoMu27_v", ev.triggerBits) );   
      }
      if(!hasMTrigger) continue;

      //select one offline muon
      std::vector<Particle> leptons = selector.flaggedLeptons(ev);     
      leptons = selector.selLeptons(leptons,SelectionTool::TIGHT,SelectionTool::MVA90,minLeptonPt,2.1);
      if(leptons.size()!=1) continue;
      if(leptons[0].id()!=13) continue;

      //select jets
      btvSF.addBTagDecisions(ev);
      if(!ev.isData) btvSF.updateBTagDecisions(ev);      
      std::vector<Jet> allJets = selector.getGoodJets(ev,30.,2.4,leptons,{});
      if(allJets.size()<minJetMultiplicity) continue;

      //met
      TLorentzVector met(0,0,0,0);
      met.SetPtEtaPhiM(ev.met_pt,0,ev.met_phi,0.);

      //mlm = missing longitudinal momentum
      TLorentzVector mlm(0,0,0,0);
      mlm.SetPxPyPzE(0,0,-leptons[0].Pz(),0);

      //me = missing energy (neutrino energy-momentum 4vector)
      TLorentzVector me = met + mlm; 

      //event weight
      float evWgt(1.0);
      
      //data specific: check event rates after selection
      if(ev.isData){
        std::map<Int_t,Float_t>::iterator rIt=lumiPerRun.find(ev.run);
        if(rIt!=lumiPerRun.end()){
          int runBin=std::distance(lumiPerRun.begin(),rIt);
          float lumi=1./rIt->second;
          ht.fill("ratevsrun",runBin,lumi,"inc");
        }else{
          cout << "[Warning] Unable to find run=" << ev.run << endl;
        }
      }

      //MC specific: compute event weight
      if (!ev.isData) {

        float normWgt(normH? normH->GetBinContent(1) : 1.0);        
        TString period = lumi.assignRunPeriod();
        double puWgt(lumi.pileupWeight(ev.g_pu,period)[0]);
        EffCorrection_t selSF = lepEffH.getOfflineCorrection(leptons[0], period);
        EffCorrection_t l1prefireProb=l1PrefireWR.getCorrection(allJets,{});

        evWgt  = normWgt*puWgt*selSF.first*l1prefireProb.first;
        evWgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);        
      }
      

      ht.fill("nvtx",       ev.nvtx,        evWgt, "inc");

      //calculate invariant mass of the system
      TLorentzVector lnjets = leptons[0]+me;
      //prepare variable for no neutrino plot
      float mlnjets_no_neutrino(0);
      for(size_t ij=0; ij<allJets.size(); ij++)
	{
	  lnjets+=allJets[ij];
	  mlnjets = lnjets.M();
	  mlnjets_no_neutrino = (lnjets-me).M();
	}
      ht.fill("mlnjets",mlnjets,evWgt,"invariant_mass");
      mlnjets_vect.push_back(mlnjets);

      //lepton-b systems
      for(size_t ij=0; ij<allJets.size(); ij++) 
        {
          int idx=allJets[ij].getJetIndex();
          bool passBtag(ev.j_btag[idx]>0);
          if(!passBtag) continue;

          float mlb( (leptons[0]+allJets[ij]).M() );
          std::vector<TString> tags={"inc",leptons[0].charge()>0 ? "plus" : "minus"};
          ht.fill("mlb",mlb,evWgt,tags);
        }

      //roman pots
      int nprotons23(0), nprotons123(0);
      int nprotons03(0), nprotons103(0);
      int npdiff0(0), npdiff1(0);
      int ntrks( isLowPUrun ? ev.nppstrk : ev.nfwdtrk );
      int v_low_cutoff(25); //change to mean of nvtx.
      for (int ift=0; ift<ntrks; ift++) {

        //single pot reconstruction
        if(!isLowPUrun && ev.fwdtrk_method[ift]!=0) continue;

        //only near (pixels) detectors
        const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123 && pot_raw_id!=03 && pot_raw_id!=103) continue;            
        //count nr of protons in each pot  
        nprotons23 += (pot_raw_id==23);
        nprotons123 += (pot_raw_id==123);
        nprotons03 += (pot_raw_id==03);
        nprotons103 += (pot_raw_id==103);
	npdiff0 += (pot_raw_id==03);
	npdiff0 -= (pot_raw_id==23);
	npdiff1 += (pot_raw_id==103);
	npdiff1 -= (pot_raw_id==123);

        float xi= (isLowPUrun ? 0.               : ev.fwdtrk_xi[ift]);
        float x=  (isLowPUrun ? ev.ppstrk_x[ift] :  0. );
	std::vector<TString> tags={"inc",ev.nvtx<v_low_cutoff ? "v_low" : "not_v_low"};

	//fill proton histograms
        ht.fill("csi",     xi,                    evWgt,tags);
        ht.fill("x",       x,                     evWgt,tags);
        ht.fill("nprotons",nprotons23+nprotons123,evWgt,tags);
        ht.fill("nprotons",nprotons23,            evWgt,tags);
	ht.fill("nprotons",nprotons123,           evWgt,tags);
        ht.fill("nprotons_new",nprotons03+nprotons103,evWgt,tags);
        ht.fill("nprotons_new",nprotons03,            evWgt,tags);
	ht.fill("nprotons_new",nprotons103,           evWgt,tags);
	ht.fill("npdiff0", npdiff0,               evWgt,tags);
	ht.fill("npdiff1", npdiff1,               evWgt,tags);
      }
 
      //proton energy loss
      float nrp23(0);
      float nrp123(0);
      //loop over forward trackers to store the number of protons in the
      //two RP:s we're considering (TODO: do this in an earlier loop)
      for (int ift=0; ift<ntrks; ift++) {
        //only near (pixels) detectors
	  const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123) continue;

	nrp23 += (pot_raw_id==23);
	nrp123 += (pot_raw_id==123);

      }
      
      //select events where one proton is captured by each RP
      if (nrp23!=1) continue;
      if (nrp123!=1) continue;
      
      //store xi for each RP
      float xi_23(0);
      float xi_123(0);
      for (int ift=0; ift<ntrks; ift++) {

        //only near (pixels) detectors
        const unsigned short pot_raw_id = (isLowPUrun ? ev.ppstrk_pot[ift] : ev.fwdtrk_pot[ift]);
        if (pot_raw_id!=23 && pot_raw_id!=123) continue;
	if (pot_raw_id==23) {
	  xi_23 = (isLowPUrun ? 0. : ev.fwdtrk_xi[ift]);
	  }
	if (pot_raw_id==123) {
	xi_123 = (isLowPUrun ? 0. : ev.fwdtrk_xi[ift]);
        }
	
      }
      //calculate proton energy according to P.Meiring's Eq. (9)
      //assuming 13 TeV collisions. Unit: TeV
      float lost_proton_energy = sqrt(13*xi_23*xi_123);
      lost_proton_energy_vect.push_back(lost_proton_energy);
      //fill 1D difference hists. Convert mlnjets to TeV
      ht.fill("Ecentral_minus_Eprotons", mlnjets/1000 - lost_proton_energy, 1);
      //fill slices
      if (0.0 <= mlnjets/1000 && mlnjets/1000 < 0.2){
	Ecentral_minus_Eprotons_00_02_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	signal_minus_bg_00_02_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	nevents_00_02_TeV++;
      }
      else if (0.2 <= mlnjets/1000 && mlnjets/1000 < 0.4) {
	Ecentral_minus_Eprotons_02_04_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	signal_minus_bg_02_04_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	nevents_02_04_TeV++;
      }
      else if (0.4 <= mlnjets/1000 && mlnjets/1000 < 0.6) {
        Ecentral_minus_Eprotons_04_06_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	signal_minus_bg_04_06_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	nevents_04_06_TeV++;
      }
      else if (0.6 <= mlnjets/1000 && mlnjets/1000 < 0.8) {
	Ecentral_minus_Eprotons_06_08_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	signal_minus_bg_06_08_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	nevents_06_08_TeV++;
      }
      else if (0.8 <= mlnjets/1000 && mlnjets/1000 < 1.0) {
	Ecentral_minus_Eprotons_08_10_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	signal_minus_bg_08_10_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	nevents_08_10_TeV++;
      }
      else if (1.0 <= mlnjets/1000 && mlnjets/1000 < 1.2) {
	Ecentral_minus_Eprotons_10_12_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	signal_minus_bg_10_12_TeV->Fill(mlnjets/1000 - lost_proton_energy, 1);
	nevents_10_12_TeV++;
      }

      signal_minus_bg->Fill(mlnjets/1000 - lost_proton_energy);
      ht.fill("Ecentral_minus_Eprotons_no_neutrino", mlnjets_no_neutrino/1000 - lost_proton_energy, 1, "");
      //fill 2D hist. Convert mlnjets to TeV
      protons_vs_CM_energy->Fill(mlnjets/1000, lost_proton_energy);
    }

  //Randomize mlnjets_vect, match random CM energies to lost proton energies to
  //simulate background
  //std::random_shuffle(mlnjets_vect.begin(), mlnjets_vect.end());
  std::random_shuffle(rand_proton_energy_vect.begin(), rand_proton_energy_vect.end());
  for (size_t i=0; i<mlnjets_vect.size(); i++) {
    //fill 1D difference hist. Convert mlnjets to TeV
    Ecentral_minus_Eprotons_bg->Fill(((mlnjets_vect[i]/1000) - rand_proton_energy_vect[i]), 1);
    //fill 2D hist. Convert mlnjets to TeV
    protons_vs_CM_energy_bg->Fill(mlnjets_vect[i]/1000, rand_proton_energy_vect[i]);
  }
  
  float lost_size(lost_proton_energy_vect.size());
  float mlnjets_size(mlnjets_vect.size());
  float scale_1(lost_size/mlnjets_size);
  Ecentral_minus_Eprotons_bg->Scale(scale_1);

  //fill bg slices
  //for (size_t i=0; i<nevents_00_02_TeV; i++) {
  for (size_t i=0; i<mlnjets_size; i++) {
    if (0.0 <= mlnjets_vect[i]/1000 && mlnjets_vect[i]/1000 < 0.2){
      Ecentral_minus_Eprotons_bg_00_02_TeV->Fill(mlnjets_vect[i]/1000 - rand_proton_energy_vect[i], 1);
      nevents_00_02_bg_TeV++;
    }
    //}
  //for (size_t i=0; i<nevents_02_04_TeV; i++) {
    if (0.2 <= mlnjets_vect[i+nevents_00_02_TeV]/1000 && mlnjets_vect[i+nevents_00_02_TeV]/1000 < 0.4){
      Ecentral_minus_Eprotons_bg_02_04_TeV->Fill(mlnjets_vect[i+nevents_00_02_TeV]/1000 - rand_proton_energy_vect[i+nevents_00_02_TeV], 1);
      nevents_02_04_bg_TeV++;
    }
    //}
    //for (size_t i=0; i<nevents_04_06_TeV; i++) {
    if (0.4 <= mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV]/1000 && mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV]/1000 < 0.6){
      Ecentral_minus_Eprotons_bg_04_06_TeV->Fill(mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV]/1000 - rand_proton_energy_vect[i+nevents_00_02_TeV+nevents_02_04_TeV], 1);
      nevents_04_06_bg_TeV++;
    }
    //}
    //for (size_t i=0; i<nevents_06_08_TeV; i++) {
    if (0.6 <= mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV]/1000 && mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV]/1000 < 0.8){
      Ecentral_minus_Eprotons_bg_06_08_TeV->Fill(mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV]/1000 - rand_proton_energy_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV], 1);
      nevents_06_08_bg_TeV++;
    }
    //}
    //for (size_t i=0; i<nevents_08_10_TeV; i++) {
    if (0.8<= mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV]/1000 && mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV]/1000 < 1.0){
      Ecentral_minus_Eprotons_bg_08_10_TeV->Fill(mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV]/1000 - rand_proton_energy_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV], 1);
      nevents_08_10_bg_TeV++;
    }
    //}
    //for (size_t i=0; i<nevents_10_12_TeV; i++) {
    if (1.0 <= mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV+nevents_08_10_TeV]/1000 && mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV+nevents_08_10_TeV]/1000 < 1.2){
      Ecentral_minus_Eprotons_bg_10_12_TeV->Fill(mlnjets_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV+nevents_08_10_TeV]/1000 - rand_proton_energy_vect[i+nevents_00_02_TeV+nevents_02_04_TeV+nevents_04_06_TeV+nevents_06_08_TeV+nevents_08_10_TeV], 1);
      nevents_10_12_bg_TeV++;
    }
  }
  //Scale bg sliced histo's
  
  if(nevents_00_02_bg_TeV!=0 && nevents_00_02_TeV!=0){
    float scale_00_02(nevents_00_02_TeV/nevents_00_02_bg_TeV);
    Ecentral_minus_Eprotons_bg_00_02_TeV->Scale(scale_00_02);
  }
  if(nevents_02_04_bg_TeV!=0 && nevents_02_04_TeV!=0){
    float scale_02_04(nevents_02_04_TeV/nevents_02_04_bg_TeV);
    Ecentral_minus_Eprotons_bg_02_04_TeV->Scale(scale_02_04);
  }
  if(nevents_04_06_bg_TeV!=0 && nevents_04_06_TeV!=0){
    float scale_04_06(nevents_04_06_TeV/nevents_04_06_bg_TeV);
    Ecentral_minus_Eprotons_bg_04_06_TeV->Scale(scale_04_06);
  }
  if(nevents_06_08_bg_TeV!=0 && nevents_06_08_TeV!=0){
    float scale_06_08(nevents_06_08_TeV/nevents_06_08_bg_TeV);
    Ecentral_minus_Eprotons_bg_06_08_TeV->Scale(scale_06_08);
  }
  if(nevents_08_10_bg_TeV!=0 && nevents_08_10_TeV!=0){
    float scale_08_10(nevents_08_10_TeV/nevents_08_10_bg_TeV);
    Ecentral_minus_Eprotons_bg_08_10_TeV->Scale(scale_08_10);
  }
  if(nevents_10_12_bg_TeV!=0 && nevents_10_12_TeV!=0){
    float scale_10_12(nevents_10_12_TeV/nevents_10_12_bg_TeV);
    Ecentral_minus_Eprotons_bg_10_12_TeV->Scale(scale_10_12);
  }

  //cout << "float_06_08 = " << nevents_06_08_TeV << endl;
  //cout << "float_06_08 = " << nevents_06_08 << endl;
  //cout << "float_06_08_bg = " << nevents_06_08_bg_TeV << endl;
  //cout << "float_06_08_bg = " << nevents_06_08_bg << endl;
  //cout << "scale_06_08 = " << scale_06_08 << endl;

  //Substract normalised backgrounds from selection
  signal_minus_bg_00_02_TeV->Add(Ecentral_minus_Eprotons_bg_00_02_TeV, -1);
  signal_minus_bg_02_04_TeV->Add(Ecentral_minus_Eprotons_bg_02_04_TeV, -1);
  signal_minus_bg_04_06_TeV->Add(Ecentral_minus_Eprotons_bg_04_06_TeV, -1);
  signal_minus_bg_06_08_TeV->Add(Ecentral_minus_Eprotons_bg_06_08_TeV, -1);
  signal_minus_bg_08_10_TeV->Add(Ecentral_minus_Eprotons_bg_08_10_TeV, -1);
  signal_minus_bg_10_12_TeV->Add(Ecentral_minus_Eprotons_bg_10_12_TeV, -1);
  signal_minus_bg->Add(Ecentral_minus_Eprotons_bg, -1);


  //Add hists to histtool
  ht.addHist("signal_minus_bg", signal_minus_bg);
  ht.addHist("Ecentral_minus_Eprotons_bg", Ecentral_minus_Eprotons_bg);
  ht.addHist("Eprotons_vs_Ecentral", protons_vs_CM_energy);
  ht.addHist("Eprotons_vs_Ecentral_bg", protons_vs_CM_energy_bg);
  ht.addHist("Ecentral_minus_Eprotons_bg_00_02_TeV", Ecentral_minus_Eprotons_bg_00_02_TeV);
  ht.addHist("Ecentral_minus_Eprotons_bg_02_04_TeV", Ecentral_minus_Eprotons_bg_02_04_TeV);
  ht.addHist("Ecentral_minus_Eprotons_bg_04_06_TeV", Ecentral_minus_Eprotons_bg_04_06_TeV);
  ht.addHist("Ecentral_minus_Eprotons_bg_06_08_TeV", Ecentral_minus_Eprotons_bg_06_08_TeV);
  ht.addHist("Ecentral_minus_Eprotons_bg_08_10_TeV", Ecentral_minus_Eprotons_bg_08_10_TeV);
  ht.addHist("Ecentral_minus_Eprotons_00_02_TeV", Ecentral_minus_Eprotons_00_02_TeV);
  ht.addHist("Ecentral_minus_Eprotons_02_04_TeV", Ecentral_minus_Eprotons_02_04_TeV);
  ht.addHist("Ecentral_minus_Eprotons_04_06_TeV", Ecentral_minus_Eprotons_04_06_TeV);
  ht.addHist("Ecentral_minus_Eprotons_06_08_TeV", Ecentral_minus_Eprotons_06_08_TeV);
  ht.addHist("Ecentral_minus_Eprotons_08_10_TeV", Ecentral_minus_Eprotons_08_10_TeV);
  ht.addHist("signal_minus_bg_00_02_TeV", signal_minus_bg_00_02_TeV);
  ht.addHist("signal_minus_bg_02_04_TeV", signal_minus_bg_02_04_TeV);
  ht.addHist("signal_minus_bg_04_06_TeV", signal_minus_bg_04_06_TeV);
  ht.addHist("signal_minus_bg_06_08_TeV", signal_minus_bg_06_08_TeV);
  ht.addHist("signal_minus_bg_08_10_TeV", signal_minus_bg_08_10_TeV);
  ht.addHist("signal_minus_bg_10_12_TeV", signal_minus_bg_10_12_TeV);

  cout << "lost = " <<  lost_proton_energy_vect.size() << endl;
  cout << "rand = " << rand_proton_energy_vect.size() << endl;
  cout << "mlnjets = " << mlnjets_vect.size() << endl;

  //Write 2D hists to file
  auto output_1 = new TCanvas("Eprotons_vs_Ecentral.root");
  protons_vs_CM_energy->Draw();
  linear_line->Draw("Same");  
  TFile out_protons_vs_CM_energy("Eprotons_vs_Ecentral.root","RECREATE");
  output_1->Write();
  out_protons_vs_CM_energy.Close();

  auto output_2 = new TCanvas("Eprotons_vs_Ecentral_bg.root");
  protons_vs_CM_energy_bg->Draw();
  linear_line->Draw("Same");
  TFile out_protons_vs_CM_energy_bg("Eprotons_vs_Ecentral_bg.root","RECREATE");
  output_2->Write();
  out_protons_vs_CM_energy_bg.Close();

  //close input file
  f->Close();
  
  //save histos to file  
  fOut->cd();
  for (auto& it : ht.getPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    if(it.second->GetEntries()==0) continue;
    it.second->SetDirectory(fOut); it.second->Write(); 
  }  
  fOut->Close();
}

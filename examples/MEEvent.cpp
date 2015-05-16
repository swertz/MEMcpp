#include "TLorentzVector.h"

#include "MEEvent.h"

MEEvent::MEEvent(){
  //hst_TTbar = new TH1D("test_TTbar", "test_TTbar", BINNING, START, STOP); 
}

void MEEvent::SetVectors(const TLorentzVector ep, const TLorentzVector mum, const TLorentzVector b, const TLorentzVector bbar, const TLorentzVector met){
  p3 = ep;
  p5 = mum;
  p4 = b;
  p6 = bbar;
  Met = met;

  /*if(hst_TTbar)
    hst_TTbar->Reset();*/
}

MEEvent::~MEEvent(){
  //delete hst_TTbar; hst_TTbar = NULL;
}

/*TH1D* MEEvent::GetTTbar(){
  return hst_TTbar;
}*/

/*void MEEvent::writeHists(void){
  TCanvas *c = new TCanvas(TString("Mttbar")+"_"+SSTR(count_wgt)+"_"+SSTR(count_perm),"Canvas for plotting",600,600);
  c->cd();
  hst_TTbar->Draw();
  c->Write();
  c->Print(TString("plots/Mttbar_")+SSTR(count_wgt)+"_"+SSTR(count_perm)+".png");
  delete c; c = 0;
}*/

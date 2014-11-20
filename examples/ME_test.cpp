/*
root -l examples/THDM_width.C\(\"delphes_output.root\"\)
*/

//------------------------------------------------------------------------------


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TLorentzVector.h"
#include "TMath.h"

#include <sstream>
#include <iomanip>
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TFile.h"
#include "TObject.h"



#include "ExRootTreeReader.h"
#include "ExRootTreeWriter.h"
#include "ExRootTreeBranch.h"
#include "ExRootResult.h"

#include "DelphesClasses.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TSystem.h"

#include <iomanip>

#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/thdm_cpp/src/HelAmps_2HDM4MG5.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/thdm_cpp/SubProcesses/P0_Sigma_2HDM4MG5_gg_bbxepem/CPPProcess.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/thdm_cpp/src/rambo.h"
#include "/home/fynu/amertens/scratch/MadGraph/madgraph5/thdm_cpp/src/Parameters_2HDM4MG5.h"


void InRestFrame(TLorentzVector v1, TLorentzVector v2){
  Float_t gamma = v1.E()/v1.M();
  Float_t betax = v1.Px()/(gamma*v1.M());
  Float_t betay = v1.Py()/(gamma*v1.M());
  Float_t betaz = v1.Pz()/(gamma*v1.M());

  v2.Boost(-betax,-betay,-betaz);
}


Float_t Theta_Star(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3){

  Float_t gamma = v3.E()/v3.M();
  Float_t betax = -v3.Px()/(gamma*v3.M());
  Float_t betay = -v3.Py()/(gamma*v3.M());
  Float_t betaz = -v3.Pz()/(gamma*v3.M());

  v1.Boost(betax,betay,betaz);
  v2.Boost(betax,betay,betaz);

  Float_t Theta =  v1.Angle(v2.Vect());

  return Theta;
}


Float_t Delta_Phi(TLorentzVector b1, TLorentzVector b2, TLorentzVector l1, TLorentzVector l2, TLorentzVector Zs){
  Float_t gamma = Zs.E()/Zs.M();
  Float_t betax = -Zs.Px()/(gamma*Zs.M());
  Float_t betay = -Zs.Py()/(gamma*Zs.M());
  Float_t betaz = -Zs.Pz()/(gamma*Zs.M());

  b1.Boost(betax,betay,betaz);
  b2.Boost(betax,betay,betaz);
  l1.Boost(betax,betay,betaz);
  l2.Boost(betax,betay,betaz);

  TVector3 plan1 = b1.Vect().Cross(b2.Vect());
  TVector3 plan2 = l1.Vect().Cross(l2.Vect());

  Float_t CosDelta = ( (plan1 * plan2)  / TMath::Sqrt((plan1*plan1)*(plan2*plan2)));
  // Float_t SinDelta = ( (plan1.Cross(plan2)).Mag() / (plan1.Mag()*plan2.Mag()));
  Float_t Phi12 = plan1.DeltaPhi(plan2);

  Float_t Delta;

  if (Phi12 > 0) Delta = acos (CosDelta);
  else Delta = - acos (CosDelta);

  //cout << "phi12 : " << Phi12 << " cos : " << CosDelta << " Delta : " << Delta << endl;
  return Delta;
}

int main(int argc, char *argv[])
//const char *inputFile,const char *outputFile
{

  TString inputFile(argv[1]);
  TString outputFile(argv[2]);

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElec = treeReader->UseBranch("Electron");
  TClonesArray *branchGen = treeReader->UseBranch("Particle");

  // Book histograms
  TH1 *histllMass = new TH1F("ll_Mass", "M_{inv}(l_{1},l_{2})", 120, 0.0, 120.0);
  TH1 *histbbMass = new TH1F("bb_Mass", "M_{inv}(b_{1},b_{2})", 400, 0.0, 1200.0);
  TH1 *histllbbMass = new TH1F("llbb_Mass", "M_{inv}(b_{1}, b_{2}, #mu_{1}, #mu_{2})", 400, 100.0, 1700.0);
  TH1 *histME = new TH1D("ME4THDM", "ME4THDM", 60, 0.0, 30);


  TH1 *hist_cosS_Z = new TH1D("cosS_Z","cosS_Z",8,-1,1);
  TH1 *hist_cosS_A = new TH1D("cosS_A","cosS_A",8,-1,1);
  TH1 *hist_cosS_H = new TH1D("cosS_H","cosS_H",8,-1,1);
  TH1 *hist_Thetal_gen = new TH1D("Theta_l","Theta_l",20,-1,1);
  TH1 *hist_Thetab_gen = new TH1D("Theta_b","Theta_b",20,-1,1);
  TH1 *hist_ThetaS_gen = new TH1D("Theta_S","Theta_S",20,-1,1);
  TH1 *hist_DeltaPhi_gen = new TH1D("Delta_Phi","Delta_Phi",20,-TMath::Pi(),TMath::Pi());

  TH1 *hist_ThetaS_Zgen = new TH1D("ThetaS_Zgen","ThetaS_Zgen",20,0,TMath::Pi());

  TH2 *histME_bbMass = new TH2D("ME4THDM_vs_bbMass", "ME4THDM_vs_bbMass", 1000, -30, 30 ,400, 0.0, 1200.0);
  TH2 *histME_llbbMass = new TH2D("ME4THDM_vs_llbbMass", "ME4THDM_vs_llbbMass", 1000, -30, 30,400, 100.0, 1700.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
  
    // If event contains at least 1 jet

    vector<TLorentzVector> bjets;
    vector<TLorentzVector> muons;
    vector<TLorentzVector> elecs;

    TLorentzVector gen_lm, gen_lp, gen_b1, gen_b2;

    GenParticle *gen;

    for (Int_t i = 0; i < branchGen->GetEntries(); i++)
      {
      gen = (GenParticle *) branchGen->At(i);
      if (gen->Status == 3)
        {
        if (gen->PID == -11 || gen->PID == -13) gen_lm = gen->P4();
        else if (gen->PID == 11 || gen->PID == 13) gen_lp = gen->P4();
        else if (gen->PID == -5) gen_b1 = gen->P4();
        else if (gen->PID == 5) gen_b2 = gen->P4();
        }
      }

    TLorentzVector beam;
    beam.SetPxPyPzE(0,0,4000,4000);

    Double_t Thetab, Thetal, ThetaS, DeltaPhi;

    Thetab = cos(Theta_Star(gen_b1,gen_lm+gen_lp,gen_b1+gen_b2));
    Thetal = cos(Theta_Star(gen_lm,gen_b1+gen_b2,gen_lm+gen_lp));
    ThetaS = cos(Theta_Star(gen_b1+gen_b2,beam,gen_lm+gen_lp+gen_b1+gen_b2));
    DeltaPhi = Delta_Phi(gen_lm+gen_lp, beam, gen_lp, gen_lm, gen_lm+gen_lp+gen_b1+gen_b2);

    hist_Thetal_gen->Fill(Thetal);
    hist_Thetab_gen->Fill(Thetab);
    hist_ThetaS_gen->Fill(ThetaS);
    hist_DeltaPhi_gen->Fill(DeltaPhi);

    TLorentzVector jet_4v;
    TLorentzVector muon_4v;
    TLorentzVector elec_4v;
    Jet* jet;
    Muon* muon;
    Electron* elec;

    for ( Int_t i = 0; i < branchJet->GetEntries(); i++){
      jet = (Jet *) branchJet->At(i);
      if (jet->PT > 30 && jet->BTag == 1 && abs(jet->Eta) < 2.4){
        jet_4v.SetPtEtaPhiM(jet->PT, jet->Eta, jet->Phi, jet->Mass); 
        bjets.push_back(jet_4v);
	}
      }

    for ( Int_t i = 0; i < branchMuon->GetEntries(); i++){
      muon = (Muon *) branchMuon->At(i);
      if (muon->PT > 20 && abs(muon->Eta) < 2.4){
        muon_4v.SetPtEtaPhiM(muon->PT, muon->Eta, muon->Phi, 0.0);   
        muons.push_back(muon_4v);
        }
      }

    for ( Int_t i = 0; i < branchElec->GetEntries(); i++){
      elec = (Electron *) branchElec->At(i);
      if (elec->PT > 20 && abs(elec->Eta) < 2.5){
        elec_4v.SetPtEtaPhiM(elec->PT, elec->Eta, elec->Phi, 0.0);
        elecs.push_back(elec_4v);
        }
      }

    if (muons.size() > 1 || elecs.size() > 1){

      TLorentzVector ll;
      if (muons.size() > 1) { ll = muons[0] + muons[1]; histllMass->Fill(ll.M());}
      else if (elecs.size() > 1)  { ll = elecs[0] + elecs[1]; histllMass->Fill(ll.M());}

  if (bjets.size() > 1 && ll.M() > 76 && ll.M() < 106 && muons.size() > 1){
        TLorentzVector bb = bjets[0] + bjets[1];
        TLorentzVector llbb = bjets[0] + bjets[1] + ll;
        histbbMass->Fill(bb.M());
        histllbbMass->Fill(llbb.M());

	hist_cosS_Z->Fill(cos(Theta_Star(muons[0],bb,ll)));
	hist_cosS_A->Fill(cos(Theta_Star(bjets[0],ll,bb)));
	hist_cosS_H->Fill(cos(Theta_Star(bb,beam,llbb)));

	cout << "create gluons initial states : " ;
        //TLorentzVector *g1;
	//TLorentzVector *g2;
	Double_t totEcm = llbb.M();
	Double_t gEcm = totEcm/2.0;
	
	cout << gEcm << " GeV" << endl;

        //g1->SetPxPyPzE(0.0,0.0,150,150);
        //g2->SetPxPyPzE(0.0,0.0,-150,150);

	cout << "changing the frame" << endl;

	InRestFrame(llbb, muons[0]);
	InRestFrame(llbb, muons[1]);
	InRestFrame(llbb, bjets[0]);
	InRestFrame(llbb, bjets[1]);


	//InRestFrame(llbb, g2);

	// Create a process object
        CPPProcess process;

        // Read param_card and set parameters
        process.initProc("/home/fynu/amertens/scratch/MadGraph/madgraph5/thdm_cpp/Cards/param_card.dat");


	cout << "Px : " << bjets[0].Px() << " " << bjets[0].Py() << " " << bjets[0].Pz() <<" " << bjets[0].E() <<" " << endl;	
	//cout << "g2 : " << g2->Px() << " " << g2->Py() << " " << g2->Pz() <<" " << g2->E() <<" " << endl;
	// momentum vector definition
	vector<double*> p(1, new double[4]);
	p[0][0] = gEcm; p[0][1] = 0.0; p[0][2] = 0.0; p[0][3] = gEcm;
	p.push_back(new double[4]);
	p[1][0] = gEcm; p[1][1] = 0.0; p[1][2] = 0.0; p[1][3] = -gEcm;
	p.push_back(new double[4]);
	p[2][0] = bjets[0].E(); p[2][1] = bjets[0].Px(); p[2][2] = bjets[0].Py(); p[2][3] = bjets[0].Pz();
	p.push_back(new double[4]);
	p[3][0] = bjets[1].E(); p[3][1] = bjets[1].Px(); p[3][2] = bjets[1].Py(); p[3][3] = bjets[1].Pz();
	p.push_back(new double[4]);
	p[4][0] = muons[0].E(); p[4][1] = muons[0].Px(); p[4][2] = muons[0].Py(); p[4][3] = muons[0].Pz();
	p.push_back(new double[4]);
	p[5][0] = muons[1].E(); p[5][1] = muons[1].Px(); p[5][2] = muons[1].Py(); p[5][3] = muons[1].Pz();

	  // Set momenta for this event
	process.setMomenta(p);

	process.setmh2mh3(300,500);


	// Evaluate matrix element
	process.sigmaKin();

	std::cout<<"ME computation : " ;
	const double* matrix_elements = process.getMatrixElements();

	cout << matrix_elements[0] << endl;
	
	histME->Fill(TMath::Log10(matrix_elements[0]));
	histME_bbMass->Fill(TMath::Log10(matrix_elements[0]), bb.M() );
	histME_llbbMass->Fill(TMath::Log10(matrix_elements[0]), llbb.M());
        }
      }
    }
  // Show resulting histograms
  TCanvas *c1 = new TCanvas("c1","multipads",900,700);
  c1->Divide(4);
  c1->cd(1);
  histllMass->Draw();
  c1->cd(2);
  histbbMass->Draw();
  c1->cd(3);
  histllbbMass->Draw();
  c1->cd(4);
  histME->Draw();

  TFile f(outputFile,"recreate");
  histllMass->Write();
  histbbMass->Write();
  histllbbMass->Write();
  histME->Write();
  hist_cosS_Z->Write();
  hist_cosS_A->Write();
  hist_cosS_H->Write();
  hist_Thetal_gen->Scale(1.0/hist_Thetal_gen->GetEntries());
  hist_Thetal_gen->Write();
  hist_Thetab_gen->Scale(1.0/hist_Thetab_gen->GetEntries());
  hist_Thetab_gen->Write();
  hist_ThetaS_gen->Scale(1.0/hist_ThetaS_gen->GetEntries());
  hist_ThetaS_gen->Write();
  hist_DeltaPhi_gen->Scale(1.0/hist_DeltaPhi_gen->GetEntries());
  hist_DeltaPhi_gen->Write();

  histME_bbMass->Write();
  histME_llbbMass->Write();
  hist_ThetaS_Zgen->Scale(1.0/hist_ThetaS_Zgen->GetEntries());
  hist_ThetaS_Zgen->Write();
  f.Write();

}

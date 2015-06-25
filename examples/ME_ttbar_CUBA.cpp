#include <string>
#include <iostream>

#include "classes/DelphesClasses.h"

#include "TStopwatch.h"
#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "MEWeight.h"
#include "utils.h"

#define BINNING 1750
#define START 250
#define STOP 2000

using namespace std;

int main(int argc, char *argv[])
{
  std::string inputFile(argv[1]);
  std::string outputFile(argv[2]);
  std::string fileTF(argv[3]);
  int start_evt = atoi(argv[4]);
  int end_evt = atoi(argv[5]);

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile.c_str());

  double MadWeight, MadWeight_Error;
  chain.SetBranchAddress("Weight_TT", &MadWeight);
  chain.SetBranchAddress("Weight_TT_Error", &MadWeight_Error);

  // Get pointers to branches used in this analysis
  TClonesArray *branchGen = NULL;
  chain.SetBranchAddress("Particle", &branchGen);

  cout << "Entries:" << chain.GetEntries() << endl;

  TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
  TTree* outTree = chain.CloneTree(0);

  double Weight_TT_cpp, Weight_TT_Error_cpp;
  bool Weighted_TT_cpp;
  double time;
  outTree->Branch("Weight_TT_cpp", &Weight_TT_cpp);
  outTree->Branch("Weight_TT_Error_cpp", &Weight_TT_Error_cpp);
  outTree->Branch("Weighted_TT_cpp", &Weighted_TT_cpp);
  outTree->Branch("Weight_TT_cpp_time", &time);

  double mTruth, mAverage, mMaxL;
  outTree->Branch("MTTbar_MCTruth", &mTruth);
  outTree->Branch("MTTbar_DMEMAverage", &mAverage);
  outTree->Branch("MTTbar_DMEMMaxL", &mMaxL);

  if(end_evt >= chain.GetEntries())
    end_evt = chain.GetEntries()-1;
  
  MEWeight* myWeight = new MEWeight("/home/fynu/swertz/scratch/Madgraph/madgraph5/cpp_ttbar_epmum_NOPOW/Cards/param_card.dat", "cteq6l1", fileTF);
  myWeight->AddTF("electron", "Binned_Egen_DeltaE_Norm_ele");
  myWeight->AddTF("muon", "Binned_Egen_DeltaE_Norm_muon");
  myWeight->AddTF("jet", "Binned_Egen_DeltaE_Norm_jet");

  TH1D* truth_TTbar = new TH1D("MTruth_TTbar", "M_{tt}  Truth", 1750, 250, 2000);

  for(int entry = start_evt; entry <= end_evt ; ++entry){
    // Load selected branches with data from specified event
    chain.GetEntry(entry);

    TLorentzVector gen_ep, gen_mum, gen_b, gen_bbar, gen_Met;

    GenParticle *gen;

    //int count_ep=0, count_mum=0;

    for (int i = 0; i < branchGen->GetEntries(); i++){
      gen = (GenParticle*) branchGen->At(i);
      //cout << "Status=" << gen->Status << ", PID=" << gen->PID << ", E=" << gen->P4().E() << endl;
      if (gen->Status == 1){
        if (gen->PID == -11){
          gen_ep = gen->P4();
          //count_ep++;
        }else if (gen->PID == 13){
          gen_mum = gen->P4();
          //count_mum++;
        }
        else if (gen->PID == 12) gen_Met += gen->P4();
        else if (gen->PID == -14) gen_Met += gen->P4();
        else if (gen->PID == 5) gen_b = gen->P4();
        else if (gen->PID == -5) gen_bbar = gen->P4();
      }
    }

    //if(count_ep != 1 || count_mum != 1)
    //  continue;
    //gen_Met.SetPz(0.);
  
    cout << "From MadGraph:" << endl;
    cout << "Electron" << endl;
    cout << gen_ep.E() << "," << gen_ep.Px() << "," << gen_ep.Py() << "," << gen_ep.Pz() << endl;
    cout << "b quark" << endl;
    cout << gen_b.E() << "," << gen_b.Px() << "," << gen_b.Py() << "," << gen_b.Pz() << endl;
    cout << "Muon" << endl;
    cout << gen_mum.E() << "," << gen_mum.Px() << "," << gen_mum.Py() << "," << gen_mum.Pz() << endl;
    cout << "Anti b quark" << endl;
    cout << gen_bbar.E() << "," << gen_bbar.Px() << "," << gen_bbar.Py() << "," << gen_bbar.Pz() << endl;
    cout << "MET" << endl;
    cout << gen_Met.E() << "," << gen_Met.Px() << "," << gen_Met.Py() << "," << gen_Met.Pz() << endl << endl;

    Weight_TT_cpp = 0.;
    Weight_TT_Error_cpp = 0.;
    time = 0.;
  
    TStopwatch chrono;
    chrono.Start();
    
    mTruth = (gen_ep + gen_b + gen_mum + gen_bbar + gen_Met).M();
    truth_TTbar->Fill(mTruth);
    cout << "M_ttbar truth: " << mTruth << endl;

    for(int permutation = 1; permutation <= 2; permutation++){

      //permutation = 1;

      double weight = 0;
      double error = 0;

      if(permutation == 1)
        myWeight->SetEvent(gen_ep, gen_mum, gen_b, gen_bbar, gen_Met);
      if(permutation == 2)
        myWeight->SetEvent(gen_ep, gen_mum, gen_bbar, gen_b, gen_Met);

      weight = myWeight->ComputeWeight(error)/2.; 

      Weight_TT_cpp += weight;
      Weight_TT_Error_cpp += SQ(error/2.);

      //break;
    }

    time = chrono.CpuTime();
    
    Weight_TT_Error_cpp = sqrt(Weight_TT_Error_cpp);
    Weighted_TT_cpp = true;

    cout << "====> Event " << entry << ": weight = " << Weight_TT_cpp << " +- " << Weight_TT_Error_cpp << endl;
    cout << "      CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << endl;
    cout << "      MadWeight: " << MadWeight << " +- " << MadWeight_Error << endl << endl;

    mAverage = myWeight->GetTempAverage();
    mMaxL = myWeight->GetTempMaxLikelihood();
    cout << "      Most likely M_ttbar : " << mMaxL << endl;
    cout << "      Average M_ttbar     : " << mAverage << endl << endl;
    myWeight->AddAndResetTempHist();

    outTree->Fill();

    //count_wgt++;
  }
  
  outFile->cd();

  truth_TTbar->Scale(1./truth_TTbar->Integral());
  truth_TTbar->Write();
  myWeight->WriteHist();

  outTree->Write();
  delete truth_TTbar; truth_TTbar = nullptr;
  delete myWeight; myWeight = nullptr;
  delete outFile; outFile = nullptr;
}


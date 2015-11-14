#include <string>
#include <iostream>
#include <iomanip>
#include <memory>
#include <map>

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/PtEtaPhiE4D.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>

#include "SubProcesses/P1_Sigma_TopEffTh_gg_epvebemvexbx/cpp_TopEffTh_pp_ttx_epem_NP2_QED1.h"
typedef cpp_TopEffTh_pp_ttx_epem_NP2_QED1 MyTopEffThProcess;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> LorentzVectorM;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double>> LorentzVectorE;

#include "MEWeight.h"
#include "TopEffTh/TopEffTh.h"

extern const std::vector<std::string> hypothesisNames { "SM", "OtG", "OG", "OphiG", "O81qq", "O83qq", "O8ut", "O8dt" };

using namespace std;

int main(int argc, char *argv[])
{
  std::string inputFile(argv[1]);
  std::string outputFile(argv[2]);
  std::string fileTF(argv[3]);
  int start_evt = atoi(argv[4]);
  int end_evt = atoi(argv[5]);

  cout << "Reading from file  " << inputFile.c_str() << endl;
  cout << "Will write to file " << outputFile.c_str() << endl;

  // Create chain of root trees
  unique_ptr<TChain> chain( new TChain("t") );
  chain->Add(inputFile.c_str());

  LorentzVectorM *lep1(nullptr), *lep2(nullptr), *bjet1(nullptr), *bjet2(nullptr);
  float MET_met, MET_phi;
  bool IsEE, IsEMu, IsMuMu;
  int leadLepPID;

  chain->SetBranchAddress("lep1_p4", &lep1);
  chain->SetBranchAddress("lep2_p4", &lep2);
  chain->SetBranchAddress("bjet1_p4", &bjet1);
  chain->SetBranchAddress("bjet2_p4", &bjet2);
  chain->SetBranchAddress("MET_met", &MET_met);
  chain->SetBranchAddress("MET_phi", &MET_phi);
  chain->SetBranchAddress("IsEE", &IsEE);
  chain->SetBranchAddress("IsEMu", &IsEMu);
  chain->SetBranchAddress("IsMuMu", &IsMuMu);
  chain->SetBranchAddress("leadLepPID", &leadLepPID);

  cout << "Entries:" << chain->GetEntries() << endl;

  std::unique_ptr<TFile> outFile( new TFile(outputFile.c_str(), "RECREATE") );
  TTree *outTree( chain->CloneTree(0) );

  std::vector<double> weights(hypothesisNames.size(), 0.);
  std::vector<double> weightErrors(hypothesisNames.size(), 0.);
  for(size_t i = 0; i < hypothesisNames.size(); ++i){
    std::string baseName = "Weight_" + hypothesisNames[i];
    outTree->Branch(baseName.c_str(), &weights[i]);
    outTree->Branch((baseName + "_Err").c_str(), &weightErrors[i]);
  }
  double time;
  outTree->Branch("MEMpp_time_s", &time);

  if(end_evt >= chain->GetEntries())
    end_evt = chain->GetEntries()-1;

  // Create CPPProcess and MEWeight objects
  Parameters_TopEffTh myParams("/home/fynu/swertz/tests_MEM/MEMcpp/MatrixElements/cpp_TopEffTh_pp_ttx_epem_NP2_QED1/Cards/param_card.dat");
  MyTopEffThProcess myProcess(myParams);
  unique_ptr<MEWeight> myWeight( new MEWeight(myProcess, "cteq6l1", fileTF) );
  
  myWeight->AddTF("electron", "Binned_Egen_DeltaE_Norm_ele");
  myWeight->AddTF("muon", "Binned_Egen_DeltaE_Norm_muon");
  myWeight->AddTF("jet", "Binned_Egen_DeltaE_Norm_jet");

  for(int entry = start_evt; entry <= end_evt ; entry++){
    // Load selected branches with data from specified event
    chain->GetEntry(entry);

    LorentzVectorE lep1_e(lep1->Pt(), lep1->Eta(), lep1->Phi(), lep1->E());
    LorentzVectorE lep2_e(lep2->Pt(), lep2->Eta(), lep2->Phi(), lep2->E());
    LorentzVectorE bjet1_e(bjet1->Pt(), bjet1->Eta(), bjet1->Phi(), bjet1->E());
    LorentzVectorE bjet2_e(bjet2->Pt(), bjet2->Eta(), bjet2->Phi(), bjet2->E());
    LorentzVectorE met_e(MET_met, 0., MET_phi, MET_met);

    cout << "Input event " << entry << ":\n";
    cout << "  Lepton 1: " << lep1_e << endl;
    cout << "  Lepton 2: " << lep2_e << endl;
    cout << "  B-jet 1 : " << bjet1_e << endl;
    cout << "  B-jet 2 : " << bjet2_e << endl;
    cout << "  MET     : " << met_e << endl;

    time = 0.;
  
    TStopwatch chrono;
    chrono.Start();

    LepType myLepType(LepType::Unknown);
    if(IsEE){
      myLepType = LepType::EE;
    }else if(IsMuMu){
      myLepType = LepType::MuMu;
    }else if(IsEMu){
      if(std::abs(leadLepPID) == 11)
        myLepType = LepType::EMu;
      else
        myLepType = LepType::MuE;
    }else{
      cout << "ERROR: Something went wrong: unknown leptonic decay type.\n";
      continue;
    }

    myWeight->SetEvent(lep1_e, lep2_e, bjet1_e, bjet2_e, met_e, myLepType);
    myWeight->ComputeWeight(weights, weightErrors);

    time = chrono.CpuTime();
    
    cout << "\n====> Event " << entry << ":" << endl;
    for(size_t i = 0; i < hypothesisNames.size(); ++i)
      cout << "      Hypothesis " << hypothesisNames[i] << ": " << weights[i] << " +- " << weightErrors[i] << " (" << setprecision(1) << fixed << 100*std::abs(weightErrors[i]/weights[i]) << "\%)" << endl << scientific << setprecision(6);
    cout << "      CPU time : " << chrono.CpuTime() << "  Real-time : " << chrono.RealTime() << endl << endl;

    outTree->Fill();
  }
 
  outFile->cd();
  outTree->Write();

  return 0;
}


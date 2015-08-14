void add(TString target, TString source){
    gSystem->Load("libDelphes");

    TChain chain = TChain("Delphes");
    chain.Add(source);

    TFile outFile(target, "RECREATE");
    TTree *outTree = chain.CloneTree(0);

    for(unsigned long i =0; i < chain.GetEntries(); ++i){
        chain.GetEntry(i);
        outTree->Fill();
    }

    outTree->Write();
    outFile.Close();
}

void lhco_from_root(const char* in, const char* out){
	gSystem->Load("~/storage/Delphes/Delphes-3.1.2/libDelphes.so");

	TChain chain("Delphes");
	chain.Add(in);

    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    TClonesArray *branchGen = treeReader->UseBranch("Particle");

    ofstream fout(out);

    for(Int_t entry = 0; entry < chain.GetEntries(); ++entry){
	  //if(!entry%1000) 
	  cout << "Reading entry " << entry << endl;
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);

      TLorentzVector gen_ep, gen_mum, gen_b, gen_bbar, gen_Met;

      GenParticle *gen;

      for (Int_t i = 0; i < branchGen->GetEntries(); i++){
        gen = (GenParticle *) branchGen->At(i);
        if (gen->Status == 1){
          if (gen->PID == -11) gen_ep = gen->P4();
          else if (gen->PID == 13) gen_mum = gen->P4();
      	  else if (gen->PID == 12) gen_Met += gen->P4();
      	  else if (gen->PID == -14) gen_Met += gen->P4();
      	  else if (gen->PID == 5) gen_b = gen->P4();
      	  else if (gen->PID == -5) gen_bbar = gen->P4();
		}
	 }
	
	 fout << "0 " << entry+1 << " 0" << endl;
	 fout << "1 1 " << gen_ep.Eta() << " " << gen_ep.Phi()+TMath::Pi() << " " << gen_ep.Pt() << " 0.000511 1 0 0 0 0" << endl;
	 fout << "2 2 " << gen_mum.Eta() << " " << gen_mum.Phi()+TMath::Pi() << " " << gen_mum.Pt() << " 0.105658 -1 0 0 0 0" << endl;
	 fout << "3 4 " << gen_b.Eta() << " " << gen_b.Phi()+TMath::Pi() << " " << gen_b.Pt() << " " << gen_b.M() << " 0 2 0 0 0" << endl;
	 fout << "4 4 " << gen_bbar.Eta() << " " << gen_bbar.Phi()+TMath::Pi() << " " << gen_bbar.Pt() << " " << gen_bbar.M() << " 0 2 0 0 0" << endl;
	 fout << "5 6 0 " << gen_Met.Phi()+TMath::Pi() << " " << gen_Met.Pt() << " 0 0 0 0 0 0" << endl;
	}

}
	


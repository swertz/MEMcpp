#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace std;

void load_weights(const char *inputFile, vector< vector<double> > &weights);

int weights_join(const char *inputFile_weights, const char *inputFile_tree, const char *outputFile_tree){

	//gSystem->Load("libDelphes");

	cout << "Loading data..." << endl;
	
	TFile inputfile(inputFile_tree, "READ");
	TTree *Tree_input = (TTree*) inputfile.Get("Delphes");

	int allentries = Tree_input->GetEntries(); // nombre total d'événements

	cout << "Setting variables..." << endl;

	///// output variables /////////////////////////////////////////
	
	// on crée un fichier root comme output:
	TFile outputfile(outputFile_tree, "RECREATE"); // recreate permet de ne pas demander de permission lors de l'écrasement du fichier

	// on définit notre tree et les branches qui vont contenir les variables d'intérêt:
	TTree *Tree_output = Tree_input->CloneTree(0);

	double O_Weight;
	double O_Weight_Error;
	bool O_Weighted;

	Tree_output->Branch("Weight_TT", &O_Weight);
	Tree_output->Branch("Weight_TT_Error", &O_Weight_Error);
	Tree_output->Branch("Weighted_TT", &O_Weighted);
	
	///// loading weights /////////////////////////////////////////
	
  cout << "Loading weights... "  << endl;
	cout.flush();

	vector< vector<double> > weights;

	load_weights(inputFile_weights, weights);

	cout << weights.size() << " weights loaded." << endl;

	cout << "Processing " << allentries << " events..." << endl;

	int nbr_weights = weights.size();

	// loop over all events
	for(int i=0; i < allentries; i++){
		
		Tree_input->GetEntry(i);
		
		if(i < nbr_weights){
      if( weights[i][0] != i ){
        cout << "Numbers don't match!" << endl;
        exit(1);
      }
			O_Weight = weights[i][1];
			O_Weight_Error = weights[i][2];
      if(O_Weight > 0)
			  O_Weighted = true;
      else
        O_Weighted = false;
		}else{
			O_Weight = 0.;
			O_Weighted = false;
		}

		Tree_output->Fill();
	}

	cout << "Saving data..." << endl;

	// on enregistre l'arbre:
	inputfile.Close();
	
	Tree_output->Write();
	outputfile.Close();

  return 0;
}

// this function reads the output file containing the weights,
// sorts them (MW writes them in the wrong order) and stores them
// in a dynamic array

bool weightSorter(vector<double> i, vector<double> j) { return i[0] < j[0]; }

void load_weights(const char *inputFile, vector< vector<double> > &weights){

	string dummy;
	double tmp_weight, tmp_weight_error;
	int tmp_number;

	ifstream file_weights(inputFile);
	if(!file_weights){
		cerr << "Error opening input file " << inputFile << endl;
		exit(1);
	}
	
	// get first two comment lines:
	file_weights.ignore(2147483647, '\n');
	file_weights.ignore(2147483647, '\n');
	
	// loading weights
	while( !file_weights.eof() ){
    if( weights.size() % 1000 == 0 )
      cout << "Reading line " << weights.size() << endl;
		file_weights >> tmp_number >> dummy >> dummy >> tmp_weight >> tmp_weight_error; // get only number and weight & error
    vector<double> tmp_vec;
		tmp_vec.push_back(tmp_number-1);
		tmp_vec.push_back(tmp_weight);
		tmp_vec.push_back(tmp_weight_error);
    weights.push_back(tmp_vec);
	}
	weights.pop_back(); // there's one too many
  cout << "Reading line " << weights.size()-1 << endl;

	// ordering them:
  sort(weights.begin(), weights.end(), weightSorter);
	
	// making sure there's no gap:
	if( weights.size()-1 != weights.back().at(0) ){
		cout << "There are some missing weights!" << endl;
		
    for(unsigned int k=0; k<weights.size()-2; k++){
		  
      vector<double> tmp_vec;
			
      if(weights[k][0]+1 != weights[k+1][0]){
				
        cout << "Weight nr. " << k+1 << " is missing." << endl;
        
        tmp_vec.clear();
        tmp_vec.push_back(0);
        tmp_vec.push_back(0);
        tmp_vec.push_back(0);
				weights.push_back(tmp_vec);
				
        for(unsigned int j=weights.size()-1; j>k+1; j--)
					weights[j] = weights[j-1];
				
        tmp_vec.clear();
        tmp_vec.push_back(k+1);
        tmp_vec.push_back(-1);
        tmp_vec.push_back(0);
				weights[k+1] = tmp_vec;
			
      }
		
    }
	
  }

	for(unsigned int k=0; k<weights.size(); k++){
    if( k % 1000 == 0)
		  cout << "entry " << k <<": weight " << weights.at(k).at(0) << ": " << weights.at(k).at(1) << " +- " << weights.at(k).at(2) << endl;
  }
  cout << "entry " << weights.size()-1 <<": weight " << weights.back().at(0) << ": " << weights.back().at(1) << " +- " << weights.back().at(2) << endl;

	file_weights.close();
}


#include <cmath>
#include <vector>

#include "/home/fynu/swertz/memoire/scripts/utils.C"

void load_weights(const char *inputFile, vector<double> &weights);

int weights_join(const char *inputFile_weights, const char *inputFile_tree, const char *outputFile_tree){

	//gSystem->Load("~/storage/Delphes/Delphes-3.0.12/libDelphes.so");
	gSystem->Load("libDelphes");

	cout << "Loading data..." << endl;
	
	TFile inputfile(inputFile_tree, "READ");
	TTree *Tree_input = (TTree*) inputfile.Get("Delphes");

	int allentries = Tree_input->GetEntries(); // nombre total d'événements

	cout << "Setting variables..." << endl;

	///// input variables /////////////////////////////////////////

	
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
	
	///// output variables /////////////////////////////////////////
	cout << "Loading weights... ";
	cout.flush();

	vector<double> weights;
	vector<double> weightErrors;

	load_weights(inputFile_weights, weights, weightErrors);

	/*if( allentries < weights.size() ){
		cerr << endl << "The number of events in the three categories don't add up to total number of events. Exiting..." << endl;
		exit(1);
	}*/
	cout << weights.size() << " weights loaded." << endl;

	cout << "Processing " << allentries << " events..." << endl;

	int nbr_weights = weights.size();

	// loop over all events
	for(int i=0; i < allentries; i++){
		
		Tree_input->GetEntry(i);
		
		if(i < nbr_weights){
			O_Weight = weights[i];;
			O_Weight_Error = weightErrors[i];
			O_Weighted = true;
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
}

// this function reads the output file containing the weights,
// sorts them (MW writes them in the wrong order) and stores them
// in a dynamic array

void load_weights(const char *inputFile, vector<double> &weights, vector<double> &weightErrors){

	vector<int> number;
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
		file_weights >> tmp_number >> dummy >> dummy >> tmp_weight >> tmp_weight_error; // get only number and weight & error
		weights.push_back(tmp_weight);
		weightErrors.push_back(tmp_weight_error);
		number.push_back(tmp_number);
	}
	number.pop_back(); // there's one too many
	weights.pop_back();
	weightErrors.pop_back();
	if(weights.size() != number.size()){
		cerr << "Something's gone wrong in pushing back the vectors!" << endl;
		exit(1);
	}

	// ordering them:
	for(int k=1; k < weights.size(); k++){
		tmp_number = number[k];
		tmp_weight = weights[k];
		tmp_weight_error = weightErrors[k];
		int pos = k;
		while(tmp_number < number[pos-1]){
			number[pos] = number[pos-1];
			weights[pos] = weights[pos-1];
			weightErrors[pos] = weightErrors[pos-1];
			pos--;
			if(pos == 0) break;
		}
		if(pos != k){
			number[pos] = tmp_number;
			weights[pos] = tmp_weight;
			weightErrors[pos] = tmp_weight_error;
		}
	}
	
	// making sure there's no gap:
	if(weights.size()-1 != number[weights.size()-1]){
		cout << endl << "There are some missing weights!" << endl;
		for(int k=0; k<weights.size()-1; k++){
			if(number[k]+1 != number[k+1]){
				cout << "Weight nr. " << k+1 << " is missing." << endl;
				number.push_back(0);
				weights.push_back(0.);
				weightErrors.push_back(0.);
				for(int j=weights.size()-1; j>k+1; j--){
					weights[j] = weights[j-1];
					weightErrors[j] = weightErrors[j-1];
					number[j] = number[j-1];
				}
				number[k+1] = k+1;
				weights[k+1] = 0.;
			}
		}
	}

	for(int k=0; k<weights.size(); k++)
		cout << "weight " << number[k] << ": " << weights[k] << " +- " << weightErrors[k] << endl;

	file_weights.close();
}


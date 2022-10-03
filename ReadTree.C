
#include <iostream>              // main Header that defines the standard input/output stream objects:
#include "TString.h"             // provide Short String Optimization (SSO) so that short strings
#include "TClonesArray.h"        // Memory for the objects stored in the array is allocated only once in the lifetime of the clones array. All objects must be of the same class. For the rest this class has the same properties as TObjArray
#include "TProfile.h"            // A chain is a collection of files containing TTree objects.
#include "TChain.h"              // A chain is a collection of files containing TTree objects.
#include "TTree.h"               //library
#include "TFile.h"               //library
#include "TMath.h"               //library
#include "TH2.h"                 //library
#include "TH1.h"                 //library
#include "TParticle.h"           //library
#include "TParticlePDG.h"        //library

using namespace std;

bool IsLeft(const double& dphi); //function-shortcut
bool IsRight(const double& dphi); //function-shortcut
int GetPIDCode(int pdg);        //function-shortcut
double DeltaPhi(double phia, double phib, double rangeMin=-TMath::Pi()/2, double rangeMax=3*TMath::Pi()/2);

void ReadTree() {

	std::vector<TString> inputFiles;   //making a vector name input files, unlimited size

	const char* name = "Monash";       //name if the output
	const char* treeName = "TT";       //name if the TTREE

	int nCores = 1;                    //number of files
	for (int i = 0; i < nCores; ++i)   //adding files to the vector inputFiles
		inputFiles.push_back(Form("/Users/omar/Documents/Pythia/Simulations/%s_Output/AnalysisResults_%d.root",name,i));

	//@ create an instance of the TChain class
	TChain* m_chain = new TChain(treeName);           //making an object m_chain from the class TChain and name it with the output tree name

	//@ loop through the input files and add them to the chain
	for(auto inputFile : inputFiles) {              //start a loop over input file, automatics
		m_chain->Add(inputFile);                    //adding the input file into the tchain object with the name
		std::cout << "Added file: " << inputFile << std::endl; //output text
	}

	if(!m_chain) { //condition if the the was bull then show the error message
		throw std::runtime_error("Calling execute while the event loop was not initialized.");
	}

	std::cout << "Total number of Events : " << m_chain->GetEntries() << std::endl;         //getting the total number of events with the function in Tchain class

	TTree* tree = m_chain;                                                                  //initializing the TTREE object by puting it equal to tchain

	int v0m, cl1, mpi;                                                                      //define three integer variables, v0m (mult)+cl1(mult in transverese region) + number of multi parton interaction
	TClonesArray* arr = new TClonesArray("TParticle");                                      //making object arr form the class TClonesArray,

	tree->SetBranchAddress("nMPI", &mpi);                                                   //rebuilding the ttree branches with the same name as store in TTREE
	tree->SetBranchAddress("fV0M", &v0m);                                                   //rebuilding the ttree branches with the same name as store in TTREE
	tree->SetBranchAddress("fCL1", &cl1);                                                   //rebuilding the ttree branches with the same name as store in TTREE
	tree->SetBranchAddress("tracks", &arr);                                                 //rebuilding the ttree branches with the same name as store in TTREE

	int entries = tree->GetEntries();                                                       //getting numbe rof events from tree not the tchain
	cout << "Number of Events : " << entries << endl;                                       //print

	const int nPtbins = 15;                                                                 //Pt bins
	double Ptbins[nPtbins+1] = {
		0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 
		1.20, 1.4, 1.6, 1.8, 2.0, 2.5, 
		3.0, 3.5, 4.0, 5.0};

	double minNT = -0.5;                                                                //maximum and minimum bins for mult
	double maxNT = 80.5;

	TH2F* NTminvsMPI_h = new TH2F("NTminvsMPI_h","; NTmin; MPI", maxNT-minNT, minNT, maxNT, maxNT-minNT, minNT, maxNT);  //mult vs mpi
	TProfile* NTminvsMPI_p = new TProfile("NTminvsMPI_p","; NTmin; #LT MPI #GT", maxNT-minNT, minNT, maxNT, minNT, maxNT);  //tprofile of above histogram
	TH2F* NTmaxvsMPI_h = new TH2F("NTmaxvsMPI_h","; NTmax; MPI", maxNT-minNT, minNT, maxNT, maxNT-minNT, minNT, maxNT);     //max vs mpi
	TProfile* NTmaxvsMPI_p = new TProfile("NTmaxvsMPI_p","; NTmax; #LT MPI #GT", maxNT-minNT, minNT, maxNT, minNT, maxNT); //tprofile of above histogram
	TH2F* NTvsMPI_h = new TH2F("NTvsMPI_h","; NT; MPI", maxNT-minNT, minNT, maxNT, maxNT-minNT, minNT, maxNT);             //total vs mpi
	TProfile* NTvsMPI_p = new TProfile("NTvsMPI_p","; NT; #LT MPI #GT", maxNT-minNT, minNT, maxNT, minNT, maxNT);           //profile

	TH1F* hV0M = new TH1F("hV0M","Multiplicity in the V0M; #it{N}_{ch}; Entries;", 101, -0.5, 100.5);                       //hist for v0 mult
	TH1F* hCL1 = new TH1F("hCL1","Multiplicity in |#eta|<1; #it{N}_{ch}; Entries;", 101, -0.5, 100.5);                      //hist for v0 hCL1
	TH1F* hDeltaPhi = new TH1F("hDeltaPhi","; #Delta #varphi; Entries;", 200, -2.2*TMath::Pi(), 2.2*TMath::Pi());           //delta phi distribution
	hDeltaPhi->Sumw2();                                                                                                     //Create structure to store sum of squares of weights.

	TString NTnames[4] = {"NTmin","NTmax","NT","NT"};                                                                       //names for???
	TString names[3] = {"Toward","Away","Transverse"};                                                                      //names for regions

	TH1F* NT_h[4];                                                                                                          //making 4 dimensional histograms, correspond to 4 classes

	TH2F* pTpiNT_h[3];                      //for NT total                                                                                 //pt and mult for pions
	TH2F* pTkNT_h[3];                                                                                                       //pt and mult for kaons
	TH2F* pTpNT_h[3];                                                                                                       //pt and mult for kaons
	TH2F* pTk0sNT_h[3];                                                                                                     //pt and mult for kaons0
	TH2F* pTLNT_h[3];                                                                                                       //pt and mult for Lambda
	TH2F* pTXiNT_h[3];                                                                                                      //pt and mult for xi
	TH2F* pTONT_h[3];                                                                                                       //pt and mult for omega
	TH2F* pTD0NT_h[3];                                                                                                      //pt and mult for D0
	TH2F* pTLcNT_h[3];                                                                                                      //pt and mult for Lc

	TH2F* pTpiNTmin_h[3];                  //for NT_min total
	TH2F* pTkNTmin_h[3];
	TH2F* pTpNTmin_h[3];
	TH2F* pTk0sNTmin_h[3];
	TH2F* pTLNTmin_h[3];
	TH2F* pTXiNTmin_h[3];
	TH2F* pTONTmin_h[3];
	TH2F* pTD0NTmin_h[3];
	TH2F* pTLcNTmin_h[3];

	for (int i = 0; i < 4; ++i)
	{
		NT_h[i] = new TH1F(Form("%s_h",NTnames[i].Data()),Form("; %s; Entries;",NTnames[i].Data()), maxNT-minNT, minNT, maxNT);     //????
		NT_h[i]->Sumw2();
	}

	for (int i = 0; i < 3; ++i){                                                                            //initilazing the defined parameters

		pTpiNTmin_h[i] = new TH2F(Form("pTpiNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTpiNTmin_h[i]->Sumw2();
		pTkNTmin_h[i] = new TH2F(Form("pTkNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTkNTmin_h[i]->Sumw2();
		pTpNTmin_h[i] = new TH2F(Form("pTpNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTpNTmin_h[i]->Sumw2();
		pTk0sNTmin_h[i] = new TH2F(Form("pTk0sNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTk0sNTmin_h[i]->Sumw2();
		pTLNTmin_h[i] = new TH2F(Form("pTLNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTLNTmin_h[i]->Sumw2();
		pTXiNTmin_h[i] = new TH2F(Form("pTXiNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTXiNTmin_h[i]->Sumw2();
		pTONTmin_h[i] = new TH2F(Form("pTONTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTONTmin_h[i]->Sumw2();
		pTD0NTmin_h[i] = new TH2F(Form("pTD0NTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTD0NTmin_h[i]->Sumw2();
		pTLcNTmin_h[i] = new TH2F(Form("pTLcNTmin_%s_h",names[i].Data()), Form("pT Spectra (%s) vs. NTmin; #it{N}_{T}; Entries;",names[i].Data()), maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTLcNTmin_h[i]->Sumw2();

		pTpiNT_h[i] = new TH2F(Form("pTpiNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTpiNT_h[i]->Sumw2();
		pTkNT_h[i] = new TH2F(Form("pTkNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTkNT_h[i]->Sumw2();
		pTpNT_h[i] = new TH2F(Form("pTpNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTpNT_h[i]->Sumw2();
		pTk0sNT_h[i] = new TH2F(Form("pTk0sNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTk0sNT_h[i]->Sumw2();
		pTLNT_h[i] = new TH2F(Form("pTLNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTLNT_h[i]->Sumw2();
		pTXiNT_h[i] = new TH2F(Form("pTXiNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTXiNT_h[i]->Sumw2();
		pTONT_h[i] = new TH2F(Form("pTONT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTONT_h[i]->Sumw2();
		pTD0NT_h[i] = new TH2F(Form("pTD0NT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTD0NT_h[i]->Sumw2();
		pTLcNT_h[i] = new TH2F(Form("pTLcNT_%s_h",names[i].Data()),"; #it{p}_{T} (GeV/#it{c}); Entries;", maxNT-minNT, minNT, maxNT, nPtbins, Ptbins);
		pTLcNT_h[i]->Sumw2();

	}

	int Nevents = 0;                                                                        //number of events
	for (int i = 0; i < entries; i++)                                                       //loop over events
	{

		if (i%1000000==0.0) { cout << i << endl; }

		tree->GetEntry(i);                                                                  //loging into the branhces of ith events

		hV0M->Fill(v0m);                                                                    //fill histogram in the initialize branch
		hCL1->Fill(cl1);                                                                    //fill histogram in the initialize branch

		int lines = arr->GetEntries();
		cout << "Numer of Particles in Event " << i << " =  " << lines << endl;             //

		//@ Here is the implementation of the RT analysis
		int index_leading = -1.0;
		double pt_leading = -1.0;
		double phi_leading = -1.0;
		
		//@ This for loop find the leading particle
		//@ This means, tis pT and index
		for (int j = 0; j < lines; ++j ) {

			TParticle* particle = (TParticle*)arr->At(j);

			double pt = particle->Pt();
			double eta = particle->Eta();
			double phi = particle->Phi();
			double charge = particle->GetPDG()->Charge();

			if (TMath::Abs(charge) == 0.0) { continue; }
			if (TMath::Abs(eta) > 0.8) { continue; }

			if (pt_leading < pt){
				pt_leading = pt;
				phi_leading = phi;
				index_leading = j;
			}

		} //@ Leading Particle loop

		if (pt_leading < 5.0) { continue; }
		Nevents++;
		//std::cout << "pT leading= " << pt_leading << endl;

		int NT = 0;
		int NTmin = 0;
		int NTmax = 0;
		int NTequal = 0;
		int NTleft = 0;
		int NTright = 0;

		for (int j = 0; j < lines; ++j ) 
		{

			TParticle* particle = (TParticle*)arr->At(j);
			if (!particle) { continue; }

			//@ Avoids auto-correlation 
			if (index_leading==j) { continue; }

			double pt = particle->Pt();
			double eta = particle->Eta();
			double phi = particle->Phi();
			double charge = particle->GetPDG()->Charge();

			if (TMath::Abs(charge) == 0.0) { continue; }

			double dphi = DeltaPhi(phi_leading,phi); 
			hDeltaPhi->Fill(dphi);

			double pi = TMath::Pi();
			if (TMath::Abs(dphi) < pi/3.) { continue; } //@ Toward region 
			else if (TMath::Abs(dphi-pi) < pi/3.) { continue; } //@ Away region
			else {  
				if ( IsLeft(dphi) ) { NTleft++; } 
				if ( IsRight(dphi) ) { NTright++; } 
			}

		} //@ Particle loop

		//@ How many events have equal number of particles on the left and right sides
		if (NTleft == NTright) { NTequal = NTleft; NT_h[2]->Fill(NTequal); }   
		//@ Check if either the left or right side has the minimum number of particles
		if (NTleft < NTright) { NTmin = NTleft; NTmax = NTright; NT_h[0]->Fill(NTmin); NT_h[1]->Fill(NTmax); }
		if (NTleft > NTright) { NTmin = NTright; NTmax = NTleft; NT_h[0]->Fill(NTmin); NT_h[1]->Fill(NTmax); }

		if (NTleft == NTright) { NT = 2.0*NTequal; }
		else { NT = NTleft + NTright; }
		NT_h[3]->Fill(NT);

		NTvsMPI_h->Fill(NT,mpi);
		NTvsMPI_p->Fill(NT,mpi);

		for (int j = 0; j < lines; ++j ) 
		{

			TParticle* particle = (TParticle*)arr->At(j);
			if (!particle) { continue; }

			//@ Avoid auto-correlations
			if (index_leading==j) { continue; }

			double pt = particle->Pt();
			double phi = particle->Phi();
			int pid_code = GetPIDCode(particle->GetPdgCode());
			if (pid_code < 0) { continue; }

			double dphi = DeltaPhi(phi_leading,phi); 
			double pi = TMath::Pi();

			int index_region = -1;
			if (TMath::Abs(dphi) < pi/3.) { //@ Toward region 
				index_region = 0;
				if (pid_code==0) pTpiNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==1) pTkNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==2) pTpNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==3) pTk0sNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==4) pTLNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==5) pTXiNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==6) pTONT_h[index_region]->Fill(NT,pt); 
				if (pid_code==7) pTD0NT_h[index_region]->Fill(NT,pt); 
				if (pid_code==8) pTLcNT_h[index_region]->Fill(NT,pt); 
			}
			else if (TMath::Abs(dphi-pi) < pi/3.) { //@ Away region
				index_region = 1;
				if (pid_code==0) pTpiNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==1) pTkNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==2) pTpNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==3) pTk0sNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==4) pTLNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==5) pTXiNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==6) pTONT_h[index_region]->Fill(NT,pt); 
				if (pid_code==7) pTD0NT_h[index_region]->Fill(NT,pt); 
				if (pid_code==8) pTLcNT_h[index_region]->Fill(NT,pt); 
			}
			else{ //@ Transverse region
				index_region = 2; 
				if (pid_code==0) pTpiNT_h[index_region]->Fill(NT,pt); 
				if (pid_code==1) pTkNT_h[index_region]->Fill(NT,pt);
				if (pid_code==2) pTpNT_h[index_region]->Fill(NT,pt);
				if (pid_code==3) pTk0sNT_h[index_region]->Fill(NT,pt);
				if (pid_code==4) pTLNT_h[index_region]->Fill(NT,pt);
				if (pid_code==5) pTXiNT_h[index_region]->Fill(NT,pt);
				if (pid_code==6) pTONT_h[index_region]->Fill(NT,pt);
				if (pid_code==7) pTD0NT_h[index_region]->Fill(NT,pt);
				if (pid_code==8) pTLcNT_h[index_region]->Fill(NT,pt);
			}
		} //@ Particle loop

		//@ Skip events that have equal number of particles
		//@ either of the left or right regions
		if (NTleft == NTright) { continue; } 

		NTminvsMPI_h->Fill(NTmin,mpi);
		NTminvsMPI_p->Fill(NTmin,mpi);
		NTmaxvsMPI_h->Fill(NTmax,mpi);
		NTmaxvsMPI_p->Fill(NTmax,mpi);

		for (int j = 0; j < lines; ++j ) 
		{

			TParticle* particle = (TParticle*)arr->At(j);
			if (!particle) { continue; }

			if (index_leading==j) { continue; }

			double pt = particle->Pt();
			double phi = particle->Phi();
			int pid_code = GetPIDCode(particle->GetPdgCode());
			if (pid_code < 0) { continue; }

			double dphi = DeltaPhi(phi_leading,phi); 
			double pi = TMath::Pi();

			int index_region = -1;
			if (TMath::Abs(dphi) < pi/3.) { 
				index_region = 0; 
				if (pid_code==0) pTpiNTmin_h[index_region]->Fill(NTmin,pt);
				if (pid_code==1) pTkNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==2) pTpNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==3) pTk0sNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==4) pTLNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==5) pTXiNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==6) pTONTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==7) pTD0NTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==8) pTLcNTmin_h[index_region]->Fill(NTmin,pt); 
			}
			if (TMath::Abs(dphi-pi) < pi/3.) { 
				index_region = 1; 
				if (pid_code==0) pTpiNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==1) pTkNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==2) pTpNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==3) pTk0sNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==4) pTLNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==5) pTXiNTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==6) pTONTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==7) pTD0NTmin_h[index_region]->Fill(NTmin,pt); 
				if (pid_code==8) pTLcNTmin_h[index_region]->Fill(NTmin,pt); 
			}
			if (NTleft==NTmin) {
				index_region = 2;
				if ( IsLeft(dphi) ) { 
					if (pid_code==0) pTpiNTmin_h[index_region]->Fill(NTmin,pt);
					if (pid_code==1) pTkNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==2) pTpNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==3) pTk0sNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==4) pTLNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==5) pTXiNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==6) pTONTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==7) pTD0NTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==8) pTLcNTmin_h[index_region]->Fill(NTmin,pt); 
				}
			}
			if (NTright==NTmin) {
				index_region = 2;
				if ( IsRight(dphi) ) { 
					if (pid_code==0) pTpiNTmin_h[index_region]->Fill(NTmin,pt);
					if (pid_code==1) pTkNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==2) pTpNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==3) pTk0sNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==4) pTLNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==5) pTXiNTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==6) pTONTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==7) pTD0NTmin_h[index_region]->Fill(NTmin,pt); 
					if (pid_code==8) pTLcNTmin_h[index_region]->Fill(NTmin,pt); 
				} 
			}
		} //@ Particle loop

	} //@ Event loop

	TFile* fOut = new TFile("./OutputTree.root", "recreate");
	fOut->cd();

	NTminvsMPI_h->Write();
	NTminvsMPI_p->Write();
	NTmaxvsMPI_h->Write();
	NTmaxvsMPI_p->Write();
	NTvsMPI_h->Write();
	NTvsMPI_p->Write();
	hV0M->Write();
	hCL1->Write();
	for (int i = 0; i < 4; ++i)
		NT_h[i]->Write();

	hDeltaPhi->Write();

	for (int i = 0; i < 3; ++i)
	{
		pTpiNTmin_h[i]->Write();
		//	pTkNTmin_h[i]->Write();
		pTpNTmin_h[i]->Write();
		//	pTk0sNTmin_h[i]->Write();
		pTLNTmin_h[i]->Write();
		pTXiNTmin_h[i]->Write();
		pTONTmin_h[i]->Write();
		pTD0NTmin_h[i]->Write();
		pTLcNTmin_h[i]->Write();
	}

	for (int i = 0; i < 3; ++i)
	{
		pTpiNT_h[i]->Write();
		//	pTkNT_h[i]->Write();
		pTpNT_h[i]->Write();
		//	pTk0sNT_h[i]->Write();
		pTLNT_h[i]->Write();
		pTXiNT_h[i]->Write();
		pTONT_h[i]->Write();
		pTD0NT_h[i]->Write();
		pTLcNT_h[i]->Write();
	}


	fOut->Close();
	delete fOut;

}

double DeltaPhi(double phia, double phib, double rangeMin, double rangeMax)
{
	double dphi = -999;
	//	double pi = TMath::Pi();

	dphi = phia - phib;

	return dphi;
}

int GetPIDCode(int pdgCode)
{
	int pidCode = -1;

	switch (TMath::Abs(pdgCode)) {
		case 211:
			pidCode = 0; // pion
			break;
		case 321:
			pidCode = 1; // kaon
			break;
		case 2212:
			pidCode = 2; // proton
			break;
		case 310:
			pidCode = 3; // K0s
			break;
		case 3122:
			pidCode = 4; // Lambda
			break;
		case 3312:
			pidCode = 5; // Xi-
			break;
		case 3334:
			pidCode = 6; // Omega-
			break;
		case 421:
			pidCode = 7; // D0
			break;
		case 4122:
			pidCode = 8; // Lambda c
			break;
		default:
			break;
	};

	return pidCode;

}

bool IsLeft(const double& dphi) 
{
	bool isLeft = kFALSE;
	if ( -2.0*TMath::Pi()/3. <= dphi && dphi <= -1.0*TMath::Pi()/3.0 ) { isLeft = kTRUE; }				 

	return isLeft;
}

bool IsRight(const double& dphi) 
{
	bool isRight = kFALSE;
	if ( TMath::Pi()/3.0 <= dphi && dphi <= 2.0*TMath::Pi()/3.0 ) { isRight = kTRUE; }

	return isRight;
}

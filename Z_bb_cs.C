#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TImage.h"
#include "TROOT.h"
#include "fstream"
#include "string"
#include "sstream"
#include "iostream"
#include "iomanip"
#include "vector"
#include "TLegend.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TBranch.h"
#include "THStack.h"
#include "TLatex.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TObject.h"
#include <array>
#include "TVectorD.h"
#include "TVectorT.h"
#include "TMath.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TAttText.h"
#include "TClonesArray.h"
// #include "../../NeroProducer/Core/interface/BareJets.hpp"

using namespace std;

Bool_t Two_Gen_b_Quarks_Matched_to_Z( std::vector<Double_t> *genPdgId , std::vector<Double_t> *genMotherPdgId )
{
	Bool_t two_gen_b_quarks_matched_to_Z ;
	int n_b_matched_to_Z = 0 ;
	// Count how many matches you get
	for (int i = 0; i < genPdgId->size(); ++i)
	{
		if ((abs(genPdgId->at(i))==5) && ( genMotherPdgId->at(i) == 23 ))
		{
			n_b_matched_to_Z = n_b_matched_to_Z + 1 ;
		}
	}
	// Determnie if signal or background
	if (n_b_matched_to_Z < 2 )
	{
		two_gen_b_quarks_matched_to_Z = kFALSE ;
	}
	else if ( n_b_matched_to_Z == 2 )
	{
		two_gen_b_quarks_matched_to_Z = kTRUE ;
		// cout << "Signal event found" << endl;
	}
	else if ( n_b_matched_to_Z > 2 )
	{
		two_gen_b_quarks_matched_to_Z = kFALSE ;
		cout << " !!! More than 2 b quarks matched to a Z !!! " << endl ;
	}
	else 
	{
		cout << "*** Error in finding the Z->bb signal ***" << endl ;
	}
return two_gen_b_quarks_matched_to_Z;
}

THStack* Plot_Signal_and_Backgrounf(TH1D* h_sig , TH1D* h_bkg, TString variable )
{
	THStack* s_sig_and_bkg = new THStack("h_sig_bkg_" + variable , variable) ;
	TLegend* leg = new TLegend(0.7, 0.9, 0.90, 0.7) ;
	h_sig->SetLineColor(kRed) ;
	h_bkg->SetLineColor(kBlue) ;
	leg->AddEntry(h_bkg, "Background" , "l") ;
	leg->AddEntry(h_sig, "Signal" , "l") ;
	TCanvas* c = new TCanvas() ; 
	s_sig_and_bkg->Add(h_sig) ;
	s_sig_and_bkg->Add(h_bkg) ;
	s_sig_and_bkg->Draw("nostack, hist") ;
	leg->Draw();
	c->SaveAs("s_sig_and_bkg_.png") ;
	c->SaveAs("s_sig_and_bkg_.pdf") ;
return s_sig_and_bkg;
}


int Get_Cut_Parameters(std::vector<TTree*> *v_trees , TString variable)
{
	cout << "Determing the " << variable <<  " cut" << endl ;
	TTree* tree = new TTree();

	// Make histograms to store the signal and backgrounds regions
	Int_t nbins = 9 ; Double_t xmin = 0 ; Double_t xmax = 10 ;
	TH1D* h_sig = new TH1D("h_sig", "h_sig", nbins, xmin, xmax);
	TH1D* h_bkg = new TH1D("h_bkg", "h_bkg", nbins, xmin, xmax);



	// loop over the trees
	for (int i = 0; i < v_trees->size(); ++i)
		{
			tree = v_trees->at(i) ;	

			std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
			std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
			std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
			TClonesArray *jetP4 = new TClonesArray();
			TClonesArray *genP4 = new TClonesArray();
			float mcWeight;

			tree->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId) ;
			cout << "Located Branch: jetMotherPdgId"<< endl;
			tree->SetBranchAddress("genPdgId", &genPdgId) ;
			cout << "Located Branch: genPdgId"<< endl;
			tree->SetBranchAddress("genMotherPdgId", &genMotherPdgId) ;
			cout << "Located Branch: genMotherPdgId"<< endl;
			tree->SetBranchAddress("jetP4", &jetP4) ;
			cout << "Located Branch: jetP4"<< endl;
			tree->SetBranchAddress("genP4", &genP4) ;
			cout << "Located Branch: gen_b_P4"<< endl;
			tree->SetBranchAddress("mcWeight",&mcWeight);
			cout << "Located Branch: mcWeight" << endl;
			
			Double_t nentries = tree->GetEntries() ;
			// Decrease the number of entries for testing to 0.05%
			nentries = nentries * 0.05 ;
			cout << "Looping over " << nentries << " events in dataset " << i << endl ;
			for (int j = 0; j < nentries; ++j)
			{
				tree->GetEntry(j) ;
				int x = 0 ;
				if ( variable == "njets")
					{	
						x = jetP4->GetEntries() ;
					}
				// cout << "njets = " << njets << endl ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
				// Check to see if there are 2 b quarks that are matched to a Z
				if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId ) == kTRUE)
				{
					//fill the signal histogram
					h_sig->Fill( x , mcWeight ) ;
				}
				else
				{
					// fill the background histogram
					h_bkg->Fill( x , mcWeight ) ;
				}
			}
		}
		// Get the integrals of the signal and background
		Double_t sig_integral = 0 ;  
		Double_t bkg_integral = 0 ; 
		sig_integral = h_sig->Integral(); cout << "Signal Integral = " << sig_integral << " Signal Entries = " << h_sig->GetEntries() << endl ;
		bkg_integral = h_bkg->Integral(); cout << "Background Integral = " << bkg_integral <<  " Background Entries = " << h_bkg->GetEntries() <<  endl ;
		THStack* s_sig_and_bkg_njets =  Plot_Signal_and_Backgrounf( h_sig , h_bkg, "njets" ) ;
return 0 ;
}


// Main Function
void Z_bb_cs()
{
	std::vector<TTree*> *v_trees = new std::vector<TTree*> ;
	TFile *f_TTJets = new TFile("~/cms/hist/TTJets.root") ;
   	cout << "Found file TTJets.root" << endl;
   	TTree *t_TTJets = (TTree*) f_TTJets->FindObjectAny("events") ; v_trees->push_back(t_TTJets) ;

   	TFile *f_WJetsToLNu = new TFile("~/cms/hist/WJetsToLNu.root") ;
   	cout << "Found file WJetsToLNu.root" << endl;
   	TTree* t_WJetsToLNu = (TTree*) f_WJetsToLNu->FindObjectAny("events") ; v_trees->push_back(t_WJetsToLNu) ;

   	TFile *f_ZZTo2Q2Nu = new TFile("~/cms/hist/ZZTo2Q2Nu.root") ;
   	cout << "Found file ZZTo2Q2Nu.root" << endl;
   	TTree* t_ZZTo2Q2Nu = (TTree*) f_ZZTo2Q2Nu->FindObjectAny("events") ; v_trees->push_back(t_ZZTo2Q2Nu) ;

   	TFile *f_WZTo1L1Nu2Q = new TFile("~/cms/hist/WZTo1L1Nu2Q.root") ;
   	cout << "Found file WZTo1L1Nu2Q.root" << endl ;
   	TTree* t_WZTo1L1Nu2Q = (TTree*) f_WZTo1L1Nu2Q->FindObjectAny("events") ; v_trees->push_back(t_WZTo1L1Nu2Q) ;

   	TFile *f_WWToNuQQ = new TFile("~/cms/hist/WWToNuQQ.root") ;
   	cout << "Found file WWToNuQQ.root" << endl;
   	TTree* t_WWToNuQQ = (TTree*) f_WWToNuQQ->FindObjectAny("events") ; v_trees->push_back(t_WWToNuQQ) ;

   	// Plot the nJets distributions in the signal and background region to see what is a good cut to apply
   	//return the desired njets cut value
   	int njets_cut = Get_Cut_Parameters( v_trees , "njets" );
}
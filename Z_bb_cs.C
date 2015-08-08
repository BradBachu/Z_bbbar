#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
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
#include <math.h>
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
	h_sig->SetFillColor(kRed) ;
	h_bkg->SetLineColor(kBlue) ;
	h_bkg->SetFillColor(kBlue) ;
	leg->AddEntry(h_bkg, "Background" , "l") ;
	leg->AddEntry(h_sig, "Signal" , "l") ;
	TCanvas* c = new TCanvas() ; 
	s_sig_and_bkg->Add(h_sig) ;
	s_sig_and_bkg->Add(h_bkg) ;
	s_sig_and_bkg->Draw("nostack, hist") ;
	leg->Draw();
	c->SaveAs("s_sig_and_bkg_"+variable+".png") ;
	c->SaveAs("s_sig_and_bkg_"+variable+".pdf") ;
return s_sig_and_bkg;
}

Bool_t Two_Jets(TClonesArray *jetP4)
{
	Bool_t two_jets;
	int n_jets = jetP4->GetEntries() ;
	if (n_jets > 1)
	{
		two_jets = kTRUE ;
	}
	else
	{
		two_jets = kFALSE ;
	}
return two_jets ;
}

THStack* Plot_Dijet_Mass(TH1D* h_sig , TH1D* h_bkg, Double_t cut1 , Double_t cut2)
{
	THStack* s_sig_and_bkg = new THStack(TString::Format("Dijet Mass/_%.3f/_%3f",cut1, cut2) , TString::Format("Dijet Mass _%.3f/_%3f",cut1, cut2)) ;
	TLegend* leg = new TLegend(0.7, 0.9, 0.90, 0.7) ;
	h_sig->SetLineColor(kRed) ;
	h_bkg->SetLineColor(kBlue) ;
	leg->AddEntry(h_bkg, "Background" , "l") ;
	leg->AddEntry(h_sig, "Signal" , "l") ;
	TCanvas* c = new TCanvas() ; 
	s_sig_and_bkg->Add(h_sig) ;
	s_sig_and_bkg->Add(h_bkg) ;
	s_sig_and_bkg->Draw() ;
	leg->Draw();
	c->SaveAs(TString::Format("h_sig_dijet_mass_%.3f_%.3f.png",cut1,cut2)) ;
	c->SaveAs(TString::Format("h_sig_dijet_mass_%.3f_%.3f.pdf",cut1,cut2)) ;
	// c->SaveAs("s_sig_and_bkg_"+".pdf") ;
return s_sig_and_bkg;
}

// Get the Z pt fromt the reconstruction from jets
Double_t Get_dijet_mass( TClonesArray* jetP4 )
{
	Double_t n_jets = jetP4->GetEntries() ;
	Double_t mass = 0 ;
	Double_t true_Z_mass = 91.1876 ;
	Double_t Best_Z_mass = 0 ;
	Double_t Z_pt = 0;
	if ( Two_Jets(jetP4) == kFALSE ) {return Z_pt ;}

	int i_best_jet_1 = 0 ;
	int i_best_jet_2 = 0 ;

	// cout << "n_jets = " << n_jets << endl;

	for (int i = 0; i < n_jets; ++i)
	{
		TLorentzVector *jet1 = new TLorentzVector();
		jet1 = dynamic_cast<TLorentzVector*>(jetP4->At(i)) ;
		for (int j = i+1 ; j < n_jets; ++j)
		{

			TLorentzVector* jet2 = new TLorentzVector(); 
			jet2= dynamic_cast<TLorentzVector*>(jetP4->At(j)) ;
			TLorentzVector* Z = new TLorentzVector();
			*Z = *jet1 + *jet2 ;
			mass = Z->M() ;
			// cout << " Z mass = " << mass << " for pair i = " << i << " j = " << j << endl ;
			if ( abs(mass - true_Z_mass) < abs(Best_Z_mass - true_Z_mass) )
			{
				//this means that this is the best combination of jets that came from a Z
				// store the indexes 
				Best_Z_mass = mass ;
				i_best_jet_1 = i ;
				i_best_jet_2 = j ;
			}
			// cout << "Best Mass = " << Best_Z_mass<< " at i = " << i << " j = " << j << endl;
		}
	}
	// cout << "Got best Z mass" << endl;

	//now we should have the best Jets to reconstruct the Z from
	TLorentzVector *best_jet1 = dynamic_cast<TLorentzVector*>(jetP4->At(i_best_jet_1)) ;
	TLorentzVector *best_jet2 = dynamic_cast<TLorentzVector*>(jetP4->At(i_best_jet_2)) ;
	TLorentzVector *Z = new TLorentzVector() ;
	*Z = *best_jet1 + *best_jet2 ;
	Double_t Z_mass = Z->M() ;

return Z_mass ;
}



Double_t Get_Cut_Parameters(std::vector<TTree*> *v_trees , TString variable , Int_t nbins , Double_t xmin , Double_t xmax)
{
	cout << "Determing the " << variable <<  " cut" << endl ;
	TTree* tree = new TTree();

	// Make histograms to store the signal and backgrounds regions
	// Int_t nbins = 9 ; Double_t xmin = 0 ; Double_t xmax = 10 ;
	TH1D* h_sig = new TH1D("h_sig", "h_sig", nbins, xmin, xmax);
	TH1D* h_bkg = new TH1D("h_bkg", "h_bkg", nbins, xmin, xmax);



	// loop over the trees
	for (int i = 0; i < v_trees->size(); ++i)
		{
			tree = v_trees->at(i) ;	

			std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
			std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
			std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
			TClonesArray *jetP4 = new TClonesArray() ;
			TClonesArray *genP4 = new TClonesArray() ; 
			TClonesArray *metP4 = new TClonesArray() ;
			float mcWeight;

			tree->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId) ;
			cout << "Located Branch: jetMotherPdgId"<< endl ;
			tree->SetBranchAddress("genPdgId", &genPdgId) ;
			cout << "Located Branch: genPdgId"<< endl ;
			tree->SetBranchAddress("genMotherPdgId", &genMotherPdgId) ;
			cout << "Located Branch: genMotherPdgId"<< endl ;
			tree->SetBranchAddress("jetP4", &jetP4) ;
			cout << "Located Branch: jetP4"<< endl ;
			tree->SetBranchAddress("genP4", &genP4) ;
			cout << "Located Branch: gen_b_P4"<< endl;
			tree->SetBranchAddress("mcWeight",&mcWeight) ;
			cout << "Located Branch: mcWeight" << endl ;
			tree->SetBranchAddress("metP4", &metP4) ;
			cout << "Located Branch: metP4" << endl ;
			
			Double_t nentries = tree->GetEntries() ;
			// Decrease the number of entries for testing to 0.05%
			nentries = nentries * 0.05 ;
			cout << "Looping over " << nentries << " events in dataset " << i << endl ;
			for (int j = 0; j < nentries; ++j)
			{
				tree->GetEntry(j) ;
				int x = 0 ;
				if      ( variable == "njets") {	x = jetP4->GetEntries() ;}
				else if ( variable == "met")   { TLorentzVector* MetP4 = dynamic_cast<TLorentzVector*>(metP4->At(0)) ;    x = MetP4->Pt() ;}  
				else {cout << " !!! Incorrect choice of variable !!!" << endl;}
				// cout << variable << " = " << x << endl ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
				// Check to see if there are 2 b quarks that are matched to a Z
				if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId ) == kTRUE)
				{
					//fill the signal histogram
					h_sig->Fill( Get_dijet_mass(jetP4) , mcWeight ) ;
				}
				else
				{
					// fill the background histogram
					h_bkg->Fill( Get_dijet_mass(jetP4) , mcWeight ) ;
				}
			}
		}
		// Get the integrals of the signal and background
		Double_t sig_integral = 0 ;  
		Double_t bkg_integral = 0 ; 
		sig_integral = h_sig->Integral(); cout << "Signal Integral = " << sig_integral << " Signal Entries = " << h_sig->GetEntries() << endl ;
		bkg_integral = h_bkg->Integral(); cout << "Background Integral = " << bkg_integral <<  " Background Entries = " << h_bkg->GetEntries() <<  endl ;
		THStack* s_sig_and_bkg_njets =  Plot_Signal_and_Backgrounf( h_sig , h_bkg, variable ) ;
return 0 ;
}

Double_t GetSensitivity( Double_t sig_integral , Double_t bkg_integral )
{
	Double_t sensitivity = 0 ;
	sensitivity = (sig_integral / sqrt( bkg_integral ) );
	return sensitivity;
}

// The goal of this function is to scan over my cuts and determine what are the best values to place to decrease sensitivity
void Minimize_Sensitivity( std::vector<TTree*> *v_trees, 
										    std::vector<TString> *v_cut_names, 
										    std::vector<Double_t> *v_cut_min, 
										    std::vector<Double_t> *v_cut_max , 
										    std::vector<Int_t> *v_cut_increments)
{	

	cout << "Plotting the Sensitivity as a function of njets and met" << endl;
 	// Make the 2D hist to store the sensitivitiy each time
 	// X-Axis = njets
	// Y-Axis = met
 	TH2D* h_2D_sensitivity = new TH2D("h_2D_sensitivity", "Sensitivitiy", 
 										v_cut_increments->at(0) ,  v_cut_min->at(0) , v_cut_max->at(0) , 
 										v_cut_increments->at(1) , v_cut_min->at(1) , v_cut_max->at(1) ) ;

 	Double_t cut1 = 0 ;
 	Double_t cut2 = 0 ;
 	Double_t increment1 = 0 ;
 	Double_t increment2 = 0 ;
 	increment1 = ( ( v_cut_max->at(0) - v_cut_min->at(0) ) / v_cut_increments->at(0) ) ;
 	increment2 = ( ( v_cut_max->at(1) - v_cut_min->at(1) ) / v_cut_increments->at(1) ) ;
	//loop over the regions. This requires n nests
	 for (int x = 0; x < v_cut_increments->at(0); ++x)
	 {
	 	//njets
	 	cut1 = v_cut_min->at(0) + (x *  increment1 ) ;
	 
	 	for (int y = 0; y < v_cut_increments->at(1); ++y)
	 	{	
	 		//met
		 	cut2 = v_cut_min->at(1) + (y *  increment2 ) ;
		 	cout << "Getting sensitivity for njets >= " << cut1 << " and met >=" << cut2 << endl;

			// Make histograms to store the signal and backgrounds regions per cut criteria
			Int_t nbins = 100 ; Double_t xmin = 60 ; Double_t xmax = 160 ;
			TH1D* h_sig_dijet_mass = new TH1D("h_sig", "h_sig", nbins, xmin, xmax);
			TH1D* h_bkg_dijet_mass = new TH1D("h_bkg", "h_bkg", nbins, xmin, xmax);
			TTree* tree = new TTree();

			// loop over the trees
			for (int i = 0; i < v_trees->size(); ++i)
				{
					tree = v_trees->at(i) ;	

					std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
					std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
					std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
					TClonesArray *jetP4 = new TClonesArray() ;
					TClonesArray *genP4 = new TClonesArray() ; 
					TClonesArray *metP4 = new TClonesArray() ;
					float mcWeight;

					tree->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId) ;
					tree->SetBranchAddress("genPdgId", &genPdgId) ;
					tree->SetBranchAddress("genMotherPdgId", &genMotherPdgId) ;
					tree->SetBranchAddress("jetP4", &jetP4) ;
					tree->SetBranchAddress("genP4", &genP4) ;
					tree->SetBranchAddress("mcWeight",&mcWeight) ;
					tree->SetBranchAddress("metP4", &metP4) ;
					Double_t nentries = tree->GetEntries() ;
					nentries = nentries;
					for (int j = 0; j < nentries; ++j)
					{
						tree->GetEntry(j) ;
						int x = 0 ;
						// Apply cuts
						Int_t njets = 0; Double_t met = 0 ;
						njets = jetP4->GetEntries() ;
						TLorentzVector* MetP4 = dynamic_cast<TLorentzVector*>(metP4->At(0))  ;
						met = MetP4->Pt() ;
						if (!(njets >= cut1 )) continue ;
						if (!(met >= cut2 )) continue ;
						if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
						// Check to see if there are 2 b quarks that are matched to a Z
						if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId ) == kTRUE) 
						{
							//fill the signal histogram
							h_sig_dijet_mass->Fill( Get_dijet_mass(jetP4) , mcWeight ) ;
						}
						else
						{
							// fill the background histogram
							h_bkg_dijet_mass->Fill( Get_dijet_mass(jetP4) , mcWeight ) ;
						}
					}
				}
			THStack* s_sig_and_bkg_dijetmass =  Plot_Dijet_Mass( h_sig_dijet_mass , h_bkg_dijet_mass , cut1, cut2 ) ;
			Double_t sig_integral = 0 ;  
			Double_t bkg_integral = 0 ; 
			sig_integral = h_sig_dijet_mass->Integral() ;
			bkg_integral = h_bkg_dijet_mass->Integral() ;
			Double_t sensitivity = 0 ;
			sensitivity = GetSensitivity( sig_integral , bkg_integral ) ;
			cout << "Filling hist " << x+1 << " " << y+1 << " with sensitivity " << sensitivity << endl;
			if (sensitivity > 0)
			{
				h_2D_sensitivity->Fill( x+1 , y+1 , sensitivity*1000);
			}
			else
			{
				h_2D_sensitivity->SetBinContent( x+1 , y+1 , 0);
			}
		}
	}
	h_2D_sensitivity->SetMinimum(0);
	TCanvas* c1 = new TCanvas() ;
	h_2D_sensitivity->Draw() ;
	c1->SaveAs("surface1.png") ;
	c1->SaveAs("surface1.pdf") ;
	TCanvas* c2 = new TCanvas() ;
	h_2D_sensitivity->Draw("E") ;
	c2->SaveAs("surface2.png") ;
	c2->SaveAs("surface2.pdf") ;
	TCanvas* c3 = new TCanvas() ;
	h_2D_sensitivity->Draw("LEGO") ;
	c3->SaveAs("surface3.png") ;
	c3->SaveAs("surface3.pdf") ;
	TCanvas* c4 = new TCanvas() ;
	h_2D_sensitivity->Draw("SURF") ;
	c4->SaveAs("surface4.png") ;
	c4->SaveAs("surface4.pdf") ;
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
   	// Double_t njets_cut = Get_Cut_Parameters( v_trees , "njets" , 9, 0 , 10 ) ;
   	// Double_t met_cut = Get_Cut_Parameters(v_trees , "met" , 100, 0, 1000) ;

   	// information on each cut should be packed in the form Name,Min,Max,Increments
   	std::vector<TString>  *v_cut_names = new std::vector<TString> ; v_cut_names->push_back("njets") ; v_cut_names->push_back("met") ;
   	std::vector<Double_t> *v_cut_min = new std::vector<Double_t> ;
   	std::vector<Double_t> *v_cut_max = new std::vector<Double_t> ;
   	std::vector<Int_t> *v_cut_increments = new std::vector<Int_t> ;
   	v_cut_names->push_back("njets") ;  v_cut_min->push_back(0) ; v_cut_max->push_back(10) ; v_cut_increments->push_back(10) ;
   	v_cut_names->push_back("met") ;  v_cut_min->push_back(0) ;  v_cut_max->push_back(1000); v_cut_increments->push_back(20) ;

   	Minimize_Sensitivity(v_trees, v_cut_names, v_cut_min, v_cut_max , v_cut_increments) ;

}
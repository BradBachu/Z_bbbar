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
#include "TGraphAsymmErrors.h"
#include <math.h>
// #include "../../NeroProducer/Core/interface/BareJets.hpp"

using namespace std;

Bool_t Two_Gen_b_Quarks_Matched_to_Z( std::vector<Double_t> *genPdgId , std::vector<Double_t> *genMotherPdgId  , TClonesArray *genP4 , Double_t b_pt_cut)
{
	Bool_t two_gen_b_quarks_matched_to_Z ;
	int n_b_matched_to_Z = 0 ;
	// Count how many matches you get
	for (unsigned int i = 0; i < genPdgId->size(); ++i)
	{
		cout << i <<" : genPdgId = " << genPdgId->at(i) << " and genMotherPdgId = " << genMotherPdgId->at(i) << endl ;
		if ((abs(genPdgId->at(i))==5) && ( genMotherPdgId->at(i) == 23 ))
		{
			TLorentzVector* genbP4 = dynamic_cast<TLorentzVector*> (genP4->At(i));
			if (b_pt_cut != 0)
			{
				if (genbP4->Pt() < b_pt_cut) continue ;
			}
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

THStack* Plot_Signal_and_Background(TH1D* h_sig , TH1D* h_bkg, TString variable )
{
	THStack* s_sig_and_bkg = new THStack("h_sig_bkg_" + variable , variable) ;
	TLegend* leg = new TLegend(0.7, 0.9, 0.90, 0.7) ;
	h_sig->SetLineColor(kRed) ;
	// h_sig->SetFillColor(kRed) ;
	h_bkg->SetLineColor(kBlue) ;
	// h_bkg->SetFillColor(kBlue) ;
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

// Used to get the Z pt from gen level. this in this case from two b's
Double_t Get_Z_pT_gen_level( std::vector<Double_t> *genPdgId , std::vector<Double_t> *genMotherPdgId , TClonesArray* genP4 )
{
	// loop through and find what index of the b matched to the Z
	std::vector<Int_t> *v_matched_b_to_Z = new std::vector<Int_t> ;
	for (unsigned int i = 0; i < genPdgId->size(); ++i)
	{
		if( (abs(genPdgId->at(i)) == 5 ) && (genMotherPdgId->at(i) == 23))
		{
			v_matched_b_to_Z->push_back(i);
		}
	}
	// now that we have the b's that came from a Z we can reconstruct them to get the Z
	if (v_matched_b_to_Z->size() > 2) 
	{
		cout << "More than 2 b's matched to a Z !!!! Something is wrong" << endl ;
	}
	//now we should have the b's to reconstruct the Z from
	TLorentzVector *genb1 = dynamic_cast<TLorentzVector*>(genP4->At(v_matched_b_to_Z->at(0))) ;
	TLorentzVector *genb2 = dynamic_cast<TLorentzVector*>(genP4->At(v_matched_b_to_Z->at(1))) ;
	TLorentzVector *Z = new TLorentzVector() ;
	*Z = *genb1 + *genb2 ;
	Double_t Z_Pt = Z->Pt() ;

return Z_Pt ;
}

Double_t Get_dR_between_bs( std::vector<Double_t> *genPdgId , std::vector<Double_t> *genMotherPdgId , TClonesArray* genP4)
{
	// loop through and find what index of the b matched to the Z
	std::vector<Int_t> *v_matched_b_to_Z = new std::vector<Int_t> ;
	for (unsigned int i = 0; i < genPdgId->size(); ++i)
	{
		if( (abs(genPdgId->at(i)) == 5 ) && (genMotherPdgId->at(i) == 23))
		{
			v_matched_b_to_Z->push_back(i);
		}
	}
	// now that we have the b's that came from a Z we ca take the dR from them
	if (v_matched_b_to_Z->size() > 2) 
	{
		cout << "More than 2 b's matched to a Z !!!! Something is wrong" << endl ;
	}
	TLorentzVector *genb1 = dynamic_cast<TLorentzVector*>(genP4->At(v_matched_b_to_Z->at(0))) ;
	TLorentzVector *genb2 = dynamic_cast<TLorentzVector*>(genP4->At(v_matched_b_to_Z->at(1))) ;
	Double_t dR = 0 ;
	dR = genb1->DrEtaPhi(*genb2) ;
return dR ;
}


void Signal_Z_bb_jets(std::vector<TTree*> *v_trees)
{
	cout << "Understanding the signal region" << endl ;
	TTree* tree = new TTree();

	// Make histograms to store the signal and backgrounds regions
	// Int_t nbins = 9 ; Double_t xmin = 0 ; Double_t xmax = 10 ;
	TH1D* h_bb_eta = new TH1D("h_bb_eta", "h_bb_eta", 16, -4, 4);
	TH1D* h_jets_eta = new TH1D("h_jets_eta", "h_jets_eta", 16, -4, 4);
	// TH1D* h_bb_delta_eta = new TH1D("h_bb_delta_eta", "h_bb_delta_eta", nbins, xmin, xmax);
	// TH1D* h_jets_delta_eta = new TH1D("h_jets_delta_eta", "h_jets_delta_eta", nbins, xmin, xmax);
	TH1D* h_bb_pt = new TH1D("h_bb_pt", "h_bb_pt", 100, 0, 1000) ;
	TH1D* h_jets_pt = new TH1D("h_jets_pt", "h_jets_pt", 100, 0, 1000) ;

	TH2D* h_bb_eta_phi = new TH2D("h_bb_eta_phi", "h_bb_eta_phi", 40, -180, 180,  16, 4 , 4 ) ;
	TH2D* h_jet_eta_phi = new TH2D("h_jet_eta_phi", "h_jet_eta_phi",  40, -180, 180,  16, 4 , 4) ;

	TH1D* h_dRsquared = new TH1D("h_dRsquared", "h_dRsquared", 20, 0, 1 ) ;

	// loop over the trees
	for (unsigned int i = 0; i < v_trees->size(); ++i)
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
			
			// Decrease the number of entries for testing to 0.05%
			Double_t nentries = tree->GetEntries() ;
			nentries = nentries*0.2 ;
			// determine the min jet pt 
			Double_t lowest_jetpt = 500 ;
			for (int j = 0; j < nentries; ++j)
			{
				tree->GetEntry(j) ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
				if (jetP4->GetEntries() != 1) continue;
				TLorentzVector *jet = new TLorentzVector();
				jet = dynamic_cast<TLorentzVector*>(jetP4->At(0)) ;
				if (jet->Pt() < lowest_jetpt)
				{
					lowest_jetpt = jet->Pt() ;
				} 
			}
			cout << "Lowest jet Pt = " << lowest_jetpt << endl;
			cout << "Looping over " << nentries << " events in dataset " << i << endl ;
			for (int j = 0; j < nentries; ++j)
			{
				tree->GetEntry(j) ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
			
				// Check to see if there are 2 b quarks that are matched to a Z
				if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId , genP4 , 0) == kFALSE) continue ;
				// look for signal with 1 jet
				if (jetP4->GetEntries() != 1) continue;
				//apply cut to b
				std::vector<int>* v_matched_b_to_Z = new std::vector<int> ;
				for (unsigned int k = 0; k < genPdgId->size(); ++k)
				{
					if ((abs(genPdgId->at(k)) == 5) && (genMotherPdgId->at(k) == 23)) // it is a b matched to a Z so plot as signal
					{
						v_matched_b_to_Z->push_back(k) ;
						//get the 4 vector
						TLorentzVector *genb = new TLorentzVector();
						genb = dynamic_cast<TLorentzVector*>(genP4->At(k)) ;
						// add pt cut to the gen b
						if ( genb->Pt() < 30 ) continue ;
						h_bb_pt->Fill(genb->Pt() , mcWeight);
						h_bb_eta->Fill(genb->Eta() , mcWeight);
						h_bb_eta_phi->Fill( genb->Phi() , genb->Eta() , mcWeight );
					}
				}
				// since the jet vector should only be of size 1 we can select that jet
				TLorentzVector *jet = new TLorentzVector();
				jet = dynamic_cast<TLorentzVector*>(jetP4->At(0)) ;
				h_jets_eta->Fill(jet->Eta(), mcWeight) ;
				h_jets_pt->Fill(jet->Pt(), mcWeight) ;
				h_jet_eta_phi->Fill(jet->Phi() , jet->Eta() , mcWeight) ;

				// Get info for dRsquared
				TLorentzVector *genb1 = new TLorentzVector();
				genb1 = dynamic_cast<TLorentzVector*>(genP4->At(v_matched_b_to_Z->at(0))) ;
				TLorentzVector *genb2 = new TLorentzVector();
				genb2 = dynamic_cast<TLorentzVector*>(genP4->At(v_matched_b_to_Z->at(1))) ;
				if (( genb1->Pt() < lowest_jetpt ) || (genb2->Pt() < lowest_jetpt)) continue ;
				Double_t dRsquared = 0 ;
				dRsquared =  pow(genb1->Phi() - genb2->Phi(), 2) + pow(genb1->Eta() - genb2->Eta() , 2) ;
				h_dRsquared->Fill(dRsquared,mcWeight) ;
				v_matched_b_to_Z->clear() ;

			}
		}

		h_bb_eta->SetLineColor(kRed);
		h_jets_eta->SetLineColor(kBlue);
		TLegend* leg1 = new TLegend(0.7, 0.9, 0.90, 0.7) ;
		leg1->AddEntry(h_bb_eta, "b quarks", "l") ;
		leg1->AddEntry(h_jets_eta, "Jets" , "l") ;
		
		h_bb_pt->SetLineColor(kRed);
		h_jets_pt->SetLineColor(kBlue);
		TLegend* leg2 = new TLegend(0.7, 0.9, 0.90, 0.7) ;
		leg2->AddEntry(h_bb_pt, "b quarks", "l") ;
		leg2->AddEntry(h_jets_pt, "Jet" , "l") ;

		h_bb_eta_phi->SetLineColor(kRed) ;
		h_jet_eta_phi->SetLineColor(kBlue) ;
		TLegend* leg3 = new TLegend(0.7, 0.9, 0.90, 0.7) ;
		leg3->AddEntry(h_bb_eta_phi, "b quarks", "l") ;
		leg3->AddEntry(h_jet_eta_phi, "Jet" , "l") ;

		TCanvas* c = new TCanvas() ;
		c->Divide(2,3) ;
		
		c->cd(1);
		h_bb_eta->Draw("hist");
		h_jets_eta->Draw("same hist");
		leg1->Draw() ;
		
		c->cd(2) ;
		h_bb_pt->Draw("hist");
		h_jets_pt->Draw("same, hist");
		leg2->Draw() ;

		c->cd(3) ;
		h_bb_eta_phi->Draw("PSR SURF1") ;
		leg3->Draw() ;

		c->cd(4) ;
		h_dRsquared->Draw("hist");

		c->cd(5) ;


		c->cd(6) ;
		h_jet_eta_phi->Draw("CONTZ") ;

		c->SaveAs("Signal_bb_jets.png") ;

		TCanvas* c2 = new TCanvas() ;
		h_dRsquared->Draw("hist");
		c2->SaveAs("dRsquared.png") ;

}

Bool_t Gen_b_from_Z_cut( std::vector<Double_t> *genPdgId ,std::vector<Double_t> *genMotherPdgId , TClonesArray *genP4 ) 
{
	Bool_t passed_cut = kFALSE;
	Int_t count = 0 ;
	// make sure the b is matched to a Z
	for (unsigned int i = 0; i < genPdgId->size(); ++i)
	{
		if (( abs(genPdgId->at(i)) == 5 ) && (genMotherPdgId->at(i) == 23 ))
		{
			TLorentzVector* genb = dynamic_cast<TLorentzVector*> (genP4->At(i)) ;
			// place the pt cut on the b's
			if ( genb->Pt() < 30 ) continue ;
			// Place the Eta cut on the b's
			if (abs(genb->Eta()) < 2.4) continue ;
			count = count + 1 ;
		}
	}
	if (count == 2) 
	{
		passed_cut = kTRUE ;
	}
	else
	{
		passed_cut = kFALSE ;
	}
return passed_cut ;
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
	for (unsigned int i = 0; i < v_trees->size(); ++i)
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
			nentries = nentries*0.01 ;
			cout << "Looping over " << nentries << " events in dataset " << i << endl ;
			for (int j = 0; j < nentries; ++j)
			{
				Double_t x = 0 ;
				tree->GetEntry(j) ;
				if      ( variable == "njets") 
				{	
					Int_t njets = jetP4->GetEntries() ;
					// loop over jets and determine how many of them are within the eta range
					for (int i = 0; i < njets; ++i)
					{
						TLorentzVector* jet = dynamic_cast<TLorentzVector*> ( jetP4->At(i) ) ;
						if ( abs(jet->Eta()) < 2.4 )
						{
							x = x+1 ;
						}
					}
				}
				else if ( variable == "met")   { TLorentzVector* MetP4 = dynamic_cast<TLorentzVector*>(metP4->At(0)) ;    x = MetP4->Pt() ;}  
				else {cout << " !!! Incorrect choice of variable !!!" << endl;}
				// cout << variable << " = " << x << endl ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}

				// Check to see if there are 2 b quarks that are matched to a Z
				if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId , genP4 , 0) == kTRUE)
				{
					//fill the signal histogram
					// h_sig->Fill( Get_dijet_mass(jetP4) , mcWeight ) ;
					if ( Gen_b_from_Z_cut( genPdgId , genMotherPdgId , genP4 ) == kFALSE ) continue ;
					h_sig->Fill( x , mcWeight ) ;
				}
				else 
				{
					// fill the background histogram
					// h_bkg->Fill( Get_dijet_mass(jetP4) , mcWeight ) ;
					if ( Gen_b_from_Z_cut( genPdgId , genMotherPdgId , genP4 ) == kFALSE ) continue ;
					h_bkg->Fill( x , mcWeight ) ;
				}
			}
		}
		// Get the integrals of the signal and background
		Double_t sig_integral = 0 ;  
		Double_t bkg_integral = 0 ; 
		sig_integral = h_sig->Integral(); cout << "Signal Integral = " << sig_integral << " Signal Entries = " << h_sig->GetEntries() << endl ;
		bkg_integral = h_bkg->Integral(); cout << "Background Integral = " << bkg_integral <<  " Background Entries = " << h_bkg->GetEntries() <<  endl ;
		THStack* s_sig_and_bkg_njets =  Plot_Signal_and_Background( h_sig , h_bkg, variable ) ;
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
			for (unsigned int i = 0; i < v_trees->size(); ++i)
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
					nentries = nentries*0.1;
					for (int j = 0; j < nentries; ++j)
					{
						tree->GetEntry(j) ;
						if (j%10000 == 0)
						{
							cout << "10000 events" << endl;
						}
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
						if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId , genP4 , 0 ) == kTRUE) 
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





// Use dR to determine if the jets are matching to the b's that came from the Z
Int_t Get_n_jets_matched_to_gen_b_from_Z( std::vector<Double_t>* genPdgId , std::vector<Double_t>* genMotherPdgId , TClonesArray* genP4 , TClonesArray* jetP4 )
{
	Int_t n_jets_matched_to_b_from_Z = 0 ;
	// Start with a jet
	std::vector<Int_t> *b_index = new std::vector<Int_t> ;
	// start looping through the jets
	for (int i = 0; i < jetP4->GetEntries(); ++i)
	{
		TLorentzVector* jet = dynamic_cast<TLorentzVector*> ( jetP4->At(i) ) ;
		// now loop through the b's from Z and determine which is within dR
		for (unsigned int j = 0; j < genPdgId->size(); ++j)
		{
			if ((abs(genPdgId->at(j)) == 5) && (genMotherPdgId->at(j) == 23))
			{
				TLorentzVector* genb = dynamic_cast<TLorentzVector*> (genP4->At(j)) ;
				// calclate the dR squared with this jet and a b quark
				Double_t dRsquared = 0 ;
				dRsquared = pow( genb->Eta() - jet->Eta()  ,2) + pow( genb->Phi() - jet->Phi()  , 2) ;
				// check if dRsquared > 0.3^2
				if (dRsquared > pow(0.3,2) ) continue ;
				// Make sure it is not claimed already
				if (b_index->size() == 0) // then it is the first b so we are safe to match it
				{
					b_index->push_back(j) ;
				}
				else
				{
					// loop thorugh the vector and make sure the b quark was not cliamed already by another jet
					for (unsigned int k = 0; k < b_index->size(); ++k)
					{
						if (  j != b_index->at(k) )
						{
							b_index->push_back(j) ;
						}
					}
				}
			}
		}
	}	

	n_jets_matched_to_b_from_Z = b_index->size() ;
	b_index->clear() ;

return n_jets_matched_to_b_from_Z;
}



void N_jets_and_Zbb(std::vector<TTree*> *v_trees )
{
	cout << "Looking at how njets varies with Zpt" << endl ;
	TTree* tree = new TTree();

	// Make histograms to store the signal and backgrounds regions
	// Int_t nbins = 9 ; Double_t xmin = 0 ; Double_t xmax = 10 ;
	TH2D* h_nbjets_Z_pt = new TH2D( "h_nbjets_Z_pt" , "h_nbjets_Z_pt" , 20, 0, 1000, 5, 0, 5 );

	// loop over the trees
	for (unsigned int i = 0; i < v_trees->size(); ++i)
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
			nentries = nentries ;
			cout << "Looping over " << nentries << " events in dataset " << i << endl ;
			for (int j = 0; j < nentries; ++j)
			{
				Double_t x = 0 ;
				tree->GetEntry(j) ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
				// Variable to keep track of jets that are matched to b's that come from Z's
				Int_t n_jets_matched_to_b_from_Z = 0 ;
				// make sure the event is signal i.e Z->bb
				if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId , genP4 , 0) == kFALSE) continue ;
				// Figure out the Z_pt 
				Double_t Z_pt = Get_Z_pT_gen_level(genPdgId , genMotherPdgId ,  genP4 ) ;
				// now figure out how many jets were matched to a b 
				// since this is signal we should ideally have 2 but we are seeing cases where we get 1
				n_jets_matched_to_b_from_Z = Get_n_jets_matched_to_gen_b_from_Z( genPdgId , genMotherPdgId , genP4 , jetP4 ) ;
				// Fill my 2D histogram
				h_nbjets_Z_pt->Fill( Z_pt , n_jets_matched_to_b_from_Z  , mcWeight) ;
			}
		}
	h_nbjets_Z_pt->SetStats(0) ;
	TCanvas* c3 = new TCanvas() ;
	c3->Divide(2,1) ;
	c3->cd(1);
	h_nbjets_Z_pt->Draw("LEGO") ;
	c3->cd(2);
	h_nbjets_Z_pt->Draw("COLZ") ;
	c3->SaveAs("h_nbjets_Z_pt.png") ;

}

Int_t Get_n_b_in_fidacc( std::vector<Double_t>* genPdgId , std::vector<Double_t>* genMotherPdgId , TClonesArray* genP4 )
{
	Int_t n_b_from_Z_in_fidacc = 0 ;
	// start looping through the b's matched to the Z's
	for (unsigned int i = 0; i < genPdgId->size() ; ++i)
	{
		if (!((abs(genPdgId->at(i))==5) && (genMotherPdgId->at(i) == 23))) continue ;
		TLorentzVector* genb = dynamic_cast<TLorentzVector*> (genP4->At(i)) ;
		cout << "Gen b from Z with eta = " << genb->Eta() << "; " ;
		if (abs(genb->Eta())<2.5) 
		{
			n_b_from_Z_in_fidacc = n_b_from_Z_in_fidacc + 1 ;
		}
		cout << "  " << endl ;
	}	

return n_b_from_Z_in_fidacc;
}

std::vector<Double_t>* Get_n_jets_matched_to_gen_b_from_Z_in_fidacc( std::vector<Double_t>* genPdgId , std::vector<Double_t>* genMotherPdgId , TClonesArray* genP4 , TClonesArray* jetP4 )
{
	std::vector<Double_t>* v_njets_btag_value = new std::vector<Double_t> ;
	Double_t n_jets_matched_to_b_from_Z_in_fidacc = 0 ;
	// Start with a b in the accetped region
	std::vector<Int_t> *b_index = new std::vector<Int_t> ;
	// start looping through the b to find the ones that came from a Z
	for (unsigned int i = 0; i < genPdgId->size(); ++i)
	{
		// make sure it came from a Z (we already looking at signal events)
		if (!((abs(genPdgId->at(i)) == 5) && (genMotherPdgId->at(i) == 23))) continue ;
		// Now I have a b that came from a Z
		TLorentzVector* genb = dynamic_cast<TLorentzVector*> ( genP4->At(i) ) ;
		// Apply the fiducial accpetance to the b
		if ( abs( genb->Eta() ) > 2.5 ) continue ;
		// now loop through the jets from the Z->b and determine which is within dR < 0.3
		for (int j = 0; j < jetP4->GetEntries(); ++j)
		{
			TLorentzVector* jet = dynamic_cast<TLorentzVector*> (jetP4->At(j)) ;
			// calclate the dR squared with this jet and a b quark
			Double_t dR = 0 ;
			dR = genb->DrEtaPhi( *jet );
			// check if dRsquared > 0.3^2
			if (dR > 0.5 ) continue ;
			// cout << "Jet matched to b with dRsquared = " << dRsquared << endl ;
			// Make sure it is not claimed already
			if (b_index->size() == 0) // then it is the first jet so we are safe to match it
			{
				b_index->push_back(j) ;
			}
			else
			{
				// loop thorugh the vector and make sure the jet was not cliamed already by another b quark 
				for (unsigned int k = 0; k < b_index->size(); ++k)
				{
					if (  j != b_index->at(k) )
					{
						b_index->push_back(j) ;
					}
				}
			}
		}
	}	

	n_jets_matched_to_b_from_Z_in_fidacc = b_index->size() ;
	v_njets_btag_value->push_back(n_jets_matched_to_b_from_Z_in_fidacc) ;
	b_index->clear() ;
	// cout << "N jets matched to b quarks from a Z = " << n_jets_matched_to_b_from_Z_in_fidacc << endl ;
return v_njets_btag_value;
}


std::vector<Double_t>* Get_n_btagged_jets_matched_to_Z_in_fidacc( 	std::vector<Double_t>* genPdgId , 
																	std::vector<Double_t>* genMotherPdgId , 
																	TClonesArray* genP4 , TClonesArray* jetP4 , 
																	std::vector<float> *btagger , 
																	TString working_point, 
																	TH1D* h_tagger_value )
{
	std::vector<Double_t> *v_njets_btag_value = new std::vector<Double_t> ;
	Double_t n_jets_matched_to_b_from_Z_in_fidacc = 0 ;
	Double_t Working_Point = 0 ;
	if (working_point == "No_Tag")
	{
		Working_Point = 0 ;
	}
	else if(working_point == "Loose" )
	{
		Working_Point = 0.432 ;
	}
	else if (working_point == "Medium")
	{
		Working_Point = 0.841 ;
	}	
	else if (working_point == "Tight")
	{
		Working_Point = 0.941 ;
	}
	else 
	{
		Working_Point = 0 ;
		cout << " !!! Error in working points !!!" << endl;
	}
	// Start with a b in the accetped region
	std::vector<Int_t> *b_index = new std::vector<Int_t> ;
	// start looping through the b
	cout << "Start looping through the b quarks that came from a Z to match to jet "  << endl ;
	for (unsigned int i = 0; i < genPdgId->size(); ++i)
	{
		// make sure it came from a Z (we already looking at signal events)
		if (!((abs(genPdgId->at(i)) == 5) && (genMotherPdgId->at(i) == 23))) continue ;
		TLorentzVector* genb = dynamic_cast<TLorentzVector*> ( genP4->At(i) ) ;
		// Apply the fiducial accpetance to the b
		if ( abs( genb->Eta() ) > 2.5 ) continue ;
		cout << " B quark from Z identified at index " << i << " with Eta = " << genb->Eta() << endl ;
		// now loop through the jets from Z and determine which is within dR
		for (int j = 0; j < jetP4->GetEntries(); ++j)
		{
			cout << "Looking for jet within dR < 0.5 to this b quark" << endl ;
			TLorentzVector* jet = dynamic_cast<TLorentzVector*> (jetP4->At(j)) ;
			if ( Working_Point != 0 )
			{
				if (btagger->at(j) < Working_Point) continue;
			}
			// calclate the dR squared with this jet and a b quark
			Double_t dR = 0 ;
			dR = genb->DrEtaPhi( *jet );
			cout << "Jet at index " << j << " has dR = " << dR << endl ;
			// check if dRsquared > 0.3^2
			if (dR > 0.5 ) continue ;
			cout << "Jet was within range" << endl ;
			// Make sure it is not claimed already
			if (b_index->size() == 0) // then it is the first jet so we are safe to match it
			{
				b_index->push_back(j) ;
				h_tagger_value->Fill(btagger->at(j)) ;
			}
			else
			{
				// loop thorugh the vector and make sure the jet was not cliamed already by another b quark 
				Double_t not_claimed = 0 ;
				for (unsigned int k = 0; k < b_index->size(); ++k)
				{
					if ( j != b_index->at(k) ) 
					{
						not_claimed = not_claimed + 1 ;
					}
				}
				if ( not_claimed == b_index->size() ) // then all the indexes are taken already so it is safe to add the jet as matched to the b
				{
					b_index->push_back(j) ;
					h_tagger_value->Fill(btagger->at(j)) ;
				}
			}
		}
	}	

	n_jets_matched_to_b_from_Z_in_fidacc = Double_t(b_index->size()) ;
	v_njets_btag_value->push_back(n_jets_matched_to_b_from_Z_in_fidacc) ;
	b_index->clear() ;
	// cout << "N jets matched to b quarks from a Z = " << n_jets_matched_to_b_from_Z_in_fidacc << endl ;
return v_njets_btag_value;
}


void Plot_n_jets_in_Signal(std::vector<TTree*> *v_trees)
{
	cout << "Looking at number of jets in the signal region"  <<endl ;
	TTree* tree = new TTree();

	// X will be the Z_pt and y will be njets
	TH2D* h_2D_njets_Z_pt = new TH2D("h_2D_njets_Z_pt", "Njets as a function of Z Pt" , 20, 0, 1000, 10, 0, 7);
	// loop over the trees
	for (unsigned int i = 0; i < v_trees->size(); ++i)
		{
			tree = v_trees->at(i) ;	

			std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
			std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
			std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
			std::vector<float> *btagger = new std::vector<float> ;	
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
			tree->SetBranchAddress("CombinedInclusiveSecondaryVertexV2", &btagger) ;

			Double_t nentries = tree->GetEntries() ;
			// Decrease the number of entries for testing to 0.05%
			nentries = nentries  ;
			// cout << "Looping over " << nentries << " events in dataset " << i << endl ;
			for (int j = 0; j < nentries; ++j)
			{
				Double_t x = 0 ;
				tree->GetEntry(j) ;
				if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
				// make sure the event is signal i.e Z->bb
				if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId , genP4 , 0) == kFALSE) continue ;
				// Figure out the Z_pt 
				Double_t Z_pt = Get_Z_pT_gen_level(genPdgId , genMotherPdgId ,  genP4 ) ;
				// determine what bin this corresponds to in my Z_pt binning
				// Int_t bin_no = int(Z_pt / 50) + 1 ;
				Int_t nJets = 0 ;
				nJets = jetP4->GetEntries() ;
				// Fill my 2D histogram with the Z_pt and the njets distributions
				h_2D_njets_Z_pt->Fill( Z_pt , nJets , mcWeight) ;
			}
		}
		h_2D_njets_Z_pt->GetXaxis()->SetTitle("Z Pt") ;
		h_2D_njets_Z_pt->GetYaxis()->SetTitle("# of Jets") ;
		h_2D_njets_Z_pt->SetStats(0) ;
		TCanvas* c = new TCanvas("c", "N jets as a function of Z Pt", 1000, 750) ;
		c->SetLogy();
		c->Divide(2,1) ;
		c->cd(1) ;
		h_2D_njets_Z_pt->Draw("LEGO") ;
		c->cd(2) ;
		h_2D_njets_Z_pt->Draw("COLZ  TEXT") ;
		c->SaveAs("Njets_Z_pt.png") ;
		cout << "Finished the njets as a function of Z-pt" << endl ;
}

std::vector< std::vector<Double_t>* >* Get_v_dR__ij__b_jet(std::vector<Double_t>* genPdgId ,std::vector<Double_t>* genMotherPdgId, TClonesArray*  genP4 , TClonesArray* jetP4 ) 
{
	// make matrix to store all the dR's
	std::vector< std::vector<Double_t>* > *v_dR__ij__b_jet = new std::vector< std::vector<Double_t>* > ;
	// Loop over the b-quarks
	for (unsigned int i = 0; i < genPdgId->size() ; ++i)
	{
		// Get a b-quark that came from a Z
		if (!( abs(genPdgId->at(i)) == 5 && genMotherPdgId->at(i) == 23 )) continue ;
		TLorentzVector* genb = dynamic_cast<TLorentzVector*> (genP4->At(i)) ;
		// now loop over the jets
		std::vector<Double_t>* v_dR_b_jet = new std::vector<Double_t> ;
		for (unsigned int j = 0; j < jetP4->GetEntries(); ++j)
		{
			TLorentzVector* jet = dynamic_cast<TLorentzVector*> ( jetP4->At(j) ) ;
			// get dR with the b and this jet
			Double_t dR = 0 ;
			dR = genb->DrEtaPhi( *jet ) ;
			v_dR_b_jet->push_back(dR) ;
		}
		// now push this vector into the matrix
		v_dR__ij__b_jet->push_back( v_dR_b_jet ) ;
	}
return v_dR__ij__b_jet ;
}

void Output_Matrix( std::vector<std::vector<Double_t>* > *v_dR__ij__b_jet)
{
	for (unsigned int i = 0; i < v_dR__ij__b_jet->size(); ++i)
	{
		// Get a single vector and loop over, 
		for (unsigned int j = 0; j < v_dR__ij__b_jet->at(i)->size(); ++j)
		{
			cout << " " << v_dR__ij__b_jet->at(i)->at(j) << " " ;
		}
	cout << " " << endl ;
	}
}

std::vector<Double_t>* Get_v_dR1_i1_dR2_i2( std::vector<std::vector<Double_t>* > *v_dR__ij__b_jet)
{
	// Loop through the matrix and get best dR and row value of best dR
	Double_t best_dR = 10000 ; // set to very large value 
	Double_t best_dR_i = 0 ;
	Double_t best_dR_j = 0 ;
	// Find the best dR
	for (unsigned int i = 0; i < v_dR__ij__b_jet->size(); ++i)
	{
		for (unsigned int j= 0; j < v_dR__ij__b_jet->at(i)->size(); ++j)
		{
			if ( v_dR__ij__b_jet->at(i)->at(j) < best_dR ) 
			{
				best_dR = v_dR__ij__b_jet->at(i)->at(j) ;
				best_dR_i = i ;
				best_dR_j = j ;
			}
		}
	}
	cout << "Best dR = " << best_dR << " at i = " << best_dR_i <<  " and j = " << best_dR_j << endl ;
	// Now find the second best dR
	Double_t second_best_dR = 10000;
	Double_t second_best_dR_i = 0 ;
	Double_t second_best_dR_j = 0 ;
	for (unsigned int i = 0; i < v_dR__ij__b_jet->size(); ++i)
	{
		// cout << " i = " << i << endl ;
		for (unsigned int j = 0; j < v_dR__ij__b_jet->at(i)->size(); ++j)
		{
			// cout << "j = " << j << endl ;
			if ( (v_dR__ij__b_jet->at(i)->at(j) < second_best_dR) && ( j != best_dR_j) )
			{
				second_best_dR = v_dR__ij__b_jet->at(i)->at(j) ;
				second_best_dR_i = i ;
				second_best_dR_j = j ;
			}
		}
	}
	cout << "Second best dR = " << second_best_dR <<  " at i = " << second_best_dR_i << " and j = " << second_best_dR_j << endl ;
	// Now put the values we want in the vector
	std::vector<Double_t> *v_dR1_i1_dR2_i2 = new std::vector<Double_t> ;
	v_dR1_i1_dR2_i2->push_back(best_dR) ; v_dR1_i1_dR2_i2->push_back(best_dR_j ) ;
	v_dR1_i1_dR2_i2->push_back(second_best_dR) ; v_dR1_i1_dR2_i2->push_back(second_best_dR_j ) ;
return v_dR1_i1_dR2_i2 ;
}

void Get_Efficiency_func(std::vector<TTree*> *v_trees , TString Jet_Option, Bool_t bMatched)
{
	cout << "Looking at Efficiency with intermediate bins and jet option " << Jet_Option <<endl ;
	TTree* tree = new TTree();
	TString s1 ;
	if (bMatched == kTRUE)
	{
		s1 = "b_matching" ;
	}
	else
	{
		s1 = "no_b_matching" ;
	}

	Double_t Z_pt_bins[10] = {0 , 50, 100, 150, 200, 250, 300, 400, 500, 1000} ;

	cout << "Creating necessary histograms" << endl;
	// Make hists for n accepted as a function of Z_pt
	TH1D* h_accepted_funtion_Z_pt = new TH1D("h_accepted_funtion_Z_pt" , " Accepted Z#rightarrow b#bar{b} at Gen level " + s1+ " " +Jet_Option , 9 , Z_pt_bins);
	TH1D* h_accepted_funtion_Z_pt_pos = new TH1D("h_accepted_funtion_Z_pt_pos" , " Accepted Z#rightarrow b#bar{b} at Gen level with Pos weight " + s1+ " " +Jet_Option, 9 , Z_pt_bins);
	TH1D* h_accepted_funtion_Z_pt_neg = new TH1D("h_accepted_funtion_Z_pt_neg" , " Accepted Z#rightarrow b#bar{b} at Gen level with Neg weight "+ s1+ " "+ Jet_Option , 9 , Z_pt_bins);
	h_accepted_funtion_Z_pt->SetLineColor(kBlue) ;
	h_accepted_funtion_Z_pt_neg->SetLineColor(kBlue);
	h_accepted_funtion_Z_pt_pos->SetLineColor(kBlue) ;
	// Make hists for n signal events where the jets were matched
	TH1D* h_bmatched_function_Z_pt = new TH1D("h_bmatched_function_Z_pt" , " Accepted Z#rightarrow 2 b Jets "+ s1+ Jet_Option , 9 , Z_pt_bins) ;
	TH1D* h_bmatched_function_Z_pt_pos = new TH1D("h_bmatched_function_Z_pt_pos" , " Accepted Z#rightarrow 2 b Jets with Pos weight " + s1+ " "+ Jet_Option , 9 , Z_pt_bins) ;
	TH1D* h_bmatched_function_Z_pt_neg = new TH1D("h_bmatched_function_Z_pt_neg" , " Accepted Z#rightarrow 2 b Jets with Neg Weight "+ s1+ " " + Jet_Option, 9 , Z_pt_bins) ;
	h_bmatched_function_Z_pt->SetLineColor(kRed) ;
	h_bmatched_function_Z_pt_pos->SetLineColor(kRed) ;
	h_bmatched_function_Z_pt_neg->SetLineColor(kRed) ;

	TH1D* h_accepted_function_dR =  new TH1D( "h_accepted_function_dR" , "Accepted Z#rightarrow b#bar{b} at Gen Level "+ s1+ " "+ Jet_Option , 10, 0, 1 );
	TH1D* h_accepted_function_dR_pos = new TH1D("h_accepted_function_dR_pos" , " Accepted Z#rightarrow b#bar{b} at Gen Level with postive weight "+ s1+ " "+ Jet_Option , 10, 0, 1);
	TH1D* h_accepted_function_dR_neg = new TH1D(" h_accepted_funtion_Z_pt_neg" , " Accepted Z#rightarrow b#bar{b} at Gen Level with negative weight " + s1+ " "+ Jet_Option, 10, 0, 1) ;
	h_accepted_function_dR->SetLineColor(kBlue) ;
	h_accepted_function_dR_pos->SetLineColor(kBlue) ;
	h_accepted_function_dR_neg->SetLineColor(kBlue) ;


	TH1D* h_bmatched_function_dR = new TH1D("h_bmatched_function_dR" , "Accepted Z#rightarrow 2 b Jets "+ s1+ " " + Jet_Option, 10, 0, 1);
	TH1D* h_bmatched_function_dR_pos = new TH1D("h_bmatched_function_dR_pos" , "Accepted Z#rightarrow 2 b Jets with Pos weight "+ s1+ " "+ Jet_Option , 10, 0, 1);
	TH1D* h_bmatched_function_dR_neg = new TH1D("h_bmatched_function_dR_neg" , "Accepted Z#rightarrow 2 b Jets with Neg weight "+ s1+ " " + Jet_Option, 10, 0, 1);

	h_bmatched_function_dR->SetLineColor(kRed) ;
	h_bmatched_function_dR_pos->SetLineColor(kRed) ;
	h_bmatched_function_dR_neg->SetLineColor(kRed) ;

	TH1D* h_best_btagger_value = new TH1D("h_best_btagger_value" , "Best dR " , 10, 0, 1) ;
	h_best_btagger_value->SetFillColor(kPink) ;
	TH1D* h_second_best_btagger_value = new TH1D("h_second_best_btagger_value" , "Second best dR" , 10, 0, 1) ;
	h_second_best_btagger_value->SetFillColor(kGreen) ;

	THStack* s_btagger = new THStack("s_btagger" , " B-tagger " ) ;

	TH1D* h_tagger_value = new TH1D("h_tagger_value" , " Value of b-tagger on b-matched Jets" + s1 +"_" + Jet_Option  , 20 , 0, 1) ;
	cout << "Begin looping over Trees" << endl ;
	// loop over the trees
	for (unsigned int i = 0; i < v_trees->size(); ++i)
	{
		tree = v_trees->at(i) ;	

		std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
		std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
		std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
		std::vector<float> *btagger = new std::vector<float> ;	
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
		tree->SetBranchAddress("CombinedInclusiveSecondaryVertexV2", &btagger) ;

		Double_t nentries = tree->GetEntries() ;
		// Decrease the number of entries for testing to 0.05%
		nentries = nentries*0.01;
		// cout << "Looping over " << nentries << " events in dataset " << i << endl ;
		for ( unsigned int j = 0; j < 50; ++j)
		{
			cout << "Looking at event = " << j << endl ;
			// std::vector<Double_t> *v_njets_btag_value = new std::vector<Double_t> ;
			tree->GetEntry(j) ;
			if (mcWeight > 0) {mcWeight = 1;}else{mcWeight = -1;}
			
			// make sure the event is signal i.e Z->bb with pt cut on b quarks
			if (Two_Gen_b_Quarks_Matched_to_Z( genPdgId , genMotherPdgId , genP4 , 0 ) == kFALSE)
			{
				cout << "Not signal event"<< endl ;
				continue ;
			}
			cout << "Event is SIGNAL" << endl ;
			// Now figure out how many b quarks from Z are within the detector acceptance
			Double_t n_b_in_fidacc = 0 ;
			n_b_in_fidacc = Get_n_b_in_fidacc(genPdgId , genMotherPdgId ,  genP4 );
			if (n_b_in_fidacc < 2) 
			{
				cout << "At least 1 b-quark was outside of fiducial acceptance. Go to next event " << endl ;
				continue ; // only signal events pass	
			}
			cout << "Both B-quarks are in fiducial acceptance" << endl ;
			// Figure out the Z_pt for histogram
			Double_t Z_pt = Get_Z_pT_gen_level(genPdgId , genMotherPdgId ,  genP4 ) ;
			// Figure out the dR between the b quarks for histogram
			Double_t dR_bquarks = Get_dR_between_bs(genPdgId , genMotherPdgId ,  genP4) ; 
			
			// fill the accepted histogram
			h_accepted_funtion_Z_pt->Fill( Z_pt , mcWeight ) ;
			h_accepted_function_dR->Fill(dR_bquarks , mcWeight ) ;
			//---------------------------------------------
			if (mcWeight == 1) 
			{
				h_accepted_funtion_Z_pt_pos->Fill(Z_pt , 1) ;
				h_accepted_function_dR_pos->Fill(dR_bquarks, 1) ;
			} 
			else 
			{ 
				h_accepted_funtion_Z_pt_neg->Fill(Z_pt, 1); 
				h_accepted_function_dR_neg->Fill(dR_bquarks, 1) ;
			}
			//--------------------------------------------

			// now figure out how many jets were matched to a b if the b is already in the accepted region
			// Double_t n_jets_matched_to_b_from_Z_in_fidacc = 0 ;
			// if (bMatched == kTRUE)
			// {
			// 	std::vector<Double_t> *v_njets_btag_value = Get_n_btagged_jets_matched_to_Z_in_fidacc( genPdgId ,genMotherPdgId ,  genP4 ,jetP4 , btagger  ,  Jet_Option , h_tagger_value) ;
			// 	n_jets_matched_to_b_from_Z_in_fidacc = v_njets_btag_value->at(0) ;
			// 	// Fill_btagger_hist()
			// }
			// else
			// {
			// 	n_jets_matched_to_b_from_Z_in_fidacc = jetP4->GetEntries() ;
			// }

			// // now deal with the numerator
			// if (n_jets_matched_to_b_from_Z_in_fidacc >= 2)
			// {
			// 	h_bmatched_function_Z_pt->Fill(Z_pt , mcWeight) ;
			// 	h_bmatched_function_dR->Fill(dR_bquarks, mcWeight) ;
			// 	if (mcWeight == 1) 
			// 	{
			// 		h_bmatched_function_Z_pt_pos->Fill(Z_pt , 1) ;
			// 		h_bmatched_function_dR_pos->Fill(dR_bquarks, 1) ;
			// 	} 
			// 	else 
			// 	{
			// 		h_bmatched_function_Z_pt_neg->Fill(Z_pt , 1) ;
			// 		h_bmatched_function_dR_neg->Fill(dR_bquarks, 1) ;
			// 	}
			// }

			// Try tighter requirement on the jet matching
			// Find every combination of dR with a signal b and jet
			// then select best two dR's 
			if (bMatched == kTRUE)
			{
				cout << "Option selected to use b-matching" << endl ;
				cout << "Number of jets in this event = " << jetP4->GetEntries() << endl ;
				// Produce matrix of dR's
				cout << "Making dR matrix" << endl;
				std::vector< std::vector<Double_t>* >* v_dR__ij__b_jet = Get_v_dR__ij__b_jet( genPdgId , genMotherPdgId,  genP4 , jetP4 ) ;
				// now make an output for self checking that looks like the matrix format
				cout << "Outputing Matrix" << endl ;
				Output_Matrix( v_dR__ij__b_jet ) ;
				// Now get the best 2 dR's in the matrix and the indexs of the jet
				cout << "Getting best 2 dRs" << endl ;
				std::vector<Double_t>* v_dR1_i1_dR2_i2 = Get_v_dR1_i1_dR2_i2(v_dR__ij__b_jet) ;
				// The requirement is that dR < 0.5. 
				cout << "Checking if they are in the required dR < 0.5" << endl ;
				// Check if they are both in the region
				if (( v_dR1_i1_dR2_i2->at(0) < 0.5 ) &&  (v_dR1_i1_dR2_i2->at(2) < 0.5 ))
				{
					cout << "Best 2 dR's were in the required dR < 0.5" << endl ;
					h_bmatched_function_Z_pt->Fill(Z_pt , mcWeight) ;
					h_bmatched_function_dR->Fill(dR_bquarks, mcWeight) ;
					if (mcWeight == 1) 
					{
						h_bmatched_function_Z_pt_pos->Fill(Z_pt , 1) ;
						h_bmatched_function_dR_pos->Fill(dR_bquarks, 1) ;
					} 
					else 
					{
						h_bmatched_function_Z_pt_neg->Fill(Z_pt , 1) ;
						h_bmatched_function_dR_neg->Fill(dR_bquarks, 1) ;
					}
					cout << "Filling the btagger values" << endl ;
					// Fill the btagger values to see how the distribution looks like
					cout << "Best dR has btagger value        = "  << btagger->at(v_dR1_i1_dR2_i2->at(1)) << endl ;
					cout << "Second best dR has btagger value = " << btagger->at(v_dR1_i1_dR2_i2->at(3)) << endl ;
					h_best_btagger_value->Fill( btagger->at(v_dR1_i1_dR2_i2->at(1)) , mcWeight ) ;
					h_second_best_btagger_value->Fill(btagger->at(v_dR1_i1_dR2_i2->at(3)) , mcWeight) ;
				}
				else
				{
					cout << "Selected dR's were not within the required dR" << endl ;
					continue ;
				}
			}
			else
			{
				cout << "Option selected to not use b-matching" << endl ;
				Double_t n_jets_matched_to_b_from_Z_in_fidacc = 0 ;
				n_jets_matched_to_b_from_Z_in_fidacc = jetP4->GetEntries() ;
				// now deal with the numerator
				if (n_jets_matched_to_b_from_Z_in_fidacc >= 2)
				{
					h_bmatched_function_Z_pt->Fill(Z_pt , mcWeight) ;
					h_bmatched_function_dR->Fill(dR_bquarks, mcWeight) ;
					if (mcWeight == 1) 
					{
						h_bmatched_function_Z_pt_pos->Fill(Z_pt , 1) ;
						h_bmatched_function_dR_pos->Fill(dR_bquarks, 1) ;
					} 
					else 
					{
						h_bmatched_function_Z_pt_neg->Fill(Z_pt , 1) ;
						h_bmatched_function_dR_neg->Fill(dR_bquarks, 1) ;
					}
				}
			}
		}
	}

		// TCanvas* c = new TCanvas("c_Efficiency_Zpt" + s1 +"_" + Jet_Option, "c_Efficiency_Zpt" + s1 +"_" + Jet_Option , 1000, 2000) ;
		// c->Divide(2,3) ;
		
		// c->cd(1) ;
		// h_accepted_funtion_Z_pt->Draw("hist") ;
		// h_bmatched_function_Z_pt->Draw("same hist") ;
		
		// c->cd(2) ;

		// TGraphAsymmErrors* g_Efficiency = new TGraphAsymmErrors( h_bmatched_function_Z_pt , h_accepted_funtion_Z_pt , "n"  );
		// g_Efficiency->Draw() ;
		
		// c->cd(3);
		// h_accepted_funtion_Z_pt_pos->Draw("hist") ;
		// h_bmatched_function_Z_pt_neg->Draw("same hist") ;
		
		// c->cd(4) ;
		// TGraphAsymmErrors* g_Efficiency_pos = new TGraphAsymmErrors( h_bmatched_function_Z_pt_pos , h_accepted_funtion_Z_pt_pos , "n" ) ;
		// g_Efficiency_pos->Draw() ;
	
		// c->cd(5) ;
		// h_accepted_funtion_Z_pt_neg->Draw() ;
		// h_bmatched_function_Z_pt_neg->Draw("same") ;
		
		// c->cd(6) ;
		// TGraphAsymmErrors* g_Efficiency_neg = new TGraphAsymmErrors(h_bmatched_function_Z_pt_neg , h_accepted_funtion_Z_pt_neg, "n") ;
		// g_Efficiency_neg->Draw() ;

		// c->SaveAs("Efficiency_Zpt"+ s1 +"_" + Jet_Option+".png") ;
		// c->Clear() ;
		// // TString name = "Efficiency" + Jet_Option + "bmatched.png" ;
		// // c->SaveAs(name) ;

		// TCanvas *c_Eff_dR_b_quarks = new TCanvas("c_Eff_dR_b_quarks" + s1 +"_" + Jet_Option, "c_Efficiency_Eff_dR_b_quarks" + s1 +"_" + Jet_Option, 1000, 2000) ;
		// c_Eff_dR_b_quarks->Divide(2,3) ;

		// c_Eff_dR_b_quarks->cd(1);
		// h_accepted_function_dR->Draw("hist") ;
		// h_bmatched_function_dR->Draw("hist same") ;
		
		// c_Eff_dR_b_quarks->cd(2);
		// TGraphAsymmErrors* g_Efficiency_dR = new TGraphAsymmErrors(  h_bmatched_function_dR , h_accepted_function_dR, "n" ) ;
		// g_Efficiency_dR->Draw() ;
		
		// c_Eff_dR_b_quarks->cd(3);
		// h_accepted_function_dR_pos->Draw("hist") ;
		// h_bmatched_function_dR_pos->Draw("hist same") ;
		
		// c_Eff_dR_b_quarks->cd(4);
		// TGraphAsymmErrors* g_Efficiency_dR_pos = new TGraphAsymmErrors(h_bmatched_function_dR_pos , h_accepted_function_dR_pos, "n") ;
		// g_Efficiency_dR_pos->Draw() ;
		
		// c_Eff_dR_b_quarks->cd(5);
		// h_accepted_function_dR_neg->Draw("hist") ;
		// h_bmatched_function_dR_neg->Draw("hist same") ;
		
		// c_Eff_dR_b_quarks->cd(6);
		// TGraphAsymmErrors* g_Efficiency_dR_neg = new TGraphAsymmErrors(h_bmatched_function_dR_neg, h_accepted_function_dR_neg ,"n") ;
		// g_Efficiency_dR_neg->Draw() ;

		// c_Eff_dR_b_quarks->SaveAs("Efficiency_dR "+ s1 +"_" + Jet_Option + ".png") ;
		// c_Eff_dR_b_quarks->Clear() ;

		TCanvas* c_tagger = new TCanvas("c_tagger" , "Tagger Value" , 600, 800);
		// h_tagger_value->Draw() ;
		s_btagger->Add(h_best_btagger_value) ;
		s_btagger->Add(h_second_best_btagger_value) ;
		s_btagger->Draw() ;
		c_tagger->SaveAs("Btagger_Value_ "+ s1 +"_" + Jet_Option + ".png") ;

}



// Main Function
void Z_bb_cs()
{
	std::vector<TTree*> *v_trees = new std::vector<TTree*> ;
	// TFile *f_TTJets = new TFile("~/cms/hist/TTJets.root") ;
 //   	cout << "Found file TTJets.root" << endl;
 //   	TTree *t_TTJets = (TTree*) f_TTJets->FindObjectAny("events") ; v_trees->push_back(t_TTJets) ;

 //   	TFile *f_WJetsToLNu = new TFile("~/cms/hist/WJetsToLNu.root") ;
 //   	cout << "Found file WJetsToLNu.root" << endl;
 //   	TTree* t_WJetsToLNu = (TTree*) f_WJetsToLNu->FindObjectAny("events") ; v_trees->push_back(t_WJetsToLNu) ;

   	TFile *f_ZZTo2Q2Nu = new TFile("~/cms/hist/ZZTo2Q2Nu.root") ;
   	cout << "Found file ZZTo2Q2Nu.root" << endl;
   	TTree* t_ZZTo2Q2Nu = (TTree*) f_ZZTo2Q2Nu->FindObjectAny("events") ; v_trees->push_back(t_ZZTo2Q2Nu) ;

   	TFile *f_WZTo1L1Nu2Q = new TFile("~/cms/hist/WZTo1L1Nu2Q.root") ;
   	cout << "Found file WZTo1L1Nu2Q.root" << endl ;
   	TTree* t_WZTo1L1Nu2Q = (TTree*) f_WZTo1L1Nu2Q->FindObjectAny("events") ; v_trees->push_back(t_WZTo1L1Nu2Q) ;

   	// TFile *f_WWToNuQQ = new TFile("~/cms/hist/WWToNuQQ.root") ;
   	// cout << "Found file WWToNuQQ.root" << endl;
   	// TTree* t_WWToNuQQ = (TTree*) f_WWToNuQQ->FindObjectAny("events") ; v_trees->push_back(t_WWToNuQQ) ;

   	// Plot the nJets distributions in the signal and background region to see what is a good cut to apply
   	//return the desired njets cut value
   	// Double_t njets_cut = Get_Cut_Parameters( v_trees , "njets" , 9, 0 , 10 ) ;
   	// Double_t met_cut = Get_Cut_Parameters(v_trees , "met" , 100, 0, 1000) ;

   	// information on each cut should be packed in the form Name,Min,Max,Increments
   	// std::vector<TString>  *v_cut_names = new std::vector<TString> ; v_cut_names->push_back("njets") ; v_cut_names->push_back("met") ;
   	// std::vector<Double_t> *v_cut_min = new std::vector<Double_t> ;
   	// std::vector<Double_t> *v_cut_max = new std::vector<Double_t> ;
   	// std::vector<Int_t> *v_cut_increments = new std::vector<Int_t> ;
   	// v_cut_names->push_back("njets") ;  v_cut_min->push_back(0) ; v_cut_max->push_back(10) ; v_cut_increments->push_back(10) ;
   	// v_cut_names->push_back("met") ;  v_cut_min->push_back(0) ;  v_cut_max->push_back(1000); v_cut_increments->push_back(20) ;

   	//Minimize_Sensitivity(v_trees, v_cut_names, v_cut_min, v_cut_max , v_cut_increments) ;

	// Signal_Z_bb_jets(v_trees);

	// N_jets_and_Zbb(v_trees );

	// Plot_All_Efficiency(v_trees) ;

	// Plot_dR_between_b_and_jet(v_trees );

	// Plot_n_jets_in_Signal(v_trees) ;
	Get_Efficiency_func( v_trees , "No_Tag" , kTRUE) ;
	// Get_Efficiency_func( v_trees , "No_Tag" , kFALSE) ;
	// Get_Efficiency_func(v_trees , "Loose" , kTRUE) ;
	// Get_Efficiency_func(v_trees , "Medium" , kTRUE) ;
	// Get_Efficiency_func(v_trees , "Tight" , kTRUE) ;
	// Get_Efficiency_func( v_trees , "No_Tag" , kFALSE) ;

}
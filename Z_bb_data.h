#ifndef Z_bb_data_H
#define Z_bb_data_H

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


Double_t Get_Reco_Mass_Data(std::vector<Int_t>* v_two_b_jets_selected , TClonesArray* jetP4)
{
	Double_t mass = 0  ;
	TLorentzVector* jet1 = dynamic_cast<TLorentzVector*> (jetP4->At( v_two_b_jets_selected->at(0) )) ;
	TLorentzVector* jet2 = dynamic_cast<TLorentzVector*> (jetP4->At( v_two_b_jets_selected->at(1) )) ;
	TLorentzVector *Z = new TLorentzVector() ;
	 *Z= *jet1 + *jet2 ;
	mass = Z->M() ;
return mass ;
}

void Plot_hist_Reco_Mass( TH1D* h_reco_mass , TString s1)
{
	h_reco_mass->SetTitle("Reconsructed Mass from 2 bjets") ;
	h_reco_mass->GetXaxis()->SetTitle( "Mass/GeV " ) ;
	h_reco_mass->GetYaxis()->SetTitle("Events/Gev") ;
	TCanvas* c = new TCanvas("c", "c_mass" , 600, 800) ;
	h_reco_mass->SetMarkerStyle(8) ;
	h_reco_mass->SetMarkerColor(kBlack) ;
	h_reco_mass->Draw("EP") ;
	c->SaveAs("/tmp/bbachu/h_mass"+s1+".png") ;
}

void Plot_hist_btagger(TH1D* h_tagger_value , TString s1 ) 
{

}

// This function loops over the jet tagger given and returns the indices of the best jets that pass the b tagging requirement
Bool_t Get_best_2_btagged_jets( std::vector<float>* btagger , Double_t Working_Point , std::vector<Int_t>* v_two_b_jets , TClonesArray *jetP4 , Double_t jet_pt_cut) 
{
	// Note: for an event to reach this stage it must already have more than 2 jets
	Int_t best_jet = 0 ; Double_t best_jet_tag = 0 ;
	Int_t second_best_jet = 0 ; Double_t second_best_jet_tag = 0 ;
	Int_t b_jet_count = 0 ;
	Bool_t two_b_jets = kFALSE ;
	for (int i = 0; i < btagger->size(); ++i)
	{
		// first check that the btagger passes the required working point and the jet_pt requirement
		TLorentzVector* jet1 = dynamic_cast<TLorentzVector*> (jetP4->At(i)) ;
		cout << "Value of btagger = " << btagger->at(i) << endl;
		if (( btagger->at(i) < Working_Point ) || (jet1->Pt() < jet_pt_cut)) continue ; // then it is not a b-jet
		cout << "Btagger was above working point of "<< Working_Point << endl ;
		b_jet_count = b_jet_count+1 ;
		cout << "Number of btagged jets = " << b_jet_count << endl ;
		// if it is the first jet that passes the working point then it is automatically the best jet
		if ( i == 0)
		{
			best_jet = 0 ;
			best_jet_tag = btagger->at(i) ;
		}
		// if it is not the first jet then we know that there should exist a best and a second best
		else
		{
			if ( btagger->at(i) > best_jet_tag ) // then this jet should be the best jet
			{
				//demote the best jet to second best jet
				second_best_jet = best_jet ;
				second_best_jet_tag = best_jet_tag ;
				// promote this jet as the best jet
				best_jet = i ;
				best_jet_tag = btagger->at(i) ;
			}
			else if ( btagger->at(i) > second_best_jet_tag )
			{
				// set this as the second best tagged jet
				second_best_jet = i ;
				second_best_jet_tag = btagger->at(i) ;
			}
		}
	}
	if ( b_jet_count > 1)
	{
		cout << "More than 2 jets were btagged. Taking the jet indices " << best_jet << " and " << second_best_jet  << endl ;
		two_b_jets = kTRUE ;
		v_two_b_jets->push_back(best_jet) ;
		v_two_b_jets->push_back(second_best_jet) ; 
	}

return two_b_jets ;
}


Bool_t Get_jets_closest_to_Z_mass( std::vector<float> *btagger ,  TClonesArray* jetP4 , Double_t Working_Point , std::vector<Int_t> *v_two_b_jets  , Double_t jet_pt_cut) 
{

	Double_t true_Z_mass = 91.1876 ;
	Double_t Best_Z_mass = 0 ;
	Int_t b_jet_count  = 0 ;
	Bool_t two_b_jets = kFALSE ;
	Int_t n_jets = btagger->size() ;
	Double_t mass = 0 ;
	int i_best_jet_1 = 0 ;
	int i_best_jet_2 = 0 ;

	for (int i = 0; i < n_jets; ++i)
	{
		TLorentzVector *jet1 = new TLorentzVector();
		jet1 = dynamic_cast<TLorentzVector*>(jetP4->At(i)) ;
		if ((btagger->at(i) < Working_Point) && (jet1->Pt() < jet_pt_cut)) continue ;
		b_jet_count = b_jet_count + 1 ;
		for (int j = i+1 ; j < n_jets; ++j)
		{
			TLorentzVector* jet2 = new TLorentzVector(); 
			jet2= dynamic_cast<TLorentzVector*>(jetP4->At(j)) ;
			if ((btagger->at(j) < Working_Point)||(jet2->Pt() < jet_pt_cut)) continue ;
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
	// cout << "Final mass selected = " << Best_Z_mass << endl ;
	if ( b_jet_count > 1  )
	{
		// cout << "Pushing back these indices into the vector " << i_best_jet_1 << " and "  << i_best_jet_2 << endl;
		v_two_b_jets->push_back(i_best_jet_1) ;
		v_two_b_jets->push_back(i_best_jet_2) ;
		two_b_jets = kTRUE ;
		// cout << "More than 2 bjets" << endl;
	}
	else
	{
		// cout << "Jets failed requirements for btagging " << endl;
	}
return two_b_jets ;
}


void Draw_all_Data(TString name , TH1D* h_loose ,TH1D* h_medium ,TH1D* h_tight )
{
	h_loose->SetMarkerColor(kBlue) ;
	h_medium->SetMarkerColor(kRed) ;
	h_tight->SetMarkerColor(kPink) ;

	TLegend* leg = new TLegend(0.7, 0.9, 0.90, 0.7) ;
	leg->AddEntry(h_loose , "Loose" , "p") ;
	leg->AddEntry(h_medium , "Medium" , "p") ;
	leg->AddEntry(h_tight , "Tight" , "p") ;

	TCanvas* c = new  TCanvas("c" , name , 600, 800) ;
	h_loose->Draw("EP ") ;
	h_medium->Draw("EP same") ;
	h_tight->Draw("EP same") ;
	leg->Draw("same");
	c->SaveAs("/tmp/bbachu."+name+".png") ; 

}


TH1D* Get_Hist_Reco_Mass_Data(std::vector<TTree*> *v_Data_trees ,  TString Jet_Selection_Option, TString Working_Point_Choice , Double_t jet_pt_cut)
{
	cout << "Plotting the reconstructed Mass from the bb-jets in data"  <<endl ;
	Double_t Working_Point = 0 ;
	if (Working_Point_Choice == "No_Tag")
	{
		Working_Point = 0 ;
	}
	else if(Working_Point_Choice == "Loose" )
	{
		Working_Point = 0.432	 ;
	}
	else if (Working_Point_Choice == "Medium")
	{
		Working_Point = 0.841 ;
	}	
	else if (Working_Point_Choice == "Tight")
	{
		Working_Point = 0.941 ;
	}
	else 
	{
		Working_Point = 0 ;
		cout << " !!! Error in working points !!!" << endl;
	}

	TTree* tree = new TTree();
	// INCOPERATE WORKING POINT IN THIS STRING NAME
	TString s1  = Jet_Selection_Option ;
	TString s2 =  TString::Format("_%.3f_%.3f",Working_Point, jet_pt_cut) ;
	TString s = s1 + s2 + Working_Point_Choice;

	TH1D* h_tagger_value = new TH1D("h_tagger_value" , " Value of b-tagger on b-matched Jets" + s1 +"_" + Jet_Selection_Option  , 20 , 0, 1) ;
	
	TH1D* h_reco_mass = new TH1D("h_reco_mass" , "h_reco_mass" , 100, 0, 200) ;

	cout << "Begin looping over Trees" << endl ;
	// loop over the trees
	for (unsigned int i = 0; i < v_Data_trees->size(); ++i)
	{
		tree = v_Data_trees->at(i) ;	

		std::vector<float> *btagger = new std::vector<float> ;	
		TClonesArray *jetP4 = new TClonesArray() ;

		float mcWeight;

		tree->SetBranchAddress("jetP4", &jetP4) ;
		tree->SetBranchAddress("CombinedInclusiveSecondaryVertexV2", &btagger) ;

		Double_t nentries = tree->GetEntries() ;
		// Decrease the number of entries for testing to 0.05%
		nentries = nentries*0.5;

		std::vector<Int_t> *v_2_b_jets_selected = new std::vector<Int_t> ;
		// cout << "Looping over " << nentries << " events in dataset " << i << endl ;

		for ( unsigned int j = 0; j < nentries; ++j)
		{
			tree->GetEntry(j) ;
			cout << "Entry = " << j << endl ;

			// Make sure that there are at least 2 jets in the event to reconstruct the mass from data
			if ( jetP4->GetEntries() < 2)
			{
				cout << "Less than 2 jets in this event" << endl ;
			}
			// Now pick out the best 2 jets that pass the jet selection criteria
			// The event will be skipped if more than 2 bjets are not identified
			if ( Jet_Selection_Option == "Best_btag" )
			{
				if (Get_best_2_btagged_jets( btagger ,  Working_Point ,  v_2_b_jets_selected ,jetP4 , jet_pt_cut) == kFALSE ) continue ;
			}
			else if ( Jet_Selection_Option == "Best_Z_Mass" )
			{
				if (Get_jets_closest_to_Z_mass( btagger , jetP4,  Working_Point , v_2_b_jets_selected , jet_pt_cut) == kFALSE) continue ;
			}
			else
			{
				cout << " !!! Wrong jet selection option selected !!!" << endl ;
			}
			// Now that I hve the incdices of the best 2 jets that fullfill my reqirements I can plot the reconstructed Z_mass
			Double_t Reco_mass = Get_Reco_Mass_Data( v_2_b_jets_selected , jetP4) ;
			// cout << "Produced Z mass from 2 jets = " << Reco_mass << endl ;
			h_reco_mass->Fill( Reco_mass ) ;

			// The Distribution of the btagger will also be useful information 
			h_tagger_value->Fill( btagger->at( v_2_b_jets_selected->at(0) ) ) ;
			h_tagger_value->Fill( btagger->at( v_2_b_jets_selected->at(1) ) ) ;

			v_2_b_jets_selected->clear() ;
		}
		
		Plot_hist_Reco_Mass( h_reco_mass , s ) ;
		Plot_hist_btagger( h_tagger_value , s ) ;
	}

return h_reco_mass ;
}

#endif
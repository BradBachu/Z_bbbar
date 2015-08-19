#ifndef Z_bb_MC_h
#define Z_bb_MC_h

#include "TAxis.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TImage.h"
#include "TROOT.h"
#include "TChain.h"
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




void Get_Efficiency_func(std::vector<TChain*> *v_trees , TString Jet_Option, Bool_t bMatched)
{
	cout << "Looking at Efficiency with intermediate bins and jet option " << Jet_Option <<endl ;
	// TTree* tree = new TTree();
	TChain *chain = new TChain() ;
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
	// Make hists for n accepted as a function of Z_pt used for EFFICIENCY plots
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
		chain = v_trees->at(i) ;	
		std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
		std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
		std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
		std::vector<float> *btagger = new std::vector<float> ;	
		TClonesArray *jetP4 = new TClonesArray() ;
		TClonesArray *genP4 = new TClonesArray() ; 
		TClonesArray *metP4 = new TClonesArray() ;
		float mcWeight;
		chain->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId) ;
		chain->SetBranchAddress("genPdgId", &genPdgId) ;
		chain->SetBranchAddress("genMotherpdgId", &genMotherPdgId) ;
		chain->SetBranchAddress("jetP4", &jetP4) ;
		chain->SetBranchAddress("genP4", &genP4) ;
		chain->SetBranchAddress("mcWeight",&mcWeight) ;
		chain->SetBranchAddress("metP4", &metP4) ;
		chain->SetBranchAddress("CombinedInclusiveSecondaryVertexV2", &btagger) ;

		Double_t nentries = chain->GetEntries() ;
		// Decrease the number of entries for testing to 0.05%
		nentries = nentries*0.01 ;
		// cout << "Looping over " << nentries << " events in dataset " << i << endl ;

		for ( unsigned int j = 0; j < 50; ++j)
		{
			cout << "Looking at event = " << j << endl ;
			// std::vector<Double_t> *v_njets_btag_value = new std::vector<Double_t> ;
			chain->GetEvent(j) ;
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
	


		TCanvas* c = new TCanvas("c_Efficiency_Zpt" + s1 +"_" + Jet_Option, "c_Efficiency_Zpt" + s1 +"_" + Jet_Option , 1000, 2000) ;
		c->Divide(2,3) ;
		
		c->cd(1) ;
		h_accepted_funtion_Z_pt->Draw("hist") ;
		h_bmatched_function_Z_pt->Draw("same hist") ;
		
		c->cd(2) ;

		TGraphAsymmErrors* g_Efficiency = new TGraphAsymmErrors( h_bmatched_function_Z_pt , h_accepted_funtion_Z_pt , "n"  );
		g_Efficiency->Draw() ;
		
		c->cd(3);
		h_accepted_funtion_Z_pt_pos->Draw("hist") ;
		h_bmatched_function_Z_pt_neg->Draw("same hist") ;
		
		c->cd(4) ;
		TGraphAsymmErrors* g_Efficiency_pos = new TGraphAsymmErrors( h_bmatched_function_Z_pt_pos , h_accepted_funtion_Z_pt_pos , "n" ) ;
		g_Efficiency_pos->Draw() ;
	
		c->cd(5) ;
		h_accepted_funtion_Z_pt_neg->Draw() ;
		h_bmatched_function_Z_pt_neg->Draw("same") ;
		
		c->cd(6) ;
		TGraphAsymmErrors* g_Efficiency_neg = new TGraphAsymmErrors(h_bmatched_function_Z_pt_neg , h_accepted_funtion_Z_pt_neg, "n") ;
		g_Efficiency_neg->Draw() ;

		c->SaveAs("Efficiency_Zpt"+ s1 +"_" + Jet_Option+".png") ;
		c->Clear() ;
		// TString name = "Efficiency" + Jet_Option + "bmatched.png" ;
		// c->SaveAs(name) ;

		TCanvas *c_Eff_dR_b_quarks = new TCanvas("c_Eff_dR_b_quarks" + s1 +"_" + Jet_Option, "c_Efficiency_Eff_dR_b_quarks" + s1 +"_" + Jet_Option, 1000, 2000) ;
		c_Eff_dR_b_quarks->Divide(2,3) ;

		c_Eff_dR_b_quarks->cd(1);
		h_accepted_function_dR->Draw("hist") ;
		h_bmatched_function_dR->Draw("hist same") ;
		
		c_Eff_dR_b_quarks->cd(2);
		TGraphAsymmErrors* g_Efficiency_dR = new TGraphAsymmErrors(  h_bmatched_function_dR , h_accepted_function_dR, "n" ) ;
		g_Efficiency_dR->Draw() ;
		
		c_Eff_dR_b_quarks->cd(3);
		h_accepted_function_dR_pos->Draw("hist") ;
		h_bmatched_function_dR_pos->Draw("hist same") ;
		
		c_Eff_dR_b_quarks->cd(4);
		TGraphAsymmErrors* g_Efficiency_dR_pos = new TGraphAsymmErrors(h_bmatched_function_dR_pos , h_accepted_function_dR_pos, "n") ;
		g_Efficiency_dR_pos->Draw() ;
		
		c_Eff_dR_b_quarks->cd(5);
		h_accepted_function_dR_neg->Draw("hist") ;
		h_bmatched_function_dR_neg->Draw("hist same") ;
		
		c_Eff_dR_b_quarks->cd(6);
		TGraphAsymmErrors* g_Efficiency_dR_neg = new TGraphAsymmErrors(h_bmatched_function_dR_neg, h_accepted_function_dR_neg ,"n") ;
		g_Efficiency_dR_neg->Draw() ;

		c_Eff_dR_b_quarks->SaveAs("Efficiency_dR "+ s1 +"_" + Jet_Option + ".png") ;
		c_Eff_dR_b_quarks->Clear() ;

		TCanvas* c_tagger = new TCanvas("c_tagger" , "Tagger Value" , 600, 800);
		// h_tagger_value->Draw() ;
		s_btagger->Add(h_best_btagger_value) ;
		s_btagger->Add(h_second_best_btagger_value) ;
		s_btagger->Draw() ;
		c_tagger->SaveAs("Btagger_Value_ "+ s1 +"_" + Jet_Option + ".png") ;
	}
}


#endif
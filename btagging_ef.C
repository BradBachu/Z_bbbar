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

//get the jets with the highest pt
std::vector<int> Highest_pt_pair_index(TClonesArray* jetarray )
{
	int size_of_array = jetarray->GetEntries() ; //cout << "Number of jets = " << size_of_array << endl ;
	int pt1_index = 0;
	int pt2_index = 0;
	float highest_jet_pt = 0;
	float second_highest_jet_pt = 0;

	//loop throught the jets in the array
	for (int j = 0; j < size_of_array; ++j)
	{
		// cout << "j = " << j << endl;
		TLorentzVector* jet = dynamic_cast<TLorentzVector*>(jetarray->At(j)) ;
		// TLorentzVector* jet = (TLorentzVector*) jetarray->At(j) ;
		float jet_pT = jet->Pt();
		if (j == 0) //the fist pt is the highest and second highest
		{
			pt1_index = 0 ; pt2_index = 0 ;
			highest_jet_pt = jet_pT ; second_highest_jet_pt = jet_pT ;

		}
		else if (j == 1) //now one has to be higher than the other
		{ 
			if (jet_pT > highest_jet_pt ) // now this is the highest pt jet
			{
				//first demote the hightest jet to second highest
				second_highest_jet_pt = highest_jet_pt ;
				pt2_index = pt1_index ;
				highest_jet_pt = jet_pT ; 
				pt1_index = 1 ;

			}
			else // well then we will set it to be the second highest pt
			{
				second_highest_jet_pt = jet_pT ;
				pt2_index = 1;
			}
		}
		else //now we have more than 2 jets so we need to check if they are any good
		{
			if ((jet_pT > second_highest_jet_pt) && ( jet_pT < highest_jet_pt)) //then it should become the new second hightest jet
			{
				second_highest_jet_pt = jet_pT ;
				pt2_index = j ; 
			}
			if (jet_pT > highest_jet_pt) //then it should become the highest jet and the previous one should be demoted
			{
				second_highest_jet_pt = highest_jet_pt ;
				pt2_index = pt1_index ;
				highest_jet_pt = jet_pT ;
				pt1_index = j ;
			}
		}
	// cout << "highest is at index " << pt1_index << " with p = "<< highest_jet_pt << endl ;
	// cout << "second is at index " << pt2_index << " with p = " <<second_highest_jet_pt << endl;
	}
	std::vector<int> highest_pt_pair_index;
	highest_pt_pair_index.push_back(pt1_index) ; 
	highest_pt_pair_index.push_back(pt2_index) ;
return highest_pt_pair_index;
}

Bool_t two_b_quarks(  std::vector<Double_t> *genPdgId) 
{
	Bool_t two_b_quarks ;
	int n_b_quarks = 0 ;
	for (int i = 0; i < genPdgId->size(); ++i)
	{
		if ( abs(genPdgId->at(i)) == 5 )
		{
			n_b_quarks = n_b_quarks + 1 ;
		}
	}
	if (n_b_quarks > 1)
	{
		two_b_quarks = kTRUE ;
		// cout << "Has 2 b-quarks"<< endl;
	}
	else
	{
		two_b_quarks = kFALSE ;
		// cout <<"Doesnt have 2 b-quarks"<< endl;
	}
return two_b_quarks ;
}

Bool_t Z_pT_Range(TClonesArray* jetarray , Double_t Z_pT_Min, Double_t Z_pT_Max)
{
	Bool_t In_range;
	std::vector<int> highest_pt_pair_index =  Highest_pt_pair_index(jetarray);
	TLorentzVector* jet1 = dynamic_cast<TLorentzVector*>(jetarray->At(highest_pt_pair_index.at(0))) ;
	TLorentzVector* jet2 = dynamic_cast<TLorentzVector*>(jetarray->At(highest_pt_pair_index.at(1))) ;
	highest_pt_pair_index.clear() ;
	TLorentzVector Z = *jet1 + *jet2 ;
	Double_t Z_pt = Z.Pt();
	if (  (Z_pt > Z_pT_Max) || (Z_pt < Z_pT_Min) ) 
	{
		In_range = kFALSE ;
		// cout << "Not in Range"<<endl;
	}
	else 
	{
		In_range = kTRUE ;
		// cout << "In Range "<< endl;
	}
return In_range ;
}

Bool_t Z_pT_Range_gen_level(std::vector<Double_t> *genPdgId , TClonesArray *genP4, std::vector<Double_t> *genMotherPdgId , Double_t Z_pT_Min, Double_t  Z_pT_Max)
{
	std::vector<int> matched_idex;
	Bool_t In_range;
	//loop through the genPdgId's and Mothers to see which b-quarks are from Z's
	for (int i = 0; i < genPdgId->size(); ++i)
	{
		if ( abs(genPdgId->at(i)) == 5 && genMotherPdgId->at(i) == 23 )	
		{
			matched_idex.push_back(i);
		}
	}
	TLorentzVector* gen1_b_P4 = dynamic_cast<TLorentzVector*>(genP4->At(matched_idex.at(0))) ;
	TLorentzVector* gen2_b_P4 = dynamic_cast<TLorentzVector*>(genP4->At(matched_idex.at(1))) ;
	TLorentzVector Z = *gen1_b_P4 + *gen2_b_P4 ;
	if ( (Z.Pt() > Z_pT_Min) && (Z.Pt() < Z_pT_Max) )
	{
		In_range = kTRUE ;
	}
	else
	{
		In_range = kFALSE;
	}
return In_range;
}

int matched_b_to_Z(std::vector<Double_t> *genPdgId, std::vector<Double_t> *genMotherPdgId)
{
	int count = 0 ;
	for (int i = 0; i < genPdgId->size(); ++i)
	{
		if ( (abs(genPdgId->at(i)) ==5) && (genMotherPdgId->at(i)== 23))
		{
			count = count + 1;
		}
	}
return count;
}

Double_t Get_Acceptance_in_ZpT_Range(Double_t Z_pT_Min , Double_t Z_pT_Max)
{
	TFile *f =new TFile("../../nero.root","READ");
    TTree *tree = (TTree*)f->FindObjectAny("events");
    cout<< "Getting acceptance measurements in range "<<Z_pT_Min << " " << Z_pT_Max << endl;
	
	std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
	std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
	std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
	TClonesArray *jetP4 = new TClonesArray();
	TClonesArray *genP4 = new TClonesArray();

	tree->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId) ;
	tree->SetBranchAddress("genPdgId", &genPdgId) ;
	tree->SetBranchAddress("genMotherPdgId", &genMotherPdgId) ;
	tree->SetBranchAddress("jetP4", &jetP4) ;
	tree->SetBranchAddress("genP4", &genP4) ;

	//create variables to keep track of bquark acceptance ratio
	Double_t total_gen_b = 0 ;
	Double_t detector_acceptance_gen_b = 0 ;

	for (int i = 0; i < tree->GetEntries(); ++i)
	{
		tree->GetEntry(i) ;
		//make sure there are at least 2 b-quarks
		if (two_b_quarks(genPdgId) == kFALSE) continue;
		//make sure there are at least 2 jets
		if (matched_b_to_Z(genPdgId, genMotherPdgId) < 2 ) continue;
		// if ( jetP4->GetEntries() == 0 || jetP4->GetEntries() == 1) continue;
		// make sure it is in the Z_pt range
		// if (Z_pT_Range(jetP4 ,  Z_pT_Min,  Z_pT_Max) == kFALSE) continue;
		if (Z_pT_Range_gen_level(genPdgId , genP4, genMotherPdgId ,  Z_pT_Min,  Z_pT_Max) == kFALSE) continue;
		// cout  << "Main criteria passed" << endl;
		//loop over vector in the entry 
		for (int j = 0; j < genPdgId->size(); ++j)
		{
			// if it is a b-quark and came from a Z
			if (abs((genPdgId->at(j) == 5)) && (genMotherPdgId->at(j) == 22))
			{
				// cout << "b-quark came from Z"<< endl;
				total_gen_b = total_gen_b + 1 ;
				TLorentzVector* gen_b_P4 = dynamic_cast<TLorentzVector*>(genP4->At(j)) ;
				//check if it is in the fidicual acceptance
				if ( abs(gen_b_P4->Eta() ) < 2.4 )
					{
						// cout << "b-quark in acceptance"<< endl;
						detector_acceptance_gen_b = detector_acceptance_gen_b + 1 ;
					}
			}
		}
	}
		cout << "total = "<< total_gen_b << endl;
		cout << "accepted = " << detector_acceptance_gen_b << endl;
		Double_t ratio;
		if (total_gen_b == 0)
		{
			ratio = 0 ;
		}
		else
		{
		ratio = detector_acceptance_gen_b / total_gen_b ;
		}
	
	cout << "ratio = " << ratio << endl;
return ratio ;
}

void Acceptance_function_Z_pT(Double_t Z_pt_bins[], int size)
{
	TH1F* h_acceptance_ratio = new TH1F("h_acceptance_ratio", "Acceptance of Z#Rightarrow b #bar{b}", size-1 , Z_pt_bins);
	for (int i = 0; i < size-1; ++i)
	{
		h_acceptance_ratio->SetBinContent( (i+1) ,  Get_Acceptance_in_ZpT_Range(Z_pt_bins[i] , Z_pt_bins[i+1]) );
	}
	h_acceptance_ratio->GetXaxis()->SetTitle("Z p_{#perp}");
	h_acceptance_ratio->GetYaxis()->SetTitle("Acceptance");
	TCanvas *c = new TCanvas();
	h_acceptance_ratio->Draw("HIST ");
	c->SaveAs("Acceptance.png");
	c->SaveAs("Acceptance.pdf");
}

std::vector<Double_t> Discrimination_limits(TString btagger)
{
	std::vector<Double_t> v;
	Double_t min; Double_t max;
	if (btagger == "JetProbability")                      { min =  0 ; max = 3.5 ; }
	if (btagger == "JetBProbability")                     { min =  0 ; max = 10 ; }
	if (btagger == "SimpleSecondaryVertexHighEff")        { min =  1 ; max = 5 ; }
	if (btagger == "SimpleSecondaryVertexHighPur")        { min =  1 ; max = 5 ; }
	if (btagger == "CombinedSecondaryVertex")             { min =  0 ; max = 1 ; }
	if (btagger == "CombinedSecondaryVertexV2")           { min =  0 ; max = 1 ; }
	if (btagger == "CombinedSecondaryVertexSoftLepton")   { min =  0 ; max = 1 ; }
	if (btagger == "CombinedInclusiveSecondaryVertexV2")  { min =  0 ; max = 1 ; }
	if (btagger == "CombinedMVA")                         { min =  0 ; max = 1 ; }
	//pushback the starting point of the discriminator
	v.push_back(min);
	//pushback the range
	v.push_back(max - min) ;
return v;
}

//make a function of Z pT where the Zpt is constructed from the mass of the 2 leading quarks
// the pt of the Z is the sum of the 4 vectors of the two leding jets in an event
std::vector< std::vector<Double_t> >  Function_of_Z_pt( TString btagger , Double_t n_coordinates, Double_t Z_pT_Min , Double_t Z_pT_Max )
{

	//Create the vectors to store the Efficiency and Fake rate
	std::vector<Double_t> v_eff; v_eff.clear() ;
	std::vector<Double_t> v_fake; v_fake.clear();
	//Create vector to store the vectors of coordinates
	std::vector<std::vector<Double_t> > v_coordinate; v_coordinate.clear() ;

	// TLorentzVector* jetp4 = dynamic_cast<TLorentzVector*> obj ;

	TFile *f =new TFile("../../nero.root","READ");
    TTree *tree = (TTree*)f->FindObjectAny("events");
	
	std::vector<float> *var_btagger = new std::vector<float> ;	
	std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
	std::vector<Double_t> *genPdgId = new std::vector<Double_t>;
	std::vector<Double_t> *genMotherPdgId = new std::vector<Double_t>;
	TClonesArray *jetp4 = new TClonesArray();

	tree->SetBranchAddress(btagger, &var_btagger) ;
	tree->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId) ;
	tree->SetBranchAddress("genPdgId", &genPdgId) ;
	tree->SetBranchAddress("genMotherPdgId", &genMotherPdgId) ;
	tree->SetBranchAddress("jetP4", &jetp4) ;

	TLorentzVector* jet1 = new TLorentzVector() ;
	TLorentzVector* jet2 = new TLorentzVector() ;
	std::vector<Double_t> v_discrimination_limits = Discrimination_limits(btagger);
	Double_t discrimination_strength = v_discrimination_limits.at(0) ;
	Double_t increment = v_discrimination_limits.at(1) / (n_coordinates-1) ;
	for (int i_coordinate = 0; i_coordinate < n_coordinates; ++i_coordinate)
	{
		Double_t n_gen_bjets = 0 ;
		Double_t n_gen_light_jets = 0 ;
		Double_t n_btagged_jets_given_bgen = 0 ;
		Double_t n_btagged_given_lightgen = 0 ;
		int count = 0 ;
		Double_t n_entries = tree->GetEntries(); 
		for (int i = 0; i < n_entries ; ++i)
		{
			tree->GetEntry(i) ;
			std::vector<int> highest_pt_pair_index ;
			if (!( jetMotherPdgId->size() == jetp4->GetEntries() ))
			{
			cout << "Size of jetMotherPdgId = " << jetMotherPdgId->size() << ". Size of jetsP4 = " << jetp4->GetEntries() << " at entry " << i<<endl;
			}
			if ( jetp4->GetEntries() == 0 || jetp4->GetEntries() == 1) continue;
			//determine which jets have the highest pt
			highest_pt_pair_index =  Highest_pt_pair_index(jetp4);
			int jet1_index = highest_pt_pair_index.at(0) ;
			int jet2_index = highest_pt_pair_index.at(1) ;
			TLorentzVector* jet1 = dynamic_cast<TLorentzVector*>(jetp4->At(jet1_index)) ;
			TLorentzVector* jet2 = dynamic_cast<TLorentzVector*>(jetp4->At(jet2_index)) ;
			highest_pt_pair_index.clear() ;
			TLorentzVector Z = *jet1 + *jet2 ;
			Double_t Z_pt = Z.Pt();
			if (  (Z_pt > Z_pT_Max) || (Z_pt < Z_pT_Min) ) continue ;	
			count = count + 1 ;
			for( int j = 0 ; j < jetMotherPdgId->size() ; ++j)
			{
				// Eff = P(pass tag | b gen)
				if (abs(jetMotherPdgId->at(j)) == 5)
				{
					// cout << "Gen ID identified as b quark " << endl;
					n_gen_bjets = n_gen_bjets  + 1 ;
					//now check if it was btagged
					if ( var_btagger->at(j) > discrimination_strength )
					{
						// cout << "Actual b-jet tagged as b-jet" << endl;
						n_btagged_jets_given_bgen  = n_btagged_jets_given_bgen + 1 ;
					}					
				}
				// Fake = P(pass tag | light gen)
				if (!(abs(jetMotherPdgId->at(j)) == 5))
				{
					// cout << "Gen ID not b quark"<< endl;
					n_gen_light_jets = n_gen_light_jets + 1 ;
					// now check if it was btagged
					if ( var_btagger->at(j) > discrimination_strength )
					{
						// cout << "Not b-jet identified as b-jet"<< endl;
						n_btagged_given_lightgen = n_btagged_given_lightgen + 1 ;
					}
				}				
			}
		}
		discrimination_strength = discrimination_strength +  increment ;
		Double_t efficiency = n_btagged_jets_given_bgen / n_gen_bjets ; 
		v_eff.push_back(efficiency) ;
		Double_t fake_rate = n_btagged_given_lightgen / n_gen_light_jets ; 
		v_fake.push_back(fake_rate) ;		
		cout << "Got coordinate " << i_coordinate << endl;
		if ( Z_pT_Min == 500 && Z_pT_Max == 1000 )
		{
			cout<< "(Efficiency, Fake rate) : ( " << efficiency << " , " << fake_rate << " ) with count " << count << endl;
		}
	}
	v_coordinate.push_back(v_fake) ;
	v_coordinate.push_back(v_eff) ;
return v_coordinate ;
}




//Use to import a ntuple from a directory
std::vector< std::vector<Double_t> >  v_Eff_Fak_Pair( TString btagger , Double_t n_coordinates )
{
	//Create the vectors to store the Efficiency and Fake rate
	std::vector<Double_t> v_eff; v_eff.clear() ;
	std::vector<Double_t> v_fake; v_fake.clear();
	//Create vector to store the vectors
	std::vector<std::vector<Double_t> > v_coordinate; v_coordinate.clear() ;

	std::vector<Double_t> *jetMotherPdgId  = new std::vector<Double_t>;
	std::vector<float> *var_btagger = new std::vector<float> ;	

	TFile *f =new TFile("../../nero.root","READ");
    TTree *tree = (TTree*)f->FindObjectAny("events");
	
	tree->SetBranchAddress(btagger, &var_btagger);
	tree->SetBranchAddress("jetMotherPdgId", &jetMotherPdgId);

	Double_t discrimination_strength = 0 ;
	Double_t increment = 1 / (n_coordinates-1) ;
	// cout << "n_coordinates = " << n_coordinates << endl;
	// cout << "Discrimination increment = " << increment << endl;
	Double_t n_entries = tree->GetEntries(); 

	for (int i_coordinate = 0; i_coordinate < n_coordinates; ++i_coordinate)
	{

		Double_t n_gen_bjets = 0 ;
		Double_t n_gen_light_jets = 0 ;
		Double_t n_btagged_jets_given_bgen = 0 ;
		Double_t n_btagged_given_lightgen = 0 ;

		// cout << "Discrimination Strength = "  << discrimination_strength << endl ;

		for (Double_t i = 0; i < n_entries ; ++i)
		{
			tree->GetEntry(i) ;
			// Now loop through the jets
			for( int j = 0 ; j < jetMotherPdgId->size() ; ++j)
			{
				// Eff = P(pass tag | b gen)
				if (abs(jetMotherPdgId->at(j)) == 5)
				{
					// cout << "Gen ID identified as b quark " << endl;
					n_gen_bjets = n_gen_bjets  + 1 ;
					//now check if it was btagged
					if ( var_btagger->at(j) > discrimination_strength )
					{
						// cout << "Actual b-jet tagged as b-jet" << endl;
						n_btagged_jets_given_bgen  = n_btagged_jets_given_bgen + 1 ;
					}					
				}
				// Fake = P(pass tag | light gen)
				if (!(abs(jetMotherPdgId->at(j)) == 5))
				{
					// cout << "Gen ID not b quark"<< endl;
					n_gen_light_jets = n_gen_light_jets + 1 ;
					// now check if it was btagged
					if ( var_btagger->at(j) > discrimination_strength )
					{
						// cout << "Not b-jet identified as b-jet"<< endl;
						n_btagged_given_lightgen = n_btagged_given_lightgen + 1 ;
					}
				}				
			}
		}
		discrimination_strength = discrimination_strength +  increment ;
		// cout << "n_btagged_jets_given_bgen = " <<  n_btagged_jets_given_bgen << "      n_btagged_given_lightgen = " << n_btagged_given_lightgen << endl;
		// cout << "n_gen_bjets = " << n_gen_bjets << "       n_gen_light_jets = " << n_gen_light_jets << endl;
		Double_t efficiency = n_btagged_jets_given_bgen / n_gen_bjets ; 
		v_eff.push_back(efficiency) ;
		Double_t fake_rate = n_btagged_given_lightgen / n_gen_light_jets ; 
		v_fake.push_back(fake_rate) ;		
		// cout << "Efficieny = " << efficiency << " Fake Rate = " << fake_rate << endl;
	}
	v_coordinate.push_back(v_fake) ;
	v_coordinate.push_back(v_eff) ;

	return v_coordinate ;
}


//loop through the btagger increments to get poDouble_ts for ROC curve
TGraph* Make_ROC_curve( TString btagger , Double_t n_coordinates, Double_t Z_pT_Min, Double_t Z_pT_Max)
{
	cout << "Making ROC curve for " << btagger << endl;
	// get the coordinates 
	std::vector< std::vector<Double_t> >  v_coordinate =  Function_of_Z_pt(  btagger , n_coordinates  , Z_pT_Min , Z_pT_Max);
	std::vector<double> vx = v_coordinate.at(0) ;
	std::vector<double> vy = v_coordinate.at(1) ;
	int size_a = vx.size(); //cout << "Vector size_a = " << size_a << endl;
	int size_b = vy.size(); //cout << "Vector size_b = " << size_b << endl;

	Double_t x[size_a] ;
	Double_t y[size_a] ;
	for (int k = 0; k < n_coordinates ; ++k)
	{
		x[k] = vx.at(k) ; 
		y[k] = vy.at(k) ;
		// cout << "Coordinate " << k << " ( " << x[k] << " , " << y[k] << " ) " << endl;
	}
 	TGraph* ROC_curve = new TGraph(  size_a , x , y ) ;
	vx.clear(); vy.clear();

 return ROC_curve;
}


void Make_ROC_Curves_for_all_taggers(TString btagger_list[], Int_t size, Double_t n_coordinates, Double_t Z_pT_Min, Double_t Z_pT_Max)
{
	TString low = TString::Format("%4.2f", Z_pT_Min);
	TString high = TString::Format("%4.2f", Z_pT_Max);
	TString title = "b-tagger ROC Curves " + low + " <  Z p_{#perp}  < " + high ;
	TMultiGraph* all_ROC_curves = new TMultiGraph("ROC", title);
	all_ROC_curves->SetTitle( title + "; False Positive Rate ; True Positive Rate");
	// all_ROC_curves->GetXaxis()->SetTitle("False Positive Rate") ;
	// all_ROC_curves->GetYaxis()->SetTitle("True Positive Rate") ;
	// all_ROC_curves->GetYaxis()->SetLimits(0,1);
	// all_ROC_curves->GetXaxis()->SetLimits(0,1);
	// Make a ROC curve for each tagger and add it to the multi graph
	for ( int i = 0 ; i < size ; ++i ) 
	{
		// declare a TGraph for a tagger
		TGraph* g1 =  Make_ROC_curve( btagger_list[i] , n_coordinates , Z_pT_Min, Z_pT_Max);
		g1->SetMarkerColor(1 + i) ;
		if( i +1 == 10)
		{
			g1->SetLineColor(40);
		}
		else
		{
		g1->SetLineColor(1+i);
		}
		g1->SetMarkerStyle(24 + i) ;
		g1->SetFillStyle(0);
		
		g1->SetTitle(btagger_list[i]);
		g1->SetName(btagger_list[i]);
		all_ROC_curves->Add(g1);
	}	
	TCanvas* c = new TCanvas("c","c",500,600);
	// TF1* f_guess = new TF1("f_guess","x",0,1);
	// f_guess->Draw();
	all_ROC_curves->Draw("ALP");
	TLegend* leg =	c->BuildLegend(0.6, 0.5 , 0.92, 0.14,  "b-taggers");
	TString name;
	name = "ROC_curves_"+ low + "_"+high ;
	c->SaveAs(name+".png");
	c->SaveAs(name+".pdf");
	TCanvas* c1 = new TCanvas("c1", "c1", 500, 600);
	all_ROC_curves->Draw("a fb pl3d");
	c1->SaveAs(name+"3D.png");
	c1->SaveAs(name+"3D.pdf");
}

void btagging_ef()
{
	Double_t n_coordinates = 11 ;
	Double_t Z_pt_bins[4] = {0, 250, 500, 1000} ;
	// Double_t Z_pt_bins[2] = {500, 1000} ;

	Double_t size =  9 ;
	TString btagger_list[9] = { "JetProbability", 
								"JetBProbability", 
								"SimpleSecondaryVertexHighEff" , 
								"SimpleSecondaryVertexHighPur", 
								"CombinedSecondaryVertex", 
								"CombinedSecondaryVertexV2", 
								"CombinedSecondaryVertexSoftLepton",
								"CombinedInclusiveSecondaryVertexV2", 
								"CombinedMVA"/*, 
								"TrackCountingHighEff" ,
								"TrackCountingHighPur" */
							  } ;


	// for (int i = 0; i < 3; ++i)
	// {
	// 	Make_ROC_Curves_for_all_taggers(btagger_list , size, n_coordinates, Z_pt_bins[i], Z_pt_bins[i+1]);
	// }
	Acceptance_function_Z_pT( Z_pt_bins,  4);
} 




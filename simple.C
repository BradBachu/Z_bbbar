void simple(){

  TFile *f =new TFile("../../nero.root","READ");
  // nero->cd();

  TTree *tree = (TTree*)f->FindObjectAny("events");

  for (int i = 0; i<tree->GetEntries(); ++i)
  {
    tree->GetEntry(i);
  }


  std::cout << "finished looping" << std::endl;

}

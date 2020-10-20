#include "ShowerVetoMaxDumb.h"
#include "ShowerVetoMax.h"
#include "ShowerMaxVetoDumb.h"
#include "ShowerMaxVeto.h"
#include "ShowerGenerateSelectDumb.h"
#include "ShowerGenerateSelect.h"
#include "ShowerMaxVetoSmart.h"
#include "TH1D.h"
#include "TFile.h"

class InputParser {

public:

  InputParser (int &argc, char **argv) {
    for (int i = 1; i < argc; ++i) arglist.push_back(string(argv[i]));
  }

  const string& getOption(const string &opt) const {
    vector<string>::const_iterator itr = find(arglist.begin(),
      arglist.end(), opt);
    if (itr != arglist.end() && ++itr != arglist.end()) return *itr;
    return "";
  }

  bool hasOption(const string &opt) const {
    return find(arglist.begin(), arglist.end(), opt) != arglist.end();
  }

private:
  vector<string> arglist;
};

int main(int argc, char* argv[]) {

  // Parser object for command line input.
  InputParser ip(argc, argv);

  Shower* shower;
  int seed;
  if (ip.hasOption("-seed")) {
    seed = stoi(ip.getOption("-seed"));
  }
  else {
    cout << "Provide seed" << endl;
    exit(0);
  }

  string showerName; 
  if (ip.hasOption("-shower")) {
    showerName = ip.getOption("-shower");
    if (showerName == "vetomaxdumb") {
      shower = new ShowerVetoMaxDumb(seed);
    }
    else if (showerName == "vetomax") {
      shower = new ShowerVetoMax(seed);
    }
    else if (showerName == "maxvetodumb") {
      shower = new ShowerMaxVetoDumb(seed);
    }
    else if (showerName == "maxveto") {
      shower = new ShowerMaxVeto(seed);
    }
    else if (showerName == "generateselectdumb") {
      shower = new ShowerGenerateSelectDumb(seed);
    }
    else if (showerName == "generateselect") {
      shower = new ShowerGenerateSelect(seed);
    }
    else if (showerName == "maxvetosmart") {
      shower = new ShowerMaxVetoSmart(seed);
    }
    else {
      cout << "Unknown shower" << endl;
      exit(0);
    }
  }
  else {
    cout << "Provide shower" << endl;
    exit(0);
  }
  
  int nEvents = 1e8;
  int multMax = 500;

  TH1D* counts = new TH1D("counts", "counts", multMax-1, 1.5, multMax+0.5);
  TH1D* times = new TH1D("times", "times", multMax-1, 1.5, multMax+0.5);
  TH1D* calls = new TH1D("calls", "calls", multMax-1, 1.5, multMax+0.5);
  
  double average = 0;
  for (int i=0; i<nEvents; i++) {
    if (i!=0 && i%(nEvents/100) == 0) {cout << i << endl;}

    shower->run();
    counts->Fill(shower->getMultiplicity());
    times->Fill(shower->getMultiplicity(), shower->getRunningTime());
    calls->Fill(shower->getMultiplicity(), shower->getFunctionCalls());
  }

  string fileName = "hists/" + showerName + std::to_string(seed) + ".root";
  TFile f(fileName.c_str(), "recreate");
  counts->Write();
  times->Write();
  calls->Write();

  cout << counts->GetMean() << endl;
  TH1D* timesDivide = (TH1D*)times->Clone("timesDivide");
  timesDivide->Divide(times, counts);

  TH1D* callsDivide = (TH1D*)calls->Clone("callsDivide");
  callsDivide->Divide(calls, counts);

  timesDivide->Write();
  callsDivide->Write();
}
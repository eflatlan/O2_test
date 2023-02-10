#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)


#include "HMPIDTools.h"
#include "DataFormatsHMP/Cluster.h"
#include "DataFormatsHMP/Digit.h"
#include "DataFormatsHMP/Trigger.h"
#include "HMPIDReconstruction/Clusterer.h"
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TROOT.h> // gRoot

#include <TRandom.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <fairlogger/Logger.h>
#include <fstream>
#include <vector>

#include "CommonDataFormat/InteractionRecord.h"

// C++ header files and libraries
#include <chrono>
#include <ctime>
#include <fstream>
#include <gsl/gsl>
#include <iostream>
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string.h>
#include <ctype.h>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include "randutils.hpp"
#include <random>
#endif




int range;    // = numTriggers - 1;
int numEvent; // = rng.(0, range)


using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger,
    o2::hmpid::Clusterer;
using std::vector, std::cout, std::cin, std::endl;
using std::this_thread::sleep_for;

using o2::InteractionRecord;




#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;


void changeCanvas() {

  auto folderName = (gSystem->GetWorkingDirectory());

  auto gwd = gSystem->GetWorkingDirectory();
  std::string path = static_cast<string>(gwd);

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));
  string fname = "";
  vector<string> fileInfo;


  bool fileFound = false;
  for (const auto &entry : fs::directory_iterator(path)) {
    const auto &pathName = static_cast<string>(entry.path());
    // std::cout << pathName << std::endl;


    const auto fileName = pathName.substr(pathName.find_last_of("\\/")+1, pathName.length());
    std::cout << " fileName " << fileName << std::endl;

    const auto &fileEnd =
        (fileName.substr(fileName.length() - 5));
    std::cout << " fileEnd " << fileEnd << std::endl;
    if (fileEnd != ".root")
      continue;
    
    
    auto fName = fileName;
    
    for(auto &c : fName)  c = std::tolower(static_cast<unsigned char>(c));



    auto tmp = pathName.substr(0, pathName.find_last_of("\\/"));
    auto folderName = tmp.substr(tmp.length() - 6);

     
    std::unique_ptr<TDirectory> tDir = std::make_unique<TDirectory>();
    std::unique_ptr<TTree> mTree;

    if (true/*folderName.length() == 6*/) {
      if (true/*std::all_of(folderName.begin(), folderName.end(), ::isdigit)*/) {

        std::cout << " folderName " << folderName << std::endl;
        fileFound = true;

	std::unique_ptr<TFile> file = std::make_unique<TFile>(fileName.c_str(), "READ");
        //auto& file2 = tDir->OpenFile(fileName.c_str());
	std::unique_ptr<TFile> fileOut = std::make_unique<TFile>(Form("fileName%s.c_str()", "mod"), "UPDATE");
	
	if(!file){
	  std::cout << "file Nullptr: " << endl;
	  continue;
	}
	assert(file && !file->IsZombie());

	auto fCl = file->Clone();
	auto t = file->Get("ccdb_object");
	if(t==nullptr) continue;
	auto tClone = t->Clone();
	tClone->Print(); 

	fCl->Print(); 
	if(t == nullptr) {
	  cout << "Nullptr t " << endl;
	  continue; 
	}

	auto tCan = std::make_unique<TCanvas>();	
	auto tPad = static_cast<TPad*>(tCan->cd());
	gStyle->SetOptStat("");
	tClone->Draw("Colz");
	tCan->SaveAs(Form("mod%s",fileName.c_str()));
	auto a = fileName.substr(fileName.length()-4);
	tCan->SaveAs(Form("mod%s.png",a.c_str()));
	auto a2 = fileName.substr(fileName.length()-4, fileName.length());
	tClone->SaveAs(Form("mod__%s", fileName.c_str()));

	/*
	auto tObj = fCl->Get("ccdb_object");
	tObj->Draw("Colz");
	tCan->SaveAs(Form("mod_%s",fileName.c_str()));
	auto a2 = fileName.substr(fileName.length()-4);
	tCan->SaveAs(Form("mod_%s.png",a2.c_str()));
	* /

        if(fileFound) {
          std::cout << "Found file digit " << endl;
	}*/
        // readDigits();
      } 
    }
    // cout << " Digits file << with {} Clusters, {} Triggers", ,
    // clusters.size(), clusterTriggers.size());
  }
  if (!fileFound) {
    cout << "No fitting file found!";
    return;
    std::exit(0);
  }

}

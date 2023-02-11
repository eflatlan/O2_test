#if !defined(__CLING__) || defined(__ROOTCLING__)
//#if !defined(__CINT__) || defined(__MAKECINT__)


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
float ratio;

using o2::hmpid::Cluster, o2::hmpid::Digit, o2::hmpid::Trigger,
    o2::hmpid::Clusterer;
using std::vector, std::cout, std::cin, std::endl;
using std::this_thread::sleep_for;

using o2::InteractionRecord;


vector<Digit> oneEventDigits; vector<Cluster> oneEventClusters;
vector<vector<Digit>> nestedVectorDig; vector<vector<Cluster>> nestedVectorClu;

std::vector<o2::hmpid::Digit> mDigitsFromFile,
    *mDigitsFromFilePtr = &mDigitsFromFile;
std::vector<o2::hmpid::Trigger> mTriggersFromFile,
    *mTriggersFromFilePtr = &mTriggersFromFile;


    
// void SaveFolder(string inpFile);
void changeFont();

vector<double> lastTimes, firstTimes, timeOfEvents;

bool mReadFile = false;
float mSigmaCut[7] = {4, 4, 4, 4, 4, 4, 4};

std::unique_ptr<TFile> mFile, mFileOut; ///< input file containin the tree
std::unique_ptr<TTree> mTree, mTreeOut; ///< input tree

vector<Triggers> mTriggersOut;
vector<Digits> mDigitsOut;

std::unique_ptr<Clusterer> mRec; // ef: changed to smart-pointer

void initFileIn(const std::string &fileName);
void initFileOut(const std::string &fileName);



vector<TriggerTimeInf> triggerInfoVec;*/

#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;



void readClusters(float _ratio) 
{
  ratio = _ratio;
  auto folderName = (gSystem->GetWorkingDirectory());

  auto gwd = gSystem->GetWorkingDirectory();
  std::string path = static_cast<string>(gwd);

  auto runNumber = (gwd.substr(gwd.length() - 9, gwd.length()));

  std::array<std::unique_ptr<TH1F>, 7> digCharge, hMipCharge, digPerEvent,
      digChargeSmall;
  std::array<std::unique_ptr<TH2F>, 7>  mapCharge4, digMapEntyAvg, mapEntCount, , oneEventDigMap, oneEventCluMap;

  vector<Cluster> clusters;
  vector<Trigger> clusterTriggers;
  vector<Digit> digits;
  string fname = "";
  vector<string> fileInfo;
  int numTriggers;

  bool fileFound = false;
  for (const auto &entry : fs::directory_iterator(path)) {
    const auto &pathName = static_cast<string>(entry.path());
    // std::cout << pathName << std::endl;
    if (pathName.size() <= 5)
      continue;

    const auto fileName = pathName.substr(pathName.find_last_of("\\/")+1, pathName.length());
    std::cout << " fileName " << fileName << std::endl;

    if (fileName.size() <= 5)
      continue;

    const auto &fileEnd =
        (fileName.substr(fileName.length() - 5));
    std::cout << " fileEnd " << fileEnd << std::endl;
    if (fileEnd != ".root")
      continue;
    
    
    auto fName = fileName;
    
    for(auto &c : fName)  c = std::tolower(static_cast<unsigned char>(c));

    if(fName.find("hmp") == std::string::npos || fName.find("dig") == std::string::npos)
      continue;

    std::cout << " dig2clus " << pathName << std::endl;


    auto tmp = pathName.substr(0, pathName.find_last_of("\\/"));
    auto folderName = tmp.substr(tmp.length() - 6);
    std::cout << " folderName " << folderName << std::endl;

    std::cout << " folderName.find(\\/) " << folderName.find("\\/") << std::endl;
    std::cout << " std::string::npos " << std::string::npos << std::endl;
    if(folderName.find("\\/") != std::string::npos){
      fname = folderName;
    } else{
      fname = "simulation";
    }

    if (true/*folderName.length() == 6*/) {
      if (true/*std::all_of(folderName.begin(), folderName.end(), ::isdigit)*/) {

        std::cout << " folderName " << folderName << std::endl;
        fileFound = true;

        initFileIn(fileName);
        initFileOut(fileName);

        iterate(mTriggersFromFile, mDigitsFromFile,  mTriggersOut ,mDigitsOut);

        if(fileFound) {
        std::cout << "Found file digit = " << digits.size(); }
        numTriggers = clusterTriggers.size();
        char *fn = strdup(folderName.c_str());
        std::cout << " fname " << fn << std::endl;
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



  for (int i = 0; i < 7; i++) {

    // ========== Map of One Event Digit Hits  [2D Hist 160x144] =========================
    const char *canoneEventDigMap = Form("One Event Digits %i; x [cm];y [cm]", i);
    oneEventDigMap[i].reset(
        new TH2F(canoneEventDigMap, canoneEventDigMap, 160, 0, 159, 144, 0, 143));
    oneEventDigMap[i]->SetTitleSize(oneEventDigMap[i]->GetTitleSize("x") * 1.3, "xy");

    // ========== Map of One Event Cluster Hits  [2D Hist 160x144] =========================
    const char *canOneEventClu = Form("One Event Clusters %i; x [cm];y [cm]", i);
    oneEventCluMap[i].reset(
        new TH2F(canOneEventClu, canOneEventClu, 160, 0, 159, 144, 0, 143));
    oneEventCluMap[i]->SetTitleSize(oneEventCluMap[i]->GetTitleSize("x") * 1.3, "xy");

    // ================= Average Digit Charge Map by number of Entries per Channel ( [2D histogram 144x160]) ===================
    // canvas digMapEntyAvg
    const char *canEntAvg =
        Form("Average Charge Entry Per Channel %i;  x [cm];y [cm]", i);
    digMapEntyAvg[i].reset(
        new TH2F(canEntAvg, canEntAvg, 160, 0, 159, 144, 0, 143));
    digMapEntyAvg[i]->SetTitleSize(digMapEntyAvg[i]->GetTitleSize("x") * 1.3,
                                   "xyz");


    // ================= Nonzero Digit Entries ( [2D histogram 144x160 mapEntCount]) ===================
    // canvas canvDigitCnt
    const char *canEntCnt =
        Form("Digit-Entires Per Channel %i;  x [cm];y [cm]", i);
    mapEntCount[i].reset(
        new TH2F(canEntCnt, canEntCnt, 160, 0, 159, 144, 0, 143));
    // digMap[i].reset(new TH2F(canDigMap, canDigMap, 160*0.8, 0, 159*0.8,
    // 144*0.8, 0, 160*0.8));
    mapEntCount[i]->SetTitleSize(mapEntCount[i]->GetTitleSize("x") * 1.3,
                                 "xyz");
  }


  std::unique_ptr<Trigger> pTgr;
  std::unique_ptr<Cluster> pClu, pCluEvt, cluster;
  const int digSize = digits.size();
  Printf("digit size = %i", digSize);
  int padChX = 0, padChY = 0, module = 0;

  std::array<float, 7> avgDig;


  for (const auto &dig : digits) {
    const auto& digId = dig.getPadID();

    Digit::pad2Absolute(digId, &module, &padChX, &padChY);
    const auto &charge = dig.getQ();


    if (charge == 0) {
      cout << " charge == 0 " << module << endl;
      dig0Charge[module]->Fill(padChX, padChY, 1.);
    }

    if (charge > 0) {
      chargeAvgByEntries[module][padChX][padChY] += charge;
      chargeAvgCount[module][padChX][padChY]++;

      // digMapEntyAvg[module]->Fill(padChX, padChY, charge);
      mapEntCount[module]->Fill(padChX, padChY, 1.0);
    }

    // digMap[module]->Fill(padChX, padChY, padDigOff[module][padChX][padChY]);
    // padDigits[module][padChX][padChY] += dig.getQ();
  }


  }


  // const char* runLabel = Form("%i  Duration = %s", fname, f1);
  const char *runLabel = Form("%s", folder);


  // location of Pads in Canvas to follow HMPID-modules
  static constexpr int posArr[] = {9, 8, 6, 5, 4, 2, 1};


    // ========== Map of One Event Digit Hits  [2D Hist 160x144] =========================


  const char* oneEvDigStr = Form("One Event Digit%s", folder);
  auto oneEventDigitCanvas = std::make_unique<TCanvas>(oneEvDigStr, oneEvDigStr, 1200, 1200);
  oneEventDigitCanvas->Divide(3, 3);

  for(const auto& dig : oneEventDigits){
    const auto& digId = dig.getPadID();
    //cout << "filling digits " << digId << endl;
    Digit::pad2Absolute(digId, &module, &padChX, &padChY);
    oneEventDigMap[module]->Fill(padChX, padChY, dig.getQ());
  } 

  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    TPad *pad = static_cast<TPad *>(oneEventDigitCanvas->cd(pos));
    oneEventDigMap[iCh]->SetStats(kFALSE);
    oneEventDigMap[iCh]->SetMarkerStyle(3);
    oneEventDigMap[iCh]->Draw("Colz");
  }
  oneEventDigitCanvas->SaveAs(Form("One Event Digit Map_%s_.png", folder));
  oneEventDigitCanvas->Show();


    // ========== Map of One Event Cluster Hits  [2D Hist 160x144] =========================

  const char* oneEvCluStr = Form("One Event Cluster%s", folder);
  auto oneEventClusterCanvas = std::make_unique<TCanvas>(oneEvCluStr, oneEvCluStr, 1200, 1200);
  oneEventClusterCanvas->Divide(3, 3);

  for(auto& clu : oneEventClusters){
    auto chamber = clu.ch();
    auto x = clu.x();
    auto y = clu.y();
    auto charge = clu.q();
    oneEventCluMap[chamber]->Fill(x, y, charge);
  }


  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    TPad *pad = static_cast<TPad *>(oneEventClusterCanvas->cd(pos));
    oneEventCluMap[iCh]->SetStats(kFALSE);
    oneEventCluMap[iCh]->SetMarkerStyle(3);
    oneEventCluMap[iCh]->Draw("Colz");//    oneEventCluMap[iCh]->DrawCopy("Colz");
  }
  oneEventClusterCanvas->SaveAs(Form("One Event Cluster Map_%s_.png", folder));
  oneEventClusterCanvas->Show();


  // ================= Average Digit Charge Map by number of Entries per Channel (canvDigitAvgEvt [2D histogram 144x160]) ===================
  std::unique_ptr<TCanvas> canvDigitAvgEvt, canvDigitCnt;
  canvDigitAvgEvt.reset(new TCanvas(Form("Event Channel Average %s", folder),
                          Form("Event Channel Average %s", folder), 1200, 2000));
  canvDigitAvgEvt->Divide(3, 3);
  canvDigitAvgEvt->cd(3);

  // ================= Nonzero Digit Entries (canvDigitCnt [2D histogram 144x160 mapEntCount]) ===================
  canvDigitCnt.reset(new TCanvas(Form("Event Channel Count %s", folder),
                          Form("Event Channel Count %s", folder), 1200, 2000));
  canvDigitCnt->Divide(3, 3);
  canvDigitCnt->cd(3);

  int maxCnt[7] = {0};
  for (int iCh = 0; iCh < 7; iCh++) {
    const auto &pos = posArr[iCh];
    TPad *pad = static_cast<TPad *>(canvDigitAvgEvt->cd(pos));
    for (int x = 0; x < 160; x++) {
      for (int y = 0; y < 144; y++) {
        if (chargeAvgCount[iCh][x][y] > maxCnt[iCh]) {
          maxCnt[iCh] = chargeAvgCount[iCh][x][y];
        }
        auto chAverage =
            chargeAvgByEntries[iCh][x][y] / chargeAvgCount[iCh][x][y];
        if (chargeAvgCount[iCh][x][y] != 0) {
          digMapEntyAvg[iCh]->Fill(x, y, chAverage);
        } else {
          digMapEntyAvg[iCh]->Fill(x, y, 0.);
        }
      }
    }
    cout << "maxCnt " << iCh << maxCnt[iCh] << endl;

    pad->SetLeftMargin(+.035 + pad->GetLeftMargin());
    pad->SetBottomMargin(.085 + pad->GetBottomMargin());
    pad->SetRightMargin(.085 + pad->GetRightMargin());
    pad->SetLogz(1);
    digMapEntyAvg[iCh]->SetStats(kFALSE);
    digMapEntyAvg[iCh]->Draw("Colz");

    TPad *pad2 = static_cast<TPad *>(canvDigitCnt->cd(pos));
    pad2->SetLogz(1);
    pad2->SetLeftMargin(pad->GetLeftMargin());
    pad2->SetBottomMargin(pad->GetBottomMargin());
    pad2->SetRightMargin(pad->GetRightMargin());
    mapEntCount[iCh]->SetStats(kFALSE);
    mapEntCount[iCh]->Draw("Colz");
  }
  canvDigitAvgEvt->SaveAs(Form("AverageByEntries_%s_.png", folder));
  canvDigitCnt->SaveAs(Form("CountEntries_%s_.png", folder));




  sleep_for(5000ms);
  // to edit TCanvases::
  bool userInput = false;
  while(!userInput){
    sleep_for(5000ms);   
    string uInputString;
    sleep_for(50ms);
    cin >> uInputString;
    if(uInputString == "C"){
      sleep_for(10000ms);

      while(uInputString != "Q"){
        cin >> uInputString;
        sleep_for(10000ms);
      }
    }
    sleep_for(1000ms);
  }

}

void initFileOut(const std::string &filename) {

  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  // Create the TFIle
  mFileOut = std::make_unique<TFile>(Form("out_%f_%s",ratio,filename.c_str()), "RECREATE");
  assert(mFileOut && !mFileOut->IsZombie());

  tit = TString::Format("HMPID Clusters File Decoding");
  mTreeOut = std::make_unique<TTree>("o2hmp", tit);
  mTreeOut->Branch("InteractionRecords", &mTriggersOut);
  mTreeOut->Branch("HMPIDDigits", &mDigitsOut);

  if (!mTreeOut) {
    LOG(warn)
        << "HMPID DigitToClusterSpec::init() : Did not find o2sim tree in "
        << filename.c_str() << endl;
    return;
    std::exit(0);
  }

  mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
}



void initFileIn(const std::string &filename) {

  long mDigitsReceived, mClustersReceived, mTriggersReceived = 0;
  // Create the TFIle
  mFile = std::make_unique<TFile>(filename.c_str(), "OLD");
  assert(mFile && !mFile->IsZombie());

  mTree.reset((TTree *)mFile->Get("o2sim"));
  if (!mTree) {
    mTree.reset((TTree *)mFile->Get("o2hmp"));
  }

  if (!mTree) {

    LOG(warn)
        << "HMPID DigitToClusterSpec::init() : Did not find o2sim tree in "
        << filename.c_str() << endl;
    return;
    std::exit(0);
  }

  if ((mTree->GetBranchStatus("HMPDigit")) == 1) {
    mTree->SetBranchAddress("HMPDigit", &mDigitsFromFilePtr);
  } else if ((mTree->GetBranchStatus("HMPIDDigits")) == 1) {
    mTree->SetBranchAddress("HMPIDDigits", &mDigitsFromFilePtr);
  } else {
    LOG(warn) << "HMPID DigitToClusterSpec::init() : Error in branches!"
              << endl;
    return;
    std::exit(0);
  }

  mTree->SetBranchAddress("InteractionRecords", &mTriggersFromFilePtr);
  mTree->Print("toponly");
}


void iterate(vector<Trigger>& triggersIn, vector<Digit>& digitsIn, vector<Trigger>& triggersOut, vector<Digit>& digitsOut)
{

  std::random_device rd;
  std::uniform_real_distribution<> range(0, 1.0);

  for(const auto& trig : triggersIn){
    std::mt19937 mt(rd());
    const int& firstTrigger = trig.getFirstEntry();
    const int& lastTrigger = trig.getLastEntry();
    const int& startIndexDigit = digitsOut.size();
    int numberOfDigits = 0;
    if (mt > ratio) {
      for(int j = firstTrigger; j < lastTrigger, j++){
        digitsOut.push_back(digitsIn[j]);
        numberOfDigits++;
      }
      mTriggers.push_back(trig);
      mTriggers.back().setDataRange(startIndexDigit, numberOfDigits);
    }
}

numEvent = static_cast<int>(range(mt));
}
void changeFont() {
  /*
  std::unique_ptr<TStyle> mStyle;
  mStyle.reset(new TStyle("canvasStyle", "Canvas Root Styles"));
  */
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.925);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.25);
  gStyle->SetStatFontSize(0.065);
  gStyle->SetLegendTextSize(0.08); //
  gStyle->SetTitleSize(.055, "xzy");
  gStyle->SetTitleOffset(.925, "xz"); //.95
  gStyle->SetTitleOffset(1, "y");     // 1.1
  gStyle->SetTitleFontSize(.05);
  // gStyle->SetTitleFont(16, "xz");
  gStyle->SetLabelOffset(0.0065, "y");
  gStyle->SetLabelFont(22, "xyz");
  gStyle->SetLabelSize(.055, "xyz"); //.0525 // verdi av akser
  // mStyle->SetStyle("canvasStyle");
}
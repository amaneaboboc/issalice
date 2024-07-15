//este modificat cum trebuie
#include "TrackN.h"
#include "EventN.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TMath.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TCut.h>
#include <TLegend.h>
#include <TSpline.h>

#include "ParticleFilter.h"
#include "HistogramManager.h"

//#include "./treeClass/Include.h"

#include <iostream>
#include <fstream>
#include <map>

using namespace std;

Double_t GetPhi(Double_t fX, Double_t fY)
{
    return fX == 0.0 && fY == 0.0 ? 0.0 : TMath::ATan2(fY, fX);
}

Double_t GetEta(Double_t fX, Double_t fY, Double_t fZ)
{
    Double_t sq_roots = TMath::Sqrt(fX*fX + fY*fY + fZ*fZ);
    if(sq_roots==0){return 0.0;}
    return 0.5*TMath::Log( (sq_roots + fZ) / (sq_roots - fZ));
}

Double_t GetEnergy(Double_t m, Double_t px, Double_t py, Double_t pz)
{
return TMath::Sqrt(TMath::Power(m,2)+ TMath::Power(px,2)+ TMath::Power(py,2)+TMath::Power(pz,2));
}


TChain* CreateChainLocal(Int_t nFilesMax, const Char_t* filename, const Char_t* treeName);

void getMult(const Char_t* inFileName, const Char_t* outFileName, Bool_t useChain)
{
    
    TTree* inTree = 0;
    if (useChain){
        cout<<"Running analysis on a chain of files"<<endl;
        inTree = CreateChainLocal(0, inFileName, "teposevent");
        if (!inTree){
            cout << "Chain: " << inFileName << " does not exist!" << endl;
            return;
        }
    } else{
        cout<<"Running analysis on a single file"<<endl;
        TFile* inFile = TFile::Open(inFileName);
        if(!inFile) {
            cout << "File: " << inFileName << " does not exist!" << endl;
            return;
        }
        inTree = (TTree*)inFile->Get("Tampt");
    }

   
    Int_t centr = 10;

    Int_t nParticleFilters = 6;

    ParticleFilter particlefilter[nParticleFilters];

    particlefilter[0].addCondition(2,321,321,322); //conditie K+
    
    particlefilter[0].setName("K+");
    
    particlefilter[0].setTitle("K+");
   
    particlefilter[1].addCondition(2,-321,-321,-322); //conditie K+
    
    particlefilter[1].setName("K-");
    
    particlefilter[1].setTitle("K-");
  
    particlefilter[2].addCondition(2,211,211,212); //conditie pi+
    
    particlefilter[2].setName("pi+");
    
    particlefilter[2].setTitle("pi+");
    
    particlefilter[3].addCondition(2,-211,-211,-212); //conditie pi+
    
    particlefilter[3].setName("pi-");
    
    particlefilter[3].setTitle("pi-");
 
    particlefilter[4].addCondition(2,2212,2212,2213); //conditie p
    
    particlefilter[4].setName("p");
    
    particlefilter[4].setTitle("p");
   
    particlefilter[5].addCondition(2,-2212,-2212,-2213); //conditie pbar
    
    particlefilter[5].setName("pBar");
    
    particlefilter[5].setTitle("pBar");
    
    TH1I* h_events[centr];
    
    for(Int_t i = 0 ; i< centr; i++){
    
    h_events[i] = new TH1I(Form("h_events_%d",i),"; counts; No.",1,0,1);
    }

    EventN* eventData = 0;
    inTree->SetBranchAddress("event", &eventData);
    
    TClonesArray* trackArray = 0;
    inTree->SetBranchAddress("particle", &trackArray);
      
    Int_t nEv = inTree->GetEntries();
    cout << "Number of events: " << nEv << endl;

}

//__________________________________________________________
TChain* CreateChainLocal(Int_t nFilesMax, const Char_t* filename, const Char_t* treeName)
{
    
    TChain* chain = new TChain(treeName);
    
    // Open the input stream
    ifstream in;
    in.open(filename);
    
    Int_t nFiles = 0;
    
    // Read the input list of files and add them to the chain
    TString file;
    while(in.good() && (nFiles<nFilesMax || nFilesMax==0)) {
        in >> file;
        if (!file.Contains("root"))
            continue; // protection
        
        nFiles++;
        chain->Add(file.Data());
    }
    
    in.close();
    
    return chain;
}
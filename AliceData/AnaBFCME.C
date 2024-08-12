#include "Event.h"
#include "TrackN.h"

//#include "treeLocal/EventQcut.h"
//#include "treeLocal/TrackQcut.h"

#include "ParticleFilter.h"
#include "HistogramManager.h"

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
#include <TSpline.h>
#include <THnSparse.h>
#include <TRandom.h>

#include <iostream>
#include <fstream>
using namespace std;


TChain* CreateChainLocal(Int_t nFilesMax, const Char_t* filename, const Char_t* treeName);

void AnaBFESE(const Char_t* inFileName, const Char_t* outFileName, const Char_t* inCorrMaps, Bool_t useChain)
{

   TFile* inSp = TFile::Open(inCorrMaps);
    if (!inSp){
        cout<<"Calibration file for spline does not exist!"<< endl;
        return;
    }

    TTree* inTree = 0;
    if (useChain){
        cout<<"Running analysis on a chain of files"<<endl;
        inTree = CreateChainLocal(0, inFileName, "tree");
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
        inTree = (TTree*)inFile->Get("tree");
    }
    
      Int_t nParticleFilters = 2;

    ParticleFilter particlefilter[nParticleFilters];

    particlefilter[0].addCondition(1,2,321,322); //conditie K+
    
    particlefilter[0].setName("n+");
    
    particlefilter[0].setTitle("n+");

    particlefilter[1].addCondition(1,3,321,322); //conditie K+
    
    particlefilter[1].setName("n-");
    
    particlefilter[1].setTitle("n-");

    TH1I* h_events = new TH1I("h_events","h_events",1,0,1);

    HistogramManager histogramManager;

    std::vector<Int_t> particleCuts;
    
    for(int i = 0; i < nParticleFilters; i++){  
    
    	TString name(particlefilter[i].getName()); //particlefilter[0].getName()
    	
    	ParticleSingleHistos* single_histos = new ParticleSingleHistos(name); 

        //ParticleSingleHistos single_histos(name);
        //new(&single_histos[j]) ParticleSingleHistos(name);
        
        single_histos->createHistograms();

        histogramManager.addHistoInSet(i, single_histos);
	    
    }
    
       for(int j = 0; j < nParticleFilters; j++){
    
        for(int k = 0; k < nParticleFilters; k++){
        

            TString pairname = particlefilter[j].getName()+particlefilter[k].getName();
    
            ParticlePairHistos* pairHisto = new ParticlePairHistos(pairname);
            pairHisto->createHistograms();
//        PairDerivedHistos* derivedHisto = new PairDerivedHistos(Form("name%d",j* nParticleFilters + k));
            PairDerivedHistos* derivedHisto = new PairDerivedHistos(pairname);
            derivedHisto->createHistograms();
    
            histogramManager.addHistoInSet(j* nParticleFilters + k,pairHisto);
            histogramManager.addHistoInSet(j* nParticleFilters + k,derivedHisto);        
    	}
    
    }

 Event* eventData = 0;

    inTree->SetBranchAddress("event", &eventData);
    
    TClonesArray* trackArray = 0;
    inTree->SetBranchAddress("track", &trackArray);
    
    const Int_t nEvents = inTree->GetEntries();
    cout << "Number of events: " << nEvents << endl;
    

    Int_t currentRun = 0;
    
    for(Int_t n = 0; n < 1000; n++) {
    
     inTree->GetEntry(n);
        
        if(eventData->run != currentRun){
            cout << "New run: " << eventData->run << endl;
            currentRun = eventData->run;
            //continue;
        }
    
        if((n+1)%1000 == 0)
            cout << "Analysis - Event: " << n+1 << "/" << nEvents << endl;

             Int_t nTrk = 0;
        
        particleCuts.clear();    

        Int_t nTracks = trackArray->GetEntries();

        for(Int_t jj = 0; jj < nTracks; jj++) {
	
            TrackN* trk = (TrackN*)trackArray->At(jj);
	
            if (!trk){
                delete trk;
                continue;
            }
        double pt = trk->pt;

        if(trk->q == 1){ trk->pt = trk->pt * (1.0f / spEffPos[binQ2C]->Eval(trk->pt));}
        if(trk->q == -1){ trk->pt = trk->pt * (1.0f/spEffNeg[binQ2C]->Eval(trk->pt));}

        double trkEta = trk->eta;

        if (pt < 0.2 || pt > 2.0) continue; 

        if(trkEta < -1.0 || trkEta > 1.0) continue;

        if (trk->q == 0) continue;

        particleCuts.push_back(jj);

        }

        h_events->Fill(0);

      for(Int_t k = 0; k < particleCuts.size(); k++) {
                          
	      TrackN* trkl = (TrackN*)trackArray->At(particleCuts[k]);
            
	      for(int i = 0; i < nParticleFilters; i++){
            
            if(particlefilter[i].accept(trkl)){

        histogramManager.HistogramSet[i]->fill(trkl,1.0);
            }
	      }
    }

 for(Int_t jl = 0; jl < particleCuts.size(); jl++){

 	    TrackN* trk1 = (TrackN*)trackArray->At(particleCuts[jl]);
 
 	    for(int ii = 0; ii < nParticleFilters; ii++){
 
 	    if(particlefilter[ii].accept(trk1)){
 

            for(Int_t jk = 0; jk < particleCuts.size(); jk++) {

            TrackN* trk2 = (TrackN*)trackArray->At(particleCuts[jk]);
            if(jl == jk)continue;
            
            for(int ij = 0; ij < nParticleFilters; ij++){
            
            if (particlefilter[ij].accept(trk2)){
            
            histogramManager.HistogramPairSet[ii*nParticleFilters + ij]->fill(trk1,trk2,1.0);

	                }
	            }
	        }//pair2
	    
            }
            }

        }//pair1


    }

TFile* fOut = new TFile(outFileName,"RECREATE");

histogramManager.Write();

h_events->Write();

fOut->Close();

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
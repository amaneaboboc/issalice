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

//#include "histogramming.h"
//#include "routine.h"
//#include "Filter.h"
#include "ParticleFilter.h"
#include "HistogramManager.h"
//#include "Configuration.h"
//#include "BalanceFunctionCalculator.h"

#include <iostream>
#include <fstream>

using namespace std;

TChain* CreateChainLocal(Int_t nFilesMax, const Char_t* filename, const Char_t* treeName);


void AnaBFMult(const Char_t* inFileName, const Char_t* outFileName, const Char_t* inFileSpName, Bool_t useChain)
{


    TFile* inSp = TFile::Open(inFileSpName);
    if (!inSp){
        cout<<"Calibration file for spline does not exist!"<< endl;
        return;
    }
    
    TSpline3* sp3MultV0A = (TSpline3*)inSp->Get("spMultV0A");
    
    TTree* inTree = 0;
    if (useChain){
        cout<<"Running analysis on a chain of files"<<endl;
        inTree = CreateChainLocal(0, inFileName, "Tree");
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
        inTree = (TTree*)inFile->Get("Tree");
    }

    Int_t centr = 10;
    
    TH1I* h_events[centr];
    
    for(Int_t i = 0; i < centr;i++){
    h_events[i] = new TH1I(Form("h_events_%d",i),"; events; Counts",1,0,1);
    }

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
    
    
    HistogramManager histogramManager;
      
    std::vector<Int_t> particleCuts;
    
    for(int i = 0; i < nParticleFilters; i++){  
    
    	for(int j = 0; j < centr; j++){
    	
    	TString name(particlefilter[i].getName()); //particlefilter[0].getName()
   	name.Append("_"+to_string(j));
    	
    	ParticleSingleHistos* single_histos = new ParticleSingleHistos(name); 

        //ParticleSingleHistos single_histos(name);
        //new(&single_histos[j]) ParticleSingleHistos(name);
        
        single_histos->createHistograms();

        histogramManager.addHistoInSet(i*centr+j, single_histos);
	}
    }
    
  
        //pair and derived histo initialization
    

 /*   
    for(int j = 0; j < nParticleFilters; j++){
    
    for(int k = 0; k < nParticleFilters; k++){
    
   // for(int l = 0; l < centr; l++){
    
    TString pairname = particlefilter[j].getName()+particlefilter[k].getName();

    //pairname.Append("_"+to_string(l));
    
        ParticlePairHistos* pairHistos = new ParticlePairHistos(pairname);
        
        pairHistos->createHistograms();
//        PairDerivedHistos* derivedHistos = new PairDerivedHistos(Form("name%d",j* nParticleFilters + k));
	//PairDerivedHistos* derivedHistos = new PairDerivedHistos(particlefilter[j].getName() + particlefilter[k].getName()+Form("_%d",l));
        //derivedHistos->createHistograms();
    
        histogramManager.addHistoInSet(j*nParticleFilters+ k  ,pairHistos);
        //histogramManager.addHistoInSet(l*centr + j* nParticleFilters + k,derivedHistos);
 		
 	    }       
    	}
    */
    //}
    
        //pair and derived histo initialization
    
    for(int j = 0; j < nParticleFilters; j++){
    
    for(int k = 0; k < nParticleFilters; k++){
    
    for(int l = 0; l < centr; l++){
    
    TString pairname = particlefilter[j].getName()+particlefilter[k].getName();
    pairname.Append("_"+to_string(l));
    
        ParticlePairHistos* pairHisto = new ParticlePairHistos(pairname);
        pairHisto->createHistograms();
//        PairDerivedHistos* derivedHisto = new PairDerivedHistos(Form("name%d",j* nParticleFilters + k));
PairDerivedHistos* derivedHisto = new PairDerivedHistos(pairname);
        derivedHisto->createHistograms();
    
        histogramManager.addHistoInSet(centr*(j* nParticleFilters + k)+l,pairHisto);
        histogramManager.addHistoInSet(centr*(j* nParticleFilters + k)+l,derivedHisto);
        
            }
        
    	}
    
    }

    EventN* eventData = 0;
    inTree->SetBranchAddress("event", &eventData);
    
    TClonesArray* trackArray = 0;
    inTree->SetBranchAddress("particle", &trackArray);
      
    Int_t nEv = inTree->GetEntries();
    cout << "Number of events: " << nEv << endl;
        

        for(Int_t n = 0; n < nEv; n++) { //nEv
        
        inTree->GetEntry(n);
        
        if((n+1)%100000  == 0)
            cout << "Analysis - Event: " << n+1 << "/" << nEv << endl;
  
        Int_t nTracks = trackArray->GetEntries();
        
              //cout<<"new Evt"<<endl;
              //cout<<"nr of tracks "<<nTracks<<endl;
              
            particleCuts.clear();    
  	    
  	    Double_t multch2 = 0;
  	         
            for(Int_t jj = 0; jj < nTracks; jj++) {    

            TrackN* trk = (TrackN*)trackArray->At(jj);
 
            //if (i == jj-1)   continue;
            
            if (!trk){
                delete trk;
                continue;
            }

            double pt = TMath::Sqrt(trk->px*trk->px + trk->py*trk->py);
            
            int charge = trk->q;
            
            double phiNew = trk->phi + TMath::Pi();
            double trkEta = trk->eta;
            int pdg = trk->pdgCode;
            double trkY = 0.5 * log((trk->e + trk->pz)/(trk->e - trk->pz));
          
            if(charge == 0) continue;
            
            if(pt > 0.2 && pt < 3.0 && trkY>3.0 && trkY<5.0) {
              multch2++;
              }
            
            if (pt < 0.2 || pt > 2.0) continue;
            
            //if(trkEta < -1.0 || trkEta > 1.0) continue; //!fillY &&
            
            if (trkY < -5.0 || trkY > 5.0) continue;
              
            if(TMath::Abs(pdg) < 90) continue;
           
            particleCuts.push_back(jj);
            
	    }//track

            if(multch2 < 2) continue;
            
        Double_t percMultV0A = 100.*(1. - sp3MultV0A->Eval(multch2));

        if(percMultV0A < 0){
            cout<<"Negative percentile: put it to 0"<<endl;
            percMultV0A = 0;
        }
       // cout<<"val spline: "<<percMultV0A<<endl;
        
        //percentile_dist->Fill(percMultV0A);
        //fill percentile distribution

        Short_t multp = -1;
        if (percMultV0A < 10.)
            multp = 0;
        else if (percMultV0A >= 10. && percMultV0A < 20.)
            multp = 1;
        else if (percMultV0A >= 20. && percMultV0A < 30.)
            multp = 2;
        else if (percMultV0A >= 30. && percMultV0A < 40.)
            multp = 3;
        else if (percMultV0A >= 40. && percMultV0A < 50.)
            multp = 4;
        else if (percMultV0A >= 50. && percMultV0A < 60.)
            multp = 5;
        else if (percMultV0A >= 60. && percMultV0A < 70.)
            multp = 6;
        else if (percMultV0A >= 70. && percMultV0A < 80.)
            multp = 7;
        else if (percMultV0A >= 80. && percMultV0A < 90.)
            multp = 8;
        else if (percMultV0A >= 90.)
            multp = 9;
            
            //histogramManager.HistogramSet[multp]->fillEvents(1);

	    h_events[multp]->Fill(0);
	    
	    for(Int_t k = 0; k < particleCuts.size(); k++) {
                          
	    TrackN* trkl = (TrackN*)trackArray->At(particleCuts[k]);
            
	    for(int i = 0; i < nParticleFilters; i++){
 
            if (particlefilter[i].accept(trkl)){
            
            histogramManager.HistogramSet[i*centr+multp]->fill(trkl,1.0);
            

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
            
            histogramManager.HistogramPairSet[centr * (ii*nParticleFilters + ij) + multp]->fill(trk1,trk2,1.0);

	    }
	    }
	   }//pair2
	    
}
}

}//pair1
       
      
}//event


//-------------- single particle derived histos (n1n1) , particle pair derived histos (R2, shft)

//scaling by number of events

for(Int_t i = 0; i< nParticleFilters; i++){

for(Int_t k = 0; k<centr; k++){
Int_t evts = histogramManager.HistogramSet[k]->GetEvents();
//histogramManager.HistogramSet[i*centr+k]->Scale(evts);
}

for(Int_t j = 0;j <nParticleFilters; j++){

for(Int_t k = 0; k<centr; k++){
Int_t evts = histogramManager.HistogramSet[k]->GetEvents();
//histogramManager.HistogramPairSet[centr * (i*nParticleFilters + j) + k]->Scale(evts);
}

}

}


std::vector<ParticleSingleHistos*> vect = histogramManager.get_sgHist();
std::vector<ParticlePairHistos*> vect_pair = histogramManager.get_pairHist();

for( Int_t ii = 0; ii < nParticleFilters; ii++){

    for(Int_t ij = 0; ij < nParticleFilters; ij++){
	
	for(Int_t k = 0; k<centr; k++){		
	
    histogramManager.HistogramDerivedSet[centr*(ii*nParticleFilters + ij) + k]->calculatePairDerivedHistograms(vect[ii*centr+k],vect[ij*centr+k],vect_pair[centr * (ii*nParticleFilters + ij) + k],1.0);
    
   // histogramManager.HistogramDerivedSet[ii*nParticleFilters + ij]->calculate_BalFct(particlefilter[ii].getName(),particlefilter[ij].getName());//,vect_pair[ij*nParticleFilters+ii]
	}
    }

}





TString strExt(outFileName);

//strExt.Append("_"+to_string(centr));
strExt.Append(".root");

cout<<"File output name: "<< strExt <<endl;
cout<<"---------- End of Analysis ----------"<<endl;

TFile* outfile =  new TFile(strExt,"RECREATE");
             
histogramManager.Write();

for(Int_t j = 0; j<centr;j++){
h_events[j]->Write();
}

outfile->Close();

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

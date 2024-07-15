//este modificat cum trebuie
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

#include "Track.h"
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

void Ana(const Char_t* inFileName, const Char_t* outFileName, const Char_t* inFileSpName, Bool_t useChain)
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
        inTree = (TTree*)inFile->Get("teposevent");
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
    
    HistogramManager histogramManager;
    
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
    Int_t np=0;
    
    Float_t px[5000], py[5000], pz[5000], e[5000];
    Int_t id[5000], np;
    
    inTree->SetBranchAddress("px",&px);
    inTree->SetBranchAddress("py",&py);
    inTree->SetBranchAddress("pz",&pz);
    inTree->SetBranchAddress("e",&e);
    inTree->SetBranchAddress("id",&id);
    //inTree->SetBranchAddress("ist",&ist);
    inTree->SetBranchAddress("np",&np);

    Int_t nEv = inTree->GetEntries();
    cout << "Number of events: " << nEv << endl;
         
         
    for(Int_t n = 0; n < nEv; n++) { //nEv
        
            if((n+1)%100000  == 0)
            cout << "Analysis - Event: " << n+1 << "/" << nEv << endl;
        
        
        inTree->GetEntry(n);
     
         //TClonesArray* fPartArray = new TClonesArray("TrackN", 5000);

         Track track[5000];
         
         Int_t nAddTrk=0;
         
     	for(Int_t jj = 0; jj < np; jj++) {
     	
    	    if(ist[jj] != 0) continue;
    	    
     	    double pt = TMath::Sqrt(px[jj]*px[jj] + py[jj]*py[jj]);
            
            int charge = m[];
            
            double phiNew = GetPhi(px[jj],py[jj]);//+ 3.1415;
            double trkEta = GetEta(px[jj],py[jj],pz[jj]);
            int pdg = id[jj];
            Double_t trkY = 0.5 * TMath::Log((GetEnergy(e[jj],px[jj],py[jj],pz[jj]) + pz[jj])/(GetEnergy(e[jj],px[jj],py[jj],pz[jj]) - pz[jj]));
     	    //cout<<"trkEta "<<trkEta<<" trkY "<<trkY<<endl;

            if(trkY > -5 && trkY < -3)
            {
                nAddTrk++;
            } //multiplicity estimator

            if(charge == 0) continue;
            
            if (pt < 0.2 || pt > 2.0) continue;
            
            //if(trkEta < -1.0 || trkEta > 1.0) continue; //!fillY &&
            
            if (trkY < -1.0 || trkY > 1.0) continue;
     
     	     //TrackN* track = new((*fPartArray)[nAddTrk]) TrackN();

     		
     	        track[nAddTrk].pdgCode = pdg;
     	        track[nAddTrk].phi = phiNew;
     	        track[nAddTrk].status = ist[jj];
                track[nAddTrk].px = px[jj];
                track[nAddTrk].py = py[jj];
                track[nAddTrk].pz = pz[jj];
                track[nAddTrk].e = GetEnergy(e[jj],px[jj],py[jj],pz[jj]);
                track[nAddTrk].phi = phiNew;
                track[nAddTrk].eta = trkEta;
                track[nAddTrk].q = charge;
                
     	        /*
                track->status = ist[jj];
                track->px = px[jj];
                track->py = py[jj];
                track->pz = pz[jj];
                track->e = GetEnergy(e[jj],px[jj],py[jj],pz[jj]);
                track->phi = phiNew;
                track->eta = trkEta;
                track->q = charge;
     	*/
     	}
     	
     	if(nAddTrk < 2) continue;
        //cout<<"nr trackuri dsa:"<< nAddTrk <<endl;
        Double_t percMultV0A = 100.*(1. - sp3MultV0A->Eval(nAddTrk));

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

     	h_events[multp]->Fill(0);
     
    //cout<<"nr trackuri:"<< fPartArray->GetEntries()<<endl;
  
    Int_t nTracks = nAddTrk;//fPartArray->GetEntries();
  
    for(Int_t j = 0; j < nTracks; j++){
    
    //TrackN* trk = (TrackN*)fPartArray->At(j);
     	   
     	   //TrackN trk = track[j];

                  
	    for(int i = 0; i < nParticleFilters; i++){
 
            if (particlefilter[i].accept(track[j])){
            
            cout<<"track filt: "<<track[j].pdgCode<<endl;
            histogramManager.HistogramSet[i*centr+multp]->fill(track[j],1.0);
            
           }            
	}    
 
    }

/*
    
    	    for(Int_t jl = 0; jl < nTracks; jl++){

 	   // TrackN* trk1 = (TrackN*)fPartArray->At(jl);
 
 	    for(int ii = 0; ii < nParticleFilters; ii++){
 
 	    if(particlefilter[ii].accept(trk1)){
 

            for(Int_t jk = 0; jk < nTracks; jk++) {

            //TrackN* trk2 = (TrackN*)fPartArray->At(jk);
            
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
    */
    
    //fPartArray->Clear("C");
        
    }// event loop
    
    
    TCanvas* c = new TCanvas();


TString strExt(outFileName);

//strExt.Append("_"+to_string(centr));
strExt.Append(".root");

cout<<"File output name: "<< strExt <<endl;
cout<<"---------- End of Analysis ----------"<<endl;

TFile* outfile =  new TFile(strExt,"RECREATE");
             
histogramManager.Write();          

for(Int_t i = 0; i< centr;i++){
h_events[i]->Write();
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

/*
In root terminal:
   .L libMyTree.so
   .L doAnaFlow.C+
   doAnaFlow(2, 0, "chain.txt", "histAnaFlowHM.root", 1)

pentru chain in directorul newCR:
find `pwd` | grep treePythiaNewCRPP.root > chain.txt
*/


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


#include <iostream>
#include <fstream>
using namespace std;


TChain* CreateChainLocal(Int_t nFilesMax, const Char_t* filename, const Char_t* treeName);


void doAnaFlowQCgap(Float_t nHarm, const Char_t* inFileName, const Char_t* outFileName, Bool_t isMB, Bool_t isLowM, Bool_t useChain)
{
    
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
    
 
    EventN* eventData = 0;
    inTree->SetBranchAddress("event", &eventData);
    
    TClonesArray* trackArray = 0;
    inTree->SetBranchAddress("particle", &trackArray);

        
    Int_t nEv = inTree->GetEntries();
    cout << "Number of events: " << nEv << endl;
    
   
    const Int_t nPt = 25;
    Double_t binsPt[nPt+1] = {0., 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10., 13., 16., 20};
    
    const Int_t nPtV0Ca = 14;
    Double_t binsPtV0Ca[nPt+1] = {0., 0.2, 0.5, 0.8, 1.1, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10., 15., 20};
    

    TH2D* fEtaPhiAll = new TH2D("fEtaPhiAll", "; #eta; #phi", 88, -1.1, 1.1, 72, 0, 2*TMath::Pi());
    TH1D* fPtAll = new TH1D("fPtAll","; p_{T} (GeV/c); Counts", nPt, binsPt);
    TH1I* fMult = new TH1I("fMult", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    //pi+K0+lambda+xi+omega
    TH1D* fPtOm = new TH1D("fPtOm","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiOm = new TH2D("fEtaPhiOm","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultOm = new TH1I("fMultOm", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    TH1D* fPtXi = new TH1D("fPtXi","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiXi = new TH2D("fEtaPhiXi","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultXi = new TH1I("fMultXi", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    TH1D* fPtL = new TH1D("fPtL","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiL = new TH2D("fEtaPhiL","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultL = new TH1I("fMultL", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    TH1D* fPtK0 = new TH1D("fPtK0","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiK0 = new TH2D("fEtaPhiK0","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultK0 = new TH1I("fMultK0", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    TH1D* fPtPi = new TH1D("fPtPi","; p_{T} (GeV/c); Counts", nPt, binsPt);
    TH2D* fEtaPhiPi = new TH2D("fEtaPhiPi","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultPi = new TH1I("fMultPi", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    
    TProfile* fSinAll = new TProfile("fSinAll", "; multiplicity; sin(n*#phi)", 600, 0, 600);
    TProfile* fCosAll = new TProfile("fCosAll", "; multiplicity; cos(n*#phi)", 600, 0, 600);

    TProfile* fQC2AllQ = new TProfile("fQC2AllQ", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC4AllQ = new TProfile("fQC4AllQ", "; multiplicity; QC{4}", 600, 0, 600);
    TProfile* fQC6AllQ = new TProfile("fQC6AllQ", "; multiplicity; QC{6}", 600, 0, 600);
    
    TProfile* fQC2AllQGap1 = new TProfile("fQC2AllQGap1", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap2 = new TProfile("fQC2AllQGap2", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap3 = new TProfile("fQC2AllQGap3", "; multiplicity; QC{2}", 600, 0, 600);
    
    TProfile* fQC4AllQGap = new TProfile("fQC4AllQGap", "; multiplicity; QC{4}", 600, 0, 600);
    TProfile* fQC2AllQGap12 = new TProfile("fQC2AllQGap12", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap13 = new TProfile("fQC2AllQGap13", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap14 = new TProfile("fQC2AllQGap14", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap23 = new TProfile("fQC2AllQGap23", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap24 = new TProfile("fQC2AllQGap24", "; multiplicity; QC{2}", 600, 0, 600);
    TProfile* fQC2AllQGap34 = new TProfile("fQC2AllQGap34", "; multiplicity; QC{2}", 600, 0, 600);
    
    
    TProfile* fVnAllg1A = new TProfile("fVnAllg1A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnAllg1C = new TProfile("fVnAllg1C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
        
    TProfile* fVnPig1A = new TProfile("fVnPig1A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnPig1C = new TProfile("fVnPig1C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
        
    TProfile* fVnKg1A = new TProfile("fVnKg1A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnKg1C = new TProfile("fVnKg1C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    
    TProfile* fVnPg1A = new TProfile("fVnPg1A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnPg1C = new TProfile("fVnPg1C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    
    TProfile* fVnLg1A = new TProfile("fVnLg1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnLg1C = new TProfile("fVnLg1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fVnK0g1A = new TProfile("fVnK0g1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnK0g1C = new TProfile("fVnK0g1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fVnXig1A = new TProfile("fVnXig1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnXig1C = new TProfile("fVnXig1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fVnOmg1A = new TProfile("fVnOmg1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnOmg1C = new TProfile("fVnOmg1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    
    
    TProfile* fVnAllg2A = new TProfile("fVnAllg2A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnAllg2C = new TProfile("fVnAllg2C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
        
    TProfile* fVnPig2A = new TProfile("fVnPig2A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnPig2C = new TProfile("fVnPig2C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
        
    TProfile* fVnKg2A = new TProfile("fVnKg2A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnKg2C = new TProfile("fVnKg2C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    
    TProfile* fVnPg2A = new TProfile("fVnPg2A", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    TProfile* fVnPg2C = new TProfile("fVnPg2C", "; p_{T} (GeV/c); v_{n}", nPt, binsPt);
    
    TProfile* fVnLg2A = new TProfile("fVnLg2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnLg2C = new TProfile("fVnLg2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fVnK0g2A = new TProfile("fVnK0g2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnK0g2C = new TProfile("fVnK0g2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fVnXig2A = new TProfile("fVnXig2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnXig2C = new TProfile("fVnXig2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fVnOmg2A = new TProfile("fVnOmg2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fVnOmg2C = new TProfile("fVnOmg2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
        
    TProfile* fQRes = new TProfile("fQRes", ";mode; resolution", 4, 0, 4);
      
    
    TH1F* fMultC = new TH1F("fMultC", "; multiplicity; Counts", 600, 0, 600);
    TH1F* fMultF = new TH1F("fMultF", "; multiplicity; Counts", 600, 0, 600);
    TH1F* fMultB = new TH1F("fMultB", "; multiplicity; Counts", 600, 0, 600);
    
    
    TH2D* fEtaPhiAllF = new TH2D("fEtaPhiAllF", "; #eta; #phi", 88, 2.9, 4.1, 72, 0, 2*TMath::Pi());
    TH2D* fEtaPhiAllB = new TH2D("fEtaPhiAllB", "; #eta; #phi", 88, -4.1, -2.9, 72, 0, 2*TMath::Pi());
    
    
    
    TProfile* fUQAllg1A = new TProfile("fUQAllg1A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQAllg1C = new TProfile("fUQAllg1C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQPig1A = new TProfile("fUQPig1A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPig1C = new TProfile("fUQPig1C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQKg1A = new TProfile("fUQKg1A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQKg1C = new TProfile("fUQKg1C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQPg1A = new TProfile("fUQPg1A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPg1C = new TProfile("fUQPg1C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQLg1A = new TProfile("fUQLg1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQLg1C = new TProfile("fUQLg1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQK0g1A = new TProfile("fUQK0g1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQK0g1C = new TProfile("fUQK0g1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQXig1A = new TProfile("fUQXig1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQXig1C = new TProfile("fUQXig1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQOmg1A = new TProfile("fUQOmg1A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQOmg1C = new TProfile("fUQOmg1C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TH1I* fMultg1A = new TH1I("fMultg1A", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg1C = new TH1I("fMultg1C", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    
    TProfile* fUQAllg1Apos = new TProfile("fUQAllg1Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQAllg1Cpos = new TProfile("fUQAllg1Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQPig1Apos = new TProfile("fUQPig1Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPig1Cpos = new TProfile("fUQPig1Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQKg1Apos = new TProfile("fUQKg1Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQKg1Cpos = new TProfile("fUQKg1Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQPg1Apos = new TProfile("fUQPg1Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPg1Cpos = new TProfile("fUQPg1Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQLg1Apos = new TProfile("fUQLg1Apos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQLg1Cpos = new TProfile("fUQLg1Cpos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQXig1Apos = new TProfile("fUQXig1Apos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQXig1Cpos = new TProfile("fUQXig1Cpos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQOmg1Apos = new TProfile("fUQOmg1Apos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQOmg1Cpos = new TProfile("fUQOmg1Cpos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TH1I* fMultg1Apos = new TH1I("fMultg1Apos", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg1Cpos = new TH1I("fMultg1Cpos", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    
    TProfile* fUQAllg1Aneg = new TProfile("fUQAllg1Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQAllg1Cneg = new TProfile("fUQAllg1Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQPig1Aneg = new TProfile("fUQPig1Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPig1Cneg = new TProfile("fUQPig1Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQKg1Aneg = new TProfile("fUQKg1Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQKg1Cneg = new TProfile("fUQKg1Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQPg1Aneg = new TProfile("fUQPg1Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPg1Cneg = new TProfile("fUQPg1Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQLg1Aneg = new TProfile("fUQLg1Aneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQLg1Cneg = new TProfile("fUQLg1Cneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQXig1Aneg = new TProfile("fUQXig1Aneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQXig1Cneg = new TProfile("fUQXig1Cneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQOmg1Aneg = new TProfile("fUQOmg1Aneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQOmg1Cneg = new TProfile("fUQOmg1Cneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TH1I* fMultg1Aneg = new TH1I("fMultg1Aneg", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg1Cneg = new TH1I("fMultg1Cneg", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    TProfile* fResUQg1 = new TProfile("fResUQg1", "; centrality percentile; resolution", 3, 0, 3);

    
    
    
    TProfile* fUQAllg2A = new TProfile("fUQAllg2A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQAllg2C = new TProfile("fUQAllg2C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQPig2A = new TProfile("fUQPig2A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPig2C = new TProfile("fUQPig2C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQKg2A = new TProfile("fUQKg2A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQKg2C = new TProfile("fUQKg2C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQPg2A = new TProfile("fUQPg2A", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPg2C = new TProfile("fUQPg2C", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQLg2A = new TProfile("fUQLg2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQLg2C = new TProfile("fUQLg2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQK0g2A = new TProfile("fUQK0g2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQK0g2C = new TProfile("fUQK0g2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQXig2A = new TProfile("fUQXig2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQXig2C = new TProfile("fUQXig2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQOmg2A = new TProfile("fUQOmg2A", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQOmg2C = new TProfile("fUQOmg2C", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TH1I* fMultg2A = new TH1I("fMultg2A", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg2C = new TH1I("fMultg2C", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    
    TProfile* fUQAllg2Apos = new TProfile("fUQAllg2Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQAllg2Cpos = new TProfile("fUQAllg2Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQPig2Apos = new TProfile("fUQPig2Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPig2Cpos = new TProfile("fUQPig2Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQKg2Apos = new TProfile("fUQKg2Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQKg2Cpos = new TProfile("fUQKg2Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQPg2Apos = new TProfile("fUQPg2Apos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPg2Cpos = new TProfile("fUQPg2Cpos", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQLg2Apos = new TProfile("fUQLg2Apos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQLg2Cpos = new TProfile("fUQLg2Cpos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQXig2Apos = new TProfile("fUQXig2Apos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQXig2Cpos = new TProfile("fUQXig2Cpos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQOmg2Apos = new TProfile("fUQOmg2Apos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQOmg2Cpos = new TProfile("fUQOmg2Cpos", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
        
    TH1I* fMultg2Apos = new TH1I("fMultg2Apos", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg2Cpos = new TH1I("fMultg2Cpos", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    TProfile* fUQAllg2Aneg = new TProfile("fUQAllg2Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQAllg2Cneg = new TProfile("fUQAllg2Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQPig2Aneg = new TProfile("fUQPig2Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPig2Cneg = new TProfile("fUQPig2Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
        
    TProfile* fUQKg2Aneg = new TProfile("fUQKg2Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQKg2Cneg = new TProfile("fUQKg2Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQPg2Aneg = new TProfile("fUQPg2Aneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    TProfile* fUQPg2Cneg = new TProfile("fUQPg2Cneg", "; p_{T} (GeV/c); uQ", nPt, binsPt);
    
    TProfile* fUQLg2Aneg = new TProfile("fUQLg2Aneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQLg2Cneg = new TProfile("fUQLg2Cneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQXig2Aneg = new TProfile("fUQXig2Aneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQXig2Cneg = new TProfile("fUQXig2Cneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TProfile* fUQOmg2Aneg = new TProfile("fUQOmg2Aneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    TProfile* fUQOmg2Cneg = new TProfile("fUQOmg2Cneg", "; p_{T} (GeV/c); v_{n}", nPtV0Ca, binsPtV0Ca);
    
    TH1I* fMultg2Aneg = new TH1I("fMultg2Aneg", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg2Cneg = new TH1I("fMultg2Cneg", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    TH1I* fMultg1 = new TH1I("fMultg1", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg1pos = new TH1I("fMultg1pos", "; multiplicity; Counts", 600, -0.5, 599.5);
    TH1I* fMultg1neg = new TH1I("fMultg1neg", "; multiplicity; Counts", 600, -0.5, 599.5);
    
    
    TProfile* fResUQg2 = new TProfile("fResUQg2", "; centrality percentile; resolution", 3, 0, 3);
    TProfile* fResUQg2pos = new TProfile("fResUQg2pos", "; centrality percentile; resolution", 3, 0, 3);
    TProfile* fResUQg2neg = new TProfile("fResUQg2neg", "; centrality percentile; resolution", 3, 0, 3);
    
    
    
    
    for(Int_t n = 0; n < nEv; n++) {
    
        inTree->GetEntry(n);
    
        if((n+1)%500000 == 0)
            cout << "Analysis - Event: " << n+1 << "/" << nEv << endl;
        

        Double_t Qx = 0, Qy = 0;
        Double_t Q2x = 0, Q2y = 0;
        Double_t Q3x = 0, Q3y = 0;
        Int_t mult = 0, multXi = 0, multOm = 0, multL = 0, multK0 = 0, multPi = 0;

        Double_t QxEtaPos1 = 0, QyEtaPos1 = 0;
        Double_t QxEtaNeg1 = 0, QyEtaNeg1 = 0;
        Int_t multEtaPos1 = 0, multEtaNeg1 = 0;
        
        Double_t QxEtaPos1pos = 0, QyEtaPos1pos = 0;
        Double_t QxEtaNeg1pos = 0, QyEtaNeg1pos = 0;
        Int_t multEtaPos1pos = 0, multEtaNeg1pos = 0;
        
        Double_t QxEtaPos1neg = 0, QyEtaPos1neg = 0;
        Double_t QxEtaNeg1neg = 0, QyEtaNeg1neg = 0;
        Int_t multEtaPos1neg = 0, multEtaNeg1neg = 0;
        
        
        Double_t QxEtaPos2 = 0, QyEtaPos2 = 0;
        Double_t QxEtaNeg2 = 0, QyEtaNeg2 = 0;
        Int_t multEtaPos2 = 0, multEtaNeg2 = 0;
        
        Double_t QxEtaPos2pos = 0, QyEtaPos2pos = 0;
        Double_t QxEtaNeg2pos = 0, QyEtaNeg2pos = 0;
        Int_t multEtaPos2pos = 0, multEtaNeg2pos = 0;
        
        Double_t QxEtaPos2neg = 0, QyEtaPos2neg = 0;
        Double_t QxEtaNeg2neg = 0, QyEtaNeg2neg = 0;
        Int_t multEtaPos2neg = 0, multEtaNeg2neg = 0;
        
        
        Double_t QxG1 = 0, QyG1 = 0, QxG2 = 0, QyG2 = 0, QxG3 = 0, QyG3 = 0, QxG4 = 0, QyG4 = 0;
        Int_t multG1 = 0, multG2 = 0, multG3 = 0, multG4 = 0;
        
        Double_t QxEta1 = 0, QyEta1 = 0;
        Double_t QxEta1pos = 0, QyEta1pos = 0;
        Double_t QxEta1neg = 0, QyEta1neg = 0;
        Int_t multEta1 = 0, multEta1pos = 0, multEta1neg = 0;
        
        
        Double_t QxQC2EtaPos1 = 0, QyQC2EtaPos1 = 0, QxQC2EtaPos2 = 0, QyQC2EtaPos2 = 0, QxQC2EtaPos3 = 0, QyQC2EtaPos3 = 0;
        Double_t QxQC2EtaNeg1 = 0, QyQC2EtaNeg1 = 0, QxQC2EtaNeg2 = 0, QyQC2EtaNeg2 = 0, QxQC2EtaNeg3 = 0, QyQC2EtaNeg3 = 0;
        Int_t multQC2EtaPos1 = 0, multQC2EtaPos2 = 0, multQC2EtaPos3 = 0, multQC2EtaNeg1 = 0, multQC2EtaNeg2 = 0, multQC2EtaNeg3 = 0;
        
                
        Int_t nTracks = trackArray->GetEntries();

        for(Int_t jj = 0; jj < nTracks; jj++) {
    
            TrackN* trk = (TrackN*)trackArray->At(jj);
    
            if (!trk){
                delete trk;
                continue;
            }
            
            Double_t pt = TMath::Sqrt(trk->px*trk->px + trk->py*trk->py);
            Int_t charge = trk->q;
            
            
            if (pt < 0.2 || pt > 3. || charge == 0)
                continue;
            
            
            Double_t phiNew = trk->phi + TMath::Pi();
            Double_t trkEta = trk->eta;
   
                          
            Double_t sinHarm = TMath::Sin(nHarm*phiNew);
            Double_t cosHarm = TMath::Cos(nHarm*phiNew);
            
            
            if (TMath::Abs(trkEta) < 2.5){
                                        
                Double_t sinHarm2 = TMath::Sin(2.*nHarm*phiNew);
                Double_t cosHarm2 = TMath::Cos(2.*nHarm*phiNew);
                    
                Double_t sinHarm3 = TMath::Sin(3.*nHarm*phiNew);
                Double_t cosHarm3 = TMath::Cos(3.*nHarm*phiNew);
        
                Qx += cosHarm;
                Qy += sinHarm;
        
                Q2x += cosHarm2;
                Q2y += sinHarm2;
                    
                Q3x += cosHarm3;
                Q3y += sinHarm3;
                    
                mult++;
                            
                if (trkEta > 0.5 && trkEta < 1.){
                    
                    QxEtaPos1 += cosHarm;
                    QyEtaPos1 += sinHarm;
                    multEtaPos1++;
                    
                    if (charge > 0){
                        QxEtaPos1pos += cosHarm;
                        QyEtaPos1pos += sinHarm;
                        multEtaPos1pos++;
                    } else if (charge < 0){
                        QxEtaPos1neg += cosHarm;
                        QyEtaPos1neg += sinHarm;
                        multEtaPos1neg++;
                    }
                    
                }
                
                if (trkEta > -1. && trkEta < -0.5){
                    
                    QxEtaNeg1 += cosHarm;
                    QyEtaNeg1 += sinHarm;
                    multEtaNeg1++;
                    
                    if (charge > 0){
                        QxEtaNeg1pos += cosHarm;
                        QyEtaNeg1pos += sinHarm;
                        multEtaNeg1pos++;
                    } else if (charge < 0){
                        QxEtaNeg1neg += cosHarm;
                        QyEtaNeg1neg += sinHarm;
                        multEtaNeg1neg++;
                    }
                    
                }
                
                if (TMath::Abs(trkEta) < 1.){
                    
                    QxEta1 += cosHarm;
                    QyEta1 += sinHarm;
                    multEta1++;
                    
                    if (charge > 0){
                        QxEta1pos += cosHarm;
                        QyEta1pos += sinHarm;
                        multEta1pos++;
                    } else if (charge < 0){
                        QxEta1neg += cosHarm;
                        QyEta1neg += sinHarm;
                        multEta1neg++;
                    }
                    
                }
                
                
                
                if (trkEta > -2.5 && trkEta < -1.55){
                    QxG1 += cosHarm;
                    QyG1 += sinHarm;
                    multG1++;
                }
                    
                if (trkEta > -1.15 && trkEta < -0.2){
                    QxG2 += cosHarm;
                    QyG2 += sinHarm;
                    multG2++;
                }
                
                if (trkEta > 0.2 && trkEta < 1.15){
                    QxG3 += cosHarm;
                    QyG3 += sinHarm;
                    multG3++;
                }
                 
                if (trkEta > 1.55 && trkEta < 2.5){
                    QxG4 += cosHarm;
                    QyG4 += sinHarm;
                    multG4++;
                }
                
                                
                if (trkEta < -1.5){
                    QxQC2EtaNeg3 += cosHarm;
                    QyQC2EtaNeg3 += sinHarm;
                    multQC2EtaNeg3++;
                }
                
                if (trkEta < -1.){
                    QxQC2EtaNeg2 += cosHarm;
                    QyQC2EtaNeg2 += sinHarm;
                    multQC2EtaNeg2++;
                }
                
                if (trkEta < -0.5){
                    QxQC2EtaNeg1 += cosHarm;
                    QyQC2EtaNeg1 += sinHarm;
                    multQC2EtaNeg1++;
                }
                
                if (trkEta > 0.5){
                    QxQC2EtaPos1 += cosHarm;
                    QyQC2EtaPos1 += sinHarm;
                    multQC2EtaPos1++;
                }
                
                if (trkEta > 1.){
                    QxQC2EtaPos2 += cosHarm;
                    QyQC2EtaPos2 += sinHarm;
                    multQC2EtaPos2++;
                }
                
                if (trkEta > 1.5){
                    QxQC2EtaPos3 += cosHarm;
                    QyQC2EtaPos3 += sinHarm;
                    multQC2EtaPos3++;
                }
                                
            }
            
            
            if (trkEta > 3. && trkEta < 5.){
                
                QxEtaPos2 += cosHarm;
                QyEtaPos2 += sinHarm;
                multEtaPos2++;
                
                if (charge > 0){
                    QxEtaPos2pos += cosHarm;
                    QyEtaPos2pos += sinHarm;
                    multEtaPos2pos++;
                } else if (charge < 0){
                    QxEtaPos2neg += cosHarm;
                    QyEtaPos2neg += sinHarm;
                    multEtaPos2neg++;
                }
                
                fEtaPhiAllF->Fill(trkEta, phiNew);
            }
            
            if (trkEta > -5. && trkEta < -3.){
                QxEtaNeg2 += cosHarm;
                QyEtaNeg2 += sinHarm;
                multEtaNeg2++;
                
                if (charge > 0){
                    QxEtaNeg2pos += cosHarm;
                    QyEtaNeg2pos += sinHarm;
                    multEtaNeg2pos++;
                } else if (charge < 0){
                    QxEtaNeg2neg += cosHarm;
                    QyEtaNeg2neg += sinHarm;
                    multEtaNeg2neg++;
                }
                
                fEtaPhiAllB->Fill(trkEta, phiNew);
            }
    
        }
        
        
        
        if (!isMB){
            if (isLowM){
                if (multEta1 >= 50)
                    continue;
            } else {
                if (multEta1 < 80)
                    continue;
            }
        }
        
        
        
        fMultg1A->Fill(multEtaNeg1);
        fMultg1Apos->Fill(multEtaNeg1pos);
        fMultg1Aneg->Fill(multEtaNeg1neg);
        
        fMultg1C->Fill(multEtaPos1);
        fMultg1Cpos->Fill(multEtaPos1pos);
        fMultg1Cneg->Fill(multEtaPos1neg);
        
        fMultg2C->Fill(multEtaPos2);
        fMultg2Cpos->Fill(multEtaPos2pos);
        fMultg2Cneg->Fill(multEtaPos2neg);
        
        fMultg2A->Fill(multEtaNeg2);
        fMultg2Apos->Fill(multEtaNeg2pos);
        fMultg2Aneg->Fill(multEtaNeg2neg);
        
        fMultg1->Fill(multEta1);
        fMultg1pos->Fill(multEta1pos);
        fMultg1neg->Fill(multEta1neg);
        
        
        
        fMultC->Fill(mult);
        fMultF->Fill(multEtaPos2);
        fMultB->Fill(multEtaNeg2);
        

        if (multQC2EtaPos1 > 0 && multQC2EtaNeg1 > 0){
            Double_t nmG1 = multQC2EtaPos1 * multQC2EtaNeg1;
            Double_t uvG1 = (QxQC2EtaPos1*QxQC2EtaNeg1 + QyQC2EtaPos1*QyQC2EtaNeg1)/nmG1;
            fQC2AllQGap1->Fill(mult, uvG1);
        }
        
        if (multQC2EtaPos2 > 0 && multQC2EtaNeg2 > 0){
            Double_t nmG2 = multQC2EtaPos2 * multQC2EtaNeg2;
            Double_t uvG2 = (QxQC2EtaPos2*QxQC2EtaNeg2 + QyQC2EtaPos2*QyQC2EtaNeg2)/nmG2;
            fQC2AllQGap2->Fill(mult, uvG2);
        }
        
        if (multQC2EtaPos3 > 0 && multQC2EtaNeg3 > 0){
            Double_t nmG3 = multQC2EtaPos3 * multQC2EtaNeg3;
            Double_t uvG3 = (QxQC2EtaPos3*QxQC2EtaNeg3 + QyQC2EtaPos3*QyQC2EtaNeg3)/nmG3;
            fQC2AllQGap3->Fill(mult, uvG3);
        }
        
        
        if (mult > 5){
            
            Double_t nm = mult;
            Double_t nm1 = nm*(nm-1.);
            Double_t nm2 = nm1*(nm-2.);
            Double_t nm3 = nm2*(nm-3.);
            Double_t nm4 = nm3*(nm-4.);
            Double_t nm5 = nm4*(nm-5.);
      
            Double_t Qsqx=Qx*Qx-Qy*Qy;
            Double_t Qsqy=2.*Qx*Qy;
            Double_t QQ = Qx*Qx+Qy*Qy;
            Double_t QsqQsq = Qsqx*Qsqx+Qsqy*Qsqy;
            Double_t Q2Q2 = Q2x*Q2x+Q2y*Q2y;
            Double_t QsqQ2 = Qsqx*Q2x+Qsqy*Q2y;
            Double_t uv = (QQ-nm) / nm1;
            Double_t u2v2 = (Q2Q2-nm) / nm1;
            Double_t uuv2 = (QsqQ2 - nm -2.*nm1*uv - nm1*u2v2 )/ nm2;
            Double_t uuvv = (QsqQsq - nm -2.*nm1 -4.*nm1*(nm-1)*uv - 2.*nm2*uuv2 - nm1*u2v2 ) / nm3;
            Double_t uuuvvv = ((Qx*Qx + Qy*Qy)*(Qx*Qx + Qy*Qy)*(Qx*Qx + Qy*Qy) + 9.*(Q2x*Q2x + Q2y*Q2y)*(Qx*Qx + Qy*Qy) - 6.*(Q2x*Qx*Qx*Qx*Qx + Q2x*Qy*Qy*Qx*Qx + 2.*Q2y*Qx*Qx*Qx*Qy + 2.*Q2y*Qx*Qy*Qy*Qy - Q2x*Qy*Qy*Qy*Qy))/nm5 + 4.*(Q3x*Qx*Qx*Qx + Q3x*Qx*Qy*Qy + 3.*Q3y*Qx*Qx*Qy - Q3y*Qy*Qy*Qy - 3.*Q3x*Q2x*Qx - 3.*Q2y*Q3x*Qy + 3.*Q3y*Q2x*Qy + 3.*Q3y*Q2y*Qx)/nm5 + 2.*(9.*(nm-4)*(Q2x*Qx*Qx - Qy*Qy*Q2x + 2.*Q2y*Qx*Qy) + 2.*(Q3x*Q3x + Q3y*Q3y))/nm5 - 9.*((Qx*Qx + Qy*Qy)*(Qx*Qx + Qy*Qy) + Q2x*Q2x + Q2y*Q2y)/nm3/(nm-5.) + 18.*(Qx*Qx + Qy*Qy)/nm1/(nm-3.)/(nm-4.) - 6.*nm/nm3;
      
            fQC2AllQ->Fill(mult, uv);
            fQC4AllQ->Fill(mult, uuvv);
            fQC6AllQ->Fill(mult, uuuvvv);
        
        }
        
        
        if (multG1 != 0 && multG2 != 0 && multG3 != 0 && multG4 != 0){
          
            Double_t qqqqg = QxG1*QxG2*QxG3*QxG4 - QxG1*QxG2*QyG3*QyG4 + QxG1*QyG2*QxG3*QyG4 + QxG1*QyG2*QyG3*QxG4 + QyG1*QxG2*QxG3*QyG4 + QyG1*QxG2*QyG3*QxG4 - QyG1*QyG2*QxG3*QxG4 + QyG1*QyG2*QyG3*QyG4;
          
            fQC4AllQGap->Fill(mult, qqqqg/(multG1*multG2*multG3*multG4));
            fQC2AllQGap12->Fill(mult, (QxG1*QxG2 + QyG1*QyG2)/(multG1*multG2));
            fQC2AllQGap13->Fill(mult, (QxG1*QxG3 + QyG1*QyG3)/(multG1*multG3));
            fQC2AllQGap14->Fill(mult, (QxG1*QxG4 + QyG1*QyG4)/(multG1*multG4));
            fQC2AllQGap23->Fill(mult, (QxG2*QxG3 + QyG2*QyG3)/(multG2*multG3));
            fQC2AllQGap24->Fill(mult, (QxG2*QxG4 + QyG2*QyG4)/(multG2*multG4));
            fQC2AllQGap34->Fill(mult, (QxG3*QxG4 + QyG3*QyG4)/(multG3*multG4));

        }
    
       
        if (multEtaPos1 > 0 && multEtaNeg1 > 0){
            Double_t resGap1 = (QxEtaPos1*QxEtaNeg1 + QyEtaPos1*QyEtaNeg1) / (multEtaPos1*multEtaNeg1);
            fQRes->Fill(0., resGap1);
        }
        
        if (multEta1 > 0 && multEtaPos2 > 0 && multEtaNeg2 > 0){
            
            Double_t resGapCF = (QxEta1*QxEtaPos2 + QyEta1*QyEtaPos2) / (multEta1*multEtaPos2);
            Double_t resGapCB = (QxEta1*QxEtaNeg2 + QyEta1*QyEtaNeg2) / (multEta1*multEtaNeg2);
            Double_t resGapFB = (QxEtaPos2*QxEtaNeg2 + QyEtaPos2*QyEtaNeg2) / (multEtaPos2*multEtaNeg2);
            
            fQRes->Fill(1., resGapCF);
            fQRes->Fill(2., resGapCB);
            fQRes->Fill(3., resGapFB);
        }
        
        
        Double_t resUQ1 = QxEtaPos1*QxEtaNeg1 + QyEtaPos1*QyEtaNeg1;
        Double_t resUQ1pos = QxEtaPos1pos*QxEtaNeg1pos + QyEtaPos1pos*QyEtaNeg1pos;
        Double_t resUQ1neg = QxEtaPos1neg*QxEtaNeg1neg + QyEtaPos1neg*QyEtaNeg1neg;
        
        fResUQg1->Fill(0., resUQ1);
        fResUQg1->Fill(1., resUQ1pos);
        fResUQg1->Fill(2., resUQ1neg);
        
        
        Double_t resUQCF = QxEta1*QxEtaPos2 + QyEta1*QyEtaPos2;
        Double_t resUQCB = QxEta1*QxEtaNeg2 + QyEta1*QyEtaNeg2;
        Double_t resUQFB = QxEtaPos2*QxEtaNeg2 + QyEtaPos2*QyEtaNeg2;
        
        fResUQg2->Fill(0., resUQCF);
        fResUQg2->Fill(1., resUQCB);
        fResUQg2->Fill(2., resUQFB);
        
        
        Double_t resUQCFpos = QxEta1pos*QxEtaPos2pos + QyEta1pos*QyEtaPos2pos;
        Double_t resUQCBpos = QxEta1pos*QxEtaNeg2pos + QyEta1pos*QyEtaNeg2pos;
        Double_t resUQFBpos = QxEtaPos2pos*QxEtaNeg2pos + QyEtaPos2pos*QyEtaNeg2pos;
        
        fResUQg2pos->Fill(0., resUQCFpos);
        fResUQg2pos->Fill(1., resUQCBpos);
        fResUQg2pos->Fill(2., resUQFBpos);
        
        
        Double_t resUQCFneg = QxEta1neg*QxEtaPos2neg + QyEta1neg*QyEtaPos2neg;
        Double_t resUQCBneg = QxEta1neg*QxEtaNeg2neg + QyEta1neg*QyEtaNeg2neg;
        Double_t resUQFBneg = QxEtaPos2neg*QxEtaNeg2neg + QyEtaPos2neg*QyEtaNeg2neg;
        
        fResUQg2neg->Fill(0., resUQCFneg);
        fResUQg2neg->Fill(1., resUQCBneg);
        fResUQg2neg->Fill(2., resUQFBneg);
        

        
        
        Int_t multN = 0;
        
        for(Int_t jjj = 0; jjj < nTracks; jjj++) {
            
            TrackN* tk = (TrackN*)trackArray->At(jjj);
            
            if (!tk){
                delete tk;
                continue;
            }
            
            Double_t ptN = TMath::Sqrt(tk->px*tk->px + tk->py*tk->py);
            //Double_t yN = 0.5*TMath::Log((tk->e + tk->pz) / (tk->e - tk->pz));
            Double_t etaN = tk->eta;
            
            
            if (ptN < 0.2 || TMath::Abs(etaN) >= 2.5)
                continue;
            
            
            Double_t phiN = tk->phi + TMath::Pi();
                    
            Double_t sinHa = TMath::Sin(nHarm*phiN);
            Double_t cosHa = TMath::Cos(nHarm*phiN);
            
            if (ptN < 3.){
                fSinAll->Fill(mult, sinHa);
                fCosAll->Fill(mult, cosHa);
            }
            
            
            
            if (ptN >= 20. || TMath::Abs(etaN) >= 1.)
                continue;
            
            
            Int_t chargeN = tk->q;
            Int_t absidN = TMath::Abs(tk->pdgCode);
                        
            
            if (chargeN != 0){
            	fPtAll->Fill(ptN);
            	fEtaPhiAll->Fill(etaN, phiN);
                multN++;
            }
             
             
            if(absidN == 3334){
                fPtOm->Fill(ptN);
                fEtaPhiOm->Fill(etaN, phiN);
                multOm++;
            }
            
            if (absidN == 3312){
                fPtXi->Fill(ptN);
                fEtaPhiXi->Fill(etaN, phiN);
                multXi++;
            }
                
            if (absidN == 3122){
                fPtL->Fill(ptN);
                fEtaPhiL->Fill(etaN, phiN);
                multL++;
            }
                
            if (absidN == 310){
                fPtK0->Fill(ptN);
                fEtaPhiK0->Fill(etaN, phiN);
                multK0++;
            }
                
            if (absidN == 211){
                fPtPi->Fill(ptN);
                fEtaPhiPi->Fill(etaN, phiN);
                multPi++;
            }
      
            
            
            Double_t uqCg2 = cosHa*QxEtaPos2 + sinHa*QyEtaPos2;
            Double_t uqAg2 = cosHa*QxEtaNeg2 + sinHa*QyEtaNeg2;
            
            if (chargeN != 0 ){
                fUQAllg2C->Fill(ptN, uqCg2);
                fUQAllg2A->Fill(ptN, uqAg2);
            }
            
            if (absidN == 211){
                fUQPig2C->Fill(ptN, uqCg2);
                fUQPig2A->Fill(ptN, uqAg2);
            }
            
            if (absidN == 321){
                fUQKg2C->Fill(ptN, uqCg2);
                fUQKg2A->Fill(ptN, uqAg2);
            }
            
            if (absidN == 2212){
                fUQPg2C->Fill(ptN, uqCg2);
                fUQPg2A->Fill(ptN, uqAg2);
            }
            
            if(absidN == 3334){
                fUQOmg2C->Fill(ptN, uqCg2);
                fUQOmg2A->Fill(ptN, uqAg2);
            }
            
            if (absidN == 3312){
                fUQXig2C->Fill(ptN, uqCg2);
                fUQXig2A->Fill(ptN, uqAg2);
            }
            
            if (absidN == 3122){
                fUQLg2C->Fill(ptN, uqCg2);
                fUQLg2A->Fill(ptN, uqAg2);
            }
            
            if (absidN == 310){
                fUQK0g2C->Fill(ptN, uqCg2);
                fUQK0g2A->Fill(ptN, uqAg2);
            }
            
            
            if (chargeN > 0){
                
                Double_t uqCg2pos = cosHa*QxEtaPos2pos + sinHa*QyEtaPos2pos;
                Double_t uqAg2pos = cosHa*QxEtaNeg2pos + sinHa*QyEtaNeg2pos;
                
                fUQAllg2Cpos->Fill(ptN, uqCg2pos);
                fUQAllg2Apos->Fill(ptN, uqAg2pos);
                
                if (absidN == 211){
                    fUQPig2Cpos->Fill(ptN, uqCg2pos);
                    fUQPig2Apos->Fill(ptN, uqAg2pos);
                }
                
                if (absidN == 321){
                    fUQKg2Cpos->Fill(ptN, uqCg2pos);
                    fUQKg2Apos->Fill(ptN, uqAg2pos);
                }
                
                if (absidN == 2212){
                    fUQPg2Cpos->Fill(ptN, uqCg2pos);
                    fUQPg2Apos->Fill(ptN, uqAg2pos);
                }
                
                if(absidN == 3334){
                    fUQOmg2Cpos->Fill(ptN, uqCg2pos);
                    fUQOmg2Apos->Fill(ptN, uqAg2pos);
                }
                
                if (absidN == 3312){
                    fUQXig2Cpos->Fill(ptN, uqCg2pos);
                    fUQXig2Apos->Fill(ptN, uqAg2pos);
                }
                
                if (absidN == 3122){
                    fUQLg2Cpos->Fill(ptN, uqCg2pos);
                    fUQLg2Apos->Fill(ptN, uqAg2pos);
                }
                
            } else if (chargeN < 0){
                
                Double_t uqCg2neg = cosHa*QxEtaPos2neg + sinHa*QyEtaPos2neg;
                Double_t uqAg2neg = cosHa*QxEtaNeg2neg + sinHa*QyEtaNeg2neg;
                
                fUQAllg2Cneg->Fill(ptN, uqCg2neg);
                fUQAllg2Aneg->Fill(ptN, uqAg2neg);
                
                if (absidN == 211){
                    fUQPig2Cneg->Fill(ptN, uqCg2neg);
                    fUQPig2Aneg->Fill(ptN, uqAg2neg);
                }
                
                if (absidN == 321){
                    fUQKg2Cneg->Fill(ptN, uqCg2neg);
                    fUQKg2Aneg->Fill(ptN, uqAg2neg);
                }
                
                if (absidN == 2212){
                    fUQPg2Cneg->Fill(ptN, uqCg2neg);
                    fUQPg2Aneg->Fill(ptN, uqAg2neg);
                }
                
                if(absidN == 3334){
                    fUQOmg2Cneg->Fill(ptN, uqCg2neg);
                    fUQOmg2Aneg->Fill(ptN, uqAg2neg);
                }
                
                if (absidN == 3312){
                    fUQXig2Cneg->Fill(ptN, uqCg2neg);
                    fUQXig2Aneg->Fill(ptN, uqAg2neg);
                }
                
                if (absidN == 3122){
                    fUQLg2Cneg->Fill(ptN, uqCg2neg);
                    fUQLg2Aneg->Fill(ptN, uqAg2neg);
                }
                
            }
            
            
            
            Double_t harmGap1C = cosHa * QxEtaPos1 + sinHa * QyEtaPos1;
            Double_t harmGap1A = cosHa * QxEtaNeg1 + sinHa * QyEtaNeg1;
            
            if (etaN > 0.5){
                
                fUQAllg1C->Fill(ptN, harmGap1A);
                
                if (chargeN != 0)
                    fUQAllg1C->Fill(ptN, harmGap1A);
                
                if (absidN == 211)
                    fUQPig1C->Fill(ptN, harmGap1A);
                                
                if (absidN == 321)
                    fUQKg1C->Fill(ptN, harmGap1A);
                
                if (absidN == 2212)
                    fUQPg1C->Fill(ptN, harmGap1A);
                
                if(absidN == 3334)
                    fUQOmg1C->Fill(ptN, harmGap1A);
  
                if (absidN == 3312)
                    fUQXig1C->Fill(ptN, harmGap1A);
                
                if (absidN == 3122)
                    fUQLg1C->Fill(ptN, harmGap1A);
   
                if (absidN == 310)
                    fUQK0g1C->Fill(ptN, harmGap1A);

                
                if (chargeN > 0){
                    
                    Double_t uqCg1pos = cosHa*QxEtaNeg1pos + sinHa*QyEtaNeg1pos;
                    
                    fUQAllg1Cpos->Fill(ptN, uqCg1pos);
                    
                    if (absidN == 211)
                        fUQPig1Cpos->Fill(ptN, uqCg1pos);

                    if (absidN == 321)
                        fUQKg1Cpos->Fill(ptN, uqCg1pos);

                    if (absidN == 2212)
                        fUQPg1Cpos->Fill(ptN, uqCg1pos);
                    
                    if(absidN == 3334)
                        fUQOmg1Cpos->Fill(ptN, uqCg1pos);
      
                    if (absidN == 3312)
                        fUQXig1Cpos->Fill(ptN, uqCg1pos);
                    
                    if (absidN == 3122)
                        fUQLg1Cpos->Fill(ptN, uqCg1pos);
                    
                } else if (chargeN < 0){

                    Double_t uqCg1neg = cosHa*QxEtaNeg1neg + sinHa*QyEtaNeg1neg;
                    
                    fUQAllg1Cneg->Fill(ptN, uqCg1neg);
                    
                    if (absidN == 211)
                        fUQPig1Cneg->Fill(ptN, uqCg1neg);
                    
                    if (absidN == 321)
                        fUQKg1Cneg->Fill(ptN, uqCg1neg);
                    
                    if (absidN == 2212)
                        fUQPg1Cneg->Fill(ptN, uqCg1neg);
                    
                    if(absidN == 3334)
                        fUQOmg1Cneg->Fill(ptN, uqCg1neg);
      
                    if (absidN == 3312)
                        fUQXig1Cneg->Fill(ptN, uqCg1neg);
                    
                    if (absidN == 3122)
                        fUQLg1Cneg->Fill(ptN, uqCg1neg);
                    
                }
        
            }
            
            
            if (etaN < -0.5){
                
                fUQAllg1A->Fill(ptN, harmGap1C);
                
                if (chargeN != 0)
                    fUQAllg1A->Fill(ptN, harmGap1C);
                
                if (absidN == 211)
                    fUQPig1A->Fill(ptN, harmGap1C);
                                
                if (absidN == 321)
                    fUQKg1A->Fill(ptN, harmGap1C);
                
                if (absidN == 2212)
                    fUQPg1A->Fill(ptN, harmGap1C);
                
                if(absidN == 3334)
                    fUQOmg1A->Fill(ptN, harmGap1C);
  
                if (absidN == 3312)
                    fUQXig1A->Fill(ptN, harmGap1C);
                
                if (absidN == 3122)
                    fUQLg1A->Fill(ptN, harmGap1C);
   
                if (absidN == 310)
                    fUQK0g1A->Fill(ptN, harmGap1C);

                
                if (chargeN > 0){
                    
                    Double_t uqAg1pos = cosHa*QxEtaPos1pos + sinHa*QyEtaPos1pos;
                    
                    fUQAllg1Apos->Fill(ptN, uqAg1pos);
                    
                    if (absidN == 211)
                        fUQPig1Apos->Fill(ptN, uqAg1pos);

                    if (absidN == 321)
                        fUQKg1Apos->Fill(ptN, uqAg1pos);

                    if (absidN == 2212)
                        fUQPg1Apos->Fill(ptN, uqAg1pos);
                    
                    if(absidN == 3334)
                        fUQOmg1Apos->Fill(ptN, uqAg1pos);
      
                    if (absidN == 3312)
                        fUQXig1Apos->Fill(ptN, uqAg1pos);
                    
                    if (absidN == 3122)
                        fUQLg1Apos->Fill(ptN, uqAg1pos);
                    
                } else if (chargeN < 0){

                    Double_t uqAg1neg = cosHa*QxEtaPos1neg + sinHa*QyEtaPos1neg;
                    
                    fUQAllg1Aneg->Fill(ptN, uqAg1neg);
                    
                    if (absidN == 211)
                        fUQPig1Aneg->Fill(ptN, uqAg1neg);
                    
                    if (absidN == 321)
                        fUQKg1Aneg->Fill(ptN, uqAg1neg);
                    
                    if (absidN == 2212)
                        fUQPg1Aneg->Fill(ptN, uqAg1neg);
                    
                    if(absidN == 3334)
                        fUQOmg1Aneg->Fill(ptN, uqAg1neg);
      
                    if (absidN == 3312)
                        fUQXig1Aneg->Fill(ptN, uqAg1neg);
                    
                    if (absidN == 3122)
                        fUQLg1Aneg->Fill(ptN, uqAg1neg);
                    
                }
        
            }
     

                    
            if (multEtaPos2 > 0){
                
                Double_t vnCg2 = (cosHa * QxEtaPos2 + sinHa * QyEtaPos2) / multEtaPos2;
                
                if (chargeN != 0)
                    fVnAllg2C->Fill(ptN, vnCg2);
                                
                if (absidN == 211)
                    fVnPig2C->Fill(ptN, vnCg2);
                
                if (absidN == 321)
                    fVnKg2C->Fill(ptN, vnCg2);
                
                if (absidN == 2212)
                    fVnPg2C->Fill(ptN, vnCg2);
                
                if(absidN == 3334)
                    fVnOmg2C->Fill(ptN, vnCg2);
                
                if (absidN == 3312)
                    fVnXig2C->Fill(ptN, vnCg2);
                
                if (absidN == 3122)
                    fVnLg2C->Fill(ptN, vnCg2);
                
                if (absidN == 310)
                    fVnK0g2C->Fill(ptN, vnCg2);
                
            }
            
            
            if (multEtaNeg2 > 0){
                
                Double_t vnAg2 = (cosHa * QxEtaNeg2 + sinHa * QyEtaNeg2) / multEtaNeg2;
                
                if (chargeN != 0)
                    fVnAllg2A->Fill(ptN, vnAg2);
                
                if (absidN == 211)
                    fVnPig2A->Fill(ptN, vnAg2);
                
                if (absidN == 321)
                    fVnKg2A->Fill(ptN, vnAg2);
                
                if (absidN == 2212)
                    fVnPg2A->Fill(ptN, vnAg2);
                
                if(absidN == 3334)
                    fVnOmg2A->Fill(ptN, vnAg2);
                
                if (absidN == 3312)
                    fVnXig2A->Fill(ptN, vnAg2);
                
                if (absidN == 3122)
                    fVnLg2A->Fill(ptN, vnAg2);
                
                if (absidN == 310)
                    fVnK0g2A->Fill(ptN, vnAg2);
                
            }
            
            
            
            if (etaN > 0.5 && multEtaNeg1 > 0){
                
                Double_t vnCg1 = harmGap1A / multEtaNeg1;
                
                if (chargeN != 0)
                    fVnAllg1C->Fill(ptN, vnCg1);
                
                if (absidN == 211)
                    fVnPig1C->Fill(ptN, vnCg1);
                
                if (absidN == 321)
                    fVnKg1C->Fill(ptN, vnCg1);
                
                if (absidN == 2212)
                    fVnPg1C->Fill(ptN, vnCg1);
                
                if(absidN == 3334)
                    fVnOmg1C->Fill(ptN, vnCg1);
                
                if (absidN == 3312)
                    fVnXig1C->Fill(ptN, vnCg1);
                
                if (absidN == 3122)
                    fVnLg1C->Fill(ptN, vnCg1);
                
                if (absidN == 310)
                    fVnK0g1C->Fill(ptN, vnCg1);
                
            }
            
            if (etaN < -0.5 && multEtaPos1 > 0){
                
                Double_t vnAg1 = harmGap1C / multEtaPos1;
                
                if (chargeN != 0)
                    fVnAllg1A->Fill(ptN, vnAg1);
                
                if (absidN == 211)
                    fVnPig1A->Fill(ptN, vnAg1);
                
                if (absidN == 321)
                    fVnKg1A->Fill(ptN, vnAg1);
                
                if (absidN == 2212)
                    fVnPg1A->Fill(ptN, vnAg1);
                
                if(absidN == 3334)
                    fVnOmg1A->Fill(ptN, vnAg1);
                
                if (absidN == 3312)
                    fVnXig1A->Fill(ptN, vnAg1);
                
                if (absidN == 3122)
                    fVnLg1A->Fill(ptN, vnAg1);
                
                if (absidN == 310)
                    fVnK0g1A->Fill(ptN, vnAg1);
                
            }
            
        }
        
        
        fMult->Fill(multN);
        fMultOm->Fill(multOm);
        fMultXi->Fill(multXi);
        fMultL->Fill(multL);
        fMultK0->Fill(multK0);
        fMultPi->Fill(multPi);
        
    }
    
    
    TFile* out = new TFile(outFileName, "RECREATE");

    fEtaPhiAll->Write();
    fPtAll->Write();
    fMult->Write();
        
    fPtOm->Write();
    fEtaPhiOm->Write();
    fMultOm->Write();
        
    fPtXi->Write();
    fEtaPhiXi->Write();
    fMultXi->Write();
        
    fPtL->Write();
    fEtaPhiL->Write();
    fMultL->Write();
        
    fPtK0->Write();
    fEtaPhiK0->Write();
    fMultK0->Write();

    fPtPi->Write();
    fEtaPhiPi->Write();
    fMultPi->Write();

        
    fSinAll->Write();
    fCosAll->Write();

    fQC2AllQ->Write();
    fQC4AllQ->Write();
    fQC6AllQ->Write();
    
    fQC2AllQGap1->Write();
    fQC2AllQGap2->Write();
    fQC2AllQGap3->Write();
    
    fQC4AllQGap->Write();
    fQC2AllQGap12->Write();
    fQC2AllQGap13->Write();
    fQC2AllQGap14->Write();
    fQC2AllQGap23->Write();
    fQC2AllQGap24->Write();
    fQC2AllQGap34->Write();

    
    fVnAllg1A->Write();
    fVnAllg1C->Write();
    
    fVnPig1A->Write();
    fVnPig1C->Write();
    
    fVnKg1A->Write();
    fVnKg1C->Write();
    
    fVnPg1A->Write();
    fVnPg1C->Write();
    
    fVnLg1A->Write();
    fVnLg1C->Write();
    
    fVnK0g1A->Write();
    fVnK0g1C->Write();
    
    fVnXig1A->Write();
    fVnXig1C->Write();
    
    fVnOmg1A->Write();
    fVnOmg1C->Write();
    
    
    fVnAllg2A->Write();
    fVnAllg2C->Write();
    
    fVnPig2A->Write();
    fVnPig2C->Write();
    
    fVnKg2A->Write();
    fVnKg2C->Write();
    
    fVnPg2A->Write();
    fVnPg2C->Write();
    
    fVnLg2A->Write();
    fVnLg2C->Write();
    
    fVnK0g2A->Write();
    fVnK0g2C->Write();
    
    fVnXig2A->Write();
    fVnXig2C->Write();
    
    fVnOmg2A->Write();
    fVnOmg2C->Write();

    fQRes->Write();
    
    
    fMultC->Write();
    fMultF->Write();
    fMultB->Write();
    
    
    fEtaPhiAllF->Write();
    fEtaPhiAllB->Write();
    
    
    fUQAllg1A->Write();
    fUQAllg1C->Write();
        
    fUQPig1A->Write();
    fUQPig1C->Write();
        
    fUQKg1A->Write();
    fUQKg1C->Write();
    
    fUQPg1A->Write();
    fUQPg1C->Write();
    
    fUQLg1A->Write();
    fUQLg1C->Write();
    
    fUQK0g1A->Write();
    fUQK0g1C->Write();
    
    fUQXig1A->Write();
    fUQXig1C->Write();
    
    fUQOmg1A->Write();
    fUQOmg1C->Write();
    
    fMultg1A->Write();
    fMultg1C->Write();
    
    
    fUQAllg1Apos->Write();
    fUQAllg1Cpos->Write();
        
    fUQPig1Apos->Write();
    fUQPig1Cpos->Write();
        
    fUQKg1Apos->Write();
    fUQKg1Cpos->Write();
    
    fUQPg1Apos->Write();
    fUQPg1Cpos->Write();
    
    fUQLg1Apos->Write();
    fUQLg1Cpos->Write();
    
    fUQXig1Apos->Write();
    fUQXig1Cpos->Write();
    
    fUQOmg1Apos->Write();
    fUQOmg1Cpos->Write();
    
    fMultg1Apos->Write();
    fMultg1Cpos->Write();
    
    
    fUQAllg1Aneg->Write();
    fUQAllg1Cneg->Write();
        
    fUQPig1Aneg->Write();
    fUQPig1Cneg->Write();
        
    fUQKg1Aneg->Write();
    fUQKg1Cneg->Write();
    
    fUQPg1Aneg->Write();
    fUQPg1Cneg->Write();
    
    fUQLg1Aneg->Write();
    fUQLg1Cneg->Write();
    
    fUQXig1Aneg->Write();
    fUQXig1Cneg->Write();
    
    fUQOmg1Aneg->Write();
    fUQOmg1Cneg->Write();
    
    fMultg1Aneg->Write();
    fMultg1Cneg->Write();
    
    fResUQg1->Write();

    
    fUQAllg2A->Write();
    fUQAllg2C->Write();
        
    fUQPig2A->Write();
    fUQPig2C->Write();
        
    fUQKg2A->Write();
    fUQKg2C->Write();
    
    fUQPg2A->Write();
    fUQPg2C->Write();
    
    fUQLg2A->Write();
    fUQLg2C->Write();
    
    fUQK0g2A->Write();
    fUQK0g2C->Write();
    
    fUQXig2A->Write();
    fUQXig2C->Write();
    
    fUQOmg2A->Write();
    fUQOmg2C->Write();
    
    fMultg2A->Write();
    fMultg2C->Write();
    
    
    fUQAllg2Apos->Write();
    fUQAllg2Cpos->Write();
        
    fUQPig2Apos->Write();
    fUQPig2Cpos->Write();
        
    fUQKg2Apos->Write();
    fUQKg2Cpos->Write();
    
    fUQPg2Apos->Write();
    fUQPg2Cpos->Write();
    
    fUQLg2Apos->Write();
    fUQLg2Cpos->Write();
    
    fUQXig2Apos->Write();
    fUQXig2Cpos->Write();
    
    fUQOmg2Apos->Write();
    fUQOmg2Cpos->Write();
    
    fMultg2Apos->Write();
    fMultg2Cpos->Write();
    
    
    fUQAllg2Aneg->Write();
    fUQAllg2Cneg->Write();
        
    fUQPig2Aneg->Write();
    fUQPig2Cneg->Write();
        
    fUQKg2Aneg->Write();
    fUQKg2Cneg->Write();
    
    fUQPg2Aneg->Write();
    fUQPg2Cneg->Write();
    
    fUQLg2Aneg->Write();
    fUQLg2Cneg->Write();
    
    fUQXig2Aneg->Write();
    fUQXig2Cneg->Write();
    
    fUQOmg2Aneg->Write();
    fUQOmg2Cneg->Write();
    
    fMultg2Aneg->Write();
    fMultg2Cneg->Write();
    
    
    fMultg1->Write();
    fMultg1pos->Write();
    fMultg1neg->Write();
    
    
    fResUQg2->Write();
    fResUQg2pos->Write();
    fResUQg2neg->Write();
     
    out->Close();
    
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

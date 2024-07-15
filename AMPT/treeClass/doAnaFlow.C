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


void doAnaFlow(Float_t nHarm, Bool_t isMB, const Char_t* inFileName, const Char_t* outFileName, Bool_t useChain)
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
    TH1I* fMult = new TH1I("fMult", "; multiplicity; Counts", 500, -0.5, 499.5);
    
    
    //pi+K0+lambda+xi+omega
    TH1D* fPtOm = new TH1D("fPtOm","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiOm = new TH2D("fEtaPhiOm","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultOm = new TH1I("fMultOm", "; multiplicity; Counts", 500, -0.5, 499.5);
    
    TH1D* fPtXi = new TH1D("fPtXi","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiXi = new TH2D("fEtaPhiXi","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultXi = new TH1I("fMultXi", "; multiplicity; Counts", 500, -0.5, 499.5);
    
    TH1D* fPtL = new TH1D("fPtL","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiL = new TH2D("fEtaPhiL","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultL = new TH1I("fMultL", "; multiplicity; Counts", 500, -0.5, 499.5);
    
    TH1D* fPtK0 = new TH1D("fPtK0","; p_{T} (GeV/c); Counts", nPtV0Ca, binsPtV0Ca);
    TH2D* fEtaPhiK0 = new TH2D("fEtaPhiK0","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultK0 = new TH1I("fMultK0", "; multiplicity; Counts", 500, -0.5, 499.5);
    
    TH1D* fPtPi = new TH1D("fPtPi","; p_{T} (GeV/c); Counts", nPt, binsPt);
    TH2D* fEtaPhiPi = new TH2D("fEtaPhiPi","; #eta; #varphi", 40, -1., 1., 72, 0, 2*TMath::Pi());
    TH1I* fMultPi = new TH1I("fMultPi", "; multiplicity; Counts", 500, -0.5, 499.5);
    
    
    
    TProfile* fSinAll = new TProfile("fSinAll", "; multiplicity; sin(n*#phi)", 300, 0, 300);
    TProfile* fCosAll = new TProfile("fCosAll", "; multiplicity; cos(n*#phi)", 300, 0, 300);

    TProfile* fQC2AllQ = new TProfile("fQC2AllQ", "; multiplicity; QC{2}", 300, 0, 300);
    TProfile* fQC4AllQ = new TProfile("fQC4AllQ", "; multiplicity; QC{4}", 300, 0, 300);
    TProfile* fQC6AllQ = new TProfile("fQC6AllQ", "; multiplicity; QC{6}", 300, 0, 300);

    
    
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
      
    
    TH1F* fMultC = new TH1F("fMultC", "; multiplicity; Counts", 500, 0, 500);
    TH1F* fMultF = new TH1F("fMultF", "; multiplicity; Counts", 500, 0, 500);
    TH1F* fMultB = new TH1F("fMultB", "; multiplicity; Counts", 500, 0, 500);
    
    
    TH2D* fEtaPhiAllF = new TH2D("fEtaPhiAllF", "; #eta; #phi", 88, 2.9, 4.1, 72, 0, 2*TMath::Pi());
    TH2D* fEtaPhiAllB = new TH2D("fEtaPhiAllB", "; #eta; #phi", 88, -4.1, -2.9, 72, 0, 2*TMath::Pi());
    
    
    
    
    for(Int_t n = 0; n < nEv; n++) {
    
        inTree->GetEntry(n);
    
        if((n+1)%500000 == 0)
            cout << "Analysis - Event: " << n+1 << "/" << nEv << endl;
        
        
        if (!isMB){
            if (eventData->multChEta1 < 100)
                continue;
        }

        
        
        Double_t Qx = 0, Qy = 0;
        Double_t Q2x = 0, Q2y = 0;
        Double_t Q3x = 0, Q3y = 0;
        Int_t mult = 0, multXi = 0, multOm = 0, multL = 0, multK0 = 0, multPi = 0;

        Double_t QxEtaPos1 = 0, QyEtaPos1 = 0;
        Double_t QxEtaNeg1 = 0, QyEtaNeg1 = 0;
        Int_t multEtaPos1 = 0, multEtaNeg1 = 0;
        
        Double_t QxEtaPos2 = 0, QyEtaPos2 = 0;
        Double_t QxEtaNeg2 = 0, QyEtaNeg2 = 0;
        Int_t multEtaPos2 = 0, multEtaNeg2 = 0;
        
                
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
            
            
            if (TMath::Abs(trkEta) < 2.){
                                        
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
                            
                    if (trkEta > 0.5){
                        QxEtaPos1 += cosHarm;
                        QyEtaPos1 += sinHarm;
                        multEtaPos1++;
                    } else if (trkEta < -0.5){
                        QxEtaNeg1 += cosHarm;
                        QyEtaNeg1 += sinHarm;
                        multEtaNeg1++;
                    }
                
            }
            
            
            if (trkEta > 3. && trkEta < 4.){
                QxEtaPos2 += cosHarm;
                QyEtaPos2 += sinHarm;
                multEtaPos2++;
                
                fEtaPhiAllF->Fill(trkEta, phiNew);
            }
            
            if (trkEta > -4. && trkEta < -3.){
                QxEtaNeg2 += cosHarm;
                QyEtaNeg2 += sinHarm;
                multEtaNeg2++;
                
                fEtaPhiAllB->Fill(trkEta, phiNew);
            }
    
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
    
       
        if (multEtaPos1 > 0 && multEtaNeg1 > 0){
            Double_t resGap1 = (QxEtaPos1 * QxEtaNeg1 + QyEtaPos1 * QyEtaNeg1) / (multEtaPos1 * multEtaNeg1);
            fQRes->Fill(0., resGap1);
        }
        
        if (mult > 0 && multEtaPos2 > 0 && multEtaNeg2 > 0){
            
            Double_t resGapCF = (Qx*QxEtaPos2 + Qy*QyEtaPos2) / (mult*multEtaPos2);
            Double_t resGapCB = (Qx*QxEtaNeg2 + Qy*QyEtaNeg2) / (mult*multEtaNeg2);
            Double_t resGapFB = (QxEtaPos2 * QxEtaNeg2 + QyEtaPos2 * QyEtaNeg2) / (multEtaPos2 * multEtaNeg2);
            
            fQRes->Fill(1., resGapCF);
            fQRes->Fill(2., resGapCB);
            fQRes->Fill(3., resGapFB);
        }
        
        
        fMultC->Fill(mult);
        fMultF->Fill(multEtaPos2);
        fMultB->Fill(multEtaNeg2);

        
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
            
            
            if (ptN < 0.2 || ptN > 20. || TMath::Abs(etaN) >= 1.)
                continue;
            
            
            Double_t phiN = tk->phi + TMath::Pi();
            Int_t chargeN = tk->q;
            Int_t absidN = TMath::Abs(tk->pdgCode);
            
            
            Double_t sinHa = TMath::Sin(nHarm*phiN);
            Double_t cosHa = TMath::Cos(nHarm*phiN);
            
            fSinAll->Fill(mult, sinHa);
            fCosAll->Fill(mult, cosHa);
            
            
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
      
            
            
            
            Double_t harmGap1C = cosHa * QxEtaPos1 + sinHa * QyEtaPos1;
            Double_t harmGap1A = cosHa * QxEtaNeg1 + sinHa * QyEtaNeg1;
            
            if (multEtaPos2 > 0){
                
                Double_t vnCg2 = (cosHa * QxEtaPos2 + sinHa * QyEtaPos2) / multEtaPos2;
                
                if (chargeN != 0 )
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
                
                if (chargeN != 0 )
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
                
                if (chargeN != 0 )
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
                
                if (chargeN != 0 )
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

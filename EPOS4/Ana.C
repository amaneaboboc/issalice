// este modificat cum trebuie
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TCut.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TSpline.h>
#include <TTree.h>

#include "HistogramManager.h"
#include "ParticleFilter.h"
#include "Track.h"
//#include "./treeClass/Include.h"

#include <fstream>
#include <iostream>
#include <map>

using namespace std;

Double_t GetPhi(Double_t fX, Double_t fY) {
  return fX == 0.0 && fY == 0.0 ? 0.0 : TMath::ATan2(fY, fX);
}

Double_t GetEta(Double_t fX, Double_t fY, Double_t fZ) {
  Double_t sq_roots = TMath::Sqrt(fX * fX + fY * fY + fZ * fZ);
  if (sq_roots == 0) {
    return 0.0;
  }
  return 0.5 * TMath::Log((sq_roots + fZ) / (sq_roots - fZ));
}

Double_t GetEnergy(Double_t m, Double_t px, Double_t py, Double_t pz) {
  return TMath::Sqrt(TMath::Power(m, 2) + TMath::Power(px, 2) +
                     TMath::Power(py, 2) + TMath::Power(pz, 2));
}

/*
Int_t sarcina(Int_t pdg) {
  Int_t charge = 0;

  switch (pdg) {
    case 12:  // electron
      charge = -1;
      break;
    case -12:  // electron
      charge = 1;
      break;
    case 14:  // muon
      charge = -1;
      break;
    case -14:  // muon
      charge = 1;
      break;
    case 16:  // tau
      charge = -1;
      break;
    case -16:  // tau
      charge = 1;
      break;
    case 120:  // pi+
      charge = 1;
      break;
    case -120:  // pi-
      charge = -1;
      break;
    case 130:  // K+
      charge = 1;
      break;
    case -130:  // K-
      charge = -1;
      break;
    case 131:  // K*+
      charge = 1;
      break;
    case -131:  // K*-
      charge = -1;
      break;
    case 451:  // B*c+
      charge = 1;
      break;
    case -451:
      charge = -1;
      break;
    case 1120:  // proton
      charge = 1;
      break;
    case -1120:  // antiproton
      charge = -1;
      break;
    case 1130:  // Sigma+
      charge = 1;
      break;
    case -1130:  // anti Sigma+
      charge = -1;
      break;
    case 2230:  // Sigma-
      charge = -1;
      break;
    case -2230:  // anti Sigma-
      charge = 1;
      break;
    case 2330:  // Xi-
      charge = -1;
      break;
    case -2330:  // Xi-
      charge = +1;
      break;
    case 2331:  // Xi*- trecut
      charge = -1;
      break;
    case 3140:  // Xi(c)+
      charge = 1;
      break;
    case 1340:  // Xi'(c)+
      charge = 1;
      break;
    case 3331:  // Omega-
      charge = -1;
      break;
    case -3331:  // Omega+
      charge = 1;
      break;
    case 2140:  // Lamda(c)
      charge = 1;
      break;
    default:
      charge = 0;
  };

   Particule neutre gasite in lista de id
case 20: //Kshort - (pdg:310) trecut
  charge = 0;
  break;
case -20: //Klong  trecut
  charge = 0;
  break;
case 110: //pi0
  charge = 0;
   break;
case 1330: //Xi0
  charge = 0;
  break;
case 1331: //Xi*0 trecut
  charge = 0;
  break;
case 231: //K*0
  charge = 0;
  break;
case -231: //K*0b
  charge = 0;
  break;
case 3240: //Xi(c)0
  charge = 0;
  break;
case 2340: //Xi'(c)0
  charge = 0;
  break;
case 2130: //Lambda
  charge = 0;
  break;
case -2130: //anti-Lambda
  charge = 0;
  break;

  return charge;
}
  */

void print_map(const std::map<int, int> m) {
  for (const auto &[key, value] : m)
    std::cout << "code : charge << [" << key << "] = " << value << "; ";
  cout << '\n';
}

TChain *CreateChainLocal(Int_t nFilesMax, const Char_t *filename,
                         const Char_t *treeName);

void Ana(const Char_t *inFileName, const Char_t *outFileName,
         const Char_t *inFileSpName, Bool_t useChain) {
  TFile *inSp = TFile::Open(inFileSpName);
  if (!inSp) {
    cout << "Calibration file for spline does not exist!" << endl;
    return;
  }

  TSpline3 *sp3MultV0A = (TSpline3 *)inSp->Get("spMultV0A");

  TTree *inTree = 0;
  if (useChain) {
    cout << "Running analysis on a chain of files" << endl;
    inTree = CreateChainLocal(0, inFileName, "teposevent");
    if (!inTree) {
      cout << "Chain: " << inFileName << " does not exist!" << endl;
      return;
    }
  } else {
    cout << "Running analysis on a single file" << endl;
    TFile *inFile = TFile::Open(inFileName);
    if (!inFile) {
      cout << "File: " << inFileName << " does not exist!" << endl;
      return;
    }
    inTree = (TTree *)inFile->Get("teposevent");
  }

  std::map<int, int> m{
      {12, -1},   {-12, 1},    {14, -1},   {-14, 1},    {16, -1},  {-16, 1},
      {120, 1},   {-120, -1},  {130, 1},   {-130, -1},  {131, 1},  {-131, -1},
      {451, 1},   {-451, -1},  {1120, 1},  {-1120, -1}, {1130, 1}, {-1130, -1},
      {2230, 1},  {-2230, -1}, {2331, -1}, {3140, 1},   {1340, 1}, {3331, -1},
      {-3331, 1}, {2140, 1}};

  cout << "print charge of 120 " << m[120] << endl;

  print_map(m);

  Int_t centr = 10;

  Int_t nParticleFilters = 6;

  ParticleFilter particlefilter[nParticleFilters];

  particlefilter[0].addCondition(2, 130, 130, 131);  // conditie K+

  particlefilter[0].setName("K+");

  particlefilter[0].setTitle("K+");

  particlefilter[1].addCondition(2, -130, -130, -131);  // conditie K+

  particlefilter[1].setName("K-");

  particlefilter[1].setTitle("K-");

  particlefilter[2].addCondition(2, 120, 120, 121);  // conditie pi+

  particlefilter[2].setName("pi+");

  particlefilter[2].setTitle("pi+");

  particlefilter[3].addCondition(2, -120, -120, -121);  // conditie pi+

  particlefilter[3].setName("pi-");

  particlefilter[3].setTitle("pi-");

  particlefilter[4].addCondition(2, 1120, 1120, 1121);  // conditie p

  particlefilter[4].setName("p");

  particlefilter[4].setTitle("p");

  particlefilter[5].addCondition(2, -1120, -1120, -1121);  // conditie pbar

  particlefilter[5].setName("pBar");

  particlefilter[5].setTitle("pBar");

  TH1I *h_events[centr];

  for (Int_t i = 0; i < centr; i++) {
    h_events[i] = new TH1I(Form("h_events_%d", i), "; counts; No.", 1, 0, 1);
  }

  HistogramManager histogramManager;

  for (int i = 0; i < nParticleFilters; i++) {
    for (int j = 0; j < centr; j++) {
      TString name(particlefilter[i].getName());  // particlefilter[0].getName()
      name.Append("_" + to_string(j));

      ParticleSingleHistos *single_histos = new ParticleSingleHistos(name);

      // ParticleSingleHistos single_histos(name);
      // new(&single_histos[j]) ParticleSingleHistos(name);

      single_histos->createHistograms();

      histogramManager.addHistoInSet(i * centr + j, single_histos);
    }
  }

  for (int j = 0; j < nParticleFilters; j++) {
    for (int k = 0; k < nParticleFilters; k++) {
      for (int l = 0; l < centr; l++) {
        TString pairname =
            particlefilter[j].getName() + particlefilter[k].getName();
        pairname.Append("_" + to_string(l));

        ParticlePairHistos *pairHisto = new ParticlePairHistos(pairname);
        pairHisto->createHistograms();
        //        PairDerivedHistos* derivedHisto = new
        //        PairDerivedHistos(Form("name%d",j* nParticleFilters +
        //        k));
        PairDerivedHistos *derivedHisto = new PairDerivedHistos(pairname);
        derivedHisto->createHistograms();

        histogramManager.addHistoInSet(centr * (j * nParticleFilters + k) + l,
                                       pairHisto);
        histogramManager.addHistoInSet(centr * (j * nParticleFilters + k) + l,
                                       derivedHisto);
      }
    }
  }

  Float_t px[5000], py[5000], pz[5000], e[5000];
  Int_t id[5000], ist[5000], np;

  inTree->SetBranchAddress("px", &px);
  inTree->SetBranchAddress("py", &py);
  inTree->SetBranchAddress("pz", &pz);
  inTree->SetBranchAddress("e", &e);
  inTree->SetBranchAddress("id", &id);
  inTree->SetBranchAddress("ist", &ist);
  inTree->SetBranchAddress("np", &np);

  Int_t nEv = inTree->GetEntries();
  cout << "Number of events: " << nEv << endl;

  for (Int_t n = 0; n < nEv; n++) {  // nEv

    if ((n + 1) % 100000 == 0)
      cout << "Analysis - Event: " << n + 1 << "/" << nEv << endl;

    inTree->GetEntry(n);

    std::vector<Track> track;

    Int_t nAddTrk = 0;

    for (Int_t jj = 0; jj < np; jj++) {
      if (ist[jj] != 0) continue;

      double pt = TMath::Sqrt(px[jj] * px[jj] + py[jj] * py[jj]);
      int charge = 0;

      if (m.find(id[jj]) != m.end()) {
        charge = m[id[jj]];
      }

      else {
        // cout<<id[jj]<<endl;
        charge = 0;
      }
      double phiNew = GetPhi(px[jj], py[jj]);  //+ 3.1415;
      double trkEta = GetEta(px[jj], py[jj], pz[jj]);
      int pdg = id[jj];
      Double_t trkY =
          0.5 * TMath::Log((GetEnergy(e[jj], px[jj], py[jj], pz[jj]) + pz[jj]) /
                           (GetEnergy(e[jj], px[jj], py[jj], pz[jj]) - pz[jj]));
      // cout<<"trkEta "<<trkEta<<" trkY "<<trkY<<endl;

      if (charge == 0) continue;

      if (pt < 0.2 || pt > 2.0) continue;

      // if(trkEta < -1.0 || trkEta > 1.0) continue; //!fillY &&

      if (trkY < -1.0 || trkY > 1.0) continue;

      // TrackN* track = new((*fPartArray)[nAddTrk]) TrackN();
      Track *trk = new Track();

      trk->pdgCode = pdg;
      trk->phi = phiNew;
      trk->px = px[jj];
      trk->py = py[jj];
      trk->pz = pz[jj];
      trk->e = GetEnergy(e[jj], px[jj], py[jj], pz[jj]);
      trk->status = ist[jj];
      trk->eta = trkEta;
      trk->q = charge;

      track.push_back(*trk);
      nAddTrk++;

      delete trk;
    }

    if (nAddTrk < 2) continue;
    // cout<<"nr trackuri dsa:"<< nAddTrk <<endl;
    Double_t percMultV0A = 100. * (1. - sp3MultV0A->Eval(nAddTrk));

    if (percMultV0A < 0) {
      cout << "Negative percentile: put it to 0" << endl;
      percMultV0A = 0;
    }

    // cout<<"val spline: "<<percMultV0A<<endl;

    // percentile_dist->Fill(percMultV0A);
    // fill percentile distribution

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

    Int_t nTracks = nAddTrk;

    for (Int_t jl = 0; jl < nTracks; jl++) {
      for (int i = 0; i < nParticleFilters; i++) {
        if (particlefilter[i].accept(track[jl])) {
          histogramManager.HistogramSet[i * centr + multp]->fill(track[jl],
                                                                 1.0);

          for (Int_t jk = 0; jk < nTracks; jk++) {
            for (int j = 0; j < nParticleFilters; j++) {
              if (particlefilter[j].accept(track[jk])) {
                if (jl == jk) continue;

                histogramManager
                    .HistogramPairSet[centr * (i * nParticleFilters + j) +
                                      multp]
                    ->fill(track[jl], track[jk], 1.0);

              }  // filter2 accept

            }  // filters2

          }  // track2

        }  // filter1 accept

      }  // filters

    }  // tracks1

    track.clear();

  }  // event loop

  TCanvas *c = new TCanvas();

  TString strExt(outFileName);

  // strExt.Append("_"+to_string(centr));
  strExt.Append(".root");

  cout << "File output name: " << strExt << endl;
  cout << "---------- End of Analysis ----------" << endl;

  TFile *outfile = new TFile(strExt, "RECREATE");

  histogramManager.Write();

  for (Int_t i = 0; i < centr; i++) {
    h_events[i]->Write();
  }

  outfile->Close();
}

//__________________________________________________________
TChain *CreateChainLocal(Int_t nFilesMax, const Char_t *filename,
                         const Char_t *treeName) {
  TChain *chain = new TChain(treeName);

  // Open the input stream
  ifstream in;
  in.open(filename);

  Int_t nFiles = 0;

  // Read the input list of files and add them to the chain
  TString file;
  while (in.good() && (nFiles < nFilesMax || nFilesMax == 0)) {
    in >> file;
    if (!file.Contains("root")) continue;  // protection

    nFiles++;
    chain->Add(file.Data());
  }

  in.close();

  return chain;
}

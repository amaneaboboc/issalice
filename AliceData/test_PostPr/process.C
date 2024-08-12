#include "HistogramManager.h"
#include "ParticleFilter.h"
#include "histogramming.h"

void process()
{
    
    TFile* f =  new TFile("out.root");
    
    if(f->IsZombie()) {
        cout<<"check input file"<<endl;
        return;
    }
    
    Int_t centr = 7;

    Int_t nParticleFilters = 2;

    Int_t nBinsQ2 = 0;

    ParticleFilter particlefilter[nParticleFilters];

    particlefilter[0].addCondition(2, 130, 130, 131);  // conditie K+

    particlefilter[0].setName("n+");

    particlefilter[0].setTitle("n+");

    particlefilter[1].addCondition(2, -130, -130, -131);  // conditie K+

    particlefilter[1].setName("n-");

    particlefilter[1].setTitle("n-");
    
    TH1I *h_events = (TH1I*)f->Get("h_events");
    //TH1I* h_events[centr][10];
        
   // for(Int_t i = 0; i < centr;i++){
   //     for(Int_t j = 0; j < 10 ; j++){
   //     h_events[i][j] = (TH1I*)f->Get(Form("h_events_%d_%d",i,j));
   //     }
  //}

    
    HistogramManager histogramManager;

    for (int i = 0; i < nParticleFilters; i++) {
   //   for (int j = 0; j < centr; j++) {
        TString name(particlefilter[i].getName());  // particlefilter[0].getName()
        //name.Append("_" + to_string(j));
        //name.Append("_" + to_string(nBinsQ2));

        ParticleSingleHistos *single_histos = new ParticleSingleHistos(name);

        // ParticleSingleHistos single_histos(name);
        // new(&single_histos[j]) ParticleSingleHistos(name);

        //name.Append("_n1_phiEta");
        single_histos->readHistograms(f, name);

//          single_histos->h_n1_phiY->Draw();
          
        //histogramManager.addHistoInSet(i * centr + j, single_histos);
        histogramManager.addHistoInSet(i, single_histos);
     // }
    }
    
    for (int j = 0; j < nParticleFilters; j++) {
      for (int k = 0; k < nParticleFilters; k++) {
       // for (int l = 0; l < centr; l++) {
          TString pairname =
              particlefilter[j].getName() + particlefilter[k].getName();
          //pairname.Append("_" + to_string(l));
          //pairname.Append("_" + to_string(nBinsQ2));

          PairDerivedHistos *derivedHisto = new PairDerivedHistos(pairname);
          derivedHisto->createHistograms();
            
          //pairname.Append("_n2_DetaDphi");
          ParticlePairHistos *pairHisto = new ParticlePairHistos(pairname);
          pairHisto->readHistograms(f,pairname);
          //        PairDerivedHistos* derivedHisto = new
          //        PairDerivedHistos(Form("name%d",j* nParticleFilters +
          //        k));
          
          //histogramManager.addHistoInSet(centr * (j * nParticleFilters + k) + l,pairHisto);
          //histogramManager.addHistoInSet(centr * (j * nParticleFilters + k) + l,derivedHisto);


          histogramManager.addHistoInSet(j * nParticleFilters + k,
                                         pairHisto);
          histogramManager.addHistoInSet(j * nParticleFilters + k,
                                         derivedHisto);
       // }
      }
    }
    

    for(Int_t i = 0; i< nParticleFilters; i++){

    //for(Int_t k = 0; k<centr; k++){

      Double_t evts = h_events->GetEntries();

//    Int_t evts = histogramManager.HistogramSet[k]->GetEvents();
      histogramManager.HistogramSet[i]->Scale(evts);
        
    for(Int_t j = 0; j <nParticleFilters; j++){

    Double_t evts = h_events->GetEntries();
    //histogramManager.HistogramPairSet[centr * (i*nParticleFilters + j) + k]->Scale(1000);
        histogramManager.HistogramPairSet[i*nParticleFilters + j]->Scale(evts);

          //  }

        }

    }

    std::vector<ParticleSingleHistos*> vect = histogramManager.get_sgHist();
    std::vector<ParticlePairHistos*> vect_pair = histogramManager.get_pairHist();

    for( Int_t ii = 0; ii < nParticleFilters; ii++){

        for(Int_t ij = 0; ij < nParticleFilters; ij++){

       // for(Int_t k = 0; k<centr; k++){

       // histogramManager.HistogramDerivedSet[centr*(ii*nParticleFilters + ij) + k]->calculatePairDerivedHistograms(vect[ii*centr+k],vect[ij*centr+k],vect_pair[centr * (ii*nParticleFilters + ij) + k],1.0);

        histogramManager.HistogramDerivedSet[ii*nParticleFilters + ij]->calculatePairDerivedHistograms(vect[ii],vect[ij],vect_pair[ii*nParticleFilters+ij],0.0);//,vect_pair[ij*nParticleFilters+ii]
      //  }
        }

    }


    cout << "---------- End of Analysis ----------" << endl;

    TFile *outfile = new TFile("hist_PbPb_ESE_022.root", "RECREATE");

    histogramManager.Write();

    for (Int_t i = 0; i < centr; i++) {
      //h_events[i]->Write();
    }

    outfile->Close();
    
}

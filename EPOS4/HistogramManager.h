#include "TH1.h"
#include "TMath.h"
#include "Task.h"
#include "Track.h"
#include "histogramming.h"

class ParticleSingleHistos {
 public:
  ////////////////////////////////////////////////////////////////////////////
  // Data Members - HistogramGroup
  ////////////////////////////////////////////////////////////////////////////
  TString name;

  bool fillEta;
  bool fillY;
  bool fillP2;
  bool useEffCorrection;
  int efficiencyOpt;

  unsigned int nBins_n1;
  float min_n1;
  float max_n1;
  unsigned int nBins_pt;
  float min_pt;
  float max_pt;
  float scale_pt;
  unsigned int nBins_phi;
  float min_phi;
  float max_phi;
  float scale_phi;
  unsigned int nBins_eta;
  float min_eta;
  float max_eta;
  float range_eta;
  unsigned int nBins_y;
  float min_y;
  float max_y;
  float range_y;
  unsigned int nBins_phiEta;
  unsigned int nBins_phiEtaPt;
  unsigned int nBins_phiY;
  unsigned int nBins_phiYPt;

  //! Primary histograms

  TH1D *h_n1;
  TH1D *h_n1_pt;
  // TH1 * h_n1_ptXS;  // 1/pt dN/dptdy

  TH2D *h_n1_phiEta;
  // TH2 * h_spt_phiEta;

  TH2D *h_n1_phiY;
  // TH2 * h_spt_phiY;

  TH1D *h_events;

  ParticleSingleHistos(const TString &_name)
      :  // HistogramGroup(_name),
        fillEta(1),
        fillY(0),
        nBins_n1(100),
        min_n1(0),
        max_n1(1000),
        nBins_pt(100),
        min_pt(0.2),
        max_pt(2.0),
        nBins_phi(72),
        min_phi(0),
        max_phi(TMath::TwoPi()),
        nBins_eta(20),
        min_eta(-1.0),
        max_eta(1.0),
        range_eta(0),
        nBins_y(0),
        min_y(0),
        max_y(0),
        range_y(0),
        h_n1(nullptr),
        // h_n1_eTotal(nullptr),
        h_n1_pt(nullptr),
        // h_n1_ptXS(nullptr),
        h_n1_phiEta(nullptr),
        h_n1_phiY(nullptr),
        // h_spt_phiEta(nullptr),
        // h_n1_phiY(nullptr),
        // h_spt_phiY(nullptr),
        // h_pdgId(nullptr)
        h_events(nullptr) {
    // appendClassName("ParticleSingleHistos");
    name = _name;
  }

  ~ParticleSingleHistos() {
    // deleteHistograms();
  }

  TString getName() {
    if (1)  // Conditie filtru
      return name;
    else
      return "UnknownType";
  }

  void createHistograms() {
    const TString &bn = getName();

    nBins_n1 = 100;
    min_n1 = 0;
    max_n1 = 1000;

    nBins_pt = 100;
    min_pt = 0.2;
    max_pt = 2.0;
    scale_pt = max_pt - min_pt;

    nBins_phi = 72;
    min_phi = 0.0;
    max_phi = TMath::TwoPi();
    scale_phi = max_phi - min_phi;

    nBins_eta = 20;
    min_eta = -1.0;
    max_eta = 1.0;
    range_eta = max_eta - min_eta;

    nBins_y = 20;
    min_y = -1.0;
    max_y = 1.0;
    range_y = max_y - min_y;

    //        nBins_y = configuration.getValueInt(ppn,"nBins_y");
    //        min_y   = configuration.getValueDouble(ppn,"Min_y");
    //        max_y   = configuration.getValueDouble(ppn,"Max_y");
    //        range_y = max_y - min_y;

    fillEta = false;
    fillY = true;
    fillP2 = false;

    h_n1 = createHistogram(createName(bn, "n1"), nBins_n1, min_n1, max_n1,
                           "n_1", "N");

    h_events =
        createHistogram(createName(bn, "events"), 1, 0, 1, "No. Events", "bin");

    if (fillEta) {
      h_n1_phiEta = createHistogram(createName(bn, "n1_phiEta"), nBins_eta,
                                    min_eta, max_eta, nBins_phi, min_phi,
                                    max_phi, "#eta", "#varphi", "N");
    }

    if (fillY) {
      h_n1_phiY =
          createHistogram(createName(bn, "n1_phiY"), nBins_y, min_y, max_y,
                          nBins_phi, min_phi, max_phi, "y", "#varphi", "N");
    }
  }

  void fill(Track trk, double weight) {
    double pt = TMath::Sqrt(trk.px * trk.px + trk.py * trk.py);
    int charge = trk.q;

    double phiNew = trk.phi + TMath::Pi();
    double trkEta = trk.eta;
    int pdg = trk.pdgCode;
    double trkY = 0.5 * log((trk.e + trk.pz) / (trk.e - trk.pz));

    // std::cout<<pt<<endl;

    // here fill the needed histograms for 1 particle filters
    h_n1->Fill(pt);

    // h_n1_phiEta->Fill(trkEta,phiNew);
    h_n1_phiY->Fill(trkY, phiNew);
  }

  // once per event fill multiplicity of the global histogram
  void fillMultiplicity(double nAccepted, double weight) {
    h_n1->Fill(nAccepted, weight);
  }

  void fillEvents(double nAccepted) { h_events->Fill(nAccepted); }

  void toWrite() {
    h_n1->Write();
    // h_n1_phiEta->Write();
    h_n1_phiY->Write();
    // h_events->Write();
  }

  Int_t GetEvents() { return h_events->GetEntries(); }

  void Scale(Int_t nEv) {
    Double_t norm;

    // cout<<"nr of events "<<nEv<<endl;

    norm = 1 / Double_t(nEv);

    // h_n1_phiEta->Scale(norm);
    h_n1_phiY->Scale(norm);
  }
};

// ------------------Particle Pair Histograms -----------------

class ParticlePairHistos {
 public:
  TString name;

  bool fillEta;
  bool fillY;
  bool fillP2;
  bool useEffCorrection;
  int efficiencyOpt;

  unsigned int nBins_n2;
  double min_n2;
  double max_n2;
  unsigned int nBins_pt;
  double min_pt;
  double max_pt;
  double scale_pt;
  unsigned int nBins_phi;
  double min_phi;
  double max_phi;
  double range_phi;
  double scale_phi;
  unsigned int nBins_eta;
  double min_eta;
  double max_eta;
  double range_eta;
  double scale_eta;
  double etabinwidth;
  unsigned int nBins_y;
  double min_y;
  double max_y;
  double range_y;
  double scale_y;
  unsigned int nBins_phiEta;
  unsigned int nBins_phiEtaPt;
  unsigned int nBins_phiY;
  unsigned int nBins_phiYPt;
  unsigned int nBins_Deta;
  double min_Deta;
  double max_Deta;
  double range_Deta;
  unsigned int nBins_Dy;
  double min_Dy;
  double max_Dy;
  double range_Dy;
  unsigned int nBins_Dphi;
  double min_Dphi;
  double max_Dphi;
  double range_Dphi;
  double width_Dphi;
  unsigned int nBins_Dphi_shft;
  double min_Dphi_shft;
  double max_Dphi_shft;

  //! Primary histograms

  TH1D *h_n2;
  // TH2D* h_n2_Dpt;
  // TH2D* h_n1_ptXS;  // 1/pt dN/dptdy

  TH2D *h_n2_DetaDphi;
  TH2D *h_n2_DyDphi;
  // TH2 * h_spt_phiEta;

  ParticlePairHistos(const TString &_name)
      :  // HistogramGroup(_name),
        fillEta(1),
        fillY(0),
        nBins_n2(100),
        min_n2(0),
        max_n2(1000),
        nBins_pt(100),
        min_pt(0.2),
        max_pt(2.0),
        nBins_phi(72),
        min_phi(0),
        max_phi(TMath::TwoPi()),
        range_phi(0),
        scale_phi(0),
        nBins_eta(20),
        min_eta(-1.0),
        max_eta(1.0),
        range_eta(0),
        scale_eta(0),
        etabinwidth(0),
        nBins_y(0),
        min_y(0),
        max_y(0),
        range_y(0),
        scale_y(0),
        nBins_Deta(0),
        min_Deta(0),
        max_Deta(0),
        range_Deta(0),
        nBins_Dy(0),
        min_Dy(0),
        max_Dy(0),
        nBins_Dphi(0),
        min_Dphi(0),
        max_Dphi(0),
        range_Dphi(0),
        width_Dphi(0),
        min_Dphi_shft(0),
        max_Dphi_shft(0),
        h_n2(nullptr),
        h_n2_DetaDphi(nullptr) {
    // appendClassName("ParticleSingleHistos");
    name = _name;
  }

  void createHistograms() {
    const TString &bn = getName();

    nBins_n2 = 1000;
    min_n2 = 0;
    max_n2 = 1000;

    nBins_pt = 100;
    min_pt = 0.2;
    max_pt = 2.0;
    scale_pt = max_pt - min_pt;

    nBins_phi = 72;
    min_phi = 0.0;
    max_phi = TMath::TwoPi();
    range_phi = max_phi - min_phi;
    scale_phi = double(nBins_phi) / range_phi;
    width_Dphi = range_phi / double(nBins_phi);

    nBins_eta = 20;
    min_eta = -1.0;
    max_eta = 1.0;
    range_eta = max_eta - min_eta;
    scale_eta = double(nBins_eta) / range_eta;

    etabinwidth = range_eta / double(nBins_eta);

    nBins_y = 20;
    min_y = -1.0;
    max_y = 1.0;
    range_y = max_y - min_y;
    scale_y = double(nBins_y) / range_y;

    //      nBins_y = configuration.getValueInt(ppn,"nBins_y");
    //      min_y   = configuration.getValueDouble(ppn,"Min_y");
    //      max_y   = configuration.getValueDouble(ppn,"Max_y");
    //      range_y = max_y - min_y;

    fillEta = true;
    fillY = true;
    fillP2 = false;

    nBins_Deta = 2 * nBins_eta - 1;
    min_Deta = -range_eta;
    max_Deta = range_eta;

    nBins_Dy = 2 * nBins_y - 1;
    min_Dy = -range_y;
    max_Dy = range_y;

    nBins_Dphi = nBins_phi;
    min_Dphi = 0.0;             //-width_Dphi/2.;
    max_Dphi = TMath::TwoPi();  // - width_Dphi/2.;
    nBins_Dphi_shft = int(double(nBins_Dphi) / 4.0);
    min_Dphi_shft = min_Dphi - range_Dphi * double(nBins_Dphi_shft);
    max_Dphi_shft = max_Dphi - range_Dphi * double(nBins_Dphi_shft);

    h_n2 = createHistogram(createName(bn, "n2"), nBins_n2, min_n2, max_n2,
                           "n_{2}", "Yield");

    h_n2_DetaDphi = createHistogram(
        createName(bn, "n2_DetaDphi"), nBins_Deta, min_Deta, max_Deta,
        nBins_Dphi, min_Dphi, max_Dphi, "#Delta#eta", "#Delta#phi", "N_{2}");

    if (fillY) {
      h_n2_DyDphi = createHistogram(
          createName(bn, "n2_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#phi", "N_{2}");
    }
  }

  /*
    inline int getPtBinFor(float v) const
    {
    int index = 0; // indicates a value out of bounds
    if (v<min_pt || v>=max_pt) return index;
    index = 1+int(scale_pt*(v-min_pt));
    return index;
    }
  */
  inline int getPhiBinFor(float v) const {
    int index = 0;  // indicates a value out of bounds
    if (v < min_phi || v >= max_phi) return index;
    index = 1 + int(scale_phi * (v - min_phi));
    return index;
  }

  inline int getEtaBinFor(float v) const {
    int index = 0;  // indicates a value out of bounds
    if (v < min_eta || v >= max_eta) return index;
    index = 1 + int(scale_eta * (v - min_eta));
    return index;
  }

  inline int getYBinFor(float v) const {
    int index = 0;  // indicates a value out of bounds
    if (v < min_y || v >= max_y) return index;
    index = 1 + int(scale_y * (v - min_y));
    return index;
  }

  void fill(Track trk1, Track trk2, double weight) {
    int iPhi1, iPt1, iEta1, iY1;
    int iPhi2, iPt2, iEta2, iY2;
    // Double_t phi1,phi2;

    int iGDeltaEtaDeltaPhi;
    int iGDeltaYDeltaPhi;

    double phi1 = trk1.phi;
    if (phi1 < 0.0) {
      phi1 += TMath::TwoPi();
    }
    iPhi1 = getPhiBinFor(phi1);

    double eta1 = trk1.eta;
    iEta1 = fillEta ? getEtaBinFor(eta1) : 0;

    double y1 = 0.5 * log((trk1.e + trk1.pz) / (trk1.e - trk1.pz));
    iY1 = fillY ? getYBinFor(y1) : 0;

    double phi2 = trk2.phi;
    if (phi2 < 0.0) {
      phi2 += TMath::TwoPi();
    }
    iPhi2 = getPhiBinFor(phi2);

    double eta2 = trk2.eta;
    iEta2 = fillEta ? getEtaBinFor(eta2) : 0;

    double y2 = 0.5 * log((trk2.e + trk2.pz) / (trk2.e - trk2.pz));
    iY2 = fillY ? getEtaBinFor(y2) : 0;

    // iEta2 = getEtaBinFor(trk2->eta);//static_cast<int>(( (trk2->eta) -
    // min_eta)/etabinwidth); if(trk2->phi <= 0){phi2=trk2->phi;phi2
    // += 6.2831853;} iPhi2 = static_cast<int>((trk2->phi+TMath::Pi()) *
    // double(nBins_phi)/range_phi);

    // cout<<"phi2 dupa: "<<phi2<<endl;
    // if (iPt1==0  || iPt2==0)  return;
    if (iPhi1 == 0 || iPhi1 == 0) return;
    if (iEta1 == 0 && iY1 == 0) return;
    if (iEta2 == 0 && iY2 == 0) return;

    int iDeltaEta = iEta1 - iEta2 + nBins_eta - 1;
    int iDeltaPhi = iPhi1 - iPhi2;
    if (iDeltaPhi < 0) iDeltaPhi += nBins_phi;
    int iDeltaY = iY1 - iY2 + nBins_y - 1;

    // iGDeltaEtaDeltaPhi = h_n2_DetaDphi->GetBin(iDeltaEta+1,iDeltaPhi+1);
    // h_n2_DetaDphi->AddBinContent(iGDeltaEtaDeltaPhi,weight);
    // h_n2_DetaDphi->SetEntries(h_n2_DetaDphi->GetEntries()+1);

    iGDeltaYDeltaPhi = h_n2_DyDphi->GetBin(iDeltaY + 1, iDeltaPhi + 1);
    h_n2_DyDphi->AddBinContent(iGDeltaYDeltaPhi, weight);
    h_n2_DyDphi->SetEntries(h_n2_DyDphi->GetEntries() + 1);
  }

  void toWrite() {
    h_n2->Write();
    // h_n2_DetaDphi->Write();
    h_n2_DyDphi->Write();
  }

  void Scale(Int_t nEv) {
    Double_t norm;

    cout << "nr of events " << nEv << endl;

    norm = 1 / Double_t(nEv);

    // h_n2_DetaDphi->Scale(norm);
    h_n2_DyDphi->Scale(norm);
  }

  TString getName() {
    if (1)  // Conditie filtru
      return name;
    else
      return "UnknownType";
  }
};

//--------------Pair Derived Histograms-----------

class PairDerivedHistos {
 public:
  TString name;

  bool fillEta;
  bool fillY;
  bool fillP2;
  bool useEffCorrection;
  int efficiencyOpt;

  unsigned int nBins_n2;
  double min_n2;
  double max_n2;
  unsigned int nBins_pt;
  double min_pt;
  double max_pt;
  double scale_pt;
  unsigned int nBins_phi;
  double min_phi;
  double max_phi;
  double range_phi;
  unsigned int nBins_eta;
  double min_eta;
  double max_eta;
  double range_eta;
  double etabinwidth;
  unsigned int nBins_y;
  double min_y;
  double max_y;
  double range_y;
  double scale_y;
  unsigned int nBins_phiEta;
  unsigned int nBins_phiEtaPt;
  unsigned int nBins_phiY;
  unsigned int nBins_phiYPt;
  unsigned int nBins_Deta;
  double min_Deta;
  double max_Deta;
  double range_Deta;
  unsigned int nBins_Dy;
  double min_Dy;
  double max_Dy;
  double range_Dy;
  unsigned int nBins_Dphi;
  double min_Dphi;
  double max_Dphi;
  double range_Dphi;
  double width_Dphi;
  unsigned int nBins_Dphi_shft;
  double min_Dphi_shft;
  double max_Dphi_shft;

  // TH2D* h_n2_Dpt;
  // TH2D* h_n1_ptXS;  // 1/pt dN/dptdy

  TH2D *h_n1n1_DetaDphi;
  TH2D *h_rho2_DetaDphi;
  TH2D *h_A2_DetaDphi;
  TH2D *h_B2_DetaDphi;
  TH2D *h_C2_DetaDphi;
  TH2D *h_R2_DetaDphi;

  TH2D *h_n1n1_DetaDphi_shft;
  TH2D *h_rho2_DetaDphi_shft;
  TH2D *h_A2_DetaDphi_shft;
  TH2D *h_B2_DetaDphi_shft;
  TH2D *h_C2_DetaDphi_shft;
  TH2D *h_R2_DetaDphi_shft;
  // TH2 * h_spt_phiEta;

  TH2D *h_n1n1_DyDphi;
  TH2D *h_rho2_DyDphi;
  TH2D *h_A2_DyDphi;
  TH2D *h_B2_DyDphi;
  TH2D *h_C2_DyDphi;
  TH2D *h_R2_DyDphi;

  TH2D *h_rho2_DyDphi_shft;
  TH2D *h_A2_DyDphi_shft;
  TH2D *h_B2_DyDphi_shft;
  TH2D *h_C2_DyDphi_shft;
  TH2D *h_R2_DyDphi_shft;

  PairDerivedHistos(const TString &_name)
      :  // HistogramGroup(_name),
        /*  fillEta(1),
      fillY(0),
      nBins_n2(100),
      min_n2(0),
      max_n2(1000),
      nBins_pt(100),
      min_pt(0.2),
      max_pt(2.0),
      nBins_phi(72),
      min_phi(0),
      max_phi(TMath::TwoPi()),
      range_phi(0),
      nBins_eta(20),
      min_eta(-1.0),
      max_eta(1.0),
      range_eta(0),
      etabinwidth(0),
      nBins_y(0),
      min_y(0),
      max_y(0),
      nBins_Deta(0),
      min_Deta(0),
      max_Deta(0),
      range_Deta(0),
      nBins_Dphi(0),
      min_Dphi(0),
      max_Dphi(0),
      range_Dphi(0),
      min_Dphi_shft(0),
      max_Dphi_shft(0),
      range_y(0),*/
        h_n1n1_DetaDphi(nullptr),
        h_rho2_DetaDphi(nullptr),
        h_A2_DetaDphi(nullptr),
        h_B2_DetaDphi(nullptr),
        h_C2_DetaDphi(nullptr),
        h_R2_DetaDphi(nullptr) {
    // appendClassName("ParticleSingleHistos");
    name = _name;
  }

  //! Primary histograms

  void createHistograms() {
    const TString &bn = getName();

    nBins_n2 = 1000;
    min_n2 = 0;
    max_n2 = 1000;

    nBins_pt = 100;
    min_pt = 0.2;
    max_pt = 2.0;
    scale_pt = max_pt - min_pt;

    nBins_phi = 72;
    min_phi = 0.0;
    max_phi = TMath::TwoPi();
    range_phi = max_phi - min_phi;

    nBins_eta = 20;
    min_eta = -1.0;
    max_eta = 1.0;
    range_eta = max_eta - min_eta;

    nBins_y = 20;
    min_y = -1.0;
    max_y = 1.0;
    range_y = max_y - min_y;
    scale_y = double(nBins_y) / range_y;

    etabinwidth = range_eta / double(nBins_eta);

    double scale_phi = max_phi - min_phi;
    double width_Dphi = scale_phi / nBins_phi;

    fillEta = false;
    fillY = true;
    // fillP2  = false;

    nBins_Deta = 2 * nBins_eta - 1;
    min_Deta = -range_eta;
    max_Deta = range_eta;

    nBins_Dy = 2 * nBins_y - 1;
    min_Dy = -range_y;
    max_Dy = range_y;

    nBins_Dphi = nBins_phi;
    min_Dphi = 0.0;             //-width_Dphi/2.;
    max_Dphi = TMath::TwoPi();  // - width_Dphi/2.;
    nBins_Dphi_shft = int(double(nBins_Dphi) / 4.0);

    width_Dphi = range_phi / double(nBins_phi);

    min_Dphi_shft = min_Dphi - width_Dphi * double(nBins_Dphi_shft);
    max_Dphi_shft = max_Dphi - width_Dphi * double(nBins_Dphi_shft);

    if (fillEta) {
      h_n1n1_DetaDphi =
          createHistogram(createName(bn, "n1n1_DetaDphi"), nBins_Deta, min_Deta,
                          max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
                          "#Delta#eta", "#Delta#varphi", "<n_{1}><n_{1}>");
      h_rho2_DetaDphi =
          createHistogram(createName(bn, "rho2_DetaDphi"), nBins_Deta, min_Deta,
                          max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
                          "#Delta#eta", "#Delta#varphi", "#rho_{2}");
      h_A2_DetaDphi =
          createHistogram(createName(bn, "A2_DetaDphi"), nBins_Deta, min_Deta,
                          max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
                          "#Delta#eta", "#Delta#varphi", "A_{2}");
      h_B2_DetaDphi =
          createHistogram(createName(bn, "B2_DetaDphi"), nBins_Deta, min_Deta,
                          max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
                          "#Delta#eta", "#Delta#varphi", "B_{2}");
      h_C2_DetaDphi =
          createHistogram(createName(bn, "C2_DetaDphi"), nBins_Deta, min_Deta,
                          max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
                          "#Delta#eta", "#Delta#varphi", "C_{2}");
      // h_D2_DetaDphi        = createHistogram(createName(bn,"D2_DetaDphi"),
      // nBins_Deta, min_Deta, max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
      // "#Delta#eta","#Delta#varphi", "D_{2}");
      h_R2_DetaDphi =
          createHistogram(createName(bn, "R2_DetaDphi"), nBins_Deta, min_Deta,
                          max_Deta, nBins_Dphi, min_Dphi, max_Dphi,
                          "#Delta#eta", "#Delta#varphi", "R_{2}");

      h_rho2_DetaDphi_shft = createHistogram(
          createName(bn, "rho2_DetaDphi_shft"), nBins_Deta, min_Deta, max_Deta,
          nBins_Dphi, min_Dphi_shft, max_Dphi_shft, "#Delta#eta",
          "#Delta#varphi", "#rho_{2}");
      h_A2_DetaDphi_shft = createHistogram(
          createName(bn, "A2_DetaDphi_shft"), nBins_Deta, min_Deta, max_Deta,
          nBins_Dphi, min_Dphi_shft, max_Dphi_shft, "#Delta#eta",
          "#Delta#varphi", "A_{2}");
      h_B2_DetaDphi_shft = createHistogram(
          createName(bn, "B2_DetaDphi_shft"), nBins_Deta, min_Deta, max_Deta,
          nBins_Dphi, min_Dphi_shft, max_Dphi_shft, "#Delta#eta",
          "#Delta#varphi", "B_{2}");
      h_C2_DetaDphi_shft = createHistogram(
          createName(bn, "C2_DetaDphi_shft"), nBins_Deta, min_Deta, max_Deta,
          nBins_Dphi, min_Dphi_shft, max_Dphi_shft, "#Delta#eta",
          "#Delta#varphi", "C_{2}");
      // h_D2_DetaDphi_shft   =
      // createHistogram(createName(bn,"D2_DetaDphi_shft"),    nBins_Deta,
      // min_Deta,  max_Deta,  nBins_Dphi,  min_Dphi_shft, max_Dphi_shft,
      // "#Delta#eta","#Delta#varphi", "D_{2}");
      h_R2_DetaDphi_shft = createHistogram(
          createName(bn, "R2_DetaDphi_shft"), nBins_Deta, min_Deta, max_Deta,
          nBins_Dphi, min_Dphi_shft, max_Dphi_shft, "#Delta#eta",
          "#Delta#varphi", "R_{2}");
    }

    if (fillY) {
      h_n1n1_DyDphi = createHistogram(
          createName(bn, "n1n1_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#varphi", "<n_{1}><n_{1}>");
      h_rho2_DyDphi = createHistogram(
          createName(bn, "rho2_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#varphi", "#rho_{2}>");
      h_A2_DyDphi = createHistogram(
          createName(bn, "A2_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#varphi", "A_{2}");
      h_B2_DyDphi = createHistogram(
          createName(bn, "B2_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#varphi", "B_{2}");
      h_C2_DyDphi = createHistogram(
          createName(bn, "C2_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#varphi", "C_{2}");
      h_R2_DyDphi = createHistogram(
          createName(bn, "R2_DyDphi"), nBins_Dy, min_Dy, max_Dy, nBins_Dphi,
          min_Dphi, max_Dphi, "#Delta y", "#Delta#varphi", "R_{2}");

      h_rho2_DyDphi_shft =
          createHistogram(createName(bn, "rho2_DyDphi_shft"), nBins_Dy, min_Dy,
                          max_Dy, nBins_Dphi, min_Dphi_shft, max_Dphi_shft,
                          "#Delta y", "#Delta#varphi", "#rho_{2}>");
      h_A2_DyDphi_shft =
          createHistogram(createName(bn, "A2_DyDphi_shft"), nBins_Dy, min_Dy,
                          max_Dy, nBins_Dphi, min_Dphi_shft, max_Dphi_shft,
                          "#Delta y", "#Delta#varphi", "A_{2}");
      h_B2_DyDphi_shft =
          createHistogram(createName(bn, "B2_DyDphi_shft"), nBins_Dy, min_Dy,
                          max_Dy, nBins_Dphi, min_Dphi_shft, max_Dphi_shft,
                          "#Delta y", "#Delta#varphi", "B_{2}");
      h_C2_DyDphi_shft =
          createHistogram(createName(bn, "C2_DyDphi_shft"), nBins_Dy, min_Dy,
                          max_Dy, nBins_Dphi, min_Dphi_shft, max_Dphi_shft,
                          "#Delta y", "#Delta#varphi", "C_{2}");
      h_R2_DyDphi_shft =
          createHistogram(createName(bn, "R2_DyDphi_shft"), nBins_Dy, min_Dy,
                          max_Dy, nBins_Dphi, min_Dphi_shft, max_Dphi_shft,
                          "#Delta y", "#Delta#varphi", "R_{2}");
    }
  }

  void calculatePairDerivedHistograms(ParticleSingleHistos *part1BaseHistos,
                                      ParticleSingleHistos *part2BaseHistos,
                                      ParticlePairHistos *pairHistos,
                                      double bincorrection) {
    double yield2;
    // int ijNormalization = 0;
    if (fillEta) {
      yield2 = part2BaseHistos->h_n1_phiEta->Integral();

      cout << "yield2: " << yield2 << endl;

      reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi(
          part1BaseHistos->h_n1_phiEta, part2BaseHistos->h_n1_phiEta,
          h_n1n1_DetaDphi, nBins_Deta, nBins_Dphi);
      h_rho2_DetaDphi->Add(pairHistos->h_n2_DetaDphi);

      h_C2_DetaDphi->Add(h_rho2_DetaDphi, h_n1n1_DetaDphi, 1.0, -1.0);
      h_B2_DetaDphi->Add(h_rho2_DetaDphi);
      h_B2_DetaDphi->Scale(1.0 / yield2);
      h_A2_DetaDphi->Add(h_C2_DetaDphi);
      h_A2_DetaDphi->Scale(1.0 / yield2);
      calculateR2_H2H2H2(h_rho2_DetaDphi, h_n1n1_DetaDphi, h_R2_DetaDphi, 0,
                         1.0, 1.0);

      shiftY(*h_rho2_DetaDphi, *h_rho2_DetaDphi_shft, nBins_Dphi_shft);
      shiftY(*h_A2_DetaDphi, *h_A2_DetaDphi_shft, nBins_Dphi_shft);
      shiftY(*h_B2_DetaDphi, *h_B2_DetaDphi_shft, nBins_Dphi_shft);
      shiftY(*h_C2_DetaDphi, *h_C2_DetaDphi_shft, nBins_Dphi_shft);
      shiftY(*h_R2_DetaDphi, *h_R2_DetaDphi_shft, nBins_Dphi_shft);
    }

    double yield2_y;
    if (fillY) {
      yield2_y = part2BaseHistos->h_n1_phiY->Integral();

      cout << "yield2_y: " << yield2 << endl;

      reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi(
          part1BaseHistos->h_n1_phiY, part2BaseHistos->h_n1_phiY, h_n1n1_DyDphi,
          nBins_Dy, nBins_Dphi);
      h_rho2_DyDphi->Add(pairHistos->h_n2_DyDphi);
      h_C2_DyDphi->Add(h_rho2_DyDphi, h_n1n1_DyDphi, 1.0, -1.0);
      h_B2_DyDphi->Add(h_rho2_DyDphi);
      h_B2_DyDphi->Scale(1.0 / yield2_y);
      h_A2_DyDphi->Add(h_C2_DyDphi);
      h_A2_DyDphi->Scale(1.0 / yield2_y);
      calculateR2_H2H2H2(h_rho2_DyDphi, h_n1n1_DyDphi, h_R2_DyDphi, 0, 1.0,
                         1.0);

      shiftY(*h_rho2_DyDphi, *h_rho2_DyDphi_shft, nBins_Dphi_shft);
      shiftY(*h_A2_DyDphi, *h_A2_DyDphi_shft, nBins_Dphi_shft);
      shiftY(*h_B2_DyDphi, *h_B2_DyDphi_shft, nBins_Dphi_shft);
      shiftY(*h_C2_DyDphi, *h_C2_DyDphi_shft, nBins_Dphi_shft);
      shiftY(*h_R2_DyDphi, *h_R2_DyDphi_shft, nBins_Dphi_shft);
    }
    /*
    TCanvas* r1 = new TCanvas();

    part1BaseHistos->h_n1_phiEta->Draw("surf1,fb");

    //TCanvas* r = new TCanvas();

    //h_rho2_DetaDphi->Draw("surf1,fb");

    TCanvas* cd = new TCanvas();

    h_C2_DetaDphi->Draw("surf1,fb");

    TCanvas* ci = new TCanvas();

    h_A2_DetaDphi_shft->Draw("surf1,fb");
*/
  }

  void toWrite() {
    // h_n1n1_DetaDphi->Write();
    // h_rho2_DetaDphi->Write();
    // h_R2_DetaDphi_shft->Write();
    // h_C2_DetaDphi_shft->Write();
    // h_B2_DetaDphi_shft->Write();
    // h_A2_DetaDphi_shft->Write();

    // h_R2_DyDphi_shft->Write();
    // h_C2_DyDphi_shft->Write();
    // h_B2_DyDphi_shft->Write();
    // h_A2_DyDphi_shft->Write();
  }

  TString getName() {
    if (1)  // Conditie filtru
      return name;
    else
      return "UnknownType";
  }
};

// Histogram manager--------------------------------------

class HistogramManager {
 public:
  vector<ParticleSingleHistos *> HistogramSet;
  vector<ParticlePairHistos *> HistogramPairSet;
  vector<PairDerivedHistos *> HistogramDerivedSet;

  vector<TString> names;
  // vector<HistogramManager> sets;

  void addHistoInSet(unsigned int index, ParticleSingleHistos *group) {
    unsigned int nSets = HistogramSet.size();
    if (index >= nSets) {
      TString s("");
      s += int(index);
    }
    HistogramSet.push_back(group);
  }

  void addHistoInSet(unsigned int index, ParticlePairHistos *group) {
    unsigned int nSets = HistogramSet.size();
    if (index >= nSets) {
      TString s("");
      s += int(index);
    }
    HistogramPairSet.push_back(group);
  }

  void addHistoInSet(unsigned int index, PairDerivedHistos *group) {
    unsigned int nSets = HistogramSet.size();
    if (index >= nSets) {
      TString s("");
      s += int(index);
    }
    HistogramDerivedSet.push_back(group);
  }

  void addSet(const TString &name) {
    if (name.Length() != 0) {
      names.push_back(name);
    }
  }

  void getHistoInSet(unsigned int index, ParticleSingleHistos *hist) {
    hist = HistogramSet[index];
  }

  void Write() {
    int nSets = HistogramSet.size();

    for (int i = 0; i < nSets; i++) {
      HistogramSet[i]->toWrite();
    }

    int nPairSets = HistogramPairSet.size();

    for (int i = 0; i < nPairSets; i++) {
      HistogramPairSet[i]->toWrite();

      HistogramDerivedSet[i]->toWrite();
    }
  }

  const std::vector<ParticleSingleHistos *> &get_sgHist() const {
    return HistogramSet;
  }

  const std::vector<ParticlePairHistos *> &get_pairHist() const {
    return HistogramPairSet;
  }
};

// Histogram group----------------------------------------

class HistogramGroup : HistogramManager {};

#ifndef histogramming_h
#define histogramming_h
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"

TH1D *createHistogram(const TString &name, int n, double min_x, double max_x,
                      const TString &title_x, const TString &title_y) {
  // cout << "Creating  1D histo " << name << " nBins:" << n << " min_x:" <<
  // min_x << " max_x:" << max_x << endl;
  TH1D *h = new TH1D(name, name, n, min_x, max_x);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  return h;
}

TH1D *createHistogram(const TString &name, int n, double *bins,
                      const TString &title_x, const TString &title_y) {
  // std::cout << "Creating  1D histo " << name << " with " << n << " non
  // uniform nBins:" << endl;
  TH1D *h = new TH1D(name, name, n, bins);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  return h;
}

TH2D *createHistogram(const TString &name, int n_x, double min_x, double max_x,
                      int n_y, double min_y, double max_y,
                      const TString &title_x, const TString &title_y,
                      const TString &title_z) {
  // std::cout << "Creating  2D histo " << name << " n_x:" << n_x << " min_x:"
  // << min_x << " max_x:" << max_x << " n_y:" << n_y << " min_y:" << min_y <<
  // " max_y:" << max_y<< endl;
  TH2D *h;
  h = new TH2D(name, name, n_x, min_x, max_x, n_y, min_y, max_y);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  if (title_z.Sizeof() > 0) h->GetZaxis()->SetTitle(title_z);
  return h;
}

TH2D *createHistogram(const TString &name, int n_x, double *xbins, int n_y,
                      double min_y, double max_y, const TString &title_x,
                      const TString &title_y, const TString &title_z)

{
  // std::cout << "Creating  2D histo " << name << " with " << n_x << " vs " <<
  // n_y << " non uniform nBins:" << endl;
  TH2D *h;
  h = new TH2D(name, name, n_x, xbins, n_y, min_y, max_y);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  if (title_z.Sizeof() > 0) h->GetZaxis()->SetTitle(title_z);
  return h;
}

//!
// Create 3D histogram
//!
TH3D *createHistogram(const TString &name, int n_x, double min_x, double max_x,
                      int n_y, double min_y, double max_y, int n_z,
                      double min_z, double max_z, const TString &title_x,
                      const TString &title_y, const TString &title_z,
                      const TString &title_w = "w")

{
  // std::cout << "Creating  3D histo " << name
  //   << " n_x:" << n_x << " min_x:" << min_x << " max_x:" << max_x << " X:"
  //   << title_x
  // << " n_y:" << n_y << " min_y:" << min_y << " max_y:" << max_y << " Y:" <<
  // title_y
  // << " n_z:" << n_z << " min_z:" << min_z << " max_z:" << max_z << " Z:" <<
  // title_z
  // << " W:" << title_w  << endl;
  TH3D *h;
  h = new TH3D(name, name, n_x, min_x, max_x, n_y, min_y, max_y, n_z, min_z,
               max_z);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  if (title_z.Sizeof() > 0) h->GetZaxis()->SetTitle(title_z);
  return h;
}

//!
// Create 1D profile histogram
//!
TProfile *createProfile(const TString &name, int n_x, double min_x,
                        double max_x, const TString &title_x,
                        const TString &title_y) {
  // cout << "Creating  1D profile " << name << " n_x:" << n_x << " min_x:" <<
  // min_x << " max_x:" << max_x << endl;
  TProfile *h = new TProfile(name, name, n_x, min_x, max_x);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  return h;
}

//!
// Create 1D profile histogram
//!
TProfile *createProfile(const TString &name, int n_x, double *bins,
                        const TString &title_x, const TString &title_y) {
  // cout << "Creating  1D profile " << name << " w/ variabe size bins" <<
  // endl;
  TProfile *h = new TProfile(name, name, n_x, bins);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  return h;
}

//!
// Create 2D profile histogram
//!
TProfile2D *createProfile(const TString &name, int n_x, double min_x,
                          double max_x, int n_y, double min_y, double max_y,
                          const TString &title_x, const TString &title_y,
                          const TString &title_z) {
  // cout << "Creating  2D profile " << name
  // << " n_x:" << n_x << " min_x:" << min_x << " max_x:" << max_x
  // << " n_y:" << n_y << " min_y:" << min_y << " max_y:" << max_y
  // << endl;
  TProfile2D *h =
      new TProfile2D(name, name, n_x, min_x, max_x, n_y, min_y, max_y);
  if (title_x.Sizeof() > 0) h->GetXaxis()->SetTitle(title_x);
  if (title_y.Sizeof() > 0) h->GetYaxis()->SetTitle(title_y);
  if (title_z.Sizeof() > 0) h->GetZaxis()->SetTitle(title_z);
  return h;
}

const TString createName(const TString &s0, const TString &s1, int option = 0) {
  TString separator;
  switch (option) {
    default:
    case 0:
      separator = "_";
      break;
    case 1:
      separator = " ";
      break;
    case 2:
      separator = "";
      break;
  }
  TString name = s0;
  name += separator;
  name += s1;
  return name;
}

const TString createName(const TString &s0, const TString &s1,
                         const TString &s2, int option = 0) {
  TString separator;
  switch (option) {
    default:
    case 0:
      separator = "_";
      break;
    case 1:
      separator = " ";
      break;
    case 2:
      separator = "";
      break;
  }
  TString name = s0;
  name += separator;
  name += s1;
  name += s2;
  return name;
}

const TString createName(const TString &s0, const TString &s1,
                         const TString &s2, const TString &s3, int option = 0) {
  TString separator;
  switch (option) {
    default:
    case 0:
      separator = "_";
      break;
    case 1:
      separator = " ";
      break;
    case 2:
      separator = "";
      break;
  }
  TString name = s0;
  name += separator;
  name += s1;
  name += separator;
  name += s2;
  name += separator;
  name += s3;

  return name;
}

const TString createName(const TString &s0, const TString &s1,
                         const TString &s2, const TString &s3,
                         const TString &s4, int option = 0) {
  TString separator;
  switch (option) {
    default:
    case 0:
      separator = "_";
      break;
    case 1:
      separator = " ";
      break;
    case 2:
      separator = "";
      break;
  }
  TString name = s0;
  name += separator;
  name += s1;
  name += separator;
  name += s2;
  name += separator;
  name += s3;
  name += separator;
  name += s4;
  return name;
}

TH1D * loadH1(TFile  inputFile,const TString & histoName)
{
    
    TH1D* h = (TH1D*)inputFile.Get(histoName);
    
    if(!h)
    {
        //cout<<"Histo not found"<<'\n';
    }
    return h;
    
}

TH2D * loadH2(TFile * inputFile,const TString & histoName)
{
    TH2D* h = (TH2D*)inputFile->Get(histoName);
 
    if(!h)
    {
        //cout<<"Histo not found"<<'\n';
    }
    
    return h;

}



#endif

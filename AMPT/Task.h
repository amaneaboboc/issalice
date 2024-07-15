

//!
//! Calculate the external product of n_1(eta,phi) by n_2(eta,phi) and project
//! onto n1n1_12(Deta,Dphi) First create a 4D external product of h_1 and h_2
//! using a 1D vector than project it onto DeltaEta(x-axis) DeltaPhi (y-axis)
//! Errors are neglected and must be accounted for using a sub-sample analysis.
//!
void
reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi (const TH2 *h_1, TH2 *h_2, TH2 *h_12,
                                         int nDeta, int nDphi)
{

  // if (reportStart(__FUNCTION__))
  //  ;
  // if (!ptrExist(__FUNCTION__,h_1,h_2,h_12)) return;
  // if (!sameDimensions(__FUNCTION__,h_1,h_2)) return;

  vector<double> numerator (nDeta * nDphi, 0.0);
  vector<double> denominator (nDeta * nDphi, 0.0);
  int nEta = h_1->GetNbinsX ();
  int nPhi = h_1->GetNbinsY ();
  nDeta = h_12->GetNbinsX ();
  nDphi = h_12->GetNbinsY ();

  int index;

  // if (reportDebug(__FUNCTION__))
  //   {
  /*
    cout << endl;
    cout << "         nEta:" << nEta << endl;
    cout << "         nPhi:" << nPhi << endl;
    cout << "        nDeta:" << nDeta << endl;
    cout << "        nDphi:" << nDphi << endl;
    cout << "  nDeta*nDphi:" << nDeta*nDphi << endl;
    */
  // }

  int iDeta, iDphi;
  double v1, v2, v, r;
  for (int iPhi1 = 1; iPhi1 <= nPhi; iPhi1++)
    {
      for (int iEta1 = 1; iEta1 <= nEta; iEta1++)
        {
          v1 = h_1->GetBinContent (iEta1, iPhi1);
          for (int iPhi2 = 1; iPhi2 <= nPhi; iPhi2++)
            {
              for (int iEta2 = 1; iEta2 <= nEta; iEta2++)
                {
                  v2 = h_2->GetBinContent (iEta2, iPhi2);
                  v = v1 * v2;
                  iDphi = iPhi1 - iPhi2;
                  if (iDphi < 0)
                    iDphi += nDphi; // dPhi+=1;
                  iDeta = iEta1 - iEta2 + nEta - 1;
                  // if (iDphi<0 || )
                  index = iDphi * nDeta + iDeta;
                  if (index < 0)
                    {

                      // cout << "<F>
                      // HistogramCollection::reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi()
                      // index<0" << endl;
                      // throws
                      // HistogramException("LOGIC","index<0","HistogramCollection::reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi()");
                    }
                  if (index >= nDeta * nDphi)
                    {
                      // cout << "<F>
                      // HistogramCollection::reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi()
                      // index>nDeta*nDphi" << endl;
                      // throw
                      // HistogramException("LOGIC","index>=nDeta*nDphi","HistogramCollection::reduce_n1EtaPhiN1EtaPhiOntoN1N1DetaDphi()");
                    }
                  numerator[index] += v; // cout << "numerator" <<
                                         // numerator[index]<< "  .. ";
                  denominator[index] += 1.0; // cout<<" denominator "<<
                                             // denominator[index]<<endl;
                }
            }
        }
    }
  // if (reportDebug(__FUNCTION__)) cout << "Compute ratio and fill target
  // histogram." << endl;
  double zero = 0.0;
  for (iDphi = 0; iDphi < nDphi; iDphi++)
    {
      for (iDeta = 0; iDeta < nDeta; iDeta++)
        {
          index = iDphi * nDeta + iDeta;
          double num = numerator[index];
          double den = denominator[index];
          // cout << " iDphi: " << iDphi << " iDeta:" << iDeta << " index:" <<
          // index << " num:" << num << " denom:" << den << endl; if (den>0) r
          // = num/den;
          if (den > 0)
            r = num;
          //          if(den>0) r = den;
          else
            r = zero;
          h_12->SetBinContent (iDeta + 1, iDphi + 1, r);
          h_12->SetBinError (iDeta + 1, iDphi + 1, zero);
        }
    }
  // if (reportEnd(__FUNCTION__))
  //  ;
}

// Calculate External Product n1_1 x n1_2 and store into n1n1_12
double
calculateN1N1_H1H1H2 (const TH1 *h_1, const TH1 *h_2, TH2 *h_12, double a1,
                      double a2)
{

  // if (reportStart(__FUNCTION__))
  ;
  // if (!ptrExist(__FUNCTION__,h_1,h_2,h_12)) return -1.0;
  int n1 = h_1->GetNbinsX ();
  int n2 = h_2->GetNbinsX ();
  int n3x = h_12->GetNbinsX ();
  int n3y = h_12->GetNbinsY ();

  if (n3x != n1 || n3y != n2)
    {

      // if (reportDebug(__FUNCTION__))
      //   {
      /*
        cout << endl;
        cout << "Incompatible histo dimensions" << endl;
        cout << "H1: " << h_1->GetName()    << " nBins:" << n1 << endl;
        cout << "H2: " << h_2->GetName()    << " nBins:" << n2 << endl;
        cout << "H3: " << h_12->GetName()   << " nBins_x:" << n3x << "
        nBins_y:" << n3y << endl;
        */
      return -1.0;
      //}
    }
  double v1, ev1, v2, ev2, v, ev, r1, r2;
  double sum = 0.;
  double norm = 0.;
  for (int i1 = 1; i1 <= n1; ++i1)
    {
      v1 = a1 * h_1->GetBinContent (i1);
      ev1 = a1 * h_1->GetBinError (i1);
      for (int i2 = 1; i2 <= n2; ++i2)
        {
          v2 = a2 * h_2->GetBinContent (i2);
          ev2 = a2 * h_2->GetBinError (i2);
          v = v1 * v2;
          if (v > 0)
            {
              r1 = ev1 / v1;
              r2 = ev2 / v2;
              ev = v * sqrt (r1 * r1 + r2 * r2);
            }
          else
            ev = 0.;
          h_12->SetBinContent (i1, i2, v);
          h_12->SetBinError (i1, i2, ev);
          sum += v;
          norm += 1.;
        }
    }
  // return average across bins
  return sum / norm;
}

// Calculate R2 = N2/N1/N1 - 1
void
calculateR2_H2H2H2 (const TH2 *n2_12, const TH2 *n1n1_12, TH2 *r2_12,
                    bool ijNormalization, double a1, double a2)
{

  //  if (reportStart(__FUNCTION__))
  //    ;

  //  if (!ptrExist(__FUNCTION__,n2_12,n1n1_12,r2_12)) return;
  //  if (!sameDimensions(__FUNCTION__,n2_12,n1n1_12,r2_12)) return;
  int n2_12_n_x = n2_12->GetNbinsX ();
  int n2_12_n_y = n2_12->GetNbinsY ();

  double v1, ev1, v2, ev2, v, ev, re1, re2;
  for (int i_x = 1; i_x <= n2_12_n_x; ++i_x)
    {
      for (int i_y = 1; i_y <= n2_12_n_y; ++i_y)
        {
          v1 = a1 * n2_12->GetBinContent (i_x, i_y);
          ev1 = a1 * n2_12->GetBinError (i_x, i_y);
          v2 = a2 * n1n1_12->GetBinContent (i_x, i_y);
          ev2 = a2 * n1n1_12->GetBinError (i_x, i_y);
          if (v1 > 0 && v2 > 0) //   && ev1/v1<0.5  && ev2/v2<0.5)
            {
              if (ijNormalization) // account for the fact only half the pairs
                                   // were counted
                v = 2 * v1 / v2;
              else // all pairs counted - no need to multiply by 2
                v = v1 / v2;
              re1 = ev1 / v1;
              re2 = ev2 / v2;
              ev = v * sqrt (re1 * re1 + re2 * re2);
              v -= 1.;
            }
          else
            {
              v = 0.;
              ev = 0;
            }
          r2_12->SetBinContent (i_x, i_y, v);
          r2_12->SetBinError (i_x, i_y, ev);
        }
    }
}

/// shift the given source to the target vertically by nbins
void
shiftY (const TH2 &source, TH2 &target, int nbins)
{

  // if (reportStart(__FUNCTION__))
  //  ;

  int i_x, i_y;
  int n_x = source.GetNbinsX ();
  int n_y = source.GetNbinsY ();

  // shift the 1st area
  for (i_x = 1; i_x <= n_x; ++i_x)
    {
      for (i_y = 1; i_y <= n_y - nbins; ++i_y)
        {
          double v = source.GetBinContent (i_x, i_y);
          double ev = source.GetBinError (i_x, i_y);
          target.SetBinContent (i_x, i_y + nbins, v);
          target.SetBinError (i_x, i_y + nbins, ev);
        }
      for (i_y = n_y - nbins + 1; i_y <= n_y; ++i_y)
        {
          double v = source.GetBinContent (i_x, i_y);
          double ev = source.GetBinError (i_x, i_y);
          target.SetBinContent (i_x, i_y - (n_y - nbins), v);
          target.SetBinError (i_x, i_y - (n_y - nbins), ev);
        }
    }
}

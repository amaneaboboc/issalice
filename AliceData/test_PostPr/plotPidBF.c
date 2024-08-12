#include "routine.h"
#include "makePads.C"
class BalFct
{
public:
    TH2D h_n1_2;
    TH2D h_A2_abBar;
    TH2D h_A2_aBarbBar;
    TH2D h_R2_abBar;
    TH2D h_R2_aBarbBar;
    
    TString str;
    TString str_cen;
    
    ~BalFct(){
        h_R2_abBar.Clear();
        h_R2_aBarbBar.Clear();
    }
    
    TString GetFile(TString & file){
        
        cout<<"fisierul de input "<<file<<endl;
        return str = file;
    }
    
    
    TH2* CalcBF(TString str_cen, TString str_p1, TString str_p2){
        
        TFile* inFile = new TFile(str);

        if(!inFile){
            cout<<"file not found"<<endl;
        }
        cout<<"str: "<<str<<endl;
                
        
        map<string, string> particlesData;
        particlesData.insert(pair<string, string>("pi+","pi-"));
        particlesData.insert(pair<string, string>("K+","K-"));
        particlesData.insert(pair<string, string>("p","pBar"));
        
        particlesData.insert(pair<string, string>("pi-","pi+"));
        particlesData.insert(pair<string, string>("K-","K+"));
        particlesData.insert(pair<string, string>("pBar","p"));
        
        string caca = (string)str_p1;
        string caca2 = (string)str_p2;
        
        TString str_antip1(particlesData[caca]);
        TString str_antip2(particlesData[caca2]);
        
        cout<<"return pair of pi+ "<< particlesData[caca]<<endl;

        TH2D* h_n1_2 = (TH2D*)inFile->Get(str_p1+"_"+str_cen+"_n1_phiY");

        TH2D* h_R2_abBar =(TH2D*)inFile->Get(str_p1+str_antip2+"_"+str_cen+"_R2_DyDphi_shft");
        TH2D* h_R2_aBarbBar =(TH2D*)inFile->Get(str_antip1+str_antip2+"_"+str_cen+"_R2_DyDphi_shft");
        h_R2_abBar->Add(h_R2_aBarbBar,-1.0);
        
        TH2D* h_R2_aBarb =(TH2D*)inFile->Get(str_antip1+str_p2+"_"+str_cen+"_R2_DyDphi_shft");
        TH2D* h_R2_ab =(TH2D*)inFile->Get(str_p1+str_p2+"_"+str_cen+"_R2_DyDphi_shft");
        
        h_R2_aBarb->Add(h_R2_aBarbBar,-1.0);
        
        //media B
        h_R2_abBar->Add(h_R2_aBarb,1.0);
        h_R2_abBar->Scale(0.5);
        
        double scalingfactor =  average_rho(h_n1_2);
        h_R2_abBar->Scale(scalingfactor);
        
        Double_t wx = h_R2_abBar->GetXaxis()->GetBinWidth(1);
        Double_t wy = h_R2_abBar->GetYaxis()->GetBinWidth(1);
        h_R2_abBar->Scale(1.0/(wx*wy));
        
        Double_t max_Deta = 2.0;

        TCutG *gcut = new TCutG("gcut",4);
        gcut->SetPoint(0,-max_Deta,-TMath::PiOver2());
        gcut->SetPoint(1,max_Deta,-TMath::PiOver2());
        gcut->SetPoint(2,max_Deta,TMath::PiOver2());
        gcut->SetPoint(3,-max_Deta,TMath::PiOver2());
//
        TH1* hist1_proj = h_R2_abBar->ProjectionY("R2_1_x");
        
        return h_R2_abBar;
        
    }
    
    TH1D* GetIntegrals(TString str_p1, TString str_p2){
        
        const Int_t nBinsPerc = 10;
        Double_t binsPerc[nBinsPerc+1] = {0.0,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};

        Double_t err_mult[nBinsPerc];
        Double_t integral_A2_plmi[nBinsPerc];
        Double_t integral_A2_away[nBinsPerc];
        Double_t integral_A2_near[nBinsPerc];
        
        Double_t mult_class[nBinsPerc];
        Double_t err_mult_A_plmi[nBinsPerc];
        
        for(Int_t j = 0; j < nBinsPerc; j++)
        {
          mult_class[j] = (binsPerc[j] + binsPerc[j+1])/2;
            //cout<< mult_class[j]<<endl;
            err_mult[j] = 0 ;
        }
        
        TFile* inFile = new TFile(str);

        if(!inFile){
            cout<<"file not found"<<endl;
        }
//        cout<<"str: "<<str<<endl;
        for(Int_t i = 0; i < nBinsPerc;i++){
        
        //TString str_p1("K+");// adauga pentru alta particula
        
        //TString str_p2("K-");
            
            map<string, string> particlesData;
            particlesData.insert(pair<string, string>("pi+","pi-"));
            particlesData.insert(pair<string, string>("K+","K-"));
            particlesData.insert(pair<string, string>("p","pBar"));
            particlesData.insert(pair<string, string>("pi-","pi+"));
            particlesData.insert(pair<string, string>("K-","K+"));
            particlesData.insert(pair<string, string>("pBar","p"));
            
            string caca = (string)str_p1;
            string caca2 = (string)str_p2;
            
            TString str_antip1(particlesData[caca]);
            TString str_antip2(particlesData[caca2]);
            
        TString str_cen(to_string(i));
        
        TH2D* h_n1_2 = (TH2D*)inFile->Get(str_p1+"_"+str_cen+"_n1_phiY");

        TH2D* h_A2_abBar =(TH2D*)inFile->Get(str_p1+str_antip2+"_"+str_cen+"_R2_DyDphi_shft");
        TH2D* h_A2_aBarbBar =(TH2D*)inFile->Get(str_antip1+str_antip2+"_"+str_cen+"_R2_DyDphi_shft");
        h_A2_abBar->Add(h_A2_aBarbBar,-1.0);
            
        TH2D* h_n1_1 = (TH2D*)inFile->Get(str_p2+"_"+str_cen+"_n1_phiY");

        TH2D* h_A2_aBarb =(TH2D*)inFile->Get(str_antip1+str_p2+"_"+str_cen+"_R2_DyDphi_shft");
        TH2D* h_A2_ab =(TH2D*)inFile->Get(str_p1+str_p2+"_"+str_cen+"_R2_DyDphi_shft");
        h_A2_aBarb->Add(h_A2_ab,-1.0);
            
        double scalingfactor =  average_rho(h_n1_2);
        h_A2_abBar->Scale(scalingfactor);
            
        scalingfactor =  average_rho(h_n1_1);
        h_A2_aBarb->Scale(scalingfactor);
            
            h_A2_abBar->Add(h_A2_aBarb,1);
            h_A2_abBar->Scale(0.5);
            
        Double_t wx = h_A2_abBar->GetXaxis()->GetBinWidth(1);
        Double_t wy = h_A2_abBar->GetYaxis()->GetBinWidth(1);

//        h_A2_abBar->Scale(wx*wy);
            TH2D* h_A2_abBar_corr = (TH2D*) h_A2_abBar->Clone();

            acceptance(h_A2_abBar,h_A2_abBar_corr,1.0);
        
        integral_A2_plmi[i] = h_A2_abBar->IntegralAndError(1,h_A2_abBar->GetNbinsX(),1,h_A2_abBar->GetNbinsY(),err_mult_A_plmi[i],"");
    }
        
        TH1D* A2_int_all = new TH1D("test","test",nBinsPerc,binsPerc);
//= new TGraphErrors(nBinsPerc,mult_class,integral_A2_plmi,err_mult,err_mult_A_plmi);
        for(Int_t i=0; i< nBinsPerc;i++){
            
            A2_int_all->SetBinContent(i,integral_A2_plmi[i]);
            A2_int_all->SetBinError(i,err_mult_A_plmi[i]);
            
        }
//
//        TCanvas* inte = new TCanvas();
//
//        A2_int_all->Draw();
//
//        A2_int_all->GetXaxis()->SetTitle("Multiplicity class");
//        A2_int_all->GetYaxis()->SetTitle("Integral");
//
//        A2_int_all->SetLineColor(1);
//        A2_int_all->SetMarkerColor(1);
//        A2_int_all->SetMarkerStyle(20);
//        A2_int_all->Draw("AP");
//
//        A2_int_all->GetXaxis()->SetTitleSize(0.045);
//        A2_int_all->GetYaxis()->SetTitleSize(0.045);
//
//        A2_int_all->GetYaxis()->SetLabelSize(0.045);
//        A2_int_all->GetXaxis()->SetLabelSize(0.045);
        
        return A2_int_all;
    }
    
};

void PlotInt(TH1* hist){
    
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    
//    hist->GetXaxis()->SetRangeUser(-1.8,1.8);
    hist->SetMaximum(1.2);
    hist->SetMinimum(0.0);

    hist->GetYaxis()->SetTitle("I_{B}");

    gStyle->SetOptStat(0);
    hist->GetXaxis()->SetNdivisions(505, "X");
    hist->GetYaxis()->SetNdivisions(505, "Y");
    hist->GetZaxis()->SetNdivisions(505, "Z");
}

void Plotting2D(TH2* hist){
    
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetZaxis()->CenterTitle();

    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetZaxis()->SetTitleSize(0.06);
    
    hist->GetXaxis()->SetLabelSize(0.06);
    hist->GetYaxis()->SetLabelSize(0.06);
    hist->GetZaxis()->SetLabelSize(0.06);

    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetZaxis()->SetTitleOffset(1.2);

    hist->GetXaxis()->SetRangeUser(-1.8,1.8);
    
    hist->GetYaxis()->SetTitle("#Delta#varphi(rad)");

    hist->SetMinimum(0);
    hist->SetMaximum(0.115);

    gStyle->SetOptStat(0);
    hist->GetXaxis()->SetNdivisions(505, "X");
    hist->GetYaxis()->SetNdivisions(505, "Y");
    hist->GetZaxis()->SetNdivisions(505, "Z");
    
}

void SetStyle(Bool_t graypalette=kFALSE);


void plotPidBF(){
    
    //    SetStyle(kFALSE);
    
    TString inFileP("hist_Monash_022.root");
    TString inFileE("hist_Epos_022.root");
        
    BalFct EPOS1, EPOS2;
    BalFct PYTH1, PYTH2;
    
    EPOS1.GetFile(inFileE);
    EPOS2.GetFile(inFileE);
    PYTH1.GetFile(inFileP);
    PYTH2.GetFile(inFileP);
    
    TString cen_1("0");
    TString cen_2("4");
    TString cen_3("8");
    
    TString pion_pl("pi+");
    TString pion_mi("pi-");
    
    TString kaon_pl("K+");
    TString kaon_mi("K-");
    
    TString prot("p");
    TString anti_prot("pBar");
    
    TH2* hist_epos = EPOS1.CalcBF(cen_2,"pi+","pi+");
    TH2* hist_pyth = PYTH1.CalcBF(cen_2,"pi+","pi+");
    
    SetStyle();
    
    TCanvas* c = new TCanvas("canv", "canv");
    hist_epos->Draw("surf1,fb");
    
    hist_epos->GetZaxis()->SetTitle("B^{pp}(#Deltay,#Delta#varphi)");
    
    Plotting2D(hist_epos);
    
    c->SetTheta(30.);
    c->SetPhi(-60.);
    
    TLatex* text = new TLatex(0.09, 0.9, "#splitline{EPOS 4 p-p}{0.2 < p_{T} < 2.0 GeV/c}");
    text->SetNDC();
    text->SetTextSize(0.06);
    text->SetTextColor(kBlack);
    text->Draw();
    
    TLatex* textM = new TLatex(0.75, 0.85, "0-10%");
    textM->SetNDC();
    textM->SetTextSize(0.07);
    textM->SetTextColor(kBlack);
    textM->Draw();
    
    TCanvas* c1 = new TCanvas("canv_pyth", "canv_pyth");
    hist_pyth->Draw("surf1,fb");
    
    hist_pyth->GetZaxis()->SetTitle("B^{pp}(#Deltay,#Delta#varphi)");
    
    Plotting2D(hist_pyth);
    
    c1->SetTheta(30.);
    c1->SetPhi(-60.);
    
    TLatex* text1 = new TLatex(0.09, 0.9, "#splitline{PYTHIA 8 p-p}{0.2 < p_{T} < 2.0 GeV/c}");
    text1->SetNDC();
    text1->SetTextSize(0.06);
    text1->SetTextColor(kBlack);
    text1->Draw();
    
    textM->Draw();
    
    TH1D* h_int_pipi_pyth = PYTH1.GetIntegrals("pi+","pi+");
    TH1D* h_int_KK_pyth = PYTH1.GetIntegrals("K+","K+");
    TH1D* h_int_pp_pyth = PYTH1.GetIntegrals("p","p");
    
    TH1D* h_int_pipi_epos = EPOS1.GetIntegrals("pi+","pi+");
    TH1D* h_int_KK_epos = EPOS1.GetIntegrals("K+","K+");
    TH1D* h_int_pp_epos = EPOS1.GetIntegrals("p","p");
    
    TCanvas* c_int = new TCanvas("integral","integral",600,500);
    
    h_int_pipi_pyth->Draw("P");
    h_int_KK_pyth->Draw("same");
    h_int_pp_pyth->Draw("same");
    
    h_int_pipi_epos->Draw("same");
    h_int_KK_epos->Draw("same");
    h_int_pp_epos->Draw("same");
    
    h_int_pipi_epos->SetMarkerColor(kBlue);
    h_int_KK_epos->SetMarkerColor(kRed);
    h_int_pp_epos->SetMarkerColor(kBlack);
    h_int_pipi_epos->SetLineColor(kWhite);
    h_int_KK_epos->SetLineColor(kWhite);
    h_int_pp_epos->SetLineColor(kWhite);
    
    
    h_int_pipi_epos->SetMarkerStyle(20);
    h_int_KK_epos->SetMarkerStyle(20);
    h_int_pp_epos->SetMarkerStyle(20);
    
    h_int_pipi_pyth->SetMarkerColor(kBlue);
    h_int_KK_pyth->SetMarkerColor(kRed);
    h_int_pp_pyth->SetMarkerColor(kBlack);
    h_int_pipi_pyth->SetLineColor(kWhite);
    h_int_KK_pyth->SetLineColor(kWhite);
    h_int_pp_pyth->SetLineColor(kWhite);
    
    h_int_pipi_pyth->SetMarkerStyle(24);
    h_int_KK_pyth->SetMarkerStyle(24);
    h_int_pp_pyth->SetMarkerStyle(24);
    
    h_int_pipi_pyth->GetXaxis()->CenterTitle();
    h_int_pipi_pyth->GetYaxis()->CenterTitle();
    
    h_int_pipi_pyth->GetXaxis()->SetTitleSize(0.06);
    h_int_pipi_pyth->GetYaxis()->SetTitleSize(0.06);
    
    h_int_pipi_pyth->GetXaxis()->SetLabelSize(0.06);
    h_int_pipi_pyth->GetYaxis()->SetLabelSize(0.06);
    
    h_int_pipi_pyth->GetXaxis()->SetTitleOffset(1.1);
    h_int_pipi_pyth->GetYaxis()->SetTitleOffset(1.1);
    
    h_int_pipi_pyth->SetMaximum(1.2);
    h_int_pipi_pyth->SetMinimum(0.0);
    
    h_int_pipi_pyth->GetXaxis()->SetTitle("Multiplicity class(%)");
    h_int_pipi_pyth->GetYaxis()->SetTitle("I_{B}");
    
    h_int_pipi_pyth->GetXaxis()->SetNdivisions(505, "X");
    h_int_pipi_pyth->GetYaxis()->SetNdivisions(505, "Y");
    
    TLatex* tex = new TLatex(0.18, 0.75, "#splitline{0.2 < p_{T} < 2.0}{|y|<1}");
    tex->SetNDC();
    tex->SetTextSize(0.06);
    tex->SetTextColor(kBlack);
    tex->Draw();
    

    
    //    TFile* fout = new TFile("hist_Integrals.root","RECREATE");
    //    hist1->Write();
    //    hist2->Write();
    //    Int_pyth->Write();
    //    Int_epos->Write();
    //    Int_pyth_pi->Write();
    //    Int_epos_pi->Write();
    //    Int_pyth_p->Write();
    //    Int_epos_p->Write();
    //    fout->Close();
    
    
    TFile* fEpos_unId = new TFile("Int_unId_epos.root");
    TFile* fPyth_unId = new TFile("Int_unId_monash.root");

    TGraph* graph_Epos = (TGraph*)fEpos_unId->Get("Graph");
    TGraph* graph_Pyth = (TGraph*)fPyth_unId->Get("Graph");
    
    const Int_t nBins =7;
    Double_t nCent[nBins+1] = {0., 5., 10., 20.0, 40.0, 60., 80., 100.};
    
    TH1D* int_Epos = new TH1D("int_Epos",";Multiplicity class(%);I_{B}",nBins,nCent);
    TH1D* int_Pyth = new TH1D("int_Pyth",";Multiplicity class(%);I_{B}",nBins,nCent);// the histogram (you should set the number of bins, the title etc)
    auto nPoints = graph_Epos->GetN(); // number of points in your TGraph
    for(int i=0; i < nPoints; ++i) {
       double x,y;
        graph_Epos->GetPoint(i, x, y);
        int_Epos->SetBinContent(i,y);
        graph_Pyth->GetPoint(i, x, y);
        int_Pyth->SetBinContent(i,y);
    }

    TCanvas* c4 = new TCanvas("c1", "c1", 900, 600);
    TPad** pads = makePads("pads", 2, 1);
    
    TH2F* hdummy1 = new TH2F("hdummy1", ";Multiplicity class(%);B (2, |#Delta#eta| > 0.8)", 100, -0.1, 4.1, 100, -0.01, 0.21);
    TH2F* hdummy2 = new TH2F("hdummy2", ";Multiplicity class(%);B (2, |#Delta#eta| > 0.8)", 100, -0.1, 4.1, 100, -0.01, 0.31);
    
    hdummy1->SetStats(0);
    hdummy2->SetStats(0);
    
    pads[0]->cd();
    
    int_Pyth->Draw("P");
    PlotInt(int_Pyth);
    int_Epos->Draw("same");

    int_Pyth->SetMarkerStyle(20);
    int_Epos->SetMarkerStyle(20);

    pads[1]->cd();

    h_int_pipi_epos->Draw();
    PlotInt(h_int_pipi_epos);

    h_int_KK_epos->Draw("same");
    h_int_pp_epos->Draw("same");
    
    h_int_pipi_pyth->Draw("same");
    h_int_KK_pyth->Draw("same");
    h_int_pp_pyth->Draw("same");
    
    TLatex* tex1 = new TLatex(0.38, 0.85, "#splitline{0.2 < p_{T} < 2.0}{|y|<1}");
    tex1->SetNDC();
    tex1->SetTextSize(0.06);
    tex1->SetTextColor(kBlack);
    tex1->Draw();
    
    auto legend = new TLegend(0.67,0.67,0.88,0.89);
    legend->SetHeader("EPOS 4","C"); // option "C" allows to center the header
    legend->AddEntry(h_int_pipi_epos,"#pi#pi","p");
    legend->AddEntry(h_int_KK_epos,"KK","p");
    legend->AddEntry(h_int_pp_epos,"p+#bar{p}","p");
    legend->Draw();
    
    auto legend1 = new TLegend(0.47,0.67,0.68,0.89);
    legend1->SetHeader("PYTHIA 8","C"); // option "C" allows to center the header
    legend1->AddEntry(h_int_pipi_pyth,"#pi#pi","p");
    legend1->AddEntry(h_int_KK_pyth,"KK","p");
    legend1->AddEntry(h_int_pp_pyth,"p+#bar{p}","p");
    legend1->Draw();
    
//    int_Pyth->GetXaxis()->CenterTitle();
//    int_Pyth->GetYaxis()->CenterTitle();
    
//    int_Pyth->GetXaxis()->SetTitleSize(0.06);
//    int_Pyth->GetYaxis()->SetTitleSize(0.06);
    
//    int_Pyth->GetXaxis()->SetLabelSize(0.06);
//    int_Pyth->GetYaxis()->SetLabelSize(0.06);
    
//    hist->GetXaxis()->SetRangeUser(-1.8,1.8);
//    int_Pyth->SetMaximum(1.2);
//    int_Pyth->SetMinimum(0.0);

//    int_Pyth->GetXaxis()->SetNdivisions(505, "X");
//    int_Pyth->GetYaxis()->SetNdivisions(505, "Y");
    
}

//______________________________________________________________________
void SetStyle(Bool_t graypalette) {
    cout << "Setting style!" << endl;
    
    gStyle->Reset("Plain");
    gStyle->SetOptTitle(0);
    //gStyle->SetOptStat(0);
    if(graypalette) gStyle->SetPalette(8,0);
    else gStyle->SetPalette(kBird);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    //gStyle->SetPadTickY(1);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kGreen);
    //gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.65,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.09,"y");
    gStyle->SetTitleSize(0.15,"xz");
    gStyle->SetTitleOffset(1.25,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    gStyle->SetTitleW(0.8);
    //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");
    
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(kWhite);
    //  gStyle->SetFillColor(kWhite);
    gStyle->SetLegendFont(42);
    
    gStyle->SetNdivisions(505, "Z");
//    gStyle->SetNdivisions(505, "Z");
    gStyle->SetEndErrorSize(0);
    
}

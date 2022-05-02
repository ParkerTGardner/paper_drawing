//blh
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
//#include "TGraphMultiErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include <vector>
#include "math.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include <numeric>
#include "TLatex.h"
//#include "coordinateTools.h"

using TMath::ATan;
using TMath::Exp;


double VOut(TH2D* h1, TH2D* h2, int YPlo, int YPhi,int vn)
{

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlo,YPhi)->Clone();
    TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlo,YPhi)->Clone();
    histfit1->Divide(histfit2);
    histfit1->Scale(h2->GetMaximum());
    histfit1->SetMarkerStyle(20);
    std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";//+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";//+[6]*TMath::Cos(6*x)+[7]*TMath::Cos(7*x)))";
    TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
    func1.SetParameter(0, histfit1->GetMaximum());
    func1.SetParameter(1, 0.1);
    func1.SetParameter(2, 0.1);
    func1.SetParameter(3, 0.1);
    //func1.SetParameter(4, 0.1);
    //func1.SetParameter(5, 0.1);
    histfit1->Fit(&func1, "q 0");
    histfit1->Fit(&func1, "q 0");
    histfit1->Fit(&func1, "m q 0");
    histfit1->Fit(&func1, "m q 0");
    histfit1->Fit(&func1, "m q E 0");
    histfit1->Fit(&func1, "m E q 0");
    return func1.GetParameter(vn);
}
//GetParError
double VEOut(TH2D* h1, TH2D* h2, int YPlo, int YPhi, int vn)
{   

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlo,YPhi)->Clone();
    TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlo,YPhi)->Clone();
    histfit1->Divide(histfit2);
    histfit1->Scale(h2->GetMaximum());
    histfit1->SetMarkerStyle(20);
    std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";//+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";//+[6]*TMath::Cos(6*x)+[7]*TMath::Cos(7*x)))";
    TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
    func1.SetParameter(0, histfit1->GetMaximum());
    func1.SetParameter(1, 0.1);
    func1.SetParameter(2, 0.1);
    func1.SetParameter(3, 0.1);
    //func1.SetParameter(4, 0.1);
    //func1.SetParameter(5, 0.1);
    histfit1->Fit(&func1, "q 0");
    histfit1->Fit(&func1, "q 0");
    histfit1->Fit(&func1, "m q 0");
    histfit1->Fit(&func1, "m q 0");
    histfit1->Fit(&func1, "m q E 0");
    histfit1->Fit(&func1, "m E q 0");
    return func1.GetParError(vn);
}




void diff_reg_pu123()
{
//all_data_paper.root
    TFile *f_PU   = new TFile("all_PU.root");
    TFile *f_data = new TFile("all_data_paper.root");

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int trackbin 	= 9;

    const int ptbin 	= 1;

    const int   trackbinbounds[trackbin]         = { 0,20,30,40,50,60,68,74,80};
    const int   trackbinboundsUpper[trackbin]    = {20,30,40,50,60,68,74,80,1000};


    const float ptbinbounds_lo[ptbin]       = {3};
    const float ptbinbounds_hi[ptbin]       = {30};
    float BW1 = 0.3;
    int ptwant 	= 1;

    int pt_lo 	= ptbinbounds_lo[ptwant-1];
    int pt_hi 	= ptbinbounds_hi[ptwant-1];

    int Xscale 	= 800;
    int Yscale 	= 800;


    TH1D* hJ1_data;
    hJ1_data = (TH1D*)f_data->Get("hJet_Pass")->Clone("j1_reco_co");

    double reco_co_Lx_1[trackbin] = {1,1,1,1,1,1,1,1,1};
    double reco_co_Ly_1[trackbin];
    double reco_co_Lxel_1[trackbin] = {10,5,5,5,5,5,4,4,4};
    double reco_co_Lxeh_1[trackbin] = {10,5,5,5,5,5,4,4,10};
    double reco_co_Lye_1[trackbin];

    double reco_co_Lx_2[trackbin] = {1,1,1,1,1,1,1,1,1};//{10,25,35,45,55,65,74,82,90};
    double reco_co_Ly_2[trackbin];
    double reco_co_Lxel_2[trackbin] = {10,5,5,5,5,5,4,4,4};
    double reco_co_Lxeh_2[trackbin] = {10,5,5,5,5,5,4,4,10};
    double reco_co_Lye_2[trackbin];

    double reco_co_Lx_3[trackbin] = {1,1,1,1,1,1,1,1,1};//{10,25,35,45,55,65,74,82,90};
    double reco_co_Ly_3[trackbin];
    double reco_co_Lxel_3[trackbin] = {10,5,5,5,5,5,4,4,4};;
    double reco_co_Lxeh_3[trackbin] = {10,5,5,5,5,5,4,4,10};
    double reco_co_Lye_3[trackbin];

    double reco_co_cLx_1[trackbin];
    double reco_co_cLy_1[trackbin];
    double reco_co_cLxel_1[trackbin];
    double reco_co_cLxeh_1[trackbin];
    double reco_co_cLye_1[trackbin];

    double reco_co_cLx_2[trackbin];
    double reco_co_cLy_2[trackbin];
    double reco_co_cLxel_2[trackbin];
    double reco_co_cLxeh_2[trackbin];
    double reco_co_cLye_2[trackbin];

    double reco_co_cLx_3[trackbin];
    double reco_co_cLy_3[trackbin];
    double reco_co_cLxel_3[trackbin];
    double reco_co_cLxeh_3[trackbin];
    double reco_co_cLye_3[trackbin];

    for(int wtrk = 0; wtrk < trackbin; wtrk++){
        //hSigS_60_to_70_and_3_to_30;1
        TH2D* h1;
        h1 = (TH2D*)f_data->Get(Form("hSigS_Cor_%d_to_%d_and_%d_to_%d_w_PU_1",trackbinbounds[wtrk], trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        h1->Scale(1/(hJ1_data->GetBinContent(wtrk+1)));
        h1->Scale(1./(BW1));

        TH2D* h2;
        h2 = (TH2D*)f_data->Get(Form("hBckS_Cor_%d_to_%d_and_%d_to_%d_w_PU_1",trackbinbounds[wtrk], trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();

        reco_co_Ly_1[wtrk]  = VOut(h1,h2,28,41, 1);       
        reco_co_Lye_1[wtrk] = VEOut(h1,h2,28,41, 1);       

        reco_co_Ly_2[wtrk]  = VOut(h1,h2,28,41, 2);       
        reco_co_Lye_2[wtrk] = VEOut(h1,h2,28,41, 2);       

        reco_co_Ly_3[wtrk]  = VOut(h1,h2,28,41, 3);       
        reco_co_Lye_3[wtrk] = VEOut(h1,h2,28,41, 3);       

    }


    for (int i = 0; i < trackbin; i++) {
        reco_co_cLx_1[i] = reco_co_Lx_1[i];
        reco_co_cLy_1[i] = reco_co_Ly_1[i];
        reco_co_cLxel_1[i] = reco_co_Lxel_1[i];
        reco_co_cLxeh_1[i] = reco_co_Lxeh_1[i];
        reco_co_cLye_1[i] = (sqrt (2) )*(reco_co_Lye_1[i]);
    }
    for (int i = 0; i < trackbin; i++) {
        reco_co_cLx_2[i] = reco_co_Lx_2[i];
        reco_co_cLy_2[i] = reco_co_Ly_2[i];
        reco_co_cLxel_2[i] = reco_co_Lxel_2[i];
        reco_co_cLxeh_2[i] = reco_co_Lxeh_2[i];
        reco_co_cLye_2[i] = (sqrt (2) )*(reco_co_Lye_2[i]);
    }
    for (int i = 0; i < trackbin; i++) {
        reco_co_cLx_3[i] = reco_co_Lx_3[i];
        reco_co_cLy_3[i] = reco_co_Ly_3[i];
        reco_co_cLxel_3[i] = reco_co_Lxel_3[i];
        reco_co_cLxeh_3[i] = reco_co_Lxeh_3[i];
        reco_co_cLye_3[i] = (sqrt (2) )*(reco_co_Lye_3[i]);
    }


    // NUMBER 2 NUMBER 2 NUMBER 2 **************************************************************************************
    // NUMBER 2 NUMBER 2 NUMBER 2 **************************************************************************************
    // NUMBER 2 NUMBER 2 NUMBER 2 **************************************************************************************
    // NUMBER 2 NUMBER 2 NUMBER 2 **************************************************************************************


    TH1D* hJ1_f2_PU;
    hJ1_f2_PU = (TH1D*)f_PU->Get("hJet_Pass")->Clone("j2_data_co");

    double data_co_Lx_1[trackbin] = {10,25,35,45,55,65,74,82,90};
    double data_co_Ly_1[trackbin];
    double data_co_Lxel_1[trackbin] = {10,5,5,5,5,5,4,4,4};
    double data_co_Lxeh_1[trackbin] = {10,5,5,5,5,5,4,4,10};
    double data_co_Lye_1[trackbin];

    double data_co_Lx_2[trackbin] = {10,25,35,45,55,65,74,82,90};
    double data_co_Ly_2[trackbin];
    double data_co_Lxel_2[trackbin] = {10,5,5,5,5,5,4,4,4};
    double data_co_Lxeh_2[trackbin] = {10,5,5,5,5,5,4,4,10};
    double data_co_Lye_2[trackbin];

    double data_co_Lx_3[trackbin] = {10,25,35,45,55,65,74,82,90};
    double data_co_Ly_3[trackbin];
    double data_co_Lxel_3[trackbin] = {10,5,5,5,5,5,4,4,4};;
    double data_co_Lxeh_3[trackbin] = {10,5,5,5,5,5,4,4,10};
    double data_co_Lye_3[trackbin];

    double data_co_cLx_1[trackbin];
    double data_co_cLy_1[trackbin];
    double data_co_cLxel_1[trackbin];
    double data_co_cLxeh_1[trackbin];
    double data_co_cLye_1[trackbin];

    double data_co_cLx_2[trackbin];
    double data_co_cLy_2[trackbin];
    double data_co_cLxel_2[trackbin];
    double data_co_cLxeh_2[trackbin];
    double data_co_cLye_2[trackbin];

    double data_co_cLx_3[trackbin];
    double data_co_cLy_3[trackbin];
    double data_co_cLxel_3[trackbin];
    double data_co_cLxeh_3[trackbin];
    double data_co_cLye_3[trackbin];

    for(int wtrk = 0; wtrk < trackbin; wtrk++){
        TH2D* h1;
        h1 = (TH2D*)f_PU->Get(Form("hSigS_Cor_%d_to_%d_and_%d_to_%d_w_PU_2",trackbinbounds[wtrk], trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();
        h1->Scale(1/(hJ1_f2_PU->GetBinContent(wtrk+1)));
        h1->Scale(1./(BW1));

        TH2D* h2;
        h2 = (TH2D*)f_PU->Get(Form("hBckS_Cor_%d_to_%d_and_%d_to_%d_w_PU_2",trackbinbounds[wtrk], trackbinboundsUpper[wtrk], pt_lo, pt_hi))->Clone();

        data_co_Ly_1[wtrk]  = VOut(h1,h2,28,41, 1);       
        data_co_Lye_1[wtrk] = VEOut(h1,h2,28,41, 1);       

        data_co_Ly_2[wtrk]  = VOut(h1,h2,28,41, 2);       
        data_co_Lye_2[wtrk] = VEOut(h1,h2,28,41, 2);       

        data_co_Ly_3[wtrk]  = VOut(h1,h2,28,41, 3);       
        data_co_Lye_3[wtrk] = VEOut(h1,h2,28,41, 3);       

    }


    for (int i = 0; i < trackbin; i++) {
        data_co_cLx_1[i] = data_co_Lx_1[i];
        data_co_cLy_1[i] = data_co_Ly_1[i];
        data_co_cLxel_1[i] = data_co_Lxel_1[i];
        data_co_cLxeh_1[i] = data_co_Lxeh_1[i];
        data_co_cLye_1[i] = (sqrt (2) )*(data_co_Lye_1[i]);
    }
    for (int i = 0; i < trackbin; i++) {
        data_co_cLx_2[i] = data_co_Lx_2[i];
        data_co_cLy_2[i] = data_co_Ly_2[i];
        data_co_cLxel_2[i] = data_co_Lxel_2[i];
        data_co_cLxeh_2[i] = data_co_Lxeh_2[i];
        data_co_cLye_2[i] = (sqrt (2) )*(data_co_Lye_2[i]);
    }
    for (int i = 0; i < trackbin; i++) {
        data_co_cLx_3[i] = data_co_Lx_3[i];
        data_co_cLy_3[i] = data_co_Ly_3[i];
        data_co_cLxel_3[i] = data_co_Lxel_3[i];
        data_co_cLxeh_3[i] = data_co_Lxeh_3[i];
        data_co_cLye_3[i] = (sqrt (2) )*(data_co_Lye_3[i]);
    }


    auto reco_co_grL1 = new TGraphAsymmErrors(trackbin,reco_co_cLx_1,reco_co_cLy_1,reco_co_cLxel_1,reco_co_cLxeh_1,reco_co_cLye_1,reco_co_cLye_1);
    auto reco_co_grL2 = new TGraphAsymmErrors(trackbin,reco_co_cLx_2,reco_co_cLy_2,reco_co_cLxel_2,reco_co_cLxeh_2,reco_co_cLye_2,reco_co_cLye_2);
    auto reco_co_grL3 = new TGraphAsymmErrors(trackbin,reco_co_cLx_3,reco_co_cLy_3,reco_co_cLxel_3,reco_co_cLxeh_3,reco_co_cLye_3,reco_co_cLye_3);

    auto reco_co_mg1  = new TMultiGraph();


    auto data_co_grL1 = new TGraphAsymmErrors(trackbin,data_co_cLx_1,data_co_cLy_1,data_co_cLxel_1,data_co_cLxeh_1,data_co_cLye_1,data_co_cLye_1);
    auto data_co_grL2 = new TGraphAsymmErrors(trackbin,data_co_cLx_2,data_co_cLy_2,data_co_cLxel_2,data_co_cLxeh_2,data_co_cLye_2,data_co_cLye_2);
    auto data_co_grL3 = new TGraphAsymmErrors(trackbin,data_co_cLx_3,data_co_cLy_3,data_co_cLxel_3,data_co_cLxeh_3,data_co_cLye_3,data_co_cLye_3);










    TH1D* hDiff_v1 = new TH1D("hDiff_v1","hDiff_v1", 9,0,9);
    TH1D* hDiff_v2 = new TH1D("hDiff_v2","hDiff_v2", 9,0,9);
    TH1D* hDiff_v3 = new TH1D("hDiff_v3","hDiff_v3", 9,0,9);

    TH1D* hErr_v1 = new TH1D("hErr_v1","hErr_v1", 9,0,9);
    TH1D* hErr_v2 = new TH1D("hErr_v2","hErr_v2", 9,0,9);
    TH1D* hErr_v3 = new TH1D("hErr_v3","hErr_v3", 9,0,9);

    for(int i=0;i<trackbin;i++){
        hDiff_v1->SetBinContent(i+1, fabs(data_co_cLy_1[i] - reco_co_cLy_1[i]));
        hDiff_v2->SetBinContent(i+1, fabs(data_co_cLy_2[i] - reco_co_cLy_2[i]));
        hDiff_v3->SetBinContent(i+1, fabs(data_co_cLy_3[i] - reco_co_cLy_3[i]));

        hErr_v1->SetBinContent(i+1, TMath::Sqrt(pow(data_co_cLye_1[i],2) + pow(reco_co_cLye_1[i],2)));
        hErr_v2->SetBinContent(i+1, TMath::Sqrt(pow(data_co_cLye_2[i],2) + pow(reco_co_cLye_2[i],2)));
        hErr_v3->SetBinContent(i+1, TMath::Sqrt(pow(data_co_cLye_3[i],2) + pow(reco_co_cLye_3[i],2)));
    }
    TFile* fS_tempA = new TFile("subtraction_reg_PU_2.root", "recreate");

    hDiff_v1->Write();
    hDiff_v2->Write();
    hDiff_v3->Write();

    hErr_v1 ->Write();
    hErr_v2 ->Write();
    hErr_v3 ->Write();

    fS_tempA->Close();



}


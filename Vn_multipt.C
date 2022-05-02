//blh
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TGraphMultiErrors.h"
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
	std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";//+[4]*TMath::Cos(4*x)))";//+[5]*TMath::Cos(5*x)))";
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
	std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";//+[4]*TMath::Cos(4*x)))";
	TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
	func1.SetParameter(0, histfit1->GetMaximum());
	func1.SetParameter(1, 0.1);
	func1.SetParameter(2, 0.1);
	func1.SetParameter(3, 0.1);
	//func1.SetParameter(4, 0.1);
	histfit1->Fit(&func1, "q 0");
	histfit1->Fit(&func1, "q 0");
	histfit1->Fit(&func1, "m q 0");
	histfit1->Fit(&func1, "m q 0");
	histfit1->Fit(&func1, "m q E 0");
	histfit1->Fit(&func1, "m E q 0");
	return func1.GetParError(vn);
}

void Fitter(TH2D* h1, TH2D* h2, int YPlo, int YPhi,int color, int boolsame)
{

	TH1::SetDefaultSumw2(kTRUE);
	TH2::SetDefaultSumw2(kTRUE);

	TH1D *histfit1 = (TH1D*) h1->ProjectionY("",YPlo,YPhi)->Clone();
	TH1D *histfit2 = (TH1D*) h2->ProjectionY("",YPlo,YPhi)->Clone();
	histfit1->Divide(histfit2);
	histfit1->Scale(h2->GetMaximum());
	histfit1->SetMarkerStyle(20);
	std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))";//+[4]*TMath::Cos(4*x)))";
	TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
	func1.SetParameter(0, histfit1->GetMaximum());
	func1.SetParameter(1, 0.1);
	func1.SetParameter(2, 0.1);
	func1.SetParameter(3, 0.1);
	//func1.SetParameter(4, 0.1);
	histfit1->Fit(&func1, "q 0");
	histfit1->Fit(&func1, "q 0");
	histfit1->Fit(&func1, "m q 0");
	histfit1->Fit(&func1, "m q 0");
	histfit1->Fit(&func1, "m q E 0");
	histfit1->Fit(&func1, "m E q 0");
	histfit1->SetStats(kFALSE);
	histfit1->SetMarkerColor(color);
	histfit1->SetMarkerStyle(33);
	histfit1->SetMarkerSize(1);
	if(boolsame == 1){
		histfit1->DrawCopy("PE");
	}
	if(boolsame == 2){
		histfit1->DrawCopy("PE SAME");
	}
	func1.SetLineColor(color);
	func1.DrawCopy("SAME");
}



void Vn_multipt()
{

	const int trackbin 	= 9;
	const int ptbin 	= 3;
	const int vnbin 	= 3;
	const int pubin 	= 3;
	//Main Files
	TFile *f1 = new TFile("bigpapa.root");
	TFile *f2 = new TFile("may1_mc.root");

	const int   MC_trackbinbounds[trackbin]         = { 0,20,30,40,50,58,66,72,79};
	const int   MC_trackbinboundsUpper[trackbin]    = {20,30,40,50,58,66,72,79,1000};

	const int   trackbinbounds[trackbin]         = { 0,20,30,40,50,60,68,74,80};
	const int   trackbinboundsUpper[trackbin]    = {20,30,40,50,60,68,74,80,1000};
	const int ptbinbounds_lo[ptbin]       = {0 ,3 ,5};
	const int  ptbinbounds_hi[ptbin]       = {30,30,30};
	float BW1 = 0.3;
	int ptwant 	= 1;
	int Xscale 	= 800;
	int Yscale 	= 800;
	TCanvas* c1 = new TCanvas("c1","c1", Xscale, Yscale);
	TCanvas* c2 = new TCanvas("c2","c2", Xscale, Yscale);

	// MC MC MC
	TH1D* hJ1_MC;
	hJ1_MC = (TH1D*)f2->Get("hJet_PassW")->Clone("jMC");

	double MC_Lx_1[ptbin][trackbin];// = {10,25,35,45,55,65,74,82};
	double MC_Ly_1[ptbin][trackbin];
	double MC_Lxel_1[ptbin][trackbin];
	double MC_Lxeh_1[ptbin][trackbin];
	double MC_Lye_1[ptbin][trackbin];

	double MC_Lx_2[ptbin][trackbin];// = {10,25,35,45,55,65,74,82};
	double MC_Ly_2[ptbin][trackbin];
	double MC_Lxel_2[ptbin][trackbin];
	double MC_Lxeh_2[ptbin][trackbin];
	double MC_Lye_2[ptbin][trackbin];

	double MC_Lx_3[ptbin][trackbin] = {10,25,35,45,55,65,74,82};
	double MC_Ly_3[ptbin][trackbin];
	double MC_Lxel_3[ptbin][trackbin];
	double MC_Lxeh_3[ptbin][trackbin];
	double MC_Lye_3[ptbin][trackbin];

	double MC_cLx_1[ptbin][trackbin];
	double MC_cLy_1[ptbin][trackbin];
	double MC_cLxel_1[ptbin][trackbin];
	double MC_cLxeh_1[ptbin][trackbin];
	double MC_cLye_1[ptbin][trackbin];

	double MC_cLx_2[ptbin][trackbin];
	double MC_cLy_2[ptbin][trackbin];
	double MC_cLxel_2[ptbin][trackbin];
	double MC_cLxeh_2[ptbin][trackbin];
	double MC_cLye_2[ptbin][trackbin];

	double MC_cLx_3[ptbin][trackbin];
	double MC_cLy_3[ptbin][trackbin];
	double MC_cLxel_3[ptbin][trackbin];
	double MC_cLxeh_3[ptbin][trackbin];
	double MC_cLye_3[ptbin][trackbin];

	for(int wtrk = 0; wtrk < trackbin; wtrk++){
		for(int wppt = 0; wppt < ptbin; wppt++){
			int special_trk = wtrk;

			TH1D* h0;
			h0 = (TH1D*)f2->Get(Form("hBinDist_gen_W_%d",special_trk+1))->Clone();
			//h0 = (TH2D*)f1->Get(Form("hBinDist_cor_%d",wtrk+1))->Clone();
			double new_x = h0->GetMean();
			MC_Lx_1[wppt][wtrk]= new_x;
			MC_Lx_2[wppt][wtrk]= new_x;
			MC_Lx_3[wppt][wtrk]= new_x;

			TH2D* h1;
			h1 = (TH2D*)f2->Get(Form("hSigS_%d_to_%d_and_%d_to_%d_w_PU_1",MC_trackbinbounds[wtrk], MC_trackbinboundsUpper[wtrk], ptbinbounds_lo[wppt], ptbinbounds_hi[wppt]))->Clone();

			h1->Scale(1/(hJ1_MC->GetBinContent(special_trk+1)));
			//h1->Scale(1/(hJ1_f1->GetBinContent(wtrk+1)));
			h1->Scale(1./(BW1));

			TH2D* h2;
			h2 = (TH2D*)f2->Get(Form("hBckS_%d_to_%d_and_%d_to_%d_w_PU_1",MC_trackbinbounds[wtrk], MC_trackbinboundsUpper[wtrk], ptbinbounds_lo[wppt], ptbinbounds_hi[wppt]))->Clone();

			MC_Ly_1[wppt][wtrk]  = VOut(h1,h2,28,41, 1);       
			MC_Lye_1[wppt][wtrk] = VEOut(h1,h2,28,41, 1);       

			MC_Ly_2[wppt][wtrk]  = VOut(h1,h2,28,41, 2);       
			MC_Lye_2[wppt][wtrk] = VEOut(h1,h2,28,41, 2);       

			MC_Ly_3[wppt][wtrk]  = VOut(h1,h2,28,41, 3);       
			MC_Lye_3[wppt][wtrk] = VEOut(h1,h2,28,41, 3);       
		}
	}
	// DATA DATA DATA
	TH1D* hJ1_f1;
	hJ1_f1 = (TH1D*)f1->Get("hJet_Pass")->Clone("j1");

	double Lx_1[ptbin][trackbin] = {10,25,35,45,55,65,74,82};
	double Ly_1[ptbin][trackbin];
	double Lxel_1[ptbin][trackbin];
	double Lxeh_1[ptbin][trackbin];
	double Lye_1[ptbin][trackbin];

	double Lx_2[ptbin][trackbin] = {10,25,35,45,55,65,74,82};
	double Ly_2[ptbin][trackbin];
	double Lxel_2[ptbin][trackbin];
	double Lxeh_2[ptbin][trackbin];
	double Lye_2[ptbin][trackbin];

	double Lx_3[ptbin][trackbin] = {10,25,35,45,55,65,74,82};
	double Ly_3[ptbin][trackbin];
	double Lxel_3[ptbin][trackbin];
	double Lxeh_3[ptbin][trackbin];
	double Lye_3[ptbin][trackbin];

	double cLx_1[ptbin][trackbin];
	double cLy_1[ptbin][trackbin];
	double cLxel_1[ptbin][trackbin];
	double cLxeh_1[ptbin][trackbin];
	double cLye_1[ptbin][trackbin];

	double cLx_2[ptbin][trackbin];
	double cLy_2[ptbin][trackbin];
	double cLxel_2[ptbin][trackbin];
	double cLxeh_2[ptbin][trackbin];
	double cLye_2[ptbin][trackbin];

	double cLx_3[ptbin][trackbin];
	double cLy_3[ptbin][trackbin];
	double cLxel_3[ptbin][trackbin];
	double cLxeh_3[ptbin][trackbin];
	double cLye_3[ptbin][trackbin];

	for(int wtrk = 0; wtrk < trackbin; wtrk++){
		for(int wppt = 0; wppt < ptbin; wppt++){
			int special_trk = wtrk;
			//if(wtrk <5 ) special_trk = wtrk;
			//if(wtrk >=5 ) special_trk = wtrk+4;

			TH1D* h0;
			h0 = (TH1D*)f1->Get(Form("hBinDist_cor_%d",special_trk+1))->Clone();
			//h0 = (TH2D*)f1->Get(Form("hBinDist_cor_%d",wtrk+1))->Clone();
			double new_x = h0->GetMean();
			Lx_1[wppt][wtrk]= new_x;
			Lx_2[wppt][wtrk]= new_x;
			Lx_3[wppt][wtrk]= new_x;

			TH2D* h1;
			h1 = (TH2D*)f1->Get(Form("hSigS_Cor_%d_to_%d_and_%d_to_%d_w_PU_1",trackbinbounds[wtrk], trackbinboundsUpper[wtrk], ptbinbounds_lo[wppt], ptbinbounds_hi[wppt]))->Clone();

			h1->Scale(1/(hJ1_f1->GetBinContent(special_trk+1)));
			//h1->Scale(1/(hJ1_f1->GetBinContent(wtrk+1)));
			h1->Scale(1./(BW1));

			TH2D* h2;
			h2 = (TH2D*)f1->Get(Form("hBckS_Cor_%d_to_%d_and_%d_to_%d_w_PU_1",trackbinbounds[wtrk], trackbinboundsUpper[wtrk], ptbinbounds_lo[wppt], ptbinbounds_hi[wppt]))->Clone();

			Ly_1[wppt][wtrk]  = VOut(h1,h2,28,41, 1);       
			Lye_1[wppt][wtrk] = VEOut(h1,h2,28,41, 1);       

			Ly_2[wppt][wtrk]  = VOut(h1,h2,28,41, 2);       
			Lye_2[wppt][wtrk] = VEOut(h1,h2,28,41, 2);       

			Ly_3[wppt][wtrk]  = VOut(h1,h2,28,41, 3);       
			Lye_3[wppt][wtrk] = VEOut(h1,h2,28,41, 3);       
		}
	}

	c1->cd();

	for (int i = 0; i < trackbin; i++) {
		for (int j = 0; j < ptbin; j++) {
			cLx_1[j][i] = Lx_1[j][i];
			cLy_1[j][i] = Ly_1[j][i];
			cLxel_1[j][i] = Lxel_1[j][i];
			cLxeh_1[j][i] = Lxeh_1[j][i];
			cLye_1[j][i] = (sqrt(2))*(Lye_1[j][i]);
		}
	}
	for (int i = 0; i < trackbin; i++) {
		for (int j = 0; j < ptbin; j++) {
			cLx_2[j][i] = Lx_2[j][i];
			cLy_2[j][i] = Ly_2[j][i];
			cLxel_2[j][i] = Lxel_2[j][i];
			cLxeh_2[j][i] = Lxeh_2[j][i];
			cLye_2[j][i] = (sqrt(2))*(Lye_2[j][i]);
		}
	}
	for (int i = 0; i < trackbin; i++) {
		for (int j = 0; j < ptbin; j++) {
			cLx_3[j][i] = Lx_3[j][i];
			cLy_3[j][i] = Ly_3[j][i];
			cLxel_3[j][i] = Lxel_3[j][i];
			cLxeh_3[j][i] = Lxeh_3[j][i];
			cLye_3[j][i] = (sqrt(2))*(Lye_3[j][i]);
		}
	}

	//SYSYSYSYSYSYSYSYSYSYSYSYSYSYYSYSYSYSYYSYSYSYSYYS
	//SYSYSYSYSYSYSYSYSYSYSYSYSYSYYSYSYSYSYYSYSYSYSYYS
	//SYSYSYSYSYSYSYSYSYSYSYSYSYSYYSYSYSYSYYSYSYSYSYYS
	//SYSYSYSYSYSYSYSYSYSYSYSYSYSYYSYSYSYSYYSYSYSYSYYS

	double cLye_1_sys[ptbin][trackbin];
	double cLye_2_sys[ptbin][trackbin];
	double cLye_3_sys[ptbin][trackbin];

	TFile *f_jPt[ptbin];
	//FILES FOR Jet PT SYStematics
	f_jPt[0] = new TFile("sys_ptBin_jetpt/sys_jetpt_0_30_83.root");
	f_jPt[1] = new TFile("sys_ptBin_jetpt/sys_jetpt_3_30_83.root");
	f_jPt[2] = new TFile("sys_ptBin_jetpt/sys_jetpt_5_30_83.root");

	TFile *f_jEP[ptbin];
	//FILES FOR Jet EP SYStematics
	f_jEP[0] = new TFile("sys_ptBin_etaphi/sys_etaphi_0_30_83.root");
	f_jEP[1] = new TFile("sys_ptBin_etaphi/sys_etaphi_3_30_83.root");
	f_jEP[2] = new TFile("sys_ptBin_etaphi/sys_etaphi_5_30_83.root");

	TFile *f_PU[ptbin][pubin];
	//FILES FOR Jet EP SYStematics
	f_PU[0][0] = new TFile("sys_ptBin_PU/sys_jetpt_0_30_83_PU1.root");
	f_PU[0][1] = new TFile("sys_ptBin_PU/sys_jetpt_0_30_83_PU2.root");
	f_PU[0][2] = new TFile("sys_ptBin_PU/sys_jetpt_0_30_83_PU3.root");
	f_PU[1][0] = new TFile("sys_ptBin_PU/sys_jetpt_3_30_83_PU1.root");
	f_PU[1][1] = new TFile("sys_ptBin_PU/sys_jetpt_3_30_83_PU2.root");
	f_PU[1][2] = new TFile("sys_ptBin_PU/sys_jetpt_3_30_83_PU3.root");
	f_PU[2][0] = new TFile("sys_ptBin_PU/sys_jetpt_5_30_83_PU1.root");
	f_PU[2][1] = new TFile("sys_ptBin_PU/sys_jetpt_5_30_83_PU2.root");
	f_PU[2][2] = new TFile("sys_ptBin_PU/sys_jetpt_5_30_83_PU3.root");


	TH1::SetDefaultSumw2(kTRUE);
	TH2::SetDefaultSumw2(kTRUE);

	TH1D* h_jPt[ptbin][vnbin];
	TH1D* h_jEP[ptbin][vnbin];
	TH1D* h_PU[ptbin][vnbin][pubin];

	for(int pt_i = 0; pt_i < ptbin; pt_i++){
		for(int vn_i = 0; vn_i < vnbin; vn_i++){
			h_jPt[pt_i][vn_i] =(TH1D*)f_jPt[pt_i]->Get(Form("hDiff_v%d",vn_i+1))->Clone(Form("pt_%d_%d_%d",pt_i,vn_i,vn_i+pt_i));
			h_jEP[pt_i][vn_i] =(TH1D*)f_jEP[pt_i]->Get(Form("hDiff_v%d",vn_i+1))->Clone(Form("EP_%d_%d_%d",pt_i,vn_i,vn_i+pt_i));
			for(int pu_i=0; pu_i<pubin;pu_i++){
				h_PU[pt_i][vn_i][pu_i]  = (TH1D*)f_PU[pt_i][pu_i]->Get(Form("hDiff_v%d",vn_i+1))->Clone(Form("PU_%d_%d_%d",pt_i,vn_i,vn_i+pt_i));
			}
		}
	}

	for (int pt_i = 0; pt_i < ptbin; pt_i++) {
		for (int i = 0; i < trackbin; i++) {
			cLye_1_sys[pt_i][i]=1.0 * sqrt(
					pow(h_jPt[pt_i][0]->GetBinContent(i+1),2) 
					+ pow(
						pow(
							(1/3)*(  
								pow( (h_PU[pt_i][0][0]->GetBinContent(i+1)),2)
								+pow( (h_PU[pt_i][0][1]->GetBinContent(i+1)),2)
								+pow( (h_PU[pt_i][0][2]->GetBinContent(i+1)),2)
							      )
							,0.5)
						,2)	
					+ pow(h_jEP[pt_i][0]->GetBinContent(i+1),2) 
					);

			cLye_2_sys[pt_i][i]=1.0 * sqrt( 
					pow(h_jPt[pt_i][1]->GetBinContent(i+1),2) 
					+ pow(
						pow(
							(1/3)*(  
								pow( (h_PU[pt_i][1][0]->GetBinContent(i+1)),2)
								+pow( (h_PU[pt_i][1][1]->GetBinContent(i+1)),2)
								+pow( (h_PU[pt_i][1][2]->GetBinContent(i+1)),2)
							      )
							,0.5)
						,2)	
					+ pow(h_jEP[pt_i][1]->GetBinContent(i+1),2) 
					);

			cLye_3_sys[pt_i][i]=1.0 * sqrt( 
					pow(h_jPt[pt_i][2]->GetBinContent(i+1),2) 
					+ pow(
						pow(
							(1/3)*(  
								pow( (h_PU[pt_i][2][0]->GetBinContent(i+1)),2)
								+pow( (h_PU[pt_i][2][1]->GetBinContent(i+1)),2)
								+pow( (h_PU[pt_i][2][2]->GetBinContent(i+1)),2)
							      )
							,0.5)
						,2)	
					+ pow(h_jEP[pt_i][2]->GetBinContent(i+1),2) 
					);
		}
	}

	TGraphMultiErrors* MC_grL1[ptbin];
	TGraphMultiErrors* MC_grL2[ptbin];
	TGraphMultiErrors* MC_grL3[ptbin];

	TGraphMultiErrors* grL1[ptbin];
	TGraphMultiErrors* grL2[ptbin];
	TGraphMultiErrors* grL3[ptbin];

	int v1color[ptbin] = {2,3,4};
	int v2color[ptbin] = {2,3,4};
	int v3color[ptbin] = {2,3,4};
	int v1markerC[ptbin] = {2,3,4};
	int v2markerC[ptbin] = {2,3,4};
	int v3markerC[ptbin] = {2,3,4};

	int v1style[ptbin] = {1,1,1};
	int v2style[ptbin] = {1,1,1};
	int v3style[ptbin] = {1,1,1};

	int v1marker[ptbin] = {20,21,22};
	int v2marker[ptbin] = {20,21,22};
	int v3marker[ptbin] = {20,21,22};

	int MC_v1style[ptbin] = {3,3,3};
	int MC_v2style[ptbin] = {3,3,3};
	int MC_v3style[ptbin] = {3,3,3};

	int MC_v1marker[ptbin] = {24,25,26};
	int MC_v2marker[ptbin] = {24,25,26};
	int MC_v3marker[ptbin] = {24,25,26};

	auto legendTrack = new TLegend(0.53, 0.785, 0.80, 0.89);
	legendTrack->SetLineColor(0);
	legendTrack->SetFillStyle(0);
	legendTrack->SetLineWidth(0);

	double zero[trackbin] = {0,0,0,0,0,0,0,0,0};

	auto mg1 = new TMultiGraph();

	bool bool_pt = 1;
	int low_range = 0;
	int high_range = 0;
	if(bool_pt){ 
		low_range = 0;	
		high_range = ptbin;
	}else{ 
		low_range = 1;
		high_range = 2;
	}


	for( int ii = low_range; ii < high_range; ii++){
		MC_grL1[ii] = new TGraphMultiErrors("","",trackbin,MC_Lx_1[ii],MC_Ly_1[ii],zero,zero,MC_Lye_1[ii],MC_Lye_1[ii]);
		MC_grL2[ii] = new TGraphMultiErrors("","",trackbin,MC_Lx_2[ii],MC_Ly_2[ii],zero,zero,MC_Lye_2[ii],MC_Lye_2[ii]);
		MC_grL3[ii] = new TGraphMultiErrors("","",trackbin,MC_Lx_3[ii],MC_Ly_3[ii],zero,zero,MC_Lye_3[ii],MC_Lye_3[ii]);

		MC_grL1[ii]->SetLineColor(v1color[ii]);
		MC_grL1[ii]->SetMarkerStyle(MC_v1marker[ii]);
		MC_grL1[ii]->SetMarkerColor(v1markerC[ii]);
		MC_grL1[ii]->SetLineStyle(MC_v1style[ii]);

		MC_grL2[ii]->SetLineColor(v2color[ii]);
		MC_grL2[ii]->SetMarkerStyle(MC_v2marker[ii]);
		MC_grL2[ii]->SetMarkerColor(v2markerC[ii]);
		MC_grL2[ii]->SetLineStyle(MC_v2style[ii]);

		MC_grL3[ii]->SetLineColor(v3color[ii]);
		MC_grL3[ii]->SetMarkerStyle(MC_v3marker[ii]);
		MC_grL3[ii]->SetMarkerColor(v3markerC[ii]);
		MC_grL3[ii]->SetLineStyle(MC_v3style[ii]);

		//mg1->Add(grL1[ii],"pc");
		mg1->Add(MC_grL2[ii],"pc");
		//mg1->Add(grL3[ii],"pc");
	}

	for( int ii = low_range; ii < high_range; ii++){
		//grL1[ii] = new TGraphMultiErrors("","",trackbin,Lx_1[ii],Ly_1[ii],zero,zero,Lye_1[ii],Lye_1[ii]);
		grL2[ii] = new TGraphMultiErrors("","",trackbin,Lx_2[ii],Ly_2[ii],zero,zero,Lye_2[ii],Lye_2[ii]);
		//grL3[ii] = new TGraphMultiErrors("","",trackbin,Lx_3[ii],Ly_3[ii],zero,zero,Lye_3[ii],Lye_3[ii]);

		//grL1[ii]->SetLineColor(v1color[ii]);
		//grL1[ii]->SetMarkerStyle(v1marker[ii]);
		//grL1[ii]->SetMarkerColor(v1markerC[ii]);
		//grL1[ii]->SetLineStyle(v1style[ii]);

		grL2[ii]->SetLineColor(v2color[ii]);
		grL2[ii]->SetMarkerStyle(v2marker[ii]);
		grL2[ii]->SetMarkerColor(v2markerC[ii]);
		grL2[ii]->SetLineStyle(v2style[ii]);

		//grL3[ii]->SetLineColor(v3color[ii]);
		//grL3[ii]->SetMarkerStyle(v3marker[ii]);
		//grL3[ii]->SetMarkerColor(v3markerC[ii]);
		//grL3[ii]->SetLineStyle(v3style[ii]);

		//legendTrack->AddEntry(grL1[ii],Form("V_{1}{2}, %.1f < j_{T} < 3.0 GeV",ptbinbounds_lo[ii]/10.0),"pl");
		//legendTrack->AddEntry(grL2[ii],Form("V_{2}{2}, %.1f < j_{T} < 3.0 GeV",ptbinbounds_lo[ii]/10.0),"pl");
		//legendTrack->AddEntry(grL3[ii],Form("V_{3}{2}, %.1f < j_{T} < 3.0 GeV",ptbinbounds_lo[ii]/10.0),"pl");

		//mg1->Add(grL1[ii],"pc");
		mg1->Add(grL2[ii],"pc");
		//mg1->Add(grL3[ii],"pc");
	}




	double sys_x[trackbin];
	for(int i =0; i<trackbin; i++){
		sys_x[i]=0.5;
	}

	TGraphMultiErrors* grL1_0_sys = new TGraphMultiErrors("","",trackbin,Lx_1[0],Ly_1[0],sys_x,sys_x,cLye_1_sys[0],cLye_1_sys[0]);
	grL1_0_sys->SetFillColor(2);
	grL1_0_sys->SetFillStyle(3001);
	//if(bool_pt) mg1->Add(grL1_0_sys,"2");

	TGraphMultiErrors* grL2_0_sys = new TGraphMultiErrors("","",trackbin,Lx_2[0],Ly_2[0],sys_x,sys_x,cLye_2_sys[0],cLye_2_sys[0]);
	grL2_0_sys->SetFillColor(2);
	grL2_0_sys->SetFillStyle(3001);
	if(bool_pt) mg1->Add(grL2_0_sys,"2");

	TGraphMultiErrors* grL3_0_sys = new TGraphMultiErrors("","",trackbin,Lx_3[0],Ly_3[0],sys_x,sys_x,cLye_3_sys[0],cLye_3_sys[0]);
	grL3_0_sys->SetFillColor(2);
	grL3_0_sys->SetFillStyle(3001);
	//if(bool_pt) mg1->Add(grL3_0_sys,"2");

	TGraphMultiErrors* grL1_1_sys = new TGraphMultiErrors("","",trackbin,Lx_1[1],Ly_1[1],sys_x,sys_x,cLye_1_sys[1],cLye_1_sys[1]);
	grL1_1_sys->SetFillColor(2);
	grL1_1_sys->SetFillStyle(3001);
	//mg1->Add(grL1_1_sys,"2");

	TGraphMultiErrors* grL2_1_sys = new TGraphMultiErrors("","",trackbin,Lx_2[1],Ly_2[1],sys_x,sys_x,cLye_2_sys[1],cLye_2_sys[1]);
	grL2_1_sys->SetFillColor(2);
	grL2_1_sys->SetFillStyle(3001);
	mg1->Add(grL2_1_sys,"2");

	TGraphMultiErrors* grL3_1_sys = new TGraphMultiErrors("","",trackbin,Lx_3[1],Ly_3[1],sys_x,sys_x,cLye_3_sys[1],cLye_3_sys[1]);
	grL3_1_sys->SetFillColor(2);
	grL3_1_sys->SetFillStyle(3001);
	//mg1->Add(grL3_1_sys,"2");

	TGraphMultiErrors* grL1_2_sys = new TGraphMultiErrors("","",trackbin,Lx_1[2],Ly_1[2],sys_x,sys_x,cLye_1_sys[2],cLye_1_sys[2]);
	grL1_2_sys->SetFillColor(2);
	grL1_2_sys->SetFillStyle(3001);
	//if(bool_pt) mg1->Add(grL1_2_sys,"2");

	TGraphMultiErrors* grL2_2_sys = new TGraphMultiErrors("","",trackbin,Lx_2[2],Ly_2[2],sys_x,sys_x,cLye_2_sys[2],cLye_2_sys[2]);
	grL2_2_sys->SetFillColor(2);
	grL2_2_sys->SetFillStyle(3001);
	if(bool_pt) mg1->Add(grL2_2_sys,"2");

	TGraphMultiErrors* grL3_2_sys = new TGraphMultiErrors("","",trackbin,Lx_3[2],Ly_3[2],sys_x,sys_x,cLye_3_sys[2],cLye_3_sys[2]);
	grL3_2_sys->SetFillColor(2);
	grL3_2_sys->SetFillStyle(3001);
	//if(bool_pt) mg1->Add(grL3_2_sys,"2");




	double l05y[2] = {(.05*.05),(.05*.05)};
	double l05x[2] = {0,110};

	auto l05 = new TGraph(2, l05x,l05y);

	l05->SetLineColor(7);
	l05->SetLineWidth(3);
	l05->SetLineStyle(2);

	double l10y[2] = {(.1*.1),(.1*.1)};
	double l10x[2] = {0,110};

	auto l10 = new TGraph(2, l10x,l10y);

	l10->SetLineColor(7);
	l10->SetLineWidth(3);
	l10->SetLineStyle(7);

	double l15y[2] = {(.15*.15),(.15*.15)};
	double l15x[2] = {0,110};

	auto l15 = new TGraph(2, l15x,l15y);

	l15->SetLineColor(7);
	l15->SetLineWidth(3);
	l15->SetLineStyle(9);

	mg1->Add(l05,"L");
	mg1->Add(l10,"L");
	mg1->Add(l15,"L");

	c1->cd();

	c1->SetTickx(1);
	c1->SetTicky(1);
	c1->SetLeftMargin(0.15);
	c1->SetTopMargin(0.05);
	mg1->Draw("AP");

	mg1->GetXaxis()->SetLimits(0,110);
	mg1->GetXaxis()->SetTitle("N_{ch}^{j}");
	mg1->GetYaxis()->SetTitle("V_{n#Delta}");
	mg1->GetYaxis()->SetTitleOffset(1.55);
	mg1->GetXaxis()->SetTitleOffset(1);
	mg1->GetYaxis()->CenterTitle(true);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetXaxis()->SetTitleSize(.04);
	mg1->GetYaxis()->SetTitleSize(.04);
	//mg1->GetYaxis()->SetTitleOffset(.04);
	gPad->Modified();
	//auto legend21 = new TLegend(0.16, 0.39, 0.46, 0.5);

	//legend21->Draw();

	auto legend22 = new TLegend(0.5, 0.11, 0.85, 0.38);
	legend22->AddEntry("","CMS pp 13 TeV, 147 fb^{-1} ","");
	legend22->AddEntry("","Anti k_{t} R=0.8","");
	legend22->AddEntry("","p_{T}^{jet} > 500 GeV","");
	legend22->AddEntry("","|#eta^{jet}| < 1.6","");
	legend22->SetLineColor(0);
	legend22->SetFillStyle(0);
	legend22->SetLineWidth(0);
	legend22->Draw();

	auto legend23 = new TLegend(0.53, 0.785, 0.80, 0.89);
	legend23->AddEntry(l15,"v_{2}{2} = 15%","l");
	legend23->AddEntry(l10,"v_{2}{2} = 10%","l");
	legend23->AddEntry(l05,"v_{2}{2} =  5%","l");
	legend23->SetLineColor(0);
	legend23->SetFillStyle(0);
	legend23->SetLineWidth(0);
	legend23->Draw();
	legendTrack->Draw();

	TH1D* hj550;
	hj550 = (TH1D*)f1->Get("hJet_Pass550")->Clone("j550");
	hj550->SetMarkerColor(3);
	hj550->SetMarkerStyle(21);

	TH1D* hj500;
	hj500 = (TH1D*)f1->Get("hJet_Pass500")->Clone("j500");
	hj500->SetMarkerColor(2);
	hj500->SetMarkerStyle(20);

	c2->cd();
	hj550->Draw("p");
	hj500->Draw("p SAME");
	auto legendJET = new TLegend(0.5, 0.11, 0.85, 0.38);
	legendJET->AddEntry(hj550, "550 Jets", "p");
	legendJET->AddEntry(hj500, "500 Jets", "p");
	legendJET->Draw();


}
/*
 *
 *
 double tempX1[trackbin];
 memcpy(tempX1, Lx_1[ii], sizeof(tempX1));	
 double tempX2[trackbin];
 std::copy(std::begin(Lx_1[ii]), std::end(Lx_1[ii]), std::begin(tempX2));
 *
 for(int jjj=0; jjj<trackbin; jjj++){
 cout << tempX1[jjj] << " , ";
 }
 cout << "tempX2:  " <<endl;
 for(int jjj=0; jjj<trackbin; jjj++){
 cout << tempX2[jjj] << " , ";
 }

 cout << "lx_1 all:  " <<endl;
 for(int jjj=0; jjj<ptbin; jjj++){
 for(int qqq=0; qqq<trackbin; qqq++){
 cout << Lx_1[jjj][qqq] << " , ";
 }
 cout << " " << endl;
 }

 */


/*
//goal is to create a temp array that copies existing array, does not preserve after-copy changes to existing array
//only idneitcal at time of copy and both arrays can be independelty changed afterward.

//target is an existing array

double temp1[elem_N]; //way one to initialize array
memcpy(temp1, target, sizeof(temp1)); //way one to copy array

std::array<double, elem_N> temp2;	//way two to initialize array
std::copy(std::begin(target), std::end(target), std::begin(temp2)); //way two to copy array
 */




/*
//DATA DATA DATA
grL1[ii]=(TGraphMultiErrors*) new TGraphMultiErrors("","",trackbin,tempX1,tempX1,tempX1,tempX1,tempX1,tempX1);
//grL1[ii]->SetMarkerColor(2);
//grL1[ii]->SetMarkerStyle(20);
//grL1[ii]->Draw("P");
//grL1[ii]->GetXaxis()->SetRangeUser(0,100);
grL2[ii]= new TGraphMultiErrors("","",trackbin,cLx_2[ii],cLy_2[ii],cLxel_2[ii],cLxeh_2[ii],cLye_2[ii],cLye_2[ii]);
//grL2[ii]->SetTitle("v_{2}{2}");
grL2[ii]->SetMarkerColor(4);
grL2[ii]->SetMarkerStyle(20);
grL3[ii] = new TGraphMultiErrors("","",trackbin,cLx_3[ii],cLy_3[ii],cLxel_3[ii],cLxeh_3[ii],cLye_3[ii],cLye_3[ii]);
//grL3[ii]->SetTitle("v_{3}{2}");
grL3[ii]->SetMarkerColor(3);
grL3[ii]->SetMarkerStyle(20);

grL1[ii]->SetLineWidth(1);
grL2[ii]->SetLineWidth(1);
grL3[ii]->SetLineWidth(1);

grL1[ii]->SetLineColor(2);
grL2[ii]->SetLineColor(4);
grL3[ii]->SetLineColor(3);
mg1->Add(grL1[ii],"pc");
//		mg1->Add(grL2[ii],"pc");
//		mg1->Add(grL3[ii],"pc");
}
 */

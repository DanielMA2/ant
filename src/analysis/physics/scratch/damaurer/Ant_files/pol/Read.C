#define analysis_cxx
#include "TH2.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <math.h>
#include "TLine.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TClonesArray.h"
#include "PParticle.h"

Double_t Pol_Bremsstrahlung_fit(Double_t *x,Double_t *par) {
     Double_t fitval = par[0]*((1/x[0])+par[1]*TMath::Gaus(x[0],0.301,0.0377,kFALSE)-par[2]);
     return fitval;
}

Double_t Unpol_Bremsstrahlung_fit(Double_t *x,Double_t *par) {
     Double_t fitval = par[0]*(1/x[0]);
     return fitval;
}

Double_t peak_fit(Double_t *x,Double_t *par) {
     Double_t fitval = par[0]*TMath::Gaus(x[0],par[1],par[2]);
     return fitval;
}


void Optics(TH1 *hist){
	
	hist->SetLineColor(4);
	hist->SetMarkerStyle(3);
	hist->SetMarkerColor(4);
	hist->SetMarkerSize(1);

}

void Read(){
TFile *data = new TFile("pol_brems_gp_ppi0_100kEvs_ant_output.root");

bool polBeam = true;
 
TH1 *h_TaggerTime = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/Tagger/h_TaggerTime");
TH1 *h_nClusters = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/Tagger/h_nClusters");
TH1 *h_VetoEnergies = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/Tagger/h_VetoEnergies");
TH1 *h_InitialBeamE = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/Tagger/h_InitialBeamE");

TH1 *h_PolarAngles = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/PolarAngles/h_PolarAngles");
TH1 *h_proton_PolarAngles = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/PolarAngles/h_ProtonPolarAngles");
TH1 *h_photon_PolarAngles = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/PolarAngles/h_PhotonPolarAngles");

TH1 *h_AzimuthAngles = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/AzimuthAngles/h_AzimuthAngles");
TH1 *h_proton_AzimuthAngles = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/AzimuthAngles/h_ProtonAzimuthAngles");
TH1 *h_photon_AzimuthAngles = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/AzimuthAngles/h_PhotonAzimuthAngles");

TH1 *h_PhotonIm = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/Invariant mass/h_PhotonIm");
TH1 *h_ProtonIm = (TH1*)data->Get("scratch_damaurer_MC_ppi0_pgg/Invariant mass/h_ProtonIm");

TH2 *h_PhotonETheta = (TH2*)data->Get("scratch_damaurer_MC_ppi0_pgg/Particle energies/h_PhotonETheta");
TH2 *h_ProtonETheta = (TH2*)data->Get("scratch_damaurer_MC_ppi0_pgg/Particle energies/h_ProtonETheta");
TH2 *h_PhotonEPhi = (TH2*)data->Get("scratch_damaurer_MC_ppi0_pgg/Particle energies/h_PhotonEPhi");
TH2 *h_ProtonEPhi = (TH2*)data->Get("scratch_damaurer_MC_ppi0_pgg/Particle energies/h_ProtonEPhi");

//Setting the optics of the histograms

Optics(h_InitialBeamE);
Optics(h_TaggerTime);
Optics(h_AzimuthAngles);
Optics(h_proton_AzimuthAngles);
Optics(h_photon_AzimuthAngles);
Optics(h_PhotonIm);
Optics(h_ProtonIm);
Optics(h_VetoEnergies);
Optics(h_nClusters);
Optics(h_PolarAngles);
Optics(h_photon_PolarAngles);
Optics(h_proton_PolarAngles);

//Plotting & saving basic histograms from tutorial

TCanvas *c1 = new TCanvas("c1","polBeam_basic_plots");
c1->Divide(3);
c1->cd(1);
gStyle->SetOptStat(1);
h_TaggerTime->SetTitle("Tagger time");
h_TaggerTime->Draw();
c1->cd(2);
gStyle->SetOptStat(1);
h_nClusters->SetTitle("n_Clusters");
h_nClusters->Draw();
c1->cd(3);
gStyle->SetOptStat(1);
h_VetoEnergies->SetTitle("Deposited Veto energies of all candidates");
h_VetoEnergies->Draw();
c1->SaveAs("new_hists/basic_hists.pdf");

//Plotting & saving Polar angle histograms

TLine *line = new TLine(20,0,20,4500);
line->SetLineColorAlpha(kGreen, 0.35);

TLine *line2 = new TLine(160,0,160,4500);
line2->SetLineColorAlpha(kGreen, 0.35);

TCanvas *c2 = new TCanvas("c2","polBeam_Polar_angle_hists");
c2->Divide(3);
c2->cd(1);
gStyle->SetOptStat(1);
h_PolarAngles->SetTitle("Particle polar angles");
//h_PolarAngles->GetYaxis()->SetTitle("# Entries");
h_PolarAngles->Draw();
//h_PolarAngles->GetYaxis()->SetRange(0.,6000.);
line->Draw();
line2->Draw();
c2->cd(2);
gStyle->SetOptStat(1);
h_photon_PolarAngles->SetTitle("Photon polar angles");
//h_PolarAngles->GetYaxis()->SetTitle("# Entries");
h_photon_PolarAngles->Draw();
//h_PolarAngles->GetYaxis()->SetRange(0.,6000.);
line->Draw();
line2->Draw();
c2->cd(3);
gStyle->SetOptStat(1);
h_proton_PolarAngles->SetTitle("Proton polar angles");
//h_PolarAngles->GetYaxis()->SetTitle("# Entries");
h_proton_PolarAngles->Draw();
//h_PolarAngles->GetYaxis()->SetRange(0.,6000.);
line->Draw();
line2->Draw();
c2->SaveAs("new_hists/polBeam_Polar_angle_hists.pdf");

//Plotting & saving azimuth angle histograms

TCanvas *c3 = new TCanvas("c3","polBeam_Phi_angle_hists");
c3->Divide(3);
c3->cd(1);
h_AzimuthAngles->SetTitle("Azimuthangles of all particles");
h_AzimuthAngles->Draw();
c3->cd(2);
h_photon_AzimuthAngles->SetTitle("Photon azimuthangles");
h_photon_AzimuthAngles->Draw();
c3->cd(3);
h_proton_AzimuthAngles->SetTitle("Proton azimuthangles");
h_proton_AzimuthAngles->Draw();
c3->SaveAs("new_hists/polBeam_Azimuth_angle_hists.pdf");

//Plotting & saving particle invariant mass histograms

TCanvas *c4 = new TCanvas("c4","polBeam_IM_hists");
c4->Divide(2);
c4->cd(1);
h_PhotonIm->SetTitle("Photons invariant masses");
TF1 *func1 = new TF1("fit1","gaus",100,160);
func1->SetParameter(1,130);
func1->SetLineColor(kGreen);
gStyle->SetOptFit();
h_PhotonIm->Fit(func1,"R");
h_PhotonIm->Draw();
c4->cd(2);
h_ProtonIm->SetTitle("Protons invariant masses");
TF1 *func2 = new TF1("fit2","gaus",900,1000);
func2->SetParameter(1,130);
func2->SetLineColor(kGreen);
gStyle->SetOptFit();
h_ProtonIm->Fit(func2,"R");
h_ProtonIm->Draw();
c4->SaveAs("new_hists/polBeam_Im_hists.pdf");

//Plotting & saving particle energy vs theta histograms

TCanvas *c5 = new TCanvas("c5","polBeam_Energy_Theta_hists");
c5->Divide(2,2);
c5->cd(1);
h_PhotonETheta->SetTitle("Photons - Energy vs Theta");
h_PhotonETheta->Draw("pcolz");
c5->cd(2);
h_ProtonETheta->SetTitle("Protons - Energy vs Theta");
h_ProtonETheta->Draw("pcolz");
c5->cd(3);
h_PhotonEPhi->SetTitle("Photons - Energy vs phi");
h_PhotonEPhi->Draw("pcolz");
c5->cd(4);
h_ProtonEPhi->SetTitle("Protons - Energy vs phi");
h_ProtonEPhi->Draw("pcolz");
c5->SaveAs("new_hists/polBeam_EvsAngles_hists.pdf");

//Plotting & saving inital photonbeam

TCanvas *c6 = new TCanvas("c6","Initial Photon beam");
h_InitialBeamE->SetTitle("Initial Photonbeam - Energy distribution");
TF1 *func3;

if(polBeam){
	func3 = new TF1("fit3",Pol_Bremsstrahlung_fit,0.16,0.75,3);
}
else{
	func3 = new TF1("fit3",Unpol_Bremsstrahlung_fit,0.16,0.75,1);
}

func3->SetLineColor(kGreen);
gStyle->SetOptFit();
h_InitialBeamE->Fit(func3,"R");
h_InitialBeamE->Draw();
c6->SaveAs("new_hists/polBeam_InitialBeamE.pdf");

}

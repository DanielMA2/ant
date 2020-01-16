#define analysis_cxx
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <math.h>
#include "TLine.h"
#include "TBox.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TClonesArray.h"
#include "PParticle.h"
#include "TLegend.h"
#include "TLegendEntry.h"

void PaintBin (TH1 *h, Int_t bin, Int_t color) 
{
   printf("%d %d %f\n", bin, color, h->GetBinContent(bin));
   TBox *b = new TBox(h->GetBinLowEdge(bin),
                      h->GetMinimum(),
                      h->GetBinWidth(bin)+h->GetBinLowEdge(bin),
                      h->GetBinContent(bin));
   b->SetFillColor(color);
   b->Draw();

}

void test(){

        Int_t number_of_bins = 10000;
	Int_t lower_edge = 0;
	Int_t upper_edge = 6;
	int steps = (int)(number_of_bins/(upper_edge-lower_edge));
	TH1* hist=new TH1F("hist","hist",number_of_bins,lower_edge,upper_edge);

        hist->GetXaxis()->SetTitle("Cases");
        hist->GetYaxis()->SetTitle("#");

        hist->SetBinContent(1*steps,10);
        hist->SetBinContent(2*steps,20);
        hist->SetBinContent(3*steps,25);

	TCanvas *c1 = new TCanvas("c1","hist");
        gStyle->SetOptStat(0);
  	c1->cd();
   	hist->Draw();
        
        auto legend = new TLegend(0.6,0.7,0.9,0.9);
        legend->SetHeader("Legend","C"); 

        PaintBin (hist, 1*steps, kRed);
        PaintBin (hist, 2*steps, kBlue);
        PaintBin (hist, 3*steps, kGreen);

        TLegendEntry* l1 = legend->AddEntry((TObject*)0,"1) ...","l");
        l1->SetTextColor(kRed);
        l1->SetLineColor(kRed);

        TLegendEntry* l2 = legend->AddEntry((TObject*)0,"2) ...","l");
        l2->SetTextColor(kBlue);
        l2->SetLineColor(kBlue);

        TLegendEntry* l3 = legend->AddEntry((TObject*)0,"3) ...","l");
        l3->SetTextColor(kGreen);
        l3->SetLineColor(kGreen);

        legend->Draw();
   	c1->SaveAs("hist.pdf");

}

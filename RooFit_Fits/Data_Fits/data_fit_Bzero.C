using namespace RooFit;

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include "ACCSEL.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooProduct.h>
#include <RooFormulaVar.h>
#include <RooExponential.h>
#include <RooGenericPdf.h>  
#include <RooProdPdf.h>
#include <RooBifurGauss.h>
#include <RooBinning.h>
#include <RooChebychev.h>  

#include <RooFitResult.h>
#include <RooFit.h>
#include <RooCmdArg.h>
#include <RooCurve.h>
#include <RooExtendPdf.h>


// Auxiliar Function to Compute Black Dotted Lines to the Height
double getYatMass(RooPlot* frame, double mass) {
    for (int i = 0; i < frame->numItems(); ++i) {
        RooCurve* curve = dynamic_cast<RooCurve*>(frame->getObject(i));
        if (!curve) continue;
        int n = curve->GetN();
        double* x = curve->GetX();
        double* y = curve->GetY();
        for (int j = 0; j < n - 1; ++j) {
            if (x[j] <= mass && mass <= x[j+1]) {
                double slope = (y[j+1] - y[j]) / (x[j+1] - x[j]);
                return y[j] + slope * (mass - x[j]);
            }
        }
    }
    return 0.0;
}


TString root_file_name = "XYZ.root";
TString mc_root_file_name = "MC_XYZ.root";

// B0 Particle
void total_data_fit_Bd() {
    const int nbins_plot = 100; // Number of bins for the plotted data

    double min_signal = 5.1;
    double max_signal = 5.6;

    double xlow = 5.0;
    double xhigh = 5.9;
    const double bin_width_plot = (xhigh - xlow) / nbins_plot;

    double mc_mean = 5.28056;

    double mc_sigma1_rt = 0.01344;
    double mc_sigma2_rt = 0.03301;
    double mc_c1_rt = 0.6584;

    double mc_sigma1_wt = 0.07205;
    double mc_sigma2_left_wt = 0.10056;
    double mc_sigma2_right_wt = 0.23594;
    double mc_c1_wt = 0.3179;

    double mc_N_rt = 3239.96;
    double mc_N_wt = 3156.04;
    double mc_f1 = mc_N_rt / (mc_N_rt + mc_N_wt);

    // Load real data tree
    TFile* file = TFile::Open(root_file_name);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open real data file." << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get("tree");
    if (!tree) {
        std::cerr << "Error: TTree not found in file." << std::endl;
        return;
    }

    // RooFit: variable and data
    RooRealVar B_mass("B_mass", "B_mass", xlow, xhigh);
    B_mass.setRange("gaussRange", min_signal, max_signal);

    RooBinning mainBins(nbins_plot, xlow, xhigh);
    B_mass.setBinning(mainBins, "mainBins");

    RooDataSet dataset("dataset", "Unbinned dataset from TTree", tree, RooArgSet(B_mass));

    // Signal: RT (2 Gauss) + WT (Gauss + asymmetric Gauss)
    // Shared mean (FIXED to MC)
    RooRealVar mean("mean", "Mean", mc_mean);
    mean.setConstant(true);

    // Resolution Scale
    RooRealVar Cs("Cs", "Resolution scale", 1.0, 0.3, 3.0);

    // Right Tag (RT): two standard Gaussians, fixed sigmas and fixed fraction
    RooRealVar sigma1_rt("sigma1_rt", "RT sigma1", mc_sigma1_rt);          
    sigma1_rt.setConstant(true);
    
    RooRealVar sigma2_rt("sigma2_rt", "RT sigma2", mc_sigma2_rt);          
    sigma2_rt.setConstant(true);
    
    RooRealVar c1_rt("c1_rt", "RT fraction of G1", mc_c1_rt);      
    c1_rt.setConstant(true);

    RooProduct sigma1_rt_eff("sigma1_rt_eff", "sigma1_rt_eff", RooArgList(sigma1_rt, Cs));
    RooProduct sigma2_rt_eff("sigma2_rt_eff", "sigma2_rt_eff", RooArgList(sigma2_rt, Cs));

    RooGaussian rt_g1("rt_g1", "RT Gaussian 1", B_mass, mean, sigma1_rt_eff);
    RooGaussian rt_g2("rt_g2", "RT Gaussian 2", B_mass, mean, sigma2_rt_eff);
    RooAddPdf   rt_pdf("rt_pdf", "RT: two Gaussians", RooArgList(rt_g1, rt_g2), RooArgList(c1_rt));


    // Wrong Tag (WT): one standard Gaussian + one asymmetric Gaussian (bifurcated), fixed params
    RooRealVar sigma1_wt("sigma1_wt", "WT sigma1 (Gauss)", mc_sigma1_wt);        
    sigma1_wt.setConstant(true);

    RooRealVar sigma2L_wt("sigma2L_wt", "WT sigma2 left (bif)", mc_sigma2_left_wt);
    sigma2L_wt.setConstant(true);

    RooRealVar sigma2R_wt("sigma2R_wt", "WT sigma2 right (bif)", mc_sigma2_right_wt);  
    sigma2R_wt.setConstant(true);

    RooRealVar c1_wt("c1_wt","WT fraction of Gauss1",mc_c1_wt);            
    c1_wt.setConstant(true);

    RooProduct sigma1_wt_eff("sigma1_wt_eff", "sigma1_wt_eff", RooArgList(sigma1_wt, Cs));
    RooProduct sigma2L_wt_eff("sigma2L_wt_eff", "sigma2L_wt_eff", RooArgList(sigma2L_wt, Cs));
    RooProduct sigma2R_wt_eff("sigma2R_wt_eff", "sigma2R_wt_eff", RooArgList(sigma2R_wt, Cs));

    RooGaussian wt_g1("wt_g1", "WT Gaussian 1", B_mass, mean, sigma1_wt_eff);
    RooBifurGauss wt_g2("wt_g2", "WT asymmetric Gaussian", B_mass, mean, sigma2L_wt_eff, sigma2R_wt_eff);
    RooAddPdf wt_pdf("wt_pdf", "WT: Gauss + bifurGauss", RooArgList(wt_g1, wt_g2), RooArgList(c1_wt));

    // RT/WT mixture with floating fraction f1 (RT weight)
    RooRealVar f1("f1", "Fraction of RT in signal", mc_f1);  
    f1.setConstant(true);

    // Full signal shape (normalized): f1 * RT  +  (1-f1) * WT
    RooAddPdf signal("signal", "Signal = f1*RT + (1-f1)*WT", RooArgList(rt_pdf, wt_pdf), RooArgList(f1));

    RooRealVar Nsig("Nsig", "Signal Yield", 261, 0, 300000);
    RooExtendPdf signal_ext("signal_ext", "Extended Signal", signal, Nsig);

    /*
    // Background: Exponential model
    RooRealVar lambda("lambda", "Lambda", -2.72, -6.0, -0.1);
    RooExponential expo("expo", "Background", B_mass, lambda);

    RooRealVar Nbkg("Nbkg", "Background Yield", 1614, 0, 1700000);
    RooExtendPdf expo_ext("expo_ext", "Extended Exponential Background", expo, Nbkg);
    */

    RooRealVar a0("a0","a0",  0.0, -2.0, 2.0);
    RooRealVar a1("a1","a1",  0.0, -2.0, 2.0);
    RooRealVar a2("a2","a2",  0.0, -2.0, 2.0);
    RooRealVar a3("a3","a3",  0.0, -2.0, 2.0);
    RooRealVar a4("a4","a4",  0.0, -2.0, 2.0);
    RooChebychev bkg_pdf("bkg_pdf","Cheb4", B_mass, RooArgList(a0,a1,a2,a3,a4));

    RooRealVar Nbkg("Nbkg","Background Yield", 1.5e6, 0, 1e8);
    RooExtendPdf bkg_ext("bkg_ext","Extended background", bkg_pdf, Nbkg);

    // Combined model
    RooAddPdf model("model", "Signal + Background", RooArgList(signal_ext, bkg_ext));

    // Fit (Extended Maximum Likelihood Method)
    RooFitResult* result = model.fitTo(dataset, Save());
    
    // Compute background-only integrals in signal and sideband regions
    B_mass.setRange("signalRegion", min_signal, max_signal);
    B_mass.setRange("lowSideband", xlow, min_signal);
    B_mass.setRange("highSideband", max_signal, xhigh);

    double frac_bkg_signal = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal(); // Background in Signal Region
    double frac_bkg_low    = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("lowSideband"))->getVal(); // Background in Left Noise Region
    double frac_bkg_high   = bkg_pdf.createIntegral(B_mass, NormSet(B_mass), Range("highSideband"))->getVal(); // Background in Right Noise Region

    double total_bkg_yield = Nbkg.getVal(); // Total background

    double bkg_in_signal = total_bkg_yield * frac_bkg_signal; // Ammount of Noise in Signal Region 
    double bkg_out_signal = total_bkg_yield * (frac_bkg_low + frac_bkg_high); // Ammount of Noise in Sidebands
    double f_b = bkg_in_signal / bkg_out_signal; // Calculating F_b

    double frac_sig_in_signal = signal.createIntegral(B_mass, NormSet(B_mass), Range("signalRegion"))->getVal(); // Signal in Signal Region
    double sig_yield_in_region = Nsig.getVal() * frac_sig_in_signal; // Ammount of Signal in Signal Region


    // Opening and Checking MC File
    TFile *file_mc = TFile::Open(mc_root_file_name);
    if (!file_mc || file_mc->IsZombie()) {
        std::cerr << "Error: Could not open MC file." << std::endl;
        return;
    }

    TTree *treemc = nullptr;
    file_mc->GetObject("Bfinder/ntKstar", treemc);
    if (!treemc) {
        std::cerr << "Error: MC TTree not found!" << std::endl;
        return;
    }

    // Apply same cuts
    TString cut_mc = Form("Bnorm_svpvDistance_2D>2.5177 && ((%s)||(Bgen==41000)) && (%s) && (%s) && (%s)",
                        isMCsignal.Data(),
                        ACCcuts_ppRef.Data(),
                        SELcuts_ppRef.Data(),
                        TRGmatching.Data());

    int nbins_mc = 150;

    TH1F *hist_mc = new TH1F("hist_mc", "MC Bmass in Signal Region;Bmass [GeV/c^{2}];Entries", nbins_mc, min_signal, max_signal);

    treemc->Draw("Bmass >> hist_mc", cut_mc + Form(" && Bmass > %.4f && Bmass < %.4f", min_signal, max_signal), "goff");

    double mc_yield_in_signal = hist_mc->Integral();

    delete hist_mc;
    file_mc->Close();
    // --------------------------------------------------------------------------------------------------
    
    double f_s = sig_yield_in_region / mc_yield_in_signal; // Calculating F_s

    // Canvas with two pads
    TCanvas* c = new TCanvas("c", "Bd Fit with Pulls", 800, 800);
    c->Divide(1, 2);

    // Top pad (fit)
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetPad(0.0, 0.15, 1.0, 1.0);
    p1->SetBottomMargin(0.02);
    p1->Draw();
    p1->cd();

    RooPlot* frame = B_mass.frame(Range(xlow, xhigh), Bins(nbins_plot));

    dataset.plotOn(frame, Binning(B_mass.getBinning("mainBins")), MarkerStyle(20), MarkerSize(1.2), Name("data"), DataError(RooAbsData::Poisson));

    // Background component
    model.plotOn(frame, Components(bkg_ext), LineColor(kRed), LineStyle(kDashed), LineWidth(2), Name("background"));

    // RT and WT components (shape-only subcomponents of the signal)
    model.plotOn(frame, Components(rt_pdf), LineColor(kGreen + 2), LineStyle(kDashed), LineWidth(2), Name("RT"));
    model.plotOn(frame, Components(wt_pdf), LineColor(kMagenta + 2), LineStyle(kDashed), LineWidth(2), Name("WT"));

    // Total model
    model.plotOn(frame, LineColor(kBlue), LineWidth(2), Name("global"));

    frame->SetTitle("");
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetXaxis()->SetLabelSize(0);  // hide x labels on top pad
    frame->GetYaxis()->SetTitle(Form("Events / ( %.4f )", bin_width_plot));
    frame->Draw();

    // Vertical dashed lines at signal-region edges (heights taken from drawn curve)
    double y_low  = getYatMass(frame, min_signal);
    double y_high = getYatMass(frame, max_signal);

    TLine* line_low  = new TLine(min_signal, 0, min_signal, y_low);
    TLine* line_high = new TLine(max_signal, 0, max_signal, y_high);
    for (TLine* l : {line_low, line_high}) {
        l->SetLineColor(kBlack);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw("same");
    }

    // Chi2 after plotting on frame
    int nParams = result->floatParsFinal().getSize();
    double chi2 = frame->chiSquare("global", "data", nParams);


    // Legend (same place), on TOP pad
    p1->cd();
    TLegend* legend = new TLegend(0.48, 0.77, 0.88, 0.98);
    legend->SetTextFont(42);
    legend->SetTextSize(0.025);
    legend->SetBorderSize(1);
    legend->SetLineColor(kBlack);
    legend->SetFillStyle(1001);
    legend->SetFillColor(kWhite);
    legend->AddEntry(frame->findObject("data"), "Data (B^{0} )", "lep");
    legend->AddEntry(frame->findObject("background"), "Background (4th degree poly)", "l");
    legend->AddEntry(frame->findObject("RT"),         "RT Signal", "l");
    legend->AddEntry(frame->findObject("WT"),         "WT Signal", "l");
    legend->AddEntry(frame->findObject("global"),     "Total Fit", "l");
    legend->Draw();


    // TPaveText (same place), on TOP pad
    p1->cd();
    TPaveText* pave = new TPaveText(0.63, 0.20, 0.88, 0.77, "NDC");
    pave->SetTextAlign(12);
    pave->SetTextFont(42);
    pave->SetTextSize(0.025);
    pave->SetFillColor(0);
    pave->SetBorderSize(1);

    
    pave->AddText(Form("Mean (fixed) = %.5f", mean.getVal()));
    // RT fixed params
    pave->AddText(Form("#sigma_{1,RT} (fixed) = %.5f", sigma1_rt.getVal()));
    pave->AddText(Form("#sigma_{2,RT} (fixed) = %.5f", sigma2_rt.getVal()));
    pave->AddText(Form("c_{1,RT} (fixed) = %.3f", c1_rt.getVal()));
    // WT fixed params
    pave->AddText(Form("#sigma_{1,WT} (fixed) = %.5f", sigma1_wt.getVal()));
    pave->AddText(Form("#sigma_{2L,WT} (fixed) = %.5f", sigma2L_wt.getVal()));
    pave->AddText(Form("#sigma_{2R,WT} (fixed) = %.5f", sigma2R_wt.getVal()));
    pave->AddText(Form("c_{1,WT} (fixed) = %.3f", c1_wt.getVal()));
    // Floating RT/WT mixture
    pave->AddText(Form("f_{1} (fixed) = %.3f", f1.getVal()));
    pave->AddText(Form("C_{s} = %.4f #pm %.4f", Cs.getVal(), Cs.getError()));

    pave->AddText(Form("N_{sig} = %.1f #pm %.1f", Nsig.getVal(), Nsig.getError()));
    pave->AddText(Form("a_{0} = %.4f #pm %.4f", a0.getVal(), a0.getError()));
    pave->AddText(Form("a_{1} = %.4f #pm %.4f", a1.getVal(), a1.getError()));
    pave->AddText(Form("a_{2} = %.4f #pm %.4f", a2.getVal(), a2.getError()));
    pave->AddText(Form("a_{3} = %.4f #pm %.4f", a3.getVal(), a3.getError()));
    pave->AddText(Form("a_{4} = %.4f #pm %.4f", a4.getVal(), a4.getError()));
    pave->AddText(Form("N_{bkg} = %.1f #pm %.1f", Nbkg.getVal(), Nbkg.getError()));
    pave->AddText(Form("#chi^{2}/ndf = %.5f", chi2));
    pave->Draw();
    

    // f_b / f_s box (same place), on TOP pad
    p1->cd();
    TPaveText* pave_fb_fs = new TPaveText(0.36, 0.87, 0.48, 0.98, "NDC");
    pave_fb_fs->SetTextAlign(12);
    pave_fb_fs->SetTextFont(42);
    pave_fb_fs->SetTextSize(0.025);
    pave_fb_fs->SetFillColor(0);
    pave_fb_fs->SetBorderSize(1);
    pave_fb_fs->AddText(Form("f_{b} = %.3f", f_b));
    pave_fb_fs->AddText(Form("f_{s} = %.3f", f_s));
    pave_fb_fs->Draw();

    // Bottom pad (pulls)
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetPad(0.0, 0.0, 1.0, 0.15);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.25);
    p2->Draw();
    p2->cd();

    RooPlot* pullFrame = B_mass.frame();
    RooHist* pullHist = frame->pullHist("data", "global"); 
    pullHist->SetMarkerSize(0.6);
    pullFrame->addPlotable(pullHist, "XP");

    pullFrame->SetTitle("");
    pullFrame->GetYaxis()->SetTitle("Pull");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitleSize(0.10);
    pullFrame->GetYaxis()->SetTitleOffset(0.40);
    pullFrame->GetYaxis()->SetLabelSize(0.08);
    pullFrame->GetXaxis()->SetTitle("m_{J/#Psi K^{*}} [GeV/c^{2}]");
    pullFrame->GetXaxis()->SetTitleSize(0.10);
    pullFrame->GetXaxis()->SetTitleOffset(1.0);
    pullFrame->GetXaxis()->SetLabelSize(0.08);
    pullFrame->SetMinimum(-3.5);
    pullFrame->SetMaximum(3.5);
    pullFrame->Draw("AP");

    TLine* zeroLine = new TLine(xlow, 0, xhigh, 0);
    zeroLine->SetLineColor(kBlue);
    zeroLine->SetLineStyle(1);
    zeroLine->SetLineWidth(1);
    zeroLine->Draw("same");

    TString name_of_the_file = "Bd_Total_Fit_RT_and_WT.pdf";
    c->SaveAs(name_of_the_file);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << " " << std::endl;
    std::cout << "R3 = " << bkg_in_signal << " events" << std::endl;
    std::cout << "R1 + R2 = " << bkg_out_signal << " events" << std::endl;
    std::cout << "f_{b} = " << f_b << std::endl;

    std::cout << " " << std::endl;

    std::cout << "S_{data} " << sig_yield_in_region << " events" << std::endl;
    std::cout << "S_{MC} = " << mc_yield_in_signal << " events" << std::endl;
    std::cout << "f_{s} = " << f_s << std::endl;
    std::cout << " " << std::endl; 

    std::cout << "Fit saved to " << name_of_the_file << std::endl;

    delete c;
    delete tree;
    delete line_low;
    delete line_high;
}



void data_fit_Bzero() {
    total_data_fit_Bd();
}
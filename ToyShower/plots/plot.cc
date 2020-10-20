// g++ plot.cc -o plot `root-config --cflags --glibs` && ./plot

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "iostream"

using namespace std;

int main() {
  TCanvas* c = new TCanvas("c", "canvas", 1600, 1600);
  //TCanvas* clog = new TCanvas("clog", "canvaslog", 1600, 1600);
  //clog->SetLogy();

  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  TLegend *legend = new TLegend(0.17, 0.63, 0.45, 0.89); 
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);

  // -------------------------------------------------------------------------------
  TFile* fvetomaxdumb = new TFile("vetomaxdumb.root");
  TH1D* vetomaxdumbcounts = (TH1D*)fvetomaxdumb->Get("counts")->Clone();
  TH1D* vetomaxdumbtimes  = (TH1D*)fvetomaxdumb->Get("times")->Clone();
  TH1D* vetomaxdumbcalls  = (TH1D*)fvetomaxdumb->Get("calls")->Clone();

  TH1D* vetomaxdumbtimesdivide = (TH1D*)vetomaxdumbtimes->Clone("timesdivide");
  vetomaxdumbtimesdivide->Divide(vetomaxdumbtimes, vetomaxdumbcounts);
  TH1D* vetomaxdumbcallsdivide = (TH1D*)vetomaxdumbcalls->Clone("callsdivide");
  vetomaxdumbcallsdivide->Divide(vetomaxdumbcalls, vetomaxdumbcounts);

  vetomaxdumbtimesdivide->SetLineColor(1);
  vetomaxdumbtimesdivide->SetLineWidth(2);
  vetomaxdumbcallsdivide->SetLineColor(1);
  vetomaxdumbcallsdivide->SetLineWidth(2);
  vetomaxdumbcounts->SetLineColor(1);
  vetomaxdumbcounts->SetLineWidth(2);


  // -------------------------------------------------------------------------------
  TFile* fmaxvetodumb = new TFile("maxvetodumb.root");
  TH1D* maxvetodumbcounts = (TH1D*)fmaxvetodumb->Get("counts")->Clone();
  TH1D* maxvetodumbtimes  = (TH1D*)fmaxvetodumb->Get("times")->Clone();
  TH1D* maxvetodumbcalls  = (TH1D*)fmaxvetodumb->Get("calls")->Clone();

  TH1D* maxvetodumbtimesdivide = (TH1D*)maxvetodumbtimes->Clone("timesdivide");
  maxvetodumbtimesdivide->Divide(maxvetodumbtimes, maxvetodumbcounts);
  TH1D* maxvetodumbcallsdivide = (TH1D*)maxvetodumbcalls->Clone("callsdivide");
  maxvetodumbcallsdivide->Divide(maxvetodumbcalls, maxvetodumbcounts);

  maxvetodumbtimesdivide->SetLineColor(38);
  maxvetodumbtimesdivide->SetLineWidth(2);
  maxvetodumbcallsdivide->SetLineColor(38);
  maxvetodumbcallsdivide->SetLineWidth(2);
  maxvetodumbcounts->SetLineColor(38);
  maxvetodumbcounts->SetLineWidth(2);


  // -------------------------------------------------------------------------------
  TFile* fgenerateselectdumb = new TFile("generateselectdumb.root");
  TH1D* generateselectdumbcounts = (TH1D*)fgenerateselectdumb->Get("counts")->Clone();
  TH1D* generateselectdumbtimes  = (TH1D*)fgenerateselectdumb->Get("times")->Clone();
  TH1D* generateselectdumbcalls  = (TH1D*)fgenerateselectdumb->Get("calls")->Clone();

  TH1D* generateselectdumbtimesdivide = (TH1D*)generateselectdumbtimes->Clone("timesdivide");
  generateselectdumbtimesdivide->Divide(generateselectdumbtimes, generateselectdumbcounts);
  TH1D* generateselectdumbcallsdivide = (TH1D*)generateselectdumbcalls->Clone("callsdivide");
  generateselectdumbcallsdivide->Divide(generateselectdumbcalls, generateselectdumbcounts);

  generateselectdumbtimesdivide->SetLineColor(41);
  generateselectdumbtimesdivide->SetLineWidth(2);
  generateselectdumbcallsdivide->SetLineColor(41);
  generateselectdumbcallsdivide->SetLineWidth(2);
  generateselectdumbcounts->SetLineColor(41);
  generateselectdumbcounts->SetLineWidth(2);


  // -------------------------------------------------------------------------------
  TFile* fvetomax = new TFile("vetomax.root");
  TH1D* vetomaxcounts = (TH1D*)fvetomax->Get("counts")->Clone();
  TH1D* vetomaxtimes  = (TH1D*)fvetomax->Get("times")->Clone();
  TH1D* vetomaxcalls  = (TH1D*)fvetomax->Get("calls")->Clone();

  TH1D* vetomaxtimesdivide = (TH1D*)vetomaxtimes->Clone("timesdivide");
  vetomaxtimesdivide->Divide(vetomaxtimes, vetomaxcounts);
  TH1D* vetomaxcallsdivide = (TH1D*)vetomaxcalls->Clone("callsdivide");

  vetomaxcallsdivide->Divide(vetomaxcalls, vetomaxcounts);
  vetomaxtimesdivide->SetLineColor(46);
  vetomaxtimesdivide->SetLineWidth(2);
  vetomaxcallsdivide->SetLineColor(46);
  vetomaxcallsdivide->SetLineWidth(2);
  vetomaxcounts->SetLineColor(46);
  vetomaxcounts->SetLineWidth(2);

  // -------------------------------------------------------------------------------
  TFile* fmaxveto = new TFile("maxveto.root");
  TH1D* maxvetocounts = (TH1D*)fmaxveto->Get("counts")->Clone();
  TH1D* maxvetotimes  = (TH1D*)fmaxveto->Get("times")->Clone();
  TH1D* maxvetocalls  = (TH1D*)fmaxveto->Get("calls")->Clone();

  TH1D* maxvetotimesdivide = (TH1D*)maxvetotimes->Clone("timesdivide");
  maxvetotimesdivide->Divide(maxvetotimes, maxvetocounts);
  TH1D* maxvetocallsdivide = (TH1D*)maxvetocalls->Clone("callsdivide");
  maxvetocallsdivide->Divide(maxvetocalls, maxvetocounts);

  maxvetotimesdivide->SetLineColor(30);
  maxvetotimesdivide->SetLineWidth(2);
  maxvetocallsdivide->SetLineColor(30);
  maxvetocallsdivide->SetLineWidth(2);
  maxvetocounts->SetLineColor(30);
  maxvetocounts->SetLineWidth(2);

  // -------------------------------------------------------------------------------
  TFile* fgenerateselect = new TFile("generateselect.root");
  TH1D* generateselectcounts = (TH1D*)fgenerateselect->Get("counts")->Clone();
  TH1D* generateselecttimes  = (TH1D*)fgenerateselect->Get("times")->Clone();
  TH1D* generateselectcalls  = (TH1D*)fgenerateselect->Get("calls")->Clone();

  TH1D* generateselecttimesdivide = (TH1D*)generateselecttimes->Clone("timesdivide");
  generateselecttimesdivide->Divide(generateselecttimes, generateselectcounts);
  TH1D* generateselectcallsdivide = (TH1D*)generateselectcalls->Clone("callsdivide");
  generateselectcallsdivide->Divide(generateselectcalls, generateselectcounts);

  generateselecttimesdivide->SetLineColor(6);
  generateselecttimesdivide->SetLineWidth(2);
  generateselectcallsdivide->SetLineColor(6);
  generateselectcallsdivide->SetLineWidth(2);
  generateselectcounts->SetLineColor(6);
  generateselectcounts->SetLineWidth(2);

  THStack* hsTimes = new THStack();
  THStack* hsTimes2 = new THStack();
  THStack* hsCalls = new THStack();

  // -------------------------------------------------------------------------------
  hsTimes->Add(vetomaxdumbtimesdivide);
  hsCalls->Add(vetomaxdumbcallsdivide);
  legend->AddEntry(vetomaxdumbtimesdivide, "#font[132]{Veto-Max}");

  hsTimes->Draw("nostack HIST L");
  hsTimes->GetXaxis()->SetRangeUser(251, 500);
  hsTimes->GetXaxis()->SetTitle("#font[132]{Multiplicity}");
  hsTimes->GetYaxis()->SetTitle("#font[132]{Average Time [s]}");
  hsTimes->SetMaximum(1.5);
  hsTimes->Draw("nostack HIST L");
  legend->Draw();
  c->Print("times1.eps", "eps");

  hsCalls->Draw("nostack HIST L");
  hsCalls->GetXaxis()->SetRangeUser(251, 500);
  hsCalls->GetXaxis()->SetTitle("#font[132]{Multiplicity}");
  hsCalls->GetYaxis()->SetTitle("#font[132]{Function Calls}");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->Print("calls1.eps");

  // -------------------------------------------------------------------------------
  hsTimes->Add(maxvetodumbtimesdivide);
  hsCalls->Add(maxvetodumbcallsdivide);
  legend->AddEntry(maxvetodumbtimesdivide, "#font[132]{Max-Veto}"); 
  hsTimes->Draw("nostack HIST L");
  legend->Draw();
  c->Print("times2.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->Print("calls2.eps");

  // -------------------------------------------------------------------------------
  hsTimes->Draw("nostack HIST L");
  legend->Draw();
  c->Print("times3.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  hsCalls->SetMaximum(4e7);
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(1);
  c->Print("calls3.eps");
  
  // -------------------------------------------------------------------------------
  hsTimes->Add(generateselectdumbtimesdivide);
  hsCalls->Add(generateselectdumbcallsdivide);
  legend->AddEntry(generateselectdumbtimesdivide, "#font[132]{Generate-Select}"); 
  hsTimes->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(0);
  c->Print("times4.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(1);
  c->Print("calls4.eps");

  // -------------------------------------------------------------------------------
  hsTimes->Add(vetomaxtimesdivide);   
  hsCalls->Add(vetomaxcallsdivide);
  legend->AddEntry(vetomaxtimesdivide, "#font[132]{Veto-Max (Scale Saving)}"); 
  hsTimes->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(0);
  c->Print("times5.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(1);
  c->Print("calls5.eps");

  // -------------------------------------------------------------------------------
  hsTimes2->Add(vetomaxtimesdivide);
  hsTimes2->Draw("nostack HIST L");
  hsTimes2->GetXaxis()->SetRangeUser(251, 500);
  hsTimes2->GetXaxis()->SetTitle("#font[132]{Multiplicity}");
  hsTimes2->GetYaxis()->SetTitle("#font[132]{Average Time [s]}");
  hsTimes2->SetMaximum(0.021);
  hsTimes2->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(0);
  c->Print("times6.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(1);
  c->Print("calls6.eps");

  // -------------------------------------------------------------------------------
  hsTimes2->Add(maxvetotimesdivide);
  hsCalls->Add(maxvetocallsdivide);
  legend->AddEntry(maxvetotimesdivide, "#font[132]{Max-Veto (Scale Saving)}"); 
  hsTimes2->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(0);
  c->Print("times7.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(1);
  c->Print("calls7.eps");

  // -------------------------------------------------------------------------------
  hsTimes2->Add(generateselecttimesdivide);
  hsCalls->Add(generateselectcallsdivide);
  legend->AddEntry(generateselecttimesdivide, "#font[132]{Generate-Select (Roulette Wheel)}"); 
  hsTimes2->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(0);
  c->Print("times8.eps", "eps");
  hsCalls->Draw("nostack HIST L");
  legend->Draw();
  c->SetLogy(1);
  c->Print("calls8.eps");
}
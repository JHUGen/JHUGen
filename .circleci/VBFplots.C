#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <cstdlib>
using namespace std;

double pi = TMath::Pi();
const int nvariables = 15;
const TString variables[nvariables] = {"Gencosthetastar_VBF", "GenhelcosthetaV1_VBF", "GenhelcosthetaV2_VBF", "Genhelphi_VBF", "GenphistarV1_VBF",
                                       "GenQ_V1", "GenQ_V2", "GenDijetMass", "GenDRjet",
                                       "leadingpT", "subleadingpT", "leadingeta", "subleadingeta", "leadingphi", "subleadingphi"};
const int ntreevariables = nvariables-6;

const int nfiles = 4;
vector<TString> files[nfiles];
TString names[nfiles];
bool filesfilled = false;
void setupfiles();


double mins[nvariables] = {-1, -1, -1, -pi, -pi,
                           0, 0, 0, 0,
                           0, 0, -6, -6, -pi, -pi};
double maxes[nvariables] = {1, 1, 1, pi, pi,
                            1000, 1000, 2000, 12,
                            250, 250, 6, 6, pi, pi};

bool withptcut = false;

void VBFplots()
{
    TString dir = "$CIRCLE_ARTIFACTS/VBF/";
    if (withptcut) dir += "/withptcut";

    setupfiles();
    THStack *h[nvariables];
    TH1F *hh[nvariables][nfiles];
    for (int i = 0; i < nvariables; i++)
    {
        h[i] = new THStack(variables[i], variables[i]);
        for (int j = 0; j < nfiles; j++)
        {
            hh[i][j] = new TH1F(TString(files[j][0]).ReplaceAll("/","") += variables[i], "h", 100, mins[i], maxes[i]);
            hh[i][j]->SetLineColor(j+1);
        }
    }

    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->SetFillStyle(0);

    for (int j = 0; j < nfiles; j++)
    {
        TChain *t = new TChain("SelectedTree", "SelectedTree");
        for (unsigned int k = 0; k < files[j].size(); k++)
        {
            int code = t->Add(files[j][k]);
            if (code < 1)
                cout << files[j][k] << endl;
        }

        float x[nvariables];
        for (int i = 0; i < ntreevariables; i++)
            t->SetBranchAddress(variables[i], &(x[i]));

        vector<double>  *pt = 0;
        vector<double> *eta = 0;
        vector<double> *phi = 0;
        t->SetBranchAddress("GenAssociatedParticlePt", &pt);
        t->SetBranchAddress("GenAssociatedParticleEta", &eta);
        t->SetBranchAddress("GenAssociatedParticlePhi", &phi);

        long length = t->GetEntries();
        if (length == 0) continue;
        for (int l = 0; l < length; l++)
        {
            t->GetEntry(l);
            bool passptcut = true;
            for (unsigned int i = 0; i < pt->size(); i++)
                if (pt->at(i) < 15)
                    passptcut = false;

            if (pt->at(0) <= pt->at(1))
                cout << "bad pt " << l << " " << pt->at(0) << " " << pt->at(1) << endl;
            x[9] = pt->at(0);
            x[10] = pt->at(1);
            x[11] = eta->at(0);
            x[12] = eta->at(1);
            x[13] = phi->at(0);
            x[14] = phi->at(1);

            if (passptcut || !withptcut)
                for (int i = 0; i < nvariables; i++)
                    hh[i][j]->Fill(x[i]);

            if ((l+1) % 10000 == 0 || l+1 == length)
                cout << l+1 << " / " << length << endl;
        }
        for (int i = 0; i < nvariables; i++)
        {
            hh[i][j]->Scale(1.0/hh[i][j]->Integral());
            h[i]->Add(hh[i][j]);
        }

        leg->AddEntry(hh[0][j], names[j], "l");
    }

    TCanvas *c1 = new TCanvas();
    gSystem->mkdir(dir, true);
    for (int i = 0; i < nvariables; i++)
    {
        h[i]->Draw("nostack");
        leg->Draw();
        c1->SaveAs(TString(dir) += TString("/") += TString(variables[i]) += ".png");
        c1->SaveAs(TString(dir) += TString("/") += TString(variables[i]) += ".eps");
        c1->SaveAs(TString(dir) += TString("/") += TString(variables[i]) += ".root");
        c1->SaveAs(TString(dir) += TString("/") += TString(variables[i]) += ".pdf");
    }
}




void setupfiles()
{
    if (filesfilled) return;

    files[0].push_back("$CIRCLE_ARTIFACTS/VBF/SM.root");
    names[0] = "VBF SM";

    files[1].push_back("$CIRCLE_ARTIFACTS/VBF/PS.root");
    names[1] = "VBF PS";

    files[2].push_back("$CIRCLE_ARTIFACTS/VBF/a2.root");
    names[2] = "VBF a2";

    files[3].push_back("$CIRCLE_ARTIFACTS/HJJ/SM.root");
    names[3] = "HJJ SM";

    filesfilled = true;
};

{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
   gROOT->Reset();
#include "Riostream.h"

   ifstream in;
// we assume a file basic.dat in the current directory
// this file has 3 columns of float data
//   in.open("Events_11-15.txt_corr");
   in.open("Output.txt");
   
   Double_t pt,event,dphi, ptnt , phnt , ptt , pht , pth , phh , ppnt , ppt , pph, ptcostht, ptcosthh, ptsintht, ptsinthh,dptcosthth;
   Int_t nlines = 0;
   TFile *f = new TFile("pythia4.root","RECREATE");
   TH1D *dphithLQ = new TH1D("dphithLQ","\dphi (tau-string) Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,-180,180);
   TH1D *ptcosthtLQ = new TH1D("ptcosthtLQ","px=pt Cos(theta) tau Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,0,100);
   TH1D *ptcosthhLQ = new TH1D("ptcosthhLQ","px=pt Cos(theta) string Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,0,100);
   TH1D *ptsinthtLQ = new TH1D("ptsinthtLQ","py=pt Sin(theta) tau Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,0,100);
   TH1D *ptsinthhLQ = new TH1D("ptsinthhLQ","py=pt Sin(theta) string Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,0,100);

//   TH1D *dptcosththLQ = new TH1D("dptcosththLQ","\dpx (tau-string) Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,0,100);

//  TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","pt");

   while (1) {
      in >> event  >> ptnt >> phnt >> ptt >> pht >> pth >> phh >> ppnt >> ppt >> pph;
      if (!in.good()) break;
//      if (nlines < 5) printf("pt=%8f",pt);
      dphi=pht-phh;
      dphithLQ->Fill(dphi);
      dphithLQ.GetXaxis()->SetTitle("\dphi");
      dphithLQ.GetYaxis()->SetTitle("");
//      ntuple->Fill(pt);

     ptcostht=ptt*cos(ppt*3.14159/180);
     ptcosthh=pth*cos(pph*3.14159/180);
     ptsintht=ptt*sin(ppt*3.14159/180);
     ptsinthh=pth*sin(pph*3.14159/180);
     
     dptcosthth=ptcostht-ptcosthh;
     ptcosthtLQ->Fill(ptcostht);
     ptcosthhLQ->Fill(ptcosthh);
      ptcosthtLQ.GetXaxis()->SetTitle("\px tau(GeV)");
      ptcosthhLQ.GetXaxis()->SetTitle("\px string(GeV)");
     ptsinthtLQ->Fill(ptsintht);

     ptsinthhLQ->Fill(ptsinthh);
     
//     dptcosththLQ->Fill(dptcosthth);
      nlines++;
   }
   printf(" found %d pointsn",nlines);

   in.close();

   f->Write();
}


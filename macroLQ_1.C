{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
   gROOT->Reset();
#include "Riostream.h"

   ifstream in;
// we assume a file basic.dat in the current directory
// this file has 3 columns of float data
   in.open("Events_11-16.txt");
   
   Double_t pt,event,ptx,phx,pty,phy,ptz,phz;
   Int_t nlines = 0;
   TFile *f = new TFile("pythia4.root","RECREATE");
   TH1D *ptxLQ = new TH1D("ptxLQ","\dphi (tau-string) Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,-50,50);
   TH1D *phxLQ = new TH1D("phxLQ","px=pt Cos(theta) tau Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,-50,50);
   TH1D *ptyLQ = new TH1D("ptyLQ","px=pt Cos(theta) string Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,-50,50);
   TH1D *phyLQ = new TH1D("phyLQ","py=pt Sin(theta) tau Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,-50,50);



   while (1) {
      in >> event  >> ptx >> phx >> pty >> phy >> ptz >> phz;
      if (!in.good()) break;

      ptxLQ->Fill(ptx);
      ptxLQ.GetXaxis()->SetTitle("px tau(GeV)");
      ptxLQ.GetYaxis()->SetTitle("");

     phxLQ->Fill(phx);
     phxLQ.GetXaxis()->SetTitle("\px string(GeV)");

     ptyLQ->Fill(pty);
     ptyLQ.GetXaxis()->SetTitle("\py tau(GeV)");

     phyLQ->Fill(phy);
     phyLQ.GetXaxis()->SetTitle("\py string(GeV)");
      nlines++;
   }
   printf(" found %d pointsn",nlines);

   in.close();

   f->Write();
}


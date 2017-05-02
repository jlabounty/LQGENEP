{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
   gROOT->Reset();
#include "Riostream.h"
#include <stdio.h>
#include <math.h>

#define PI 3.14159265

   ifstream in;
// we assume a file basic.dat in the current directory
// this file has 3 columns of float data
   in.open("Events_11-29.txt");
   
   Double_t pt,event, ptMiss,px,py,pz,phi_nutau, phi_tau,phi_had,phimiss,ptmissphi;
   Int_t nlines = 0;
   TFile *f = new TFile("pythia4.root","RECREATE");
   TH1D *ptMissLQ = new TH1D("ptMissLQ","Missing mass Spectrum for 500 LQ ep events with 10 on 250 Gev",50,0,80);
   TH1D *ptmissphiLQ = new TH1D("ptmissphiLQ","\dphi_miss Spectrum for 1000 LQ ep events with 10 on 250 Gev",50,-180,180);

   while (1) {
      in >> event  >> ptMiss >> px >> py >> pz >> phi_nutau >> phi_tau >> phi_had;
      ptmissphi=atan(py/px)*180/PI;
      phimiss=abs((atan(py/px)*180/PI)-phi_tau); 
      if (!in.good()) break;
      ptMissLQ->Fill(ptMiss);
      ptMissLQ.GetXaxis()->SetTitle("\p_T Miss(GeV)");

      ptmissphiLQ->Fill(ptmissphi);
      ptmissphiLQ.GetXaxis()->SetTitle("\phi Miss");

           
      nlines++;
   }
   printf(" found %d pointsn",nlines);

   in.close();

   f->Write();
}


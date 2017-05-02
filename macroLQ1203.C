{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
   gROOT->Reset();
#include "Riostream.h"
#include <stdio.h>
#include <math.h>

#define PI 3.14159265

   ifstream in;
//       open(unit=30,file='DATA/27on820_miss.dat'
//       open(unit=31, file='DATA/27on820_jet.dat'
//       open(unit=32, file='DATA/27on820_e.dat'
//       open(unit=33, file='DATA/27on820_tau.dat'
//       open(unit=34, file='DATA/27on820_nue.dat'

//        write(30,10)lqgpar1(4),ptid,pxid,pyid,pzid
//        write(31,11)lqgpar1(4),ptj,pxj,pyj,phj,ppj
//        write(32,11)lqgpar1(4),pte,pxe,pye,phe,ppe
//        write(33,11)lqgpar1(4),ptt,pxt,pyt,pht,ppt
//        write(34,11)lqgpar1(4),ptne,pxne,pyne,phne,ppne

   in.open("fort.8");
   
   Double_t event,pxid,pyid,pxt,pxt,pyt,pht,ppt,ptmissx,ptmissy,ptmiss2,ptmiss;
   Int_t nlines = 0;
   TFile *f = new TFile("pythia4.root","RECREATE");
   TH1D *ptMissLQ = new TH1D("ptMissLQ","Missing mass Spectrum for 500K LQ ep events with 27.5 on 820 Gev",40,0,80);

   while (1) {
      in >> event>>pxid>>pyid>>pxt>>pxt>>pyt>>pht>>ppt>>ptmissx>>ptmissy>>ptmiss2>>ptmiss;
      if (!in.good()) break;
      ptmissx=(pxid)*(pxid);
      ptmissy=(pyid)*(pyid);


      ptmiss2=(ptmissx+ptmissy);
      ptmiss=pow(ptmiss2,0.5);
      ptMissLQ->Fill(ptmiss);
      ptMissLQ.GetXaxis()->SetTitle("\p_T Miss(GeV)");

           
      nlines++;
   }
   printf(" found %d pointsn",nlines);

   in.close();

   f->Write();
}


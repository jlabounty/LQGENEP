*-----------------
* File: LQguser.f
*-----------------
*
      Program LQGUSER
*******************************************
* Generates 1000 events of the subprocess * 
*         e+ qi --> S_1/2L --> mu+ qj     *    
*           at the HERA energies          *
*                                         *
* LQ mass = 200 GeV                       *
* LQ couplings = 0.3                      *
* first generation quark involved         *
******************************************* 
      Implicit None
*
      double precision mass,cross_max,
     >  ptnt,phnt,ptt,pht,pth,phh,
     >  ppnt,ppt,pph,
     >  ptx,pty,ptz,phx,phy,phz,
     >  ptid,phid,ppid,pxid,pyid,pzid,lam

      integer lqtype,Nevt,qi,qj,iproc,id,
     > Nevtdown,Nevtup
      logical first
      data first /.true./
      save first
*
C...LQGENEP run setup parameters
      double precision BEAMPAR,LQGPAR3
      integer LQGPAR1,LQGPAR2
      COMMON/LQGDATA/BEAMPAR(3),LQGPAR1(10),LQGPAR2(10),LQGPAR3(20)

C...LQGENEP event informations
      double precision LQGKPAR,LQGDDX
      integer LQGPID
      COMMON/LQGEVT/LQGKPAR(3),LQGDDX(3),LQGPID(3)
*
      integer NeV,i
*
c       Nevt=10000
       OPEN(UNIT=22,FILE='inputfile',STATUS='OLD')

         Read(22,*)Nevt,Mass,beampar(2),beampar(3),
     >    qi,qj
         close(22)




c      Mass=1936.D0
c      qi=1
c      qj=1
      iproc=0
c      cross_max=2.d-5
      cross_max=7.d-15
*
* Charge of lepton beam (+1/-1 positron/electron) 
c      beampar(1)=+1.  ! positron
      beampar(1)=-1.  ! electron
* Energy of lepton beam 
c      beampar(2)=10.     ! EIC1
c      beampar(2)=20.    ! EIC2 
c       beampar(2)=27.5  ! HERA
* Energy of proton beam 
c      beampar(3)=250.   ! EIC1
c      beampar(3)=325.   ! EIC2
c      beampar(3)=820.   ! HERA
 
* Number of events to generate 
      lqgpar1(1)=Nevt
* Number of events which will be fully listed in the output 
      lqgpar1(2)=0
* Histogram flag (if 1 some histograms will be produced) 
      lqgpar1(3)=1
* Current event number (first generated event will be lqgpar1(4)+1) 
      lqgpar1(4)=0
* Pythia process number (should be > 400, if = 0 it will be set
* to 401, first value available for external processes)  
      lqgpar1(5)=iproc
* LQ type
      lqgpar2(1)=14
* generation of the input quark in the s-channel process 
      lqgpar2(2)=qi
* generation of the output quark in the s-channel process 
      lqgpar2(3)=qj
* generation of the output lepton      
      lqgpar2(4)=3
* LQ mass in GeV
      lqgpar3(1)=Mass
* Initial state s-channel coupling
      lqgpar3(2)=0.3
* Final state s-channel coupling (in case of process eq -> LQ -> eq
*                               the two couplings should be the same) 
      lqgpar3(3)=0.3
* x range low limit
      lqgpar3(4)=0.011
* x range high limit
      lqgpar3(5)=1.
* y range low limit
      lqgpar3(6)=0.1
* y range high limit
      lqgpar3(7)=1.       
* minimum allowed Q^2 in Gev^2
      lqgpar3(8)=1000.
* Parton distribution type according to PDFLIB      
      lqgpar3(9)=1.       
* Parton distribution group according to PDFLIB      
      lqgpar3(10)=4.       
* Parton distribution set according to PDFLIB      
      lqgpar3(11)=32.       
* Max cross section
      lqgpar3(12)=cross_max
* Eventually switch off initial state QCD and QED radiation
* setting to 0 the following Pythia parameters
c     MSTP(61)=0
c     MSTP(71)=0
* Eventually switch off multiple interaction
c     MSTP(81)=0
* Eventually switch off fragmentation and decay
c     MSTP(111)=0

* LQGENEP Initialization
       call LQGENEP(Nevt,0)


      Nev=lqgpar1(1)

* LQGENEP generation loop
      do i=1,Nevt
cc       if (i.NE.Nevt) then
c        if ((i.GT.(Nevt-1000)).and.
c     >    (i.LE.Nevt)) then
cc       CALL PYEVNT
cc      else
       call LQGENEP(Nevt,1)
cc      endif
      enddo       

* LQGENEP termination
       call LQGENEP(Nevt,2)
      stop
      end

*-----------------
* File: lqgenep.f
*-----------------
*
           subroutine LQGENEP(Nevt,flag)
C------------------------------------------
C...Main program for leptoquark generation 
C...in electron-proton scattering
C------------------------------------------
C...All real arithmetic in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      Integer flag, itau,id,ii,ij

C...LQGENEP run setup parameters
      double precision BEAMPAR,LQGPAR3,
     >  ptnt,phnt,ptt,pht,pth,phh,
     >  ptx,pty,ptz,phx,phy,phz,
     > ptid, phid,ppid,pxid,pyid,pzid

      integer LQGPAR1,LQGPAR2
      COMMON/LQGDATA/BEAMPAR(3),LQGPAR1(10),LQGPAR2(10),LQGPAR3(20)

C...LQGENEP event informations
       double precision LQGKPAR,LQGDDX
       integer LQGPID
       COMMON/LQGEVT/LQGKPAR(3),LQGDDX(3),LQGPID(3)

C...Pythia declarations. 
C...Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers
C...(left- and righthanded SUSY, excited fermions).
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)

C...EXTERNAL statement links PYDATA on most machines.
      EXTERNAL PYDATA
*
C...Pythia Commonblocks.
C...The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...Process information.
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
C...Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,
     >/PYINT2/,/PYMSSM/
C...Local Array.
      DIMENSION NCHN(12),QVEC(4)
      DATA NCHN/12*0/

C...Internal used common
C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar

C# LQGpdfC.inc #
      character*20 parm(20)
      double precision pdfsf(20)
      common /LQGpdfC/ pdfsf 
      
C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc # 
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2

C# LQGout.inc #
      double precision DXSEC(3),pvalence
      integer q_o,q_s,genproc,genprtyp,sch,uch,gproc(8)
      common /LQGout/ DXSEC,pvalence,q_o,q_s,genproc,genprtyp,
     >sch,uch,gproc


C...processes
      integer ibea
      Character*9 chbea(2)
      Character*12 chprod(2,3)
      data chbea /'e+ qi -> ','e- qi -> '/
      data chprod /' -> e+ qj',' -> e- qj',
     >             ' -> mu+ qj',' -> mu- qj',
     >             ' -> tau+ qj',' -> tau- qj'/
C...LQ type
      Character*7 LQCHA(14)
      DATA LQCHA /'S_0L','S_0R','~S_0R','S_1L',
     >            'V_1/2L','V_1/2R','~V_1/2L',
     >            'V_0L','V_0R','~V_0R','V_1L',
     >            'S_1/2L','S_1/2R','~S_1/2L'/                        

C...Local declarations
      Real pxtot,pytot,pztot,etot,ecmtot,
     >     pxsum,pysum,pzsum,esum,ecmsum,
     >     pxlow,pxhi,pylow,pyhi,pzlow,pzhi,elow,ehi,ecmlow,ecmhi 

C...ASCII Output File
      INTEGER asciiN
      PARAMETER (asciiN=29)
      CHARACTER*256 out_file

C-----------------------------------------------------------------

      Integer Nwds_HBOOK
      Parameter (Nwds_HBOOK=100000)
      Real          HMEM
      Common /PAWC/ HMEM(Nwds_HBOOK)

C Give output file name

      out_file='./TestOut.txt'
C      out_file='./outdir/TestOut.txt'

C...FLAG=0 -> First section: inizialization  
      If(flag.eq.0)then

C.. LQGENEP banner
       call LQGBAN

C.. Hbook inizialization       
       if(LQGPAR1(3).gt.0)Call HLIMIT(Nwds_HBOOK)
       if(LQGPAR1(3).lt.0)Call HLIMIT(-Nwds_HBOOK)

C...beams properties.
       echar=beampar(1)
       Ebeam=beampar(2)
       Pbeam=beampar(3)
       S=4.d0*Ebeam*Pbeam

C...LQ properties
       MLQ=LQGPAR3(1)
       G1=LQGPAR3(2)
       G2=LQGPAR3(3)
       LQTYPE=LQGPAR2(1)       

C... outcoming lepton
       l_o=LQGPAR2(4)

C... incoming and outcoming quark generation
       q_i=LQGPAR2(2)       
       q_j=LQGPAR2(3)

C... kinematic ranges
       xmin=LQGPAR3(4)              
       xmax=LQGPAR3(5)              
       ymin=LQGPAR3(6)              
       ymax=LQGPAR3(7) 
       Zmin=0.d0
       Zmax=1.d0
       Q2min=LQGPAR3(8)             

C... print LQGENEP generation run settings
       call LQGPRI1

C... structure function
       parm(1)='NPTYPE'
       parm(2)='NGROUP'
       parm(3)='NSET'
       pdfsf(1)= LQGPAR3(9)
       pdfsf(2)= LQGPAR3(10) 
       pdfsf(3)= LQGPAR3(11)
       call PDFSET(PARM,pdfsf)

C... Pythia initialization 
       P(1,1)=0D0
       P(1,2)=0D0
       P(1,3)=-Ebeam
       P(2,1)=0D0
       P(2,2)=0D0
       P(2,3)=Pbeam

C...Evaluate limits for total momentum and energy histos 
       if(LQGPAR1(3).gt.0)then
        pxtot=sngl(P(1,1)+P(2,1))
        pytot=sngl(P(1,2)+P(2,2))
        pztot=sngl(P(1,3)+P(2,3))
        etot=sngl(dabs(P(1,3))+dabs(P(2,3)))
        ecmtot=sqrt(etot*etot-(pxtot*pxtot+pytot*pytot+pztot*pztot))
        if(pxtot.gt.0)then
         pxlow=pxtot-0.01*pxtot
         pxhi=pxtot+0.01*pxtot
        else
         pxlow=pxtot-1.
         pxhi=pxtot+1.
        endif
        if(pytot.gt.0)then
         pylow=pytot-0.01*pytot
         pyhi=pytot+0.01*pytot
        else
         pylow=pytot-1.
         pyhi=pytot+1.
        endif
        if(pztot.gt.0)then
         pzlow=pztot-0.01*pztot
         pzhi=pztot+0.01*pztot
        else
         pzlow=pztot-1.
         pzhi=pztot+1.
        endif
        if(etot.gt.0)then
         elow=etot-0.01*etot
         ehi=etot+0.01*etot
        else
         elow=etot-1.
         ehi=etot+1.
        endif                        
        if(ecmtot.gt.0)then
         ecmlow=ecmtot-0.01*ecmtot
         ecmhi=ecmtot+0.01*ecmtot
        else
         ecmlow=ecmtot-1.
         ecmhi=ecmtot+1.
        endif        
       endif
C...Initialize Pythia
       if(LQGPAR1(5).eq.0)then
        Isub=401
       else
        Isub=LQGPAR1(5)
       endif
       sch=0.d0
       uch=0.d0
       call vzero(gproc,8)
       sigmax=LQGPAR3(12)
       if(echar.gt.0)then
        ibea=1
       else
        ibea=2
       endif
       CALL PYUPIN(ISUB,
     >      CHBEA(ibea)//LQCHA(LQTYPE)//CHPROD(ibea,l_o),sigmax)
       MSEL=0
       MSUB(ISUB)=1
*
       if(beampar(1).GT.0)then
        CALL PYINIT('USER','e+','p',0D0)
       else
        CALL PYINIT('USER','e-','p',0D0)
       endif


       if(LQGPAR1(3).ne.0)then
C...Book histos.
        call hropen(69,'lqgenep','lqgenep.histo','N',1024,ierr) 
        CALL hbook1(1000,'x gen',50,sngl(xmin),sngl(xmax),0.)
        CALL hbook1(1001,'x gen s-ch.',50,sngl(xmin),sngl(xmax),0.)
        CALL hbook1(1002,'x gen u-ch',50,sngl(xmin),sngl(xmax),0.)
        CALL hbook1(2000,'y gen',50,sngl(ymin),sngl(ymax),0.) 
        CALL hbook1(2001,'y gen s-ch',50,sngl(ymin),sngl(ymax),0.) 
        CALL hbook1(2002,'y gen u-ch',50,sngl(ymin),sngl(ymax),0.) 
        CALL hbook1(3000,'Q2 gen',50,0.,6.,0.) 
        call hbook1(5001,'sum px',100,pxlow,pxhi,0.)
        call hbook1(5002,'sum py',100,pylow,pyhi,0.) 
        call hbook1(5003,'sum pz',100,pzlow,pzhi,0.)
       call hbook1(5004,'sum e',100,elow,ehi,0.) 
        call hbook1(5000,'center of mass energy',
     >  100,ecmlow,ecmhi,0.) 
       endif
       write(6,*)
      endif

C...Open File
      OPEN(asciiN, file=out_file)

C      WRITE(29,*)'============================================'
      Integer HeaderTest = 0
      if(HeaderTest .LE. 1) then
      WRITE(29,*)' Pythia Output File'
      WRITE(29,*)'============================================'
      WRITE(29,*)'I, ievent, genevent, subprocess, nucleon,
     & targetparton, xtargparton, beamparton, xbeamparton,
     & thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu,
     & leptonphi, s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R,
     &  sigma_rad, SigRadCor, EBrems, photonflux, nrTracks'
      WRITE(29,*)'============================================'
C      WRITE(29,30) 'I', 'KS', 'KF(ID)', 'ORIG'
C30    FORMAT(A5, A5, A7, A5)
C      WRITE(29,30) 'I', 'KS', 'KF(ID)', 'ORIG', 'px'
C     &           , 'py', 'pz', 'E', 'm'
C30    FORMAT(A5, A5, A7, A5, A10, A10, A10, A10, A10)
      WRITE(29,30) 'I', 'KS', 'KF(ID)', 'ORIG', 'D1', 'D2', 'px'
     &           , 'py', 'pz', 'E', 'm', 'vx', 'vy', 'vz'
30    FORMAT(A5, A5, A7, A5, A5, A5, A10, A10, A10, A10, A10, 
     &           A10, A10, A10)
      WRITE(29,*)'============================================'
      endif
      HeaderTest = 999
C-----------------------------------------------------------------
C...FLAG=1 -> Second section: event generation  

      if(flag.eq.1)Then
       CALL PYEVNT
c       print*,"The no. of event is",LQGPAR1(4)
       CALL PYHEPC(1)

C...s-u channel
       if(genproc.eq.1)sch=sch+1     
       if(genproc.eq.2)uch=uch+1     

C...process type
      gproc(genprtyp)=gproc(genprtyp)+1

C...Fill event informations common
       LQGKPAR(1)=X
       LQGKPAR(2)=Y
       LQGKPAR(3)=Q2
       LQGDDX(1)=(DXSEC(2)+DXSEC(3))*1.d-9
       LQGDDX(2)=DXSEC(3)*1.d-9
       LQGDDX(3)=DXSEC(2)*1.d-9
       LQGPID(1)=q_s
       LQGPID(2)=q_o
       if(genproc.eq.1)then
        LQGPID(3)=1
       elseif(genproc.eq.2)then
        LQGPID(3)=2
       endif

**swadhin: HADRONIC TAU DECAY
*Here particle that are not decayed are called (except neutrinos) and their pT are added. The sum is called pT Miss.


        do 60 J=1,N

          if ((K(J,1).EQ.11).and.
     >        (K(J,2).EQ.15)) then      ! find tau and get it's line number
           idt=J
           idtd=K(J,4)                  ! tau decay's to what line number
          endif

          do 45 l=1,N              
          if  (l.EQ.idt) then      
                                    
             ptt=PYP(l,10)
             pxt=P(l,1)
             pyt=P(l,2)
             pht=PYP(l,16)
             ppt=PYP(l,14)
           endif
45        enddo

          if ((K(J,1).EQ.1).and.
     >        (K(J,2).EQ.16)
     >    .and.(K(J,3).EQ.idt)) then      ! find tau neutrino and get it's line number
           idtnu=J
          endif

          if ((J.EQ.idtd)      ! line number is the decay of tau
     >   .and.(K(J,2).NE.-12)  ! this decay is not nu_ebar
     >   .and.(K(J,2).NE.-14)) then ! this decay is not nu_mubar

          ptid=0.d0
          pxid=0.d0
          pyid=0.d0
          pzid=0.d0

          do 50 I=1,N

          if ((K(I,1).LT.11)    ! Anything that doesn't decay = Final Product
     >   .and.(K(I,2).NE.12)  ! Final Product NOT AN Electron neutrino
     >   .and.(K(I,2).NE.14)  ! Final Product NOT AN Muon neutrino
     >   .and.(K(I,2).NE.16)  ! Final Product NOT AN Tau neutrino
     >   .and.(K(I,2).NE.-12) ! Final Product NOT AN Electron neutrino
     >   .and.(K(I,2).NE.-14) ! Final Product NOT AN  Muon neutrino
     >   .and.(K(I,2).NE.-16)) then ! Final Product NOT A Tau neutrino

          ptid=ptid+PYP(I,10)
          pxid=pxid+P(I,1)
          pyid=pyid+P(I,2)
          pzid=pzid+P(I,3)
          endif
50        enddo

        !!Delta Phi < 20 cut
 
        if (pyid.GT.0.) then
        phimiss=acos(pxid/sqrt(pxid*pxid+pyid*pyid))*180/(3.14159)
        endif

        if (pyid.LT.0.) then
        phimiss=-acos(pxid/sqrt(pxid*pxid+pyid*pyid))*180/(3.14159)

        endif

        if(abs(phimiss-pht).GT.180.) then
        dphimiss= 360-abs(phimiss-pht)
        endif

        if(abs(phimiss-pht).GT.180.) then
        dphimiss= abs(phimiss-pht)
        endif

C        write(*,*) dphimiss

        if (dphimiss.LT.20.) then  
          write(8,10)LQGPAR1(4),pxid,pyid,pxt,
     >    pyt,pht,ppt  
        endif

10      FORMAT(I8,6(1PE14.6))
           endif
60      enddo

C...List first few events.
       LQGPAR1(4)=LQGPAR1(4)+1
       if(LQGPAR1(4).LE.LQGPAR1(2))         CALL PYLIST(2)
C.. Swadhin
       if(mod(LQGPAR1(4),1000).eq.0)then
        write(6,1000) LQGPAR1(4)
1000    format('>>>>>> ',I8,
     >         ' events succesfully generated    <<<<<<')
       endif

C Write to output file
C Josh

C      WRITE(29,*) 0, LQGPAR1(4), 1, MSTI(1), MSTI(12), MSTI(16), 
C     &           PARI(34), MSTI(15), PARI(33), PARI(53)
       WRITE(29,*) '   0          1    1   95 2212         21  0.000000
     & 21 0.000000     0.000000      0.03587194790        
     & 0.00016131040 0.00000022484      718.31915283203      
     & 382.32009887695
     &   0.00773132309        0.00000000000           0.000000000
     &   0.000000000           0.000000000           0.000000000
     &   0.000000000           0.000000000           0.000000000
     &   0.000000000           0.000000000           0.000000000
     &   1.422321899             670'
      WRITE(29,*)'============================================'
       DO 100 L=1,N
C         WRITE(29,40) L, K(L,1), K(L,2), K(L,3)
C40       FORMAT(I5.2, I5.2, I7.2, I5.2)
          if(L == 3) then
          WRITE(29,*)'============================================'
          endif
C          if((K(L,1) .NE. 21) .AND. (KSprev .EQ. 21))
C     &           WRITE(29,*)'============================================'
C          WRITE(29,40) L, K(L,1), K(L,2), K(L,3), P(L,1), P(L,2)
C     &                 , P(L,3), P(L,4), P(L,5)
C40        FORMAT(I5.2, I5.2, I7.2, I5.2, F10.5, F10.5, F10.5
C     &           , F10.5, F10.5)
          WRITE(29,40) L, K(L,1), K(L,2), K(L,3), K(L,4), K(L,5),
     &                 P(L,1), P(L,2)
     &                 , P(L,3), P(L,4), P(L,5), V(L,1), V(L,2), V(L,3)
40        FORMAT(I5, I5, I7, I5, I5, I5, F10.5, F10.5, F10.5
     &           , F10.5, F10.5, F15.5, F15.5, F15.5)
          KSprev = K(L,1)
100    continue
       WRITE(29,*)'=============== Event finished ==============='
       CALL PYHIST

       if(LQGPAR1(3).ne.0)then
C...Fill histos
        CALL HF1(1000,sngl(x),1.) 
        if(genproc.eq.1)CALL HF1(1001,sngl(x),1.) 
        if(genproc.eq.2)CALL HF1(1002,sngl(x),1.) 
        CALL HF1(2000,sngl(y),1.) 
        if(genproc.eq.1)CALL HF1(2001,sngl(y),1.) 
        if(genproc.eq.2)CALL HF1(2002,sngl(y),1.) 
        CALL HF1(3000,log10(sngl(Q2)),1.) 

C... final energy and momentum checks
        px_sum=0.
        py_sum=0.
        pz_sum=0.
        e_sum=0.
        cme=0.
        do 222 i=1,N
         if(K(I,1).le.10)then
          px_sum=px_sum+P(I,1)
          py_sum=py_sum+P(I,2)
          pz_sum=pz_sum+P(I,3)
          e_sum=e_sum+P(I,4)
         endif
222     enddo
        cme=sqrt(e_sum**2-px_sum**2-py_sum**2-pz_sum**2) 
        call hf1(5001,sngl(px_sum),1.)
        call hf1(5002,sngl(py_sum),1.)
        call hf1(5003,sngl(pz_sum),1.)
        call hf1(5004,sngl(e_sum),1.)
        call hf1(5000,sngl(cme),1.)    
       endif
      endif

C-----------------------------------------------------------------
C...FLAG=2 -> Third section: Termination  
      if(flag.eq.2)Then
        write(6,*)

C...Pythia final table.
        CALL PYSTAT(1)

        write(6,*)
C...LQGENEP final statistics.
        CALL LQGPRI2

C...Closing Histograms.
        if(LQGPAR1(3).ne.0)then
         Call HCDIR('//lqgenep',' ')
         CALL HROUT(0,ICYCLE,' ')
         CALL HREND('lqgenep')
        endif
      endif
      END
*
C*********************************************************************
 
      SUBROUTINE PYUPEV(ISUB,SIGEV)
C-------------------------------------------
C...Pythia routine for user external process
C------------------------------------------- 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP

C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar

C# LQGpdfC.inc #
      double precision pdfsf(20)
      common /LQGpdfC/ pdfsf 
      
C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc # 
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2

C# LQGout.inc #
      double precision DXSEC(3),pvalence
      integer q_o,q_s,genproc,genprtyp,sch,uch,gproc(8)
      common /LQGout/ DXSEC,pvalence,q_o,q_s,genproc,genprtyp,
     >sch,uch,gproc
C...
      CHARACTER CHAF*16
      COMMON /PYDAT4/CHAF(500,2)  
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYUPPR/NUP,KUP(20,7),NFUP,IFUP(10,2),PUP(20,5),Q2UP(0:10)
      SAVE /PYDAT1/,/PYUPPR/
C...Local arrays and parameters.
      DIMENSION XPPs(-25:25),XPPu(-25:25),XPE(-25:25),TERM(20)
      DATA PI/3.141592653589793D0/
* conversion from pb to mb
      DATA CONV/1.D-9/
c      DATA CONV/1.D0/

C...LQGENEP parameters
      double precision BEAMPAR,LQGPAR3
      integer LQGPAR1,LQGPAR2
      COMMON/LQGDATA/BEAMPAR(3),LQGPAR1(10),LQGPAR2(10),LQGPAR3(20)

C...LQ names according to Aachen convention
      Character*7 LQCHA(14)
      DATA LQCHA /'S_0L','S_0R','~S_0R','S_1L',
     >            'V_1/2L','V_1/2R','~V_1/2L',
     >            'V_0L','V_0R','~V_0R','V_1L',
     >            'S_1/2L','S_1/2R','~S_1/2L'/                  
C...
      DATA KLQ /39/
*      
      sigev=0.d0
      irej=0

* sigma
      X=pyr(0)*(Xmax-Xmin)+Xmin
      Y=pyr(0)*(Ymax-Ymin)+Ymin
      Z=1
      Srad=S*z
      Q2=S*X*Y*Z
*      
C... Evaluate double differential cross section
      call LQGDDXS
*
      dxdydz=(Xmax-Xmin)*(Ymax-Ymin)*(Zmax-Zmin)
      Sigev=(DXSEC(2)+DXSEC(3))*conv*dxdydz
 
* fill Pythia variables for the generated process
* e beam         
      ECM2XZ=S*X*Z
      ECMXZ=sqrt(ECM2XZ)
      NUP=5
      KUP(1,1)=1
      KUP(1,2)=(echar)*-11
      KUP(1,3)=0
      KUP(1,4)=0
      KUP(1,5)=0
      KUP(1,6)=0
      KUP(1,7)=0
      PUP(1,1)=0.
      PUP(1,2)=0.
      PUP(1,4)=Z*sqrt(S)/2.d0
      PUP(1,3)=PUP(1,4)            
      PUP(1,5)=0.            
* p beam
      KUP(2,1)=1
      KUP(2,2)=q_s
      KUP(2,3)=0
      KUP(2,4)=0
      KUP(2,5)=0
      KUP(2,6)=0
      KUP(2,7)=0
      if(q_s.gt.0)then
       KUP(2,6)=5
      else
       KUP(2,7)=5
      endif
      PUP(2,1)=0.
      PUP(2,2)=0.
      PUP(2,4)=X*sqrt(S)/2.d0
      PUP(2,3)=-PUP(2,4)            
      PUP(2,5)=0.            
* LQ
      KUP(3,1)=2
      KUP(3,2)=KLQ      
      CHAF(pycomp(KLQ),1)=LQCHA(LQTYPE)      
      KUP(3,3)=0
      KUP(3,4)=0
      KUP(3,5)=0
      KUP(3,6)=0
      KUP(3,7)=0
      PUP(3,1)=PUP(2,1)+PUP(1,1)
      PUP(3,2)=PUP(2,2)+PUP(1,2)
      PUP(3,3)=PUP(2,3)+PUP(1,3)            
      PUP(3,4)=sqrt(ECM2XZ+PUP(3,1)**2+PUP(3,2)**2+PUP(3,3)**2)            
      PUP(3,5)=sqrt(ECM2XZ)            

* final state in sub-system cm.
*   final state lepton
      theta=acos(1.d0-2.d0*Y)
      PHI=2D0*PI*PYR(0)
      rtshat=ECMXZ
      KUP(4,1)=1
      KUP(4,2)=echar*-(11+2*(l_o-1))
      KUP(5,1)=1
      KUP(5,2)=q_o
      PUP(4,5)=PYMASS(KUP(4,2))
      PUP(5,5)=PYMASS(KUP(5,2))
      PUP44=0.5D0*(RTSHAT**2+PUP(4,5)**2-PUP(5,5)**2)/RTSHAT      
      PUP54=RTSHAT-PUP44
      KUP(4,3)=3
      KUP(4,4)=0
      KUP(4,5)=0
      KUP(4,6)=0
      KUP(4,7)=0
      if(irej.eq.1.and.PUP44**2-PUP(4,5)**2.lt.0)then
       PMOD=1.d0
      else
       PMOD=sqrt(PUP44**2-PUP(4,5)**2)      
      endif
      PUP(4,1)=PMOD*sin(theta)*cos(phi)
      PUP(4,2)=PMOD*sin(theta)*sin(phi)
      PUP43=PMOD*cos(theta)      
      PUP44=PUP(3,5)/2.d0
*   final state quark
      KUP(5,3)=3
      KUP(5,4)=0
      KUP(5,5)=0
      KUP(5,6)=0
      KUP(5,7)=0
      if(q_o.gt.0)then
       KUP(5,4)=2
      else
       KUP(5,5)=2
      endif
      PUP(5,1)=0.
      PUP(5,2)=0.
      if(irej.eq.1.and.PUP54**2-PUP(5,5)**2.lt.0)then
       PMOD=1.d0
      else
       PMOD=sqrt(PUP54**2-PUP(5,5)**2)      
      endif
      PUP(5,1)=-PUP(4,1)
      PUP(5,2)=-PUP(4,2)
      PUP53=-PUP43      
*   Longitudinal boost to cm frame
      beta=(z-x)/(z+x)
      gamma=0.5d0*(z+x)/sqrt(x*z)
      PUP(4,3)=GAMMA*(PUP43+BETA*PUP44)
      PUP(4,4)=GAMMA*(PUP44+BETA*PUP43)
      PUP(5,3)=GAMMA*(PUP53+BETA*PUP54)
      PUP(5,4)=GAMMA*(PUP54+BETA*PUP53)
* 
      NFUP=1
      IFUP(1,1)=4
      IFUP(1,2)=5
      Q2UP(0)=Q2
      Q2UP(1)=Q2
      RETURN
      END
*

      SUBROUTINE LQGDDXS
C----------------------------------------------
C...Evaluate double differential cross section 
C...           d^2 sigma / dx dy
C----------------------------------------------
*
      implicit none
*
C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar
            
C# LQGpdfC.inc #
      double precision pdfsf(20)
      common /LQGpdfC/ pdfsf
            
C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min
            
C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc #
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2
            
C# LQGout.inc #
      double precision DXSEC(3),pvalence
      integer q_o,q_s,genproc,genprtyp,sch,uch,gproc(8)
      common /LQGout/ DXSEC,pvalence,q_o,q_s,genproc,genprtyp,
     >sch,uch,gproc

      double precision DSIGMADXDY(4)
      double precision rand(1),sfrac1,sfrac2,ufrac1,ufrac2
      double precision pvalences_u,pvalences_d,pvalenceu_u,pvalenceu_d
*
cc--------------------------------------------------------
C...Leptoquark types ranges from 1 to 14.
C 
C     1->S_0 LEFT
C     2->S_0 RIGHT
C     3->~S_0 RIGHT
C     4->S_1 LEFT
C     5->V_1/2 LEFT
C     6->V_1/2 RIGHT
C     7->~V_1/2 LEFT
C     8->V_0 LEFT
C     9->V_0 RIGHT
C     10->~V_0 RIGHT
C     11->V_1 LEFT
C     12->S_1/2 LEFT
C     13->S_1/2 RIGHT   
C     14->~S_1/2 LEFT
C
C       DSIGMADXDY(4) - Double differential cross section
C                       DSIGMADXDY(1) = Standard Model term (SM processes)
C                       DSIGMADXDY(2) = Interference term between SM and LQ
C                       DSIGMADXDY(3) = LQ term - u channel
C                       DSIGMADXDY(4) = LQ term - s channel
C ========================================================================
C INPUT PARAMETERS:
C                   X - standard DIS x variable 
C                   Y - standard DIS y variable
C                   S - Center of mass energy
C                 MLQ - Leptoquark mass
C                  G1 - initial state coupling
C                  G2 - final state coupling
C                 l_o - generation of the outcoming lepton
C               echar - charge of the incoming lepton
C                 q_i - generation of initial state quark
C                 q_j - generation of the final state quark
C              LQTYPE - Leptoquark type (see table above)
C
C OUTPUT PARAMETERS:
C                       DXSEC = double differential cross section (pb):
C                            DXSEC(1)= LQ-SM interference term
C                            DXSEC(2)= LQ term - u channel 
C                            DXSEC(3)= LQ term - s channel
C                       q_o = output quark (from LQ decay):
C                            1 down        -1 antidown 
C                            2 up          -2 antiup
C                            3 strange     -3 antistrange 
C                            4 charm       -4 anticharm
C                            5 bottom      -5 antibottom
C                            6 top         -6 antitop
C            
C----------------------------------------------------------
      double precision pyr
      double precision C_R_P,C_L_P,C_R_E,C_L_E
     &                ,C_R_U,C_L_U,C_R_D,C_L_D
     &                ,B_RR_U,B_RL_U
     &                ,B_LR_U,B_LL_U
     &                ,B_RR_D,B_RL_D
     &                ,B_LR_D,B_LL_D
      double precision CCCf2u,DDDf2u,EEEf2u,FFFf2u,
     &                 CCCf2d,DDDf2d,EEEf2d,FFFf2d
      double precision CCCf0u,DDDf0u,EEEf0u,FFFf0u,
     &                 CCCf0d,DDDf0d,EEEf0d,FFFf0d
      double precision CCCv2u,DDDv2u,EEEv2u,FFFv2u,
     &                 CCCv2d,DDDv2d,EEEv2d,FFFv2d
      double precision CCCv0u,DDDv0u,EEEv0u,FFFv0u,
     &                 CCCv0d,DDDv0d,EEEv0d,FFFv0d

      double precision weakmix,Mz2,Gz2,Gf,alpha,A,P,pi
      parameter (weakmix=0.2315, Mz2=(91.187)**2,Gz2=(2.490)**2)
      parameter (Gf=0.0000116639,alpha=1.0/137.036, A=27.5, P=820.0)
      parameter (pi=3.141592653589793d0)

      double precision UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL
      double precision UPQs(3),DNQs(3),UPQBs(3),DNQBs(3)
      double precision UPQu(3),DNQu(3),UPQBu(3),DNQBu(3)
      double precision scales,scaleu,LambdaL2,LambdaR2,Ms2
     &                ,aaa,bbb,ggg
     &                ,LambdaL2_1,LambdaL2_2
     &                ,LambdaR2_1,LambdaR2_2
     &                ,GAM

c......LH 0, RH 1
      integer LRindex(14)
      data LRindex/0,1,1,0,0,1,0,0,1,1,0,0,1,0/
      integer SVindex(14)
      data SVindex/1,1,1,1,2,2,2,3,3,3,3,4,4,4/

* protection
      if(y.eq.1)y=1.d0-1.d-13
*
      DSIGMADXDY(1)=0.d0
      DSIGMADXDY(2)=0.d0
      DSIGMADXDY(3)=0.d0
      DSIGMADXDY(4)=0.d0
      q_o=0
      pvalence=0.
      pvalences_u=0.
      pvalenceu_u=0.
      pvalences_d=0.
      pvalenceu_d=0.
*
      C_R_P = weakmix
      C_L_P = -0.5+weakmix
      C_R_U = -2.0*weakmix/3.0
      C_L_U = 0.5-2.0*weakmix/3.0
      C_R_D = weakmix/3.0
      C_L_D = -0.5+weakmix/3.0
      C_R_E = weakmix
      C_L_E = -0.5+weakmix

      Ms2=(MLQ)**2
      Q2=Srad*X*Y
      if(Q2.lt.Q2min)goto 999
* u channel densities
      scaleu=sqrt(Srad*X*(1.-Y))    
      CALL STRUCTM(X,scaleu,UPV
     &               ,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)                   
      if(UPV+USEA.gt.0)
     & pvalenceu_u=UPV/(UPV+USEA)
      if(DNV+DSEA.gt.0)
     & pvalenceu_d=DNV/(DNV+DSEA)
      if(echar.eq.1)then
* case e+, mu+, tau+ 
       UPQu(1)=UPV+USEA
       UPQBu(1)=USEA
       UPQu(2)=CHM
       UPQBu(2)=CHM
       UPQu(3)=TOP
       UPQBu(3)=TOP
       DNQu(1)=DNV+DSEA
       DNQBu(1)=DSEA
       DNQu(2)=STR
       DNQBu(2)=STR
       DNQu(3)=BOT
       DNQBu(3)=BOT
      elseif(echar.eq.-1)then
* case e-, mu-, tau- 
       UPQu(1)=USEA
       UPQBu(1)=UPV+USEA
       UPQu(2)=CHM
       UPQBu(2)=CHM
       UPQu(3)=TOP
       UPQBu(3)=TOP
       DNQu(1)=DSEA
       DNQBu(1)=DNV+DSEA
       DNQu(2)=STR
       DNQBu(2)=STR
       DNQu(3)=BOT
       DNQBu(3)=BOT
      endif          
* s channel densities
      scales=sqrt(Srad*X)    
      CALL STRUCTM(X,scales,UPV
     &               ,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GL)        
      if(UPV+USEA.gt.0)
     & pvalences_u=UPV/(UPV+USEA)
      if(DNV+DSEA.gt.0)
     & pvalences_d=DNV/(DNV+DSEA)           
      if(echar.eq.1)then
* case e+, mu+, tau+ 
       UPQs(1)=UPV+USEA
       UPQBs(1)=USEA
       UPQs(2)=CHM
       UPQBs(2)=CHM
       UPQs(3)=TOP
       UPQBs(3)=TOP
       DNQs(1)=DNV+DSEA
       DNQBs(1)=DSEA
       DNQs(2)=STR
       DNQBs(2)=STR
       DNQs(3)=BOT
       DNQBs(3)=BOT
      elseif(echar.eq.-1)then
* case e-, mu-, tau- 
       UPQs(1)=USEA
       UPQBs(1)=UPV+USEA
       UPQs(2)=CHM
       UPQBs(2)=CHM
       UPQs(3)=TOP
       UPQBs(3)=TOP
       DNQs(1)=DSEA
       DNQBs(1)=DNV+DSEA
       DNQs(2)=STR
       DNQBs(2)=STR
       DNQs(3)=BOT
       DNQBs(3)=BOT
      endif
*         
      aaa = Q2*(Q2+Mz2)/((Q2+Mz2)**2+Mz2*Gz2)
      bbb = sqrt(2.0)*Gf*Mz2/(pi*alpha)
    
      B_RR_U = -2./3. + aaa*bbb*C_R_P*C_R_U
      B_RL_U = -2./3. + aaa*bbb*C_R_P*C_L_U
      B_LR_U = -2./3. + aaa*bbb*C_L_P*C_R_U
      B_LL_U = -2./3. + aaa*bbb*C_L_P*C_L_U
      B_RR_D = 1./3. + aaa*bbb*C_R_P*C_R_D
      B_RL_D = 1./3. + aaa*bbb*C_R_P*C_L_D
      B_LR_D = 1./3. + aaa*bbb*C_L_P*C_R_D
      B_LL_D = 1./3. + aaa*bbb*C_L_P*C_L_D

      IF (SVindex(LQTYPE).EQ.1) THEN

       CCCf2u=Q2*(1.-Y)**2*(UPV+USEA)
     &                /(4.0*pi*alpha)
       DDDf2u=Q2*USEA/(4.0*pi*alpha)
       EEEf2u=(Q2-X*Srad)**2*Y**2*(UPQu(q_j))
     &            /(64.0*(pi*alpha)**2)
       FFFf2u=Q2**2*UPQBs(q_i)/(64.0*(pi*alpha)**2)

       CCCf2d=Q2*(1.-Y)**2*(DNV+DSEA)
     &                /(4.0*pi*alpha)
       DDDf2d=Q2*DSEA/(4.0*pi*alpha)
       EEEf2d=(Q2-X*Srad)**2*Y**2*(DNQu(q_j))
     &            /(64.0*(pi*alpha)**2)
       FFFf2d=Q2**2*DNQBs(q_i)/(64.0*(pi*alpha)**2)

      ELSE IF (SVindex(LQTYPE).EQ.4) THEN

       CCCf0u=Q2*(1.-Y)**2*USEA/(4.0*pi*alpha)
       DDDf0u=Q2*(UPV+USEA)/(4.0*pi*alpha)
       EEEf0u=(Q2-X*Srad)**2*Y**2*UPQBu(q_j)
     &            /(64.0*(pi*alpha)**2)
       FFFf0u=Q2**2*(UPQs(q_i))/(64.0*(pi*alpha)**2)

       CCCf0d=Q2*(1.-Y)**2*DSEA/(4.0*pi*alpha)
       DDDf0d=Q2*(DNV+DSEA)/(4.0*pi*alpha)
       EEEf0d=(Q2-X*Srad)**2*Y**2*DNQBu(q_j)
     &            /(64.0*(pi*alpha)**2)
       FFFf0d=Q2**2*(DNQs(q_i))/(64.0*(pi*alpha)**2)

      ELSE IF (SVindex(LQTYPE).EQ.2) THEN

       CCCv2u=Q2*(1.-Y)**2*USEA/(2.0*pi*alpha)
       DDDv2u=Q2*(UPV+USEA)/(2.0*pi*alpha)
       EEEv2u=Q2**2*UPQBs(q_i)/(16.0*(pi*alpha)**2)
       FFFv2u=(Q2-X*Srad)**2*Y**2*(UPQu(q_j))
     &            /(16.0*(pi*alpha)**2*(1.0-Y)**2)

       CCCv2d=Q2*(1.-Y)**2*DSEA/(2.0*pi*alpha)
       DDDv2d=Q2*(DNV+DSEA)/(2.0*pi*alpha)
       EEEv2d=Q2**2*DNQBs(q_i)/(16.0*(pi*alpha)**2)
       FFFv2d=(Q2-X*Srad)**2*Y**2*(DNQu(q_j))
     &            /(16.0*(pi*alpha)**2*(1.0-Y)**2)

      ELSE IF (SVindex(LQTYPE).EQ.3) THEN

       CCCv0u=Q2*(1.-Y)**2*(UPV+USEA)
     &            /(2.0*pi*alpha)
       DDDv0u=Q2*USEA/(2.0*pi*alpha)
       EEEv0u=Q2**2*(UPQs(q_i))/(16.0*(pi*alpha)**2)
       FFFv0u=(Q2-X*Srad)**2*Y**2*UPQBu(q_j)
     &            /(16.0*(pi*alpha)**2*(1.0-Y)**2)

       CCCv0d=Q2*(1.-Y)**2*(DNV+DSEA)
     &              /(2.0*pi*alpha)
       DDDv0d=Q2*DSEA/(2.0*pi*alpha)
       EEEv0d=Q2**2*(DNQs(q_i))/(16.0*(pi*alpha)**2)
       FFFv0d=(Q2-X*Srad)**2*Y**2*DNQBu(q_j)
     &            /(16.0*(pi*alpha)**2*(1.0-Y)**2)
      ENDIF

      DSIGMADXDY(1)=(B_RL_U**2+B_LR_U**2+
     &              (B_RR_U**2+B_LL_U**2)*(1.-Y)**2)
     &               *(UPV+USEA+CHM+TOP)+
     &              (B_RL_D**2+B_LR_D**2+
     &              (B_RR_D**2+B_LL_D**2)*(1.-Y)**2)
     &               *(DNV+DSEA+STR+BOT)+
     &              (B_RR_U**2+B_LL_U**2+
     &              (B_RL_U**2+B_LR_U**2)*(1.-Y)**2)
     &               *(USEA+CHM+TOP) +
     &              (B_RR_D**2+B_LL_D**2+
     &              (B_RL_D**2+B_LR_D**2)*(1.-Y)**2)
     &               *(DSEA+STR+BOT)

      IF (LRindex(LQTYPE).EQ.0) THEN
         LambdaL2_1=G1
         LambdaL2_2=G2
         if(G2.ne.0)then
          LambdaL2=LambdaL2_1*LambdaL2_2
         else
          LambdaL2=G1*G1
         endif
         LambdaR2_1=0.0
         LambdaR2_2=0.0
         LambdaR2=LambdaR2_1*LambdaR2_2
      ELSE IF (LRindex(LQTYPE).EQ.1) THEN
         LambdaR2_1=G1
         LambdaR2_2=G2
         if(G2.ne.0)then
          LambdaR2=LambdaR2_1*LambdaR2_2
         else
          LambdaR2=G1*G1
         endif
         LambdaL2_1=0.0
         LambdaL2_2=0.0
         LambdaL2=LambdaL2_1*LambdaL2_2
      ENDIF

      IF (LQTYPE.EQ.1.or.LQTYPE.EQ.2) THEN
         GAM=MLQ*(sqrt(2.0)*LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &           MLQ*(sqrt(2.0)*LambdaL2_2+LambdaR2_2)**2
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=LambdaR2*B_RR_U+LambdaL2*B_LL_U
         DSIGMADXDY(2)=BBB*CCCf2u/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf2u*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
         if(q_i.ne.3)then
          DSIGMADXDY(3)=EEEf2u*GGG/(Q2-X*Srad-Ms2)**2
         else
* Top quark in the final state
          DSIGMADXDY(3)=0.d0
         endif         
         if(q_j.ne.3)then
          DSIGMADXDY(4)=FFFf2u*GGG/AAA
         else
* Top quark in the final state
          DSIGMADXDY(4)=0.d0
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(4)+DSIGMADXDY(3))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ ub -> l+ ub
          genproc=1
          genprtyp=3
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=-2
          if(q_j.eq.2)q_o=-4
          if(q_j.eq.3)q_o=-6
          if(q_i.eq.1)q_s=-2
          if(q_i.eq.2)q_s=-4
          if(q_i.eq.3)q_s=-6          
         else
* u channel: e+ u -> l+ u
          genproc=2
          genprtyp=5
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=2
          if(q_i.eq.2)q_o=4
          if(q_i.eq.3)q_o=6
          if(q_j.eq.1)q_s=2
          if(q_j.eq.2)q_s=4
          if(q_j.eq.3)q_s=6          
         endif
      ELSE IF (LQTYPE.EQ.3) THEN
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=LambdaR2*B_RR_D+LambdaL2*B_LL_D
         DSIGMADXDY(2)=BBB*CCCf2d/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf2d*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
         DSIGMADXDY(3)=EEEf2d*GGG/(Q2-X*Srad-Ms2)**2
         DSIGMADXDY(4)=FFFf2d*GGG/AAA
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(4)+DSIGMADXDY(3))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ db -> l+ db
          genproc=1
          genprtyp=4
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=-1
          if(q_j.eq.2)q_o=-3
          if(q_j.eq.3)q_o=-5
          if(q_i.eq.1)q_s=-1
          if(q_i.eq.2)q_s=-3
          if(q_i.eq.3)q_s=-5          
         else
* u channel: e+ d -> l+ d
          genproc=2
          genprtyp=6
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=1
          if(q_i.eq.2)q_o=3
          if(q_i.eq.3)q_o=5
          if(q_j.eq.1)q_s=1
          if(q_j.eq.2)q_s=3
          if(q_j.eq.3)q_s=5
         endif
      ELSE IF (LQTYPE.EQ.4) THEN
         GAM=MLQ*(sqrt(2.0)*LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(sqrt(2.0)*LambdaL2_2+LambdaR2_2)**2
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=LambdaR2*B_RR_D+2.0*LambdaL2*B_LL_D
         DSIGMADXDY(2)=BBB*CCCf2d/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf2d*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+2.0*LambdaL2)**2
         DSIGMADXDY(3)=EEEf2d*GGG/(Q2-X*Srad-Ms2)**2
         DSIGMADXDY(4)=FFFf2d*GGG/AAA
         sfrac1=DSIGMADXDY(4)
         ufrac1=DSIGMADXDY(3)
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=LambdaR2*B_RR_U+LambdaL2*B_LL_U
         DSIGMADXDY(2)=DSIGMADXDY(2)+BBB*CCCf2u/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf2u*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
* no Top quark in the final state, otherwise no u-channel contribution
         if(q_i.ne.3)then
          DSIGMADXDY(3)=DSIGMADXDY(3)+EEEf2u*GGG/(Q2-X*Srad-Ms2)**2
         endif
* no Top quark in the final state, otherwise no s-channel contribution
         if(q_j.ne.3)then
          DSIGMADXDY(4)=DSIGMADXDY(4)+FFFf2u*GGG/AAA
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac2=DSIGMADXDY(4)-sfrac1
         ufrac2=DSIGMADXDY(3)-ufrac1
         sfrac1=sfrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         sfrac2=sfrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac1=ufrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac2=ufrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ db -> l+ db
          genproc=1
          genprtyp=4
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=-1
          if(q_j.eq.2)q_o=-3
          if(q_j.eq.3)q_o=-5
          if(q_i.eq.1)q_s=-1
          if(q_i.eq.2)q_s=-3
          if(q_i.eq.3)q_s=-5          
         elseif(rand(1).lt.(sfrac1+sfrac2))then
* s channel: e+ ub -> l+ ub
          genproc=1
          genprtyp=3
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=-2
          if(q_j.eq.2)q_o=-4
          if(q_j.eq.3)q_o=-6
          if(q_i.eq.1)q_s=-2
          if(q_i.eq.2)q_s=-4
          if(q_i.eq.3)q_s=-6
         elseif(rand(1).lt.(sfrac1+sfrac2+ufrac1))then
* u channel: e+ d -> l+ d
          genproc=2
          genprtyp=6
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=1
          if(q_i.eq.2)q_o=3
          if(q_i.eq.3)q_o=5
          if(q_j.eq.1)q_s=1
          if(q_j.eq.2)q_s=3
          if(q_j.eq.3)q_s=5          
         else
* u channel: e+ u -> l+ u
          genproc=2
          genprtyp=5
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=2
          if(q_i.eq.2)q_o=4
          if(q_i.eq.3)q_o=6
          if(q_j.eq.1)q_s=2
          if(q_j.eq.2)q_s=4
          if(q_j.eq.3)q_s=6
         endif
      ELSE IF (LQTYPE.EQ.5) THEN
         AAA=(X*Srad-Ms2)
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         BBB=LambdaR2*B_RL_D+LambdaL2*B_LR_D
         GGG=(LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=BBB*CCCv2d*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv2d/(Q2-X*Srad-Ms2)
         DSIGMADXDY(4)=EEEv2d*
     &         ((LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         DSIGMADXDY(3)=FFFv2d*(LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ db -> l+ db
          genproc=1
          genprtyp=4
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=-1
          if(q_j.eq.2)q_o=-3
          if(q_j.eq.3)q_o=-5
          if(q_i.eq.1)q_s=-1
          if(q_i.eq.2)q_s=-3
          if(q_i.eq.3)q_s=-5
         else
* u channel: e+ d -> l+ d
          genproc=2
          genprtyp=6
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=1
          if(q_i.eq.2)q_o=3
          if(q_i.eq.3)q_o=5
          if(q_j.eq.1)q_s=1
          if(q_j.eq.2)q_s=3
          if(q_j.eq.3)q_s=5          
         endif
      ELSE IF (LQTYPE.EQ.6) THEN
         AAA=(X*Srad-Ms2)
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         BBB=LambdaR2*B_RL_D+LambdaL2*B_LR_D
         GGG=(LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=BBB*CCCv2d*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv2d/(Q2-X*Srad-Ms2)
         DSIGMADXDY(4)=EEEv2d*
     &         ((LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         DSIGMADXDY(3)=FFFv2d*(LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
* output quark
         sfrac1=DSIGMADXDY(4)
         ufrac1=DSIGMADXDY(3)
         BBB=LambdaR2*B_RL_U+LambdaL2*B_LR_U
         GGG=(LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=DSIGMADXDY(2)+BBB*CCCv2u*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv2u/(Q2-X*Srad-Ms2)
* no Top quark in the final state, otherwise no s-channel contribution
         if(q_j.ne.3)then
          DSIGMADXDY(4)=DSIGMADXDY(4)+EEEv2u*
     &         ((LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         endif
* no Top quark in the final state, otherwise no u-channel contribution
         if(q_i.ne.3)then
          DSIGMADXDY(3)=DSIGMADXDY(3)+FFFv2u*(LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac2=DSIGMADXDY(4)-sfrac1
         ufrac2=DSIGMADXDY(3)-ufrac1
         sfrac1=sfrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         sfrac2=sfrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac1=ufrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac2=ufrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ db -> l+ db
          genproc=1
          genprtyp=4
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=-1
          if(q_j.eq.2)q_o=-3
          if(q_j.eq.3)q_o=-5
          if(q_i.eq.1)q_s=-1
          if(q_i.eq.2)q_s=-3
          if(q_i.eq.3)q_s=-5          
         elseif(rand(1).lt.(sfrac1+sfrac2))then
* s channel: e+ ub -> l+ ub
          genproc=1
          genprtyp=3
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=-2
          if(q_j.eq.2)q_o=-4
          if(q_j.eq.3)q_o=-6
          if(q_i.eq.1)q_s=-2
          if(q_i.eq.2)q_s=-4
          if(q_i.eq.3)q_s=-6          
         elseif(rand(1).lt.(sfrac1+sfrac2+ufrac1))then
* u channel: e+ d -> l+ d
          genproc=2
          genprtyp=6
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=1
          if(q_i.eq.2)q_o=3
          if(q_i.eq.3)q_o=5
          if(q_j.eq.1)q_s=1
          if(q_j.eq.2)q_s=3
          if(q_j.eq.3)q_s=5          
         else
* u channel: e+ u -> l+ u
          genproc=2
          genprtyp=5
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=2
          if(q_i.eq.2)q_o=4
          if(q_i.eq.3)q_o=6
          if(q_j.eq.1)q_s=2
          if(q_j.eq.2)q_s=4
          if(q_j.eq.3)q_s=6
         endif
      ELSE IF (LQTYPE.EQ.7) THEN
         AAA=(X*Srad-Ms2)
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         BBB=LambdaR2*B_RL_U+LambdaL2*B_LR_U
         GGG=(LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=BBB*CCCv2u*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv2u/(Q2-X*Srad-Ms2)
* no Top quark in the final state, otherwise no s-channel contribution
         if(q_j.ne.3)then
          DSIGMADXDY(4)=EEEv2u*
     &         ((LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         else    
          DSIGMADXDY(4)=0.d0
         endif
* no Top quark in the final state, otherwise no s-channel contribution
         if(q_i.ne.3)then
          DSIGMADXDY(3)=FFFv2u*(LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
         else    
          DSIGMADXDY(3)=0.d0
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ ub -> l+ ub
          genproc=1
          genprtyp=3
          if(q_i.eq.1)then
           if(echar.lt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=-2
          if(q_j.eq.2)q_o=-4
          if(q_j.eq.3)q_o=-6
          if(q_i.eq.1)q_s=-2
          if(q_i.eq.2)q_s=-4
          if(q_i.eq.3)q_s=-6          
         else
* u channel: e+ u -> l+ u
          genproc=2
          genprtyp=5
          if(q_j.eq.1)then
           if(echar.gt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=2
          if(q_i.eq.2)q_o=4
          if(q_i.eq.3)q_o=6
          if(q_j.eq.1)q_s=2
          if(q_j.eq.2)q_s=4
          if(q_j.eq.3)q_s=6          
         endif
      ELSE IF (LQTYPE.EQ.8.or.LQTYPE.EQ.9) THEN
         AAA=(X*Srad-Ms2)
         GAM=MLQ*(sqrt(2.0)*LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(sqrt(2.0)*LambdaL2_2+LambdaR2_2)**2
         BBB=-LambdaR2*B_RR_D-LambdaL2*B_LL_D
         GGG=(2.0*LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=BBB*CCCv0d*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv0d/(Q2-X*Srad-Ms2)
         DSIGMADXDY(4)=EEEv0d*
     &         ((LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         DSIGMADXDY(3)=FFFv0d*(LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ d -> l+ d
          genproc=1
          genprtyp=2
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=1
          if(q_j.eq.2)q_o=3
          if(q_j.eq.3)q_o=5
          if(q_i.eq.1)q_s=1
          if(q_i.eq.2)q_s=3
          if(q_i.eq.3)q_s=5
         else
* u channel: e+ db -> l+ db
          genproc=2
          genprtyp=8
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=-1
          if(q_i.eq.2)q_o=-3
          if(q_i.eq.3)q_o=-5
          if(q_j.eq.1)q_s=-1
          if(q_j.eq.2)q_s=-3
          if(q_j.eq.3)q_s=-5
         endif
      ELSE IF (LQTYPE.EQ.10) THEN
         AAA=(X*Srad-Ms2)
         GAM=MLQ*(sqrt(2.0)*LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(sqrt(2.0)*LambdaL2_2+LambdaR2_2)**2
         BBB=-LambdaR2*B_RR_U-2.0*LambdaL2*B_LL_U
         GGG=(2.0*LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=BBB*CCCv0u*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv0u/(Q2-X*Srad-Ms2)
         if(q_j.ne.3)then
          DSIGMADXDY(4)=EEEv0u*
     &         ((4.0*LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         else
* Top quark in the final state: no s-channel contribution
          DSIGMADXDY(4)=0.d0
         endif
         if(q_i.ne.3)then
          DSIGMADXDY(3)=
     &          FFFv0u*(4.0*LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*2.0*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
         else
* Top quark in the final state: no u-channel contribution
          DSIGMADXDY(3)=0.d0
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ u -> l+ u
          genproc=1
          genprtyp=1
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=2
          if(q_j.eq.2)q_o=4
          if(q_j.eq.3)q_o=6
          if(q_i.eq.1)q_s=2
          if(q_i.eq.2)q_s=4
          if(q_i.eq.3)q_s=6
         else
* u channel: e+ ub -> l+ ub
          genproc=2
          genprtyp=7
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=-2
          if(q_i.eq.2)q_o=-4
          if(q_i.eq.3)q_o=-6
          if(q_j.eq.1)q_s=-2
          if(q_j.eq.2)q_s=-4
          if(q_j.eq.3)q_s=-6
         endif
      ELSE IF (LQTYPE.EQ.11) THEN
         AAA=(X*Srad-Ms2)
         GAM=MLQ*(sqrt(2.0)*LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(sqrt(2.0)*LambdaL2_2+LambdaR2_2)**2
         BBB=-LambdaR2*B_RR_D-LambdaL2*B_LL_D
         GGG=(2.0*LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=BBB*CCCv0d*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv0d/(Q2-X*Srad-Ms2)
         DSIGMADXDY(4)=EEEv0d*
     &         ((LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         DSIGMADXDY(3)=FFFv0d*(LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
* output quark
         sfrac1=DSIGMADXDY(4)
         ufrac1=DSIGMADXDY(3)
         BBB=-LambdaR2*B_RR_U-2.0*LambdaL2*B_LL_U
         GGG=(2.0*LambdaL2+LambdaR2)/(24.0*pi)
         DSIGMADXDY(2)=DSIGMADXDY(2)+BBB*CCCv0u*AAA/(AAA**2+
     &          (Ms2*GGG)**2)+
     &          BBB*DDDv0u/(Q2-X*Srad-Ms2)
         if(q_j.ne.3)then
* no Top quark in the final state otherwise no s-channel contribution
          DSIGMADXDY(4)=DSIGMADXDY(4)+EEEv0u*
     &         ((4.0*LambdaL2**2+LambdaR2**2)*(1.-Y)**2+
     &         2.0*2.0*LambdaL2*LambdaR2*Y**2)/
     &         (AAA**2+Ms2*GAM**2/(24.0*pi)**2)
         endif
         if(q_i.ne.3)then
* no Top quark in the final state otherwise no u-channel contribution
          DSIGMADXDY(3)=DSIGMADXDY(3)+
     &          FFFv0u*(4.0*LambdaL2**2+LambdaR2**2+
     &          2.0*Y**2*2.0*LambdaL2*LambdaR2)/
     &          (Q2-X*Srad-Ms2)**2
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac2=DSIGMADXDY(4)-sfrac1
         ufrac2=DSIGMADXDY(3)-ufrac1
         sfrac1=sfrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         sfrac2=sfrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac1=ufrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac2=ufrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ d -> l+ d
          genproc=1
          genprtyp=2
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=1
          if(q_j.eq.2)q_o=3
          if(q_j.eq.3)q_o=5
          if(q_i.eq.1)q_s=1
          if(q_i.eq.2)q_s=3
          if(q_i.eq.3)q_s=5
         elseif(rand(1).lt.(sfrac1+sfrac2))then
* s channel: e+ u -> l+ u
          genproc=1
          genprtyp=1
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=2
          if(q_j.eq.2)q_o=4
          if(q_j.eq.3)q_o=6
          if(q_i.eq.1)q_s=2
          if(q_i.eq.2)q_s=4
          if(q_i.eq.3)q_s=6
         elseif(rand(1).lt.(sfrac1+sfrac2+ufrac1))then
* u channel: e+ db -> l+ db
          genproc=2
          genprtyp=8
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=-1
          if(q_i.eq.2)q_o=-3
          if(q_i.eq.3)q_o=-5
          if(q_j.eq.1)q_s=-1
          if(q_j.eq.2)q_s=-3
          if(q_j.eq.3)q_s=-5
         else
* u channel: e+ ub -> l+ ub
          genproc=2
          genprtyp=7
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=-2
          if(q_i.eq.2)q_o=-4
          if(q_i.eq.3)q_o=-6
          if(q_j.eq.1)q_s=-2
          if(q_j.eq.2)q_s=-4
          if(q_j.eq.3)q_s=-6
         endif
      ELSE IF (LQTYPE.EQ.12) THEN
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=-LambdaR2*B_RL_U-LambdaL2*B_LR_U
         DSIGMADXDY(2)=BBB*CCCf0u/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf0u*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
         if(q_i.ne.3)then
          DSIGMADXDY(3)=EEEf0u*GGG/(Q2-X*Srad-Ms2)**2
         else
* Top quark in the final state no u-channel contribution
          DSIGMADXDY(3)=0.d0
         endif
         if(q_j.ne.3)then
          DSIGMADXDY(4)=FFFf0u*GGG/AAA
         else
* Top quark in the final state no s-channel contribution
          DSIGMADXDY(4)=0.d0
         endif
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         if(DSIGMADXDY(3).gt.0)then
          sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(3)+DSIGMADXDY(4))
         else
          sfrac1=0.d0
         endif
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ u -> l+ u
          genproc=1
          genprtyp=1
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=2
          if(q_j.eq.2)q_o=4
          if(q_j.eq.3)q_o=6
          if(q_i.eq.1)q_s=2
          if(q_i.eq.2)q_s=4
          if(q_i.eq.3)q_s=6
         else
* u channel: e+ ub -> l+ ub
          genproc=2
          genprtyp=7
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=-2
          if(q_i.eq.2)q_o=-4
          if(q_i.eq.3)q_o=-6
          if(q_j.eq.1)q_s=-2
          if(q_j.eq.2)q_s=-4
          if(q_j.eq.3)q_s=-6
         endif
      ELSE IF (LQTYPE.EQ.13) THEN
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=-LambdaR2*B_RL_U-LambdaL2*B_LR_U
         DSIGMADXDY(2)=BBB*CCCf0u/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf0u*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
         if(q_i.ne.3)then
          DSIGMADXDY(3)=EEEf0u*GGG/(Q2-X*Srad-Ms2)**2
         else
* Top quark in the final state no u-channel contribution 
          DSIGMADXDY(3)=0.d0
         endif
         if(q_j.ne.3)then
          DSIGMADXDY(4)=FFFf0u*GGG/AAA
         else
* Top quark in the final state no s-channel contribution 
          DSIGMADXDY(3)=0.d0
         endif
* output quark
         sfrac1=DSIGMADXDY(4)
         ufrac1=DSIGMADXDY(3)
         BBB=-LambdaR2*B_RL_D-LambdaL2*B_LR_D
         DSIGMADXDY(2)=DSIGMADXDY(2)+BBB*CCCf0d/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf0d*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
         DSIGMADXDY(3)=DSIGMADXDY(3)+EEEf0d*GGG/(Q2-X*Srad-Ms2)**2
         DSIGMADXDY(4)=DSIGMADXDY(4)+FFFf0d*GGG/AAA
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac2=DSIGMADXDY(4)-sfrac1
         ufrac2=DSIGMADXDY(3)-ufrac1
         sfrac1=sfrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         sfrac2=sfrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac1=ufrac1/(DSIGMADXDY(3)+DSIGMADXDY(4))
         ufrac2=ufrac2/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ u -> l+ u
          genproc=1
          genprtyp=1
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_u
          endif
          if(q_j.eq.1)q_o=2
          if(q_j.eq.2)q_o=4
          if(q_j.eq.3)q_o=6
          if(q_i.eq.1)q_s=2
          if(q_i.eq.2)q_s=4
          if(q_i.eq.3)q_s=6
         elseif(rand(1).lt.(sfrac1+sfrac2))then
* s channel: e+ d -> l+ d
          genproc=1
          genprtyp=2
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=1
          if(q_j.eq.2)q_o=3
          if(q_j.eq.3)q_o=5
          if(q_i.eq.1)q_s=1
          if(q_i.eq.2)q_s=3
          if(q_i.eq.3)q_s=5
         elseif(rand(1).lt.(sfrac1+sfrac2+ufrac1))then
* u channel: e+ ub -> l+ ub
          genproc=2
          genprtyp=7
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_u
          endif
          if(q_i.eq.1)q_o=-2
          if(q_i.eq.2)q_o=-4
          if(q_i.eq.3)q_o=-5
          if(q_j.eq.1)q_s=-2
          if(q_j.eq.2)q_s=-4
          if(q_j.eq.3)q_s=-5
         else
* u channel: e+ db -> l+ db
          genproc=2
          genprtyp=8
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=-1
          if(q_i.eq.2)q_o=-3
          if(q_i.eq.3)q_o=-5
          if(q_j.eq.1)q_s=-1
          if(q_j.eq.2)q_s=-3
          if(q_j.eq.3)q_s=-5
         endif
      ELSE IF (LQTYPE.EQ.14) THEN
         GAM=MLQ*(LambdaL2_1+LambdaR2_1)**2
         if(l_o.gt.1)GAM=GAM+
     &            MLQ*(LambdaL2_2+LambdaR2_2)**2
         AAA=(X*Srad-Ms2)**2+
     &        Ms2*GAM**2/((16.0*pi)**2)
         BBB=-LambdaR2*B_RL_D-LambdaL2*B_LR_D
         DSIGMADXDY(2)=BBB*CCCf0d/(Q2-X*Srad-Ms2)+
     &             BBB*DDDf0d*(X*Srad-Ms2)/AAA
         GGG=(LambdaR2+LambdaL2)**2
         DSIGMADXDY(3)=EEEf0d*GGG/(Q2-X*Srad-Ms2)**2
         DSIGMADXDY(4)=FFFf0d*GGG/AAA
* output quark
         if(DSIGMADXDY(3)+DSIGMADXDY(4).eq.0)then
          goto 999
         endif          
         sfrac1=DSIGMADXDY(4)/(DSIGMADXDY(3)+DSIGMADXDY(4))
c         call rm48(rand,1)
         rand(1)=pyr(0)
         if(rand(1).lt.sfrac1)then
* s channel: e+ d -> l+ d
          genproc=1
          genprtyp=2
          if(q_i.eq.1)then
           if(echar.gt.0)pvalence=pvalences_d
          endif
          if(q_j.eq.1)q_o=1
          if(q_j.eq.2)q_o=3
          if(q_j.eq.3)q_o=5
          if(q_i.eq.1)q_s=1
          if(q_i.eq.2)q_s=3
          if(q_i.eq.3)q_s=5
         else
* u channel: e+ db -> l+ db
          genproc=2
          genprtyp=8
          if(q_j.eq.1)then
           if(echar.lt.0)pvalence=pvalenceu_d
          endif
          if(q_i.eq.1)q_o=-1
          if(q_i.eq.2)q_o=-3
          if(q_i.eq.3)q_o=-5
          if(q_j.eq.1)q_s=-1
          if(q_j.eq.2)q_s=-3
          if(q_j.eq.3)q_s=-5
         endif
      ENDIF
      DSIGMADXDY(1)=DSIGMADXDY(1)*pi*alpha**2/(Srad*(X*Y)**2)
      DSIGMADXDY(2)=DSIGMADXDY(2)*pi*alpha**2/(Srad*(X*Y)**2)
      DSIGMADXDY(3)=DSIGMADXDY(3)*pi*alpha**2/(Srad*(X*Y)**2)
      DSIGMADXDY(4)=DSIGMADXDY(4)*pi*alpha**2/(Srad*(X*Y)**2)
      DSIGMADXDY(1)=DSIGMADXDY(1)*0.38938*(10.0)**9
      DSIGMADXDY(2)=DSIGMADXDY(2)*0.38938*(10.0)**9
      DSIGMADXDY(3)=DSIGMADXDY(3)*0.38938*(10.0)**9
      DSIGMADXDY(4)=DSIGMADXDY(4)*0.38938*(10.0)**9
*
999   continue
*
      if(echar.lt.0)q_o=-q_o
      if(echar.lt.0)q_s=-q_s
      if(l_o.eq.1)then
* if lepton generation = 1 => load interference term  
       DXSEC(1)=DSIGMADXDY(2)
      else
       DXSEC(1)=0.d0
      endif
      DXSEC(2)=DSIGMADXDY(3)
      DXSEC(3)=DSIGMADXDY(4)
* 
      RETURN
      END
*
      subroutine LQGBAN
C-------------------------
C...Print LQGENEP banner 
C-------------------------
      implicit none

      write(6,*)
      write(6,*) ' *************************************************'
      write(6,*) '-                                                 -'
      write(6,*) '-                      LQGENEP                    -'
      write(6,*) '-                                                 -'
      write(6,*) '-                                                 -'
      write(6,*) '-                                                 -'
      write(6,*) '-    LeptoQuark GENerator for Electron-Proton     -'
      write(6,*) '-                   scattering                    -'      
      write(6,*) '-                                                 -'
      write(6,*) '-                                                 -'
      write(6,*) '-                                                 -'
      write(6,*) '- Author: L.Bellagamba                            -'
      write(6,*) '- e-mail lorenzo.bellagamba@bo.infn.it            -'
      write(6,*) '-                                                 -'
      write(6,*) '- Version: 1.0                                    -'
      write(6,*) '- Date: 01.03.2001                                -'
      write(6,*) ' *************************************************'
      write(6,*)
      write(6,*)
      write(6,*)
      return
      end
*
      subroutine LQGPRI1
C------------------------
C...Print Run requests
C------------------------
*
      implicit none
*
C...LQGENEP parameters
      double precision BEAMPAR,LQGPAR3
      integer LQGPAR1,LQGPAR2
      COMMON/LQGDATA/BEAMPAR(3),LQGPAR1(10),LQGPAR2(10),LQGPAR3(20)
*       
      write(6,*) '>>>>>>>>>>>>>>  LQGENEP RUN REQUEST  <<<<<<<<<<<<<<'
      write(6,*)
      write(6,101) LQGPAR1(1)
      write(6,1011) BEAMPAR(1)
      write(6,1012) BEAMPAR(2)
      write(6,1013) BEAMPAR(3)
      write(6,102) LQGPAR2(1)
      write(6,103) sngl(LQGPAR3(1))
      write(6,104) sngl(LQGPAR3(2))
      write(6,105) sngl(LQGPAR3(3))
      write(6,106) LQGPAR2(2)
      write(6,107) LQGPAR2(3)
      write(6,108) LQGPAR2(4)
      write(6,*)
      write(6,109) sngl(LQGPAR3(4)),sngl(LQGPAR3(5))
      write(6,110) sngl(LQGPAR3(6)),sngl(LQGPAR3(7))
      write(6,111) sngl(LQGPAR3(8))
      write(6,*)
      write(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<'
101   format('   Number of events.................. ',I8)
1011  format('   Electron beam charge.............. ',F3.0)
1012  format('   Electron beam energy.............. ',F7.2)
1013  format('   Proton beam energy................ ',F7.2)
102   format('   LQ type...........................  ',I2)
103   format('   LQ mass (GeV)..................... ',F7.2)
104   format('   LQ production coupling (s-ch.).... ',F7.4)
105   format('   LQ decay coupling (s-ch.)......... ',F7.4)
106   format('   struck quark generation (s-ch.)... ',I2)
107   format('   output quark generation (s-ch.)... ',I2)
108   format('   output lepton generation.......... ',I2)
109   format('   X generation range................ ',F6.3,' - ',F6.3)
110   format('   Y generation range................ ',F6.3,' - ',F6.3)
111   format('   minimum allowed Q2 (GeV^2)........ ',F7.2)
      return
      end
*
      subroutine LQGPRI2
C------------------------------------------
C... Print LQGENEP final statistics 
C------------------------------------------
*
      Implicit none
*
      double precision DXSEC(3),pvalence
      integer q_o,q_s,genproc,genprtyp,sch,uch,gproc(8)
      common /LQGout/ DXSEC,pvalence,q_o,q_s,genproc,genprtyp,
     >sch,uch,gproc
*
C...LQGENEP run setup parameters
      double precision BEAMPAR,LQGPAR3
      integer LQGPAR1,LQGPAR2
      COMMON/LQGDATA/BEAMPAR(3),LQGPAR1(10),LQGPAR2(10),LQGPAR3(20)

C...LQGENEP event informations
      double precision LQGKPAR,LQGDDX
      integer LQGPID
      COMMON/LQGEVT/LQGKPAR(3),LQGDDX(3),LQGPID(3)

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j  

C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar           

      Character*40 CHANS1,CHANS2,CHANU1,CHANU2
*       
      CALL LQGCHAN(CHANS1,CHANS2,CHANU1,CHANU2)
*
      write(6,*) '>>>>>>>>>>>  LQGENEP FINAL STATISTICS  <<<<<<<<<<<'
      write(6,*)
      write(6,101) LQGPAR1(4)
      write(6,104) sch
      write(6,106) uch
      write(6,*)
      write(6,*)
      write(6,*) '            Details of the generation             '
      write(6,*) '           ---------------------------            '
      write(6,*)
      write(6,*) '=================================================='
      write(6,*) 'I                                       I        I'
      write(6,*) 'I                process                I events I'
      write(6,*) 'I                                       I        I'
      write(6,*) '=================================================='
      write(6,*) 'I                                       I        I'      
      if(lqtype.eq.1.or.lqtype.eq.2)then
        write(6,111)
        write(6,110) CHANS1,gproc(3)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(5)
      elseif(lqtype.eq.3)then
        write(6,111)
        write(6,110) CHANS1,gproc(4)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(6)
      elseif(lqtype.eq.4)then
        write(6,111)
        write(6,110) CHANS1,gproc(4)
        write(6,110) CHANS2,gproc(3)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(6)
        write(6,110) CHANU2,gproc(5)
      elseif(lqtype.eq.5)then
        write(6,111)
        write(6,110) CHANS1,gproc(4)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(6)
      elseif(lqtype.eq.6)then
        write(6,111)
        write(6,110) CHANS1,gproc(4)
        write(6,110) CHANS2,gproc(3)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(6)
        write(6,110) CHANU2,gproc(5)
      elseif(lqtype.eq.7)then
        write(6,111)
        write(6,110) CHANS1,gproc(3)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(5)
      elseif(lqtype.eq.8.or.lqtype.eq.9)then
        write(6,111)
        write(6,110) CHANS1,gproc(2)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(8)
      elseif(lqtype.eq.10)then
        write(6,111)
        write(6,110) CHANS1,gproc(1)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(7)
      elseif(lqtype.eq.11)then
        write(6,111)
        write(6,110) CHANS1,gproc(2)
        write(6,110) CHANS2,gproc(1)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(8)
        write(6,110) CHANU2,gproc(7)
      elseif(lqtype.eq.12)then
        write(6,111)
        write(6,110) CHANS1,gproc(1)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(7)
      elseif(lqtype.eq.13)then
        write(6,111)
        write(6,110) CHANS1,gproc(2)
        write(6,110) CHANS2,gproc(1)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(8)
        write(6,110) CHANU2,gproc(7)
      elseif(lqtype.eq.14)then
        write(6,111)
        write(6,110) CHANS1,gproc(2)
        write(6,*) 'I        -----------------------        I        I'
        write(6,112)
        write(6,110) CHANU1,gproc(8)
      endif                                                                 
      write(6,*) 'I                                       I        I'        
      write(6,*) '=================================================='
      write(6,*)
      write(6,*) '           ---------------------------            '
      write(6,*)
      write(6,*)
      write(6,*) '>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<'
101   format('   Number of generated events.............',I8)
C 102   format('   Total cross section (mb)...............',I2)
C 103   format('   s-channel cross section (mb)....... ',F7.2)
104   format('   Number of s-channel generated events...',I8)
C 105   format('   u-channel cross section (mb)....... ',F7.2)
106   format('   Number of u-channel generated events...',I8)
110   format(1X,'I ',A38,'I',1X,I6,' I')
111   format(1X,'I ','s-channel:',28X,'I',7X,' I')
112   format(1X,'I ','u-channel:',28X,'I',7X,' I')
      return
      end      
*
      Subroutine LQGCHAN(CHANS1,CHANS2,CHANU1,CHANU2)
C-----------------------------------------------------
C...Set up the character variables containing 
C...         the generated processes
C-----------------------------------------------------
*
      Implicit None
*
      Character*38 CHANS1,CHANS2,CHANU1,CHANU2
*
C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j       
*
      Character*7 LQCHA(14)
      DATA LQCHA /'S_0L','S_0R','~S_0R','S_1L',
     >            'V_1/2L','V_1/2R','~V_1/2L',
     >            'V_0L','V_0R','~V_0R','V_1L',
     >            'S_1/2L','S_1/2R','~S_1/2L'/                        
      character*5 uch(3), dch(3)
      character*5 ubch(3), dbch(3)
      character*5 lchp(3),lchm(3)
      data uch,ubch /' u ',' c ',' t ', ' ubar',' cbar',' tbar'/
      data dch,dbch /' d ',' s ',' b ', ' dbar',' sbar',' bbar'/
      data lchp,lchm /'e+','mu+','tau+','e-','mu-','tau-'/
      character*5 l_in
      character*5 q_in
      character*5 l_ou
      character*5 q_ou
*      
      if(echar.gt.0) then
       l_in ='e+'
       l_ou=lchp(l_o)
      else
       l_in ='e-'
       l_ou=lchm(l_o)
      endif
      IF(lqtype.eq.1.or.lqtype.eq.2.or.lqtype.eq.7)Then
       if(echar.gt.0) then
        q_in=ubch(q_i)
        q_ou=ubch(q_j)
       else
        q_in=uch(q_i)
        q_ou=uch(q_j)
       endif
       CHANS1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANS2=' '
       if(echar.gt.0) then
        q_in=uch(q_j)
        q_ou=uch(q_i)
       else
        q_in=ubch(q_i)
        q_ou=ubch(q_j)
       endif
       CHANU1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANU2=' '
      ELSEIF(lqtype.eq.3.or.lqtype.eq.5)Then
       if(echar.gt.0) then
        q_in=dbch(q_i)
        q_ou=dbch(q_j)
       else
        q_in=dch(q_i)
        q_ou=dch(q_j)
       endif
       CHANS1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANS2=' '
       if(echar.gt.0) then
        q_in=dch(q_j)
        q_ou=dch(q_i)
       else
        q_in=dbch(q_j)
        q_ou=dbch(q_i)        
       endif
       CHANU1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANU2=' '
      ELSEIF(lqtype.eq.4.or.lqtype.eq.6)Then
       if(echar.gt.0) then
        q_in=dbch(q_i)
        q_ou=dbch(q_j)
       else
        q_in=dch(q_i)
        q_ou=dch(q_j)
       endif
       CHANS1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       if(echar.gt.0) then
        q_in=ubch(q_i)
        q_ou=ubch(q_j)
       else
        q_in=uch(q_i)
        q_ou=uch(q_j)
       endif
       CHANS2=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou
       if(echar.gt.0) then
        q_in=dch(q_j)
        q_ou=dch(q_i)
       else
        q_in=dbch(q_j)
        q_ou=dbch(q_i)
       endif
       CHANU1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       if(echar.gt.0) then
        q_in=uch(q_j)
        q_ou=uch(q_i)
       else
        q_in=ubch(q_j)
        q_ou=ubch(q_i)
       endif
       CHANU2=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
      ELSEIF(lqtype.eq.8.or.lqtype.eq.9.or.lqtype.eq.14)Then
       if(echar.gt.0) then
        q_in=dch(q_i)
        q_ou=dch(q_j)
       else
        q_in=dbch(q_i)
        q_ou=dbch(q_j)
       endif
       CHANS1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANS2=' '
       if(echar.gt.0) then
        q_in=dbch(q_j)
        q_ou=dbch(q_i)
       else
        q_in=dch(q_i)
        q_ou=dch(q_j)
       endif
       CHANU1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANU2=' '
      ELSEIF(lqtype.eq.10.or.lqtype.eq.12)Then
       if(echar.gt.0) then
        q_in=uch(q_i)
        q_ou=uch(q_j)
       else
        q_in=ubch(q_i)
        q_ou=ubch(q_j)
       endif
       CHANS1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANS2=' '
       if(echar.gt.0) then
        q_in=ubch(q_j)
        q_ou=ubch(q_i)
       else
        q_in=uch(q_i)
        q_ou=uch(q_j)
       endif
       CHANU1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       CHANU2=' '       
      ELSEIF(lqtype.eq.11.or.lqtype.eq.13)Then
       if(echar.gt.0) then
        q_in=dch(q_i)
        q_ou=dch(q_j)
       else
        q_in=dbch(q_i)
        q_ou=dbch(q_j)
       endif
       CHANS1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       if(echar.gt.0) then
        q_in=uch(q_i)
        q_ou=uch(q_j)
       else
        q_in=ubch(q_i)
        q_ou=ubch(q_j)
       endif
       CHANS2=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       if(echar.gt.0) then
        q_in=dbch(q_j)
        q_ou=dbch(q_i)
       else
        q_in=dch(q_j)
        q_ou=dch(q_i)
       endif
       CHANU1=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
       if(echar.gt.0) then
        q_in=ubch(q_j)
        q_ou=ubch(q_i)
       else
        q_in=uch(q_j)
        q_ou=uch(q_i)
       endif
       CHANU2=l_in//q_in//' -> '//LQCHA(lqtype)//' -> '//l_ou//q_ou       
      ENDIF      
      Return
      End
*
cc      Subroutine LQGTESTL
C-------------------------------------------
C...   Test routine to generate the process
C...   e q_1 -> S_1/2^L -> mu q_1
C...   Mass of the LQ = 250 GeV
C...   beams at HERA energies
C-------------------------------------------
ccc      Implicit none 
C...Pythia parameters
C...Parameters. 
cc      double precision PARP,PARI
cc      integer MSTP,MSTI
cc      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...LQGENEP run setup parameters
cc       double precision BEAMPAR,LQGPAR3
cc       integer LQGPAR1,LQGPAR2
cc       COMMON/LQGDATA/BEAMPAR(3),LQGPAR1(10),LQGPAR2(10),LQGPAR3(20)


C...LQGENEP event informations
cc       double precision LQGKPAR,LQGDDX
cc       integer LQGPID
cc       COMMON/LQGEVT/LQGKPAR(3),LQGDDX(3),LQGPID(3)
*

cc       integer NeV,i
*
cc       beampar(1)=1.
cc       beampar(2)=27.5
cc       beampar(3)=820.
cc       lqgpar1(1)=5000
cc       lqgpar1(2)=10
cc       lqgpar1(3)=1
cc       lqgpar1(4)=0
cc       lqgpar1(5)=0
cc       lqgpar2(1)=12
cc       lqgpar2(2)=1
cc       lqgpar2(3)=1
cc       lqgpar2(4)=2
cc       lqgpar3(1)=250.
cc       lqgpar3(2)=0.3
cc       lqgpar3(3)=0.3
cc       lqgpar3(4)=0.
cc       lqgpar3(5)=1.
cc       lqgpar3(6)=0.
cc       lqgpar3(7)=1.cc       
cc       lqgpar3(8)=500.
cc       lqgpar3(9)=1.cc       
cc       lqgpar3(10)=4.cc       
cc       lqgpar3(11)=32.cc       
* Max cross section
cc       lqgpar3(12)=3.d-6
*
* switch off initial state QCD and QED radiation
ccc       MSTP(61)=0
* switch off final state QCD and QED radiation
ccc       MSTP(71)=0
* switch off multiple interaction
ccc       MSTP(81)=0
* switch off fragmentation and decay
ccc       MSTP(111)=0

* LQGENEP Initialization
cc       call LQGENEP(0)
cc       Nev=lqgpar1(1)

* LQGENEP generation loop
cc       do i=1,Nev
cc        call LQGENEP(1)
cc       enddo

* LQGENEP termination
cc       call LQGENEP(2)

cc       return
cc       end
*
      Subroutine LQGTESTH
***************************
C-------------------------------------------
C...   Test routine to generate the process
C...   e+ q_2 -> ~S_0^R -> mu+ q_1
C...   Mass of the LQ = 600 GeV
C...   beams at HERA energies
C-------------------------------------------
      implicit none 
C...Pythia parameters
C...Parameters. 
      double precision PARP,PARI
      integer MSTP,MSTI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

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
       beampar(1)=1.
       beampar(2)=27.5
       beampar(3)=820.
       lqgpar1(1)=5000
       lqgpar1(2)=10
       lqgpar1(3)=1
       lqgpar1(4)=0
       lqgpar1(5)=0
       lqgpar2(1)=3
       lqgpar2(2)=2
       lqgpar2(3)=1
       lqgpar2(4)=3
       lqgpar3(1)=600.
       lqgpar3(2)=0.3
       lqgpar3(3)=0.3
       lqgpar3(4)=0.
       lqgpar3(5)=1.
       lqgpar3(6)=0.
       lqgpar3(7)=1.       
       lqgpar3(8)=500.
       lqgpar3(9)=1.       
       lqgpar3(10)=4.       
       lqgpar3(11)=32.       
* Max cross section
       lqgpar3(12)=2.d-11
*
* switch off initial state QCD and QED radiation
c       MSTP(61)=0
* switch off final state QCD and QED radiation
c       MSTP(71)=0
* switch off multiple interaction
c       MSTP(81)=0
* switch off fragmentation and decay
c       MSTP(111)=0

* LQGENEP Initialization
cc       call LQGENEP(0)
cc       Nev=lqgpar1(1)

* LQGENEP generation loop
cc       do i=1,Nev
cc        call LQGENEP(1)
cc       enddo       

* LQGENEP termination
cc       call LQGENEP(2)

       return
       end
*
      Subroutine LQGXNWA(XNWA,IERR)
C-------------------------------------------------------      
C...Evaluates Narrow Width Approximation Cross Section  
C...                                                      
C...  output argument: XNWA = NWA cross section           
C...                   IERR = 0 -> OK , 1 -> LQ higher 
C...                                    than center of 
C...                                    mass energy 
C-------------------------------------------------------      
      Implicit None
*
      double precision XNWA
      integer IERR
      double precision xx,cutfac0,cutfac1,ycut,factor
      double precision upv,dnv,usea,dsea,str,chm,bot,top,gl
      double precision xsu,xsd,xsub,xsdb
      double precision Br_rat
      
C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar

C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc # 
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2
*
      double precision pi
      parameter (pi=3.141592653589793d0)
*
      IERR=0
      XNWA=0.d0
      xx=MLQ**2/S
      if (xx.gt.1d0.or.xx.lt.1d-6) then
       write(6,*)'#################################################'
       write(6,*)'LQGXNWA Error:'
       write(6,*)'LQ mass not compatible with center of mass energy'
       write(6,*)'#################################################'
       IERR=1
       return
      endif
      if (Q2min.le.1.) then 
       cutfac0=1.d0
       cutfac1=2.d0
      else
       ycut=Q2min/xx/s
       if (ycut.ge.1d0) then
        XNWA=0.d0
       else
        cutfac0=1.d0-ycut
        cutfac1=2.*(cutfac0)**3
       endif
      endif
      factor=(pi*G1**2)/(4.d0*s)*(0.3894)
      call structm(xx,mlq,upv,dnv,usea,dsea,str,chm,bot,top,gl)
      if (echar.ge.0) then
       xsu=factor*(upv+usea)/xx
       xsd=factor*(dnv+dsea)/xx
       xsub=factor*usea/xx
       xsdb=factor*dsea/xx
      else
       xsu=factor*usea/xx
       xsd=factor*dsea/xx
       xsub=factor*(upv+usea)/xx
       xsdb=factor*(dnv+dsea)/xx
      endif
      IF(l_o.eq.1)then
* only electron and eventually neutrino final states considered      
       IF(lqtype.eq.1)XNWA = cutfac0*xsub*0.5
       IF(lqtype.eq.2)XNWA = cutfac0*xsub
       IF(lqtype.eq.3)XNWA = cutfac0*xsdb
       IF(lqtype.eq.4)XNWA = cutfac0*(xsub*0.5+2.*xsdb)
       IF(lqtype.eq.12)XNWA = cutfac0*xsu
       IF(lqtype.eq.13)XNWA = cutfac0*(xsu+xsd)
       IF(lqtype.eq.14)XNWA = cutfac0*xsd
       IF(lqtype.eq.5)XNWA = cutfac1*(xsdb)
       IF(lqtype.eq.6)XNWA = cutfac1*(xsub+xsdb)
       IF(lqtype.eq.7)XNWA = cutfac1*(xsub)*0.5
       IF(lqtype.eq.8)XNWA = cutfac1*(xsd)*0.5
       IF(lqtype.eq.9)XNWA = cutfac1*(xsd)
       IF(lqtype.eq.10)XNWA = cutfac1*(xsu)
       IF(lqtype.eq.11)XNWA = cutfac1*(2.*xsu+xsd*0.5)
      ELSE
* electron and muon (or tau) possible final states
       Br_rat=G2**2/(G1**2+G2**2)      
       IF(lqtype.eq.1)XNWA = Br_rat*cutfac0*xsub*0.5
       IF(lqtype.eq.2)XNWA = Br_rat*cutfac0*xsub
       IF(lqtype.eq.3)XNWA = Br_rat*cutfac0*xsdb
       IF(lqtype.eq.4)XNWA = Br_rat*cutfac0*(xsub*0.5+2.*xsdb)
       IF(lqtype.eq.12)XNWA = Br_rat*cutfac0*xsu
       IF(lqtype.eq.13)XNWA = Br_rat*cutfac0*(xsu+xsd)
       IF(lqtype.eq.14)XNWA = Br_rat*cutfac0*xsd
       IF(lqtype.eq.5)XNWA = Br_rat*cutfac1*(xsdb)
       IF(lqtype.eq.6)XNWA = Br_rat*cutfac1*(xsub+xsdb)
       IF(lqtype.eq.7)XNWA = Br_rat*cutfac1*(xsub)*0.5
       IF(lqtype.eq.8)XNWA = Br_rat*cutfac1*(xsd)*0.5
       IF(lqtype.eq.9)XNWA = Br_rat*cutfac1*(xsd)
       IF(lqtype.eq.10)XNWA = Br_rat*cutfac1*(xsu)
       IF(lqtype.eq.11)XNWA = Br_rat*cutfac1*(2.*xsu+xsd*0.5)
      ENDIF
      return
      end
*
      subroutine LQGXINT(XINT,XINTE,IERR)
C-------------------------------------------------------      
C... Cross Section evaluation by Double Differential      
C... Cross Section integration                            
C...                                                      
C... Output argument: XINT =  Integrated Cross Section    
C...                  XINTE = Error on Cross Section      
C...                  IERR =  0 -> OK , >0 -> problems in  
C...                            cross section evaluation    
C-------------------------------------------------------      
      Implicit None
*
      Double precision XINT,XINTE
      Integer IERR
*
      double precision relerr,result
      external DSDXDY
      integer minpts,maxpts,iwk,nfnevl,ifail
      parameter (minpts=10000)
      parameter (maxpts=5000000)
      parameter (iwk=1000000)
      double precision wk(iwk),a(2),b(2),eps
      parameter (eps=1d-3) 

C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

      IERR=0
      XINT=0.d0
      XINTE=0.d0
*
      a(1)=XMIN
      b(1)=XMAX
      a(2)=YMIN
      b(2)=YMAX
*      
      call dadmul(DSDXDY,2,a,b,minpts,maxpts,eps,
     +               wk,iwk,result,relerr,nfnevl,ifail)
*
      if(ifail.ne.0)then
       Write(6,*)'LQGXINT: Error in Cross Section integration'
       IERR=ifail
       return
      endif     
*
      XINT=result
      XINTE=abs(result)*relerr
*
      return
      end
*
      Double precision function dsdxdy(n,xx)
*
      IMPLICIT None
      
      integer n
      double precision xx(*)

C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc # 
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2

C# LQGout.inc #
      double precision DXSEC(3),pvalence
      integer q_o,q_s,genproc,genprtyp,sch,uch,gproc(8)
      common /LQGout/ DXSEC,pvalence,q_o,q_s,genproc,genprtyp,
     >sch,uch,gproc

* conversion from pb to mb
      double precision CONV
      DATA CONV/1.D-9/
*
      DSDXDY=0.
      X=xx(1)
      Y=xx(2)
      Q2=x*y*s        
      dsdxdy=0.d0
      if(Q2.gt.Q2min) then             
       call LQGDDXS
       DSDXDY=(DXSEC(2)+DXSEC(3))*conv
      endif
      return
      end
*
      Subroutine LQGXHMA(XSEC,XSECE,IERR)
C---------------------------------------------------------      
C...    Evaluates High Mass Approximation Cross Section     
C...                                                        
C...  Output arguments: XSEC  = HMA cross section (mb)      
C...                    XSECE = error on cross section (mb) 
C...                    IERR = 0 -> OK , >0 -> problems in   
C...                             cross section evaluation      
C...                                                        
C---------------------------------------------------------      
*
      implicit none
*     
      double precision XSEC,XSECE
      integer IERR
      external qdens

* declaration for DADMUL variables
      integer minpts,maxpts,iwk
      parameter (minpts=10000)
      parameter (maxpts=5000000)
      parameter (iwk=1000000)
      double precision wk(iwk),a(2),b(2)
      double precision eps
      parameter (eps=1d-3)  
      double precision result,relerr
      integer nfnevl,ifail
      double precision pi
      parameter (pi=3.141592653589793d0)
*
C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc # 
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2

*
      IERR=0
      XSEC=0.d0
      XSECE=0.d0
      a(1)=xmin
      b(1)=xmax
      a(2)=ymin
      b(2)=ymax
*      
      call dadmul(qdens,2,a,b,minpts,maxpts,eps,
     +            wk,iwk,result,relerr,nfnevl,ifail)         
*
      if(ifail.ne.0)then
       Write(6,*)'LQGXHMA: Error in Cross Section integration'
       IERR=ifail
       return
      endif     
      XSEC=S/32./pi*G1*G1*G2*G2/MLQ/MLQ/MLQ/MLQ*result*0.38938
      XSECE=relerr*XSEC
      return
      end
*
      double precision function qdens(n,xx)
*
      implicit none
*
      integer n
      double precision xx(*)
      double precision UPVs,DNVs,USEAs,DSEAs,STRs,
     +CHMs,BOTs,TOPs,GLs
      double precision UPVu,DNVu,USEAu,DSEAu,STRu,
     +CHMu,BOTu,TOPu,GLu
*
      double precision sh,u
      double precision UPQs(3),DNQs(3),UPQBs(3),DNQBs(3)
      double precision UPQu(3),DNQu(3),UPQBu(3),DNQBu(3)
*
C# LQGKinC.inc #
      double precision xmax,xmin,ymax,ymin,zmax,zmin,Q2min
      common /LQGKinC/ xmax,xmin,ymax,ymin,zmax,zmin,Q2min

C# LQGproc.inc #
      double precision Mlq,G1,G2
      Integer LQtype,l_o,q_i,q_j
      common /LQGproc/ Mlq,G1,G2,LQtype,l_o,q_i,q_j 

C# LQGKinV.inc # 
      double precision S,Srad,x,y,z,Q2
      common /LQGKinV/ S,Srad,x,y,z,Q2

C# LQGpar1.inc #
      integer echar
      double precision ebeam,pbeam
      common /LQGbeam/ebeam,pbeam,echar
*
      X=xx(1)
      Y=xx(2)
*      
      qdens=0.d0
      Q2=S*x*y
      sh=sqrt(S*X)
      u=sqrt(S*X*(1-Y))
      if(Q2.gt.Q2min)Then
       CALL STRUCTM(X,sh,UPVs
     &            ,DNVs,USEAs,DSEAs,STRs,CHMs,BOTs,TOPs,GLs) 
      else
       UPVs=0.
       DNVs=0.
       USEAs=0.
       DSEAs=0.
       STRs=0.
       CHMs=0.
       BOTs=0.
       TOPs=0.
       GLs=0.        
      endif
      if(Q2.gt.Q2min)Then
       CALL STRUCTM(X,u,UPVu
     &            ,DNVu,USEAu,DSEAu,STRu,CHMu,BOTu,TOPu,GLu)
      else
       UPVu=0.
       DNVu=0.
       USEAu=0.
       DSEAu=0.
       STRu=0.
       CHMu=0.
       BOTu=0.
       TOPu=0.
       GLu=0.        
      endif 
*
      if(echar.eq.1)then
* case e+, mu+, tau+ 
       UPQu(1)=UPVu+USEAu
       UPQBu(1)=USEAu
       UPQu(2)=CHMu
       UPQBu(2)=CHMu
       UPQu(3)=TOPu
       UPQBu(3)=TOPu
       DNQu(1)=DNVu+DSEAu
       DNQBu(1)=DSEAu
       DNQu(2)=STRu
       DNQBu(2)=STRu
       DNQu(3)=BOTu
       DNQBu(3)=BOTu
      elseif(echar.eq.-1)then
* case e-, mu-, tau- 
       UPQu(1)=USEAu
       UPQBu(1)=UPVu+USEAu
       UPQu(2)=CHMu
       UPQBu(2)=CHMu
       UPQu(3)=TOPu
       UPQBu(3)=TOPu
       DNQu(1)=DSEAu
       DNQBu(1)=DNVu+DSEAu
       DNQu(2)=STRu
       DNQBu(2)=STRu
       DNQu(3)=BOTu
       DNQBu(3)=BOTu
      endif          
* s channel densities
      if(echar.eq.1)then
* case e+, mu+, tau+ 
       UPQs(1)=UPVs+USEAs
       UPQBs(1)=USEAs
       UPQs(2)=CHMs
       UPQBs(2)=CHMs
       UPQs(3)=TOPs
       UPQBs(3)=TOPs
       DNQs(1)=DNVs+DSEAs
       DNQBs(1)=DSEAs
       DNQs(2)=STRs
       DNQBs(2)=STRs
       DNQs(3)=BOTs
       DNQBs(3)=BOTs
      elseif(echar.eq.-1)then
* case e-, mu-, tau- 
       UPQs(1)=USEAs
       UPQBs(1)=UPVs+USEAs
       UPQs(2)=CHMs
       UPQBs(2)=CHMs
       UPQs(3)=TOPs
       UPQBs(3)=TOPs
       DNQs(1)=DSEAs
       DNQBs(1)=DNVs+DSEAs
       DNQs(2)=STRs
       DNQBs(2)=STRs
       DNQs(3)=BOTs
       DNQBs(3)=BOTs
      endif
*
      if(LQTYPE.eq.1)Then
       qdens=UPQBs(q_i)/2.+UPQu(q_j)*(1-y)*(1-y)/2.
      elseif(LQTYPE.eq.2)Then
       qdens=UPQBs(q_i)/2.+UPQu(q_j)*(1-y)*(1-y)/2.
      elseif(LQTYPE.eq.3)Then
       qdens=DNQBs(q_i)/2.+DNQu(q_j)*(1-y)*(1-y)/2.
      elseif(LQTYPE.eq.4)Then
       if(q_j.eq.3)UPQBs(q_i)=0.
       if(q_i.eq.3)UPQu(q_j)=0.
       qdens=(UPQBs(q_i)+4.*DNQBs(q_i))/2.
     >       +(UPQu(q_j)+4.*DNQu(q_j))*(1.-Y)*(1.-Y)/2.
      elseif(LQTYPE.eq.5)Then
       qdens=DNQBs(q_i)*(1-y)*(1-y)*2.+DNQu(q_j)*2.
      elseif(LQTYPE.eq.6)Then
       if(q_j.eq.3)UPQBs(q_i)=0.
       if(q_i.eq.3)UPQu(q_j)=0.
       qdens=(UPQBs(q_i)+DNQBs(q_i))*(1-y)*(1-y)*2.
     >       +(UPQu(q_j)+DNQu(q_j))*2.
      elseif(LQTYPE.eq.7)Then
       qdens=UPQBs(q_i)*(1-y)*(1-y)*2.+UPQu(q_j)*2.
      elseif(LQTYPE.eq.8)Then
       qdens=DNQs(q_i)*2.*(1.-Y)*(1.-Y)+DNQBu(q_j)*2.
      elseif(LQTYPE.eq.9)Then
       qdens=DNQs(q_i)*2.*(1.-Y)*(1.-Y)+DNQBu(q_j)*2.
      elseif(LQTYPE.eq.10)Then  
       qdens=UPQs(q_i)*2.*(1.-Y)*(1.-Y)+UPQBu(q_j)*2.
      elseif(LQTYPE.eq.11)Then  
       if(q_j.eq.3)UPQs(q_i)=0.
       if(q_i.eq.3)UPQBu(q_j)=0.
       qdens=(4.*UPQs(q_i)+DNQs(q_i))*2.*(1.-Y)*(1.-Y)
     >      +(4.*UPQBu(q_j)+DNQBu(q_j))*2.
      elseif(LQTYPE.EQ.12)Then
       qdens=UPQs(q_i)/2.+UPQBu(q_j)*(1.-Y)*(1.-Y)/2.
      elseif(LQTYPE.EQ.13)Then
       if(q_j.eq.3)UPQs(q_i)=0.
       if(q_i.eq.3)UPQBu(q_j)=0.
       qdens=(UPQs(q_i)+DNQs(q_i))/2.
     >       +(UPQBu(q_j)+DNQBu(q_j))*(1.-Y)*(1.-Y)/2.
      elseif(LQTYPE.EQ.14)Then
       qdens=DNQs(q_i)/2.+DNQBu(q_j)*(1.-Y)*(1.-Y)/2.
      endif
      return
      end
*

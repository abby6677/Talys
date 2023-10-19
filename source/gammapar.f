      subroutine gammapar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 28, 2021
c | Task  : Gamma ray parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*6  gamchar
      character*20 wtablechar
      character*80 key
      character*90 gamfile
      integer      Zix,Nix,Z,A,N,ia,irad,l,nen,it,i
      real         eg1,sg1,gg1,eg2,sg2,gg2,egamref,enum,denom,ee,et,ft,
     +             factor,fe1(numTqrpa),fstrength,temp,dtemp,fe1t,fm1,
     +             fmax,dE,Eq,Emid,wt,wtld(6)
c
c ***************** Default giant resonance parameters *****************
c
c We use indices to indicate the type of radiation, the multipolarity
c and the number of the peak (which sometimes is two). They are
c (0(M) or 1(E), l, number), i.e. egr(0,1,1) means a constant for
c M1-radiation and egr(1,2,1) a constant for E2-radiation.
c
c 1. Read experimental E1 values from GR parameter file
c
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c ZZ,Z   : charge number of residual nucleus
c AA,A   : mass number of residual nucleus
c N      : neutron number of residual nucleus
c flagpsfglobal: flag for global photon strength functions only
c gamchar: help variable
c gamfile: giant resonance parameter file
c path   : directory containing structure files to be read
c ia     : mass number from GR table
c eg1,...: help variables
c egr    : energy of GR
c ggr    : width of GR
c gg1    : width of GR
c gg2    : width of GR
c sgr    : strength of GR
c sg1    : strength of GR
c sg2    : strength of GR
c
c GDR parameters from the table can always be overruled by a value
c given in the input file.
c
      nTqrpa=1
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      N=A-Z
      if (.not.flagpsfglobal) then
        gamchar=trim(nuc(Z))//'.gdr'
        gamfile=trim(path)//'gamma/gdr/'//gamchar
        inquire (file=gamfile,exist=lexist)
        if (.not.lexist) goto 20
        open (unit=2,file=gamfile,status='old')
   10   read(2,'(4x,i4,6f8.2)',end=20) ia,eg1,sg1,gg1,eg2,sg2,gg2
        if (A.ne.ia) goto 10
        if (egr(Zix,Nix,1,1,1).eq.0.) egr(Zix,Nix,1,1,1)=eg1
        if (sgr(Zix,Nix,1,1,1).eq.0.) sgr(Zix,Nix,1,1,1)=sg1
        if (ggr(Zix,Nix,1,1,1).eq.0.) ggr(Zix,Nix,1,1,1)=gg1
        if (egr(Zix,Nix,1,1,2).eq.0.) egr(Zix,Nix,1,1,2)=eg2
        if (sgr(Zix,Nix,1,1,2).eq.0.) sgr(Zix,Nix,1,1,2)=sg2
        if (ggr(Zix,Nix,1,1,2).eq.0.) ggr(Zix,Nix,1,1,2)=gg2
   20   close (unit=2)
      endif
c
c 1. Default GR parameterization compiled by Kopecky in RIPL
c    IAEA-TECDOC-1034, August 1998.
c
c onethird: 1/3
c pi      : pi
c kgr     : constant for gamma-ray strength function
c pi2h2c2 : 1/(pi*pi*clight*clight*hbar**2) in mb**-1.MeV**-2
c
c E1-radiation
c
      if (egr(Zix,Nix,1,1,1).eq.0.) egr(Zix,Nix,1,1,1)=
     +  31.2*A**(-onethird)+20.6*A**(-(1./6.))
      if (ggr(Zix,Nix,1,1,1).eq.0.) ggr(Zix,Nix,1,1,1)=
     +  0.026*(egr(Zix,Nix,1,1,1)**1.91)
      if (sgr(Zix,Nix,1,1,1).eq.0.) sgr(Zix,Nix,1,1,1)=
     +  1.2*120.*Z*N/(A*pi*ggr(Zix,Nix,1,1,1))
c
c E2-radiation
c
      if (egr(Zix,Nix,1,2,1).eq.0.) egr(Zix,Nix,1,2,1)=
     +  63.*A**(-onethird)
      if (ggr(Zix,Nix,1,2,1).eq.0.) ggr(Zix,Nix,1,2,1)=
     +  6.11-0.012*A
      if (sgr(Zix,Nix,1,2,1).eq.0.) sgr(Zix,Nix,1,2,1)=
     +  1.4e-4*(Z**2)*egr(Zix,Nix,1,2,1)/
     +  (A**onethird*ggr(Zix,Nix,1,2,1))
c
c E3-6 radiation
c
c gammax: number of l-values for gamma multipolarity
c
      do 110 l=3,gammax
        if (egr(Zix,Nix,1,l,1).eq.0.) egr(Zix,Nix,1,l,1)=
     +    egr(Zix,Nix,1,l-1,1)
        if (ggr(Zix,Nix,1,l,1).eq.0.) ggr(Zix,Nix,1,l,1)=
     +    ggr(Zix,Nix,1,l-1,1)
        if (sgr(Zix,Nix,1,l,1).eq.0.) sgr(Zix,Nix,1,l,1)=
     +    sgr(Zix,Nix,1,l-1,1)*8.e-4
  110 continue
c
c Check whether number of giant resonances is two (in which case
c all parameters must be specified)
c
c irad    : variable to indicate M(=0) or E(=1) radiation
c ngr     : number of GR
c strength: model for E1 gamma-ray strength function
c
      do 130 irad=0,1
        do 130 l=1,gammax
          if (egr(Zix,Nix,irad,l,2).ne.0..and.
     +      sgr(Zix,Nix,irad,l,2).ne.0..and.
     +      egr(Zix,Nix,irad,l,2).ne.0.)  ngr(Zix,Nix,irad,l)=2
  130 continue
      if (strength.le.2.or.strength.eq.5) goto 300
c
c
c ***************** HFbcs or HFB QRPA strength functions ***************
c
c For Goriely's HFbcs or HFB QRPA strength function we overwrite the
c E1 strength function with tabulated results, if available.
c
c numgamqrpa   : number of energies for QRPA strength function
c eqrpa        : energy grid for QRPA strength function
c fqrpa,fe1    : tabulated QRPA strength function
c adjust       : subroutine for energy-dependent parameter adjustment
c factor       : multiplication factor
c gamadjust    : logical for energy-dependent gamma adjustment
c etable,ftable: constant to adjust tabulated strength functions
c nTqrpa       : number of temperatures for QRPA
c Tqrpa        : temperature for QRPA
c dtemp        : temperature increment
c qrpaexist    : flag for existence of tabulated QRPA strength functions
c
      gamchar=trim(nuc(Z))//'.psf'
      if (strength.eq.3) gamfile=trim(path)//'gamma/hfbcs/'//gamchar
      if (strength.eq.4) gamfile=trim(path)//'gamma/hfb/'//gamchar
      if (strength.eq.6) gamfile=trim(path)//'gamma/hfbt/'//gamchar
      if (strength.eq.7) gamfile=trim(path)//'gamma/rmf/'//gamchar
      if (strength.eq.8) gamfile=trim(path)//'gamma/gogny/'//gamchar
      if (strength.eq.9) then
        if (flagpsfglobal) then
          gamfile=trim(path)//'gamma/smlo2019global/'//gamchar
        else
          gamfile=trim(path)//'gamma/smlo2019/'//gamchar
        endif
      endif
      if (strength.eq.10) gamfile=trim(path)//'gamma/hfbt_bsk27/'//
     &  gamchar
      inquire (file=gamfile,exist=lexist)
      if (.not.lexist) goto 210
      open (unit=2,file=gamfile,status='old')
  220 read(2,'(10x,i4)',err=210,end=210) ia
      read(2,*)
      if (ia.ne.A) then
        do 230 nen=1,numgamqrpa+1
          read(2,*)
  230   continue
        goto 220
      endif
      if (strength.eq.6.or.strength.eq.7.or.strength.eq.9.or.
     +  strength.eq.10) nTqrpa=11
      do 240 nen=1,numgamqrpa
        read(2,'(f9.3,20es12.3)') ee,(fe1(it),it=1,nTqrpa)
        if (gamadjust(Zix,Nix)) then
          key='etable'
          call adjust(ee,key,Zix,Nix,0,0,factor)
          et=etable(Zix,Nix,1,1)+etableadjust(Zix,Nix,1,1)+factor-1.
          key='ftable'
          call adjust(ee,key,Zix,Nix,0,0,factor)
          ft=ftable(Zix,Nix,1,1)*ftableadjust(Zix,Nix,1,1)+factor-1.
        else
          et=etable(Zix,Nix,1,1)+etableadjust(Zix,Nix,1,1)
          ft=ftable(Zix,Nix,1,1)*ftableadjust(Zix,Nix,1,1)
        endif
        eqrpa(Zix,Nix,nen,1,1)=ee+et
        do 250 it=1,nTqrpa
          fqrpa(Zix,Nix,nen,it,1,1)=onethird*pi2h2c2*fe1(it)*ft
  250   continue
  240 continue
      if (nTqrpa.gt.1) then
        dtemp=0.2
        temp=-dtemp
        do 260 it=1,nTqrpa
          temp=temp+dtemp
          Tqrpa(it)=temp
  260   continue
      endif
      qrpaexist(Zix,Nix,1,1)=.true.
  210 close (unit=2)
c
c ############################# M radiation ############################
c
c M1 radiation: strength function is related to that of E1
c
c egamref   : help variable
c enum,denom: help variables
c strengthM1: model for M1 gamma-ray strength function
c
  300 if (strengthM1.le.2.or.strengthM1.eq.4) then
        if (egr(Zix,Nix,0,1,1).eq.0.) egr(Zix,Nix,0,1,1)=
     +    41.*A**(-onethird)
        if (ggr(Zix,Nix,0,1,1).eq.0.) ggr(Zix,Nix,0,1,1)=4.
        if (sgr(Zix,Nix,0,1,1).eq.0.) then
          egamref=7.
          enum=kgr(1)*egamref*ggr(Zix,Nix,0,1,1)**2
          denom=(egamref**2-egr(Zix,Nix,0,1,1)**2)**2+
     +      (ggr(Zix,Nix,0,1,1)*egamref)**2
          if (strengthM1.eq.1) then
            factor=1.58e-9*A**0.47
          else
            factor=fstrength(Zix,Nix,0.,egamref,1,1)/(0.0588*A**0.878)
          endif
          sgr(Zix,Nix,0,1,1)=factor*denom/enum
        endif
      endif
      if (strengthM1.eq.3) then
c spin-flip mode
        if (egr(Zix,Nix,0,1,1).eq.0.) egr(Zix,Nix,0,1,1)=18./A**(1./6.)
        if (ggr(Zix,Nix,0,1,1).eq.0.) ggr(Zix,Nix,0,1,1)=4.
        if (sgr(Zix,Nix,0,1,1).eq.0.) sgr(Zix,Nix,0,1,1)=0.03*A**(5./6.)
c scissors mode
        if (epr(Zix,Nix,0,1,1).eq.0.) epr(Zix,Nix,0,1,1)=5./A**(0.1)
        if (gpr(Zix,Nix,0,1,1).eq.0.) gpr(Zix,Nix,0,1,1)=1.5
        if (tpr(Zix,Nix,0,1,1).eq.0.) tpr(Zix,Nix,0,1,1)=
     +    1.0e-2*abs(beta2(Zix,Nix,0))*A**(0.9)
      endif
c Includes Spin-flip from RIPL & Scissors mode from Kawano
      if (strengthM1.eq.4) then
        if (epr(Zix,Nix,0,1,1).eq.0.) epr(Zix,Nix,0,1,1)=
     +    80.*abs(beta2(Zix,Nix,0))/A**(1./3.)
        if (gpr(Zix,Nix,0,1,1).eq.0.) gpr(Zix,Nix,0,1,1)=1.5
        if (tpr(Zix,Nix,0,1,1).eq.0.) tpr(Zix,Nix,0,1,1)=
     +    42.4*beta2(Zix,Nix,0)**2/gpr(Zix,Nix,0,1,1)
      endif
c
c Tabulated M1 strength
c
c fm1: tabulated M1 strength function
c ft: tabulated strength function
c
      if (strengthM1.eq.8.or.strengthM1.eq.10) then
        gamchar=trim(nuc(Z))//'.psf'
        if (strengthM1.eq.8)
     &    gamfile=trim(path)//'gamma/gognyM1/'//gamchar
        if (strengthM1.eq.10)
     &    gamfile=trim(path)//'gamma/bsk27M1/'//gamchar
        inquire (file=gamfile,exist=lexist)
        if (.not.lexist) goto 350
        open (unit=2,file=gamfile,status='old')
  320   read(2,'(10x,i4)',end=310) ia
        read(2,*)
        if (ia.ne.A) then
          do 330 nen=1,numgamqrpa+1
            read(2,*)
  330     continue
          goto 320
        endif
        do 340 nen=1,numgamqrpa
          read(2,'(f9.3,20es12.3)') ee,fm1
          if (gamadjust(Zix,Nix)) then
            key='etable'
            call adjust(ee,key,Zix,Nix,0,0,factor)
            et=etable(Zix,Nix,0,1)+etableadjust(Zix,Nix,0,1)+factor-1.
            key='ftable'
            call adjust(ee,key,Zix,Nix,0,0,factor)
            ft=ftable(Zix,Nix,0,1)*ftableadjust(Zix,Nix,0,1)+factor-1.
          else
            et=etable(Zix,Nix,0,1)+etableadjust(Zix,Nix,0,1)
            ft=ftable(Zix,Nix,0,1)*ftableadjust(Zix,Nix,0,1)
          endif
          eqrpa(Zix,Nix,nen,0,1)=ee+et
          fqrpa(Zix,Nix,nen,1,0,1)=onethird*pi2h2c2*fm1*ft
  340   continue
c
c Avoid having an increasing function for extrapolation purposes
c
        if (fqrpa(Zix,Nix,numgamqrpa,1,0,1).gt.
     +    fqrpa(Zix,Nix,numgamqrpa-1,1,0,1))
     +  fqrpa(Zix,Nix,numgamqrpa,1,0,1)=
     +  fqrpa(Zix,Nix,numgamqrpa-1,1,0,1)
        qrpaexist(Zix,Nix,0,1)=.true.
  310   close (unit=2)
c
c add some scissors mode contribution to spherical QRPA calculation 
c
        if (strengthM1.eq.10) then
          if (epr(Zix,Nix,0,1,1).eq.0.) epr(Zix,Nix,0,1,1)=3.
          if (gpr(Zix,Nix,0,1,1).eq.0.) gpr(Zix,Nix,0,1,1)=2.
          if (tpr(Zix,Nix,0,1,1).eq.0.) tpr(Zix,Nix,0,1,1)=
     &      7.e-5*beta2(Zix,Nix,0)**2*A**2
        endif
      endif
c
c M2-6 radiation
c
  350 do 360 l=2,gammax
        if (egr(Zix,Nix,0,l,1).eq.0.) egr(Zix,Nix,0,l,1)=
     +    egr(Zix,Nix,0,l-1,1)
        if (ggr(Zix,Nix,0,l,1).eq.0.) ggr(Zix,Nix,0,l,1)=
     +    ggr(Zix,Nix,0,l-1,1)
        if (sgr(Zix,Nix,0,l,1).eq.0.) sgr(Zix,Nix,0,l,1)=
     +    sgr(Zix,Nix,0,l-1,1)*8.e-4
  360 continue
c
c External strength functions
c
c fe1t: tabulated strength function
c
      do 410 irad=0,1
        do 420 l=1,gammax
          if (Exlfile(Zix,Nix,irad,l)(1:1).ne.' ') then
            nen=0
            open (unit=2,file=Exlfile(Zix,Nix,irad,l),status='old')
  430       read(2,*,err=440,end=500) ee,fe1t
            if (gamadjust(Zix,Nix)) then
              key='etable'
              call adjust(ee,key,Zix,Nix,0,0,factor)
              et=etable(Zix,Nix,irad,l)+etableadjust(Zix,Nix,irad,l)+
     +          factor-1.
              key='ftable'
              call adjust(ee,key,Zix,Nix,0,0,factor)
              ft=ftable(Zix,Nix,irad,l)*ftableadjust(Zix,Nix,irad,l)+
     +          factor-1.
            else
              et=etable(Zix,Nix,irad,l)+etableadjust(Zix,Nix,irad,l)
              ft=ftable(Zix,Nix,irad,l)*ftableadjust(Zix,Nix,irad,l)
            endif
            nen=nen+1
            eqrpa(Zix,Nix,nen,irad,l)=ee+et
            fqrpa(Zix,Nix,nen,1,irad,l)=onethird*pi2h2c2*fe1t*ft
  440       if (nen.lt.numgamqrpa) goto 430
  500       if (nen.gt.0) qrpaexist(Zix,Nix,irad,l)=.true.
            close (unit=2)
          endif
  420   continue
  410 continue
c
c Reading of width of tabulated PSF
c
      if (flagngfit.and.k0.eq.1.and..not.flagpsfglobal.and.
     +  (strength.eq.8.or.strength.eq.9)) then
        if (.not.gamadjust(Zix,Nix)) then
          wtablechar=trim(nuc(Z))//'.wtable'
          if (strength.eq.8) then
            gamfile=trim(path)//'gamma/gogny/'//wtablechar
          else
            gamfile=trim(path)//'gamma/smlo2019/'//wtablechar
          endif
          inquire (file=gamfile,exist=lexist)
          if (.not.lexist) goto 480
          open (unit=2,file=gamfile,status='old')
  460     read(2,'(4x,i4,6f10.5)',end=480) ia,(wtld(i),i=1,6)
          if (A.ne.ia) goto 460
          wt=wtld(ldmodel(Zix,Nix))
          if (wt.gt.0.) wtable(Zix,Nix,1,1)=wt
  480     close (unit=2)
        endif
      endif
c
c Adjustment of width of tabulated PSF
c
      do 510 irad=0,1
        do 520 l=1,gammax
          wt=wtable(Zix,Nix,irad,l)*wtableadjust(Zix,Nix,irad,l)
          Emid=0.
          fmax=0.
          do 530 nen=1,numgamqrpa
            if (fqrpa(Zix,Nix,nen,1,irad,l).gt.fmax) then
              fmax=fqrpa(Zix,Nix,nen,1,irad,l)
              Emid=eqrpa(Zix,Nix,nen,irad,l)
            endif
  530     continue
          do 540 nen=1,numgamqrpa
            Eq=eqrpa(Zix,Nix,nen,irad,l)
            dE=Eq-Emid
            eqrpa(Zix,Nix,nen,irad,l)=Emid+dE*wt
  540     continue
  520   continue
  510 continue
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

      subroutine omppar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 31, 2020
c | Task  : Optical model parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*1  omptype,Rpart
      character*8  ompchar
      character*72 optmodfile
      character*132 ompfile
      character*132 string
      integer      Zix,Nix,Z,N,A,k,iz,ia,i,nomp,ii,nen,kk,j,NE,Zt,At,
     +             Zbeg,Zend,Abeg,Aend,Rnum
      real         e,eferm,dum,Er,Eeps,Eripl(numen),dEripl,Ebeg,Efin
c
c ************** Read optical model parameters from database ***********
c
c Zix         : charge number index for residual nucleus
c Nix         : neutron number index for residual nucleus
c ZZ,Z        : charge number of residual nucleus
c NN,N        : neutron number of residual nucleus
c AA,A        : mass number of residual nucleus
c omptype     : type of optical model (spherical or coupled)
c flaglocalomp: flag for local (y) or global (n) optical model
c ompchar     : help variable
c ompfile     : optical model parameter file
c optmodfileN : optical model parameter file for neutrons
c optmodfileP : optical model parameter file for protons
c path        : directory containing structure files to be read
c ia          : mass number from level file
c nomp        : number of particles for optical model parameters
c ef          : Fermi energy
c rc0         : Coulomb radius
c rv0,av0     : real volume radius, diffuseness
c v1,v2,v3    : components for V
c w1,w2       : components for W
c rvd0,avd0   : real surface radius, diffuseness
c d1,d2,d3    : components for Wd
c rvso0,avso0 : real spin-orbit radius, diffuseness
c vso1,vso2   : components for Vso
c wso1,wso2   : components for Wso
c flagdisp    : flag for dispersive optical model
c flagjlm     : flag for using semi-microscopic JLM OMP
c disp        : flag for dispersive optical model
c colltype    : type of collectivity (D, V or R)
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      omptype=' '
      if (.not.flaglocalomp) goto 100
      do 10 k=1,2
        ompchar=parsym(k)//'-'//trim(nuc(Z))//'.omp'
        if (k.eq.1) then
          if (optmodfileN(Zix)(1:1).ne.' ') then
            ompfile=optmodfileN(Zix)
          else
            ompfile=trim(path)//'optical/neutron/'//ompchar
          endif
        else
          if (optmodfileP(Zix)(1:1).ne.' ') then
            ompfile=optmodfileP(Zix)
          else
            ompfile=trim(path)//'optical/proton/'//ompchar
          endif
        endif
        omptype=' '
        inquire (file=ompfile,exist=lexist)
        if (.not.lexist) goto 10
        open (unit=2,file=ompfile,status='old')
   20   read(2,'(4x,2i4,3x,a1)',end=60) ia,nomp,omptype
        if (A.ne.ia) then
          do 30 i=1,nomp
            do 30 ii=1,4
              read(2,'()')
   30     continue
          goto 20
        endif
        omptype=' '
        do 40 i=1,nomp
          read(2,'(4x,f7.2,f8.3)') ef(Zix,Nix,k),rc0(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f9.6,f6.1,f7.1)') rv0(Zix,Nix,k),
     +      av0(Zix,Nix,k),v1(Zix,Nix,k),v2(Zix,Nix,k),v3(Zix,Nix,k),
     +      w1(Zix,Nix,k),w2(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f7.2)') rvd0(Zix,Nix,k),
     +      avd0(Zix,Nix,k),d1(Zix,Nix,k),d2(Zix,Nix,k),d3(Zix,Nix,k)
          read(2,'(2f8.3,f6.1,f10.4,f6.1,f7.1)') rvso0(Zix,Nix,k),
     +      avso0(Zix,Nix,k),vso1(Zix,Nix,k),vso2(Zix,Nix,k),
     +      wso1(Zix,Nix,k),wso2(Zix,Nix,k)
          if (flagdisp) then
            if (nomp.eq.1.or.flagjlm) then
              disp(Zix,Nix,k)=.false.
            else
              disp(Zix,Nix,k)=.true.
            endif
          else
            goto 60
          endif
   40   continue
   60   close (unit=2)
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,k)=0.85*d1(Zix,Nix,k)
   10 continue
c
c ************************* Global optical model ***********************
c
c 1. Neutrons
c
c ompglobal: flag for use of global optical model
c onethird : 1/3
c twothird : 2/3
c
c Test if local OMP has been assigned.
c
  100 if (rv0(Zix,Nix,1).eq.0.) then
        if (flagdisp) disp(Zix,Nix,1)=.false.
        ompglobal(Zix,Nix,1)=.true.
        ef(Zix,Nix,1)=-11.2814+0.02646*A
        rv0(Zix,Nix,1)=1.3039-0.4054*A**(-onethird)
        av0(Zix,Nix,1)=0.6778-1.487e-4*A
        v1(Zix,Nix,1)=59.30-21.0*real(N-Z)/A-0.024*A
        v2(Zix,Nix,1)=7.228e-3-1.48e-6*A
        v3(Zix,Nix,1)=1.994e-5-2.0e-8*A
        w1(Zix,Nix,1)=12.195+0.0167*A
        w2(Zix,Nix,1)=73.55+0.0795*A
        rvd0(Zix,Nix,1)=1.3424-0.01585*A**onethird
        avd0(Zix,Nix,1)=0.5446-1.656e-4*A
        d1(Zix,Nix,1)=16.0-16.0*real(N-Z)/A
        d2(Zix,Nix,1)=0.0180+3.802e-3/(1.+exp((A-156.)/8.0))
        d3(Zix,Nix,1)=11.5
        vso1(Zix,Nix,1)=5.922+0.0030*A
        vso2(Zix,Nix,1)=0.0040
        rvso0(Zix,Nix,1)=1.1854-0.647*A**(-onethird)
        avso0(Zix,Nix,1)=0.59
        wso1(Zix,Nix,1)=-3.1
        wso2(Zix,Nix,1)=160.
        rc0(Zix,Nix,1)=0.
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,1)=0.85*d1(Zix,Nix,1)
      endif
c
c 2. Protons
c
      if (rv0(Zix,Nix,2).eq.0.) then
        if (flagdisp) disp(Zix,Nix,2)=.false.
        ompglobal(Zix,Nix,2)=.true.
        ef(Zix,Nix,2)=-8.4075+0.01378*A
        rv0(Zix,Nix,2)=1.3039-0.4054*A**(-onethird)
        av0(Zix,Nix,2)=0.6778-1.487e-4*A
        v1(Zix,Nix,2)=59.30+21.0*real(N-Z)/A-0.024*A
        v2(Zix,Nix,2)=7.067e-3+4.23e-6*A
        v3(Zix,Nix,2)=1.729e-5+1.136e-8*A
        w1(Zix,Nix,2)=14.667+0.009629*A
        w2(Zix,Nix,2)=73.55+0.0795*A
        rvd0(Zix,Nix,2)=1.3424-0.01585*A**onethird
        avd0(Zix,Nix,2)=0.5187+5.205e-4*A
        d1(Zix,Nix,2)=16.0+16.0*real(N-Z)/A
        d2(Zix,Nix,2)=0.0180+3.802e-3/(1.+exp((A-156.)/8.0))
        d3(Zix,Nix,2)=11.5
        vso1(Zix,Nix,2)=5.922+0.0030*A
        vso2(Zix,Nix,2)=0.0040
        rvso0(Zix,Nix,2)=1.1854-0.647*A**(-onethird)
        avso0(Zix,Nix,2)=0.59
        wso1(Zix,Nix,2)=-3.1
        wso2(Zix,Nix,2)=160.
        rc0(Zix,Nix,2)=1.198+0.697*A**(-twothird)+12.994*A**(-5./3.)
c
c Reduce d1 parameter (of Wd) in case of coupled-channels, unless
c already specified in OMP parameterization.
c
        if (colltype(Zix,Nix).ne.'S'.and.omptype.ne.'C')
     +    d1(Zix,Nix,2)=0.85*d1(Zix,Nix,2)
      endif
c
c ************** Optical model parameters from RIPL ********************
c
      if (Zix.le.numZph.and.Nix.le.numNph) then
        if (flagriplomp) then
c
c Copy RIPL parameter files to working directory
c
          ompfile=trim(path)//'optical/ripl/om-parameter-u.dat'
          open (unit=41,file=ompfile,status='unknown')
          open (unit=42,file='om-parameter-u.dat',status='unknown')
  110     read(41,'(a)',end=120) string
          write(42,'(a)') trim(string)
          goto 110
  120     close (unit=42)
          close (unit=41)
          ompfile=trim(path)//'optical/ripl/gs-mass-sp.dat'
          open (unit=41,file=ompfile,status='unknown')
          open (unit=42,file='gs-mass-sp.dat',status='unknown')
  130     read(41,'(a)',end=140) string
          write(42,'(a)') trim(string)
          goto 130
  140     close (unit=42)
          close (unit=41)
c
c Set energy grid for OMP table
c
          Er=0.
          nen=0
          dEripl=0.001
  150     Er=Er+dEripl
          Eeps=Er+1.e-4
          if (Er.gt.enincmax+12.) goto 160
          if (nen.eq.numen) goto 160
          nen=nen+1
          Eripl(nen)=Er
          if (Eeps.gt.0.01) dEripl=0.01
          if (Eeps.gt.0.1) dEripl=0.1
          if (Eeps.gt.4.) dEripl=0.2
          if (Eeps.gt.10.) dEripl=0.5
          if (Eeps.gt.30.) dEripl=1.
          if (Eeps.gt.100.) dEripl=2.
          goto 150
  160     NE=nen
          do 194 k=1,6
            if (riplomp(k).gt.0.and.Zix.eq.parZ(k).and.Nix.eq.parN(k).
     +        and.optmod(Zix,Nix,k).eq.' ') then
c
c Test if RIPL number is correct and inside mass and energy ranges 
c
              Zt=ZZ(0,0,k)
              At=AA(0,0,k)
              ompfile=trim(path)//'optical/ripl/om-index.txt'
              open (unit=31,file=ompfile,status='old')
              read(31,'(//)')
  170         read(31,'(i4,3x,a1,23x,i2,1x,i2,2x,i3,1x,i3,1x,
     +          f5.1,1x,f5.1)',end=180) 
     +          Rnum,Rpart,Zbeg,Zend,Abeg,Aend,Ebeg,Efin
              if (riplomp(k).eq.Rnum) then
                if (parsym(k).eq.Rpart.and.Zt.ge.Zbeg.and.Zt.le.Zend.
     +            and.At.ge.Abeg.and.At.le.Aend) then
                  Eompbeg1(k,1)=Ebeg
                  Eompend1(k,1)=Efin
                  Eompbeg0(k,1)=Ebeg*0.8
                  Eompend0(k,1)=Efin*1.2
                  close (unit=31)
                  goto 190
                endif
              endif
              goto 170
  180         close (unit=31)
              write(*,'(" TALYS-warning: RIPL OMP ",i4," for ",i3,a2,
     +          " out of range")') riplomp(k),At,nuc(Zt)
              if (.not.flagriplrisk) stop
c
c Create input file for RIPL OMP retrieval code
c
  190         open (unit=42,file='ominput.inp',status='unknown')
              write(42,*) NE
              write(42,*) (Eripl(i),i=1,NE)
              write(42,*) Zt,At,riplomp(k),-2
              close (unit=42)
c
c Run RIPL OMP retrieval code
c
              write(*,*) "RIPL OMP for ",Zt,At
              call riplomp_mod(NE,Zt,At,Eripl,riplomp(k),-2,0)
              close (unit=35)
c
c Retrieve  RIPL OMP parameters on energy grid
c
              e=0.
              open(42,file='omp-table.dat',status='unknown')
  192         read(42,'(a)') string
              kk=index(string,'Ef =')
              if (kk > 0) read(string(kk+4:132),*) ef(Zix,Nix,k)
              read(string(1:7),*,err=192,end=192) e
              if (e.gt.0.) then
                backspace 42
                do i=1,NE 
                  read(42,*) eomp(Zix,Nix,k,i),
     +              (vomp(Zix,Nix,k,i,j),j=1,3),dum,
     +              (vomp(Zix,Nix,k,i,j),j=4,18)
                enddo
                close(42)
                omplines(Zix,Nix,k)=NE
              else
                goto 192
              endif
            endif
  194     enddo
        endif
c
c ************** Optical model file from user input file ***************
c
c numNph    : maximal number of neutrons away from the initial
c             compound nucleus for multiple pre-equilibrium emission
c numZph    : maximal number of protons away from the initial
c             compound nucleus for multiple pre-equilibrium emission
c optmod    : file with optical model parameters
c optmodfile: file with optical model parameters
c omplines  : number of lines in optical model file
c numomp    : maximum number of lines in optical model file
c eomp      : energies on optical model file
c vomp      : optical model parameters from file
c
        do 210 k=1,6
          optmodfile=optmod(Zix,Nix,k)
          if (optmodfile(1:1).ne.' ') then
            open (unit=2,file=optmodfile,status='old')
            eferm=0.
            read(2,'(a)',end=300,err=300) string
            read(string,*,end=240,err=300) iz,ia,omplines(Zix,Nix,k),
     +        eferm
  240       if (eferm.gt.0.) ef(Zix,Nix,k)=eferm
            if (omplines(Zix,Nix,k).gt.numomp) then
              write(*,'(" TALYS-error: number of lines in OMP file > ",
     +          i4)') numomp
              stop
            endif
            eomp(Zix,Nix,k,0)=0.
            do 220 nen=1,omplines(Zix,Nix,k)
              read(2,*,err=300,end=220) eomp(Zix,Nix,k,nen),
     +          (vomp(Zix,Nix,k,nen,ii),ii=1,19)
  220       continue
            close (unit=2)
            if (omplines(Zix,Nix,k).eq.1) then
              omplines(Zix,Nix,k)=2
              eomp(Zix,Nix,k,2)=eomp(Zix,Nix,k,1)
              do 230 ii=1,19
                vomp(Zix,Nix,k,2,ii)=vomp(Zix,Nix,k,1,ii)
  230         continue
            endif
          endif
  210   continue
      endif
c
c Set OMP parameters for extension up to 1 GeV
c
c Zindex  : charge number index for residual nucleus
c Nindex  : neutron number index for residual nucleus
c Ejoin   : joining energy for high energy OMP
c enincmax: maximum incident energy
c optical : subroutine for determination of optical potential
c V0      : V at zero MeV
c Vjoin   : V at joining energy
c Wjoin   : W at joining energy
c
      do 310 k=1,2
        w3(Zix,Nix,k)=25.-0.0417*A
        w4(Zix,Nix,k)=250.
        if (Zix.eq.Zindex(0,0,k).and.Nix.eq.Nindex(0,0,k).and.
     +    enincmax.gt.Ejoin(k)) then
          e=0.
          call optical(Zix,Nix,k,e)
          V0(k)=v
          e=Ejoin(k)
          call optical(Zix,Nix,k,e)
          Vjoin(k)=v
          Wjoin(k)=w
        endif
  310 continue
c
c Set energy ranges of alternative optical models
c
c Eompbeg0: lower energy of KD03 OMP
c Eompbeg1: lower energy of alternative OMP
c Eompend1: upper energy of alternative OMP
c Eompend0: upper energy of KD03 OMP
c
c Eompbeg0 <=  Eompbeg1 <=  Eompend1 <=  Eompend0
c
c
c Deuteron OMPs
c
      do 410 i=2,5
        Eompbeg0(3,i)=0.
        Eompbeg1(3,i)=0.
        Eompend1(3,i)=200.
        Eompend0(3,i)=300.
  410 continue
      Eompend1(3,2)=90.
      Eompend0(3,2)=150.
      Eompend1(3,3)=100.
      Eompend0(3,3)=150.
c
c Alpha OMPs
c
      do 420 i=2,8
        Eompbeg0(6,i)=0.
        Eompbeg1(6,i)=0.
        Eompend1(6,i)=200.
        Eompend0(6,i)=300.
  420 continue
      Eompend1(6,2)=25.
      Eompend0(6,2)=50.
      return
  300 write(*,'(" TALYS-error: Format error in ",a72)') optmodfile
      stop
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

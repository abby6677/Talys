      subroutine gammaout(Zcomp,Ncomp)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : January 8, 2021
c | Task  : Output of gamma-ray strength functions, transmission
c |         coefficients and cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1  radtype
      character*12 psffile
      character*25 model
      integer      Zcomp,Ncomp,i,Z,N,A,l,nen,irad,Npsf,nenb,nene,MM
      real         Epsf(numen),e,fstrength
c
c ***** Gamma-ray strength functions and transmission coefficients *****
c
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c ZZ,Z      : charge number of residual nucleus
c NN,N      : neutron number of residual nucleus
c AA,A      : mass number of residual nucleus
c nuc       : symbol of nucleus
c gamgam    : total radiative width
c dgamgam   : uncertainty in gamgam
c gamgamth  : theoretical total radiative width
c D0        : experimental s-wave resonance spacing in eV
c dD0       : uncertainty in D0
c D0theo    : theoretical s-wave resonance spacing
c swaveth   : theoretical strength function for s-wave
c gnorm     : gamma normalization factor
c strength  : model for E1 gamma-ray strength function
c model     : string for gamma-ray strength function
c strengthM1: model for M1 gamma-ray strength function
c etable    : constant to adjust tabulated strength functions
c ftable    : constant to adjust tabulated strength functions
c wtable    : constant to adjust tabulated strength functions
c flagupbend: flag for low-energy upbend of photon strength function
c upbendc   : properties of the low-energy upbend of given multipolarity
c nTqrpa    : number of temperatures for QRPA
c gammax    : number of l-values for gamma multipolarity
c sgr       : strength of GR
c ngr       : number of GR
c egr       : energy of GR
c ggr       : width of GR
c epr       : energy of PR
c gpr       : width of PR
c tpr       : strength of PR
c kgr       : constant for gamma-ray strength function
c ebegin    : first energy point of energy grid
c eend      : last energy point of energy grid
c egrid,e   : outgoing energy grid
c fstrength : gamma-ray strength function
c Tjl       : transmission coefficients as a function of particle type,
c             energy, spin and l-value
c filepsf   : flag for photon strength functions on separate files
c
      Z=ZZ(Zcomp,Ncomp,0)
      N=NN(Zcomp,Ncomp,0)
      A=AA(Zcomp,Ncomp,0)
      if (filepsf) then
        i=0
        do 10 nen=ebegin(0),eend(0)
          i=i+1
          Epsf(i)=egrid(nen)
   10   continue
        Npsf=i
        if (Npsf.gt.0) then
          if (Epsf(Npsf).lt.30.) then
            nenb=int(2.*(Epsf(Npsf)+0.5))
            nene=60
            do 20 nen=nenb,nene
              i=i+1
              if (i.le.numen) Epsf(i)=0.5*nen
   20       continue
            Npsf=min(i,numen)
          endif
        endif
      endif
      write(*,'(/" ########## GAMMA STRENGTH FUNCTIONS, TRANSMISSION",
     +  " COEFFICIENTS AND CROSS SECTIONS ##########")')
      write(*,'(/" Gamma-ray information for Z=",i3," N=",i3,
     +  " (",i3,a2,") "/)') Z,N,A,nuc(Z)
      write(*,'(" S-wave strength function parameters:"/)')
      write(*,'(" Exp. total radiative width=",f10.5," eV +/-",f8.5,
     +  " Theor. total radiative width for l=0:",f15.5," eV")')
     +  gamgam(Zcomp,Ncomp),dgamgam(Zcomp,Ncomp),gamgamth(Zcomp,Ncomp,0)
      write(*,'(53x," Theor. total radiative width for l=0:",f15.5,
     +  " eV (gnorm=",f6.3,")")') gamgamth(Zcomp,Ncomp,0)*gnorm,gnorm
      write(*,'(53x," Theor. total radiative width for l=1:",f15.5,
     +  " eV")') gamgamth(Zcomp,Ncomp,1)
      write(*,'(" Exp. D0                   =",f10.2," eV +/-",f8.2,
     +  " Theor. D0                   =",f15.5," eV")')
     +  D0(Zcomp,Ncomp),dD0(Zcomp,Ncomp),D0theo(Zcomp,Ncomp)
      write(*,'(53x," Theor. D1                   =",f15.5," eV")')
     +  D1theo(Zcomp,Ncomp)
      write(*,'(" Theor. S-wave strength f. =",f10.5,"E-4")')
     +  1.e4*swaveth(Zcomp,Ncomp)
      write(*,'(" Average resonance energy  =",f13.2," eV")') 
     +  Eavres*1.e6
      write(*,'(" Normalization factor      =",f10.5)') gnorm
      write(*,'(/" Incident energy: E[MeV]=",f8.3)') Einc
      if (strength.eq.1) model="Kopecky-Uhl              "
      if (strength.eq.2) model="Brink-Axel               "
      if (strength.eq.3) model="Goriely HFbcs tables     "
      if (strength.eq.4) model="Goriely HFB tables       "
      if (strength.eq.5) model="Goriely Hybrid model     "
      if (strength.eq.6) model="Goriely T-dep. HFB Tables"
      if (strength.eq.7) model="Goriely T-dep. RMF Tables"
      if (strength.eq.8) model="Gogny D1M HFB+QRPA Tables"
      if (strength.eq.9) model="IAEA-CRP SMLO 2019 Tables"
      if (strength.eq.10) model="BSk27+QRPA 2018 Tables   "
      write(*,'(/" Gamma-ray strength function model for E1: ",a25)')
     +  model
      if (strengthM1.eq.1) model="RIPL-1                   "
      if (strengthM1.eq.2) model="RIPL-2                   "
      if (strengthM1.eq.3) model="IAEA GSF CRP (2018)      "
      if (strengthM1.eq.4) model="RIPL-2+Scissors Kawano   "
      if (strengthM1.eq.8) model="Gogny D1M HFB+QRPA Tables"
      if (strengthM1.eq.10) model="BSk27+QRPA 2018 Tables   "
      write(*,'(/" Gamma-ray strength function model for M1: ",a25)')
     +  model
      if (strength.eq.3.or.strength.eq.4.or.strength.ge.6) then
        write(*,'(/" Adjustable parameters for E1: etable=",f10.5,
     +    " ftable=",f10.5," wtable=",f10.5," number of T=",i3)')
     +    etable(Zcomp,Ncomp,1,1),ftable(Zcomp,Ncomp,1,1),
     +    wtable(Zcomp,Ncomp,1,1),nTqrpa
      endif
      if (strengthM1.eq.8.or.strengthM1.eq.10) then
        write(*,'(/" Adjustable parameters for M1: etable=",f10.5,
     +    " ftable=",f10.5," wtable=",f10.5)')
     +    etable(Zcomp,Ncomp,0,1),ftable(Zcomp,Ncomp,0,1),
     +    wtable(Zcomp,Ncomp,0,1)
      endif
      if (flagupbend) then
        write(*,'(/" Inclusion of an E1 upbend C x U*",
     +    "/ (1+exp(E-eta)) with  C=",es10.2," eta=",f8.2)')
     +    upbend(Zcomp,Ncomp,1,1,1),upbend(Zcomp,Ncomp,1,1,2)
        write(*,'(/" Inclusion of an M1 upbend C exp(-F*|beta2|) ",
     +    "exp(-eta*E) with C=",es10.2," eta=",f8.2," F=",f8.2)')
     +    upbend(Zcomp,Ncomp,0,1,1),upbend(Zcomp,Ncomp,0,1,2),
     +    upbend(Zcomp,Ncomp,0,1,3)
      endif
      do 30 l=1,gammax
        write(*,'(/" Normalized gamma-ray strength functions and ",
     +    "transmission coefficients for l=",i2,/)') l
        write(*,'(" Giant resonance parameters :"/)')
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'(" sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3," and ",f8.3,
     +      "    PR: sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3)')
     +      l,sgr(Zcomp,Ncomp,0,l,1),l,
     +      sgr(Zcomp,Ncomp,1,l,1),sgr(Zcomp,Ncomp,1,l,2),
     +      l,tpr(Zcomp,Ncomp,0,l,1),l,tpr(Zcomp,Ncomp,1,l,1)
        else
          write(*,'(" sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3,
     +      "    PR: sigma0(M",i1,") =",f8.3,"       sigma0(E",i1,
     +      ") =",f8.3)')
     +      l,sgr(Zcomp,Ncomp,0,l,1),l,sgr(Zcomp,Ncomp,1,l,1),
     +      l,tpr(Zcomp,Ncomp,0,l,1),l,tpr(Zcomp,Ncomp,1,l,1)
        endif
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'("      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3," and ",f8.3,
     +    "    PR:      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3)')
     +      l,egr(Zcomp,Ncomp,0,l,1),l,
     +      egr(Zcomp,Ncomp,1,l,1),egr(Zcomp,Ncomp,1,l,2),
     +      l,epr(Zcomp,Ncomp,0,l,1),l,epr(Zcomp,Ncomp,1,l,1)
        else
          write(*,'("      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3,
     +    "    PR:      E(M",i1,") =",f8.3,"            E(E",i1,
     +      ") =",f8.3)')
     +      l,egr(Zcomp,Ncomp,0,l,1),l,egr(Zcomp,Ncomp,1,l,1),
     +      l,epr(Zcomp,Ncomp,0,l,1),l,epr(Zcomp,Ncomp,1,l,1)
        endif
        if (ngr(Zcomp,Ncomp,1,l).eq.2) then
          write(*,'("  gamma(M",i1,") =",f8.3,"        gamma(E",i1,
     +      ") =",f8.3," and ",f8.3,
     +      "    PR:  gamma(M",i1,") =",f8.3,"        gamma(E",i1,") =",
     +      f8.3)')
     +      l,ggr(Zcomp,Ncomp,0,l,1),l,
     +      ggr(Zcomp,Ncomp,1,l,1),ggr(Zcomp,Ncomp,1,l,2),
     +      l,gpr(Zcomp,Ncomp,0,l,1),l,gpr(Zcomp,Ncomp,1,l,1)
        else
          write(*,'("  gamma(M",i1,") =",f8.3,"        gamma(E",i1,
     +      ") =",f8.3,
     +      "    PR:  gamma(M",i1,") =",f8.3,"        gamma(E",i1,") =",
     +      f8.3)')
     +      l,ggr(Zcomp,Ncomp,0,l,1),l,ggr(Zcomp,Ncomp,1,l,1),
     +      l,gpr(Zcomp,Ncomp,0,l,1),l,gpr(Zcomp,Ncomp,1,l,1)
        endif
        write(*,'("      k(M",i1,") =",es14.5,"      k(E",i1,") =",
     +    es14.5/)') l,kgr(l),l,kgr(l)
        write(*,'("      E       f(M",i1,")        f(E",i1,")",
     +    "        T(M",i1,")        T(E",i1,")"/)')  l,l,l,l
        do 40 nen=ebegin(0),eend(0)
          e=egrid(nen)
          write(*,'(1x,f8.3,4es13.5)') e,
     +      fstrength(Zcomp,Ncomp,0.,e,0,l)*gnorm,
     +      fstrength(Zcomp,Ncomp,0.,e,1,l)*gnorm,Tjl(0,nen,0,l),
     +      Tjl(0,nen,1,l)
   40   continue
c
c Output on separate files
c
        if (filepsf.and.l.eq.1) then
          do 50 irad=0,1
            psffile='psf000000.E1'
            write(psffile(4:9),'(2i3.3)') Z,A
            if (irad.eq.0) then
              radtype='M'
              MM=strengthM1
            else
              radtype='E'
              MM=strength
            endif
            psffile(11:11)=radtype
            open (unit=1,file=psffile,status='replace')
            write(1,'("#     E      f(",a1,"1)     Model: ",i2)') 
     +        radtype,MM
            do 60 nen=1,Npsf
              e=Epsf(nen)
              write(1,'(1x,f8.3,es13.5)') e,
     +          fstrength(Zcomp,Ncomp,0.,e,irad,1)*gnorm
   60       continue
            close (1)
   50     continue
        endif
   30 continue
c
c **************** Cross sections for inverse channels *****************
c
c xsreac: reaction cross section
c
      write(*,'(/" Photoabsorption cross sections"/)')
      write(*,'("     E      reaction"/)')
      do 110 nen=ebegin(0),eend(0)
        write(*,'(1x,f8.3,es12.4)') egrid(nen),xsreac(0,nen)
  110 continue
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

      subroutine densityout(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 27, 2020
c | Task  : Output of level density
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*8      str(-1:1)
      character*12     ldfile,ldstring
      character*13     ldfileout
      character*25     model
      character*30     collstring
      integer          Zix,Nix,Z,N,A,ldmod,ibar,odd,J,ploop,parity,nex,
     +                 NL,NT,i
      real             aldmatch,SS,P,Eex,ignatyuk,spincut,ald,Krot,Kvib,
     +                 Kcoll,chi2D0,Dratio,dEx,sigma,Tnuc
      double precision densitytot,densitytotP,density,dens,Ncum,chi2,
     +                 avdev,ldtot,ldtotP
c
c ********************** Level density parameters **********************
c
c Zix        : charge number index for residual nucleus
c Nix        : neutron number index for residual nucleus
c ZZ,Z       : charge number of residual nucleus
c NN,N       : neutron number of residual nucleus
c AA,A       : mass number of residual nucleus
c S,SS       : separation energy per particle
c nuc        : symbol of nucleus
c ldmodel    : level density model
c model      : string for level density model
c flagcol    : flag for collective enhancement of level density
c ldexist    : flag for existence of level density table
c flagfission: flag for fission
c ibar       : fission barrier
c nfisbar    : number of fission barrier parameters
c alev       : level density parameter
c aldmatch   : function to determine effective level density parameter
c ldfile     : level density file
c ldstring   : string for level density file
c
      Z=ZZ(Zix,Nix,0)
      N=NN(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      SS=S(Zix,Nix,1)
      ldmod=ldmodel(Zix,Nix)
      write(*,'(/" Level density parameters for Z=",i3," N=",i3,
     +  " (",i3,a2,") "/)')  Z,N,A,nuc(Z)
      if (ldmod.eq.1) model="Gilbert-Cameron          "
      if (ldmod.eq.2) model="Back-shifted Fermi Gas   "
      if (ldmod.eq.3) model="Generalized superfluid   "
      if (ldmod.eq.4) model="Goriely tables           "
      if (ldmod.eq.5) model="Hilaire-Goriely tables   "
      if (ldmod.eq.6) model="Hilaire-Goriely Gogny    "
      write(*,'(" Model: ",a25)') model
      if (ldmod.ge.4.and..not.ldexist(Zix,Nix,0))
     +  write(*,'(" Tables not available ")')
      if (flagcol(Zix,Nix).and..not.ldexist(Zix,Nix,0)) then
        write(*,'(" Collective enhancement: yes"/)')
      else
        write(*,'(" Collective enhancement: no"/)')
      endif
      if (flagfission) then
        write(*,'(21x," g.s.     Fission barriers ")')
        write(*,'(29x,3(5x,i1,4x))') (ibar,ibar=1,nfisbar(Zix,Nix))
        write(*,'()')
      endif
      write(*,'(" a(Sn)           :",f10.5)') alev(Zix,Nix)
      if (flagfission) write(*,'(" a-effective     :",10x,3f10.5)')
     +  (aldmatch(Zix,Nix,SS,ibar),ibar=1,nfisbar(Zix,Nix))
c
c Theoretical D0
c
c D0    : experimental s-wave resonance spacing in eV
c dD0   : uncertainty in D0
c D0theo: theoretical s-wave resonance spacing
c
      write(*,'(" Experimental D0 :",f18.2," eV +- ",f15.5)')
     +  D0(Zix,Nix),dD0(Zix,Nix)
      write(*,'(" Theoretical D0  :",f18.2," eV")') D0theo(Zix,Nix)
      write(*,'(" Theoretical D1  :",f18.2," eV")') D1theo(Zix,Nix)
      write(*,'(" Av. res. energy :",f18.2," eV")') Eavres*1.e6
c
c Other parameters
c
c alimit      : asymptotic level density parameter
c gammald     : gamma-constant for asymptotic level density parameter
c pair,P      : total pairing correction
c deltaW      : shell correction in nuclear mass
c Nlast       : last discrete level
c Nlow        : lowest discrete level for temperature matching
c Ntop        : highest discrete level for temperature matching
c Exmatch     : matching point for Ex
c T           : nuclear temperature
c E0          : constant of temperature formula
c Pshift      : adjustable pairing shift
c Ucrit       : critical U
c Econd       : condensation energy
c Tcrit       : critical temperature
c scutoffdisc : spin cutoff factor for discrete level region
c spincut     : spin cutoff factor
c ignatyuk    : function for energy dependent level density parameter a
c beta2       : deformation parameter
c Krotconstant: normalization constant for rotational enhancement
c Ufermi      : energy of Fermi distribution for damping of ground-state
c               rotational effects
c cfermi      : width of Fermi distribution for damping of ground-state
c               rotational effects
c Ufermibf    : energy of Fermi distribution for damping of barrier
c               rotational effects
c cfermibf    : width of Fermi distribution for damping of barrier
c               rotational effects
c
      P=pair(Zix,Nix)
      write(*,'(" Asymptotic a    :",f10.5)') alimit(Zix,Nix)
      write(*,'(" Damping gamma   :",f10.5)') gammald(Zix,Nix)
      write(*,'(" Pairing energy  :",f10.5)') P
      write(*,'(" Shell correction:",4f10.5)') (deltaW(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'(" Last disc. level:",4(7x,i3))') (Nlast(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'(" Nlow            :",4(7x,i3))') (Nlow(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'(" Ntop            :",4(7x,i3))') (Ntop(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      if (ldmod.eq.1) then
        write(*,'(" Matching Ex     :",4f10.5)') (Exmatch(Zix,Nix,ibar),
     +    ibar=0,nfisbar(Zix,Nix))
        write(*,'(" Temperature     :",4f10.5)') (T(Zix,Nix,ibar),
     +    ibar=0,nfisbar(Zix,Nix))
        write(*,'(" E0              :",4f10.5)') (E0(Zix,Nix,ibar),
     +    ibar=0,nfisbar(Zix,Nix))
      endif
      write(*,'(" Adj. pair shift :",4f10.5)') (Pshift(Zix,Nix,ibar),
     +  ibar=0,nfisbar(Zix,Nix))
      if (ldmod.eq.3) then
        write(*,'(" Critical energy :",4f10.5)') (Ucrit(Zix,Nix,ibar),
     +    ibar=0,nfisbar(Zix,Nix))
        write(*,'(" Condensation en.:",4f10.5)') (Econd(Zix,Nix,ibar),
     +    ibar=0,nfisbar(Zix,Nix))
        write(*,'(" Critical temp.  :",f10.5)') Tcrit(Zix,Nix)
      endif
      write(*,'(" Discrete sigma  :",4f10.5)')
     +  (sqrt(scutoffdisc(Zix,Nix,ibar)),ibar=0,nfisbar(Zix,Nix))
      write(*,'(" Sigma (Sn)      :",4f10.5)')
     +  (sqrt(spincut(Zix,Nix,ignatyuk(Zix,Nix,SS,ibar),SS,ibar)),
     +  ibar=0,nfisbar(Zix,Nix))
      write(*,'(" Rhotot(Sn=",f5.2,"):",1p,4e10.3)')
     +  SS,(densitytot(Zix,Nix,SS,ibar,ldmod),
     +  ibar=0,nfisbar(Zix,Nix))
      if (flagcol(Zix,Nix).and..not.ldexist(Zix,Nix,1)) then
        write(*,'(" beta2           :",f10.5)') beta2(Zix,Nix,0)
        write(*,'(" Krotconstant    :",4f10.5)')
     +    (Krotconstant(Zix,Nix,ibar),ibar=0,nfisbar(Zix,Nix))
        write(*,'(" Ufermi          :",f10.5)') Ufermi
        write(*,'(" cfermi          :",f10.5)') cfermi
        if (ibar.gt.0) then
          write(*,'(" Ufermibf        :",f10.5)') Ufermibf
          write(*,'(" cfermibf        :",f10.5)') cfermibf
        endif
      endif
c
c ********************** Total level density ***************************
c
c odd: odd (1) or even (0) nucleus
c
      odd=mod(A,2)
      do 110 ibar=0,nfisbar(Zix,Nix)
        if (ibar.eq.0) then
          write(*,'(/" Level density per parity for ground state")')
        else
          write(*,'(/" Level density per parity for fission barrier",
     +      i3)') ibar
        endif
        write(*,'(" (Total level density also per parity)"/)')
        if (flagcol(Zix,Nix).and..not.ldexist(Zix,Nix,ibar)) then
          write(*,'("    Ex     a    sigma   total ",
     +      9("  JP= ",f4.1),"    Krot      Kvib      Kcoll")')
     +      (real(J+0.5*odd),J=0,8)
        else
          write(*,'("    Ex     a    sigma   total ",
     +      9("  JP= ",f4.1)/)') (real(J+0.5*odd),J=0,8)
        endif
c
c Tabulated level densities
c
c ploop      : help variable
c parity     : parity
c flagparity : flag for non-equal parity distribution
c edens      : energy grid for tabulated level densities
c densitytotP: total level density per parity
c nendens    : number of energies for level density grid
c ctable     : constant to adjust tabulated level densities
c ptable     : constant to adjust tabulated level densities
c
        if (ldmod.ge.4.and.ldexist(Zix,Nix,ibar)) then
          if (ldmod.eq.4) then
            ploop=1
          else
            ploop=-1
          endif
          do 120 parity=1,ploop,-2
            if (flagparity.and.parity.eq.1)
     +        write(*,'(/" Positive parity"/)')
            if (parity.eq.-1) write(*,'(/" Negative parity"/)')
            do 130 nex=1,nendens(Zix,Nix)
              Eex=edens(nex)
              write(*,'(1x,f6.2,14x,11es10.3)') Eex,
     +          densitytotP(Zix,Nix,Eex,parity,ibar,ldmod),
     +          (density(Zix,Nix,Eex,real(J+0.5*odd),parity,ibar,ldmod),
     +          J=0,8)
  130     continue
  120     continue
          write(*,'(/" Normalization:")')
          write(*,'("        ctable=",f10.5)') ctable(Zix,Nix,ibar)
          write(*,'("        ptable=",f10.5)') ptable(Zix,Nix,ibar)
        else
c
c Analytical level densities
c
c Eex       : excitation energy
c ald       : level density parameter
c aldcrit   : critical level density parameter
c colenhance: subroutine for collective enhancement
c Krot      : rotational enhancement factor
c Kvib      : vibrational enhancement factor
c Kcoll     : total collective enhancement
c pardis    : parity distribution
c densitytot: total level density
c density   : level density
c
          do 140 nex=1,nendens(Zix,Nix)
            Eex=edens(nex)
            ald=ignatyuk(Zix,Nix,Eex,ibar)
            if (ldmod.eq.3.and.Eex.lt.Ucrit(Zix,Nix,ibar)-
     +        P-Pshift(Zix,Nix,ibar)) ald=aldcrit(Zix,Nix,ibar)
            if (flagcol(Zix,Nix).and..not.ldexist(Zix,Nix,ibar)) then
              call colenhance(Zix,Nix,Eex,ald,ibar,Krot,Kvib,Kcoll)
              write(collstring,'(3es10.3)') Krot,Kvib,Kcoll
            else
              collstring=' '
            endif
            write(*,'(1x,f6.2,2f7.3,10es10.3,a30)') Eex,ald,
     +        sqrt(spincut(Zix,Nix,ald,Eex,ibar)),
     +        densitytotP(Zix,Nix,Eex,1,ibar,ldmod),
     +        (density(Zix,Nix,Eex,real(J+0.5*odd),1,ibar,ldmod),
     +        J=0,8),collstring
  140     continue
        endif
  110 continue
c
c Cumulative number of discrete levels vs. integrated level density
c
c NL,NT      : help variables
c filedensity: flag for level densities on separate files
c chi2D0     : chi-square of D0
c Dratio     : ratio D0theo/Dexp
c nucmass    : mass of nucleus
c delta0     : systematical pairing energy
c aldcrit    : critical level density parameter
c edis       : energy of level
c dEx        : energy bin for integration
c Ncum       : number of cumulative levels (integral of level density)
c C          : integration constant
c nlevmax2   : maximum number of levels
c sigma      : square root of spin cutoff factor
c chi2       : chi-square
c avdev      : average deviation
c
c Output given in general output file and on separate files.
c
      do 210 ibar=0,nfisbar(Zix,Nix)
        NL=Nlow(Zix,Nix,ibar)
        NT=Ntop(Zix,Nix,ibar)
        write(*,'(/" Discrete levels versus total level density"/)')
        write(*,'("   Energy Level   N_cumulative"/)')
        if (filedensity) then
          ldfile='ld000000.tot'
          if (ibar.gt.0) write(ldfile(10:12),'(i3.3)') ibar
          write(ldfile(3:8),'(2i3.3)') Z,A
          open (unit=1,file=ldfile,status='replace')
          if (ibar.gt.0) then
            ldstring=' Barrier    '
            write(ldstring(10:10),'(i1)') ibar
          else
            ldstring='            '
          endif
          write(1,'("# Level density for ",i3,a2,a12)') A,nuc(Z),
     +      ldstring
          if (flagcol(Zix,Nix).and..not.ldexist(Zix,Nix,ibar)) then
            write(1,'("# Model: ",a25,"Collective enhancement: yes",
     +        " beta2:",f12.5)') model,beta2(Zix,Nix,ibar)
          else
            write(1,'("# Model: ",a25,"Collective enhancement: no")')
     +        model
          endif
          write(1,'("# Experimental D0 :",f18.2," eV +- ",f15.5)')
     +      D0(Zix,Nix),dD0(Zix,Nix)
          if (dD0(Zix,Nix).eq.0.) then
            chi2D0=0.
          else
            chi2D0=(abs(D0theo(Zix,Nix)-D0(Zix,Nix))/dD0(Zix,Nix))**2
          endif
          if (D0(Zix,Nix).eq.0.) then
            Dratio=0.
          else
            Dratio=D0theo(Zix,Nix)/D0(Zix,Nix)
          endif
          write(1,'("# Theoretical D0  :",f18.2," eV Chi2:",es12.5,
     +      " Dtheo/Dexp: ",f9.5)') D0theo(Zix,Nix),chi2D0,Dratio
          write(1,'("# Nlow: ",i3," Ntop: ",i3,
     +      "           Mass in a.m.u.  : ",f10.6)') Nlow(Zix,Nix,ibar),
     +      Ntop(Zix,Nix,ibar),nucmass(Zix,Nix)
          if (ldmod.le.3) then
            write(1,'("# a(Sn)           :",f10.5,
     +        "   Asymptotic a    :",f10.5)')  alev(Zix,Nix),
     +        alimit(Zix,Nix)
            write(1,'("# Shell correction:",f10.5,
     +        "   Damping gamma   :",f10.5)') deltaW(Zix,Nix,ibar),
     +        gammald(Zix,Nix)
            write(1,'("# Pairing energy  :",f10.5,
     +        "   Separation en.  :",f10.5)') P,SS
            write(1,'("# Disc. sigma     :",f10.5,
     +        "   Sigma (Sn)      :",f10.5)')
     +        sqrt(scutoffdisc(Zix,Nix,ibar)),
     +        sqrt(spincut(Zix,Nix,ignatyuk(Zix,Nix,SS,ibar),SS,ibar))
            if (ldmod.eq.3) then
              write(1,'("# Adj. pair shift :",f10.5,
     +          "   delta0          :",f10.5,"   Crit. a:",f10.5)')
     +          Pshift(Zix,Nix,ibar),delta0(Zix,Nix),
     +          aldcrit(Zix,Nix,ibar)
            else
              write(1,'("# Adj. pair shift :",f10.5)')
     +          Pshift(Zix,Nix,ibar)
            endif
            if (ldmod.eq.1) then
              write(1,'("# Matching Ex     :",f10.5,
     +          "   Temperature     :",f10.5,"   E0:",f10.5)')
     +          Exmatch(Zix,Nix,ibar),T(Zix,Nix,ibar),E0(Zix,Nix,ibar)
            endif
            if (ldmod.eq.2) write(1,'("#")')
            if (ldmod.eq.3) then
              write(1,'("# Critical energy :",f10.5,
     +          "   Condensation en.:",f10.5,"   Crit. T:",f10.5)')
     +          Ucrit(Zix,Nix,ibar),Econd(Zix,Nix,ibar),Tcrit(Zix,Nix)
            endif
          else
            write(1,'("# ctable          :",f10.5)')
     +        ctable(Zix,Nix,ibar)
            write(1,'("# ptable          :",f10.5)')
     +        ptable(Zix,Nix,ibar)
            write(1,'("#")')
            write(1,'("#")')
            write(1,'("#")')
            write(1,'("#")')
          endif
          write(1,'("# Energy  Level N_cumulative    Total l.d. ",
     +      "       a            Sigma  ")')
        endif
        Ncum=real(NL)
        chi2=0.
        avdev=0.
        do 220 i=NL+1,nlevmax2(Zix,Nix)-1
          if (edis(Zix,Nix,i+1).eq.0.) goto 220
          Eex=0.5*(edis(Zix,Nix,i)+edis(Zix,Nix,i-1))
          dEx=edis(Zix,Nix,i)-edis(Zix,Nix,i-1)
          dens=densitytot(Zix,Nix,Eex,ibar,ldmod)
          Ncum=Ncum+dens*dEx
          if (Ncum.lt.1.e5) then
            Eex=edis(Zix,Nix,i)
            write(*,'(1x,f8.4,i4,f12.3)') Eex,i,Ncum
            if (filedensity) then
              if (ldmod.le.3) then
                ald=ignatyuk(Zix,Nix,Eex,ibar)
                if (ldmod.eq.3.and.Eex.lt.Ucrit(Zix,Nix,ibar)-P-
     +            Pshift(Zix,Nix,ibar)) ald=aldcrit(Zix,Nix,ibar)
                dens=densitytot(Zix,Nix,Eex,ibar,ldmod)
                sigma=sqrt(spincut(Zix,Nix,ald,Eex,ibar))
                write(1,'(f8.4,i4,4f14.3)') Eex,i,Ncum,dens,ald,sigma
              else
                write(1,'(f8.4,i4,2f14.3)') Eex,i,Ncum,dens
              endif
            endif
          endif
          if (i.le.NT) then
            chi2=chi2+(Ncum-i)**2/i
            avdev=avdev+abs(Ncum-i)/(NT-NL)
          endif
  220   continue
        if (filedensity) then
          write(1,'("# Chi-square per point for levels between ",
     +      "Nlow and Ntop: ",es12.5," Average deviation: ",es12.5)')
     +      chi2,avdev
          close (1)
        endif
  210 continue
c
c Level densities per parity on separate files
c
      if (filedensity) then
        ibar=0
        str(-1)='Negative'
        str(1)='Positive'
        do 310 parity=1,-1,-2
          ldfileout='nld000000.tab'
          write(ldfileout(4:9),'(2i3.3)') Z,A
          open (unit=1,status='unknown',file=ldfileout)
          write(1,'(20x,96("*"))')
          write(1,'(20x,"*  Z=",i3," A=",i3,": ",a8,"-Parity ",
     +      "Spin-dependent Level Density [MeV-1] for ",a2,i3,
     +        " and ldmodel=",i2,"  *")') Z,A,str(parity),nuc(Z),A,ldmod
          write(1,'(20x,96("*"))')
          if (mod(A,2).eq.0) then
            write(1,'(" U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT  ",
     +        31("   J=",i2.2,2x,:))') (J,J=0,29)
          else
            write(1,'(" U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT  ",
     +        31("  J=",i2.2,"/2",1x,:))') (J,J=1,59,2)
          endif
          Ncum=0.
          dEx=edens(1)
          do 320 nex=1,nendens(Zix,Nix)
            Eex=edens(nex)
            Tnuc=sqrt(Eex/alev(Zix,Nix))
            if (nex.gt.1) dEx=Eex-edens(nex-1)
            dens=densitytotP(Zix,Nix,Eex,parity,ibar,ldmod)
            Ncum=Ncum+dens*dEx
            ldtot=0.
            do 330 J=0,29
              ldtot=ldtot+(2.*J+1)*
     +          density(Zix,Nix,Eex,real(J+0.5*odd),parity,ibar,ldmod)
  330       continue
            ldtotP=densitytotP(Zix,Nix,Eex,parity,ibar,ldmod)
            write(1,'(1x,f6.2,f7.3,1x,1p,33e9.2)') Eex,Tnuc,Ncum,
     +        ldtotP,ldtot,
     +        (density(Zix,Nix,Eex,real(J+0.5*odd),1,ibar,ldmod),J=0,29)
  320     continue
          write(1,*)
  310   continue
        close(1)
      endif
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

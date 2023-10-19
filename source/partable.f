      subroutine partable(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 24, 2021
c | Task  : Write model parameters per nucleus to separate file
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,A,ibar,l
c
c ****************************** Z and A of nucleus ********************
c
c Zix : charge number index for residual nucleus
c Nix : neutron number index for residual nucleus
c ZZ,Z: charge number of residual nucleus
c AA,A: mass number of residual nucleus
c nuc : symbol of nucleus
c
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      write(51,'("##")')
      write(51,'("## Parameters for ",i3,a2)')  A,nuc(Z)
c
c ********************** Level density parameters **********************
c
c alev         : level density parameter
c aadjust....  : adjustable factors for level density parameters
c                (default 1.)
c gammald      : gamma-constant for asymptotic level density parameter
c pair         : total pairing correction
c D0           : experimental s-wave resonance spacing in eV
c ibar         : fission barrier
c nfisbar      : number of fission barrier parameters
c Pshift       : adjustable pairing shift
c deltaW       : shell correction in nuclear mass
c ldmodel      : level density model
c T            : nuclear temperature
c E0           : constant of temperature formula
c Exmatch      : matching point for Ex
c Ntop         : highest discrete level for temperature matching
c Nlow         : lowest discrete level for temperature matching
c s2adjust     : adjustable constant (Z,A,barrier-dependent) for spin
c                cutoff parameter
c flagcolall   : flag for collective enhancement of level density
c                for all nuclides
c Krotconstant : normalization constant for rotational enhancement
c ctable,ptable: constant to adjust tabulated level densities
c phmodel      : particle-hole state density model
c g            : single-particle level density parameter
c gp           : single-particle proton level density parameter
c gn           : single-particle neutron level density parameter
c
      write(51,'("##")')
      write(51,'("## Level density")')
      write(51,'("##")')
      write(51,'("a              ",2i4,f10.5)') Z,A,alev(Zix,Nix)
      write(51,'("aadjust        ",2i4,f10.5)') Z,A,aadjust(Zix,Nix)
      write(51,'("gammald        ",2i4,f10.5)') Z,A,gammald(Zix,Nix)
      write(51,'("pair           ",2i4,f10.5)') Z,A,pair(Zix,Nix)
      do 10 ibar=0,nfisbar(Zix,Nix)
        write(51,'("Pshift         ",2i4,f10.5,i4)') Z,A,
     +    Pshift(Zix,Nix,ibar),ibar
        write(51,'("Pshiftadjust   ",2i4,f10.5,i4)') Z,A,
     +    pshiftadjust(Zix,Nix,ibar),ibar
        write(51,'("deltaW         ",2i4,f10.5,i4)') Z,A,
     +    deltaW(Zix,Nix,ibar),ibar
        if (ldmodel(Zix,Nix).eq.1) then
          write(51,'("T              ",2i4,f10.5,i4)') Z,A,
     +      T(Zix,Nix,ibar),ibar
          write(51,'("E0             ",2i4,f10.5,i4)') Z,A,
     +      E0(Zix,Nix,ibar),ibar
          write(51,'("Exmatch        ",2i4,f10.5,i4)') Z,A,
     +      Exmatch(Zix,Nix,ibar),ibar
          write(51,'("Tadjust        ",2i4,f10.5,i4)') Z,A,
     +      Tadjust(Zix,Nix,ibar),ibar
          write(51,'("E0adjust       ",2i4,f10.5,i4)') Z,A,
     +      E0adjust(Zix,Nix,ibar),ibar
          write(51,'("Exmatchadjust  ",2i4,f10.5,i4)') Z,A,
     +      Exmatchadjust(Zix,Nix,ibar),ibar
        endif
        write(51,'("Ntop           ",2i4,2i4)') Z,A,Ntop(Zix,Nix,ibar),
     +   ibar
        write(51,'("Nlow           ",2i4,2i4)') Z,A,Nlow(Zix,Nix,ibar),
     +   ibar
        write(51,'("s2adjust       ",2i4,f10.5,i4)') Z,A,
     +    s2adjust(Zix,Nix,ibar),ibar
        write(51,'("ctable         ",2i4,f10.5,i4)') Z,A,
     +    ctable(Zix,Nix,ibar),ibar
        write(51,'("ptable         ",2i4,f10.5,i4)') Z,A,
     +    ptable(Zix,Nix,ibar),ibar
        write(51,'("ctableadjust   ",2i4,f10.5,i4)') Z,A,
     +    ctableadjust(Zix,Nix,ibar),ibar
        write(51,'("ptableadjust   ",2i4,f10.5,i4)') Z,A,
     +    ptableadjust(Zix,Nix,ibar),ibar
        if (flagcolall) write(51,'("Krotconstant   ",2i4,f10.5,i4)')
     +    Z,A,Krotconstant(Zix,Nix,ibar),ibar
   10 continue
      if (D0(Zix,Nix).ne.0.)
     +   write(51,'("D0             ",2i4,es12.5)') Z,A,
     +   D0(Zix,Nix)*0.001
      if (phmodel.eq.1) then
        write(51,'("g              ",2i4,f10.5)') Z,A,g(Zix,Nix)
        write(51,'("gp             ",2i4,f10.5)') Z,A,gp(Zix,Nix)
        write(51,'("gn             ",2i4,f10.5)') Z,A,gn(Zix,Nix)
        write(51,'("gnadjust       ",2i4,f10.5)') Z,A,gnadjust(Zix,Nix)
        write(51,'("gpadjust       ",2i4,f10.5)') Z,A,gpadjust(Zix,Nix)
        write(51,'("gadjust        ",2i4,f10.5)') Z,A,gadjust(Zix,Nix)
      endif
c
c ************************ Gamma-ray parameters ************************
c
c gamgam       : total radiative width in eV
c gammax       : number of l-values for gamma multipolarity
c sgr          : strength of GR
c egr          : energy of GR
c ggr          : width of GR
c egradjust....: adjustable factors for giant resonance parameters
c strength     : model for E1 gamma-ray strength function
c etable,ftable: constant to adjust tabulated strength functions
c ngr          : number of GR
c
      write(51,'("##")')
      write(51,'("## Gamma-ray")')
      write(51,'("##")')
      write(51,'("gamgam         ",2i4,f10.5)') Z,A,gamgam(Zix,Nix)
      write(51,'("gamgamadjust   ",2i4,f10.5)') Z,A,
     +  gamgamadjust(Zix,Nix)
      do 110 l=1,gammax
        if (strength.le.2.or.strength.eq.5) then
          write(51,'("sgr            ",2i4,f8.3," E",i1)') Z,A,
     +      sgr(Zix,Nix,1,l,1),l
          write(51,'("egr            ",2i4,f8.3," E",i1)') Z,A,
     +      egr(Zix,Nix,1,l,1),l
          write(51,'("ggr            ",2i4,f8.3," E",i1)') Z,A,
     +      ggr(Zix,Nix,1,l,1),l
          write(51,'("sgradjust      ",2i4,f8.3," E",i1)') Z,A,
     +      sgradjust(Zix,Nix,1,l,1),l
          write(51,'("egradjust      ",2i4,f8.3," E",i1)') Z,A,
     +      egradjust(Zix,Nix,1,l,1),l
          write(51,'("ggradjust      ",2i4,f8.3," E",i1)') Z,A,
     +      ggradjust(Zix,Nix,1,l,1),l
        else
          write(51,'("etable         ",2i4,f10.5," E",i1)') Z,A,
     +      etable(Zix,Nix,1,l),l
          write(51,'("ftable         ",2i4,f10.5," E",i1)') Z,A,
     +      ftable(Zix,Nix,1,l),l
          write(51,'("wtable         ",2i4,f10.5," E",i1)') Z,A,
     +      wtable(Zix,Nix,1,l),l
          write(51,'("etableadjust   ",2i4,f10.5," E",i1)') Z,A,
     +      etableadjust(Zix,Nix,1,l),l
          write(51,'("ftableadjust   ",2i4,f10.5," E",i1)') Z,A,
     +      ftableadjust(Zix,Nix,1,l),l
          write(51,'("wtableadjust   ",2i4,f10.5," E",i1)') Z,A,
     +      wtableadjust(Zix,Nix,1,l),l
        endif
        if (strengthM1.eq.8.or.strengthM1.eq.10) then
          write(51,'("etable         ",2i4,f10.5," M",i1)') Z,A,
     +      etable(Zix,Nix,0,l),l
          write(51,'("ftable         ",2i4,f10.5," M",i1)') Z,A,
     +      ftable(Zix,Nix,0,l),l
          write(51,'("wtable         ",2i4,f10.5," M",i1)') Z,A,
     +      wtable(Zix,Nix,0,l),l
          write(51,'("etableadjust   ",2i4,f10.5," M",i1)') Z,A,
     +      etableadjust(Zix,Nix,0,l),l
          write(51,'("ftableadjust   ",2i4,f10.5," M",i1)') Z,A,
     +      ftableadjust(Zix,Nix,0,l),l
          write(51,'("wtableadjust   ",2i4,f10.5," M",i1)') Z,A,
     +      wtableadjust(Zix,Nix,0,l),l
        endif
        if (strengthM1.ge.3) then
          write(51,'("upbendc        ",2i4,es12.5," M",i1)') Z,A,
     +      upbend(Zix,Nix,0,l,1),l
        endif
        if (ngr(Zix,Nix,1,l).eq.2) then
          write(51,'("sgr            ",2i4,f8.3," E",i1," 2")') Z,A,
     +      sgr(Zix,Nix,1,l,2),l
          write(51,'("egr            ",2i4,f8.3," E",i1," 2")') Z,A,
     +      egr(Zix,Nix,1,l,2),l
          write(51,'("ggr            ",2i4,f8.3," E",i1," 2")') Z,A,
     +      ggr(Zix,Nix,1,l,2),l
          write(51,'("sgradjust      ",2i4,f8.3," E",i1," 2")') Z,A,
     +      sgradjust(Zix,Nix,1,l,2),l
          write(51,'("egradjust      ",2i4,f8.3," E",i1," 2")') Z,A,
     +      egradjust(Zix,Nix,1,l,2),l
          write(51,'("ggradjust      ",2i4,f8.3," E",i1," 2")') Z,A,
     +      ggradjust(Zix,Nix,1,l,2),l
        endif
        write(51,'("sgr            ",2i4,f8.3," M",i1)') Z,A,
     +    sgr(Zix,Nix,0,l,1),l
        write(51,'("egr            ",2i4,f8.3," M",i1)') Z,A,
     +    egr(Zix,Nix,0,l,1),l
        write(51,'("ggr            ",2i4,f8.3," M",i1)') Z,A,
     +    ggr(Zix,Nix,0,l,1),l
        write(51,'("sgradjust      ",2i4,f8.3," M",i1)') Z,A,
     +    sgradjust(Zix,Nix,0,l,1),l
        write(51,'("egradjust      ",2i4,f8.3," M",i1)') Z,A,
     +    egradjust(Zix,Nix,0,l,1),l
        write(51,'("ggradjust      ",2i4,f8.3," M",i1)') Z,A,
     +    ggradjust(Zix,Nix,0,l,1),l
  110 continue
c
c ************************** Fission parameters ************************
c
c flagfission: flag for fission
c nfisbar    : number of fission barrier parameters
c fbarrier   : height of fission barrier
c fbaradjust : adjustable factors for fission parameters
c fwidth     : width of fission barrier
c bdamp      : fission partial damping parameter
c fismodelx  : fission model
c widthc2    : width of class2 states
c Rtransmom  : normalization constant for moment of inertia for
c              transition states
c Rclass2mom : normalization constant for moment of inertia for
c              class 2 states
c
      if (flagfission) then
        write(51,'("##")')
        write(51,'("## Fission parameters")')
        write(51,'("##")')
        do 210 ibar=1,nfisbar(Zix,Nix)
          if (fismodelx(Zix,Nix).eq.5) then
            if (ibar.eq.1) then
              write(51,'("betafiscor     ",2i4,f10.5)') Z,A,
     +          betafiscor(Zix,Nix)
              write(51,'("vfiscor        ",2i4,f10.5)') Z,A,
     +          vfiscor(Zix,Nix)
              write(51,'("betafiscoradjust ",2i4,f10.5)') Z,A,
     +          betafiscoradjust(Zix,Nix)
              write(51,'("vfiscoradjust  ",2i4,f10.5)') Z,A,
     +          vfiscoradjust(Zix,Nix)
            endif
            write(51,'("bdamp          ",2i4,f10.5,i3)') Z,A,
     +        bdamp(Zix,Nix,ibar),ibar
            write(51,'("bdampadjust    ",2i4,f10.5,i3)') Z,A,
     +        bdampadjust(Zix,Nix,ibar),ibar
          else
            write(51,'("fisbar         ",2i4,f10.5,i3)') Z,A,
     +        fbarrier(Zix,Nix,ibar),ibar
            write(51,'("fishw          ",2i4,f10.5,i3)') Z,A,
     +        fwidth(Zix,Nix,ibar),ibar
            write(51,'("fisbaradjust   ",2i4,f10.5,i3)') Z,A,
     +        fbaradjust(Zix,Nix,ibar),ibar
            write(51,'("fishwadjust    ",2i4,f10.5,i3)') Z,A,
     +        fwidthadjust(Zix,Nix,ibar),ibar
          endif
          if (ibar.lt.nfisbar(Zix,Nix))
     +      write(51,'("class2width    ",2i4,f10.5,i3)') Z,A,
     +      widthc2(Zix,Nix,ibar),ibar
          write(51,'("Rtransmom      ",2i4,f10.5,i3)') Z,A,
     +      Rtransmom(Zix,Nix,ibar),ibar
          write(51,'("Rclass2mom     ",2i4,f10.5,i3)') Z,A,
     +      Rclass2mom(Zix,Nix,ibar),ibar
  210   continue
      endif
      write(51,'("##--------------------------------------------")')
      return
      end
Copyright (C)  2019 A.J. Koning, S. Hilaire and S. Goriely

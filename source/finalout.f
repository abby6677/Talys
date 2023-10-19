        subroutine finalout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 2, 2021
c | Task  : Output of final results
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*6   discfile,contfile
      character*9   totfile,reactionstring(0:6)
      character*10  binfile
      character*11  fisfile
      character*12  rpfile,isofile,xsfile
      character*19  gamfile
      character*80  string
      character*126 stringtot
      integer       istat,nen,type,Acomp,Zcomp,Ncomp,Z,A,nex,Zix,Nix,
     +              npart,ia,ih,it,id,ip,in,ident,idc,i1,i2
      real          Egam
c
c Write model parameters to separate file
c
c flagpartable  : flag for output of model parameters on separate file
c ldmodel       : level density model
c alphad        : alpha-constant for asymptotic level density parameter
c betald        : beta-constant for asymptotic level density parameter
c gammashell1   : gamma-constant for asymptotic level density parameter
c gammashell2   : gamma-constant for asymptotic level density parameter
c pairconstant  : constant for pairing energy systematics
c Pshiftconstant: global constant for pairing shift
c Rspincut      : adjustable constant (global) for spin cutoff factor
c Rspincutff    : parameter (global) for FF spin cutoff factor
c cglobal,...   : global constant to adjust tabulated level densities
c flagcolall    : flag for collective enhancement of level density
c                 for all nuclides
c Ufermi        : energy of Fermi distribution for damping of
c               : ground-state rotational effects
c cfermi        : width of Fermi distribution for damping of
c               : ground-state rotational effects
c Ufermibf      : energy of Fermi distribution for damping of barrier
c               : rotational effects
c cfermibf      : width of Fermi distribution for damping of barrier
c               : rotational effects
c phmodel       : particle-hole state density model
c Kph           : constant for single-particle level density parameter
c                 (g=A/Kph)
c gnorm         : gamma normalization factor
c xscaptherm    : thermal capture cross section
c M2constant    : overall constant for matrix element in exciton model
c M2limit       : constant for asymptotic value for matrix element
c M2shift       : constant for energy shift for matrix element
c Rpinu,Rnupi...: ratio for two-component matrix element
c Rgamma        : adjustable parameter for pre-equilibrium gamma decay
c Esurf         : well depth for surface interaction
c parsym        : symbol of particle
c xspreeqtotps  : preequilibrium cross section per particle type for
c                 pickup and stripping
c xspreeqtotki  : preequilibrium cross section per particle type for
c                 knockout and inelastic
c xspreeqtotbu  : preequilibrium cross section per particle type for
c                 breakup
c Cstrip        : adjustable parameter for stripping/pick-up reactions
c Cknock        : adjustable parameter for knockout reactions
c Cbreak        : adjustable parameter for breakup reactions
c Cnubar1       : adjustable parameter for nubar constant value
c Cnubar2       : adjustable parameter for nubar energy slope
c Tmadjust      : adjustable parameter for PFNS temperature
c Fsadjust      : adjustable parameter for PFNS scission fraction
c flagjlm       : flag for using semi-microscopic JLM OMP
c alphaomp      : alpha optical model (1=normal, 2= McFadden-Satchler,
c                 3-5= folding potential, 6,8= Avrigeanu, 7=Nolte)
c aradialcor    : adjustable parameter for shape of DF alpha potential
c adepthcor     : adjustable parameter for depth of DF alpha potential
c v1adjust..    : adjustable factors for OMP (default 1.)
c k0            : index of incident particle
c parinclude    : logical to include outgoing particle
c flagcomp      : flag for compound nucleus calculation
c Rprime        : potential scattering radius
c Sstrength     : s,p,d,etc-wave strength function
c D0theo        : theoretical s-wave resonance spacing
c
      if (flagpartable) then
        write(51,'("##")')
        write(51,'("## General parameters")')
        write(51,'("##")')
        write(51,'("## Level density")')
        write(51,'("##")')
        if (ldmodel(0,0).le.3) then
          write(51,'("alphald        ",f10.5)') alphald(0,0)
          write(51,'("betald         ",f10.5)') betald(0,0)
          write(51,'("gammashell1    ",f10.5)') gammashell1(0,0)
          write(51,'("gammashell2    ",f10.5)') gammashell2
          write(51,'("pairconstant   ",f10.5)') pairconstant
          write(51,'("pshiftconstant ",f10.5)') Pshiftconstant(0,0)
          write(51,'("Rspincut       ",f10.5)') Rspincut
        endif
        write(51,'("Rspincutff     ",f10.5)') Rspincutff
        write(51,'("cglobal        ",es12.5)') cglobal
        write(51,'("pglobal        ",es12.5)') pglobal
        if (flagcolall) then
          write(51,'("Ufermi         ",f10.5)') Ufermi
          write(51,'("cfermi         ",f10.5)') cfermi
          write(51,'("Ufermibf       ",f10.5)') Ufermibf
          write(51,'("cfermibf       ",f10.5)') cfermibf
        endif
        if (phmodel.eq.1) write(51,'("Kph            ",f10.5)') Kph
        write(51,'("##")')
        write(51,'("## Gamma-ray")')
        write(51,'("##")')
        write(51,'("gnorm          ",f10.5)') gnorm
        write(51,'("xscaptherm     ",es12.5)') xscaptherm(-1)
        write(51,'("##")')
        write(51,'("## Pre-equilibrium")')
        write(51,'("##")')
        write(51,'("M2constant     ",f10.5)') M2constant
        write(51,'("M2limit        ",f10.5)') M2limit
        write(51,'("M2shift        ",f10.5)') M2shift
        write(51,'("Rpipi          ",f10.5)') Rpipi
        write(51,'("Rnunu          ",f10.5)') Rnunu
        write(51,'("Rpinu          ",f10.5)') Rpinu
        write(51,'("Rnupi          ",f10.5)') Rnupi
        write(51,'("Rgamma         ",f10.5)') Rgamma
        write(51,'("Esurf          ",f10.5)') Esurf
        do type=1,6
          write(51,'("Cstrip        ",a1,f10.5)') parsym(type),
     +      Cstrip(type)
          write(51,'("Cknock        ",a1,f10.5)') parsym(type),
     +      Cknock(type)
          write(51,'("Cbreak        ",a1,f10.5)') parsym(type),
     +      Cbreak(type)
        enddo
        write(51,'("##")')
        write(51,'("## Fission")')
        write(51,'("##")')
        write(51,'("Cnubar1        ",f10.5)') Cnubar1
        write(51,'("Cnubar2        ",f10.5)') Cnubar2
        write(51,'("Tmadjust       ",f10.5)') Tmadjust
        write(51,'("Fsadjust       ",f10.5)') Fsadjust
        write(51,'("##")')
        write(51,'("## Optical model")')
        write(51,'("##")')
        if (flagjlm) then
          write(51,'("lvadjust       ",f10.5)') lvadjust
          write(51,'("lwadjust       ",f10.5)') lwadjust
          write(51,'("lv1adjust      ",f10.5)') lv1adjust
          write(51,'("lw1adjust      ",f10.5)') lw1adjust
          write(51,'("lvsoadjust     ",f10.5)') lvsoadjust
          write(51,'("lwsoadjust     ",f10.5)') lwsoadjust
        endif
        if (alphaomp.ge.3.and.alphaomp.le.5) then
          write(51,'("aradialcor     ",f10.5)') aradialcor
          write(51,'("adepthcor      ",f10.5)') adepthcor
        endif
        do type=1,6
          if (.not.flagjlm.or.type.gt.2) then
            write(51,'("v1adjust      ",a1,f10.5)') parsym(type),
     +        v1adjust(type)
            write(51,'("v2adjust      ",a1,f10.5)') parsym(type),
     +        v2adjust(type)
            write(51,'("v3adjust      ",a1,f10.5)') parsym(type),
     +        v3adjust(type)
            write(51,'("v4adjust      ",a1,f10.5)') parsym(type),
     +        v4adjust(type)
            write(51,'("rvadjust      ",a1,f10.5)') parsym(type),
     +        rvadjust(type)
            write(51,'("avadjust      ",a1,f10.5)') parsym(type),
     +        avadjust(type)
            write(51,'("w1adjust      ",a1,f10.5)') parsym(type),
     +        w1adjust(type)
            write(51,'("w2adjust      ",a1,f10.5)') parsym(type),
     +        w2adjust(type)
            write(51,'("w3adjust      ",a1,f10.5)') parsym(type),
     +        w3adjust(type)
            write(51,'("w4adjust      ",a1,f10.5)') parsym(type),
     +        w4adjust(type)
            write(51,'("rwadjust      ",a1,f10.5)') parsym(type),
     +        rwadjust(type)
            write(51,'("awadjust      ",a1,f10.5)') parsym(type),
     +        awadjust(type)
            write(51,'("rvdadjust     ",a1,f10.5)') parsym(type),
     +        rvdadjust(type)
            write(51,'("avdadjust     ",a1,f10.5)') parsym(type),
     +        avdadjust(type)
            write(51,'("d1adjust      ",a1,f10.5)') parsym(type),
     +        d1adjust(type)
            write(51,'("d2adjust      ",a1,f10.5)') parsym(type),
     +        d2adjust(type)
            write(51,'("d3adjust      ",a1,f10.5)') parsym(type),
     +        d3adjust(type)
            write(51,'("rwdadjust     ",a1,f10.5)') parsym(type),
     +        rwdadjust(type)
            write(51,'("awdadjust     ",a1,f10.5)') parsym(type),
     +        awdadjust(type)
            write(51,'("rvsoadjust    ",a1,f10.5)') parsym(type),
     +        rvsoadjust(type)
            write(51,'("avsoadjust    ",a1,f10.5)') parsym(type),
     +        avsoadjust(type)
            write(51,'("vso1adjust    ",a1,f10.5)') parsym(type),
     +        vso1adjust(type)
            write(51,'("vso2adjust    ",a1,f10.5)') parsym(type),
     +        vso2adjust(type)
            write(51,'("wso1adjust    ",a1,f10.5)') parsym(type),
     +        wso1adjust(type)
            write(51,'("wso2adjust    ",a1,f10.5)') parsym(type),
     +        wso2adjust(type)
            write(51,'("rcadjust      ",a1,f10.5)') parsym(type),
     +        rcadjust(type)
          endif
        enddo
        if (k0.eq.1.and.(parinclude(0).or.flagcomp).and.
     +    Rprime.ne.0.) then
          write(51,'("##")')
          write(51,'("## Resonance parameters")')
          write(51,'("## Z   A     S0        R      xs(therm)    D0",
     +      "         a         P        Sn")')
          write(51,'("##",2i4,7es10.3)') Ztarget,Atarget,
     +      Sstrength(0)*1.e4,Rprime,xscaptherm(-1),D0theo(0,0),
     +      alev(0,0),pair(0,0),S(0,0,1)
        endif
      endif
      close (unit=51)
c
c ****************** Integrated binary cross sections ******************
c
c flagomponly: flag to execute ONLY an optical model calculation
c binfile    : file for binary output
c stringtot  : string from file
c flagexc    : flag for output of excitation functions
c numinc     : number of incident energies
c parname    : name of particle
c
      if (flagomponly) return
      if (flagnatural.or..not.flagexc) return
      write(*,'(/" ########## EXCITATION FUNCTIONS ###########"/)')
      totfile='total.tot'
      open (unit=1,file=totfile,status='old',iostat=istat)
      if (istat.eq.0) then
        write(*,'(" 1. Total (binary) cross sections"/)')
        write(*,'("   Energy    Non-elastic  Elastic     Total  ",
     +    "  Comp. el.  Shape el.  Reaction   Comp non-el",
     +    "  Direct   Pre-equil."/)')
        read(1,'(////)')
        do 10 nen=1,numinc
          read(1,'(a112)',iostat=istat) stringtot
          if (istat.ne.0) goto 10
          write(*,'(1x,a112)') stringtot
   10   continue
        close (unit=1)
      endif
      binfile='binary.tot'
      open (unit=1,file=binfile,status='old',iostat=istat)
      if (istat.eq.0) then
        write(*,'(/" 2. Binary non-elastic cross sections",
     +    " (non-exclusive)"/)')
        write(*,'("   Energy   ",7(2x,a8,1x)/)')
     +    (parname(type),type=0,6)
        read(1,'(////)')
        do 20 nen=1,numinc
          read(1,'(a112)',iostat=istat) stringtot
          if (istat.ne.0) goto 20
          write(*,'(1x,a112)') stringtot
   20   continue
        close (unit=1)
      endif
c
c ************** Total particle production cross sections **************
c
c parskip    : logical to skip outgoing particle
c flagfission: flag for fission
c
      write(*,'(/" 3. Total particle production cross sections")')
      do 110 type=0,6
        if (parskip(type)) goto 110
        totfile=' prod.tot'
        write(totfile(1:1),'(a1)') parsym(type)
        open (unit=1,file=totfile,status='old',iostat=istat)
        if (istat.eq.0) then
          write(*,'(/1x,a8, " production"/)') parname(type)
          write(*,'(" Energy   Cross section Multiplicity"/)')
          read(1,'(////)')
          do 120 nen=1,numinc
            read(1,'(a80)',iostat=istat) string
            if (istat.ne.0) goto 120
            write(*,'(1x,a80)') string
  120     continue
          close (unit=1)
        endif
  110 continue
      if (flagfission) then
        fisfile='fission.tot'
        open (unit=1,file=fisfile,status='old',iostat=istat)
        if (istat.eq.0) then
          write(*,'(/" 3b. Total fission cross sections "/)')
          write(*,'(" Energy   Cross section "/)')
          read(1,'(////)')
          do 130 nen=1,numinc
            read(1,'(a80)',iostat=istat) string
            if (istat.ne.0) goto 130
            write(*,'(1x,a80)') string
  130     continue
          close (unit=1)
        endif
      endif
c
c ******************** Residual production cross sections **************
c
c Acomp     : mass number index for compound nucleus
c maxA      : maximal number of nucleons away from the initial
c                  compound nucleus
c Zcomp     : charge number index for compound nucleus
c maxZ      : maximal number of protons away from the initial compound
c             nucleus
c Ncomp     : neutron number index for compound nucleus
c maxN      : maximal number of neutrons away from the initial
c             compound nucleus
c rpexist   : flag for existence of residual production cross section
c ZZ,Z      : charge number of residual nucleus
c AA,A      : mass number of residual nucleus
c Qres      : Q-value for residual nucleus
c Ethresh   : threshold incident energy for residual nucleus
c Nlast     : last discrete level
c rpisoexist: flag for existence of isomeric residual production cross
c             section
c nuc       : symbol of nucleus
c tau       : lifetime of state in seconds
c
      write(*,'(/" 4. Residual production cross sections")')
      do 210 Acomp=0,maxA
        do 215 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 215
          if (.not.rpexist(Zcomp,Ncomp)) goto 215
          Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
          rpfile='rp000000.tot'
          write(rpfile(3:8),'(2i3.3)') Z,A
          open (unit=1,file=rpfile,status='old',iostat=istat)
          if (istat.eq.0) then
            write(*,'(/" Production of Z=",i3,
     +        " A=",i3," (",i3,a2,") - Total"/)') Z,A,A,nuc(Z)
            write(*,'(" Q-value    =",f11.6)') Qres(Zcomp,Ncomp,0)
            write(*,'(" E-threshold=",f11.6/)') Ethresh(Zcomp,Ncomp,0)
            write(*,'(" Energy   Cross section"/)')
            read(1,'(////)')
            do 220 nen=1,numinc
              read(1,'(a80)',iostat=istat) string
              if (istat.ne.0) goto 220
              write(*,'(1x,a80)') string
  220       continue
            close (unit=1)
          endif
          do 230 nex=0,Nlast(Zcomp,Ncomp,0)
            if (.not.rpisoexist(Zcomp,Ncomp,nex)) goto 230
            isofile='rp000000.L00'
            write(isofile(3:8),'(2i3.3)') Z,A
            write(isofile(11:12),'(i2.2)') nex
            open (unit=1,file=isofile,status='old',iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Production of Z=",i3," A=",i3," (",i3,a2,
     +        ") - Level",i3,"  (lifetime:",es12.5," sec.)"/)')
     +          Z,A,A,nuc(Z),nex,tau(Zcomp,Ncomp,nex)
              write(*,'(" Q-value    =",f11.6)') Qres(Zcomp,Ncomp,nex)
              write(*,'(" E-threshold=",f11.6/)')
     +          Ethresh(Zcomp,Ncomp,nex)
              write(*,'(" Energy   Cross section Branching ratio"/)')
              read(1,'(////)')
              do 240 nen=1,numinc
                read(1,'(a80)',iostat=istat) string
                if (istat.ne.0) goto 240
                write(*,'(1x,a80)') string
  240         continue
              close (unit=1)
            endif
  230     continue
  215   continue
  210 continue
c
c ********************** Fission cross sections ************************
c
c fisexist: flag for existence of fission cross section
c
      if (flagfission) then
        write(*,'(/" 4b. Fission cross sections per fissioning",
     +    " nuclide")')
        do 310 Acomp=0,maxA
          do 315 Zcomp=0,maxZ
            Ncomp=Acomp-Zcomp
            if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 315
            if (.not.fisexist(Zcomp,Ncomp)) goto 315
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            rpfile='rp000000.fis'
            write(rpfile(3:8),'(2i3.3)') Z,A
            open (unit=1,file=rpfile,status='old',iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Fission cross section for Z=",i3,
     +          " A=",i3," (",i3,a2,")"/)') Z,A,A,nuc(Z)
              write(*,'("  Energy     Cross section"/)')
              read(1,'(////)')
              do 320 nen=1,numinc
                read(1,'(a80)',iostat=istat) string
                if (istat.ne.0) goto 320
                write(*,'(1x,a80)') string
  320         continue
              close (unit=1)
            endif
  315     continue
  310   continue
      endif
c
c ******************** Reactions to discrete states ********************
c
c flagdisc    : flag for output of discrete state cross sections
c Zindex,Zix  : charge number index for residual nucleus
c Nindex,Nix  : neutron number index for residual nucleus
c flagchannels: flag for exclusive channels calculation
c Ltarget     : excited level of target
c jdis        : spin of level
c cparity     : parity of level (character)
c parlev      : parity of level
c edis        : energy of level
c
      if (flagdisc) then
        do 410 type=0,6
          if (type.eq.k0) then
            reactionstring(type)='Inelastic'
          else
            reactionstring(type)='  ( , )  '
            write(reactionstring(type)(4:4),'(a1)') parsym(k0)
            write(reactionstring(type)(6:6),'(a1)') parsym(type)
          endif
  410   continue
        write(*,'(/" 5. Binary reactions to discrete levels",
     +    " and continuum")')
        Zcomp=0
        Ncomp=0
        do 420 type=0,6
          if (parskip(type)) goto 420
          Zix=Zindex(Zcomp,Ncomp,type)
          Nix=Nindex(Zcomp,Ncomp,type)
          Z=ZZ(Zcomp,Ncomp,type)
          A=AA(Zcomp,Ncomp,type)
          if (flagchannels) then
            contfile='  .tot'
            write(contfile(1:2),'(2a1)') parsym(k0),parsym(type)
            open (unit=1,file=contfile,status='old',iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Total exclusive ",a9," scattering")')
     +          reactionstring(type)
              write(*,'(/"  Energy     Total      Discrete   ",
     +          "Continuum     (",a1,",g",a1,")"/)') parsym(k0),
     +          parsym(type)
              read(1,'(////)')
              do 430 nen=1,numinc
                read(1,'(a80)',iostat=istat) string
                if (istat.ne.0) goto 430
                write(*,'(1x,a80)') string
  430         continue
              close (unit=1)
            endif
          endif
          do 440 nex=0,Nlast(Zix,Nix,0)
            if (type.eq.k0.and.nex.eq.Ltarget) goto 440
            discfile='  .L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
            write(discfile(5:6),'(i2.2)') nex
            open (unit=1,file=discfile,status='old',iostat=istat)
            if (istat.eq.0) then
              write(*,'(/1x,a9," to level",i3," of ",i3,a2,":",f4.1,a1,
     +          " at",f8.5," MeV"/)') reactionstring(type),nex,A,nuc(Z),
     +          jdis(Zix,Nix,nex),cparity(parlev(Zix,Nix,nex)),
     +          edis(Zix,Nix,nex)
              write(*,'("  Energy     Total      Direct",
     +          "      Compound"/)')
              read(1,'(////)')
              do 450 nen=1,numinc
                read(1,'(a80)',iostat=istat) string
                if (istat.ne.0) goto 450
                write(*,'(1x,a80)') string
  450         continue
              close (unit=1)
            endif
  440     continue
          if (flagchannels) then
            contfile='  .con'
            write(contfile(1:2),'(2a1)') parsym(k0),parsym(type)
            open (unit=1,file=contfile,status='old',iostat=istat)
            if (istat.eq.0) then
              write(*,'(/" Continuum exclusive ",a9," scattering")')
     +          reactionstring(type)
              write(*,'(/"  Energy     Total      Continuum   (",a1,
     +          ",g",a1,")"/)') parsym(k0),parsym(type)
              read(1,'(////)')
              do 460 nen=1,numinc
                read(1,'(a80)',iostat=istat) string
                if (istat.ne.0) goto 460
                write(*,'(1x,a80)') string
  460         continue
              close (unit=1)
            endif
          endif
  420   continue
      endif
c
c ******************** Exclusive channels cross sections ***************
c
c npart     : number of particles in outgoing channel
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))
c numia,....: maximal number of ejectile in channel description
c chanopen  : flag to open channel with first non-zero cross section
c idnum     : counter for exclusive channel
c idchannel : identifier for exclusive channel
c reacstring: string for exclusive reaction channel
c Qexcl     : Q-value for exclusive channel
c Ethrexc   : threshold incident energy for exclusive channel
c
c 1. Exclusive cross sections
c
      if (flagchannels) then
        write(*,'(/" 6. Exclusive cross sections")')
        do 510 npart=0,maxchannel
        do 511 ia=0,numia
        do 512 ih=0,numih
        do 513 it=0,numit
        do 514 id=0,numid
        do 515 ip=0,numip
        do 516 in=0,numin
          if (in+ip+id+it+ih+ia.ne.npart) goto 516
          if (.not.chanopen(in,ip,id,it,ih,ia)) goto 516
          ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
          do 520 idc=0,idnum
            if (idchannel(idc).eq.ident) then
              xsfile='xs000000.tot'
              write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
              open (unit=1,file=xsfile,status='old',iostat=istat)
              if (istat.eq.0) then
                write(*,'(/"     Emitted particles     reaction")')
                write(*,'("   n   p   d   t   h   a")')
                write(*,'(6i4,7x,a17)')  in,ip,id,it,ih,ia,
     +            reacstring(idc)
                write(*,'(/" Q-value    =",f11.6)') Qexcl(idc,0)
                write(*,'(" E-threshold=",f11.6/)') Ethrexcl(idc,0)
                write(*,'(" Energy   Cross section Gamma c.s. ",
     +            "c.s./res.prod.cs"/)')
                read(1,'(////)')
                do 530 nen=1,numinc
                  read(1,'(a100)',iostat=istat) stringtot
                  if (istat.ne.0) goto 530
                  write(*,'(1x,a100)') stringtot
  530           continue
                close (unit=1)
              endif
              Zcomp=ip+id+it+2*ih+2*ia
              Ncomp=in+id+2*it+ih+2*ia
              isofile='xs000000.tot'
              write(isofile(3:8),'(6i1)') in,ip,id,it,ih,ia
              do 540 nex=0,Nlast(Zcomp,Ncomp,0)
                write(isofile(11:12),'(i2.2)') nex
                open (unit=1,file=isofile,status='old',iostat=istat)
                if (istat.eq.0) then
                  write(*,'(/" Level",i3," (lifetime:",
     +              es12.5," sec.)"/)') nex,tau(Zcomp,Ncomp,nex)
                  write(*,'(" Q-value    =",f11.6)') Qexcl(idc,nex)
                  write(*,'(" E-threshold=",f11.6/)') Ethrexcl(idc,nex)
                  write(*,'(" Energy   Cross section  Branching"/)')
                  read(1,'(////)')
                  do 550 nen=1,numinc
                    read(1,'(a80)',iostat=istat) string
                    if (istat.ne.0) goto 550
                    write(*,'(1x,a80)') string
  550             continue
                  close (unit=1)
                endif
  540         continue
            endif
  520     continue
  516   continue
  515   continue
  514   continue
  513   continue
  512   continue
  511   continue
  510   continue
      endif
c
c 2. Exclusive fission cross sections
c
        if (flagfission) then
          write(*,'(/" 6b. Exclusive fission cross sections")')
          do 610 npart=0,maxchannel
          do 611 ia=0,numia
          do 612 ih=0,numih
          do 613 it=0,numit
          do 614 id=0,numid
          do 615 ip=0,numip
          do 616 in=0,numin
            if (in+ip+id+it+ih+ia.ne.npart) goto 616
            if (.not.chanopen(in,ip,id,it,ih,ia)) goto 616
            ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
            do 620 idc=0,idnum
              if (idchannel(idc).eq.ident) then
                xsfile='xs000000.fis'
                write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
                open (unit=1,file=xsfile,status='old',iostat=istat)
                if (istat.eq.0) then
                  write(*,'(/"     Emitted particles     reaction")')
                  write(*,'("   n   p   d   t   h   a")')
                  write(*,'(6i4,7x,a17)') in,ip,id,it,ih,ia,
     +              reacstring(idc)
                  write(*,'(/" Q-value    =",f11.6)') Qexcl(idc,0)
                  write(*,'(" E-threshold=",f11.6/)') Ethrexcl(idc,0)
                  write(*,'(" Energy   Cross section c.s./res.pr.cs"/)')
                  read(1,'(////)')
                  do 630 nen=1,numinc
                    read(1,'(a80)',iostat=istat) string
                    if (istat.ne.0) goto 630
                    write(*,'(1x,a80)') string
  630             continue
                  close (unit=1)
                endif
              endif
  620       continue
  616     continue
  615     continue
  614     continue
  613     continue
  612     continue
  611     continue
  610     continue
        endif
c
c ************************* Gamma-ray intensities **********************
c
c flaggamdis: flag for output of discrete gamma-ray intensities
c numlev    : maximum number of included discrete levels
c gamexist  : flag for existence of gamma production cross section
c Egam      : gamma energy
c
      if (flaggamdis) then
        write(*,'(/" 7. Gamma-ray intensities")')
        do 710 Zcomp=0,maxZ
          do 715 Ncomp=0,maxN
            Z=ZZ(Zcomp,Ncomp,0)
            A=AA(Zcomp,Ncomp,0)
            do 720 i1=0,numlev
              do 725 i2=0,i1
                if (.not.gamexist(Zcomp,Ncomp,i1,i2)) goto 725
                gamfile='gam000000L00L00.tot'
                write(gamfile(4:9),'(2i3.3)') Z,A
                write(gamfile(11:12),'(i2.2)') i1
                write(gamfile(14:15),'(i2.2)') i2
                open (unit=1,file=gamfile,status='old',iostat=istat)
                if (istat.eq.0) then
                  Egam=edis(Zcomp,Ncomp,i1)-edis(Zcomp,Ncomp,i2)
                  write(*,'(/1x,i3,a2,": Initial state",i3," (",f4.1,a1,
     +              " at",f8.4,") ---> Final state",i3," (",f4.1,a1,
     +              " at",f8.4,")  Egamma= ",f8.4/)') A,nuc(Z),i1,
     +              jdis(Zcomp,Ncomp,i1),
     +              cparity(parlev(Zcomp,Ncomp,i1)),
     +              edis(Zcomp,Ncomp,i1),i2,jdis(Zcomp,Ncomp,i2),
     +              cparity(parlev(Zcomp,Ncomp,i2)),
     +              edis(Zcomp,Ncomp,i2),Egam
                  write(*,'(" Energy   Cross section"/)')
                  read(1,'(////)')
                  do 730 nen=1,numinc
                    read(1,'(a80)',iostat=istat) string
                    if (istat.ne.0) goto 730
                    write(*,'(1x,a80)') string
  730             continue
                  close (unit=1)
                endif
  725         continue
  720       continue
  715     continue
  710   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

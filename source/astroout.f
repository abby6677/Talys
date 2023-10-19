      subroutine astroout
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely, Stephane Hilaire, Arjan Koning
c | Date  : December 12, 2016
c | Task  : Output of astrophysical reaction rates
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numZ1,numN1
      parameter    (numZ1=numZ+1,numN1=numN+1)
      logical      lexist
      character*4  machar
      character*90 mafile
      character*13 astrofile
      character*1  targetparity,yesno
      integer      Acomp,Zcomp,Ncomp,i,Z,A,ires,iresprod,iwriterp,ia,
     +             zrespro(numZ1*numN1),arespro(numZ1*numN1),nex,
     +             maxAastro
      real         rateastrorp(numT,numZ1*numN1),rateastroeps,
     +             ratio,macs,dmacs,xs,dxs,branch
c
c ************************ Output of reaction rates ********************
c
c numZ1      : numZ + 1
c numN1      : numN + 1
c Acomp      : mass number index for compound nucleus
c maxAastro  : maximal number of nucleons away from initial compound
c              nucleus for astrophysical calculations
c Zcomp      : proton number index for compound nucleus
c maxZastro  : maximal number of protons away from initial compound
c              nucleus for astrophysical calculations
c Ncomp      : neutron number index for compound nucleus
c maxNastro  : maximal number of neutrons away from initial compound
c              nucleus for astrophysical calculations
c nTmax      : effective number of temperatures for Maxwellian
c rateastro  : thermonuclear reaction rate factor
c macsastro  : thermonuclear reaction cross section
c flagastrogs: flag for calculation of astrophysics reaction rate with
c              target in ground state only
c ZZ,Z       : charge number of residual nucleus
c AA,A       : mass number of residual nucleus
c nuc        : symbol of nucleus
c T9         : Temperature grid in 10**9 K
c partf      : integrated partition function
c edis       : energy of level
c astrofile: file with astro results
c targetparity: parity of target
c branch     : branching ratio to a given excited state
c
      write(*,'(/" 8. Thermonuclear reaction rates")')
      maxAastro=maxZastro+maxNastro
      do 10 Acomp=0,maxAastro
        do 20 Zcomp=0,maxZastro
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxNastro) goto 20
          do 30 i=1,nTmax
            if (rateastro(Zcomp,Ncomp,i).gt.0.) goto 40
   30     continue
          goto 20
   40     Z=ZZ(Zcomp,Ncomp,0)
          A=AA(Zcomp,Ncomp,0)
         if (nTmax.eq.1) then
            write(*,'(/" Reaction rate for Z=",i3," A=",i3," (",i3,a2,
     +        ") at <E>=",f8.5," MeV (Excited States Contribution : ",
     +        a1,")"/)') Z,A,A,nuc(Z),astroE,yesno(.not.flagastrogs)
          else
            write(*,'(/" Reaction rate for Z=",
     +        i3," A=",i3," (",i3,a2,")"/)') Z,A,A,nuc(Z)
          endif
          write(*,'("    T        G(T)        Rate       MACS "/)')
          do 50 i=1,nTmax
            write(*,'(1x,f8.4,3es12.5)') T9(i),partf(i),
     +        rateastro(Zcomp,Ncomp,i),macsastro(Zcomp,Ncomp,i)
   50     continue
          if (flagastroex) then
            do 200 nex=0,Nlast(Zcomp,Ncomp,0)
              if (nex.gt.0.and.tau(Zcomp,Ncomp,nex).eq.0.) goto 200
              write(*,'(/" Reaction rate for Z=",
     +          i3," A=",i3," (",i3,a2,") to the excited states L",
     +          i2.2," at E=",f12.5," MeV",/)') Z,A,A,nuc(Z),nex,
     +          edis(Zcomp,Ncomp,nex)
              write(*,'("    T       Rate         MACS      ",
     +          "Branching"/)')
              do 210 i=1,nTmax
                branch=0.
                if (rateastro(Zcomp,Ncomp,i).gt.0.) branch=
     +            rateastroex(Zcomp,Ncomp,i,nex)/
     +            rateastro(Zcomp,Ncomp,i)
                write(*,'(1x,f8.4,3es12.5)') T9(i),
     +            rateastroex(Zcomp,Ncomp,i,nex),
     +            macsastroex(Zcomp,Ncomp,i,nex),branch
  210         continue
  200       continue
          endif
          if (flagracap.and.Zcomp.eq.0.and.Acomp.eq.0) then
            write(*,'(/"    T      Rate(Eq)    Rate(DC)  ",
     +        "  MACS(Eq)    MACS(DC)  "/)')
            do 60 i=1,nTmax
              write(*,'(1x,f8.4,4es12.5)') T9(i),
     +          rateastro(Zcomp,Ncomp,i)-rateastroracap(i),
     +          rateastroracap(i),
     +          macsastro(Zcomp,Ncomp,i)-macsastroracap(i),
     +          macsastroracap(i)
   60       continue
          endif
   20   continue
   10 continue
c
c Comparison with experimental MACS at 30 keV
c
c machar : part of filename for MACS
c mafile : file with MACS
c dmacs: uncertainty of MACS
c rateastrorp: rate for astro residual production
c
      if (astroE.ge.0.029.and.astroE.le.0.031) then
        Z=Ztarget
        A=Atarget
        ratio=0.
        macs=0.
        dmacs=0.
        machar='z   '
        write(machar(2:4),'(i3.3)') Z
        mafile=trim(path)//'gamma/macs/'//machar
        inquire (file=mafile,exist=lexist)
        if (lexist) then
          open (unit=2,file=mafile,status='old')
   70     read(2,'(4x,i4,2(es12.4))',end=80) ia,xs,dxs
          if (A.ne.ia) goto 70
          macs=xs
          dmacs=dxs
   80     close (unit=2)
        endif
        if (macs.gt.0.) ratio=macsastro(0,0,1)/macs
        astrofile='macs.g'
        open (unit=1,file=astrofile,status='replace')
        write(1,'("# Z   A     MACS(mb)    Exp(mb)     dEXP",
     +    "(mb)   MACS/Exp")')
        write(1,'(2i4,4es12.5)') Z,A,macsastro(0,0,1),macs,dmacs,
     +    ratio
        close (unit=1)
      endif
c
c Write results to separate files
c
c rateastroeps: cutoff value
c zresprod    : help variable
c Atarget     : mass number of target nucleus
c Ztarget     : charge number of target nucleus
c Starget     : symbol of target nucleus
c parsym      : symbol of particle
c k0          : index of incident particle
c flagfission : flag for fission
c
      astrofile='astrorate.g'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",g) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",g)")')
     +    Atarget,Starget,parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS")')
      do 110 i=1,nTmax
        write(1,'(f8.4,2es12.5)') T9(i),rateastro(0,0,i),
     +    macsastro(0,0,i)
  110 continue
      close (unit=1)
c
c output partial rates(n,g) to given excited states in a specific file
c
      if (.not.flagastroex) goto 240
      do 220 nex=1,Nlast(0,0,0)
        if (tau(0,0,nex).ne.0.) goto 230
  220 continue
      goto 240
  230 do 250 nex=0,Nlast(0,0,0)
        if (nex.eq.0.or.tau(0,0,nex).ne.0.) then
          astrofile='astrorate.g'
          write(astrofile(12:13),'(i2.2)') nex
          open (unit=1,file=astrofile,status='replace')
          if (nTmax.eq.1) then
            write(1,'("# Reaction rate for ",i3,a2,"(",a1,",g) at <E>=",
     +        f8.5," MeV (Exc. States Cont. : ",a1,
     +        ") to final level L",i2.2," at E=",f12.5," MeV")')
     +        Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs),
     +        nex,edis(0,0,nex)
          else
            write(1,'("# Reaction rate for Z=",i3," A=",i3,2x,i3,a2,"(",
     +        a1,",g) to the excited state L",i2.2," at E=",f12.5,
     +        " MeV")') Ztarget,Atarget,Atarget,Starget,parsym(k0),
     +        nex,edis(0,0,nex)
          endif
          write(1,'("#   T9      Rate         MACS      Branching")')
          do 260 i=1,nTmax
            write(1,'(1x,f8.4,3es12.5)') T9(i),
     +        rateastroex(0,0,i,nex),
     +        macsastroex(0,0,i,nex),
     +        rateastroex(0,0,i,nex)/rateastro(0,0,i)
  260     continue
          close (unit=1)
        endif
  250 continue
  240 continue
      astrofile='astrorate.p'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",p) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",p)")')
     +    Atarget,Starget,parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS")')
      do 120 i=1,nTmax
        write(1,'(f8.4,2es12.5)') T9(i),rateastro(1,0,i),
     +    macsastro(1,0,i)
  120 continue
      close (unit=1)
c
c output partial rates(p,g) to given excited states in a specific file
c
      if (.not.flagastroex) goto 340
      do 320 nex=1,Nlast(1,0,0)
        if (tau(1,0,nex).ne.0.) goto 330
  320 continue
      goto 340
  330 do 350 nex=0,Nlast(1,0,0)
        if (nex.eq.0.or.tau(1,0,nex).ne.0.) then
          astrofile='astrorate.p'
          write(astrofile(12:13),'(i2.2)') nex
          open (unit=1,file=astrofile,status='replace')
          if (nTmax.eq.1) then
            write(1,'("# Reaction rate for ",i3,a2,"(",a1,",p) at <E>=",
     +        f8.5," MeV (Exc. States Cont. : ",a1,
     +        ") to final level L",i2.2," at E=",f12.5," MeV")')
     +        Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs),
     +        nex,edis(1,0,nex)
          else
            write(1,'("# Reaction rate for Z=",i3," A=",i3,2x,i3,a2,"(",
     +        a1,",p) to the excited state L",i2.2," at E=",f12.5,
     +        " MeV")') Ztarget,Atarget,Atarget,Starget,parsym(k0),
     +        nex,edis(1,0,nex)
          endif
          write(1,'("#   T       Rate         MACS      Branching")')
          do 360 i=1,nTmax
            write(1,'(1x,f8.4,3es12.5)') T9(i),
     +        rateastroex(1,0,i,nex),
     +        macsastroex(1,0,i,nex),
     +        rateastroex(1,0,i,nex)/rateastro(1,0,i)
  360     continue
          close (unit=1)
        endif
  350 continue
  340 continue
      astrofile='astrorate.a'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a) at <E>=",
     +    f8.5," MeV (Excited States Contribution : ",a1,")")')
     +    Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs)
      else
        write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a)")')
     +  Atarget,Starget,parsym(k0)
      endif
      write(1,'("#    T       Rate       MACS")')
      do 130 i=1,nTmax
        write(1,'(f8.4,2es12.5)') T9(i),rateastro(2,2,i),
     +    macsastro(2,2,i)
  130 continue
      close (unit=1)
c
c output partial rates(p,g) to given excited states in a specific file
c
c ires   : counter
c iresprod   : counter
c iwriterp   : counter
c zrespro: Z of residual product
c arespro: A of residual product
c
      if (.not.flagastroex) goto 440
      do 420 nex=1,Nlast(2,2,0)
        if (tau(2,2,nex).ne.0.) goto 430
  420 continue
      goto 440
  430 do 450 nex=0,Nlast(2,2,0)
        if (nex.eq.0.or.tau(2,2,nex).ne.0.) then
          astrofile='astrorate.a'
          write(astrofile(12:13),'(i2.2)') nex
          open (unit=1,file=astrofile,status='replace')
          if (nTmax.eq.1) then
            write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a) at <E>=",
     +        f8.5," MeV (Exc. States Cont. : ",a1,
     +        ") to final level L",i2.2," at E=",f12.5," MeV")')
     +        Atarget,Starget,parsym(k0),astroE,yesno(.not.flagastrogs),
     +        nex,edis(2,2,nex)
          else
            write(1,'("# Reaction rate for Z=",i3," A=",i3,2x,i3,a2,"(",
     +        a1,",a) to the excited state L",i2.2," at E=",f12.5,
     +        " MeV")') Ztarget,Atarget,Atarget,Starget,parsym(k0),
     +        nex,edis(2,2,nex)
          endif
          write(1,'("#   T       Rate         MACS      Branching")')
          do 460 i=1,nTmax
            write(1,'(1x,f8.4,3es12.5)') T9(i),
     +        rateastroex(2,2,i,nex),
     +        macsastroex(2,2,i,nex),
     +        rateastroex(2,2,i,nex)/rateastro(2,2,i)
  460     continue
          close (unit=1)
        endif
  450 continue
  440 continue
      if (flagfission) then
        astrofile='astrorate.f'
        open (unit=1,file=astrofile,status='replace')
        if (nTmax.eq.1) then
          write(1,'("# Reaction rate for ",i3,a2,"(",a1,",f) at <E>=",
     +      f8.5," MeV (Excited States Contribution : ",a1,")")')
     +      Atarget,Starget,parsym(k0),astroE,
     +      yesno(.not.flagastrogs)
        else
          write(1,'("# Reaction rate for ",i3,a2,"(",a1,",a)")')
     +    Atarget,Starget,parsym(k0)
        endif
        write(1,'("#    T       Rate       MACS")')
        do 140 i=1,nTmax
          write(1,'(f8.4,2es12.5)') T9(i),rateastrofis(i),
     +      macsastrofis(i)
  140   continue
        close (unit=1)
      endif
      rateastroeps=1.e-10
      iresprod=0
      do 150 Acomp=0,maxAastro
        do 150 Zcomp=0,maxZastro
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxNastro) goto 150
          if (Zcomp.eq.0.and.Ncomp.eq.0) goto 150
          if (Zcomp.eq.0.and.Ncomp.eq.1) goto 150
          if (Zcomp.eq.1.and.Ncomp.eq.0) goto 150
          if (Zcomp.eq.2.and.Ncomp.eq.2) goto 150
          iwriterp=0
          do 160 i=1,nTmax
            if (rateastro(Zcomp,Ncomp,i).gt.rateastroeps) iwriterp=1
  160     continue
          if (iwriterp.eq.1) then
            iresprod=iresprod+1
            do 170 i=1,nTmax
              rateastrorp(i,iresprod)=rateastro(Zcomp,Ncomp,i)
  170       continue
            zrespro(iresprod)=ZZ(Zcomp,Ncomp,0)
            arespro(iresprod)=AA(Zcomp,Ncomp,0)
          endif
  150 continue
      astrofile='astrorate.tot'
      targetparity='+'
      if (targetP.eq.-1) targetparity='-'
      open (unit=1,file=astrofile,status='replace')
      if (nTmax.eq.1) then
        write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +    " reactions  Jp(GS)=",f5.1,a1," at <E>=",f8.5,
     +    " MeV (Excited States Contribution : ",a1,")")')
     +    Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +    targetparity,astroE,yesno(.not.flagastrogs)
      else
        if (.not.flagastrogs) then
          if (nonthermlev.eq.-1) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " Tot: fully thermalized target")')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity
          endif
          if (nonthermlev.gt.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " GS: thermalized target on the GS",
     +        " excluding Level nb=",i2)')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity,nonthermlev
          endif
          if (nonthermlev.eq.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " Isom: thermalized target in Level nb=",i2,
     +        " and excluding Level nb=",i2)')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity,Ltarget,nonthermlev
          endif
        else
          if (Ltarget.eq.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " GS: non-thermalized target in its ground state")')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity
          endif
          if (Ltarget.gt.0) then
            write(1,'("# Reaction rate for ",2i4,a2,"+",a1,3x,i4,
     +        " reactions  Jp(GS)=",f5.1,a1,
     +        " Isom: non-thermalized target in its excited state = L",
     +        i2.2)')
     +        Ztarget,Atarget,Starget,parsym(k0),iresprod+5,targetspin,
     +        targetparity,Ltarget
          endif
        endif
      endif
      write(1,'("     T9      G(T)    (",a1,",g)",i3,a2,"  (",a1,
     +  ",n)",i3,a2,"  (",a1,",p)",i3,a2,"  (",a1,",a)",i3,a2,1x,
     +  "   Fission  ",
     +  200(6x,i3,a2,1x))')
     +  parsym(k0),AA(0,0,0),nuc(ZZ(0,0,0)),
     +  parsym(k0),AA(0,1,0),nuc(ZZ(0,1,0)),
     +  parsym(k0),AA(1,0,0),nuc(ZZ(1,0,0)),
     +  parsym(k0),AA(2,2,0),nuc(ZZ(2,2,0)),
     +  (arespro(ires),nuc(zrespro(ires)),ires=1,iresprod)
      do 180 i=1,nTmax
        write(1,'(f8.4,200es12.5)') T9(i),partf(i),rateastro(0,0,i),
     +    rateastro(0,1,i),rateastro(1,0,i),rateastro(2,2,i),
     +    rateastrofis(i),(rateastrorp(i,ires),ires=1,iresprod)
  180 continue
      close (unit=1)
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

      subroutine angleout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 21, 2020
c | Task  : Output of discrete angular distributions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*21 discfile,legfile
      integer      LL,iang,Zix,Nix,i,type
c
c **************** Elastic scattering angular distribution *************
c
c 1. Legendre coefficients
c
c flaglegendre: flag for output of Legendre coefficients
c J2end       : 2 * end of J summation
c tleg        : total Legendre coefficient
c k0          : index of incident particle
c LL          : counter for l value
c Ltarget     : excited level of target
c dleg        : direct reaction Legendre coefficient
c cleg        : compound nucleus Legendre coefficient
c tlegnor     : total Legendre coefficient normalized to 1
c fileelastic : flag for elastic angular distribution on separate file
c parsym      : symbol of particle
c Atarget     : mass number of target nucleus
c Ztarget     : charge number of target nucleus
c Starget     : symbol of target nucleus
c Einc        : incident energy in MeV
c cleg0       : Legendre coefficient normalized to the first one
c
      write(*,'(/" 8. Discrete state angular distributions")')
      if (flaglegendre) then
        write(*,'(/" 8a1. Legendre coefficients for elastic",
     +    " scattering"/)')
        write(*,'("   L       Total           Direct",
     +    "         Compound       Normalized"/)')
        do 10 LL=0,J2end
          write(*,'(1x,i3,4es16.5)') LL,tleg(k0,Ltarget,LL),
     +      dleg(k0,Ltarget,LL),cleg(k0,Ltarget,LL),
     +      tlegnor(k0,Ltarget,LL)
   10   continue
c
c Write results to separate file
c
        if (fileelastic) then
          if (flagblock) then
            legfile='  leg.L00'
            write(legfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(legfile(8:9),'(i2.2)') Ltarget
            if (.not.legexist(k0,k0,Ltarget)) then
              legexist(k0,k0,Ltarget)=.true.
              open (unit=1,file=legfile,status='unknown')
            else
              open (unit=1,file=legfile,status='unknown',
     +          position='append')
            endif
          else
            legfile='          leg.L00'
            write(legfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(legfile(3:10),'(f8.3)') Einc
            write(legfile(3:6),'(i4.4)') int(Einc)
            write(legfile(16:17),'(i2.2)') Ltarget
            open (unit=1,file=legfile,status='unknown')
          endif
          write(1,'("# ",a1," + ",i3,a2,
     +      " Elastic scattering Legendre coefficients")') parsym(k0),
     +      Atarget,Starget
          write(1,'("# E-incident = ",f10.5)') Einc
          write(1,'("# ")')
          write(1,'("# # coeff.   =",i4)') J2end+1
          write(1,'("#  L       Total           Direct",
     +      "        Compound       Normalized    ENDF-6")')
          do 20 LL=0,J2end
            write(1,'(i3,5es16.5)') LL,tleg(k0,Ltarget,LL),
     +        dleg(k0,Ltarget,LL),cleg(k0,Ltarget,LL),
     +        tlegnor(k0,Ltarget,LL),cleg0(k0,Ltarget,LL)
   20     continue
          close (unit=1)
        endif
      endif
c
c 2. Angular distributions
c
c nangle   : number of angles
c angle    : angle
c discad   : discrete state angular distribution
c directad : direct angular distribution
c compad   : compound angular distribution
c natstring: string extension for file names
c iso      : counter for isotope
c ruth     : elastic/Rutherford ratio
c
      write(*,'(/" 8a2. Elastic scattering angular distribution"/)')
      if (k0.eq.1) then
        write(*,'(" Angle        Total          Direct",
     +    "         Compound"/)')
        do 30 iang=0,nangle
          write(*,'(1x,f5.1,3es16.5)') angle(iang),
     +      discad(k0,Ltarget,iang),directad(k0,Ltarget,iang),
     +      compad(k0,Ltarget,iang)
   30   continue
c
c Write results to separate file
c
        if (fileelastic) then
          if (flagblock) then
            discfile='nnang.L00'
            write(discfile(8:9),'(i2.2)') Ltarget
            if (.not.angexist(k0,k0,Ltarget)) then
              angexist(k0,k0,Ltarget)=.true.
              open (unit=1,file=discfile,status='unknown')
            else
              open (unit=1,file=discfile,status='unknown',
     +          position='append')
            endif
          else
            discfile='nn        ang.L00'
            write(discfile(3:10),'(f8.3)') Einc
            write(discfile(3:6),'(i4.4)') int(Einc)
            write(discfile(16:17),'(i2.2)') Ltarget
            open (unit=1,file=discfile,status='unknown')
          endif
          write(1,'("# ",a1," + ",i3,a2,
     +      " Elastic scattering angular distribution")') parsym(k0),
     +      Atarget,Starget
          write(1,'("# E-incident = ",f10.5)') Einc
          write(1,'("# ")')
          write(1,'("# # angles   =",i4)') nangle+1
          write(1,'("# Angle       xs            Direct",
     +      "         Compound")')
          do 40 iang=0,nangle
            write(1,'(f5.1,3es16.5)') angle(iang),
     +        discad(k0,Ltarget,iang),directad(k0,Ltarget,iang),
     +        compad(k0,Ltarget,iang)
   40     continue
          close (unit=1)
        endif
      else
        write(*,'(" Angle        Total            Direct",
     +    "       Compound       c.s/Rutherford"/)')
        do 50 iang=0,nangle
          write(*,'(1x,f5.1,4es16.5)') angle(iang),
     +      max(discad(k0,Ltarget,iang),directad(k0,Ltarget,iang)),
     +      directad(k0,Ltarget,iang),compad(k0,Ltarget,iang),ruth(iang)
   50   continue
c
c Write results to separate file
c
        if (fileelastic) then
          if (flagblock) then
            discfile='  ang.L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(discfile(8:9),'(i2.2)') Ltarget
            if (.not.angexist(k0,k0,Ltarget)) then
              angexist(k0,k0,Ltarget)=.true.
              open (unit=1,file=discfile,status='unknown')
            else
              open (unit=1,file=discfile,status='unknown',
     +          position='append')
            endif
          else
            discfile='          ang.L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(discfile(3:10),'(f8.3)') Einc
            write(discfile(3:6),'(i4.4)') int(Einc)
            write(discfile(16:17),'(i2.2)') Ltarget
            open (unit=1,file=discfile,status='unknown')
          endif
          write(1,'("# ",a1," + ",i3,a2,
     +      " Elastic scattering angular distribution")') parsym(k0),
     +      Atarget,Starget
          write(1,'("# E-incident = ",f10.5)') Einc
          write(1,'("# ")')
          write(1,'("# # angles   =",i4)') nangle+1
          write(1,'("# Angle       xs            Direct",
     +      "         Compound    c.s./Rutherford")')
          do 60 iang=0,nangle
            write(1,'(f5.1,4es16.5)') angle(iang),
     +        max(discad(k0,Ltarget,iang),directad(k0,Ltarget,iang)),
     +        directad(k0,Ltarget,iang),compad(k0,Ltarget,iang),
     +        ruth(iang)
   60     continue
          close (unit=1)
        endif
      endif
c
c ************** Inelastic scattering angular distributions ************
c
c Zindex,Zix: charge number index for residual nucleus
c Nindex,Nix: neutron number index for residual nucleus
c nlev      : number of levels for nucleus
c xsdisc    : total cross section for discrete state
c fileangle : designator for angular distributions on separate file
c
      Zix=Zindex(0,0,k0)
      Nix=Nindex(0,0,k0)
c
c 1. Legendre coefficients
c
      if (flaglegendre) then
        write(*,'(/" 8b1. Legendre coefficients for inelastic",
     +    " scattering")')
        do 110 i=0,nlev(Zix,Nix)
          if (i.eq.Ltarget) goto 110
          if (xsdisc(k0,i).eq.0.) goto 110
          write(*,'(/"    Level ",i2/)') i
          write(*,'("   L       Total           Direct",
     +      "         Compound       Normalized"/)')
          do 120 LL=0,J2end
            write(*,'(1x,i3,4es16.5)') LL,tleg(k0,i,LL),dleg(k0,i,LL),
     +        cleg(k0,i,LL),tlegnor(k0,i,LL)
  120     continue
c
c Write results to separate file
c
          if (fileangle(i)) then
            if (flagblock) then
              legfile='  leg.L00'
              write(legfile(1:2),'(2a1)') parsym(k0),parsym(k0)
              write(legfile(8:9),'(i2.2)') i
              if (.not.legexist(k0,k0,i)) then
                legexist(k0,k0,i)=.true.
                open (unit=1,file=legfile,status='unknown')
              else
                open (unit=1,file=legfile,status='unknown',
     +            position='append')
              endif
            else
              legfile='          leg.L00'
              write(legfile(1:2),'(2a1)') parsym(k0),parsym(k0)
              write(legfile(3:10),'(f8.3)') Einc
              write(legfile(3:6),'(i4.4)') int(Einc)
              write(legfile(16:17),'(i2.2)') i
              open (unit=1,file=legfile,status='unknown')
            endif
            write(1,'("# ",a1," + ",i3,a2,
     +        " Inelastic scattering Legendre coefficients"," - Level",
     +        i3)')  parsym(k0),Atarget,Starget,i
            write(1,'("# E-incident = ",f10.5)') Einc
            write(1,'("# ")')
            write(1,'("# # coeff.   =",i4)') J2end+1
            write(1,'("#  L       Total           Direct",
     +        "        Compound       Normalized    ENDF-6")')
            do 130 LL=0,J2end
              write(1,'(i3,5es16.5)') LL,tleg(k0,i,LL),dleg(k0,i,LL),
     +          cleg(k0,i,LL),tlegnor(k0,i,LL),cleg0(k0,i,LL)
  130       continue
            close (unit=1)
          endif
  110   continue
      endif
c
c 2. Angular distributions
c
      write(*,'(/" 8b2. Inelastic angular distributions")')
      do 140 i=0,nlev(Zix,Nix)
        if (i.eq.Ltarget) goto 140
        if (xsdisc(k0,i).eq.0.) goto 140
        write(*,'(/"    Level ",i2/)') i
        write(*,'(" Angle       Total         Direct       Compound"/)')
        do 150 iang=0,nangle
          write(*,'(1x,f5.1,3es15.5)') angle(iang),discad(k0,i,iang),
     +      directad(k0,i,iang),compad(k0,i,iang)
  150   continue
c
c Write results to separate file
c
        if (fileangle(i)) then
          if (flagblock) then
            discfile='  ang.L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(discfile(8:9),'(i2.2)') i
            if (.not.angexist(k0,k0,i)) then
              angexist(k0,k0,i)=.true.
              open (unit=1,file=discfile,status='unknown')
            else
              open (unit=1,file=discfile,status='unknown',
     +          position='append')
            endif
          else
            discfile='          ang.L00'
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(discfile(3:10),'(f8.3)') Einc
            write(discfile(3:6),'(i4.4)') int(Einc)
            write(discfile(16:17),'(i2.2)') i
            open (unit=1,file=discfile,status='unknown')
          endif
          write(1,'("# ",a1," + ",i3,a2,
     +      " Inelastic scattering angular distribution",
     +      " - Level",i3)') parsym(k0),Atarget,Starget,i
          write(1,'("# E-incident = ",f10.5)') Einc
          write(1,'("# ")')
          write(1,'("# # angles   =",i4)') nangle+1
          write(1,'("# Angle      xs           Direct       Compound")')
          do 160 iang=0,nangle
            write(1,'(f5.1,3es15.5)') angle(iang),discad(k0,i,iang),
     +        directad(k0,i,iang),compad(k0,i,iang)
  160     continue
          close (unit=1)
        endif
  140 continue
c
c ********** Non-inelastic scattering angular distributions ************
c
c parskip  : logical to skip outgoing particle
c xsdisctot: total cross section summed over discrete states
c
      write(*,'(/" 8c. Angular distributions for other reactions")')
      do 210 type=0,6
        if (parskip(type)) goto 210
        if (type.eq.k0) goto 210
        if (xsdisctot(type).eq.0.) goto 210
        Zix=Zindex(0,0,type)
        Nix=Nindex(0,0,type)
c
c 1. Legendre coefficients
c
        if (flaglegendre) then
          write(*,'(/" 8c1. Legendre coefficients for (",
     +      a1,",",a1,")")') parsym(k0),parsym(type)
          do 220 i=0,nlev(Zix,Nix)
            if (xsdisc(type,i).eq.0.) goto 220
            write(*,'(/"    Level ",i2/)') i
            write(*,'("   L       Total           Direct         ",
     +        "Compound       Normalized"/)')
            do 230 LL=0,J2end
              write(*,'(1x,i3,4es16.5)') LL,tleg(type,i,LL),
     +          dleg(type,i,LL),cleg(type,i,LL),tlegnor(type,i,LL)
  230       continue
c
c Write results to separate file
c
            if (fileangle(i)) then
              if (flagblock) then
                legfile='  leg.L00'
                write(legfile(1:2),'(2a1)') parsym(k0),parsym(type)
                write(legfile(8:9),'(i2.2)') i
                if (.not.legexist(k0,type,i)) then
                  legexist(k0,type,i)=.true.
                  open (unit=1,file=legfile,status='unknown')
                else
                  open (unit=1,file=legfile,status='unknown',
     +              position='append')
                endif
              else
                legfile='          leg.L00'
                write(legfile(1:2),'(2a1)') parsym(k0),parsym(type)
                write(legfile(3:10),'(f8.3)') Einc
                write(legfile(3:6),'(i4.4)') int(Einc)
                write(legfile(16:17),'(i2.2)') i
                open (unit=1,file=legfile,status='unknown')
              endif
              write(1,'("# ",a1," + ",i3,a2," (",a1,",",a1,
     +          ") Legendre coefficients"," - Level",i3)') parsym(k0),
     +          Atarget,Starget,parsym(k0),parsym(type),i
              write(1,'("# E-incident = ",f10.5)') Einc
              write(1,'("# ")')
              write(1,'("# # coeff.   =",i4)') J2end+1
              write(1,'("#  L       Total           Direct",
     +          "        Compound       Normalized    ENDF-6")')
              do 240 LL=0,J2end
                write(1,'(i3,5es16.5)') LL,tleg(type,i,LL),
     +            dleg(type,i,LL),cleg(type,i,LL),tlegnor(type,i,LL),
     +            cleg0(type,i,LL)
  240       continue
            close (unit=1)
          endif
  220     continue
        endif
c
c 2. Angular distributions
c
        write(*,'(/" 8c2. (",a1,",",a1,") angular distributions")')
     +    parsym(k0),parsym(type)
        do 250 i=0,nlev(Zix,Nix)
          if (xsdisc(type,i).eq.0.) goto 250
          write(*,'(/"    Level ",i2/)') i
          write(*,'(" Angle      Total          Direct",
     +      "        Compound"/)')
          do 260 iang=0,nangle
            write(*,'(1x,f5.1,3es15.5)') angle(iang),
     +        discad(type,i,iang),directad(type,i,iang),
     +        compad(type,i,iang)
  260     continue
c
c Write results to separate file
c
          if (fileangle(i)) then
            if (flagblock) then
              discfile='  ang.L00'
              write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
              write(discfile(8:9),'(i2.2)') i
              if (.not.angexist(k0,type,i)) then
                angexist(k0,type,i)=.true.
                open (unit=1,file=discfile,status='unknown')
              else
                open (unit=1,file=discfile,status='unknown',
     +            position='append')
              endif
            else
              discfile='          ang.L00'
              write(discfile(1:2),'(2a1)') parsym(k0),parsym(type)
              write(discfile(3:10),'(f8.3)') Einc
              write(discfile(3:6),'(i4.4)') int(Einc)
              write(discfile(16:17),'(i2.2)') i
              open (unit=1,file=discfile,status='unknown')
            endif
            write(1,'("# ",a1," + ",i3,a2," (",a1,",",a1,
     +        ") angular distributions - Level",i3)') parsym(k0),
     +        Atarget,Starget,parsym(k0),parsym(type),i
            write(1,'("# E-incident = ",f10.5)') Einc
            write(1,'("# ")')
            write(1,'("# # angles   =",i4)') nangle+1
            write(1,'("# Angle      xs           Direct",
     +        "        Compound")')
            do 270 iang=0,nangle
              write(1,'(f5.1,3es15.5)') angle(iang),
     +          discad(type,i,iang),directad(type,i,iang),
     +          compad(type,i,iang)
  270       continue
            close (unit=1)
          endif
  250   continue
  210 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

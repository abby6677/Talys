      subroutine spectraout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 21, 2020
c | Task  : Output of particle spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*10 Efile
      character*21 specfile
      integer      type,nen
c
c ***************************** Spectra ********************************
c
c parskip    : logical to skip outgoing particle
c xsparticle : total particle production cross section
c parname    : name of particle
c ebegin     : first energy point of energy grid
c eend       : last energy point of energy grid
c espec      : outgoing energy grid
c xssumout   : cross section summed over mechanisms
c xsdiscout  : total smoothed cross section for discrete state
c xspreeqout : preequilibrium cross section per particle type and
c              outgoing energy
c xsmpreeqout: multiple pre-equilibrium emission spectrum
c xscompout  : compound emission cross section
c flagrecoil : flag for calculation of recoils
c flaglabddx : flag for calculation of DDX in LAB system
c iejlab     : number of ejectile lab bins
c Eejlab     : center of ejectile lab bin
c xsejlab    : LAB ejectile spectrum
c xsejlabint : LAB energy-integrated spectrum
c
      write(*,'(/" 7. Composite particle spectra")')
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (xsparticle(type).eq.0.) goto 10
        write(*,'(/" Spectra for outgoing ",a8/)') parname(type)
        if (k0.le.2.and.type.le.2) then
          write(*,'("  Energy   Total       Direct    Pre-equil.",
     +      "  Mult. preeq  Compound"/)')
          do 20 nen=ebegin(type),eendout(type)
            write(*,'(f8.3,5es12.5)') espec(type,nen),
     +        xssumout(type,nen),xsdiscout(type,nen),
     +        xspreeqout(type,nen),xsmpreeqout(type,nen),
     +        xscompout(type,nen)
   20     continue
        else
          write(*,'("  Energy   Total       Direct    Pre-equil.",
     +      "  Mult. preeq  Compound    Stripping   Knock-out",
     +      "   Break-up"/)')
          do 25 nen=ebegin(type),eendout(type)
            write(*,'(f8.3,8es12.5)') espec(type,nen),
     +        xssumout(type,nen),xsdiscout(type,nen),
     +        xspreeqout(type,nen),xsmpreeqout(type,nen),
     +        xscompout(type,nen),xspreeqpsout(type,nen),
     +        xspreeqkiout(type,nen),xspreeqbuout(type,nen)
   25     continue
        endif
        if (flagrecoil.and.flaglabddx) then
          write(*,'(/" LAB spectra for outgoing ",a8/)') parname(type)
          write(*,'("  Energy   Cross section"/)')
          do 30 nen=1,iejlab(type)
            write(*,'(f8.3,es12.5)') Eejlab(type,nen),
     +        xsejlab(type,nen)
   30     continue
          write(*,'(/" Energy-integrated cross section:",es12.5/)')
     +      xsejlabint(type)
        endif
c
c Write results to separate file
c
c filespectrum: designator for spectrum on separate file
c natstring   : string extension for file names
c iso         : counter for isotope
c Einc        : incident energy in MeV
c specfile    : file with spectrum
c parsym      : symbol of particle
c preeqratio  : pre-equilibrium ratio
c
        if (filespectrum(type)) then
          if (flagblock) then
            specfile=' spec.tot'//natstring(iso)
            write(specfile(1:1),'(a1)') parsym(type)
            if (.not.spexist1(type)) then
              spexist1(type)=.true.
              open (unit=1,file=specfile,status='unknown')
            else
              open (unit=1,file=specfile,status='unknown',
     +          position='append')
            endif
          else
            specfile=' spec0000.000.tot'//natstring(iso)
            write(specfile(1:1),'(a1)') parsym(type)
            write(specfile(6:13),'(f8.3)') Einc
            write(specfile(6:9),'(i4.4)') int(Einc)
            open (unit=1,file=specfile,status='unknown')
          endif
          write(1,'("# ",a1," + ",i3,a2,": ",a8," spectrum")')
     +      parsym(k0),Atarget,Starget,parname(type)
          write(1,'("# E-incident = ",f10.5)') Einc
          write(1,'("# E-average  = ",f8.3)') Eaverage(type)
          write(1,'("# # energies =",i6)') eendout(type)-ebegin(type)+1
          if (k0.le.2.and.type.le.2) then
            write(1,'("# E-out    Total       Direct    Pre-equil.",
     +        "  Mult. preeq  Compound   PE ratio   ")')
            do 40 nen=ebegin(type),eendout(type)
              write(1,'(f8.3,6es12.5)')
     +          espec(type,nen),xssumout(type,nen),xsdiscout(type,nen),
     +          xspreeqout(type,nen),xsmpreeqout(type,nen),
     +          xscompout(type,nen),preeqratio(type,nen)
   40       continue
          else
            write(1,'("# E-out    Total       Direct    Pre-equil.",
     +        "  Mult. preeq  Compound    PE ratio   BU ratio   ",
     +        " Stripping   Knock-out   Break-up")')
            do 45 nen=ebegin(type),eendout(type)
              write(1,'(f8.3,10es12.5)')
     +          espec(type,nen),xssumout(type,nen),xsdiscout(type,nen),
     +          xspreeqout(type,nen),xsmpreeqout(type,nen),
     +          xscompout(type,nen),preeqratio(type,nen),
     +          buratio(type,nen),xspreeqpsout(type,nen),
     +          xspreeqkiout(type,nen),xspreeqbuout(type,nen)
   45       continue
          endif
          close (unit=1)
          if (flagrecoil.and.flaglabddx) then
            if (flagblock) then
              specfile=' spec.lab'//natstring(iso)
              write(specfile(1:1),'(a1)') parsym(type)
              if (.not.spexist2(type)) then
                spexist2(type)=.true.
                open (unit=1,file=specfile,status='unknown')
              else
                open (unit=1,file=specfile,status='unknown',
     +            position='append')
              endif
            else
              specfile=' spec0000.000.lab'//natstring(iso)
              write(specfile(1:1),'(a1)') parsym(type)
              write(specfile(6:13),'(f8.3)') Einc
              write(specfile(6:9),'(i4.4)') int(Einc)
              open (unit=1,file=specfile,status='unknown')
            endif
            write(1,'("# ",a1," + ",i3,a2,": ",a8,
     +        " spectrum in LAB frame")') parsym(k0),Atarget,
     +        Starget,parname(type)
            write(1,'("# E-incident = ",f10.5)') Einc
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') iejlab(type)
            write(1,'("# E-out    Total")')
            do 50 nen=1,iejlab(type)
              write(1,'(f8.3,es12.5)') Eejlab(type,nen),
     +          xsejlab(type,nen)
   50       continue
            close (unit=1)
          endif
        endif
   10 continue
      write(*,'(/" Average emission energies"/)')
      do 110 type=0,6
        if (parskip(type)) goto 110
        write(*,'(1x,a8,4x,f8.3)') parname(type),Eaverage(type)
        if (filespectrum(type)) then
          Efile='Eaverage.'//parsym(type)
          if (nin.eq.numinclow+1) then
            open (unit=1,file=Efile,status='replace')
            write(1,'("# ",a1," + ",i3,a2," Average ",a8,
     +        " emission energy")') parsym(k0),Atarget,Starget,
     +        parname(type)
            write(1,'("# Q-value    =",es12.5)') Q(type)
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') numinc
            write(1,'("#    E       E-average")')
            do 210 nen=1,numinclow
              write(1,'(3es12.5)') eninc(nen),0.,0.
  210       continue
          else
            open (unit=1,file=Efile,status='old')
            do 230 nen=1,nin+4
              read(1,*,end=240,err=240)
  230       continue
          endif
          write(1,'(2es12.5)') Einc,Eaverage(type)
  240     close (unit=1)
        endif
  110 continue
      return
      end
Copyright (C)  2019 A.J. Koning, S. Hilaire and S. Goriely

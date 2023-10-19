      subroutine totalout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 16, 2016
c | Task  : Output of total cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*18 totfile
      integer      nen
c
c ********************* Total cross sections ***************************
c
c Einc        : incident energy in MeV
c nin         : counter for incident energy
c flaginitpop : flag for initial population distribution
c xsinitpop   : initial population cross section
c eninccm     : center-of-mass incident energy in MeV
c k0          : index of incident particle
c xstotinc    : total cross section (neutrons only) for incident channel
c xselasinc   : total elastic cross section (neutrons only) for incident
c               channel
c xsreacinc   : reaction cross section for incident channel
c xscompel    : compound elastic cross section
c xsnonel     : non-elastic cross section
c xsdirdiscsum: total direct cross section
c xspreeqsum  : total preequilibrium cross section summed over particles
c flaggiant   : flag for collective contribution from giant resonances
c xsgrsum     : sum over giant resonance cross sections
c xscompnonel : total compound non-elastic cross section
c xselastot   : total elastic cross section (shape + compound)
c
      if (Einc.ge.0.001) then
        write(*,'(/" ########### REACTION SUMMARY FOR E=",f10.5,
     +    " ###########"/)') Einc
      else
        write(*,'(/" ########### REACTION SUMMARY FOR E=",es12.5,
     +    " ###########"/)') Einc
      endif
      if (flaginitpop) then
        write(*,'(" 1. Initial population cross section =",es12.5/)')
     +    xsinitpop
        return
      endif
      write(*,'(" Center-of-mass energy: ",f8.3/)') eninccm
      write(*,'(" 1. Total (binary) cross sections"/)')
      if (k0.eq.1) write(*,'(" Total           =",es12.5)') xstotinc
      if (k0.eq.1) write(*,'("   Shape elastic   =",es12.5)')
     +  xselasinc
      write(*,'("   Reaction        =",es12.5)') xsreacinc
      write(*,'("     Compound elastic=",es12.5)') xscompel
      write(*,'("     Non-elastic     =",es12.5)') xsnonel
      write(*,'("       Direct          =",es12.5)') xsdirdiscsum
      write(*,'("       Pre-equilibrium =",es12.5)') xspreeqsum
      if (flagracap) write(*,'("       Direct Capture  =",es12.5)')
     +  xsracape
      if (flaggiant) write(*,'("       Giant resonance =",es12.5)')
     +  xsgrsum
      write(*,'("       Compound non-el =",es12.5)') xscompnonel
      if (k0.eq.1) write(*,'("     Total elastic   =",es12.5)')
     +  xselastot
c
c Write results to separate file
c
c filetotal: flag for total cross sections on separate file
c natstring: string extension for file names
c iso      : counter for isotope
c numinclow: number of incident energies below Elow
c parsym   : symbol of particle
c k0       : index of incident particle
c Atarget  : mass number of target nucleus
c Starget  : symbol of target nucleus
c Ztarget  : charge number of target nucleus
c numinc   : number of incident energies
c
      if (filetotal) then
        totfile='total.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,file=totfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("#     E      Non-elastic    Elastic     Total",
     +      "     Comp. el.   Shape el.   Reaction",
     +      "   Comp. nonel   Direct    Pre-equil.  Dir. Capt.")')
          do 10 nen=1,numinclow
            write(1,'(11es12.5)') eninc(nen),fxsnonel(nen),
     +        fxselastot(nen),fxstotinc(nen),fxscompel(nen),
     +        fxselasinc(nen),fxsreacinc(nen),fxscompnonel(nen),
     +        fxsdirdiscsum(nen),fxspreeqsum(nen),fxsracape(nen)
   10     continue
        else
          open (unit=1,file=totfile,status='old')
          do 20 nen=1,nin+4
            read(1,*,end=30,err=30)
   20     continue
        endif
        write(1,'(11es12.5)') Einc,xsnonel,xselastot,
     +    xstotinc,xscompel,xselasinc,xsreacinc,xscompnonel,
     +    xsdirdiscsum,xspreeqsum,xsracape
   30   close (unit=1)
c
c Total cross sections (i.e. from OMP) only
c
        totfile='totalxs.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,file=totfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2," Total cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("#    E      Cross section")')
          do 40 nen=1,numinclow
            write(1,'(2es12.5)') eninc(nen),fxstotinc(nen)
   40     continue
        else
          open (unit=1,file=totfile,status='old')
          do 50 nen=1,nin+4
            read(1,*,end=60,err=60)
   50     continue
        endif
        write(1,'(2es12.5)') Einc,xstotinc
   60   close (unit=1)
c
c Elastic cross sections only
c
        totfile='elastic.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,file=totfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2," Elastic cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("#    E      Cross section")')
          do 70 nen=1,numinclow
            write(1,'(2es12.5)') eninc(nen),fxselastot(nen)
   70     continue
        else
          open (unit=1,file=totfile,status='old')
          do 80 nen=1,nin+4
            read(1,*,end=90,err=90)
   80     continue
        endif
        write(1,'(2es12.5)') Einc,xselastot
   90   close (unit=1)
c
c Nonelastic cross sections only
c
        totfile='nonelastic.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,file=totfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2," Nonelastic cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("#    E      Cross section")')
          do 100 nen=1,numinclow
            write(1,'(2es12.5)') eninc(nen),fxsnonel(nen)
  100     continue
        else
          open (unit=1,file=totfile,status='old')
          do 110 nen=1,nin+4
            read(1,*,end=120,err=120)
  110     continue
        endif
        write(1,'(2es12.5)') Einc,xsnonel
  120   close (unit=1)
c
c Reaction cross sections only
c
c xsreacinc     : reaction cross section for incident channel
c
        totfile='reaction.tot'//natstring(iso)
        if (nin.eq.numinclow+1) then
          open (unit=1,file=totfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2," Reaction cross sections")')
     +      parsym(k0),Atarget,Starget
          write(1,'("# ")')
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("#    E      Cross section")')
          do 210 nen=1,numinclow
            write(1,'(2es12.5)') eninc(nen),fxsreacinc(nen)
  210     continue
        else
          open (unit=1,file=totfile,status='old')
          do 220 nen=1,nin+4
            read(1,*,end=230,err=230)
  220     continue
        endif
        write(1,'(2es12.5)') Einc,xsreacinc
  230   close (unit=1)
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

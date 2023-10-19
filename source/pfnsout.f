      subroutine pfnsout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 24, 2021
c | Task  : Output of fission neutrons and spectra
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*90 pfnsfile
      integer      type,nen
c
c Write results to separate files
c
c pfnsfile   : file for PFNS
c natstring : string extension for file names
c iso       : counter for isotope
c Einc0     : incident energy in MeV
c parsym    : symbol of particle
c
c Write PFNS, PFGS, etc.
c
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (pfnsmodel.eq.1.and.type.ne.1) goto 10
        write(*,'(/" +++ Prompt fission ",a8," spectrum +++")')
     +    parname(type)
        pfnsfile='pfxs0000.000.fis'//natstring(iso)
        pfnsfile(3:3)=parsym(type)
        if (Einc0.lt.0.001) then
          write(pfnsfile(5:12),'(es8.2)') Einc0
        else
          write(pfnsfile(5:12),'(f8.3)') Einc0
          write(pfnsfile(5:8),'(i4.4)') int(Einc0)
        endif
        open (unit=1,file=pfnsfile,status='replace')
        write(1,'("# ",a1," + ",i3,a2,": Prompt fission ",a8,
     +    " spectrum ")') parsym(k0),Atarget,Starget,parname(type)
        if (Einc0.lt.0.001) then
          write(1,'("# E-incident = ",es8.2," MeV")') Einc0
        else
          write(1,'("# E-incident = ",f8.3," MeV")') Einc0
        endif
        write(1,'("# Number of energies:",i6)') NEpfns
        write(*,'(/" E-average         = ",f8.3," MeV")') Eavpfns(type)
        write(*,'(" Weighted E-average= ",f8.3," MeV")') 
     +    Epfnsaverage(type)
        write(1,'("# E-average  = ",f8.3," MeV")') Eavpfns(type)
        write(*,'(/"       E-out         spectrum    Maxwell ratio",
     +    "   spectrum_CM"/)')
        write(1,'("#      E-out          spectrum    Maxwell ratio",
     +    "   spectrum_CM")')
        do 110 nen=1,NEpfns
          write(*,'(4es15.4)') Epfns(nen),pfns(type,nen),
     +      maxpfns(type,nen),pfnscm(type,nen)
          write(1,'(4es15.4)') Epfns(nen),pfns(type,nen),
     +      maxpfns(type,nen),pfnscm(type,nen)
 110    continue
        close (unit=1)
   10 continue
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

      subroutine nubarout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 23, 2021
c | Task  : Output of average number of fission neutrons
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*90 nufile
      integer      type,nen
c
c Write results to separate files
c
c nufile    : file for nubar
c yieldfile : file with fission yields
c natstring : string extension for file names
c iso       : counter for isotope
c Einc0     : incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c
c nu per number, P(nu) and nubar as function of Z and A
c
      if (nin0.eq.numinclow+1.and.numinclow.gt.0) then
        do nen=1,numinclow
          do type=0,6
            fnubar(nen,type)=nubar(type)
          enddo
        enddo
      endif
      write(*,'(/" +++ AVERAGE NUMBER OF PROMPT FISSION NEUTRONS +++")')
      do 10 type=0,6
        if (parskip(type)) goto 10
        if (nubar(type).eq.0.) goto 10
c
c Write nubar
c
c nubarexist : flag for existence of nubar file
c numinc     : number of incident energies
c eninc      : incident energy in MeV
c
        nufile='nubarx.tot'//natstring(iso)//'       '
        nufile(6:6)=parsym(type)
        if (.not.nubarexist(type)) then
          nubarexist(type)=.true.
          open (unit=1,file=nufile,status='replace')
          write(1,'("# ",a1," + ",i3,a2,
     +      ": Average prompt ",a8," multiplicity (nubar-prompt)")')
     +    parsym(k0),Atarget,Starget,parname(type)
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("# ")')
          write(1,'("# E-in           nubar")')
          do 40 nen=1,numinclow
            write(1,'(2es12.5)') eninc(nen),fnubar(nen,type)
   40     continue
          do 50 nen=numinclow+1,nin0-1
            write(1,'(2es12.5)') eninc(nen),0.
   50     continue
        else
          open (unit=1,file=nufile,status='old')
          do 60 nen=1,nin0+4
            read(1,*,end=70,err=70)
   60     continue
        endif
        write(1,'(2es12.5)') Einc0,nubar(type)
        write(*,'(2es12.5)') Einc0,nubar(type)
   70   close (unit=1)
   10 continue
      return
      end
Copyright (C)  2021 A.J. Koning, S. Hilaire and S. Goriely

      subroutine massdisout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : October 30, 2020
c | Task  : Output of fission fragment yields
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*12 isostring(0:1)
      character*90 fpyieldfile,fpfile
      integer      iz,ia,in,nen,nex
c
c ****************** Output of fission yields **************************
c
      isostring(0)='ground state'
      isostring(1)='isomer      '
      write(*,'(/" ++++++++++ FISSION YIELDS ++++++++++"/)')
c
c Write results to separate files
c
c yieldfile : file with fission yields
c natstring : string extension for file names
c iso       : counter for isotope
c Einc0     : incident energy in MeV
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c Ztarget   : charge number of target nucleus
c yieldApre : pre-neutron emission mass yield
c yieldApost: post-neutron emission corrected mass yield
c
      fpyieldfile='yieldA0000.000.fis'//natstring(iso)
      if (Einc0.lt.0.001) then
        write(fpyieldfile(7:14),'(es8.2)') Einc0
      else
        write(fpyieldfile(7:14),'(f8.3)') Einc0
        write(fpyieldfile(7:10),'(i4.4)') int(Einc0)
      endif
      open (unit=1,file=fpyieldfile,status='replace')
      write(1,'("# ",a1," + ",i3,a2,": Fission yields")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# E-incident = ",es12.5)') Einc0
      write(1,'("# Number of nuclides: ",i3)') Atarget
      write(1,'("# ")')
      write(*,'(" Fission yields as function of A"/)')
      write(*,'("   A       FP yield       FF yield ",
     +  "         FP xs          FF xs")')
      write(1,'("#  A       FP yield       FF yield ",
     +  "         FP xs          FF xs")')
      do 10 ia=1,Atarget
        write(*,'(i4,4es15.4)') ia,yieldApost(ia),yieldApre(ia),
     +    xsApost(ia),xsApre(ia)
        write(1,'(i4,4es15.4)') ia,yieldApost(ia),yieldApre(ia),
     +    xsApost(ia),xsApre(ia)
  10  continue
      write(*,'(/" Tot",4es15.4)') yieldtotpost,yieldtotpre,
     +  xstotpost,xstotpre
      write(1,'(/"#Tot",4es15.4)') yieldtotpost,yieldtotpre,
     +  xstotpost,xstotpre
      close (unit=1)
c
c Write ff/fp production
c
c fpexist    : flag for existence of fission product
c fpfile     : file with fission product
c yieldZApre : pre-neutron emission isotopic yield
c yieldZApost: post-neutron emission corrected isotopic yield
c
      write(*,'(/" Fission yields as function of Z, A"/)')
      write(*,'("    Z    A iso     FP yield       FF yield",
     +  "         FP xs          FF xs    Isom. Ratio")')
      fpyieldfile='yieldZA0000.000.fis'//natstring(iso)
      if (Einc0.lt.0.001) then
        write(fpyieldfile(8:15),'(es8.2)') Einc0
      else
        write(fpyieldfile(8:15),'(f8.3)') Einc0
        write(fpyieldfile(8:11),'(i4.4)') int(Einc0)
      endif
      open (unit=2,file=fpyieldfile,status='replace')
      write(2,'("# ",a1," + ",i3,a2,": Z, A Fission yields")')
     +  parsym(k0),Atarget,Starget
      write(2,'("# E-incident = ",es12.5)') Einc0
      write(2,'("# ")')
      write(2,'("# ")')
      write(2,'("#   Z    A iso      FP yield       FF yield",
     +  "         FP xs          FF xs    Isom. Ratio")')
      do 210 ia=1,Atarget
        do 220 iz=1,Ztarget
          in=ia-iz
          if (in.lt.1.or.in.gt.Ninit) goto 220
          if (xsZApre(iz,in).lt.fpeps.and.
     +      xsZApost(iz,in).lt.fpeps.and..not.fpexist(iz,in))
     +      goto 220
          fpfile='fp000000.tot'//natstring(iso)
          write(fpfile(3:8),'(2i3.3)') iz,ia
          if (.not.fpexist(iz,in)) then
            open (unit=1,file=fpfile,status='replace')
            write(1,'("# ",a1," + ",i3,a2,": Fission product yield of ",
     +        i3,a2)') parsym(k0),Atarget,Starget,ia,nuc(iz)
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') numinc
            write(1,'("# ")')
            write(1,'("# E-incident    FP yield    FF yield ",
     +        "     FP xs       FF xs")')
            do 230 nen=1,nin0-1
              write(1,'(5es12.5)') eninc(nen),0.,0.,0.,0.
  230       continue
          else
            open (unit=1,file=fpfile,status='old')
            do 240 nen=1,nin0+4
              read(1,*,end=250,err=250)
  240       continue
          endif
          if (xsZApre(iz,in).ge.fpeps.and.xsZApost(iz,in).ge.fpeps) then
            write(*,'(2i5,i3,4es15.4)') iz,ia,-1,yieldZApost(iz,in),
     +        yieldZApre(iz,in),xsZApost(iz,in),xsZApre(iz,in)
            write(2,'(2i5,i3,4es15.4)') iz,ia,-1,yieldZApost(iz,in),
     +        yieldZApre(iz,in),xsZApost(iz,in),xsZApre(iz,in)
          endif
          write(1,'(5es12.5)') Einc0,yieldZApost(iz,in),
     +      yieldZApre(iz,in),xsZApost(iz,in),xsZApre(iz,in)
  250     close (unit=1)
          if (xsfpex(iz,in,1).gt.0.) then
            do 260 nex=0,1
              write(fpfile(10:12),'("L",i2.2)') nex
              if (.not.fpexist(iz,in)) then
                open (unit=1,file=fpfile,status='replace')
                write(1,'("# ",a1," + ",i3,a2,": Fission product yield",
     +            " of ",i3,a2,1x,a12)') parsym(k0),Atarget,Starget,ia,
     +            nuc(iz),isostring(nex)
                write(1,'("# ")')
                write(1,'("# # energies =",i6)') numinc
                write(1,'("# ")')
                write(1,'("# E-incident     FP yield      FP xs ",
     +            "  Ratio ")')
                do 270 nen=1,nin0-1
                  write(1,'(5es12.5)') eninc(nen),0.,0.,0.,0.
  270           continue
              else
                open (unit=1,file=fpfile,status='old')
                do 280 nen=1,nin0+4
                  read(1,*,end=290,err=290)
  280           continue
              endif
              if (xsZApre(iz,in).ge.fpeps.and.xsZApost(iz,in).ge.fpeps)
     +          then
                write(*,'(2i5,i3,2(es15.4,15x),es15.4)') iz,ia,nex,
     +            yieldfpex(iz,in,nex),xsfpex(iz,in,nex),
     +            fpratio(iz,in,nex)
                write(2,'(2i5,i3,2(es15.4,15x),es15.4)') iz,ia,nex,
     +            yieldfpex(iz,in,nex),xsfpex(iz,in,nex),
     +            fpratio(iz,in,nex)
              endif
              write(1,'(4es12.5)') Einc0,yieldfpex(iz,in,nex),
     +          xsfpex(iz,in,nex),fpratio(iz,in,nex)
  290         close (unit=1)
  260       continue
          endif
          if (.not.fpexist(iz,in)) fpexist(iz,in)=.true.
  220   continue
c
c Write cumulative ff/fp production
c
        if (xsApre(ia).lt.fpeps.and.
     +    xsApost(ia).lt.fpeps.and..not.fpaexist(ia)) goto 210
        fpfile='fp000000.tot'//natstring(iso)
        write(fpfile(6:8),'(i3.3)') ia
        if (.not.fpaexist(ia)) then
          fpaexist(ia)=.true.
          open (unit=1,file=fpfile,status='replace')
          write(1,'("# ",a1," + ",i3,a2,": Fission product yield of A=",
     +      i3)') parsym(k0),Atarget,Starget,ia
          write(1,'("# ")')
          write(1,'("# # energies =",i6)') numinc
          write(1,'("# ")')
          write(1,'("# E-incident    FP yield    FF yield",
     +      "      FP xs       FF xs")')
          do 310 nen=1,nin0-1
            write(1,'(5es12.5)') eninc(nen),0.,0.,0.,0.
  310     continue
        else
          open (unit=1,file=fpfile,status='old')
          do 320 nen=1,nin0+4
            read(1,*,end=330,err=330)
  320     continue
        endif
        write(1,'(5es12.5)') Einc0,yieldApost(ia),yieldApre(ia),
     +    xsApost(ia),xsApre(ia)
  330   close (unit=1)
  210 continue
      close (unit=2)
      write(*,'(/"Total     ",4es15.4)') yieldtotpost,yieldtotpre,
     +  xstotpost,xstotpre
      return
      end
Copyright (C)  2019 A.J. Koning, S. Hilaire and S. Goriely

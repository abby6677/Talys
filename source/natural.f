      subroutine natural
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : February 8, 2019
c | Task  : Calculation for natural elements
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numen2nat
      parameter    (numen2nat=7*numen2)
      logical      lexist,elexist,specexist,resexist,xsexist,
     +             fissionexist,fyexist
      character*10 headstring,totstring
      character*13 prodfile
      character*15 fisfile,Yfile
      character*16 resfile,xsfile
      character*20 totfile,tot0file
      character*21 discfile,specfile,fyfile
      character*29 recfile
      character*80 str(9)
      integer      i,k,k2,type,jend,iang,j,n1,n2,nen,neniso(numiso),
     +             nenen,zbeg,zend,abeg,aend,iz,ia,npart,ih,it,id,ip,in,
     +             nex
      real         en(numen2),xst(10),xstotnat(10,numen2),ang(0:numang),
     +             xsprodnat(numen2),xsyieldnat(numen2),xsnat(0:numen2),
     +             xs1nat(0:numen2),xs2nat(0:numen2),xs,y,xs1,xs2,
     +             enspec(numiso,0:numen2nat),E,entmp,Ea,Eb,Efac,
     +             xsspec(numiso,0:numen2nat,5),enspecnat(0:numen2nat),
     +             xsspecnat(0:numen2nat,5),Nis,act,frac,
     +             actnat(0:numtime),Nisnat(0:numtime),Ynat(0:numtime),
     +             fracnat(0:numtime)
c
c **************** Create runs and directories per isotope *************
c
c numen2nat    : number of energies for natural element
c isonum       : number of isotopes
c iso          : counter for isotope
c talysinput   : subroutine for user input and defaults
c talysinitial : subroutine for initialization of nuclear structure
c talysreaction: subroutine with reaction models
c
c For mono-isotopic nuclides we are done.
c
      if (isonum.eq.1) return
c
c Do a full TALYS calculation for each isotope
c
      do iso=2,isonum
        call talysinput
        call talysinitial
        call talysreaction
      enddo
c
c ******** Merge output files in results for natural elements **********
c
c 1. Total cross sections
c
c numen2    : maximum number of outgoing energies
c en        : incident energy
c xsprodnat : production cross sections for natural element
c xsyieldnat: yields for natural element
c xst       : help variable
c jend      : end of jloop
c xstotnat  : total cross sections for natural element
c filetotal : flag for total cross sections on separate file
c totfile   : file with total cross sections
c totstring : string for totfile
c headstring: string for first part of filename
c tot0file  : file with total cross sections
c numinc    : number of incident energies
c natstring : string extension for file names
c abun      : isotopic abundance
c parsym    : symbol of particle
c nuc       : symbol of nucleus
c
      do 10 k=1,numen2
        en(k)=0.
        xsprodnat(k)=0.
        xsyieldnat(k)=0.
        do 10 k2=1,10
          xst(k2)=0.
          xstotnat(k2,k)=0.
   10 continue
      if (filetotal) then
        do 20 i=1,isonum
          totfile='total.tot'//natstring(i)
          inquire (file=totfile,exist=lexist)
          if (lexist) then
            open (2,file=totfile,status='old')
            read(2,'(////)',end=50,err=50)
            do 30 k=1,numinc
              read(2,'(11e11.4)',end=50,err=50) en(k),
     +          (xst(k2),k2=1,10)
              do 40 k2=1,10
                xstotnat(k2,k)=xstotnat(k2,k)+abun(i)*xst(k2)
   40         continue
   30       continue
   50       close (unit=2)
          endif
   20   continue
        open (3,file='total.tot',status='replace')
        write(3,'("# ",a1," + nat-",a2," Total cross sections")')
     +    parsym(k0),Starget
        write(3,'("# ")')
        write(3,'("# ")')
        write(3,'("# # energies =",i6)') numinc
        write(3,'("#    E      Non-elastic  Elastic     Total",
     +    "     Comp. el.  Shape el.  Reaction",
     +    " Comp. nonel   Direct   Pre-equil.   Dir. Capt.")')
        do 60 k=1,numinc
          write(3,'(11es11.4)') en(k),(xstotnat(k2,k),k2=1,10)
   60   continue
        close (unit=3)
        do 70 k=1,numen2
          en(k)=0.
          do 70 k2=1,10
            xst(k2)=0.
            xstotnat(k2,k)=0.
   70   continue
        if (k0.eq.1) then
          jend=4
        else
          jend=3
        endif
        do 72 j=1,jend
          do 74 i=1,isonum
            if (j.eq.1) totstring='totalxs'
            if (j.eq.2) totstring='nonelastic'
            if (j.eq.3) totstring='reaction'
            if (j.eq.4) totstring='elastic'
            tot0file=trim(totstring)//'.tot'
            totfile=trim(tot0file)//natstring(i)
            inquire (file=totfile,exist=lexist)
            if (lexist) then
              open (2,file=totfile,status='old')
              read(2,'(////)',end=78,err=78)
              do 76 k=1,numinc
                read(2,'(2es12.5)',end=78,err=78) en(k),xst(j)
                  xstotnat(j,k)=xstotnat(j,k)+abun(i)*xst(j)
   76         continue
   78         close (unit=2)
            endif
   74     continue
          if (j.eq.1) headstring='Total     '
          if (j.eq.2) headstring='Nonelastic'
          if (j.eq.3) headstring='Reaction  '
          if (j.eq.4) headstring='Elastic   '
          open (3,file=tot0file,status='replace')
          write(3,'("# ",a1," + nat-",a2," ",a10," cross sections")')
     +      parsym(k0),Starget,headstring
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# # energies =",i6)') numinc
          write(3,'("#    E      Cross section")')
          do 80 k=1,numinc
            write(3,'(2es12.5)') en(k),xstotnat(j,k)
   80     continue
          close (unit=3)
   72   continue
c
c 2. Particle production cross sections
c
c prodfile: file with total particle production cross sections
c parname : name of particle
c
        do 110 type=1,6
          do 120 k=1,numinc
            xsprodnat(k)=0.
            xsyieldnat(k)=0.
  120     continue
          do 130 i=1,isonum
            prodfile=' prod.tot'//natstring(i)
            write(prodfile(1:1),'(a1)') parsym(type)
            inquire (file=prodfile,exist=lexist)
            if (lexist) then
              open (2,file=prodfile,status='old')
              read(2,'(////)',end=150,err=150)
              do 140 k=1,numinc
                read(2,'(3es12.5)',end=150,err=150) en(k),xs,y
                xsprodnat(k)=xsprodnat(k)+abun(i)*xs
                xsyieldnat(k)=xsyieldnat(k)+abun(i)*y
  140         continue
  150         close (unit=2)
            endif
  130     continue
          open (3,file=prodfile(1:9),status='replace')
          write(3,'("# ",a1," + nat-",a2," Total ",a8," production")')
     +      parsym(k0),Starget,parname(type)
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# # energies =",i6)') numinc
          write(3,'("#    E         xs         Yield")')
          do 160 k=1,numinc
            write(3,'(3es12.5)') en(k),xsprodnat(k),xsyieldnat(k)
  160     continue
          close (unit=3)
  110   continue
      endif
c
c 3. Elastic scattering angular distribution
c
c fileelastic: flag for elastic angular distribution on separate file
c xs1nat     : help variable
c xs2nat     : help variable
c nangle     : number of angles
c xsnat,...  : cross section for natural element
c elexist    : logical to determine existence of elastic scattering file
c discfile   : file with elastic scattering angular distribution
c angle,ang  : angle
c eninc      : incident energy in MeV
c
      if (fileelastic) then
        do 210 k=1,numinc
          do 220 iang=0,nangle
            xsnat(iang)=0.
            xs1nat(iang)=0.
            xs2nat(iang)=0.
  220     continue
          elexist=.false.
          do 230 i=1,isonum
            discfile='nn        ang.L00'//natstring(i)
            write(discfile(1:2),'(2a1)') parsym(k0),parsym(k0)
            write(discfile(3:10),'(f8.3)') eninc(k)
            write(discfile(3:6),'(i4.4)') int(eninc(k))
            inquire (file=discfile,exist=lexist)
            if (lexist) then
              elexist=.true.
              open (2,file=discfile,status='old')
              read(2,'(////)',end=250,err=250)
              do 240 iang=0,nangle
                read(2,'(f5.1,3e16.5)',end=250,err=250) ang(iang),xs,
     +            xs1,xs2
                xsnat(iang)=xsnat(iang)+abun(i)*xs
                xs1nat(iang)=xs1nat(iang)+abun(i)*xs1
                xs2nat(iang)=xs2nat(iang)+abun(i)*xs2
  240         continue
  250         close (unit=2)
            endif
  230     continue
          if (elexist) then
            open (3,file=discfile(1:16),status='replace')
            write(3,'("# ",a1," + nat-",a2," Elastic scattering",
     +        " angular distribution")') parsym(k0),Starget
            write(3,'("# E-incident = ",f8.3)') eninc(k)
            write(3,'("# ")')
            write(3,'("# # angles   =",i4)') nangle+1
            write(3,'("#   E         xs            Direct",
     +        "         Compound")')
            do 260 iang=0,nangle
              write(3,'(f5.1,3es16.5)') ang(iang),xsnat(iang),
     +          xs1nat(iang),xs2nat(iang)
  260       continue
            close (unit=3)
          endif
  210   continue
      endif
c
c 4. Composite particle spectra
c
c filespectrum: designator for spectrum on separate file
c specexist   : logical to determine existence of spectrum file
c enspecnat   : emission energy for natural element
c enspec      : emission energy
c specfile    : file with composite particle spectra
c xsspecnat   : differential cross section for natural element
c xsspec      : differential cross section
c neniso      : number of emission energies per isotope
c Efac,....   : help variables
c ebegin      : first energy point of energy grid
c eendout     : last energy point of energy grid
c
      do 310 k=1,numinc
        do 320 type=1,6
          if (.not.filespectrum(type)) goto 320
          specexist=.false.
          do 330 k2=1,numen2nat
            enspecnat(k2)=0.
            do 340 j=1,5
              xsspecnat(k2,j)=0.
  340       continue
  330     continue
          do 350 i=1,isonum
            do 360 k2=1,numen2
              enspec(i,k2)=0.
              do 370 j=1,5
                xsspec(i,k2,j)=0.
  370         continue
  360       continue
c
c In general, the secondary spectra for the various isotopes are
c all different. Therefore, we first read the secondary energy
c grids and cross sections into memory.
c
            neniso(i)=0
            specfile=parsym(type)//'spec0000.000.tot'//natstring(i)
            write(specfile(6:13),'(f8.3)') eninc(k)
            write(specfile(6:9),'(i4.4)') int(eninc(k))
            inquire (file=specfile,exist=lexist)
            if (lexist) then
              specexist=.true.
              open (2,file=specfile,status='old')
              read(2,'(////)',end=390,err=390)
              do 380 k2=1,numen2
                read(2,'(f8.3,5es12.5)',end=390,err=390)
     +            enspec(i,k2),(xsspec(i,k2,j),j=1,5)
                neniso(i)=neniso(i)+1
  380         continue
  390         close (unit=2)
            endif
  350     continue
c
c Make one unifying energy grid by removing double energies
c and sorting the remaining energies.
c
          if (specexist) then
            nenen=0
            do 400 i=1,isonum
              do 410 k2=1,neniso(i)
                E=enspec(i,k2)
                do 420 nen=1,nenen
                  if (E.eq.enspecnat(nen)) goto 410
  420           continue
                nenen=nenen+1
                enspecnat(nenen)=E
  410         continue
  400       continue
            do 430 n1=1,nenen
              do 440 n2=n1,nenen
                if (enspecnat(n1).le.enspecnat(n2)) goto 440
                entmp=enspecnat(n1)
                enspecnat(n1)=enspecnat(n2)
                enspecnat(n2)=entmp
  440         continue
  430       continue
c
c Interpolation and construction of natural spectra.
c
            do 450 nen=1,nenen
              do 460 i=1,isonum
                do 470 k2=1,neniso(i)-1
                  Ea=enspec(i,k2)
                  Eb=enspec(i,k2+1)
                  if (enspecnat(nen).ge.Ea.and.enspecnat(nen).lt.Eb)
     +              then
                    Efac=(enspecnat(nen)-Ea)/(Eb-Ea)
                    do 480 j=1,5
                      xs=xsspec(i,k2,j)+Efac*(xsspec(i,k2+1,j)-
     +                  xsspec(i,k2,j))
                      xsspecnat(nen,j)=xsspecnat(nen,j)+abun(i)*xs
  480               continue
                    goto 460
                  endif
  470           continue
  460         continue
  450       continue
            open (3,file=specfile(1:17),status='replace')
            write(3,'("# ",a1," + nat-",a2,": ",a8," spectrum")')
     +        parsym(k0),Starget,parname(type)
            write(3,'("# E-incident = ",f8.3)') eninc(k)
            write(3,'("# ")')
            write(3,'("# # energies =",i6)') nenen
            write(3,'("#  E-out    Total       Direct    Pre-equil.",
     +        "  Mult. preeq  Compound")')
            do 490 nen=1,nenen
              write(3,'(f8.3,5es12.5)') enspecnat(nen),
     +          (xsspecnat(nen,j),j=1,5)
  490       continue
            close (unit=3)
          endif
  320   continue
  310 continue
c
c 5. Residual production cross sections
c
c fileresidual: flag for residual production cross sections  on
c               separate file
c zbeg        : start of Z loop
c zend        : end of Z loop
c abeg        : start of A loop
c aend        : end of A loop
c Ztarget     : charge number of target nucleus
c Starget     : symbol of target nucleus
c isotope     : isotope number of residual nucleus
c numZ        : maximal number of protons away from the initial compound
c               nucleus
c numN        : maximal number of neutrons away from the initial
c               compound nucleus
c
      if (fileresidual) then
        zbeg=max(Ztarget-numZ-2,1)
        zend=Ztarget+2
        abeg=max(isotope(1)-numN-2,1)
        aend=isotope(isonum)+4
        do 510 iz=zbeg,zend
          do 510 ia=abeg,aend
c
c Total
c
c resexist: logical to determine existence of residual production file
c resfile : file with residual production cross sections
c
            do 520 k=1,numinc
              xsnat(k)=0.
  520       continue
            resexist=.false.
            do 530 i=1,isonum
              resfile='rp000000.tot'//natstring(i)
              write(resfile(3:8),'(2i3.3)') iz,ia
              inquire (file=resfile,exist=lexist)
              if (lexist) then
                resexist=.true.
                open (2,file=resfile,status='old')
                read(2,'(////)',end=550,err=550)
                do 540 k=1,numinc
                  read(2,'(2es12.5)',end=550,err=550) en(k),xs
                  xsnat(k)=xsnat(k)+abun(i)*xs
  540           continue
  550           close (unit=2)
              endif
  530       continue
          if (resexist) then
              open (3,file=resfile(1:12),status='replace')
              write(3,'("# ",a1," + nat-",a2,": Production of ",i3,a2,
     +          " - Total")') parsym(k0),Starget,ia,nuc(iz)
              write(3,'("# ")')
              write(3,'("# ")')
              write(3,'("# # energies =",i6)') numinc
              write(3,'("#    E         xs")')
              do 560 k=1,numinc
                write(3,'(2es12.5)') en(k),xsnat(k)
  560         continue
              close (unit=3)
            endif
c
c Per ground state and isomer
c
            write(resfile(10:12),'("L00")')
            inquire (file=resfile,exist=lexist)
            if (.not.lexist) goto 510
            do 570 nex=0,numlev
              do 580 k=1,numinc
                xsnat(k)=0.
  580         continue
              resexist=.false.
              do 590 i=1,isonum
                resfile='rp000000.L00'//natstring(i)
                write(resfile(3:8),'(2i3.3)') iz,ia
                write(resfile(11:12),'(i2.2)') nex
                inquire (file=resfile,exist=lexist)
                if (lexist) then
                  resexist=.true.
                  open (2,file=resfile,status='old')
                  read(2,'(////)',end=610,err=610)
                  do 600 k=1,numinc
                    read(2,'(2es12.5)',end=610,err=610) en(k),xs
                    xsnat(k)=xsnat(k)+abun(i)*xs
  600             continue
  610             close (unit=2)
                endif
  590         continue
              if (resexist) then
                open (3,file=resfile(1:12),status='replace')
                write(3,'("# ",a1," + nat-",a2,": Production of ",i3,a2,
     +            " - Level",i3)') parsym(k0),Starget,ia,nuc(iz),nex
                write(3,'("# ")')
                write(3,'("# ")')
                write(3,'("# # energies =",i6)') numinc
                write(3,'("#    E         xs")')
                do 620 k=1,numinc
                  write(3,'(2es12.5)') en(k),xsnat(k)
  620           continue
                close (unit=3)
              endif
  570       continue
  510   continue
      endif
c
c 6. Exclusive channel cross sections
c
c filechannels: flag for exclusive channel cross sections on
c               separate file
c npart       : number of particles in outgoing channel
c maxchannel  : maximal number of outgoing particles in individual
c               channel description (e.g. this is 3 for (n,2np))
c numin,....  : maximal number of ejectile in channel description
c xsexist     : flag for existence of cross section file
c xsfile      : file with channel cross sections
c
      if (filechannels) then
        do 710 npart=0,maxchannel
        do 711 ia=0,numia
        do 712 ih=0,numih
        do 713 it=0,numit
        do 714 id=0,numid
        do 715 ip=0,numip
        do 716 in=0,numin
          if (in+ip+id+it+ih+ia.ne.npart) goto 716
          do 720 k=1,numinc
            xsnat(k)=0.
  720     continue
          xsexist=.false.
          do 730 i=1,isonum
            xsfile='xs000000.tot'//natstring(i)
            write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
            inquire (file=xsfile,exist=lexist)
            if (lexist) then
              xsexist=.true.
              open (2,status='old',file=xsfile)
              read(2,'(////)',end=750,err=750)
              do 740 k=1,numinc
                read(2,'(2es12.5)',end=750,err=750) en(k),xs
                xsnat(k)=xsnat(k)+abun(i)*xs
  740         continue
  750         close (unit=2)
            endif
  730     continue
          if (xsexist) then
            open (3,status='replace',file=xsfile(1:12))
            write(3,'("# ",a1," + nat-",a2)') parsym(k0),Starget
            write(3,'("# ")')
            write(3,'("# ")')
            write(3,'("# # energies =",i6)') numinc
            write(3,'("#    E         xs")')
            do 760 k=1,numinc
              write(3,'(2es12.5)') en(k),xsnat(k)
  760       continue
            close (unit=3)
          endif
  716   continue
  715   continue
  714   continue
  713   continue
  712   continue
  711   continue
  710   continue
      endif
c
c 7. Fission cross sections
c
c filefission : flag for fission cross sections on separate file
c fissionexist: logical to determine existence of fission file
c fisfile     : file with fission cross sections
c
      if (filefission) then
        do 810 k=1,numinc
          xsnat(k)=0.
  810   continue
        fissionexist=.false.
        do 820 i=1,isonum
          fisfile='fission.tot'//natstring(i)
          inquire (file=fisfile,exist=lexist)
          if (lexist) then
            fissionexist=.true.
            open (2,status='old',file=fisfile)
            read(2,'(////)',end=840,err=840)
            do 830 k=1,numinc
              read(2,'(2es12.5)',end=840,err=840) en(k),xs
              xsnat(k)=xsnat(k)+abun(i)*xs
  830       continue
  840       close (unit=2)
          endif
  820   continue
        if (fissionexist) then
          open (3,status='unknown',file='fission.tot')
          write(3,'("# ",a1," + nat-",a2," Total fission cross ",
     +      "section")') parsym(k0),Starget
          write(3,'("# ")')
          write(3,'("# ")')
          write(3,'("# # energies =",i6)') numinc
          write(3,'("#     E         xs")')
          do 850 k=1,numinc
            write(3,'(f8.3,es12.5)') en(k),xsnat(k)
  850     continue
          close (unit=3)
        endif
      endif
c
c 8. Fission yields
c
c flagmassdis: flag for calculation of fission fragment mass yields
c Atarget    : mass number of target nucleus
c flagffevap : flag for calculation of particle evaporation from
c              fission fragment mass yields
c
c Excitation function per fission product
c
      if (flagmassdis) then
        do 910 iz=1,Ztarget
          do 910 ia=1,Atarget
            do 920 k=1,numinc
              xs1nat(k)=0.
              xs2nat(k)=0.
  920       continue
            resexist=.false.
            do 930 i=1,isonum
              resfile='fp000000.tot'//natstring(i)
              write(resfile(3:8),'(2i3.3)') iz,ia
              inquire (file=resfile,exist=lexist)
              if (lexist) then
                resexist=.true.
                open (2,status='old',file=resfile)
                read(2,'(////)',end=950,err=950)
                do 940 k=1,numinc
                  read(2,'(e10.3,es12.4,3x,es12.4)',end=950,err=950)
     +              en(k),xs1,xs2
                  xs1nat(k)=xs1nat(k)+abun(i)*xs1
                  xs2nat(k)=xs2nat(k)+abun(i)*xs2
  940           continue
  950           close (unit=2)
              endif
  930       continue
            if (resexist) then
              open (3,status='unknown',file=resfile(1:12))
              write(3,'("# ",a1," + nat-",a2,": ff yield of ",i3,a2)')
     +          parsym(k0),Starget,ia,nuc(iz)
              write(3,'("# ")')
              write(3,'("# ")')
              write(3,'("# # energies =",i6)') numinc
              write(3,'("#    E         xs")')
              if (flagffevap) then
                write(3,'("# E-incident   FF Yield   FP yield")')
                do 960 nen=1,numinc
                  write(3,'(3es12.5)') eninc(nen),
     +              xs1nat(nen),xs2nat(nen)
  960           continue
              else
                write(3,'("# E-incident   FF Yield")')
                do 970 nen=1,numinc
                  write(3,'(2es12.5)') eninc(nen),xs1nat(nen)
  970           continue
              endif
              close (unit=3)
            endif
  910   continue
c
c Mass distribution per incident energy
c
c fyexist: flag for existence of fission yields
c fyfile: file with fission yields
c
        do 1010 k=1,numinc
          do 1020 ia=1,Atarget
            xs1nat(ia)=0.
            xs2nat(ia)=0.
 1020     continue
          fyexist=.false.
          do 1030 i=1,isonum
            fyfile='yield0000.000.fis'//natstring(i)
            write(fyfile(6:13),'(f8.3)') eninc(k)
            write(fyfile(6:9),'(i4.4)') int(eninc(k))
            inquire (file=fyfile,exist=lexist)
            if (lexist) then
              fyexist=.true.
              open (2,status='old',file=fyfile)
              read(2,'(////)',end=1050,err=1050)
              do 1040 ia=1,isotope(i)
                read(2,'(3x,2e15.4)',end=1050,err=1050) xs1,xs2
                xs1nat(ia)=xs1nat(ia)+abun(i)*xs1
                xs2nat(ia)=xs2nat(ia)+abun(i)*xs2
 1040         continue
 1050         close (unit=2)
            endif
 1030     continue
          if (fyexist) then
            open (3,status='unknown',file=fyfile(1:16))
            write(3,'("# ",a1," +  nat-",a2,": mass yields")')
     +        parsym(k0),Starget
            write(3,'("# E-incident = ",f8.3)') eninc(k)
            write(3,'("# ")')
            write(3,'("# ")')
            write(3,'("# Mass    Yield   Corrected yield")')
            do 1060 ia=1,Atarget
              write(3,'(i3,3x,es12.4,3x,es12.4)') ia,xs1nat(ia),
     +          xs2nat(ia)
 1060       continue
            close (unit=3)
          endif
 1010   continue
      endif
c
c 9. Recoil spectra
c
c recfile: file with recoil spectra
c
c flagrecoil : flag for calculation of recoils
c
      if (flagrecoil) then
        zbeg=max(Ztarget-numZ-2,1)
        zend=Ztarget+2
        abeg=max(isotope(1)-numN-2,1)
        aend=isotope(isonum)+4
        do 1110 iz=zbeg,zend
          do 1110 ia=abeg,aend
            do 1120 k=1,numinc
              specexist=.false.
              do 1130 k2=1,numen2nat
                enspecnat(k2)=0.
                xsspecnat(k2,1)=0.
 1130         continue
              do 1140 i=1,isonum
                do 1150 k2=1,numen2
                  enspec(i,k2)=0.
                  xsspec(i,k2,1)=0.
 1150           continue
c
c In general, the recoil spectra for the various isotopes are
c all different. Therefore, we first read the recoil energy
c grids and cross sections into memory.
c
                neniso(i)=0
                recfile='rec000000spec0000.000.tot'//natstring(i)
                write(recfile(4:9),'(2i3.3)') iz,ia
                write(recfile(14:21),'(f8.3)') eninc(k)
                write(recfile(14:17),'(i4.4)') int(eninc(k))
                inquire (file=recfile,exist=lexist)
                if (lexist) then
                  specexist=.true.
                  open (2,status='old',file=recfile)
                  read(2,'(////)',end=1170,err=1170)
                  do 1160 k2=1,numen2
                    read(2,'(f8.3,es12.5)',end=1170,err=1170)
     +                enspec(i,k2),xsspec(i,k2,1)
                    neniso(i)=neniso(i)+1
 1160             continue
 1170             close (unit=2)
                endif
 1140         continue
c
c Make one unifying energy grid by removing double energies
c and sorting the remaining energies.
c
c nenen: energy counter
c entmp : help variable
c
              if (specexist) then
                nenen=0
                do 1200 i=1,isonum
                  do 1210 k2=1,neniso(i)
                    E=enspec(i,k2)
                    do 1220 nen=1,nenen
                      if (E.eq.enspecnat(nen)) goto 1210
 1220               continue
                    nenen=nenen+1
                    enspecnat(nenen)=E
 1210             continue
 1200           continue
                do 1230 n1=1,nenen
                  do 1240 n2=n1,nenen
                    if (enspecnat(n1).le.enspecnat(n2)) goto 1240
                    entmp=enspecnat(n1)
                    enspecnat(n1)=enspecnat(n2)
                    enspecnat(n2)=entmp
 1240             continue
 1230           continue
c
c Interpolation and construction of natural spectra.
c
                do 1250 nen=1,nenen
                  do 1260 i=1,isonum
                    do 1270 k2=1,neniso(i)-1
                      Ea=enspec(i,k2)
                      Eb=enspec(i,k2+1)
                      if (enspecnat(nen).ge.Ea.and.enspecnat(nen).lt.Eb)
     +                  then
                        Efac=(enspecnat(nen)-Ea)/(Eb-Ea)
                        xs=xsspec(i,k2,1)+Efac*(xsspec(i,k2+1,1)-
     +                    xsspec(i,k2,1))
                        xsspecnat(nen,1)=xsspecnat(nen,1)+abun(i)*xs
                        goto 1260
                      endif
 1270               continue
 1260             continue
 1250           continue
                open (3,status='unknown',file=recfile(1:24))
                write(3,'("# ",a1," + nat-",a2,": recoil spectrum for",
     +            i3,a2)') parsym(k0),Starget,ia,nuc(iz)
                write(3,'("# E-incident = ",f8.3)') eninc(k)
                write(3,'("# ")')
                write(3,'("# # energies =",i6)') nenen
                write(3,'("#  E-out    Cross section ")')
                do 1280 nen=1,nenen
                  write(3,'(f8.3,es12.5)') enspecnat(nen),
     +              xsspecnat(nen,1)
 1280           continue
                close (unit=3)
              endif
 1120     continue
 1110   continue
      endif
c
c 10. Isotope production
c
c Merge output files in results for natural elements
c
c flagprod   : flag for isotope production
c maxZ       : maximal number of protons away from the initial
c              compound nucleus
c Z          : charge number of nucleus
c Zinit      : charge number of initial compound nucleus
c maxN       : maximal number of neutrons away from the initial
c              compound nucleus
c N          : neutron number of nucleus
c Ninit      : neutron number of initial compound nucleus
c A          : mass number of nucleus
c Nis        : number of produced isotopes
c Nisnat     : number of produced isotopes for natural element
c Nisomer    : number of isomers for this nuclide
c prodexist  : logical to determine existence of residual production
c               file
c natstring  : string extension for file names
c Yfile      : file with production yields
c state      : state of final nuclide
c numtime    : maximum number of time points
c activitynat: activity of produced element in MBq
c activity   : yield of produced isotope in MBq
c yieldnat   : activity of produced element in MBq/(mA.h)
c yield      : yield of produced isotope in MBq/(mA.h)
c act        : activity
c actnat     : activity for natural element
c Nisonat    : number of isotopes produced after irradiation
c              for natural target
c Nisototnat : number of elemental isotopes produced after
c              irradiation for natural target
c Nisorelnat : fraction of number of produced isotopes per element
c              for natural target
c Nelrel     : relative amount of produced element
c fracnat    : fraction for natural element
c Td         : half life per time unit
c Tp         : irradiation time with maximal yield per time unit
c prate      : production rate per isotope
c lambda     : decay rate per isotope
c Ynat       : yield for natural element
c Tgrid      : time
c
      if (flagprod) then
        zbeg=max(Ztarget-numZ-2,1)
        zend=Ztarget+2
        abeg=max(isotope(1)-numN-2,1)
        aend=isotope(isonum)+4
        do 1310 iz=zbeg,zend
          do 1310 ia=abeg,aend
            do 1320 nex=-1,min(numlev,99)
              do k=0,numtime
                actnat(k)=0.
                Nisnat(k)=0.
                Ynat(k)=0.
                fracnat(k)=0.
              enddo
              resexist=.false.
              do 1330 i=1,isonum
                Yfile='Y000000.tot'//natstring(i)
                write(Yfile(2:7),'(2i3.3)') iz,ia
                if (nex.ge.0) write(Yfile(9:11),'("L",i2.2)') nex
                inquire (file=Yfile,exist=lexist)
                if (lexist) then
                  resexist=.true.
                  open (2,status='old',file=Yfile)
                  do k=1,9
                    read(2,'(a80)') str(k)
                  enddo
                  do 1340 k=1,numtime
                    read(2,'(f8.1,3es15.5,f15.5)',end=1350,err=1350) 
     +                xs,act,Nis,Y,frac
                    actnat(k)=actnat(k)+abun(i)*act
                    Nisnat(k)=Nisnat(k)+abun(i)*Nis
                    Ynat(k)=Ynat(k)+abun(i)*Y
                    fracnat(k)=fracnat(k)+abun(i)*frac
 1340             continue
 1350             close (unit=2)
                endif
 1330         continue
              if (resexist) then
                open (3,status='unknown',file=Yfile(1:11))
                write(3,'("# Reaction: ",a1," + nat",a2,a58)')
     +             parsym(k0),nuc(Ztarget),str(1)(22:80)
                do k=2,9
                  write(3,'(a80)') str(k)
                enddo
                do 1360 k=1,numtime
                  write(3,'(f8.1,3es15.5,f15.5)') Tgrid(k),
     +              actnat(k),Nisnat(k),Ynat(k),fracnat(k)
 1360           continue
              endif
 1320       continue
 1310     continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

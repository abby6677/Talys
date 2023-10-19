      subroutine channelsout
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 21, 2020
c | Task  : Output of exclusive reaction channels
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*3  isostring
      character*8  spstring
      character*12 isofile,gamfile
      character*16 xsfile
      character*21 spfile
      integer      npart,ia,ih,it,id,ip,in,nex,ident,idc,Zcomp,Ncomp,NL,
     +             nen,Ngam,i1,i2,type
      real         emissum,xs
c
c ****************** Output of channel cross sections ******************
c
c npart     : number of particles in outgoing channel
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))
c numia,....: maximal number of ejectile in channel description
c chanopen  : flag to open channel with first non-zero cross section
c idc,ident : help variables
c idnum     : counter for exclusive channel
c idchannel : identifier for exclusive channel
c xschannel : channel cross section
c xseps     : limit for cross sections
c reacstring: string for exclusive reaction channel
c Zcomp     : charge number index for compound nucleus
c Ncomp     : neutron number index for compound nucleus
c Nlast,NL  : last discrete level
c tau       : lifetime of state in seconds
c edis      : energy of level
c xschaniso : channel cross section per isomer
c exclbranch: exclusive channel yield per isomer
c
      write(*,'(/" 6. Exclusive cross sections"/)')
      write(*,'(" 6a. Total exclusive cross sections "/)')
      write(*,'("     Emitted particles     cross section reaction",
     +  "         level    isomeric    isomeric    lifetime",
     +  " relative yield")')
      write(*,'("    n   p   d   t   h   a",40x,"cross section",
     +  "   ratio")')
      do 10 npart=0,maxchannel
      do 11 ia=0,numia
      do 12 ih=0,numih
      do 13 it=0,numit
      do 14 id=0,numid
      do 15 ip=0,numip
      do 16 in=0,numin
        if (in+ip+id+it+ih+ia.ne.npart) goto 16
        if (.not.chanopen(in,ip,id,it,ih,ia)) goto 16
        ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
        do 20 idc=0,idnum
          if (idchannel(idc).eq.ident) then
            if (xschannel(idc).lt.xseps) goto 16
            write(*,'(1x,6i4,3x,es12.5,2x,a17,40x,es12.5)') 
     +        in,ip,id,it,ih,ia,xschannel(idc),reacstring(idc),
     +        yieldchannel(idc)
            Zcomp=ip+id+it+2*ih+2*ia
            Ncomp=in+id+2*it+ih+2*ia
            NL=Nlast(Zcomp,Ncomp,0)
            do 30 nex=1,NL
              if (tau(Zcomp,Ncomp,nex).ne.0.) goto 40
   30       continue
            goto 16
   40       write(*,'(61x,"0    ",es12.5,f9.5)')
     +        xschaniso(idc,0),exclbranch(idc,0)
            do 50 nex=1,NL
              if (tau(Zcomp,Ncomp,nex).ne.0.) then
                write(*,'(59x,i3,4x,es12.5,f9.5,2x,es12.5,
     +            " sec. ")') nex,xschaniso(idc,nex),
     +            exclbranch(idc,nex),tau(Zcomp,Ncomp,nex)
              endif
   50       continue
          endif
   20   continue
   16 continue
   15 continue
   14 continue
   13 continue
   12 continue
   11 continue
   10 continue
c
c Write results on separate files
c
c isostring    : string to designate target isomer
c filechannels : flag for exclusive channel cross sections on
c                separate file
c Liso         : isomeric number of target
c chanexist    : flag for existence of exclusive cross section
c natstring    : string extension for file names
c iso          : counter for isotope
c parsym       : symbol of particle
c k0           : index of incident particle
c Atarget      : mass number of target nucleus
c Ztarget      : charge number of target nucleus
c Starget      : symbol of target nucleus
c Qexcl        : Q-value for exclusive channel
c Ethrexc      : threshold incident energy for exclusive channel
c flagcompo    : flag for output of cross section components
c numinc       : number of incident energies
c numinclow    : number of incident energies below Elow
c eninc,Einc   : incident energy in MeV
c nin          : counter for incident energy
c xsratio      : ratio of exclusive cross section over residual
c                production cross section (for exclusive gamma ray
c                intensities)
c Fdir         : direct population fraction per nucleus
c Fpreeq       : preequilibrium population fraction per nucleus
c Fcomp        : compound population fraction per nucleus
c flagendf     : flag for information for ENDF-6 file
c fxsgamchannel: gamma channel cross section
c fxsgamdischan: discrete gamma channel cross section
c
      isostring='   '
      if (filechannels) then
        if (Liso.ge.1) isostring='('//isochar(Liso)//')'
        do 110 npart=0,maxchannel
        do 111 ia=0,numia
        do 112 ih=0,numih
        do 113 it=0,numit
        do 114 id=0,numid
        do 115 ip=0,numip
        do 116 in=0,numin
          if (in+ip+id+it+ih+ia.ne.npart) goto 116
          if (.not.chanopen(in,ip,id,it,ih,ia)) goto 116
          ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
          do 120 idc=0,idnum
            if (idchannel(idc).eq.ident) then
              if (xschannel(idc).lt.xseps.and.npart.ne.0.and.
     +          .not.chanexist(in,ip,id,it,ih,ia)) goto 116
              Zcomp=ip+id+it+2*ih+2*ia
              Ncomp=in+id+2*it+ih+2*ia
c
c A. Total
c
              xsfile='xs000000.tot'//natstring(iso)
              write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
              if (.not.chanexist(in,ip,id,it,ih,ia)) then
                chanexist(in,ip,id,it,ih,ia)=.true.
                open (unit=1,file=xsfile,status='unknown')
                write(1,'("# ",a1," + ",i3,a2,a3,": ",a17," Total")')
     +            parsym(k0),Atarget,Starget,isostring,reacstring(idc)
                write(1,'("# Q-value    =",es12.5)') Qexcl(idc,0)
                write(1,'("# E-threshold=",es12.5)') Ethrexcl(idc,0)
                write(1,'("# # energies =",i6)') numinc
                if (flagcompo) then
                  write(1,'("#     E          xs       gamma xs  ",
     +              "xs/res.prod.xs              Direct  ",
     +              "Preequilibrium Compound")')
                else
                  write(1,'("#     E          xs       gamma xs  ",
     +              "xs/res.prod.xs")')
                endif
                if (flagcompo) then
                  do 130 nen=1,numinclow
                    write(1,'(4es12.5,12x,3es12.5)') eninc(nen),
     +                fxschannel(nen,idc),fxsgamchannel(nen,idc),
     +                fxsratio(nen,idc),
     +                Fdir(Zcomp,Ncomp)*fxschannel(nen,idc),
     +                Fpreeq(Zcomp,Ncomp)*fxschannel(nen,idc),
     +                Fcomp(Zcomp,Ncomp)*fxschannel(nen,idc)
                      fxschannel(nen,idc)=0.
                      fxsgamchannel(nen,idc)=0.
                      fxsratio(nen,idc)=0.
  130             continue
                  do 140 nen=numinclow+1,nin-1
                    write(1,'(4es12.5,12x,3es12.5)')
     +                eninc(nen),0.,0.,0.,0.,0.,0.
  140             continue
                else
                  do 150 nen=1,numinclow
                    write(1,'(4es12.5)') eninc(nen),
     +                fxschannel(nen,idc),fxsgamchannel(nen,idc),
     +                fxsratio(nen,idc)
                      fxschannel(nen,idc)=0.
                      fxsgamchannel(nen,idc)=0.
                      fxsratio(nen,idc)=0.
  150             continue
                  do 160 nen=numinclow+1,nin-1
                    write(1,'(4es12.5)') eninc(nen),0.,0.,0.
  160             continue
                endif
              else
                open (unit=1,file=xsfile,status='old')
                do 170 nen=1,nin+4
                  read(1,*,end=180,err=180)
  170           continue
              endif
              if (flagcompo) then
                write(1,'(4es12.5,12x,3es12.5)') Einc,
     +            xschannel(idc),xsgamchannel(idc),xsratio(idc),
     +            Fdir(Zcomp,Ncomp)*xschannel(idc),
     +            Fpreeq(Zcomp,Ncomp)*xschannel(idc),
     +            Fcomp(Zcomp,Ncomp)*xschannel(idc)
              else
                write(1,'(4es12.5)') Einc,xschannel(idc),
     +            xsgamchannel(idc),xsratio(idc)
              endif
  180         close (unit=1)
c
c B. Ground state and isomers
c
c chanisoexist: flag for existence of exclusive isomeric cross section
c
              NL=Nlast(Zcomp,Ncomp,0)
              do 210 nex=1,NL
                if (tau(Zcomp,Ncomp,nex).ne.0.) goto 220
  210         continue
              goto 300
  220         do 230 nex=0,NL
                if (nex.eq.0.or.tau(Zcomp,Ncomp,nex).ne.0.) then
                  isofile='xs000000.L00'
                  write(isofile(3:8),'(6i1)') in,ip,id,it,ih,ia
                  write(isofile(11:12),'(i2.2)') nex
                  if (.not.chanisoexist(in,ip,id,it,ih,ia,nex)) then
                    chanisoexist(in,ip,id,it,ih,ia,nex)=.true.
                    open (unit=1,file=isofile,status='unknown')
                    write(1,'("# ",a1," + ",i3,a2,a3,": ",a17," Level",
     +               i3)') parsym(k0),Atarget,Starget,isostring,
     +               reacstring(idc),nex
                    write(1,'("# Q-value    =",es12.5," Elevel=",
     +                f11.6)') Qexcl(idc,nex),edis(Zcomp,Ncomp,nex)
                    write(1,'("# E-threshold=",es12.5)')
     +                Ethrexcl(idc,nex)
                    write(1,'("# # energies =",i6)') numinc
                    write(1,'("#     E          xs       Branching")')
                    do 240 nen=1,numinclow
                      write(1,'(3es12.5)') eninc(nen),
     +                  fxschaniso(nen,idc,nex),fexclbranch(nen,idc,nex)
                        fxschaniso(nen,idc,nex)=0.
                        fexclbranch(nen,idc,nex)=0.
  240               continue
                    do 250 nen=numinclow+1,nin-1
                      write(1,'(3es12.5)') eninc(nen),0.,0.
  250               continue
                  else
                    open (unit=1,file=isofile,status='old')
                    do 260 nen=1,nin+4
                      read(1,*,end=270,err=270)
  260               continue
                  endif
                  write(1,'(3es12.5)') Einc,
     +              xschaniso(idc,nex),exclbranch(idc,nex)
  270             close (unit=1)
                endif
  230         continue
c
c C. Discrete gamma-rays
c
c gamchanexist: flag for existence of exclusive discrete gamma-rays
c Ngam        : counter for discrete gamma rays
c numlev      : maximum number of included discrete levels
c
  300         if (flagendf) then
                gamfile='xs000000.gam'
                write(gamfile(3:8),'(6i1)') in,ip,id,it,ih,ia
                if (.not.gamchanexist(in,ip,id,it,ih,ia)) then
                  gamchanexist(in,ip,id,it,ih,ia)=.true.
                  open (unit=1,file=gamfile,status='unknown')
                  write(1,'("# ",a1," + ",i3,a2,a3,": ",a17,
     +              " Discrete gamma-rays")') parsym(k0),Atarget,
     +              Starget,isostring,reacstring(idc)
                  write(1,'("# Q-value    =",es12.5)') Qexcl(idc,0)
                  write(1,'("# E-threshold=",es12.5)') Ethrexcl(idc,0)
                  write(1,'("# # energies =",i6)') numinc
                  write(1,'("# L1 L2    xs         Estart      Eend")')
                  do 310 nen=1,numinclow
                    Ngam=0
                    do 320 i1=1,numlev
                      do 320 i2=0,i1
                        if (fxsgamdischan(nen,idc,i1,i2).gt.0.)
     +                    Ngam=Ngam+1
  320               continue
                    write(1,'(es12.5,i5)') eninc(nen),Ngam
                    do 330 i1=1,numlev
                      do 330 i2=0,i1
                        if (fxsgamdischan(nen,idc,i1,i2).gt.0.)
     +                    write(1,'(2i3,es12.5,2f11.6)') i1,i2,
     +                    fxsgamdischan(nen,idc,i1,i2),
     +                    edis(Zcomp,Ncomp,i1),edis(Zcomp,Ncomp,i2)
  330                 continue
  310             continue
                  do 340 nen=numinclow+1,nin-1
                    write(1,'(es12.5,i5)') eninc(nen),0
  340             continue
                else
                  open (unit=1,file=gamfile,status='old')
                  do 350 nen=1,5
                    read(1,*,end=400,err=400)
  350             continue
                  do 360 nen=1,nin-1
                    read(1,'(12x,i5)',end=400,err=400) Ngam
                    do 370 i1=1,Ngam
                      read(1,*,end=400,err=400)
  370               continue
  360             continue
                endif
                Ngam=0
                do 380 i1=1,numlev
                  do 380 i2=0,i1
                    if (xsgamdischan(idc,i1,i2).gt.0.) Ngam=Ngam+1
  380           continue
                write(1,'(es12.5,i5)') Einc,Ngam
                do 390 i1=1,numlev
                  do 390 i2=0,i1
                    if (xsgamdischan(idc,i1,i2).gt.0.)
     +                write(1,'(2i3,es12.5,2f11.6)') i1,i2,
     +                xsgamdischan(idc,i1,i2),edis(Zcomp,Ncomp,i1),
     +                edis(Zcomp,Ncomp,i2)
  390           continue
  400           close (unit=1)
              endif
            endif
  120     continue
  116   continue
  115   continue
  114   continue
  113   continue
  112   continue
  111   continue
  110   continue
      endif
c
c *************** Output of fission channel cross sections *************
c
c flagfission : flag for fission
c chanfisexist: flag for existence of exclusive fission cross section
c fisstring   : string for exclusive fission reaction channel
c xsfischannel: fission channel cross section
c xsabs       : absorption cross section
c channelsum  : sum over exclusive channel cross sections
c xsngnsum    : sum over total (projectile,gamma-ejectile) cross
c               sections
c xsnonel     : non-elastic cross section
c
      if (flagfission) then
        write(*,'(/" 6a2. Exclusive fission cross sections "/)')
        write(*,'("     Emitted particles     cross section reaction")')
        write(*,'("    n   p   d   t   h   a")')
        do 410 npart=0,maxchannel
        do 411 ia=0,numia
        do 412 ih=0,numih
        do 413 it=0,numit
        do 414 id=0,numid
        do 415 ip=0,numip
        do 416 in=0,numin
          if (in+ip+id+it+ih+ia.ne.npart) goto 416
          if (.not.chanopen(in,ip,id,it,ih,ia)) goto 416
          if (nin.eq.numinclow+1)
     +      chanfisexist(in,ip,id,it,ih,ia)=.false.
          ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
          do 420 idc=0,idnum
            if (idchannel(idc).eq.ident) goto 430
  420     continue
          goto 416
  430     if (xsfischannel(idc).lt.xseps) goto 416
          write(*,'(1x,6i4,3x,es12.5,2x,a17)') in,ip,id,it,ih,ia,
     +      xsfischannel(idc),fisstring(idc)
  416   continue
  415   continue
  414   continue
  413   continue
  412   continue
  411   continue
  410   continue
      endif
      write(*,'(/" Absorption cross section                 :",f12.5)')
     +  xsabs
      write(*,'(/" Sum over exclusive channel cross sections:",f12.5)')
     +  channelsum
      write(*,'(" (n,gn) + (n,gp) +...(n,ga) cross sections:",f12.5)')
     +  xsngnsum
      write(*,'(" Total                                    :",f12.5)')
     +  channelsum+xsngnsum
      if (flaginitpop) then
        write(*,'(" Initial population cross section         :",f12.5)')
     +    xsinitpop
      else
        write(*,'(" Non-elastic cross section                :",f12.5)')
     +    xsnonel
      endif
c
c Write results on separate files
c
c filefission: flag for fission cross sections on separate file
c
      if (filefission) then
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
              if (xsfischannel(idc).lt.xseps.and.npart.ne.0.and.
     +          .not.chanexist(in,ip,id,it,ih,ia)) goto 516
              Zcomp=ip+id+it+2*ih+2*ia
              Ncomp=in+id+2*it+ih+2*ia
              xsfile='xs000000.fis'
              write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
              if (.not.chanfisexist(in,ip,id,it,ih,ia)) then
                chanfisexist(in,ip,id,it,ih,ia)=.true.
                open (unit=1,file=xsfile,status='unknown')
                write(1,'("# ",a1," + ",i3,a2,a3,": ",a17," Fission")')
     +            parsym(k0),Atarget,Starget,isostring,fisstring(idc)
                write(1,'("# Q-value    =",es12.5)') Qexcl(idc,0)
                write(1,'("# E-threshold=",es12.5)') Ethrexcl(idc,0)
                write(1,'("# # energies =",i6)') numinc
                write(1,'("#    E         xs")')
                do 530 nen=1,nin-1
                  write(1,'(2es12.5)') eninc(nen),0.
  530           continue
              else
                open (unit=1,file=xsfile,status='old')
                do 540 nen=1,nin+4
                  read(1,*,end=550,err=550)
  540           continue
              endif
              write(1,'(2es12.5)') Einc,xsfischannel(idc)
  550         close (unit=1)
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
c ********** Check of total particle production cross section **********
c
c flagcheck : flag for output of numerical checks
c xsparticle: total particle production cross section
c xsparcheck: total particle production cross section
c
      if (flagcheck) then
        write(*,'(/" Check of particle production cross sections")')
        write(*,'(" (Only applies if non-elastic cross section is",
     +    " exhausted by exclusive cross sections)"/)')
        do 560 type=1,6
          if (parskip(type)) goto 560
          write(*,'(1x,a8,"=",es12.5,
     +      "    Summed exclusive cross sections=",es12.5)')
     +      parname(type),xsparticle(type),xsparcheck(type)
  560 continue
      endif
c
c *************** Output of channel cross section spectra **************
c
c flagspec   : flag for output of spectra
c parskip    : logical to skip outgoing particle
c ebegin     : first energy point of energy grid
c eendhigh   : last energy point for energy grid for any particle
c egrid      : outgoing energy grid
c xschannelsp: channel cross section spectra
c preeqratio : pre-equilibrium ratio
c emissum    : integrated emission spectrum
c xschancheck: integrated channel spectra
c xsdisctot  : total cross section summed over discrete states
c gmult      : continuum gamma multiplicity
c Eavchannel : channel average energy
c Qexcl      : Q-value for exclusive channel
c eninccm    : center-of-mass incident energy in MeV
c Especsum   : total emission energy
c
      if (flagspec) then
        write(*,'(/" 6b. Exclusive spectra ")')
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
            if (idchannel(idc).eq.ident) goto 630
  620     continue
          goto 616
  630     if (xschannel(idc).lt.xseps) goto 616
          write(*,'(/"      Emitted particles     ",
     +      "cross section reaction      gamma cross section")')
          write(*,'("    n   p   d   t   h   a")')
          write(*,'(1x,6i4,3x,es12.5,2x,a17,es12.5)') in,ip,id,it,
     +      ih,ia,xschannel(idc),reacstring(idc),xsgamchannel(idc)
          write(*,'(/"  Outgoing spectra"/)')
          write(*,'("   Energy  ",7(a8,4x)/)') (parname(type),type=0,6)
          do 640 nen=ebegin(0),eendhigh
            write(*,'(1x,f8.3,7es12.5)') egrid(nen),
     +        (xschannelsp(idc,type,nen),type=0,6)
  640     continue
          if (filechannels) then
            spstring='sp000000'
            write(spstring(3:8),'(6i1)') in,ip,id,it,ih,ia
            if (flagblock) then
              spfile=spstring//'.tot'
              if (.not.spchanexist(in,ip,id,it,ih,ia)) then
                spchanexist(in,ip,id,it,ih,ia)=.true.
                open (unit=1,file=spfile,status='unknown')
              else
                open (unit=1,file=spfile,status='unknown',
     +            position='append')
              endif
            else
              spfile = spstring//'E0000.000.tot'
              write(spfile(10:17), '(f8.3)') Einc
              write(spfile(10:13), '(i4.4)') int(Einc)
              open (unit=1,file=spfile,status='unknown')
            endif
            write(1,'("# ",a1," + ",i3,a2,a3,": ",a17," Spectra")')
     +        parsym(k0),Atarget,Starget,isostring,reacstring(idc)
            write(1,'("# E-incident = ",f10.5)') Einc
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') eendhigh-ebegin(0)+1
            write(1,'("#  E-out  ",7(2x,a8,2x))')
     +        (parname(type),type=0,6)
            if (npart.eq.0) then
              do 650 nen=ebegin(0),eendhigh
                write(1,'(f8.3,7es12.5)') egrid(nen),
     +            xschannelsp(idc,0,nen),(xsngnspec(type,nen),type=1,6)
  650         continue
            else
              do 660 nen=ebegin(0),eendhigh
                write(1,'(f8.3,7es12.5)') egrid(nen),
     +            (xschannelsp(idc,type,nen),type=0,6)
  660         continue
            endif
            close (unit=1)
          endif
          if (flagcheck) then
            emissum=xschancheck(idc)
            xs=0.
            if (npart.eq.1.and.in.eq.1) xs=xsdisctot(1)
            if (npart.eq.1.and.ip.eq.1) xs=xsdisctot(2)
            if (npart.eq.1.and.id.eq.1) xs=xsdisctot(3)
            if (npart.eq.1.and.it.eq.1) xs=xsdisctot(4)
            if (npart.eq.1.and.ih.eq.1) xs=xsdisctot(5)
            if (npart.eq.1.and.ia.eq.1) xs=xsdisctot(6)
            emissum=emissum+xs
            write(*,'(/"  E-av    ",7(f7.3,5x))')
     +        (Eavchannel(idc,type),type=0,6)
            write(*,'("  multi   ",f7.3,6(11x,i1))')
     +        gmult(idc),in,ip,id,it,ih,ia
            write(*,'("  Total   ",7(f7.3,5x))')
     +        gmult(idc)*Eavchannel(idc,0),
     +        in*Eavchannel(idc,1),ip*Eavchannel(idc,2),
     +        id*Eavchannel(idc,3),it*Eavchannel(idc,4),
     +        ih*Eavchannel(idc,5),ia*Eavchannel(idc,6)
            write(*,'(/" Available energy:",f10.5)')
     +        Qexcl(idc,0)+eninccm
            write(*,'(" Emission energy :",f10.5)') Especsum(idc)
            write(*,'(/" Check of integrated emission spectra:")')
            if (npart.eq.0) then
              write(*,'(" Cross section (x multiplicity)    =",
     +          es12.5)') gmult(idc)*xschannel(idc)
            else
              write(*,'(" Cross section (x multiplicity)    =",
     +          es12.5)') npart*xschannel(idc)
            endif
            if (npart.eq.1) then
              write(*,'(" Integrated spectra + discrete c.s.=",
     +          es12.5)') emissum
              write(*,'(" Integrated spectra                  =",
     +          es12.5)') xschancheck(idc)
              write(*,'(" Discrete cross sections             =",
     +          es12.5)') xs
            else
              write(*,'(" Integrated spectra                =",
     +          es12.5)') emissum
            endif
          endif
  616   continue
  615   continue
  614   continue
  613   continue
  612   continue
  611   continue
  610   continue
c
c *********** Output of fission channel cross section spectra **********
c
c eend          : last energy point of energy grid
c xsfischannelsp: fission channel cross section spectra
c xsfischancheck: integrated fission channel spectra
c
        if (flagfission) then
          write(*,'(/" 6b2. Exclusive fission spectra ")')
          do 710 npart=0,maxchannel
          do 711 ia=0,numia
          do 712 ih=0,numih
          do 713 it=0,numit
          do 714 id=0,numid
          do 715 ip=0,numip
          do 716 in=0,numin
            if (in+ip+id+it+ih+ia.ne.npart) goto 716
            if (.not.chanopen(in,ip,id,it,ih,ia)) goto 716
            ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
            do 720 idc=0,idnum
              if (idchannel(idc).eq.ident) goto 730
  720       continue
            goto 716
  730       if (xsfischannel(idc).le.xseps) goto 716
            write(*,'(/"      Emitted particles     ",
     +        "cross section reaction")')
            write(*,'("    n   p   d   t   h   a")')
            write(*,'(1x,6i4,3x,es12.5,2x,a17)') in,ip,id,it,ih,ia,
     +        xsfischannel(idc),fisstring(idc)
            write(*,'(/"  Outgoing spectra"/)')
            write(*,'("   Energy  ",7(a8,4x)/)')
     +        (parname(type),type=0,6)
            do 760 nen=ebegin(0),eend(0)
              write(*,'(1x,f8.3,7es12.5)') egrid(nen),
     +          (xsfischannelsp(idc,type,nen),type=0,6)
  760       continue
            if (filechannels) then
              spstring='sp000000'
              write(spstring(3:8),'(6i1)') in,ip,id,it,ih,ia
              if (flagblock) then
                spfile=spstring//'.fis'
                if (.not.spfischanexist(in,ip,id,it,ih,ia)) then
                  spfischanexist(in,ip,id,it,ih,ia)=.true.
                  open (unit=1,file=spfile,status='unknown')
                else
                  open (unit=1,file=spfile,status='unknown',
     +              position='append')
                endif
              else
                spfile = spstring//'E0000.000.fis'
                write(spfile(10:17), '(f8.3)') Einc
                write(spfile(10:13), '(i4.4)') int(Einc)
                open (unit=1,file=spfile,status='unknown')
              endif
              write(1,'("# ",a1," + ",i3,a2,a3,": ",a17," Spectra")')
     +          parsym(k0),Atarget,Starget,isostring,fisstring(idc)
              write(1,'("# E-incident = ",f10.5)') Einc
              write(1,'("# ")')
              write(1,'("# # energies =",i6)') eendhigh-ebegin(0)+1
              write(1,'("#  E-out  ",7(2x,a8,2x))')
     +          (parname(type),type=0,6)
              do 770 nen=ebegin(0),eendhigh
                write(1,'(f8.3,7es12.5)') egrid(nen),
     +            (xsfischannelsp(idc,type,nen),type=0,6)
  770         continue
              close (unit=1)
            endif
            if (flagcheck) then
              write(*,'(/" Check of integrated emission spectra:")')
              write(*,'(" Cross section (x multiplicity)=",es12.5)')
     +          npart*xsfischannel(idc)
              write(*,'(" Integrated spectra            =",es12.5)')
     +          xsfischancheck(idc)
            endif
  716     continue
  715     continue
  714     continue
  713     continue
  712     continue
  711     continue
  710     continue
        endif
c
c ************* Check of total particle production spectra *************
c
c xspreeq    : preequilibrium cross section per particle type and
c              outgoing energy
c xsmpreeq   : multiple pre-equilibrium emission spectrum
c xsgr       : smoothed giant resonance cross section
c xscomp     : compound emission spectrum
c xsspeccheck: total particle production spectra
c
        if (flagcheck) then
          write(*,'(/" Check of total particle production spectra"/)')
          write(*,'(" (Only applies if non-elastic cross section",
     +      " is exhausted by exclusive cross sections)")')
          do 810 type=0,6
            if (parskip(type)) goto 810
            write(*,'(/" Summed exclusive ",a8," spectra"/)')
     +        parname(type)
            write(*,'("  Energy   Summed     Composite   Difference"/)')
            do 820 nen=ebegin(type),eend(type)
              xs=xspreeq(type,nen)+xsmpreeq(type,nen)+xsgr(type,nen)+
     +          xscomp(type,nen)
              write(*,'(1x,f8.3,2es12.5,2x,es12.5)') egrid(nen),
     +          xsspeccheck(type,nen),xs,xsspeccheck(type,nen)-xs
  820       continue
  810     continue
        endif
      endif
c
c ***************** Output of channel recoil spectra *******************
c
c flagrecoil: flag for calculation of recoils
c maxenrec  : number of recoil energies
c Erec      : recoil energy
c specrecoil: recoil spectrum
c filerecoil: flag for recoil spectra on separate file
c
      if (flagrecoil) then
        write(*,'(/" 6c. Exclusive recoil spectra ")')
        do 910 npart=0,maxchannel
        do 911 ia=0,numia
        do 912 ih=0,numih
        do 913 it=0,numit
        do 914 id=0,numid
        do 915 ip=0,numip
        do 916 in=0,numin
          if (in+ip+id+it+ih+ia.ne.npart) goto 916
          if (.not.chanopen(in,ip,id,it,ih,ia)) goto 916
          ident=100000*in+10000*ip+1000*id+100*it+10*ih+ia
          do 920 idc=0,idnum
            if (idchannel(idc).eq.ident) goto 930
  920     continue
          goto 916
  930     if (xschannel(idc).lt.xseps) goto 916
          Zcomp=ip+id+it+2*ih+2*ia
          Ncomp=in+id+2*it+ih+2*ia
          write(*,'(/"      Emitted particles     ",
     +      "cross section reaction")')
          write(*,'("    n   p   d   t   h   a")')
          write(*,'(1x,6i4,3x,es12.5,2x,a17)') in,ip,id,it,ih,ia,
     +      xschannel(idc),reacstring(idc)
          write(*,'(/" Recoil spectrum"/)')
          write(*,'("   Energy  Cross section"/)')
          do 940 nen=0,maxenrec
            write(*,'(1x,f8.3,es12.5)') Erec(Zcomp,Ncomp,nen),
     +        specrecoil(Zcomp,Ncomp,nen)*xsratio(idc)
  940     continue
          if (filechannels.and.filerecoil) then
            spstring='sp000000'
            write(spstring(3:8),'(6i1)') in,ip,id,it,ih,ia
            if (flagblock) then
              spfile=spstring//'.rec'
              if (.not.recchanexist(in,ip,id,it,ih,ia)) then
                recchanexist(in,ip,id,it,ih,ia)=.true.
                open (unit=1,file=spfile,status='unknown')
              else
                open (unit=1,file=spfile,status='unknown',
     +            position='append')
              endif
            else
              spfile = spstring//'E0000.000.rec'
              write(spfile(10:17), '(f8.3)') Einc
              write(spfile(10:13), '(i4.4)') int(Einc)
              open (unit=1,file=spfile,status='unknown')
            endif
            write(1,'("# ",a1," + ",i3,a2,a3,": ",a17,
     +        " Recoil Spectrum")') parsym(k0),Atarget,Starget,
     +        isostring,reacstring(idc)
            write(1,'("# E-incident = ",f10.5)') Einc
            write(1,'("# ")')
            write(1,'("# # energies =",i6)') maxenrec+1
            write(1,'("#  E-out Cross section")')
            do 950 nen=0,maxenrec
              write(1,'(f8.3,es12.5)') Erec(Zcomp,Ncomp,nen),
     +          specrecoil(Zcomp,Ncomp,nen)*xsratio(idc)
  950       continue
            close (unit=1)
          endif
  916   continue
  915   continue
  914   continue
  913   continue
  912   continue
  911   continue
  910   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

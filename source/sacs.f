        subroutine sacs
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 2, 2015
c | Task  : Statistical analysis of cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer numP
      parameter (numP=1000000)
      character*16 xsfile
      integer       npart,ia,ih,it,id,ip,in,ident,idc,nen,istat,Nxs
      real          Exs(0:numP),xs(0:numP),Emax,xsmax,Ea,Eb,xsa,xsb,
     +              Eh1,Eh2,xshalf,Ewidth
c
c ********* Read exclusive cross sections ******************************
c
c parsym    : symbol of particle
c k0        : index of incident particle
c Atarget   : mass number of target nucleus
c Starget   : symbol of target nucleus
c npart     : number of particles in outgoing channel
c maxchannel: maximal number of outgoing particles in individual
c             channel description (e.g. this is 3 for (n,2np))
c numia,....: maximal number of ejectile in channel description
c chanopen  : flag to open channel with first non-zero cross section
c idnum     : counter for exclusive channel
c idchannel : identifier for exclusive channel
c Ztarget   : charge number of target nucleus
c Ethrexc   : threshold incident energy for exclusive channel
c Emax      : energy at maximum
c xsmax     : maximum cross sections
c Ea,Eb,....: help variables
c Eh1       : help variable
c Eh2       : help variable
c xshalf    : cross section at half maximum
c Ewidth    : full width at half maximum
c
      open (unit=1,file='sacs.dat',status='replace')
      write(1,'("# ",a1," + ",i3,a2,
     +  ": Statistical analysis of cross sections")')
     +  parsym(k0),Atarget,Starget
      write(1,'("# Z   A     channel    E-thresh.   Emax    ",
     +  " xsmax    width ")')
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
            xsfile='xs000000.tot'
            write(xsfile(3:8),'(6i1)') in,ip,id,it,ih,ia
            open (unit=3,file=xsfile,status='old',iostat=istat)
            if (istat.eq.0) then
              Exs(0)=0.
              xs(0)=0.
              Emax=0.
              xsmax=0.
              read(3,'(///,14x,i6,/)') Nxs
              do 30 nen=1,Nxs
                read(3,*) Exs(nen),xs(nen)
   30         continue
              close (unit=3)
              do 40 nen=1,Nxs
                if (xs(nen).gt.xsmax) then
                  Emax=Exs(nen)
                  xsmax=xs(nen)
                endif
   40         continue
              xshalf=0.5*xsmax
              Eh1=0.
              Eh2=0.
              Ewidth=0.
              do 50 nen=1,Nxs
                if (xs(nen-1).lt.xshalf.and.xs(nen).ge.xshalf) then
                  Ea=Exs(nen-1)
                  Eb=Exs(nen)
                  xsa=xs(nen-1)
                  xsb=xs(nen)
                  call pol1(xsa,xsb,Ea,Eb,xshalf,Eh1)
                endif
                if (xs(nen-1).gt.xshalf.and.xs(nen).le.xshalf) then
                  Ea=Exs(nen-1)
                  Eb=Exs(nen)
                  xsa=xs(nen-1)
                  xsb=xs(nen)
                  call pol1(xsa,xsb,Ea,Eb,xshalf,Eh2)
                  goto 60
                endif
   50         continue
   60         if (Eh1.gt.0..and.Eh2.gt.0.) Ewidth=Eh2-Eh1
              write(1,'(2i4,1x,a15,2f8.3,es12.5,f10.3)')
     +          Ztarget,Atarget,xsfile,Ethrexcl(idc,0),Emax,xsmax,Ewidth
            endif
          endif
   20   continue
   16 continue
   15 continue
   14 continue
   13 continue
   12 continue
   11 continue
   10 continue
      close (unit=1)
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

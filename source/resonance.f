        subroutine resonance
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 12, 2016
c | Task  : Reconstruction and broadening of resonance information
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer numP
      parameter (numP=1000000)
      logical      lexist
      character*9  afile
      character*12 xsfile,rpfile,pfile
      character*72 resfile
      character*72 outfile
      character*132 infile
      character*80 string,headstring(5)
      integer      iT,Ntemp,MF,MT,i,nlin,j,k,NR,NPr,NP0,NT,Ntot,Nh,
     +             ifile,Zix,Nix,nen
      real         temp,x(3),y(3),E(numP),xs(numP),Et(numP),xst(numP),
     +             fact0,am,convfac,fex,ee,ee0,ee1,smacs
      double precision fact,sum,xxs,term,term0,dE,rate
c
c prepro   : PREPRO routines for making cross section pointwise
c Liso     : isomeric number of target
c path     : directory containing structure files to be read
c reslib   : library with resonance parameters
c Tres     : temperature for broadening low energy cross sections
c flaggroup: flag for output of low energy groupwise cross sections
c infile   : input file
c iT       : counter
c Nh       : help variable
c NP0      : number of points
c NPr      : number of points
c Ntemp    : number of temperatures
c
      if (Liso.eq.0) then
        afile='000.res. '
      else
        afile='000i.res.'
        write(afile(4:4),'(a1)') isochar(Liso)
      endif
      write(afile(1:3),'(i3.3)') Atarget
      resfile='n-'//trim(nuc(Ztarget))//trim(afile)//trim(reslib)
      infile=trim(path)//'resfiles/'//trim(reslib)//'/'
     +  //trim(resfile)
      inquire (file=infile,exist=lexist)
      if (.not.lexist) then
        write(*,'(" TALYS-warning: ",a," does not exist")') trim(infile)
        return
      endif
      open(unit=41,file=infile,status='unknown')
      open(unit=42,file=resfile,status='unknown')
      do
        read(41,'(a)',end=5) string
        write(42,'(a)') string
      enddo
    5 close (41)
      close (42)
      if (flagastro) then
        fact0=3.951776e+17
        Zix=parZ(k0)
        Nix=parN(k0)
        am=redumass(Zix,Nix,k0)
        convfac=2.45484e+5/sqrt(am)
        open(unit=2,file='astrorateres.g',status='replace')
        write(2,'("# Reaction rate for ",i3,a2,"(",a1,",g)")')
     +    Atarget,Starget,parsym(k0)
        write(2,'("#    T     kT[keV]     Rate       MACS")')
        Ntemp=nTmax
      else
        Ntemp=1
      endif
      do 10 iT=1,Ntemp
        if (flagastro) then
          temp=T9(iT)*1.e9
        else
          temp=Tres
        endif
        if (flaggroup) then
          outfile=resfile(1:11)//'group'
        else
          if (temp.gt.0.) then
            outfile=resfile(1:11)//'point'
          else
            outfile=resfile(1:11)//'point0'
          endif
        endif
        call prepro(resfile,outfile,temp,flaggroup)
        open (unit=1,file=outfile,status='old')
   20   read(1,'(a80)',err=100,end=100) string
        read(string(71:72),'(i2)') MF
        if (MF.ne.3) goto 20
        read(string(73:75),'(i3)',end=100,err=100) MT
        read(1,'(a80)',err=100,end=100) string
        read(string(45:55),'(i11)',err=100) NR
        read(string(56:66),'(i11)',err=100) NPr
        nlin=1+(NR-1)/3
        do 110 i=1,nlin
          read(1,'(a80)',err=100) string
  110   continue
        nlin=1+(NPr-1)/3
        k=0
        do 120 i=1,nlin
          read(1,'(6e11.6)',err=100) (x(j),y(j),j=1,3)
          do 130 j=1,3
            k=k+1
            E(k)=x(j)*1.e-6
            xs(k)=y(j)*1.e3
  130     continue
  120   continue
        NPr=NPr-1
        NP0=NPr
        do 140 k=NP0,1,-1
          if (xs(k).gt.0.) then
            NPr=k
            goto 150
          endif
  140   continue
  150   read(1,'(a80)',err=100,end=100) string
        if (MT.eq.1) xsfile='totalxs.tot'
        if (MT.eq.2) xsfile='elastic.tot'
        if (MT.eq.18) xsfile='fission.tot'
        rpfile='rp000000.tot'
        if (MT.eq.102) then
            xsfile='xs000000.tot'
          write(rpfile(3:5),'(i3.3)') Ztarget
          write(rpfile(6:8),'(i3.3)') Atarget+1
        endif
        if (MT.eq.103)  then
          xsfile='xs010000.tot'
          write(rpfile(3:5),'(i3.3)') Ztarget-1
          write(rpfile(6:8),'(i3.3)') Atarget
        endif
        if (MT.eq.107) then
          xsfile='xs000001.tot'
          write(rpfile(3:5),'(i3.3)') Ztarget-2
          write(rpfile(6:8),'(i3.3)') Atarget-3
        endif
c
c Read TALYS output files
c
c Et: energy
c
        do 210 ifile=1,2
          if (ifile.eq.1) then
            pfile=xsfile
          else
            if (MT.lt.102) goto 210
            pfile=rpfile
          endif
          inquire (file=pfile,exist=lexist)
          if (.not.lexist) goto 10
          open (unit=2,file=pfile,status='old')
          do i=1,5
            read(2,'(a80)') headstring(i)
          enddo
          read(headstring(4)(15:80),*) NT
          do k=1,NT
            read(2,*,err=155,end=155) Et(k),xst(k)
          enddo
  155     close (2)
          Nh=NT+1
          do k=1,NT
            if (Et(k).gt.E(NPr)) then
              Nh=k
              goto 160
            endif
          enddo
  160     Ntot=NPr+NT-Nh+1
          write(headstring(4)(15:20),'(i6)') Ntot
c
c Write final output files
c
c xxs: cross section
c
          open (unit=2,file=pfile,status='replace')
            do i=1,5
            write(2,'(a80)') headstring(i)
          enddo
          do k=1,NPr
            write(2,'(2es12.5)') E(k),xs(k)
          enddo
          do k=Nh,NT
            write(2,'(2es12.5)') Et(k),xst(k)
          enddo
          close(2)
  210   continue
        goto 20
        if (flagastro) then
          fact=3.7335e+10/(sqrt(am*(temp**3)))
          sum=0.
          term0=0.
          do 220 nen=1,Ntot
            ee=real(Et(nen)*specmass(Zix,Nix,k0))
            ee1=ee
            if (nen.gt.1) then
              ee0=real(Et(nen-1)*specmass(Zix,Nix,k0))
            else
              ee0=ee
            endif
            xxs=xst(nen)/1000.
            fex=11.605*ee/temp
            if (fex.gt.80.) goto 220
            term=fact*ee*xxs*exp(-fex)
            dE=(ee1-ee0)
            if (nen.gt.1) sum=sum+(term+term0)/2.*dE
            term0=term
  220     continue
          rate=sum
          smacs=rate/convfac/sqrt(temp)
          write(2,'(f8.4,f9.3,2es12.5)') temp,temp*0.086173e3,
     +      rate,smacs
        endif
  100 close(1)
   10 continue
      if (flagastro) close(2)
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

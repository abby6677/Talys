       subroutine thermalxs
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 8, 2019
c | Task  : Cross sections at thermal energies
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*8  thchar
      character*90 thfile
      integer      Z,A,ia,isoT,isoR,i
      real         xs,xsp,xsalpha,ald,Spair,sumcap,sump,sumalpha
c
c ********** Resonance spacings and total radiative widths *************
c
c ZZ,Z   : charge number of residual nucleus
c AA,A   : mass number of residual nucleus
c thchar : help variable
c thfile : thermal cross section file
c path   : directory containing structure files to be read
c
c Read experimental values from thermal cross section file
c Values from the table can always be overruled by values given in the
c input file.
c
c 1. Inquire whether file is present
c
      Z=ZZ(0,0,1)
      A=AA(0,0,1)
      thchar=trim(nuc(Z))//'.ther'
      thfile=trim(path)//'thermal/'//thchar
      inquire (file=thfile,exist=lexist)
      if (.not.lexist) goto 110
      open (unit=2,file=thfile,status='old')
c
c 2. Search for the isotope under consideration and read information
c
c ia           : mass number from resonance table
c xscaptherm   : thermal capture cross section
c xsptherm     : thermal (n,p) cross section
c xsalphatherm : thermal (n,a) cross section
c xsalpha      : (n,a) cross section
c xsp          : (n,p) cross section
c xs,....      : help variables
c
      xs=0.
      xsp=0.
      xsalpha=0.
   10 read(2,'(4x,3i4,3(e9.2,9x))',end=20) ia,isoT,isoR,xs,xsp,xsalpha
      if (A.ne.ia) goto 10
      if (xs.ne.0.) xscaptherm(isoR)=xs
      if (xsp.ne.0.) xsptherm(isoR)=xsp
      if (xsalpha.ne.0.) xsalphatherm(isoR)=xsalpha
      goto 10
   20 close (unit=2)
      sumcap=0.
      sump=0.
      sumalpha=0.
      do 30 i=0,numisom
        sumcap=sumcap+xscaptherm(i)
        sump=sump+xsptherm(i)
        sumalpha=sumalpha+xsalphatherm(i)
   30 continue
      if (xscaptherm(-1).eq.0.) xscaptherm(-1)=sumcap
      if (xsptherm(-1).eq.0.) xsptherm(-1)=sump
      if (xsalphatherm(-1).eq.0.) xsalphatherm(-1)=sumalpha
c
c 2. Systematics
c
c Kopecky's value for (n,gamma) cross section at thermal energy.
c (J. Kopecky, M.G. Delfini, H.A.J. van der Kamp and D. Nierop:
c Revisions and extensions of neutron capture cross-sections in
c the European Activation File EAF-3, ECN-C--92-051, July 1992.)
c
c alev,ald: level density parameter
c Spair   : help variable
c S       : separation energy per particle
c pair    : pairing energy
c
  110 if (xscaptherm(-1).eq.0.) then
        ald=alev(0,0)
        Spair=S(0,0,1)-pair(0,0)
        Spair=max(Spair,1.)
        xscaptherm(-1)=1.5e-3*(ald*Spair)**3.5
      endif
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

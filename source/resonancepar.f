      subroutine resonancepar(Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 16, 2018
c | Task  : S-wave resonance parameters
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*6  reschar
      character*90 resfile
      integer      Zix,Nix,Z,A,ia,Nrrf
      real         D0f,dD0f,gamgamf,dgamgamf
c
c ********** Resonance spacings and total radiative widths *************
c
c Eavres : average resonance energy
c Zix    : charge number index for residual nucleus
c Nix    : neutron number index for residual nucleus
c ZZ,Z   : charge number of residual nucleus
c AA,A   : mass number of residual nucleus
c reschar: help variable
c resfile: resonance parameter file
c path   : directory containing structure files to be read
c
c Read experimental values from resonance parameter file
c Resonance parameters from the table can always be overruled
c by a value given in the input file.
c
c 1. Inquire whether file is present
c
      Eavres=0.01
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      reschar=trim(nuc(Z))//'.res'
      resfile=trim(path)//'resonances/'//reschar
      inquire (file=resfile,exist=lexist)
      if (.not.lexist) goto 30
      open (unit=2,file=resfile,status='old')
c
c 2. Search for the isotope under consideration and read information
c
c ia     : mass number from resonance table
c D0     : experimental s-wave resonance spacing in eV
c dD0    : uncertainty in D0
c gamgam : experimental total radiative width in eV
c dgamgam: uncertainty in gamgam
c D0f    : help variable
c dD0f   : help variable
c gamgamf: experimental total radiative width in eV
c dgamgamf: uncertainty in gamgam
c Nrr,Nrrf: number of resonances
c
   10 read(2,'(4x,i4,2e9.2,10x,2f9.5,i4)',end=20) ia,D0f,dD0f,gamgamf,
     +  dgamgamf,Nrrf
      if (A.ne.ia) goto 10
      if (dD0f.ne.0..and.D0(Zix,Nix).eq.0.) dD0(Zix,Nix)=dD0f*1000.
      if (D0f.ne.0..and.D0(Zix,Nix).eq.0.) D0(Zix,Nix)=D0f*1000.
      if (dgamgamf.ne.0..and.gamgam(Zix,Nix).eq.0.) dgamgam(Zix,Nix)=
     +  dgamgamf
      if (gamgamf.ne.0..and.gamgam(Zix,Nix).eq.0.)
     +  gamgam(Zix,Nix)=gamgamf
      if (Nrrf.ne.0) Nrr(Zix,Nix)=Nrrf
      if (Zix.eq.0.and.Nix.eq.0) then
        if (Nrrf.gt.0.and.D0(Zix,Nix).gt.0.) 
     +    Eavres=0.5*((Nrrf-1)*D0(Zix,Nix))*1.e-6
      endif
   20 close (unit=2)
c
c 3. Tabulated value for total radiative width or systematics
c    (Kopecky, 2002)
c
c gamkopecky  : radiative width in eV by spline fit of Kopecky
c gamgamadjust: adjustable factor for radiative parameters
c             (default 1.)
c
   30 if (gamgam(Zix,Nix).eq.0.) then
        if (A.ge.40.and.A.le.250) then
          gamgam(Zix,Nix)=gamkopecky(A)
        else
          gamgam(Zix,Nix)=min(1593./(A*A),10.)
        endif
      endif
      gamgam(Zix,Nix)=gamgamadjust(Zix,Nix)*gamgam(Zix,Nix)
      if (Zix.eq.0.and.Nix.eq.0.and.gnorm.ne.-1.)
     +  gnorm=gamgamadjust(Zix,Nix)*gnorm
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

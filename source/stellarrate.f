      subroutine stellarrate
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely and Stephane Hilaire
c | Date  : November 10, 2016
c | Task  : Calculate reaction rate for a Maxwell-Boltzmann
c |           distribution
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zix,Nix,Ncomp,Zcomp,nloop,iloop,nen,i,nex,nexloop
      real             fact0,am,convfac,e,fex,ee0,ee1,smacs
      double precision fact,sum,xs,term,dE,rate
c
c ************************ Partition Function **************************
c
c fact0      : reaction rate factor for photons
c Zix        : charge number index for target nucleus
c parZ       : charge number of particle
c k0         : index of incident particle
c Nix        : neutron number index for target nucleus
c parN       : neutron number of particle
c redumass,am: reduced mass
c partfunc   : subroutine to calculate partition function
c
c Initialization
c
      fact0=3.951776e+17
      Zix=parZ(k0)
      Nix=parN(k0)
      am=redumass(Zix,Nix,k0)
c
c Calculation of the T-dependent partition function
c
      call partfunc
c
c Calculation of thermonuclear reaction rate factor
c
c convfac     : conversion factor for the average velocity at
c               temperature T9
c redumass    : reduced mass
c parZ        : charge number of particle
c parN        : neutron number of particle
c k0          : index of incident particle
c Ncomp       : neutron number index for compound nucleus
c maxNastro   : maximal number of neutron away from initial compound
c               nucleus for astrophysical calculations
c Zcomp       : proton number index for compound nucleus
c maxZastro   : maximal number of protons away from initial compound
c               nucleus for astrophysical calculations
c flagfission : flag for fission
c nTmax       : number of temperatures
c T9          : Temperature grid in 10**9 K
c fact        : reaction rate factor for particles
c sum,fex     : help variables
c numinc      : number of incident energies
c nexloop     : counter
c eninc,e     : incident energy in MeV
c xs          : cross section expressed in barn
c xsastro     : cross section for astrophysical calculation
c xsastrofis  : astrophysical fission cross section
c term        : help variable
c dE          : incident energy bin
c flagastrogs : flag for calculation of astrophysics reaction rate with
c               target in ground state only
c partf       : integrated partition function
c rateastro   : thermonuclear reaction rate factor
c rateastroex : thermonuclear reaction rate to a given excited state in
c               the residual nucleus
c macsastro   : thermonuclear reaction cross section
c macsastroex : thermonuclear reaction cross section to a given excited
c               state in the residual nucleus
c rateastrofis: thermonuclear reaction rate factor for fission
c macsastrofis: thermonuclear reaction cross section for fission
c smacs       : total conversion factor for the average velocity at
c               temperature T9
c
      convfac=2.45484e+5/sqrt(redumass(parZ(k0),parN(k0),k0))
      do 10 Ncomp=0,maxNastro
        do 20 Zcomp=0,maxZastro
          nloop=1
          if (flagfission.and.Ncomp.eq.0.and.Zcomp.eq.0) nloop=2
          if (flagracap.and.Ncomp.eq.0.and.Zcomp.eq.0) nloop=3
          do 30 iloop=0,nloop
            nexloop=0
            if (iloop.eq.0) nexloop=Nlast(Zcomp,Ncomp,0)
            do 40 nex=0,nexloop
              if (nex.gt.0.and.tau(Zcomp,Ncomp,nex).eq.0.) goto 40
              do 50 i=1,nTmax
                fact=3.7335e+10/(sqrt(am*(T9(i)**3)))
                sum=0.
                do 60 nen=1,numinc
                  e=real(eninc(nen)*specmass(Zix,Nix,k0))
                  if (nen.lt.numinc) then
                    ee1=real(eninc(nen+1)*specmass(Zix,Nix,k0))
                  else
                    ee1=e
                  endif
                  if (nen.gt.1) then
                    ee0=real(eninc(nen-1)*specmass(Zix,Nix,k0))
                  else
                    ee0=e
                  endif
                  if (iloop.eq.0) xs=xsastroex(Zcomp,Ncomp,nen,nex)/
     +              1000.
                  if (iloop.eq.1) xs=xsastro(Zcomp,Ncomp,nen)/1000.
                  if (iloop.eq.2) xs=xsastrofis(nen)/1000.
                  if (iloop.eq.3) xs=xsracap(nen)/1000.
                  fex=11.605*e/T9(i)
                  if (fex.gt.80.) goto 60
                  if (k0.ne.0) then
                    term=fact*e*xs*exp(-fex)
                  else
                    term=fact0*e**2*xs/(exp(fex)-1.)
                  endif
                  dE=(ee1-ee0)/2.
                  sum=sum+term*dE
  60            continue
                if (flagastrogs.or.iloop.eq.3) then
                  rate=sum
                else
                  rate=sum/partf(i)
                endif
                smacs=convfac*sqrt(T9(i))
                if (iloop.eq.0) then
                  rateastroex(Zcomp,Ncomp,i,nex)=rate
                  macsastroex(Zcomp,Ncomp,i,nex)=rate/smacs
                endif
                if (iloop.eq.1) then
                  rateastro(Zcomp,Ncomp,i)=rate
                  macsastro(Zcomp,Ncomp,i)=rate/smacs
                endif
                if (iloop.eq.2) then
                  rateastrofis(i)=rate
                  macsastrofis(i)=rate/smacs
                endif
                if (iloop.gt.2) then
                  rateastroracap(i)=rate
                  macsastroracap(i)=rate/smacs
                  rateastro(Zcomp,Ncomp,i)=rateastro(Zcomp,Ncomp,i)+rate
                  macsastro(Zcomp,Ncomp,i)=macsastro(Zcomp,Ncomp,i)+
     +              rate/smacs
                endif
  50          continue
  40        continue
  30      continue
  20    continue
  10  continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

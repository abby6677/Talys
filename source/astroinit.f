      subroutine astroinit
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Stephane Goriely
c | Date  : November 10, 2016
c | Task  : Initialization of astrophysics quantities
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer in,iz,nen,i,nex
c
c **** Initialization of arrays of astrophysics interest ***************
c
c maxNastro   : maximal number of neutrons away from the initial
c               compound nucleus for astrophysical calculations
c maxZastro   : maximal number of protons away from the initial
c               compound nucleus for astrophysical calculations
c numT        : number of temperatures
c numinc      : number of incident energies
c xsastro     : cross section for astrophysical calculation
c xsastroex   : cross section for astrophysical calculation to a given excited state
c rateastro   : thermonuclear reaction rate factor
c rateastroex : thermonuclear reaction rate factor to a given excited state
c macsastro   : Maxwellian-averaged thermonuclear reaction cross section
c partf       : integrated partition function
c rateastrofis: thermonuclear reaction rate factor for fission
c macsastro   : thermonuclear reaction cross section
c macsastroex : thermonuclear reaction cross section to a given excited state
c macsastrofis: thermonuclear reaction cross section for fission
c xsastrofis  : astrophysical fission cross section
c
      maxZastro=numZastro
      maxNastro=numNastro
      do 10 in=0,numNastro
        do 20 iz=0,numZastro
          do 30 nen=1,numinc
            xsastro(iz,in,nen)=0.
            do 40 nex=0,numlev
              xsastroex(iz,in,nen,nex)=0.
   40       continue
   30     continue
          do 50 i=1,numT
            rateastro(iz,in,i)=0.
            macsastro(iz,in,i)=0.
            do 60 nex=0,numlev
              rateastroex(iz,in,i,nex)=0.
              macsastroex(iz,in,i,nex)=0.
   60       continue
   50     continue
   20   continue
   10 continue
      do 110 i=1,numT
        partf(i)=0.
        rateastrofis(i)=0.
        rateastroracap(i)=0.
        macsastrofis(i)=0.
        macsastroracap(i)=0.
  110 continue
      do 120 nen=1,numenin
        xsastrofis(nen)=0.
  120 continue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

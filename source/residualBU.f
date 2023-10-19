      subroutine residualBU
c
c +---------------------------------------------------------------------
c | Author: Marilena Avrigeanu
c | Date  : January 14, 2021
c | Task  : Residual production cross sections for break-up
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Acomp,Zcomp,Ncomp,nex
      real    branch(0:numlev)
c
c ************************ Cross sections ******************************
c
c Zcomp      : charge number index
c Ncomp      : neutron number index
c xspopnuc   : population cross section per nucleus
c Nlast      : last discrete level
c tau        : lifetime of state in seconds
c branch   : branching ratio for isomeric cross section
c xspopex    : population cross section summed over spin and parity
c
      do 10 Acomp=0,maxA
        do 20 Zcomp=0,maxZ
          Ncomp=Acomp-Zcomp
          if (Ncomp.lt.0.or.Ncomp.gt.maxN) goto 20
          if (xspopnuc(Zcomp,Ncomp).eq.0.) goto 20
          do 30 nex=0,Nlast(Zcomp,Ncomp,0)
            branch(nex)=0.
            if (nex.eq.0.or.tau(Zcomp,Ncomp,nex).ne.0.) then
              branch(nex)=xspopex(Zcomp,Ncomp,nex)/xspopnuc(Zcomp,Ncomp)
            endif
   30     continue
          xspopnuc(Zcomp,Ncomp)=xspopnuc(Zcomp,Ncomp)+
     +      xsBFnuc(Zcomp,Ncomp)  
          do 40 nex=0,Nlast(Zcomp,Ncomp,0)
            if (nex.eq.0.or.tau(Zcomp,Ncomp,nex).ne.0.) then
              xspopex(Zcomp,Ncomp,nex)=branch(nex)*xspopnuc(Zcomp,Ncomp)
            endif
   40     continue
   20   continue
   10 continue
      return
      end
Copyright (C)  2021 A.J. Koning, S. Hilaire and S. Goriely

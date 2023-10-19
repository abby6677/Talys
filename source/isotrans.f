      subroutine isotrans(Z,N)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 17, 2020
c | Task  : Correction factors for isospin forbidden transitions
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical flagiso
      integer Z,N,type
      real    ff(-1:numpar)
c
c ******** Correction factors for isospin forbidden transitions ********
c
c Z      : charge number of nucleus
c N      : neutron number of nucleus
c Zinit  : charge number of initial compound nucleus
c Ninit  : neutron number of initial compound nucleus
c k0     : index for incident particle
c primary: flag to designate primary (binary) reaction
c fiso   : correction factor for isospin forbidden transitions
c fisom  : correction factor for isospin forbidden transitions
c          in multiple emission
c
      do type=-1,6
        ff(type)=1.
      enddo
      flagiso=.false.
      if (Z.eq.N) then
        flagiso=.true.
        if (k0.eq.0) then
          ff(1)=2.
          ff(2)=2.
          ff(6)=5.
        endif
        if (k0.eq.1) ff(0)=2.
        if (k0.eq.2) ff(0)=2.
        if (k0.eq.6) ff(0)=5.
      endif
      if (Z.eq.N-1.or.Z.eq.N+1) then
        flagiso=.true.
        if (k0.eq.0) then
          ff(1)=1.5
          ff(2)=1.5
          ff(6)=1.5
        endif
        if (k0.eq.1) ff(0)=1.5
        if (k0.eq.2) ff(0)=1.5
        if (k0.eq.6) ff(0)=1.5
      endif
      if (flagiso) then
        if (primary) then
          do type=-1,6
            if (fiso(type).eq.-1.) fiso(type)=ff(type)
          enddo
        else
          do type=-1,6
            if (fisom(type).eq.-1.) fisom(type)=ff(type)
          enddo
        endif
      else
        if (primary) then
          do type=-1,6
            if (fiso(type).eq.-1.) fiso(type)=1.
          enddo
        else
          do type=-1,6
            if (fisom(type).eq.-1.) fisom(type)=1.
          enddo
        endif
      endif
      end
Copyright (C)  2020 A.J. Koning, S. Hilaire and S. Goriely

      function twkbphaseint(efis,ibar,Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Guillaume Scamps
c | Date  : August 12, 2019
c | Task  : Interpolation of direct WKB penetrability
c +---------------------------------------------------------------------
c
c       
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ibar,Zix,Nix,nen
      real    twkbphaseint,efis,Ea,Eb,Ta,Tb,Tf,Ewkb(0:nbinswkb)
c
c ************************* Interpolation ******************************
c
c Ewkb : energy
c Ta   : transmission coefficient
c Tb   : transmission coefficient
c
      do nen=0,nbinswkb
        Ewkb(nen)=Uwkb(Zix,Nix,nen)
      enddo
      if (efis.gt.Ewkb(nbinswkb)) then
        twkbphaseint=1.
      else
        call locate(Ewkb,0,nbinswkb,efis,nen)
        Ea=Ewkb(nen)
        Eb=Ewkb(nen+1)
        Ta=Twkbphase(Zix,Nix,nen,ibar)
        Tb=Twkbphase(Zix,Nix,nen+1,ibar)
        call pol1(Ea,Eb,Ta,Tb,efis,Tf)
        twkbphaseint=Tf
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

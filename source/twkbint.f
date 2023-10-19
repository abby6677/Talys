      function twkbint(efis,ibar,Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Guillaume Scamps
c | Date  : October 11, 2009
c | Task  : Interpolation of WKB penetrability
c +---------------------------------------------------------------------
c
c      SCAMPS: Modification in order to do logarithmic interpolation
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ibar,Zix,Nix,nen
      real    twkbint,efis,Ea,Eb,Ta,Tb,Tf,Ewkb(0:nbinswkb)
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
        twkbint=1.
      else
        call locate(Ewkb,0,nbinswkb,efis,nen)
        Ea=Ewkb(nen)
        Eb=Ewkb(nen+1)
        Ta=Twkb(Zix,Nix,nen,ibar)
        Tb=Twkb(Zix,Nix,nen+1,ibar)
        If (Ta.gt.0.and.Tb.gt.0) then
          Ta=log(Ta)
          Tb=log(Tb)
          call pol1(Ea,Eb,Ta,Tb,efis,Tf)
          twkbint=exp(Tf)
        else
          call pol1(Ea,Eb,Ta,Tb,efis,Tf)
          twkbint=Tf
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

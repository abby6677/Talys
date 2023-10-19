      function twkbtransint(efis,ibar,Zix,Nix)
c
c +---------------------------------------------------------------------
c | Author: Guillaume Scamps
c | Date  : August 13, 2019
c | Task  : Interpolation of direct WKB penetrability
c +---------------------------------------------------------------------
c       
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer ibar,Zix,Nix,nen
      real    Twkbtransint,efis,Ea,Eb,Ta,Tb,Tf,Ewkb(0:nbinswkb)
c
c ************************* Interpolation ******************************
c
c Ewkb : energy
c Ta   : transmission coefficient
c Tb   : transmission coefficient
c
      if (ibar.gt.2) stop 'pb tdirwkint ' 
      do nen=0,nbinswkb
        Ewkb(nen)=Uwkb(Zix,Nix,nen)
      enddo
      if (efis.gt.Ewkb(nbinswkb)) then
        Twkbtransint=1.
      else
        call locate(Ewkb,0,nbinswkb,efis,nen)
        Ea=Ewkb(nen)
        Eb=Ewkb(nen+1)
        Ta=Twkbtrans(Zix,Nix,nen,ibar)
        Tb=Twkbtrans(Zix,Nix,nen+1,ibar)
        if ( Ta.gt.0. and. Tb.gt.0 ) then
          Ta=log(Ta)
          Tb=log(Tb)
          call pol1(Ea,Eb,Ta,Tb,efis,Tf)
          Twkbtransint=exp(Tf)
        else
          call pol1(Ea,Eb,Ta,Tb,efis,Tf)
          Twkbtransint=Tf
        endif
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

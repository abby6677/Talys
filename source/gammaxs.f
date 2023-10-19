      function gammaxs(Zcomp,Ncomp,Egamma)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 15, 2018
c | Task  : Gamma ray cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zcomp,Ncomp,irad,l
      real    gammaxs,Egamma,xsgdr,fstrength,xsqd,quasideuteron
c
c ************** Calculate photo-absorption cross section **************
c
c Zcomp        : charge number index for compound nucleus
c Ncomp        : neutron number index for compound nucleus
c Einc         : incident energy in MeV
c Egamma       : gamma energy
c xsgdr        : photo-absorption cross section from GDR part
c fstrength    : gamma ray strength function
c kgr          : constant for gamma-ray strength function
c xsqd         : photo-absorption cross section from QD part
c quasideuteron: function for quasi-deuteron cross section
c
c 1. GDR part
c
      xsgdr=0.
      do 10 irad=0,1
        do 10 l=1,gammax
          xsgdr=xsgdr+fstrength(Zcomp,Ncomp,Einc,Egamma,irad,l)/
     +      kgr(l)*Egamma
   10 continue
c
c 2. QD part
c
      xsqd=quasideuteron(Egamma)
c
c Total absorption cross section
c
      gammaxs=xsgdr+xsqd
      return
      end
Copyright (C)  2017 A.J. Koning, S. Hilaire and S. Goriely

      function kalbachBU(type,Ein,ang)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 12, 2016
c | Task  : Kalbach systematics for break-up
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Ap,Npr,Ab,Za,Zb
      real    kalbachBU,Ein,ang,K,term,abu,K1,K2,K3,ang0,wang,Tc
c
c ************************ Kalbach systematics *************************
c
c kalbachBU : Kalbach function for break-up
c Ein       : incident energy
c ang       : angle
c K         : constant of Kalbach systematics
c Za        : charge of projectile
c Npr       : neutron number of projectile
c Zb        : charge of ejectile
c K1        : constant of Kalbach systematics
c K2        : constant of Kalbach systematics
c K3        : constant of Kalbach systematics
c parA      : mass number of particle
c k0        : index of incident particle
c parZ      : charge number of particle
c parN      : neutron number of particle
c term      : help variable
c abu       : Kalbach parameter
c Deff      : effective target-projectile separation
c Sab       : separation energy for projectile
c Ecent     : centroid energy for emission spectrum
c ang0      : angular barrier
c wang      : angular width
c Tc        : function for Coulomb dip
c
c  Systematics of Kalbach: Phys. Rev. Cxx, xxxxx, (2017)
c
      kalbachBU=1./fourpi
      K=1.4
      Ap=parA(k0)
      Npr=parN(k0)
      Za=parZ(k0)
      Ab=parA(type)
      Zb=parZ(type)
      abu=4.7+Ab
      if (Ab.eq.Ap-1) then
        if (Sab.gt.0.) then
          term=1.+exp((12.*Sab-Ecent)/(0.84*Sab))
          abu=4.*Ab+Zb-2.+0.029*Ecent+7.6/Ap/term
        endif
      endif
      K1=1.8
      K2=0.26
      K3=0.035
      if (Npr.gt.0.) then
        term=1+exp((K2-Ca/Ein)/(K3*Za/Npr))
        ang0=K1/(Ap**1.5)/term
        wang=min(0.09,ang0/3.)
        term=1.+exp((ang0-ang)/wang)
        Tc=1./term
        kalbachBU=(abu*abu+1.)/twopi*exp(-abu*cos(ang))*Tc
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

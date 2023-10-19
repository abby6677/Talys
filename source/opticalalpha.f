      subroutine opticalalpha(Zix,Nix,eopt)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning, Vlad and Marilena Avrigeanu
c | Date  : December 12, 2014
c | Task  : Other optical potential for alphas
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer Zix,Nix,Z,A
      real    eopt,A13,rb,e1opt,e2opt,e3opt,e4opt
c
c Alternative options for alpha OMP
c
c S. Goriely: inclusion of the alpha OMP of Mc Fadden & Satchler
c for alphaomp=2, folding model for alphaomp=3,4,5
c alphaomp=6 --> Avrigeanu et al. potential [PRC90,044612(2014)]
c
c Overwrite some of the previous values.
c
c v,rv,...: optical model parameters
c A13: A**1/3
c
      call ompadjust(eopt,6)
      Z=ZZ(Zix,Nix,0)
      A=AA(Zix,Nix,0)
      A13=real(A)**onethird
c
c alphaomp=2 --> Global OMP of Mc Fadden & Satchler,
c Nucl. Phys. A 84, 177 (1966).
c
      if (alphaomp.eq.2) then
        v=185.0
        rv=1.40
        av=0.52
        w=25.0
        rw=rv
        aw=av
        vd=0.
        wd=0.
        vso=0.
        wso=0.
        rc=1.3
      endif
c
c alphaomp=6 --> Global OMP of Avrigeanu et al.,
c Phys. Rev. C 90, 044612 (2014).
c
c e1opt: energy for alpha OMP
c e2opt: energy for alpha OMP
c e3opt: energy for alpha OMP
c e4opt: energy for alpha OMP
c
      if (alphaomp.eq.6) then
        rb=2.66+1.36*A13!Norenberg, HIC, N-H, 1980, p.8
        e2opt=(2.59+10.4/A)*Z/rb
        e1opt=-3.03-0.76*A13+1.24*e2opt
        e3opt=22.2+0.181*Z/A13
        e4opt=29.1-0.22*Z/A13
        if (eopt.le.e3opt) then
          v=165.0 + 0.733*Z/A13 - 2.64*eopt
        else
          v=116.5 + 0.337*Z/A13 - 0.453*eopt
        endif
        v=max(v,-100.)
        if (eopt.le.25.0) then
          rv=1.18 + 0.012*eopt
        else
          rv=1.48
        endif
        if (eopt.le.e2opt) then
          av=0.631+(0.016-0.001*e2opt)*Z/A13
        elseif (eopt.le.e4opt) then
          av=0.631+(0.016-0.001*eopt)*Z/A13
        else
          av=0.684-0.016*Z/A13-(0.0026-0.00026*Z/A13)*eopt
        endif
        av=max(av,0.1)
c
        w=amax1(2.73-2.88*A13+1.11*eopt,0.0)
        rw=1.34
        aw=0.50
c
        vd=0.
        if (eopt.le.e1opt) then
          wd=4.0
        elseif (eopt.le.e2opt) then
          wd=22.2+4.57*A13-7.446*e2opt+6.0*eopt
        else
          wd=22.2+4.57*A13-1.446*eopt
        endif
        wd=max(wd,0.)
        if (A.le.152.or.A.ge.190) then
          rwd=1.52
        else
          rwd=amax1(1.74-0.01*eopt,1.52)
        endif
        awd=0.729-0.074*A13
c
        vso=0.
        wso=0.
        rc=1.3
      endif
c
c alphaomp=7 --> Global OMP of Nolte et al.,
c Phys. Rev. C 36, 1312 (1987).
c
      if (alphaomp.eq.7) then
        v=101.1+6.051*Z/A13-0.248*eopt
        rv=1.245
        av=0.817-0.0085*A13
        w=26.82-1.706*A13+0.006*eopt
        rw=1.57
        aw=0.692-0.02*A13
        vd=0.
        wd=0.
        vso=0.
        wso=0.
        rc=1.3
      endif
c
c alphaomp=8 --> Global OMP of Avrigeanu et al.,
c Phys. Rev. C 49, 2136 (1994).
c
      if (alphaomp.eq.8) then
        v=101.1+6.051*Z/A13-0.248*eopt
        rv=1.245
        av=0.817-0.0085*A13
        if (eopt.le.73.) then
          w=12.64-1.706*A13+0.20*eopt
        else
          w=26.82-1.706*A13+0.006*eopt
        endif
        rw=1.57
        aw=0.692-0.02*A13
        vd=0.
        wd=0.
        vso=0.
        wso=0.
        rc=1.3
      endif
c
c Possible adjustment of parameters
c
c Fv1,....     : adjustable factors for OMP (default 1.)
c
      v=Fv1*v
      rv=Frv*rv
      av=Fav*av
      w=Fw1*w
      rw=Frw*rw
      aw=Faw*aw
      vd=Fd1*vd
      rvd=Frvd*rvd
      avd=Favd*avd
      wd=Fd1*wd
      rwd=Frwd*rwd
      awd=Fawd*awd
      vso=Fvso1*vso
      rvso=Frvso*rvso
      avso=Favso*avso
      wso=Fwso1*wso
      rwso=Frwso*rwso
      awso=Fawso*awso
      rc=Frc*rc
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

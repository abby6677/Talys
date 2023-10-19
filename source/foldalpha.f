      subroutine foldalpha(Zix,Nix,E)
c
c Global alpha potential = real double folding potential + imaginary WS potential
c   dispersive relations also included if alphaomp=5
c Ref: Demetriou, Grama and Goriely (2002) NPA707, 253
c
      include "talys.cmb"
      character*80 key
      real rhomom(numjlm),radmom(numjlm)
      real va(1000),rva(1000)
      real z,a,ee,rb,e2exp,factor1,factor2
      real vopr(numjlm),wopr(numjlm)
      real alphav(10)
      real e,a2,a3,a13,as,es,aaw,raw,te,expj0,ww,dwv,dws,
     &  efermia,rww,aww,www,ws,rws,aws,wws,
     &  vdiv,vdis,am1,am2,am3,rr,c,hh,aa1,bb1,v5d,eref,dvolj
      integer Zix,Nix,i,k,nu,nradrho,klo,khi,kk
c
c ************************ Alpha Optical potential *************************
c
c Zix        : charge  number index of target nucleus
c Nix        : neutron number index of target nucleus
c E          : incident energy
c Z          : charge number of target nucleus
c A          : mass   number of target nucleus
c a13        : mass**1/3
c a3         : mass**3
c e2exp      : lowest inelastic threshold
c radmom     : radial grid for the potential, identical as the one used for JLM
c rhomom     : total density (neutron + proton) at a given radius radmom
c alphaomp   : model for the alpha OMP
c            : 3=OMP I, 4=OMP II, 5=OMP III of Demetriou et al. (2002)
c vopr       : Final real part of the alpha optical potential
c wopr       : Final imaginary part of the alpha optical potential
c
      a=AA(Zix,Nix,0)
      z=ZZ(Zix,Nix,0)
      ee=E
      do i=1,numjlm
        radmom(i)=radjlm(Zix,Nix,i)
        rhomom(i)=rhojlmn(Zix,Nix,i,1)+rhojlmp(Zix,Nix,i,1)
      enddo
      e2exp=edis(Zix,Nix,1)
      a2=a**2
      a3=a**3
      a13=a**(1./3.)
c
c
c--------------------------------------------------------------------------------------
c alphav: optical model part
c alphav =  1: reduced radius r0=r/a13 of real pot
c           2: diffuseness a of woods-saxon
c           3: reduced radius r0s=r/a13 of imag. pot.
c           4: diffuseness as of imag. pot.
c           5: fraction of surface-peaked im. pot. relative to volume part
c           6: depth vso of spin-orbit pot.
c           7: depth v of real pot.
c           8: depth w of imag. pot.
c rww        : radius of imaginary volume potential
c aww        : diffuseness of imaginary volume potential
c www        : depth of imaginary volume potential
c rws        : radius of imaginary surface potential
c aws        : diffuseness of imaginary surface potential
c wws        : depth of imaginary surface potential
c te         : energy-dependent term of the imaginary volume integral
c as         : diffuseness of the functional form of the imaginary volume integral
c es         : energy threshold of the functional form of the imaginary volume integral
c expj0      : saturation value of the imaginary volume integral
c--------------------------------------------------------------------------------------
c
      v5d=0.
      goto (10,20,30,200) alphaomp-2
   10 continue
c--------------------------------------------------------------------------------------
c Case alphaomp=3
c Imaginary potential OMP I of Demetriou et al. (2002)
c Woods-Saxon E-dependent volume imaginary potential
c--------------------------------------------------------------------------------------
      alphav(1)=1.25
      raw=0.85281+0.02202*a-2.14551e-4*a2+7.92942e-7*a3-
     1     9.94658e-10*a2*a2
      aaw=-0.13526+0.02029*a-1.98441e-4*a2+7.35104e-7*a3-
     2     9.15272e-10*a2*a2
      if (aaw.le.1.e-2) aaw=1.e-2
      as=7.65867-7.5669*e2exp+2.50486*e2exp**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(ee-es)/as))
      expj0=77.
      if (a.le.90.) expj0=135.-0.644*a
      ww=3./pi/raw**3*expj0*te
      alphav(3)=raw
      alphav(4)=aaw
      alphav(8)=ww
      alphav(5)=0.
      dwv=0.
      dws=0.
      goto 50
c
   20 continue
c--------------------------------------------------------------------------------------
c Case alphaomp=4
c Imaginary potential OMP II of Demetriou et al. (2002)
c Woods-Saxon E-dependent volume+surface imaginary potential
c--------------------------------------------------------------------------------------
      alphav(1)=1.25
      raw=1.47385-0.00134615*a
      aaw=0.29
      alphav(5)=0.9
      as=7.65867-7.5669*e2exp+2.50486*e2exp**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(ee-es)/as))
      expj0=77.
      if(a.le.90.) expj0=135.-0.644*a
      ww=expj0*te/(raw**3/3.+7.603*alphav(5)*aaw*raw**2/a13)/pi
c     ww=expj0*te/((1.-alphav(5))*raw**3/3.+7.603*alphav(5)*aaw*
c    &  raw**2/a13)/pi
      alphav(3)=raw
      alphav(4)=aaw
      alphav(8)=ww
      dwv=0.
      dws=0.
      goto 50
c
   30 continue
c--------------------------------------------------------------------------------------------------
c Case alphaomp=5
c Imaginary potential OMP III of Demetriou et al. (2002)
c Woods-Saxon E-dependent volume+surface imaginary potential plus Dispersion contribution
c  using prescription and functions of Capote et al. J. Phys. G 27 (2001) B15
c
c efermia    : Fermi energy
c am1        : mass excess of the target-projectile system
c am2        : mass excess of the target+projectile system
c am3        : mass excess of the projectile
c dwv        : dispersive contribution to the real potential arising from the volume imaginary pot
c dws        : dispersive contribution to the real potential arising from the surface imaginary pot
c--------------------------------------------------------------------------------------------------
c
c v5d: OMP component
c
      if (expmexc(Zix-2,Nix-2).ne.0..and.expmexc(Zix+2,Nix+2).ne.0.)
     &  then
        am1=expmexc(Zix+2,Nix+2)
        am2=expmexc(Zix-2,Nix-2)
      else
        am1=thmexc(Zix+2,Nix+2)
        am2=thmexc(Zix-2,Nix-2)
      endif
      am3=(parmass(6)-4.)*amu
      efermia=-0.5*(am1-am2+2.*am3)
c
      alphav(1)=1.25
      raw=1.47385-0.00134615*a
      aaw=0.30
      alphav(5)=0.35
      v5d=alphav(5)
      c=0.005
      if(ee.lt.13.) c=-0.165*ee+2.15
      alphav(5)=alphav(5)*exp(-c*abs(ee-efermia))
c
      as=7.65867-7.5669*e2exp+2.50486*e2exp**2
      es=0.0854*a+1.1307*as
      te=1./(1.+exp(-(ee-es)/as))
      expj0=75.
      if(a.le.96.) expj0=135.-0.644*a
c volume integral
      ww=expj0*te/((1.-alphav(5))*raw**3/3.+7.603*alphav(5)*aaw*
     &  raw**2/a13)/pi
      alphav(3)=raw
      alphav(4)=aaw
      alphav(8)=ww
c
   50 continue

c final radius, diffuseness and depth of the imaginary WS-type potential
c They are multiplied by the OMP adjustment keywords of TALYS
c (default 1.)
c   volume term
      rww=rvadjust(6)*alphav(3)*a13
      aww=avadjust(6)*alphav(4)
      if (alphaomp.eq.4) then
        www=w1adjust(6)*alphav(8)
      else
        www=w1adjust(6)*alphav(8)*(1.-alphav(5))
      endif
c   surface term
      rws=rwdadjust(6)*1.09*rww
      aws=awdadjust(6)*1.6*aww
      wws=d1adjust(6)*alphav(8)*alphav(5)
c
c dispersive contributions: dwv = volume; dws = surface; dvolj = volume integral
c new method for dispersive relation following Mahaux, Ngo and Satchler, NPA 449 (1986) 354
c
c raw: radius
c dvolj: volume integral
c eref: reference energy
c
      if (alphaomp.eq.5) then
        eref=150.
        call mahaux(a,ee,eref,expj0,as,es,v5d,raw,aaw,dwv,dws,dvolj,
     &    efermia)
      endif
c
c--------------------------------------------------------------------------------------------------
c determination of the real folding potential through the product of the fourier transforms
c rb         : maximum radius value
c nu         : number of radial grid point
c rva        : radial coordinate for the double folding potential
c va         : double folding potential
c nradrho    : number of grid point used in radmom and rhomom
c vdiv       : dispersive contribution to the real potential arising from the volume imaginary pot
c vdis       : dispersive contribution to the real potential arising from the surface imaginary pot
c ws         : surface component of the imaginary potential
c rc         : coulomb radius based on elton's formula
c--------------------------------------------------------------------------------------------------
      rb=radjlm(Zix,Nix,numjlm)
      nu=1000
      nradrho=numjlm

      call afold(z,a,ee,rb,nu,rva,va,radmom,rhomom,nradrho)
c
c Depth and shape are multiplied by the OMP adjustment keywords of TALYS
c (default 1.)
c
c ompadjustp: flag for local optical model parameter adjustment
c adjust    : subroutine for energy-dependent parameter adjustment
c factor1,2 : multiplication factor
c
      if (ompadjustp(6)) then
        key='aradialcor'
        call adjust(E,key,Zix,Nix,0,0,factor1)
        key='adepthcor'
        call adjust(E,key,Zix,Nix,0,0,factor2)
      else
        factor1=1.
        factor2=1.
      endif
      do i=1,nu
c       rva(i)=factor1*aradialcor*rva(i)
c
csg Correction of the radius dependence of the real part after an analysis of the 
c   of the (a,g) and (a,n) data of deformed nuclei : 27/4/2018 (Brussels)
c   increase of rva by 3% for deformed nuclei but only below typically 18 MeV

        rva(i)=factor1*aradialcor*rva(i)*(1.+beta2(Zix,Nix,0)/15.
     &    /(1.+exp((E-18.)/2.)))
        va(i)=factor2*adepthcor*va(i)
      enddo
c
c interpolation of the alpha optical potential va(rva) on the given radial grid radjlm
c
c aa1: help variable
c bb1: help variable
c aaw: help variable
c hh : radius difference
c
      do 100 k=1,numjlm
        rr=radjlm(Zix,Nix,k)
        klo=1
        khi=nu
        if (rr.le.rva(klo)) then
           khi=klo+1
           goto 110
        endif
        if (rr.ge.rva(khi)) then
           klo=khi-1
           goto 110
        endif
  120   if (khi-klo.gt.1) then
          kk=(khi+klo)/2.
          if (rva(kk).gt.rr) then
            khi=kk
          else
          klo=kk
          endif
          goto 120
        endif
  110   hh=rva(khi)-rva(klo)
        aa1=(rva(khi)-rr)/hh
        bb1=(rr-rva(klo))/hh
        vopr(k)=aa1*va(klo)+bb1*va(khi)
c dispersive contributions to real potential for alphaomp=5  (OMP III)
        if (alphaomp.lt.5) then
          vdiv=0.
          vdis=0.
        else
          if (abs((rr-rww)/aww).lt.88.) then
            vdiv=-dwv/(1.+exp((rr-rww)/aww))
          else
            vdiv=0.
          endif
          if (abs((rr-rws)/aws).lt.88.) then
            vdis=-4.*dws*exp((rr-rws)/aws)/(1.+exp((rr-rws)/aws))**2
          else
            vdis=0.
          endif
        endif
        vopr(k)=vopr(k)+vdiv+vdis
c final imaginary potential
        wopr(k)=0.
        if (abs((rr-rww)/aww).lt.88.)
     &    wopr(k)=-www/(1.+exp((rr-rww)/aww))
        if (abs((rr-rws)/aws).lt.88.) then
          ws=-4.*wws*exp((rr-rws)/aws)/(1.+exp((rr-rws)/aws))**2
        else
          ws=0.
        endif
        wopr(k)=wopr(k)+ws
c-----------------------------------------------------------------
c potjlm   1 : final real central potential
c potjlm   2 : final imaginary central potential
c-----------------------------------------------------------------
        potjlm(Zix,Nix,k,1)=vopr(k)
        potjlm(Zix,Nix,k,2)=wopr(k)
        potjlm(Zix,Nix,k,3)=0.
        potjlm(Zix,Nix,k,4)=0.
        potjlm(Zix,Nix,k,5)=0.
        potjlm(Zix,Nix,k,6)=0.
  100 continue
c calculate the coulomb radius based on elton's formula
      rc=1.123+2.352*(a**(-.666666))-2.07*(a**(-1.333333))
  200 continue
      return
      end

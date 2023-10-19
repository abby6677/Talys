      subroutine mom(Zix,Nix,pZ,e)
c
c +---------------------------------------------------------------------
c | Author: Eric Bauge and Arjan Koning
c | Date  : April 4, 2012
c | Task  : Microscopic JLM OMP
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*80 key
      integer      Zix,Nix,Z,N,A,i,j
      double precision       pZ,e,lv,lw,lv1,lw1,lvso,lwso
      real         factor,rhomomn(numjlm,6),rhomomp(numjlm,6),
     +             vpot(numjlm,6),rcjlm,alam
c
c *********************** Parameterization *****************************
c
c pZ  : product of charges
c ZZ,Z: charge number of residual nucleus
c NN,N: neutron number of residual nucleus
c AA,A: mass number of residual nucleus
c
c Parameters from prc 63, 024607 2001
c lv: real potential depth normalization factor
c lw: imag potential depth normalization factor
c lv1: real isovector potential depth normalization factor
c lw1: imag isovector potential depth normalization factor
c lvso: real spin orbit potential depth normalization factor
c lwso: imaginary spin orbit potential depth normalization factor
c ompadjustp: logical for energy-dependent OMP adjustment
c adjust    : subroutine for energy-dependent parameter adjustment
c factor    : multiplication factor
c lvadjust   : adjustable parameter
c alam    : adjustable parameter
c potjlm : JLM potential depth values
c normjlm: JLM potential normalization factors
c
      A=AA(Zix,Nix,0)
      Z=ZZ(Zix,Nix,0)
      N=A-Z
      lv=.0008*log(e*1000.)+.00018*(log(e*1000.))**2.+.951
      lw=(1.24-(1.0/(1.+exp(((e)-4.5)/2.9))))
     &     *(1.+.06*exp(-((e-14.)/3.7)**2 ))
     &     *(1.-0.09*exp(-((e-80.)/78.)**2))
      if(e.gt.80.)lw=lw*(1+(e-80.)/(400.))
      lv1=1.5-(0.65/(1.+exp((e-1.3)/3.0)))
c
c Possible modification of JLMB potential with jlmmode
c  cf Goriely & Delaroche: PLB653, 158 (2007)
c
      alam=0.44
      if (jlmmode.eq.1) alam=1.10*exp(-0.4*e**0.25)
      if (jlmmode.ge.2) alam=1.375*exp(-0.2*e**0.5)
      if (jlmmode.eq.3) lw=lw*2.
c
      lw1=(1.1+(alam/(1.+(exp(((e)-40.)/50.9))**4 )))
     &     *(1.-.065*exp(-((40.-e)/13.)**2))
     &     *(1.-.083*exp(-((200.-e)/80.)**2))
c
      lvso=40.+exp(-e*0.013)*130.
      lwso=-0.2*(e-20)
      if (ompadjustp(1)) then
        key='lvadjust'
        call adjust(real(e),key,0,0,0,0,factor)
        lv=factor*lvadjust*lv
        key='lwadjust'
        call adjust(real(e),key,0,0,0,0,factor)
        lw=factor*lwadjust*lw
        key='lv1adjust'
        call adjust(real(e),key,0,0,0,0,factor)
        lv1=factor*lv1adjust*lv1
        key='lw1adjust'
        call adjust(real(e),key,0,0,0,0,factor)
        lw1=factor*lw1adjust*lw1
        key='lvsoadjust'
        call adjust(real(e),key,0,0,0,0,factor)
        lvso=factor*lvsoadjust*lvso
        key='lwsoadjust'
        call adjust(real(e),key,0,0,0,0,factor)
        lwso=factor*lwsoadjust*lwso
      endif
c
c Calculate and write the potentials in ecis format
c
c rcjlm: coulomb radius for JLM
c rhomomn    : total neutron density at a given radius radmom
c rhomomp    : total proton density at a given radius radmom
c vpot       : optical potential
c
      do i=1,numjlm
        do j=1,6
          rhomomn(i,j)=rhojlmn(Zix,Nix,i,j)
          rhomomp(i,j)=rhojlmp(Zix,Nix,i,j)
        enddo
      enddo
      call momjlmecis(Z,N,pZ,e,lv1,lw1,1.25d0,1.35d0,
     +  rhomomn,rhomomp,vpot,rcjlm)
      do i=1,numjlm
        do j=1,6
          potjlm(Zix,Nix,i,j)=min(vpot(i,j),1.e30)
          potjlm(Zix,Nix,i,j)=max(potjlm(Zix,Nix,i,j),-1.e30)
        enddo
      enddo
      normjlm(Zix,Nix,1)=lv
      normjlm(Zix,Nix,2)=lw
      normjlm(Zix,Nix,3)=lv1
      normjlm(Zix,Nix,4)=lw1
      normjlm(Zix,Nix,5)=lvso
      normjlm(Zix,Nix,6)=lwso
      rc=rcjlm
      return
      end

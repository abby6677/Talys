      subroutine tdirbarrier(Zcomp,Ncomp,J2,parity,ibar,ibar2,trfis,
     +  rhof,Eex)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Marieke Duijvestijn
c | Date  : October 31, 2019
c | Task  : Fission transmission coefficient for one barrier
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,J2,parity,ibar,J,itr,j2trans,
     +                 pitrans,i,ibar2,ibar3
      real             Eex,etrans,Eeff,twkbint,elow,
     +                 emid,eup,dE1,dE2,Twkbphaseint
      double precision trfis,rhof,trfisone,rho1,rho2,rho3,r1log,r2log,
     +                 r3log,rho,trfisonetwo,trfistwo,trfisthree
      external Twkbphaseint
c
c ********** Fission transmission coefficient for one barrier **********
c
c Zcomp : charge number index for compound nucleus
c Ncomp : neutron number index for compound nucleus
c J2    : 2 * J
c parity: parity     
c ibar  : fission barrier
c trfis : transmission coefficient
c rhof  : integrated level density
c Eex   : excitation energy of entrance bin
c
      J=J2/2
      trfis=0.
      rhof=0.
c
c Correct LDM barrier height with ground state shell correction
c
c fismodelx: fission model
c nfisbar  : number of fission barrier parameters
c fbarrier : height of fission barrier
c deltaW   : shell correction in nuclear mass
c fwidth   : width of fission barrier
c
      if (.not.flagfispartdamp) return
c
c 1. Discrete states
c
c nfistrrot     : number of rotational transition states for barrier
c etrans,j2trans: help variables
c efistrrot     : energy of rotational transition states
c jfistrrot     : spin of rotational transition states
c pitrans,Eeff  : help variables
c pfistrrot     : parity of rotational transition states
c primary       : flag to designate primary (binary) reaction
c trfisone      : help variable
c twkbint       : WKB penetrability
c numhill       : maximum number of Hill-Wheeler points
c tfisA         : transmission coefficient for Hill-Wheeler magnitude
c rhofisA       : integrated level density corresponding to tfisA
c
      if (abs(ibar-ibar2).eq.2) then
        ibar3=(ibar2+ibar)/2
      else
        ibar3=ibar
      endif
      do 10 itr=1,nfistrrot(Zcomp,Ncomp,ibar3)
        etrans=efistrrot(Zcomp,Ncomp,ibar3,itr)
        if (Eex.lt.etrans) goto 10
        j2trans=int(2.*jfistrrot(Zcomp,Ncomp,ibar3,itr))
        pitrans=pfistrrot(Zcomp,Ncomp,ibar3,itr)
        Eeff=Eex-etrans
        if ((J2.eq.j2trans).and.(parity.eq.pitrans)) then
          if (abs(ibar-ibar2).eq.1) then
            trfisone=twkbint(Eeff,ibar,Zcomp,Ncomp)
            trfistwo=twkbint(Eeff,ibar2,Zcomp,Ncomp)
            if ( twkbphaseint(Eeff,ibar,Zcomp,Ncomp) .gt. 0 ) then
              trfis=trfis+trfisone*trfistwo/
     +        (1+(1.-trfisone)*(1.-trfistwo))
            else
              trfis=trfis+trfisone*trfistwo
            endif
          elseif (abs(ibar-ibar2).eq.2) then
            trfisone=twkbint(Eeff,ibar,Zcomp,Ncomp)
              trfistwo=twkbint(Eeff,ibar2,Zcomp,Ncomp)
            ibar3=(ibar2+ibar)/2
              trfisthree=twkbint(Eeff,ibar3,Zcomp,Ncomp)
            if ( twkbphaseint(Eeff,ibar,Zcomp,Ncomp) .gt. 0 ) then
              trfisonetwo=trfisone*trfisthree/
     +            (1+(1.-trfisone)*(1.-trfisthree))
            else
              trfisonetwo=trfisone*trfisthree
            endif
            if ( twkbphaseint(Eeff,ibar3,Zcomp,Ncomp) .gt. 0 ) then
                trfis=trfis+trfisonetwo*trfistwo/
     +            (1+(1.-trfisonetwo)*(1.-trfistwo))
            else
              trfis=trfis+trfisonetwo*trfistwo
            endif
          endif
          rhof=rhof+1.
        endif
   10 continue
c
c 2. Continuum
c
c fecont       : start of continuum energy
c nbintfis     : number of integration bins
c elow,eup,emid: help variables
c eintfis      : excitation energy for fission   
c dE1,dE2      : help variables
c rho1-3       : help variables
c r1log,...    : help variables
c rhofis       : integrated level density  
c rho          : integrated level density
c
      if (abs(ibar-ibar2).eq.2) then
        ibar3=(ibar2+ibar)/2
      else
        ibar3=ibar
      endif
      if (Eex.ge.fecont(Zcomp,Ncomp,ibar3)) then
        do 20 i=1,nbintfis(ibar3)-2,2
          elow=eintfis(i,ibar3)
          if (elow.gt.Eex) goto 20
          emid=min(eintfis(i+1,ibar3),Eex)
          eup=min(eintfis(i+2,ibar3),Eex)
          dE1=emid-elow
          dE2=eup-emid
          if (abs(ibar-ibar2).eq.1) then       
            rho1=rhofis(i,J,parity,ibar3)*(1.+1.d-10)
            rho2=rhofis(i+1,J,parity,ibar3)
            rho3=rhofis(i+2,J,parity,ibar3)*(1.+1.d-10)
          else
            rho1=min(rhofis(i,J,parity,ibar)*(1.+1.d-10),
     +        rhofis(i,J,parity,ibar2)*(1.+1.d-10))
            rho2=min(rhofis(i+1,J,parity,ibar),
     +        rhofis(i+1,J,parity,ibar2))
            rho3=min(rhofis(i+2,J,parity,ibar)*(1.+1.d-10),
     +        rhofis(i+2,J,parity,ibar2)*(1.+1.d-10))
        endif
        r1log=log(rho1)
        r2log=log(rho2)
        r3log=log(rho3)
        if (r2log.ne.r1log.and.r2log.ne.r3log) then
          rho=(rho1-rho2)/(r1log-r2log)*dE1
     +      +(rho2-rho3)/(r2log-r3log)*dE2
        else
          rho=rho2*(dE1+dE2)
        endif 
        Eeff=Eex-emid
        if (abs(ibar-ibar2).eq.1) then       
          trfisone=twkbint(Eeff,ibar,Zcomp,Ncomp)
          trfistwo=twkbint(Eeff,ibar2,Zcomp,Ncomp)
          if ( twkbphaseint(Eeff,ibar,Zcomp,Ncomp) .gt. 0 ) then
            trfis=trfis+rho*trfisone*trfistwo/
     +        (1+(1.-trfisone)*(1.-trfistwo))
          else
            trfis=trfis+rho*trfisone*trfistwo
          endif      
        elseif (abs(ibar-ibar2).eq.2) then
          trfisone=twkbint(Eeff,ibar,Zcomp,Ncomp)
          trfistwo=twkbint(Eeff,ibar2,Zcomp,Ncomp)
          ibar3=(ibar2+ibar)/2
            trfisthree=twkbint(Eeff,ibar3,Zcomp,Ncomp)
          if ( twkbphaseint(Eeff,ibar,Zcomp,Ncomp) .gt. 0 ) then
            trfisonetwo=rho*trfisone*trfisthree/
     +        (1+(1.-trfisone)*(1.-trfisthree))
          else
            trfisonetwo=rho*trfisone*trfisthree
          endif
          if ( twkbphaseint(Eeff,ibar3,Zcomp,Ncomp) .gt. 0 ) then
            trfis=trfis+rho*trfisonetwo*trfistwo/
     +        (1+(1.-trfisonetwo)*(1.-trfistwo))
          else
            trfis=trfis+rho*trfisonetwo*trfistwo
          endif
        endif
          rhof=rhof+rho
   20   continue
      endif
      return
      end
Copyright (C) 2004  A.J. Koning, S. Hilaire and M.C. Duijvestijn

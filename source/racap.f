      subroutine racap
c
c +---------------------------------------------------------------------
c | Author: Stephane Goriely, Xu Yi, and Arjan Koning
c | Date  : January 6, 2017
c | Task  : Radiative capture model
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          i,J,par,parity,Zix,Nix,nex,iopt,prac(numlev2)
      real             jrac(numlev2),qq,erac(numlev2),
     +                 exfin(numlev2),jlmracap1(numjlm),vncap1,rvncap1,
     +                 avncap1,xsall(3),spfacst(numlev2),
     +                 xspex(numlev2),xsp(numlev2,0:numJph,2)
      double precision pZ
c
c ********************** General initializations ***********************
c
c iopt    : index
c jlmracap: JLMB potential for the radial grid: [0.1, 20], step of 0.1
c jlmracap1: JLMB potential for the radial grid: [0.1, 20], step of 0.1
c
c Calculation the initial wood saxon potential parameter
c
      iopt=racopt
      Zix=parZ(k0)
      Nix=parN(k0)
      call optical(Zix,Nix,k0,Einc)
      vncap1=v
      rvncap1=rv
      avncap1=av
c
c iopt=1 - Woods-Saxon potential
c iopt=2 - No.2, no use
c iopt=3 - JLMB, Jeukenne et al. (1977) and Burge PRC, 1998
c iopt=4 - Folding potential M3Y (now not available)
c
c Calculation the initial state JLMB potential
c
      pZ=dble(real(Zix)*ZTarget)
      call mom(Zix,Nix,pZ,dble(Einc))
      do i=1,numjlm
        jlmracap1(i)=normjlm(Zix,Nix,1)*
     &    (potjlm(Zix,Nix,i,1)+
     &    normjlm(Zix,Nix,3)*
     &    (rhojlmn(Zix,Nix,i,1)+
     &    rhojlmp(Zix,Nix,i,1))/
     &    (rhojlmn(Zix,Nix,i,1)-
     &    rhojlmp(Zix,Nix,i,1))
     &    *potjlm(Zix,Nix,i,3))
      enddo
c
c Set the experimental level position, spin and parity.
c
c erac: ebergy
c prac: parity
c
      do i=1,numlev2
        spfacst(i)=0.
      enddo
      do i=1,nlevexpracap
        jrac(i)=jdis(0,0,i)
        prac(i)=parlev(0,0,i)
        erac(i)=edis(0,0,i)
        spfacst(i)=spectfac(0,0,i-1)
      enddo
      do i= nlevexpracap,maxex(0,0)+1
        spfacst(i)=spectfac(0,0,i-1)
      enddo
      qq=min(S(0,0,1),S(0,0,2))
c
c End of initial setting and begin to racap calculation.
c
c rvncap1: radius
c avncap1: diffuseness
c vncap1: potential
c jrac: spin
c spfacst: spin
c xsall: cross section
c exfin: energy
c qq: help variable
c xspex: help variable
c
      call racapcalc(Einc,Ztarget,Atarget,beta2(Zix,Nix,0),
     &  jdis(Zix,Nix,0),parlev(Zix,Nix,0),gsspin(Zix,Nix),
     &  gsparity(Zix,Nix),expmass(Zix,Nix),thmass(Zix,Nix),
     &  Zix,parA(k0),parspin(k0),parmass(k0),jdis(0,0,0),parlev(0,0,0),
     &  gsspin(0,0),gsparity(0,0),expmass(0,0),thmass(0,0),qq,
     &  numlev2,nlevexpracap,erac,jrac,prac,exfin,
     &  numjlm,jlmracap1,jlmracap2,numJ,numdensracap,chglposj,chglnegj,
     &  vncap1,rvncap1,avncap1,vncap2,rvncap2,avncap2,xsall,xsracape,
     &  xspex,xsp,iopt,pi,e2,amu,hbarc,nlevracap(0,0),spfacst,ispect)
c
c Racap output cross section,transfer the unit of the cross section to millibarns.
c
      xsracape=xsracape*1000.0
      xsracap(nin)=xsracape
      xsracapEM(nin,1,1)=xsall(1)*1000.0
      xsracapEM(nin,1,2)=xsall(2)*1000.0
      xsracapEM(nin,0,1)=xsall(3)*1000.0
      xsracapedisc=0.
      xsracapecont=0.
      do i=1,nlevracap(0,0)
        if (i.le.nlevexpracap) then
c          spectfac(0,0,i-1)=spfacst(i)
          xsracapedisc=xsracapedisc+xspex(i)*1000.0
          xsracappopex(i-1)=xspex(i)*1000.0
          do J=0,numJ
          do parity=-1,1,2
            if (parity.eq.-1) par=1
            if (parity.eq.1) par=2
            xsracappop(i-1,J,parity)=xsp(i,J,par)*1000.0
          enddo
          enddo
        else
          xsracapecont=xsracapecont+xspex(i)*1000.0
          do nex=nlevexpracap,maxex(0,0)
c            spectfac(0,0,nex)=spfacst(i)
            if (exfin(i).lt.Ex(0,0,nex)+deltaEx(0,0,nex)/2..and.
     +        exfin(i).ge.Ex(0,0,nex)-deltaEx(0,0,nex)/2.) then
              xsracappopex(nex)=xsracappopex(nex)+xspex(i)*1000.0
              do J=0,numJ
              do parity=-1,1,2
                if (parity.eq.-1) par=1
                if (parity.eq.1) par=2
                xsracappop(nex,J,parity)=xsracappop(nex,J,parity)+
     +            xsp(i,J,par)*1000.0
              enddo
              enddo
            endif
          enddo
        endif
      enddo
      return
      end

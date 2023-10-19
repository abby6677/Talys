      subroutine rpevap
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 1, 2015
c | Task  : Evaporation of residual products
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer iz,ia,nex,type,nen,npar,in,ih,it,id,ip,
     +        ident,idc,iaa,inn
      real    xsexcpart(0:numpar,numin),sumpost,sumpfns,Esumpfns,E,Eav,
     +        maxwell,summax,dE,xsc,yZA,yA
c
c ********************** Loop over fission fragments *******************
c
c Do a full TALYS calculation for each residual product
c
      flagrpruns=.true.
c     do type=0,6
c       do nen=0,numen2
c         pfns(type,nen)=0.
c         maxpfns(type,nen)=0.
c       enddo
c     enddo
c     fiseps=Rfiseps*xsfistot
      write(*,'(/" ########## Start of loop over residual products"/)')
      do 20 ia=1,Atarget0
        do 30 iz=1,Ztarget0
          in=ia-iz
          if (in.lt.1.or.in.gt.numneu) goto 30
          if (xspopnuc0(iz,ia).lt.5.) goto 30
          Arp=ia
          Zrp=iz
          call evaptalys
c
c Add fission product cross sections
c
c         do 110 Zix=0,maxZ
c           do 120 Nix=0,maxN
c             xsfpZApost(iz,in)=xsfpZApost(iz,in)+xspopnuc(Zix,Nix)
c             xsfpApost(ia)=xsfpApost(ia)+xspopnuc(Zix,Nix)
c             do 130 nex=0,Nlast(Zix,Nix,0)
c               if (nex.eq.0.or.tau(Zix,Nix,nex).ne.0.)
c    +            xsfpex(iz,in,nex)=xsfpex(iz,in,nex)+
c    +            xspopex(Zix,Nix,nex)
c 130         continue
c 120       continue
c 110     continue
c
c Add prompt fission particle and gamma production and spectra
c
c         if (xsinitpop.gt.0.) then
c           do 202 npar=1,numin
c             do 204 type=0,6
c               xsexcpart(type,npar)=0.
c 204         continue
c             do 206 iaa=0,numia
c             do 206 ih=0,numih
c             do 206 it=0,numit
c             do 206 id=0,numid
c             do 206 ip=0,numip
c             do 206 inn=0,numin
c               if (inn+ip+id+it+ih+iaa.ne.npar) goto 206
c               ident=100000*inn+10000*ip+1000*id+100*it+10*ih+iaa
c               do 208 idc=0,idnum
c                 if (idchannel(idc).eq.ident) then
c                   xsc=xschannel(idc)
c                   if (inn.gt.0) xsexcpart(1,inn)=xsexcpart(1,inn)+xsc
c                   if (ip.gt.0) xsexcpart(2,ip)=xsexcpart(2,ip)+xsc
c                   if (id.gt.0) xsexcpart(3,id)=xsexcpart(3,id)+xsc
c                   if (it.gt.0) xsexcpart(4,it)=xsexcpart(4,it)+xsc
c                   if (ih.gt.0) xsexcpart(5,ih)=xsexcpart(5,ih)+xsc
c                   if (iaa.gt.0) xsexcpart(6,iaa)=xsexcpart(6,iaa)+xsc
c                 endif
c 208           continue
c 206         continue
c 202       continue
c           yA=yieldApre(ia)
c           yZA=yieldZApre(iz,in)
c           do 210 type=0,6
c             nubar(type)=nubar(type)+yZA*multiplicity(type)
c             if (yA.gt.0.) then
c               nupre(type,ia)=nupre(type,ia)+yZA/yA*multiplicity(type)
c               do 220 npar=1,numin
c                 Pdisnu(type,npar)=Pdisnu(type,npar)+
c    +              xsexcpart(type,npar)/xsinitpop
c 220           continue
c               do 230 nen=0,numen2
c                 pfns(type,nen)=pfns(type,nen)+yZA*xssumout(type,nen)
c 230           continue
c             endif
c 210       continue
c         endif
   30   continue
   20 continue
      write(*,'(/" ########## End of loop over residual products"/)')
c
c Average energy and relation to Maxwellian
c
c     do 290 type=0,6
c       sumpfns=0.
c       Esumpfns=0.
c       do 292 nen=1,numen2
c         dE=espec(type,nen)-espec(type,nen-1)
c         sumpfns=sumpfns+pfns(type,nen)*dE
c         Esumpfns=Esumpfns+espec(type,nen)*pfns(type,nen)*dE
c         maxpfns(type,nen)=0.
c 292   continue
c       if (sumpfns.gt.0.) then
c         Eavpfns(type)=Esumpfns/sumpfns
c       else
c         Eavpfns(type)=0.
c       endif
c       summax=0.
c       do 294 nen=1,numen2
c         dE=espec(type,nen)-espec(type,nen-1)
c         E=espec(type,nen)
c         Eav=Eavpfns(type)
c         maxwell=sqrt(E)*exp(-E/Eav)
c         summax=summax+maxwell*dE
c 294   continue
c       do 296 nen=1,numen2
c         E=espec(type,nen)
c         Eav=Eavpfns(type)
c         maxwell=sqrt(E)*exp(-E/Eav)
c         if (maxwell.gt.0..and.sumpfns.gt.0.)
c    +      maxpfns(type,nen)=pfns(type,nen)/maxwell*summax/sumpfns
c 296   continue
c 290 continue
c     sumpost=0.
c     do 310 ia=1,Atarget0
c       sumpost=sumpost+xsfpApost(ia)
c 310 continue
c     sumpost=0.5*sumpost
c     if (sumpost.gt.0.) then
c       xsfptotpost=0.
c       yieldtotpost=0.
c       do 320 iz=1,Ztarget0
c         do 330 ia=iz+1,Atarget0
c           if (xsfpZApost(iz,ia).eq.0.) goto 330
c           in=ia-iz
c           if (in.gt.numneu) goto 330
c           yieldZApost(iz,in)=xsfpZApost(iz,in)/sumpost
c           yieldApost(ia)=yieldApost(ia)+yieldZApost(iz,in)
c           yieldnpost(in)=yieldnpost(in)+yieldZApost(iz,in)
c           xsfptotpost=xsfptotpost+xsfpZApost(iz,in)
c           yieldtotpost=yieldtotpost+yieldZApost(iz,in)
c 330     continue
c 320   continue
c     endif
c
c Reset variables to those of original target.
c
c flagffruns: flag to denote that run is for residual product
c
      flagrpruns=.false.
      call talysinput
      flagmain=.false.
      call talysinitial
      flagmain=.true.
      if (.not.flagompall) call basicxs(0,0)
      if (parinclude(0)) call gamma(0,0)
      if (enincmax.ge.epreeq.or.flagracap) then
        call preeqinit
        call excitoninit
      endif
      if (flagracap) call racapinit
      if (flagcomp) call compoundinit
      if (flagastro) call astroinit
c
c Output
c
c     call massdisout
c     call nubarout
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

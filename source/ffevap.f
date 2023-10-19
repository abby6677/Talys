      subroutine ffevap
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Jean-Francois Lemaitre
c | Date  : December 11, 2021
c | Task  : Evaporation of fission fragments
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*8  ffstring
      character*14 fffile
      integer iz,ia,Zix,Nix,izp,inp,iap,nex,nex0,type,nen,npar,in,ih,it,
     +        id,ip,ident,idc,iaa,inn,nen2,izH,iaH,inH,i,ZCN,ACN
      real    xsexcpart(0:numpar,0:numin),sumpost,sumpfns,Esumpfns,E,
     +        Eav,fiseps,maxwell,summax,dE,xsc,yZA,yA,spec,sum,fac,
     +        Ekintot,Ekinff,Eb,Ee,sqrtEkinff,sqrtE,Eex,dEex,fac1,fac2,
     +        gauss,Pmultiff(numelem,numneu,0:numpar,0:numnu),Emaxff,
     +        gau(0:numpop),term,speccm,xsb,xse,sumpfnscm,Ecm1,Ecm2,
     +        xs1,xs2
c
c ********************** Loop over fission fragments *******************
c
c Do a full TALYS calculation for each fission fragment and incident
c energy
c
      flagffruns=.true.
      do type=0,6
        Epfnsaverage(type)=0.
        do nen=0,numpfns
          pfns(type,nen)=0.
          pfnscm(type,nen)=0.
          maxpfns(type,nen)=0.
        enddo
      enddo
      fiseps=Rfiseps*xsfistot
      Pmultiff=0.
      write(*,'(/" ########## Start of loop over fission fragments"/)')
c
c Use Viola systematics for first order guess of kinetic energy of FF
c
      ZCN=Ztarget0+parZ(k0)
      ACN=Atarget0+parA(k0)
      Ekintot=0.1189*Ztarget0*Ztarget0/((Atarget0+1)**onethird)+7.3
      do 20 ia=1,ACN
        do 30 iz=1,ZCN
          in=ia-iz
          if (in.lt.1.or.in.gt.numneu) goto 30
          if (xsZApre(iz,in).le.fiseps) goto 30
          Ekinff=real(ACN-ia)/real(ia)*Ekintot/ACN
c
c Okumura model
c
          ffstring='ff000000'
          write(ffstring(3:5),'(i3.3)') iz
          write(ffstring(6:8),'(i3.3)') ia
          if (fymodel.ge.4) then
            fffile=ffstring//'.ex'
            if (fymodel.eq.4) then
              open (unit=1,file=fffile,status='replace')
              write(1,'(i6,2i4,f12.5,2es12.5,f12.5,a)') numpop,0,1,
     +          Einc+S(0,0,1)+Q(1),yieldZApre(iz,in),xsZApre(iz,in),
     +          TKE(iz,in)," (numpop numJ numparity Ex Yff xsff Ekin)"
              Eex=Excff(iz,in)
              dEex=dExcff(iz,in)
              Emaxff=min(Eex+4.*dEex,80.)
              fac1=1./(dEex*sqrttwopi)
              fac2=1./(2.*dEex**2)
              dE=Emaxff/numpop
              sum=0.
              do nex=0,numpop
                E=dE*nex
                gau(nex)=fac1*exp(-((E-Eex)**2)*fac2)
                sum=sum+dE*gau(nex)
              enddo
              do nex=0,numpop
                E=dE*nex
                gauss=gau(nex)/sum
                write(1,'(f10.5,es12.5)') E,dE*gauss*xsZApre(iz,in)
              enddo
              close(1)
            endif
            Ekinff=real(ACN-ia)/real(ia)*TKE(iz,in)/ACN
          endif
          sqrtEkinff=sqrt(Ekinff)
          Aff=ia
          Zff=iz
          call evaptalys
c
c Add fission product cross sections
c
c iap: mass number
c izp: charge number
c inn: neutron number
c Eav: average energy
c Ecm: C.M. energy
c Ee: energy
c Ekinff: kinetic energy of F.F.
c
          do 110 Zix=0,maxZ
            do 120 Nix=0,maxN
              izp=iz-Zix
              if (izp.lt.1) goto 120
              inp=in-Nix
              if (inp.lt.1) goto 120
              iap=izp+inp
              xsZApost(izp,inp)=xsZApost(izp,inp)+xspopnuc(Zix,Nix)
              xsApost(iap)=xsApost(iap)+xspopnuc(Zix,Nix)
              do 130 nex=0,Nlast(Zix,Nix,0)
                if (nex.eq.0.or.tau(Zix,Nix,nex).ne.0.) then
                  nex0=min(nex,1)
                  xsfpex(izp,inp,nex0)=xsfpex(izp,inp,nex0)+
     +              xspopex(Zix,Nix,nex)
                endif
  130         continue
  120       continue
  110     continue
c
c Add prompt fission particle and gamma production and spectra
c
c npar: counter
c iaa: counter
c yA : pre-neutron emission mass yield
c xsexcpart: partial cross section
c yZA : pre-neutron emission isotopic yield
c
          if (xsinitpop.gt.0.) then
            do 201 npar=0,numin
              do 202 type=0,6
                xsexcpart(type,npar)=0.
  202         continue
              do 203 iaa=0,numia
              do 204 ih=0,numih
              do 205 it=0,numit
              do 206 id=0,numid
              do 207 ip=0,numip
              do 208 inn=0,numin
                if (inn+ip+id+it+ih+iaa.ne.npar) goto 208
                ident=100000*inn+10000*ip+1000*id+100*it+10*ih+iaa
                do 209 idc=0,idnum
                  if (idchannel(idc).eq.ident) then
                    xsc=xschannel(idc)
                    xsexcpart(1,inn)=xsexcpart(1,inn)+xsc
                    xsexcpart(2,ip)=xsexcpart(2,ip)+xsc
                    xsexcpart(3,id)=xsexcpart(3,id)+xsc
                    xsexcpart(4,it)=xsexcpart(4,it)+xsc
                    xsexcpart(5,ih)=xsexcpart(5,ih)+xsc
                    xsexcpart(6,iaa)=xsexcpart(6,iaa)+xsc
                  endif
  209           continue
  208         continue
  207         continue
  206         continue
  205         continue
  204         continue
  203         continue
  201       continue
            yA=yieldApre(ia)
            yZA=yieldZApre(iz,in)
            do 210 type=0,6
              if (parskip(type)) goto 210
              nubar(type)=nubar(type)+yZA*multiplicity(type)
              if (yA.gt.0.) then
                nuA(type,ia)=nuA(type,ia)+yZA/yA*multiplicity(type)
                nuZA(type,iz,in)=nuZA(type,iz,in)+yZA*multiplicity(type)
                do 220 npar=0,numin
                  Pmultiff(iz,in,type,npar)=
     +              xsexcpart(type,npar)/xsinitpop
  220           continue
                term=yZa*ACN/(ACN-ia)
                Epfnsaverage(type)=Epfnsaverage(type)+
     +            0.5*term*Eaverage(type)
                EaverageZA(type,iz,in)=Eaverage(type)
                EaverageA(type,ia)=EaverageA(type,ia)+
     +            yZA/yA*Eaverage(type)
                fffile=ffstring//'.'//parsym(type)//'spec'
                open (unit=1,file=fffile,status='unknown')
                write(1,'("# Prompt fission ",a8," spectrum of ",
     +            i3,a2)') parname(type),ia,nuc(iz)
                write(1,'("# Number of energies:",i6)') NEpfns+1
                write(1,'("#    E-out spectrum_CM spectrum_lab")')
                do 230 nen=1,NEpfns
                  E=Epfns(nen)
c
c C.M.
c
                  speccm=0.
                  do 235 nen2=0,eend(type)-1
                    Eb=espec(type,nen2)
                    Ee=espec(type,nen2+1)
                    xsb=xssumout(type,nen2)
                    xse=xssumout(type,nen2+1)
                    if (Eb.eq.0..and.Ee.eq.0.) goto 235
                    if (E.ge.Eb.and.E.le.Ee) then
                      call pol1(Eb,Ee,xsb,xse,E,speccm)
                      goto 237
                    endif
  235             continue
c
c Lab
c
  237             sum=0.
                  sqrtE=sqrt(E)
                  Eb=(sqrtE-sqrtEkinff)**2
                  Ee=(sqrtE+sqrtEkinff)**2
                  fac=1./(4.*sqrtEkinff)
c
c Trapezoidal rule to integrate
c
                  if (abs(Ee-Eb).lt.1.e-5) goto 250
                  do 240 nen2=1,eend(type)-1
                    Ecm1=espec(type,nen2)
                    Ecm2=espec(type,nen2+1)
                    if (Ecm1.le.Eb.and.Ecm2.le.Eb) goto 240
                    if (Ecm1.ge.Ee.and.Ecm2.ge.Ee) goto 240
                    xs1=0.
                    xs2=0.
                    if (Ecm1.gt.0.) xs1=xssumout(type,nen2)/sqrt(Ecm1)
                    if (Ecm2.gt.0.) xs2=xssumout(type,nen2+1)/sqrt(Ecm2)
                    term=0.5*(xs1+xs2)
                    if (Ecm1.lt.Eb.and.Ecm2.ge.Eb) dE=Ecm2-Eb
                    if (Ecm1.gt.Eb.and.Ecm2.lt.Ee) dE=Ecm2-Ecm1
                    if (Ecm1.le.Ee.and.Ecm2.gt.Ee) dE=Ee-Ecm1
                    sum=sum+term*dE
  240             continue
  250             spec=sum*fac
                  if (type.eq.0) spec=speccm
                  pfns(type,nen)=pfns(type,nen)+spec
                  pfnscm(type,nen)=pfnscm(type,nen)+speccm
                  write(1,'(f10.5,2es12.5)') E,speccm,spec
  230           continue
                close(1)
              endif
  210       continue
          endif
   30   continue
   20 continue
      write(*,'(/" ########## End of loop over fission fragments"/)')
c
c Calculate P(nu)
c
      do 281 ia=1,ACN/2
        do 282 iz=1,ZCN-1
          in=ia-iz
          if (in.lt.1.or.in.gt.numneu) goto 282
          if (xsZApre(iz,in).lt.fiseps.and..not.fpexist(iz,in))
     +      goto 282
          do 283 type=0,6
            if (parskip(type)) goto 283
            sum=0.
            do 284 npar=0,numin
              sum=sum+Pmultiff(iz,in,type,npar)
  284       continue
            Pmultiff(iz,in,type,npar)=Pmultiff(iz,in,type,npar)/sum
            izH=ZCN-iz
            iaH=ACN-ia
            inH=iaH-izH
            sum=0.
            do 285 npar=0,numin
              sum=sum+Pmultiff(izH,inH,type,npar)
  285       continue
            Pmultiff(izH,inH,type,npar)=Pmultiff(izH,inH,type,npar)/sum
            do 286 npar=0,numin
              do 287 i=0,npar
                Pdisnu(type,npar)=Pdisnu(type,npar)+yieldZApre(iz,in)*
     +            Pmultiff(iz,in,type,i)*Pmultiff(izH,inH,type,npar-i)
  287         continue
  286       continue
  283     continue
  282   continue
  281 continue
      do 288 type=0,6
        if (parskip(type)) goto 288
        sum=0.
        do 289 npar=0,numin
          sum=sum+Pdisnu(type,npar)
  289   continue
        Pdisnuav(type)=0.
        if (sum.gt.0.) then
          do 290 npar=0,numin
            Pdisnu(type,npar)=Pdisnu(type,npar)/sum
            Pdisnuav(type)=Pdisnuav(type)+npar*Pdisnu(type,npar)
  290     continue
        endif
  288 continue
c
c Average energy and relation to Maxwellian
c
c Ekintot: total kinetic energy
c Esumpfns: integrated PFNS
c maxwell: Maxwell distribution
c sumpfns: integrated PFNS
c sumpost: sum over post-neutron FP's
c sqrtE: square root of energy
c sqrtEkinff: square root of kinetic energy of FF's
c summax : integral over Maxwellian
c spec  : spectrum
c
      do 291 type=0,6
        if (parskip(type)) goto 291
        sumpfns=0.
        sumpfnscm=0.
        Esumpfns=0.
        Eavpfns(type)=0.
        do 292 nen=1,NEpfns
          sumpfns=sumpfns+pfns(type,nen)*dEpfns(nen)
          sumpfnscm=sumpfnscm+pfnscm(type,nen)*dEpfns(nen)
  292   continue
        if (sumpfns.gt.0.) then
          do 293 nen=1,NEpfns
            pfns(type,nen)=pfns(type,nen)/sumpfns
            if (sumpfnscm.gt.0.)
     +        pfnscm(type,nen)=pfnscm(type,nen)/sumpfnscm
            Esumpfns=Esumpfns+Epfns(nen)*pfns(type,nen)*dEpfns(nen)
            maxpfns(type,nen)=0.
  293     continue
          Eavpfns(type)=Esumpfns
          Eav=Eavpfns(type)
          if (Eav.gt.0.) then
            summax=0.
            do 294 nen=1,NEpfns
              E=Epfns(nen)
              maxwell=sqrt(E)*exp(-E/Eav)
              summax=summax+maxwell*dEpfns(nen)
  294       continue
            do 296 nen=1,NEpfns
              E=Epfns(nen)
              maxwell=sqrt(E)*exp(-E/Eav)
              if (maxwell.gt.0.)
     +          maxpfns(type,nen)=pfns(type,nen)/maxwell*summax
  296       continue
          endif
        endif
  291 continue
      sumpost=0.
      do 310 ia=1,ACN
        sumpost=sumpost+xsApost(ia)
  310 continue
      sumpost=0.5*sumpost
      if (sumpost.gt.0.) then
        xstotpost=0.
        yieldtotpost=0.
        do 320 iz=1,ZCN
          do 330 ia=iz+1,ACN
            in=ia-iz
            if (in.gt.numneu) goto 330
            if (xsZApost(iz,in).eq.0.) goto 330
            yieldZApost(iz,in)=xsZApost(iz,in)/sumpost
            yieldApost(ia)=yieldApost(ia)+yieldZApost(iz,in)
            xstotpost=xstotpost+xsZApost(iz,in)
            yieldtotpost=yieldtotpost+yieldZApost(iz,in)
            if (xsfpex(iz,in,1).gt.0.) then
              do 340 nex=0,1
                yieldfpex(iz,in,nex)=xsfpex(iz,in,nex)/sumpost
                if (yieldZApost(iz,in).gt.0.) fpratio(iz,in,nex)=
     +            yieldfpex(iz,in,nex)/yieldZApost(iz,in)
  340         continue
            endif
  330     continue
  320   continue
      endif
c
c Reset variables to those of original target.
c
c flagffruns: flag to denote that run is for fission fragment
c
      flagffruns=.false.
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
      call massdisout
      call nubarout
      call nudisout
      if (flagspec) call pfnsout
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely

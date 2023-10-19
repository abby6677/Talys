      subroutine tfission(Zcomp,Ncomp,nex,J2,parity)
c
c +---------------------------------------------------------------------
c | Author: Stephane Hilaire and Pascal Romain
c | Date  : October 11, 2019
c | Task  : Fission transmission coefficients
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          Zcomp,Ncomp,nex,J2,parity,J,iloop,ihill,ic2,jc2,
     +                 pc2
      real             Eex,term1,term2,damper,expo,ec2,wo2,diffnrj,
     +                 wo2damp,boost,boostmax,Ecut,Ecut1,Ecut2,term11,
     +                 term21,term12,term22,damper1,damper2,wo2damp1,
     +                 wo2damp2,Twkbtransint
      double precision tf,tfb1,rnfb1,tfb2,rnfb2,tfb3,rnfb3,tfii,addnrj,
     +                 tf12,tsum123,tfiii,density,tdir,trans,rnorm,
     +                 sumt2,sumt3,ta12,ta13,ta23,ta32,tdir12,tdir13,
     +                 tdir21,tdir23,tgam2,tgam3,ti2,ti3,trans2,trans3
      external Twkbtransint
c
c ********** Calculation of fission transmission coefficients **********
c
c Zcomp          : charge number index for compound nucleus
c Ncomp          : neutron number index for compound nucleus
c J2             : 2 * J
c parity         : parity of target
c dExinc         : excitation energy bin for mother nucleus
c deltaEx        : excitation energy bin for population arrays
c tf             : help variable
c Eex            : excitation energy
c Exinc          : excitation energy of entrance bin
c tfisdown,tfisup: fission transmission coefficients
c tfis           : fission transmission coefficients
c numhill        : maximum number of Hill-Wheeler points
c tfisA          : transmission coefficient for Hill-Wheeler magnitude
c rhofisA        : integrated level density corresponding to tfisA
c
c The fission transmission coefficients decrease very rapidly with
c excitation energy. Therefore, we calculate them at the endpoints and
c at the middle of each excitation energy bin. With this information, we
c can do logarithmic integration in subroutine compound.
c
c J and parity are in loops outside this subroutine
c
      J=J2/2
      dExinc=deltaEx(Zcomp,Ncomp,nex)
      do 10 iloop=1,3
        tf=0.
        if (iloop.eq.1) then
          Eex=max(Exinc-0.5*dExinc,0.)
          tfisdown(J,parity)=0.
        endif
        if (iloop.eq.2) then
          Eex=Exinc
          tfis(J,parity)=0.
          gamfis(J,parity)=0.
          taufis(J,parity)=0.
          denfis(J,parity)=0.
          do 20 ihill=0,numhill
            tfisA(J,parity,ihill)=0.
            rhofisA(J,parity,ihill)=1.
   20     continue
        endif
        if (iloop.eq.3) then
          Eex=min(Exinc+0.5*dExinc,Exmax(Zcomp,Ncomp))
          tfisup(J,parity)=0.
        endif
c
c 1. One barrier
c
c nfisbar   : number of fission barrier parameters
c t1barrier : subroutine for fission transmission coefficient for one
c             barrier
c tfb1      : help variable
c tfb2      : help variable
c tfb3      : help variable
c tfbii     : help variable
c tfbiii    : help variable
c tsum123   : help variable
c rnfb1: help variable
c rnfb2: help variable
c rnfb3: help variable
c
        if (nfisbar(Zcomp,Ncomp).eq.1) then
          call t1barrier(Zcomp,Ncomp,J2,parity,1,tfb1,rnfb1,Eex,iloop)
          tf=tfb1
        endif
c
c 2. Two barriers
c
c transeps: absolute limit for transmission coefficient
c flagfispartdamp: flag for fission partial damping
c
        if (nfisbar(Zcomp,Ncomp).eq.2) then
          if (flagfispartdamp) call tdirbarrier(Zcomp,Ncomp,J2,parity,
     +      1,2,tdir,rnfb1,Eex,iloop)
          call t1barrier(Zcomp,Ncomp,J2,parity,1,tfb1,rnfb1,Eex,iloop)
          if (tfb1.lt.transeps) goto 30
          call t1barrier(Zcomp,Ncomp,J2,parity,2,tfb2,rnfb2,Eex,iloop)
          if (tfb2.lt.transeps) goto 30
          if (flagfispartdamp) then
            trans=Twkbtransint(Eex,1,Zcomp,Ncomp)
            tf=tfb1*tfb2/(tfb1+tfb2)*trans + tdir*(1.-trans)
          else
            tf=tfb1*tfb2/(tfb1+tfb2)
          endif
c
c ****************** Special treatment for class2 states ***************
c
c Ecut       : help variable
c Emaxclass2 : maximum energy for class2 states
c widthc2    : width of class2 states
c flagclass2 : flag for class2 states in fission
c term1,term2: help variables
c damper     : damping factor
c efisc2rot  : energy of rotational class2 states
c nfisc2rot  : number of rotational class2 states per set
c nc2        : help variable
c ic2    : help variable
c ec2: help variable
c wo2: help variable
c jc2: help variable
c pc2: help variable
c addnrj    : added energy
c jfisc2rot  : spin of rotational class2 states
c pfisc2rot  : parity of rotational class2 states
c
          Ecut=Emaxclass2(Zcomp,Ncomp,1)+0.5*widthc2(Zcomp,Ncomp,1)
          if (flagclass2.and.(Eex.le.Ecut)) then
            term1=-Eex+0.5*(efisc2rot(Zcomp,Ncomp,1,1) +
     +        efisc2rot(Zcomp,Ncomp,1,nfisc2rot(Zcomp,Ncomp,1)))
            term2=efisc2rot(Zcomp,Ncomp,1,nfisc2rot(Zcomp,Ncomp,1))-
     +        efisc2rot(Zcomp,Ncomp,1,1)
            damper=1.
            if (term2.gt.0.) then
              expo=24.*term1/term2
              if (abs(expo).le.80.) damper=1./(1.+exp(expo))
            endif
            wo2=0.5*widthc2(Zcomp,Ncomp,1)
            wo2damp=wo2*damper
            tfii=0.
            do 40 ic2=nfisc2rot(Zcomp,Ncomp,1),1,-1
              ec2=efisc2rot(Zcomp,Ncomp,1,ic2)
              diffnrj=abs(Eex-ec2)
              addnrj=diffnrj/wo2damp
              boost=1./(1.+addnrj**2)
              if (boost.ge.0.25) then
                jc2=int(2.*jfisc2rot(Zcomp,Ncomp,1,ic2))
                pc2=pfisc2rot(Zcomp,Ncomp,1,ic2)
                if ((jc2.eq.J2).and.(pc2.eq.parity)) then
                  boostmax=4./(tfb1+tfb2)
                  boost=boostmax*boost
                  tfii=tfii+tf*boost
                endif
              endif
   40       continue
            if (tfii.gt.0.) tf=tfii
          endif
        endif
c
c 3. Three barriers
c
        if (nfisbar(Zcomp,Ncomp).eq.3) then
          call t1barrier(Zcomp,Ncomp,J2,parity,1,tfb1,rnfb1,Eex,iloop)
          if (tfb1.lt.transeps) goto 30
          call t1barrier(Zcomp,Ncomp,J2,parity,2,tfb2,rnfb2,Eex,iloop)
          if (tfb2.lt.transeps) goto 30
          call t1barrier(Zcomp,Ncomp,J2,parity,3,tfb3,rnfb3,Eex,iloop)
          if (tfb3.lt.transeps) goto 30
          if (flagfispartdamp) then
            call tdirbarrier(Zcomp,Ncomp,J2,parity,1,2,
     +        tdir12,rnfb1,Eex,iloop)
            call tdirbarrier(Zcomp,Ncomp,J2,parity,2,3,
     +        tdir23,rnfb1,Eex,iloop)
            call tdirbarrier(Zcomp,Ncomp,J2,parity,1,3,
     +        tdir13,rnfb1,Eex,iloop)
            trans2=Twkbtransint(Eex,1,Zcomp,Ncomp)
            trans3=Twkbtransint(Eex,2,Zcomp,Ncomp)
            tdir21=(1-trans2)*tdir12
            tdir12=(1-trans2)*tdir12
            tdir23=(1-trans3)*tdir23
            tdir13=(1-trans2)*(1-trans3)*tdir13
            Ta12=tfb1*trans2
            Ta23=tfb2*trans3
            Ta13=trans3*tdir12
            Ta32=trans2*tfb2
            tgam2=0.
            tgam3=0.
            sumT2= tfb1+tdir23+Ta23+Tgam2
            sumT3= tdir21+tfb3+Ta32+Tgam3
            Ti2=Ta12*(tdir23/sumT2+Ta23*tfb3/(sumT2*sumT3))
            Ti3=Ta13*(tfb3/sumT3+Ta32*tdir23/(sumT2*sumT3))
            if (abs((Ta23*Ta32)/(sumT2*sumT3)-1.).le.1.e-8 )  then
              Rnorm=  sumT2*sumT3/( (tfb1+tdir23+Tgam2)*sumT3
     +          + Ta23*(tdir21+tfb3+Tgam3) )
            else
              Rnorm=1./(1.-(Ta23*Ta32)/(sumT2*sumT3))
            endif
            tf=tdir13+Rnorm*(Ti2+Ti3)
          else
            tf12=tfb1*tfb2/(tfb1+tfb2)
            tsum123=tf12+tfb3
            tf=tf12*tfb3/tsum123
          endif
c
c *********** Special treatment for class2 and class3 states ***********
c
c Ecut1: cutoff energy
c Ecut2: cutoff energy
c term11: help variable
c term21: help variable
c term22: help variable
c damper1 : energy damping function
c damper2 : energy damping function
c ec2   : rotational energy
c diffnrj: energy difference
c boost: energy boost
c boostmax: maximum energy boost
c wo2damp : energy damping function
c wo2damp1: energy damping function
c wo2damp2: energy damping function
c tfii    : help variable
c tfiii   : help variable
c tf12    : help variable
c
          Ecut1=Emaxclass2(Zcomp,Ncomp,1)+0.5*widthc2(Zcomp,Ncomp,1)
          Ecut2=Emaxclass2(Zcomp,Ncomp,2)+0.5*widthc2(Zcomp,Ncomp,2)
          Ecut=max(Ecut1,Ecut2)
          if (flagclass2.and.(Eex.le.Ecut)) then
            term11=-Eex+0.5*(efisc2rot(Zcomp,Ncomp,1,1) +
     +        efisc2rot(Zcomp,Ncomp,1,nfisc2rot(Zcomp,Ncomp,1)))
            term21=efisc2rot(Zcomp,Ncomp,1,nfisc2rot(Zcomp,Ncomp,1))-
     +        efisc2rot(Zcomp,Ncomp,1,1)
            term12=-Eex+0.5*(efisc2rot(Zcomp,Ncomp,2,1) +
     +        efisc2rot(Zcomp,Ncomp,2,nfisc2rot(Zcomp,Ncomp,2)))
            term22=efisc2rot(Zcomp,Ncomp,2,nfisc2rot(Zcomp,Ncomp,2))-
     +        efisc2rot(Zcomp,Ncomp,2,1)
            damper1=1.
            if (term21.gt.0.) then
              expo=24.*term11/term21
              if (abs(expo).le.80.) damper1=1./(1.+exp(expo))
            endif
            damper2=1.
            if (term22.gt.0.) then
              expo=24.*term12/term22
              if (abs(expo).le.80.) damper2=1./(1.+exp(expo))
            endif
            wo2=0.5*widthc2(Zcomp,Ncomp,1)
            wo2damp1=wo2*damper1
            tfii=0.
            do 50 ic2=nfisc2rot(Zcomp,Ncomp,1),1,-1
              ec2=efisc2rot(Zcomp,Ncomp,1,ic2)
              diffnrj=abs(Eex-ec2)
              addnrj=diffnrj/wo2damp1
              boost=1./(1.+addnrj**2)
              if (boost.ge.0.25) then
                jc2=int(2.*jfisc2rot(Zcomp,Ncomp,1,ic2))
                pc2=pfisc2rot(Zcomp,Ncomp,1,ic2)
                if ((jc2.eq.J2).and.(pc2.eq.parity)) then
                  boostmax=4./(tfb1+tfb2)
                  boost=boostmax*boost
                  tfii=tfii+tfb1*tfb2/(tfb1+tfb2)*boost
                endif
              endif
   50       continue
            if (tfii.gt.0.) then
               tf12=tfii
               tsum123=tf12+tfb3
               tf=tf12*tfb3/tsum123
            endif
            wo2=0.5*widthc2(Zcomp,Ncomp,2)
            wo2damp2=wo2*damper2
            tfiii=0.
            do 60 ic2=nfisc2rot(Zcomp,Ncomp,2),1,-1
              ec2=efisc2rot(Zcomp,Ncomp,2,ic2)
              diffnrj=abs(Eex-ec2)
              addnrj=diffnrj/wo2damp2
              boost=1./(1.+addnrj**2)
              if (boost.ge.0.25) then
                jc2=int(2.*jfisc2rot(Zcomp,Ncomp,2,ic2))
                pc2=pfisc2rot(Zcomp,Ncomp,2,ic2)
                if ((jc2.eq.J2).and.(pc2.eq.parity)) then
                  boostmax=4./tsum123
                  boost=boostmax*boost
                  tfiii=tfiii+tf*boost
                endif
              endif
   60       continue
            if (tfiii.gt.0.) tf=tfiii
          endif
        endif
c
c Optional adjustment factors
c
c Fnorm    : multiplication factor
c
   30   tf=tf*Fnorm(-1)
        if (iloop.eq.1) tfisdown(J,parity)=tf
        if (iloop.eq.2) tfis(J,parity)=tf
        if (iloop.eq.3) tfisup(J,parity)=tf
   10 continue
c
c Partial fission widths and lifetimes
c
c gamfis  : fission width
c taufis  : fission lifetime
c denfis  : fission level density
c
      denfis(J,parity)=density(Zcomp,Ncomp,Exinc,real(J),parity,0,
     +  ldmodel(Zcomp,Ncomp))
      if (denfis(J,parity).gt.0.) 
     +  gamfis(J,parity)=tfis(J,parity)/(twopi*denfis(J,parity))
      if (gamfis(J,parity).gt.0.) 
     +  taufis(J,parity)=hbar/gamfis(J,parity)
      return
      end
Copyright (C)  2019 A.J. Koning, S. Hilaire and S. Goriely

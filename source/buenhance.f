      subroutine buenhance(Z00,N00,BFsum,Ebu)
c +---------------------------------------------------------------------
c
c | Author   : Marilena Avrigeanu
c | Date     : April 2020
c | Task     : Inelastic breakup enhancement brought by breakup neutrons
c              and protons to any (Z,N) residual nucleus from d+Atarge
c              interaction, Eq. 2 from Phys.Rev.C 94,014606(2016) and
c              References therein
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,Z,N,nen,Z00,N00,Zix,Nix
      real BFenhance(2),enhance03(2),conv03(2),Ebreak(2),BCout(2)
      real sumtest0(2,numenout),sumtest(2),sumGauss(2)
      real Emaxn,E0n03,term03,Gauss,width,w03,En,BCin,Bdeut,BFsum,Ebu
c
c ************************** Avrigeanu model ***************************
c
c DEUTERON Break-up model by M. and V. Avrigeanu, PRC95,024607(2017)
c Inelastic breakup enhancement, M.Avrigeanu et al., PRC94,014606(2016)
c
c numenout    : maximal number of breakup nucleons outgoing energies
c               included in talys.cmb as well as in reacinitial
c Ebu         : incident energy in MeV
c Atarget     : mass number of target nucleus
c BCin,BCout  : effective Coulomb barrier
c Ztarget     : charge number of target nucleus
c Bdeut       : deuteon binding energy
c Emaxn       : maximum breakup nucleons energy in Laboratory System
c Ebreak      : Centroid energy of the breakup nucleon energy
c               distributions, Kalbach 2003, in Laboratory System
c width       : Full width at half maximum of the breakup nucleon
c               energy distribution, Kalbach 2003
c Z00         : charge number for residual nucleus
c N00         : neutron number for residual nucleus
c Z           : charge number for residual nucleus
c N           : neutron number for residual nucleus
c ENHratio    : breakup nucleons enhancing reaction cross sections
c               ratios, Eq. (2) from PRC94,014606(2016)          
c               n + Atarget  sig(n,Z,A,Eout)/sig_Total(n,Enout);
c               p + Atarget  sig(p,Z,A,Eout)/sig_Reaction(p,Epout)
c BFenhance   : inelastic breakup enhancement brought by breakup
c               nucleons to any (Z,N) residual nucleus from d+Atarge
c               interaction
c BFsum       : total inelastic breakup enhancement to (Z,N) residual
c               nucleus population
c
      BFsum = 0.
c
      do type=1,2
        Ebreak(type)=0.
        conv03(type)=0.
        enhance03(type)=0.
        sumtest(type)=0.
        sumGauss(type)=0.
        BFenhance(type)=0.
      enddo
c
      do 103 type=1,2
      do 103 nen=1,numenout
        sumtest0(type,nen)=0.
103   continue
c
      Z=Z00
      N=N00
      Bdeut = 2.225
      BCin=Ztarget/9.5
      BCout(1)=0.
      BCout(2)=BCin
      Emaxn = Ebu*(Atarget+1.)/(Atarget+2.)-Bdeut*(Atarget+1.)/Atarget
c
c      Breakup threshold
c
      if(Emaxn.le.0.d+00) then
        Emaxn=0.
        write(8,*)" deuteron energy lower than Bd, no breakup"
        go to 499
      endif
c
c----------------------     breakup nucleon type  ----------------------
c
      do 401, type=1,2
c
c
        Ebreak(type) = 0.5*(Atarget+1.)/(Atarget+2.)*Ebu+
     +      0.5*(Atarget+1.)/Atarget*(-Bdeut-BCin+2.*BCout(type))
        width = 1.15 + 0.12*Ebu - Atarget/140.
C
        if(width.le.0.d+00) then
          width=0.
          go to 499
        endif
c
        if(Ebreak(type).le.0.d+00) then
          Ebreak(type)=0.01
        endif
c
        E0n03=Ebreak(type)
        w03=width
c
        conv03(type)     = 0.
        sumGauss(type)   = 0.
        BFenhance(type)  = 0.
c
c----------------------     breakup nucleon energy   -------------------
c
c xsBF        : nucleon inelastic breakup c.s. (see breakupAVR)
c ebubin      : set = 0.1 (see npxsratios), outgoing breakup nucleon
c               energy bin for integration in Eq.(2), PRC94,014606.
c ebegin      : first energy point of energy grid
c eend        : last energy point of energy grid
c En,egrid    : outgoing energy
c term03      : help variable for  breakup enhancement calculations
c conv03      : help variable for  breakup enhancement calculations
c enhance03   : help variable for  breakup enhancement calculations
c sumtest0    : help variables for ckhecking the calculation procedure
c sumtest     : help variables for ckhecking the calculation procedure
c sumGauss    : help variables for ckhecking the calculation procedure
c
         Zix=Zinit-Z
         Nix=Ninit-N
         DO 301, nen=1,numenout
           En=ebubin*nen     
          if(En.le.Emaxn) then
             sumtest0(type,nen) = GAUSS(En,E0n03,w03)
             term03 = ENHratio(type,Zix,Nix,nen)*
     +         GAUSS(En,E0n03,w03)*ebubin
c
c        BREAKUP THRESHOLD:  En>Emax BU
c
          else
             sumtest0(type,nen) = 0.
            term03 = 0.
          endif
c
           sumGauss(type)=sumGauss(type)+sumtest0(type,nen)*ebubin
          conv03(type) = conv03(type) + term03
301     continue
c
c--------------------   END   breakup nucleon energy   -----------------
c
        enhance03(type) = xsBF(type)*conv03(type)
         if(sumGauss(type).eq.0.d+00) go to 401
        BFenhance(type) = enhance03(type)/sumGauss(type)
        BFsum= BFsum + BFenhance(TYPE)
c
401   continue
c
c----------------------     END   breakup nucleon type  ----------------
c
499      continue
c
ccc        write(8,*)'   '
ccc        write(8,*)"   Ed    type    BFenhance     Emaxn   ",
ccc     1            " E0_SL      w03   sumtest     xsBUnuc"
ccc        write(8,*)"   Ed    type    BFenhance  "
c
        do 506 type=1,2
          if(sumGauss(type).eq.0.d+00)cycle
         do 505 nen=1,numenout
           sumtest0(type,nen)=xsBUnuc(type)*
     +    sumtest0(type,nen)/sumGauss(type)
           sumtest(type)=sumtest(type)+sumtest0(type,nen)*ebubin
505       continue
c
ccc        write(8,999)Ebu,type,BFenhance(type),Emaxn,Ebreak(type),
ccc     +              width,sumtest(type),xsBUnuc(type)
        write(8,999)Ebu,type,BFenhance(type)
506     continue
c
ccc         write(8,*)' ... ine buenhance:  Z=',Z,' A=',Z+N,nuc(Z),
ccc     +                                '    BFsum=',BFsum
         write(8,*)'                 BFsum=',
     +              BFsum,'  for  Z=',Z,' A=',Z+N,nuc(Z)
c
999     Format(2x,f8.3,1x,I3,3x,f12.5,3x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,
     +       4x,f8.3,3x,f8.3,4x,f8.3,4x,f8.3,4x,f8.3,4x,f8.3)
c
c
      return
      end

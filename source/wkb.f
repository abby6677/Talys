      subroutine wkb(Z,A,Zix,Nix,nbar)
c
c +---------------------------------------------------------------------
c | Author: Roberto Capote, Arjan Koning, Stephane Goriely, G. Scamps
c | Date  : January 8, 2021
c | Task  : Initialization of WKB approximation for fission
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer nsmooth
      integer Z,A,nbar,nbarr,i,j,n1,Find_Extrem,Zix,Nix,ii
      real    rmiu,centr,heigth,width,ucentr,uheigth,uwidth,uexc0,uexc1,
     +        tdir,dE,uexc,phase(2*numbar),tff(2*numbar),Vwell,Vwell2,
     +        Vtop,Vtop2,b1,b2
      character*9 filewkb
c
c ******* Finding maxima and minima of the deformation energy curve ****
c
c Find_Extrem: function to find extrema
c nbarr: number of fission barriers
c nsmooth: smoothing parameter
c centr: center of fitted parabola
c Nsmooth: smoothing parameter
c logmax: logical
c logmin: logical
c filewkb: WKB file
c heigth: barrier height
c ucentr: energy of center
c iext : counter
c
      nsmooth=5
      rmiu = 0.054*A**(5./3.)
c initializes iiextr() and nextr
      nextr = Find_Extrem(nsmooth)
      nextr=min(nextr,2*numbar-1)
      nbarr = nextr/2 + 1
      nbar=nbarr
c     nwell = nextr/2
c     Fitting parabola
      filewkb='         '
      write(filewkb,'("wkb",2i3.3)') Z,A
      if (flagfisout) then
        open (unit=22,file=filewkb,status='replace')
        write(22,'("Z=",i3," A=",i3)') Z,A
      endif
      do j=1,2*numbar
        Vheight(j)=0.
        Vwidth(j)=0.
      enddo
      do j=1,nextr
        CALL ParabFit(iiextr(j),3,rmiu,betafis,vfis,
     &  centr,  heigth,  width, ucentr, uheigth, uwidth)
        if(width.LT.0.05d0) CYCLE ! Skipping very narrow peaks
        if (flagfisout) then
          write(22,*) ' Def: ',betafis(iiextr(j)),
     &      ' (',centr,' +/- ',ucentr,')'
          write(22,*) ' Heigth :',vfis(iiextr(j)),
     &      ' (',heigth,' +/- ',uheigth,')'
          write(22,*) ' Width :', width, ' +/- ', uwidth
          write(22,*) '*******************************************'
        endif
c--------------------------------------------------
c       Initializes parabola's parameters
c
c       The real height of the barrier is used (not the fitted
c       parabola height)
c
c height: height of barrier
c
c
        Vheight(j) = vfis(iiextr(j))
        Vwidth(j) = width
        if (mod(j,2).eq.1) then
          ii=(j+1)/2
          fbarrier(Zix,Nix,ii)=vfis(iiextr(j))
          fwidth(Zix,Nix,ii)=width
        endif
c       The real position of the barrier is used (not the fitted
c       parabola centr)
        Vpos(j) = betafis(iiextr(j))
c--------------------------------------------------
      enddo
      Vpos(nextr+1) = 100.d0
C
c     PAUSE
C==================================================================
C     Calculating the limit for the energy loop
C     (Depends on barriers' height and depth)
C
c     ustep = 0.01d0
      uexc0 = 1.d0
      uexc1 = 7.d0
      uexc0 = Vheight(2)
      if(nextr.eq.5) uexc0 = min(uexc0,Vheight(4))
      uexc1 =     max(Vheight(1),Vheight(3))
      if(nextr.eq.5) uexc1 =  max(uexc1,Vheight(5))
C==================================================================
      if (flagfisout)
     &  write(22,'(1x,A5,2x,18(A10,1x))')  ' Uexc', '   Tdir   ',
     &  'TaTb/Ta+Tb','    Ta    ','    Tb    ', '    Tc    '
C
C     ENERGY LOOP
C
      b1=bdamp(Zix,Nix,1)
      b2=bdamp(Zix,Nix,2)
      n1=3*nbinswkb/4
      uexc=0.
      Uwkb(Zix,Nix,0)=0.
      do j=1,nbar
        Twkb(Zix,Nix,0,j)=0.
        Twkbdir(Zix,Nix,0,j)=0.
        Twkbtrans(Zix,Nix,0,j)=0.
        Twkbphase(Zix,Nix,0,j)=0.
      enddo
      do i=1,nbinswkb
        if (i.le.n1) then
          dE=uexc1/n1
        else
          dE=0.2
        endif
        uexc=uexc+dE
        CALL WKBFIS(Z,A-Z,uexc, tff, phase, tdir)
        Uwkb(Zix,Nix,i)=uexc
        if (flagfispartdamp) then
          if (nbar.gt.1) then
            Vwell=0.
            Vwell2=0.
            if (iiextr(2).gt.0) Vwell=vfis(iiextr(2))
            if (iiextr(4).gt.0) Vwell2=vfis(iiextr(4))
            Vtop=min(Vheight(1),Vheight(3))
            Vtop2=min(Vheight(3),Vheight(5))
            Twkbdir(Zix,Nix,i,1)=tdir
            if (uexc.lt.Vwell) then
              Twkbtrans(Zix,Nix,i,1)=0.
            else if (uexc.gt.Vtop) then
              Twkbtrans(Zix,Nix,i,1)=1.
            else
              Twkbtrans(Zix,Nix,i,1)= (uexc**2-Vwell**2)/
     +          ( (Vtop**2-Vwell**2)*exp(-(uexc-Vtop)/b1) )
            endif
            if (uexc.lt.Vwell2) then
              Twkbtrans(Zix,Nix,i,2)=0.
            else if (uexc.gt.Vtop2) then
              Twkbtrans(Zix,Nix,i,2)=1.
            else
              Twkbtrans(Zix,Nix,i,2)= (uexc**2-Vwell2**2)/
     +          ( (Vtop2**2-Vwell2**2)*exp(-(uexc-Vtop2)/b2) )
            endif
c           write(16,*)uexc,Twkbtrans(Zix,Nix,i,1),Twkbtrans(Zix,Nix,i,2)
          else
            Twkbdir(Zix,Nix,i,1)=0.
            Twkbtrans(Zix,Nix,i,1)=1.
            Twkbtrans(Zix,Nix,i,2)=1.
          endif
        endif
        do j=1,nbar
          Twkb(Zix,Nix,i,j)=tff(2*j-1)
          if (flagfispartdamp) Twkbphase(Zix,Nix,i,j)=phase(2*j)
          if (flagfisout)
     +      write(22,'("i=",i2," E=",f8.3," j=",i2," T=",es12.3)')
     +      i,uexc,j,Twkb(Zix,Nix,i,j)
        enddo
      enddo
      if (flagfisout) close (22)
      END
      subroutine wkbfis(iz,in,uexcit, tff, phase, tdir)
C
C     Roberto Capote, IAEA/NDS
C     e-mail: r.capotenoy@iaea.org, 24 August 2006, v1.00
C             rcapotenoy@yahoo.com
C
C     Mihaela Sin, University of Bucharest
C     e-mail: mihaela.sin@gmail.com
C
C     For theory behind see
C     M. Sin, R. Capote, A. Ventura et al, Phys.Rev.C74, 014608 (2006)
C
C     WKBFIS calculates WKB momentum integrals by Gauss-Legendre method
C     for any given one-dimensional barrier and a given excitation
C     energy.
C
C     * uexcit is the excitation energy
C     * tfis is the fission transmission coefficient (see eq.(4) PRC)
C     * tff(extrema), phase(extrema) are the transmission coefficients
C       and phase integrals for all extrema
C       tff(1)= Ta, tff(3)= Tb, tff(5)= Tc (for triple humped barrier)
C
c     IMPLICIT none
      include "talys.cmb"
      INTEGER iz,in
c     include 'vdeform.inc'
      real uexcit, tdir, tff(2*numbar), phase(2*numbar)
      real tdirv(2*numbar)
C-------------------------------------------------------------------
c tdirv : WKB variable
c tdir  : WKB variable
c tff : fission transmission coefficient
c
      real Uexc, Smiu
      INTEGER K
      COMMON /VARGS/ Uexc, Smiu, K
c     real pi
c     PARAMETER (pi = 3.14159259d0)
      real dmom,abserr,rmiu,epsa,epsb,dummy,eps,deps,vdef,phasecal
      INTEGER j,ieps
c
c     FUNCTIONS
c
c vdef : deformed potential
c Vdef : deformed potential
c rmiu: parameter for WKB
c uexc : excitation energy
c uexc0: excitation energy
c uexc1: excitation energy
c uheigth: height of parabola
c uwidth: height of parabola
c uexcit: excitation energy
c
      real   Fmoment, FmomentParab, GaussLegendre41, FindIntersect
      EXTERNAL FmomentParab, Fmoment
      Uexc   = uexcit
      rmiu = 0.054d0*(iz+in)**(5.d0/3.d0)
C     Below is an square root of (MIU divided by 2)
      smiu = sqrt(0.5d0*rmiu)
      DO k=1, 2*numbar ! barriers and wells
        phase(k) = 0.d0
        tff(k)   = 0.d0
        tdirv(k) = 0.d0
      enddo
C-----------------------------------------------------------------------
C     Momentum integrals are calculated
      if (flagfisout) then
        write(22,*)
        write(22,*) ' Excitation Energy : ',uexcit
      endif
      DO k=1, nextr ! barriers and wells
          if(mod(k,2).eq.1) then
C         BARRIERS
            if(uexcit.ge.Vheight(k)) then
C
C           For excitation energies above barriers the transmission
C           coefficient is calculated by Hill-Wheeler formula
C
c dmom: momentum integral
c abserr: absolute error
c
            if (Vwidth(k).gt.0) then
              dmom = pi*(Vheight(k) - uexcit)/Vwidth(k)
            else
              dmom=-50.
            endif
            phase(k)   = min(dmom,50.)
            tff(k) = 1.d0/(1.d0 + DEXP(2.d0*dmom))
            if (flagfisout)
     &      write(22,'(1x,A6,I2,A10,d10.3,A3,d10.3,A15)')
     &      ' BARR ',k,'  Mom.Int=',dmom,' T=',tff(k),' (Hill-Wheeler)'

            else
C
C           For excitation energies below barriers the transmission
C           coefficient is calculated by WKB method
C
c epsa: help variable
c epsb: help variable
c deps: help variable
c phasecal: help variable
c dummy: help variable
c ieps: counter
c GaussLegendre41: function for Gauss Legendre integration
c Fmoment: function for integrand for Gauss-Legendre integration
c FmomentParab: function for integrand for Gauss-Legendre integration
c
            epsa = FindIntersect(uexcit,iiextr(k-1),iiextr(k),.false.)
            epsb = FindIntersect(uexcit,iiextr(k),iiextr(k+1),.false.)
C
C           Calculating phase integral for real shape
C
            dmom = GaussLegendre41(Fmoment,epsa,epsb,abserr)
            if(flagfisout.and.dmom.gt.0.d0.and.abserr.gt.dmom*0.03) then
              write(22,*) ' WARNING: For extremum ',k,
     &                   ' phase integral is not accurate (',
     &        abserr/dmom*100.d0,' %)'
            endif
            phase(k) = min(dmom,50.)
            tff(k) = 1.d0/(1.d0 + DEXP(2.d0*phase(k)))
C
            if (flagfisout) then
         write(22,'(1x,A6,I2, A10,f7.4,A4,f7.4,1x,A9,d10.3,A3,d10.3)')
     &      ' BARR ',k,' at beta2 ',epsa,' to ',epsb,
     &      ' Mom.Int=',phase(k),' T=',tff(k)
         deps=(epsb-epsa)/50.
         eps=epsa-deps
         phasecal=0.
         do ieps=1,50
           eps=eps+deps
           if (mod(ieps,5).eq.0)
     &     write(22,'(10x,"eps=",f6.3," Vdef-E=",f8.3)') eps,
     &       vdef(eps)-Uexc
           if (vdef(eps).ge.Uexc.and.vdef(eps+deps).ge.Uexc)
     &     phasecal=phasecal+2.*smiu*((vdef(eps)-Uexc)**0.5+
     &       (vdef(eps+deps)-Uexc)**0.5)/2.*deps
         enddo
         if (phasecal.lt.30.)
     &     write(22,'(10x," Kcal=",f10.4," Tcal=",es12.4)')
     &     phasecal,1./(1.+exp(2.*phasecal))
         endif
            endif
          else
C         WELLS
            if(uexcit.LE.Vheight(k)) then
C
C           Excitation energies below the well
C
            phase(k) = 0.d0
            else
C
C           For excitation energies above the well the transmission
C           coefficient is calculated by WKB method
C
            epsa = FindIntersect(uexcit,iiextr(k-1),iiextr(k),.true.)
            epsb = FindIntersect(uexcit,iiextr(k),iiextr(k+1),.true.)
C
C           Calculating phase integral for real shape
C
            dmom = GaussLegendre41(Fmoment,epsa,epsb,abserr)
         if (flagfisout.and.dmom.gt.0.d0.and.abserr.gT.dmom*0.03)
     &        write(22,*) ' WARNING: For extremum ',k,
     &        ' phase integral is not accurate (',
     &        abserr/dmom*100.d0,' %)'
            phase(k) = min(dmom,50.)
            if (flagfisout)
     &    write(22,'(1x,A6,I2, A10,f7.4,A4,f7.4,1x,A9,d10.3,A3,d10.3)')
     &      ' WELL ',k,' at beta2 ',epsa,' to ',epsb,
     &      ' Mom.Int=',phase(k)
            endif
          endif
      ENDDO
C
C     Fission transmission for double/triple humped barrier
C
C     Iteration over barriers
C
      if (nextr.gt.0) tdirv(nextr) = tff(nextr)
      do k=nextr-2,1,-2
         dmom = (1.d0 - tff(k))*(1.d0 - tdirv(k+2))
         if(k.gt.1) then
           tdirv(k) = tff(k)*tdirv(k+2) /
     &     (1.d0 + dmom)
c     &     (1.d0 + 2.d0*sqrt(dmom)*cos(2.d0*phase(k+1)) + dmom)
         else
           tdirv(1) = tff(1)*tdirv(3) /
     &     (1.d0 + dmom)
c     &     (1.d0 + 2.d0*sqrt(dmom)*cos(2.d0*phase(2)) + dmom)
         endif
      enddo
        tdir = tdirv(1)
        dummy = 1.d0/(1.d0/tff(1) + 1.d0/tff(3))
        if(nextr.gt.3)
     &  dummy = 1.d0/(1.d0/tff(1) + 1.d0/tff(3) + 1.d0/tff(5))

      if (flagfisout) write(22,'(1x,f5.2,2x,21(d10.3,1x))')
     &   uexc,tdir,dummy,(tff(j),j=1,nextr),(phase(j),j=1,nextr)
      RETURN
      END
      function VdefParab(EPS)
C
C     This function calculates parabolic shape of the deformation energy
C
C     Called by gaussian integration
C
c     IMPLICIT NONE
      include "talys.cmb"
      real VdefParab
      real EPS
c     INCLUDE 'vdeform.inc'
c VdefParab: function for parabolic shape
c Smiu: help variable
c Uexc: excitation energy
c
      real Uexc, Smiu
      INTEGER K
      COMMON /VARGS/ Uexc, Smiu, K
      VdefParab = Vheight(K) + (-1)**K*(SMIU*Vwidth(K)*(EPS-Vpos(K)))**2
      return
      end
      function Vdef(EPS)
C
C     This function calculates real shape of the deformation energy
C     by linear interpolation to obtain the value of the barrier Vdef
C     at deformation EPS (needed to integrte this function)
C
C     Called by gaussian integration
C
c     IMPLICIT NONE
c EPS : deformation
c
      include "talys.cmb"
      real Vdef
      real EPS
c     INCLUDE 'vdeform.inc'
C
C     Local variables
c ei: energy
c eip: energy
c vi: potential
c vip: potential
c
      INTEGER idef
      real vi, vip, ei, eip
C
C     The four lines below is a very simple (and inefficient) search
C     Should be replaced by efficient search routine to locate the
C     element idef of the array betafis()
      idef=1
      do while (EPS.GT.betafis(idef) .and. idef.LE.nbeta)
        idef = idef + 1
      enddo
      if (idef.ne.1) idef = idef - 1
      vi  = vfis(idef)
      vip = vfis(idef+1)
      ei  = betafis(idef)
      eip = betafis(idef+1)
      if(ei.eq.eip) then
C       Special case treated here to avoid division by zero
C       We assume that in this case vi = vip
        Vdef = vi
        return
      endif
      Vdef = vi + (EPS-ei)/(eip-ei)*(vip-vi)
      return
      end
      FUNCTION FindIntersect(uexc,ja,jb,iswell)
C
C     Calculates solutions (beta2(i)) of the equation
C         V(beta2) = Excitation Energy
C
C     If solution not found, then assign first or last point
C     of the interval depending on whether we are solving the equation
C     in the right(left) side of the barrier(well) case
C
c FindIntersect: function to find intersection
c ja: parameter for intersection
c jb: parameter for intersection
c iswell: logical
c
c
      include "talys.cmb"
c     IMPLICIT NONE
      INTEGER ja, jb
      real uexc
      real FindIntersect
      LOGICAL iswell
c     INCLUDE 'vdeform.inc'
C
C     Local variables
c is0: index
c
      INTEGER j, is0, is1
      real slope
c
c slope: slope
c
      is0 = -1
      IF(ABS(uexc-vfis(ja)).EQ.uexc-vfis(ja)) is0 = 1
      DO j=ja,jb
C      checking for a sign change in Uexc-Vdef
       is1 = -1
       IF(ABS(uexc-vfis(j)).EQ.uexc-vfis(j)) is1 = 1
       IF(is1.EQ.is0) CYCLE
C      Sign of (Uexc-vfis(j)) changed, calculating the precise value
C      of the deformation EPS at which Uexc = Vdef
       FindIntersect = betafis(j-1) + (betafis(j)-betafis(j-1))*
     >       (uexc-vfis(j-1))/(vfis(j)-vfis(j-1))
       RETURN
      ENDDO
C
C     Below is the analysis if intersection not found in [ja,jb]
C
      slope = vfis(jb) - vfis(ja)
      IF(iswell) then
C       WELLS
        IF(slope.ge.0) then
          FindIntersect = betafis(jb) ! ascending
        ELSE
          FindIntersect = betafis(ja) ! descending
        ENDIF
      ELSE
C       BARRIERS
        IF(slope.ge.0) then
          FindIntersect = betafis(ja) ! ascending
        ELSE
          FindIntersect = betafis(jb) ! descending
        ENDIF
      ENDIF
      RETURN
      END
      FUNCTION Find_Extrem (Nsmooth)
C
C     Find all extremums of the smoothed deformation energy curve
C
C     Nsmooth - Number of smoothing points
C
c
      include "talys.cmb"
c     IMPLICIT NONE
      integer Find_Extrem
      INTEGER Nsmooth
c     INCLUDE 'vdeform.inc'
      LOGICAL logmin, logmax
      INTEGER j,k
      INTEGER iext
      iext=0
      Find_Extrem = 0
      iiextr(0)=1
CC----------------------------------------------------------------------
C     determination of the minima   and maxima
      do j=nsmooth+1 , nbeta-nsmooth

c Modification SCAMPS to prevent problem (ex Fm324):
        if (vfis(j).eq.vfis(j-1)) cycle
        logmax=.true.
        do k=j-nsmooth,j+nsmooth
          if (k.eq.j) cycle
          if (vfis(k).gt.vfis(j)) logmax=.false.
        enddo
        if (logmax) then
          iext=iext+1
          iiextr(iext)=j
        endif
        if (iext.eq.2*numbar-1) exit
        logmin=.true.
        do k=j-nsmooth,j+nsmooth
          if (k.eq.j.or.k.lt.1) cycle
          if (vfis(k).lt.vfis(j)) logmin=.false.
        enddo
        if (logmin) then
          iext=iext+1
          iiextr(iext)=j
        endif
        if (iext.eq.2*numbar-1) exit
      enddo
      Find_Extrem = iext
      iiextr(iext+1)= nbeta
      return
      end

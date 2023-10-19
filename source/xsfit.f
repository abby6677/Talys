      subroutine xsfit(Z,A)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 8, 2021
c | Task  : Adjusted parameters to fit cross sections
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      lexist
      character*90 xsfile
      integer      Z,A,iz,ia,i
      real         x(10)
c
c ***************** Read parameters for (n,a) cross sections ***********
c
      x=0.
      if (flagnafit) then
        xsfile=trim(path)//'best/na.par'
        inquire (file=xsfile,exist=lexist)
        if (.not.lexist) goto 20
        open (unit=2,file=xsfile,status='old')
   10   read(2,'(2i4,5f12.5)',err=20,end=20) iz,ia,(x(i),i=1,5)
        if (Z.eq.iz.and.A.eq.ia) then
          rvadjust(6)=x(1)
          avadjust(6)=x(2)
          ompadjustp(6)=.true.
          Cstrip(6)=x(3)
          ctable(2,2,0)=x(4)
          ptable(2,2,0)=x(5)
        else
          goto 10
        endif
   20   close (unit=2)
      endif
      if (flagnnfit) then
        xsfile=trim(path)//'best/nn.par'
        inquire (file=xsfile,exist=lexist)
        if (.not.lexist) goto 120
        open (unit=2,file=xsfile,status='old')
  110   read(2,'(2i4,10f12.5)',err=120,end=120) iz,ia,(x(i),i=1,10)
        if (Z.eq.iz.and.A.eq.ia) then
          rvadjust(2)=x(1)
          avadjust(2)=x(2)
          rwadjust(2)=x(1)
          awadjust(2)=x(2)
          ompadjustp(2)=.true.
          gnadjust(0,0)=x(3)
          gpadjust(0,0)=x(4)
          ctable(0,1,0)=x(5)
          ptable(0,1,0)=x(6)
          ctable(0,2,0)=x(7)
          ptable(0,2,0)=x(8)
          ctable(1,0,0)=x(9)
          ptable(1,0,0)=x(10)
        else
          goto 110
        endif
  120   close (unit=2)
      endif
      return
      end
Copyright (C)  2021 A.J. Koning, S. Hilaire and S. Goriely

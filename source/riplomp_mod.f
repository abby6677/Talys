      subroutine riplomp_mod(NE,Zt,At,Eripl,index1,index2,index3)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : May 25, 2020
c | Task  : Interface for RIPL OMP 
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      USE om_retrieve
      integer NE
      real Eripl(NE)
      integer Zt,At,index1,index2,index3
      Number_Energies=NE
      Ztarget=Zt
      Atarget=At
      Energies=0.
      do i=1,NE
        Energies(i)=Eripl(i)
      enddo
      RIPL_Index=index1
      Calc_Type=index2
      Def_Index=index3
      call retrieve
      end
Copyright (C)  2020 A.J. Koning, S. Hilaire and S. Goriely

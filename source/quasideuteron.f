      function quasideuteron(Egamma)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 22, 2017
c | Task  : Quasi-deuteron model of Chadwick and Oblozinsky
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real    quasideuteron,Egamma,levinger,freedeut,fpauli
c
c ****** Calculate quasi-deuteron photo-absorption cross section *******
c
c levinger: Levinger parameter
c Egamma  : gamma energy
c freedeut: free deuteron cross section
c fpauli  : Pauli-blocking function of Chadwick
c xsqd    : photo-absorption cross section from QD part
c Ntarget : neutron number of target nucleus
c Ztarget : charge number of target nucleus
c Atarget : mass number of target nucleus
c
      levinger=6.5
      if (Egamma.gt.2.224) then
        freedeut=61.2/(Egamma**3)*(Egamma-2.224)**1.5
      else
        freedeut=0.
      endif
      if (Egamma.le.140.) then
        if (Egamma.ge.20.) then
          fpauli=8.3714e-2-9.8343e-3*Egamma+4.1222e-4*Egamma*Egamma-
     +      3.4762e-6*(Egamma**3)+9.3537e-9*(Egamma**4)
        else
          fpauli=exp(-73.3/Egamma)
        endif
      else
        fpauli=exp(-24.2348/Egamma)
      endif
      quasideuteron=levinger*Ntarget*Ztarget/real(Atarget)*
     +  freedeut*fpauli
      return
      end
Copyright (C)  2017 A.J. Koning, S. Hilaire and S. Goriely

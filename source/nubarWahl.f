      function nubarWahl(Ein)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : September 12, 2021
c | Task  : Wahl systematics for nubar
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      real    nubarWahl,Ein,Fz,Fn,term1,term2,Bn,Th
c
c ************************ Wahl systematics ****************************
c
c nubar  : average number of fission neutrons
c Ein    : incident energy
c Bn     : neutron separation energy
c Zinit  : charge number of initial compound nucleus
c Ainit  : mass number of initial compound nucleus
c Ninit  : neutron number of initial compound nucleus
c Fz,Fn  : help variables
c term1  : help variable
c Cnubar1: adjustable parameter for nubar constant value
c Cnubar2: adjustable parameter for nubar energy slope
c
c A.C. Wahl, "Systematics of fission-product yields", LA-13926 (2002), 
c Eq. (21)
c
      Bn=S(0,0,1)
      if (mod(Zinit,2).eq.0) then
        Fz=1.
      else
        Fz=-1.
      endif
      if (mod(Ninit,2).eq.0) then
        Fn=1.
      else
        Fn=-1.
      endif
      term1=2.286+0.147*(Zinit-92)+0.054*(Ainit-236)+0.040*(2.-Fz-Fn)
      if (Ein.gt.0.) then
        Th=11.47-0.166*Zinit*Zinit/real(Ainit)+0.093*(2.-Fz-Fn)-Bn
        term2=(0.145-0.0043*(Ainit-236))*(Ein-Th)
      else
        term2=0.
      endif
      nubarWahl=Cnubar1*term1+Cnubar2*term2
      return
      end
Copyright (C)  2021 A.J. Koning, S. Hilaire and S. Goriely

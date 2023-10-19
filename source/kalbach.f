      function kalbach(type,Ein,Eout,ang)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 10, 2015
c | Task  : Kalbach systematics
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer type,i
      real    kalbach,Ein,Eout,ang,Et1,Et3,C1,C2,C3,Cmin(numpar),
     +        Cmout(numpar),ea,Ek1,Ek3,eb,X1,X3,akal
c
c ************************ Kalbach systematics *************************
c
c kalbach   : Kalbach function
c Ein       : incident energy
c Eout      : outgoing energy
c ang       : angle
c Et1       : constant of Kalbach systematics
c Et3       : constant of Kalbach systematics
c Cmin      : constant of Kalbach systematics
c Cmout     : constant of Kalbach systematics
c C1        : constant
c C2        : constant
c X1        : energy
c X3        : energy
c Ek1       : energy
c Ek3       : energy
c C3        : constant
c Smyers    : Myers-Swiatecki separation energy
c k0        : index of incident particle
c akal      : Kalbach 'a' parameter
c fourpi    : 4.*pi
c
c  Systematics of Kalbach: Phys. Rev. C37, 2350, (1987)
c
c Since we separate the pre-equilibrium (xspreeqad) and compound
c (xscompad) angular distribution in the output we only need to take
c the forward peaked component of the Kalbach formula to calculate the
c pre-equilibrium angular distribution.
c
      data (Cmin(i),i=1,6) /1.,1.,1.,1.,1.,0./
      data (Cmout(i),i=1,6) /0.5,1.,1.,1.,1.,2./
      Et1=130.
      Et3=41.
      C1=0.04
      C2=1.8e-6
      C3=6.7e-7
c
c Isotropic distribution for photons
c
      if (k0.eq.0.or.type.eq.0) then
        kalbach=1./fourpi
      else
        ea=Ein+Smyers(k0)
        Ek1=min(ea,Et1)
        Ek3=min(ea,Et3)
        eb=Eout+Smyers(type)
        X1=Ek1*eb/ea
        X3=Ek3*eb/ea
        akal=C1*X1+C2*(X1**3)+C3*Cmin(k0)*Cmout(type)*(X3**4)
        kalbach=0.
        if (abs(akal).le.80.) kalbach=1./fourpi*akal/sinh(akal)*
     +    exp(akal*cos(ang))
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

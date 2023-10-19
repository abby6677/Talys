      subroutine dwbaread(nen1,nen2)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning and Eric Bauge
c | Date  : August 10, 2015
c | Task  : Read ECIS results for DWBA for MSD
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer          nen1,nen2,J,nS,iang,k,itype,istat
      double precision xs
c
c ********************** Read DWBA cross sections **********************
c
c maxJmsd   : maximal spin for MSD calculation
c xsdwin    : DWBA cross section as a function of incident energy,
c             outgoing energy and angular momentum
c flagddx   : flag for output of double-differential cross sections
c nanglecont: number of angles for continuum
c xsdw      : DWBA angular distribution as a function of incident
c             energy, outgoing energy, angular momentum and angle
c nS        : counter
c
      read(10,'()')
      do 10 J=0,maxJmsd
        read(10,*) xs
        xsdwin(nen1,nen2,J,0)=real(xs)
   10 continue
      if (flagddx) then
        read(8,'()')
        read(8,'(12x,i3)') nS
        do 30 iang=0,nanglecont
          do 30 k=1,nS
            read(8,'()')
   30   continue
        do 40 J=0,maxJmsd
          read(8,'(12x,i3)',iostat=istat) nS
          if (istat.ne.0) goto 40
          do 50 iang=0,nanglecont
            do 60 k=1,nS
              read(8,'(i3,12x,e12.5)',iostat=istat) itype,xs
              if (istat.ne.0) goto 60
              if (itype.eq.0) xsdw(nen1,nen2,J,iang,0)=real(xs)
   60       continue
   50     continue
   40   continue
      endif
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

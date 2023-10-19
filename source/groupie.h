C=======================================================================
C
C     GROUPIE COMMON
C
C=======================================================================
C
C     PARAMETERS
C
C-----------------------------------------------------------------------
c-----2017/3/07 - INCREASED TO 3,000,000 FROM 600,000.
CAK   PARAMETER (MAXPOINT = 3000000) ! Data points in memory
c-----REPORT FOR EACH MAT
      PARAMETER (MAXMAT  =  1000)
c-----MAX. NUMBER OF BANDS
      PARAMETER (MAXBAND = 5)
C-----MAX. NUMBER OF GROUPS
c-----2015/12/22 - Increased to 20,000 from 1,000
CAK   PARAMETER (MAXGROUP = 20000)
c-----MULTI-GROUP OUTPUT TO ENDF/B FORMAT
c-----2015/12/22 - Increased to 20,000 from 3,000
CAK   PARAMETER (MAXOUT   = 20000)
CAK Reduced to previous values to make PREPRO routines fit inside TALYS
      PARAMETER (MAXPOINT = 600000) ! Data points in memory
      PARAMETER (MAXGROUP = 1000)
      PARAMETER (MAXOUT   = 3000)
C-----------------------------------------------------------------------
C
C     STORAGE
C
C-----------------------------------------------------------------------
c
C     PAGED STORAGE FOR
C     1 = SPECTRUM
C     2 = TOTAL
C     3 = ELASTIC
C     4 = CAPTURE
C     5 = FISSION
C     6 = OTHER - NOT IN PAGING SYSTEM, ONLY IN MEMORY ARRAYS
C
C-----------------------------------------------------------------------
      COMMON XPAGE(MAXPOINT,5),YPAGE(MAXPOINT,5)
      COMMON/REPORT/ERBTAB(6,MAXMAT),ERLIB(25,MAXBAND),NBNTAB(MAXMAT),
     1 IZATAB(MAXMAT),LZA
      COMMON/GROUPR/EGROUP(MAXGROUP+1),TOTAV(MAXGROUP+1),NGR,NGRP1,IGR
C-----PLOTTING COMMON = UNSHIELDED AND SHIELDED
      COMMON/PLOTCOM/XCPLOT(MAXGROUP,5,2)
C-----MULTI-BAND STORAGE
      COMMON/BANDID/WTBAND(6,MAXBAND,MAXGROUP),
     1              XCBAND(6,MAXBAND,MAXGROUP)
      COMMON/INTNRM/XCINT(25,6),XCNORM(25)
      COMMON/MATERR/ERMAT(25,MAXBAND),ERNOW(25,MAXBAND)
      COMMON/SIMPLE/XCFI(25,6),AVEXP(25),AVNORM(25),SIGMAB(25)
      COMMON/MISSIT/FST(25),DFST(25),DDFST(25),DEAVST(25)
      COMMON/ELPAS1/SHIELD(25),YLOW(5),YHIGH(5),YLOWP1(5),YLOWP2(5),
     1 YLOWP3(5),YHIGHP1(5),YHIGHP2(5),YHIGHP3(5)
      CHARACTER*4 POINT,REACT2,REACT3
      COMMON/ELPASC/POINT(25),REACT2(2,6),REACT3(2,6)
C-----ENDF Formatted Output.
      COMMON/GROUPOUT/XOUT(MAXOUT),YOUT(MAXOUT)

C=======================================================================
C
C     RECENT COMMON
C
C=======================================================================
C
C     PARAMETERS
C
C-----------------------------------------------------------------------
C-----Maximum Reactions
      PARAMETER (MAXNREACT = 11)
C-----Max. points in tabulated Rho
      PARAMETER (MAXRHO = 5000)
c-----2014/2/20 - reduced to 300000 from 600000
      PARAMETER (MAXRES = 300000) ! Max. # of resonances
c-----2014/2/1 - decreased from 1,000 to 100
      PARAMETER (MAXSEC =  100)   ! Max. # of sections
C-----2017/3/7 - INCREASED FROM 600,000 TO 1,200,000
C
C     WARNING - Change both MAXPTX and MAXPTXP1
C
CAK   PARAMETER (MAXPTX   = 1200000) ! Max. # energy points per page
CAK   PARAMETER (MAXPTXP1 = 1200001) ! MAXPTX + 1
CAK Reduced to previous values to make PREPRO routines fit inside TALYS
      PARAMETER (MAXPTX   = 600000) ! Max. # energy points per page
      PARAMETER (MAXPTXP1 = 600001) ! MAXPTX + 1
C-----12/20/06 - INCREASED MAXSAVE FROM 200 TO 2,000
      PARAMETER (MAXSAVE = 2000)  ! Max. # of saved iteration points
C-----------------------------------------------------------------------
C
C     STORAGE
C
C-----------------------------------------------------------------------
C-----TABULATED ENERGY DEPENDENT RHO
      COMMON/TABRHO/ERHOTB(MAXRHO),RHOTAB(MAXRHO),APTAB(MAXRHO),
     1 NRHO(MAXSEC),NBTRHO(MAXSEC),INTRHO(MAXSEC),INXRHO(4,MAXSEC),
     2 NUMRHO
C-----RESONANCE SECTION
      COMMON/SECTON1/BETA(MAXSEC),RHOX2(MAXSEC),RHOC2(MAXSEC),
     1 RHOP1(MAXSEC),GJTAB(MAXSEC),QVALUE(MAXSEC),EXCITE(4,MAXSEC),
     2 NLOW(MAXSEC),NHIGH(MAXSEC),LVALUE(MAXSEC),LVALUEC(MAXSEC),
     3 LRXTAB(MAXSEC),NAPTAB(MAXSEC),NGRTAB(2,MAXSEC),NPPTAB(MAXSEC)
      COMMON/SECTON2/EL(MAXSEC),EH(MAXSEC),POTL(MAXSEC),ADDL(MAXSEC),
     1 MODE(MAXSEC),ISECT,NSECT
c-----2016/11/15 - Added to allow L dependent fission = speed if no
      COMMON/GOFISH/LFWL(MAXSEC)
C-----2017/4/13 - INTERPOLATION LAW - now fixed dimension = 100
      COMMON/TERPCOM/INTF(100),NBTF(100)
C-----------------------------------------------------------------------
C
C     BLANK COMMON = Only 1 declaration
C
C     RESONANCE DATA
C
C     CROSS SECTION ARRAYS
C     --------------------
C     ENRES   = Resonance Energy
C     SHIFT2  = Shift Factor
C     ENODE   = Starting node energy
C     WIDNOD  = Starting node width
C     RESTAB  = Resonance parameters
C     RESJTAB = Resonance J values
C     ETAB2   = FILE 2 RESONANCE  ENERGY
C     SIG2    = FILE 2 RESONANCE  CROSS SECTION
C     ETAB3   = FILE 3 BACKGROUND ENERGY
C     SIG3    = FILE 3 BACKGROUND CROSS SECTION
C     ETAB23  = COMBINED FILE 2 + 3 ENERGY
C     SIG23   = COMBINED FILE 2 + 3 CROSS SECTION
C     SIGMID  = CROSS SECTION AT MIDPOINT OF ITERATION INTERVAL
C
C-----------------------------------------------------------------------
C-----08/29/09 - EXPANDED RESTAB FROM 6 TO 12 PARAMETERS PER RESONANCE
C-----           TO HANDLE LRF=7 REICH-MOORE.
C----- WARNING - DO NOT DECREASE BELOW 12, TO HANDLE ADLER-ADLER
c-----2017/4/6 - Combined 2 blank COMMON groups
c-----2017/4/6 - Deleted RESTABER
      COMMON ENRES(MAXRES),SHIFT2(MAXRES),ENODE(MAXRES),WIDNOD(MAXRES),
     1 RESTAB(12,MAXRES),RESJTAB(MAXRES),
     1       ETAB2 (MAXPTXP1),SIG2 (MAXNREACT,MAXPTXP1),
     1       ETAB3 (MAXPTX),  SIG3 (MAXPTX),
     1       ETAB23(MAXPTX),  SIG23(MAXPTX),SIGMID(MAXNREACT),
     1       ESAVE(MAXSAVE),SIGSAVE(MAXNREACT,MAXSAVE),NSAVE
C
C     NOTE THAT ETAB2 AND SIG2 ARE DIMENSIONED MAXPTX+1, WHILE ETAB2X
C     AND SIG2X ARE DIMENSIONED MAXPTX - TO ALLOW PAGES TO BE WRITTEN
C     OR READ WHILE STILL KEEPING 1 ENERGY POINT IN MEMORY.
 
      DIMENSION ETAB2X(MAXPTX),SIG2X(MAXNREACT,MAXPTX)
      EQUIVALENCE (ETAB2(1) ,ETAB2X(1)),(SIG2(1,1),SIG2X(1,1))

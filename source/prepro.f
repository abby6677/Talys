        subroutine prepro(infile,outfile,Tres,flaggroup)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : July 8, 2018
c | Task  : PREPRO routines for making cross section pointwise
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
c recentt: TALYS version of RECENT code of PREPRO package
c
      logical      flaggroup
      character*72 infile,outfile,inrecent,insigma1,ingroupie,outrecent,
     +             outsigma1,outgroupie
      real         Tres
c
c Call required modules of prepro package
c
c recentt : TALYS version of RECENT code of PREPRO package
c sigma1t : TALYS version of SIGMA1 code of PREPRO package
c groupiet: TALYS version of GROUPIE code of PREPRO package
c
      inrecent=infile
      outrecent=outfile(1:11)//'point0'
      call recentt(inrecent,outrecent)
      if (Tres.gt.0.) then
        insigma1=outrecent
        outsigma1=outrecent(1:11)//'point'
        call sigma1t(insigma1,outsigma1,dble(Tres))
      endif
      if (flaggroup) then
        if (Tres.gt.0.) then
          ingroupie=outsigma1
        else
          ingroupie=outrecent
        endif
        outgroupie=outfile
        call groupiet(ingroupie,outgroupie)
      endif
      return
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely
      subroutine recentt(infile,outfile)
C=======================================================================
C
C     PROGRAM RECENT
C     ==============
C     VERSION 79-1 (OCTOBER 1979)  CDC-7600
C     VERSION 80-1 (MAY 1980)      IBM, CDC AND CRAY VERSION
C     VERSION 80-2 (DECEMBER 1980) IMPROVED TREATMENT OF UNRESOLVED
C                                  REGION TO COMPUTE ALL REACTIONS AT
C                                  THE SAME TIME.
C     VERSION 81-1 (MARCH 1981)  IMPROVED BASED ON USER COMMENTS.
C     VERSION 81-2 (AUGUST 1981) ADDED MONITOR MODE. ADDED SPEED OPTION
C                                TO BYPASS BACKWARDS THINNING IF FILE 3
C                                ALLOWABLE ERROR = 0.0 (NOTE THIS OPTION
C                                WILL RESULT IN ALL TABULATED POINTS
C                                FROM THE EVALUATION BEING KEPT IN THE
C                                OUTPUT FROM THIS PROGRAM).
C     VERSION 82-1 (JANUARY 1982) IMPROVED COMPUTER COMPATIBILITY.
C     VERSION 83-1 (JANUARY 1983)*MAJOR RE-DESIGN.
C                                *PAGE SIZES INCREASED.
C                                *ELIMINATED COMPUTER DEPENDENT CODING.
C                                *NEW, MORE COMPATIBLE I/O UNIT NUMBERS.
C                                *ADDED OPTION TO KEEP ALL RECONSTRUCTED
C                                 AND BACKGROUND ENERGY POINTS.
C                                *ADDED STANDARD ALLOWABLE ERROR OPTIONS
C                                 (CURRENTLY 0.1 PER-CENT RECONSTRUCTION
C                                 AND 0.0 PER-CENT THINNING).
C     VERSION 83-2 (OCTOBER 1983) IMPROVED BASED ON USER COMMENTS.
C     VERSION 84-1 (JANUARY 1984) IMPROVED INTERVAL HALFING CONVERGENCE.
C     VERSION 85-1 (APRIL 1985)  *A BRAND NEW PROGRAM WHICH COMPLETELY
C                                 SUPERCEDES ALL PREVIOUS VERSIONS OF
C                                 THIS PROGRAM.
C                                *UPDATED FOR ENDF/B-VI FORMATS.
C                                *ADDED GENERAL REICH-MOORE FORMALISM
C                                 (WITH TWO FISSION CHANNELS).
C                                *DECREASED RUNNING TIME.
C                                *SPECIAL I/O ROUTINES TO GUARANTEE
C                                 ACCURACY OF ENERGY.
C                                *DOUBLE PRECISION TREATMENT OF ENERGY
C                                 (REQUIRED FOR NARROW RESONANCES).
C     VERSION 85-2 (AUGUST 1985) *FORTRAN-77/H VERSION
C     VERSION 86-1 (JANUARY 1986)*ENERGY DEPENDENT SCATTERING RADIUS
C     VERSION 86-2 (JUNE 1986)   *IF FIRST CHANCE FISSION (MT=19)
C                                 BACKGROUND IS PRESENT ADD RESONANCE
C                                 CONTRIBUTION OF FISSION TO IT.
C     VERSION 86-3 (OCTOBER 1986)*MULTI-LEVEL OR REICH-MOORE..CORRECT
C                                 POTENTIAL SCATTERING CROSS SECTION FOR
C                                 MISSING AND/OR FICTICIOUS (L,J)
C                                 SEQUENCES.
C     VERSION 87-1 (JANUARY 1987)*IMPROVED COMBINING FILE 2+3
C     VERSION 87-2 (MARCH 1987)  *CORRECTED ADLER-ADLER CALCULATIONS.
C     VERSION 88-1 (JULY 1988)   *UPDATED REICH-MOORE ENDF/B-VI FORMAT
C                                 TO BE THE SAME AS REICH-MOORE FORMAT
C                                 IN EARLIER VERSIONS OF ENDF/B FORMAT.
C                                *CHECK FOR PRELIMINARY ENDF/B-VI
C                                 REICH-MOORE FORMAT (NOW ABANDONED)
C                                 AND TERMINATE EXECUTION IF DATA IS
C                                 IN THIS FORMAT.
C                                *CALCULATE CHANNEL RADIUS OR SET IT
C                                 EQUAL TO THE SCATTERING RADIUS.
C                                *IMPLEMENTED HYBRID R-FUNCTION WITH THE
C                                 FOLLOWING RESTRICTIONS
C                                 - ONLY INELASTIC COMPETITION (NO
C                                   CHARGED PARTICLES)
C                                 - NO TABULATED FILE 2 BACKGROUND
C                                 - NO TABULATED OPTICAL MODEL PHASE
C                                   SHIFT
C                                *PROGRAM EXIT IF GENERAL R-MATRIX IN
C                                 THE EVALUATION (THIS FORMALISM WILL
C                                 BE IMPLEMENTED ONLY AFTER THE AUTHOR
C                                 RECEIVES REAL EVALUATIONS WHICH USE
C                                 THIS FORMALISM...UNTIL THEN IT IS
C                                 IMPOSSIBLE TO ADEQUATELY TEST THAT
C                                 THE CODING FOR THIS FORMALISM IS
C                                 CORRECT).
C                                *INCREASED MAXIMUM NUMBER OF RESONANCES
C                                 FROM 1002 TO 4008.
C                                *DOUBLE PRECISION RESONANCE REGION
C                                 LIMITS.
C                                *FILE 2 AND FILE 3 ENERGIES WHICH ARE
C                                 NEARLY EQUAL ARE TREATED AS EQUAL
C                                 (I.E., SAME TO ABOUT 9 DIGITS).
C                                *CHECK FILE 3 BACKGROUND CROSS SECTIONS
C                                 IN EDIT MODE.
C                                *OPTION...INTERNALLY DEFINE FILENAMES
C                                 (SEE SUBROUTINE FILEIO FOR DETAILS).
C     VERSION 89-1 (JANUARY 1989)*PSYCHOANALYZED BY PROGRAM FREUD TO
C                                 INSURE PROGRAM WILL NOT DO ANYTHING
C                                 CRAZY.
C                                *UPDATED TO USE NEW PROGRAM CONVERT
C                                 KEYWORDS.
C                                *CORRECTED MULTILEVEL, REICH-MOORE AND
C                                 HYBRID R-FUNCTION POTENTIAL SCATTER
C                                 TO ACCOUNT FOR REPEATED J-VALUES FOR
C                                 THE SAME TARGET SPIN AND L-VALUE.
C                                *ADDED LIVERMORE CIVIC COMPILER
C                                 CONVENTIONS.
C                                *UPDATED TO USE NEW ENDF/B-VI
C                                 CONVENTION TO ALLOW UNRESOLVED
C                                 RESONANCE CONTRIBUTION TO ALREADY
C                                 BE INCLUDED IN THE FILE 3 CROSS
C                                 SECTIONS (INFINITELY DIULUTE
C                                 CONTRIBUTION).
C     VERSION 90-1 (JUNE 1990)   *UPDATED BASED ON USER COMMENTS
C                                *ADDED FORTRAN SAVE OPTION
C                                *NEW MORE CONSISTENT ENERGY OUTPUT
C                                 ROUTINE
C     VERSION 91-1 (JULY 1991)   *NEW UNIFORM TREATMENT OF ALL RESONANCE
C                                 FORMALISMS (SEE, COMMENTS BELOW)
C                                *NEW REICH-MOORE ALGORITHM
C                                *MORE EXTENSIVE ERROR CHECKING AND
C                                 ERROR MESSAGE EXPLANATIONS
C     VERSION 92-1 (JANUARY 1992)*MAJOR RESTRUCTING TO IMPROVE ACCURACY
C                                 AND COMPUTER INDEPENDENCE.
C                                *INCREASED ENERGY POINT PAGE SIZE FROM
C                                 1002 TO 4008.
C                                *NO MORE THAN 2 ENERGY POINTS WHERE
C                                 CROSS SECTION IS ZERO AT BEGINNING
C                                 OF A SECTION FOR EACH REACTION,E.G.,
C                                 THRESHOLD FISSION.
C                                *PROCESS ONLY A PORTION OF RESONANCE
C                                 REGION - SEE EXPLANATION BELOW
C                                *ALL ENERGIES INTERNALLY ROUNDED PRIOR
C                                 TO CALCULATIONS.
C                                *COMPLETELY CONSISTENT I/O AND ROUNDING
C                                 ROUTINES - TO MINIMIZE COMPUTER
C                                 DEPENDENCE.
C     VERSION 93-1 (MARCH 1993)  *UPDATED REICH-MOORE TREATMENT TO USE
C                                 L DEPENDENT SCATTERING RADIUS (APL)
C                                 RATHER THAN SCATTERING RADIUS (AP)
C                                 (SEE, ENDF/B-VI FORMATS AND
C                                  PROCEDURES MANUAL, PAGE 2.6)
C                                *INCREASED PAGE SIZE FROM 4008 TO
C                                 20040 DATA POINTS.
C                                *INCREASED MAXIMUM NUMBER OF RESONANCES
C                                 FROM 4008 TO 20040.
C     VERSION 94-1 (JANUARY 1994)*VARIABLE ENDF/B DATA FILENAMES
C                                 TO ALLOW ACCESS TO FILE STRUCTURES
C                                 (WARNING - INPUT PARAMETER FORMAT
C                                 HAS BEEN CHANGED).
C                                *CLOSE ALL FILES BEFORE TERMINATING
C                                 (SEE, SUBROUTINE ENDIT)
C     VERSION 94-2 (AUGUST 1994) *CORRECTED ADDJ FOR ENERGY DEPENDENT
C                                 (TABULATED) SCATTERING RADIUS CASE.
C     VERSION 96-1 (JANUARY 1996) *COMPLETE RE-WRITE
C                                 *IMPROVED COMPUTER INDEPENDENCE
C                                 *ALL DOUBLE PRECISION
C                                 *ON SCREEN OUTPUT
C                                 *UNIFORM TREATMENT OF ENDF/B I/O
C                                 *IMPROVED OUTPUT PRECISION
C                                 *ALWAYS INCLUDE THERMAL VALUE
C                                 *DEFINED SCRATCH FILE NAMES
C     VERSION 97-1 (APRIL 1997)   *OPTIONAL MAKE NEGATIVE CROSS
C                                  SECTION = 0 FOR OUTPUT
C                                *INCREASED PAGE SIZE FROM 20040 TO
C                                 120000 DATA POINTS.
C                                *INCREASED MAXIMUM NUMBER OF RESONANCES
C                                 FROM 20040 TO 120000.
C     VERSION 99-1 (MARCH 1999)   *CORRECTED CHARACTER TO FLOATING
C                                  POINT READ FOR MORE DIGITS
C                                 *UPDATED TEST FOR ENDF/B FORMAT
C                                  VERSION BASED ON RECENT FORMAT CHANGE
C                                 *UPDATED CONSTANTS BASED ON CSEWG
C                                  SUBCOMMITTEE RECOMMENDATIONS
C                                 *GENERAL IMPROVEMENTS BASED ON
C                                  USER FEEDBACK
C     VERSION 99-2 (JUNE 1999)    *IMPLEMENTED NEW REICH-MOORE FORMALISM
C                                  TO ALLOW DEFINITION OF (L,J,S) FOR
C                                  EACH SEQUENCE.
C                                 *ASSUME ENDF/B-VI, NOT V, IF MISSING
C                                  MF=1, MT-451.
C     VERS. 2000-1 (FEBRUARY 2000)*GENERAL IMPROVEMENTS BASED ON
C                                  USER FEEDBACK
C     VERS. 2002-1 (MAY 2002)     *OPTIONAL INPUT PARAMETERS
C                  (SEPT. 2002)   *OUTPUT RESONANCE WITH 9 DIGITS
C                                 *TO BE C AND C++ COMPATIBLE OUTPUT
C     VERS. 2004-1 (JAN. 2004)    *ADDED INCLUDE 'recent.h'
C                                 *MADE ENDF/B-VII READY
C                                 *UPDATED FOR NEW REICH-MOORE LRF=7
C                                  PARAMETERS WITH COMPETITION
C                                 *ADDED COULOMB PENETRATION FACTORS FOR
C                                  LRF=7 COMPETITIVE CHANNELS.
C                                 *EXTENDED DEFINITIONS OF PENETRATION
C                                  FACTOR, LEVEL SHIFT FACTOR, AND
C                                  POTENTIAL SCATTERING PHASE SHIFT
C                                  ABOVE L = 5 TO INFINITY.
C                                 *ADDED QUICK CALCULATION - IF THE
C                                  INPUT ALLOWABLE ERROR IS 1.0 OR MORE
C                                  (100 % OR MORE) THERE IS NO ITERATION
C                                  TO CONVERGENCE - CROSS SECTION ARE
C                                  QUICKLY CALCULATED ONLY AT A FIXED
C                                  SET OF ENERGY POINTS, BASED ON THE
C                                  ENERGY AND WIDTH OF ALL RESONANCES.
C                                  THIS CAN BE USED TO QUICKLY "SEE"
C                                  NEW EVALUATIONS THAT MAY CONTAIN
C                                  ERRORS, THAT WOULD OTHERWISE CAUSE
C                                  THIS CODE TO RUN FOR AN EXCESSIVELY
C                                  LONG TIME.
C     VERS. 2005-1 (JUNE 2005)    *ADDED ENERGY DEPENDENT SCATTERING
C                                  RADIUS FOR ALL RESONANCE TYPES
C                                  (EARLIER ONLY BREIT-WIGNER ALLOWED).
C     VERS. 2007-1 (JAN. 2007)    *CHECKED AGAINST ALL ENDF/B-VII.
C                                 *DECOUPLED PAGE SIZE FROM MAX. # OF
C                                  RESONANCES.
C                                 *INCREASED PAGE SIZE FROM 120,000 TO
C                                  750,000 DATA POINTS.
C                                 *KEPT MAX. # OF RESONANCE AT 120,000.
C                                 *CORRECTED ALL BACKGROUND = 0 CASE
C     VERS. 2007-2 (OCT. 2007)    *NO MT=19 OUTPUT IF NO BACKGROUND,
C                                  REGARDLESS OF INPUT OPTION.
C                                 *72 CHARACTER FILE NAMES.
C     VERS. 2008-1 (FEB. 2008)    *CORRECTED NAPS ERROR - NOW DEFINE FOR
C                                  ALL TYPES OF PARAMETERS - EARLIER
C                                  ONLY DEFINED FOR B-W PARAMETERS.
C     VERS. 2008-2 (APRIL 2008)   *CORRECTED NRO/NAPS=1/1 - MUST
C                                  DEFINE RHOX2 AT EACH RESONANCE USING
C                                  SETRHO1 BEFORE ENERGY DEPENDENT
C                                  CALCULATION.
C                                 *ADDED PRECISION TO RESONANCE PROFILE
C                                  IN SUBROUTINE SUBINT
C     VERS. 2009-1 (JULY 2009)    *NEW REICH-MOORE COMPETITIVE WIDTHS -
C                                  IF CHARGED PARTICLE REACTION (MT=103
C                                  THROUGH 107) WILL ADD RESONANCE
C                                  CONTRIBUTION TO COMPETITIVE MT AND IF
C                                  PRESENT, THE GROUND LEVEL, MT = 600
C                                  THROUGH 800. IF COMPETITIVE CHANNEL
C                                  IS mt=4 (TOTAL N.N') IT WILL ALSO ADD
C                                  COMPETITIVE RESONANCE CONTRIBUTION TO
C                                  MT=50 (N,N' GROUND).
C                                 *NEW REICH-MOORE - SUM COMPETITIVE
C                                  WIDTHS IF ALL FOR THE SAME STATE (MT)
C     VERS. 2009-2 (AUG. 2009)    *RE-WRITE TO USE 12, RATHER THAN 6,
C                                  PAAMETERS PER RESONANCE.
C                                 *MAJOR RE-WRITE TO ACCOMODATE GENERAL
C                                  REICH-MOORE (LRF=7).
C                                 *COMPLETE RE-WRITE FOR ADLER-ADLER
C                                  AND HRF (N O LONGER USED IN ENDF/B)
C                                  TO USE 12 PARAMETERS PER RESNANCE.
C     VERS. 2010-1 (April 2010)   *ADDED SAMRML LOGIC TO HANDLE ALL
C                                  LRF=7 CASES.
C                                 *EXTENDED SAMRML LOGIC TO PROCESS ALL
C                                  EVALUATIONS = RESOLVED + UNRESOLVED +
C                                  TABULATED - SAMRML ONLY DOES ONE
C                                  SECTION OF RESOLVED LRF=7 DATA
C                                  WITHOUT TABULATED BACKGROUND.
C                                 *UPDATED ELASTIC POTENTIAL CALCULATION
C                                  FOR TOTAL (SLBW) AND CORRECTION FOR
C                                  MISSING SEQUENCES (MLBW, RM, HRF).
C                                 *ADDED HIDDEN (OPTIONAL) UNRESOLVED
C                                  COMPETITION LISTING (NOT ENDF/B).
C                                 *ADDED BOB MACFARLANE'S PROPOSAL - USE
C                                  LRX TO DEFINE COMPETITIVE L VALUE -
C                                  COMPETITIVE L = LRX - 1, IF LRX > 0.
C                                 *CHECKED FOR NEGATIVE WIDTHS.
C     VERS. 2012-1 (Nov.  2012)   *ADDED ENERGY DEPENDENT STEP SIZE
C                                  FOR STARTING GRID AROUND RESONANCES.
C                                 *Added CODENAME
C                                 *32 and 64 bit Compatible
C                                 *Added ERROR stops
C                                 *Check for no capture for Reich-Moore.
C     VERS. 2012-2 (Nov.  2012)   *Eliminated ERROR in NHIGH(0) index.
C     VERS. 2013-1 (Nov.  2013)   *Extended OUT9.
C     VERS. 2015-1 (Jan.  2015)   *Multiple LRF=7, General Reich-Moore
C                                  Resonance Regions.
C                                 *Added OUT10.
C                                 *Replaced ALL 3 way IF Statements.
C                                 *Replaced ALL LOGICAL by INTEGER.
C     VERS. 2016-1 (Jan.  2016)   *Do not Change LSSF during the
C                                  reconstrcution - for compatibility
C                                  with later URR treatment.
C                                 *Insured that all ERROR stops print
C                                  a message explaining why the code
C                                  stopped.
C                                 *Partial Energy Range Processing
C                                  no longer allowed - today's computers
C                                  are so fast that this option is now
C                                  out-of-date and no longer allowed.
C                                 *L-Value dependent fission = Earlier
C                                  was done only by entire isotope.
C                                 *Denser Starting Energy Grid.
C     VERS. 2017-1 (May   2017)   *Corrected ERROR in LRF=3 treatment.
C                                  This ERROR only existed in version
C                                  2016-1, which was never released to
C                                  the general public, so it will not
C                                  effect any results calculated by code
C                                  users.
C                                 *All floating input parameters changed
C                                  to character input + IN9 conversion.
C                                 *Added points to starting energy grid
C                                  to approximate the shape of each
C                                  resonance = based on comparisons of
C                                  0.01% to 0.1% results.
C                                 *Increased max. points to 1,200,000.
C                                 *LRF=7 Shift option no longer allowed
C                                  Set = 0, print WARNING and continue.
C                                 *Corrected COMMON/NAPRHO/NRO,NAPS
C                                  /NAPRHO/ mispelled - Freud found.
C     VERS. 2017-2 (Sept. 2017)   *Corrected Write statemnt at 5731.
C
C     OWNED, MAINTAINED AND DISTRIBUTED BY
C     ------------------------------------
C     THE NUCLEAR DATA SECTION
C     INTERNATIONAL ATOMIC ENERGY AGENCY
C     P.O. BOX 100
C     A-1400, VIENNA, AUSTRIA
C     EUROPE
C
C     ORIGINALLY WRITTEN BY
C     ------------------------------------
C     Dermott E. Cullen
C
C     PRESENT CONTACT INFORMATION
C     ---------------------------
C     Dermott E. Cullen
C     1466 Hudson Way
C     Livermore, CA 94550
C     U.S.A.
C     Telephone  925-443-1911
C     E. Mail    RedCullen1@Comcast.net
C     Website    RedCullen1.net/HOMEPAGE.NEW
C
C     Acknowledgement (Version 2004-1)
C     ==================================================================
C     The author thanks Nancy Larson, ORNL, for providing her SAMRML
c     code for comparison to RECENT output for Reich-Moore evaluations,
c     in particular to verify results for the new LFR=7 evaluations. I
c     also thank her for providing guidance to help me understand and
c     implement this new teatment for Reich-Moore parameters.
C
C     ACKNOWLEDGEMENT (VERSION 92-1)
C     ==================================================================
C     THE AUTHOR THANKS SOL PEARLSTEIN (BROOKHAVEN NATIONAL LAB) FOR
C     SIGNIFICANTLY CONTRIBUTING TOWARD IMPROVING THE ACCURACY AND
C     COMPUTER INDEPENDENCE OF THIS CODE - THANKS, SOL
C     ==================================================================
C
C     AUTHORS MESSAGE
C     ==================================================================
C     THE REPORT DESCRIBED ABOVE IS THE LATEST PUBLISHED DOCUMENTATION
C     FOR THIS PROGRAM. HOWEVER, THE COMMENTS BELOW SHOULD BE CONSIDERED
C     THE LATEST DOCUMENTATION INCLUDING ALL RECENT IMPROVEMENTS. PLEASE
C     READ ALL OF THESE COMMENTS BEFORE IMPLEMENTATION, PARTICULARLY
C     THE COMMENTS CONCERNING MACHINE DEPENDENT CODING.
C
C     AT THE PRESENT TIME WE ARE ATTEMPTING TO DEVELOP A SET OF COMPUTER
C     INDEPENDENT PROGRAMS THAT CAN EASILY BE IMPLEMENTED ON ANY ONE
C     OF A WIDE VARIETY OF COMPUTERS. IN ORDER TO ASSIST IN THIS PROJECT
C     IT WOULD BE APPECIATED IF YOU WOULD NOTIFY THE AUTHOR OF ANY
C     COMPILER DIAGNOSTICS, OPERATING PROBLEMS OR SUGGESTIONS ON HOW TO
C     IMPROVE THIS PROGRAM. HOPEFULLY, IN THIS WAY FUTURE VERSIONS OF
C     THIS PROGRAM WILL BE COMPLETELY COMPATIBLE FOR USE ON YOUR
C     COMPUTER.
C
C     PURPOSE
C     ==================================================================
C     THIS PROGRAM IS DESIGNED TO RECONSTRUCT THE RESONANCE CONTRIBUTION
C     TO THE CROSS SECTION IN LINEARLY INTERPOLABLE FORM, ADD IN ANY
C     LINEARLY INTERPOLABLE BACKGROUND CROSS SECTION AND OUTPUT THE
C     RESULT IN THE ENDF/B FORMAT. THE CROSS SECTIONS OUTPUT BY THIS
C     PROGRAM WILL BE LINEARLY INTERPOLABLE OVER THE ENTIRE ENERGY RANGE
C
C     THE RESONANCE CONTRIBUTION IS CALCULATED FOR TOTAL (MT=1),
C     ELASTIC (MT=2), CAPTURE (MT=102) AND FISSION (MT=18), ADDED
C     TO THE BACKGROUND (IF ANY) AND OUTPUT. IN ADDITION, IF THERE
C     IS A FIRST CHANCE FISSION (MT=19) BACKGROUND PRESENT THE RESONANCE
C     CONTRIBUTION OF FISSION WILL BE ADDED TO THE BACKGROUND AND
C     OUTPUT. IF THERE IS NO FIRST CHANCE FISSION (MT=19) BACKGROUND
C     PRESENT THE PROGRAM WILL NOT OUTPUT MT=19.
C
C     IN THE FOLLOWING FOR SIMPLICITY THE ENDF/B TERMINOLOGY--ENDF/B
C     TAPE--WILL BE USED. IN FACT THE ACTUAL MEDIUM MAY BE TAPE, CARDS,
C     DISK OR ANY OTHER MEDIUM.
C
C     PROCESSING DATA IN THE ENDF/B-VI FORMAT
C     ==================================================================
C     IT HAS NOW BEEN CONFIRMED (PRIVATE COMMUNICATION, CHARLES DUNFORD,
C     APRIL, 1991) THAT THE PROPER PROCEDURE TO FOLLOW WHEN THERE ARE
C     MISSING OR DUPLICATE J VALUES IS TO IN ALL CASES ADD A SEQUENCE
C     WITH NO RESONANCES TO ACCOUNT FOR THE CONTRIBUTION OF THE SEQUENCE
C     TO THE POTENTIAL SCATTERING CROSS SECTION.
C
C     THIS IS THE PROCEDURE WHICH WAS FOLLOWED BY ALL VERSIONS OF RECENT
C     SINCE 86-3 AND WILL CONTINUE TO BE THE PROCEDURE.
C
C     INPUT ENDF/B FORMAT AND CONVENTIONS
C     ==================================================================
C     ENDF/B FORMAT
C     -------------
C     THIS PROGRAM ONLY USES THE ENDF/B BCD OR LINE IMAGE FORMAT (AS
C     OPPOSED TO THE BINARY FORMAT) AND CAN HANDLE DATA IN ANY VERSION
C     OF THE ENDF/B FORMAT (I.E., ENDF/B-I, II,III, IV, V OR VI FORMAT).
C
C     IT IS ASSUMED THAT THE DATA IS CORRECTLY CODED IN THE ENDF/B
C     FORMAT AND NO ERROR CHECKING IS PERFORMED. IN PARTICULAR IT IS
C     ASSUMED THAT THE MAT, MF AND MT ON EACH LINE IS CORRECT. SEQUENCE
C     NUMBERS (COLUMNS 76-80) ARE IGNORED ON INPUT, BUT WILL BE
C     CORRECTLY OUTPUT ON ALL CARDS. THE FORMAT OF SECTION MF=1, MT=451
C     AND ALL SECTIONS OF MF=2 AND 3 MUST BE CORRECT. THE PROGRAM COPIES
C     ALL OTHER SECTION OF DATA AS HOLLERITH AND AS SUCH IS INSENSITIVE
C     TO THE CORRECTNESS OR INCORRECTNESS OF ALL OTHER SECTIONS.
C
C     ENDF/B FORMAT VERSION
C     ---------------------
C     THE FORMATS AND CONVENTIONS FOR READING AND INTERPRETING THE DATA
C     VARIES FROM ONE VERSION OF ENDF/B TO THE NEXT. HOWEVER, IF THE
C     HOLLERITH SECTION (MF=1, MT=451) IS PRESENT IT IS POSSIBLE FOR
C     THIS PROGRAM TO DISTINGUISH BETWEEN DATA IN THE ENDF/B-IV, V AND
C     VI FORMATS AND TO USE THE APPROPRIATE CONVENTIONS FOR EACH
C     ENDF/B VERSION (SEE, SUBROUTINE FILE1 FOR A DESCRIPTION OF HOW
C     THIS IS DONE). IF THE HOLLERITH SECTION IS NOT PRESENT THE
C     PROGRAM WILL ASSUME THE DATA IS IN THE ENDF/B-VI FORMAT AND USE
C     ALL CONVENTIONS APPROPRIATE TO ENDF/B-V. USERS ARE ENCOURAGED TO
C     INSURE THAT THE HOLLERITH SECTION (MF=1, MT=451) IS PRESENT IN
C     ALL EVALUATIONS.
C
C     INPUT OF ENERGIES
C     -----------------
C     ALL ENERGIES ARE READ IN DOUBLE PRECISION (BY SPECIAL FORTRAN I/O
C     ROUTINES) AND ARE TREATED IN DOUBLE PRECISION IN ALL CALCULATIONS.
C
C     OUTPUT ENDF/B FORMAT AND CONVENTIONS
C     ==================================================================
C     CONTENTS OF OUTPUT
C     ------------------
C     ENTIRE EVALUATIONS ARE OUTPUT, NOT JUST THE RECONSTRUCTED FILE
C     3 CROSS SECTIONS, E.G. ANGULAR AND ENERGY DISTRIBUTIONS ARE
C     ALSO INCLUDED.
C
C     DOCUMENTATION
C     -------------
C     THE FACT THAT THIS PROGRAM HAS OPERATED ON THE DATA IS DOCUMENTED
C     BY THE ADDITION OF COMMENT CARDS AT THE END OF EACH HOLLERITH
C     SECTION IN THE FORM
C
C     ***************** RECENT (VERSION 2017-2) ***************
C     RESONANCE CONTRIBUTION RECONSTRUCTED TO WITHIN   0.100 PER-CENT
C     COMBINED DATA NOT THINNED (ALL RESONANCE + BACKGROUND DATA KEPT)
C
C     THE ORDER OF ALL SIMILAR COMMENTS (FROM LINEAR, SIGMA1 AND GROUPY)
C     REPRESENTS A COMPLETE HISTORY OF ALL OPERATIONS PERFORMED ON
C     THE DATA, INCLUDING WHICH VERSION OF EACH PROGRAM WAS USED.
C
C     THESE COMMENT CARDS ARE ONLY ADDED TO EXISTING HOLLERITH SECTIONS,
C     I.E., THIS PROGRAM WILL NOT CREATE A HOLLERITH SECTION. THE FORMAT
C     OF THE HOLLERITH SECTION IN ENDF/B-V DIFFERS FROM THE THAT OF
C     EARLIER VERSIONS OF ENDF/B. BY READING AN EXISTING MF=1, MT=451
C     IT IS POSSIBLE FOR THIS PROGRAM TO DETERMINE WHICH VERSION OF
C     THE ENDF/B FORMAT THE DATA IS IN. WITHOUT HAVING A SECTION OF
C     MF=1, MT=451 PRESENT IT IS IMPOSSIBLE FOR THIS PROGRAM TO
C     DETERMINE WHICH VERSION OF THE ENDF/B FORMAT THE DATA IS IN, AND
C     AS SUCH IT IS IMPOSSIBLE FOR THE PROGRAM TO DETERMINE WHAT FORMAT
C     SHOULD BE USED TO CREATE A HOLLERITH SECTION.
C
C     REACTION INDEX
C     --------------
C     THIS PROGRAM DOES NOT USE THE REACTION INDEX WHICH IS GIVEN IN
C     SECTION MF=1, MT=451 OF EACH EVALUATION.
C
C     THIS PROGRAM DOES NOT UPDATE THE REACTION INDEX IN MF=1, MT=451.
C     THIS CONVENTION HAS BEEN ADOPTED BECAUSE MOST USERS DO NOT
C     REQUIRE A CORRECT REACTION INDEX FOR THEIR APPLICATIONS AND IT WAS
C     NOT CONSIDERED WORTHWHILE TO INCLUDE THE OVERHEAD OF CONSTRUCTING
C     A CORRECT REACTION INDEX IN THIS PROGRAM. HOWEVER, IF YOU REQUIRE
C     A REACTION INDEX FOR YOUR APPLICATIONS, AFTER RUNNING THIS PROGRAM
C     YOU MAY USE PROGRAM DICTIN TO CREATE A CORRECT REACTION INDEX.
C
C     OUTPUT FORMAT OF ENERGIES
C     -------------------------
C     IN THIS VERSION OF RECENT ALL FILE 3 ENERGIES WILL BE OUTPUT IN
C     F (INSTEAD OF E) FORMAT IN ORDER TO ALLOW ENERGIES TO BE WRITTEN
C     WITH UP TO 9 DIGITS OF ACCURACY. IN PREVIOUS VERSIONS THIS WAS AN
C     OUTPUT OPTION. HOWEVER USE OF THIS OPTION TO COMPARE THE RESULTS
C     OF ENERGIES WRITTEN IN THE NORMAL ENDF/B CONVENTION OF 6 DIGITS
C     TO THE 9 DIGIT OUTPUT FROM THIS PROGRAM DEMONSTRATED THAT FAILURE
C     TO USE THE 9 DIGIT OUTPUT CAN LEAD TO LARGE ERRORS IN THE DATA
C     JUST DUE TO TRANSLATION OF ENERGIES FROM THEIR INTERNAL (BINARY)
C     REPRESENTATION TO THE ENDF/B FORMAT.
C
C     ACCURACY OF ENERGY
C     ------------------
C     IN ORDER TO ALLOW ENERGIES TO BE ACCURATELY OUTPUT TO 9 DIGITS
C     ON SHORT WORD LENGTH COMPUTERS (E.G. IBM) ALL ENERGIES AND
C     ENERGY DEPENDENT TERMS ARE READ AND TREATED IN DOUBLE PRECISION.
C
C     OUTPUT OF RESONANCE PARAMETERS
C     ------------------------------
C     A SPECIAL CONVENTION HAS BEEN INTRODUCED REGARDING RESONANCE
C     PARAMETERS. IN ORDER TO ALLOW THE USER TO DOPPLER BROADEN AND/OR
C     SELF-SHIELD CROSS SECTIONS THE RESONANCE PARAMETERS ARE ALSO
C     INCLUDED IN THE OUTPUT WITH THE EVALUATION. IN ORDER TO AVOID THE
C     POSSIBILITY OF ADDING THE RESONANCE CONTRIBUTION A SECOND TIME
C     TWO CONVENTIONS HAVE BEEN ADOPTED TO INDICATE THAT THE RESONANCE
C     CONTRIBUTION HAS ALREADY BEEN ADDED TO THE FILE 3 CROSS SECTIONS,
C
C     (1) WHEN THE DATA IS PROCESSED BY THIS PROGRAM LRP (IN MF=1,
C     MT=451) IS SET EQUAL TO 2. THIS IS A CONVENTION WHICH HAS BEEN
C     ADOPTED AS A STANDARD CONVENTION IN ENDF/B-VI, BUT IS ONLY TO BE
C     USED FOR PROCESSED DATA, AS OPPOSED TO THE ORIGINAL EVALUATIONS.
C     IN EVALUATIONS WHICH CONTAIN MF=1, MT=451 LRP CAN BE USED TO
C     DETERMINE IF THE MATERIAL HAS BEEN PROCESSED.
C
C     (2) THE LRU FLAG IN EACH SECTION OF FILE 2 DATA IS CHANGED TO
C     LRU=LRU+3. FOR EXAMPLE WHEN READING AN ENDF/B EVALUATION LRU=0
C     (NO RESONANCES), =1 (RESOLVED) OR =2 (UNRESOLVED) INDICATES THAT
C     THE DATA IS IN THE ORIGINAL ENDF/B FORM. LRU=3 (NO RESONANCES),
C     =4 (RESOLVED) OR =5 (UNRESOLVED) INDICATES THAT THE RESONANCE
C     CONTRIBUTION HAS ALREADY BEEN ADDED TO THE FILE 3 DATA. THIS
C     SECOND CONVENTION HAS BEEN ADOPTED AS INSURANCE THAT THE RESONANCE
C     CONTRIBUTION WILL NOT BE ADDED TWICE, EVEN FOR EVALUATIONS WHICH
C     DO NOT CONTAIN MF=1, MT=451 (EVALUATIONS WHICH CONTAIN MF=1,
C     MT=451 ARE COVERED BY CONVENTION (1), DESCRIBED ABOVE).
C
C     UNIFORM TREATMENT OF RESONANCE FORMALISMS
C     ==================================================================
C     NORMALIZATION
C     =============
C     ALL OF THE RESONANCE FORMALISMS INCLUDE A FACTOR OF,
C
C     PI*(FRACTIONAL ABUNDANCE)/(K**2)
C
C     THIS FACTOR HAS BEEN REMOVED FROM THE CALCULATION OF EACH TYPE
C     OF RESONANCE FORMALISM AND IS APPLIED AS A FINAL NORMALIZATION
C     AFTER THE CALCULATION, ONLY ONE PLACE IN THIS PROGRAM.
C
C     FOR SIMPLICITY THIS TERM IS NOT INCLUDED IN THE FOLLOWING
C     DERIVATIONS - IN ALL CASES THE ACTUAL CROSS SECTION IS A PRODUCT
C     OF THE ABOVE FACTOR TIMES THE RESULTS PRESENTED BELOW.
C
C     SIMILARITIES
C     ============
C     FOR THE RESOLVED RESONANCE REGION, EXCEPT FOR SINGLE LEVEL BREIT
C     WIGNER, PARAMETERS ALL OF THE FORMALISMS DEFINE THE CROSS SECTIONS
C     IN AN EQUIVALENT FORM,
C
C     TOTAL    = 2*GJ*REAL(1 - U)
C              = 2*GJ*(1 - REAL(U))
C     ELASTIC  =   GJ*(1 - U)**2
C              =   GJ*((1 - 2*REAL(U)) + (REAL(U)**2 + IM(U)**2))
C              = 2*GJ*(1 - REAL(U)) - GJ*(1 - (REAL(U)**2 + IM(U)**2))
C
C     SINCE THE FIRST TERM IS THE TOTAL, THE SECOND TERM MUST BE
C     ABSORPTION. SO WE FIND,
C
C     ABSORPTION = GJ*(1 - (REAL(U)**2 + IM(U)**2))
C
C     IN ALL CASES U IS DEFINED IN THE FORM,
C
C     U        = EXP(-I*2*PS)*((1-X) - I*Y)
C
C     WHERE (X) AND (Y) ARE RELATED TO THE SYMMETRIC AND ANTI-SYMMETRIC
C     CONTRIBUTIONS OF THE RESONANCES, RESPECTIVELY. ONLY THE DEFINITION
C     OF (X) AND (Y) WILL BE DIFFERENT FOR EACH RESONANCE FORMALISM.
C     BELOW WE WILL SHOW THAT WHAT MIGHT APPEAR TO BE A STRANGE CHOICE
C     OF DEFINITION OF THE SIGN OF (X) AND(Y) HAS BEEN SELECTED SO THAT
C     FOR BREIT-WIGNER PARAMETERS (X) AND (Y) CORRESPOND EXACTLY TO THE
C     SYMMETRIC AND ANTI-SYMMETRIC CONTRIBUTION OF THE RESONANCES.
C
C     U        = (COS(2*PS) - I*SIN(2*PS))*((1-X) - I*Y)
C              =   ((1-X)*COS(2*PS) - Y*SIN(2*PS))
C              =-I*((1-X)*SIN(2*PS) + Y*COS(2*PS))
C
C     REAL(U)  = ((1-X)*COS(2*PS) - Y*SIN(2*PS))
C     IM(U)    =-((1-X)*SIN(2*PS) + Y*COS(2*PS))
C
C     R(U)**2  =((1-X)*COS(2*PS))**2 + (Y*SIN(2*PS))**2
C               -2*(1-X)*Y*COS(2*PS)*SIN(2*PS)
C     I(U)**2  =((1-X)*SIN(2*PS))**2 + (Y*COS(2*PS))**2
C               +2*(1-X)*Y*COS(2*PS)*SIN(2*PS)
C
C     THE TERMS 2*(1-X)*Y*COS(2*PS)*SIN(2*PS) CANCEL AND UPON USING
C     THE IDENTITY COS(2*PS)**2 + SIN(2*PS)**2 = 1,
C
C     SUM      = (1-X)**2 + (Y)**2
C
C     WE NOW HAVE ALL THE QUANTITIES THAT WE NEED TO DEFINE THE CROSS
C     SECTIONS,
C
C     ELASTIC
C     =======
C     ELASTIC  =GJ*(1 - 2*REAL(U) + (REAL(U)**2 + IM(U)**2))
C              =GJ*(1 - 2*((1-X)*COS(2*PS)-Y*SIN(2*PS))+(1-X)**2+(Y)**2)
C
C     THIS CAN BE WRITTEN AS A SUM OF 2 SQUARES,
C
C     ELASTIC  =GJ*(COS(2*PS) - (1-X))**2 + (SIN(2*PS) + Y)**2)
C
C              =GJ*((COS(2*PS))**2 - 2*(1-X)*COS(2*PS) + (1-X)**2) +
C                   (SIN(2*PS))**2 + 2*Y*SIN(2*PS)    + (Y)**2)
C
C     AGAIN USING THE IDENTITY COS(2*PS)**2 + SIN(2*PS)**2 = 1, WE CAN
C     SEE THAT THE DEFINITION AS THE SUM OF 2 SQUARES IS IDENTICAL TO
C     THE PRECEDING DEFINITION OF THE ELASTIC.
C
C     ELASTIC  =GJ*(COS(2*PS) - (1-X))**2 + (SIN(2*PS) + Y)**2)
C              =GJ*((COS(2*PS)-1) + X)**2 + (SIN(2*PS) + Y)**2)
C
C     USING THE IDENTITY (1 - COS(2*PS))) = 2*SIN(PS)**2, WE OBTAIN
C     THE FINAL FORM FOR THE ELASTIC,
C
C     ELASTIC =GJ*(2*SIN(PS)**2 - X)**2 + (SIN(2*PS) + Y)**2)
C
C     ABSORPTION
C     ==========
C     ABSORPTION = GJ*(1 - (REAL(U)**2 + IM(U)**2))
C                = GJ*(1 - ((1-X)**2   + (Y)**2)
C                = GJ*(1 - (1 - 2*X + (X)**2 + (Y)**2)
C                = GJ*(2*X - (X)**2 + (Y)**2)
C
C     SINCE PHYSICALLY THE ABSORPTION CANNOT BE NEGATIVE WE CAN SEE
C     THAT (X) MUST BE POSITIVE AND 2*X MUST BE GREATER THAN
C     (X)**2 + (Y)**2, FOR ALL OF THE FORMALISMS.
C
C     TOTAL
C     =====
C     IN THIS PROGRAM THE TOTAL CROSS SECTION IS ALWAYS DEFINED TO BE
C     THE SUM OF ITS PARTS - SO THE ABOVE DEFINITION IS NEVER EXPLICITLY
C     USED. HOWEVER, WE CAN LEARN SOMETHING BY EXAMINING THE DEFINITION,
C
C     TOTAL    = 2*GJ*REAL(1 - U)
C              = 2*GJ*(1 - (((1-X)*COS(2*PS) - Y*SIN(2*PS)))
C              = 2*GJ*((1 - COS(2*PS))*(1-X) - (1-X) + Y*SIN(2*PS))
C              = 2*GJ*(2*SIN(PS)**2*(1-X)    - (1-X) + Y*SIN(2*PS))
C
C              = 4*GJ*SIN(PS)**2 +
C                2*GJ*((X-1) - 2*X*SIN(PS)**2 +  Y*SIN(2*PS))
C
C     THE IMPORTANT POINT TO NOTE IS THAT THE DEFINITION OF THE TOTAL
C     DOES NOT EXPLICITLY CONTAIN ANY DEPENDENCE ON X**2 AND Y**2 -
C     THE LEVEL-LEVEL INTERFERENCE TERMS.
C
C     THIS IMPLIES THAT IF A GIVEN SET OF RESONANCE PARAMETERS ARE USED
C     WITH THIS DEFINITION THEY WILL PRODUCE EXACTLY THE SAME TOTAL
C     CROSS SECTION - WHETHER WE CLAIM THE PARAMETERS HAVE BEEN
C     PRODUCED USING A SINGLE OR MULTI-LEVEL FIT. THIS RESULT COULD
C     BE VERY MISLEADING, IF THIS RESULT FOR THE TOTAL IS IMPLIED TO
C     MEAN THAT ONE INTERPRETATION OR THE OTHER WILL NOT HAVE ANY
C     EFFECT ON THE INDIVIDUAL CROSS SECTIONS.
C
C     STARTING FROM EXACTLY THE SAME RESONANCE PARAMETERS, RELATIVE TO
C     THE RESULTS OBTAINED USING THE SINGLE LEVEL FORMULA, MULTI-LEVEL
C     RESULTS WILL TEND TO ALWAYS DECREASE THE ABSORPTION AND INCREASE
C     THE ELASTIC. THIS CAN BE IMMEDIATELY SEEN FROM OUR GENERAL
C     MULTI-LEVEL DEFINITION OF ABSORPTION,
C
C     ABSORPTION =GJ*(2*X - ((X)**2 + (Y)**2))
C
C     THE SINGLE LEVEL ABSORPTION IS,
C
C     ABSORPTION =GJ*(2*X)
C
C     THE DIFFERENCE BETWEEN THE TWO IS -2*GJ*(X**2 + Y**2), SO THAT
C     REGARDLESS OF HOW WE DEFINE (X) AND (Y) THE INCLUSION OF THIS
C     TERM WILL ALWAYS DECREASE ABSORPTION. SINCE THE TOTAL CROSS
C     SECTION IS THE SAME IN BOTH CASE, THIS MEANS THAT THE ELASTIC
C     HAS BEEN INCREASED BY THIS AMOUNT.
C
C     AGAIN, THESE RESULTS ARE BASED ON STARTING FROM EXACTLY THE SAME
C     PARAMETERS - IN ANY ACTUAL CASE THE PARAMETERS BASED ON A SINGLE
C     OR MULTI-LEVEL FIT WILL BE QUITE DIFFERENT - THE POINT THAT WE
C     WANT TO STRESS HERE IS THAT YOU SHOULD NEVER USE PARAMETERS
C     WHICH HAVE BEEN DEFINED BY A FIT USING ONE FORMALISM - IN THE
C     EQUATIONS FOR A DIFFERENT FORMALISM - AND ASSUME THAT THE RESULTS
C     WILL BE CONSISTENT - AND NEVER USE THE TOTAL CROSS SECTION TO
C     SEE WHETHER OR NOT A SET OF SINGLE LEVEL PARAMETERS CAN BE USED
C     WITH A MULTI-LEVEL FORMALISM.
C
C     POTENTIAL CROSS SECTION
C     =======================
C     FAR FROM RESONANCES (X) AND (Y) WILL BE SMALL AND THE ELASTIC
C     CROSS SECTION REDUCES TO,
C
C     ELASTIC =GJ*(2*SIN(PS)**2)**2     + (SIN(2*PS))**2
C             =GJ*4*(SIN(PS)**4         + SIN(2*PS)**2
C
C     USING THE IDENTITY SIN(2*PS) = 2*SIN(PS)*COS(PS)
C
C             =4*GJ*(SIN(PS)**4         + (SIN(PS)*COS(PS))**2)
C             =4*GJ*SIN(PS)**2*(SIN(PS)**2 + COS(PS)**2)
C             =4*GJ*SIN(PS)**2
C
C     WHICH IS THE POTENTIAL CROSS SECTION. NOTE THAT THIS RESULT IS
C     INDEPENDENT OF THE FORMALISM USED, AS IT MUST PHYSICALLY BE,
C     AND AS SUCH ALTHOUGH AS YET WE HAVE NOT DEFINED IT, WE CAN
C     NOW SEE THAT IN ALL CASES (PS) MUST BE THE PHASE SHIFT AND FOR
C     CONSISTENCY IT MUST BE DEFINED USING EXACTLY THE SAME DEFINITION
C     IN ALL CASES.
C
C     IN ADDITION SINCE PHYSICALLY FOR EACH L VALUE WE EXPECT TO OBTAIN
C     A POTENTIAL CROSS SECTION,
C
C     4*(2*L+1)*SIN(PS)**2
C
C     OBVIOUSLY FOR CONSISTENCY WE MUST HAVE,
C
C     (2*L+1) = (SUM OVER J) GJ
C
C     ONLY IN THIS CASE WILL THE RESULTS BE CONSISTENT - THIS POINT WILL
C     BE DISCUSSED IN DETAIL BELOW.
C
C     WHAT ARE THIS TERMS (X) AND (Y)
C     ===============================
C     (X) AND (Y) CAN BE EASILY IDENTIFIED BY CONSIDERING THE SINGLE
C     AND MULTI-LEVEL BREIT WIGNER FORMALISMS. IN THESE CASES WE WILL
C     FIND THAT,
C
C     X        = GAM(N)*GAM(T)/2/DEN
C     Y        = GAM(N)*(E-ER)/DEN
C     DEN      = ((E-ER)**2 + (GAM(T)/2)**2)
C
C     EXTREME CARE HAS TO BE USED TO PROPERLY DEFINE (Y) SUCH THAT IT
C     IS NEGATIVE FOR E LESS THAN ER AND POSITIVE FOR E GREATER THAN
C     ER. I WILL MERELY MENTION THAT THE EQUATIONS FOR ALL FORMALISMS
C     IN ENDF-102 DO NOT CONSISTENTLY USE (E - ER) - IN SOME CASES
C     THIS IS WRITTEN AS (ER - E), WHICH CAN LEAD TO AN INCORRECT
C     SIGN IN THE DEFINITION OF THE (Y) THAT WE REQUIRE.
C
C     THE INTERFERENCE TERMS CAN BE WRITTEN IN TERMS OF,
C     1) LEVEL-SELF INTERFERENCE  = THE CONTRIBUTION OF EACH LEVEL
C                                   INTERFERRING WITH ITSELF
C     2) LEVEL-LEVEL INTERFERENCE = THE CONTRIBUTION OF EACH LEVEL
C                                   INTERFERRRING WITH ALL OTHER LEVELS
C
C     WE WILL REFER TO THESE TWO AS (L-S) AND (L-L),
C
C     X**2     = (GAM(N)*(GAM(T)/2)**2/(DEN)**2      + (L-L)
C              = (GAM(N)**2*((GAM(T)/2)**2)/(DEN)**2 + (L-L)
C     Y**2     = (GAM(N))**2*((E-ER))**2/(DEN)**2    + (L-L)
C
C     X**2+Y**2= GAM(N)**2*DEN/(DEN)**2 = GAM(N)**2/DEN + (L-L)
C
C     TO SEE THE EFFECT OF INCLUDING MULTI-LEVEL INTERFERENCE WE CAN
C     CONSIDER OUR GENERAL EXPRESSION FOR ABSORPTION,
C
C     ABSORPTION =GJ*(2*X - ((X)**2 + (Y)**2))
C
C     AND NOTE THAT FOR BOTH SINGLE AND MULTI-LEVEL BREIT WIGNER THE
C     ENDF-102 SAYS TO TREAT ABSORPTION IN A SINGLE LEVEL APPROXIMATION
C     I.E., IGNORE LEVEL-LEVEL INTERFERENCE. IF ALL INTERFERENCE IS
C     IGNORED THIS IS EQUIVALENT TO COMPLETELY IGNORING X**2 + Y**2 AND
C     DEFINING,
C
C     ABSORPTION =GJ*2*X
C                =2*GJ*GAM(N)*GAM(T)/DEN
C
C     WHICH IS INCORRECT - SINCE THIS SEEMS TO INDICATE EVERYTHING IS
C     ABSORBED. IN ORDER TO OBTAIN THE CORRECT EXPRESSION WE CANNOT
C     COMPLETELY IGNORE INTERFERENCE - WE CAN IGNORE LEVEL-LEVEL
C     INTERFERENCE, BUT WE MUST INCLUDE LEVEL-SELF INTERFERENCE,
C
C     X**2+Y**2= GAM(N)**2/DEN
C
C     ABSORPTION =GJ*(2*X - ((X)**2 + (Y)**2))
C                =GJ*GAM(N)*(GAM(T)-GAM(N))/DEN
C                =GJ*GAM(N)*GAM(A)/DEN
C
C     SUMMARY
C     =======
C     AN IMPORTANT POINT TO NOTE IS THE DEFINITION OF (X) AND (Y)
C     WHICH IN ALL CASES WILL CORRESPOND TO THE SYMMETRIC AND
C     ANTI-SYMMETRIC CONTRIBUTION OF THE RESONANCES. IN PARTICULAR
C     DEFINING (U) IN TERMS OF (1-X) INSTEAD OF (X) IS EXTREMELY
C     IMPORTANT. NOTE, THAT THE DEFINITION OF THE ELASTIC AND
C     ABSORPTION ONLY INVOLVE (X), NOT (1-X). FAR FROM RESONANCES
C     (X) CAN BE EXTREMELY SMALL, THEREFORE (1-X) WILL BE VERY CLOSE
C     TO (1). IF THE CALCULATION PROCEEDS BY FIRST CALCULATING (1-X)
C     AND THEN DEFINING (X) BY SUBTRACTING (1), EXTREME ROUND-OFF
C     PROBLEMS CAN RESULT. THESE PROBLEMS CAN BE AVOIDED BY IN ALL
C     CASES DEFINING (X) DIRECTLY, WITHOUT ANY DIFFERENCES.
C
C     IN EACH FORMALISM THE DEFINITION OF (X) AND (Y) MAY BE DIFFERENT
C     BUT ONCE WE HAVE DEFINED (X) AND (Y) WE CAN IMMEDIATELY WRITE
C     THE CROSS SECTIONS USING A UNIFORM DEFINITION,
C
C     ELASTIC =GJ*(2*SIN(PS)**2 - X)**2 + (SIN(2*PS) + Y)**2)
C
C     ABSORPTION =-GJ*(2*X + (X)**2 + (Y)**2)
C
C     AND DEFINE THE TOTAL AS THE SUM OF THESE 2 PARTS.
C
C     RELATIONSHIP TO SINGLE LEVEL
C     ============================
C     HOW DO THE SINGLE AND MULTI-LEVEL FORMALISMS COMPARE. TO SEE,
C     STARTING FROM OUR GENERAL DEFINITION OF THE ELASTIC IN THE FORM,
C
C     ELASTIC =GJ*(2*SIN(PS)**2 + X)**2 + (SIN(2*PS) + Y)**2)
C             =GJ*(4*SIN(PS)**4 - 4*X*SIN(PS)**2 + X**2
C                + SIN(2*PS)**2 + 2*Y*SIN(2*PS)  + Y**2)
C
C             =4*GJ*SIN(PS)**2 +
C                GJ*(X**2 + Y**2
C                   -4*X*SIN(PS)**2
C                   +2*Y*SIN(2*PS))
C
C     AND OUR SPECIFIC DEFINITIONS OF (X) AND (Y) FOR MULTI-LEVEL BREIT-
C     WIGNER PARAMETERS,
C
C     X        = GAM(N)*GAM(T)/2/DEN
C     Y        = GAM(N)*(E-ER)/DEN
C     DEN      = ((E-ER)**2 + (GAM(T)/2)**2)
C
C     X**2+Y**2= GAM(N)**2/DEN + (L-L)
C
C     WE CAN RECOGNIZE X**2 AND Y**2 AS THE INTERFERENCE - (L-S) + (L-L)
C     TERMS IN THE MULTI-LEVEL FORMALISM. IN ORDER TO OBTAIN THE SINGLE
C     LEVEL EQUATION WE CAN ASSUME THAT EACH LEVEL DOES NOT INTERFERE
C     WITH ANY OTHER LEVEL - THEREFORE THE (L-L) CONTRIBUTION IS ZERO.
C
C     ELASTIC =4*GJ*SIN(PS)**2 +
C                GJ*GAM(N)*(GAM(N)
C                           -2*GAM(T)*SIN(PS)**2
C                           +2*(E-ER)*SIN(2*PS))/DEN
C
C     WHICH IS THE FORM THAT IT APPEARS IN ENDF-102, EXCEPT FOR TWO
C     TYPOGRAPHICAL ERRORS IN THE SECOND TERM,
C
C     -2*GAM(T)*SIN(PS)**2
C
C     WHICH IN ENDF-102 IS WRITTEN,
C
C     -2*(GAM(T)-GAM(N))*SIN(2*PS)**2
C
C     PROGRAM CONVENTIONS
C     ==================================================================
C     MINIMUM INPUT DATA
C     ------------------
C     FOR EACH MATERIAL TO BE PROCESSED THE MINIMUM INPUT DATA ARE THE
C     RESONANCE PARAMETERS IN FILE 2. IF THERE ARE NO FILE 2 PARAMETERS
C     IN A GIVEN MATERIAL THE ENTIRE MATERIAL WILL SIMPLY BE COPIED.
C     NEITHER THE HOLLERITH SECTION (MF=1, MT=451) NOR THE BACKGROUND
C     CROSS SECTION (SECTIONS OF MF=3) NEED BE PRESENT FOR THIS PROGRAM
C     TO EXECUTE PROPERLY. HOWEVER, SINCE THE CONVENTIONS USED IN
C     INTERPRETING THE RESONANCE PARAMETERS DEPENDS ON ENDF/B VERSION
C     USERS ARE STRONGLY RECOMMENDED TO INSURE THAT MF=1, MT=451 IS
C     PRESENT IN EACH MATERIAL TO ALLOW THE PROGRAM TO DETERMINE THE
C     ENDF/B FORMAT VERSION.
C
C     RESONANCE PARAMETERS
C     --------------------
C     RESONANCE PARAMETERS MAY BE REPRESENTED USING ANY COMBINATION
C     OF THE REPRESENTATIONS ALLOWED IN ENDF/B,
C     (1) RESOLVED DATA
C         (A) SINGLE LEVEL BREIT-WIGNER
C         (B) MULTI-LEVEL BREIT-WIGNER
C         (C) ADLER-ADLER
C         (D) REICH-MOORE
C         (E) HYBRID R-FUNCTION
C     (2) UNRESOLVED DATA
C         (A) ALL PARAMETERS ENERGY INDEPENDENT
C         (B) FISSION PARAMETERS ENERGY DEPENDENT
C         (C) ALL PARAMETERS ENERGY DEPENDENT
C
C     THE FOLLOWING RESOLVED DATA FORMALISMS ARE NOT TREATED BY THIS
C     VERSION OF THE CODE AND WILL ONLY BE IMPLEMENTED AFTER EVALUATIONS
C     USING THESE FORMALISMS ARE AVAILABLE TO THE AUTHOR OF THIS CODE
C     FOR TESTING IN ORDER TO INSURE THAT THEY CAN BE HANDLED PROPERLY
C         (A) GENERAL R-MATRIX
C
C     CALCULATED CROSS SECTIONS
C     -------------------------
C     THIS PROGRAM WILL USE THE RESONANCE PARAMETERS TO CALCULATE THE
C     TOTAL, ELASTIC, CAPTURE AND POSSIBLY FISSION CROSS SECTIONS. THE
C     COMPETITIVE WIDTH WILL BE USED IN THESE CALCULATIONS, BUT THE
C     COMPETITIVE CROSS SECTION ITSELF WILL NOT BE CALCULATED. THE
C     ENDF/B CONVENTION IS THAT ALTHOUGH A COMPETITIVE WIDTH MAY BE
C     GIVEN, THE COMPETITIVE CROSS SECTION MUST BE SEPARATELY TABULATED
C     AS A SECTION OF FILE 3 DATA.
C
C     RESOLVED REGION
C     ---------------
C     IN THE RESOLVED REGION THE RESOLVED PARAMETERS ARE USED TO
C     CALCULATE COLD (0 KELVIN), LINEARLY INTERPOLABLE, ENERGY DEPENDENT
C     CROSS SECTIONS.
C
C     SCATTERING RADIUS
C     -----------------
C     FOR SINGLE OR MULTI LEVEL BREIT-WIGNER PARAMETERS THE SCATTERING
C     RADIUS MAY BE SPECIFIED IN EITHER ENERGY INDEPENDENT (CONSTANT)
C     OR ENERGY DEPENDENT FORM (A TABLE OF ENERGY VS. RADIUS AND AN
C     ASSOCIATED INTERPOLATION LAW). IN ALL OTHER CASE ONLY AN ENERGY
C     INDEPENDENT SCATTERING RADIUS IS ALLOWED.
C
C     FOR ANY ONE MATERIAL (I.E. MAT) IF ENERGY DEPENDENT SCATTERING
C     RADII ARE GIVEN THE TOTAL NUMBER OF INTERPOLATION REGIONS AND
C     TABULATED VALUES FOR THE ENTIRE MATERIAL CANNOT EXCEED,
C     200 - INTERPOLATION REGIONS
C     500 - TABULATED VALUES
C     IF THESE LIMITS ARE EXCEEDED THE PROGRAM WILL PRINT AN ERROR
C     MESSAGE AND TERMINATE.
C
C     IF YOU REQUIRE A LARGER NUMBER OF INTERPOLATION REGION AND/OR
C     TABULATED VALUES,
C     (1) INTERPOLATION REGIONS - INCREASE THE DIMENSION OF NBTRHO AND
C     INTRHO IN COMMON/TABRHO/ THROUGHOUT THE PROGRAM AND CHANGE MAXSEC
C     IN SUBROUTINE RDAP (MAXSEC = MAXIMUM NUMBER OF INTERPOLATION
C     REGIONS).
C     (2) TABULATED VALUES - INCREASE THE DIMENSION OF ERHOTB, RHOTAB
C     AND APTAB IN COMMON/TABRHO/ THROUGHOUT THE PROGRAM AND CHANGE
C     MAXRHO IN SUBROUTINE RDAP (MAXRHO = MAXIMUM NUMBER OF TABULATED
C     VALUES).
C
C     RESOLVED REICH-MOORE AND MULTI-LEVEL BREIT-WIGNER PARAMETERS
C     ------------------------------------------------------------
C     CROSS SECTIONS FOR REICH-MOORE PARAMETERS ARE CALCULATED ACCORDING
C     TO THE EQUATION (1) - (8) OF SECTION D.1.3 OF ENDF-102. IN ORDER
C     TO CALCULATE CROSS SECTIONS FROM MULTI-LEVEL PARAMETERS IN A
C     REASONABLE AMOUNT OF TIME THIS PROGRAM EXPRESSES THE CROSS SECTION
C     IN TERMS OF A SINGLE SUM OVER RESONANCES (SEE, ENDF-102, SECTION
C     D.1.2, EQUATIONS 6-7), RATHER THAN AS A DOUBLE SUM (SEE, ENDF-102
C     SECTION D.1.2, EQUATION 1-2). IN ORDER FOR THE ENDF-102 EQUATIONS
C     TO BE CORRECT THE PARAMETERS MUST MEET THE FOLLOWING CONDITIONS,
C
C     (1) FOR EACH L STATE ALL PHYSICALLY POSSIBLE J SEQUENCES MUST BE
C         PRESENT. ONLY IN THIS CASE WILL THE CONTRIBUTIONS OF THE
C         INDIVIDUAL J SEQUENCES ADD UP TO PRODUCE THE CORRECT POTENTIAL
C         SCATTERING CONTRIBUTION FOR THE L STATE (SEE, ENDF-102,
C         SECTION D.1.2, EQUATIONS 6-7). IF ANY J SEQUENCE IS MISSING
C         THE PROGRAM WILL PRINT A WARNING AND ADD THE J SEQUENCE WITH
C         NO RESONANCE PARAMETERS IN ORDER TO ALLOW THE POTENTIAL
C         SCATTERING TO BE CALCULATED CORRECTLY (THIS IS EQUIVALENT TO
C         ASSUMING THAT THE EVALUATOR REALIZES THAT ALL J SEQUENCES MUST
C         BE AND ARE PRESENT AND THAT THE EVALUATION STATES THAT THERE
C         ARE NO RESONANCES WITH CERTAIN PHYSICALLY POSSIBLE J VALUES...
C         IN THIS CASE POTENTIAL CONTRIBUTION MUST STILL BE CONSIDERED).
C
C         EXAMPLE
C         =======
C         AN EXAMPLE OF WHERE THIS OCCURS AND IS IMPORTANT TO CONSIDER
C         IS U-238 IN ENDF/B-IV AND V LIBRARIES WHERE FOR L=1 THERE IS
C         ONLY A J=1/2 SEQUENCE. NOT INCLUDING THE J=3/2 SEQUENCE LEADS
C         TO UNDERESTIMATING THE POTENTIAL SCATTERING AND PRODUCES
C         MINIMA IN THE ELASTIC CROSS SECTION WHICH ARE AN ORDER OF
C         MAGNITUDE LOWER THAN THE CROSS SECTIONS OBTAINED BE INCLUDING
C         THE J=3/2 SEQUENCE.
C
C     (2) FOR A GIVEN TARGET SPIN AND L VALUE THERE MAY BE 2 POSSIBLE
C         MEANS OF OBTAINING THE SAME J VALUE. WHEN THIS OCCURS IN
C         ORDER TO CALCULATE THE CORRECT POTENTIAL SCATTERING CROSS
C         SECTION IT IS IMPORTANT TO INCLUDE THE EFFECT OF BOTH
C         POSSIBLE J SEQUENCES, EVEN THOUGH FROM THE ENDF/B DATA IT IS
C         NOT POSSIBLE TO DETERMINE WHICH OF THE 2 POSSIBLE SEQUENCES
C         ANY GIVEN RESONANCE BELONGS TO. IN THIS CASE THIS PROGRAM
C         TREAT ALL RESONANCES WITH THE SAME J VALUE AS BELONGING TO
C         THE SAME J SEQUENCE (TO ALLOW INTERFERENCE) AND WILL ADD AN
C         ADDITIONAL J SEQUENCE WITH NO RESONANCES IN ORDER TO ALLOW
C         THE POTENTIAL CROSS SECTION TO BE CALCULATED CORRECTLY. WHEN
C         THIS OCCURS A WARNING MESSAGE IS PRINTED, BUT BASED ON THE
C         ENDF/B DATA THERE IS NOTHING WRONG WITH THE DATA AND THERE IS
C         NOTHING THAT THE USER CAN DO TO CORRECT OR IN ANY WAY MODIFY
C         THE DATA TO ELIMINATE THE PROBLEM.
C
C         EXAMPLE
C         =======
C         FOR A TARGET SPIN =1 AND L=1 THE 2 RANGES OF PHYSICALLY
C         POSSIBLE J ARE 1/2, 3/2, 5/2 AND 1/2, 3/2. BY CHECKING THE
C         ENDF/B DATA IT IS POSSIBLE TO INSURE THAT THE 3 POSSIBLE
C         J VALUES (1/2, 3/2, 5/2) ARE PRESENT AND TO INCLUDE ALL 3
C         J SEQUENCES IN THE CALCULATIONS. HOWEVER, UNLESS ALL 5
C         POSSIBLE J SEQUENCES ARE INCLUDED THE STATISTICAL WEIGHTS
C         OF THE J SEQUENCES WILL NOT SUM UP TO 2*L+1 AND THE
C         POTENTIAL CROSS SECTION WILL BE UNDERESTIMATED. IN THIS
C         EXAMPLE THE SUM OF THE 3 J SEQUENCES 1/2, 3/2, 5/2 IS 2,
C         RATHER THAN 3 AS IT SHOULD BE FOR L=1, AND THE CONTRIBUTION
C         OF THE L=1 RESONANCES TO THE POTENTIAL SCATTERING CROSS
C         SECTION WILL ONLY BE 2/3 OF WHAT IT SHOULD BE, UNLESS THE
C         OTHER 2 J SEQUENCES (WITH DUPLICATE J VALUES) ARE INCLUDED
C         IN THE CALCULATION.
C
C     (3) EACH RESONANCE MUST HAVE AN ASSIGNED, PHYSICALLY POSSIBLE
C         J VALUE. PHYSICALLY IMPOSSIBLE OR AVERAGE J VALUES CANNOT BE
C         UNIQUELY INTERPRETED USING THE EQUATIONS IN ENDF-102 AND
C         THEIR USE WILL USUALLY RESULT IN PHYSICALLY UNRELIABLE CROSS
C         SECTIONS. THIS PROGRAM WILL CHECK ALL J VALUES AND IF ANY ARE
C         ARE FOUND TO BE PHYSICALLY IMPOSSIBLE (BASED ON TARGET SPIN
C         AND L VALUE) AN ERROR MESSAGE WILL BE PRINTED TO INDICATE THAT
C         THE RECONSTRUCTED CROSS SECTIONS WILL BE UNRELIABLE AND THE
C         PROGRAM WILL CONTINUE. IN AN ATTEMPT TO CALCULATE THE CORRECT
C         POTENTIAL SCATTERING CROSS SECTION THIS PROGRAM WILL SUBTRACT
C         THE POTENTIAL SCATTERING CONTRIBUTION DUE TO ALL FICTICIOUS J
C         SEQUENCES AND ADD THE CONTRIBUTION OF ALL PHYSICALLY POSSIBLE
C         J SEQUENCES (AS DESCRIBED ABOVE).
C
C         WARNING (LET THE USER BEWARE)
C         =============================
C         (A) IT CANNOT BE STRESSED ENOUGH THAT CROSS SECTIONS OBTAINED
C             USING PHYSICALLY IMPOSSIBLE J VALUES FOR REICH-MOORE AND
C             MULTI-LEVEL BREIT-WIGNER RESONANCE PARAMETERS WILL RESULT
C             IN UNRELIABLE CROSS SECTIONS. THE DECISION TO HAVE THIS
C             PROGRAM CONTINUE TO PROCESS WHEN THIS CONDITION IS FOUND
C             IS BASED ON AN ATTEMPT TO ALLOW THE USER TO AT LEAST HAVE
C             SOME RESULTS (HOWEVER BAD THEY MAY BE) IF THERE IS NO
C             OTHER EVALUATED DATA AVAILABLE.
C         (B) EVEN THOUGH THE REICH-MOORE AND MULTI-LEVEL EQUATIONS ARE
C             DEFINED AS ABSOLUTE OR SQUARED CONTRIBUTIONS WHICH MUST
C             ALL BE PHYSICALLY POSSIBLE, ATTEMPTING TO CORRECT THE
C             POTENTIAL CROSS SECTION (AS DESCRIBED ABOVE) CAN LEAD TO
C             NEGATIVE ELASTIC CROSS SECTIONS. THIS IS BECAUSE BASED ON
C             THE INFORMATION AVAILABLE IN THE EVALUATION IT IS NOT
C             NOT POSSIBLE TO CORRECTLY ACCOUNT FOR THE INTERFERENCE
C             BETWEEN THE RESONANCE AND POTENTIAL CONTRIBUTIONS FOR EACH
C             J SEQUENCE.
C
C     UNRESOLVED RESONANCE REGION
C     ---------------------------
C     IN THE UNRESOLVED RESONANCE REGION THE UNRESOLVED PARAMETERS
C     ARE USED TO CALCULATE INFINITELY DILUTE AVERAGE CROSS SECTIONS.
C     NOTE, IT IS IMPORTANT TO UNDERSTAND THAT FROM THE DEFINITION OF
C     THE UNRESOLVED PARAMETERS IT IS NOT POSSIBLE TO UNIQUELY CALCULATE
C     ENERGY DEPENDENT CROSS SECTIONS. ONLY AVERAGES OR DISTRIBUTIONS
C     MAY BE CALCULATED.
C
C     UNRESOLVED INTERPOLATION
C     ------------------------
C     IN THE UNRESOLVED RESONANCE REGION CROSS SECTIONS AT EACH ENERGY
C     ARE CALCULATED BY INTERPOLATING PARAMETERS. THIS IS THE CONVENTION
C     USED IN ENDF/B-IV AND EARLIER VERSIONS OF ENDF/B. THE ENDF/B-V
C     CONVENTION OF INTERPOLATING CROSS SECTIONS, NOT PARAMETERS, HAS
C     BEEN ABANDONED AS IMPRACTICAL SINCE IT CAN LEAD TO THE SITUATION
C     WHERE EXACTLY THE SAME PHYSICAL DATA CAN LEAD TO DIFFERENT RESULTS
C     DEPENDING ON WHICH OF THE THREE ENDF/B UNRESOLVED PARAMTER FORMATS
C     IS USED. FOR EXAMPLE, GIVEN A SET OF ENERGY INDEPENDENT UNRESOLVED
C     PARAMETERS IT IS POSSIBLE TO CODE THESE PARAMETERS IN EACH OF THE
C     THREE ENDF/B UNRESOLVED PARAMETER FORMATS. SINCE PHYSICALLY WE
C     ONLY HAVE ONE SET OF PARAMETERS WE WOULD EXPECT THE RESULTS TO BE
C     INDEPENDENT OF HOW THEY ARE REPRESENTED IN ENDF/B. UNFORTUNATELY
C     USING THE ENDF/B-V CONVENTION TO INTERPOLATE CROSS SECTIONS CAN
C     LEAD TO THREE COMPLETELY DIFFERENT RESULTS. IN CONTRAST USING THE
C     ENDF/B-IV AND EARLIER CONVENTION OF INTERPOLATING PARAMETERS LEADS
C     TO COMPLETELY CONSISTENT RESULTS.
C
C     INTERNAL REPRESENTATION OF UNRESOLVED PARAMETERS
C     ------------------------------------------------
C     ANY OF THE THREE POSSIBLE REPRESENTATIONS OF UNRESOLVED PARAMETERS
C     CAN BE UNIQUELY REPRESENTED IN THE ALL PARAMETERS ENERGY DEPENDENT
C     REPRESENTATIONS WITH THE APPROPRIATE (ENDF/B VERSION DEPENDENT)
C     INTERPOLATION LAW. THIS IS DONE BY THE PROGRAM WHILE READING THE
C     UNRESOLVED PARAMETERS AND ALL SUBSEQUENT CALCULATIONS NEED ONLY
C     CONSIDER THE ALL PARAMETERS ENERGY DEPENDENT REPRESENTATION.
C
C     RESONANCE RECONSTRUCTION STARTING ENERGY GRID
C     ---------------------------------------------
C     AS IN ANY ITERATIVE METHOD THE WAY TO SPEED CONVERGENCE IS TO TRY
C     TO START CLOSE TO THE ANSWER. THIS PROGRAM ATTEMPTS TO DO THIS BY
C     STARTING FROM AN ENERGY GRID WHICH IS A GOOD APPROXIMATION TO A
C     SIMPLE BREIT-WIGNER LINE SHAPE,
C
C     SIGMA(X)=1.0/(1.0+X*X)
C
C     WHERE X IS THE DISTANCE FROM THE PEAK IN HALF-WIDTHS
C
C     SUBROUTINE SUBINT HAS A BUILT-IN TABLE OF NODES WHICH ARE THE
C     HALF-WIDTH MULTIPLES TO APPROXIMATE THE SIMPLE BREIT-LINE SHAPE
C     TO WITHIN 1 PER-CENT OVER THE ENTIRE INTERVAL 0 TO 500 HALF-WIDTHS
C
C     BETWEEN ANY TWO RESOLVED RESONANCES THE STARTING GRID IS BASED ON
C     THE HALF-WIDTHS OF THE TWO RESONANCES. FROM THE LOWER ENERGY
C     RESONANCE UP TO THE MID-POINT BETWEEN THE RESONANCES (MID-POINT
C     IS DEFINED HERE AS AN EQUAL NUMBER OF HALF-WIDTHS FROM EACH
C     RESONANCE) THE HALF-WIDTH OF THE LOWER ENERGY RESONANCE IS USED.
C     FROM THE MID-POINT UP TO THE HIGHER ENERGY RESONANCE THE HALF-
C     WIDTH OF THE UPPER ENERGY RESONANCE IS USED.
C
C     WITH THIS ALOGORITHM CLOSELY SPACED RESONANCES WILL HAVE ONLY
C     A FEW STARTING NODES PER RESONANCE (E.G. U-235). WIDELY SPACED
C     RESONANCES WILL HAVE MORE NODES PER RESONANCE (E.G. U-238). FOR
C     A MIX OF S, P, D ETC. RESONANCES THIS ALOGORITHM GUARANTEES AN
C     ADEQUTE DESCRIPTION OF THE PROFILE OF EVEN EXTREMELY NARROW
C     RESONANCES (WHICH MAY IMMEDIATELY CONVERGENCE TO THE ACCURACY
C     REQUESTED, THUS MINIMIZING ITERATION).
C
C     BACKGROUND CROSS SECTIONS
C     -------------------------
C     THE PROGRAM WILL SEARCH FOR BACKGROUND CROSS SECTIONS FOR TOTAL
C     (MT=1), ELASTIC (MT=2), FISSION (MT=18), FIRST CHANCE FISSION
C     (MT=19) AND CAPTURE (MT=102).
C
C     (1) THE BACKGROUND CROSS SECTIONS (FILE 3) CAN BE PRESENT OR NOT
C         PRESENT FOR EACH REACTION.
C     (2) IF FOR A GIVEN REACTION THE BACKGROUND CROSS SECTION IS
C         PRESENT, IT WILL BE ADDED TO THE RESONANCE CONTRIBUTION AND
C         THE RESULT WILL BE OUTPUT.
C     (3) IF FOR A GIVEN REACTION THE BACKGROUND IS NOT PRESENT THE
C         PROGRAM WILL,
C         (A) IF THE INPUT TO THE PROGRAM SPECIFIES NO OUTPUT FOR
C             REACTIONS WITH NO BACKGROUND THERE WILL BE NO OUTPUT.
C         (B) IF THE INPUT TO THE PROGRAM SPECIFIES OUTPUT FOR REACTIONS
C             WITH NO BACKGROUND,
C             (I) THE RESONANCE CONTRIBUTION TO TOTAL, ELASTIC OR
C                 CAPTURE WILL BE OUTPUT.
C             (II) IF ALL FISSION RESONANCE PARAMETERS ARE ZERO THE
C                  FISSION CROSS SECTION (MT=18) WILL NOT BE OUTPUT.
C                  OTHERWISE THE RESONANCE CONTRIBUTION OF THE FISSION
C                  (MT=18) WILL BE OUTPUT.
C             (III) THERE WILL BE NO OUTPUT FOR FIRST CHANCE FISSION
C                   (MT=19).
C
C     COMBINING RESONANCES AND BACKGROUND CROSS SECTIONS
C     --------------------------------------------------
C     IN ORDER TO BE COMBINED WITH THE RESONANCE CONTRIBUTION THE
C     BACKGROUND CROSS SECTIONS MUST BE GIVEN AT 0 KELVIN TEMPERATURE
C     AND MUST BE LINEARLY INTERPOLABLE. IF THESE CONDITIONS ARE MET
C     THE RESONANCE AND BACKGROUND CONTRIBUTIONS WILL BE ADDED TOGETHER
C     AND OUTPUT. IF THESE CONDITIONS ARE NOT MET THE BACKGROUND CROSS
C     SECTION WILL BE IGNORED AND ONLY THE RESONANCE CONTRIBUTION WILL
C     BE OUTPUT. IF THE BACKGROUND HAS NOT BEEN ADDED TO THE RESONANCE
C     CONTRIBUTION AFTER THIS PROGRAM FINISHES THE USER CAN MAKE THE
C     RESONANCE AND BACKGROUND CONTRIBUTIONS COMPATIBLE BY,
C
C     (1) IF THE BACKGROUND IS NOT LINEARLY INTERPOABLE, LINEARIZE THE
C         BACKGROUND (E.G., USE PROGRAM LINEAR).
C     (2) IF THE BACKGROUND IS NOT GIVEN AT 0 KELVIN, DOPPLER BROADEN
C         THE RESONANCE (NOT BACKGROUND) CONTRIBUTION TO THE SAME
C         TEMPERATURE AS THE BACKGROUND (E.G., USE PROGRAM SIGMA1).
C
C     ONCE THE RESONANCE AND BACKGROUND CONTRIBUTIONS HAVE BEEN MADE
C     COMPATIBLE THEY CAN BE ADDED TOGETHER (E.G., USE PROGRAM MIXER).
C
C     THE RECONSTRUCTION OF THE RESONANCE CONTRIBUTION TO THE CROSS
C     SECTION CAN BE QUITE EXPENSIVE (IN TERMS OF COMPUTER TIME). SINCE
C     THE RECONSTRUCTION IS PERFORMED BEFORE THE BACKGROUND CROSS
C     SECTIONS ARE READ, THE ABOVE CONVENTIONS HAVE BEEN ADOPTED IN
C     ORDER TO AVOID LOSE OF COMPUTER TIME INVOLVED IN RECONSTRUCTING
C     THE RESONANCE CONTRIBUTION.
C
C     COMMON ENERGY GRID
C     ------------------
C     THIS PROGRAM WILL RECONSTRUCT THE RESONANCE CONTRIBUTION TO THE
C     TOTAL, ELASTIC, FISSION AND CAPTURE CROSS SECTIONS ALL ON THE
C     SAME ENERGY GRID. EACH REACTION WILL THEN BE COMBINED WITH ITS
C     BACKGROUND CROSS SECTION (IF ANY) AND OUTPUT WITHOUT ANY FURTHER
C     THINNING. IF THERE ARE NO BACKGROUND CROSS SECTIONS, OR IF THE
C     BACKGROUND CROSS SECTION FOR ALL FOUR REACTIONS ARE GIVEN ON A
C     COMMON ENERGY GRID, THE OUTPUT FROM THIS PROGRAM WILL BE ON A
C     COMMON ENERGY GRID FOR ALL FOUR REACTIONS.
C
C     THERMAL ENERGY
C     --------------
C     IF THE RESONANCE REGION SPANS THERMAL ENERGY (0.0253 EV) THIS
C     POINT IS ALWAYS INCLUDED IN THE COMMON ENERGY GRID USED FOR ALL
C     REACTIONS AND WILL ALWAYS APPEAR IN THE OUTPUT DATA.
C
C     SECTION SIZE
C     ------------
C     SINCE THIS PROGRAM USES A LOGICAL PAGING SYSTEM THERE IS NO LIMIT
C     TO THE NUMBER OF POINTS IN ANY SECTION, E.G., THE TOTAL CROSS
C     SECTION MAY BE REPRESENTED BY 200,000 DATA POINTS.
C
C     SELECTION OF DATA
C     -----------------
C     THE PROGRAM SELECTS MATERIALS TO BE PROCESSED BASED EITHER ON
C     MAT (ENDF/B MAT NO.) OR ZA. THE PROGRAM ALLOWS UP TO 100 MAT OR
C     ZA RANGES TO BE SPECIFIED. THE PROGRAM WILL ASSUME THAT THE
C     ENDF/B TAPE IS IN EITHER MAT OR ZA ORDER, WHICHEVER CRITERIA IS
C     USED TO SELECT MATERIALS, AND WILL TERMINATE WHEN A MAT OR ZA
C     IS FOUND THAT IS ABOVE THE RANGE OF ALL REQUESTS.
C
C     ALLOWABLE ERROR
C     ---------------
C     THE RECONSTRUCTION OF LINEARLY INTERPOLABLE CROSS SECTIONS FROM
C     RESONANCE PARAMETERS CANNOT BE PERFORMED EXACTLY. HOWEVER IT CAN
C     BE PERFORMED TO VIRTUALLY ANY REQUIRED ACCURACY AND MOST
C     IMPORTANTLY CAN BE PERFORMED TO A TOLERANCE THAT IS SMALL COMPARED
C     TO THE UNCERTAINTY IN THE CROSS SECTIONS THEMSELVES. AS SUCH THE
C     CONVERSION OF CROSS SECTIONS TO LINEARLY INTERPOLABLE FORM CAN BE
C     PERFORMED WITH ESSENTIALLY NO LOSS OF INFORMATION.
C
C     THE ALLOWABLE ERROR MAY BE ENERGY INDEPENDENT (CONSTANT) OR ENERGY
C     DEPENDENT. THE ALLOWABLE ERROR IS DESCRIBED BY A TABULATED
C     FUNCTION OF UP TO 20 (ENERGY,ERROR) PAIRS AND LINEAR INTERPOLATION
C     BETWEEN TABULATED POINTS. IF ONLY ONE TABULATED POINT IS GIVEN THE
C     ERROR WILL BE CONSIDERED CONSTANT OVER THE ENTIRE ENERGY RANGE.
C     WITH THIS ENERGY DEPENDENT ERROR ONE MAY OPTIMIZE THE OUTPUT FOR
C     ANY GIVEN APPLICATION BY USING A SMALL ERROR IN THE ENERGY RANGE
C     OF INTEREST AND A LESS STRINGENT ERROR IN OTHER ENERGY RANGES,
C     E.G., 0.1 PER-CENT FROM 0 UP TO THE LOW EV RANGE AND A LESS
C     STRINGENT TOLERANCE AT HIGHER ENERGIES.
C
C     DEFAULT ALLOWABLE ERROR
C     -----------------------
C     IN ORDER TO INSURE CONVERENCE OF THE RESONANCE RECONSTRUCTION THE
C     ALLOWABLE ERROR MUST BE POSITIVE. IF THE USER INPUTS AN ERROR FOR
C     RESONANCE RECONSTRUCTION THAT IS NOT POSITIVE IT WILL BE SET TO
C     THE DEFAULT VALUE (CURRENTLY 0.1 PER-CENT) AND INDICATED AS SUCH
C     IN THE OUTPUT LISTING.
C
C     INTERVAL HALVING ALGORITHM
C     -------------------------
C     THIS PROGRAM WILL START BY CALCULATING THE CROSS SECTIONS AT THE
C     ENERGIES CORRESPONDING TO THE PEAK OF EACH RESONANCE, AS WELL AS
C     A FIXED NUMBER OF HALF-WIDTHS ON EACH SIDE OF EACH RESONANCE.
C     STARTING FROM THIS BASIC GRID OF POINTS THE PROGRAM WILL CONTINUE
C     TO HALF EACH INTERVAL UNTIL THE CROSS SECTIONS FOR ALL REACTIONS
C     AT THE CENTER OF THE INTERVAL CAN BE DEFINED BY LINEAR
C     INTERPOLATION FROM THE ENDS OF THE INTERVAL TO WITHIN THE USER
C     SPECIFIED ACCURACY CRITERIA.
C
C     DISTANT RESONANCE TREATMENT
C     ---------------------------
C     THE OPTION TO TREAT DISTANT RESONANCES, WHICH WAS AVAILABLE IN
C     EARLIER VERSIONS OF THIS PROGRAM, IS NO LONGER AVAILABLE, BECAUSE
C     IT WAS FOUND TO PRODUCE UNRELIABLE RESULTS. IN THIS VERSION OF
C     THE PROGRAM ALL RESONANCES ARE TREATED EXACTLY.
C
C     PROGRAM OPERATION
C     ==================================================================
C     EDIT MODE
C     ---------
C     IT IS SUGGESTED THAT BEFORE RUNNING THIS PROGRAM TO RECONSTRUCT
C     CROSS SECTIONS FROM RESONANCE PARAMETERS (WHICH CAN BE QUITE
C     EXPENSIVE) THE USER FIRST RUN THE PROGRAM IN THE EDIT MODE (SEE,
C     DESCRIPTION OF INPUT PARAMETERS BELOW). IN THE EDIT MODE THE
C     PROGRAM WILL READ, LIST AND EXTENSIVELY CHECK THE CONSISTENCY OF
C     ALL RESONANCE PARAMETERS AND ENDF/B DEFINED RESONANCE FLAGS. THIS
C     IS A VERY INEXPENSIVE MEANS OF CHECKING ALL DATA BEFORE INVESTING
C     A LARGE AMOUNT OF MONEY IN RECONSTRUCTING CROSS SECTIONS. ANY AND
C     ALL DIGNOSTICS RECEIVED FROM THE EDIT WILL SUGGEST HOW TO CORRECT
C     THE EVALUATED DATA TO MAKE IT CONSISTENT BEFORE RECONSTRUCTING
C     CROSS SECTIONS. IN ORDER TO OBTAIN MEANINGFUL RESULTS FROM THE
C     RECONSTRUCTION ALL SUGGESTED CHANGES TO THE EVALUATION SHOULD BE
C     PERFORMED BEFORE TRYING RECONSTRUCTION (OTHERWISE THE RESULT OF
C     RECONSTRUCTION WILL NOT BE RELIABLE).
C
C     RECONSTRUCTION MODE
C     -------------------
C     FOR EACH REQUESTED MATERIAL
C     ---------------------------
C     IF SECTION MF=1, MT=451 IS PRESENT COMMENTS WILL BE ADD TO
C     DOCUMENT THAT THE MATERIAL HAS BEEN PROCESSED. MF=1, MT=451 WILL
C     ALSO BE USED TO DETERMINE THE VERSION OF THE ENDF/B FORMAT WHICH
C     WILL ALLOW THE PROGRAM TO USE THE APPROPRIATE CONVENTIONS.
C
C     ALL OF THE FILE 2 RESONANCE PARAMETERS ARE FIRST READ AND THE
C     LINEARLY INTERPOLABLE CONTRIBUTION OF THE RESONANCE PARAMETERS
C     TO THE TOTAL, ELASTIC, CAPTURE AND FISSION CROSS SECTIONS IS
C     CALCULATED SIMULTANEOUSLY USING A COMMON ENERGY GRID FOR ALL
C     FOUR REACTIONS.
C
C     AFTER THE RESONANCE CONTRIBUTION HAS BEEN RECONSTRUCTED EACH OF
C     THE FIVE REACTIONS (MT=1, 2, 18, 19, 102) IS CONSIDERED SEPARATELY
C     FOR COMBINATION WILL THE BACKGROUND CROSS SECTION, IF ANY, AS
C     DESCRIBED ABOVE.
C
C     OUTPUT WILL INCLUDE THE ENTIRE EVALUATION, INCLUDING RESONANCES
C     PARAMETERS WITH LRU MODIFIED (AS DESCRIBED ABOVE) TO INDICATE
C     THAT THE RESONANCE CONTRIBUTION HAS ALREADY BEEN ADDED TO THE
C     FILE 3 CROSS SECTIONS.
C
C     THE CYCLE OF RECONSTRUCTING THE RESONANCE CONTRIBUTION AND ADDING
C     THE BACKGROUND WILL BE REPEATED FOR EACH MATERIAL REQUESTED.
C
C-----2016/3/10 - This option is no longer allowed - today's computers
C                 are so mjuch faster that this option is no longer
C                 needed.
C     PROCESS ONLY A PORTION OF RESONANCE REGION
C     ==================================================================
C     MODERN EVALUATIONS MAY BE EXTREMELY LARGE AND IT MAY NOT BE
C     POSSIBLE TO PROCESS AN ENTIRE EVALUATION (I.E., ADD THE RESONANCE
C     CONTRIBUTION) DURING A SINGLE COMPUTER RUN.
C
C     ALSO IN THE CASE WHERE YOU ARE ONLY INTERESTED IN THE CROSS
C     SECTIONS OVER A SMALL ENERGY RANGE, YOU MAY NOT WANT TO PROCESS
C     AN ENTIRE EVALUATION, E.G., IF YOU ONLY WANT TO KNOW WHAT THE
C     CROSS SECTIONS ARE NEAR THERMAL ENERGY, 0.0253 EV.
C
C     IN ORDER TO ALLOW AN EVALUATION TO BE PROCESSED USING A NUMBER OF
C     SHORTER COMPUTER RUNS AN OPTION HAS BEEN ADDED TO THIS PROGRAM TO
C     ALLOW THE USER TO SPECIFY THE ENERGY RANGE TO BE PROCESSED.
C
C     USING THIS OPTION YOU MAY START AT THE LOWEST ENERGY (ZERO UP TO
C     SOME ENERGY) AND USE THE RESULTS OF THIS RUN AS INPUT TO THE
C     NEXT RUN, WHERE YOU CAN SPECIFY THE NEXT ENERGY RANGE. THIS
C     CYCLE CAN BE REPEATED UNTIL YOU HAVE PROCESSED THE ENTIRE
C     EVALUATION.
C
C     WARNING - THIS OPTION SHOULD BE USED WITH EXTREME CARE - THIS
C     OPTION HAS BEEN RELUCTANTLY ADDED - RELUCTANTLY BECAUSE IT CAN
C     BE EXTREMELY DANGEROUS TO USE THIS OPTION UNLESS YOU CAREFULLY
C     CHECKED WHAT YOU ARE DOING.
C
C     THE OPTION SHOULD ONLY BE USED AS FOLLOWS,
C     1) YOU MUST PROCESS USING ENERGY RANGES STARTING AT LOW ENERGY
C        AND WORKING YOUR WAY TOWARD HIGH ENERGY, E.G.,
C         0.0   TO  3.0+3
C         3.0+3 TO 10.0+3
C        10.0+3 TO 80.0+3, ETC.
C     2) FOR THE LAST ENERGY RANGE THE LOWER ENERGY LIMIT MUST BE
C        NON-ZERO (WHERE TO START) AND THE UPPER ENERGY LIMIT MUST
C        BE ZERO (NO LIMIT)
C        80.0+3 TO  0.0
C
C     IF YOU ARE ONLY INTERESTED IN THE CROSS SECTION OVER A NARROW
C     ENERGY INTERVAL AND DO NOT INTENT TO MAKE ANY OTHER USE OF THE
C     RESULTS, YOU CAN IGNORE THESE WARNINGS AND MERELY SPECIFY ANY
C     ENERGY INTERVAL OVER WHICH YOU WISH CALCULATIONS TO BE
C     PERFORMED.
C
C     NORMALLY WHEN THIS PROGRAM PROCESSES AN EVALUATION IT WILL SET
C     FLAGS IN THE EVALUATION TO PREVENT THE SAME RESONANCE
C     CONTRIBUTION FROM BEING ADDED TO THE CROSS SECTION MORE THAN
C     ONCE, SHOULD YOU USE THE OUTPUT FROM THIS PROGRAM AS INPUT TO
C     THE PROGRAM.
C
C     WHEN PROCESSING ONLY PORTIONS OF THE RESONANCE REGION THIS
C     PROGRAM CANNOT SET THESE FLAGS TO PROTECT AGAINST ADDING THE
C     RESONANCE CONTRIBUTION MORE THAN ONCE - WHICH MAKES USE OF
C     THIS OPTION EXTREMELY DANGEROUS.
C
C     ONLY YOU CAN CHECK TO MAKE SURE THAT YOU HAVE CORRECTLY
C     INCLUDED EACH ENERGY RANGE ONLY ONCE - SEE THE COMMENT LINES
C     AT THE END OF SECTION, MF=1, MT=451, FOR A COMPLETE RECORD
C     OF EACH RUN USING THIS PROGRAM. THIS SECTION WILL CONTAIN
C     LINES OF THE FORM
C
C     ***************** PROGRAM RECENT (VERSION 2017-2) *************
C     ONLY PROCESS  0.00000+ 0 TO  3.00000+ 3 EV
C     ***************** PROGRAM RECENT (VERSION 2017-2) *************
C     ONLY PROCESS  3.00000+ 3 TO  1.00000+ 4 EV
C     ***************** PROGRAM RECENT (VERSION 2017-2) *************
C     ONLY PROCESS  1.00000+ 4 TO  8.00000+ 4 EV
C     ***************** PROGRAM RECENT (VERSION 2017-2) *************
C     ONLY PROCESS  8.00000+ 4 TO  2.00000+ 7 EV
C
C     YOU SHOULD CHECK TO INSURE THAT THERE ARE NO OVERLAPPING ENERGY
C     RANGES OR MISSING ENERGY RANGES.
C
C     WHEN YOU INDICATE BY INPUT THAT YOU ARE ABOUT TO PROCESS THE
C     LAST ENERGY RANGE (SEE ABOVE, LOWER ENERGY LIMIT = NON-ZERO,
C     UPPER ENERGY LIMIT = ZERO), THIS PROGRAM WILL ASSUME THAT
C     YOU HAVE NOW COMPLETED ALL PROCESSING - AND ONLY THEN WILL
C     IT SET FLAGS IN THE EVALUATION TO PREVENT THE RESONANCE
C     CONTRIBUTION FROM BEING ADDED MORE THAN ONCE. FOR THIS REASON
C     YOU CANNOT PROCESS STARTING WITH ENERGY INTERVALS AT HIGH
C     ENERGY AND WORKING TOWARD LOW ENERGY - YOU MUST START AT LOW
C     ENERGY AND WORK TOWARD HIGH ENERGY.
C-----2016/3/10 - This option is no longer allowed - today's computers
C
C     I/O FILES
C     ==================================================================
C     INPUT FILES
C     -----------
C     UNIT  DESCRIPTION
C     ----  -----------
C       2   INPUT LINE (BCD - 80 CHARACTERS/RECORD)
C      10   ORIGINAL ENDF/B DATA (BCD - 80 CHARACTERS/RECORD)
C
C     OUTPUT FILES
C     ------------
C     UNIT  DESCRIPTION
C     ----  -----------
C       3   OUTPUT REPORT (BCD - 120 CHARACTERS/RECORD)
C      11   FINAL ENDF/B DATA (BCD - 80 CHARACTERS/RECORD)
C
C     SCRATCH FILES
C     -------------
C     UNIT  DESCRIPTION
C     ----  -----------
C      12   SCRATCH FILE FOR DATA RECONSTRUCTED FROM RESONANCE
C           PARAMETERS (BINARY - 100200 WORDS/RECORD)
C      14   SCRATCH FILE FOR COMBINED FILE 2 AND 3 DATA
C           (BINARY - 40080 WORDS/RECORD)
C
C     OPTIONAL STANDARD FILE NAMES (SEE SUBROUTINE FILEIO)
C     ==================================================================
C     UNIT  FILE NAME
C     ----  ----------
C       2   RECENT.INP
C       3   RECENT.LST
C      10   ENDFB.IN
C      11   ENDFB.OUT
C      12   (SCRATCH)
C      14   (SCRATCH)
C
C     INPUT CARDS
C     ==================================================================
C     LINE  COLS.  FORMAT  DESCRIPTION
C     ----  -----  ------  -----------
C       1    1-11    I11   RETRIEVAL CRITERIA (0=MAT, 1=ZA)
C                          THIS OPTION DEFINED WHETHER COLUMNS 1-22 OF
C                          SUBSEQUENT INPUT CARDS SHOULD BE INTERPRETED
C                          TO BE MAT OR ZA RANGES.
C           12-22   E11.4  FILE 2 MINIMUM ABSOLUTE CROSS SECTION
C                          (IF 1.0E-10 OR LESS IS INPUT THE PROGRAM
C                          WILL USE 1.0E-10)
C           23-33    I11   TREATMENT OF REACTIONS FOR WHICH BACKGROUND
C                          CROSS SECTION IS NOT GIVEN.
C                          = 0 - IGNOR (I.E. NO OUTPUT)
C                          = 1 - OUTPUT RESONANCE CONTRIBUTION.
C                          THIS OPTION IS USEFUL WITH PARTIAL EVALUATION
C                          (E.G. ENDF/B-V DOSIMETRY LIBRARY) WHERE ONLY
C                          ONE OR MORE OF THE REACTIONS ARE OF ACTUAL
C                          INTEREST.
C                          WARNING...THE USE OF THIS FIELD HAS BEEN
C                          CHANGED. THIS FIELD WAS PREVIOUSLY USED TO
C                          DEFINE THE PRECISION OF THE CALCULATION AND
C                          OUTPUT. THE FORMER DEFINITION OF THIS FIELD
C                          WAS...
C                          MINIMUM ENERGY SPACING FLAG
C                          = 0 - 6 DIGIT MINIMUM ENERGY SPACING.
C                                STANDARD 6 DIGIT E11.4 OUTPUT.
C                          = 1 - 9 DIGIT MINIMUM ENERGY SPACING.
C                                STANDARD 6 DIGIT E11.4 OUTPUT.
C                          = 2 - 9 DIGIT MINIMUM ENERGY SPACING.
C                                VARIABLE 9 DIGIT F FORMAT OUTPUT.
C                          FROM EXPERIENCE IT HAS BEEN FOUND THAT
C                          FAILURE TO SET THIS OPTION TO 2 CAN RESULT
C                          IN LARGE ERRORS IN THE FINAL DATA. THEREFORE
C                          INTERNALLY THIS OPTION IS SET TO 2.
C           34-44    I11   OPERATING MODE
C                          = 0 - CACULATE. MINIMUM OUTPUT LISTING
C                          = 1 - CACULATE. LIST ALL RESONANCE PARAMETERS
C                          = 2 - EDIT MODE. NO CALCULATION. LIST ALL
C                                RESONANCE PARAMETERS.
C                          NOTE, THE EDIT MODE (=2) IS THE SUGGESTED
C                          MODE TO FIRST TEST THE CONSISTENCY OF THE
C                          EVALUATED DATA, BEFORE RECONSTRUCTING CROSS
C                          SECTIONS (SEE, COMMENTS ABOVE).
C           45-55    I11   NEGATIVE CROSS SECTIOIN TREATMENT
C                          = 0 - O.K. - NO CHANGE
C                          = 1 - SET = 0
C           56-66    I11   MONITOR MODE SELECTOR
C                          = 0 - NORMAL OPERATION
C                          = 1 - MONITOR PROGRESS OF RECONSTRUCTION OF
C                                FILE 2 DATA AND COMBINING FILE 2 AND
C                                FILE 3 DATA. EACH TIME A PAGE OF DATA
C                                POINTS IS WRITTEN TO A SCRATCH FILE
C                                PRINT OUT THE TOTAL NUMBER OF POINTS
C                                ON SCRATCH AND THE LOWER AND UPPER
C                                ENERGY LIMITS OF THE PAGE (THIS OPTION
C                                MAY BE USED IN ORDER TO MONITOR THE
C                                EXECUTION SPEED OF LONG RUNNING JOBS).
C       2    1-72    A72   ENDF/B INPUT DATA FILENAME
C                          (STANDARD OPTION = ENDFB.IN)
C       3    1-72    A72   ENDF/B OUTPUT DATA FILENAME
C                          (STANDARD OPTION = ENDFB.OUT)
C     4-N    1-11    I11   MINIMUM MAT OR ZA (SEE COLS. 1-11, LINE 1)
C           12-22    I11   MAXIMUM MAT OR ZA (SEE COLS. 1-11, LINE 1)
C                          UP TO 100 MAT OR ZA RANGES MAY BE SPECIFIED,
C                          ONE RANGE PER LINE. THE LIST IS TERMINATED
C                          BY A BLANK LINE. IF THE THE UPPER LIMIT OF
C                          ANY REQUEST IS LESS THAN THE LOWER LIMIT THE
C                          UPPER LIMIT WILL BE SET EQUAL TO THE LOWER
C                          LIMIT. IF THE FIRST REQUEST LINE IS BLANK IT
C                          WILL TERMINATE THE REQUEST LIST AND CAUSE ALL
C                          DATA TO BE RETRIEVED (SEE EXAMPLE INPUT).
C----- 2016/3/10 - Partial Processing no longer allowed.
C                  If these fields are not blank the code will STOP
C                  with a WARNING that this is no longer allowed.
C           23-33   E11.4  LOWER ENERGY LIMIT FOR PROCESSING.
C           34-44   E11.4  UPPER ENERGY LIMIT FOR PROCESSING.
C                         *THE LOWER AND UPPER ENERGY LIMITS MUST BE
C                          ZERO, OR BLANK, UNLESS YOU WISH TO ONLY
C                          PROCESS A PORTION OF RESONANCE REGIONS.
C                         *THESE ENERGY LIMITS ARE ONLY READ FROM THE
C                          FIRST MAT/ZA REQUEST LINE
C                         *IF BOTH ARE ZERO (OR BLANK) THE ENTIRE
C                          RESONANCE REGION FOR EACH MATERIAL WILL BE
C                          PROCESSED
C                         *IF LIMITS ARE INPUT ONLY THAT PORTION OF THE
C                          RESONANCE REGION FOR EACH MATERIAL WHICH
C                          LIES BETWEEN THESE LIMITS WILL BE PROCESSED
C                         *SEE INSTRUCTIONS ABOVE BEFORE USING THIS
C                          OPTION.
C----- 2016/3/10 - Partial Processing no longer allowed.
C     VARY   1-11   E11.4  ENERGY FOR FILE 2 ERROR LAW     (  SEE   )
C           12-22   E11.4  ERROR FOR FILE 2 ERROR LAW      (COMMENTS)
C                                                          ( BELOW  )
C
C     NOTE, THIS VERSION OF THE PROGRAM DOES NOT THIN THE COMBINED FILE
C     FILE 2 + 3 DATA. AS SUCH THE ERROR LAW FOR COMBINING FILE 2 + 3
C     WHICH WAS REQUIRED IN EARLIER VERSIONS OF THIS CODE ARE NO LONGER
C     REQUIRED.
C
C     THE FILE 2 ERROR LAW MAY BE ENERGY INDEPENDENT (DEFINED BY A
C     SINGLE ERROR) OR ENERGY DEPENDENT (DEFINED BY UP TO 20 ENERGY,
C     ERROR PAIRS). FOR THE ENERGY DEPENDENT CASE LINEAR INTERPOLATION
C     WILL BE USED TO DEFINE THE ERROR AT ENERGIES BETWEEN THOSE AT
C     WHICH THE ERROR IS TABULATED. THE ERROR LAW IS TERMINATED BY A
C     BLANK LINE. IF ONLY ONE ENERGY, ERROR PAIR IS GIVEN THE LAW WILL
C     BE CONSIDERED TO BE ENERGY INDEPENDENT. IF MORE THAN ONE PAIR
C     IS GIVEN IT BE CONSIDERED TO BE ENERGY DEPENDENT (NOTE, THAT
C     FOR A CONSTANT ERROR THE ENERGY INDEPENDENT FORM WILL RUN FASTER.
C     HOWEVER, FOR SPECIFIC APPLICATIONS AN ENERGY DEPENDENT ERROR MAY
C     BY USED TO MAKE THE PROGRAM RUN CONSIDERABLE FASTER).
C
C     ALL ENERGIES MUST BE IN ASCENDING ENERGY ORDER. FOR CONVERGENCE
C     OF THE FILE 2 RECONSTRUCTION ALGORITHM ALL THE ERRORS MUST BE
C     POSITIVE. IF ERROR IS NOT POSITIVE IT WILL BE SET EQUAL TO THE
C     STANDARD OPTION (CURRENTLY 0.001, CORRRESPONDING TO 0.1 PER-CENT).
C     IF THE FIRST LINE OF THE ERROR LAW IS BLANK IT WILL TERMINATE THE
C     ERROR LAW AND THE ERROR WILL BE TREATED AS ENERGY INDEPENDENT,
C     EQUAL TO THE STANDARD OPTION (CURRENTLY, 0.1 PER-CENT). SEE,
C     EXAMPLE INPUT 4.
C
C     EXAMPLE INPUT NO. 1
C     -------------------
C     CONSIDER ALL URANIUM ISOTOPES AND TH-232. CONSIDER CROSS SECTIONS
C     WHICH ARE LARGER THAN 1.0E-8 BARNS IN ABSOLUTE VALUE. ONLY OUTPUT
C     REACTIONS FOR WHICH A BACKGROUND IS GIVEN. LIST ALL PARAMETERS AND
C     CALCULATE CROSS SECTIONS. MONITOR THE EXECUTION PROGRESS OF THE
C     PROGRAM. BETWEEN 0 AND 100 EV USE 0.1 PER-CENT ACCURACY. BETWEEN
C     100 EV AND 1 KEV VARY THE ACCURACY FROM 0.1 TO 1 PER-CENT. ABOVE
C     1 KEV USE 1 PER-CENT ACCURACY.
C
C     EXPLICITLY SPECIFY THE STANDARD FILENAMES.
C
C     THE FOLLOWING 11 INPUT CARDS ARE REQUIRED.
C
C          1 1.00000-08          0          1          0         1
C ENDFB.IN
C ENDFB.OUT
C      92000      92999
C      90232                  (UPPER LIMIT AUTOMATICALLY SET TO 90232)
C                             (END REQUEST LIST)
C 0.00000+ 0 1.00000-03
C 1.00000+02 1.00000-03
C 1.00000+03 1.00000-02
C 1.00000+09 1.00000-02
C                             (END FILE 2 ERROR LAW)
C
C     EXAMPLE INPUT NO. 2
C     -------------------
C     CONSIDER ALL URANIUM ISOTOPES AND TH-232. CONSIDER CROSS SECTIONS
C     WHICH ARE LARGER THAN 1.0E-8 BARNS IN ABSOLUTE VALUE. ONLY OUTPUT
C     REACTIONS FOR WHICH A BACKGROUND IS GIVEN. CROSS SECTIONS WILL BE
C     CALCULATED, BUT PARAMETERS WILL NOT BE LISTED. THE PROGRESS OF THE
C     PROGRAM WILL NOT BE MONITORED. USE 0.1 PER-CENT ACCURACY FOR ALL
C     ENERGIES. SINCE 0.1 PER-CENT IS THE STANDARD OPTION FOR THE ERROR
C     LAW THE FIRST ERROR LAW LINE MAY BE LEFT BLANK.
C
C     LEAVE THE DEFINITION OF THE FILENAMES BLANK - THE PROGRAM WILL
C     THEN USE THE STANDARD FILENAMES.
C
C     THE FOLLOWING 7 INPUT CARDS ARE REQUIRED.
C
C          1 1.00000-08          0          0          0         0
C
C
C      92000      92999
C      90232                  (UPPER LIMIT AUTOMATICALLY SET TO 90232)
C                             (END REQUEST LIST)
C                             (USE STANDARD OPTION FOR ERROR LAW)
C
C     EXAMPLE INPUT NO. 3
C     -------------------
C     THE SAME AS EXAMPLE INPUT NO. 2, ONLY IN THIS CASE ONLY CALCULATE
C     CROSS SECTIONS OVER THE ENERGY RANGE 0.01 TO 0.1 EV - ACROSS THE
C     THERMAL ENERGY RANGE. NOTE, THE ONLY DIFFERENCE BETWEEN THE INPUT
C     PARAMETERS IN THIS CASE AND IN EXAMPLE NO. 2, IS THAT ON THE
C     SECOND INPUT LINE WE HAVE ADDED THE ENERGY RANGE 0.01 TO 0.1 EV.
C     USE \PREPRO94\LINEAR\ENDFB.OUT AS INPUT AND ENDFB.OUT AS OUTPUT -
C     SINCE ENDFB.OUT IS THE STANDARD OUTPUT FILENAME THE NAME CAN BE
C     EITHER INCLUDED IN THE INPUT OR LEFT BLANK.
C
C     THE FOLLOWING 7 INPUT CARDS ARE REQUIRED.
C
C          1 1.00000-08          0          0          0         0
C \PREPRO94\LINEAR\ENDFB.OUT
C ENDFB.OUT
C      92000      92999 1.00000- 2 1.00000- 1
C      90232                  (UPPER LIMIT AUTOMATICALLY SET TO 90232)
C                             (END REQUEST LIST)
C                             (USE STANDARD OPTION FOR ERROR LAW)
C
C     EXAMPLE INPUT NO. 4
C     -------------------
C     RECONSTRUCT ALL DATA. OUTPUT ALL REACTIONS, REGARDING OF WHETHER
C     OR NOT THERE IS A BACKGROUND CROSS SECTION. DO NOT MONITOR THE
C     PROGRESS OF THE PROGRAM. RECONSTRUCT CROSS SECTIONS TO 1 PER-CENT
C     ACCURACY. USE \ENDFB6\LINEAR\ZA092238 AS INPUT AND
C     \ENDFB6\RECENT\ZA092238 AS OUTPUT.
C
C     THE FOLLOWING 6 INPUT CARDS ARE REQUIRED.
C
C          0 0.0                 1          0          0         0
C \ENDFB6\ZA092238
C \ENDFB6\RECENT\ZA092238
C                       (RETRIEVE ALL DATA, END REQUEST LIST)
C            1.00000- 2
C                       (END FILE 2 ERROR LAW)
C
C     EXAMPLE INPUT NO. 5
C     -------------------
C     RECONSTRUCT ALL DATA. ONLY OUTPUT REACTIONS FOR WHICH A BACKGROUND
C     CROSS SECTION IS GIVEN. DO NOT MONITOR THE PROGRESS OF THE PROGRAM
C     RECONSTRUCT CROSS SECTIONS TO 0.1 PER-CENT ACCURACY. USE ENDFB.IN
C     AS INPUT AND ENDFB.OUT AS OUTPUT.
C
C     THIS CORRESPONDS TO USING ALL OF THE STANDARD OPTONS BUILT-IN TO
C     THE PROGRAM AND ALL INPUT CARDS MAY BE BLANK.
C
C     IN THIS CASE THE FOLLOWING 5 INPUT CARDS ARE REQUIRED.
C     (ZEROES ARE INDICATED ON THE FIRST LINE, BELOW, ONLY TO INDICATE
C     WHERE THE LINE IS. THE ACTUAL INPUT LINE CAN BE COMPLETELY BLANK).
C
C          0 0.0                 0          0          0         0
C                       (USE STANDARD INPUT FILENAME = ENDFB.IN)
C                       (USE STANDARD OUTPUT FILENAME = ENDFB.OUT)
C                       (RETRIEVE ALL DATA, END REQUEST LIST)
C                       (0.1 ERROR, END FILE 2 ERROR LAW)
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
C-----08/08/2012 DEFINE CODE NAME
      CHARACTER*8 CODENAME
      COMMON/NAMECODE/CODENAME
      CHARACTER*4 LINE
CAK
      CHARACTER*72 infile,outfile
CAK
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/MAXIE/NEDSEC,NEDRES,NEDNOD
      COMMON/PAGER/NPAGE,NPAGP1,NPAGM1
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/COPC/LINE(17)
      COMMON/COPI/MFIELD(3)
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/MINNIE/EMIN,EMAX,DEMIN
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C-----08/08/2012 DEFINE CODE NAME
      CODENAME = 'RECENT  '
C-----INITIALIZE TIMER
CAK   CALL TIMER
      CALL TIMERpre
C
C     DEFINE ALL I/O UNITS.
C
CAK   CALL FILEIO
      CALL FILEIOr
C-----DEFINE THE NUMBER OF POINTS IN EACH PAGE OF DATA.
C-----01/04/07 - SWITCHED FROM MAXRES TO MAXPTX.
      NPAGE=MAXPTX
      NPAGP1=NPAGE+1
      NPAGM1=NPAGE-1
C-----DEFINE NUMBER OF SECTIONS, RESONANCES AND NODES.
      NEDSEC=0
      NEDRES=0
      NEDNOD=0
C-----DEFINE MINIMUM AND MAXIMUM ENERGIES OF INTEREST AND MINIMUM
C-----ALLOWABLE ENERGY INTERVAL FOR ANY ONE RESONANCE REGION.
      EMIN=1.0D-05
      EMAX=1.0D+10
      DEMIN=1.0D-05
C-----DEFINE COMMONLY USED CONSTANTS.
      ZERO=0.0D+00
      HALF=5.0D-01
      ONE=1.0D+00
      TWO=2.0D+00
      THREE=3.0D+00
      FOUR=4.0D+00
      EIGHT=8.0D+00
      PI = DACOS(-1.0D+00)
      PI2=TWO*PI
C-----OUTPUT PROGRAM IDENTIFICATION.
      WRITE(OUTP,130)
      WRITE(*   ,130)
C-----READ AND CHECK ALL INPUT PARAMETERS.
CAK   CALL READIN
      CALL READINr(infile,outfile)
C
C     READ, LIST AND OUTPUT ENDF/B TAPE LABEL.
C
      CALL COPYL
      WRITE(OUTP,140) LINE,MFIELD(1)
      WRITE(*   ,140) LINE,MFIELD(1)
C-----INITIALIZE FISSILE FLAG OFF AND MF=1, MT=451 SECTION PRESENT FLAG
C-----TO OFF.
   10 LFI=0
      MT451=0
C-----FIND NEXT REQUESTED MAT AND TERMINATE IF RETURNED MAT IS NOT
C-----POSITIVE.
CAK   CALL NXTMAT
      CALL NXTMATr
      IF(MATH.GT.0) GO TO 40
C-----END OF RUN. TEND LINE HAS ALREADY BEEN OUTPUT.
   20 CONTINUE
      WRITE(OUTP,160) MAXSEC,MAXRES,MAXRES,NEDSEC,NEDNOD,NEDRES
      IF(NEDSEC.GT.MAXSEC.OR.NEDNOD.GT.MAXRES.OR.NEDRES.GT.MAXRES)
     1 WRITE(OUTP,170)
      WRITE(OUTP,180)
      WRITE(*   ,160) MAXSEC,MAXRES,MAXRES,NEDSEC,NEDNOD,NEDRES
      IF(NEDSEC.GT.MAXSEC.OR.NEDNOD.GT.MAXRES.OR.NEDRES.GT.MAXRES)
     1 WRITE(*   ,170)
      WRITE(*   ,180)
CAK
      CLOSE (OTAPE)
      RETURN
CAK
      CALL ENDIT
C
C     FIND FILE 1 OR 2.
C
   30 CALL CONTI
      IF(MTH.gt.0) go to 40
      CALL CONTO               ! SEND
      IF(MATH.lt.0) go to 20
      IF(MATH.eq.0) go to 10
      go to 30
c  50 IF(MFH-1) 90,60,70       ! MF/MT=1/451 Comments
c  60 IF(MTH-451) 80,130,90
   40 IF(MFH.lt.1) go to 70
      IF(MFH.gt.1) go to 50
      IF(MTH.lt.451) go to 60
      IF(MTH.eq.451) go to 110
      go to 70
   50 IF(MFH.lt.2) go to 70
      IF(MFH.eq.2) go to 120
      go to 90
C-----COPY SECTION.
   60 CALL CONTO
      CALL COPYS
      GO TO 30
C-----COPY FILE.
   70 CALL CONTO
   80 CALL COPYF
      GO TO 30
C-----COPY MAT
   90 CALL CONTO
  100 CALL COPYM
      GO TO 10
C-----MF=1, MT=451 FOUND. ADD COMMENT CARDS AND DEFINE ENDF/B FORMAT
C-----VERSION - TO USE ALL CONVENTIONS ASSOICATED WITH CORRECT FORMAT
C-----VERSION).
  110 MT451=1
CAK   CALL FILE1(IMDONE)
      CALL FILE1r(IMDONE)
C-----COPY REMAINDER OF MAT IF NO RESONANCE PARAMETERS OR MAT HAS
C-----ALREADY BEEN PROCESSED. OTHERWISE COPY TO END OF FILE 1 (MF=1).
      IF(IMEDIT.EQ.2) GO TO 80
      IF(IMDONE.le.0) go to 100
      go to 80
C-----FILE 2 FOUND. PRINT WARNING MESSAGE IF NO MF=1, MT=451.
  120 IF(MT451.LE.0) WRITE(OUTP,150)
C-----PROCESS ALL OF FILE 2 (LIST ALL PARAMETERS AND/OR PROCESS INTO
C-----POINTWISE FORM).
      CALL CONTO
CAK   CALL FILE2
      CALL FILE2r
C-----EDIT FILE 3 OR COMBINE FILE 2 AND 3 CONTRIBUTIONS AND OUTPUT.
C-----COPY REMAINDER OF MAT.
CAK   CALL FILE3
      CALL FILE3r
C-----ENTIRE MATERIAL HAS PROCESSED. PROCEED TO NEXT MAT.
      GO TO 10
  130 FORMAT(' Calculate Cross Sections from Resonance Parameters',
     1 ' (RECENT 2017-2)'/1X,78('='))
  140 FORMAT(1X,78('=')/' ENDF/B Tape Label'/1X,78('=')/1X,16A4,A2,I4)
  150 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 ' No Section MF=1, MT=451.'/
     2 ' (1) Cannot Determine ENDF/B Format Version. Will'/
     3 '     Assume and Use ALL ENDF/B-VI Conventions.'/
     4 ' (2) Cannot Determine if there are Resonance Parameters(LRP).'/
     6 ' (3) Cannot Determine if Material is Fissile (LFI).'/
     7 ' Will Read Data to Answer Points (2) AND (3).')
  160 FORMAT(
     1 ' Core Allocation and Requirements'/1X,78('=')/
     2 9X,'Sections   Nodes Parameter'/28X,'Storage'/1X,78('=')/
     3 ' Allocated',I7,I8,I10/
     4 ' Required ',I7,I8,I10/1X,78('='))
  170 FORMAT(1X,7('WARNING...'),'WARNING'/
     1       ' WARNING...Before Using this Program to Calculate'/
     2       '           Cross Sections the Core Allocation MUST'/
     3       '           be Increased to at Least the Requirements'/
     4       '           Described Above. Failure to do this Will'/
     5       '           Cause the Program to Abort During Eexecution'/
     6 1X,78('='))
  180 FORMAT(' End of Run'/1X,78('='))
      END
C     SUBROUTINE READIN
      SUBROUTINE READINr(infile,outfile)
C=======================================================================
C
C     READ AND CHECK ALL INPUT PARAMETERS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*4 MESS1,MESS2,MESS3,MESS4,MESS5,MESS6
CAK
      CHARACTER*72 infile,outfile
CAK
      CHARACTER*72 NAMEIN,NAMEOUT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/IOSTATUS/ISTAT1,ISTAT2
      COMMON/MATZA/MODGET,NMATZA,MATMIN(101),MATMAX(101)
      COMMON/OKERR2/ERRXC2,ERXC20,ENER2(1000),ER2(1000),MAXER2,NOERR
      COMMON/FIXPOT/MYPOT
      COMMON/MINSIG/SIGMIN
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
      COMMON/NAMEX/NAMEIN,NAMEOUT
      DIMENSION MESS1(2),MESS2(9,2),MESS3(10,3),MESS4(9,2),
     1 MESS5(5,2),MESS6(2)
      DATA EPMAX/2.0D+07/
C-----DEFINE STANDARD MINIMUM ALLOWABLE ERROR (PRESENTLY 0.1 PER-CENT).
      DATA ERRMIN/1.0D-03/
C-----DEFINE MINIMUM ABSOLUTE CROSS SECTION OF INTEREST (CROSS SECTION
C-----MAY BE POSITIVE OR NEGATIVE, BUT IF BOTH ENDS OF ANY INTERVAL ARE
C-----CLOSER TO ZERO THAN SIGLOW THE INTERVAL WILL NOT BE FURTHER
C-----SUB-DIVIDED).
      DATA SIGLOW/1.0D-10/
C-----DEFINE ALL OUTPUT MESSAGES.
      DATA MESS1/' MAT','  ZA'/
      DATA MESS2/
     1 '  No',' Out','put ','    ','    ','    ','    ','    ','    ',
     2 '    ',' Out','put ','(Res','onan','ce C','ontr','ibut','ion)'/
      DATA MESS3/
     1 ' Cal','cula','te. ','Mini','mum ','Outp','ut L','isti','ng. ',
     2 '    ',
     3 ' Cal','cula','te. ','List',' Res','onan','ce P','aram','eter',
     4 's.  ',
     5 ' Edi','t Mo','de. ','List',' Res','onan','ce P','aram','eter',
     6 's.  '/
      DATA MESS4/
     1 '  No',' Cha','nge ','(All','ow N','egat','ive ','Outp','ut) ',
     2 '  Ma','ke =',' 0 (','No N','egat','ive ','Outp','ut) ','    '/
      DATA MESS5/'(Sta','ndar','d Op','tion',')   ',
     1 '    ','    ','    ','    ','    '/
      DATA MESS6/' Off','  On'/
C-----INITIALIZE TO ITERATE
      NOERR = 0
C-----READ AND PRINT INTERPRETATION OF FIRST LINE OF INPUT PARAMETERS.
CAK   IF(ISTAT1.EQ.1) GO TO 20
C-----2017/5/6 - Changed all floating point to character.
CAK   READ(INP,10,END=20,ERR=20)
CAK  1 MODGET,(FIELD(j,1),j=1,11),IMBACK,IMEDIT,MAKEPLUS,MONITR
CAK10 FORMAT(I11,11A1,2I11,I11,I11)
CAK   CALL IN9(SIGMIN,FIELD(1,1))
C-----2017/5/6 - Changed all floating point to character.
CAK   GO TO 30
C-----DEFINE DEFAULT VALUES
CAK20 ISTAT1   = 1
      MODGET   = 0
CAK   SIGMIN   = 0.0
      SIGMIN   =1.e-10
      IMBACK   = 1
      IMEDIT   = 1
      MAKEPLUS = 1
      MONITR   = 1
   30 IF(MODGET.NE.0) MODGET=1
      IF(IMBACK.NE.0) IMBACK=1
      IF(MAKEPLUS.NE.0) MAKEPLUS=1
C-----INITIALIZE TO ADD MISSING OR DUPLICATE (L,S,J) SEQUENCES.
      MYPOT=1
      IF(IMEDIT.LT.10) GO TO 40
C-----AS REQUESTED, TURN ON OVERRIDE NOT TO ADD MISSING OR DUPLICATE
C-----(L,S,J) SEQUENCES.
      MYPOT=0
      IMEDIT=IMEDIT-10*(IMEDIT/10)
   40 IF(IMEDIT.LT.0.OR.IMEDIT.GT.2) IMEDIT=2
C-----IF IN EDIT MODE SET OTAPE = 0 TO INDICATE NO ENDF/B OUTPUT.
      IF(IMEDIT.EQ.2) OTAPE=0
      IF(MONITR.LE.0) MONITR=0
      IF(MONITR.GT.0) MONITR=1
      MIN=2
      IF(SIGMIN.GT.0.0D+0) GO TO 50
      SIGMIN=SIGLOW
      MIN=1
   50 CALL OUT9(SIGMIN,FIELD(1,1))
      WRITE(OUTP,360) MESS1(MODGET+1),(FIELD(M,1),M=1,11),
     1 (MESS5(J,MIN),J=1,5),(MESS2(J,IMBACK+1),J=1,9)
      WRITE(OUTP,370) (MESS3(J,IMEDIT+1),J=1,10),
     1                (MESS4(J,MAKEPLUS+1),J=1,9),
     2 MESS6(MONITR+1)
      WRITE(*   ,360) MESS1(MODGET+1),(FIELD(M,1),M=1,11),
     1 (MESS5(J,MIN),J=1,5),(MESS2(J,IMBACK+1),J=1,9)
      WRITE(*   ,370) (MESS3(J,IMEDIT+1),J=1,10),
     1                (MESS4(J,MAKEPLUS+1),J=1,9),
     2 MESS6(MONITR+1)
C-----PRINT WARNING MESSAGE IF OVERRIDE NOT TO ADD MISSING OR DUPLICATE
C-----(L,S,J) SEQUENCES HAS BEEN TURNED ON.
      IF(MYPOT.EQ.0) WRITE(OUTP,480)
C
C     READ FILENAMES - IF BLANK USE STANDARD FILENAMES
C
C-----INPUT DATA.
CAK   IF(ISTAT1.EQ.1) GO TO 70
CAK   READ(INP,60,END=70,ERR=70) NAMEIN
      NAMEIN=infile
CAK60 FORMAT(A72)
      IF(NAMEIN.EQ.' ') NAMEIN = 'ENDFB.IN'
C-----OUTPUT DATA.
CAK   READ(INP,60,END=80,ERR=80) NAMEOUT
      NAMEOUT=outfile
      IF(NAMEOUT.EQ.' ') NAMEOUT = 'ENDFB.OUT'
      GO TO 90
C-----USE DEFAULT FILENAMES
CAK70 NAMEIN  = 'ENDFB.IN'
   80 NAMEOUT = 'ENDFB.OUT'
CAK   ISTAT1 = 1
C-----PRINT FINAL FILENAMES
   90 WRITE(OUTP,100) NAMEIN,NAMEOUT
      WRITE(*   ,100) NAMEIN,NAMEOUT
  100 FORMAT(1X,78('=')/
     1 ' ENDF/B Input and Output Data Filenames'/1X,A72/
     2 1X,A72)
C
C     OPEN ENDF/B DATA FILES
C
CAK   CALL FILIO2
      CALL FILIO2r
C***** DEBUG - ACTIVATE FOR UNRESOLVED COMPETITION LISTING
C     OPEN(22,FILE='RECENT.COMPETE')
C***** DEBUG - ACTIVATE FOR UNRESOLVED COMPETITION LISTING
C
C     TERMINATE IF ERROR OPENING ENDF/B DATA FILE
C
      IF(ISTAT2.EQ.1) THEN
      WRITE(OUTP,110) NAMEIN
      WRITE(   *,110) NAMEIN
  110 FORMAT(//' ERROR - open ENDF/B data file'/1X,A72//)
      CALL ENDERROR
      ENDIF
C-----READ SELECTION RANGES (EITHER MAT OR ZA). IF MAXIMUM IS LESS
C-----THAN MINIMUM SET IT EQUAL TO MINIMUM.
      IF(MODGET.EQ.0) WRITE(OUTP,380)
      IF(MODGET.EQ.1) WRITE(OUTP,390)
      IF(MODGET.EQ.0) WRITE(*   ,380)
      IF(MODGET.EQ.1) WRITE(*   ,390)
CAK   IF(ISTAT1.EQ.1) GO TO 130
C-----2017/5/6 - Changed all floating point to character.
CAK   READ(INP,120) MATMIN(1),MATMAX(1),((FIELD(KK,KKK),KK=1,11),
CAK  1 KKK=1,2)
CAK120 FORMAT(2I11,22A1)
C-----2017/5/6 - Changed all floating point to character.
C     GO TO 140
C-----DEFINE DEFAULT VALUES
CAK Initialize MATMIN and MATMAX
      MATMIN(1)=1
      MATMAX(1)=9999
      DO NMATZA=2,101
        MATMIN(NMATZA)=0
        MATMAX(NMATZA)=0
      ENDDO
CAK
CAK130 ISTAT1    = 1
CAK   MATMIN(1) = 0
CAK   MATMAX(1) = 0
      GO TO 160
C-----------------------------------------------------------------------
C
C     2016/3/10 - Partial Range Processing no longer allowed.
C
C-----------------------------------------------------------------------
C-----CONVERT ENERGY RANGE TO PROCESS FROM CHARACTERS TO FLOATING POINT.
CAK140 CALL IN9(EPART1,FIELD(1,1))
CAK   CALL IN9(EPART2,FIELD(1,2))
C-----DEFINE WHETHER ALL OR A PORTION OF THE RESONANCE REGION WILL BE
C-----PROCESSED.
      NPART=0
      IF(EPART1.GT.0.0) NPART=1
      IF(EPART2.GT.0.0) NPART=2
      IF(EPART2.LE.0.0) EPART2=EPMAX
      IF(NPART.LE.0) GO TO 160
      CALL OUT9(EPART1,FIELD(1,1))
      CALL OUT9(EPART2,FIELD(1,2))
      WRITE(OUTP,150) ((FIELD(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,150) ((FIELD(M,KK),M=1,11),KK=1,2)
  150 FORMAT(///' ERROR - You Have Tried to Turned on the Override to',
     1                                              ' ONLY Process a'/
     2          '         Portion of the Resonance Region ',
     3                                     11A1,' to ',11A1,' eV.'/
     2          '         This option is no longer allowed.'/
     3          '         Execution Terminated.'///)
      CALL ENDERROR
C-----------------------------------------------------------------------
C
C     Requested MAT/ZA Ranges.
C
C-----------------------------------------------------------------------
C-----IF NO MAT/ZA RANGE USE STANDARD OPTION = ALL.
  160 IF(MATMIN(1).GT.0.OR.MATMAX(1).GT.0) GO TO 170
      MATMAX(1)=9999
      MODGET=0
      NMATZA=2
C-----Default = ALL
      WRITE(OUTP,410) MATMIN(1),MATMAX(1)
      WRITE(*   ,410) MATMIN(1),MATMAX(1)
      go to 180
C-----Define Range
  170 IF(MATMAX(1).LT.MATMIN(1)) MATMAX(1)=MATMIN(1)
      WRITE(OUTP,400) MATMIN(1),MATMAX(1)
      WRITE(*   ,400) MATMIN(1),MATMAX(1)
C-----PROCESS REMAINING DATA REQUESTS.
  180 DO 190 NMATZA=2,101
CAK   IF(ISTAT1.EQ.1) GO TO 210
CAK   READ(INP,120,END=200,ERR=200) MATMIN(NMATZA),MATMAX(NMATZA)
c-----Check input and define defaults.
      IF(MATMIN(NMATZA).LE.0.AND.MATMAX(NMATZA).LE.0) GO TO 210
      IF(MATMAX(NMATZA).LT.MATMIN(NMATZA)) MATMAX(NMATZA)=MATMIN(NMATZA)
      WRITE(OUTP,400) MATMIN(NMATZA),MATMAX(NMATZA)
      WRITE(*   ,400) MATMIN(NMATZA),MATMAX(NMATZA)
  190 CONTINUE
      GO TO 310
C-----------------------------------------------------------------------
C
C     READ AND PRINT FILE 2 ERROR LAW. ERROR MUST BE POSITIVE AND
C     ENERGIES MUST BE IN ASCENDING ORDER. IF ERROR IS ZERO SET TO
C     STANDARD OPTION. ERROR LAW IS TERMINATED BY BLANK LINE. IF
C     FIRST LINE IS BLANK TERMINATE ERROR LAW AND SET ERROR TO
C     ENERGY INDEPENDENT STANDARD OPTION (CURRENTLY 0.1 PER-CENT).
C
C-----------------------------------------------------------------------
CAK200 ISTAT1 = 1
  210 NMATZA=NMATZA-1
CAK   IF(ISTAT1.EQ.1) GO TO 230
C-----2017/5/6 - Changed all floating point to character.
CAK   READ(INP,220,END=230,ERR=230) ((FIELD(j,k),j=1,11),k=1,2)
CAK220 FORMAT(22A1)
CAK   CALL IN9(ENER2(1),FIELD(1,1))
CAK   CALL IN9(ER2  (1),FIELD(1,2))
C-----2017/5/6 - Changed all floating point to character.
      GO TO 240
CAK230 ISTAT1   = 1
      ENER2(1) = 0.0
      ER2(1)   = 0.0
  240 IF(ENER2(1).LT.0.0) ENER2(1)=0.0
      CALL OUT9(ENER2(1),FIELD(1,1))
      IF(ER2(1).GT.0.0) GO TO 250
      ER2(1)=ERRMIN
      PERCNT=100.0*ERRMIN
      CALL OUT9(ER2(1),FIELD(1,2))
      WRITE(OUTP,430) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,430) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
      IF(ENER2(1).GT.0.0) GO TO 260
      MAXER2=2
      GO TO 300
  250 PERCNT=100.0*ER2(1)
      CALL OUT9(ER2(1),FIELD(1,2))
      WRITE(OUTP,420) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,420) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
  260 DO 290 MAXER2=2,1000
CAK   IF(ISTAT1.EQ.1) GO TO 300
C-----2017/5/6 - Changed all floating point to character.
CAK   READ(INP,220,END=230,ERR=230) ((FIELD(j,k),j=1,11),k=1,2)
CAK   CALL IN9(ENER2(MAXER2),FIELD(1,1))
CAK   CALL IN9(ER2  (MAXER2),FIELD(1,2))
C-----2017/5/6 - Changed all floating point to character.
      IF(ENER2(MAXER2).LE.0.0.AND.ER2(MAXER2).LE.0.0) GO TO 300
      CALL OUT9(ENER2(MAXER2),FIELD(1,1))
      IF(ER2(MAXER2).LE.0.0) GO TO 270
      PERCNT=100.0*ER2(MAXER2)
      CALL OUT9(ER2(MAXER2),FIELD(1,2))
      WRITE(OUTP,440) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,440) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
      GO TO 280
  270 ER2(MAXER2)=ERRMIN
      PERCNT=100.0*ERRMIN
      CALL OUT9(ER2(MAXER2),FIELD(1,2))
      WRITE(OUTP,450) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,450) ((FIELD(M,I),M=1,11),I=1,2),PERCNT
  280 IF(ENER2(MAXER2).LT.ENER2(MAXER2-1)) GO TO 350
  290 CONTINUE
      GO TO 340
C-----INITIALIZE FILE 2 RECONSTRUCTION LAW INDICES AND ALLOWABLE ERROR.
  300 MAXER2=MAXER2-1
      ERRXC2=ER2(1)
C-----09/21/02 - DECREASED FROM 0.10 TO 0.02
C     ERXC20=0.10*ERRXC2
      ERXC20=0.02*ERRXC2
C-----02/14/04 - CHECK FOR NO ITERATION = LARGE ALLOWABLE ERROR
      DO I=1,MAXER2
      IF(ER2(I).LT.1.0) RETURN
      ENDDO
      NOERR = 1 ! SET FOR NO ITERATION
      RETURN
C
C     ERROR MESSAGE SECTION. PRINT ERROR MESSAGE AND TERMINATE.
C
C-----OVER 100 MAT OR ZA RANGES.
  310 WRITE(OUTP,320)
      WRITE(*   ,320)
  320 FORMAT(///' ERROR - Over 100 Ranges----Execution Terminated'///)
  330 CALL ENDERROR
C-----OVER 1000 ENTRIES IN ERROR LAW.
  340 WRITE(OUTP,460)
      WRITE(*   ,460)
      GO TO 330
C-----ERROR LAW ENERGIES NOT IN ASCENDING ORDER.
  350 WRITE(OUTP,470)
      WRITE(*   ,470)
      GO TO 330
  360 FORMAT(
     1 ' Retrieval Criteria-----------',7X,A4/
     2 ' File 2 Mimimum Cross Section-',11A1,1X,5A4/
     3 ' Reactions with No Background-',9A4)
  370 FORMAT(
     1 ' Calculate/Edit Mode----------',1X,10A4/
     2 ' Negative Cross Sections------',1X,9A4/
     3 ' Monitor Mode-----------------',7X,A4)
  380 FORMAT(1X,78('=')/' Requested MAT Ranges'/1X,78('=')/
     1 4X,'Mimimum',4X,'Maximum'/1X,78('='))
  390 FORMAT(1X,78('=')/' Requested ZA Ranges'/1X,78('=')/
     1 4X,'Mimimum',4X,'Maximum'/1X,78('='))
  400 FORMAT(2I11)
  410 FORMAT(2I11,' (Default Option)')
  420 FORMAT(1X,78('=')/' Allowable Uncertainty'/1X,78('=')/
     1 6X,'Energy',' Uncertainty',3X,'per-cent'/1X,78('=')/
     2 1X,11A1,1X,11A1,F11.3)
  430 FORMAT(1X,78('=')/' File 2 Reconstruction Error'/1X,78('=')/
     1 6X,'Energy',6X,'Error',3X,'per-cent'/1X,78('=')/
     3 1X,11A1,1X,11A1,F11.3,' (Default Option)')
  440 FORMAT(1X,11A1,1X,11A1,F11.3)
  450 FORMAT(1X,11A1,1X,11A1,F11.3,' (Default Option)')
  460 FORMAT(' Over 1000 Ranges--Execution Terminated')
  470 FORMAT(' Energies MUST be in Ascending Order----',
     1 'Execution Terminated')
  480 FORMAT(1X,78('=')/
     1 ' WARNING - You Have Turned on the Override NOT to Add Missing'/
     2 ' or Duplicate (L,S,J) Sequences to Account for their Potential'/
     3 ' Contribution. This is NOT Recommended, Based on the Decision'/
     4 ' of the National Nuclear Data Center, Brookhaven National Lab,'/
     5 ' Private Communication, Charles Dunford, (April 1991)'/
     6  1X,78('='))
      END
CAK   SUBROUTINE NXTMAT
      SUBROUTINE NXTMATr
C=======================================================================
C
C     FIND NEXT REQUESTED MATERIAL BASED EITHER ON ZA OR MAT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 ZABCD
      CHARACTER*4 FMT5,FMTHOL
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/MATZA/MODGET,NMATZA,MATMIN(101),MATMAX(101)
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/WHATZA/IZANOW,MATNOW,TEMP3,IVERSE,INT45
      COMMON/HOLFMT/FMTHOL,ZABCD(10)
      DIMENSION IZAMIN(101),IZAMAX(101)
      EQUIVALENCE (MATMIN(1),IZAMIN(1)),(MATMAX(1),IZAMAX(1))
      DATA FMT5/' V'/
C-----READ NEXT LINE AND CHECK FOR END OF ENDF/B TAPE.
   10 CALL CONTI
      IF(MTH.gt.0) go to 20
      IF(MATH.lt.0) go to 60
      go to 10
C-----DEFINE FIXED POINT ZA.
   20 IZANOW=C1H
C-----COMPARE MAT OR ZA TO SELECTION CRITERIA.
      IMHIGH=0
      DO 50 IMATZA=1,NMATZA
      IF(MODGET.NE.0) GO TO 30
      IF(MATH.lt.MATMIN(IMATZA)) go to 40
      IF(MATH.eq.MATMIN(IMATZA)) go to 70
      IF(MATH.le.MATMAX(IMATZA)) go to 70
      go to 50
   30 IF(IZANOW.lt.IZAMIN(IMATZA)) go to 40
      IF(IZANOW.eq.IZAMIN(IMATZA)) go to 70
      IF(IZANOW.le.IZAMAX(IMATZA)) go to 70
   40 IMHIGH=1
   50 CONTINUE
C-----THIS MATERIAL HAS NOT BEEN REQUESTED. IF BEYOND RANGE OF ALL
C-----REQUESTS RUN IF COMPLETED. IF NOT SKIP TO NEXT MATERIAL.
      IF(IMHIGH.LE.0) GO TO 60
C-----SKIP TO MATERIAL END (MEND) LINE.
      CALL SKIPM
      GO TO 10
C-----END OF RUN. RETURN NEGATIVE MATH AS INDICATOR. OUTPUT TAPE END
C-----(TEND) RECORD.
   60 MATH=-1
      MFH=0
      MTH=0
      CALL OUTT
      WRITE(OUTP,90)
      WRITE(*   ,90)
      RETURN
C-----THIS MATERIAL REQUESTED. INITIALIZE OUTPUT SEQUENCE NUMBER,
C-----ENDF/B FORMAT VERSION NUMBER, ASSUME ENDF/B-VI FORMAT AND
C-----INITIALIZE FILE 3 TEMPERATURE TO ZERO (IF MF=1, MT-451 IS
C-----PRESENT THE ENDF/B VERSION NUMBER WILL BE RE-DEFINED BASED
C-----ON THE FORMAT OF MF=1, MT=451).
   70 NOSEQ=1
      MATNOW=MATH
      FMTHOL=FMT5
      IVERSE=6
      TEMP3=0.0
      INT45=2
C-----DEFINE BCD EQUIVALENT OF ZA.
      CALL ZAHOL(IZANOW,ZABCD)
C-----IDENTIFY MAT BEING PROCESSED.
      WRITE(OUTP,80) ZABCD,MATNOW
      WRITE(*   ,80) ZABCD,MATNOW
      RETURN
   80 FORMAT(1X,78('*')/' Processing ',10A1,' MAT=',I5/
     1 1X,78('*'))
   90 FORMAT(1X,78('*')/' End of ENDF/B Input Data'/1X,78('*'))
      END
CAK   SUBROUTINE FILE1(IMDONE)
      SUBROUTINE FILE1r(IMDONE)
C=======================================================================
C
C     ADD COMMENTS AT THE END OF FILE 1, SECTION 451 TO INDICATE
C     THAT THIS MATERIAL HAS BEEN PROCESSED BY PROGRAM RECENT AND
C     TO SPECIFY THE MAXIMUM ALLOWABLE ERROR.
C
C     DEFINE FORMAT TO BE ENDF/B-IV, V OR VI.
C
C     THE ENDF/B FORMAT CAN BE DETERMINED FROM THE SECOND LINE.
C     ENDF/B-IV = N1 > 0, N2 = 0,LINE COUNT (POSITIVE)
C     ENDF/B-V  = N1 = N2 =0
C     ENDF/B-VI =      N2 = VERSION NUMBER (6 OR MORE)
C     ENDF/B-VII= SAME AS VI, BUT FORMAT FIELD SAYS 7
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 ZABCD,FIELD,PROGDOC1
      CHARACTER*4 FMTHOL,FMTTAB
      CHARACTER*66 PROGDOC
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/MINSIG/SIGMIN
      COMMON/OKERR2/ERRXC2,ERXC20,ENER2(1000),ER2(1000),MAXER2,NOERR
      COMMON/WHATZA/IZANOW,MATNOW,TEMP3,IVERSE,INT45
      COMMON/HOLFMT/FMTHOL,ZABCD(10)
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
      DIMENSION FMTTAB(4),PROGDOC(9),PROGDOC1(66,9)
      EQUIVALENCE (PROGDOC(1),PROGDOC1(1,1))
      DATA FMTTAB/'IV',' V','VI','VII'/
C-----DOCUMENTATION TO ADD TO ENDF/B OUTPUT - EACH LINE IS 66
C-----CHARACTERS LONG - FIELDS 12345678901 ARE FILLED IN WITH
C-----11 CHARACTERS DURING EXECUTION.
C               1         2         3         4         5         6
C       12345678901234567890123456789012345678901234567890123456789012
C       3456
      DATA PROGDOC/
     1 ' ***************** Program RECENT (VERSION 2017-2) ***********',
     2 ' Only Process12345678901 to12345678901 eV                     ',
     3 ' for All Data Greater than12345678901 barns in Absolute Value ',
     4 ' Data Linearized to within an Accuracy of12345678901 per-cent ',
     5 ' Data Linearized Using Energy Dependent Uncertainty           ',
     6 '      Energy    Accuracy                                      ',
     7 '        (eV)  (per-cent)                                      ',
     8 ' ----------- -----------                                      ',
     9 ' 12345678901 12345678901                                      '/
C-----FILL IN REMAINDER OF FIRST LINE
      PROGDOC1(63,1) = '*'
      PROGDOC1(64,1) = '*'
      PROGDOC1(65,1) = '*'
      PROGDOC1(66,1) = '*'
C-----INITIALIZE MAT ALREADY PROCESSED FLAG TO OFF.
      IMDONE=1
C-----SAVE FLAGS WHICH INDICATES WHETHER OR NOT MAT IS FISSIONABLE AND
C-----IF RESONANCE PARAMETERS ARE PRESENT.
      LRP=L1H
      LFI=L2H
C-----HEAD LINE OF SECTION HAS BEEN READ. READ SECOND LINE AND TEST
C-----FOR ENDF/B-IV FORMAT.
      CALL CARDI(C1A,C2A,L1A,L2A,N1A,N2A)
C-----IV N1 > 0, N2 = 0
      IF(N1A.LE.0.OR.N2A.NE.0) GO TO 30
C
C     ENDF/B-IV FORMAT. UPDATE NUMBER OF COMMENT CARDS AND OUTPUT
C     TWO CARDS SET FLAGS FOR ENDF/B-IV FORMAT.
C
C-----IF THERE ARE RESONANCE PARAMETERS AND THEIR CONTRIBUTION HAS NOT
C-----BEEN ADDED TO THE BACKGROUND SET LRP TO INDICATE THAT IT HAS NOW
C-----BEEN ADDED. ONLY ADD COMMENT CARDS IF MAT WILL BE PROCESSED.
      IF(LRP.NE.1) GO TO 10
      L1H=2
      N1OUT=N1A+3
      IF(MAXER2.GT.1) N1OUT=N1OUT+MAXER2+3
      GO TO 20
   10 N1OUT=N1A
   20 CALL CONTO
      CALL CARDO(C1A,C2A,L1A,L2A,N1OUT,N2A)
      N1X=N1A
      FMTHOL=FMTTAB(1)
      IVERSE=4
      INT45=5
      GO TO 90
C-----NOT ENDF/B-IV. READ THIRD LINE AND TEST FOR ENDF/B-V FORMAT.
   30 CALL CARDI(C1B,C2B,L1B,L2B,N1B,N2B)
      IF(N2A.GT.0) GO TO 60
C
C     ENDF/B-V FORMAT. UPDATE NUMBER OF COMMENT CARDS AND OUTPUT
C     THREE CARDS AND SET FLAGS FOR ENDF/B=V FORMAT.
C
C-----IF THERE ARE RESONANCE PARAMETERS AND THEIR CONTRIBUTION HAS NOT
C-----BEEN ADDED TO THE BACKGROUND SET LRP TO INDICATE THAT IT HAS NOW
C-----BEEN ADDED. ONLY ADD COMMENT CARDS IF MAT WILL BE PROCESSED.
      IF(LRP.NE.1) GO TO 40
      L1H=2
      N1OUT=N1B+3
      IF(MAXER2.GT.1) N1OUT=N1OUT+MAXER2+3
      GO TO 50
   40 N1OUT=N1B
   50 CALL CONTO
      CALL CARDO(C1A,C2A,L1A,L2A,N1A,N2A)
      CALL CARDO(C1B,C2B,L1B,L2B,N1OUT,N2B)
      N1X=N1B
      FMTHOL=FMTTAB(2)
      IVERSE=5
      INT45=2
      GO TO 90
C
C     ENDF/B-VI FORMAT. UPDATE NUMBER OF COMMENT CARDS AND OUTPUT
C     THREE CARDS AND SET FLAGS FOR ENDF/B=VI FORMAT.
C
C-----READ FOURTH LINE.
   60 CALL CARDI(C1C,C2C,L1C,L2C,N1C,N2C)
C-----IF THERE ARE RESONANCE PARAMETERS AND THEIR CONTRIBUTION HAS NOT
C-----BEEN ADDED TO THE BACKGROUND SET LRP TO INDICATE THAT IT HAS NOW
C-----BEEN ADDED. ONLY ADD COMMENT CARDS IF MAT WILL BE PROCESSED.
      IF(LRP.NE.1) GO TO 70
      L1H=2
      N1OUT=N1C+3
      IF(MAXER2.GT.1) N1OUT=N1OUT+MAXER2+3
      GO TO 80
   70 N1OUT=N1C
   80 CALL CONTO
      CALL CARDO(C1A,C2A,L1A,L2A,N1A,N2A)
      CALL CARDO(C1B,C2B,L1B,L2B,N1B,N2B)
      N1X=N1C
      FMTHOL=FMTTAB(3)
      IVERSE=6
      INT45=2
      TEMP3=C1C
      INPART=N1B/10
C-----SET DERIVED MATERIAL FLAG.
      L1C=1
      CALL CARDO(C1C,C2C,L1C,L2C,N1OUT,N2C)
C-----IF NO RESONANCE PARAMETERS OR THERE CONTRIBUTION HAS ALREADY BEEN
C-----ADDED TO THE BACKGROUND COPY MAT. OTHERWISE INSERT COMMENTS.
   90 IF(LRP.NE.1) GO TO 130
C
C     SKIP OUTPUT IF IN EDIT (I.E. NO OUTPUT) MODE.
C
      IF(OTAPE.LE.0) GO TO 140
C-----COPY TO END OF HOLLERITH.
      DO 100 N=1,N1X
  100 CALL COPY1
C
C     ADD COMMENTS TO DOCUMENT WHAT WAS DONE TO DATA
C
C-----OUTPUT PROGRAM VERSION I.D.
      CALL HOLLYO(PROGDOC1(1,1))
C
C     DESCRIBE RESONANCE RECONSTRUCTION CRITERIA.
C
C-----OUTPUT MINIMUM CROSS SECTION
      CALL OUT9(SIGMIN,PROGDOC1(27,3))
      CALL HOLLYO(PROGDOC1(1,3))
      IF(MAXER2.GT.1) GO TO 110
C-----OUTPUT ENERGY INDEPENDENT ERROR USED FOR RECONSTRCUTION.
      PERCNT=100.0*ER2(1)
      CALL OUT9(PERCNT,PROGDOC1(42,4))
      CALL HOLLYO(PROGDOC1(1,4))
      GO TO 140
C-----OUTPUT 4 LINE TITLE
  110 CALL HOLLYO(PROGDOC1(1,5))
      CALL HOLLYO(PROGDOC1(1,6))
      CALL HOLLYO(PROGDOC1(1,7))
      CALL HOLLYO(PROGDOC1(1,8))
      DO 120 I=1,MAXER2
      PERCNT=100.0*ER2(I)
      CALL OUT9(ENER2(I),PROGDOC1( 2,9))
      CALL OUT9(PERCNT  ,PROGDOC1(14,9))
  120 CALL HOLLYO(PROGDOC1(1,9))
      GO TO 140
C
C     END OF HOLLERITH OUTPUT TO FILE1
C
C-----INDICATE THAT MAT NEED NOT BE PROCESSED (EITHER NO PARAMETERS OR
C-----THEIR CONTRIBUTION HAS ALREADY BEEN ADDED TO THE BACKGROUND).
  130 IMDONE=0
C-----IDENTIFY ENDF/B FORMAT.
  140 WRITE(OUTP,160) FMTHOL
      WRITE(*   ,160) FMTHOL
C-----DEFINE WHETHER OR NOT MATERIAL IS FISSILE.
      IF(LFI.EQ.0) WRITE(OUTP,170)
      IF(LFI.GT.0) WRITE(OUTP,180)
      IF(LFI.EQ.0) WRITE(*   ,170)
      IF(LFI.GT.0) WRITE(*   ,180)
C-----DEFINE WHETHER OR NOT THERE ARE RESONANCE PARAMETERS AND WHETHER
C-----THERE CONTRIBUTION HAS ALREADY BEEN ADDED TO THE BACKGROUND CROSS
C-----SECTION.
      IF(LRP.EQ.0) WRITE(OUTP,190)
      IF(LRP.EQ.1) WRITE(OUTP,200)
      IF(LRP.EQ.2.AND.IMEDIT.NE.2) WRITE(OUTP,210)
      IF(LRP.EQ.2.AND.IMEDIT.EQ.2) WRITE(OUTP,220)
      IF(LRP.EQ.0) WRITE(*   ,190)
      IF(LRP.EQ.1) WRITE(*   ,200)
      IF(LRP.EQ.2.AND.IMEDIT.NE.2) WRITE(*   ,210)
      IF(LRP.EQ.2.AND.IMEDIT.EQ.2) WRITE(*   ,220)
C-----FOR ENDF/B-VI FORMATTED DATA DEFINE PROJECTILE AND BACKGROUND
C-----TEMPERATURE.
      IF(IVERSE.NE.6) GO TO 150
      IF(INPART.EQ.1) WRITE(OUTP,230) INPART
      IF(INPART.NE.1) WRITE(OUTP,240) INPART
      CALL OUT9(TEMP3,FIELD(1,1))
      WRITE(OUTP,250) (FIELD(M,1),M=1,11)
  150 RETURN
  160 FORMAT(' Based on the Format and Contents of MF=1, MT=451'/
     1 ' (1) ENDF/B-',A2,' Format.')
  170 FORMAT(' (2) Material is NOT Fissile (LFI=0).')
  180 FORMAT(' (2) Material is Fissile (LFI=1).')
  190 FORMAT(' (3) No Resonance Parameters Given (LRP=0).',
     1 ' MAT Copied.')
  200 FORMAT(' (3) Resonance Parameters are Given (LRP=1).')
  210 FORMAT(' (3) Resonance Data Already Added to Background Cross'/
     1       '     Section (LRP=2). MAT Copied.')
  220 FORMAT(' (3) Resonance Data Already Added to Background Cross'/
     1       '     Section (LRP=2).')
  230 FORMAT(' (4) Projectile ZA =',I7,' (Neutron).')
  240 FORMAT(' (4) Projectile ZA =',I7,' (Expect Neutron = 1).')
  250 FORMAT(' (5) Temperature of Background ',11A1,' Kelvin.')
      END
CAK   SUBROUTINE FILE2
      SUBROUTINE FILE2r
C=======================================================================
C
C     READ ALL FILE2 DATA AND PROCESS INTO POINTWISE FORM.
C
C=======================================================================
      INCLUDE 'implicit.h'
      REAL*4 SECONDS
      CHARACTER*1 FIELD
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/UNITS/ISCR2,ISCR23
      COMMON/SUBS/ESUB(1000),ENODP,ENODM,WIDP,WIDM,ISUB,NSUB
      COMMON/PAGER/NPAGE,NPAGP1,NPAGM1
      COMMON/POINTN/NPOINT,KPOINT
      COMMON/OKERR2/ERRXC2,ERXC20,ENER2(1000),ER2(1000),MAXER2,NOERR
      COMMON/MINSIG/SIGMIN
      COMMON/MINNIE/EMIN,EMAX,DEMIN
      COMMON/EDGES/ERANGE(50),NODLOW(50),NODHI(50),NRANGE,IRANGE
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FIELDC/FIELD(11,12)
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
      COMMON/NEG1COM/NEG1    ! 2014/11/02 - add to FILE2 and FILE3
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DATA IUSE2/0/
      DATA TINY/1.0D-09/
c-----reset LRF=7 used and IMFISSY flags for each evaluation
      call same0
C-----READ ALL FILE 2 PARAMETERS.
      CALL READ2
c-----DFEFINE REACTIONS TO CALCULATE
      call Answer1
C-----INITIALIZE POINT COUNTS
C-----NPOINT = NUMBER OF POINTS ON SCRATCH. POINT TOTAL AT END.
C-----KPOINT = NUMBER OF POINTS IN CORE (0 TO MAXRES).
C-----KPTP1  = INDEX TO POSITION INTO WHICH NEXT POINT WILL BE STORED.
C-----NBASE  = NUMBER OF POINTS GENERATED UP TO END OF LAST ENERGY.
C-----RANGE).
      NPOINT=0
      KPOINT=0
      KPTP1=1
      NBASE=0
      NEG1 =0    ! added 2014/11/02 - to initialize NEG1
C-----IF THERE ARE NO SECTIONS WHOSE RESONANCE CONTRIBUTION MUST BE
C-----ADDED TO BACKGROUND CROSS SECTION THERE IS NOTHING ELSE TO DO.
      IF(NSECT.GT.0) GO TO 10
      WRITE(OUTP,500)
      WRITE(*   ,500)
      GO TO 430
C-----IF IN EDIT MODE RETURN AFTER READING ALL FILE 2 DATA.
   10 IF(IMEDIT.EQ.2) GO TO 430
C
C     CALCULATE RESONANCE CONTRIBUTION TO CROSS SECTION.
C
      WRITE(OUTP,490)
      WRITE(*   ,490)
C-----INITIALIZE SAVED POINT COUNT.
      NSAVE=0
C
C     SET UP LOOP OVER ENERGY RANGES (E.G. RESOLVED AND UNRESOLVED)
C
      DO 420 IRANGE=2,NRANGE
      IRM1=IRANGE-1
C-----INITIALIZE FLAG TO DEFINE TYPE OF RANGE (RESOLVED, UNRESOLVED
C-----OR BOTH).
      IR=0
      IU=0
C-----DEFINE WHICH SECTIONS CONTRIBUTE TO THIS ENERGY RANGE.
      DO 30 I=1,NSECT
      IF(EL(I).LT.ERANGE(IRANGE).AND.EH(I).GT.ERANGE(IRM1)) GO TO 20
C-----TURN OFF SECTION.
      MODE(I)=-IABS(MODE(I))
      GO TO 30
C-----TURN ON SECTION.
   20 MODE(I)=IABS(MODE(I))
C-----DETERMINE TYPE OF RESONANCE REGION (I.E. RESOLVED, UNRESOLVED,
C-----BOTH OR NONE).
      IF(MODE(I).LT.11) IR=1
      IF(MODE(I).GE.11) IU=1
   30 CONTINUE
C-----IDENTIFY TYPE OF ENERGY REGION.
      NREG=1+IR+2*IU
C
C     SET UP LOOP OVER ENERGY NODES IN CURRENT ENERGY RANGE.
C
C-----DEFINE WHICH NODES LIE IN THIS ENERGY RANGE AND SET UP LOOP OVER
C-----NODES.
      NODE1=NODLOW(IRANGE)+1
      NODE2=NODHI(IRANGE)
      DO 280 KNODE=NODE1,NODE2
C-----DEFINE INTERVAL END POINTS.
      ENODM=ENODE(KNODE-1)
      WIDM=WIDNOD(KNODE-1)
      ENODP=ENODE(KNODE)
      WIDP=WIDNOD(KNODE)
C-----DEFINE SUB-INTERVALS.
      CALL SUBINT
C
C     SPECIAL TREATMENT FOR FIRST ENERGY POINT. IF FIRST RESONANCE
C     REGION STARTS ABOVE 1.0E-5 EV START TABLE AT LOWER LIMIT OF
C     RESONANCE REGION WITH ONE POINT WITH ZERO CROSS SECTION FOLLOWED
C     BY A POINT WITH THE SAME ENERGY AND THE CROSS SECTION CALCULATED
C     AT THE LOWER LIMIT OF THE RESONANCE REGION. IF FIRST RESONANCE
C     STARTS AT 1.0E-5 START TABLE WITH CROSS SECTION CALCULATED AT
C     LOWER ENERGY LIMIT OF RESONANCE REGION.
C
      IF(KPOINT.LE.0) GO TO 40
C-----CALCULATE CROSS SECTION AT FIRST POINT OF EACH ENERGY REGION.
      IF(KNODE.le.NODE1) go to 60
      go to 70
C-----INITIALIZE CROSS SECTION TABLE TO ONE POINT WITH ZERO CROSS
C-----SECTION AT THE BEGINNING OF THE TABLE UNLESS THE TABLE STARTS
C-----AT 1.0-5 EV OR LESS (IN WHICH CASE TABLE WILL START WITH NON-ZERO
C-----CROSS SECTIONS).
   40 IF(ERANGE(1).LE.1.001*EMIN) GO TO 60
      KPOINT=1
      KPTP1=2
      ETAB2(1)=ERANGE(1)
      DO 50 IR=1,NREACT
   50 SIG2(IR,1)=0.0
C-----INITIALIZE INDEX TO FIRST SUB-INTERVAL AND INDICATE THAT THE
C-----THE SUB-INTERVAL HAS NOT YET BEEN SUB-DIVIDED.
   60 ISUB=1
C-----CALCULATE CROSS SECTION AT NEXT ENERGY (EITHER FIRST OR SECOND
C-----ENERGY).
      ETAB2(KPTP1)=ESUB(1)
      CALL SIGMA(ETAB2(KPTP1),SIG2(1,KPTP1))
C-----SAVE STARTING POINT (NO CONVERGENCE TESTS UNTIL SECOND POINT IS
C-----GENERATED).
      GO TO 160
C
C     INTERNAL SUB-INTERVAL. SET UP LOOP OVER SUB-INTERVALS REMAINING
C     SUB-INTERVALS (NOTE, AT THIS POINT THE ENERGY INTERVAL UP TO AND
C     INCLUDING THE ENERGY OF THE FIRST SUB-INTERVAL (ISUB=1) HAS
C     ALREADY BEEN CALCULATED, EITHER BY CALCULATING THE FIRST POINT
C     OR AS THE LAST POINT OF THE PRECEDING ENERGY INTERVAL).
C
C-----INITIALIZE INDEX TO SECOND SUB-INTERVAL AND INDICATE THAT THE
C-----SUB-INTERVAL HAS NOT YET BEEN SUB-DIVIDED.
   70 ISUB=2
C-----CALCULATE CROSS SECTION AT END OF SUB-INTERVAL.
   80 ETAB2(KPTP1)=ESUB(ISUB)
      CALL SIGMA(ETAB2(KPTP1),SIG2(1,KPTP1))
C-----SKIP IF NO ITERATION
      IF(NOERR.GT.0) GO TO 160
C
C     ITERATION AND CONVERGENCE TESTS.
C
C-----DEFINE ENERGY AT MIDPOINT AND TEST FOR CONVERGENCE BASED ON SHORT
C-----ENERGY INTERVAL.
   90 EA=ETAB2(KPOINT)
      EB=ETAB2(KPTP1)
      EMID=HALF*(EA+EB)
      CALL INCORE9(EMID)
      IF(EMID.LE.EA.OR.EMID.GE.EB) GO TO 160
      IF((EMID-EA).LE.EA*TINY) GO TO 160
C-----DEFINE CROSS SECTION AT MIDPOINT.
      CALL SIGMA(EMID,SIGMID)
C-----IF AN ENERGY DEPENDENT ERROR LAW IS USED DEFINE ERROR AT EMID.
  100 IF(MAXER2.GT.1) CALL ERROK2(EMID)
C-----DEFINE CONTRIBUTION OF EACH ENDPOINT TO MIDPOINT.
      DE=EB-EA
      WTA=(EB-EMID)/DE
      WTB=(EMID-EA)/DE
C-----TEST EACH REACTION FOR CONVERGENCE. TO BE ACCEPTABLE ALL REACTIONS
C-----MUST CONVERGE (I.E., IF EVEN ONE REACTION FAIL ONE CONVEREGENCE
C-----TEST CONTINUE SUB-DIVIDING AND INTERATING TO CONVERGENCE).
      DO 130 IR=1,NREACT
C-----DEFINE EXACT END POINT CROSS SECTIONS.
      SIGA=SIG2(IR,KPOINT)
      SIGB=SIG2(IR,KPTP1)
C-----CONVERGENCE IF CROSS SECTION AT BOTH ENDS OF INTERVAL ARE LESS
C-----THAN MINIMUM CROSS SECTION OF INTEREST.
      IF(DABS(SIGA).LE.SIGMIN.AND.DABS(SIGB).LE.SIGMIN) GO TO 130
C-----TO PREVENT INFINITE ITERATION TOWARD ZERO ONLY APPLY CHANGE TEST
C-----IF CROSS SECTION AT BOTH ENDS OF INTERVAL ARE POSITIVE.
      IF(SIGA.LE.0.0.OR.SIGB.LE.0.0) GO TO 110
C-----TEST FOR CROSS SECTION CHANGE.
      IF(SIGA.GT.1.4*SIGB.OR.SIGB.GT.1.4*SIGA) GO TO 140
C-----DEFINE EXACT AND LINEARLY INTERPOLATED MID-POINT CROSS SECTIONS.
  110 SIGM=SIGMID(IR)
      SIGLIN=WTA*SIGA+WTB*SIGB
C-----09/22/02 - ADDED TEST FOR MAXIMUM
C-----TEST FOR ITERATION TOWARD MIN/MAX (MINIMUM IS INDICATED IF
C-----THE EXACT CROSS SECTION AT THE MID-POINT IS LESS THAN THE
C-----EXACT CROSS SECTION AT BOTH ENDS OF THE INTERVAL).
      IF(SIGM.LT.SIGA.AND.SIGM.LT.SIGB) GO TO 120
      IF(SIGM.GT.SIGA.AND.SIGM.GT.SIGB) GO TO 120
C-----NORMAL CONVERGENCE TEST.
      IF(DABS(SIGM-SIGLIN).le.DABS(ERRXC2*SIGM)) go to 130
      go to 140
C-----ITERATING TOWARD MIN/MAX. USE MORE STRINGENT CONVERGENCE TEST.
  120 IF(DABS(SIGM-SIGLIN).gt.DABS(ERXC20*SIGM)) go to 140
C-----END OF CONVERGENCE TEST LOOP. REACTION HAS PASSED ALL TESTS AND
C-----FOR THIS REACTION MID-POINT IS NOT REQUIRED.
  130 CONTINUE
C-----CONVERGENCE FOR ALL REACTIONS. KEEP END POINT.
      GO TO 160
C
C     NO CONVERGENCE. SAVE VALUE FROM END OF INTERVAL AND SHORTEN
C     INTERVAL.
C
  140 IF(NSAVE.LT.MAXSAVE) NSAVE=NSAVE+1
      ESAVE(NSAVE)=ETAB2(KPTP1)
      ETAB2(KPTP1)=EMID
      DO 150 IR=1,NREACT
      SIGSAVE(IR,NSAVE)=SIG2(IR,KPTP1)
  150 SIG2(IR,KPTP1)=SIGMID(IR)
      GO TO 90
C
C     CONVERGENCE.
C
C-----SAVE END OF INTERVAL BY INCREASING KPOINT BY ONE. IF IN CORE PAGE
C-----IS FULL UNLOAD IT TO SCRATCH, MOVE LAST POINT TO TABLE BEGINNING
C-----AND RE-INITIALIZE IN CORE POINT INDICES.
  160 IF(KPOINT.LT.NPAGE) GO TO 230
      IF(NPOINT.EQ.0.AND.IUSE2.GT.0) REWIND ISCR2
C
C     DEFINE START OF EACH REACTION
C
      IF(NEG1.NE.0) GO TO 200
      NEG1 = 1
      DO 190 II=1,NREACT
      IF(NEGTAB(II).NE.0) GO TO 190
      DO 170 I=1,NPAGE
      IF(SIG2X(II,I).GT.0.0D+0) GO TO 180
  170 CONTINUE
      NEG1 = 0
      GO TO 190
  180 NEGTAB(II) = I + NPOINT
      IF(I.GT.1) NEGTAB(II) = I - 1
  190 CONTINUE
C
C     OUTPUT TO SCRATCH.
C
  200 WRITE(ISCR2) ETAB2X,SIG2X
      NPOINT=NPOINT+NPAGE
C-----IF REQUESTED PRINT MESSAGE EVERYTIME A PAGE IS OUTPUT TO SCRATCH.
      IF(MONITR.EQ.0) GO TO 210
      CALL TIMER1(SECONDS)
      CALL OUT9(ETAB2(    1),FIELD(1,1))
      CALL OUT9(ETAB2(NPAGE),FIELD(1,2))
      WRITE(OUTP,510) NPOINT,((FIELD(M,I),M=1,11),I=1,2),SECONDS
      WRITE(*   ,510) NPOINT,((FIELD(M,I),M=1,11),I=1,2),SECONDS
  210 ETAB2(1)=ETAB2(NPAGP1)
      DO 220 IR=1,NREACT
  220 SIG2(IR,1)=SIG2(IR,NPAGP1)
      KPTP1=1
  230 KPOINT=KPTP1
      KPTP1=KPOINT+1
C
C     IF THERE ARE NO REMAINING SAVED POINTS THE CURRENT SUB-INTERVAL
C     IS NOW FINISHED. IF THERE ARE REMAINING SAVED POINTS USE THEM
C     TO DEFINE NEXT ITERATION INTERVAL.
C
      IF(NSAVE.lt.1) go to 270
      IF(NSAVE.gt.1) go to 250
C-----ONLY 1 POINT SAVED...USE IT TO DEFINE ENDPOINT AND BRANCH BACK
C-----TO DEFINE MIDPOINT...THEN ITERATE.
      ETAB2(KPTP1)=ESAVE(1)
      DO 240 IR=1,NREACT
  240 SIG2(IR,KPTP1)=SIGSAVE(IR,1)
      NSAVE=0
      GO TO 90
C-----MORE THAN 1 POINT SAVED...USE LAST 2 SAVED POINTS TO DEFINE
C-----MIDPOINT AND ENDPOINT...THEN ITERATE.
  250 NSAVM1=NSAVE-1
      ETAB2(KPTP1)=ESAVE(NSAVM1)
      EMID=ESAVE(NSAVE)
      DO 260 IR=1,NREACT
      SIG2(IR,KPTP1)=SIGSAVE(IR,NSAVM1)
  260 SIGMID(IR)=SIGSAVE(IR,NSAVE)
      NSAVE=NSAVE-2
C-----CONVERGENCE IF ENERGY INTERVAL IS TOO SMALL.
      EA=ETAB2(KPOINT)
      EB=ETAB2(KPTP1)
      IF((EMID-EA).LE.EA*TINY) GO TO 160
C-----USE 3 CURRENTLY KNOWN POINTS TO TEST FOR CONVERGENCE.
      GO TO 100
C
C     END OF SUB-INTERVAL LOOP. CONTINUE UNTIL ALL SUB-INTERVALS HAVE
C     BEEN USED.
C
  270 ISUB=ISUB+1
      IF(ISUB.LE.NSUB) GO TO 80
C
C     END OF NODE LOOP.
C
  280 CONTINUE
C
C     END OF 1 RESONANCE REGION. IS THIS THE END OF THE ENTIRE RESONANCE
C     REGION (E.G., RESOLVED AND UNRESOLVED).
C
      IF(IRANGE.LT.NRANGE) GO TO 410
C
C     END OF ENTIRE RESONANCE REGION. ADD LAST POINT WITH ZERO CROSS
C     SECTION AND LEAVE DATA IN CORE OR MOVE TO SCRATCH.
C
C-----IF IN CORE PAGE IS FULL UNLOAD IT TO SCRATCH.
      IF(KPOINT.LT.NPAGE) GO TO 340
      IF(NPOINT.LE.0.AND.IUSE2.GT.0) REWIND ISCR2
C
C     DEFINE START OF EACH REACTION.
C
      IF(NEG1.NE.0) GO TO 320
      NEG1 = 1
      DO 310 II=1,NREACT
      IF(NEGTAB(II).NE.0) GO TO 310
      DO 290 I=1,NPAGE
      IF(SIG2X(II,I).GT.0.0D+0) GO TO 300
  290 CONTINUE
      NEG1 = 0
      GO TO 310
  300 NEGTAB(II) = I + NPOINT
      IF(I.GT.1) NEGTAB(II) = I - 1
  310 CONTINUE
C
C     OUTPUT TO SCRATCH
C
  320 WRITE(ISCR2) ETAB2X,SIG2X
      NPOINT=NPOINT+NPAGE
C-----IF REQUESTED PRINT MESSAGE EVERYTIME A PAGE IS OUTPUT TO SCRATCH.
      IF(MONITR.EQ.0) GO TO 330
      CALL TIMER1(SECONDS)
      CALL OUT9(ETAB2(    1),FIELD(1,1))
      CALL OUT9(ETAB2(NPAGE),FIELD(1,2))
      WRITE(OUTP,510) NPOINT,((FIELD(M,I),M=1,11),I=1,2),SECONDS
      WRITE(*   ,510) NPOINT,((FIELD(M,I),M=1,11),I=1,2),SECONDS
  330 KPTP1=1
  340 KPOINT=KPTP1
      KPTP1=KPOINT+1
      ETAB2(KPOINT)=ERANGE(IRANGE)
      DO 350 IR=1,NREACT
  350 SIG2(IR,KPOINT)=0.0
C-----------------------------------------------------------------------
C
C     END OF RESONANCE REGION
C
C-----------------------------------------------------------------------
C
C     DEFINE START OF EACH REACTION
C
      IF(NEG1.NE.0) GO TO 390
      NEG1 = 1
      DO 380 II=1,NREACT
      IF(NEGTAB(II).NE.0) GO TO 380
      DO 360 I=1,KPOINT
      IF(SIG2X(II,I).GT.0.0D+0) GO TO 370
  360 CONTINUE
      NEG1 = 0
      GO TO 380
  370 NEGTAB(II) = I + NPOINT
      IF(I.GT.1) NEGTAB(II) = I - 1
  380 CONTINUE
C
C     DEFINE TOTAL POINT COUNT.
C
  390 NPOINT=NPOINT+KPOINT
C-----ALL POINTS HAVE BEEN CALCULATED. LEAVE IN CORE OR SAVE ON SCRATCH.
      IF(NPOINT.LE.NPAGE) GO TO 400
C-----OUTPUT TO SCRATCH
      WRITE(ISCR2) ETAB2X,SIG2X
      END FILE ISCR2
      REWIND ISCR2
      IUSE2=1
C-----IF REQUESTED PRINT MESSAGE EVERYTIME A PAGE IS OUTPUT TO SCRATCH.
      IF(MONITR.EQ.0) GO TO 400
      CALL TIMER1(SECONDS)
      CALL OUT9(ETAB2(     1),FIELD(1,1))
      CALL OUT9(ETAB2(KPOINT),FIELD(1,2))
      WRITE(OUTP,510) NPOINT,((FIELD(M,I),M=1,11),I=1,2),SECONDS
      WRITE(*   ,510) NPOINT,((FIELD(M,I),M=1,11),I=1,2),SECONDS
C-----SET IN CORE POINT COUNT TO ZERO (NCOUNT IS NOW THE TOTAL POINT
C-----COUNT).
  400 KPOINT=0
C-----PRINT SUMMARY OF ENERGY RANGE.
  410 NKPONT=NPOINT+KPOINT
      NBASE=NKPONT-NBASE
      CALL OUT9(ERANGE(IRM1  ),FIELD(1,1))
      CALL OUT9(ERANGE(IRANGE),FIELD(1,2))
      IF(NREG.EQ.1) WRITE(OUTP,440)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.2) WRITE(OUTP,450)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.3) WRITE(OUTP,460)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.4) WRITE(OUTP,470)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.1) WRITE(*   ,440)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.2) WRITE(*   ,450)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.3) WRITE(*   ,460)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      IF(NREG.EQ.4) WRITE(*   ,470)
     1 ((FIELD(M,I),M=1,11),I=1,2),NBASE
      NBASE=NKPONT
C-----END OF ENERGY RANGE LOOP.
  420 CONTINUE
C
C     ALL CROSS SECTION POINTS HAVE NOW BEEN CALCULATED. PRINT SUMMARY
C     OF ENTIRE RESONANCE REGION.
C
C-----OUTPUT SUMMARY OF TOTAL NUMBER OF POINTS GENERATED.
      WRITE(OUTP,480) NPOINT
      WRITE(*   ,480) NPOINT
  430 RETURN
  440 FORMAT(2X,22A1,I11,' Not in Any Resonance Region'/
     1 35X,' (WARNING - Not Expected - Check Data)')
  450 FORMAT(2X,22A1,I11,' Resolved')
  460 FORMAT(2X,22A1,I11,' Unresolved')
  470 FORMAT(2X,22A1,I11,' Overlapping Resolved/Unresolved'/
     1 35X,' (WARNING - This is Illegal in ENDF/B)')
  480 FORMAT(1X,78('=')/' Entire Resonance Region',I11,' Points'/
     1 1X,78('='))
  490 FORMAT(1X,78('=')/' Reconstructing Cross Sections from',
     1 ' Resonance Parameters'/1X,78('=')/
     2 '        E-Low     E-High   Points  ',
     3 ' Type of Resonance Region'/
     3 '        (eV)       (eV)   Generated',
     4 ' Messages'/1X,78('='))
  500 FORMAT(1X,78('=')/' No Resonance Contribution to Add to',
     1 ' Background Cross Section. MAT Copied.')
  510 FORMAT(25X,I10,11A1,' to',11A1,' eV',F11.2,' Sec.')
      END
      SUBROUTINE ERROK2(E)
C=======================================================================
C
C     DEFINE ALLOWABLE ERROR FOR RECONSTRUCTION OF ENERGY DEPENDENT
C     CROSS SECTIONS FROM FILE 2 RESONANCE PARAMETERS. THE ERROR LAW
C     CAN BE ENERGY INDEPENDENT (CONSTANT) OR ENERGY DEPENDENT
C     (GIVEN BY A LINEARLY INTERPOLABLE TABLE IN ENERGY VS. ERROR).
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/OKERR2/ERRXC2,ERXC20,ENER2(1000),ER2(1000),MAXER2,NOERR
C-----INITIALIZE INDEX TO INTERPOLATION TABLE.
      DATA MINER2/2/
C-----INTERPOLATE TO FIND FIRST TABULATED ENERGY ABOVE E.
      IF(E.le.ENER2(1)) go to 80
      DO 10 NOWER2=MINER2,MAXER2
      IF(E.lt.ENER2(NOWER2)) go to 20
      IF(E.eq.ENER2(NOWER2)) go to 70
   10 CONTINUE
C-----EXTEND ERROR AS CONSTANT ABOVE TABULATED RANGE.
      GO TO 90
   20 NOWER1=NOWER2-1
      IF(E.eq.ENER2(NOWER1)) go to 60
      IF(E.gt.ENER2(NOWER1)) go to 50
      DO 30 NOWER2=2,MAXER2
      IF(E.lt.ENER2(NOWER2)) go to 40
      IF(E.eq.ENER2(NOWER2)) go to 70
   30 CONTINUE
      GO TO 90
C-----DEFINE INDEX AND INTERPOLATE ERROR.
   40 NOWER1=NOWER2-1
   50 ERRXC2=((ENER2(NOWER2)-E)*ER2(NOWER1)+
     1 (E-ENER2(NOWER1))*ER2(NOWER2))/(ENER2(NOWER2)-ENER2(NOWER1))
C-----09/21/02 - DECREASED FROM 0.10 TO 0.02
C     ERXC20=0.10*ERRXC2
      ERXC20=0.02*ERRXC2
      MINER2=NOWER2
      RETURN
C-----EXACT MATCH TO TABULATED ENERGY. DEFINE INDEX AND ERROR.
   60 ERRXC2=ER2(NOWER1)
C-----09/21/02 - DECREASED FROM 0.10 TO 0.02
C     ERXC20=0.10*ERRXC2
      ERXC20=0.02*ERRXC2
      MINER2=NOWER1
      IF(MINER2.LE.1) MINER2=2
      RETURN
C-----EXACT MATCH TO TABULATED ENERGY. DEFINE INDEX AND ERROR.
   70 ERRXC2=ER2(NOWER2)
C-----09/21/02 - DECREASED FROM 0.10 TO 0.02
C     ERXC20=0.10*ERRXC2
      ERXC20=0.02*ERRXC2
      MINER2=NOWER2
      RETURN
C-----EXTEND ERROR AS CONSTANT BELOW TABULATED RANGE.
   80 ERRXC2=ER2(1)
C-----09/21/02 - DECREASED FROM 0.10 TO 0.02
C     ERXC20=0.10*ERRXC2
      ERXC20=0.02*ERRXC2
      MINER2=2
      RETURN
C-----EXTEND ERROR AS CONSTANT ABOVE TABULATED RANGE.
   90 ERRXC2=ER2(MAXER2)
C-----09/21/02 - DECREASED FROM 0.10 TO 0.02
C     ERXC20=0.10*ERRXC2
      ERXC20=0.02*ERRXC2
      MINER2=MAXER2
      RETURN
      END
      SUBROUTINE READ2
C=======================================================================
C
C     READ ALL FILE2 DATA. DEFINE ENERGY RANGES AND ENERGY NODES WITHIN
C     EACH ENERGY RANGE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 ZABCDI,FIELD
      CHARACTER*28 FISLST,REGLST,RTYPE,UTYPE,HOLNRO,HOLNAP
      COMMON/WHATZA/IZANOW,MATNOW,TEMP3,IVERSE,INT45
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/NAPRHO/NRO,NAPS
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/ETABS/NODES
      COMMON/RANGER/LOW,LHI
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/MINNIE/EMIN,EMAX,DEMIN
      COMMON/EDGES/ERANGE(50),NODLOW(50),NODHI(50),NRANGE,IRANGE
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/LRUNOW/LRUIN
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FIELDC/FIELD(11,12)
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION ZABCDI(10),FISLST(3),REGLST(4),RTYPE(8),
     1 UTYPE(4),HOLNRO(3),HOLNAP(3)
      DATA FISLST/
     1 ' (No Fission Widths)        ',
     2 ' (Fission Widths Given)     ',
     3 ' (ERROR - MUST be 0 or 1)   '/
      DATA REGLST/
     1 ' (No Parameters)            ',
     2 ' (Resolved Region)          ',
     3 ' (Unresolved Region)        ',
     4 ' (ERROR - See Below)        '/
      DATA RTYPE/
     1 ' (Single Level Breit-Wigner)',
     2 ' (Multi-Level Breit-Wigner) ',
     3 ' (Reich-Moore)              ',
     4 ' (Adler-Adler)              ',
     5 ' (General R-Matrix)         ',
     6 ' (Hybrid R-Function)        ',
     7 ' (Reich-Moore + Competition)',
     8 ' (ERROR - see Below)        '/
      DATA UTYPE/
     1 ' (All Energy Independent)   ',
     2 ' (Fission Energy Dependent) ',
     3 ' (All Energy Dependent)     ',
     4 ' (ERROR - see Below)        '/
      DATA HOLNRO/
     1 '(Energy Independent)        ',
     2 '(Energy Dependent)          ',
     3 '(ERROR - MUST be 0 or 1)    '/
      DATA HOLNAP/
     1 '(Calculate)                 ',
     2 '(= Scattering Radius)       ',
     3 '(ERROR - MUST be 0 to 2)    '/
C-----DEFINE THERMAL ENERGY FOR INCLUSION AS A NODE IF THE RESONANCE
C-----REGION SPANS THERMAL ENERGY.
      DATA ETHERM/2.53D-02/
C-----DEFINE MINIMUM ENERGY OF RESONANCE REGION (ALLOWING FOR ROUNDOFF).
      EMINL=0.999*EMIN
C-----FILE 2 FOUND. DEFINE ZA, AWR AND NIS (NUMBER OF ISOTOPES).
      ZA=C1H
      AWR=C2H
      NIS=N1H
C-----INITIALIZE FISSILE FLAG OFF (INDICATES ALL FISSION WIDTHS READ SO
C-----FAR ARE ZERO - IF THIS FLAG IS NOT RESET WHILE READING PARAMETERS
C-----THE MAT WILL BE CONSIDERED TO BE NON-FISSILE AND FISSION CROSS
C-----SECTIONS WILL NOT BE OUTPUT).
      LFWX  =0
      LFWSUM=0
C-----PRINT ZA, ATOMIC WEIGHT AND NUMBER OF ISOTOPES.
      IZA=ZA
      CALL ZAHOL(IZA,ZABCDI)
      CALL OUT9(AWR,FIELD(1,1))
      WRITE(OUTP,410) ZABCDI,(FIELD(M,1),M=1,11),NIS
      WRITE(*   ,410) ZABCDI,(FIELD(M,1),M=1,11),NIS
C
C     INITIALIZE RESONANCE PARAMETER TABLE INDICES.
C
C-----INITIALIZE SECTION COUNT (A SECTION IS DATA FOR ONE ISOTOPE,
C-----ENERGY RANGE AND L VALUE) AND THE NUMBER OF ENERGY RANGES WITH
C-----ENERGY DEPENDENT SCATTERING RADIUS. INITIALIZE ALL SECTION
C-----INDICES TO INDICATE NO ENERGY DEPENDENT SCATTERING RADIUS.
      NSECT=0
      NUMRHO=0
      DO 10 I=1,MAXSEC
   10 NRHO(I)=0
C-----INITIALIZE COUNT OF ENERGY NODES.
      NODES=0
C-----INITIALIZE INDEX TO NEXT RESONANCE TO READ.
      LHI=0
C
C     READ DATA FOR EACH ISOTOPE.
C
C-----INITIALIZE SUM OF FRACTIONAL ABUNDANCES.
      ABNSUM=0.0
      DO 320 IS=1,NIS
C-----SAVE CURRENT LFWX FLAG AND SET LFWX=0 FOR EACH ISOTOPE.
      IF(LFWX.gt.0) LFWSUM=LFWX
      LFWX=0
C-----DEFINE ISOTOPE ZA, ABUNDANCE, FLAG DEFINING WHETHER OR NOT FISSION
C-----WIDTHS ARE GIVEN AND NUMBER OF ENERGY RANGES.
      CALL CARDIO(ZAI,ABN,L1,LFW,NER,N2)
C-----INCREMENT SUM OF FRACTONAL ABUNDANCE.
      ABNSUM=ABNSUM+ABN
      IZAI=ZAI
      CALL ZAHOL(IZAI,ZABCDI)
      LFWO=LFW+1
      IF(LFW.NE.0.AND.LFW.NE.1) LFWO=3
      CALL OUT9(ABN,FIELD(1,1))
      WRITE(OUTP,460) IS,ZABCDI,(FIELD(M,1),M=1,11),LFW,FISLST(LFWO),NER
      WRITE(*   ,460) IS,ZABCDI,(FIELD(M,1),M=1,11),LFW,FISLST(LFWO),NER
C
C     PRINT WARNING IF SUM OF FRACTIONAL ABUNDANCES IS NOT CLOSE TO
C     1.0
C
      IF(IS.NE.NIS) GO TO 20
      IF(DABS(ABNSUM-1.0D+00).LE.0.001) GO TO 20
      WRITE(OUTP,420)
      IF(NIS.EQ.1) WRITE(OUTP,430) ABNSUM
      IF(NIS.GT.1) WRITE(OUTP,440) ABNSUM
      WRITE(OUTP,450)
      WRITE(*   ,420)
      IF(NIS.EQ.1) WRITE(*   ,430) ABNSUM
      IF(NIS.GT.1) WRITE(*   ,440) ABNSUM
      WRITE(*   ,450)
C
C     READ DATA FOR EACH ENERGY RANGE.
C
   20 DO 300 IE=1,NER
C-----DEFINE ENERGY RANGE, TYPE OF RESONANCE RANGE (NO PARAMETERS/
C-----RESOLVED/UNRESOLVED) AND TYPE OF PARAMETER.
      CALL CARDI(ELX,EHX,LRUIN,LRF,N1,N2)
c-----Save original formalism = it ma be changed for testing
      LRFIN = LRF
c***** DEBUG
c-----Activate to change LRF=3 to LRF=2
c-----This is the one and only line to change
c     if(LRF.eq.3) LRF = 2
c***** DEBUG
C
C     IF ENDS OF RESONANCE REGION ARE VERY CLOSE TO THE ENDS OF ANY
C     OTHER RESONANCE REGION INSURE THEY ARE EXACTLY THE SAME (THIS
C     PROCEDURE WILL AVOID MICRO OVERLAPS OR HOLES BETWEEN RESONANCE
C     REGIONS).
C
      IF(NSECT.LE.0) GO TO 40
      DO 30 MSECT=1,NSECT
C-----FIRST CHECK FOR SAME ENERGY REGION (COMPARE LOWER TO LOWER AND
C-----UPPER TO UPPER ENERGY LIMITS).
      IF(DABS(ELX-EL(MSECT)).LE.DEMIN*EL(MSECT)) ELX=EL(MSECT)
      IF(DABS(EHX-EH(MSECT)).LE.DEMIN*EH(MSECT)) EHX=EH(MSECT)
C-----NEXT CHECK FOR ADJOINING ENERGY REGIONS (COMPARE LOWER TO UPPER
C-----AND UPPER TO LOWER ENERGY LIMITS).
      IF(DABS(ELX-EH(MSECT)).LE.DEMIN*EH(MSECT)) ELX=EH(MSECT)
      IF(DABS(EHX-EL(MSECT)).LE.DEMIN*EL(MSECT)) EHX=EL(MSECT)
   30 CONTINUE
C-----SAVE FLAG TO DETERMINE WHETHER OR NOT SCATTERING RADIUS IS
C-----ENERGY DEPENDENT (ONLY LEGAL FOR BREIT-WIGNER PARAMETERS).
   40 NRO=N1
C-----SAVE FLAG TO DETERMINE HOW TO USE THE CHANNEL AND SCATTERING RADII
      NAPS=N2
C
C     DEFINE LRU FOR INTERNAL USE
C     LRUIN=0 -     NO RESONANCE PARAMETERS
C          =1 -     RESOLVED PARAMETERS.
C          =2 -     UNRESOLVED PARAMETERS.
C          =3 - 5 - SAME AS 0 - 2, ONLY THE RESONANCE CONTRIBUTION HAS
C                   ALREADY BEEN ADDED TO THE FILE 3 CROSS SECTIONS
C                   (USUALLY INDICATING THAT THIS MAT HAS ALREADY BEEN
C                   PROCESSED BY THIS OR A SIMILAR PROGRAM).
C     IF LRUIN IS 3 TO 5 THE PARAMETERS WILL BE READ BUT IGNORED IN ALL
C     SUBSEQUENT RESONANCE REGION CALCULATIONS.
C     SAVE LRU AS INPUT (LRUIN), DEFINE LRU FOR OUTPUT INDICATING THAT
C     THE RESONANCE CONTRIBUTION HAS ALREADY BEEN ADDED TO CROSS SECTION
C     (LRUOUT) AND DEFINE LRU FOR READING (LRU= 0, 1 OR 2).
C
C     PRINT ENERGY LIMITS AND TYPE OF RESONANCE REGION (UNRESOLVED OR
C     UNRESOLVED).
C
      LRUOUT=LRUIN
C-----DO NOT CHANGE PARAMETER TYPE FOR ENDF/B-VI DATA (FILE 1 IS USED
C-----TO DEFINE THAT MATERIAL HAS BEEN PROCESSED).
      IF(IVERSE.EQ.6) GO TO 50
      IF(LRUOUT.GE.0.AND.LRUOUT.LE.2) LRUOUT=LRUOUT+3
   50 LRU=LRUIN
      IF(LRU.GT.2) LRU=LRU-3
C-----DEFINE AND PRINT INTERPRETATION OF TYPE OF RESONANCE REGION.
      LRULST=LRUIN+1
      IF(LRUIN.GE.3.AND.LRUIN.LE.5) LRULST=LRULST-3
      IF(LRUIN.LT.0.OR.LRUIN.GT.5) LRUIN=4
      CALL OUT9(ELX,FIELD(1,1))
      CALL OUT9(EHX,FIELD(1,2))
      WRITE(OUTP,470) ((FIELD(M,I),M=1,11),I=1,2),LRUIN,REGLST(LRULST)
      WRITE(*   ,470) ((FIELD(M,I),M=1,11),I=1,2),LRUIN,REGLST(LRULST)
C
C     DEFINE AND PRINT INTERPRETATION OF TYPE OF RESOLVED PARAMETERS.
C
      IF(LRU.NE.1) GO TO 60
      LRFLST=LRF
      IF(LRF.LT.1.OR.LRF.GT.7) LRFLST=8
      WRITE(OUTP,480) LRF,RTYPE(LRFLST)
      WRITE(*   ,480) LRF,RTYPE(LRFLST)
      GO TO 70
C
C     DEFINE AND PRINT INTERPRETATION OF TYPE OF UNRESOLVED PARAMETERS.
C
   60 IF(LRU.NE.2) GO TO 70
      LRFLST=4
      IF(LRF.EQ.1.AND.LFW.EQ.0) LRFLST=1
      IF(LRF.EQ.1.AND.LFW.EQ.1) LRFLST=2
      IF(LRF.EQ.2) LRFLST=3
      WRITE(OUTP,490) LRF,UTYPE(LRFLST)
      WRITE(*   ,490) LRF,UTYPE(LRFLST)
C
C     PRINT INTERPRETATION OF NRO (ENERGY INDEPENDENT OR DEPENDENT
C     SCATTERING RADIUS) AND NAPS (CALCULATE CHANNEL RADIUS OR DEFINE
C     IT AS EQUAL TO THE SCATTERING RADIUS).
C
C-----DEFINE WHETHER OR SCATTERING RADIUS IS ENERGY DEPENDENT.
   70 KNRO=NRO+1
      IF(NRO.LT.0.OR.NRO.GT.1) KNRO=3
      WRITE(OUTP,500) NRO,HOLNRO(KNRO)
      WRITE(*   ,500) NRO,HOLNRO(KNRO)
C-----DEFINE WHETHER CHANNEL RADIUS WILL BE CALCULATED OR SET EQUAL TO
C-----THE SCATTERING RADIUS.
      KNAPS=NAPS+1
      IF(KNAPS.EQ.3) KNAPS=2
      IF(NAPS.LT.0.OR.NAPS.GT.2) KNAPS=3
      WRITE(OUTP,510) NAPS,HOLNAP(KNAPS)
      WRITE(*   ,510) NAPS,HOLNAP(KNAPS)
C
C     PRINT WARNING IF RESONANCE CONTRIBUTION HAS ALREADY BEEN ADDED TO
C     BACKGROUND CROSS SECTIONS.
C
      IF(LRUIN.GE.3.AND.LRUIN.LE.5) WRITE(OUTP,520) LRUIN
      IF(LRUIN.GE.3.AND.LRUIN.LE.5) WRITE(*   ,520) LRUIN
C
C     ERROR IF RESONANCE REGION BOUNDARIES ARE NOT IN ASCENDING ENERGY
C     ORDER OR NOT IN EXPECTED ENERGY RANGE.
C
C-----IF RESONANCE REGION BOUNDARIES ARE NOT IN ASCENDING ENERGY ORDER
C-----PRINT MESSAGE AND TERMINATE.
      IF(ELX.GE.EHX) GO TO 80
C-----IF RESONANCE REGION BOUNDARIES ARE NOT IN EXPECTED ENERGY RANGE
C-----PRINT MESSAGE AND TRUNCATE RANGE TO EXPECTED RANGE.
      IF(ELX.GE.EMINL.AND.ELX.LE.EMAX.AND.EHX.GE.EMINL.AND.EHX.LE.EMAX)
     1 GO TO 100
      IF(ELX.LT.EMIN) ELX=EMIN
      IF(ELX.GT.EMAX) ELX=EMAX
      IF(EHX.LT.EMIN) EHX=EMIN
      IF(EHX.GT.EMAX) EHX=EMAX
      CALL OUT9(EMIN,FIELD(1,1))
      CALL OUT9(EMAX,FIELD(1,2))
      WRITE(OUTP,600) ((FIELD(M,I),M=1,11),I=1,2)
      WRITE(*   ,600) ((FIELD(M,I),M=1,11),I=1,2)
      CALL OUT9(ELX,FIELD(1,1))
      CALL OUT9(EHX,FIELD(1,2))
      WRITE(OUTP,610) ((FIELD(M,I),M=1,11),I=1,2)
      WRITE(*   ,610) ((FIELD(M,I),M=1,11),I=1,2)
      IF(ELX.LT.EHX) GO TO 100
   80 WRITE(OUTP,90)
      WRITE(*   ,90)
   90 FORMAT(///' ERROR - Resonance Region Energy Limits MUST be in'/
     1       '            Ascending Energy Order.'/
     2       '            Execution Terminated.'///)
      CALL ENDERROR
C
C     ERROR STOP IF LRU, LRF OR LRF, LFW COMBINATION IS ILLEGAL.
C     UNLESS THESE PARAMETERS ARE CORRECT THE PROGRAM CANNOT DETERMINE
C     THE FORMAT OF THE ENDF/B DATA.
C
  100 CALL CARDO(ELX,EHX,LRUOUT,LRF,N1,N2)
C-----IF NECESSARY READ ENERGY DEPENDENT SCATTERING RADIUS.
C-----06/24/05 - MOVED CALL REAP AFTER ABOVE OUTPUT LINE
      IF(NRO.EQ.1) CALL RDAP
C-----TERMINATE IF ILLEGAL LRU.
      IF(LRUIN.GE.0.AND.LRUIN.LE.5) GO TO 120
      WRITE(OUTP,110) LRUIN
      WRITE(*   ,110) LRUIN
  110 FORMAT(///' ERROR - Illegal LRU=',I5,' (Expect 0 to 5).'/
     2          '         Cannot Determine Format of ENDF/B Data.'/
     3          '         Execution Terminated.'///)
      CALL ENDERROR
C-----COPY SECTION WITH NO PARAMETERS.
  120 IF(LRU.EQ.0) GO TO 180
C
C     TERMINATE IF ILLEGAL RESOLVED REGION LRF.
C
      IF(LRU.NE.1) GO TO 140
C-----02/14/04 - CHANGE TO 7 FOR NEW REIC H=MOORE
      IF(LRF.GE.1.AND.LRF.LE.7) GO TO 160
      WRITE(OUTP,130) LRF
      WRITE(*   ,130) LRF
  130 FORMAT(///' ERROR - Illegal LRF=',I5,' (Expect 1 to 7).'/
     1          '         Cannot Determine Format of ENDF/B Data'/
     2          '         Execution Terminated'///)
      CALL ENDERROR
C
C     TERMINATE IF ILLEGAL UNRESOLVED LRF, LFW COMBINATION.
C
  140 IF(LRU.NE.2) GO TO 160
      LRFLST=4
      IF(LRF.EQ.1.AND.LFW.EQ.0) LRFLST=1
      IF(LRF.EQ.1.AND.LFW.EQ.1) LRFLST=2
      IF(LRF.EQ.2) LRFLST=3
      IF(LRFLST.NE.4) GO TO 160
      WRITE(OUTP,150) LRF,LFW
      WRITE(*   ,150) LRF,LFW
  150 FORMAT(///' ERROR - Illegal LRF=',I5,' LFW=',I5,' Combination'/
     2          '         (Expect LRF=1,LFW=0 OR LRF=1,LFW=1 OR LRF=2)'/
     3          '         Cannot Determine Format of ENDF/B Data'/
     4          '         Execution Terminated'///)
      CALL ENDERROR
C
C     ERROR IF ILLEGAL VALUE FOR NRO (ENERGY INDEPENDENT OR DEPENDENT
C     SCATTERING RADIUS) AND/OR NAPS (CALCULATE CHANNEL RADIUS OR SET
C     IT EQUAL TO THE SCATTERING RADIUS).
C
  160 IF(KNRO.NE.3) GO TO 170
C-----ILLEGAL NRO VALUE.
      WRITE(OUTP,580) NRO
      WRITE(*   ,580) NRO
      NRO=0
  170 IF(KNAPS.NE.3) GO TO 190
C-----ILLEGAL NAPS VALUE.
      WRITE(OUTP,590) NAPS
      WRITE(*   ,590) NAPS
      NAPS=0
      GO TO 190
C
C     READ AND SAVE ALL PARAMETERS.
C
C-----COPY SECTION WITH NO PARAMETERS.
  180 CALL CARDIO(SPI,AP,L1,L2,N1,N2)
      CALL OUT9(SPI,FIELD(1,1))
      CALL OUT9(AP ,FIELD(1,2))
      WRITE(OUTP,530) ((FIELD(M,I),M=1,2),I=1,2)
      WRITE(*   ,530) ((FIELD(M,I),M=1,2),I=1,2)
      GO TO 300
C-----SAVE INITIAL RESONANCE TABLE PARAMETERS. IF RESONANCE CONTRIBUTION
C-----FOR THIS SECTION HAS ALREADY BEEN ADDED TO THE FILE 3 CROSS
C-----SECTION THEY WILL BE RESTORED AFTER THE DATA HAS BEEN READ
C-----(THIS WILL EFFECTIVELY IGNORE THE SECTION).
  190 NSECTI=NSECT
      NODESI=NODES
      LHII=LHI
C-----ARE PARAMETERS RESOLVED OR UNRESOLVED.
      IF(LRU.EQ.2) GO TO 260
C-----READ ALL RESOLVED RESONANCE PARAMETERS.
      GO TO (200,200,210,220,230,240,250),LRF
  200 CALL RDBW
      GO TO 270
C-----FOR ALL VERSIONS OF ENDF/B FORMAT ASSUME GENERAL REICH-MOORE
C-----TREATMENT WITH TWO FISSION CHANNELS (NOTE, ENDF/B-VI REICH-MOORE
C-----FORMAT HAS BEEN UPDATED TO BE EXACTLY THE SAME AS THE FORMAT IN
C-----EARLIER VERSIONS OF ENDF/B).
  210 CALL RDRM1
      GO TO 270
  220 CALL RDAA
      GO TO 270
  230 CALL RDGRM
      GO TO 270
  240 CALL RDHRF
      GO TO 270
C-----NEW (2003) REICH-MOORE WITH COMPETITION
  250 CALL RDRML
      GO TO 270
C-----READ UNRESOLVED PARAMETERS (ANY OF THE 3 REPRESENTATIONS).
  260 CALL RDUR
C
C     IF RESONANCE CONTRIBUTION OF CURRENT SECTION ALREADY BEEN ADDED
C     TO CROSS SECTIONS, SKIP THIS SECTION.
C
  270 IF(LRUIN.NE.1.AND.LRUIN.NE.2) GO TO 280
      GO TO 290
C-----RESTORE INITIAL RESONANCE TABLE PARAMETERS (I.E., IGNOR CURRENT
C-----SECTION).
  280 NSECT=NSECTI
      NODES=NODESI
      LHI=LHII
      GO TO 300
C-----NO. USE LOWER AND UPPER ENERGY LIMIT AS NODES (USE ZERO WIDTH).
  290 CALL NOODLE(ELX,ZERO,EMIN,EMAX)
      CALL NOODLE(EHX,ZERO,EMIN,EMAX)
C-----END OF ENERGY RANGE LOOP.
  300 CONTINUE
C-----END OF ISOTOPE LOOP. TEST CONSISTENCY OF LFW AND LFI FLAGS AND
C-----FISSION DATA READ.
      IF(LFW.EQ.0.AND.LFWX.GT.0) WRITE(OUTP,540) LFW
      IF(LFW.GT.0.AND.LFWX.EQ.0) WRITE(OUTP,550) LFW
      IF(LFW.EQ.0.AND.LFWX.GT.0) WRITE(*   ,540) LFW
      IF(LFW.GT.0.AND.LFWX.EQ.0) WRITE(*   ,550) LFW
C-----ONLY CHECK LFI IF SECTION MF=1, MT-451 IS PRESENT (SINCE LFI IS
C-----DEFINED IN MF=1, MT=451).
      IF(MT451.LE.0) GO TO 310
      IF(LFI.EQ.0.AND.LFWX.GT.0) WRITE(OUTP,560) LFI
      IF(LFI.GT.0.AND.LFWX.EQ.0) WRITE(OUTP,570) LFI
      IF(LFI.EQ.0.AND.LFWX.GT.0) WRITE(*   ,560) LFI
      IF(LFI.GT.0.AND.LFWX.EQ.0) WRITE(*   ,570) LFI
C-----DEFINE CUMULATIVE LFWX FLAG.
  310 if(LFWX.gt.0) LFWSUM = 1
      IF(LFWSUM.GT.0.OR.LFWX.GT.0) LFWX=1
  320 CONTINUE
C-----ALL PARAMETERS READ. COPY TO END OF FILE 2.
      CALL COPYF
C
C     ALL RESONANCE PARAMETERS HAVE BEEN READ. IF THERE ARE NO SECTIONS
C     WHOSE RESONANCE CONTRIBUTION MUST BE ADDED TO THE CROSS SECTIONS
C     THERE IS NOTHING TO DO. OTHERWISE, (1) DEFINE TABLE OF ENERGY
C     RANGES (USUALLY 1 OR 2, I.E. RESOLVED AND/OR UNRESOLVED), (2)
C     DEFINE WHICH NODES ARE WITHIN EACH ENERGY RANGE.
C
      IF(NSECT.LE.0) GO TO 400
C
C     DEFINE ENERGY RANGES.
C
C-----DEFINE TABLE OF ENDS OF ALL RESONANCE REGIONS.
      IRANGE=0
      DO 350 IS=1,NSECT
      ELTNOW=EL(IS)
      EHTNOW=EH(IS)
      IF(ELTNOW.LT.EMIN) ELTNOW=EMIN
      IF(EHTNOW.GT.EMAX) EHTNOW=EMAX
      IF(ELTNOW.GE.EHTNOW) GO TO 350
C-----IGNOR REPEATED ENERGY RANGES.
      IF(IRANGE.LT.2) GO TO 340
      DO 330 I=1,IRANGE,2
      IF(ELTNOW.EQ.ERANGE(I).AND.EHTNOW.EQ.ERANGE(I+1)) GO TO 350
  330 CONTINUE
C-----SAVE BOTH ENDS OF ENERGY REGION.
  340 IRANGE=IRANGE+1
      ERANGE(IRANGE)=ELTNOW
      IRANGE=IRANGE+1
      ERANGE(IRANGE)=EHTNOW
  350 CONTINUE
C-----SORT RESONANCE REGION BOUNDARIES INTO ASCENDING SORTS.
      CALL SORTD(ERANGE,IRANGE)
C-----ELIMINATE DUPLICATE ENERGY LIMITS.
      NRANGE=1
      DO 360 IS=2,IRANGE
      IF(DABS(ERANGE(IS)-ERANGE(NRANGE)).LE.DABS(DEMIN*ERANGE(NRANGE)))
     1 GO TO 360
      NRANGE=NRANGE+1
      ERANGE(NRANGE)=ERANGE(IS)
  360 CONTINUE
C-----ADD THERMAL ENERGY (0.0253 EV) IF RESONANCE REGION SPANS THERMAL
C-----ENERGY.
      CALL NOODLE(ETHERM,ZERO,ERANGE(1),ERANGE(NRANGE))
C
C     DEFINE WHICH NODES CONTRIBUTE TO EACH RESONANCE REGION.
C
      IN1=1
      IN2=1
      DO 390 IRANGE=2,NRANGE
      IF(IN1.GE.NODES) GO TO 380
      DO 370 IN2=IN1,NODES
      IF(ENODE(IN2).ge.ERANGE(IRANGE)) go to 380
  370 CONTINUE
      IN2=NODES
  380 NODLOW(IRANGE)=IN1
      NODHI(IRANGE)=IN2
  390 IN1=IN2
  400 RETURN
  410 FORMAT(1X,78('=')/' Listing of All Resonance Parameters'/
     1 1X,78('=')/
     1       ' Element or Material------------------',10A1/
     2       ' Atomic Weight Ratio------------------',11A1/
     3       ' Number of Isotopes-------------------',I11)
  420 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' For ENDF/B Evaluations the Fractional Abundance is Defined'/
     2 ' to be the Fraction of Each Isotope in the Evaluation.'/
     3 ' It is NOT the Naturally Occurring Isotopic Abundance'/
     4 ' Unless the Evaluation Contains Data for ALL Isotopes of'/
     5 ' an Element. For Any Evaluation the Sum of the Fractional'/
     6 ' Abundances for All Isotopes Included in the Evaluation'/
     7 ' Should be 1.0.')
  430 FORMAT(' For this Evaluation there is Only One Isotope'/
     1 ' and the Fractional Abundance is',F10.5,' Rather Than 1.0.')
  440 FORMAT(' For this Evaluation the Sum of the Fractional'/
     1 ' Abundances is',F10.5,' Rather than 1.0.')
  450 FORMAT(
     1 ' Either You are Doing this on Purpose for Some Special'/
     2 ' Application or Else the ENDF/B Data is Incorrect.'/
     3 ' This Program will Continue with the Calculations Asuuming'/
     4 ' that the Abundances Read from the Evaluation are Correct.'/
     5 ' Check ENDF/B Data and Correct the Fractional Abundances.')
  460 FORMAT(1X,78('=')/
     1       ' Isotope Number-----------------------',I11/
     2       ' Isotope------------------------------',10A1/
     3       ' Fractional Abundance-----------------',11A1/
     4       ' LFW (Fission Widths)-----------------',I11,A28/
     4       ' Number of Energy Ranges--------------',I11)
  470 FORMAT(1X,78('=')/
     1 ' Lower Limit of the Energy Range------',11A1,' eV'/
     2 ' Upper Limit of the Energy Range------',11A1,' eV'/
     3 ' LRU (Type of Region)-----------------',I11,A28)
  480 FORMAT(' LRF (Type of Resolved Parameters)----',I11,A28)
  490 FORMAT(' LRF (Energy Dependence of Widths)----',I11,A28)
  500 FORMAT(' NRO (Scattering Radius)--------------',I11,1X,A28)
  510 FORMAT(' NAPS (Channel Radius)----------------',I11,1X,A28)
  520 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' Note, LRU=',I2,' Indicates that the Resonance Contribution'/
     2 ' has Already been Added to the Cross Sections. This Section'/
     3 ' will be Read, but Ignored in All Resonance Calculations.')
  530 FORMAT(1X,78('=')/' No Parameters'/1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Effective Scattering Radius (A+)-----',11A1)
  540 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' WARNING - LFW=',I2,' (Indicates NO Fission Widths Given).'/
     1 '           but Fission Data Read are NOT Zero. Will Ignor'/
     2 '           LFW and Assume Data Read is Correct.'/1X,78('='))
  550 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' WARNING - LFW=',I2,' (Indicates Fission Widths Given).'/
     1 '           but Fission Data Read are ALL Zero. Will Ignor'/
     2 '           LFW and Assume Data Read is Correct.'/1X,78('='))
  560 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' WARNING - LFI=',I2,' from MF=1, MT=451 (Indicates Material',
     2 ' is NOT Fissile)'/
     1 '           but Fission Data Read are NOT Zero. Will Ignor',
     2 ' LFI and Assume'/
     2 '           Data Read is Correct.'/1X,78('='))
  570 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' WARNING - LFI=',I2,' from MF=1, MT=451 (Indicates Material',
     2 ' is Fissile)'/
     1 '           but Fission Data Read are ALL Zero. Will Ignor',
     2 ' LFI and Assume'/
     2 '           Data Read is Correct.'/1X,78('='))
  580 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' WARNING - NRO=',I5,' (Illegal Value)'/
     2 10X,' in Order to Continue this Program will Set NRO=0'/
     3 ' WARNING - This May Result in Errors - Check ENDF/B Data.')
  590 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' WARNING - NAPS=',I5,' (Illegal Value)'/
     2 10X,' In Order to Continue this Program will Set NAPS=0'/
     3 ' WARNING - This May Result in Errors - Check ENDF/B Data.')
  600 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 ' Program Expects Resonance Region Between ',11A1,' and',
     2 11A1,' eV.')
  610 FORMAT(
     3 ' Program Will Truncate Resonance Region to',11A1,' and',
     2 11A1,' eV.'/
     4 ' Check Evaluated Data.'/1X,78('='))
      END
CAK   SUBROUTINE FILE3
      SUBROUTINE FILE3r
C=======================================================================
C
C     OUTPUT DATA IN ENDF/B FORMAT. ALL SECTIONS OF FILE 3 DATA WILL BE
C     PROCESSED BY THIS SUBROUTINE. ALL SECTIONS EXCEPT TOTAL, ELASTIC,
C     FISSION AND CAPTURE WILL BE COPIED BY THIS ROUTINE. FOR TOTAL,
C     ELASTIC, FISSION AND CAPTURE IF THERE ARE BACKGROUND CROSS SECTION
C     THE RESONANCE AND BACKGROUND CROSS SECTIONS WILL BE COMBINED AND
C     OUTPUT. IF THERE IS NO BACKGROUND CROSS SECTION ONLY THE RESONANCE
C     CONTRIBUTION WILL BE OUTPUT. IN EITHER CASE THE DATA (RESONANCE OR
C     RESONANCE PLUS BACKGROUND CONTRIBUTIONS) WILL NOT BE THINNED, I.E.
C     THE OUTPUT ENERGY GRID WILL BE THE UNION OF THE RESONANCE AND
C     BACKGROUND ENERGY GRIDS.
C
C     AFTER ALL OF FILE 3 HAS BEEN PROCESSED THE REMAINDER OF THE MAT
C     WILL BE COPIED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      REAL*4 SECONDS
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*8 REACTS,RSORT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/UNITS/ISCR2,ISCR23
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/POINTN/NPOINT,KPOINT
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/EDGES/ERANGE(50),NODLOW(50),NODHI(50),NRANGE,IRANGE
      COMMON/PAGER/NPAGE,NPAGP1,NPAGM1
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/LASTE/ELAST
      COMMON/FLAGS/MINUS3,IMPLUS
      COMMON/WHATZA/IZANOW,MATNOW,TEMP3,IVERSE,INT45
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
      COMMON/NEG1COM/NEG1    ! 2014/11/02 - add to FILE2 and FILE3
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
c-----07/26/09 - ADDED Q-VALUE FOR COMPETITIVE MF=3 OUTPUT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION MREACT(20),MTADD(20),REACTS(20),QADD(11),NEGADD(11),
     1 NBTO(1),INTO(1),IHEAD(7),ISAVE1(7),SAVE1(2),MTCP(103:107)
      CHARACTER*8 BCDCP(103:107)
C-----3/29/07 - ADDED MISSING EQUIVALENCE
      EQUIVALENCE (L1H,IHEAD(1))
C     CHARGED.. 103 104 105 106 107
      DATA BCDCP/
     1 '   (n,p)',
     1 '   (n,d)',
     1 '   (n,t)',
     1 ' (n,He3)',
     1 '   (n,a)'/
c       12345678
      DATA MTCP/600,650,700,750,800/
      DATA MREACT/1,2, 4, 4,  3,15*0/
      DATA MTADD/ 1,2,18,19,102,15*0/
      DATA REACTS/
     1 '   Total',
     2 ' Elastic',
     3 ' Fission',
     4 ' (n,f) 1',
     5 ' Capture',
     6 15*' Others '/
      DATA INTO/2/
      DATA IZERO/0/
      DATA IUSE23/0/
      DATA ISAVE1/0,0,0,0,0,0,0/
      DATA SAVE1/0.0D+00,0.0D+00/
      DATA OKDIFF/1.0D-09/
      DATA IONE/1/
      DATA ZEROD/0.0D+00/
c-----------------------------------------------------------------------
C
C     DEFINE REACTIONS TO OUTPUT.
C
c-----------------------------------------------------------------------
C-----SET UP INITIALLY FOR SIMPLE CASE = TOTAL, ELASTIC, CAPTURE
      MREACT(1) = 1     ! total
      MREACT(2) = 2     ! elastic
      MREACT(3) = 3     ! capture
      NEGADD(1) = NEGTAB(1)
      NEGADD(2) = NEGTAB(2)
      NEGADD(3) = NEGTAB(3)
      MTADD (1) = 1
      MTADD (2) = 2
      MTADD (3) = 102
      REACTS(1) = '   Total'
      REACTS(2) = ' Elastic'
      REACTS(3) = ' Capture'
      QADD  (1) = 0.0D+0
      QADD  (2) = 0.0D+0
      QADD  (3) = 0.0D+0
      NROUT     = 3
C
      IF(LRF7.ne.0) GO TO 10  ! LRF7 = Counter of LRF=7 sections
C
C-----SIMPLE CASE - IF NECESSARY ADD 2 FISSION CHANNELS (MT=18 AND 19)
      IF(NREACT.EQ.3) GO TO 70
      MREACT(3) = 4     ! (n,f)
      MREACT(4) = 4     ! (n,f) first chance
      MREACT(5) = 3     ! capture
      NEGADD(3) = NEGTAB(4)
      NEGADD(4) = NEGTAB(4)
      NEGADD(5) = NEGTAB(3)
      MTADD (3) = 18
      MTADD (4) = 19
      MTADD (5) = 102
      REACTS(3) = ' Fission'
      REACTS(4) = ' (n,f) 1'
      REACTS(5) = ' Capture'
      QADD  (3) = 0.0D+0
      QADD  (4) = 0.0D+0
      QADD  (5) = 0.0D+0
      NROUT     = 5
      GO TO 70
C
C     GENERAL LRF=7 CASE
C
   10 NROUT  = NREACT
      KREACT = NREACT
C-----ASSUME 1 AND 2 ARE TOTAL AND ELASTIC - DEFINE ALL OTHERS.
      DO 40 J=3,NREACT
      MREACT(J) = J
      NEGADD(J) = NEGTAB (J)
      MTADD (J) = MTREACT(J)
      REACTS(J) = '   Other'
      QADD  (J) = QREACT (J)
      IF(MTADD(J).EQ.18) GO TO 20
      DO K=103,107
      IF(MTADD(J).EQ.K) GO TO 30
      ENDDO
      if(MTADD(J).eq.102) REACTS(J) = ' Capture'
      if(MTADD(J).eq. 19) REACTS(J) = ' (n,f) 1'
      if(MTADD(J).ge. 50.and.
     1   MTADD(J).le. 91) REACTS(J) = '   Level'
      GO TO 40
c-----MT=18 - add MT=19, second chance fission
   20 REACTS(J) = ' Fission'
      KREACT = KREACT + 1
      MREACT(KREACT) = MREACT(J)
      NEGADD(KREACT) = NEGTAB(J)
      MTADD (KREACT) = 19
      REACTS(KREACT) = ' (n,f) 1'
      QADD  (KREACT) = QREACT(J)
      GO TO 40
c-----MT = 103 rhrough 107 - add charged particle level
   30 REACTS(J) = BCDCP(K)
      KREACT = KREACT + 1
      MREACT(KREACT) = MREACT(J)
      NEGADD(KREACT) = NEGTAB(J)
      MTADD (KREACT) = MTCP(K)
      REACTS(KREACT) = '   Level'
      QADD  (KREACT) = QREACT(J)
   40 CONTINUE
      NREACT = KREACT
      NROUT  = KREACT
c
C     SORT INTO MT ORDER (Total and elastic already in order).
c
      DO 60 J =  3,NREACT
      DO 50 JJ=J+1,NREACT
      IF(MTADD(JJ).GE.MTADD(J)) GO TO 50
      MM1   = MTADD (JJ)
      MM2   = MREACT(JJ)
      MM3   = NEGADD(JJ)
      RSORT = REACTS(JJ)
      QSORT = QADD  (JJ)
      MTADD (JJ) = MTADD (J)
      NEGADD(JJ) = NEGADD(J)
      MREACT(JJ) = MREACT(J)
      REACTS(JJ) = REACTS(J)
      QADD  (JJ) = QADD  (J)
      MTADD (J)  = MM1
      MREACT(J)  = MM2
      NEGADD(J)  = MM3
      REACTS(J)  = RSORT
      QADD  (J)  = QSORT
   50 CONTINUE
   60 CONTINUE
C
c-----Print summary.
   70 write(*,80) (MTADD(kk),kk=1,NROUT)
      write(3,80) (MTADD(kk),kk=1,NROUT)
   80 format(' Summary of Output MT #'/1x,78('-')/(10i5))
C
C-----IF NO POINTS WERE GENERATED IN THE RESONANCE CALCULATION MERELY
C-----COPY REMAINDER OF MAT UNLESS IN EDIT MODE.
      IF(IMEDIT.EQ.2) GO TO 90
      IF(NPOINT.LE.0) GO TO 820
      WRITE(OUTP,840)
      WRITE(*   ,840)
      GO TO 100
   90 WRITE(OUTP,850)
      WRITE(*   ,850)
C-----INITIALIZE END OF FILE 3 FLAG OFF.
  100 MF3END=0
c-----------------------------------------------------------------------
c
c     Output results
c
c-----------------------------------------------------------------------
c-----initialize to ALL cross sections negative
      NEG1 = 0
      do i=1,NROUT
      NEGTAB(i) = 0
      enddo
C-----LOOP OVER RESONANCE REACTIONS.
      DO 750 IREACT=1,NROUT
C-----INITIALIZE TO OUTPUT ALL POINTS (MAY BE LESS FOR COMPETITION)
      NOUT2 = NPOINT
      MSKIP = 0
C-----DEFINE INDEX TO CURRENT REACTION.
      INDEX=MREACT(IREACT)
      IF(INDEX.LE.0) GO TO 750
      IMESS=1
C-----INITIALIZE NEGATIVE CROSS SECTION POINT COUNT.
      MINUS3=0
      IMPLUS=0
C-----INITIALIZE LAST ENERGY READ FROM FILE 3 FOR ASCENDING ENERGY TEST.
      ELAST=0.0
C-----TEST FOR CURRENT STATUS
C-----MF3END=-1 - STILL IN FILE 3. RESTORE LAST LINE READ.
C-----      = 0 - STILL IN FILE 3. READ NEXT LINE.
C-----      = 1 - END OF FILE 3 READ (NO BACKGROUND)
C-----      = 2 - PAST END OF FILE 3 READ (NO BACKGROUND)
      IF(MF3END.eq.0) go to 120
      IF(MF3END.gt.0) go to 640
C-----RESTORE LAST LINE READ.
      C1H=SAVE1(1)
      C2H=SAVE1(2)
      DO 110 I=1,7
  110 IHEAD(I)=ISAVE1(I)
      MF3END=0
      GO TO 150
C-----READ FIRST LINE OF NEXT SECTION OF FILE 3, FEND OR MEND LINE.
  120 CALL CONTI
      IF(MFH.eq.3) go to 150
      IF(MFH.gt.3) go to 130
C-----END OF FILE 3. SET FLAG. IF MEND LINE TREAT AS PAST END OF FILE 3.
      IF(MATH.LE.0) GO TO 130
      MF3END=1
      GO TO 640
C-----PAST THE END OF FILE 3. SET FLAG AND SAVE LINE.
  130 MF3END=2
      SAVE1(1)=C1H
      SAVE1(2)=C2H
      DO 140 I=1,7
  140 ISAVE1(I)=IHEAD(I)
      GO TO 640
C-----------------------------------------------------------------------
C
C     STILL IN FILE 3. SEE IF RESONANCE CONTRIBUTION MUST BE ADDED TO
C     CURRENT CROSS SECTION. DEPENDING ON CURRENT MT FROM FILE 3
C     EITHER,
C     (1) COPY SECTION - NOT UP TO REQUIRED SECTION YET.
C     (2) COMBINE RESONANCE AND BACKGROUND CONTRBUTION -SECTION MATCH.
C     (3) OUTPUT RESONANCE CONTRIBUTION - NO BACKGROUND.
C
C-----------------------------------------------------------------------
  150 CONTINUE
      IF(MTH.eq.MTADD(IREACT)) go to 160
      IF(MTH.gt.MTADD(IREACT)) go to 620
C-----COPY SECTION.
      CALL CONTO
      CALL COPYS
      GO TO 120
C-----------------------------------------------------------------------
C
C     BACKGROUND DATA PRESENT.
C
C-----------------------------------------------------------------------
C-----READ SECTION LEADER LINE AND INTERPOLATION LAW.
  160 CALL CARDI(C1,C2,L1,L2,N1,N2)
      CALL TERPI(NBTF(1),INTF(1),N1)
C-----DEFINE BACKGROUND TEMPERATURE.
      IF(IVERSE.NE.6) GO TO 170
      TEMP=TEMP3
      GO TO 180
  170 TEMP=C1
      IF(MTH.EQ.1) TEMP1=TEMP
      IF(L2.NE.0) TEMP=TEMP1
C-----IF BACKGROUND TEMPERATURE IS NOT 0 KELVIN PRINT WARNING MESSAGE
C-----AND DO NOT COMBINE RESONANCE AND BACKGROUND CONTRIBUTIONS (SKIP
C-----BACKGROUND AND OUTPUT ONLY THE RESONANCE CONTRIBUTION).
  180 IF(TEMP.LE.0.0) GO TO 190
      IMESS=2
      CALL SKIPS
      GO TO 670
C-----CHECK INTERPOLATION LAW.
  190 DO 200 I=1,N1
      IF(INTF(I).NE.2) GO TO 210
  200 CONTINUE
      GO TO 220
C-----INTERPLATION LAW IS NOT LINEAR-LINEAR. PRINT WARNING MESSAGE
C-----AND DO NOT COMBINE RESONANCE AND BACKGROUND CONTRIBUTIONS (SKIP
C-----BACKGROUND AND OUTPUT ONLY THE RESONANCE CONTRIBUTION).
  210 IMESS=3
      CALL SKIPS
      GO TO 670
C-----SECTION HEADER AND LEADER CARDS HAVE BEEN READ. TEMPERATURE AND
C-----INTERPOLATION LAW CHECKED. IF IN EDIT MOODE SKIP SECTION AND
C-----CONTINUE.
  220 IF(IMEDIT.NE.2) GO TO 230
      WRITE(OUTP,880) REACTS(IREACT),MTADD(IREACT),N2
      WRITE(*   ,880) REACTS(IREACT),MTADD(IREACT),N2
      CALL SKIPS
      GO TO 750
C-----IF ALL FISSION WIDTHS ARE ZERO COPY BACKGROUND (NO RESONANCE
C-----CONTRIBUTION...BOTH FISSION AND FIRST CHANCE FISSION).
  230 IF(LFWX.NE.0) GO TO 240
      IF(MTH.NE.18.AND.MTH.NE.19) GO TO 240
      WRITE(OUTP,970) REACTS(IREACT),MTADD(IREACT),IZERO,N2,N2
      WRITE(*   ,970) REACTS(IREACT),MTADD(IREACT),IZERO,N2,N2
C-----OUTPUT SECTION HEADER AND LEADER CARDS AND INTERPOLATION LAW.
      CALL CONTO
      CALL CARDO(C1,C2,L1,L2,IONE,N2)
      NBTO(1)=N2
      CALL TERPO(NBTO,INTO,1)
C-----COPY REMAINDER OF SECTION.
      CALL COPYS
      GO TO 750
C-----------------------------------------------------------------------
C
C     COMBINE RESONANCE AND BACKGROUND CONTRIBUTIONS.
C
C-----------------------------------------------------------------------
C-----DEFINE RESONANCE AND TABULATED POINT COUNTS.
  240 NPT2=NPOINT
      NPT3=N2
      LEFT2=NPOINT
      LEFT3=N2
C-----INITIALIZE COMBINED POINT COUNTS.
      IPT23=0
      NPT23=0
C-----INITIALIZE NON-ZERO CROSS SECTION FLAG.
      KMPLUS=0
C-----INITIALIZE RESONANCE AND TABULATED POINT INDICES.
      IPT2=1
      IPT3=1
      KPT2=IPT2
      KPT3=IPT3
C-----IF REQUIRED POSITION FILE 2 SCRATCH FILE FOR READING.
      IF(NPOINT.GT.NPAGE) REWIND ISCR2
C-----READ FIRST PAGE OF TABULATED DATA.
      IF(NPT3.GT.NPAGE) NPT3=NPAGE
      CALL POINTI(ETAB3,SIG3,NPT3)
      LEFT3=LEFT3-NPT3
C-----IF ENTIRE BACKGROUND FITS IN CORE AND IT IS ZERO AT ALL ENERGIES
C-----SIMPLY ADD END POINT ENERGIES AND OUTPUT (NO NEED TO PERFORM
C-----ADDITION).
      IF(LEFT3.GT.0) GO TO 320
      DO 250 I=1,NPT3
      IF(DABS(SIG3(I)).NE.0.0D+00) GO TO 320
  250 CONTINUE
C-----------------------------------------------------------------------
C
C     BACKGROUND IS ZERO. ADDITION NOT REQUIRED.
C
C-----------------------------------------------------------------------
C-----IF REQUIRED ADD END POINTS.
      NPT23=NPOINT
      E3LST=ETAB3(1)
      E3=ETAB3(NPT3)
      IF(E3LST.LT.ERANGE(1)) NPT23=NPT23+1
      IF(E3.GT.ERANGE(NRANGE)) NPT23=NPT23+1
      WRITE(OUTP,890) REACTS(IREACT),MTADD(IREACT),NPOINT,N2,NPT23
      WRITE(*   ,890) REACTS(IREACT),MTADD(IREACT),NPOINT,N2,NPT23
C-----OUTPUT SECTION HEAD AND LEADER CARDS AND INTERPOLATION LAW.
C-----USE Q VALUE AS INPUT WITH BACKGROUND = NO C2 CHANGE HERE
      CALL CONTO
      CALL CARDO(C1,C2,L1,L2,IONE,NPT23)
      NBTO(1)=NPT23
      CALL TERPO(NBTO,INTO,1)
C-----IF REQUIRED INSERT FIRST BACKGROUND POINT.
      IPT23=0
      IF(E3LST.GE.ERANGE(1)) GO TO 260
      IPT23=1
      ETAB23(1)=E3LST
      SIG23(1)=0.0
C-----COPY FILE 2 POINTS TO FILE 3 ARRAYS AND OUTPUT A PAGE AT A TIME.
  260 NPT2=1
  270 IPT2=NPT2+NPAGM1
      IF(IPT2.GT.NPOINT) IPT2=NPOINT
      IF(NPOINT.GT.NPAGE) READ(ISCR2) ETAB2X,SIG2X
      JPT2=(IPT2-NPT2)+1
      DO 290 I=1,JPT2
      IPT23=IPT23+1
      IF(IPT23.LE.NPAGE) GO TO 280
C-----IF REQUEST MAKE ALL NEGATIVE CROSS SECTIONS = 0
      IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,NPAGE
      IF(SIG23(KP).LT.0.0D+00) THEN
      MINUS3=MINUS3+1
      SIG23(KP)=0.0
      ENDIF
      ENDDO
      ENDIF
      CALL POINTO(ETAB23,SIG23,NPAGE)
      IPT23=1
  280 ETAB23(IPT23)=ETAB2(I)
      SIG23(IPT23)=SIG2(INDEX,I)
  290 NPT2=NPT2+NPAGE
      IF(NPT2.LE.NPOINT) GO TO 270
C-----IF REQUIRED INSERT LAST BACKGROUND POINT.
      IF(E3.LE.ERANGE(NRANGE)) GO TO 310
      IPT23=IPT23+1
      IF(IPT23.LE.NPAGE) GO TO 300
C-----IF REQUEST MAKE ALL NEGATIVE CROSS SECTIONS = 0
      IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,NPAGE
      IF(SIG23(KP).LT.0.0D+00) THEN
      MINUS3=MINUS3+1
      SIG23(KP)=0.0
      ENDIF
      ENDDO
      ENDIF
      CALL POINTO(ETAB23,SIG23,NPAGE)
      IPT23=1
  300 ETAB23(IPT23)=E3
      SIG23(IPT23)=0.0
C-----IF REQUEST MAKE ALL NEGATIVE CROSS SECTIONS = 0
  310 IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,IPT23
      IF(SIG23(KP).LT.0.0D+00) THEN
      MINUS3=MINUS3+1
      SIG23(KP)=0.0
      ENDIF
      ENDDO
      ENDIF
C-----OUTPUT LAST PAGE OF POINTS.
      CALL POINTO(ETAB23,SIG23,IPT23)
C-----COPY SEND LINE.
      CALL COPYS
      GO TO 740
C-----------------------------------------------------------------------
C
C     RESONANCE AND BACKGROUND CONTRIBUTIONS MUST BE ADDED TOGETHER.
C
C-----------------------------------------------------------------------
  320 E3=ETAB3(IPT3)
      XC3=SIG3(IPT3)
      E3LST=E3
      XC3LST=XC3
C-----INSURE FIRST PAGE OF RESONANCE DATA IS IN CORE AND DEFINE FIRST
C-----POINT.
      IF(NPOINT.LE.NPAGE) GO TO 330
      NPT2=NPAGE
      READ(ISCR2) ETAB2X,SIG2X
  330 LEFT2=LEFT2-NPT2
      E2=ETAB2(1)
      XC2=SIG2(INDEX,1)
      E2LST=E2
      XC2LST=XC2
C-----SELECT LOWEST ENERGY.
  340 IF(E2.eq.E3) go to 470
      IF(E2.gt.E3) go to 350
C-----TREAT SMALL DIFFERENCES AS EQUALITY (SAME TO ABOUT 9 DIGITS).
      IF(DABS(E2-E3).LE.OKDIFF*E2) GO TO 470
      GO TO 360
  350 IF(DABS(E3-E2).LE.OKDIFF*E3) GO TO 460
      GO TO 430
C
C     FILE 2 ENERGY IS LOWER. INTERPOLATE FILE 3 DATA, OR DEFINE TO
C     BE ZERO IF FILE 3 DOES NOT SPAN ENERGY RANGE.
C
  360 E23=E2
      IF(KPT3.GT.1) GO TO 370
      XC23=XC2
      GO TO 390
  370 IF(E3.GT.E3LST) GO TO 380
      XC23=XC2+XC3
      GO TO 390
  380 XC23=XC2+((E3-E2)*XC3LST+(E2-E3LST)*XC3)/(E3-E3LST)
C-----DEFINE NEXT FILE 2 POINT (USE EITHER NEXT POINT - FROM CORE OR
C-----SCRATCH - OR EXTEND CROSS SECTION BEYOND TABULATED RANGE AS
C-----ZERO).
  390 E2LST=E2
      XC2LST=XC2
      IPT2=IPT2+1
      KPT2=KPT2+1
      IF(IPT2.LE.NPT2) GO TO 420
      IF(LEFT2.GT.0) GO TO 410
C-----EXTEND CROSS SECTION AS ZERO (IF CROSS SECTION IS NOT YET ZERO
C-----SET CROSS SECTION TO ZERO AT CURRENT ENERGY. FOR ALL FOLLOWING
C-----POINTS LEAVE CROSS SECTION EQUAL TO ZERO AND SET CURRENT ENERGY
C-----EQUAL TO CURRENT ENERGY FROM FILE 3.)
      IF(XC2.LE.0.0) GO TO 400
      XC2=0.0
      GO TO 550
  400 E2=E3
      GO TO 550
C-----LOAD NEXT PAGE FROM SCRATCH AND RE-INITIALIZE IN CORE INDEX.
  410 IF(NPT2.GT.LEFT2) NPT2=LEFT2
      READ(ISCR2) ETAB2X,SIG2X
      LEFT2=LEFT2-NPT2
      IPT2=1
C-----USE NEXT POINT FROM CORE.
  420 E2=ETAB2(IPT2)
      XC2=SIG2(INDEX,IPT2)
      GO TO 550
C
C     FILE 3 ENERGY IS LOWER. INTERPOLATE FILE 2 DATA, OR DEFINE TO
C     BE ZERO IF FILE 2 DOES NOT SPAN ENERGY RANGE.
C
  430 E23=E3
      IF(KPT2.GT.1) GO TO 440
      XC23=XC3
      GO TO 510
  440 IF(E2.GT.E2LST) GO TO 450
      XC23=XC2+XC3
      GO TO 510
  450 XC23=XC3+((E2-E3)*XC2LST+(E3-E2LST)*XC2)/(E2-E2LST)
      GO TO 510
C
C     ENERGIES ARE VERY CLOSE AND WILL BE TREATED AS EQUAL. SELECT THE
C     SMALLER ENERGY FOR OUTPUT.
C
  460 E2=E3
C
C     ENERGIES ARE EQUAL.
C
  470 E23=E2
      XC23=XC2+XC3
C-----DEFINE NEXT FILE 2 POINT (USE EITHER NEXT POINT - FROM CORE OR
C-----SCRATCH - OR EXTEND CROSS SECTION BEYOND TABULATED RANGE AS
C-----ZERO).
      E2LST=E2
      XC2LST=XC2
      IPT2=IPT2+1
      KPT2=KPT2+1
      IF(IPT2.LE.NPT2) GO TO 500
      IF(LEFT2.GT.0) GO TO 490
C-----EXTEND CROSS SECTION AS ZERO (IF CROSS SECTION IS NOT YET ZERO
C-----SET CROSS SECTION TO ZERO AT CURRENT ENERGY. FOR ALL FOLLOWING
C-----POINTS LEAVE CROSS SECTION EQUAL TO ZERO AND SET CURRENT ENERGY
C-----EQUAL TO CURRENT ENERGY FROM FILE 3.)
      IF(XC2.LE.0.0) GO TO 480
      XC2=0.0
      GO TO 510
  480 E2=E3
      GO TO 510
C-----LOAD NEXT PAGE FROM SCRATCH AND RE-INITIALIZE IN CORE INDEX.
  490 IF(NPT2.GT.LEFT2) NPT2=LEFT2
      READ(ISCR2) ETAB2X,SIG2X
      LEFT2=LEFT2-NPT2
      IPT2=1
C-----USE NEXT POINT FROM CORE.
  500 E2=ETAB2(IPT2)
      XC2=SIG2(INDEX,IPT2)
C-----DEFINE NEXT FILE 3 POINT (USE EITHER NEXT POINT - FROM CORE OR
C-----FILE - OR EXTEND CROSS SECTION BEYOND TABULATED RANGE AS ZERO).
  510 E3LST=E3
      XC3LST=XC3
      IPT3=IPT3+1
      KPT3=KPT3+1
      IF(IPT3.LE.NPT3) GO TO 540
      IF(LEFT3.GT.0) GO TO 530
C-----EXTEND CROSS SECTION AS ZERO (IF CROSS SECTION IS NOT YET ZERO
C-----SET CROSS SECTION TO ZERO AT CURRENT ENERGY. FOR ALL FOLLOWING
C-----POINTS LEAVE CROSS SECTION EQUAL TO ZERO AND SET CURRENT ENERGY
C-----EQUAL TO CURRENT ENERGY FROM FILE 2.)
      IF(XC3.LE.0.0) GO TO 520
      XC3=0.0
      GO TO 550
  520 E3=E2
      GO TO 550
C-----LOAD NEXT PAGE FROM INPUT FILE AND RE-INITIALIZE IN CORE INDEX.
  530 IF(NPT3.GT.LEFT3) NPT3=LEFT3
      CALL POINTI(ETAB3,SIG3,NPT3)
      LEFT3=LEFT3-NPT3
      IPT3=1
C-----USE NEXT POINT FROM CORE.
  540 E3=ETAB3(IPT3)
      XC3=SIG3(IPT3)
C-----STORE COMBINED POINT.
  550 IPT23=IPT23+1
      IF(IPT23.LE.NPAGE) GO TO 570
      IF(NPT23.LE.0.AND.IUSE23.GT.0) REWIND ISCR23
      NPT23=NPT23+NPAGE
      WRITE(ISCR23) ETAB23,SIG23
C-----IF REQUESTED PRINT MESSAGE EVERYTIME A PAGE IS OUTPUT TO SCRATCH.
      IF(MONITR.EQ.0) GO TO 560
      CALL TIMER1(SECONDS)
      CALL OUT9(ETAB23(    1),FIELD(1,1))
      CALL OUT9(ETAB23(NPAGE),FIELD(1,2))
      WRITE(OUTP,1000) NPT23,((FIELD(M,I),M=1,11),I=1,2),SECONDS
      WRITE(*   ,1000) NPT23,((FIELD(M,I),M=1,11),I=1,2),SECONDS
  560 IPT23=1
C
C     SAVE NEXT COMBINED POINT.
C
  570 ETAB23(IPT23)=E23
      SIG23(IPT23)=XC23
C
C     ONLY ALLOW UP TO 2 ENERGY POINTS WITH ZERO CROSS SECTION AT THE
C     BEGINNING OF THE TABLE.
C
      IF(KMPLUS.NE.0) GO TO 590
      IF(XC23.NE.0.0) GO TO 580
      IF(IPT23.LE.2) GO TO 590
C-----KEEP SHIFTING POINTS FORWARD AND RESETTING POINT COUNT.
      ETAB23(2)=ETAB23(IPT23)
      SIG23(2)=SIG23(IPT23)
      IPT23=2
      GO TO 590
  580 KMPLUS=1
C-----CONTINUE COMBINING POINTS IF ANYMORE LEFT.
  590 IF(IPT2.LE.NPT2.OR.LEFT2.GT.0) GO TO 340
C-----WHEN OUT OF FILE 2 POINTS INSURE CURRENT FILE 2 ENERGY (WITH ZERO
C-----CROSS SECTION) IS EQUAL TO FILE 3 ENERGY. NOTE, WHEN ENERGIES ARE
C-----EQUAL NEXT E2 IS DEFINED BEFORE NEXT E3. THEREFORE WHEN E2 IS
C-----BEING EXTEND TO HIGHER ENERGY WITH ZERO CROSS SECTION E2 DEFINED
C-----ABOVE WILL BE SET EQUAL TO THE PRECEDING E3. THE FOLLOWING
C-----STATEMENT WILL INSURE THAT IN THIS CASE E2 IS EXTEND AS EQUAL TO
C-----THE CURRENT E3.
      IF(XC2.EQ.0.0) E2=E3
      IF(IPT3.LE.NPT3.OR.LEFT3.GT.0) GO TO 340
C-----END OF SECTION. DEFINE COMBINED POINT COUNT AND IF REQUIRED SET UP
C-----SCRATCH FILE TO READ.
      NPT23=NPT23+IPT23
      IF(NPT23.LE.NPAGE) GO TO 600
      WRITE(ISCR23) ETAB23,SIG23
      END FILE ISCR23
      REWIND ISCR23
      IUSE23=1
C-----IF REQUESTED PRINT MESSAGE EVERYTIME A PAGE IS OUTPUT TO SCRATCH.
      IF(MONITR.EQ.0) GO TO 600
      CALL TIMER1(SECONDS)
      CALL OUT9(ETAB23(    1),FIELD(1,1))
      CALL OUT9(ETAB23(IPT23),FIELD(1,2))
      WRITE(OUTP,1000) NPT23,((FIELD(M,I),M=1,11),I=1,2),SECONDS
      WRITE(*   ,1000) NPT23,((FIELD(M,I),M=1,11),I=1,2),SECONDS
C-----OUTPUT SECTION HEAD AND LEADER CARDS AND INTERPOLATION LAW.
  600 WRITE(OUTP,870) REACTS(IREACT),MTADD(IREACT),NPOINT,N2,NPT23
      WRITE(*   ,870) REACTS(IREACT),MTADD(IREACT),NPOINT,N2,NPT23
      CALL CONTO
C-----USE Q VALUE FROM BACKGROUND = NO C2 CHANGE HERE
      CALL CARDO(C1,C2,L1,L2,IONE,NPT23)
      NBTO(1)=NPT23
      CALL TERPO(NBTO,INTO,1)
C-----OUTPUT DATA A PAGE AT A TIME.
      IPTB=1
  610 IF(NPT23.GT.NPAGE) READ(ISCR23) ETAB23,SIG23
      IPT23=IPTB+NPAGM1
      IF(IPT23.GT.NPT23) IPT23=NPT23
      JPT23=(IPT23-IPTB)+1
C-----IF REQUEST MAKE ALL NEGATIVE CROSS SECTIONS = 0
      IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,JPT23
      IF(SIG23(KP).LT.0.0D+00) THEN
      MINUS3=MINUS3+1
      SIG23(KP)=0.0
      ENDIF
      ENDDO
      ENDIF
      CALL POINTO(ETAB23,SIG23,JPT23)
      IPTB=IPTB+NPAGE
      IF(IPTB.LE.NPT23) GO TO 610
C-----COPY SEND LINE.
      CALL COPYS
      GO TO 740
C-----------------------------------------------------------------------
C
C     BEYOND REQUIRED SECTION (I.E. NO BACKGROUND). SAVE CURRENT LINE
C     AND SET FLAG.
C
C-----------------------------------------------------------------------
  620 SAVE1(1)=C1H
      SAVE1(2)=C2H
      DO 630 I=1,7
  630 ISAVE1(I)=IHEAD(I)
      MF3END=-1
C-----------------------------------------------------------------------
C
C     NO BACKGROUND. IF FISSION CROSS SECTION AND ALL FISSION WIDTHS
C     ARE ZERO SIMPLY SKIP THE REACTION (I.E., NO OUTPUT). OTHERWISE
C     USE INPUT PARAMETERS TO DECIDE TO EITHER HAVE NO OUTPUT OR TO
C     OUTPUT RESONANCE CONTRIBUTION.
C
C-----------------------------------------------------------------------
C-----DO NOT OUTPUT FISSION CROSS SECTION UNLESS FISSION WIDTHS ARE NOT
C-----ZERO.
  640 IF(MTADD(IREACT).EQ.18.OR.MTADD(IREACT).EQ.19) THEN
      IF(LFWX.LE.0) GO TO 750
      ENDIF
C-----DO NOT OUTPUT COMPETITIVE LEVEL IF NO BACKGROUND
      IF(MTADD(IREACT).GE.600) GO TO 750 ! CHARGED PARTICLES
      IF(MTADD(IREACT).EQ. 50) GO TO 750 ! (N,N')
C-----IF IN EDIT MODE PRINT MESSAGE AND CONTINUE.
      IF(IMEDIT.NE.2) GO TO 650
      WRITE(OUTP,860) REACTS(IREACT),MTADD(IREACT)
      WRITE(*   ,860) REACTS(IREACT),MTADD(IREACT)
      GO TO 750
C-----USE INPUT OPTION TO DECIDE WHETHER OR NOT TO OUTPUT SECTION.
C-----07/10/16 - CORRECTED FOR NO MT=19 OUTPUT IF NO BACKGHROUND
  650 IF(MTADD(IREACT).NE.19) GO TO 660
C-----NO OUTPUT. PRINT MESSAGE EXPLAINING WHY.
      WRITE(OUTP,950) REACTS(IREACT),MTADD(IREACT),IZERO,IZERO,IZERO
      WRITE(*   ,950) REACTS(IREACT),MTADD(IREACT),IZERO,IZERO,IZERO
      GO TO 750
  660 IF(IMBACK.EQ.1) GO TO 670
C-----NO OUTPUT. PRINT MESSAGE EXPLAINING WHY.
      WRITE(OUTP,940) REACTS(IREACT),MTADD(IREACT),IZERO,IZERO,IZERO
      WRITE(*   ,940) REACTS(IREACT),MTADD(IREACT),IZERO,IZERO,IZERO
      GO TO 750
C-----OUTPUT RESONANCE CONTRIBUTION. DEFINE OUTPUT SECTION MAT/MF/MT.
  670 MATH=MATNOW
      MFH=3
      MTH=MTADD(IREACT)
C-----USE Q-VALUE DEFINED BY RESONANCE PARAMETER INPUT
      C2OUT = QADD(IREACT)
C
C     DEFINE POSSIBLE OUTSET FROM START OF SECTION TO ALLOW FOR
C     THRESHOLD
C
      IF(NEGADD(IREACT).LE.0) NEGADD(IREACT) = 1
      MSKIP = NEGADD(IREACT) - 1
      MSKIP = 3*(MSKIP/3)
      NOUT2 = NPOINT - MSKIP
C
C     PRINT DESCRIPTION OF SECTION.
C
      IF(IMESS.eq.2) go to 680
      IF(IMESS.gt.2) go to 690
C-----SECTION OF O.K. TO BE COMBINED WITH BACKGROUND.
      WRITE(OUTP,900) REACTS(IREACT),MTADD(IREACT),NOUT2,IZERO,NOUT2
      WRITE(*   ,900) REACTS(IREACT),MTADD(IREACT),NOUT2,IZERO,NOUT2
      GO TO 700
C-----BACKGROUND TEMPERATURE IS NOT 0 KELVIN.
  680 CALL OUT9(TEMP,FIELD(1,1))
      IF(IMEDIT.NE.2) THEN
      WRITE(OUTP,910) REACTS(IREACT),MTADD(IREACT),NOUT2,IZERO,
     1 NOUT2,(FIELD(M,1),M=1,11)
      WRITE(*   ,910) REACTS(IREACT),MTADD(IREACT),NOUT2,IZERO,
     1 NOUT2,(FIELD(M,1),M=1,11)
      WRITE(OUTP,1010) REACTS(IREACT),MTADD(IREACT),
     1 N2,(FIELD(M,1),M=1,11)
      WRITE(*   ,1010) REACTS(IREACT),MTADD(IREACT),
     1 N2,(FIELD(M,1),M=1,11)
      ENDIF
      GO TO 700
C-----INTERPOLATION LAW IS NOT LINEAR.
  690 IF(IMEDIT.NE.2) THEN
      WRITE(OUTP,920) REACTS(IREACT),MTADD(IREACT),NOUT2,IZERO,
     1 NOUT2,(NBTF(I),INTF(I),I=1,N1)
      WRITE(*   ,920) REACTS(IREACT),MTADD(IREACT),NOUT2,IZERO,
     1 NOUT2,(NBTF(I),INTF(I),I=1,N1)
      WRITE(OUTP,1020) REACTS(IREACT),MTADD(IREACT),N2,
     1 (NBTF(I),INTF(I),I=1,N1)
      WRITE(*   ,1020) REACTS(IREACT),MTADD(IREACT),N2,
     1 (NBTF(I),INTF(I),I=1,N1)
      WRITE(OUTP,930)
      WRITE(*   ,930)
      ELSE
      WRITE(OUTP,1030)
      WRITE(*   ,1030)
      ENDIF
C-----OUTPUT SECTION HEAD CARDS AND INTERPOLATION LAW.
  700 CALL CARDO(ZA,AWR,IZERO,IZERO,IZERO,IZERO)
      CALL CARDO(ZEROD,C2OUT,IZERO,IZERO,IONE,NOUT2)
      NBTO(1)=NOUT2
      CALL TERPO(NBTO,INTO,1)
C-----RESONANCE CONTRIBUTION A PAGE AT A TIME.
      IF(NPOINT.GT.NPAGE) REWIND ISCR2
      IPOINT=1
  710 IF(NPOINT.GT.NPAGE) READ(ISCR2) ETAB2X,SIG2X
      KPOINT=IPOINT+NPAGM1
      IF(KPOINT.GT.NPOINT) KPOINT=NPOINT
      KOUT=(KPOINT-IPOINT)+1
      DO 720 K=1,KOUT
  720 SIG3(K)=SIG2(INDEX,K)
C-----IF REQUEST MAKE ALL NEGATIVE CROSS SECTIONS = 0
      IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,KOUT
      IF(SIG3(KP).LT.0.0D+00) THEN
      MINUS3=MINUS3+1
      SIG3(KP)=0.0
      ENDIF
      ENDDO
      ENDIF
      IF(MSKIP.LE.0) THEN
C-----NO THRESHOLD OR NOW ABOVE IT
      CALL POINTO(ETAB2,SIG3,KOUT)
      ELSE
C-----SKIP TO THRESHOLD
      IF(MSKIP.GT.KOUT) GO TO 730
      KKOUT = KOUT  - MSKIP
      NOUT1 = MSKIP + 1
      CALL POINTO(ETAB2(NOUT1),SIG3(NOUT1),KKOUT)
  730 MSKIP = MSKIP - KOUT
      IF(MSKIP.LE.0) MSKIP = 0
      ENDIF
      IPOINT=IPOINT+NPAGE
      IF(IPOINT.LE.NPOINT) GO TO 710
C-----ADD SEND LINE.
      CALL OUTS(MATH,MFH)
C
C     END OF SECTION LOOP
C
C-----PRINT WARNING IF CROSS SECTION IS NEGATIVE.
  740 IF(MINUS3.GT.0) WRITE(OUTP,980) MINUS3
      IF(MINUS3.GT.0) WRITE(*   ,980) MINUS3
C-----PRINT WARNING MESSAGE IF CROSS SECTION IS NOT POSITIVE AT ANY ENER
C-----END OF REACTION LOOP.
      IF(IMPLUS.LE.0) WRITE(OUTP,990)
      IF(IMPLUS.LE.0) WRITE(*   ,990)
  750 CONTINUE
      WRITE(OUTP,960)
      WRITE(*   ,960)
C-----------------------------------------------------------------------
C
C     ALL RESONANCE REACTIONS HAVE BEEN OUTPUT. TEST FOR CURRENT STATUS
C     MF3END=-1 - STILL IN FILE 3. RESTORE LAST LINE READ, OUTPUT AND
C                 COPY TO END OF FILE 3.
C           = 0 - STILL IN FILE 3. COPY TO END OF FILE 3.
C           = 1 - END OF FILE 3 READ. OUTPUT FEND LINE.
C           = 2 - PAST END OF FILE READ. OUTPUT FEND LINE. RESTORE LAST
C                 LINE AND OUTPUT.
C     IF ALL CASES COPY REMAINDER OF MAT.
C
C-----------------------------------------------------------------------
      IF(MF3END.lt.0) go to 760
      IF(MF3END.eq.0) go to 780
      IF(MF3END.lt.1) go to 780
      IF(MF3END.eq.1) go to 790
      go to 800
C-----STILL IN FILE 3. RESTORE LAST LINE READ AND OUTPUT.
  760 C1H=SAVE1(1)
      C2H=SAVE1(2)
      DO 770 I=1,7
  770 IHEAD(I)=ISAVE1(I)
      CALL CONTO
C-----STILL IN FILE 3. COPY TO END OF FILE 3.
  780 CALL COPYF
      GO TO 820
C-----END OF FILE 3 REACHED. OUTPUT FEND.
  790 CALL OUTF(MATH)
      GO TO 820
C-----PAST END OF FILE 3. OUTPUT FEND. RESTORE AND OUTPUT LAST LINE.
  800 CALL OUTF(MATH)
      C1H=SAVE1(1)
      C2H=SAVE1(2)
      DO 810 I=1,7
  810 IHEAD(I)=ISAVE1(I)
      CALL CONTO
C-----IF LAST LINE WAS MEND THERE IS NOTHING ELSE TO DO.
      IF(MATH.LE.0) GO TO 830
C-----COPY REMAINDER OF MAT.
  820 CALL COPYM
C-----PRINT RUNNING TIME
  830 CALL TIMEMAT
      RETURN
  840 FORMAT(1x,78('=')/' Combining File 2 and File 3 Data'/1X,78('=')/
     1 ' Reaction   MT  File 2  File 3  Combined'/
     2 '                Points  Points    Points Comments'/1X,78('='))
  850 FORMAT(1X,78('=')/' Edit of File 3 Data'/1X,78('=')/
     1 ' Reaction   MT  File 3'/
     2 '                Points Comments'/1X,78('='))
  860 FORMAT(1X,A8,I5,8X,' No Background.')
  870 FORMAT(1X,A8,I5,2I8,I10)
  880 FORMAT(1X,A8,I5,I8)
  890 FORMAT(1X,A8,I5,2I8,I10,' Background is Zero at ALL Energies.')
  900 FORMAT(1X,A8,I5,2I8,I10,' No Background. Output Resonance Part.')
  910 FORMAT(1X,A8,I5,2I8,I10,' Background ',11A1,' Kelvin.'/
     2 1X,7('WARNING...'),'WARNING'/
     3 35X,' Cannot Combine 0 Kelvin Resonance and Hot Background.'/
     4 35X,' Will Only Output Resonance Contribution.')
  920 FORMAT(1X,A8,I5,2I8,I10,' Background NOT Linear-Linear'/
     1 35X,'   NBT   INT'/(35X,2I6))
  930 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 35X,' Cannot Combine Resonance and Background.'/
     2 35X,' Background MUST be Linearized.'/
     3 35X,' Will ONLY Output Resonance Contribution.')
  940 FORMAT(1X,A8,I5,2I8,I10,' No Background. No Output (as per Input',
     1 ' Parameter)')
  950 FORMAT(1X,A8,I5,2I8,I10,' No Background. Never Output for first',
     1 ' chance fission')
  960 FORMAT(1X,78('='))
  970 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 1X,A8,I5,2I8,I10,' All Fission Widths are Zero.'/
     3 35X,' WilL Output Background Cross Section.')
  980 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 35X,' Above Cross Section is Negative at',I6,' Energies.')
  990 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 35X,' Above Cross Section is NOT Positive at ANY Energy')
 1000 FORMAT(25X,I10,11A1,' to',11A1,' eV',F11.2,' Sec.')
 1010 FORMAT(1X,A8,I5,I8,' Background Temperature=',11A1,
     1  ' Kelvin.'/
     2 1X,7('WARNING...'),'WARNING'/
     3 35X,' Cannot Combine 0 Kelvin Resonance and Hot Background.'/
     4 35X,' Will Only Output Resonance Contribution.')
 1020 FORMAT(1X,A8,I5,I8,' Background NOT Linear-Linear'/
     1 15X,'   NBT   INT'/(15X,2I6))
 1030 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 35X,' Cannot Combine Resonance and Background.'/
     2 35X,' Background MUST be Linearized.'/
     3 35X,' Will ONLY Output Resonance Contribution.')
      END
      SUBROUTINE RDBW
C=======================================================================
C
C     READ BREIT-WIGNER SINGLE OR MULTI LEVEL DATA FOR ONE ENERGY RANGE.
C     EACH L STATE WILL BE TREATED AS A SEPERATE SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*40 LRXHOL
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/NAPRHO/NRO,NAPS
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RANGER/LOW,LHI
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
      DIMENSION LRXHOL(3)
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DATA LRXHOL/
     1 '(No Competitive Widths)                 ',
     2 '(WARNING - Defined LC=LRX-1)            ',
     3 '(ERROR - < 0 Illegal - Will Assume = 0) '/
c       1234567890123456789012345678901234567890
      DATA C3/1.23D-01/
      DATA C5/8.0D-02/
C-----UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITTEE RECOMMENDATION
c***** Updated 9/01/04 - ENDF-102, appendix H
      DATA C4/1.00866491578D+00/
c***** Updated 9/01/04
      DATA C6/2.196807122623D-03/
C-----INITIALIZE NEGATIVE WIDTH COUNT
      NEGSUM = 0
C-----DEFINE TARGET SPIN, SPIN UP SCATTERING RADIUS AND NUMBER OF L
C-----STATES.
      CALL CARDIO(SPI,AP,L1,L2,NLS,N2)
      CALL OUT9(SPI,FIELD(1,1))
      IF(NRO.EQ.0) GO TO 10
      WRITE(OUTP,150) (FIELD(M,1),M=1,11),NLS
      WRITE(*   ,150) (FIELD(M,1),M=1,11),NLS
      GO TO 20
   10 CALL OUT9(AP ,FIELD(1,2))
      WRITE(OUTP,160) ((FIELD(M,I),M=1,11),I=1,2),NLS
      WRITE(*   ,160) ((FIELD(M,I),M=1,11),I=1,2),NLS
C
C     L STATE LOOP.
C
   20 DO 130 IL=1,NLS
C-----DEFINE ATOMIC WEIGHT RATIO, Q VALUE FOR COMPETITIVE REACTION,
C-----L VALUE, COMPETITIVE REACTION FLAG AND NUMBER OF RESONANCES.
      CALL CARDIO(AWRI,QX,LNOW,LRX,NRS6X,NRS)
C-----CHECK LRX
      LRXOUT = 1               ! LRX = 0
      IF(LRX.GT.0) LRXOUT = 2  ! LRX > 0
      IF(LRX.LT.0) LRXOUT = 3  ! LRX < 0
      CALL OUT9(AWRI,FIELD(1,1))
      CALL OUT9(QX  ,FIELD(1,2))
      WRITE(OUTP,170) ((FIELD(M,I),M=1,11),I=1,2),
     1 LNOW,LRX,LRXHOL(LRXOUT)
      WRITE(*   ,170) ((FIELD(M,I),M=1,11),I=1,2),
     1 LNOW,LRX,LRXHOL(LRXOUT)
      IF(LRX.LT.0) LRX = 0   ! SET AFTER OUTPUT LISTING
      IF(LRX.GT.0) THEN
      WRITE(OUTP,180) LRX-1
      WRITE(*   ,180) LRX-1
      ENDIF
      WRITE(OUTP,190) NRS
      WRITE(*   ,190) NRS
C
C     CONSTRUCT TABLE OF PERMISSIBLE J-VALUES
C     (THIS TABLE WILL BE USED TO CHECK FOR MISSING OR ILLEGAL J-VALUES)
C
      CALL SETJ(LNOW)
C-----INITIALIZE FLAGS TO CHECK COMPETITIVE WIDTHS.
      LRXIN=0
      LRXNEG=0
C-----INCREMENT SECTION COUNT AND INDICES TO RESONANCE PARAMETER TABLE.
      CALL LIMIT1(1)
C-----DEFINE ALL PARAMETERS FOR A SECTION.
      AWRICM=AWRI/(AWRI+ONE)
C-----EITHER CALCULATE CHANNEL RADIUS OR DEFINE TO BE EQUAL TO THE
C-----SCATTERING RADIUS.
      IF(NAPS.NE.0) GO TO 30
      APX=C3*((C4*AWRI)**(ONE/THREE))+C5
      GO TO 40
   30 APX=AP
   40 AK1=C6*AWRICM
      AK2=AK1*AK1
      BETA(NSECT)=PI*ABN/AK2
      RHOX2(NSECT)=APX*APX*AK2
      RHOP1(NSECT)=AP*AK1
      EL(NSECT)=ELX
      EH(NSECT)=EHX
      QVALUE(NSECT)=QX/AWRICM
      NLOW(NSECT)=LOW
      NHIGH(NSECT)=LHI
      LVALUE(NSECT)=LNOW
c-----2016/11/16 - Added for L-state dependent fission.
      LFWL(NSECT) = 0
c-----10/05/28 - added for potential correction
      POTL(NSECT) = 4*LNOW + 2  ! 2(2l+1) total weight of all sequences
      ADDL(NSECT) = 0.0D+0      ! initialize weight of missing sequences
      LRXTAB(NSECT)=LRX
      NAPTAB(NSECT)=NAPS
      MODE(NSECT)=LRF
C
C     06/09/05 - THIS WAS ALREADY DONE FOR BREIT-WIGNER - NOW ALL.
C
C     IF ENERGY DEPENDENT SCATTERING RADIUS DEFINE INDEX TO TABULATED
C     DATA (OTHERWISE INDEX HAS ALREADY BEEN INITIALIZED TO 0).
C     CALCULATE ENERGY DEPENDENT RHOP FROM SCATTERING RADII.
C
      IF(NRO.GT.0) THEN
      NRHO(NSECT)=NUMRHO
      INX3=INXRHO(3,NUMRHO)
      INX4=INXRHO(4,NUMRHO)
      DO I=INX3,INX4
      RHOTAB(I)=APTAB(I)*AK1
      ENDDO
      ENDIF
C
C     READ RESONANCE PARAMETERS (6 PER RESONANCE).
C
C     (1) ENERGY
C     (2) J
C     (3) TOTAL WIDTH
C     (4) NEUTRON WIDTH
C     (5) CAPTURE WIDTH
C     (6) FISSION WIDTH
C
C     SOME OF THESE PARAMETERS WILL BE CONVERTED TO THE FORM REQUIRED
C     FOR CALCULATIONS.
C
C     (2) STATISTICAL WEIGHT, GJ (J NOT NEEDED)
C     (3) COMPETITIVE WIDTH DIVIDED BY PENETRATION FACTOR (TOTAL WIDTH
C         NOT NEEDED)
C     (4) NEUTRON WIDTH DIVIDED BY PENETRATION FACTOR (NEUTRON WIDTH IS
C         ONLY USED IN A FORM DIVIVED BY THE PENETRATION FACTOR).
C
      CALL LISPAR6(LOW,NRS6X)
      IF(IMEDIT.EQ.0) GO TO 50
      IF(LRF.EQ.1) WRITE(OUTP,200)
      IF(LRF.EQ.2) WRITE(OUTP,210)
      IF(LRF.EQ.1) WRITE(*   ,200)
      IF(LRF.EQ.2) WRITE(*   ,210)
C-----FOR EACH SECTION (ISOTOPE, ENERGY RANGE, L VALUE) SORT
C-----RESONANCES IN ASCENDING (J, E) ORDER.
   50 CALL ORDER(LOW,LHI)
      AJNOW=RESTAB(2,LOW)
C-----CHECK FOR LEGAL J-VALUE AND DEFINE STATISTICAL WEIGHT.
      CALL CHECKJ(AJNOW,LNOW,GJ,LRF)
      DO 120 JR=LOW,LHI
c***** DEBUG
c-----Change LRF=3 to LRF=2 if input (LRFIN) was actually = 3
      if(LRFIN.eq.3) then
      W3 = DABS(RESTAB(3,JR))     ! elastic
      W4 = DABS(RESTAB(4,JR))     ! capture
      W5 = DABS(RESTAB(5,JR))     ! fission
      W6 = DABS(RESTAB(6,JR))
      RESTAB(6,JR) = W5 + W6      ! fission
      RESTAB(5,JR) = W4           ! capture
      RESTAB(4,JR) = W3           ! elastic
      RESTAB(3,JR) = W3+W4+W5+W6  ! total
      endif
c***** DEBUG
C-----IF FISSION WIDTH IS NOT ZERO TURN ON FISSILE FLAG.
      IF(DABS(RESTAB(6,JR)).NE.0.0) LFWX=1
      IF(DABS(RESTAB(6,JR)).NE.0.0) LFWL(NSECT) = LFWL(NSECT) + 1
C-----INSERT DIVIDING LINE BETWEEN DIFFERENT J VALUES.
      IF(DABS(AJNOW-RESTAB(2,JR)).LE.0.01) GO TO 60
      AJNOW=RESTAB(2,JR)
C-----CHECK FOR LEGAL J-VALUE AND DEFINE STATISTICAL WEIGHT.
      CALL CHECKJ(AJNOW,LNOW,GJ,LRF)
      IF(IMEDIT.NE.0) WRITE(OUTP,230)
C-----SEE IF TOTAL EQUALS ELASTC + CAPTURE + FISSION.
   60 GAMC=RESTAB(3,JR)-(RESTAB(4,JR)+RESTAB(5,JR)+RESTAB(6,JR))
C-----IF DIFFERENCE IS VERY SMALL COMPARED TO TOTAL SET IT TO ZERO.
      IF(DABS(GAMC).LE.0.00001*RESTAB(3,JR)) GAMC=0.0
C-----CHECK FOR NEGATIVE DIFFERENCES AND SIGNIFICANT POSITIVE DIFFERENCE
      IF(GAMC.LT.0.0) LRXNEG=LRXNEG+1
      IF(GAMC.NE.0.0) LRXIN=LRXIN+1
C-----LIST RESONANCE PARAMETERS.
C-----CHECK FOR NEGATIVE WIDTHS
      NEGRES = 0
      DO I=3,6
      IF(RESTAB(I,JR).LT.0.0D+0) NEGRES = NEGRES + 1
      ENDDO
      IF(IMEDIT.EQ.0.AND.NEGRES.EQ.0) GO TO 90
      DO 70 I=1,6
   70 CALL OUT9(RESTAB(I,JR),FIELD(1,I))
      CALL OUT9(GAMC        ,FIELD(1,7))
      WRITE(OUTP,220) (FIELD(M,1),M=1,11),RESTAB(2,JR),
     1 ((FIELD(M,I),M=1,11),I=3,7)
      IF(NEGRES.GT.0) THEN
      NEGSUM = NEGSUM + NEGRES
      WRITE(   *,220) (FIELD(M,1),M=1,11),RESTAB(2,JR),
     1 ((FIELD(M,I),M=1,11),I=3,7)
      WRITE(   *,80) NEGRES
      WRITE(OUTP,80) NEGRES
   80 FORMAT(' ERROR - ',I1,' Negative Widths')
      ENDIF
C-----REPLACE J VALUE BY STATISTICAL WEIGHT.
   90 RESJTAB(JR) = RESTAB(2,JR)
      RESTAB(2,JR)= GJ
C-----04/21/07 - IF NECESSARY DEFINE RHOX2
      IF(NRO.EQ.1.AND.NAPS.EQ.1) THEN
      ERABS=DABS(RESTAB(1,JR))
      CALL SETRHO1(ERABS,NSECT)
      ENDIF
C-----DEFINE COMPETITIVE WIDTH DIVIDED BY PENETRATION FACTOR.
C-----(WILL REPLACE TOTAL WIDTH IN RESTAB).
      IF(LRX.LE.0) GO TO 100
C-----IF NEGATIVE OR ZERO COMPETITIVE WIDTH SET IT TO ZERO.
      IF(GAMC.LE.0.0) GO TO 100
      RHOZ2=DABS(DABS(RESTAB(1,JR))+QVALUE(NSECT))*RHOX2(NSECT)
C-----08/30/04 - CHANGED TO ||Er| + Q| - from |Er + Q|
c-----           These are the SAME if Er > 0
C-----04/09/10 - CHANGED FROM FACTS2 TO FACTS3 - SHIFT NOT USED.
      CALL FACTS3(LRX-1,RHOZ2,PENFACC)
C-----02/14/04 - PROTECT AGAINST ZERO ENERGY RESONANCES
      IF(PENFACC.NE.0.0D+00) THEN
      RESTAB(3,JR)=GAMC/PENFACC
      ELSE
      RESTAB(3,JR)=0.0D+00
      ENDIF
      GO TO 110
  100 RESTAB(3,JR)=0.0
C-----DEFINE NEUTRON WIDTH DIVIDED BY PENETRATION FACTOR.
  110 RHO2X=DABS(RESTAB(1,JR))*RHOX2(NSECT)
      CALL FACTS2(LVALUE(NSECT),RHO2X,SHIFT2(JR),PENFAC)
C-----02/14/04 - PROTECT AGAINST ZERO ENERGY RESONANCES
      IF(PENFAC.NE.0.0D+00) THEN
      RESTAB(4,JR)=RESTAB(4,JR)/PENFAC
      ELSE
      RESTAB(4,JR)=0.0D+00
      ENDIF
  120 CONTINUE
C-----PRINT ERROR MESSAGES IF WIDTHS DO NOT ADD UP OR COMPETITIVE
C-----WIDTHS ARE NEGATIVE, OR LRX=0 AND ALL COMPETITIVE WIDTHS ARE
C-----ZERO.
      IF(LRX.EQ.0.AND.LRXIN.GT.0) WRITE(OUTP,240) LRXIN
      IF(LRX.GT.0.AND.LRXNEG.GT.0) WRITE(OUTP,250) LRXNEG
      IF(LRX.GT.0.AND.LRXIN.EQ.0) WRITE(OUTP,260)
C-----IF NO POSITIVE COMPETITIVE WIDTHS SET FLAG OFF.
      IF(LRXIN.EQ.0) LRXTAB(NSECT)=0
C
C     FOR MULTI-LEVEL PARAMETERS PRINT WARNING IF ANY J SEQUENCES ARE
C     MISSING AND ADD MISSING J SEQUENCE WITH NO RESONANCES IN ORDER TO
C     ALLOW POTENTIAL SCATTERING TO BE CORRECTLY CALCULATED.
C
      IF(LRF.EQ.2) THEN
      LSECT = NSECT
      CALL MISSINGJ(LNOW,LSECT)
      ENDIF
C-----END OF L STATE LOOP.
  130 CONTINUE
C
C     STOP IF NEGATIVE WIDTHS
C
      IF(NEGSUM.GT.0) THEN
      WRITE(   *,140) NEGSUM
      WRITE(OUTP,140) NEGSUM
  140 FORMAT(///
     1 ' ERROR - Execution Terminated because of',i8,
     1 ' Negative Widths.'///)
      CALL ENDERROR
      ENDIF
      RETURN
  150 FORMAT(1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Energy Dependent Scattering Radius---(LISTED ABOVE)'/
     3       ' Number of L Values-------------------',I11)
  160 FORMAT(1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Effective Scattering Radius (A+)-----',11A1/
     3       ' Number of L Values-------------------',I11)
  170 FORMAT(1X,78('=')/
     1       ' Atomic Weight Ratio of Isotope-------',11A1/
     2       ' Competitive Reaction Q Value---------',11A1,' eV'/
     3       ' Angular Momentum (L)-----------------',I11/
     4       ' LRX (Competitive Widths)-------------',I11,1X,A40)
  180 FORMAT(
     1       ' LC (Competitive L Value: LRX-1)------',I11)
  190 FORMAT(
     5       ' Number of Resonances-----------------',I11)
  200 FORMAT(1X,78('=')/' Single Level Breit-Wigner',
     1 ' Resonance Parameters'/1X,78('=')/
     1 '      Energy  J Value      Total    Neutron',
     2 '    Capture    Fission Competition'/
     2 8X,'(eV)',9X,4(7X,'(eV)'),8X,'(eV)'/1X,78('='))
  210 FORMAT(1X,78('=')/' Multi-Level Breit-Wigner',
     1 ' Resonance Parameters'/1X,78('=')/
     1 '      Energy  J Value      Total    Neutron',
     2 '    Capture    Fission Competition'/
     2 8X,'(eV)',9X,4(7X,'(eV)'),8X,'(eV)'/1X,78('='))
  220 FORMAT(1X,11A1,F7.2,2X,44A1,1X,11A1)
  230 FORMAT(1X,78('-'))
  240 FORMAT(1X,78('=')/1X,7('WARNING...'),'WARNING'/
     1 ' LRX=0 IndicateS ALL Competitive Widths Should be Zero and'/
     2 ' Data Should have Total Width = (Elastic+Capture+Fission).'/
     3 I5,' Resonances do NOT Satisfy this Condition.'/
     4 ' This Program will Assume that LRX is Correct and that ALL'/
     5 ' Competitive Widths are Zero and will Define the Total Width'/
     6 ' to be the Sum of Elastic+Capture+Fission for ALL Resonances.'/
     7 ' Check Evaluated Data. Correct Either Sum of Widths or LRX.')
  250 FORMAT(1X,78('=')/1X,7('WARNING...'),'WARNING'/
     1 ' LRX=1 Indicates Competitive Widths Present.'/
     2 I5,' Competitive Widths are Negative.'/
     3 ' Will Set ALL Negative Competitive Widths Equal to Zero.'/
     4 ' Check Evaluated Data. Correct Either Sum of Widths or LRX.')
  260 FORMAT(1X,78('=')/1X,7('WARNING...'),'WARNING'/
     1 ' LRX=1 Indicates Competitive Widths Present.'/
     2 ' There are NO Positive Competitive Widths.'/
     3 ' Will Assume Competitive Widths are ALL Zero.'/
     4 ' Check Evaluated Data. Correct Either Sum of Widths or LRX.')
      END
      SUBROUTINE GJWAIT(SPI,AJ,GJ)
C=======================================================================
C
C     DEFINE STATISTICAL WEIGHT BASED ON S AND J.
C
C     NOTE - THIS SHOULD BE THE ONLY PLACE IN THE PROGRAM THAT DEFINES
C            STATISTICAL WEIGHTS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      GJ=(TWO*AJ+ONE)/(FOUR*SPI+TWO)
      RETURN
      END
      SUBROUTINE SETJ(LNOW)
C=======================================================================
C
C     USE CHANNEL SPIN AND L VALUE TO DEFINE TABLE OF LEGAL J VALUES
C     AND STATISTICAL WEIGHTS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/JFIX/GJSET(50),GJBAD,TABJ(50),JOK(50),JMISS,NUMJ,IGJBAD
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
C-----INITIALIZE SUM OF STATISTICAL WEIGHTS.
      CALL SUMJ(GJBAD,LNOW,1)
C-----CHANNEL SPIN I+0.5 FIRST
      NUMJ=1
      ELL=LNOW
      CHSPIN=SPI+0.5D+00
      TABJ(NUMJ)=DABS(CHSPIN-ELL)
      ENDJ=CHSPIN+ELL
   10 IF(TABJ(NUMJ).EQ.ENDJ) GO TO 20
      TABJ(NUMJ+1)=TABJ(NUMJ)+1.0D+00
      NUMJ=NUMJ+1
      GO TO 10
C-----CHANNEL SPIN I-0.5 NEXT.
   20 IF(SPI.EQ.0.0) GO TO 40
      NUMJ=NUMJ+1
      CHSPIN=DABS(SPI-0.5D+00)
      TABJ(NUMJ)=DABS(CHSPIN-ELL)
      ENDJ=CHSPIN+ELL
   30 IF(TABJ(NUMJ).EQ.ENDJ) GO TO 40
      TABJ(NUMJ+1)=TABJ(NUMJ)+1.0
      NUMJ=NUMJ+1
      GO TO 30
C-----SORT J VALUES INTO ASCENDING ORDER, ELIMINATE DUPLICATES AND
C-----DEFINE STATISTIC WEIGHTS.
   40 IF(NUMJ.GT.1) CALL SORTS(TABJ,NUMJ)
      NUMIN=0
C-----INITIALIZE ERROR MESSAGE COUNT
      MYERROR = 0
      DO 80 I=1,NUMJ
      IF(I.EQ.1) GO TO 60
      DO 50 II=1,NUMIN
      IF(DABS(TABJ(II)-TABJ(I)).LT.0.01) GO TO 70
   50 CONTINUE
   60 NUMIN=NUMIN+1
      TABJ(NUMIN)=TABJ(I)
      JOK(NUMIN)=0
C-----DEFINE STATISTICAL WEIGHT.
      AJ=TABJ(NUMIN)
      CALL GJWAIT(SPI,AJ,GJSET(NUMIN))
      GO TO 80
C-----DUPLICATE J VALUE. ONLY ADD POTENTIAL IF OPTION IS TURNED ON
C-----AND THESE ARE NOT SINGLE LEVEL PARAMETERS.
   70 IF(LRF.EQ.1) GO TO 80
C-----PRINT WARNING BEFORE FIRST ERROR MESSAGE
      IF(MYERROR.EQ.0) WRITE(OUTP,90)
      MYERROR = MYERROR + 1
      WRITE(OUTP,100) LNOW,TABJ(I)
   80 CONTINUE
C-----PRINT FINAL LINE IF ANY ERROR MESSAGES
      IF(MYERROR.GT.0) WRITE(OUTP,110)
      NUMJ=NUMIN
      RETURN
   90 FORMAT(1X,78('-')/1X,7('WARNING...'),'WARNING')
  100 FORMAT(
     1 ' L=',I2,' J =',F7.3,' Corresponds to 2 Resonance Sequences.')
  110 FORMAT(1X,78('-'))
      END
      SUBROUTINE SUMJ(GJNOW,LNOW,MYWAY)
C=======================================================================
C
C     EITHER,
C     MYWAY = 1 - INITIALIZE SUM OF STATISICAL WEIGHTS (GJ)
C           = 2 - ADD GJ TO SUM
C           = 3 - COMPARE SUM TO (2*L+1) AND PRINT WARNING IF THEY
C                 DO NOT AGREE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/FIXPOT/MYPOT
      COMMON/JFIX/GJSET(50),GJBAD,TABJ(50),JOK(50),JMISS,NUMJ,IGJBAD
      GO TO (10,20,30),MYWAY
C-----INITIALIZE SUM.
   10 GJSUM=0.0
      RETURN
C-----ADD TO SUM
   20 GJSUM=GJSUM+GJNOW
      RETURN
C-----COMPARE SUM TO (2*L+1). PRINT ERROR MESSAGE IF MORE THAN 1
C-----PER-CENT DISAGREEMENT.
   30 GJWANT=2*LNOW+1
      IF(DABS(GJWANT-GJSUM).LE.0.01*GJWANT) GO TO 50
      WRITE(OUTP,60) LNOW,GJWANT,GJSUM
C-----IF REQUESTED ADD MISSING SEQUENCES
      IF(MYPOT.EQ.0) GO TO 40
      WRITE(OUTP,70)
      IGJBAD = 1
      GJBAD  = GJWANT - GJSUM
      RETURN
C-----OTHERWISE, DO NOT ADD SEQUENCES
   40 WRITE(OUTP,80)
   50 IGJBAD = 0
      GJBAD  = 0.0
      RETURN
   60 FORMAT(1X,78('-')/1X,7('WARNING...'),'WARNING'/
     1 ' FOR L =',I3,' Expect Sum of Statistical Weights (GJ) to Equal'/
     2 ' (2*L + 1) =',F7.3/
     3 ' Found     =',F7.3)
   70 FORMAT(
     4 ' Corrective Action Will be Taken to Correctly Calculate'/
     4 ' the Potential Scattering Cross Section - This Procedure is'/
     5 ' Based on the Decision of the National Nuclear Data Center,'/
     6 ' Brookhaven National Laboratory, Private Communication,'/
     7 ' Charles Dunford, (April 1991)')
   80 FORMAT(
     1 ' No Corrective Action Taken - Based on Input Parameter')
      END
      SUBROUTINE CHECKJ(AJNOW,LNOW,GJ,MYWAY)
C=======================================================================
C
C     DEFINE GJ AND CHECK J VALUE AGAINST TABLE OF LEGAL VALUES.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/JFIX/GJSET(50),GJBAD,TABJ(50),JOK(50),JMISS,NUMJ,IGJBAD
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/FIXPOT/MYPOT
C-----99/06/18 - UPDATED FOR NEW REICH-MOORE CONVENTION
      AJPLUS = DABS(AJNOW)
C-----DEFINE STATISTICAL WEIGHT.
      CALL GJWAIT(SPI,AJPLUS,GJ)
C-----ADD TO SUM OF WEIGHTS.
      CALL SUMJ(GJ,LNOW,2)
C
C     DO NOT CHECK SINGLE LEVEL J VALUE - AVERAGE J VALUE ALLOWED.
C
      IF(MYWAY.NE.2) GO TO 30
C-----CHECK J VALUE.
      DO 10 I=1,NUMJ
      IF(DABS(AJPLUS-TABJ(I)).LE.0.01) GO TO 20
   10 CONTINUE
      WRITE(OUTP,40) LNOW,AJPLUS
      GO TO 30
C-----INDICATE THAT J SEQUENCE HAS BEEN FOUND.
   20 JOK(I)=1
   30 RETURN
   40 FORMAT(1X,78('-')/1X,7('WARNING...'),'WARNING'/
     1 ' Based on Target Spin and L=',I2,' J =',F7.3,
     2 ' is NOT Physically Possible.'/
     3 ' The Use of Fictitious J Values is NOT Allowed in ENDF/B.'/
     5 ' Cross Sections Obtained from these Parameters will be',
     6 ' Unreliable.'/
     7 ' Correct Evaluation Before Attempting Reconstruction.'/
     8 1X,78('-'))
      END
      SUBROUTINE MISSINGJ(LNOW,LSECT)
C=======================================================================
C
C     IF NECESSARY ADD SECTION FOR MISSING J VALUES, ONLY IF MYPOT = 1.
C
C     CHECK SUM OF STATISTICAL WEIGHTS.
C
C     NEW SECTION WILL HAVE MOST PARAMETERS THE SAME AS LSECT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/JFIX/GJSET(50),GJBAD,TABJ(50),JOK(50),JMISS,NUMJ,IGJBAD
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/RANGER/LOW,LHI
      COMMON/FIXPOT/MYPOT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      JMISS=0
      DO 20 J=1,NUMJ
      IF(JOK(J).GT.0) GO TO 20
C-----PRINT TITLE BEFORE FIRST MISSING J VALUE.
      IF(JMISS.GT.0) GO TO 10
      JMISS=1
      WRITE(OUTP,50) LNOW
C-----MISSING J-VALUE. INCREMENT BAD J SEQUENCE COUNT AND ADD
C-----STATISTICAL WEIGHT.
   10 WRITE(OUTP,60) TABJ(J)
   20 CONTINUE
      IF(JMISS.LE.0) GO TO 30
      WRITE(OUTP,70)
   30 CALL SUMJ(GJNOW,LNOW,3)
      IF(IGJBAD.EQ.0.OR.MYPOT.LE.0) GO TO 40
c-----10/05/28 - define correction for missing sequences
      ADDL(LSECT) = 2.0D+0*GJBAD
   40 RETURN
   50 FORMAT(1X,78('-')/1X,7('WARNING...'),'WARNING'/
     1 ' for L=',I2,' The Following J Sequences are Missing.')
   60 FORMAT(' J =',F8.2)
   70 FORMAT(1X,78('-'))
      END
      SUBROUTINE RDAP
C=======================================================================
C
C     READ, WRITE AND SAVE ENERGY DEPENDENT SCATTERING RADIUS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*4 TYPINT
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/LASTE/ELAST
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
C-----2/14/04 - ADD INCLUDE
      INCLUDE 'recent.h'
      DIMENSION TYPINT(4,6)
      DATA TYPINT/
     1 '(His','togr','am) ','    ',
     2 '(Lin',' X-L','in Y',')   ',
     3 '(Log',' X-L','in Y',')   ',
     4 '(Lin',' X-L','og Y',')   ',
     5 '(Log',' X-L','og Y',')   ',
     6 '(ERR','OR) ','    ','    '/
C-----READ TAB1 LEADER LINE AND INTERPOLATION LAW.
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      CALL CARDO(C1,C2,L1,L2,N1,N2)
C-----INCREMENT NUMBER OF ENERGY RANGES WITH ENERGY DEPENDENT SCATTERING
C-----RADIUS AND DEFINE INDICES TO INTERPOLATION LAW AND TABULATED DATA.
      NUMRHO=NUMRHO+1
      IF(NUMRHO.GT.1) GO TO 10
      INXRHO(1,1)=1
      INXRHO(3,1)=1
      GO TO 20
   10 INXRHO(1,NUMRHO)=INXRHO(2,NUMRHO-1)+1
      INXRHO(3,NUMRHO)=INXRHO(4,NUMRHO-1)+1
   20 INXRHO(2,NUMRHO)=INXRHO(1,NUMRHO)+N1-1
      INXRHO(4,NUMRHO)=INXRHO(3,NUMRHO)+N2-1
      INX1=INXRHO(1,NUMRHO)
      INX2=INXRHO(2,NUMRHO)
      INX3=INXRHO(3,NUMRHO)
      INX4=INXRHO(4,NUMRHO)
C-----TEST FOR TABLE OVERFLOW.
      IF(INX2.LE.MAXSEC.AND.INX4.LE.MAXRHO) GO TO 40
      WRITE(OUTP,30) INX2,MAXSEC,INX4,MAXRHO
   30 FORMAT(///' ERROR'/
     1 ' Energy Dependent Scattering Radius Table Overflow'/
     2 ' Interpolation Regions---',I6,' (Cannot Exceed',I6,')'/
     3 ' Tabulated Energies------',I6,' (Cannot Exceed',I6,')'/
     4 ' Increase Dimension in COMMON/TABRHO/ and Re-Run'/
     5 ' Execution Terminated.'///)
      CALL ENDERROR
C-----READ AND CHECK INTERPOLATION LAW.
   40 CALL TERPI(NBTRHO(INX1),INTRHO(INX1),N1)
      CALL TERPO(NBTRHO(INX1),INTRHO(INX1),N1)
C-----PRINT TITLE IN EDIT MODE.
      IF(IMEDIT.GT.0) WRITE(OUTP,160) N1,N2
      IF(N1.LE.0) GO TO 90
      NBTLST=0
      IERAP=0
      DO 80 I=INX1,INX2
      II=INTRHO(I)
      IF(II.LT.1.OR.II.GT.5) GO TO 50
      IF(NBTRHO(I).le.NBTLST) go to 60
      go to 70
   50 II=6
   60 IERAP=1
   70 IF(IMEDIT.GT.0) WRITE(OUTP,170) I,NBTRHO(I),II,
     1 (TYPINT(J,II),J=1,4)
   80 NBTLST=NBTRHO(I)
      IF(IERAP.EQ.0.AND.N2.GT.1.AND.NBTRHO(INX2).EQ.N2) GO TO 120
C-----ERROR IN SCATTERING RADIUS INTERPOLATION LAW.
   90 WRITE(OUTP,100)
      WRITE(*   ,100)
  100 FORMAT(///' ERROR'/
     1 ' Energy Dependent Scattering Radius Interpolation Law Error'/
     3 ' Number of Regions (N1) MUST be Positive'/
     4 ' Number of Energies (N2) MUST be 2 or More'/
     5 ' Interpolation Region Boundaries MUST be in Ascending Order'/
     6 ' Interpolation Law MUST be 1 to 5'/
     7 ' Last Region Boundary MUST be Equal to the Number of Points',
     8 ' (N2)'/' Execution Terminated.'///)
      IF(IMEDIT.EQ.0) then
      WRITE(OUTP,110)
      WRITE(*   ,110)
      ENDIF
  110 FORMAT(' Suggest You Re-Run in Edit Mode to Determine Error')
      CALL ENDERROR
C-----READ AND CHECK TABULATED ENERGY DEPENDENT SCATTERING RADIUS.
  120 ELAST=0.0
      CALL POINTI(ERHOTB(INX3),APTAB(INX3),N2)
      CALL POINTO(ERHOTB(INX3),APTAB(INX3),N2)
      IF(IMEDIT.GT.0) WRITE(OUTP,180)
      IERAP=0
      DO 130 I=INX3,INX4
      IF(APTAB(I).LE.0.0) IERAP=1
      IF(IMEDIT.EQ.0) GO TO 130
      CALL OUT9(ERHOTB(I),FIELD(1,1))
      CALL OUT9(APTAB (I),FIELD(1,2))
      WRITE(OUTP,190) I,((FIELD(M,J),M=1,11),J=1,2)
  130 CONTINUE
      IF(IERAP.EQ.0) GO TO 150
      WRITE(OUTP,140)
      WRITE(*   ,140)
  140 FORMAT(///' ERROR - Energy Dependent Scattering Law Error.'/
     3 '         Energies MUST be Monotonically Increasing.'/
     4 '         Scattering Radius MUST be Positive.'/
     5 '         Execution Terminated.'///)
      CALL ENDERROR
  150 RETURN
  160 FORMAT(1X,78('=')/
     1 ' Energy Dependent Scattering Radius Interpolation Law'/
     2 1X,78('=')/
     3 I5,' Interpolation Ranges'/
     4 I5,' Tabulated Values'/
     5 1X,78('=')/' Interpolation Law'/1X,78('=')/
     6 ' Index  Boundary       Law'/1X,78('='))
  170 FORMAT(I6,2I10,1X,4A4)
  180 FORMAT(1X,78('=')/
     1 ' Energy Dependent Scattering Radius'/1X,78('=')/
     2 ' Index      Energy      Radius'/
     3 '             (eV)      (Fermi)'/1X,78('='))
  190 FORMAT(I6,1X,11A1,1X,11A1)
      END
      SUBROUTINE RDRM1
C=======================================================================
C
C     READ REICH-MOORE DATA FOR ONE ENERGY RANGE.
C     EACH L STATE WILL BE TREATED AS A SEPERATE SECTION.
C
C     THIS ROUTINE USES THE GENERAL REICH-MOORE FORMALISM WITH TWO
C     FISSION CHANNELS AS DEFINED IN ENDF/B-IV. THIS ROUTINE WILL BE
C     USED TO READ DATA IN ANY VERSION OF THE ENDF/B FORMAT (NOTE,
C     THE ENDF/B-VI REICH-MOORE FORMAT HAS NOW BEEN UPDATED TO BE
C     EXACTLY THE SAME AS EARLIER VERSIONS OF THE ENDF/B FORMAT).
C
C     CHECK FOR PRELIMINARY ENDF/B-VI FORMAT(NOW ABANDONED). TERMINATE
C     EXECUTION IF DATA IS IN PRELIMINARY ENDF/B-VI FORMAT.
C
C     FIELD DEFINITIONS FOR EACH RESONANCE
C
C     FIELD          PRELIMINARY          CURRENT
C                    ENDF/B-VI FORMAT     ENDF/B-VI FORMAT
C     =====          ================     ================
C       1            ENERGY               ENERGY
C       2            J                    J
C       3            TOTAL WIDTH          ELASTIC WIDTH
C       4            ELASTIC WIDTH        CAPTURE WIDTH
C       5            CAPTURE WIDTH        FISSION WIDTH 1
C       6            NOT USED             FISSION WIDTH 2
C                    (FISSION NOT
C                    ALLOWED)
C
C     IF THIRD FIELD (PRELIMINARY TOTAL) IS EQUAL TO THE SUM OF THE
C     FOURTH (PRELIMINARY ELASTIC) AND FIFTH (PRELIMINARY CAPTURE)
C     AND SIXTH FIELD IS ZERO (FISSION NOT ALLOWED IN PRELIMINARY
C     FORMAT) FOR ALL RESONANCES THIS PROGRAM WILL ASSUME THAT THE
C     DATA IS IN THE PRELIMINARY REICH-MOORE FORMAT AND TERMINATE
C     EXECUTION WITH A WARNING MESSAGE THAT THE DATA MUST BE CONVERTED
C     TO THE CURRENT REICH-MOORE FORMAT BEFORE IT CAN BE PROCESSED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/NAPRHO/NRO,NAPS
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RANGER/LOW,LHI
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DATA C3/1.23D-01/
      DATA C5/8.0D-02/
C-----UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITTEE RECOMMENDATION
c***** Updated 9/01/04 - ENDF-102, appendix H
c     DATA C4/1.008664904D+00/
      DATA C4/1.00866491578D+00/
c***** Updated 9/01/04
      DATA C6/2.196807122623D-03/
C-----DEFINE CONSTANT TO ALLOW UP TO 1 PER-CENT DIFFERENCE BETWEEN
C-----TOTAL AND ELASTIC PLUS CAPTURE FOR TEST OF PRELIMINARY ENDF/B-VI
C-----FORMAT.
      DATA ALLOW/0.01/
C-----INITIALIZE NEGATIVE WIDTH COUNT
      NEGSUM = 0
C-----INITIALIZE NO CAPTURE
      NOCAPSUM = 0
C-----INITIALIZE FLAG TO INDICATE THAT DATA IS IN THE PRELIMINARY
C-----ENDF/B-VI FORMAT.
      IMOLD=1
C-----DEFINE TARGET SPIN, SCATTERING RADIUS AND NUMBER OF L STATES.
      CALL CARDIO(SPI,AP,L1,L2,NLS,N2)
      CALL OUT9(SPI,FIELD(1,1))
      CALL OUT9(AP ,FIELD(1,2))
      WRITE(OUTP,180) ((FIELD(M,I),M=1,11),I=1,2),NLS
      WRITE(*   ,180) ((FIELD(M,I),M=1,11),I=1,2),NLS
C
C     L STATE LOOP.
C
      DO 130 IL=1,NLS
C-----DEFINE ATOMIC WEIGHT RATIO, L DEPENDENT SCATTERING RADIUS, L VALUE
C-----AND NUMBER OF RESONANCES.
      CALL CARDIO(AWRI,APL,LNOW,L2,NRS6X,NRS)
      AWRICM = AWRI/(AWRI+ONE)
      CALL OUT9(AWRI,FIELD(1,1))
      CALL OUT9(APL ,FIELD(1,2))
      WRITE(OUTP,190) ((FIELD(M,I),M=1,11),I=1,2),LNOW,NRS
      WRITE(*   ,190) ((FIELD(M,I),M=1,11),I=1,2),LNOW,NRS
C-----PRINT WARNING IF L DEPENDENT SCATTERING RADIUS IS NOT POSITIVE.
      IF(APL.GT.0.0) GO TO 10
C-----IF L DEPENDENT RADIUS IS NOT POSITIVE USE SCATTERING RADIUS.
C-----(SEE, ENDF/B-VI FORMATS AND PROCEDURES MANUAL).
      APL=AP
      CALL OUT9(APL ,FIELD(1,2))
      WRITE(OUTP,200) (FIELD(KK,2),KK=1,11)
      WRITE(*   ,200) (FIELD(KK,2),KK=1,11)
C
C     CONSTRUCT TABLE OF PERMISSIBLE J-VALUES
C     (THIS TABLE WILL BE USED TO CHECK FOR MISSING OR ILLEGAL J-VALUES)
C
   10 CALL SETJ(LNOW)
C-----INCREMENT SECTION COUNT AND INDICES TO RESONANCE PARAMETER TABLE.
      CALL LIMIT1(1)
C-----DEFINE ALL PARAMETERS FOR A SECTION.
C-----EITHER CALCULATE CHANNEL RADIUS OR DEFINE TO BE EQUAL TO THE
C-----SCATTERING RADIUS.
      IF(NAPS.NE.0) GO TO 20
      APX=C3*((C4*AWRI)**(ONE/THREE))+C5
      GO TO 30
   20 APX=APL
C
C     APX = IS USED IN PENETRABILITIES AND SHIFT FACTORS.
C         = EITHER BASED ON FORMULA OR EQUAL TO APL, DEPENDING ON NAPS.
C     APL = IS USED IN THE HARD SPHERE PHASE SHIFTS.
C     (SEE, ENDF/B-VI FORMATS AND PROCEDURES MANUAL, PAGE 2.6)
C
   30 AK1=C6*AWRICM
      AK2=AK1*AK1
      BETA(NSECT)=PI*ABN/AK2
      RHOX2(NSECT)=APX*APX*AK2
      RHOP1(NSECT)=APL*AK1
      EL(NSECT)=ELX
      EH(NSECT)=EHX
      NLOW(NSECT)=LOW
      NHIGH(NSECT)=LHI
      LVALUE(NSECT)=LNOW
C-----Added for L-State Dependent fission.
      LFWL(NSECT)  = 0
C-----10/05/28 - added for potential correction
      POTL(NSECT) = 4*LNOW + 2  ! 2(2l+1) total weight of all sequences
      ADDL(NSECT) = 0.0D+0      ! initialize weight of missing sequences
c----- 02/08/08 - ADDED NAPS DEFINITION
      NAPTAB(NSECT)=NAPS
      MODE(NSECT)=LRF
C
C     06/09/05 - THIS WAS ALREADY DONE FOR BREIT-WIGNER - NOW ALL.
C
C     IF ENERGY DEPENDENT SCATTERING RADIUS DEFINE INDEX TO TABULATED
C     DATA (OTHERWISE INDEX HAS ALREADY BEEN INITIALIZED TO 0).
C     CALCULATE ENERGY DEPENDENT RHOP FROM SCATTERING RADII.
C
      IF(NRO.GT.0) THEN
      NRHO(NSECT)=NUMRHO
      INX3=INXRHO(3,NUMRHO)
      INX4=INXRHO(4,NUMRHO)
      DO I=INX3,INX4
      RHOTAB(I)=APTAB(I)*AK1
      ENDDO
      ENDIF
C
C     READ RESONANCE PARAMETERS (6 PER RESONANCE).
C
C     (1) ENERGY
C     (2) J
C     (3) NEUTRON WIDTH
C     (4) CAPTURE WIDTH
C     (5) FIRST FISSION WIDTH
C     (6) SECOND FISSION WIDTH
C
C     SOME OF THESE PARAMETERS WILL CONVERTED TO THE FORM REQUIRED
C     FOR CALCULATIONS.
C
C     (2) STATISTICAL WEIGHT, GJ (J IS NOT NEEDED)
C     (3) 1/2 NEUTRON WIDTH DIVIDED BY PENETRATION FACTOR (NEUTRON WIDTH
C         IS ONLY USED DIVIDED BY THE PENETRATION FACTOR)
C     (4) 1/2 THE CAPTURE WIDTH (THIS IS THE ONLY FORM IN WHICH THE
C         CAPTURE WIDTH IS USED)
C     (5) SIGNED SQUARE ROOT OF 1/2 FIRST FISSION WIDTH (ONLY FORM USED)
C     (6) SIGNED SQUARE ROOT OF 1/2 SECOND FISSION WIDTH (ONLY FORM USED
C
      CALL LISPAR6(LOW,NRS6X)
      IF(IMEDIT.NE.0) WRITE(OUTP,210)
C-----FOR EACH SECTION (ISOTOPE, ENERGY RANGE, L VALUE) SORT
C-----RESONANCES IN ASCENDING (J, E) ORDER.
      CALL ORDER(LOW,LHI)
      AJNOW=RESTAB(2,LOW)
C-----CHECK FOR LEGAL J-VALUE AND DEFINE STATISTICAL WEIGHT.
      CALL CHECKJ(AJNOW,LNOW,GJ,2)
      DO 110 JR=LOW,LHI
C-----INSERT DIVIDING LINE BETWEEN DIFFERENT J VALUES.
      IF(DABS(AJNOW-RESTAB(2,JR)).LE.0.01) GO TO 40
      AJNOW=RESTAB(2,JR)
C-----CHECK FOR LEGAL J-VALUE AND DEFINE STATISTICAL WEIGHT.
      CALL CHECKJ(AJNOW,LNOW,GJ,2)
      IF(IMEDIT.NE.0) WRITE(OUTP,230)
C-----CHECK FOR NEGATIVE WIDTHS OR NO CAPTURE
   40 NEGRES = 0
      NOCAPT = 0
      DO I=3,4
      IF(RESTAB(I,JR).LT.0.0D+0) NEGRES = NEGRES + 1
      ENDDO
      IF(RESTAB(4,JR).LE.0.0D+0) NOCAPT = 1
      IF(IMEDIT.EQ.0.AND.NEGRES.EQ.0) GO TO 80
C-----LIST RESONANCE PARAMETERS.
      DO 50 I=1,6
   50 CALL OUT9(RESTAB(I,JR),FIELD(1,I))
      WRITE(OUTP,220) (FIELD(M,1),M=1,11),RESTAB(2,JR),
     1 ((FIELD(M,I),M=1,11),I=3,6)
C-----NEGATIVE WIDTHS?
      IF(NEGRES.GT.0) THEN
      NEGSUM = NEGSUM + NEGRES
      WRITE(OUTP,220) (FIELD(M,1),M=1,11),RESTAB(2,JR),
     1 ((FIELD(M,I),M=1,11),I=3,6)
      WRITE(   *,60) NEGRES
      WRITE(OUTP,60) NEGRES
   60 FORMAT(' ERROR - ',I1,' Negative Widths')
      ENDIF
C-----NO CAPTURE?
      IF(NOCAPT.GT.0) THEN
      NOCAPSUM = NOCAPSUM + NOCAPT
      WRITE(OUTP,220) (FIELD(M,1),M=1,11),RESTAB(2,JR),
     1 ((FIELD(M,I),M=1,11),I=3,6)
      WRITE(   *,70)
      WRITE(OUTP,70)
   70 FORMAT(' ERROR - Capture MUST be Positive')
      ENDIF
C-----CHECK FOR PRELIMINARY ENDF/B-VI FORMAT.
   80 IF(IMOLD.LE.0) GO TO 100
C-----NOT PRELIMINARY FORMAT IF THIRD FIELD IS NOT THE SUM OF THE
C-----FOURTH AND FIFTH FIELDS.
      IF(DABS(RESTAB(3,JR)-(RESTAB(4,JR)+RESTAB(5,JR))).GT.
     1 DABS(ALLOW*RESTAB(3,JR))) GO TO 90
C-----NOT PRELIMINARY FORMAT IF SIXTH FIELD IS NOT ZERO.
      IF(DABS(RESTAB(6,JR)).EQ.0.0) GO TO 100
C-----DATA IS NOT IN THE PRELIMINARY ENDF/B-VI FORMAT. TURN OFF FLAG
C-----AND STOP TESTING FOR PRELIMINARY FORMAT.
   90 IMOLD=0
C-----REPLACE J VALUE BY STATISTICAL WEIGHT.
  100 RESJTAB(JR) = RESTAB(2,JR)
      RESTAB(2,JR)= GJ
C-----DEFINE SIGNED SQUARE ROOT OF 1/2 FISSION WIDTHS.
      GAMF1=DSQRT(DABS(0.5d0*RESTAB(5,JR)))
      IF(RESTAB(5,JR).LT.0.0) GAMF1=-GAMF1
      RESTAB(5,JR)=GAMF1
      GAMF2=DSQRT(DABS(0.5d0*RESTAB(6,JR)))
      IF(RESTAB(6,JR).LT.0.0) GAMF2=-GAMF2
      RESTAB(6,JR)=GAMF2
C-----IF FISSION WIDTHS ARE NOT ZERO TURN ON FISSILE FLAG.
      IF(GAMF1.NE.0.0.OR.GAMF2.NE.0.0) LFWX=1
      IF(GAMF1.NE.0.0.OR.GAMF2.NE.0.0) LFWL(NSECT) = LFWL(NSECT) + 1
      IF(GAMF1.NE.0.0.and.GAMF2.NE.0.0) IMBOTH = IMBOTH + 1
C-----DEFINE 1/2 NEUTRON WIDTH DIVIDED BY PENETRATION FACTOR.
C-----04/21/07 - IF NECESSARY DEFINE RHOX2
      ERABS=DABS(RESTAB(1,JR))
      IF(NRO.EQ.1.AND.NAPS.EQ.1) CALL SETRHO1(ERABS,NSECT)
      RHO2X=ERABS*RHOX2(NSECT)
      CALL FACTS3(LVALUE(NSECT),RHO2X,PENFAC)
C-----02/14/04 - PROTECT AGAINST ZERO ENERGY RESONANCES
      IF(PENFAC.NE.0.0D+00) THEN
      RESTAB(3,JR)=0.5d0*RESTAB(3,JR)/PENFAC
      ELSE
      RESTAB(3,JR)=0.0D+00
      ENDIF
C-----DEFINE 1/2 CAPTURE WIDTH.
      RESTAB(4,JR)=0.5d0*RESTAB(4,JR)
  110 CONTINUE
C
C     PRINT WARNING IF ANY J SEQUENCES ARE MISSING AND ADD MISSING J
C     SEQUENCE WITH NO RESONANCES IN ORDER TO ALLOW POTENTIAL SCATTERING
C     TO BE CORRECTLY CALCULATED.
C
      LSECT = NSECT
      CALL MISSINGJ(LNOW,LSECT)
c
c     2016/11/19 - Added Summary
c
      iii = (NHIGH(NSECT) - NLOW(NSECT)) + 1
      write(OUTP,120) NSECT,LVALUE(NSECT),iii,LFWL(NSECT),IMBOTH
      write(   *,120) NSECT,LVALUE(NSECT),iii,LFWL(NSECT),IMBOTH
  120 FORMAT(1X,78('-')/' Summary of the Above Section'/1x,78('-')/
     1       ' Section Number.........--------------',I11/
     2       ' L Value................--------------',I11/
     3       ' Total Number of Resonances-----------',I11/
     4       ' Number With Non-Zero Fission Width---',I11/
     5       ' Number With 2 Non-Zero Fission Widths',I11)
c
c     2016/11/19 - Added Summar
c
C-----END OF L STATE LOOP.
  130 CONTINUE
C-----TERMINATE EXECUTION IF DATA IS IN THE PRELIMINARY ENDF/B-VI
C-----FORMAT.
      IF(IMOLD.LE.0) GO TO 150
      WRITE(OUTP,140)
      WRITE(*   ,140)
  140 FORMAT(///' ERROR'/
     1 ' Data are in the Preliminary ENDF/B-VI Format, which is No'/
     2 ' Longer Used. Data MUST be Reformated to Comform to the'/
     3 ' Current ENDF/B-VI Reich-Moore Format Before it can be'/
     4 ' Processed by this Program. Execution Terminated.'///)
      CALL ENDERROR
C
C     STOP IF NEGATIVE WIDTHS
C
  150 IF(NEGSUM.GT.0) THEN
      WRITE(   *,160) NEGSUM
      WRITE(OUTP,160) NEGSUM
  160 FORMAT(///
     1 ' ERROR - Execution Terminated because of',i8,
     2 ' Negative Widths.'///)
      CALL ENDERROR
      ENDIF
C
C     STOP IF NO CAPTURE
C
      IF(NOCAPSUM.GT.0) THEN
      WRITE(   *,170) NOCAPSUM
      WRITE(OUTP,170) NOCAPSUM
  170 FORMAT(///
     1 ' Execution Terminated because of No Capture in',i8,
     2 ' Resonances.'///)
      CALL ENDERROR
      ENDIF
      RETURN
  180 FORMAT(1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Scattering Radius--------------------',11A1/
     3       ' Number of L Values-------------------',I11)
  190 FORMAT(1X,78('=')/
     1       ' Atomic Weight Ratio of Isotope-------',11A1/
     2       ' L Dependent Scattering Radius--------',11A1/
     3       ' Angular Momentum (L)-----------------',I11/
     4       ' Number of Resonances-----------------',I11)
  200 FORMAT(1X,7('WARNING...'),'WARNING'/
     1 ' L Dependent Scattering Radius in the Evaluation is Zero.'/
     2 ' Have Defined it to be Equal to the Scattering Radius',11A1/
     3 ' (see, ENDF/B-VI Formats and Procedures Manual, page 2.11)')
  210 FORMAT(1X,78('=')/' Reich-Moore Resonance Parameters'/1X,78('=')/
     1 '      Energy  J Value    Neutron    Capture  Fission-1',
     2 '  Fission-2'/8X,'(eV)',9X,4(7X,'(eV)')/1X,78('='))
  220 FORMAT(1X,11A1,F7.2,2X,44A1)
  230 FORMAT(1X,78('-'))
      END
      SUBROUTINE RDAA
C=======================================================================
C
C     READ ADLER-ADLER DATA FOR ONE ENERGY RANGE.
C     SINCE ADLER-ADLER PARAMETERS ARE INDEPENDENT OF L AND J THE
C     ENTIRE ENERGY RANGE WILL BE TREATED AS ONE SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*8 REACTR
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RANGER/LOW,LHI
      COMMON/WHATZA/IZANOW,MATNOW,TEMP3,IVERSE,INT45
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
      COMMON/NAPRHO/NRO,NAPS
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION REACTR(3)
      DATA REACTR/
     1 'Total   ',
     2 'Fission ',
     3 'Capture '/
C-----UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITTEE RECOMMENDATION
      DATA C6/2.196807122623D-03/
C-----DEFINE TARGET SPIN, SPIN UP SCATTERING RADIUS AND NUMBER OF L
C-----STATES.
      CALL CARDIO(SPI,AP,L1,L2,NLS,N2)
      CALL OUT9(SPI,FIELD(1,1))
      CALL OUT9(AP ,FIELD(1,2))
      WRITE(OUTP,120) ((FIELD(M,I),M=1,11),I=1,2),NLS
      WRITE(*   ,120) ((FIELD(M,I),M=1,11),I=1,2),NLS
C-----DEFINE ATOMIC WEIGHT RATIO, TYPE AND NUMBER OF BACKGROUND
C-----PARAMETERS.
      CALL CARDIO(AWRI,XN,LI,L2,NRS6X,NRS)
      CALL OUT9(AWRI,FIELD(1,1))
      WRITE(OUTP,130) (FIELD(M,1),M=1,11),LI,NRS
      WRITE(*   ,130) (FIELD(M,1),M=1,11),LI,NRS
C-----INCREMENT SECTION COUNT AND INDICES TO RESONANCE PARAMETER TABLE
C-----(RESONANCE TABLE INDICES WILL ONLY REFER TO BACKGROUND).
      CALL LIMIT1(1)
C-----DEFINE ALL PARAMETERS FOR A SECTION.
      AK1=C6*AWRI/(AWRI+ONE)
      BETA(NSECT)=PI*ABN/(AK1*AK1)
      RHOP1(NSECT)=TWO*AP*AK1
      EL(NSECT)=ELX
      EH(NSECT)=EHX
      NLOW(NSECT)=LOW
C-----MODE DEPENDS ON VERSION OF ENDF/B FORMAT.
      MODE(NSECT)=LRF
C
C     06/09/05 - THIS WAS ALREADY DONE FOR BREIT-WIGNER - NOW ALL.
C
C     IF ENERGY DEPENDENT SCATTERING RADIUS DEFINE INDEX TO TABULATED
C     DATA (OTHERWISE INDEX HAS ALREADY BEEN INITIALIZED TO 0).
C     CALCULATE ENERGY DEPENDENT RHOP FROM SCATTERING RADII.
C
      IF(NRO.GT.0) THEN
      NRHO(NSECT)=NUMRHO
      INX3=INXRHO(3,NUMRHO)
      INX4=INXRHO(4,NUMRHO)
      DO I=INX3,INX4
      RHOTAB(I)=APTAB(I)*AK1
      ENDDO
      ENDIF
      IF(IVERSE.LT.5) MODE(NSECT)=8
      LVALUE(NSECT)=LI
c----- 02/08/08 - ADDED NAPS DEFINITION
      NAPTAB(NSECT)=NAPS
C
C     READ BACKGROUND CROSS SECTIONS.
C
C     LI = 1 = TOTAL
C        = 2 = FISSION
C        = 3 = TOTAL + FISSION
C        = 4 = CAPTURE
C        = 5 = TOTAL + CAPTURE
C        = 6 = FISSION + CAPTURE
C        = 7 = ALL 3
C
C     FOR CONVENIENCE ALWAYS ALLOW ROOM FOR 3 SETS OF BACKGROUND
C     CROSS SECTIONS AND POSITION BACKGROUND CROSS SECTIONS TO
C     STANDARD POSITIONS (TOTAL, FISSION AND CAPTURE, IN THAT ORDER).
C
C-----INIIALIZE ALL OFF.
      LHI=LOW+2
      DO I=LOW,LHI
      DO K=1,6
      RESTAB(K,I) = 0.0
      ENDDO
      ENDDO
      IF(NRS.LE.0) GO TO 80
      GO TO (10,20,30,40,50,60,70),LI
   10 CALL LISTIO(RESTAB(1,LOW  ),6) ! TOTAL
      GO TO 80
   20 CALL LISTIO(RESTAB(1,LOW+1),6) ! FISSION
      GO TO 80
   30 CALL LISTIO(RESTAB(1,LOW  ),6) ! TOTAL
      CALL LISTIO(RESTAB(1,LOW+1),6) ! FISSION
      GO TO 80
   40 CALL LISTIO(RESTAB(1,LOW+2),6) ! CAPTURE
      GO TO 80
   50 CALL LISTIO(RESTAB(1,LOW  ),6) ! TOTAL
      CALL LISTIO(RESTAB(1,LOW+2),6) ! CAPTURE
      GO TO 80
   60 CALL LISTIO(RESTAB(1,LOW+1),6) ! FISSION
      CALL LISTIO(RESTAB(1,LOW+2),6) ! CAPTURE
      GO TO 80
C-----TOTAL, FISSION AND CAPTURE. NOTHING TO D0.
   70 CALL LISTIO(RESTAB(1,LOW  ),6) ! TOTAL
      CALL LISTIO(RESTAB(1,LOW+1),6) ! FISSION
      CALL LISTIO(RESTAB(1,LOW+2),6) ! CAPTURE
C-----SET FISSION FLAG IF ANY FISSION BACKGROUND
   80 DO I=1,6
      IF(DABS(RESTAB(I,LOW+1)).NE.0.0) LFWX=1
      IF(DABS(RESTAB(I,LOW+1)).NE.0.0) LFWL(NSECT) = LFWL(NSECT) + 1
      ENDDO
C
C     LIST BACKGROUND CROSS SECTIONS IN STANDARD POSITION.
C
      IF(IMEDIT.LE.0) GO TO 90
      WRITE(OUTP,140)
C-----TOTAL
      DO I=1,6
      CALL OUT9(RESTAB(I,LOW  ),FIELD(1,I))
      ENDDO
      WRITE(OUTP,150) REACTR(1),((FIELD(M,I),M=1,11),I=1,6)
C-----FISSION
      DO I=1,6
      CALL OUT9(RESTAB(I,LOW+1),FIELD(1,I))
      ENDDO
      WRITE(OUTP,150) REACTR(2),((FIELD(M,I),M=1,11),I=1,6)
C-----CAPTURE
      DO I=1,6
      CALL OUT9(RESTAB(I,LOW+2),FIELD(1,I))
      ENDDO
      WRITE(OUTP,150) REACTR(3),((FIELD(M,I),M=1,11),I=1,6)
C
C     L STATE LOOP.
C
   90 DO 110 IL=1,NLS
C-----READ L VALUE AND NUMBER OF J STATES.
      CALL CARDIO(XN,XN,LNOW,L2,NJS,N2)
      WRITE(OUTP,190) LNOW,NJS
C
C     J STATE LOOP.
C
      DO 110 IJ=1,NJS
C-----READ J VALUE AND NUMBER OF RESONANCES
      CALL CARDIO(AJ,C2,L1,L2,NRS12,NRS)
      WRITE(OUTP,200) AJ,NRS
C-----INCREMENT INDICES TO RESONANCE PARAMETER TABLE.
      CALL LIMIT1(2)
C-----READ RESONANCE PARAMETERS - 12 PARAMETERS PER RESONANCE.
      CALL LISPAR12(LOW,NRS12)
      IF(IMEDIT.NE.0) WRITE(OUTP,160)
      DO 100 JR=LOW,LHI
      IF(DABS(RESTAB(7,JR)).NE.0.0.OR.DABS(RESTAB(8,JR)).NE.0.0)
     1 LFWX = 1
      IF(DABS(RESTAB(7,JR)).NE.0.0.OR.DABS(RESTAB(8,JR)).NE.0.0)
     1 LFWL(NSECT) = 1
C
C     LIST RESONANCE PARAMETERS - NOW STORED TO 12 WORDS PER RESONANCE
C
      IF(IMEDIT.EQ.0) GO TO 100
C-----TOTAL
      DO I=1,4
      CALL OUT9(RESTAB(I  ,JR),FIELD(1,I))
      ENDDO
      WRITE(OUTP,170) REACTR(1),((FIELD(M,I),M=1,11),I=1,4)
C-----FISSION
      DO I=1,4
      CALL OUT9(RESTAB(I+4,JR),FIELD(1,I))
      ENDDO
      WRITE(OUTP,170) REACTR(2),((FIELD(M,I),M=1,11),I=1,4)
C-----CAPTURE
      DO I=1,4
      CALL OUT9(RESTAB(I+8,JR),FIELD(1,I))
      ENDDO
      WRITE(OUTP,170) REACTR(3),((FIELD(M,I),M=1,11),I=1,4)
      WRITE(OUTP,180)
  100 CONTINUE
C-----END OF L AND J STATE LOOPS.
  110 CONTINUE
C-----DEFINE UPPER INDEX TO RESONANCE TABLE.
      NHIGH(NSECT)=LHI
      RETURN
  120 FORMAT(1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Effective Scattering Radius (A+)-----',11A1/
     3       ' Number of L Values-------------------',I11)
  130 FORMAT(1X,78('=')/
     1       ' Atomic Weight Ratio of Isotope-------',11A1/
     2       ' LI-----------------------------------',I11/
     3       ' Number of Background Sets------------',I11)
  140 FORMAT(1X,78('=')/' Adler-Adler Resonance Parameters'/1X,78('=')/
     1 ' Background'/1X,78('=')/' Reaction',
     1 '         A1         A2         A3         A4         B1',
     2 '         B2'/1X,78('='))
  150 FORMAT(1X,A8,66A1)
  160 FORMAT(1X,78('=')/' Resonance Parameters'/1X,78('=')/' Reaction',
     1  '  Resonance Half-Width  Symmetric Asymmetric'/
     2  9X,'   Energy     (eV)      Parameter Parameter'/
     3  9X,'    (eV)                  (GRT)     (GIT)  '/1X,78('='))
  170 FORMAT(1X,A8,44A1)
  180 FORMAT(1X,78('='))
  190 FORMAT(1X,78('=')/
     1       ' L Value------------------------------',I11/
     2       ' Number of J Values-------------------',I11)
  200 FORMAT(' J Value------------------------------',5X,F6.2/
     3       ' Number of Resonances-----------------',I11)
      END
      SUBROUTINE RDHRF
C=======================================================================
C
C     READ HYBRID R-FUNCTION DATA FOR ONE ENERGY RANGE.
C     EACH (L,S,J) STATE WILL BE TREATED AS A SEPARATE SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*4 MTHOL,ANSWER
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/NAPRHO/NRO,NAPS
      COMMON/RANGER/LOW,LHI
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/HRFTAB/QHRF(4),MTHRF(4),NHRF
      COMMON/FIELDC/FIELD(11,12)
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION MTOK(9),MTHOL(3,10),ANSWER(3,2)
      DATA ANSWER/' (No',')   ',
     1 ' (Ye','s)  ',
     2 ' (ER','ROR)'/
C-----DEFINE ALLOWED INELASTIC AND CHARGED PARTICLE MT NUMBERS.
      DATA MTOK/51,52,53,54,103,104,105,106,107/
      DATA MTHOL/
     1 'Inel','asti','c   ',
     2 'Inel','asti','c   ',
     3 'Inel','asti','c   ',
     4 'Inel','asti','c   ',
     5 'n,p ','    ','    ',
     6 'n,d ','    ','    ',
     7 'n,t ','    ','    ',
     8 'n,He','-3  ','    ',
     9 'n,Al','pha ','    ','NOT ','Allo','wed '/
      DATA C3/1.23D-01/
      DATA C5/8.0D-02/
C-----UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITTEE RECOMMENDATION
c***** Updated 9/01/04 - ENDF-102, appendix H
      DATA C4/1.00866491578D+00/
c***** Updated 9/01/04
      DATA C6/2.196807122623D-03/
C-----DEFINE TARGET SPIN, SPIN UP SCATTERING RADIUS AND NUMBER OF L
C-----STATES.
      CALL CARDIO(SPI,C2,L1,L2,NLS,N2)
      CALL OUT9(SPI,FIELD(1,1))
      WRITE(OUTP,230) (FIELD(M,1),M=1,11),NLS
      WRITE(*   ,230) (FIELD(M,1),M=1,11),NLS
C-----DEFINE NUMBER OF EACH KIND OF REACTION.
      CALL CARDIO(C1,C2,NGRE,NFRE,NIRE,NCRE)
      WRITE(OUTP,240) NGRE,NFRE,NIRE,NCRE
      WRITE(*   ,240) NGRE,NFRE,NIRE,NCRE
C-----CHECK FOR ALLOWABLE VALUES.
      IF(NGRE.LT.0.OR.NGRE.GT.1) GO TO 10
      IF(NFRE.LT.0.OR.NFRE.GT.1) GO TO 10
      IF(NIRE.LT.0.OR.NIRE.GT.4) GO TO 10
      IF(NCRE.LT.0.OR.NCRE.GT.4) GO TO 10
      NHRF=NIRE+NCRE
      IF(NHRF.LE.4) GO TO 30
   10 WRITE(OUTP,20) NHRF
      WRITE(*   ,20) NHRF
   20 FORMAT(' ERROR - Illegal Number of Reactions',I3,' (Must be 4)'/
     1       '         Execution Terminated.'///)
      CALL ENDERROR
C-----DEFINE KIND OF REACTION.
   30 CALL CARDIO(C1,C2,MTHRF(1),MTHRF(2),MTHRF(3),MTHRF(4))
C-----READ Q-VALUES.
      CALL CARDIO(C1,C2,L1,L2,N1,N2)
      CALL LISTIO(QHRF,N1)
C-----LIST MT NUMBERS AND Q-VALUES.
      IF(NHRF.LE.0) GO TO 100
      IF(IMEDIT.NE.0) WRITE(OUTP,250)
C-----CHECK AND LIST MT NUMBERS AND Q-VALUES.
      IERR=0
      DO 60 I=1,NHRF
      DO 40 J=1,9
      IF(MTHRF(I).EQ.MTOK(J)) GO TO 50
   40 CONTINUE
      J=10
      IERR=1
C-----TEMPORARILY ONLY ALLOW INELASTIC REACTIONS.
   50 IF(J.GT.4.AND.J.LE.9) IERR=2
      CALL OUT9(QHRF(I),FIELD(1,1))
      IF(IMEDIT.NE.0) WRITE(OUTP,260) MTHRF(I),
     1 (FIELD(M,1),M=1,11),(MTHOL(K,J),K=1,3)
C-----FROM HERE ON TREAT ALL INELASTIC AS THE SAME MT NUMBER (MT=51).
      IF(MTHRF(I).GE.51.AND.MTHRF(I).LE.54) MTHRF(I)=51
   60 CONTINUE
C
C     TERMINATED EXECUTION IF REACTION IS NOT ALLOWED.
C
      IF(IERR.EQ.0) GO TO 100
      IF(IERR.EQ.2) GO TO 80
C-----ILLEGAL MT NUMBER.
      WRITE(OUTP,70)
      WRITE(*   ,70)
   70 FORMAT(///' ERROR - Illegal MT Number. Execution Terminated.'///)
      CALL ENDERROR
C-----CHARGED PARTICLE REACTION (TEMPORARILY NOT ALLOWED)
   80 WRITE(OUTP,90)
      WRITE(*   ,90)
   90 FORMAT(///' ERROR -',
     1 ' Currently Program Does NOT Allow Charged Particle Reactions.'/
     2 ' Execution Terminated.'///)
      CALL ENDERROR
C
C     L STATE LOOP.
C
  100 DO 220 IL=1,NLS
C-----DEFINE ATOMIC WEIGHT RATIO, L VALUE AND NUMBER OF S-VALUES
C-----(CHANNEL SPINS).
      CALL CARDIO(AWRI,C1,LNOW,L2,NSS,N2)
C-----DEFINE ATOMIC WEIGHT DEPENDENT QUANTITIES.
      AWRICM=AWRI/(AWRI+ONE)
C-----INITIALIZE CHANNEL RADIUS TO CALCULATED VALUE (MAY BE CHANGED
C-----LATER AFTER SCATTERING RADIUS IS READ).
      APX=C3*((C4*AWRI)**(ONE/THREE))+C5
      AK1=C6*AWRICM
      AK2=AK1*AK1
      CALL OUT9(AWRI,FIELD(1,1))
      WRITE(OUTP,270) (FIELD(M,1),M=1,11),LNOW,NSS
C
C     CONSTRUCT TABLE OF PERMISSIBLE J-VALUES
C     (THIS TABLE WILL BE USED TO CHECK FOR MISSING OR ILLEGAL J-VALUES)
C
      CALL SETJ(LNOW)
C
C     S-VALUE LOOP.
C
      DO 210 IS=1,NSS
C-----DEFINE S-VALUE AND NUMBER OF J-VALUES.
      CALL CARDIO(AS,C2,L1,L2,NJS,N2)
      CALL OUT9(AS,FIELD(1,1))
      WRITE(OUTP,280) (FIELD(M,1),M=1,11),NJS
C
C     J STATE LOOP.
C
      DO 210 IJ=1,NJS
C-----DEFINE J-VALUE, CHANNEL RADIUS, FLAGS FOR BACKGROUND AND OPTICAL
C-----MODEL PHASE SHIFTS AND THE NUMBER OF RESONANCES.
      CALL CARDIO(AJ,AC,LBK,LPS,NLSJ12,NLSJ)
C-----CHECK FOR LEGAL J-VALUE AND DEFINE STATISTICAL WEIGHT.
      CALL CHECKJ(AJ,LNOW,GJ,2)
C-----IF REQUESTED SET CHANNEL RADIUS EQUAL TO SCATTERING RADIUS.
      IF(NAPS.NE.0) APX=AC
      LBKIN=LBK+1
      IF(LBK.LT.0.OR.LBK.GT.1) LBKIN=3
      LPSIN=LPS+1
      IF(LPS.LT.0.OR.LPS.GT.1) LPSIN=3
      CALL OUT9(AJ,FIELD(1,1))
      CALL OUT9(AC,FIELD(1,2))
      WRITE(OUTP,290) ((FIELD(M,I),M=1,11),I=1,2),LBK,
     1 ANSWER(1,LBKIN),ANSWER(2,LBKIN),LPS,ANSWER(1,LPSIN),
     2 ANSWER(2,LPSIN),NLSJ
C-----CHECK FOR ILLEGAL BACKGROUND OR OPTICAL MODEL PHASE SHIFT FLAG.
      IF(LBKIN.NE.2.AND.LPSIN.NE.2) GO TO 120
      WRITE(OUTP,110)
      WRITE(*   ,110)
  110 FORMAT(///' ERROR -',
     1 ' Background and Optical Phase Shift Flags MUST be 0 or 1.'/
     2 ' Execution Terminated.'///)
      CALL ENDERROR
C-----TEMPORARILY DO NOT ALLOW TABULATED BACKGROUND OR OPTICAL MODEL
C-----PHASE SHIFT.
  120 IF(LBK.EQ.0.AND.LPS.EQ.0) GO TO 140
      WRITE(OUTP,130)
c-----207/9/6 - Added missing format # 130
      WRITE(*   ,130)
  130 FORMAT(///' ERROR -',
     1 ' Currently Program Does NOT Allow Tabulated Background or'/
     2 ' Optical Phase Shift. Execution Terminated.'///)
      CALL ENDERROR
C-----ALLOW SPACE FOR 12 PARAMETERS PER RESONANCE.
  140 NRS=  NLSJ
C-----INCREMENT SECTION COUNT AND INDICES TO RESONANCE PARAMETER TABLE.
      CALL LIMIT1(1)
C-----DEFINE ALL PARAMETERS FOR A SECTION.
      BETA(NSECT)=PI*ABN/AK2
      RHOX2(NSECT)=APX*APX*AK2
      RHOP1(NSECT)=AC*AK1
      EL(NSECT)=ELX
      EH(NSECT)=EHX
      GJTAB(NSECT)=GJ
      NLOW(NSECT)=LOW
      NHIGH(NSECT)=LHI
      LVALUE(NSECT)=LNOW
      LRXTAB(NSECT)=NHRF
c----- 02/08/08 - ADDED NAPS DEFINITION
      NAPTAB(NSECT)=NAPS
      MODE(NSECT)=LRF
C
C     06/09/05 - THIS WAS ALREADY DONE FOR BREIT-WIGNER - NOW ALL.
C
C     IF ENERGY DEPENDENT SCATTERING RADIUS DEFINE INDEX TO TABULATED
C     DATA (OTHERWISE INDEX HAS ALREADY BEEN INITIALIZED TO 0).
C     CALCULATE ENERGY DEPENDENT RHOP FROM SCATTERING RADII.
C
      IF(NRO.GT.0) THEN
      NRHO(NSECT)=NUMRHO
      INX3=INXRHO(3,NUMRHO)
      INX4=INXRHO(4,NUMRHO)
      DO I=INX3,INX4
      RHOTAB(I)=APTAB(I)*AK1
      ENDDO
      ENDIF
C-----DEFINE EXCITATION ENERGIES.
      IF(NHRF.LE.0) GO TO 160
      DO 150 I=1,NHRF
  150 EXCITE(I,NSECT)=QHRF(I)/AWRICM
C
C     READ RESONANCE PARAMETERS (12 PER RESONANCE).
C
C     (1) ENERGY
C     (2) NEUTRON WIDTH
C     (3) CAPTURE WIDTH
C     (4) FISSION WIDTH
C     (5) FIRST COMPETITIVE WIDTH
C     (6) SECOND COMPETITIVE WIDTH
C     (7) THIRD COMPETITIVE WIDTH
C     (8) FOURTH COMPETITIVE WIDTH
C     (9) FIRST EXIT CHANNEL L-VALUE
C    (10) SECOND EXIT CHANNEL L-VALUE
C    (11) THIRD EXIT CHANNEL L-VALUE
C    (12) FOURTH EXIT CHANNEL L-VALUE
C
C     SOME OF THESE PARAMETERS WILL BE CONVERTED TO THE FORM REQUIRED
C     FOR CALCULATIONS.
C
C     (2) NEUTRON WIDTH DIVIDED BY PENETRATION FACTOR
C   (5-8) INELASTIC WIDTH DIVIDED BY PENETRATION FACTOR
C
  160 CALL LISPAR12(LOW,NLSJ12)
      IF(IMEDIT.NE.0) WRITE(OUTP,300)
      DO 200 JR=LOW,LHI
C-----LIST RESONANCE PARAMETERS.
      IF(IMEDIT.EQ.0) GO TO 180
      LOUT1=RESTAB( 9,JR)
      LOUT2=RESTAB(10,JR)
      LOUT3=RESTAB(11,JR)
      LOUT4=RESTAB(12,JR)
      DO 170 I=1,6
  170 CALL OUT9(RESTAB(I,JR  ),FIELD(1,I))
      CALL OUT9(RESTAB(7,JR  ),FIELD(1,7))
      CALL OUT9(RESTAB(8,JR  ),FIELD(1,8))
      WRITE(OUTP,310) ((FIELD(M,I),M=1,11),I=1,8),LOUT1,LOUT2,
     1 LOUT3,LOUT4
C-----04/21/07 - IF NECESSARY DEFINE RHOX2
  180 ERABS=DABS(RESTAB(1,JR))
      IF(NRO.EQ.1.AND.NAPS.EQ.1) CALL SETRHO1(ERABS,NSECT)
      RHO2X=ERABS*RHOX2(NSECT)
C-----DEFINE NEUTRON WIDTH DIVIDED BY PENETRATION FACTOR.
      CALL FACTS2(LVALUE(NSECT),RHO2X,SHIFT2(JR),PENFAC)
C-----02/14/04 - PROTECT AGAINST ZERO ENERGY RESONANCES
      IF(PENFAC.NE.0.0D+00) THEN
      RESTAB(2,JR)=RESTAB(2,JR)/PENFAC
      ELSE
      RESTAB(2,JR)=0.0D+00
      ENDIF
C-----DEFINE INELASTIC WIDTHS DIVIDED BY PENETRATION FACTOR.
      IF(NHRF.LE.0) GO TO 200
c----- 5 -  9 = WIDTHS
C----- 9 - 12 = L VALUES
C-----INITIALIZE INDICES TO WIDTH AND EXIT CHANNEL L-VALUE.
      JJ1=4
      JJ2=8
      DO 190 IHRF=1,NHRF
      JJ1=JJ1+1
      JJ2=JJ2+1
      IF(MTHRF(IHRF).NE.51) GO TO 190
      RHOZ2=DABS(RESTAB(1,JR)+EXCITE(IHRF,NSECT))*RHOX2(NSECT)
      LCOM=RESTAB(JJ2,JR)
c-----10/10/10 - switched from FACTS2 to FACTS3 - shift not used
      CALL FACTS3(LCOM,RHOZ2,PENFAC)
C-----02/14/04 - PROTECT AGAINST ZERO ENERGY RESONANCES
      IF(PENFAC.NE.0.0D+00) THEN
      RESTAB(JJ1,JR)=RESTAB(JJ1,JR)/PENFAC
      ELSE
      RESTAB(JJ1,JR)=0.0D+00
      ENDIF
  190 CONTINUE
C-----END OF RESONANCE LOOP.
  200 CONTINUE
C-----END OF S AND J LOOPS.
  210 CONTINUE
C
C     PRINT WARNING IF ANY J SEQUENCES ARE MISSING AND ADD MISSING J
C     SEQUENCE WITH NO RESONANCES IN ORDER TO ALLOW POTENTIAL
C     SCATTERING TO BE CORRECTLY CALCULATED.
C
      LSECT = NSECT
      CALL MISSINGJ(LNOW,LSECT)
C-----END OF L LOOP.
  220 CONTINUE
      RETURN
  230 FORMAT(1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Number of L Values-------------------',I11)
  240 FORMAT(1X,78('=')/
     1       ' Number of Capture Reactions----------',I11,' (0 TO 1)'/
     2       ' Number of Fission Reactions----------',I11,' (0 TO 1)'/
     3       ' Number of Inelastic Reactions--------',I11,' (0 TO 4)'/
     3       ' Number of Charged Particle Reactions-',I11,' (0 TO 4)')
  250 FORMAT(1X,78('='))
  260 FORMAT(' MT and Q-Value (eV)------------------',I11,11A1,1X,3A4)
  270 FORMAT(1X,78('=')/
     1       ' Atomic Weight Ratio of Isotope-------',11A1/
     2       ' L-Value------------------------------',I11/
     3       ' Number of S-Values (Channel Spin)----',I11)
  280 FORMAT(1X,78('=')/
     1       ' S-Value------------------------------',11A1/
     2       ' Number of J-Values-------------------',I11)
  290 FORMAT(1X,78('=')/
     1       ' J-Value------------------------------',11A1/
     2       ' Channel Radius (10E-12 cm)-----------',11A1/
     3       ' Background Tabulated-----------------',I11,2A4/
     4       ' Optical Model Phase SshiftTabulated--',I11,2A4/
     5       ' Number of Resonances-----------------',I11)
  300 FORMAT(1X,112('=')/' Hybrid R-Function Resonance Parameters'/
     1 1X,112('=')/
     2 '                                              Competitive',
     3 ' Widths                          Exit Channel'/
     2 '      Energy    Neutron    Capture    Fission Reaction-1',
     4 ' Reaction-2 Reaction-3 Reaction-4 L-Values'/
     6 '        (eV)       (eV)       (eV)       (eV)       (eV)',
     7 '       (eV)       (eV)       (eV) L-1 L-2 L-3 L-4'/
     8 1X,112('='))
  310 FORMAT(1X,88A1,4I4)
      END
      SUBROUTINE RDGRM
C=======================================================================
C
C     GENERAL R-MATRIX TREATMENT IS NOT YET INCLUDED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      WRITE(OUTP,10)
      WRITE(*   ,10)
      CALL ENDERROR
   10 FORMAT(///' ERROR -',
     1 ' General R-Matrix Parameters are NOT Yet Treated by this'/
     2 ' Program....Execution Terminated.'/1X,78('=')/
     3 ' In Order to Demonstrate that this Formalism is Generally'/
     4 ' Useful the Author of this Code Asks that Evaluators Send'/
     5 ' Him Evaluations Using this Formalism and Their Calculated'/
     6 ' Cross Sections for Comparison.'/1X,78('=')/
     7 ' Until Parameters and Calculated Cross Sections Using this'/
     8 ' Formalism are Available it is NOT Possible to Verify the'/
     9 ' Accuracy of Any Results Calculated by this Code.'/1X,78('='))
      END
      SUBROUTINE RDUR
C=======================================================================
C
C     READ UNRESOLVED RESONANCE PARAMETERS. THIS ROUTINE HANDLES ALL
C     POSSIBLE REPRESENTATIONS OF UNRESOLVED PARAMETERS,
C     (1) NO FISSION WIDTHS, ALL PARAMETERS ENERGY INDEPENDENT
C     (2) FISSION WIDTHS GIVEN, ONLY FISSION WIDTHS ENERGY DEPENDENT
C     (3) FISSION WIDTHS GIVEN, ALL PARAMETERS ENERGY DEPENDENT
C
C     FOR SIMPLICITY IN LATER CALCULATIONS ALL INPUT REPRESENTATIONS
C     ARE INTERNALLY CONVERTED TO THE ALL PARAMETERS ENERGY DEPENDENT
C     FORM TREATING EACH (L, J) STATE AS A SEPARATE SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      CHARACTER*4 TYPINT
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/NAPRHO/NRO,NAPS
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RANGER/LOW,LHI
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/WHATZA/IZANOW,MATNOW,TEMP3,IVERSE,INT45
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
      COMMON/LRUNOW/LRUIN
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
      COMMON/FIELDC/FIELD(11,12)
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION DUMSET(6),DUMBO(6),TYPINT(4,6)
      EQUIVALENCE (DUMSET(1),DX),(DUMSET(2),AJ),(DUMSET(3),AMUN),
     1 (DUMSET(4),GNO),(DUMSET(5),GG)
      DATA TYPINT/
     1 '(His','togr','am) ','    ',
     2 '(Lin',' X-L','in Y',')   ',
     3 '(Log',' X-L','in Y',')   ',
     4 '(Lin',' X-L','og Y',')   ',
     5 '(Log',' X-L','og Y',')   ',
     6 '(ERR','OR) ','    ','    '/
      DATA DUMSET/6*0.0D+00/
      DATA DUMBO/6*0.0D+00/
      DATA C3/1.23D-01/
      DATA C5/8.0D-02/
C-----UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITTEE RECOMMENDATION
c***** Updated 9/01/04 - ENDF-102, appendix H
c     DATA C4/1.008664904D+00/
      DATA C4/1.00866491578D+00/
c***** Updated 9/01/04
      DATA C6/2.196807122623D-03/
C-----INITIALIZE NEGATIVE WIDTH COUNT
      NEGSUM = 0
C-----PRINT TITLE FOR OUTPUT LISTING.
      IF(LRF.EQ.1.AND.LFW.EQ.0) WRITE(OUTP,260)
      IF(LRF.EQ.1.AND.LFW.EQ.1) WRITE(OUTP,270)
      IF(LRF.EQ.2) WRITE(OUTP,280)
      IF(LRF.EQ.1.AND.LFW.EQ.0) WRITE(*   ,260)
      IF(LRF.EQ.1.AND.LFW.EQ.1) WRITE(*   ,270)
      IF(LRF.EQ.2) WRITE(*   ,280)
C
C     READ BEGINNING OF EACH TYPE OF REPRESENTATION AND DEFINE ALL
C     COMMON TERMS.
C
C-----DEFINE TARGET SPIN, SPIN UP SCATTERING RADIUS AND,
C-----(LRF=1, LFW=1) - NUMBERS OF ENERGIES AT WHICH FISSION WIDTHS ARE
C-----                 GIVEN AND NUMBER OF L VALUES, OR
C-----(OTHERWISE)    - NUMBER OF L VALUES.
C
C     FOR ENDF/B-VI FORMATS CHECK L1 FIELD (LSSF) TO SEE IF THE
C     UNRESOLVED INFINITELY DILUTE CROSS SECTION HAS ALREADY BEEN
C     INCLUDED IN THE FILE 3 CROSS SECTIONS. IF IT HAS CHANGE LRUIN
C     TO EFFECTIVELY IGNORE SECTION WHILE CALCULATING RESONANCE
C     CONTRIBUTION TO THE FILE 3 CROSS SECTIONS.
C
      CALL CARDI(SPI,AP,LSSF,L2,NE,NLS)
      CALL CARDO(SPI,AP,LSSF,L2,NE,NLS)
      IF(LFW.NE.1.OR.LRF.NE.1) NLS=NE
C-----LIST SPIN AND SPIN UP SCATTERING.
      CALL OUT9(SPI,FIELD(1,1))
      CALL OUT9(AP ,FIELD(1,2))
      WRITE(OUTP,320) ((FIELD(M,I),M=1,11),I=1,2)
      WRITE(*   ,320) ((FIELD(M,I),M=1,11),I=1,2)
C-----FOR ENDF/B-VI PRINT INTERPRETATION OF LSSF.
      IF(IVERSE.NE.6) GO TO 10
      if(LSSF.eq.0) then
      WRITE(OUTP,300) LSSF
      WRITE(*   ,300) LSSF
      else
      WRITE(OUTP,310) LSSF
      WRITE(*   ,310) LSSF
      endif
C-----LIST NUMBER OF L VALUES.
   10 WRITE(OUTP,330) NLS
      WRITE(*   ,330) NLS
C
C     FOR ENDF/B-VI IF LRUIN INDICATES THAT CROSS SECTIONS NOT YET
C     ADDED TO FILE 3, BUT LSSF INDICATES THAT THEY HAVE PRINT
C     WARNING AND CHANGE LRUIN TO IGNORE THIS SECTION.
C
      IF(IVERSE.NE.6.OR.LRUIN.EQ.5.OR.LSSF.EQ.0) GO TO 20
      WRITE(OUTP,290) LSSF
      WRITE(*   ,290) LSSF
      LRUIN=5
C
C     IF ONLY FISSION WIDTH ARE ENERGY DEPENDENT READ ENERGIES AT WHICH
C     FISSION WIDTHS ARE GIVEN.
C
   20 IF(LRF.NE.1.OR.LFW.NE.1) GO TO 40
C-----INCREMENT SECTION COUNT (ALLOWING ROOM TO INSERT THE DEGREES OF
C-----FREEDOM BEFORE THE ENERGIES).
      NRS=NE+1
      CALL LIMIT1(1)
      LOWP1=LOW+1
C-----READ ENERGIES AT WHICH FISSION WIDTHS ARE GIVEN AND MOVE THEM
C-----INTO THE 1-ST POSITION FOR EACH RESONANCE.
      JRC=(LHI-LOWP1)+1
      CALL LISTIO(ENRES(LOWP1),JRC)
      DO 30 JR1=LOWP1,LHI
      CALL NOODLE(ENRES(JR1),ZERO,ELX,EHX)
   30 RESTAB(1,JR1)=ENRES(JR1)
C-----SAVE INDEX TO FIRST ENERGY TO ALLOW COPY OF ENERGIES FOR EACH
C-----SECTION.
      LOWP1X=LOWP1
C
C     L STATE LOOP.
C
   40 DO 240 IL=1,NLS
C-----DEFINE ATOMIC WEIGHT RATIO, L VALUE AND NUMBER OF J STATES.
      CALL CARDIO(AWRI,C2,LNOW,L2,NJS,NJSX)
      IF(LRF.EQ.1.AND.LFW.EQ.0) NJS=NJSX
      CALL OUT9(AWRI,FIELD(1,1))
      WRITE(OUTP,340) (FIELD(M,1),M=1,11),LNOW,NJS
      WRITE(*   ,340) (FIELD(M,1),M=1,11),LNOW,NJS
C-----DEFINE COMMON PARAMETERS FOR ALL SECTIONS WITH SAME L VALUE.
C-----EITHER CALCULATE CHANNEL RADIUS OR DEFINE TO BE EQUAL TO THE
C-----SCATTERING RADIUS.
      IF(NAPS.NE.0) GO TO 50
      APX=C3*((C4*AWRI)**(ONE/THREE))+C5
      GO TO 60
   50 APX=AP
   60 AK1=C6*AWRI/(AWRI+ONE)
      BETAD=PI*ABN/(AK1*AK1)
      RHOX2D=(APX*AK1)**2
      RHOP1D=AP*AK1
C
C     J STATE LOOP.
C
      DO 240 IJ=1,NJS
C-----READ NEXT LINE FOR (LRF=1, LFW=1 OR LRF=2, ANY LFW...NOTHING TO
C-----READ FOR LRF=1, LFW=0).
      IF(LRF.NE.2) GO TO 70
C-----DEFINE J VALUE, INTERPOLATION LAW AND NUMBER OF ENERGIES.
      CALL CARDIO(AJ,C2,INTX,L2,NET6P6,NE)
C-----DEFINE THE NUMBER OF RESONANCE LOCATIONS WHICH WILL BE USED
C-----ALLOWING FOR THE PRECEDING DEGREES OF FREEDOM.
      NRS=NE+1
      GO TO 90
   70 IF(LFW.EQ.0) GO TO 80
C-----DEFINE THE NUMBER OF DEGREES OF FREEDOM FOR FISSION (ASSUME THE
C-----SAME NUMBER OF ENERGIES, I.E. NEED NOT RE-DEFINE NE OR NRS AT
C-----THIS POINT).
      CALL CARDIO(C1,C2,L1,MUF,NEP6,N2)
      IF(IL.EQ.1.AND.IJ.EQ.1) GO TO 100
      GO TO 90
C-----INDICATE 3 LOCATIONS PER (L, J) STATE (FIRST LOCATION CONTAINS
C-----PARAMETERS AND DEGREES OF FREEDOM, THE SECOND AND THIRD CONTAIN
C-----ENERGY POINTS AT THE LOWER AND UPPER ENERGY LIMITS OF UNRESOLVED
C-----REGION.
   80 NRS=3
C-----INCREMENT SECTION COUNT AND INDICES TO RESONANCE PARAMETER TABLE
   90 CALL LIMIT1(1)
  100 LOWP1=LOW+1
C-----DEFINE ALL PARAMETERS FOR A SECTION.
      BETA(NSECT)=BETAD
      RHOX2(NSECT)=RHOX2D
      RHOP1(NSECT)=RHOP1D
      EL(NSECT)=ELX
      EH(NSECT)=EHX
      NLOW(NSECT)=LOW
      NHIGH(NSECT)=LHI
      LVALUE(NSECT)=LNOW
c----- 02/08/08 - ADDED NAPS DEFINITION
      NAPTAB(NSECT)=NAPS
      MODE(NSECT)=11
C
C     06/09/05 - THIS WAS ALREADY DONE FOR BREIT-WIGNER - NOW ALL.
C
C     IF ENERGY DEPENDENT SCATTERING RADIUS DEFINE INDEX TO TABULATED
C     DATA (OTHERWISE INDEX HAS ALREADY BEEN INITIALIZED TO 0).
C     CALCULATE ENERGY DEPENDENT RHOP FROM SCATTERING RADII.
C
      IF(NRO.GT.0) THEN
      NRHO(NSECT)=NUMRHO
      INX3=INXRHO(3,NUMRHO)
      INX4=INXRHO(4,NUMRHO)
      DO I=INX3,INX4
      RHOTAB(I)=APTAB(I)*AK1
      ENDDO
      ENDIF
C
C     SELECT PARAMETER REPRESENTATION.
C
C-----ARE ALL WIDTHS ENERGY DEPENDENT.
      IF(LRF.EQ.2) GO TO 160
C-----ARE FISSION WIDTHS GIVEN.
      IF(LFW.NE.0) GO TO 120
C
C     ALL PARAMETERS ARE ENERGY INDEPENDENT. CONVERT TO ENERGY DEPENDENT
C     FORM USING INTERPOLATION LAW 1 (HISTOGRAM).
C
C-----READ LEVEL SPACING, J, NEUTRON DEGREES OF FREEDOM, NEUTRON AND
C-----CAPTURE WIDTHS (SEE, ABOVE EQUIVALENCE TO MEMBERS OF DUMSET).
      CALL LISTIO(DUMSET,6)
C-----DEFINE INTERPOLATION LAW TO BE HISTOGRAM.
      INTX=1
C-----DEFINE NUMBER OF DEGREES OF FREEDOM FOR COMPETITION, CAPTURE
C-----AND FISSION (10) AND ELASTIC (AS INPUT).
      RESTAB(3,LOW)=10.0
      RESTAB(4,LOW)=AMUN
      RESTAB(5,LOW)=10.0
      RESTAB(6,LOW)=10.0
C-----COPY PARAMETERS AS CONSTANT BETWEEN LOWER AND UPPER ENERGY LIMITS
C-----OF THE UNRESOLVED REGION.
      ENRES(LOWP1)=ELX
      ENRES(LHI)=EHX
      RESTAB(1,LOWP1)=ELX
      RESTAB(1,LHI)=EHX
      DO 110 I=LOWP1,LHI
      RESTAB(2,I)=DX
      RESTAB(3,I)=0.0
      RESTAB(4,I)=GNO
      RESTAB(5,I)=GG
  110 RESTAB(6,I)=0.0
      GO TO 170
C
C     ENERGY DEPENDENT FISSION WIDTHS, ALL OTHER PARAMETERS ARE ENERGY
C     INDEPENDENT. TREAT EACH (L, J) STATE AS A SECTION AND CONVERT
C     DATA TO ALL PARAMETERS ENERGY DEPENDENT FORM USING IMPLIED ENDF/B
C     VERSION DEPENDENT INTERPOLATION LAW.
C
C-----READ LEVEL SPACING, J, NEUTRON DEGREES OF FREEDOM, NEUTRON AND
C-----CAPTURE WIDTHS (SEE, ABOVE EQUIVALENCE TO MEMBERS OF DUMSET).
  120 CALL LISTIO(DUMSET,6)
C-----READ ENERGY DEPENDENT FISSION WIDTHS UP TO 6 AT A TIME AND THEN
C-----MOVE THEM INTO THE 6-TH POSITION FOR EACH RESONANCE.
      JR=LOWP1
      DO 140 JR1=LOWP1,LHI,6
      JR2=JR1+5
      IF(JR2.GT.LHI) JR2=LHI
      JR3=(JR2-JR1)+1
      CALL LISTIO(DUMBO,JR3)
      DO 130 K=1,JR3
      RESTAB(6,JR)=DUMBO(K)
  130 JR=JR+1
  140 CONTINUE
C-----DEFINE INTERPOLATION LAW BASED ON ENDF/B FORMAT VERSION (E.G.,
C-----ENDF/B-IV = LOG-LOG, ENDF/B-V = LINEAR-LINEAR).
      INTX=INT45
C-----DEFINE NUMBER OF DEGREES OF FREEDOM FOR COMPETITION AND CAPTURE
C-----(10) AND ELASTIC AND FISSION (AS INPUT).
      RESTAB(3,LOW)=10.0
      RESTAB(4,LOW)=AMUN
      RESTAB(5,LOW)=10.0
      RESTAB(6,LOW)=MUF
C-----COPY ENERGIES AT WHICH FISSION WIDTHS ARE GIVEN FROM FIRST SECTION
C-----COPY SPACING, COMPETITIVE, NEUTRON AND CAPTURE WIDTHS AS CONSTANTS
      JRX=LOWP1X
      DO 150 JR=LOWP1,LHI
      ENRES(JR)=ENRES(JRX)
      RESTAB(1,JR)=RESTAB(1,JRX)
      RESTAB(2,JR)=DX
      RESTAB(3,JR)=0.0
      RESTAB(4,JR)=GNO
      RESTAB(5,JR)=GG
  150 JRX=JRX+1
      GO TO 170
C
C     ALL WIDTHS ARE ENERGY DEPENDENT.
C
C-----READ NUMBER OF DEGREES OF FREEDOM FOLLOWED BY ALL RESONANCE
C-----PARAMETERS.
  160 CALL LISTIO(RESTAB(1,LOW),6)
      CALL LISPAR6(LOW+1,6*(NRS-1))
C
C     LIST DATA CONVERTED TO ALL PARAMETERS ENERGY DEPENDENT FORM.
C
C-----SAVE INTERPOLATION LAW AND DEFINE STATISTICAL WEIGHT (THESE 2
C-----QUANTITIES ARE SAVED IN THE 2 UNUSED LOCATIONS PRECEDING THE
C-----NUMBER OF DEGREES OF FREEDOM FOR EACH REACTION).
  170 RESTAB(1,LOW)=INTX
      CALL GJWAIT(SPI,AJ,GJ)
      RESTAB(2,LOW)=GJ
C-----LIST INTERPOLATION LAW, J AND DEGREES OF FREEDOM.
      IF(IMEDIT.EQ.0) GO TO 180
      LSTINT=6
      IF(INTX.GE.1.AND.INTX.LE.5) LSTINT=INTX
      WRITE(OUTP,350) INTX,(TYPINT(I,LSTINT),I=1,4),AJ,
     1 (RESTAB(K,LOW),K=3,6)
      IF(IMEDIT.NE.2) GO TO 180
C-----IN EDIT MODE IF INTERPOLATION LAW IS INCORRECT PRINT WARNING.
      IF(LSTINT.EQ.6) WRITE(OUTP,370) INTX
      GO TO 200
C-----IN CALCULATION MODE IF INTERPOLATION LAW IS INCORRECT PRINT
C-----WARNING AND TERMINATE.
  180 IF(INTX.GE.1.AND.INTX.LE.5) GO TO 200
      WRITE(OUTP,190) INTX
      WRITE(*   ,190) INTX
  190 FORMAT(///' ERROR - Unresolved Interpolation Law=',I5,
     1 ' (MUST be 1 to 5).'/
     2          '         Cannot Interpolate Unresolved Data.'/
     3          '         Execution Terminated'///)
      CALL ENDERROR
C-----LIST PARAMETERS.
  200 DO 230 JR=LOWP1,LHI
C-----IF FISSION WIDTH IS NOT ZERO TURN ON FISSILE FLAG.
      IF(DABS(RESTAB(6,JR)).NE.0.0) LFWX=1
      IF(DABS(RESTAB(6,JR)).NE.0.0) LFWL(NSECT) = LFWL(NSECT) + 1
C-----CHECK FOR NEGATIVE WIDTHS
      NEGRES = 0
      DO I=3,6
      IF(RESTAB(I,JR).LT.0.0D+0) NEGRES = NEGRES + 1
      ENDDO
      IF(IMEDIT.EQ.0.AND.NEGRES.EQ.0) GO TO 230
      DO 210 I=1,6
  210 CALL OUT9(RESTAB(I,JR),FIELD(1,I))
      GAMT=RESTAB(3,JR)+RESTAB(4,JR)+RESTAB(5,JR)+RESTAB(6,JR)
      CALL OUT9(GAMT        ,FIELD(1,7))
      WRITE(OUTP,360) ((FIELD(M,I),M=1,11),I=1,7)
C-----NEGATIVE WIDTHS?
      IF(NEGRES.GT.0) THEN
      NEGSUM = NEGSUM + NEGRES
      WRITE(   *,360) ((FIELD(M,I),M=1,11),I=1,7)
      WRITE(   *,220) NEGRES
      WRITE(OUTP,220) NEGRES
  220 FORMAT(' ERROR - ',I1,' Negative Widths')
      ENDIF
  230 CONTINUE
C-----END OF L AND J STATE LOOPS.
  240 CONTINUE
C
C     STOP IF NEGATIVE WIDTHS
C
      IF(NEGSUM.GT.0) THEN
      WRITE(   *,250) NEGSUM
      WRITE(OUTP,250) NEGSUM
  250 FORMAT(///' ERROR -',
     1 ' Execution Terminated because of',i8,' Negative Widths.'///)
      CALL ENDERROR
      ENDIF
      RETURN
  260 FORMAT(1X,78('=')/' Energy Independent Unresolved Parameters')
  270 FORMAT(1X,78('=')/' Unresolved Fission Widths Energy Dependent')
  280 FORMAT(1X,78('=')/' All Unresolved Parameters Energy Dependent')
  290 FORMAT(1X,78('=')/ 1X,7('WARNING...'),'WARNING'/
     1 ' Note, LSSF=',I2,' Indicates that the Resonance Contribution'/
     2 ' has Already been Added to the Cross Sections. This Section'/
     3 ' will be Read, but Ignored in ALL Resonance Calculations.')
  300 FORMAT(' LSSF---------------------------------',I11,
     1 ' (Add to File 3)')
  310 FORMAT(' LSSF---------------------------------',I11,
     1 ' (Included in File 3)')
  320 FORMAT(1X,78('=')/
     1       ' Nuclear Spin of Target---------------',11A1/
     2       ' Effective Scattering Radius (A+)-----',11A1)
  330 FORMAT(' Number of L Values-------------------',I11)
  340 FORMAT(1X,78('=')/
     1       ' Atomic Weight Ratio of Isotope-------',11A1/
     2       ' L Value------------------------------',I11/
     3       ' Number of J Values-------------------',I11)
  350 FORMAT(1X,78('=')/32X,' Degrees of Freedom'/1X,78('=')/
     1 '   Interpolation  J Value  Competition    Neutron   Capture ',
     2 '    Fission'/7X,' Law'/1X,78('=')/I2,1X,3A4,A1,
     3 F9.2,6X,F6.2,1X,3(5X,F6.2)/
     3 1X,78('=')/' Resonance Parameters'/1X,78('=')/'      Energy',
     3 '     Level Competition    Neutron    Capture    Fission',
     4 '      Total'/
     5 12X,'    Spacing',5(6X,'Width')/1X,7(7X,'(eV)')/1X,78('='))
  360 FORMAT(1X,77A1)
  370 FORMAT(1X,78('=')/1X,7('WARNING...'),'WARNING'/
     1 ' Unresolved Interpolation Law=',I5,
     2 ' (MUST be 1 to 5). Cannot Interpolate Unresolved Data.'/
     3 ' Correct Evaluated Data Before Using it for Calculations.'/
     4 ' If Data is NOT Corrected this Program will Abort During',
     5 ' Calculations.'/' Execution Terminated'/1X,78('='))
      END
      SUBROUTINE LIMIT1(IPATH)
C=======================================================================
C
C     CORE ALLOCATION.
C
C     INCREMENT SECTION COUNT AND LOWER AND/OR UPPER INDICES TO
C     RESONANCE PARAMETER TABLE.
C
C     IF AVAILABLE CORE ALLOCATION IS EXCEEDED TERMINATE UNLESS PROGRAM
C     IS IN THE EDIT MODE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,ADDSEC,ADDRES
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RANGER/LOW,LHI
      COMMON/MAXIE/NEDSEC,NEDRES,NEDNOD
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DATA ADDSEC/0/
      DATA ADDRES/0/
C-----SELECT PATH.
      IF(IPATH.EQ.2) GO TO 30
C-----INCREMENT SECTION COUNT AND INSURE THAT CORE ALLOCATION WILL NOT
C-----BE EXCEEDED.
C-----11/24/2012 - ELIMINATED NHIGH(0) INDEXING ERROR.
      IF(NSECT.GT.0) THEN
      LHI = NHIGH(NSECT)
      ELSE
      LHI = 0
      ENDIF
C-----11/24/2012 - ELIMINATED NHIGH(0) INDEXING ERROR.
      NSECT=NSECT+1
      IF((NSECT+ADDSEC).GT.NEDSEC) NEDSEC=NSECT+ADDSEC
      IF(NSECT.LE.MAXSEC) GO TO 30
C-----AVAILABLE CORE EXCEEDED. TERMINATE UNLESS IN EDIT MODE.
      IF(IMEDIT.NE.2) GO TO 10
      ADDSEC=MAXSEC
      NSECT=1
      GO TO 30
   10 WRITE(OUTP,20) MAXSEC
      WRITE(*   ,20) MAXSEC
   20 FORMAT(///' ERROR - More than',I6,' Sections'/
     1          '         Re-Run Program in Edit Mode to Determine'/
     1          '         Memory Requirements.'/
     2          '         Execution Terminated'///)
      CALL ENDERROR
C
C     INCREMENT LOWER AND UPPER INDICES TO RESONANCE PARAMETER TABLE.
C     PARAMETER TABLE.
C
C-----DEFINE INDICES FOR NEXT SET OF RESONANCES TO READ AND INSURE THAT
C-----CORE ALLOCATION WILL NOT BE EXCEEDED.
   30 LOW=LHI+1
      LHI=LHI+NRS
      IF((LHI+ADDRES).GT.NEDRES) NEDRES=LHI+ADDRES
      IF(LHI.LE.MAXRES) GO TO 60
C-----AVAILABLE CORE EXCEEDED. TERMINATE UNLESS IN EDIT MODE.
      IF(IMEDIT.NE.2) GO TO 40
      ADDRES=ADDRES+(LOW-1)
      LOW=1
      LHI=NRS
      GO TO 60
   40 WRITE(OUTP,50) MAXRES
      WRITE(*   ,50) MAXRES
   50 FORMAT(///' ERROR - More than',I6,' Resonances'/
     1          '         Re-Run Program in Edit Mode to Determine'/
     2          '         Requirements.'/
     3          '         Execution Terminated'///)
      CALL ENDERROR
   60 RETURN
      END
      SUBROUTINE SIGMA(E,SIGNOW)
C=======================================================================
C
C     CALCULATE CONTRIBUTION FROM ALL SECTIONS.
C
C     NOTE - THE ALGORITHMS FOR ALL POSSIBLE TYPES OF RESONANCE
C            PARAMETERS HAVE BEEN WRITTEN IN A FORM TO REMOVE A
C            COMMON FACTOR OF,
C
C            PI*ABUNDANCE*(LAMBDA**2)
C
C            AFTER CALCULATING THE CONTRIBUTION FROM ANY TYPE OF
C            RESONANCE PARAMETER THIS FACTOR IS USED TO MULTIPLY
C            THE RESULT TO DEFINE THE FINAL CROSS SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      common/rmlfinal/rmlsigma(11)
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
c-----2017/4/13 - Deleted Lrf
      COMMON/MRMLWCOM/Npp,Ngroup,Nres,Nro7,Naps7,Kg,Minr,Maxr
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION SIGNOW(*)
C
C-----IF ANY SECTIONS HAVE ENERGY DEPENDENT SCATTERING RADIUS DEFINE
C-----RADIUS AT ENERGY E.
      IF(NUMRHO.GT.0) CALL SETRHO(E)
C-----INITIALIZE ALL CROSS SECTIONS (EXCEPT TOTAL)
      DO J=2,NREACT
      SIGNOW(J)=0.0
      ENDDO
      SIGXSUM = 0.0
C
C     ADD CONTRIBUTION FROM EACH SECTION.
C
      DO 140 ISECT=1,NSECT
C-----SELECT PARAMETER TYPE.
      MODE1=MODE(ISECT)
C-----SKIP SECTION IF IT DOES NOT CONTRIBUTE TO THIS ENERGY RANGE.
      IF(MODE1.LE.0) GO TO 140
C-----DEFINE INDICES FOR CURRENT SECTION OF RESONANCES.
      LOW=NLOW(ISECT)
      LHI=NHIGH(ISECT)
      Ngr1 = NGRTAB(1,ISECT)
      Ngr2 = NGRTAB(2,ISECT)
      Npp  = NPPTAB(ISECT)
C-----INITIALIZE CONTRIBUTION OF SECTION TO SUM.
      SIGNTI=0.0
      SIGNNI=0.0
      SIGNGI=0.0
      SIGNFI=0.0
C-----RESOLVED    = 1 TO 10
C-----UNRESOLVED = 11 TO 12 (12 NO LONGER USED)
C              1   2   3   4   5   6   7   8   9  10
      GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 90,
     1       110,120), MODE1
C
C     UNRESOLVED
C
C-----LRF = 1) SINGLE LEVEL BREIT-WIGNER.
   10 CALL SIGBW1(E)
      GO TO 130
C-----LRF = 2) MULTI-LEVEL BREIT-WIGNER.
   20 CALL SIGBWM(E)
      GO TO 130
C-----LRF = 3) REICH-MOORE.
   30 CALL SIGRM1(E)
      GO TO 130
C-----LRF = 4) ADLER-ADLER (ENDF/B-V AND LATER).
   40 CALL SIGAA5(E)
      GO TO 130
C-----LRF = 5) GENERAL R-MATRIX.
   50 CALL SIGGRM(E)
      GO TO 130
C-----LRF = 6) HYBRID R-FUNCTION.
   60 CALL SIGHRF(E)
      GO TO 130
C-----LRF = 7) NEW (2003) REICH-MOORE WITH COMPETITION
   70 CALL SIGRML(E)
      do j=1,NREACT
      SIGNOW(j) = rmlsigma(j)
      enddo
C-----SKIP NORMALIZATION
      GO TO 140
C-----LRF = 4A) ADLER-ADLER (ENDF/B-IV AND EARLIER).
   80 CALL SIGAA4(E)
      GO TO 130
C-----ILLEGAL
   90 WRITE(3,100)
      WRITE(*,100)
  100 FORMAT(///' ERROR - Illegal Resonance Parameter Mode.'/
     1          '         Execution Terminated.'///)
      CALL ENDERROR
C
C     UNRESOLVED
C
C-----UNRESOLVED (INTERPOLATE PARAMETERS).
  110 CALL SIGURP(E)
      GO TO 130
C-----UNRESOLVED (INTERPOLATE CROSS SECTIONS) - NO LONGER USED
  120 CALL SIGURS(E)
C
C     ALL
C
C-----MULTIPLY CONTRIBUTION OF SECTION BY PI*ABUNDANCE*(LAMBDA**2)
C-----(WHICH IS BETAE) AND ADD TO SUM.
  130 BETAE=BETA(ISECT)/E
      SIGNOW(2)=SIGNOW(2)+SIGNNI*BETAE
      SIGNOW(3)=SIGNOW(3)+SIGNGI*BETAE
      SIGNOW(4)=SIGNOW(4)+SIGNFI*BETAE
      SIGXSUM  =SIGXSUM  +SIGNXI*BETAE  ! OPTIONAL COMPETIRION
C-----END OF SECTION LOOP.
  140 CONTINUE
C-----DEFINE TOTAL TO EQUAL SUM OF PARTS.
      TOTSUM = 0.0
      DO J=2,NREACT
      TOTSUM = TOTSUM + SIGNOW(J)
      ENDDO
      SIGNOW(1) = TOTSUM
C***** DEBUG - ACTIVATE FOR UNRESOLVED COMPETITION LISTING
C     IF(SIGXSUM.GT.0.0) WRITE(22,2200) E,SIGNOW(2),SIGNOW(3),SIGXSUM
C2200 FORMAT(1P5D12.4)
C***** DEBUG - ACTIVATE FOR UNRESOLVED COMPETITION LISTING
      RETURN
      END
      SUBROUTINE SETRHO(E)
C=======================================================================
C
C     FOR ALL SECTIONS WHICH HAVE AN ENERGY DEPENDENT SCATTERING RADIUS
C     INTERPOLATE TABULATED VALUES TO DEFINE SCATTERINMG RADIUS AT E.
C     STORE INTERPOLATED VALUES IN THE ARRAY RHOP1 (THE REMAINDER OF
C     THE CALCULATIONS CAN THEN PROCEED IGNORING ENERGY DEPENDENCE OF
C     THE SCATTERING RADIUS).
C
C=======================================================================
      INCLUDE 'implicit.h'
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C-----SELECT SECTIONS WITH ENERGY DEPENDENT SCATTERING RADIUS.
      DO 70 KSECT=1,NSECT
      IF(NRHO(KSECT).EQ.0) GO TO 70
C-----DEFINE INDICES TO INTERPOLATION LAW AND TABULATED DATA.
      LNX=NRHO(KSECT)
      INX1=INXRHO(1,LNX)
      INX2=INXRHO(2,LNX)
      INX3=INXRHO(3,LNX)
      INX4=INXRHO(4,LNX)
C-----FIND ENERGY OR ENERGY RANGE WHERE RADII ARE GIVEN.
      DO 10 L=INX3,INX4
      IF(E.lt.ERHOTB(L)) go to 20
      IF(E.eq.ERHOTB(L)) go to 50
   10 CONTINUE
C-----EXTEND RADIUS AS CONSTANT OUTSIDE TABULATED ENERGY RANGE.
      L=INX4
      GO TO 50
   20 IF(L.EQ.INX3) GO TO 50
C-----E IS BETWEEN ERHOTB(L-1) AND ERHOTB(L).
      L=L-1
C-----DEFINE INTERPOLATION LAW. DEFINE TRUE POINT INDEX.
      INX5=L-INX3+1
      DO 30 J=INX1,INX2
      IF(INX5.GT.NBTRHO(J)) GO TO 40
   30 CONTINUE
      J=INX2+1
C-----POINT IS IN INTERPOLATION REGION J-1.
   40 J=J-1
C-----INTERPOLATE RADIUS TO ENERGY E.
      CALL RHOINT(E,RHOP1(KSECT),ERHOTB(L),RHOTAB(L),INTRHO(J))
      GO TO 60
C-----DEFINE RADIUS SET AT E (E IS AN ENERGY AT WHICH RADIUS IS
C-----TABULATED OR E IS OUTSIDE TABULATED RANGE AND WILL BE EXTENDED
C-----AS CONSTANT).
   50 RHOP1(KSECT)=RHOTAB(L)
C-----IF REQUESTED SET CHANNEL RADIUS EQUAL TO SCATTERING RADIUS.
   60 IF(NAPTAB(KSECT).EQ.1) RHOX2(KSECT)=RHOP1(KSECT)**2
   70 CONTINUE
      RETURN
      END
      SUBROUTINE SETRHO1(ERABS,KSECT)
C=======================================================================
C
C     DEFINE SCATTERING RADIUSA, ETC., AT ONE RESONANCE ENERGY.
C
C     THIS IS A SHORTER VERSION OF SETRHO.
C
C=======================================================================
      INCLUDE 'implicit.h'
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C-----DEFINE INDICES TO INTERPOLATION LAW AND TABULATED DATA.
      LNX=NRHO(KSECT)
      INX1=INXRHO(1,LNX)
      INX2=INXRHO(2,LNX)
      INX3=INXRHO(3,LNX)
      INX4=INXRHO(4,LNX)
C-----FIND ENERGY OR ENERGY RANGE WHERE RADII ARE GIVEN.
      DO 10 L=INX3,INX4
      IF(ERABS.lt.ERHOTB(L)) go to 20
      IF(ERABS.eq.ERHOTB(L)) go to 50
   10 CONTINUE
C-----EXTEND RADIUS AS CONSTANT OUTSIDE TABULATED ENERGY RANGE.
      L=INX4
      GO TO 50
   20 IF(L.EQ.INX3) GO TO 50
C-----ERABS IS BETWEEN ERHOTB(L-1) AND ERHOTB(L).
      L=L-1
C-----DEFINE INTERPOLATION LAW. DEFINE TRUE POINT INDEX.
      INX5=L-INX3+1
      DO 30 J=INX1,INX2
      IF(INX5.GT.NBTRHO(J)) GO TO 40
   30 CONTINUE
      J=INX2+1
C-----POINT IS IN INTERPOLATION REGION J-1.
   40 J=J-1
C-----INTERPOLATE RADIUS TO ENERGY ERABS.
      CALL RHOINT(ERABS,RHOP1(KSECT),ERHOTB(L),RHOTAB(L),INTRHO(J))
      GO TO 60
C-----DEFINE RADIUS SET AT ERABS (ERABS IS AN ENERGY AT WHICH RADIUS IS
C-----TABULATED OR ERABS IS OUTSIDE TABULATED RANGE AND WILL BE EXTENDED
C-----AS CONSTANT).
   50 RHOP1(KSECT)=RHOTAB(L)
C-----IF REQUESTED SET CHANNEL RADIUS EQUAL TO SCATTERING RADIUS.
   60 RHOX2(KSECT)=RHOP1(KSECT)**2
      CONTINUE
      RETURN
      END
      SUBROUTINE RHOINT(E,RHOP1,ERHOTB,RHOTAB,INTX)
C=======================================================================
C
C     INTERPOLATE ENERGY DEPENDENT SCATTERING RADIUS.
C
C     THIS ROUTINE HAS BEEN RECODED IN ORDER TO AVOID ROUND-OFF
C     PROBLEMS ON SHORT WORD LENGTH COMPUTERS, E.G. IBM-360, 370, ETC.
C     THIS ROUTINE IS NOW SLIGHTLY LESS EFFICIENT THAN IN ITS PREVIOUS
C     FORM. HOWEVER ALL INTERPOLATION IS NOW DEFINED AS A WEIGHTED SUM
C     OF TERMS, AS OPPOSED TO THE PREVIOUS FORM WHICH USED DIFFERENCES.
C
C     AN ILLEGAL INTERPOLATION CODE OR A NON-POSITIVE ENERGY WHERE LOG
C     ENERGY INTERPOLATION IS REQUIRED INDICATES EITHER AN ERROR IN THE
C     DATA AS IT APPEARS IN THE ENDF/B FORMAT, OR AN ERROR IN THIS
C     PROGRAM. THEREFORE ERRORS OF THIS TYPE WILL CAUSE THE PROGRAM TO
C     PRINT A WARNING MESSAGE AND TERMINATE EXECUTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FIELDC/FIELD(11,12)
      DIMENSION ERHOTB(*),RHOTAB(*)
C-----CHECK INTERPOLATION CODE.
      IF(INTX.LT.1.OR.INTX.GT.5) GO TO 60
C-----DEFINE ENERGIES AT THE 2 ENDS OF THE INTERVAL.
      E1=ERHOTB(1)
      E2=ERHOTB(2)
C-----CHECK FOR ZERO LENGTH INTERVAL.
      DE12=E2-E1
      IF(DE12.EQ.0.0) GO TO 100
C-----SELECT INTERPOLATION LAW.
      GO TO (10,20,30,40,50),INTX
C
C     HISTOGRAM. CONSTANT EQUAL TO VALUE AT LOWER ENERGY LIMIT.
C
   10 RHOP1=RHOTAB(1)
      RETURN
C
C     LINEAR X AND LINEAR Y.
C
   20 WT2=(E-E1)/DE12
      WT1=ONE-WT2
      RHOP1=WT1*RHOTAB(1)+WT2*RHOTAB(2)
      RETURN
C
C     LOG X AND LINEAR Y.
C
C-----INSURE ALL X VALUES ARE POSITIVE FOR LOG.
   30 IF(E1.LE.0.0.OR.E2.LE.0.0.OR.E.LE.0.0) GO TO 80
      WT2=DLOG(E/E1)/DLOG(E2/E1)
      WT1=ONE-WT2
      RHOP1=WT1*RHOTAB(1)+WT2*RHOTAB(2)
      RETURN
C
C     LINEAR X AND LOG Y.
C
   40 WT2=(E-E1)/DE12
      WT1=ONE-WT2
      RHOP1=DEXP(WT1*DLOG(RHOTAB(1))+WT2*DLOG(RHOTAB(2)))
      RETURN
C
C     LOG X AND LOG Y.
C
   50 IF(E1.LE.0.0.OR.E2.LE.0.0.OR.E.LE.0.0) GO TO 80
      WT2=DLOG(E/E1)/DLOG(E2/E1)
      WT1=ONE-WT2
      RHOP1=DEXP(WT1*DLOG(RHOTAB(1))+WT2*DLOG(RHOTAB(2)))
      RETURN
C-----ILLEGAL INTERPOLATE CODE.
   60 WRITE(OUTP,70) INTX
      WRITE(*   ,70) INTX
   70 FORMAT(///' ERROR - Interpolating Scattering Radius.'/
     2          '         Illegal Interpolation Code =',I5,
     2                                    ' (MUST be 1 to 5).'/
     3          '         Correct Evaluated Data and Re-Run Program.')
      CALL ENDERROR
C-----ILLEGAL LOG ENERGY INTERPOLATION WITH NEGATIVE VALUES.
   80 CALL OUT9(E1,FIELD(1,1))
      CALL OUT9(E2,FIELD(1,2))
      CALL OUT9(E ,FIELD(1,3))
      WRITE(OUTP,90) ((FIELD(M,J),M=1,11),J=1,3)
      WRITE(*   ,90) ((FIELD(M,J),M=1,11),J=1,3)
   90 FORMAT(///' ERROR - Interpolating Scattering Radius.'/
     1          '         Illegal Log Energy Interpolation',
     2                                  ' Using Negative Energy.'/
     3          '         Interpolation Code=',I5,
     3                                    ' (Cannot be 3 or 5).'/
     4          '         E1,E2,E=',3(11A1,1X)/
     5          '         Correct Evaluated Data and Re-Run Program.'/
     6          '         Execution Terminated.'///)
      CALL ENDERROR
C-----ZERO LENGTH ENERGY INTERVAL.
  100 CALL OUT9(E1,FIELD(1,1))
      CALL OUT9(E2,FIELD(1,2))
      WRITE(OUTP,110) ((FIELD(M,J),M=1,11),J=1,2)
      WRITE(*   ,110) ((FIELD(M,J),M=1,11),J=1,2)
  110 FORMAT(///' ERROR - Interpolating Scattering Radius.'/
     2          '         Illegal Interpolation Over Zero',
     2                                  ' Length Interval.'/
     3          '         E1,E2=',11A1,1X,11A1/
     4          '         Correct Evaluated Data and Rr-Run Program.'/
     5          '         Execution Terminated.'///)
      CALL ENDERROR
      RETURN        ! Dummy cannot be reached
      END
      SUBROUTINE SIGBW1(E)
C=======================================================================
C
C     ADD CONTRIBUTION OF ONE SECTION OF SINGLE LEVEL BREIT-WIGNER
C     PARAMETERS.
C
C     DEFINITIONS FROM ENDF-102 FORMATS AND PROCEDURES MANUAL
C     =======================================================
C     ELASTIC
C     =======
C     THE DEFINITION OF SINGLE LEVEL ELASTIC SCATTERING IN ENDF-102 IS
C     INCORRECT. HERE WE USE THE CORRECT DEFINITION,
C
C     ELASTIC = (2*L + 1)*SIN(PS)**2 +
C
C       GJ*GAM(N)*(GAM(N) - 2*GAM(T)*SIN(PS)**2 + 2*DE*SIN(2*PS))/DEN
C
C     CAPTURE
C     =======
C     CAPTURE = GAM(N)*GAM(C)/DEN
C
C     FISSION
C     =======
C     FISSION = GAM(N)*GAM(F)/DEN
C
C     DE  = (E - ER)
C     DEN = ((DE)**2 + (GAM(T)/2)**2)
C
C     SUMMED OVER ALL RESONANCES WITH THE SAME L VALUE. NOTE, THE
C     POTENTIAL CONTRIBUTION IS INCLUDED ONLY ONCE FOR EACH L VALUE
C     AND THE TERMS GJ, 2*SIN(PS)**2, 2*SIN(2*PS) ARE THE SAME FOR
C     ALL RESONANCES WITH THE SAME (L,J) VALUE. THEREFORE THE ROUTINE
C     WILL DEFINE 3 SEPERATE SUMS,
C
C     1) GAM(N)*GAM(N)/DEN
C     2) GAM(N)*GAM(T)/DEN
C     3) GAM(N)*DE/DEN
C
C     ONLY AFTER SUMMING THESE EXPRESSIONS OVER ALL RESONANCES WILL
C     THEY BE MULTIPLIER BY THE RESONANCE INDEPENDENT FACTORS TO
C     DEFINE THE CROSS SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C
C     DEFINE RESONANCE INDEPENDENT QUANTITIES.
C
      RHO2=E*RHOX2(ISECT)
      RHOP=RHOP1(ISECT)*DSQRT(E)
C-----DEFINE PENETRATION FACTOR, SHIFT FACTOR, PHASE SHIFT.
      CALL FACTS2(LVALUE(ISECT),RHO2,SF2,PF)
      CALL FACPHI(LVALUE(ISECT),RHOP,PS)
      COSPS=DCOS(PS)
      SINPS=DSIN(PS)
C-----DEFINE SIN(2*PS) AND 2*SIN(PS)**2
      SIN2PS=TWO*SINPS*COSPS
      SINPS2=TWO*SINPS*SINPS
C-----DEFINE PENETRATION FACTOR FOR COMPETITIVE WIDTH.
      IF(LRXTAB(ISECT).LE.0) GO TO 10
      RHOZ2=(E+QVALUE(ISECT))*RHOX2(ISECT)
c-----10/10/10 - switched from FACTS2 to FACTS3 - shift not used
      CALL FACTS3(LRXTAB(ISECT)-1,RHOZ2,PFC)
      GO TO 20
   10 PFC=0.0
C
C     J VALUE LOOP.
C
C-----DEFINE STATISTICAL WEIGHT FOR ALL RESONANCES WITH THE SAME (L,J).
   20 GJ=RESTAB(2,LOW)
      AJNOW = RESJTAB(LOW)
C-----INITIALIZE CONTRIBUTION OF (L,J) SEQUENCE.
      SIGNN1=ZERO
      SIGNN2=ZERO
      SIGNN3=ZERO
      SIGNGJ=ZERO
      SIGNFJ=ZERO
C
C     RESONANCE LOOP. ADD CONTRIBUTION OF ALL RESONANCES WITH SAME (L,J)
C
      DO 30 JR=LOW,LHI
C-----ONLY USE RESONANCES WITH SAME J VALUE.
      IF(DABS(AJNOW-RESJTAB(JR)).GT.0.01) GO TO 40
      ER=ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 30  ! Note, DABS to only get = 0
      GAMC=PFC*RESTAB(3,JR)
      GN=RESTAB(4,JR)
      GAMN=PF*GN
      GAMG=RESTAB(5,JR)
      GAMF=RESTAB(6,JR)
      GAMT=GAMN+GAMG+GAMF+GAMC
      GAMT2=HALF*GAMT
      DE=E-(ER+GN*(SHIFT2(JR)-SF2))
      DEN=DE*DE+GAMT2*GAMT2
      COMFAC=GAMN/DEN
      SIGNN1=SIGNN1+COMFAC*GAMN
      SIGNN2=SIGNN2+COMFAC*GAMT
      SIGNN3=SIGNN3+COMFAC*DE
      SIGNGJ=SIGNGJ+COMFAC*GAMG
      SIGNFJ=SIGNFJ+COMFAC*GAMF
C-----END OF RESONANCE LOOP.
   30 CONTINUE
      JR=LHI+1
C-----MULTIPLY CONTRIBUTION OF (L,J) SEQUENCE BY STATISTICAL WEIGHT
C-----AND ADD TO SUM...NOTE, SINPS2 IS 2*SIN(PS)**2
   40 SIGNNI=SIGNNI+GJ*(SIGNN1-SIGNN2*SINPS2+TWO*SIGNN3*SIN2PS)
      SIGNGI=SIGNGI+GJ*SIGNGJ
      SIGNFI=SIGNFI+GJ*SIGNFJ
C-----TEST FOR ANOTHER J VALUE.
      LOW=JR
      IF(LOW.LE.LHI) GO TO 20
C
C     ADD POTENTIAL SCATTERING CONTRIBUTION FOR L VALUE.
C
C-----ADD 2*(2*L+1)*[2*SIN(PS)**2] (NOTE, SINPS2 IS 2*SIN(PS)**2)
      SIGNNI=SIGNNI+POTL(ISECT)*SINPS2
      RETURN
      END
      SUBROUTINE SIGBWM(E)
C=======================================================================
C
C     ADD CONTRIBUTION OF ONE SECTION OF MULTI-LEVEL BREIT-WIGNER
C     PARAMETERS.
C
C     DEFINITIONS FROM ENDF-102 FORMATS AND PROCEDURES MANUAL
C     =======================================================
C     CAPTURE AND FISSION CROSS SECTIONS ARE CALCULATED EXACTLY AS IN
C     THE CASE OF SINGLE LEVEL BREIT-WIGNER RESONANCES (FOR DETAILS
C     SEE, SUBROUTINE SIGBW1). HERE WE WILL ONLY CONSIDER THE ELASTIC
C     CROSS SECTION.
C
C     THE ELASTIC CROSS SECTION IS DEFINED TO BE,
C
C     ELASTIC = GJ*(1 - U(N,N))**2
C
C                                     I*GAM(N)
C     U(N,N)= EXP(-I*2*PS)*(I + ====================)
C                               ((ER-E)-I*GAM(T)/2)
C
C                                I*GAM(N)*((ER-E)+I*GAM(T)/2)
C           = EXP(-I*2*PS)*(I + =======================================)
C                               ((ER-E)-I*GAM(T)/2)*((ER-E)+I*GAM(T)/2)
C
C                                I*GAM(N)*((ER-E)-GAM(N)*GAM(T)/2)
C           = EXP(-I*2*PS)*(I + =======================================)
C                               ((ER-E)**2+(GAM(T)/2)**2)
C
C     U(N,N) = (COS(2*PS) - I*SIN(2*PS))*((1 - R) - I*S)
C
C     R      = GAM(N)*GAM(T)/2/DEN
C     S      = GAM(N)*DE/DEN
C     GAM(N) = NEUTRON WIDTH
C     GAM(R) = TOTAL WIDTH
C     DE     = E - ER
C     DEN    = ((E - ER)**2 + GAM(T)**2)
C
C     SUMMED OVER RESONANCES FOR EACH (L,J) SEQUENCE.
C
C     FROM OUR DEFINITION OF A UNIFORM TREATMENT OF ALL FORMALISMS WE
C     CAN IMMEDIATELY IDENTIFY,
C
C     X =  R
C     Y =  S
C
C     AND SUBSTITUTE INTO OUR GENERAL EXPRESSION,
C
C     ELASTIC =GJ*(2*SIN(PS)**2 - R)**2 + (SIN(2*PS) + S)**2)
C
C     ABSORPTION (CAPTURE AND FISSION) ARE TREATED IN A SINGLE LEVEL
C     APPROXIMATION,
C
C     CAPTURE =GJ*GAM(N)*GAM(G)/DEN
C
C     FISSION =GJ*GAM(N)*GAM(F)/DEN
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C
C     DEFINE RESONANCE INDEPENDENT QUANTITIES.
C
      RHO2=E*RHOX2(ISECT)
      RHOP=RHOP1(ISECT)*DSQRT(E)
C-----DEFINE PENETRATION FACTOR, SHIFT FACTOR, PHASE SHIFT.
      CALL FACTS2(LVALUE(ISECT),RHO2,SF2,PF)
      CALL FACPHI(LVALUE(ISECT),RHOP,PS)
      COSPS=DCOS(PS)
      SINPS=DSIN(PS)
C-----DEFINE SIN(2*PS) AND 2*SIN(PS)**2
      SIN2PS=TWO*SINPS*COSPS
      SINPS2=TWO*SINPS*SINPS
C-----DEFINE PENETRATION FACTOR FOR COMPETITIVE WIDTH.
      IF(LRXTAB(ISECT).LE.0) GO TO 10
      RHOZ2=(E+QVALUE(ISECT))*RHOX2(ISECT)
c-----10/10/10 - switched from FACTS2 to FACTS3 - shift not used
      CALL FACTS3(LRXTAB(ISECT)-1,RHOZ2,PFC)
      GO TO 20
   10 PFC=0.0
C
C     J VALUE LOOP.
C
C-----DEFINE STATISTICAL WEIGHT FOR ALL RESONANCES WITH THE SAME (L,J).
   20 GJ=RESTAB(2,LOW)
      AJNOW = RESJTAB(LOW)
C-----INITIALIZE CONTRIBUTION OF (L,J) SEQUENCE.
      SIGNN1=ZERO
      SIGNN2=ZERO
      SIGNGJ=ZERO
      SIGNFJ=ZERO
C
C     RESONANCE LOOP. ADD CONTRIBUTION OF ALL RESONANCES WITH SAME (L,J)
C
      DO 30 JR=LOW,LHI
C-----ONLY USE RESONANCES WITH SAME J VALUE.
      IF(DABS(AJNOW-RESJTAB(JR)).GT.0.01) GO TO 40
      ER=ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 30  ! Note, DABS to only get = 0
      GAMC=PFC*RESTAB(3,JR)
      GN=RESTAB(4,JR)
      GAMN=PF*GN
      GAMG=RESTAB(5,JR)
      GAMF=RESTAB(6,JR)
      DE=E-(ER+GN*(SHIFT2(JR)-SF2))
      GAMT2=HALF*(GAMN+GAMG+GAMF+GAMC)
      DEN=DE*DE+GAMT2*GAMT2
      COMFAC=GAMN/DEN
      SIGNN1=SIGNN1+COMFAC*GAMT2
      SIGNN2=SIGNN2+COMFAC*DE
      SIGNGJ=SIGNGJ+COMFAC*GAMG
      SIGNFJ=SIGNFJ+COMFAC*GAMF
C-----END OF RESONANCE LOOP.
   30 CONTINUE
      JR=LHI+1
C-----DEFINE NEUTRON CROSS SECTION, MULTIPLY BY STATISTICAL WEIGHT AND
C-----ADD CONTRIBUTION OF (L,J) SEQUENCE TO SUM.
   40 SIGNNI=SIGNNI+GJ*((SINPS2-SIGNN1)**2+(SIN2PS+SIGNN2)**2)
      SIGNGI=SIGNGI+GJ*SIGNGJ
      SIGNFI=SIGNFI+GJ*SIGNFJ
C-----TEST FOR ANOTHER J VALUE.
      LOW=JR
      IF(LOW.LE.LHI) GO TO 20
C
C     ADD MISSING SEQUENCES TO POTENTIAL
C
C-----ADD 2*GJ*[2*SIN(PS)**2] (NOTE, SINPS2 IS 2*SIN(PS)**2)
      SIGNNI=SIGNNI+ADDL(ISECT)*SINPS2
      RETURN
      END
      SUBROUTINE SIGRM1(E)
C=======================================================================
C
C     ADD CONTRIBUTION OF ONE SECTION OF REICH-MOORE PARAMETERS.
C
C     DEFINITIONS FROM ENDF-102 FORMATS AND PROCEDURES MANUAL
C     =======================================================
C     THE CROSS SECTIONS ARE DEFINED TO BE,
C
C     TOTAL        =2*GJ*REAL(I - U(N,N))
C     ABSORPTION   =4*GJ*(REAL(RHO(N,N)) - RHO(N,N)**2)
C                  =  GJ*(I - U(N,N)**2)
C     ELASTIC      =  GJ*(I - U(N,N))**2
C     FISSION      =4*GJ*(SUM OVER C)(RHO(N,C)**2)
C     CAPTURE      = ABSORPTION - FISSION
C
C     WHICH ARE COMPLETELY DEFINED IN TERMS OF U(N,N) AND RHO(N,C),
C
C     RHO(N,C)     =I - INVERSE(I - K)
C     U(N,N)       =EXP(-I*2*PS)*(2*INVERSE(I - K) - I)
C                  =(COS(2*PS) - I*SIN(2*PS))*(2*INVERSE(I - K) - I)
C                  =(COS(2*PS) - I*SIN(2*PS))*(1 - 2*RHO(N,N))
C
C                  =COS(2*PS)*(I-2*REAL(RHO)) - I*2*SIN(2*PS)*IM(RHO)
C
C     REAL(U(N,N)) =   COS(2*PS)*(I-2*REAL(RHO))
C     IM(U(N,N))   =-2*SIN(2*PS)*IM(RHO)
C
C     MATRIX ELEMENTS
C     =======================================================
C
C                     I*SQRT(GAM(C)/2*GAM(C*)/2)
C     (I - K)    =I - ==============================================
C                     ((ER-E) -I*(GAM(R)/2)
C
C                      I*SQRT(GAM(C)/2*GAM(C*)/2)*((ER-E) + I*(GAM(R)/2)
C                =I - ==================================================
C                     ((ER-E) -I*(GAM(R)/2)*((ER-E) + I*(GAM(R)/2)
C
C                     SQRT(GAM(C)/2*GAM(C*)/2)*(GAM(R)/2)
C                =I + ==============================================
C                     ((ER-E)**2 + (GAM(R)/2)**2)
C
C                     SQRT(GAM(C)/2*GAM(C*)/2)*(ER-E)
C                =  - ==============================================
C                     ((ER-E)**2 + (GAM(R)/2)**2)
C
C     (I - K)      = (R + I) - I*S
C
C     R            = SQRT(GAM(C)/2*GAM(C*)/2)*(GAM/2)/DEN
C     S            = SQRT(GAM(C)/2*GAM(C*)/2)*(ER-E)/DEN
C     GAM(R)       = ELIMINATED RADIATIVE WIDTH
C     GAM(C)       = PARTIAL WIDE FOR CHANNEL C FOR A RESONANCE
C     DEN          = ((ER - E)**2 + (GAM/2)**2)
C
C     SUMMED OVER RESONANCES FOR EACH (L,J) SEQUENCE.
C
C     PHYSICALLY (R) IS THE SYMMETRIC CONTRIBUTION OF THE RESONANCES
C     AND IS ALWAYS POSITIVE. SIMILARLY (S) IS THE ANTI-SYMMETRIC
C     CONTRIBUTION OF THE RESONANCES AND IS NEGATIVE FOR ENERGIES
C     LESS THAN THE RESONANCE AND POSITIVE FOR ENERGIES GREATER THAN
C     THE RESONANCE.
C
C     NOTE, IN ENDF-102 THE 2 VARIABLES DESCRIBED CORRESPOND TO (R+I)
C     AND (S). IN ORDER TO IMPROVE NUMERICAL STABILITY AND SIMPLIFY
C     THE RESULTS THE FOLLOWING EQUATIONS WILL BE DEFINED IN TERMS OF
C     (R) AND (S), SO THAT THE SOLUTION CAN BE DEFINED DIRECTLY IN
C     TERMS OF THE SYMMETRIC AND ANTI-SYMMETRIC RESONANCE CONTRIBUTIONS.
C
C     SQUARES, SUCH AS RHO(N,N)**2 AND (I - U(N,N))**2, MUST BE
C     EVALUATED IN TERMS OF THE VARIABLE AND COMPLEX CONGEGATE, E.G.,
C     F(X)    =  A + B*I
C     F(X)**2 = (A + B*I)*(A - B*I) = A**2 + B**2
C
C     SO THAT,
C     RHO(N,N)**2     = (REAL(RHO(N,N))**2 + (IM(RHO(N,N))**2
C     (1 - U(N,N))**2 = (REAL(I - U(N,N))**2) + (IM(U(N,N))**2
C                     = I - 2*REAL(U(N,N))
C                     + (REAL(U(N,N))**2)+(IM(U(N,N))**2)
C                     = I - 2*REAL(U(N,N)
C                     + (U(N,N))**2
C     SOLUTION
C     =======================================================
C     INVERTING THE MATRIX
C     =======================================================
C     SOLVING FOR R' (NOT R' + I) AND S'
C     ==================================
C     I = [(R+I)+iS][(R'+I)+iS']
C     1) I = (R+I)(R'+I) - S*S'
C     2) 0 = (R+I)S' + (R'+I)S
C
C     2) (R+I)S'  = -(R'+I)S
C     2)       S' = -(R'+I)S/(R+I)
C
C     1) I = (R+I)(R'+I) + S*(R'+I)S/(R+I)
C       (R+I) = (R+I)^2(R'+I) + S^2(R'+I)
C       (R+I) = (R'+I)[(R+I)^2+S^2] ....................
C       (R+I) = [(R+I)^2+S^2] + R'[(R+I)^2+S^2]        .
C       (R+I) = [R^2+2R+I+S^2] + R'[(R+I)^2+S^2]       .
C      -[R^2+R+S^2] = R'[(R+I)^2+S^2]                  .
C      -[R(R+I)+S^2] = R'[(R+I)^2+S^2]                 .
C                                                      .
C     R' = -[R(R+I)+S^2]/[(R+I)^2+S^2]                 .
C                                                      .
C     S' = -(R'+I)S/(R+I) ..............................
C
C     S' = -S/[(R+I)^2+S^2]
C
C     IF THERE IS NO FISSION THEN R, S, RI AND SI ARE MERELY SCALARS
C     AND THE INVERSION IS ACCOMPLISHED AS DESCRIBED HERE.
C
C     IF THERE IS FISSION R, S, RI AND SI ARE 3 X 3 MATRICES AND THE
C     INVERSION IS ACCOMPLISHED BY MATRIX INVERSION.
C
C     IN EITHER CASE THE ABOVE RELATIONSHIPS ARE AT LEAST SYMBOLICALLY
C     CORRECT AND CAN BE USED TO DEFINE USEFUL PROPERTIES OF THE
C     INVERSE MATRIX. AT LEAST SYMBOLICALLY WE CAN TREAT THESE TERMS
C     AS A LINEAR SYSTEM OF 2 EQUATIONS IN 2 UNKNOWNS AND IMMEDIATELY
C     WRITE THE SOLUTION IN THE FORM,
C
C     RI      = -(R*(R+I)+S*S)/((R+I)*(R+I)+S*S)
C     SI      =  ((R+I)*S - S*R)/((R+I)*(R+I)+S*S)
C             = -S/((R+I)*(R+I)+S*S)
C
C     DEFINITION OF TERMS APPEARING IN THE CROSS SECTIONS
C     =======================================================
C     RHO(N,C)       =I - INVERSE(I - K)
C                    =I - ((RI+I)-I*SI)
C                    =-(RI - I*SI)
C     REAL(RHO(N,C)) =-RI
C     IM(RHO(N,C))   = SI
C     RHO(N,C)**2    =(RI)**2 + (SI)**2
C
C     U(N,N)         =(COS(2*PS) - I*SIN(2*PS))*(2*INVERSE(I - K) - I)
C                    =(COS(2*PS)-I*SIN(2*PS))*(2*INVERSE(I-K)-1) + I)
C                    =(COS(2*PS)-I*SIN(2*PS))*(1 - 2*RHO(N,N))
C                    =(COS(2*PS)-I*SIN(2*PS))*(1 + 2*RI - I*2*SI)
C
C                    =COS(2*PS)*(1 + 2*RI) + SIN(2*PS)*2*SI
C                    -I*(COS(2*PS)*2*SI - SIN(2*PS)*(1 + 2*RI))
C
C     REAL(U(N,N))   =COS(2*PS)*(1 + 2*RI) + SIN(2*PS)*2*SI
C     IM(U(N,N))     =-(2*COS(2*PS)*SI - SIN(2*PS)*(1 + 2*RI))
C     (U(N,N))**2    =REAL(U(N,N))**2 + IM(U(N,N))**2
C     REAL(U(N,N))**2=(COS(2*PS)**2)*(((1 + 2*RI)**2)
C                    +4*(SIN(2*PS)**2)*(SI**2)
C                    +4*COS(2*PS)*SIN(2*PS)*(1 + 2*RI)*SI
C     IM(U(N,N))**2  =(SIN(2*PS)**2)*(((1 + 2*RI)**2)
C                    +4*(COS(2*PS)**2)*(SI**2)
C                    -4*COS(2*PS)*SIN(2*PS)*(1 + 2*RI)*SI
C
C     NOTE,           4*COS(2*PS)*SIN(2*PS)*(1 + 2*RI)*SI CANCELS OUT
C                     AND COS(2*PS)**2 + SIN(2*PS)**2 = 1
C
C     (U(N,N))**2    =((1 + 2*RI)**2 +4*(SI**2)
C                    = 1 + 4*(RI + RI**2 +SI**2)
C
C     DEFINITION OF THE CROSS SECTIONS
C     =======================================================
C     TOTAL
C     =======================================================
C     IN THIS PROGRAM THE TOTAL IS DEFINED TO BE THE SUM OF ITS PART,
C
C     TOTAL = ELASTIC + ABSORPTION
C
C     HOWEVER FOR COMPLETENESS WE WILL DEFINE IT HERE.
C
C     TOTAL        =2*GJ*REAL(1 - U(N,N))
C                  =2*GJ*(1 -((COS(2*PS)*(2*RI+1)+SIN(2*PS)*2*SI)))
C                  =2*GJ*(1 - COS(2*PS)-2*(COS(2*PS)*RI-SIN(2*PS)*SI))
C                  =4*GJ*(SIN(PS)**2 - (COS(2*PS)*RI-SIN(2*PS)*SI))
C
C     HERE WE HAVE USED THE IDENTITY 1 - COS(2*PS) = 2*SIN(PS)**2
C
C     ABSORPTION
C     =======================================================
C     ABSORPTION   =4*GJ*(REAL(RHO(N,N) - RHO(N,N)**2)
C                  =4*GJ*(-RI           - (RI**2 + SI**2))
C                  =-4*GJ*(RI + RI**2 + SI**2))
C
C     FISSION
C     =======================================================
C     FISSION      =SUM RHO(N,C)**2
C
C     WHERE RHO(N,C) ARE THE OFF DIAGONAL TERMS OF THE MATRICES.
C     WRITTEN EXPLIVITLY,
C
C     FISSION      =RHO(1,2)**2 + RHO(1,3)**2
C                  =RI(1,2)**2 + RI(1,3)**2 + SI(1,2)**2 + SI(1,3)**2
C
C     CAPTURE
C     =======================================================
C     CAPTURE IS DEFINED TO BE,
C
C     CAPTURE      = ABSORPTION - FISSION
C
C     WHEN THERE IS NO FISSION CAPTURE IS EQUAL TO ABSORPTION.
C
C     ELASTIC
C     =======================================================
C     ELASTIC      =  GJ*(1 - U(N,N))**2
C
C     NOTE THAT THIS DEFINITION IS EXACTLY THE SAME AS IN THE CASE
C     OF MULTI-LEVEL BREIT-WIGNER RESONANCES, IF WE USE (RI) AND (SI)
C     INSTEAD OF (R) AND (S) (SEE, SUBROUTINE SIGBWM FOR DETAILS).
C     THEREFORE WE SHOULD NOT BE SURPRISED TO FIND A SIMILAR
C     EXPRESSION FOR THE ELASTIC CROSS SECTION,
C
C     ELASTIC=GJ*(1 - U(N,N))**2
C            =GJ*(1-2*REAL(U(N,N))+(REAL(U(N,N))**2+IM(U(N,N))**2)
C            =GJ*(1-2*(COS(2*PS)*(RI+1)-SIN(2*PS)*SI)+(RI+1)**2+SI**2)
C            =GJ*(1-2*(COS(2*PS)*(RI+1)-SIN(2*PS)*SI)+(RI+1)**2+SI**2)
C
C     BY ADDING AND SUBTRACTING TERMS THIS CAN BE WRITTEN AS THE SUM OF
C     THE SQUARE OF 2 TERMS,
C
C     ELASTIC=GJ*((COS(2*PS)-(RI+1))**2 + (SIN(2*PS)+SI)**2)
C
C            =GJ*(COS(2*PS)**2 -2*(RI+1)*COS(2*PS) + (RI+1)**2
C                +SIN(2*PS)**2 +2*SI*SIN(2*PS)     + SI**2)
C             ========================================================
C            =GJ*(1-2*(COS(2*PS)*(RI+1)-SIN(2*PS)*SI)+(RI+1)**2+SI**2)
C
C     HERE WE HAVE USED THE IDENTITY COS(2*PS)**2+SIN(2*PS)**2 = 1
C
C     THIS FORM CAN BE FURTHER SIMPLIFIED TO AVOID ROUND-OFF,
C
C     (COS(2*PS)-(RI+1))**2 = ((COS(2*PS)-1)-RI)**2
C                          = (-2*SIN(PS)**2-RI)**2
C                          = (2*SIN(PS)**2+RI)**2
C     TO FIND,
C
C     ELASTIC  = GJ*(2*SIN(PS)**2+RI)**2 + (SIN(2*PS)+SI)**2)
C
C    THIS IS THE FORM IN WHICH THE ELASTIC CROSS SECTION IS CALCULATED,
C
C    POTENTIAL CROSS SECTION
C    =======================
C    FAR FROM RESONANCES AND FOR (L,J) SEQUENCES WHICH DO NOT HAVE ANY
C    RESONANCES (RI) AND (SI) WILL BE SMALL OR ZERO AND WE FIND THAT THE
C    ELASTIC CROSS SECTION REDUCES TO THE POTENTIAL CROSS SECTION,
C
C    ELASTIC  =GJ*(4*SIN(PS)**4 + SIN(2*PS)**2)
C
C    USING THE IDENTITY SIN(2*PS) = 2*SIN(PS)*COS(PS)
C
C             =GJ*(4*SIN(PS)**4 + 4*(SIN(PS)*COS(PS))**2)
C             =4*GJ*SIN(PS)**2*(SIN(PS)**2 + COS(PS)**2)
C
C    BUT SIN(PS)**2 + COS(PS)**2 = 1, TO FIND,
C
C    ELASTIC  =4*GJ*SIN(PS)**2
C
C    IT IS IMPORTANT TO REALIZE THAT EVEN IN THE CASE WHERE THERE ARE
C    NO RESONANCES SPECIFIED FOR A (L,J) SEQUENCE THE CROSS SECTION
C    DOES NOT BECOME ZERO, SINCE THERE IS STILL A CONTRIBUTION TO THE
C    POTENTIAL CROSS SECTION.
C
C     NUMERICAL STABILITY
C     =======================================================
C     OBVIOUSLY, PHYSICALLY THE TOTAL, ELASTIC AND ABSORPTION CANNOT
C     BE NEGATIVE, AND ALL OF THE ABOVE EQUATIONS REFLECT THIS FACT.
C     HOWEVER, CARE MUST BE USED TO AVOID NUMERICAL INSTABILITY.
C
C     FROM THE DEFINITION OF ABSORPTION,
C
C     ABSORPTION   =-4*GJ*(RI + RI**2 + SI**2)
C
C     BUT BY THE ABOVE DEFINITIONS,
C
C     RI      = -(R*(R+I)+S*S)/((R+I)*(R+I)+S*S)
C     SI      = -S/((R+I)*(R+I)+S*S)
C
C     RI      = -((R+I)*(R+I)+S*S - (R+I))/((R+I)*(R+I)+*S*S)
C             = (R+I)/((R+I)*(R+I)+S*S) - I
C     RI**2   = (R+I)**2/((R+I)**2+S**2)**2)-2*(R+I)/((R+I)**2+S**2)+1
C     SI**2   =     S**2/((R+I)**2+S**2)**2)
C
C     RI**2+SI**2 = (I - 2*(R+I) + ((R+I)**2+S**2))/((R+I)**2+S**2)
C     RI          = (  +   (R+I) - ((R+I)**2+S**2))/((R+I)**2+S**2)
C     =============================================================
C     SUM         = (I -   (R+I)                  )/((R+I)**2+S**2)
C                 =-R/((R+I)**2+S**2)
C
C     ABSORPTION  =4*GJ*R/((R+I)*(R+I)+S*S)
C
C     FROM THE DEFINITION OF R WE CAN SEE THAT THIS IS MERELY THE
C     SYMMETRIC CONTRIBUTION OF THE RESONANCES, WHICH IS INHERENTLY
C     POSITIVE AND NOT SUBJECT TO ROUND-OFF ERROR. THEREFORE THIS
C     DEFINITION OF THE ABSORPTION WILL BE INHERENTLY STABLE IF WE
C     AVOID PROBLEMS BY INITIALLY DEFINING (R) RATHER THAN (R+I) AND
C     USE THIS QUANTITY TO DEFINE THE ABSORPTION (NOTE, IF WE FIRST TRY
C     TO DEFINE (R+I) IT INVOLVES ADDING I WHICH INTRODUCES ROUND-OFF
C     ERROR FAR FROM RESONANCES - IF WE SUBSEQUENTLY TRY TO DEFINE R
C     BY SUBTRACTING I WE CAN MERELY INTRODUCE MORE ROUND-OFF - HOWEVER
C     IF WE DEFINE THE QUANTITY (R) DIRECTLY WE DO NOT INTRODUCE ANY
C     ROUND-OFF).
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FISSY/LFWX,LFI,MT451,LFWSUM
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C----- 07/26/09 - EXPANDED DIMENSION FROM 3 X 3 TO 10 X 10
      DIMENSION R(10,10),S(10,10),RI(10,10),SI(10,10)
      EQUIVALENCE (RI(1,1),RI11),(SI(1,1),SI11),(RI(1,2),RI12),
     1 (SI(1,2),SI12),(RI(1,3),RI13),(SI(1,3),SI13),
     2 (R(1,1),R11),(S(1,1),S11)
C
C     DEFINE RESONANCE INDEPENDENT QUANTITIES.
C
      RHO2=E*RHOX2(ISECT)
      RHOP=RHOP1(ISECT)*DSQRT(E)
C-----DEFINE PENETRATION FACTOR, SHIFT FACTOR, PHASE SHIFT.
      CALL FACTS3(LVALUE(ISECT),RHO2,PF)
      CALL FACPHI(LVALUE(ISECT),RHOP,PS)
      COSPS=DCOS(PS)
      SINPS=DSIN(PS)
C-----DEFINE SIN(2*PS) AND 2*SIN(PS)**2
      SIN2PS=TWO*SINPS*COSPS
      SINPS2=TWO*SINPS*SINPS
C-----IS THIS SECTUION FISSIONABLE.
      IF(LFWL(ISECT).NE.0) GO TO 40
C-----------------------------------------------------------------------
C
C     NO FISSION.
C
C     J VALUE LOOP.
C
C-----------------------------------------------------------------------
C-----DEFINE STATISTICAL WEIGHT FOR ALL RESONANCES WITH THE SAME (L,J).
   10 GJ=RESTAB(2,LOW)
      AJNOW = RESJTAB(LOW)
C-----INITIALIZE ELEMENTS.
      R11=ZERO
      S11=ZERO
C
C     RESONANCE LOOP. ADD CONTRIBUTION OF ALL RESONANCES WITH SAME (L,J)
C
      DO 20 JR=LOW,LHI
C-----ONLY USE RESONANCES WITH SAME STATISTICAL WEIGHT, I.E. J VALUE.
      IF(DABS(AJNOW-RESJTAB(JR)).GT.0.01) GO TO 30
      ER = ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 20  ! Note, DABS to only get = 0
      GAMN=PF*RESTAB(3,JR)
      GAMNG2=RESTAB(4,JR)
      DE = ER - E
      DEN=DE*DE+GAMNG2*GAMNG2
      DE2=DE/DEN
      GAMNG4=GAMNG2/DEN
C-----DEFINE TERMS IN UPPER TRIANGLE OF MATRIX.
      R11=R11+GAMNG4*GAMN
      S11=S11+DE2*GAMN
C-----END OF RESONANCE LOOP.
   20 CONTINUE
      JR=LHI+1
   30 DET = (R11+ONE)**2+S11**2
      SI11=-S11/DET
      RI11=-(R11*(R11+ONE)+S11**2)/DET
C
C     ADD CONTRIBUTION OF CURRENT (L,J) SEQUENCE TO SUM.
C
C-----ABSORPTION (IN THIS CASE ONLY CAPTURE).
      SIGNG1=-FOUR*GJ*(RI11+(RI11**2+SI11**2))
C-----ELASTIC.
      SIGNNI=SIGNNI+GJ*((SINPS2+TWO*RI11)**2+(SIN2PS+TWO*SI11)**2)
C-----CAPTURE (ABSORPTION).
      SIGNGI=SIGNGI+SIGNG1
C-----TEST FOR ANOTHER J VALUE.
      LOW=JR
      IF(LOW.LE.LHI) GO TO 10
C
C     ADD MISSING SEQUENCES TO POTENTIAL
C
C-----ADD 2*GJ*[2*SIN(PS)**2] (NOTE, SINPS2 IS 2*SIN(PS)**2)
      SIGNNI=SIGNNI+ADDL(ISECT)*SINPS2
      RETURN
C-----------------------------------------------------------------------
C
C     FISSION
C
C     J VALUE LOOP.
C
C-----------------------------------------------------------------------
C-----DEFINE STATISTICAL WEIGHT FOR ALL RESONANCES WITH THE SAME (L,J).
   40 GJ=RESTAB(2,LOW)
      AJNOW = RESJTAB(LOW)
C-----INITIALIZE MATRICES.
      DO 50 I=1,3
      DO 50 J=1,3
      R (I,J)=ZERO
   50 S (I,J)=ZERO
C
C     RESONANCE LOOP. ADD CONTRIBUTION OF ALL RESONANCES WITH SAME (L,J)
C
      DO 60 JR=LOW,LHI
C-----ONLY USE RESONANCES WITH SAME STATISTICAL WEIGHT, I.E. J VALUE.
      IF(DABS(AJNOW-RESJTAB(JR)).GT.0.01) GO TO 70
      ER = ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 60  ! Note, DABS to only get = 0
      GAMN=PF*RESTAB(3,JR)
      GAMN2=DSQRT(GAMN)
      GAMNG2=RESTAB(4,JR)
      GAMFA2=RESTAB(5,JR)
      GAMFA=GAMFA2*GAMFA2
      GAMFB2=RESTAB(6,JR)
      GAMFB=GAMFB2*GAMFB2
      DE = ER - E
      DEN=DE*DE+GAMNG2*GAMNG2
      DE2=DE/DEN
      GAMNG4=GAMNG2/DEN
      GAMNFA=GAMN2*GAMFA2
      GAMNFB=GAMN2*GAMFB2
      GAMFAB=GAMFA2*GAMFB2
C-----DEFINE TERMS IN UPPER TRIANGLE OF MATRIX.
      R(1,1)=R(1,1)+GAMNG4*GAMN
      R(1,2)=R(1,2)+GAMNG4*GAMNFA
      R(1,3)=R(1,3)+GAMNG4*GAMNFB
      R(2,2)=R(2,2)+GAMNG4*GAMFA
      R(2,3)=R(2,3)+GAMNG4*GAMFAB
      R(3,3)=R(3,3)+GAMNG4*GAMFB
      S(1,1)=S(1,1)+DE2*GAMN
      S(1,2)=S(1,2)+DE2*GAMNFA
      S(1,3)=S(1,3)+DE2*GAMNFB
      S(2,2)=S(2,2)+DE2*GAMFA
      S(2,3)=S(2,3)+DE2*GAMFAB
      S(3,3)=S(3,3)+DE2*GAMFB
C-----END OF RESONANCE LOOP.
   60 CONTINUE
      JR=LHI+1
C-----MAKE MATRICES SYMMETRIC.
   70 R(2,1)=R(1,2)
      S(2,1)=S(1,2)
      R(3,1)=R(1,3)
      S(3,1)=S(1,3)
      R(3,2)=R(2,3)
      S(3,2)=S(2,3)
C-----INVERT COMPLEX MATRIX.
C----- 2/20/10 - switch to general form if needed.
c----- 1/31/17 - Not needed = checked - both give EXACTLY same answer
      CALL FROBNS3(R,S,RI,SI)
C----- 2/20/10 - switch to general form if needed.
C
C     ADD CONTRIBUTION OF CURRENT (L,J) SEQUENCE TO SUM.
C
C-----(SEE, EQULVALENCE FOR DEFINITION OF RI11, SI11, RI12, SI12, RI13,
C----- SI13).
      GJ4=FOUR*GJ
C-----ABSORPTION (IN THIS CASE CAPTURE + FISSION).
      SIGNG1=-GJ4*(RI11+(RI11**2+SI11**2))
C-----ELASTIC.
      SIGNNI=SIGNNI+GJ*((SINPS2+TWO*RI11)**2+(SIN2PS+TWO*SI11)**2)
C-----FISSION.
      SIGNF1=GJ4*(RI12**2+RI13**2+SI12**2+SI13**2)
C-----CAPTURE (ABSORPTION - FISSION).
      SIGNGI=SIGNGI+(SIGNG1-SIGNF1)
C-----FISSION.
      SIGNFI=SIGNFI+SIGNF1
C-----TEST FOR ANOTHER J VALUE.
      LOW=JR
      IF(LOW.LE.LHI) GO TO 40
C
C     ADD MISSING SEQUENCES TO POTENTIAL
C
C-----ADD 2*GJ*[2*SIN(PS)**2] (NOTE, SINPS2 IS 2*SIN(PS)**2)
      SIGNNI=SIGNNI+ADDL(ISECT)*SINPS2
      RETURN
      END
      SUBROUTINE SIGAA4(E)
C=======================================================================
C
C     ADD CONTRIBUTION OF ONE SECTION OF ADLER-ADLER PARAMETERS.
C     THIS WILL INCLUDE THE BACKGROUND CROSS SECTION PLUS ALL
C     RESONANCES, FOR ALL VALUES OF L AND J.
C
C     THIS ROUTINE USES THE ENDF/B-IV AND EARLIER DEFINITIONS OF THE
C     ADLER-ADLER PARAMETERS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DATA COMFC1/0.0D+00/
C-----DEFINE ARITHMETIC STATEMENT FUNCTION FOR CALCULATION OF BACKGROUND
      BACKGR(JR)=((((RESTAB(4,JR)/E+RESTAB(3,JR))/E+RESTAB(2,JR))/E+
     1 RESTAB(1,JR))+(RESTAB(6,JR)*E+RESTAB(5,JR))*E)
C-----DEFINE ARITHMETIC STATEMENT FUNCTION FOR CALCULATION OF RESONANCE.
      RESER(GR,GI)=COMFC1*(GR*COSPS+GI*SINPS)+COMFC2*(GI*COSPS-GR*SINPS)
C
C     DEFINE RESONANCE INDEPENDENT QUANTITIES.
C
      ESQRT=DSQRT(E)
C-----TYPE OF BACKGROUND (NOT L VALUE).
      LI=LVALUE(ISECT)
C-----ADLER-ADLER PARAMETERS ARE ONLY DEFINED FOR L=0, IN WHICH CASE
C-----THE PHASE SHIFT IS EQUAL TO RHO (NO NEED TO CALL FACPHI).
      RHOP=RHOP1(ISECT)*DSQRT(E)
      COSPS=DCOS(RHOP)
      SINPS=DSIN(RHOP)
C
C     DEFINE BACKGROUND CROSS SECTIONS.
C
      GO TO (10,20,30,40,50,60,70),LI
   10 SIGNTI=BACKGR(LOW)
      GO TO 80
   20 SIGNFI=BACKGR(LOW)
      GO TO 80
   30 SIGNTI=BACKGR(LOW)
      SIGNFI=BACKGR(LOW+1)
      GO TO 80
   40 SIGNGI=BACKGR(LOW)
      GO TO 80
   50 SIGNTI=BACKGR(LOW)
      SIGNGI=BACKGR(LOW+1)
      GO TO 80
   60 SIGNFI=BACKGR(LOW)
      SIGNGI=BACKGR(LOW+1)
      GO TO 80
   70 SIGNTI=BACKGR(LOW)
      SIGNFI=BACKGR(LOW+1)
      SIGNGI=BACKGR(LOW+2)
C
C     RESONANCE LOOP.
C
C-----DEFINE BASE ADDRESS FOR RESONANCE PARAMETERS.
   80 LOW=LOW+3
      DO 90 JR=LOW,LHI
      ER = ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 90 ! Note, DABS to only get = 0
C-----CONSTRAIN FIT TO USE SAME RESONANCE ENERGY AND TOTAL WIDTH FOR
C-----ALL REACTIONS (SEE ENDF/B-102 PROCEDURES MANUAL).
      DET=ER-E
      DWT=RESTAB(2,JR)
      DEN=DET*DET+DWT*DWT
      COMFC1=DWT/DEN
      COMFC2=DET/DEN
      SIGNTI=SIGNTI+RESER(RESTAB(3 ,JR),RESTAB( 4,JR))
      SIGNFI=SIGNFI+RESER(RESTAB(7 ,JR),RESTAB( 8,JR))
      SIGNGI=SIGNGI+RESER(RESTAB(11,JR),RESTAB(12,JR))
C-----END OF RESONANCE LOOP.
   90 CONTINUE
C-----MULTIPLY ALL BY THE SQUARE ROOT OF E AND ADD TERM TO TOTAL.
      SIGNTI=ESQRT*SIGNTI+TWO*(ONE-COSPS)
      SIGNFI=ESQRT*SIGNFI
      SIGNGI=ESQRT*SIGNGI
C-----DEFINE ELASTIC BY SUBTRACTION.
      SIGNNI=SIGNTI-(SIGNFI+SIGNGI)
      RETURN
      END
      SUBROUTINE SIGAA5(E)
C=======================================================================
C
C     ADD CONTRIBUTION OF ONE SECTION OF ADLER-ADLER PARAMETERS.
C     THIS WILL INCLUDE THE BACKGROUND CROSS SECTION PLUS ALL
C     RESONANCES, FOR ALL VALUES OF L AND J.
C
C     THIS ROUTINE USES THE ENDF/B-V AND VI DEFINITIONS OF THE
C     ADLER-ADLER PARAMETERS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C-----DEFINE ARITHMETIC STATEMENT FUNCTION FOR CALCULATION OF BACKGROUND
      BACKGR(JR)=((((RESTAB(4,JR)/E+RESTAB(3,JR))/E+RESTAB(2,JR))/E+
     1 RESTAB(1,JR))+(RESTAB(6,JR)*E+RESTAB(5,JR))*E)
C
C     DEFINE RESONANCE INDEPENDENT QUANTITIES.
C
      ESQRT=DSQRT(E)
C-----TYPE OF BACKGROUND (NOT L VALUE).
      LI=LVALUE(ISECT)
C-----ADLER-ADLER PARAMETERS ARE ONLY DEFINED FOR L=0, IN WHICH CASE
C-----THE PHASE SHIFT IS EQUAL TO RHO (NO NEED TO CALL FACPHI).
      RHOP=RHOP1(ISECT)*ESQRT
      COSPS=DCOS(RHOP)
      SINPS=DSIN(RHOP)
C-----INITIALIZE SUMS.
      SIGNT1=ZERO
      SIGNT2=ZERO
      SIGNT3=ZERO
      SIGNT4=ZERO
      SIGNF1=ZERO
      SIGNF2=ZERO
      SIGNG1=ZERO
      SIGNG2=ZERO
C
C     DEFINE BACKGROUND CROSS SECTIONS.
C
      GO TO (10,20,30,40,50,60,70),LI
   10 SIGNTI=BACKGR(LOW)
      GO TO 80
   20 SIGNFI=BACKGR(LOW)
      GO TO 80
   30 SIGNTI=BACKGR(LOW)
      SIGNFI=BACKGR(LOW+1)
      GO TO 80
   40 SIGNGI=BACKGR(LOW)
      GO TO 80
   50 SIGNTI=BACKGR(LOW)
      SIGNGI=BACKGR(LOW+1)
      GO TO 80
   60 SIGNFI=BACKGR(LOW)
      SIGNGI=BACKGR(LOW+1)
      GO TO 80
   70 SIGNTI=BACKGR(LOW)
      SIGNFI=BACKGR(LOW+1)
      SIGNGI=BACKGR(LOW+2)
C
C     RESONANCE LOOP.
C
C-----DEFINE BASE ADDRESS FOR RESONANCE PARAMETERS.
   80 LOW=LOW+3
      DO 90 JR=LOW,LHI
      ER = ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 90 ! Note, DABS to only get = 0
C-----CONSTRAIN FIT TO USE SAME RESONANCE ENERGY AND TOTAL WIDTH FOR
C-----ALL REACTIONS (SEE ENDF/B-102 PROCEDURES MANUAL).
      DET=ER-E
      DWT=RESTAB(2,JR)
      GRT=RESTAB(3,JR)
      GIT=RESTAB(4,JR)
      DEN=DET*DET+DWT*DWT
      COMFC1=DWT/DEN
      COMFC2=DET/DEN
      SIGNT1=SIGNT1+COMFC1*GRT
      SIGNT2=SIGNT2+COMFC1*GIT
      SIGNT3=SIGNT3+COMFC2*GIT
      SIGNT4=SIGNT4+COMFC2*GRT
      SIGNF1=SIGNF1+COMFC1*RESTAB( 7,JR)
      SIGNF2=SIGNF2+COMFC2*RESTAB( 8,JR)
      SIGNG1=SIGNG1+COMFC1*RESTAB(11,JR)
      SIGNG2=SIGNG2+COMFC2*RESTAB(12,JR)
C-----END OF RESONANCE LOOP.
   90 CONTINUE
C-----ADD POTENTIAL TO TOTAL, MULTIPLY BY THE SQUARE ROOT OF E.
      SIGNTI=TWO*(ONE-COSPS)+ESQRT*(SIGNTI+
     1 (SIGNT1+SIGNT3)*COSPS+(SIGNT2-SIGNT4)*SINPS)
      SIGNFI=ESQRT*(SIGNFI+SIGNF1+SIGNF2)
      SIGNGI=ESQRT*(SIGNGI+SIGNG1+SIGNG2)
C-----DEFINE ELASTIC BY SUBTRACTION.
      SIGNNI=SIGNTI-(SIGNFI+SIGNGI)
      RETURN
      END
      SUBROUTINE SIGHRF(E)
C=======================================================================
C
C     ADD CONTRIBUTION OF ONE SECTION (L,S,J) OF HYBRID R-FUNCTION
C     PARAMETERS.
C
C     EACH (L,S,J) STATE WILL BE TREATED AS A SEPARATE SECTION.
C
C     THE EQUATIONS DEFINED IN ENDF-102 SAY TO USE R-MATRIX FOR THE
C     ELASTIC AND EXACTLY SINGLE LEVEL BREIT-WIGNER FOR ALL OTHER
C     REACTIONS. THIS LEADS TO AN ODD MIXTURE TO EQUATIONS IN THE
C     FOLLOWING TREATMENT.
C
C     CAPTURE AND FISSION CROSS SECTIONS ARE CALCULATED EXACTLY AS IN
C     THE CASE OF SINGLE LEVEL BREIT-WIGNER RESONANCES (FOR DETAILS
C     SEE, SUBROUTINE SIGBW1). HERE WE WILL ONLY CONSIDER THE ELASTIC
C     CROSS SECTION.
C
C     TERMS USED TO DEFINE ELASTIC CROSS SECTION
C     ==========================================
C     WE DEFINE,
C
C     RP   = RLSJ*PLSJ
C     RR0P = RR0LSJ*PLSJ
C     RI0P = RI0LSJ*PLSJ
C
C            GAM(N)
C     RP   = ========================        + (RR0P + I*RI0P)
C            2*((ER-E) - I*(GAM(R)/2)
C
C            GAM(N)*(((ER-E)+I*GAM(R)/2)
C          = ===========================     + (RR0P + I*RI0P)
C            2*((ER-E)**2+(GAM(R)/2)**2)
C
C            GAM(N)*(ER-E)+I*GAM(N)*GAM(R)/2
C          = =============================== + (RR0P + I*RI0P)
C            2*((ER-E)**2+(GAM(R)/2)**2)
C
C     RP   = RRP   + I*RIP
C
C     RRP  =(GAM(N)*(ER-E))/DEN1         + RR0P
C     RIP  =(GAM(N)*GAM(R)/2)/DEN1       + RI0P
C     DEN1 =2*((ER-E)**2+(GAM(R)/2)**2)
C
C     SUMMED OVER RESONANCES FOR EACH (L,S,J) SEQUENCE. FOLLOWING THE
C     SUM OVER RESONANCES WE DEFINE,
C
C            1 + I*PLSJ*RLSJ   1 + I*RP   (1 - RIP) + I*RRP
C     XLSJ = =============== = ========   =================
C            1 - I*PLSJ*RLSJ   1 - I*RP   (1 + RIP) - I*RRP
C
C            ((1 - RIP) + I*RRP)*((1 + RIP) + I*RRP)
C          = =======================================
C            ((1 + RIP)**2 + (RRP)**2)
C
C            1 - RIP**2 - RRP**2 + I*2*RRP
C          = =============================
C            ((1 + RIP)**2 + (RRP)**2)
C
C                  2*(RIP+ RIP**2 + RRP**2 - I*RRP)
C          = 1 -  =================================
C                 ((1 + RIP)**2 + (RRP)**2)
C
C     XLSJ = (1 - R) - I*IS
C
C     R    = 2*(RIP+((RIP)**2+(RRP)**2)/DEN2
C     S    =-2*RRP/DEN2
C     DEN2 =(1 + RIP)**2 + (RRP)**2
C
C     DEFINITION OF THE ELASTIC CROSS SECTION
C     =======================================
C     THE ELASTIC CROSS SECTION IS DEFINED TO BE,
C
C     ELASTIC      = GJ*(1 - ULSJ)**2
C
C     ULSJ         = EXP(-I*2*PS)*XLSJ
C                  = (COS(2*PS) - I*SIN(2*PS))*((R + 1) - I*S)
C
C     ELASTIC  = GJ*((2*SIN(PS)**2+R)**2 + (SIN(2*PS)+S)**2)
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C-----TEMPORARILY DEFINE REAL AND IMAGINARY BACKGROUND TO BE ZERO.
      DATA R0R/0.0D+00/
      DATA R0I/0.0D+00/
C
C     DEFINE RESONANCE INDEPENDENT QUANTITIES.
C
      GJ=GJTAB(ISECT)
      RHO2=E*RHOX2(ISECT)
      RHOP=RHOP1(ISECT)*DSQRT(E)
C-----DEFINE PENETRATION FACTOR, SHIFT FACTOR, PHASE SHIFT.
      CALL FACTS2(LVALUE(ISECT),RHO2,SF2,PF)
      CALL FACPHI(LVALUE(ISECT),RHOP,PS)
      COSPS=DCOS(PS)
      SINPS=DSIN(PS)
C-----DEFINE SIN(2*PS) AND 2*SIN(PS)**2
      SIN2PS=TWO*SINPS*COSPS
      SINPS2=TWO*SINPS*SINPS
C-----DEFINE PENETRATION FACTORS FOR COMPETITIVE WIDTHS, IF ANY.
      LRX=LRXTAB(ISECT)
C-----INITIALIZE CONTRIBUTION OF (L,S,J) SEQUENCE.
      RR=PF*R0R
      RI=PF*R0I
      SIGNGJ=ZERO
      SIGNFJ=ZERO
C
C     RESONANCE LOOP. ADD CONTRIBUTION OF ALL RESONANCES WITH SAME
C     (L,S,J).
C
      DO 30 JR=LOW,LHI
      ER=ENRES(JR)
C-----02/14/04 - SKIP ZERO ENERGY RESONANCES
      IF(DABS(ER).LE.0.0) GO TO 30 ! Note, DABS to only get = 0
C-----DEFINE NEUTRON, CAPTURE AND FISSION WIDTHS. INITIALIZE ELIMINATED
C-----WIDTH TO CAPTURE (COMPETING REACTIONS, IF ANY, WILL BE ADDED
C-----BELOW).
      DE = ER - E
      GN=RESTAB(2,JR)
      GAMN=PF*GN
      GAMG=RESTAB(3,JR)
      GAMF=RESTAB(4,JR)
      GAMER=GAMG
C-----DEFINE COMPETITIVE WIDTHS, IF ANY AND ADD TO ELIMINATED WIDTH.
      IF(LRX.LE.0) GO TO 20
C----- 5 -  8 = WIDTHS
C----- 9 - 12 = L VALUES
C-----INITIALIZE INDEX TO WIDTH AND EXIT CHANNEL L-VALUE.
      LL1=4
      LL2=8
C-----DEFINE COMPETITIVE WIDTHS AND ADD TO ELIMINATED WIDTH.
      DO 10 I=1,LRX
      LL1=LL1+1
      LL2=LL2+1
      LCOM=RESTAB(LL2,JR)
      RHOZ2=(E+EXCITE(I,ISECT))*RHOX2(ISECT)
c-----10/10/10 - switched from FACTS2 to FACTS3 - shift not used
      CALL FACTS3(LCOM,RHOZ2,PFC)
      GAMC=PFC*RESTAB(LL1,JR)
   10 GAMER=GAMER+GAMC
C-----ADD CONTRIBUTION OF RESONANCE.
   20 GAMER2=HALF*GAMER
      DE=E-(ER+GN*(SHIFT2(JR)-SF2))
C-----R-MATRIX TREATMENT FOR ELASTIC.
      DEN=DE*DE+GAMER2*GAMER2
      COMFAC=GAMN/DEN
      RR=RR+COMFAC*DE
      RI=RI+COMFAC*GAMER
C-----BREIT-WIGNER TREATMENT FOR CAPTURE AND FISSION.
      GAMT=GAMER+GAMN+GAMF
      GAMT2=HALF*GAMT
      DEN=DE*DE+GAMT2*GAMT2
      COMFAC=GAMN/DEN
      SIGNGJ=SIGNGJ+COMFAC*GAMG
      SIGNFJ=SIGNFJ+COMFAC*GAMF
C-----END OF RESONANCE LOOP.
   30 CONTINUE
C-----DEFINE 2 TERMS FOR UNIFORM TREATMENT.
      DEN2=(ONE+RI)**2+RR**2
      RI=TWO*(RI+RI**2+RR**2)/DEN2
      RR=-TWO*RR/DEN2
C-----MULTIPLY CONTRIBUTION OF (L,S,J) SEQUENCE BY STATISTICAL WEIGHT
C-----AND ADD TO SUM (NOTE, SINPS2 = 2*SIN(PS)**2).
      SIGNNI=SIGNNI+GJ*((SINPS2+RI)**2+(SIN2PS+RR)**2)
      SIGNGI=SIGNGI+GJ*SIGNGJ
      SIGNFI=SIGNFI+GJ*SIGNFJ
C
C     ADD MISSING SEQUENCES TO POTENTIAL
C
C-----ADD 2*GJ*[2*SIN(PS)**2] (NOTE, SINPS2 IS 2*SIN(PS)**2)
      SIGNNI=SIGNNI+ADDL(ISECT)*SINPS2
      RETURN
      END
      SUBROUTINE SIGGRM(E)
C=======================================================================
C
C     GENERAL R-MATRIX FORMALISM
C     ==========================
C
C     GENERAL R-MATRIX FORMALISM HAS NOT YET BEEN IMPLEMENTED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DATA IMOK/0/
      IF(IMOK.NE.0) GO TO 20
      E=0.0
      WRITE(OUTP,10)
      WRITE(*   ,10)
   10 FORMAT(///' ERROR - General R-Matrix Formalism has NOT yet',
     1                                        ' been Implemented'/
     2          '         Execution Terminated'///)
      CALL ENDERROR
   20 RETURN
      END
      SUBROUTINE SIGURP(E)
C=======================================================================
C
C     DEFINE CROSS SECTIONS FOR A SECTION OF UNRESOLVED RESONANCE
C     (LRU=2), ENERGY DEPENDENT PARAMETERS (LRF=2).
C
C     CROSS SECTIONS ARE DEFINED BY INTERPOLATING PARAMETERS (NOT CROSS
C     SECTIONS).
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/RANGER/LOW,LHI
      COMMON/PARAMS/RHO2,RHOP,BETAE,SF2,PF,PS,COSPS,SINPS,SIGNTI,
     1 SIGNNI,SIGNGI,SIGNFI,SIGNXI
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION DUMSET(6)
      EQUIVALENCE (DUMSET(2),DX),(DUMSET(3),GXX),(DUMSET(4),GNOX),
     1 (DUMSET(5),GGX),(DUMSET(6),GFX)
C-----DEFINE L VALUE.
      LNOW=LVALUE(ISECT)
C-----DEFINE INTERPOLATION LAW, STATISTICAL WEIGHT AND DEGREES OF
C-----FREEDOM FOR EACH REACTION (CAPTURE DEGREES OF FREEDOM IS ASSUMED
C-----TO BE INFINITY).
      INTX=RESTAB(1,LOW)
      GJ=RESTAB(2,LOW)
      MUX=RESTAB(3,LOW)
      MUN=RESTAB(4,LOW)
      MUF=RESTAB(6,LOW)
      AMUN=MUN
C-----FIND ENERGY OR ENERGY RANGE WHERE PARAMETERS ARE GIVEN.
      LOWP1=LOW+1
      DO 10 L=LOWP1,LHI
      IF(E.lt.ENRES(L)) go to 20
      IF(E.eq.ENRES(L)) go to 30
   10 CONTINUE
C-----EXTEND PARAMETERS AS CONSTANT OUTSIDE THEIR TABULATED ENERGY
C-----RANGE.
      L=LHI
      GO TO 30
   20 IF(L.EQ.LOWP1) GO TO 30
C-----INTERPOLATE PARAMETERS TO ENERGY E.
      CALL TERPUP(E,DUMSET,L,INTX)
      GO TO 40
C-----DEFINE PARAMETER SET AT E (E IS AN ENERGY AT WHICH PARAMETERS
C-----IS TABULATED).
   30 DX=RESTAB(2,L)
      GXX=RESTAB(3,L)
      GNOX=RESTAB(4,L)
      GGX=RESTAB(5,L)
      GFX=RESTAB(6,L)
C-----CALCULATE PENETRABILITY (VL) AND PHASE SHIFT(PS)
   40 E2=DSQRT(E)
      RHO2=E*RHOX2(ISECT)
      RHOC=E2*RHOP1(ISECT)
CAK   CALL UNFAC(LNOW,RHO2,RHOC,AMUN,VL,PS)
      CALL UNFACpre(LNOW,RHO2,RHOC,AMUN,VL,PS)
C-----DEFINE NEUTRON WIDTH.
      GNX=GNOX*VL*E2
C-----CALCULATE FLUCTUATION INTEGRALS (RN, RC AND RF).
C-----9/20/10 - ADDED RX
      CALL GNRL3(GNX,GGX,GFX,GXX,MUN,MUF,MUX,RN,RC,RF,RX)
C-----DEFINE COMMON FACTOR FOR ALL REACTIONS.
      TEMP=PI2*GJ*GNX/DX
C-----DEFINE CROSS SECTIONS FOR THIS (L, J) STATE.
      SIGNNI=(RN*GNX-TWO*(DSIN(PS)**2))*TEMP
      SIGNGI=RC*GGX*TEMP
      SIGNFI=RF*GFX*TEMP
      SIGNXI=RX*GXX*TEMP ! OPTIONAL COMPETITION
C-----ADD L COMPONENT OF POTENTIAL SCATTERING DURING PASS THROUGH
C-----FIRST J STATE OF EACH L VALUE (INDICATED BY EITHER A NEW MODE
C-----MODE=RESONANCE REPRESENTATION, OR NEW L VALUE).
      IF(ISECT.EQ.1) GO TO 50
      IF(MODE(ISECT-1).EQ.MODE(ISECT).AND.
     1 LVALUE(ISECT-1).EQ.LVALUE(ISECT)) GO TO 60
   50 SIGNNI=SIGNNI+(EIGHT*FLOAT(LNOW)+FOUR)*((DSIN(PS))**2)
   60 RETURN
      END
      SUBROUTINE SIGURS(E)
C=======================================================================
C
C     CROSS SECTIONS ARE DEFINED BY INTERPOLATING CROSS SECTIONS (NOT
C     PARAMETERS). THIS ROUTINE IS PRESENTLY NOT IMPLEMENTED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DATA IMOK/0/
      IF(IMOK.NE.0) GO TO 20
      E=0.0
      WRITE(OUTP,10)
      WRITE(*   ,10)
   10 FORMAT(///' ERROR - Unresolved Parameter Interpolation has',
     1                                ' not yet been imlemented.'/
     2          '         Execution Terminated')
      CALL ENDERROR
   20 RETURN
      END
      SUBROUTINE FROBNS3(AM1,B,CM1,D)
C=======================================================================
C
C     3 X 3 MATRIX
C     ------------
C     STARTING FROM THE COMPLEX MATRIX (A+I)+I*B, WHERE A AND B ARE
C     SYMMETRIC, INVERT THE COMPLEX MATRIX TO DEFINE INVERSE (C+I)+I*D.
C
C     THE FROBENIUS-SCHUR METHOD IS USED TO INVERT THE COMPLEX MATRIX.
C
C     STARTING FROM THE DEFINITION OF THE COMPLEX MATRIX AND ITS INVERSE
C
C     IDENTITY = ((A+I)+I*B)*((C+I)+I*D)
C
C     WE OBTAIN 2 REAL MATRIX EQUATIONS
C
C     I       = (A+I)*(C+I)-B*D
C     0       = B*(C+I)+(A+I)*D
C
C    -A       = (A+I)*C-B*D
C    -B       = B*C+(A+I)*D
C
C     FROM THE SECOND EQUATION,
C
C     (A+I)*D = -(B+B*C) = -B*(C+I)
C     D       =-INVERSE((A+I))*(B*(C+I))
C
C     SUBSTITUTING THIS EXPRESSION FOR (D) INTO THE FIRST EQUATION.
C
C    -A       = (A+I)*C+B*INVERSE((A+I))*(B*(C+I))
C    -A       =((A+I)+B*INVERSE((A+I))*B)*C+B*INVERSE((A+I))*B
C
C     C       =-INVERSE((A+I)+B*INVERSE((A+I))*B)*(A+B*INVERSE((A+I))*B)
C     D       =-INVERSE((A+I))*B*(C+I)
C
C     THIS METHOD REQUIRES ONLY 2 MATRIX INVERSIONS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
C----- 07/26/09 - EXPANDED DIMENSION FROM 3 X 3 TO 10 X 10
      DIMENSION A(10,10),B(10,10),C(10,10),D(10,10),AI(10,10),
     1 AIB(10,10),AM1(10,10),CM1(10,10),C2(10,10),C3(10,10)
C-----DEFINE A = (A-I)+I
      DO 20 I=1,3
      DO 10 J=1,3
   10 A(I,J)=AM1(I,J)
   20 A(I,I)=A(I,I)+ONE
C-----DEFINE THE INVERSE OF A = 3 X 3
      CALL COFACT3(A,AI)
C-----DEFINE INVERSE(A)*B
      DO 40 I=1,3
      DO 40 J=1,3
      SUM=ZERO
      DO 30 K=1,3
   30 SUM=SUM+AI(I,K)*B(K,J)
   40 AIB(I,J)=SUM
C-----STORE A+B*INVERSE(A)*B IN C AND A'+B*INVERSE(A)*B IN C2.
      DO 60 I=1,3
      DO 60 J=1,3
      SUM=ZERO
      DO 50 K=1,3
   50 SUM=SUM+B(I,K)*AIB(K,J)
      C(I,J)=A(I,J)+SUM
   60 C2(I,J)=AM1(I,J)+SUM
C-----DEFINE THE INVERSE OF A+B*INVERSE(A)*B
      CALL COFACT3(C,C3)
C-----DEFINE C' AND C = C' + I
      DO 90 I=1,3
      DO 80 J=1,3
      SUM=ZERO
      DO 70 K=1,3
   70 SUM=SUM+C3(I,K)*C2(K,J)
      C(I,J)=-SUM
   80 CM1(I,J)=-SUM
   90 C(I,I)=C(I,I)+ONE
C-----DEFINE D.
      DO 110 I=1,3
      DO 110 J=1,3
      SUM=ZERO
      DO 100 K=1,3
  100 SUM=SUM+AIB(I,K)*C(K,J)
  110 D(I,J)=-SUM
      RETURN
      END
      SUBROUTINE COFACT3(A,AI)
C=======================================================================
C
C     USE METHOD OF COFACTORS TO INVERT SYMMETRIC 3 X 3 MATRIX (A).
C     SINCE (A) IS SYMMETRIC, ITS INVERSE (AI) WILL ALSO BE
C     SYMMETRIC.
C
C=======================================================================
      INCLUDE 'implicit.h'
C-----07/26/09 - EXPANDED DIMENSION FROM 3 X 3 TO 10 X 10
      DIMENSION A(10,10),AI(10,10)
C-----DEFINE DETERMINANT (ASSUMING SYMMETRIC).
      DET1=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      DET2=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      DET3=A(2,1)*A(3,2)-A(2,2)*A(3,1)
      DET=A(1,1)*DET1+A(1,2)*DET2+A(1,3)*DET3
C-----DEFINE INVERSE (ASSUMING SYMMETRIC).
      AI(1,1)=DET1/DET
      AI(2,1)=DET2/DET
      AI(3,1)=DET3/DET
      AI(1,2)=AI(2,1)
      AI(2,2)=(A(1,1)*A(3,3)-A(3,1)*A(1,3))/DET
      AI(3,2)=(A(1,2)*A(3,1)-A(1,1)*A(3,2))/DET
      AI(1,3)=AI(3,1)
      AI(2,3)=AI(3,2)
      AI(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))/DET
      RETURN
      END
      SUBROUTINE ORDER(LOW,LHI)
C=======================================================================
C
C     FOR EACH SECTION (I.E., SAME ISOTOPE, ENERGY RANGE, L VALUE) SORT
C     RESONANCES INTO ASCENDING (J, E) ORDER.
C
C     08/29/09 - UPDATED FOR 12 PARAMETERS PER RESONANCE
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 SWITCH
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
C-----SORT IS NOT REQUIRED IF ONLY 1 RESONANCE.
      IF(LOW.EQ.LHI) GO TO 50
      LTOP=LHI
      LOWP1=LOW+1
C
C     J, E SORT.
C
      DO 40 IR=LOWP1,LHI
      JRM1=LOW
      SWITCH=0
      DO 30 JR=LOWP1,LTOP
      IF(RESTAB(2,JRM1).lt.RESTAB(2,JR)) go to 30
      IF(RESTAB(2,JRM1).gt.RESTAB(2,JR)) go to 10
      IF(ENRES(JRM1).le.ENRES(JR)) go to 30
C-----EXCHANGE RESONANCES INTO J, E ORDER.
   10 SWITCH=1
      DRES=ENRES(JRM1)
      ENRES(JRM1)=ENRES(JR)
      ENRES(JR)=DRES
C----- 08/29/09 - UPDATED FOR 12 PARAMETERS PER RESONANCE
      DO 20 I=1,12
      A=RESTAB(I,JRM1)
      RESTAB(I,JRM1)=RESTAB(I,JR)
   20 RESTAB(I,JR)=A
   30 JRM1=JR
C-----STOP SORT IF ALL RESONANCES ARE ALREADY IN J, E ORDER.
      IF(SWITCH.LE.0) GO TO 50
C-----RESONANCE WITH LARGEST J, E IS NOW IN LAST LOCATION (LTOP).
C-----SHORTEN LENGTH OF TABLE TO CONTINUE SORT OF REMAINING RESONANCES.
   40 LTOP=LTOP-1
   50 RETURN
      END
      SUBROUTINE SORTS(X,LX)
C=======================================================================
C
C     SORT AN ARRAY INTO ASCENDING FLOATING POINT ORDER.
C
C     ARGUMENTS
C     ---------
C     X      = ARRAY TO SORT (DIMENSION LX)
C     LX     = NUMBER OF ELEMENTS TO SORT
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 SWITCH
      DIMENSION X(LX)
C-----IF LESS THAN TWO (2) ELEMENTS NO SORT IS REQUIRED.
      IF(LX.LT.2) RETURN
C-----SET INNER LOOP INDICES
      LTOP=LX+1
C-----SET UP OUTER LOOP
      DO 30 IN=2,LX
C-----INITIALIZE EXCHANGE SWITCH OFF.
      SWITCH=0
C-----SET UPPER INDEX TO INNER LOOP
      LTOP=LTOP-1
C-----SET LARGEST ELEMENT INDICATOR TO FIRST ELEMENT
      LBIG=1
C-----SET UP INNER LOOP
      DO 20 J=2,LTOP
C-----COMPARE ELEMENTS
      IF(X(LBIG).GT.X(J)) GO TO 10
C-----ELEMENTS ARE IN NUMERICAL ORDER. RESET INDEX TO LARGER ELEMENT.
      LBIG=J
      GO TO 20
C-----ELEMENTS ARE NOT IN NUMERICAL ORDER. SET INTERCHANGE SWITCH.
   10 SWITCH=1
C-----END OF INNER LOOP
   20 CONTINUE
C-----ARE ALL ELEMENTS ALREADY IN ORDER......
      IF(SWITCH.LE.0) RETURN
C-----NO. MOVE LARGEST ELEMENT TO TOP OF REMAINING TABLE
      DUMMY=X(LBIG)
      X(LBIG)=X(LTOP)
      X(LTOP)=DUMMY
   30 CONTINUE
      RETURN
      END
      SUBROUTINE SORTD(X,LX)
C=======================================================================
C
C     SORT AN ARRAY INTO ASCENDING DOUBLE PRECISION ORDER.
C
C     ARGUMENTS
C     ---------
C     X      = ARRAY TO SORT (DIMENSION LX)
C     LX     = NUMBER OF ELEMENTS TO SORT
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 SWITCH
      DIMENSION X(LX)
C-----IF LESS THAN TWO (2) ELEMENTS NO SORT IS REQUIRED.
      IF(LX.LT.2) RETURN
C-----SET INNER LOOP INDICES
      LTOP=LX+1
C-----SET UP OUTER LOOP
      DO 30 IN=2,LX
C-----INITIALIZE EXCHANGE SWITCH OFF.
      SWITCH=0
C-----SET UPPER INDEX TO INNER LOOP
      LTOP=LTOP-1
C-----SET LARGEST ELEMENT INDICATOR TO FIRST ELEMENT
      LBIG=1
C-----SET UP INNER LOOP
      DO 20 J=2,LTOP
C-----COMPARE ELEMENTS
      IF(X(LBIG).GT.X(J)) GO TO 10
C-----ELEMENTS ARE IN NUMERICAL ORDER. RESET INDEX TO LARGER ELEMENT.
      LBIG=J
      GO TO 20
C-----ELEMENTS ARE NOT IN NUMERICAL ORDER. SET INTERCHANGE SWITCH.
   10 SWITCH=1
C-----END OF INNER LOOP
   20 CONTINUE
C-----ARE ALL ELEMENTS ALREADY IN ORDER......
      IF(SWITCH.LE.0) RETURN
C-----NO. MOVE LARGEST ELEMENT TO TOP OF REMAINING TABLE
      DUMMY=X(LBIG)
      X(LBIG)=X(LTOP)
      X(LTOP)=DUMMY
   30 CONTINUE
      RETURN
      END
      SUBROUTINE SUBINT
C=======================================================================
C
C     11/01/14 - UPDATED FOR ENERGY DEPENDENT STEPS
C     1) 1 UP TO 0.1 KEV
C     2) 2 UP TO   1 KEV
C     3) 4 ABOVE   1 KEV
C     WARNING - TO WORK CORRECTLY (ITEMS=256) MUST BE A MULTIPLE OF 4.
C     ==================================================================
C     DEFINE SUBINTERVALS BETWEEN NODES ACCORDING TO THE WIDTH OF
C     THE TWO ADJACENT RESONANCES. THIS ROUTINE HAS A BUILT-IN TABLE
C     OF MULTIPLES OF HALF-WIDTHS REQUIRED TO FIT A SIMPLE BREIT-WIGNER
C     LINE SPACE TO WITHIN 1.0 PER-CENT ACCURACY OVER THE RANGE 0 TO
C     500 HALF-WIDTHS,
C
C     SIGMA(X)=1.0/(1.0+X*X)
C
C     WHERE X IS THE DISTANCE PEAK IN HALF-WIDTHS
C
C     STARTING FROM THE LOWER ENERGY END OF AN INTERVAL AND USING
C     SUB-INTERVALS BASED UPON THE WIDTH AT THE LOWER ENERGY END, NODES
C     ARE INSERTED AT SUCESSIVELY HIGHER ENERGIES UNTIL A NODE IS
C     CLOSER (IN HALF-WIDTH UNITS) TO THE UPPER ENERGY END OF THE
C     INTERVAL. THEN USING SUB-INTERVALS BASED UPON THE WIDTH AT THE
C     UPPER ENERGY, NODES ARE INSERTED UP TO THE UPPER ENERGY END OF
C     THE INTERVAL.
C
C     WITH THIS ALOGORITHM CLOSELY SPACED RESONANCES WILL HAVE ONLY
C     A FEW SUB-INTERVALS PER INTERVAL (E.G. U-235). WIDELY SPACED
C     RESONANCES WILL HAVE MORE SUB-INTERVALS PER INTERVAL (E.G. U-238).
C     FOR A MIX OF S, P, D ETC. RESONANCES THIS ALOGORITHM GUARANTEES
C     AN ADEQUTE DESCRIPTION OF THE PROFILE OF EVEN EXTREMELY NARROW
C     RESONANCES (WHICH MAY IMMEDIATELY CONVERGENCE TO THE ACCURACY
C     REQUESTED, THUS MINIMIZING ITERATION).
C
C=======================================================================
      INCLUDE 'implicit.h'
c----- 2016/11/24 - Added 18 (48 - 30)
c----- 2016/11/24 - Added 2 more to insure TOTAL is multiple of 4
      PARAMETER (ITIMES = 276)  ! 4 X 69
      COMMON/SUBS/ESUB(1000),ENODP,ENODM,WIDP,WIDM,ISUB,NSUB
C-----10/02/02 - GREATLY INCREASED # OF INTERVALS
c-----added 18 = 3 lines
      DIMENSION TIMES(0:ITIMES)
C-----DEFINE HALF WIDTH SPACINGS OF SUBINTERVALS.
      DATA TIMES/
     1 0.000D+00, 1.000D-03, 2.000D-03, 3.000D-03, 4.000D-03, 5.000D-03,
     2 6.000D-03, 7.000D-03, 8.000D-03, 9.000D-03, 1.000D-02, 2.000D-02,
     3 3.000D-02, 4.000D-02, 5.000D-02, 6.000D-02, 7.000D-02, 8.000D-02,
     4 9.000D-02, 1.000D-01, 1.100D-01, 1.200D-01, 1.300D-01, 1.400D-01,
     4 1.500D-01, 1.600D-01, 1.700D-01, 1.800D-01, 1.900D-01, 2.000D-01,
     4 2.100D-01, 2.200D-01, 2.300D-01, 2.400D-01, 2.500D-01, 2.600D-01,
     5 2.800D-01, 3.000D-01, 3.200D-01, 3.400D-01, 3.600D-01, 3.800D-01,
     6 4.000D-01, 4.200D-01, 4.400D-01, 4.600D-01, 4.800D-01, 5.000D-01,
c-----08/04/23 - added 18 = 3 lines
     7 5.200D-01, 5.400D-01, 5.600D-01, 5.800D-01, 6.000D-01, 6.200D-01,
     7 6.400D-01, 6.600D-01, 6.800D-01, 7.000D-01, 7.200D-01, 7.400D-01,
     7 7.600D-01, 7.800D-01, 8.000D-01, 8.200D-01, 8.400D-01, 8.600D-01,
     7 8.800D-01, 9.000D-01, 9.200D-01, 9.400D-01, 9.600D-01, 9.800D-01,
     7 1.000D+00, 1.020D+00, 1.040D+00, 1.060D+00, 1.080D+00, 1.100D+00,
c-----08/04/23 - added 30 = 5 lines
     9 1.120D+00, 1.140D+00, 1.160D+00, 1.180D+00, 1.200D+00, 1.220D+00,
     9 1.240D+00, 1.260D+00, 1.280D+00, 1.300D+00, 1.320D+00, 1.340D+00,
     9 1.360D+00, 1.380D+00, 1.400D+00, 1.420D+00, 1.440D+00, 1.460D+00,
     9 1.480D+00, 1.500D+00, 1.520D+00, 1.540D+00, 1.560D+00, 1.580D+00,
     9 1.600D+00, 1.620D+00, 1.640D+00, 1.660D+00, 1.680D+00, 1.700D+00,
     9 1.720D+00, 1.740D+00, 1.760D+00, 1.780D+00, 1.800D+00, 1.820D+00,
     9 1.840D+00, 1.860D+00, 1.880D+00, 1.900D+00, 1.920D+00, 1.940D+00,
     9 1.960D+00, 1.980D+00, 2.000D+00, 2.020D+00, 2.060D+00, 2.100D+00,
     9 2.140D+00, 2.180D+00, 2.220D+00, 2.240D+00, 2.260D+00, 2.300D+00,
     3 2.350D+00, 2.400D+00, 2.450D+00, 2.500D+00, 2.600D+00, 2.700D+00,
     4 2.800D+00, 2.900D+00, 3.000D+00, 3.100D+00, 3.200D+00, 3.300D+00,
     5 3.400D+00, 3.600D+00, 3.800D+00, 4.000D+00, 4.200D+00, 4.400D+00,
     6 4.600D+00, 4.800D+00, 5.000D+00, 5.200D+00, 5.400D+00, 5.600D+00,
     7 5.800D+00, 6.000D+00, 6.200D+00, 6.400D+00, 6.500D+00, 6.800D+00,
     8 7.000D+00, 7.500D+00, 8.000D+00, 8.500D+00, 9.000D+00, 9.500D+00,
     9 1.000D+01, 1.050D+01, 1.100D+01, 1.150D+01, 1.200D+01, 1.250D+01,
     A 1.300D+01, 1.350D+01, 1.400D+01, 1.450D+01, 1.500D+01, 1.550D+01,
     1 1.600D+01, 1.700D+01, 1.800D+01, 1.900D+01, 2.000D+01, 2.100D+01,
     2 2.200D+01, 2.300D+01, 2.400D+01, 2.500D+01, 2.600D+01, 2.700D+01,
     3 2.800D+01, 2.900D+01, 3.000D+01, 3.100D+01, 3.200D+01, 3.300D+01,
c----- 2016/11/24 - Added 18 (48 - 30)
     1 3.400D+01, 3.500D+01, 3.600D+01, 3.700D+01, 3.800D+01, 3.900D+01,
     2 4.000D+01, 4.100D+01, 4.200D+01, 4.300D+01, 4.400D+01, 4.500D+01,
     3 4.600D+01, 4.700D+01, 4.800D+01, 4.900D+01, 5.000D+01, 5.100D+01,
     4 5.200D+01, 5.300D+01, 5.400D+01, 5.500D+01, 5.600D+01, 5.700D+01,
     5 5.800D+01, 6.000D+01, 6.200D+01, 6.400D+01, 6.600D+01, 6.800D+02,
     6 7.000D+01, 7.200D+01, 7.400D+01, 7.600D+01, 7.800D+01, 8.000D+02,
     7 8.400D+01, 8.800D+01, 9.200D+01, 9.600D+01, 1.000D+02, 1.040D+02,
     8 1.080D+02, 1.120D+02, 1.180D+02, 1.232D+02, 1.260D+02, 1.300D+02,
     9 1.382D+02, 1.550D+02, 1.600D+02, 1.739D+02, 1.800D+02, 1.951D+02,
     A 2.000D+02, 2.100D+02, 2.189D+02, 2.300D+02, 2.456D+02, 2.500D+02,
     1 2.600D+02, 2.756D+02, 3.092D+02, 3.200D+02, 3.469D+02, 3.600D+02,
c----- 2016/11/24 - Added 2 to be sure total is multiple of 4
     2 3.892D+02, 4.000D+02, 4.200D+02, 4.367D+02, 4.600D+02, 4.800D+02,
     3 5.000D+02, 5.500D+02, 6.000D+02, 7.000D+02, 8.000D+02, 9.000D+02,
     3 1.000D+03/
C
C     LOWER ENERGY END OF INTERVAL.
C
C-----11/01/14 - ENERGY DEPENDENT STEP SIZE
      NDSUB = 1                          ! up to 100 eV
      IF(ENODM.GT.100.0D+0)   NDSUB = 2  ! 100 to 1,000 eV
      IF(ENODM.GT.1000.0D+0)  NDSUB = 4  ! over 1,000 eV
C-----START FIRST SUBINTERVAL AT BEGINNING OF INTERVAL.
      NSUB=1
      ESUB(1)=ENODM
C-----IF LOWER WIDTH IS ZERO DO NOT INSERT ANY MORE SUB-INTERVALS NEAR
C-----LOWER ENERGY END OF INTERVAL.
      IF(WIDM.gt.0.0D+0) go to 20
C-----IF UPPER WIDTH IS ZERO DO NOT INSERT ANY MORE SUB-INTERVALS NEAR
C-----UPPER ENERGY END OF INTERVAL (START ITERATION FROM LOWER AND UPPER
C-----ENERGY LIMITS...UNRESOLVED RESONANCE REGION).
      IF(WIDP.gt.0.0D+0) go to 10
      NSUB=2
      ESUB(2)=ENODP
      RETURN
C-----LOWER WIDTH IS ZERO BUT UPPER LIMIT IS NOT. DEFINE INDEX TO INSERT
C-----AS MANY SUB-INTERVALS AS POSSIBLE.
   10 ISUB=ITIMES
      GO TO 70
C-----LOWER WIDTH IS NOT ZERO. IF UPPER WIDTH IS ZERO, SET MIDPOINT TO
C-----BE EXACTLY EQUAL TO UPPER LIMIT (AVOID ROUND-OFF IN DEFINITION).
   20 IF(WIDP.gt.0.0D+0) go to 30
      EMID=ENODP
      GO TO 40
   30 EMID=(WIDM*ENODP+WIDP*ENODM)/(WIDM+WIDP)
      CALL INCORE9(EMID)
C-----INSERT SUBINTERVALS FROM LOWER ENERGY LIMIT UP TO MIDPOINT
C-----BETWEEN RESONANCES (BASED ON HALF-WIDTH UNITS).
   40 DO 50 ISUB=NDSUB,ITIMES,NDSUB
      ENOW=ENODM+TIMES(ISUB)*WIDM
      CALL INCORE9(ENOW)
      IF(ENOW.ge.EMID) go to 60
      IF(ENOW.le.ESUB(NSUB)) go to 50
      NSUB=NSUB+1
      ESUB(NSUB)=ENOW
   50 CONTINUE
      ISUB=ITIMES
C-----ADD MIDPOINT AS NODE.
   60 NSUB=NSUB+1
      ESUB(NSUB)=EMID
C
C     UPPER ENERGY END OF INTERVAL.
C
C-----IF UPPER WIDTH IS ZERO USE END POINT AS LAST NODE (EITHER ADD
C-----AN ADDITIONAL NODE OR SET LAST NODE ENERGY TO END POINT ENERGY).
      IF(WIDP.gt.0.0D+0) go to 70
      IF(ENODP.le.ESUB(NSUB)) go to 100
      go to 90
C-----INSERT SUBINTERVALS FROM MIDPOINT UP TO NEXT RESONANCE.
   70 IISUB=ISUB
      DO 80 ISUB=IISUB,NDSUB,-NDSUB
      ENOW=ENODP-TIMES(ISUB)*WIDP
      CALL INCORE9(ENOW)
      IF(ENOW.ge.ENODP) go to 90
      IF(ENOW.le.ESUB(NSUB)) go to 80
      NSUB=NSUB+1
      ESUB(NSUB)=ENOW
   80 CONTINUE
C-----ADD ENDPOINT AS NODE.
   90 NSUB=NSUB+1
  100 ESUB(NSUB)=ENODP
      RETURN
      END
      SUBROUTINE NOODLE(ERES,WID,ELOW,EHIGH)
C=======================================================================
C
C     CORE ALLOCATION FOR NODES.
C
C     SAVE ENERGY (ERES) AND WIDTH (WID) AS NODE IN ASCENDING ENERGY
C     ORDER IF IT IS IN THE ENERGY RANGE (ELOW,EHIGH). OTHERWISE IGNOR.
C
C     IF CORE ALLOCATION IS EXCEEDED TERMINATE UNLESS IN EDIT MODE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,ADDNOD
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/ETABS/NODES
      COMMON/MAXIE/NEDSEC,NEDRES,NEDNOD
      COMMON/IWATCH/IMEDIT,MAKEPLUS,MONITR,IMBACK
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DATA ADDNOD/0/
C-----ONLY USE NODES IN RESONANCE REGION.
      IF(ERES.LT.ELOW.OR.ERES.GT.EHIGH) GO TO 90
C-----IF NO OTHER NODES YET SAVE FIRST NODE.
      IF(NODES.LE.0) GO TO 40
C-----COMPARE TO ALL NODES SAVED SO FAR.
      DO 10 I=1,NODES
      IF(ERES.lt.ENODE(I)) go to 20
      IF(ERES.gt.ENODE(I)) go to 10
C-----SAME AS PREVIOUS ENERGY. ADD WIDTHS AND RETURN.
      WIDNOD(I)=WIDNOD(I)+WID
      GO TO 90
   10 CONTINUE
C-----NEW NODE TO ADD TO END OF TABLE.
      GO TO 40
C-----MOVE ALL FOLLOWING NODES BACK ONE LOCATION IN TABLE AND INSERT
C-----NODE IN ENERGY ORDER.
   20 NN=NODES+1
      IF((NN+ADDNOD).GT.NEDNOD) NEDNOD=NN+ADDNOD
      IF(NN.GT.MAXRES) GO TO 60
      DO 30 J=I,NODES
      ENODE(NN)=ENODE(NN-1)
      WIDNOD(NN)=WIDNOD(NN-1)
   30 NN=NN-1
      ENODE(I)=ERES
      WIDNOD(I)=WID
      NODES=NODES+1
      GO TO 90
C-----ADD NEW NODE TO END OF TABLE.
   40 NODES=NODES+1
      IF((NODES+ADDNOD).GT.NEDNOD) NEDNOD=NODES+ADDNOD
      IF(NODES.GT.MAXRES) GO TO 60
   50 ENODE(NODES)=ERES
      WIDNOD(NODES)=WID
      GO TO 90
C-----TOO MANY NODES. PRINT ERROR MESSAGE AND TERMINATE UNLESS IN EDIT
C-----MODE.
   60 IF(IMEDIT.NE.2) GO TO 70
      ADDNOD=ADDNOD+MAXRES
      NODES=1
      GO TO 50
   70 WRITE(OUTP,80) MAXRES
      WRITE(*   ,80) MAXRES
   80 FORMAT(///' ERROR - More than',I5,' Nodes'/
     1          '         Re-Run Program in Edit Mode to Determine',
     2                                          ' Memory Requirements'/
     3          '         Execution Terminated'///)
      CALL ENDERROR
   90 RETURN
      END
      SUBROUTINE FACTS2(L,RHO2,SF2,PF)
C=======================================================================
C
C     CALCULATE 1/2 THE SHIFT FACTOR (SF2) AND PENETRATION FACTOR (PF)
C     FROM L VALUE AND RHO2 (RHO2=RHO*RHO).
C
C     THE SHIFT FACTOR IS USED TO DEFINE THE ENERGY DEPENDENT RESONANCE
C     PEAK ENERGY IN THE FORM
C
C     ERP=ER+0.5*GAMN*(SF(ABS(ER))-SF(E))
C
C     IN ORDER TO AVOID HAVING TO MULTIPLY BY 0.5 EVERYTIME THE SHIFT
C     FACTOR IS USED THIS ROUTINE WILL RETURN 1/2 THE SHIFT FACTOR.
C     THIS IS DONE SIMPLY BY DEFINING ALL CONSTANTS WHICH ARE USED TO
C     DERIVE THE SHIFT FACTOR TO BE 1/2 THERE NORMAL DEFINITION (THUS
C     COMPLETELY AVOIDING THE MULTIPLICATION BY 1/2).
C
C     FOR NEGATIVE ENERGY PENETRATION AND SHIFT FACTORS ARE SET EQUAL
C     TO ZERO (SEE SUMMARY OF CSEWG MEETING, MAY 16-17, 1979
C     RESONANCE SUB-COMMITTEE REPORT). THIS CAN OCCUR FOR THE
C     COMPETITIVE WIDTH AT ENERGIES BELOW THE COMPETITIVE THRESHOLD.
C
C     04/29/04 - ADDED EXTENSION ABOVE L = 5 TO INFINITY
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      DATA DNINE/9.0D+00/
      DATA CA/1.5D+00/
      DATA CB/2.25D+02/
      DATA C3/4.5D+01/
      DATA C4/6.0D+00/
      DATA C5/3.375D+02/
      DATA C6/3.0D+00/
      DATA C7/1.1025D+04/
      DATA C8/1.575D+03/
      DATA C9/1.35D+02/
      DATA C10/1.0D+01/
      DATA C11/2.2050D+04/
      DATA C12/2.3625D+03/
      DATA C13/1.35D+02/
      DATA C14/5.0D+00/
      DATA C15/8.93025D+05/
      DATA C16/9.9225D+04/
      DATA C17/6.3D+03/
      DATA C18/3.15D+02/
      DATA C19/1.5D+01/
      DATA C20/2.2325625D+06/
      DATA C21/1.9845D+05/
      DATA C22/9.45D+03/
      DATA C23/3.15D+02/
      DATA C24/7.5D+00/
      IF(RHO2.GT.ZERO) GO TO 10
      SF2=ZERO
      PF=ZERO
      RETURN
   10 RHO=DSQRT(RHO2)
      IF(L.GT.5) GO TO 80
      GO TO (20,30,40,50,60,70),L+1
C-----S-WAVE (L=0)
   20 SF2=ZERO
      PF=RHO
      RETURN
C-----P-WAVE (L=1)
   30 DEN=ONE+RHO2
      PF=RHO2*RHO/DEN
      SF2=-HALF/DEN
      RETURN
C-----D-WAVE (L=2)
   40 RHO4=RHO2*RHO2
      DEN=THREE*RHO2+RHO4+DNINE
      PF=RHO4*RHO/DEN
      SF2=-(DNINE+CA*RHO2)/DEN
      RETURN
C-----F-WAVE (L=3)
   50 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      DEN=CB+C3*RHO2+C4*RHO4+RHO6
      PF=RHO6*RHO/DEN
      SF2=-(C5+C3*RHO2+C6*RHO4)/DEN
      RETURN
C-----G-WAVE (L=4)
   60 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      RHO8=RHO4*RHO4
      DEN=C7+C8*RHO2+C9*RHO4+C10*RHO6+RHO8
      PF=RHO8*RHO/DEN
      SF2=-(C11+C12*RHO2+C13*RHO4+C14*RHO6)/DEN
      RETURN
C-----(L=5)
   70 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      RHO8=RHO4*RHO4
      RHO10=RHO4*RHO6
      DEN=C15+C16*RHO2+C17*RHO4+C18*RHO6+C19*RHO8+RHO10
      PF=RHO10*RHO/DEN
      SF2=-(C20+C21*RHO2+C22*RHO4+C23*RHO6+C24*RHO8)/DEN
      RETURN
C
C     (L>5)
C
C-----FIRST DEFINE L = 5 = EXACTLY AS ABOVE, BUT HERE WE MUST START
C-----WITH THE SHIFT FACTOR - ABOVE 1/2 THE SHIFT FACTOR IS USED -
C-----AT THE END THIS IS CORRECTED TO RETURN 1/2 THE SHIFT FACTOR.
   80 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      RHO8=RHO4*RHO4
      RHO10=RHO4*RHO6
      DEN=C15+C16*RHO2+C17*RHO4+C18*RHO6+C19*RHO8+RHO10
      PF=RHO10*RHO/DEN
      SF=-2.0D+0*(C20+C21*RHO2+C22*RHO4+C23*RHO6+C24*RHO8)/DEN
C-----THEN RECURSION TO L
      DO LL=6,L
      FLL  = LL
      PAR  = FLL-SF
      DEN  = PAR**2+PF**2
      PF   = RHO2*PF/DEN
      SF   = RHO2*PAR/DEN - FLL
      ENDDO
C-----RETURN 1/2 SHIFT FACTOR
      SF2  = SF/2.0D+0
      RETURN
      END
      SUBROUTINE FACTS3(L,RHO2,PF)
C=======================================================================
C
C     SAME AS FACTS2 WITHOUT THE SHIFT FACTOR
C     REICH-MOORE DOES NOT USE SHIFT FACTOR
C     USE THIS DURING READ AND SIGMA CALCULATIONS
C
C     CALCULATE 1/2 THE SHIFT FACTOR (SF2) AND PENETRATION FACTOR (PF)
C     FROM L VALUE AND RHO2 (RHO2=RHO*RHO).
C
C     THE SHIFT FACTOR IS USED TO DEFINE THE ENERGY DEPENDENT RESONANCE
C     PEAK ENERGY IN THE FORM
C
C     ERP=ER+0.5*GAMN*(SF(ABS(ER))-SF(E))
C
C     IN ORDER TO AVOID HAVING TO MULTIPLY BY 0.5 EVERYTIME THE SHIFT
C     FACTOR IS USED THIS ROUTINE WILL RETURN 1/2 THE SHIFT FACTOR.
C     THIS IS DONE SIMPLY BY DEFINING ALL CONSTANTS WHICH ARE USED TO
C     DERIVE THE SHIFT FACTOR TO BE 1/2 THERE NORMAL DEFINITION (THUS
C     COMPLETELY AVOIDING THE MULTIPLICATION BY 1/2).
C
C     FOR NEGATIVE ENERGY PENETRATION AND SHIFT FACTORS ARE SET EQUAL
C     TO ZERO (SEE SUMMARY OF CSEWG MEETING, MAY 16-17, 1979
C     RESONANCE SUB-COMMITTEE REPORT). THIS CAN OCCUR FOR THE
C     COMPETITIVE WIDTH AT ENERGIES BELOW THE COMPETITIVE THRESHOLD.
C
C     04/29/04 - ADDED EXTENSION ABOVE L = 5 TO INFINITY
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      DATA DNINE/9.0D+00/
      DATA CB/2.25D+02/
      DATA C3/4.5D+01/
      DATA C4/6.0D+00/
      DATA C7/1.1025D+04/
      DATA C8/1.575D+03/
      DATA C9/1.35D+02/
      DATA C10/1.0D+01/
      DATA C15/8.93025D+05/
      DATA C16/9.9225D+04/
      DATA C17/6.3D+03/
      DATA C18/3.15D+02/
      DATA C19/1.5D+01/
      DATA C20/2.2325625D+06/
      DATA C21/1.9845D+05/
      DATA C22/9.45D+03/
      DATA C23/3.15D+02/
      DATA C24/7.5D+00/
      IF(RHO2.GT.ZERO) GO TO 10
      PF=ZERO
      RETURN
   10 RHO=DSQRT(RHO2)
      IF(L.GT.5) GO TO 80
      GO TO (20,30,40,50,60,70),L+1
C-----S-WAVE (L=0)
   20 PF=RHO
      RETURN
C-----P-WAVE (L=1)
   30 DEN=ONE+RHO2
      PF=RHO2*RHO/DEN
      RETURN
C-----D-WAVE (L=2)
   40 RHO4=RHO2*RHO2
      DEN=THREE*RHO2+RHO4+DNINE
      PF=RHO4*RHO/DEN
      RETURN
C-----F-WAVE (L=3)
   50 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      DEN=CB+C3*RHO2+C4*RHO4+RHO6
      PF=RHO6*RHO/DEN
      RETURN
C-----G-WAVE (L=4)
   60 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      RHO8=RHO4*RHO4
      DEN=C7+C8*RHO2+C9*RHO4+C10*RHO6+RHO8
      PF=RHO8*RHO/DEN
      RETURN
C-----(L=5)
   70 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      RHO8=RHO4*RHO4
      RHO10=RHO4*RHO6
      DEN=C15+C16*RHO2+C17*RHO4+C18*RHO6+C19*RHO8+RHO10
      PF=RHO10*RHO/DEN
      RETURN
C
C     (L>5)
C
C-----FIRST DEFINE L = 5 = EXACTLY AS ABOVE, BUT HERE WE MUST START
C-----WITH THE SHIFT FACTOR - ABOVE 1/2 THE SHIFT FACTOR IS USED -
C-----AT THE END THIS IS CORRECTED TO RETURN 1/2 THE SHIFT FACTOR.
   80 RHO4=RHO2*RHO2
      RHO6=RHO4*RHO2
      RHO8=RHO4*RHO4
      RHO10=RHO4*RHO6
      DEN=C15+C16*RHO2+C17*RHO4+C18*RHO6+C19*RHO8+RHO10
      PF=RHO10*RHO/DEN
      SF=-2.0D+0*(C20+C21*RHO2+C22*RHO4+C23*RHO6+C24*RHO8)/DEN
C-----THEN RECURSION TO L
      DO LL=6,L
      FLL  = LL
      PAR  = FLL-SF
      DEN  = PAR**2+PF**2
      PF   = RHO2*PF/DEN
      SF   = RHO2*PAR/DEN - FLL
      ENDDO
      RETURN
      END
      SUBROUTINE FACPHI(L,RHOP,PS)
C=======================================================================
C
C     CALCULATE PHASE SHIFT (PS) FROM L VALUE AND RHOP.
C
C     04/29/04 - ADDED EXTENSION ABOVE L = 5 TO INFINITY
C     02/20/14 - No Lower Value Cutoff (compatible with SAMRML)
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      DATA CA/1.5D+01/
      DATA CB/6.0D+00/
      DATA C3/1.05D+02/
      DATA C4/1.0D+01/
      DATA C5/4.5D+01/
      DATA C6/9.45D+02/
      DATA C7/4.2D+02/
C-----COPIED FROM FACT2S FOR L = 5 PENETRATION AND SHIFT FACTORS
      DATA C15/8.93025D+05/
      DATA C16/9.9225D+04/
      DATA C17/6.3D+03/
      DATA C18/3.15D+02/
      DATA C19/1.5D+01/
      DATA C20/2.2325625D+06/
      DATA C21/1.9845D+05/
      DATA C22/9.45D+03/
      DATA C23/3.15D+02/
      DATA C24/7.5D+00/
      IF(L.GT.5) GO TO 90
      GO TO (10,20,30,40,50,60),L+1
C-----S-WAVE (L=0)
   10 PS=RHOP
      RETURN
C-----P-WAVE (L=1)
   20 PS=RHOP-DATAN(RHOP)
      GO TO 80
C-----D-WAVE (L=2)
   30 RATIO=THREE/(THREE-RHOP*RHOP)
      GO TO 70
C-----F-WAVE (L=3)
   40 RHOP2=RHOP*RHOP
      RATIO=(CA-RHOP2)/(CA-CB*RHOP2)
      GO TO 70
C-----G-WAVE (L=4)
   50 RHOP2=RHOP*RHOP
      RATIO=(C3-C4*RHOP2)/(C3-(C5-RHOP2)*RHOP2)
      GO TO 70
C-----(L=5)
   60 RHOP2=RHOP*RHOP
      RATIO=(C6-(C3-RHOP2)*RHOP2)/(C6-(C7-CA*RHOP2)*RHOP2)
   70 PS=RHOP-DATAN(RATIO*RHOP)
   80 RETURN
C
C     (L>5)
C
C-----FIRST DEFINE L = 5
   90 RHOP2=RHOP*RHOP
C-----PENETRATION AND LEVEL SHIFT FACTORS
      RHOP4=RHOP2*RHOP2
      RHOP6=RHOP4*RHOP2
      RHOP8=RHOP4*RHOP4
      RHOP10=RHOP4*RHOP6
      DEN=C15+C16*RHOP2+C17*RHOP4+C18*RHOP6+C19*RHOP8+RHOP10
      PFX=RHOP10*RHOP/DEN
      SF=-2.0D+0*(C20+C21*RHOP2+C22*RHOP4+C23*RHOP6+C24*RHOP8)/DEN
C-----PHASE SHIFT
      RHOP2=RHOP*RHOP
      RATIO=(C6-(C3-RHOP2)*RHOP2)/(C6-(C7-CA*RHOP2)*RHOP2)
      PS=RHOP-DATAN(RATIO*RHOP)
C-----THEN RECURSION TO L
      DO LL=6,L
      FLL  = LL
      PAR  = FLL-SF
      PS   = PS - DATAN(PFX/PAR)
      DEN  = PAR**2+PFX**2
      PFX   = RHOP2*PFX/DEN
      SF   = RHOP2*PAR/DEN - FLL
      ENDDO
      RETURN
      END
CAK   SUBROUTINE UNFAC(L,RHO2,RHOC,AMUN,VL,PS)
      SUBROUTINE UNFACpre(L,RHO2,RHOC,AMUN,VL,PS)
C=======================================================================
C
C     CALCULATE THE PENETRABILITY FACTOR (VL) AND PHASE SHIFT (PS)
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      DATA DNINE/9.0D+00/
      LL=L+1
      GO TO (10,20,30),LL
C-----S-WAVE
   10 VL=AMUN
      PS=RHOC
      GO TO 40
C-----P-WAVE
   20 VL=AMUN*RHO2/(ONE+RHO2)
      PS=RHOC-DATAN(RHOC)
      GO TO 40
C-----D-WAVE
   30 RHO4=RHO2*RHO2
      VL=AMUN*RHO4/(DNINE+THREE*RHO2+RHO4)
      PS=RHOC-DATAN(THREE*RHOC/(THREE-RHOC*RHOC))
   40 RETURN
      END
      SUBROUTINE GNRL3(GNX,GGX,GFX,GXX,MUN,MUF,MUX,RN,RC,RF,RX)
C-----9/20/10 - ADDED RX to above argument list
C=======================================================================
C
C     CALCULATE UNRESOLVED RESONANCE FLUCTUATION FUNCTION
C     (ORIGINAL CODING FROM AVERAGE4 BY MULKI BHAT)
C     (NEW WEIGHTING SCHEME FROM MC SQUARED)
C
C     THIS ROUTINE HAS BEEN MODIFIED TO CALCULATE ELASTIC, CAPTURE
C     AND FISSION FLUCTUATION FUNCTIONS ALL DURING ONE CALL (AS
C     OPPOSED TO THE ORIGINAL VERSION WHICH CALCULATED EACH REACTION
C     SEPARATELY).
C
C     GNX, GGX, GFX AND GXX ARE THE WIDTHS FOR ELASTIC, CAPTURE,
C     FISSION AND COMPETITION. MUN, MUF AND MUX ARE THE NUMBER OF
C     DEGREES OF FREEDOM FOR ELASTIC, FISSION AND COMPETITION (INFINITE
C     NUMBER OF DEGREES ASSUMED FOR CAPTURE). RN, RC AND RF ARE THE
C     CALCULATED FLUCTUATION INTEGRALS FOR ELASTIC, CAPTURE AND FISSION
C
C     THE NUMBER OF DEGREES OF FREEDOM FOR EACH DISTRIBUTION (ELASTIC,
C     FISSION OR COMPETITION) MAY BE 1 TO 4. IF THE NUMBER OF DEGREES
C     OF FREEDOM FOR ANY DISTRIBUTION IS LESS THAN 1 OR MORE THAN 4
C     IT WILL BE TREATED AS AN INFINITE NUMBER OF DEGREES OF FREEDOM
C     (WHICH INFERS THAT THE WIDTHS ARE NOT DISTRIBUTED, BUT ARE RATHER
C     ALL EQUAL TO THE AVERAGE VALUE). THIS LAST CASE IS SIMULATED BY
C     DEFINING AN ADDITIONAL 10 POINT QUADRATURE IN WHICH THE WEIGHT
C     FOR ONE POINT IS 1.0 AND THE WEIGHT FOR ALL OTHER POINTS IS ZERO.
C     FOR THE ONE POINT OF WEIGHT 1.0 THE AVERAGE WIDTH WILL BE USED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION XX(5,10),WW(5,10)
      DATA XX/
     A 3.0013465D-03,1.3219203D-02,1.0004488D-03,1.3219203D-02,1.0D+0,
     1 7.8592886D-02,7.2349624D-02,2.6197629D-02,7.2349624D-02,0.0D+0,
     2 4.3282415D-01,1.9089473D-01,1.4427472D-01,1.9089473D-01,0.0D+0,
     3 1.3345267D+00,3.9528842D-01,4.4484223D-01,3.9528842D-01,0.0D+0,
     4 3.0481846D+00,7.4083443D-01,1.0160615D+00,7.4083443D-01,0.0D+0,
     5 5.8263198D+00,1.3498293D+00,1.9421066D+00,1.3498293D+00,0.0D+0,
     6 9.9452656D+00,2.5297983D+00,3.3150885D+00,2.5297983D+00,0.0D+0,
     7 1.5782128D+01,5.2384894D+00,5.2607092D+00,5.2384894D+00,0.0D+0,
     8 2.3996824D+01,1.3821772D+01,7.9989414D+00,1.3821772D+01,0.0D+0,
     9 3.6216208D+01,7.5647525D+01,1.2072069D+01,7.5647525D+01,0.0D+0/
      DATA WW/
     A 1.1120413D-01,3.3773418D-02,3.3376214D-04,1.7623788D-03,1.0D+0,
     1 2.3546798D-01,7.9932171D-02,1.8506108D-02,2.1517749D-02,0.0D+0,
     2 2.8440987D-01,1.2835937D-01,1.2309946D-01,8.0979849D-02,0.0D+0,
     3 2.2419127D-01,1.7652616D-01,2.9918923D-01,1.8797998D-01,0.0D+0,
     4 1.0967668D-01,2.1347043D-01,3.3431475D-01,3.0156335D-01,0.0D+0,
     5 3.0493789D-02,2.1154965D-01,1.7766657D-01,2.9616091D-01,0.0D+0,
     6 4.2930874D-03,1.3365186D-01,4.2695894D-02,1.0775649D-01,0.0D+0,
     7 2.5827047D-04,2.2630659D-02,4.0760575D-03,2.5171914D-03,0.0D+0,
     8 4.9031965D-06,1.6313638D-05,1.1766115D-04,8.9630388D-10,0.0D+0,
     9 1.4079206D-08,0.0000000D+00,5.0989546D-07,0.0000000D+00,0.0D+0/
C-----INITIALIZE FLUCTUATION INTEGRALS FOR ELASTIC, CAPTURE AND FISSION.
      RN=0.0
      RC=0.0
      RF=0.0
      RX=0.0  ! OPTIONAL COMPETITION
C-----INTEGRALS ARE ZERO IF NEUTRON AND CAPTURES WIDTHS ARE NOT
C-----POSITIVE.
      IF(GNX.LE.0.0.OR.GGX.LE.0.0) RETURN
C-----INSURE NUMBER OF DEGREES OF FREEDOM FOR ELASTIC WIDTHS IN O.K.
      IF(MUN.LT.1.OR.MUN.GT.4) MUN=5
C-----IS THERE FISSION.
      IF(GFX.lt.0.0D+0) go to 120
      IF(GFX.gt.0.0D+0) go to 50
C
C     NOT FISSILE.
C
C-----NO FISSION. IS THERE A COMPETITIVE WIDTH.
      IF(GXX.GT.0.0) GO TO 20
C-----NO COMPETITIVE WIDTH. 1-D QUADRATURE.
      DO 10 J=1,10
      XJ=XX(MUN,J)
      FACTOR=WW(MUN,J)*XJ/(GNX*XJ+GGX)
      RN=RN+XJ*FACTOR
   10 RC=RC+FACTOR
      RETURN
C-----COMPETITIVE WIDTH. 2-D QUADRATURE. INSURE NUMBER OF DEGREES OF
C-----FREEDOM FOR COMPETITIVE WIDTH IS O.K.
   20 IF(MUX.LT.1.OR.MUX.GT.4) MUX=5
      DO 40 J=1,10
      XJ=XX(MUN,J)
      WJXJ=WW(MUN,J)*XJ
      EFFJ=GNX*XJ+GGX
      DO 30 K=1,10
      XK=XX(MUX,K)
      FACTOR=WW(MUX,K)*WJXJ/(EFFJ+GXX*XK)
      RN=RN+XJ*FACTOR
      RX=RX+XK*FACTOR  ! OPTONAL COMPETITION
   30 RC=RC+FACTOR
   40 CONTINUE
      RETURN
C
C     FISSILE.
C
C-----FISSION PRESENT. INSURE NUMBER OF DEGREES OF FREEDOM FOR FISSION
C-----IS O.K.
   50 IF(MUF.LT.1.OR.MUF.GT.4) MUF=5
C-----IS THERE A COMPETITIVE WIDTH.
      IF(GXX.GT.0.0) GO TO 80
C-----NO COMPETITIVE WIDTH. 2-D QUADRATURE.
      DO 70 J=1,10
      XJ=XX(MUN,J)
      WJXJ=WW(MUN,J)*XJ
      EFFJ=GNX*XJ+GGX
      DO 60 K=1,10
      XK=XX(MUF,K)
      FACTOR=WW(MUF,K)*WJXJ/(EFFJ+GFX*XK)
      RN=RN+XJ*FACTOR
      RC=RC+FACTOR
   60 RF=RF+XK*FACTOR
   70 CONTINUE
      RETURN
C-----COMPETITIVE WIDTH. 3-D QUADRATURE. INSURE NUMBER OF DEGREES OF
C-----FREEDOM FOR COMPETITIVE WIDTH IS O.K.
   80 IF(MUX.LT.1.OR.MUX.GT.4) MUX=5
      DO 110 J=1,10
      XJ=XX(MUN,J)
      WJXJ=WW(MUN,J)*XJ
      EFFJ=GNX*XJ+GGX
      DO 100 K=1,10
      XK=XX(MUF,K)
      WKWJXJ=WW(MUF,K)*WJXJ
      EFFJK=EFFJ+GFX*XK
      DO 90 I=1,10
      XL=XX(MUX,I)
      FACTOR=WW(MUX,I)*WKWJXJ/(EFFJK+GXX*XL)
      RN=RN+XJ*FACTOR
      RX=RX+XL*FACTOR  ! OPTIONAL COMPETITION
      RC=RC+FACTOR
   90 RF=RF+XK*FACTOR
  100 CONTINUE
  110 CONTINUE
  120 RETURN
      END
      SUBROUTINE TERPUP(E,DUMSET,L,INTX)
C=======================================================================
C
C     INTERPOLATE UNRESOLVED RESONANCE PARAMETERS. RESONANCE SPACING,
C     COMPETITIVE, NEUTRON, CAPTURE AND FISSION WIDTHS ARE ALL
C     INTERPOLATED AT THE SAME TIME AND RETURNED IN LOCATIONS DUMSET(2)
C     THROUGH DUMSET(6), RESPECTIVELY.
C
C     THIS ROUTINE HAS BEEN RECODED IN ORDER TO AVOID ROUND-OFF
C     PROBLEMS ON SHORT WORD LENGTH COMPUTERS, E.G. IBM-360, 370, ETC.
C     THIS ROUTINE IS NOW SLIGHTLY LESS EFFICIENT THAN IN ITS PREVIOUS
C     FORM. HOWEVER ALL INTERPOLATION IS NOW DEFINED AS A WEIGHTED SUM
C     OF TERMS, AS OPPOSED TO THE PREVIOUS FORM WHICH USED DIFFERENCES.
C
C     THIS ROUTINE IS ONLY USED TO INTERPOLATE PARAMETERS IN THE
C     UNRESOLVED RESONANCE REGION. IN ORDER TO INSURE THAT IN ALL CASES
C     INTERPOLATION CAN BE PERFORMED, IF THE PARAMETER IS NOT POSITIVE
C     AND LOG INTERPOLATION OF THE CROSS SECTION IS SPECIFIED (I =4 OR
C     5), THE INTERPOLATION OF THE CROSS SECTION WILL BE SWITCHED TO
C     LINEAR (I=2 OR 3) AND A WARNING MESSAGE WILL BE PRINTED.
C
C     AN ILLEGAL INTERPOLATION CODE OR A NON-POSITIVE ENERGY WHERE LOG
C     ENERGY INTERPOLATION IS REQUIRED INDICATES EITHER AN ERROR IN THE
C     DATA AS IT APPEARS IN THE ENDF/B FORMAT, OR AN ERROR IN THIS
C     PROGRAM. THEREFORE ERRORS OF THIS TYPE WILL CAUSE THE PROGRAM TO
C     PRINT A WARNING MESSAGE AND TERMINATE EXECUTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/CONJOB/PI,PI2,ZERO,HALF,ONE,TWO,THREE,FOUR,EIGHT
      COMMON/FIELDC/FIELD(11,12)
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION DUMSET(6)
C-----CHECK INTERPOLATION CODE.
      IF(INTX.LT.1.OR.INTX.GT.5) GO TO 170
C-----DEFINE INDEX TO LOWER ENERGY LIMIT.
      LM1=L-1
C-----DEFINE ENERGIES AT THE 2 ENDS OF THE INTERVAL.
      E1=ENRES(LM1)
      E2=ENRES(L)
C-----CHECK FOR ZERO LENGTH INTERVAL.
      DE12=E2-E1
      IF(DE12.EQ.0.0) GO TO 210
C-----SELECT INTERPOLATION LAW.
      GO TO (10,110,140,30,70),INTX
C
C     HISTOGRAM. CONSTANT EQUAL TO VALUE AT LOWER ENERGY LIMIT.
C
   10 DO 20 I=2,6
   20 DUMSET(I)=RESTAB(I,LM1)
      RETURN
C
C     LINEAR X AND LOG Y.
C
   30 WT2=(E-E1)/DE12
      WT1=ONE-WT2
      DO 60 I=2,6
C-----IF Y IS NOT POSITIVE OR CONSTANT SET FLAG AND USE LINEAR Y
C-----INTERPOLATION.
      IF(RESTAB(I,LM1).GT.0.0.AND.RESTAB(I,L).GT.0.0) GO TO 40
      IF(RESTAB(I,LM1).EQ.RESTAB(I,L)) GO TO 50
C-----LOG INTERPOLATION OF PARAMETERS IS NOT POSSIBLE. SWITCH TO LINEAR.
      GO TO 120
   40 DUMSET(I)=DEXP(WT1*DLOG(RESTAB(I,LM1))+WT2*DLOG(RESTAB(I,L)))
      GO TO 60
   50 DUMSET(I)=RESTAB(I,LM1)
   60 CONTINUE
      RETURN
C
C     LOG X AND LOG Y.
C
   70 IF(E1.LE.0.0.OR.E2.LE.0.0.OR.E.LE.0.0) GO TO 190
      WT2=DLOG(E/E1)/DLOG(E2/E1)
      WT1=ONE-WT2
      DO 100 I=2,6
C-----IF Y IS NOT POSITIVE OR CONSTANT SET FLAG AND USE LINEAR Y
C-----INTERPOLATIOM.
      IF(RESTAB(I,LM1).GT.0.0.AND.RESTAB(I,L).GT.0.0) GO TO 80
      IF(RESTAB(I,LM1).EQ.RESTAB(I,L)) GO TO 90
C-----LOG INTERPOLATION OF PARAMETERS IS NOT POSSIBLE. SWITCH TO LINEAR.
      GO TO 150
   80 DUMSET(I)=DEXP(WT1*DLOG(RESTAB(I,LM1))+WT2*DLOG(RESTAB(I,L)))
      GO TO 100
   90 DUMSET(I)=RESTAB(I,LM1)
  100 CONTINUE
      RETURN
C
C     LINEAR X AND LINEAR Y.
C
  110 WT2=(E-E1)/DE12
      WT1=ONE-WT2
  120 DO 130 I=2,6
  130 DUMSET(I)=WT1*RESTAB(I,LM1)+WT2*RESTAB(I,L)
      RETURN
C
C     LOG X AND LINEAR Y.
C
C-----INSURE ALL X VALUES ARE POSITIVE FOR LOG.
  140 IF(E1.LE.0.0.OR.E2.LE.0.0.OR.E.LE.0.0) GO TO 190
      WT2=DLOG(E/E1)/DLOG(E2/E1)
      WT1=ONE-WT2
  150 DO 160 I=2,6
  160 DUMSET(I)=WT1*RESTAB(I,LM1)+WT2*RESTAB(I,L)
      RETURN
C-----ILLEGAL INTERPOLATE CODE.
  170 WRITE(OUTP,180) INTX
      WRITE(*   ,180) INTX
  180 FORMAT(///' ERROR - Interpolating Unresolved Parameters.'/
     1          '         Illegal Interpolation Code =',I5,
     2                                    ' (MUST be 1 to 5).'/
     3          '         Correct Evaluated Data and Re-Run Program.'/
     4          '         Execution Terminated.'///)
      CALL ENDERROR
C-----ILLEGAL LOG ENERGY INTERPOLATION WITH NEGATIVE VALUES.
  190 CALL OUT9(E1,FIELD(1,1))
      CALL OUT9(E2,FIELD(1,2))
      CALL OUT9(E ,FIELD(1,3))
      WRITE(OUTP,200) ((FIELD(M,J),M=1,11),J=1,3)
      WRITE(*   ,200) ((FIELD(M,J),M=1,11),J=1,3)
  200 FORMAT(///' ERROR - Interpolating Unresolved Parameters.'/
     2          '         Illegal Log Energy Interpolation Using',
     2                                        ' Negative Energy.'/
     3          '         Interpolation Code=',I5,' (Cannot be 3, 5).'/
     4          '         E1,E2,E=',3(11A1,1X)/
     5          '         Correct Evaluated Data and Re-Run Program.'/
     6          '         Execution Terminated.'///)
      CALL ENDERROR
C-----ZERO LENGTH ENERGY INTERVAL.
  210 CALL OUT9(E1,FIELD(1,1))
      CALL OUT9(E2,FIELD(1,2))
      WRITE(OUTP,220) ((FIELD(M,J),M=1,11),J=1,2)
      WRITE(*   ,220) ((FIELD(M,J),M=1,11),J=1,2)
  220 FORMAT(///' ERROR - Interpolating Unresolved Parameters.'/
     1          '         Illegal Interpolation Over Zero',
     2                                 ' Length Interval.'/
     3          '         E1,E2=',11A1,1X,11A1/
     4          '         Correct Evaluated Data and Re-Run Program.'/
     5          '         Execution Terminated.'///)
      CALL ENDERROR
      RETURN
      END
      SUBROUTINE LISPAR6(IXLOW,IX)
C=======================================================================
C
C     SAME AS LISPAR - BUT 6 PARAMETERS PER RESONANCE
C
C     READ AND LIST ALL RESONANCE PARAMETERS
C
C     PARAMETERS ARE READ IN DOUBLE PRECISION.
C     THE ENERGY AND WIDTH ARE SAVED IN DOUBLE PRECISION TO DEFINE THE
C     ENERGY GRID.
C
C     IXLOW   = STARTING INDEX TO RESONANCE PARAMETER TABLE
C     IX      = NUMBER OF WORDS TO READ
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION XD1(12)
      DATA ZEROD/0.0D+00/
C
C     READ PARAMETERS FOR EACH RESONANCE.
C
      KXLOW=IXLOW-1
      DO 100 I1=1,IX,6
C-----READ IN DOUBLE PRECISION.
      CALL LISTIO9(XD1,6)
C-----COPY TO PARAMETER ARRAY.
      KXLOW=KXLOW+1
      ENRES(KXLOW)=XD1(1)
      DO 10 I=1,6
   10 RESTAB(I,KXLOW)=XD1(I)
C
C     SAVE RESONANCE ENERGY AS A NODE.
C
      IF(LRU.EQ.2) GO TO 90
C-----RESOLVED.
      GO TO (20,20,30,40,50,60,70),LRF
C-----BREIT-WIGNER.
   20 WIDTOT=XD1(3)
      GO TO 80
C-----REICH-MOORE
   30 WIDTOT=XD1(3)+XD1(4)+DABS(XD1(5))+DABS(XD1(6))
      GO TO 80
C-----ADLER-ADLER
   40 WIDTOT=XD1(2)
      GO TO 80
C-----GENERAL-R-MATRIX (NOT IMPLEMENTED).
   50 WIDTOT=ZEROD
      GO TO 80
C-----HYBRID-R-MATRIX
   60 WIDTOT=XD1(2)+XD1(3)+XD1(4)+XD1(5)+XD1(6)+XD1(7)+XD1(8)
      GO TO 80
C-----NEW (2003) REICH-MOORE
   70 WIDTOT = DABS(XD1(2))+DABS(XD1(3))+
     1         DABS(XD1(4))+DABS(XD1(5))
C-----DEFINE RESOLVED REGION NODE.
   80 CALL NOODLE(XD1(1),WIDTOT,ELX,EHX)
      GO TO 100
C
C     DEFINE UNRESOLVED REGION NODE.
C
   90 CALL NOODLE(XD1(1),ZEROD,ELX,EHX)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE LISPAR12(IXLOW,IX)
C=======================================================================
C
C     SAME AS LISPAR - BUT 12 PARAMETERS PER RESONANCE
C
C     READ AND LIST ALL RESONANCE PARAMETERS
C
C     PARAMETERS ARE READ IN DOUBLE PRECISION.
C     THE ENERGY AND WIDTH ARE SAVED IN DOUBLE PRECISION TO DEFINE THE
C     ENERGY GRID.
C
C     IXLOW   = STARTING INDEX TO RESONANCE PARAMETER TABLE
C     IX      = NUMBER OF WORDS TO READ
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
C-----2/14/04 - ADDED INCLUDE
      INCLUDE 'recent.h'
      DIMENSION XD1(12)
      DATA ZEROD/0.0D+00/
C
C     READ PARAMETERS FOR EACH RESONANCE.
C
      KXLOW=IXLOW-1
      DO 100 I1=1,IX,12
C-----READ IN DOUBLE PRECISION.
      CALL LISTIO9(XD1,12)
C-----COPY TO PARAMETER ARRAY.
      KXLOW=KXLOW+1
      ENRES(KXLOW)=XD1(1)
      DO 10 I=1,12
   10 RESTAB(I,KXLOW)=XD1(I)
C
C     SAVE RESONANCE ENERGY AS A NODE.
C
      IF(LRU.EQ.2) GO TO 90
C-----RESOLVED.
      GO TO (20,20,30,40,50,60,70),LRF
C-----BREIT-WIGNER.
   20 WIDTOT=XD1(3)
      GO TO 80
C-----REICH-MOORE
   30 WIDTOT=XD1(3)+XD1(4)+DABS(XD1(5))+DABS(XD1(6))
      GO TO 80
C-----ADLER-ADLER
   40 WIDTOT=XD1(2)
      GO TO 80
C-----GENERAL-R-MATRIX (NOT IMPLEMENTED).
   50 WIDTOT=ZEROD
      GO TO 80
C-----HYBRID-R-MATRIX
   60 WIDTOT=XD1(2)+XD1(3)+XD1(4)+XD1(5)+XD1(6)+XD1(7)+XD1(8)
      GO TO 80
C-----NEW (2003) REICH-MOORE
   70 WIDTOT = DABS(XD1(2))+DABS(XD1(3))+
     1         DABS(XD1(4))+DABS(XD1(5))
C-----DEFINE RESOLVED REGION NODE.
   80 CALL NOODLE(XD1(1),WIDTOT,ELX,EHX)
      GO TO 100
C
C     DEFINE UNRESOLVED REGION NODE.
C
   90 CALL NOODLE(XD1(1),ZEROD,ELX,EHX)
  100 CONTINUE
      RETURN
      END
CAK   SUBROUTINE FILEIO
      SUBROUTINE FILEIOr
C=======================================================================
C
C     DEFINE ALL I/O UNITS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*72 NAMEIN,NAMEOUT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/IOSTATUS/ISTAT1,ISTAT2
      COMMON/UNITS/ISCR2,ISCR23
      COMMON/NAMEX/NAMEIN,NAMEOUT
C-----DEFINE ALL I/O UNITS.
      INP=2
      OUTP=3
      ITAPE=10
      OTAPE=11
      ISCR2=12
      ISCR23=14
C-----DEFINE ALL FILE NAMES.
      OPEN(OUTP,FILE='RECENT.LST',STATUS='UNKNOWN')
      CALL SCRATCH1(ISCR2 ,'RECENT.001  ')
      CALL SCRATCH1(ISCR23,'RECENT.002  ')
      OPEN(INP,FILE='RECENT.INP',STATUS='OLD',ERR=10)
      ISTAT1 = 0
      RETURN
   10 ISTAT1 = 1
      RETURN
CAK   ENTRY FILIO2
      ENTRY FILIO2r
C=======================================================================
C
C     DEFINE ENDF/B DATA I/O UNITS AND OPTIONALLY DEFINE FILE NAMES.
C
C=======================================================================
      OPEN(OTAPE,FILE=NAMEOUT,STATUS='UNKNOWN')
      OPEN(ITAPE,FILE=NAMEIN,STATUS='OLD',ERR=20)
      ISTAT2 = 0
      RETURN
   20 ISTAT2 = 1
      RETURN
      END
      SUBROUTINE COULOMB(RHO, LNOW, ETA, PENCOUL)
C=======================================================================
C
C     CALCULATE PENETRATION FACTOR INCLUDING COULOMB
C     ----------------------------------------------
C     04/09/10 - LIBERATED (I.E., STOLEN) FROM NANCY LARSON'S SAMMY CODE
C                MODIFIED TO ONLY CALCULATE PENETRATION FACTOR -
C                NO SHIFT FACTOR OR PHASE SHIFT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      PARAMETER (MAXL=100)
      COMMON/COOLCOM/F(MAXL),FPR(MAXL),G(MAXL),GPR(MAXL),
     1 SIGCOOL(MAXL)
C
C     INITIALIZE AND TEST
C
      IFAIL   = 0
      PENCOUL = 0.0D+0
      IF (LNOW.LT.0) RETURN
      LLMAX = LNOW + 2
      IF (LLMAX.GT.100) RETURN
C
C     CALCULATE PENETRATION FACTOR
C
      IF (RHO.LT.0.02D0) THEN
C-----SMALL RHO
      CALL COOL1 (ETA, RHO, LNOW, LLMAX, PENCOUL)
      ELSE
C-----LARGE RHO
      CALL COOL2 (RHO, ETA, LNOW, LLMAX, PENCOUL, IFAIL, IEXP)
C-----IF LARGE RHO BRANCH FAILS, TRY SMALL RHO
      IF (IFAIL.NE.0) THEN
      ELSE IF (IFAIL.GT.0) THEN
      CALL COOL1 (ETA, RHO, LNOW, LLMAX, PENCOUL)
      ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE COOL1 (EETA,RRHO,LNOW,LLMAX,PENCOUL)
C=======================================================================
C
C     GERRY HALE'S ROUTINE FOR SMALL RHO
C
C=======================================================================
      INCLUDE 'implicit.h'
      PARAMETER (MAXL=100)
      COMMON/COOLCOM/F(MAXL),FPR(MAXL),G(MAXL),GPR(MAXL),
     1 SIGCOOL(MAXL)
      DATA FIVE /5.0D0/
      DATA TEN /10.0D0/
C
      ETA = EETA
      RHO = RRHO
      LMAX = LLMAX
C
C
      IF (ETA.GT.TEN*RHO .AND. ETA.GT.FIVE) THEN
C ------------------------------------------------
C ***    HERE FOR ETA>>RHO
      CALL BIGETA (ETA, RHO, LNOW, LLMAX, PENCOUL)
C ------------------------------------------------
C
      ELSE
C
C ------------------------------------------------
C ***    HERE FOR ETA AND RHO NOT SO VERY DIFFERENT
C
      IF (ETA.GE.FIVE) THEN
C *** GENERATE U, UPR, FIRST FOR RHOI.NE.RHO
C *** I.E., USE ASYMPTOTIC FORMULA FOR LARGE ETA TO GIVE U=G0(2*ETA,ETA)
      CALL ASYMP1 (ETA, RHOI, U, UPR)
C *** NOW GENERATE TAYLOR SERIES EXPANSION OF G0(RHO) AROUND G0(RHOI),
C *** WITH RHOI REDEFINED IF NECESSARY TO OBTAIN CONVERGENCE
      CALL TAYLOR (ETA, RHO, U, UPR, RHOI)
C
      ELSE
C ***       GENERATE SIGCOOL TO USE IN ASYMPTOTIC EXPANSION FOR RHOI
      CALL XSIGLL (ETA, SIGMA0, LMAX)
C ***       USE ASYMPTOTIC FORMULA TO GIVE G0(RHOI,ETA)
      CALL ASYMP2 (ETA, RHO, U, UPR, RHOI, SIGMA0)
      IF (DABS(U).LE.1.0E25) THEN
      CALL TAYLOR (ETA, RHO, U, UPR, RHOI)
      ENDIF
      ENDIF
C
C ***    NOW FIND PENCOUL, ETC...
      G0 = U
      G0PR = UPR
      IF (DABS(G0).GT.1.D+25) THEN
      LLMAX  = 0
      G  (1) = G0
      GPR(1) = G0PR
      F  (1) = 0.0
      FPR(1) = 0.0
      ELSE
      CALL GETFG (ETA, RHO, LLMAX, LNOW, G0, G0PR)
      ENDIF
      LPLUS1  = LNOW + 1
      ASQ     = F(LPLUS1)**2 + G(LPLUS1)**2
      PENCOUL = RHO/ASQ
C
C ------------------------------------------------
C
      ENDIF
      RETURN
      END
      SUBROUTINE COOL2 (XX, ETA1, LNOW, LLMAX, PENCOUL, IFAIL, IEXP)
C=======================================================================
C
C     NANCY LARSON'S ROUTINE FOR LARGE RHO
C
C=======================================================================
      INCLUDE 'implicit.h'
      PARAMETER (MAXL=100)
      COMMON/COOLCOM/F(MAXL),FPR(MAXL),G(MAXL),GPR(MAXL),
     1 SIGCOOL(MAXL)
      DATA ZERO/0.0D0/
      DATA ONE /1.0D0/
      DATA TWO /2.0D0/
      DATA TEN/10.0D0/
      DATA TEN2/100.0D0/
      DATA HALF/0.5D0/
      DATA TM30/1.0D-30/
      DATA ABORT/2.0D4/
C      DATA RT2EPI /0.79788 45608 02865 35587 98921 19868 76373 D0/
C *** THIS CONSTANT IS  DSQRT(TWO/PI)&  USE Q0 FOR IBM REAL*16& D0 FOR
C ***  REAL*8 + CDC DOUBLE P&  E0 FOR CDC SINGLE P; AND TRUNCATE VALUE.
C
      ACCUR = 1.0D-16
C ***            CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
C
C
      IFAIL = 0
      IEXP  = 1
      ETA   = ETA1
      GJWKB = ZERO
      PACCQ = ONE
      ACC   = ACCUR
      ACC4  = ACC*TEN2*TEN2*0.1D0
C                          *0.1D0 ADDED BY ROS 8 APR 98
      ACCH  = DSQRT(ACC)
      IF (XX.LE.ACCH) GO TO 50
      X   = XX
      XLM = ZERO
C
      E2MM1 = ETA*ETA + XLM*XLM + XLM
C           =   ETA^2 + LL(LL+1)
C
      if(X*(X-TWO*ETA) .LT. XLM*XLM + XLM) then
      IXLTURN = 1
      else
      IXLTURN = 0
      endif
C
      XLL   = DFLOAT(LLMAX)
C ***         XLL IS MAX LAMBDA VALUE, OR 0.5 SMALLER FOR J,Y BESSELS
C ***         DETERMINE STARTING ARRAY ELEMENT (1) FROM 0
      L1 = LLMAX + 1
C
C
C
C *** EVALUATE CF1 = FX  = FPRIME(XL,ETA,X)/FX(XL,ETA,X)
C
      XI  = ONE/X
      FCL = ONE
      PK  = XLL + ONE
      PX  = PK  + ABORT
   10 CONTINUE
      EK  = ETA / PK
      FX  = (EK+PK*XI)*FCL + (FCL-ONE)*XI
      PK1 =  PK + ONE
      IF (DABS(ETA*X+PK*PK1).LE.ACC) THEN
C ***    TEST ENSURES B1.NE.ZERO FOR NEGATIVE ETA; FIXUP IS EXACT.
      FCL  = (ONE+EK*EK)/(ONE+(ETA/PK1)**2)
      PK   = TWO + PK
      GO TO 10
      ENDIF
      D  = ONE / ( (PK+PK1)*(XI+EK/PK1) )
      DF = -FCL*(ONE+EK*EK)*D
      IF (FCL.NE.ONE) FCL = -ONE
      IF (D.LT.ZERO ) FCL = -FCL
      FX  = FX + DF
C
C
C
C ***   BEGIN CF1 LOOP ON PK = K = LAMBDA + 1
C
      P = ONE
   20 CONTINUE
      PK  = PK1
      PK1 = PK1 + ONE
      EK  = ETA / PK
      TK  = (PK+PK1)*(XI+EK/PK1)
      D   =  TK - D*(ONE+EK*EK)
      IF (DABS(D).LE.ACCH) THEN
      WRITE (6,30) D, DF, ACCH, PK, EK, ETA, X
   30 FORMAT(/' CF1 ACCURACY LOSS: D,DF,ACCH,K,ETA/K,ETA,X = ',
     *         1P7E9.2/)
      P = P + ONE
      IF (P.GT.TWO) GO TO 60
      ENDIF
      D = ONE/D
      IF (D.LT.ZERO) FCL = -FCL
      DF = DF*(D*TK-ONE)
      FX = FX + DF
      IF (PK.GT.PX) GO TO 60
      IF (DABS(DF).GE.DABS(FX)*ACC) GO TO 20
C
C DOWNWARD RECURRENCE TO LAMBDA=XLM. ARRAY G, IF PRESENT,STORES RL
C
      IF (LLMAX.GT.0) THEN
      FCL = FCL*TM30
      FPL = FCL*FX
      FPR(L1) = FPL
      F (L1) = FCL
      XL  = XLL
      RL  = ONE
      EL  = ZERO
      DO LP=1,LLMAX
      EL     = ETA/XL
      RL     = DSQRT(ONE+EL*EL)
      SL     =  EL + XL*XI
      L      =  L1 - LP
      FCL1   = (FCL *SL+FPL)/RL
      FPL    =  FCL1*SL-FCL *RL
      FCL    =  FCL1
      F(L)  =  FCL
      FPR(L) = FPL
      G(L+1)= RL
      XL     = XL - ONE
      ENDDO
      IF (FCL.EQ.ZERO) FCL = ACC
      FX  = FPL/FCL
C ***    NOW WE HAVE REACHED LAMBDA = 0
C ***    EVALUATE CF2 = P + I.Q  AGAIN USING STEED'S ALGORITHM
C ***    SEE TEXT FOR COMPACT COMPLEX CODE FOR SP CDC OR NON-ANSI IBM
      ENDIF
C
C
      if (IXLTURN.ne.0) then
c-----2014/11/02 - corrected third argument to floating point
c-----             it was MAX(XLM,ZERO) here, but not in Subroutine
      XLMMAX = DMAX1(XLM,ZERO)
      CALL JWKB (X, ETA, XLMMAX, FJWKB, GJWKB, IEXP)
      endif
C
C
      IF (IEXP.GT.1 .OR. GJWKB.GT.ONE/(ACCH*TEN2)) THEN
C
C ***ARRIVE HERE IF [G(XLM).GT.10**6] OR [(IEXP.GT.250+XLTURN)=.TRUE.]
C ***ARRIVE HERE IF [G(XLM).GT.10**6] OR [(IEXP.GT.70 &XLTURN)=.TRUE.]
C ***IN OTHER WORDS, WHERE VALUES ARE EXTREME
      W   = FJWKB
      GAM = GJWKB*W
      P   = FX
      Q   = ONE
C
      ELSE
C
      IXLTURN = 0
      TA =  TWO*ABORT
      PK =  ZERO
      WI =  ETA + ETA
      P  =  ZERO
      Q  =  ONE - ETA*XI
      AR = -E2MM1
      AI =  ETA
      BR =  TWO*(X-ETA)
      BI =  TWO
      DR =  BR/(BR*BR+BI*BI)
      DI = -BI/(BR*BR+BI*BI)
      DP = -XI*(AR*DI+AI*DR)
      DQ =  XI*(AR*DR-AI*DI)
   40 CONTINUE
      P  = P  + DP
      Q  = Q  + DQ
      PK = PK + TWO
      AR = AR + PK
      AI = AI + WI
      BI = BI + TWO
      D  = AR*DR - AI*DI + BR
      DI = AI*DR + AR*DI + BI
      C  = ONE/(D*D+DI*DI)
      DR =  C*D
      DI = -C*DI
      A  = BR*DR - BI*DI - ONE
      B  = BI*DR + BR*DI
      C  = DP*A  - DQ*B
      DQ = DP*B  + DQ*A
      DP = C
      IF (PK.GT.TA) GO TO 70
      IF (DABS(DP)+DABS(DQ).GE.(DABS(P)+DABS(Q))*ACC) GO TO 40
      PACCQ = HALF*ACC/DMIN1(DABS(Q),ONE)
      IF (DABS(P).GT.DABS(Q)) PACCQ = PACCQ*DABS(P)
      IF (Q.LE.ACC4*DABS(P)) GO TO 80
      GAM = (FX-P)/Q
C ***    SOLVE FOR FCM=FX AT LAMBDA=XLM, THEN FIND NORM FACTOR W=W/FCM
      W   = ONE/ DSQRT((FX-P)*GAM+Q)
C
      ENDIF
C
C
      FCM = DSIGN(W,FCL)
C ***                    = (SIGN OF FCL) * DABS(W)
      F(1) = FCM
      IF (IXLTURN.eq.0) THEN
      GCL = FCM*GAM
      ELSE
      GCL = GJWKB
      ENDIF
      G(1) = GCL
      GPL =  GCL * (P-Q/GAM)
      GPR(1) = GPL
      FPR(1) = FCM * FX
C
C
C *** UPWARD RECURRENCE FROM G(1),GPR(1)  STORED VALUE IS RL
C *** RENORMALISE F,FPR AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C ***    XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
      W    = W/DABS(FCL)
      DO L=1,LLMAX
      XL = XL + ONE
      EL = ETA/XL
      RL = G(L+1)
      SL = EL + XL*XI
      GCL1     = ( SL*GCL-GPL )/RL
      GPL      =   RL*GCL - SL*GCL1
      GCL      = GCL1
      G (L+1) = GCL1
      GPR(L+1) = GPL
      FPR(L+1) = W * FPR(L+1)
      F (L+1) = W * F (L+1)
      ENDDO
C
C GENERATE PENETRABILITY, SHIFT FACTOR, AND PHASE SHIFT, & DERIVATIVES
C NOTE THAT IEXP = 1 MEANS "NORMAL" VERSION WORKED...
      L = LNOW + 1
      IF (IEXP.GT.1) THEN
      IF (IEXP.LT.150) THEN
      ASQ = G(L)**2
      AAA = TEN ** (-IEXP*2)
      PENCOUL = XX/ASQ * AAA
      ENDIF
C        NOTE WE'VE ALREADY SET DEFAULTS TO ZERO SO NO NEED TO REPEAT
      ELSE
      ASQ = F(L)**2 + G(L)**2
      PENCOUL = XX/ASQ
      ENDIF
C
      RETURN
C-----------------------------------------------------------------------
C
C     ERROR CONDITION
C
C-----------------------------------------------------------------------
   50 IFAIL = -1
      RETURN
   60 IFAIL =  1
      RETURN
   70 IFAIL =  2
      RETURN
   80 IFAIL =  3
      RETURN
      END
      SUBROUTINE JWKB (XX, ETA1, XL, FJWKB, GJWKB, IEXP)
C=======================================================================
C
C *** COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS FOR XL.GE.0
C *** AS MODIFIED BY BIEDENHARN ET.AL. PHYS REV 97 (1955) 542-554
C *** CALLS DMAX1, DSQRT, DLOG, DEXP, DATAN2, DFLOAT, INT
C *** BARNETT FEB 1981
C
C=======================================================================
      INCLUDE 'implicit.h'
      DATA ZERO/0.0D0/
      DATA HALF/0.5D0/
      DATA ONE /1.0D0/
      DATA TEN/10.0D0/
C
      DATA ALOGE /0.434294481903251816667932D0/
C *** ALOGE IS LOG-BASE-10 OF E ==> E=DEXP(1.0D0),ALOGE=DLOG10(E)
      DATA SIX35 /0.1714285714285714285714285714285714285714285D0/
C *** SIX35 IS 6.0/35.0
C
      X     = XX
      ETA   = ETA1
      GH2   = X* (ETA+ETA-X)
      XLL1  = DMAX1 ( XL*XL + XL, ZERO )
      IF (GH2+XLL1.LE.ZERO) RETURN
      HLL = XLL1 + SIX35
      HL  = DSQRT(HLL)
      SL  = ETA/HL + HL/X
      RL2 = ONE + ETA*ETA/HLL
      GH  = DSQRT(GH2+HLL)/X
      PHI = X*GH - HALF*( HL*DLOG((GH+SL)**2/RL2) - DLOG(GH) )
      IF (ETA.NE.ZERO) PHI = PHI - ETA*DATAN2(X*GH,X-ETA)
C
      PHI10 = -PHI*ALOGE
      IEXP  = INT(PHI10)
      IF (IEXP.GT.70) THEN
      GJWKB = TEN**(PHI10-DFLOAT(IEXP))
      ELSE
      GJWKB = DEXP(-PHI)
      IEXP = 0
      ENDIF
      FJWKB = HALF/(GH*GJWKB)
      RETURN
      END
      SUBROUTINE ASYMP1 (EETA, RHOI, U, UPR)
C=======================================================================
C
C *** PURPOSE -- CALCULATE ASYMPTOTIC EXPANSION OF G0 & G0PR FROM EQS.
C ***            14.5.12B AND 14.5.13B IN ABROMOWITZ & STEGUN
C ***            NOTE THAT THESE FORMULAE ARE VALID FOR RHOI=2*ETA
C
C=======================================================================
      INCLUDE 'implicit.h'
      ETA = EETA
      CETA = (ETA)**0.3333333333333333D0
      SETA = DSQRT(CETA)
      TEMP = 1.0D0/CETA**2
      U = 1.223404016D0*SETA*
     *     (1.0D0 + TEMP**2*(.04959570165 D0 + TEMP*(-.008888888889 D0 +
     *              TEMP**2*(.002455199181D0 + TEMP*(-.0009108958061D0 +
     *              TEMP**2*.0002534684115D0)))))
      UPR = -.7078817734D0*
     *     (1.0D0 + TEMP*(-.1728260369  D0 + TEMP**2*(.0003174603174D0 +
     *              TEMP*(-.003581214850D0 + TEMP**2*(.0003117824680D0
     *            - TEMP*.0009073966427D0)))))/SETA
      RHOI  =  2.0D0*ETA
      RETURN
      END
      SUBROUTINE XSIGLL (EETA, SIGMA0, LMAX)
C=======================================================================
C
C PURPOSE -- GENERATE SIGCOOL(LL) FOR LL=1 THRU LMAX+1(IE L=0 THRU LMAX)
C
C=======================================================================
      INCLUDE 'implicit.h'
      PARAMETER (MAXL=100)
      COMMON/COOLCOM/F(MAXL),FPR(MAXL),G(MAXL),GPR(MAXL),
     1 SIGCOOL(MAXL)
      DIMENSION BER(5)
      DATA MMMXXX /100000/
      DATA SMALL /0.000001D0/
      DATA ZERO  /0.0D0/
      DATA ONE   /1.0D0/
      DATA TWO   /2.0D0/
      DATA PI    /3.141592653589793238462643D0/
      DATA EULER /0.577215664901532860606512D0/
      DATA BER /0.1666666666666666666666666666666666667D0,
     *         -0.0333333333333333333333333333333333333D0,
     *          0.0238095238095238095238095238095238095D0,
     *         -0.0333333333333333333333333333333333333D0,
     *          0.0757575757575757575757575757575757576D0/
C     DATA BER [ 1/6, -1/30, 1/42, -1/30, 5/66 ]
C
      ETA = EETA
      PETA =  DABS(ETA)
C
      IF (PETA.GE.3.0D0) THEN
C
C ***    HERE (DABS(ETA).GE.3) SO USE TABLE BER TO ESTIMATE SIGMA0
      SUM = ZERO
      DO I=1,5
      XI = I
      M = 2*I - 1
      XM = M
      SUM = SUM + BER(I) / (TWO*XI*XM*(PETA**M))
      ENDDO
      SIGMA0 = PI/4.0D0 + PETA*(DLOG(PETA)-ONE) - SUM
C ***    EQ. 14.6.16 GIVES ALL BUT "-SUM" OF THIS EQUATION ???
C
      ELSE
C
C ***    HERE DABS(ETA).LT.3 SO GENERATE SIGMA0 MORE ACCURATELY, FROM
C ***       EQ. 14.5.6 PG 650 ABRAMOWITZ & STEGUN
      SUMAS = ZERO
      DO IS=1,MMMXXX
      S = IS
      TEMP1 = PETA/S
      IF (S.LE.TWO*PETA) THEN
C ***          THIS IS EXACT FROM EQ. 6.1.27 FOR ARG(GAMMA(1+PETA/S))
C ***             USING EQ. 6.3.2 FOR PSI(1)
      AS = TEMP1 - DATAN(TEMP1)
      ELSE
C ***    THIS IS APPROXIMATION FOR ATAN(TEMP1) VIA TAYLOR EXPANSION
C ***             FOR PETA/S -> 0
      AS = ZERO
      K = 0
      DO J=1,MMMXXX
      M = J + J + 1
      XM = M
      ADD = (TEMP1**M)/XM
      IF (K.EQ.0) THEN
      AS = AS + ADD
      K = 1
      ELSE
      AS = AS - ADD
      K = 0
      ENDIF
      IF (DABS(ADD/AS).LE.SMALL) GO TO 10
      ENDDO
   10 CONTINUE
      ENDIF
      SUMAS = SUMAS + AS
      IF (DABS(AS/SUMAS).LE.SMALL) GO TO 20
      ENDDO
   20 CONTINUE
      SIGMA0 = - EULER*PETA + SUMAS
C
      ENDIF
      IF (ETA.LT.ZERO) SIGMA0 = -SIGMA0
      SIGCOOL(1) = SIGMA0
C
C *** NOW SET SIGCOOL(LL) FOR ALL L
      IF (LMAX.GT.0) THEN
      DO LL=1,LMAX
      XL = DFLOAT(LL)
      SIGCOOL(LL+1) = SIGCOOL(LL) + DATAN(ETA/XL)
      ENDDO
      ENDIF
C
      RETURN
      END
      SUBROUTINE ASYMP2 (EETA, RRHO, U, UPR, RHOI, SIGMA0)
C=======================================================================
C
C *** PURPOSE -- CALCULATE U, UPR FOR LARGE RHOI BUT FINITE ETA
C ***            USING EQS. 14.5.1-8 ON PG 540 ABROMOWITZ & STEGUN
C ***            NOTE THAT RHOI IS CHOSEN TO BE "LARGE ENOUGH" SO THAT
C ***            FORMULA IS VALID.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION ZOLD(4), ZNEW(4), Z(4), BIGZ(4)
      DATA DEL/100.0D0/
      DATA EPSLON/0.000001D0/
      DATA ZERO  /0.0D0/
      DATA ONE   /1.0D0/
      DATA TWO   /2.0D0/
      DATA TEN  /10.0D0/
C
      ETA = EETA
      RHO = RRHO
      RHOI = DMAX1 (RHO*TWO, TEN, TEN*ETA)
CX      RHOI = DMAX1 (RHO, TEN, TEN*ETA)
C
   10 CONTINUE
C *** INITIALIZE
      XN = ZERO
      ZOLD(1) = ONE
      ZOLD(2) = ZERO
      ZOLD(3) = ZERO
      ZOLD(4) = ONE - ETA/RHOI
      DO I=1,4
      Z(I) = ZOLD(I)
      BIGZ(I) = DABS(Z(I))
      ENDDO
C
      JCHECK = 0
      DO N=1,100
      TEMP = TWO*(XN+ONE)*RHOI
      AN = (TWO*XN+ONE)*ETA/TEMP
      BN = (ETA*ETA-XN*(XN+ONE))/TEMP
      XN = XN + ONE
      ZNEW(1) = AN*ZOLD(1) - BN*ZOLD(2)
      ZNEW(2) = AN*ZOLD(2) + BN*ZOLD(1)
      ZNEW(3) = AN*ZOLD(3) - BN*ZOLD(4) - ZNEW(1)/RHOI
      ZNEW(4) = AN*ZOLD(4) + BN*ZOLD(3) - ZNEW(2)/RHOI
      ICHECK = 0
      DO I=1,4
      Z(I)    = Z(I) + ZNEW(I)
      ZOLD(I) = ZNEW(I)
      TEMP2   = DABS(Z(I))
      BIGZ(I) = DMAX1 (BIGZ(I), DABS(Z(I)))
      IF (BIGZ(I)/TEMP2 .GT. DEL) GO TO 20
      IF (DABS(ZNEW(I)/Z(I)) .LE. EPSLON) ICHECK = ICHECK + 1
      ENDDO
      W = Z(1)*Z(4) - Z(2)*Z(3)
      IF (DABS(W).GT.TEN) GO TO 20
      IF (ICHECK.EQ.4) THEN
      JCHECK = JCHECK + 1
      IF (JCHECK.GE.4) GO TO 30
      ELSE
      JCHECK = 0
      ENDIF
      ENDDO
C
C
   20 CONTINUE
C OOPS.  RHOI ISN'T BIG ENOUGH FOR THE ASYMPTOTIC FORMULA TO CONVERGE.
C        DOUBLE RHOI & TRY AGAIN.
      RHOI = RHOI*2.0D0
      GO TO 10
C
   30 CONTINUE
C HERE THE FORMULA FOR Z'S HAS CONVERGED.
C ERGO CALCULATE U, UPR, RHOI:
      PHI = RHOI - ETA*DLOG(TWO*RHOI) + SIGMA0
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
      G0   = Z(1)*COSPHI - Z(2)*SINPHI
      G0PR = Z(3)*COSPHI - Z(4)*SINPHI
      U = G0
      UPR = G0PR
      RETURN
      END
      SUBROUTINE TAYLOR (EETA, RRHO, U, UPR, RHOI)
C=======================================================================
C
C *** PURPOSE -- DO TAYLOR SERIES INTEGRATION OF G0, G0PR FOR ARG=RHO,
C ***            STARTING FROM THE (NOW KNOWN) VALUES AT ARG=RHOI
C
C
C *** FIND SOLUTION OF DIFFERENTIAL EQUATION
C         U" + (1-2*ETA/RHO) U = 0    (COULOMB EQN FOR L=0, 14.1.1 A&S)
C     VIA TAYLOR EXPANSION
C         U(RHO) = U + U' D + U" D**2/2 + U"' D**3/6 + U"" D**4/4! + ...
C     WHERE RIGHT-HAND-SIDE IS EVALUATED AT RHOI
C     AND WHERE D = DELTA = RHO-RHOI
C
C     REWRITING U(RHO) = SUM  A(N)  WHERE N=1 TO INFINITY, WITH
C
C          A(1) = U
C          A(2) = U'  D
C          A(3) = U"  D**2 / 2!
C          A(4) = U"' D**3 / 3!
C          A(5) = U"" D**4 / 4!
C          ...
C          A(N+1) = U[N] D**N / N!
C
C
C     SUBSTITUTING EQUATION (1) INTO A(3) GIVES
C     A(3) = {-(1-2 ETA/RHOI) U } D**2/2
C          = - (1-2 ETA/RHOI) D**2 / 2 * A(1)
C     A(4) = - { 2 ETA/RHOI**2 U + (1-2 ETA/RHOI) U'} D**3/ 3!
C          = - { 2 ETA/RHOI**2 A(1) + (1-2 ETA/RHOI) A(2) /D} D**3/3!
C          = - {   2 ETA/RHOI**2  D**3   /3! A(1)
C                + (1-2 ETA/RHOI) D**2 2!/3! A(2)  }
C
C          = - {     1/RHOI       D    /3     A(3)
C                + (1-2 ETA/RHOI) D**2 /(3*2) A(2)
C                +                D**3 /(3*2) A(1) }  USING (3) ABOVE
C
C    A(5) = {4ETA/RHOI**3 U - 4ETA/RHOI**2 U'-(1-2 ETA/RHOI) U"}D**4/4!
C         = - {   2/RHOI        D**3   /4! A(2)
C              + (1-2 ETA/RHOI) D**2 2!/4! A(3)
C              +  2/RHOI        D    3!/4! A(4) }  USE A(4) TO REPLACE A
C
C     A(N) = - {    1/RHOI       D    (N-3) /(N-1)        A(N-1)
C               + (1-2 ETA/RHOI) D**2       /((N-1)(N-2)  A(N-2)
C               +                D**3       /((N-1)(N-2)  A(N-3) }
C
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION A(100)
      DATA EPSLON/1.0D-6/
      DATA BIGGER/1.0D10/
      DATA BIGGST/1.0D30/
      DATA DEL   /100.D0/
      DATA ZERO  /0.0D0/
      DATA ONE   /1.0D0/
      DATA TWO   /2.0D0/
C
      ETA = EETA
      RHO = RRHO
      DELTA = RHO - RHOI
      IF (DELTA.EQ.ZERO) RETURN
C
   10 CONTINUE
      A(1) = U
      A(2) = DELTA*UPR
      A(3) = - DELTA*DELTA/TWO*(ONE-TWO*ETA/RHOI)*A(1)
      NSTART = 4
   20 CONTINUE
      JCHECK = 0
      SUM   = ZERO
      SUMPR = ZERO
      BIG   = ZERO
      BIGPR = ZERO
      DO N=1,100
      XN = DFLOAT(N-1)
      IF (N.GE.NSTART) THEN
      A(N) = - ( DELTA*(XN-ONE)*(XN-TWO  )*A(N-1) +
     *               (  RHOI-TWO*ETA)*(DELTA**2)*A(N-2) +
     *                                (DELTA**3)*A(N-3)    ) /
     *               (RHOI*(XN-ONE)*XN)
      IF (A(N).GT.BIGGER) GO TO 30
      ENDIF
      SUM   = SUM   +    A(N)
      SUMPR = SUMPR + XN*A(N)
      IF (SUM  .GE.BIGGST) GO TO 30
      IF (SUMPR.GE.BIGGST) GO TO 30
      BIG   = DMAX1 (BIG  , DABS(SUM)  )
      BIGPR = DMAX1 (BIGPR, DABS(SUMPR))
      IF (SUM.EQ.ZERO .OR. SUMPR.EQ.ZERO) THEN
      JCHECK = 0
      ELSE
C        ELSE IF (SUM.NE.ZERO .AND. SUMPR.NE.ZERO) THEN
      IF (DABS(BIG  /SUM  ).GE.DEL) GO TO 30
      IF (DABS(BIGPR/SUMPR).GE.DEL) GO TO 30
      IF (DABS(A(N) /SUM  ).GE.EPSLON .OR.
     *          DABS(XN*A(N)/SUMPR).GE.EPSLON) THEN
      JCHECK = 0
      ELSE
      JCHECK = JCHECK + 1
      IF (JCHECK.GE.4) GO TO 40
      ENDIF
      ENDIF
      ENDDO
      N = 100
C
   30 CONTINUE
C *** THE SERIES DID NOT CONVERGE FOR U(RHO) USING TAYLOR EXPANSION
C ***   AROUND ARG=RHOI.  ERGO, FIND U(RHOI+X), WHERE X=DELTA/2, USING
C ***   TAYLOR EXPANSION AROUND ARG=RHOI.  ASSUMING THIS CONVERGES,
C ***   REDEFINE RHOI TO BE RHOI + X AND TRY AGAIN TO GET U(RHO).
      NSTART = MAX0 (NSTART,N+1)
      M = NSTART - 1
      DELTA = DELTA/TWO
      TEMP = TWO
      DO N=1,M
      TEMP = TEMP/TWO
      A(N) = A(N)*TEMP
      ENDDO
      GO TO 20
C
   40 CONTINUE
C *** HERE WE KNOW SUM = U(RHOI+DELTA).  REDEFINE RHOI & DELTA AND TRY
C ***   AGAIN TO FIND U(RHO) AS TAYLOR EXPANSION AROUND ARG=RHOI.
      U = SUM
      UPR = SUMPR/DELTA
      RHOI = RHOI + DELTA
      DELTA = RHO - RHOI
      IF (DABS(DELTA).GE.EPSLON) GO TO 10
      RETURN
      END
      SUBROUTINE GETFG (ETA, RHO, LLMAX, LNOW, G0, G0PR)
C=======================================================================
C
C PURPOSE -- CALCULATE G & F FOR ALL L, GIVEN G(L=0) AND DERIV G(L=0)
C
C=======================================================================
      INCLUDE 'implicit.h'
      PARAMETER (MAXL=100)
      COMMON/COOLCOM/F(MAXL),FPR(MAXL),G(MAXL),GPR(MAXL),
     1 SIGCOOL(MAXL)
      DATA ONE   /1.0D0/
      DATA TWO   /2.0D0/
      DATA THREE /3.0D0/
      LMAX = LLMAX
C
C *** SET G(L)
      LIMIT =  MAX0 (3, LMAX+1)
      G  (1) = G0
      GPR(1) = G0PR
      G  (2) = ((ETA+ONE/RHO)*G(1)-GPR(1))/DSQRT(ETA**2+ONE)
C ***        FROM ABRAMOWITZ & STEGUN EQ. 14.2.2
      DO L=3,LIMIT
      XL = L - 1
      TEMP1 = DSQRT(XL*XL+ETA*ETA)
      G(L) = (TWO*XL-ONE)/TEMP1 *
     *          (ETA/(XL-ONE)+XL/RHO)*G(L-1) -
     *          XL/TEMP1 * DSQRT(ONE+(ETA/(XL-ONE))**2)*G(L-2)
C ***    FROM ABRAMOWITZ & STEGUN EQ. 14.2.3 REWRITTEN
      IF (DABS(G(L)).GT.1.0D+12 .AND. L.GT.LNOW) THEN
      LIMIT = L
      GO TO 10
      ENDIF
      ENDDO
C
   10 CONTINUE
C *** FIND MAXIMUM L VALUE TO USE IN FIGURING F(LIMIT)
C *** I.E. FIND J SUCH THAT G(J) < 1.E-4 * G(J-1) THREE
C *** TIMES IN A ROW
      GM2 = G(LIMIT-1)
      GM1 = G(LIMIT  )
      IL = - 1
      DO J=LIMIT,10000
      XL = J
      TEMP1 = DSQRT(XL*XL+ETA*ETA)
      GM = (TWO*XL-ONE)/TEMP1*(ETA/(XL-ONE)+XL/RHO)*GM1 -
     *          XL/TEMP1* DSQRT(ONE+(ETA/(XL-ONE))**2)*GM2
C ***    AGAIN EQ. 14.2.3
C GMH    IF ( (G(LIMIT)/GM)**2 .GT. 1.0D-8 ) IL = -2
      IF ( DABS(G(LIMIT)/GM) .GT. 1.0D-4 ) IL = -2
      IF (IL.GT.0) GO TO 20
      IL = IL + 1
      GM2 = GM1
      GM1 = GM
      ENDDO
C
   20 CONTINUE
C *** FIGURE APPROXIMATE F(LIMIT+3), F(LIMIT+4), ... F(J-1) IN REVERSE
C ***    ORDER.  DO NOT STORE THESE NUMBERS ANYWHERE PERMANENT
      XL = J
      FP1 = XL/GM/DSQRT(XL**2+ETA**2)
C ***       FROM ABRAMOWITZ & STEGUN EQ. 14.2.5 WITH
C ***                         F(L)G(L-1) ASSUMED TO BE NEGLIGIBLE
      FP2 = 0
      L = J - 1
      DO LL=1,J-3-LIMIT
      L = L - 1
      XL = L
      TEMP2 = DSQRT( (XL+ONE)**2 + ETA**2 )
      FP = ( (TWO*XL+THREE)*(ETA/(XL+TWO)+ (XL+ONE)/RHO )*FP1 -
     *        (XL+ONE) * DSQRT(ONE+(ETA/(XL+TWO))**2)*FP2)/ TEMP2
C ***    FROM ABRAMOWITZ & STEGUN EQ. 14.2.3 AGAIN, IN REVERSE
      FP2 = FP1
      FP1 = FP
      ENDDO
C
C *** NOW FIGURE F(1), F(2), ..., F(LIMIT), F(LIMIT+1),
C ***      AGAIN IN REVERSE ORDER.  STORE THESE IN ARRAY F(L)
      DO LL=1,LIMIT+2
      L = L - 1
      XL = L
      TEMP2 = DSQRT( (XL+ONE)**2 + ETA**2 )
      FP = ( (TWO*XL+THREE)*(ETA/(XL+TWO) + (XL+ONE)/RHO )*FP1 -
     *        (XL+ONE) * DSQRT(ONE+(ETA/(XL+TWO))**2)*FP2)/ TEMP2
      F(L+1) = FP
      FP2 = FP1
      FP1 = FP
      ENDDO
C
      CONTINUE
C *** GENERATE DERIVATIVE FUNCTIONS FPR & GPR FOR L=1,LIMIT
C ***      (REMEMBER, WE ALREADY HAVE GPR(1))
      FPR(1) = (ONE/RHO+ETA) * F(1) - DSQRT(ONE+ETA**2)*F(2)
      DO L=2,LIMIT
      XL = L
      FPR(L) = ( XL/RHO + ETA/XL ) * F(L)   -
     *       DSQRT(ONE+(ETA/XL)**2)   * F(L+1)
      TEMP1 = ETA/(XL-ONE)
      GPR(L) = DSQRT(ONE+TEMP1**2)  *G(L-1)
     *           - ((XL-ONE)/RHO+TEMP1)*G(L)
      ENDDO
C
      IF (LIMIT.LE.LMAX) LLMAX = LIMIT - 1
      RETURN
      END
      SUBROUTINE BIGETA (EETA,RRHO,LNOW,LLMAX,PENCOUL)
C=======================================================================
C
C *** FORMULAE 14.6.7-8 PAGE 542 IN ABRAMOWITZ & STEGUN, FOR ETA >> RHO
C
C=======================================================================
      INCLUDE 'implicit.h'
      PARAMETER (MAXL=100)
      COMMON/COOLCOM/F(MAXL),FPR(MAXL),G(MAXL),GPR(MAXL),
     1 SIGCOOL(MAXL)
      DATA HALF /0.5D0/
      DATA ONE  /1.0D0/
      DATA TWO  /2.0D0/
      DATA EULER /0.577215664901532860606512D0/
      DATA PI    /3.141592653589793238462643D0/
C
      ETA = EETA
      RHO = RRHO
      Q = TWO*RHO*ETA
      ZHALF = DSQRT(Q)
      Z = ZHALF * TWO
C
C *** EVALUATE I_0(Z) FROM 9.6.12 A&S
      SUM = ONE
      A   = Q
      DO K=1,100
      IF (SUM+A.EQ.SUM) GO TO 20
      SUM = SUM + A
      A = A * Q/DFLOAT(K+1)**2
      ENDDO
      WRITE(3,10)
      WRITE(*,10)
   10 FORMAT(///' ERROR - STOP IN BIGETA IN RML/MRML08.F    # 1'///)
      CALL ENDERROR
   20 CONTINUE
      AI0 = SUM
C
C *** EVALUATE K_0(Z) FROM 9.6.13 A&S
      SUM = - (DLOG(ZHALF)+EULER) * AI0
      A = Q
      B = ONE
      DO K=1,100
      IF (SUM+A*B.EQ.SUM) GO TO 40
      SUM = SUM + A*B
      B = B + ONE/DFLOAT(K+1)
      A = A * Q/DFLOAT(K+1)**2
      ENDDO
      WRITE(*,30)
      WRITE(3,30)
   30 FORMAT(///' ERROR - STOP IN BIGETA IN RML/MRML08.F    # 2'///)
      CALL ENDERROR
   40 CONTINUE
      AK0 = SUM
C
C *** EVALUATE I_1(Z) FROM 9.6.10 A&S
      SUM = ONE
      A   = Q/TWO
      DO K=1,100
      IF (SUM+A.EQ.SUM) GO TO 60
      SUM = SUM + A
      A = A * Q/(DFLOAT(K+1)*DFLOAT(K+2))
      ENDDO
      WRITE(*,50)
      WRITE(3,50)
   50 FORMAT(///' ERROR - STOP IN BIGETA IN RML/MRML08.F    # 3'///)
      CALL ENDERROR
   60 CONTINUE
      AI1 = SUM * ZHALF
C
C *** EVALUATE K_1(Z) FROM 9.6.11 A&S
      SUM = ONE/Z + (DLOG(ZHALF)+EULER)*AI1 - ZHALF*HALF
      A = ZHALF * Q/TWO
      B = ONE
      C = HALF**2
      DO K=1,100
      IF (SUM-A*(B+C).EQ.SUM) GO TO 80
      SUM = SUM - A*(B+C)
      C = HALF/DFLOAT(K+2)
      B = B + ONE/DFLOAT(K+1)
      A = A * Q/(DFLOAT(K+1)*DFLOAT(K+2))
      ENDDO
      WRITE(*,70)
      WRITE(3,70)
   70 FORMAT(///' ERROR - STOP IN BIGETA IN RML/MRML08.F    # 4'///)
      CALL ENDERROR
   80 CONTINUE
      AK1 = SUM
C
C *** NOW HAVE I_0(Z), I_1(Z), K_0(Z), K_1(Z)
C *** ERGO, CAN GET F,G,FP,GP FOR L = 0, FROM EQ. 14.6.8 A&S
C
      C = DSQRT(PI*RHO)
      D = DSQRT(TWO*PI*ETA)
      F  (1) = C * AI1
      FPR(1) = D * AI0
      C = TWO*C/PI
      D = TWO*D/PI
      G  (1) =   C * AK1
      GPR(1) = - D * AK0
C
C *** STORE RESULTS BECAUSE THEY'LL BE CHANGED IN GETFG
      G0 = G(1)
      G0PR = GPR(1)
C
      IF (LNOW.GT.0) THEN
C ***    OBTAIN VALUES FOR ALL L'S
      CALL GETFG (ETA, RHO, LLMAX, LNOW, G0, G0PR)
      ENDIF
      A = PI*ETA
      B = DEXP(-A)
      N = LNOW + 1
      IF (F(N)*B*B+G(N).EQ.G(N)) THEN
      PENCOUL = ( RHO*B /G(N))/G(N) * B
C        NO NEED TO ZERO THE OTHERS BECAUSE THAT'S ALREADY BEEN DONE
      ELSE
      F(N) = F(N)*B
      G(N) = G(N)/B
      FPR(N) = FPR(N)*B
      GPR(N) = GPR(N)/B
c-----put GETPS inline
      LPLUS1  = LNOW + 1
      ASQ     = F(LPLUS1)**2 + G(LPLUS1)**2
      A       = DSQRT(ASQ)
      PENCOUL = RHO/ASQ
      ENDIF
      RETURN
      END
      SUBROUTINE RDRML
C=======================================================================
C
C
C     NEW (2003) REICH-MOORE WITH COMPETITION
C
C     READ REICH-MOORE DATA FOR ONE ENERGY RANGE.
C     ALL DATA WILL BE TREATED AS A SINGLE SECTION.
C
C     THIS ROUTINE USES THE GENERAL REICH-MOORE FORMALISM WITH TWO
C     FISSION CHANNELS AS DEFINED IN ENDF/B-IV. THIS ROUTINE WILL BE
C     USED TO READ DATA IN ANY VERSION OF THE ENDF/B FORMAT (NOTE,
C     THE ENDF/B-VI REICH-MOORE FORMAT HAS NOW BEEN UPDATED TO BE
C     EXACTLY THE SAME AS EARLIER VERSIONS OF THE ENDF/B FORMAT).
C
C     CHECK FOR PRELIMINARY ENDF/B-VI FORMAT(NOW ABANDONED). TERMINATE
C     EXECUTION IF DATA IS IN PRELIMINARY ENDF/B-VI FORMAT.
C
C     FIELD DEFINITIONS FOR EACH RESONANCE
C
C     FIELD          PRELIMINARY          CURRENT
C                    ENDF/B-VI FORMAT     ENDF/B-VI FORMAT
C     =====          ================     ================
C       1            ENERGY               ENERGY
C       2            J                    J
C       3            TOTAL WIDTH          ELASTIC WIDTH
C       4            ELASTIC WIDTH        CAPTURE WIDTH
C       5            CAPTURE WIDTH        FISSION WIDTH 1
C       6            NOT USED             FISSION WIDTH 2
C                    (FISSION NOT
C                    ALLOWED)
C
C     IF THIRD FIELD (PRELIMINARY TOTAL) IS EQUAL TO THE SUM OF THE
C     FOURTH (PRELIMINARY ELASTIC) AND FIFTH (PRELIMINARY CAPTURE)
C     AND SIXTH FIELD IS ZERO (FISSION NOT ALLOWED IN PRELIMINARY
C     FORMAT) FOR ALL RESONANCES THIS PROGRAM WILL ASSUME THAT THE
C     DATA IS IN THE PRELIMINARY REICH-MOORE FORMAT AND TERMINATE
C     EXECUTION WITH A WARNING MESSAGE THAT THE DATA MUST BE CONVERTED
C     TO THE CURRENT REICH-MOORE FORMAT BEFORE IT CAN BE PROCESSED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
c-----Set FIXED dimensions for necessary arrays
      CALL Setdim
c-----Read and store particle-pair information
      CALL READPART
c-----Read and store resonance parameters
      CALL READRESP
C-----Check quantum numbers for consistency
      CALL CHEKQUAN
C-----Generate Zke, Zkte, Zkfe, Zeta7 for use in calculating k, rho,& et
      CALL Fxradi
C-----Generate energy-independent channel constants; count parameters
      CALL Betset
      RETURN
      END
      SUBROUTINE same1(Awrx,Nrox,Napsx)
c-----2014/2/20 - ereslow, ereshigh passed through common
C=======================================================================
C
C     DEFINE PARAMETERS READ BY RECENT FOR SAMRML
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/INDATS/ZA,AWR,ZAI,ABN,SPI,AP,AWRI,QX,NRS,NIS,LFW,NER,
     1 LRU,LRF,LRFIN,NLS
      COMMON/INDATD/ELX,EHX
      COMMON/NAPRHO/NRO,NAPS
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
c
c     First check for ERRORS.
c
      if(NIS.ne.1) then
      write(*,10)
      write(3,10)
   10 format(///' ERROR - Only one isotope allowed for LRF=7 data.'///)
      CALL ENDERROR
      endif
c
c     O.K. - set parameters.
c
      LRF7     = LRF7+1 ! LRF7 = # of LRF=7 Energy Ranges.
      Awrx     = AWR
      Nrox     = NRO
      Napsx    = NAPS
      ereslow  = ELX
      ereshigh = EHX
      RETURN
      ENTRY same2(NLSX)
C=======================================================================
C
C     DEFINE PARAMETERS READ BY SAMRML FOR RECENT
C
C=======================================================================
      NLS  = NLSX
      RETURN
      ENTRY same0
C=======================================================================
C
C     Reset LRF=7 used flag for each evaluation.
C
C=======================================================================
      LRF7    = 0   ! # of LRF=7 Sections
      IMFISSY = 0
c-----2014/2/20 - added for multi-region
      KresSum = 0
      KgrSum  = 0
      NppSum  = 0
      RETURN
      END
      SUBROUTINE same3(NRESX,Npp)
C=======================================================================
C
C     DEFINE PARAMETERS READ BY SAMRML FOR RECENT
C
c     1 dummy section of Lrf = 7
c
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'recent.h'
      COMMON/INDATD/ELX,EHX
      COMMON/MAXIE/NEDSEC,NEDNOD,NEDRES
      COMMON/ETABS/NODES
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
      NSECT       = NSECT + 1
      EL(NSECT)   = ELX
      EH(NSECT)   = EHX
      NGRTAB(1,NSECT) = Ngr1
      NGRTAB(2,NSECT) = Ngr2
      NPPTAB(  NSECT) = Npp
      MODE(NSECT) = 7
      NEDSEC      = 1
      NEDNOD      = NODES
      NEDRES      = NRESX
      if(NSECT.eq.1) then
      ii1 = 1
      else
      ii1 = NHIGH(NSECT-1) + 1
      endif
      NLOW(NSECT)  = ii1
      NHIGH(NSECT) = ii1 + NRESX - 1
      return
      END
      SUBROUTINE SIGRML(Enow)
c=======================================================================
c
c     Calculate ALL cross sections at energy ENOW
c
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
c-----define RECENT energy for SAMRML
      Su = Enow
C-----Generate energy-dependent pieces of cross section (no derivs)
      CALL Abpart
C-----Form the cross section Crss (no Deriv)
      CALL Crosss
C-----Final answer
      CALL Answers
      return
      END
      SUBROUTINE READPART
c=======================================================================
C
C *** PURPOSE -- Read particle-pair definitions for Lrf=7
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
c-----limit of resolved region for SAMRML
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
      common/usemt/dumrml(100),Quse(11),Mtuse(11)
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
      character*1 FIELD(11,12)
      DATA Aneutron /1.00866491578d0/
      DATA Proton   /1.00727646688d0/
      DATA Deuteron /2.01355321271d0/
      DATA Alpha    /4.0015061747d0/
      DATA He3      /3.00000000000d0/
      DATA Triton   /3.0000000000d0/
      DATA One  /1.0d0/
C
C *** NOTE: Above the last two ATOMIC WEIGHTS are not correct!
C
      CALL CARDIO(C1, C2, Lx, Lx, Nls, Lx)
c-----initialize - only on first region
      if(NppSum.le.0) then
      do i=1,11
      Quse    (i) = 0.0d+0
      Mtuse   (i) = 0
      QuseSum (i) = 0.0d+0
      MtuseSum(i) = 0
      ENDDO
      endif
c-----define parameters read by SAMRML for RECENT
      call same2(Nls)
c-----and for SAMRML
      Ngroup = Nls
      CALL CARDIO(C1, C2, Npp, L2, N1, N2)
C-----Spi and Ap fields are not used for LRF=7
      write(*,10) Ngroup,Npp
      write(3,10) Ngroup,Npp
   10 format(1x,78('=')/
     3 ' Number of J Values-------------------',I11/
     4 ' Number of Particle-Pairs-------------',I11)
C
C ***    Read information re particle-pairs
      if(IMEDIT.ne.0) then
      Write(3,20)
   20 format(1x,78('=')/' Particle Pair Information'/1x,78('=')/
     1       '   I',6X,'MassA',6X,'MassB',9X,'ZA',9X,'ZB',6X,'SpinA',
     1         6X,'SPinB'/
     2         4X,10X,'Q',5X,'Pentra',6X,'Shift',9X,'MT',4X,'ParityA',
     3         4X,'ParityB'/1x,78('-'))
      endif
      DO I=1,Npp
      CALL LISTIO(dumrml(1),12)
      Ema    (I)   = dumrml(1)
      Emb    (I)   = dumrml(2)
      Za           = dumrml(3)
      Zb           = dumrml(4)
      Spina  (I)   = dumrml(5)
      Spinb  (I)   = dumrml(6)
      Qqq    (I)   = dumrml(7)
      XLpent       = dumrml(8)
      xIshift      = dumrml(9)
      aMT          = dumrml(10)
      Pa     (I)   = dumrml(11)
      Pb     (I)   = dumrml(12)
      Ishift (I)   = xIshift + 0.01
      Lpent  (I)   = xLpent + 0.01
      Kza    (I)   = Za
      Kzb    (I)   = Zb
      Mt7     (I)  = aMt
      Quse   (I)   = Qqq(I)
      Mtuse  (I)   = Mt7(I)
c-----Accumulate MT from ALL regions
      if(NppSum.gt.0) then
      do k=1,NppSum
      if(MtuseSum(k).eq.Mtuse(I)) go to 30
      enddo
      endif
c-----New MT
      NppSum = NppSum + 1
      MtuseSum(NppSum) = Mtuse(I)
      QuseSum (NppSum) = Quse(I)
   30 continue
      call OUT9(Ema(I)  ,FIELD(1,1))
      call OUT9(Emb(I)  ,FIELD(1,2))
      call OUT9(Spina(I),FIELD(1,3))
      call OUT9(Spinb(I),FIELD(1,4))
      call OUT9(Qqq(I)  ,FIELD(1,5))
      call OUT9(Pa(I)   ,FIELD(1,6))
      call OUT9(Pb(I)   ,FIELD(1,7))
      if(IMEDIT.ne.0) then
      write( 3,40) I,((FIELD(II,JJ),II=1,11),JJ=1,2), Kza(I),
     1 Kzb(I),((FIELD(II,JJ),II=1,11),JJ=3,4),(FIELD(II,5),II=1,11),
     2 Lpent(I),Ishift(I),Mt7(I),((FIELD(II,JJ),II=1,11),JJ=6,7)
   40 FORMAT(  I4,22A1    ,I11,I11,22A1/
     1         4X,11A1   ,I11,I11,I11,22A1)
      endif
      ENDDO
C
C     Fill in the defaults
C
      DO I=1,Npp
      IF (Mt7(I).EQ.102) THEN
C ***       Mt7 =102 => one "particle" is gamma
      ELSE IF (Mt7(I).EQ.2) THEN
C ***       Mt7 =  2 => neutron
      Ema(I) = One
      Spina(I) = 0.5d0
      Kza(I) = 0
      Lpent(I) = 1
      ELSE IF (Mt7(I).EQ.18) THEN
C ***       Mt7 = 18 => fission
      Lpent(I) = 0
      ELSE IF (Mt7(I).GT.50 .AND. Mt7(I).LT.99) THEN
C ***  50 < Mt7 < 99 => inelastic
      Ema(I) = One
      Spina(I) = 0.5d0
      Kza(I) = 0
      Lpent(I) = 1
      ELSE IF (Mt7(I).EQ.103) THEN
C ***       Mt7 =103 => proton
      Ema(I) = Proton/Aneutron
      Spina(I) = 0.5d0
      Kza(I) = 1
      Lpent(I) = 1
      ELSE IF (Mt7(I).EQ.104) THEN
C ***       Mt7 =104 => deuteron
      Ema(I) = Deuteron/Aneutron
      Spina(I) = 1.0d0
      Kza(I) = 1
      Lpent(I) = 1
      ELSE IF (Mt7(I).EQ.105) THEN
C ***       Mt7 =105 => triton
      Ema(I) = Triton/Aneutron
      Spina(I) = 0.5d0
      Kza(I) = 1
      Lpent(I) = 1
      ELSE IF (Mt7(I).EQ.106) THEN
C ***       Mt7 =106 => He3
      Ema(I) = He3/Aneutron
      Spina(I) = 0.5d0
      Kza(I) = 2
      Lpent(I) = 1
      ELSE IF (Mt7(I).EQ.107) THEN
C ***       Mt7 =107 => alpha
      Ema(I) = Alpha/Aneutron
      Spina(I) = 0.0d0
      Pa(I) = 0.0d0
      Kza(I) = 2
      Lpent(I) = 1
      ELSE
      ENDIF
      ENDDO
      RETURN
      END
      SUBROUTINE READRESP
c=======================================================================
C
C *** Purpose -- Read the channel information and the resonance paramete
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
      common/usemt/dumrml(100),Quse(11),Mtuse(11)
      character*1 FIELD(11,12)
      character*4 MTBCD
      DATA MTBCD/'MT='/
c-----2014/2/20 - multi-region - sum groups and resonances
      Ngr1   = KgrSum + 1
      Ngr2   = KgrSum + Ngroup
c-----Insure groups are in legal range.
      if(Ngr2.gt.MaxNgroup) then
      write(*,10) Ngr2,MaxNgroup
      write(3,10) Ngr2,MaxNgroup
   10 format(//' ERROR - Ngr2=',i5,' Exceeeds MaxNgroup=',i5/
     1         '         Increase MaxNgroup and try again.'/
     2         '         Execution Terminated.'///)
      CALL ENDERROR
      endif
      KgrSum = Ngr2
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      DO Igroup=Ngr1,Ngr2
C
      CALL CARDIO(Spin(Igroup), Parity(Igroup),L1, L2, N1, IchanP1)
      call OUT9(Spin(Igroup)  ,FIELD(1,1))
      call OUT9(Parity(Igroup),FIELD(1,2))
      write( 3,20) Igroup, ((FIELD(II,JJ),II=1,11),JJ=1,2), IchanP1
      write( *,20) Igroup, ((FIELD(II,JJ),II=1,11),JJ=1,2), IchanP1
   20 format(1x,78('=')/
     1       ' J Value #----------------------------',I11/
     1       ' J Value------------------------------',11A1/
     1       ' Parity-------------------------------',11A1/
     1       ' Number of Spin-Pairs-----------------',I11)
C
      Ichan         = IchanP1 - 1
      Nchan(Igroup) = Ichan
C ***    Read channel parameters
      Ix     = 0
      I      = 0
      Kgamma = 0
      if(IMEDIT.ne.0) then
      Write(3,30)
   30 FORMAT(1x,78('-')/
     1       ' Spin-Pair#',4X,'L Value',3X,'    Spin',8X,'Bnd',
     1       8X,'Rde',8X,'Rdt')
      endif
c-----Ich = Channel = width count - DO NOT OFFSET BY SUM
      DO Ich=1,IchanP1
      CALL LISTIO(dumrml(1),6)
      Pp  = dumrml(1)
      eL  = dumrml(2)
      Spx = dumrml(3)
      Bnd = dumrml(4)
      Rde = dumrml(5)
      Rdt = dumrml(6)
      Mpp = Pp
      if(IMEDIT.ne.0) then
      call OUT9(eL ,FIELD(1,1))
      call OUT9(Spx,FIELD(1,2))
      call OUT9(Bnd,FIELD(1,3))
      call OUT9(Rde,FIELD(1,4))
      call OUT9(Rdt,FIELD(1,5))
      write( 3,40) Mpp,((FIELD(II,JJ),II=1,11),JJ=1,5)
   40 FORMAT(i11,55A1)
      endif
      Ippx = Pp
      IF (Ippx.EQ.1) THEN
      Kgamma =  Ich
      ELSE
      I = I + 1
      Ipp(I,Igroup) = Ippx
      IF (Ipp(I,Igroup).EQ.2) Ix = Ix + 1
C                  Pp=2 =>incident channel, except really it's Mt7(Ipp)=
C                         that defines what's an incident channel...
      Lspin  (I,Igroup) = eL
      Chspin7(I,Igroup) = Spx
      Bndry  (I,Igroup) = Bnd
      Rdeff  (I,Igroup) = Rde
      Rdtru  (I,Igroup) = Rdt
      ENDIF
      ENDDO
      IF (Ix.EQ.0) THEN
c-----2016/3/8 - Added ERROR message
      write(*,50)
      write(3,50)
   50 format(//' ERROR - At least one Spin-Pair# MUST = 2'/
     1         '         Check and Correct Data and Try Again.'/
     2         '         Execution Terminated.'///)
c-----2016/3/8 - Added ERROR message
      CALL ENDERROR
      ENDIF
C
C ***    Read number of resonances for this spin group
c
c      Read resonance parameters
c
      write(3,60)
   60 format(1x,78('-'))
      CALL CARDIO(C1, C2, L1, Nresg(Igroup), N1, N2)
      write( 3,70) Nresg(Igroup)
      write( *,70) Nresg(Igroup)
   70 format(' Number of Resonances-----------------',i11)
C
C ***    Read and store the resonance parameters for this spin group
      IF (Nresg(Igroup).EQ.0) THEN
c-----skip 1 line
      CALL LISTIO(dumrml(1),1)
      ELSE
c-----read resonances
      do i=1,MaxNchan
      if(Ipp(i,igroup).le.0) go to 80
      enddo
      i = MaxNchan + 1
   80 i = i - 1
      if(IMEDIT.ne.0) then
      write(3,90) (MTBCD,Mtuse(Ipp(k,igroup)),k=1,i)
   90 format(1x,78('=')/' Reich-Moore Resonance Parameters'/1x,78('=')/
     1       '     Energy'/'       (eV)',
     1       5X,'MT=102',11(5x,A3,i3)/1x,78('='))
      WRITE(3,100)
  100 FORMAT(1X,78('='))
      endif
      DO Ires=1,Nresg(Igroup)
      KresSum = KresSum + 1
c-----Check available storage
      if(KresSum.gt.MaxNres) then
      write(*,110) KresSum,MaxNres
      write(3,110) KresSum,MaxNres
  110 format(//' ERROR - Number of resonances',i8,' exceeds maximum',
     1 ' allowed',i8/
     3         '         Increase MaxNres and try again.'/
     2         '         Execution Terminated.'///)
      call ENDERROR
      endif
c-----2014/2/20 - Initialize = energy + max. # of widths.
      do k=1,MaxNchan+1
      dumrml(k) = 0.0d+0
      enddo
      CALL LISTIO(dumrml(1),Ichan+2)
      Eres  (KresSum)  = dumrml(1)
      Gamgam(KresSum)  = dumrml(2)
      widtot        = dabs(dumrml(2))
      if(Ichan.gt.0) then
      do I=1,Ichan
      Gamma(I,KresSum) = dumrml(I+2)
      widtot        = widtot + dabs(dumrml(I+2))
      call OUT9(Gamma(I,KresSum),FIELD(1,I+2))
      enddo
      endif
      call OUT9(Eres(KresSum)  ,FIELD(1,1))
      call OUT9(Gamgam(KresSum),FIELD(1,2))
      if(IMEDIT.ne.0) then
      write(3,120) ((FIELD(K,I),K=1,11),I=1,Ichan+2)
  120 format(132A1)
      endif
c-----use resonance energy as a node.
      call NOODLE(Eres(KresSum),widtot,ereslow,ereshigh)
c
      IF (Kgamma.NE.1) THEN
      X = Gamma(Kgamma+1,KresSum)
      If (Kgamma.GT.1) THEN
      DO I=Kgamma,2,-1
      Gamma(I,KresSum) = Gamma(I-1,KresSum)
      ENDDO
      ENDIF
      Gamma(1,KresSum) = Gamgam(KresSum)
      Gamgam(KresSum)  = X
      ENDIF
      ENDDO
      ENDIF
      ENDDO
C
c-----sum # of resonances
c-----2014/2/20 - multi-region sum
      if(Ngr1.gt.1) Nresg(Ngr1) = Nresg(Ngr1-1) + Nresg(Ngr1)
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      IF (Ngr2.GT.Ngr1) THEN
      DO Igroup=Ngr1+1,Ngr2
      Nresg(Igroup) = Nresg(Igroup-1) + Nresg(Igroup)
      ENDDO
      ENDIF
c-----Nres = sum of # of resonances
      Nres = Nresg(Ngr2)
c-----KresSum = one-by-one count of resonances
      IF (KresSum.NE.Nresg(Ngr2)) then
      write(3,130) Ngroup,Nresg(Ngr2),KresSum
      write(*,130) Ngroup,Nresg(Ngr2),KresSum
  130 format(' ERROR - Nresg/KresSum=',3i8,' MUST be equal.'/
     1       '         Execution Terminated.'///)
      CALL ENDERROR
      endif
c-----------------------------------------------------------------------
c
c     Summary of sequences
c
c-----------------------------------------------------------------------
      write(*,140)
      write(3,140)
  140 format(1x,78('=')/' Summary of Sequences (MT #s used)'/1x,78('-'))
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      DO Igroup=Ngr1,Ngr2
      do i=1,MaxNchan
      if(Ipp(i,igroup).le.0) go to 150
      enddo
      i = MaxNchan + 1
  150 i = i - 1
      write(*,160) igroup,(Mtuse(Ipp(k,igroup)),k=1,i)
      write(3,160) igroup,(Mtuse(Ipp(k,igroup)),k=1,i)
  160 format(10i5)
      ENDDO
      write(*,170)
      write(3,170)
  170 format(1x,78('='))
c
c     set up RECENT to use 1 section
c
      call same3(Nres,Npp)
      RETURN
      END
      SUBROUTINE CHEKQUAN
c=======================================================================
C
C *** PURPOSE -- Report quantum number information etc when particle-pai
C ***            definitions are given
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
      DATA Zero /0.0d0/
      DATA Half /0.5d0/
      DATA One  /1.0d0/
      DATA Two  /2.0d0/
C-----------------------------------------------------------------------
C
C     2017/4/6 - RML (LRF=7) Shift NO Longer Allowed - set=0, continue
C
C-----------------------------------------------------------------------
      DO J=1,Npp
      IF (Lpent(J).LT.0 .OR. Lpent(J).GT.1) THEN
      write(3,10) J,Lpent(J)
   10 format(' WARNING..J/Lpent(J)=',2i5,' Expect either 0 or 1')
      ENDIF
      IF (Ishift(J).ne.0) then
      write(3,20) Ishift(J)
      write(*,20) Ishift(J)
   20 format(1x,78('-')/
     1 ' WARNING...Shift=',i5,' LRF=7 Shift Option NOT ALLOWED in',
     2 ' ENDF.'/
     2 '                 ',5x,' Will Ignore and assume Shift = 0.')
      Ishift(J) = 0
      ENDIF
      ENDDO
C
C
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      DO J=Ngr1,Ngr2
C
      IF (dMOD(Spin(J),Half).NE.Zero) THEN
      Write (3,30) J, Spin(J)
      Write (*,30) J, Spin(J)
   30 FORMAT (' ERROR - Quantum Numbers for Spin Group '/
     1        '         Number ', I3,' spin =', F10.5,
     2                                 '(MUST be multiple of 1/2)'/
     3        '         Execution Terminated.'///)
      CALL ENDERROR
      ENDIF
C
      Goj(J) = Zero
      Nent(J) = 0
      Next(J) = 0
      Nchanj = Nchan(J)
      DO N=1,Nchanj
      IF (Ipp(N,J).EQ.2) THEN
      Nent(J) = Nent(J) + 1
      IF (Goj(J).EQ.Zero) THEN
      Goj(J) = (Two*DABS(Spin (      J) )+One)/
     *                   ( (Two*DABS(Spina(Ipp(N,J)))+One)
     *                   * (Two*DABS(Spinb(Ipp(N,J)))+One) )
      ENDIF
      ELSE IF (Ipp(N,J).GT.2) THEN
      Next(J) = Next(J) + 1
      ELSE
      ENDIF
      IF (Lspin(N,J).LT.0) THEN
      write(3,40) N,J,Lspin(N,J)
   40 format(' WARNING..N/J/Lspin(N,J)=',3i5,' Is Negative')
      ENDIF
      IF (Chspin7(N,J).EQ.Zero .OR. Spin(J).EQ.Zero) THEN
      ELSE
      X = One
      I = MOD(Lspin(N,J),2)
      IF (I.NE.0) X = -X
      IF (Chspin7(N,J).LT.Zero) X = -X
      IF (Spin(J).LT.Zero) X = -X
      IF (X.LT.Zero) THEN
      Write (3,50) J,N, Spin(J), Lspin(N,J), Chspin7(N,J)
   50 FORMAT (' *** Parity problem ***', /,
     *               '     Group and channel #', 2I5, /,
     *               '     Spin, L, Chspin7 =', F7.1, I4, F7.1)
      ENDIF
      ENDIF
      Smin = DABS( DABS(Spina(Ipp(N,J))) - DABS(Spinb(Ipp(N,J))) )
      Smax =       DABS(Spina(Ipp(N,J))) + DABS(Spinb(Ipp(N,J)))
      IF (DABS(Chspin7(N,J)).LT.Smin) THEN
      Write (3,80) J, N
      Write (3,60) Chspin7(N,J), Spina(Ipp(N,J)),
     *             Spinb(Ipp(N,J)), Smin
c               Write (*,11500) J, N
c               Write (*,11501) Chspin7(N,J), Spina(Ipp(N,J)),
c    *             Spinb(Ipp(N,J)), Smin
   60 FORMAT (' ******* |Chspin|<|Spina-Spinb|=Smin',3X,4F5.1)
      ELSE IF (DABS(Chspin7(N,J)).GT.Smax) THEN
      Write (3,80) J, N
      Write (3,70) Chspin7(N,J), Spina(Ipp(N,J)),
     *             Spinb(Ipp(N,J)), Smax
   70 FORMAT (' ******* |Chspin|>|Spina+Spinb|=Smax',3X,4F5.1)
      ENDIF
   80 FORMAT (/, ' ****** Error in quantum numbers for',
     *         1x, 'Group Number', I3, ' and Channel Number', I2)
      IF (Lpent(Ipp(N,J)).NE.0) THEN
      Smin = DABS(DABS(Dfloat(Lspin(N,J))-DABS(Chspin7(N,J))))
      Smax = Dfloat(Lspin(N,J)) + DABS(Chspin7(N,J))
      IF (DABS(Spin(J)).LT.Smin) THEN
      Write (3,80) J, N
      Write (3,90) Spin(J), Lspin(N,J), Chspin7(N,J),Smin
   90 FORMAT(' ****** |Spin|<|Lspin-Chspin|=Smin', 3X,4F5.1)
      ELSE IF (DABS(Spin(J)).GT.Smax) THEN
      Write (3,80) J, N
      Write (3,100) Spin(J), Lspin(N,J), Chspin7(N,J),Smin
  100 FORMAT(' ****** |Spin|>|Lspin-Chspin|=Smin', 3X,4F5.1)
      ENDIF
      ENDIF
      ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE Setdim
c=======================================================================
c
c     SAMRML routine with fixed dimensions
c
c     Call once for first Lrf=7 case,
c     to initialize ALL parameters.
c
c=======================================================================
c-----ADD ALL to insure initialization
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
c-----energy limits of resolved range
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
c-----(define parameters read by RECENT for SAMRML)---------
c-----2014/2/20 - ereslow, ereshigh passed through common
      call same1(Awr7,Nro7,Naps7)
c-----(initializememory storage parameters)-------
      RETURN
      END
      SUBROUTINE Fxradi
c=======================================================================
C
C *** Purpose -- fix the following parameters:
C ***            Zke , where k   = Zke * sqrt(E)
C ***            Zkte, where Rho = ka = Zkte * sqrt(E) in penetrability
C ***            Zkfe, where Rho = ka = Zkfe * sqrt(E) in phi
C ***            Zeta7, where eta = Zeta7/sqrt(E) for charged particles
C
C ***            Also fix masses
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
C
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
      DATA Ten /10.0d0/
C
c-----initialize
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      do ii=Ngr1,Ngr2
      do jj=1,MaxNchan
      Zke  (jj,ii) = 0.0d+0
      Zkfe (jj,ii) = 0.0d+0
      Zkte (jj,ii) = 0.0d+0
      Zeta7(jj,ii) = 0.0d+0
      Echan(jj,ii) = 0.0d+0
      enddo
      enddo
C
      Emneut = 1.00866491578d0
C            = mass of neutron in amu
      Hbarrr = 6.582118890d-16
C            = Planck's constant over 2 Pi in (eV s)
      Amuevv = 931.494013d+06
C            = atomic mass unit in eV
      Cspeed = 2.99792458d+08
C            = speed of light in m/s
      Fininv = 1.0d0/137.03599976d0
C            = alpha = fine structure constant
C            = e^2 / (hbar c) dimensionless
      Ff = 1.0d+15
      Twomhb = dSQRT(2.0d0*Emneut*Amuevv)/(Hbarrr*Ff*Cspeed)
      Etac   =           Fininv * Amuevv /(Hbarrr*Ff*Cspeed) * Emneut
C
      DO Ippx=1,Npp
      IF (Ema(Ippx).EQ.Zero) THEN
      Ema(Ippx) = One
      ENDIF
      IF (Emb(Ippx).EQ.Zero) THEN
      Emb(Ippx) = Awr7
      ENDIF
      ENDDO
C
C
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      DO Kgroup=Ngr1,Ngr2
      NchanK = Nchan(Kgroup)
      Factor = Emb(2) + Ema(2)
      Alabcm = Emb(2)/Factor
      Factor = Alabcm/Ema(2)
      DO Ichan=1,NchanK
      Ippx = Ipp(Ichan,Kgroup)
      Aa = Emb(Ippx) + Ema(Ippx)
      Aa = Emb(Ippx)/Aa
      IF (Qqq(Ippx).NE.Zero .AND. Echan(Ichan,Kgroup).EQ.Zero)THEN
      Echan(Ichan,Kgroup) = - Qqq(Ippx) / Alabcm
      ENDIF
      Redmas = Aa * Ema(Ippx)
      Z = Twomhb * dSQRT(Redmas*Factor)
      Zke (Ichan,Kgroup) = Z
      Zkfe(Ichan,Kgroup) = Z*Rdeff(Ichan,Kgroup)*Ten
      Zkte(Ichan,Kgroup) = Z*Rdtru(Ichan,Kgroup)*Ten
      Dcoulomb = Kzb(Ippx)*Kza(Ippx)
      IF(Dcoulomb.gt.0.0d+0)
     1 Zeta7(Ichan,Kgroup) = Etac*Dcoulomb*Redmas/Zke(Ichan,Kgroup)
      ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE Abpart
c=======================================================================
C
C *** Purpose -- Generate Alphar and Alphai => for cross section
C ***            and Upr and Upi = Energy-dependent pieces of Pr & Pi
C ***            Also generate Pr and Pi = partial of R wrt U-parameters
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/RANGER/LOW,LHI
C
      DATA One  /1.0d0/
C
C *** Generate Alphar =  (DEL E) / ( (DEL E)**2 + (Gamgam/2)**2 )
C ***      AND Alphai = Gamgam/2 / ( Ditto )
C
c     DO N=1,Nres
      DO N=LOW,LHI
      Difen (N) = Eres(N) - Su
      Aa = Difen(N)**2 + Gbetpr(3,N)
      Xden  (N) = One/Aa
      Alphar(N) = Difen(N)*Xden(N)
      Alphai(N) = Gbetpr(2,N)*Xden(N)
      ENDDO
      RETURN
      END
      SUBROUTINE Setr
c=======================================================================
C
C *** PURPOSE -- GENERATE Linv = 1/(S-B+IP)
C ***                    Rootp = sqrt(P)
C ***                     Rmat = SUM Beta7*Beta7/((DEL E)-i(GAMGAM/2))
C ***                     Ymat = Linv - Rmat
C ***            Also return Lrmat = 1 if no R-matrix contribution
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
C
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
      DATA Two  /2.0d0/
      DATA Tiny /1.d-8/
C
      Nnntot = Nchann + 1
      DO I=1,Nchann
      Nnntot = Nnntot - 1
      IF (Su.GE.Zero .AND. Su.LE.Echan(Nnntot,Kg)) THEN
      Nchann = Nnntot - 1
      ELSE
      GO TO 10
      ENDIF
      ENDDO
   10 CONTINUE
C
C *** INITIALIZE Rmat = R-MATRIX
C
      KL = 0
      DO K=1,Nchann
      DO L=1,K
      KL = KL + 1
      Rmat(1,KL) = Zero
      Rmat(2,KL) = Zero
      ENDDO
      ENDDO
C
      IF (Maxr.GE.Minr .AND. Minr.GT.0) THEN
      DO Ires=Minr,Maxr
      KL = 0
      DO K=1,Nchann
      DO L=1,K
      KL = KL + 1
      IF (Su.GT.Echan(K,Kg) .AND. Su.GT.Echan(L,Kg) .AND.
     *              Beta7(KL,Ires).NE.Zero) THEN
      Rmat(1,KL) = Rmat(1,KL) + Alphar(Ires)*Beta7(KL,Ires)
      Rmat(2,KL) = Rmat(2,KL) + Alphai(Ires)*Beta7(KL,Ires)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
C *** Check if Rmat is Zero; if so, set Lrmat=1
      KL = 0
      DO K=1,Nchann
      DO L=1,K
      KL = KL + 1
      IF (Rmat(1,KL).NE.Zero) GO TO 20
      IF (Rmat(2,KL).NE.Zero) GO TO 20
      ENDDO
      ENDDO
      Lrmat = 1
   20 CONTINUE
C
C
C *** GENERATE Rootp,PH,(H+Rmat) matrices
C ***          Rootp = SQRT(P)
C ***             PH = 1/(S-B+IP)
C ***           Ymat = (1/(S-B+IP) - Rmat)
C
      KL = 0
      DO K=1,Nchann
      DO L=1,K
      KL = KL + 1
      Ymat(1,KL) = - Rmat(1,KL)
      Ymat(2,KL) = - Rmat(2,KL)
      ENDDO
      ENDDO
C
      II = 0
c
c     Loop over channels
c
      DO 50 I=1,Nchann
      Ipx = Ipp(I,Kg)
      II = II + I
      Rootp (I) =  One
      Elinvr(I) =  Zero
      Elinvi(I) = -One
c-----------------------------------------------------------------------
c
c     Define penetrability.
c
c-----------------------------------------------------------------------
c-----assume below threshold
      P = 1.0d+00
c
c     Only if incident energy is above threshold and Lpent is ON
c
      IF (Su.le.Echan(I,Kg)) go to 40   ! above threshold?
      IF (Lpent(Ipx).NE.1) go to 30     ! penetration turned on
c
c     Penetrability is Calculated
c
      Lsp  = Lspin(I,Kg)
c-----2012/2/19 - added DABS
      Ex   = dSQRT (DABS(Su-Echan(I,Kg)))
      Rho  = Zkte(I,Kg)*Ex
      Rhof = Zkfe(I,Kg)*Ex
C-----coulomb or not?
      IF (Zeta7(I,Kg).EQ.Zero) THEN
c-----no coulomb
      Rho2 = Rho**2  ! SAMRML uses Rho, not Rhof
      CALL FACTS3 (Lsp, Rho2, P)
      ELSE
c-----coulomb
      Eta = Zeta7(I,Kg)/Ex
      CALL COULOMB (Rho, Lsp, Eta, P)
      ENDIF
      if(P.le.1.0d-35) P = 0.0d+0
c-----------------------------------------------------------------------
c
c     no shift
c
c-----------------------------------------------------------------------
      call FACPHI(Lsp, Rhof, Ps) ! SAMRML uses Rhof, not Rho
      Sinphi(I) = Dsin(Ps)
      Cosphi(I) = Dcos(Ps)
      Sinsqr(I) = Sinphi(I)**2             ! sin^2(PS)
      Sin2ph(I) = Two*Sinphi(I)*Cosphi(I)  ! sin(2*PS)
      Rootp(I)  = dSQRT(P)
c
c     Add results to matrices
c
      IF(.NOT.((One-P*Rmat(2,II).EQ.One .OR. P.LT.Tiny))) THEN
c-----Most cases = P is not close to zero
      Hr = 0.0
      Hi = 0.0
      if(P.ne.0.0d+0) Hi = -1.0d+0/P
      Elinvr(I) = Hr
      Elinvi(I) = Hi
      Ymat(1,II) = Hr + Ymat(1,II)
      Ymat(2,II) = Hi + Ymat(2,II)
      ELSE
C-----Rare cases = P is very small but non-Zero
      Ymat(1,II) = P*Ymat(1,II)
      Ymat(2,II) = P*Ymat(2,II) - One
      Rmat(1,II) = P*Rmat(1,II)
      Rmat(2,II) = P*Rmat(2,II)
      IF (Nchann.GT.1) THEN
      IF (I.GT.1) THEN
      DO J=1,I-1
      Ji = (I*(I-1))/2 + J
      Ymat(1,Ji) = Rootp(I)*Ymat(1,Ji)
      Ymat(2,Ji) = Rootp(I)*Ymat(2,Ji)
      Rmat(1,Ji) = Rootp(I)*Rmat(1,Ji)
      Rmat(2,Ji) = Rootp(I)*Rmat(2,Ji)
      ENDDO
      ENDIF
      IF (I.LT.Nchann) THEN
      DO J=I+1,Nchann
      Ji = (J*(J-1))/2 + I
      Ymat(1,Ji) = Rootp(I)*Ymat(1,Ji)
      Ymat(2,Ji) = Rootp(I)*Ymat(2,Ji)
      Rmat(1,Ji) = Rootp(I)*Rmat(1,Ji)
      Rmat(2,Ji) = Rootp(I)*Rmat(2,Ji)
      ENDDO
      ENDIF
c-----end channels > 1 loop
      ENDIF
c-----end of rare case - assume P = 1
      Rootp(I) = One
      Elinvr(I) = Zero
      Elinvi(I) = -One
      ENDIF
      go to 50
C-----Penetrability is NOT Calculated
   30 Ymat(2,II) = Ymat(2,II) - One
      go to 50
c-----Below threshold - assume P = 0
   40 Rootp (I) = Zero
      Elinvr(I) = One
      Elinvi(I) = Zero
c-----end of channel loop
   50 CONTINUE
C
      IF (Lrmat.NE.1) THEN
C ***    Check if one channel is irrelevant; if so, set to unity
      KL = 0
      DO K=1,Nchann
      KL = KL + K
      IF (Ymat(1,KL).EQ.Zero .AND. Ymat(2,KL).EQ.Zero) THEN
      Ymat(1,KL) = One
      ENDIF
      ENDDO
      ENDIF
      RETURN
      END
      SUBROUTINE Crosss
c=======================================================================
C
C *** PURPOSE -- Form the cross sections Sigma(Ksigma) and the
C ***   ( partial derivatives of the cross section with respect to
C ***   the resonance parameters ) = Dsigma(Ksigma,Ipar)
C
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
      common/rmlcross/Sigma(100)
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
C-----note - Fourpi = 4*pi/100
      DATA Fourpi /0.1256637061435917295385057D0/
C
C *** Initialize
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      do ii=1,NppSum
      Sigma(ii) = 0.0d+0
      enddo
c
c     an implicit loop is here faster than an explicit loop.
c
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      Kg = Ngr1
C
C *** DO LOOP OVER GROUPS (IE SPIN-PARITY GROUPS)  -
C ***   GOES TO END OF ROUTINE
C
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      if(Kg.le.1) then
      Maxr   = 0
      else
      Maxr   = Nresg(Kg-1)
      endif
C
C ***    Initialize for this group
c-----2014/2/20 - Extended for Multiple Energy Ranges.
   10 do ii=1,NppSum
      Crss(ii) = 0.0d+0
      enddo
C
      Minr = Maxr + 1
      Maxr = Nresg(Kg)
      Nchann = Nchan(Kg)
      Nentnn = Nent(Kg)
      Nextnn = Next(Kg)
C
      Lrmat = 0
C ***    Set R-Matrix and other necessary arrays
c
      CALL Setr
C
      IF (Lrmat.EQ.1) THEN
C ***       Calculate Xq & Xxxx matrices if trivial
      CALL Zeror
      ELSE
C
C ***       INVERT Ymat; note that Xqr is "Dummy" here
      CALL Yinvrs
C ***       GENERATE XQ & Xxxx matrices
      CALL Setxqx
      ENDIF
C
      Dgoj = Goj(Kg)
C ***    Generate cross section pieces
c
      CALL Sectio
C
C ***    Done calculating contribution for this group; ergo, add to tota
      DO Ip=1,Npp
      Sigma(Ip) = Crss(Ip) + Sigma(Ip)
      ENDDO
C
c----- end of implicit loop
      Kg = Kg + 1
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      if(Kg.le.Ngr2  ) go to 10
C
C *** Normalize properly
      DO Ip=1,Npp
C-----note - Fourpi = 4*pi/100
      Sigma(Ip) = Sigma(Ip)*Fourpi/Su
      ENDDO
C
C
      RETURN
      END
      SUBROUTINE Answer1
c=======================================================================
c
c     define Mt7's to Output.
c
c=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
c-----limit of resolved region for SAMRML
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
      common/usemt/dumrml(100),Quse(11),Mtuse(11)
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
      common/FISSY/LFWX,LFI,MT451,LFWSUM
c
      if(LRF7.gt.0) go to 10  ! LRF7 = # of LRF=7 sections
c
c     NO LRF=7 = Simple case of other formalisms.
c
      MTREACT(1) = 1
      MTREACT(2) = 2
      MTREACT(3) = 102
      QREACT (1) = 0.0D+0
      QREACT (2) = 0.0D+0
      QREACT (3) = 0.0D+0
      NREACT     = 3
      if(LFWSUM.gt.0) then
      MTREACT(4) = 18
      QREACT (4) = 0.0D+0
      NREACT     = 4
      endif
      go to 20
c
c     General LRF = 7 case
c
c-----this NOW handles fission before and after LRF=7 section.
   10 IMFISSY = 0
      if(LFWSUM.gt.0) IMFISSY = 1
c-----does LRF=7 section have fission?
      nofiss = 0
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      do i=1,NppSum
      if(MtuseSum(i).eq.18) nofiss = 1 ! yes
      enddo
c-----assume MT =1, 2 and 102 are first 3 output
      MTREACT(1) = 1
      MTREACT(2) = 2
      MTREACT(3) = 102
      QREACT (1) = 0.0D+0
      QREACT (2) = 0.0D+0
      QREACT (3) = 0.0D+0
      NREACT     = 3
c-----if other range has fission and current does not, add fission.
      if(nofiss.ne.0) IMFISSY = 0
      if(IMFISSY.eq.1.and.nofiss.eq.0) then
      MTREACT(4) = 18
      QREACT (4) = 0.0D+0
      NREACT     = 4
      endif
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      if(NPPSum.le.2) go to 20
      do i=1,NppSum
      if(MtuseSum(i).ne.2.and.MtuseSum(i).ne.102) then
      NREACT = NREACT + 1
      MTREACT(NREACT) = MtuseSum(i)
      QREACT (NREACT) = QuseSum (i)
      endif
      enddo
c
c     All cases - print summary
c
   20 write(*,30) (MTREACT(i),i=1,NREACT)
      write(3,30) (MTREACT(i),i=1,NREACT)
   30 format(1x,78('=')/' Summary of Calculated MT #'/1x,78('-')/(10i5))
      return
      END
      SUBROUTINE Answers
C=======================================================================
c
c     Transform from Samrml to RECENT results
c
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlcross/Sigma(100)     ! Samrml results
      common/rmlfinal/rmlsigma(11)   ! RECENT results
      common/outmt/QREACT(11),MTREACT(11),NEGTAB(11),NREACT,IMFISSY,
     1 LRF7
c
c     WARNING - rmlsigma is NOT initialized = MUST define ALL.
c
c-----Initialize ALL.
      do i=1,11
      rmlsigma(i) = 0.0d+0
      enddo
c----=total = elastic + non-elastic
      rmlsigma(1) = Sigma(1) + Sigma(2)
c-----elastic
      rmlsigma(2) = Sigma(1)
c-----capture (subtract others, if any)
      if(Npp.le.2) then
c-----only capture
      rmlsigma(3) = Sigma(2)
      if(IMFISSY.eq.1) rmlsigma(4) = 0.0
      return
      endif
c-----otherwise, subtract others to define capture
      subrml = 0.0
      do i=3,Npp
      subrml = subrml + Sigma(i)
      enddo
      rmlsigma(3) = Sigma(2) - subrml
c-----define remaining
      KPP = 3
      if(IMFISSY.eq.1) then   ! if needed, define fission = 0
      rmlsigma(4) = 0.0
      KPP = 4
      endif
      do i=3,Npp
      KPP = KPP + 1
      rmlsigma(KPP) = Sigma(i)
      enddo
      RETURN
      END
C ======================================================================
C
C NIST Guide to Available Math Software.
C Source for module Coulfg from package COULOMB.
C Retrieved from TIBER on Fri Apr 3 15:08:25 1998.
C ======================================================================
C *** Converted to double precision.   8 Apr 98.   ROS
C *** Specialized to avoid excess computations in SAMMY.  29 Dec 99. NML
C *** Latest update 19 June 2002 NML
C
C ======================================================================
      SUBROUTINE Zeror
C=======================================================================
C
C     Initialize arrays
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
      DATA Zero /0.0d0/
      N = Nchann
      JK = 0
      DO J=1,N
      DO K=1,J
      JK = JK + 1
      Xxxxr(JK) = Zero
      Xxxxi(JK) = Zero
      ENDDO
      ENDDO
      DO J=1,N
      DO K=1,N
      Xqr(K,J) = Zero
      Xqi(K,J) = Zero
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE Scale3 (V1, V2, V3, Rmat)
C=======================================================================
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Rmat(2,6)
      DATA Zero /0.0d0/
      Aa = 1.D10
      V1 = Zero
      V2 = Zero
      V3 = Zero
C
      IF (DABS(Rmat(1,1)).GE.Aa .OR. DABS(Rmat(2,1)).GE.Aa) THEN
      Bb = DABS(Rmat(1,1))
      Cc = DABS(Rmat(2,1))
      IF (Cc.GT.Bb) Bb = Cc
      V1 = dSQRT(Bb)
      Rmat(1,1) = Rmat(1,1)/Bb
      Rmat(2,1) = Rmat(2,1)/Bb
      Rmat(1,2) = Rmat(1,2)/V1
      Rmat(2,2) = Rmat(2,2)/V1
      Rmat(1,4) = Rmat(1,4)/V1
      Rmat(2,4) = Rmat(2,4)/V1
      ENDIF
C
      IF (DABS(Rmat(1,3)).GE.Aa .OR. DABS(Rmat(2,3)).GE.Aa) THEN
      Bb = DABS(Rmat(1,3))
      Cc = DABS(Rmat(2,3))
      IF (Cc.GT.Bb) Bb = Cc
      V2 = dSQRT(Bb)
      Rmat(1,2) = Rmat(1,2)/V2
      Rmat(2,2) = Rmat(2,2)/V2
      Rmat(1,3) = Rmat(1,3)/Bb
      Rmat(2,3) = Rmat(2,3)/Bb
      Rmat(1,5) = Rmat(1,5)/V2
      Rmat(2,5) = Rmat(2,5)/V2
      ENDIF
C
      IF (DABS(Rmat(1,6)).GE.Aa .OR. DABS(Rmat(2,6)).GE.Aa) THEN
      Bb = DABS(Rmat(1,6))
      Cc = DABS(Rmat(2,6))
      IF (Cc.GT.Bb) Bb = Cc
      V3 = dSQRT(Bb)
      Rmat(1,4) = Rmat(1,4)/V3
      Rmat(2,4) = Rmat(2,4)/V3
      Rmat(1,5) = Rmat(1,5)/V3
      Rmat(2,5) = Rmat(2,5)/V3
      Rmat(1,6) = Rmat(1,6)/Bb
      Rmat(2,6) = Rmat(2,6)/Bb
      ENDIF
C
      RETURN
      END
      SUBROUTINE Unscl3 (V1, V2, V3, Rmat, Rinv)
C=======================================================================
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Rmat(2,6), Rinv(2,6)
      DATA Zero /0.0d0/
C
      IF (V1.GT.Zero) THEN
      Bb = V1**2
      Rmat(1,1) = Rmat(1,1)*Bb
      Rmat(2,1) = Rmat(2,1)*Bb
      Rmat(1,2) = Rmat(1,2)*V1
      Rmat(2,2) = Rmat(2,2)*V1
      Rmat(1,4) = Rmat(1,4)*V1
      Rmat(2,4) = Rmat(2,4)*V1
      Rinv(1,1) = Rinv(1,1)/Bb
      Rinv(2,1) = Rinv(2,1)/Bb
      Rinv(1,2) = Rinv(1,2)/V1
      Rinv(2,2) = Rinv(2,2)/V1
      Rinv(1,4) = Rinv(1,4)/V1
      Rinv(2,4) = Rinv(2,4)/V1
      ENDIF
C
      IF (V2.GT.Zero) THEN
      Bb = V2**2
      Rmat(1,2) = Rmat(1,2)*V2
      Rmat(2,2) = Rmat(2,2)*V2
      Rmat(1,3) = Rmat(1,3)*Bb
      Rmat(2,3) = Rmat(2,3)*Bb
      Rmat(1,5) = Rmat(1,5)*V2
      Rmat(2,5) = Rmat(2,5)*V2
      Rinv(1,2) = Rinv(1,2)/V2
      Rinv(2,2) = Rinv(2,2)/V2
      Rinv(1,3) = Rinv(1,3)/Bb
      Rinv(2,3) = Rinv(2,3)/Bb
      Rinv(1,5) = Rinv(1,5)/V2
      Rinv(2,5) = Rinv(2,5)/V2
      ENDIF
C
      IF (V3.GT.Zero) THEN
      Bb = V3**2
      Rmat(1,4) = Rmat(1,4)*V3
      Rmat(2,4) = Rmat(2,4)*V3
      Rmat(1,5) = Rmat(1,5)*V3
      Rmat(2,5) = Rmat(2,5)*V3
      Rmat(1,6) = Rmat(1,6)*Bb
      Rmat(2,6) = Rmat(2,6)*Bb
      Rinv(1,4) = Rinv(1,4)/V3
      Rinv(2,4) = Rinv(2,4)/V3
      Rinv(1,5) = Rinv(1,5)/V3
      Rinv(2,5) = Rinv(2,5)/V3
      Rinv(1,6) = Rinv(1,6)/Bb
      Rinv(2,6) = Rinv(2,6)/Bb
      ENDIF
C
      RETURN
      END
      SUBROUTINE Xspfa (Bp, N, Kpvt, Info)
C=======================================================================
C
C *** modified October 14, 1993, by NML to use for complex arrays
C
C     Xspfa FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN
C     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING.
C
C     TO SOLVE  A*X = B , FOLLOW Xspfa BY Xspsl.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW Xspfa BY Xspsl.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW Xspfa BY Xspdi.
C     TO COMPUTE  INVERSE(A) , FOLLOW Xspfa BY Xspdi... which I don't ha
C
C     ON ENTRY
C
C        Bp      REAL*8 (2,(N*(N+1)/2))
C                Bp(1,*) is real part, Bp(2,*) is imaginary part.
C                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE
C                COLUMNS OF THE UPPER TRIangle ARE STORED SEQUENTIALLY
C                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 .
C                SEE COMMENTS BELOW FOR DETAILS.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     OUTPUT
C
C        Bp      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHIch
C                WERE USED TO OBTAIN IT STORED IN PACKED FORM.
C                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U)
C                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT
C                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE
C                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL
C                WITH 1 BY 1 AND 2 BY 2 BLOCKS.
C
C        Kpvt    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        Info    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS
C                     Not AN ERROR CONDITION FOR THIS SBROUTINE,
C                     BUT IT DOES INDICATE THAT SSPSL OR SSPDI MAY
C                     DIVIDE BY ZERO IF CALLED.
C
C     PACKED STORAGE
C
C          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER
C          TRIangle OF A SYMMETRIC MATRIX.
C
C                K = 0
C                DO J=1,N
C                   DO I=1,J
C                      K = K + 1
C                      Bp(K) = A(I,J)
C                   ENDDO
C                ENDDO
C
C     SuB-ROUTINES AND FNCTIONS
C       BLAS Xaxpy, Xswap, Ixamax
C       FORTRAN DABS, dMAX1, dSQRT
C
C     INTERNAL VARIABLES
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Bp(2,*), Kpvt(*)
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
C
C
C     INITIALIZE
C
C     Alpha IS USED IN CHOOSING PIVOT BLOCK SIZE.
      Alpha = (One + dSQRT(17.0d0))/8.0d0
C
      Info = 0
C
C     MAIN LOOP ON K, WHIch GOES FROM N TO 1.
C
      K = N
      Ik = (N*(N-1))/2
   10 CONTINUE
C
C        LEAVE THE LOOP IF K=0 OR K=1.
C
C     ...EXIT
      IF (K.EQ.0) GO TO 20
      IF (K.LE.1) THEN
      Kpvt(1) = 1
      IF (Bp(1,1).EQ.Zero .AND. Bp(2,1).EQ.Zero) Info = 1
C     ......EXIT
      GO TO 20
      ENDIF
C
C
C        THIS SECTION OF CODE DETERMINES THE KIND OF
C        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
C        Kstep WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
C        ISWAP WILL BE SET TO 1 IF AN INTERCHANGE IS REQUIRED.
C
      Km1 = K - 1
      Kk = Ik + K
      Absakk = Bp(1,Kk)**2 + Bp(2,Kk)**2
C
C        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN COLUMN K.
C
      Imax = Ixamax (K-1, Bp(1,Ik+1), 1)
      Imk = Ik + Imax
      Colmax = Bp(1,Imk)**2 + Bp(2,Imk)**2
C
      IF (Absakk .GE. Alpha*Colmax) THEN
C
      Kstep = 1
      ISWAP = 0
C
      ELSE
C
C ***       DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN ROW Imax.
      Rowmax = Zero
      Imaxp1 = Imax + 1
      Im = (Imax*(Imax-1))/2
      Imj = Im + 2*Imax
      DO J=Imaxp1,K
      Aa = Bp(1,Imj)**2 + Bp(2,Imj)**2
      Rowmax = dMAX1 (Rowmax, Aa)
      Imj = Imj + J
      ENDDO
      IF (Imax.NE.1) THEN
      Jmax = Ixamax (Imax-1, Bp(1,IM+1), 1)
      Jmim = Jmax + IM
      Aa = Bp(1,Jmim)**2 + Bp(2,Jmim)**2
      Rowmax = dMAX1(Rowmax, Aa)
      ENDIF
      Imim = Imax + Im
      Aa = Bp(1,Imim)**2 + Bp(2,Imim)**2
      IF (Aa .GE. Alpha*Rowmax) THEN
      Kstep = 1
      ISWAP = 1
      ELSE
      IF (Absakk .GE. Alpha*Colmax*(Colmax/Rowmax)) THEN
      Kstep = 1
      ISWAP = 0
      ELSE
      Kstep = 2
      if(Imax .NE. Km1) then
      ISWAP = 1
      else
      ISWAP = 0
      endif
      ENDIF
      ENDIF
C
      ENDIF
C
      IF (dMAX1(Absakk,Colmax) .EQ. Zero) THEN
C
C           COLUMN K IS ZERO.  SET Info AND ITERATE THE LOOP.
      Kpvt(K) = K
      Info = K
      ELSE
C
      IF (Kstep.NE.2) THEN
C              1 X 1 PIVOT BLOCK.
      if(ISWAP.ne.0) then
C                 PERFORM AN INTERCHANGE.
      CALL Xswap (Imax, Bp(1,IM+1), 1, Bp(1,Ik+1), 1)
      Imj = Ik + Imax
      DO Jj=Imax,K
      J = K + Imax - Jj
      Jk = Ik + J
      T = Bp(1,Jk)
      Bp(1,Jk) = Bp(1,Imj)
      Bp(1,Imj) = T
      T = Bp(2,Jk)
      Bp(2,Jk) = Bp(2,Imj)
      Bp(2,Imj) = T
      Imj = Imj - (J - 1)
      ENDDO
      ENDIF
C
C              PERFORM THE ELIMINATION.
      Ij = Ik - (K-1)
      DO Jj=1,Km1
      J = K - Jj
      Jk = Ik + J
      Aa = Bp(1,Kk)**2 + Bp(2,Kk)**2
      Dmulk  = -(Bp(1,Jk)*Bp(1,Kk)+Bp(2,Jk)*Bp(2,Kk))/Aa
      Dmulki =  (Bp(1,Jk)*Bp(2,Kk)-Bp(2,Jk)*Bp(1,Kk))/Aa
      T  = Dmulk
      TI = Dmulki
      CALL Xaxpy (J, T, TI, Bp(1,Ik+1), 1, Bp(1,Ij+1), 1)
      Bp(1,Jk) = Dmulk
      Bp(2,Jk) = Dmulki
      Ij = Ij - (J - 1)
      ENDDO
C
C              SET THE PIVOT ARRAY.
C
      Kpvt(K) = K
      IF (ISWAP.ne.0) Kpvt(K) = Imax
C
      ELSE
C           ELSE IF (Kstep.EQ.2)
C
C              2 X 2 PIVOT BLOCK.
      Km1k = Ik + K - 1
      Ikm1 = Ik - (K-1)
      IF (ISWAP.ne.0) THEN
C
C                 PERFORM AN INTERCHANGE.
      CALL Xswap (Imax, Bp(1,IM+1), 1, Bp(1,Ikm1+1), 1)
      Imj = Ikm1 + Imax
      DO Jj=Imax,Km1
      J = Km1 + Imax - Jj
      Jkm1 = Ikm1 + J
      T = Bp(1,Jkm1)
      Bp(1,Jkm1) = Bp(1,Imj)
      Bp(1,Imj) = T
      T = Bp(2,Jkm1)
      Bp(2,Jkm1) = Bp(2,Imj)
      Bp(2,Imj) = T
      Imj = Imj - (J - 1)
      ENDDO
      T = Bp(1,Km1k)
      Bp(1,Km1k) = Bp(1,Imk)
      Bp(1,Imk) = T
      T = Bp(2,Km1k)
      Bp(2,Km1k) = Bp(2,Imk)
      Bp(2,Imk) = T
      ENDIF
C
C              PERFORM THE ELIMINATION.
      Km2 = K - 2
      IF (Km2.NE.0) THEN
      Aa = Bp(1,Km1k)**2 + Bp(2,Km1k)**2
      Ak  = (Bp(1,Kk)*Bp(1,Km1k)+Bp(2,Kk)*Bp(2,Km1k))/Aa
      Aki = (Bp(2,Kk)*Bp(1,Km1k)-Bp(1,Kk)*Bp(2,Km1k))/Aa
      Km1km1 = Ikm1 + K - 1
      Akm1  = ( Bp(1,Km1km1)*Bp(1,Km1k) +
     *                      Bp(2,Km1km1)*Bp(2,Km1k) ) /Aa
      Akm1i = ( Bp(2,Km1km1)*Bp(1,Km1k) -
     *                      Bp(1,Km1km1)*Bp(2,Km1k) ) /Aa
      Denom  = One - (Ak*Akm1-Aki*Akm1i)
      Denomi =     - (Ak*Akm1i+Aki*Akm1)
      Dd = Denom**2 + Denomi**2
      Ij = Ik - (K-1) - (K-2)
      DO Jj=1,Km2
      J = Km1 - Jj
      Jk = Ik + J
      Bk  = (Bp(1,Jk)*Bp(1,Km1k)+Bp(2,Jk)*Bp(2,Km1k))/Aa
      Bki = (Bp(2,Jk)*Bp(1,Km1k)-Bp(1,Jk)*Bp(2,Km1k))/Aa
      Jkm1 = Ikm1 + J
      Bkm1  = ( Bp(1,Jkm1)*Bp(1,Km1k) +
     *                         Bp(2,Jkm1)*Bp(2,Km1k) ) /Aa
      Bkm1i = ( Bp(2,Jkm1)*Bp(1,Km1k) -
     *                         Bp(1,Jkm1)*Bp(2,Km1k) ) /Aa
      Xx  = Akm1*Bk - Akm1i*Bki - Bkm1
      XxI = Akm1*Bki + Akm1i*Bk - Bkm1i
      Dmulk  = (Xx*Denom+XxI*Denomi)/Dd
      Dmulki = (XxI*Denom-Xx*Denomi)/Dd
      Xx  = Ak*Bkm1 - Aki*Bkm1i - Bk
      XxI = Ak*Bkm1i + Aki*Bkm1 - Bki
      Dmlkm1 = (Xx*Denom+XxI*Denomi)/Dd
      Dmlkmi = (XxI*Denom-Xx*Denomi)/Dd
      T  = Dmulk
      TI = Dmulki
      CALL Xaxpy (J, T, TI, Bp(1,Ik+1), 1, Bp(1,Ij+1), 1)
      T  = Dmlkm1
      TI = Dmlkmi
      CALL Xaxpy (J, T, TI, Bp(1,Ikm1+1), 1,Bp(1,Ij+1),1)
      Bp(1,Jk) = Dmulk
      Bp(2,Jk) = Dmulki
      Bp(1,Jkm1) = Dmlkm1
      Bp(2,Jkm1) = Dmlkmi
      Ij = Ij - (J-1)
      ENDDO
      ENDIF
C
C              SET THE PIVOT ARRAY.
      Kpvt(K) = 1 - K
      IF (ISWAP.ne.0) Kpvt(K) = -Imax
      Kpvt(K-1) = Kpvt(K)
      ENDIF
      ENDIF
      Ik = Ik - (K-1)
      IF (Kstep.EQ.2) Ik = Ik - (K-2)
      K = K - Kstep
      GO TO 10
   20 CONTINUE
      RETURN
      END
      INTEGER*4 FUNCTION Ixamax (N, Sx, Incx)
C=======================================================================
C
C     FINDS THE INDEX OF ELEMENT HAVING Maximum squared value
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Sx(2,*)
C
      Ixamax = 0
      IF (N.LT.1) RETURN
      Ixamax = 1
      IF (N.EQ.1) RETURN
C
      IF (Incx.NE.1) THEN
C
C        CODE FOR INCREMENT Not EQUAL TO 1
      Ix = 1
      Smax = Sx(1,1)**2 + Sx(2,1)**2
      Ix = Ix + Incx
      DO I=2,N
      Aa = Sx(1,Ix)**2 + Sx(2,Ix)**2
      IF (Aa.GT.Smax) THEN
      Ixamax = I
      Smax = Aa
      ENDIF
      Ix = Ix + Incx
      ENDDO
C
      ELSE
C
C        CODE FOR INCREMENT EQUAL TO 1
      Smax = Sx(1,1)**2 + Sx(2,1)**2
      DO I=2,N
      Aa = Sx(1,I)**2 + Sx(2,I)**2
      IF (Aa.GT.Smax) THEN
      Ixamax = I
      Smax = Aa
      ENDIF
      ENDDO
C
      ENDIF
      RETURN
      END
      SUBROUTINE Xaxpy (N, Sa, Sai, Sx, Incx, Sy, Incy)
C=======================================================================
C
C     Complex CONSTANT TIMES A VECTOR PLUS Another VECTOR.
C     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Sx(2,1), Sy(2,1)
C
      IF (N.LE.0) RETURN
      IF (Sa.EQ.0.0D0 .AND. Sai.EQ.0.0D0) RETURN
C
      IF (Incx.NE.1 .OR. Incy.NE.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          Not EQUAL TO 1
      Ix = 1
      Iy = 1
      IF (Incx.LT.0) Ix = (-N+1)*Incx + 1
      IF (Incy.LT.0) Iy = (-N+1)*Incy + 1
      DO I=1,N
      Aa = Sy(1,Iy) + Sa*Sx(1,Ix) - Sai*Sx(2,Ix)
      Sy(2,Iy) = Sy(2,Iy) + Sa*Sx(2,Ix) + Sai*Sx(1,Ix)
      Sy(1,Iy) = Aa
      Ix = Ix + Incx
      Iy = Iy + Incy
      ENDDO
C
C
      ELSE
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C        CLEAN-UP LOOP
      M = MOD (N, 4)
      IF (M.NE.0) THEN
      DO I=1,M
      Aa = Sy(1,I) + Sa*Sx(1,I) - Sai*Sx(2,I)
      Sy(2,I) = Sy(2,I) + Sa*Sx(2,I) + Sai*Sx(1,I)
      Sy(1,I) = Aa
      ENDDO
      IF (N.LT.4) RETURN
      ENDIF
      Mp1 = M + 1
      DO I=Mp1,N,4
      Aa        = Sy(1,I  ) + Sa*Sx(1,I  ) - Sai*Sx(2,I  )
      Sy(2,I  ) = Sy(2,I  ) + Sa*Sx(2,I  ) + Sai*Sx(1,I  )
      Sy(1,I  ) = Aa
      Aa        = Sy(1,I+1) + Sa*Sx(1,I+1) - Sai*Sx(2,I+1)
      Sy(2,I+1) = Sy(2,I+1) + Sa*Sx(2,I+1) + Sai*Sx(1,I+1)
      Sy(1,I+1) = Aa
      Aa        = Sy(1,I+2) + Sa*Sx(1,I+2) - Sai*Sx(2,I+2)
      Sy(2,I+2) = Sy(2,I+2) + Sa*Sx(2,I+2) + Sai*Sx(1,I+2)
      Sy(1,I+2) = Aa
      Aa        = Sy(1,I+3) + Sa*Sx(1,I+3) - Sai*Sx(2,I+3)
      Sy(2,I+3) = Sy(2,I+3) + Sa*Sx(2,I+3) + Sai*Sx(1,I+3)
      Sy(1,I+3) = Aa
      ENDDO
C
      ENDIF
      RETURN
      END
      REAL*8 FUNCTION Xdot (Xdoti, N, Sx, Incx, Sy, Incy)
C=======================================================================
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Sx(2,1), Sy(2,1)
      DATA Zero /0.0d0/
C
      Stemp  = Zero
      Stempi = Zero
      Xdot   = Zero
      Xdoti  = Zero
      IF (N.LE.0) RETURN
      IF (Incx.NE.1 .OR. Incy.NE.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
      Ix = 1
      Iy = 1
      IF (Incx.LT.0) Ix = (-N+1)*Incx + 1
      IF (Incy.LT.0) Iy = (-N+1)*Incy + 1
      DO I=1,N
      Stemp  = Stemp  + Sx(1,Ix)*Sy(1,Iy) - Sx(2,Ix)*Sy(2,Iy)
      Stempi = Stempi + Sx(2,Ix)*Sy(1,Iy) + Sx(1,Ix)*Sy(2,Iy)
      Ix = Ix + Incx
      Iy = Iy + Incy
      ENDDO
      Xdot  = Stemp
      Xdoti = Stempi
C
      ELSE
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C        CLEAN-UP LOOP
      M = MOD (N, 5)
      IF (M.NE.0) THEN
      DO I=1,M
      Stemp  = Stemp  + Sx(1,I)*Sy(1,I) - Sx(2,I)*Sy(2,I)
      Stempi = Stempi + Sx(2,I)*Sy(1,I) + Sx(1,I)*Sy(2,I)
      ENDDO
      IF (N.LT.5) THEN
      Xdot  = Stemp
      Xdoti = Stempi
      RETURN
      ENDIF
      ENDIF
      Mp1 = M + 1
      DO I=Mp1,N,5
      Stemp  = Stemp  + Sx(1,I  )*Sy(1,I  ) - Sx(2,I  )*Sy(2,I  )
     *                      + Sx(1,I+1)*Sy(1,I+1) - Sx(2,I+1)*Sy(2,I+1)
     *                      + Sx(1,I+2)*Sy(1,I+2) - Sx(2,I+2)*Sy(2,I+2)
     *                      + Sx(1,I+3)*Sy(1,I+3) - Sx(2,I+3)*Sy(2,I+3)
     *                      + Sx(1,I+4)*Sy(1,I+4) - Sx(2,I+4)*Sy(2,I+4)
      Stempi = Stempi + Sx(2,I  )*Sy(1,I  ) + Sx(1,I  )*Sy(2,I  )
     *                      + Sx(2,I+1)*Sy(1,I+1) + Sx(1,I+1)*Sy(2,I+1)
     *                      + Sx(2,I+2)*Sy(1,I+2) + Sx(1,I+2)*Sy(2,I+2)
     *                      + Sx(2,I+3)*Sy(1,I+3) + Sx(1,I+3)*Sy(2,I+3)
     *                      + Sx(2,I+4)*Sy(1,I+4) + Sx(1,I+4)*Sy(2,I+4)
      ENDDO
      Xdot  = Stemp
      Xdoti = Stempi
C
      ENDIF
      RETURN
      END
      SUBROUTINE Xswap (N, SX, Incx, SY, Incy)
C=======================================================================
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Sx(2,1),Sy(2,1)
C
      IF (N.LE.0) RETURN
      IF (Incx.NE.1 .OR. Incy.NE.1) THEN
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
      Ix = 1
      Iy = 1
      IF (Incx.LT.0) Ix = (-N+1)*Incx + 1
      IF (Incy.LT.0) Iy = (-N+1)*Incy + 1
      DO I=1,N
      Stemp    = Sx(1,Ix)
      Sx(1,Ix) = Sy(1,Iy)
      Sy(1,Iy) = Stemp
      Stemp    = Sx(2,Ix)
      Sx(2,Ix) = Sy(2,Iy)
      Sy(2,Iy) = Stemp
      Ix = Ix + Incx
      Iy = Iy + Incy
      ENDDO
C
      ELSE
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C        CLEAN-UP LOOP
      M = MOD (N, 3)
      IF (M.NE.0) THEN
      DO I=1,M
      Stemp   = Sx(1,I)
      Sx(1,I) = Sy(1,I)
      Sy(1,I) = Stemp
      Stemp   = Sx(2,I)
      Sx(2,I) = Sy(2,I)
      Sy(2,I) = Stemp
      ENDDO
      IF (N.LT.3) RETURN
      ENDIF
      Mp1 = M + 1
      DO I=Mp1,N,3
      Stemp     = Sx(1,I  )
      Sx(1,I  ) = Sy(1,I  )
      Sy(1,I  ) = Stemp
      Stemp     = Sx(2,I  )
      Sx(2,I  ) = Sy(2,I  )
      Sy(2,I  ) = Stemp
      Stemp     = Sx(1,I+1)
      Sx(1,I+1) = Sy(1,I+1)
      Sy(1,I+1) = Stemp
      Stemp     = Sx(2,I+1)
      Sx(2,I+1) = Sy(2,I+1)
      Sy(2,I+1) = Stemp
      Stemp     = Sx(1,I+2)
      Sx(1,I+2) = Sy(1,I+2)
      Sy(1,I+2) = Stemp
      Stemp     = Sx(2,I+2)
      Sx(2,I+2) = Sy(2,I+2)
      Sy(2,I+2) = Stemp
      ENDDO
C
      ENDIF
      RETURN
      END
      SUBROUTINE Xspsl (Bp, N, Kpvt, B)
C=======================================================================
C
C     Xspsl SOLVES THE complex SYMMETRIC SYSTEM
C                 A * X = B
C     USING THE FACTORS COMPUTED BY Xspfa.
C
C     ON ENTRY
C
C        Bp      REAL*8 (2,N*(N+1)/2)
C                THE OUTPUT FROM SSPFA.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        Kpvt    INTEGER(N)
C                THE PIVOT VECTOR FROM SSPFA.
C
C        B       REAL(2,N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A Division BY ZERO MAY OCCUR IF  Xspco  HAS SET RCOND .EQ. 0.0
C        OR  Xspfa  HAS SET Info .NE. 0  .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL Xspfa (Bp, N, Kpvt, Info)
C           IF (Info .NE. 0) GO TO ...
C           DO J=1,P
C              CALL Xspsl (Bp, N, Kpvt, C(1,J))
C           ENDDO
C
C     SBROUTINES AND FNCTIONS
C
C     BLAS Xaxpy, Xdot
C     FORTRAN IABS
C
C     INTERNAL VARIABLES.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION Bp(2,1), B(2,1), Kpvt(*)
      DATA One /1.0d0/
C
C     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
C     D INVERSE TO B.
C
      K = N
      Ik = (N*(N-1))/2
   10 CONTINUE
      IF (K.NE.0) THEN
      Kk = Ik + K
C
      IF (Kpvt(K).GE.0) THEN
C
C           1 X 1 PIVOT BLOCK.
      IF (K.NE.1) THEN
      Kp = Kpvt(K)
      IF (Kp.NE.K) THEN
C                 INTERCHANGE.
      Temp    = B(1,K)
      B(1,K)  = B(1,Kp)
      B(1,Kp) = Temp
      Temp    = B(2,K)
      B(2,K)  = B(2,Kp)
      B(2,Kp) = Temp
      ENDIF
C              APPLY THE TRANSFORMATION.
      CALL Xaxpy (K-1, B(1,K), B(2,K), Bp(1,Ik+1), 1, B(1,1),1)
      ENDIF
C           APPLY D INVERSE.
      Aa = Bp(1,Kk)**2 + Bp(2,Kk)**2
      Ab     = (B(1,K)*Bp(1,Kk)+B(2,K)*Bp(2,Kk))/Aa
      B(2,K) = (B(2,K)*Bp(1,Kk)-B(1,K)*Bp(2,Kk))/Aa
      B(1,K) = Ab
      K = K - 1
      Ik = Ik - K
C
      ELSE
C        ELSE IF (Kpvt(K).LT.0)
C
C           2 X 2 PIVOT BLOCK.
      Ikm1 = Ik - (K-1)
      IF (K.NE.2) THEN
      Kp = IABS(Kpvt(K))
      IF (Kp.NE.K-1) THEN
C                 INTERCHANGE.
      Temp     = B(1,K-1)
      B(1,K-1) = B(1,Kp)
      B(1,Kp)  = Temp
      Temp     = B(2,K-1)
      B(2,K-1) = B(2,Kp)
      B(2,Kp)  = Temp
      ENDIF
C              APPLY THE TRANSFORMATION.
      CALL Xaxpy (K-2, B(1,K  ), B(2,K  ), Bp(1,Ik  +1),
     *                       1, B(1,1), 1)
      CALL Xaxpy (K-2, B(1,K-1), B(2,K-1), Bp(1,Ikm1+1),
     *                       1, B(1,1), 1)
      ENDIF
C           APPLY D INVERSE.
      Km1k = Ik + K - 1
      Kk = Ik + K
      Aa = Bp(1,Km1k)**2 + Bp(2,Km1k)**2
      Ak  = (Bp(1,Kk)*Bp(1,Km1k)+Bp(2,Kk)*Bp(2,Km1k))/Aa
      Aki = (Bp(2,Kk)*Bp(1,Km1k)-Bp(1,Kk)*Bp(2,Km1k))/Aa
      Km1km1 = Ikm1 + K - 1
      Akm1  = (Bp(1,Km1km1)*Bp(1,Km1k)+Bp(2,Km1km1)*Bp(2,Km1k))/Aa
      Akm1i = (Bp(2,Km1km1)*Bp(1,Km1k)-Bp(1,Km1km1)*Bp(2,Km1k))/Aa
      Bk  = (B(1,K)*Bp(1,Km1k)+B(2,K)*Bp(2,Km1k))/Aa
      Bki = (B(2,K)*Bp(1,Km1k)-B(1,K)*Bp(2,Km1k))/Aa
      Bkm1  = (B(1,K-1)*Bp(1,Km1k)+B(2,K-1)*Bp(2,Km1k))/Aa
      Bkm1i = (B(2,K-1)*Bp(1,Km1k)-B(1,K-1)*Bp(2,Km1k))/Aa
      Denom  = Ak*Akm1 - Aki*Akm1i - One
      Denomi = Ak*Akm1i + Aki*Akm1
      Dd = Denom**2 + Denomi**2
      Ab  = Akm1*Bk - Akm1i*Bki - Bkm1
      Abi = Akm1i*Bk + Akm1*Bki - Bkm1i
      B(1,K) = (Ab*Denom+Abi*Denomi)/Dd
      B(2,K) = (Abi*Denom-Ab*Denomi)/Dd
      Ab  = Ak*Bkm1 - Aki*Bkm1i - Bk
      Abi = Aki*Bkm1 + Ak*Bkm1i - Bki
      B(1,K-1) = (Ab*Denom+Abi*Denomi)/Dd
      B(2,K-1) = (Abi*Denom-Ab*Denomi)/Dd
      K = K - 2
      Ik = Ik - (K+1) - K
C
      ENDIF
C
      GO TO 10
      ENDIF
C
C     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
      K = 1
      Ik = 0
   20 CONTINUE
      IF (K.LE.N) THEN
      IF (Kpvt(K).GE.0) THEN
C
C           1 X 1 PIVOT BLOCK.
      IF (K.NE.1) THEN
C              APPLY THE TRANSFORMATION.
      B(1,K) = B(1,K) + Xdot (Xdoti, K-1, Bp(1,Ik+1),
     *                                          1, B(1,1), 1)
      B(2,K) = B(2,K) + Xdoti
      Kp = Kpvt(K)
      IF (Kp.NE.K) THEN
C                 INTERCHANGE.
      Temp    = B(1,K )
      B(1,K ) = B(1,Kp)
      B(1,Kp) = Temp
      Temp    = B(2,K )
      B(2,K ) = B(2,Kp)
      B(2,Kp) = Temp
      ENDIF
      ENDIF
      Ik = Ik + K
      K = K + 1
C
      ELSE
C
C           2 X 2 PIVOT BLOCK.
      IF (K.NE.1) THEN
C              APPLY THE TRANSFORMATION.
      B(1,K  ) = B(1,K  ) + Xdot (Xdoti, K-1, Bp(1,Ik  +1),
     *                                           1, B(1,1), 1)
      B(2,K  ) = B(2,K  ) + Xdoti
      Ikp1 = Ik + K
      B(1,K+1) = B(1,K+1) + Xdot (Xdoti, K-1, Bp(1,Ikp1+1),
     *                                           1, B(1,1),1)
      B(2,K+1) = B(2,K+1) + Xdoti
      Kp = IABS(Kpvt(K))
      IF (Kp.NE.K) THEN
C                 INTERCHANGE.
      Temp    = B(1,K )
      B(1,K ) = B(1,Kp)
      B(1,Kp) = Temp
      Temp    = B(2,K )
      B(2,K ) = B(2,Kp)
      B(2,Kp) = Temp
      ENDIF
      ENDIF
      Ik = Ik + K + K + 1
      K  = K  + 2
C
      ENDIF
      GO TO 20
      ENDIF
      RETURN
      END
      SUBROUTINE Setxqx
C=======================================================================
C
C *** purpose -- form XQ & XXXX matrices, where
C ***            XQ   = Yinv * Rmat       and
C ***            XXXX = P/L + sqrt(P)/L         (1/L-R)**-1          sqr
C ***                 =       sqrt(P)/(S-B+IP) * Yinv       * Rmat * sqr
C ***                 =       sqrt(P)/L        * XQ                * sqr
C
C ***            Note that the matrix W defined in SAMMY manual is given
C ***                 by W(c,c') = delta(c,c') + 2i XXXX(c,c')
C ***                 as in Eq. (III.D.4) in SAMMY manual R3
C
C ***         ie W    = I + 2i XXXX
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
C
      do i=1,Nchann
      do j=1,Nchann
      Xqr(i,j) = 0.0d+0
      Xqi(i,j) = 0.0d+0
      enddo
      enddo
C
C *** Xqr(k,i) = (L**-1-R)**-1 * R ... note asymmetry
      DO I=1,Nchann
      DO J=1,Nchann
      Ij = Ijkl(J,I)
      DO K=1,Nchann
      Jk = Ijkl(K,J)
      Xqr(K,I) = Xqr(K,I) + Yinv(1,Ij)*Rmat(1,Jk) -
     *                               Yinv(2,Ij)*Rmat(2,Jk)
      Xqi(K,I) = Xqi(K,I) + Yinv(1,Ij)*Rmat(2,Jk) +
     *                               Yinv(2,Ij)*Rmat(1,Jk)
      ENDDO
      ENDDO
      ENDDO
C
C *** Xxxx = sqrt(P)/L  * xq * sqrt(P) ... symmetric
      IJ = 0
      DO I=1,Nchann
      Plr = Rootp(I)*Elinvr(I)
      Pli = Rootp(I)*Elinvi(I)
      DO J=1,I
      Ij = Ij + 1
      Xxxxr(Ij) = Rootp(J)* (Xqr(J,I)*Plr-Xqi(J,I)*Pli)
      Xxxxi(Ij) = Rootp(J)* (Xqi(J,I)*Plr+Xqr(J,I)*Pli)
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE Sectio
C=======================================================================
C
C *** purpose -- generate pieces of cross sections (except for "4 pi / E
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
      DATA Two  /2.0d0/
C
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      DO Jj=1,NppSum
      Crss(Jj) = Zero
      ENDDO
C
C     Entrance channel, Ipp=2
C ***    elastic  Crss(1) = g*0.25* sum(entrance chs c,c')
C ***                                     times |(1-U(c,c'))| **2 / Zz
C ***                     = g* [ sin(phi)**2 * (1-2XXXXi)
C ***                            - sin(2phi)*XXXXr
C ***                            + (XXXXr**2 + XXXXi**2) ] / Zz
      Ii = 0
      Ij = 0
      DO I=1,Nentnn
      Zz = Zke(I,Kg)**2
      Ii = Ii + I
      Termn = Sinsqr(I)*( One - Two * Xxxxi(Ii) )
     *                    - Sin2ph(I)*Xxxxr(Ii)
      Termn = Termn / Zz
      DO J=1,I
      Ij = Ij + 1
      Ar = ( Xxxxr(Ij)**2 + Xxxxi(Ij)**2 )/Zz
      IF (I.NE.J) Ar = Ar + Ar
      Termn = Termn + Ar
      ENDDO
      Crss(1) = Termn + Crss(1)
      ENDDO
      Crss(1) = Crss(1)*Dgoj
C *** End of Ipp=2 term (elastic)
C
C *** Ipp=1 term, sort-of
C ***    absorption = g*0.25 * sum(inc c)
C ***                  [ 1 -  sum(inc c') |U(c,c')| **2 ] / Zz
C ***               = - g* (XXXXr**2 + XXXXi**2) / Zz
      Ii = 0
      Ij = 0
      DO I=1,Nentnn
      Ii = Ii + I
      Zz = Zke(I,Kg)**2
      Terma =  Xxxxi(Ii) / Zz
      DO J=1,I
      Ij = Ij + 1
      Ar = (- Xxxxr(Ij)**2 - Xxxxi(Ij)**2) / Zz
      IF (I.NE.J) Ar = Ar + Ar
      Terma = Terma + Ar
      ENDDO
      Crss(2) = Terma + Crss(2)
      ENDDO
      Crss(2) = Crss(2)*Dgoj
C *** End of absorption term
C
C *** All other channels, classed by particle-pair number
C ***    reaction ch c'= g*0.25 * sum(inc c) |U(c,c')|**2 / Zz
C ***                  = g* (XxxxR**2 + XxxxI**2) / Zz
      DO Jj=1,Nextnn
      J = Jj + Nentnn
      IF (J.LE.Nchann) THEN
      Ip = Ipp(J,Kg)
      DO I=1,Nentnn
      Zz = Zke(I,Kg)**2
      Ij = (J*(J-1))/2 + I
Cq             Ij = Ijkl(I,J) but I < J always
      Crss(Ip) = Crss(Ip) + (Xxxxr(Ij)**2+Xxxxi(Ij)**2)/Zz
      ENDDO
      ENDIF
      ENDDO
      IF (Npp.GT.2) THEN
      DO Ip=3,Npp
      Crss(Ip) = Crss(Ip)*Dgoj
      ENDDO
      ENDIF
C *** End of other reaction types
C
C
      RETURN
      END
      INTEGER*4 FUNCTION Ijkl (M,N)
C=======================================================================
C
C     Define Ijkl INDEX
C
C=======================================================================
      INCLUDE 'implicit.h'
      IF (M.LE.N) THEN
      Ijkl = (N*(N-1))/2 + M
      ELSE
      Ijkl = (M*(M-1))/2 + N
      ENDIF
      RETURN
      END
      SUBROUTINE Betset
C=======================================================================
C
C *** PURPOSE -- GENERATE Betapr, Gbetpr, Beta7, U
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
      common/comresol/ereslow,ereshigh,QuseSum(11),KresSum,KgrSum,
     1 Ngr1,Ngr2,MtuseSum(11),NppSum
C
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
      DATA Half /0.5d0/
C
C
C *** Convert to Betapr parameters
C
      IF (Nres.GT.0) THEN
C
c-----2014/2/20 - Extended for Multiple Energy Ranges.
      if(Ngr1.le.1) then
      Minres = 1
      else
      Minres = Nresg(Ngr1-1) + 1
      endif
      DO Ig=Ngr1,Ngr2
C
      Mmaxc = Nchan(Ig)
      DO J=1,Mmaxc
      Lsp = Lspin(J,Ig)
      DO Ires=Minres,Nresg(Ig)
      Betapr(J,Ires) = Zero
      Rho = Zero
      P = One
      IF (Gamma(J,Ires).EQ.Zero) THEN
      ELSE
      IF (Lpent(Ipp(J,Ig)).NE.0) THEN
c-----2012/2/19 - Removed DABS
      Ex = DABS (      Eres(Ires) -Echan(J,Ig) )
      IF (Ex.NE.Zero) THEN
      Ex = dSQRT(Ex)
      Rho = Zkte(J,Ig)*Ex
      IF (Zeta7(J,Ig).EQ.Zero) THEN
c-----no coulomb
      Rho2 = Rho**2
      CALL FACTS3 (Lsp, Rho2, P)
      if(P.le.0.0d+0) P = 1.0d+0
      ELSE
c-----coulomb
      Eta = Zeta7(J,Ig)/Ex
      CALL COULOMB (Rho, Lsp, Eta, P)
      if(P.le.0.0d+0) P = 1.0d+0
      ENDIF
      ENDIF
      ENDIF
      Betapr(J,Ires) = dSQRT(Half*DABS(Gamma(J,Ires))/P)
      IF (Gamma(J,Ires).LT.Zero) Betapr(J,Ires) =
     *                                      - Betapr(J,Ires)
      ENDIF
      ENDDO
      ENDDO
C
C ***       Generate Beta7 parameters
      DO Ires=Minres,Nresg(Ig)
      KL = 0
      DO K=1,Mmaxc
      DO L=1,K
      KL = KL + 1
      Beta7(KL,Ires) = Betapr(L,Ires)*Betapr(K,Ires)
      ENDDO
      ENDDO
      ENDDO
C
      DO Ires=Minres,Nresg(Ig)
      Gbetpr(2,Ires) = Half*DABS(Gamgam(Ires))
      Gbetpr(1,Ires) = dSQRT(Gbetpr(2,Ires))
      IF (Gamgam(Ires).LT.Zero) Gbetpr(1,Ires)= -Gbetpr(1,Ires)
      Gbetpr(3,Ires) = Gbetpr(2,Ires)**2
      ENDDO
C
      Minres = Nresg(Ig) + 1
      ENDDO
C
      ENDIF
C
C
      RETURN
      END
      SUBROUTINE Yinvrs
C=======================================================================
C
C *** PURPOSE -- INVERT Ymat TO GIVE Yinv
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
c
      if(Nchann.gt.3) go to 40
      go to (10,20,30),Nchann
c----- 1 Channel
   10 CALL Onech
      return
c----- 2 Channels
   20 CALL Twoch
      return
c----- 3 Channels
   30 CALL Three
      return
c----- >3 Channels
   40 CALL Yfour
      return
      END
      SUBROUTINE Onech
C=======================================================================
C
C *** PURPOSE -- Invert Ymat to give Yinv, for the one-channel case
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
C
      IF (Ymat(1,1)+Ymat(2,1).EQ.Ymat(1,1)) THEN
      Yinv(1,1) = One/Ymat(1,1)
      Yinv(2,1) = -(Ymat(2,1)/Ymat(1,1))/Ymat(1,1)
C
      ELSE IF (Ymat(1,1)+Ymat(2,1).EQ.Ymat(2,1)) THEN
      Yinv(1,1) = (Ymat(1,1)/Ymat(2,1))/Ymat(2,1)
      Yinv(2,1) = -One/Ymat(2,1)
C
      ELSE IF (Ymat(1,1).EQ.Zero) THEN
      Yinv(1,1) = Zero
      Yinv(2,1) = -One/Ymat(2,1)
C
      ELSE
      Aa = Ymat(1,1)**2 + Ymat(2,1)**2
      Yinv(1,1) = Ymat(1,1)/Aa
      Yinv(2,1) = -Ymat(2,1)/Aa
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE Twoch
C=======================================================================
C
C *** PURPOSE -- Invert Ymat to give Yinv, for the two-channel case
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
C
      IF (Ymat(1,1).NE.Zero .OR. Ymat(1,2).NE.Zero .OR.
     *    Ymat(1,3).NE.Zero) THEN
C
C
C ***    First step must be to scale such that one large number
C ***       does not cause this to blow up...
      V = Zero
      K = 0
      DO I=1,3
      IF (DABS(Ymat(1,I)).GT.V) K = I
      IF (DABS(Ymat(1,I)).GT.V) V = DABS(Ymat(1,I))
      IF (DABS(Ymat(2,I)).GT.V) K = I
      IF (DABS(Ymat(2,I)).GT.V) V = DABS(Ymat(2,I))
      ENDDO
      IF (V.GT.Zero) THEN
      IF (K.EQ.2) THEN
      DO I=1,3
      Ymat(1,I) = Ymat(1,I)/V
      Ymat(2,I) = Ymat(2,I)/V
      ENDDO
      ELSE
      Ymat(1,K) = Ymat(1,K)/V
      Ymat(2,K) = Ymat(2,K)/V
      Ymat(1,2) = Ymat(1,2)/dSQRT(V)
      Ymat(2,2) = Ymat(2,2)/dSQRT(V)
      ENDIF
      ENDIF
C
      Bbr = Ymat(1,1)*Ymat(1,3) -
     *      Ymat(2,1)*Ymat(2,3) - Ymat(1,2)**2 +
     *      Ymat(2,2)**2
      Bbi = Ymat(1,1)*Ymat(2,3) +
     *      Ymat(2,1)*Ymat(1,3) - 2.*Ymat(1,2)*
     *      Ymat(2,2)
C
      IF (Bbr+Bbi.NE.Bbi) THEN
      IF (Bbr+Bbi.NE.Bbr) THEN
      Aa = One/(Bbr*Bbr+Bbi*Bbi)
      Aar = Bbr*Aa
      Aai = -Bbi*Aa
      ELSE
      Aar = One/Bbr
      Aai = -(Bbi/Bbr)/Bbr
      ENDIF
      ELSE
      Aar = (Bbr/Bbi)/Bbi
      Aai = -One/Bbi
      ENDIF
C
      Yinv(1,1) =  Aar*Ymat(1,3) - Aai*Ymat(2,3)
      Yinv(2,1) =  Aar*Ymat(2,3) + Aai*Ymat(1,3)
      Yinv(1,2) = -Aar*Ymat(1,2) + Aai*Ymat(2,2)
      Yinv(2,2) = -Aar*Ymat(2,2) - Aai*Ymat(1,2)
      Yinv(1,3) =  Aar*Ymat(1,1) - Aai*Ymat(2,1)
      Yinv(2,3) =  Aar*Ymat(2,1) + Aai*Ymat(1,1)
C
      IF (V.NE.Zero) THEN
      IF (K.EQ.2) THEN
      DO I=1,3
      Yinv(1,I) = Yinv(1,I)/V
      Yinv(2,I) = Yinv(2,I)/V
      ENDDO
      ELSE
      Yinv(1,K) = Yinv(1,K)/V
      Yinv(2,K) = Yinv(2,K)/V
      Yinv(1,2) = Yinv(1,2)/dSQRT(V)
      Yinv(2,2) = Yinv(2,2)/dSQRT(V)
      ENDIF
      ENDIF
C
C
      ELSE IF (Ymat(2,2).NE.Zero) THEN
C ***    Here when real part of Ymat is Zero everywhere, imaginary part
C ***       dense
      Yinv(1,1) = Zero
      Yinv(1,2) = Zero
      Yinv(1,3) = Zero
      V = Ymat(2,1)*Ymat(2,3) - Ymat(2,2)**2
      Yinv(2,1) = -Ymat(2,3)/V
      Yinv(2,2) =  Ymat(2,2)/V
      Yinv(2,3) = -Ymat(2,1)/V
C
C
      ELSE
C ***    Here when only real part of Ymat is Zero everywhere; imaginary
C ***       part is Zero off-diagonal and non-Zero on diagonal
      Yinv(1,1) =  Zero
      Yinv(1,2) =  Zero
      Yinv(1,3) =  Zero
      Yinv(2,1) = -One /(Ymat(2,1))
      Yinv(2,2) =  Zero
      Yinv(2,3) = -One /(Ymat(2,3))
C
C
      ENDIF
      RETURN
      END
      SUBROUTINE Three
C=======================================================================
C
C *** Purpose -- Invert Ymat to give Yinv, for three-channel case
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
      DATA Two  /2.0d0/
C
      CALL Scale3 (V1, V2, V3, Ymat(1,1))
C
      Izz = 0
      DO I=1,6
      IF (Ymat(1,I).NE.Zero) Izz = 1
      ENDDO
C
      IF (Ymat(2,2).NE.Zero) Izz = 1
      IF (Ymat(2,4).NE.Zero) Izz = 1
      IF (Ymat(2,5).NE.Zero) Izz = 1
C
C
      IF (Izz.EQ.0) THEN
C
C ***    Here if only imaginary parts of diagonal terms are nonZero
      DO I=1,6
      Yinv(1,I) = Zero
      Yinv(2,I) = Zero
      ENDDO
      Yinv(2,1) = -One/Ymat(2,1)
      Yinv(2,3) = -One/Ymat(2,3)
      Yinv(2,6) = -One/Ymat(2,6)
C
      ELSE
C
C ***    Here real parts of diagonal terms are nonZero
      Fc1r = Ymat(1,3)*Ymat(1,6) -
     *          Ymat(2,3)*Ymat(2,6) - Ymat(1,5)**2 +
     *          Ymat(2,5)**2
      Fc1i = Ymat(2,3)*Ymat(1,6) +
     *          Ymat(1,3)*Ymat(2,6) - Two*Ymat(1,5)*
     *          Ymat(2,5)
      Fc2r = Ymat(1,2)*Ymat(1,6) -
     *          Ymat(2,2)*Ymat(2,6) -
     *          Ymat(1,4)*Ymat(1,5) +
     *          Ymat(2,4)*Ymat(2,5)
      Fc2i = Ymat(2,2)*Ymat(1,6) +
     *          Ymat(1,2)*Ymat(2,6) -
     *          Ymat(1,4)*Ymat(2,5) -
     *          Ymat(2,4)*Ymat(1,5)
      Fc3r = Ymat(1,2)*Ymat(1,5) -
     *          Ymat(2,2)*Ymat(2,5) -
     *          Ymat(1,4)*Ymat(1,3) +
     *          Ymat(2,4)*Ymat(2,3)
      Fc3i = Ymat(2,2)*Ymat(1,5) +
     *          Ymat(1,2)*Ymat(2,5) -
     *          Ymat(2,4)*Ymat(1,3) -
     *          Ymat(1,4)*Ymat(2,3)
C
      Aar = Ymat(1,1)*Fc1r -  Ymat(1,2)*Fc2r +
     *         Ymat(1,4)*Fc3r - (Ymat(2,1)*Fc1i
     *           -Ymat(2,2)*Fc2i+Ymat(2,4)*Fc3i)
      Aai = Ymat(1,1)*Fc1i -  Ymat(1,2)*Fc2i +
     *         Ymat(1,4)*Fc3i + (Ymat(2,1)*Fc1r
     *           -Ymat(2,2)*Fc2r+Ymat(2,4)*Fc3r)
C
      IF (Aar+Aai.NE.Aai) THEN
      IF (Aar+Aai.NE.Aar) THEN
      Aa = Aar**2 + Aai**2
      Dcr = Aar/Aa
      Dci = -Aai/Aa
      ELSE
      Dcr = One/Aar
      Dci = -(Aai/Aar)/Aar
      ENDIF
      ELSE
      Dcr = (Aar/Aai)/Aai
      Dci = -One/Aai
      ENDIF
C
      Yinv(1,1) = Fc1r*Dcr - Fc1i*Dci
      Yinv(2,1) = Fc1r*Dci + Fc1i*Dcr
      Yinv(1,2) = -Fc2r*Dcr + Fc2i*Dci
      Yinv(2,2) = -Fc2r*Dci - Fc2i*Dcr
      Yinv(1,3) = (Ymat(1,1)*Ymat(1,6)
     *                    -Ymat(1,4)**2
     *                    -Ymat(2,1)*Ymat(2,6)
     *                    +Ymat(2,4)**2               )  *  Dcr
     *              - (    Ymat(2,1)*Ymat(1,6)
     *                -Two*Ymat(1,4)*Ymat(2,4)
     *                    +Ymat(1,1)*Ymat(2,6)   )  *  Dci
      Yinv(2,3) = (Ymat(1,1)*Ymat(1,6)
     *                    -Ymat(1,4)**2
     *                    -Ymat(2,1)*Ymat(2,6)
     *                    +Ymat(2,4)**2             )   *   Dci
     *              + (    Ymat(2,1)*Ymat(1,6)
     *                -Two*Ymat(1,4)*Ymat(2,4)
     *                    +Ymat(1,1)*Ymat(2,6) )   *   Dcr
      Yinv(1,4) = Fc3r*Dcr - Fc3i*Dci
      Yinv(2,4) = Fc3r*Dci + Fc3i*Dcr
      Yinv(1,5) = -(Ymat(1,1)*Ymat(1,5)
     *                     -Ymat(1,2)*Ymat(1,4)
     *                     -Ymat(2,1)*Ymat(2,5)
     *                     +Ymat(2,2)*Ymat(2,4)  ) * Dcr
     *               + (    Ymat(1,1)*Ymat(2,5)
     *                     -Ymat(1,2)*Ymat(2,4)
     *                     +Ymat(2,1)*Ymat(1,5)
     *                     -Ymat(2,2)*Ymat(1,4)  ) * Dci
      Yinv(2,5) = -(Ymat(1,1)*Ymat(1,5)
     *                     -Ymat(1,2)*Ymat(1,4)
     *                     -Ymat(2,1)*Ymat(2,5)
     *                     +Ymat(2,2)*Ymat(2,4)  ) * Dci
     *                 - (  Ymat(1,1)*Ymat(2,5)
     *                     -Ymat(1,2)*Ymat(2,4)
     *                     +Ymat(2,1)*Ymat(1,5)
     *                     -Ymat(2,2)*Ymat(1,4)  ) * Dcr
      Yinv(1,6) = (Ymat(1,1)*Ymat(1,3)
     *                    -Ymat(1,2)**2
     *                    -Ymat(2,1)*Ymat(2,3)
     *                    +Ymat(2,2)**2              ) * Dcr
     *             - (     Ymat(2,1)*Ymat(1,3)
     *                -Two*Ymat(1,2)*Ymat(2,2)
     *                    +Ymat(1,1)*Ymat(2,3)  ) * Dci
      Yinv(2,6) = (Ymat(1,1)*Ymat(1,3)
     *                    -Ymat(1,2)**2
     *                    -Ymat(2,1)*Ymat(2,3)
     *                    +Ymat(2,2)**2              ) * Dci
     *             + (     Ymat(2,1)*Ymat(1,3)
     *                -Two*Ymat(1,2)*Ymat(2,2)
     *                    +Ymat(1,1)*Ymat(2,3)  ) * Dcr
      CALL Unscl3 (V1, V2, V3, Ymat(1,1), Yinv(1,1))
      ENDIF
      RETURN
      END
      SUBROUTINE Yfour
C=======================================================================
C
C *** PURPOSE -- calculate Ymat**-1, for any number of channels.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'rmlcom.h'
      common/rmlflags/Dgoj,Nchann,Lrmat,Nentnn,Nextnn
      DIMENSION Dummy(2,MaxNchan),Kpvt(100)
      equivalence (Xqr(1,1),Dummy(1,1))
      DATA Zero /0.0d0/
      DATA One  /1.0d0/
C
      CALL Xspfa (Ymat(1,1), Nchann, Kpvt(1), Info)
      IF (Info.NE.0) Write (3,10) Info
   10 FORMAT (' Problem in Xspfa with Info=', I5)
      Kj = 0
      DO K=1,Nchann
      DO J=1,Nchann
      Dummy(1,J) = Zero
      Dummy(2,J) = Zero
      ENDDO
      Dummy(1,K) = One
      CALL Xspsl (Ymat(1,1), Nchann, Kpvt(1), Dummy(1,1))
      DO J=1,K
      Kj = Kj + 1
      Yinv(1,KJ) = Dummy(1,J)
      Yinv(2,KJ) = Dummy(2,J)
      ENDDO
      ENDDO
      RETURN
      END
      subroutine sigma1t(infile,outfile,Tres)
C=======================================================================
C
C     PROGRAM SIGMA1
C     ==============
C     VERSION 73-1 (MARCH 1973)
C     VERSION 76-1 (FEBRUARY 1976)
C     VERSION 76-2 (OCTOBER 1976)
C     VERSION 77-1 (JANUARY 1977)
C     VERSION 78-1 (JULY 1978)
C     VERSION 79-1 (JULY 1979)    CDC-7600 AND CRAY-1 VERSION.
C     VERSION 80-1 (MAY 1980)     IBM, CDC AND CRAY VERSION
C     VERSION 80-2 (DECEMBER 1980)IMPROVED BASED ON USER COMMENTS.
C     VERSION 81-1 (MARCH 1981)   DOUBLE PRECISION IBM VERSION
C     VERSION 81-2 (AUGUST 1981)  IMPROVED IBM SPEED AND STABILITY
C     VERSION 82-1 (JANUARY 1982) IMPROVED COMPUTER COMPATIBILITY
C     VERSION 83-1 (JANUARY 1983)*MAJOR RE-DESIGN.
C                                *PAGE SIZE INCREASED - 1002 TO 2004.
C                                *ELIMINATED COMPUTER DEPENDENT CODING.
C                                *NEW, MORE COMPATIBLE I/O UNIT NUMBER.
C                                *ADDED STANDARD ALLOWABLE ERROR OPTION
C                                 (CURRENTLY 0.1 PER-CENT).
C                                *UNRESOLVED RESONANCE REGION COPIED.
C                                *1/V EXTENSION OF CROSS SECTIONS
C                                 OUTSIDE OF TABULATED ENERGY RANGE AND
C                                 INTO UNRESOLVED ENERGY RANGE.
C     VERSION 83-2 (OCTOBER 1983)*IMPROVED BASED ON USER COMMENTS.
C     VERSION 84-1 (APRIL 1984)  *IMPROVED NUMERICAL STABILITY.
C                                *PARTIAL EVALUATION TREATMENT.
C     VERSION 85-1 (APRIL 1985)  *ITERATE TO CONVERGENCE (USING THE SAME
C                                 ENERGY GRID FOR HOT CROSS SECTION AS
C                                 COLD CROSS SECTIONS WAS FOUND TO BE
C                                 INACCURATE).
C                                *NEW FASTER HIGH ENERGY BROADENING.
C                                *UPDATED FOR ENDF/B-VI FORMATS.
C                                *SPECIAL I/O ROUTINES TO GUARANTEE
C                                 ACCURACY OF ENERGY.
C                                *DOUBLE PRECISION TREATMENT OF ENERGY
C                                 (REQUIRED FOR NARROW RESONANCES).
C     VERSION 85-2 (AUGUST 1985) *FORTRAN-77/H VERSION
C     VERSION 86-1 (JANUARY 1986)*ENERGY DEPENDENT SCATTERING RADIUS
C     VERSION 88-1 (JULY 1988)   *OPTION...INTERNALLY DEFINE ALL I/O
C                                 FILE NAMES (SEE, SUBROUTINE FILEIO
C                                 FOR DETAILS).
C                                *IMPROVED BASED ON USER COMMENTS.
C     VERSION 89-1 (JANUARY 1989)*PSYCHOANALYZED BY PROGRAM FREUD TO
C                                 INSURE PROGRAM WILL NOT DO ANYTHING
C                                 CRAZY.
C                                *UPDATED TO USE NEW PROGRAM CONVERT
C                                 KEYWORDS.
C                                *ADDED LIVERMORE CIVIC COMPILER
C                                 CONVENTIONS.
C     VERSION 90-1 (JUNE 1990)   *UPDATED BASED ON USER COMMENTS
C                                *ADDED FORTRAN SAVE OPTION
C                                *NEW MORE CONSISTENT ENERGY OUTPUT
C                                 ROUTINES
C     VERSION 91-1 (JULY 1991)   *WARNING...INPUT PARAMETER FORMAT
C                                 HAS BEEN CHANGED - SEE BELOW FOR
C                                 DETAILS.
C                                *ADDED CHARGED PARTICLE PROJECTILES
C                                *OUTPUT ENERGY RANGE IS ALWAYS AT
C                                 LEAST AS LARGE AS INPUT ENERGY RANGE.
C                                *NO 1/V EXTENSION OF CROSS SECTIONS
C                                 FROM UNRESOLVED ENERGY RANGE.
C     VERSION 92-1 (JANUARY 1992)*INSURE MINIMUM AND MAXIMUM CROSS
C                                 SECTIONS ARE ALWAYS KEPT (NOT THINNED)
C                                *MT=19 (FIRST CHANCE FISSION) TREATED
C                                 THE SAME AS FISSION.
C                                *VARIABLE MINIMUM CROSS SECTION OF
C                                 INTEREST - TO ALLOW SMALL CROSS
C                                 SECTIONS NEAR THRESHOLDS TO BE
C                                 TREATED PROPERLY.
C                                *ALL ENERGIES INTERNALLY ROUNDED PRIOR
C                                 TO CALCULATIONS.
C                                *COMPLETELY CONSISTENT I/O AND ROUNDING
C                                 ROUTINES - TO MINIMIZE COMPUTER
C                                 DEPENDENCE.
C     VERSION 92-2 (JULY 1992)   *CORRECTED BUG ASSOCIATED WITH
C                                 THRESHOLD REACTIONS.
C                                *UNRESOLVED REGION COPIED WITHOUT
C                                 THINNING (IT SHOULD BE EXACTLY THE
C                                 SAME AT ALL TEMPERATURES).
C                                *NO THINNING OF REACTIONS (MT) THAT
C                                 WERE NOT BROADENED.
C     VERSION 93-1 (APRIL 1993)  *INCREASED PAGE SIZE FROM 2004
C                                 TO 24000 ENERGY PONTS.
C     VERSION 94-1 (JANUARY 1994)*VARIABLE ENDF/B DATA FILENAMES
C                                 TO ALLOW ACCESS TO FILE STRUCTURES
C                                 (WARNING - INPUT PARAMETER FORMAT
C                                 HAS BEEN CHANGED)
C                                *CLOSE ALL FILES BEFORE TERMINATING
C                                 (SEE, SUBROUTINE ENDIT)
C     VERSION 96-1 (JANUARY 1996) *COMPLETE RE-WRITE
C                                 *IMPROVED COMPUTER INDEPENDENCE
C                                 *ALL DOUBLE PRECISION
C                                 *ON SCREEN OUTPUT
C                                 *UNIFORM TREATMENT OF ENDF/B I/O
C                                 *IMPROVED OUTPUT PRECISION
C                                 *DEFINED SCRATCH FILE NAMES
C                                 *ALWAYS INCLUDE THERMAL VALUE
C     VERSION 97-1 (APRIL 1997)   *OPTIONALLY SET NEGATIVE CROSS
C                                  SECTIONS = 0 ON INPUT AND
C                                  OUTPUT.
C                                 *INCREASED PAGE SIZE FROM 24000
C                                  TO 60000 ENERGY POINTS.
C     VERSION 99-1 (MARCH 1999)   *CORRECTED CHARACTER TO FLOATING
C                                  POINT READ FOR MORE DIGITS
C                                 *UPDATED TEST FOR ENDF/B FORMAT
C                                  VERSION BASED ON RECENT FORMAT CHANGE
C                                 *TREAT LOW ENERGY INITIAL CROSS
C                                  SECTIONS AS LOG-LOG INTERPOLABLE
C                                 *CONSTANT (RATHER THAN 1/V) EXTENSION
C                                  TO HIGHER ENERGY.
C                                 *UPDATED CONSTANTS BASED ON CSEWG
C                                  SUBCOMMITTEE RECOMMENDATIONS
C                                 *GENERAL IMPROVEMENTS BASED ON
C                                  USER FEEDBACK
C     VERSION 99-2 (JUNE 1999)    *EXTENDED RANGE OF INTEGRALS FROM 4
C                                  TO 5 UNITS ON EACH SIDE OF ENERGY
C                                  POINT TO ALLOW FOR LARGER VARIATION
C                                  IN THE LOCAL CROSS SECTION
C                                 *ASSUME ENDF/B-VI, NOT V, IF MISSING
C                                  MF=1, MT-451.
C     VERSION 99-3 (OCTOBER 1999))*IMPROVED ERFC FUNCTION DEFINITION.
C                                  I THANK BOB MACFARLANE (LANL) FOR
C                                  SUPPLYING A MORE ACCURATE ERFC
C                                  FUNCTION.
C     VERS. 2000-1 (FEBRUARY 2000)*CORRECTED LOW ENERGY INTERPOLATION
C                                  FOR NON-POSITIVE CROSS SECTIONS
C                                 *GENERAL IMPROVEMENTS BASED ON
C                                  USER FEEDBACK
C     VERS. 2002-1 (MAY 2002)     *OPTIONAL INPUT PARAMETERS
C     VERS. 2004-1 (JAN. 2004)    *OPTIONALLY IGNORE UNRESOLVED REGION
C                                 *CORRECTED PROBLEM AT THE RESOLVED/
C                                  UNRESOLVED ENERGY BOUNDARY.
C                                 *CORRECTED HIGH ENERGY CONSTANT CROSS
C                                  SECTION EXTENSION.
C                                 *TIGHTER CRITERIA FOR INITIAL ENERGY
C                                  POINT SPACING
C                                 *TEMPERATURE DEPENDENT ENERGY POINT
c                                  SPACING.
C                                 *ADDED NEW REICH-MOORE (LRF=7) TO
C                                  FILE2 TO ALLOW COPY TO FIND ANY
C                                  FOLLOWING UNRESOLVED PARAMETERS
C     VERS. 2005-1 (JUNE 2005)    *CORRECTED ERROR IN EHOT3 EQUIVALENCE
C                                  TO EHOT - THIS ONLY EFFECTS VERY BIG
C                                  OUTPUT FILES.
C     VERS. 2007-1 (JAN. 2007)    *CHECKED AGAINST ALL ENDF/B-VII.
C                                 *INCREASED PAGE SIZE FROM 60,000
C                                  TO 360,000 ENERGY POINTS.
C     VERS. 2008-1 (APRIL 2008)   *1/2 INITIAL ENERGY POINT SPACING
C                                 *72 CHARACTER FILE NAMES.
C     VERS. 2010-1 (Apr. 2010)    *ASSUME LOW ENERGY LOG-LOG VARIATION
C                                  UP TO 1/A (eV) FOR ALL BUT TOTAL AND
C                                  ELASTIC.
C                                 *CHANGED DEFAULT UNCERTAINTY TO 0.01%
C                                  FROM 0.1%
C                                 *ALLOW MULTIPLE, ADJACENT UNRESOLVED
C                                  RESONANCE REGIONS = COMBINE INTO ONE
C                                  LARGER ENERGY RANGE TO COPY.
C                                 *DO NOT BROADEN SECTIONS THAT START
C                                  ABOVE 1 MILLION KT - PREVIOUSLY IT
C                                  WAS ASSUMED TOTAL, ELASTIC, CAPTURE
C                                  AND FISSION, AND LARGE SECTIONS (OVER
C                                  10,000 ENERGY POINTS) WOULD BROADEN.
C     VERS. 2012-1 (Aug. 2012)    *CHANGE COPY CRITERIA TO HANDLE NEW
C                                  (N,N') DATA = THRESHOLD MAY BE VERY
C                                  HIGH (OLD CRITERIA) BUT INCLUDES MANY
C                                  TABULATED ENERGY POINTS (NEW ADDED
C                                  CRITERIA).
C                                 *ADDED STOP IF INCIDENT PARTICLE DATA
C                                  CANNOT BE DOPPLER BROADENED, E.G.,
C                                  PHOTON INCIDENT.
C                                 *Added CODENAME
C                                 *32 and 64 bit Compatible
C                                 *Added ERROR stop
C     VERS. 2013-1 (Nov. 2013)    *Added NO broadening above 10 MeV -
C                                  this is to handle newer evaluations
C                                  that extend to higher energies and
C                                  may do "strange" things to stop one
C                                  MT and then include it as part of
C                                  a sum at higher energies, e.g. this
C                                  change will copy ALL points above
C                                  10 MeV, thus avoiding problems near
C                                  transistion energies at 20. 30, etc.
C                                  MeV or higher energies.
C     VERS. 2015-1 (Jan. 2015)    *Replaced ALL 3 way IF Statements.
C                                 *Replaced ALL LOGICAL by INTEGER.
C                                 *Extended OUT9.
C     VERS. 2017-1 (May  2017)    *For MF=2 only use MT=151 = Defines
C                                  Unresolved Resonance Region (URR).
C                                  Ignore - NJOY created MT=152 and 153.
C                                 *Increased page size to 1,2000,000.
C                                 *All floating input parameters changed
C                                  to character input + IN9 conversion.
C                                 *Added NRO = energy dependent scatter
C                                  radius to copying FILE2 parameters
C                                  to define unresolved energy range.
C
C     OWNED, MAINTAINED AND DISTRIBUTED BY
C     ------------------------------------
C     THE NUCLEAR DATA SECTION
C     INTERNATIONAL ATOMIC ENERGY AGENCY
C     P.O. BOX 100
C     A-1400, VIENNA, AUSTRIA
C     EUROPE
C
C     ORIGINALLY WRITTEN BY
C     ------------------------------------
C     Dermott E. Cullen
C
C     PRESENT CONTACT INFORMATION
C     ---------------------------
C     Dermott E. Cullen
C     1466 Hudson Way
C     Livermore, CA 94550
C     U.S.A.
C     Telephone  925-443-1911
C     E. Mail    RedCullen1@Comcast.net
C     Website    RedCullen1.nedt/HOMEPAGE.NEW
C
C     Acknowledgement 2004
C     --------------------
C     Currently almost all improvements to this code are based upon
C     feedback from code users who report problems. This feedback
C     benefits ALL users of this code, and ALL users are encouraged
C     to report problems.
C
C     Improvements on the 2004 version of this code based on user
C     feedback including,
C     1) Bret Beck  - reported a problem at the resolved/unresolved
C                     energy boundary.
C     2) S. Ganesan - reported a problem for small temperature changes.
C
C     AUTHORS MESSAGE
C     ---------------
C     THE REPORT DESCRIBED ABOVE IS THE LATEST PUBLISHED DOCUMENTATION
C     FOR THIS PROGRAM. HOWEVER, THE COMMENTS BELOW SHOULD BE CONSIDERED
C     THE LATEST DOCUMENTATION INCLUDING ALL RECENT IMPROVEMENTS. PLEASE
C     READ ALL OF THESE COMMENTS BEFORE IMPLEMENTATION, PARTICULARLY
C     THE COMMENTS CONCERNING MACHINE DEPENDENT CODING.
C
C     AT THE PRESENT TIME WE ARE ATTEMPTING TO DEVELOP A SET OF COMPUTER
C     INDEPENDENT PROGRAMS THAT CAN EASILY BE IMPLEMENTED ON ANY ONE
C     OF A WIDE VARIETY OF COMPUTERS. IN ORDER TO ASSIST IN THIS PROJECT
C     IT WOULD BE APPECIATED IF YOU WOULD NOTIFY THE AUTHOR OF ANY
C     COMPILER DIAGNOSTICS, OPERATING PROBLEMS OR SUGGESTIONS ON HOW TO
C     IMPROVE THIS PROGRAM. HOPEFULLY, IN THIS WAY FUTURE VERSIONS OF
C     THIS PROGRAM WILL BE COMPLETELY COMPATIBLE FOR USE ON YOUR
C     COMPUTER.
C
C     PURPOSE
C     -------
C     THIS PROGRAM IS DESIGNED TO DOPPLER BROADEN NEUTRON INDUCED
C     CROSS SECTIONS. EACH SECTION OF CROSS SECTIONS (FILE 3) IS READ
C     FROM THE ENDF/B FORMAT. THE DATA IS DOPPLER BROADENED, THINNED
C     AND OUTPUT IN THE ENDF/B FORMAT.
C
C     IN THE FOLLOWING DISCUSSION FOR SIMPLICITY THE ENDF/B TERMINOLOGY
C     ---ENDF/B TAPE---WILL BE USED. IN FACT THE ACTUAL MEDIUM MAY BE
C     TAPE, CARDS, DISK OR ANY OTHER MEDIUM.
C
C     ENDF/B FORMAT
C     -------------
C     THIS PROGRAM ONLY USES THE ENDF/B BCD OR CARD IMAGE FORMAT (AS
C     OPPOSED TO THE BINARY FORMAT) AND CAN HANDLE DATA IN ANY VERSION
C     OF THE ENDF/B FORMAT (I.E., ENDF/B-I, II, III, IV OR V FORMAT).
C
C     IT IS ASSUMED THAT THE DATA IS CORRECTLY CODED IN THE ENDF/B
C     FORMAT AND NO ERROR CHECKING IS PERFORMED. IN PARTICULAR IT IS
C     ASSUMED THAT THE MAT, MF AND MT ON EACH CARD IS CORRECT. SEQUENCE
C     NUMBERS (COLUMNS 76-80) ARE IGNORED ON INPUT, BUT WILL BE
C     CORRECTLY OUTPUT ON ALL CARDS. THE FORMAT OF SECTION MF=1, MT=451
C     AND ALL SECTIONS OF MF=3 MUST BE CORRECT. THE PROGRAM COPIES ALL
C     OTHER SECTION OF DATA AS HOLLERITH AND AS SUCH IS INSENSITIVE TO
C     THE CORRECTNESS OR INCORRECTNESS OF ALL OTHER SECTIONS.
C
C     ALL CROSS SECTIONS THAT ARE USED BY THIS PROGRAM MUST BE TABULATED
C     AND LINEARLY INTERPOLABLE IN ENERGY AND CROSS SECTION (ENDF/B
C     INTERPOLATION LAW 2). FILE 3 CROSS SECTIONS MAY BE MADE LINEARLY
C     INTERPOLABLE BY USING PROGRAM LINEAR (UCRL-50400, VOL.17, PART A).
C     FILE 2 RESONANCE PARAMETERS MAY BE USED TO RECONSTRUCT ENERGY
C     DEPENDENT CROSS SECTIONS AND ADD IN FILE 3 BACKGROUND CROSS
C     SECTIONS TO DEFINE LINEARLY INTERPOLABLE CROSS SECTIONS BY USING
C     PROGRAM RECENT (UCRL-50400, VOL. 17, PART C). IF THIS PROGRAM
C     FINDS THAT THE FILE 3 CROSS SECTIONS ARE NOT LINEARLY INTERPOLABLE
C     THIS PROGRAM WILL TERMINATE EXECUTION.
C
C     UNRESOLVED RESONANCE REGION
C     ---------------------------
C     IN THE UNRESOLVED RESONANCE REGION IT IS NOT POSSIBLE TO EXACTLY
C     DEFINE THE ENERGY DEPENDENCE OF THE CROSS SECTIONS. THE AVERAGE
C     WIDTHS AND SPACINGS GIVEN IN ENDF/B ARE ONLY ADEQUATE TO DEFINE
C     AVERAGE VALUES OF THE CROSS SECTIONS. THEREFORE ALL CROSS SECTIONS
C     IN THE ENDF/B FORMAT FOR THE UNRESOLVED REGION ARE REALLY AVERAGE
C     VALUES WHICH CANNOT BE DOPPLER BROADENED USING THE SIGMA1 METHOD
C     (WHICH REQUIRES TABULATED, LINEARLY INTERPOLABLE, ENERGY DEPENDENT
C     CROSS SECTIONS.
C
C     THEREFORE,
C     (1) ALL TABULATED POINTS WITHIN THE UNRESOLVED RESONANCE REGION
C     WILL BE COPIED, WITHOUT MODIFICATION OR BROADENING. ADOPTION OF
C     THIS CONVENTION WILL ALLOW SUBSEQUENT PROGRAMS TO PROPERLY DEFINE
C     SELF-SHIELDED, DOPPLER BROADENED CROSS SECTIONS IN THE UNRESOLVED
C     RESONANCE REGION.
C     (2) CROSS SECTIONS WILL BE EXTENDED AS 1/V ABOVE THE UPPER ENERGY
C     LIMIT OF THE RESOLVED RESONANCE REGION AND BELOW THE LOWER ENERGY
C     LIMIT OF THE CONTINUUUM REGION (I.E. INTO THE UNRESOLVED
C     RESONANCE REGION). THIS CONVENTION WILL GUARANTEE A SMOOTH
C     BEHAVIOR CLOSE TO THE UNRESOLVED RESONANCE REGION BOUNDARIES.
C
C     OUTPUT FORMAT
C     -------------
C     IN THIS VERSION OF SIGMA1 ALL FILE 3 ENERGIES WILL BE OUTPUT IN
C     F (INSTEAD OF E) FORMAT IN ORDER TO ALLOW ENERGIES TO BE WRITTEN
C     WITH UP TO 9 DIGITS OF ACCURACY. IN PREVIOUS VERSIONS THIS WAS AN
C     OUTPUT OPTION. HOWEVER USE OF THIS OPTION TO COMPARE THE RESULTS
C     OF ENERGIES WRITTEN IN THE NORMAL ENDF/B CONVENTION OF 6 DIGITS
C     TO THE 9 DIGIT OUTPUT FROM THIS PROGRAM DEMONSTRATED THAT FAILURE
C     TO USE THE 9 DIGIT OUTPUT CAN LEAD TO LARGE ERRORS IN THE DATA
C     JUST DUE TO TRANSLATION OF THE ENERGIES TO THE ENDF/B FORMAT.
C
C     CONTENTS OF OUTPUT
C     ------------------
C     ENTIRE EVALUATIONS ARE OUTPUT, NOT JUST THE BROADENED FILE 3
C     CROSS SECTIONS, E.G. ANGULAR AND ENERGY DISTRIBUTIONS ARE ALSO
C     INCLUDED.
C
C     DOCUMENTATION
C     -------------
C     THE FACT THAT THIS PROGRAM HAS OPERATED ON THE DATA IS DOCUMENTED
C     BY THE ADDITION OF THREE COMMENTS CARDS AT THE END OF EACH
C     HOLLERITH SECTION IN THE FORM
C
C     ***************** PROGRAM SIGMA1 (2017-1) ***************
C     DATA DOPPLER BROADENED TO 300.0   KELVIN AND
C     DATA THINNED TO WITHIN AN ACCURACY OF  0.1 PER-CENT
C
C     THE ORDER OF ALL SIMILAR COMMENTS (FROM LINEAR,RECENT AND GROUPY)
C     REPRESENTS A COMPLETE HISTORY OF ALL OPERATIONS PERFORMED ON
C     THE DATA.
C
C     THESE COMMENT CARDS ARE ONLY ADDED TO EXISTING HOLLERITH SECTIONS,
C     I.E., THIS PROGRAM WILL NOT CREATE A HOLLERITH SECTION. THE FORMAT
C     OF THE HOLLERITH SECTION IN ENDF/B-V DIFFERS FROM THE THAT OF
C     EARLIER VERSIONS OF ENDF/B. BY READING AN EXISTING MF=1, MT=451
C     IT IS POSSIBLE FOR THIS PROGRAM TO DETERMINE WHICH VERSION OF
C     THE ENDF/B FORMAT THE DATA IS IN. WITHOUT HAVING A SECTION OF
C     MF=1, MT=451 PRESENT IT IS IMPOSSIBLE FOR THIS PROGRAM TO
C     DETERMINE WHICH VERSION OF THE ENDF/B FORMAT THE DATA IS IN, AND
C     AS SUCH IT IS IMPOSSIBLE FOR THE PROGRAM TO DETERMINE WHAT FORMAT
C     SHOULD BE USED TO CREATE A HOLLERITH SECTION.
C
C     REACTION INDEX
C     --------------
C     THIS PROGRAM DOES NOT USE THE REACTION INDEX WHICH IS GIVEN IN
C     SECTION MF=1, MT=451 OF EACH EVALUATION.
C
C     THIS PROGRAM DOES NOT UPDATE THE REACTION INDEX IN MF=1, MT=451.
C     THIS CONVENTION HAS BEEN ADOPTED BECAUSE MOST USERS DO NOT
C     REQUIRE A CORRECT REACTION INDEX FOR THEIR APPLICATIONS AND IT WAS
C     NOT CONSIDERED WORTHWHILE TO INCLUDE THE OVERHEAD OF CONSTRUCTING
C     A CORRECT REACTION INDEX IN THIS PROGRAM. HOWEVER, IF YOU REQUIRE
C     A REACTION INDEX FOR YOUR APPLICATIONS, AFTER RUNNING THIS PROGRAM
C     YOU MAY USE PROGRAM DICTIN TO CREATE A CORRECT REACTION INDEX.
C
C     SECTION SIZE
C     ------------
C     SINCE THIS PROGRAM USES A LOGICAL PAGING SYSTEM THERE IS NO LIMIT
C     TO THE NUMBER OF POINTS IN ANY SECTION, E.G., THE TOTAL CROSS
C     SECTION MAY BE REPRESENTED BY 200,000 DATA POINTS.
C
C     SELECTION OF DATA
C     -----------------
C     THE PROGRAM SELECTS MATERIALS TO BE BROADENED BASED EITHER ON
C     MAT (ENDF/B MAT NO.) OR ZA. THE PROGRAM ALLOWS UP TO 100 MAT OR
C     ZA RANGES TO BE SPECIFIED. THE PROGRAM WILL ASSUME THAT THE
C     ENDF/B TAPE IS IN EITHER MAT OR ZA ORDER, WHICHEVER CRITERIA IS
C     USED TO SELECT MATERIALS, AND WILL TERMINATE WHEN A MAT OR ZA
C     IS FOUND THAT IS ABOVE THE RANGE OF ALL REQUESTS.
C
C     ENERGY GRID OF BROADENED DATA
C     -----------------------------
C     THE ENERGY GRID FOR THE DOPPLER BROADENED CROSS SECTIONS IS
C     SELECTED TO INSURE THAT THE BROADENED DATA IS LINEAR-LINEAR
C     INTERPOLABLE. AS SUCH THE ENERGY GRID FOR THE BROADENED DATA
C     MAY NOT BE THE SAME AS THE ENERGY GRID FOR THE ORIGINAL
C     UNBROADENED DATA. GENERALLY AFTER BROADENING THERE WILL BE
C     FEWER DATA POINTS IN THE RESONANCE REGION, BUT AT LOW ENERGY
C     THERE MAY BE MORE POINTS, DUE TO THE 1/V LOW ENERGY EFFECT
C     CREATED BY DOPPLER BROADENING.
C
C     EFFECTIVE TEMERATURE INCREASE
C     -----------------------------
C     IF THE ORIGINAL DATA IS NOT AT ZERO KELVIN THE PROGRAM WILL
C     BROADEN THE DATA BY THE EFFECTIVE TEMPERATURE DIFFENCE TO THE
C     FINAL TEMPERATURE. IF THE DATA IS ALREADY AT A TEMPERATURE THAT
C     IS HIGHER THAN THE FINAL TEMPERATURE DOPPLER BROADENING IS
C     NATURALLY NOT PERFORMED AND THE TEMPERATURE IN THE SECTION IS LEFT
C     AT ITS ORIGINAL VALUE.
C
C     MULTIPLE FINAL TEMPERATURES
C     ---------------------------
C     THE PRESENT VERSION ONLY DOPPLER BROADENS TO ONE FINAL TEMPERATURE
C     (IF THERE IS SUFFICIENT INTEREST EXPRESSED BY USERS FUTURE
C     VERSION MAY BROADEN TO MULTIPLE TEMPERATURES. PLEASE
C     CONTACT THE AUTHOR IF YOU ARE INTERESTED IN A MULTIPLE
C     TEMPERATURE OPTION).
C
C     PROGRAM OPERATION
C     -----------------
C     EACH SECTION OF FILE 3 DATA IS CONSIDERED SEPERATELY. THE DATA
C     IS READ AND DOPPLER BROADENED A PAGE AT A TIME (ONE PAGE IS
C     60000 DATA POINTS). UP TO THREE PAGES OF DATA MAY BE IN THE CORE
C     AT ANY GIVEN TIME, THE PAGE BEING BROADENED, THE PAGE BELOW IT
C     IN ENERGY AND THE PAGE ABOVE IT IN ENERGY. AFTER A PAGE HAS BEEN
C     BROADENED IT IS THINNED, IF THE ENTIRE SECTION CONTAINS ONLY
C     ONE PAGE OR LESS, IT WILL STILL BE CORE RESIDENT AND WILL BE
C     WRITTEN DIRECTLY FROM CORE TO THE OUTPUT TAPE. IF THE BROADENED,
C     THINNED SECTION IS LARGER THAN A PAGE, AFTER A PAGE HAS BEEN
C     BROADENED AND THINNED IT IS WRITTEN TO A SCRATCH FILE. AFTER THE
C     ENTIRE SECTION HAS BEEN BROADENED AND THINNED THE DATA IS READ
C     FROM SCRATCH TO CORE, ONE PAGE AT A TIME, THE OUTPUT TO THE OUTPUT
C     TAPE.
C
C     ALLOWABLE ERROR
C     ---------------
C     AFTER DOPPLER BROADENING THE CROSS SECTION IN THE RESONANCE REGION
C     WILL GENERALLY BE MUCH SMOOTHER THAN THE UNBROADENED DATA AND CAN
C     BE REPRESENTED TO THE SAME ACCURACY BY A SMALLER NUMBER OF ENERGY
C     POINTS. THEREFORE AFTER DOPPLER BROADENING THE DATA CAN BE THINNED
C     WITH ESSENTIALLY NO LOSE OF INFORMATION.
C
C     THE ALLOWABLE ERROR MAY BE ENERGY INDEPENDENT (CONSTANT) OR ENERGY
C     DEPENDENT. THE ALLOWABLE ERROR IS DESCRIBED BY A TABULATED
C     FUNCTION OF UP TO 20 (ENERGY,ERROR) PAIRS AND LINEAR INTERPOLATION
C     BETWEEN TABULATED POINTS. IF ONLY ONE TABULATED POINT IS GIVEN THE
C     ERROR WILL BE CONSIDERED CONSTANT OVER THE ENTIRE ENERGY RANGE.
C     WITH THIS ENERGY DEPENDENT ERROR ONE MAY OPTIMIZE THE OUTPUT FOR
C     ANY GIVEN APPLICATION BY USING A SMALL ERROR IN THE ENERGY RANGE
C     OF INTEREST AND A LESS STRINGENT ERROR IN OTHER ENERGY RANGES.
C
C     INPUT FILES
C     -----------
C     UNIT  DESCRIPTION
C     ----  -----------
C        2  INPUT CARDS (BCD - 80 CHARACTERS/RECORD)
C       10  ORIGINAL ENDF/B DATA (BCD - 80 CHARACTERS/RECORD)
C
C     OUTPUT FILES
C     ------------
C     UNIT  DESCRIPTION
C     ----  -----------
C        3  OUTPUT REPORT (BCD - 120 CHARACTERS/RECORD)
C       11  FINAL ENDF/B DATA (BCD - 80 CHARACTERS/RECORD)
C
C     SCRATCH FILES
C     -------------
C     UNIT  DESCRIPTION
C     ----  -----------
C       12  SCRATCH FILE FOR BROADENED DATA
C           (BINARY - 180000 WORDS/RECORD - DOUBLE PRECISION/
C                      42000 WORDS/RECORD - SINLGE PRECISION)
C
C     OPTIONAL STANDARD FILE NAMES (SEE SUBROUTINE FILEIO)
C     ----------------------------------------------------
C     UNIT  FILE NAME
C     ----  ----------
C       2   SIGMA1.INP
C       3   SIGMA1.LST
C      10   ENDFB.IN
C      11   ENDFB.OUT
C      12   (SCRATCH)
C
C     INPUT CARDS
C     -----------
C     CARD  COLS.  DESCRIPTION
C     ----  -----  -----------
C        1   1-11  SELECTION CRITERIA (0=MAT, 1=ZA)
C           12-22  MONITOR MODE SELECTOR
C                  = 0 - NORMAL OPERATION
C                  = 1 - MONITOR PROGRESS OF DOPPLER BROADENING OF DATA.
C                        EACH TIME A PAGE OF DATA POINTS IS WRITTEN TO
C                        THE SCRATCH FILE PRINT OUT THE TOTAL NUMBER OF
C                        POINTS ON SCRATCH AND THE LOWER AND UPPER
C                        ENERGY LIMITS OF THE PAGE (THIS OPTION MAY BE
C                        USED IN ORDER TO MONITOR THE EXECUTION SPEED
C                        OF LONG RUNNING JOBS).
C           23-33  KELVIN TEMPERATURE
C           34-44  MINIMUM CROSS SECTION OF INTEREST
C                  (DEFAULT VALUE = 1.0E-10 BARNS).
C           45-55  NEGATIVE CROSS SECTION TREATMENT
C                  = 0 - O.K.
C                  = 1 - SET = 0
C           56-66  UNRESOLVED RESONANCE REGION TREATMENT
C                  = 0 - COPY (NO BROADENING)
C                  = 1 - IGNORE (BROADEN)
C        2   1-72  ENDF/B INPUT DATA FILENAME
C                  (STANDARD OPTION = ENDFB.IN)
C        3   1-72  ENDF/B OUTPUT DATA FILENAME
C                  (STANDARD OPTION = ENDFB.OUT)
C      4-N   1-11  LOWER MAT OR ZA LIMIT
C           12-22  UPPER MAT OR ZA LIMIT
C                  UP TO 100 MAT OR ZA RANGES MAY BE SPECIFIED, ONE
C                  RANGE PER CARD. THE LIST OF RANGES IS TERMINATED BY
C                  A BLANK CARD. IF THE UPPER LIMIT IS LESS THAN THE
C                  LOWER LIMIT THE UPPER LIMIT WILL BE SET EQUAL TO THE
C                  LOWER LIMIT. IF THE FIRST REQUEST CARD IS BLANK IT
C                  WILL TERMINATE THE LIST OF REQUESTS AND CAUSE ALL
C                  DATA TO BE RETRIEVED (SEE EXAMPLE INPUT).
C      VARY  1-11  ENERGY FOR ERROR LAW
C           12-22  ERROR FOR ERROR LAW
C                  THE ACCEPTABLE LINEARIZING ERROR CAN BE GIVEN AS AN
C                  ENERGY DEPENDENT FUNCTION SPECIFIED BY UP TO 20
C                  (ENERGY,ERROR) PAIRS AND LINEAR INTERPOLATION
C                  TABULATE POINTS. ENERGIES MUST BE IN ASCENDING ORDER.
C                  THE ERROR LAW IS TERMINATED BY A BLANK CARD. IF THE
C                  FIRST ERROR LAW CARD IS BLANK IT WILL TERMINATE THE
C                  ERROR LAW AND THE ERROR WILL BE TREATED AS ENERGY
C                  INDEPENDENT, EQUAL TO ZERO, WHICH INDICATES THAT THE
C                  BROADENED DATA SHOULD NOT BE THINNED.
C
C     EXAMPLE INPUT NO. 1
C     -------------------
C     BROADEN ALL URANIUM ISOTOPES AND THORIUM-232 TO 300 KELVIN. FROM
C     0 TO 100 EV THIN OUTPUT DATA TO 0.1 PER-CENT ACCURACY. FROM 100 EV
C     TO 1 KEV VARY THE ERROR BETWEEN 0.1 AND 1 PER-CENT. ABOVE 1 KEV
C     USE 1 PER-CENT ACCURACY.
C
C     EXPLICITLY SPECIFY THE STANDARD FILENAMES.
C
C     THE FOLLOWING 11 CARDS ARE REQUIRED
C
C          1          0 3.00000+ 2
C ENDFB.IN
C ENDFB.OUT
C      92000      92999
C      90232               (UPPER LIMIT WILL AUTOMATICALLY BE DEFINED)
C                          (BLANK CARD INDICATES END OF REQUEST LIST)
C 0.00000+ 0 1.00000-03
C 1.00000+ 2 1.00000-03
C 1.00000+ 3 1.00000-02
C 1.00000+ 9 1.00000-02
C                          (BLANK CARD INDICATES END OF ERROR LAW)
C
C     EXAMPLE INPUT NO. 2
C     -------------------
C     BROADEN ALL DATA TO 300 KELVIN AND DO NOT THIN THE BROADEN DATA.
C     ALL OF THE STANDARD OPTION MAY BE INVOKED MERELY BY SPECIFYING
C     THE KELVIN TEMPERATURE ON THE FIRST CARD. ALL OTHER FIELDS MAY
C     BE LEFT BLANK.
C
C     LEAVE THE DEFINITION OF THE FILENAMES BLANK - THE PROGRAM WILL
C     THEN USE STANDARD FILENAMES.
C
C     THE FOLLOWING 5 CARDS ARE REQUIRED
C
C                       3.00000+ 2
C                       (USE STANDARD FILENAME = ENDFB.IN)
C                       (USE STANDARD FILENAME = ENDFB.OUT)
C                       (RETRIEVE ALL DATA, TERMINATE REQUEST LIST)
C                       (0.0 ALLOWABLE ERROR, TERMINATE ERROR LAW)
C
C     EXAMPLE INPUT NO. 3
C     -------------------
C     THE SAME AS ABOVE, ONLY DEFINE THE MINIMUM CROSS SECTION OF
C     INTEREST TO BE 1.0E-30 BARNS (INSTEAD OF THE DEFAULT VALUE OF
C     1.0E-10).
C
C     READ ENDF/B DATA FROM \ENDFB6\RECENT\ZA092238 AND WRITE ENDF/B
C     DATA TO \ENDFB\SIGMA1\ZA092238
C
C     THE FOLLOWING 5 CARDS ARE REQUIRED
C
C                       3.00000+ 2 1.00000-30
C \ENDFB6\RECENT\ZA092238
C \ENDFB6\SIGMA1\ZA092238
C                       (RETRIEVE ALL DATA, TERMINATE REQUEST LIST)
C                       (0.0 ALLOWABLE ERROR, TERMINATE ERROR LAW)
C
C=======================================================================
      INCLUDE 'implicit.h'
C-----08/08/2012 DEFINE CODE NAME
      CHARACTER*8 CODENAME
      COMMON/NAMECODE/CODENAME
      INTEGER*4 OTAPE,OUTP,COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,
     1 HOT3M1,UNRES1,UNRES2,UREVAL,UREACT
      CHARACTER*1 FIELD6
      CHARACTER*4 CARD
CAK
      CHARACTER*72 infile,outfile
CAK
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
CAK   COMMON/IWATCH/MONITR,MAKEPLUS,MYUNRES
      COMMON/IWATCHs/MONITR,MAKEPLUS,MYUNRES
      COMMON/COPI/MFIELD(3)
      COMMON/COPC/CARD(17)
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HOTS/ALPHA,HOTSY1,HOTSY2,TEMPK,TEMPEF,N2TAPI,N2TAPO
      COMMON/MATTOT/MATIN,MATOUT
      COMMON/CONTAC/OVPI,OVPI2,ATOP
      COMMON/EXTEND/DTMAX,MESS
      COMMON/THRESH/ETHRES,EMIN
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
CAK   COMMON/PAGER/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/PAGERs/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/RESOLV/EULOW,EUHIGH,UREVAL,UREACT,UNRES1,UNRES2,IL,IH
CAK   COMMON/MAXIE/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/MAXIEs/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/PARTIN/ATWTP,INPART
      COMMON/OKERRT/ERXC3T,KERR3T,MAXERT,ENER3T(21),ER3T(21)
CAK   COMMON/FIELDC/FIELD6(11,6)
      COMMON/FIELDCs/FIELD6(11,6)
C-----08/20/2013 - added high energy cutoff
      COMMON/HIGHESTE/EHIGHEST,EMILLION
      INCLUDE 'sigma1.h'
C-----------------------------------------------------------------------
C
C     INITIALIZE
C
C-----------------------------------------------------------------------
C-----08/08/2012 DEFINE CODE NAME
      CODENAME = 'SIGMA1  '
C-----DEFINE PI TO DOUBLE PRECISION ACCURACY.
      PI = DACOS(-1.0D+00)
C-----08/20/2013 - added high energy cutoff
      EHIGHEST = 10.1D+6 ! this is eV (10.1 MeV, just above 10 MeV).
      EMILLION =  1.0D+6 ! this is multiple of AE/KT
C-----INITIALIZE TIMER
CAK   CALL TIMER
      CALL TIMERpre
C-----INITIALIZE ERFC TABLES
      CALL TABERFC
C-----------------------------------------------------------------------
C
C     DEFINE ALL I/O UNITS.
C
C-----------------------------------------------------------------------
CAK   CALL FILEIO
      CALL FILEIOs
C-----INITIALIZE LAST MAT READ.
      MATNOW=0
C-----DEFINE WHERE TO START HIGH ENERGY APPROXIMATE TREATMENT (NOTE,
C-----HOTSY2 CORRESPONDS TO SQRT(A*E/KT). THEREFORE SETTING HOTSY2 EQUAL
C-----TO 1,000 CORRESPONDS TO STARTING THE HIGH ENERGY APPROXIMATION AT
C-----A*E/KT = 1,000,000, E.G., FOR U-238 DOPPLER BROADENED TO ROOM
C-----TEMPERATURE, A=238, KT=0.0253, A/KT IS APPROXIMATELY 10,000 AND
C-----THE APPROXIMATION STARTS AT ABOUT 100 EV...SAVING A GREAT DEAL OF
C-----TIME WITH ESSENTIALLY NO LOSS OF ACCURACY IN THE RESOLVED
C-----RESONANCE REGION ABOVE 100 EV).
      HOTSY1=8.0
      HOTSY2=1000.0
C-----DEFINE STARTING LOCATION FOR COLD DATA POINTS (COLD DATA WILL
C-----ALWAYS START AT FIRST LOCATION IN TABLES).
      COLD1=1
C-----DEFINE PAGE SIZE (BY CHANGING DIMENSION, COMMON STATEMENTS
C-----AND THIS NUMBER THE SIZE OF THE CORE TABLES MAY BE CHANGED).
      NPAGE=MAXPAGE
      NPT2=2*NPAGE
      NPT3=3*NPAGE
      NP1P1=NPAGE+1
      NP2P1=NPT2+1
C-----DEFINE MINIMUM ENERGY OF INTEREST (IN EV).
      EMIN=1.0E-05
C-----DEFINE MINIMUM ENERGY OF THRESHOLD (IN EV).
      ETHRES=1.0
C-----DEFINE REQUIRED CONSTANTS.
      ATOP=5.0
      OVPI=1.0/DSQRT(PI)
      OVPI2=2.0*OVPI
C-----INITIALIZE MAXIMUM ALLOWABLE TEMPERATURE STEP (TO AVOID
C-----CROSS SECTION EXTENSION).
      DTMAX=0.0
C-----INITIALIZE TOTAL NUMBER OF FILE3 POINTS READ AND WRITTEN.
      N2TAPI=0
      N2TAPO=0
C-----------------------------------------------------------------------
C
C     READ ALL INPUT AND CREATE OUTPUT REPORT.
C
C-----------------------------------------------------------------------
C-----LIST TITLE FOR OUTPUT.
      WRITE(OUTP,180)
      WRITE(*   ,180)
C-----READ AND CHECK ALL INPUT PARAMETERS.
CAK   CALL READIN
      CALL READINs(infile,outfile,Tres)
C-----DEFINE MAXIMUM ALLOWABLE ENERGY SPACING TO APPROXIMATE 1/V.
C-----AT LOW AND HIGH ENERGY.
      DEMAXLOW=ESPACE(ER3T(1))
      DEMAXHI =ESPACE(ER3T(MAXERT))
C-----------------------------------------------------------------------
C
C     POSITION TO BEGINNING OF NEXT REQUESTED SECTION. IF BEGINNING OF
C     EVALUATION (NEW MAT) PRINT SUMMARY OF LAST EVALUATION (IF ANY)
C     AND INITIALIZE COUNTS AND FLAGS FOR NEW EVALUATION.
C
C-----------------------------------------------------------------------
C-----COPY TAPE LABEL TO FROM INPUT TO OUTPUT FILE.
      CALL COPYL
C-----LIST TAPE LABEL.
      WRITE(OUTP,190) CARD,MFIELD(1)
      WRITE(*   ,190) CARD,MFIELD(1)
C-----FIND NEXT REQUESTED MATERIAL.
CAK10 CALL NXTMAT
   10 CALL NXTMATs
      CALL CONTO
C-----CHECK FOR BEGINNING OF NEW MATERIAL.
      IF(MATH.eq.MATNOW) go to 60
      IF(MATNOW.LE.0) GO TO 40
C-----WRITE TOTALS FOR LAST MATERIAL (EITHER WITH OR WITHOUT UNRESOLVED
C-----REGION ENERGY LIMITS).
      IF(UREVAL.LE.0) GO TO 20
      CALL OUT9(EULOW ,FIELD6(1,1))
      CALL OUT9(EUHIGH,FIELD6(1,2))
      WRITE(OUTP,150) ((FIELD6(M,I),M=1,11),I=1,2),MATIN,MATOUT
      WRITE(*   ,150) ((FIELD6(M,I),M=1,11),I=1,2),MATIN,MATOUT
      GO TO 30
   20 WRITE(OUTP,160) MATIN,MATOUT
      WRITE(*   ,160) MATIN,MATOUT
C-----PRINT RUNNING TIME
   30 CALL TIMEMAT
C-----INCREMENT TOTAL POINT COUNTS FOR TAPE.
      N2TAPI=N2TAPI+MATIN
      N2TAPO=N2TAPO+MATOUT
C-----CHECK FOR END OF RUN.
   40 IF(MATH.le.0) go to 120
C-----SAVE CURRENT MAT NUMBER AND INITIALIZE COUNT OF THE NUMBER OF
C-----FILE 3 POINTS READ AND WRITTEN.
      MATNOW=MATH
      MATIN=0
      MATOUT=0
C-----INITIALIZE TO NEUTRON AS INCIDENT PARTICLE.
c-----3/22/2012 - INPART = 1000*Z + A = 1 for neutron, 0 for photon
      INPART=1
      CALL WHATIN
C-----INITIALIZE TO NO UNRESOLVED RESONANCE REGION AT THE BEGINNING
C-----OF EACH NEW EVALUATION. SECTION HEAD CARD ALREADY READ AND WRITTEN
      UREVAL=0
      EULOW=-1000.0
      EUHIGH=-1000.0
      GO TO 60
C-----READ HEAD CARD OF NEXT ENDF/B FILE AND CHECK FOR END OF
C-----MATERIAL.
   50 CALL CONTI
      CALL CONTO
      IF(MATH.le.0) go to 10
C-----------------------------------------------------------------------
C
C     PROCESS FILE 1 (COMMENTS), FILE 2 (RESONANCE PARAMETERS) AND
C     FILE 3 (CROSS SECTIONS).
C
C-----------------------------------------------------------------------
C-----FIND FILE 1, SECTION 451 AND ADD COMMENTS TO INDICATE THAT
C-----THIS MATERIAL HAS BEEN PROCESSED.
   60 IF(MFH.lt.1) go to 50
      IF(MFH.gt.1) go to 70
      IF(MTH.ne.451) go to 90
C-----ADD COMMENTS.
CAK   CALL FILE1
      CALL FILE1s
      GO TO 100
C-----FIND FILE 2 AND DEFINE THE ENERGY RANGE OF THE UNRESOLVED
C-----RESONANCE REGION.
   70 IF(MFH.lt.2) go to 100
      IF(MFH.gt.2) go to 80
C-----2016/11/25 - Added to only allow MT=151, resonance parameters
      IF(MTH.ne.151) go to 90
C-----IF REQUESTED, IGNORE UNRESOLVED RESONANCE REGION.
      IF(MYUNRES.NE.0) GO TO 100
CAK   CALL FILE2
      CALL FILE2s
      GO TO 100
C-----COPY UP TO FILE 3.
   80 IF(MFH.lt.3) go to 100
      IF(MFH.gt.3) go to 110
C-----A SECTION OF FILE 3 SECTION FOUND. DOPPLER BROADEN ONE SECTION.
CAK   CALL FILE3
      CALL FILE3s
C-----------------------------------------------------------------------
C
C     COPY TO END OF SECTION, FILE OR MATERIAL.
C
C-----------------------------------------------------------------------
C-----COPY TO SECTION END (SEND) CARD.
   90 CALL COPYS
      GO TO 50
C-----COPY TO NEXT FILE.
  100 CALL COPYF
      GO TO 50
C-----COPY REMAINDER OF MATERIAL.
  110 CALL COPYM
      GO TO 10
C-----------------------------------------------------------------------
C
C     END OF RUN
C
C-----------------------------------------------------------------------
C-----PRINT WARNING MESSAGE IF NO DATA WAS FOUND THAT SATISIFED
C-----REQUESTS.
  120 IF(MATNOW.GT.0) GO TO 130
      WRITE(OUTP,210)
      WRITE(*   ,210)
      CALL ENDERROR
C-----LIST TOTAL NUMBER OF FILE 3 POINTS READ AND WRITTEN.
  130 WRITE(OUTP,170) N2TAPI,N2TAPO
      WRITE(*   ,170) N2TAPI,N2TAPO
C-----WRITE WARNING MESSAGE IF CROSS SECTION EXTENSION USED.
      IF(DTMAX.LE.0.0) GO TO 140
      CALL OUT9(DTMAX,FIELD6(1,1))
      WRITE(OUTP,200) (FIELD6(M,1),M=1,11)
      WRITE(*   ,200) (FIELD6(M,1),M=1,11)
C-----END ENDF/B FORMAT OUTPUT FILE.
CAK
  140 CLOSE (OTAPE)
      RETURN
CAK
CAK  140 CALL ENDIT
      GO TO 140  ! CANNOT GET TO HERE.
  150 FORMAT(1X,79('-')/' Unresolved ',11A1,' to',11A1,
     1 ' eV',4X,'MAT Totals',2I9/1X,79('-'))
  160 FORMAT(1X,79('-')/' No Unresolved Region',23X,'MAT Totals',2I9/
     1 1X,79('-'))
  170 FORMAT(1X,79('-')/43X,'Tape Totals',2I9/1X,79('-'))
  180 FORMAT(' Doppler Broaden ENDF/B Cross Sections (SIGMA1 2017-1)'/
     1 1X,79('-'))
  190 FORMAT(1X,79('-')/' ENDF/B Tape Label'/1X,79('-')/1X,16A4,A2,I4/
     1 1X,79('-')/
     2 ' Projectile',4X,'  MAT  MT ENDF/B',5X,'Kelvin',5X,
     3 'Q-Value   Points   Points'/
     3 '       Material',9X,' Format',9X,'In',10X,'eV','       In',
     3 '      Out'/1X,79('-'))
  200 FORMAT(///' Extension'/1X,9('-')/
     1 ' Cross Section Extension Can be Avoided by'/
     2 ' Thinning Data or Doppler Broadening in Steps'/
     3 ' of Less than',11A1,' Kelvin')
  210 FORMAT(' WARNING - No Data Found that Satisfied Retrieval',
     1 ' Criteria.'/12X,
     2 ' Therefore No Data was Broadened or Written to Output File.'/
     3 1X,79('-'))
      END
      SUBROUTINE WHATIN
C=======================================================================
C
C     DEFINE MASS OF INCIDENT PROJECTILE RELATIVE TO THE MASS OF THE
C     NEUTRON.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*4 FMTHOL,PROHOL,PROTAB
      COMMON/PARTIN/ATWTP,INPART
C     COMMON/HOLFMT/FMTHOL,PROHOLK
      COMMON/HOLFMTs/FMTHOL,PROHOL
      DIMENSION NPART(6),ATWTN(6),PROTAB(6)
C-----DEFINE ZA AND MASS OF ALLOWABLE PARTICLES.
      DATA NPART/
     1     1, 1001, 1002, 1003, 2003, 2004/
C-----NEUTRON MASS UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITTEE
C-----RECOMMENDATION (NANCY LARSON)
c-----04/09/20 - Updated according to ENDF-102, Appendix H.
      DATA ATWTN/
     1 1.00866491578D+00,     ! Neutron
     2 1.00727646688D+00,     ! Proton
     3 2.01355321271D+00,     ! Deuteron
     4 3.016049268D+00,       ! Triton
     5 3.01493223469D+00,     ! He3
     6 4.0015061747D+00/      ! Alpha
c----Old values
c    1 1.008664904D+00,
c    2 1.007825D+00,
c    3 2.014102D+00,
c    4 3.016050D+00,
c    5 3.016030D+00,
c    6 4.002603D+00/
      DATA PROTAB/
     1 'n   ','p   ','d   ','t   ','He3 ','He4 '/
C-----LOOK UP ZA OF PROJECTILE.
      DO 10 I=1,6
      IF(INPART.EQ.NPART(I)) GO TO 30
   10 CONTINUE
C-----3/22/2012 - NOT FOUND = ERROR STOP
      WRITE(*,20) INPART
      WRITE(3,20) INPART
   20 FORMAT(///' ERROR...Incident Particle ZA=',I6,
     1 ' CANNOT be Doppler Broadened.'/
     2          '         Correct DATA and RE-TRY.'///)
      CALL ENDERROR
C-----FOUND. DEFINE MASS RELATIVE TO NEUTRON.
   30 ATWTP=ATWTN(I)/ATWTN(1)
      PROHOL=PROTAB(I)
      RETURN
      END
CAK   SUBROUTINE FILE1
      SUBROUTINE FILE1s
C=======================================================================
C
C     ADD COMMENTS AT THE END OF FILE1, SECTION 451 TO INDICATE
C     THAT THIS MATERIAL HAS BEEN PROCESSED BY PROGRAM SIGMA1 AND
C     TO SPECIFY THE TEMPERATURE AND MAXIMUM ALLOWABLE ERROR.
C
C     DEFINE FORMAT TO BE ENDF/B-IV, V OR VI.
C
C     THE ENDF/B FORMAT CAN BE DETERMINED FROM THE SECOND CARD.
C     ENDF/B-IV = N1 > 0, N2 = 0, CARD COUNT (POSITIVE)
C     ENDF/B-V  = N1 = N2 = 0
C     ENDF/B-IV =      N2 = NUMBER NUMBER (6 OR MORE)
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 PROGDOC1
      CHARACTER*4 FMTTAB,FMTHOL,PROHOL
      CHARACTER*66 PROGDOC
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/OKERRT/ERXC3T,KERR3T,MAXERT,ENER3T(21),ER3T(21)
      COMMON/HOTS/ALPHA,HOTSY1,HOTSY2,TEMPK,TEMPEF,N2TAPI,N2TAPO
C     COMMON/HOLFMT/FMTHOL,PROHOLK
      COMMON/HOLFMTs/FMTHOL,PROHOL
      COMMON/TEMPO/TEMP3,IVERSE
      COMMON/PARTIN/ATWTP,INPART
CAK   COMMON/PARAMS/XCMIN
      COMMON/PARAMSs/XCMIN
      DIMENSION FMTTAB(3),PROGDOC(9),PROGDOC1(66,9)
      EQUIVALENCE (PROGDOC(1),PROGDOC1(1,1))
      DATA FMTTAB/'IV  ','V   ','VI  '/
C-----DOCUMENTATION TO ADD TO ENDF/B OUTPUT - EACH LINE IS 66
C-----CHARACTERS LONG - FIELDS 12345678901 ARE FILLED IN WITH
C-----11 CHARACTERS DURING EXECUTION.
C               1         2         3         4         5         6
C       12345678901234567890123456789012345678901234567890123456789012
C       3456
      DATA PROGDOC/
     1 ' ***************** Program SIGMA1 (VERSION 2017-1) ***********',
     2 ' Data Doppler Broadened to12345678901 Kelvin                  ',
     3 ' for All Data Greater than12345678901 barns in Absolute Value ',
     4 ' Data Linearized to Within an Accuracy of12345678901 per-cent ',
     5 ' Data Linearized Using Energy Dependent Uncertainty           ',
     6 '      Energy    Accuracy                                      ',
     7 '        (eV)  (per-cent)                                      ',
     8 ' ----------- -----------                                      ',
     9 ' 12345678901 12345678901                                      '/
C-----FILL IN REMAINDER OF FIRST LINE.
      PROGDOC1(63,1) = '*'
      PROGDOC1(64,1) = '*'
      PROGDOC1(65,1) = '*'
      PROGDOC1(66,1) = '*'
C-----HEAD CARD OF SECTION HAS BEEN READ AND WRITTEN. READ NEXT CARD
C-----AND DETERMINE IF THIS IS THE ENDF/B-IV, V OR VI FORMAT.
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      IVERSE=4
C-----CHECK FOR ENDF/B-IV.
C-----IV N1 > 0, N2 = 0
      IF(N1.GT.0.AND.N2.EQ.0) GO TO 10
C-----NOT ENDF/B-IV. READ THIRD CARD.
      N2X=N2
      CALL CARDO(C1,C2,L1,L2,N1,N2)
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      IVERSE=5
C-----CHECK FOR ENDF/B-V FORMAT.
      IF(N2X.LE.0) GO TO 10
      N1X=N1
C-----ENDF/B-VI FORMAT. READ FOURTH CARD.
      CALL CARDO(C1,C2,L1,L2,N1,N2)
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      IVERSE=6
C-----DEFINE INCIDENT PARTICLE AND MASS RELATIVE TO NEUTRON.
      INPART=N1X/10
c---- 3/22/2012 = 0 = WRONG - that mean photon - use = 1 = neutron
      IF(INPART.LE.0) INPART=0  ! Default to Photon
      CALL WHATIN
C-----DEFINE INITIAL TEMPERATURE FROM C1 FIELD.
      TEMP3=C1
C-----IF FINAL TEMPERATURE IS HIGHER THAN INITIAL TEMPERATURE DEFINE
C-----FINAL TEMPERATURE.
      IF(TEMPK.GT.C1) C1=TEMPK
C-----SET DERIVED MATERIAL FLAG IF OUTPUT IS FOR POSITIVE TEMPERATURE.
      IF(C1.GT.0.0) L1=1
C-----DEFINE ENDF/B FORMAT NUMBER.
   10 FMTHOL=FMTTAB(IVERSE-3)
C-----INCREASE COMMENT CARD COUNT AND COPY TO END OF HOLLERITH.
      IF(MAXERT.LE.1) N1OUT=N1+4
      IF(MAXERT.GT.1) N1OUT=N1+7+MAXERT
      CALL CARDO(C1,C2,L1,L2,N1OUT,N2)
      DO 20 N=1,N1
   20 CALL COPY1
C-----------------------------------------------------------------------
C
C     ADD COMMENTS TO DOCUMENT WHAT WAS DONE TO DATA
C
C-----------------------------------------------------------------------
C-----OUTPUT PROGRAM NAME AND VERSION I.D.
      CALL HOLLYO(PROGDOC1(1,1))
C-----OUTPUT FINAL TEMPERATURE IN KELVIN
      CALL OUT9(TEMPK,PROGDOC1(27,2))
      CALL HOLLYO(PROGDOC1(1,2))
C-----OUTPUT MINIMUM CROSS SECTION
      CALL OUT9(XCMIN,PROGDOC1(27,3))
      CALL HOLLYO(PROGDOC1(1,3))
      IF(MAXERT.GT.1) GO TO 30
C-----ENERGY INDEPENDENT UNCERTAINTY
      PERCNT=100.0*ER3T(1)
      CALL OUT9(PERCNT,PROGDOC1(42,4))
      CALL HOLLYO(PROGDOC1(1,4))
      RETURN
C-----ADD FOUR COMMENT CARDS PLUS ENERGY DEPENDENT UNCERTAINTY
   30 CALL HOLLYO(PROGDOC1(1,5))
      CALL HOLLYO(PROGDOC1(1,6))
      CALL HOLLYO(PROGDOC1(1,7))
      CALL HOLLYO(PROGDOC1(1,8))
      DO 40 I=1,MAXERT
      PERCNT=100.0*ER3T(I)
      CALL OUT9(ENER3T(I),PROGDOC1( 2,9))
      CALL OUT9(PERCNT   ,PROGDOC1(14,9))
   40 CALL HOLLYO(PROGDOC1(1,9))
      RETURN
      END
CAK   SUBROUTINE FILE2
      SUBROUTINE FILE2s
C=======================================================================
C
C     READ RESONANCE PARAMETERS IN ORDER TO DEFINE THE ENERGY RANGE
C     OF THE UNRESOLVED RESONANCE REGION.
C
C     INSURE THAT UNRESOLVED REGION IS THE SAME FOR ALL ISOTOPES AND
C     SPIN STATES. IF NOT, PRINT WARNING AND USE MINIMUM UNRESOLVED
C     REGION FOR ALL ISOTOPES AND SPIN STATES.
C
C     UNRESOLVED ARGUMENTS
C     --------------------
C     UREVAL = EVALUATION UNRESOLVED RESONANCE REGION INDICATOR.
C            = 0 - EVALUATION DOES NOT HAVE UNRESOLVED REGION.
C            = 1 - EVALUATION DOES HAVE UNRESOLVED REGION.
C     UREACT = REACTION UNRESOLVED RESONANCE REGION INDICATOR.
C            = UREVAL - IF REACTION (MT) IS TOTAL, ELASTIC, CAPTURE
C              OR FISSION.
C            = 0 - OTHERWISE.
C     UNRES1 = INDEX TO LAST DATA POINT BELOW LOWER ENERGY LIMIT OF
C              UNRESOLVED RESONANCE REGION.
C     UNRES2 = INDEX TO FIRST DATA POINT ABOVE UPPER ENERGY LIMIT OF
C              UNRESOLVED RESONANCE REGION.
C     EULOW  = LOWER ENERGY LIMIT OF UNRESOLVED RESONANCE REGION.
C     EUHIGH = UPPER ENERGY LIMIT OF UNRESOLVED RESONANCE REGION.
C     IL     = COUNT OF THE NUMBER OF DATA POINTS AT THE LOWER ENERGY
C              LIMIT OF THE UNRESOLVED RESONANCE REGION (IF NOT AT LEAST
C              TWO PROGRAM WILL INSERT EXTRA POINTS).
C     IH     = COUNT OF THE NUMBER OF DATA POINTS AT THE UPPER ENERGY
C              LIMIT OF THE UNRESOLVED RESONANCE REGION (IF NOT AT LEAST
C              TWO PROGRAM WILL INSERT EXTRA POINTS).
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*4 FMTHOL,PROHOL
      INTEGER*4 OUTP,OTAPE,UNRES1,UNRES2,UREVAL,UREACT
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RESOLV/EULOW,EUHIGH,UREVAL,UREACT,UNRES1,UNRES2,IL,IH
CAK   COMMON/WHATZA/IZA
      COMMON/WHATZAs/IZA
CAK   COMMON/HOLFMT/FMTHOL,PROHOL
      COMMON/HOLFMTs/FMTHOL,PROHOL
C-----INITIALIZE ERROR FLAG OFF.
      IBAD=0
C-----HEAD RECORD ALREADY READ. DEFINE NUMBER OF ISOTOPES.
      NIS=N1H
C-----DO LOOP OVER ALL ISOTOPES
      DO 270 IS=1,NIS
      CALL CARDIO(C1H,C2H,L1H,LFW,NER,N2H)
C-----DO LOOP OVER ALL ENERGY RANGES
      DO 260 JER=1,NER
      CALL CARDIO(EL,EH,LRU,LRF,N1H,N2H)
C
C     2017/5/16 - Added energy dependent scattering radius
C
      NRO = N1H
      if(NRO.eq.1) then
c-----Energy dependent scattering radius = copy TAB1 record
      CALL CARDIO(C1,C2,L1,L2,N1,N2)
      do i = 1,N1,3     ! Interpolation law (NBT,INT) Pairs
      CALL COPY1
      enddo
      do i = 1,N2,3     ! Scattering radius (X,Y) Pairs
      CALL COPY1
      enddo
      endif
C
C-----DEFINE LRU FOR INTERNAL USE AS ORIGINAL LRU (BEFORE RECENT).
      IF(LRU.GT.3) LRU=LRU-3
C-----Select resolved or unresolved.
      IF(LRU.eq.1) go to 60    ! Resolved
      IF(LRU.gt.1) go to 10    ! unresolved
C
C     NO RESONANCE PARAMETERS PRESENT
C
C-----COPY SECTION WITH NO RESONANCE PARAMETERS.
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,N2H)
      GO TO 260
C-----------------------------------------------------------------------
C
C     UNRESONANCE
C
C-----------------------------------------------------------------------
C
C     ALLOW FOR MULTIPLE, ADJACENT UNRESOLVED REGIONS.
C
C-----CHECK FOR ONLY ONE UNRESOLVED RESONANCE REGION.
   10 IF(UREVAL.LE.0) GO TO 50
C-----MULTIPLE UNRESOLVED - FIRST CHECK FOR SAME ENERGY RANGE.
      IF(DABS(EULOW-EL) .LE.0.0001*DABS(EULOW).AND.
     1   DABS(EUHIGH-EH).LE.0.0001*DABS(EUHIGH)) GO TO 60
C-----MULTIPLE UNRESOLVED - NEXT CHECK FOR ADJACENT ENERGY RANGE.
      IF(DABS(EUHIGH-EL).GT.0.0001*DABS(EUHIGH)) GO TO 30
C-----EUHIGH OF LAST = EL OF NEXT = ADJACENT RANGES - EXTEND RANGE UP
      WRITE(OUTP,20) EULOW,EUHIGH,EL,EH
      WRITE(*   ,20) EULOW,EUHIGH,EL,EH
   20 FORMAT(1X,79('-')/
     1       ' WARNING - Combining Adjacent Unresolved Ranges'/
     1       '          ',1pd11.4,' eV to ',1pd11.4,' eV and'/
     1       '          ',1pd11.4,' eV to ',1pd11.4,' eV'/1X,79('-'))
      EUHIGH = EH
      GO TO 60
   30 IF(DABS(EULOW-EH).GT.0.0001*DABS(EULOW)) GO TO 40
C-----EULOW OF LAST = EH OF NEXT = ADJACENT RANGES - EXTEND RANGE DOWN
      WRITE(OUTP,20) EL,EH,EULOW,EUHIGH
      WRITE(*   ,20) EL,EH,EULOW,EUHIGH
      EULOW = EL
      GO TO 60
C-----MULTIPLE UNRESOLVED INCOMPATIBLE RANGES
   40 IBAD=1
      GO TO 60
C-----SET FLAG TO INDICATE THE PRESENCE OF AN UNRESOLVED RESONANCE
C-----REGION AND DEFINE ENERGY LIMITS OF THE UNRESOLVED REGION.
   50 UREVAL=1
      EULOW=EL
      EUHIGH=EH
C-----------------------------------------------------------------------
C
C     RESONANCE PARAMETERS PRESENT
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C     LRU= 1 - RESOLVED
C     LRF= 1 - SLBW, = 2 - MLBW, = 3 - REICH-MOORE, = 4 - ADLER-ADLER
C            - NEW REICH-MOORE = 7
C
C     LRU= 2 - UNRESOLVED
C     LRF= 1 - ENERGY INDEPENDENT WIDTHS (EXCEPT POSSIBLY FISSION)
C        = 2 - ENERGY   DEPENDENT WIDTHS
C
C-----------------------------------------------------------------------
   60 IF(LRU.NE.1) GO TO 70    ! Resolved?
      IF(LRF.EQ.1.OR.          ! Single Level Breit-Wigner
     1   LRF.EQ.2.OR.          ! Multi-Level  Breit-Wigner
     2   LRF.EQ.3) GO TO 80    ! Reich=Moore
      IF(LRF.EQ.4) GO TO 130   ! Adler-Adler
      IF(LRF.EQ.7) GO TO 160   ! New Reich-Moore
C-----ILLEGAL - IGNORE REMAINDER OF FILE 2
      GO TO 280
   70 IF(LRU.NE.2) GO TO 280   ! Unresolved?
      IF(LRF.EQ.1) GO TO 180   ! Energy Independent Widths
      IF(LRF.EQ.2) GO TO 240   ! Energy   Dependent Widths
C-----ILLEGAL - IGNORE REMAINDER OF FILE 2
      GO TO 280
C-----------------------------------------------------------------------
C
C     BREIT-WIGNER (SINGLE OR MULTI-LEVEL) OR REICH-MOORE FORMALISM
C
C-----------------------------------------------------------------------
C-----IF SCATTERING RADIUS IS ENERGY DEPENDENT COPY TABULATED VALUES.
   80 IF(N1H.EQ.0) GO TO 110
      IF(LRF.NE.1.AND.LRF.NE.2) GO TO 110
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,N2H)
C-----SKIP SCATTERING RADIUS INTERPOLATION LAW.
      DO 90 I=1,N1H,3
   90 CALL COPY1
C-----SKIP SCATTERING RADIUS TABULATED VALUES.
      DO 100 I=1,N2H,3
  100 CALL COPY1
C-----READ NEXT CARD.
  110 CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
C-----LOOP OVER ALL L STATES
      DO 120 ILS=1,NLS
C-----READ NEXT CARD.
      CALL CARDIO(C1H,C2H,L1H,L2H,NRS6,NRS)
C-----COPY RESONANCE PARAMETERS.
      DO 120 IRS=1,NRS
  120 CALL COPY1
      GO TO 260
C-----------------------------------------------------------------------
C
C     ADLER-ADLER FORMALISM
C
C-----------------------------------------------------------------------
C-----READ NEXT CARD.
  130 CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
C-----READ BACKGROUND CORRECTIONS.
      CALL CARDIO(C1H,C2H,L1H,L2H,NX6,N2H)
C-----COPY BACKGROUND CORRECTION CONSTANTS.
      DO 140 I=1,NX6,6
  140 CALL COPY1
C-----LOOP OVER L STATES
      DO 150 I=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,NJS,N2H)
C-----LOOP OVER J STATES
      DO 150 J=1,NJS
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,NLJ)
C-----COPY ALL RESONANCE DATA
      DO 150 K=1,NLJ
      CALL COPY1
  150 CALL COPY1
      GO TO 260
C-----------------------------------------------------------------------
C
C     NEW REICH-MOORE FORMALISM
C
C-----------------------------------------------------------------------
C-----DEFINE NUMBER OF J STATES
  160 CALL CARDIO(C1,C2,L1,L2,NJS,N2)
C-----DEFINE NUMBER OF PARTICLE-PAIRS
      CALL CARDIO(C1,C2,L1,L2,NPP12,N2)
C-----COPY PARTICLE-PAIR DATA
      DO N=1,NPP12,6
      CALL COPY1
      ENDDO
C-----LOOP OVER J STATES
      DO 170 IJ=1,NJS
C-----J, PARITY, AND NUMBER OF CHANNELS
      CALL CARDIO(C1,C2,L1,L2,NCH6,N2)
C-----COPY CHANNEL DATA
      DO N=1,NCH6,6
      CALL COPY1
      ENDDO
C-----DEFINE NUMBER OF RESONANCES
      CALL CARDIO(C1,C2,L1,L2,N1,NRS)
C-----COPY RESONANCE PARAMETERS
      DO N=1,NRS
      CALL COPY1
      ENDDO
  170 CONTINUE
      GO TO 280
C-----------------------------------------------------------------------
C
C     UNRESOLVED WITH ENERGY INDEPENDENT WIDTHS.
C
C-----------------------------------------------------------------------
C-----TEST IF FISSION WIDTHS GIVEN
  180 IF(LFW.gt.0) go to 200
C-----FISSION WIDTHS NOT GIVEN
      CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
C-----LOOP OVER ALL L-STATES
      DO 190 ILS=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,NJS)
      DO 190 N=1,NJS
  190 CALL COPY1
      GO TO 260
C-----FISSION WIDTHS GIVEN (LFW=1)
  200 CALL CARDIO(C1H,C2H,L1H,L2H,NE,NLS)
C-----COPY FISSION WIDTH ENERGY POINTS
      DO 210 I=1,NE,6
  210 CALL COPY1
C-----LOOP OVER L-STATES
      DO 230 I=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,NJS,N2H)
C-----LOOP OVER J STATES
      DO 230 J=1,NJS
      CALL CARDIO(C1H,C2H,L1H,L2H,NEP6,N2H)
      DO 220 K=1,NEP6,6
  220 CALL COPY1
  230 CONTINUE
      GO TO 260
C-----------------------------------------------------------------------
C
C     UNRESOLVED WITH ENERGY DEPENDENT WIDTHS.
C
C-----------------------------------------------------------------------
C-----READ NEXT CARD.
  240 CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
C-----DO LOOP OVER L-STATES
      DO 250 I=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,NJS,N2H)
      DO 250 J=1,NJS
      CALL CARDIO(C1H,C2H,L1H,L2H,NE6P6,N2H)
C-----COPY NUMBER OF DEGREES OF FREEDOM AND PARAMETERS.
      DO 250 K=1,NE6P6,6
  250 CALL COPY1
  260 CONTINUE
  270 CONTINUE
C-----------------------------------------------------------------------
C
C     END OF RESONANCE REGION (FILE 2) DATA.
C
C-----------------------------------------------------------------------
C-----PRINT WARNING IF UNRESOLVED RESONANCE REGION IS NOT UNIQUE.
  280 IF(UREVAL.LE.0) GO TO 300
      IF(EULOW.LT.EUHIGH) GO TO 290
      UREVAL=0
      WRITE(OUTP,320) IZA,MATH,MTH,FMTHOL
      WRITE(*   ,320) IZA,MATH,MTH,FMTHOL
      GO TO 300
  290 IF(IBAD.GT.0) WRITE(OUTP,310) IZA,MATH,MTH,FMTHOL
      IF(IBAD.GT.0) WRITE(*   ,310) IZA,MATH,MTH,FMTHOL
  300 RETURN
  310 FORMAT(1X,I6,2I5,2X,A2,
     1 ' WARNING - Unresolved Resonance Energy Range NOT the Same'/
     2 29X,' for ALL Isotopes and Spin States. This Violates ENDF/B'/
     3 29X,' Conventions. The Unresolved Resonance Energy Range'/
     4 29X,' will be Considered to Extend from the Maximum Lower'/
     5 29X,' Energy Limit up to the Mimimum Upper Energy Limit'/
     6 29X,' of ALL the Unresolved Ranges Defined for this MAT.')
  320 FORMAT(1X,I6,2I5,2X,A2,
     1 ' WARNING - Unresolved Resonance Energy Range NOT the Same'/
     2 29X,' for ALL Isotopes and Spin States. This Violates ENDF/B'/
     3 29X,' Conventions. Cannot Locate ANY Energy Range that'/
     4 29X,' ONLY Contains Unresolved Resonances. Will Ignor'/
     5 29X,' Unresolved Resonance Region.')
      END
CAK   SUBROUTINE FILE3
      SUBROUTINE FILE3s
C=======================================================================
C
C     THIS ROUTINE IS DESIGNED TO DOPPLER BROADEN ONE ENDF/B SECTION
C     OF DATA (I.E., ONE REACTION). DATA IS READ, DOPPLER BROADENED AND
C     OUTPUT IN THE ENDF/B FORMAT. IF THE SECTION CONTAINS 180000 FEWER
C     POINTS THE ENTIRE OPERATION IS PERFORMED IN CORE. IF THE SECTION
C     CONTAINS MORE 180000 POINTS THE DATA WILL BE BROADENED A PAGE
C     (1 PAGE = 60000 POINTS) AT A TIME, WRITTEN TO SCRATCH AND AFTER
C     THE ENTIRE SECTION HAS BEEN BROADENED IT WILL BE READ BACK FROM
C     SCRATCH AND OUTPUT IN THE ENDF/B FORMAT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,TOOHI,COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,
     1 HOT3,HOT3M1,UNRES1,UNRES2,UREVAL,UREACT
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
CAK   COMMON/PAGER/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/PAGERs/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/HOTS/ALPHA,HOTSY1,HOTSY2,TEMPK,TEMPEF,N2TAPI,N2TAPO
      COMMON/LOGLOGE/ELOGLOG
      COMMON/HEADER/ZA,AWR,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/C1,Q,L1,L2,N1,N2,MAT,MF,MT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/FILLER/N2LEFT,TOOHI,ITHRES,LOAD1,LOAD2
      COMMON/EXTEND/DTMAX,MESS
      COMMON/RESOLV/EULOW,EUHIGH,UREVAL,UREACT,UNRES1,UNRES2,IL,IH
      COMMON/LASTE/ELAST
      COMMON/TEMPO/TEMP3,IVERSE
      COMMON/SLIM/ISTART,NOTHIN,ITHIN1,ITHIN2,ITHIN3,MTEND
      COMMON/PARTIN/ATWTP,INPART
CAK   COMMON/MAXIE/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/MAXIEs/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
c----- 2012 - number of points read in any MT
      COMMON/READN2/N2READ
      INCLUDE 'sigma1.h'
C-----DEFINE BOLTZMANN CONSTANT IN EV/DEGREE KELVIN
C-----(ASSUMING ENERGY WILL BE EXPRESSED IN EV AND THE ATOMIC
C-----WEIGHT RATIO WILL BE EXPRESSED IN NEUTRON MASS UNITS,
C-----AS OPPOSED TO ATOMIC MASS UNITS (AMU)).
C-----UPDATED NOV. 12, 1998 AS PER CSEWG SUBCOMMITEE RECOMMENDATION
C***** DEVUG
      DATA BOLTZM/8.617385D-05/
C     DATA BOLTZM/8.6164D-05/
C***** DEVUG
C-----INITIALIZE TEMPERATURE READ FROM TOTAL CROSS SECTION (MT=1).
      DATA TEMP1/0.0D+00/
C-----INITIALIZE MAXIMUM ENERGY SPACING FOR LOW ENERGY START.
      DEMAXC=DEMAXLOW
      DEMAXH=DEMAXLOW
C-----------------------------------------------------------------------
C
C     READ TAB1 LEADER CARD AND INTERPOLATION LAW. CHECK INTERPOLATION
C     LAW AND IF DATA IS NOT LINEARLY INTERPOLABLE TERMINATE EXECUTION.
C
C-----------------------------------------------------------------------
C-----READ TAB1 LEAD CARD.
      CALL CARDI(C1,Q,L1,L2,N1,N2)
C-----READ INTERPOLATION LAW.
      CALL TERPI(NBTF,INTF,N1)
C-----INSURE DATA ARE LINEAR-LINEAR INTERPOLABLE.
      DO 10 K=1,N1
      IF(INTF(K).NE.2) GO TO 20
   10 CONTINUE
      GO TO 30
C-----DATA ARE NOT LINEAR-LINEAR. TERMINATE EXECUTION.
   20 WRITE(OUTP,130) MATH,MFH,MTH,(NBTF(I),INTF(I),I=1,N1)
      WRITE(OUTP,140)
      CALL ENDERROR
C-----------------------------------------------------------------------
C
C     INITIALIZE ALL COUNTS AND FLAGS FOR SECTION.
C
C-----------------------------------------------------------------------
C-----INITIALIZE FLAG NOT TO PRINT EXTENSION MESSAGE WITH SECTION.
   30 MESS=1
C-----INITIALIZE FLAG TO DOPPLER BROADEN SECTION.
      TOOHI=-1
C-----TURN OFF DOPPLER BROADENING FLAG IF SECTION IS MU-BAR, XI OR
C-----GAMMA.
c----- 2012 - number of points read in any MT
      N2READ = N2
      IF(MTH.GE.251.AND.MTH.LE.253) TOOHI=1
C-----------------------------------------------------------------------
C
C     DEFINE ORIGINAL TEMPERATURE. FOR ENDF/B-V AND EARLIER VERSIONS
C     C1 OF THE TAB1 LEAD CARD IS EITHER THE TEMPERATURE (L2=0) OR THE
C     Q-VALUE CORRESPONDING TO THE REACTION WITH MT=N2. IN THE LATTER
C     CASE ASSUME THE SAME TEMPERATURE AS FOR THE TOTAL. FOR ENDF/B-VI
C     THE INITIAL TEMPERATURE HAS ALREADY BEEN READ FROM FILE 1 (TEMP3).
C
C-----------------------------------------------------------------------
C-----CHECK FOR ENDF/B-VI FORMAT.
      IF(IVERSE.NE.6) GO TO 40
C-----ENDF/B-VI FORMAT. INITIAL TEMPERATURE HAS ALREADY BEEN DEFINED
C-----FROM FILE 1 (TEMP3).
      TEMPEF=TEMPK-TEMP3
      GO TO 60
C-----ENDF/B-V OR EARLIER. USE EITHER C1 FIELD OR TEMPERATURE FROM
C-----TOTAL.
   40 IF(L2.LE.0) GO TO 50
      TEMP3=TEMP1
      TEMPEF=TEMPK-TEMP3
      GO TO 60
   50 TEMP3=C1
      IF(MTH.EQ.1) TEMP1=TEMP3
      TEMPEF=TEMPK-TEMP3
C-----IN ENDF/B-V OR EARLIER DEFINE TEMPERATURE IN C1 FIELD.
      IF(TEMPEF.GT.0.0.AND.TOOHI.LT.0) C1=TEMPK
C-----SHOULD DATA BE DOPPLER BROADENED....
   60 IF(TEMPEF.gt.0.0D+0) go to 70
C-----NO.
      TOOHI=1
      GO TO 80
C-----YES.
   70 ALPHA=(AWR/ATWTP)/(BOLTZM*TEMPEF)
C-----DEFINE CUTOFF ENERGY FOR LOW ENERGY LOG-LOG ASSUMPTION.
      IF(MTH.LE.2) THEN
      ELOGLOG = 0.0D+0             ! total an elastic
      ELSE
c-----this is the same as older versions, but easier and faster
c-----E * a   < 40,   a = atwt/(kT) ~ 40*atwt (room kt~ 1/40 eV)
c-----E *atwt < 1
c-----E       < 1/atwt
      ELOGLOG = 1.0D+0*(ATWTP/AWR) ! capture, fission, etc.
      ENDIF
C-----INITIALIZE FLAG TO INDICATE BEGINNING OF SECTION.
   80 ISTART=1
C-----INITIALIZE TOTAL NUMBER OF POINTS IN SECTION AND NUMBER OF POINTS
C-----LEFT TO READ.
      N2IN=N2
      N2LEFT=N2
C-----INITIALIZE BROADENING INDICES TO FIRST TWO PAGES.
      HOT1=1
      HOT2=NPT2
C-----INITIALIZE END OF SECTION FLAG.
      MTEND=0
C-----INITIALIZE THINNING INDICES.
      ITHIN1=1
      ITHIN2=2
      ITHIN3=2
C-----INITIALIZE COUNT OF POINTS ON SCRATCH.
      N2SCR=0
C-----INITIALIZE LAST ENERGY READ FOR ASCENDING ENERGY TEST.
      ELAST=0.0D+00
C-----------------------------------------------------------------------
C
C     INITIALIZE UNRESOLVED RESONANCE REGION PARAMETERS FOR THIS SECTION
C
C-----------------------------------------------------------------------
C-----IF THERE IS AN UNRESOLVED RESONANCE REGION AND THIS SECTION IS
C-----TOTAL, ELASTIC, FISSION OR CAPTURE INITIAL INDICES.
      UREACT=UREVAL
      UNRES1=10000000
      UNRES2=10000000
      IL=0
      IH=0
      IF(UREACT.LE.0) GO TO 90
      IF(MTH.NE.1.AND.MTH.NE.2.AND.MTH.NE.18.AND.MTH.NE.19.AND.
     1 MTH.NE.102) UREACT=0
C-----------------------------------------------------------------------
C
C     LOAD DATA.
C
C-----------------------------------------------------------------------
C-----LOAD NEXT PAGE OF DATA AT ORIGINAL TEMPERATURE.
   90 CALL FILLUP
C-----SHOULD DATA BE DOPPLER BROADENED....
      IF(TOOHI.LE.0) GO TO 100
C-----NO. IF END OF SECTION BRANCH FOR FINAL THIN AND OUTPUT.
      IF(MTEND.ne.0) go to 120
C-----OTHERWISE THIN DATA AND LOAD NEXT PAGES.
      CALL THINIT
      HOT1=NP1P1
      GO TO 90
C-----------------------------------------------------------------------
C
C     BROADEN AND THIN DATA.
C
C-----------------------------------------------------------------------
C-----YES. INITIALIZE INDICES TO FIRST AND LAST POINTS TO USE IN DOPPLER
C-----BROADENING INTEGRAL TO POINT TO FIRST AND LAST POINTS THAT ARE
C-----LOADED IN CORE.
  100 COLD1P=COLD1
      COLD2P=COLD2
C-----IF UNRESOLVED RESONANCE REGION DATA IS IN CORE TRUNCATE RANGE
C-----OF INTEGRATION AT UNRESOLVED RESONANCE REGION ENERGY BOUNDARY.
C-----FOR ENERGIES BELOW UNRESOLVED REGION SET UPPER INTEGRATION LIMIT
C-----TO LOWER ENERGY LIMIT OF UNRESOLVED REGION.
      IF(HOT1.LE.UNRES1.AND.UNRES1.LE.COLD2) COLD2P=UNRES1
C-----FOR ENERGIES ABOVE UNRESOLVED REGION SET LOWER INTEGRATION LIMIT
C-----TO UPPER ENERGY LIMIT OF UNRESOLVED REGION.
      IF(HOT1.GE.UNRES2.AND.UNRES2.GT.0) COLD1P=UNRES2
C-----BROADEN DATA.
      CALL BROADN
C-----HAS THE LAST POINT BEEN BROADENED....
      IF(MTEND.ne.0) go to 120
C-----------------------------------------------------------------------
C
C     MORE POINTS TO BROADEN. SET UP DATA AND INDICES FOR NEXT PAGE
C     TO BE BROADENED.
C
C-----------------------------------------------------------------------
C-----SHIFT PAGES TWO AND THREE OF COLD DATA FORWARD ONE PAGE IN CORE.
      KK=0
      DO 110 K=NP1P1,COLD2
      KK=KK+1
      ECOLD(KK)=ECOLD(K)
      YCOLD(KK)=YCOLD(K)
      XCCOLD(KK)=XCCOLD(K)
  110 DCOLD(KK)=DCOLD(K)
C-----SET INDEX TO FIRST POINT TO DOPPLER BROADEN TO BEGINNING OF
C-----SECOND PAGE.
      HOT1=NP1P1
C-----IF ANY POINTS ARE WITHIN THE UNRESOLVED RESONANCE REGION SHIFT
C-----INDICES TO UNRESOLVED RESONANCE REGION FORWARD ONE PAGE.
      IF(UREACT.LE.0) GO TO 90
      UNRES1=UNRES1-NPAGE
      UNRES2=UNRES2-NPAGE
      GO TO 90
C-----------------------------------------------------------------------
C
C     END OF SECTION. THIN REMAINDER OF SECTION AND OUTPUT SECTION
C     (FROM CORE OR SCRATCH).
C
C-----------------------------------------------------------------------
  120 MTEND=1
      CALL THINIT
      CALL COPOUT
      RETURN
  130 FORMAT(///'  Interpolation Law is NOT Linear-Linear'/
     1 '  MAT/MF/MT=',I5,I3,I4/
     2 '  ---------------------------------------'/'      NBT  INT'/
     3 '  ---------------------------------------'/(I9,I5))
  140 FORMAT(//'  Execution Terminated'///)
      END
      SUBROUTINE FILLUP
C=======================================================================
C
C     LOAD NEXT PAGE OR PAGES OF CROSS SECTIONS AT THE ORIGINAL
C     TEMPERATURE INTO CORE. INSURE THAT THE MAXIMUM ENERGY
C     SPACING REQUIRED FOR LINEAR-LINEAR INTERPOLATABLE DOPPLER
C     BROADENED DATA IS NOT EXCEEDED.
C
C     IF THERE IS AN UNRESOLVED RESONANCE REGION INSURE THAT THERE IS
C     AT LEAST TWO POINTS AT THE LOWER AND UPPER ENERGY LIMITS OF THE
C     UNRESOLVED REGION AND DEFINE INDICES TO LOWER AND UPPER ENERGY
C     LIMITS OF THE UNRESOLVED REGION. DO NOT ADD ADDITIONAL ENERGY
C     POINTS WITHIN THE UNRESOLVED RESONANCE REGION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 TOOHI,COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,
     1 UNRES1,UNRES2,UREVAL,UREACT
      COMMON/HEADER/C1H,AWR,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/HOTS/ALPHA,HOTSY1,HOTSY2,TEMPK,TEMPEF,N2TAPI,N2TAPO
      COMMON/LOGLOGE/ELOGLOG
      COMMON/SLIM/ISTART,NOTHIN,ITHIN1,ITHIN2,ITHIN3,MTEND
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
CAK   COMMON/PAGER/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/PAGERs/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/FILLER/N2LEFT,TOOHI,ITHRES,LOAD1,LOAD2
      COMMON/CONTAC/OVPI,OVPI2,ATOP
      COMMON/THRESH/ETHRES,EMIN
      COMMON/FILLXY/EIN(1002),XCIN(1002)
      COMMON/RESOLV/EULOW,EUHIGH,UREVAL,UREACT,UNRES1,UNRES2,IL,IH
CAK   COMMON/MAXIE/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/MAXIEs/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/POINT1/EIN1,EUSE1
CAK   COMMON/IWATCH/MONITR,MAKEPLUS,MYUNRES
      COMMON/IWATCHs/MONITR,MAKEPLUS,MYUNRES
      COMMON/TEMPO/TEMP3,IVERSE
      COMMON/HIGHESTE/EHIGHEST,EMILLION
      COMMON/EXTEND/DTMAX,MESS
c----- 2012 - number of points read in any MT
      COMMON/READN2/N2READ
      INCLUDE 'sigma1.h'
      DATA XCKM1/0.0D+00/
      DATA EKM1/0.0D+00/
C-----THERMAL ENERGY - ALWAYS KEEP
      DATA ETHERMAL/2.53D-02/
C-----------------------------------------------------------------------
C
C     AT BEGINNING OF SECTION SET INDICES TO FORCE IMMEDIATE READ AND TO
C     LOAD UP TO THREE PAGES (UP TO 180000 POINTS) OF DATA. THIS IS NOT
C     THE BEGINNING OF THE SECTION SET INDICES TO ONLY LOAD THE THIRD
C     PAGE (60000 POINTS) OF DATA.
C
C-----------------------------------------------------------------------
C-----IS THIS THE BEGINNING OF THE SECTION.
      IF(HOT1.GT.1) GO TO 10
C-----BEGINNING OF SECTION. SET INDICES TO FORCE READING WHEN FIRST DATA
C-----POINT IS REQUESTED FROM ENDF/B DATA FILE.
      IPASS=0
      NFILL=1002
      IFILL=NFILL+1
C-----INITIALIZE INDICES TO READ UP TO THREE PAGES OF DATA POINTS.
      LOAD1=1
      LOAD2=NPT3
C-----INITIALIZE SIGN OF CROSS SECTION TO POSITIVE.
      SIGNT=1.0
C-----INITIALIZE INDEX TO FIRST POINT TO BROADEN.
      HOT3=0
      HOT3M1=-1
      GO TO 20
C-----SET INDEX TO ONLY LOAD LAST PAGE OF DATA POINTS.
   10 LOAD1=NP2P1
C-----INITIALIZE LOAD INDEX TO FIRST LOCATION.
   20 K=LOAD1
C-----IF NO MORE POINTS IN CORE READ NEXT PAGE (IF ANY).
   30 IF(IFILL.LE.NFILL) GO TO 50
      IF(N2LEFT.GT.NFILL) GO TO 40
      IF(N2LEFT.LE.0) GO TO 310
      NFILL=N2LEFT
   40 N2LEFT=N2LEFT-NFILL
      CALL POINTI(EIN,XCIN,NFILL)
C-----IF REQUESTED MAKE ALL NEGATIVE CROSS SECTIONS = 0
      IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,NFILL
      IF(XCIN(KP).LT.0.0D+00) XCIN(KP)=0.0D+00
      ENDDO
      ENDIF
      IFILL=1
C-----SAVE FIRST ENERGY OF SECTION.
      IF(IPASS.NE.0) GO TO 50
      IPASS=1
      EIN1=EIN(1)
C-----------------------------------------------------------------------
C
C     AT THE BEGINNING OF SECTION DECIDE WHETHER OR NOT TO DOPPLER
C     BROADEN. DO NOT DOPPLER BROADEN ANY SECTION IN WHICH THE
C     THRESHOLD IS HIGHER THAN 1,000,000* KT/A (I.E. SECTIONS WHICH
C     HAVE A HIGH ENERGY THRESHOLD, ABOVE WHICH BROADENING WILL HAVE
C     A NEGLIGABLE EFFECT).
C
C-----------------------------------------------------------------------
   50 IF(TOOHI.eq.0) go to 150
      IF(TOOHI.gt.0) go to 90
C-----SKIP POINTS BELOW THRESHOLD.
      DO 60 IFILL=1,NFILL
      IF(DABS(XCIN(IFILL)).GT.0.0) GO TO 70
   60 CONTINUE
C-----ENTIRE PAGE IS ZERO. SET INDEX AND READ NEXT PAGE.
      IFILL=NFILL+1
      GO TO 30
C-----DEFINE INDEX TO THRESHOLD.
   70 IFILL=IFILL-1
      IF(IFILL.LE.0) IFILL=1
C-----INITIALIZE PARAMETERS AND THEN DECIDE WHETHER OR NOT TO DOPPLER
C-----BROADEN.
      AX=ALPHA*EIN(IFILL)
      TOOHI=1
      ITHRES=0
C-----NO DOPPLER BROADENING IF THRESHOLD OVER 1,000,000*KT/A.
c-----2012 - AND LESS THAN 1,000 POINTS - ADDED FOR NEW (N,N')
C-----2013 - ADDED ABSOLUTE CUTOFF AT 10 MEV.
      IF(EIN(IFILL).gt.EHIGHEST) go to 80  ! less than 10 MeV
      IF(AX        .lt.EMILLION.or.        ! less than 1,000,000 AE/KT
     1   N2READ    .ge.1000)     go to 100 ! more than 1,000 points
C-----SECTION WILL ONLY BE COPIED. TURN OFF UNRESOLVED REGION FLAG.
   80 UREACT=0
      MESS  =3
C-----DATA WILL NOT BE DOPPLER BROADENED. MAXIMUM ENERGY SPACING IS
C-----NOT REQUIRED. JUST COPY POINTS DIRECTLY INTO HOT CROSS SECTION
C-----ARRAY.
   90 HOT3M1=HOT3
      HOT3=HOT3+1
      EHOT(HOT3)=EIN(IFILL)
      XCHOT(HOT3)=XCIN(IFILL)
      IFILL=IFILL+1
      GO TO 300
C-----SECTION WILL BE DOPPLER BROADENED.
  100 TOOHI=0
C-----IF FIRST POINT IS NEGATIVE (E.G. BACKGROUND CROSS SECTION) REVERSE
C-----SIGN OF CROSS SECTION.
      IF(XCIN(IFILL).LT.0.0) SIGNT=-1.0
C-----------------------------------------------------------------------
C
C     DEFINE FIRST DATA POINT. USE EITHER FIRST TABULATED POINT OR IF
C     THIS IS A THRESHOLD REACTION CALCULATE POSITION OF NEW LABORATORY
C     EFFECTIVE THRESHOLD DUE TO DOPPLER BROADENING AND USE THIS AS
C     FIRST POINT.
C
C-----------------------------------------------------------------------
      IF(EIN(IFILL).LE.ETHRES) GO TO 120
C-----------------------------------------------------------------------
C
C     DO NOT EXTEND CROSS SECTION FROM UNRESOLVED RESONANCE REGION.
C
C-----------------------------------------------------------------------
      IF(UREACT.LE.0) GO TO 110
      IF(EIN(IFILL).GE.EULOW) GO TO 120
C-----------------------------------------------------------------------
C
C     INSERT ONE OR TWO (IF FIRST POINT CROSS SECTION IS NOT ZERO)
C     POINTS BELOW THRESHOLD.
C
C-----------------------------------------------------------------------
  110 ITHRES=1
      V=DSQRT(AX)-ATOP
      ENEXT=V*V/ALPHA
      IF(ENEXT.LT.ETHRES.OR.V.LE.0.0) ENEXT=EMIN
      ECOLD(K)=ENEXT
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      XCCOLD(K)=0.0
      IF(K.EQ.LOAD1) EUSE1=ECOLD(K)
      K=K+1
      IF(XCIN(IFILL).eq.0.0D+0) go to 120
      ECOLD(K)=EIN(IFILL)
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      XCCOLD(K)=0.0
      IF(K.EQ.LOAD1) EUSE1=ECOLD(K)
      GO TO 130
C-----USE FIRST TABULATED POINT.
  120 ECOLD(K)=EIN(IFILL)
      XCCOLD(K)=XCIN(IFILL)
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      IFILL=IFILL+1
      IF(K.EQ.LOAD1) EUSE1=ECOLD(K)
C-----------------------------------------------------------------------
C
C     FIRST POINT DEFINED. RE-INITIALIZE UNRESOLVED RESONANCE REGION
C     INDICES IF SECTION STARTS ABOVE LOWER ENERGY OF UNRESOLVED REGION.
C
C-----------------------------------------------------------------------
  130 IF(UREACT.LE.0) GO TO 270
C-----IF SECTION STARTS ABOVE UPPER ENERGY LIMIT OF UNRESOLVED REGION
C-----TURN OFF UNRESOLVED REGION.
      IF(ECOLD(LOAD1).LT.EUHIGH) GO TO 140
      UREACT=0
      GO TO 270
C-----IF SECTION STARTS ABOVE LOWER ENERGY LIMIT OF UNRESOLVED REGION
C-----INDICATE THAT ADDITIONAL POINTS ARE NOT REQUIRED AT THE LOWER
C-----ENERGY LIMIT OF UNRESOLVED REGION.
  140 IF(ECOLD(LOAD1).LT.EULOW) GO TO 270
      IL=2
      UNRES1=1
      GO TO 270
C-----------------------------------------------------------------------
C
C     FIRST POINT DEFINED. LOAD ALL OTHER POINTS.
C
C     CHECK FOR UNRESOLVED RESONANCE REGION. IF THERE IS AN UNRESOLVED
C     RESONANCE REGION DEFINE INDICES TO THE LAST POINT BELOW THE
C     UNRESOLVED RESONANCE REGION (UNRES1) AND THE FIRST POINT ABOVE THE
C     UNRESOLVED RESONANCE REGION (UNRES2). INSURE THAT THERE ARE AT
C     LEAST TWO POINTS AT THE LOWER AND UPPER ENERGY LIMITS OF THE
C     UNRESOLVED RESONANCE REGION. WHEN THE LOWER LIMIT OF THE
C     UNRESOLVED REGION IS FOUND SET THE INDEX TO THE UPPER ENERGY LIMIT
C     TO A LARGE NUMBER TO INDICATE THAT THE UPPER ENERGY LIMIT IS STILL
C     BEYOND THE UPPER LIMIT OF THE POINTS LOADED INTO CORE SO FAR. ONCE
C     THE UPPER ENERGY LIMIT OF THE UNRESOLVED REGION IS FOUND THE INDEX
C     TO THE UPPER ENERGY LIMIT WILL BE PROPER DEFINED.
C
C-----------------------------------------------------------------------
  150 IF(UREACT.LE.0) GO TO 230
C-----INSURE THERE ARE AT LEAST TWO POINTS AT LOWER ENERGY LIMIT OF
C-----UNRESOLVED REGION.
      IF(EIN(IFILL).lt.EULOW) go to 230
      IF(EIN(IFILL).gt.EULOW) go to 160
      IF(IL.EQ.0) UNRES1=K
      IPATH=1
      XCULOW=XCIN(IFILL)
      GO TO 180
  160 IF(IL.eq.1) go to 170
      IF(IL.gt.1) go to 190
C-----THERE ARE NO POINTS AT LOWER ENERGY LIMIT. INTERPOLATE AND
C-----INSERT POINT.
      UNRES1=K
      XCULOW=((EIN(IFILL)-EULOW)*XCKM1+(EULOW-EKM1)*XCIN(IFILL))/
     1 (EIN(IFILL)-EKM1)
  170 IPATH=0
  180 ECOLD(K)=EULOW
      XCCOLD(K)=XCULOW
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      IL=IL+1
      GO TO 260
C-----INSURE THERE ARE AT LEAST TWO POINTS AT UPPER ENERGY LIMIT OF
C-----UNRESOLVED REGION.
  190 IF(EIN(IFILL).lt.EUHIGH) go to 240
      IF(EIN(IFILL).gt.EUHIGH) go to 200
      XCUHI=XCIN(IFILL)
      IPATH=1
      GO TO 220
  200 IF(IH.eq.1) go to 210
      IF(IH.gt.1) go to 230
C-----THERE ARE NO POINTS AT UPPER ENERGY LIMIT. INTERPOLATE AND
C-----INSERT FIRST POINT.
      XCUHI=((EIN(IFILL)-EUHIGH)*XCKM1+(EUHIGH-EKM1)*XCIN(IFILL))/
     1 (EIN(IFILL)-EKM1)
  210 IPATH=0
  220 ECOLD(K)=EUHIGH
      XCCOLD(K)=XCUHI
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      UNRES2=K
      IH=IH+1
      GO TO 260
C-----IF CROSS SECTION HAS CHANGED SIGN INTERPOLATE AND INSERT ENERGY
C-----POINT AT ZERO CROSS SECTION UNLESS CROSS SECTION AT LAST POINT
C-----WAS ZERO.
  230 IF(SIGNT*XCIN(IFILL).GE.0.0) GO TO 240
      SIGNT=-SIGNT
      IF(XCKM1.eq.0.0D+0) go to 240
      IPATH=0
      XCCOLD(K)=0.0
      ECOLD(K)=(XCIN(IFILL)*EKM1-XCKM1*EIN(IFILL))/(XCIN(IFILL)-XCKM1)
      IF(ECOLD(K).LT.EKM1) ECOLD(K)=EKM1
      IF(ECOLD(K).GT.EIN(IFILL)) ECOLD(K)=EIN(IFILL)
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      GO TO 260
C-----CHECK MAXIMUM ENERGY SPACING.
  240 IF(EIN(IFILL).LE.ENEXT) GO TO 250
      IPATH=0
      ECOLD(K)=ENEXT
      IF(ECOLD(K).LT.EKM1) ECOLD(K)=EKM1
      IF(ECOLD(K).GT.EIN(IFILL)) ECOLD(K)=EIN(IFILL)
C-----------------------------------------------------------------------
C
C     USE LOG-LOG INTERPOLATION AT LOW ENERGY
C
C-----------------------------------------------------------------------
      ITERP=2
      IF(ENEXT.LT.ELOGLOG) THEN
      IF(XCIN(IFILL).GT.0.0.AND.XCKM1.GT.0.0) ITERP = 5
      ENDIF
      XCCOLD(K)=TERPIT(ENEXT,EIN(IFILL),EKM1,XCIN(IFILL),XCKM1,ITERP)
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
      GO TO 260
C-----ENERGY SPACING IS O.K. ACCEPT NEXT TABULATED POINT.
  250 IPATH=1
      ECOLD(K)=EIN(IFILL)
      XCCOLD(K)=XCIN(IFILL)
      ARG=ALPHA*ECOLD(K)
      YCOLD(K)=DSQRT(ARG)
C-----NEXT INPUT POINT IS ACCEPTABLE.
  260 IFILL=IFILL+IPATH
C-----------------------------------------------------------------------
C
C     DEFINE SPEED-LIKE TERMS (SEE UCRL-50400, VOL. 17, PART C) AND
C     SLOPE BETWEEN POINTS.
C
C-----------------------------------------------------------------------
  270 IF(K.LE.1) GO TO 290
      DX=YCOLD(K)-YCOLD(K-1)
      DSCOLD=XCCOLD(K)-XCCOLD(K-1)
      IF(DX.eq.0.0D+0) go to 280
      DCOLD(K-1)=DSCOLD/(DX*(YCOLD(K)+YCOLD(K-1)))
      GO TO 290
  280 DCOLD(K-1)=0.0
C-----------------------------------------------------------------------
C
C     DEFINE NEXT ALLOWABLE ENERGY INTERVAL. IF THERE IS AN UNRESOLVED
C     RESONANCE REGION INSURE THAT THERE ARE AT LEAST TWO POINTS AT THE
C     LOWER AND UPPER ENERGY LIMITS OF THE UNRESOLVED REGION. DO NOT ADD
C     ENERGY POINTS WITHIN THE UNRESOLVED ENERGY REGION.
C
C-----------------------------------------------------------------------
C-----SAVE LAST POINT FOR INTERPOLATION.
  290 EKM1=ECOLD(K)
      XCKM1=XCCOLD(K)
  300 ENEXT=DEMAXC*ECOLD(K)
C-----RELAX ENERGY STEP ABOVE 1 EV.
      IF(EKM1.GE.1.0D+00) DEMAXC=DEMAXHI
C-----------------------------------------------------------------------
C
C     INSERT THERMAL POINT BY SETTING ENEXT = ETHERMAL
C     WHEN THERMAL ENERGY IS CROSSED.
C
C-----------------------------------------------------------------------
      IF(ECOLD(K).LT.ETHERMAL.AND.
     1   ENEXT   .GT.ETHERMAL) ENEXT=ETHERMAL
      K=K+1
      IF(K.LE.LOAD2) GO TO 30
C-----------------------------------------------------------------------
C
C     ALL POINTS REQUESTED, OR ALL REMAINING POINTS, HAVE BEEN LOADED.
C
C-----------------------------------------------------------------------
C-----IS THIS THE END OF THE DATA TABLE.
      IF(N2LEFT.LE.0.AND.IFILL.GT.NFILL) GO TO 320
C-----NO.
      GO TO 330
C-----END OF DATA TABLE. DEFINE NUMBER OF POINTS IN CORE AND SET INDEX
C-----TO BROADEN ALL REMAINING POINTS.
  310 LOAD2=K-1
  320 HOT2=LOAD2
      MTEND=-1
C-----DEFINE INDEX TO LAST DATA POINT IN CORE.
  330 COLD2=LOAD2
      DCOLD(COLD2)=0.0
      RETURN
      END
      SUBROUTINE BROADN
C=======================================================================
C
C     GIVEN A CROSS SECTION THAT IS DESCRIBED BY A TABLE OF CROSS
C     SECTION VS. E AND LINEAR-LINEAR INTERPOLATION BETWEEN POINTS
C     THIS ROUTINE WILL EXACTLY DOPPLER BROADENING THE CROSS
C     SECTION AT A PORTION OF THE ENERGIES (ONE, TWO OR THREE PAGES).
C
C     INPUT PARAMETERS
C     COLD1 =INDEX TO FIRST POINT LOADED IN CORE.
C     COLD2 =INDEX TO LAST POINT LOADED IN CORE.
C     COLD1P=INDEX TO FIRST POINT TO USE IN DOPPLER INTEGRATION
C     COLD2P=INDEX TO LAST POINT TO USE IN DOPPLER INTEGRATION
C            THE PAIRS (COLD1,COLD2) AND (COLD1P,COLD2P) WILL DIFFER
C            ONLY IF THERE IS AN UNRESOLVED RESONANCE PRESENT.
C     HOT1  =INDEX TO FIRST POINT TO BROADEN.
C     HOT2  =INDEX TO LAST POINT TO BROADEN.
C     YCOLD =TABLE OF SPEEDS CORRESPONDING TO XCCOLD.
C     XCCOLD =TABLE OF CROSS SECTIONS AT INITIAL TEMPERATURE.
C     DCOLD=SLOPE BETWEEN POINTS IN (YCOLD,XCCOLD) TABLE.
C     EHOT  =TABLE OF ENERGIES CORRESPONDING TO XCHOT (AND XCCOLD).
C     XCHOT =TABLE OF CROSS SECTIONS BROADENED TO TEMPK.
C     ALPHA =DOPPLER WIDTH (11505.3*AWR/TEMPEF)
C
C     THE TABLE OF POINTS (COLD1P,COLD2P) IS USED TO BROADEN THE
C     TABLE OF POINTS (HOT1,HOT2).
C
C     IF THERE ARE 180000 (3 PAGES) OR FEWER DATA POINTS ALL DOPPLER
C     BROADENED CROSS SECTIONS WILL BE CALCULATED IN ONE PASS THROUGH
C     THIS ROUTINE. IF THERE ARE OVER 180000 POINTS THE DATA WILL BE
C     DOPPLER BROADENED A PAGE AT A TIME. AN EXCEPTION IS THAT AT THE
C     BEGINNING OF THE SECTION THE FIRST TWO PAGES WILL BE BROADENED
C     IN ONE PASS AND AT THE END OF THE SECTION THE LAST TWO PAGES WILL
C     BE DOPPLER BROADENED.
C
C     NORMALLY ALL OF THE POINTS LOADED INTO CORE WILL BE USED IN THE
C     DOPPLER INTEGRATION, I.E. COLD1P=COLD1 AND COLD2P=COLD1P. HOWEVER
C     IF THE THREE PAGES OF DATA THAT ARE IN CORE CONTAIN ANY PORTION OF
C     THE UNRESOLVED RESONANCE REGION THE DOPPLER INTEGRATION WILL BE
C     CUT OFF AT THE EDGE OF THE UNRESOLVED REGION AND THE CROSS SECTION
C     WILL BE EXTENDED FROM THAT POINT AS 1/V. IN THIS CASE FOR DATA
C     POINTS AT ENERGIES LESS THAN THE UNRESOLVED RESONANCE REGION
C     COLD2P WILL POINT TO THE LOWER ENERGY LIMIT OF THE UNRESOLVED
C     REGION AND COLD1P WILL POINT TO THE FIRST POINT IN CORE. SIMILARLY
C     FOR POINTS ABOVE THE UPPER ENERGY LIMIT OF THE UNRESOLVED REGION
C     COLD1P WILL POINT TO THE UPPER ENERGY LIMIT OF THE UNRESOLVED
C     REGION AND COLD2P WILL POINT TO THE LAST POINT LOADED IN CORE.
C     WITHIN THE UNRESOLVED ENERGY RANGE THE ORIGINAL DATA POINT WILL
C     JUST BE COPIED WITHOUT DOPPLER BROADENING AND COLD1P AND COLD2P
C     ARE NOT USED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,
     1 HEATER,HIHEAT,UREVAL,UREACT,UNRES1,UNRES2,HEATM1,HEATM2,HOTEND
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
CAK   COMMON/PAGER/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/PAGERs/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/HOTS/ALPHA,HOTSY1,HOTSY2,TEMPK,TEMPEF,N2TAPI,N2TAPO
      COMMON/CONTAC/OVPI,OVPI2,ATOP
      COMMON/INDICE/HEATER,LOHEAT,HIHEAT
      COMMON/EXTEND/DTMAX,MESS
CAK   COMMON/MAXIE/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/MAXIEs/DEMAXC,DEMAXH,DEMAXLOW,DEMAXHI
      COMMON/RESOLV/EULOW,EUHIGH,UREVAL,UREACT,UNRES1,UNRES2,IL,IH
      COMMON/OKERRC/ERXC3C,ERXC30,KERR3C,MAXERC,ENER3C(21),ER3C(21)
      COMMON/SLIM/ISTART,NOTHIN,ITHIN1,ITHIN2,ITHIN3,MTEND
CAK   COMMON/PARAMS/XCMIN
      COMMON/PARAMSs/XCMIN
C-----08/20/2013 - added high energy cutoff
      COMMON/HIGHESTE/EHIGHEST,EMILLION
      INCLUDE 'sigma1.h'
      DATA SIGY/0.0D+00/
      DATA XCIHM1/0.0D+00/
      DATA E1MEV/1.0D+06/
      DATA ENEXT/0.0D+00/
      DATA HALF/5.0D-01/
      DATA EIHM1/0.0D+00/
C-----THERMAL ENERGY - ALWAYS KEEP
      DATA ETHERMAL/2.53D-02/
C-----01/16/04 - INITIALIZE TEMPERATURE DEPENDENT ENERGY SPACING
      ENEXT2 = 0.0
C-----------------------------------------------------------------------
C
C     TEST FOR CROSS SECTION EXTENSION
C
C-----------------------------------------------------------------------
C-----TEST LOWER LIMIT FOR ALL PAGES BEYOND SECOND.
      IF(HOT1.EQ.COLD1) GO TO 10
      DX=YCOLD(NP1P1)-YCOLD(1)
      IF(DX.GE.ATOP) GO TO 10
C-----DEFINE MAXIMUM ALLOWABLE TEMPERATURE CHANGE TO AVOID
C-----CROSS SECTION EXTENSION.
      DX=DX/ATOP
      DTMAX2=TEMPEF*DX*DX
      IF(DTMAX.LE.0.0.OR.DTMAX2.LT.DTMAX) DTMAX=DTMAX2
C-----SET FLAG TO PRINT OUT MESSAGE WITH SECTION
      MESS=2
C-----TEST UPPER LIMIT FOR ALL PAGES BUT LAST ONE.
   10 IF(HOT2.EQ.COLD2) GO TO 20
      DX=YCOLD(NPT3)-YCOLD(NPT2)
      IF(DX.GT.ATOP) GO TO 20
C-----DEFINE MAXIMUM ALLOWABLE TEMPERATURE CHANGE TO AVOID
C-----CROSS SECTION EXTENSION.
      DX=DX/ATOP
      DTMAX2=TEMPEF*DX*DX
      IF(DTMAX.LE.0.0.OR.DTMAX2.LT.DTMAX) DTMAX=DTMAX2
C-----SET FLAG TO PRINT OUT MESSAGE WITH SECTION
      MESS=2
C-----------------------------------------------------------------------
C
C     SET UP LOOP TO SELECT POINTS FROM WHICH TO BEGIN ITERATION.
C
C-----------------------------------------------------------------------
C-----TURN OFF BEGINNING OF ENERGY RANGE FLAG.
   20 HOTEND=0
      NLOW=HOT1
      IF(NLOW.GT.1) NLOW=NLOW-1
C-----INITIALIZE INDICES TO INTERIOR COLD POINTS IN ENERGY RANGE TO
C-----USE FOR INTEGRATION.
      LOHEAT=COLD1P+1
      HIHEAT=COLD2P-1
      DO 340 IHEAT=HOT1,HOT2
      IHEATM=IHEAT-1
C-----DEFINE COLD ENERGY AND CROSS SECTION IN SCALAR FORM.
      EIH=ECOLD(IHEAT)
      XCIH=XCCOLD(IHEAT)
C-----08/20/2013 - Copy High Energy Points (above 10 MeV).
      IF(EIH.GE.EHIGHEST) GO TO 30
C-----------------------------------------------------------------------
C
C     UNRESOLVED REGION WILL BE COPIED AND NOT BROADENED. WHEN THE
C     UPPER ENERGY LIMIT OF THE UNRESOLVED REGION IS REACHED RESET
C     INTEGRATION LIMITS TO EXTEND FROM UPPER ENERGY OF UNRESOLVED
C     REGION TO UPPER ENERGY OF DATA IN CORE.
C
C-----------------------------------------------------------------------
      IF(IHEAT.le.UNRES1) go to 50
      IF(IHEAT.eq.UNRES2) go to 40
      IF(IHEAT.gt.UNRES2) go to 50
C-----COPY POINTS WITHIN UNRESOLVED REGION.
   30 HOT3M1=HOT3
      HOT3=HOT3+1
      EHOT(HOT3)=EIH
      XCHOT(HOT3)=XCIH
C-----IF OVER 2 PAGES OF BROADENED DATA THIN IT.
      IF(HOT3.GT.NPT2) CALL THINIT
      GO TO 320
C-----UPPER LIMIT OF UNRESOLVED REGION REACHED. DEFINE RANGE OF
C-----OF INTEGRATION FROM UPPER LIMIT OF UNRESOLVED RANGE UP TO LAST
C-----POINT IN CORE.
   40 COLD1P=UNRES2
      COLD2P=COLD2
C-----RE-DEFINE INDICES TO INTERIOR COLD POINTS IN ENERGY RANGE TO
C-----USE FOR INTEGRATION.
      LOHEAT=COLD1P+1
      HIHEAT=COLD2P-1
C-----------------------------------------------------------------------
C
C     ALWAYS USE FIRST POINT AS NODE.
C
C-----------------------------------------------------------------------
C-----IS THIS FIRST POINT IN ENERGY RANGE.
   50 IF(IHEAT.NE.COLD1P) GO TO 60
C-----INITIALIZE SIGN OF SLOPE.
      DSIGNT=1.0
      IF(DCOLD(COLD1P).LT.0.0) DSIGNT=-1.0
C-----SET FLAG TO INDICATE BEGINNING OF ENERGY RANGE.
      HOTEND=1
      GO TO 100
C-----------------------------------------------------------------------
C
C     ALWAYS USE LAST POINT AS NODE.
C
C-----------------------------------------------------------------------
   60 IF(IHEAT.GE.HOT2.OR.IHEAT.GE.COLD2P) GO TO 100
C-----------------------------------------------------------------------
C
C     USE PRECEDING POINT IF IT WAS MAXIMUM OR MINIMUM AND ENERGY HAS
C     NOT YET BEEN USED.
C
C-----------------------------------------------------------------------
      IF(DSIGNT*DCOLD(IHEATM).GE.0.0) GO TO 70
      DSIGNT=-DSIGNT
C-----ONLY USE POINT IF SAME ENERGY HAS NOT YET BEEN USED.
      IF(EIHM1.le.EHOT(HOT3)) go to 330
      go to 90
C-----------------------------------------------------------------------
C
C     ONLY USE POINT IF SAME ENERGY HAS NOT YET BEEN USED (BROADENING
C     ELIMINATES ALL DISCONTINUITIES).
C
C-----------------------------------------------------------------------
   70 IF(EHOT(HOT3).GE.EIH) GO TO 330
C
C     USE ALL POINTS ABOVE 1 MEV.
C
      IF(EIH.GE.E1MEV) GO TO 100
C
C     USE ALL NON-POSITIVE CROSS SECTION POINTS.
C
      IF(XCIH.LE.0.0) GO TO 100
C
C     USE ALL POINTS WHERE COLD CROSS SECTION CHANGES BY MORE THAN A
C     FACTOR OF 1.5.
C
C-----01/16/04 - REDUCED FROM 2 TO 1.5
C-----           1.1 INCREASES RUNNING TIME WITHOUT IMPROVING RESULTS
      IF(XCIH.LE.1.5*SIGY.AND.SIGY.LE.1.5*XCIH) GO TO 80
C-----COLD CROSS SECTION HAS CHANGED BY FACTOR. USE EITHER CROSS
C-----SECTION OR ENERGY INTERVAL.
      IF(EIH.le.ENEXT) go to 100
      go to 110
C-----------------------------------------------------------------------
C
C     IF POINT IS NOT WITHIN ALLOWABLE ENERGY SPACING INTERPOLATE TO
C     INSERT ENERGY POINT.
C
C-----------------------------------------------------------------------
   80 IF(EIH.GT.ENEXT) GO TO 110
C
C     USE AT LEAST EVERY TENTH POINT.
C
C-----01/16/04 - TEMPERATURE DEEPENDENT ENERGY SPACING TEST.
      IF(EIH.GT.ENEXT2) GO TO 100
      IF((IHEAT-NLOW).GE.10) GO TO 100
C
C     DO NOT USE POINT AS NODE.
C
      GO TO 330
C-----------------------------------------------------------------------
C
C     DEFINE PARAMETERS AT ENERGY POINT TO BROADEN.
C
C-----------------------------------------------------------------------
C-----DEFINE BROADENED CROSS SECTION AT SAME ENERGY THAT COLD CROSS
C-----SECTION IS GIVEN (PRECEDING POINT).
   90 HOT3M1=HOT3
      HOT3=HOT3+1
      EHOT(HOT3)=EIHM1
      Y=YCOLD(IHEATM)
      SIGY=XCIHM1
      HEATM1=IHEATM
      HEATM2=IHEATM
      NHIGH=IHEATM
      GO TO 120
C-----DEFINE BROADENED CROSS SECTION AT SAME ENERGY THAT COLD CROSS
C-----SECTION IS GIVEN (CURRENT POINT).
  100 HOT3M1=HOT3
      HOT3=HOT3+1
      EHOT(HOT3)=EIH
      Y=YCOLD(IHEAT)
      SIGY=XCIH
      HEATM1=IHEAT
      HEATM2=IHEAT
      NHIGH=IHEAT
      GO TO 120
C-----DEFINE BROADENED CROSS SECTION AT ENERGY AT WHICH COLD CROSS
C-----SECTION NOT GIVEN (DEFINE CROSS SECTION BY ENERGY INTERPOLATIONS).
  110 HOT3M1=HOT3
      HOT3=HOT3+1
      EHOT(HOT3)=ENEXT
      ARG=ALPHA*EHOT(HOT3)
      Y=DSQRT(ARG)
      SIGY=((EIH-ENEXT)*XCIHM1+(ENEXT-EIHM1)*XCIH)/(EIH-EIHM1)
      HEATM1=IHEAT
      HEATM2=IHEATM
      NHIGH=IHEATM
C-----------------------------------------------------------------------
C
C     BROADEN DATA.
C
C-----------------------------------------------------------------------
C-----INITIALIZE COUNT OF SAVED DATA POINTS.
  120 ISAVE=0
C-----SELECT DOPPLER BROADENING METHOD.
      IF(Y.GT.HOTSY1) GO TO 130
      CALL BROADL(Y,SIGY,XCHOT(HOT3),HEATM1,HEATM2)
      GO TO 150
  130 IF(Y.GT.HOTSY2) GO TO 140
      CALL BROADH(Y,SIGY,XCHOT(HOT3),HEATM1,HEATM2)
      GO TO 150
  140 CALL BROADS(Y,SIGY,XCHOT(HOT3),HEATM1,HEATM2)
  150 CONTINUE
C-----CHECK FOR BEGINNING OF ENERGY RANGE (ONLY BROADEN ONE POINT AND
C-----THEN GO TO END OF LOOP TO SET UP FOR FIRST ENERGY INTERVAL).
      IF(HOTEND.GT.0) GO TO 320
C-----DO NOT INTERPOLATE ZERO LENGTH INTERVALS.
  160 IF(EHOT(HOT3).LE.EHOT(HOT3M1)) GO TO 300
C-----------------------------------------------------------------------
C
C     CHECK FOR CONVERGENCE.
C
C-----------------------------------------------------------------------
C-----DO NOT INTERPOLATE IF CROSS SECTIONS AT BOTH ENDS OF INTERVAL
C-----ARE ABSOLUTELY LESS THAN MINIMUM CROSS SECTION OF INTEREST.
      IF(DABS(XCHOT(HOT3)).LT.XCMIN.AND.DABS(XCHOT(HOT3M1)).LT.XCMIN)
     1 GO TO 300
C-----------------------------------------------------------------------
C
C     PERFORM CALCULATION AT MIDPOINT.
C
C-----------------------------------------------------------------------
C-----DEFINE MID-POINT ENERGY AND SPEED-LIKE TERM.
      EMID=HALF*(EHOT(HOT3)+EHOT(HOT3M1))
C-----ROUND.
      CALL INCORE9(EMID)
C-----CHECK FOR ZERO LENGTH INTERVALS.
      IF(EMID.LE.EHOT(HOT3M1).OR.EMID.GE.EHOT(HOT3)) GO TO 300
      ARG=ALPHA*EMID
      YMID=DSQRT(ARG)
C-----DEFINE INDICES TO ENERGY INTERVAL OR POINT.
      DO 170 NUSE=NLOW,IHEAT
      IF(YMID.lt.YCOLD(NUSE)) go to 190
      IF(YMID.eq.YCOLD(NUSE)) go to 180
  170 CONTINUE
      NUSE=IHEAT
C-----ENERGY POINT.
  180 SIGMID=XCCOLD(NUSE)
      HEATM1=NUSE
      HEATM2=NUSE
      GO TO 200
C-----ENERGY INTERVAL.
  190 HEATM1=NUSE
      HEATM2=NUSE-1
      WT1    = (ECOLD(NUSE)-EMID)/(ECOLD(NUSE)-ECOLD(HEATM2))
      WT2    = 1.0D+00 - WT1
      SIGMID = WT1*XCCOLD(HEATM2) + WT2*XCCOLD(NUSE)
  200 CONTINUE
C-----SELECT BROADENING METHOD.
      IF(YMID.GT.HOTSY1) GO TO 210
      CALL BROADL(YMID,SIGMID,XCMID,HEATM1,HEATM2)
      GO TO 230
  210 IF(YMID.GT.HOTSY2) GO TO 220
      CALL BROADH(YMID,SIGMID,XCMID,HEATM1,HEATM2)
      GO TO 230
  220 CALL BROADS(YMID,SIGMID,XCMID,HEATM1,HEATM2)
C-----------------------------------------------------------------------
C
C     CHECK FOR CONVERGENCE.
C
C-----------------------------------------------------------------------
C-----DO NOT FURTHER SUB-DIVIDE ENERGY RANGES WITH NEGATIVE CROSS
C-----SECTIONS (SUCH DATA HAS NO PHYSICAL MEANING AND ADDITIONAL
C-----TIME SHOULD NOT BE SPEND ON IT).
  230 IF(XCMID.LT.0.0) GO TO 290
C-----NO CONVERGENCE IF CROSS SECTION CHANGES BY MORE THAN 10 PER-CENT.
C-----01/16/04 - REDUCED FROM 1.4 TO 1.1
      IF(XCHOT(HOT3M1).GT.1.1*XCHOT(HOT3).OR.
     1 XCHOT(HOT3)    .GT.1.1*XCHOT(HOT3M1)) GO TO 280
C-----DEFINE CROSS SECTION AT MID-POINT BY LINEAR INTERPOLATION.
      WT1   = (EHOT(HOT3)-EMID)/(EHOT(HOT3)-EHOT(HOT3M1))
      WT2   = 1.0 - WT1
      XCLIN = WT1*XCHOT(HOT3M1) + WT2*XCHOT(HOT3)
C-----IF REQUESTED, DEFINE ENERGY DEPENDENT CONVERGENCE CRITERIA.
      IF(KERR3C.NE.0) CALL ERROKC(EMID)
C-----------------------------------------------------------------------
C
C     USE MORE STRINGENT CONVERGENCE CRITERIA IF ITERATING TOWARD.....
C
C-----------------------------------------------------------------------
      IF(XCMID.eq.XCHOT(HOT3)) go to 250
      IF(XCMID.gt.XCHOT(HOT3)) go to 240
C-----MINIMUM.
      IF(XCMID.lt.XCHOT(HOT3M1)) go to 260
      go to 250
C-----MAXIMUM.
  240 IF(XCMID.gt.XCHOT(HOT3M1)) go to 260
C-----STANDARD LINEAR INTERPOLATION.
  250 IF(DABS(XCMID-XCLIN).LE.DABS(ERXC3C*XCMID)) GO TO 290
      GO TO 270
C-----MORE STRINGENT LINEAR INTERPOLATION TOWARD MINIMUM/MAXIMUM.
  260 IF(DABS(XCMID-XCLIN).LE.DABS(ERXC30*XCMID)) GO TO 290
C-----LOW CROSS SECTION CONVERGENCE TEST.
  270 IF(DABS(XCMID).LT.XCMIN.AND.DABS(XCHOT(HOT3M1)).LT.XCMIN)
     1 GO TO 290
C-----NO CONVERGENCE. SAVE POINT AND HALF INTERVAL.
  280 IF(ISAVE.LT.MAXSAVE) ISAVE=ISAVE+1
      ESAVE(ISAVE)=EHOT(HOT3)
      XCSAVE(ISAVE)=XCHOT(HOT3)
      EHOT(HOT3)=EMID
      XCHOT(HOT3)=XCMID
      GO TO 160
C-----CONVERGENCE. IF THINNING WILL BE PERFORMED SAVE MIDPOINT AND
C-----END POINT. IF NO THINNING ONLY SAVE END POINT.
  290 IF(NOTHIN.GT.0) GO TO 300
      EHOT(HOT3+1)=EHOT(HOT3)
      XCHOT(HOT3+1)=XCHOT(HOT3)
      EHOT(HOT3)=EMID
      XCHOT(HOT3)=XCMID
      HOT3M1=HOT3
      HOT3=HOT3+1
C-----CONVERGENCE. IF OVER 2 PAGES OF BROADENED DATA THIN IT.
  300 IF(HOT3.GT.NPT2) CALL THINIT
C-----IF END OF INTERVAL PROCEED TO NEXT INTERVAL. OTHERWISE USE LAST
C-----POINT GENERATED.
      IF(ISAVE.LE.0) GO TO 310
      HOT3M1=HOT3
      HOT3=HOT3+1
      EHOT(HOT3)=ESAVE(ISAVE)
      XCHOT(HOT3)=XCSAVE(ISAVE)
      ISAVE=ISAVE-1
      GO TO 160
C-----END OF INTERVAL. SET UP FOR NEXT INTERVAL BY PLACING END OF LAST
C-----INTERVAL AT BEGINNING OF NEXT INTERVAL.
  310 NLOW=HEATM1
      ENEXT=DEMAXH*EHOT(HOT3)
C-----01/16/04 - DEFINE SPACING FOR TEMPERATURE DEPENDENT TEST.
      ENEXT2 = DSQRT(ALPHA*EHOT(HOT3)) + ATOP
      ENEXT2 = ENEXT2**2/ALPHA
C-----RELAX ENERGY STEP ABOVE 1 EV.
      IF(EHOT(HOT3).GE.1.0D+00) DEMAXH=DEMAXHI
C-----------------------------------------------------------------------
C
C     INSERT THERMAL POINT BY SETTING ENEXT = ETHERMAL
C     WHEN THERMAL ENERGY IS CROSSED.
C
C-----------------------------------------------------------------------
      IF(EHOT(HOT3).LT.ETHERMAL.AND.
     1   ENEXT     .GT.ETHERMAL) ENEXT=ETHERMAL
C-----CONTINUE IN LOOP IF PRECEDING (NOT CURRENT) POINT WAS USED OR
C-----IF STILL ITERATING IN ENERGY RANGE DUE TO ENERGY POINT SPACING.
      IF(NHIGH.NE.IHEAT) GO TO 70
C-----DEFINE BEGINNING OF NEXT TABULATED ENERGY INTERVAL.
  320 HOTEND=0
      NLOW=IHEAT
      SIGY=XCIH
      ENEXT=DEMAXH*EHOT(HOT3)
C-----RELAX ENERGY STEP ABOVE 1 EV.
      IF(EHOT(HOT3).GE.1.0D+00) DEMAXH=DEMAXHI
C-----------------------------------------------------------------------
C
C     INSERT THERMAL POINT BY SETTING ENEXT = ETHERMAL
C     WHEN THERMAL ENERGY IS CROSSED.
C
C-----------------------------------------------------------------------
      IF(EHOT(HOT3).LT.ETHERMAL.AND.
     1   ENEXT     .GT.ETHERMAL) ENEXT=ETHERMAL
C-----SAVE LAST COLD ENERGY AND CROSS SECTION.
  330 EIHM1=EIH
      XCIHM1=XCIH
  340 CONTINUE
      RETURN
      END
      SUBROUTINE BROADL(Y,SIGY,XCY,HEATM1,HEATM2)
C=======================================================================
C
C     LOW ENERGY DOPPLER BROADENING ROUTINE. THIS ROUTINE WILL
C     BE USED TO DOPPLER BROADEN ALL CROSS SECTIONS AT ENERGIES
C     WHERE AE/KT IS LESS THAN OR EQUAL TO 16. ANY POINT WITH A
C     HIGHER ENERGY WILL BE PASSED ON TO ROUTINE BROADH.
C
C     FOR AE/KT LESS THAN OR EQUAL TO 16 BOTH EXPONENTIALS IN THE
C     DOPPLER BROADENING KERNEL MUST BE CONSIDERED. FOR HIGHER ENERGIES
C     THE SECOND EXPONENTIAL MAY BE IGNORED, WHICH SIMPLIFIES THE
C     DOPPLER BROADENING.
C
C     THE ROUTINE HAS BEEN DESIGNED WITH NO SUBROUTINE CALLS
C     IN ORDER TO MINIMIZE RUNNING TIME. THE ARITHMETIC
C     STATEMENT FUNCTIONS RATION(A) AND ERFC(R,EXPERF) WILL
C     BE COMPILED AS IN LINE CODING BY VIRTUALLY ANY FORTRAN
C     COMPILER, AND AS SUCH DO NOT REPRESENT FUNCTION CALLS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HEATER,HIHEAT,
     1 HEATM1,HEATM2,HOT3,HOT3M1
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
      COMMON/CONTAC/OVPI,OVPI2,ATOP
      COMMON/INDICE/HEATER,LOHEAT,HIHEAT
      INCLUDE 'sigma1.h'
      DATA ZERO/0.0D+00/
      DATA HALF/0.5D+00/
      DATA ONE/1.0D+00/
      DATA TWO/2.0D+00/
      DATA FOUR/4.0D+00/
C-----------------------------------------------------------------------
C
C     SET UP LOOP TO DOPPLER BROADEN CROSS SECTIONS.
C
C-----------------------------------------------------------------------
C-----DEFINE ALL REQUIRED CONSTANTS FOR POINT.
      Y2=Y+Y
      YY=Y*Y
      YYA=YY+HALF
C-----INITIALIZE INTEGRALS.
      DUMMY1=ZERO
      IF(Y2.LE.ATOP) GO TO 10
      EXPY2=ZERO
      ERFCY2=ZERO
      GO TO 20
   10 EXPY2=DEXP(-Y2*Y2)
      IX2 = 1.0D+03*Y2
      Z = Y2 - ERFCX(IX2)
      ERFCY2 = EXPY2*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
C-----------------------------------------------------------------------
C
C     INTEGRATE OVER ALL ENERGY INTERVALS ABOVE CURRENT ENERGY POINT.
C
C-----------------------------------------------------------------------
C-----INITIALIZE
   20 F2A=YYA*(ONE-ERFCY2)+OVPI2*Y
      SIGA=SIGY
C-----SET UP LOOP OVER INTERVALS ABOVE CURRENT ONE + 1 (FOR EXTENSION)
      DO 60 HEATER=HEATM2,HIHEAT+1
C-----EXTEND AS CONSTANT
      IF(HEATER.GT.HIHEAT) THEN
      B=Y+ATOP
      BMY=ATOP
      SIGB=SIGA
      ELSE
C-----USE TABULATED INTERVAL
      B=YCOLD(HEATER+1)
      SIGB=XCCOLD(HEATER+1)
C-----ONLY EXTEND RANGE OF INTEGRATION UP TO ATOP UNITS ABOVE Y.
      BMY=B-Y
      IF(BMY.LE.ATOP) GO TO 30
      XXA=YCOLD(HEATER)**2
      XXB=B*B
      BMY=ATOP
      B=Y+ATOP
      BS=B*B
      SIGB=((BS-XXB)*SIGA+(XXA-BS)*SIGB)/(XXA-XXB)
      ENDIF
C-----DEFINE CONTRIBUTION OF FIRST INTEGRAL.
   30 EXPBM=DEXP(-BMY*BMY)
      IX2 = 1.0D+03*BMY
      Z = BMY - ERFCX(IX2)
      ERFCBM = EXPBM*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
C-----DEFINE CONTRIBUTION OF SECOND INTEGRAL.
      BPY=B+Y
      IF(BPY.LE.ATOP) GO TO 40
      EXPBP=ZERO
      ERFCBP=ZERO
      F2B=YYA*ERFCBM+OVPI*BPY*EXPBM
      GO TO 50
   40 EXPBP=DEXP(-BPY*BPY)
      IX2 = 1.0D+03*BPY
      Z = BPY - ERFCX(IX2)
      ERFCBP = EXPBP*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
C-----ADD CONTRIBUTION FROM CURRENT INTERVAL.
      F2B=YYA*(ERFCBM-ERFCBP)+OVPI*(B*(EXPBM-EXPBP)+Y*(EXPBM+EXPBP))
   50 DUMMY1=DUMMY1+(SIGA+SIGB)*(F2A-F2B)
C-----TEST FOR END OF RANGE OF INTEGRATION.
      IF(BMY.GE.ATOP) GO TO 70
C-----SAVE VALUES FROM LAST INTEGRAL.
      F2A=F2B
   60 SIGA=SIGB
C-----------------------------------------------------------------------
C
C     INTEGRATE OVER ALL ENERGY INTERVALS BELOW CURRENT ENERGY POINT.
C
C-----------------------------------------------------------------------
C-----NO INTERVALS BELOW FIRST POINT.
   70 IF(HEATM1.LE.COLD1P) GO TO 120
C-----RE-INITIALIZE INTEGRALS TO ZERO DISTANCE VALUES
      F2B=YYA*(ONE+ERFCY2)-OVPI2*Y
      SIGB=SIGY
C-----SET UP LOOP OVER INTERVALS BELOW CURRENT POINT.
      HEATER=HEATM1
      DO 110 LL=LOHEAT,HEATM1
      HEATER=HEATER-1
      A=YCOLD(HEATER)
      SIGA=XCCOLD(HEATER)
C-----ONLY EXTEND RANGE OF INTEGRATION DOWN TO ATOP UNITS BELOW Y.
      YMA=Y-A
      IF(YMA.LE.ATOP) GO TO 80
      XXB=YCOLD(HEATER+1)**2
      XXA=A*A
      YMA=ATOP
      A=Y-ATOP
      AS=A*A
      SIGA=((AS-XXB)*SIGA+(XXA-AS)*SIGB)/(XXA-XXB)
C-----DEFINE CONTRIBUTION OF FIRST INTEGRAL.
   80 EXPAM=DEXP(-YMA*YMA)
      IX2 = 1.0D+03*YMA
      Z = YMA - ERFCX(IX2)
      ERFCAM = EXPAM*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
C-----DEFINE CONTRIBUTION OF SECOND INTEGRAL.
      YPA=Y+A
      IF(YPA.LE.ATOP) GO TO 90
      EXPAP=ZERO
      ERFCAP=ZERO
      F2A=YYA*ERFCAM-OVPI*YPA*EXPAM
      GO TO 100
   90 EXPAP=DEXP(-YPA*YPA)
      IX2 = 1.0D+03*YPA
      Z = YPA - ERFCX(IX2)
      ERFCAP = EXPAP*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
C-----ADD CONTRIBUTION FROM CURRENT INTERVAL.
      F2A=YYA*(ERFCAM+ERFCAP)-OVPI*(A*(EXPAM-EXPAP)+Y*(EXPAM+EXPAP))
  100 DUMMY1=DUMMY1+(SIGA+SIGB)*(F2B-F2A)
C-----TEST FOR END OF RANGE OF INTEGRATION.
      IF(YMA.GE.ATOP) GO TO 130
C-----SAVE VALUES FROM LAST INTEGRAL.
      F2B=F2A
  110 SIGB=SIGA
C
C     CROSS SECTION EXTENSION AS 1/V BELOW TABULATED ENERGY RANGE.
C
C-----CONTINUE CROSS SECTION AS 1/V FROM 0 TO GENERAL LOWER LIMIT
      ADDER =       XCCOLD(COLD1P)*YCOLD(COLD1P)*TWO*
     1 (Y*(ERFCAM-ERFCAP)-OVPI*(EXPAM-EXPAP))
      IF(ADDER.GT.0.0D+00) DUMMY1 = DUMMY1 + ADDER
      GO TO 130
C-----Y IS AT LOWER END OF TABULATED RANGE. CONTINUE CROSS SECTION AS
C-----1/V FROM 0 TO LOWER LIMIT Y.
C-----CONTINUE CROSS SECTION AS 1/V TO YCOLD=0.0
  120 ADDER =       XCCOLD(COLD1P)*YCOLD(COLD1P)*TWO*
     1 (Y*(ONE-ERFCY2)-OVPI*(ONE-EXPY2))
      IF(ADDER.GT.0.0D+00) DUMMY1 = DUMMY1 + ADDER
C-----------------------------------------------------------------------
C
C     DEFINE BROADENED CROSS SECTION.
C
C-----------------------------------------------------------------------
  130 XCY=DUMMY1/(FOUR*YY)
      RETURN
      END
      SUBROUTINE BROADH(Y,SIGY,XCY,HEATM1,HEATM2)
C=======================================================================
C
C     HIGH ENERGY DOPPLER BROADENING ROUTINE. THIS ROUTINE WILL
C     BE USED TO DOPPLER BROADEN ALL CROSS SECTIONS AT ENERGIES
C     WHERE AE/KT IS GREATER THAN 16. ANY POINT WITH A LOWER
C     ENERGY WILL ALREADY HAVE BEEN DOPPLER BROADENED BY
C     ROUTINE BROADL.
C
C     FOR AE/KT LESS THAN OR EQUAL TO 16 BOTH EXPONENTIALS IN THE
C     DOPPLER BROADENING KERNEL MUST BE CONSIDERED. FOR HIGHER ENERGIES
C     THE SECOND EXPONENTIAL MAY BE IGNORED, WHICH SIMPLIFIES THE
C     DOPPLER BROADENING.
C
C     THE ROUTINE HAS BEEN DESIGNED WITH NO SUBROUTINE CALLS
C     IN ORDER TO MINIMIZE RUNNING TIME. THE ARITHMETIC
C     STATEMENT FUNCTIONS RATION(A) AND ERFC(R,EXPERF) WILL
C     BE COMPILED AS IN LINE CODING BY VIRTUALLY ANY FORTRAN
C     COMPILER, AND AS SUCH DO NOT REPRESENT FUNCTION CALLS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HEATER,HIHEAT,
     1 HEATM1,HEATM2,HOT3,HOT3M1
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
      COMMON/CONTAC/OVPI,OVPI2,ATOP
      COMMON/INDICE/HEATER,LOHEAT,HIHEAT
      INCLUDE 'sigma1.h'
      DATA ZERO/0.0D+00/
      DATA HALF/0.5D+00/
      DATA QUART3/0.75D+00/
      DATA ONEP5/1.5D+00/
      DATA TWO/2.0D+00/
      DATA TWOP5/2.5D+00/
C-----------------------------------------------------------------------
C
C     SET UP LOOP TO DOPPLER BROADEN CROSS SECTIONS.
C
C-----------------------------------------------------------------------
C-----DEFINE ALL REQUIRED CONSTANTS FOR POINT.
      YY=Y*Y
      YYA=YY+HALF
      YYB=TWOP5*YY+QUART3
      YYC=OVPI2*Y
      YYD=YYC*(YY+TWO)
C-----INITIALIZE INTEGRALS.
      DUMMY1=ZERO
      DUMMY2=ZERO
C-----------------------------------------------------------------------
C
C     INTEGRATE OVER ALL ENERGY INTERVALS ABOVE CURRENT ENERGY POINT.
C
C-----------------------------------------------------------------------
C-----INITIALIZE INTEGRALS.
      A=Y
      SIGA=SIGY
      F2A=YYA+YYC
      F4A=YYB+YYD
C-----SET UP LOOP OVER INTERVALS ABOVE CURRENT ONE + 1 (FOR EXTENSION).
      DO 20 HEATER=HEATM2,HIHEAT+1
      IF(HEATER.GT.HIHEAT) THEN
      BMY=ATOP
      B=Y+ATOP
      ELSE
      B=YCOLD(HEATER+1)
      BMY=B-Y
C-----ONLY EXTEND RANGE OF INTEGRATION UP TO ATOP UNITS ABOVE Y.
      IF(BMY.LE.ATOP) GO TO 10
      BMY=ATOP
      B=Y+ATOP
      ENDIF
C-----DEFINE CONSTANTS FOR THIS POINT.
   10 BB1P5=B*B+ONEP5
      BPY=B+Y
      EXPBM=DEXP(-BMY*BMY)
      IX2 = 1.0D+03*BMY
      Z = BMY - ERFCX(IX2)
      ERFCBM = EXPBM*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
      F2B=YYA*ERFCBM+OVPI*BPY*EXPBM
      YYMAA=(Y-A)*(Y+A)
C-----ADD CONTRIBUTION FROM CURRENT INTERVAL.
      DUMMY1=DUMMY1+SIGA*(F2A-F2B)
      IF(HEATER.LE.HIHEAT) THEN
      F4B=(YYA*YYMAA+YYB)*ERFCBM+OVPI*(BPY*(BB1P5+YYMAA)+Y)*EXPBM
      DUMMY2=DUMMY2+DCOLD(HEATER)*(F4A-F4B)
      ENDIF
C-----TEST FOR END OF RANGE OF INTEGRATION.
      IF(BMY.GE.ATOP) GO TO 30
C-----SAVE VALUES FROM LAST INTERGAL.
      A=B
      SIGA=XCCOLD(HEATER+1)
      F2A=F2B
      YYMAA=(Y-A)*(Y+A)
   20 F4A=(YYA*YYMAA+YYB)*ERFCBM+OVPI*(BPY*(BB1P5+YYMAA)+Y)*EXPBM
C-----------------------------------------------------------------------
C
C     INTEGRATE OVER ALL ENERGY INTERVALS BELOW CURRENT ENERGY POINT.
C
C-----------------------------------------------------------------------
C-----NO INTERVALS BELOW CURRENT FIRST POINT.
   30 IF(HEATM1.LE.COLD1P) GO TO 60
C-----RE-INITIALIZE INTEGRALS TO ZERO DISTANCE VALUES
      B=Y
      SIGB=SIGY
      F2B=YYA-YYC
      F4B=YYB-YYD
C-----SET UP LOOP OVER INTERVALS BELOW CURRENT POINT.
      HEATER=HEATM1
      DO 50 LL=LOHEAT,HEATM1
      HEATER=HEATER-1
      A=YCOLD(HEATER)
      YMA=Y-A
C-----ONLY EXTEND RANGE OF INTEGRATION DOWN TO ATOP UNITS BELOW Y.
      IF(YMA.LE.ATOP) GO TO 40
      YMA=ATOP
      A=Y-ATOP
C-----DEFINE CONSTANTS FOR THIS POINT.
   40 AA1P5=A*A+ONEP5
      YPA=Y+A
      EXPAM=DEXP(-YMA*YMA)
      IX2 = 1.0D+03*YMA
      Z = YMA - ERFCX(IX2)
      ERFCAM = EXPAM*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
      F2A=YYA*ERFCAM-OVPI*YPA*EXPAM
      YYMBB=(Y-B)*(Y+B)
      F4A=(YYA*YYMBB+YYB)*ERFCAM-OVPI*(YPA*(AA1P5+YYMBB)+Y)*EXPAM
C-----ADD CONTRIBUTION FROM CURRENT INTERVAL.
      DUMMY1=DUMMY1+SIGB*(F2B-F2A)
      DUMMY2=DUMMY2+DCOLD(HEATER)*(F4B-F4A)
C-----TEST FOR END OF RANGE OF INTEGRATION.
      IF(YMA.GE.ATOP) GO TO 70
C-----SAVE VALUES FROM LAST INTEGRAL.
      B=A
      SIGB=XCCOLD(HEATER)
      F2B=F2A
      YYMBB=(Y-B)*(Y+B)
   50 F4B=(YYA*YYMBB+YYB)*ERFCAM-OVPI*(YPA*(AA1P5+YYMBB)+Y)*EXPAM
C
C     EXTEND CROSS SECTION AS 1/V BELOW TABULATED ENERGY RANGE
C
C-----CONTINUE CROSS SECTION AS 1/V FROM 0 TO GENERAL LOWER LIMIT
      ADDER =       XCCOLD(COLD1P)*YCOLD(COLD1P)*(Y*ERFCAM-OVPI*EXPAM)
      IF(ADDER.GT.0.0D+00) DUMMY1 = DUMMY1 + ADDER
      GO TO 70
C-----Y IS AT LOWER END OF TABULATED RANGE. CONTINUE CROSS SECTION AS
C-----1/V FROM 0 TO LOWER LIMIT Y
   60 ADDER =       XCCOLD(COLD1P)*YCOLD(COLD1P)*(Y-OVPI)
      IF(ADDER.GT.0.0D+00) DUMMY1 = DUMMY1 + ADDER
C-----------------------------------------------------------------------
C
C     DEFINE BROADENED CROSS SECTION.
C
C-----------------------------------------------------------------------
   70 XCY=HALF*(DUMMY1+DUMMY2)/YY
      RETURN
      END
      SUBROUTINE BROADS(Y,SIGY,XCY,HEATM1,HEATM2)
C=======================================================================
C
C     HIGH ENERGY DOPPLER BROADENING ROUTINE. THIS ROUTINE WILL
C     BE USED TO DOPPLER BROADEN ALL CROSS SECTIONS AT ENERGIES
C     WHERE AE/KT IS GREATER THAN ONE MILLION. ANY POINT WITH A LOWER
C     ENERGY WILL ALREADY HAVE BEEN DOPPLER BROADENED BY ROUTINE
C     BROADL OR BROADH.
C
C     FOR AE/KT LESS THAN OR EQUAL TO 16 BOTH EXPONENTIALS IN THE
C     DOPPLER BROADENING KERNEL MUST BE CONSIDERED. FOR HIGHER ENERGIES
C     THE SECOND EXPONENTIAL MAY BE IGNORED, WHICH SIMPLIFIES THE
C     DOPPLER BROADENING.
C
C     FOR AE/KT GREATER THAN 1,000,000 IT IS POSSIBLE TO ASSUME THE
C     TERM (X/Y)**2 IN THE DOPPLER INTEGRAL IS JUST UNITY (SINCE X ONLY
C     VARIES FROM Y-4 TO Y+4, FOR Y GREATER THAN 1000 THIS TERM IS
C     ESSENTIALLY CONSTANT SINCE IT VARIES FROM (996/1000)**2 TO
C     (1004/1000)**2.= 1.0 +/- 0.008.
C
C     THE ROUTINE HAS BEEN DESIGNED WITH NO SUBROUTINE CALLS
C     IN ORDER TO MINIMIZE RUNNING TIME. THE ARITHMETIC
C     STATEMENT FUNCTIONS RATION(A) AND ERFC(R,EXPERF) WILL
C     BE COMPILED AS IN LINE CODING BY VIRTUALLY ANY FORTRAN
C     COMPILER, AND AS SUCH DO NOT REPRESENT FUNCTION CALLS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HEATER,HIHEAT,
     1 HEATM1,HEATM2,HOT3,HOT3M1
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
      COMMON/CONTAC/OVPI,OVPI2,ATOP
      COMMON/INDICE/HEATER,LOHEAT,HIHEAT
      INCLUDE 'sigma1.h'
      DATA ZERO/0.0D+00/
      DATA HALF/0.5D+00/
      DATA ONE/1.0D+00/
      DATA TWOP5/2.5D+00/
C-----------------------------------------------------------------------
C
C     SET UP LOOP TO DOPPLER BROADEN CROSS SECTIONS.
C
C-----------------------------------------------------------------------
C-----DEFINE ALL REQUIRED CONSTANTS FOR POINT.
      YY = Y*Y
      YYA=OVPI2*Y
C-----INITIALIZE INTEGRALS.
      DUMMY1=ZERO
      DUMMY2=ZERO
C-----------------------------------------------------------------------
C
C     INTEGRATE OVER ALL ENERGY INTERVALS ABOVE CURRENT ENERGY POINT.
C
C-----------------------------------------------------------------------
C-----INITIALIZE INTEGRALS.
      A=Y
      SIGA=SIGY
      ERFCAM=ONE
      F4A=TWOP5+YYA
C-----SET UP LOOP OVER INTERVALS ABOVE CURRENT ONE + 1 (FOR EXTENSION)
      DO 20 HEATER=HEATM2,HIHEAT+1
      IF(HEATER.GT.HIHEAT) THEN
C-----EXTEND AS CONSATANT
      BMY=ATOP
      B=Y+ATOP
      ELSE
C-----USE TABULATED INTERVAL
      B=YCOLD(HEATER+1)
      BMY=B-Y
C-----ONLY EXTEND RANGE OF INTEGRATION UP TO ATOP UNITS ABOVE Y.
      IF(BMY.LE.ATOP) GO TO 10
      BMY=ATOP
      B=Y+ATOP
      ENDIF
C-----DEFINE CONSTANTS FOR THIS POINT.
   10 EXPBM=DEXP(-BMY*BMY)
      IX2 = 1.0D+03*BMY
      Z = BMY - ERFCX(IX2)
      ERFCBM = EXPBM*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
C-----ADD CONTRIBUTION FROM CURRENT INTERVAL.
      F2B=OVPI*(B+Y)*EXPBM
      DUMMY1=DUMMY1+SIGA*(ERFCAM-ERFCBM)
      IF(HEATER.LE.HIHEAT) THEN
      F4B=((Y-A)*(Y+A)+TWOP5)*ERFCBM+F2B
      DUMMY2=DUMMY2+DCOLD(HEATER)*(F4A-F4B)
      ENDIF
C-----TEST FOR END OF RANGE OF INTEGRATION.
      IF(BMY.GE.ATOP) GO TO 30
C-----SAVE VALUES FROM LAST INTERGAL.
      A=B
      SIGA=XCCOLD(HEATER+1)
      ERFCAM=ERFCBM
   20 F4A=((Y-A)*(Y+A)+TWOP5)*ERFCBM+F2B
C-----------------------------------------------------------------------
C
C     INTEGRATE OVER ALL ENERGY INTERVALS BELOW CURRENT ENERGY POINT.
C
C-----------------------------------------------------------------------
C-----NO INTERVALS BELOW CURRENT FIRST POINT.
   30 IF(HEATM1.LE.COLD1P) GO TO 60
C-----RE-INITIALIZE INTEGRALS TO ZERO DISTANCE VALUES
      B=Y
      SIGB=SIGY
      ERFCBM=ONE
      F4B=TWOP5-YYA
C-----SET UP LOOP OVER INTERVALS BELOW CURRENT POINT.
      HEATER=HEATM1
      DO 50 LL=LOHEAT,HEATM1
      HEATER=HEATER-1
      A=YCOLD(HEATER)
      YMA=Y-A
C-----ONLY EXTEND RANGE OF INTEGRATION DOWN TO ATOP UNITS BELOW Y.
      IF(YMA.LE.ATOP) GO TO 40
      YMA=ATOP
      A=Y-ATOP
C-----DEFINE CONSTANTS FOR THIS POINT.
   40 EXPAM=DEXP(-YMA*YMA)
      IX2 = 1.0D+03*YMA
      Z = YMA - ERFCX(IX2)
      ERFCAM = EXPAM*(((F3(IX2)*Z+F2(IX2))*Z+F1(IX2))*Z+F0(IX2))
      F2A=OVPI*(Y+A)*EXPAM
      F4A=((Y-B)*(Y+B)+TWOP5)*ERFCAM-F2A
C-----ADD CONTRIBUTION FROM CURRENT INTERVAL.
      DUMMY1=DUMMY1+SIGB*(ERFCBM-ERFCAM)
      DUMMY2=DUMMY2+DCOLD(HEATER)*(F4B-F4A)
C-----TEST FOR END OF RANGE OF INTEGRATION.
      IF(YMA.GE.ATOP) GO TO 70
C-----SAVE VALUES FROM LAST INTEGRAL.
      B=A
      SIGB=XCCOLD(HEATER)
      ERFCBM=ERFCAM
   50 F4B=((Y-B)*(Y+B)+TWOP5)*ERFCAM-F2A
C
C     CONTINUE CROSS SECTION AT 1/V BELOW TABULATED ENERGY RANGE.
C
C-----CONTINUE CROSS SECTION AS 1/V FROM 0 TO GENERAL LOWER LIMIT
      ADDER = XCCOLD(COLD1P)*YCOLD(COLD1P)*(Y*ERFCAM-OVPI*EXPAM)/YY
      IF(ADDER.GT.0.0D+00) DUMMY1 = DUMMY1 + ADDER
      GO TO 70
C-----Y IS AT LOWER END OF TABULATED RANGE. CONTINUE CROSS SECTION AS
C-----1/V FROM 0 TO LOWER LIMIT Y
   60 ADDER =       XCCOLD(COLD1P)*YCOLD(COLD1P)*(Y-OVPI)/YY
      IF(ADDER.GT.0.0D+00) DUMMY1 = DUMMY1 + ADDER
C-----------------------------------------------------------------------
C
C     DEFINE BROADENED CROSS SECTION.
C
C-----------------------------------------------------------------------
   70 XCY=HALF*(DUMMY1+DUMMY2)
      RETURN
      END
      SUBROUTINE THINIT
C=======================================================================
C
C     GIVEN A FUNCTION REPRESENTED BY A TABLE OF ENERGIES EHOT AND
C     CROSS SECTIONS XCHOT AND LINEAR-LINEAR INTERPOLATION BETWEEN
C     TABULATED POINTS, THIS ROUTINE WILL REMOVE ALL EXTRANEOUS
C     POINTS THAT LIE WITHIN A GIVEN ACCURACY OF THE FUNCTION BASED
C     UPON INTERPOLATION FROM THE SURROUNDING POINTS.
C
C     DURING ONE PASS THIS ROUTINE WILL THIN ANY UNCONVERGED POINTS
C     FROM THE LAST BROADENED PAGE PLUS THE CURRENT BROADENED PAGE.
C     AT THE END OF THE CURRENT PAGE THE LAST CONVERGED POINT AND ALL
C     UNCONVERGED POINTS WILL BE LEFT IN THE CORE TO BE THINNED WITH
C     THE NEXT PAGE (UNLESS THE END OF THE TABLE IS REACHED IN WHICH
C     CASE THE LAST POINT IS CONVERGED). IF THERE IS MORE THAN A FULL
C     PAGE OF THINNED POINTS IT WILL BE COPIED TO SCRATCH AND ALL
C     REMAINING POINTS WILL BE SHIFTED FORWARD ONE PAGE IN THE DOPPLER
C     BROADENED DATA PAGES.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,SCR,COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,
     1 HOT3M1,UREVAL,UREACT,UNRES1,UNRES2,TOOHI
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/UNITS/SCR
      COMMON/UNITSs/SCR
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
CAK   COMMON/PAGER/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/PAGERs/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/OKERRT/ERXC3T,KERR3T,MAXERT,ENER3T(21),ER3T(21)
      COMMON/SLIM/ISTART,NOTHIN,ITHIN1,ITHIN2,ITHIN3,MTEND
CAK   COMMON/IWATCH/MONITR,MAKEPLUS,MYUNRES
      COMMON/IWATCHs/MONITR,MAKEPLUS,MYUNRES
      COMMON/RESOLV/EULOW,EUHIGH,UREVAL,UREACT,UNRES1,UNRES2,IL,IH
      COMMON/FILLER/N2LEFT,TOOHI,ITHRES,LOAD1,LOAD2
CAK   COMMON/FIELDC/FIELD6(11,6)
      COMMON/FIELDCs/FIELD6(11,6)
      COMMON/HIGHESTE/EHIGHEST,EMILLION
      INCLUDE 'sigma1.h'
      DATA IMUSED/0/
C-----THERMAL ENERGY - ALWAYS KEEP
      DATA ETHERMAL/2.53D-02/
C-----------------------------------------------------------------------
C
C     IF MAXIMUM ALLOWABLE ERROR IS NOT POSITIVE NO THINNING CAN BE
C     PERFORMED ON THE TABLE. DO NOT THIN REACTIONS (MT) THAT WERE
C     NOT BROADENED.
C
C-----------------------------------------------------------------------
      IF(NOTHIN.LE.0.0.AND.TOOHI.LE.0) GO TO 10
C-----NO THINNING. SET THINNED POINT INDEX TO LAST BROADENED POINT.
      ITHIN1=HOT3
      GO TO 190
C-----------------------------------------------------------------------
C
C     THINNING WILL BE PERFORMED.
C
C-----------------------------------------------------------------------
C-----IF NO POINTS FOR THINNING SKIP TO END OF THINNING SECTION.
   10 IF(ITHIN2.GT.HOT3) GO TO 180
C-----------------------------------------------------------------------
C
C     INITIALIZE SIGN OF DERIVATIVE AT BEGINNING OF SECTION.
C
C-----------------------------------------------------------------------
      IF(ISTART.le.0) go to 20
      ISTART=0
      DSIGNT=1.0
      IF(XCHOT(1).GT.XCHOT(2)) DSIGNT=-1.0
C-----------------------------------------------------------------------
C
C     SET UP LOOP OVER BROADENED POINTS.
C
C-----------------------------------------------------------------------
   20 DO 170 M=ITHIN2,HOT3
      MM1=M-1
      DXC=XCHOT(M)-XCHOT(MM1)
C-----------------------------------------------------------------------
C
C     KEEP ALL MAXIMA AND MINIMA.
C
C-----------------------------------------------------------------------
C-----PRECEDING POINT WAS MAXIMUM OR MINIMUM IF SIGN OF CHANGE IN
C-----CROSS SECTION HAS REVERSED.
      IF(DXC*DSIGNT.ge.0.0D+0) go to 30
      DSIGNT=-DSIGNT
C-----IF ENERGIES OF TWO POINTS ARE THE SAME TREAT AS DISCONTINUITY.
      IF(EHOT(M).LE.EHOT(MM1)) GO TO 50
C-----ENERGY NOT THE SAME. SAVE MAXIMUM OR MINIMUM IF IT HAS NOT
C-----ALREADY BEEN SAVED.
      IF(MM1.ge.ITHIN3) go to 130
C-----MAXIMUM OR MINIMUM ALREADY SAVED.
C-----------------------------------------------------------------------
C
C     KEEP PRECEDING AND PRESENT POINTS IF SIGN OF CROSS SECTION HAS
C     CHANGED.
C
C-----------------------------------------------------------------------
   30 IF(XCHOT(M)*XCHOT(MM1).ge.0.0D+0) go to 40
C-----IF PRECEDING POINT HAS NOT BEEN KEPT KEEP PRECEDING AND CURRENT
C-----POINTS, OTHERWISE KEEP ONLY CURRENT POINT.
      IF(M.le.ITHIN3) go to 150
      go to 140
C-----------------------------------------------------------------------
C
C     KEEP DISCONTINUITY.
C
C-----------------------------------------------------------------------
C-----IS ENERGY OF TWO POINTS THE SAME.
   40 IF(EHOT(M).GT.EHOT(MM1)) GO TO 80
C-----YES. CHECK FOR SAME CROSS SECTION.
   50 DXC=DABS(XCHOT(M)-XCHOT(MM1))
C-----YES. CHECK FOR BEGINNING OF THINNING INTERVAL.
      IF(M.NE.ITHIN3) GO TO 60
C-----BEGINNING OF INTERVAL (M-1 ALREADY SAVED). IF CROSS SECTIONS ARE,
C-----(1) SAME -SKIP M AND UPDATE INDEX DEFINING M+1 BEGINNING INTERVAL
C-----(2) DIFFERENT - KEEP POINT M.
      IF(DXC.le.0.0D+0) go to 160
      go to 70
C-----NOT BEGINNING OF INTERVAL (M-1 NOT SAVED). IF CROSS SECTIONS ARE,
C-----(1) SAME - KEEP M, SKIP M-1.
C-----(2) DIFFERENT - KEEP POINTS M-1 AND M.
   60 IF(DXC.gt.0.0D+0) go to 140
   70 IF(M.le.ITHIN3) go to 150
      go to 140
C-----------------------------------------------------------------------
C
C     2013/9/18 - KEEP ALL POINTS ABOVE HIGH ENERGY CUTOFF.
C
C-----------------------------------------------------------------------
   80 IF(EHOT(M).LT.EHIGHEST) GO TO 90
      IF(M.le.ITHIN3) go to 150
      go to 140
C-----------------------------------------------------------------------
C
C     KEEP THERMAL POINT AND PRECEDING POINT, IF NOT ALREADY KEPT.
C
C-----------------------------------------------------------------------
   90 IF(EHOT(M).NE.ETHERMAL) GO TO 100
      IF(M.le.ITHIN3) go to 150
      go to 140
C-----------------------------------------------------------------------
C
C     KEEP ALL POINTS IN THE UNRESOLVED RESONANCE REGION.
C
C-----------------------------------------------------------------------
  100 IF(UREACT.LE.0) GO TO 110
      IF(EHOT(M).LT.EULOW.OR.EHOT(M).GT.EUHIGH) GO TO 110
      IF(M.le.ITHIN3) go to 150
      go to 140
C-----------------------------------------------------------------------
C
C     DEFINE SLOPE OF STRAIGHT LINE THAT WILL PASS WITHIN THE ALLOWABLE
C     ERROR OF EACH POINT. KEEP ELIMINATING POINTS UNTIL ONE OR MORE
C     POINTS CANNOT BE APPROXIMATED TO WITHIN THE ALLOWABLE ERROR. AT
C     THAT POINT KEEP THE LAST PRECEDING POINT (I.E., KEEP THE LAST
C     POINT THAT PASASED THE TEST).
C
C-----------------------------------------------------------------------
C-----DEFINE ENERGY INTERVAL BETWEEN CURRENT POINT AND LAST CONVERGED
C-----POINT.
  110 DE=EHOT(M)-EHOT(ITHIN1)
      SLOPE=(XCHOT(M)-XCHOT(ITHIN1))/DE
C-----INITIALIZE MAXIMUM AND MINIMUM ALLOWABLE SLOPE AT FIRST POINT OF
C-----INTERVAL.
      IF(M.NE.ITHIN3) GO TO 120
      IF(KERR3T.NE.0) CALL ERROKT(EHOT(M))
      DSLOPE=ERXC3T*XCHOT(M)/DE
      SLPMAX=SLOPE+DSLOPE
      SLPMIN=SLOPE-DSLOPE
      GO TO 170
C-----AFTER FIRST POINT OF INTERVAL SEE IF SLOPE TO CURRENT POINT PASSES
C-----WITHIN THE ALLOWABLE ERROR OF ALL PRECEDING POINTS IN CURRENT
C-----INTERVAL.
  120 IF(SLOPE.GT.SLPMAX.OR.SLOPE.LT.SLPMIN) GO TO 130
C-----CAN ELIMINATE CURRENT POINT. UPDATE SLOPE LIMITS.
      IF(KERR3T.NE.0) CALL ERROKT(EHOT(M))
      DSLOPE=ERXC3T*XCHOT(M)/DE
      SLP1=SLOPE+DSLOPE
      IF(SLP1.LT.SLPMAX) SLPMAX=SLP1
      SLP2=SLOPE-DSLOPE
      IF(SLP2.GT.SLPMIN) SLPMIN=SLP2
      GO TO 170
C-----NEED TO KEEP LAST PRECEDING POINT (LAST ONE TO PASS TEST).
  130 ITHIN1=ITHIN1+1
      EHOT(ITHIN1)=EHOT(MM1)
      XCHOT(ITHIN1)=XCHOT(MM1)
C-----RE-DEFINE INDEX TO BEGINNING OF NEXT INTERVAL.
      ITHIN3=M
      GO TO 80
C-----NEED TO KEEP LAST PRECEDING AND CURRENT POINTS.
  140 ITHIN1=ITHIN1+1
      EHOT(ITHIN1)=EHOT(MM1)
      XCHOT(ITHIN1)=XCHOT(MM1)
C-----NEED TO KEEP CURRENT POINT.
  150 ITHIN1=ITHIN1+1
      EHOT(ITHIN1)=EHOT(M)
      XCHOT(ITHIN1)=XCHOT(M)
C-----RE-DEFINE INDEX TO BEGINNING OF NEXT INTERVAL.
  160 ITHIN3=M+1
C-----END OF THINNING LOOP.
  170 CONTINUE
C-----IF LAST POINT OF ARRAY WAS NOT SAVED SAVE IT NOW.
      IF(ITHIN3.GT.HOT3) GO TO 180
      ITHIN1=ITHIN1+1
      EHOT(ITHIN1)=EHOT(HOT3)
      XCHOT(ITHIN1)=XCHOT(HOT3)
C-----RE-DEFINE NUMBER OF BROADENED POINTS IN CORE.
  180 HOT3=ITHIN1
      HOT3M1=HOT3-1
  190 CONTINUE
C-----HAS LAST POINT OF CROSS SECTION BEEN DOPPLER BROADENED.
      IF(MTEND.GT.0) GO TO 230
C-----------------------------------------------------------------------
C
C     MORE POINTS REMAIN TO DOPPLER BROADEN. IF MORE THAN ONE PAGE
C     OF THINNED POINTS COPY ONE PAGE TO SCRATCH AND MOVE OTHER POINTS
C     ONE PAGE FORWARD IN CORE. OTHERWISE NOTHING MORE TO DO.
C
C-----------------------------------------------------------------------
C-----IF THERE IS A FULL PAGE OF THINNED POINTS UNLOAD THEM TO SCRATCH.
      IF(ITHIN1.LE.NPAGE) GO TO 220
C-----POSITION SCRATCH BEFORE FIRST WRITE.
      IF(N2SCR.EQ.0.AND.IMUSED.GT.0) REWIND SCR
      IMUSED=1
C-----COPY PAGE TO SCRATCH AND INCREMENT SCRATCH POINT COUNT.
      N2SCR=N2SCR+NPAGE
      WRITE(SCR) EHOT1,XCHOT1
      IF(MONITR.LE.0) GO TO 200
      CALL OUT9(EHOT1(1)    ,FIELD6(1,1))
      CALL OUT9(EHOT1(NPAGE),FIELD6(1,2))
      WRITE(OUTP,290) N2SCR,((FIELD6(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,290) N2SCR,((FIELD6(M,KK),M=1,11),KK=1,2)
C-----MOVE REMAINING POINTS FORWARD IN CORE.
  200 K=0
      DO 210 J=NP1P1,ITHIN1
      K=K+1
      EHOT(K)=EHOT(J)
  210 XCHOT(K)=XCHOT(J)
C-----DEFINE INDEX TO LAST THINNED POINT.
      ITHIN1=K
      HOT3=K
      HOT3M1=HOT3-1
C-----DEFINE THINNING INDICES FOR NEXT TIME THAT THIS ROUTINE WILL BE
C-----CALLED.
  220 ITHIN2=HOT3+1
      ITHIN3=HOT3+1
      RETURN
C-----------------------------------------------------------------------
C
C     END OF REACTION. IF ALL DATA IS CORE RESIDENT LEAVE IT IN CORE.
C     OTHERWISE COPY ALL TO SCRATCH AND INDICATE NONE REMAINING IN
C     CORE (AT THIS POINT THERE MAY BE UP TO THREE PAGES OF CORE
C     RESIDENT DATA). IF SCRATCH FILE IS USED POSITION IT TO BE READ.
C
C     FIRST PAGE OUTPUT.
C
C-----------------------------------------------------------------------
  230 CONTINUE
C 220 IF(N2SCR.LE.0) GO TO 270
      IF(N2SCR.LE.0) GO TO 280
      IF(MONITR.LE.0) GO TO 240
      ITOP=ITHIN1
      IF(ITOP.GT.NPAGE) ITOP=NPAGE
      N2SCRP=N2SCR+ITOP
      CALL OUT9(EHOT(1)   ,FIELD6(1,1))
      CALL OUT9(EHOT(ITOP),FIELD6(1,2))
      WRITE(OUTP,290) N2SCRP,((FIELD6(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,290) N2SCRP,((FIELD6(M,KK),M=1,11),KK=1,2)
  240 WRITE(SCR) EHOT1,XCHOT1
C-----------------------------------------------------------------------
C
C     SECOND PAGE OUTPUT.
C
C-----------------------------------------------------------------------
      IF(ITHIN1.LE.NPAGE) GO TO 270
      IF(MONITR.LE.0) GO TO 250
      ITOP=ITHIN1
      IF(ITOP.GT.NPT2) ITOP=NPT2
      N2SCRP=N2SCR+ITOP
      CALL OUT9(EHOT2(1)   ,FIELD6(1,1))
      CALL OUT9(EHOT (ITOP),FIELD6(1,2))
      WRITE(OUTP,290) N2SCRP,((FIELD6(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,290) N2SCRP,((FIELD6(M,KK),M=1,11),KK=1,2)
  250 WRITE(SCR) EHOT2,XCHOT2
C-----------------------------------------------------------------------
C
C     THIRD PAGE OUTPUT.
C
C-----------------------------------------------------------------------
      IF(ITHIN1.LE.NPT2) GO TO 270
      IF(MONITR.LE.0) GO TO 260
      N2SCRP=N2SCR+ITHIN1
      CALL OUT9(EHOT3(1)     ,FIELD6(1,1))
      CALL OUT9(EHOT (ITHIN1),FIELD6(1,2))
      WRITE(OUTP,290) N2SCRP,((FIELD6(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,290) N2SCRP,((FIELD6(M,KK),M=1,11),KK=1,2)
  260 WRITE(SCR) EHOT3,XCHOT3
C-----------------------------------------------------------------------
C
C     ALL POINTS ARE NOW ON SCRATCH.
C
C-----------------------------------------------------------------------
C-----INCREMENT SCRATCH POINT COUNT AND INDICATE NONE REMAIN IN CORE.
  270 N2SCR=N2SCR+ITHIN1
      ITHIN1=0
C-----END FILE AND POSITION SCRATCH TO BE READ.
      END FILE SCR
      REWIND SCR
C-----DEFINE FINAL NUMBER OF POINTS TO OUTPUT.
  280 N2TOT=ITHIN1+N2SCR
      RETURN
  290 FORMAT(33X,I9,11A1,' to',11A1,' eV Finished')
      END
      SUBROUTINE COPOUT
C=======================================================================
C
C     THIS ROUTINE IS DESIGNED TO COPY THE SECTION OF BROADENED AND/OR
C     THINNED DATA FROM THE CORE OR SCRATCH FILE TO THE RESULT FILE IN
C     THE ENDF/B FORMAT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OTAPE,OUTP,SCR,COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,
     1 HOT3M1,TOOHI
      CHARACTER*1 ZABCD,FIELD6
      CHARACTER*4 MESSAG,FMTHOL,PROHOL
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/UNITS/SCR
      COMMON/UNITSs/SCR
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/C1,Q,L1,L2,N1,N2,MAT,MF,MT
CAK   COMMON/WHATZA/IZA
      COMMON/WHATZAs/IZA
      COMMON/INDEX/COLD1,COLD2,COLD1P,COLD2P,HOT1,HOT2,HOT3,HOT3M1,N2IN,
     1 N2TOT,N2SCR
CAK   COMMON/PAGER/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/PAGERs/NPAGE,NPT2,NPT3,NP1P1,NP2P1
      COMMON/FILLER/N2LEFT,TOOHI,ITHRES,LOAD1,LOAD2
      COMMON/THRESH/ETHRES,EMIN
      COMMON/FLAGS/MINUS3,IMPLUS
      COMMON/MATTOT/MATIN,MATOUT
      COMMON/TEMPO/TEMP3,IVERSE
      COMMON/EXTEND/DTMAX,MESS
CAK   COMMON/HOLFMT/FMTHOL,PROHOL
      COMMON/HOLFMTs/FMTHOL,PROHOL
      COMMON/POINT1/EIN1,EUSE1
CAK   COMMON/FIELDC/FIELD6(11,6)
      COMMON/FIELDCs/FIELD6(11,6)
CAK   COMMON/IWATCH/MONITR,MAKEPLUS,MYUNRES
      COMMON/IWATCHs/MONITR,MAKEPLUS,MYUNRES
      INCLUDE 'sigma1.h'
c----- 2012 - deleted E3UT(3), XC3OUT(3)
      DIMENSION MESSAG(3,3),NBTO(1),INTO(1),ZABCD(10)
C-----DEFINE CROSS SECTION EXTENSION MESSAGE.
      DATA MESSAG/'    ','    ','    ',
     1            'Exte','nsio','n   ',
     2            'Copi','ed  ','    '/
C-----DEFINE LINEAR-LINEAR INTERPOLATION LAW.
      DATA INTO/2/
C-----------------------------------------------------------------------
C
C     DEFINE TABLE SIZE, OUTPUT BEGINNING OF TAB1 RECORD AND DEFINE
C     WHERE DATA IS (I.E. ON SCRATCH OR IN CORE).
C
C-----------------------------------------------------------------------
C-----2013/3/17 - Output 'Copied' IF NOT BRODENED.
      IF(TOOHI.GT.0) MESS = 3
C-----INITIALIZE NEGATIVE CROSS SECTION FLAG OFF.
      MINUS3=0
      IMPLUS=0
C-----OUTPUT TAB1 LEAD CARD (SECTION HEAD CARD ALREADY OUTPUT IN MAIN)
      CALL CARDO(C1,Q,L1,L2,N1,N2TOT)
C-----OUTPUT INTERPOLATION LAW.
      NBTO(1)=N2TOT
      CALL TERPO(NBTO,INTO,1)
C-----DEFINE NUMBER OF DATA POINTS IN CORE (ALL IF CORE RESIDENT, OR
C-----THREE PAGES IF DATA IS ON SCRATCH).
      IOUT=N2TOT
      IF(N2SCR.GT.0) IOUT=NPT3
C-----------------------------------------------------------------------
C
C     SET UP LOOP OVER PAGES OF DATA.
C
C-----------------------------------------------------------------------
      LOOP=1
C-----IF REQUIRED LOAD UP TO THREE THINNED PAGES INTO CORE.
   10 IF(N2SCR.LE.0) GO TO 30
      READ(SCR) EHOT1,XCHOT1
      IF(N2SCR.LE.NPAGE) GO TO 20
      READ(SCR) EHOT2,XCHOT2
      IF(N2SCR.LE.NPT2) GO TO 20
      READ(SCR) EHOT3,XCHOT3
C-----DEFINE NUMBER OF POINTS IN CORE AND DECREMENT COUNT OF POINTS
C-----ON SCRATCH.
   20 IF(N2SCR.LT.IOUT) IOUT=N2SCR
      N2SCR=N2SCR-NPT3
C-----IF THRESHOLD REPLACED DUE TO DOPPLER BROADENING INSURE THAT IF
C-----FIRST TABULATED ENERGY IS ABOVE MINIMUM ENERGY OF INTEREST THE
C-----CROSS SECTION IS ZERO AT THE NEW THRESHOLD.
   30 IF(ITHRES.le.0) go to 40
      ITHRES=0
      IF(EHOT(1).GT.1.1*EMIN) XCHOT(1)=0.0
C-----------------------------------------------------------------------
C
C     SELECT NORMAL OUTPUT OR OUTPUT WITH ADDITIONAL POINT AT LOWEST
C     ENERGY READ AS INPUT.
C
C-----------------------------------------------------------------------
   40 IF(EUSE1.LE.EIN1) GO TO 50
      XCHOT(1) = 0.0
C-----INSURE THIS CAN ONLY BE DONE ONCE.
      EUSE1=-10.0
      EIN1 = 10.0
C-----------------------------------------------------------------------
C
C     NORMAL OUTPUT ROUTE.
C
C-----------------------------------------------------------------------
C-----IF REQUESTED MAKE ALL NEGATIVE CROSS SECTIONS = 0
   50 IF(MAKEPLUS.EQ.1) THEN
      DO KP=1,IOUT
      IF(XCHOT(KP).LT.0.0D+00) XCHOT(KP)=0.0D+00
      ENDDO
      ENDIF
      CALL POINTO(EHOT,XCHOT,IOUT)
C-----------------------------------------------------------------------
C
C     END OF PAGE LOOP.
C
C-----------------------------------------------------------------------
      LOOP=LOOP+NPT3
      IF(LOOP.LE.N2TOT) GO TO 10
C-----------------------------------------------------------------------
C
C     OUTPUT SECTION REPORT AND INCREMENT POINT COUNTS FOR MATERIAL.
C     PRINT WARNING IF OUTPUT CONTAINS ANY NEGATIVE CROSS SECTIONS.
C
C-----------------------------------------------------------------------
      CALL ZAHOL(IZA,ZABCD)
      CALL OUT9(TEMP3,FIELD6(1,1))
      CALL OUT9(Q    ,FIELD6(1,2))
      WRITE(OUTP,60) PROHOL,ZABCD,MATH,MTH,FMTHOL,
     1 ((FIELD6(M,L),M=1,11),L=1,2),N2IN,N2TOT,
     2 (MESSAG(KK,MESS),KK=1,3)
      WRITE(*   ,60) PROHOL,ZABCD,MATH,MTH,FMTHOL,
     1 ((FIELD6(M,L),M=1,11),L=1,2),N2IN,N2TOT,
     2 (MESSAG(KK,MESS),KK=1,3)
      MATIN=MATIN+N2IN
      MATOUT=MATOUT+N2TOT
C-----PRINT WARNING MESSAGE IF THIS SECTION CONTAINS NEGATIVE CROSS
C-----SECTIONS.
      IF(MINUS3.GT.0) WRITE(OUTP,70) MINUS3
      IF(MINUS3.GT.0) WRITE(*   ,70) MINUS3
C-----PRINT WARNING IF CROSS SECTION IS NOT POSITIVE AT ANY ENERGY.
      IF(IMPLUS.LE.0) WRITE(OUTP,80)
      IF(IMPLUS.LE.0) WRITE(*   ,80)
      RETURN
   60 FORMAT(1X,A4,10A1,I5,I4,2X,A2,2X,2(1X,11A1),2I9,1X,2A4,A1)
   70 FORMAT(19X,'WARNING - Above Cross Section Negative at',I9,
     1 ' Energies')
   80 FORMAT(19X,'WARNING - Above Cross Section NOT',
     1 ' Positive at ANY Energy')
      END
CAK   SUBROUTINE READIN
      SUBROUTINE READINs(infile,outfile,Tres)
C=======================================================================
C
C     READ AND CHECK ALL INPUT PARAMETERS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 FIELD6
      CHARACTER*4 MESS1,MESS2
CAK
      CHARACTER*72 infile,outfile
CAK
      CHARACTER*72 NAMEIN,NAMEOUT
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/IOSTATUS/ISTAT1,ISTAT2
CAK   COMMON/MATZA/MODGET,NMATZA,MATMIN(101),MATMAX(101)K
      COMMON/MATZAs/MODGET,NMATZA,MATMIN(101),MATMAX(101)
      COMMON/OKERRT/ERXC3T,KERR3T,MAXERT,ENER3T(21),ER3T(21)
      COMMON/OKERRC/ERXC3C,ERXC30,KERR3C,MAXERC,ENER3C(21),ER3C(21)
      COMMON/HOTS/ALPHA,HOTSY1,HOTSY2,TEMPK,TEMPEF,N2TAPI,N2TAPO
      COMMON/SLIM/ISTART,NOTHIN,ITHIN1,ITHIN2,ITHIN3,MTEND
CAK   COMMON/IWATCH/MONITR,MAKEPLUS,MYUNRES
      COMMON/IWATCHs/MONITR,MAKEPLUS,MYUNRES
CAK   COMMON/PARAMS/XCMIN
      COMMON/PARAMSs/XCMIN
      COMMON/NAMEX/NAMEIN,NAMEOUT
CAK   COMMON/FIELDC/FIELD6(11,6)
      COMMON/FIELDCs/FIELD6(11,6)
      DIMENSION MESS1(2),MESS2(2)
      DATA MESS1/' MAT','  ZA'/
      DATA MESS2/' Off','  On'/
C-----DEFINE DEFAULT OPTION FOR MINIMUM CROSS SECTION.
      DATA XCMIN1/1.0D-10/
C-----------------------------------------------------------------------
C
C     READ AND CHECK FIRST LINE OF INPUT PARAMETERS.
C
C-----------------------------------------------------------------------
CAK   IF(ISTAT1.EQ.1) GO TO 20
c----- 2017 - Changed all floating point to character
CAK   READ(INP,10,END=20,ERR=20) MODGET,MONITR,
CAK  1 ((FIELD6(j,i),j=1,11),i=1,2),MAKEPLUS,MYUNRES
CAK10 FORMAT(2I11,22A1,2I11)
CAK   CALL IN9(TEMPK,FIELD6(1,1))
CAK   CALL IN9(XCMIN,FIELD6(1,2))
c----- 2017 - Changed all floating point to character
CAK   GO TO 30
C-----DEFINE DEFAULT VALUES
CAK20 ISTAT1   = 1
      MODGET   = 0
      MONITR   = 0
C-----09/01/2013 - Changed from 300 to 293.6 K
CAK   TEMPK    = 293.6d+0
      TEMPK    = Tres
      XCMIN    = 0.0
      MAKEPLUS = 1
      MYUNRES  = 0
c-----RETRIEVAL MODE
   30 IF(MODGET.NE.0) MODGET=1
c-----MONITOR MODE
      IF(MONITR.NE.0) MONITR=1
C-----TEMPERATURE
      CALL OUT9(TEMPK,FIELD6(1,1))
      IF(XCMIN.LE.0.0) GO TO 40
C-----MINIMUM CROSS SECTION
      CALL OUT9(XCMIN,FIELD6(1,2))
      WRITE(OUTP,280) MESS1(MODGET+1),MESS2(MONITR+1),
     1 ((FIELD6(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,280) MESS1(MODGET+1),MESS2(MONITR+1),
     1 ((FIELD6(M,KK),M=1,11),KK=1,2)
      GO TO 50
C-----USE DEFAULT OPTION FOR MINIMUM CROSS SECTION.
   40 XCMIN=XCMIN1
      CALL OUT9(XCMIN,FIELD6(1,2))
      WRITE(OUTP,290) MESS1(MODGET+1),MESS2(MONITR+1),
     1 ((FIELD6(M,KK),M=1,11),KK=1,2)
      WRITE(*   ,290) MESS1(MODGET+1),MESS2(MONITR+1),
     1 ((FIELD6(M,KK),M=1,11),KK=1,2)
C-----NEGATIVE CROSS SECTION TREATMENT
   50 IF(MAKEPLUS.NE.0) MAKEPLUS=1
      IF(MAKEPLUS.EQ.0) THEN
      WRITE(OUTP,300)
      WRITE(*   ,300)
      ELSE
      WRITE(OUTP,310)
      WRITE(*   ,310)
      ENDIF
C-----UNRESOLVED TREATMENT
      IF(MYUNRES.NE.0) MYUNRES = 1
      IF(MYUNRES.EQ.0) THEN
      WRITE(OUTP,320)
      WRITE(*   ,320)
      ELSE
      WRITE(OUTP,330)
      WRITE(*   ,330)
      ENDIF
C-----------------------------------------------------------------------
C
C     READ FILENAMES - IF BLANK USE STANDARD FILENAMES
C
C-----------------------------------------------------------------------
C-----INPUT DATA.
CAK   IF(ISTAT1.EQ.1) GO TO 70
CAK   READ(INP,60,END=70,ERR=70) NAMEIN
      NAMEIN=infile
CAK60 FORMAT(A72)
      IF(NAMEIN.EQ.' ') NAMEIN = 'ENDFB.IN'
C-----OUTPUT DATA.
CAK   READ(INP,60,END=80,ERR=80) NAMEOUT
      NAMEOUT=outfile
      IF(NAMEOUT.EQ.' ') NAMEOUT = 'ENDFB.OUT'
      GO TO 90
C-----USE DEFAULT NAMES
CAK70 NAMEIN  = 'ENDFB.IN'
   80 NAMEOUT = 'ENDFB.OUT'
CAK   ISTAT1 = 1
C-----PRINT FINAL FILENAME
   90 WRITE(OUTP,100) NAMEIN,NAMEOUT
      WRITE(*   ,100) NAMEIN,NAMEOUT
  100 FORMAT(1X,79('-')/
     1 ' ENDF/B Input and Output Data Filenames'/1X,79('-')/
     2 1X,A72/1X,A72)
C-----------------------------------------------------------------------
C
C     OPEN ENDF/B DATA FILES
C
C-----------------------------------------------------------------------
CAK   CALL FILIO2
      CALL FILIO2s
C-----------------------------------------------------------------------
C
C     TERMINATE IF ERROR OPENING FILE
C
C-----------------------------------------------------------------------
      IF(ISTAT2.EQ.1) THEN
      WRITE(OUTP,110) NAMEIN
      WRITE(   *,110) NAMEIN
  110 FORMAT(//' ERROR - Opening ENDF/B formatted file'/1X,A72//)
      CALL ENDERROR
      ENDIF
C-----------------------------------------------------------------------
C
C     READ SELECTION RANGES (EITHER MAT OR ZA). IF MAXIMUM IS LESS
C     THAN MINIMUM SET IT EQUAL TO MINIMUM. IF FIRST CARD IS BLANK
C     RETRIEVE ALL DATA.
C
C-----------------------------------------------------------------------
      IF(MODGET.EQ.0) WRITE(OUTP,340)
      IF(MODGET.EQ.1) WRITE(OUTP,350)
      IF(MODGET.EQ.0) WRITE(*   ,340)
      IF(MODGET.EQ.1) WRITE(*   ,350)
CAK   IF(ISTAT1.EQ.1) GO TO 140
CAK   READ(INP,120,END=130,ERR=130) MATMIN(1),MATMAX(1)
CAK120 FORMAT(2I11)
      IF(MATMIN(1).GT.0.OR.MATMAX(1).GT.0) GO TO 150
      GO TO 140
CAK130 ISTAT1    = 1
  140 MATMAX(1) = 9999
      MODGET=0
      NMATZA=2
      WRITE(OUTP,370) MATMIN(1),MATMAX(1)
      WRITE(*   ,370) MATMIN(1),MATMAX(1)
      GO TO 180
c-----Check input and define defaults.
  150 IF(MATMAX(1).LT.MATMIN(1)) MATMAX(1)=MATMIN(1)
      WRITE(OUTP,360) MATMIN(1),MATMAX(1)
      WRITE(*   ,360) MATMIN(1),MATMAX(1)
      DO 160 NMATZA=2,101
CAK   READ(INP,120,END=170,ERR=170) MATMIN(NMATZA),MATMAX(NMATZA)
c-----Check input and define defaults.
      IF(MATMIN(NMATZA).LE.0.AND.MATMAX(NMATZA).LE.0) GO TO 180
      IF(MATMAX(NMATZA).LT.MATMIN(NMATZA)) MATMAX(NMATZA)=MATMIN(NMATZA)
      WRITE(OUTP,360) MATMIN(NMATZA),MATMAX(NMATZA)
      WRITE(*   ,360) MATMIN(NMATZA),MATMAX(NMATZA)
  160 CONTINUE
      WRITE(OUTP,380)
      WRITE(*   ,380)
      CALL ENDERROR
C-----------------------------------------------------------------------
C
C     READ AND LIST FILE 3 ERROR LAW. ENERGIES MUST BE IN ASCENDING
C     ORDER. IT IS O.K. IF ERROR IS ZERO (WHICH INDICATES THAT
C     THINNING SHOULD NOT BE PERFORMED). ERROR LAW IS TERMINATED BY
C     BLANK CARD. IF FIRST CARD IS BLANK, TERMINATE ERROR LAW AND
C     DEFINE ALLOWABLE ERROR TO BE 0.0 (I.E., NO THINNING).
C
C-----------------------------------------------------------------------
CAK170 ISTAT1 = 1
  180 NMATZA=NMATZA-1
CAK   IF(ISTAT1.EQ.1) GO TO 200
c----- 2017 - Changed all floating point to character
CAK   READ(INP,190,END=200,ERR=200) ((FIELD6(j,i),j=1,11),i=1,2)
CAK190 FORMAT(22A1)
CAK   CALL IN9(ENER3T(1),FIELD6(1,1))
CAK   CALL IN9(ER3T  (1),FIELD6(1,2))
c----- 2017 - Changed all floating point to character
      GO TO 210
CAK200 ISTAT1    = 1
      ENER3T(1) = 0.0
      ER3T(1)   = 0.0
  210 IF(ENER3T(1).LE.0.0) ENER3T(1)=0.0
      IF(ER3T(1).LE.0.0) ER3T(1)=0.0
      IF(ENER3T(1).GT.0.0.OR.ER3T(1).GT.0.0) GO TO 220
C-----USE DEFAULT VALUES.
      MAXERT=2
      ENER3T(1)=0.0
      ER3T(1)=0.0001  ! 1/16/09 - CHANGED DEFAULT TO 0.01%
      PERCNT=100.0*ER3T(1)
      ERRMAX=0.0
      CALL OUT9(ENER3T(1),FIELD6(1,1))
      CALL OUT9(ER3T(1)  ,FIELD6(1,2))
      WRITE(OUTP,420) ((FIELD6(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,420) ((FIELD6(M,I),M=1,11),I=1,2),PERCNT
      GO TO 250
C-----USE VALUES AS READ.
  220 IF(ER3T(1).LE.0.0) ER3T(1)=0.001
      PERCNT=100.0*ER3T(1)
      ERRMAX=ER3T(1)
      CALL OUT9(ENER3T(1),FIELD6(1,1))
      CALL OUT9(ER3T(1)  ,FIELD6(1,2))
      WRITE(OUTP,410) ((FIELD6(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,410) ((FIELD6(M,I),M=1,11),I=1,2),PERCNT
      DO 230 MAXERT=2,21
c----- 2017 - Changed all floating point to character
CAK   READ(INP,190,END=200,ERR=200) ((FIELD6(j,i),j=1,11),i=1,2)
CAK   CALL IN9(ENER3T(MAXERT),FIELD6(1,1))
CAK   CALL IN9(ER3T  (MAXERT),FIELD6(1,2))
c----- 2017 - Changed all floating point to character
      IF(ENER3T(MAXERT).LE.0.0.AND.ER3T(MAXERT).LE.0.0) GO TO 250
      IF(ER3T(MAXERT).LE.0.0) ER3T(MAXERT)=0.001
      PERCNT=100.0*ER3T(MAXERT)
      IF(ER3T(MAXERT).GT.ERRMAX) ERRMAX=ER3T(MAXERT)
      CALL OUT9(ENER3T(MAXERT),FIELD6(1,1))
      CALL OUT9(ER3T(MAXERT)  ,FIELD6(1,2))
      WRITE(OUTP,390) ((FIELD6(M,I),M=1,11),I=1,2),PERCNT
      WRITE(*   ,390) ((FIELD6(M,I),M=1,11),I=1,2),PERCNT
      IF(ENER3T(MAXERT).LT.ENER3T(MAXERT-1)) GO TO 270
  230 CONTINUE
      WRITE(OUTP,400)
      WRITE(*   ,400)
  240 CALL ENDERROR
  250 MAXERT=MAXERT-1
      KERR3T=0
      IF(MAXERT.GT.1) KERR3T=1
C-----SET C = T LAW
      MAXERC=MAXERT
      KERR3C=KERR3T
      DO 260 I=1,MAXERT
      ENER3C(I)=ENER3T(I)
  260 ER3C  (I)=ER3T  (I)
C-----INITIALIZE TO FIRST VALUES
      ERXC3T=ER3T(1)
      ERXC3C=ER3C(1)
      ERXC30=0.1d0*ERXC3C
C-----------------------------------------------------------------------
C
C     DATA WILL ALWAYS BE THINNED
C
C-----------------------------------------------------------------------
      NOTHIN=0
C***** THINNING ON/OFF
C
C     TO TURN OFF THINNING ACTIVATE THE FOLLOWING
C
C     NOTHIN=1
C***** THINNING ON/OFF
      RETURN
  270 WRITE(OUTP,430)
      WRITE(*   ,430)
      GO TO 240
  280 FORMAT(' Retrieval Criteria-----------',7X,A4/
     1 ' Monitor Mode-----------------',7X,A4/
     2 ' Temperature------------------',11A1,' KELVIN'/
     3 ' Mimimum Cross Section--------',11A1)
  290 FORMAT(' Retrieval Criteria-----------',7X,A4/
     1 ' Monitor Mode-----------------',7X,A4/
     2 ' Temperature------------------',11A1,' Kelvin'/
     3 ' Minimum Cross Section--------',11A1,' (Default Option)')
  300 FORMAT(' Negative Cross Section-------',
     1 ' No Change (Allow Negative Outout)')
  310 FORMAT(' Negative Cross Section-------',
     1 ' Make = 0 (No Negative Outout)')
  320 FORMAT(' Unresolved Resonances Region-',' Copy (No Broadening)')
  330 FORMAT(' Unresolved Resonances Region-',' Ignore (Broaden)')
  340 FORMAT(1X,79('-')/' Requested MAT Ranges'/1X,79('-')/
     1 5X,'Mimimum',4X,'Maximum'/1X,79('-'))
  350 FORMAT(1X,79('-')/' Requested ZA Ranges'/1X,79('-')/
     1 5X,'Minimum',4X,'Maximum'/1X,79('-'))
  360 FORMAT(1X,2I11)
  370 FORMAT(1X,2I11,' (Default Option)')
  380 FORMAT(1X,79('-')/' Over 100 Ranges----Execution Terminated')
  390 FORMAT(1X,11A1,1X,11A1,F11.4)
  400 FORMAT(1X,79('-')/' Over 20 Ranges----Execution Terminated')
  410 FORMAT(1X,79('-')/' Allowable Uncertainty'/1X,79('-')/
     1 6X,'Energy',1X,'Uncertainty',3X,'per-cent'/1X,79('-')/
     2 1X,11A1,1X,11A1,F11.4)
  420 FORMAT(1X,79('-')/' Allowable Uncertainty'/1X,79('-')/
     1 6X,'Energy',1X,'Unceryainty',3X,'per-cent'/1X,79('-')/
     2 1X,11A1,1X,11A1,F11.4,' (Default Option)')
  430 FORMAT(' Energies MUST be in Ascending Order----',
     1 'Execution Terminated')
      END
      REAL*8 FUNCTION ESPACE(ERROR)
C=======================================================================
C
C     DEFINE ENERGY SPACING TO APPROXIMATE 1/V TO WITHIN ERROR
C
C     THE BELOW TABLE HANDLES THE FRACTIONAL ERROR RANGE
C     1.0D-02 TO 1.0D-06 = 1.0 TO 0.0001 PER-CENT
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION ERRTAB(5),DETAB(5)
      DATA ERRTAB/
     1 1.0D-02,
     2 1.0D-03,
     3 1.0D-04,
     4 1.0D-05,
     5 1.0D-06/
C-----04/23/08 - HALVED SPACING AGAIN
      DATA DETAB/
     1 1.04888250D+00,
     2 1.01410820D+00,
     3 1.00423800D+00,
     4 1.00134370D+00,
     5 1.00041050D+00/
C-----12/22/06 - HALVED SPACING AGAIN
c     DATA DETAB/
c    1 1.09776500D+00,
c    2 1.02821650D+00,
c    3 1.00847700D+00,
c    4 1.00268750D+00,
c    5 1.00082100D+00/
C-----01/16/04 - HALVED SPACING
c     DATA DETAB/
c    1 1.19553000D+00,
c    2 1.05643300D+00,
c    3 1.01695400D+00,
c    4 1.00537500D+00,
c    5 1.00164200D+00/
C-----12/20/06 - REDUCED SPACING BY FACTOR OF 10 (INSERTED 0 AFTER .)
c     DATA DETAB/
c    1 1.019553000D+00,
c    2 1.005643300D+00,
c    3 1.001695400D+00,
c    4 1.000537500D+00,
c    5 1.000164200D+00/
C-----12/20/06 - REDUCED SPACING BY FACTOR OF 100 (INSERTED 00 AFTER .)
C     DATA DETAB/
C    1 1.0019553000D+00,
C    2 1.0005643300D+00,
C    3 1.0001695400D+00,
C    4 1.0000537500D+00,
C    5 1.0000164200D+00/
      IF(ERROR.LT.ERRTAB(1)) GO TO 10
      ESPACE=DETAB(1)
      RETURN
   10 DO 20 I=1,5
      IF(ERROR.ge.ERRTAB(I)) go to 30
   20 CONTINUE
      I=5
   30 ESPACE=DETAB(I)
      RETURN
      END
CAK   SUBROUTINE NXTMAT
      SUBROUTINE NXTMATs
C=======================================================================
C
C     FIND NEXT REQUESTED MATERIAL BASED EITHER ON ZA OR MAT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*4 FMTHOL,PROHOL
CAK   COMMON/MATZA/MODGET,NMATZA,MATMIN(101),MATMAX(101)
      COMMON/MATZAs/MODGET,NMATZA,MATMIN(101),MATMAX(101)
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
C     COMMON/WHATZA/IZAK
      COMMON/WHATZAs/IZA
CAK   COMMON/HOLFMT/FMTHOL,PROHOL
      COMMON/HOLFMTs/FMTHOL,PROHOL
      COMMON/TEMPO/TEMP3,IVERSE
      DIMENSION IZAMIN(101),IZAMAX(101)
      EQUIVALENCE (MATMIN(1),IZAMIN(1)),(MATMAX(1),IZAMAX(1))
C-----READ NEXT CARD AND CHECK FOR END OF ENDF/B TAPE.
   10 CALL CONTI
      IF(MTH.gt.0) go to 20
      IF(MATH.lt.0) go to 60
      go to 10
C-----DEFINE FIXED POINT ZA.
   20 IZA=C1H
C-----COMPARE MAT OR ZA TO SELECTION CRITERIA.
      LOW=0
      DO 50 IMATZA=1,NMATZA
      IF(MODGET.NE.0) GO TO 30
      IF(MATH.lt.MATMIN(IMATZA)) go to 40
      IF(MATH.eq.MATMIN(IMATZA)) go to 70
      IF(MATH.le.MATMAX(IMATZA)) go to 70
      go to 50
   30 IF(IZA.lt.IZAMIN(IMATZA)) go to 40
      IF(IZA.eq.IZAMIN(IMATZA)) go to 70
      IF(IZA.le.IZAMAX(IMATZA)) go to 70
   40 LOW=1
   50 CONTINUE
C-----THIS MATERIAL HAS NOT BEEN REQUESTED. IF BEYOND RANGE OF ALL
C-----REQUESTS RUN IF COMPLETED. IF NOT SKIP TO NEXT MATERIAL.
      IF(LOW.LE.0) GO TO 60
C-----SKIP TO MATERIAL END (MEND) CARD.
      CALL SKIPM
      GO TO 10
C-----END OF RUN. RETURN NEGATIVE MATH AS INDICATOR.
   60 MATH=-1
      MFH=0
      MTH=0
C-----THIS MATERIAL REQUESTED. INITIALIZE OUTPUT SEQUENCE NUMBER,
C-----ENDF/B FORMAT VERSION TO BLANK, ENDF/B FORMAT TO VERSION V AND
C-----INITIAL TEMPERATURE TO ZERO.
   70 NOSEQ=1
      FMTHOL='VI  '
      IVERSE=6
      TEMP3=0.0
      RETURN
      END
      SUBROUTINE ERROKT(E)
C=======================================================================
C
C     DEFINE ALLOWABLE ERROR FOR FILE 3 LINEARIZED CROSS SECTIONS.
C     THE ERROR LAW CAN BE ENERGY INDEPENDENT (CONSTANT) OR ENERGY
C     DEPENDENT (GIVEN BY A LINEARLY INTERPOLABLE TABLE IN ENERGY
C     VS. ERROR).
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/OKERRT/ERXC3T,KERR3T,MAXERT,ENER3T(21),ER3T(21)
C-----INITIALIZE INDEX TO INTERPOLATION TABLE.
      DATA MINER3/2/
C-----ENERGY DEPENDENT. WITHIN ENERGY RANGE OF ERROR LAW USE LINEAR
C-----INTERPOLATION. OUTSIDE RANGE EXTEND ERROR AS CONSTANT FROM
C-----CLOSEST END OF TABLE.
      IF(E.le.ENER3T(1)) go to 80
      DO 10 NOWER3=MINER3,MAXERT
      IF(E.lt.ENER3T(NOWER3)) go to 20
      IF(E.eq.ENER3T(NOWER3)) go to 70
   10 CONTINUE
C-----EXTEND AS CONSTANT TO HIGHER ENERGIES.
      GO TO 90
   20 NM1=NOWER3-1
      IF(E.eq.ENER3T(NM1)) go to 60
      IF(E.gt.ENER3T(NM1)) go to 50
      DO 30 NOWER3=2,MAXERT
      IF(E.lt.ENER3T(NOWER3)) go to 40
      IF(E.eq.ENER3T(NOWER3)) go to 70
   30 CONTINUE
      GO TO 90
C-----INTERPOLATE BETWEEN ENERGIES.
   40 NM1=NOWER3-1
   50 MINER3=NOWER3
      ERXC3T=((ENER3T(NOWER3)-E)*ER3T(NM1)+(E-ENER3T(NM1))*ER3T(NOWER3))
     1 /(ENER3T(NOWER3)-ENER3T(NM1))
      RETURN
C-----EXACT ENERGY MATCH.
   60 MINER3=NM1
      IF(MINER3.LE.1) MINER3=2
      ERXC3T=ER3T(NM1)
      RETURN
C-----EXACT ENERGY MATCH.
   70 MINER3=NOWER3
      ERXC3T=ER3T(NOWER3)
      RETURN
C-----EXTEND AS CONSTANT TO LOWER ENERGIES.
   80 MINER3=2
      ERXC3T=ER3T(1)
      RETURN
C-----EXTEND AS CONSTANT TO HIGHER ENERGIES.
   90 MINER3=MAXERT
      ERXC3T=ER3T(MAXERT)
      RETURN
      END
      SUBROUTINE ERROKC(E)
C=======================================================================
C
C     DEFINE ALLOWABLE ERROR FOR FILE 3 RECONSTRUCTED CROSS SECTIONS.
C     THE ERROR LAW CAN BE ENERGY INDEPENDENT (CONSTANT) OR ENERGY
C     DEPENDENT (GIVEN BY A LINEARLY INTERPOLABLE TABLE IN ENERGY
C     VS. ERROR).
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/OKERRC/ERXC3C,ERXC30,KERR3C,MAXERC,ENER3C(21),ER3C(21)
C-----INITIALIZE INDEX TO INTERPOLATION TABLE.
      DATA MINERD/2/
C-----ENERGY DEPENDENT. WITHIN ENERGY RANGE OF ERROR LAW USE LINEAR
C-----INTERPOLATION. OUTSIDE RANGE EXTEND ERROR AS CONSTANT FROM
C-----CLOSEST END OF TABLE.
      IF(E.le.ENER3C(1)) go to 80
      DO 10 NOWERD=MINERD,MAXERC
      IF(E.lt.ENER3C(NOWERD)) go to 20
      IF(E.eq.ENER3C(NOWERD)) go to 70
   10 CONTINUE
C-----EXTEND TO HIGHER ENERGIES AS CONSTANT.
      GO TO 90
   20 NM1=NOWERD-1
      IF(E.eq.ENER3C(NM1)) go to 60
      IF(E.gt.ENER3C(NM1)) go to 50
      DO 30 NOWERD=2,MAXERC
      IF(E.lt.ENER3C(NOWERD)) go to 40
      IF(E.eq.ENER3C(NOWERD)) go to 70
   30 CONTINUE
      GO TO 90
C-----INTERPOLATE BETWEEN ENERGIES.
   40 NM1=NOWERD-1
   50 MINERD=NOWERD
      ERXC3C=((ENER3C(NOWERD)-E)*ER3C(NM1)+(E-ENER3C(NM1))*ER3C(NOWERD))
     1 /(ENER3C(NOWERD)-ENER3C(NM1))
      ERXC30=0.1d0*ERXC3C
      RETURN
C-----EXACT ENERGY MATCH.
   60 MINERD=NM1
      IF(MINERD.LE.1) MINERD=2
      ERXC3C=ER3C(NM1)
      ERXC30=0.1d0*ERXC3C
      RETURN
C-----EXACT ENERGY MATCH.
   70 MINERD=NOWERD
      ERXC3C=ER3C(NOWERD)
      ERXC30=0.1d0*ERXC3C
      RETURN
C-----EXTEND TO LOWER ENERGIES AS CONSTANT.
   80 MINERD=2
      ERXC3C=ER3C(1)
      ERXC30=0.1d0*ERXC3C
      RETURN
C-----EXTEND TO HIGHER ENERGIES AS CONSTANT.
   90 MINERD=MAXERC
      ERXC3C=ER3C(MAXERC)
      ERXC30=0.1d0*ERXC3C
      RETURN
      END
CAK   SUBROUTINE FILEIO
      SUBROUTINE FILEIOs
C=======================================================================
C
C     DEFINE ALL I/O UNITS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,SCR
      CHARACTER*72 NAMEIN,NAMEOUT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/IOSTATUS/ISTAT1,ISTAT2
CAK   COMMON/UNITS/SCR
      COMMON/UNITSs/SCR
      COMMON/NAMEX/NAMEIN,NAMEOUT
C-----DEFINE ALL I/O UNIT NUMBERS.
      INP=2
      OUTP=3
      ITAPE=10
      OTAPE=11
      SCR=12
C-----DEFINE ALL FILE NAMES.
      OPEN(OUTP,FILE='SIGMA1.LST',STATUS='UNKNOWN')
      CALL SCRATCH1(SCR,'SIGMA1.001  ')
      OPEN(INP,FILE='SIGMA1.INP',STATUS='OLD',ERR=10)
      ISTAT1 = 0
      RETURN
   10 ISTAT1 = 1
      RETURN
CAK   ENTRY FILIO2
      ENTRY FILIO2s
C=======================================================================
C
C     DEFINE ENDF/B DATA I/O UNITS AND OPTIONAL DEFINE FILE NAMES.
C
C=======================================================================
      OPEN(OTAPE,FILE=NAMEOUT,STATUS='UNKNOWN')
      OPEN(ITAPE,FILE=NAMEIN,STATUS='OLD',ERR=20)
      ISTAT2 = 0
      RETURN
   20 ISTAT2 = 1
      RETURN
      END
      REAL*8 FUNCTION TERFC(x)
C=======================================================================
C
c     compute the complementary error function.
c     from the slatec library fnlib ERFC.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION erfcs(21),erfccs(59),erc2cs(49)
      data erfcs/
     1  -.4904612123 4691808039 9845440333 76 d-1      ,
     2  -.1422612051 0371364237 8247418996 31 d+0      ,
     3  +.1003558218 7599795575 7546767129 33 d-1      ,
     4  -.5768764699 7674847650 8270255091 67 d-3      ,
     5  +.2741993125 2196061034 4221607914 71 d-4      ,
     6  -.1104317550 7344507604 1353812959 05 d-5      ,
     7  +.3848875542 0345036949 9613114981 74 d-7      ,
     8  -.1180858253 3875466969 6317518015 81 d-8      ,
     9  +.3233421582 6050909646 4029309533 54 d-10     ,
     A  -.7991015947 0045487581 6073747085 95 d-12     ,
     1  +.1799072511 3961455611 9672454866 34 d-13     ,
     2  -.3718635487 8186926382 3168282094 93 d-15     ,
     3  +.7103599003 7142529711 6899083946 66 d-17     ,
     4  -.1261245511 9155225832 4954248533 33 d-18     ,
     5  +.2091640694 1769294369 1705002666 66 d-20     ,
     6  -.3253973102 9314072982 3641600000 00 d-22     ,
     7  +.4766867209 7976748332 3733333333 33 d-24     ,
     8  -.6598012078 2851343155 1999999999 99 d-26     ,
     9  +.8655011469 9637626197 3333333333 33 d-28     ,
     B  -.1078892517 7498064213 3333333333 33 d-29     ,
     1  +.1281188399 3017002666 6666666666 66 d-31     /
c
      data erc2cs/
     1  -.6960134660 2309501127 3915082619 7 d-1      ,
     2  -.4110133936 2620893489 8221208466 6 d-1      ,
     3  +.3914495866 6896268815 6114370524 4 d-2      ,
     4  -.4906395650 5489791612 8093545077 4 d-3      ,
     5  +.7157479001 3770363807 6089414182 5 d-4      ,
     6  -.1153071634 1312328338 0823284791 2 d-4      ,
     7  +.1994670590 2019976350 5231486770 9 d-5      ,
     8  -.3642666471 5992228739 3611843071 1 d-6      ,
     9  +.6944372610 0050125899 3127721463 3 d-7      ,
     A  -.1371220902 1043660195 3460514121 0 d-7      ,
     1  +.2788389661 0071371319 6386034808 7 d-8      ,
     2  -.5814164724 3311615518 6479105031 6 d-9      ,
     3  +.1238920491 7527531811 8016881795 0 d-9      ,
     4  -.2690639145 3067434323 9042493788 9 d-10     ,
     5  +.5942614350 8479109824 4470968384 0 d-11     ,
     6  -.1332386735 7581195792 8775442057 0 d-11     ,
     7  +.3028046806 1771320171 7369724330 4 d-12     ,
     8  -.6966648814 9410325887 9586758895 4 d-13     ,
     9  +.1620854541 0539229698 1289322762 8 d-13     ,
     B  -.3809934465 2504919998 7691305772 9 d-14     ,
     1  +.9040487815 9788311493 6897101297 5 d-15     ,
     2  -.2164006195 0896073478 0981204700 3 d-15     ,
     3  +.5222102233 9958549846 0798024417 2 d-16     ,
     4  -.1269729602 3645553363 7241552778 0 d-16     ,
     5  +.3109145504 2761975838 3622741295 1 d-17     ,
     6  -.7663762920 3203855240 0956671481 1 d-18     ,
     7  +.1900819251 3627452025 3692973329 0 d-18     ,
     8  -.4742207279 0690395452 2565599996 5 d-19     ,
     9  +.1189649200 0765283828 8068307845 1 d-19     ,
     C  -.3000035590 3257802568 4527131306 6 d-20     ,
     1  +.7602993453 0432461730 1938527709 8 d-21     ,
     2  -.1935909447 6068728815 6981104913 0 d-21     ,
     3  +.4951399124 7733378810 0004238677 3 d-22     ,
     4  -.1271807481 3363718796 0862198988 8 d-22     ,
     5  +.3280049600 4695130433 1584165205 3 d-23     ,
     6  -.8492320176 8228965689 2479242239 9 d-24     ,
     7  +.2206917892 8075602235 1987998719 9 d-24     ,
     8  -.5755617245 6965284983 1281950719 9 d-25     ,
     9  +.1506191533 6392342503 5414405119 9 d-25     ,
     D  -.3954502959 0187969531 0428569599 9 d-26     ,
     1  +.1041529704 1515009799 8464505173 3 d-26     ,
     2  -.2751487795 2787650794 5017890133 3 d-27     ,
     3  +.7290058205 4975574089 9770368000 0 d-28     ,
     4  -.1936939645 9159478040 7750109866 6 d-28     ,
     5  +.5160357112 0514872983 7005482666 6 d-29     ,
     6  -.1378419322 1930940993 8964480000 0 d-29     ,
     7  +.3691326793 1070690422 5109333333 3 d-30     ,
     8  -.9909389590 6243654206 5322666666 6 d-31     ,
     9  +.2666491705 1953884133 2394666666 6 d-31     /
c
      data erfccs/
     1  +.7151793102 0292477450 3697709496 d-1        ,
     2  -.2653243433 7606715755 8893386681 d-1        ,
     3  +.1711153977 9208558833 2699194606 d-2        ,
     4  -.1637516634 5851788416 3746404749 d-3        ,
     5  +.1987129350 0552036499 5974806758 d-4        ,
     6  -.2843712412 7665550875 0175183152 d-5        ,
     7  +.4606161308 9631303696 9379968464 d-6        ,
     8  -.8227753025 8792084205 7766536366 d-7        ,
     9  +.1592141872 7709011298 9358340826 d-7        ,
     A  -.3295071362 2528432148 6631665072 d-8        ,
     1  +.7223439760 4005554658 1261153890 d-9        ,
     2  -.1664855813 3987295934 4695966886 d-9        ,
     3  +.4010392588 2376648207 7671768814 d-10       ,
     4  -.1004816214 4257311327 2170176283 d-10       ,
     5  +.2608275913 3003338085 9341009439 d-11       ,
     6  -.6991110560 4040248655 7697812476 d-12       ,
     7  +.1929492333 2617070862 4205749803 d-12       ,
     8  -.5470131188 7543310649 0125085271 d-13       ,
     9  +.1589663309 7626974483 9084032762 d-13       ,
     B  -.4726893980 1975548392 0369584290 d-14       ,
     1  +.1435873376 7849847867 2873997840 d-14       ,
     2  -.4449510561 8173583941 7250062829 d-15       ,
     3  +.1404810884 7682334373 7305537466 d-15       ,
     4  -.4513818387 7642108962 5963281623 d-16       ,
     5  +.1474521541 0451330778 7018713262 d-16       ,
     6  -.4892621406 9457761543 6841552532 d-17       ,
     7  +.1647612141 4106467389 5301522827 d-17       ,
     8  -.5626817176 3294080929 9928521323 d-18       ,
     9  +.1947443382 2320785142 9197867821 d-18       ,
     C  -.6826305642 9484207295 6664144723 d-19       ,
     1  +.2421988887 2986492401 8301125438 d-19       ,
     2  -.8693414133 5030704256 3800861857 d-20       ,
     3  +.3155180346 2280855712 2363401262 d-20       ,
     4  -.1157372324 0496087426 1239486742 d-20       ,
     5  +.4288947161 6056539462 3737097442 d-21       ,
     6  -.1605030742 0576168500 5737770964 d-21       ,
     7  +.6063298757 4538026449 5069923027 d-22       ,
     8  -.2311404251 6979584909 8840801367 d-22       ,
     9  +.8888778540 6618855255 4702955697 d-23       ,
     D  -.3447260576 6513765223 0718495566 d-23       ,
     1  +.1347865460 2069650682 7582774181 d-23       ,
     2  -.5311794071 1250217364 5873201807 d-24       ,
     3  +.2109341058 6197831682 8954734537 d-24       ,
     4  -.8438365587 9237891159 8133256738 d-25       ,
     5  +.3399982524 9452089062 7359576337 d-25       ,
     6  -.1379452388 0732420900 2238377110 d-25       ,
     7  +.5634490311 8332526151 3392634811 d-26       ,
     8  -.2316490434 4770654482 3427752700 d-26       ,
     9  +.9584462844 6018101526 3158381226 d-27       ,
     E  -.3990722880 3301097262 4224850193 d-27       ,
     1  +.1672129225 9444773601 7228709669 d-27       ,
     2  -.7045991522 7660138563 8803782587 d-28       ,
     3  +.2979768402 8642063541 2357989444 d-28       ,
     4  -.1262522466 4606192972 2422632994 d-28       ,
     5  +.5395438704 5424879398 5299653154 d-29       ,
     6  -.2380992882 5314591867 5346190062 d-29       ,
     7  +.1099052830 1027615735 9726683750 d-29       ,
     8  -.4867713741 6449657273 2518677435 d-30       ,
     9  +.1525877264 1103575676 3200828211 d-30       /
c
      data sqrtpi / 1.772453850 9055160272 9816748334 115d0 /
c-----------------------------------------------------------------------
c
c     Initialize on first call
c
c-----------------------------------------------------------------------
      DATA IFIRST/1/
      if (IFIRST.ne.0) then
      eta=0.1d0*D1MACH(3)
      nterf =INITDS(erfcs ,21,eta)
      nterfc=INITDS(erfccs,59,eta)
      nterc2=INITDS(erc2cs,49,eta)
      xsml=-DSQRT(-DLOG(sqrtpi*D1MACH(3)))
      txmax=DSQRT(-DLOG(sqrtpi*D1MACH(1)))
      xmax=txmax-0.5d0*DLOG(txmax)/txmax-0.01d0
      sqeps=DSQRT(2.0d0*D1MACH(3))
      IFIRST = 0
      endif
c
      if (x.le.xsml) then
      TERFC=2.0d0
      else if (x.le.xmax) then
      y=abs(x)
      if (y.le.1.0d0) then
      if (y.lt.sqeps) TERFC=1.0d0-2.0d0*x/sqrtpi
      if (y.ge.sqeps) TERFC=1.0d0-x*(1.0d0
     &        +DCSEVL(2.d0*x*x-1.d0,erfcs,nterf))
      else
      y=y*y
      if (y.le.4.d0) TERFC=exp(-y)/abs(x)*(0.5d0
     &        +DCSEVL((8.d0/y-5.d0)/3.d0,erc2cs,nterc2))
      if (y.gt.4.d0) TERFC=exp(-y)/abs(x)*(0.5d0
     &        +DCSEVL(8.d0/y-1.d0,erfccs,nterfc))
      if (x.lt.0.d0) TERFC=2.0d0-TERFC
      endif
      else
c        call mess('derfc','x so big erfc underflows',' ')
      TERFC=0.d0
      endif
      return
      END
      REAL*8 FUNCTION D1MACH(i)
C=======================================================================
C
c     return floating point machine dependent constants.
c     used by slatec library routines.
c
c      D1MACH( 1) = b**(emin-1), the smallest positive magnitude.
c      D1MACH( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c      D1MACH( 3) = b**(-t), the smallest relative spacing.
c      D1MACH( 4) = b**(1-t), the largest relative spacing.
c      D1MACH( 5) = log10(b)
c
c     assume double precision numbers are represented in the t-digit,
c     base-b form
c
c              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c     where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
c     emin .le. e .le. emax.
c
c     the values of b, t, emin and emax are provided in i1mach as
c     follows:
c
c      i1mach(10) = b, the base.
c      i1mach(14) = t, the number of base-b digits.
c      i1mach(15) = emin, the smallest exponent e.
c      i1mach(16) = emax, the largest exponent e.
c
c     to alter this function for a particular environment, the desired
c     set of data statements should be activated by removing the c from
c     column 1.  also, the values of D1MACH(1) - D1MACH(4) should be
c     checked for consistency with the local operating system.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION dmach(5)
      dmach(1)=2.22507385d-100
      dmach(2)=1.79769313d+100
      dmach(3)=1.11022302d-16
      dmach(4)=2.22044605d-16
      dmach(5)=0.30103001d+0
      D1MACH=dmach(i)
      return
      END
      INTEGER*4 FUNCTION INITDS(os,nos,eta)
C=======================================================================
C
c     determine the number of terms needed in an orthogonal
c     polynomial series so that it meets a specified accuracy.
c     from the slatec library fnlib INITDS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION os(*)
      err=0.0D+0
      do i=nos,1,-1
      err=err+DABS(    os(i) )
      if(err.le.eta) go to 10
      enddo
      i = 1
   10 INITDS=i
      return
      END
      REAL*8 FUNCTION DCSEVL(x,cs,n)
C=======================================================================
C
c     evaluate a chebyshev series.
c     from the slatec library fnlib DCSEVL.
C
C     WARNING - local variable onepl NEVER USED?
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION cs(*)
c
      b1=0.0d0
      b0=0.0d0
      twox=2.0d0*x
      do i=1,n
      b2=b1
      b1=b0
      ni=n+1-i
      b0=twox*b1-b2+cs(ni)
      enddo
      DCSEVL=0.5d0*(b0-b2)
      return
      END
      SUBROUTINE TABERFC
C=======================================================================
C
C     TABULATE ERFC AND CUBIC SPLINE COEFFICIENTS
C
C=======================================================================
      INCLUDE 'implicit.h'
      INCLUDE 'sigma1.h'
      DO 10 IX=-1,10002
      X = IX
      X = 1.0D-03*X
      ERFCX  (IX) = X
   10 ERFCTAB(IX) = TERFC(X)*DEXP(X*X)
C
C     FIRST DEFINE COMPLICATED CUBIC SPLINE COEFFICIENTS
C
C     F(X) = (X-X2)(X-X3)(X-X4)F1/[(X1-X2)(X1-X3)(X1-X4)] +
C            (X-X1)(X-X3)(X-X4)F2/[(X2-X1)(X2-X3)(X2-X4)] +
C            (X-X1)(X-X2)(X-X4)F3/[(X3-X1)(X3-X2)(X3-X4)] +
C            (X-X1)(X-X2)(X-X3)F4/[(X4-X1)(X4-X2)(X4-X3)]
C     DEFINE THE 4 COEFFICIENTS
C     ERFC1 = F1/[(X1-X2)(X1-X3)(X1-X4)]
C     ERFC2 = F2/[(X2-X1)(X2-X3)(X2-X4)]
C     ERFC3 = F3/[(X3-X1)(X3-X2)(X3-X4)]
C     ERFC4 = F4/[(X4-X1)(X4-X2)(X4-X3)]
C     X1-X2 = -D
C     X1-X3 = -D*2
C     X1-X4 = -D*3
c     X2-X1 = +D
C     X2-X3 = -D
C     X2-X4 = -D*2
C     X3-X1 = +D*2
C     X3-X2 = +D
C     X3-X4 = -D
C     X4-X1 = +D*3
C     X4-X2 = +D*2
C     X4-X3 = +D
C     ERFC1 = F1/[(-D   )(-D*2 )(-D*3 )] = -F1/[6*D^3]
C     ERFC2 = F2/[(+D   )(-D   )(-D*2 )] = +F2/[2*D^3]
C     ERFC3 = F3/[(+D*2 )(+D   )(-D   )] = -F3/[2*D^3]
C     ERFC4 = F4/[(+D*3 )(+D*2 )(+D   )] = +F4/[6*D^3]
C
      D   = 1.0D-03
      DDD = D**3
      DO 20 IX=0,10000
      ERFC1(IX) = -ERFCTAB(IX-1)/(6.0D+00*DDD)
      ERFC2(IX) =  ERFCTAB(IX  )/(2.0D+00*DDD)
      ERFC3(IX) = -ERFCTAB(IX+1)/(2.0D+00*DDD)
   20 ERFC4(IX) =  ERFCTAB(IX+2)/(6.0D+00*DDD)
C
C     NOW CONVERT TO SIMPLE CUBIC SPLINE COEFFICIENTS
C
C     F3 = [ERFC1+ERFC2+ERFC3+ERC4]
C     F2 =-D*[3*ERFC1+2*ERFC2+ERFC3]
C     F1 = D^2*[2*ERFC1-ERFC2-2*ERFC3-ERFC4]
C     F0 = 2*D^3*ERFC2
C     Z = X-X2
C     ERFC = (((F3*Z+F2)*Z+F1)*Z+F0)
C
      D = 1.0D-03
      DO 30 IX2=0,10000
      F3(IX2) = ERFC1(IX2)+ERFC2(IX2)+ERFC3(IX2)+ERFC4(IX2)
      F2(IX2) =-D*(3.0D+00*ERFC1(IX2)+
     1             2.0D+00*ERFC2(IX2)+
     1                     ERFC3(IX2))
      F1(IX2) =(D*D)*(2.0D+00*ERFC1(IX2) -
     1                        ERFC2(IX2) -
     2                2.0D+00*ERFC3(IX2) -
     3                        ERFC4(IX2))
   30 F0(IX2) =2.0*D*D*D*ERFC2(IX2)
      RETURN
      END
      subroutine groupiet(infile,outfile)
C=======================================================================
C
C     PROGRAM GROUPIE
C     ===============
C     VERSION 76-1 (NOVEMBER 1976)
C     VERSION 79-1 (OCTOBER 1979) CDC-7600 AND CRAY-1 VERSION.
C     VERSION 80-1 (MAY 1980) IBM, CDC AND CRAY VERSION
C     VERSION 81-1 (JANUARY 1981) EXTENSION TO 3000 GROUPS
C     VERSION 81-2 (MARCH 1981) IMPROVED SPEED
C     VERSION 81-3 (AUGUST 1981) BUILT-IN 1/E WEIGHTING SPECTRUM
C     VERSION 82-1 (JANUARY 1982) IMPROVED COMPUTER COMPATIBILITY
C     VERSION 83-1 (JANUARY 1983)*MAJOR RE-DESIGN.
C                                *ELIMINATED COMPUTER DEPENDENT CODING.
C                                *NEW, MORE COMPATIBLE I/O UNIT NUMBERS.
C                                *NEW MULTI-BAND LIBRARY BINARY FORMAT.
C     VERSION 83-2 (OCTOBER 1983) ADDED OPTION TO ALLOW SIGMA-0 TO BE
C                                 DEFINED EITHER AS MULTIPLES OF
C                                 UNSHIELDED TOTAL CROSS SECTION IN EACH
C                                 GROUP, OR POWERS OF 10 IN ALL GROUPS.
C     VERSION 84-1 (APRIL 1984)   ADDED MORE BUILT IN MULTIGROUP ENERGY
C                                 STRUCTURES.
C     VERSION 85-1 (APRIL 1985)  *UPDATED FOR ENDF/B-VI FORMATS.
C                                *SPECIAL I/O ROUTINES TO GUARANTEE
C                                 ACCURACY OF ENERGY.
C                                *DOUBLE PRECISION TREATMENT OF ENERGY
C                                 (REQUIRED FOR NARROW RESONANCES).
C                                *MINIMUM TOTAL CROSS SECTION TREATMENT
C     VERSION 85-2 (AUGUST 1985) *FORTRAN-77/H VERSION
C     VERSION 86-1 (JANUARY 1986)*ENDF/B-VI FORMAT
C     VERSION 86-2 (JUNE 1986)   *BUILT-IN MAXWELLIAN, 1/E AND FISSION
C                                 WEIGHTING SPECTRUM.
C     VERSION 88-1 (JULY 1988)   *OPTION...INTERNALLY DEFINE ALL I/O
C                                 FILE NAMES (SEE, SUBROUTINES FILIO1
C                                 FILIO2 FOR DETAILS).
C                                *IMPROVED BASED ON USER COMMENTS.
C     VERSION 89-1 (JANUARY 1989)*PSYCHOANALYZED BY PROGRAM FREUD TO
C                                 INSURE PROGRAM WILL NOT DO ANYTHING
C                                 CRAZY.
C                                *UPDATED TO USE NEW PROGRAM CONVERT
C                                 KEYWORDS.
C                                *ADDED LIVERMORE CIVIC COMPILER
C                                 CONVENTIONS.
C     VERSION 91-1 (JUNE 1991)   *INCREASED PAGE SIZE FROM 1002 TO 5010
C                                 POINTS
C                                *UPDATED BASED ON USER COMMENTS
C                                *ADDED FORTRAN SAVE OPTION
C                                *COMPLETELY CONSISTENT ROUTINE TO READ
C                                 FLOATING POINT NUMBERS.
C     VERSION 92-1 (JANUARY 1992)*ADDED RESONANCE INTEGRAL CALCULATION -
C                                 UNSHIELDED AND/OR SHIELDED - FOR
C                                 DETAILS SEE BELOW
C                                *INCREASED NUMBER OF ENERGY POINTS
C                                 IN BUILT-IN SPECTRA - TO IMPROVE
C                                 ACCURACY.
C                                *ALLOW SELECTION OF ZA/MF/MT OR
C                                 MAT/MF/MT RANGES - ALL DATA NOT
C                                 SELECTED IS SKIPPED ON INPUT AND
C                                 NOT WRITTEN AS OUTPUT.
C                                *COMPLETELY CONSISTENT I/O ROUTINES -
C                                 TO MINIMIZE COMPUTER DEPENDENCE.
C                                *NOTE, CHANGES IN INPUT PARAMETER
C                                 FORMAT - FOR ZA/MF/MT OR MAT/MF/MT
C                                 RANGES.
C     VERSION 92-2 (JUNE 1992)   *MULTIBAND PARAMETERS OUTOUT AS
C                                 CHARACTER (RATHER THAN BINARY) FILE.
C     VERSION 93-1 (APRIL 1993)  *INCREASED PAGE SIZE FROM 5010 TO
C                                 30000 POINTS
C                                *ELIMINATED COMPUTER DEPENDENCE.
C     VERSION 94-1 (JANUARY 1994)*VARIABLE ENDF/B DATA FILENAMES
C                                 TO ALLOW ACCESS TO FILE STRUCTURES
C                                 (WARNING - INPUT PARAMETER FORMAT
C                                 HAS BEEN CHANGED)
C                                *CLOSE ALL FILES BEFORE TERMINATING
C                                 (SEE, SUBROUTINE ENDIT)
C     VERSION 95-1 (JANUARY 1994)*CORRECTED MAXWELLIAN WEIGHTING
C                                *CHANGING WEIGHTING SPECTRUM FROM
C                                 0.1 TO 0.001 % UNCERTAINTY
C     VERSION 96-1 (JANUARY 1996) *COMPLETE RE-WRITE
C                                 *IMPROVED COMPUTER INDEPENDENCE
C                                 *ALL DOUBLE PRECISION
C                                 *ON SCREEN OUTPUT
C                                 *UNIFORM TREATMENT OF ENDF/B I/O
C                                 *IMPROVED OUTPUT PRECISION
C                                 *DEFINED SCRATCH FILE NAMES
C                                 *UP TO 1000 GROUP MULTI-BAND
C                                  CALCULATION (PREVIOUSLY 175)
C                                 *MAXIMUM NUMBER OF GROUPS REDUCED
C                                  FROM 3,000 TO 1,000
C                                 *UP TO 1000 MATERIALS
C                                  (PREVIOUSLY 100)
C                                 *CORRECTED USE OF MAXWELLIAN +
C                                  1/E + FISSION SPECTRUM
C                                 *ONLY 2 BAND VERSION DISTRIBUTED
C                                  (CONTACT AUTHOR FOR DETAILS)
C                                 *DEFINED SCRATCH FILE NAMES
C     VERSION 99-1 (MARCH 1999)   *CORRECTED CHARACTER TO FLOATING
C                                  POINT READ FOR MORE DIGITS
C                                 *UPDATED TEST FOR ENDF/B FORMAT
C                                  VERSION BASED ON RECENT FORMAT CHANGE
C                                 *GENERAL IMPROVEMENTS BASED ON
C                                  USER FEEDBACK
C     VERSION 99-2 (JUNE 1999)    *ASSUME ENDF/B-VI, NOT V, IF MISSING
C                                  MF=1, MT-451.
C     VERS. 2000-1 (FEBRUARY 2000)*ADDED MF=10, ACTIVATION CROSS SECTION
C                                  PROCESSING.
C                                 *GENERAL IMPROVEMENTS BASED ON
C                                  USER FEEDBACK
C     VERS. 2002-1 (FEBRUARY 2002)*ADDED TART 700 GROUP STRUCTURE
C                                 *ADDED VARIABLE SIGMA0 INPUT OPTION
C                  (MAY 2002)     *OPTIONAL INPUT PARAMETERS
C                  (NOV. 2002)    *ADDED SAND-II EXTENDED DOWN TO
C                                  1.0D-5 EV.
C                  (JUNE 2003)    *CORRECTED SAND-II 620 AND 640 GROUP
C                                  ENERGY BOUNDARIES DEFINITIONS.
C     VERS. 2004-1 (SEPT. 2004)  *INCREASED PAGE SIZE FROM 30000 TO
C                                 120000 POINTS
C                                *ADDED "OTHER" AS ADDITIONAL REACTION
C                                 TO IMPROVE MULTI-BAND FITTING
C                                *ADDED ITERATION FOR "BEST" PARTIAL
C                                 PARAMETERS.
C                                *DO NOT SKIP LOW TOTAL ENERGY RANGES
C                                 WHEN DEFINING AVERAGE CROSS SECTIONS -
C                                 THIS MAKES OUTPUT COMPATIBLE WITH
C                                 ANY STANDARD AVERAGING PROCEDURE
C     VERS. 2005-1 (JAN. 2005)   *ADDED OPTION TO CHANGE TEMPERATURE OF
C                                 BUILT-IN STANDARD SPECTRUM.
C     VERS. 2007-1 (JAN. 2007)   *CHECKED AGAINST ALL ENDF/B-VII.
C                                *INCREASED PAGE SIZE FROM 120,000 TO
C                                 600,000 POINTS
C     VERS. 2008-1 (JAN. 2008)   *72 CHARACTER FILE NAMES.
C                                *GENERAL UPDATES
C     VERS. 2010-1 (Apr. 2010)   *INCREASED WEIGHTING SPECTRUM TO 30,000
C                                 FROM 3,000 ENERGY POINTS.
C                                *ADDED OUTPUT TO PLOT/COMPARE SHIELDED
C                                 AND UNSHIELDED CROSS SECTIONS.
C     VERS. 2011-1 (June 2011)   *Corrected TART 700 groups to extend up
C                                 to 1 GeV (1,000 MeV) - previously it
C                                 was ERRONEOUSLY cutoff at 20 MeV.
C     VERS. 2011-2 (Nov. 2011)   *Corrected TART 616 groups lowest
C                                 energy from 1.0D-4 eV to 1.0D-5 eV.
C                                *Added TART 666 to 200 MeV (for TENDL).
C                                *Optional high energy cross section
C                                 extension above tabulated energy range
C                                 (either = 0 = standard, or constant)
C                                 WARNING - ENDF/B standard convention
C                                 is that the cross section = 0 where it
C                                 is not explicitly defined - extension
C                                 = 0 is standard, constant is NOT, so
C                                 constant extension is NOT RECOMMENDED.
C     VERS. 2012-1 (Aug. 2012)   *Added CODENAME
C                                *32 and 64 bit Compatible
C                                *Added ERROR stop.
C     VERS. 2013-1 (Nov. 2013)   *Extended OUT9.
C                                *Uses OUTG, not OUT10 for energies.
C     VERS. 2015-1 (Jan. 2015)   *Corrected SPECTM - handle ALL included
C                                 group structures, i.e., even those
C                                 that start above thremal range by
C                                 ALWAYS constructing weigthing spectrum
C                                 to be AT LEAST 1.0D-5 eV to 20 MeV.
C                                *Extended OUTG
C                                *Replaced ALL 3 way IF Statements.
C                                *Generalized TART Group Strructures.
C                                *Generalized SAND-II Group Structures.
C                                *Extended SAND-II to 60, 150, 200 MeV.
C     VERS. 2015-2 (Mar. 2015)   *Deleted 1P from formats reading input
C                                 parameters, causing incorrect scaling
C                                *Changed ALL data to "D" instead of
C                                 "E" to insure it is REAL*8 and avoid
C                                 Truncation ERRORS.
C     VERS. 2015-3 (July 2015)   *Insure no 10 digit output - not
C                                 needed for multi-group and this makes
C                                 listings simpler.
C                                *Corrected High Energy Extension =
C                                 Can effect highest energy group.
C     VERS. 2016-1 (July 2016)   *Added UKAEA 1102 Group Structure.
C                                *Increased storage to accommodate
C                                 much larger group structures =
C                                 up to 20,000 Groups.
C                                *Added output listing of the complete
C                                 input parameters for URRFIT, including
C                                 the NJOY parameters LSSF and ICOMP.
C                                *Changed multiple IF statements to
C                                 accommodate compiler optimizer
C                                *Cosmetic changes based on FREUD
C                                 psychoanalysis.
C                                *Updated multi-band treatment to
C                                 explcitly handle small shielding
C                                 limit - without this update the small
C                                 limit becomes numerically unstable.
C     VERS. 2017-1 (May  2017)   *Increased max. points to 3,000,000.
C                                *METHODB was incorrecctly named
C                                 METHOD  in one routine = corrected.
C                                *Default multi-band is method #2 =
C                                 conserve <x>, <1/(x+<x>>, <1/x>.
C                                *Definition of built-in group structure
C                                 using SUBROUTINE GROPE is identical
C                                 for GROUPIE and VIRGIN.
C                                *All floating input parameters changed
C                                 to character input + IN9 conversion.
C                                *Output report identfies MF now that
C                                 this code does more than just MF=3.
C                                *Added NRO = energy dependent scatter
C                                 radius to copying FILE2 parameters
C                                 to define unresolved energy range.
C
C     2015-2 Acknowledgment
C     =====================
C     I thank Chuck Whitmer (TerraPower,WA) and Andrej Trkov (NDS,IAEA)
C     for reporting the errors that led to the 2015-2 Improvements in
C     this code.
C
C     I thank Jean-Christophe Sublet (UKAEA) for contributing MAC
C     executables and Bojan Zefran (IJS, Slovenia) for contributing
C     LINUX (32 or 63 bit) executables. And most of all I must thank
C     Andrej Trkov (NDS, IAEA) for overseeing the entire PREPRO project
C     at IAEA, Vienna. This was a truly International team who worked
C     together to produce PREPRO 2015-2.
C
C     OWNED, MAINTAINED AND DISTRIBUTED BY
C     ------------------------------------
C     THE NUCLEAR DATA SECTION
C     INTERNATIONAL ATOMIC ENERGY AGENCY
C     P.O. BOX 100
C     A-1400, VIENNA, AUSTRIA
C     EUROPE
C
C     ORIGINALLY WRITTEN BY
C     ------------------------------------
C     Dermott E. Cullen
C
C     PRESENT CONTACT INFORMATION
C     ---------------------------
C     Dermott E. Cullen
C     1466 Hudson Way
C     Livermore, CA 94550
C     U.S.A.
C     Telephone  925-443-1911
C     E. Mail    RedCullen1@Comcast.net
C     Website    RedCullen1.net/HOMEPAGE.NEW
C
C     AUTHORS MESSAGE
C     ---------------
C     THE REPORT DESCRIBED ABOVE IS THE LATEST PUBLISHED DOCUMENTATION
C     FOR THIS PROGRAM. HOWEVER, THE COMMENTS BELOW SHOULD BE CONSIDERED
C     THE LATEST DOCUMENTATION INCLUDING ALL RECENT IMPROVEMENTS. PLEASE
C     READ ALL OF THESE COMMENTS BEFORE IMPLEMENTATION, PARTICULARLY
C     THE COMMENTS CONCERNING MACHINE DEPENDENT CODING.
C
C     AT THE PRESENT TIME WE ARE ATTEMPTING TO DEVELOP A SET OF COMPUTER
C     INDEPENDENT PROGRAMS THAT CAN EASILY BE IMPLEMENTED ON ANY ONE
C     OF A WIDE VARIETY OF COMPUTERS. IN ORDER TO ASSIST IN THIS PROJECT
C     IT WOULD BE APPECIATED IF YOU WOULD NOTIFY THE AUTHOR OF ANY
C     COMPILER DIAGNOSTICS, OPERATING PROBLEMS OR SUGGESTIONS ON HOW TO
C     IMPROVE THIS PROGRAM. HOPEFULLY, IN THIS WAY FUTURE VERSIONS OF
C     THIS PROGRAM WILL BE COMPLETELY COMPATIBLE FOR USE ON YOUR
C     COMPUTER.
C
C     PURPOSE
C     -------
C     THIS PROGRAM IS DESIGNED TO CALCULATE ANY COMBINATION OF
C     THE FOLLOWING QUANTITIES FROM LINEARLY INTERPOLABLE TABULATED
C     CROSS SECTIONS IN THE ENDF/B FORMAT
C
C     (1) UNSHIELDED GROUP AVERAGED CROSS SECTIONS
C     (2) BONDARENKO SELF-SHIELDED GROUP AVERAGED CROSS SECTIONS
C     (3) MULTI-BAND PARAMETERS
C
C     IN THE FOLLOWING FOR SIMPLICITY THE ENDF/B TERMINOLOGY--ENDF/B
C     TAPE--WILL BE USED. IN FACT THE ACTUAL MEDIUM MAY BE TAPE, CARDS,
C     DISK OR ANY OTHER MEDIUM.
C
C     ENDF/B FORMAT
C     -------------
C     THIS PROGRAM ONLY USES THE ENDF/B BCD OR CARD IMAGE FORMAT (AS
C     OPPOSED TO THE BINARY FORMAT) AND CAN HANDLE DATA IN ANY VERSION
C     OF THE ENDF/B FORMAT (I.E., ENDF/B-I, II,III, IV OR V FORMAT).
C
C     IT IS ASSUMED THAT THE DATA IS CORRECTLY CODED IN THE ENDF/B
C     FORMAT AND NO ERROR CHECKING IS PERFORMED. IN PARTICULAR IT IS
C     ASSUMED THAT THE MAT, MF AND MT ON EACH CARD IS CORRECT. SEQUENCE
C     NUMBERS (COLUMNS 76-80) ARE IGNORED ON INPUT, BUT WILL BE
C     CORRECTLY OUTPUT ON ALL CARDS. THE FORMAT OF SECTION MF=1, MT=451
C     AND ALL SECTIONS OF MF= 3 MUST BE CORRECT. THE PROGRAM COPIES ALL
C     OTHER SECTION OF DATA AS HOLLERITH AND AS SUCH IS INSENSITIVE TO
C     THE CORRECTNESS OR INCORRECTNESS OF ALL OTHER SECTIONS.
C
C     ALL FILE 3 CROSS SECTIONS THAT ARE USED BY THIS PROGRAM MUST BE
C     LINEARLY INTERPOLABLE IN ENERGY AND CROSS SECTION (ENDF/B
C     INTERPOLATION LAW 2). FILE 3 BACKGROUND CROSS SECTIONS MAY BE MADE
C     LINEARLY INTERPOLABLE USING PROGRAM LINEAR (UCRL-50400, VOL. 17,
C     PART A). THE RESONANCE CONTRIBUTION MAY BE ADDED TO THE BACKGROUND
C     CROSS SECTIONS USING PROGRAM RECENT (UCRL-50400, VOL. 17, PART B).
C     IF THIS PROGRAM FINDS THAT THE FILE 3 CROSS SECTIONS ARE NOT
C     LINEARLY INTERPOLABLE THIS PROGRAM WILL TERMINATE EXECUTION.
C
C     CONTENTS OF OUTPUT
C     ------------------
C     IF ENDF/B FORMATTED OUTPUT IS REQUESTED ENTIRE EVALUATIONS ARE
C     OUTPUT, NOT JUST THE MULTI-GROUPED FILE 3 CROSS SECTIONS, E.G.
C     ANGULAR AND ENERGY DISTRIBUTIONS ARE ALSO INCLUDED.
C
C     DOCUMENTATION
C     -------------
C     THE FACT THAT THIS PROGRAM HAS OPERATED ON THE DATA IS DOCUMENTED
C     BY THE ADDITION OF THREE COMMENT CARDS AT THE END OF EACH
C     HOLLERITH SECTION TO DESCRIBE THE GROUP STRUCTURE AND WEIGHTING
C     SPECTRUM, E.G.
C
C     ********************** PROGRAM GROUPIE (2017-1) ***************
C     UNSHIELDED GROUP AVERAGES USING   69 GROUPS (WIMS)
C     MAXWELLIAN, 1/E AND FISSION WEIGHTING SPECTRUM
C
C     THE ORDER OF ALL SIMILAR COMMENTS (FROM LINEAR, RECENT AND SIGMA1)
C     REPRESENTS A COMPLETE HISTORY OF ALL OPERATIONS PERFORMED ON
C     THE DATA.
C
C     THESE COMMENT CARDS ARE ONLY ADDED TO EXISTING HOLLERITH SECTIONS,
C     I.E., THIS PROGRAM WILL NOT CREATE A HOLLERITH SECTION. THE FORMAT
C     OF THE HOLLERITH SECTION IN ENDF/B-V DIFFERS FROM THE THAT OF
C     EARLIER VERSIONS OF ENDF/B. BY READING AN EXISTING MF=1, MT=451
C     IT IS POSSIBLE FOR THIS PROGRAM TO DETERMINE WHICH VERSION OF
C     THE ENDF/B FORMAT THE DATA IS IN. WITHOUT HAVING A SECTION OF
C     MF=1, MT=451 PRESENT IT IS IMPOSSIBLE FOR THIS PROGRAM TO
C     DETERMINE WHICH VERSION OF THE ENDF/B FORMAT THE DATA IS IN, AND
C     AS SUCH IT IS IMPOSSIBLE FOR THE PROGRAM TO DETERMINE WHAT FORMAT
C     SHOULD BE USED TO CREATE A HOLLERITH SECTION.
C
C     REACTION INDEX
C     --------------
C     THIS PROGRAM DOES NOT USE THE REACTION INDEX WHICH IS GIVEN IN
C     SECTION MF=1, MT=451 OF EACH EVALUATION.
C
C     THIS PROGRAM DOES NOT UPDATE THE REACTION INDEX IN MF=1, MT=451.
C     THIS CONVENTION HAS BEEN ADOPTED BECAUSE MOST USERS DO NOT
C     REQUIRE A CORRECT REACTION INDEX FOR THEIR APPLICATIONS AND IT WAS
C     NOT CONSIDERED WORTHWHILE TO INCLUDE THE OVERHEAD OF CONSTRUCTING
C     A CORRECT REACTION INDEX IN THIS PROGRAM. HOWEVER, IF YOU REQUIRE
C     A REACTION INDEX FOR YOUR APPLICATIONS, AFTER RUNNING THIS PROGRAM
C     YOU MAY USE PROGRAM DICTIN TO CREATE A CORRECT REACTION INDEX.
C
C     SECTION SIZE
C     ------------
C     SINCE THIS PROGRAM USES A LOGICAL PAGING SYSTEM THERE IS NO LIMIT
C     TO THE NUMBER OF POINTS IN ANY SECTION, E.G., THE TOTAL CROSS
C     SECTION MAY BE REPRESENTED BY 200,000 DATA POINTS.
C
C     SELECTION OF DATA
C     -----------------
C     THE PROGRAM SELECTS MATERIALS TO BE PROCESSED BASED EITHER ON
C     MAT (ENDF/B MAT NO.) OR ZA. THE PROGRAM ALLOWS UP TO 100 MAT OR
C     ZA RANGES TO BE SPECIFIED. THE PROGRAM WILL ASSUME THAT THE
C     ENDF/B TAPE IS IN EITHER MAT OR ZA ORDER, WHICHEVER CRITERIA IS
C     USED TO SELECT MATERIALS, AND WILL TERMINATE WHEN A MAT OR ZA
C     IS FOUND THAT IS ABOVE THE RANGE OF ALL REQUESTS.
C
C     ENERGY ORDER AND UNITS
C     ----------------------
C     ALL ENERGIES (FOR CROSS SECTIONS, WEIGHTING SPECTRUM OR GROUP
C     BOUNDARIES) MUST BE IN UNITS OF EV AND MUST BE IN ASCENDING
C     NUMERICAL ORDER.
C
C     ENERGY GRID
C     -----------
C     ALTHOUGH ALL REACTIONS MUST TO LINEARLY INTERPOLABLE, THEY DO NOT
C     ALL HAVE TO USE THE SAME ENERGY GRID. EACH REACTION CAN BE GIVEN
C     BY AN INDEPENDENT ENERGY GRID. THIS PROGRAM WILL PROCEED FROM
C     THE LOWEST TO HIGHEST ENERGY SELECTING EACH ENERGY INTERVAL OVER
C     WHICH ALL DATA, FOR ANY GIVEN CALCULATION, ARE ALL LINEARLY
C     INTERPOLABLE.
C
C     GROUP STRUCTURE
C     ---------------
C     THIS PROGRAM IS DESIGNED TO USE AN ARBITRARY ENERGY GROUP
C     STRUCTURE WHERE THE ENERGIES ARE IN EV AND ARE IN INCREASING
C     ENERGY ORDER. THE MAXIMUM NUMBER OF GROUPS IS 20,000.
C
C     THE USER MAY INPUT AN ARBITRARY GROUP STRUCTURE OR THE USER MAY
C     USE USE ONE OF THE SEVEN BUILT-IN GROUP STRUCTURES.
C     (0) 175 GROUP (TART STRUCTURE)
C     (1)  50 GROUP (ORNL STRUCTURE)
C     (2) 126 GROUP (ORNL STRUCTURE)
C     (3) 171 GROUP (ORNL STRUCTURE)
C     (4) 620 GROUP (SAND-II STRUCTURE, UP TO 18 MEV)
C     (5) 640 GROUP (SAND-II STRUCTURE, UP TO 20 MEV)
C     (6)  69 GROUP (WIMS STRUCTURE)
C     (7)  68 GROUP (GAM-I STRUCTURE)
C     (8)  99 GROUP (GAM-II STRUCTURE)
C     (9)  54 GROUP (MUFT STRUCTURE)
C    (10)  28 GROUP (ABBN STRUCTURE)
C    (11) 616 GROUP (TART STRUCTURE TO 20 MeV)
C    (12) 700 GROUP (TART STRUCTURE TO 1 GEV)
C    (13) 665 GROUP (SAND-II STRUCTURE, 1.0D-5 eV, UP TO 18 MEV)
C    (14) 685 GROUP (SAND-II STRUCTURE, 1.0D-5 eV, UP TO 20 MEV)
C    (15) 666 GROUP (TART STRUCTURE TO 200 MeV)
C    (16) 725 GROUP (SAND-II STRUCTURE, 1.0D-5 eV, UP TO  60 MEV)
C    (17) 755 GROUP (SAND-II STRUCTURE, 1.0D-5 eV, UP TO 150 MEV)
C    (18) 765 GROUP (SAND-II STRUCTURE, 1.0D-5 eV, UP TO 200 MEV)
C    (19)1102 GROUP (UKAEA   STRUCTURE, 1.0D-5 eV, UP TO   1 GeV)
C
C     GROUP AVERAGES
C     --------------
C     THIS PROGRAM DEFINES GROUP AVERAGED CROSS SECTIONS AS...
C
C               (INTEGRAL E1 TO E2) (SIGMA(E)*S(E)*WT(E)*DE)
C     AVERAGE = -----------------------------------------
C               (INTEGRAL E1 TO E2) (S(E)*WT(E)*DE)
C     WHERE...
C
C     AVERAGE  = GROUP AVERAGED CROSS SECTION
C     E1, E2   = ENERGY LIMITS OF THE GROUP
C     SIGMA(E) = ENERGY DEPENDENT CROSS SECTION FOR ANY GIVEN REACTION
C     S(E)     = ENERGY DEPENDENT WEIGHTING SPECTRUM
C     WT(E)    = ENERGY DEPENDENT SELF-SHIELDING FACTOR.
C
C     ENERGY DEPENDENT WEIGHTING SPECTRUM
C     -----------------------------------
C     THE ENERGY DEPENDENT WEIGHTING SPECTRUM IS GIVEN BY AN ARBITRARY
C     TABULATED LINERLY INTERPOLABLE FUNCTION WHICH CAN BE DESCRIBED
C     BY AN ARBITRARY NUMBER OF POINTS. THIS ALLOWS THE USER TO
C     SPECIFY ANY DESIRED WEIGHTING SPECTRUM TO ANY GIVEN DEGREE OF
C     ACCURACY. REMEMBER THAT THE PROGRAM WILL ASSUME THAT THE SPECTRUM
C     IS LINEARLY INTERPOLABLE BETWEEN TABULATED POINTS. THEREFORE THE
C     USER SHOULD USE ENOUGH POINTS TO INSURE AN ADEQUATE REPRESENTATION
C     OF THE SPECTRUM BETWEEN TABULATED DATA POINTS.
C
C     THE PRESENT VERSION OF THE CODE HAS THREE BULIT-IN WEIGHTING
C     SPECTRA,
C
C     (1) CONSTANT
C     (2) 1/E
C     (3) MAXWELLIAN = E*EXP(-E/KT)/KT                (0.0 TO 4*KT)
C         1/E        = C1/E                           (4*KT TO 67 KEV)
C         FISSION    = C2*EXP(-E/WA)*SINH(SQRT(E*WB)) (ABOVE 67 KEV)
C
C         KT     = 0.253 EV (293 KELVIN)
C         WA     = 9.65D+5
C         WB     = 2.29D-6
C         C1, C2 = DEFINED TO MAKE SPECTRUM CONTINUOUS
C
C         FISSION SPECTRUM CONSTANTS FROM
C         A.F.HENRY, NUCLEAR REACTOR ANALYSIS, P. 11, MIT PRESS (1975)
C
C     UNSHIELDED GROUP AVERAGES
C     -------------------------
C     FOR UNSHIELDED AVERAGES THE SELF-SHIELDING FACTOR (WT(E)) IS SET
C     TO UNITY. THIS PROGRAM ALLOWS UP TO 20,000 GROUPS.
C
C     SELF-SHIELDED GROUP AVERAGES
C     ----------------------------
C     IF SELF-SHIELDED AVERAGES AND/OR MULTI-BAND PARAMETERS ARE
C     CALCULATED THIS PROGRAM ALLOWS UP TO 20,000 GROUPS. SELF-SHIELDED
C     AVERAGES AND/OR MULTI-BAND PARAMETERS ARE CALCULATED FOR THE
C     TOTAL, ELASTIC, CAPTURE AND FISSION.
C
C     FOR THE TOTAL, ELASTIC, CAPTURE AND FISSION THE PROGRAM USES A
C     WEIGHTING FUNCTION THAT IS A PRODUCT OF THE ENERGY DEPENDENT
C     WEIGHTING SPECTRUM TIMES A BONDERENKO TYPE SELF-SHIELDING FACTOR.
C
C     WT(E) = S(E)/(TOTAL(E)+SIGMA0)**N
C
C     WHERE...
C
C     S(E)     - ENERGY DEPENDENT WEIGHTING SPECTRUM (DEFINED BY
C                TABULATED VALUES AND LINEAR INTERPOLATION BETWEEN
C                TABULATED VALUES).
C     TOTAL(E) - ENERGY DEPENDENT TOTAL CROSS SECTION FOR ONE MATERIAL
C                (DEFINED BY TABULATED VALUES AND LINEAR INTERPOLATION
C                BETWEEN TABULATED VALUES).
C     SIGMA0   - CROSS SECTION TO REPRESENT THE EFFECT OF ALL OTHER
C                MATERIALS AND LEAKAGE (DEFINED WITHIN EACH GROUP TO BE
C                A MULTIPLE OF THE UNSHIELDED TOTAL CROSS SECTION WITHIN
C                THAT GROUP OR POWERS OF 10 - INPUT OPTION).
C     N        - A POSITIVE INTEGER (0, 1, 2 OR 3).
C
C     THE PROGRAM WILL USE ONE ENERGY DEPENDENT WEIGHTING SPECTRUM S(E)
C     AND 25 DIFFERENT BONDERENKO TYPE SELF-SHIELDING FACTORS (25 SIGMA0
C     AND N COMBINATIONS) TO DEFINE 25 DIFFERENT AVERAGE CROSS SECTIONS,
C     FOR EACH REACTION, WITHIN EACH GROUP.
C
C     THE 25 WEIGHTING FUNCTIONS USED ARE....
C     (1)   - UNSHIELDED CROSS SECTIONS (N=0)
C     (2-22)- PARTIALLY SHIELDED CROSS SECTIONS (N=1 ,VARIOUS SIGMA0)
C             THE VALUES OF SIGMA0 USED WILL BE EITHER,
C             (A) THE VALUES OF SIGMA0 THAT ARE USED VARY FROM 1024
C             TIMES THE UNSHIELDED TOTAL CROSS SECTIONS IN STEPS OF 1/2
C             DOWN TO 1/1024 TIMES THE UNSHIELDED TOTAL CROSS SECTION
C             (A RANGE OF OVER 1 MILLION, CENTERED ON THE UNSHIELDED
C             TOTAL CROSS SECTION WITHIN EACH GROUP).
C             (B) THE SAME CONSTANT VALUES OF SIGMA0 IN EACH GROUP. THE
C             VALUES OF SIGMA0 USED INCLUDE 40000, 20000, 10000, 7000,
C             4000, 2000, 1000, 700, 400, 200, 100, 70, 40, 20, 10, 7,
C             4, 2, 1, 0.7, 0.4 (A RANGE OF 100,000 SPANNING MORE THAN
C             THE RANGE OF SIGMA0 VALUES THAT MAY BE ENCOUNTERED IN
C             ACTUAL APPLICATIONS)
C     (23)  - TOTALLY SHIELDED FLUX WEIGHTED CROSS SECTION
C             (N=1, SIGMA0=0)
C     (24)  - TOTALLY SHIELDED CURRENT WEIGHTED CROSS SECTION
C             (N=2, SIGMA0=0)
C     (25)  - TOTALLY SHIELDED COSINE SQUARED WEIGHTED CROSS SECTION
C             (N=3, SIGMA0=0)
C
C     FOR ALL OTHER REACTIONS (EXCEPT TOTAL, ELASTIC, CAPTURE AND
C     FISSION) THE PROGRAM WILL USE THE ENERGY DEPENDENT WEIGHTING
C     SPECTRUM S(E) TO DEFINE THE UNSHIELDED (BONDERENKO N=0)
C     AVERAGED CROSS SECTION WITHIN EACH GROUP.
C
C     CALCULATION OF RESONANCE INTEGRALS
C     ----------------------------------
C     IN A PURE ELASTIC ISOTROPICALLY SCATTERING MATERIAL WITH A
C     CONSTANT CROSS SECTION THE SPECTRUM WILL BE 1/E AND THERE WILL
C     BE NO SELF-SHIELDING.
C
C     IN THIS CASE IF THE CROSS SECTION VARIES WITH ENERGY THE
C     SPECTRUM WILL STILL BE 1/E AND THE SELF-SHIELDING FACTOR WILL
C     BE EXACTLY 1/SIG-TOT(E) - WHERE SIG-TOT(E) = SIG-EL(E), SINCE
C     THERE IS ONLY SCATTERING.
C
C     IF WE HAVE AN INFINITELY DILUTE AMOUNT OF A MATERIAL UNIFORMLY
C     MIXED WITH A PURE ELASTIC ISOTROPICALLY SCATTERING MATERIAL WITH
C     A CONSTANT CROSS SECTION THE STANDARD DEFINITION OF THE RESONANCE
C     INTEGRAL CAN BE USED TO DEFINE REACTION RATES FOR EACH REACTION.
C
C     THE RESONANCE INTEGRAL IS DEFINED AS,
C
C     RI      = (INTEGRAL E1 TO E2) (SIGMA(E)*S(E)*WT(E)*DE)
C
C     WHERE NORMALLY,
C     S(E)    = 1/E
C     WT(E)   = 1    - NO SELF-SHIELDING
C
C     FROM THE ABOVE DEFINITION OF GROUP AVERAGED CROSS SECTIONS THE
C     RESONANCE INTEGRAL IS,
C
C     RI      = AVERAGE * (INTEGRAL E1 TO E2) (S(E)*WT(E)*DE)
C
C     FOR A 1/E SPECTRUM AND NO SELF-SHIELDING THIS REDUCES TO,
C
C     RI      = AVERAGE* LOG(E2/E1)
C
C     IN ANY OTHER SITUATION, INCLUDING ABSORPTION AND/OR ENERGY
C     DEPENDENT CROSS SECTIONS, THE SPECTRUM WILL NOT BE 1/E -
C     ABSORPTION WILL TEND TO DECREASE THE SPECTRUM PROGRESSIVELY
C     MORE AT LOWER ENERGIES - ENERGY DEPENDENCE OF THE CROSS SECTION
C     WILL LEAD TO SELF-SHIELDING.
C
C     HERE WE WILL NOT ATTEMPT TO PERFORM A DETAILED SPECTRUM
C     CALCULATION TO ACCOUNT FOR ABSORPTION.
C
C     HOWEVER, WE WILL EXTEND THE DEFINITION OF THE RESONANCE INTEGRAL
C     TO ACCOUNT FOR SELF-SHIELDING EFFECTS BY ALLOWING FOR INCLUSION
C     OF SELF-SHIELDING EFFECTS IN THE DEFINITION OF GROUP AVERAGES
C     AND THEN DEFINING THE RESONANCE INTEGRAL AS,
C
C     RI      = AVERAGE* LOG(E2/E1)
C
C     IN ORDER TO CALCULATE RESONANCE INTEGRALS YOU MUST FOLLOW THESE
C     STEPS,
C
C     1) SELECT A 1/E SPECTRUM - ON FIRST LINE OF INPUT PARAMETERS.
C     2) SELECT THE ENERGY BOUNDARIES - NORMALLY ONLY 1 GROUP FROM
C        0.5 EV UP TO 20 MEV - HOWEVER, YOU ARE FREE TO SELECT ANY
C        ENERGY RANGE THAT YOU WISH - YOU MAY EVEN SELECT MORE THAN
C        1 GROUP MERELY BY SPECIFYING MORE THAN 1 GROUP AS INPUT -
C        THIS CAN BE USED TO DEFINE THE CONTRIBUTIONS TO THE RESONANCE
C        INTEGRAL FROM INDIVIDUAL ENERGY RANGES.
C     3) SELECT THIS OPTION FOR THE UNSHIELDED AND/OR SHIELDED OUTPUT
C        LISTING - ON THE SECOND LINE OF INPUT PARAMETERS.
C
C     WHEN THIS OPTION IS USED THE PROGRAM WILL CALCULATE GROUP AVERAGED
C     CROSS SECTIONS - AS DEFINED ABOVE - PRIOR TO OUTPUT THE RESULTS
C     WILL MERELY BE MULTIPLIED BY THE WIDTH OF THE GROUP ASSUMING YOU
C     HAVE SELECTED A 1/E SPECTRUM - THERE IS NO CHECK ON THIS - THE
C     PROGRAM MERELY MULTIPLIES THE GROUP AVERAGED CROSS SECTIONS BY,
C
C     LOG(E2/E1) - WHERE E2 AND E1 ARE THE GROUP ENERGY BOUNDARIES.
C
C     WARNING - IT IS UP TO YOU TO INSURE THAT YOU FOLLOW EXACTLY THE
C               STEPS OUTLINED ABOVE IF YOU WISH TO OBTAIN MEANINGFUL
C               RESULTS.
C
C     NOTE - OUTPUT IN THE ENDF/B FORMAT IS ALWAYS GROUP AVERAGED CROSS
C            SECTIONS, REGARDLESS OF WHETHER YOU ASK FOR AVERAGED CROSS
C            SECTIONS OR RESONANCE INTEGRALS - THIS IS BECAUSE DATA IN
C            THE ENDF/B FORMAT IS EXPLICITLY DEFINED TO BE CROSS
C            SECTIONS.
C
C            RESONANCE INTEGRAL OUTPUT CAN ONLY BE OBTAINED IN THE
C            LISTING FORMATS.
C
C     MINIMUM TOTAL CROSS SECTION TREATMENT
C     -------------------------------------
C     SINCE THE BONDARENKO SELF-SHIELDING DEPENDS ON 1/TOTAL CROSS
C     SECTION, THE ALGORITHM WILL BECOME NUMERICALLY UNSTABLE IF THE
C     TOTAL CROSS SECTION IS NEGATIVE (AS OCCURS IN MANY ENDF/B
C     EVALUATIONS). IF THE TOTAL IS LESS THAN SOME MINIMUM ALLOWABLE
C     VALUE (DEFINE BY OKMIN, PRESENTLY 1 MILLI-BARN) AN ERROR MESSAGE
C     WILL BE PRINTED AND FOR THE SELF-SHIELDING CALCULATION ALL ENERGY
C     INTERVALS IN WHICH THE TOTAL IS LESS THAN THE MINIMUM WILL BE
C     IGNORED.
C
C     NOTE, FOR THE UNSHIELDED CALCULATIONS ALL CROSS SECTIONS WILL BE
C     CONSIDERED WHETHER THEY ARE POSITIVE OR NEGATIVE. THEREFORE IF
C     THE TOTAL CROSS SECTION IS NEGATIVE OR LESS THAN THE MINIMUM
C     VALUE THERE MAY BE AN INCONSISTENCY BETWEEN THE UNSHIELDED AND
C     THE SELF-SHIELDED CROSS SECTIONS. IF THE TOTAL CROSS SECTION IS
C     NEGATIVE AND SELF-SHIELDED CROSS SECTIONS ARE CALCULATED THE
C     PROGRAM WILL PRINT AN ERROR MESSAGE INDICATING THAT THE SELF-
C     SHIELDED RESULTS ARE UNRELIABLE AND SHOULD NOT BE USED. THEREFORE
C     IN THIS CASE THE PROGRAM WILL NOT ATTEMPT TO MODIFY THE UNSHIELDED
C     RESULTS TO ELIMINATE THE EFFECT OF NEGATIVE CROSS SECTIONS, SINCE
C     THE UNSHIELDED RESULTS ARE THE ONLY ONES WHICH TRULY REFLECT THE
C     ACTUAL INPUT.
C
C     RESOLVED RESONANCE REGION
C     -------------------------
C     IN THE RESOLVED RESONANCE REGION (ACTUALLY EVERYWHERE BUT IN THE
C     UNRESOLVED RESONANCE REGION) THE CROSS SECTIONS OUTPUT BY LINEAR-
C     RECENT-SIGMA1 WILL BE ACTUAL ENERGY DEPENDENT CROSS SECTIONS AND
C     THE CALCULATIONS BY THIS PROGRAM WILL YIELD ACTUAL SHIELDED AND
C     UNSHIELDED CROSS SECTIONS.
C
C     UNRESOLVED RESONANCE REGION
C     ---------------------------
C     IN THE UNRESOLVED RESONANCE REGION PROGRAM RECENT USES THE
C     UNRESOLVED RESONANCE PARAMETERS TO CALCULATE INFINITELY DILUTE
C     AVERAGE CROSS SECTIONS. THIS PROGRAM WILL MERELY READ THIS
C     INFINITELY DILUTE DATA AS IF IT WERE ENERGY DEPENDENT DATA AND
C     GROUP AVERAGE IT. AS SUCH THIS PROGRAM WILL PRODUCE THE CORRECT
C     UNSHIELDED CROSS SECTION IN THE UNRESOLVED RESONANCE REGION, BUT
C     IT WILL NOT PRODUCE THE CORRECT SELF-SHIELDING EFFECTS.
C
C     ACCURACY OF RESULTS
C     -------------------
C     ALL INTEGRALS ARE PERFORMED ANALYTICALLY. THEREFORE NO ERROR IS
C     INTRODUCED DUE TO THE USE OF TRAPAZOIDAL OR OTHER INTEGRATION
C     SCHEME. THE TOTAL ERROR THAT CAN BE ASSIGNED TO THE RESULTING
C     AVERAGES IS JUST THAT DUE TO THE ERROR IN THE CROSS SECTIONS
C     AND ENERGY DEPENDENT WEIGHTING SPECTRUM. GENERALLY SINCE THE
C     THE ENERGY DEPENDENT WEIGHTING SPECTRUM APPEARS IN BOTH THE
C     NUMERATOR AND THE DENOMINATOR THE AVERAGES RAPIDLY BECOME
C     INSENSITIVE TO THE WEIGHTING SPECTRUM AS MORE GROUPS ARE USED.
C     SINCE THE WEIGHTING SPECTRUM IS LOADED IN THE PAGING SYSTEM THE
C     USER CAN DESCRIBE THE SPECTRUM TO ANY REQUIRED ACCURACY USING
C     ANY NUMBER OF ENERGY VS. SPECTRUM PAIRS.
C
C     MULTI-BAND PARAMETERS
C     ---------------------
C     MULTI-BAND PARAMETERS ARE CALCULATED FOR THE TOTAL, ELASTIC,
C     CAPTURE AND FISSION REACTIONS. WITH THE NUMBER OF GROUPS THAT
C     ARE NORMALLY USED (SEE BUILT IN GROUP STRUCTURES) ALL OTHER
C     REACTIONS RESULT IN A NEGLIGABLE AMOUNT OF SELF-SHIELDING. AS
C     SUCH THEIR EQUIVALENT BAND CROSS SECTION WILL MERELY BE THEIR
C     UNSHIELDED VALUE WITHIN EACH BAND.
C
C     FOR ANY GIVEN EVALUATION, WITHIN ANY GIVEN GROUP THIS PROGRAM
C     WILL GENERATE THE MINIMUM NUMBER OF BANDS REQUIRED WITHIN THAT
C     GROUP. AS OUTPUT TO THE COMPUTER READABLE DISK FILE THE BAND
C     PARAMETERS FOR EACH EVALUATION WILL BE FORMATTED TO HAVE THE
C     SAME NUMBER OF BANDS IN ALL GROUPS (WITH ZERO WEIGHT FOR SOME
C     BANDS WITHIN ANY GROUP). THE USER MAY DECIDE TO HAVE OUTPUT
C     EITHER WITH THE MINIMUM NUMBER OF BANDS REQUIRED FOR EACH
C     EVALUATION (E.G. 2 BANDS FOR HYDROGEN AND 4 BANDS FOR U-233) OR
C     THE SAME NUMBER OF BANDS FOR ALL EVALUATIONS (E.G. 4 BANDS FOR
C     BOTH HYDROGEN AND U-233).
C
C     FOR 2 OR FEWER BANDS THE PROGRAM USES AN ANALYTIC EXPRESSION
C     TO DEFINE ALL MULTI-BAND PARAMETERS. FOR MORE THAN 2 BANDS THE
C     PROGRAM PERFORMS A NON-LINEAR FIT TO SELECT THE MULTI-BAND
C     PARAMETERS THAT MINIMIZE THE MAXIMUM FRACTIONAL ERROR AT ANY
C     POINT ALONG THE ENTIRE SELF-SHIELDING CURVE. THE NUMBER OF BANDS
C     REQUIRED WITHIN ANY GIVEN GROUP IS DEFINED BY INSURING THAT THE
C     MULTI-BAND PARAMETERS CAN BE USED TO ACCURATELY DEFINE SELF-
C     SHIELDED CROSS SECTIONS ALONG THE ENTIRE SELF-SHIELDING CURVE
C     FROM SIGMA0 = 0 TO INFINITY. THE USER MAY DEFINE THE ACCURACY
C     REQUIRED.
C
C     ENDF/B FORMATTED UNSHIELDED AVERAGES
C     ------------------------------------
C     UNSHIELDED MULTI-GROUP AVERAGED CROSS SECTIONS FOR ALL REACTIONS
C     MAY BE OBTAINED IN THE ENDF/B FORTRAN IN EITHER HISTOGRAM
C     (INTERPOLATION LAW 1) OR LINEARLY INTERPOLABLE (INTERPOLATION
C     LAW 2) FORM. SEE INPUT BELOW FOR DETAILS.
C
C     MIXTURES OF MATERIALS AND RESONANCE OVERLAP
C     -------------------------------------------
C     THE SELF-SHIELDED CROSS SECTIONS FOR THE INDIVIDUAL CONSTITUENTS
C     OF ANY MIXTURE CAN BE CALCULATED BY THIS PROGRAM BY REALIZING THAT
C     THIS PROGRAM ESSENTIALLY ONLY USES THE TOTAL CROSS SECTION AS A
C     WEIGHTING FUNCTION TO ACCOUNT FOR SELF-SHIELDING EFFECTS. FOR A
C     MIXTURE IT IS THEREFORE ONLY NECESSARY TO USE THE TOTAL CROSS
C     SECTION FOR THE MIXTURE IN PLACE OF THE ACTUAL TOTAL CROSS SECTION
C     FOR EACH CONSTITUENT AND TO RUN THIS PROGRAM. THIS CAN BE DONE BY
C     FIRST RUNNING PROGRAM MIXER TO CALCULATE THE ENERGY DEPENDENT
C     TOTAL CROSS SECTION FOR ANY COMPOSITE MIXTURE. NEXT, SUBSTITUTE
C     THIS COMPOSITE TOTAL CROSS SECTION FOR THE ACTUAL TOTAL CROSS
C     SECTION OF EACH CONSTITUENT (IN EACH ENDF/B FORMATTED EVALUATION).
C     FINALLY, RUN THIS PROGRAM TO CALCULATE THE SELF-SHIELDED CROSS
C     SECTION FOR EACH CONSTITUENT, PROPERLY ACCOUNTING FOR RESONANCE
C     OVERLAP BETWEEN THE RESONANCES OF ALL OF THE CONSTITUENTS OF THE
C     MIXTURE. DURING THE SAME RUN THESE SELF-SHIELDED CROSS SECTIONS
C     CAN IN TURN BE USED TO CALCULATE FULLY CORRELATED MULT-BAND
C
C     MULTI-BAND PARAMETER OUTPUT FORMAT
C     ----------------------------------
C     FOR VERSIONS 92-2 AND LATER VERSIONS THE MULTI-BAND PARAMETERS
C     ARE OUTPUT IN A SIMPLE CHARACTER FORMAT, THAT CAN BE TRANSFERRED
C     AND USED ON VIRTUALLY ANY COMPUTER.
C
C     THE BINARY FORMAT USED IN EARLIER VERSIONS OF THIS CODE IS NO
C     LONGER USED.
C
C     CONTACT THE AUTHOR IF YOU WOULD LIKE TO RECEIVE A SIMPLE PROGRAM
C     TO READ THE CHARACTER FORMATTED MULTI-BAND PARAMETER FILE AND
C     CREATE A BINARY, RANDOM ACCESS FILE FOR USE ON VIRTUALLY ANY
C     COMPUTER.
C
C     THE FORMAT OF THE CHARACTER FILE IS,
C
C     RECORD   COLUMNS   FORMAT   DESCRIPTION
C        1       1-72     18A4    LIBRARY DESCRIPTION (AS READ)
C        2       1-11      I11    MATERIAL ZA
C               12-22      I11    NUMBER GROUPS
C               23-33      I11    NUMBER OF BANDS
C               34-44     E11.4   TEMPERATURE (KELVIN)
C               45-55    1X,10A1  HOLLERITH DESCRIPTION OF ZA
C        3       1-11     E11.4   ENERGY (EV) - GROUP BOUNDARY.
C               12-22     E11.4   TOTAL      (FIRST BAND)
C               23-33     E11.4   ELASTIC
C               34-44     E11.4   CAPTURE
C               35-55     E11.4   FISSION
C        4       1-11     -----   BLANK
C               12-22     E11.4   TOTAL      (SECOND BAND)
C               23-33     E11.4   ELASTIC
C               34-44     E11.4   CAPTURE
C               35-55     E11.4   FISSION
C
C     LINES 3 AND 4 ARE REPEATED FOR EACH GROUP. THE LAST LINE FOR EACH
C     MATERIAL (ZA) IS,
C
C        N       1-11     E11.4   ENERGY (EV) - UPPER ENERGY LIMIT OF
C                                               LAST GROUP.
C
C     FOR EXAMPLE, A 175 GROUP, 2 BAND FILE, FOR EACH MATERIAL WILL
C     CONTAIN 352 LINES = 1 HEADER LINE, 175 * 2 LINES OF PARAMETERS,
C                         AND 1 FINAL LINE WITH THE UPPER ENERGY LIMIT
C                         OF THE LAST GROUP.
C
C     INPUT FILES
C     -----------
C     UNIT  DESCRIPTION
C     ----  -----------
C       2   INPUT DATA (BCD - 80 CHARACTERS/RECORD)
C      10   ORIGINAL ENDF/B DATA (BCD - 80 CHARACTERS/RECORD)
C
C     OUTPUT FILES
C     ------------
C     UNIT  DESCRIPTION
C     ----  -----------
C      31   MULTI-BAND PARAMETERS CHARACTER FILE - OPTIONAL
C           (BCD - 80 CHARACTERS/RECORD)
C      32   SELF-SHIELDED CROSS SECTION LISTING - OPTIONAL
C           (BCD - 120 CHARACTERS/RECORD)
C      33   MULTI-BAND PARAMETER LISTING - OPTIONAL
C           (BCD - 120 CHARACTERS/RECORD)
C      34   UNSHIELDED CROSS SECTION LISTING - OPTION
C           (BCD - 120 CHARACTERS/RECORD)
C       3   OUTPUT REPORT (BCD - 80 CHARACTERS/RECORD)
C      11   MULTI-GROUP ENDF/B DATA - OPTIONAL
C           (BCD - 80 CHARACTERS/RECORD)
C
C     SCRATCH FILES
C     -------------
C     UNIT  FILENAME  DESCRIPTION
C     ----  --------  -----------
C       8   ENERGY DEPENDENT WEIGHTING SPECTRUM
C           (BINARY - 40080 WORDS/BLOCK)
C       9   TOTAL CROSS SECTION
C           (BINARY - 40080 WORDS/BLOCK)
C      12   ELASTIC CROSS SECTION - ONLY FOR SELF-SHIELDING CALCULATION
C           (BINARY - 40080 WORDS/BLOCK)
C      13   CAPTURE CROSS SECTION - ONLY FOR SELF-SHIELDING CALCULATION
C           (BINARY - 40080 WORDS/BLOCK)
C      14   FISSION CROSS SECTION - ONLY FOR SELF-SHIELDING CALCULATION
C           (BINARY - 40080 WORDS/BLOCK)
C
C     OPTIONAL STANDARD FILE NAMES (SEE SUBROUTINES FILIO1 AND FILIO2)
C     ----------------------------------------------------------------
C     UNIT  FILE NAME
C     ----  ----------
C       2   GROUPIE.INP
C       3   GROUPIE.LST
C       8   (SCRATCH)
C       9   (SCRATCH)
C      10   ENDFB.IN
C      11   ENDFB.OUT
C      12   (SCRATCH)
C      13   (SCRATCH)
C      14   (SCRATCH)
C      31   MULTBAND.TAB
C      32   SHIELD.LST
C      33   MULTBAND.LST
C      34   UNSHIELD.LST
C
C      I/O UNITS USED
C      --------------
C      UNITS 2, 3 8, 9 AND 10 WILL ALWAYS BE USED.
C      UNITS 31 THROUGH 34 AND 11 ARE OPTIONALLY USED DEPENDING ON THE
C      OUTPUT REQUESTED.
C      UNITS 12, 13 AND 14 WILL ONLY BE USED IF SELF-SHIELDED OR
C      MULTIBAND OUTPUT IS REQUESTED.
C
C     INPUT CARDS
C     -----------
C     CARD  COLS.  FORMAT  DESCRIPTION
C     ----  -----  ------  -----------
C       1    1-11    I11   SELECTION CRITERIA (0=MAT, 1=ZA)
C       1   12-22    I11   NUMBER OF GROUPS.
C                          =.GT.0 - ARBITRARY GROUP BOUNDARIES ARE READ
C                                   FROM INPUT FILE (N GROUPS REQUIRE
C                                   N+1 GROUP BOUNDARIES). CURRENT
C                                   PROGRAM MAXIMUM IS 20,000 GROUPS.
C                                   BUILT-IN OPTIONS INCLUDE....
C                          =  0   - TART    175 GROUPS
C                          = -1   - ORNL     50 GROUPS
C                          = -2   - ORNL    126 GROUPS
C                          = -3   - ORNL    171 GROUPS
C                          = -4   - SAND-II 620 (665) GROUPS TO 18 MEV
C                          = -5   - SAND-II 640 (685) GROUPS TO 20 MEV
C                          = -6   - WIMS     69 GROUPS
C                          = -7   - GAM-I    68 GROUPS
C                          = -8   - GAM-II   99 GROUPS
C                          = -9   - MUFT     54 GROUPS
C                          =-10   - ABBN     28 GROUPS
C                          =-11   - TART    616 GROUPS TO 20 MEV
C                          =-12   - TART    700 GROUPS TO 1 GEV
C                          =-13   - SAND-II 665 GROUPS TO 18 MEV
C                          =-14   - SAND-II 685 GROUPS TO 20 MEV
C                          =-15   - TART    666 GROUPS TO 200 MEV
C                          =-16   - SAND-II 725 GROUPS TO 60 MEV
C                          =-17   - SAND-II 755 GROUPS TO 150 MEV
C                          =-18   - SAND-II 765 GROUPS TO 200 MEV
C                          =-19   - UKAEA  1102 GROUPS TO   1 GeV
C       1   23-33    I11   MULTI-BAND SELECTOR
C                          =  0 - NO MULTI-BAND CALCULATIONS
C                          =  1 - 2 BAND. CONSERVE AV(TOT), AV(1/TOT)
C                                 AND AV(1/TOT**2)
C                          =  2 - 2 BAND. CONSERVE AV(TOT), AV(1/TOT)
C                                 AND AV(1/(TOT+SIGMA0)) WHERE
C                                 SIGMA0 = AV(TOT) IN EACH GROUP
C                          = 3-5- MULTI-BAND FIT. CONSERVE AV(TOT) AND
C                                 MINIMIZE FRACTIONAL ERROR FOR ENTIRE
C                                 SELF-SHIELDING CURVE (SIGMA0 = 0 TO
C                                 INFINITY)
C                          IF THE SELECTOR IS POSITIVE (1 TO 5) THE
C                          MINIMUM NUMBER OF BANDS WILL BE OUTPUT FOR
C                          EACH ISOTOPE INDEPENDENTLY. IF THE SELECTOR
C                          IS NEGATIVE (-1 TO -5) THE SAME NUMBER OF
C                          BANDS (ABS(SELECTOR)) WILL BE OUTPUT FOR
C                          ALL ISOTOPES.
C       1   34-44    I11   NUMBER OF POINTS USED TO DESCRIBE ENERGY
C                          DEPENDENT WEIGHTING SPECTRUM S(E).
C                          = -2    - MAXWELLIAN - UP TO 0.1 EV
C                                    1/E        - 0.1 EV TO 67 KEV
C                                    FISSION    - ABOVE 67 KEV
C 05/01/20-----------------ADDED OPTION TO ALLOW TEMPERATURE OF THE
C                          MAXWELLIAN TO BE CHANGED - SEE INPUT LINE 4,
C                          COLUMNS 55 - 66.
C                          = -1    - 1/E
C                          = 0 OR 1- ENERGY INDEPENDENT (SO CALLED FLAT
C                                    WEIGHTING SPECTRUM).
C                          = .GT.1 - READ THIS MANY POINTS FROM INPUT
C                                    TO DESCRIBE WEIGHTING SPECTRUM.
C                                    NO LIMIT TO THE NUMBER OF POINTS
C                                    USED TO DESCRIBE WEIGHTING.
C       1   45-55   E11.4  MULTI-BAND CONVERGENCE CRITERIA.
C                          ONLY USED FOR 3 OR MORE BANDS. THE NUMBER OF
C                          BANDS IN EACH GROUPS IS SELECTED TO INSURE
C                          THAT THE ENTIRE SELF-SHIELDING CURVE CAN BE
C                          REPRODUCED TO WITHIN THIS FRACTIONAL ERROR.
C                          = .LT. 0.0001 - USE STANDARD 0.001
C                                          (0.1 PER-CENT)
C                          = .GE. 0.0001 - USE AS CONVERGENCE CRITERIA
C       1   56-66    I11   SIGMA-0 DEFINITION SELECTOR.
C                          < 0 - 21 VALUES OF SIGMA0 ARE READ INPUT AND
C                                INTERPRETED AS FIXED VALUES = SAME AS
C                                = 1 DESCRIPTION BELOW
C                                INPUT VALUES MUST ALL BE,
C                                1) GREATER THAN 0
C                                2) IN DESCENDING VALUE ORDER
C                          = 0 - SIGMA-0 WILL BE DEFINED AS A MULTIPLE
C                                OF THE UNSHIELDED TOTAL CROSS SECTION
C                                IN EACH GROUP (VALUES OF 1/1024 TO
C                                1024 IN STEPS OF A FACTOR OF 2 WILL
C                                BE USED AS THE MULTIPLIER).
C                          = 1 - SIGMA-0 WILL BE DEFINED AS THE SAME
C                                NUMBER OF BARNS IN EACH GROUP (VALUES
C                                40000 TO 0.4 BARNS WILL BE USED. WITHIN
C                                EACH DECADE VALUES OF 10, 7, 4, 2, 1
C                                BARNS WILL BE USED).
C       1   67-70    I4    High energy extension = definition of cross
C                          section above highest tabulated energy.
C                          = 0 = cross section = 0 (standard ENDF/B)
C                          = 1 = cross section = constant (equal to
C                                value at highest tabulated energy).
C     2-4    1-66 6E11.4   IF SIGMA-0 DEFINITION SELECTOR < 0, THE NEXT
C                          4 LINES OF INPUT ARE THE 22 VALUES OF SIGMA0,
C                          6 PER LINE.
C       2    1-72    A72   ENDF/B INPUT DATA FILENAME
C                          (STANDARD OPTION = ENDFB.IN)
C       3    1-72    A72   ENDF/B OUTPUT DATA FILENAME
C                          (STANDARD OPTION = ENDFB.OUT)
C
C     THE FOURTH INPUT CARD IS USED TO SELECT ALL DESIRED OUTPUT MODES.
C     EACH OUTPUT DEVICE MAY BE TURNED OFF (0) OR ON (1). THEREFORE
C     THEREFORE EACH OF THE FOLLOWING INPUT PARAMETERS MAY BE EITHER
C     ZERO TO INDICATE NO OUTPUT OR NON-ZERO TO INDICATE OUTPUT.
C
C       4     1-11   I11   SELF-SHIELDED CROSS SECTION LISTING
C                          = 1 - CROSS SECTIONS
C                          = 2 - RESONANCE INTEGRALS
C       4    12-22   I11   MULTI-BAND PARAMETER LISTING
C       4    23-33   I11   MULTI-BAND PARAMETERS COMPUTER READABLE
C       4    34-44   I11   UNSHIELDED CROSS SECTIONS IN ENDF/B FORMAT
C                          = 1 - HISTOGRAM FORMAT (INTERPOLATION LAW 1)
C                          = 2 - LINEAR-LINEAR (INTERPOLATION LAW 2)
C       4    45-55   I11   UNSHIELDED CROSS SECTIONS LISTING
C                          = 1 - CROSS SECTIONS
C                          = 2 - RESONANCE INTEGRALS
C 05/01/20 - ADDED THE BELOW OPTION
C       4    56-66   E11.4 IF THE STANDARD BUILT-IN SPECTRA IS USED,
C                          INPUT LINE 1, COLUMNS 34-44 = 2, THIS FIELD
C                          CAN BE USED TO OPTIONALLY CHANGE TEMPERATURE
C                          OF THE MAXWELLIAN.
C                          INPUT IS IN EV (0.0253 EV = ROOM TEMPERATURE)
C                          = 0 - USE DEFAULT 0.0253 EV, ROOM TEMPERATURE
C                          > 0 - USE THIS AS THE TEMPERATURE
C                          RESTRICTION - TEMPERATURE CANNOT EXCEED
C                          1000 EV.
C
C       5     1-80   18A4  LIBRARY IDENTIFICATION. ANY TEXT THAT THE
C                          USER WISHES TO IDENTIFY THE MULTI-BAND
C                          PARAMETERS. THIS LIBRARY IDENTIFICATION IS
C                          WRITTEN INTO THE COMPUTER READABLE MULTI-BAND
C                          DATA FILE.
C
C      6-N    1- 6    I6   LOWER MAT OR ZA LIMIT
C             7- 8    I2   LOWER MF LIMIT
C             9-11    I3   LOWER MT LIMIT
C            12-17   I11   UPPER MAT OR ZA LIMIT
C            18-19    I2   UPPER MF LIMIT
C            20-22    I3   UPPER MT LIMIT
C                          UP TO 100 RANGES MAY BE SPECIFIED, ONE RANGE
C                          PER LINE. THE LIST OF RANGES IS TERMINATED
C                          BY A BLANK CARD. IF THE UPPER MAT OR ZA
C                          LIMIT IS LESS THAN THE LOWER LIMIT THE UPPER
C                          IS SET EQUAL TO THE LOWER LIMIT. IF THE UPPER
C                          MF OR MT LIMIT IS ZERO IT WILL BE SET EQUAL
C                          TO ITS MAXIMUM VALUE, 99 OR 999, RESPECTIVELY
C                          IF THE FIRST REQUEST LINE IS BLANK IT WILL
C                          TERMINATE THE LIST OF REQUESTS AND CAUSE ALL
C                          DATA TO BE RETRIEVED (SEE EXAMPLE INPUT).
C
C      VARY   1-66  6E11.4 ENERGY GROUP BOUNDARIES. ONLY REQUIRED IF
C                          THE NUMBER OF GROUPS INDICATED ON THE FIRST
C                          INPUT CARD IS POSITIVE. ALL ENERGIES MUST
C                          BE IN ASCENDING ENERGY IN EV. THE PRESENT
C                          LIMITS ARE 1 TO 20,000 GROUPS. FOR N GROUPS
C                          N+1 BOUNDARIES WILL BE READ FROM THE
C                          INPUT FILE, E.G. IF THE FIRST INPUT CARD
C                          INDICATES 20 GROUPS, 21 ENERGY BOUNDARIES
C                          WILL BE READ FROM THE INPUT FILE.
C
C      VARY   1-66  6E11.4 ENERGY DEPENDENT WEIGHTING SPECTRUM. ONLY
C                          REQUIRED IF THE NUMBER OF POINTS INDICATED
C                          ON FIRST CARD IS MORE THAN ONE. DATA IS
C                          GIVEN IN (ENERGY, WEIGHT) PAIRS, UP TO 3
C                          PAIRS PER CARD, USING ANY NUMBER OF CARDS
C                          REQUIRED. ENERGIES MUST BE IN ASCENDING
C                          ORDER IN EV. THE SPECTRUM VALUES MUST BE
C                          NON-NEGATIVE. THE ENERGY RANGE OF SPECTRUM
C                          MUST AT LEAST SPAN THE ENERGY RANGE OF THE
C                          ENERGY GROUPS. SINCE SPECTRUM IS STORED IN
C                          PAGING SYSTEM THERE IS NO LIMIT TO NUMBER
C                          OF POINTS THAT CAN BE USED TO DESCRIBE THE
C                          WEIGHTING SPECTRUM.
C
C     EXAMPLE INPUT NO. 1
C     -------------------
C     REQUEST DATA BY MAT AND PROCESS ALL DATA (ALL MAT BETWEEN 1 AND
C     9999). USE THE TART 175 GROUP STRUCTURE, GENERATE 2 BAND
C     PARAMETERS (THE FOR ALL ISOTOPES) TO 0.1 PER-CENT ACCURACY
C     IN THE SELF-SHIELDING CURVE. OUTPUT ALL  LISTING, COMPUTER
C     READABLE AND ENDF/B FORMAT GROUP AVERAGES.
C
C     EXPLICITLY SPECIFY THE STANDARD FILENAMES.
C
C     THE FOLLOWING 7 INPUT LINES ARE REQUIRED.
C
C          0          0         -2          0 1.00000-03          0
C ENDFB.IN
C ENDFB.OUT
C          1          1          1          1          1
C TART 175 GROUP, 2 BAND LIBRARY TO 0.1 PER-CENT ACCURACY
C     1 1  1  9999 0  0
C                       (BLANK CARD TERMINATES REQUEST LIST)
C
C     EXAMPLE INPUT NO. 2
C     -------------------
C     THE SAME EXAMPLE 1, AS ABOVE, ONLY THE ENDF/B DATA WILL BE READ
C     FROM \ENDFB6\SIGMA1\K300\ZA092238 (U-238 AT 300 KELVIN) AND
C     WRITTEN TO \ENDFB6\GROUPIE\K300\ZA092238
C
C     THE FOLLOWING 7 INPUT LINES ARE REQUIRED.
C
C          0          0         -2          0 1.00000-03          0
C \ENDFB6\SIGMA1\K300\ZA092238
C \ENDFB6\GROUPIE\K300\ZA092238
C          1          1          1          1          1
C TART 175 GROUP, 2 BAND LIBRARY TO 0.1 PER-CENT ACCURACY
C     1 1  1  9999 0  0
C                       (BLANK CARD TERMINATES REQUEST LIST)
C
C     EXAMPLE INPUT NO. 3
C     -------------------
C     PROCESS ALL DATA. USE 1/E WEIGHTING IN ORDER TO CALCULATE
C     UNSHIELDED ONE GROUP CROSS SECTIONS OVER THE ENERGY RANGE 0.5 EV
C     TO 1 MEV (NOTE THAT THE RESULTS ARE SIMPLY PROPORTIONAL TO THE
C     RESONANCE INTEGRAL FOR EACH REACTION). OUTPUT UNSHIELDED LISTING.
C
C     LEAVE THE DEFINITION OF THE FILENAMES BLANK - THE PROGRAM WILL
C     THEN USE STANDARD FILENAMES.
C
C     THE FOLLOWING 7 INPUT CARDS ARE REQUIRED.
C
C          0          0          1         -1                     0
C                       (USE STANDARD FILENAME = ENDFB.IN)
C                       (USE STANDARD FILENAME = ENDFB.OUT)
C          0          0          0          0          1
C RESONANCE INTEGRAL CALCULATION (FROM 0.5 EV TO 1 MEV)
C                       (RETRIEVE ALL DATA, TERMINATE REQUEST LIST)
C 5.00000-01 1.00000+06
C
C=======================================================================
      INCLUDE 'implicit.h'
C-----08/08/2012 DEFINE CODE NAME
      CHARACTER*8 CODENAME
      COMMON/NAMECODE/CODENAME
      CHARACTER*1 ZABCD,FIELDX
      CHARACTER*4 FMTHOL,TUNITS,TMPHOL,CARD
CAK
      CHARACTER*72 infile,outfile
CAK
      INTEGER*4 OUTP,OTAPE,OTAPE2,TYPSIG
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
CAK   COMMON/FILLER/ISECT
      COMMON/FILLERg/ISECT
      COMMON/COPC/CARD(17)
      COMMON/LEADER/TEMP,Q,L1,L2,N1,N2,MAT,MF,MT
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
CAK   COMMON/TEMPO/TMPTAB(3),TEMP1,NTEMP
      COMMON/TEMPOg/TMPTAB(3),TEMP1,NTEMP
      COMMON/TEMPC/TINOUT(4)
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
      COMMON/ELPASD/TUNITS(2,4),TMPHOL(3)
      COMMON/ELPASZ/ZABCD(10)
      COMMON/ELPAS2/EGB,TMPK
CAK   COMMON/PAGER/NPAGE,NPAGM1
      COMMON/PAGERg/NPAGE,NPAGM1
      COMMON/COPI/MFIELD(3)
      COMMON/GRABER/NOSELF,LSECT
CAK   COMMON/HOLFMT/FMTHOL
      COMMON/HOLFMTg/FMTHOL
      COMMON/MINOK/OKMIN(5),REDMIN(5)
      COMMON/TYPER/NTYPE
      COMMON/VERSES/TEMP3,IVERSE
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      COMMON/LOGCOM/MYLOG(8)
      COMMON/RESCOM/RES1,RES2,URES1,URES2,LSSF,ICOMP,INPART
      INCLUDE 'groupie.h'
c-----2016/5/21 - Increased NBT and INT dimension from 20 to 100
      DIMENSION NPTUSE(5),NPTTOT(5),NBT(100),INT(100),MTNEED(5,3)
C-----DEFINE REQUIRED REACTIONS FOR NEUTRONS AND PHOTONS.
      DATA MTNEED/0,  1  ,2,102, 18,
     1            0,501,522,502,504,
     1            0,  0,  0,  0,  0/
C-----DEFINE LAST ZA READ FROM ENDF/B FORMAT FILE.
      DATA LASTZA/0/
C=======================================================================
C
C     INITIALIZATION.
C
C=======================================================================
C-----08/08/2012 DEFINE CODE NAME
      CODENAME = 'GROUPIE '
C-----INITIALIZE TIME
CAK   CALL TIMER
      CALL TIMERpre
C-----DEFINE ALL I/O UNITS AND OPTIONALLY DEFINE FILE NAMES.
CAK   CALL FILIO1
      CALL FILIO1s
C-----LOAD ALL DATA INTO COMMON
      CALL BLOCK
C-----DEFINE PAGE SIZE.
      NPAGE=MAXPOINT
      NPAGM1=NPAGE-1
C-----DEFINE MINIMUM ACCEPTABLE SPECTRUM AND CROSS SECTION VALUES AND
C-----INITIALIZE POINT COUNTS FOR ALL REACTIONS.
      DO 10 ISECT=1,5
      OKMIN(ISECT)=0.0
   10 NPTTOT(ISECT)=0
C-----INITIALIZE COUNT OF SECTIONS
      DO I=1,8
      MYLOG(I) = 0
      ENDDO
C-----DEFINE MINIMUM ALLOWABLE CROSS SECTION FOR TOTAL AND ELASTIC
      OKMIN(2)=0.001
      OKMIN(3)=0.001
C-----INITIALIZE COUNT OF EVALUATIONS FOR WHICH MULTIBAND PARAMETERS
C-----HAVE BEEN GENERATED.
      LZA=0
C-----Initialize Resonance Region Boundaries.
      INPART= 1     ! Initialize to neutron (in case not ENDF6 format)
      RES1  = 0.0d0
      RES2  = 0.0d0
      URES1 = 0.0d0
      URES2 = 0.0d0
      LSSF  = 0
      ICOMP = 0
C=======================================================================
C
C     READ AND CHECK ALL INPUT PARAMETERS AND COPY TAPE LABEL.
C
C=======================================================================
CAK   CALL READIN
      CALL READINg(infile,outfile)
C-----READ ENDF/B TAPE LABEL AND IF THERE IS MULTI-GROUP OUTPUT IN THE
C-----ENDF/B FORMAT WRITE LABEL TO MULTI-GROUP OUTPUT FILE.
      CALL COPYL
C-----LIST TAPE LABEL.
      WRITE(OUTP,650) CARD,MFIELD(1)
      WRITE(*   ,650) CARD,MFIELD(1)
C=======================================================================
C
C     SEARCH FOR NEXT REQUESTED MATERIAL. TEST FOR END OF DATA.
C
C=======================================================================
CAK20 CALL NXTMAT
   20 CALL NXTMATg
      IF(MATH.LE.0) GO TO 420
C
C     SECTION HEAD CARD FOUND. LOCATE FILE 1, SECTION 451 OR FILE 3.
C
      IF(MFH.NE.1.OR.MTH.NE.451) GO TO 30
C-----FILE 1 , SECTION 451. ADD COMMENTS AND COPY REMAINDER OF FILE 1.
CAK   CALL FILE1
      CALL FILE1g
      GO TO 20
C
C     FILE2 - Define Resonance Region Boundaries
C
   30 IF(MFH.NE.2.OR.MTH.NE.151) GO TO 40
      CALL FILE2
      GO TO 20
C=======================================================================
C
C     SEARCH FOR MF = 3, 10 OR 23 - OTHERWISE COPY.
C
C=======================================================================
   40 IF(MFH.eq. 3) go to 50   ! Neutron Data
      IF(MFH.eq.10) go to 50   ! Activation Data
      IF(MFH.eq.23) go to 50   ! Photon data
C-----NOT FILE 1, 2, 3 OR 23. COPY SECTION (MT).
      CALL CONTOG
      CALL COPYSG    ! Copy to end of Section (MT=0)
      GO TO 20
C=======================================================================
C
C     FILE 3, 10 OR 23 FOUND. LOCATE REQUIRED REACTIONS.
C
C=======================================================================
C-----SAVE MF FOR OUTPUT REPORT.
   50 MFHOUT = MFH
C-----DEFINE TEMPERATURE FROM FIRST SECTION OF EACH MAT.
      IF(IZA.EQ.LASTZA) GO TO 60
      LASTZA=IZA
C-----DEFINE CHEMICAL SYMBOL FOR THIS MAT.
      CALL ZAHOL(IZA,ZABCD)
C-----SET FLAG TO INDICATE NEW MAT...SELECT TEMPERATURE FROM FIRST
C-----SECTION USED, WHEN NEXT CARD IS READ.
   60 MPASS=0
C-----INITIALIZE POINT COUNTS FOR NEXT MATERIAL.
      DO 70 ISECT=2,5
      NPTUSE(ISECT)=0
   70 NPTAB(ISECT)=0
      NPTUSE(1)=0
C-----INITIALIZE MINIMUM VALUES CROSS SECTIONS READ.
      DO 80 I=1,5
   80 REDMIN(I)=2.0*OKMIN(I)
C-----SET FLAG TO INDICATE IF SELF-SHIELDING CALCULATION MUST BE
C-----PERFORMED AND HAS NOT YET BEEN DONE ( 0 -YES, 1 -NO).
      IMDONE=0
      IF(NOSELF.NE.0) IMDONE=1
C=======================================================================
C
C     BEGINNING OF REACTION LOOP.
C
C=======================================================================
C-----IS SECTION REQUIRED FOR SELF-SHIELDING CALCULATION.
   90 MYTYPE=1
      IF(MFNOW.EQ.23) MYTYPE=2
      IF(MFNOW.EQ.10) MYTYPE=3
      DO 100 ISECT=2,5
      IF(MTH.eq.MTNEED(ISECT,MYTYPE)) go to 130
  100 CONTINUE
C-----THIS IS NOT TOTAL, ELASTIC, CAPTURE OR FISSION. SET POINT
C-----COUNT INDEX.
      NTYPE=1
C-----IF UNSHIELDING OUTPUT PERFORM CALCULATION. OTHERWISE SKIP SECTION.
      IF(OTAPE.GT.0.OR.LIST3.GT.0) GO TO 140
C=======================================================================
C
C     REACTION LOOP
C
C=======================================================================
C-----SKIP TO BEGINNING OF NEXT SECTION.
  110 CALL SKIPS
  120 CALL CONTIG
      IF(MTH.gt.0) go to 90
C-----CHECK FOR END OF FILE 3, 10 OR 23.
      IF(MFH.le.0) go to 250
      go to 120
C=======================================================================
C
C     REQUIRED REACTION FOUND. DEFINE INDEX TO SAVE POINT COUNT
C     OF THIS SECTION (NTYPE). IF ONLY DOING UNSHIELDED CALCULATION
C     SET INDEX (ISECT) TO LOAD DATA INTO PAGE 2.
C
C=======================================================================
  130 NTYPE=ISECT
      IF(NOSELF.le.0) go to 150
C-----LOAD DATA INTO DEFAULT PAGE (2 - IF NO SELF-SHIELDING,
C----- 4 - OTHERWISE).
  140 ISECT=LSECT
C=======================================================================
C
C     LOOP OVER TABLES = 1 TABLE, EXCEPT FOR MF=10
C
C=======================================================================
  150 NS=1
      IF(MFH.EQ.10) NS = N1H
C-----04/11/00 - MOVED CONTOG OUTPUT HERE - FROM TAB1
      IF(OTAPE.GT.0) CALL CONTOG
      DO 240 LOOP=1,NS
C
C     READ AND CHECK INTERPOLATION LAW. THEN LOAD SECTION INTO PAGING
C     SECTION.
C
      CALL CARDI(TEMP,Q,L1,L2,N1,N2)
C-----DEFINE TEMPERATURE FROM FIRST SECTION READ, OR FILE 1 IF ENDF/B-VI
C-----FORMAT.
      IF(MPASS.GT.0) GO TO 180
      MPASS=1
      TEMP1=TEMP
      IF(IVERSE.GE.6) TEMP1=TEMP3
C-----SELECT OUTPUT TEMPERATURE UNITS.
      DO 160 NTEMP=1,3
      IF(TEMP1.LE.TMPTAB(NTEMP)) GO TO 170
  160 CONTINUE
      NTEMP=4
  170 TMPK=TEMP1*TINOUT(NTEMP)
      TMPHOL(2)=TUNITS(1,NTEMP)
      TMPHOL(3)=TUNITS(2,NTEMP)
  180 CALL TERPI(NBT,INT,N1)
      NPTAB(ISECT)=N2
C-----FOR TOTAL, ELASTIC, CAPTURE OR FISSION SAVE POINT COUNT. FOR ALL
C-----OTHERS ADD UP THE TOTAL NUMBER OF POINTS.
      IF(NTYPE.NE.1) GO TO 190
      NPTUSE(NTYPE)=NPTUSE(NTYPE)+N2
      GO TO 200
  190 NPTUSE(NTYPE)=N2
C-----INSURE INTERPOLATION LAW IS LINEAR-LINEAR.
  200 DO 210 I=1,N1
      IF(INT(I).NE.2) GO TO 220
  210 CONTINUE
      GO TO 230
C-----INTERPOLATION LAW IS NOT LINEAR-LINEAR. WRITE ERROR MESSAGE AND
C-----TERMINATE.
  220 WRITE(OUTP,540) MATH,MFH,MTH,N1
      WRITE(OUTP,550) (NBT(I),INT(I),I=1,N1)
      WRITE(OUTP,560)
      CALL ENDERROR
C-----READ SECTION AND LOAD INTO PAGING SYSTEM.
  230 CALL PAGIN5(ITAPE,ISECT,NTYPE)
C
C     UNSHIELDED CROSS SECTION CALCULATION.
C
C-----UNSHIELDED TOTAL AVERAGES ARE ALWAYS PERFORMED. IF ENDF/B OR LIST
C-----FORMATTED OUTPUT IS REQUESTED UNSHIELDED AVERAGES WILL ALSO BE
C-----PERFORMED FOR ALL OTHER REACTIONS.
      IF(MTH.NE.1.AND.OTAPE.LE.0.AND.LIST3.LE.0) GO TO 240
C-----CALCULATE UNSHIELDED GROUP AVERAGES.
      CALL GROUPU
C-----IF NOT TOTAL, ELASTIC, CAPTURE OR FISSION RESET SAVED POINT
C-----COUNT TO ZERO.
      IF(MTH.NE.MTNEED(ISECT,MYTYPE)) NPTAB(ISECT)=0
C=======================================================================
C
C     END OF TABLE LOOP.
C
C=======================================================================
  240 CONTINUE
C-----04/11/00 - MOVED SEND HERE FROM TAB1.
      IF(OTAPE.GT.0) CALL OUTS(MATH,MFH)
C-----IF ALL REQUIRED SECTIONS READ PERFORM SELF-SHIELDING AND
C-----MULTI-BAND CALCULATION (LAST SECTION REQUIRED FOR SELF-SHIELDING
C-----CALCULATION IS CAPTURE, MT=102).
      IF(MTH.lt.102) go to 110
C
C     SELF-SHIELDING AND MULTI-BAND CALCULATIONS.
C
C-----IF SELF-SHIELDING CALCULATION ALREADY DONE OR NOT REQUIRED SKIP
C-----THIS SECTION.
  250 IF(IMDONE.GT.0) GO TO 270
      IMDONE=1
C-----DEFINE NUMBER OF SECTIONS ACTUALLY READ.
      NSECT=5
      IF(NPTAB(5).GT.0) GO TO 260
      NSECT=4
      IF(NPTAB(4).GT.0) GO TO 260
      NSECT=3
      IF(NPTAB(3).GT.0) GO TO 260
      NSECT=2
C-----PERFORM SELF-SHIELDING AND MULTI-BAND CALCULATIONS.
  260 CALL GROUPS
C-----IF STILL IN FILE 3, 10 OR 23 BRANCH BACK TO PROCESS NEXT SECTION.
  270 IF(MFH.GT.0) GO TO 110
C=======================================================================
C
C     END OF MF - IF END MF=3 PRINT RESONANCE REGION SUMMARY.
C
C=======================================================================
      IF(MFHOUT.eq.3) then
      if(RES2.gt.0.0d0) then
      CALL OUT9G( RES1,FIELDX(1,1))
      CALL OUT9G( RES2,FIELDX(1,2))
      WRITE(3,280) ZABCD,((FIELDX(k,kk),k=1,11),kk=1,2)
      WRITE(*,280) ZABCD,((FIELDX(k,kk),k=1,11),kk=1,2)
  280 format(1x,78('-')/' Resonance Region for ',10A1/
     1 ' Resolved Resonance Region.....',11A1,' to ',11A1, ' eV')
      else
      WRITE(3,290) ZABCD
      WRITE(*,290) ZABCD
  290 format(1x,78('-')/' Resonance Region for ',10A1/
     1 ' Resolved Resonance Region..... NONE')
      endif
c-----Unresolved
      if(URES2.gt.0.0d0) then
      CALL OUT9G( ATWT,FIELDX(1,3))
      WRITE(3,300) IZA,MATNOW,(FIELDX(k,3),k=1,11)
      WRITE(*,300) IZA,MATNOW,(FIELDX(k,3),k=1,11)
  300 format(' Complete Input Parameters for URRFIT are listed below'/
     1 ' ZA............................',I11/
     2 ' MAT...........................',I11/
     3 ' Atomic Weight Ratio (ATWT)....',11A1)
      CALL OUT9G(URES1,FIELDX(1,5))
      CALL OUT9G(URES2,FIELDX(1,6))
      WRITE(3,310) ((FIELDX(k,kk),k=1,11),kk=5,6)
      WRITE(*,310) ((FIELDX(k,kk),k=1,11),kk=5,6)
  310 format(
     1 ' Unresolved Resonance Region...',11A1,' to ',11A1, ' eV')
      WRITE(3,320) ICOMP,LSSF
      WRITE(*,320) ICOMP,LSSF
  320 format(
     1 ' Unresolved Competition (ICOMP)',I11/
     2 ' Is Unresolved Tabulated (LSSF)',I11/
     3       1x,78('-'))
      else
      WRITE(3,330)
      WRITE(*,330)
  330 format(
     1 ' Unresolved Resonance Region... NONE'/
     2       1x,78('-'))
      endif
      endif
C=======================================================================
C
C     END OF MAT. PRINT SUMMARY OF MATERIAL AND ERROR MESSAGES.
C
C=======================================================================
C-----SUMMARY OF MATERIAL.
      CALL OUT9G(TEMP1,FIELDX(1,1))
      WRITE(OUTP,640) ZABCD,MATNOW,MFHOUT,FMTHOL,(FIELDX(M,1),M=1,11),
     1 (NPTUSE(I),I=2,5),NPTUSE(1)
      WRITE(*   ,640) ZABCD,MATNOW,MFHOUT,FMTHOL,(FIELDX(M,1),M=1,11),
     1 (NPTUSE(I),I=2,5),NPTUSE(1)
C-----PRINT WARNING MESSAGE IF SELF-SHIELDING CALCULATION WAS REQUESTED
C-----AND THE TOTAL CROSS SECTION IS NOT PRESENT.
      IF(MFHOUT.NE.3.OR.NPTUSE(2).GT.0.OR.NOSELF.NE.0) GO TO 340
      WRITE(OUTP,620)
      WRITE(*   ,620)
C-----PRINT WARNING IF MINIMUM IS LESS THAN ALLOWABLE VALUE FOR ANY
C-----CROSS SECTION.
  340 DO 400 I=1,NSECT
      IF(REDMIN(I).GE.OKMIN(I)) GO TO 400
      CALL OUT9G(REDMIN(I),FIELDX(1,1))
      IF(I.lt.2) go to 350
      IF(I.eq.2) go to 360
      IF(I.lt.4) go to 370
      IF(I.eq.4) go to 380
      go to 390
C-----OTHER.
  350 WRITE(OUTP,570) (FIELDX(M,1),M=1,11)
      GO TO 400
C-----TOTAL.
  360 WRITE(OUTP,580) (FIELDX(M,1),M=1,11)
      GO TO 400
C-----ELASTIC.
  370 WRITE(OUTP,590) (FIELDX(M,1),M=1,11)
      GO TO 400
C-----CAPTURE.
  380 WRITE(OUTP,600) (FIELDX(M,1),M=1,11)
      GO TO 400
C-----FISSION.
  390 WRITE(OUTP,610) (FIELDX(M,1),M=1,11)
  400 CONTINUE
C-----INCREMENT POINT COUNTS FOR ENTIRE FILE.
      DO 410 I=1,5
  410 NPTTOT(I)=NPTTOT(I)+NPTUSE(I)
C-----BRANCH BACK TO START NEXT MATERIAL.
      GO TO 20
C=======================================================================
C
C     END OF RUN. IF MULTIBAND PARAMETERS WERE GENERATED WRITE REPORT
C     FOR EACH EVALUATION INDICATING THE MAXIMUM ERROR IN ANY GROUP
C     THAT WILL RESULT IF MULTIBAND PARAMETERS ARE USED TO DEFINE
C     SELF-SHIELDED CROSS SECTIONS AT ANY POINT ALONG THE ENTIRE
C     SELF-SHIELDED CURVE FOR SIGMA0=0 TO INFINITY. ALSO WRITE SUMMARY
C     OF MAXIMUM ERROR FOR ANY NUMBER OF BANDS, FOR EACH COMBINATION
C     OF (SIGMAO,N).
C
C=======================================================================
C-----PRINT WARNING MESSAGE IF NO DATA SATISFIED REQUESTS.
  420 IF(LASTZA.GT.0) GO TO 430
      WRITE(OUTP,750)
      WRITE(*   ,750)
      GO TO 530
C-----OUTPUT TOTAL POINT COUNTS.
  430 WRITE(OUTP,660) (NPTTOT(I),I=2,5),NPTTOT(1)
      WRITE(*   ,660) (NPTTOT(I),I=2,5),NPTTOT(1)
C-----NO MULTIBAND LIBRARY SUMMARY IF MULTIBAND PARAMETERS WERE
C-----NOT GENERATED.
      IF(LZA.LE.0) GO TO 530
C-----FOR EACH EVALUATION LIST ERROR VS. BANDS TABLE.
      WRITE(OUTP,690)
      DO 490 LLZA=1,LZA
      NBNEED=NBNTAB(LLZA)
      GO TO (440,450,460,470,480),NBNEED
  440 WRITE(OUTP,740) IZATAB(LLZA),ERBTAB(1,LLZA)
      GO TO 490
  450 WRITE(OUTP,700) IZATAB(LLZA),(ERBTAB(I,LLZA),I=1,NBNEED)
      GO TO 490
  460 WRITE(OUTP,710) IZATAB(LLZA),(ERBTAB(I,LLZA),I=1,NBNEED)
      GO TO 490
  470 WRITE(OUTP,720) IZATAB(LLZA),(ERBTAB(I,LLZA),I=1,NBNEED)
      GO TO 490
  480 WRITE(OUTP,730) IZATAB(LLZA),(ERBTAB(I,LLZA),I=1,NBNEED)
  490 CONTINUE
C-----LIST MAXIMUM ERROR VS. SIGMA0 FOR THE ENTIRE FILE.
      WRITE(OUTP,670)
      DO 510 I=1,25
      DO 500 NB=1,NBMAX
  500 ERLIB(I,NB)=100.0*ERLIB(I,NB)
  510 WRITE(OUTP,680) I,(ERLIB(I,J),J=1,NBMAX)
C-----------------------------------------------------------------------
C
C     PRINT SUMMARY OF TYPES OF CONVERGENCE FOR MULTI-BAND PARAMETERS
C
C-----------------------------------------------------------------------
      MYLOGSUM =0
      do k=1,8
      MYLOGSUM = MYLOGSUM + MYLOG(k)
      enddo
      WRITE(*,520) MYLOG,MYLOGSUM,LZA*NGR
      WRITE(3,520) MYLOG,MYLOGSUM,LZA*NGR
  520 FORMAT(1X,78('-')/' Summary of Multi-Band Results'/1X,78('-')/
     1 ' No Self-Shielding (Conserve 1 Moment)..........',I11/
     2 ' Little Self-Shielding (Conserve 2 Moments).....',I11/
     3 ' General Self-Shielding (Conserve 3 Moments)....',I11/
     4 ' Strict Convergence.............................',I11/
     5 ' Not Strict Convergence.........................',I11/
     6 ' Soft Convergence...............................',I11/
     7 ' Very Soft Convergence..........................',I11/
     8 ' No Convergence.................................',I11/
     9 ' Sum............................................',I11/
     9 ' Check (Evaluations x Groups)...................',I11)
C-----IF ENDF/B FORMAT OUTPUT, OUTPUT FEND. MEND AND TEND LINES, AS
C-----REQUIRED.
  530 MATH=-1
      MFH=0
      MTH=0
      NOSEQ=0
      CALL CONTOG
C-----END ENDF/B FORMATTED FILE AND LISTING FILES.
      WRITE(OUTP,630)
      WRITE(*   ,630)
      OTAPE = 11        ! FORCE ON-LINE RUNING TIME REPORT.
CAK
      CLOSE (OTAPE)
      RETURN
CAK
      CALL ENDIT
      GO TO 420    ! cannot get to here.
  540 FORMAT(///' MAT/MF/MT/N1=',I5,I3,I4,I5,
     1 ' Interpolation Law is Not Llnear-Linear'/
     2 '   NBT   INT'/1X,11('-'))
  550 FORMAT(2I6)
  560 FORMAT(//'  **** Execution Terminated ***')
  570 FORMAT(' WARNING....Minimum Other   in Above',
     1 ' Material is',
     2 1X,11A1,' barns.')
  580 FORMAT(' WARNING....Minimum Total   in Above',
     1 ' Material is',
     2 1X,11A1,' barns.')
  590 FORMAT(' WARNING....Minimum Elastic in Above',
     1 ' Material is',
     2 1X,11A1,' barns.')
  600 FORMAT(' WARNING....Minimum Capture in Above',
     1 ' Material is',
     2 1X,11A1,' barns.')
  610 FORMAT(' WARNING....Minimum Fission in Above',
     1 ' Material is',
     2 1X,11A1,' barns.')
  620 FORMAT(' WARNING...Total Cross Section is Not Given for Above',
     1 ' Material.'/
     2 ' Self-Shielding and Multiband Calculation Skipped',
     3 ' (i.e. Impossible)')
  630 FORMAT(1X,78('-')/' Other Points = Points Other Than Total,',
     1 ' Elastic, Capture or Fission.'/
     3 ' This Will Always be Zero Unless Unshielded Output is',
     4 ' Requested.'/1X,78('-'))
  640 FORMAT(10A1,I5,I3,2X,A1,11A1,5I8)
  650 FORMAT(1X,78('-')/' ENDF/B Tape Label'/1X,78('-')/1X,16A4,A2,I4/
     1 1X,78('-')/' Isotope    MAT MF Fmt Kelvin   ',
     2 '   Total Elastic Capture Fission   Other'/13X,'     ',14X,
     3 '  Points  Points  Points  Points  Points'/1X,78('-'))
  660 FORMAT(1X,78('-')/25X,'Totals ',5I8)
  670 FORMAT(///8X,'Maximum per-cent Error vs. Sigma-0 for Entire',
     1 ' FILE'/1X,78('-')/' Index *    1 Band  *  2 Bands  ',
     2 '*  3 Bands  *  4 Bands  *  5 Bands'/1X,78('-'))
  680 FORMAT(I6,1X,'*',F10.2,2X,4('*',F9.2,2X))
  690 FORMAT(///8X,'Maximum per-cent Error During Sigma-0',
     1 ' INTERPOLATION'/1X,78('-')/'    ZA *    1 Band  *  2 Bands  ',
     2 '*  3 Bands  *  4 Bands  *  5 Bands'/1X,78('-'))
  700 FORMAT(I6,1X,'*',F10.2,2X,'*',F9.2,2X,3('*',11X))
  710 FORMAT(I6,1X,'*',F10.2,2X,2('*',F9.2,2X),2('*',11X))
  720 FORMAT(I6,1X,'*',F10.2,2X,3('*',F9.2,2X),'*')
  730 FORMAT(I6,1X,'*',F10.2,2X,4('*',F9.2,2X))
  740 FORMAT(I6,1X,'*',F10.2,2X,4('*',11X))
  750 FORMAT(/' WARNING....No Data Found That Satisfied Retrieval',
     1 ' Criteria.'/11X,
     2 ' Therefore No Data was Group Averaged or Written to Output',
     3 'FILE.'/1X,78('-'))
      END
CAK   SUBROUTINE FILE1
      SUBROUTINE FILE1g
C=======================================================================
C
C     ADD COMMENTS AT THE END OF FILE 1, SECTION 451 TO INDICATE
C     THAT THIS MATERIAL HAS BEEN PROCESSED BY PROGRAM GROUPIE AND
C     TO SPECIFY THE NUMBER OF GROUPS USED.
C
C     DEFINE FORMAT TO BE ENDF/B-4, 5 OR 6.
C
C     THE ENDF/B FORMAT CAN BE DETERMINED FROM THE SECOND CARD.
C     ENDF/B-4  = N1 > 0, N2 = 0, CARD COUNT (POSITIVE)
C     ENDF/B-5  = N1 = N2 = 0
C     ENDF/B-6         N2 = VERSION NUMBER (6 OR MORE)
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 PROGDOC1
      CHARACTER*4 FMTTAB,FMTHOL
      CHARACTER*66 PROGDOC
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/TEMP,Q,L1,L2,N1,N2,MAT,MF,MT
CAK   COMMON/HOLFMT/FMTHOL
      COMMON/HOLFMTg/FMTHOL
      COMMON/SPIM/IMSP
      COMMON/VERSES/TEMP3,IVERSE
      COMMON/RESCOM/RES1,RES2,URES1,URES2,LSSF,ICOMP,INPART
      INCLUDE 'groupie.h'
      DIMENSION FMTTAB(4),PROGDOC(6),PROGDOC1(66,6)
      EQUIVALENCE (PROGDOC(1),PROGDOC1(1,1))
      DATA FMTTAB/'4','5','6','7'/
C-----DOCUMENTATION TO ADD TO ENDF/B OUTPUT - EACH LINE IS 66
C-----CHARACTERS LONG - FIELDS 12345678901 ARE FILLED IN WITH
C-----11 CHARACTERS DURING EXECUTION.
C                1         2         3         4         5         6
C       12345678901234567890123456789012345678901234567890123456789012
C       3456
      DATA PROGDOC/
     1 ' ***************** Program GROUPIE (VERSION 2017-1)***********',
     2 ' Unshielded Group Averages Using 12345 Groups                 ',
     3 ' Weighting Spectrum: Maxwellian, 1/E and Fission Spectrum     ',
     4 ' Weighting Spectrum: 1/E Spectrum                             ',
     5 ' Weighting Spectrum: Flat (Constant) Spectrum                 ',
     6 ' Weighting Spectrum: Input Spectrum                           '/
C-----WARNING - Only 3 of the above 6 lines are OUTPUT.
C
C-----FILL IN REMAINDER OF FIRST LINE.
      PROGDOC1(63,1) = '*'
      PROGDOC1(64,1) = '*'
      PROGDOC1(65,1) = '*'
      PROGDOC1(66,1) = '*'
      INPART= 1     ! Initialize to neutron (in case not ENDF6 format)
C-----HEAD CARD OF SECTION HAS BEEN READ. WRITE IT AND READ NEXT CARD
C-----AND DETERMINE IF THIS IS THE ENDF/B-IV, V OR VI FORMAT.
      CALL CONTOG
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      IVERSE=4
C-----CHECK FOR ENDF/B-IV.
C-----IV N1 > 0, N2 = 0
      IF(N1.GT.0.AND.N2.EQ.0) GO TO 10
C-----NOT ENDF/B-IV. READ THIRD CARD.
      N2IN=N2
      CALL CARDO(C1,C2,L1,L2,N1,N2)
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      IVERSE=5
C-----CHECK FOR ENDF/B-V FORMAT.
      IF(N2IN.LE.0) GO TO 10
      N1IN = N1
C-----ENDF/B-VI FORMAT. READ FOURTH CARD.
      CALL CARDO(C1,C2,L1,L2,N1,N2)
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      IVERSE=6
      TEMP3=C1
      INPART = N1IN/10
C-----SET DERIVED MATERIAL FLAG.
      L1=1
C-----DEFINE ENDF/B FORMAT NUMBER.
   10 FMTHOL=FMTTAB(IVERSE-3)
C-----SKIP OUTPUT IF NO ENDF/B FORMATTED OUTPUT.
      IF(OTAPE.LE.0) GO TO 30
C-----INCREASE COMMENT CARD COUNT AND COPY TO END OF HOLLERITH.
      N1P3=N1+3
      CALL CARDO(C1,C2,L1,L2,N1P3,N2)
      DO 20 N=1,N1
   20 CALL COPY1G
C
C     ADD THREE COMMENT LINES.
C
C-----PROGRAM NAME AND VERSION
      CALL HOLLYO(PROGDOC1(1,1))
C-----NUMBER OF GROUPS
      CALL INTOUT(NGR,PROGDOC1(34,2),5)
      CALL HOLLYO(PROGDOC1(1,2))
C-----WEIGHTING SPECTRUM.
      IF(IMSP.EQ.1) CALL HOLLYO(PROGDOC1(1,3))
      IF(IMSP.EQ.2) CALL HOLLYO(PROGDOC1(1,4))
      IF(IMSP.EQ.3) CALL HOLLYO(PROGDOC1(1,5))
      IF(IMSP.EQ.4) CALL HOLLYO(PROGDOC1(1,6))
C-----COPY TO END OF FILE (MF=0)
   30 CALL COPYFG    ! Copy to end of FILE (MF=0)
      RETURN
      END
      SUBROUTINE FILE2
C=======================================================================
C
C     SIGMA1 FILE2 ADAPTED FOR USE BY GROUPIE.
C     ========================================
C     READ RESONANCE PARAMETERS IN ORDER TO DEFINE THE ENERGY RANGE
C     OF THE RESOLVED AND UNRESOLVED RESONANCE REGION.
C
C     NO TESTS FOR INCONSISTENCY = STOP ON ERROR
C
C     RES1    = Lower Energy limit of Resolved.
C     RES2    = Upper Energy limit of Resolved.
C     URES1   = Lower Energy limit of Unresolved.
C     URES2   = Upper Energy limit of Unresolved.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*4 FMTHOL
      INTEGER*4 OUTP,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LEADER/C1,C2,L1,L2,N1,N2,MAT,MF,MT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/RESCOM/RES1,RES2,URES1,URES2,LSSF,ICOMP,INPART
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
CAK   COMMON/HOLFMT/FMTHOL
      COMMON/HOLFMTg/FMTHOL
C-----Initialize
      INPART= 1     ! Initialize to neutron (in case not ENDF6 format)
      RES1  = 0.0d0
      RES2  = 0.0d0
      URES1 = 0.0d0
      URES2 = 0.0d0
c-----GROUPIE/SIGMA1 Differ
      CALL CONTOG
C-----HEAD RECORD ALREADY READ. DEFINE NUMBER OF ISOTOPES.
      NIS=N1H
C-----DO LOOP OVER ALL ISOTOPES
      DO 240 IS=1,NIS
      CALL CARDIO(C1H,C2H,L1H,LFW,NER,N2H)
C-----DO LOOP OVER ALL ENERGY RANGES
      DO 230 JER=1,NER
      CALL CARDIO(EL,EH,LRU,LRF,N1H,N2H)
c
c     2017/5/16 - added NRO - Energy dependent scattering radius
c
      NRO = N1H
      if(NRO.eq.1) then
c-----Energy dependent scattering radius = copy TAB1 record
      CALL CARDIO(C1,C2,L1,L2,N1,N2)
      do i=1,N1,3                  ! Interpolation law (NBT,INT) Pairs
      CALL COPY1G
      enddo
      do I=1,N2,3                  ! Scattering radius (X,Y) Pairs
      CALL COPY1G
      enddo
      endif
c
C-----DEFINE LRU FOR INTERNAL USE AS ORIGINAL LRU (BEFORE RECENT).
      IF(LRU.GT.3) LRU=LRU-3
C-----Select resolved or unresolved or NONE.
      IF(LRU.eq.1) go to 10    ! Resolved
      IF(LRU.gt.1) go to 20    ! unresolved
C
C     NO RESONANCE PARAMETERS PRESENT
C
C-----COPY SECTION WITH NO RESONANCE PARAMETERS.
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,N2H)
      GO TO 230
C-----------------------------------------------------------------------
C
C     Resolved.
C
C-----------------------------------------------------------------------
C-----DEFINE UNITED ENERGY RANGE.
   10 IF(RES2.le.0.0d0) then
      RES1 = EL
      RES2 = EH
      else
      IF(EL.lt.RES1) RES1 = EL
      IF(EH.gt.RES2) RES2 = EH
      endif
      go to 30
C-----------------------------------------------------------------------
C
C     Unresolved.
C
C-----------------------------------------------------------------------
C-----DEFINE UNITED ENERGY RANGE.
   20 IF(URES2.le.0.0d0) then
      URES1 = EL
      URES2 = EH
      else
      IF(EL.lt.URES1) URES1 = EL
      IF(EH.gt.URES2) URES2 = EH
      endif
C-----------------------------------------------------------------------
C
C     RESONANCE PARAMETERS PRESENT
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C     LRU= 1 - RESOLVED
C     LRF= 1 - SLBW, = 2 - MLBW, = 3 - REICH-MOORE, = 4 - ADLER-ADLER
C            - NEW REICH-MOORE = 7
C
C     LRU= 2 - UNRESOLVED
C     LRF= 1 - ENERGY INDEPENDENT WIDTHS (EXCEPT POSSIBLY FISSION)
C        = 2 - ENERGY   DEPENDENT WIDTHS
C
C-----------------------------------------------------------------------
   30 IF(LRU.NE.1) GO TO 40    ! Resolved?
      IF(LRF.EQ.1.OR.          ! Single Level Breit-Wigner
     1   LRF.EQ.2.OR.          ! Multi-Level  Breit-Wigner
     2   LRF.EQ.3) GO TO 50    ! Reich=Moore
      IF(LRF.EQ.4) GO TO 100   ! Adler-Adler
      IF(LRF.EQ.7) GO TO 130   ! New Reich-Moore
C-----ILLEGAL - IGNORE REMAINDER OF FILE 2
      GO TO 250
   40 IF(LRU.NE.2) GO TO 250   ! Unresolved?
      IF(LRF.EQ.1) GO TO 150   ! Energy Independent Widths
      IF(LRF.EQ.2) GO TO 210   ! Energy   Dependent Widths
C-----ILLEGAL - IGNORE REMAINDER OF FILE 2
      GO TO 250
C-----------------------------------------------------------------------
C
C     BREIT-WIGNER (SINGLE OR MULTI-LEVEL) OR REICH-MOORE FORMALISM
C
C-----------------------------------------------------------------------
C-----IF SCATTERING RADIUS IS ENERGY DEPENDENT COPY TABULATED VALUES.
   50 IF(N1H.EQ.0) GO TO 80
      IF(LRF.NE.1.AND.LRF.NE.2) GO TO 80
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,N2H)
C-----SKIP SCATTERING RADIUS INTERPOLATION LAW.
      DO 60 I=1,N1H,3
   60 CALL COPY1G
C-----SKIP SCATTERING RADIUS TABULATED VALUES.
      DO 70 I=1,N2H,3
   70 CALL COPY1G
C-----READ NEXT CARD.
   80 CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
C-----LOOP OVER ALL L STATES
      DO 90 ILS=1,NLS
C-----READ NEXT CARD.
      CALL CARDIO(C1H,C2H,L1H,L2H,NRS6,NRS)
C-----COPY RESONANCE PARAMETERS.
      DO 90 IRS=1,NRS
   90 CALL COPY1G
      GO TO 230
C-----------------------------------------------------------------------
C
C     ADLER-ADLER FORMALISM
C
C-----------------------------------------------------------------------
C-----READ NEXT CARD.
  100 CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
C-----READ BACKGROUND CORRECTIONS.
      CALL CARDIO(C1H,C2H,L1H,L2H,NX6,N2H)
C-----COPY BACKGROUND CORRECTION CONSTANTS.
      DO 110 I=1,NX6,6
  110 CALL COPY1G
C-----LOOP OVER L STATES
      DO 120 I=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,NJS,N2H)
C-----LOOP OVER J STATES
      DO 120 J=1,NJS
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,NLJ)
C-----COPY ALL RESONANCE DATA
      DO 120 K=1,NLJ
      CALL COPY1G
  120 CALL COPY1G
      GO TO 230
C-----------------------------------------------------------------------
C
C     NEW REICH-MOORE FORMALISM
C
C-----------------------------------------------------------------------
C-----DEFINE NUMBER OF J STATES
  130 CALL CARDIO(C1,C2,L1,L2,NJS,N2)
C-----DEFINE NUMBER OF PARTICLE-PAIRS
      CALL CARDIO(C1,C2,L1,L2,NPP12,N2)
C-----COPY PARTICLE-PAIR DATA
      DO N=1,NPP12,6
      CALL COPY1G
      ENDDO
C-----LOOP OVER J STATES
      DO 140 IJ=1,NJS
C-----J, PARITY, AND NUMBER OF CHANNELS
      CALL CARDIO(C1,C2,L1,L2,NCH6,N2)
C-----COPY CHANNEL DATA
      DO N=1,NCH6,6
      CALL COPY1G
      ENDDO
C-----DEFINE NUMBER OF RESONANCES
      CALL CARDIO(C1,C2,L1,L2,N1,NRS)
C-----COPY RESONANCE PARAMETERS
      DO N=1,NRS
      CALL COPY1G
      ENDDO
  140 CONTINUE
      GO TO 230
C-----------------------------------------------------------------------
C
C     UNRESOLVED WITH ENERGY INDEPENDENT WIDTHS (LRF = 1)
C
C-----------------------------------------------------------------------
C-----TEST IF FISSION WIDTHS GIVEN
  150 IF(LFW.gt.0) go to 170
C
C     Case A: FISSION WIDTHS NOT GIVEN (LFW = 0/ LRF = 1)
C
      CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
      LSSF = L1H
C-----LOOP OVER ALL L-STATES
      DO 160 ILS=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,N1H,NJS)
      DO 160 N=1,NJS
  160 CALL COPY1G
      GO TO 230
C
C     Case B: FISSION WIDTHS GIVEN (LFW = 1/ LRF = 1)
C
  170 CALL CARDIO(C1H,C2H,L1H,L2H,NE,NLS)
      LSSF = L1H
C-----COPY FISSION WIDTH ENERGY POINTS
      DO 180 I=1,NE,6
  180 CALL COPY1G
C-----LOOP OVER L-STATES
      DO 200 I=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,NJS,N2H)
C-----LOOP OVER J STATES
      DO 200 J=1,NJS
      CALL CARDIO(C1H,C2H,L1H,L2H,NEP6,N2H)
      DO 190 K=1,NEP6,6
  190 CALL COPY1G
  200 CONTINUE
      GO TO 230
C-----------------------------------------------------------------------
C
C     Case C: UNRESOLVED WITH ALL ENERGY DEPENDENT WIDTHS (LRF = 2)
C             Independent of LFW
C
C-----------------------------------------------------------------------
C-----READ NEXT CARD.
  210 CALL CARDIO(C1H,C2H,L1H,L2H,NLS,N2H)
      LSSF = L1H
C-----DO LOOP OVER L-STATES
      DO 220 I=1,NLS
      CALL CARDIO(C1H,C2H,L1H,L2H,NJS,N2H)
      DO 220 J=1,NJS
      CALL CARDIO(C1H,C2H,L1H,L2H,NE6P6,N2H)
C-----COPY NUMBER OF DEGREES OF FREEDOM AND PARAMETERS.
      DO 220 K=1,NE6P6,6
  220 CALL COPY1G
C-----END OF ENERGY RANGE LOOP
  230 CONTINUE
C-----END OF ISOTOPE LOOP
  240 CONTINUE
C-----------------------------------------------------------------------
C
C     FINISHED O.K.
C
C-----------------------------------------------------------------------
C-----COPY TO END OF FILE (MF=0)
      CALL COPYFG
      RETURN
C-----------------------------------------------------------------------
C
C     ERROR.
C
C-----------------------------------------------------------------------
  250 WRITE(3,260)
      WRITE(*,260)
  260 FORMAT(///' ERROR - Reading MF/MT=2/151 Resonance Parameters.'/
     1          '         Execution Terminated.'///)
      CALL ENDERROR
      RETURN          ! Dummy RETURN to satisfy some compilers.
      END
      SUBROUTINE PAGIN5(NTAPE,ITYPE,NTYPE)
C=======================================================================
C
C     READ TABLE OF DATA IN THE ENDF/B FORMAT AND LOAD IT INTO THE
C     PAGING SYSTEM. THIS ROUTINE IS USED TO READ EITHER THE ENERGY
C     DEPENDENT WEIGHTING SPECTRUM FROM THE INPUT FILE OR A SECTION
C     OF CROSS SECTIONS FROM THE ENDF/B FORMAT FILE. IN EITHER CASE
C     THIS ROUTINE WILL ONLY READ (X,Y) PAIRS IN 6E11.4 FORMAT (I.E.
C     IF READING FROM THE ENDF/B FORMAT THE SECTION HEAD CARDS AND
C     INTERPOLATION LAW MUST BE READ BEFORE CALLING THIS ROUTINE).
C
C     ARGUMENTS
C     ---------
C     NTAPE = LOGICAL NUMBER UNIT OF INPUT FILE.
C     ITYPE = PAGING SYSTEM STORAGE INDEX THAT DEFINES WHERE TO STORE
C             DATA.
C           = 1 ENERGY DEPENDENT WEIGHTING SPECTRUM
C           = 2 TOTAL CROSS SECTION (IF SELF-SHIELDING CALCULATION).
C               ALL CROSS SECTIONS (IF ONLY UNSHIELDED CALCULATION).
C           = 3 ELASTIC CROSS SECTION
C           = 4 CAPTURE CROSS SECTION.
C               ALL CROSS SECTIONS EXCEPT TOTAL, ELASTIC OR FISSION
C               (IF SELF-SHIELDING CALCULATION).
C           = 5 FISSION CROSS SECTION
C     NTYPE = DEFINES THE TYPE OF DATA
C           = 1 ENERGY DEPENDENT WEIGHTING SPECTRUM.
C               ALL CROSS SECTIONS (EXCEPT TOTAL, ELASTIC, CAPTURE
C               OR FISSION).
C           = 2 TOTAL CROSS SECTION
C           = 3 ELASTIC CROSS SECTION
C           = 4 CAPTURE CROSS SECTION
C           = 5 FISSION CROSS SECTION
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELDX
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/ZA,AWRIN,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
CAK   COMMON/PAGER/NPAGE,NPAGM1
      COMMON/PAGERg/NPAGE,NPAGM1
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
      COMMON/LASTE/ELAST
      COMMON/MINOK/OKMIN(5),REDMIN(5)
      COMMON/TABASE/NPT
      COMMON/RESCOM/RES1,RES2,URES1,URES2,LSSF,ICOMP,INPART
      COMMON/INPAGE/IXYLOW(5),IXYHI(5),ISCR(5)
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      INCLUDE 'groupie.h'
C-----PRINT TITLE FOR WEIGHTING SPECTRUM.
      IF(ITYPE.EQ.1) WRITE(OUTP,120)
C-----DEFINE SCRATCH UNIT.
      NSCR=ISCR(ITYPE)
C-----READ FROM STANDARD FILE NAME.
      ISAVE=ITAPE
      ITAPE=NTAPE
C-----INITIALIZE LAST ENERGY FOR ASCENDING ENERGY TEST.
      ELAST=0.0
C-----SET UP LOOP OVER PAGES.
      N2X=NPTAB(ITYPE)
      DO 70 NPT=1,N2X,NPAGE
C-----READ NEXT PAGE.
      NNPT=NPT+NPAGM1
      IF(NNPT.GT.N2X) NNPT=N2X
      IHIGH=NNPT-NPT+1
      CALL POINTI(XPAGE(1,ITYPE),YPAGE(1,ITYPE),IHIGH)
c-----------------------------------------------------------------------
c
c     Define ICOMP for URR Treatment
c     1) Only beginning of table
c     2) Only MF=3
c     3) Ony threshold below URR Max. energy
c     4) MT =  51 to  91 = Competition
c     5) MT = 103 to 107 = Absorption
c
c-----------------------------------------------------------------------
      IF(NPT.ne.1.or.MFH.ne.3) go to 10
      if(XPAGE(1,ITYPE).ge.URES2) go to 10
      if(MTH.ge.51.and.MTH.le.91) then    ! Competition
      ICOMP = MTH
      if(MTH.gt.51) ICOMP = 4             ! More than 1 level = use MT=4
      endif
      if(MTH.ge.103.and.MTH.le.107) then         ! Absorption
      if(ICOMP.lt.1000) ICOMP = 1000*MTH + ICOMP ! Only allow 1
      endif
C
C     CHECK FOR MINIMUM ALLOWABLE VALUE.
C
C-----(ALLOW MU-BAR, XI AND ETA TO BE NEGATIVE).
   10 IF(MTH.GE.251.AND.MTH.LE.253) GO TO 30
      DO 20 I=1,IHIGH
C-----IGNORE POINTS OUTSIDE MULTIGROUP ENERGY RANGE.
      IF(XPAGE(I,ITYPE).LT.EGROUP(    1)) GO TO 20
      IF(XPAGE(I,ITYPE).GT.EGROUP(NGRP1)) GO TO 20
      IF(YPAGE(I,ITYPE).LT.REDMIN(NTYPE)) REDMIN(NTYPE)=YPAGE(I,ITYPE)
   20 CONTINUE
C-----PRINT WEIGHTING SPECTRUM.
   30 IF(ITYPE.NE.1) GO TO 60
      DO 50 I=1,IHIGH,3
      II=I+2
      IF(II.GT.IHIGH) II=IHIGH
      J=0
      DO 40 III=I,II
      J=J+1
      CALL OUT9G(XPAGE(III,1),FIELDX(1,J))
      J=J+1
   40 CALL OUT9G(YPAGE(III,1),FIELDX(1,J))
   50 WRITE(OUTP,130) ((FIELDX(M,JJ),M=1,11),JJ=1,J)
C-----IF OVER ONE PAGE OF DATA MOVE DATA TO SCRATCH FILE.
   60 IF(N2X.LE.NPAGE) GO TO 70
      IF(NPT.EQ.1) REWIND NSCR
      CALL OBLOCK(NSCR,XPAGE(1,ITYPE),YPAGE(1,ITYPE),NPAGE)
   70 CONTINUE
C-----IS CROSS SECTION IN CORE OR ON SCRATCH.
      IF(N2X.GT.NPAGE) GO TO 80
C-----IN CORE. SET INDICES TO ENTIRE TABLE.
      IXYLOW(ITYPE)=0
      IXYHI(ITYPE)=N2X
      GO TO 90
C-----ON SCRATCH. LOAD FIRST PAGE AND SET INDICES TO FIRST PAGE.
   80 END FILE NSCR
      REWIND NSCR
      CALL IBLOCK(NSCR,XPAGE(1,ITYPE),YPAGE(1,ITYPE),NPAGE)
      IXYLOW(ITYPE)=0
      IXYHI(ITYPE)=NPAGE
C-----IF INPUT WEIGHTING SPECTRUM IS NEGATIVE PRINT MESSAGE AND
C-----TERMINATE.
   90 IF(ITYPE.NE.1.OR.(REDMIN(1).GE.OKMIN(1))) GO TO 100
      CALL OUT9G(REDMIN(1),FIELDX(1,1))
      WRITE(OUTP,110) (FIELDX(M,1),M=1,11)
      CALL ENDERROR
C-----RESTORE STANDARD UNIT NUMBER.
  100 ITAPE=ISAVE
      RETURN
  110 FORMAT(/' Minimum Weighting Spectrum Value is'/
     1 1X,11A1,'. it MUST be Positive for Calculations.'/
     2 ' Execution Terminated.')
  120 FORMAT(1X,78('-')/' Weighting Spectrum (Linearly Interpolable)'/
     1 1X,78('-')/'   Energy-eV   Spectrum',2X,
     2 '  Energy-eV   Spectrum',2X,'  Energy-eV   Spectrum'/
     3 1X,78('-'))
  130 FORMAT(1X,3(22A1,2X))
      END
      SUBROUTINE GROUPU
C=======================================================================
C
C     THIS ROUTINE IS DESIGNED TO COMPUTE UNSHIELDED GROUP
C     AVERAGED CROSS SECTIONS FOR ONE REACTION. BOTH THE ENERGY
C     DEPENDENT WEIGHTING SPECTRUM AND THE CROSS SECTION MUST BE
C     LINEAR INTERPOLABLE BETWEEN ADJACENT TABULATED POINTS. ALL
C     INTEGRALS ARE PERFORMED ANALYTICALLY.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,OTAPE2,TYPSIG
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
      COMMON/COMXTEND/MYXTEND
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
CAK   COMMON/FILLER/ISECT
      COMMON/FILLERg/ISECT
      COMMON/MINOK/OKMIN(5),REDMIN(5)
      COMMON/TYPER/NTYPE
      INCLUDE 'groupie.h'
c-----2016/3/27 - AVN increased to 20,001 from 1001
      DIMENSION JSECT(2),LPTAB(5),AVN(20001)
C-----INITIALIZE INDICES TO SOURCE SPECTRUM.
      DATA JSECT/1,0/
      DATA LPTAB/5*1/
C
C     FIND FIRST CROSS SECTION POINT OF INTEREST (EITHER THE FIRST
C     TABULATED POINT OR THE THRESHOLD. THE THRESHOLD IS HEREIN
C     DEFINED TO BE THE HIGHEST ENERGY POINT AT WHICH THE CROSS
C     SECTION IS ZERO. IN ENDF/B THIS IS NOT NECESSARILY THE FIRST
C     TABULATED DATA POINT). DEFINE ENERGY AT THIS POINT. SAVE
C     INDEX TO FIRST POINT OF INTEREST (LPTAB) AND INITIALIZE POINT
C     INDEX TO PRECEEDING POINT (IPTAB).
C
      JJ=NPTAB(ISECT)
      DO 10 JP=1,JJ
      J=JP
      CALL XYPAGE(J,ISECT,XA(2),YA(2))
      IF(YA(2).gt.0.0D+0) go to 20
   10 XALAST=XA(2)
      GO TO 180
   20 IPTAB(1)=0
      J=J-1
      IF(J.GT.0) XA(2)=XALAST
      IF(J.LT.1) J=1
      LPTAB(ISECT)=J
      IPTAB(ISECT)=J-1
C
C     OUTPUT SHOULD NOT INCLUDE GROUPS WHOLLY OUTSIDE THE ENERGY
C     RANGE OF THE REACTION. SKIP ALL GROUPS THAT PRECEED FIRST
C     CROSS SECTION POINT OF INTEREST (IF ANY). SAVE INDEX TO FIRST
C     GROUP OF INTEREST (IGRLOW) TO DEFINE WHERE OUTPUT BEGINS.
C
      DO 30 IGR=2,NGRP1
      IF(XA(2).LT.EGROUP(IGR)) GO TO 40
   30 CONTINUE
      GO TO 180
   40 IGRLOW=IGR-1
      JSECT(2)=ISECT
C
C     DEFINE CROSS SECTION AND SPECTRUM VALUES OVER FIRST INTERVAL
C     OF FIRST GROUP OF INTEREST (FROM LOWER GROUP BOUNDARY TO FIRST
C     TABULATED POINT ABOVE LOWER GROUP BOUNDARY). IF THRESHOLD IS
C     WITHIN GROUP DEFINE CROSS SECTION BETWEEN LOWER GROUP BOUNDARY
C     AND THRESHOLD TO BE ZERO.
C
      XB=EGROUP(IGRLOW)
      DO 90 I=1,2
      IISECT=JSECT(I)
C-----FIND FIRST DATA POINT ABOVE LOWER ENERGY BOUNDARY OF GROUP.
   50 IF(IPTAB(IISECT).LT.NPTAB(IISECT)) GO TO 60
      XA(I)=2.0D+20
      GO TO 90
   60 IPTAB(IISECT)=IPTAB(IISECT)+1
      CALL XYPAGE(IPTAB(IISECT),IISECT,XA(I),YA(I))
      IF(XA(I).gt.XB) go to 70
C-----POINT IS STILL BELOW ENERGY LIMIT. SAVE ENERGY AND CROSS SECTION
C-----FOR INTERPOLATION.
      XC=XA(I)
      YB(I)=YA(I)
      GO TO 50
C-----POINT IS ABOVE LOWER ENERGY LIMIT. IF THIS IS FIRST POINT OF
C-----INTEREST DEFINE PRECEDING INTERVAL WITH ZERO CROSS SECTION.
C-----OTHERWISE DEFINE CROSS SECTION AT LOWER ENERGY BOUNDARY OF
C-----GROUP BY LINEAR INTERPOLATION BETWEEN LAST TWO POINTS.
   70 IF(IPTAB(IISECT).GT.LPTAB(IISECT)) GO TO 80
      YA(I)=0.0
      YB(I)=0.0
      IPTAB(IISECT)=LPTAB(IISECT)-1
      GO TO 90
   80 YB(I)=((XB-XC)*YA(I)+(XA(I)-XB)*YB(I))/(XA(I)-XC)
   90 CONTINUE
C
C     GROUP LOOP. DEFINE UPPER ENERGY LIMIT OF GROUP. INITIALIZE
C     INTEGRAL OF CROSS SECTION TIMES SPECTRUM (XCINT1) AND INTEGRAL
C     OF SPECTRUM (XCNRM1).
C
      DO 170 IGR=IGRLOW,NGR
      XE=EGROUP(IGR+1)
      XCINT1=0.0
      XCNRM1=0.0
C
C     DEFINE BEGINNING OF NEXT INTERVAL TO BE THE SAME AS THE END
C     OF THE LAST INTERVAL. SELECT LONGEST INTERVAL (XC TO XB)
C     WITHIN CURRENT ENERGY GROUP OVER WHICH BOTH SPECTRUM AND
C     CROSS SECTION ARE LINEARLY INTERPOLABLE (THIS IS MINIMUM ENERGY
C     OF XE, XA(1) AND XA(2)).
C
  100 XC=XB
      YC(1)=YB(1)
      YC(2)=YB(2)
      XB=XE
      IF(XA(1).LT.XB) XB=XA(1)
      IF(XA(2).LT.XB) XB=XA(2)
C
C     NEXT CALCULATION WILL BE INTEGRAL BETWEEN ENERGIES XC AND XB.
C     INTERPOLATE OVER THE ENERGY INTERVAL XC TO XA TO DEFINE THE
C     CROSS SECTION AND SPECTRUM AT ENERGY XB. IF THE ENTIRE INTERVAL
C     HAS THEN BEEN USED (I.E. XA.LE.XB) ADVANCE TO THE NEXT ENERGY
C     INTERVAL. SKIP ENERGY INTERVALS WITH NEGATIVE TOTAL CROSS SECTION
C     AND SKIP ALL EMPTY INTERVALS (ZERO LENGTH INTERVALS DUE TO
C     DISCONTINUITIES IN THE CROSS SECTION AND/OR SPECTRUM).
C
      DO 130 I=1,2
      IISECT=JSECT(I)
      IF(XA(I).gt.XB) go to 120
C-----NEXT INTEGRAL WILL EXTEND TO END OF INTERVAL XC TO XA. DEFINE
C-----INTERPOLATED VALUE AS END OF INTERVAL AND THEN SELECT NEXT
C-----INTERVAL. IF BEYOND THE TABULATED ENERGY RANGE EXTEND AS CONSTANT
C-----(I.E. KEEP SAME YA AND DEFINE NEW ENERGY AT 2.0D+20).
      YB(I)=YA(I)
      IF(IPTAB(IISECT).LT.NPTAB(IISECT)) GO TO 110
c-----2015/8/6 - Insert another point at end with either
c-----           same or 0 cross section.
      IF(IPTAB(IISECT).EQ.NPTAB(IISECT)) then
      IPTAB(IISECT)=IPTAB(IISECT)+1
      if(MYXTEND.LE.0) YA(I) = 0.0d0
      go to 130
      ENDIF
C-----Beyond tabulated values.
      XA(I)=2.0D+20           ! DEFINE ENERGY WELL ABOVE GROUPS
      IF(MYXTEND.LE.0) THEN   ! USE LAST TABULATED VALUE OR ZERO?
      YA(I) = 0.0d0
      YB(I) = 0.0d0
      GO TO 130               ! Cross section = 0 - skip calculation
      ENDIF
      GO TO 120               ! Use cross section extension
c-----Next tabulated point.
  110 IPTAB(IISECT)=IPTAB(IISECT)+1
      CALL XYPAGE(IPTAB(IISECT),IISECT,XA(I),YA(I))
      GO TO 130
C-----NEXT INTEGRAL WILL ONLY BE OVER A PORTION OF THE INTERVAL XC TO
C-----XA. INTERPOLATE BETWEEN XC AND XA TO DEFINE VALUE AT XB.
  120 YB(I)=((XB-XC)*YA(I)+(XA(I)-XB)*YC(I))/(XA(I)-XC)
  130 CONTINUE
      DE=XB-XC
      IF(DE.LE.0.0) GO TO 140
C-----------------------------------------------------------------------
C
C     04/09/12 - USE ENTIRE ENERGY RANGE EVEN IF TOTAL IS NOT > 0
C
C-----------------------------------------------------------------------
C
C     COMPUTE CONTRIBUTION TO INTEGRAL OF CROSS SECTION TIMES
C     SPECTRUM AND INTEGRAL OF SPECTRUM.
C
      AVSOUR=(YB(1)+YC(1))
      XCNRM1=XCNRM1+DE*AVSOUR
      XCINT1=XCINT1+DE*(AVSOUR*(YB(2)+YC(2))+
     1 (YB(1)-YC(1))*(YB(2)-YC(2))/3.0)
C
C     TEST FOR END OF GROUP. END OF INTEGRAL OVER GROUP IF LAST ENERGY
C     INTERVAL (XC TO XB) EXTENDS TO THE UPPER BOUNDARY OF THE GROUP
C     (XE). IF NOT END OF GROUP CONTINUE WITH INTEGRALS. IF END OF
C     GROUP DEFINE AVERAGE VALUE AS RATIO OF INTEGRALS. IF THIS IS THE
C     TOTAL CROSS SECTION SAVE UNSHIELDED GROUP AVERAGE VALUE FOR LATER
C     SELF-SHIELDING CALCULATION (SEE GROUPS).
C
  140 IF(XB.LT.XE) GO TO 100
      IF(XCNRM1.GT.0.0) GO TO 150
      AVN(IGR)=0.0
      GO TO 160
  150 AVN(IGR)=0.5*XCINT1/XCNRM1
  160 IF(NTYPE.EQ.2) TOTAV(IGR)=AVN(IGR)
C
C     END OF GROUP LOOP.
C
  170 CONTINUE
      GO TO 190
C
C     GROUP STRUCTURE AND CROSS SECTION ENERGY RANGES DO NOT OVERLAP.
C     DEFINE ONE GROUP (THE HIGHEST ENERGY GROUP) WITH A CROSS SECTION
C     EQUAL TO ZERO.
C
  180 IGRLOW=NGR
      AVN(NGR)=0.0
C
C     SECTION HAS BEEN PROCESSED. PERFORM ALL REQUESTED OUTPUT.
C
C-----OUTPUT IN ENDF/B FORMAT (IF REQUESTED)
  190 IF(OTAPE.GT.0) CALL TAB1(EGROUP,AVN,IGRLOW)
C-----OUTPUT IN LISTING FORMAT (IF REQUESTED)
      IF(LIST3.GT.0) CALL LISTAV(EGROUP,AVN,IGRLOW)
      RETURN
      END
      SUBROUTINE GROUPS
C=======================================================================
C
C     THIS ROUTINE IS DESIGNED TO COMPUTE SHIELDED AND
C     UNSHIELDED GROUP AVERAGED CROSS SECTIONS FOR TOTAL,
C     ELASTIC, CAPTURE AND FISSION CROSS SECTIONS. BOTH THE ENERGY
C     DEPENDENT WEIGHTING SPECTRUM AND ALL CROSS SECTIONS MUST BE LINEAR
C     INTERPOLABLE BETWEEN TABULATED VALUES. ALL INTEGRALS ARE
C     PERFORMED ANALYTICALLY.
C
C     ALTHOUGH ALL INTEGRALS ARE PERFORMED ANALYTICALLY WHEN PERFORMING
C     THE INTEGRAL (E1 TO E2) (SIGMAI(E)*S(E)/(SIGMAT(E)+SIGMAO)**N))
C     SPECIAL CARE MUST BE TAKEN WHEN THE TOTAL CROSS SECTION (SIGMAT)
C     IS CONSTANT OR ALMOST CONSTANT OVER ANY ENERGY INTERVAL OF
C     INTEGRATION. IN THIS CASE THE ANALYTICAL INTERGALS THAT ARE
C     OBTAINED CONTAIN A GREAT DEAL OF CROSS CANCELLATION OF TERMS. IN
C     ORDER TO ACCURATELY HANDLE THIS CASE AN EXPANSION IS USED TO
C     EVALUATE THE EXPRESSION...
C
C     (ALOG((1+R)/(1-R))-2.0*R)/(R**3)
C
C     WHERE R IS THE RATIO OF THE CHANGE IN TOTAL CROSS SECTION ACROSS
C     THE ENERGY INTERVAL TO TWICE THE AVERAGE TOTAL CROSS SECTION IN
C     THE INTERVAL....
C
C     R=(SIGMAT(E2)-SIGMAT(E1))/(SIGMAT(E2)+SIGMAT(E1))
C
C     DEFINING X=R**2 THE APPROPRIATE EXPANSION IS
C
C     2.0*(1/3+X/5+X**2/7+X**3/9+X**4/11+X**5/13+X**6/15+...
C
C     THIS EXPANSION IS TRUNCATED AFTER R**8 (ERROR OF ORDER R**10).
C     THEREFORE FOR ABS(R) LESS THAN 0.01 THIS EXPANSION WILL YIELD
C     16 DIGIT ACCURACY (ERROR OF ORDER 10**(-20)).
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OTAPE2,TYPSIG,OPS
      CHARACTER*1 ZABCD,FIELDX
      CHARACTER*4 STAR,BLANK,TUNITS,TMPHOL
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
      COMMON/OPUS/OPS(5)
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
      COMMON/ELPASD/TUNITS(2,4),TMPHOL(3)
      COMMON/ELPASZ/ZABCD(10)
      COMMON/ELPAS2/EGB,TMPK
      COMMON/MINOK/OKMIN(5),REDMIN(5)
      COMMON/ZABAD/IZABAD
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      COMMON/COMXTEND/MYXTEND
      INCLUDE 'groupie.h'
      DIMENSION ITSIG(12),ERBIG(5),YBB(5),YCC(5),XCFIS(11)
C-----DEFINE INDICES TO SIGMA0 SHIELDED VALUES TO LIST IN OUTPUT.
      DATA ITSIG/1,4,6,8,10,12,14,16,18,20,23,24/
      DATA STAR/'*'/
      DATA IPASS/0/
      DATA BLANK/'    '/
      DATA ZERO/0.0D+0/
      DATA ONE /1.0D+0/
      DATA TWO /2.0D+0/
C
C     FIRST TIME THROUGH INITIALIZE SECOND HALF OF SIGMA0
C     MULTIPLIER TABLE AND CONSTANTS FOR EXPANSION OF LOG
C     AS DESCRIBED ABOVE.
C
      IF(IPASS.GT.0) GO TO 20
      IPASS=1
C-----DEFINE SECOND HALF OF SIGMA0 TABLE TO BE RECIPROCAL OF FIRST HALF.
      II=22
      DO 10 I=2,11
      SIGMAB(II)=1.0/SIGMAB(I)
   10 II=II-1
C-----DEFINE CONSTANTS FOR USE IN EXPANSION OF LOG (AS DESCRIBED ABOVE).
      F3 =2.0D+0/3.0D+0
      F5 =2.0D+0/5.0D+0
      F7 =2.0D+0/7.0D+0
      F9 =2.0D+0/9.0D+0
      F11=2.0D+0/11.0D+0
      F13=2.0D+0/13.0D+0
      F15=2.0D+0/15.0D+0
C-----THIS IS MINIMUM R**2 BELOW WHICH THE EXPANSION WILL BE USED.
      RATMIN=1.0D-4
C-----SIGMAB(23) MUST BE ZERO IF AV(1/TOT**2) AND AV(1/TOT**3) ARE TO BE
C-----CALCULATED PROPERLY (WARNING....DO NOT CHANGE)
      SIGMAB(23)=0.0
C
C     PRINT TITLE IF USING COMPACT, SIMPLIFIED LISTING FOR 1 GROUP
C     (SPECTRUM AVERAGED) CALCULATIONS.
C
      CALL TOP1
C
C     CANNOT PERFORM SELF-SHIELDING AND MULTIBAND CALCULATION IF
C     TOTAL CROSS SECTION IS NOT GIVEN.
C
   20 IF(NPTAB(2).LE.0) RETURN
C-----INITIALIZE
      DO I=1,25
      DO J=1,6
      XCINT(I,J) = 0.0D+0
      ENDDO
      ENDDO
C
C     OUTPUT IDENTIFICATION INFORMATION FOR THIS MATERIAL.
C
      IZABAD=0
C-----CONVERT KELVIN TEMPERATURE TO OUTPUT FORM.
      CALL OUT9G(TMPK,FIELDX(1,4))
C
C     PRINT NORMAL TITLE IF MORE THAN 1 GROUP.
C
      CALL TOP1
      IF(LIST2.LE.0) GO TO 40
      IF(NSECT.LT.5) GO TO 30
      WRITE(LIST2,700) MATNOW,(FIELDX(M,4),M=1,11),TMPHOL(2),
     1 TMPHOL(3),ZABCD,
     2 (REACT2(1,I),REACT2(2,I),I=2,5)
     3,(REACT2(1,I),REACT2(2,I),I=2,6)
      GO TO 40
   30 WRITE(LIST2,710) MATNOW,(FIELDX(M,4),M=1,11),TMPHOL(2),
     1 TMPHOL(3),ZABCD,
     2 (REACT2(1,I),REACT2(2,I),I=2,4),
     3 (REACT2(1,I),REACT2(2,I),I=2,4),REACT2(1,6),REACT2(2,6)
C-----INITIALIZE ERROR VECTOR FOR THIS MATERIAL.
   40 IF(NBMAX.LT.1) GO TO 60
      DO 50 NB=1,NBMAX
      DO 50 I=1,25
   50 ERMAT(I,NB)=0.0
   60 NBNEED=1
C
C     DEFINE SPECTRUM AND ALL CROSS SECTION VALUES OVER FIRST
C     INTERVAL OF FIRST GROUP. IF CROSS SECTION EXTENDS TO LOWER
C     ENERGY LIMIT OF GROUP DEFINE CROSS SECTION AT LOWER ENERGY
C     GROUP BOUNDARY BY LINEAR INTERPOLATION. IF THRESHOLD OF ANY
C     CROSS SECTION IS ABOVE LOWER ENERGY LIMIT OF FIRST GROUP
C     DEFINE CROSS SECTION TO BE ZERO OVER AN INTERVAL THAT EXTENDS
C     FROM THE LOWER GROUP BOUNDARY TO THE THRESHOLD.
C
      XB=EGROUP(1)
      DO 110 ISECTP=1,NSECT
      KSECT=ISECTP
C-----FIND FIRST DATA POINT ABOVE LOWER ENERGY BOUNDARY OF FIRST GROUP.
      IPTAB(KSECT)=0
   70 IF(IPTAB(KSECT).LT.NPTAB(KSECT)) GO TO 80
      XA(KSECT)=2.0D+20
      GO TO 110
   80 IPTAB(KSECT)=IPTAB(KSECT)+1
      CALL XYPAGE(IPTAB(KSECT),KSECT,XA(KSECT),YA(KSECT))
      IF(XA(KSECT).gt.XB) go to 90
C-----POINT IS STILL BELOW LOWER ENERGY LIMIT. SAVE ENERGY AND CROSS
C-----SECTION FOR INTERPOLATION.
      XC=XA(KSECT)
      YB(KSECT)=YA(KSECT)
      GO TO 70
C-----POINT IS ABOVE LOWER ENERGY LIMIT. IF THIS IS FIRST POINT
C-----DEFINE PRECEEDING INTERVAL WITH ZERO CROSS SECTION. OTHERWISE
C-----DEFINE CROSS SECTION AT LOWER ENERGY LIMIT BY LINEAR
C-----INTERPOLATION BETWEEN LAST TWO POINTS.
   90 IF(IPTAB(KSECT).GT.1) GO TO 100
      YA(KSECT)=0.0
      YB(KSECT)=0.0
      IPTAB(KSECT)=0
      GO TO 110
  100 YB(KSECT)=((XB-XC)*YA(KSECT)+(XA(KSECT)-XB)*YB(KSECT))/
     1 (XA(KSECT)-XC)
  110 CONTINUE
C
C     GROUP LOOP. DEFINE LOWER AND UPPER ENERGY LIMITS OF GROUP.
C     INITIALIZE INTEGRAL OF CROSS SECTION TIMES SPECTRUM TIMES
C     SELF-SHIELDING FACTOR (XCINT) AND INTEGRAL OF SPECTRUM TIMES
C     SELF-SHIELDING FACTOR (XCNORM). DEFINE SIGMA0 VALUES AS
C     MULTIPLES OF UNSHIELDED TOTAL CROSS SECTION FOR CURRENT
C     GROUP OR THE SAME BARNS VALUES IN EACH GROUP (SHIELD).
C     INITIALIZE CROSS SECTION LIMITS FOR EACH REACTION (YLOW,YHIGH).
C
      DO 460 IGR=1,NGR
C-----INITIALIZE TOTAL OR ELASTIC NOT POSITIVE FLAG
      NOTPLUS = 0
      EGB=XB
      XE=EGROUP(IGR+1)
      DO 120 L=1,25
  120 XCNORM(L)=0.0
      DO 130 KSECT=2,NSECT
      YLOW(KSECT)=YB(KSECT)
      YHIGH(KSECT)=YB(KSECT)
      DO 130 L=1,25
  130 XCINT(L,KSECT)=0.0
C
C     DEFINE SIGMA-0 TO BE EITHER,
C     (1) MULTIPLE OF UNSHIELDED TOTAL IN EACH GROUP, OR,
C     (2) THE SAME NUMBER OF BARNS IN EACH GROUP.
C
C     IN ALL CASES SIGMAB IS THE SIGMA-0 DIVIDED BY THE UNSHIELDED TOTAL
C
      AVXCGP=TOTAV(IGR)
      IF(TYPSIG.NE.0) GO TO 150
      DO 140 L=2,23
  140 SHIELD(L)=SIGMAB(L)*AVXCGP
      GO TO 190
  150 IF(AVXCGP.GT.0.0) GO TO 170
      DO 160 L=2,23
  160 SIGMAB(L)=0.0
      GO TO 190
  170 DO 180 L=2,23
  180 SIGMAB(L)=SHIELD(L)/AVXCGP
C
C     DEFINE BEGINNING OF NEXT INTERVAL TO BE THE SAME AS THE END OF
C     THE LAST INTERVAL. SELECT LONGEST ENERGY INTERVAL (XC TO XB)
C     WITHIN CURRENT GROUP OVER WHICH THE SPECTRUM AND ALL CROSS
C     SECTIONS ARE LINEARLY INTERPOLABLE (THIS IS THE MINIMUM ENERGY
C     OF XE AND THE XA ARRAY).
C
  190 XC=XB
      XB=XE
      DO 200 KSECT=1,NSECT
      YC(KSECT)=YB(KSECT)
      IF(XA(KSECT).LT.XB) XB=XA(KSECT)
  200 CONTINUE
C
C     NEXT CALCULATION WILL BE INTEGRAL BETWEEN ENERGIES XC AND XB.
C     INTERPOLATE OVER THE ENERGY INTERVAL XC TO XA TO DEFINE THE
C     SPECTRUM AND ALL CROSS SECTIONS AT ENERGY XB. IF THE ENTIRE
C     ENERGY INTERVAL HAS BEEN USED UP (XA.LE.XB) ADVANCE TO THE
C     NEXT ENERGY INTERVAL. SKIP ENERGY INTERVALS WITH NEGATIVE CROSS
C     SECTION AND SKIP ALL EMPTY INTERVALS (ZERO LENGTH DUE TO
C     DISCONTINUTY IN THE SPECTRUM AND/OR CROSS SECTIONS).
C
      DO 230 KSECT=1,NSECT
      IF(XA(KSECT).gt.XB) go to 220
C-----NEXT INTEGRAL WILL EXTEND TO END OF THE ENERGY INTERVAL XC TO
C-----XA. DEFINE INTERPOLATED VALUE AS END INTERVAL AND THEN SELECT
C-----THE NEXT ENERGY INTERVAL. IF BEYOND TABULATED ENERGY RANGE EXTEND
C-----AS CONSTANT (I.E. KEEP SAME VALUE FOR YA AND DEFINE NEW ENERGY
C----- XA TO 2.0D+20 EV).
      YB(KSECT)=YA(KSECT)
      IF(IPTAB(KSECT).LT.NPTAB(KSECT)) GO TO 210
c-----2015/8/6 - Insert another point at end with either
c-----           same or 0 cross section
      IF(IPTAB(KSECT).eq.NPTAB(KSECT)) then
      IPTAB(KSECT)=IPTAB(KSECT)+1
      IF(MYXTEND.LE.0) YA(KSECT) = 0.0d0
      GO TO 230
      ENDIF
c-----Beyond tabulated data.
      XA(KSECT)=2.0D+20
      IF(MYXTEND.LE.0) THEN
      YA(KSECT) = 0.0d0
      YB(KSECT) = 0.0d0
      GO TO 230           ! Cross section = 0 - skip calculation
      ENDIF
      GO TO 220           ! Use extension
c-----Next tabulated point.
  210 IPTAB(KSECT)=IPTAB(KSECT)+1
      CALL XYPAGE(IPTAB(KSECT),KSECT,XA(KSECT),YA(KSECT))
      GO TO 230
C-----NEXT INTEGRAL WILL ONLY BE OVER A PORTION OF THE ENERGY INTERVAL.
C-----INTERPOLATE BETWEEN XC AND XA TO DEFINE VALUE AT XB.
  220 YB(KSECT)=((XB-XC)*YA(KSECT)+(XA(KSECT)-XB)*YB(KSECT))/
     1 (XA(KSECT)-XC)
  230 CONTINUE
      DE=XB-XC
      IF(DE.LE.0.0) GO TO 320
C-----DEFINE ALL END POINTS IN SCRATCH ARRAY.
      DO 240 KSECT=1,NSECT
      YBB(KSECT)=YB(KSECT)
  240 YCC(KSECT)=YC(KSECT)
C-----------------------------------------------------------------------
C
C     04/09/12 - SMALL TOTAL TREATMENT NO LONGER USED.
C
C-----------------------------------------------------------------------
C-----SET FLAG IF TOTAL OR ELASTIC IS NOT POSITIVE.
      IF(YB(2).LE.0.0D+0.OR.YC(2).LE.0.0D+0) NOTPLUS = 1
      IF(YB(3).LE.0.0D+0.OR.YC(3).LE.0.0D+0) NOTPLUS = 1
C
C     KEEP TRACK OF MINIMUM AND MAXMIMUM CROSS SECTION FOR REACTIONS.
C
      DO 260 KSECT=2,NSECT
      IF(YLOW(KSECT).lt.YBB(KSECT)) go to 250
      IF(YLOW(KSECT).eq.YBB(KSECT)) go to 260
      YLOW(KSECT)=YBB(KSECT)
      GO TO 260
  250 IF(YHIGH(KSECT).LT.YBB(KSECT)) YHIGH(KSECT)=YBB(KSECT)
  260 CONTINUE
C
C     COMPUTE EFFECT OF SELF-SHIELDING FACTOR 1/(TOTAL(E)+SIGMA0)**N
C     FOR 25 COMBINATIONS OF (SIGMA0,N). DEFINE INTEGRAL OF SPECTRUM
C     TIMES SELF-SHIELDING FACTOR (XCNORM) FOR ALL 25 VALUES OF
C     (SIGMA0,N). THE SELF-SHIELDED TOTAL CROSS SECTIONS CAN LATER BE
C     DEFINED IN TERMS OF COMBINATIONS OF THESE INTEGRALS. THE EFFECT
C     OF THE SELF-SHIELDING FACTOR ON THE INDIVIDUAL REACTION CROSS
C     SECTIONS MAY BE DIVIDED INTO THREE COMPONENTS..
C
C     (1) FST   - EFFECT ON AVERAGE OF SPECTRUM TIMES AVERAGE OF
C                 REACTION CROSS SECTION.
C     (2) DFST  - EFFECT ON AVERAGE OF SPECTRUM TIMES CHANGE IN
C                 REACTION CROSS SECTION AND AVERAGE REACTION CROSS
C                 SECTION TIMES CHANGE IN SPECTRUM.
C     (3) DDFST - EFFECT ON CHANGE IN SPECTRUM TIMES CHANGE IN REACTION
C                 CROSS SECTION.
C
C     THESE FUNCTIONS ARE CALCULATED BELOW IN A DIMENSIONLESS FORM THAT
C     IS ONLY A FUNCTION OF THE RATIO OF THE CHANGE IN THE TOTAL CROSS
C     SECTION TO TWICE THE AVERAGE TOTAL CROSS SECTION.
C
C     RATIO=(SIGMAT(E2)-SIGMAT(E1))/(SIGMAT(E2)+SIGMAT(E1))
C
C     NOTE THAT AS LONG AS THE TOTAL CROSS SECTION IS POSITIVE THIS
C     RATIO IS RESTRICTED TO THE RANGE -1.0 TO 1.0.
C
C-----DEFINE COMPONENTS OF SPECTRUM (AVERAGE AND CHANGE).
      AVSOUR=YBB(1)+YCC(1)
      DSOUR=YBB(1)-YCC(1)
C-----DEFINE CHANGE IN TOTAL CROSS SECTION (SAME FOR ALL SIGMA0).
      DST1=YBB(2)-YCC(2)
C-----N=1 AND 22 VALUES OF SIGMA0 (I.E. 1/(TOTAL(E)+SIGMA0)).
      DO 290 L=2,23
      YBP=YBB(2)+SHIELD(L)
      YCP=YCC(2)+SHIELD(L)
      AVST1=YBP+YCP
C-----ALLOW FOR NEGATIVE TOTAL
      IF(AVST1.GT.0.0D+0) THEN
      DEAVST(L)= DE/AVST1
      RATIO    = DST1/AVST1
      ELSE
      DEAVST(L)= ZERO
      RATIO=     ZERO
      ENDIF
C
      RATIO2=RATIO**2
      IF(RATIO2.LT.RATMIN) GO TO 270
      XX=(ONE+RATIO)/(ONE-RATIO)
      XX=DLOG(XX)/RATIO
      FST(L)=XX
      XX=(TWO-XX)/RATIO
      DFST(L)=XX
      DDFST(L)=-XX/RATIO
      GO TO 280
  270 XX=(((((F15*RATIO2+F13)*RATIO2+F11)*RATIO2+F9)*RATIO2+F7)*
     1                                RATIO2+F5)*RATIO2+F3
      DDFST(L)=XX
      DFST(L)=-RATIO*XX
      FST(L)=TWO+RATIO2*XX
  280 CONTINUE
  290 CONTINUE
C-----SIGMA0=0, N=2 AND 3 (I.E. 1/TOTAL(E)**2 AND 1/TOTAL(E)**3). THE
C-----FOLLOWING CODING ASSUMES THAT SHIELD(23) IS 0.0. THEREFORE AVST1
C-----AND RATIO NEED NOT BE CALCULATED AGAIN IN THIS SECTION (THEY ARE
C-----THE SAME AS CALCULATED IN THE LAST PASS OF THE ABOVE DO LOOP).
      FST(24)=2.0/(1.0-RATIO2)
      DFST(24)=RATIO*(DDFST(23)-FST(24))
      DDFST(24)=FST(24)-2.0*DDFST(23)
      FST(25)=0.5*(FST(24)**2)
      DFST(25)=-RATIO*FST(25)
      DDFST(25)=DDFST(23)+RATIO2*FST(25)
C
C     NORMALIZE DIMENSIONLESS INTEGRALS AND THEN DEFINE INTEGRAL OF
C     SPECTRUM TIMES SELF-SHIELDING FACTOR (XCNORM) TO USE AS
C     NORMALIZATION FOR ALL REACTIONS.
C
      DEAVST(24)=DEAVST(23)/AVST1
      DEAVST(25)=DEAVST(24)/AVST1
      XCNORM(1)=XCNORM(1)+DE*AVSOUR
      DO 300 L=2,25
      FST(L)=DEAVST(L)*FST(L)
      DFST(L)=DEAVST(L)*DFST(L)
      DDFST(L)=DEAVST(L)*DDFST(L)
  300 XCNORM(L)=XCNORM(L)+(AVSOUR*FST(L)+DSOUR*DFST(L))
C
C     FOR EACH REACTION (TOTAL, ELASTIC, CAPTURE AND FISSION) DEFINE
C     INTEGRAL OF THE SPECTRUM TIMES CROSS SECTION TIMES SELF-SHIELDING
C     FACTOR FOR ALL 25 COMBINATION OF (SIGMA0,N).
C
      DO 310 KSECT=2,NSECT
      AVSTN=YBB(KSECT)+YCC(KSECT)
      DSTN=YBB(KSECT)-YCC(KSECT)
      A1=AVSTN*AVSOUR
      A2=DSTN*AVSOUR+DSOUR*AVSTN
      A3=DSTN*DSOUR
      XCINT(1,KSECT)=XCINT(1,KSECT)+DE*(A1+A3/3.0)
      DO 310 L=2,25
  310 XCINT(L,KSECT)=XCINT(L,KSECT)+(A1*FST(L)+A2*DFST(L)+A3*DDFST(L))
C
C     TEST FOR END OF GROUP. IF NOT, CONTINUE WITH INTEGRALS.
C     IF END OF GROUP, NORMALIZE ALL INTEGRALS AND DEFINE COMMON
C     FACTORS. DEFINE SELF-SHIELDED AVERAGED CROSS SECTIONS AS
C     RATIO OF INTEGRALS (REACTIONS TO FLUX).
C
  320 IF(XB.LT.XE) GO TO 190
C-----CONVERT TOTAL AND REACTION INTEGRALS TO AVERAGES.
      DO 340 I=1,25
C-----IF NORMALIZATION = 0 SET VALUES = 0
      IF(XCNORM(I).GT.0.0) GO TO 330
      DO KSECT=2,NSECT
      XCINT(I,KSECT)=0.0
      ENDDO
      GO TO 340
C-----AT THIS POINT THE NUMERATOR (XCINT) AND DENOMINATOR (XCNORM)
C-----DIFFER FROM THE EXACT INTEGRALS BY SCALAR FACTORS.....
C-----XCINT(I),I=1,23 - IS 4 TIMES THE CORRECT INTEGRAL
C-----XCINT(24)       - IS 2 TIMES THE CORRECT INTEGRAL
C-----XCINT(25)       - IS 1 TIMES THE CORRECT INTEGRAL
C-----XCNORM(I),I=1,23- IS 2 TIMES THE CORRECT INTEGRAL
C-----XCNORM(24)      - IS 1 TIMES THE CORRECT INTEGRAL
C-----XCNORM(25)      - IS 1/2 TIMES THE CORRECT INTEGRAL.
C-----RENORMALIZE XCNORM TO OBTAIN THE CORRECT RATIOS.
  330 XCNORM(I)=2.0*XCNORM(I)
      DO KSECT=2,NSECT
      XCINT(I,KSECT)=XCINT(I,KSECT)/XCNORM(I)
      ENDDO
  340 CONTINUE
C
C     CHECK TO INSURE TOTAL DOES NOT INCREASE WITH SIGMA0
C
      IF(NOTPLUS.EQ.0) THEN
      DO KK=2,25
      IF(XCINT(KK,2).GT.1.000001*XCINT(KK-1,2).OR.
     1   XCINT(KK,2).GT.1.000001*XCINT(1,2)) THEN
      WRITE(3,350) IGR
      WRITE(*,350) IGR
  350 FORMAT(' IGR=',I5)
      WRITE(3,360) (II,XCINT(II,2),XCINT(II,2)/XCINT(1,2),II=1,25)
      WRITE(*,360) (II,XCINT(II,2),XCINT(II,2)/XCINT(1,2),II=1,25)
  360 FORMAT(I8,1PD20.12,0PF20.12,' BEFORE SUM')
      ENDIF
      ENDDO
      ENDIF
C
C     DEFINE THE REST = UNSHIELDED [TOTAL-(ELASTIC+CAPTURE+FISSION)]
C
      DO I=1,25
      THEREST=0.0
      DO KSECT=3,NSECT
      THEREST=THEREST+XCINT(I,KSECT)
      ENDDO
      XCINT(I,6) = XCINT(I,2) - THEREST
c-----2016/7/3 - Added lower limit for "Other" -
C-----           Insure TOTAL is EXACTLY SUM OF PARTS
      if(XCINT(I,6).lt.1.0d-6*XCINT(I,2)) then
      XCINT(I,6) = 0.0d0
      XCINT(I,2) = THEREST
      endif
      ENDDO
C
C     IF NON-POSITIVE TOTAL OR ELASTIC DEFINE ALL = UNSHIELDED
C
C-----IF LITTLE SELF-SHIELDING SET TO NONE
      IF(XCINT(23,2).GT.0.99999D+0*XCINT(1,2)) NOTPLUS = 1
      IF(NOTPLUS.NE.0) THEN
      DO KSECT=2,6
      DO I=2,25
      XCINT(I,KSECT) = XCINT(1,KSECT)
      ENDDO
      ENDDO
      ENDIF
C
C     DEFINE F FACTORS AND OUTPUT.
C
C-----DEFINE LOWER GROUP BOUNDARY.
      CALL OUT9G(EGB,FIELDX(1,1))
C-----LOOP OVER TOTAL AND PARTIALS.
      DO 450 KSECT=2,6
C-----ASSUME SMALL "OTHER" IS NOISE AND IGNORE
C-----2016/5/21 - Changed multiple IF statement to accommodate compiler
C-----            optimizer.
C     IF(KSECT.EQ.6.AND.
C    1   XCINT(1,KSECT).LE.1.0D-6*XCINT(1,2)) GO TO 450
c-----The following replaces the above
      IF(KSECT.EQ.6) THEN        ! Only execute the next line if KSECT=6
      IF(XCINT(1,KSECT).LE.1.0D-6*XCINT(1,2)) GO TO 450
      ENDIF
C-----CROSS SECTION.
      XCINTS=XCINT(1,KSECT)
      CALL OUT9G(XCINTS,FIELDX(1,2))
C-----CROSS SECTION OR RESONANCE INTEGRAL.
      IF(OPS(1).EQ.2) XCINTS=DLOG(EGROUP(IGR+1)/EGROUP(IGR))*XCINTS
      CALL OUT9G(XCINTS,FIELDX(1,3))
      IF(XCINT(1,KSECT).GT.0.0) GO TO 380
      XCFI(1,KSECT)=1.0
      DO 370 L=2,25
  370 XCFI(L,KSECT)=0.0
      GO TO 400
  380 DO 390 L=2,25
  390 XCFI(L,KSECT)=XCINT(L,KSECT)/XCINT(1,KSECT)
C
C     2010/3/12 - SAVE PLOTTAB OUTPUT
C
  400 XCPLOT(IGR,KSECT-1,1) = XCINT( 1,KSECT) ! UNSHIELDED
      XCPLOT(IGR,KSECT-1,2) = XCINT(23,KSECT) !   SHIELDED
C
C     OUTPUT.
C
      DO 410 L=1,11
      LOUT = ITSIG(L+1)
  410 XCFIS(L)=XCFI(LOUT,KSECT)
      IF(KSECT.EQ.2) THEN
      DO KKK=1,11
      IF(XCFIS(KKK).GT.1.000001D+0) THEN
      WRITE(3,420) (LLL,XCFIS(LLL),LLL=1,11)
      WRITE(*,420) (LLL,XCFIS(LLL),LLL=1,11)
      WRITE(3,420) (LLL,XCFI(LLL,2),LLL=1,25)
      WRITE(*,420) (LLL,XCFI(LLL,2),LLL=1,25)
  420 FORMAT(I8,0PF20.12)
      ENDIF
      ENDDO
      ENDIF
C-----INCLUDE GROUP NUMBER AND LOWER ENERGY WITH TOTAL.
      IF(KSECT.NE.2) GO TO 430
      IF(LIST1.GT.0) WRITE(LIST1,650) IGR,(FIELDX(M,1),M=1,11),
     1 REACT3(1,KSECT),REACT3(2,KSECT),
     2 (FIELDX(M,3),M=1,11),XCFIS
      GO TO 440
  430 IF(LIST1.GT.0) WRITE(LIST1,660) REACT3(1,KSECT),REACT3(2,KSECT),
     1 (FIELDX(M,3),M=1,11),XCFIS
  440 CONTINUE
  450 CONTINUE
C
C     MULTI-BAND CALCULATION.
C
      IF(METHODB.GT.0) CALL BANDIT
C
C     END OF GROUP LOOP
C
  460 CONTINUE
C
C     PLOTTAB OUTPUT
C
C-----TOTAL
      WRITE(16,470) ZABCD
  470 FORMAT(10A1,'   Total',' Unshielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,1,1),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1               (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
  480 FORMAT(22A1)
      ENDDO
      WRITE(16,490)
  490 FORMAT(30X,'(BLANK LINE)')
      WRITE(16,500) ZABCD
  500 FORMAT(10A1,'   Total',' Shielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,1,2),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1               (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
C-----ELASTIC
      WRITE(16,510) ZABCD
  510 FORMAT(10A1,' Elastic',' Unshielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,2,1),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1              (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
      WRITE(16,520) ZABCD
  520 FORMAT(10A1,' Elastic',' Shielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,2,2),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1              (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
C-----CAPTURE
      IF(NSECT.GE.4) THEN
      WRITE(16,530) ZABCD
  530 FORMAT(10A1,' Capture',' Unshielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,3,1),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1              (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
      WRITE(16,540) ZABCD
  540 FORMAT(10A1,' Capture',' Shielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,3,2),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1              (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
      ENDIF
C-----FISSION
      IF(NSECT.GE.5) THEN
      WRITE(16,550) ZABCD
  550 FORMAT(10A1,' Fission',' Unshielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,4,1),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1              (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
      WRITE(16,560) ZABCD
  560 FORMAT(10A1,' Fission',' Shielded')
      DO IGG=1,NGR
      CALL OUT9G(EGROUP(IGG    ),FIELDX(1,1))
      CALL OUT9G(XCPLOT(IGG,4,2),FIELDX(1,2))
      CALL OUT9G(EGROUP(IGG+1  ),FIELDX(1,3))
      WRITE(16,480) (FIELDX(I,1),I=1,11),(FIELDX(I,2),I=1,11),
     1              (FIELDX(I,3),I=1,11),(FIELDX(I,2),I=1,11)
      ENDDO
      WRITE(16,490)
      ENDIF
C
C     FINISH LAST PAGE OF LISTINGS FOR THIS MATERIAL.
C
      CALL OUT9G(XE,FIELDX(1,1))
C-----FINISH PAGE OF SELF-SHIELDED CROSS SECTIONS.
      IF(LIST1.LE.0.OR.NGR.LE.1) GO TO 570
      WRITE(LIST1,650) NGRP1,(FIELDX(M,1),M=1,11)
      WRITE(LIST1,670) BLANK
C-----FINISH PAGE OF MULTI-BAND PARAMETERS.
  570 IF(LIST2.LE.0) GO TO 600
      IF(NSECT.LT.5) GO TO 580
      WRITE(LIST2,680) NGRP1,(FIELDX(M,1),M=1,11),STAR
      GO TO 590
  580 WRITE(LIST2,690) NGRP1,(FIELDX(M,1),M=1,11),STAR
  590 WRITE(LIST2,670) BLANK
C
C     UPDATE MAXMIMUM ERROR FOR THE ENTIRE LIBRARY FOR EACH NUMBER
C     OF BANDS. DEFINE MAXIMUM ERROR FOR CURRENT MATERIAL FOR EACH
C     NUMBER OF BANDS AND OUTPUT THE RESULTS.
C
  600 IF(METHODB.EQ.0) GO TO 640
      DO 620 NB=1,NBNEED
      ERBIG(NB)=0.0
      DO 610 I=1,25
      IF(ERMAT(I,NB).GT.ERLIB(I,NB)) ERLIB(I,NB)=ERMAT(I,NB)
      IF(I.GT.23) GO TO 610
      IF(ERMAT(I,NB).GT.ERBIG(NB)) ERBIG(NB)=ERMAT(I,NB)
  610 CONTINUE
  620 ERBIG(NB)=100.0*ERBIG(NB)
C-----SAVE PARAMETERS FOR FINAL OUTPUT REPORT.
      LZA=LZA+1
      IF(LZA.GT.MAXMAT) THEN
      WRITE(3,720) MAXMAT
      WRITE(*,720) MAXMAT
      CALL ENDERROR
      ENDIF
      IZATAB(LZA)=IZA
      NBNTAB(LZA)=NBNEED
      DO 630 I=1,NBNEED
  630 ERBTAB(I,LZA)=ERBIG(I)
  640 CONTINUE
C
C     MULTI-BAND LIBRARY.
C
      IF(OTAPE2.GT.0) CALL BANOUT
      RETURN
  650 FORMAT(I5,11A1,1X,A4,A1,11A1,11F8.5)
  660 FORMAT(       17X,A4,A1,11A1,11F8.5)
  670 FORMAT(73X,A1)
  680 FORMAT(    I5,11A1,58X,A1)
  690 FORMAT( 9X,I5,11A1,48X,A1)
  700 FORMAT(' MAT',I5,11X,11A1,1X,
     1 A4,A3,' Multi-Band Parameters',12X,'*',
     2 7X,' Unshielded Group Averages',12X,10A1/73X,'*'/
     1 '  No.   Group-eV Band  Weight',4(3X,2A4),' *',5(3X,2A4)/73X,'*')
  710 FORMAT(' MAT',I5,11X,11A1,1X,
     1 A4,A3,' Multi-Band Parameters',12X,'*',
     2 7X,' Unshielded Group Averages',12X,10A1/73X,'*'/
     1 10X,'  No.   Group-eV Band  Weight',3(3X,2A4),' *',4(3X,2A4)/73X
     2,'*')
  720 FORMAT(//' ERROR - More than',I6,' Materials'//)
      END
      SUBROUTINE TOP1
C=======================================================================
C
C     PRINT HEADING FOR SHIELDED OUTPUT LISTING.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OTAPE2,TYPSIG,OPS
      CHARACTER*1 ZABCD,FIELDX
      CHARACTER*4 TUNITS,TMPHOL
      CHARACTER*8 SIGHL1,SIGHL2,RICHAR
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/OPUS/OPS(5)
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
      COMMON/ELPASD/TUNITS(2,4),TMPHOL(3)
      COMMON/ELPASZ/ZABCD(10)
      COMMON/ELPAS2/EGB,TMPK
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      COMMON/SIGGIE/SIGHL1(12),SIGHL2(12)
      INCLUDE 'groupie.h'
      DATA RICHAR/'R.I.    '/
      DATA IPASS/0/
C
C     NOTHING TO DO IF NO SHIELDED OUTPUT LISTING.
C
      IF(LIST1.LE.0) RETURN
C
C     ON FIRST CALL SELECT EITHER COMPACT, SIMPLIFIED LISTING IF ONLY
C     1 GROUP (SPECTRUM AVERAGES) OR NORMAL LISTING IF MORE THAN
C     1 GROUP.
C
      IF(IPASS.GT.0) GO TO 30
      IPASS=1
C
C     IF OUTPUT IS RESONANCE INTEGRALS CHANGE TITLE.
C
      IF(OPS(1).NE.2) GO TO 10
      SIGHL1(1)=RICHAR
      SIGHL2(1)=RICHAR
   10 IF(NGR.GT.1) RETURN
C
C     USE COMPACT, SIMPLIFIED FORMAT FOR 1 GROUP RESULTS.
C
      CALL OUT9G(TMPK,     FIELDX(1,1))
      CALL OUT9G(EGROUP(1),FIELDX(1,2))
      CALL OUT9G(EGROUP(2),FIELDX(1,3))
      IF(TYPSIG.NE.0) GO TO 20
      IF(OPS(1).EQ.1) WRITE(LIST1,90) (FIELDX(M,1),M=1,11),
     1 TMPHOL(2),TMPHOL(3),((FIELDX(M,J),M=1,11),J=2,3),SIGHL1
      IF(OPS(1).EQ.2) WRITE(LIST1,100) (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),((FIELDX(M,J),M=1,11),J=2,3),SIGHL1
      RETURN
C-----IDENTIFY CROSS SECTIONS OR RESONANCE INTEGRALS.
   20 IF(OPS(1).EQ.1) WRITE(LIST1,110) (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),((FIELDX(M,J),M=1,11),J=2,3),SIGHL2
      IF(OPS(1).EQ.2) WRITE(LIST1,120) (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),((FIELDX(M,J),M=1,11),J=2,3),SIGHL2
      RETURN
C-----IDENTIFY CROSS SECTIONS OR RESONANCE INTEGRALS.
C
C     USE NORMAL OUTPUT LISTING IF MORE THAN 1 GROUP.
C
C-----CONVERT KELVIN TEMPERATURE TO OUTPUT FORM.
   30 IF(NGR.LE.1) RETURN
      CALL OUT9G(TMPK,FIELDX(1,4))
      IF(TYPSIG.NE.0) GO TO 40
C-----IDENTIFY CROSS SECTIONS OR RESONANCE INTEGRALS.
      IF(OPS(1).EQ.1) WRITE(LIST1,50) MATNOW,ZABCD,
     1 (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),ZABCD,SIGHL1
      IF(OPS(1).EQ.2) WRITE(LIST1,60) MATNOW,ZABCD,
     1 (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),ZABCD,SIGHL1
      RETURN
C-----IDENTIFY CROSS SECTIONS OR RESONANCE INTEGRALS.
   40 IF(OPS(1).EQ.1) WRITE(LIST1,70) MATNOW,ZABCD,
     1 (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),ZABCD,SIGHL2
      IF(OPS(1).EQ.2) WRITE(LIST1,80) MATNOW,ZABCD,
     1 (FIELDX(M,4),M=1,11),
     1 TMPHOL(2),TMPHOL(3),ZABCD,SIGHL2
      RETURN
   50 FORMAT(' MAT',I5,1X,10A1,14X,11A1,1X,
     1 A4,A2,' Self-Shielded Cross Sections',30X,10A1//
     2 21X,'  Unshielded (Sigma-0 = the Unshielded Total Cross',
     3 ' Section in Each Group times Below Multipliers)'/
     4 '  No.   Group-eV React',3X,12A8)
   60 FORMAT(' MAT',I5,1X,10A1,11X,11A1,1X,
     1 A4,A2,' Self-Shielded Resonance Integrals',28X,10A1//
     2 21X,'  Unshielded (Sigma-0 = the Unshielded Total Cross',
     3 ' Section in Each Group times Below Multipliers)'/
     4 '  No.   Group-eV React',3X,12A8)
   70 FORMAT(' MAT',I5,1X,10A1,14X,11A1,1X,
     1 A4,A2,' Self-Shielded Cross Sections',30X,10A1//
     2 21X,'  Unshielded (Sigma-0 = the barns Values Listed',
     3 'Below)'/
     4 '  No.   Group-eV React',3X,12A8)
   80 FORMAT(' MAT',I5,1X,10A1,11X,11A1,1X,
     1 A4,A2,' Self-Shielded Resonance Integrals',28X,10A1//
     2 21X,'  Unshielded (Sigma-0 = the barns Values Listed',
     3 ' Below)'/
     4 '  No.   Group-eV React',3X,12A8)
   90 FORMAT(' Self-Shielded Cross Sections'/
     1 1X,11A1,1X,A4,A2/
     2 1X,11A1,' to',11A1,' eV'//
     3 21X,'  Unshielded (Sigma-0 = the Unshielded Total Cross',
     4 ' Section in Each Group times Below Multipliers)'/
     5 ' Material   MAT React',2X,12A8)
  100 FORMAT(' Self-Shielded Resonance Integrals'/
     1 1X,11A1,1X,A4,A2/
     2 1X,11A1,' to',11A1,' eV'//
     3 21X,'  Unshielded (Sigma-0 = the Unshielded Total Cross',
     4 ' Section in Each Group times Below Multipliers)'/
     5 ' Material   MAT React',2X,12A8)
  110 FORMAT(' Self-Shielded Cross Sections'/
     1 1X,11A1,1X,A4,A2/
     2 1X,11A1,' to',11A1,' eV'//
     3 21X,'  Unshielded (Sigma-0 = the barns Values Listed',
     5 'Below)'/' Material   MAT React',2X,12A8)
  120 FORMAT(' Self-Shielded Resonance Integrals'/
     1 1X,11A1,1X,A4,A2/
     2 1X,11A1,' to',11A1,' eV'//
     3 21X,'  Unshielded (Sigma-0 = the barns Values Listed',
     5 ' Below)'/' Material   MAT React',2X,12A8)
      END
      SUBROUTINE TAB1(E,AV,IGRLOW)
C=======================================================================
C
C     OUTPUT MULTI-GROUP CROSS SECTIONS IN ENDF/B FORMAT IN EITHER
C     HISTOGRAM OR LINEAR-LINEAR INTERPOLABLE FORM.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/LEADER/TEMP,Q,L1,L2,N1,N2,MAT,MF,MT
      COMMON/LAWYER/MYLAWO
      COMMON/FLAGS/MINUS3,IMPLUS
      INCLUDE 'groupie.h'
      DIMENSION E(*),AV(*),NBTO(1),INTO(2)
C-----DEFINE HISTOGRAM AND LINEAR-LINEAR INTERPOLATION LAW.
      DATA INTO/1,2/
C-----INITIALIZE INDEX TO OUTPUT ARRAY.
      II=0
      MINUS3=0
      IMPLUS=0
C-----SELECT HISTOGRAM OR LINEAR-LINEAR INTERPOLABLE FORM.
      IF(MYLAWO.GT.1) GO TO 20
C
C     HISTOGRAM OUTPUT.
C
C-----THERE WILL BE ONE POINT OUTPUT FOR EACH GROUP, PLUS A FINAL
C-----ENERGY AT UPPER ENERGY LIMIT OF LAST GROUP.
      N2OUT=NGR-IGRLOW+2
C-----OUTPUT LEADER CARD AND INTERPOLATION LAW (ONE INTERPOLATION REGION
C-----USING INTERPOLATION LAW 1).
      N1XX=1
      CALL CARDO(TEMP,Q,L1,L2,N1XX,N2OUT)
      NBTO(1)=N2OUT
      CALL TERPO(NBTO,INTO,1)
C-----OUTPUT ONE DATA POINT FOR EACH DATA POINT, PLUS AN EXTRA ENERGY
C-----AS THE UPPER ENERGY LIMIT OF THE LAST GROUP
      DO 10 ILOW=IGRLOW,NGR
      II=II+1
      XOUT(II)=E(ILOW)
      YOUT(II)=AV(ILOW)
      IF(II.LT.MAXOUT) GO TO 10
      CALL POINTO9(XOUT,YOUT,II)
      II=0
   10 CONTINUE
C-----ADD UPPER ENERGY LIMIT OF LAST GROUP AND OUTPUT.
      II=II+1
      XOUT(II)=E(NGR+1)
      YOUT(II)=0.0
      CALL POINTO9(XOUT,YOUT,II)
      GO TO 60
C
C     LINEAR-LINEAR INTERPOLABLE OUTPUT.
C
C-----THERE WILL BE TWO POINTS OUTPUT FOR EACH GROUP.
   20 N2OUT=2*(NGR-IGRLOW+1)
C-----IF CROSS SECTION HAS A THRESHOLD AND DOES NOT SPAN GROUP
C-----STRUCTURE DEFINE AN ADDITIONAL POINT TO DEFINE CROSS SECTION
C-----TO BE ZERO AT LOWER ENERGY LIMIT OF FIRST GROUP WITH NON-ZERO
C-----AVERAGE CROSS SECTION (GROUP IGRLOW).
      IF(IGRLOW.LE.1.OR.AV(IGRLOW).LE.0.0) GO TO 30
      XOUT(1)=E(IGRLOW)
      YOUT(1)=0.0
C-----INCREASE OUTPUT POINT COUNT FOR POINT AT THRESHOLD.
      N2OUT=N2OUT+1
      II=1
C-----OUTPUT LEADER CARD AND INTERPOLATION LAW (ONE INTERPOLATION REGION
C-----USING INTERPOLATION LAW 2).
   30 L1XX=1
      CALL CARDO(TEMP,Q,L1,L2,L1XX,N2OUT)
      NBTO(1)=N2OUT
      CALL TERPO(NBTO,INTO(2),1)
C-----CREATE TWO DATA POINTS FOR EACH ENERGY GROUP (POINTS AT LOWER
C-----AND UPPER ENERGY LIMITS OF GROUP. CROSS SECTION SAME AT BOTH
C-----POINTS).
      DO 50 ILOW=IGRLOW,NGR
      II=II+1
      XOUT(II)=E(ILOW)
      YOUT(II)=AV(ILOW)
      IF(II.LT.MAXOUT) GO TO 40
      CALL POINTO9(XOUT,YOUT,II)
      II=0
   40 II=II+1
      XOUT(II)=E(ILOW+1)
      YOUT(II)=AV(ILOW)
      IF(II.LT.MAXOUT) GO TO 50
      CALL POINTO9(XOUT,YOUT,II)
      II=0
   50 CONTINUE
      IF(II.GT.0) CALL POINTO9(XOUT,YOUT,II)
   60 RETURN
      END
      SUBROUTINE LISTAV(E,AV,IGRLOW)
C=======================================================================
C
C     LIST UNSHIELDED GROUP AVERAGE CROSS SECTION FOR ONE REACTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OTAPE2,OPS
      CHARACTER*1 ZABCD,FIELD1,FIELD2
      CHARACTER*4 QBLANK
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
      COMMON/ELPASZ/ZABCD(10)
      COMMON/OPUS/OPS(5)
      INCLUDE 'groupie.h'
      DIMENSION E(*),AV(*),KEY(4),FIELD1(11,6),FIELD2(11,6)
      DATA QBLANK/'    '/
      DATA IPASS/0/
      DATA LSTMAT/-9999/
C-----DEFINE PARAMETERS FOR PAGE LAYOUT AT FIRST CALL.
      IF(IPASS.GT.0) GO TO 10
      IPASS=1
C-----PRINT TITLE IF USING ONLY 1 GROUP (SPECTRUM AVERAGES).
      CALL TOP3
C
C     SELECT COMPACT OR NORMAL OUTPUT FORMAT.
C
   10 IF(NGRP1.GT.2) GO TO 30
C-----SKIP ZERO OUTPUT.
      IF(AV(1).LE.0.0) GO TO 20
C-----SELECT CROSS SECTION OR RESONANCE INTEGRAL.
      AVUSE=AV(1)
      IF(OPS(5).EQ.2) AVUSE=DLOG(E(2)/E(1))*AVUSE
      CALL OUT9G(AVUSE,FIELD2(1,1))
      IF(MATNOW.NE.LSTMAT)
     1 WRITE(LIST3,130) ZABCD,MATNOW,MTH,(FIELD2(M,1),M=1,11)
      IF(MATNOW.EQ.LSTMAT)
     1 WRITE(LIST3,140) MTH,(FIELD2(M,1),M=1,11)
      LSTMAT=MATH
   20 RETURN
C
C     NORMAL PAGE LAYOUT.
C
C-----DEFINE NUMBER OF VALUES TO PRINT.
   30 LGR=(NGRP1-IGRLOW)+1
C-----DEVIDE INTO 4 COLUMNS.
      LCOL=(LGR+3)/4
      LCOLM1=LCOL-1
      LPAGE=4*LCOL
C-----SET UP LOOP OVER PAGES OF OUTPUT.
      DO 100 I1=IGRLOW,NGRP1,LPAGE
C-----LIST TITLE FOR PAGE.
      CALL TOP3
C-----DEFINE INDEX TO END OF FIRST COLUMN.
      I2=I1+LCOLM1
      IF(I2.GT.NGRP1) I2=NGRP1
C-----SET  UP LOOP OVER LINES OF OUTPUT ON THIS PAGE (UP TO 53 LINES).
      DO 90 I3=I1,I2
C-----SELECT UP TO FOUR GROUPS THAT WILL APPEAR ON THIS LINE.
      KEY(1)=I3
      CALL OUT9G(E(I3),FIELD1(1,1))
      IF(I3.GE.NGRP1) GO TO 40
C-----SELECT CROSS SECTION OR RESONANCE INTEGRAL.
      AVUSE=AV(I3)
      IF(OPS(5).EQ.2) AVUSE=DLOG(E(I3+1)/E(I3))*AVUSE
      CALL OUT9G(AVUSE,FIELD2(1,1))
   40 KK=I3
      DO 50 K=2,4
      KK=KK+LCOL
      IF(KK.GT.NGRP1) GO TO 60
      KEY(K)=KK
      CALL OUT9G(E(KK),FIELD1(1,K))
      IF(KK.GE.NGRP1) GO TO 50
C-----SELECT CROSS SECTION OR RESONANCE INTEGRAL.
      AVUSE=AV(KK)
      IF(OPS(5).EQ.2) AVUSE=DLOG(E(I3+1)/E(I3))*AVUSE
      CALL OUT9G(AVUSE,FIELD2(1,K))
   50 CONTINUE
      K=5
   60 K=K-1
C-----FOR UPPER ENERGY OF LAST GROUP DO NOT PRINT AVERAGE VALUE.
      IF(KEY(K).GE.NGRP1) GO TO 70
      WRITE(LIST3,110) (KEY(L),(FIELD1(M,L),M=1,11),
     1 (FIELD2(M,L),M=1,11),L=1,K)
      GO TO 90
   70 IF(K.LE.1) GO TO 80
      KM1=K-1
      WRITE(LIST3,110) (KEY(L),(FIELD1(M,L),M=1,11),
     1 (FIELD2(M,L),M=1,11),L=1,KM1),
     2 KEY(K),(FIELD1(M,K),M=1,11)
      GO TO 90
   80 WRITE(LIST3,110) NGRP1,(FIELD1(M,1),M=1,11)
   90 CONTINUE
C-----SPACE TO THE BOTTOM OF THE PAGE.
  100 WRITE(LIST3,120) QBLANK
      RETURN
  110 FORMAT(I6,22A1,3(I7,22A1))
  120 FORMAT(A1)
  130 FORMAT(1X,10A1,I5,I4,11A1)
  140 FORMAT(16X,I4,11A1)
      END
      SUBROUTINE TOP3
C=======================================================================
C
C     PRINT HEADING FOR UNSHIELDED OUTPUT LISTING.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OTAPE2,OPS
      CHARACTER*1 ZABCD,FIELD1
      CHARACTER*4 TUNITS,TMPHOL
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
      COMMON/ELPASD/TUNITS(2,4),TMPHOL(3)
      COMMON/ELPASZ/ZABCD(10)
      COMMON/ELPAS2/EGB,TMPK
      COMMON/OPUS/OPS(5)
      INCLUDE 'groupie.h'
      DIMENSION FIELD1(11,3)
      DATA IPASS/0/
C
C     NOTHING TO DO IF NO UNSHIELDED OUTPUT LISTING.
C
      IF(LIST3.LE.0) RETURN
C
C     ON FIRST CALL SELECT EITHER SPECIAL, COMPACT FORMAT IF ONLY
C     1 GROUP OR NORMAL OUTPUT LISTING.
C
      IF(IPASS.GT.0) GO TO 10
      IPASS=1
C
C     SPECIAL, COMPACT FORMAT IF ONLY 1 GROUP.
C
      IF(NGR.GT.1) RETURN
      CALL OUT9G(TMPK     ,FIELD1(1,1))
      CALL OUT9G(EGROUP(1),FIELD1(1,2))
      CALL OUT9G(EGROUP(2),FIELD1(1,3))
C-----IDENTIFY EITHER CROSS SECTIONS OR RESONANCE INTEGRALS.
      IF(OPS(5).EQ.1) WRITE(LIST3,40) (FIELD1(M,1),M=1,11),
     1 TMPHOL(2),TMPHOL(3),((FIELD1(M,J),M=1,11),J=2,3)
      IF(OPS(5).EQ.2) WRITE(LIST3,50) (FIELD1(M,1),M=1,11),
     1 TMPHOL(2),TMPHOL(3),((FIELD1(M,J),M=1,11),J=2,3)
      RETURN
C
C     IF MORE THAN 1 GROUP PRINT NORMAL HEADING.
C
   10 IF(NGR.LE.1) RETURN
C-----LIST TITLE FOR PAGE.
      CALL OUT9G(TMPK     ,FIELD1(1,1))
C-----IDENTIFY EITHER CROSS SECTIONS OR RESONANCE INTEGRALS.
      IF(OPS(5).EQ.1) WRITE(LIST3,20) MTH,MATNOW,(FIELD1(M,1),M=1,11),
     1 TMPHOL(2),TMPHOL(3),ZABCD
      IF(OPS(5).EQ.2) WRITE(LIST3,30) MTH,MATNOW,(FIELD1(M,1),M=1,11),
     1 TMPHOL(2),TMPHOL(3),ZABCD
      RETURN
   20 FORMAT(53X,'MT=',I3/' MAT',I5,
     1 24X,11A1,1X,A4,A2,' Unshielded Cross Sections',
     2 25X,10A1//1X,'  No.  Group-eV   Average',
     3 3(4X,'  No.  Group-eV   Average')/)
   30 FORMAT(53X,'MT=',I3/' MAT',I5,
     1 22X,11A1,1X,A4,A2,' Unshielded Resonance Integrals',
     2 22X,10A1//1X,'  No.  Group-eV   R.I.   ',
     3 3(4X,'  No.  Group-eV   R.I.   ')/)
   40 FORMAT(
     1 ' Unshielded Cross Sections'/
     1 1X,11A1,1X,A4,A2/
     3 1X,11A1,' to',11A1,' eV'//
     2 '   Material  MAT  MT Average'/)
   50 FORMAT(
     1 ' Unshielded ResonAnce Integrals'/
     1 1X,11A1,1X,A4,A2/
     3 1X,11A1,' to',11A1,' eV'//
     2 '   Material  MAT  MT R.I.'/)
      END
      SUBROUTINE XYPAGE(I,ITYPE,X,Y)
C=======================================================================
C
C     RETRIEVE AND RETURN SELECTED X AND Y VALUES FROM THE PAGING
C     SYSTEM.
C
C     I     = DATA POINT INDEX (I=1 UP TO TABLE SIZE)
C     ITYPE = PAGING SYSTEM ARRAY INDEX WHICH DEFINES TYPE OF DATA.
C           = 1 ENERGY DEPENDENT WEIGHTING SPECTRUM.
C           = 2 TOTAL CROSS SECTION.
C           = 3 ELASTIC CROSS SECTION.
C           = 4 CAPTURE CROSS SECTION
C           = 5 FISSION CROSS SECTION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/PAGER/NPAGE,NPAGM1
      COMMON/PAGERg/NPAGE,NPAGM1
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
      COMMON/INPAGE/IXYLOW(5),IXYHI(5),ISCR(5)
      INCLUDE 'groupie.h'
C-----INSURE POINT INDEX IS IN LEGAL RANGE.
      IF(I.GT.0.AND.I.LE.NPTAB(ITYPE)) GO TO 10
C-----ILLEGAL POINT INDEX.
      WRITE(OUTP,40) I,NPTAB(ITYPE)
      X=0.0
      Y=0.0
      RETURN
C-----IF DATA IS NOT IN CORE LOAD CORRECT PAGE.
   10 IF(I.GT.IXYLOW(ITYPE)) GO TO 30
      IXYHI(ITYPE)=0
      NSCR=ISCR(ITYPE)
      REWIND NSCR
   20 IXYLOW(ITYPE)=IXYHI(ITYPE)
      IXYHI(ITYPE)=IXYHI(ITYPE)+NPAGE
      CALL IBLOCK(ISCR(ITYPE),XPAGE(1,ITYPE),YPAGE(1,ITYPE),NPAGE)
   30 IF(I.GT.IXYHI(ITYPE)) GO TO 20
C-----DEFINE REQUIRED POINT.
      ICORE=I-IXYLOW(ITYPE)
      X=XPAGE(ICORE,ITYPE)
      Y=YPAGE(ICORE,ITYPE)
      RETURN
   40 FORMAT(15H XYPAGE..Index=,I6,' (MUST be 1 TO',I6,')')
      END
      SUBROUTINE BANDIT
C=======================================================================
C
C     DEFINE MULTI-BAND WEIGHTS AND CROSS SECTIONS FOR TOTAL, ELASTIC,
C     CAPTURE AND FISSION. MULTI-BAND PARAMETERS ARE SELECTED TO
C     CONSERVE MOMENTS OF THE CROSS SECTIONS BASED UPON BONDERENKO
C     SELF-SHIELDED CROSS SECTIONS. IN ALL CASES THE AVERAGE SHIELDED
C     AND UNSHIELDED CROSS SECTIONS WILL BE CONSERVED. ADDITIONAL
C     MOMENTS WILL BE CONSERVED AS MORE CROSS SECTION BANDS ARE
C     USED.
C
C     THIS ROUTINE WILL START WITH ONE BAND AND INCREASE THE NUMBER OF
C     BANDS UNTIL EITHER ALL BONDERENKO SELF-SHIELDED CROSS SECTIONS
C     CAN BE APPROXIMATED TO WITHIN SOME ACCURACY (ERROK) OR THE MAXIMUM
C     DESIRED NUMBER OF BAND (NBMAX) IS REACHED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,OTAPE2,TYPSIG
      CHARACTER*1 FIELDX,ZABCD
      CHARACTER*4 BUMMER
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
      COMMON/ELPAS2/EGB,TMPK
CAK   COMMON/FILLER/ISECT
      COMMON/FILLERg/ISECT
      COMMON/ZABAD/IZABAD
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      COMMON/ELPASZ/ZABCD(10)
      COMMON/LOGCOM/MYLOG(8)
      INCLUDE 'groupie.h'
      DIMENSION XCB(5),WTB(5)
      DATA IERR/0/
C-----INITIALIZE ALL PARAMETERS (ONLY USED FOR 3 OR MORE BANDS).
      DATA BUMMER/'ERR '/
      DATA ZERO   /0.00D+00/
      DATA QUARTER/0.25D+00/
      DATA HALF   /0.50D+00/
      DATA ONE    /1.00D+00/
      DATA TWO    /2.00D+00/
C
C     INITIALIZE TO 2 BANDS BEING USED
C
      NBUSE=2
      NBNEED=NBUSE
      NBAND=NBUSE
C-----INITIALIZE ALL BAND VALUES FOR GROUP
      DO I=1,6
      DO J=1,5
      XCBAND(I,J,IGR) = 0.0D+0
      ENDDO
      ENDDO
C
C     EXTEND CROSS SECTION LIMITS OUTWARD
C     DEFINE 3 RANGES FOR SUCCESSIVE PASSES
C     #1:   1%
C     #2:  10%
C     #3: 100%
C
      DO K=2,NSECT
      YLOWP1(K)  = 0.99*YLOW(K)
      YLOWP2(K)  = 0.90*YLOW(K)
      YLOWP3(K)  = 0.50*YLOW(K)
      YLOW (K)   = 0.50*YLOW(K)
      YHIGHP1(K) = 1.01*YHIGH(K)
      YHIGHP2(K) = 1.10*YHIGH(K)
      YHIGHP3(K) = 2.00*YHIGH(K)
      YHIGH(K)   = 2.00*YHIGH(K)
      ENDDO
C-----------------------------------------------------------------------
C
C     IF TOTAL = 0 OR NO SHIELDING USE 1 BAND
C
C-----------------------------------------------------------------------
      SIGT1   = XCINT( 1,2)
      SIGT2   = XCINT(23,2)
      IF(SIGT1.LE.ZERO)  GO TO 240
      IF(SIGT2.ge.SIGT1) GO TO 240
C-----------------------------------------------------------------------
C
C     LITTLE SHIELDING = CONSERVE 2 MOMENTS.
C
C-----------------------------------------------------------------------
c-----SMALL 0.01% TEST
      IF(SIGT2.gt.0.9999 *SIGT1) THEN
      MYLOGX = 2
      ABAND  = ONE/SIGT2
      BB2    = (SIGT1 - SIGT2)/SIGT1
      if(BB2.le.0.0d0) go to 240
      BBAND  = DSQRT(BB2)/SIGT2
      DBAND  = 0.0d0                          ! for below tests
      WT1    = HALF
      WT2    = HALF
      BANDT1 = ONE/(ABAND + BBAND)
      BANDT2 = ONE/(ABAND - BBAND)
      go to 10
      endif
C-----------------------------------------------------------------------
C
C     USE 2 BANDS
C
C     METHODB #1: CONSERVE 1/TOT AND 1/TOT**2
C     METHODB #2: CONSERVE 1/TOT AND 1/(TOT + <TOT>)
C
C-----------------------------------------------------------------------
      MYLOGX = 3
      if(METHODB.eq.1) then
      SIGT3 = XCINT(24,2)
      ABAND = (SIGT1 - SIGT3)/(TWO*(SIGT1 - SIGT2)*SIGT3)
      else
      SIGT3 = XCINT(12,2)
      ABAND  = (SIGT1 - SIGT2)/(TWO*(SIGT1 - SIGT3)*SIGT2)
      endif
C
C     THE FOLLOWING IS THE SAME FOR BOTH METHODS
C
      BB2   = (SIGT2*ABAND*(SIGT1*ABAND-TWO)+ONE)/(SIGT1*SIGT2)
      if(BB2.le.0.0d0) go to 240
      BBAND = DSQRT(BB2)
      DBAND = (ONE-ABAND*SIGT2)/(TWO*SIGT2*BBAND) ! Note: divide by B
      WT1   = HALF + DBAND
      WT2   = HALF - DBAND
      BANDT1 = ONE/(ABAND + BBAND)
      BANDT2 = ONE/(ABAND - BBAND)
C-----------------------------------------------------------------------
C
C     TEST MOMENTS
C
C-----------------------------------------------------------------------
   10 APPROX1 = WT1*BANDT1 + WT2*BANDT2
      APPROX2 = ONE/(WT1/BANDT1 + WT2/BANDT2)
      if(MYLOGX.eq.2) then
      SIGT3   = XCINT(12,2)                   ! for below tests
      APPROX3 = SIGT3
      else                ! for below tests
      if(METHODB.eq.1) then
      APPROX3 = (WT1/BANDT1+WT2/BANDT2)/(WT1/BANDT1**2+WT2/BANDT2**2)
      else
      BB1     = BANDT1+SIGT1
      BB2     = BANDT2+SIGT1
      APPROX3 = (WT1*BANDT1/BB1 + WT2*BANDT2/BB2)/
     1          (WT1       /BB1 + WT2       /BB2)
      endif
      endif
      ERR1 = (SIGT1 - APPROX1)/SIGT1
      ERR2 = (SIGT2 - APPROX2)/SIGT2
      ERR3 = (SIGT3 - APPROX3)/SIGT3
      IF(DABS(ERR1).GT.1.0D-6.OR.
     1   DABS(ERR2).GT.1.0D-6.OR.
     1   DABS(ERR3).GT.1.0D-6) THEN
c-----2017/4/14 - ERROR message turned off
      WRITE(3,20) IGR,SIGT1,     SIGT2,     SIGT3,   METHODB,
     1                APPROX1,   APPROX2,   APPROX3,
     1                100.0*ERR1,100.0*ERR2,100*ERR3
      WRITE(*,20) IGR,SIGT1,     SIGT2,     SIGT3,   METHODB,
     1                APPROX1,   APPROX2,   APPROX3,
     1                100.0*ERR1,100.0*ERR2,100*ERR3
   20 FORMAT(I6,1P3D12.5,' moments   Method #',i1/
     1       6X,1P3D12.5,' approximations'/
     1       6X,1P3D12.5,' % differences')
      CALL ENDERROR
c-----2017/4/14 - ERROR message turned off
      ENDIF
C-----TEST D, A AND B
      IF(DABS(DBAND).GT.HALF.OR.BBAND.GT.ABAND) THEN
      WRITE(3,30) IGR,DBAND,ABAND,BBAND,METHODB
      WRITE(*,30) IGR,DBAND,ABAND,BBAND,METHODB
   30 FORMAT(' ERROR - IGR/D/A/B=',I6,1P3D20.12,' Method #',i1)
      CALL ENDERROR
      ENDIF
C-----------------------------------------------------------------------
C
C     DEFINE TOTAL BAND PARAMETERS
C
C-----------------------------------------------------------------------
      WTBAND(1,1,IGR) = WT1
      WTBAND(1,2,IGR) = WT2
      XCBAND(2,1,IGR) = BANDT1
      XCBAND(2,2,IGR) = BANDT2
C-----------------------------------------------------------------------
C
C     NEW LINEAR EQUATIONS.
C
C     SIG1 = <SIG1> - p/P1  BOTH TOTAL AND PARTIALS
c     SIG2 = <SIG1> + p/P2
C
C     T^2 = [(1/4-D^2)*<SIGT1>+2DT][<SIGT1>-<SIGT2>]
C     T^2 + 2bT + c = 0
C     Switching signs of both b and c
C     b = D[<sigt1>-<sigt2>]                > 0
C     c = (1/4-D^2)<sigt1>[<sigt1>-<sigt2>] > 0
C     T = b+[b^2 + c]^1/2  Only the + sign has a positive solution
C
C     P   = T[<SIGP1>-<SIGP2>]/[<SIGT1>-<SIGT2>]
C
C-----------------------------------------------------------------------
      DSIGT12 = SIGT1 - SIGT2
      if(DSIGT12.eq.0.0d0) go to 240            ! No Shielding
      BTRY = DBAND*DSIGT12
      CTRY = (QUARTER - DBAND**2)*SIGT1*DSIGT12
      ROOT = DSQRT(BTRY**2 + CTRY)
      DTOT = BTRY+ROOT
C-----------------------------------------------------------------------
C
C     LOOP OVER PARTIALS
C
C-----------------------------------------------------------------------
C
C     FIRST TRY STANDARD
C
      TOTNORM = DTOT/DSIGT12
      DO 70 MSECT=3,6
      SIGP1 = XCINT( 1,MSECT)
      SIGP2 = XCINT(23,MSECT)
C-----NOTHING TO DO IF NO CROSS SECTION.
      IF(SIGP1.LE.ZERO) GO TO 50
      CBAND  = TOTNORM*(SIGP1 - SIGP2)
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
C
C     TEST MOMENTS
C
      SIXP1  = WT1*BANDP1 + WT2*BANDP2
      SIXP2  = (WT1*BANDP1/BANDT1 + WT2*BANDP2/BANDT2)/
     1         (WT1       /BANDT1 + WT2       /BANDT2)
      ERR1 = DABS(SIGP1 - SIXP1)
      ERR2 = DABS(SIGP2 - SIXP2)
      if(DABS(SIGP1).ne.0.0d0.and.DABS(SIGP2).ne.0.0d0) then
      if(ERR1.gt.1.0d-6*DABS(SIXP1).or.
     1   ERR2.gt.1.0d-6*DABS(SIXP2)) then
c-----2017/4/14 - ERROR message turned off
      write(3,40) IGR,MSECT,SIGP1,     SIGP2,
     1                      SIXP1,     SIXP2,
     2                      100.0*ERR1,100.0*ERR2
      write(*,40) IGR,MSECT,SIGP1,     SIGP2,
     1                      SIXP1,     SIXP2,
     2                      100.0*ERR1,100.0*ERR2
   40 format(2i4,1p2d12.5,' partials'/
     1        8x,1p2d12.5,' approximations'/
     1        8x,1p2d12.5,' % differences')
c     CALL ENDERROR
c-----2017/4/14 - ERROR message turned off
      endif
      endif
C
C     INSURE PARTIAL CROSS SECTIONS ARE POSITIVE.
C
      if(CBAND.ge.0.0d0) then
      if(CBAND.gt. 0.999*WT1*SIGP1) then
      CBAND  =     0.999*WT1*SIGP1
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
      endif
      else
      if(CBAND.lt.-0.999*WT2*SIGP1) then
      CBAND  =    -0.999*WT2*SIGP1
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
      endif
      endif
C-----SKIP RANGE TEST FOR "OTHER"
      IF(MSECT.EQ.6) GO TO 60
C-----INSURE PARTIALS ARE IN LEGAL RANGE (STRICT HERE).
      IF(BANDP1.LT.YLOW(MSECT).OR.BANDP1.GT.YHIGH(MSECT)) GO TO 80
      IF(BANDP2.LT.YLOW(MSECT).OR.BANDP2.GT.YHIGH(MSECT)) GO TO 80
      GO TO 60
C
C     USE NO SHIELDING FOR PARTIAL.
C
   50 BANDP1=SIGP1
      BANDP2=SIGP1
C
C     PARTIAL PARAMETERS ARE O.K. - SAVE THEM.
C
   60 XCBAND(MSECT,1,IGR)=BANDP1
      XCBAND(MSECT,2,IGR)=BANDP2
C-----END OF PARTIAL LOOP
   70 CONTINUE
C-----ALL PARTIALS ARE O.K.
      MYLOG(MYLOGX) = MYLOG(MYLOGX) + 1
      GO TO 260
C-----------------------------------------------------------------------
C
C     NEXT TRY ITERATION
C
C-----------------------------------------------------------------------
   80 DBSTART = DBAND
      DBSTEP = 1.0D-05
      IF(DBAND.GT.0.0D+0) DBSTEP = -DBSTEP
C-----LOOP FOR STRICT/LOOSE CONVERGENCE
C     MYLOOP = 1:   1%
C              2:  10%
C              3: 100%
C              4: POSITIVE
      DO 200 MYLOOP=1,4
      MYCHANGE = 0
C-----LOOP TO VARY DBAND
      DO 190 IDB = -50000,50000
      DBAND = IDB*DBSTEP
      BTRY   = DBAND*DSIGT12
      CTRY   = (QUARTER - DBAND**2)*SIGT1*DSIGT12
      ROOT   = DSQRT(BTRY**2 + CTRY)
      DTOT   = BTRY+ROOT
      WT1    = HALF + DBAND
      WT2    = HALF - DBAND
      BANDT1 = SIGT1 - DTOT/WT1
      BANDT2 = SIGT1 + DTOT/WT2
C-----SKIP IF TOTAL PARAMETERS AREE NOT O.K.
      IF(WT1.LT.ZERO.OR.WT2.LT.ZERO) GO TO 190
      IF(BANDT1.LT.YLOW(2).OR.BANDT1.GT.YHIGH(2)) GO TO 190
      IF(BANDT2.LT.YLOW(2).OR.BANDT2.GT.YHIGH(2)) GO TO 190
      TOTNORM = DTOT/DSIGT12
      DO 180 MSECT=3,6
      SIGP1 = XCINT( 1,MSECT)
      SIGP2 = XCINT(23,MSECT)
C-----NOTHING TO DO IF NO CROSS SECTION.
      IF(SIGP1.LE.ZERO) GO TO 160
      CBAND = TOTNORM*(SIGP1 - SIGP2)
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
C-----SKIP RANGE TEST FOR "OTHER"
      IF(MSECT.EQ.6) GO TO 170
C-----INSURE PARTIALS ARE IN LEGAL RANGE
      GO TO (90,110,130,150), MYLOOP
C-----PASS 1: STRICT 1%
   90 IF(BANDP1.GE.YLOWP1(MSECT).AND.BANDP1.LE.YHIGHP1(MSECT)) GO TO 100
      GO TO 190
  100 IF(BANDP2.GE.YLOWP1(MSECT).AND.BANDP2.LE.YHIGHP1(MSECT)) GO TO 170
      GO TO 190
C-----PASS 2: NOT STRICT 10%
  110 IF(BANDP1.GE.YLOWP2(MSECT).AND.BANDP1.LE.YHIGHP2(MSECT)) GO TO 120
      GO TO 190
  120 IF(BANDP2.GE.YLOWP2(MSECT).AND.BANDP2.LE.YHIGHP2(MSECT)) GO TO 170
      GO TO 190
C-----PASS 3: SOFT 100%
  130 IF(BANDP1.GE.YLOWP3(MSECT).AND.BANDP1.LE.YHIGHP3(MSECT)) GO TO 140
      GO TO 190
  140 IF(BANDP2.GE.YLOWP3(MSECT).AND.BANDP2.LE.YHIGHP3(MSECT)) GO TO 170
      GO TO 190
C-----PASS 4: VERY SOFT
  150 IF(BANDP1.GE.ZERO.AND.BANDP2.GE.ZERO) GO TO 170
      GO TO 190
C
C     USE NO SHIELDING FOR PARTIAL.
C
  160 BANDP1=SIGP1
      BANDP2=SIGP1
C
C     PARTIAL PARAMETERS ARE O.K. - SAVE THEM.
C
  170 XCBAND(MSECT,1,IGR)=BANDP1
      XCBAND(MSECT,2,IGR)=BANDP2
C-----END OF PARTIAL LOOP
  180 CONTINUE
C-----ALL PARTIALS ARE O.K. - IF CHANGED, COPY NEW TOTALS
      MYLOGX = 3+MYLOOP
      MYLOG(MYLOGX) = MYLOG(MYLOGX) + 1
      IF(MYCHANGE.NE.0) THEN
      WTBAND(1,1,IGR) = WT1
      WTBAND(1,2,IGR) = WT2
      XCBAND(2,1,IGR) = BANDT1
      XCBAND(2,2,IGR) = BANDT2
      ENDIF
      GO TO 260
C-----END OF LOOP TO CHANGE DBAND
  190 MYCHANGE = 1
C-----END OF STRICT/LOOSE LOOP
  200 CONTINUE
C-----------------------------------------------------------------------
C
C     CHANGING DBAND DID NOT HELP - RESET TO ORIGINAL
C     AND INSURE PARTIALS ARE NOT NEGATIVE.
C
C-----------------------------------------------------------------------
      MYLOGX   = 8
      MYLOG(8) = MYLOG(8) + 1
      DBAND = DBSTART
      WT1   = HALF + DBAND
      WT2   = HALF - DBAND
      BTRY = DBAND*DSIGT12
      CTRY = (QUARTER - DBAND**2)*SIGT1*DSIGT12
      ROOT = DSQRT(BTRY**2 + CTRY)
      DTOT= BTRY+ROOT
      BANDT1 = SIGT1 - DTOT/WT1
      BANDT2 = SIGT1 + DTOT/WT2
      WTBAND(1,1,IGR) = WT1
      WTBAND(1,2,IGR) = WT2
      XCBAND(2,1,IGR) = BANDT1
      XCBAND(2,2,IGR) = BANDT2
      TOTNORM = DTOT/DSIGT12
C-----LOOP OVER REACTIONS
      DO 230 MSECT=3,6
      SIGP1 = XCINT( 1,MSECT)
      SIGP2 = XCINT(23,MSECT)
C-----NOTHING TO DO IF NO CROSS SECTION.
      IF(SIGP1.LE.ZERO) GO TO 210
      CBAND = TOTNORM*(SIGP1 - SIGP2)
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
C-----NO NEGATIVES
      IF(BANDP1.LT.ZERO) THEN
      CBAND =  WT1*SIGP1
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
      ENDIF
      IF(BANDP2.LT.ZERO) THEN
      CBAND = -WT2*SIGP1
      BANDP1 = SIGP1 - CBAND/WT1
      BANDP2 = SIGP1 + CBAND/WT2
      ENDIF
      IF(BANDP1.LT.ZERO) BANDP1 = ZERO
      IF(BANDP2.LT.ZERO) BANDP2 = ZERO
      GO TO 220
  210 BANDP1=SIGP1
      BANDP2=SIGP1
  220 XCBAND(MSECT,1,IGR)=BANDP1
      XCBAND(MSECT,2,IGR)=BANDP2
  230 CONTINUE
      GO TO 260
C-----------------------------------------------------------------------
C
C     NO SHIELDING - USE 1 BAND
C
C-----------------------------------------------------------------------
  240 MYLOGX   = 1
      MYLOG(1) = MYLOG(1) + 1
      WTBAND(1,1,IGR)=0.5d0
      WTBAND(1,2,IGR)=0.5d0
      DO 250 MSECT=2,6
      XCBAND(MSECT,1,IGR)=XCINT(1,MSECT)
      XCBAND(MSECT,2,IGR)=XCINT(1,MSECT)
  250 CONTINUE
C-----------------------------------------------------------------------
C
C     DEFINE TOTAL = SUM OF PARTS
C
C     WARNING - This SUM allows for competition - MSECT = 6
C
C-----------------------------------------------------------------------
  260 XCBAND(2,1,IGR)=ZERO
      XCBAND(2,2,IGR)=ZERO
      DO 270 MSECT=3,6
      XCBAND(2,1,IGR)=XCBAND(2,1,IGR)+XCBAND(MSECT,1,IGR)
  270 XCBAND(2,2,IGR)=XCBAND(2,2,IGR)+XCBAND(MSECT,2,IGR)
C-----------------------------------------------------------------------
C
C     DEFINE ERROR IN FIT
C
C-----------------------------------------------------------------------
      ERRMAX=ZERO
      XCB(1)=XCBAND(2,1,IGR)/XCINT(1,2)
      XCB(2)=XCBAND(2,2,IGR)/XCINT(1,2)
      WTB(1)=WTBAND(1,1,IGR)
      WTB(2)=WTBAND(1,2,IGR)
      DO 290 I=2,23
      XCTOP=0.0
      AVNORM(I)=ZERO
      DO 280 NB=1,NBUSE
      FACTOR=WTB(NB)/(XCB(NB)+SIGMAB(I))
      XCTOP=XCTOP+FACTOR*XCB(NB)
  280 AVNORM(I)=AVNORM(I)+FACTOR
      IF(AVNORM(I).gt.0.0d0.and.XCFI(I,2).gt.0.0d0) then
      XCAV=XCTOP/AVNORM(I)
      ERNOW(I,NBUSE)=DABS(ONE-XCAV/XCFI(I,2))
      IF(ERNOW(I,NBUSE).GT.ERRMAX) ERRMAX=ERNOW(I,NBUSE)
      else
      ERNOW(I,NBUSE) = 0.0d0
      endif
  290 CONTINUE
      AVNORM(24)=WTB(1)/(XCB(1)**2)
      AVNORM(25)=AVNORM(24)/XCB(1)
      DO 300 NB=2,NBUSE
      FACTOR=WTB(NB)/(XCB(NB)**2)
      AVNORM(24)=AVNORM(24)+FACTOR
  300 AVNORM(25)=AVNORM(25)+FACTOR/XCB(NB)
C-----------------------------------------------------------------------
C
C     MULTI-BAND PARAMETERS ARE PHYSICALLY ACCEPTABLE AND YIELD LOWEST
C     ERROR OF ANY NUMBER OF BANDS CONSIDERED SO FAR.
C
C-----------------------------------------------------------------------
      DO 310 I=2,23
      IF(ERNOW(I,NBAND).GT.ERMAT(I,NBAND)) ERMAT(I,NBAND)=ERNOW(I,NBAND)
  310 CONTINUE
C-----------------------------------------------------------------------
C
C     ITERATION IS COMPLETED. SET ALL REMAINING WEIGHTS TO ZERO AND
C     CROSS SECTIONS TO UNSHIELDED AVERAGES. SAVE LAST ACCEPTABLE
C     SET OF WEIGHTS AND CROSS SECTIONS FOR DETERMNATION OF PARTIAL
C     CROSS SECTIN BAND PARAMETERS.
C
C-----------------------------------------------------------------------
      IF(NBAND.GE.NBMAX) GO TO 330
      NBP1=NBAND+1
      DO 320 NB=NBP1,NBMAX
      WTBAND(1,NB,IGR)=ZERO
      DO 320 MSECT=2,6
  320 XCBAND(MSECT,NB,IGR)=XCINT(1,MSECT)
C-----------------------------------------------------------------------
C
C     CHECK THAT TOTAL IS SUM OF PARTS IN EACH BAND
C
C-----------------------------------------------------------------------
  330 IF(XCINT(1,2).LE.ZERO) GO TO 610
      DO 370 NB=1,NBAND
      SUMMER=0.0
      IF(XCBAND(2,NB,IGR).GT.ZERO) THEN
      DO 340 KSECT=3,6
  340 SUMMER=SUMMER+XCBAND(KSECT,NB,IGR)
      DIFF=(XCBAND(2,NB,IGR)-SUMMER)/XCBAND(2,NB,IGR)
      IF(DABS(DIFF).LE.1.0D-03) GO TO 370
      GO TO 350
      ELSE
      DIFF   = ZERO
      SUMMER = ZERO
      ENDIF
  350 WRITE(OUTP,360) IGR,NB,WTBAND(1,NB,IGR),100.0*DIFF,SUMMER,
     1 (XCBAND(K,NB,IGR),K=2,6)
      WRITE(   *,360) IGR,NB,WTBAND(1,NB,IGR),100.0*DIFF,SUMMER,
     1 (XCBAND(K,NB,IGR),K=2,6)
  360 FORMAT(' Band Sum ERROR',I6,I2,1PE11.4,0PF10.2,' %',1P8E11.4)
  370 CONTINUE
C-----------------------------------------------------------------------
C
C     COMPARE EXACT CROSS SECTION MOMENTS TO THOSE RECONSTRUCTED
C     FROM THE BAND PARAMETERS.
C
C-----------------------------------------------------------------------
C-----DEFINE NORMALIZATION (DENOMINATOR) FOR ALL SELF-SHIELDED
C-----CROSS SECTIONS.
      DO 380 I=2,25
  380 AVNORM(I)=ZERO
      DO 400 NB=1,NBAND
      DO 390 I=2,23
  390 AVNORM(I)=AVNORM(I)+WTBAND(1,NB,IGR)/
     1 (XCBAND(2,NB,IGR)+SHIELD(I))
      SA2=XCBAND(2,NB,IGR)**2
      SA3=SA2*XCBAND(2,NB,IGR)
      AVNORM(24)=AVNORM(24)+WTBAND(1,NB,IGR)/SA2
  400 AVNORM(25)=AVNORM(25)+WTBAND(1,NB,IGR)/SA3
C-----DEFINE NUMERATOR FOR EACH REACTION AND EACH SELF-SHIELDED
C-----CROSS SECTION.
      DO 600 ISECT=2,6
      IF(XCINT(1,ISECT).LE.ZERO) GO TO 600
      DO 410 I=1,25
  410 AVEXP(I)=0.0
      DO 430 NB=1,NBAND
      ADDNB=WTBAND(1,NB,IGR)*XCBAND(ISECT,NB,IGR)
      AVEXP(1)=AVEXP(1)+ADDNB
      DO 420 I=2,23
  420 AVEXP(I)=AVEXP(I)+ADDNB/(XCBAND(2,NB,IGR)+SHIELD(I))
      SA2=XCBAND(2,NB,IGR)**2
      SA3=SA2*XCBAND(2,NB,IGR)
      AVEXP(24)=AVEXP(24)+ADDNB/SA2
  430 AVEXP(25)=AVEXP(25)+ADDNB/SA3
C-----DEFINE NORMALIZED RECONSTRUCTED SELF-SHIELDED
C-----CROSS SECTIONS.
      DO 440 I=2,25
  440 AVEXP(I)=AVEXP(I)/AVNORM(I)
C-----DEFINE DIFFERENCE BETWEEN EXACT AND RECONSTRUCTED VALUES.
      DO 460 I=1,25
      IF(XCINT(I,ISECT).GT.0.0) GO TO 450
      ERNOW(I,NBAND)=ZERO
      GO TO 460
  450 ERNOW(I,NBAND)=DABS(ONE-AVEXP(I)/XCINT(I,ISECT))
  460 CONTINUE
C-----UPDATE MAXIMUM ERROR FOR SELF-SHIELDED TOTAL CROSS SECTIONS.
      IF(ISECT.NE.2) GO TO 480
      IF(ERNOW(24,NBAND).GT.ERMAT(24,NBAND))
     1 ERMAT(24,NBAND)=ERNOW(24,NBAND)
      IF(ERNOW(25,NBAND).GT.ERMAT(25,NBAND))
     1 ERMAT(25,NBAND)=ERNOW(25,NBAND)
      DO 470 I=1,25
C-----UPDATE MAXIMUM ERROR IF MULTI-BANDS ARE NOT USED
C-----(I.E. ONE BAND/GROUP)
      IF(XCFI(I,2).gt.ZERO) then
      ERNOW(I,1)=DABS(ONE-ONE/XCFI(I,2))
      IF(ERNOW(I,1).GT.ERMAT(I,1)) ERMAT(I,1)=ERNOW(I,1)
      else
      ERNOW(I,1) = 0.0d0
      endif
  470 CONTINUE
C-----------------------------------------------------------------------
C
C     CHECK DATA AND PRINT IF ANY ERRORS
C
C-----------------------------------------------------------------------
C-----CHECK WEIGHTS BETWEEN ZERO AND ONE AND BAND CROSS SECTIONS
C-----BETWEEN MINIMUM AND MAXIMUM FOR EACH REACTION.
  480 DO 490 NB=1,NBMAX
      IF(WTBAND(1,NB,IGR).LT.ZERO.OR.WTBAND(1,NB,IGR).GT.ONE)
     1 GO TO 500
C-----NO LIMIT TEST FOR OTHER
      IF(ISECT.GE.5) GO TO 490
C-----NO LIMIT TEST FOR CROSS SECTION SMALL COMPARED TO TOTAL
      IF(XCINT(1,ISECT).LE.0.005*XCINT(1,2)) GO TO 490
      IF(XCBAND(ISECT,NB,IGR).LT.YLOW (ISECT).OR.
     1   XCBAND(ISECT,NB,IGR).GT.YHIGH(ISECT)) GO TO 500
  490 CONTINUE
C-----FOR EACH CROSS SECTION CHECK TWO CONSERVED MOMENTS.
      IF(ERNOW(1,NBAND).LE.ER10PC.AND.ERNOW(23,NBAND).LE.ERROK)
     1 GO TO 600
C-----------------------------------------------------------------------
C
C     ONE OR MORE ERRORS IN PARAMETERS. RE-CHECK, FLAG AND LIST ALL
C     PARAMETERS.
C
C-----------------------------------------------------------------------
  500 XES=XE
C-----NO ERROR MESSAGE FOR OTHER
      IF(ISECT.GT.5) GO TO 600
      WRITE(OUTP,790) IGR,ZABCD,REACT2(1,ISECT),REACT2(2,ISECT),
     1 EGB,XES,YLOW(2),YHIGH(2),YLOW(ISECT),YHIGH(ISECT),XCINT(1,2),
     2             XCINT(23,2),XCINT(1,ISECT),
     3 XCINT(23,ISECT),AVEXP(1),          AVEXP(23),
     4 100.0*ERNOW(1,NBAND),
     5 100.0*ERNOW(23,NBAND)
C-----CHECK ERRORS.
      IERR=0
      DO 510 I=1,3
  510 POINT(I)='    '
      IF(ERNOW(1,NBAND).LE.ER10PC) GO TO 520
      IERR=1
      POINT(1)=BUMMER
  520 IF(ERNOW(23,NBAND).LE.ERROK) GO TO 530
      IERR=1
      POINT(2)=BUMMER
  530 IF(IERR.GT.0) WRITE(OUTP,830) (POINT(I),I=1,2)
C-----CHECK BAND WEIGHTS.
      WRITE(OUTP,800) (WTBAND(1,NB,IGR),NB=1,NBMAX)
      IF(ISECT.GT.2) GO TO 560
      IERR=0
      DO 550 NB=1,NBMAX
      IF(WTBAND(1,NB,IGR).LT.ZERO.OR.WTBAND(1,NB,IGR).GT.ONE)
     1 GO TO 540
      POINT(NB)='    '
      GO TO 550
  540 IERR=1
      POINT(NB)=BUMMER
      IF(IZABAD.LE.0) WRITE(OUTP,780) IZA,IGR
      IZABAD=1
  550 CONTINUE
      IF(IERR.GT.0) WRITE(OUTP,830) (POINT(NB),NB=1,NBMAX)
  560 WRITE(OUTP,810) (XCBAND(2,NB,IGR),NB=1,NBMAX)
C-----CHECK PARTIAL BAND CROSS SECTIONS.
      WRITE(OUTP,820) (XCBAND(ISECT,NB,IGR),NB=1,NBMAX)
      IERR=0
      DO 580 NB=1,NBMAX
      IF(XCBAND(ISECT,NB,IGR).LT.YLOW (ISECT).OR.
     1   XCBAND(ISECT,NB,IGR).GT.YHIGH(ISECT)) GO TO 570
      POINT(NB)='    '
      GO TO 580
  570 IERR=1
      POINT(NB)=BUMMER
  580 CONTINUE
      IF(IERR.GT.0) WRITE(OUTP,830) (POINT(NB),NB=1,NBMAX)
      WRITE(OUTP,840)
      IF(XCBAND(ISECT,1,IGR).LT.0.0.OR.
     1   XCBAND(ISECT,2,IGR).LT.0.0) THEN
      WRITE(3,590) ISECT,IGR,(XCBAND(ISECT,k,IGR),k=1,2)
      WRITE(*,590) ISECT,IGR,(XCBAND(ISECT,k,IGR),k=1,2)
  590 FORMAT(' ERROR - Negative Partials..ISECT/IGR=',2i5/
     1       '         Partialls=',1P2D12.5)
      CALL ENDERROR
      ENDIF
C-----------------------------------------------------------------------
C
C     END OF CHECK LOOP
C
C-----------------------------------------------------------------------
  600 CONTINUE
C-----------------------------------------------------------------------
C
C     LIST MULTI-BAND PARAMETERS.
C
C-----------------------------------------------------------------------
  610 IF(LIST2.LE.0) GO TO 770
      LINOUT=1
      DO 620 NB=2,NBAND
      IF(WTBAND(1,NB,IGR).GT.ZERO) LINOUT=LINOUT+1
  620 CONTINUE
C-----CONVERT BAND CROSS SECTIONS TO HOLLERITH. CONVERT UNUSED TO BLANK.
      DO 760 NB=1,LINOUT
      DO 650 ISECT=2,5
      IF(ISECT.GT.NSECT) GO TO 630
      CALL OUT9G(XCBAND(ISECT,NB,IGR),FIELDX(1,ISECT))
      GO TO 650
C-----UNUSED =  BLANK
  630 DO 640 M=1,11
  640 FIELDX(M,ISECT)=' '
  650 CONTINUE
      IF(NB.GT.1) GO TO 740
C-----CONVERT UNSHIELDED AVERAGES TO HOLLERITH.
      DO 680 ISECT=2,5
      ISP4=ISECT+4
      IF(ISECT.GT.NSECT) GO TO 660
      CALL OUT9G(XCINT(1,ISECT),FIELDX(1,ISP4))
      GO TO 680
C-----UNUSED = BLANK
  660 DO 670 M=1,11
  670 FIELDX(M,ISP4)=' '
  680 CONTINUE
C-----LIST REMAINDER WITH FIRST BAND IF REMAINDER IS AT LEAST
C-----0.05 PER-CENT OF THE UNSHIELDED TOTAL CROSS SECTION.
      REMAIN=XCINT(1,2)
      DO 690 ISECT=3,NSECT
  690 REMAIN=REMAIN-XCINT(1,ISECT)
      IF(REMAIN.LT.0.0005*XCINT(1,2)) GO TO 700
      CALL OUT9G(REMAIN,FIELDX(1,10))
      GO TO 720
C-----UNUSED = BLANK
  700 DO 710 M=1,11
  710 FIELDX(M,10)=' '
C-----LIST PARAMETERS FOR FIRST BAND AND UNSHIELDED GROUP
C-----AVERAGES.
  720 IF(NSECT.LT.5) GO TO 730
      WRITE(LIST2,850) IGR,(FIELDX(M,1),M=1,11),NB,
     1 WTBAND(1,NB,IGR),
     2 ((FIELDX(M,I),M=1,11),I=2,10)
      GO TO 760
  730 WRITE(LIST2,870) IGR,(FIELDX(M,1),M=1,11),NB,
     1 WTBAND(1,NB,IGR),
     2 ((FIELDX(M,I),M=1,11),I=2,4),
     3 ((FIELDX(M,I),M=1,11),I=6,8),
     4 (FIELDX(M,10),M=1,11)
      GO TO 760
C-----LIST PARAMETERS FOR SECOND TO N-TH BANDS.
  740 IF(NSECT.LT.5) GO TO 750
      WRITE(LIST2,860) NB,WTBAND(1,NB,IGR),
     1 ((FIELDX(M,I),M=1,11),I=2,5)
      GO TO 760
  750 WRITE(LIST2,880) NB,WTBAND(1,NB,IGR),
     1 ((FIELDX(M,I),M=1,11),I=2,4)
  760 CONTINUE
  770 RETURN
  780 FORMAT(' ERROR - Bad Weights in',I6,' Group',I6)
  790 FORMAT(1X,78('-')/
     1  1X,'ERROR in Group',I5,1X,10A1,2A4/
     2  1X,'Energy Range    ',1P2E11.4,' eV'/
     3  1X,'Total Limits    ',1P2E11.4/
     4  1X,'Partial Limits  ',1P2E11.4/
     5  1X,'Total Averages  ',1P2E11.4/
     6  1X,'Partial Averages',1P2E11.4/
     7  1X,'Approximation   ',1P2E11.4/
     8  1X,'Error           ',0P2F11.2,' %')
  800 FORMAT(1X,'Weights         ',1P5E11.4)
  810 FORMAT(1X,'Total Bands     ',1P5E11.4)
  820 FORMAT(1X,'Partial Bands   ',1P5E11.4)
  830 FORMAT(17X,5(6X,A4))
  840 FORMAT(1X,78('-'))
  850 FORMAT(I5,11A1,    I5,F8.5,44A1,' *',55A1)
  860 FORMAT(16X,        I5,F8.5,44A1,' *')
  870 FORMAT(11X,I5,11A1,I5,F8.5,33A1,' *',44A1)
  880 FORMAT(27X,        I5,F8.5,33A1,' *')
      END
      SUBROUTINE GROPE(IWANT,EGROUP,NGROUP)
C=======================================================================
C
C     THIS ROUTINE IS DESIGNED TO DEFINE ONE OF THE BUILT IN GROUP
C     STRUCTURES.
C
C     ENERGIES ARE RETURNED IN ASCENDING ORDER IN EV.
C
C     ARGUMENTS
C     ---------
C     IWANT   = ENERGY GROUP STRUCTURE SELECTOR.
C             = 0  175 GROUPS (TART)
C             = 1   50 GROUPS (ORNL)
C             = 2  126 GROUPS (ORNL)
C             = 3  171 GROUPS (ORNL)
C             = 4  620 GROUPS (SAND-II, 1.0D-4, UP TO 18 MEV)
C             = 5  640 GROUPS (SAND-II, 1.0D-4, UP TO 20 MEV)
C             = 6   69 GROUPS (WIMS)
C             = 7   68 GROUPS (GAM-I)
C             = 8   99 GROUPS (GAM-II)
C             = 9   54 GROUPS (MUFT)
C             =10   28 GROUPS (ABBN)
C             =11  616 GROUPS (TART TO 20 MEV)
C             =12  700 GROUPS (TART TO 1 GEV)
C             =13  665 GROUPS (SAND-II, 1.0D-5 eV, UP TO 18 MEV)
C             =14  685 GROUPS (SAND-II, 1.0D-5 eV, UP TO 20 MEV)
C             =15  666 GROUPS (TART TO 200 MEV)
C             =16  725 GROUPS (SAND-II, 1.0D-5 eV, UP TO  60 MEV)
C             =17  755 GROUPS (SAND-II, 1.0D-5 eV, UP TO 150 MEV)
C             =18  765 GROUPS (SAND-II, 1.0D-5 eV, UP TO 200 MEV)
C             =19 1102 GROUPS (UKAEA  , 1.0D-5 eV, UP TO   1 GeV)
C             = OTHERWISE ZERO GROUP BOUNDARIES ARE RETURNED
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION EGROUP(*),
     1 TART(176),ORNLA(51),ORNLB(127),ORNLC(172),
     1 WIMS(70),GAMI(69),GAMII(100),XMUFT(55),ABBN(29),
     2 TART1(40),TART2(40),TART3(40),
     3 TART4(40),TART5(16),GR50A(45),GR50B(6),GR126A(45),GR126B(45),
     4 GR126C(37),GR171A(45),GR171B(45),GR171C(45),GR171D(37),
     7 WIMSA(45),WIMSB(25),
     8 GAMIA(45),GAMIB(24),GAMIIA(45),GAMIIB(45),GAMIIC(10)
      EQUIVALENCE (TART(1),TART1(1)),(TART(41),TART2(1)),
     1 (TART(81),TART3(1)),(TART(121),TART4(1)),(TART(161),TART5(1)),
     2 (ORNLA(1),GR50A(1)),(ORNLA(46),GR50B(1)),(ORNLB(1),GR126A(1)),
     3 (ORNLB(46),GR126B(1)),(ORNLB(91),GR126C(1)),(ORNLC(1),GR171A(1)),
     4 (ORNLC(46),GR171B(1)),(ORNLC(91),GR171C(1)),
     5 (ORNLC(136),GR171D(1))
      EQUIVALENCE (WIMS(1),WIMSA(1)),
     1 (WIMS(46),WIMSB(1)),(GAMI(1),GAMIA(1)),(GAMI(46),GAMIB(1)),
     2 (GAMII(1),GAMIIA(1)),(GAMII(46),GAMIIB(1)),(GAMII(91),GAMIIC(1))
C
C     DEFINE TART 175 GROUP STRUCTURE
C
      DATA TART1/ 1.3068D-03, 5.2271D-03, 2.0908D-02, 3.2669D-02,
     1 4.7044D-02, 8.3215D-02, 1.3068D-01, 1.8817D-01,
     2 2.5613D-01, 3.3453D-01, 4.2339D-01, 5.1230D-01,
     3 7.5270D-01, 1.1761D+00, 1.5106D+00, 2.0908D+00,
     4 2.7411D+00, 3.5335D+00, 4.7044D+00, 5.6578D+00,
     5 6.7367D+00, 8.3215D+00, 9.6199D+00, 1.1012D+01,
     6 1.3068D+01, 1.4683D+01, 1.5812D+01, 1.7584D+01,
     7 1.8817D+01, 2.0746D+01, 2.2769D+01, 2.4170D+01,
     8 2.5613D+01, 2.7097D+01, 2.8623D+01, 2.9402D+01,
     9 3.0991D+01, 3.3453D+01, 3.6009D+01, 3.8659D+01/
      DATA TART2/ 4.0478D+01, 4.2339D+01, 4.3285D+01, 4.6186D+01,
     1 4.7174D+01, 4.9181D+01, 5.1230D+01, 5.3321D+01,
     2 5.6536D+01, 5.7628D+01, 6.0968D+01, 6.3247D+01,
     3 6.5568D+01, 6.6744D+01, 7.0335D+01, 7.1553D+01,
     4 7.5270D+01, 7.7800D+01, 7.9080D+01, 8.1673D+01,
     5 8.4307D+01, 8.8337D+01, 9.1076D+01, 9.3857D+01,
     6 9.6680D+01, 9.8107D+01, 1.0245D+02, 1.0990D+02,
     7 1.1761D+02, 1.2558D+02, 1.3381D+02, 1.4059D+02,
     8 1.5106D+02, 1.6008D+02, 1.6936D+02, 1.7890D+02,
     9 1.8870D+02, 1.9876D+02, 2.0908D+02, 2.7411D+02/
      DATA TART3/ 3.2669D+02, 3.8105D+02, 4.7044D+02, 4.9908D+02,
     1 5.6578D+02, 6.0425D+02, 6.3666D+02, 7.1558D+02,
     2 8.3215D+02, 9.1767D+02, 1.0585D+03, 1.3068D+03,
     3 1.5812D+03, 1.8817D+03, 2.2084D+03, 2.5613D+03,
     4 2.9402D+03, 3.3453D+03, 3.7765D+03, 4.2339D+03,
     5 5.7628D+03, 7.5270D+03, 1.0245D+04, 1.5106D+04,
     6 2.0908D+04, 2.6462D+04, 3.2669D+04, 3.9530D+04,
     7 4.7044D+04, 5.7615D+04, 7.0020D+04, 8.3215D+04,
     8 9.8909D+04, 1.3068D+05, 1.8195D+05, 2.0746D+05,
     9 2.4170D+05, 2.7097D+05, 2.9402D+05, 3.3453D+05/
      DATA TART4/ 3.7765D+05, 4.2339D+05, 5.1230D+05, 6.3247D+05,
     1 7.5270D+05, 8.8337D+05, 1.0245D+06, 1.1761D+06,
     2 1.3381D+06, 1.5106D+06, 1.6936D+06, 1.8870D+06,
     3 2.0908D+06, 2.3051D+06, 2.5299D+06, 2.7411D+06,
     4 3.0108D+06, 3.2669D+06, 3.5335D+06, 3.8105D+06,
     5 4.0688D+06, 4.3960D+06, 4.7044D+06, 4.9908D+06,
     6 5.3525D+06, 5.6578D+06, 6.0425D+06, 6.3666D+06,
     7 6.7367D+06, 7.1558D+06, 7.5479D+06, 7.9096D+06,
     8 8.3215D+06, 8.7867D+06, 9.1767D+06, 9.6648D+06,
     9 1.0120D+07, 1.0585D+07, 1.1012D+07, 1.1547D+07/
      DATA TART5/ 1.1993D+07, 1.2499D+07, 1.3068D+07, 1.3542D+07,
     1 1.3863D+07, 1.4134D+07, 1.4407D+07, 1.4683D+07,
     2 1.5186D+07, 1.5754D+07, 1.6334D+07, 1.6923D+07,
     3 1.7523D+07, 1.8134D+07, 1.8755D+07, 2.0000D+07/
C
C     DEFINE ORNL 50, 126 AND 171 GROUP STRUCTURES.
C
C-----DEFINE 50 GROUP ENERGY BOUNDARIES.
      DATA GR50A/
     1 1.0000D-05, 6.8256D-01, 1.1250D+00, 1.8550D+00, 3.0590D+00,
     2 5.0430D+00, 8.3150D+00, 1.3710D+01, 2.2600D+01, 3.7270D+01,
     3 6.1440D+01, 1.0130D+02, 1.6700D+02, 2.7540D+02, 3.5360D+02,
     4 4.5400D+02, 5.8290D+02, 7.4850D+02, 9.6110D+02, 1.2340D+03,
     5 1.5850D+03, 2.0350D+03, 2.7470D+03, 3.3550D+03, 4.3070D+03,
     6 5.5310D+03, 7.1020D+03, 9.1190D+03, 1.1710D+04, 1.5030D+04,
     7 1.9300D+04, 2.4790D+04, 3.1830D+04, 4.0870D+04, 5.2480D+04,
     8 6.7380D+04, 8.6520D+04, 1.1110D+05, 1.4260D+05, 1.8320D+05,
     9 2.3520D+05, 3.0200D+05, 3.8770D+05, 4.9790D+05, 8.2080D+05/
      DATA GR50B/
     1 1.3530D+06, 2.2310D+06, 3.6790D+06, 6.0650D+06, 1.0000D+07,
     2 1.9970D+07/
C-----DEFINE 126 GROUP ENERGY BOUNDARIES.
      DATA GR126A/
     1 1.0000D-05, 1.0000D-01, 4.1399D-01, 1.1254D+00, 2.3724D+00,
     2 5.0435D+00, 1.0677D+01, 2.2603D+01, 3.7267D+01, 4.7851D+01,
     3 6.1442D+01, 1.0130D+02, 1.6702D+02, 2.1445D+02, 2.7536D+02,
     4 4.5400D+02, 7.4852D+02, 9.6112D+02, 1.2341D+03, 1.5846D+03,
     5 2.0347D+03, 2.2487D+03, 2.4852D+03, 2.6126D+03, 2.7465D+03,
     6 3.0354D+03, 3.3546D+03, 3.7074D+03, 4.3074D+03, 5.5308D+03,
     7 7.1017D+03, 9.1188D+03, 1.1709D+04, 1.5034D+04, 1.9305D+04,
     8 2.1875D+04, 2.3579D+04, 2.4788D+04, 2.6058D+04, 2.7000D+04,
     9 2.8500D+04, 3.1828D+04, 3.4307D+04, 4.0868D+04, 4.6309D+04/
      DATA GR126B/
     1 5.2475D+04, 5.6562D+04, 6.7379D+04, 7.2000D+04, 7.9500D+04,
     2 8.2500D+04, 8.6517D+04, 9.8037D+04, 1.1109D+05, 1.1679D+05,
     3 1.2277D+05, 1.2907D+05, 1.3569D+05, 1.4264D+05, 1.4996D+05,
     4 1.5764D+05, 1.6573D+05, 1.7422D+05, 1.8316D+05, 1.9255D+05,
     5 2.0242D+05, 2.1280D+05, 2.2371D+05, 2.4724D+05, 2.7324D+05,
     6 2.8725D+05, 2.9452D+05, 2.9720D+05, 2.9850D+05, 3.0197D+05,
     7 3.3373D+05, 3.6883D+05, 4.0762D+05, 4.5049D+05, 4.9787D+05,
     8 5.2340D+05, 5.5023D+05, 5.7844D+05, 6.0810D+05, 6.3928D+05,
     9 6.7206D+05, 7.0651D+05, 7.4274D+05, 7.8082D+05, 8.2085D+05/
      DATA GR126C/
     1 8.6294D+05, 9.0718D+05, 9.6164D+05, 1.0026D+06, 1.1080D+06,
     2 1.1648D+06, 1.2246D+06, 1.2873D+06, 1.3534D+06, 1.4227D+06,
     3 1.4957D+06, 1.5724D+06, 1.6530D+06, 1.7377D+06, 1.8268D+06,
     4 1.9205D+06, 2.0190D+06, 2.1225D+06, 2.2313D+06, 2.3069D+06,
     5 2.3653D+06, 2.3852D+06, 2.4660D+06, 2.5924D+06, 2.7253D+06,
     6 2.8650D+06, 3.0119D+06, 3.1664D+06, 3.6788D+06, 4.4933D+06,
     7 5.4881D+06, 6.0653D+06, 6.7032D+06, 8.1873D+06, 1.0000D+07,
     8 1.2214D+07, 1.7333D+07/
C-----DEFINE 171 GROUP ENERGY BOUNDARIES.
      DATA GR171A/
     1 1.0000D-05, 1.0000D-01, 4.1399D-01, 5.3158D-01, 6.8256D-01,
     2 8.7642D-01, 1.1254D+00, 1.4450D+00, 1.8554D+00, 2.3724D+00,
     3 3.0590D+00, 3.9279D+00, 5.0435D+00, 6.4760D+00, 8.3153D+00,
     4 1.0677D+01, 1.3710D+01, 1.7603D+01, 2.2603D+01, 2.9203D+01,
     5 3.7267D+01, 4.7851D+01, 6.1442D+01, 7.8893D+01, 1.0130D+02,
     6 1.3007D+02, 1.6702D+02, 2.1445D+02, 2.7536D+02, 3.5358D+02,
     7 4.5400D+02, 5.8295D+02, 7.4852D+02, 9.6112D+02, 1.2341D+03,
     8 1.5846D+03, 2.0347D+03, 2.2487D+03, 2.4852D+03, 2.6126D+03,
     9 2.7465D+03, 3.0354D+03, 3.3546D+03, 3.7074D+03, 4.3074D+03/
      DATA GR171B/
     1 5.5308D+03, 7.1017D+03, 9.1188D+03, 1.1709D+04, 1.5034D+04,
     2 1.9305D+04, 2.1875D+04, 2.3579D+04, 2.4176D+04, 2.4788D+04,
     3 2.6058D+04, 2.7000D+04, 2.8500D+04, 3.1828D+04, 3.4307D+04,
     4 4.0868D+04, 4.6309D+04, 5.2475D+04, 5.6562D+04, 6.7379D+04,
     5 7.2000D+04, 7.9500D+04, 8.2500D+04, 8.6517D+04, 9.8037D+04,
     6 1.1109D+05, 1.1679D+05, 1.2277D+05, 1.2907D+05, 1.3569D+05,
     7 1.4264D+05, 1.4996D+05, 1.5764D+05, 1.6573D+05, 1.7422D+05,
     8 1.8316D+05, 1.9255D+05, 2.0242D+05, 2.1280D+05, 2.2371D+05,
     9 2.3518D+05, 2.4724D+05, 2.7324D+05, 2.8725D+05, 2.9452D+05/
      DATA GR171C/
     1 2.9720D+05, 2.9850D+05, 3.0197D+05, 3.3373D+05, 3.6883D+05,
     2 3.8774D+05, 4.0762D+05, 4.5049D+05, 4.9787D+05, 5.2340D+05,
     3 5.5023D+05, 5.7844D+05, 6.0810D+05, 6.3928D+05, 6.7206D+05,
     4 7.0651D+05, 7.4274D+05, 7.8082D+05, 8.2085D+05, 8.6294D+05,
     5 9.0718D+05, 9.6164D+05, 1.0026D+06, 1.1080D+06, 1.1648D+06,
     6 1.2246D+06, 1.2873D+06, 1.3534D+06, 1.4227D+06, 1.4957D+06,
     7 1.5724D+06, 1.6530D+06, 1.7377D+06, 1.8268D+06, 1.9205D+06,
     8 2.0190D+06, 2.1225D+06, 2.2313D+06, 2.3069D+06, 2.3457D+06,
     9 2.3653D+06, 2.3852D+06, 2.4660D+06, 2.5924D+06, 2.7253D+06/
      DATA GR171D/
     1 2.8650D+06, 3.0119D+06, 3.1664D+06, 3.3287D+06, 3.6788D+06,
     2 4.0657D+06, 4.4933D+06, 4.7237D+06, 4.9659D+06, 5.2205D+06,
     3 5.4881D+06, 5.7695D+06, 6.0653D+06, 6.3763D+06, 6.5924D+06,
     4 6.7032D+06, 7.0469D+06, 7.4082D+06, 7.7880D+06, 8.1873D+06,
     5 8.6071D+06, 9.0484D+06, 9.5123D+06, 1.0000D+07, 1.0513D+07,
     6 1.1052D+07, 1.1618D+07, 1.2214D+07, 1.2840D+07, 1.3499D+07,
     7 1.3840D+07, 1.4191D+07, 1.4550D+07, 1.4918D+07, 1.5683D+07,
     8 1.6487D+07, 1.7333D+07/
C
C     DEFINE WIMS 69 GROUP STRUCTURE.
C
      DATA WIMSA/
     1 0.0001D0  , 0.0050D0  , 0.0100D0  , 0.0150D0  , 0.0200D0  ,
     2 0.0250D0  , 0.0300D0  , 0.0350D0  , 0.0420D0  , 0.0500D0  ,
     3 0.0580D0  , 0.0670D0  , 0.0800D0  , 0.1000D0  , 0.1400D0  ,
     4 0.1800D0  , 0.2200D0  , 0.2500D0  , 0.2800D0  , 0.3000D0  ,
     5 0.3200D0  , 0.3500D0  , 0.4000D0  , 0.5000D0  , 0.6250D0  ,
     6 0.7800D0  , 0.8500D0  , 0.9100D0  , 0.9500D0  , 0.9720D0  ,
     7 0.9960D0  , 1.0200D0  , 1.0450D0  , 1.0710D0  , 1.0970D0  ,
     8 1.1230D0  , 1.1500D0  , 1.3000D0  , 1.5000D0  , 2.1000D0  ,
     9 2.6000D0  , 3.3000D0  , 4.0000D0  , 9.8770D0  , 15.968D0  /
      DATA WIMSB/
     1 27.7000D0  , 48.0520D0   , 75.5014D0  , 148.728D0 , 367.262D0 ,
     2 906.898D0  , 1425.10D0   , 2239.45D0  , 3519.10D0 , 5530.00D0 ,
     3 9118.00D0  , 15030.0D0   , 24780.0D0  , 40850.0D0 , 67340.0D0 ,
     4 111000.0D0 , 183000.0D0  , 302000.0D0 , 500000.0D0, 821000.0D0,
     5 1353000.0D0, 2231000.0D0 , 3679000.0D0,
     6 6065000.0D0, 10000000.0D0/
C
C     DEFINE GAM-I 68 GROUP STRUCTURE.
C
      DATA GAMIA/
     1 0.414D0,    0.532D0,    0.683D0,    0.876D0,    1.125D0,
     2 1.44D0,     1.86D0,     2.38D0,     3.06D0,     3.93D0,
     3 5.04D0,     6.48D0,     8.32D0,    10.68D0,    13.7D0,
     4 17.6D0,     22.6D0,     29.0D0,     37.3D0,     47.9D0,
     5 61.4D0,     78.9D0,    101.0D0,    130.0D0,    167.0D0,
     6 215.0D0,    275.0D0,    354.0D0,    454.0D0,    583.0D0,
     7 748.0D0,    961.0D0,   1230.0D0,   1590.0D0,   2040.0D0,
     8 2610.0D0,   3360.0D0,   4310.0D0,   5530.0D0,   7100.0D0,
     9 9120.0D0,  11700.0D0,  15000.0 D0, 19300.0D0,  24800.0D0/
      DATA GAMIB/
     1 31800.0D0 , 40900.0D0 , 52500.0D0 , 67400.0D0 , 86500.0D0 ,
     2 111000.0D0, 143000.0D0, 183000.0D0, 235000.0D0, 302000.0D0,
     3 388000.0D0, 498000.0D0, 639000.0D0, 821000.0D0, 1.0500D+06,
     4 1.3600D+06, 1.7400D+06, 2.2300D+06, 2.8700D+06, 3.6800D+06,
     5 4.7200D+06, 6.0700D+06, 7.7900D+06, 1.0000D+07/
C
C     DEFINE GAM-II 99 GROUP STRUCTURE.
C
      DATA GAMIIA/
     1 0.414D0,    0.532D0,    0.683D0,    0.876D0,    1.125D0,
     2 1.44D0,     1.86D0,     2.38D0,     3.06D0,     3.93D0,
     3 5.04D0,     6.48D0,     8.32D0,    10.68D0,    13.7D0,
     4 17.6D0,     22.6D0,     29.0D0,     37.3D0,     47.9D0,
     5 61.4D0,     78.9D0,    101.0D0,    130.0D0,    167.0D0,
     6 215.0D0,    275.0D0,    354.0D0,    454.0D0,    583.0D0,
     7 748.0D0,    961.0D0,   1230.0D0,   1590.0D0,   2040.0D0,
     8 2610.0D0,   3360.0D0,   4310.0D0,   5530.0D0,   7100.0D0,
     9 9120.0D0,  11700.0D0,  15000.0D0,  19300.0D0,  24800.0D0/
      DATA GAMIIB/
     1 31800.00D0,  40900.0D0,  52500.0D0,  67400.0D0,  86500.0D0,
     2 111000.0D0, 128000.0D0, 136000.0D0, 150000.0D0, 166000.0D0,
     3 183000.0D0, 202000.0D0, 224000.0D0, 247000.0D0, 273000.0D0,
     4 302000.0D0, 334000.0D0, 369000.0D0, 408000.0D0, 450000.0D0,
     5 498000.0D0, 550000.0D0, 608000.0D0, 672000.0D0, 743000.0D0,
     6 821000.0D0, 907000.0D0, 1.0000D+06, 1.1100D+06, 1.2200D+06,
     7 1.3500D+06, 1.5000D+06, 1.6500D+06, 1.8300D+06, 2.0200D+06,
     8 2.2300D+06, 2.4700D+06, 2.7300D+06, 3.0100D+06, 3.3300D+06,
     9 3.6800D+06, 4.0700D+06, 4.4900D+06, 4.9600D+06, 5.4900D+06/
      DATA GAMIIC/
     1 6.0700D+06, 6.7000D+06, 7.4100D+06, 8.1900D+06, 9.0500D+06,
     2 1.0000D+07, 1.1100D+07, 1.2200D+07, 1.3500D+07, 1.4900D+07/
C
C     DEFINE MUFT 54 GROUP STRUCTURE.
C
      DATA XMUFT/
     1   0.625D+0,   0.835D+0,   1.125D+0,   1.440D+0,   1.855D+0,
     2    2.38D+0,    3.06D+0,    3.97D+0,    5.10D+0,    6.50D+0,
     3    8.32D+0,    10.7D+0,    13.7D+0,    17.6D+0,    22.6D+0,
     4    29.0D+0,    37.2D+0,    47.8D+0,    61.3D+0,    78.7D+0,
     5   101.0D+0,   130.0D+0,   167.0D+0,   275.0D+0,   454.0D+0,
     6   750.0D+0,  1230.0D+0,  2030.0D+0,  3350.0D+0,  5530.0D+0,
     7  9120.0D+0, 15000.0D+0, 24800.0D+0, 40900.0D+0, 67400.0D+0,
     8 86500.0D+0, 1.11000D+5, 1.43000D+5, 1.83000D+5, 2.35000D+5,
     9 3.02000D+5, 3.87000D+5, 4.98000D+5, 6.39000D+5, 8.21000D+5,
     A 1.0500D+06, 1.3500D+06, 1.7400D+06, 2.2300D+06, 2.8600D+06,
     1 3.6800D+06, 4.7200D+06, 6.0700D+06, 7.7900D+06, 1.0000D+07/
C
C    ABBN GROUP STRUCTURE (28 GROUPS - NARROW GROUP NEAR THERMAL,
C    DUMMY GROUP FROM 0.0256 TO 0.215 EV, 26 GROUPS UP TO 15 MEV).
C
      DATA ABBN/  2.500D-02, 2.560D-02,
     1 2.150D-01, 4.650D-01, 1.000D+00, 2.150D+00, 4.650D+00, 1.000D+01,
     1 2.150D+01, 4.650D+01, 1.000D+02, 2.150D+02, 4.650D+02, 1.000D+03,
     1 2.150D+03, 4.650D+03, 1.000D+04, 2.150D+04, 4.650D+04, 1.000D+05,
     1 2.000D+05, 4.000D+05, 8.000D+05, 1.400D+06, 2.500D+06, 4.000D+06,
     1 6.500D+06, 1.050D+07, 1.500D+07/
C
C     CHECK FOR ALLOWABLE VALUES - IF NOT, RETURN NO GROUPS.
C
      IF(IWANT.GE.0.AND.IWANT.LE.19) GO TO 10
      NGROUP=0
      RETURN
C
C     SELECT GROUP STRUCTURE
C
   10 II=IWANT+1
c              1   2   3   4   5   6   7   8   9  10
      GO TO ( 20, 40, 60, 80,100,110,120,140,160,180,
     1       200,220,230,250,260,240,270,280,290,300),II
C-----TART 176 GROUPS
   20 DO 30 I=1,176
   30 EGROUP(I)=TART(I)
      NGROUP=175
      RETURN
C-----ORNL 50 GROUPS
   40 DO 50 I=1,51
   50 EGROUP(I)=ORNLA(I)
      NGROUP=50
      RETURN
C-----ORNL 126 GROUPS
   60 DO 70 I=1,127
   70 EGROUP(I)=ORNLB(I)
      NGROUP=126
      RETURN
C-----ORNL 171 GROUPS
   80 DO 90 I=1,172
   90 EGROUP(I)=ORNLC(I)
      NGROUP=171
      RETURN
C
C     SAND-II - Original - 1.0D-4 eV to 18 or 20 MeV.
C
C-----SAND 620 GROUPS (1.0D-4 eV to 18 MeV)
  100 EGRPMIN = 1.0D-4
      EGRPMAX = 1.8D+7
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C-----SAND 640 GROUPS (1.0D-4 eV to 20 MeV)
  110 EGRPMIN = 1.0D-4
      EGRPMAX = 2.0D+7
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C-----WIMS 69 GROUPS
  120 DO 130 I=1,70
  130 EGROUP(I)=WIMS(I)
      NGROUP=69
      RETURN
C-----GAM-I 68 GROUPS
  140 DO 150 I=1,69
  150 EGROUP(I)=GAMI(I)
      NGROUP=68
      RETURN
C-----GAM-II 99 GROUPS
  160 DO 170 I=1,100
  170 EGROUP(I)=GAMII(I)
      NGROUP=99
      RETURN
C-----MUFT 54 GROUPS
  180 DO 190 I=1,55
  190 EGROUP(I)=XMUFT(I)
      NGROUP=54
      RETURN
C-----ABBN 28 GROUPS
  200 DO 210 I=1,29
  210 EGROUP(I)=ABBN(I)
      NGROUP=28
      RETURN
C-----TART 616 GROUPS UP TO 20 MeV
  220 EGRPMAX = 2.0D+7
      CALL TARTGEN(EGROUP,NGROUP,EGRPMAX)
      RETURN
C-----TART 700 GROUPS UP TO 1 GeV
  230 EGRPMAX = 1.0D+9
      CALL TARTGEN(EGROUP,NGROUP,EGRPMAX)
      RETURN
C-----TART 666 GROUPS UP TO 200 MeV
  240 EGRPMAX = 2.0D+8
      CALL TARTGEN(EGROUP,NGROUP,EGRPMAX)
      RETURN
C
C     SAND-II - Extended down to 1.0D-5 up to 18 or 20 MeV.
C
C-----SAND 665 GROUPS (1.0D-5 eV to 18 MeV)
  250 EGRPMIN = 1.0D-5
      EGRPMAX = 1.8D+7
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C-----SAND 685 GROUPS (1.0D-5 eV to 20 MeV)
  260 EGRPMIN = 1.0D-5
      EGRPMAX = 2.0D+7
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C
C     SAND-II - Extended up to 60, 150 or 200 MeV.
C
C-----SAND 725 GROUPS (1.0D-5 eV to 60 MeV)
  270 EGRPMIN = 1.0D-5
      EGRPMAX = 6.0D+7
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C-----SAND 755 GROUPS (1.0D-5 eV to 150 MeV)
  280 EGRPMIN = 1.0D-5
      EGRPMAX = 1.5D+8
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C-----SAND 765 GROUPS (1.0D-5 eV to 200 MeV)
  290 EGRPMIN = 1.0D-5
      EGRPMAX = 2.0D+8
      CALL SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
      RETURN
C
C     UKAEA 1102 GROUPS (1.0D-5 eV to 1 GeV)
C
  300 EGRPMAX = 1.0D+9
      CALL UKAEAGEN(EGROUP,NGROUP,EGRPMAX)
      RETURN
      END
      SUBROUTINE TARTGEN(EGROUP,NGROUP,EGRPMAX)
C===============================================================
C
C     DEFINE TART GROUP STRUCTURE FROM 1.0D-5 eV TO UP
C     ANY UPPER LIMIT DEFINED BY EGRPMAX INPUT.
C
C===============================================================
      INCLUDE 'implicit.h'
      DIMENSION EGROUP(*)
c-----2011/11 - Corrected lower energy to 1.0D-5 from 1.0D-04.
      DATA EGRPMIN/1.0D-05/
      DATA TEN/10.0D+00/
      DATA FIFTY/50.0D+00/
C-----DEFAULT IS 700 GROUPS
      NGROUP=700
C-----50 PER ENERGY DECADE
      DE=DEXP(DLOG(TEN)/FIFTY)
      ENOW=EGRPMIN
C-----DEFINE GROUPS UP TO MAXIMUM ENDL ENERGY
      DO 10 I=1,NGROUP+1
      EGROUP(I)=ENOW
      IF(EGROUP(I).eq.EGRPMAX) go to 30
      IF(EGROUP(I).gt.EGRPMAX) go to 20
   10 ENOW=ENOW*DE
      RETURN
C-----DEFINE END OF GROUPS AT MAXIMUM ENDL ENERGY
   20 EGROUP(I)=EGRPMAX
   30 NGROUP=I-1
      RETURN
      END
      SUBROUTINE SAND2GEN(EGROUP,NGROUP,EGRPMIN,EGRPMAX)
C=======================================================================
C
C     GENERAL SAND-II GROUP STRUCTURE TO SELECT ALL GROUPS
C     BETWEEN EGRPMIN AND EGRPMAX (INPUT PARAMETERS).
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION EGROUP(*),SAND(766),
     1 SANDA(45),SANDB(45),SANDC(45),SANDD(45),SANDE(45),SANDF(45),
     2 SANDG(45),SANDH(45),SANDI(45),SANDJ(45),SANDK(45),SANDL(45),
     3 SANDM(45),SANDN(45),SANDO(36),SANDP(20),SANDQ(40),SANDR(30),
     4 SANDS(10)
c
      EQUIVALENCE (SAND(  1),SANDA(1)),    ! 45
     1            (SAND( 46),SANDB(1)),    ! 45
     1            (SAND( 91),SANDC(1)),    ! 45
     1            (SAND(136),SANDD(1)),    ! 45
     1            (SAND(181),SANDE(1)),    ! 45
     1            (SAND(226),SANDF(1)),    ! 45
     1            (SAND(271),SANDG(1)),    ! 45
     1            (SAND(316),SANDH(1)),    ! 45
     1            (SAND(361),SANDI(1)),    ! 45
     1            (SAND(406),SANDJ(1)),    ! 45
     1            (SAND(451),SANDK(1)),    ! 45
     1            (SAND(496),SANDL(1)),    ! 45
     1            (SAND(541),SANDM(1)),    ! 45
     1            (SAND(586),SANDN(1)),    ! 45
     1            (SAND(631),SANDO(1)),    ! 36
     1            (SAND(667),SANDP(1)),    ! 20
     1            (SAND(687),SANDQ(1)),    ! 40
     1            (SAND(727),SANDR(1)),    ! 30
     1            (SAND(757),SANDS(1))     ! 10
c                       766
C
C     DEFINE SAND-II 45 GROUP EXTENSION DOWN TO 10D-5 EV
C
      DATA SANDA/
     1 1.0000D-05, 1.0500D-05, 1.1000D-05, 1.1500D-05, 1.2000D-05,
     2 1.2750D-05, 1.3500D-05, 1.4250D-05, 1.5000D-05, 1.6000D-05,
     3 1.7000D-05, 1.8000D-05, 1.9000D-05, 2.0000D-05, 2.1000D-05,
     4 2.2000D-05, 2.3000D-05, 2.4000D-05, 2.5500D-05, 2.7000D-05,
     5 2.8000D-05, 3.0000D-05, 3.2000D-05, 3.4000D-05, 3.6000D-05,
     6 3.8000D-05, 4.0000D-05, 4.2500D-05, 4.5000D-05, 4.7500D-05,
     7 5.0000D-05, 5.2500D-05, 5.5000D-05, 5.7500D-05, 6.0000D-05,
     8 6.3000D-05, 6.6000D-05, 6.9000D-05, 7.2000D-05, 7.6000D-05,
     9 8.0000D-05, 8.4000D-05, 8.8000D-05, 9.2000D-05, 9.6000D-05/
C
C     DEFINE SAND-II 620 (665) AND 620 (685) GROUP STRUCTURES.
C
      DATA SANDB/
     1 1.0000D-04, 1.0500D-04, 1.1000D-04, 1.1500D-04, 1.2000D-04,
     2 1.2750D-04, 1.3500D-04, 1.4250D-04, 1.5000D-04, 1.6000D-04,
     3 1.7000D-04, 1.8000D-04, 1.9000D-04, 2.0000D-04, 2.1000D-04,
     4 2.2000D-04, 2.3000D-04, 2.4000D-04, 2.5500D-04, 2.7000D-04,
     5 2.8000D-04, 3.0000D-04, 3.2000D-04, 3.4000D-04, 3.6000D-04,
     6 3.8000D-04, 4.0000D-04, 4.2500D-04, 4.5000D-04, 4.7500D-04,
     7 5.0000D-04, 5.2500D-04, 5.5000D-04, 5.7500D-04, 6.0000D-04,
     8 6.3000D-04, 6.6000D-04, 6.9000D-04, 7.2000D-04, 7.6000D-04,
     9 8.0000D-04, 8.4000D-04, 8.8000D-04, 9.2000D-04, 9.6000D-04/
      DATA SANDC/
     1 1.0000D-03, 1.0500D-03, 1.1000D-03, 1.1500D-03, 1.2000D-03,
     2 1.2750D-03, 1.3500D-03, 1.4250D-03, 1.5000D-03, 1.6000D-03,
     3 1.7000D-03, 1.8000D-03, 1.9000D-03, 2.0000D-03, 2.1000D-03,
     4 2.2000D-03, 2.3000D-03, 2.4000D-03, 2.5500D-03, 2.7000D-03,
     5 2.8000D-03, 3.0000D-03, 3.2000D-03, 3.4000D-03, 3.6000D-03,
     6 3.8000D-03, 4.0000D-03, 4.2500D-03, 4.5000D-03, 4.7500D-03,
     7 5.0000D-03, 5.2500D-03, 5.5000D-03, 5.7500D-03, 6.0000D-03,
     8 6.3000D-03, 6.6000D-03, 6.9000D-03, 7.2000D-03, 7.6000D-03,
     9 8.0000D-03, 8.4000D-03, 8.8000D-03, 9.2000D-03, 9.6000D-03/
      DATA SANDD/
     1 1.0000D-02, 1.0500D-02, 1.1000D-02, 1.1500D-02, 1.2000D-02,
     2 1.2750D-02, 1.3500D-02, 1.4250D-02, 1.5000D-02, 1.6000D-02,
     3 1.7000D-02, 1.8000D-02, 1.9000D-02, 2.0000D-02, 2.1000D-02,
     4 2.2000D-02, 2.3000D-02, 2.4000D-02, 2.5500D-02, 2.7000D-02,
     5 2.8000D-02, 3.0000D-02, 3.2000D-02, 3.4000D-02, 3.6000D-02,
     6 3.8000D-02, 4.0000D-02, 4.2500D-02, 4.5000D-02, 4.7500D-02,
     7 5.0000D-02, 5.2500D-02, 5.5000D-02, 5.7500D-02, 6.0000D-02,
     8 6.3000D-02, 6.6000D-02, 6.9000D-02, 7.2000D-02, 7.6000D-02,
     9 8.0000D-02, 8.4000D-02, 8.8000D-02, 9.2000D-02, 9.6000D-02/
      DATA SANDE/
     1 1.0000D-01, 1.0500D-01, 1.1000D-01, 1.1500D-01, 1.2000D-01,
     2 1.2750D-01, 1.3500D-01, 1.4250D-01, 1.5000D-01, 1.6000D-01,
     3 1.7000D-01, 1.8000D-01, 1.9000D-01, 2.0000D-01, 2.1000D-01,
     4 2.2000D-01, 2.3000D-01, 2.4000D-01, 2.5500D-01, 2.7000D-01,
     5 2.8000D-01, 3.0000D-01, 3.2000D-01, 3.4000D-01, 3.6000D-01,
     6 3.8000D-01, 4.0000D-01, 4.2500D-01, 4.5000D-01, 4.7500D-01,
     7 5.0000D-01, 5.2500D-01, 5.5000D-01, 5.7500D-01, 6.0000D-01,
     8 6.3000D-01, 6.6000D-01, 6.9000D-01, 7.2000D-01, 7.6000D-01,
     9 8.0000D-01, 8.4000D-01, 8.8000D-01, 9.2000D-01, 9.6000D-01/
      DATA SANDF/
     1 1.0000D+00, 1.0500D+00, 1.1000D+00, 1.1500D+00, 1.2000D+00,
     2 1.2750D+00, 1.3500D+00, 1.4250D+00, 1.5000D+00, 1.6000D+00,
     3 1.7000D+00, 1.8000D+00, 1.9000D+00, 2.0000D+00, 2.1000D+00,
     4 2.2000D+00, 2.3000D+00, 2.4000D+00, 2.5500D+00, 2.7000D+00,
     5 2.8000D+00, 3.0000D+00, 3.2000D+00, 3.4000D+00, 3.6000D+00,
     6 3.8000D+00, 4.0000D+00, 4.2500D+00, 4.5000D+00, 4.7500D+00,
     7 5.0000D+00, 5.2500D+00, 5.5000D+00, 5.7500D+00, 6.0000D+00,
     8 6.3000D+00, 6.6000D+00, 6.9000D+00, 7.2000D+00, 7.6000D+00,
     9 8.0000D+00, 8.4000D+00, 8.8000D+00, 9.2000D+00, 9.6000D+00/
      DATA SANDG/
     1 1.0000D+01, 1.0500D+01, 1.1000D+01, 1.1500D+01, 1.2000D+01,
     2 1.2750D+01, 1.3500D+01, 1.4250D+01, 1.5000D+01, 1.6000D+01,
     3 1.7000D+01, 1.8000D+01, 1.9000D+01, 2.0000D+01, 2.1000D+01,
     4 2.2000D+01, 2.3000D+01, 2.4000D+01, 2.5500D+01, 2.7000D+01,
     5 2.8000D+01, 3.0000D+01, 3.2000D+01, 3.4000D+01, 3.6000D+01,
     6 3.8000D+01, 4.0000D+01, 4.2500D+01, 4.5000D+01, 4.7500D+01,
     7 5.0000D+01, 5.2500D+01, 5.5000D+01, 5.7500D+01, 6.0000D+01,
     8 6.3000D+01, 6.6000D+01, 6.9000D+01, 7.2000D+01, 7.6000D+01,
     9 8.0000D+01, 8.4000D+01, 8.8000D+01, 9.2000D+01, 9.6000D+01/
      DATA SANDH/
     1 1.0000D+02, 1.0500D+02, 1.1000D+02, 1.1500D+02, 1.2000D+02,
     2 1.2750D+02, 1.3500D+02, 1.4250D+02, 1.5000D+02, 1.6000D+02,
     3 1.7000D+02, 1.8000D+02, 1.9000D+02, 2.0000D+02, 2.1000D+02,
     4 2.2000D+02, 2.3000D+02, 2.4000D+02, 2.5500D+02, 2.7000D+02,
     5 2.8000D+02, 3.0000D+02, 3.2000D+02, 3.4000D+02, 3.6000D+02,
     6 3.8000D+02, 4.0000D+02, 4.2500D+02, 4.5000D+02, 4.7500D+02,
     7 5.0000D+02, 5.2500D+02, 5.5000D+02, 5.7500D+02, 6.0000D+02,
     8 6.3000D+02, 6.6000D+02, 6.9000D+02, 7.2000D+02, 7.6000D+02,
     9 8.0000D+02, 8.4000D+02, 8.8000D+02, 9.2000D+02, 9.6000D+02/
      DATA SANDI/
     1 1.0000D+03, 1.0500D+03, 1.1000D+03, 1.1500D+03, 1.2000D+03,
     2 1.2750D+03, 1.3500D+03, 1.4250D+03, 1.5000D+03, 1.6000D+03,
     3 1.7000D+03, 1.8000D+03, 1.9000D+03, 2.0000D+03, 2.1000D+03,
     4 2.2000D+03, 2.3000D+03, 2.4000D+03, 2.5500D+03, 2.7000D+03,
     5 2.8000D+03, 3.0000D+03, 3.2000D+03, 3.4000D+03, 3.6000D+03,
     6 3.8000D+03, 4.0000D+03, 4.2500D+03, 4.5000D+03, 4.7500D+03,
     7 5.0000D+03, 5.2500D+03, 5.5000D+03, 5.7500D+03, 6.0000D+03,
     8 6.3000D+03, 6.6000D+03, 6.9000D+03, 7.2000D+03, 7.6000D+03,
     9 8.0000D+03, 8.4000D+03, 8.8000D+03, 9.2000D+03, 9.6000D+03/
      DATA SANDJ/
     1 1.0000D+04, 1.0500D+04, 1.1000D+04, 1.1500D+04, 1.2000D+04,
     2 1.2750D+04, 1.3500D+04, 1.4250D+04, 1.5000D+04, 1.6000D+04,
     3 1.7000D+04, 1.8000D+04, 1.9000D+04, 2.0000D+04, 2.1000D+04,
     4 2.2000D+04, 2.3000D+04, 2.4000D+04, 2.5500D+04, 2.7000D+04,
     5 2.8000D+04, 3.0000D+04, 3.2000D+04, 3.4000D+04, 3.6000D+04,
     6 3.8000D+04, 4.0000D+04, 4.2500D+04, 4.5000D+04, 4.7500D+04,
     7 5.0000D+04, 5.2500D+04, 5.5000D+04, 5.7500D+04, 6.0000D+04,
     8 6.3000D+04, 6.6000D+04, 6.9000D+04, 7.2000D+04, 7.6000D+04,
     9 8.0000D+04, 8.4000D+04, 8.8000D+04, 9.2000D+04, 9.6000D+04/
      DATA SANDK/
     1 1.0000D+05, 1.0500D+05, 1.1000D+05, 1.1500D+05, 1.2000D+05,
     2 1.2750D+05, 1.3500D+05, 1.4250D+05, 1.5000D+05, 1.6000D+05,
     3 1.7000D+05, 1.8000D+05, 1.9000D+05, 2.0000D+05, 2.1000D+05,
     4 2.2000D+05, 2.3000D+05, 2.4000D+05, 2.5500D+05, 2.7000D+05,
     5 2.8000D+05, 3.0000D+05, 3.2000D+05, 3.4000D+05, 3.6000D+05,
     6 3.8000D+05, 4.0000D+05, 4.2500D+05, 4.5000D+05, 4.7500D+05,
     7 5.0000D+05, 5.2500D+05, 5.5000D+05, 5.7500D+05, 6.0000D+05,
     8 6.3000D+05, 6.6000D+05, 6.9000D+05, 7.2000D+05, 7.6000D+05,
     9 8.0000D+05, 8.4000D+05, 8.8000D+05, 9.2000D+05, 9.6000D+05/
      DATA SANDL/
     1 1.0000D+06, 1.1000D+06, 1.2000D+06, 1.3000D+06, 1.4000D+06,
     2 1.5000D+06, 1.6000D+06, 1.7000D+06, 1.8000D+06, 1.9000D+06,
     3 2.0000D+06, 2.1000D+06, 2.2000D+06, 2.3000D+06, 2.4000D+06,
     4 2.5000D+06, 2.6000D+06, 2.7000D+06, 2.8000D+06, 2.9000D+06,
     5 3.0000D+06, 3.1000D+06, 3.2000D+06, 3.3000D+06, 3.4000D+06,
     6 3.5000D+06, 3.6000D+06, 3.7000D+06, 3.8000D+06, 3.9000D+06,
     7 4.0000D+06, 4.1000D+06, 4.2000D+06, 4.3000D+06, 4.4000D+06,
     8 4.5000D+06, 4.6000D+06, 4.7000D+06, 4.8000D+06, 4.9000D+06,
     9 5.0000D+06, 5.1000D+06, 5.2000D+06, 5.3000D+06, 5.4000D+06/
      DATA SANDM/
     1 5.5000D+06, 5.6000D+06, 5.7000D+06, 5.8000D+06, 5.9000D+06,
     2 6.0000D+06, 6.1000D+06, 6.2000D+06, 6.3000D+06, 6.4000D+06,
     3 6.5000D+06, 6.6000D+06, 6.7000D+06, 6.8000D+06, 6.9000D+06,
     4 7.0000D+06, 7.1000D+06, 7.2000D+06, 7.3000D+06, 7.4000D+06,
     5 7.5000D+06, 7.6000D+06, 7.7000D+06, 7.8000D+06, 7.9000D+06,
     6 8.0000D+06, 8.1000D+06, 8.2000D+06, 8.3000D+06, 8.4000D+06,
     7 8.5000D+06, 8.6000D+06, 8.7000D+06, 8.8000D+06, 8.9000D+06,
     8 9.0000D+06, 9.1000D+06, 9.2000D+06, 9.3000D+06, 9.4000D+06,
     9 9.5000D+06, 9.6000D+06, 9.7000D+06, 9.8000D+06, 9.9000D+06/
      DATA SANDN/
     1 1.0000D+07, 1.0100D+07, 1.0200D+07, 1.0300D+07, 1.0400D+07,
     2 1.0500D+07, 1.0600D+07, 1.0700D+07, 1.0800D+07, 1.0900D+07,
     3 1.1000D+07, 1.1100D+07, 1.1200D+07, 1.1300D+07, 1.1400D+07,
     4 1.1500D+07, 1.1600D+07, 1.1700D+07, 1.1800D+07, 1.1900D+07,
     5 1.2000D+07, 1.2100D+07, 1.2200D+07, 1.2300D+07, 1.2400D+07,
     6 1.2500D+07, 1.2600D+07, 1.2700D+07, 1.2800D+07, 1.2900D+07,
     7 1.3000D+07, 1.3100D+07, 1.3200D+07, 1.3300D+07, 1.3400D+07,
     8 1.3500D+07, 1.3600D+07, 1.3700D+07, 1.3800D+07, 1.3900D+07,
     9 1.4000D+07, 1.4100D+07, 1.4200D+07, 1.4300D+07, 1.4400D+07/
      DATA SANDO/
     1 1.4500D+07, 1.4600D+07, 1.4700D+07, 1.4800D+07, 1.4900D+07,
     2 1.5000D+07, 1.5100D+07, 1.5200D+07, 1.5300D+07, 1.5400D+07,
     3 1.5500D+07, 1.5600D+07, 1.5700D+07, 1.5800D+07, 1.5900D+07,
     4 1.6000D+07, 1.6100D+07, 1.6200D+07, 1.6300D+07, 1.6400D+07,
     5 1.6500D+07, 1.6600D+07, 1.6700D+07, 1.6800D+07, 1.6900D+07,
     6 1.7000D+07, 1.7100D+07, 1.7200D+07, 1.7300D+07, 1.7400D+07,
     7 1.7500D+07, 1.7600D+07, 1.7700D+07, 1.7800D+07, 1.7900D+07,
     8 1.8000D+07/
C------DEFINE EXTENSION OF SAND-II STRUCTURE FROM 18 TO 20 MEV USING
C------20 ADDITIONAL GROUPS EACH 0.1 MEV WIDTH.
      DATA SANDP/
     1 1.8100D+07, 1.8200D+07, 1.8300D+07, 1.8400D+07, 1.8500D+07,
     1 1.8600D+07, 1.8700D+07, 1.8800D+07, 1.8900D+07, 1.9000D+07,
     1 1.9100D+07, 1.9200D+07, 1.9300D+07, 1.9400D+07, 1.9500D+07,
     1 1.9600D+07, 1.9700D+07, 1.9800D+07, 1.9900D+07, 2.0000D+07/
C------DEFINE EXTENSION OF SAND-II STRUCTURE FROM 20 TO 60 MEV USING
C------20 ADDITIONAL GROUPS EACH 0.5 MEV WIDTH.
C------10 ADDITIONAL GROUPS EACH 1.0 MEV WIDTH.
C------10 ADDITIONAL GROUPS EACH 2.0 MEV WIDTH.
      DATA SANDQ/
     1 20.500D+06, 21.000D+06, 21.500D+06, 22.000D+06, 22.500D+06,
     1 23.000D+06, 23.500D+06, 24.000D+06, 24.500D+06, 25.000D+06,
     1 25.500D+06, 26.000D+06, 26.500D+06, 27.000D+06, 27.500D+06,
     1 28.000D+06, 28.500D+06, 29.000D+06, 29.500D+06, 30.000D+06,
     1 31.000D+06, 32.000D+06, 33.000D+06, 34.000D+06, 35.000D+06,
     1 36.000D+06, 37.000D+06, 38.000D+06, 39.000D+06, 40.000D+06,
     1 42.000D+06, 44.000D+06, 46.000D+06, 48.000D+06, 50.000D+06,
     1 52.000D+06, 54.000D+06, 56.000D+06, 58.000D+06, 60.000D+06/
C------DEFINE EXTENSION OF SAND-II STRUCTURE FROM 60 TO 150 MEV USING
C------20 ADDITIONAL GROUPS EACH 2.0 MEV WIDTH.
C------10 ADDITIONAL GROUPS EACH 5.0 MEV WIDTH.
      DATA SANDR/
     1 62.000D+06, 64.000D+06, 66.000D+06, 68.000D+06, 70.000D+06,
     1 72.000D+06, 74.000D+06, 76.000D+06, 78.000D+06, 80.000D+06,
     1 82.000D+06, 84.000D+06, 86.000D+06, 88.000D+06, 90.000D+06,
     1 92.000D+06, 94.000D+06, 96.000D+06, 98.000D+06,100.000D+06,
     1105.000D+06,110.000D+06,115.000D+06,120.000D+06,125.000D+06,
     1130.000D+06,135.000D+06,140.000D+06,145.000D+06,150.000D+06/
C------DEFINE EXTENSION OF SAND-II STRUCTURE FROM 150 TO 200 MEV USING
C------10 ADDITIONAL GROUPS EACH 5.0 MEV WIDTH.
      DATA SANDS/
     1 155.000D+06,160.000D+06,165.000D+06,170.000D+06,175.000D+06,
     1 180.000D+06,185.000D+06,190.000D+06,195.000D+06,200.000D+06/
C
C     USE ALL GROUPS BETWEEN EGRPMIN AND EGRPMAX.
C
      ig = 0
      do i=1,766
      if(SAND(i).ge.EGRPMIN.and.SAND(i).le.EGRPMAX) then
      ig = ig + 1
      EGROUP(ig) = SAND(i)
      endif
      enddo
      NGROUP = ig - 1 ! # of groups is one less than # group boundaries
      return
      END
      SUBROUTINE UKAEAGEN(EGROUP,NGROUP,EGRPMAX)
C===============================================================
C
C     DEFINE UKAEA GROUP STRUCTURE FROM 1.0D-5 EV TO
C     ANY UPPER LIMIT ABOVE 30 MEV DEFINED BY EGRPMAX INPUT.
C
C     WARNING - As used by GROUPIE   EGRPMAX = 1 GeV.
C               There is no test for EGRPMAX > 30 MeV (Caveat Emptor).
C
C===============================================================
      INCLUDE 'implicit.h'
      DIMENSION EGROUP(*)
      DATA TEN/  10.0D+00/
      DATA FIFTY/50.0D+00/
      DATA EGRPMIN1/1.0D-05/
      DATA EGRPMIN2/5.5D-01/
      DATA EGRPMIN3/1.0D+01/
      DATA EGRPMIN4/5.0D+06/
      DATA E10MEV  /1.0D+07/
      DATA E30MEV  /3.0D+07/
C
C     GROUPS SPLIT INTO EQUAL E (energy) AND U (lethargy) SECTIONS
C      237 EQUAL U FROM 1.0E-05 TO 5.5D-01
C      378 EQUAL E FROM 5.5D-01 TO 1.0D+01
C      285 EQUAL U FROM 1.0D+01 TO 5.0D+06
C      126 EQUAL E FROM 5.0D+06 TO ~3.02D+07
C       76 EQUAL U FROM ~3.02D+07 TO 1.0D+09
C     1102 TOTAL GROUPS
C
C-----CUMULATIVE GROUPS USED IN BUILD BELOW
      NGROUP1= 237
      NGROUP2= 615
      NGROUP3= 900
      NGROUP4=1026
      NGROUP5=1102
      NGROUP= 1102
C-----50 GROUPS PER ENERGY DECADE (Same as TART Group Structure).
      DEU=DEXP(DLOG(TEN)/FIFTY)
C-----------------------------------------------------------------------
C
C     DEFINE GROUPS SEQUENTIALLY OVER 5 LOOPS
C     OPTION TO EXIT AT OTHER MAXIMUM ENERGY ABOVE 30 MEV
C
C-----------------------------------------------------------------------
c
c     10E-5 eV to 0.55 eV: 50 Groups per energy decade (same as TART)
c
      ENOW=EGRPMIN1
      DO I=1,NGROUP1
      EGROUP(I)=ENOW
      ENOW=ENOW*DEU
      ENDDO
c
c     0.55 eV to 10 eV: 0.025 eV spacing.
c
      ENOW=EGRPMIN2
      DO I=NGROUP1+1,NGROUP2
      EGROUP(I)=ENOW
      ENOW=ENOW+2.5D-02
      ENDDO
c
c     10 eV to 5 MeV: 50 Groups per energy decade (same as TART)
c
      ENOW=EGRPMIN3
      DO I=NGROUP2+1,NGROUP3
      EGROUP(I)=ENOW
      ENOW=ENOW*DEU
      ENDDO
c
c     5 MeV to 30 MeV: 200 keV spacing.
c
      ENOW=EGRPMIN4
      DO I=NGROUP3+1,NGROUP4
      EGROUP(I)=ENOW
      ENOW=ENOW+2.0D+05
      ENDDO
c
c     30 MeV to 1 GeV: 50 Groups per energy decade (same as TART)
c
c-----First TART group boundary over 30 MeV
      EGRPMIN5 = E10MEV
   10 EGRPMIN5 = EGRPMIN5*DEU
      if(EGRPMIN5.le.E30MEV) go to 10
      ENOW=EGRPMIN5
      DO I=NGROUP4+1,NGROUP5+1
      EGROUP(I)=ENOW
      IF(EGROUP(I).gt.EGRPMAX) go to 20
      IF(EGROUP(I).eq.EGRPMAX) go to 30
      ENOW=ENOW*DEU
      ENDDO
      I = NGROUP5+1
C-----DEFINE END OF GROUPS AT MAXIMUM ENDL ENERGY
   20 EGROUP(I)=EGRPMAX
   30 NGROUP=I-1
      RETURN
      END
      SUBROUTINE BANOUT
C=======================================================================
C
C     OUTPUT MULTI-BAND PARAMETERS FOR NEXT MATERIAL. PACK TABLE
C     TO OUTPUT ONLY THE NUMBER OF BANDS REQUIRED PER GROUP AND WORDS
C     PER BAND AND THEN OUTPUT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OTAPE2,TYPSIG
      CHARACTER*1 ZABCD,FIELDX
      CHARACTER*72 LIBID
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
CAK   COMMON/TEMPO/TMPTAB(3),TEMP1,NTEMP
      COMMON/TEMPOg/TMPTAB(3),TEMP1,NTEMP
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
      COMMON/GROUPC/LIBID
      COMMON/ELPASZ/ZABCD(10)
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      INCLUDE 'groupie.h'
C
C     ON FIRST CALL OUTPUT LIBRARY I.D.
C
      DATA IPASS/0/
      IF(IPASS.GT.0) GO TO 10
      IPASS=1
      WRITE(OTAPE2,50) LIBID
C
C     OUTPUT MULTIBAND PARAMETERS.
C
C-----OUTPUT,
C-----1) ZA
C-----2) NUMBER OF GROUPS
C-----3) NUMBER OF BANDS
C-----4) TEMPERATURE (KELVIN)
C-----5) ZA I.D. IN HOLLERITH.
C
C     ONLY OUTPUT THE NUMBER OF BANDS NEEDED, BUT ALWAYS
C     AT LEAST 2.
C
   10 NBOUT=NBNEED
      IF(NBOUT.LT.2) NBOUT=2
      CALL OUT9G(TEMP1,FIELDX(1,1))
      WRITE(OTAPE2,60) IZA,NGR,NBOUT,(FIELDX(M,1),M=1,11),ZABCD
C
C     OUTPUT PARAMETERS FOR EACH GROUP/BAND
C
C-----LOOP OVER GROUPS.
      DO 40 IG=1,NGR
C-----FORMAT LOWER GROUP ENERGY BOUNDARY FOR OUTPUT.
      CALL OUT9G(EGROUP(IG)     ,FIELDX(1,   1))
C-----LOOP OVER BANDS
      DO 30 IB=1,NBOUT
C-----FORMAT GROUP/BAND WEIGHT FOR OUTPUT.
      CALL OUT9G(WTBAND(1,IB,IG),FIELDX(1,   2))
C-----FORMAT GROUP/BAND CROSS SECTIONS FOR OUTPUT.
      IOUT=2
      DO 20 IS=2,NSECT
      IOUT=IOUT+1
   20 CALL OUT9G(XCBAND(IS,IB,IG),FIELDX(1,IOUT))
C-----ENERGY OUTPUT WITH FIRST BAND
      IF(IB.EQ.1) WRITE(OTAPE2,70) ((FIELDX(M,K),M=1,11),K=1,IOUT)
C-----NO ENERGY OUTPUT WITH OTHER BANDS
      IF(IB.NE.1) WRITE(OTAPE2,80) ((FIELDX(M,K),M=1,11),K=2,IOUT)
   30 CONTINUE
   40 CONTINUE
C-----OUTPUT LAST GROUP BOUNDARY.
      CALL OUT9G(EGROUP(NGRP1),FIELDX(1,1))
      WRITE(OTAPE2,70) (FIELDX(M,1),M=1,11)
      RETURN
   50 FORMAT(A72)
   60 FORMAT(3I11,11A1,1X,10A1)
   70 FORMAT(11A1,11A1,55A1)
   80 FORMAT(11X ,11A1,55A1)
      END
      SUBROUTINE OVERE(ESPECT,SPECT,N,ELOW,EHIGH)
C=======================================================================
C
C     CONSTRUCT TABLE OF 1/E AS DEFINED BY A TABLE OF N POINTS
C     LOGARITHMICALLY SPACED BETWEEN THE ENERGIES ELOW AND EHIGH.
C     VALUES ARE RETURNED IN ASCENDING ENERGY ORDER.
C
C     THE MAXIMUM PER-CENT DIFFERENCE BETWEEN THE EXACT 1/E SHAPE
C     (MAGNITUDE DOES NOT MATTER) AND THE TABULATED APPROXIMATION
C     (ASSUMING LINEAR VARIATION BETWEEN TABULATED POINTS) DEPENDS
C     ONLY ON THE NUMBER OF POINTS PER ENERGY DECADE. THE MAXIMUM
C     ERROR VERSUS POINTS PER DECADE WILL BE...
C
C     POINTS PER DECADE    MAXIMUM PER-CENT ERROR
C     -----------------    ----------------------
C                 10           0.67
C                 20           0.17
C                 30           0.074
C                 40           0.042
C                 50           0.027
C                 60           0.019
C                 70           0.014
C                 80           0.011
C                 90           0.009
C                100           0.007
C
C     FOR EXAMPLE IN ORDER TO CONSTRUCT A TABULATED LINEARLY
C     INTERPOLABLE APPROXIMATION TO 1/E FROM 1.0D-5 EV UP TO 10 MEV
C     (12 DECADES OF ENERGY) TO WITHIN AN ACCURACY OF 0.1 PER-CENT
C     CAN BE ACCOMPLISHED BY USING 30 POINTS PER DECADE, OR A TOTAL
C     OF 361 POINTS (12 X 30 + 1 AT END OF ENERGY INTERVAL). THIS
C     WOULD YIELD A MAXIMUM ERROR OF 0.074 PER-CENT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION ESPECT(*),SPECT(*)
C-----DEFINE MULTIPLIER TO SPACE N POINTS BETWEEN ELOW AND EHIGH
C-----USING EQUAL LOG (I.E., LETHARGY) SPACING.
      FN=N-1
      DE=DEXP(DLOG(EHIGH/ELOW)/FN)
C-----DEFINE FIRST POINT AT LOWER ENERGY LIMIT.
      ESPECT(1)=ELOW
      SPECT(1)=1.0/ELOW
C-----DEFINE REMAINING POINTS AS LOGARITHMICALLY SPACED, I.E. EACH
C-----ENERGY POINT IS A MULTILPE OF THE PRECEEDING ENERGY POINT.
      DO 10 J=2,N
      ESPECT(J)=DE*ESPECT(J-1)
   10 SPECT(J)=1.0/ESPECT(J)
C-----ELIMINATE ROUNDOFF PROBLEMS BY DEFINING LAST POINT AS EXACTLY
C-----AT UPPER ENERGY LIMIT.
      ESPECT(N)=EHIGH
      SPECT(N)=1.0/EHIGH
      RETURN
      END
      SUBROUTINE SPECTM(ESPECT,SPECT,N,ETAB,TEMPM,WA,WB)
C=======================================================================
C
C     CONSTRUCT A TABULATED WEIGHTING SPECTRUM COMPOSED OF,
C     (1) MAXWELLIAN - FROM ETAB(1) TO ETAB(2)
C                      S(E) = E*EXP(-E/TEMPM)/TEMPM
C     (2) 1/E        - FROM ETAB(2) TO ETAB(3)
C                      S(E) = C1/E
C     (3) FISSION    - FROM ETAB(3) TO ETAB(4)
C                      S(E) = C2*EXP(-E/WA)*SINH(SQRT(E*WB))
C     (4) FLAT       - FROM ETAB(4) TO ETAB(5)
C                      S(E) = C3 - HIGH ENERGY FLAT.
C
C     NOTE, THIS IS A FLUX WEIGHTING SPECTRUM - NOT NEUTRON
C     DENSITY. A MAXWELLIAN DEFINES THE NEUTRON DENSITY,
C
C     N(E) = C*SQRT(E)*DEXP(-B*E)*D(E)
C
C     FOR THE CORRECT FLUX WEIGHTING MULTIPLY THIS BY NEUTRON
C     SPEED, WHICH INTRODUCES AN ADDITIONAL FACTOR OF SQRT(E).
C     SO THAT THE CORRECT THERMAL WEIGHTING IS,
C
C     FLUX(E) = C*E*EXP(-B*E)*D(E)
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION ESPECT(*),SPECT(*),ETAB(5)
C-----DEFINE ALLOWABLE ERROR.
C-----05/24/09 - DECREASED FROM 1.0D-5 TO 1.0D-6
      DATA ERROKS/1.0D-6/
      DATA ONE/1.0D+0/
      DATA HALF/5.0D-1/
C-----INITIALIZE COEFFICIENTS TO MAKE SPECTRUM CONTINUOUS.
      C1=1.0
      C2=1.0
C-----START WITH A MINIMUM OF 50 POINTS IN EACH PART OF SPECTRUM.
      I=50
      XN=I-1
C-----INITIALIZE NUMBER OF POINTS IN SPECTRUM.
      N=0
C
C     CONSTRUCT MAXWELLIAN.
C
      IF(ETAB(2).LE.ETAB(1)) GO TO 50
      DU=DLOG(ETAB(2)/ETAB(1))/XN
      DE=DEXP(DU)
      ENOW=ETAB(1)
   10 N=N+1
      ESPECT(N)=ENOW
      SPECT(N)=ENOW*DEXP(-ENOW/TEMPM)/TEMPM
      IF(N.EQ.1) GO TO 30
   20 EMID=HALF*(ESPECT(N)+ESPECT(N-1))
      SMID=EMID*DEXP(-EMID/TEMPM)/TEMPM
      STERP=HALF*(SPECT(N)+SPECT(N-1))
      IF(ABS(STERP-SMID).LE.ERROKS*ABS(SMID)) GO TO 30
      ESPECT(N)=EMID
      SPECT(N)=SMID
      GO TO 20
   30 IF(ESPECT(N).GE.ETAB(2)) GO TO 40
      ENOW=ESPECT(N)*DE
      IF(ENOW.GE.ETAB(2)) ENOW=ETAB(2)
      GO TO 10
C
C     CONSTRUCT 1/E
C
   40 C1=ESPECT(N)*SPECT(N)
   50 IF(ETAB(3).LE.ETAB(2)) GO TO 100
      DU=DLOG(ETAB(3)/ETAB(2))/XN
      DE=DEXP(DU)
      ENOW=ESPECT(N)*DE
   60 N=N+1
      ESPECT(N)=ENOW
      SPECT(N)=C1/ENOW
      IF(N.EQ.1) GO TO 80
   70 EMID=HALF*(ESPECT(N)+ESPECT(N-1))
      SMID=C1/EMID
      STERP=HALF*(SPECT(N)+SPECT(N-1))
      IF(ABS(STERP-SMID).LE.ERROKS*ABS(SMID)) GO TO 80
      ESPECT(N)=EMID
      SPECT(N)=SMID
      GO TO 70
   80 IF(ESPECT(N).GE.ETAB(3)) GO TO 90
      ENOW=ESPECT(N)*DE
      IF(ENOW.GE.ETAB(3)) ENOW=ETAB(3)
      GO TO 60
C
C     CONSTRUCT FISSION SPECTRUM.
C
   90 R=DSQRT(ESPECT(N)*WB)
      EXPR=DEXP(R)
      C2=SPECT(N)/(DEXP(-ESPECT(N)/WA)*(EXPR-ONE/EXPR))
  100 IF(ETAB(4).LE.ETAB(3)) GO TO 140
      DU=DLOG(ETAB(4)/ETAB(3))/XN
      DE=DEXP(DU)
      ENOW=ESPECT(N)*DE
  110 N=N+1
      ESPECT(N)=ENOW
      R=DSQRT(ENOW*WB)
      EXPR=DEXP(R)
      SPECT(N)=C2*(DEXP(-ENOW/WA)*(EXPR-ONE/EXPR))
      IF(N.EQ.1) GO TO 130
  120 EMID=HALF*(ESPECT(N)+ESPECT(N-1))
      R=DSQRT(EMID*WB)
      EXPR=DEXP(R)
      SMID=C2*(DEXP(-EMID/WA)*(EXPR-ONE/EXPR))
      STERP=HALF*(SPECT(N)+SPECT(N-1))
      IF(ABS(STERP-SMID).LE.ERROKS*ABS(SMID)) GO TO 130
      ESPECT(N)=EMID
      SPECT(N)=SMID
      GO TO 120
  130 IF(ESPECT(N).GE.ETAB(4)) GO TO 140
      ENOW=ESPECT(N)*DE
      IF(ENOW.GE.ETAB(4)) ENOW=ETAB(4)
      GO TO 110
C
C     CONSTRUCT CONSTANT
C
  140 N=N+1
      ESPECT(N)=ETAB(5)
      SPECT(N)=SPECT(N-1)
      RETURN
      END
CAK   SUBROUTINE READIN
      SUBROUTINE READINg(infile,outfile)
C=======================================================================
C
C     INITIALIZE PARAMETERS AND THEN READ AND CHECK INPUT DATA.
C     INITIALIZE ALL OUTPUT FILES THAT WILL BE USED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,OTAPE2,OPS,TYPSIG
      CHARACTER*1 FIELDX
      CHARACTER*4 MESS1,NO,YES,OPSC
      CHARACTER*16 MYLAW,MYGRUP
      CHARACTER*24 TYPEL
      CHARACTER*36 HWEIGH
      CHARACTER*28 MESSBAND
      CHARACTER*72 LIBID
CAK
      CHARACTER*72 infile,outfile
CAK
      CHARACTER*8 SIGHL1,SIGHL2
      CHARACTER*72 NAMEIN,NAMEOUT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/IOSTATUS/ISTAT1,ISTAT2
      COMMON/REALLY/XE,XA(5),YA(5),YB(5),YC(5),NPTAB(5),IPTAB(5),NSECT
      COMMON/HEADER/ZA,AWRIN,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LAWYER/MYLAWO
C     COMMON/TEMPO/TMPTAB(3),TEMP1,NTEMPK
      COMMON/TEMPOg/TMPTAB(3),TEMP1,NTEMP
      COMMON/TEMPC/TINOUT(4)
      COMMON/TRASH/ERROK,ER10PC,MGROUP,METHODB,NBAND,NBMAX,
     1 NBNEED,TYPSIG
      COMMON/COMXTEND/MYXTEND
      COMMON/GROUPC/LIBID
CAK   COMMON/MATZA/MODGET,NMATZA,MATMIN(101),MFMIN(101),MTMIN(101),
CAK  1 MATMAX(101),MFMAX(101),MTMAX(101)
      COMMON/MATZAg/MODGET,NMATZA,MATMIN(101),MFMIN(101),MTMIN(101),
     1 MATMAX(101),MFMAX(101),MTMAX(101)
      COMMON/GRABER/NOSELF,LSECT
      COMMON/SPIM/IMSP
      COMMON/OPUS/OPS(5)
      COMMON/NAMEX/NAMEIN,NAMEOUT
      COMMON/INPAGE/IXYLOW(5),IXYHI(5),ISCR(5)
CAK   COMMON/FIELDC/FIELDX(11,22)
      COMMON/FIELDCg/FIELDX(11,22)
      COMMON/SIGGIE/SIGHL1(12),SIGHL2(12)
      INCLUDE 'groupie.h'
      DIMENSION MYLAW(3),MYGRUP(21),TYPEL(3),HWEIGH(4),MESSBAND(3),
     1 MESS1(2),MYGRUI(21),OPSC(5),ETAB(5)
C-----DEFINE DESCRIPTION FOR ENDF/B FORMATTED OUTPUT INTERPOLATION LAW.
      DATA MYLAW/'                ',
     1           '(Histogram)     ',
     2           '(Linear Linear) '/
      DATA MESSBAND/
     1 '                            ',
     2 'Conserve 1/Total^2          ',
     2 'Conserve 1/[Total+<Total>]  '/
c       1234567890123456789012345678901234567890
C-----DEFINE TWO POSSIBLE SELECTION CRITERIA.
      DATA MESS1/' MAT','  ZA'/
C-----DEFINE TWO STATES FOR EACH OUTPUT DEVICE.
      DATA NO/'  No'/
      DATA YES/' Yes'/
C-----NONE, CROSS SECTION OR RESONANCE INTEGRAL OUTPUT.
      DATA TYPEL/
     1 '                        ',
     1 '(Cross Sections)        ',
     2 '(Resonance Integrals)   '/
C-----DEFINE DESCRIPTION OF BUILT IN GROUP STRUCTURES.
      DATA MYGRUP/
     1 '(From Input)    ','(TART)          ',
     2 '(ORNL)          ','(ORNL)          ',
     3 '(ORNL)          ','(SAND-II)       ',
     4 '(SAND-II EXTEND)','(WIMS)          ',
     5 '(GAM-I)         ','(GAM-II)        ',
     6 '(MUFT)          ','(ABBN)          ',
     7 '(TART)          ','(TART)          ',
     8 '(SAND-II)       ','(SAND-II EXTEND)',
     9 '(TART)          ','(SAND-II 60 MeV)',
     A '(SAND-II 150MeV)','(SAND-II 200MeV)',
     B '(UKAEA    1 GeV)'/
      DATA MYGRUI/
     1    0,175, 50,126,171,620,640, 69, 68, 99,
     2   74, 28,616,700,665,685,666,725,755,765,
     3 1102/
C-----DEFINE DESCRIPTION OF WEIGHTING SPECTRUM.
      DATA HWEIGH/
     1 '(Maxwellian-1/E-Fission Spectrum)   ',
     2 '(1/E)                               ',
     3 '(Flat Weighted)                     ',
     4 '(From Input)                        '/
C-----DEFINE MINIMUM ALLOWABLE CONVERGENCE CRITERIA.
      DATA ERRMIN/1.0D-04/
C-----DEFINE STANDARD ALLOWABLE CONVERGENCE CRITERIA.
      DATA ERRUSE/1.0D-03/
C-----INITIALIZE CURRENT ENDF/B MAT, MF AND MT.
      MATH=0
      MFH=0
      MTH=0
C-----INITIALIZE ERROR IN INPUT DATA FLAG OFF.
      KERR=0
C-----DEFINE LIMITS FOR SELECTION OF OUTPUT TEMPERATURE.
      BOLTZM=8.6348D-05
      TMPTAB(1)=0.9999D+00/BOLTZM
      TMPTAB(2)=9.999D+02/BOLTZM
      TMPTAB(3)=9.999D+05/BOLTZM
C-----DEFINE CONVERSION FACTORS FROM KELVIN TO OUTPUT UNITS
C-----(KELVIN, EV, KEV OR MEV).
      TINOUT(1)=1.0D+0
      TINOUT(2)=BOLTZM
      TINOUT(3)=BOLTZM/1.0D+03
      TINOUT(4)=BOLTZM/1.0D+06
C-----INITIALIZE EVALUATION COUNT FOR OUTPUT LISTING.
      LZA=0
      WRITE(OUTP,540)
      WRITE(*   ,540)
C
C     READ ALL INPUT PARAMETERS.
C
CAK   IF(ISTAT1.EQ.1) GO TO 20
c-----2017/5/6 - Changed floating point to character.
CAK   READ(INP,10,END=20,ERR=20)
CAK  1 MODGET,MGROUP,METHODB,MSPECT,(FIELDX(k,1),k=1,11),TYPSIG,MYXTEND
CAK10 FORMAT(4I11,11A1,I11,I4)
CAK   CALL IN9(ERRIN,FIELDX(1,1))
c-----2017/5/6 - Changed floating point to character.
CAK   GO TO 30
C-----DEFINE DEFAULT VALUES
CAK20 ISTAT1    = 1
      MODGET    = 0
CAK   MGROUP    = 0
      MGROUP    = -12
      METHODB   = 2         ! <X> <1/(X+<X>) <1/X>
CAK   MSPECT    = 0
      MSPECT    = -1
CAK   ERRIN     = 0.0
      ERRIN     = 1.e-3
      TYPSIG    = 0
      MYXTEND   = 0
C-----02/02/02 - UPDATE TO ALLOW INPUT SIGMA0 VALUES
C
C     IF TYPSIG = -1, READ 21 VALUES OF FIXED SIGMA0
C
   30 IF(TYPSIG.GE.0) GO TO 50
      TYPSIG = 1
c-----2017/5/6 - Changed floating point to character.
      READ(INP,40)    ((FIELDX(j,i),j=1,11),i=2,22)
   40 FORMAT(66A1)
      do i=2,22
      CALL IN9(SHIELD(i),FIELDX(1,i))
      enddo
c-----2017/5/6 - Changed floating point to character.
      WRITE(OUTP,760) (SHIELD(I),I=2,22)
      WRITE(   *,760) (SHIELD(I),I=2,22)
C-----THEY MUST BE IN DESCENDING ORDER
      MYERR = 0
      DO I=2,22
      IF(SHIELD(I).LE.SHIELD(I+1)) THEN
      WRITE(OUTP,770) I,SHIELD(I),SHIELD(I+1)
      WRITE(   *,770) I,SHIELD(I),SHIELD(I+1)
      MYERR = MYERR + 1
      ENDIF
      ENDDO
      IF(MYERR.NE.0) THEN
      WRITE(OUTP,780)
      WRITE(   *,780)
      CALL ENDERROR
      ENDIF
C-----CHANGE TITLES IDENTIFYING SIGMA0 VALUES
      DO I=2,10
      CALL OUT8(SHIELD(2*I),SIGHL2(I))
      ENDDO
C
C     READ FILENAMES - IF BLANK USE STANDARD FILENAMES
C
C-----INPUT DATA.
CAK50 IF(ISTAT1.EQ.1) GO TO 70
CAK   READ(INP,60,END=70,ERR=70) NAMEIN
CAK60 FORMAT(A72)
   50 NAMEIN=infile
      IF(NAMEIN.EQ.' ') NAMEIN = 'ENDFB.IN'
C-----OUTPUT DATA.
CAK   READ(INP,60,END=80,ERR=80) NAMEOUT
      NAMEOUT=outfile
      IF(NAMEOUT.EQ.' ') NAMEOUT = 'ENDFB.OUT'
      GO TO 90
C-----USE DEFAULT FILENAMES
CAK70 NAMEIN  = 'ENDFB.IN'
   80 NAMEOUT = 'ENDFB.OUT'
CAK   ISTAT1 = 1
C-----PRINT FINAL FILENAMES
   90 WRITE(OUTP,100) NAMEIN,NAMEOUT
      WRITE(*   ,100) NAMEIN,NAMEOUT
  100 FORMAT(
     1 ' ENDF/B Input and Output Data Filenames'/1X,A72/
     2 1X,A72/1X,78('-'))
C
C     READ OUTPUT OPTIONS.
C
CAK   IF(ISTAT1.EQ.1) GO TO 120
c-----2017/5/6 - Changed floating point to character.
CAK   READ(INP,110,END=120,ERR=120) OPS,(FIELDX(j,1),j=1,11)
CAK110 FORMAT(5I11,11A1)
CAK   CALL IN9(TEMPMIN,FIELDX(1,1))
c-----2017/5/6 - Changed floating point to character.
CAK   GO TO 130
C-----DEFINE DEFAULT VALUES
CAK120 ISTAT1 = 1
      OPS(1) = 1
      OPS(2) = 1
      OPS(3) = 1
      OPS(4) = 1
      OPS(5) = 1
      TEMPMIN = 0.0
CAK130 IF(ISTAT1.EQ.1) GO TO 150
CAK   READ(INP,140) LIBID
CAK140 FORMAT(A72)
      LIBID  = ' '
      GO TO 160
  150 ISTAT1 = 1
      LIBID  = ' '
  160 IF(MODGET.NE.0) MODGET=1
C
C     IN PRODUCTION VERSION ONLY ALLOW 0 OR 2 BANDS
C
C-----ONLY ALLOW 2 BANDS
      IF(METHODB.LT.0) METHODB = 0
      IF(METHODB.GT.2) METHODB = 2
C-----IF NOT CROSS SECTIONS OR RESONANCE INTEGRALS SELECTED TURN OFF
      IF(OPS(1).LT.0.OR.OPS(1).GT.2) OPS(1)=0
      IF(OPS(5).LT.0.OR.OPS(5).GT.2) OPS(5)=0
      IMOPS1=OPS(1)+1
      IMOPS5=OPS(5)+1
C-----IF ILLEGAL INTERPOLATION LAW FOR ENDF/B FORMATTED OUTPUT, TURN
C-----OFF ENDF/B FORMATTED OUTPUT.
      IF(IABS(OPS(4)).GT.2) OPS(4)=0
C-----DEFINE MESSAGE TO DESCRIBE ENDF/B FORMATTED INTERPOLATION LAW.
      MYLAWO=OPS(4)
C-----IF MULTI-BAND PARAMETERS NOT CALCULATED OUTPUT IS IMPOSSIBLE.
      IF(METHODB.NE.0) GO TO 170
      OPS(2)=0
      OPS(3)=0
      GO TO 180
C-----IF NO MULTI-BAND OUTPUT DO NOT CALCULATE PARAMETERS.
  170 IF(OPS(2).NE.0.OR.OPS(3).NE.0) GO TO 180
      METHODB=0
  180 NBMAX=METHODB
      IF(NBMAX.EQ.1) NBMAX=2
C
C     DEFINE WEIGHTING SPECTRUM
C
C-----INSURE THAT STANDARD SPECTRUM OPTION IS IN STANDARD FORM.
      IF(MSPECT.GE.0) GO TO 190
C-----DEFINE SPECTRUM TO BE 1/E OR MAXWELLIAN, 1/E AND FISSION SPECTRUM
      IF(MSPECT.LT.-2) MSPECT=-2
c-----05/24/09 - increased from 3,000 to 30,000
      IF(MSPECT.EQ.-1) THEN
      IMSP=2        ! 1/E
      MSPTAB=30000
      ELSE
      IMSP=1        ! Maxwellian + 1/E + Fission
      MSPTAB=16724
      ENDIF
      GO TO 200
C-----DEFINE SPECTRUM TO FLAT WEIGHTED OR READ FROM INPUT.
  190 IF(MSPECT.LE.1) MSPECT=0
      MSPTAB=MSPECT
      IF(MSPTAB.LT.1) MSPTAB=2
      IMSP=3
      IF(MSPECT.GT.0) IMSP=4
C
C     INITIALIZE ERROR TO STANDARD OPTION. ONLY USE INPUT ERROR IF IT
C     IS IN THE LEGAL RANGE (I.E. GREATER THAN 0.0001).
C
  200 ERROK=ERRUSE
      IF(ERRIN.GE.ERRMIN) ERROK=ERRIN
      ERRPC=100.0*ERROK
      ER10PC=0.1*ERROK
C-----INSURE SIGMA-0 DEFINITION SELECTOR IS STANDARD VALUE.
      IF(TYPSIG.LE.0) TYPSIG=0
      IF(TYPSIG.GT.0) TYPSIG=1
C-----CONVERT OUTPUT DEVICE SELECTORS TO HOLLERITH (FOR OUTPUT LISTING)
C-----AND DEFINE ALL UNUSED LOGICAL NUMBERS TO BE ZERO.
      IPICK=0
      DO 270 I=1,5
      IF(IABS(OPS(I)).GT.0) GO TO 260
      OPSC(I)=NO
      GO TO (210,220,230,240,250),I
  210 LIST1=0
C-----INITIALIZE TO NO PLOTTED OUTPUT (SHIELDED/UNSHIELDED)
      IPLOT=0
      GO TO 270
  220 LIST2=0
      GO TO 270
  230 OTAPE2=0
      GO TO 270
  240 OTAPE=0
      GO TO 270
  250 LIST3=0
      GO TO 270
  260 OPSC(I)=YES
      IPICK=1
  270 CONTINUE
C-----OPTIONALLY DEFINE FILENAMES.
CAK   CALL FILIO2
      CALL FILIO2g
C
C     TERMINATE IF ERROR OPENING ENDF/B DATA FILE
C
      IF(ISTAT2.EQ.1) THEN
      WRITE(OUTP,280) NAMEIN
      WRITE(   *,280) NAMEIN
  280 FORMAT(//' ERROR - Opening ENDF/B data file'/1X,A72//)
      CALL ENDERROR
      ENDIF
C
C     IF NO SELF-SHIELDED OR MULTI-BAND OUTPUT SET FLAG TO BYPASS
C     THESE CALCULATIONS (I.E. ONLY DO UNSHIELDED MULTI-GROUP
C     CALCULATIONS).
C
      IF(OTAPE2.NE.0.OR.LIST1.NE.0.OR.LIST2.NE.0) GO TO 290
C-----NO SELF-SHIELDING CALCULATION. UP TO 20,000 GROUPS ALLOWED.
C-----CHECK FOR LEGAL NUMBER OF GROUPS.
      NOSELF=1
C-----LOAD ALL CROSS SECTIONS INTO SECOND PAGE.
      LSECT=2
      IF(MGROUP.GE.-19.AND.MGROUP.LE.MAXGROUP) GO TO 300
      WRITE(OUTP,700) MGROUP,MAXGROUP
      GO TO 530
C-----SELF-SHIELDING CALCULATION. UP TO 20,000 GROUPS ALLOWED.
C-----CHECK FOR LEGAL NUMBER OF GROUPS.
  290 NOSELF=0
C-----LOAD ALL OTHER (OTHER= NOT TOTAL, ELASTIC, CAPTURE OR FISSION)
C-----CROSS SECTIONS INTO PAGE 4 (THIS PAGE WILL BE LAST ONE USED
C-----WHEN TOTAL, ELASTIC, CAPTURE AND FISSION ARE READ.
      LSECT=4
      IF(MGROUP.GE.-19.AND.MGROUP.LE.MAXGROUP) GO TO 300
      WRITE(OUTP,690) MGROUP,MAXGROUP
      GO TO 530
  300 IMGR=IABS(MGROUP)+2
      IF(MGROUP.GT.0) IMGR=1
      MYGRUI(1)=MGROUP
C
C     LIST INTERPRETATION OF INPUT PARAMETERS.
C
      CALL OUT9G(ERROK,FIELDX(1,1))
      WRITE(OUTP,550) MESS1(MODGET+1),MYGRUI(IMGR),
     1 MYGRUP(IMGR),NBMAX,MESSBAND(METHODB+1),MSPTAB,HWEIGH(IMSP),
     2 (FIELDX(M,1),M=1,11),ERRPC
      WRITE(*   ,550) MESS1(MODGET+1),MYGRUI(IMGR),
     1 MYGRUP(IMGR),NBMAX,MESSBAND(METHODB+1),MSPTAB,HWEIGH(IMSP),
     2 (FIELDX(M,1),M=1,11),ERRPC
      IF(TYPSIG.EQ.0) THEN
      WRITE(OUTP,640)
      WRITE(*   ,640)
      ELSE
      WRITE(OUTP,650)
      WRITE(*   ,650)
      ENDIF
      IF(MYXTEND.EQ.0) THEN   ! high energy cross section extension.
      WRITE(OUTP,660)
      WRITE(*   ,660)
      ELSE
      MYXTEND = 1
      WRITE(OUTP,670)
      WRITE(*   ,670)
      ENDIF
C-----MUST USE MULTIPLES OF UNSHIELDED TOTAL IF CALCULATING
C-----ULTI-BAND PARAMETERS
      IF(TYPSIG.NE.0.AND.NBMAX.GT.0) THEN
      WRITE(OUTP,310)
      WRITE(*   ,310)
  310 FORMAT(1X,78('-')/
     1 ' ERROR - Multi-Band Parameters Require that Sigma-0'/
     2 '         be Multiples of the Unshielded Total. You MUST'/
     3 '         Turn off either Multi-bands or Sigma-0 Definition.'/
     3 '         EXECUTION TERMINATED.'///)
      CALL ENDERROR
      ENDIF
      WRITE(OUTP,560) OPSC(1),TYPEL(IMOPS1),
     1 (OPSC(J),J=2,4),MYLAW(MYLAWO+1),
     2 OPSC(5),TYPEL(IMOPS5)
      WRITE(*   ,560) OPSC(1),TYPEL(IMOPS1),
     1 (OPSC(J),J=2,4),MYLAW(MYLAWO+1),
     2 OPSC(5),TYPEL(IMOPS5)
C-----IF USING STANDARD BUILT-IN SPECTRA PRINT MAXWELLIAN TEMPERATURE
      IF(MSPECT.EQ.-2) THEN
      IF(TEMPMIN.GT.1000.0D+0) TEMPMIN = 1000.0D+0
      IF(TEMPMIN.LE.0.0) TEMPMIN = 0.0253 ! DEFAULT ROOM TEMPERATURE
      WRITE(OUTP,570) TEMPMIN
      WRITE(   *,570) TEMPMIN
      ENDIF
      WRITE(OUTP,580) LIBID
      WRITE(   *,580) LIBID
C-----IF NO OUTPUT DEVICES SELECTED TERMINATE RUN.
      IF(IPICK.NE.0) GO TO 320
      WRITE(OUTP,720)
      GO TO 530
C
C     READ SELECTION RANGES (EITHER MAT OR ZA). IF MAXIMUM IS LESS THAN
C     MINIMUM SET IT EQUAL TO MINIMUM. IF FIRST CARD IS BLANK RETRIEVE
C     ALL DATA.
C
  320 IF(MODGET.EQ.0) WRITE(OUTP,590)
      IF(MODGET.EQ.1) WRITE(OUTP,600)
      IF(MODGET.EQ.0) WRITE(*   ,590)
      IF(MODGET.EQ.1) WRITE(*   ,600)
CAK   IF(ISTAT1.EQ.1) GO TO 340
CAK   READ(INP,330,END=340,ERR=340)
CAK  1 MATMIN(1),MFMIN(1),MTMIN(1),MATMAX(1),MFMAX(1),MTMAX(1)
  330 FORMAT(I6,I2,I3,I6,I2,I3)
      GO TO 350
C-----DEFINE DEFAULT VALUES
  340 ISTAT1    = 1
      MATMIN(1) = 0
      MFMIN(1)  = 0
      MTMIN(1)  = 0
      MATMAX(1) = 0
      MFMAX(1)  = 0
      MTMAX(1)  = 0
  350 IF(MATMIN(1).GT.0.OR.MATMAX(1).GT.0) GO TO 360
      MATMAX(1)=9999
      IF(MFMIN(1).LT.0) MFMIN(1)=0
      IF(MFMAX(1).LE.0) MFMAX(1)=99
      IF(MTMIN(1).LT.0) MTMIN(1)=0
      IF(MTMAX(1).LE.0) MTMAX(1)=999
      MODGET=0
      NMATZA=2
      WRITE(OUTP,620) MATMIN(1),MFMIN(1),MTMIN(1),MATMAX(1),MFMAX(1),
     1 MTMAX(1)
      WRITE(*   ,620) MATMIN(1),MFMIN(1),MTMIN(1),MATMAX(1),MFMAX(1),
     1 MTMAX(1)
      GO TO 390
  360 IF(MATMAX(1).LT.MATMIN(1)) MATMAX(1)=MATMIN(1)
      IF(MFMIN(1).LT.0) MFMIN(1)=0
      IF(MFMAX(1).LE.0) MFMAX(1)=99
      IF(MTMIN(1).LT.0) MTMIN(1)=0
      IF(MTMAX(1).LE.0) MTMAX(1)=999
      WRITE(OUTP,610) MATMIN(1),MFMIN(1),MTMIN(1),MATMAX(1),MFMAX(1),
     1 MTMAX(1)
      WRITE(*   ,610) MATMIN(1),MFMIN(1),MTMIN(1),MATMAX(1),MFMAX(1),
     1 MTMAX(1)
      DO 370 NMATZA=2,101
      IF(ISTAT1.EQ.1) GO TO 390
      READ(INP,330,END=380,ERR=380)
     1 MATMIN(NMATZA),MFMIN(NMATZA),MTMIN(NMATZA),
     2 MATMAX(NMATZA),MFMAX(NMATZA),MTMAX(NMATZA)
      IF(MATMIN(NMATZA).LE.0.AND.MATMAX(NMATZA).LE.0) GO TO 390
      IF(MATMAX(NMATZA).LT.MATMIN(NMATZA)) MATMAX(NMATZA)=MATMIN(NMATZA)
      IF(MFMIN(NMATZA).LT.0) MFMIN(NMATZA)=0
      IF(MFMAX(NMATZA).LE.0) MFMAX(NMATZA)=99
      IF(MTMIN(NMATZA).LT.0) MTMIN(NMATZA)=0
      IF(MTMAX(NMATZA).LE.0) MTMAX(NMATZA)=999
      WRITE(OUTP,610) MATMIN(NMATZA),MFMIN(NMATZA),MTMIN(NMATZA),
     1 MATMAX(NMATZA),MFMAX(NMATZA),MTMAX(NMATZA)
  370 WRITE(*   ,610) MATMIN(NMATZA),MFMIN(NMATZA),MTMIN(NMATZA),
     1 MATMAX(NMATZA),MFMAX(NMATZA),MTMAX(NMATZA)
      WRITE(OUTP,630)
      WRITE(*   ,630)
      GO TO 530
  380 ISTAT1 = 1
  390 NMATZA=NMATZA-1
C
C     DEFINE GROUP STRUCTURE.
C
      IF(MGROUP.eq.0) go to 400
      IF(MGROUP.gt.0) go to 410
C
C     SELECT BUILT IN GROUP STRUCTURE.
C
C             = 0  175 GROUPS (TART)
C             = 1   50 GROUPS (ORNL)
C             = 2  126 GROUPS (ORNL)
C             = 3  171 GROUPS (ORNL)
C             = 4  620 GROUPS (SAND-II, UP TO 18 MEV)
C             = 5  640 GROUPS (SAND-II, UP TO 20 MEV)
C             = 6   69 GROUPS (WIMS)
C             = 7   68 GROUPS (GAM-I)
C             = 8   99 GROUPS (GAM-II)
C             = 9   54 GROUPS (MUFT)
C             =10   28 GROUPS (ABBN)
C             =11  616 GROUPS (TART TO 20 MEV)
C             =12  700 GROUPS (TART TO 1 GEV)
C             =13  665 GROUPS (SAND-II, 1.0D-5 eV, UP TO 18 MEV)
C             =14  685 GROUPS (SAND-II, 1.0D-5 eV, UP TO 20 MEV)
C             =15  666 GROUPS (TART TO 200 MEV)
C             =16  725 GROUPS (SAND-II, 1.0D-5 eV, UP TO  60 MEV)
C             =17  755 GROUPS (SAND-II, 1.0D-5 eV, UP TO 150 MEV)
C             =18  765 GROUPS (SAND-II, 1.0D-5 eV, UP TO 200 MEV)
C             =19 1102 GROUPS (UKAEA  , 1.0D-5 eV, UP TO   1 GeV)
C
      MGROUP=-MGROUP
  400 CALL GROPE(MGROUP,EGROUP,NGR)
      NGRP1=NGR+1
      GO TO 440
C-----READ AND PRINT ARBITRARY GROUP STRUCTURE.
  410 NGR=MGROUP
      NGRP1=NGR+1
c-----2017/5/6 - Changed floating point to character.
      CALL LISTIV(INP,EGROUP,NGRP1)
c-----2017/5/6 - Changed floating point to character.
      WRITE(OUTP,740)
      DO 430 I=1,NGRP1,6
      II=I+5
      IF(II.GT.NGRP1) II=NGRP1
      J=0
      DO 420 III=I,II
      J=J+1
  420 CALL OUT9G(EGROUP(III),FIELDX(1,J))
  430 WRITE(OUTP,750) ((FIELDX(M,JJ),M=1,11),JJ=1,J)
C-----INSURE GROUP BOUNDARIES ARE IN ASCENDING ENERGY ORDER.
  440 DO 450 I=1,NGR
      IF(EGROUP(I).LT.EGROUP(I+1)) GO TO 450
      KERR=1
      CALL OUT9G(EGROUP(I)  ,FIELDX(1,1))
      CALL OUT9G(EGROUP(I+1),FIELDX(1,2))
      WRITE(OUTP,680) I,((FIELDX(M,KKK),M=1,11),KKK=1,2)
  450 CONTINUE
      IF(KERR.gt.0) go to 530
C-----------------------------------------------------------------------
C
C     DEFINE ENERGY DEPENDENT WEIGHTING SPECTRUM. USE EITHER LINEARLY
C     INTERPOLABLE APPROXIMATION TO 1/E, CONSTANT OR READ ABRITRARY
C     WEIGHTING FUNCTION.
C
C-----------------------------------------------------------------------
      IF(MSPECT.eq.0) go to 470
      IF(MSPECT.gt.0) go to 480
C
C     DEFINE LINEARLY INTERPOLABLE TABULATED SPECTRUM EITHER,
C     (1) 1/E, OR
C     (2) MAXWELLIAN, 1/E, FISSION SPECTRUM
C
      ELOW=EGROUP(1)
      EHIGH=EGROUP(NGRP1)
C
C     2013/7/7 - Corrected Spectrum range to include ALL parts defined
C                below in OVERE and SPECTM.
C
      IF(ELOW .gt.1.0D-05) ELOW  = 1.0D-05   ! 1.0D-5 eV
      IF(EHIGH.lt.2.0D+07) EHIGH = 2.0D+07   ! 20 MeV
C
C     SELECT 1/E OR MAXWELLIAN, 1/E AND FISSION SPECTRUM.
C
      IF(MSPECT.EQ.-2) GO TO 460
C-----CONSTRUCT LINEARLY INTERPOLABLE APPROXIMATION TO 1/E.
      NPTAB(1)=MSPTAB
      IXYLOW(1)=0
      IXYHI(1)=MSPTAB
      CALL OVERE(XPAGE(1,1),YPAGE(1,1),NPTAB(1),ELOW,EHIGH)
      GO TO 490
C
C     CONSTRUCT LINEARLY INTERPOLABLE APPROXIMATION TO MAXWELLIAN, 1/E
C     FISSION AND CONSTANT SPECTRUM.
C
  460 TEMPM=TEMPMIN
      WA=9.65D+5
      WB=2.29D-6
      ETAB(1)=ELOW
      ETAB(2)=4.0*TEMPM
      ETAB(3)=6.7D+04
      ETAB(4)=1.0D+07
      ETAB(5)=EHIGH      ! Note - EHIGH is at least 20 MeV.
      NPTAB(1)=MSPTAB
      CALL SPECTM(XPAGE(1,1),YPAGE(1,1),NPTAB(1),ETAB,TEMPM,WA,WB)
      IXYLOW(1)=0
      IXYHI(1)=NPTAB(1)
      GO TO 490
C
C     DEFINE CONSTANT WEIGHTING SPECTRUM THAT SPANS THE ENERGY GROUP
C
  470 NPTAB(1)=2
      IXYLOW(1)=0
      IXYHI(1)=2
      YPAGE(1,1)=1.0
      YPAGE(2,1)=1.0
      XPAGE(1,1)=EGROUP(1)
      XPAGE(2,1)=EGROUP(NGRP1)
      GO TO 490
C
C     LOAD ARBITRARY TABULATED SPECTRUM INTO PAGING SYSTEM.
C
  480 NPTAB(1)=MSPECT
      CALL PAGIN5(INP,1,1)
C-----INSURE SPECTRUM SPANS ENERGY RANGE OF GROUPS.
      CALL XYPAGE(1,1,ELOWD,SHIGH)
      CALL XYPAGE(MSPECT,1,EHIGHD,SHIGH)
      IF(ELOWD.LE.EGROUP(1).AND.EHIGHD.GE.EGROUP(NGRP1)) GO TO 490
      CALL OUT9G(ELOWD        ,FIELDX(1,1))
      CALL OUT9G(EHIGHD       ,FIELDX(1,2))
      CALL OUT9G(EGROUP(1)    ,FIELDX(1,3))
      CALL OUT9G(EGROUP(NGRP1),FIELDX(1,4))
      WRITE(OUTP,710) ((FIELDX(M,I),M=1,11),I=1,4)
      GO TO 530
C-----OPEN EVALUATED DATA FILE AND IF REQUIRED CREATE MULTI-GROUP
C-----DATA FILE (BOTH USE THE ENDF/B FORMAT).
  490 NSECT=4
C-----INITIALIZE MAXIMUM ERROR VECTOR.
      IF(NBMAX.LT.1) GO TO 510
      DO 500 NB=1,5
      DO 500 I=1,25
  500 ERLIB(I,NB)=0.0
C-----ALL INPUT DATA IS O.K.
  510 RETURN
C
C     ERROR IN INPUT DATA. PRINT MESSAGE AND TERMINATE.
C
  520 CALL ENDERROR
  530 WRITE(OUTP,730)
      GO TO 520
  540 FORMAT(' Multi-Group and Multi-Band Parameters from',
     1 ' ENDF/B Data (GROUPIE 2017-1)'/1X,78('-'))
  550 FORMAT(' Retrieval Criteria------------',7X,A4/
     1       ' Number of Energy Groups-------',I11,1X,A16/
     2       ' Maximum Number of Bands/Group-',I11,1X,A28/
     3       ' Points in Weighting Spectrum--',I11,1X,A36/
     4       ' Band Selection Criteria-------',11A1,
     5 ' (',F9.4,' per-cent)')
  560 FORMAT(1X,78('-')/
     1 ' Self-Shielded Listing---------',7X,A4,1X,A24/
     2 ' Multi-Band Listing------------',7X,A4/
     3 ' Multi-Band Library File-------',7X,A4/
     4 ' Unshielded Averages in ENDF/B-',7X,A4,1X,A16/
     5 ' Unshielded Averages Listing---',7X,A4,1X,A24)
  570 FORMAT(
     1 ' Maxwellian Temperature--------',1PE11.4,' eV')
  580 FORMAT(
     6 1X,78('-')/' Multi-Band Library Identification'/1X,78('-')/
     7 1X,A72/1X,78('-'))
  590 FORMAT(' MAT/MF/MT Ranges'/1X,78('-')/
     1 '    Minimum       Maximum'/
     2 '    MAT MF  MT    MAT MF  MT'/1X,78('-'))
  600 FORMAT(' ZA/MF/MT Ranges'/1X,78('-')/
     1 '    Minimum       Maximum'/
     2 '     ZA MF  MT     ZA MF  MT'/1X,78('-'))
  610 FORMAT(I7,I3,I4,I7,I3,I4)
  620 FORMAT(I7,I3,I4,I7,I3,I4,' (Default Option)')
  630 FORMAT(//' ERROR - Over 100 ranges----Execution Terminated')
  640 FORMAT(' Sigma-0 Definition------------',
     1 ' Multiple of Unshielded Total in Each Group')
  650 FORMAT(' Sigma-0 Definition------------',
     1 ' Same barns Values in Each Group')
  660 FORMAT(' High E Cross Section Extension',
     1 ' = Zero (Standard Convention)')
  670 FORMAT(' High E Cross Section Extension',
     1 ' = Constant (WARNING - not Standard)')
  680 FORMAT(/' Group Boundaries are NOT in Ascending Energy Order'/
     1 ' Group',I5,1X,11A1,' to ',11A1,' eV')
  690 FORMAT(/' Number of Groups=',I5,' (MUST be -19 to ',I5,')')
  700 FORMAT(/' Number of Groups=',I5,' (MUST be -19 to ',I5,')')
  710 FORMAT(/' Energy Spectrum from',11A1,' to ',11A1,' eV'/
     1 ' Does NOT Span the Group Range from',11A1,' to ',
     1 11A1,' eV')
  720 FORMAT(/' No Output Options Selected---Execution Terminated')
  730 FORMAT(' Execution Terminated')
  740 FORMAT(1X,78('-')/' Group Energy Boundaries'/1X,78('-')/
     1 6(3X,'Energy-eV')/1X,78('-'))
  750 FORMAT(6(1X,11A1))
  760 FORMAT(' Sigma0 values read as input'/1X,78('-')/
     1 1X,1P6E11.4/1X,1P6E11.4/1X,1P6E11.4/1X,1P3E11.4/1X,78('-'))
  770 FORMAT(' ERROR - Sigma0 NOT in DESCENDING value order',I5,
     1 1P2E11.4)
  780 FORMAT(1X,78('-')/' Execution terminated due to error in input'//)
      END
CAK   SUBROUTINE NXTMAT
      SUBROUTINE NXTMATg
C=======================================================================
C
C     FIND NEXT REQUESTED MATERIAL BASED EITHER ON ZA OR MAT.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*4 FMTHOL
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
CAK   COMMON/WHATZA/ATWT,IZA,MATNOW,MFNOW
      COMMON/WHATZAg/ATWT,IZA,MATNOW,MFNOW
CAK   COMMON/HOLFMT/FMTHOL
      COMMON/HOLFMTg/FMTHOL
      COMMON/VERSES/TEMP3,IVERSE
      COMMON/RESCOM/RES1,RES2,URES1,URES2,LSSF,ICOMP,INPART
C-----READ NEXT CARD AND CHECK FOR END OF ENDF/B TAPE.
   10 CALL CONTIG
      IF(MFH.gt.0) go to 20
      IF(MATH.lt.0) go to 30
      go to 10
C-----DEFINE FIXED POINT ZA.
   20 IZA  = C1H
      ATWT = C2H
C-----DEFINED WHETHER OR NOT THIS SECTION IS REQUESTED.
      MMM = MYWANT(IZA,MATH,MFH,MTH)
      IF(MMM.lt.0) go to 30
      IF(MMM.gt.0) go to 40
C-----NOT REQUESTED. SKIP SECTION.
      CALL SKIPS
      GO TO 10
C-----END OF RUN. RETURN NEGATIVE MATH AS INDICATOR.
   30 MATH=-1
      MFH=0
      MTH=0
C-----THIS MATERIAL REQUESTED. IF A NEW MAT INITIALIZE PARAMETERS.
   40 IF(MATH.EQ.MATNOW) GO TO 50
      NOSEQ=1
      IVERSE=6
      FMTHOL=' '
      MATNOW=MATH
C-----Initialize Resonance Region Boundaries.
      INPART= 1     ! Initialize to neutron (in case not ENDF6 format)
      RES1  = 0.0d0
      RES2  = 0.0d0
      URES1 = 0.0d0
      URES2 = 0.0d0
   50 MFNOW=MFH
      RETURN
      END
      INTEGER*4 FUNCTION MYWANT(IZA,MAT,MF,MT)
C=======================================================================
C
C     DEFINE WHETHER OR NOT THIS SECTION IS REQUESTED.
C
C=======================================================================
      INCLUDE 'implicit.h'
CAK   COMMON/MATZA/MODGET,NMATZA,MATMIN(101),MFMIN(101),MTMIN(101),
CAK  1 MATMAX(101),MFMAX(101),MTMAX(101)
      COMMON/MATZAg/MODGET,NMATZA,MATMIN(101),MFMIN(101),MTMIN(101),
     1 MATMAX(101),MFMAX(101),MTMAX(101)
      DIMENSION IZAMIN(101),IZAMAX(101)
      EQUIVALENCE (MATMIN(1),IZAMIN(1)),(MATMAX(1),IZAMAX(1))
      LOW=0
      IF(MODGET.NE.0) LOW=1
      DO 30 IMATZA=1,NMATZA
      IF(MODGET.NE.0) GO TO 10
C
C     MAT RANGES.
C
      IF(MAT.GT.MATMAX(IMATZA)) GO TO 30
      LOW=1
      IF(MAT.LT.MATMIN(IMATZA)) GO TO 30
      GO TO 20
C
C     ZA RANGES.
C
   10 IF(IZA.LT.IZAMIN(IMATZA).OR.IZA.GT.IZAMAX(IMATZA)) GO TO 30
C
C     MF/MT RANGES.
C
C-----COMPARE MF AND MT SELECTION CRITERIA.
   20 IF(MF.LT.MFMIN(IMATZA).OR.MF.GT.MFMAX(IMATZA)) GO TO 30
      IF(MT.LT.MTMIN(IMATZA).OR.MT.GT.MTMAX(IMATZA)) GO TO 30
C
C     THIS SECTION IS REQUESTED.
C
      GO TO 50
   30 CONTINUE
C
C     THIS SECTION HAS NOT BEEN REQUESTED. IF BEYOND RANGE OF ALL
C     REQUESTS RUN IF COMPLETED. IF NOT SKIP TO NEXT SECTION.
C
      IF(LOW.LE.0) GO TO 40
C-----NOT REQUESTED.
      MYWANT=0
      RETURN
C-----BEYOND THE RANGE OF ALL REQUESTS.
   40 MYWANT=-1
      RETURN
C-----REQUESTED.
   50 MYWANT=IMATZA
      RETURN
      END
      SUBROUTINE CONTIG
C=======================================================================
C
C     READ ONE ENDF/B CONTROL CARD.
C
C=======================================================================
      INCLUDE 'implicit.h'
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
C-----READ CARD.
   10 CALL CONTI
      IF(MTH .gt.0) go to 20     ! MTH > 0 = Start of SECTION
      IF(MFH .le.0) go to 30     ! MFH or MATH = 0 = End
      IF(MATH.eq.0) go to 30
      IF(MATH.lt.0) go to 50     ! MATH < 0 = TEND
      go to 10                   ! Skip SEND
C-----SKIP SECTION IF NOT REQUESTED.
   20 IZATRY=C1H
      MMM = MYWANT(IZATRY,MATH,MFH,MTH)
      IF(MMM.lt.0) go to 50
      IF(MMM.eq.0) go to 40
C-----REQUESTED OR END OF MF OR MAT.
   30 RETURN
C-----NOT REQUESTED. SKIP TO NEXT SECTION.
   40 CALL SKIPS
      GO TO 10
C-----END OF REQUESTED DATA.
   50 MATH=-1
      MFH=0
      MTH=0
      RETURN
      END
      SUBROUTINE CONTOG
C=======================================================================
C
C     WRITE ONE ENDF/B CONTROL RECORD.
C
C     PRECEED BY FEND, MEND AS NEEDED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      CHARACTER*4 CARDX(17)
      DIMENSION MFIELD(3)
      IF(OTAPE.LE.0) RETURN
C-----ADD FEND AND MEND LINES AS NEEDED.
      IF(LSTMAT.LE.0) GO TO 10
C-----ADD FEND LINE IF NEW MAT OR MF.
      IF(MATH.EQ.LSTMAT.AND.MFH.EQ.LSTMF) GO TO 10
      CALL OUTF(LSTMAT)
C-----ADD MEND LINE IF NEW MAT.
      IF(MATH.EQ.LSTMAT) GO TO 10
      CALL OUTM
      NOSEQ=1
   10 CALL CONTO
C-----SAVE LAST MAT/MF OUTPUT.
      LSTMAT = MATH
      LSTMF  = MFH
      RETURN
      ENTRY COPYFG
C=======================================================================
C
C     SAME AS COPYF, EXCEPT READ, BUT DO NOT OUTPUT FEND.
C     PLUS SAVE LAST MAT, MF OUTPUT.
C
C=======================================================================
   20 READ(ITAPE,30) CARDX,MFIELD
      IF(MFIELD(2).le.0) GO TO 40       ! If MF=0 return without write
      IF(OTAPE.LE.0) GO TO 20
      IF(MFIELD(3).le.0) then
      CALL OUTS(MFIELD(1),MFIELD(2))
      GO TO 20
      ENDIF
      WRITE(OTAPE,30) CARDX,MFIELD,NOSEQ
      NOSEQ = NXTSEQ(NOSEQ)
C-----SAVE LAST MAT/MF OUTPUT.
      LSTMAT = MFIELD(1)
      LSTMF  = MFIELD(2)
      GO TO 20
   30 FORMAT(16A4,A2,I4,I2,I3,I5)
   40 RETURN
      ENTRY COPYSG
C=======================================================================
C
C     SAME AS COPYS PLUS SAVE LAST MAT, MF OUTPUT
C
C=======================================================================
   50 READ(ITAPE,30) CARDX,MFIELD
      IF(OTAPE.LE.0) GO TO 60
      IF(MFIELD(3).gt.0) then
      WRITE(OTAPE,30) CARDX,MFIELD,NOSEQ
      NOSEQ = NXTSEQ(NOSEQ)
      else
      CALL OUTS(MFIELD(1),MFIELD(2))
      ENDIF
C-----SAVE LAST MAT/MF OUTPUT.
   60 LSTMAT = MFIELD(1)
      LSTMF  = MFIELD(2)
      IF(MFIELD(3).gt.0) go to 50
      RETURN
      ENTRY COPY1G
C=======================================================================
C
C     SAME AS COPY1 PLUS SAVE LAST MAT, MF OUTPUT
C
C=======================================================================
      READ(ITAPE,30) CARDX,MFIELD
      IF(OTAPE.LE.0) RETURN
      IF(MFIELD(3).gt.0) then
      WRITE(OTAPE,30) CARDX,MFIELD,NOSEQ
      NOSEQ = NXTSEQ(NOSEQ)
      else
      CALL OUTS(MFIELD(1),MFIELD(2))
      ENDIF
C-----SAVE LAST MAT/MF OUTPUT.
      LSTMAT = MFIELD(1)
      LSTMF  = MFIELD(2)
      RETURN
      END
      SUBROUTINE IBLOCK(ISCR,X,Y,IX)
C=======================================================================
C
C     READ A BLOCK OF WORDS FROM SCRATCH FILE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION X(IX),Y(IX)
      READ(ISCR) X,Y
      RETURN
      END
      SUBROUTINE OBLOCK(ISCR,X,Y,IX)
C=======================================================================
C
C     WRITE A BLOCK OF WORDS TO SCRATCH FILE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION X(IX),Y(IX)
      WRITE(ISCR) X,Y
      RETURN
      END
      SUBROUTINE BLOCK
C=======================================================================
C
C     INITIALIZE LABELLED COMMON.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*8 SIGHL1,SIGHL2,SIGHX1,SIGHX2
      CHARACTER*4 TUNITS,TUNIS,TMPHOL,REAC2,REAC3
      COMMON/ELPASD/TUNITS(2,4),TMPHOL(3)
      COMMON/SIGGIE/SIGHL1(12),SIGHL2(12)
      INCLUDE 'groupie.h'
      DIMENSION SIGMAZ(25),SHIELX(25),REAC2(2,6),REAC3(2,6),TUNIS(2,4),
     1 SIGHX1(12),SIGHX2(12)
C-----DEFINE SIGMA0 TABLE AS MULTIPLES OF UNSHIELDED GROUP
C-----AVERAGED CROSS SECTION IN EACH GROUP.
      DATA SIGMAZ/
     1   0.0D+00,1024.0D+00, 512.0D+00, 256.0D+00, 128.0D+00,
     3  64.0D+00,  32.0D+00,  16.0D+00,   8.0D+00,   4.0D+00,
     3   2.0D+00,   1.0D+00, 13*0.0D+00/
C-----DEFINE SIGMAO TABLE AS ABSOLUTE NUMBERS.
      DATA SHIELX/
     1   0.0D+00,40000.0D+00,20000.0D+00,10000.0D+00, 7000.0D+00,
     2            4000.0D+00, 2000.0D+00, 1000.0D+00,  700.0D+00,
     3             400.0D+00,  200.0D+00,  100.0D+00,   70.0D+00,
     4              40.0D+00,   20.0D+00,   10.0D+00,    7.0D+00,
     5               4.0D+00,    2.0D+00,    1.0D+00,    0.7D+00,
     6               0.4D+00,    0.0D+00,    0.0D+00,    0.0D+00/
C-----DEFINE TITLES FOR ABOVE VALUES OF SIGMAO FOR LISTING.
      DATA SIGHX1/'Averages','     256',
     1            '      64','      16',
     1            '       4','       1',
     2            '     1/4','    1/16',
     2            '    1/64','   1/256',
     3            '  1/Tot ','1/Tot**2'/
      DATA SIGHX2/'Averages','   10000',
     1            '    4000','    1000',
     1            '     400','     100',
     2            '      40','      10',
     2            '       4','       1',
     3            '  1/Tot ','1/Tot**2'/
C-----DEFINE HOLLERITH EQUIVALENT OF REACTIONS.
      DATA REAC2/
     1 '    ','    ',
     2 '   T','otal',
     3 ' Ela','stic',
     4 ' Cap','ture',
     5 ' Fis','sion',
     6 ' Oth','er  '/
      DATA REAC3/
     1 '    ','    ',
     2 'Tota','l   ',
     3 'Elas','t   ',
     4 'Capt','.   ',
     5 'Fiss','.   ',
     6 'Othe','r   '/
C-----DEFINE TEMPERATURE OUTOUT UNITS.
      DATA TUNIS/'Kelv','in  ','eV  ',' ','keV ',' ','MeV ',' '/
      DO 10 I=1,25
      SIGMAB(I)=SIGMAZ(I)
   10 SHIELD(I)=SHIELX(I)
      DO 20 I=1,2
      DO 20 J=1,6
   20 REACT2(I,J)=REAC2(I,J)
      DO 30 I=1,2
      DO 30 J=1,6
   30 REACT3(I,J)=REAC3(I,J)
      DO 40 I=1,2
      DO 40 J=1,4
   40 TUNITS(I,J)=TUNIS(I,J)
      DO 50 I=1,12
      SIGHL1(I) = SIGHX1(I)
   50 SIGHL2(I) = SIGHX2(I)
      RETURN
      END
      SUBROUTINE INTOUT(INTX,FIELD,LENGTH)
C=======================================================================
C
C     CONVERT INTEGER TO CHARACTERS FOR OUTPUT
C
C     WARNING - ONLY CONSIDERS POSITIVE INTEGERS
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 DIGITS,FIELD
      DIMENSION DIGITS(0:9),FIELD(LENGTH)
      DATA DIGITS/'0','1','2','3','4','5','6','7','8','9'/
C-----INITIALIZE TO BLANK
      DO 10 I=1,LENGTH
   10 FIELD(I)=' '
C-----FILL IN LAST DIGIT TO FIRST
      II=INTX
      DO 20 I=LENGTH,1,-1
      IF(II.LE.0) GO TO 30
      KK=II/10
      LL=II-10*KK
      FIELD(I)=DIGITS(LL)
   20 II=KK
   30 RETURN
      END
CAK   SUBROUTINE FILIO1
      SUBROUTINE FILIO1s
C=======================================================================
C
C     DEFINE ALL I/O UNITS AND OPTIONALLY DEFINE FILE NAMES.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE,OTAPE2
      CHARACTER*72 NAMEIN,NAMEOUT
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
CAK   COMMON/UNITS/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/UNITSg/OTAPE2,LIST1,LIST2,LIST3,IPLOT
      COMMON/IOSTATUS/ISTAT1,ISTAT2
      COMMON/NAMEX/NAMEIN,NAMEOUT
      COMMON/INPAGE/IXYLOW(5),IXYHI(5),ISCR(5)
C-----DEFINE ALL I/O UNITS. THEY ARE,
C-----INP     = INPUT FILE
C-----OUTP    = OUTPUT LISTING
C-----ITAPE   = ENDF/B FORMAT DATA TO BE AVERAGED
C-----OTAPE   = ENDF/B FORMAT MULTIGROUP OUTPUT (OPTIONAL)
C-----OTAPE2  = MULTIBAND PARAMETERS FILE (OPTIONAL)
C-----LIST1   = SELF-SHIELDED CROSS SECTION LISTING (OPTIONAL)
C-----LIST2   = MULTIBAND PARAMETERS LISTING (OPTIONAL)
C-----LIST3   = UNSHIELDED CROSS SECTION LISTING (OPTIONAL)
C-----ISCR(1) = SCRATCH FILE FOR WEIGHTING SPECTRUM
C-----ISCR(2) = SCRATCH FILE FOR TOTAL CROSS SECTION
C-----ISCR(3) = SCRATCH FILE FOR ELASTIC CROSS SECTION
C-----ISCR(4) = SCRATCH FILE FOR CAPTURE CROSS SECTION
C-----ISCR(5) = SCRATCH FILE FOR FISSION CROSS SECTION
C-----IPLOT   = PLOTTAB OUTPUT
      INP=2
      OUTP=3
      ITAPE=10
      OTAPE=11
      OTAPE2=31
      LIST1=32
      LIST2=33
      LIST3=34
      ISCR(1)=8
      ISCR(2)=9
      ISCR(3)=12
      ISCR(4)=13
      ISCR(5)=14
      IPLOT  =16
C-----DEFINE ALL FILE NAMES.
      OPEN(OUTP,FILE='GROUPIE.LST',STATUS='UNKNOWN')
      CALL SCRATCH1(ISCR(1),'GROUPIE.001 ')
      CALL SCRATCH1(ISCR(2),'GROUPIE.002 ')
      CALL SCRATCH1(ISCR(3),'GROUPIE.003 ')
      CALL SCRATCH1(ISCR(4),'GROUPIE.004 ')
      CALL SCRATCH1(ISCR(5),'GROUPIE.005 ')
C-----2010/3/12 - INITIALIZE PLOTTING FILE
      OPEN(IPLOT,FILE='PLOTTAB.CUR')
      OPEN(INP,FILE='GROUPIE.INP',STATUS='OLD',ERR=10)
      ISTAT1 = 0
      RETURN
   10 ISTAT1 = 1
      RETURN
CAK   ENTRY FILIO2
      ENTRY FILIO2g
C
C     DEFINE FILES (BASED ON REQUESTED OUTPUT).
C
      IF(OTAPE.GT.0) OPEN(OTAPE,FILE=NAMEOUT,STATUS='UNKNOWN')
      IF(OTAPE2.GT.0) OPEN(OTAPE2,FILE='MULTBAND.TAB',STATUS='UNKNOWN')
      IF(LIST1.GT.0) OPEN(LIST1,FILE='SHIELD.LST',STATUS='UNKNOWN')
      IF(LIST2.GT.0) OPEN(LIST2,FILE='MULTBAND.LST',STATUS='UNKNOWN')
      IF(LIST3.GT.0) OPEN(LIST3,FILE='UNSHIELD.LST',STATUS='UNKNOWN')
      IF(ITAPE.GT.0) OPEN(ITAPE,FILE=NAMEIN,STATUS='OLD',ERR=20)
      ISTAT2 = 0
      RETURN
   20 ISTAT2 = 1
      RETURN
      END
      SUBROUTINE OUT8(Z,ZCHAR)
C=======================================================================
C
C     FORMAT SIGMA0 FOR OUTPUT IN 8 COLUMNS
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 ZCHAR   ,DIGITS
      DIMENSION   ZCHAR(8),DIGITS(0:9)
      DATA DIGITS/'0','1','2','3','4','5','6','7','8','9'/
      DO I=1,8
      ZCHAR(I) = ' '
      ENDDO
C-----INTEGER IF > 1
      IF(Z.LT.1.0) GO TO 20
C
C     DIGITS LAST TO FIRST
C
      II = Z
      DO 10 I=8,2,-1
      JJ = II/10
      KK = II - 10*JJ
      ZCHAR(I) = DIGITS(KK)
      IF(JJ.LE.0) RETURN
   10 II = JJ
      RETURN
C
C     . FOLLOWED BY 6 DIGITS
C
   20 ZCHAR(2) = '.'
      DO I=3,8
      ZCHAR(I) = '0'
      ENDDO
      II = 1.0D+6*Z
      DO 30 I=8,3,-1
      JJ = II/10
      KK = II - 10*JJ
      ZCHAR(I) = DIGITS(KK)
      IF(JJ.LE.0) RETURN
   30 II = JJ
      RETURN
      END
C=======================================================================
C
C     ENDFIO.F includes routines to read and write ALL of the types
C     of ENDF records, e.g., CONT, TAB1, TAB2, LIST, etc.
C
C     The objective is to remove ALL direct reads and writes from
C     the ENDF Pre-Processing Codes (PREPRO) codes.
C
C     Note, that I (Red Cullen) have used these routines for so long
C     that some of the comments and variables still refer to CARDS -
C     as in way back when everything was based on 80 column computer
C     CARDS. I have purposely left this terminology in place, as a
C     reminder of how far we have come.
C
C     Version 2017-1 (May 2017)
C     =========================
C     *Treat End-of-File as TEND file during CONTI read.
C     *Added INNEXT to maintain INTEGER*8 logic for 9 to 10 digit outout
C     *Added POINTIV and LISTIV - added V = Variable input unit - to
C      read tables of PREPRO Input Parameters (rather than ENDF data).
C
C=======================================================================
C
C     OWNED, MAINTAINED AND DISTRIBUTED BY
C     ------------------------------------
C     THE NUCLEAR DATA SECTION
C     INTERNATIONAL ATOMIC ENERGY AGENCY
C     P.O. BOX 100
C     A-1400, VIENNA, AUSTRIA
C     EUROPE
C
C     ORIGINALLY WRITTEN BY
C     ------------------------------------
C     Dermott E. Cullen
C
C     PRESENT CONTACT INFORMATION
C     ---------------------------
C     Dermott E. Cullen
C     1466 Hudson Way
C     Livermore, CA 94550
C     U.S.A.
C     Telephone  925-443-1911
C     E. Mail    RedCullen1@Comcast.net
C     Website    RedCullen1.net/HOMEPAGE.NEW
C
C=======================================================================
      SUBROUTINE CARDIO(C1,C2,L1,L2,N1,N2)
C=======================================================================
C
C     READ AND WRITE ENDF/B CARD.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CALL CARDI(C1,C2,L1,L2,N1,N2)
      CALL CARDO(C1,C2,L1,L2,N1,N2)
      RETURN
      END
      SUBROUTINE CONTI
C=======================================================================
C
C     READ ONE ENDF/B CONTROL LINE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 FIELD2
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      DIMENSION FIELD2(11,2)
C-----READ FLOATING POINT FIELDS AS CHARACTERS
      READ(ITAPE,20,END=10) FIELD2,L1H,L2H,N1H,N2H,MATH,MFH,MTH
C-----TRANSLATE FROM CHARACTERS TO FLOATING POINT.
      CALL IN9(C1H,FIELD2(1,1))
      CALL IN9(C2H,FIELD2(1,2))
C-----ELIMINATE -0
      IF(IABS(L1H).LE.0) L1H=0
      IF(IABS(L2H).LE.0) L2H=0
      IF(IABS(N1H).LE.0) N1H=0
      IF(IABS(N2H).LE.0) N2H=0
      RETURN
c
c     ERROR - treat EOF as TEND line
c
   10 MATH = -1
      MFH  =  0
      MTH  =  0
      RETURN
   20 FORMAT(22A1,4I11,I4,I2,I3)
      END
      SUBROUTINE CONTO
C=======================================================================
C
C     WRITE ONE ENDF/B CONTROL RECORD.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD2
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      DIMENSION FIELD2(11,2)
C-----LAST MAT NUMBER - USED TO RESET SEQUENCE NUMBER
      DATA LASTMAT/-100000/
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----IF SEND, FEND OR MEND OUTPUT IN STANDARD FORM.
      IF(MTH.GT.0) GO TO 10
      CALL OUTS(MATH,MFH)
      RETURN
C-----CONVERT FLOATING POINT NUMBERS TO STANDARD OUTPUT FORM.
   10 CALL OUT9(C1H,FIELD2(1,1))
      CALL OUT9(C2H,FIELD2(1,2))
C-----ELIMINATE -0
      IF(IABS(L1H).LE.0) L1H=0
      IF(IABS(L2H).LE.0) L2H=0
      IF(IABS(N1H).LE.0) N1H=0
      IF(IABS(N2H).LE.0) N2H=0
C-----IF NEW MAT RESET SEQUENCE NUMBER
      IF(MATH.eq.LASTMAT) go to 20
      LASTMAT=MATH
      NOSEQ=1
C-----OUTPUT LINE IMAGE.
   20 IF(NOSEQ.LE.0) NOSEQ=1
      WRITE(OTAPE,30) FIELD2,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   30 FORMAT(22A1,4I11,I4,I2,I3,I5)
      END
      SUBROUTINE CARDI(C1,C2,L1,L2,N1,N2)
C=======================================================================
C
C     READ ENDF/B LINE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 FIELD2
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/LEADER/C1X,C2X,L1X,L2X,N1X,N2X,MAT,MF,MT
      DIMENSION FIELD2(11,2)
      READ(ITAPE,10) FIELD2,L1,L2,N1,N2,MAT,MF,MT
      CALL IN9(C1,FIELD2(1,1))
      CALL IN9(C2,FIELD2(1,2))
C-----ELIMINATE -0
      IF(IABS(L1).LE.0) L1=0
      IF(IABS(L2).LE.0) L2=0
      IF(IABS(N1).LE.0) N1=0
      IF(IABS(N2).LE.0) N2=0
      RETURN
   10 FORMAT(22A1,4I11,I4,I2,I3)
      END
      SUBROUTINE CARDO(C1,C2,L1,L2,N1,N2)
C=======================================================================
C
C     WRITE ENDF/B LINE.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD2
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      DIMENSION FIELD2(11,2)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
      IF(MTH.GT.0) GO TO 10
      CALL OUTS(MATH,MFH)
      RETURN
C-----CONVERT FLOATING POINT NUMBERS TO STANDARD OUTPUT FORM.
   10 CALL OUT9(C1,FIELD2(1,1))
      CALL OUT9(C2,FIELD2(1,2))
C-----ELIMINATE -0
      IF(IABS(L1).LE.0) L1=0
      IF(IABS(L2).LE.0) L2=0
      IF(IABS(N1).LE.0) N1=0
      IF(IABS(N2).LE.0) N2=0
C-----OUTPUT LINE.
      IF(NOSEQ.LE.0) NOSEQ=1
      WRITE(OTAPE,20) FIELD2,L1,L2,N1,N2,MATH,MFH,MTH,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   20 FORMAT(22A1,4I11,I4,I2,I3,I5)
      END
      SUBROUTINE TERPI(NBT,INT,N1)
C=======================================================================
C
C     READ TAB1 INTERPOLATION LAW.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DIMENSION NBT(N1),INT(N1)
      READ(ITAPE,10) (NBT(I),INT(I),I=1,N1)
      RETURN
   10 FORMAT(6I11)
      END
      SUBROUTINE TERPO(NBT,INT,N1)
C=======================================================================
C
C     WRITE TAB1 INTERPOLATION LAW.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      DIMENSION NBT(N1),INT(N1)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----LOOP OVER RANGES - UP TO 3 PER LINE
      DO 40 I1=1,N1,3
      I2=I1+2
      IF(I2.GT.N1) I2=N1
C-----OUTPUT LINE
      IOUT=(I2-I1)+1
      GO TO (10,20,30),IOUT
   10 WRITE(OTAPE,50) NBT(I1),INT(I1),MATH,MFH,MTH,NOSEQ
      GO TO 40
   20 WRITE(OTAPE,60) (NBT(II),INT(II),II=I1,I2),MATH,MFH,MTH,NOSEQ
      GO TO 40
   30 WRITE(OTAPE,70) (NBT(II),INT(II),II=I1,I2),MATH,MFH,MTH,NOSEQ
   40 NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   50 FORMAT(2I11,44X,I4,I2,I3,I5)
   60 FORMAT(4I11,22X,I4,I2,I3,I5)
   70 FORMAT(6I11    ,I4,I2,I3,I5)
      END
      SUBROUTINE POINTI(X,Y,IXY)
C=======================================================================
C
C     READ A PAGE OF DATA POINTS AND INSURE THAT THE ENERGIES ARE IN
C     ASCENDING ORDER. IF ENERGIES ARE NOT, TERMINATE EXECUTION.
C
C     WARNING - BEFORE STARTING TO READ EACH TABLE OF POINTS,
C               ELAST MUST BE INITIALIZED = 0, TO ALLOW ENERGY
C               ORDER TEST.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/LASTE/ELAST
      DIMENSION X(IXY),Y(IXY),FIELD6(11,6)
      DATA OKDIFF/1.0D-09/
C-----SET UP LOOP OVER LINES.
      DO 20 I=1,IXY,3
      II=I+2
      IF(II.GT.IXY) II=IXY
      IN=2*(II-I)+2
C-----READ ENERGY AS HOLLERITH AND CROSS SECTION AS FLOATING POINT.
      READ(ITAPE,50) FIELD6
      J=I
C-----CONVERT ENERGY TO FLOATING POINT.
      DO 10 K=1,IN,2
      CALL IN9(X(J),FIELD6(1,K))
      CALL IN9(Y(J),FIELD6(1,K+1))
   10 J=J+1
   20 CONTINUE
C-----CHECK ENERGY ORDER.
      DO 40 I=1,IXY
      IF(X(I).GE.ELAST) GO TO 40
C-----ALLOW FOR SMALL DIFFERENCES (ABOUT THE SAME TO 9 DIGITS).
      IF(DABS(ELAST-X(I)).LE.OKDIFF*ELAST) GO TO 30
      CALL OUT9(ELAST,FIELD6(1,1))
      CALL OUT9(X(I) ,FIELD6(1,2))
      WRITE(OUTP,60) MATH,MFH,MTH,
     1 I-1,(FIELD6(M,1),M=1,11),
     2 I  ,(FIELD6(M,2),M=1,11)
C-----WHEN SMALL DIFFERENCES OCCUR INSURE THAT ENERGIES ARE NOT IN
C-----DESCENDING ORDER.
   30 X(I)=ELAST
   40 ELAST=X(I)
      RETURN
   50 FORMAT(66A1)
   60 FORMAT(2X,78('-')/I5,I3,I4/
     1 ' Energies Not in Ascending Energy Order'/
     2 '  Index      Energy'/
     3 I7,1X,11A1      /I7,1X,11A1      /
     4 19X,' Execution Terminated.'/2X,78('-'))
      END
      SUBROUTINE POINTIV(IUNIT,X,Y,IXY)
C=======================================================================
C
C     POINTI + V = VARIABLE INPUT UNIT (IUNIT).
C
C     WARNING - NO X Order Test (as in POINTI).
C
C     Used to READ PREPRO Input Parameters.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 FIELD6
      DIMENSION X(IXY),Y(IXY),FIELD6(11,6)
C-----SET UP LOOP OVER LINES.
      DO 20 I=1,IXY,3
      II=I+2
      IF(II.GT.IXY) II=IXY
      IN=2*(II-I)+2
C-----READ ENERGY AS HOLLERITH AND CROSS SECTION AS FLOATING POINT.
      READ(IUNIT,50) FIELD6
      J=I
C-----CONVERT ENERGY TO FLOATING POINT.
      DO 10 K=1,IN,2
      CALL IN9(X(J),FIELD6(1,K))
      CALL IN9(Y(J),FIELD6(1,K+1))
   10 J=J+1
   20 CONTINUE
      RETURN
   50 FORMAT(66A1)
      END
      SUBROUTINE POINTO(X,Y,IXY)
C=======================================================================
C
C     WRITE IXY DATA POINTS. FORMAT OF ENERGIES WILL VARY TO ALLOW
C     MAXIMUM PRECISION OF BETWEEN 6 AND 9 DIGITS ACCURACY.
C
C     CROSS SECTIONS WILL ALWAYS BE OUTPUT E11.4 FORMAT.
C
C     PHOTON DATA WILL ALWAYS BE IN E11.4 FORMAT.
C
C     ENERGIES WILL BE OUTPUT IN EITHER STANDARD E11.4 FORMAT OR A
C     VARIABLE F FORMAT (VARIABLE FROM F11.8 TO F11.0) TO GIVE THE
C     MAXIMUM NUMBER OF DIGITS OF ACCURACY. AS OUTPUT BY THIS ROUTINE
C     STANDARD FORM E11.4 FORMAT GIVES 6 DIGITS OF ACCURACY. THE
C     VARIABLE FORM F FORMAT WILL GIVE 6 TO 9 DIGITS ACCURACY. AS
C     LONG AS THE EXPONENT OF AN ENERGY IN E11.4 FORMAT IS BETWEEN -3
C     AND +8, MORE DIGITS WILL BE INCLUDED IF THE NUMBER IS OUTPUT IN
C     VARIABLE F FORMAT. IN PARTICULAR A FULL 9 DIGITS WILL BE OUTPUT
C     FOR ALL ENERGIES BETWEEN 1 EV AND 100 MEV. BETWEEN 1 MILLI-EV
C     AND 1 EV THE NUMBER OF DIGITS WILL VARY FROM 6 TO 8.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/FLAGS/MINUS3,IMPLUS
      DIMENSION X(IXY),Y(IXY),FIELD6(11,6)
      DATA ZEROD/0.0D+00/
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----NOTHING TO DO IF NO POINTS
      IF(IXY.LE.0) RETURN
C-----SET UP LOOP OVER LINES (UP TO 3 POINTS PER LINE).
      DO 50 I1=1,IXY,3
      I2=I1+2
      IF(I2.GT.IXY) I2=IXY
C
C     OUTPUT ONE LINE WITH ENERGY IN F OR E FORMAT AND CROSS SECTION
C     IN E FORMAT.
C
C-----CONVERT DATA TO NORMAL FORM.
      K=0
      DO 10 II=I1,I2
C-----Negative Imaginary Anomalous Scattering is O.K.
      IF(MFH.ne.27.or.MTH.ne.506) then
C-----COUNT IF CROSS SECTION IS NEGATIVE.
      IF(Y(II).LT.ZEROD) MINUS3=MINUS3+1
      ENDIF
C-----SET FLAG IF POSITIVE.
      IF(Y(II).GT.ZEROD) IMPLUS=1
      K=K+1
c-----2013/1/12 - changed ENERGY to OUT10 from OUT9.
      CALL OUT10(X(II),FIELD6(1,K))
      K=K+1
C-----CHANGED CROSS SECTION TO 9 DIGIT OUTPUT
   10 CALL OUT9(Y(II),FIELD6(1,K))
C-----OUTPUT ONE LINE.
      IOUT=(I2-I1)+1
      GO TO (20,30,40),IOUT
   20 WRITE(OTAPE,60) ((FIELD6(M,II),M=1,11),II=1,2),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 50
   30 WRITE(OTAPE,70) ((FIELD6(M,II),M=1,11),II=1,4),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 50
   40 WRITE(OTAPE,80) ((FIELD6(M,II),M=1,11),II=1,6),
     1 MATH,MFH,MTH,NOSEQ
   50 NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   60 FORMAT(22A1,44X,I4,I2,I3,I5)
   70 FORMAT(44A1,22X,I4,I2,I3,I5)
   80 FORMAT(66A1    ,I4,I2,I3,I5)
      END
      SUBROUTINE POINTO9(X,Y,IXY)
C=======================================================================
C
C     Identical to POINTO, but uses OUT9G rather than OUT10 =
C     Used by GROUPIE - 10 digits never needed for  multi-group.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/FLAGS/MINUS3,IMPLUS
      DIMENSION X(IXY),Y(IXY),FIELD6(11,6)
      DATA ZEROD/0.0D+00/
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----NOTHING TO DO IF NO POINTS
      IF(IXY.LE.0) RETURN
C-----SET UP LOOP OVER LINES (UP TO 3 POINTS PER LINE).
      DO 50 I1=1,IXY,3
      I2=I1+2
      IF(I2.GT.IXY) I2=IXY
C
C     OUTPUT ONE LINE WITH ENERGY IN F OR E FORMAT AND CROSS SECTION
C     IN E FORMAT.
C
C-----CONVERT DATA TO NORMAL FORM.
      K=0
      DO 10 II=I1,I2
C-----Negative Imaginary Anomalous Scattering is O.K.
      IF(MFH.ne.27.or.MTH.ne.506) then
C-----COUNT IF CROSS SECTION IS NEGATIVE.
      IF(Y(II).LT.ZEROD) MINUS3=MINUS3+1
      ENDIF
C-----SET FLAG IF POSITIVE.
      IF(Y(II).GT.ZEROD) IMPLUS=1
      K=K+1
c-----2013/1/12 - changed ENERGY to OUT10 from OUT9.
c-----2014/3/27 - added POINTO9 using OUT9 - this is the only
c-----            difference between POINTO and POINTO9.
c-----2015/7/30 - Changed OUT9 to OIUT9G.
      CALL OUT9G(X(II),FIELD6(1,K))
      K=K+1
C-----CHANGED CROSS SECTION TO 9 DIGIT OUTPUT
   10 CALL OUT9G(Y(II),FIELD6(1,K))
C-----OUTPUT ONE LINE.
      IOUT=(I2-I1)+1
      GO TO (20,30,40),IOUT
   20 WRITE(OTAPE,60) ((FIELD6(M,II),M=1,11),II=1,2),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 50
   30 WRITE(OTAPE,70) ((FIELD6(M,II),M=1,11),II=1,4),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 50
   40 WRITE(OTAPE,80) ((FIELD6(M,II),M=1,11),II=1,6),
     1 MATH,MFH,MTH,NOSEQ
   50 NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   60 FORMAT(22A1,44X,I4,I2,I3,I5)
   70 FORMAT(44A1,22X,I4,I2,I3,I5)
   80 FORMAT(66A1    ,I4,I2,I3,I5)
      END
      SUBROUTINE LISTIO(X,IX)
C=======================================================================
C
C     READ AND WRITE LIST DATA.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION X(IX)
      CALL LISTI(X,IX)
      CALL LISTO(X,IX)
      RETURN
      END
      SUBROUTINE LISTIO9(X,IX)
C=======================================================================
C
C     READ AND WRITE LIST DATA.
C
C     SAME AS LISTIO, WITHOUT OUT10
C
C=======================================================================
      INCLUDE 'implicit.h'
      DIMENSION X(IX)
      CALL LISTI(X,IX)
      CALL LISTO9(X,IX)
      RETURN
      END
      SUBROUTINE LISTSKIP(IX)
C=======================================================================
C
C     SKIP LIST DATA.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 DUMMY
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DATA DUMMY/' '/
      DO 10 I=1,IX,6
   10 READ(ITAPE,20) DUMMY
C-----USE DUMMY TO PREVENT COMPILER WARNING
      IF(DUMMY.NE.' ') I=1
      RETURN
   20 FORMAT(A1)
      END
      SUBROUTINE LISTI(X,IX)
C=======================================================================
C
C     READ LIST DATA.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DIMENSION X(IX),FIELD6(11,6)
C
C     READ LIST RECORD PARAMETERS
C
C-----SET UP LOOP OVER CARDS.
      DO 20 I1=1,IX,6
      I2=I1+5
      IF(I2.GT.IX) I2=IX
C-----READ AS CHARACTERS
      READ(ITAPE,30) FIELD6
C-----CONVERT FROM CHARACTERS TO FLOATING POINT
      K=0
      DO 10 L=I1,I2
      K=K+1
   10 CALL IN9(X(L),FIELD6(1,K))
   20 CONTINUE
      RETURN
   30 FORMAT(66A1)
      END
      SUBROUTINE LISTIV(IUNIT,X,IX)
C=======================================================================
C
C     READ LIST DATA LISTI + V = With Variable Input unit (IUNIT)
C
C     Used to Read PREPRO Input Parameters.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 FIELD6
      DIMENSION X(IX),FIELD6(11,6)
C
C     READ LIST RECORD PARAMETERS
C
C-----SET UP LOOP OVER CARDS.
      DO 20 I1=1,IX,6
      I2=I1+5
      IF(I2.GT.IX) I2=IX
C-----READ AS CHARACTERS
      READ(IUNIT,30) FIELD6
C-----CONVERT FROM CHARACTERS TO FLOATING POINT
      K=0
      DO 10 L=I1,I2
      K=K+1
   10 CALL IN9(X(L),FIELD6(1,K))
   20 CONTINUE
      RETURN
   30 FORMAT(66A1)
      END
      SUBROUTINE LISTO(X,IX)
C=======================================================================
C
C     WRITE LIST DATA
C
C     2013/1/12 - changed to OUT10 from OUT9.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      DIMENSION X(IX),FIELD6(11,6)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----NOTHING TO DO IF NO POINTS TO OUTPUT
      IF(IX.LE.0) RETURN
C-----SET UP LOOP OVER CARDS.
      DO 80 I1=1,IX,6
      I2=I1+5
      IF(I2.GT.IX) I2=IX
C-----CONVERT DATA TO NORMAL FORM.
      K=0
      DO 10 L=I1,I2
      K=K+1
c-----2013/1/12 - changed to OUT10 from OUT9.
   10 CALL OUT10(X(L),FIELD6(1,K))
C-----OUTPUT ONE LINE.
      GO TO (20,30,40,50,60,70),K
   20 WRITE(OTAPE,90)   (FIELD6(M,1),M=1,11),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   30 WRITE(OTAPE,100) ((FIELD6(M,II),M=1,11),II=1,2),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   40 WRITE(OTAPE,110) ((FIELD6(M,II),M=1,11),II=1,3),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   50 WRITE(OTAPE,120) ((FIELD6(M,II),M=1,11),II=1,4),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   60 WRITE(OTAPE,130) ((FIELD6(M,II),M=1,11),II=1,5),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   70 WRITE(OTAPE,140) ((FIELD6(M,II),M=1,11),II=1,6),
     1 MATH,MFH,MTH,NOSEQ
   80 NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   90 FORMAT(11A1,55X,I4,I2,I3,I5)
  100 FORMAT(22A1,44X,I4,I2,I3,I5)
  110 FORMAT(33A1,33X,I4,I2,I3,I5)
  120 FORMAT(44A1,22X,I4,I2,I3,I5)
  130 FORMAT(55A1,11X,I4,I2,I3,I5)
  140 FORMAT(66A1    ,I4,I2,I3,I5)
      END
      SUBROUTINE LISTO9(X,IX)
C=======================================================================
C
C     WRITE LIST DATA
C
C     SAME AS LISTO, BUT USE OUT9G NOT OUT10.
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*1 FIELD6
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      DIMENSION X(IX),FIELD6(11,6)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----NOTHING TO DO IF NO POINTS TO OUTPUT
      IF(IX.LE.0) RETURN
C-----SET UP LOOP OVER CARDS.
      DO 80 I1=1,IX,6
      I2=I1+5
      IF(I2.GT.IX) I2=IX
C-----CONVERT DATA TO NORMAL FORM.
      K=0
      DO 10 L=I1,I2
      K=K+1
c-----USE OUT9G - NOT OUT10.
   10 CALL OUT9G(X(L),FIELD6(1,K))
C-----OUTPUT ONE LINE.
      GO TO (20,30,40,50,60,70),K
   20 WRITE(OTAPE,90)   (FIELD6(M,1),M=1,11),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   30 WRITE(OTAPE,100) ((FIELD6(M,II),M=1,11),II=1,2),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   40 WRITE(OTAPE,110) ((FIELD6(M,II),M=1,11),II=1,3),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   50 WRITE(OTAPE,120) ((FIELD6(M,II),M=1,11),II=1,4),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   60 WRITE(OTAPE,130) ((FIELD6(M,II),M=1,11),II=1,5),
     1 MATH,MFH,MTH,NOSEQ
      GO TO 80
   70 WRITE(OTAPE,140) ((FIELD6(M,II),M=1,11),II=1,6),
     1 MATH,MFH,MTH,NOSEQ
   80 NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   90 FORMAT(11A1,55X,I4,I2,I3,I5)
  100 FORMAT(22A1,44X,I4,I2,I3,I5)
  110 FORMAT(33A1,33X,I4,I2,I3,I5)
  120 FORMAT(44A1,22X,I4,I2,I3,I5)
  130 FORMAT(55A1,11X,I4,I2,I3,I5)
  140 FORMAT(66A1    ,I4,I2,I3,I5)
      END
      SUBROUTINE LINEIN
C=======================================================================
C
C     READ A LINE, INCLUDING MAT, MF, MT
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*4 CARD
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/COPC/CARD(17)
      READ(ITAPE,10) CARD,MATH,MFH,MTH
      RETURN
   10 FORMAT(16A4,A2,I4,I2,I3)
      END
      SUBROUTINE LINEOUT
C=======================================================================
C
C     WRITE A LINE
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*4 CARD
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/COPC/CARD(17)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----USE STANDARD FORM FOR END LINES
      IF(MTH.GT.0) GO TO 10
      CALL OUTS(MATH,MFH)
      RETURN
   10 WRITE(OTAPE,20) CARD,MATH,MFH,MTH,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   20 FORMAT(16A4,A2,I4,I2,I3,I5)
      END
      SUBROUTINE COPYT
C=======================================================================
C
C     COPY TO TEND, MEND, FEND OR SEND RECORDS.
C     ENTRY POINTS ARE,
C     COPYT = COPY TO TEND RECORD
C     COPYM = COPY TO MEND RECORD
C     COPYF = COPY TO FEND RECORD
C     COPYS = COPY TO SEND RECORD
C     COPYL = COPY TAPE LABEL
C     COPY1 = COPY ONE LINE
C
C     COPY TO TEND RECORD
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      CHARACTER*4 CARD
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/COPC/CARD(17)
      COMMON/COPI/MFIELD(3)
   10 READ(ITAPE,150) CARD,MFIELD
      IF(OTAPE.LE.0) GO TO 30
      IF(MFIELD(3).GT.0) GO TO 20
      CALL OUTS(MFIELD(1),MFIELD(2))
      GO TO 30
   20 WRITE(OTAPE,150) CARD,MFIELD,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
C-----NEED MAT < 0
   30 IF(MFIELD(1).ge.0) go to 10
      RETURN
C=======================================================================
C
C     COPY TO MEND RECORD
C
C=======================================================================
      ENTRY COPYM
   40 READ(ITAPE,150) CARD,MFIELD
      IF(OTAPE.LE.0) GO TO 60
      IF(MFIELD(3).GT.0) GO TO 50
      CALL OUTS(MFIELD(1),MFIELD(2))
      GO TO 60
   50 WRITE(OTAPE,150) CARD,MFIELD,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
C-----NEED MAT <= 0
   60 IF(MFIELD(1).gt.0) go to 40
      RETURN
C=======================================================================
C
C     COPY TO FEND RECORD
C
C=======================================================================
      ENTRY COPYF
   70 READ(ITAPE,150) CARD,MFIELD
      IF(OTAPE.LE.0) GO TO 90
      IF(MFIELD(3).GT.0) GO TO 80
      CALL OUTS(MFIELD(1),MFIELD(2))
      GO TO 90
   80 WRITE(OTAPE,150) CARD,MFIELD,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
C-----NEED MF <= 0
   90 IF(MFIELD(2).gt.0) go to 70
      RETURN
C=======================================================================
C
C     COPY TO SEND RECORD
C
C=======================================================================
      ENTRY COPYS
  100 READ(ITAPE,150) CARD,MFIELD
      IF(OTAPE.LE.0) GO TO 120
      IF(MFIELD(3).GT.0) GO TO 110
      CALL OUTS(MFIELD(1),MFIELD(2))
      GO TO 120
  110 WRITE(OTAPE,150) CARD,MFIELD,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
C-----NEED MT <= 0
  120 IF(MFIELD(3).gt.0) go to 100
      RETURN
C=======================================================================
C
C     COPY TAPE LABEL
C
C=======================================================================
      ENTRY COPYL
      READ(ITAPE,150) CARD,MFIELD
      IF(OTAPE.LE.0) RETURN
C-----MF/MT/NOSEQ = 0 ON TEND LINE
      NOSEQ=0
C-----USE STANDARD TAPE NUMBER IF INPUT = 0
      IF(MFIELD(1).EQ.0) MFIELD(1)=7000      ! 2014/4/20 6000 to 7000
      MFIELD(2)=0
      MFIELD(3)=0
      WRITE(OTAPE,150) CARD,MFIELD,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
      RETURN
C=======================================================================
C
C     COPY ONE LINE
C
C=======================================================================
      ENTRY COPY1
      READ(ITAPE,150) CARD,MFIELD
      IF(OTAPE.LE.0) RETURN
      IF(MFIELD(3).GT.0) GO TO 130
      CALL OUTS(MFIELD(1),MFIELD(2))
      GO TO 140
  130 WRITE(OTAPE,150) CARD,MFIELD,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
  140 RETURN
  150 FORMAT(16A4,A2,I4,I2,I3,I5)
      END
      SUBROUTINE OUTT
C=======================================================================
C
C     OUTPUT TEND, MEND, FEND OR SEND RECORDS.
C     ENTRY POINTS ARE,
C     OUTT = OUTPUT TEND RECORD
C     OUTM = OUTPUT MEND RECORD
C     OUTF = OUTPUT FEND RECORD
C     OUTS = OUTPUT SEND RECORD
C
C     OUTPUT TEND RECORD
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
      WRITE(OTAPE,10)
      RETURN
C=======================================================================
C
C     OUTPUT MEND RECORD
C
C=======================================================================
      ENTRY OUTM
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
      WRITE(OTAPE,20)
C-----RESET NOSEQ AFTER MEND OUTPUT
      NOSEQ=1
      RETURN
C=======================================================================
C
C     OUTPUT FEND RECORD
C
C=======================================================================
      ENTRY OUTF(MATOUT)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
      WRITE(OTAPE,30) MATOUT
      NOSEQ=1
      RETURN
C=======================================================================
C
C     OUTPUT SEND RECORD
C
C=======================================================================
      ENTRY OUTS(MATOUT,MFOUT)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
C-----NOSEQ = 0 ON FEND/MEND/TEND LINE
      IF(MFOUT.LE.0) then
      WRITE(OTAPE,30) MATOUT
      else
      WRITE(OTAPE,40) MATOUT,MFOUT  ! NOSEQ = 99999
      endif
      NOSEQ=1
      RETURN
   10 FORMAT(66X,'  -1 0  0    0')
   20 FORMAT(66X,'   0 0  0    0')
   30 FORMAT(66X,I4, ' 0  0    0')
   40 FORMAT(66X,I4,I2,'  099999')
      END
      SUBROUTINE SKIPT
C=======================================================================
C
C     SKIP TO TEND, MEND, FEND OR SEND RECORDS.
C     ENTRY POINTS ARE,
C     SKIPT = SKIP TO TEND RECORD
C     SKIPM = SKIP TO MEND RECORD
C     SKIPF = SKIP TO FEND RECORD
C     SKIPS = SKIP TO SEND RECORD
C     SKIPL = SKIP TAPE LABEL
C     SKIP1 = SKIP ONE LINE
C
C     SKIP TO TEND RECORD
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      COMMON/COPI/MFIELD(3)
   10 READ(ITAPE,50) MFIELD
C-----NEED MAT < 0
      IF(MFIELD(1).ge.0) go to 10
      RETURN
C=======================================================================
C
C     SKIP TO MEND RECORD
C
C=======================================================================
      ENTRY SKIPM
   20 READ(ITAPE,50) MFIELD
C-----NEED MAT <= 0
      IF(MFIELD(1).gt.0) go to 20
      RETURN
C=======================================================================
C
C     SKIP TO FEND RECORD
C
C=======================================================================
      ENTRY SKIPF
   30 READ(ITAPE,50) MFIELD
C-----NEED MF <= 0
      IF(MFIELD(2).gt.0) go to 30
      RETURN
C=======================================================================
C
C     SKIP TO SEND RECORD
C
C=======================================================================
      ENTRY SKIPS
   40 READ(ITAPE,50) MFIELD
C-----NEED MT <= 0
      IF(MFIELD(3).gt.0) go to 40
      RETURN
C=======================================================================
C
C     SKIP TAPE LABEL
C
C=======================================================================
      ENTRY SKIPL
      READ(ITAPE,50) MFIELD
      RETURN
C=======================================================================
C
C     SKIP ONE LINE
C
C=======================================================================
      ENTRY SKIP1
      READ(ITAPE,50) MFIELD
      RETURN
   50 FORMAT(66X,I4,I2,I3,I5)
      END
      SUBROUTINE HOLLYI(LINE66)
C=======================================================================
C
C     READ A LINE OF 66 CHARACTERS
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 LINE66
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      READ(ITAPE,10) LINE66
      RETURN
   10 FORMAT(66A1,I4,I2,I3,I5)
      END
      SUBROUTINE HOLLYO(LINE66)
C=======================================================================
C
C     WRITE A LINE OF 66 CHARACTERS
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*1 LINE66
      INTEGER*4 OUTP,OTAPE
      COMMON/HEADER/C1H,C2H,L1H,L2H,N1H,N2H,MATH,MFH,MTH,NOSEQ
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DIMENSION LINE66(66)
C-----NO OUTPUT IF OUTPUT UNIT IS TURNED OFF
      IF(OTAPE.LE.0) RETURN
      WRITE(OTAPE,10) LINE66,MATH,MFH,MTH,NOSEQ
      NOSEQ=NXTSEQ(NOSEQ)
      RETURN
   10 FORMAT(66A1,I4,I2,I3,I5)
      END
      FUNCTION NXTSEQ(NNSEQ)
C=======================================================================
C
C     DEFINE NEXT SEQUENCE NUMBER FOR ENDF/B OUTPUT. ALLOW FOR
C     MORE THAN 100000 LINES PER EVALUATION BY RESETTING NUMBER
C     TO 1 EVERY TIME 100000 IS REACHED.
C
C=======================================================================
      INCLUDE 'implicit.h'
      NN=NNSEQ+1
      IF(NN.EQ.100000) NN=1
      NXTSEQ=NN
      RETURN
      END
      SUBROUTINE ZAHOL(ZA,ZABCD)
C=======================================================================
C
C     GIVEN ANY ZA (1000*Z+A) THIS ROUTINE WILL DEFINE A 10 CHARACTER
C     EQUIVALENT IN THE FORM,
C
C     10/07/04 - ADDED Z = 104 THROUGH 110.
C
C     CHARACTER POSITION (1 THROUGH 10)
C              1
C     1234567890
C
C     ZZZ-SS-AAA
C
C     ZZZ  - CHARACTER REPRESENTATION OF Z
C     SS   - CHEMICAL SYMBOL FOR ELEMENT
C     AAA  - CHARACTER REPRESENTATION FOR A OR NAT, IF A = 0
C
C     Z IS RIGHT ADJUSTED TO END IN CHARACTER 3
C     A IS LEFT ADJUSTED TO START IN CHARACTER 8
C
C     EXAMPLE, ZA = 6012 IS RETURNED AS,
C
C     CHARACTER POSITION (1 THROUGH 10)
C              1
C     1234567890
C
C       6-C -12
C
C=======================================================================
      INCLUDE 'implicit.h'
      INTEGER*4 ZA,Z,A
      CHARACTER*1 DUM1,DUM2,ZABCD,ZATAB,DIGITS,NEUTRON,FISSPRO,PHOTON
      DIMENSION ZATAB(2,110),DUM1(2,54),DUM2(2,56),ZABCD(10),
     1 DIGITS(10),NEUTRON(10),FISSPRO(10),PHOTON(10)
      EQUIVALENCE (ZATAB(1,1),DUM1(1,1)),(ZATAB(1,55),DUM2(1,1))
      DATA DIGITS/'0','1','2','3','4','5','6','7','8','9'/
      DATA DUM1/
     1 'H',' ','H','e','L','i','B','e','B',' ','C',' ',
     2 'N',' ','O',' ','F',' ','N','e','N','a','M','g',
     3 'A','l','S','i','P',' ','S',' ','C','l','A','r',
     4 'K',' ','C','a','S','c','T','i','V',' ','C','r',
     5 'M','n','F','e','C','o','N','i','C','u','Z','n',
     6 'G','a','G','e','A','s','S','e','B','r','K','r',
     7 'R','b','S','r','Y',' ','Z','r','N','b','M','o',
     8 'T','c','R','u','R','h','P','d','A','g','C','d',
     9 'I','n','S','n','S','b','T','e','I',' ','X','e'/
      DATA DUM2/
     1 'C','s','B','a','L','a','C','e','P','r','N','d',
     2 'P','m','S','m','E','u','G','d','T','b','D','y',
     3 'H','o','E','r','T','m','Y','b','L','u','H','f',
     4 'T','a','W',' ','R','e','O','s','I','r','P','t',
     5 'A','u','H','g','T','l','P','b','B','i','P','o',
     6 'A','t','R','n','F','r','R','a','A','c','T','h',
     7 'P','a','U',' ','N','p','P','u','A','m','C','m',
     8 'B','k','C','f','E','s','F','m','M','d','N','o',
     9 'L','r','R','f','D','b','S','g','B','h','H','a',
     A 'M','t','D','s'/
      DATA NEUTRON/' ','N','e','u','t','r','o','n',' ',' '/
      DATA PHOTON /' ','P','h','o','t','o','n',' ',' ',' '/
      DATA FISSPRO/'F','i','s','s','.','P','r','o','d','.'/
C
C     SPECIAL TREATMENT FOR ENDL NEUTRON, FISSION PRODUCTS
C     AND PHOTON
C
C-----NEUTRON?
      IF(ZA.NE.1) GO TO 20
      DO 10 I=1,10
   10 ZABCD(I)=NEUTRON(I)
      RETURN
C-----FISSION PRODUCT?
   20 IF(ZA.NE.99120.AND.ZA.NE.99125) GO TO 40
      DO 30 I=1,10
   30 ZABCD(I)=FISSPRO(I)
      RETURN
C-----PHOTON?
   40 IF(ZA.EQ.0) THEN
      DO 50 I=1,10
   50 ZABCD(I)=PHOTON(I)
      RETURN
      ENDIF
C
C     NORMAL TREATMENT
C
C-----BLANK OUT ZABCD TO START.
      DO 60 I=1,10
   60 ZABCD(I)=' '
C-----DEFINE Z AND A SEPARATELY.
      Z=ZA/1000
      A=ZA-1000*Z
C-----DEFINE SYMBOL FOR ELEMENT.
      ZABCD(4)='-'
      ZABCD(7)='-'
      IF(Z.GT.0.AND.Z.LE.110) GO TO 70
      ZABCD(5)='?'
      ZABCD(6)='?'
      IF(Z.LT.0.OR.Z.GT.999) GO TO 100
      GO TO 80
   70 ZABCD(5)=ZATAB(1,Z)
      ZABCD(6)=ZATAB(2,Z)
C-----DEFINE Z LAST DIGIT TO FIRST.
   80 II=3
      DO 90 I=1,3
      NEXTZ=Z/10
      KZ=Z-10*NEXTZ
      ZABCD(II)=DIGITS(KZ+1)
      Z=NEXTZ
      IF(Z.LE.0) GO TO 100
   90 II=II-1
  100 IF(A.GT.0) GO TO 110
C-----NATURAL ISOTOPIC MIXTURE.
      ZABCD(8) ='N'
      ZABCD(9) ='a'
      ZABCD(10)='t'
      GO TO 140
C-----DEFINE A FIRST DIGIT TO LAST.
  110 IDIV=100
      IMON=0
      II=7
      DO 130 I=1,3
      IA=A/IDIV
      IF(IA.EQ.0.AND.IMON.EQ.0) GO TO 120
      IMON=1
      II=II+1
      ZABCD(II)=DIGITS(IA+1)
  120 A=A-IDIV*IA
  130 IDIV=IDIV/10
  140 RETURN
      END
      REAL*8 FUNCTION TERPIT(X,X1,X2,Y1,Y2,INTERP)
C=======================================================================
C
C     INTERPOLATION ACCORDING TO ENDF/B LAWS 1 THROUGH 6.
C
C     INTERPOLATE BETWEEN (X1,Y1) AND (X2,Y2) TO DEFINE
C     Y AT X.
C
C     WARNING - THIS ROUTINE DOES NOT CHECK THE CONSISTENCY
C     BETWEEN DATA AND THE ENDF/B INTERPOLATION LAW - THEREFORE
C     TO AVOID ERRORS DURING EXECUTION THE USER MUST CHECK
C     CONSISTENCY BEFORE CALLING THIS ROUTINE.
C
C     CONSISTENCY = INTERPOLATION LAW = 1 THROUGH 6
C                 = NO DISCONTINUITIES, X1 = X2
C                 = ONLY POSITIVE VALUES IF LOG INTERPOLATION
C
C     06/02/10 - ADDED CONSISTENCY CHECKS FOR ALL PARAMETERS THAT
C                ARE NON-POSTITIVE AND REQUIRE TAKING THEIR LOG -
C                IN ALL SUCH CASE THIS ROUTINE SWITCHES TO LINEAR
C                (INTERP=2) INTERPOLATION.
C
C=======================================================================
      INCLUDE 'implicit.h'
      DATA ONED /1.0D+00/
      DATA ZEROD/0.0D+00/
C
C     FOR X1 = X2 OR Y1 = Y2 USE Y1.
C
      IF(X1.EQ.X2.OR.Y1.EQ.Y2) GO TO 10
C
C     SELECT INTERPOLATION METHOD.
C
C     IN ALL CASES THE RESULT (Y) IS THE WEIGHTED SUM OF
C     CONTRIBUTIONS FROM THE 2 ENDS OF THE INTERVAL.
C
      GO TO (10,20,30,40,50,60),INTERP
C-----1) HISTOGRAM - OR X1=X2 OR Y1=Y2 DEFINE Y = Y1
   10 TERPIT=Y1
      RETURN
C-----2) LIN X VS. LIN Y.
   20 WT2=(X-X1)/(X2-X1)
      WT1=ONED-WT2
      TERPIT=WT2*Y2+WT1*Y1
      RETURN
C-----3) LOG X VS. LIN Y.
   30 IF(X.LE.ZEROD.OR.X1.LE.ZEROD.OR.X2.LE.ZEROD) GO TO 20
      WT2=DLOG(X/X1)/DLOG(X2/X1)
      WT1=ONED-WT2
      TERPIT=WT2*Y2+WT1*Y1
      RETURN
C-----4) LIN X VS. LOG Y.
   40 IF(Y1.LE.ZEROD.OR.Y2.LE.ZEROD) GO TO 20
      WT2=(X-X1)/(X2-X1)
      WT1=ONED-WT2
      TERPIT=DEXP(WT2*DLOG(Y2)+WT1*DLOG(Y1))
      RETURN
C-----5) LOG X VS. LOG Y.
   50 IF(X.LE.ZEROD.OR.X1.LE.ZEROD.OR.X2.LE.ZEROD) GO TO 20
      IF(Y1.LE.ZEROD.OR.Y2.LE.ZEROD) GO TO 20
      WT2=DLOG(X/X1)/DLOG(X2/X1)
      WT1=ONED-WT2
      TERPIT=DEXP(WT2*DLOG(Y2)+WT1*DLOG(Y1))
      RETURN
C-----6) CHARGED PARTICLE THRESHOLDS...WARNING = THIS ASSUMES T = 0.0.
C-----06/02/09 - ORIGINAL DID NOT INCLUDE (A/E) TERM.
C
C     SIG = (A/E)*EXP[-B/SQRT(E - T)]
C     E*SIG = A*EXP[-B/SQRT(E-T)]
C     LOG(E*SIG) = LOG(A) - B/SQRT(E-T)
C     LOG(E*SIG) = WT2*LOG(X2*Y2) + WT1*LOG(X1*Y1)
C
C-----06/02/09 = USE LINEAR NEAR X OR Y <= 0
   60 IF(X.LE.ZEROD.OR.X1.LE.ZEROD.OR.X2.LE.ZEROD) GO TO 20
      IF(Y1.LE.ZEROD.OR.Y2.LE.ZEROD) GO TO 20
C-----OTHERWISE WEIGHT FOR E*SIG IS,
C-----WT2 = (1/SQRT(E)-1/SQRT(E1))/(1/SQRT(E2)-1/SQRT(E1))
      WT2=(ONED/DSQRT( X) - ONED/DSQRT(X1))/
     1    (ONED/DSQRT(X2) - ONED/DSQRT(X1))
      WT1=ONED-WT2
C-----LOG(E*SIG) = WT2*LOG(X2*Y2) + WT1*LOG(X1*Y1)
      TERPIT=DEXP(WT2*DLOG(X2*Y2)+WT1*DLOG(X1*Y1))/X
      RETURN
      END
      SUBROUTINE INCORE9(ZIN)
C=======================================================================
C
C     PURPOSE
C     =======
C     ROUND NUMBER TO FROM 5 TO 9 DIGITS OF ACCURACY.
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C                  No test for 4 or more digit exponent
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused
C                  problems on several different types of computers.
C
C     ARGUMENTS
C     =========
C     ZIN      = NUMBER OF BE ROUNDED (INPUT/OUTPUT)
C
C     METHOD
C     ======
C     COLUMNS            12345678901     ACCURACY
C     -------------------------------------------
C     0 TO 10^-9          1.2345E-12     5 DIGITS
C     10^-9 TO 10^-4      1.23456E-8     6 DIGITS
C     10^-4 TO 10^-3      .000123456     6 DIGITS
C     10^-3 TO 10^-2      .001234567     7 DIGITS
C     10^-2 TO 10^-1      .012345678     8 DIGITS
C     10^-1 TO 1          .123456789     9 DIGITS
C     1 TO 10^9           12345.6789     9 DIGITS
C     10^9 TO 10^10       1.23456E+9     6 DIGITS
C     10^10 >             1.2345E+12     5 DIGITS
C
C=======================================================================
      INCLUDE 'implicit.h'
      REAL*8 IN
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DIMENSION TENS(-99:99),ROUNDER(-99:99)
C-----ON FIRST CALL INITIALIZE POWERS OF 10
      DATA IPASS/0/
      IF(IPASS.NE.0) GO TO 50
      IPASS=1
      INMAN8 = 100000000
      INMAN9 = 1000000000
      TENS(0)=1.0D+00
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DO 10 I=1,99
      TENS( I)=TENS(I-1)*10.0D+00
      TENS(-I)=TENS(1-I)/10.0D+00
   10 ROUNDER(I)  = 5.001D-05
      DO 20 I=-99,-10
   20 ROUNDER(I)  = 5.001D-05
      DO 30 I=-9,-4
   30 ROUNDER(I)  = 5.001D-06
      ROUNDER(-3) = 5.001D-07
      ROUNDER(-2) = 5.001D-08
      ROUNDER(-1) = 5.001D-09
      ROUNDER( 0) = 5.001D-09
      DO 40 I=1,8
   40 ROUNDER(I)  = 5.001D-09
      ROUNDER(9)  = 5.001D-06
C
C     NO ROUNDING NECESSARY FOR ZERO - RETURN
C     OTHERWISE DEFINE SIGN AND ABSOLUTE VALUE.
C
   50 IF(ZIN.eq.0.0D+0) go to 160 ! no rounding if = 0
      IF(ZIN.gt.0.0D+0) go to 60
C-----NEGATIVE.
      ZSIGN=-1.0D+00
      Z=-ZIN
      GO TO 70
C-----POSITIVE.
   60 ZSIGN=1.0D+00
      Z=ZIN
C
C     DEFINE EXPONENT AND NORMALIZED MANTISSA
C
   70 IEXP=DLOG10(Z)
      IF(Z.LT.1.0D+00) IEXP = IEXP - 1
      IF(iabs(IEXP).gt.99) go to 160    ! no 2 digit exponent
      ZN=Z*TENS(-IEXP) + ROUNDER(IEXP)
      IF(ZN.eq.1.0D+00) go to 160 ! no rounding powers of 10
      IF(ZN.gt.1.0D+00) go to 80
      IEXP=IEXP-1                 ! < 1
      ZN=10.0D+00*ZN
      IF(iabs(IEXP).gt.99) go to 160    ! no 2 digit exponent
      Z = ZN*TENS(IEXP)
      GO TO 90
   80 IF(ZN.eq.10.0D+00) go to 160 ! no rounding powers of 10
      IF(ZN.lt.10.0D+00) go to 90
      IEXP=IEXP+1                ! > 10
      ZN=ZN/10.0D+00
      IF(iabs(IEXP).gt.99) go to 160    ! no 2 digit exponent
      Z = ZN*TENS(IEXP)
C
C     ZN IS NOW IN NORMAL FORM 1.23456789...
C
C-----------------------------------------------------------------------
C
C     TEST FOR SPECIAL RANGES = VERY LOW PROBABILITY
C
C-----------------------------------------------------------------------
   90 IF(Z.GE.1.0D+00) GO TO 110
C
C     IF EXTREMELY LOW ENERGY RANGE < 10^-10 USE 5 DIGITS
C
      IF(Z.LT.1.0D-09) GO TO 120
      IF(Z.GE.1.0D-04) GO TO 100
C-----10^-10 TO 10^-4 = 6 DIGITS
      IN = ZN*TENS(5)
      KEXP = IEXP-5
      GO TO 140
C-----10^-4 TO 1: 6 TO 9 DIGITS
  100 II = 9 + IEXP
      IF(iabs(II).gt.99) go to 160    ! no 2 digit exponent
      IN = ZN*TENS(II)
      KEXP = IEXP-II
      GO TO 140
C
C     HIGH ENERGY RANGE CHECK > 10^9
C
  110 IF(Z.LT.1.0D+09) GO TO 130
      IF(Z.GE.1.0D+10) GO TO 120
C
C     10^9 TO 10^10 = 6 DIGITS
C
      IN = ZN*TENS(5)
      KEXP = IEXP-5
      GO TO 140
C
C     EXTREME LOW AND HIGH ENERGY RANGE - USE 5 DIGITS
C
  120 IN = ZN*TENS(4)
      KEXP = IEXP-4
      GO TO 140
C-----------------------------------------------------------------------
C
C     NORMAL RANGE - 1 TO < 10^10 - USE 9 DIGITS = HIGH PROBABILITY
C
C-----------------------------------------------------------------------
  130 IN = ZN*TENS(8)
      KEXP = IEXP-8
C
C     IN IS NOW IN 9 DIGIT FORM 123456789
C     IF 10 DIGIT, DUE TO ROUNDING - DECREASE BY 10 AND INCREASE IEXP
C
c-----2014/4/14 - changed to INMAN9, instead of integer strings.
  140 IF(IN.lt.INMAN9) go to 150
      IN  = INMAN8
      IEXP = IEXP + 1
C
C     FLOAT 9 DIGIT AND RESTORE EXPONENT
C
  150 Z   = IN
      IF(iabs(KEXP).gt.99) go to 160    ! no 2 digit exponent
      ZIN = ZSIGN*Z*TENS(KEXP)
      RETURN
C
C     NO ROUNDING NECESSARY FOR 0
C
  160 RETURN
      END
      SUBROUTINE OUT9(ZIN,FIELD)
C=======================================================================
C
C     PURPOSE
C     =======
C     FORMAT NUMBER FOR OUTPUT TO INCLUDE AS MANY DIGITS OF
C     ACCURACY AS POSSIBLE.
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C                  No test for 4 or more digit exponent
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused
C                  problems on several different types of computers.
C     10/25/2014 - Changed from D to E exponential form to improve
C                  compatibility between computer languages
C
C     040923 - CHANGED ZLOW FROM 1.0D-4 TO 1.0D-3
C              NEAR 1.0D-3 PRECISION GOES FROM 999,999 TO 1,000,000
C              FOR SMOOTHEST TRANSITION FROM 6 TO 7 DIGITS
C
C     ARGUMENTS
C     =========
C     Z        = FLOATING POINT NUMBER OF BE OUTPUT (INPUT)
C     FIELD    = 11A1 CHARACTERS TO OUTPUT          (OUTPUT)
C
C     METHOD
C     ======
C     COLUMNS            12345678901     ACCURACY
C     -------------------------------------------
C     0 TO 10^-9          1.2345E-12     5 DIGITS
C     10^-9 TO 10^-3      1.23456E-8     6 DIGITS
C     10^-3 TO 10^-2      .001234567     7 DIGITS
C     10^-2 TO 10^-1      .012345678     8 DIGITS
C     10^-1 TO 1          .123456789     9 DIGITS
C     1 TO 10^9           12345.6789     9 DIGITS
C     10^9 TO 10^10       1.23456E+9     6 DIGITS
C     10^10 >             1.2345E+12     5 DIGITS
C
C     OUTPUT WILL BE IN 11 COLUMN FORMAT
C
C     WARNING - THIS IS NOT A GENERAL ROUNDING ROUTINE WHICH WILL WORK
C               FOR ROUNDING TO ANY NUMBER OF DIGITS - IT WILL ONLY
C               WORK PROPERLY FOR THE RANGES INDICATED ABOVE, FOR
C               11 COLUMN OUTPUT.
C
C=======================================================================
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      save
c-----2017/5/3 - Add INNEXT = must be INTEGER*8
c-----INTEGER*8 for long integers
      integer*8 INMANT,INMAN9,INMAN8,INMAN7,INMAN6,INMAN5,INMAN4,INNEXT
c-----2017/5/3 - Add INNEXT = must be INTEGER*8
      CHARACTER*1 FIELD,DIGITS,ZEROH
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DIMENSION DIGITS(0:9),FIELD(11),ZEROH(11),
     1 TENS(-99:99),ROUNDER(-99:99)
      DATA DIGITS/
     1 '0','1','2','3','4','5','6','7','8','9'/
C-----RETURN FOR = 0
      DATA ZEROH/
     1 ' ','0','.','0',' ',' ',' ',' ',' ',' ',' '/
C-----LOWER TRANSITION POINT FROM 7 TO 9 DIGIT OUTPUT
      DATA ZLOW/1.0D-03/
C-----UPPER TRANSITION POINT FROM 9 TO 7 DIGIT OUTPUT
      DATA ZHIGH/1.0D+09/
      DATA TENTH/1.0D-01/
c-----------------------------------------------------------------------
c
c     ON FIRST CALL INITIALIZE POWERS OF 10
c
c-----------------------------------------------------------------------
      DATA IPASS/0/
      IF(IPASS.NE.0) GO TO 50
      IPASS=1
c-----INTEGER*8 for long integers
      INMAN4  = 10000
      INMAN5  = 10*INMAN4
      INMAN6  = 10*INMAN5
      INMAN7  = 10*INMAN6
      INMAN8  = 10*INMAN7
      INMAN9  = 10*INMAN8
      TENS(0)=1.0D+00
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DO 10 I=1,99
      TENS( I)=TENS(I-1)*10.0D+00
      TENS(-I)=TENS(1-I)/10.0D+00
   10 ROUNDER(I)  = 5.001D-05
      DO 20 I=-99,-10
   20 ROUNDER(I)  = 5.001D-05
      DO 30 I=-9,-4
   30 ROUNDER(I)  = 5.001D-06
      ROUNDER(-3) = 5.001D-07
      ROUNDER(-2) = 5.001D-08
      ROUNDER(-1) = 5.001D-09
      ROUNDER( 0) = 5.001D-09
      DO 40 I=1,8
   40 ROUNDER(I)  = 5.001D-09
      ROUNDER(9)  = 5.001D-06
c-----------------------------------------------------------------------
C
C     IMMEDIATELY RETURN 0.00000+00.
C
c-----------------------------------------------------------------------
C-----02/14/04 - ADDED, JUST IN CASE
   50 IF(DABS(ZIN).le.0.0D+0) GO TO 60
      IF(DABS(ZIN).lt.1.0D-99) GO TO 60    ! Only 2 digit exponents
      IF(DABS(ZIN).gt.1.0D+99) GO TO 60
c----- ZIN = 0 handled above
      IF(ZIN.lt.0.0d+0) go to 80
      go to 90
c-----Return 0
   60 DO 70 I=1,11
   70 FIELD(I)=ZEROH(I)
      RETURN
c-----------------------------------------------------------------------
C
C     DEFINE SIGN OF MANTISSA AND ABSOLUTE MANTISSA
C
c-----------------------------------------------------------------------
C-----NEGATIVE.
   80 FIELD(1)='-'
      Z=-ZIN
      GO TO 100
C-----POSITIVE.
   90 FIELD(1)=' '
      Z=ZIN
c-----------------------------------------------------------------------
c
c     DEFINE EXPONENT AND NORMALIZED MANTISSA
c
c-----------------------------------------------------------------------
  100 IEXP=DLOG10(Z)
      IF(Z.LT.1.0D+00) IEXP = IEXP - 1
c-----11/22/2013 - Decide here on F or E Output
c-----11/26/2013 - Restrict to 2 digit exponent - may change by 1 or 2
      if(IABS(IEXP).gt.97) go to 60
      if(ZIN.lt.0.0D+0) then
      ZN=Z*TENS(-IEXP) + ROUNDER(IEXP)       ! Standard
      else
      if(Z.LE.ZLOW.OR.Z.GE.ZHIGH) then       ! F or E Format?
      ZN=Z*TENS(-IEXP) + TENTH*ROUNDER(IEXP) ! Extra Digit
      else
      ZN=Z*TENS(-IEXP) + ROUNDER(IEXP)       ! Standard
      endif
      endif
      IF(ZN.eq.1.0D+00) go to 120            ! Test rounding underflow
      IF(ZN.gt.1.0D+00) go to 110            ! Test rounding underflow
      IEXP=IEXP-1         ! MUST be < 1
      GO TO 120
  110 IF(ZN.lt.10.0D+00) go to 120           ! Test rounding overflow
      IEXP=IEXP+1         ! MUST be >= 10
      ZN=ZN/10.0D+00
  120 Z = ZN*TENS(IEXP)
c-----------------------------------------------------------------------
C
C     SELECT F OR E FORMAT
C
c-----------------------------------------------------------------------
      IF(Z.LE.ZLOW.OR.Z.GE.ZHIGH) GO TO 150
c-----------------------------------------------------------------------
C
C     F FORMAT
C
C     12345678901
C      X.XXXXXXXX = 9 DIGITS
C      .001234567
C      123456789.
C
c-----------------------------------------------------------------------
C-----DEFINE 6 TO 9 DIGIT MANTISSA WITH ROUNDING
      IPOWER=8-IEXP
      IF(IEXP.LT.0) IPOWER=8
      INMANT=Z*TENS(IPOWER)
C-----CHECK FOR OVERFLOW DUE TO ROUNDING
      if(INMANT.lt.INMAN9) go to 130
      INMANT=INMAN8
      IEXP=IEXP+1
C-----DECIMAL POINT.
  130 IDOT=3+IEXP
      IF(IDOT.LE.2) THEN
C----- IF < 1, MOVE DECIMAL POINT TO COLUMN 2 AND ADD A DIGIT
      IDOT=2
      INMANT=Z*10.0D+00*TENS(IPOWER)
      ENDIF
      FIELD(IDOT)='.'
C-----MANTISSA - LAST DIGIT TO FIRST.
      II=11
      DO 140 I=2,11
      IF(II.EQ.IDOT) GO TO 140
      INNEXT=INMANT/10
      I3=INMANT-10*INNEXT
      FIELD(II)=DIGITS(I3)
      INMANT=INNEXT
  140 II=II-1
      RETURN
c-----------------------------------------------------------------------
C
C     E FORMAT
C
C     12345678901
C      X.XXXXE+NN = 5 DIGITS
C      X.XXXXXE+N = 6 DIGITS
C
C     11/22/2013 - If not negative, use first column
C     12345678901
C     X.XXXXXE+NN = 6 DIGITS
C     X.XXXXXXE+N = 7 DIGITS
c
c-----------------------------------------------------------------------
C-----Negative?
  150 IF(ZIN.lt.0.0d+0) go to 170
c
c     POsitive. Use first column - Decimal point is always in column 2
c
      FIELD(2)='.'
      KDOT    = 2
      ISTART  = 1
      IF(IABS(IEXP).GE.10) GO TO 160
      ID=8                               ! 1 Digit exponent
      INMANT=(1.0D+06)*ZN
      IF(INMANT.lt.INMAN7) go to 190
      INMANT=INMAN6
      IEXP=IEXP+1
      IF(IABS(IEXP).LT.10) GO TO 190
  160 ID=7                               ! 2 Digit exponent
      INMANT=(1.0D+05)*ZN
C-----CHECK FOR OVERFLOW DUE TO ROUNDING
      IF(INMANT.lt.INMAN6) go to 190
      INMANT=INMAN5
      IEXP=IEXP+1
c
c     Negative Number - Cannot use first column
c
C-----DECIMAL POINT IS ALWAYS IN COLUMN 3
  170 FIELD(3)='.'
      KDOT    = 3
      ISTART  = 2
      IF(IABS(IEXP).GE.10) GO TO 180
      ID=8                                 ! 1 Digit Exponent
      INMANT=(1.0D+05)*ZN
      IF(INMANT.lt.INMAN6) go to 190
      INMANT=INMAN5
      IEXP=IEXP+1
      IF(IABS(IEXP).LT.10) GO TO 190
  180 ID=7                                 ! 2 Digit Exponent
      INMANT=(1.0D+04)*ZN
      IF(INMANT.lt.INMAN5) go to 190
      INMANT=INMAN4
      IEXP=IEXP+1
C
C     DEFINE MANTISSA
C
  190 IEXPS=ID+1
      II=ID
      DO 200 I=ISTART,ID
      IF(II.EQ.KDOT) GO TO 200
      INNEXT=INMANT/10
      I3=INMANT-10*INNEXT
      FIELD(II)=DIGITS(I3)
      INMANT=INNEXT
  200 II=II-1
C
C     E
C
      FIELD(IEXPS) = 'E'
      IEXPS = IEXPS + 1
C
C     SIGN OF EXPONENT
C
      IF(IEXP.ge.0) go to 210
      IEXP=-IEXP
      FIELD(IEXPS)='-'
      GO TO 220
  210 FIELD(IEXPS)='+'
C
C     EXPONENT
C
  220 IF(IEXP.lt.10) go to 230
      KEXP=IEXP/10
      FIELD(10)=DIGITS(KEXP)
      IEXP=MOD(IEXP,10)
  230 FIELD(11)=DIGITS(IEXP)
c
c     If using column 1 but mantissa ends in 0, move mantissa right.
c
      if(KDOT.eq.2.and.FIELD(ID).eq.'0') then
      do k=ID-1,1,-1
      FIELD(k+1)=FIELD(k)
      enddo
      FIELD(1) = ' '
      endif
      RETURN
      END
      SUBROUTINE OUT9G(ZIN,FIELD)
C=======================================================================
C
C     PURPOSE
C     =======
C     20915/7/28 - same as OUT9, but column 1 is always '-' or blank.
C                  Used by GROUPIE to simplify listings.
C
C     FORMAT NUMBER FOR OUTPUT TO INCLUDE AS MANY DIGITS OF
C     ACCURACY AS POSSIBLE.
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C                  No test for 4 or more digit exponent
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused
C                  problems on several different types of computers.
C     10/25/2014 - Changed from D to E exponential form to improve
C                  compatibility between computer languages
C
C     040923 - CHANGED ZLOW FROM 1.0D-4 TO 1.0D-3
C              NEAR 1.0D-3 PRECISION GOES FROM 999,999 TO 1,000,000
C              FOR SMOOTHEST TRANSITION FROM 6 TO 7 DIGITS
C
C     ARGUMENTS
C     =========
C     Z        = FLOATING POINT NUMBER OF BE OUTPUT (INPUT)
C     FIELD    = 11A1 CHARACTERS TO OUTPUT          (OUTPUT)
C
C     METHOD
C     ======
C     COLUMNS            12345678901     ACCURACY
C     -------------------------------------------
C     0 TO 10^-9          1.2345E-12     5 DIGITS
C     10^-9 TO 10^-3      1.23456E-8     6 DIGITS
C     10^-3 TO 10^-2      .001234567     7 DIGITS
C     10^-2 TO 10^-1      .012345678     8 DIGITS
C     10^-1 TO 1          .123456789     9 DIGITS
C     1 TO 10^9           12345.6789     9 DIGITS
C     10^9 TO 10^10       1.23456E+9     6 DIGITS
C     10^10 >             1.2345E+12     5 DIGITS
C
C     OUTPUT WILL BE IN 11 COLUMN FORMAT
C
C     WARNING - THIS IS NOT A GENERAL ROUNDING ROUTINE WHICH WILL WORK
C               FOR ROUNDING TO ANY NUMBER OF DIGITS - IT WILL ONLY
C               WORK PROPERLY FOR THE RANGES INDICATED ABOVE, FOR
C               11 COLUMN OUTPUT.
C
C=======================================================================
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      save
c-----INTEGER*8 for long integers
c-----2017/5/3 - Add INNEXT = must be INTEGER*8
      integer*8 INMANT,INMAN9,INMAN8,INMAN7,INMAN6,INMAN5,INMAN4,INNEXT
c-----2017/5/3 - Add INNEXT = must be INTEGER*8
      CHARACTER*1 FIELD,DIGITS,ZEROH
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DIMENSION DIGITS(0:9),FIELD(11),ZEROH(11),
     1 TENS(-99:99),ROUNDER(-99:99)
      DATA DIGITS/
     1 '0','1','2','3','4','5','6','7','8','9'/
C-----RETURN FOR = 0
      DATA ZEROH/
     1 ' ','0','.','0',' ',' ',' ',' ',' ',' ',' '/
C-----LOWER TRANSITION POINT FROM 7 TO 9 DIGIT OUTPUT
      DATA ZLOW/1.0D-03/
C-----UPPER TRANSITION POINT FROM 9 TO 7 DIGIT OUTPUT
      DATA ZHIGH/1.0D+09/
c-----------------------------------------------------------------------
c
c     ON FIRST CALL INITIALIZE POWERS OF 10
c
c-----------------------------------------------------------------------
      DATA IPASS/0/
      IF(IPASS.NE.0) GO TO 50
      IPASS=1
c-----INTEGER*8 for long integers
      INMAN4  = 10000
      INMAN5  = 10*INMAN4
      INMAN6  = 10*INMAN5
      INMAN7  = 10*INMAN6
      INMAN8  = 10*INMAN7
      INMAN9  = 10*INMAN8
      TENS(0)=1.0D+00
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DO 10 I=1,99
      TENS( I)=TENS(I-1)*10.0D+00
      TENS(-I)=TENS(1-I)/10.0D+00
   10 ROUNDER(I)  = 5.001D-05
      DO 20 I=-99,-10
   20 ROUNDER(I)  = 5.001D-05
      DO 30 I=-9,-4
   30 ROUNDER(I)  = 5.001D-06
      ROUNDER(-3) = 5.001D-07
      ROUNDER(-2) = 5.001D-08
      ROUNDER(-1) = 5.001D-09
      ROUNDER( 0) = 5.001D-09
      DO 40 I=1,8
   40 ROUNDER(I)  = 5.001D-09
      ROUNDER(9)  = 5.001D-06
c-----------------------------------------------------------------------
C
C     IMMEDIATELY RETURN 0.00000+00.
C
c-----------------------------------------------------------------------
C-----02/14/04 - ADDED, JUST IN CASE
   50 IF(DABS(ZIN).le.0.0D+0) GO TO 60
      IF(DABS(ZIN).lt.1.0D-99) GO TO 60    ! Only 2 digit exponents
      IF(DABS(ZIN).gt.1.0D+99) GO TO 60
c----- ZIN = 0 handled above
      IF(ZIN.lt.0.0d+0) go to 80
      go to 90
c-----Return 0
   60 DO 70 I=1,11
   70 FIELD(I)=ZEROH(I)
      RETURN
c-----------------------------------------------------------------------
C
C     DEFINE SIGN OF MANTISSA AND ABSOLUTE MANTISSA
C
c-----------------------------------------------------------------------
C-----NEGATIVE.
   80 FIELD(1)='-'
      Z=-ZIN
      GO TO 100
C-----POSITIVE.
   90 FIELD(1)=' '
      Z=ZIN
c-----------------------------------------------------------------------
c
c     DEFINE EXPONENT AND NORMALIZED MANTISSA
c
c-----------------------------------------------------------------------
  100 IEXP=DLOG10(Z)
      IF(Z.LT.1.0D+00) IEXP = IEXP - 1
c-----11/22/2013 - Decide here on F or E Output
c-----11/26/2013 - Restrict to 2 digit exponent - may change by 1 or 2
      if(IABS(IEXP).gt.97) go to 60
      ZN=Z*TENS(-IEXP) + ROUNDER(IEXP)       ! Standard
      IF(ZN.eq.1.0D+00) go to 120            ! Test rounding underflow
      IF(ZN.gt.1.0D+00) go to 110            ! Test rounding underflow
      IEXP=IEXP-1         ! MUST be < 1
      GO TO 120
  110 IF(ZN.lt.10.0D+00) go to 120           ! Test rounding overflow
      IEXP=IEXP+1         ! MUST be >= 10
      ZN=ZN/10.0D+00
  120 Z = ZN*TENS(IEXP)
c-----------------------------------------------------------------------
C
C     SELECT F OR E FORMAT
C
c-----------------------------------------------------------------------
      IF(Z.LE.ZLOW.OR.Z.GE.ZHIGH) GO TO 150
c-----------------------------------------------------------------------
C
C     F FORMAT
C
C     12345678901
C      X.XXXXXXXX = 9 DIGITS
C      .001234567
C      123456789.
C
c-----------------------------------------------------------------------
C-----DEFINE 6 TO 9 DIGIT MANTISSA WITH ROUNDING
      IPOWER=8-IEXP
      IF(IEXP.LT.0) IPOWER=8
      INMANT=Z*TENS(IPOWER)
C-----CHECK FOR OVERFLOW DUE TO ROUNDING
      if(INMANT.lt.INMAN9) go to 130
      INMANT=INMAN8
      IEXP=IEXP+1
C-----DECIMAL POINT.
  130 IDOT=3+IEXP
      IF(IDOT.LE.2) THEN
C----- IF < 1, MOVE DECIMAL POINT TO COLUMN 2 AND ADD A DIGIT
      IDOT=2
      INMANT=Z*10.0D+00*TENS(IPOWER)
      ENDIF
      FIELD(IDOT)='.'
C-----MANTISSA - LAST DIGIT TO FIRST.
      II=11
      DO 140 I=2,11
      IF(II.EQ.IDOT) GO TO 140
      INNEXT=INMANT/10
      I3=INMANT-10*INNEXT
      FIELD(II)=DIGITS(I3)
      INMANT=INNEXT
  140 II=II-1
      RETURN
c-----------------------------------------------------------------------
C
C     E FORMAT
C
C     12345678901
C      X.XXXXE+NN = 5 DIGITS
C      X.XXXXXE+N = 6 DIGITS
C
C     11/22/2013 - If not negative, use first column
C     12345678901
C     X.XXXXXE+NN = 6 DIGITS
C     X.XXXXXXE+N = 7 DIGITS
c
c-----------------------------------------------------------------------
c
c     Do not cannot use first column
c
C-----DECIMAL POINT IS ALWAYS IN COLUMN 3
  150 FIELD(3)='.'
      KDOT    = 3
      ISTART  = 2
      IF(IABS(IEXP).GE.10) GO TO 160
      ID=8                                 ! 1 Digit Exponent
      INMANT=(1.0D+05)*ZN
      IF(INMANT.lt.INMAN6) go to 170
      INMANT=INMAN5
      IEXP=IEXP+1
      IF(IABS(IEXP).LT.10) GO TO 170
  160 ID=7                                 ! 2 Digit Exponent
      INMANT=(1.0D+04)*ZN
      IF(INMANT.lt.INMAN5) go to 170
      INMANT=INMAN4
      IEXP=IEXP+1
C
C     DEFINE MANTISSA
C
  170 IEXPS=ID+1
      II=ID
      DO 180 I=ISTART,ID
      IF(II.EQ.KDOT) GO TO 180
      INNEXT=INMANT/10
      I3=INMANT-10*INNEXT
      FIELD(II)=DIGITS(I3)
      INMANT=INNEXT
  180 II=II-1
C
C     E
C
      FIELD(IEXPS) = 'E'
      IEXPS = IEXPS + 1
C
C     SIGN OF EXPONENT
C
      IF(IEXP.ge.0) go to 190
      IEXP=-IEXP
      FIELD(IEXPS)='-'
      GO TO 200
  190 FIELD(IEXPS)='+'
C
C     EXPONENT
C
  200 IF(IEXP.lt.10) go to 210
      KEXP=IEXP/10
      FIELD(10)=DIGITS(KEXP)
      IEXP=MOD(IEXP,10)
  210 FIELD(11)=DIGITS(IEXP)
      RETURN
      END
      SUBROUTINE OUT10(ZIN,FIELD)
C=======================================================================
C
C     PURPOSE
C     =======
C     FORMAT NUMBER FOR OUTPUT TO INCLUDE AS MANY DIGITS OF
C     ACCURACY AS POSSIBLE.
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C                  No test for 4 or more digit exponent
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused
C                  problems on several different types of computers.
C     10/25/2014 - Changed from D to E exponential form to improve
C                  compatibility between computer languages
C
C     040923 - CHANGED ZLOW FROM 1.0D-4 TO 1.0D-3
C              NEAR 1.0D-3 PRECISION GOES FROM 999,999 TO 1,000,000
C              FOR SMOOTHEST TRANSITION FROM 6 TO 7 DIGITS
C
C     ARGUMENTS
C     =========
C     Z        = FLOATING POINT NUMBER OF BE OUTPUT (INPUT)
C     FIELD    = 11A1 CHARACTERS TO OUTPUT          (OUTPUT)
C
C     METHOD
C     ======
C     COLUMNS            12345678901     ACCURACY
C     -------------------------------------------
C     0 TO 10^-9         1.23456E-12     6 DIGITS
C     10^-9 TO 10^-3     1.234567E-8     7 DIGITS
C     10^-3 TO 10^-2     .0012345678     8 DIGITS
C     10^-2 TO 10^-1     .0123456789     9 DIGITS
C     10^-1 TO 1         .1234567891    10 DIGITS
C     1 TO 10^9          12345.67891    10 DIGITS
C     10^9 TO 10^10      1.234567E+9     7 DIGITS
C     10^10 >            1.23456E+12     6 DIGITS
C
C     OUTPUT WILL BE IN 11 COLUMN FORMAT
C
C     WARNING - THIS IS NOT A GENERAL ROUNDING ROUTINE WHICH WILL WORK
C               FOR ROUNDING TO ANY NUMBER OF DIGITS - IT WILL ONLY
C               WORK PROPERLY FOR THE RANGES INDICATED ABOVE, FOR
C               11 COLUMN OUTPUT.
C
C=======================================================================
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      save
c-----INTEGER*8 for long integers
c-----2017/5/3 - Add INNEXT = must be INTEGER*8
      integer*8 INMANT,INMAN10,INMAN9,INMAN7,INMAN6,INMAN5,INMAN4,INNEXT
c-----2017/5/3 - Add INNEXT = must be INTEGER*8
      CHARACTER*1 FIELD,DIGITS,ZEROH
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DIMENSION DIGITS(0:9),FIELD(11),ZEROH(11),
     1 TENS(-99:99),ROUNDER(-99:99)
      DATA DIGITS/
     1 '0','1','2','3','4','5','6','7','8','9'/
C-----RETURN FOR = 0
      DATA ZEROH/
     1 ' ','0','.','0',' ',' ',' ',' ',' ',' ',' '/
C-----LOWER TRANSITION POINT FROM 7 TO 9 DIGIT OUTPUT
      DATA ZLOW/1.0D-03/
C-----UPPER TRANSITION POINT FROM 9 TO 7 DIGIT OUTPUT
      DATA ZHIGH/1.0D+09/
      DATA TENTH/1.0D-01/
c-----------------------------------------------------------------------
c
c     ON FIRST CALL INITIALIZE POWERS OF 10
c
c-----------------------------------------------------------------------
      DATA IPASS/0/
      IF(IPASS.NE.0) GO TO 50
      IPASS=1
c-----INTEGER*8 for long integers
      INMAN4  = 10000
      INMAN5  = 10*INMAN4
      INMAN6  = 10*INMAN5
      INMAN7  = 10*INMAN6
      INMAN9  = 100*INMAN7
      INMAN10 = 10*INMAN9  ! WARNING - this exceeds 32 bit integer
      TENS(0)=1.0D+00
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DO 10 I=1,99
      TENS( I)=TENS(I-1)*10.0D+00
      TENS(-I)=TENS(1-I)/10.0D+00
   10 ROUNDER(I)  = 5.001D-05
      DO 20 I=-99,-10
   20 ROUNDER(I)  = 5.001D-05
      DO 30 I=-9,-4
   30 ROUNDER(I)  = 5.001D-06
      ROUNDER(-3) = 5.001D-07
      ROUNDER(-2) = 5.001D-08
      ROUNDER(-1) = 5.001D-09
      ROUNDER( 0) = 5.001D-09
      DO 40 I=1,8
   40 ROUNDER(I)  = 5.001D-09
      ROUNDER(9)  = 5.001D-06
c-----------------------------------------------------------------------
C
C     IMMEDIATELY RETURN 0.00000+00.
C
c-----------------------------------------------------------------------
C-----02/14/04 - ADDED, JUST IN CASE
   50 IF(DABS(ZIN).le.0.0D+0) GO TO 60
      IF(DABS(ZIN).lt.1.0D-99) GO TO 60    ! Only 2 digit expoents
      IF(DABS(ZIN).gt.1.0D+99) GO TO 60
c----- ZIN = 0 handled above
      IF(ZIN.lt.0.0d+0) go to 80
      go to 90
   60 DO 70 I=1,11
   70 FIELD(I)=ZEROH(I)
      RETURN
c-----------------------------------------------------------------------
C
C     DEFINE SIGN OF MANTISSA AND ABSOLUTE MANTISSA
C
c-----------------------------------------------------------------------
c-----< 0 = use OUT9
   80 CALL OUT9(ZIN,FIELD)
      RETURN
C-----POSITIVE.
   90 FIELD(1)=' '
      Z=ZIN
c-----------------------------------------------------------------------
c
c     DEFINE EXPONENT AND NORMALIZED MANTISSA
c
c-----------------------------------------------------------------------
      IEXP=DLOG10(Z)
c-----11/26/2013 - Restrict to 2 digit exponent - may change by 1 or 2
      if(IABS(IEXP).gt.97) go to 60
      IF(Z.LT.1.0D+00) IEXP = IEXP - 1
c-----11/22/2013 - Decide here on F or E Output
      if(ZIN.lt.0.0D+0) then
      ZN=Z*TENS(-IEXP) + ROUNDER(IEXP)       ! Standard
      else
      ZN=Z*TENS(-IEXP) + TENTH*ROUNDER(IEXP) ! Extra Digit
      endif
      IF(ZN.eq.1.0D+00) go to 110
      IF(ZN.gt.1.0D+00) go to 100
      IEXP=IEXP-1               ! MUST be < 1
      go to 110
  100 IF(ZN.lt.10.0D+00) go to 110
      IEXP=IEXP+1               ! MUST be >= 10
      ZN=ZN/10.0D+00
  110 Z = ZN*TENS(IEXP)
c-----------------------------------------------------------------------
C
C     SELECT F OR E FORMAT
C
c-----------------------------------------------------------------------
      IF(Z.LE.ZLOW.OR.Z.GE.ZHIGH) GO TO 140
c-----------------------------------------------------------------------
C
C     F FORMAT
C
C     12345678901
C     X.XXXXXXXXX = 10 DIGITS
C     .0012345678
C     1234567891.
C
c-----------------------------------------------------------------------
C-----DEFINE 7 TO 10 DIGIT MANTISSA WITH ROUNDING
      IPOWER=9-IEXP
      IF(IEXP.LT.0) IPOWER=9
      INMANT=Z*TENS(IPOWER)
C-----CHECK FOR OVERFLOW DUE TO ROUNDING
      IF(INMANT.lt.INMAN10) go to 120
      INMANT=INMAN9
      IEXP=IEXP+1
C-----DECIMAL POINT.
  120 IDOT=2+IEXP
      IF(IDOT.LE.1) THEN
C----- IF < 1, MOVE DECIMAL POINT TO COLUMN 1 AND ADD A DIGIT
      IDOT=1
      INMANT=Z*10.0D+00*TENS(IPOWER)
      ENDIF
      FIELD(IDOT)='.'
C-----MANTISSA - LAST DIGIT TO FIRST.
      II=11
      DO 130 I=1,11
      IF(II.EQ.IDOT) GO TO 130
      INNEXT=INMANT/10
      I3=INMANT-10*INNEXT
      FIELD(II)=DIGITS(I3)
      INMANT=INNEXT
  130 II=II-1
c
c     If 10 digit ends in 0, shift right = leave column 1 blank
c
      if(FIELD(1).ne.'-'.and.FIELD(11).eq.'0') then
      do k=10,1,-1
      FIELD(k+1)=FIELD(k)
      enddo
      FIELD(1) = ' '
      endif
      RETURN
c-----------------------------------------------------------------------
C
C     E FORMAT
C
C     12345678901
C      X.XXXXE+NN = 5 DIGITS
C      X.XXXXXE+N = 6 DIGITS
c
c     11/22/2013 - If first column is blank use it.
C
C     12345678901
C     X.XXXXXE+NN = 6 DIGITS
C     X.XXXXXXE+N = 7 DIGITS
C
C==============================================================
C-----Negative?
  140 IF(ZIN.lt.0.0d+0) go to 160
c
c     POsitive. Use first column - Decimal point is always in column 2
c
      FIELD(2)='.'
      KDOT    = 2
      ISTART  = 1
      IF(IABS(IEXP).GE.10) GO TO 150
      ID=8                               ! 1 Digit exponent
      INMANT=(1.0D+06)*ZN
      IF(INMANT.lt.INMAN7) go to 180
      INMANT=INMAN6
      IEXP=IEXP+1
      IF(IABS(IEXP).LT.10) GO TO 180
  150 ID=7                               ! 2 Digit exponent
      INMANT=(1.0D+05)*ZN
C-----CHECK FOR OVERFLOW DUE TO ROUNDING
      IF(INMANT.lt.INMAN6) go to 180
      INMANT=INMAN5
      IEXP=IEXP+1
c
c     Negative Number - Cannot use first column
c
C-----DECIMAL POINT IS ALWAYS IN COLUMN 3
  160 FIELD(3)='.'
      KDOT    = 3
      ISTART  = 2
      IF(IABS(IEXP).GE.10) GO TO 170
      ID=8                                 ! 1 Digit Exponent
      INMANT=(1.0D+05)*ZN
      IF(INMANT.lt.INMAN6) go to 180
      INMANT=INMAN5
      IEXP=IEXP+1
      IF(IABS(IEXP).LT.10) GO TO 180
  170 ID=7                                 ! 2 Digit Exponent
      INMANT=(1.0D+04)*ZN
      IF(INMANT.lt.INMAN5) go to 180
      INMANT=INMAN4
      IEXP=IEXP+1
C
C     DEFINE MANTISSA
C
  180 IEXPS=ID+1
      II=ID
      DO 190 I=ISTART,ID
      IF(II.EQ.KDOT) GO TO 190
      INNEXT=INMANT/10
      I3=INMANT-10*INNEXT
      FIELD(II)=DIGITS(I3)
      INMANT=INNEXT
  190 II=II-1
C
C     E
C
      FIELD(IEXPS) = 'E'
      IEXPS = IEXPS + 1
C
C     SIGN OF EXPONENT
C
      IF(IEXP.ge.0) go to 200
      IEXP=-IEXP
      FIELD(IEXPS)='-'
      GO TO 210
  200 FIELD(IEXPS)='+'
C
C     EXPONENT
C
  210 IF(IEXP.lt.10) go to 220
      KEXP=IEXP/10
      FIELD(10)=DIGITS(KEXP)
      IEXP=MOD(IEXP,10)
  220 FIELD(11)=DIGITS(IEXP)
c
c     If using column 1 but mantissa ends in 0, move mantissa right.
c
      if(KDOT.eq.2.and.FIELD(ID).eq.'0') then
      do k=ID-1,1,-1
      FIELD(k+1)=FIELD(k)
      enddo
      FIELD(1) = ' '
      endif
      RETURN
      END
      SUBROUTINE IN9(E,FIELD)
C=======================================================================
C
C     PURPOSE
C     =======
C     CONVERT FROM HOLLERITH TO FLOATING POINT.
C     12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENT (I Thank Viktor Zerkin)
C                  More than 3 digit exponent = ERROR - return 0.0
C     06/12/2013 - BACK tO 2 digit exponent - 3 digits caused
C                  problems on several different types of computers.
C
C     ARGUMENTS
C     =========
C     E       = FLOATING POINT NUMBER (OUTPUT)
C     FIELD   = 11A1 CHARACTER STRING (INPUT)
C
C     METHOD
C     ======
C     FIELD IS A STRING OF 11 CHARACTERS.
C     IT IS CONVERTED INTO A FLOATING POINR NUMBER (E)
C
C=======================================================================
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      save
      CHARACTER*1 MESS,DIGIT,FIELD,IFIELD
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENTS
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DIMENSION FIELD(11),TENS(-99:99),DIGIT(0:9),MESS(11),XDIG(0:9)
      DATA MESS/' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/
      DATA DIGIT/'0','1','2','3','4','5','6','7','8','9'/
      DATA XDIG/0.0D+00,1.0D+00,2.0D+00,3.0D+00,4.0D+00,
     1          5.0D+00,6.0D+00,7.0D+00,8.0D+00,9.0D+00/
C-----ON FIRST CALL DEFINE POWERS OF 10
      DATA IPASS/0/
      IF(IPASS.NE.0) GO TO 20
      IPASS=1
      TENS(0)=1.0D+00
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENTS
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      DO 10 I=1,99
      TENS( I)=TENS(I-1)*10.0D+00
   10 TENS(-I)=TENS(1-I)/10.0D+00
C
C     TRANSLATE MANTISSA.
C
C-----SKIP LEADING BLANK CHARACTERS.
   20 DO 30 I=1,11
      IF(FIELD(I).NE.' ') GO TO 40
   30 CONTINUE
C-----FIELD IS COMPLETELY BLANK. RETURN ZERO.
      E=0.0D+00
      RETURN
C-----INITIALIZE SIGN TO PLUS AND THEN CHECK FOR LEADING MINUS SIGN.
   40 SIGN=1.0D+00
C-----06/02/09 - ADDED TEST FOR LEADING + SIGN.
      IF(FIELD(I).NE.'-') GO TO 50
      SIGN=-1.0D+00
      I=I+1
c-----added leading +
      GO TO 60
   50 IF(FIELD(I).EQ.'+') I = I + 1
C-----INITIALIZE FIXED POINT INPUT FIELD AND POSITION OF DECIMAL POINT.
   60 X9IN=0.0
      IPT=-20
      IMZERO=0
C-----SCAN REMAINDER OF MANTISSA.
      DO 120 J=I,11
      IFIELD=FIELD(J)
C-----SCAN FOR DIGIT OR DECIMAL POINT (WHICH ARE PART OF MANTISSA).
      DO 70 K=0,9
      IF(IFIELD.EQ.DIGIT(K)) GO TO 90
   70 CONTINUE
      IF(IFIELD.NE.'.') GO TO 80
      IPT=0
      GO TO 120
C-----SCAN FOR BLANK (WHICH ENDS MANTISSA).
   80 IF(IFIELD.EQ.' ') GO TO 130
C-----SCAN FOR e, E, d, D, - OR + (WHICH BEGINS EXPONENT).
      IF(IFIELD.EQ.'e'.OR.IFIELD.EQ.'E') GO TO 160
      IF(IFIELD.EQ.'d'.OR.IFIELD.EQ.'D') GO TO 160
      IF(IFIELD.EQ.'-') GO TO 190
      IF(IFIELD.EQ.'+') GO TO 170
C-----ERROR. CANNOT IDENTIFY CHARACTER.
      GO TO 240
C-----DIGIT FOUND. SAVE TRAILING ZEROES AFTER DECIMAL POINT.
   90 IF(IPT.LT.0) GO TO 110
      IF(K.NE.0) GO TO 100
C-----SAVE TRAILING ZEROES.
      IMZERO=IMZERO+1
      GO TO 120
  100 IF(IMZERO.LE.0) GO TO 110
C-----INSERT ZEROES BEFORE NEXT NUMBER.
      X9IN=10.0D+00*X9IN
      IPT=IPT+1
      IMZERO=IMZERO-1
      GO TO 100
C-----DIGIT FOUND. INCREMENT FIXED POINT EQUIVALENT AND DECIMAL POINT
C-----OFFSET.
  110 X9IN=10.0D+00*X9IN+XDIG(K)
      IPT=IPT+1
  120 CONTINUE
C-----ENTIRE FIELD TRANSLATED (NO EXPONENT). CONVERT TO FLOATING POINT.
      GO TO 150
C-----BLANK FOUND (END OF MANTISSA). SCAN REMAINDER OF FIELD FOR
C-----EXPONENT.
  130 I=J+1
      IF(I.GT.11) GO TO 150
      DO 140 J=I,11
      IFIELD=FIELD(J)
      IF(IFIELD.EQ.' ') GO TO 140
      IF(IFIELD.EQ.'e'.OR.IFIELD.EQ.'E') GO TO 160
      IF(IFIELD.EQ.'d'.OR.IFIELD.EQ.'D') GO TO 160
      IF(IFIELD.EQ.'-') GO TO 190
      IF(IFIELD.EQ.'+') GO TO 170
C-----ERROR. CANNOT IDENTIFY CHARACTER.
      GO TO 240
  140 CONTINUE
C-----ENTIRE FIELD TRANSLATED (NO EXPONENT). CONVERT TO FLOATING POINT.
  150 E=X9IN
      IF(IPT.GT.0) E=E/TENS(IPT)
      E=SIGN*E
      RETURN
C
C     TRANSLATE EXPONENT.
C
C-----BEGINNING OF EXPONENT FOUND (E OR D). CHECK FOR FOLLOWING - OR +.
  160 J=J+1
      IFIELD=FIELD(J)
      IF(IFIELD.EQ.'-') GO TO 190
      IF(IFIELD.NE.'+') GO TO 180
C----- + FOUND. INITIALIZE EXPONENT SIGN.
  170 J=J+1
  180 ISIGN=1
      GO TO 200
C----- - FOUND. INITIALIZE EXPONENT SIGN.
  190 J=J+1
      ISIGN=-1
C-----INITIALIZE EXPONENT AND SCAN REMAINING CHARACTERS FOR EXPONENT.
  200 IEXP=0
      DO 230 I=J,11
      IFIELD=FIELD(I)
      IF(IFIELD.EQ.' ') GO TO 230
      DO 210 K=0,9
      IF(IFIELD.EQ.DIGIT(K)) GO TO 220
  210 CONTINUE
C-----ERROR. CANNOT IDENTIFY CHARACTER.
      GO TO 240
C-----DIGIT FOUND. INCREMENT EXPONENT.
C-----OFFSET.
  220 IEXP=10*IEXP+K
  230 CONTINUE
C-----ENTIRE FIELD TRANSLATED (WITH EXPONENT). CONVERT TO FLOATING
C-----POINT.
      E=X9IN
      IEXP=ISIGN*IEXP
      IF(IPT.GT.0) IEXP=IEXP-IPT
C-----12/18/2012 - EXTENDED FOR 3 DIGIT EXPONENTS
C-----06/12/2013 - BACK tO 2 digit exponent - 3 digits caused problems
      IF(iabs(IEXP).GT.99) GO TO 240
      E=SIGN*E*TENS(IEXP)
      RETURN
C
C     ERROR CONDITIONS.
C
C-----ILLEGAL CHARACTER.
  240 MESS(J)='*'
      write(*,250) FIELD,MESS
  250 FORMAT(1X,11A1/1X,11A1/' ERROR in Input Data...Translated as 0.')
      E=0.0D+00
      MESS(J)=' '
      RETURN
      END
      SUBROUTINE SCRATCH1(ISCR,FILE12)
C=======================================================================
C
C     PURPOSE
C     =======
C     OPEN SCRATCH FILE EITHER WITH OR WITHOUT FILENAME.
C
C     ARGUMENTS
C     =========
C     ISCR     = I/O NUMBER NUMBER (INTEGER)
C     FILENAME = FILE NAME (CHRACTER*12)
C
C     VERSIONS
C     ========
C     THERE ARE 2 VERSIONS OF THIS ROUTINE,
C     scratcha.f = OPEN WITH FILENAME
C     scratchb.f = OPEN WITHOUT FILENAME
C
C     YOU NEED MERELY LOAD THE CODES WITH WHICHEVER IS BEST
C     FOR YOU - SEE BELOW FOR DETAILS
C
C     THERE IS ACTUALLY ONLY ONE VERSION WITH PART OF THE
C     CODE COMMENTED OUT, TO ONLY OPEN FILES EITHER WITH OR
C     WITHOUT FILE NAMES - SEE BELOW.
C
C     COMPILER DEPENDENCE
C     ===================
C     1) MOST COMPILERS ALLOW SCRATCH FILES TO BE OPENED
C        WITH OR WITHOUT NAMES = RECOMMENDED OPEN WITH
C     2) ON IBM-PC ABSOFT AND MICROSOFT DO NOT ALLOW NAMES
C     3) ON IBM-PC LAHEY REQUIRES NAMES (GETS CONFUSED IF
C        THERE ARE MULTIPLE SCRATCH FILES WITHOUT NAMES)
C
C     THE CONVENIENCE OF NAMES
C     ========================
C     1) MANY SYSTEMS WILL CREATE SCRATCH FILES WITH RANDOM
C        NAMES, AND NOT DELETE THEM WHEN EXECUTION ENDS
C     2) THIS CAN LEAD TO AN ACCUMULATION OF MANY UNUSED
C        SCRATCH FILES
C     3) WITH NAMES, SOME SYSTEMS WILL DELETE AND SOME WON'T
C        WHEN EXECUTION ENDS - IS YOU END UP WITH EITHER
C        NO SCRATCH FILES OR A FIXED NUMBER THAT DOESN'T
C        ACCUMULATE OVER MANY RUNS
C     4) THAT'S WHY USING NAMES IS RECOMMENDED
C
C     WARNING
C     =======
C     IN ORDER TO WORK PROPERLY THE FILENAME MUST BE
C     EXACTLY 12 CHARACTER LONG = MAXIMUM ALLOWED FOR
C     IBM-PC FILENAME COMPATIBILITY
C     MAXIMUM = 8 CHARACTER FILE NAME
C               1 CHARACTER .
C               3 CHARACTER EXTENSION
C
C     IF THE ACTUAL NAME IS LESS THAN 12 CHARACTERS BE SURE TO
C     PAD THE NAME OUT TO 12 CHARACTERS WHEN THIS ROUTINE IS
C     CALLED, E.G., IF THE NAME IS - EXAMPLE.X, USE,
C     'EXAMPLE.X   '
C      123456789012
C
C     NOTE - THE 3 TRAILING BLANKS TO PAD THE NAME TO 12 CHARACTERS.
C
C=======================================================================
      INCLUDE 'implicit.h'
      CHARACTER*12 FILE12
C-----USE FILENAME TO PREVENT COMPILER WARNING
      IF(FILE12.EQ.'            ') RETURN
C
C     BELOW IS THE ONLY DIFFERENCE BETWEEN SCRATCHA AND SCRATCHB
C     SCRATCHA = OPEN WITH    FILENAME
C     SCRSTCHB = OPEN WITHOUT FILENAME
C
C***** DEBUG
C-----SCRATCHA: OPEN FILE WITH    FILENAME
C     OPEN(ISCR,FILE=FILE12,FORM='UNFORMATTED',STATUS='SCRATCH')
C***** DEBUG
C-----SCRATCHB: OPEN FILE WITHOUT FILENAME
      OPEN(ISCR,            FORM='UNFORMATTED',STATUS='SCRATCH')
C***** DEBUG
      RETURN
      END
      SUBROUTINE ENDIT
C=======================================================================
C
C     PRINT EXECUTION TIME AND TERMINATE = NORMAL FINISH
C
C=======================================================================
CAK   CALL TIMER
      CALL TIMERpre
      STOP
      ENTRY ENDERROR
C=======================================================================
C
C     ENTRY POINT TO STOP ON ERROR.
C
C=======================================================================
      CALL TIMEERR
      STOP
      END
CAK   SUBROUTINE TIMER
      SUBROUTINE TIMERpre
C=======================================================================
C
C     TOTAL EXECUTION TIME
C
C     WARNING - ALL TIMES ARE IN SINGLE PRECISION.
C               DO NOT ADD ANY DOUBLE OR REAL*8 STATEMENTS.
C
C=======================================================================
      IMPLICIT REAL*4 (A-H,O-Z)  ! note - REAL*4 not 8
      IMPLICIT INTEGER*4 (I-N)
      SAVE
      CHARACTER*8 CODENAME
      COMMON/NAMECODE/CODENAME
      INTEGER*4 OUTP,OTAPE
      COMMON/ENDFIO/INP,OUTP,ITAPE,OTAPE
      DATA TSTART/0.0/
      DATA IPASS/0/
C-----DEFINE CURRENT TIME
      CALL TIMEIT(TNOW)
C-----ON FIRST PASS DEFINE STARTING TIME
      IF(IPASS.EQ.0) TSTART=TNOW
      IPASS=IPASS+1
C-----PRINT EVERY PASS EXCEPT FIRST ONE
      IF(IPASS.LE.1) RETURN
      WRITE(OUTP,10) CODENAME,TNOW-TSTART
      WRITE(OUTP,20)
C-----OUTPUT TO SCREEN ONLY IF ENDF/B OUTPUT IS PRODUCED
      IF(OTAPE.GT.0) THEN
      WRITE(*   ,10) CODENAME,TNOW-TSTART
      WRITE(*   ,20)
      ENDIF
   10 FORMAT(1X,78('=')/1X,A8,' Total Execution Time',F20.2,' Seconds')
   20 FORMAT(1X,78('='))
      RETURN
C=======================================================================
C
C     ENTRY TIME TO RETURN ELAPSED TIME
C
C     WARNING - TIMER MUST BE CALLED FIRST TO DEFINE TSTART
C
C=======================================================================
      ENTRY TIMER1(SECONDS)
      CALL TIMEIT(TNOW)
      SECONDS = TNOW - TSTART
      RETURN
C=======================================================================
C
C     ENTRY TIME TO RETURN TIME FOR EACH MAT.
C     IDENTICAL TO TIMER, BUT MAT, NOT CODENAME.
C
C=======================================================================
      ENTRY TIMEMAT
      CALL TIMEIT(TNOW)
      WRITE(OUTP,30) TNOW-TSTART
      WRITE(OUTP,20)
      IF(OTAPE.GT.0) THEN
      WRITE(*   ,30) TNOW-TSTART
      WRITE(*   ,20)
      ENDIF
   30 FORMAT(1X,78('=')/1X,'     MAT',' Total Execution Time',F20.2,
     1 ' Seconds')
      RETURN
C=======================================================================
C
C     ENTRY TIME TO RETURN TIME FOR ERROR STOP.
C     IDENTICAL TO TIMER, BUT ERROR, NOT CODENAME.
C
C=======================================================================
      ENTRY TIMEERR
      CALL TIMEIT(TNOW)
      WRITE(OUTP,40) TNOW-TSTART
      WRITE(OUTP,20)
      IF(OTAPE.GT.0) THEN
      WRITE(*   ,40) TNOW-TSTART
      WRITE(*   ,20)
      ENDIF
   40 FORMAT(1X,78('=')/1X,'   ERROR',' Total Execution Time',F20.2,
     1 ' Seconds')
      RETURN
      END
      SUBROUTINE TIMEIT(SECONDS)
C=======================================================================
C
C     IBM-PC - TIME SINCE START OF PROBLEM (BASE TIME) IN SECONDS
C     UNIX     This works on ALL of these.
C     LINUX
C     MAC
C
C     WARNING - ALL TIMES ARE IN SINGLE PRECISION -
C               DO NOT ADD DOUBLE PRECISION OR REAL*8 STATEMENTS.
C
C=======================================================================
      IMPLICIT REAL*4 (A-H,O-Z) ! note REAL*4, not 8
      SAVE
      DIMENSION TARRAY(2)
      SECONDS=ETIME(TARRAY)
      SECONDS=(TARRAY(1)+TARRAY(2))
      RETURN
      END

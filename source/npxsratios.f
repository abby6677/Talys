      subroutine npxsratios
c
c +---------------------------------------------------------------------------
c | Author: Marilena Avrigeanu
c | Date  : April, 2020
c | Task  : Reads from 'structure' breakup n & p c.s. ratios requested
c           by the inelastic breakup enhancement calculation,
c           Eq. 2 from Phys. Rev. C 94,014606 (2016)
c +---------------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include      "talys.cmb"
      logical      lexist
      character*3  Zstring
      character*90 enhBUnxs,enhBUpxs
      integer      Zcomp,Acomp,Ncomp,Z,A,N,type,ien,Ioutn,iep,Ioutp
      integer      Nix,Zix
      real         ENHN,ENHP,enout,epout
      real*8       NUCLEUN,NUCLEUP,NUCLEUNr,NUCLEUPr,Eoutn,Eoutp
c
c ***  n,p reaction cross sections ratios for breakup enhancement  ***
c
c path          : directory containing structure files to be read
c lenpath       : length of pathname
c enhBUnxs      : breakup neutron enhancing ratios file
c enhBUpxs      : breakup proton enhancing ratios file
c ENHratioN/P   : cross sections ratios for breakup enhancement,
c                 Eq. (2), PRC 94,014606 (2016)                                  
c                 n + Atarget:  sig(n,Z,A,Enout)/sig_Total(n,Enout)
c                 p + Atarget:  sig(p,Z,A,Epout)/sig_Reaction(p,Epout)
c                 included in STRUCTURE directory as follows:
c                 structure\TENDL_npxs\neutrons\ENHratioN.da
c                 structure\TENDL_npxs\protons\ENHratioP.da
c numenout      : maximal number of breakup nucleons outgoing energies
c                 included in talys.cmb as well as in reacinitial
c Zix           : charge number index for residual nucleus
c Nix           : neutron number index for residual nucleus
c
      write(8,*)' '
      write(8,*) '  start npxsratios  for   d + ',Atarget,nuc(Ztarget)
      write(8,*)' '
      write(8,*)' numZ numN numA numenout',numZ,NumN,numZ+numN,numenout
c     
      do 115 Ien=1,numenout
        do 115 Nix=0,numN
          do 115 Zix=0,numZ
            do 115 type=0,numpar                         
                     ENHratio(type,Zix,Nix,Ien)=0.          
115      continue
c
c     Inquire whether breakup enhancing ratios file are presen
c
       Zstring='000'
       write(Zstring,'(i3.3)') Ztarget
       do type=1,2
        if (type.eq.1) then
c         enhBUnxs=trim(path)//'TENDL_npxs/neutrons/ENHratioN.dat'
          enhBUnxs=trim(path)//'TENDL_npxs/neutrons/ENHratioN.'//Zstring
        inquire (file=enhBUnxs,exist=lexist)
          if (.not.lexist) go to 409
        open (unit=25,status='unknown',file=enhBUnxs)
c         enhBUpxs=trim(path)//'TENDL_npxs/protons/ENHratioP.dat'
          enhBUpxs=trim(path)//'TENDL_npxs/protons/ENHratioP.'//Zstring
        inquire (file=enhBUpxs,exist=lexist)
        if (.not.lexist) go to 409
          open (unit=26,status='unknown',file=enhBUpxs)
        endif
      enddo     
c
      ebubin = 0.1
c          
c ************     Inelastic breakup enhancing ratios     ************
c
c ebubin       : outgoing breakup nucleon energy bin for integration
c                in Eq.(2), PRC 94,014606
c Atarget      : mass number of target nucleus
c Ztarget      : charge number of target nucleus
c type=1,2     : breakup neutron, proton from d + Atarget interaction
c Acomp        : mass number index for compound nucleus
c Zcomp        : charge number index for compound nucleus
c Ncomp        : neutron number index for compound nucleus     
c ZZ,Z         : charge number of residual nucleus
c AA,A         : mass number of residual nucleus
c N            : neutron number of residual nucleus
c NUCLEUN      : residual nuclei of interest for breakup enhancemen
c                brought by type=1 (breakup neutrons) are those in
c                the area:  Zcomp=1:5
c                           Ncomp=0:4
c                           Acomp=1:9
c NUCLEUP      : residual nuclei of interest for breakup enhancemen
c                brought by type=2 (breakup protons) are those in
c                the area:  Zcomp=0:4
c                           Ncomp=1:5
c                           Acomp=1:9
c Ioutn        : number of breakup neutron outgoing energies
c Ioutp        : number of breakup proton outgoing energies
c enout        : outgoing energy for breakup neutron
c epout        : outgoing energy for breakup proton
c NUCLEUNr     : help variable
c NUCLEUPr     : help variable
c Eoutn,Eoutp  : help variables
c
c ************   Read breakup neutrons enhancing ratios   ************
c
      write(8,*) 'residual nuclei for BU enhancement calculation'
      do 33 type=1,2
       if(type.eq.1)then
         do 15 Zcomp=1,5
         do 13 Acomp=1,9
           Ncomp=Acomp-Zcomp
           if (Ncomp.lt.0.or.Ncomp.gt.4) goto 11
           Z=ZZ(Zcomp,Ncomp,0)
           A=AA(Zcomp,Ncomp,0)       
           N=AA(Zcomp,Ncomp,0)-ZZ(Zcomp,Ncomp,0)
           NUCLEUN=Ztarget*1.0D+9+Atarget*1.0D+6+Z*1.0D+3+A
5            read(25,*,end=409)NUCLEUNr,Eoutn
           Ioutn=int(Eoutn)
           if(NUCLEUNr.EQ.NUCLEUN)then
             write(8,*)' NUCLEUN=',NUCLEUN,' Z=',Z,' A=',A,
     +                     ' N=',N,' * NUCLEUNr=',NUCLEUNr
             Zix=Zinit-Z
             Nix=Ninit-N
             do ien=1,Ioutn
               enout=ebubin*ien
               read(25,*)enout,ENHN
               ENHratio(type,Zix,Nix,ien)=ENHN
             enddo
           else
             do ien=1,Ioutn+1
               read(25,*)
             enddo
             go to 5
           endif
11           continue     
13         rewind (25)
         continue
15         continue
c          
       else                                 
         write(8,*)'               ***'
c
c ************   Read breakup protons enhancing ratios   ************
c
         do 29 Zcomp=0,4
         do 27 Acomp=1,9
           Ncomp=Acomp-Zcomp
           if (Ncomp.le.0.or.Ncomp.gt.5) goto 25
           Z=ZZ(Zcomp,Ncomp,0)
           A=AA(Zcomp,Ncomp,0)       
           N=AA(Zcomp,Ncomp,0)-ZZ(Zcomp,Ncomp,0)
           NUCLEUP=Ztarget*1.0D+9+Atarget*1.0D+6+Z*1.0D+3+A           
21            read(26,*,end=409)NUCLEUPr,Eoutp
           Ioutp=int(Eoutp)
             if(NUCLEUPr.EQ.NUCLEUP)then
             write(8,*)' NUCLEUP=',NUCLEUP,' Z=',Z,' A=',A,
     +                     ' N=',N,' * NUCLEUPr=',NUCLEUPr
              Zix=Zinit-Z
              Nix=Ninit-N
              do iep=1,Ioutp
                   epout=ebubin*Iep
               read(26,*)epout,ENHP
               ENHratio(type,Zix,Nix,Iep)=ENHP
             enddo
               else
             do ien=1,Ioutp+1
               read(26,*)
             enddo
            go to 21
           endif
25           continue     
27         rewind(26)
         continue
29         continue
       endif
33     continue
c
          write(8,*)'  '
          write(8,*)' *** Test npxsratios *** '
c
      do 319 type=1,2
      do 317 Zcomp=1,5
      do 315 Acomp=1,9
       Ncomp=Acomp-Zcomp                    
       if (Ncomp.lt.0.or.Ncomp.gt.4) goto 313
       Z=ZZ(Zcomp,Ncomp,0)
       A=AA(Zcomp,Ncomp,0)       
       N=AA(Zcomp,Ncomp,0)-ZZ(Zcomp,Ncomp,0)
       Zix=Zinit-Z
       Nix=Ninit-N
       write(8,*)' type Z A N ENHratio(type,Z,N,320)',
     +               type,Z,A,N,ENHratio(type,Zix,Nix,320)
313       continue
315     continue
317     continue
319     continue     
c
      write(8,*)' *** END test npxsratios *** '
      write(8,*)'  '
c
      go to 601
c
c               
409     write(8,*)' '
      write(8,*)'subroutine npxsratios: missing TENDL_npxs library in',
     +          ' TALYS structure or missing xs files from TENDL-2019 '
      write(8,*)'***     NO BUenhancement taken into account     ***'
      write(8,*)'  '
c
      write(*,*)' '
      write(*,*)'subroutine npxsratios: missing TENDL_npxs library in',
     +          ' TALYS structure or missing xs files from TENDL-2019 '
      write(*,*)'***     NO BUenhancement taken into account     ***'
      write(*,*)' '
c
c     
      do 599 Ien=1,numenout
        do 599 Nix=0,numN
          do 599 Zix=0,numZ
            do 599 type=0,numpar                         
                     ENHratio(type,Zix,Nix,Ien)=0.          
599      continue
c
c
601     CONTINUE
      return
      end

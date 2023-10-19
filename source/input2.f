      subroutine input2
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 11, 2021
c | Task  : Read input for second set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      logical      fcol,flagassign
      character*1  ch
      character*80 word(40),key,value,cval
      integer      type,Zix,Nix,col(0:numZ,0:numN),i,ip,i2,iz,ia,ldmod,
     +             nex,class,ival,ibar,irad,lval,igr
      real         val,sfthall,sfexpall,sfth,sfexp
c
c ************* Defaults for second set of input variables *************
c
c outtype     : type of outgoing particles
c maxZ,numZ   : maximal number of protons away from the initial
c               compound nucleus
c maxN,numN   : maximal number of neutrons away from the initial
c               compound nucleus
c nbins0      : number of continuum excitation energy bins
c segment     : number of segments to divide emission energy grid
c nlevmax     : maximum number of included discrete levels for target
c Ltarget     : excited level of target
c nlevmaxres  : maximum number of included discrete levels for residual
c               nucleus
c nlevbin     : number of excited levels for binary nucleus
c k0          : index of incident particle
c Lisoinp     : user assignment of target isomer number
c isomer      : definition of isomer in seconds
c core        : even-even core for weakcoupling (-1 or 1)
c gammax      : number of l-values for gamma multipolarity
c nangle      : number of angles
c numang      : maximum number of angles
c nanglecont  : number of angles for continuum
c maxenrec    : number of recoil energies
c massmodel   : model for theoretical nuclear mass
c disctable   : table with discrete levels
c flagmicro   : flag for completely microscopic Talys calculation
c ldmodelall  : level density model for all nuclides
c ldmodelCN   : level density model for compound nucleus
c flagcolall  : flag for collective enhancement of level density for all
c               nuclides
c fislim      : mass above which nuclide fissions
c wmode       : designator for width fluctuation model
c preeqmode   : designator for pre-equilibrium model
c mpreeqmode  : designator for multiple pre-equilibrium model
c phmodel     : particle-hole state density model
c nlev        : number of excited levels for nucleus
c ldmodel     : level density model
c skipCN      : flag to skip compound nucleus in evaporation chain
c col         : help variable
c flagcol     : flag for collective enhancement of level density
c numlev      : maximum number of included discrete levels
c flagomponly : flag to execute ONLY an optical model calculation
c flagequi    : flag to use equidistant excitation bins instead of
c               logarithmic bins
c flagequispec: flag to use equidistant bins for emission spectra
c flagpopMeV  : flag to use initial population per MeV instead of
c               histograms
c flagmassdis : flag for calculation of fission fragment mass yields
c flagracap   : flag for radiative capture model
c ldmodelracap: level density model for direct radiative capture
c spectfacexp : experimental spectroscopic factor
c spectfacth  : theoretical spectroscopic factor
c maxZrp      : maximal number of protons away from the initial
c               compound nucleus before new residual evaporation
c maxNrp      : maximal number of neutrons away from the initial
c               compound nucleus before new residual evaporation
c
      do 10 type=0,6
        outtype(type)=' '
   10 continue
      maxZ=numZ-2
      maxN=numN-2
      nbins0=40
      if (flagffruns.or.flagrpruns) then
        segment=2
        flagequispec=.true.
      else
        segment=1
        flagequispec=.false.
      endif
      nlevmax=max(30,Ltarget)
      nlevmaxres=10
      do 20 type=0,6
        if (type.le.2.or.type.eq.6) then
          nlevbin(type)=10
        else
          nlevbin(type)=5
        endif
   20 continue
      nlevbin(k0)=nlevmax
      Lisoinp=-1
      isomer=1.
      core=-1
      gammax=2
      nangle=numang
      nanglecont=18
      maxenrec=numenrec
      massmodel=2
      disctable=1
      phmodel=1
      if (flagmicro) then
        ldmodelall=5
      else
        ldmodelall=1
      endif
      ldmodelCN=0
      if (Atarget.gt.fislim) then
        flagcolall=.true.
      else
        flagcolall=.false.
      endif
      preeqmode=2
      wmode=1
      mpreeqmode=2
      sfthall=1.
      sfexpall=0.347
      if (mod(Atarget,2).ne.0) sfexpall=1.
      do 30 Nix=0,numN
        do 30 Zix=0,numZ
          nlev(Zix,Nix)=0
          ldmodel(Zix,Nix)=0
          skipCN(Zix,Nix)=0
          col(Zix,Nix)=0
          flagcol(Zix,Nix)=flagcolall
          spectfacth(Zix,Nix)=0.
          do 30 nex=0,numlev
            spectfacexp(Zix,Nix,nex)=0.
   30 continue
      flagomponly=.false.
      flagequi=.true.
      flagpopMeV=.false.
      flagmassdis=.false.
      flagracap=.false.
      ldmodelracap=1
      maxZrp=numZ-2
      maxNrp=numN-2
c
c **************** Read second set of input variables ******************
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input
c              line
c inline     : input line
c word       : words on input line
c key        : keyword
c value      : value or string
c ch         : character
c
c The keyword is identified and the corresponding values are read.
c Erroneous input is immediately checked. The keywords and number of
c values on each line are retrieved from the input.
c
      do 110 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c Test for keywords
c
c parsym   : symbol of particle
c Zinit    : charge number of initial compound nucleus
c Ninit    : neutron number of initial compound nucleus
c getvalues: subroutine to assign values to keywords
c fcol     : flag for collective enhancement
c sfexp    : variable for spectrocopic factor
c sfexpall : variable for spectrocopic factor
c sfth     : variable for spectrocopic factor
c sfthall  : variable for spectrocopic factor
c
        if (key.eq.'ejectiles') then
          ip=-1
          do 210 i2=2,40
            ch=word(i2)(1:1)
            do 220 type=0,6
              if (ch.eq.parsym(type)) then
                ip=ip+1
                if (ip.le.6) outtype(ip)=ch
                goto 210
              endif
  220       continue
            if (ip.eq.-1) goto 300
  210     continue
          goto 110
        endif
        if (key.eq.'maxz') then
          read(value,*,end=300,err=300) maxZ
          goto 110
        endif
        if (key.eq.'maxn') then
          read(value,*,end=300,err=300) maxN
          goto 110
        endif
        if (key.eq.'bins') then
          read(value,*,end=300,err=300) nbins0
          goto 110
        endif
        if (key.eq.'segment') then
          read(value,*,end=300,err=300) segment
          goto 110
        endif
        if (key.eq.'maxlevelstar') then
          read(value,*,end=300,err=300) nlevmax
          nlevmax=max(nlevmax,Ltarget)
          nlevbin(k0)=nlevmax
          goto 110
        endif
        if (key.eq.'maxlevelsres') then
          read(value,*,end=300,err=300) nlevmaxres
          goto 110
        endif
        if (key.eq.'liso') then
          read(value,*,end=300,err=300) Lisoinp
          goto 110
        endif
        if (key.eq.'isomer') then
          read(value,*,end=300,err=300) isomer
          goto 110
        endif
        if (key.eq.'core') then
          read(value,*,end=300,err=300) core
          goto 110
        endif
        if (key.eq.'gammax') then
          read(value,*,end=300,err=300) gammax
          goto 110
        endif
        if (key.eq.'angles') then
          read(value,*,end=300,err=300) nangle
          goto 110
        endif
        if (key.eq.'anglescont') then
          read(value,*,end=300,err=300) nanglecont
          goto 110
        endif
        if (key.eq.'maxenrec') then
          read(value,*,end=300,err=300) maxenrec
          goto 110
        endif
        if (key.eq.'massmodel') then
          read(value,*,end=300,err=300) massmodel
          goto 110
        endif
        if (key.eq.'disctable') then
          read(value,*,end=300,err=300) disctable
          goto 110
        endif
        if (key.eq.'ldmodel') then
          read(value,*,end=300,err=300) ldmod
          read(word(3),*,end=230,err=300) iz
          read(word(4),*,end=300,err=300) ia
  230     if (word(3).eq.' ') then
            ldmodelall=ldmod
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
              goto 1000
            else
              ldmodel(Zix,Nix)=ldmod
            endif
          endif
          goto 110
        endif
        if (key.eq.'ldmodelcn') then
          read(value,*,end=300,err=300) ldmodelCN
          goto 110
        endif
        if (key.eq.'colenhance') then
          if (ch.eq.'n') fcol=.false.
          if (ch.eq.'y') fcol=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          read(word(3),*,end=240,err=300) iz
          read(word(4),*,end=300,err=300) ia
  240     if (word(3).eq.' ') then
            flagcolall=fcol
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
              goto 1000
            else
              flagcol(Zix,Nix)=fcol
              col(Zix,Nix)=1
            endif
          endif
          goto 110
        endif
        if (key.eq.'skipcn') then
          read(word(2),*,end=300,err=300) iz
          read(word(3),*,end=300,err=300) ia
          Zix=Zinit-iz
          Nix=Ninit-ia+iz
          if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
            goto 1000
          else
            skipCN(Zix,Nix)=1
          endif
          goto 110
        endif
        if (key.eq.'widthmode') then
          read(value,*,end=300,err=300) wmode
          goto 110
        endif
        if (key.eq.'preeqmode') then
          read(value,*,end=300,err=300) preeqmode
          goto 110
        endif
        if (key.eq.'mpreeqmode') then
          read(value,*,end=300,err=300) mpreeqmode
          goto 110
        endif
        if (key.eq.'phmodel') then
          read(value,*,end=300,err=300) phmodel
          goto 110
        endif
        if (key.eq.'nlevels') then
          class=2
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) nlev(Zix,Nix)=ival
          goto 110
        endif
        if (key.eq.'maxlevelsbin') then
          class=7
          call getvalues(class,word,Zix,Nix,type,
     +      ibar,irad,lval,igr,val,ival,cval,flagassign)
          if (flagassign) nlevbin(type)=ival
          goto 110
        endif
        if (key.eq.'omponly') then
          if (ch.eq.'n') flagomponly=.false.
          if (ch.eq.'y') flagomponly=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'equidistant') then
          if (ch.eq.'n') flagequi=.false.
          if (ch.eq.'y') flagequi=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'equispec') then
          if (ch.eq.'n') flagequispec=.false.
          if (ch.eq.'y') flagequispec=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'popmev') then
          if (ch.eq.'n') flagpopMeV=.false.
          if (ch.eq.'y') flagpopMeV=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'massdis') then
          if (ch.eq.'n') flagmassdis=.false.
          if (ch.eq.'y') flagmassdis=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'racap') then
          if (ch.eq.'n') flagracap=.false.
          if (ch.eq.'y') flagracap=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 300
          goto 110
        endif
        if (key.eq.'ldmodelracap') then
          read(value,*,end=300,err=300) ldmodelracap
          goto 110
        endif
        if (key.eq.'maxzrp') then
          read(value,*,end=300,err=300) maxZrp
          goto 110
        endif
        if (key.eq.'maxnrp') then
          read(value,*,end=300,err=300) maxNrp
          goto 110
        endif
        if (key.eq.'sfth') then
          read(value,*,end=300,err=300) sfth
          read(word(3),*,end=250,err=300) iz
          read(word(4),*,end=300,err=300) ia
  250     if (word(3).eq.' ') then
            sfthall=sfth
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN) then
              goto 1000
            else
              spectfacth(Zix,Nix)=sfth
            endif
          endif
          goto 110
        endif
        if (key.eq.'sfexp') then
          nex=-1
          Zix=0
          Nix=0
          read(value,*,end=300,err=300) sfexp
          read(word(3),*,end=260,err=300) iz
          read(word(4),*,end=260,err=300) ia
          read(word(5),*,end=260,err=300) nex
  260     if (word(3).eq.' ') then
            sfexpall=sfexp
          else
            Zix=Zinit-iz
            Nix=Ninit-ia+iz
            if (Zix.lt.0.or.Zix.gt.numZ.or.Nix.lt.0.or.Nix.gt.numN)
     +        goto 1000
            if (word(3).eq.' '.or.word(4).eq.' ') nex=iz
            if (nex.eq.-1) then
              sfexpall=sfexp
            else
              if (nex.lt.0.or.nex.gt.numlev) then
                write(*,'(" TALYS-error: 0 <= nex <= numlev ",a80)')
     +            inline(i)
                stop
              else
                spectfacexp(Zix,Nix,nex)=sfexp
              endif
            endif
          endif
          goto 110
        endif
        goto 110
 1000   write(*,'(" TALYS-warning: Z,N index out of range,",
     +    " keyword ignored: ",a80)') inline(i)
  110 continue
c
c Set level density models and spectroscopic factors per nucleus
c
      if (ldmodelCN.gt.0) then
        ldmodel(0,0)=ldmodelCN
      else
        ldmodelCN=ldmodelall
      endif
      do 310 Nix=0,numN
        do 310 Zix=0,numZ
          if (ldmodel(Zix,Nix).eq.0) ldmodel(Zix,Nix)=ldmodelall
          if (col(Zix,Nix).eq.0) flagcol(Zix,Nix)=flagcolall
          if (spectfacth(Zix,Nix).eq.0.) spectfacth(Zix,Nix)=sfthall
          do 310 nex=0,numlev
            if (spectfacexp(Zix,Nix,nex).eq.0.)
     +        spectfacexp(Zix,Nix,nex)=sfexpall
  310 continue
      return
  300 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

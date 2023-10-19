      subroutine checkkeyword
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 29, 2021
c | Task  : Check for errors in keywords
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer      numkey
      parameter    (numkey=388)
      integer      i,j
      character*80 keyword(numkey),word(40),key
c
c Although it is difficult to prevent the user from all possible input
c errors, we can check for the use of wrong keywords and for unphysical
c values for most of the input variables.
c
c *********************** Check for wrong keywords *********************
c
c keyword: keyword
c numkey : number of keywords
c
c TALYS will stop if a keyword is incorrect
c
      data (keyword(i),i=1,numkey) /
     +  ' ', 'a', 'aadjust', 'abundance', 'adddiscrete', 'addelastic',
     +  'adepthcor', 'alimit', 'alphald', 'alphaomp', 'angles',
     +  'anglescont', 'anglesrec', 'aradialcor', 'area', 'astro',
     +  'astroe', 'astroex', 'astrogs', 'astrot', 'asys', 
     +  'autorot', 'avadjust', 'avadjustf', 'avdadjust', 'avdadjustf',
     +  'avsoadjust', 'avsoadjustf', 'awadjust', 'awadjustf',
     +  'awdadjust',
     +  'awdadjustf', 'awsoadjust', 'awsoadjustf', 'axtype', 'bdamp', 
     +  'bdampadjust', 'best', 'bestbranch', 'bestend', 'bestpath', 
     +  'beta2', 'betafiscor', 'betafiscoradjust', 'betald', 'bins', 
     +  'block', 'branch', 
     +  'breakupmodel', 'cbreak', 'cfermi', 'cfermibf', 'cglobal', 
     +  'channelenergy', 'channels', 'cknock', 'class2', 'class2file', 
     +  'class2width', 'cnubar1', 'cnubar2', 'colenhance', 'colldamp',
     +  'components', 'compound',
     +  'core', 'coulomb', 'cpang', 'cstrip', 'ctable', 
     +  'ctableadjust', 'ctmglobal', 'd0', 'd1adjust', 'd2adjust', 
     +  'd3adjust', 'ddxmode', 'deformfile', 'deltaw', 'densfile', 
     +  'deuteronomp', 
     +  'disctable', 'dispersion', 'e0', 'e0adjust', 'e1file', 
     +  'eciscalc', 'eciscompound', 'ecisdwba', 'ecissave', 'eback', 
     +  'ebeam', 'egr', 'egradjust', 'ejectiles', 'ejoin', 
     +  'electronconv', 'element', 'elow', 'elwidth', 'emsdmin', 
     +  'endf', 'endfdetail', 'endfecis', 'energy', 'epr', 
     +  'epradjust', 'equidistant', 'equispec', 'estop', 'esurf', 
     +  'etable', 'etableadjust', 'exmatch', 'exmatchadjust', 'expmass',
     +  'ffevaporation', 'ffmodel',
     +  'ffspin', 'fileangle', 'filechannels', 'fileddxa', 'fileddxe', 
     +  'filedensity', 'filediscrete', 'fileelastic', 'filefission', 
     +  'filegamdis', 
     +  'filepsf', 'filerecoil', 'fileresidual', 'filespectrum', 
     +  'filetotal', 
     +  'fisbar', 'fisbaradjust', 'fisfeed', 'fishw', 'fishwadjust', 
     +  'fismodel', 'fismodelalt', 'fispartdamp', 'fiso', 'fisom',
     +  'fission', 'fsadjust',
     +  'ftable', 'ftableadjust', 'fullhf', 'fymodel', 'g', 'gadjust',
     +  'gamgam', 'gamgamadjust', 'gammald', 'gammashell1', 
     +  'gammashell2', 
     +  'gammax', 'gefran', 'ggr', 'ggradjust', 'giantresonance', 
     +  'gn', 'gnadjust', 'gnorm', 'gp', 'gpadjust',
     +  'gpr', 'gpradjust', 'group', 'gshell', 'hbstate', 
     +  'hbtransfile', 'ibeam','incadjust', 'inccalc', 'integral', 
     +  'isomer', 'jlmmode', 'jlmomp', 'kph', 'krotconstant', 
     +  'kvibmodel', 'labddx', 'ldmodel', 'ldmodelcn', 'ldmodelracap',
     +  'levelfile', 'liso', 'localomp', 'ltarget', 'lurr', 'lv1adjust', 
     +  'lvadjust', 'lvsoadjust', 'lw1adjust', 'lwadjust', 'lwsoadjust',
     +  'm1file', 'm2constant', 'm2limit', 'm2shift', 'mass', 
     +  'massdir', 'massdis', 'massexcess', 'massmodel', 'massnucleus', 
     +  'maxband', 'maxchannel', 'maxenrec', 'maxlevelsbin', 
     _  'maxlevelsres',
     +  'maxlevelstar', 'maxn', 'maxnrp', 'maxrot', 'maxz', 
     +  'maxzrp', 'micro', 'mpreeqmode', 'msdbins', 'multipreeq', 
     +  'nafit', 'ngfit', 'nlevels', 'nlow', 'nnfit', 'nonthermlev', 
     +  'ntop', 'nulldev', 
     +  'ompenergyfile', 'omponly', 'onestep', 'optmod', 'optmodall', 
     +  'optmodfilen', 'optmodfilep', 'outangle', 'outbasic', 
     +  'outbinspectra', 
     +  'outcheck','outdecay', 'outdensity', 'outdirect', 'outdiscrete', 
     +  'outdwba', 'outecis', 'outexcitation', 'outfission', 'outfy', 
     +  'outgamdis', 'outgamma', 'outinverse', 'outkd', 'outlegendre', 
     +  'outlevels', 
     +  'outmain', 'outomp', 'outpopulation', 'outpreequilibrium',
     +  'outspectra', 
     +  'outtransenergy', 'pair', 'pairconstant', 'pairmodel', 'parity',
     +  'partable', 'pfnsmodel', 'pglobal', 'phmodel', 'popeps', 
     +  'popmev', 'preeqcomplex', 'preeqmode', 'preeqspin', 
     +  'preeqsurface', 'preequilibrium', 
     +  'production', 'projectile', 'pshift', 'pshiftadjust', 
     +  'pshiftconstant', 
     +  'psfglobal', 'ptable', 'ptableadjust', 'racap', 'radialfile', 
     +  'radialmodel', 'radiounit', 'rcadjust', 'rclass2mom', 
     +  'reaction',
     +  'recoil', 'recoilaverage', 'relativistic', 'rescuefile', 
     +  'reslib',
     +  'resonance', 'rfiseps', 'rgamma', 'rho', 'riplomp', 'riplrisk', 
     +  'rnunu', 'rnupi', 'rotational', 'rpevap', 'rpinu', 'rpipi', 
     +  'rprime', 'rspincut', 'rspincutff', 'rtransmom', 'rvadjust', 
     +  'rvadjustf', 'rvdadjust', 'rvdadjustf', 'rvsoadjust', 
     +  'rvsoadjustf', 'rwadjust', 
     +  'rwadjustf', 'rwdadjust', 'rwdadjustf', 'rwsoadjust', 
     +  'rwsoadjustf',
     +  's2adjust', 'sacs', 'segment', 'sfexp', 'sfth', 
     +  'sgr', 'sgradjust', 'shellmodel', 'skipcn', 'soswitch', 
     +  'soukho',
     +  'spherical', 'spincutmodel', 'spr', 'spradjust', 'statepot', 
     +  'strength', 'strengthm1', 'strucpath', 'sysreaction', 't', 
     +  'tadjust', 'tcool', 'tirrad', 'tjadjust', 'tmadjust', 
     +  'transeps', 'transpower', 'tres', 'twocomponent', 'ufermi', 
     +  'ufermibf', 'upbend', 'upbendc', 'upbende', 'upbendf',  'urr',
     +  'urrnjoy', 'v1adjust', 'v2adjust', 'v3adjust', 'v4adjust', 
     +  'vfiscor', 'vfiscoradjust', 'vso1adjust', 'vso2adjust', 
     +  'vinfadjust', 'w1adjust', 
     +  'w2adjust', 'w3adjust', 'w4adjust', 'wso1adjust', 'wso2adjust', 
     +  'widthfluc', 'widthmode', 'wtable', 'wtableadjust', 
     +  'xsalphatherm', 'xscaptherm', 'xseps', 'xsptherm', 
     +  'yieldfile', 'yieldunit'/
c
c A keyword can be de-activated by putting a # in front of it.
c All first words of the input lines are checked against the list
c of keywords.
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input
c              line
c inline     : input line
c word       : words on input line
c key        : keyword
c
c The keyword is identified.
c
      do 10 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        if (key(1:1).eq.'#') goto 10
        do 20 j=1,numkey
          if (keyword(j).eq.key) goto 10
   20   continue
        write(*,'(/" TALYS-error: Wrong keyword: ",a20)') key
        stop
   10 continue
      return
      end
Copyright (C)  2019 A.J. Koning, S. Hilaire and S. Goriely

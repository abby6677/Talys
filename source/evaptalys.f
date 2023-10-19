      subroutine evaptalys
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : November 1, 2015
c | Task  : Execute TALYS subroutines for secondary evaporation
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
c
c These subroutines are executed for each nuclide as a fission fragment
c or as a remaining resdiual product of a high-energy calculation.
c
c talysinput  : subroutine for user input and defaults
c talysinitial: subroutine for initialization of nuclear structure
c
      call talysinput
      call talysinitial
c
c basicxs     : subroutine for basic cross sections and transmission
c               coefficients
c parinclude  : logical to include outgoing particle
c gamma       : subroutine for gamma cross section and transmission
c               coefficients
c enincmax    : maximum incident energy
c epreeq      : on-set incident energy for preequilibrium calculation
c preeqinit   : subroutine for initialization of general
c               pre-equilibrium parameters
c excitoninit : subroutine for initialization of exciton model
c               parameters
c flagcomp    : flag for compound nucleus calculation
c compoundinit: subroutine for initialization of compound model
c               parameters
c
      call basicxs(0,0)
      if (parinclude(0)) call gamma(0,0)
      if (enincmax.ge.epreeq) then
        call preeqinit
        call excitoninit
      endif
      if (flagcomp) call compoundinit
c
c energies   : subroutine for energies
c reacinitial: subroutine for initialization of arrays for various
c              cross sections
c
      call energies
      call reacinitial
c
c Optical model
c
c incident: subroutine for main settings and basic cross sections for
c           incident energy
c exgrid  : subroutine to set excitation energy grid
c
      call incident
      call exgrid(0,0)
c
c Pre-equilibrium reactions
c
c flagpreeq : flag for pre-equilibrium calculation
c preeq     : subroutine for preequilibrium reactions
c population: subroutine for processing of pre-equilibrium spectra
c             into population bins
c
      if (flagpreeq) then
        call preeq
        call population
      endif
c
c Multiple emission
c
c multiple: subroutine for multiple emission
c
      call multiple
c
c Exclusive channels
c
c flagchannels: flag for exclusive channels calculation
c channels    : subroutine for exclusive reaction channels
c
      if (flagchannels) call channels
c
c Collecting total cross sections, spectra, angular distributions, etc.
c
c totalxs : subroutine for total cross sections
c flagspec: flag for output of spectra
c spectra : subroutine for creation of particle spectra
c residual: subroutine for residual production cross sections
c output  : subroutine for output
c
      call totalxs
      if (flagspec.or.flagddx) call spectra
      call residual
      call output
      return
      end
Copyright (C)  2015 A.J. Koning, S. Hilaire and S. Goriely

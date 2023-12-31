      block data gamdata
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : August 10, 2015
c | Task  : Kopecky's spline fit of radiative widths
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer indx
c
c ***************** Gamma_gamma tabulated values ***********************
c
c Spline fit to experimental data made by J. Kopecky, 2002.
c
c gamkopecky: radiative width in eV by spline fit of Kopecky
c
      data (gamkopecky(indx),indx=1,250) /
     +  1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,
     +  1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,
     +  1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,
     +  1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,1.15,
     +  1.1,1.08,1.02,1.005,0.95,1.0,1.0,1.15,1.3,1.5,
     +  1.7,1.95,2.2,2.4,2.6,2.8,2.8,2.65,2.3,1.9,
     +  1.6,1.3,1.0,0.8,0.65,0.5,0.46,0.42,0.4,0.36,
     +  0.33,0.32,0.315,0.315,0.305,0.295,0.285,0.275,0.265,0.25,
     +  0.24,0.23,0.22,0.21,0.20,0.19,0.185,0.175,0.17,0.165,
     +  0.16,0.155,0.15,0.145,0.14,0.135,0.132,0.128,0.125,0.12,
     +  0.1195,0.119,0.1185,0.118,0.1175,0.117,0.1165,0.116,0.114,0.112,
     +  0.11,0.109,0.1085,0.108,0.1075,0.107,0.106,0.106,0.106,0.105,
     +  0.105,0.104,0.103,0.103,0.102,0.101,0.1,0.097,0.094,0.091,
     +  0.088,0.085,0.082,0.079,0.076,0.073,0.072,0.071,0.07,0.068,
     +  0.066,0.065,0.063,0.061,0.06,0.062,0.064,0.066,0.068,0.072,
     +  0.075,0.078,0.081,0.083,0.087,0.091,0.093,0.096,0.098,0.101,
     +  0.103,0.106,0.11,0.107,0.104,0.101,0.095,0.092,0.089,0.084,
     +  0.081,0.078,0.075,0.073,0.071,0.068,0.065,0.064,0.062,0.06,
     +  0.0615,0.063,0.0645,0.066,0.07,0.074,0.078,0.088,0.0815,0.082,
     +  0.083,0.0845,0.086,0.095,0.11,0.128,0.15,0.18,0.21,0.32,
     +  0.5,0.75,1.05,1.38,1.5,0.75,0.3,0.12,0.06,0.05,
     +  0.045,0.042,0.04,0.039,0.038,0.037,0.036,0.035,0.034,0.033,
     +  0.033,0.0325,0.0325,0.032,0.0315,0.031,0.0305,0.03,0.03,0.03,
     +  0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,
     +  0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03/
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely

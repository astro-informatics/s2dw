
                            S2DW CODE
         Steerable scale discretised wavelets on the sphere
  ----------------------------------------------------------------

DESCRIPTION
     The S2DW code provides functionality to perform the scale
     discretised wavelet transform on the sphere developed in our
     paper: Exact reconstruction with directional wavelets on the
     sphere. Routines are provided to compute wavelet and scaling
     coefficients from the spherical harmonic coefficients of a signal
     on the sphere and to synthesise the spherical harmonic
     coefficients of the original signal from its wavelet and scaling
     coefficients. The reconstruction of the spherical harmonic
     coefficients of the original signal is exact to numerical
     precision. Typically, maximum reconstruction errors are of the
     order 10^(-12) or smaller. Please see our paper for further
     details of the wavelet transform and a discussion of typical
     reconstruction errors and execution times of this implementation.

VERSION
     Release 1.1

AUTHORS
     J. D. McEwen (http://www.jasonmcewen.org) and Y. Wiaux

REFERENCE
     Y. Wiaux, J. D. McEwen, P. Vandergheynst, and O. Blanc
     Exact reconstruction with directional wavelets on the
     sphere. Mon. Not. Roy. Astron. Soc., 388(2):770-788, 2008
     (arXiv:0712.3519)

DOCUMENTATION
     See doc/index.html

REQUIREMENTS
     Library: 
       FFTW (http://www.fftw.org/)
       CFITSIO (http://heasarc.gsfc.nasa.gov/docs/software/fitsio/)
     Utility programs:
       S2 (http://www.jasonmcewen.org/codes.html)
       HEALPix (http://healpix.jpl.nasa.gov/)

INSTALLATION
     See doc/index.html

DOWNLOAD
     http://www.jasonmcewen.org/codes.html

SUPPORT
     See http://cosmocoffee.info/ forums

NOTE
     The package is still under development
     Please report problems/bugs by email to Jason McEwen

LICENSE
     S2DW package to compute the scale discretised wavelet transform on
     the sphere
     Copyright (C) 2008  Jason McEwen and Yves Wiaux

     This program is free software; you can redistribute it and/or
     modify it under the terms of the GNU General Public License
     as published by the Free Software Foundation; either version 2
     of the License, or (at your option) any later version.

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details (LICENSE.txt).

     You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, 
     MA  02110-1301, USA.


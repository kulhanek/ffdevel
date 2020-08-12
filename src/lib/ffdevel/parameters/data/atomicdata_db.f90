! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
!
! FFDevel is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! FFDevel is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FFDevel. If not, see <http://www.gnu.org/licenses/>.
! ==============================================================================

module ffdev_atomicdata_db

use ffdev_constants
use ffdev_atomicdata_dat

! ==============================================================================
! source: Consistent van der Waals Radii for the Whole Main Group, all in A
! DOI: 10.1021/jp8111556

integer,parameter   :: VDW_RADII_MAXZ = 18

real(DEVDP)     :: vdw_radii(1:VDW_RADII_MAXZ) = (/ &
 1.10, 1.40,  &
 1.81, 1.53, 1.92, 1.70, 1.55, 1.52, 1.47, 1.54, &
 2.27, 1.73, 1.84, 2.10, 1.80, 1.80, 1.75, 1.88  &
/)

! ==============================================================================
! source: https://en.wikipedia.org/wiki/Effective_nuclear_charge
! hydrogen - decreased from 1.000 -> 0.9000

integer,parameter   :: ZEFF_CLEMENTI_MAXZ = 18

real(DEVDP)     :: zeff_clementi(1:ZEFF_CLEMENTI_MAXZ) = (/ &
 0.900, 1.688, &
 1.279, 1.912, 2.421, 3.136, 3.834, 4.453, 5.100, 5.758, &
 2.507, 3.308, 4.066, 4.285, 4.886, 5.482, 6.116, 6.764  &
/)

! ==============================================================================
! ionization potential of elements, all in eV
! source: https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
! ref: CRC

integer,parameter   :: IPEA_MAXZ = 18

! ionization potenatial
real(DEVDP)     :: atomicdata_ip0(1:IPEA_MAXZ) = (/ &
 13.59844, 24.58741,  &
  5.39172,  9.32270,  8.29803, 11.26030, 14.53414, 13.61806, 17.42282, 21.56460, &
  5.13908,  7.64624,  5.98577,  8.15169, 10.48669, 10.36001, 12.96764, 15.75962  &
/)

! ionization potential +
real(DEVDP)     :: atomicdata_ipp(1:IPEA_MAXZ) = (/ &
   0.0000, 54.41778, &
 75.64018, 18.21116, 25.15484, 24.38332, 29.6013, 35.11730, 34.97082, 40.96328, &
 47.28640, 15.03528, 18.82856, 16.34585, 19.7694, 23.33790, 23.81400, 27.62967  &
/)

! ==============================================================================

! linearized atomic densities of spherically symmetrized atoms

integer,parameter   :: ATDENS_MAX_Z_ALL = 86
integer             :: ATDENS_MAX_Z

real(DEVDP)     :: atomicdata_rho_b012(1:ATDENS_MAX_Z_ALL,3)

! like-only atoms --------------------------------------------------------------
! calculated by Horton

integer,parameter   :: atomicdata_rho_PBE0_def2QZVPP_maxZ = 18

real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_b012(1:atomicdata_rho_PBE0_def2QZVPP_maxZ * 3) = (/ &
    3.700590,     1.142950,     0.000000, &
    5.366070,     1.799430,     0.000000, &
    1.388400,     0.449738,     0.000000, &
    2.344620,     0.617900,     0.044300, &
    3.270980,     0.940463,     0.058098, &
    3.795113,     1.014338,     0.084911, &
    4.339434,     0.994576,     0.063807, &
    4.732243,     1.048393,     0.083294, &
    5.104189,     1.064493,     0.119258, &
    5.426750,     1.142565,     0.105795, &
    1.544900,     0.513380,     0.000000, &
    2.001590,     0.517720,     0.077590, &
    2.483738,     0.584371,     0.071295, &
    2.943062,     0.693044,     0.043570, &
    3.326283,     0.643320,     0.032266, &
    3.706425,     0.684930,     0.033126, &
    4.046629,     0.710852,     0.042020, &
    4.376460,     0.746810,     0.052440  &
/)

! ==============================================================================

end module ffdev_atomicdata_db

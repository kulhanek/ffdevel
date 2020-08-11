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
! source: https://en.wikipedia.org/wiki/Electron_affinity_(data_page)
! ref: CRC

integer,parameter   :: IPEA_MAXZ = 18

! electron affinities - ! FIXME - what to do with negative EA
real(DEVDP)     :: atomicdata_ipm(1:IPEA_MAXZ) = (/ &
 0.75420, -0.50000, &
 0.61805, -0.50000, 0.27972, 1.26212, -0.07000, 1.46111, 3.40119, -1.20000, &
 0.54793, -0.40000, 0.43283, 1.38952,  0.74661, 2.07710, 3.61273, -1.00000 &
/)

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

real(DEVDP)     :: atomicdata_rho_bm(1:ATDENS_MAX_Z_ALL)
real(DEVDP)     :: atomicdata_rho_b0(1:ATDENS_MAX_Z_ALL)
real(DEVDP)     :: atomicdata_rho_bp(1:ATDENS_MAX_Z_ALL)

real(DEVDP)     :: atomicdata_rho_am(1:ATDENS_MAX_Z_ALL)
real(DEVDP)     :: atomicdata_rho_a0(1:ATDENS_MAX_Z_ALL)
real(DEVDP)     :: atomicdata_rho_ap(1:ATDENS_MAX_Z_ALL)

! like-only atoms --------------------------------------------------------------
! calculated by Horton

integer         :: atomicdata_rho_PBE0_def2QZVPP_maxZ = 18
real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_bm(1:18) = (/ &
  2.5576,   2.6798, &
  0.9387,   1.7710,   2.3017,   2.8901,   3.3317,   3.7041,   4.1621,   2.9105, &
  1.0315,   1.5615,   1.9226,   2.3014,   2.6817,   3.0216,   3.3783,   2.4285  &
/)
real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_am(1:18) = (/ &
 -1.6272,  -1.4519, &
 -4.3363,  -2.2427,  -1.2817,  -0.3971,   0.1482,   0.4974,   0.9154,  -0.6041, &
 -4.1687,  -2.5828,  -1.6265,  -0.8177,  -0.1558,   0.3373,   0.7945,  -0.7770  &
/)
real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_b0(1:18) = (/ &
  3.7006,   5.3661, &
  1.3884,   2.3446,   3.2275,   3.7990,   4.3814,   4.6959,   5.0948,   5.4267, &
  1.5449,   2.0016,   2.4213,   2.9353,   3.3520,   3.6896,   4.0451,   4.3765  &
/)
real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_a0(1:18) = (/ &
 -1.2153,   0.3398, &
 -3.8653,  -1.6793,  -0.2552,   0.4146,   0.9908,   1.2061,   1.4900,   1.6313, &
 -3.4829,  -2.0237,  -1.0388,  -0.0482,   0.5722,   1.0230,   1.4331,   1.7463  &
/)
real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_bp(1:18) = (/ &
  0.0000,   7.1655, &
  8.0572,   3.0068,   4.4144,   4.8648,   5.4188,   5.9755,   6.2974,   6.6751, &
  7.3439,   2.5969,   3.2707,   3.6823,   4.0010,   4.4799,   4.8010,   5.1757  &
/)
real(DEVDP)     :: atomicdata_rho_PBE0_def2QZVPP_ap(1:18) = (/ &
  0.0000,   0.6189, &
  0.5422,  -1.4942,   0.8243,   1.2137,   1.6719,   2.0694,   2.2212,   2.3825, &
  3.1808,  -1.7262,  -0.2019,   0.7068,   1.1421,   1.7653,   2.1040,   2.4752  &
/)

! ==============================================================================

end module ffdev_atomicdata_db

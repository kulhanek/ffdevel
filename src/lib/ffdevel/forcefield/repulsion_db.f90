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

module ffdev_repulsion_db

use ffdev_constants
use ffdev_repulsion_dat

! ==============================================================================
! ionization potential of elements, all in eV
! source: https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
! ref: CRC

real(DEVDP)     :: ionization_potential(1:18) = (/ &
13.59844, 24.58741,  5.39172,  9.32270,  8.29803, 11.26030, 14.53414, 13.61806, 17.42282, 21.5646, &
 5.13908,  7.64624,  5.98577,  8.15169, 10.48669, 10.36001, 12.96764, 15.75962 &
/)

! ==============================================================================

! calculated by gaussian, electron density symmetrized

integer         :: repulsion_PBE0_def2QZVPP_maxZ = 18
real(DEVDP)     :: repulsion_PBE0_def2QZVPP_b0(1:18) = (/ &
  2.8662,   4.9313,   1.6671,   2.2593,   2.6735,   3.1934,   3.7611,   4.1103,   4.6049,   5.0782, &
  1.8265,   2.2678,   2.2114,   2.6141,   2.9995,   3.2670,   3.5915,   4.0491 &
/)
real(DEVDP)     :: repulsion_PBE0_def2QZVPP_a0(1:18) = (/ &
 -3.6253,   0.0329,  -3.4905,  -1.7308,  -0.6980,   0.2912,   1.2173,   1.7185,   2.4225,   3.0735, &
 -1.4312,  -0.0893,  -0.1280,   0.9381,   1.7120,   2.2309,   2.7892,   3.7182 &
/)
real(DEVDP)     :: repulsion_PBE0_def2QZVPP_bm(1:18) = (/ &
  1.8156,   3.1076,   0.9885,   1.4491,   1.8984,   2.5109,   3.0250,   3.4741,   3.9990,   3.5671, &
  1.0793,   1.2725,   1.6271,   2.2451,   2.5645,   2.8504,   3.1708,   3.1367 &
/)
real(DEVDP)     :: repulsion_PBE0_def2QZVPP_am(1:18) = (/ &
 -3.7372,  -0.5703,  -4.4464,  -2.8069,  -1.5045,  -0.1585,   0.8239,   1.5624,   2.3400,   2.2042, &
 -2.7989,  -2.0164,  -0.8956,   0.8836,   1.4854,   2.0484,   2.6019,   3.1434 &
/)
real(DEVDP)     :: repulsion_PBE0_def2QZVPP_bp(1:18) = (/ &
  0.0000,   6.1748,   3.1610,   3.0346,   3.8767,   4.2199,   4.6646,   5.1174,   5.4225,   5.8875, &
  2.1916,   2.9533,   3.2884,   3.2224,   3.5468,   3.9306,   4.1884,   4.5738 &
/)
real(DEVDP)     :: repulsion_PBE0_def2QZVPP_ap(1:18) = (/ &
  0.0000,  -1.0271, -12.4032,  -1.3689,   0.8785,   1.3050,   1.8458,   2.3144,   2.6083,   3.1406, &
-10.9977,   0.4484,   1.4350,   1.4073,   2.1393,   2.8722,   3.2817,   3.9813 &
/)

! ==============================================================================

end module ffdev_repulsion_db

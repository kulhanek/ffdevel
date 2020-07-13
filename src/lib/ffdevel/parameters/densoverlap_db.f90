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

module ffdev_densoverlap_db

use ffdev_constants
use ffdev_densoverlap_dat

! ==============================================================================

! calculated by orca

integer         :: densoverlap_PBE0_def2QZVPP_maxZ = 18
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_b0(1:18) = (/ &
  2.8530,   4.9073,   1.6671,   2.2593,   2.5524,   3.8065,   3.7506,   4.0987,   4.8721,   5.0781, &
  1.8265,   2.2677,   2.6670,   2.5424,   2.9917,   3.3898,   3.5202,   4.0491 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_a0(1:18) = (/ &
 -3.6546,  -0.0341,  -3.4906,  -1.7308,  -0.6302,   0.7411,   1.1856,   1.7406,   2.2396,   3.0736, &
 -1.4312,  -0.0894,   0.3405,   1.0696,   1.6845,   2.1861,   2.8202,   3.7182 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_bm(1:18) = (/ &
  1.7963,   3.0749,   0.9885,   1.4425,   1.9154,   2.4981,   3.1766,   3.3945,   3.9990,   3.9895, &
  1.0793,   1.3564,   1.5862,   2.2356,   2.4456,   2.8942,   3.1573,   2.6989 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_am(1:18) = (/ &
 -3.7791,  -0.6588,  -4.4463,  -2.7657,  -1.4830,  -0.2017,   0.7506,   1.6219,   2.3401,   2.8516, &
 -2.7989,  -1.9267,  -0.7590,   0.8438,   1.5534,   1.9728,   2.5544,   2.6976 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_bp(1:18) = (/ &
  0.0000,   6.1702,   3.1610,   3.0345,   3.8768,   4.5417,   4.5846,   5.1174,   5.6238,   5.9821, &
  2.1915,   2.9533,   3.2883,   3.0399,   3.5963,   3.9270,   4.2638,   4.6352 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_ap(1:18) = (/ &
  0.0000,  -1.0374, -12.4033,  -1.3689,   0.8787,   1.4751,   2.0065,   2.3144,   2.5853,   2.9977, &
-10.9979,   0.4483,   1.4350,   1.8124,   2.0562,   2.8596,   3.2719,   3.9094 &
/)

! ==============================================================================

end module ffdev_densoverlap_db

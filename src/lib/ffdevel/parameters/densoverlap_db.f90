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

! calculated by gaussian

integer         :: densoverlap_PBE0_def2QZVPP_maxZ = 18
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_b0(1:18) = (/ &
  2.8466,   4.5366,   1.1915,   1.9748,   2.9905,   2.9653,   3.5561,   4.0574,   4.0421,   4.4957, &
  1.2871,   1.6338,   2.3148,   2.2904,   2.7310,   3.1763,   3.2862,   3.7102 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_a0(1:18) = (/ &
 -3.6617,  -0.8222,  -6.1412,  -3.0073,  -0.6720,  -0.1259,   0.5768,   0.9957,   1.3163,   1.6035, &
 -5.6002,  -3.5562,  -1.5003,  -0.3900,   0.3829,   1.1497,   1.8194,   2.4244 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_bm(1:18) = (/ &
  1.5194,   1.4307,   0.6787,   1.4911,   2.0029,   1.8219,   2.6112,   2.9218,   2.7754,   1.4911, &
  0.7403,   0.6679,   1.7622,   1.6706,   2.0499,   2.3916,   2.6597,   1.8549 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_am(1:18) = (/ &
 -4.4059,  -3.8045,  -6.5686,  -3.8780,  -2.6969,  -2.2901,  -0.7692,  -0.5707,  -0.3557,  -1.1648, &
 -6.3000,  -5.3788,  -2.6076,  -1.8753,  -1.0743,  -0.4710,   0.7372,   0.2790 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_bp(1:18) = (/ &
  0.0000,   6.5008,  10.6752,   2.6575,   3.9098,   4.0446,   5.2910,   5.2744,   5.7705,   6.1833, &
  7.6799,   2.1765,   2.9032,   3.6605,   4.0319,   3.9057,   4.3498,   4.5297 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_ap(1:18) = (/ &
  0.0000,  -0.3776, -31.8672,  -2.9457,   0.8864,   1.6583,   2.1160,   2.6076,   2.8086,   2.8220, &
-34.6914,  -3.5482,  -0.4928,   1.5281,   1.8540,   2.5289,   3.0924,   3.9454 &
/)

! ==============================================================================

end module ffdev_densoverlap_db

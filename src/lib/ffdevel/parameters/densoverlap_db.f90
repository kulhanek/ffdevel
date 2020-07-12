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

integer         :: densoverlap_HF_DKH_ANORCC_maxZ = 18
real(DEVDP)     :: densoverlap_HF_DKH_ANORCC_b0(1:18) = (/ &
  2.9809,   4.9440,   1.6297,   2.1939,   2.6507,   3.7370,   3.8368,   4.0485,   5.0292,   5.0731, &
  1.7924,   2.1563,   2.6213,   2.5356,   2.9846,   3.4189,   3.4736,   3.9110 &
/)
real(DEVDP)     :: densoverlap_HF_DKH_ANORCC_a0(1:18) = (/ &
 -3.4630,  -0.1335,  -3.5656,  -1.8557,  -0.3557,   0.5534,   1.2660,   1.4751,   2.3744,   2.7897, &
 -1.4173,  -0.3225,   0.2748,   1.0914,   1.6823,   2.2653,   2.6864,   3.2422 &
/)
real(DEVDP)     :: densoverlap_HF_DKH_ANORCC_bm(1:18) = (/ &
  1.5595,   1.3942,   0.9121,   0.8532,   1.7486,   2.1494,   0.8466,   2.4868,   3.0382,   3.1373, &
  0.6317,   0.6628,   1.4253,   1.8880,   1.8756,   2.2758,   2.6368,   2.7814 &
/)
real(DEVDP)     :: densoverlap_HF_DKH_ANORCC_am(1:18) = (/ &
 -4.3113,  -3.9367,  -4.7399,  -5.0213,  -1.6450,  -1.2272,  -3.6053,  -0.5125,   0.1373,   1.5367, &
 -3.2175,  -4.4631,  -1.1810,  -0.2900,  -0.2757,   0.4973,   1.1468,   2.5274 &
/)
real(DEVDP)     :: densoverlap_HF_DKH_ANORCC_bp(1:18) = (/ &
  0.0000,   6.1447,   6.8097,   3.0236,   3.9154,   4.6829,   5.4320,   5.4961,   6.0861,   6.1309, &
  2.1721,   2.9367,   3.2214,   3.5658,   3.7659,   4.0860,   4.4140,   4.8815 &
/)
real(DEVDP)     :: densoverlap_HF_DKH_ANORCC_ap(1:18) = (/ &
  0.0000,  -1.5314,  -2.9304,  -1.3572,   1.0874,   1.8332,   2.5369,   3.0301,   3.3798,   3.7211, &
-11.9831,   0.4941,   1.3703,   2.0510,   2.4282,   3.3373,   3.7854,   4.6686 &
/)

! ==============================================================================

end module ffdev_densoverlap_db

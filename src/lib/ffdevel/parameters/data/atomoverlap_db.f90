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

module ffdev_atomoverlap_db

use ffdev_constants
use ffdev_atomoverlap_dat

! ==============================================================================

! like-only atoms --------------------------------------------------------------
! calculated by Horton

integer         :: densoverlap_PBE0_def2QZVPP_maxZ = 18
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_b0(1:18) = (/ &
  3.6464,   5.6182,   2.1287,   2.8275,   3.1490,   3.7455,   4.3807,   4.7596,   5.2731,   5.7001, &
  2.2057,   2.7548,   2.6870,   3.0485,   3.5056,   3.7726,   4.1119,   4.5279 &
/)
real(DEVDP)     :: densoverlap_PBE0_def2QZVPP_a0(1:18) = (/ &
 -3.3637,  -0.7481,  -3.2328,  -1.5942,  -1.0525,  -0.1065,   0.7363,   1.1655,   1.7232,   2.0859, &
 -1.5928,  -0.1990,  -0.2275,   0.4291,   1.2303,   1.6119,   2.0659,   2.6541 &
/)

integer         :: wfoverlap_PBE0_def2QZVPP_maxZ = 18
real(DEVDP)     :: wfoverlap_PBE0_def2QZVPP_b0(1:18) = (/ &
  3.6528,   5.7293,   1.9433,   2.7361,   3.1064,   3.7210,   4.3807,   4.7738,   5.3179,   5.7630, &
  1.9973,   2.5316,   2.5266,   2.9709,   3.4405,   3.7307,   4.0899,   4.5697 &
/)
real(DEVDP)     :: wfoverlap_PBE0_def2QZVPP_a0(1:18) = (/ &
 -0.0357,   1.4696,   0.9733,   2.1648,   2.5187,   3.0141,   3.4550,   3.6600,   3.9725,   4.1492, &
  1.8825,   3.0207,   3.2401,   3.8623,   4.3549,   4.6047,   4.8759,   5.3788 &
/)

! ==============================================================================

end module ffdev_atomoverlap_db

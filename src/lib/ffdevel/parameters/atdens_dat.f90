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

module ffdev_atdens_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

integer,parameter   :: ATDENS_HF_UGBS       = 1
integer,parameter   :: ATDENS_CC_UGBS       = 2
integer,parameter   :: ATDENS_HF_DKH_ANORCC = 3
integer,parameter   :: ATDENS_CC_DKH_ANORCC = 4

integer             :: atdens_source = ATDENS_HF_DKH_ANORCC

! atom densities ---------------------------------------------------------------

integer,parameter   :: ATDENS_MAX_Z = 86

! ==============================================================================

real(DEVDP)     :: atdens_HF_UGBS_b(1:86) = (/  &
 3.7790,  5.4950,  2.0630,  2.7756,  3.5105,  3.7221,  4.4116,  5.0522,  4.9972,  5.5710, &
 1.9765,  2.4461,  2.4855,  2.8736,  3.3728,  3.7569,  3.9158,  4.3394,  1.8049,  2.0923, &
 2.1631,  2.2384,  1.9466,  2.2685,  2.3643,  2.3072,  2.3339,  2.4904,  2.5066,  2.7373, &
 3.2279,  2.8367,  3.2496,  3.2315,  3.7190,  4.0200,  1.7793,  2.0097,  2.1371,  2.1630, &
 1.9280,  2.5274,  2.3813,  2.2776,  2.2487,  4.1119,  2.2476,  2.6243,  2.4803,  2.7154, &
 3.0467,  3.3461,  3.5456,  3.7296,  1.7083,  1.9236,  1.9844,  1.9211,  1.9347,  1.9326, &
 1.7659,  1.7725,  1.9554,  2.4825,  2.2017,  2.1997,  1.9957,  1.7874,  1.7939,  2.0198, &
 2.1702,  2.4591,  2.3944,  2.3406,  2.5850,  2.3060,  2.4838,  2.4917,  2.4901,  2.6335, &
 2.9619,  3.3629,  2.9591,  3.4743,  3.6900,  3.4836 &
/)

real(DEVDP)     :: atdens_HF_UGBS_a(1:86) = (/  &
-1.1453,  0.4309, -2.0449, -0.7484, -0.4090,  0.5174,  0.9975,  1.3293,  1.5321,  1.7863, &
-2.0320, -0.9133,  0.0920,  0.1719,  0.6045,  0.7326,  1.3344,  1.6471, -1.7673, -0.9173, &
-1.0749, -0.9473, -2.1726, -1.5930, -0.9963, -1.6556, -1.6320, -1.3085, -1.3056, -0.5940, &
-0.1228,  0.1810,  0.5936,  0.9089,  1.4086,  1.6172, -1.5587, -0.8757, -0.5650, -0.8042, &
-1.9435, -0.5892, -0.9196, -1.5044, -1.6767,  1.8139, -1.7196, -0.4776,  0.3832,  0.3537, &
 0.6806,  0.9385,  1.6708,  1.7176, -1.4506, -0.6551, -1.1852, -0.7289, -0.7236, -0.8035, &
-1.4294, -1.4275, -0.8605,  0.6070, -1.0997, -1.1589, -0.8829, -1.6732, -1.6690, -0.9136, &
-0.6503, -0.0061, -0.4971, -0.6772, -0.4272, -1.4232, -0.9576, -0.9204, -0.9410, -0.3290, &
-0.0120,  0.6045,  0.7199,  1.4508,  1.6885,  1.6156 &
/)

logical     :: atdens_HF_UGBS_avail(1:86) = (/  &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true. &
/)

! ==============================================================================

real(DEVDP)     :: atdens_HF_DKH_ANORCC_b(1:86) = (/  &
 3.7779,  5.5439,  2.0456,  2.7635,  3.0463,  4.2651,  4.3863,  5.0200,  5.6352,  5.6131, &
 1.9884,  2.4710,  2.7125,  2.8968,  3.3855,  3.8113,  3.9610,  4.3533,  1.7156,  2.1053, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  3.2626,  3.6764,  3.9910,  4.0824,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000 &
/)

real(DEVDP)     :: atdens_HF_DKH_ANORCC_a(1:86) = (/  &
-1.1428,  0.4977, -2.1124, -0.8276, -0.0360,  0.2028,  0.9616,  1.3211,  1.6026,  1.8473, &
-1.9967, -0.8334, -0.9050,  0.2330,  0.6335,  0.8893,  1.2227,  1.6779, -2.1748, -0.9353, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.6022,  1.0680,  1.3492,  1.7412,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000, &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000 &
/)

logical     :: atdens_HF_DKH_ANORCC_avail(1:86) = (/  &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
 .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true.,  .true., &
.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
.false., .false.,  .true.,  .true.,  .true.,  .true., .false., .false., .false., .false., &
.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
.false., .false., .false., .false., .false., .false., .false., .false., .false., .false., &
.false., .false., .false., .false., .false., .false. &
/)

! ------------------------------------------------------------------------------

end module ffdev_atdens_dat

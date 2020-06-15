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

integer,parameter   :: ATDENS_CCSD_UGBS = 1

integer             :: atdens_source = ATDENS_CCSD_UGBS

! atom densities ---------------------------------------------------------------

integer,parameter   :: ATDENS_MAX_Z = 86

! CCSD(Full)/UGBS

real(DEVDP)     :: atdens_CCSD_UGBS_b(1:ATDENS_MAX_Z) = (/  &
 3.7790,  5.3741,  2.0756,  2.9083,  2.9689,  3.6980,  4.2326,  4.2047,  4.8179,  5.3434,  &
 1.9757,  2.5588,  2.4174,  2.8571,  3.2696,  3.6437,  3.9208,  4.2437,  1.6666,  2.2007,  &
 2.3615,  2.4331,  2.2014,  2.4184,  2.5628,  2.6566,  2.4805,  2.6305,  2.6181,  2.9051,  &
 2.8131,  2.8424,  3.2301,  3.5458,  3.8806,  3.9938,  1.6469,  2.0686,  2.2876,  2.3910,  &
 2.4568,  2.5881,  2.4525,  2.4171,  2.4339,  3.9963,  2.3601,  2.7213,  2.4814,  2.7218,  &
 3.0361,  3.2820,  3.5456,  3.7077,  1.4922,  1.8731,  1.8765,  0.0000,  0.0000,  0.0000,  &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000  &
/)

real(DEVDP)     :: atdens_CCSD_UGBS_a(1:ATDENS_MAX_Z) = (/  &
-1.1453,  0.3410, -2.0266, -0.6422, -0.6747,  0.0190,  0.8117,  0.9757,  1.2757,  1.5976,  &
-2.0709, -0.6953, -0.5308, -0.4649,  0.4559,  0.5986,  1.0944,  1.5491, -2.5370, -0.6817,  &
-0.5918, -0.5439, -1.5041, -1.2384, -0.6529, -0.5415, -1.3916, -0.9798, -1.1011, -0.3416,  &
-0.7885, -0.0354,  0.5389,  0.8042,  1.2042,  1.5901, -2.5358, -0.8695, -0.4101, -0.2212,  &
-0.4298, -0.4592, -0.9656, -1.2059, -1.1332,  1.6917, -1.4895, -0.3646, -0.0477,  0.2468,  &
 0.6500,  0.8748,  1.6084,  1.6937, -2.8590, -1.1649, -1.6011,  0.0000,  0.0000,  0.0000,  &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  &
 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000  &
/)

! ------------------------------------------------------------------------------

end module ffdev_atdens_dat

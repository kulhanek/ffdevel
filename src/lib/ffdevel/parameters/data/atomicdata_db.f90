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
! ionization potential of elements, all in eV
! source: https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
! ref: CRC

integer,parameter   :: IONIZATION_POTENTIAL_MAXZ = 18

real(DEVDP)     :: ionization_potential(1:IONIZATION_POTENTIAL_MAXZ) = (/ &
13.59844, 24.58741,  &
 5.39172,  9.32270,  8.29803, 11.26030, 14.53414, 13.61806, 17.42282, 21.56460, &
 5.13908,  7.64624,  5.98577,  8.15169, 10.48669, 10.36001, 12.96764, 15.75962  &
/)

! ==============================================================================

! linearized atomic densities of spherically symmetrized atoms

integer,parameter   :: ATOMICDATA_RHO_MAX_Z = 86

real(DEVDP)     :: atomicdata_rho_bmii(1:ATOMICDATA_RHO_MAX_Z)
real(DEVDP)     :: atomicdata_rho_b0ii(1:ATOMICDATA_RHO_MAX_Z)
real(DEVDP)     :: atomicdata_rho_bpii(1:ATOMICDATA_RHO_MAX_Z)

real(DEVDP)     :: atomicdata_rho_amii(1:ATOMICDATA_RHO_MAX_Z)
real(DEVDP)     :: atomicdata_rho_a0ii(1:ATOMICDATA_RHO_MAX_Z)
real(DEVDP)     :: atomicdata_rho_apii(1:ATOMICDATA_RHO_MAX_Z)

! like-only atoms --------------------------------------------------------------
! calculated by Horton



! ==============================================================================

end module ffdev_atomicdata_db

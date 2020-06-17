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

module ffdev_buried_dat

use ffdev_constants
use ffdev_variables

! buried data ------------------------------------------------------------------

! buried atom info
type BURRIED_ATOM_TYPE
    real(DEVDP)         :: expave   ! average exposition, 0.0 - buried, 1.0 - exposed
    real(DEVDP)         :: expsig
    real(DEVDP)         :: weight   ! 0.0 - buried, 1.0 - exposed
    integer             :: num
end type BURRIED_ATOM_TYPE

type(BURRIED_ATOM_TYPE),allocatable :: buried_atoms(:)    ! ntypes

integer,parameter   :: SURF_SESA = 1
integer,parameter   :: SURF_SASA = 2

integer             :: surface_mode = SURF_SESA
real(DEVDP)         :: ProbeR      =  1.30d0             ! solvent probe radius
real(DEVDP)         :: BuriedExp0  =  0.25d0
real(DEVDP)         :: BuriedBeta  = 50.00d0

! ------------------------------------------------------------------------------

end module ffdev_buried_dat

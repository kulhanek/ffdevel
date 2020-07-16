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

module ffdev_topology_utils

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! function ffdev_topology_find_nbtype
! ==============================================================================

integer function ffdev_topology_find_nbtype_by_aindex(top,ai,aj)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ai,aj
    ! --------------------------------------------
    integer         :: ti,tj
    ! --------------------------------------------------------------------------

    ! convert to types
    ti = top%atoms(ai)%typeid
    tj = top%atoms(aj)%typeid

    ffdev_topology_find_nbtype_by_aindex = ffdev_topology_find_nbtype_by_tindex(top,ti,tj)

end function ffdev_topology_find_nbtype_by_aindex

! ==============================================================================
! function ffdev_topology_find_nbtype_by_tindex
! ==============================================================================

integer function ffdev_topology_find_nbtype_by_tindex(top,ti,tj)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ti,tj
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,top%nnb_types
        if(  ( (top%nb_types(i)%ti .eq. ti) .and. (top%nb_types(i)%tj .eq. tj) ) .or. &
             ( (top%nb_types(i)%ti .eq. tj) .and. (top%nb_types(i)%tj .eq. ti) ) ) then
             ffdev_topology_find_nbtype_by_tindex = i
             return
        end if
    end do

    ! not found
    ffdev_topology_find_nbtype_by_tindex = 0
    return

end function ffdev_topology_find_nbtype_by_tindex

! ==============================================================================
! function ffdev_topology_find_nbtype_by_types
! ==============================================================================

integer function ffdev_topology_find_nbtype_by_types(top,sti,stj)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    character(*)    :: sti,stj
    ! --------------------------------------------
    integer         :: i,ti,tj
    ! --------------------------------------------------------------------------

    ti = 0
    tj = 0

    do i=1,top%natom_types
        if( top%atom_types(i)%name .eq. sti ) ti = i
        if( top%atom_types(i)%name .eq. stj ) tj = i
    end do

    ffdev_topology_find_nbtype_by_types = ffdev_topology_find_nbtype_by_tindex(top,ti,tj)

end function ffdev_topology_find_nbtype_by_types

! ------------------------------------------------------------------------------

end module ffdev_topology_utils


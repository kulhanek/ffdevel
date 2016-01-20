! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_genpoints_utils

use ffdev_sizes
use ffdev_constants

implicit none
contains

!===============================================================================
! subroutine genpoints_genrot
!===============================================================================

subroutine genpoints_genrot(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_genpoints_dat

    implicit none
    type(TOPOLOGY)          :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,nrotors
        call genpoints_genrot_mask(top,rotors(i)%ai,rotors(i)%aj)
        call genpoints_genrot_for_dihedral(geo,rotors(i)%ai,rotors(i)%aj,rotors(i)%angle)
    end do

end subroutine genpoints_genrot

!===============================================================================
! subroutine genpoints_genrot_mask
!===============================================================================

subroutine genpoints_genrot_mask(top,ai,aj)

    use ffdev_topology
    use ffdev_genpoints_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ai,aj
    ! --------------------------------------------
    integer         :: stacktop, i, is, js
    ! --------------------------------------------------------------------------

    ! reset mask
    AtomMask(:) = .false.

    AtomMask(ai) = .true. ! do not traverse in this direction
    AtomMask(aj) = .true.

    ! put first atom on the stack
    stacktop = 1
    ProcessingStack(stacktop) = aj

    do while( stacktop .gt. 0 )
        ! get atom from the stack
        is = ProcessingStack(stacktop)
        stacktop = stacktop - 1

        ! for all atom neighbours
        do i=1,top%atoms(is)%nbonds
            js = top%atoms(is)%bonded(i)
            if( AtomMask(js) ) cycle    ! already processed

            AtomMask(js) = .true.
            stacktop = stacktop + 1
            ProcessingStack(stacktop) = js
        end do
    end do

end subroutine genpoints_genrot_mask

!===============================================================================
! subroutine genpoints_genrot_for_dihedral
!===============================================================================

subroutine genpoints_genrot_for_dihedral(geo,ai,aj,dangle)

    use ffdev_geometry
    use ffdev_genpoints_maths
    use ffdev_genpoints_dat

    implicit none
    type(GEOMETRY)          :: geo
    integer                 :: ai,aj
    real(DEVDP)             :: dangle
    ! --------------------------------------------
    real(DEVDP)             :: dir(3)
    type(TRANSFORMATION)    :: trans
    ! --------------------------------------------------------------------------

    call trans%set_identity()

    dir(:) = - geo%crd(:,ai)
    call trans%move_by(dir)

    dir(:) = geo%crd(:,aj) - geo%crd(:,ai)
    call normalize(dir)
    call trans%rotate_by(dangle,dir)

    dir(:) = geo%crd(:,ai)
    call trans%move_by(dir)

    call trans%apply_to(geo%crd,AtomMask)

end subroutine genpoints_genrot_for_dihedral

!===============================================================================

end module ffdev_genpoints_utils

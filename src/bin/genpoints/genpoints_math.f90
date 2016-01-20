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
! WARRANTY without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FFDevel. If not, see <http://www.gnu.org/licenses/>.
! ==============================================================================

module ffdev_genpoints_maths

use ffdev_sizes
use ffdev_constants
implicit none

! ------------------------------------------------------------------------------

type TRANSFORMATION
    real(DEVDP)     :: Matrix(4,4)  ! transformation matrix

    contains
        ! executive methods
        procedure   :: set_zero
        procedure   :: set_identity
        procedure   :: mult_from_right
        procedure   :: move_by
        procedure   :: rotate_by
        procedure   :: apply_to
end type TRANSFORMATION

contains

!===============================================================================
! subroutine set_zero
!===============================================================================

subroutine set_zero(trans)

    implicit none
    class(TRANSFORMATION)   :: trans
    ! --------------------------------------------
    integer                 :: i,j
    ! --------------------------------------------------------------------------

    do i=1,4
        do j=1,4
            trans%Matrix(i,j) = 0.0d0
        end do
    end do

end subroutine set_zero

!===============================================================================
! subroutine set_identity
!===============================================================================

subroutine set_identity(trans)

    implicit none
    class(TRANSFORMATION)   :: trans
    ! --------------------------------------------
    integer                 :: i,j
    ! --------------------------------------------------------------------------

    do i=1,4
        do j=1,4
            if( i .eq. j ) then
                trans%Matrix(i,j) = 1.0d0
            else
                trans%Matrix(i,j) = 0.0d0
            end if
        end do
    end do

end subroutine set_identity

!===============================================================================
! subroutine set_identity
!===============================================================================

subroutine mult_from_right(trans1,trans2)

    implicit none
    class(TRANSFORMATION)   :: trans1
    type(TRANSFORMATION)    :: trans2
    ! --------------------------------------------
    integer                 :: i,j,k
    real(DEVDP)             :: pomvec(4), d
    ! --------------------------------------------------------------------------

    do i=1,4
        ! backup row
        do k=1,4
            pomvec(k) = trans1%Matrix(i,k)
        end do
        ! mult
        do j=1,4
            d = 0.0d0
            do k=1,4
                d = d + pomvec(k)*trans2%Matrix(k,j)
            end do
            trans1%Matrix(i,j) = d
        end do
    end do

end subroutine mult_from_right

!===============================================================================
! subroutine move_by
!===============================================================================

subroutine move_by(trans,mov)

    implicit none
    class(TRANSFORMATION)   :: trans
    real(DEVDP)             :: mov(3)
    ! --------------------------------------------
    type(TRANSFORMATION)    :: tmp
    ! --------------------------------------------------------------------------

    call tmp%set_identity()
    tmp%Matrix(4,1:3) = mov(1:3)
    call trans%mult_from_right(tmp)

end subroutine move_by

!===============================================================================
! subroutine move_by
!===============================================================================

subroutine rotate_by(trans,alpha,dir)

    implicit none
    class(TRANSFORMATION)   :: trans
    real(DEVDP)             :: alpha
    real(DEVDP)             :: dir(3)
    ! --------------------------------------------
    type(TRANSFORMATION)    :: tmp
    real(DEVDP)             :: d,pd
    ! --------------------------------------------------------------------------

    d  = sqrt(dir(1)**2 + dir(2)**2 + dir(3)**2)
    pd = sqrt(dir(1)**2 + dir(2)**2)

    call tmp%set_identity()

    if( pd .gt. 0.0d0 ) then                 ! rotate around z
        tmp%Matrix(1,1) =  dir(1) / pd
        tmp%Matrix(2,2) =  tmp%Matrix(1,1)
        tmp%Matrix(1,2) = -dir(2) / pd
        tmp%Matrix(2,1) = -tmp%Matrix(1,2)

        call trans%mult_from_right(tmp)
        call tmp%set_identity()
    end if

    if( d .gt. 0.0d0 ) then                 ! rotate around y
        tmp%Matrix(1,1) = dir(3) / d
        tmp%Matrix(3,3) = tmp%Matrix(1,1)
        tmp%Matrix(1,3) = pd / d
        tmp%Matrix(3,1) = -tmp%Matrix(1,3)

        call trans%mult_from_right(tmp)
        call tmp%set_identity()
    end if

! rotate around z about angle aplha

    tmp%Matrix(1,1) = cos(alpha)
    tmp%Matrix(2,2) = tmp%Matrix(1,1)
    tmp%Matrix(1,2) = sin(alpha)
    tmp%Matrix(2,1) = -tmp%Matrix(1,2)

    call trans%mult_from_right(tmp)
    call tmp%set_identity()

! backward transformation

    if( pd .gt. 0.0d0 ) then                 ! rotate around y
        tmp%Matrix(1,1) =  dir(3) / d
        tmp%Matrix(3,3) =  tmp%Matrix(1,1)
        tmp%Matrix(1,3) = -pd / d
        tmp%Matrix(3,1) = -tmp%Matrix(1,3)

        call trans%mult_from_right(tmp)
        call tmp%set_identity()
    end if

    if( d .gt. 0.0d0 ) then                ! rotate around z
        tmp%Matrix(1,1) =  dir(1) / pd
        tmp%Matrix(2,2) =  tmp%Matrix(1,1)
        tmp%Matrix(1,2) =  dir(2) / pd
        tmp%Matrix(2,1) = -tmp%Matrix(1,2)

        call trans%mult_from_right(tmp)
    end if

end subroutine rotate_by

!===============================================================================
! subroutine move_by
!===============================================================================

subroutine apply_to(trans,points,mask)

    implicit none
    class(TRANSFORMATION)   :: trans
    real(DEVDP)             :: points(:,:)
    logical                 :: mask(:)
    ! --------------------------------------------
    integer                 :: i,j
    real(DEVDP)             :: pomvec(4)
    ! --------------------------------------------------------------------------

    do i=1,size(points,2)
        if( .not. mask(i) ) cycle ! shall we rotate atom?

        pomvec(1:3) = points(1:3,i)
        pomvec(4) = 1.0d0

        ! prepare for final data
        points(:,i) = 0.0d0

        ! calculate transformation
        do j=1,4
            points(1:3,i) = points(1:3,i) + pomvec(j)*trans%Matrix(j,1:3)
        end do
    end do

end subroutine apply_to

!===============================================================================
! subroutine normalize
!===============================================================================

subroutine normalize(vec)

    implicit none
    real(DEVDP) :: vec(3)
    ! --------------------------------------------
    real(DEVDP) :: s
    ! --------------------------------------------------------------------------

    s = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    if( s .gt. 0.0d0 ) then
        vec(:) = vec(:) / s
    end if

end subroutine normalize

! ------------------------------------------------------------------------------

end module ffdev_genpoints_maths

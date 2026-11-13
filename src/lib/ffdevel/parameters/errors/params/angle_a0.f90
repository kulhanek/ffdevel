! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_angle_a0

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTAngleA0

    logical         :: BoltzmannWeights         ! extra Boltzmann weights from target energies

    contains
        ! executive methods
        procedure   :: init_errfce              => ffdev_err_angle_a0_init
        procedure   :: load_errfce              => ffdev_err_angle_a0_ctrl
        procedure   :: set_title_errfce         => ffdev_err_angle_a0_set_title
        procedure   :: calc_errfce              => ffdev_err_angle_a0_error
end type TypeEFTAngleA0

contains

! ==============================================================================
! subroutine ffdev_err_angle_a0_init
! ==============================================================================

subroutine ffdev_err_angle_a0_init(err_item)

    implicit none
    class(TypeEFTAngleA0)   :: err_item
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%init_errfce()

    err_item%BoltzmannWeights = .false.

end subroutine ffdev_err_angle_a0_init

! ==============================================================================
! subroutine ffdev_err_angle_a0_ctrl
! ==============================================================================

subroutine ffdev_err_angle_a0_ctrl(err_item,fin)

    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    class(TypeEFTAngleA0)       :: err_item
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%load_errfce(fin)

    if( prmfile_get_logical_by_key(fin,'bweight', err_item%BoltzmannWeights)) then
        write(DEV_OUT,170) prmfile_onoff(err_item%BoltzmannWeights)
    else
        write(DEV_OUT,175) prmfile_onoff(err_item%BoltzmannWeights)
    end if

170  format ('Extra Boltzmann weight (bweight)       = ',a12)
175  format ('Extra Boltzmann weight (bweight)       = ',a12,'                  (default)')

end subroutine ffdev_err_angle_a0_ctrl

!===============================================================================
! Subroutine:  ffdev_err_angle_a0_set_title
!===============================================================================

subroutine ffdev_err_angle_a0_set_title(err_item)

    implicit none
    class(TypeEFTAngleA0)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'angle_a0'

end subroutine ffdev_err_angle_a0_set_title

! ==============================================================================
! subroutine ffdev_err_angle_a0_error
! ==============================================================================

subroutine ffdev_err_angle_a0_error(err_item,opterr)

    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTAngleA0)    :: err_item
    logical                 :: opterr
    ! --------------------------------------------
    integer                 :: i,j,q,ai,aj,ak
    real(DEVDP)             :: err,seterrangles,totw
    real(DEVDP)             :: d0,dt,w
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    seterrangles = 0.0
    totw = 0

    do i=1,nsets
        do q=1,sets(i)%top%nangles
            if( err_item%OnlyFFOpt ) then
                if( .not. sets(i)%top%angle_types(sets(i)%top%angles(q)%at)%ffoptactive ) cycle
            end if
            ai = sets(i)%top%angles(q)%ai
            aj = sets(i)%top%angles(q)%aj
            ak = sets(i)%top%angles(q)%ak

            d0 = sets(i)%top%angle_types(sets(i)%top%angles(q)%at)%a0 * DEV_R2D

            do j=1,sets(i)%ngeos
                if( .not. sets(i)%geo(j)%trg_crd_loaded ) cycle

                dt = ffdev_geometry_get_angle(sets(i)%geo(j)%trg_crd,ai,aj,ak) * DEV_R2D
                err = d0 - dt
                w = sets(i)%geo(j)%weight
                if( err_item%BoltzmannWeights .and. (.not. sets(i)%isref) ) then
                    w = w * exp( -sets(i)%geo(j)%trg_energy / (DEV_Rgas*sets(i)%temperature))
                end if
                seterrangles = seterrangles + w * err**2
                totw = totw + w
            end do
        end do
    end do

    if( totw .gt. 0 ) then
        err_item%RepValue = sqrt(seterrangles/totw)
    end if


end subroutine ffdev_err_angle_a0_error

! ------------------------------------------------------------------------------

end module ffdev_err_angle_a0

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

module ffdev_err_rmsd

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTRMSD

    logical         :: MassWeighted

    contains
        ! executive methods
        procedure   :: init_errfce              => ffdev_err_rmsd_init
        procedure   :: load_errfce              => ffdev_err_rmsd_ctrl
        procedure   :: set_title_errfce         => ffdev_err_rmsd_set_title
        procedure   :: calc_errfce              => ffdev_err_rmsd_error
        procedure   :: print_set_summary_errfce => ffdev_err_rmsd_summary
end type TypeEFTRMSD

contains

!===============================================================================
! Subroutine:  ffdev_err_rmsd_init
!===============================================================================

subroutine ffdev_err_rmsd_init(err_item)

    implicit none
    class(TypeEFTRMSD)  :: err_item
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%init_errfce()

    err_item%MassWeighted = .true.

end subroutine ffdev_err_rmsd_init

! ==============================================================================
! subroutine ffdev_err_rmsd_ctrl
! ==============================================================================

subroutine ffdev_err_rmsd_ctrl(err_item,fin)

    use ffdev_utils
    use prmfile

    implicit none
    class(TypeEFTRMSD)  :: err_item
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%load_errfce(fin)

    if( prmfile_get_logical_by_key(fin,'mw', err_item%MassWeighted)) then
        write(DEV_OUT,140) prmfile_onoff(err_item%MassWeighted)
    else
        write(DEV_OUT,145) prmfile_onoff(err_item%MassWeighted)
    end if

140  format ('RMSD mass weighted (mw)                = ',a12)
145  format ('RMSD mass weighted (mw)                = ',a12,'                  (default)')

end subroutine ffdev_err_rmsd_ctrl

!===============================================================================
! Subroutine:  ffdev_err_rmsd_set_title
!===============================================================================

subroutine ffdev_err_rmsd_set_title(err_item)

    implicit none
    class(TypeEFTRMSD)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'RMSD'

end subroutine ffdev_err_rmsd_set_title

! ==============================================================================
! subroutine ffdev_err_rmsd_error
! ==============================================================================

subroutine ffdev_err_rmsd_error(err_item,opterr)

    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_geometry_utils
    use ffdev_errors_dat

    implicit none
    class(TypeEFTRMSD)  :: err_item
    logical             :: opterr
    ! --------------------------------------------
    integer             :: i,j
    real(DEVDP)         :: rmsd,seterrrmsd
    real(DEVDP)         :: swsum
    ! --------------------------------------------------------------------------

    err_item%ErrFceValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    seterrrmsd = 0.0
    swsum = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_crd_loaded .and. sets(i)%geo(j)%trg_crd_optimized ) then

                rmsd = ffdev_geometry_utils_get_rmsd(sets(i)%geo(j)%natoms,sets(i)%geo(j)%z, &
                                                     sets(i)%geo(j)%trg_crd,sets(i)%geo(j)%crd, &
                                                     err_item%MassWeighted)

                seterrrmsd = seterrrmsd + sets(i)%geo(j)%weight*rmsd**2
                swsum = swsum + 1
            end if
        end do
    end do

    if( swsum .gt. 0 ) then
        err_item%ErrFceValue = sqrt(seterrrmsd/real(swsum))
    end if

end subroutine ffdev_err_rmsd_error

! ==============================================================================
! subroutine ffdev_err_rmsd_summary
! ==============================================================================

subroutine ffdev_err_rmsd_summary(err_item,set,printsum)

    use ffdev_targetset_dat
    use ffdev_geometry_utils

    implicit none
    class(TypeEFTRMSD)  :: err_item
    type(TARGETSET)     :: set
    logical             :: printsum
    ! --------------------------------------------
    real(DEVDP)         :: rmsd,serr
    integer             :: j,num
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return

    if( printsum .eqv. .false. ) then
        do j=1,set%ngeos
            if( set%geo(j)%trg_crd_loaded .and. set%geo(j)%trg_crd_optimized ) then
                printsum = .true.
                return
            end if
        end do
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    serr = 0.0d0
    num = 0

    do j=1,set%ngeos
        if( set%geo(j)%trg_crd_loaded .and. set%geo(j)%trg_crd_optimized ) then

            rmsd = ffdev_geometry_utils_get_rmsd(set%geo(j)%natoms,set%geo(j)%z, &
                                                 set%geo(j)%trg_crd,set%geo(j)%crd,err_item%MassWeighted)

            serr = serr + set%geo(j)%weight * rmsd**2
            num = num + 1
            write(DEV_OUT,30) j, set%geo(j)%weight, rmsd
        end if
    end do

    if( num .gt. 0 ) then
        serr = sqrt(serr / real(num))
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40)  serr

    10 format('# ID   Weight       RMSD')
    20 format('# ---- ------ ----------')
    30 format(I6,1X,F6.3,1X,F10.3)
    40 format('# Final weighted error = ',F10.3)

end subroutine ffdev_err_rmsd_summary

! ------------------------------------------------------------------------------

end module ffdev_err_rmsd

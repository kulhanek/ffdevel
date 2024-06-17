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

module ffdev_err_aimr0

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_aimr0_init
! ==============================================================================

subroutine ffdev_err_aimr0_init

    use ffdev_err_aimr0_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableAIMR0Error        = .false.
    PrintAIMR0ErrorSummary  = .false.
    AIMR0ErrorWeight        = 1.0
    AIMR0BurriedOnly        = .false.

end subroutine ffdev_err_aimr0_init

! ==============================================================================
! subroutine ffdev_err_aimr0_error
! ==============================================================================

subroutine ffdev_err_aimr0_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_aimr0_dat
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat
    use ffdev_xdm_dat
    use ffdev_atomicdata_db

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer                 :: i,z
    real(DEVDP)             :: r0,r0opt,r0aim,v0,vaim,diff,err,tote,totw,w
    ! --------------------------------------------------------------------------

    error%aimr0 = 0.0d0

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        if( params(i)%realm .ne. REALM_VDW_R0 ) cycle

        if( params(i)%ti .ne. params(i)%tj ) then
            call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_aimr0_error!')
        end if

        if( AIMR0BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) cycle
        end if

        z = types(params(i)%ti)%z
        if( (z .lt. 1) .or. (z .gt. VDW_RADII_MAXZ) ) then
            call ffdev_utils_exit(DEV_ERR,1,'z .gt. VDW_RADII_MAXZ in ffdev_err_aimr0_error!')
        end if

        r0    = params(i)%value
        r0opt = atomicdata_r0opt(z)
        v0    = xdm_atoms(params(i)%ti)%v0ave
        vaim  = xdm_atoms(params(i)%ti)%vave

        r0aim = 2.0d0 * (vaim / v0)**(1.0d0/3.0d0) * r0opt

        diff = r0 - r0aim
        err = diff**2
        w = 1.0d0

        tote = tote + err
        totw = totw + w
    end do

    if( totw .gt. 0.0d0 ) then
        error%aimr0 = tote / totw
    end if

end subroutine ffdev_err_aimr0_error

! ==============================================================================
! subroutine ffdev_err_aimr0_summary
! ==============================================================================

subroutine ffdev_err_aimr0_summary()

    use ffdev_err_aimr0_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat
    use ffdev_xdm_dat
    use ffdev_atomicdata_db

    implicit none
    integer                 :: i,z
    real(DEVDP)             :: r0,r0opt,r0aim,v0,vaim,diff,err,tote,totw,w
    character(len=MAX_PATH) :: flag
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        if( params(i)%realm .ne. REALM_VDW_R0 ) cycle

        if( params(i)%ti .ne. params(i)%tj ) then
            call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_aimr0_error!')
        end if

        if( AIMR0BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) cycle
        end if

        z = types(params(i)%ti)%z
        if( (z .lt. 1) .or. (z .gt. VDW_RADII_MAXZ) ) then
            call ffdev_utils_exit(DEV_ERR,1,'z .gt. VDW_RADII_MAXZ in ffdev_err_aimr0_error!')
        end if

        r0    = params(i)%value
        r0opt = atomicdata_r0opt(z)
        v0    = xdm_atoms(params(i)%ti)%v0ave
        vaim  = xdm_atoms(params(i)%ti)%vave

        r0aim = 2.0d0 * (vaim / v0)**(1.0d0/3.0d0) * r0opt

        diff = r0 - r0aim

        err = diff**2
        w = 1.0d0
        flag = ''

        if( AIMR0BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) then
                err = 0.0d0
                w = 0.0d0
                flag = '   EXPOSED'
            else
                flag = 'BURIED'
            end if
        end if

        write(DEV_OUT,30) i,types(params(i)%ti)%name,types(params(i)%tj)%name,r0opt,v0,vaim,r0aim,r0, &
                          diff,err,trim(flag)

        tote = tote + err
        totw = totw + w
    end do

    if( totw .gt. 0.0d0 ) then
        tote = tote / totw
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40) tote
    write(DEV_OUT,45) tote*AIMR0ErrorWeight

 5 format('# NB R0 Penalties')
10 format('# ID TypA TypB R0Opt (DB)      V0       Vave         Raim      R0 (FF)     Diff       Error        Flag  ')
20 format('# -- ---- ---- ----------  ---------- ----------  ---------- ---------- ---------- ---------- ---------- ')
30 format(I4,1X,A4,1X,A4,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,A)
40 format('# Final penalty      =                          ',F10.5)
45 format('# Final penalty w/w  =                          ',F10.5)

end subroutine ffdev_err_aimr0_summary

! ------------------------------------------------------------------------------

end module ffdev_err_aimr0



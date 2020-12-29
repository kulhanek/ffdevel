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

module ffdev_err_nbr0

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_nbr0_init
! ==============================================================================

subroutine ffdev_err_nbr0_init

    use ffdev_err_nbr0_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableNBR0Error        = .false.
    PrintNBR0ErrorSummary  = .false.
    NBR0ErrorWeight        = 1.0
    NBR0BurriedOnly        = .true.

end subroutine ffdev_err_nbr0_init

! ==============================================================================
! subroutine ffdev_err_nbr0_error
! ==============================================================================

subroutine ffdev_err_nbr0_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_nbr0_dat
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer                 :: i
    real(DEVDP)             :: r0,rc,diff,err,tote,totw,w
    ! --------------------------------------------------------------------------

    error%nbr0 = 0.0d0

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        if( params(i)%realm .ne. REALM_VDW_R0 ) cycle

        if( params(i)%ti .ne. params(i)%tj ) then
            call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_nbr0_error!')
        end if

        r0 = params(i)%value
        rc = ffdev_atomicdata_rcii(params(i)%ti,damp_fa,damp_fb)

        diff = r0 - rc
        err = diff**2
        w = 1.0d0

        if( NBR0BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) cycle
        end if

        tote = tote + err
        totw = totw + w
    end do

    if( totw .gt. 0.0d0 ) then
        error%nbr0 = tote / totw
    end if

end subroutine ffdev_err_nbr0_error

! ==============================================================================
! subroutine ffdev_err_nbr0_summary
! ==============================================================================

subroutine ffdev_err_nbr0_summary()

    use ffdev_err_nbr0_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat

    implicit none
    integer                 :: i
    real(DEVDP)             :: r0,rc,diff,err,tote,totw,w
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
            call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_nbr0_summary!')
        end if

        r0 = params(i)%value
        rc = ffdev_atomicdata_rcii(params(i)%ti,damp_fa,damp_fb)

        diff = r0 - rc
        err = diff**2
        w = 1.0d0
        flag = ''

        if( NBR0BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) then
                err = 0.0d0
                w = 0.0d0
                flag = '   EXPOSED'
            else
                flag = 'BURIED'
            end if
        end if

        write(DEV_OUT,30) i,types(params(i)%ti)%name,types(params(i)%tj)%name,r0,rc, &
                          diff,err,trim(flag)

        tote = tote + err
        totw = totw + w
    end do

    if( totw .gt. 0.0d0 ) then
        tote = tote / totw
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40) tote
    write(DEV_OUT,45) tote*NBR0ErrorWeight

 5 format('# NB R0 Penalties')
10 format('# ID TypA TypB   R0 (FF)    RC (DB)     Diff       Error       Flag  ')
20 format('# -- ---- ---- ---------- ---------- ---------- ---------- ----------')
30 format(I4,1X,A4,1X,A4,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,A)
40 format('# Final penalty      =                          ',F10.5)
45 format('# Final penalty w/w  =                          ',F10.5)

end subroutine ffdev_err_nbr0_summary

! ------------------------------------------------------------------------------

end module ffdev_err_nbr0



! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2024 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_nbc6

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_nbc6_init
! ==============================================================================

subroutine ffdev_err_nbc6_init

    use ffdev_err_nbc6_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableNBC6Error        = .false.
    PrintNBC6ErrorSummary  = .false.
    NBC6ErrorWeight        = 1.0
    NBC6BurriedOnly        = .false.
    NBC6Sqrt16             = .true.
    NBC6Eff                = 0

end subroutine ffdev_err_nbc6_init

! ==============================================================================
! subroutine ffdev_err_nbc6_error
! ==============================================================================

subroutine ffdev_err_nbc6_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_nbc6_dat
    use ffdev_parameters_dat
    use ffdev_disp_dat
    use ffdev_buried_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer                 :: i,j
    real(DEVDP)             :: r0,eps,c6mm,c6ex,diff,err,tote,totw,w
    logical                 :: i_enabled,j_enabled
    ! --------------------------------------------------------------------------

    error%nbc6 = 0.0d0

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nparams

        if( params(i)%realm .ne. REALM_VDW_R0 ) cycle
        if( params(i)%ti .ne. params(i)%tj ) then
            call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_nbc6_error I (r0)!')
        end if

        if( NBC6BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) cycle
        end if

        r0  = params(i)%value
        eps = 0.0d0

        i_enabled = params(i)%enabled
        j_enabled = .false.

        do j=1,nparams
            if( params(j)%realm .ne. REALM_VDW_EPS ) cycle
            if( params(j)%ti .ne. params(j)%tj ) then
                call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_nbc6_error II (eps)!')
            end if
            if( params(i)%ti .eq. params(j)%ti ) then
                eps = params(j)%value
                j_enabled = params(j)%enabled
                exit
            end if
        end do

        if( .not. (i_enabled .or. j_enabled) ) cycle    ! at least one parameter, either r0 or eps must be enabled

        c6mm = 2.0d0 * eps * r0**6

        select case(NBC6Eff)
            case(0)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6
            case(1)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / disp_pairs(params(i)%ti,params(i)%ti)%rc**2
            case(2)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / disp_pairs(params(i)%ti,params(i)%ti)%rc**2 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c10 / disp_pairs(params(i)%ti,params(i)%ti)%rc**4
            case(3)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / r0**2
            case(4)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / r0**2
                if( types(params(i)%ti)%z .eq. 1 ) then
                    c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 / 2.0d0
                end if
            case default
                stop 'ffdev_err_nbc6_error'
        end select

        diff = c6mm - c6ex
        err = diff**2

        if( NBC6Sqrt16 ) then
            err = err ** (1.0d0/6.0d0)
        end if

        w = 1.0d0

        tote = tote + err
        totw = totw + w
    end do

    if( totw .gt. 0.0d0 ) then
        error%nbc6 = tote / totw
    end if

end subroutine ffdev_err_nbc6_error

! ==============================================================================
! subroutine ffdev_err_nbc6_summary
! ==============================================================================

subroutine ffdev_err_nbc6_summary()

    use ffdev_err_nbc6_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_disp_dat
    use ffdev_buried_dat

    implicit none
    integer                 :: i,j
    real(DEVDP)             :: c6mm,c6ex,diff,err,tote,totw,w,r0,eps
    character(len=MAX_PATH) :: flag
    logical                 :: i_enabled,j_enabled
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)

    if( NBC6Sqrt16 ) then
        write(DEV_OUT,6)
    end if

    select case(NBC6Eff)
        case(0)
            write(DEV_OUT,71)
        case(1)
            write(DEV_OUT,72)
        case(2)
            write(DEV_OUT,73)
        case(3)
            write(DEV_OUT,74)
        case(4)
            write(DEV_OUT,75)
        case default
            stop 'ffdev_err_nbc6_summary'
    end select

    write(DEV_OUT,10)
    write(DEV_OUT,20)

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nparams

       if( params(i)%realm .ne. REALM_VDW_R0 ) cycle
        if( params(i)%ti .ne. params(i)%tj ) then
            call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_nbc6_summary I (r0)!')
        end if

        r0  = params(i)%value
        eps = 0.0d0

        i_enabled = params(i)%enabled
        j_enabled = .false.

        do j=1,nparams
            if( params(j)%realm .ne. REALM_VDW_EPS ) cycle
            if( params(j)%ti .ne. params(j)%tj ) then
                call ffdev_utils_exit(DEV_ERR,1,'ti .ne. tj in ffdev_err_nbc6_summary II (eps)!')
            end if
            if( params(i)%ti .eq. params(j)%ti ) then
                eps = params(j)%value
                j_enabled = params(j)%enabled
                exit
            end if
        end do

        if( .not. (i_enabled .or. j_enabled) ) cycle    ! at least one parameter, either r0 or eps must be enabled

        c6mm = 2.0d0 * eps * r0**6

        select case(NBC6Eff)
            case(0)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6
            case(1)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / disp_pairs(params(i)%ti,params(i)%ti)%rc**2
            case(2)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / disp_pairs(params(i)%ti,params(i)%ti)%rc**2 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c10 / disp_pairs(params(i)%ti,params(i)%ti)%rc**4
            case(3)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / r0**2
            case(4)
                c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 + &
                       disp_pairs(params(i)%ti,params(i)%ti)%c8 / r0**2
                if( types(params(i)%ti)%z .eq. 1 ) then
                    c6ex = disp_pairs(params(i)%ti,params(i)%ti)%c6 / 2.0d0
                end if
        end select

        diff = c6mm - c6ex

        err = diff**2

        if( NBC6Sqrt16 ) then
            err = err ** (1.0d0/6.0d0)
        end if

        w = 1.0d0
        flag = ''

        if( NBC6BurriedOnly ) then
            if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) then
                err = 0.0d0
                w = 0.0d0
                flag = '   EXPOSED'
            else
                flag = 'BURIED'
            end if
        end if

        write(DEV_OUT,30) i,types(params(i)%ti)%name,types(params(i)%tj)%name,eps,r0,c6mm,c6ex, &
                          diff,err,trim(flag)

        tote = tote + err
        totw = totw + w
    end do

    if( totw .gt. 0.0d0 ) then
        tote = tote / totw
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40) tote
    write(DEV_OUT,45) tote*NBC6ErrorWeight

 5 format('# NB C6 Penalties')
 6 format('# err**(1/6) mode')
71 format('# C6eff = C6')
72 format('# C6eff = C6 + C8/Rvdw**2')
73 format('# C6eff = C6 + C8/Rvdw**2 + C10/Rvdw**4')
74 format('# C6eff = C6 + C8/R0**2')
75 format('# C6eff = C6 + C8/R0**2, C6eff = C6/2.0 (hydrogens)')
10 format('# ID TypA TypB  eps (FF)    R0 (FF)    C6 (FF)    C6 (DB)     Diff       Error       Flag  ')
20 format('# -- ---- ---- ---------- ---------- ---------- ---------- ---------- ---------- ----------')
30 format(I4,1X,A4,1X,A4,1X,F10.5,1X,F10.5,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.5,1X,A)
40 format('# Final penalty      =                                                ',F10.5)
45 format('# Final penalty w/w  =                                                ',F10.5)

end subroutine ffdev_err_nbc6_summary

! ------------------------------------------------------------------------------

end module ffdev_err_nbc6



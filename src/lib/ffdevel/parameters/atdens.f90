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

module ffdev_atdens

use ffdev_atdens_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_atdens_update_db
! ==============================================================================

subroutine ffdev_atdens_update_db

    use ffdev_utils
    use ffdev_atdens_db

    implicit none
    ! --------------------------------------------------------------------------

    atdens_bp(:) = 0.0d0
    atdens_b0(:) = 0.0d0
    atdens_bm(:) = 0.0d0

    atdens_ap(:) = 0.0d0
    atdens_a0(:) = 0.0d0
    atdens_am(:) = 0.0d0

    select case(atdens_source)
        case(ATDENS_HF_UGBS)
!            atdens_bm = atdens_HF_UGBS_bm
!            atdens_am = atdens_HF_UGBS_am
!
!            atdens_b0 = atdens_HF_UGBS_b0
!            atdens_a0 = atdens_HF_UGBS_a0
!
!            atdens_bp = atdens_HF_UGBS_bp
!            atdens_ap = atdens_HF_UGBS_ap

        case(ATDENS_HF_DKH_ANORCC)
            atdens_bm = atdens_HF_DKH_ANORCC_bm
            atdens_am = atdens_HF_DKH_ANORCC_am

            atdens_b0 = atdens_HF_DKH_ANORCC_b0
            atdens_a0 = atdens_HF_DKH_ANORCC_a0

            atdens_bp = atdens_HF_DKH_ANORCC_bp
            atdens_ap = atdens_HF_DKH_ANORCC_ap

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atdens_update_db')
    end select

end subroutine ffdev_atdens_update_db

! ==============================================================================
! subroutine ffdev_atdens_print
! ==============================================================================

subroutine ffdev_atdens_print

    use ffdev_utils
    use smf_periodic_table_dat
    use ffdev_parameters_dat
    use prmfile

    implicit none
    integer :: z, i
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Linearized Atomic Densities', '=')

    write(DEV_OUT,*)
    write(DEV_OUT,5) trim(ffdev_atdens_source_to_string(atdens_source))
    write(DEV_OUT,6) prmfile_onoff(atdens_mod_by_charge)

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

     do i=1,ntypes
        z = types(i)%z
        write(DEV_OUT,30,ADVANCE='NO') i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z))
        if( atdens_bm(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') atdens_bm(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( atdens_b0(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') atdens_b0(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( atdens_bp(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') atdens_bp(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if

        if( atdens_bm(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') atdens_am(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( atdens_b0(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') atdens_a0(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( atdens_bp(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') atdens_ap(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( atdens_mod_by_charge ) then
            write(DEV_OUT,42,ADVANCE='NO') types(i)%aveq
        else
            write(DEV_OUT,42,ADVANCE='NO') 0.0d0
        end if
        write(DEV_OUT,40) ffdev_atdens_b(i)
    end do

  5 format('# Atom density source  : ',A)
  6 format('# Modulation by charge : ',A)
 10 format('# ID Type  Z  El   B-     B0     B+     A-     A0     A+      <Q>    Bx  ')
 20 format('# -- ---- --- -- ------ ------ ------ ------ ------ ------  ------ ------')
 30 format(I4,1X,A4,1X,I3,1X,A2)
 40 format(1X,F6.3)
 42 format(2X,F6.3)
 50 format(7X)

end subroutine ffdev_atdens_print

! ==============================================================================
! function ffdev_atdens_b
! ==============================================================================

real(DEVDP) function ffdev_atdens_b(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    real(DEVDP) :: q, b1, b2
    ! --------------------------------------------------------------------------

    ffdev_atdens_b = 1.0 ! default b

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. ATDENS_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atdens_b')
    end if

    ! effective charge of type
    q = types(gti)%aveq

    ! no modulation by charge or zero charge or no extrapolation/interpolation data
    if( (.not. atdens_mod_by_charge) .or. (q .eq. 0.0d0) .or. &
        ( (atdens_bp(z) .eq. 0.0d0) .and. (atdens_bm(z) .eq. 0.0d0) ) ) then
        ffdev_atdens_b = atdens_b0(z)
        if( ffdev_atdens_b .eq. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'No atodens_b data for given Z in ffdev_atdens_b')
        end if
        return
    end if

    ! modulation by charge
        b1 = atdens_b0(z)
    if( q .gt. 0.0d0 ) then
        ! positive mode ( q > 0 )
        if( atdens_bp(z) .ne. 0.0d0 ) then
            ! interpolation
            b2 = atdens_bp(z) - atdens_b0(z)
        else
            ! extrapolation
            b2 = atdens_b0(z) - atdens_bm(z)
        end if
    else
        ! negative mode ( q < 0 )
        if( atdens_bm(z) .ne. 0.0d0 ) then
            ! interpolation
            b2 = atdens_b0(z) - atdens_bm(z)
        else
            ! extrapolation
            b2 = atdens_bp(z) - atdens_b0(z)
        end if
    end if

    ffdev_atdens_b = b1 + b2*q

end function ffdev_atdens_b

! ==============================================================================
! function ffdev_atdens_a
! ==============================================================================

real(DEVDP) function ffdev_atdens_a(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    real(DEVDP) :: q, a1, a2
    ! --------------------------------------------------------------------------

    ffdev_atdens_a = 0.0 ! default a

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. ATDENS_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atdens_a')
    end if

    ! effective charge of type
    q = types(gti)%aveq

    ! no modulation by charge or zero charge or no extrapolation/interpolation data
    if( (.not. atdens_mod_by_charge) .or. (q .eq. 0.0d0) .or. &
        ( (atdens_bp(z) .eq. 0.0d0) .and. (atdens_bm(z) .eq. 0.0d0) ) ) then
        ffdev_atdens_a = atdens_a0(z)
        if( ffdev_atdens_a .eq. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'No atodens_b data for given Z in ffdev_atdens_a')
        end if
        return
    end if

    ! modulation by charge
    a1 = atdens_a0(z)
    if( q .gt. 0.0d0 ) then
        ! positive mode ( q > 0 )
        if( atdens_bp(z) .ne. 0.0d0 ) then
            ! interpolation
            a2 = atdens_ap(z) - atdens_a0(z)
        else
            ! extrapolation
            a2 = atdens_a0(z) - atdens_am(z)
        end if
    else
        ! negative mode ( q < 0 )
        if( atdens_bm(z) .ne. 0.0d0 ) then
            ! interpolation
            a2 = atdens_a0(z) - atdens_am(z)
        else
            ! extrapolation
            a2 = atdens_ap(z) - atdens_a0(z)
        end if
    end if

    ffdev_atdens_a = a1 + a2*q

end function ffdev_atdens_a

! ==============================================================================
! function ffdev_atdens_rc
! ==============================================================================

real(DEVDP) function ffdev_atdens_rc(gti,dens)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    real(DEVDP) :: dens
    ! --------------------------------------------
    real(DEVDP) :: b, a
    ! --------------------------------------------------------------------------

    b = ffdev_atdens_b(gti)
    a = ffdev_atdens_a(gti)

    if( b .eq. 0.0d0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'No atodens_b data for given Z in ffdev_atdens_rc')
    end if

    ffdev_atdens_rc = (a - dens)/b

end function ffdev_atdens_rc

! ==============================================================================
! subroutine ffdev_atdens_source_from_string
! ==============================================================================

integer function ffdev_atdens_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('HF/UGBS')
            ffdev_atdens_source_from_string = ATDENS_HF_UGBS
        case('CC/UGBS')
            ffdev_atdens_source_from_string = ATDENS_CC_UGBS
        case('HF-DKH/ANORCC')
            ffdev_atdens_source_from_string = ATDENS_HF_DKH_ANORCC
        case('CC-DKH/ANORCC')
            ffdev_atdens_source_from_string = ATDENS_CC_DKH_ANORCC
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atdens_source_from_string!')
    end select

end function ffdev_atdens_source_from_string

! ==============================================================================
! subroutine ffdev_atdens_source_to_string
! ==============================================================================

character(80) function ffdev_atdens_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(ATDENS_HF_UGBS)
            ffdev_atdens_source_to_string = 'HF/UGBS'
        case(ATDENS_CC_UGBS)
            ffdev_atdens_source_to_string = 'CC/UGBS'
        case(ATDENS_HF_DKH_ANORCC)
            ffdev_atdens_source_to_string = 'HF-DKH/ANORCC'
        case(ATDENS_CC_DKH_ANORCC)
            ffdev_atdens_source_to_string = 'CC-DKH/ANORCC'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atdens_source_to_string!')
    end select

end function ffdev_atdens_source_to_string

! ------------------------------------------------------------------------------

end module ffdev_atdens

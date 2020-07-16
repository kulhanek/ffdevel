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

module ffdev_densoverlap

use ffdev_densoverlap_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_densoverlap_update_db
! ==============================================================================

subroutine ffdev_densoverlap_update_db

    use ffdev_utils
    use ffdev_densoverlap_db

    implicit none
    ! --------------------------------------------------------------------------

    densoverlap_bp(:) = 0.0d0
    densoverlap_b0(:) = 0.0d0
    densoverlap_bm(:) = 0.0d0

    densoverlap_ap(:) = 0.0d0
    densoverlap_a0(:) = 0.0d0
    densoverlap_am(:) = 0.0d0

    select case(densoverlap_source)
        case(DO_PBE0_def2QZVPP)
            densoverlap_bm(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_bm
            densoverlap_am(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_am

            densoverlap_b0(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_b0
            densoverlap_a0(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_a0

            densoverlap_bp(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_bp
            densoverlap_ap(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_ap

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_densoverlap_update_db')
    end select

end subroutine ffdev_densoverlap_update_db

! ==============================================================================
! subroutine ffdev_densoverlap_print
! ==============================================================================

subroutine ffdev_densoverlap_print

    use ffdev_utils
    use smf_periodic_table_dat
    use ffdev_parameters_dat
    use prmfile

    implicit none
    integer :: z, i
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Linearized Electron Density Overlaps', '=')

    write(DEV_OUT,*)
    write(DEV_OUT,5) trim(ffdev_densoverlap_source_to_string(densoverlap_source))
    write(DEV_OUT,6) prmfile_onoff(densoverlap_mod_by_charge)

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

     do i=1,ntypes
        z = types(i)%z
        write(DEV_OUT,30,ADVANCE='NO') i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z))
        if( densoverlap_bm(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') densoverlap_bm(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( densoverlap_b0(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') densoverlap_b0(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( densoverlap_bp(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') densoverlap_bp(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if

        if( densoverlap_bm(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') densoverlap_am(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( densoverlap_b0(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') densoverlap_a0(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( densoverlap_bp(z) .ne. 0 ) then
            write(DEV_OUT,40,ADVANCE='NO') densoverlap_ap(z)
        else
            write(DEV_OUT,50,ADVANCE='NO')
        end if
        if( densoverlap_mod_by_charge ) then
            write(DEV_OUT,42,ADVANCE='NO') types(i)%aveq
        else
            write(DEV_OUT,42,ADVANCE='NO') 0.0d0
        end if
        write(DEV_OUT,40) ffdev_densoverlap_bii(i)
    end do

  5 format('# Electron density overlap source : ',A)
  6 format('# Modulation by charge            : ',A)
 10 format('# ID Type  Z  El   B-     B0     B+     A-     A0     A+      <Q>    Bx  ')
 20 format('# -- ---- --- -- ------ ------ ------ ------ ------ ------  ------ ------')
 30 format(I4,1X,A4,1X,I3,1X,A2)
 40 format(1X,F6.3)
 42 format(2X,F6.3)
 50 format(7X)

end subroutine ffdev_densoverlap_print

! ==============================================================================
! function ffdev_densoverlap_bii
! ==============================================================================

real(DEVDP) function ffdev_densoverlap_bii(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    real(DEVDP) :: q, b1, b2
    ! --------------------------------------------------------------------------

    ffdev_densoverlap_bii = 1.0 ! default b

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_densoverlap_bii')
    end if

    ! effective charge of type
    q = types(gti)%aveq

    ! no modulation by charge or zero charge or no extrapolation/interpolation data
    if( (.not. densoverlap_mod_by_charge) .or. (q .eq. 0.0d0) .or. &
        ( (densoverlap_bp(z) .eq. 0.0d0) .and. (densoverlap_bm(z) .eq. 0.0d0) ) ) then
        ffdev_densoverlap_bii = densoverlap_b0(z)
        if( ffdev_densoverlap_bii .eq. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'No atodens_b data for given Z in ffdev_densoverlap_bii')
        end if
        return
    end if

    ! modulation by charge
        b1 = densoverlap_b0(z)
    if( q .gt. 0.0d0 ) then
        ! positive mode ( q > 0 )
        if( densoverlap_bp(z) .ne. 0.0d0 ) then
            ! interpolation
            b2 = densoverlap_bp(z) - densoverlap_b0(z)
        else
            ! extrapolation
            b2 = densoverlap_b0(z) - densoverlap_bm(z)
        end if
    else
        ! negative mode ( q < 0 )
        if( densoverlap_bm(z) .ne. 0.0d0 ) then
            ! interpolation
            b2 = densoverlap_b0(z) - densoverlap_bm(z)
        else
            ! extrapolation
            b2 = densoverlap_bp(z) - densoverlap_b0(z)
        end if
    end if

    ffdev_densoverlap_bii = b1 + b2*q

end function ffdev_densoverlap_bii

! ==============================================================================
! function ffdev_densoverlap_aii
! ==============================================================================

real(DEVDP) function ffdev_densoverlap_aii(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    real(DEVDP) :: q, a1, a2
    ! --------------------------------------------------------------------------

    ffdev_densoverlap_aii = 0.0 ! default a

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_densoverlap_aii')
    end if

    ! effective charge of type
    q = types(gti)%aveq

    ! no modulation by charge or zero charge or no extrapolation/interpolation data
    if( (.not. densoverlap_mod_by_charge) .or. (q .eq. 0.0d0) .or. &
        ( (densoverlap_bp(z) .eq. 0.0d0) .and. (densoverlap_bm(z) .eq. 0.0d0) ) ) then
        ffdev_densoverlap_aii = densoverlap_a0(z)
        if( ffdev_densoverlap_aii .eq. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'No atodens_b data for given Z in ffdev_densoverlap_aii')
        end if
        return
    end if

    ! modulation by charge
    a1 = densoverlap_a0(z)
    if( q .gt. 0.0d0 ) then
        ! positive mode ( q > 0 )
        if( densoverlap_bp(z) .ne. 0.0d0 ) then
            ! interpolation
            a2 = densoverlap_ap(z) - densoverlap_a0(z)
        else
            ! extrapolation
            a2 = densoverlap_a0(z) - densoverlap_am(z)
        end if
    else
        ! negative mode ( q < 0 )
        if( densoverlap_bm(z) .ne. 0.0d0 ) then
            ! interpolation
            a2 = densoverlap_a0(z) - densoverlap_am(z)
        else
            ! extrapolation
            a2 = densoverlap_ap(z) - densoverlap_a0(z)
        end if
    end if

    ffdev_densoverlap_aii = a1 + a2*q

end function ffdev_densoverlap_aii

! ==============================================================================
! function ffdev_densoverlap_rcii
! ==============================================================================

real(DEVDP) function ffdev_densoverlap_rcii(gti,dens)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    real(DEVDP) :: dens
    ! --------------------------------------------
    real(DEVDP) :: b, a
    ! --------------------------------------------------------------------------

    b = ffdev_densoverlap_bii(gti)
    a = ffdev_densoverlap_aii(gti)

    if( b .eq. 0.0d0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'No atodens_bii data for given gti in ffdev_densoverlap_rcii')
    end if

    ffdev_densoverlap_rcii = (a - dens)/b

end function ffdev_densoverlap_rcii

! ==============================================================================
! function ffdev_densoverlap_bij
! ==============================================================================

real(DEVDP) function ffdev_densoverlap_bij(gti,gtj)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti,gtj
    ! --------------------------------------------
    integer     :: zi,zj
    ! --------------------------------------------------------------------------

    ffdev_densoverlap_bij = 1.0 ! default b

    ! get Z
    zi = types(gti)%z
    if( (zi .le. 0) .and. (zi .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Zi is out-of-range in ffdev_densoverlap_bij')
    end if

    zj = types(gti)%z
    if( (zj .le. 0) .and. (zj .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Zj is out-of-range in ffdev_densoverlap_bij')
    end if

    ! FIXME

end function ffdev_densoverlap_bij

! ==============================================================================
! function ffdev_densoverlap_aij
! ==============================================================================

real(DEVDP) function ffdev_densoverlap_aij(gti,gtj)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti,gtj
    ! --------------------------------------------
    integer     :: zi,zj
    real(DEVDP) :: q, a1, a2
    ! --------------------------------------------------------------------------

    ffdev_densoverlap_aij = 0.0 ! default a

    ! get Z
    zi = types(gti)%z
    if( (zi .le. 0) .and. (zi .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Zi is out-of-range in ffdev_densoverlap_bij')
    end if

    zj = types(gti)%z
    if( (zj .le. 0) .and. (zj .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Zj is out-of-range in ffdev_densoverlap_bij')
    end if

    ! FIXME

end function ffdev_densoverlap_aij

! ==============================================================================
! function ffdev_densoverlap_rcij
! ==============================================================================

real(DEVDP) function ffdev_densoverlap_rcij(gti,gtj,dens)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti,gtj
    real(DEVDP) :: dens
    ! --------------------------------------------
    real(DEVDP) :: b, a
    ! --------------------------------------------------------------------------

    b = ffdev_densoverlap_bij(gti,gtj)
    a = ffdev_densoverlap_aij(gti,gtj)

    if( b .eq. 0.0d0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'No atodens_bij data for given gti,gtj in ffdev_densoverlap_rcij')
    end if

    ffdev_densoverlap_rcij = (a - dens)/b

end function ffdev_densoverlap_rcij

! ==============================================================================
! subroutine ffdev_densoverlap_source_from_string
! ==============================================================================

integer function ffdev_densoverlap_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('PBE0/def2-QZVPP')
            ffdev_densoverlap_source_from_string = DO_PBE0_def2QZVPP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_densoverlap_source_from_string!')
    end select

end function ffdev_densoverlap_source_from_string

! ==============================================================================
! subroutine ffdev_densoverlap_source_to_string
! ==============================================================================

character(80) function ffdev_densoverlap_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(DO_PBE0_def2QZVPP)
            ffdev_densoverlap_source_to_string = 'PBE0/def2-QZVPP'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_densoverlap_source_to_string!')
    end select

end function ffdev_densoverlap_source_to_string

! ------------------------------------------------------------------------------

end module ffdev_densoverlap

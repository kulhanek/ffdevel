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

module ffdev_atomicdata

use ffdev_atomicdata_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_atomicdata_update_db
! ==============================================================================

subroutine ffdev_atomicdata_update_db

    use ffdev_utils
    use ffdev_atomicdata_db

    implicit none
    ! --------------------------------------------------------------------------

    atomicdata_rho_bmii(:) = 0.0d0
    atomicdata_rho_b0ii(:) = 0.0d0
    atomicdata_rho_bpii(:) = 0.0d0

    atomicdata_rho_amii(:) = 0.0d0
    atomicdata_rho_a0ii(:) = 0.0d0
    atomicdata_rho_apii(:) = 0.0d0

    select case(bii_source)
        case(AD_BII_IP)
            ! nothing to do
        case(AD_BII_PBE0_def2QZVPP)
            ! FIXME
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_update_db')
    end select

end subroutine ffdev_atomicdata_update_db

! ==============================================================================
! subroutine ffdev_atomicdata_print
! ==============================================================================

subroutine ffdev_atomicdata_print

    use ffdev_utils
    use smf_periodic_table_dat
    use ffdev_parameters_dat
    use ffdev_atomicdata_db

    implicit none
    integer     :: i,z
    real(DEVDP) :: bmii,b0ii,bpii,amii,a0ii,apii,rvdw,ip
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Atomic Database', '=')

    write(DEV_OUT,*)
    write(DEV_OUT,5) trim(ffdev_atomicdata_bii_source_to_string(bii_source))
    write(DEV_OUT,6) trim(ffdev_atomicdata_bii_mods_to_string(bii_mods))
    write(DEV_OUT,7) trim(ffdev_atomicdata_rcii_source_to_string(rcii_source))

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

     do i=1,ntypes
        z = types(i)%z

        bmii = 0.0
        b0ii = 0.0
        bpii = 0.0
        amii = 0.0
        a0ii = 0.0
        apii = 0.0
        if( (z .ge. 1) .and. (z .le. ATOMICDATA_RHO_MAX_Z ) ) then
            bmii = atomicdata_rho_bmii(z)
            b0ii = atomicdata_rho_b0ii(z)
            bpii = atomicdata_rho_bpii(z)
            amii = atomicdata_rho_amii(z)
            a0ii = atomicdata_rho_a0ii(z)
            apii = atomicdata_rho_apii(z)
        end if

        rvdw = 0.0
        if( (z .ge. 1) .and. (z .le. VDW_RADII_MAXZ ) ) then
            rvdw = vdw_radii(z)
        end if

        ip = 0.0
        if( (z .ge. 1) .and. (z .le. IONIZATION_POTENTIAL_MAXZ ) ) then
            ip = ionization_potential(z)
        end if

        write(DEV_OUT,30) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          bmii, b0ii, bpii, amii, a0ii, apii, rvdw, ip

    end do
    write(DEV_OUT,20)

  5 format('# Bii data source  : ',A)
  6 format('# Bii data mods    : ',A)
  7 format('# Rcii data source : ',A)

 10 format('# ID Type  Z  El |   B-     B0     B+   |   A-     A0     A+   | R(vdW) | IP(eV) |')
 20 format('# -- ---- --- -- | ------ ------ ------ | ------ ------ ------ | ------ | ------ |')
 30 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3,3X,F6.3,1X,F6.3)

end subroutine ffdev_atomicdata_print

! ==============================================================================
! function ffdev_atomicdata_bii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_bii(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_xdm_dat
    use ffdev_atomicdata_db

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: maxz,z
    real(DEVDP) :: q
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_bii = 0.0 ! 0.0 -> undefined

    select case(bii_source)
        case(AD_BII_IP)
            maxz = IONIZATION_POTENTIAL_MAXZ
        case(AD_BII_PBE0_def2QZVPP)
            maxz = ATOMICDATA_RHO_MAX_Z
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii 0')
    end select

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. maxz) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_bii')
    end if
    q = types(gti)%aveq

    if( bii_mods .ne. AD_BII_MOD_BY_CHRG ) then
        ! charge unmodified bii
        select case(bii_source)
            case(AD_BII_IP)
                ! in A^-1 => 1.0 / DEV_AU2A
                ffdev_atomicdata_bii = 2.0d0 * sqrt(2.0d0 * ionization_potential(z) * DEV_eV2AU) / DEV_AU2A
            case(AD_BII_PBE0_def2QZVPP)
                ffdev_atomicdata_bii = atomicdata_rho_b0ii(z)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii I')
        end select

        select case(bii_mods)
            case(AD_BII_RAW)
                ! nothing to be here
            case(AD_BII_MOD_BY_XDM)
                if( .not. xdm_data_loaded ) then
                    call ffdev_utils_exit(DEV_ERR,1,'no XDM data loaded for bii in ffdev_atomicdata_bii')
                end if
                ! XDM mod, DOI: 10.1021/ct1001494
                ffdev_atomicdata_bii = ffdev_atomicdata_bii * &
                                    (xdm_atoms(gti)%v0ave / xdm_atoms(gti)%vave)**(1.0d0/3.0d0)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii II')
        end select
        ! exit
        return
    end if

    ! charge modified bii
    ! FIXME

end function ffdev_atomicdata_bii

! ==============================================================================
! function ffdev_atomicdata_aii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_aii(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomicdata_db

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: maxz,z
    real(DEVDP) :: q
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_aii = 0.0 ! 0.0 -> undefined

    select case(bii_source)
        case(AD_BII_IP)
            call ffdev_utils_exit(DEV_ERR,1,'aii is not available from IP in ffdev_atomicdata_aii')
        case(AD_BII_PBE0_def2QZVPP)
            maxz = ATOMICDATA_RHO_MAX_Z
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_aii 0')
    end select

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. maxz) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_aii')
    end if
    q = types(gti)%aveq

    if( bii_mods .ne. AD_BII_MOD_BY_CHRG ) then
        ! charge unmodified bii
        select case(bii_source)
            case(AD_BII_PBE0_def2QZVPP)
                ffdev_atomicdata_aii = atomicdata_rho_a0ii(z)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_aii I')
        end select

        select case(bii_mods)
            case(AD_BII_RAW)
                ! nothing to be here
            case(AD_BII_MOD_BY_XDM)
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported aii modification by XDM in ffdev_atomicdata_aii')
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_aii II')
        end select
        ! exit
        return
    end if

    ! charge modified aii
    ! FIXME

end function ffdev_atomicdata_aii

! ==============================================================================
! function ffdev_atomicdata_rcii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_rcii(gti,fa,fb)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_disp_dat
    use ffdev_atomicdata_db

    implicit none
    integer     :: gti
    real(DEVDP) :: fa,fb
    ! --------------------------------------------
    integer     :: maxz, z
    real(DEVDP) :: b, a
    real(DEVDP) :: q
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_rcii = 0.0 ! 0.0 -> undefined

    select case(rcii_source)
        case(AD_RCII_VDW)
            maxz = VDW_RADII_MAXZ
    ! -------------
        case(AD_RCII_XDM)
            ! workaround - largest Z we have
            maxz = ATOMICDATA_RHO_MAX_Z
    ! -------------
        case(AD_RCII_PBE0_def2QZVPP)
            maxz = ATOMICDATA_RHO_MAX_Z
    ! -------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_rcii 0')
    end select

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. maxz) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_rcii')
    end if
    q = types(gti)%aveq

    select case(rcii_source)
        case(AD_RCII_VDW)
            ffdev_atomicdata_rcii = fa * vdw_radii(z) + fb
    ! -------------
        case(AD_RCII_XDM)
            if( .not. disp_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1,'no DISP data loaded for rcc in ffdev_atomicdata_rcii')
            end if
            ffdev_atomicdata_rcii = fa * disp_pairs(gti,gti)%rc + fb
    ! -------------
        case(AD_RCII_PBE0_def2QZVPP)
            b = ffdev_atomicdata_bii(gti)
            a = ffdev_atomicdata_aii(gti)
            if( b .eq. 0.0d0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'No ffdev_atomicdata_bii data for given gti in ffdev_atomicdata_rcii')
            end if
            ffdev_atomicdata_rcii = (a - fa)/b
    ! -------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_rcii II')
    end select

end function ffdev_atomicdata_rcii

! ==============================================================================
! subroutine ffdev_atomicdata_bii_source_from_string
! ==============================================================================

integer function ffdev_atomicdata_bii_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('IP')
            ffdev_atomicdata_bii_source_from_string = AD_BII_IP
        case('PBE0/def2-QZVPP')
            ffdev_atomicdata_bii_source_from_string = AD_BII_PBE0_def2QZVPP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atomicdata_bii_source_from_string!')
    end select

end function ffdev_atomicdata_bii_source_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_bii_source_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_bii_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AD_BII_IP)
            ffdev_atomicdata_bii_source_to_string = 'IP'
        case(AD_BII_PBE0_def2QZVPP)
            ffdev_atomicdata_bii_source_to_string = 'PBE0/def2-QZVPP'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii_source_to_string!')
    end select

end function ffdev_atomicdata_bii_source_to_string

! ==============================================================================
! subroutine ffdev_atomicdata_bii_mods_from_string
! ==============================================================================

integer function ffdev_atomicdata_bii_mods_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('RAW')
            ffdev_atomicdata_bii_mods_from_string = AD_BII_RAW
        case('XDM')
            ffdev_atomicdata_bii_mods_from_string = AD_BII_MOD_BY_XDM
        case('CHARGE')
            ffdev_atomicdata_bii_mods_from_string = AD_BII_MOD_BY_CHRG
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atomicdata_bii_mods_from_string!')
    end select

end function ffdev_atomicdata_bii_mods_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_bii_mods_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_bii_mods_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AD_BII_RAW)
            ffdev_atomicdata_bii_mods_to_string = 'RAW'
        case(AD_BII_MOD_BY_XDM)
            ffdev_atomicdata_bii_mods_to_string = 'XDM'
        case(AD_BII_MOD_BY_CHRG)
            ffdev_atomicdata_bii_mods_to_string = 'CHARGE'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii_mods_to_string!')
    end select

end function ffdev_atomicdata_bii_mods_to_string

! ==============================================================================
! subroutine ffdev_atomicdata_rcii_source_from_string
! ==============================================================================

integer function ffdev_atomicdata_rcii_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('VDW')
            ffdev_atomicdata_rcii_source_from_string = AD_RCII_VDW
        case('XDM')
            ffdev_atomicdata_rcii_source_from_string = AD_RCII_XDM
        case('PBE0/def2-QZVPP')
            ffdev_atomicdata_rcii_source_from_string = AD_RCII_PBE0_def2QZVPP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atomicdata_rcii_source_from_string!')
    end select

end function ffdev_atomicdata_rcii_source_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_rcii_source_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_rcii_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AD_RCII_VDW)
            ffdev_atomicdata_rcii_source_to_string = 'VDW'
        case(AD_RCII_XDM)
            ffdev_atomicdata_rcii_source_to_string = 'XDM'
        case(AD_RCII_PBE0_def2QZVPP)
            ffdev_atomicdata_rcii_source_to_string = 'PBE0/def2-QZVPP (bii/aii)'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_rcii_source_to_string!')
    end select

end function ffdev_atomicdata_rcii_source_to_string

! ------------------------------------------------------------------------------

end module ffdev_atomicdata

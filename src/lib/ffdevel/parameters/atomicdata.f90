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

    atomicdata_rho_b012(:,:) = 0.0d0

    ATDENS_MAX_Z = 0

    select case(atdens_source)
        case(AD_ATDENS_PBE0_def2QZVPP)
            ATDENS_MAX_Z = atomicdata_rho_PBE0_def2QZVPP_maxZ
            atomicdata_rho_b012(1:atomicdata_rho_PBE0_def2QZVPP_maxZ,1:3) = &
                transpose(reshape(atomicdata_rho_PBE0_def2QZVPP_b012,(/3,atomicdata_rho_PBE0_def2QZVPP_maxZ/)))
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_update_db')
    end select

    if( .not. atomicdata_b0opt_initialized ) then
        atomicdata_b0opt(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = atomicdata_rho_b012(1:atomicdata_rho_PBE0_def2QZVPP_maxZ,1)
        atomicdata_b0opt_initialized = .true.
    end if

    atomicdata_vdw_r0free(1:VDW_RADII_MAXZ) = 2.0d0*vdw_radii(1:VDW_RADII_MAXZ)
    ! FIXME
    atomicdata_vdw_pbfree(1:VDW_RADII_MAXZ) = 3.2

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
    real(DEVDP) :: b0,b1,b2,rvdw,bii,rcii,r0free,pbfree
    real(DEVDP) :: ip0,ipp
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Atomic Database', '=')

! ------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(ffdev_atomicdata_atdens_source_to_string(atdens_source))
    write(DEV_OUT,110)
    write(DEV_OUT,120)

    do i=1,ntypes
        z = types(i)%z

        b0 = 0.0
        b1 = 0.0
        b2 = 0.0
        if( (z .ge. 1) .and. (z .le. ATDENS_MAX_Z ) ) then
            b0 = atomicdata_rho_b012(z,1)
            b1 = atomicdata_rho_b012(z,2)
            b2 = atomicdata_rho_b012(z,3)
        end if

        write(DEV_OUT,130) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          b0, b1, b2

    end do
    write(DEV_OUT,120)

! ------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,220)

    do i=1,ntypes
        z = types(i)%z

        ip0 = 0.0
        ipp = 0.0

        if( (z .ge. 1) .and. (z .le. IPEA_MAXZ ) ) then
            ip0 = atomicdata_ip0(z)
            ipp = atomicdata_ipp(z)
        end if

        ! FIXME - what to do with negative EA
        b0 = 2.0d0 * sqrt(2.0d0 * ip0 * DEV_eV2AU) / DEV_AU2A
        b1 = 2.0d0 * sqrt(2.0d0 * ipp * DEV_eV2AU) / DEV_AU2A - b0

        write(DEV_OUT,230) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          ip0, ipp, b0, b1
    end do
    write(DEV_OUT,220)

! ------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,300)
    write(DEV_OUT,310)
    write(DEV_OUT,320)

    do i=1,ntypes
        z = types(i)%z

        rvdw = 0.0
        if( (z .ge. 1) .and. (z .le. VDW_RADII_MAXZ ) ) then
            rvdw = 2.0d0 * vdw_radii(z)
        end if

        write(DEV_OUT,330) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          rvdw

    end do
    write(DEV_OUT,320)

! ------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,600)
    write(DEV_OUT,610)
    write(DEV_OUT,320)

    do i=1,ntypes
        z = types(i)%z

        pbfree = 0.0
        r0free = 0.0
        if( (z .ge. 1) .and. (z .le. VDW_RADII_MAXZ ) ) then
            pbfree = atomicdata_vdw_pbfree(z)
            r0free = atomicdata_vdw_r0free(z)
        end if

        write(DEV_OUT,630) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          pbfree, r0free

    end do
    write(DEV_OUT,320)

! ------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,410) trim(ffdev_atomicdata_bii_source_to_string(bii_source))
    write(DEV_OUT,420) trim(ffdev_atomicdata_bii_mods_to_string(bii_mods))
    write(DEV_OUT,430) trim(ffdev_atomicdata_rcii_source_to_string(rcii_source))

    write(DEV_OUT,*)
    write(DEV_OUT,500) damp_fa,damp_fb
    write(DEV_OUT,510)
    write(DEV_OUT,520)

    do i=1,ntypes
        z = types(i)%z

        bii = ffdev_atomicdata_bii(i)
        rcii = ffdev_atomicdata_rcii(i,damp_fa,damp_fb)

        write(DEV_OUT,530) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          types(i)%aveq, bii, rcii

    end do
    write(DEV_OUT,520)

100 format('# Bii/Aii from atomic densities: ',A)
110 format('# ID Type  Z  El |   B0    B1*q  B2*q^2 |')
120 format('# -- ---- --- -- | ------ ------ ------ |')
130 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3,1X,F6.3)

200 format('# Bii/Aii from ionization potentials and electron affinities')
210 format('# ID Type  Z  El |  IP0    IP+  -->  B0    B1*q  |')
220 format('# -- ---- --- -- | ------ ------ | ------ ------ |')
230 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3,3X,F6.3,1X,F6.3)

300 format('# van der Waals radii')
310 format('# ID Type  Z  El | 2Rvdw  |')
320 format('# -- ---- --- -- | ------ |')
330 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3)

600 format('# AIM Atomic Data')
610 format('# ID Type  Z  El | PBFREE | R0FREE |')
630 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3)

410 format('# Bii data source  : ',A)
420 format('# Bii data mods    : ',A)
430 format('# Rcii data source : ',A)

500 format('# Current data: damp_fa = ',F6.3,'; damp_fb= ',F6.3)
510 format('# ID Type  Z  El  Qeff  |   Bii   Rcii  |')
520 format('# -- ---- --- -- ------ | ------ ------ |')
530 format(I4,1X,A4,1X,I3,1X,A2,1X,F6.3,3X,F6.3,1X,F6.3)

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
    real(DEVDP) :: q, b0, b1, b2
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_bii = 0.0 ! 0.0 -> undefined

    select case(bii_source)
        case(AD_BII_IPEA)
            maxz = IPEA_MAXZ
        case(AD_BII_ATDENS)
            maxz = ATDENS_MAX_Z
        case(AD_BII_B0OPT)
            maxz = B0OPT_MAX_Z
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii 0')
    end select

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. maxz) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_bii')
    end if
    q = types(gti)%aveq

    select case(bii_source)
        case(AD_BII_IPEA)
            ! in A^-1 => 1.0 / DEV_AU2A
            b0 = 2.0d0 * sqrt(2.0d0 * atomicdata_ip0(z) * DEV_eV2AU) / DEV_AU2A
            b1 = 2.0d0 * sqrt(2.0d0 * atomicdata_ipp(z) * DEV_eV2AU) / DEV_AU2A - b0
            b2 = 0.0d0
    ! ------------
        case(AD_BII_ATDENS)
            b0 = atomicdata_rho_b012(z,1)
            b1 = atomicdata_rho_b012(z,2)
            b2 = atomicdata_rho_b012(z,3)
    ! ------------
        case(AD_BII_B0OPT)
            b0 = atomicdata_b0opt(z)
            b1 = 0.0
            b2 = 0.0
    ! ------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii I')
    end select

    select case(bii_mods)
        case(AD_BII_RAW)
            ffdev_atomicdata_bii = b0
            return
    ! ------------
        case(AD_BII_MOD_BY_XDM)
            if( .not. xdm_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1,'no XDM data loaded for bii in ffdev_atomicdata_bii')
            end if
            ! XDM mod, DOI: 10.1021/ct1001494
            ffdev_atomicdata_bii = b0 * &
                                (xdm_atoms(gti)%v0ave / xdm_atoms(gti)%vave)**(1.0d0/3.0d0)
            return
    ! ------------
        case(AD_BII_MOD_BY_CHRG)
            ffdev_atomicdata_bii = b0 + b1*q + b2*q**2
    ! ------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii II')
    end select

end function ffdev_atomicdata_bii

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
    real(DEVDP) :: q
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_rcii = 0.0 ! 0.0 -> undefined

    select case(rcii_source)
        case(AD_RCII_VDW)
            maxz = VDW_RADII_MAXZ
    ! -------------
        case(AD_RCII_DISP)
            ! workaround - largest Z we have
            maxz = ATDENS_MAX_Z_ALL
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
            ffdev_atomicdata_rcii = fa * 2.0d0 * vdw_radii(z) + fb
    ! -------------
        case(AD_RCII_DISP)
            if( .not. disp_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1,'no DISP data loaded for rcc in ffdev_atomicdata_rcii')
            end if
            ffdev_atomicdata_rcii = fa * disp_pairs(gti,gti)%rc + fb
    ! -------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_rcii II')
    end select

end function ffdev_atomicdata_rcii

! ------------------------------------------------------------------------------

real(DEVDP) function ffdev_atomicdata_get_effZ(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------------------------------------

    if( ntypes .eq. 0 ) then
        ffdev_atomicdata_get_effZ = 0.0
        return
    end if

    select case(Zeff_mode)
        case(AD_ZEFF_MAX)
            ffdev_atomicdata_get_effZ = types(gti)%z
            return
    ! --------------------
        case(AD_ZEFF_VALENCE)
            ffdev_atomicdata_get_effZ = ffdev_atomicdata_get_valence_effZ(gti)
            return
    ! --------------------
        case(AD_ZEFF_OPT)
            ffdev_atomicdata_get_effZ = types(gti)%Zeff
            return
    ! --------------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_get_effZ')
    end select

end function ffdev_atomicdata_get_effZ

! ------------------------------------------------------------------------------

real(DEVDP) function ffdev_atomicdata_get_min_Zeff(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomicdata_db

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. ZEFF_CLEMENTI_MAXZ) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_get_min_Zeff')
    end if

    ffdev_atomicdata_get_min_Zeff = zeff_clementi(z)

end function ffdev_atomicdata_get_min_Zeff

! ------------------------------------------------------------------------------

real(DEVDP) function ffdev_atomicdata_get_valence_effZ(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    z = types(gti)%z
    ! determine number of valence electrons
    if( z .le. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_get_valence_effZ')
    end if
    if( z .le. 2 ) then
        ! H-He
        ffdev_atomicdata_get_valence_effZ = z
        return
    end if
    if( z .le. 10 ) then
        ! Li-Ne
        ffdev_atomicdata_get_valence_effZ = z - 2
        return
    end if
    if( z .le. 18 ) then
        ! Na-Ar
        ffdev_atomicdata_get_valence_effZ = z - 10
        return
    end if
    if( z .le. 36 ) then
        ! K-Kr
        ffdev_atomicdata_get_valence_effZ = z - 18
        return
    end if
   if( z .le. 54 ) then
        ! Rb-Xe
        ffdev_atomicdata_get_valence_effZ = z - 36
        return
   end if
   if( z .le. 86 ) then
        ! Cs-Rn
        ffdev_atomicdata_get_valence_effZ = z - 54
        return
    end if
    call ffdev_utils_exit(DEV_ERR,1,'Not defined in ffdev_atomicdata_get_valence_effZ')

end function ffdev_atomicdata_get_valence_effZ

! ==============================================================================
! subroutine ffdev_atomicdata_bii_source_from_string
! ==============================================================================

integer function ffdev_atomicdata_bii_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('IPEA')
            ffdev_atomicdata_bii_source_from_string = AD_BII_IPEA
        case('ATDENS')
            ffdev_atomicdata_bii_source_from_string = AD_BII_ATDENS
        case('B0OPT')
            ffdev_atomicdata_bii_source_from_string = AD_BII_B0OPT
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
        case(AD_BII_IPEA)
            ffdev_atomicdata_bii_source_to_string = 'IPEA'
        case(AD_BII_ATDENS)
            ffdev_atomicdata_bii_source_to_string = 'ATDENS'
        case(AD_BII_B0OPT)
            ffdev_atomicdata_bii_source_to_string = 'B0OPT'
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
! subroutine ffdev_atomicdata_atdens_source_from_string
! ==============================================================================

integer function ffdev_atomicdata_atdens_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('PBE0/def2-QZVPP')
            ffdev_atomicdata_atdens_source_from_string = AD_ATDENS_PBE0_def2QZVPP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) &
                                  //'" in ffdev_atomicdata_atdens_source_from_string!')
    end select

end function ffdev_atomicdata_atdens_source_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_atdens_source_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_atdens_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AD_ATDENS_PBE0_def2QZVPP)
            ffdev_atomicdata_atdens_source_to_string = 'PBE0/def2-QZVPP'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_atdens_source_to_string!')
    end select

end function ffdev_atomicdata_atdens_source_to_string

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
        case('DISP')
            ffdev_atomicdata_rcii_source_from_string = AD_RCII_DISP
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
        case(AD_RCII_DISP)
            ffdev_atomicdata_rcii_source_to_string = 'DISP'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_rcii_source_to_string!')
    end select

end function ffdev_atomicdata_rcii_source_to_string

! ==============================================================================
! subroutine ffdev_atomicdata_zeff_mode_from_string
! ==============================================================================

integer function ffdev_atomicdata_zeff_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('MAX')
            ffdev_atomicdata_zeff_mode_from_string = AD_ZEFF_MAX
        case('VALENCE')
            ffdev_atomicdata_zeff_mode_from_string = AD_ZEFF_VALENCE
        case('OPT')
            ffdev_atomicdata_zeff_mode_from_string = AD_ZEFF_OPT
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atomicdata_zeff_mode_from_string!')
    end select

end function ffdev_atomicdata_zeff_mode_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_zeff_mode_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_zeff_mode_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AD_ZEFF_MAX)
            ffdev_atomicdata_zeff_mode_to_string = 'MAX'
        case(AD_ZEFF_VALENCE)
            ffdev_atomicdata_zeff_mode_to_string = 'VALENCE'
        case(AD_ZEFF_OPT)
            ffdev_atomicdata_zeff_mode_to_string = 'OPT'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_zeff_mode_to_string!')
    end select

end function ffdev_atomicdata_zeff_mode_to_string

! ------------------------------------------------------------------------------

end module ffdev_atomicdata

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

    atomicdata_rho_bm(:) = 0.0d0
    atomicdata_rho_b0(:) = 0.0d0
    atomicdata_rho_bp(:) = 0.0d0

    atomicdata_rho_am(:) = 0.0d0
    atomicdata_rho_a0(:) = 0.0d0
    atomicdata_rho_ap(:) = 0.0d0

    ATDENS_MAX_Z = 0

    select case(atdens_source)
        case(AD_ATDENS_PBE0_def2QZVPP)
            ATDENS_MAX_Z = atomicdata_rho_PBE0_def2QZVPP_maxZ
            atomicdata_rho_bm(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = &
                       atomicdata_rho_PBE0_def2QZVPP_bm(1:atomicdata_rho_PBE0_def2QZVPP_maxZ)
            atomicdata_rho_b0(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = &
                       atomicdata_rho_PBE0_def2QZVPP_b0(1:atomicdata_rho_PBE0_def2QZVPP_maxZ)
            atomicdata_rho_bp(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = &
                       atomicdata_rho_PBE0_def2QZVPP_bp(1:atomicdata_rho_PBE0_def2QZVPP_maxZ)

            atomicdata_rho_am(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = &
                       atomicdata_rho_PBE0_def2QZVPP_am(1:atomicdata_rho_PBE0_def2QZVPP_maxZ)
            atomicdata_rho_a0(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = &
                       atomicdata_rho_PBE0_def2QZVPP_a0(1:atomicdata_rho_PBE0_def2QZVPP_maxZ)
            atomicdata_rho_ap(1:atomicdata_rho_PBE0_def2QZVPP_maxZ) = &
                       atomicdata_rho_PBE0_def2QZVPP_ap(1:atomicdata_rho_PBE0_def2QZVPP_maxZ)
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
    real(DEVDP) :: bm,b0,bp,am,a0,ap,rvdw,aii,bii,rcii
    real(DEVDP) :: ipm,ip0,ipp
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

        bm = 0.0
        b0 = 0.0
        bp = 0.0
        am = 0.0
        a0 = 0.0
        ap = 0.0
        if( (z .ge. 1) .and. (z .le. ATDENS_MAX_Z ) ) then
            bm = atomicdata_rho_bm(z)
            b0 = atomicdata_rho_b0(z)
            bp = atomicdata_rho_bp(z)
            am = atomicdata_rho_am(z)
            a0 = atomicdata_rho_a0(z)
            ap = atomicdata_rho_ap(z)
        end if

        write(DEV_OUT,130) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          am, a0, ap, bm, b0, bp

    end do
    write(DEV_OUT,120)

! ------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,220)

    do i=1,ntypes
        z = types(i)%z

        ipm = 0.0
        ip0 = 0.0
        ipp = 0.0

        if( (z .ge. 1) .and. (z .le. IPEA_MAXZ ) ) then
            ! FIXME - what to do with negative EA
            ipm = abs(atomicdata_ipm(z))
            ip0 = atomicdata_ip0(z)
            ipp = atomicdata_ipp(z)
        end if

        ! FIXME - what to do with negative EA
        bm = 2.0d0 * sqrt(2.0d0 * ipm * DEV_eV2AU) / DEV_AU2A
        b0 = 2.0d0 * sqrt(2.0d0 * ip0 * DEV_eV2AU) / DEV_AU2A
        bp = 2.0d0 * sqrt(2.0d0 * ipp * DEV_eV2AU) / DEV_AU2A

        write(DEV_OUT,230) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          ipm, ip0, ipp, bm, b0, bp
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
            rvdw = vdw_radii(z)
        end if

        write(DEV_OUT,330) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          rvdw

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

        aii = 0.0
        if( bii_source .eq. AD_BII_ATDENS ) then
            aii = ffdev_atomicdata_aii(i)
        end if
        bii = ffdev_atomicdata_bii(i)
        rcii = ffdev_atomicdata_rcii(i,damp_fa,damp_fb)

        write(DEV_OUT,530) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          types(i)%aveq, aii, bii, rcii

    end do
    write(DEV_OUT,520)

100 format('# Bii/Aii from atomic densities: ',A)
110 format('# ID Type  Z  El |   A-     A0     A+   |   B-     B0     B+   |')
120 format('# -- ---- --- -- | ------ ------ ------ | ------ ------ ------ |')
130 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3,1X,F6.3,3X,F6.3,1X,F6.3,1X,F6.3)

200 format('# Bii/Aii from ionization potentials and electron affinities')
210 format('# ID Type  Z  El |   EA    IP0    IP+  -->  B-     B0     B+   |')
220 format('# -- ---- --- -- | ------ ------ ------ | ------ ------ ------ |')
230 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3,1X,F6.3,3X,F6.3,1X,F6.3,1X,F6.3)

300 format('# van der Waals radii')
310 format('# ID Type  Z  El |  Rvdw  |')
320 format('# -- ---- --- -- | ------ |')
330 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3)

410 format('# Bii data source  : ',A)
420 format('# Bii data mods    : ',A)
430 format('# Rcii data source : ',A)

500 format('# Current data: damp_fa = ',F6.3,'; damp_fb= ',F6.3)
510 format('# ID Type  Z  El  Qeff  |   Aii    Bii   Rcii  |')
520 format('# -- ---- --- -- ------ | ------ ------ ------ |')
530 format(I4,1X,A4,1X,I3,1X,A2,1X,F6.3,3X,F6.3,1X,F6.3,1X,F6.3)

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
    real(DEVDP) :: q, bm, b0, bp
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_bii = 0.0 ! 0.0 -> undefined

    select case(bii_source)
        case(AD_BII_IPEA)
            maxz = IPEA_MAXZ
        case(AD_BII_ATDENS)
            maxz = ATDENS_MAX_Z
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
            ! FIXME - what to do with negative EA
            bm = 2.0d0 * sqrt(2.0d0 * abs(atomicdata_ipm(z)) * DEV_eV2AU) / DEV_AU2A
            b0 = 2.0d0 * sqrt(2.0d0 * atomicdata_ip0(z) * DEV_eV2AU) / DEV_AU2A
            bp = 2.0d0 * sqrt(2.0d0 * atomicdata_ipp(z) * DEV_eV2AU) / DEV_AU2A
    ! ------------
        case(AD_BII_ATDENS)
            bm = atomicdata_rho_bm(z)
            b0 = atomicdata_rho_b0(z)
            bp = atomicdata_rho_bp(z)
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
            ! FIXME
            stop
    ! ------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_bii II')
    end select

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
    real(DEVDP) :: q, am, a0, ap
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_aii = 0.0 ! 0.0 -> undefined

    select case(bii_source)
        case(AD_BII_IPEA)
            call ffdev_utils_exit(DEV_ERR,1,'aii is not available from IP/EA in ffdev_atomicdata_aii')
        case(AD_BII_ATDENS)
            maxz = ATDENS_MAX_Z
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_aii 0')
    end select

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. maxz) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_aii')
    end if
    q = types(gti)%aveq

    select case(bii_source)
        case(AD_BII_IPEA)
            call ffdev_utils_exit(DEV_ERR,1,'aii is not available from IP/EA in ffdev_atomicdata_aii')
    ! ------------
        case(AD_BII_ATDENS)
            am = atomicdata_rho_am(z)
            a0 = atomicdata_rho_a0(z)
            ap = atomicdata_rho_ap(z)
    ! ------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_aii I')
    end select

    select case(bii_mods)
        case(AD_BII_RAW)
            ffdev_atomicdata_aii = a0
            return
        case(AD_BII_MOD_BY_XDM)
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported aii modification by XDM in ffdev_atomicdata_aii')
    ! ------------
        case(AD_BII_MOD_BY_CHRG)
            ! FIXME
            stop
    ! ------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_aii II')
    end select


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
        case(AD_RCII_DISP)
            ! workaround - largest Z we have
            maxz = ATDENS_MAX_Z_ALL
    ! -------------
        case(AD_RCII_ATDENS)
            maxz = ATDENS_MAX_Z
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
        case(AD_RCII_DISP)
            if( .not. disp_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1,'no DISP data loaded for rcc in ffdev_atomicdata_rcii')
            end if
            ffdev_atomicdata_rcii = fa * disp_pairs(gti,gti)%rc + fb
    ! -------------
        case(AD_RCII_ATDENS)
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

! ------------------------------------------------------------------------------

real(DEVDP) function ffdev_atomicdata_get_effZ(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    select case(eff_core)
        case(AD_EFF_CORE_NONE)
            ffdev_atomicdata_get_effZ = types(gti)%z
            return
    ! --------------------
        case(AD_EFF_CORE_OPT)
            ffdev_atomicdata_get_effZ = types(gti)%Zeff
            return
    ! --------------------
        case(AD_EFF_CORE_MAX)
            z = types(gti)%z
            ! determine number of valence electrons
            if( z .le. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_get_effZ')
            end if
            if( z .le. 2 ) then
                ! H-He
                ffdev_atomicdata_get_effZ = z
                return
            end if
            if( z .le. 10 ) then
                ! Li-Ne
                ffdev_atomicdata_get_effZ = z - 2
                return
            end if
            if( z .le. 18 ) then
                ! Na-Ar
                ffdev_atomicdata_get_effZ = z - 10
                return
            end if
            ! FIXME
            call ffdev_utils_exit(DEV_ERR,1,'Not defined in ffdev_atomicdata_get_effZ')
    ! --------------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_get_effZ')
    end select

end function ffdev_atomicdata_get_effZ

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
        case('ATDENS')
            ffdev_atomicdata_rcii_source_from_string = AD_RCII_ATDENS
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
        case(AD_RCII_ATDENS)
            ffdev_atomicdata_rcii_source_to_string = 'ATDENS'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_rcii_source_to_string!')
    end select

end function ffdev_atomicdata_rcii_source_to_string

! ==============================================================================
! subroutine ffdev_atomicdata_eff_core_from_string
! ==============================================================================

integer function ffdev_atomicdata_eff_core_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('NONE')
            ffdev_atomicdata_eff_core_from_string = AD_EFF_CORE_NONE
        case('MAX')
            ffdev_atomicdata_eff_core_from_string = AD_EFF_CORE_MAX
        case('OPT')
            ffdev_atomicdata_eff_core_from_string = AD_EFF_CORE_OPT
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atomicdata_eff_core_from_string!')
    end select

end function ffdev_atomicdata_eff_core_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_eff_core_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_eff_core_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AD_EFF_CORE_NONE)
            ffdev_atomicdata_eff_core_to_string = 'NONE'
        case(AD_EFF_CORE_MAX)
            ffdev_atomicdata_eff_core_to_string = 'MAX'
        case(AD_EFF_CORE_OPT)
            ffdev_atomicdata_eff_core_to_string = 'OPT'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_eff_core_to_string!')
    end select

end function ffdev_atomicdata_eff_core_to_string

! ------------------------------------------------------------------------------

end module ffdev_atomicdata

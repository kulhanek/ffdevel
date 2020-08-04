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

    densoverlap_b0ii(:) = 0.0d0
    densoverlap_a0ii(:) = 0.0d0

    select case(atomoverlap_source)
        case(AO_PBE0_def2QZVPP)
            densoverlap_b0ii(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_b0
            densoverlap_a0ii(1:densoverlap_PBE0_def2QZVPP_maxZ) = densoverlap_PBE0_def2QZVPP_a0

            wfoverlap_b0ii(1:densoverlap_PBE0_def2QZVPP_maxZ) = wfoverlap_PBE0_def2QZVPP_b0
            wfoverlap_a0ii(1:densoverlap_PBE0_def2QZVPP_maxZ) = wfoverlap_PBE0_def2QZVPP_a0
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
    use prmfile

    implicit none
    integer :: i,z
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Linearized Electron Density/Wavefunction Overlaps', '=')

    write(DEV_OUT,*)
    write(DEV_OUT,5) trim(ffdev_atomicdata_source_to_string(atomoverlap_source))

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

     do i=1,ntypes
        z = types(i)%z
        write(DEV_OUT,30) i, adjustl(types(i)%name), types(i)%z, adjustl(pt_symbols(types(i)%z)), &
                          densoverlap_b0ii(z), densoverlap_a0ii(z), wfoverlap_b0ii(z), wfoverlap_a0ii(z)
    end do
    write(DEV_OUT,20)

  5 format('# Data source : ',A)

 10 format('# ID Type  Z  El | B0(D0) A0(D0) | B0(WF) A0(WF)')
 20 format('# -- ---- --- -- | ------ ------ | ------ ------')
 30 format(I4,1X,A4,1X,I3,1X,A2,3X,F6.3,1X,F6.3,3X,F6.3,1X,F6.3)

end subroutine ffdev_atomicdata_print

! ==============================================================================
! function ffdev_atomicdata_do_bii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_do_bii(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_xdm_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_do_bii = 1.0 ! default b

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_do_bii')
    end if

    ffdev_atomicdata_do_bii = densoverlap_b0ii(z)

    select case(atomoverlap_mods)
        case(AO_MODS_PLAIN)
            return
        case(AO_MODS_BY_XDM)
            ! XDM mod, DOI: 10.1021/ct1001494
            ffdev_atomicdata_do_bii = ffdev_atomicdata_do_bii * &
                                (xdm_atoms(gti)%v0ave / xdm_atoms(gti)%vave)**(1.0d0/3.0d0)
            return
        case default
            call ffdev_utils_exit(DEV_ERR,1,'No implemented in ffdev_atomicdata_do_bii')
    end select

end function ffdev_atomicdata_do_bii

! ==============================================================================
! function ffdev_atomicdata_do_aii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_do_aii(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_do_aii = 0.0 ! default a

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_do_aii')
    end if

    ffdev_atomicdata_do_aii = densoverlap_a0ii(z)

end function ffdev_atomicdata_do_aii

! ==============================================================================
! function ffdev_atomicdata_wo_bii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_wo_bii(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_xdm_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_wo_bii = 1.0 ! default b

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_wo_bii')
    end if

    ffdev_atomicdata_wo_bii = wfoverlap_b0ii(z)

    select case(atomoverlap_mods)
        case(AO_MODS_PLAIN)
            return
        case(AO_MODS_BY_XDM)
            ! XDM mod, DOI: 10.1021/ct1001494
            ffdev_atomicdata_wo_bii = ffdev_atomicdata_wo_bii * &
                                (xdm_atoms(gti)%v0ave / xdm_atoms(gti)%vave)**(1.0d0/3.0d0)
            return
        case default
            call ffdev_utils_exit(DEV_ERR,1,'No implemented in ffdev_atomicdata_wo_bii')
    end select

end function ffdev_atomicdata_wo_bii

! ==============================================================================
! function ffdev_atomicdata_wo_aii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_wo_aii(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    ! --------------------------------------------
    integer     :: z
    ! --------------------------------------------------------------------------

    ffdev_atomicdata_wo_aii = 0.0 ! default a

    ! get Z
    z = types(gti)%z
    if( (z .le. 0) .and. (z .gt. DENSOVERLAP_MAX_Z) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_wo_aii')
    end if

    ffdev_atomicdata_wo_aii = wfoverlap_a0ii(z)

end function ffdev_atomicdata_wo_aii

! ==============================================================================
! function ffdev_atomicdata_rcii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_rcii(gti,dens)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: gti
    real(DEVDP) :: dens
    ! --------------------------------------------
    real(DEVDP) :: b, a
    ! --------------------------------------------------------------------------

    b = ffdev_atomicdata_do_bii(gti)
    a = ffdev_atomicdata_do_aii(gti)

    if( b .eq. 0.0d0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'No atodens_bii data for given gti in ffdev_atomicdata_rcii')
    end if

    ffdev_atomicdata_rcii = (a - dens)/b

end function ffdev_atomicdata_rcii

! ==============================================================================
! function ffdev_atomicdata_ip_bii
! ==============================================================================

real(DEVDP) function ffdev_atomicdata_ip_bii(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_ip_db
    use ffdev_xdm_dat

    implicit none
    integer         :: gti
    ! --------------------------------------------
    integer         :: z
    ! --------------------------------------------------------------------------

    z = types(gti)%z

    if( (z .gt. ionization_potential_maxZ) .or. (z .le. 0) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_atomicdata_ip_bii!')
    end if

    ! in A^-1 => 1.0 / DEV_AU2A
    ffdev_atomicdata_ip_bii = 2.0d0 * sqrt(2.0d0 * ionization_potential(z) * DEV_eV2AU) / DEV_AU2A

    select case(atomoverlap_mods)
        case(AO_MODS_PLAIN)
            return
        case(AO_MODS_BY_XDM)
            ! XDM mod, DOI: 10.1021/ct1001494
            ffdev_atomicdata_ip_bii = ffdev_atomicdata_ip_bii * &
                                (xdm_atoms(gti)%v0ave / xdm_atoms(gti)%vave)**(1.0d0/3.0d0)
            return
        case default
            call ffdev_utils_exit(DEV_ERR,1,'No implemented in ffdev_atomicdata_ip_bii')
    end select

end function ffdev_atomicdata_ip_bii

! ==============================================================================
! subroutine ffdev_atomicdata_source_from_string
! ==============================================================================

integer function ffdev_atomicdata_source_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('PBE0/def2-QZVPP')
            ffdev_atomicdata_source_from_string = AO_PBE0_def2QZVPP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_atomicdata_source_from_string!')
    end select

end function ffdev_atomicdata_source_from_string

! ==============================================================================
! subroutine ffdev_atomicdata_source_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_source_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(AO_PBE0_def2QZVPP)
            ffdev_atomicdata_source_to_string = 'PBE0/def2-QZVPP'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_source_to_string!')
    end select

end function ffdev_atomicdata_source_to_string

! ------------------------------------------------------------------------------

end module ffdev_atomicdata
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

module ffdev_ip_db

use ffdev_constants
use ffdev_variables

! ==============================================================================
! ionization potential of elements, all in eV
! source: https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)
! ref: CRC

integer         :: ionization_potential_maxZ = 18

real(DEVDP)     :: ionization_potential(1:18) = (/ &
13.59844, 24.58741,  5.39172,  9.32270,  8.29803, 11.26030, 14.53414, 13.61806, 17.42282, 21.5646, &
 5.13908,  7.64624,  5.98577,  8.15169, 10.48669, 10.36001, 12.96764, 15.75962 &
/)

contains

! ==============================================================================
! function ffdev_bfac_from_ip
! ==============================================================================

real(DEVDP) function ffdev_bfac_from_ip(gti)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer         :: gti
    ! --------------------------------------------
    integer         :: z
    ! --------------------------------------------------------------------------

    z = types(gti)%z

    if( (z .gt. ionization_potential_maxZ) .or. (z .le. 0) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_bfac_from_ip!')
    end if

    ! in A^-1 => 1.0 / DEV_AU2A
    ffdev_bfac_from_ip = 2.0d0 * sqrt(2.0d0 * ionization_potential(z) * DEV_eV2AU) / DEV_AU2A

end function ffdev_bfac_from_ip

! ==============================================================================
! function ffdev_bfac_from_ip_xdm
! ==============================================================================

real(DEVDP) function ffdev_bfac_from_ip_xdm(gti)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_xdm_dat

    implicit none
    integer         :: gti
    ! --------------------------------------------
    integer         :: z
    ! --------------------------------------------------------------------------

    z = types(gti)%z

    if( (z .gt. ionization_potential_maxZ) .or. (z .le. 0) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z is out-of-range in ffdev_bfac_from_ip_xdm!')
    end if

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded in ffdev_bfac_from_ip_xdm!')
    end if

    ! in A^-1 => 1.0 / DEV_AU2A
    ffdev_bfac_from_ip_xdm = 2.0d0 * sqrt(2.0d0 * ionization_potential(z) * DEV_eV2AU) / DEV_AU2A
    ! XDM mod, DOI: 10.1021/ct1001494
    ffdev_bfac_from_ip_xdm = ffdev_bfac_from_ip_xdm * &
                        (xdm_atoms(gti)%v0ave / xdm_atoms(gti)%vave)**(1.0d0/3.0d0)

end function ffdev_bfac_from_ip_xdm

! ==============================================================================

end module ffdev_ip_db
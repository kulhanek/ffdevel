! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_errors_control

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_errors_ctrl
! ==============================================================================

subroutine ffdev_errors_ctrl(fin)

    use ffdev_err_bonds_control
    use ffdev_err_angles_control
    use ffdev_err_dihedrals_control
    use ffdev_err_impropers_control
    use ffdev_err_nbdists_control
    use ffdev_err_energy_control
    use ffdev_err_rmsd_control
    use ffdev_err_freqs_control
    use ffdev_err_bonds
    use ffdev_err_angles
    use ffdev_err_dihedrals
    use ffdev_err_impropers
    use ffdev_err_nbdists
    use ffdev_err_energy
    use ffdev_err_rmsd
    use ffdev_err_freqs
    use ffdev_errors_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FFERROR', ':')

    ! clear what should be calculated
    errors_calc_freq = .false.
    errors_calc_hess = .false.
    errors_calc_grad = .false.
    errors_calc_ene  = .false.

    ! reset setup by external request
    if( prmfile_open_section(fin,'setdefault') ) then
        ! reset all setup
        call ffdev_err_energy_init()
        call ffdev_err_bonds_init()
        call ffdev_err_angles_init()
        call ffdev_err_dihedrals_init()
        call ffdev_err_impropers_init()
        call ffdev_err_nbdists_init()
        call ffdev_err_freqs_init()
        call ffdev_err_rmsd_init()
    end if

    ! read setup
    call ffdev_err_energy_ctrl(fin)
    call ffdev_err_bonds_ctrl(fin)
    call ffdev_err_angles_ctrl(fin)
    call ffdev_err_dihedrals_ctrl(fin)
    call ffdev_err_impropers_ctrl(fin)
    call ffdev_err_nbdists_ctrl(fin)
    call ffdev_err_freqs_ctrl(fin)
    call ffdev_err_rmsd_ctrl(fin)

end subroutine ffdev_errors_ctrl

! ------------------------------------------------------------------------------

end module ffdev_errors_control

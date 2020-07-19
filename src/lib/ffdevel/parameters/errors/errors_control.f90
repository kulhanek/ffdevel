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
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_errors_ctrl
! ==============================================================================

subroutine ffdev_errors_ctrl(fin)

    use prmfile
    use ffdev_utils

    use ffdev_errors_dat
    use ffdev_errors

    use ffdev_parameters_dat

    use ffdev_err_bonds_control
    use ffdev_err_angles_control
    use ffdev_err_dihedrals_control
    use ffdev_err_impropers_control
    use ffdev_err_nbdists_control
    use ffdev_err_energy_control
    use ffdev_err_rmsd_control
    use ffdev_err_ihess_control
    use ffdev_err_sapt_control
    use ffdev_err_chrgpnl_control
    use ffdev_err_zerograd_control
    use ffdev_err_probe_control
    use ffdev_err_papnl_control
    use ffdev_err_pbpnl_control
    use ffdev_err_qnb_control
    use ffdev_err_mue_control

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FFERROR', ':')

    ! clear what should be calculated
    errors_calc_hess    = .false.
    errors_calc_grad    = .false.
    errors_calc_ene     = .false.
    errors_calc_sapt   = .false.

    ! reset setup by external request or default setup
    if( ResetAllSetup ) then
        call ffdev_errors_init_all()
        write(DEV_OUT,10)
    else if ( prmfile_open_section(fin,'setdefault') ) then
        call ffdev_errors_init_all()
        write(DEV_OUT,20)
    end if

    ! read setup
    call ffdev_err_energy_ctrl(fin)
    call ffdev_err_sapt_ctrl(fin)
    call ffdev_err_probe_ctrl(fin)
    call ffdev_err_bonds_ctrl(fin)
    call ffdev_err_angles_ctrl(fin)
    call ffdev_err_dihedrals_ctrl(fin)
    call ffdev_err_impropers_ctrl(fin)
    call ffdev_err_nbdists_ctrl(fin)
    call ffdev_err_ihess_ctrl(fin)
    call ffdev_err_rmsd_ctrl(fin)
    call ffdev_err_chrgpnl_ctrl(fin)
    call ffdev_err_zerograd_ctrl(fin)
    call ffdev_err_papnl_ctrl(fin)
    call ffdev_err_pbpnl_ctrl(fin)
    call ffdev_err_qnb_ctrl(fin)
    call ffdev_err_mue_ctrl(fin)

 10 format('>>> INFO: All errors disabled by default (resetallsetup=on)!')
 20 format('>>> INFO: All errors disabled by the explicit request ([setdefault])!')

end subroutine ffdev_errors_ctrl

! ------------------------------------------------------------------------------

end module ffdev_errors_control

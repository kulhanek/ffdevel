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

module ffdev_errors_dat

use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

integer,parameter       :: EE_ABS       = 1 ! absolute
integer,parameter       :: EE_REL       = 2 ! relative
integer,parameter       :: EE_LOG       = 3 ! log
integer,parameter       :: EE_ABSLOG    = 4 ! abs+log

! ------------------------------------------------------------------------------
! error function

type FFERROR_TYPE
    real(DEVDP)         :: total
    real(DEVDP)         :: energy
    real(DEVDP)         :: bonds
    real(DEVDP)         :: angles
    real(DEVDP)         :: dihedrals
    real(DEVDP)         :: impropers
    real(DEVDP)         :: rmsd
    real(DEVDP)         :: nbdists
    real(DEVDP)         :: ihess_bonds
    real(DEVDP)         :: ihess_angles
    real(DEVDP)         :: ihess_dihedrals
    real(DEVDP)         :: ihess_impropers
    real(DEVDP)         :: sapt0_ele
    real(DEVDP)         :: sapt0_rep
    real(DEVDP)         :: sapt0_disp
    real(DEVDP)         :: chrgpnl
end type FFERROR_TYPE

type(FFERROR_TYPE)      :: FFError

! ------------------------------------------------------------------------------
! internal setup
logical                 :: errors_calc_ene      = .false.
logical                 :: errors_calc_sapt0    = .false.
logical                 :: errors_calc_grad     = .false.
logical                 :: errors_calc_hess     = .false.

! ------------------------------------------------------------------------------
integer                 :: ProgramIndex                 = 0

! ------------------------------------------------------------------------------

end module ffdev_errors_dat

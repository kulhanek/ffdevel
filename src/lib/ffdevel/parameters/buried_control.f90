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

module ffdev_buried_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_buried_ctrl
! ==============================================================================

subroutine ffdev_buried_ctrl(fin,exec)

    use ffdev_buried_dat
    use ffdev_buried
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    character(MAX_PATH)         :: buffer
    real(DEVDP)                 :: rnum
    integer                     :: imode
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'buried') ) then
        write(DEV_OUT,115) ffdev_buried_surf_mode_to_string(surface_mode)
        write(DEV_OUT,125) ProbeR
        write(DEV_OUT,135) BuriedExp0
        write(DEV_OUT,145) BuriedBeta
        return
    end if

    if( prmfile_get_string_by_key(fin,'mode', buffer)) then
        imode = ffdev_buried_surf_mode_from_string(buffer)
        write(DEV_OUT,110) ffdev_buried_surf_mode_to_string(imode)
        if( exec ) surface_mode = imode
    else
        write(DEV_OUT,115) ffdev_buried_surf_mode_to_string(surface_mode)
    end if

    if( prmfile_get_real8_by_key(fin,'proberad', rnum)) then
        write(DEV_OUT,120) rnum
        if( exec ) ProbeR = rnum
    else
        write(DEV_OUT,125) ProbeR
    end if

    if( prmfile_get_real8_by_key(fin,'exp0', rnum)) then
        write(DEV_OUT,130) rnum
        if( exec ) BuriedExp0 = rnum
    else
        write(DEV_OUT,135) BuriedExp0
    end if

    if( prmfile_get_real8_by_key(fin,'beta', rnum)) then
        write(DEV_OUT,140) rnum
        if( exec ) BuriedBeta = rnum
    else
        write(DEV_OUT,145) BuriedBeta
    end if

 10 format('=== [buried] ===================================================================')

110  format ('Surface mode (mode)                    = ',a12)
115  format ('Surface mode (mode)                    = ',a12,'                  (default)')

120  format ('Probe radius (proberad)                = ',f21.8)
125  format ('Probe radius (proberad)                = ',f21.8,'         (default)')

130  format ('Weight exp0 (exp0)                     = ',f21.8)
135  format ('Weight exp0 (exp0)                     = ',f21.8,'         (default)')

140  format ('Weight beta (beta)                     = ',f21.8)
145  format ('Weight beta (beta)                     = ',f21.8,'         (default)')

end subroutine ffdev_buried_ctrl

! ------------------------------------------------------------------------------

end module ffdev_buried_control

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

module ffdev_densoverlap_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_densoverlap_ctrl
! ==============================================================================

subroutine ffdev_densoverlap_ctrl(fin,exec)

    use ffdev_utils
    use prmfile
    use ffdev_densoverlap_dat
    use ffdev_densoverlap

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    character(MAX_PATH)         :: buffer
    integer                     :: lsource
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'densoverlap') ) then
        write(DEV_OUT,115) ffdev_densoverlap_source_to_string(densoverlap_source)
        write(DEV_OUT,55) prmfile_onoff(densoverlap_mod_by_charge)
        call ffdev_densoverlap_update_db
        return
    end if

    if( prmfile_get_string_by_key(fin,'source', buffer)) then
        lsource = ffdev_densoverlap_source_from_string(buffer)
        write(DEV_OUT,110) ffdev_densoverlap_source_to_string(lsource)
        if( exec ) densoverlap_source = lsource
    else
        write(DEV_OUT,115) ffdev_densoverlap_source_to_string(densoverlap_source)
    end if

    if( prmfile_get_logical_by_key(fin,'modbycharge', densoverlap_mod_by_charge)) then
        write(DEV_OUT,50) prmfile_onoff(densoverlap_mod_by_charge)
    else
        write(DEV_OUT,55) prmfile_onoff(densoverlap_mod_by_charge)
    end if

    call ffdev_densoverlap_update_db
    call ffdev_densoverlap_print

 10 format('== [densoverlap] =================================================================')

110  format ('Density overlap source (source)        = ',a29)
115  format ('Density overlap source (source)        = ',a29,' (default)')
 50  format ('Modulation by charge (modbycharge)     = ',a12)
 55  format ('Modulation by charge (modbycharge)     = ',a12,'                  (default)')

end subroutine ffdev_densoverlap_ctrl

! ------------------------------------------------------------------------------

end module ffdev_densoverlap_control

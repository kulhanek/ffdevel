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

module ffdev_atomicdata_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_atomicdata_ctrl
! ==============================================================================

subroutine ffdev_atomicdata_ctrl(fin,exec)

    use ffdev_utils
    use prmfile
    use ffdev_atomicdata_dat
    use ffdev_atomicdata

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    character(MAX_PATH)         :: buffer
    integer                     :: lbuff
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'atomicdata') ) then
        write(DEV_OUT,115) ffdev_atomicdata_bii_source_to_string(bii_source)
        write(DEV_OUT,125) ffdev_atomicdata_bii_mods_to_string(bii_mods)
        write(DEV_OUT,135) ffdev_atomicdata_rcii_source_to_string(rcii_source)
        write(DEV_OUT,145) ffdev_atomicdata_effz_mode_to_string(effz_mode)
        call ffdev_atomicdata_update_db
        return
    end if

    if( prmfile_get_string_by_key(fin,'bii_source', buffer)) then
        lbuff = ffdev_atomicdata_bii_source_from_string(buffer)
        write(DEV_OUT,110) ffdev_atomicdata_bii_source_to_string(lbuff)
        if( exec ) bii_source = lbuff
    else
        write(DEV_OUT,115) ffdev_atomicdata_bii_source_to_string(bii_source)
    end if

    if( prmfile_get_string_by_key(fin,'bii_mods', buffer)) then
        lbuff = ffdev_atomicdata_bii_mods_from_string(buffer)
        write(DEV_OUT,120) ffdev_atomicdata_bii_mods_to_string(lbuff)
        if( exec )  bii_mods = lbuff
    else
        write(DEV_OUT,125) ffdev_atomicdata_bii_mods_to_string(bii_mods)
    end if

    if( prmfile_get_string_by_key(fin,'rcii_source', buffer)) then
        lbuff = ffdev_atomicdata_rcii_source_from_string(buffer)
        write(DEV_OUT,130) ffdev_atomicdata_rcii_source_to_string(lbuff)
        if( exec ) rcii_source = lbuff
    else
        write(DEV_OUT,135) ffdev_atomicdata_rcii_source_to_string(rcii_source)
    end if

    if( prmfile_get_string_by_key(fin,'effz_mode', buffer)) then
        lbuff = ffdev_atomicdata_effz_mode_from_string(buffer)
        write(DEV_OUT,140) ffdev_atomicdata_effz_mode_to_string(lbuff)
        if( exec ) effz_mode = lbuff
    else
        write(DEV_OUT,145) ffdev_atomicdata_effz_mode_to_string(effZ_mode)
    end if

    call ffdev_atomicdata_update_db
    if( .not. exec ) then
        call ffdev_atomicdata_print
    end if

 10 format('== [atomicdata] ================================================================')

110  format ('Bii source  (bii_source)               = ',a29)
115  format ('Bii source  (bii_source)               = ',a29,' (default)')
120  format ('Bii modification mode (bii_mods)       = ',a29)
125  format ('Bii modification mode (bii_mods)       = ',a29,' (default)')
130  format ('Rcii source (rcii_source)              = ',a29)
135  format ('Rcii source (rcii_source)              = ',a29,' (default)')
140  format ('Effective charge of nuclei (effz_mode) = ',a29)
145  format ('Effective charge of nuclei (effz_mode) = ',a29,' (default)')

end subroutine ffdev_atomicdata_ctrl

! ------------------------------------------------------------------------------

end module ffdev_atomicdata_control

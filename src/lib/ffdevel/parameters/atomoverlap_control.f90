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

module ffdev_atomoverlap_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_atomoverlap_ctrl
! ==============================================================================

subroutine ffdev_atomoverlap_ctrl(fin,exec)

    use ffdev_utils
    use prmfile
    use ffdev_atomoverlap_dat
    use ffdev_atomoverlap

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    character(MAX_PATH)         :: buffer
    integer                     :: lsource
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'atomoverlap') ) then
        write(DEV_OUT,115) ffdev_atomoverlap_source_to_string(atomoverlap_source)
        call ffdev_atomoverlap_update_db
        return
    end if

    if( prmfile_get_string_by_key(fin,'source', buffer)) then
        lsource = ffdev_atomoverlap_source_from_string(buffer)
        write(DEV_OUT,110) ffdev_atomoverlap_source_to_string(lsource)
        if( exec ) atomoverlap_source = lsource
    else
        write(DEV_OUT,115) ffdev_atomoverlap_source_to_string(atomoverlap_source)
    end if

    call ffdev_atomoverlap_update_db
    call ffdev_atomoverlap_print

 10 format('== [atomoverlap] ===============================================================')

110  format ('Atom overlap source (source)           = ',a29)
115  format ('Atom overlap source (source)           = ',a29,' (default)')

end subroutine ffdev_atomoverlap_ctrl

! ------------------------------------------------------------------------------

end module ffdev_atomoverlap_control

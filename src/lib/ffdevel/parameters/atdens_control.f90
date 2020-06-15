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

module ffdev_atdens_control

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_atdens_ctrl
! ==============================================================================

subroutine ffdev_atdens_ctrl(fin,exec)

    use ffdev_utils
    use prmfile
    use ffdev_atdens_dat
    use ffdev_atdens

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: exec
    character(MAX_PATH)         :: buffer
    integer                     :: lsource
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'atdens') ) then
        write(DEV_OUT,115) ffdev_atdens_source_to_string(atdens_source)
        return
    end if

    if( prmfile_get_string_by_key(fin,'source', buffer)) then
        lsource = ffdev_atdens_source_from_string(buffer)
        write(DEV_OUT,110) ffdev_atdens_source_to_string(lsource)
        if( exec ) atdens_source = lsource
    else
        write(DEV_OUT,115) ffdev_atdens_source_to_string(atdens_source)
    end if

 10 format('=== [atdens] ===================================================================')

110  format ('Atom density source (source)           = ',a12)
115  format ('Atom density source (source)           = ',a12,'                  (default)')

end subroutine ffdev_atdens_ctrl

! ------------------------------------------------------------------------------

end module ffdev_atdens_control

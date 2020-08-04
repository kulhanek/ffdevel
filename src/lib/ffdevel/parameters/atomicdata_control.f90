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
    integer                     :: lsource
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'atomicdata') ) then
        write(DEV_OUT,115) ffdev_atomicdata_source_to_string(atomoverlap_source)
        write(DEV_OUT,125) ffdev_atomicdata_mods_to_string(atomoverlap_mods)
        call ffdev_atomicdata_update_db
        return
    end if

    if( prmfile_get_string_by_key(fin,'source', buffer)) then
        lsource = ffdev_atomicdata_source_from_string(buffer)
        write(DEV_OUT,110) ffdev_atomicdata_source_to_string(lsource)
        if( exec ) atomoverlap_source = lsource
    else
        write(DEV_OUT,115) ffdev_atomicdata_source_to_string(atomoverlap_source)
    end if

    if( prmfile_get_string_by_key(fin,'mods', buffer)) then
        atomoverlap_mods = ffdev_atomicdata_mods_from_string(buffer)
        write(DEV_OUT,120) ffdev_atomicdata_mods_to_string(atomoverlap_mods)
    else
        write(DEV_OUT,125) ffdev_atomicdata_mods_to_string(atomoverlap_mods)
    end if

    call ffdev_atomicdata_update_db
    call ffdev_atomicdata_print

 10 format('== [atomicdata] ================================================================')

110  format ('Atom overlap source (source)           = ',a29)
115  format ('Atom overlap source (source)           = ',a29,' (default)')
120  format ('Bii modification mode (mods)           = ',a29)
125  format ('Bii modification mode (mods)           = ',a29,' (default)')

end subroutine ffdev_atomicdata_ctrl

! ==============================================================================
! subroutine ffdev_err_papnl_control_source_to_string
! ==============================================================================

character(80) function ffdev_atomicdata_mods_to_string(source)

    use ffdev_utils
    use ffdev_atomicdata_dat

    implicit none
    integer  :: source
    ! --------------------------------------------------------------------------

    select case(source)
        case(AO_MODS_PLAIN)
            ffdev_atomicdata_mods_to_string = 'PLAIN - Raw data'
        case(AO_MODS_BY_XDM)
            ffdev_atomicdata_mods_to_string = 'XDM - modified by XDM volumes'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_atomicdata_mods_to_string!')
    end select

end function ffdev_atomicdata_mods_to_string

! ==============================================================================
! subroutine ffdev_atomicdata_mods_from_string
! ==============================================================================

integer function ffdev_atomicdata_mods_from_string(string)

    use ffdev_utils
    use ffdev_atomicdata_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('PLAIN')
            ffdev_atomicdata_mods_from_string = AO_MODS_PLAIN
        case('XDM')
            ffdev_atomicdata_mods_from_string = AO_MODS_BY_XDM
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) &
                                            // '" in ffdev_atomicdata_mods_from_string!')
    end select

end function ffdev_atomicdata_mods_from_string

! ------------------------------------------------------------------------------

end module ffdev_atomicdata_control

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

module ffdev_mmd3_dontrol

use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_mmd3_ctrl
! ==============================================================================

subroutine ffdev_mmd3_ctrl(fin)

    use ffdev_mmd3
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    logical                     :: rst
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'{MMD3}', ':')

    ! open MMD3 group
    if( .not. prmfile_open_group(fin,'MMD3') ) then
        ! no setup -> use default
            ! damping setup
            select case(mmd3_damping)
                case(MMD3_BM)
                    write(DEV_OUT,*)
                    write(DEV_OUT,100)
                    write(DEV_OUT,115) mmd3_bj_a1
                    write(DEV_OUT,125) mmd3_bj_a2
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'mmd3_damping not implemented in ffdev_mmd3_ctrl!')
            end select
        return
    end if

! damping function setup -------------------------------------------------
    select case(MMD3_BM)
        case(MMD3_BM)
            write(DEV_OUT,*)
            write(DEV_OUT,100)

        ! BM dampning setup
            if( prmfile_open_section(fin,'BM') ) then
                if( prmfile_get_real8_by_key(fin,'a1',mmd3_bj_a2) ) then
                    write(DEV_OUT,110) mmd3_bj_a2
                else
                    write(DEV_OUT,115) mmd3_bj_a1
                end if
                if( prmfile_get_real8_by_key(fin,'a2',mmd3_bj_a2) ) then
                    write(DEV_OUT,120) mmd3_bj_a2
                else
                    write(DEV_OUT,125) mmd3_bj_a2
                end if
            else
                write(DEV_OUT,115) mmd3_bj_a1
                write(DEV_OUT,125) mmd3_bj_a2
            end if
        case default
            call ffdev_utils_exit(DEV_OUT,1,'mmd3_damping not implemented in ffdev_mmd3_ctrl!')
    end select

100 format('=== [BM] =======================================================================')
110 format('BJ damping parameter a1 (a1) = ',F10.6)
115 format('BJ damping parameter a1 (a1) = ',F10.6,' (default)')
120 format('BJ damping parameter a2 (a2) = ',F10.6)
125 format('BJ damping parameter a2 (a2) = ',F10.6,' (default)')

end subroutine ffdev_mmd3_ctrl

! ==============================================================================

end module ffdev_mmd3_dontrol

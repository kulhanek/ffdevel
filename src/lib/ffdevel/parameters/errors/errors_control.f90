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

    use ffdev_errors_dat
    use ffdev_utils
    use prmfile

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    character(MAX_PATH) :: errfcename
    logical             :: retval
    integer             :: alloc_stat, i
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FFERROR', ':')

    ! by default all setup is reset if this part is reached
    if( NumOfErrorFces .gt. 0 ) then
        do i=1,NumOfErrorFces
            deallocate(ErrorFceList(i)%ErrFce)
        end do
        deallocate(ErrorFceList)
        NumOfErrorFces = 0
    end if

    ! load new error setup
    ! count number of sections in the group
    NumOfErrorFces = prmfile_count_group(fin)
    if( NumOfErrorFces .le. 0 ) return

    ! allocate
    allocate(ErrorFceList(NumOfErrorFces), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for error functions!')
    end if

    retval = prmfile_first_section(fin)
    i = 1
    do while( retval )
        ! create error function
        retval = prmfile_get_section_name(fin,errfcename)

        write(DEV_OUT,*)
        write(DEV_OUT,10) trim(errfcename)

        call ffdev_errors_ctrl_create_errfce(ErrorFceList(i)%ErrFce,errfcename)

        ! setup error function
        call ErrorFceList(i)%ErrFce%init_errfce()
        call ErrorFceList(i)%ErrFce%load_errfce(fin)
        call ErrorFceList(i)%ErrFce%set_title_errfce()

        ! next
        retval = prmfile_next_section(fin)
        i = i + 1
    end do

10 format('# === [',A,'] ===')

end subroutine ffdev_errors_ctrl

! ==============================================================================
! subroutine ffdev_errors_ctrl_create_errfce
! ==============================================================================

subroutine ffdev_errors_ctrl_create_errfce(errfce,errfcename)

    use ffdev_errors_dat
    use ffdev_utils

! geometry based
    use ffdev_err_bonds
    use ffdev_err_angles
    use ffdev_err_dihedrals
    use ffdev_err_impropers
    use ffdev_err_nbdists
    use ffdev_err_rmsd

! energy based
    use ffdev_err_energy
    use ffdev_err_ihess
    use ffdev_err_mue
    use ffdev_err_zerograd

! parameters
    use ffdev_err_l1reg
    use ffdev_err_l2reg
    use ffdev_err_bond_r0
    use ffdev_err_angle_a0

    implicit none
    class(ErrorFceType),pointer :: errfce
    character(MAX_PATH)         :: errfcename
    ! --------------------------------------------
    integer             :: alloc_stat
    ! --------------------------------------------------------------------------

    select case(trim(errfcename))
    ! geometry based
        case('bonds')
            allocate(TypeEFTBonds::errfce, stat = alloc_stat)
        case('angles')
            allocate(TypeEFTAngles::errfce, stat = alloc_stat)
        case('dihedrals')
            allocate(TypeEFTDihedrals::errfce, stat = alloc_stat)
        case('impropers')
            allocate(TypeEFTImpropers::errfce, stat = alloc_stat)
        case('nbdists')
            allocate(TypeEFTImpropers::errfce, stat = alloc_stat)
        case('rmsd')
            allocate(TypeEFTImpropers::errfce, stat = alloc_stat)

    ! energy based
        case('energy')
            allocate(TypeEFTEnergy::errfce, stat = alloc_stat)
        case('ihess')
            allocate(TypeEFTIhess::errfce, stat = alloc_stat)
        case('mue')
            allocate(TypeEFTMUE::errfce, stat = alloc_stat)
        case('zerogrd')
            allocate(TypeEFTZeroGrad::errfce, stat = alloc_stat)

    ! parameters
        case('l1reg')
            allocate(TypeEFTL1Reg::errfce, stat = alloc_stat)
        case('l2reg')
            allocate(TypeEFTL2Reg::errfce, stat = alloc_stat)
        case('bond_r0')
            allocate(TypeEFTBondR0::errfce, stat = alloc_stat)
        case('angle_a0')
            allocate(TypeEFTAngleA0::errfce, stat = alloc_stat)

    ! not found
        case default
            call ffdev_utils_exit(DEV_ERR,1, &
                       'The error function '''//trim(errfcename)//''' is not implemented!')
    end select

    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1, &
                   'Unable to allocate memory for the error function: '''//trim(errfcename)//'''!')
    end if

end subroutine ffdev_errors_ctrl_create_errfce

! ------------------------------------------------------------------------------

end module ffdev_errors_control

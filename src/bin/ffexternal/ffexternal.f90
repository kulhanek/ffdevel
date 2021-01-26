! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
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

program ffdev_external_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use ffdev_gradient
    use ffdev_gradient_utils
    use ffdev_hessian
    use ffdev_hessian_utils
    use ffdev_timers

    implicit none

    character(len=MAX_PATH) :: topname      ! input topology name
    character               :: layer        ! requested layer
    character(len=MAX_PATH) :: inpname      ! input data
    character(len=MAX_PATH) :: outname      ! output data
    character(len=MAX_PATH) :: arg
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    integer                 :: mode,i         ! what should be calculated
    logical                 :: do_numerical_g
    logical                 :: do_numerical_h
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('FF External')

    call ffdev_timers_init_top
    call ffdev_timers_init

    ! test number of input arguments
    if( command_argument_count() .lt. 5  ) then
        call print_usage()
        write(DEV_OUT,*) 'Number of arguments: ',command_argument_count()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (4 or more are required)!')
    end if

    do_numerical_g = .false.
    do_numerical_h = .false.

    i = 1
    do while (i .le. command_argument_count() )
        call get_command_argument(i, arg)
        select case(trim(arg))
            case('-p')
                i = i + 1
                if( i .le. command_argument_count() ) call get_command_argument(i, topname)
                i = i + 1
            case('-ng')
                do_numerical_g = .true.
                i = i + 1
            case('-nh')
                do_numerical_h = .true.
                i = i + 1
            case default
                exit
        end select
    end do

    if( i .lt. 3  ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified after topology!')
    end if

    call get_command_argument(i+0, layer)
    call get_command_argument(i+1, inpname)
    call get_command_argument(i+2, outname)

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(topname)
    write(DEV_OUT,110) trim(layer)
    write(DEV_OUT,120) trim(inpname)
    write(DEV_OUT,130) trim(outname)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,topname)
    call ffdev_topology_info(top)

    ! load coordinates
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input Data','=')

    call ffdev_geometry_init(geo)
    call ffdev_geometry_load_ginp(geo,inpname,mode)
    call ffdev_geometry_info_input(geo)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo)

    ! finalize topology and geometry setup
    call ffdev_topology_finalize_setup(top)
    select case(mode)
        case(0)
            write(DEV_OUT,140)
            ! calculate energy and gradient
            call ffdev_energy_all(top,geo)
        case(1)
            write(DEV_OUT,141)
            call ffdev_gradient_allocate(geo)
            if( do_numerical_g ) then
                call ffdev_gradient_num_all(top,geo)
            else
                call ffdev_gradient_all(top,geo)
            end if
        case(2)
            write(DEV_OUT,142)
            call ffdev_gradient_allocate(geo)
            call ffdev_hessian_allocate(geo)
            if( do_numerical_h ) then
                call ffdev_hessian_num_by_grds_all(top,geo)
            else
                call ffdev_hessian_all(top,geo)
            end if
    end select

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Energies','=')
    call ffdev_geometry_info_ene(geo)

    call ffdev_geometry_save_gout(geo,outname,mode)

    call ffdev_timers_finalize(.true.)

    call ffdev_utils_footer('FF External')

100 format('Simplified topology : ',A)
110 format('Requested layer     : ',A)
120 format('Input data          : ',A)
130 format('Output data         : ',A)
140 format('Requested mode         = energy')
141 format('Requested mode         = gradient')
142 format('Requested mode         = hessian')

contains

!===============================================================================
! subroutine:  print_usage
!===============================================================================

subroutine print_usage()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(DEV_OUT,10)
    write(DEV_OUT,*)

    return

10 format('    ffexternal -p <stopology> [-ng] [-nh] Layer InputFile OutputFile MsgFile FChkFile MatElFile')

end subroutine print_usage

!===============================================================================

end program ffdev_external_program

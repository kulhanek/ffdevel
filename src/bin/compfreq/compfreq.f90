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

program ffdev_compfreq_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_topology
    use ffdev_geometry
    use ffdev_hessian
    use ffdev_hessian_utils

    implicit none
    character(len=MAX_PATH) :: crdname1     ! input coordinate name #1
    character(len=MAX_PATH) :: crdname2     ! input coordinate name #2
    character(len=MAX_PATH) :: arg
    type(GEOMETRY)          :: geo1
    type(GEOMETRY)          :: geo2
    double precision        :: rmsd    
    integer                 :: alloc_stat,i
    logical                 :: by_modes
    logical                 :: by_geo
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Compare Frequencies')

    ! test number of input arguments
    if( command_argument_count() .lt. 2 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (at least two expected)!')
    end if

    call get_command_argument(1, crdname1)
    call get_command_argument(2, crdname2)

    by_geo = .false.
    by_modes = .false.

    do i=3,command_argument_count()
        call get_command_argument(i, arg)
        select case(trim(arg))
            case('bygeo')
                by_geo = .true.
            case('bymodes')
                by_modes = .true.
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unrecognized argument ('//trim(arg)//')')
        end select
    end do

    ! by_geo is default comparison mode
    if( (by_modes .eqv. .false.) .and. (by_geo .eqv. .false.) ) then
        by_geo = .true.
    end if

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,110) trim(crdname1)
    write(DEV_OUT,120) trim(crdname2)

    ! load coordinates #1 ------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input coordinates #1','=')

    call ffdev_geometry_init(geo1)
    call ffdev_geometry_load_point(geo1,crdname1)
    call ffdev_geometry_info_input(geo1)

    if( .not. geo1%trg_hess_loaded ) then
        call ffdev_utils_exit(DEV_OUT,1,'Point #1 does not contain Hessian!')
    end if

    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(geo1)

    ! calculate frequencies
    call ffdev_hessian_allocate(geo1)
    call ffdev_hessian_allocate_freq(geo1)
    geo1%hess = geo1%trg_hess
    call ffdev_hessian_calc_freqs(geo1)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Frequencies #1', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo1)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo1)

    ! load coordinates #2 ------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input coordinates #2','=')

    call ffdev_geometry_init(geo2)
    call ffdev_geometry_load_point(geo2,crdname2)
    call ffdev_geometry_info_input(geo2)

    if( .not. geo2%trg_hess_loaded ) then
        call ffdev_utils_exit(DEV_OUT,1,'Point #2 does not contain Hessian!')
    end if

    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(geo2)

    ! calculate frequencies
    call ffdev_hessian_allocate(geo2)
    call ffdev_hessian_allocate_freq(geo2)
    geo2%hess = geo2%trg_hess
    call ffdev_hessian_calc_freqs(geo2)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Frequencies #2', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo2)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo2)

    ! final check --------------------------------------------------------------
    if( geo1%natoms .ne. geo2%natoms ) then
        call ffdev_utils_exit(DEV_OUT,1,'Points do not have the same number of atoms!')
    end if

    do i=1,geo1%natoms
        if( geo1%z(i) .ne. geo2%z(i) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Points do not have the same order/type of atoms!')
        end if
    end do

    allocate(geo1%freq_t2s_map(3*geo1%natoms), geo1%freq_t2s_angles(3*geo1%natoms), &
             geo1%freq_t2s_rmsd(3*geo1%natoms), stat = alloc_stat)

    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate array for maps!')
    end if


! by_geo comparison
    if( by_geo ) then

        ! compare freqs items
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'COMPARISON BY GEO','#')

        ! superimpose structures
        call ffdev_hessian_superimpose_freqs_by_geo(geo1%natoms,geo1%z,geo1%crd,geo2%crd,geo2%nmodes,.false.,rmsd)
        geo1%freq_t2s_rmsd(:) = rmsd

        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Superimposed coordinates #2', ':')
        write(DEV_OUT,130) rmsd
        write(DEV_OUT,*)
        call ffdev_geometry_print_xyz(geo2)

        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Superimposed vectors #2', ':')
        write(DEV_OUT,*)
        call ffdev_hessian_print_freqs(DEV_OUT,geo2)

        ! find vector mappings
        call ffdev_hessian_find_mapping_by_angle(geo1%natoms,geo1%nmodes,geo2%nmodes,geo1%freq_t2s_map,geo1%freq_t2s_angles)

        ! compare freqs items
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'COMPARISON','=')

        call ffdev_hessian_print_mapping(.true.,geo1%natoms,geo1%freq,geo1%nmodes,geo2%freq,geo2%nmodes,&
                                         geo1%freq_t2s_map,geo1%freq_t2s_angles,geo1%freq_t2s_rmsd)
    end if

! by_modes comparison
    if( by_modes ) then

        ! compare freqs items
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'COMPARISON BY MODES','#')

        ! superimpose and calculate mapping
        do i=7,3*geo1%natoms
            ! superimpose
            call ffdev_hessian_superimpose_freqs_by_mode(i,geo1%natoms,geo1%crd,geo2%crd, &
                                                 geo1%nmodes,geo2%nmodes,rmsd)
            geo1%freq_t2s_rmsd(i) = rmsd

            ! and then find mapping
            call ffdev_hessian_find_mapping_by_mode(i,geo1%natoms,geo1%nmodes, &
                                            geo2%nmodes,geo1%freq_t2s_map, &
                                            geo1%freq_t2s_angles)
        end do

        ! compare freqs items
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'COMPARISON','=')

        call ffdev_hessian_print_mapping(.true.,geo1%natoms,geo1%freq,geo1%nmodes,geo2%freq,geo2%nmodes,&
                                         geo1%freq_t2s_map,geo1%freq_t2s_angles,geo1%freq_t2s_rmsd)

    end if

    call ffdev_utils_footer('Compare Frequencies')

110 format('Input coordinates #1 : ',A)
120 format('Input coordinates #2 : ',A)
130 format('RMSD(#2->#1) = ',F12.6)

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

10 format('    compfreq <pst#1> <pst2#2> [bymodes] [bygeo]')

end subroutine print_usage

!===============================================================================

end program ffdev_compfreq_program

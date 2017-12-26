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

program ffdev_enefilter_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use prmfile
    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use smf_xyzfile
    use smf_xyzfile_type

    implicit none
    character(len=MAX_PATH)     :: topname
    character(len=MAX_PATH)     :: itrjname
    character(len=MAX_PATH)     :: otrjname
    logical                     :: rst
    character(PRMFILE_MAX_PATH) :: string
    integer                     :: i, probe_size, snap
    real(DEVDP)                 :: enetreshold, ene
    type(TOPOLOGY)              :: top
    type(GEOMETRY)              :: geo
    type(XYZFILE_TYPE)          :: fin
    type(XYZFILE_TYPE)          :: fout
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('vdWEneFilter')

    ! test number of input arguments
    if( command_argument_count() .ne. 5 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (five expected)!')
    end if

    ! input data
    write(DEV_OUT,*)
    call get_command_argument(1, topname)
    write(DEV_OUT,100) trim(topname)

    call get_command_argument(2, string)
    read(string,*,end=200,err=200) probe_size
    write(DEV_OUT,110) probe_size

    call get_command_argument(3, string)
    read(string,*,end=210,err=210) enetreshold
    write(DEV_OUT,120) enetreshold

    call get_command_argument(4, itrjname)
    write(DEV_OUT,130) trim(itrjname)

    call get_command_argument(5, otrjname)
    write(DEV_OUT,140) trim(otrjname)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,trim(topname))

    ! switch to probe mode and finalizae topology
    top%probe_size = probe_size
    call ffdev_topology_finalize_setup(top)
    call ffdev_topology_info(top)

    ! open input/output files
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Filter Snapshosts','=')

    call init_xyz(fin)
    call open_xyz(DEV_TRAJ,itrjname,fin,'OLD')

    call init_xyz(fout)
    call allocate_xyz(fout,top%natoms)
    call open_xyz(DEV_OTRAJ,otrjname,fout,'UNKNOWN')

    write(DEV_OUT,*)
    write(DEV_OUT,150)
    write(DEV_OUT,160)

    ! read xyz stream
    call ffdev_geometry_init(geo)
    snap = 1
    do while (is_next_xyz_record(DEV_TRAJ,fin))
        ! load snapshot and perform integrity check
        call read_xyz(DEV_TRAJ,fin)
        call ffdev_geometry_load_xyz_snapshot(geo,fin)
        call ffdev_geometry_check_z(top,geo)
        ! calculate energy
        call ffdev_energy_all(top,geo)
        ene = geo%total_ene
        if( ene .le. enetreshold ) then
            write(DEV_OUT,170) snap, ene, 'OK'
            ! save snapshot
            call ffdev_geometry_save_xyz_snapshot(geo,fout)
            call write_xyz(DEV_OTRAJ,fout)
        else
            write(DEV_OUT,170) snap, ene, ' OUT-OF-RANGE'
            ! skip snapshot
        end if
        snap = snap + 1
    end do

    ! close all streams
    call close_xyz(DEV_TRAJ,fin)
    call close_xyz(DEV_OTRAJ,fout)

    call ffdev_utils_footer('vdWEneFilter')

    stop

100 format('Simplified topology   : ',A)
110 format('Probe size (natoms)   : ',I6)
120 format('Energy treshold       : ',F10.6,' [kcal/mol]')
130 format('Input XYZ trajectory  : ',A)
140 format('Output XYZ trajectory : ',A)
150 format('# snapshot   energy   flag')
160 format('# -------- ---------- ----')
170 format(2X,I8,1X,F10.4,1X,A)

200 call ffdev_utils_exit(DEV_OUT,1,'Probe size must be an integer number!')
210 call ffdev_utils_exit(DEV_OUT,1,'Energy treshold must be an real number!')

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

10 format('    vdwenefilter <stop> <probesize> <enetreshold> <inxyz> <outxyz>')

end subroutine print_usage

!===============================================================================

end program ffdev_enefilter_program

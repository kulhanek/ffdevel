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

program ffdev_genpoints_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables
    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use prmfile
    use ffdev_genpoints_control
    use ffdev_genpoints_dat
    use smf_xyzfile
    use smf_xyzfile_type
    use smf_periodic_table_dat
    use ffdev_geoopt_control

    implicit none
    character(len=MAX_PATH) :: ctrlname      ! input control file name
    type(PRMFILE_TYPE)      :: fin
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    type(XYZFILE_TYPE)      :: traj
    integer                 :: alloc_stat, i, j
    real(DEVDP)             :: torang
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Generate Points')

    ! check if control file was provided
    if( command_argument_count() .ne. 1 ) then
        call print_usage
        call ffdev_utils_exit(DEV_ERR,1,'No input file specified on the command line!')
    end if

    call get_command_argument(1, ctrlname)

    ! process control file -----------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Reading control file', ':')
    write(DEV_OUT,'(a,a)') 'Control file name : ',trim(ctrlname)

    call prmfile_init(fin)

    if( .not. prmfile_read(fin,ctrlname) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Specified control file cannot be opened!')
    end if

    if( .not. prmfile_open_group(fin,'MAIN') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Specified control file does not contain {MAIN} group!')
    end if

    call ffdev_genpoints_ctrl_files(fin)
    call ffdev_genpoints_ctrl_points(fin)

    if( OptimizePoints ) then
        call ffdev_geoopt_ctrl_minimize(fin)
    end if

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(fin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_ERR,1,'Unprocessed lines found in the control file!')
    end if

    ! release the file
    call prmfile_clear(fin)

    ! --------------------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input data', ':')

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,GenTopName)
    call ffdev_topology_info(top)

    ! load coordinates
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates','=')

    call ffdev_geometry_init(geo)
    call ffdev_geometry_load_point(geo,GenCrdName)
    call ffdev_geometry_info_input(geo)
    write(DEV_OUT,*)
    call ffdev_geometry_info_point_header()
    call ffdev_geometry_info_point(geo)

    ! check coordinates and topology
    if( top%natoms .ne. geo%natoms ) then
        call ffdev_utils_exit(DEV_ERR,1,'Topology and coordinates are not compatible (different number of atoms)!')
    end if

    ! finalize topology and geometry setup
    call ffdev_topology_finalize_setup(top)

    if( RotorListLoaded ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Rotors','=')
        call read_rotors(top)
    end if

    ! --------------------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Generate points', ':')

    ! fix some dimensions
    select case(GeneratorMethod)
        case(GENPOINTS_NMODES_METHOD)
            if( top%natoms .le. 2 ) then
                call ffdev_utils_exit(DEV_ERR,1,'More than two atoms are required for the nmodes generator!')
            end if
            NPoints = (3*top%natoms - 6)*2*NModePoints
            if( MaxPoints .lt. NPoints ) then
                MaxPoints = NPoints
                write(DEV_OUT,120) MaxPoints
            end if
    end select

    NPoints = MaxPoints
    allocate( Points(NPoints), AtomMask(top%natoms), ProcessingStack(top%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
         call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate arrays for points and atom mask!')
    end if

    select case(GeneratorMethod)
        case(GENPOINTS_SYSTEMATIC_METHOD)
            call genpoints_systematic
        case(GENPOINTS_STEPBYSTEP_METHOD)
            call genpoints_stepbystep
        case(GENPOINTS_STOCHASTIC_METHOD)
            call genpoints_stochastic
        case(GENPOINTS_STOCHASTIC_OPT_METHOD)
            call genpoints_stochastic_opt
        case(GENPOINTS_STOCHASTIC_BY_STEPS_METHOD)
            call genpoints_stochastic_by_steps
        case(GENPOINTS_NMODES_METHOD)
            call genpoints_nmodes
    end select

    ! --------------------------------------------------------------------------
    if( OptimizePoints ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Optimize points', ':')
        call genpoints_optimize_points
    end if

    ! --------------------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Finalize', ':')

    write(DEV_OUT,100) trim(GenOutName)
    write(DEV_OUT,105) CurrPtsIndex

    call init_xyz(traj)
    call allocate_xyz(traj,top%natoms)
    do i=1,top%natoms
        traj%symbols(i) = pt_symbols(top%atom_types(top%atoms(i)%typeid)%z)
    end do

    ! ---------------------------
    call open_xyz(DEV_TRAJ,GenOutName,traj,'UNKNOWN')
    do i=1,CurrPtsIndex
        write(traj%comment,'(F20.6)') Points(i)%energy
        traj%cvs = Points(i)%crd
        call write_xyz(DEV_TRAJ,traj)
    end do
    call close_xyz(DEV_TRAJ,traj)

    ! ---------------------------
    if( (MinEnergyPtsIndex .gt. 0) .or. (MinOptEnergyPtsIndex .gt. 0) ) then
        write(DEV_OUT,110) trim(GenFinName)
        call open_xyz(DEV_TRAJ,GenFinName,traj,'UNKNOWN')
        if( OptimizePoints ) then
            write(traj%comment,'(F20.6)') Points(MinOptEnergyPtsIndex)%energy
            traj%cvs = Points(MinOptEnergyPtsIndex)%crd
        else
            write(traj%comment,'(F20.6)') Points(MinEnergyPtsIndex)%energy
            traj%cvs = Points(MinEnergyPtsIndex)%crd
        end if
        call write_xyz(DEV_TRAJ,traj)
        call close_xyz(DEV_TRAJ,traj)

        call free_xyz(traj)
    end if

    ! ---------------------------
    ! profile
    write(DEV_OUT,200) trim(ProfileName)
    call ffdev_utils_open(DEV_PROFILE,ProfileName,'U')

    do i=1,CurrPtsIndex
        write(DEV_PROFILE,205,ADVANCE='NO') i

        ! calculate rotor angle for given point
        do j=1,NRotors
            torang = ffdev_geometry_get_dihedral(Points(i)%crd,Rotors(j)%ait,Rotors(j)%ai,Rotors(j)%aj,Rotors(j)%ajt)
            write(DEV_PROFILE,210,ADVANCE='NO') torang*DEV_R2D
        end do

        ! write energy
        write(DEV_PROFILE,215) Points(i)%energy, Points(i)%energy-Points(MinEnergyPtsIndex)%energy
    end do

    close(DEV_PROFILE)

    ! --------------------------------------------------------------------------

    call ffdev_utils_footer('Generate Points')

100 format('File name with generated points     = ',A)
105 format('Number of generated points          = ',I6)
110 format('File name with global minimum point = ',A)
120 format('>>> INFO: Increasing maximum number of points to : ',I6)

200 format('File name with energy profile       = ',A)
205 format(I6)
210 format(1X,F10.4)
215 format(1X,F10.4,1X,F10.4)

contains

!===============================================================================
! subroutine print_usage
!===============================================================================

subroutine print_usage()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(DEV_OUT,10)
    write(DEV_OUT,*)

    return

10 format('    genpoints <ctrlinput>')

end subroutine print_usage

!===============================================================================
! subroutine read_rotors(top)
!===============================================================================

subroutine read_rotors(top)

    implicit none
    type(TOPOLOGY)      :: top
    ! --------------------------------------------
    type(PRMFILE_TYPE)  :: rotfin
    integer             :: i, ai, aj, alloc_stat, mb, j
    logical             :: first
    ! --------------------------------------------------------------------------

    call prmfile_init(rotfin)

    if( .not. prmfile_read(rotfin,GenRotName) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Specified rotor bond file cannot be opened!')
    end if

    if( .not. prmfile_open_group(rotfin,'MAIN') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Specified rotor bond file does not contain {MAIN} group!')
    end if

    if( .not. prmfile_open_section(rotfin,'rotors') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [rotors] section!')
    end if

    ! number of lines
    nrotors = prmfile_count_section(rotfin)

    if( nrotors .le. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'At least one rotor must be specified in the rotor bond file!')
    end if

    ! allocate data
    allocate( rotors(nrotors), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
         call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate list for rotors!')
    end if

    ! read individual rotors
    i = 1
    do while( prmfile_get_int_int(rotfin,ai,aj) )
        if( i .gt. nrotors ) then
            call ffdev_utils_exit(DEV_ERR,1,'Mismatch in the rotor bond file parser!')
        end if
        if( (ai .le. 0) .or. (ai .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'ai out-of-limits!')
        end if
        if( (aj .le. 0) .or. (aj .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'aj out-of-limits!')
        end if
        rotors(i)%ai = ai
        rotors(i)%aj = aj
        i = i + 1
    end do

    ! check if everything was read
    if( prmfile_count_ulines(rotfin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(rotfin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_ERR,1,'Unprocessed lines found in the rotor bond file!')
    end if

    ! release the file
    call prmfile_clear(rotfin)


    ! generate terminals
    do i=1,nrotors
        first = .true.
        mb = 0
        do j = 1,top%atoms(rotors(i)%ai)%nbonds
            if( top%atoms(rotors(i)%ai)%bonded(j) .eq. rotors(i)%aj ) cycle ! skip
            if( first .eqv. .true. ) then
                rotors(i)%ait = top%atoms(rotors(i)%ai)%bonded(j)
                mb = top%atoms( top%atoms(rotors(i)%ai)%bonded(j) )%nbonds
            end if
            if( mb .le. top%atoms( top%atoms(rotors(i)%ai)%bonded(j) )%nbonds ) then
                rotors(i)%ait = top%atoms(rotors(i)%ai)%bonded(j)
                mb = top%atoms( top%atoms(rotors(i)%ai)%bonded(j) )%nbonds
            end if
        end do

        first = .true.
        mb = 0
        do j = 1,top%atoms(rotors(i)%aj)%nbonds
            if( top%atoms(rotors(i)%aj)%bonded(j) .eq. rotors(i)%ai ) cycle ! skip
            if( first .eqv. .true. ) then
                rotors(i)%ajt = top%atoms(rotors(i)%aj)%bonded(j)
                mb = top%atoms( top%atoms(rotors(i)%aj)%bonded(j) )%nbonds
            end if
            if( mb .le. top%atoms( top%atoms(rotors(i)%aj)%bonded(j) )%nbonds ) then
                rotors(i)%ajt = top%atoms(rotors(i)%aj)%bonded(j)
                mb = top%atoms( top%atoms(rotors(i)%aj)%bonded(j) )%nbonds
            end if
        end do
    end do

    write(DEV_OUT,100) nrotors

    write(DEV_OUT,*)
    write(DEV_OUT,110)
    write(DEV_OUT,120)

    do i=1,nrotors
        write(DEV_OUT,130) rotors(i)%ait, top%atoms(rotors(i)%ait)%name, &
                           top%atoms(rotors(i)%ait)%residx,top%atoms(rotors(i)%ait)%resname, &
                           rotors(i)%ai, top%atoms(rotors(i)%ai)%name, &
                           top%atoms(rotors(i)%ai)%residx,top%atoms(rotors(i)%ai)%resname, &
                           rotors(i)%aj, top%atoms(rotors(i)%aj)%name, &
                           top%atoms(rotors(i)%aj)%residx,top%atoms(rotors(i)%aj)%resname, &
                           rotors(i)%ajt, top%atoms(rotors(i)%ajt)%name, &
                           top%atoms(rotors(i)%ajt)%residx,top%atoms(rotors(i)%ajt)%resname

    end do

100 format('Number of rotor bonds = ',I6)
110 format('# Indx Name RIndx RNam    Indx  Name RIndx RNam    Indx  Name RIndx RNam    Indx  Name RIndx RNam')
120 format('# ---- ---- ----- ---- | ------ ---- ----- ---- ~ ------ ---- ----- ---- | ------ ---- ----- ----')
130 format(I6,1X,A4,1X,I5,1X,A4,3X,I6,1X,A4,1X,I5,1X,A4,3X,I6,1X,A4,1X,I5,1X,A4,3X,I6,1X,A4,1X,I5,1X,A4)

end subroutine read_rotors

!===============================================================================
! subroutine genpoints_stepbystep
!===============================================================================

subroutine genpoints_stepbystep

    use ffdev_genpoints_utils

    implicit none
    integer         :: est_npoints, i, itmp, lticks
    integer         :: nticks, alloc_status, nruns, counter, rotor
    type(GEOMETRY)  :: tmpgeo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(A)') 'Method = stepbystep'

    ! energy of the first structure
    call ffdev_energy_all(top,geo)
    MinEnergy = geo%total_ene
    MinEnergyPtsIndex = 1

    nticks = INT(2.0d0*DEV_PI / TickAngle)
    est_npoints = nticks * nrotors ! only one * !!!!

    write(DEV_OUT,*)
    write(DEV_OUT,10) MinEnergy
    write(DEV_OUT,20) est_npoints
    write(DEV_OUT,30) MaxPoints
    write(DEV_OUT,35) nticks

    ! add the first point to stack
    CurrPtsIndex = 1
    Points(CurrPtsIndex)%energy = MinEnergy
    allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
             stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
    end if
    Points(CurrPtsIndex)%crd(:,:) = geo%crd(:,:)
    Points(CurrPtsIndex)%angles(:) = 0.0d0
    Points(CurrPtsIndex)%active_angle = 0

    ! init temporary geometry
    call ffdev_geometry_init(tmpgeo)

    RejectedPoints = 0
    nruns = 1

    write(DEV_OUT,*)
    write(DEV_OUT,'(A)',advance='NO') 'Progress: |'
    call flush(DEV_OUT)

    counter = MaxPoints / 68 + 1

    ! generate points
    do while( (CurrPtsIndex .lt. MaxPoints) .and. (nruns .lt. est_npoints) )

        ! copy initial structure
        call ffdev_geometry_copy(tmpgeo,geo)

        ! generate rotation vector
        itmp = nruns
        do i=1,NRotors
            Rotors(i)%angle = 0.0d0
        end do

        lticks = mod(nruns,nticks)
!        if( lticks .gt. nticks/2 ) then
!            lticks = nticks/2 - lticks
!        end if
        rotor   = nruns / nticks + 1
        if( rotor .le. NRotors ) then
            Rotors(rotor)%angle = lticks * TickAngle
        end if

        ! generate point
        call genpoints_genrot(top,tmpgeo)

        ! get energy
        call ffdev_energy_all(top,tmpgeo)

        nruns = nruns + 1

        ! is energy acceptable
        if( (tmpgeo%total_ene - MinEnergy) .gt. MaxEnergy ) then
            RejectedPoints = RejectedPoints + 1
            cycle
        end if

        CurrPtsIndex = CurrPtsIndex + 1

        if( mod(CurrPtsIndex,counter) .eq. 0 ) then
            write(DEV_OUT,'(A)',advance='NO') '*'
            call flush(DEV_OUT)
        end if

        ! update min energy
        if( MinEnergy .gt. tmpgeo%total_ene ) then
            MinEnergy = tmpgeo%total_ene
            MinEnergyPtsIndex = CurrPtsIndex
        end if

        ! record new point
        Points(CurrPtsIndex)%energy = tmpgeo%total_ene
        allocate(Points(CurrPtsIndex)%crd(3,top%natoms), Points(CurrPtsIndex)%angles(NRotors), &
                 stat = alloc_status)
        if( alloc_status .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
        end if
        Points(CurrPtsIndex)%crd(:,:) = tmpgeo%crd(:,:)
        Points(CurrPtsIndex)%angles(:) = Rotors(:)%angle
        Points(CurrPtsIndex)%active_angle = rotor
    end do

    write(DEV_OUT,'(A)') '|'

    write(DEV_OUT,*)
    write(DEV_OUT,40) RejectedPoints
    write(DEV_OUT,50) CurrPtsIndex
    write(DEV_OUT,60) MinEnergy,MinEnergyPtsIndex

 10 format('First geometry energy      = ',F18.7)
 20 format('Estimated number of points = ',I10)
 30 format('Maximum number of points   = ',I10)
 35 format('Ticks per rotor bond       = ',I10)

 40 format('Rejected points            = ',I10)
 50 format('Accepted points            = ',I10)
 60 format('Minimum energy             = ',F18.7,' found for point ',I10)

end subroutine genpoints_stepbystep

!===============================================================================
! subroutine genpoints_systematic
!===============================================================================

subroutine genpoints_systematic

    use ffdev_genpoints_utils

    implicit none
    integer         :: est_npoints, i, itmp, lticks
    integer         :: nticks, alloc_status, nruns, counter
    type(GEOMETRY)  :: tmpgeo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(A)') 'Method = systematic'

    ! energy of the first structure
    call ffdev_energy_all(top,geo)
    MinEnergy = geo%total_ene
    MinEnergyPtsIndex = 1

    nticks = INT(2.0d0*DEV_PI / TickAngle)
    est_npoints = nticks**nrotors

    write(DEV_OUT,*)
    write(DEV_OUT,10) MinEnergy
    write(DEV_OUT,20) est_npoints
    write(DEV_OUT,30) MaxPoints
    write(DEV_OUT,35) nticks

    ! add the first point to stack
    CurrPtsIndex = 1
    Points(CurrPtsIndex)%energy = MinEnergy
    allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
             stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
    end if
    Points(CurrPtsIndex)%crd(:,:) = geo%crd(:,:)
    Points(CurrPtsIndex)%angles(:) = 0.0d0
    Points(CurrPtsIndex)%active_angle = 0

    ! init temporary geometry
    call ffdev_geometry_init(tmpgeo)

    RejectedPoints = 0
    nruns = 1

    write(DEV_OUT,*)
    write(DEV_OUT,'(A)',advance='NO') 'Progress: |'
    call flush(DEV_OUT)

    counter = MaxPoints / 68 + 1

    ! generate points
    do while( (CurrPtsIndex .lt. MaxPoints) .and. (nruns .lt. est_npoints) )

        ! copy initial structure
        call ffdev_geometry_copy(tmpgeo,geo)

        ! generate rotation vector
        itmp = nruns
        do i=1,NRotors
            lticks = mod(itmp,nticks)
            itmp   = itmp / nticks
!            if( lticks .gt. nticks/2 ) then
!                lticks = nticks/2 - lticks
!            end if
            Rotors(i)%angle = lticks * TickAngle
        end do

        ! generate point
        call genpoints_genrot(top,tmpgeo)

        ! get energy
        call ffdev_energy_all(top,tmpgeo)

        nruns = nruns + 1

        ! is energy acceptable
        if( tmpgeo%total_ene - MinEnergy .gt. MaxEnergy ) then
            RejectedPoints = RejectedPoints + 1
            cycle
        end if

        CurrPtsIndex = CurrPtsIndex + 1

        if( mod(CurrPtsIndex,counter) .eq. 0 ) then
            write(DEV_OUT,'(A)',advance='NO') '*'
            call flush(DEV_OUT)
        end if

        ! update min energy
        if( MinEnergy .gt. tmpgeo%total_ene ) then
            MinEnergy = tmpgeo%total_ene
            MinEnergyPtsIndex = CurrPtsIndex
        end if

        ! record new point
        Points(CurrPtsIndex)%energy = tmpgeo%total_ene
        allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
                 stat = alloc_status)
        if( alloc_status .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
        end if
        Points(CurrPtsIndex)%crd(:,:) = tmpgeo%crd(:,:)
        Points(CurrPtsIndex)%angles(:) = Rotors(:)%angle
        Points(CurrPtsIndex)%active_angle = 0
    end do

    write(DEV_OUT,'(A)') '|'

    write(DEV_OUT,*)
    write(DEV_OUT,40) RejectedPoints
    write(DEV_OUT,50) CurrPtsIndex
    write(DEV_OUT,60) MinEnergy,MinEnergyPtsIndex

 10 format('First geometry energy      = ',F18.7)
 20 format('Estimated number of points = ',I10)
 30 format('Maximum number of points   = ',I10)
 35 format('Ticks per rotor bond       = ',I10)

 40 format('Rejected points            = ',I10)
 50 format('Accepted points            = ',I10)
 60 format('Minimum energy             = ',F18.7,' found for point ',I10)

end subroutine genpoints_systematic

!===============================================================================
! subroutine genpoints_stochastic
!===============================================================================

subroutine genpoints_stochastic

    use ffdev_genpoints_utils

    implicit none
    integer         :: i
    integer         :: alloc_status, nruns, counter
    real(DEVDP)     :: angle
    type(GEOMETRY)  :: tmpgeo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(A)') 'Method = stochastic'

    ! energy of the first structure
    call ffdev_energy_all(top,geo)
    MinEnergy = geo%total_ene
    MinEnergyPtsIndex = 1

    write(DEV_OUT,*)
    write(DEV_OUT,10) MinEnergy
    write(DEV_OUT,30) MaxPoints

    ! add the first point to stack
    CurrPtsIndex = 1
    Points(CurrPtsIndex)%energy = MinEnergy
    allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
             stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
    end if
    Points(CurrPtsIndex)%crd(:,:) = geo%crd(:,:)
    Points(CurrPtsIndex)%angles(:) = 0.0d0
    Points(CurrPtsIndex)%active_angle = 0

    ! init temporary geometry
    call ffdev_geometry_init(tmpgeo)

    RejectedPoints = 0
    nruns = 1

    write(DEV_OUT,*)
    write(DEV_OUT,'(A)',advance='NO') 'Progress: |'
    call flush(DEV_OUT)

    counter = MaxPoints / 68 + 1

    ! generate points
    do while( CurrPtsIndex .lt. MaxPoints )

        ! copy initial structure
        call ffdev_geometry_copy(tmpgeo,geo)

        ! generate random rotation vector
        do i=1,NRotors
            call random_number(angle)
            Rotors(i)%angle = DEV_PI*(2.0*angle - 1.0d0)
        end do

        ! generate point
        call genpoints_genrot(top,tmpgeo)

        ! get energy
        call ffdev_energy_all(top,tmpgeo)

        nruns = nruns + 1

        ! is energy acceptable
        if( tmpgeo%total_ene - MinEnergy .gt. MaxEnergy ) then
            RejectedPoints = RejectedPoints + 1
            cycle
        end if

        CurrPtsIndex = CurrPtsIndex + 1

        if( mod(CurrPtsIndex,counter) .eq. 0 ) then
            write(DEV_OUT,'(A)',advance='NO') '*'
            call flush(DEV_OUT)
        end if

        ! update min energy
        if( MinEnergy .gt. tmpgeo%total_ene ) then
            MinEnergy = tmpgeo%total_ene
            MinEnergyPtsIndex = CurrPtsIndex
        end if

        ! record new point
        Points(CurrPtsIndex)%energy = tmpgeo%total_ene
        allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
                 stat = alloc_status)
        if( alloc_status .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
        end if
        Points(CurrPtsIndex)%crd(:,:) = tmpgeo%crd(:,:)
        Points(CurrPtsIndex)%angles(:) = Rotors(:)%angle
        Points(CurrPtsIndex)%active_angle = 0
    end do

    write(DEV_OUT,'(A)') '|'

    write(DEV_OUT,*)
    write(DEV_OUT,40) RejectedPoints
    write(DEV_OUT,50) CurrPtsIndex
    write(DEV_OUT,60) MinEnergy,MinEnergyPtsIndex

 10 format('First geometry energy      = ',F18.7)
 30 format('Maximum number of points   = ',I10)

 40 format('Rejected points            = ',I10)
 50 format('Accepted points            = ',I10)
 60 format('Minimum energy             = ',F18.7,' found for point ',I10)

end subroutine genpoints_stochastic

!===============================================================================
! subroutine genpoints_stochastic_opt
!===============================================================================

subroutine genpoints_stochastic_opt

    use ffdev_genpoints_utils
    use ffdev_geoopt
    use ffdev_gradient_utils

    implicit none
    integer         :: i
    integer         :: alloc_status, nruns, counter
    real(DEVDP)     :: angle
    type(GEOMETRY)  :: tmpgeo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(A)') 'Method = stochastic-opt'

    ! energy of the first structure
    call ffdev_energy_all(top,geo)
    MinEnergy = geo%total_ene
    MinEnergyPtsIndex = 1

    write(DEV_OUT,*)
    write(DEV_OUT,10) MinEnergy
    write(DEV_OUT,30) MaxPoints

    ! add the first point to stack
    CurrPtsIndex = 1
    Points(CurrPtsIndex)%energy = MinEnergy
    allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
             stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
    end if
    Points(CurrPtsIndex)%crd(:,:) = geo%crd(:,:)
    Points(CurrPtsIndex)%angles(:) = 0.0d0
    Points(CurrPtsIndex)%active_angle = 0

    ! init temporary geometry
    call ffdev_geometry_init(tmpgeo)

    RejectedPoints = 0
    nruns = 1

    write(DEV_OUT,*)

    counter = MaxPoints / 68 + 1

    ! generate points
    do while( CurrPtsIndex .lt. MaxPoints )

        ! copy initial structure
        call ffdev_geometry_copy(tmpgeo,geo)
        call ffdev_gradient_allocate(tmpgeo)

        ! generate random rotation vector
        do i=1,NRotors
            call random_number(angle)
            Rotors(i)%angle = DEV_PI*(2.0*angle - 1.0d0)
        end do

        ! generate point
        call genpoints_genrot(top,tmpgeo)

        ! optimize
        call ffdev_geoopt_run(DEV_OUT,top,tmpgeo)

        ! get energy
        call ffdev_energy_all(top,tmpgeo)

        nruns = nruns + 1

        ! is energy acceptable
        if( tmpgeo%total_ene - MinEnergy .gt. MaxEnergy ) then
            RejectedPoints = RejectedPoints + 1
            cycle
        end if

        CurrPtsIndex = CurrPtsIndex + 1

        ! update min energy
        if( MinEnergy .gt. tmpgeo%total_ene ) then
            MinEnergy = tmpgeo%total_ene
            MinEnergyPtsIndex = CurrPtsIndex
        end if

        ! record new point
        Points(CurrPtsIndex)%energy = tmpgeo%total_ene
        allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
                 stat = alloc_status)
        if( alloc_status .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
        end if
        Points(CurrPtsIndex)%crd(:,:) = tmpgeo%crd(:,:)
        Points(CurrPtsIndex)%angles(:) = Rotors(:)%angle
        Points(CurrPtsIndex)%active_angle = 0
    end do


    write(DEV_OUT,*)
    write(DEV_OUT,40) RejectedPoints
    write(DEV_OUT,50) CurrPtsIndex
    write(DEV_OUT,60) MinEnergy,MinEnergyPtsIndex

 10 format('First geometry energy      = ',F18.7)
 30 format('Maximum number of points   = ',I10)

 40 format('Rejected points            = ',I10)
 50 format('Accepted points            = ',I10)
 60 format('Minimum energy             = ',F18.7,' found for point ',I10)

end subroutine genpoints_stochastic_opt

!===============================================================================
! subroutine genpoints_stochastic_by_steps
!===============================================================================

subroutine genpoints_stochastic_by_steps

    use ffdev_genpoints_utils

    implicit none
    integer         :: est_npoints, i, lticks
    integer         :: nticks, alloc_status, nruns, counter
    real(DEVDP)     :: angle
    type(GEOMETRY)  :: tmpgeo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(A)') 'Method = stochastic-by-steps'

    ! energy of the first structure
    call ffdev_energy_all(top,geo)
    MinEnergy = geo%total_ene
    MinEnergyPtsIndex = 1

    nticks = INT(2.0d0*DEV_PI / TickAngle)
    est_npoints = nticks**nrotors

    write(DEV_OUT,*)
    write(DEV_OUT,10) MinEnergy
    write(DEV_OUT,30) MaxPoints
    write(DEV_OUT,35) nticks

    ! add the first point to stack
    CurrPtsIndex = 1
    Points(CurrPtsIndex)%energy = MinEnergy
    allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
             stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
    end if
    Points(CurrPtsIndex)%crd(:,:) = geo%crd(:,:)
    Points(CurrPtsIndex)%angles(:) = 0.0d0
    Points(CurrPtsIndex)%active_angle = 0

    ! init temporary geometry
    call ffdev_geometry_init(tmpgeo)

    RejectedPoints = 0
    nruns = 1

    write(DEV_OUT,*)
    write(DEV_OUT,'(A)',advance='NO') 'Progress: |'
    call flush(DEV_OUT)

    counter = MaxPoints / 68 + 1

    ! generate points
    do while( (CurrPtsIndex .lt. MaxPoints) .and. (nruns .lt. est_npoints) )

        ! copy initial structure
        call ffdev_geometry_copy(tmpgeo,geo)

        ! generate random rotation vector
        do i=1,NRotors
            call random_number(angle)
            lticks = angle*nticks
            if( lticks .gt. nticks/2 ) then
                lticks = nticks/2 - lticks
            end if
            Rotors(i)%angle = lticks * TickAngle
        end do

        ! generate point
        call genpoints_genrot(top,tmpgeo)

        ! get energy
        call ffdev_energy_all(top,tmpgeo)

        nruns = nruns + 1

        ! is energy acceptable
        if( tmpgeo%total_ene - MinEnergy .gt. MaxEnergy ) then
            RejectedPoints = RejectedPoints + 1
            cycle
        end if

        CurrPtsIndex = CurrPtsIndex + 1

        if( mod(CurrPtsIndex,counter) .eq. 0 ) then
            write(DEV_OUT,'(A)',advance='NO') '*'
            call flush(DEV_OUT)
        end if

        ! update min energy
        if( MinEnergy .gt. tmpgeo%total_ene ) then
            MinEnergy = tmpgeo%total_ene
            MinEnergyPtsIndex = CurrPtsIndex
        end if

        ! record new point
        Points(CurrPtsIndex)%energy = tmpgeo%total_ene
        allocate(Points(CurrPtsIndex)%crd(3,top%natoms),Points(CurrPtsIndex)%angles(NRotors), &
                 stat = alloc_status)
        if( alloc_status .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
        end if
        Points(CurrPtsIndex)%angles(:) = Rotors(:)%angle
        Points(CurrPtsIndex)%crd(:,:) = tmpgeo%crd(:,:)
        Points(CurrPtsIndex)%active_angle = 0
    end do

    write(DEV_OUT,'(A)') '|'

    write(DEV_OUT,*)
    write(DEV_OUT,40) RejectedPoints
    write(DEV_OUT,50) CurrPtsIndex
    write(DEV_OUT,60) MinEnergy,MinEnergyPtsIndex

 10 format('First geometry energy      = ',F18.7)
 30 format('Maximum number of points   = ',I10)
 35 format('Ticks per rotor bond       = ',I10)

 40 format('Rejected points            = ',I10)
 50 format('Accepted points            = ',I10)
 60 format('Minimum energy             = ',F18.7,' found for point ',I10)

end subroutine genpoints_stochastic_by_steps

!===============================================================================
! subroutine genpoints_nmodes
!===============================================================================

subroutine genpoints_nmodes

    use ffdev_genpoints_utils
    use ffdev_hessian
    use ffdev_hessian_utils

    implicit none
    integer         :: i,j,k,l,m,n,alloc_status
    real(DEVDP)     :: angle
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(A)') 'Method = nmodes'

    ! do we have hessian? ------------
    if( .not. geo%trg_hess_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'Initial point does not contain Hessian!')
    end if

    ! calculate frequencies ----------
    call ffdev_hessian_allocate(geo)
    call ffdev_hessian_allocate_freq(geo)
    geo%hess = geo%trg_hess
    call ffdev_hessian_calc_freqs(geo)

    ! print frequencies
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Initial Point Frequencies', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo)

    ! calculate deviations
    CurrPtsIndex = 0

    ! FIXME
    ! determine zero frequency window
    do i=7,3*top%natoms
        do j=-NModePoints,NModePoints
            if( j .eq. 0 ) cycle
            CurrPtsIndex = CurrPtsIndex + 1
            angle =j*DEV_PI/real(NModePoints)
            ! copy coordinates
            allocate(Points(CurrPtsIndex)%crd(3,top%natoms), stat = alloc_status)
            if( alloc_status .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for point coordinates!')
            end if
            Points(CurrPtsIndex)%crd = geo%crd
            ! perturbe coordinates
            do k=1,top%natoms
                do l=1,3
                    m = (i-1) / 3 + 1
                    n = mod(i-1,3) + 1
                    Points(CurrPtsIndex)%crd(l,k) = Points(CurrPtsIndex)%crd(l,k) + NModeAmplitude*sin(angle)*geo%hess(l,k,n,m)
                end do
            end do
        end do
    end do

end subroutine genpoints_nmodes

!===============================================================================
! subroutine genpoints_optimize_points
!===============================================================================

subroutine genpoints_optimize_points

    use ffdev_geometry
    use ffdev_gradient_utils
    use ffdev_geoopt

    implicit none
    type(GEOMETRY)  :: tmpgeo
    integer         :: alloc_stat, i
    ! --------------------------------------------------------------------------

    call ffdev_geometry_init(tmpgeo)
    call ffdev_geometry_allocate(tmpgeo,top%natoms)
    call ffdev_gradient_allocate(tmpgeo)

    if( HoldCV ) then
        tmpgeo%nrst = NRotors
        allocate( tmpgeo%rst(tmpgeo%nrst), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
             call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate rst array!')
        end if
        do i=1,tmpgeo%nrst
            allocate( tmpgeo%rst(i)%ai(4), stat = alloc_stat)
            if( alloc_stat .ne. 0 ) then
                 call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate rst%ai array!')
            end if
            tmpgeo%rst(i)%ai(1) = Rotors(i)%ait
            tmpgeo%rst(i)%ai(2) = Rotors(i)%ai
            tmpgeo%rst(i)%ai(3) = Rotors(i)%aj
            tmpgeo%rst(i)%ai(4) = Rotors(i)%ajt
            tmpgeo%rst(i)%cvtype = 'D'
         end do
    end if

    do i=1,CurrPtsIndex
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Optimize', ':')
        write(DEV_OUT,10) i
        write(DEV_OUT,*)

        tmpgeo%crd(:,:) = Points(i)%crd(:,:)

        if( HoldCV ) then
            do j=1,tmpgeo%nrst
                tmpgeo%rst(j)%trg_value = &
                      ffdev_geometry_get_dihedral(tmpgeo%crd,Rotors(j)%ait,Rotors(j)%ai,Rotors(j)%aj,Rotors(j)%ajt)
                tmpgeo%rst(j)%value = tmpgeo%rst(j)%trg_value
            end do
        end if

        call ffdev_geoopt_run(DEV_OUT,top,tmpgeo)

        if( i .eq. 1 ) then
            MinOptEnergy = tmpgeo%total_ene
            MinOptEnergyPtsIndex = i
        end if

        if( MinOptEnergy .gt. tmpgeo%total_ene ) then
            MinOptEnergy = tmpgeo%total_ene
            MinOptEnergyPtsIndex = i
        end if

        Points(i)%crd(:,:) = tmpgeo%crd(:,:)
        Points(i)%energy = tmpgeo%total_ene
    end do

    call ffdev_geometry_destroy(tmpgeo)

    write(DEV_OUT,*)
    write(DEV_OUT,60) MinEnergy, MinEnergyPtsIndex
    write(DEV_OUT,65) MinOptEnergy, MinOptEnergyPtsIndex

 10 format('Point #',I10)
 60 format('Minimum energy                     = ',F16.7,' found for point ',I10)
 65 format('Minimum energy after optimization  = ',F16.7,' found for point ',I10)

end subroutine genpoints_optimize_points

!===============================================================================

end program ffdev_genpoints_program

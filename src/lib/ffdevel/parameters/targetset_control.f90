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

module ffdev_targetset_control

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_targetset_ctrl
! ==============================================================================

subroutine ffdev_targetset_ctrl(fin,allow_nopoints)

    use ffdev_targetset_dat
    use prmfile
    use ffdev_utils
    use ffdev_geometry
    use ffdev_hessian_utils

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: allow_nopoints
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string,topin,key,geoname,sweight,field,wmode
    integer                     :: i,j,k,l,alloc_status,minj,probesize,npoints,ai
    logical                     :: data_avail,rst,shift2zero,stream,has_sapt
    real(DEVDP)                 :: minenergy,weight
    logical                     :: unique_probe_types
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'{TARGETS}', ':')

    ! open TARGETS group
    if( .not. prmfile_open_group(fin,'TARGETS') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open {TARGETS} group!')
    end if

    ! load [setup] configuration
    call ffdev_targetset_ctrl_setup(fin)

    ! get number of sets
    nsets = 0
    rst = prmfile_first_section(fin)
    do while( rst )
        ! open section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to get section name!')
        end if

        if( string .eq. 'SET' ) then
            nsets = nsets + 1
        end if
        rst = prmfile_next_section(fin)
    end do

    if( nsets .lt. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'At least one target set must be defined!')
    end if

    ! allocate sets
    allocate(sets(nsets), stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for target sets!')
    end if

    ! pre-init
    sets(:)%isref       = .false.
    sets(:)%name        = ''

    rst = prmfile_first_section(fin)

    ! load sets ----------------------------------
    i = 1
    do while( rst )

        ! open set section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to get section name!')
        end if

        if( string .ne. 'SET' ) then
            rst = prmfile_next_section(fin)
            cycle
        end if

        write(DEV_OUT,*)
        write(DEV_OUT,1) i

!----------------------
! topology
!----------------------
        ! get name - optional
        if( .not. prmfile_get_string_by_key(fin,'name',sets(i)%name) ) then
            sets(i)%name = ''
        else
            write(DEV_OUT,5) trim(sets(i)%name)
        end if

        ! get topology name
        if( .not. prmfile_get_string_by_key(fin,'topology',topin) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to get topology name (topology) for [SET]!')
        end if
        write(DEV_OUT,12) trim(topin)

        probesize = 0
        if( prmfile_get_integer_by_key(fin,'probesize',probesize) ) then
            write(DEV_OUT,20) probesize
        else
            write(DEV_OUT,25) probesize
        end if

        unique_probe_types = .true.
        if( prmfile_get_logical_by_key(fin,'unique_probe_types', unique_probe_types)) then
            write(DEV_OUT,30) prmfile_onoff(unique_probe_types)
        else
            write(DEV_OUT,35) prmfile_onoff(unique_probe_types)
        end if

        ! load topology and print info
        write(DEV_OUT,*)
        call ffdev_topology_init(sets(i)%top)
        call ffdev_topology_load(sets(i)%top,topin)
        call ffdev_topology_switch_to_probe_mode(sets(i)%top,probesize,unique_probe_types)
        call ffdev_topology_info(sets(i)%top)

        write(DEV_OUT,*)
        sets(i)%isref = .false.
        if( prmfile_get_logical_by_key(fin,'isreference', sets(i)%isref)) then
            write(DEV_OUT,71) prmfile_onoff(sets(i)%isref)
        else
            write(DEV_OUT,75) prmfile_onoff(sets(i)%isref)
        end if

        if( sets(i)%isref ) then
            if( len(trim(sets(i)%name)) .eq. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Name is not specified for a reference set!')
            end if
        end if

        sets(i)%nrefs = 0
        nullify(sets(i)%refs)
        nullify(sets(i)%natomsrefs)

        if( prmfile_init_field_by_key(fin,'references') ) then
            if( sets(i)%isref ) then
                call ffdev_utils_exit(DEV_ERR,1,'References cannot be defined for a reference set!')
            end if

            do while ( prmfile_get_field_by_key(fin,field) )
                sets(i)%nrefs = sets(i)%nrefs + 1
            end do
            allocate(sets(i)%refs(sets(i)%nrefs),sets(i)%natomsrefs(sets(i)%nrefs), stat = alloc_status)
            if( alloc_status .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for references!')
            end if
            rst = prmfile_init_field_by_key(fin,'references')
            k = 1
            do while ( prmfile_get_field_by_key(fin,field) )
                ! find set
                sets(i)%refs(k) = 0
                do l=1,nsets
                    if( sets(l)%isref .and. (trim(sets(l)%name) .eq. trim(field)) ) then
                        sets(i)%refs(k) = l
                        exit
                    end if
                end do
                if( sets(i)%refs(k) .eq. 0 ) then
                    call ffdev_utils_exit(DEV_ERR,1,'Unable to find reference (' // trim(field) // ')!')
                end if
                k = k + 1
            end do
            rst = prmfile_get_string_by_key(fin,'references', string)
            write(DEV_OUT,70) trim(string)
        ! check size consistency
            k = 0
            do l=1,sets(i)%nrefs
                sets(i)%natomsrefs(l) = sets(sets(i)%refs(l))%top%natoms
                k = k + sets(i)%natomsrefs(l)
            end do
            if( k .ne. sets(i)%top%natoms ) then
                call ffdev_utils_exit(DEV_ERR,1,'Size (natoms) inconsistency detected between the set and reference systems!')
            end if
        ! check order consistency - needed for charges
            ai = 1
            do l=1,sets(i)%nrefs
                do k=1,sets(sets(i)%refs(l))%top%natoms
                    if( trim(sets(sets(i)%refs(l))%top%atoms(k)%name) .ne. trim(sets(i)%top%atoms(ai)%name) ) then
                        write(string,'(A,I3,A,I3)') 'Inconsistency in order of references (name)! Ref atom: ', &
                                                     k, ' Set atom: ', ai
                        call ffdev_utils_exit(DEV_ERR,1,string)
                    end if
                    if( trim(sets(sets(i)%refs(l))%top%atoms(k)%resname) .ne. trim(sets(i)%top%atoms(ai)%resname) ) then
                        write(string,'(A,I3,A,I3)') 'Inconsistency in order of references (resname)! Ref atom: ', &
                                                    k, ' Set atom: ', ai
                        call ffdev_utils_exit(DEV_ERR,1,string)
                    end if
                    ai = ai + 1
                end do
            end do

        end if

        shift2zero = .false.
        if( sets(i)%nrefs .eq. 0 ) then
            if( .not. prmfile_get_logical_by_key(fin,'shift2zero',shift2zero) ) then
                write(DEV_OUT,65) prmfile_onoff(shift2zero)
            else
                write(DEV_OUT,60) prmfile_onoff(shift2zero)
            end if
        else
            ! we cannot manipulate with energy in reference mode
            write(DEV_OUT,65) prmfile_onoff(shift2zero)
        end if

        sets(i)%optgeo = .false.
        if( prmfile_get_logical_by_key(fin,'optgeo', sets(i)%optgeo)) then
            write(DEV_OUT,40) prmfile_onoff(sets(i)%optgeo)
        else
            write(DEV_OUT,45) prmfile_onoff(sets(i)%optgeo)
        end if

        sets(i)%keepoptgeo = .false.
        if( prmfile_get_logical_by_key(fin,'keepoptgeo', sets(i)%keepoptgeo)) then
            write(DEV_OUT,50) prmfile_onoff(sets(i)%keepoptgeo)
        else
            write(DEV_OUT,55) prmfile_onoff(sets(i)%keepoptgeo)
        end if

        sets(i)%nofreq = DoNotCalcFreqs
        if( prmfile_get_logical_by_key(fin,'nofreq', sets(i)%nofreq)) then
            write(DEV_OUT,110) prmfile_onoff(sets(i)%nofreq)
        else
            write(DEV_OUT,115) prmfile_onoff(sets(i)%nofreq)
        end if

        sets(i)%savepts  = SavePts
        sets(i)%savexyzr = SaveXYZR
        if( prmfile_get_logical_by_key(fin,'savepts', sets(i)%savepts)) then
            write(DEV_OUT,80) prmfile_onoff(sets(i)%savepts)
        else
            write(DEV_OUT,85) prmfile_onoff(sets(i)%savepts)
        end if
        if( prmfile_get_logical_by_key(fin,'savexyzr', sets(i)%savexyzr)) then
            write(DEV_OUT,120) prmfile_onoff(sets(i)%savexyzr)
        else
            write(DEV_OUT,125) prmfile_onoff(sets(i)%savexyzr)
        end if

        wmode = 'auto'
        if( prmfile_get_string_by_key(fin,'wmode', wmode)) then
            write(DEV_OUT,90) trim(wmode)
        else
            write(DEV_OUT,95) trim(wmode)
        end if

        ! get name of initial driving profile
        if( .not. prmfile_get_string_by_key(fin,'initial_drvene',sets(i)%initial_drvene) ) then
            sets(i)%initial_drvene = ''
            write(DEV_OUT,18) trim('-none-')
        else
            write(DEV_OUT,18) trim(sets(i)%initial_drvene)
        end if

        ! get name of initial driving structures
        if( .not. prmfile_get_string_by_key(fin,'initial_drvxyz',sets(i)%initial_drvxyz) ) then
            sets(i)%initial_drvxyz = ''
            write(DEV_OUT,19) trim('-none-')
        else
            write(DEV_OUT,19) trim(sets(i)%initial_drvxyz)
        end if

        ! get name of resulting topology
        if( .not. prmfile_get_string_by_key(fin,'final_topology',sets(i)%final_stop) ) then
            sets(i)%final_stop = ''
            write(DEV_OUT,15) trim('-none-')
        else
            write(DEV_OUT,15) trim(sets(i)%final_stop)
        end if

        ! get name of resulting driving profile
        if( .not. prmfile_get_string_by_key(fin,'final_drvene',sets(i)%final_drvene) ) then
            sets(i)%final_drvene = ''
            write(DEV_OUT,16) trim('-none-')
        else
            write(DEV_OUT,16) trim(sets(i)%final_drvene)
        end if

        ! get name of resulting driving structures
        if( .not. prmfile_get_string_by_key(fin,'final_drvxyz',sets(i)%final_drvxyz) ) then
            sets(i)%final_drvxyz = ''
            write(DEV_OUT,17) trim('-none-')
        else
            write(DEV_OUT,17) trim(sets(i)%final_drvxyz)
        end if

!----------------------
! points
!----------------------

        ! count number of points
        sets(i)%ngeos = 0
        rst = prmfile_first_line(fin)
        do while( rst )
            rst = prmfile_get_field_on_line(fin,key)
            select case(trim(key))
                case('point')
                    sets(i)%ngeos = sets(i)%ngeos + 1
                case('points')
                    rst = prmfile_get_field_on_line(fin,geoname)
                    if( .not. rst ) then
                        call ffdev_utils_exit(DEV_ERR,1,'Name of file with stream of points is not specified!')
                    end if
                    sets(i)%ngeos = sets(i)%ngeos +  ffdev_geometry_num_of_points(geoname)
                case default
                    ! do nothing
            end select
            rst = prmfile_next_line(fin)
        end do

        if( sets(i)%ngeos .le. 0 ) then
            if( allow_nopoints .eqv. .false. ) then
                call ffdev_utils_exit(DEV_ERR,1,'No points for current set!')
            end if
        end if

        if( sets(i)%isref .and. (sets(i)%ngeos .ne. 1) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Reference set can contain only one point!')
        end if

        if( sets(i)%ngeos .gt. 0 ) then
            ! allocate points
            allocate(sets(i)%geo(sets(i)%ngeos), stat = alloc_status)
            if( alloc_status .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for set points!')
            end if
        end if

        write(DEV_OUT,*)
        write(DEV_OUT,200) sets(i)%ngeos

        if( sets(i)%ngeos .gt. 0 ) then
            write(DEV_OUT,*)
            if( shift2zero .or. sets(i)%nrefs .gt. 0 .or. sets(i)%isref ) then
                call ffdev_geometry_info_point_header_ext(.false.)
            else
                call ffdev_geometry_info_point_header()
            end if
        end if

        ! load points
        rst = prmfile_first_line(fin)
        j = 1
        do while( rst )

            rst = prmfile_get_field_on_line(fin,key)

            ! process only points
            if( (trim(key) .ne. 'point') .and. (trim(key) .ne. 'points') ) then
                rst = prmfile_next_line(fin)
                cycle
            end if

            rst = prmfile_get_field_on_line(fin,geoname)
            if( .not. rst ) then
                call ffdev_utils_exit(DEV_ERR,1,'Point name not specified!')
            end if

            weight = 1.0d0
            rst = prmfile_get_field_on_line(fin,sweight)
            if( rst ) then
                read(sweight,*,end=100,err=100) weight
100             continue
            end if

            select case(trim(key))
                case('point')
                    npoints = 1
                    call ffdev_utils_open(DEV_GEO,geoname,'O')
                    stream = .false.
                case('points')
                    npoints = ffdev_geometry_num_of_points(geoname)
                    call ffdev_utils_open(DEV_GEO,geoname,'O')
                    stream = .true.
            end select

            do k=1,npoints

                if( j .gt. sets(i)%ngeos ) then
                    call ffdev_utils_exit(DEV_ERR,1,'Mismatch in the number of data points!')
                end if

                ! load point and print its info
                call ffdev_geometry_init(sets(i)%geo(j))
                sets(i)%geo(j)%id = j
                sets(i)%geo(j)%weight = weight
                ! always load trg_crd as the geo optimization can be switched on any later
                sets(i)%geo(j)%trg_crd_loaded  = .true.
                sets(i)%geo(j)%name = geoname
                call ffdev_geometry_load_1point(sets(i)%geo(j),stream)

                ! calculate traget freqs and normal mode if needed
                if(  (.not. sets(i)%nofreq) .and. sets(i)%geo(j)%trg_hess_loaded ) then
                    call ffdev_hessian_calc_trg_freqs(sets(i)%geo(j))
                    sets(i)%geo(j)%trg_freq_loaded =  .true.
                end if

                if( shift2zero .or. sets(i)%nrefs .gt. 0 .or. sets(i)%isref ) then
                    call ffdev_geometry_info_point_ext(sets(i)%geo(j),.false.)
                else
                    call ffdev_geometry_info_point(sets(i)%geo(j))
                end if

                ! overwrite weights
                select case(trim(wmode))
                    case('auto')
                        ! nothing to do
                    case('ire2')
                        if( sets(i)%geo(j)%trg_energy .ne. 0d0 ) then
                            sets(i)%geo(j)%weight = 1.0d0 / sets(i)%geo(j)%trg_energy**2
                        end if
                    case default
                        call ffdev_utils_exit(DEV_ERR,1,'Unsupported wmode ''' // trim(wmode) // '''!')
                end select

                call ffdev_geometry_check_z(sets(i)%top,sets(i)%geo(j))

                ! reference point needs energy
                if( sets(i)%isref ) then
                    if( sets(i)%geo(j)%trg_ene_loaded .eqv. .false. ) then
                        call ffdev_utils_exit(DEV_ERR,1,'Reference point must contain ENERGY definition!')
                    end if
                end if

                ! do we have data for point
                data_avail = sets(i)%geo(j)%trg_ene_loaded .or. &
                             sets(i)%geo(j)%trg_grd_loaded .or. &
                             sets(i)%geo(j)%trg_hess_loaded .or. &
                             sets(i)%optgeo .or. sets(i)%isref

                if( .not. data_avail ) then
                    call ffdev_utils_exit(DEV_ERR,1,'No training data are available in the point!')
                end if
                j = j + 1
            end do

            ! close the file
            close(DEV_GEO)

            rst = prmfile_next_line(fin)
        end do

        ! update points with references
        if( sets(i)%nrefs .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,305)
            write(DEV_OUT,*)
            call ffdev_geometry_info_point_header_ext(.true.)
            do j=1,sets(i)%ngeos
                if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle
                do k=1,sets(i)%nrefs
                    sets(i)%geo(j)%trg_energy = sets(i)%geo(j)%trg_energy - sets( sets(i)%refs(k) )%geo(1)%trg_energy
                end do
                call ffdev_geometry_info_point_ext(sets(i)%geo(j),.true.)
            end do
        end if

        ! do energy statistics
        minj = 0
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_ene_loaded ) then
                minenergy = sets(i)%geo(j)%trg_energy
                minj = j
                exit
            end if
        end do
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle
            if( minenergy .gt. sets(i)%geo(j)%trg_energy ) then
                minenergy = sets(i)%geo(j)%trg_energy
                minj = j
            end if
        end do
        sets(i)%mineneid = minj
        write(DEV_OUT,*)
        write(DEV_OUT,300) sets(i)%mineneid,minenergy

        if( (sets(i)%mineneid .gt. 0) .and. (sets(i)%nrefs .eq. 0) ) then
            if( shift2zero ) then
                write(DEV_OUT,*)
                write(DEV_OUT,308)
                write(DEV_OUT,*)
                call ffdev_geometry_info_point_header_ext(.true.)
                do j=1,sets(i)%ngeos
                    if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle
                    sets(i)%geo(j)%trg_energy = sets(i)%geo(j)%trg_energy - minenergy
                    call ffdev_geometry_info_point_ext(sets(i)%geo(j),.true.)
                end do
            end if
        end if

        i = i + 1

        rst = prmfile_next_section(fin)

    end do

    ! build SAPT list for sets with SAPT data
    do i=1,nsets
        has_sapt = .false.
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_sapt_loaded ) then
                has_sapt = .true.
                exit
            end if
        end do
        if( has_sapt ) then
            ! gen SAPT list
            if( sets(i)%top%probe_size .gt. 0 ) then
                call ffdev_topology_gen_sapt_list_for_probes(sets(i)%top)
            else if( sets(i)%nrefs .ge. 2 ) then
                call ffdev_topology_gen_sapt_list_for_refs(sets(i)%top,sets(i)%nrefs,sets(i)%natomsrefs)
            end if
        end if
    end do

  1 format('=== [SET] #',I2.2,' ==================================================================')
  5 format('Set name (name)                         = ',A)
 12 format('Input topology name (topology)          = ',A)
 15 format('Final topology name (final_topology)    = ',A)
 16 format('Final driving profile (final_drvene)    = ',A)
 17 format('Final driving structures (final_drvxyz) = ',A)
 18 format('Initial drv profile (initial_drvene)    = ',A)
 19 format('Initial drv structures (initial_drvxyz) = ',A)

 20 format('Probe size (probesize)                  = ',I12)
 25 format('Probe size (probesize)                  = ',I12,'                  (default)')

 30 format('Unique probe (unique_probe_types)       = ',A12)
 35 format('Unique probe (unique_probe_types)       = ',A12,'                  (default)')

 40 format('Optimize geometry (optgeo)              = ',A12)
 45 format('Optimize geometry (optgeo)              = ',A12,'                  (default)')

 50 format('Keep optimized geometry (keepoptgeo)    = ',A12)
 55 format('Keep optimized geometry (keepoptgeo)    = ',A12,'                  (default)')

 60 format('Shift minimum to zero (shift2zero)      = ',A12)
 65 format('Shift minimum to zero (shift2zero)      = ',A12,'                  (default)')

 71 format('Set as reference set (isreference)      = ',A12)
 75 format('Set as reference set (isreference)      = ',A12,'                  (default)')

 70 format('Reference sets (references)             = ',A)

 80 format('Save final point geometry (savepts)     = ',a12)
 85 format('Save final point geometry (savepts)     = ',a12,'                  (default)')

120 format('Save final xyzr geometry (savexyzr)     = ',a12)
125 format('Save final xyzr geometry (savexyzr)     = ',a12,'                  (default)')

 90 format('Point weights mode (wmode)              = ',a12)
 95 format('Point weights mode (wmode)              = ',a12,'                  (default)')

110 format('Do not calculate frequencies (nofreq)   = ',a12)
115 format('Do not calculate frequencies (nofreq)   = ',a12,'                  (default)')

200 format('Number of target points                 = ',I6)
300 format('Minimum energy point #',I5.5,' has energy ',F20.4)
305 format('Substracting energy of reference states ...')
308 format('Shifting minimum to zero ...')

end subroutine ffdev_targetset_ctrl

! ==============================================================================
! subroutine ffdev_targetset_ctrl_setup
! ==============================================================================

subroutine ffdev_targetset_ctrl_setup(fin)

    use ffdev_targetset_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    integer             :: estat
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'setup') ) then
        write(DEV_OUT,105) prmfile_onoff(DoNotCalcFreqs)
        write(DEV_OUT,135) prmfile_onoff(UseOptGeometry)
        write(DEV_OUT,145) prmfile_onoff(KeepPtsNames)
        write(DEV_OUT,85)  prmfile_onoff(SavePts)
        write(DEV_OUT,95)  SavePtsPath
        write(DEV_OUT,115) prmfile_onoff(SaveXYZR)
        write(DEV_OUT,125) SaveXYZRPath
        write(DEV_OUT,155) prmfile_onoff(SaveSumLogs)
        write(DEV_OUT,165) SaveSumLogsPath

        ! create directories
        if( SavePts ) then
            call execute_command_line('mkdir -p ' // trim(SavePtsPath), exitstat = estat )
            if( estat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to create SavePtsPath!')
            end if
        end if
        if( SaveXYZR ) then
            call execute_command_line('mkdir -p ' // trim(SaveXYZRPath), exitstat = estat  )
            if( estat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to create SaveXYZRPath!')
            end if
        end if
        if( SaveSumLogs ) then
            call execute_command_line('mkdir -p ' // trim(SaveSumLogsPath), exitstat = estat  )
            if( estat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to create SaveSumLogsPath!')
            end if
        end if
        return
    end if

    if( prmfile_get_logical_by_key(fin,'nofreqs', DoNotCalcFreqs)) then
        write(DEV_OUT,100) prmfile_onoff(DoNotCalcFreqs)
    else
        write(DEV_OUT,105) prmfile_onoff(DoNotCalcFreqs)
    end if
    if( prmfile_get_logical_by_key(fin,'useoptgeo', UseOptGeometry)) then
        write(DEV_OUT,130) prmfile_onoff(UseOptGeometry)
    else
        write(DEV_OUT,135) prmfile_onoff(UseOptGeometry)
    end if
    if( prmfile_get_logical_by_key(fin,'keepptsnames', KeepPtsNames)) then
        write(DEV_OUT,140) prmfile_onoff(KeepPtsNames)
    else
        write(DEV_OUT,145) prmfile_onoff(KeepPtsNames)
    end if

    if( prmfile_get_logical_by_key(fin,'savepts', SavePts)) then
        write(DEV_OUT,80) prmfile_onoff(SavePts)
    else
        write(DEV_OUT,85) prmfile_onoff(SavePts)
    end if
    if( prmfile_get_string_by_key(fin,'saveptspath', SavePtsPath)) then
        write(DEV_OUT,90) SavePtsPath
    else
        write(DEV_OUT,95) SavePtsPath
    end if

    if( prmfile_get_logical_by_key(fin,'savexyzr', SaveXYZR)) then
        write(DEV_OUT,110) prmfile_onoff(SaveXYZR)
    else
        write(DEV_OUT,115) prmfile_onoff(SaveXYZR)
    end if
    if( prmfile_get_string_by_key(fin,'savexyzpath', SaveXYZRPath)) then
        write(DEV_OUT,120) SaveXYZRPath
    else
        write(DEV_OUT,120) SaveXYZRPath
    end if

    if( prmfile_get_logical_by_key(fin,'savesumlogs', SaveSumLogs)) then
        write(DEV_OUT,150) prmfile_onoff(SaveSumLogs)
    else
        write(DEV_OUT,155) prmfile_onoff(SaveSumLogs)
    end if
    if( prmfile_get_string_by_key(fin,'savesumlogspath', SaveSumLogsPath)) then
        write(DEV_OUT,160) SaveSumLogsPath
    else
        write(DEV_OUT,165) SaveSumLogsPath
    end if

    ! create directories
    if( SavePts ) then
        call execute_command_line('mkdir -p ' // trim(SavePtsPath), exitstat = estat )
        if( estat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to create SavePtsPath!')
        end if
    end if
    if( SaveXYZR ) then
        call execute_command_line('mkdir -p ' // trim(SaveXYZRPath), exitstat = estat )
        if( estat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to create SaveXYZRPath!')
        end if
    end if
    if( SaveSumLogs ) then
        call execute_command_line('mkdir -p ' // trim(SaveSumLogsPath), exitstat = estat  )
        if( estat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to create SaveSumLogsPath!')
        end if
    end if

 10 format('=== [setup] ====================================================================')

100  format ('Do not calculate frequencies (nofreqs)    = ',a12)
105  format ('Do not calculate frequencies (nofreqs)    = ',a12,'               (default)')

130  format ('Use optimized geometry (useoptgeo)        = ',a12)
135  format ('Use optimized geometry (useoptgeo)        = ',a12,'               (default)')

140  format ('Keep point names (keepptsnames)           = ',a12)
145  format ('Keep point names (keepptsnames)           = ',a12,'               (default)')

 80  format ('Save geometry as points (savepts)         = ',a12)
 85  format ('Save geometry as points (savepts)         = ',a12,'               (default)')

 90  format ('Save points path (saveptspath)            = ',a25)
 95  format ('Save points path (saveptspath)            = ',a25,'  (default)')

110  format ('Save geometry as xyzr (savexyzr)          = ',a12)
115  format ('Save geometry as xyzr (savexyzr)          = ',a12,'               (default)')

120  format ('Save xyzr path (savexyzrpath)             = ',a25)
125  format ('Save xyzr path (savexyzrpath)             = ',a25,'  (default)')

150  format ('Save error summary (savesumlogs)          = ',a12)
155  format ('Save error summary (savesumlogs)          = ',a12,'               (default)')

160  format ('Save error summary path (savesumlogspath) = ',a25)
165  format ('Save error summary path (savesumlogspath) = ',a25,'  (default)')

end subroutine ffdev_targetset_ctrl_setup

! ==============================================================================
! subroutine ffdev_targetset_ctrl_optgeo
! ==============================================================================

subroutine ffdev_targetset_ctrl_optgeo_set_default

    use ffdev_targetset_dat

    implicit none
    ! --------------------------------------------------------------------------

    ! setup default parameters
    GlbOptGeometryEnabled   = .false.   ! force all geometry optimization in each error evaluation
    GlbOptGeometryDisabled  = .false.   ! disable all geometry optimization in each error evaluation
    GlbShowOptProgress      = .false.   ! print geometry optimization progress
    GlbKeepOptGeometry      = .false.   ! keep geometry from previous geometry optimization
    ResetKeptOptGeoAt       = 0         ! reset initial geometry to trg for keepoptgeo = on

end subroutine ffdev_targetset_ctrl_optgeo_set_default

! ==============================================================================
! subroutine ffdev_targetset_ctrl_optgeo
! ==============================================================================

subroutine ffdev_targetset_ctrl_optgeo(fin)

    use ffdev_targetset_dat
    use prmfile
    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( ResetAllSetup ) then
        call ffdev_targetset_ctrl_optgeo_set_default()
    end if

    if( .not. prmfile_open_section(fin,'optgeo') ) then
        write(DEV_OUT,35) prmfile_onoff(GlbOptGeometryEnabled)
        write(DEV_OUT,37) prmfile_onoff(GlbOptGeometryDisabled)
        write(DEV_OUT,45) prmfile_onoff(GlbShowOptProgress)
        write(DEV_OUT,55) prmfile_onoff(GlbKeepOptGeometry)
        write(DEV_OUT,65) ResetKeptOptGeoAt
        return
    end if

    if( prmfile_get_logical_by_key(fin,'enabled', GlbOptGeometryEnabled)) then
        write(DEV_OUT,30) prmfile_onoff(GlbOptGeometryEnabled)
    else
        write(DEV_OUT,35) prmfile_onoff(GlbOptGeometryEnabled)
    end if
    if( prmfile_get_logical_by_key(fin,'disabled', GlbOptGeometryDisabled)) then
        write(DEV_OUT,36) prmfile_onoff(GlbOptGeometryDisabled)
    else
        write(DEV_OUT,37) prmfile_onoff(GlbOptGeometryDisabled)
    end if
    if( prmfile_get_logical_by_key(fin,'optprogress', GlbShowOptProgress)) then
        write(DEV_OUT,40) prmfile_onoff(GlbShowOptProgress)
    else
        write(DEV_OUT,45) prmfile_onoff(GlbShowOptProgress)
    end if
    if( prmfile_get_logical_by_key(fin,'keepoptgeo', GlbKeepOptGeometry)) then
        write(DEV_OUT,50) prmfile_onoff(GlbKeepOptGeometry)
    else
        write(DEV_OUT,55) prmfile_onoff(GlbKeepOptGeometry)
    end if
    if( prmfile_get_integer_by_key(fin,'resetat', ResetKeptOptGeoAt)) then
        write(DEV_OUT,60) ResetKeptOptGeoAt
    else
        write(DEV_OUT,65) ResetKeptOptGeoAt
    end if

 10 format('=== [optgeo] ===================================================================')

 30  format ('Force all geo optimization (enabled)   = ',a12)
 35  format ('Force all geo optimization (enabled)   = ',a12,'                  (default)')

 36  format ('Disable all geo optimization (disabled)= ',a12)
 37  format ('Disable all geo optimization (disabled)= ',a12,'                  (default)')

 40  format ('Show geo opt progress (optprogress)    = ',a12)
 45  format ('Show geo opt progress (optprogress)    = ',a12,'                  (default)')

 50  format ('Keep optimized geometry (keepoptgeo)   = ',a12)
 55  format ('Keep optimized geometry (keepoptgeo)   = ',a12,'                  (default)')

 60  format ('Reset kept geometry after (resetat)    = ',i12)
 65  format ('Reset kept geometry after (resetat)    = ',i12,'                  (default)')

end subroutine ffdev_targetset_ctrl_optgeo

! ------------------------------------------------------------------------------

end module ffdev_targetset_control

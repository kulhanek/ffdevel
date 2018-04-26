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

contains

! ==============================================================================
! subroutine ffdev_targetset_ctrl
! ==============================================================================

subroutine ffdev_targetset_ctrl(fin,allow_nopoints)

    use ffdev_targetset_dat
    use prmfile
    use ffdev_utils
    use ffdev_geometry

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: allow_nopoints
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string,topin,key,geoname,sweight,field
    integer                     :: i,j,k,l,alloc_status,minj,probesize,nb_mode,comb_rules
    logical                     :: data_avail,rst,shift2zero
    real(DEVDP)                 :: minenergy,weight
    logical                     :: unique_probe_types
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'{TARGETS}', ':')

    ! open TARGETS group
    if( .not. prmfile_open_group(fin,'TARGETS') ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to open {TARGETS} group!')
    end if
    
    ! load [setup] configuration
    call ffdev_targetset_ctrl_setup(fin)

    ! get number of sets
    nsets = 0
    rst = prmfile_first_section(fin)
    do while( rst )
        ! open section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get section name!')
        end if

        if( string .eq. 'SET' ) then
            nsets = nsets + 1
        end if
        rst = prmfile_next_section(fin)
    end do
        
    if( nsets .lt. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'At least one target set must be defined!')
    end if

    ! allocate sets
    allocate(sets(nsets), stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for target sets!')
    end if

    rst = prmfile_first_section(fin)

    ! load sets ----------------------------------
    i = 1
    do while( rst )

        ! open set section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get section name!')
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
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get topology name (topology) for [SET]!')
        end if
        write(DEV_OUT,12) trim(topin)

        ! get name of resulting topology
        if( .not. prmfile_get_string_by_key(fin,'final',sets(i)%final_stop) ) then
            sets(i)%final_stop = ''
            write(DEV_OUT,15) trim('-none-')
        else
            write(DEV_OUT,15) trim(sets(i)%final_stop)
        end if
        
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
        sets(i)%optgeo = OptimizeGeometry
        if( prmfile_get_logical_by_key(fin,'optgeo', sets(i)%optgeo)) then
            write(DEV_OUT,40) prmfile_onoff(sets(i)%optgeo)
        else
            write(DEV_OUT,45) prmfile_onoff(sets(i)%optgeo)
        end if 
        
        sets(i)%keepoptgeo = KeepOptimizedGeometry
        if( prmfile_get_logical_by_key(fin,'keepoptgeo', sets(i)%keepoptgeo)) then
            write(DEV_OUT,50) prmfile_onoff(sets(i)%keepoptgeo)
        else
            write(DEV_OUT,55) prmfile_onoff(sets(i)%keepoptgeo)
        end if  
        
        sets(i)%nrefs = 0 
        nullify(sets(i)%refs)
        
        if( prmfile_init_field_by_key(fin,'references') ) then
            do while ( prmfile_get_field_by_key(fin,field) )
                sets(i)%nrefs = sets(i)%nrefs + 1
            end do
            allocate(sets(i)%refs(sets(i)%nrefs), stat = alloc_status)
            if( alloc_status .ne. 0 ) then
                call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for references!')
            end if
            rst = prmfile_init_field_by_key(fin,'references')
            k = 1
            do while ( prmfile_get_field_by_key(fin,field) )
                ! find set
                sets(i)%refs(k) = 0
                do l=1,nsets
                    if( trim(sets(l)%name) .eq. trim(field) ) then
                        sets(i)%refs(k) = l
                        exit
                    end if
                end do
                if( sets(i)%refs(k) .eq. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to find reference (' // trim(field) // ')!')
                end if
                k = k + 1
            end do
            rst = prmfile_get_string_by_key(fin,'references', string)
            write(DEV_OUT,70) trim(string)
            ! check size consitency
            k = 0
            do l=1,sets(i)%nrefs
                k = k + sets(sets(i)%refs(l))%top%natoms
            end do
            if( k .ne. sets(i)%top%natoms ) then
                call ffdev_utils_exit(DEV_OUT,1,'Size (natoms) inconsistency between reference systems and the set detected!')
            end if
        end if
        
        sets(i)%offset = 0.0d0
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
                case default
                    ! do nothing
            end select
            rst = prmfile_next_line(fin)
        end do

        if( sets(i)%ngeos .le. 0 ) then
            if( allow_nopoints .eqv. .false. ) then
                call ffdev_utils_exit(DEV_OUT,1,'No points for current set!')
            end if
        end if
        
        if( (len(trim(sets(i)%name)) .gt. 0) .and. (sets(i)%ngeos .ne. 1) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Named set can contain only one point!')
        end if

        if( sets(i)%ngeos .gt. 0 ) then
            ! allocate points
            allocate(sets(i)%geo(sets(i)%ngeos), stat = alloc_status)
            if( alloc_status .ne. 0 ) then
                call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for set points!')
            end if
        end if

        write(DEV_OUT,*)
        write(DEV_OUT,200) sets(i)%ngeos

        if( sets(i)%ngeos .gt. 0 ) then
            write(DEV_OUT,*)
            if( shift2zero ) then
                call ffdev_geometry_info_point_header()
            else
                call ffdev_geometry_info_point_header_ext()
            end if
        end if

        ! load points
        rst = prmfile_first_line(fin)
        j = 1
        do while( rst )

            rst = prmfile_get_field_on_line(fin,key)

            ! process only points
            if( trim(key) .ne. 'point' ) then
                rst = prmfile_next_line(fin)
                cycle
            end if

            rst = prmfile_get_field_on_line(fin,geoname)
            if( .not. rst ) then
                call ffdev_utils_exit(DEV_OUT,1,'Point name not specified!')
            end if

            weight = 1.0d0
            rst = prmfile_get_field_on_line(fin,sweight)
            if( rst ) then
                read(sweight,*,end=100,err=100) weight
100             continue
            end if

            if( j .gt. sets(i)%ngeos ) then
                call ffdev_utils_exit(DEV_OUT,1,'Mismatch in the number of data points!')
            end if

            ! load point and print its info
            call ffdev_geometry_init(sets(i)%geo(j))
            sets(i)%geo(j)%id = j
            sets(i)%geo(j)%weight = weight
            sets(i)%geo(j)%trg_crd_loaded  = sets(i)%optgeo
            call ffdev_geometry_load_point(sets(i)%geo(j),geoname)
            if( shift2zero ) then
                call ffdev_geometry_info_point(sets(i)%geo(j))
            else
                call ffdev_geometry_info_point_ext(sets(i)%geo(j))
            end if

            call ffdev_geometry_check_z(sets(i)%top,sets(i)%geo(j))

            ! do we have data for point
            data_avail = sets(i)%geo(j)%trg_ene_loaded .or. &
                         sets(i)%geo(j)%trg_grd_loaded .or. &
                         sets(i)%geo(j)%trg_hess_loaded .or. &
                         sets(i)%optgeo .or. ( len(sets(i)%name) .gt. 0)

            if( .not. data_avail ) then
                call ffdev_utils_exit(DEV_OUT,1,'No training data are available in the point!')
            end if
            j = j + 1

            rst = prmfile_next_line(fin)
        end do

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
        if( sets(i)%mineneid .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,300) sets(i)%mineneid,minenergy
            if( shift2zero ) then
                call ffdev_geometry_info_point_header_ext()
                do j=1,sets(i)%ngeos
                    if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle
                    sets(i)%geo(j)%trg_energy = sets(i)%geo(j)%trg_energy - minenergy
                    write(DEV_OUT,*)
                    call ffdev_geometry_info_point_ext(sets(i)%geo(j))
                end do
            end if
        end if

        i = i + 1
        rst = prmfile_next_section(fin)

    end do

  1 format('=== [SET] #',I2.2,' ==================================================================')
  5 format('Set name (name)                         = ',A) 
 12 format('Input topology name (topology)          = ',A)
 15 format('Final topology name (final)             = ',A)
 
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
  
 70 format('Reference sets (references)             = ',A)

200 format('Number of target points                 = ',I6)
300 format('Minimum energy point #',I5.5,' has energy ',F20.4)

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
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'setup') ) then
        write(DEV_OUT,25) prmfile_onoff(ApplyCombinationRules)
        write(DEV_OUT,35) prmfile_onoff(OptimizeGeometry)
        write(DEV_OUT,45) prmfile_onoff(ShowOptimizationProgress)
        write(DEV_OUT,55) prmfile_onoff(KeepOptimizedGeometry)
        return
    end if
    
    if( prmfile_get_logical_by_key(fin,'comb_rules', ApplyCombinationRules)) then
        write(DEV_OUT,20) prmfile_onoff(ApplyCombinationRules)
    else
        write(DEV_OUT,25) prmfile_onoff(ApplyCombinationRules)
    end if
    if( prmfile_get_logical_by_key(fin,'optgeo', OptimizeGeometry)) then
        write(DEV_OUT,30) prmfile_onoff(OptimizeGeometry)
    else
        write(DEV_OUT,35) prmfile_onoff(OptimizeGeometry)
    end if
    if( prmfile_get_logical_by_key(fin,'optprogress', ShowOptimizationProgress)) then
        write(DEV_OUT,40) prmfile_onoff(ShowOptimizationProgress)
    else
        write(DEV_OUT,45) prmfile_onoff(ShowOptimizationProgress)
    end if 
    if( prmfile_get_logical_by_key(fin,'keepoptgeo', KeepOptimizedGeometry)) then
        write(DEV_OUT,50) prmfile_onoff(KeepOptimizedGeometry)
    else
        write(DEV_OUT,55) prmfile_onoff(KeepOptimizedGeometry)
    end if       
      
 10 format('=== [setup] ====================================================================')

 20  format ('Apply combination rules (comb_rules)   = ',a12)
 25  format ('Apply combination rules (comb_rules)   = ',a12,'                  (default)')

 30  format ('Optimize geometry (optgeo)             = ',a12)
 35  format ('Optimize geometry (optgeo)             = ',a12,'                  (default)')
 
 40  format ('Show geo opt progress (optprogress)    = ',a12)
 45  format ('Show geo opt progress (optprogress)    = ',a12,'                  (default)') 

 50  format ('Keep optimized geometry (keepoptgeo)   = ',a12)
 55  format ('Keep optimized geometry (keepoptgeo)   = ',a12,'                  (default)') 
  
end subroutine ffdev_targetset_ctrl_setup

! ------------------------------------------------------------------------------

end module ffdev_targetset_control

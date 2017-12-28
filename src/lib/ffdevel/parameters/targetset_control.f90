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
    character(PRMFILE_MAX_PATH) :: string,topin,key,geoname,sweight
    integer                     :: i,j,alloc_status,minj,probesize,nb_mode,lj_rule
    logical                     :: data_avail,rst,shift2zero
    real(DEVDP)                 :: minenergy,weight
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'{TARGETS}', ':')

    ! open TARGETS group
    if( .not. prmfile_open_group(fin,'TARGETS') ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to open {TARGETS} group!')
    end if

    ! get number of sets
    nsets = prmfile_count_group(fin)
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
        write(DEV_OUT,*)
        write(DEV_OUT,10) i

        ! open set section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get section name!')
        end if

        if( string .ne. 'SET' ) then
            call ffdev_utils_exit(DEV_OUT,1,'Section '''//trim(string)//''' must be [SET]!')
        end if

!----------------------
! topology
!----------------------

        ! get topology name
        if( .not. prmfile_get_string_by_key(fin,'topology',topin) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get topology name (topology) for [SET]!')
        end if
        write(DEV_OUT,12) trim(topin)

        ! get name of resulting topology
        if( .not. prmfile_get_string_by_key(fin,'final',sets(i)%final_name) ) then
            sets(i)%final_name = ''
            write(DEV_OUT,15) trim('-none-')
        else
            write(DEV_OUT,15) trim(sets(i)%final_name)
        end if

        ! load topology and print info
        write(DEV_OUT,*)
        call ffdev_topology_init(sets(i)%top)
        call ffdev_topology_load(sets(i)%top,topin)

        if( .not. prmfile_get_integer_by_key(fin,'probesize',probesize) ) then
            probesize = 0
        end if
        write(DEV_OUT,17) probesize
        call ffdev_topology_switch_to_probe_mode(sets(i)%top,probesize)

        call ffdev_topology_info(sets(i)%top)

        if( prmfile_get_string_by_key(fin,'lj_rule',string) ) then
            write(DEV_OUT,*)
            select case(trim(string))
                case('LB')
                    lj_rule = LJ_RULE_LB
                    write(DEV_OUT,19) 'LB (Lorentz-Berthelot)'
                case('WH')
                    lj_rule = LJ_RULE_WH
                    write(DEV_OUT,19) 'WH (Waldman-Hagler)'
                case('KG')
                    lj_rule = LJ_RULE_KG
                    write(DEV_OUT,19) 'KG (Kong)'
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Unsupported lj_rule in ffdev_targetset_ctrl!')
            end select

            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Original LJ parameters', '#')
            call ffdev_topology_info_types(sets(i)%top,1)

            ! remix parameters
            call ffdev_topology_apply_LJ_combrule(sets(i)%top,lj_rule)

            ! new set of parameters
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Remixed LJ parameters', '#')
             call ffdev_topology_info_types(sets(i)%top,2)
        end if

        if( prmfile_get_string_by_key(fin,'nb_mode',string) ) then
            select case(trim(string))
                case('LJ')
                    nb_mode = NB_MODE_LJ
                case('BP')
                    nb_mode = NB_MODE_BP
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Unsupported nb_mode in ffdev_targetset_ctrl!')
            end select
            call ffdev_topology_switch_nbmode(sets(i)%top,nb_mode)
        end if

        write(DEV_OUT,*)
        select case(sets(i)%top%nb_mode)
            case(NB_MODE_LJ)
                write(DEV_OUT,18) 'LJ'
            case(NB_MODE_BP)
                write(DEV_OUT,18) 'BP'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported nb_mode in ffdev_targetset_ctrl!')
        end select

!----------------------
! points
!----------------------

        sets(i)%offset = 0.0d0
        if( .not. prmfile_get_logical_by_key(fin,'shift2zero',shift2zero) ) then
            shift2zero = .false.
        end if
        write(DEV_OUT,16) prmfile_onoff(shift2zero)

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
            call ffdev_geometry_load_point(sets(i)%geo(j),geoname)
            call ffdev_geometry_info_point(sets(i)%geo(j))
            call ffdev_geometry_check_z(sets(i)%top,sets(i)%geo(j))

            ! do we have data for point
            data_avail = sets(i)%geo(j)%trg_ene_loaded .or. &
                         sets(i)%geo(j)%trg_grd_loaded .or. &
                         sets(i)%geo(j)%trg_hess_loaded

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
            write(DEV_OUT,*)
            if( shift2zero ) then
                call ffdev_geometry_info_point_header_ext()
                do j=1,sets(i)%ngeos
                    if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle
                    sets(i)%geo(j)%trg_energy = sets(i)%geo(j)%trg_energy - minenergy
                    call ffdev_geometry_info_point_ext(sets(i)%geo(j))
                end do
            end if
        end if

        i = i + 1
        rst = prmfile_next_section(fin)

    end do

 10 format('=== [SET] #',I2.2,' ==================================================================')
 12 format('Input topology name (topology)     = ',A)
 15 format('Final topology name (final)        = ',A)
 16 format('Shift minimum to zero (shift2zero) = ',A)
 17 format('Probe size (probesize)             = ',I6)
 18 format('NB mode (nb_mode)                  = ',A)
 19 format('LJ combining rule (lj_rule)        = ',A)

200 format('Number of target points            = ',I6)
300 format('Minimum energy point #',I5.5,' has energy ',F20.4)

end subroutine ffdev_targetset_ctrl

! ------------------------------------------------------------------------------

end module ffdev_targetset_control

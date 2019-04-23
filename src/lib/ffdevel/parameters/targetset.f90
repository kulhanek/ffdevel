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

module ffdev_targetset

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_targetset_init_pts
! at this point parameters should be extracted from topologies
! ==============================================================================

subroutine ffdev_targetset_init_pts

    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_gradient_utils
    use ffdev_hessian_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: i,j
    ! --------------------------------------------------------------------------

    do i=1,nsets
        call ffdev_topology_finalize_setup(sets(i)%top)
        do j=1,sets(i)%ngeos
            ! always allocate gradient as geo optimization can be switched any time
            call ffdev_gradient_allocate(sets(i)%geo(j))
            ! allocate hessian as needed
            if( sets(i)%geo(j)%trg_hess_loaded ) then
                call ffdev_hessian_allocate(sets(i)%geo(j))
                call ffdev_gradient_allocate(sets(i)%geo(j))                
            end if
        end do
    end do

end subroutine ffdev_targetset_init_pts

! ==============================================================================
! subroutine ffdev_targetset_calc_all
! NOTE - it requires a program to be loaded
! ==============================================================================

subroutine ffdev_targetset_calc_all

    use ffdev_targetset_dat
    use ffdev_energy
    use ffdev_gradient
    use ffdev_hessian
    use ffdev_geoopt
    use ffdev_errors_dat
    use ffdev_hessian_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: i,j,k
    ! --------------------------------------------------------------------------

    ! apply combination rules
    if( ApplyCombinationRules ) then
        do i=1,nsets
            if( sets(i)%top%assumed_comb_rules .eq. COMB_RULE_IN ) then
                call ffdev_topology_apply_NB_comb_rules(sets(i)%top,NBCombRules)
            else
                call ffdev_topology_apply_NB_comb_rules(sets(i)%top,sets(i)%top%assumed_comb_rules)
            end if
        end do
        call ffdev_parameters_reinit_nbparams()
    end if    
    
    ! optimize geometry
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( (sets(i)%optgeo .or. GlbOptGeometryEnabled) .and. (.not. GlbOptGeometryDisabled) &
                .and. sets(i)%geo(j)%trg_crd_loaded ) then
                if( .not. (sets(i)%keepoptgeo .or. GlbKeepOptGeometry) ) then
                    sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                end if
                if( GlbShowOptProgress ) then
                    call ffdev_geoopt_run(DEV_OUT,sets(i)%top,sets(i)%geo(j))
                else
                    call ffdev_geoopt_run(DEV_NULL,sets(i)%top,sets(i)%geo(j))
                end if
                sets(i)%geo(j)%trg_crd_optimized = .true.
            else
                sets(i)%geo(j)%trg_crd_optimized = .false.
            end if
        end do
    end do
    
    ! calculate all energies, gradients, hessians    
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_hess_loaded .and. (errors_calc_hess .or. errors_calc_freq) ) then
                ! calculate hessian
                call ffdev_hessian_num_all(sets(i)%top,sets(i)%geo(j))
                ! calculate frequencies
                if( sets(i)%geo(j)%trg_freq_loaded .and. errors_calc_freq ) then
                    call ffdev_hessian_allocate_freq(sets(i)%geo(j))
                    call ffdev_hessian_calc_freqs(sets(i)%geo(j))
                end if
            else if( sets(i)%geo(j)%trg_grd_loaded .and. errors_calc_grad ) then
                call ffdev_gradient_all(sets(i)%top,sets(i)%geo(j))
            else if( sets(i)%geo(j)%trg_ene_loaded .and. errors_calc_ene ) then
                call ffdev_energy_all(sets(i)%top,sets(i)%geo(j))
            end if
        end do
    end do
    
   ! update energies 
    do i=1,nsets
        if( sets(i)%nrefs .gt. 0 ) then
            do j=1,sets(i)%ngeos
                do k=1,sets(i)%nrefs
                    sets(i)%geo(j)%bond_ene     = sets(i)%geo(j)%bond_ene - sets( sets(i)%refs(k) )%geo(1)%bond_ene
                    sets(i)%geo(j)%angle_ene    = sets(i)%geo(j)%angle_ene - sets( sets(i)%refs(k) )%geo(1)%angle_ene
                    sets(i)%geo(j)%dih_ene      = sets(i)%geo(j)%dih_ene - sets( sets(i)%refs(k) )%geo(1)%dih_ene
                    sets(i)%geo(j)%impropr_ene  = sets(i)%geo(j)%impropr_ene - sets( sets(i)%refs(k) )%geo(1)%impropr_ene
                    sets(i)%geo(j)%ele14_ene    = sets(i)%geo(j)%ele14_ene - sets( sets(i)%refs(k) )%geo(1)%ele14_ene
                    sets(i)%geo(j)%nb14_ene     = sets(i)%geo(j)%nb14_ene - sets( sets(i)%refs(k) )%geo(1)%nb14_ene
                    sets(i)%geo(j)%ele_ene      = sets(i)%geo(j)%ele_ene - sets( sets(i)%refs(k) )%geo(1)%ele_ene
                    sets(i)%geo(j)%nb_ene       = sets(i)%geo(j)%nb_ene - sets( sets(i)%refs(k) )%geo(1)%nb_ene
                    sets(i)%geo(j)%total_ene    = sets(i)%geo(j)%total_ene - sets( sets(i)%refs(k) )%geo(1)%total_ene  
                end do
            end do
        end if
    end do          

end subroutine ffdev_targetset_calc_all

! ==============================================================================
! subroutine ffdev_parameters_reinit_nbparams
! ==============================================================================

subroutine ffdev_parameters_reinit_nbparams()

    use ffdev_parameters_dat
    use ffdev_targetset_dat

    implicit none
    integer     :: i,j
    ! --------------------------------------------------------------------------

    ! update parameter values
    do i=1,nparams
        select case(params(i)%realm)
            case(REALM_VDW_EPS,REALM_VDW_R0,REALM_VDW_ALPHA,REALM_VDW_A,REALM_VDW_B,REALM_VDW_C)
                do j=1,nsets
                    if( params(i)%ids(j) .eq. 0 ) cycle
                    select case(params(i)%realm)
                        case(REALM_VDW_EPS)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%eps
                        case(REALM_VDW_R0)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%r0
                        case(REALM_VDW_ALPHA)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%alpha
                        case(REALM_VDW_A)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%A
                        case(REALM_VDW_B)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%B
                        case(REALM_VDW_C)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%C
                    end select
                end do
        end select
    end do

end subroutine ffdev_parameters_reinit_nbparams

! ==============================================================================
! subroutine ffdev_targetset_save_final_stops
! ==============================================================================

subroutine ffdev_targetset_save_final_stops

    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    do i=1,nsets
        if( len_trim(sets(i)%final_stop) .ne. 0 ) then
            write(DEV_OUT,10) i,trim(sets(i)%final_stop)
            call ffdev_topology_save(sets(i)%top,sets(i)%final_stop)
        end if
    end do

 10 format('Final topology name for set #',I2.2,' = ',A)

end subroutine ffdev_targetset_save_final_stops

! ==============================================================================
! subroutine ffdev_targetset_save_final_pts
! ==============================================================================

subroutine ffdev_targetset_save_final_pts

    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer                 :: i,j
    character(len=MAX_PATH) :: sname
    ! --------------------------------------------------------------------------

    write(DEV_OUT,10)
    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! save geometry if requested
            if( sets(i)%savegeo .and. sets(i)%geo(j)%trg_crd_optimized ) then
                write(sname,20) i,j
                sname = trim(SavePointsPath)//'/'//trim(sname)
                write(DEV_OUT,30) trim(sname)
                call ffdev_geometry_save_point(sets(i)%geo(j),sname)
            end if
        end do
    end do

 10 format('Saving final point geometries ...')
 20 format('S',I2.2,'P',I6.6,'.pst')
 30 format(7X,A)

end subroutine ffdev_targetset_save_final_pts

! ==============================================================================
! subroutine ffdev_targetset_save_initial_drvs
! ==============================================================================

subroutine ffdev_targetset_save_initial_drvs

    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    ! write files
    do i=1,nsets
        if( len_trim(sets(i)%initial_drvene) .ne. 0 ) then
            write(DEV_OUT,10) i,trim(sets(i)%initial_drvene)
            call ffdev_targetset_save_drvene(sets(i),sets(i)%initial_drvene)
        end if
        if( len_trim(sets(i)%initial_drvxyz) .ne. 0 ) then
            write(DEV_OUT,20) i,trim(sets(i)%initial_drvxyz)
            call ffdev_targetset_save_drvxyz(sets(i),sets(i)%initial_drvxyz)
        end if
    end do

 10 format('Saving initial driving data (drvene) for set #',I2.2,' = ',A)
 20 format('Saving initial driving data (drvxyz) for set #',I2.2,' = ',A)

end subroutine ffdev_targetset_save_initial_drvs

! ==============================================================================
! subroutine ffdev_targetset_save_final_drvs
! ==============================================================================

subroutine ffdev_targetset_save_final_drvs

    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    do i=1,nsets
        if( len_trim(sets(i)%final_drvene) .ne. 0 ) then
            write(DEV_OUT,10) i,trim(sets(i)%final_drvene)
            call ffdev_targetset_save_drvene(sets(i),sets(i)%final_drvene)
        end if
        if( len_trim(sets(i)%final_drvxyz) .ne. 0 ) then
            write(DEV_OUT,20) i,trim(sets(i)%final_drvxyz)
            call ffdev_targetset_save_drvxyz(sets(i),sets(i)%final_drvxyz)
        end if
    end do

    10 format('Saving final driving data (drvene) for set #',I2.2,' = ',A)
    20 format('Saving final driving data (drvxyz) for set #',I2.2,' = ',A)

end subroutine ffdev_targetset_save_final_drvs

! ==============================================================================
! subroutine ffdev_targetset_save_drvene
! ==============================================================================

subroutine ffdev_targetset_save_drvene(set,name)

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TARGETSET)         :: set
    character(len=MAX_PATH) :: name
    ! --------------------------------------------
    integer                 :: i,j,nrst
    real(DEVDP)             :: err
    ! --------------------------------------------------------------------------

    if( set%ngeos .eq. 0 ) return ! no geometries

    nrst = set%geo(1)%nrst

    call ffdev_utils_open(DEV_PROFILE,name,'U')

    ! write header
    write(DEV_PROFILE,10,ADVANCE='NO')
    do i=1,set%geo(1)%nrst
        write(DEV_PROFILE,11,ADVANCE='NO') i
    end do
    write(DEV_PROFILE,12)

    write(DEV_PROFILE,20,ADVANCE='NO')
    do i=1,set%geo(1)%nrst
        write(DEV_PROFILE,21,ADVANCE='NO')
    end do
    write(DEV_PROFILE,22)

    do i=1,set%ngeos
        write(DEV_PROFILE,30,ADVANCE='NO') i
        if( nrst .eq. set%geo(i)%nrst ) then
            call ffdev_geometry_get_rst_values(set%geo(i))
            do j=1,set%geo(i)%nrst
                select case(trim(set%geo(i)%rst(j)%cvtype))
                    case('B')
                        write(DEV_PROFILE,31,ADVANCE='NO') set%geo(i)%rst(j)%value
                    case('A','D')
                        write(DEV_PROFILE,32,ADVANCE='NO') set%geo(i)%rst(j)%value*DEV_R2D
                    case default
                        call ffdev_utils_exit(DEV_OUT,1,'Unsupported RST in ffdev_geometry_copy!')
                end select
            end do
        else
            do j=1,nrst
                write(DEV_PROFILE,21,ADVANCE='NO')
            end do
        end if
        err = set%geo(i)%total_ene - set%offset - set%geo(i)%trg_energy
        write(DEV_PROFILE,33) set%geo(i)%trg_energy, set%geo(i)%total_ene - set%offset, err
    end do

    close(DEV_PROFILE)

 10 format('# IDX')
 11 format('       CV',I2.2)
 12 format('   E(TRG,1)    E(MM,2)     E(2-1)')
 20 format('# ---')
 21 format(' ----------')
 22 format(' ---------- ---------- ----------')
 30 format(I5)
 31 format(1X,F10.4)
 32 format(1X,F10.1)
 33 format(1X,F10.2,1X,F10.2,1X,F10.2)

end subroutine ffdev_targetset_save_drvene

! ==============================================================================
! subroutine ffdev_targetset_save_drvxyz
! ==============================================================================

subroutine ffdev_targetset_save_drvxyz(set,name)

    use ffdev_targetset_dat
    use ffdev_geometry
    use smf_xyzfile
    use smf_xyzfile_type

    implicit none
    type(TARGETSET)         :: set
    character(len=MAX_PATH) :: name
    ! --------------------------------------------
    integer                 :: i
    type(XYZFILE_TYPE)      :: xyz
    ! --------------------------------------------------------------------------

    if( set%ngeos .eq. 0 ) return ! no geometries

    call init_xyz(xyz)
    call allocate_xyz(xyz,set%top%natoms)
    call open_xyz(DEV_XYZ,name,xyz,'UNKNOWN')

    do i=1,set%ngeos
        call ffdev_geometry_save_xyz_snapshot(set%geo(i),xyz)
        call write_xyz(DEV_XYZ,xyz)
    end do

    call close_xyz(DEV_XYZ,xyz)

end subroutine ffdev_targetset_save_drvxyz

! ------------------------------------------------------------------------------

end module ffdev_targetset

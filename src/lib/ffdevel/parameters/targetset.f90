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
use ffdev_variables

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
    use ffdev_utils
    use ffdev_timers
    use ffdev_geometry
    use ffdev_parallel_dat

    implicit none
    integer         :: i,j,k
    integer,save    :: evalcounter = 0      ! global counter
    integer         :: georeseted
    ! --------------------------------------------------------------------------

    call ffdev_timers_start_timer(FFDEV_TARGETSET_ALL_TIMER)

    evalcounter = evalcounter + 1

    georeseted = 0

    ! PARALLEL !!!!! ensure that NB parameters are updated before the parallel section !!!!!
    do i=1,nsets
        if( sets(i)%top%nb_params_update ) then
            call ffdev_topology_update_nb_params(sets(i)%top)
            call ffdev_topology_update_nb_pairs(sets(i)%top)
        end if
    end do

    call ffdev_timers_start_timer(FFDEV_GEOOPT_TIMER)
    if( ParallelGeoOptBlock ) then
        !$omp parallel
        !$omp single
        ! optimize geometry
        do i=1,nsets
            if( sets(i)%top%probe_size .gt. 0 ) cycle   ! do not optimize probed structures
            if( sets(i)%top%natoms .le. 1 ) cycle       ! do not optimize systems with one atom
            do j=1,sets(i)%ngeos
                sets(i)%geo(j)%trg_crd_optimized = .false.
                if( (sets(i)%optgeo .or. GlbOptGeometryEnabled) .and. (.not. GlbOptGeometryDisabled) &
                    .and. sets(i)%geo(j)%trg_crd_loaded ) then
                    if( .not. (sets(i)%keepoptgeo .or. GlbKeepOptGeometry) ) then
                        sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                    else
                        if( ResetKeptOptGeoAt .gt. 0 ) then
                            if( mod(evalcounter,ResetKeptOptGeoAt) .eq. 0 ) then
                                sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                                georeseted = georeseted + 1
                            end if
                        end if
                    end if
                    !$omp task
                    ! no output is allowed in parallel mode
                    call ffdev_geoopt_run(0,sets(i)%top,sets(i)%geo(j))
                    sets(i)%geo(j)%trg_crd_optimized = .true.
                    !$omp end task
                end if
            end do
        end do
        !$omp end single
        !$omp end parallel
    else
        ! optimize geometry
        do i=1,nsets
            if( sets(i)%top%probe_size .gt. 0 ) cycle   ! do not optimize probed structures
            if( sets(i)%top%natoms .le. 1 ) cycle       ! do not optimize systems with one atom
            do j=1,sets(i)%ngeos
                sets(i)%geo(j)%trg_crd_optimized = .false.
                if( (sets(i)%optgeo .or. GlbOptGeometryEnabled) .and. (.not. GlbOptGeometryDisabled) &
                    .and. sets(i)%geo(j)%trg_crd_loaded ) then
                    if( .not. (sets(i)%keepoptgeo .or. GlbKeepOptGeometry) ) then
                        sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                    else
                        if( ResetKeptOptGeoAt .gt. 0 ) then
                            if( mod(evalcounter,ResetKeptOptGeoAt) .eq. 0 ) then
                                sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                                georeseted = georeseted + 1
                            end if
                        end if
                    end if
                    if( GlbShowOptProgress ) then
                        call ffdev_geoopt_run(DEV_OUT,sets(i)%top,sets(i)%geo(j))
                    else
                        call ffdev_geoopt_run(DEV_NULL,sets(i)%top,sets(i)%geo(j))
                    end if
                    sets(i)%geo(j)%trg_crd_optimized = .true.
                end if
            end do
        end do
    end if
    call ffdev_timers_stop_timer(FFDEV_GEOOPT_TIMER)

    if( georeseted .gt. 0 ) then
        write(DEV_OUT,*) '>>> Geometries (',georeseted,') reset at ', evalcounter
    end if

    call ffdev_timers_start_timer(FFDEV_PROPS_TIMER)
    if( ParallelErrorBlock ) then
        !$omp parallel
        !$omp single
        ! calculate all energies, gradients, hessians
        do i=1,nsets
            do j=1,sets(i)%ngeos
                if( sets(i)%geo(j)%trg_hess_loaded .and. errors_calc_hess ) then
                    ! calculate hessian
                    !$omp task
                    call ffdev_hessian_num_all(sets(i)%top,sets(i)%geo(j))
                    !$omp end task
                ! for zerograd error we do not need trg_grd_loaded
                ! else if( sets(i)%geo(j)%trg_grd_loaded .and. errors_calc_grad ) then
                else if( errors_calc_grad ) then
                    !$omp task
                    call ffdev_gradient_all(sets(i)%top,sets(i)%geo(j))
                    !$omp end task
                else if( (sets(i)%geo(j)%trg_ene_loaded .or. sets(i)%geo(j)%trg_probe_ene_loaded) .and. errors_calc_ene ) then
                    !$omp task
                    call ffdev_energy_all(sets(i)%top,sets(i)%geo(j))
                    !$omp end task
                end if
            end do
        end do

        ! calculate SAPT
        do i=1,nsets
            if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle
            do j=1,sets(i)%ngeos
                if( .not. ( sets(i)%geo(j)%trg_sapt_loaded .and. errors_calc_sapt ) ) cycle
                !$omp task
                call ffdev_energy_sapt(sets(i)%top,sets(i)%geo(j))
                !$omp end task
            end do
        end do
        !$omp end single
        !$omp end parallel
    else
        ! calculate all energies, gradients, hessians
        do i=1,nsets
            do j=1,sets(i)%ngeos
                if( sets(i)%geo(j)%trg_hess_loaded .and. errors_calc_hess ) then
                    ! calculate hessian
                    call ffdev_hessian_num_all(sets(i)%top,sets(i)%geo(j))
                ! for zerograd error we do not need trg_grd_loaded
                ! else if( sets(i)%geo(j)%trg_grd_loaded .and. errors_calc_grad ) then
                else if( errors_calc_grad ) then
                    call ffdev_gradient_all(sets(i)%top,sets(i)%geo(j))
                else if( sets(i)%geo(j)%trg_ene_loaded .and. errors_calc_ene ) then
                    call ffdev_energy_all(sets(i)%top,sets(i)%geo(j))
                end if
            end do
        end do

        ! calculate SAPT
        do i=1,nsets
            if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle
            do j=1,sets(i)%ngeos
                if( .not. ( (sets(i)%geo(j)%trg_sapt_loaded .or. sets(i)%geo(j)%trg_probe_ene_loaded) &
                            .and. errors_calc_sapt ) ) cycle
                call ffdev_energy_sapt(sets(i)%top,sets(i)%geo(j))
            end do
        end do
    end if

    call ffdev_timers_stop_timer(FFDEV_PROPS_TIMER)

   ! update energies
    do i=1,nsets
        if( (sets(i)%nrefs .gt. 0) .and. (sets(i)%top%probe_size .le. 0) ) then
            do j=1,sets(i)%ngeos
                do k=1,sets(i)%nrefs
                    sets(i)%geo(j)%bond_ene     = sets(i)%geo(j)%bond_ene    - sets( sets(i)%refs(k) )%geo(1)%bond_ene
                    sets(i)%geo(j)%angle_ene    = sets(i)%geo(j)%angle_ene   - sets( sets(i)%refs(k) )%geo(1)%angle_ene
                    sets(i)%geo(j)%dih_ene      = sets(i)%geo(j)%dih_ene     - sets( sets(i)%refs(k) )%geo(1)%dih_ene
                    sets(i)%geo(j)%impropr_ene  = sets(i)%geo(j)%impropr_ene - sets( sets(i)%refs(k) )%geo(1)%impropr_ene

                    sets(i)%geo(j)%ele14_ene    = sets(i)%geo(j)%ele14_ene - sets( sets(i)%refs(k) )%geo(1)%ele14_ene
                    sets(i)%geo(j)%rep14_ene    = sets(i)%geo(j)%rep14_ene - sets( sets(i)%refs(k) )%geo(1)%rep14_ene
                    sets(i)%geo(j)%dis14_ene    = sets(i)%geo(j)%dis14_ene - sets( sets(i)%refs(k) )%geo(1)%dis14_ene

                    sets(i)%geo(j)%ele_ene      = sets(i)%geo(j)%ele_ene - sets( sets(i)%refs(k) )%geo(1)%ele_ene
                    sets(i)%geo(j)%pen_ene      = sets(i)%geo(j)%pen_ene - sets( sets(i)%refs(k) )%geo(1)%pen_ene
                    sets(i)%geo(j)%ind_ene      = sets(i)%geo(j)%ind_ene - sets( sets(i)%refs(k) )%geo(1)%ind_ene
                    sets(i)%geo(j)%rep_ene      = sets(i)%geo(j)%rep_ene - sets( sets(i)%refs(k) )%geo(1)%rep_ene
                    sets(i)%geo(j)%dis_ene      = sets(i)%geo(j)%dis_ene - sets( sets(i)%refs(k) )%geo(1)%dis_ene

                    sets(i)%geo(j)%total_ene    = sets(i)%geo(j)%total_ene - sets( sets(i)%refs(k) )%geo(1)%total_ene
                end do
            end do
        end if
    end do

    call ffdev_timers_stop_timer(FFDEV_TARGETSET_ALL_TIMER)

end subroutine ffdev_targetset_calc_all

! ==============================================================================
! subroutine ffdev_targetset_reinit_nbparams
! ==============================================================================

subroutine ffdev_targetset_reinit_nbparams()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer     :: i,j
    ! --------------------------------------------------------------------------

    ! update parameter values
    do i=1,nparams
        do j=1,nsets
            if( params(i)%ids(j) .eq. 0 ) cycle
            select case(params(i)%realm)
                case(REALM_BOND_D0,REALM_BOND_K, &
                     REALM_ANGLE_A0,REALM_ANGLE_K,REALM_DIH_V,REALM_DIH_G, &
                     REALM_DIH_SCEE,REALM_DIH_SCNB,REALM_IMPR_V, &
                     REALM_IMPR_G,REALM_DIH_C,REALM_PAC)
                    ! nothing to be here
                case(REALM_VDW_EPS)
                    params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%eps
                case(REALM_VDW_R0)
                    params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%r0
                case(REALM_VDW_ALPHA)
                    params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%alpha
                case(REALM_VDW_PA)
                    params(i)%value = sets(j)%top%atom_types(params(i)%ids(j))%PA
                case(REALM_VDW_PB)
                    params(i)%value = sets(j)%top%atom_types(params(i)%ids(j))%PB
                case(REALM_VDW_RC)
                    params(i)%value = sets(j)%top%atom_types(params(i)%ids(j))%RC
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_targetset_reinit_nbparams!')
            end select
        end do
    end do

end subroutine ffdev_targetset_reinit_nbparams

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
    integer                 :: i,j,nsaved
    character(len=MAX_PATH) :: sname
    ! --------------------------------------------------------------------------

    if( UseOptGeometry ) then
        write(DEV_OUT,10)
    else
        write(DEV_OUT,15)
    end if

    nsaved = 0
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( UseOptGeometry ) then
                ! save geometry if requested - optimized geo
                if( sets(i)%savepts .and. sets(i)%geo(j)%trg_crd_optimized ) then
                    if( KeepPtsNames ) then
                        write(sname,25) trim(sets(i)%geo(j)%name)
                    else
                        write(sname,20) i,j
                    end if
                    sname = trim(SavePtsPath)//'/'//trim(sname)
                    write(DEV_OUT,30) trim(sname)
                    call ffdev_geometry_save_point(sets(i)%geo(j),sname,.false.)
                    nsaved = nsaved + 1
                end if
            else
                ! save geometry if requested - target geo
                if( sets(i)%savepts .and. sets(i)%geo(j)%trg_crd_loaded ) then
                    if( KeepPtsNames ) then
                        write(sname,25) trim(sets(i)%geo(j)%name)
                    else
                        write(sname,20) i,j
                    end if

                    sname = trim(SavePtsPath)//'/'//trim(sname)
                    write(DEV_OUT,30) trim(sname)
                    call ffdev_geometry_save_point(sets(i)%geo(j),sname,.true.)
                    nsaved = nsaved + 1
                end if
            end if
        end do
    end do

    if( nsaved .eq. 0 ) then
        write(DEV_OUT,16)
    end if

 10 format('Saving final point geometries ... using MM optimized coordinates')
 15 format('Saving final point geometries ... using target coordinates')
 16 format('       >>> INFO: No request to save pts geometries ...')
 20 format('S',I2.2,'P',I6.6,'.pst')
 25 format(A)
 30 format(7X,A)

end subroutine ffdev_targetset_save_final_pts

! ==============================================================================
! subroutine ffdev_targetset_save_final_xyzr
! ==============================================================================

subroutine ffdev_targetset_save_final_xyzr

    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer                 :: i,j,nsaved
    character(len=MAX_PATH) :: sname
    ! --------------------------------------------------------------------------

    if( UseOptGeometry ) then
        write(DEV_OUT,10)
    else
        write(DEV_OUT,15)
    end if

    nsaved = 0
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( UseOptGeometry ) then
                ! save geometry if requested - optimized geo
                if( sets(i)%savexyzr .and. sets(i)%geo(j)%trg_crd_optimized ) then
                    if( KeepPtsNames ) then
                        write(sname,25) trim(sets(i)%geo(j)%name)
                    else
                        write(sname,20) i,j
                    end if
                    sname = trim(SaveXYZRPath)//'/'//trim(sname)
                    write(DEV_OUT,30) trim(sname)
                    call ffdev_geometry_save_xyzr(sets(i)%top,sets(i)%geo(j),sname,.false.)
                end if
            else
                ! save geometry if requested - target geo
                if( sets(i)%savexyzr .and. sets(i)%geo(j)%trg_crd_loaded ) then
                    if( KeepPtsNames ) then
                        write(sname,25) trim(sets(i)%geo(j)%name)
                    else
                        write(sname,20) i,j
                    end if
                    sname = trim(SaveXYZRPath)//'/'//trim(sname)
                    write(DEV_OUT,30) trim(sname)
                    call ffdev_geometry_save_xyzr(sets(i)%top,sets(i)%geo(j),sname,.true.)
                end if
            end if
        end do
    end do

    if( nsaved .eq. 0 ) then
        write(DEV_OUT,16)
    end if

 10 format('Saving final xyzr geometries  ... using MM optimized coordinates')
 15 format('Saving final xyzr geometries  ... using target coordinates')
 16 format('       >>> INFO: No request to save xyzr geometries ...')
 20 format('S',I2.2,'P',I6.6,'.xyzr')
 25 format(A,'.xyzr')
 30 format(7X,A)

end subroutine ffdev_targetset_save_final_xyzr

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
    write(DEV_PROFILE,12,ADVANCE='NO')
    write(DEV_PROFILE,13)

    write(DEV_PROFILE,20,ADVANCE='NO')
    do i=1,set%geo(1)%nrst
        write(DEV_PROFILE,21,ADVANCE='NO')
    end do
    write(DEV_PROFILE,22,ADVANCE='NO')
    write(DEV_PROFILE,23)

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
                        call ffdev_utils_exit(DEV_ERR,1,'Unsupported RST in ffdev_geometry_copy!')
                end select
            end do
        else
            do j=1,nrst
                write(DEV_PROFILE,21,ADVANCE='NO')
            end do
        end if
        err = set%geo(i)%total_ene - set%geo(i)%trg_energy
        write(DEV_PROFILE,33,ADVANCE='NO') set%geo(i)%trg_energy, set%geo(i)%total_ene, err
        write(DEV_PROFILE,34) set%geo(i)%bond_ene, set%geo(i)%angle_ene, set%geo(i)%dih_ene, &
                              set%geo(i)%impropr_ene, set%geo(i)%dih_ene + set%geo(i)%impropr_ene, &
                              set%geo(i)%ele_ene, set%geo(i)%ele14_ene, set%geo(i)%ele_ene + set%geo(i)%ele14_ene, &
                              set%geo(i)%rep_ene, set%geo(i)%rep14_ene, set%geo(i)%rep_ene + set%geo(i)%rep14_ene, &
                              set%geo(i)%dis_ene, set%geo(i)%dis14_ene, set%geo(i)%dis_ene + set%geo(i)%dis14_ene, &
                              set%geo(i)%bond_ene + set%geo(i)%angle_ene + set%geo(i)%dih_ene + set%geo(i)%impropr_ene, &
                              set%geo(i)%ele_ene + set%geo(i)%ele14_ene + set%geo(i)%rep_ene + set%geo(i)%rep14_ene + &
                              set%geo(i)%dis_ene + set%geo(i)%dis14_ene

    end do

    close(DEV_PROFILE)

 10 format('# IDX')
 11 format('       CV',I2.2)
 12 format('   E(TRG,1)    E(MM,2)     E(2-1)')
 13 format('     Ebonds    Eangles      Etors      Eimps  Edih(t+i)        Eel      E14el    Etotele' &
           '       Erep     E14rep    Etotrep      Edisp    E14disp   Etotdisp        Ebn        Enb')
 20 format('# ---')
 21 format(' ----------')
 22 format(' ---------- ---------- ----------')
 23 format(' ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------' &
           ' ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------')
 30 format(I5)
 31 format(1X,F10.4)
 32 format(1X,F10.1)
 33 format(1X,F10.2,1X,F10.2,1X,F10.2)
 34 format(1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2, &
           1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2,1X,F10.2)

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

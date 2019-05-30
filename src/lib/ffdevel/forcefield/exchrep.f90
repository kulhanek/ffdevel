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

module ffdev_exchrep

use ffdev_geometry_dat
use ffdev_constants

private get_overlap2
private get_overlap3A
private get_overlap3B
private get_overlap4A
private get_overlap4B
private get_lda2
private get_lda3A
private get_lda3B
private get_lda4A
private get_lda4B

private grid_integrate

contains

! ==============================================================================
! subroutine ffdevel_exchrep_gen_numgrid_cache
! ==============================================================================

subroutine ffdevel_exchrep_gen_numgrid_cache(cache)

    use, intrinsic :: iso_c_binding, only: c_ptr
    use numgrid
    use ffdev_timers
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(GRID_CACHE)    :: cache
    integer             :: z1
    integer             :: z2
    ! --------------------------------------------
    integer             :: xidx
    type(c_ptr)         :: gctx
    real(DEVDP)         :: radial_precision
    integer             :: min_num_angular_points, max_num_angular_points
    real(DEVDP)         :: alpha_max(2)
    integer             :: max_l_quantum_number(2)
    real(DEVDP)         :: alpha_min(5,2)
    integer             :: az(2)
    real(DEVDP)         :: px(2),py(2),pz(2)
    integer             :: nats, cidx, npts, ipts, alloc_stat
    ! --------------------------------------------------------------------------

    ! numgrid setup
    radial_precision = 1.0d-12
    min_num_angular_points = 86
    max_num_angular_points = 302

    ! ne/ne - def2-QZVPD
    alpha_max(1) =  160676.2795500
    alpha_max(2) =  160676.2795500
    max_l_quantum_number(1) = 4
    max_l_quantum_number(2) = 4

    alpha_min(1,1) = 0.30548672205
    alpha_min(2,1) = 0.16889708453
    alpha_min(3,1) = 0.74700000
    alpha_min(4,1) = 1.52400000
    alpha_min(5,1) = 2.98300000

    alpha_min(1,2) = 0.30548672205
    alpha_min(2,2) = 0.16889708453
    alpha_min(3,2) = 0.74700000
    alpha_min(4,2) = 1.52400000
    alpha_min(5,2) = 2.98300000

    az(1) = 10 ! Ne
    az(2) = 10 ! Ne

    ! coordinates are in a.u.
    px(1) = 0.0d0
    py(1) = 0.0d0
    pz(1) = 0.0d0

    px(2) = cache%r * DEV_A2AU
    py(2) = 0.0d0
    pz(2) = 0.0d0

    ! we have two centers
    nats = 2

    ! generate grid
    call ffdev_timers_start_timer(FFDEV_POT_NB_GRID)

        cidx = 1

        gctx = numgrid_new_atom_grid(radial_precision,       &
                                   min_num_angular_points, &
                                   max_num_angular_points, &
                                   az(cidx),      &
                                   alpha_max(cidx),              &
                                   max_l_quantum_number(cidx),   &
                                   alpha_min(:,cidx))

        npts = numgrid_get_num_grid_points(gctx)

        cache%npts1 = npts
        allocate(cache%grid_1(npts,4), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_exchrep_gen_cache!')
        end if

        call numgrid_get_grid(gctx, nats, cidx-1,   &
                         px, py, pz, az,            &
                         cache%grid_1(:,1), cache%grid_1(:,2), cache%grid_1(:,3), cache%grid_1(:,4))

        call numgrid_free_atom_grid(gctx)

! --------------------------

        cidx = 2

        gctx = numgrid_new_atom_grid(radial_precision,       &
                                   min_num_angular_points, &
                                   max_num_angular_points, &
                                   az(cidx),      &
                                   alpha_max(cidx),              &
                                   max_l_quantum_number(cidx),   &
                                   alpha_min(:,cidx))

        npts = numgrid_get_num_grid_points(gctx)

        cache%npts2 = npts
        allocate(cache%grid_2(npts,4), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_exchrep_gen_cache!')
        end if

        call numgrid_get_grid(gctx, nats, cidx-1,   &
                         px, py, pz, az,            &
                         cache%grid_2(:,1), cache%grid_2(:,2), cache%grid_2(:,3), cache%grid_2(:,4))

        call numgrid_free_atom_grid(gctx)

    call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

end subroutine ffdevel_exchrep_gen_numgrid_cache

! ==============================================================================
! function ffdevel_exchrep_ene_numgrid_nocache
! ==============================================================================

real(DEVDP) function ffdevel_exchrep_ene_numgrid_nocache(mode,r,z1,pa1,pb1,pc1,pd1,z2,pa2,pb2,pc2,pd2)

    use numgrid
    use, intrinsic :: iso_c_binding, only: c_ptr
    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    integer     :: mode
    real(DEVDP) :: r
    integer     :: z1
    real(DEVDP) :: pa1,pb1,pc1,pd1
    integer     :: z2
    real(DEVDP) :: pa2,pb2,pc2,pd2
    ! --------------------------------------------
    integer                     :: xidx
    type(c_ptr)                 :: gctx
    real(DEVDP)                 :: radial_precision
    integer                     :: min_num_angular_points, max_num_angular_points
    real(DEVDP)                 :: alpha_max(2)
    integer                     :: max_l_quantum_number(2)
    real(DEVDP)                 :: alpha_min(5,2)
    integer                     :: az(2)
    real(DEVDP)                 :: px(2),py(2),pz(2)
    real(DEVDP), allocatable    :: grid_x(:)
    real(DEVDP), allocatable    :: grid_y(:)
    real(DEVDP), allocatable    :: grid_z(:)
    real(DEVDP), allocatable    :: grid_w(:)
    integer                     :: nats, cidx, npts, ipts, alloc_stat
    real(DEVDP)                 :: rsum
    real(DEVDP)                 :: epa1,epa2,epb1,epb2,epc1,epc2
    ! --------------------------------------------------------------------------

    ffdevel_exchrep_ene_numgrid_nocache = 0.0d0

    if( pb2 .eq. 0.0d0 ) return
    if( pb2 .eq. 0.0d0 ) return

    ! numgrid setup
    radial_precision = 1.0d-12
    min_num_angular_points = 86
    max_num_angular_points = 302

    ! ne/ne - def2-QZVPD
    alpha_max(1) =  160676.2795500
    alpha_max(2) =  160676.2795500
    max_l_quantum_number(1) = 4
    max_l_quantum_number(2) = 4

    alpha_min(1,1) = 0.30548672205
    alpha_min(2,1) = 0.16889708453
    alpha_min(3,1) = 0.74700000
    alpha_min(4,1) = 1.52400000
    alpha_min(5,1) = 2.98300000

    alpha_min(1,2) = 0.30548672205
    alpha_min(2,2) = 0.16889708453
    alpha_min(3,2) = 0.74700000
    alpha_min(4,2) = 1.52400000
    alpha_min(5,2) = 2.98300000

    az(1) = 10
    az(2) = 10

    ! coordinates are in a.u.
    px(1) = 0.0d0
    py(1) = 0.0d0
    pz(1) = 0.0d0

    px(2) = r * DEV_A2AU
    py(2) = 0.0d0
    pz(2) = 0.0d0

    epa1 = exp(pa1)
    epa2 = exp(pa2)

    select case(mode)
        case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_DENS3B,NB_MODE_PAULI_DENS4B)
            epb1 = 2.0d0*pb1
            epc1 = pc1
            epb2 = 2.0d0*pb2
            epc2 = pc2
        case(NB_MODE_PAULI_DENS3A,NB_MODE_PAULI_DENS4A)
            epb1 = 2.0d0*pb1
            epc1 = 2.0d0*pc1
            epb2 = 2.0d0*pb2
            epc2 = 2.0d0*pc2
        case(NB_MODE_PAULI_WAVE2,NB_MODE_PAULI_WAVE3A,NB_MODE_PAULI_WAVE3B,NB_MODE_PAULI_WAVE4A,NB_MODE_PAULI_WAVE4B)
            epb1 = pb1
            epc1 = pc1
            epb2 = pb2
            epc2 = pc2
        case(NB_MODE_PAULI_LDA2,NB_MODE_PAULI_LDA3B,NB_MODE_PAULI_LDA4B)
            epb1 = 2.0d0*pb1
            epc1 = pc1
            epb2 = 2.0d0*pb2
            epc2 = pc2
        case(NB_MODE_PAULI_LDA3A,NB_MODE_PAULI_LDA4A)
            epb1 = 2.0d0*pb1
            epc1 = 2.0d0*pc1
            epb2 = 2.0d0*pb2
            epc2 = 2.0d0*pc2
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_numgrid_nocache!')
    end select

    rsum =  0.0d0

    ! we have two centers
    nats = 2
    do cidx=1,nats

        ! generate grid
        call ffdev_timers_start_timer(FFDEV_POT_NB_GRID)

            gctx = numgrid_new_atom_grid(radial_precision,       &
                                       min_num_angular_points, &
                                       max_num_angular_points, &
                                       az(cidx),      &
                                       alpha_max(cidx),              &
                                       max_l_quantum_number(cidx),   &
                                       alpha_min(:,cidx))

            npts = numgrid_get_num_grid_points(gctx)

            allocate(grid_x(npts),grid_y(npts),grid_z(npts),grid_w(npts), &
                     stat = alloc_stat)
            if( alloc_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_exchrep_ene_nocache!')
            end if

            call numgrid_get_grid(gctx, nats, cidx-1,   &
                             px, py, pz, az,            &
                             grid_x, grid_y, grid_z, grid_w)

        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

        ! calculate integral
        call ffdev_timers_start_timer(FFDEV_POT_NB_INT)

            rsum = rsum + grid_integrate(mode,px(2),npts,grid_x(:),grid_y(:), &
                                                 grid_z(:),grid_w(:),epa1,epb1,epc1,pd1,epa2,epb2,epc2,pd2)

        call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

        ! destroy grid
        call ffdev_timers_start_timer(FFDEV_POT_NB_GRID)
            deallocate(grid_x, grid_y, grid_z, grid_w)
            call numgrid_free_atom_grid(gctx)
        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

    end do

    select case(mode)
        case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_DENS3A,NB_MODE_PAULI_DENS4A,NB_MODE_PAULI_DENS3B,NB_MODE_PAULI_DENS4B)
            ffdevel_exchrep_ene_numgrid_nocache = epa1*epa2*rsum**pauli_dens_power
        case(NB_MODE_PAULI_WAVE2,NB_MODE_PAULI_WAVE3A,NB_MODE_PAULI_WAVE4A,NB_MODE_PAULI_WAVE3B,NB_MODE_PAULI_WAVE4B)
            ffdevel_exchrep_ene_numgrid_nocache = epa1*epa2*rsum**2/px(2)
        case(NB_MODE_PAULI_LDA2,NB_MODE_PAULI_LDA3A,NB_MODE_PAULI_LDA4A,NB_MODE_PAULI_LDA3B,NB_MODE_PAULI_LDA4B)
            ffdevel_exchrep_ene_numgrid_nocache = rsum
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_nocache!')
    end select

end function ffdevel_exchrep_ene_numgrid_nocache

! ==============================================================================
! function ffdevel_exchrep_ene_numgrid_cache
! ==============================================================================

real(DEVDP) function ffdevel_exchrep_ene_numgrid_cache(cache,mode,pa1,pb1,pc1,pd1,pa2,pb2,pc2,pd2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    type(GRID_CACHE)    :: cache
    integer             :: mode
    real(DEVDP)         :: pa1,pb1,pc1,pd1
    real(DEVDP)         :: pa2,pb2,pc2,pd2
    ! --------------------------------------------
    real(DEVDP)         :: lr
    real(DEVDP)         :: rsum, eexch
    real(DEVDP)         :: epa1,epa2,epb1,epb2,epc1,epc2
    ! --------------------------------------------------------------------------

    ffdevel_exchrep_ene_numgrid_cache = 0.0d0

    if( pb2 .eq. 0.0d0 ) return
    if( pb2 .eq. 0.0d0 ) return

    lr = cache%r * DEV_A2AU

    epa1 = exp(pa1)
    epa2 = exp(pa2)

    select case(mode)
        case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_DENS3B,NB_MODE_PAULI_DENS4B)
            epb1 = 2.0d0*pb1
            epc1 = pc1
            epb2 = 2.0d0*pb2
            epc2 = pc2
        case(NB_MODE_PAULI_DENS3A,NB_MODE_PAULI_DENS4A)
            epb1 = 2.0d0*pb1
            epc1 = 2.0d0*pc1
            epb2 = 2.0d0*pb2
            epc2 = 2.0d0*pc2
        case(NB_MODE_PAULI_WAVE2,NB_MODE_PAULI_WAVE3A,NB_MODE_PAULI_WAVE3B,NB_MODE_PAULI_WAVE4A,NB_MODE_PAULI_WAVE4B)
            epb1 = pb1
            epc1 = pc1
            epb2 = pb2
            epc2 = pc2
        case(NB_MODE_PAULI_LDA2,NB_MODE_PAULI_LDA3B,NB_MODE_PAULI_LDA4B)
            epb1 = 2.0d0*pb1
            epc1 = pc1
            epb2 = 2.0d0*pb2
            epc2 = pc2
        case(NB_MODE_PAULI_LDA3A,NB_MODE_PAULI_LDA4A)
            epb1 = 2.0d0*pb1
            epc1 = 2.0d0*pc1
            epb2 = 2.0d0*pb2
            epc2 = 2.0d0*pc2
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_cache!')
    end select

    rsum = grid_integrate(mode,lr,cache%npts1,cache%grid_1(:,1),cache%grid_1(:,2), &
                    cache%grid_1(:,3),cache%grid_1(:,4),epa1,epb1,epc1,pd1,epa2,epb2,epc2,pd2)
    rsum = rsum + grid_integrate(mode,lr,cache%npts2,cache%grid_2(:,1),cache%grid_2(:,2), &
                    cache%grid_2(:,3),cache%grid_2(:,4),epa1,epb1,epc1,pd1,epa2,epb2,epc2,pd2)

    select case(mode)
        case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_DENS3A,NB_MODE_PAULI_DENS4A,NB_MODE_PAULI_DENS3B,NB_MODE_PAULI_DENS4B)
            ffdevel_exchrep_ene_numgrid_cache = epa1*epa2*rsum**pauli_dens_power
        case(NB_MODE_PAULI_WAVE2,NB_MODE_PAULI_WAVE3A,NB_MODE_PAULI_WAVE4A,NB_MODE_PAULI_WAVE3B,NB_MODE_PAULI_WAVE4B)
            ffdevel_exchrep_ene_numgrid_cache = epa1*epa2*rsum**2/lr
        case(NB_MODE_PAULI_LDA2,NB_MODE_PAULI_LDA3A,NB_MODE_PAULI_LDA4A,NB_MODE_PAULI_LDA3B,NB_MODE_PAULI_LDA4B)
            ffdevel_exchrep_ene_numgrid_cache = rsum
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_cache!')
    end select

end function ffdevel_exchrep_ene_numgrid_cache

! ==============================================================================
! function ffdevel_exchrep_ene_simgrid
! currently only LDA in for testing implemented ...
! ==============================================================================

real(DEVDP) function ffdevel_exchrep_ene_simgrid(mode,r,pa1,pb1,pc1,pd1,pa2,pb2,pc2,pd2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    integer             :: mode
    real(DEVDP)         :: r
    real(DEVDP)         :: pa1,pb1,pc1,pd1
    real(DEVDP)         :: pa2,pb2,pc2,pd2
    ! --------------------------------------------
    real(DEVDP)         :: lr
    real(DEVDP)         :: rsum,ssum
    real(DEVDP)         :: dx,dy,dz,dyz,dv
    integer             :: i, j, k
    real(DEVDP)         :: x,y,z,w1,w2,r1,r2
    ! --------------------------------------------------------------------------

    ffdevel_exchrep_ene_simgrid = 0.0d0

    if( pb2 .eq. 0.0d0 ) return
    if( pb2 .eq. 0.0d0 ) return

    lr = r * DEV_A2AU / 2.0d0

    dx = 0.1 * DEV_A2AU
    dy = 0.1 * DEV_A2AU
    dz = 0.1 * DEV_A2AU

    dv = dx*dy*dz

    rsum = 0.0d0

!$omp parallel

    !$omp do private(i,y,z,j,k,ssum,dyz,r1,r2,w1,w2), reduction(+:rsum)
    do i=1,100
        x = -5.0d0 * DEV_A2AU + i*dx ! this must be derived from i because of omp directive
        y = -3.0 * DEV_A2AU
        ssum = 0.0d0
        do j=1,60
            z = -3.0 * DEV_A2AU
            do k=1,60
                dyz = y**2 + z**2
                r1 = (x-lr)**2 + dyz
                r2 = (x+lr)**2 + dyz
                r1 = sqrt(r1)
                r2 = sqrt(r2)
                w1 = exp(pa1)*exp(-2.0d0*pb1*r1)
                w2 = exp(pa2)*exp(-2.0d0*pb2*r2)

                ssum = ssum + ((w1+w2)**pauli_lda_power - w1**pauli_lda_power - w2**pauli_lda_power)*dv

                z = z + dz
            end do
            y = y + dy
        end do
        rsum = rsum + ssum
    end do

!$omp end parallel

    ffdevel_exchrep_ene_simgrid = rsum

    !write(*,*) r, ffdevel_exchrep_ene_lda, pa1, pa2, pb1, pb2

end function ffdevel_exchrep_ene_simgrid

! ==============================================================================
! function grid_integrate
! ==============================================================================

real(DEVDP) function grid_integrate(mode,lr,npts,gx,gy,gz,gw,epa1,epb1,epc1,pd1,epa2,epb2,epc2,pd2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    integer     :: mode
    real(DEVDP) :: eexch
    real(DEVDP) :: lr
    integer     :: npts
    real(DEVDP) :: gx(:),gy(:),gz(:),gw(:)
    real(DEVDP) :: epa1,epb1,epc1,pd1
    real(DEVDP) :: epa2,epb2,epc2,pd2
    ! --------------------------------------------
    integer     :: ipts
    real(DEVDP) :: rsum
    ! --------------------------------------------------------------------------

    rsum =  0.0d0

    call ffdev_timers_start_timer(FFDEV_POT_NB_INT)

!$omp parallel

    ! calculate integral
        select case(mode)
            case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_WAVE2)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap2(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epb1,epb2)
                end do
            case(NB_MODE_PAULI_DENS3A,NB_MODE_PAULI_WAVE3A)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap3A(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epb1,epc1,epb2,epc2)
                end do
            case(NB_MODE_PAULI_DENS3B,NB_MODE_PAULI_WAVE3B)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap3B(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epb1,epc1,epb2,epc2)
                end do
            case(NB_MODE_PAULI_DENS4A,NB_MODE_PAULI_WAVE4A)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap4A(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epb1,epc1,pd1,epb2,epc2,pd2)
                end do
            case(NB_MODE_PAULI_DENS4B,NB_MODE_PAULI_WAVE4B)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap4B(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epb1,epc1,pd1,epb2,epc2,pd2)
                end do
            case(NB_MODE_PAULI_LDA2)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_lda2(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epa1,epb1,epa2,epb2)
                end do
            case(NB_MODE_PAULI_LDA3A)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_lda3A(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epa1,epb1,epc1,epa2,epb2,epc2)
                end do
            case(NB_MODE_PAULI_LDA3B)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_lda3B(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epa1,epb1,epc1,epa2,epb2,epc2)
                end do
            case(NB_MODE_PAULI_LDA4A)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_lda4A(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epa1,epb1,epc1,pd1,epa2,epb2,epc2,pd2)
                end do
            case(NB_MODE_PAULI_LDA4B)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_lda4B(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,epa1,epb1,epc1,pd1,epa2,epb2,epc2,pd2)
                end do
        case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_cache!')
        end select

!$omp end parallel

    call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

    grid_integrate = rsum

end function grid_integrate

! ------------------------------------------------------------------------------

real(DEVDP)  function get_overlap2(x,y,z,r,pb1,pb2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pb2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap2 = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = exp(-pb1*r1)
    w2 = exp(-pb2*r2)

    get_overlap2 = w1*w2

end function get_overlap2

! ------------------------------------------------------------------------------

real(DEVDP) function get_overlap3A(x,y,z,r,pb1,pc1,pb2,pc2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1
    real(DEVDP)     :: pb2,pc2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap3A = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = (r1**pc1)*exp(-pb1*r1)
    w2 = (r2**pc2)*exp(-pb2*r2)

    get_overlap3A = w1*w2

end function get_overlap3A

! ------------------------------------------------------------------------------

real(DEVDP) function get_overlap3B(x,y,z,r,pb1,pc1,pb2,pc2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1
    real(DEVDP)     :: pb2,pc2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap3B = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = exp(-pb1*r1)*(1.0d0 + pc1/r1)
    w2 = exp(-pb2*r2)*(1.0d0 + pc2/r2)

    get_overlap3B = w1*w2

end function get_overlap3B

! ------------------------------------------------------------------------------

real(DEVDP) function get_overlap4A(x,y,z,r,pb1,pc1,pd1,pb2,pc2,pd2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1,pd1
    real(DEVDP)     :: pb2,pc2,pd2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap4A = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = (r1**pc1)*exp(-pb1*r1)*(1.0d0 + pd1/r1)
    w2 = (r2**pc2)*exp(-pb2*r2)*(1.0d0 + pd2/r2)

    get_overlap4A = w1*w2

end function get_overlap4A

! ------------------------------------------------------------------------------

real(DEVDP) function get_overlap4B(x,y,z,r,pb1,pc1,pd1,pb2,pc2,pd2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1,pd1
    real(DEVDP)     :: pb2,pc2,pd2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap4B = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = exp(-pb1*r1)*(1.0d0 + pc1/r1 + pd1/(r1**2))
    w2 = exp(-pb2*r2)*(1.0d0 + pc2/r2 + pd2/(r2**2))

    get_overlap4B = w1*w2

end function get_overlap4B

! ------------------------------------------------------------------------------

real(DEVDP)  function get_lda2(x,y,z,r,epa1,pb1,epa2,pb2)

    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: epa1,pb1
    real(DEVDP)     :: epa2,pb2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_lda2 = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = epa1*exp(-pb1*r1)
    w2 = epa2*exp(-pb2*r2)

    get_lda2 = (w1+w2)**pauli_lda_power - w1**pauli_lda_power - w2**pauli_lda_power

end function get_lda2

! ------------------------------------------------------------------------------

real(DEVDP) function get_lda3A(x,y,z,r,epa1,pb1,pc1,epa2,pb2,pc2)

    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: epa1,pb1,pc1
    real(DEVDP)     :: epa2,pb2,pc2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_lda3A = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = epa1*(r1**pc1)*exp(-pb1*r1)
    w2 = epa2*(r2**pc2)*exp(-pb2*r2)

    get_lda3A = (w1+w2)**pauli_lda_power - w1**pauli_lda_power - w2**pauli_lda_power

end function get_lda3A

! ------------------------------------------------------------------------------

real(DEVDP) function get_lda3B(x,y,z,r,epa1,pb1,pc1,epa2,pb2,pc2)

    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: epa1,pb1,pc1
    real(DEVDP)     :: epa2,pb2,pc2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_lda3B = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = epa1*exp(-pb1*r1)*(1.0d0 + pc1/r1)
    w2 = epa2*exp(-pb2*r2)*(1.0d0 + pc2/r2)

    get_lda3B = (w1+w2)**pauli_lda_power - w1**pauli_lda_power - w2**pauli_lda_power

end function get_lda3B

! ------------------------------------------------------------------------------

real(DEVDP) function get_lda4A(x,y,z,r,epa1,pb1,pc1,pd1,epa2,pb2,pc2,pd2)

    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: epa1,pb1,pc1,pd1
    real(DEVDP)     :: epa2,pb2,pc2,pd2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_lda4A = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = epa1*(r1**pc1)*exp(-pb1*r1)*(1.0d0 + pd1/r1)
    w2 = epa2*(r2**pc2)*exp(-pb2*r2)*(1.0d0 + pd2/r2)

    get_lda4A = (w1+w2)**pauli_lda_power - w1**pauli_lda_power - w2**pauli_lda_power

end function get_lda4A

! ------------------------------------------------------------------------------

real(DEVDP) function get_lda4B(x,y,z,r,epa1,pb1,pc1,pd1,epa2,pb2,pc2,pd2)

    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: epa1,pb1,pc1,pd1
    real(DEVDP)     :: epa2,pb2,pc2,pd2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_lda4B = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = epa1*exp(-pb1*r1)*(1.0d0 + pc1/r1 + pd1/(r1**2))
    w2 = epa2*exp(-pb2*r2)*(1.0d0 + pc2/r2 + pd2/(r2**2))

    get_lda4B = (w1+w2)**pauli_lda_power - w1**pauli_lda_power - w2**pauli_lda_power

end function get_lda4B

! ------------------------------------------------------------------------------

end module ffdev_exchrep

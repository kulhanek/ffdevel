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

private get_dens_overlap2
private get_dens_overlap3
private get_dens_overlap3p
private get_wave_overlap2
private get_wave_overlap3
private grid_integrate

contains

! ==============================================================================
! subroutine ffdevel_exchrep_gen_cache
! ==============================================================================

subroutine ffdevel_exchrep_gen_cache(cache,z1,z2)

    use, intrinsic :: iso_c_binding, only: c_ptr
    use numgrid
    use ffdev_timers
    use ffdev_utils

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

    ! ne/ne - def2-
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

    az(1) = z1
    az(2) = z2

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

end subroutine ffdevel_exchrep_gen_cache

! ==============================================================================
! function ffdevel_exchrep_ene_nocache
! ==============================================================================

real(DEVDP) function ffdevel_exchrep_ene_nocache(mode,r,z1,pa1,pb1,pc1,pd1,pe1,z2,pa2,pb2,pc2,pd2,pe2)

    use numgrid
    use, intrinsic :: iso_c_binding, only: c_ptr
    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    integer     :: mode
    real(DEVDP) :: r
    integer     :: z1
    real(DEVDP) :: pa1,pb1,pc1,pd1,pe1
    integer     :: z2
    real(DEVDP) :: pa2,pb2,pc2,pd2,pe2
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
    real(DEVDP)                 :: eexch, rsum
    ! --------------------------------------------------------------------------

    ffdevel_exchrep_ene_nocache = 0.0d0

    if( pb2 .eq. 0.0d0 ) return
    if( pb2 .eq. 0.0d0 ) return

    ! numgrid setup
    radial_precision = 1.0d-12
    min_num_angular_points = 86
    max_num_angular_points = 302

    ! ne/ne - def2-
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

    az(1) = z1
    az(2) = z2

    ! coordinates are in a.u.
    px(1) = 0.0d0
    py(1) = 0.0d0
    pz(1) = 0.0d0

    px(2) = r * DEV_A2AU
    py(2) = 0.0d0
    pz(2) = 0.0d0

    eexch = 0.0d0
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
            select case(mode)
                case(NB_MODE_PAULI_DENS2)
                    do ipts = 1, npts
                        rsum = rsum + grid_w(ipts) * get_dens_overlap2(grid_x(ipts), grid_y(ipts), grid_z(ipts), &
                                                       px(2),pb1,pb2)
                    end do
                case(NB_MODE_PAULI_DENS3)
                    do ipts = 1, npts
                        rsum = rsum + grid_w(ipts) * get_dens_overlap3(grid_x(ipts), grid_y(ipts), grid_z(ipts), &
                                                       px(2),pb1,pc1,pb2,pc2)
                    end do
                case(NB_MODE_PAULI_DENS5P)
                    do ipts = 1, npts
                        rsum = rsum + grid_w(ipts) * get_dens_overlap5p(grid_x(ipts), grid_y(ipts), grid_z(ipts), &
                                                       px(2),pb1,pc1,pd1,pe1,pb2,pc2,pd2,pe2)
                    end do
                case(NB_MODE_PAULI_WAVE2)
                    do ipts = 1, npts
                        rsum = rsum + grid_w(ipts) * get_wave_overlap2(grid_x(ipts), grid_y(ipts), grid_z(ipts), &
                                                       px(2),pb1,pb2)
                    end do
                case(NB_MODE_PAULI_WAVE3)
                    do ipts = 1, npts
                        rsum = rsum + grid_w(ipts) * get_wave_overlap3(grid_x(ipts), grid_y(ipts), grid_z(ipts), &
                                                       px(2),pb1,pc1,pb2,pc2)
                    end do
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_nocache!')
            end select
        call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

        ! destroy grid
        call ffdev_timers_start_timer(FFDEV_POT_NB_GRID)
            deallocate(grid_x, grid_y, grid_z, grid_w)
            call numgrid_free_atom_grid(gctx)
        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

    end do

    select case(mode)
        case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_DENS3,NB_MODE_PAULI_DENS5P)
            eexch = rsum*exp(pa1)*exp(pa2)
        case(NB_MODE_PAULI_WAVE2,NB_MODE_PAULI_WAVE3)
            eexch = exp(pa1)*exp(pa2)*rsum**2/px(2)
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_nocache!')
    end select

    ffdevel_exchrep_ene_nocache = eexch

end function ffdevel_exchrep_ene_nocache

! ==============================================================================
! function ffdevel_exchrep_ene_cache
! ==============================================================================

real(DEVDP) function ffdevel_exchrep_ene_cache(cache,mode,pa1,pb1,pc1,pd1,pe1,pa2,pb2,pc2,pd2,pe2)

    use numgrid
    use, intrinsic :: iso_c_binding, only: c_ptr
    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    type(GRID_CACHE)    :: cache
    integer             :: mode
    real(DEVDP)         :: r
    integer             :: z1
    real(DEVDP)         :: pa1,pb1,pc1,pd1,pe1
    integer             :: z2
    real(DEVDP)         :: pa2,pb2,pc2,pd2,pe2
    ! --------------------------------------------
    real(DEVDP)         :: lr
    real(DEVDP)         :: rsum, eexch
    ! --------------------------------------------------------------------------

    ffdevel_exchrep_ene_cache = 0.0d0

    if( pb2 .eq. 0.0d0 ) return
    if( pb2 .eq. 0.0d0 ) return

    lr = cache%r * DEV_A2AU

    rsum = grid_integrate(mode,lr,cache%npts1,cache%grid_1(:,1),cache%grid_1(:,2), &
                    cache%grid_1(:,3),cache%grid_1(:,4),pb1,pc1,pd1,pe1,pb2,pc2,pd2,pe2)
    rsum = rsum + grid_integrate(mode,lr,cache%npts2,cache%grid_2(:,1),cache%grid_2(:,2), &
                    cache%grid_2(:,3),cache%grid_2(:,4),pb1,pc1,pd1,pe1,pb2,pc2,pd2,pe2)

    select case(mode)
        case(NB_MODE_PAULI_DENS2,NB_MODE_PAULI_DENS3,NB_MODE_PAULI_DENS5P)
            eexch = rsum*exp(pa1)*exp(pa2)
        case(NB_MODE_PAULI_WAVE2,NB_MODE_PAULI_WAVE3)
            eexch = exp(pa1)*exp(pa2)*rsum**2/lr
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_cache!')
    end select

    ffdevel_exchrep_ene_cache = eexch

end function ffdevel_exchrep_ene_cache

! ==============================================================================
! function grid_integrate
! ==============================================================================

real(DEVDP) function grid_integrate(mode,lr,npts,gx,gy,gz,gw,pb1,pc1,pd1,pe1,pb2,pc2,pd2,pe2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    integer     :: mode
    real(DEVDP) :: eexch
    real(DEVDP) :: lr
    integer     :: npts
    real(DEVDP) :: gx(:),gy(:),gz(:),gw(:)
    real(DEVDP) :: pb1,pc1,pd1,pe1
    real(DEVDP) :: pb2,pc2,pd2,pe2
    ! --------------------------------------------
    integer     :: ipts
    real(DEVDP) :: rsum
    ! --------------------------------------------------------------------------

    rsum =  0.0d0

    call ffdev_timers_start_timer(FFDEV_POT_NB_INT)

!$omp parallel

    ! calculate integral
        select case(mode)
            case(NB_MODE_PAULI_DENS2)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_dens_overlap2(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,pb1,pb2)
                end do
            case(NB_MODE_PAULI_DENS3)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_dens_overlap3(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,pb1,pc1,pb2,pc2)
                end do
            case(NB_MODE_PAULI_DENS5P)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_dens_overlap5p(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,pb1,pc1,pd1,pe1,pb2,pc2,pd2,pe2)
                end do
            case(NB_MODE_PAULI_WAVE2)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_wave_overlap2(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,pb1,pb2)
                end do
            case(NB_MODE_PAULI_WAVE3)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_wave_overlap3(gx(ipts), gy(ipts), gz(ipts), &
                                                   lr,pb1,pc1,pb2,pc2)
                end do
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_exchrep_ene_cache!')
        end select

!$omp end parallel

    call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

    grid_integrate = rsum

end function grid_integrate

! ------------------------------------------------------------------------------

real(DEVDP)  function get_dens_overlap2(x,y,z,r,pb1,pb2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pb2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_dens_overlap2 = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = exp(-2.0d0*pb1*r1)
    w2 = exp(-2.0d0*pb2*r2)

    get_dens_overlap2 = w1*w2

end function get_dens_overlap2

! ------------------------------------------------------------------------------

real(DEVDP) function get_dens_overlap3(x,y,z,r,pb1,pc1,pb2,pc2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1
    real(DEVDP)     :: pb2,pc2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_dens_overlap3 = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = r1**(2.0d0*pc1)*exp(-2.0d0*pb1*r1)
    w2 = r2**(2.0d0*pc2)*exp(-2.0d0*pb2*r2)

    get_dens_overlap3 = w1*w2

end function get_dens_overlap3

! ------------------------------------------------------------------------------

real(DEVDP) function get_dens_overlap5p(x,y,z,r,pb1,pc1,pd1,pe1,pb2,pc2,pd2,pe2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1,pd1,pe1
    real(DEVDP)     :: pb2,pc2,pd2,pe2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_dens_overlap5p = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = exp(-2.0d0*pb1*r1)*(1.0 + pc1/r1 + pd1/r1**2 + pe1/r1**3)
    w2 = exp(-2.0d0*pb2*r2)*(1.0 + pc2/r2 + pd2/r2**2 + pe2/r2**3)

    get_dens_overlap5p = w1*w2

end function get_dens_overlap5p

! ------------------------------------------------------------------------------

real(DEVDP) function get_wave_overlap2(x,y,z,r,pb1,pb2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pb2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_wave_overlap2 = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = exp(-pb1*r1)
    w2 = exp(-pb2*r2)

    get_wave_overlap2 = w1*w2

end function get_wave_overlap2

! ------------------------------------------------------------------------------

real(DEVDP) function get_wave_overlap3(x,y,z,r,pb1,pc1,pb2,pc2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: r
    real(DEVDP)     :: pb1,pc1
    real(DEVDP)     :: pb2,pc2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_wave_overlap3 = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = x**2 + dyz
    r2 = (x-r)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = r1**(pc1)*exp(-pb1*r1)
    w2 = r2**(pc2)*exp(-pb2*r2)

    get_wave_overlap3 = w1*w2

end function get_wave_overlap3

! ------------------------------------------------------------------------------

end module ffdev_exchrep

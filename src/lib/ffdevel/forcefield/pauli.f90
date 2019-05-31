! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_pauli

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdevel_pauli_gen_numgrid_cache
! ==============================================================================

subroutine ffdevel_pauli_gen_numgrid_cache(cache)

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
    real(DEVDP)         :: px(2),py(2),pz(2),lr,hr
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

    lr = cache%r * DEV_A2AU
    hr = 0.5d0 * lr

    ! coordinates are in a.u.
    px(1) = -hr
    py(1) = 0.0d0
    pz(1) = 0.0d0

    px(2) = hr
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
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_pauli_gen_cache!')
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
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_pauli_gen_cache!')
        end if

        call numgrid_get_grid(gctx, nats, cidx-1,   &
                         px, py, pz, az,            &
                         cache%grid_2(:,1), cache%grid_2(:,2), cache%grid_2(:,3), cache%grid_2(:,4))

        call numgrid_free_atom_grid(gctx)

    call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

end subroutine ffdevel_pauli_gen_numgrid_cache

! ==============================================================================
! function ffdevel_pauli_ene_numgrid_nocache
! ==============================================================================

real(DEVDP) function ffdevel_pauli_ene_numgrid_nocache(mode,r,ta,tb)

    use numgrid
    use, intrinsic :: iso_c_binding, only: c_ptr
    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    integer                     :: mode
    real(DEVDP)                 :: r
    type(NB_TYPE)               :: ta, tb
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
    real(DEVDP)                 :: lr,hr,rsum
    real(DEVDP)                 :: pa1(3),pa2(3)
    real(DEVDP)                 :: pb1(3),pb2(3)
    ! --------------------------------------------------------------------------

    ffdevel_pauli_ene_numgrid_nocache = 0.0d0

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
    lr    = r * DEV_A2AU
    hr    = 0.5d0*lr

    px(1) = -hr
    py(1) = 0.0d0
    pz(1) = 0.0d0

    px(2) = hr
    py(2) = 0.0d0
    pz(2) = 0.0d0

    pa1(1) = exp(ta%pa1)
    pa1(2) = exp(ta%pa2)
    pa1(3) = exp(ta%pa3)

    pa2(1) = exp(tb%pa1)
    pa2(2) = exp(tb%pa2)
    pa2(3) = exp(tb%pa3)

    pb1(1) = ta%pb1
    pb1(2) = ta%pb2
    pb1(3) = ta%pb3

    pb2(1) = tb%pb1
    pb2(2) = tb%pb2
    pb2(3) = tb%pb3

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
                call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_pauli_ene_nocache!')
            end if

            call numgrid_get_grid(gctx, nats, cidx-1,   &
                             px, py, pz, az,            &
                             grid_x, grid_y, grid_z, grid_w)

        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

        ! calculate integral
        call ffdev_timers_start_timer(FFDEV_POT_NB_INT)
            select case(mode)
                case(NB_MODE_PAULI_DENS,NB_MODE_PAULI_WAVE)
                    rsum = rsum + grid_integrate_overlap(lr,npts,grid_x(:),grid_y(:), &
                                    grid_z(:),grid_w(:),pa1,pb1,pa2,pb2)
                case(NB_MODE_PAULI_XFUN)
                    ! FIXME
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_nocache!')
            end select
        call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

        ! destroy grid
        call ffdev_timers_start_timer(FFDEV_POT_NB_GRID)
            deallocate(grid_x, grid_y, grid_z, grid_w)
            call numgrid_free_atom_grid(gctx)
        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

    end do

    select case(mode)
        case(NB_MODE_PAULI_DENS)
            ffdevel_pauli_ene_numgrid_nocache = rsum**pauli_dens_power
        case(NB_MODE_PAULI_WAVE)
            ffdevel_pauli_ene_numgrid_nocache = rsum**2/lr
        case(NB_MODE_PAULI_XFUN)
            ffdevel_pauli_ene_numgrid_nocache = rsum
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_nocache!')
    end select

end function ffdevel_pauli_ene_numgrid_nocache

! ==============================================================================
! function ffdevel_pauli_ene_numgrid_cache
! ==============================================================================

real(DEVDP) function ffdevel_pauli_ene_numgrid_cache(cache,mode,ta,tb)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    type(GRID_CACHE)    :: cache
    integer             :: mode
    type(NB_TYPE)       :: ta, tb
    ! --------------------------------------------
    real(DEVDP)         :: lr
    real(DEVDP)         :: rsum, eexch
    real(DEVDP)         :: pa1(3),pa2(3)
    real(DEVDP)         :: pb1(3),pb2(3)
    ! --------------------------------------------------------------------------

    ffdevel_pauli_ene_numgrid_cache = 0.0d0

    lr = cache%r * DEV_A2AU

    pa1(1) = exp(ta%pa1)
    pa1(2) = exp(ta%pa2)
    pa1(3) = exp(ta%pa3)

    pa2(1) = exp(tb%pa1)
    pa2(2) = exp(tb%pa2)
    pa2(3) = exp(tb%pa3)

    pb1(1) = ta%pb1
    pb1(2) = ta%pb2
    pb1(3) = ta%pb3

    pb2(1) = tb%pb1
    pb2(2) = tb%pb2
    pb2(3) = tb%pb3

    rsum = 0.0d0

    ! write(*,*) pa1,pb1,pa2,pb2

    call ffdev_timers_start_timer(FFDEV_POT_NB_INT)
        select case(mode)
            case(NB_MODE_PAULI_DENS,NB_MODE_PAULI_WAVE)
                rsum = rsum + grid_integrate_overlap(lr,cache%npts1,cache%grid_1(:,1),cache%grid_1(:,2), &
                                cache%grid_1(:,3),cache%grid_1(:,4),pa1,pb1,pa2,pb2)
                rsum = rsum + grid_integrate_overlap(lr,cache%npts2,cache%grid_2(:,1),cache%grid_2(:,2), &
                                cache%grid_2(:,3),cache%grid_2(:,4),pa1,pb1,pa2,pb2)
            case(NB_MODE_PAULI_XFUN)
                        ! FIXME
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_numgrid_cache!')
        end select
    call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

    select case(mode)
        case(NB_MODE_PAULI_DENS)
            ffdevel_pauli_ene_numgrid_cache = rsum**pauli_dens_power
        case(NB_MODE_PAULI_WAVE)
            ffdevel_pauli_ene_numgrid_cache = rsum**2/lr
        case(NB_MODE_PAULI_XFUN)
            ffdevel_pauli_ene_numgrid_cache = rsum
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_numgrid_cache!')
    end select

    ! write(*,*) rsum, ffdevel_pauli_ene_numgrid_cache

end function ffdevel_pauli_ene_numgrid_cache

! ==============================================================================
! function ffdevel_pauli_ene_simgrid
! currently only LDA in for testing implemented ...
! ==============================================================================

real(DEVDP) function ffdevel_pauli_ene_simgrid(mode,r,ta,tb)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    integer             :: mode
    real(DEVDP)         :: r
    type(NB_TYPE)       :: ta, tb
    ! --------------------------------------------
    real(DEVDP)         :: lr,hr
    real(DEVDP)         :: rsum,ssum
    real(DEVDP)         :: dx,dy,dz,dv
    integer             :: i, j, k
    real(DEVDP)         :: x,y,z
    real(DEVDP)         :: pa1(3),pa2(3)
    real(DEVDP)         :: pb1(3),pb2(3)
    ! --------------------------------------------------------------------------

    ffdevel_pauli_ene_simgrid = 0.0d0

    lr = r * DEV_A2AU
    hr = 0.5d0*lr

    pa1(1) = exp(ta%pa1)
    pa1(2) = exp(ta%pa2)
    pa1(3) = exp(ta%pa3)

    pa2(1) = exp(tb%pa1)
    pa2(2) = exp(tb%pa2)
    pa2(3) = exp(tb%pa3)

    pb1(1) = ta%pb1
    pb1(2) = ta%pb2
    pb1(3) = ta%pb3

    pb2(1) = tb%pb1
    pb2(2) = tb%pb2
    pb2(3) = tb%pb3

    dx = 0.1 * DEV_A2AU
    dy = 0.1 * DEV_A2AU
    dz = 0.1 * DEV_A2AU

    dv = dx*dy*dz

    rsum = 0.0d0

!$omp parallel

    !$omp do reduction(+:rsum)
    do i=1,100
        x = -5.0d0 * DEV_A2AU + i*dx ! this must be derived from i because of omp directive
        y = -3.0 * DEV_A2AU
        ssum = 0.0d0
        do j=1,60
            z = -3.0 * DEV_A2AU
            do k=1,60
                ! FIXME - well this would not be too effective due to conditions
                select case(mode)
                    case(NB_MODE_PAULI_DENS,NB_MODE_PAULI_WAVE)
                        select case(pauli_split_valence)
                            case(PAULI_SZ)
                                ssum = ssum + dv * get_overlap_sz(x,y,z,hr,pa1,pb1,pa2,pb2)
                        case default
                                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_simgrid!')
                        end select
                    case(NB_MODE_PAULI_XFUN)
                    case default
                        call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_simgrid!')
                end select

                z = z + dz
            end do
            y = y + dy
        end do
        rsum = rsum + ssum
    end do

!$omp end parallel

    select case(mode)
        case(NB_MODE_PAULI_DENS)
            ffdevel_pauli_ene_simgrid = rsum**pauli_dens_power
        case(NB_MODE_PAULI_WAVE)
            ffdevel_pauli_ene_simgrid = rsum**2/lr
        case(NB_MODE_PAULI_XFUN)
            ffdevel_pauli_ene_simgrid = rsum
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdevel_pauli_ene_simgrid!')
    end select

end function ffdevel_pauli_ene_simgrid

! ==============================================================================
! function grid_integrate_overlap
! ==============================================================================

real(DEVDP) function grid_integrate_overlap(lr,npts,gx,gy,gz,gw,pa1,pb1,pa2,pb2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    real(DEVDP) :: lr
    integer     :: npts
    real(DEVDP) :: gx(:),gy(:),gz(:),gw(:)
    real(DEVDP) :: pa1(3),pa2(3)
    real(DEVDP) :: pb1(3),pb2(3)
    ! --------------------------------------------
    integer     :: ipts
    real(DEVDP) :: rsum, hr
    ! --------------------------------------------------------------------------

    rsum =  0.0d0
    hr = 0.5d0*lr

!$omp parallel

    ! calculate integral
        select case(pauli_split_valence)
            case(PAULI_SZ)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap_sz(gx(ipts),gy(ipts),gz(ipts), &
                                                   hr,pa1,pb1,pa2,pb2)
                end do
            case(PAULI_DZ)
                !$omp do private(ipts), reduction(+:rsum)
                do ipts = 1, npts
                    rsum = rsum + gw(ipts) * get_overlap_dz(gx(ipts),gy(ipts),gz(ipts), &
                                                   hr,pa1,pb1,pa2,pb2)
                end do
        case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in grid_integrate_overlap!')
        end select

!$omp end parallel

    grid_integrate_overlap = rsum

end function grid_integrate_overlap

! ------------------------------------------------------------------------------

real(DEVDP)  function get_overlap_sz(x,y,z,hr,pa1,pb1,pa2,pb2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr
    real(DEVDP)     :: pa1(3),pa2(3)
    real(DEVDP)     :: pb1(3),pb2(3)
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap_sz = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = pa1(1)*exp(-pb1(1)*r1)
    w2 = pa2(1)*exp(-pb2(1)*r2)

    get_overlap_sz = w1*w2

end function get_overlap_sz

! ------------------------------------------------------------------------------

real(DEVDP)  function get_overlap_dz(x,y,z,hr,pa1,pb1,pa2,pb2)

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr
    real(DEVDP)     :: pa1(3),pa2(3)
    real(DEVDP)     :: pb1(3),pb2(3)
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz
    ! --------------------------------------------------------------------------

    get_overlap_dz = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = pa1(1)*exp(-pb1(1)*r1) + r1*pa1(2)*exp(-pb1(2)*r1)
    w2 = pa2(1)*exp(-pb2(1)*r2) + r2*pa2(2)*exp(-pb2(2)*r2)

    get_overlap_dz = w1*w2

end function get_overlap_dz

! ------------------------------------------------------------------------------

end module ffdev_pauli

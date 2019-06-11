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

private preform_parameters
private finalize_rsum
private grid_integrate

private get_overlap_sv
private get_overlap_rsv

private get_rhop_sv
private get_rhop_rsv

private get_kxeg_sv
private get_kxeg_rsv

contains

! ==============================================================================
! subroutine ffdev_pauli_set_default
! ==============================================================================

subroutine ffdev_pauli_set_default

    use ffdev_pauli_dat

    implicit none
    ! --------------------------------------------

    ! numgrid setup
    pauli_use_numgrid    = .true.        ! use numgrid for integration
    pauli_cache_grid     = .true.        ! cache grid for integration

    ! WFN
    pauli_wf_form        = PAULI_WF_EXP
    pauli_wf_nsto        = 1             ! number of STO orbitals
    pauli_wf2rho_power   = 2             ! wf->rho transformation, only for PDENS and XFUN

    ! PDENS
    pauli_dens_dpower    = 1.0           ! can be optimized via pauli_dp
    pauli_dens_rpower    = 0.0           ! can be optimized via pauli_rp

    ! XFUN
    pauli_xfun           = PAULI_XFUN_RHOP
    pauli_xfun_dpower    = 2.0           ! can be optimized via pauli_xd
    pauli_xfun_kpower    = 2.0d0/3.0d0   ! can be optimized via pauli_xk
    pauli_xfun_xpower    = 1.0d0/3.0d0   ! can be optimized via pauli_xx
    pauli_xfun_xfac      = 1.0           ! can be optimized via pauli_xf

end subroutine ffdev_pauli_set_default

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
    max_num_angular_points = 350

    ! Kr/Kr - def2-TZVPP
    alpha_max(1) =  600250.9757500
    alpha_max(2) =  600250.9757500
    max_l_quantum_number(1) = 3
    max_l_quantum_number(2) = 3

    alpha_min(1,1) = 0.19332399511
    alpha_min(2,1) = 0.11569704458
    alpha_min(3,1) = 0.2510000
    alpha_min(4,1) = 0.6280000

    alpha_min(1,2) = 0.19332399511
    alpha_min(2,2) = 0.11569704458
    alpha_min(3,2) = 0.2510000
    alpha_min(4,2) = 0.6280000

    az(1) = 18 ! Kr
    az(2) = 18 ! Kr

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

real(DEVDP) function ffdevel_pauli_ene_numgrid_nocache(mode,r,z1,t1,z2,t2)

    use numgrid
    use, intrinsic :: iso_c_binding, only: c_ptr
    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    integer                     :: mode
    real(DEVDP)                 :: r
    integer                     :: z1, z2
    type(NB_TYPE)               :: t1, t2
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
    integer                     :: n1, n2
    ! --------------------------------------------------------------------------

    ffdevel_pauli_ene_numgrid_nocache = 0.0d0

    ! numgrid setup
    radial_precision = 1.0d-12
    min_num_angular_points = 86
    max_num_angular_points = 350

    ! Kr/Kr - def2-TZVPP
    alpha_max(1) =  600250.9757500
    alpha_max(2) =  600250.9757500
    max_l_quantum_number(1) = 3
    max_l_quantum_number(2) = 3

    alpha_min(1,1) = 0.19332399511
    alpha_min(2,1) = 0.11569704458
    alpha_min(3,1) = 0.2510000
    alpha_min(4,1) = 0.6280000

    alpha_min(1,2) = 0.19332399511
    alpha_min(2,2) = 0.11569704458
    alpha_min(3,2) = 0.2510000
    alpha_min(4,2) = 0.6280000

    az(1) = 18 ! Kr
    az(2) = 18 ! Kr

    ! coordinates are in a.u.
    lr    = r * DEV_A2AU
    hr    = 0.5d0*lr

    px(1) = -hr
    py(1) = 0.0d0
    pz(1) = 0.0d0

    px(2) = hr
    py(2) = 0.0d0
    pz(2) = 0.0d0

! setup parameters
    call preform_parameters(mode,n1,n2)

! integrate
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
                call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_pauli_ene_numgrid_nocache!')
            end if

            call numgrid_get_grid(gctx, nats, cidx-1,   &
                             px, py, pz, az,            &
                             grid_x, grid_y, grid_z, grid_w)

        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

    ! integrate
        call ffdev_timers_start_timer(FFDEV_POT_NB_INT)
            rsum = rsum + grid_integrate(mode,lr,npts,grid_x,grid_y,grid_z,grid_w, &
                                         n1,t1,n2,t2)
        call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

        ! destroy grid
        call ffdev_timers_start_timer(FFDEV_POT_NB_GRID)
            deallocate(grid_x, grid_y, grid_z, grid_w)
            call numgrid_free_atom_grid(gctx)
        call ffdev_timers_stop_timer(FFDEV_POT_NB_GRID)

    end do

! final processing
    ffdevel_pauli_ene_numgrid_nocache = finalize_rsum(mode,rsum,lr)

end function ffdevel_pauli_ene_numgrid_nocache

! ==============================================================================
! function ffdevel_pauli_ene_numgrid_cache
! ==============================================================================

real(DEVDP) function ffdevel_pauli_ene_numgrid_cache(cache,mode,z1,t1,z2,t2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    type(GRID_CACHE)    :: cache
    integer             :: mode
    integer             :: z1, z2
    type(NB_TYPE)       :: t1, t2
    ! --------------------------------------------
    real(DEVDP)         :: lr
    real(DEVDP)         :: rsum
    integer             :: n1,n2
    ! --------------------------------------------------------------------------

    ffdevel_pauli_ene_numgrid_cache = 0.0d0

    lr = cache%r * DEV_A2AU

! setup parameters
    call preform_parameters(mode,n1,n2)

! integrate
    rsum = 0.0d0
    call ffdev_timers_start_timer(FFDEV_POT_NB_INT)
        rsum = rsum + grid_integrate(mode,lr,cache%npts1,cache%grid_1(:,1),cache%grid_1(:,2), &
                                    cache%grid_1(:,3),cache%grid_1(:,4),n1,t1,n2,t2)
        rsum = rsum + grid_integrate(mode,lr,cache%npts2,cache%grid_2(:,1),cache%grid_2(:,2), &
                                     cache%grid_2(:,3),cache%grid_2(:,4),n1,t1,n2,t2)
    call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

! final processing
    ffdevel_pauli_ene_numgrid_cache = finalize_rsum(mode,rsum,lr)

end function ffdevel_pauli_ene_numgrid_cache

! ==============================================================================
! function ffdevel_pauli_ene_simgrid
! ==============================================================================

real(DEVDP) function ffdevel_pauli_ene_simgrid(mode,r,z1,t1,z2,t2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    integer                     :: mode
    real(DEVDP)                 :: r
    integer                     :: z1, z2
    type(NB_TYPE)               :: t1, t2
    ! --------------------------------------------
    real(DEVDP)                 :: lr,hr
    real(DEVDP)                 :: rsum,ssum
    real(DEVDP)                 :: dx,dy,dz,dv
    integer                     :: i, j, k, nx, ny, nz, npts, pts
    real(DEVDP)                 :: x,y,z
    integer                     :: n1,n2,alloc_stat
    real(DEVDP), allocatable    :: grid_x(:)
    real(DEVDP), allocatable    :: grid_y(:)
    real(DEVDP), allocatable    :: grid_z(:)
    real(DEVDP), allocatable    :: grid_w(:)
    ! --------------------------------------------------------------------------

    ffdevel_pauli_ene_simgrid = 0.0d0

    lr = r * DEV_A2AU

! setup parameters
    call preform_parameters(mode,n1,n2)

! grid parameters
    nx = 100
    ny = 60
    nz = 60
    dx = 0.1 * DEV_A2AU
    dy = 0.1 * DEV_A2AU
    dz = 0.1 * DEV_A2AU

! generate simgrid
    dv = dx*dy*dz
    npts = nx*ny*nz

    allocate(grid_x(npts),grid_y(npts),grid_z(npts),grid_w(npts), &
             stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate grid buffers in ffdevel_pauli_ene_simgrid!')
    end if

    pts = 1
    x = -dx*(nx/2)
    do i=1,nx
        y = -dy*(ny/2)
        do j=1,ny
            z = -dz*(nz/2)
            do k=1,nz
                grid_x(pts) = x
                grid_y(pts) = y
                grid_z(pts) = z
                grid_w(pts) = dv
                pts = pts + 1
                z = z + dz
            end do
            y = y + dy
        end do
        x = x + dx
    end do

! integrate
    rsum = 0.0d0
    call ffdev_timers_start_timer(FFDEV_POT_NB_INT)
        rsum = rsum + grid_integrate(mode,lr,npts,grid_x,grid_y,grid_z,grid_w, &
                                     n1,t1,n2,t2)
    call ffdev_timers_stop_timer(FFDEV_POT_NB_INT)

! release grid
    deallocate(grid_x, grid_y, grid_z, grid_w)

! final processing
    ffdevel_pauli_ene_simgrid = finalize_rsum(mode,rsum,lr)

end function ffdevel_pauli_ene_simgrid

! ==============================================================================
! function finalize_rsum
! ==============================================================================

subroutine preform_parameters(mode,n1,n2)

    use ffdev_pauli_dat
    use ffdev_utils
    use ffdev_topology_dat
    use ffdev_topology

    implicit none
    integer         :: mode
    integer         :: n1,n2
    ! --------------------------------------------------------------------------

    used_wf2rho_power = pauli_wf2rho_power
    if( mode .eq. NB_MODE_PAULI_WAVE ) then
        used_wf2rho_power = 1
    end if

    n1 = pauli_wf_nsto
    n2 = pauli_wf_nsto

end subroutine preform_parameters

! ==============================================================================
! function finalize_rsum
! ==============================================================================

real(DEVDP) function finalize_rsum(mode,rsum,lr)

    use ffdev_pauli_dat
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    integer     :: mode
    real(DEVDP) :: rsum
    real(DEVDP) :: lr
    ! --------------------------------------------------------------------------

    select case(mode)
        case(NB_MODE_PAULI_DENS)
            finalize_rsum = (rsum**pauli_dens_dpower)/(lr**pauli_dens_rpower)
        case(NB_MODE_PAULI_WAVE)
            finalize_rsum = rsum**2/lr
        case(NB_MODE_PAULI_XFUN)
            finalize_rsum = rsum
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in finalize_rsum!')
    end select

end function finalize_rsum

! ==============================================================================
! function grid_integrate
! ==============================================================================

real(DEVDP) function grid_integrate(mode,lr,npts,gx,gy,gz,gw,n1,t1,n2,t2)

    use ffdev_timers
    use ffdev_topology_dat
    use ffdev_pauli_dat
    use ffdev_utils

    implicit none
    integer         :: mode
    real(DEVDP)     :: lr
    integer         :: npts
    real(DEVDP)     :: gx(:),gy(:),gz(:),gw(:)
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    integer         :: ipts
    real(DEVDP)     :: rsum, hr
    ! --------------------------------------------------------------------------

    rsum =  0.0d0
    hr = 0.5d0*lr

    !$omp parallel

    select case(mode)
        case(NB_MODE_PAULI_DENS,NB_MODE_PAULI_WAVE)
! ------------------------------------------------------
            select case(pauli_wf_form)
                case(PAULI_WF_EXP)
                    !$omp do private(ipts), reduction(+:rsum)
                    do ipts = 1, npts
                        rsum = rsum + gw(ipts) * get_overlap_sv(gx(ipts),gy(ipts),gz(ipts), &
                                                       hr,n1,t1,n2,t2)
                    end do
                case(PAULI_WF_REXP)
                    !$omp do private(ipts), reduction(+:rsum)
                    do ipts = 1, npts
                        rsum = rsum + gw(ipts) * get_overlap_rsv(gx(ipts),gy(ipts),gz(ipts), &
                                                       hr,n1,t1,n2,t2)
                    end do
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Not implemented in grid_integrate Ia!')
            end select
! ------------------------------------------------------
        case(NB_MODE_PAULI_XFUN)
            select case(pauli_xfun)
            ! ------------------------------------------------------
                case(PAULI_XFUN_RHOP)
                    select case(pauli_wf_form)
                        case(PAULI_WF_EXP)
                            !$omp do private(ipts), reduction(+:rsum)
                            do ipts = 1, npts
                                rsum = rsum + gw(ipts) * get_rhop_sv(gx(ipts),gy(ipts),gz(ipts), &
                                                               hr,n1,t1,n2,t2)
                            end do
                        case(PAULI_WF_REXP)
                            !$omp do private(ipts), reduction(+:rsum)
                            do ipts = 1, npts
                                rsum = rsum + gw(ipts) * get_rhop_rsv(gx(ipts),gy(ipts),gz(ipts), &
                                                               hr,n1,t1,n2,t2)
                            end do
                        case default
                            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in grid_integrate Ib!')
                    end select
            ! ------------------------------------------------------
                case(PAULI_XFUN_KXEG)
                    select case(pauli_wf_form)
                        case(PAULI_WF_EXP)
                            !$omp do private(ipts), reduction(+:rsum)
                            do ipts = 1, npts
                                rsum = rsum + gw(ipts) * get_kxeg_sv(gx(ipts),gy(ipts),gz(ipts), &
                                                               hr,n1,t1,n2,t2)
                            end do
                        case(PAULI_WF_REXP)
                            !$omp do private(ipts), reduction(+:rsum)
                            do ipts = 1, npts
                                rsum = rsum + gw(ipts) * get_kxeg_rsv(gx(ipts),gy(ipts),gz(ipts), &
                                                               hr,n1,t1,n2,t2)
                            end do
                        case default
                            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in grid_integrate Ic!')
                    end select
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Not implemented in grid_integrate III!')
            end select
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in grid_integrate II!')
    end select

    !$omp end parallel

    grid_integrate = rsum

end function grid_integrate

! ------------------------------------------------------------------------------

real(DEVDP)  function get_overlap_sv(x,y,z,hr,n1,t1,n2,t2)

    use ffdev_pauli_dat
    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr           ! half of the distance !!!!
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz, b1, b2
    integer         :: i
    ! --------------------------------------------------------------------------

    get_overlap_sv = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = 0.0d0
    w2 = 0.0d0

    b1 = 0.0d0
    do i=1,n1
        b1 = b1 + t1%pb(i)
        w1 = w1 + t1%pa(i)*exp(-b1*r1)
    end do
    b2 = 0.0d0
    do i=1,n2
        b2 = b2 + t2%pb(i)
        w2 = w2 + t2%pa(i)*exp(-b2*r2)
    end do

    get_overlap_sv = (w1**used_wf2rho_power)*(w2**used_wf2rho_power)

end function get_overlap_sv

! ------------------------------------------------------------------------------

real(DEVDP)  function get_rhop_sv(x,y,z,hr,n1,t1,n2,t2)

    use ffdev_pauli_dat
    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr           ! half of the distance !!!!
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz, d1, d2, b1, b2
    integer         :: i
    ! --------------------------------------------------------------------------

    get_rhop_sv = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = 0.0d0
    w2 = 0.0d0

    b1 = 0.0d0
    do i=1,n1
        b1 = b1 + t1%pb(i)
        w1 = w1 + t1%pa(i)*exp(-b1*r1)
    end do
    b2 = 0.0d0
    do i=1,n2
        b2 = b2 + t2%pb(i)
        w2 = w2 + t2%pa(i)*exp(-b2*r2)
    end do

    d1 = w1**used_wf2rho_power
    d2 = w2**used_wf2rho_power

    get_rhop_sv = 0.5d0*(d1 + d2)**pauli_xfun_dpower - d1**pauli_xfun_dpower - d2**pauli_xfun_dpower

end function get_rhop_sv

! ------------------------------------------------------------------------------

real(DEVDP)  function get_kxeg_sv(x,y,z,hr,n1,t1,n2,t2)

    use ffdev_pauli_dat
    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr           ! half of the distance !!!!
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz, d1, d2, ke, ee, b1, b2
    integer         :: i
    ! --------------------------------------------------------------------------

    get_kxeg_sv = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = 0.0d0
    w2 = 0.0d0

    b1 = 0.0d0
    do i=1,n1
        b1 = b1 + t1%pb(i)
        w1 = w1 + t1%pa(i)*exp(-b1*r1)
    end do
    b2 = 0.0d0
    do i=1,n2
        b2 = b2 + t2%pb(i)
        w2 = w2 + t2%pa(i)*exp(-b2*r2)
    end do

    d1 = w1**used_wf2rho_power
    d2 = w2**used_wf2rho_power

    ke = 0.5d0*(d1 + d2)**pauli_xfun_kpower - d1**pauli_xfun_kpower - d2**pauli_xfun_kpower
    ee = 0.5d0*(d1 + d2)**pauli_xfun_xpower - d1**pauli_xfun_xpower - d2**pauli_xfun_xpower

    get_kxeg_sv = ke - ee*pauli_xfun_xfac

end function get_kxeg_sv

! ------------------------------------------------------------------------------

real(DEVDP)  function get_overlap_rsv(x,y,z,hr,n1,t1,n2,t2)

    use ffdev_pauli_dat
    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr           ! half of the distance !!!!
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz, b1, b2
    integer         :: i
    ! --------------------------------------------------------------------------

    get_overlap_rsv = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = 0.0d0
    w2 = 0.0d0

    b1 = 0.0d0
    do i=1,n1
        b1 = b1 + t1%pb(i)
        w1 = w1 + t1%pa(i)*exp(-b1*r1)*(r1**t1%pc(i))
    end do
    b2 = 0.0d0
    do i=1,n2
        b2 = b2 + t2%pb(i)
        w2 = w2 + t2%pa(i)*exp(-b2*r2)*(r2**t2%pc(i))
    end do

    get_overlap_rsv = (w1**used_wf2rho_power)*(w2**used_wf2rho_power)

end function get_overlap_rsv

! ------------------------------------------------------------------------------

real(DEVDP)  function get_rhop_rsv(x,y,z,hr,n1,t1,n2,t2)

    use ffdev_pauli_dat
    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr           ! half of the distance !!!!
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz, d1, d2, b1, b2
    integer         :: i
    ! --------------------------------------------------------------------------

    get_rhop_rsv = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = 0.0d0
    w2 = 0.0d0

    b1 = 0.0d0
    do i=1,n1
        b1 = b1 + t1%pb(i)
        w1 = w1 + t1%pa(i)*exp(-b1*r1)*(r1**t1%pc(i))
    end do
    b2 = 0.0d0
    do i=1,n2
        b2 = b2 + t2%pb(i)
        w2 = w2 + t2%pa(i)*exp(-b2*r2)*(r2**t2%pc(i))
    end do

    d1 = w1**used_wf2rho_power
    d2 = w2**used_wf2rho_power

    get_rhop_rsv = 0.5d0*(d1 + d2)**pauli_xfun_dpower - d1**pauli_xfun_dpower - d2**pauli_xfun_dpower

end function get_rhop_rsv

! ------------------------------------------------------------------------------

real(DEVDP)  function get_kxeg_rsv(x,y,z,hr,n1,t1,n2,t2)

    use ffdev_pauli_dat
    use ffdev_topology_dat

    implicit none
    real(DEVDP)     :: x, y, z
    real(DEVDP)     :: hr           ! half of the distance !!!!
    integer         :: n1,n2
    type(NB_TYPE)   :: t1, t2
    ! --------------------------------------------
    real(DEVDP)     :: r1, r2, w1, w2, dyz, d1, d2, ke, ee, b1, b2
    integer         :: i
    ! --------------------------------------------------------------------------

    get_kxeg_rsv = 0.0

    ! geometry calculated here MUST follow grid construction (position of atom centers)
    dyz = y**2 + z**2
    r1 = (x+hr)**2 + dyz
    r2 = (x-hr)**2 + dyz

    r1 = sqrt(r1)
    r2 = sqrt(r2)

    w1 = 0.0d0
    w2 = 0.0d0

    b1 = 0.0d0
    do i=1,n1
        b1 = b1 + t1%pb(i)
        w1 = w1 + t1%pa(i)*exp(-b1*r1)*(r1**t1%pc(i))
    end do
    b2 = 0.0d0
    do i=1,n2
        b2 = b2 + t2%pb(i)
        w2 = w2 + t2%pa(i)*exp(-b2*r2)*(r2**t2%pc(i))
    end do

    d1 = w1**used_wf2rho_power
    d2 = w2**used_wf2rho_power

    ke = 0.5d0*(d1 + d2)**pauli_xfun_kpower - d1**pauli_xfun_kpower - d2**pauli_xfun_kpower
    ee = 0.5d0*(d1 + d2)**pauli_xfun_xpower - d1**pauli_xfun_xpower - d2**pauli_xfun_xpower

    get_kxeg_rsv = ke - ee*pauli_xfun_xfac

end function get_kxeg_rsv

! ------------------------------------------------------------------------------

end module ffdev_pauli

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

module ffdev_hessian

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_hessian_all
! ==============================================================================

subroutine ffdev_hessian_all(top,geo,skipnb)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_timers

    implicit none
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    logical,optional    :: skipnb
    ! -------------------------------------------
    logical             :: calcnb
    ! --------------------------------------------------------------------------

    calcnb = .true.
    if( present(skipnb) ) then
        calcnb = .not. skipnb
    end if

    call ffdev_timers_start_timer(FFDEV_POT_HESSIAN_TIMER)

    ! reset energy
    geo%bond_ene = 0.0d0
    geo%angle_ene = 0.0d0
    geo%dih_ene = 0.0d0
    geo%impropr_ene = 0.0d0

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    geo%total_ene = 0.0d0
    geo%rst_energy = 0.0d0

    ! reset gradient
    geo%grd(:,:) = 0.0d0

    ! reset hessian
    geo%hess(:,:,:,:) = 0.0d0

    stop 'incorrect hessian for dihedrals'

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_hessian_bonds(top,geo)
        call ffdev_hessian_angles(top,geo)
        call ffdev_hessian_dihedrals(top,geo)
        call ffdev_hessian_impropers(top,geo)
    end if

    ! FIXME
!    if( calcnb ) then
!        ! non-bonded terms
!        select case(top%nb_mode)
!            case(NB_MODE_LJ)
!                call ffdev_hessian_nb_lj(top,geo,+1.0d0)
!            case default
!                call ffdev_utils_exit(DEV_ERR,1,'Unsupported vdW mode in ffdev_hessian_all!')
!        end select
!    end if

    geo%total_ene = geo%bond_ene + geo%angle_ene + geo%dih_ene + geo%impropr_ene &
                  + geo%ele14_ene + geo%rep14_ene + geo%dis14_ene  &
                  + geo%ele_ene + geo%rep_ene + geo%dis_ene

    call ffdev_timers_stop_timer(FFDEV_POT_HESSIAN_TIMER)

end subroutine ffdev_hessian_all

! ==============================================================================
! subroutine ffdev_hessian_ihess_bonded
! ==============================================================================

subroutine ffdev_hessian_ihess_bonded(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_timers

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    call ffdev_timers_start_timer(FFDEV_POT_HESSIAN_TIMER)

    ! reset energy
    geo%bond_ene = 0.0d0
    geo%angle_ene = 0.0d0
    geo%dih_ene = 0.0d0
    geo%impropr_ene = 0.0d0

    geo%ele14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%dis_ene = 0.0d0

    geo%total_ene = 0.0d0
    geo%rst_energy = 0.0d0

    ! reset gradient
    geo%grd(:,:) = 0.0d0

    ! reset hessian
    geo%hess(:,:,:,:) = 0.0d0

    stop 'incorrect hessian for dihedrals'

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_hessian_bonds(top,geo)
        call ffdev_hessian_angles(top,geo)
        call ffdev_hessian_dihedrals(top,geo)
        call ffdev_hessian_impropers(top,geo)
    end if

    geo%total_ene = geo%bond_ene + geo%angle_ene + geo%dih_ene &
                  + geo%impropr_ene

    call ffdev_timers_stop_timer(FFDEV_POT_HESSIAN_TIMER)

end subroutine ffdev_hessian_ihess_bonded

! ==============================================================================
! subroutine ffdev_hessian_num_all
! ==============================================================================

subroutine ffdev_hessian_num_all(top,geo,skipnb)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use ffdev_gradient
    use ffdev_utils
    use ffdev_timers

    implicit none
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    logical,optional    :: skipnb
    ! -------------------------------------------
    logical             :: lskipnb
    type(GEOMETRY)      :: tmp_geo
    real(DEVDP)         :: d,ene1,ene2,d2inv,eb
    integer             :: i,j,k,l
    ! --------------------------------------------------------------------------

    lskipnb = .false.
    if( present(skipnb) ) then
        lskipnb = skipnb
    end if

    call ffdev_timers_start_timer(FFDEV_POT_HESSIAN_TIMER)

    d = 1.0d-5  ! differentiation parameter
    d2inv = 1.0d0/d**2

    ! calculate base energy
    call ffdev_gradient_all(top,geo,lskipnb)

    ! base energy
    eb = geo%total_ene

    ! allocate temporary geometry object
    call ffdev_geometry_allocate(tmp_geo,geo%natoms)
    tmp_geo%crd(:,:) = geo%crd(:,:)

    ! hessian by numerical differentiation - diagonal items
    do i=1,geo%natoms
        do j=1,3
            ! left
            tmp_geo%crd(j,i) = geo%crd(j,i) + d
            call ffdev_energy_all(top,tmp_geo,lskipnb)
            ene1 = tmp_geo%total_ene

            ! right
            tmp_geo%crd(j,i) = geo%crd(j,i) - d
            call ffdev_energy_all(top,tmp_geo,lskipnb)
            ene2 = tmp_geo%total_ene

            ! hessian
            geo%hess(j,i,j,i) = (ene1-2.0d0*eb+ene2)*d2inv

            ! move back
            tmp_geo%crd(j,i) = geo%crd(j,i)
        end do
    end do

    ! off diagonal elements
    do i=1,geo%natoms
        do j=1,3
            do k=1,geo%natoms
                do l=1,3
                    if( (3*(i-1)+j) .le. (3*(k-1)+l) ) cycle
                    ! left
                    tmp_geo%crd(j,i) = geo%crd(j,i) + d
                    tmp_geo%crd(l,k) = geo%crd(l,k) + d
                    call ffdev_energy_all(top,tmp_geo,lskipnb)
                    ene1 = tmp_geo%total_ene

                    ! right
                    tmp_geo%crd(j,i) = geo%crd(j,i) - d
                    tmp_geo%crd(l,k) = geo%crd(l,k) - d
                    call ffdev_energy_all(top,tmp_geo,lskipnb)
                    ene2 = tmp_geo%total_ene

                    ! hessian
                    geo%hess(j,i,l,k) = 0.5d0*(ene1-2.0d0*eb+ene2)*d2inv &
                                      - 0.5d0*(geo%hess(j,i,j,i)+geo%hess(l,k,l,k))
                    geo%hess(l,k,j,i) = geo%hess(j,i,l,k)

                    ! move back
                    tmp_geo%crd(j,i) = geo%crd(j,i)
                    tmp_geo%crd(l,k) = geo%crd(l,k)
                end do
            end do
        end do
    end do

    ! release temporary geometry object
    deallocate(tmp_geo%crd)

    call ffdev_timers_stop_timer(FFDEV_POT_HESSIAN_TIMER)

end subroutine ffdev_hessian_num_all

! ==============================================================================
! subroutine ffdev_hessian_num_by_grds_all
! ==============================================================================

subroutine ffdev_hessian_num_by_grds_all(top,geo,skipnb)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use ffdev_gradient
    use ffdev_gradient_utils
    use ffdev_utils
    use ffdev_timers

    implicit none
    type(TOPOLOGY)      :: top
    type(GEOMETRY)      :: geo
    logical,optional    :: skipnb
    ! -------------------------------------------
    logical                     :: lskipnb
    type(GEOMETRY),allocatable  :: tmp_geo1(:,:),tmp_geo2(:,:)
    real(DEVDP)                 :: d
    integer                     :: i,j,k,l
    ! --------------------------------------------------------------------------

    lskipnb = .false.
    if( present(skipnb) ) then
        lskipnb = skipnb
    end if

    call ffdev_timers_start_timer(FFDEV_POT_HESSIAN_TIMER)

    d = 1.0d-6  ! differentiation parameter

    ! calculate base energy and gradient
    call ffdev_gradient_all(top,geo,lskipnb)

    ! allocate temporary geometry objects and calculate perturbed gradients
    allocate( tmp_geo1(3,geo%natoms), tmp_geo2(3,geo%natoms))
    do i=1,geo%natoms
        do j=1,3
            call ffdev_geometry_allocate(tmp_geo1(j,i),geo%natoms)
            call ffdev_gradient_allocate(tmp_geo1(j,i))
            tmp_geo1(j,i)%crd(:,:) = geo%crd(:,:)
            tmp_geo1(j,i)%crd(j,i) = tmp_geo1(j,i)%crd(j,i) + d
            call ffdev_gradient_all(top,tmp_geo1(j,i),lskipnb)

            call ffdev_geometry_allocate(tmp_geo2(j,i),geo%natoms)
            call ffdev_gradient_allocate(tmp_geo2(j,i))
            tmp_geo2(j,i)%crd(:,:) = geo%crd(:,:)
            tmp_geo2(j,i)%crd(j,i) = tmp_geo2(j,i)%crd(j,i) - d
            call ffdev_gradient_all(top,tmp_geo2(j,i),lskipnb)
        end do
    end do

    ! hessian by central differences from gradients
    do i=1,geo%natoms
        do j=1,3
            do k=1,geo%natoms
                do l=1,3
                    geo%hess(j,i,l,k) = (tmp_geo1(l,k)%grd(j,i) - tmp_geo2(l,k)%grd(j,i))/(4.0d0*d) &
                                      + (tmp_geo1(j,i)%grd(l,k) - tmp_geo2(j,i)%grd(l,k))/(4.0d0*d)
                end do
            end do
        end do
    end do

    ! release temporary geometry object
    do i=1,geo%natoms
        do j=1,3
            deallocate(tmp_geo1(j,i)%crd,tmp_geo1(j,i)%grd,tmp_geo2(j,i)%crd,tmp_geo2(j,i)%grd)
        end do
    end do
    deallocate(tmp_geo1,tmp_geo2)

    call ffdev_timers_stop_timer(FFDEV_POT_HESSIAN_TIMER)

end subroutine ffdev_hessian_num_by_grds_all

!===============================================================================
! subroutine ffdev_hessian_bonds
!===============================================================================

subroutine ffdev_hessian_bonds(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,ib,ic
    real(DEVDP)     ::  b,db,b2,dn1,dn2,dn3
    real(DEVDP)     ::  h11,h22,h33,h12,h13,h23,h21,h31,h32
    real(DEVDP)     ::  rij(3)
    ! --------------------------------------------------------------------------

    ! reset energy
    geo%bond_ene = 0.0d0

    do ib=1,top%nbonds
        ! for each bond
        i  = top%bonds(ib)%ai
        j  = top%bonds(ib)%aj
        ic = top%bonds(ib)%bt
        ! calculate rij
        rij(:) = geo%crd(:,j) - geo%crd(:,i)

        ! calculate energy
        b2 = rij(1)**2 + rij(2)**2 + rij(3)**2
        b = sqrt ( b2 )
        db = b - top%bond_types(ic)%d0
        geo%bond_ene = geo%bond_ene + 0.5d0*top%bond_types(ic)%k*db**2

        dn1 = top%bond_types(ic)%k*db / b
        dn2 = top%bond_types(ic)%k*db / b**3
        dn3 = top%bond_types(ic)%k / b2

        ! calculate gradient
        geo%grd(:,j) = geo%grd(:,j) + rij(:)*dn1
        geo%grd(:,i) = geo%grd(:,i) - rij(:)*dn1

        ! ----------------------------------------
        ! calculate hessian - diagonals
        h11 = dn1 + (dn3 - dn2)*rij(1)**2
        h22 = dn1 + (dn3 - dn2)*rij(2)**2
        h33 = dn1 + (dn3 - dn2)*rij(3)**2

        geo%hess(1,i,1,i) = geo%hess(1,i,1,i) + h11
        geo%hess(2,i,2,i) = geo%hess(2,i,2,i) + h22
        geo%hess(3,i,3,i) = geo%hess(3,i,3,i) + h33

        geo%hess(1,j,1,j) = geo%hess(1,j,1,j) + h11
        geo%hess(2,j,2,j) = geo%hess(2,j,2,j) + h22
        geo%hess(3,j,3,j) = geo%hess(3,j,3,j) + h33

        ! calculate hessian - off-diagonals - common
        h12 = (dn3 - dn2)*rij(1)*rij(2)
        h13 = (dn3 - dn2)*rij(1)*rij(3)
        h23 = (dn3 - dn2)*rij(2)*rij(3)

        geo%hess(1,i,2,i) = geo%hess(1,i,2,i) + h12
        geo%hess(1,i,3,i) = geo%hess(1,i,3,i) + h13
        geo%hess(2,i,3,i) = geo%hess(2,i,3,i) + h23

        geo%hess(2,i,1,i) = geo%hess(2,i,1,i) + h12
        geo%hess(3,i,1,i) = geo%hess(3,i,1,i) + h13
        geo%hess(3,i,2,i) = geo%hess(3,i,2,i) + h23

        geo%hess(1,j,2,j) = geo%hess(1,j,2,j) + h12
        geo%hess(1,j,3,j) = geo%hess(1,j,3,j) + h13
        geo%hess(2,j,3,j) = geo%hess(2,j,3,j) + h23

        geo%hess(2,j,1,j) = geo%hess(2,j,1,j) + h12
        geo%hess(3,j,1,j) = geo%hess(3,j,1,j) + h13
        geo%hess(3,j,2,j) = geo%hess(3,j,2,j) + h23

       ! calculate hessian - off-diagonals - cross
        h11 = (dn2 - dn3)*rij(1)*rij(1) - dn1
        h12 = (dn2 - dn3)*rij(1)*rij(2)
        h13 = (dn2 - dn3)*rij(1)*rij(3)

        h21 = (dn2 - dn3)*rij(2)*rij(1)
        h22 = (dn2 - dn3)*rij(2)*rij(2) - dn1
        h23 = (dn2 - dn3)*rij(2)*rij(3)

        h31 = (dn2 - dn3)*rij(3)*rij(1)
        h32 = (dn2 - dn3)*rij(3)*rij(2)
        h33 = (dn2 - dn3)*rij(3)*rij(3) - dn1

        geo%hess(1,i,1,j) = geo%hess(1,i,1,j) + h11
        geo%hess(1,i,2,j) = geo%hess(1,i,2,j) + h12
        geo%hess(1,i,3,j) = geo%hess(1,i,3,j) + h13

        geo%hess(2,i,1,j) = geo%hess(2,i,1,j) + h21
        geo%hess(2,i,2,j) = geo%hess(2,i,2,j) + h22
        geo%hess(2,i,3,j) = geo%hess(2,i,3,j) + h23

        geo%hess(3,i,1,j) = geo%hess(3,i,1,j) + h31
        geo%hess(3,i,2,j) = geo%hess(3,i,2,j) + h32
        geo%hess(3,i,3,j) = geo%hess(3,i,3,j) + h33

        geo%hess(1,j,1,i) = geo%hess(1,j,1,i) + h11
        geo%hess(1,j,2,i) = geo%hess(1,j,2,i) + h21
        geo%hess(1,j,3,i) = geo%hess(1,j,3,i) + h31

        geo%hess(2,j,1,i) = geo%hess(2,j,1,i) + h12
        geo%hess(2,j,2,i) = geo%hess(2,j,2,i) + h22
        geo%hess(2,j,3,i) = geo%hess(2,j,3,i) + h32

        geo%hess(3,j,1,i) = geo%hess(3,j,1,i) + h13
        geo%hess(3,j,2,i) = geo%hess(3,j,2,i) + h23
        geo%hess(3,j,3,i) = geo%hess(3,j,3,i) + h33
    end do

end subroutine ffdev_hessian_bonds

!===============================================================================
! subroutine ffdev_hessian_angles
!===============================================================================

subroutine ffdev_hessian_angles(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,j,k,ia,ic
    real(DEVDP)     :: s,scp,angv,da,f0,f1,f2,ui,ti,usi,tsi,ene
    real(DEVDP)     :: rij(3),rkj(3),di(3),dj(3),dk(3),dn1,dn2,dn3,dn4,dn5
    real(DEVDP)     :: lh(3,3),hh(3,3)
    ! --------------------------------------------------------------------------

    ! reset energy
    geo%angle_ene = 0.0d0

    do ia=1,top%nangles
        i  = top%angles(ia)%ai
        j  = top%angles(ia)%aj
        k  = top%angles(ia)%ak
        ic = top%angles(ia)%at

        ! calculate rij and rkj
        rij(:) = geo%crd(:,i) - geo%crd(:,j)
        rkj(:) = geo%crd(:,k) - geo%crd(:,j)

        ! calculate bjiinv and bjkinv and their squares
        ui = 1.0d0 / (rij(1)**2 + rij(2)**2 + rij(3)**2 )
        ti = 1.0d0 / (rkj(1)**2 + rkj(2)**2 + rkj(3)**2 )
        usi = sqrt(ui)
        tsi = sqrt(ti)

        ! calculate scp and angv
        s = rij(1)*rkj(1) + rij(2)*rkj(2) + rij(3)*rkj(3)
        scp = s * usi * tsi
        if ( scp .gt.  1.0d0 ) then
            scp =  1.0d0
        else if ( scp .lt. -1.0d0 ) then
            scp = -1.0d0
        end if
        angv = acos(scp)

        ! calculate energy
        da = angv - top%angle_types(ic)%a0
        ene = 0.5d0*top%angle_types(ic)%k*da**2
        geo%angle_ene = geo%angle_ene + ene

        ! calculate gradient
        f0 = sin ( angv )
        if ( abs(f0) .lt. 1.e-3 ) then
            ! sin(0.1 deg) = 1.7e-3
            ! this is set for angles close to 0 deg or 180 deg by 0.1 deg
            ! the aim is to avoid division be zero
            f0 = 1.0d3
        else
            f0 = 1.0d0 / f0
        end if
        f1 = - top%angle_types(ic)%k*f0*da

        dn1 = usi*tsi
        di(:) =  dn1*rkj(:) - rij(:)*scp*ui
        dk(:) =  dn1*rij(:) - rkj(:)*scp*ti
        dj(:) = - di(:) - dk(:)

        geo%grd(:,i) = geo%grd(:,i) + f1 * di(:)
        geo%grd(:,j) = geo%grd(:,j) + f1 * dj(:)
        geo%grd(:,k) = geo%grd(:,k) + f1 * dk(:)

        ! calculate hessian
        f2 = top%angle_types(ic)%k*(f0**2 - scp*da*f0**3)
        dn1 = usi**3 * tsi
        dn2 = usi**5 * tsi
        dn3 = usi * tsi**3
        dn4 = usi * tsi**5
        dn5 = usi**3 * tsi**3

        ! hessian - diagonals - i
        lh(1,1) = - dn1*(s + 2.0d0*rij(1)*rkj(1)) + 3.0d0*rij(1)**2*s*dn2
        lh(2,2) = - dn1*(s + 2.0d0*rij(2)*rkj(2)) + 3.0d0*rij(2)**2*s*dn2
        lh(3,3) = - dn1*(s + 2.0d0*rij(3)*rkj(3)) + 3.0d0*rij(3)**2*s*dn2
        lh(1,2) = - dn1*(rij(1)*rkj(2) + rij(2)*rkj(1)) + 3.0d0*rij(1)*rij(2)*s*dn2
        lh(1,3) = - dn1*(rij(1)*rkj(3) + rij(3)*rkj(1)) + 3.0d0*rij(1)*rij(3)*s*dn2
        lh(2,3) = - dn1*(rij(2)*rkj(3) + rij(3)*rkj(2)) + 3.0d0*rij(2)*rij(3)*s*dn2

        geo%hess(1,i,1,i) = geo%hess(1,i,1,i) + f1*lh(1,1) + f2*di(1)*di(1)
        geo%hess(2,i,2,i) = geo%hess(2,i,2,i) + f1*lh(2,2) + f2*di(2)*di(2)
        geo%hess(3,i,3,i) = geo%hess(3,i,3,i) + f1*lh(3,3) + f2*di(3)*di(3)
        geo%hess(1,i,2,i) = geo%hess(1,i,2,i) + f1*lh(1,2) + f2*di(1)*di(2)
        geo%hess(1,i,3,i) = geo%hess(1,i,3,i) + f1*lh(1,3) + f2*di(1)*di(3)
        geo%hess(2,i,3,i) = geo%hess(2,i,3,i) + f1*lh(2,3) + f2*di(2)*di(3)
        geo%hess(2,i,1,i) = geo%hess(2,i,1,i) + f1*lh(1,2) + f2*di(1)*di(2)
        geo%hess(3,i,1,i) = geo%hess(3,i,1,i) + f1*lh(1,3) + f2*di(1)*di(3)
        geo%hess(3,i,2,i) = geo%hess(3,i,2,i) + f1*lh(2,3) + f2*di(2)*di(3)

        ! hessian - diagonals - j
        lh(1,1) = 2.0d0*tsi*usi - s*(dn1+dn3) - 2.0d0*(rij(1)*dn1 + rkj(1)*dn3)*(rkj(1) + rij(1)) &
                + 3.0*s*(dn2*rij(1)**2 + dn4*rkj(1)**2) + 2.0d0*rij(1)*rkj(1)*s*dn5
        lh(2,2) = 2.0d0*tsi*usi - s*(dn1+dn3) - 2.0d0*(rij(2)*dn1 + rkj(2)*dn3)*(rkj(2) + rij(2)) &
                + 3.0*s*(dn2*rij(2)**2 + dn4*rkj(2)**2) + 2.0d0*rij(2)*rkj(2)*s*dn5
        lh(3,3) = 2.0d0*tsi*usi - s*(dn1+dn3) - 2.0d0*(rij(3)*dn1 + rkj(3)*dn3)*(rkj(3) + rij(3)) &
                + 3.0*s*(dn2*rij(3)**2 + dn4*rkj(3)**2) + 2.0d0*rij(3)*rkj(3)*s*dn5

        lh(1,2) = - rij(1)*(rkj(2)+rij(2))*dn1 - rij(2)*(rkj(1)+rij(1))*dn1 + 3.0*rij(1)*rij(2)*s*dn2 &
                - (rkj(1)+rij(1))*rkj(2)*dn3 - (rkj(2)+rij(2))*rkj(1)*dn3 + s*dn5*(rij(1)*rkj(2)+rij(2)*rkj(1)) &
                + 3.0*s*dn4*rkj(1)*rkj(2)
        lh(1,3) = - rij(1)*(rkj(3)+rij(3))*dn1 - rij(3)*(rkj(1)+rij(1))*dn1 + 3.0*rij(1)*rij(3)*s*dn2 &
                - (rkj(1)+rij(1))*rkj(3)*dn3 - (rkj(3)+rij(3))*rkj(1)*dn3 + s*dn5*(rij(1)*rkj(3)+rij(3)*rkj(1)) &
                + 3.0*s*dn4*rkj(1)*rkj(3)
        lh(2,3) = - rij(2)*(rkj(3)+rij(3))*dn1 - rij(3)*(rkj(2)+rij(2))*dn1 + 3.0*rij(2)*rij(3)*s*dn2 &
                - (rkj(2)+rij(2))*rkj(3)*dn3 - (rkj(3)+rij(3))*rkj(2)*dn3 + s*dn5*(rij(2)*rkj(3)+rij(3)*rkj(2)) &
                + 3.0*s*dn4*rkj(2)*rkj(3)

        geo%hess(1,j,1,j) = geo%hess(1,j,1,j) + f1*lh(1,1) + f2*dj(1)*dj(1)
        geo%hess(2,j,2,j) = geo%hess(2,j,2,j) + f1*lh(2,2) + f2*dj(2)*dj(2)
        geo%hess(3,j,3,j) = geo%hess(3,j,3,j) + f1*lh(3,3) + f2*dj(3)*dj(3)
        geo%hess(1,j,2,j) = geo%hess(1,j,2,j) + f1*lh(1,2) + f2*dj(1)*dj(2)
        geo%hess(1,j,3,j) = geo%hess(1,j,3,j) + f1*lh(1,3) + f2*dj(1)*dj(3)
        geo%hess(2,j,3,j) = geo%hess(2,j,3,j) + f1*lh(2,3) + f2*dj(2)*dj(3)
        geo%hess(2,j,1,j) = geo%hess(2,j,1,j) + f1*lh(1,2) + f2*dj(1)*dj(2)
        geo%hess(3,j,1,j) = geo%hess(3,j,1,j) + f1*lh(1,3) + f2*dj(1)*dj(3)
        geo%hess(3,j,2,j) = geo%hess(3,j,2,j) + f1*lh(2,3) + f2*dj(2)*dj(3)

        ! hessian - diagonals - k
        lh(1,1) = - dn3*(s + 2.0d0*rij(1)*rkj(1)) + 3.0*rkj(1)**2*s*dn4
        lh(2,2) = - dn3*(s + 2.0d0*rij(2)*rkj(2)) + 3.0*rkj(2)**2*s*dn4
        lh(3,3) = - dn3*(s + 2.0d0*rij(3)*rkj(3)) + 3.0*rkj(3)**2*s*dn4
        lh(1,2) = - dn3*(rij(1)*rkj(2) + rij(2)*rkj(1)) + 3.0*rkj(1)*rkj(2)*s*dn4
        lh(1,3) = - dn3*(rij(1)*rkj(3) + rij(3)*rkj(1)) + 3.0*rkj(1)*rkj(3)*s*dn4
        lh(2,3) = - dn3*(rij(2)*rkj(3) + rij(3)*rkj(2)) + 3.0*rkj(2)*rkj(3)*s*dn4

        geo%hess(1,k,1,k) = geo%hess(1,k,1,k) + f1*lh(1,1) + f2*dk(1)*dk(1)
        geo%hess(2,k,2,k) = geo%hess(2,k,2,k) + f1*lh(2,2) + f2*dk(2)*dk(2)
        geo%hess(3,k,3,k) = geo%hess(3,k,3,k) + f1*lh(3,3) + f2*dk(3)*dk(3)
        geo%hess(1,k,2,k) = geo%hess(1,k,2,k) + f1*lh(1,2) + f2*dk(1)*dk(2)
        geo%hess(1,k,3,k) = geo%hess(1,k,3,k) + f1*lh(1,3) + f2*dk(1)*dk(3)
        geo%hess(2,k,3,k) = geo%hess(2,k,3,k) + f1*lh(2,3) + f2*dk(2)*dk(3)
        geo%hess(2,k,1,k) = geo%hess(2,k,1,k) + f1*lh(1,2) + f2*dk(1)*dk(2)
        geo%hess(3,k,1,k) = geo%hess(3,k,1,k) + f1*lh(1,3) + f2*dk(1)*dk(3)
        geo%hess(3,k,2,k) = geo%hess(3,k,2,k) + f1*lh(2,3) + f2*dk(2)*dk(3)

        ! hessian - off-diagonal - i,k
        lh(1,1) = usi*tsi - dn1*rij(1)**2 - dn3*rkj(1)**2 + rij(1)*rkj(1)*s*dn5
        lh(2,2) = usi*tsi - dn1*rij(2)**2 - dn3*rkj(2)**2 + rij(2)*rkj(2)*s*dn5
        lh(3,3) = usi*tsi - dn1*rij(3)**2 - dn3*rkj(3)**2 + rij(3)*rkj(3)*s*dn5
        lh(1,2) = - rij(1)*rij(2)*dn1 - rkj(1)*rkj(2)*dn3 + rij(1)*rkj(2)*s*dn5
        lh(1,3) = - rij(1)*rij(3)*dn1 - rkj(1)*rkj(3)*dn3 + rij(1)*rkj(3)*s*dn5
        lh(2,3) = - rij(2)*rij(3)*dn1 - rkj(2)*rkj(3)*dn3 + rij(2)*rkj(3)*s*dn5
        lh(2,1) = - rij(2)*rij(1)*dn1 - rkj(2)*rkj(1)*dn3 + rij(2)*rkj(1)*s*dn5
        lh(3,1) = - rij(3)*rij(1)*dn1 - rkj(3)*rkj(1)*dn3 + rij(3)*rkj(1)*s*dn5
        lh(3,2) = - rij(3)*rij(2)*dn1 - rkj(3)*rkj(2)*dn3 + rij(3)*rkj(2)*s*dn5

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dk(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dk(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dk(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dk(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dk(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dk(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dk(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dk(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dk(2)

        geo%hess(1,i,1,k) = geo%hess(1,i,1,k) + hh(1,1)
        geo%hess(2,i,2,k) = geo%hess(2,i,2,k) + hh(2,2)
        geo%hess(3,i,3,k) = geo%hess(3,i,3,k) + hh(3,3)
        geo%hess(1,i,2,k) = geo%hess(1,i,2,k) + hh(1,2)
        geo%hess(1,i,3,k) = geo%hess(1,i,3,k) + hh(1,3)
        geo%hess(2,i,3,k) = geo%hess(2,i,3,k) + hh(2,3)
        geo%hess(2,i,1,k) = geo%hess(2,i,1,k) + hh(2,1)
        geo%hess(3,i,1,k) = geo%hess(3,i,1,k) + hh(3,1)
        geo%hess(3,i,2,k) = geo%hess(3,i,2,k) + hh(3,2)

        geo%hess(1,k,1,i) = geo%hess(1,k,1,i) + hh(1,1)
        geo%hess(2,k,2,i) = geo%hess(2,k,2,i) + hh(2,2)
        geo%hess(3,k,3,i) = geo%hess(3,k,3,i) + hh(3,3)
        geo%hess(1,k,2,i) = geo%hess(1,k,2,i) + hh(2,1)
        geo%hess(1,k,3,i) = geo%hess(1,k,3,i) + hh(3,1)
        geo%hess(2,k,3,i) = geo%hess(2,k,3,i) + hh(3,2)
        geo%hess(2,k,1,i) = geo%hess(2,k,1,i) + hh(1,2)
        geo%hess(3,k,1,i) = geo%hess(3,k,1,i) + hh(1,3)
        geo%hess(3,k,2,i) = geo%hess(3,k,2,i) + hh(2,3)

        ! hessian - off-diagonal - i,j
        lh(1,1) = - usi*tsi + s*dn1 + rij(1)*rkj(1)*dn1 + rij(1)*(rkj(1)+rij(1))*dn1 &
                - 3.0*s*dn2*rij(1)**2 + dn3*rkj(1)**2 - rij(1)*rkj(1)*s*dn5
        lh(2,2) =- usi*tsi + s*dn1 + rij(2)*rkj(2)*dn1 + rij(2)*(rkj(2)+rij(2))*dn1 &
                - 3.0*s*dn2*rij(2)**2 + dn3*rkj(2)**2 - rij(2)*rkj(2)*s*dn5
        lh(3,3) = - usi*tsi + s*dn1 + rij(3)*rkj(3)*dn1 + rij(3)*(rkj(3)+rij(3))*dn1 &
                - 3.0*s*dn2*rij(3)**2 + dn3*rkj(3)**2 - rij(3)*rkj(3)*s*dn5

        lh(1,2) = rij(1)*(rkj(2)+rij(2))*dn1 + rkj(1)*rij(2)*dn1 - 3.0*rij(1)*rij(2)*s*dn2 &
                + rkj(1)*rkj(2)*dn3 - rij(1)*rkj(2)*s*dn5
        lh(1,3) = rij(1)*(rkj(3)+rij(3))*dn1 + rkj(1)*rij(3)*dn1 - 3.0*rij(1)*rij(3)*s*dn2 &
                + rkj(1)*rkj(3)*dn3 - rij(1)*rkj(3)*s*dn5
        lh(2,3) = rij(2)*(rkj(3)+rij(3))*dn1 + rkj(2)*rij(3)*dn1 - 3.0*rij(2)*rij(3)*s*dn2 &
                + rkj(2)*rkj(3)*dn3 - rij(2)*rkj(3)*s*dn5

        lh(2,1) = rkj(2)*rij(1)*dn1 + rij(2)*(rkj(1)+rij(1))*dn1 - 3.0*rij(2)*rij(1)*s*dn2 &
                + rkj(2)*rkj(1)*dn3 - rij(2)*rkj(1)*s*dn5
        lh(3,1) = rkj(3)*rij(1)*dn1 + rij(3)*(rkj(1)+rij(1))*dn1 - 3.0*rij(3)*rij(1)*s*dn2 &
                + rkj(3)*rkj(1)*dn3 - rij(3)*rkj(1)*s*dn5
        lh(3,2) = rkj(3)*rij(2)*dn1 + rij(3)*(rkj(2)+rij(2))*dn1 - 3.0*rij(3)*rij(2)*s*dn2 &
                + rkj(3)*rkj(2)*dn3 - rij(3)*rkj(2)*s*dn5

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dj(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dj(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dj(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dj(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dj(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dj(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dj(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dj(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dj(2)

        geo%hess(1,i,1,j) = geo%hess(1,i,1,j) + hh(1,1)
        geo%hess(2,i,2,j) = geo%hess(2,i,2,j) + hh(2,2)
        geo%hess(3,i,3,j) = geo%hess(3,i,3,j) + hh(3,3)
        geo%hess(1,i,2,j) = geo%hess(1,i,2,j) + hh(1,2)
        geo%hess(1,i,3,j) = geo%hess(1,i,3,j) + hh(1,3)
        geo%hess(2,i,3,j) = geo%hess(2,i,3,j) + hh(2,3)
        geo%hess(2,i,1,j) = geo%hess(2,i,1,j) + hh(2,1)
        geo%hess(3,i,1,j) = geo%hess(3,i,1,j) + hh(3,1)
        geo%hess(3,i,2,j) = geo%hess(3,i,2,j) + hh(3,2)

        geo%hess(1,j,1,i) = geo%hess(1,j,1,i) + hh(1,1)
        geo%hess(2,j,2,i) = geo%hess(2,j,2,i) + hh(2,2)
        geo%hess(3,j,3,i) = geo%hess(3,j,3,i) + hh(3,3)
        geo%hess(1,j,2,i) = geo%hess(1,j,2,i) + hh(2,1)
        geo%hess(1,j,3,i) = geo%hess(1,j,3,i) + hh(3,1)
        geo%hess(2,j,3,i) = geo%hess(2,j,3,i) + hh(3,2)
        geo%hess(2,j,1,i) = geo%hess(2,j,1,i) + hh(1,2)
        geo%hess(3,j,1,i) = geo%hess(3,j,1,i) + hh(1,3)
        geo%hess(3,j,2,i) = geo%hess(3,j,2,i) + hh(2,3)

        ! hessian - off-diagonal - j,k
        lh(1,1) = - usi*tsi + dn1*rij(1)**2 + s*dn3 + (rkj(1)+rij(1))*rkj(1)*dn3 &
                + rij(1)*rkj(1)*dn3 - rij(1)*rkj(1)*s*dn5 - 3.0*s*dn4*rkj(1)**2

        lh(2,2) = - usi*tsi + dn1*rij(2)**2 + s*dn3 + (rkj(2)+rij(2))*rkj(2)*dn3 &
                + rij(2)*rkj(2)*dn3 - rij(2)*rkj(2)*s*dn5 - 3.0*s*dn4*rkj(2)**2

        lh(3,3) = - usi*tsi + dn1*rij(3)**2 + s*dn3 + (rkj(3)+rij(3))*rkj(3)*dn3 &
                + rij(3)*rkj(3)*dn3 - rij(3)*rkj(3)*s*dn5 - 3.0*s*dn4*rkj(3)**2

        lh(1,2) = rij(1)*rij(2)*dn1 + (rkj(1)+rij(1))*rkj(2)*dn3 + rkj(1)*rij(2)*dn3 &
                - rij(1)*rkj(2)*s*dn5 - 3.0*rkj(1)*rkj(2)*s*dn4
        lh(1,3) = rij(1)*rij(3)*dn1 + (rkj(1)+rij(1))*rkj(3)*dn3 + rkj(1)*rij(3)*dn3 &
                - rij(1)*rkj(3)*s*dn5 - 3.0*rkj(1)*rkj(3)*s*dn4
        lh(2,3) = rij(2)*rij(3)*dn1 + (rkj(2)+rij(2))*rkj(3)*dn3 + rkj(2)*rij(3)*dn3 &
                - rij(2)*rkj(3)*s*dn5 - 3.0*rkj(2)*rkj(3)*s*dn4

        lh(2,1) = rij(2)*rij(1)*dn1 + rkj(2)*rij(1)*dn3 + (rkj(2)+rij(2))*rkj(1)*dn3 &
                - rij(2)*rkj(1)*s*dn5 - 3.0*rkj(2)*rkj(1)*s*dn4
        lh(3,1) = rij(3)*rij(1)*dn1 + rkj(3)*rij(1)*dn3 + (rkj(3)+rij(3))*rkj(1)*dn3 &
                - rij(3)*rkj(1)*s*dn5 - 3.0*rkj(3)*rkj(1)*s*dn4
        lh(3,2) = rij(3)*rij(2)*dn1 + rkj(3)*rij(2)*dn3 + (rkj(3)+rij(3))*rkj(2)*dn3 &
                - rij(3)*rkj(2)*s*dn5 - 3.0*rkj(3)*rkj(2)*s*dn4

        hh(1,1) = f1*lh(1,1) + f2*dj(1)*dk(1)
        hh(2,2) = f1*lh(2,2) + f2*dj(2)*dk(2)
        hh(3,3) = f1*lh(3,3) + f2*dj(3)*dk(3)

        hh(1,2) = f1*lh(1,2) + f2*dj(1)*dk(2)
        hh(1,3) = f1*lh(1,3) + f2*dj(1)*dk(3)
        hh(2,3) = f1*lh(2,3) + f2*dj(2)*dk(3)

        hh(2,1) = f1*lh(2,1) + f2*dj(2)*dk(1)
        hh(3,1) = f1*lh(3,1) + f2*dj(3)*dk(1)
        hh(3,2) = f1*lh(3,2) + f2*dj(3)*dk(2)

        geo%hess(1,j,1,k) = geo%hess(1,j,1,k) + hh(1,1)
        geo%hess(2,j,2,k) = geo%hess(2,j,2,k) + hh(2,2)
        geo%hess(3,j,3,k) = geo%hess(3,j,3,k) + hh(3,3)
        geo%hess(1,j,2,k) = geo%hess(1,j,2,k) + hh(1,2)
        geo%hess(1,j,3,k) = geo%hess(1,j,3,k) + hh(1,3)
        geo%hess(2,j,3,k) = geo%hess(2,j,3,k) + hh(2,3)
        geo%hess(2,j,1,k) = geo%hess(2,j,1,k) + hh(2,1)
        geo%hess(3,j,1,k) = geo%hess(3,j,1,k) + hh(3,1)
        geo%hess(3,j,2,k) = geo%hess(3,j,2,k) + hh(3,2)

        geo%hess(1,k,1,j) = geo%hess(1,k,1,j) + hh(1,1)
        geo%hess(2,k,2,j) = geo%hess(2,k,2,j) + hh(2,2)
        geo%hess(3,k,3,j) = geo%hess(3,k,3,j) + hh(3,3)
        geo%hess(1,k,2,j) = geo%hess(1,k,2,j) + hh(2,1)
        geo%hess(1,k,3,j) = geo%hess(1,k,3,j) + hh(3,1)
        geo%hess(2,k,3,j) = geo%hess(2,k,3,j) + hh(3,2)
        geo%hess(2,k,1,j) = geo%hess(2,k,1,j) + hh(1,2)
        geo%hess(3,k,1,j) = geo%hess(3,k,1,j) + hh(1,3)
        geo%hess(3,k,2,j) = geo%hess(3,k,2,j) + hh(2,3)

    end do

end subroutine ffdev_hessian_angles

!===============================================================================
! subroutine ffdev_hessian_dihedrals
!===============================================================================

subroutine ffdev_hessian_dihedrals(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip,pn
    real(DEVDP)     ::  scp1,scp2,scp3,r2ij,r2kj,r2lk,sarg,carg,dn1,dn2,dn3,dn4,dn5,dn6
    real(DEVDP)     ::  tp,y,phi,f0,f1,f2,arg,th1,th2,ibo1,ibo2
    real(DEVDP)     ::  rij(3),rkj(3),rlk(3),rnj(3),rnk(3)
    real(DEVDP)     ::  di(3),dj(3),dk(3),dl(3),lh(3,3),hh(3,3)
    ! -----------------------------------------------------------------------------

    geo%dih_ene = 0.0d0

    do ip=1,top%ndihedrals

        i  = top%dihedrals(ip)%ai
        j  = top%dihedrals(ip)%aj
        k  = top%dihedrals(ip)%ak
        l  = top%dihedrals(ip)%al

        ic = top%dihedrals(ip)%dt

        rij(:) = geo%crd(:,i) - geo%crd(:,j)
        rkj(:) = geo%crd(:,k) - geo%crd(:,j)
        rlk(:) = geo%crd(:,l) - geo%crd(:,k)

        rnj(1) =  rij(2)*rkj(3) - rij(3)*rkj(2)
        rnj(2) =  rij(3)*rkj(1) - rij(1)*rkj(3)
        rnj(3) =  rij(1)*rkj(2) - rij(2)*rkj(1)

        rnk(1) = -rkj(2)*rlk(3) + rkj(3)*rlk(2)
        rnk(2) = -rkj(3)*rlk(1) + rkj(1)*rlk(3)
        rnk(3) = -rkj(1)*rlk(2) + rkj(2)*rlk(1)

        scp1 = rij(1)*rkj(1) + rij(2)*rkj(2) + rij(3)*rkj(3)
        scp2 = rkj(1)*rlk(1) + rkj(2)*rlk(2) + rkj(3)*rlk(3)
        scp3 = rij(1)*rlk(1) + rij(2)*rlk(2) + rij(3)*rlk(3)

        r2ij = rij(1)**2 + rij(2)**2 + rij(3)**2
        r2kj = rkj(1)**2 + rkj(2)**2 + rkj(3)**2
        r2lk = rlk(1)**2 + rlk(2)**2 + rlk(3)**2

        tp  = scp3*r2kj - scp1*scp2
        ibo1 = 1.0d0 / sqrt(r2ij*r2kj - scp1**2)
        ibo2 = 1.0d0 / sqrt(r2lk*r2kj - scp2**2)

        ! calculate y and phi
        y = tp * ibo1 * ibo2

        if ( y .gt.  1.0 ) then
                y =  1.0
                phi = acos (1.0) ! const
        else if ( y .lt. -1.0 ) then
                y = -1.0
                phi = acos (-1.0) ! const
        else
            phi = acos ( y )
        end if

        ! dihedral sign
        if( rkj(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
           +rkj(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
           +rkj(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) .lt. 0) then
                    phi = -phi
        end if

        f0 = sin ( phi )
        if ( abs(f0) .lt. 1.e-12 ) f0 = 1.e-12
        f0 =  1.0d0 / f0

        f1 = 0.0d0
        f2 = 0.0d0
        do pn=1,top%dihedral_types(ic)%n
            if( .not. top%dihedral_types(ic)%enabled(pn) ) cycle

            ! calculate energy
            arg = pn*phi - top%dihedral_types(ic)%g(pn)
            sarg = sin(arg)
            carg = cos(arg)

            if( dih_cos_only ) then
                geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%v(pn)*carg
            else
                geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%v(pn)*(1.0d0+carg)
            end if

            ! outer derivatives
            f1 = f1 + top%dihedral_types(ic)%v(pn)*pn*sarg*f0
            f2 = f2 + pn    * top%dihedral_types(ic)%v(pn) * y * sarg * f0**3 &
                    - pn**2 * top%dihedral_types(ic)%v(pn)     * carg * f0**2
        end do

        ! helpers
        dn1 = ibo1    * ibo2
        dn2 = ibo1**3 * ibo2
        dn3 = ibo1    * ibo2**3
        dn4 = ibo1**5 * ibo2
        dn5 = ibo1**3 * ibo2**3
        dn6 = ibo1    * ibo2**5

        ! gradients of y
        th1 = (r2kj*scp3-scp1*scp2)*dn2
        th2 = (r2kj*scp3-scp1*scp2)*dn3

        di(1) = dn1*(rlk(1)*r2kj-rkj(1)*scp2) - (rij(1)*r2kj-rkj(1)*scp1)*th1
        di(2) = dn1*(rlk(2)*r2kj-rkj(2)*scp2) - (rij(2)*r2kj-rkj(2)*scp1)*th1
        di(3) = dn1*(rlk(3)*r2kj-rkj(3)*scp2) - (rij(3)*r2kj-rkj(3)*scp1)*th1

        dj(1) = - ( (rkj(1) + rij(1))*scp1 - rij(1)*r2kj - rkj(1)*r2ij )*th1 &
                - ( rlk(1)*scp2 - rkj(1)*r2lk )*th2 &
                + ( (rkj(1)+rij(1))*scp2 + rlk(1)*scp1 - rlk(1)*r2kj - 2.0d0*rkj(1)*scp3)*dn1
        dj(2) = - ( (rkj(2) + rij(2))*scp1 - rij(2)*r2kj - rkj(2)*r2ij )*th1 &
                - ( rlk(2)*scp2 - rkj(2)*r2lk )*th2 &
                + ( (rkj(2)+rij(2))*scp2 + rlk(2)*scp1 - rlk(2)*r2kj - 2.0d0*rkj(2)*scp3)*dn1
        dj(3) = - ( (rkj(3) + rij(3))*scp1 - rij(3)*r2kj - rkj(3)*r2ij )*th1 &
                - ( rlk(3)*scp2 - rkj(3)*r2lk )*th2 &
                + ( (rkj(3)+rij(3))*scp2 + rlk(3)*scp1 - rlk(3)*r2kj - 2.0d0*rkj(3)*scp3)*dn1

        dk(1) = - ( rkj(1)*r2ij - rij(1)*scp1 )*th1 &
                - ( -(rlk(1) - rkj(1))*scp2 + rkj(1)*r2lk - rlk(1)*r2kj )*th2 &
                + ( 2.0d0*rkj(1)*scp3 - rij(1)*scp2 - (rlk(1) - rkj(1))*scp1 - rij(1)*r2kj )*dn1
        dk(2) = - ( rkj(2)*r2ij - rij(2)*scp1 )*th1 &
                - ( -(rlk(2) - rkj(2))*scp2 + rkj(2)*r2lk - rlk(2)*r2kj )*th2 &
                + ( 2.0d0*rkj(2)*scp3 - rij(2)*scp2 - (rlk(2) - rkj(2))*scp1 - rij(2)*r2kj )*dn1
        dk(3) = - ( rkj(3)*r2ij - rij(3)*scp1 )*th1 &
                - ( -(rlk(3) - rkj(3))*scp2 + rkj(3)*r2lk - rlk(3)*r2kj )*th2 &
                + ( 2.0d0*rkj(3)*scp3 - rij(3)*scp2 - (rlk(3) - rkj(3))*scp1 - rij(3)*r2kj )*dn1

        dl(1) = dn1*(rij(1)*r2kj-rkj(1)*scp1) - (rlk(1)*r2kj-rkj(1)*scp2)*th2
        dl(2) = dn1*(rij(2)*r2kj-rkj(2)*scp1) - (rlk(2)*r2kj-rkj(2)*scp2)*th2
        dl(3) = dn1*(rij(3)*r2kj-rkj(3)*scp1) - (rlk(3)*r2kj-rkj(3)*scp2)*th2

        ! calculate gradient
        geo%grd(:,i) = geo%grd(:,i) + f1*di(:)
        geo%grd(:,j) = geo%grd(:,j) + f1*dj(:)
        geo%grd(:,k) = geo%grd(:,k) + f1*dk(:)
        geo%grd(:,l) = geo%grd(:,l) + f1*dl(:)

        ! calculate hessian - diagonal i
        lh(1,1) = - (r2kj-rkj(1)**2)*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * (rij(1)*r2kj-rkj(1)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4 &
                  - 2.0 * (rij(1)*r2kj-rkj(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2
        lh(2,2) = - (r2kj-rkj(2)**2)*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * (rij(2)*r2kj-rkj(2)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4 &
                  - 2.0 * (rij(2)*r2kj-rkj(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2
        lh(3,3) = - (r2kj-rkj(3)**2)*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * (rij(3)*r2kj-rkj(3)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4 &
                  - 2.0 * (rij(3)*r2kj-rkj(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2
        lh(1,2) = + rkj(1)*rkj(2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rij(1)*r2kj-rkj(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2
        lh(1,3) = + rkj(1)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rij(1)*r2kj-rkj(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2
        lh(2,3) = + rkj(2)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rij(2)*r2kj-rkj(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2

        geo%hess(1,i,1,i) = geo%hess(1,i,1,i) + f1*lh(1,1) + f2*di(1)*di(1)
        geo%hess(2,i,2,i) = geo%hess(2,i,2,i) + f1*lh(2,2) + f2*di(2)*di(2)
        geo%hess(3,i,3,i) = geo%hess(3,i,3,i) + f1*lh(3,3) + f2*di(3)*di(3)
        geo%hess(1,i,2,i) = geo%hess(1,i,2,i) + f1*lh(1,2) + f2*di(1)*di(2)
        geo%hess(1,i,3,i) = geo%hess(1,i,3,i) + f1*lh(1,3) + f2*di(1)*di(3)
        geo%hess(2,i,3,i) = geo%hess(2,i,3,i) + f1*lh(2,3) + f2*di(2)*di(3)
        geo%hess(2,i,1,i) = geo%hess(2,i,1,i) + f1*lh(1,2) + f2*di(1)*di(2)
        geo%hess(3,i,1,i) = geo%hess(3,i,1,i) + f1*lh(1,3) + f2*di(1)*di(3)
        geo%hess(3,i,2,i) = geo%hess(3,i,2,i) + f1*lh(2,3) + f2*di(2)*di(3)

        ! calculate hessian - diagonal j
        lh(1,1) = - (-2.0d0*scp1 + r2kj + r2ij + 4.0*rij(1)*rkj(1) - (rkj(1)+rij(1))**2 )*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * ( (rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)**2 * (r2kj*scp3-scp1*scp2)*dn4 &
                  + 2.0 * ( (rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)    * (rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2)*dn5 &
                  - (r2lk-rlk(1)**2)*(r2kj*scp3-scp1*scp2)*dn3 &
                  + 3.0 * (rlk(1)*scp2 - rkj(1)*r2lk)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * ( (rkj(1)+rij(1))*scp1 - rij(1)*r2kj - rkj(1)*r2ij) &
                        * (-2.0d0*rkj(1)*scp3 + (rkj(1)+rij(1))*scp2 + rlk(1)*scp1 - rlk(1)*r2kj) * dn2 &
                  - 2.0 * ( rlk(1)*scp2-rkj(1)*r2lk)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj)*dn3 &
                  + (2.0 * scp3 - 2.0d0*scp2 + 4.0*rkj(1)*rlk(1) - 2.0d0*(rkj(1)+rij(1))*rlk(1)) * dn1
        lh(2,2) =  - (-2.0d0*scp1 + r2kj + r2ij + 4.0*rij(2)*rkj(2) - (rkj(2)+rij(2))**2 )*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * ( (rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)**2 * (r2kj*scp3-scp1*scp2)*dn4 &
                  + 2.0 * ( (rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)    * (rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2)*dn5 &
                  - (r2lk-rlk(2)**2)*(r2kj*scp3-scp1*scp2)*dn3 &
                  + 3.0 * (rlk(2)*scp2 - rkj(2)*r2lk)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * ( (rkj(2)+rij(2))*scp1 - rij(2)*r2kj - rkj(2)*r2ij) &
                        * (-2.0d0*rkj(2)*scp3 + (rkj(2)+rij(2))*scp2 + rlk(2)*scp1 - rlk(2)*r2kj) * dn2 &
                  - 2.0 * ( rlk(2)*scp2-rkj(2)*r2lk)*(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj)*dn3 &
                  + (2.0 * scp3 - 2.0d0*scp2 + 4.0*rkj(2)*rlk(2) - 2.0d0*(rkj(2)+rij(2))*rlk(2)) * dn1
        lh(3,3) =  - (-2.0d0*scp1 + r2kj + r2ij + 4.0*rij(3)*rkj(3) - (rkj(3)+rij(3))**2 )*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * ( (rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)**2 * (r2kj*scp3-scp1*scp2)*dn4 &
                  + 2.0 * ( (rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)    * (rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2)*dn5 &
                  - (r2lk-rlk(3)**2)*(r2kj*scp3-scp1*scp2)*dn3 &
                  + 3.0 * (rlk(3)*scp2 - rkj(3)*r2lk)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * ( (rkj(3)+rij(3))*scp1 - rij(3)*r2kj - rkj(3)*r2ij) &
                        * (-2.0d0*rkj(3)*scp3 + (rkj(3)+rij(3))*scp2 + rlk(3)*scp1 - rlk(3)*r2kj) * dn2 &
                  - 2.0 * ( rlk(3)*scp2-rkj(3)*r2lk)*(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj)*dn3 &
                  + (2.0 * scp3 - 2.0d0*scp2 + 4.0*rkj(3)*rlk(3) - 2.0d0*(rkj(3)+rij(3))*rlk(3)) * dn1

        lh(1,2) = - (2.0d0*rij(1)*rkj(2)-(rkj(1)+rij(1))*(rkj(2)+rij(2))+2.0d0*rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                       *((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + rlk(1)*rlk(2)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(1)*scp2-rkj(1)*r2lk)*(rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  -  ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn2 &
                  - (rlk(1)*scp2-rkj(1)*r2lk)*(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn3 &
                  - ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn2 &
                  - (rlk(2)*scp2-rkj(2)*r2lk)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn3 &
                  + (2.0d0*rkj(1)*rlk(2)-(rkj(1)+rij(1))*rlk(2)+2.0d0*rkj(2)*rlk(1)-(rkj(2)+rij(2))*rlk(1)) * dn1
        lh(1,3) = - (2.0d0*rij(1)*rkj(3)-(rkj(1)+rij(1))*(rkj(3)+rij(3))+2.0d0*rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                       *((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + rlk(1)*rlk(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(1)*scp2-rkj(1)*r2lk)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn2 &
                  - (rlk(1)*scp2-rkj(1)*r2lk)*(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn3 &
                  - ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn2 &
                  - (rlk(3)*scp2-rkj(3)*r2lk)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn3 &
                  + (2.0d0*rkj(1)*rlk(3)-(rkj(1)+rij(1))*rlk(3)+2.0d0*rkj(3)*rlk(1)-(rkj(3)+rij(3))*rlk(1)) * dn1
        lh(2,3) = - (2.0d0*rij(2)*rkj(3)-(rkj(2)+rij(2))*(rkj(3)+rij(3))+2.0d0*rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                       *((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + rlk(2)*rlk(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(2)*scp2-rkj(2)*r2lk)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn2 &
                  - (rlk(2)*scp2-rkj(2)*r2lk)*(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn3 &
                  - ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn2 &
                  - (rlk(3)*scp2-rkj(3)*r2lk)*(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn3 &
                  + (2.0d0*rkj(2)*rlk(3)-(rkj(2)+rij(2))*rlk(3)+2.0d0*rkj(3)*rlk(2)-(rkj(3)+rij(3))*rlk(2)) * dn1

        geo%hess(1,j,1,j) = geo%hess(1,j,1,j) + f1*lh(1,1) + f2*dj(1)*dj(1)
        geo%hess(2,j,2,j) = geo%hess(2,j,2,j) + f1*lh(2,2) + f2*dj(2)*dj(2)
        geo%hess(3,j,3,j) = geo%hess(3,j,3,j) + f1*lh(3,3) + f2*dj(3)*dj(3)
        geo%hess(1,j,2,j) = geo%hess(1,j,2,j) + f1*lh(1,2) + f2*dj(1)*dj(2)
        geo%hess(1,j,3,j) = geo%hess(1,j,3,j) + f1*lh(1,3) + f2*dj(1)*dj(3)
        geo%hess(2,j,3,j) = geo%hess(2,j,3,j) + f1*lh(2,3) + f2*dj(2)*dj(3)
        geo%hess(2,j,1,j) = geo%hess(2,j,1,j) + f1*lh(1,2) + f2*dj(1)*dj(2)
        geo%hess(3,j,1,j) = geo%hess(3,j,1,j) + f1*lh(1,3) + f2*dj(1)*dj(3)
        geo%hess(3,j,2,j) = geo%hess(3,j,2,j) + f1*lh(2,3) + f2*dj(2)*dj(3)

        ! calculate hessian - diagonal k
        lh(1,1) = - (r2ij-rij(1)**2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(1)*r2ij-rij(1)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4  &
                  + 2.0d0*(rkj(1)*r2ij-rij(1)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2.0d0*scp2+r2lk+r2kj-(rlk(1)-rkj(1))**2-4.0*rkj(1)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0d0*(rkj(1)*r2ij-rij(1)*scp1)*(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn2 &
                  - 2.0d0*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn3 &
                  + (2.0d0*scp3+2.0d0*scp1-2.0d0*rij(1)*(rlk(1)-rkj(1))-4.0*rij(1)*rkj(1)) * dn1
        lh(2,2) =  - (r2ij-rij(2)**2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(2)*r2ij-rij(2)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4  &
                  + 2.0d0*(rkj(2)*r2ij-rij(2)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2.0d0*scp2+r2lk+r2kj-(rlk(2)-rkj(2))**2-4.0*rkj(2)*rlk(2))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0d0*(rkj(2)*r2ij-rij(2)*scp1)*(2.0d0*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1-rij(2)*r2kj) * dn2 &
                  - 2.0d0*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                       *(2.0d0*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1-rij(2)*r2kj) * dn3 &
                  + (2.0d0*scp3+2.0d0*scp1-2.0d0*rij(2)*(rlk(2)-rkj(2))-4.0*rij(2)*rkj(2)) * dn1
        lh(3,3) =   - (r2ij-rij(3)**2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(3)*r2ij-rij(3)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4  &
                  + 2.0d0*(rkj(3)*r2ij-rij(3)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2.0d0*scp2+r2lk+r2kj-(rlk(3)-rkj(3))**2-4.0*rkj(3)*rlk(3))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0d0*(rkj(3)*r2ij-rij(3)*scp1)*(2.0d0*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1-rij(3)*r2kj) * dn2 &
                  - 2.0d0*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                       *(2.0d0*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1-rij(3)*r2kj) * dn3 &
                  + (2.0d0*scp3+2.0d0*scp1-2.0d0*rij(3)*(rlk(3)-rkj(3))-4.0*rij(3)*rkj(3)) * dn1

        lh(1,2) = + rij(1)*rij(2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rkj(2)*r2ij-rij(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(2)-(rlk(1)-rkj(1))*(rlk(2)-rkj(2))-2*rkj(2)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                     *(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (-rij(1)*(rlk(2)-rkj(2))-rij(2)*(rlk(1)-rkj(1))+2*(-rij(1))*rkj(2)+2*(-rij(2))*rkj(1))*dn1
        lh(1,3) = + rij(1)*rij(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rkj(3)*r2ij-rij(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(3)-(rlk(1)-rkj(1))*(rlk(3)-rkj(3))-2*rkj(3)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                     *(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (-rij(1)*(rlk(3)-rkj(3))-rij(3)*(rlk(1)-rkj(1))+2*(-rij(1))*rkj(3)+2*(-rij(3))*rkj(1))*dn1
        lh(2,3) = + rij(2)*rij(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rkj(3)*r2ij-rij(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(2)*rlk(3)-(rlk(2)-rkj(2))*(rlk(3)-rkj(3))-2*rkj(3)*rlk(2))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                     *(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (-rij(2)*(rlk(3)-rkj(3))-rij(3)*(rlk(2)-rkj(2))+2*(-rij(2))*rkj(3)+2*(-rij(3))*rkj(2))*dn1

        geo%hess(1,k,1,k) = geo%hess(1,k,1,k) + f1*lh(1,1) + f2*dk(1)*dk(1)
        geo%hess(2,k,2,k) = geo%hess(2,k,2,k) + f1*lh(2,2) + f2*dk(2)*dk(2)
        geo%hess(3,k,3,k) = geo%hess(3,k,3,k) + f1*lh(3,3) + f2*dk(3)*dk(3)
        geo%hess(1,k,2,k) = geo%hess(1,k,2,k) + f1*lh(1,2) + f2*dk(1)*dk(2)
        geo%hess(1,k,3,k) = geo%hess(1,k,3,k) + f1*lh(1,3) + f2*dk(1)*dk(3)
        geo%hess(2,k,3,k) = geo%hess(2,k,3,k) + f1*lh(2,3) + f2*dk(2)*dk(3)
        geo%hess(2,k,1,k) = geo%hess(2,k,1,k) + f1*lh(1,2) + f2*dk(1)*dk(2)
        geo%hess(3,k,1,k) = geo%hess(3,k,1,k) + f1*lh(1,3) + f2*dk(1)*dk(3)
        geo%hess(3,k,2,k) = geo%hess(3,k,2,k) + f1*lh(2,3) + f2*dk(2)*dk(3)

        ! calculate hessian - diagonal l
        lh(1,1) = - (r2kj-rkj(1)**2) * (r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(1)*r2kj-rkj(1)*scp2)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * (r2kj*rij(1)-rkj(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(2,2) = - (r2kj-rkj(2)**2) * (r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(2)*r2kj-rkj(2)*scp2)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * (r2kj*rij(2)-rkj(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(3,3) = - (r2kj-rkj(3)**2) * (r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(3)*r2kj-rkj(3)*scp2)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * (r2kj*rij(3)-rkj(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3
        lh(1,2) =   rkj(1)*rkj(2)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (r2kj*rij(1)-rkj(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(1,3) = rkj(1)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (r2kj*rij(1)-rkj(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(2,3) =  rkj(2)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (r2kj*rij(2)-rkj(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3

        geo%hess(1,l,1,l) = geo%hess(1,l,1,l) + f1*lh(1,1) + f2*dl(1)*dl(1)
        geo%hess(2,l,2,l) = geo%hess(2,l,2,l) + f1*lh(2,2) + f2*dl(2)*dl(2)
        geo%hess(3,l,3,l) = geo%hess(3,l,3,l) + f1*lh(3,3) + f2*dl(3)*dl(3)
        geo%hess(1,l,2,l) = geo%hess(1,l,2,l) + f1*lh(1,2) + f2*dl(1)*dl(2)
        geo%hess(1,l,3,l) = geo%hess(1,l,3,l) + f1*lh(1,3) + f2*dl(1)*dl(3)
        geo%hess(2,l,3,l) = geo%hess(2,l,3,l) + f1*lh(2,3) + f2*dl(2)*dl(3)
        geo%hess(2,l,1,l) = geo%hess(2,l,1,l) + f1*lh(1,2) + f2*dl(1)*dl(2)
        geo%hess(3,l,1,l) = geo%hess(3,l,1,l) + f1*lh(1,3) + f2*dl(1)*dl(3)
        geo%hess(3,l,2,l) = geo%hess(3,l,2,l) + f1*lh(2,3) + f2*dl(2)*dl(3)

        ! calculate hessian - off diagonal i,j
        lh(1,1) = - (scp1-r2kj-(-rkj(1)-rij(1))*rkj(1)-2*rij(1)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + ( rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(1)*rlk(1)-rkj(1)*(-rlk(1)) ) * dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(2,2) = - (scp1-r2kj-(-rkj(2)-rij(2))*rkj(2)-2*rij(2)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + ( rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(2)*rlk(2)-rkj(1)*(-rlk(1)) ) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(3,3) = - (scp1-r2kj-(-rkj(3)-rij(3))*rkj(3)-2*rij(3)*rkj(3))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + ( rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(3)*rlk(3)-rkj(1)*(-rlk(1)) ) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        lh(1,2) = - (-2*rij(1)*rkj(2)-rkj(1)*(-rkj(2)-rij(2)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + ( -rkj(1)*(-rlk(2))-2*rkj(2)*rlk(1) ) * dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(1,3) = - (-2*rij(1)*rkj(3)-rkj(1)*(-rkj(3)-rij(3)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + ( -rkj(1)*(-rlk(3))-2*rkj(3)*rlk(1) ) * dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3
        lh(2,3) = - (-2*rij(2)*rkj(3)-rkj(2)*(-rkj(3)-rij(3)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rij(2)*r2kj-rkj(2)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + ( -rkj(2)*(-rlk(3))-2*rkj(3)*rlk(2) ) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        lh(2,1) = - (-(-rkj(1)-rij(1))*rkj(2)-2*rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + (-2*rkj(1)*rlk(2)-rkj(2)*(-rlk(1))) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(3,1) = - (-(-rkj(1)-rij(1))*rkj(3)-2*rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (-2*rkj(1)*rlk(3)-rkj(3)*(-rlk(1))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(3,2) = - (-(-rkj(2)-rij(2))*rkj(3)-2*rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (-2*rkj(2)*rlk(3)-rkj(3)*(-rlk(2))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dj(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dj(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dj(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dj(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dj(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dj(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dj(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dj(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dj(2)

        geo%hess(1,i,1,j) = geo%hess(1,i,1,j) + hh(1,1)
        geo%hess(2,i,2,j) = geo%hess(2,i,2,j) + hh(2,2)
        geo%hess(3,i,3,j) = geo%hess(3,i,3,j) + hh(3,3)
        geo%hess(1,i,2,j) = geo%hess(1,i,2,j) + hh(1,2)
        geo%hess(1,i,3,j) = geo%hess(1,i,3,j) + hh(1,3)
        geo%hess(2,i,3,j) = geo%hess(2,i,3,j) + hh(2,3)
        geo%hess(2,i,1,j) = geo%hess(2,i,1,j) + hh(2,1)
        geo%hess(3,i,1,j) = geo%hess(3,i,1,j) + hh(3,1)
        geo%hess(3,i,2,j) = geo%hess(3,i,2,j) + hh(3,2)

        geo%hess(1,j,1,i) = geo%hess(1,j,1,i) + hh(1,1)
        geo%hess(2,j,2,i) = geo%hess(2,j,2,i) + hh(2,2)
        geo%hess(3,j,3,i) = geo%hess(3,j,3,i) + hh(3,3)
        geo%hess(1,j,2,i) = geo%hess(1,j,2,i) + hh(2,1)
        geo%hess(1,j,3,i) = geo%hess(1,j,3,i) + hh(3,1)
        geo%hess(2,j,3,i) = geo%hess(2,j,3,i) + hh(3,2)
        geo%hess(2,j,1,i) = geo%hess(2,j,1,i) + hh(1,2)
        geo%hess(3,j,1,i) = geo%hess(3,j,1,i) + hh(1,3)
        geo%hess(3,j,2,i) = geo%hess(3,j,2,i) + hh(2,3)

        ! calculate hessian - off diagonal i,k
        lh(1,1) =  - (rij(1)*rkj(1)-scp1)*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (rij(1)*r2kj-rkj(1)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                   - (rkj(1)*r2ij-rij(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                   + (-rkj(3)*rlk(3)-rkj(2)*rlk(2)+rkj(1)*rlk(1)-rkj(1)*(rlk(1)-rkj(1))-rkj(3)**2-rkj(2)**2-rkj(1)**2)*dn1 &
                   - (rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(2,2) = - (rij(2)*rkj(2)-scp1)*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (rij(2)*r2kj-rkj(2)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                   - (rkj(2)*r2ij-rij(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                   + (-rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(2)*(rlk(2)-rkj(2))-rkj(1)*rlk(1)-rkj(3)**2-rkj(2)**2-rkj(1)**2)*dn1 &
                   - (rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(3,3) = - (rij(3)*rkj(3)-scp1)*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(3)*r2ij-rij(3)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (rij(3)*r2kj-rkj(3)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                   - (rkj(3)*r2ij-rij(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                   + (rkj(3)*rlk(3)-rkj(3)*(rlk(3)-rkj(3))-rkj(2)*rlk(2)-rkj(1)*rlk(1)-rkj(3)**2-rkj(2)**2-rkj(1)**2)*dn1 &
                   - (rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        lh(1,2) = - (2*rij(1)*rkj(2)-rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + (2*rkj(2)*rlk(1)-rkj(1)*(rlk(2)-rkj(2)))* dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(1,3) = - (2*rij(1)*rkj(3)-rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + (2*rkj(3)*rlk(1)-rkj(1)*(rlk(3)-rkj(3)))* dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3
        lh(2,3) = - (2*rij(2)*rkj(3)-rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + (2*rkj(3)*rlk(2)-rkj(2)*(rlk(3)-rkj(3)))* dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        lh(2,1) = - (2*rij(2)*rkj(1)-rij(1)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + (2*rkj(1)*rlk(2)-rkj(2)*(rlk(1)-rkj(1))) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(3,1) = - (2*rij(3)*rkj(1)-rij(1)*rkj(3))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (2*rkj(1)*rlk(3)-rkj(3)*(rlk(1)-rkj(1))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(3,2) = - (2*rij(3)*rkj(2)-rij(2)*rkj(3))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (2*rkj(2)*rlk(3)-rkj(3)*(rlk(2)-rkj(2))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dk(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dk(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dk(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dk(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dk(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dk(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dk(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dk(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dk(2)

        geo%hess(1,i,1,k) = geo%hess(1,i,1,k) + hh(1,1)
        geo%hess(2,i,2,k) = geo%hess(2,i,2,k) + hh(2,2)
        geo%hess(3,i,3,k) = geo%hess(3,i,3,k) + hh(3,3)
        geo%hess(1,i,2,k) = geo%hess(1,i,2,k) + hh(1,2)
        geo%hess(1,i,3,k) = geo%hess(1,i,3,k) + hh(1,3)
        geo%hess(2,i,3,k) = geo%hess(2,i,3,k) + hh(2,3)
        geo%hess(2,i,1,k) = geo%hess(2,i,1,k) + hh(2,1)
        geo%hess(3,i,1,k) = geo%hess(3,i,1,k) + hh(3,1)
        geo%hess(3,i,2,k) = geo%hess(3,i,2,k) + hh(3,2)

        geo%hess(1,k,1,i) = geo%hess(1,k,1,i) + hh(1,1)
        geo%hess(2,k,2,i) = geo%hess(2,k,2,i) + hh(2,2)
        geo%hess(3,k,3,i) = geo%hess(3,k,3,i) + hh(3,3)
        geo%hess(1,k,2,i) = geo%hess(1,k,2,i) + hh(2,1)
        geo%hess(1,k,3,i) = geo%hess(1,k,3,i) + hh(3,1)
        geo%hess(2,k,3,i) = geo%hess(2,k,3,i) + hh(3,2)
        geo%hess(2,k,1,i) = geo%hess(2,k,1,i) + hh(1,2)
        geo%hess(3,k,1,i) = geo%hess(3,k,1,i) + hh(1,3)
        geo%hess(3,k,2,i) = geo%hess(3,k,2,i) + hh(2,3)

        ! calculate hessian - off diagonal i,l
        lh(1,1) = + (rij(1)*r2kj-rkj(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)**2 + rkj(2)**2) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(2,2) = + (rij(2)*r2kj-rkj(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)**2 + rkj(1)**2) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(3,3) = + (rij(3)*r2kj-rkj(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(2)**2 + rkj(1)**2) * dn1 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3

        lh(1,2) = + (rij(1)*r2kj-rkj(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(2) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(1,3) = + (rij(1)*r2kj-rkj(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(3) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3
        lh(2,3) =  + (rij(2)*r2kj-rkj(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(2)*rkj(3) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3

        lh(2,1) = + (rij(2)*r2kj-rkj(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(2) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(3,1) = + (rij(3)*r2kj-rkj(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(3) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3
        lh(3,2) = + (rij(3)*r2kj-rkj(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(2)*rkj(3) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dl(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dl(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dl(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dl(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dl(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dl(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dl(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dl(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dl(2)

        geo%hess(1,i,1,l) = geo%hess(1,i,1,l) + hh(1,1)
        geo%hess(2,i,2,l) = geo%hess(2,i,2,l) + hh(2,2)
        geo%hess(3,i,3,l) = geo%hess(3,i,3,l) + hh(3,3)
        geo%hess(1,i,2,l) = geo%hess(1,i,2,l) + hh(1,2)
        geo%hess(1,i,3,l) = geo%hess(1,i,3,l) + hh(1,3)
        geo%hess(2,i,3,l) = geo%hess(2,i,3,l) + hh(2,3)
        geo%hess(2,i,1,l) = geo%hess(2,i,1,l) + hh(2,1)
        geo%hess(3,i,1,l) = geo%hess(3,i,1,l) + hh(3,1)
        geo%hess(3,i,2,l) = geo%hess(3,i,2,l) + hh(3,2)

        geo%hess(1,l,1,i) = geo%hess(1,l,1,i) + hh(1,1)
        geo%hess(2,l,2,i) = geo%hess(2,l,2,i) + hh(2,2)
        geo%hess(3,l,3,i) = geo%hess(3,l,3,i) + hh(3,3)
        geo%hess(1,l,2,i) = geo%hess(1,l,2,i) + hh(2,1)
        geo%hess(1,l,3,i) = geo%hess(1,l,3,i) + hh(3,1)
        geo%hess(2,l,3,i) = geo%hess(2,l,3,i) + hh(3,2)
        geo%hess(2,l,1,i) = geo%hess(2,l,1,i) + hh(1,2)
        geo%hess(3,l,1,i) = geo%hess(3,l,1,i) + hh(1,3)
        geo%hess(3,l,2,i) = geo%hess(3,l,2,i) + hh(2,3)

        ! calculate hessian - off diagonal j,k
        lh(1,1) = - (scp1-r2ij-2.0d0*rij(1)*rkj(1)+rij(1)*(rkj(1)+rij(1)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(1)*r2ij-rij(1)*scp1)*((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2lk+2.0d0*rkj(1)*rlk(1)+rlk(1)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(1)*scp2-rkj(1)*r2lk)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn2 &
                  - (rlk(1)*scp2-rkj(1)*r2lk)*(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn3 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn3 &
                  + (-2.0d0*scp3+rkj(3)*rlk(3)+rkj(2)*rlk(2)+rkj(1)*rlk(1) + (rkj(1)+rij(1))*(rlk(1)-rkj(1)) &
                     -2.0d0*rkj(1)*rlk(1)+rij(1)*rlk(1)+rkj(3)**2 &
                     -rij(3)*rkj(3) + rkj(2)**2 - rij(2)*rkj(2) + rkj(1)**2 + 2.0d0*rij(1)*rkj(1) - rij(1)*rkj(1)) * dn1
        lh(2,2) =  - (scp1-r2ij-2*rij(2)*rkj(2)-rij(2)*(-rkj(2) - rij(2)))*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(2)*r2ij-rij(2)*scp1)*(-(-rkj(2) - rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (-(-rkj(2) - rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                     *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   + (rkj(2)*r2ij-rij(2)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (-scp2-r2lk+2*rkj(2)*rlk(2)-(-rlk(2))*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                   + 3*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                   - (-(-rkj(2) - rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                     *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                   - (-(-rlk(2))*scp2-rkj(2)*r2lk)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                   - (rkj(2)*r2ij-rij(2)*scp1)*(-2*rkj(2)*scp3-(-rkj(2) - rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                   - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                     *(-2*rkj(2)*scp3-(-rkj(2) - rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                   + (-2*scp3+rkj(3)*rlk(3)+rkj(2)*rlk(2)-(-rkj(2) - rij(2))*(rlk(2)-rkj(2)) &
                      +2*rkj(2)*(-rlk(2))-rij(2)*(-rlk(2))+rkj(1)*rlk(1)+rkj(3)**2-rij(3)*rkj(3) &
                      +rkj(2)**2-2*(-rij(2))*rkj(2)-rij(2)*rkj(2)+rkj(1)**2-rij(1)*rkj(1)) * dn1
        lh(3,3) = - (scp1-r2ij-2*rij(3)*rkj(3)-rij(3)*(-rkj(3)-rij(3)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2lk+2*rkj(3)*rlk(3)-(-rlk(3))*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(-rlk(3))*scp2-rkj(3)*r2lk)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-2*scp3+rkj(3)*rlk(3)-(-rkj(3)-rij(3))*(rlk(3)-rkj(3))+2*rkj(3)*(-rlk(3))-rij(3)*(-rlk(3)) &
                     +rkj(2)*rlk(2)+rkj(1)*rlk(1)+rkj(3)**2-2*(-rij(3))*rkj(3)-rij(3)*rkj(3) &
                     +rkj(2)**2-rij(2)*rkj(2)+rkj(1)**2-rij(1)*rkj(1)) * dn1

        lh(1,2) = - (-2*rij(1)*rkj(2)-rij(2)*(-rkj(1)-rij(1)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                     *(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(2)-(-rlk(1))*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                     *(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(-rlk(1))*scp2-rkj(1)*r2lk)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*(rlk(2)-rkj(2))+2*rkj(2)*(-rlk(1))-rij(2)*(-rlk(1))-2*(-rij(2))*rkj(1)) * dn1
        lh(1,3) = - (-2*rij(1)*rkj(3)-rij(3)*(-rkj(1)-rij(1)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(3)-(-rlk(1))*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                     *(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(-rlk(1))*scp2-rkj(1)*r2lk)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*(rlk(3)-rkj(3))+2*rkj(3)*(-rlk(1))-rij(3)*(-rlk(1))-2*(-rij(3))*rkj(1)) * dn1
        lh(2,3) = - (-2*rij(2)*rkj(3)-rij(3)*(-rkj(2)-rij(2)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                     *(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(3)-(-rlk(2))*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(-rlk(2))*scp2-rkj(2)*r2lk)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-(-rkj(2)-rij(2))*(rlk(3)-rkj(3))+2*rkj(3)*(-rlk(2))-rij(3)*(-rlk(2))-2*(-rij(3))*rkj(2)) * dn1

        lh(2,1) = - (-rij(1)*(-rkj(2)-rij(2))-2*rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(1)-(rlk(1)-rkj(1))*(-rlk(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(-rlk(2))*scp2-rkj(2)*r2lk)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rkj(1)*(-rlk(2))-rij(1)*(-rlk(2))-(-rkj(2)-rij(2))*(rlk(1)-rkj(1))-2*(-rij(1))*rkj(2)) * dn1
        lh(3,1) = - (-rij(1)*(-rkj(3)-rij(3))-2*rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(1)-(rlk(1)-rkj(1))*(-rlk(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(-rlk(3))*scp2-rkj(3)*r2lk)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rkj(1)*(-rlk(3))-rij(1)*(-rlk(3))-(-rkj(3)-rij(3))*(rlk(1)-rkj(1))-2*(-rij(1))*rkj(3)) * dn1
        lh(3,2) =  - (-rij(2)*(-rkj(3)-rij(3))-2*rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(2)-(rlk(2)-rkj(2))*(-rlk(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(-rlk(3))*scp2-rkj(3)*r2lk)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (2*rkj(2)*(-rlk(3))-rij(2)*(-rlk(3))-(-rkj(3)-rij(3))*(rlk(2)-rkj(2))-2*(-rij(2))*rkj(3)) * dn1

        hh(1,1) = f1*lh(1,1) + f2*dj(1)*dk(1)
        hh(2,2) = f1*lh(2,2) + f2*dj(2)*dk(2)
        hh(3,3) = f1*lh(3,3) + f2*dj(3)*dk(3)

        hh(1,2) = f1*lh(1,2) + f2*dj(1)*dk(2)
        hh(1,3) = f1*lh(1,3) + f2*dj(1)*dk(3)
        hh(2,3) = f1*lh(2,3) + f2*dj(2)*dk(3)

        hh(2,1) = f1*lh(2,1) + f2*dj(2)*dk(1)
        hh(3,1) = f1*lh(3,1) + f2*dj(3)*dk(1)
        hh(3,2) = f1*lh(3,2) + f2*dj(3)*dk(2)

        geo%hess(1,j,1,k) = geo%hess(1,j,1,k) + hh(1,1)
        geo%hess(2,j,2,k) = geo%hess(2,j,2,k) + hh(2,2)
        geo%hess(3,j,3,k) = geo%hess(3,j,3,k) + hh(3,3)
        geo%hess(1,j,2,k) = geo%hess(1,j,2,k) + hh(1,2)
        geo%hess(1,j,3,k) = geo%hess(1,j,3,k) + hh(1,3)
        geo%hess(2,j,3,k) = geo%hess(2,j,3,k) + hh(2,3)
        geo%hess(2,j,1,k) = geo%hess(2,j,1,k) + hh(2,1)
        geo%hess(3,j,1,k) = geo%hess(3,j,1,k) + hh(3,1)
        geo%hess(3,j,2,k) = geo%hess(3,j,2,k) + hh(3,2)

        geo%hess(1,k,1,j) = geo%hess(1,k,1,j) + hh(1,1)
        geo%hess(2,k,2,j) = geo%hess(2,k,2,j) + hh(2,2)
        geo%hess(3,k,3,j) = geo%hess(3,k,3,j) + hh(3,3)
        geo%hess(1,k,2,j) = geo%hess(1,k,2,j) + hh(2,1)
        geo%hess(1,k,3,j) = geo%hess(1,k,3,j) + hh(3,1)
        geo%hess(2,k,3,j) = geo%hess(2,k,3,j) + hh(3,2)
        geo%hess(2,k,1,j) = geo%hess(2,k,1,j) + hh(1,2)
        geo%hess(3,k,1,j) = geo%hess(3,k,1,j) + hh(1,3)
        geo%hess(3,k,2,j) = geo%hess(3,k,2,j) + hh(2,3)

        ! calculate hessian - off diagonal j,l
        lh(1,1) = + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (scp2-2*rkj(1)*rlk(1)-rkj(1)*(-rlk(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-rkj(3)**2+rij(3)*rkj(3)-rkj(2)**2+rij(2)*rkj(2)-rkj(1)**2-(-rkj(1)-rij(1))*rkj(1)-rij(1)*rkj(1)) * dn1 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(2,2) = + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (scp2-2*rkj(2)*rlk(2)-rkj(2)*(-rlk(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-rkj(3)**2+rij(3)*rkj(3)-rkj(2)**2-(-rkj(2)-rij(2))*rkj(2)-rij(2)*rkj(2)-rkj(1)**2+rij(1)*rkj(1)) * dn1 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(3,3) = + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (scp2-2*rkj(3)*rlk(3)-rkj(3)*(-rlk(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-rkj(3)**2-(-rkj(3)-rij(3))*rkj(3)-rij(3)*rkj(3)-rkj(2)**2+rij(2)*rkj(2)-rkj(1)**2+rij(1)*rkj(1)) * dn1 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        lh(1,2) = + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(2)-rkj(2)*(-rlk(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*rkj(2)-2*rij(2)*rkj(1)) * dn1 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(1,3) = + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(3)-rkj(3)*(-rlk(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*rkj(3)-2*rij(3)*rkj(1)) * dn1 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(2,3) = + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(2)*rlk(3)-rkj(3)*(-rlk(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-(-rkj(2)-rij(2))*rkj(3)-2*rij(3)*rkj(2)) * dn1 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3

        lh(2,1) = + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-rkj(1)*(-rlk(2))-2*rkj(2)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-2*rij(1)*rkj(2)-rkj(1)*(-rkj(2)-rij(2))) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(3,1) = + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-rkj(1)*(-rlk(3))-2*rkj(3)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-2*rij(1)*rkj(3)-rkj(1)*(-rkj(3)-rij(3))) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3
        lh(3,2) = + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-rkj(2)*(-rlk(3))-2*rkj(3)*rlk(2))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-2*rij(2)*rkj(3)-rkj(2)*(-rkj(3)-rij(3))) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        hh(1,1) = f1*lh(1,1) + f2*dj(1)*dl(1)
        hh(2,2) = f1*lh(2,2) + f2*dj(2)*dl(2)
        hh(3,3) = f1*lh(3,3) + f2*dj(3)*dl(3)

        hh(1,2) = f1*lh(1,2) + f2*dj(1)*dl(2)
        hh(1,3) = f1*lh(1,3) + f2*dj(1)*dl(3)
        hh(2,3) = f1*lh(2,3) + f2*dj(2)*dl(3)

        hh(2,1) = f1*lh(2,1) + f2*dj(2)*dl(1)
        hh(3,1) = f1*lh(3,1) + f2*dj(3)*dl(1)
        hh(3,2) = f1*lh(3,2) + f2*dj(3)*dl(2)

        geo%hess(1,j,1,l) = geo%hess(1,j,1,l) + hh(1,1)
        geo%hess(2,j,2,l) = geo%hess(2,j,2,l) + hh(2,2)
        geo%hess(3,j,3,l) = geo%hess(3,j,3,l) + hh(3,3)
        geo%hess(1,j,2,l) = geo%hess(1,j,2,l) + hh(1,2)
        geo%hess(1,j,3,l) = geo%hess(1,j,3,l) + hh(1,3)
        geo%hess(2,j,3,l) = geo%hess(2,j,3,l) + hh(2,3)
        geo%hess(2,j,1,l) = geo%hess(2,j,1,l) + hh(2,1)
        geo%hess(3,j,1,l) = geo%hess(3,j,1,l) + hh(3,1)
        geo%hess(3,j,2,l) = geo%hess(3,j,2,l) + hh(3,2)

        geo%hess(1,l,1,j) = geo%hess(1,l,1,j) + hh(1,1)
        geo%hess(2,l,2,j) = geo%hess(2,l,2,j) + hh(2,2)
        geo%hess(3,l,3,j) = geo%hess(3,l,3,j) + hh(3,3)
        geo%hess(1,l,2,j) = geo%hess(1,l,2,j) + hh(2,1)
        geo%hess(1,l,3,j) = geo%hess(1,l,3,j) + hh(3,1)
        geo%hess(2,l,3,j) = geo%hess(2,l,3,j) + hh(3,2)
        geo%hess(2,l,1,j) = geo%hess(2,l,1,j) + hh(1,2)
        geo%hess(3,l,1,j) = geo%hess(3,l,1,j) + hh(1,3)
        geo%hess(3,l,2,j) = geo%hess(3,l,2,j) + hh(2,3)

        ! calculate hessian - off diagonal k,l
        lh(1,1) = + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2kj+2*rkj(1)*rlk(1)-rkj(1)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (-rij(3)*rkj(3)-rij(2)*rkj(2)) * dn1 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(2,2) = + (rkj(2)*r2ij-rij(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2kj+2*rkj(2)*rlk(2)-rkj(2)*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (-rij(3)*rkj(3)-rij(1)*rkj(1)) * dn1 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(3,3) = + (rkj(3)*r2ij-rij(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2kj+2*rkj(3)*rlk(3)-rkj(3)*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  + (-rij(2)*rkj(2)-rij(1)*rkj(1)) * dn1 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        lh(1,2) = + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(2)-rkj(2)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rij(2)*rkj(1)-rij(1)*rkj(2)) * dn1 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(1,3) = + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(3)-rkj(3)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rij(3)*rkj(1)-rij(1)*rkj(3)) * dn1 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(2,3) = + (rkj(2)*r2ij-rij(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(3)-rkj(3)*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (2*rij(3)*rkj(2)-rij(2)*rkj(3)) * dn1 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3

        lh(2,1) = + (rkj(2)*r2ij-rij(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(1)-rkj(1)*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) *dn3 &
                  + (2*rij(1)*rkj(2)-rij(2)*rkj(1)) * dn1 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(3,1) = + (rkj(3)*r2ij-rij(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(1)-rkj(1)*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) *dn3 &
                  + (2*rij(1)*rkj(3)-rij(3)*rkj(1)) * dn1 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3
        lh(3,2) = + (rkj(3)*r2ij-rij(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(2)-rkj(2)*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) *dn3 &
                  + (2*rij(2)*rkj(3)-rij(3)*rkj(2)) * dn1 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        hh(1,1) = f1*lh(1,1) + f2*dk(1)*dl(1)
        hh(2,2) = f1*lh(2,2) + f2*dk(2)*dl(2)
        hh(3,3) = f1*lh(3,3) + f2*dk(3)*dl(3)

        hh(1,2) = f1*lh(1,2) + f2*dk(1)*dl(2)
        hh(1,3) = f1*lh(1,3) + f2*dk(1)*dl(3)
        hh(2,3) = f1*lh(2,3) + f2*dk(2)*dl(3)

        hh(2,1) = f1*lh(2,1) + f2*dk(2)*dl(1)
        hh(3,1) = f1*lh(3,1) + f2*dk(3)*dl(1)
        hh(3,2) = f1*lh(3,2) + f2*dk(3)*dl(2)

        geo%hess(1,k,1,l) = geo%hess(1,k,1,l) + hh(1,1)
        geo%hess(2,k,2,l) = geo%hess(2,k,2,l) + hh(2,2)
        geo%hess(3,k,3,l) = geo%hess(3,k,3,l) + hh(3,3)
        geo%hess(1,k,2,l) = geo%hess(1,k,2,l) + hh(1,2)
        geo%hess(1,k,3,l) = geo%hess(1,k,3,l) + hh(1,3)
        geo%hess(2,k,3,l) = geo%hess(2,k,3,l) + hh(2,3)
        geo%hess(2,k,1,l) = geo%hess(2,k,1,l) + hh(2,1)
        geo%hess(3,k,1,l) = geo%hess(3,k,1,l) + hh(3,1)
        geo%hess(3,k,2,l) = geo%hess(3,k,2,l) + hh(3,2)

        geo%hess(1,l,1,k) = geo%hess(1,l,1,k) + hh(1,1)
        geo%hess(2,l,2,k) = geo%hess(2,l,2,k) + hh(2,2)
        geo%hess(3,l,3,k) = geo%hess(3,l,3,k) + hh(3,3)
        geo%hess(1,l,2,k) = geo%hess(1,l,2,k) + hh(2,1)
        geo%hess(1,l,3,k) = geo%hess(1,l,3,k) + hh(3,1)
        geo%hess(2,l,3,k) = geo%hess(2,l,3,k) + hh(3,2)
        geo%hess(2,l,1,k) = geo%hess(2,l,1,k) + hh(1,2)
        geo%hess(3,l,1,k) = geo%hess(3,l,1,k) + hh(1,3)
        geo%hess(3,l,2,k) = geo%hess(3,l,2,k) + hh(2,3)

    end do

end subroutine ffdev_hessian_dihedrals

!===============================================================================
! subroutine ffdev_hessian_impropers
!===============================================================================

subroutine ffdev_hessian_impropers(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip
    real(DEVDP)     ::  scp1,scp2,scp3,r2ij,r2kj,r2lk,sarg,carg,dn1,dn2,dn3,dn4,dn5,dn6
    real(DEVDP)     ::  tp,y,phi,f0,f1,f2,arg,th1,th2,ibo1,ibo2
    real(DEVDP)     ::  rij(3),rkj(3),rlk(3),rnj(3),rnk(3)
    real(DEVDP)     ::  di(3),dj(3),dk(3),dl(3),lh(3,3),hh(3,3)
    ! -----------------------------------------------------------------------------

    geo%impropr_ene = 0.0d0

    do ip=1,top%nimpropers

        i  = top%impropers(ip)%ai
        j  = top%impropers(ip)%aj
        k  = top%impropers(ip)%ak
        l  = top%impropers(ip)%al

        ic = top%impropers(ip)%dt

        rij(:) = geo%crd(:,i) - geo%crd(:,j)
        rkj(:) = geo%crd(:,k) - geo%crd(:,j)
        rlk(:) = geo%crd(:,l) - geo%crd(:,k)

        rnj(1) =  rij(2)*rkj(3) - rij(3)*rkj(2)
        rnj(2) =  rij(3)*rkj(1) - rij(1)*rkj(3)
        rnj(3) =  rij(1)*rkj(2) - rij(2)*rkj(1)

        rnk(1) = -rkj(2)*rlk(3) + rkj(3)*rlk(2)
        rnk(2) = -rkj(3)*rlk(1) + rkj(1)*rlk(3)
        rnk(3) = -rkj(1)*rlk(2) + rkj(2)*rlk(1)

        scp1 = rij(1)*rkj(1) + rij(2)*rkj(2) + rij(3)*rkj(3)
        scp2 = rkj(1)*rlk(1) + rkj(2)*rlk(2) + rkj(3)*rlk(3)
        scp3 = rij(1)*rlk(1) + rij(2)*rlk(2) + rij(3)*rlk(3)

        r2ij = rij(1)**2 + rij(2)**2 + rij(3)**2
        r2kj = rkj(1)**2 + rkj(2)**2 + rkj(3)**2
        r2lk = rlk(1)**2 + rlk(2)**2 + rlk(3)**2

        tp  = scp3*r2kj - scp1*scp2
        ibo1 = 1.0d0 / sqrt(r2ij*r2kj - scp1**2)
        ibo2 = 1.0d0 / sqrt(r2lk*r2kj - scp2**2)

        ! calculate y and phi
        y = tp * ibo1 * ibo2

        if ( y .gt.  1.0 ) then
                y =  1.0
                phi = acos (1.0) ! const
        else if ( y .lt. -1.0 ) then
                y = -1.0
                phi = acos (-1.0) ! const
        else
            phi = acos ( y )
        end if

        ! improper sign
        if( rkj(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
           -rkj(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
           +rkj(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) .lt. 0) then
                phi = -phi
        end if

        f0 = sin ( phi )
        if ( abs(f0) .lt. 1.e-12 ) f0 = 1.e-12
        f0 =  1.0d0 / f0

        f1 = 0.0d0
        f2 = 0.0d0

        ! calculate energy
        arg = 2.0d0*phi - top%improper_types(ic)%g
        sarg = sin(arg)
        carg = cos(arg)
        geo%impropr_ene = geo%impropr_ene + top%improper_types(ic)%v*(1.0d0+carg)

        ! outer derivatives
        f1 = f1 + 2.0d0*top%improper_types(ic)%v*sarg*f0
        f2 = f2 + 2.0    * top%improper_types(ic)%v * y * sarg * f0**3 &
                - 2.0d0**2 * top%improper_types(ic)%v     * carg * f0**2

        ! helpers
        dn1 = ibo1    * ibo2
        dn2 = ibo1**3 * ibo2
        dn3 = ibo1    * ibo2**3
        dn4 = ibo1**5 * ibo2
        dn5 = ibo1**3 * ibo2**3
        dn6 = ibo1    * ibo2**5

        ! gradients of y
        th1 = (r2kj*scp3-scp1*scp2)*dn2
        th2 = (r2kj*scp3-scp1*scp2)*dn3

        di(1) = dn1*(rlk(1)*r2kj-rkj(1)*scp2) - (rij(1)*r2kj-rkj(1)*scp1)*th1
        di(2) = dn1*(rlk(2)*r2kj-rkj(2)*scp2) - (rij(2)*r2kj-rkj(2)*scp1)*th1
        di(3) = dn1*(rlk(3)*r2kj-rkj(3)*scp2) - (rij(3)*r2kj-rkj(3)*scp1)*th1

        dj(1) = - ( (rkj(1) + rij(1))*scp1 - rij(1)*r2kj - rkj(1)*r2ij )*th1 &
                - ( rlk(1)*scp2 - rkj(1)*r2lk )*th2 &
                + ( (rkj(1)+rij(1))*scp2 + rlk(1)*scp1 - rlk(1)*r2kj - 2.0d0*rkj(1)*scp3)*dn1
        dj(2) = - ( (rkj(2) + rij(2))*scp1 - rij(2)*r2kj - rkj(2)*r2ij )*th1 &
                - ( rlk(2)*scp2 - rkj(2)*r2lk )*th2 &
                + ( (rkj(2)+rij(2))*scp2 + rlk(2)*scp1 - rlk(2)*r2kj - 2.0d0*rkj(2)*scp3)*dn1
        dj(3) = - ( (rkj(3) + rij(3))*scp1 - rij(3)*r2kj - rkj(3)*r2ij )*th1 &
                - ( rlk(3)*scp2 - rkj(3)*r2lk )*th2 &
                + ( (rkj(3)+rij(3))*scp2 + rlk(3)*scp1 - rlk(3)*r2kj - 2.0d0*rkj(3)*scp3)*dn1

        dk(1) = - ( rkj(1)*r2ij - rij(1)*scp1 )*th1 &
                - ( -(rlk(1) - rkj(1))*scp2 + rkj(1)*r2lk - rlk(1)*r2kj )*th2 &
                + ( 2.0d0*rkj(1)*scp3 - rij(1)*scp2 - (rlk(1) - rkj(1))*scp1 - rij(1)*r2kj )*dn1
        dk(2) = - ( rkj(2)*r2ij - rij(2)*scp1 )*th1 &
                - ( -(rlk(2) - rkj(2))*scp2 + rkj(2)*r2lk - rlk(2)*r2kj )*th2 &
                + ( 2.0d0*rkj(2)*scp3 - rij(2)*scp2 - (rlk(2) - rkj(2))*scp1 - rij(2)*r2kj )*dn1
        dk(3) = - ( rkj(3)*r2ij - rij(3)*scp1 )*th1 &
                - ( -(rlk(3) - rkj(3))*scp2 + rkj(3)*r2lk - rlk(3)*r2kj )*th2 &
                + ( 2.0d0*rkj(3)*scp3 - rij(3)*scp2 - (rlk(3) - rkj(3))*scp1 - rij(3)*r2kj )*dn1

        dl(1) = dn1*(rij(1)*r2kj-rkj(1)*scp1) - (rlk(1)*r2kj-rkj(1)*scp2)*th2
        dl(2) = dn1*(rij(2)*r2kj-rkj(2)*scp1) - (rlk(2)*r2kj-rkj(2)*scp2)*th2
        dl(3) = dn1*(rij(3)*r2kj-rkj(3)*scp1) - (rlk(3)*r2kj-rkj(3)*scp2)*th2

        ! calculate gradient
        geo%grd(:,i) = geo%grd(:,i) + f1*di(:)
        geo%grd(:,j) = geo%grd(:,j) + f1*dj(:)
        geo%grd(:,k) = geo%grd(:,k) + f1*dk(:)
        geo%grd(:,l) = geo%grd(:,l) + f1*dl(:)

        ! calculate hessian - diagonal i
        lh(1,1) = - (r2kj-rkj(1)**2)*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * (rij(1)*r2kj-rkj(1)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4 &
                  - 2.0 * (rij(1)*r2kj-rkj(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2
        lh(2,2) = - (r2kj-rkj(2)**2)*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * (rij(2)*r2kj-rkj(2)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4 &
                  - 2.0 * (rij(2)*r2kj-rkj(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2
        lh(3,3) = - (r2kj-rkj(3)**2)*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * (rij(3)*r2kj-rkj(3)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4 &
                  - 2.0 * (rij(3)*r2kj-rkj(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2
        lh(1,2) = + rkj(1)*rkj(2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rij(1)*r2kj-rkj(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2
        lh(1,3) = + rkj(1)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rij(1)*r2kj-rkj(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2
        lh(2,3) = + rkj(2)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rij(2)*r2kj-rkj(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2

        geo%hess(1,i,1,i) = geo%hess(1,i,1,i) + f1*lh(1,1) + f2*di(1)*di(1)
        geo%hess(2,i,2,i) = geo%hess(2,i,2,i) + f1*lh(2,2) + f2*di(2)*di(2)
        geo%hess(3,i,3,i) = geo%hess(3,i,3,i) + f1*lh(3,3) + f2*di(3)*di(3)
        geo%hess(1,i,2,i) = geo%hess(1,i,2,i) + f1*lh(1,2) + f2*di(1)*di(2)
        geo%hess(1,i,3,i) = geo%hess(1,i,3,i) + f1*lh(1,3) + f2*di(1)*di(3)
        geo%hess(2,i,3,i) = geo%hess(2,i,3,i) + f1*lh(2,3) + f2*di(2)*di(3)
        geo%hess(2,i,1,i) = geo%hess(2,i,1,i) + f1*lh(1,2) + f2*di(1)*di(2)
        geo%hess(3,i,1,i) = geo%hess(3,i,1,i) + f1*lh(1,3) + f2*di(1)*di(3)
        geo%hess(3,i,2,i) = geo%hess(3,i,2,i) + f1*lh(2,3) + f2*di(2)*di(3)

        ! calculate hessian - diagonal j
        lh(1,1) = - (-2.0d0*scp1 + r2kj + r2ij + 4.0*rij(1)*rkj(1) - (rkj(1)+rij(1))**2 )*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * ( (rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)**2 * (r2kj*scp3-scp1*scp2)*dn4 &
                  + 2.0 * ( (rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)    * (rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2)*dn5 &
                  - (r2lk-rlk(1)**2)*(r2kj*scp3-scp1*scp2)*dn3 &
                  + 3.0 * (rlk(1)*scp2 - rkj(1)*r2lk)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * ( (rkj(1)+rij(1))*scp1 - rij(1)*r2kj - rkj(1)*r2ij) &
                        * (-2.0d0*rkj(1)*scp3 + (rkj(1)+rij(1))*scp2 + rlk(1)*scp1 - rlk(1)*r2kj) * dn2 &
                  - 2.0 * ( rlk(1)*scp2-rkj(1)*r2lk)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj)*dn3 &
                  + (2.0 * scp3 - 2.0d0*scp2 + 4.0*rkj(1)*rlk(1) - 2.0d0*(rkj(1)+rij(1))*rlk(1)) * dn1
        lh(2,2) =  - (-2.0d0*scp1 + r2kj + r2ij + 4.0*rij(2)*rkj(2) - (rkj(2)+rij(2))**2 )*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * ( (rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)**2 * (r2kj*scp3-scp1*scp2)*dn4 &
                  + 2.0 * ( (rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)    * (rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2)*dn5 &
                  - (r2lk-rlk(2)**2)*(r2kj*scp3-scp1*scp2)*dn3 &
                  + 3.0 * (rlk(2)*scp2 - rkj(2)*r2lk)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * ( (rkj(2)+rij(2))*scp1 - rij(2)*r2kj - rkj(2)*r2ij) &
                        * (-2.0d0*rkj(2)*scp3 + (rkj(2)+rij(2))*scp2 + rlk(2)*scp1 - rlk(2)*r2kj) * dn2 &
                  - 2.0 * ( rlk(2)*scp2-rkj(2)*r2lk)*(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj)*dn3 &
                  + (2.0 * scp3 - 2.0d0*scp2 + 4.0*rkj(2)*rlk(2) - 2.0d0*(rkj(2)+rij(2))*rlk(2)) * dn1
        lh(3,3) =  - (-2.0d0*scp1 + r2kj + r2ij + 4.0*rij(3)*rkj(3) - (rkj(3)+rij(3))**2 )*(r2kj*scp3-scp1*scp2)*dn2 &
                  + 3.0 * ( (rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)**2 * (r2kj*scp3-scp1*scp2)*dn4 &
                  + 2.0 * ( (rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)    * (rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2)*dn5 &
                  - (r2lk-rlk(3)**2)*(r2kj*scp3-scp1*scp2)*dn3 &
                  + 3.0 * (rlk(3)*scp2 - rkj(3)*r2lk)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * ( (rkj(3)+rij(3))*scp1 - rij(3)*r2kj - rkj(3)*r2ij) &
                        * (-2.0d0*rkj(3)*scp3 + (rkj(3)+rij(3))*scp2 + rlk(3)*scp1 - rlk(3)*r2kj) * dn2 &
                  - 2.0 * ( rlk(3)*scp2-rkj(3)*r2lk)*(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj)*dn3 &
                  + (2.0 * scp3 - 2.0d0*scp2 + 4.0*rkj(3)*rlk(3) - 2.0d0*(rkj(3)+rij(3))*rlk(3)) * dn1

        lh(1,2) = - (2.0d0*rij(1)*rkj(2)-(rkj(1)+rij(1))*(rkj(2)+rij(2))+2.0d0*rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                       *((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + rlk(1)*rlk(2)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(1)*scp2-rkj(1)*r2lk)*(rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  -  ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn2 &
                  - (rlk(1)*scp2-rkj(1)*r2lk)*(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn3 &
                  - ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn2 &
                  - (rlk(2)*scp2-rkj(2)*r2lk)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn3 &
                  + (2.0d0*rkj(1)*rlk(2)-(rkj(1)+rij(1))*rlk(2)+2.0d0*rkj(2)*rlk(1)-(rkj(2)+rij(2))*rlk(1)) * dn1
        lh(1,3) = - (2.0d0*rij(1)*rkj(3)-(rkj(1)+rij(1))*(rkj(3)+rij(3))+2.0d0*rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                       *((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + rlk(1)*rlk(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(1)*scp2-rkj(1)*r2lk)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn2 &
                  - (rlk(1)*scp2-rkj(1)*r2lk)*(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn3 &
                  - ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn2 &
                  - (rlk(3)*scp2-rkj(3)*r2lk)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn3 &
                  + (2.0d0*rkj(1)*rlk(3)-(rkj(1)+rij(1))*rlk(3)+2.0d0*rkj(3)*rlk(1)-(rkj(3)+rij(3))*rlk(1)) * dn1
        lh(2,3) = - (2.0d0*rij(2)*rkj(3)-(rkj(2)+rij(2))*(rkj(3)+rij(3))+2.0d0*rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                       *((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(2)*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + rlk(2)*rlk(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(2)*scp2-rkj(2)*r2lk)*(rlk(3)*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - ((rkj(2)+rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn2 &
                  - (rlk(2)*scp2-rkj(2)*r2lk)*(-2.0d0*rkj(3)*scp3+(rkj(3)+rij(3))*scp2+rlk(3)*scp1-rlk(3)*r2kj) * dn3 &
                  - ((rkj(3)+rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn2 &
                  - (rlk(3)*scp2-rkj(3)*r2lk)*(-2.0d0*rkj(2)*scp3+(rkj(2)+rij(2))*scp2+rlk(2)*scp1-rlk(2)*r2kj) * dn3 &
                  + (2.0d0*rkj(2)*rlk(3)-(rkj(2)+rij(2))*rlk(3)+2.0d0*rkj(3)*rlk(2)-(rkj(3)+rij(3))*rlk(2)) * dn1

        geo%hess(1,j,1,j) = geo%hess(1,j,1,j) + f1*lh(1,1) + f2*dj(1)*dj(1)
        geo%hess(2,j,2,j) = geo%hess(2,j,2,j) + f1*lh(2,2) + f2*dj(2)*dj(2)
        geo%hess(3,j,3,j) = geo%hess(3,j,3,j) + f1*lh(3,3) + f2*dj(3)*dj(3)
        geo%hess(1,j,2,j) = geo%hess(1,j,2,j) + f1*lh(1,2) + f2*dj(1)*dj(2)
        geo%hess(1,j,3,j) = geo%hess(1,j,3,j) + f1*lh(1,3) + f2*dj(1)*dj(3)
        geo%hess(2,j,3,j) = geo%hess(2,j,3,j) + f1*lh(2,3) + f2*dj(2)*dj(3)
        geo%hess(2,j,1,j) = geo%hess(2,j,1,j) + f1*lh(1,2) + f2*dj(1)*dj(2)
        geo%hess(3,j,1,j) = geo%hess(3,j,1,j) + f1*lh(1,3) + f2*dj(1)*dj(3)
        geo%hess(3,j,2,j) = geo%hess(3,j,2,j) + f1*lh(2,3) + f2*dj(2)*dj(3)

        ! calculate hessian - diagonal k
        lh(1,1) = - (r2ij-rij(1)**2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(1)*r2ij-rij(1)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4  &
                  + 2.0d0*(rkj(1)*r2ij-rij(1)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2.0d0*scp2+r2lk+r2kj-(rlk(1)-rkj(1))**2-4.0*rkj(1)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0d0*(rkj(1)*r2ij-rij(1)*scp1)*(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn2 &
                  - 2.0d0*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn3 &
                  + (2.0d0*scp3+2.0d0*scp1-2.0d0*rij(1)*(rlk(1)-rkj(1))-4.0*rij(1)*rkj(1)) * dn1
        lh(2,2) =  - (r2ij-rij(2)**2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(2)*r2ij-rij(2)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4  &
                  + 2.0d0*(rkj(2)*r2ij-rij(2)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2.0d0*scp2+r2lk+r2kj-(rlk(2)-rkj(2))**2-4.0*rkj(2)*rlk(2))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0d0*(rkj(2)*r2ij-rij(2)*scp1)*(2.0d0*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1-rij(2)*r2kj) * dn2 &
                  - 2.0d0*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                       *(2.0d0*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1-rij(2)*r2kj) * dn3 &
                  + (2.0d0*scp3+2.0d0*scp1-2.0d0*rij(2)*(rlk(2)-rkj(2))-4.0*rij(2)*rkj(2)) * dn1
        lh(3,3) =   - (r2ij-rij(3)**2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(3)*r2ij-rij(3)*scp1)**2 * (r2kj*scp3-scp1*scp2) * dn4  &
                  + 2.0d0*(rkj(3)*r2ij-rij(3)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2.0d0*scp2+r2lk+r2kj-(rlk(3)-rkj(3))**2-4.0*rkj(3)*rlk(3))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0d0*(rkj(3)*r2ij-rij(3)*scp1)*(2.0d0*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1-rij(3)*r2kj) * dn2 &
                  - 2.0d0*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                       *(2.0d0*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1-rij(3)*r2kj) * dn3 &
                  + (2.0d0*scp3+2.0d0*scp1-2.0d0*rij(3)*(rlk(3)-rkj(3))-4.0*rij(3)*rkj(3)) * dn1

        lh(1,2) = + rij(1)*rij(2)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rkj(2)*r2ij-rij(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(2)-(rlk(1)-rkj(1))*(rlk(2)-rkj(2))-2*rkj(2)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                     *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (-rij(1)*(rlk(2)-rkj(2))-rij(2)*(rlk(1)-rkj(1))+2*(-rij(1))*rkj(2)+2*(-rij(2))*rkj(1))*dn1
        lh(1,3) = + rij(1)*rij(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rkj(3)*r2ij-rij(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(3)-(rlk(1)-rkj(1))*(rlk(3)-rkj(3))-2*rkj(3)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (-rij(1)*(rlk(3)-rkj(3))-rij(3)*(rlk(1)-rkj(1))+2*(-rij(1))*rkj(3)+2*(-rij(3))*rkj(1))*dn1
        lh(2,3) = + rij(2)*rij(3)*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rkj(3)*r2ij-rij(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(2)*rlk(3)-(rlk(2)-rkj(2))*(rlk(3)-rkj(3))-2*rkj(3)*rlk(2))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (-rij(2)*(rlk(3)-rkj(3))-rij(3)*(rlk(2)-rkj(2))+2*(-rij(2))*rkj(3)+2*(-rij(3))*rkj(2))*dn1

        geo%hess(1,k,1,k) = geo%hess(1,k,1,k) + f1*lh(1,1) + f2*dk(1)*dk(1)
        geo%hess(2,k,2,k) = geo%hess(2,k,2,k) + f1*lh(2,2) + f2*dk(2)*dk(2)
        geo%hess(3,k,3,k) = geo%hess(3,k,3,k) + f1*lh(3,3) + f2*dk(3)*dk(3)
        geo%hess(1,k,2,k) = geo%hess(1,k,2,k) + f1*lh(1,2) + f2*dk(1)*dk(2)
        geo%hess(1,k,3,k) = geo%hess(1,k,3,k) + f1*lh(1,3) + f2*dk(1)*dk(3)
        geo%hess(2,k,3,k) = geo%hess(2,k,3,k) + f1*lh(2,3) + f2*dk(2)*dk(3)
        geo%hess(2,k,1,k) = geo%hess(2,k,1,k) + f1*lh(1,2) + f2*dk(1)*dk(2)
        geo%hess(3,k,1,k) = geo%hess(3,k,1,k) + f1*lh(1,3) + f2*dk(1)*dk(3)
        geo%hess(3,k,2,k) = geo%hess(3,k,2,k) + f1*lh(2,3) + f2*dk(2)*dk(3)

        ! calculate hessian - diagonal l
        lh(1,1) = - (r2kj-rkj(1)**2) * (r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(1)*r2kj-rkj(1)*scp2)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * (r2kj*rij(1)-rkj(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(2,2) = - (r2kj-rkj(2)**2) * (r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(2)*r2kj-rkj(2)*scp2)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * (r2kj*rij(2)-rkj(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(3,3) = - (r2kj-rkj(3)**2) * (r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(3)*r2kj-rkj(3)*scp2)**2 * (r2kj*scp3-scp1*scp2) * dn6 &
                  - 2.0 * (r2kj*rij(3)-rkj(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3
        lh(1,2) =   rkj(1)*rkj(2)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (r2kj*rij(1)-rkj(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(1,3) = rkj(1)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (r2kj*rij(1)-rkj(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(2,3) =  rkj(2)*rkj(3)*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0 * (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (r2kj*rij(2)-rkj(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3

        geo%hess(1,l,1,l) = geo%hess(1,l,1,l) + f1*lh(1,1) + f2*dl(1)*dl(1)
        geo%hess(2,l,2,l) = geo%hess(2,l,2,l) + f1*lh(2,2) + f2*dl(2)*dl(2)
        geo%hess(3,l,3,l) = geo%hess(3,l,3,l) + f1*lh(3,3) + f2*dl(3)*dl(3)
        geo%hess(1,l,2,l) = geo%hess(1,l,2,l) + f1*lh(1,2) + f2*dl(1)*dl(2)
        geo%hess(1,l,3,l) = geo%hess(1,l,3,l) + f1*lh(1,3) + f2*dl(1)*dl(3)
        geo%hess(2,l,3,l) = geo%hess(2,l,3,l) + f1*lh(2,3) + f2*dl(2)*dl(3)
        geo%hess(2,l,1,l) = geo%hess(2,l,1,l) + f1*lh(1,2) + f2*dl(1)*dl(2)
        geo%hess(3,l,1,l) = geo%hess(3,l,1,l) + f1*lh(1,3) + f2*dl(1)*dl(3)
        geo%hess(3,l,2,l) = geo%hess(3,l,2,l) + f1*lh(2,3) + f2*dl(2)*dl(3)

        ! calculate hessian - off diagonal i,j
        lh(1,1) = - (scp1-r2kj-(-rkj(1)-rij(1))*rkj(1)-2*rij(1)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + ( rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(1)*rlk(1)-rkj(1)*(-rlk(1)) ) * dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(2,2) = - (scp1-r2kj-(-rkj(2)-rij(2))*rkj(2)-2*rij(2)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + ( rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(2)*rlk(2)-rkj(1)*(-rlk(1)) ) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(3,3) = - (scp1-r2kj-(-rkj(3)-rij(3))*rkj(3)-2*rij(3)*rkj(3))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + ( rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(3)*rlk(3)-rkj(1)*(-rlk(1)) ) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        lh(1,2) = - (-2*rij(1)*rkj(2)-rkj(1)*(-rkj(2)-rij(2)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + ( -rkj(1)*(-rlk(2))-2*rkj(2)*rlk(1) ) * dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(1,3) = - (-2*rij(1)*rkj(3)-rkj(1)*(-rkj(3)-rij(3)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + ( -rkj(1)*(-rlk(3))-2*rkj(3)*rlk(1) ) * dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3
        lh(2,3) = - (-2*rij(2)*rkj(3)-rkj(2)*(-rkj(3)-rij(3)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rij(2)*r2kj-rkj(2)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + ( -rkj(2)*(-rlk(3))-2*rkj(3)*rlk(2) ) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        lh(2,1) = - (-(-rkj(1)-rij(1))*rkj(2)-2*rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + (-2*rkj(1)*rlk(2)-rkj(2)*(-rlk(1))) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(3,1) = - (-(-rkj(1)-rij(1))*rkj(3)-2*rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (-2*rkj(1)*rlk(3)-rkj(3)*(-rlk(1))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(3,2) = - (-(-rkj(2)-rij(2))*rkj(3)-2*rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (-2*rkj(2)*rlk(3)-rkj(3)*(-rlk(2))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dj(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dj(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dj(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dj(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dj(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dj(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dj(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dj(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dj(2)

        geo%hess(1,i,1,j) = geo%hess(1,i,1,j) + hh(1,1)
        geo%hess(2,i,2,j) = geo%hess(2,i,2,j) + hh(2,2)
        geo%hess(3,i,3,j) = geo%hess(3,i,3,j) + hh(3,3)
        geo%hess(1,i,2,j) = geo%hess(1,i,2,j) + hh(1,2)
        geo%hess(1,i,3,j) = geo%hess(1,i,3,j) + hh(1,3)
        geo%hess(2,i,3,j) = geo%hess(2,i,3,j) + hh(2,3)
        geo%hess(2,i,1,j) = geo%hess(2,i,1,j) + hh(2,1)
        geo%hess(3,i,1,j) = geo%hess(3,i,1,j) + hh(3,1)
        geo%hess(3,i,2,j) = geo%hess(3,i,2,j) + hh(3,2)

        geo%hess(1,j,1,i) = geo%hess(1,j,1,i) + hh(1,1)
        geo%hess(2,j,2,i) = geo%hess(2,j,2,i) + hh(2,2)
        geo%hess(3,j,3,i) = geo%hess(3,j,3,i) + hh(3,3)
        geo%hess(1,j,2,i) = geo%hess(1,j,2,i) + hh(2,1)
        geo%hess(1,j,3,i) = geo%hess(1,j,3,i) + hh(3,1)
        geo%hess(2,j,3,i) = geo%hess(2,j,3,i) + hh(3,2)
        geo%hess(2,j,1,i) = geo%hess(2,j,1,i) + hh(1,2)
        geo%hess(3,j,1,i) = geo%hess(3,j,1,i) + hh(1,3)
        geo%hess(3,j,2,i) = geo%hess(3,j,2,i) + hh(2,3)

        ! calculate hessian - off diagonal i,k
        lh(1,1) =  - (rij(1)*rkj(1)-scp1)*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (rij(1)*r2kj-rkj(1)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                   - (rkj(1)*r2ij-rij(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                   + (-rkj(3)*rlk(3)-rkj(2)*rlk(2)+rkj(1)*rlk(1)-rkj(1)*(rlk(1)-rkj(1))-rkj(3)**2-rkj(2)**2-rkj(1)**2)*dn1 &
                   - (rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(2,2) = - (rij(2)*rkj(2)-scp1)*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (rij(2)*r2kj-rkj(2)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                   - (rkj(2)*r2ij-rij(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                   + (-rkj(3)*rlk(3)+rkj(2)*rlk(2)-rkj(2)*(rlk(2)-rkj(2))-rkj(1)*rlk(1)-rkj(3)**2-rkj(2)**2-rkj(1)**2)*dn1 &
                   - (rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(3,3) = - (rij(3)*rkj(3)-scp1)*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(3)*r2ij-rij(3)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (rij(3)*r2kj-rkj(3)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                   - (rkj(3)*r2ij-rij(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                   + (rkj(3)*rlk(3)-rkj(3)*(rlk(3)-rkj(3))-rkj(2)*rlk(2)-rkj(1)*rlk(1)-rkj(3)**2-rkj(2)**2-rkj(1)**2)*dn1 &
                   - (rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        lh(1,2) = - (2*rij(1)*rkj(2)-rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + (2*rkj(2)*rlk(1)-rkj(1)*(rlk(2)-rkj(2)))* dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(1,3) = - (2*rij(1)*rkj(3)-rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(rij(1)*r2kj-rkj(1)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2) * dn2 &
                  + (2*rkj(3)*rlk(1)-rkj(1)*(rlk(3)-rkj(3)))* dn1 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3
        lh(2,3) = - (2*rij(2)*rkj(3)-rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + (2*rkj(3)*rlk(2)-rkj(2)*(rlk(3)-rkj(3)))* dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        lh(2,1) = - (2*rij(2)*rkj(1)-rij(1)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2) * dn2 &
                  + (2*rkj(1)*rlk(2)-rkj(2)*(rlk(1)-rkj(1))) * dn1 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(3,1) = - (2*rij(3)*rkj(1)-rij(1)*rkj(3))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (2*rkj(1)*rlk(3)-rkj(3)*(rlk(1)-rkj(1))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(3,2) = - (2*rij(3)*rkj(2)-rij(2)*rkj(3))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2) * dn2 &
                  + (2*rkj(2)*rlk(3)-rkj(3)*(rlk(2)-rkj(2))) * dn1 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dk(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dk(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dk(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dk(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dk(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dk(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dk(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dk(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dk(2)

        geo%hess(1,i,1,k) = geo%hess(1,i,1,k) + hh(1,1)
        geo%hess(2,i,2,k) = geo%hess(2,i,2,k) + hh(2,2)
        geo%hess(3,i,3,k) = geo%hess(3,i,3,k) + hh(3,3)
        geo%hess(1,i,2,k) = geo%hess(1,i,2,k) + hh(1,2)
        geo%hess(1,i,3,k) = geo%hess(1,i,3,k) + hh(1,3)
        geo%hess(2,i,3,k) = geo%hess(2,i,3,k) + hh(2,3)
        geo%hess(2,i,1,k) = geo%hess(2,i,1,k) + hh(2,1)
        geo%hess(3,i,1,k) = geo%hess(3,i,1,k) + hh(3,1)
        geo%hess(3,i,2,k) = geo%hess(3,i,2,k) + hh(3,2)

        geo%hess(1,k,1,i) = geo%hess(1,k,1,i) + hh(1,1)
        geo%hess(2,k,2,i) = geo%hess(2,k,2,i) + hh(2,2)
        geo%hess(3,k,3,i) = geo%hess(3,k,3,i) + hh(3,3)
        geo%hess(1,k,2,i) = geo%hess(1,k,2,i) + hh(2,1)
        geo%hess(1,k,3,i) = geo%hess(1,k,3,i) + hh(3,1)
        geo%hess(2,k,3,i) = geo%hess(2,k,3,i) + hh(3,2)
        geo%hess(2,k,1,i) = geo%hess(2,k,1,i) + hh(1,2)
        geo%hess(3,k,1,i) = geo%hess(3,k,1,i) + hh(1,3)
        geo%hess(3,k,2,i) = geo%hess(3,k,2,i) + hh(2,3)

        ! calculate hessian - off diagonal i,l
        lh(1,1) = + (rij(1)*r2kj-rkj(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)**2 + rkj(2)**2) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(1)*r2kj-rkj(1)*scp2) * dn3
        lh(2,2) = + (rij(2)*r2kj-rkj(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)**2 + rkj(1)**2) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(3,3) = + (rij(3)*r2kj-rkj(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(2)**2 + rkj(1)**2) * dn1 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3

        lh(1,2) = + (rij(1)*r2kj-rkj(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(2) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(1,3) = + (rij(1)*r2kj-rkj(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(3) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3
        lh(2,3) =  + (rij(2)*r2kj-rkj(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(2)*rkj(3) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3

        lh(2,1) = + (rij(2)*r2kj-rkj(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(2) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(2)*r2kj-rkj(2)*scp2) * dn3
        lh(3,1) = + (rij(3)*r2kj-rkj(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(1)*rkj(3) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3
        lh(3,2) = + (rij(3)*r2kj-rkj(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - rkj(2)*rkj(3) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(rlk(3)*r2kj-rkj(3)*scp2) * dn3

        hh(1,1) = f1*lh(1,1) + f2*di(1)*dl(1)
        hh(2,2) = f1*lh(2,2) + f2*di(2)*dl(2)
        hh(3,3) = f1*lh(3,3) + f2*di(3)*dl(3)

        hh(1,2) = f1*lh(1,2) + f2*di(1)*dl(2)
        hh(1,3) = f1*lh(1,3) + f2*di(1)*dl(3)
        hh(2,3) = f1*lh(2,3) + f2*di(2)*dl(3)

        hh(2,1) = f1*lh(2,1) + f2*di(2)*dl(1)
        hh(3,1) = f1*lh(3,1) + f2*di(3)*dl(1)
        hh(3,2) = f1*lh(3,2) + f2*di(3)*dl(2)

        geo%hess(1,i,1,l) = geo%hess(1,i,1,l) + hh(1,1)
        geo%hess(2,i,2,l) = geo%hess(2,i,2,l) + hh(2,2)
        geo%hess(3,i,3,l) = geo%hess(3,i,3,l) + hh(3,3)
        geo%hess(1,i,2,l) = geo%hess(1,i,2,l) + hh(1,2)
        geo%hess(1,i,3,l) = geo%hess(1,i,3,l) + hh(1,3)
        geo%hess(2,i,3,l) = geo%hess(2,i,3,l) + hh(2,3)
        geo%hess(2,i,1,l) = geo%hess(2,i,1,l) + hh(2,1)
        geo%hess(3,i,1,l) = geo%hess(3,i,1,l) + hh(3,1)
        geo%hess(3,i,2,l) = geo%hess(3,i,2,l) + hh(3,2)

        geo%hess(1,l,1,i) = geo%hess(1,l,1,i) + hh(1,1)
        geo%hess(2,l,2,i) = geo%hess(2,l,2,i) + hh(2,2)
        geo%hess(3,l,3,i) = geo%hess(3,l,3,i) + hh(3,3)
        geo%hess(1,l,2,i) = geo%hess(1,l,2,i) + hh(2,1)
        geo%hess(1,l,3,i) = geo%hess(1,l,3,i) + hh(3,1)
        geo%hess(2,l,3,i) = geo%hess(2,l,3,i) + hh(3,2)
        geo%hess(2,l,1,i) = geo%hess(2,l,1,i) + hh(1,2)
        geo%hess(3,l,1,i) = geo%hess(3,l,1,i) + hh(1,3)
        geo%hess(3,l,2,i) = geo%hess(3,l,2,i) + hh(2,3)

        ! calculate hessian - off diagonal j,k
        lh(1,1) = - (scp1-r2ij-2.0d0*rij(1)*rkj(1)+rij(1)*(rkj(1)+rij(1)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3.0*(rkj(1)*r2ij-rij(1)*scp1)*((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(1)*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2lk+2.0d0*rkj(1)*rlk(1)+rlk(1)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3.0*(rlk(1)*scp2-rkj(1)*r2lk)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - ((rkj(1)+rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn2 &
                  - (rlk(1)*scp2-rkj(1)*r2lk)*(2.0d0*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1-rij(1)*r2kj) * dn3 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-2.0d0*rkj(1)*scp3+(rkj(1)+rij(1))*scp2+rlk(1)*scp1-rlk(1)*r2kj) * dn3 &
                  + (-2.0d0*scp3+rkj(3)*rlk(3)+rkj(2)*rlk(2)+rkj(1)*rlk(1) + (rkj(1)+rij(1))*(rlk(1)-rkj(1)) &
                     -2.0d0*rkj(1)*rlk(1)+rij(1)*rlk(1)+rkj(3)**2 &
                     -rij(3)*rkj(3) + rkj(2)**2 - rij(2)*rkj(2) + rkj(1)**2 + 2.0d0*rij(1)*rkj(1) - rij(1)*rkj(1)) * dn1
        lh(2,2) =  - (scp1-r2ij-2*rij(2)*rkj(2)-rij(2)*(-rkj(2) - rij(2)))*(r2kj*scp3-scp1*scp2) * dn2 &
                   + 3*(rkj(2)*r2ij-rij(2)*scp1)*(-(-rkj(2) - rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                   + (-(-rkj(2) - rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                     *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                   + (rkj(2)*r2ij-rij(2)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                   - (-scp2-r2lk+2*rkj(2)*rlk(2)-(-rlk(2))*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                   + 3*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                   - (-(-rkj(2) - rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                     *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                   - (-(-rlk(2))*scp2-rkj(2)*r2lk)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                   - (rkj(2)*r2ij-rij(2)*scp1)*(-2*rkj(2)*scp3-(-rkj(2) - rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                   - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                     *(-2*rkj(2)*scp3-(-rkj(2) - rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                   + (-2*scp3+rkj(3)*rlk(3)+rkj(2)*rlk(2)-(-rkj(2) - rij(2))*(rlk(2)-rkj(2)) &
                      +2*rkj(2)*(-rlk(2))-rij(2)*(-rlk(2))+rkj(1)*rlk(1)+rkj(3)**2-rij(3)*rkj(3) &
                      +rkj(2)**2-2*(-rij(2))*rkj(2)-rij(2)*rkj(2)+rkj(1)**2-rij(1)*rkj(1)) * dn1
        lh(3,3) = - (scp1-r2ij-2*rij(3)*rkj(3)-rij(3)*(-rkj(3)-rij(3)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2lk+2*rkj(3)*rlk(3)-(-rlk(3))*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(-rlk(3))*scp2-rkj(3)*r2lk)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-2*scp3+rkj(3)*rlk(3)-(-rkj(3)-rij(3))*(rlk(3)-rkj(3))+2*rkj(3)*(-rlk(3))-rij(3)*(-rlk(3)) &
                     +rkj(2)*rlk(2)+rkj(1)*rlk(1)+rkj(3)**2-2*(-rij(3))*rkj(3)-rij(3)*rkj(3)+rkj(2)**2-rij(2)*rkj(2) &
                     +rkj(1)**2-rij(1)*rkj(1)) * dn1

        lh(1,2) = - (-2*rij(1)*rkj(2)-rij(2)*(-rkj(1)-rij(1)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(2)-(-rlk(1))*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(-rlk(1))*scp2-rkj(1)*r2lk)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*(rlk(2)-rkj(2))+2*rkj(2)*(-rlk(1))-rij(2)*(-rlk(1))-2*(-rij(2))*rkj(1)) * dn1
        lh(1,3) = - (-2*rij(1)*rkj(3)-rij(3)*(-rkj(1)-rij(1)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(3)-(-rlk(1))*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                     *(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(-rlk(1))*scp2-rkj(1)*r2lk)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*(rlk(3)-rkj(3))+2*rkj(3)*(-rlk(1))-rij(3)*(-rlk(1))-2*(-rij(3))*rkj(1)) * dn1
        lh(2,3) = - (-2*rij(2)*rkj(3)-rij(3)*(-rkj(2)-rij(2)))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(3)*r2ij-rij(3)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                     *(r2kj*scp3-scp1*scp2) * dn4 &
                  + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (rkj(3)*r2ij-rij(3)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(3)-(-rlk(2))*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3) &
                     *r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn2 &
                  - (-(-rlk(2))*scp2-rkj(2)*r2lk)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) &
                    *(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-(-rkj(2)-rij(2))*(rlk(3)-rkj(3))+2*rkj(3)*(-rlk(2))-rij(3)*(-rlk(2))-2*(-rij(3))*rkj(2)) * dn1

        lh(2,1) = - (-rij(1)*(-rkj(2)-rij(2))-2*rij(2)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(1)-(rlk(1)-rkj(1))*(-rlk(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(-rlk(2))*scp2-rkj(2)*r2lk)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rkj(1)*(-rlk(2))-rij(1)*(-rlk(2))-(-rkj(2)-rij(2))*(rlk(1)-rkj(1))-2*(-rij(1))*rkj(2)) * dn1
        lh(3,1) = - (-rij(1)*(-rkj(3)-rij(3))-2*rij(3)*rkj(1))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(1)*r2ij-rij(1)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(1)*r2ij-rij(1)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(1)-(rlk(1)-rkj(1))*(-rlk(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) &
                    *(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn2 &
                  - (-(-rlk(3))*scp2-rkj(3)*r2lk)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rkj(1)*(-rlk(3))-rij(1)*(-rlk(3))-(-rkj(3)-rij(3))*(rlk(1)-rkj(1))-2*(-rij(1))*rkj(3)) * dn1
        lh(3,2) =  - (-rij(2)*(-rkj(3)-rij(3))-2*rij(3)*rkj(2))*(r2kj*scp3-scp1*scp2) * dn2 &
                  + 3*(rkj(2)*r2ij-rij(2)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(r2kj*scp3-scp1*scp2) * dn4 &
                  + (rkj(2)*r2ij-rij(2)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn5 &
                  + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(2)-(rlk(2)-rkj(2))*(-rlk(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn2 &
                  - (-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) &
                    *(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) &
                    *(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn2 &
                  - (-(-rlk(3))*scp2-rkj(3)*r2lk)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (2*rkj(2)*(-rlk(3))-rij(2)*(-rlk(3))-(-rkj(3)-rij(3))*(rlk(2)-rkj(2))-2*(-rij(2))*rkj(3)) * dn1

        hh(1,1) = f1*lh(1,1) + f2*dj(1)*dk(1)
        hh(2,2) = f1*lh(2,2) + f2*dj(2)*dk(2)
        hh(3,3) = f1*lh(3,3) + f2*dj(3)*dk(3)

        hh(1,2) = f1*lh(1,2) + f2*dj(1)*dk(2)
        hh(1,3) = f1*lh(1,3) + f2*dj(1)*dk(3)
        hh(2,3) = f1*lh(2,3) + f2*dj(2)*dk(3)

        hh(2,1) = f1*lh(2,1) + f2*dj(2)*dk(1)
        hh(3,1) = f1*lh(3,1) + f2*dj(3)*dk(1)
        hh(3,2) = f1*lh(3,2) + f2*dj(3)*dk(2)

        geo%hess(1,j,1,k) = geo%hess(1,j,1,k) + hh(1,1)
        geo%hess(2,j,2,k) = geo%hess(2,j,2,k) + hh(2,2)
        geo%hess(3,j,3,k) = geo%hess(3,j,3,k) + hh(3,3)
        geo%hess(1,j,2,k) = geo%hess(1,j,2,k) + hh(1,2)
        geo%hess(1,j,3,k) = geo%hess(1,j,3,k) + hh(1,3)
        geo%hess(2,j,3,k) = geo%hess(2,j,3,k) + hh(2,3)
        geo%hess(2,j,1,k) = geo%hess(2,j,1,k) + hh(2,1)
        geo%hess(3,j,1,k) = geo%hess(3,j,1,k) + hh(3,1)
        geo%hess(3,j,2,k) = geo%hess(3,j,2,k) + hh(3,2)

        geo%hess(1,k,1,j) = geo%hess(1,k,1,j) + hh(1,1)
        geo%hess(2,k,2,j) = geo%hess(2,k,2,j) + hh(2,2)
        geo%hess(3,k,3,j) = geo%hess(3,k,3,j) + hh(3,3)
        geo%hess(1,k,2,j) = geo%hess(1,k,2,j) + hh(2,1)
        geo%hess(1,k,3,j) = geo%hess(1,k,3,j) + hh(3,1)
        geo%hess(2,k,3,j) = geo%hess(2,k,3,j) + hh(3,2)
        geo%hess(2,k,1,j) = geo%hess(2,k,1,j) + hh(1,2)
        geo%hess(3,k,1,j) = geo%hess(3,k,1,j) + hh(1,3)
        geo%hess(3,k,2,j) = geo%hess(3,k,2,j) + hh(2,3)

        ! calculate hessian - off diagonal j,l
        lh(1,1) = + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (scp2-2*rkj(1)*rlk(1)-rkj(1)*(-rlk(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-rkj(3)**2+rij(3)*rkj(3)-rkj(2)**2+rij(2)*rkj(2)-rkj(1)**2-(-rkj(1)-rij(1))*rkj(1)-rij(1)*rkj(1)) * dn1 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(2,2) = + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (scp2-2*rkj(2)*rlk(2)-rkj(2)*(-rlk(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-rkj(3)**2+rij(3)*rkj(3)-rkj(2)**2-(-rkj(2)-rij(2))*rkj(2)-rij(2)*rkj(2)-rkj(1)**2+rij(1)*rkj(1)) * dn1 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(3,3) = + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (scp2-2*rkj(3)*rlk(3)-rkj(3)*(-rlk(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-rkj(3)**2-(-rkj(3)-rij(3))*rkj(3)-rij(3)*rkj(3)-rkj(2)**2+rij(2)*rkj(2)-rkj(1)**2+rij(1)*rkj(1)) * dn1 &
                  - (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        lh(1,2) = + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(2)-rkj(2)*(-rlk(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*rkj(2)-2*rij(2)*rkj(1)) * dn1 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(1,3) = + (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(1)*rlk(3)-rkj(3)*(-rlk(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(1))*scp2-rkj(1)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-2*rkj(1)*scp3-(-rkj(1)-rij(1))*scp2-(-rlk(1))*scp1+(-rlk(1))*r2kj) * dn3 &
                  + (-(-rkj(1)-rij(1))*rkj(3)-2*rij(3)*rkj(1)) * dn1 &
                  - (-(-rkj(1)-rij(1))*scp1-rij(1)*r2kj-rkj(1)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(1))*scp2-rkj(1)*r2lk) * dn3
        lh(2,3) = + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-2*rkj(2)*rlk(3)-rkj(3)*(-rlk(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-(-rkj(2)-rij(2))*rkj(3)-2*rij(3)*rkj(2)) * dn1 &
                  - (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3

        lh(2,1) = + (-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-rkj(1)*(-rlk(2))-2*rkj(2)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(2))*scp2-rkj(2)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-2*rkj(2)*scp3-(-rkj(2)-rij(2))*scp2-(-rlk(2))*scp1+(-rlk(2))*r2kj) * dn3 &
                  + (-2*rij(1)*rkj(2)-rkj(1)*(-rkj(2)-rij(2))) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(2)-rij(2))*scp1-rij(2)*r2kj-rkj(2)*r2ij) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(2))*scp2-rkj(2)*r2lk) * dn3
        lh(3,1) = + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-rkj(1)*(-rlk(3))-2*rkj(3)*rlk(1))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-2*rij(1)*rkj(3)-rkj(1)*(-rkj(3)-rij(3))) * dn1 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3
        lh(3,2) = + (-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-rkj(2)*(-rlk(3))-2*rkj(3)*rlk(2))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(-rlk(3))*scp2-rkj(3)*r2lk)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(-2*rkj(3)*scp3-(-rkj(3)-rij(3))*scp2-(-rlk(3))*scp1+(-rlk(3))*r2kj) * dn3 &
                  + (-2*rij(2)*rkj(3)-rkj(2)*(-rkj(3)-rij(3))) * dn1 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rkj(3)-rij(3))*scp1-rij(3)*r2kj-rkj(3)*r2ij) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(-rlk(3))*scp2-rkj(3)*r2lk) * dn3

        hh(1,1) = f1*lh(1,1) + f2*dj(1)*dl(1)
        hh(2,2) = f1*lh(2,2) + f2*dj(2)*dl(2)
        hh(3,3) = f1*lh(3,3) + f2*dj(3)*dl(3)

        hh(1,2) = f1*lh(1,2) + f2*dj(1)*dl(2)
        hh(1,3) = f1*lh(1,3) + f2*dj(1)*dl(3)
        hh(2,3) = f1*lh(2,3) + f2*dj(2)*dl(3)

        hh(2,1) = f1*lh(2,1) + f2*dj(2)*dl(1)
        hh(3,1) = f1*lh(3,1) + f2*dj(3)*dl(1)
        hh(3,2) = f1*lh(3,2) + f2*dj(3)*dl(2)

        geo%hess(1,j,1,l) = geo%hess(1,j,1,l) + hh(1,1)
        geo%hess(2,j,2,l) = geo%hess(2,j,2,l) + hh(2,2)
        geo%hess(3,j,3,l) = geo%hess(3,j,3,l) + hh(3,3)
        geo%hess(1,j,2,l) = geo%hess(1,j,2,l) + hh(1,2)
        geo%hess(1,j,3,l) = geo%hess(1,j,3,l) + hh(1,3)
        geo%hess(2,j,3,l) = geo%hess(2,j,3,l) + hh(2,3)
        geo%hess(2,j,1,l) = geo%hess(2,j,1,l) + hh(2,1)
        geo%hess(3,j,1,l) = geo%hess(3,j,1,l) + hh(3,1)
        geo%hess(3,j,2,l) = geo%hess(3,j,2,l) + hh(3,2)

        geo%hess(1,l,1,j) = geo%hess(1,l,1,j) + hh(1,1)
        geo%hess(2,l,2,j) = geo%hess(2,l,2,j) + hh(2,2)
        geo%hess(3,l,3,j) = geo%hess(3,l,3,j) + hh(3,3)
        geo%hess(1,l,2,j) = geo%hess(1,l,2,j) + hh(2,1)
        geo%hess(1,l,3,j) = geo%hess(1,l,3,j) + hh(3,1)
        geo%hess(2,l,3,j) = geo%hess(2,l,3,j) + hh(3,2)
        geo%hess(2,l,1,j) = geo%hess(2,l,1,j) + hh(1,2)
        geo%hess(3,l,1,j) = geo%hess(3,l,1,j) + hh(1,3)
        geo%hess(3,l,2,j) = geo%hess(3,l,2,j) + hh(2,3)

        ! calculate hessian - off diagonal k,l
        lh(1,1) = + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2kj+2*rkj(1)*rlk(1)-rkj(1)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (-rij(3)*rkj(3)-rij(2)*rkj(2)) * dn1 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(2,2) = + (rkj(2)*r2ij-rij(2)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2kj+2*rkj(2)*rlk(2)-rkj(2)*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (-rij(3)*rkj(3)-rij(1)*rkj(1)) * dn1 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(3,3) = + (rkj(3)*r2ij-rij(3)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (-scp2-r2kj+2*rkj(3)*rlk(3)-rkj(3)*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) * dn3 &
                  + (-rij(2)*rkj(2)-rij(1)*rkj(1)) * dn1 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        lh(1,2) = + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(2)-rkj(2)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rij(2)*rkj(1)-rij(1)*rkj(2)) * dn1 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(1,3) = + (rkj(1)*r2ij-rij(1)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(1)*rlk(3)-rkj(3)*(rlk(1)-rkj(1)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(2*rkj(1)*scp3-rij(1)*scp2-(rlk(1)-rkj(1))*scp1+(-rij(1))*r2kj) * dn3 &
                  + (2*rij(3)*rkj(1)-rij(1)*rkj(3)) * dn1 &
                  - (rkj(1)*r2ij-rij(1)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(1)-rkj(1))*scp2+rkj(1)*r2lk-rlk(1)*r2kj) * dn3
        lh(2,3) = + (rkj(2)*r2ij-rij(2)*scp1)*(rlk(3)*r2kj-rkj(3)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(3)-rkj(3)*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(3)*r2kj-rkj(3)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(3)*r2kj-rkj(3)*scp2)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) * dn3 &
                  + (2*rij(3)*rkj(2)-rij(2)*rkj(3)) * dn1 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rij(3)*r2kj-rkj(3)*scp1) * dn2 &
                  - (rij(3)*r2kj-rkj(3)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3

        lh(2,1) = + (rkj(2)*r2ij-rij(2)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(2)*rlk(1)-rkj(1)*(rlk(2)-rkj(2)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(2*rkj(2)*scp3-rij(2)*scp2-(rlk(2)-rkj(2))*scp1+(-rij(2))*r2kj) *dn3 &
                  + (2*rij(1)*rkj(2)-rij(2)*rkj(1)) * dn1 &
                  - (rkj(2)*r2ij-rij(2)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(2)-rkj(2))*scp2+rkj(2)*r2lk-rlk(2)*r2kj) * dn3
        lh(3,1) = + (rkj(3)*r2ij-rij(3)*scp1)*(rlk(1)*r2kj-rkj(1)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(1)-rkj(1)*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(1)*r2kj-rkj(1)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(1)*r2kj-rkj(1)*scp2)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) *dn3 &
                  + (2*rij(1)*rkj(3)-rij(3)*rkj(1)) * dn1 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rij(1)*r2kj-rkj(1)*scp1) * dn2 &
                  - (rij(1)*r2kj-rkj(1)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3
        lh(3,2) = + (rkj(3)*r2ij-rij(3)*scp1)*(rlk(2)*r2kj-rkj(2)*scp2)*(r2kj*scp3-scp1*scp2) * dn5 &
                  - (2*rkj(3)*rlk(2)-rkj(2)*(rlk(3)-rkj(3)))*(r2kj*scp3-scp1*scp2) * dn3 &
                  + 3*(rlk(2)*r2kj-rkj(2)*scp2)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj)*(r2kj*scp3-scp1*scp2) * dn6 &
                  - (rlk(2)*r2kj-rkj(2)*scp2)*(2*rkj(3)*scp3-rij(3)*scp2-(rlk(3)-rkj(3))*scp1+(-rij(3))*r2kj) *dn3 &
                  + (2*rij(2)*rkj(3)-rij(3)*rkj(2)) * dn1 &
                  - (rkj(3)*r2ij-rij(3)*scp1)*(rij(2)*r2kj-rkj(2)*scp1) * dn2 &
                  - (rij(2)*r2kj-rkj(2)*scp1)*(-(rlk(3)-rkj(3))*scp2+rkj(3)*r2lk-rlk(3)*r2kj) * dn3

        hh(1,1) = f1*lh(1,1) + f2*dk(1)*dl(1)
        hh(2,2) = f1*lh(2,2) + f2*dk(2)*dl(2)
        hh(3,3) = f1*lh(3,3) + f2*dk(3)*dl(3)

        hh(1,2) = f1*lh(1,2) + f2*dk(1)*dl(2)
        hh(1,3) = f1*lh(1,3) + f2*dk(1)*dl(3)
        hh(2,3) = f1*lh(2,3) + f2*dk(2)*dl(3)

        hh(2,1) = f1*lh(2,1) + f2*dk(2)*dl(1)
        hh(3,1) = f1*lh(3,1) + f2*dk(3)*dl(1)
        hh(3,2) = f1*lh(3,2) + f2*dk(3)*dl(2)

        geo%hess(1,k,1,l) = geo%hess(1,k,1,l) + hh(1,1)
        geo%hess(2,k,2,l) = geo%hess(2,k,2,l) + hh(2,2)
        geo%hess(3,k,3,l) = geo%hess(3,k,3,l) + hh(3,3)
        geo%hess(1,k,2,l) = geo%hess(1,k,2,l) + hh(1,2)
        geo%hess(1,k,3,l) = geo%hess(1,k,3,l) + hh(1,3)
        geo%hess(2,k,3,l) = geo%hess(2,k,3,l) + hh(2,3)
        geo%hess(2,k,1,l) = geo%hess(2,k,1,l) + hh(2,1)
        geo%hess(3,k,1,l) = geo%hess(3,k,1,l) + hh(3,1)
        geo%hess(3,k,2,l) = geo%hess(3,k,2,l) + hh(3,2)

        geo%hess(1,l,1,k) = geo%hess(1,l,1,k) + hh(1,1)
        geo%hess(2,l,2,k) = geo%hess(2,l,2,k) + hh(2,2)
        geo%hess(3,l,3,k) = geo%hess(3,l,3,k) + hh(3,3)
        geo%hess(1,l,2,k) = geo%hess(1,l,2,k) + hh(2,1)
        geo%hess(1,l,3,k) = geo%hess(1,l,3,k) + hh(3,1)
        geo%hess(2,l,3,k) = geo%hess(2,l,3,k) + hh(3,2)
        geo%hess(2,l,1,k) = geo%hess(2,l,1,k) + hh(1,2)
        geo%hess(3,l,1,k) = geo%hess(3,l,1,k) + hh(1,3)
        geo%hess(3,l,2,k) = geo%hess(3,l,2,k) + hh(2,3)

    end do

end subroutine ffdev_hessian_impropers

!===============================================================================
! subroutine ffdev_hessian_nb_lj
!===============================================================================

subroutine ffdev_hessian_nb_lj(top,geo,fac)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    real(DEVDP)     :: fac
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,rij(3)
    real(DEVDP)     :: r2a,ra,r6a,Vela,V_aa,V_ba,dva
    real(DEVDP)     :: dn1,dn2,dn3,dn4,dn5,dn6
    real(DEVDP)     :: h11,h22,h33,h12,h13,h23,h21,h31,h32,hxx,hyy
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    ! FIXME

!    do ip=1,top%nb_size
!        i = top%nb_list(ip)%ai
!        j = top%nb_list(ip)%aj
!        nbt = top%nb_list(ip)%nbt
!        aLJa  = top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**12
!        bLJa  = 2.0d0*top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**6
!        crgij = top%atoms(i)%charge*top%atoms(j)%charge*332.05221729d0
!
!        ! calculate dx, r and r2
!        rij(:) = geo%crd(:,i) - geo%crd(:,j)
!
!        r2a = rij(1)**2 + rij(2)**2 + rij(3)**2
!        r2a = 1.0d0/r2a
!        ra  = sqrt(r2a)
!        r6a = r2a*r2a*r2a
!
!        Vela = crgij*ra
!        V_aa = aLJa*r6a*r6a
!        V_ba = bLJa*r6a
!
!        ! calculate energy
!        if( top%nb_list(ip)%dt .eq. 0 ) then
!            inv_scee = 1.0d0
!            inv_scnb = 1.0d0
!            geo%ele_ene  = geo%ele_ene + Vela
!            geo%nb_ene  = geo%nb_ene + V_aa - V_ba
!            dva   = r2a*(Vela + 12.0d0*V_aa - 6.0d0*V_ba)
!        else
!            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
!            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
!            geo%ele14_ene  = geo%ele14_ene + inv_scee*Vela
!            geo%nb14_ene  = geo%nb14_ene + inv_scnb*(V_aa - V_ba)
!            dva   = r2a*(inv_scee*Vela + inv_scnb*(12.0d0*V_aa - 6.0d0*V_ba))
!        end if
!
!        ! calculate gradient
!        geo%grd(:,i) = geo%grd(:,i) - rij(:)*dva
!        geo%grd(:,j) = geo%grd(:,j) + rij(:)*dva
!
!        dn1 = ra*r2a
!        dn2 = ra*r2a*r2a
!        dn3 = r2a**7
!        dn4 = dn3*r2a
!        dn5 = r2a**4
!        dn6 = dn5*r2a
!
!        ! calculate hessian
!        ! ----------------------------------------
!        hxx = 3.0d0*inv_scee*crgij*dn2 + 168.0d0*inv_scnb*aLJa*dn4 - 48.0d0*inv_scnb*bLJa*dn6
!        hyy = crgij*inv_scee*dn1 + 12.0d0*inv_scnb*aLJa*dn3 - 6.0d0*inv_scnb*bLJa*dn5
!
!        ! calculate hessian - diagonals
!        h11 = hxx*rij(1)**2 - hyy
!        h22 = hxx*rij(2)**2 - hyy
!        h33 = hxx*rij(3)**2 - hyy
!
!        geo%hess(1,i,1,i) = geo%hess(1,i,1,i) + fac*h11
!        geo%hess(2,i,2,i) = geo%hess(2,i,2,i) + fac*h22
!        geo%hess(3,i,3,i) = geo%hess(3,i,3,i) + fac*h33
!
!        geo%hess(1,j,1,j) = geo%hess(1,j,1,j) + fac*h11
!        geo%hess(2,j,2,j) = geo%hess(2,j,2,j) + fac*h22
!        geo%hess(3,j,3,j) = geo%hess(3,j,3,j) + fac*h33
!
!        ! calculate hessian - off-diagonals - common
!        h12 = hxx*rij(1)*rij(2)
!        h13 = hxx*rij(1)*rij(3)
!        h23 = hxx*rij(2)*rij(3)
!
!        geo%hess(1,i,2,i) = geo%hess(1,i,2,i) + fac*h12
!        geo%hess(1,i,3,i) = geo%hess(1,i,3,i) + fac*h13
!        geo%hess(2,i,3,i) = geo%hess(2,i,3,i) + fac*h23
!
!        geo%hess(2,i,1,i) = geo%hess(2,i,1,i) + fac*h12
!        geo%hess(3,i,1,i) = geo%hess(3,i,1,i) + fac*h13
!        geo%hess(3,i,2,i) = geo%hess(3,i,2,i) + fac*h23
!
!        geo%hess(1,j,2,j) = geo%hess(1,j,2,j) + fac*h12
!        geo%hess(1,j,3,j) = geo%hess(1,j,3,j) + fac*h13
!        geo%hess(2,j,3,j) = geo%hess(2,j,3,j) + fac*h23
!
!        geo%hess(2,j,1,j) = geo%hess(2,j,1,j) + fac*h12
!        geo%hess(3,j,1,j) = geo%hess(3,j,1,j) + fac*h13
!        geo%hess(3,j,2,j) = geo%hess(3,j,2,j) + fac*h23
!
!       ! calculate hessian - off-diagonals - cross
!        h11 = hyy - hxx*rij(1)*rij(1)
!        h12 =     - hxx*rij(1)*rij(2)
!        h13 =     - hxx*rij(1)*rij(3)
!
!        h21 =     - hxx*rij(2)*rij(1)
!        h22 = hyy - hxx*rij(2)*rij(2)
!        h23 =     - hxx*rij(2)*rij(3)
!
!        h31 =     - hxx*rij(3)*rij(1)
!        h32 =     - hxx*rij(3)*rij(2)
!        h33 = hyy - hxx*rij(3)*rij(3)
!
!        geo%hess(1,i,1,j) = geo%hess(1,i,1,j) + fac*h11
!        geo%hess(1,i,2,j) = geo%hess(1,i,2,j) + fac*h12
!        geo%hess(1,i,3,j) = geo%hess(1,i,3,j) + fac*h13
!
!        geo%hess(2,i,1,j) = geo%hess(2,i,1,j) + fac*h21
!        geo%hess(2,i,2,j) = geo%hess(2,i,2,j) + fac*h22
!        geo%hess(2,i,3,j) = geo%hess(2,i,3,j) + fac*h23
!
!        geo%hess(3,i,1,j) = geo%hess(3,i,1,j) + fac*h31
!        geo%hess(3,i,2,j) = geo%hess(3,i,2,j) + fac*h32
!        geo%hess(3,i,3,j) = geo%hess(3,i,3,j) + fac*h33
!
!        geo%hess(1,j,1,i) = geo%hess(1,j,1,i) + fac*h11
!        geo%hess(1,j,2,i) = geo%hess(1,j,2,i) + fac*h21
!        geo%hess(1,j,3,i) = geo%hess(1,j,3,i) + fac*h31
!
!        geo%hess(2,j,1,i) = geo%hess(2,j,1,i) + fac*h12
!        geo%hess(2,j,2,i) = geo%hess(2,j,2,i) + fac*h22
!        geo%hess(2,j,3,i) = geo%hess(2,j,3,i) + fac*h32
!
!        geo%hess(3,j,1,i) = geo%hess(3,j,1,i) + fac*h13
!        geo%hess(3,j,2,i) = geo%hess(3,j,2,i) + fac*h23
!        geo%hess(3,j,3,i) = geo%hess(3,j,3,i) + fac*h33
!    end do

end subroutine ffdev_hessian_nb_lj

! ------------------------------------------------------------------------------

end module ffdev_hessian

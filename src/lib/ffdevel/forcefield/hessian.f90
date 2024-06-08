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
    use ffdev_nbmode_LJ

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

    ! reset energy
    geo%bond_ene = 0.0d0
    geo%angle_ene = 0.0d0
    geo%dih_ene = 0.0d0
    geo%impropr_ene = 0.0d0

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%pen_ene = 0.0d0
    geo%ind_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    geo%total_ene = 0.0d0
    geo%rst_energy = 0.0d0

    ! reset gradient
    geo%grd(:,:) = 0.0d0

    ! reset hessian
    geo%hess(:,:,:,:) = 0.0d0

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_hessian_bonds(top,geo)
        call ffdev_hessian_angles(top,geo)
        call ffdev_hessian_dihedrals(top,geo)
        call ffdev_hessian_impropers(top,geo)
    end if

    if( calcnb ) then
        if( top%nb_params_update ) then
            call ffdev_topology_update_nb_params(top)
            call ffdev_topology_update_nb_pairs(top)
        end if

        ! non-bonded terms
        select case(nb_mode)
            case(NB_VDW_LJ)
                call ffdev_hessian_nb_lj(top,geo)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported vdW mode in ffdev_hessian_all!')
        end select
    end if

    geo%total_ene = geo%bond_ene + geo%angle_ene + geo%dih_ene + geo%impropr_ene &
                  + geo%ele14_ene + geo%rep14_ene + geo%dis14_ene  &
                  + geo%ele_ene + geo%pen_ene + geo%ind_ene + geo%rep_ene + geo%dis_ene

end subroutine ffdev_hessian_all

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
! subroutine vector_product
!===============================================================================

subroutine vector_product(a,b,v)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    real(DEVDP)     :: a(3),b(3)
    real(DEVDP)     :: v(3)
    ! -----------------------------------------------------------------------------

    v(1) = a(2)*b(3) - a(3)*b(2)
    v(2) = a(3)*b(1) - a(1)*b(3)
    v(3) = a(1)*b(2) - a(2)*b(1)

end subroutine vector_product

!===============================================================================
! subroutine tensor_product
!===============================================================================

subroutine tensor_product(a,b,v)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    real(DEVDP)     :: a(3),b(3)
    real(DEVDP)     :: v(3,3)
    ! -----------------------------------------------------------------------------

    v(1,1) = a(1)*b(1)
    v(1,2) = a(1)*b(2)
    v(1,3) = a(1)*b(3)

    v(2,1) = a(2)*b(1)
    v(2,2) = a(2)*b(2)
    v(2,3) = a(2)*b(3)

    v(3,1) = a(3)*b(1)
    v(3,2) = a(3)*b(2)
    v(3,3) = a(3)*b(3)

end subroutine tensor_product

!===============================================================================
! subroutine tensor_transpose
!===============================================================================

subroutine tensor_transpose(a,b)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    real(DEVDP)     :: a(3,3)
    real(DEVDP)     :: b(3,3)
    ! -----------------------------------------------------------------------------

    b(1,1) = a(1,1)
    b(2,1) = a(1,2)
    b(3,1) = a(1,3)

    b(1,2) = a(2,1)
    b(2,2) = a(2,2)
    b(3,2) = a(2,3)

    b(1,3) = a(3,1)
    b(2,3) = a(3,2)
    b(3,3) = a(3,3)

end subroutine tensor_transpose

!===============================================================================
! subroutine ffdev_hessian_dihedrals
!===============================================================================

subroutine ffdev_hessian_dihedrals(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip,pn
    real(DEVDP)     ::  f(3),g(3),h(3),a(3),b(3)
    real(DEVDP)     ::  fg,hg,a2,b2,gv,scp,phi,f1,f2,arg,sarg,carg
    real(DEVDP)     ::  di(3),dj(3),dk(3),dl(3)
    real(DEVDP)     ::  tmp_tens(3,3),tmp_tens1(3,3),tmp_tens2(3,3)
    real(DEVDP)     ::  tmp_d1(3)
    real(DEVDP)     ::  diff,sc1
    ! -----------------------------------------------------------------------------

    ! source: 10.1002/(SICI)1096-987X(19960715)17:9<1132::AID-JCC5>3.0.CO;2-T

    geo%dih_ene = 0.0d0

    do ip=1,top%ndihedrals
        i  = top%dihedrals(ip)%ai
        j  = top%dihedrals(ip)%aj
        k  = top%dihedrals(ip)%ak
        l  = top%dihedrals(ip)%al
        ic = top%dihedrals(ip)%dt

        f(:) = geo%crd(:,i) - geo%crd(:,j)
        g(:) = geo%crd(:,j) - geo%crd(:,k)
        h(:) = geo%crd(:,l) - geo%crd(:,k)

        a(1) = f(2)*g(3) - f(3)*g(2)
        a(2) = f(3)*g(1) - f(1)*g(3)
        a(3) = f(1)*g(2) - f(2)*g(1)

        b(1) = h(2)*g(3) - h(3)*g(2)
        b(2) = h(3)*g(1) - h(1)*g(3)
        b(3) = h(1)*g(2) - h(2)*g(1)

        fg = dot_product(f,g)
        hg = dot_product(h,g)
        a2 = a(1)**2 + a(2)**2 + a(3)**2
        b2 = b(1)**2 + b(2)**2 + b(3)**2
        gv = sqrt( g(1)**2 + g(2)**2 + g(3)**2 )

        ! calculate scp and phi
        scp = (a(1)*b(1)+a(2)*b(2)+a(3)*b(3))/sqrt(a2*b2)
        if ( scp .gt.  1.0 ) then
                scp =  1.0
                phi = acos (1.0) ! const
        else if ( scp .lt. -1.0 ) then
                scp = -1.0
                phi = acos (-1.0) ! const
        else
            phi = acos ( scp )
        end if
        if( g(1)*(a(2)*b(3)-a(3)*b(2)) &
           +g(2)*(a(3)*b(1)-a(1)*b(3)) &
           +g(3)*(a(1)*b(2)-a(2)*b(1)) .gt. 0) then
                    phi = -phi
        end if

        f1 = 0.0d0
        f2 = 0.0d0
        select case(top%dihedral_types(ic)%mode)
            case(DIH_COS)
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

                ! calculate gradient
                    f1 = f1 - top%dihedral_types(ic)%v(pn) * pn * sarg
                    f2 = f2 - top%dihedral_types(ic)%v(pn) * pn*pn * carg
                end do
            case(DIH_GRBF)
                do pn=1,top%dihedral_types(ic)%n
                    if( .not. top%dihedral_types(ic)%enabled(pn) ) cycle

                ! calculate energy
                    diff = ffdev_geometry_get_dihedral_deviation(phi,top%dihedral_types(ic)%p(pn))

                    geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%c(pn) &
                                  * exp(-diff**2/top%dihedral_types(ic)%w2(pn))
                ! calculate gradient
                    f1 = f1 - 2.0*top%dihedral_types(ic)%c(pn)*diff &
                            * exp(-diff**2/top%dihedral_types(ic)%w2(pn)) &
                            / top%dihedral_types(ic)%w2(pn)
                    f2 = f2 - 2.0*top%dihedral_types(ic)%c(pn) &
                            * exp(-diff**2/top%dihedral_types(ic)%w2(pn)) &
                            * (top%dihedral_types(ic)%w2(pn) - 2.0d0*diff**2) &
                            / (top%dihedral_types(ic)%w2(pn)**2)
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented [ffdev_hessian_dihedrals]!')
        end select

    ! calculate gradient
        di(:) = -gv/a2*a(:)
        dj(:) = (gv/a2 + fg/(a2*gv))*a(:) - hg/(b2*gv)*b(:)
        dk(:) = (hg/(b2*gv) - gv/b2)*b(:) - fg/(a2*gv)*a(:)
        dl(:) = gv/b2*b(:)

        geo%grd(:,i) = geo%grd(:,i) + f1*di(:)
        geo%grd(:,j) = geo%grd(:,j) + f1*dj(:)
        geo%grd(:,k) = geo%grd(:,k) + f1*dk(:)
        geo%grd(:,l) = geo%grd(:,l) + f1*dl(:)


    ! calculate hessian - part I
        ! eq 28 - left part
        call tensor_product(di,di,tmp_tens)
        geo%hess(:,i,:,i) = geo%hess(:,i,:,i) + f2*tmp_tens(:,:)
        call tensor_product(di,dj,tmp_tens)
        geo%hess(:,i,:,j) = geo%hess(:,i,:,j) + f2*tmp_tens(:,:)
        call tensor_product(di,dk,tmp_tens)
        geo%hess(:,i,:,k) = geo%hess(:,i,:,k) + f2*tmp_tens(:,:)
        call tensor_product(di,dl,tmp_tens)
        geo%hess(:,i,:,l) = geo%hess(:,i,:,l) + f2*tmp_tens(:,:)

        call tensor_product(dj,di,tmp_tens)
        geo%hess(:,j,:,i) = geo%hess(:,j,:,i) + f2*tmp_tens(:,:)
        call tensor_product(dj,dj,tmp_tens)
        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + f2*tmp_tens(:,:)
        call tensor_product(dj,dk,tmp_tens)
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) + f2*tmp_tens(:,:)
        call tensor_product(dj,dl,tmp_tens)
        geo%hess(:,j,:,l) = geo%hess(:,j,:,l) + f2*tmp_tens(:,:)

        call tensor_product(dk,di,tmp_tens)
        geo%hess(:,k,:,i) = geo%hess(:,k,:,i) + f2*tmp_tens(:,:)
        call tensor_product(dk,dj,tmp_tens)
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) + f2*tmp_tens(:,:)
        call tensor_product(dk,dk,tmp_tens)
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + f2*tmp_tens(:,:)
        call tensor_product(dk,dl,tmp_tens)
        geo%hess(:,k,:,l) = geo%hess(:,k,:,l) + f2*tmp_tens(:,:)

        call tensor_product(dl,di,tmp_tens)
        geo%hess(:,l,:,i) = geo%hess(:,l,:,i) + f2*tmp_tens(:,:)
        call tensor_product(dl,dj,tmp_tens)
        geo%hess(:,l,:,j) = geo%hess(:,l,:,j) + f2*tmp_tens(:,:)
        call tensor_product(dl,dk,tmp_tens)
        geo%hess(:,l,:,k) = geo%hess(:,l,:,k) + f2*tmp_tens(:,:)
        call tensor_product(dl,dl,tmp_tens)
        geo%hess(:,l,:,l) = geo%hess(:,l,:,l) + f2*tmp_tens(:,:)


    ! calculate hessian - part II
    ! --------------------------------------------
        ! eq. 29
        ! df^2/dF^2, eq. 32
        sc1 = gv/(a2*a2)
        call vector_product(g,a,tmp_d1)
        call tensor_product(a,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,a,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1 + tmp_tens2)

        geo%hess(:,i,:,i) = geo%hess(:,i,:,i) + tmp_tens
        geo%hess(:,i,:,j) = geo%hess(:,i,:,j) - tmp_tens
        geo%hess(:,j,:,i) = geo%hess(:,j,:,i) - tmp_tens
        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens ! (-1)*(-1) = 1

    ! --------------------------------------------
        ! df^2/dG^2, eq. 44
        sc1 = 1.0d0/(2.0d0 * gv**3 * a2)
        call vector_product(g,a,tmp_d1)
        call tensor_product(tmp_d1,a,tmp_tens1)
        call tensor_product(a,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

        sc1 = fg/(gv * a2 * a2)
        call vector_product(f,a,tmp_d1)
        call tensor_product(a,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,a,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

        sc1 = - 1.0d0/(2.0d0 * gv**3 * b2)
        call vector_product(g,b,tmp_d1)
        call tensor_product(tmp_d1,b,tmp_tens1)
        call tensor_product(b,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

        sc1 = - hg/(gv * b2 * b2)
        call vector_product(h,b,tmp_d1)
        call tensor_product(b,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,b,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

    ! --------------------------------------------
        ! df^2/dH^2, eq. 33
        sc1 = -gv/(b2*b2)
        call vector_product(g,b,tmp_d1)
        call tensor_product(b,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,b,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1 + tmp_tens2)

        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1
        geo%hess(:,k,:,l) = geo%hess(:,k,:,l) - tmp_tens
        geo%hess(:,l,:,k) = geo%hess(:,l,:,K) - tmp_tens
        geo%hess(:,l,:,l) = geo%hess(:,l,:,l) + tmp_tens

    ! --------------------------------------------
        ! df^2/dFdG, eq. 38
        sc1 = 1.0d0 / (gv * a2 * a2)
        call vector_product(a,f,tmp_d1)
        call tensor_product(tmp_d1,a,tmp_tens1)
        call vector_product(a,g,tmp_d1)
        call tensor_product(a,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(gv*gv*tmp_tens1(:,:) + fg*tmp_tens2(:,:))
        call tensor_transpose(tmp_tens,tmp_tens1)

        geo%hess(:,i,:,j) = geo%hess(:,i,:,j) + tmp_tens
        geo%hess(:,i,:,k) = geo%hess(:,i,:,k) - tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) + tmp_tens

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) - tmp_tens - tmp_tens1

        geo%hess(:,j,:,i) = geo%hess(:,j,:,i) + tmp_tens1
        geo%hess(:,k,:,i) = geo%hess(:,k,:,i) - tmp_tens1
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) + tmp_tens1

     ! --------------------------------------------
        ! df^2/dGdH, eq. 39
        sc1 = - 1.0d0 / (gv * b2 * b2)
        call vector_product(b,h,tmp_d1)
        call tensor_product(tmp_d1,b,tmp_tens1)
        call vector_product(b,g,tmp_d1)
        call tensor_product(b,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(gv*gv*tmp_tens1(:,:) + hg*tmp_tens2(:,:))
        call tensor_transpose(tmp_tens,tmp_tens1)

        geo%hess(:,j,:,l) = geo%hess(:,j,:,l) + tmp_tens1
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens1
        geo%hess(:,k,:,l) = geo%hess(:,k,:,l) - tmp_tens1

        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens + tmp_tens1

        geo%hess(:,l,:,j) = geo%hess(:,l,:,j) + tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,l,:,k) = geo%hess(:,l,:,k) - tmp_tens

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
    real(DEVDP)     ::  f(3),g(3),h(3),a(3),b(3)
    real(DEVDP)     ::  fg,hg,a2,b2,gv,scp,phi,f1,f2,arg,sarg,carg
    real(DEVDP)     ::  di(3),dj(3),dk(3),dl(3)
    real(DEVDP)     ::  tmp_tens(3,3),tmp_tens1(3,3),tmp_tens2(3,3)
    real(DEVDP)     ::  tmp_d1(3)
    real(DEVDP)     ::  sc1
    ! -----------------------------------------------------------------------------

    geo%impropr_ene = 0.0d0

    do ip=1,top%nimpropers
        i  = top%impropers(ip)%ai
        j  = top%impropers(ip)%aj
        k  = top%impropers(ip)%ak
        l  = top%impropers(ip)%al
        ic = top%impropers(ip)%dt

        f(:) = geo%crd(:,i) - geo%crd(:,j)
        g(:) = geo%crd(:,j) - geo%crd(:,k)
        h(:) = geo%crd(:,l) - geo%crd(:,k)

        a(1) = f(2)*g(3) - f(3)*g(2)
        a(2) = f(3)*g(1) - f(1)*g(3)
        a(3) = f(1)*g(2) - f(2)*g(1)

        b(1) = h(2)*g(3) - h(3)*g(2)
        b(2) = h(3)*g(1) - h(1)*g(3)
        b(3) = h(1)*g(2) - h(2)*g(1)

        fg = dot_product(f,g)
        hg = dot_product(h,g)
        a2 = a(1)**2 + a(2)**2 + a(3)**2
        b2 = b(1)**2 + b(2)**2 + b(3)**2
        gv = sqrt( g(1)**2 + g(2)**2 + g(3)**2 )

        ! calculate scp and phi
        scp = (a(1)*b(1)+a(2)*b(2)+a(3)*b(3))/sqrt(a2*b2)
        if ( scp .gt.  1.0 ) then
                scp =  1.0
                phi = acos (1.0)    ! const
        else if ( scp .lt. -1.0 ) then
                scp = -1.0
                phi = acos (-1.0)   ! const
        else
            phi = acos ( scp )
        end if
        if( g(1)*(a(2)*b(3)-a(3)*b(2)) &
           +g(2)*(a(3)*b(1)-a(1)*b(3)) &
           +g(3)*(a(1)*b(2)-a(2)*b(1)) .gt. 0) then
                    phi = -phi
        end if

    ! calculate energy
        arg = 2.0*phi - top%improper_types(ic)%g
        sarg = sin(arg)
        carg = cos(arg)

        geo%impropr_ene = geo%impropr_ene + top%improper_types(ic)%v*(1.0d0 + carg)

    ! calculate gradient
        f1 = - 2.0d0 * top%improper_types(ic)%v * sarg
        f2 = - 4.0d0 * top%improper_types(ic)%v * carg

        di(:) = -gv/a2*a(:)
        dj(:) = (gv/a2 + fg/(a2*gv))*a(:) - hg/(b2*gv)*b(:)
        dk(:) = (hg/(b2*gv) - gv/b2)*b(:) - fg/(a2*gv)*a(:)
        dl(:) = gv/b2*b(:)

        geo%grd(:,i) = geo%grd(:,i) + f1*di(:)
        geo%grd(:,j) = geo%grd(:,j) + f1*dj(:)
        geo%grd(:,k) = geo%grd(:,k) + f1*dk(:)
        geo%grd(:,l) = geo%grd(:,l) + f1*dl(:)

    ! calculate hessian - part I
        ! eq 28 - left part
        call tensor_product(di,di,tmp_tens)
        geo%hess(:,i,:,i) = geo%hess(:,i,:,i) + f2*tmp_tens(:,:)
        call tensor_product(di,dj,tmp_tens)
        geo%hess(:,i,:,j) = geo%hess(:,i,:,j) + f2*tmp_tens(:,:)
        call tensor_product(di,dk,tmp_tens)
        geo%hess(:,i,:,k) = geo%hess(:,i,:,k) + f2*tmp_tens(:,:)
        call tensor_product(di,dl,tmp_tens)
        geo%hess(:,i,:,l) = geo%hess(:,i,:,l) + f2*tmp_tens(:,:)

        call tensor_product(dj,di,tmp_tens)
        geo%hess(:,j,:,i) = geo%hess(:,j,:,i) + f2*tmp_tens(:,:)
        call tensor_product(dj,dj,tmp_tens)
        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + f2*tmp_tens(:,:)
        call tensor_product(dj,dk,tmp_tens)
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) + f2*tmp_tens(:,:)
        call tensor_product(dj,dl,tmp_tens)
        geo%hess(:,j,:,l) = geo%hess(:,j,:,l) + f2*tmp_tens(:,:)

        call tensor_product(dk,di,tmp_tens)
        geo%hess(:,k,:,i) = geo%hess(:,k,:,i) + f2*tmp_tens(:,:)
        call tensor_product(dk,dj,tmp_tens)
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) + f2*tmp_tens(:,:)
        call tensor_product(dk,dk,tmp_tens)
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + f2*tmp_tens(:,:)
        call tensor_product(dk,dl,tmp_tens)
        geo%hess(:,k,:,l) = geo%hess(:,k,:,l) + f2*tmp_tens(:,:)

        call tensor_product(dl,di,tmp_tens)
        geo%hess(:,l,:,i) = geo%hess(:,l,:,i) + f2*tmp_tens(:,:)
        call tensor_product(dl,dj,tmp_tens)
        geo%hess(:,l,:,j) = geo%hess(:,l,:,j) + f2*tmp_tens(:,:)
        call tensor_product(dl,dk,tmp_tens)
        geo%hess(:,l,:,k) = geo%hess(:,l,:,k) + f2*tmp_tens(:,:)
        call tensor_product(dl,dl,tmp_tens)
        geo%hess(:,l,:,l) = geo%hess(:,l,:,l) + f2*tmp_tens(:,:)


    ! calculate hessian - part II
    ! --------------------------------------------
        ! eq. 29
        ! df^2/dF^2, eq. 32
        sc1 = gv/(a2*a2)
        call vector_product(g,a,tmp_d1)
        call tensor_product(a,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,a,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1 + tmp_tens2)

        geo%hess(:,i,:,i) = geo%hess(:,i,:,i) + tmp_tens
        geo%hess(:,i,:,j) = geo%hess(:,i,:,j) - tmp_tens
        geo%hess(:,j,:,i) = geo%hess(:,j,:,i) - tmp_tens
        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens ! (-1)*(-1) = 1

    ! --------------------------------------------
        ! df^2/dG^2, eq. 44
        sc1 = 1.0d0/(2.0d0 * gv**3 * a2)
        call vector_product(g,a,tmp_d1)
        call tensor_product(tmp_d1,a,tmp_tens1)
        call tensor_product(a,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

        sc1 = fg/(gv * a2 * a2)
        call vector_product(f,a,tmp_d1)
        call tensor_product(a,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,a,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

        sc1 = - 1.0d0/(2.0d0 * gv**3 * b2)
        call vector_product(g,b,tmp_d1)
        call tensor_product(tmp_d1,b,tmp_tens1)
        call tensor_product(b,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

        sc1 = - hg/(gv * b2 * b2)
        call vector_product(h,b,tmp_d1)
        call tensor_product(b,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,b,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1(:,:) + tmp_tens2(:,:))

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) + tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1

    ! --------------------------------------------
        ! df^2/dH^2, eq. 33
        sc1 = -gv/(b2*b2)
        call vector_product(g,b,tmp_d1)
        call tensor_product(b,tmp_d1,tmp_tens1)
        call tensor_product(tmp_d1,b,tmp_tens2)

        tmp_tens = f1*sc1*(tmp_tens1 + tmp_tens2)

        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens ! (-1)*(-1) = 1
        geo%hess(:,k,:,l) = geo%hess(:,k,:,l) - tmp_tens
        geo%hess(:,l,:,k) = geo%hess(:,l,:,K) - tmp_tens
        geo%hess(:,l,:,l) = geo%hess(:,l,:,l) + tmp_tens

    ! --------------------------------------------
        ! df^2/dFdG, eq. 38
        sc1 = 1.0d0 / (gv * a2 * a2)
        call vector_product(a,f,tmp_d1)
        call tensor_product(tmp_d1,a,tmp_tens1)
        call vector_product(a,g,tmp_d1)
        call tensor_product(a,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(gv*gv*tmp_tens1(:,:) + fg*tmp_tens2(:,:))
        call tensor_transpose(tmp_tens,tmp_tens1)

        geo%hess(:,i,:,j) = geo%hess(:,i,:,j) + tmp_tens
        geo%hess(:,i,:,k) = geo%hess(:,i,:,k) - tmp_tens
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) + tmp_tens

        geo%hess(:,j,:,j) = geo%hess(:,j,:,j) - tmp_tens - tmp_tens1

        geo%hess(:,j,:,i) = geo%hess(:,j,:,i) + tmp_tens1
        geo%hess(:,k,:,i) = geo%hess(:,k,:,i) - tmp_tens1
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) + tmp_tens1

     ! --------------------------------------------
        ! df^2/dGdH, eq. 39
        sc1 = - 1.0d0 / (gv * b2 * b2)
        call vector_product(b,h,tmp_d1)
        call tensor_product(tmp_d1,b,tmp_tens1)
        call vector_product(b,g,tmp_d1)
        call tensor_product(b,tmp_d1,tmp_tens2)

        tmp_tens = f1*sc1*(gv*gv*tmp_tens1(:,:) + hg*tmp_tens2(:,:))
        call tensor_transpose(tmp_tens,tmp_tens1)

        geo%hess(:,j,:,l) = geo%hess(:,j,:,l) + tmp_tens1
        geo%hess(:,j,:,k) = geo%hess(:,j,:,k) - tmp_tens1
        geo%hess(:,k,:,l) = geo%hess(:,k,:,l) - tmp_tens1

        geo%hess(:,k,:,k) = geo%hess(:,k,:,k) + tmp_tens + tmp_tens1

        geo%hess(:,l,:,j) = geo%hess(:,l,:,j) + tmp_tens
        geo%hess(:,k,:,j) = geo%hess(:,k,:,j) - tmp_tens
        geo%hess(:,l,:,k) = geo%hess(:,l,:,k) - tmp_tens

    end do

end subroutine ffdev_hessian_impropers

! ------------------------------------------------------------------------------

end module ffdev_hessian

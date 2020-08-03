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

module ffdev_gradient

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_gradient_all
! ==============================================================================

subroutine ffdev_gradient_all(top,geo,skipnb)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_timers

    use ffdev_nbmode_LJ
    use ffdev_nbmode_EXP_DISPBJ
    use ffdev_nbmode_EXP_DISPTT

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
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    geo%total_ene = 0.0d0
    geo%rst_energy = 0.0d0

    ! reset gradient
    geo%grd(:,:) = 0.0d0

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_gradient_bonds(top,geo)
        call ffdev_gradient_angles(top,geo)
        call ffdev_gradient_dihedrals(top,geo)
        call ffdev_gradient_impropers(top,geo)
    end if

    if( calcnb ) then
        if( top%nb_params_update ) then
            call ffdev_topology_update_nb_params(top)
        end if

        ! non-bonded terms
        select case(nb_mode)
            case(NB_VDW_LJ)
                call ffdev_gradient_nb_LJ(top,geo)

            case(NB_VDW_EXP_DISPBJ)
                call ffdev_gradient_nb_EXP_DISPBJ(top,geo)

            case(NB_VDW_EXP_DISPTT)
                call ffdev_gradient_nb_EXP_DISPTT(top,geo)

            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported vdW mode in ffdev_gradient_all!')
        end select
    end if

    geo%total_ene = geo%bond_ene + geo%angle_ene + geo%dih_ene + geo%impropr_ene &
                  + geo%ele14_ene + geo%rep14_ene + geo%dis14_ene  &
                  + geo%ele_ene + geo%pen_ene + geo%rep_ene + geo%dis_ene

end subroutine ffdev_gradient_all

! ==============================================================================
! subroutine ffdev_gradient_num_all
! ==============================================================================

subroutine ffdev_gradient_num_all(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------
    type(GEOMETRY)  :: tmp_geo
    real(DEVDP)     :: d,ene1,ene2
    integer         :: i,j
    ! --------------------------------------------------------------------------

    d = 5.0d-4  ! differentiation parameter

    ! calculate base energy
    call ffdev_energy_all(top,geo)

    ! allocate temporary geometry object
    call ffdev_geometry_allocate(tmp_geo,geo%natoms)
    tmp_geo%crd(:,:) = geo%crd(:,:)

    ! gradient by numerical differentiation
    do i=1,geo%natoms
        do j=1,3
            ! left
            tmp_geo%crd(j,i) = geo%crd(j,i) + d
            call ffdev_energy_all(top,tmp_geo)
            ene1 = tmp_geo%total_ene

            ! right
            tmp_geo%crd(j,i) = geo%crd(j,i) - d
            call ffdev_energy_all(top,tmp_geo)
            ene2 = tmp_geo%total_ene

            ! gradient
            geo%grd(j,i) = 0.5d0*(ene1-ene2)/d

            ! move back
            tmp_geo%crd(j,i) = geo%crd(j,i)
        end do
    end do

    ! release temporary geometry object
    deallocate(tmp_geo%crd)

end subroutine ffdev_gradient_num_all

!===============================================================================
! subroutine ffdev_gradient_bonds
!===============================================================================

subroutine ffdev_gradient_bonds(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,ib,ic
    real(DEVDP)     ::  b,db,dv
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
        b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
        db = b - top%bond_types(ic)%d0
        geo%bond_ene = geo%bond_ene + 0.5d0*top%bond_types(ic)%k*db**2

        ! calculate gradient
        dv = top%bond_types(ic)%k*db/b
        geo%grd(:,j) = geo%grd(:,j) + rij(:)*dv
        geo%grd(:,i) = geo%grd(:,i) - rij(:)*dv
    end do

end subroutine ffdev_gradient_bonds

!===============================================================================
! subroutine ffdev_gradient_angles
!===============================================================================

subroutine ffdev_gradient_angles(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,j,k,ia,ic
    real(DEVDP)     :: bjiinv, bjkinv, bji2inv, bjk2inv
    real(DEVDP)     :: scp,angv,da,dv,f1
    real(DEVDP)     :: rji(3),rjk(3),di(3),dk(3)
    ! --------------------------------------------------------------------------

    ! reset energy
    geo%angle_ene = 0.0d0

    do ia=1,top%nangles
        i  = top%angles(ia)%ai
        j  = top%angles(ia)%aj
        k  = top%angles(ia)%ak
        ic = top%angles(ia)%at

        ! calculate rji and rjk
        rji(:) = geo%crd(:,i) - geo%crd(:,j)
        rjk(:) = geo%crd(:,k) - geo%crd(:,j)

        ! calculate bjiinv and bjkinv and their squares
        bji2inv = 1.0d0/(rji(1)**2 + rji(2)**2 + rji(3)**2 )
        bjk2inv = 1.0d0/(rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
        bjiinv = sqrt(bji2inv)
        bjkinv = sqrt(bjk2inv)

            ! calculate scp and angv
        scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
        scp = scp * bjiinv*bjkinv
        if ( scp .gt.  1.0d0 ) then
            scp =  1.0d0
        else if ( scp .lt. -1.0d0 ) then
            scp = -1.0d0
        end if
        angv = acos(scp)

        ! calculate energy
        da = angv - top%angle_types(ic)%a0
        geo%angle_ene = geo%angle_ene + 0.5*top%angle_types(ic)%k*da**2

        ! calculate gradient
        dv = top%angle_types(ic)%k*da

        f1 = sin ( angv )
        if ( abs(f1) .lt. 1.0d-3 ) then
            ! sin(0.1 deg) = 1.7e-3
            ! this is set for angles close to 0 deg or 180 deg by 0.1 deg
            ! the aim is to avoid division be zero
            f1 = -1.0d3
        else
            f1 = -1.0d0 / f1
        end if

        di(:) = f1 * ( rjk(:)*bjiinv*bjkinv - scp*rji(:)*bji2inv )
        dk(:) = f1 * ( rji(:)*bjiinv*bjkinv - scp*rjk(:)*bjk2inv )

        geo%grd(:,i) = geo%grd(:,i) + dv*di(:)
        geo%grd(:,k) = geo%grd(:,k) + dv*dk(:)
        geo%grd(:,j) = geo%grd(:,j) - dv*( di(:) + dk(:) )
    end do

end subroutine ffdev_gradient_angles

!===============================================================================
! subroutine ffdev_gradient_dihedrals
!===============================================================================

subroutine ffdev_gradient_dihedrals(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip,pn
    real(DEVDP)     ::  scp,phi,arg,dv,diff
    real(DEVDP)     ::  a2,b2,gv,fg,hg
    real(DEVDP)     ::  f(3),g(3),h(3),a(3),b(3)
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

        dv = 0.0d0
        select case(top%dihedral_types(ic)%mode)
            case(DIH_COS)
                do pn=1,top%dihedral_types(ic)%n
                    if( .not. top%dihedral_types(ic)%enabled(pn) ) cycle
                    arg = pn*phi - top%dihedral_types(ic)%g(pn)
                    ! calculate energy
                    if( dih_cos_only ) then
                        geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%v(pn)*cos(arg)
                    else
                        geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%v(pn)*(1.0d0+cos(arg))
                    end if
                    ! calculate gradient
                    dv = dv - pn*top%dihedral_types(ic)%v(pn)*sin(arg)
                end do
            case(DIH_GRBF)
                do pn=1,top%dihedral_types(ic)%n
                    if( .not. top%dihedral_types(ic)%enabled(pn) ) cycle
                    diff = ffdev_geometry_get_dihedral_deviation(phi,top%dihedral_types(ic)%p(pn))
                    ! calculate energy
                    geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%c(pn) &
                                  * exp(-diff**2/top%dihedral_types(ic)%w2(pn))
                    ! calculate gradient
                    dv = dv - top%dihedral_types(ic)%c(pn) &
                            * exp(-diff**2/top%dihedral_types(ic)%w2(pn)) &
                            * 2.0*diff/top%dihedral_types(ic)%w2(pn)
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented [ffdev_gradient_dihedrals]!')
        end select

        ! calculate gradient
        geo%grd(:,i) = geo%grd(:,i) + dv*( -gv/a2*a(:) )
        geo%grd(:,j) = geo%grd(:,j) + dv*(  (gv/a2 + fg/(a2*gv))*a(:) - hg/(b2*gv)*b(:) )
        geo%grd(:,k) = geo%grd(:,k) + dv*(  (hg/(b2*gv) - gv/b2)*b(:) - fg/(a2*gv)*a(:) )
        geo%grd(:,l) = geo%grd(:,l) + dv*( gv/b2*b(:) )
    end do

end subroutine ffdev_gradient_dihedrals

!===============================================================================
! subroutine ffdev_gradient_impropers
!===============================================================================

subroutine ffdev_gradient_impropers(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip
    real(DEVDP)     ::  scp,phi,arg,dv
    real(DEVDP)     ::  a2,b2,gv,fg,hg
    real(DEVDP)     ::  f(3),g(3),h(3),a(3),b(3)
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

        ! calculate energy
        arg = 2*phi - top%improper_types(ic)%g
        geo%impropr_ene = geo%impropr_ene + top%improper_types(ic)%v * (1.0d0 + cos(arg))

        ! calculate gradient
        dv  = -2.0d0 * top%improper_types(ic)%v * sin(arg)

        ! calculate gradient
        geo%grd(:,i) = geo%grd(:,i) + dv*( -gv/a2*a(:) )
        geo%grd(:,j) = geo%grd(:,j) + dv*(  (gv/a2 + fg/(a2*gv))*a(:) - hg/(b2*gv)*b(:) )
        geo%grd(:,k) = geo%grd(:,k) + dv*(  (hg/(b2*gv) - gv/b2)*b(:) - fg/(a2*gv)*a(:) )
        geo%grd(:,l) = geo%grd(:,l) + dv*( gv/b2*b(:) )

    end do

end subroutine ffdev_gradient_impropers

! ------------------------------------------------------------------------------

end module ffdev_gradient

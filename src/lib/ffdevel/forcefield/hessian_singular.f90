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

module ffdev_hessian_singular

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_hessian_dihedrals
!===============================================================================

subroutine ffdev_hessian_dihedrals_singular(top,geo)

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

end subroutine ffdev_hessian_dihedrals_singular

!===============================================================================
! subroutine ffdev_hessian_impropers_singular
!===============================================================================

subroutine ffdev_hessian_impropers_singular(top,geo)

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

end subroutine ffdev_hessian_impropers_singular

! ------------------------------------------------------------------------------

end module ffdev_hessian_singular

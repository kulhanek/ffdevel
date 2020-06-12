! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_nbmode_EXP_DISPBJ

use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_EXP_DISPBJ
!===============================================================================

subroutine ffdev_energy_nb_EXP_DISPBJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    do ip=1,top%nb_size
        i    = top%nb_list(ip)%ai
        j    = top%nb_list(ip)%aj
        dt   = top%nb_list(ip)%dt

        ! electrostatics
        crgij = top%nb_list(ip)%crgij

        ! repulsion
        pa   = top%nb_list(ip)%pa
        pb   = top%nb_list(ip)%pb

        ! dispersion coefficients
        c6   = top%nb_list(ip)%c6
        c8   = top%nb_list(ip)%c8
        c10  = top%nb_list(ip)%c10

        rc6  = top%nb_list(ip)%rc6
        rc8  = top%nb_list(ip)%rc8
        rc10 = top%nb_list(ip)%rc10

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)

        r6   = r2*r2*r2
        r8   = r6*r2
        r10  = r8*r2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_ee =   crgij/r
        V_aa =   pa*exp(-pb*r)
        V_bb = - c6*r6i - c8*r8i - c10*r10i

        ! write(*,*) 'ENE  : ',ip,r,pa,pb,c6,c8,c10,rc6,rc8,rc10,V_aa,V_bb

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb
        end if
    end do

end subroutine ffdev_energy_nb_EXP_DISPBJ

!===============================================================================
! subroutine ffdev_energy_sapt_EXP_DISPBJ
!===============================================================================

subroutine ffdev_energy_sapt_EXP_DISPBJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee
    ! --------------------------------------------------------------------------

    geo%sapt_ele  = 0.0d0
    geo%sapt_rep  = 0.0d0
    geo%sapt_dis = 0.0d0

    do ip=1,top%sapt_size
        i     = top%sapt_list(ip)%ai
        j     = top%sapt_list(ip)%aj

        ! electrostatics
        crgij = top%sapt_list(ip)%crgij

        ! repulsion
        pa   = top%sapt_list(ip)%pa
        pb   = top%sapt_list(ip)%pb

        ! dispersion coefficients
        c6   = top%sapt_list(ip)%c6
        c8   = top%sapt_list(ip)%c8
        c10  = top%sapt_list(ip)%c10

        rc6  = top%sapt_list(ip)%rc6
        rc8  = top%sapt_list(ip)%rc8
        rc10 = top%sapt_list(ip)%rc10

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)

        r6   = r2*r2*r2
        r8   = r6*r2
        r10  = r8*r2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_ee =   crgij/r
        V_aa =   pa*exp(-pb*r)
        V_bb = - c6*r6i - c8*r8i - c10*r10i

        ! write(*,*) 'SAPT: ',ip,r,pa,pb,c6,c8,c10,rc6,rc8,rc10,V_aa,V_bb

        geo%sapt_ele = geo%sapt_ele + V_ee
        geo%sapt_rep = geo%sapt_rep + V_aa
        geo%sapt_dis = geo%sapt_dis + V_bb

    end do

end subroutine ffdev_energy_sapt_EXP_DISPBJ

!===============================================================================
! subroutine ffdev_energy_nbpair_EXP_DISPBJ
!===============================================================================

real(DEVDP) function ffdev_energy_nbpair_EXP_DISPBJ(nbpair,r)

    use ffdev_topology_dat

    implicit none
    type(NB_PAIR)   :: nbpair
    real(DEVDP)     :: r
    ! --------------------------------------------
    real(DEVDP)     :: pa,pb
    real(DEVDP)     :: r2,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i
    ! --------------------------------------------------------------------------

    ! repulsion
    pa   = nbpair%pa
    pb   = nbpair%pb

    ! dispersion coefficients
    c6   = nbpair%c6
    c8   = nbpair%c8
    c10  = nbpair%c10

    rc6  = nbpair%rc6
    rc8  = nbpair%rc8
    rc10 = nbpair%rc10

    ! calculate distances
    r2   = r**2

    r6   = r2*r2*r2
    r8   = r6*r2
    r10  = r8*r2

    r6i  = 1.0d0/(r6+rc6)
    r8i  = 1.0d0/(r8+rc8)
    r10i = 1.0d0/(r10+rc10)

    V_aa =   pa*exp(-pb*r)
    V_bb = - c6*r6i - c8*r8i - c10*r10i

    ffdev_energy_nbpair_EXP_DISPBJ = V_aa + V_bb

end function ffdev_energy_nbpair_EXP_DISPBJ

!===============================================================================
! subroutine ffdev_gradient_nb_EXP_DISPBJ
!===============================================================================

subroutine ffdev_gradient_nb_EXP_DISPBJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3,r2a
    real(DEVDP)     :: r2,r,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee,dva,inv_scnb,inv_scee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    do ip=1,top%nb_size
        i    = top%nb_list(ip)%ai
        j    = top%nb_list(ip)%aj
        dt   = top%nb_list(ip)%dt

        ! electrostatics
        crgij = top%nb_list(ip)%crgij

        ! repulsion
        pa   = top%nb_list(ip)%pa
        pb   = top%nb_list(ip)%pb

        ! dispersion coefficients
        c6   = top%nb_list(ip)%c6
        c8   = top%nb_list(ip)%c8
        c10  = top%nb_list(ip)%c10

        rc6  = top%nb_list(ip)%rc6
        rc8  = top%nb_list(ip)%rc8
        rc10 = top%nb_list(ip)%rc10

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a  = 1.0d0/r2
        r    = sqrt(r2)

        r6   = r2*r2*r2
        r8   = r6*r2
        r10  = r8*r2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_ee =   crgij/r
        V_aa =   pa*exp(-pb*r)
        V_bb = - c6*r6i - c8*r8i - c10*r10i

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb

            dva = r2a*(V_ee +  V_aa*pb*r &
                       - 6.0d0*c6*r6i*r6i*r6 - 8.0d0*c8*r8i*r8i*r8 - 10.0d0*c10*r10i*r10i*r10)
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb

            dva = r2a*(inv_scee*V_ee + inv_scnb*( V_aa*pb*r &
                       - 6.0d0*c6*r6i*r6i*r6 - 8.0d0*c8*r8i*r8i*r8 - 10.0d0*c10*r10i*r10i*r10))
        end if

        ! calculate gradient
        dxa1 = dva*dxa1
        dxa2 = dva*dxa2
        dxa3 = dva*dxa3

        geo%grd(1,i) = geo%grd(1,i) - dxa1
        geo%grd(2,i) = geo%grd(2,i) - dxa2
        geo%grd(3,i) = geo%grd(3,i) - dxa3

        geo%grd(1,j) = geo%grd(1,j) + dxa1
        geo%grd(2,j) = geo%grd(2,j) + dxa2
        geo%grd(3,j) = geo%grd(3,j) + dxa3
    end do

end subroutine ffdev_gradient_nb_EXP_DISPBJ

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_DISPBJ

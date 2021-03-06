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

module ffdev_nbmode_LJ

use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_LJ
!===============================================================================

subroutine ffdev_energy_nb_LJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r2a,ra,r6a,V_aa,V_bb,V_ee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    do ip=1,top%nb_size
        i     = top%nb_list(ip)%ai
        j     = top%nb_list(ip)%aj
        dt    = top%nb_list(ip)%dt

        crgij = top%nb_list(ip)%q1 * top%nb_list(ip)%q2
        aLJa  = top%nb_list(ip)%pa1
        bLJa  = top%nb_list(ip)%c6

        dxa1  = geo%crd(1,i) - geo%crd(1,j)
        dxa2  = geo%crd(2,i) - geo%crd(2,j)
        dxa3  = geo%crd(3,i) - geo%crd(3,j)

        r2    = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a   = 1.0d0/r2
        ra    = sqrt(r2a)
        r6a   = r2a*r2a*r2a

        V_ee  =   crgij*ra
        V_aa  =   aLJa*r6a*r6a
        V_bb  = - bLJa*r6a

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

end subroutine ffdev_energy_nb_LJ

!===============================================================================
! subroutine ffdev_energy_nb_LJ_EXP6_probe
!===============================================================================

subroutine ffdev_energy_nb_LJ_EXP6_probe(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,aLJb,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r,r2,r2a,ra,r6a,V_aa,V_bb,V_ee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    do ip=1,top%nb_size
        i     = top%nb_list(ip)%ai
        j     = top%nb_list(ip)%aj
        dt    = top%nb_list(ip)%dt

        crgij = top%nb_list(ip)%q1 * top%nb_list(ip)%q2
        aLJa  = top%nb_list(ip)%pa1
        aLJb  = top%nb_list(ip)%pb1
        bLJa  = top%nb_list(ip)%c6

        dxa1  = geo%crd(1,i) - geo%crd(1,j)
        dxa2  = geo%crd(2,i) - geo%crd(2,j)
        dxa3  = geo%crd(3,i) - geo%crd(3,j)

        r2    = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r     = sqrt(r2)
        ra    = 1.0d0/r
        r2a   = 1.0d0/r2

        r6a   = r2a*r2a*r2a

        V_ee  =   crgij*ra
        V_aa  =   aLJa*exp(-aLJb*r)
        V_bb  = - bLJa*r6a

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

end subroutine ffdev_energy_nb_LJ_EXP6_probe

!===============================================================================
! subroutine ffdev_energy_sapt_LJ
!===============================================================================

subroutine ffdev_energy_sapt_LJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j
    real(DEVDP)     :: aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r2a,ra,r6a,V_ee,V_aa,V_bb
    ! --------------------------------------------------------------------------

    geo%sapt_ele = 0.0d0
    geo%sapt_pen = 0.0d0
    geo%sapt_rep = 0.0d0
    geo%sapt_dis = 0.0d0

    do ip=1,top%sapt_size
        i     = top%sapt_list(ip)%ai
        j     = top%sapt_list(ip)%aj

        crgij = top%sapt_list(ip)%q1 * top%sapt_list(ip)%q2
        aLJa  = top%sapt_list(ip)%pa1
        bLJa  = top%sapt_list(ip)%c6

        dxa1  = geo%crd(1,i) - geo%crd(1,j)
        dxa2  = geo%crd(2,i) - geo%crd(2,j)
        dxa3  = geo%crd(3,i) - geo%crd(3,j)

        r2    = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a   = 1.0d0/r2
        ra    = sqrt(r2a)
        r6a   = r2a*r2a*r2a

        V_ee  =   crgij*ra
        V_aa  =   aLJa*r6a*r6a
        V_bb  = - bLJa*r6a

        geo%sapt_ele = geo%sapt_ele + V_ee
        geo%sapt_rep = geo%sapt_rep + V_aa
        geo%sapt_dis = geo%sapt_dis + V_bb
    end do

end subroutine ffdev_energy_sapt_LJ

!===============================================================================
! subroutine ffdev_energy_nbpair_LJ
!===============================================================================

subroutine ffdev_energy_nbpair_LJ(nbpair,r,nbene)

    use ffdev_topology_dat

    implicit none
    type(NB_PAIR)           :: nbpair
    real(DEVDP)             :: r
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------
    real(DEVDP)     :: aLJa,bLJa
    real(DEVDP)     :: r2,r2a,r6a,V_aa,V_bb,V_ee,V_pe
    ! --------------------------------------------------------------------------

    nbene%ele_ene = 0.0
    nbene%pen_ene = 0.0
    nbene%rep_ene = 0.0
    nbene%dis_ene = 0.0
    nbene%tot_ene = 0.0

    aLJa  = nbpair%pa1
    bLJa  = nbpair%c6

    r2    = r**2
    r2a   = 1.0d0/r2
    r6a   = r2a*r2a*r2a

    V_ee  = 0.0d0
    V_pe  = 0.0d0
    V_aa  =   aLJa*r6a*r6a
    V_bb  = - bLJa*r6a

    nbene%ele_ene = V_ee
    nbene%pen_ene = V_pe
    nbene%rep_ene = V_aa
    nbene%dis_ene = V_bb
    nbene%tot_ene = V_ee + V_pe + V_aa + V_bb

end subroutine ffdev_energy_nbpair_LJ

!===============================================================================
! subroutine ffdev_gradient_nb_LJ
!===============================================================================

subroutine ffdev_gradient_nb_LJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r2a,ra,r6a,V_ee,V_aa,V_bb,dva
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    do ip=1,top%nb_size
        i     = top%nb_list(ip)%ai
        j     = top%nb_list(ip)%aj
        dt    = top%nb_list(ip)%dt

        crgij = top%nb_list(ip)%q1 * top%nb_list(ip)%q2
        aLJa  = top%nb_list(ip)%pa1
        bLJa  = top%nb_list(ip)%c6

        dxa1  = geo%crd(1,i) - geo%crd(1,j)
        dxa2  = geo%crd(2,i) - geo%crd(2,j)
        dxa3  = geo%crd(3,i) - geo%crd(3,j)

        r2    = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a   = 1.0d0/r2
        ra    = sqrt(r2a)
        r6a   = r2a*r2a*r2a

        V_ee  =   crgij*ra
        V_aa  =   aLJa*r6a*r6a
        V_bb  = - bLJa*r6a

        ! calculate energy
        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb

            dva = r2a*(V_ee + 12.0d0*V_aa + 6.0d0*V_bb)
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb

            dva = r2a*(inv_scee*V_ee + inv_scnb*(12.0d0*V_aa + 6.0d0*V_bb))
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

end subroutine ffdev_gradient_nb_LJ

! ------------------------------------------------------------------------------

end module ffdev_nbmode_LJ

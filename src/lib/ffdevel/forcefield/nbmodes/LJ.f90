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

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_LJ_qtop
!===============================================================================

subroutine ffdev_energy_nb_LJ_qtop(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,scale2,V_aa,V_bb
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0

    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i   = top%nb_list(ip)%ai
        j   = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        aLJa  = top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**12
        bLJa  = 2.0d0*top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**6
        crgij = top%atoms(i)%charge * top%atoms(j)%charge

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        V_aa =   aLJa*r6a*r6a
        V_bb = - bLJa*r6a

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene  + scale2*crgij*ra
            geo%nb_ene   = geo%nb_ene   + V_aa + V_bb

            geo%nb_rep   = geo%nb_rep   + V_aa
            geo%nb_disp  = geo%nb_disp  + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
            geo%ele14_ene   = geo%ele14_ene + inv_scee*scale2*crgij*ra
            geo%nb14_ene    = geo%nb14_ene  + inv_scnb*(V_aa + V_bb)

            geo%nb_rep      = geo%nb_rep    + inv_scnb*V_aa
            geo%nb_disp     = geo%nb_disp   + inv_scnb*V_bb
        end if
    end do

end subroutine ffdev_energy_nb_LJ_qtop

!===============================================================================
! subroutine ffdev_energy_nb_LJ_qgeo
!===============================================================================

subroutine ffdev_energy_nb_LJ_qgeo(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,scale2,V_aa,V_bb
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0

    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i   = top%nb_list(ip)%ai
        j   = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        aLJa  = top%nb_types(nbt)%eps * top%nb_types(nbt)%r0**12
        bLJa  = 2.0d0 * top%nb_types(nbt)%eps * top%nb_types(nbt)%r0**6
        crgij = geo%sup_chrg(i) * geo%sup_chrg(j)

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        V_aa =   aLJa*r6a*r6a
        V_bb = - bLJa*r6a

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene  + scale2*crgij*ra
            geo%nb_ene   = geo%nb_ene   +  V_aa + V_bb
            geo%nb_rep   = geo%nb_rep   + V_aa
            geo%nb_disp  = geo%nb_disp  + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene   = geo%ele14_ene + inv_scee*scale2*crgij*ra
            geo%nb14_ene    = geo%nb14_ene  + inv_scnb*(V_aa + V_bb)
            geo%nb_rep      = geo%nb_rep    + inv_scnb*V_aa
            geo%nb_disp     = geo%nb_disp   + inv_scnb*V_bb
        end if
    end do

end subroutine ffdev_energy_nb_LJ_qgeo

!===============================================================================
! subroutine ffdev_energy_sapt_LJ
!===============================================================================

subroutine ffdev_energy_sapt_LJ(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,scale2
    ! --------------------------------------------------------------------------

    geo%sapt_ele = 0.0d0
    geo%sapt_rep = 0.0d0
    geo%sapt_disp = 0.0d0

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%sapt_size
        i   = top%sapt_list(ip)%ai
        j   = top%sapt_list(ip)%aj
        nbt = top%sapt_list(ip)%nbt

        ! write(*,*) ip,i,j,nbt

        aLJa  = top%nb_types(nbt)%eps * top%nb_types(nbt)%r0**12
        bLJa  = 2.0d0*top%nb_types(nbt)%eps * top%nb_types(nbt)%r0**6

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij = geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij = top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        geo%sapt_ele  = geo%sapt_ele  + scale2*crgij*ra
        geo%sapt_rep  = geo%sapt_rep  + aLJa*r6a*r6a
        geo%sapt_disp = geo%sapt_disp - bLJa*r6a
    end do

end subroutine ffdev_energy_sapt_LJ

!===============================================================================
! subroutine ffdev_gradient_nb_lj_qtop
!===============================================================================

subroutine ffdev_gradient_nb_lj_qtop(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,Vela,V_aa,V_ba,dva,scale2
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0

    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt
        aLJa  = top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**12
        bLJa  = 2.0d0*top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**6
        crgij = top%atoms(i)%charge * top%atoms(j)%charge

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        Vela = scale2*crgij*ra
        V_aa = aLJa*r6a*r6a
        V_ba = bLJa*r6a

        ! calculate energy
        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + Vela
            geo%nb_ene   = geo%nb_ene + V_aa - V_ba

            geo%nb_rep   = geo%nb_rep   + V_aa
            geo%nb_disp  = geo%nb_disp  - V_ba

            dva = r2a*(Vela + 12.0d0*V_aa - 6.0d0*V_ba)
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene   = geo%ele14_ene + inv_scee*Vela
            geo%nb14_ene    = geo%nb14_ene  + inv_scnb*(V_aa - V_ba)

            geo%nb_rep      = geo%nb_rep    + inv_scnb*V_aa
            geo%nb_disp     = geo%nb_disp   - inv_scnb*V_ba

            dva = r2a*(inv_scee*Vela + inv_scnb*(12.0d0*V_aa - 6.0d0*V_ba))
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

end subroutine ffdev_gradient_nb_lj_qtop

!===============================================================================
! subroutine ffdev_gradient_nb_lj_qgeo
!===============================================================================

subroutine ffdev_gradient_nb_lj_qgeo(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,Vela,V_aa,V_ba,dva,scale2
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0

    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt
        aLJa  = top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**12
        bLJa  = 2.0d0*top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**6
        crgij =  geo%sup_chrg(i) * geo%sup_chrg(j)

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        Vela = scale2*crgij*ra
        V_aa = aLJa*r6a*r6a
        V_ba = bLJa*r6a

        ! calculate energy
        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + Vela
            geo%nb_ene   = geo%nb_ene + V_aa - V_ba

            geo%nb_rep   = geo%nb_rep   + V_aa
            geo%nb_disp  = geo%nb_disp  - V_ba

            dva = r2a*(Vela + 12.0d0*V_aa - 6.0d0*V_ba)
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene   = geo%ele14_ene + inv_scee*Vela
            geo%nb14_ene    = geo%nb14_ene  + inv_scnb*(V_aa - V_ba)

            geo%nb_rep      = geo%nb_rep    + inv_scnb*V_aa
            geo%nb_disp     = geo%nb_disp   - inv_scnb*V_ba

            dva   = r2a*(inv_scee*Vela + inv_scnb*(12.0d0*V_aa - 6.0d0*V_ba))
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

end subroutine ffdev_gradient_nb_lj_qgeo

! ------------------------------------------------------------------------------

end module ffdev_nbmode_LJ
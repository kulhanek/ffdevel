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

module ffdev_nbmode_12_6

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_12_6
!===============================================================================

subroutine ffdev_energy_nb_12_6(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,pa,c6,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,scale2,V_aa,V_bb
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0

    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_12_6!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i   = top%nb_list(ip)%ai
        j   = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt
        pa  = exp(top%nb_types(nbt)%pa)
        c6  = top%nb_types(nbt)%c6 * disp_fa

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij =  geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij =  top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        V_aa =   pa*r6a*r6a
        V_bb = - c6*r6a

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij*ra
            geo%nb_ene   = geo%nb_ene  + V_aa + V_bb
            geo%nb_rep   = geo%nb_rep  + V_aa
            geo%nb_disp  = geo%nb_disp + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee*scale2*crgij*ra
            geo%nb14_ene  = geo%nb14_ene  + inv_scnb*(V_aa + V_bb)
            geo%nb_rep    = geo%nb_rep    + inv_scnb*V_aa
            geo%nb_disp   = geo%nb_disp   + inv_scnb*V_bb
        end if
    end do

end subroutine ffdev_energy_nb_12_6

!===============================================================================
! subroutine ffdev_energy_sapt0_12_6
!===============================================================================

subroutine ffdev_energy_sapt0_12_6(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: pa,c6,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,scale2
    ! --------------------------------------------------------------------------

    geo%sapt0_ele = 0.0d0
    geo%sapt0_rep = 0.0d0
    geo%sapt0_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_12_6!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%sapt0_size
        i   = top%sapt0_list(ip)%ai
        j   = top%sapt0_list(ip)%aj
        nbt = top%sapt0_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa)
        c6  = top%nb_types(nbt)%c6 * disp_fa

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij =  geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij =  top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)
        r6a = r2a*r2a*r2a

        geo%sapt0_ele  = geo%sapt0_ele  + scale2*crgij*ra
        geo%sapt0_rep  = geo%sapt0_rep  + pa*r6a*r6a
        geo%sapt0_disp = geo%sapt0_disp - c6*r6a
    end do

end subroutine ffdev_energy_sapt0_12_6

! ------------------------------------------------------------------------------

end module ffdev_nbmode_12_6

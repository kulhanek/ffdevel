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

module ffdev_nbmode_EXP_XDMBJ

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_EXP_XDMBJ
!===============================================================================

subroutine ffdev_energy_nb_EXP_XDMBJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,rc,rc2,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,pe
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0

    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_EXP_XDMBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt   = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa * top%nb_types(nbt)%pb)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! XDM
        c6  = xdm_pairs(agti,agtj)%c6ave
        c8  = xdm_pairs(agti,agtj)%c8ave
        c10 = xdm_pairs(agti,agtj)%c10ave

        rc  = disp_fa*xdm_pairs(agti,agtj)%rc + disp_fb

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij = geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij = top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2 = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r  = sqrt(r2)

        rc2 = rc*rc

        r6 = r2*r2*r2
        rc6 = rc2*rc2*rc2

        r8 = r6*r2
        rc8 = rc6*rc2

        r10 = r8*r2
        rc10 = rc8*rc2

        V_aa = pa*exp(-pb*r)

        r6i = 1.0d0/(r6+rc6)
        r8i = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_bb = - c6*r6i - c8*r8i - c10*r10i

        pe = (1.0d0 - exp(-disp_fc*r/xdm_pairs(agti,agtj)%rc))

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + pe * scale2*crgij/r
            geo%nb_ene  = geo%nb_ene  + V_aa + V_bb

            geo%nb_rep  = geo%nb_rep  + V_aa
            geo%nb_disp = geo%nb_disp + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + pe * inv_scee*scale2*crgij/r
            geo%nb14_ene  = geo%nb14_ene  + inv_scnb*(V_aa + V_bb)

            geo%nb_rep    = geo%nb_rep    + inv_scnb*V_aa
            geo%nb_disp   = geo%nb_disp   + inv_scnb*V_bb
        end if
    end do

end subroutine ffdev_energy_nb_EXP_XDMBJ

!===============================================================================
! subroutine ffdev_energy_sapt_EXP_XDMBJ
!===============================================================================

subroutine ffdev_energy_sapt_EXP_XDMBJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,rc,rc2,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,pe
    ! --------------------------------------------------------------------------

    geo%sapt_ele = 0.0d0
    geo%sapt_rep = 0.0d0
    geo%sapt_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_sapt_EXP_XDMBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%sapt_size
        i   = top%sapt_list(ip)%ai
        j   = top%sapt_list(ip)%aj
        nbt   = top%sapt_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa * top%nb_types(nbt)%pb)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! XDM
        c6  = xdm_pairs(agti,agtj)%c6ave
        c8  = xdm_pairs(agti,agtj)%c8ave
        c10 = xdm_pairs(agti,agtj)%c10ave

        rc  = disp_fa*xdm_pairs(agti,agtj)%rc + disp_fb

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij = geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij = top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2 = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r  = sqrt(r2)

        rc2 = rc*rc

        r6 = r2*r2*r2
        rc6 = rc2*rc2*rc2

        r8 = r6*r2
        rc8 = rc6*rc2

        r10 = r8*r2
        rc10 = rc8*rc2

        V_aa = pa*exp(-pb*r)

        r6i = 1.0d0/(r6+rc6)
        r8i = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_bb = - c6*r6i - c8*r8i - c10*r10i

        pe = (1.0d0 - exp(-disp_fc*r/xdm_pairs(agti,agtj)%rc))

        geo%sapt_ele  = geo%sapt_ele + pe * scale2*crgij/r
        geo%sapt_rep  = geo%sapt_rep + V_aa
        geo%sapt_disp = geo%sapt_disp + V_bb

    end do

end subroutine ffdev_energy_sapt_EXP_XDMBJ

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_XDMBJ

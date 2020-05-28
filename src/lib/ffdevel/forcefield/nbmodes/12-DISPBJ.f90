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

module ffdev_nbmode_12_DISPBJ

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_12_DISPBJ
!===============================================================================

subroutine ffdev_energy_nb_12_DISPBJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,rc,rc2,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    if( .not. disp_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'DISP not loaded for ffdev_energy_nb_12_DISPBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i   = top%nb_list(ip)%ai
        j   = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! dispersion coefficients
        c6  = disp_s6  * disp_pairs(agti,agtj)%c6
        c8  = disp_s8  * disp_pairs(agti,agtj)%c8
        c10 = disp_s10 * disp_pairs(agti,agtj)%c10

        select case(dispbj_mode)
            case(DISP_BJ_DRC)
                rc = damp_fa * disp_pairs(agti,agtj)%rc + damp_fb
            case(DISP_BJ_ORC)
                rc = top%nb_types(nbt)%rc
            case default
                if( .not. disp_data_loaded ) then
                    call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_gradient_nb_12_DISPBJ!')
                end if
        end select

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij = geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij = top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)

        rc2  = rc*rc

        r6   = r2*r2*r2
        rc6  = rc2*rc2*rc2

        r8   = r6*r2
        rc8  = rc6*rc2

        r10  = r8*r2
        rc10 = rc8*rc2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_ee =   scale2*crgij/r
        V_aa =   pa/(r6*r6)
        V_bb = - c6*r6i - c8*r8i - c10*r10i

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb
        end if
    end do

end subroutine ffdev_energy_nb_12_DISPBJ

!===============================================================================
! subroutine ffdev_energy_sapt_12_DISPBJ
!===============================================================================

subroutine ffdev_energy_sapt_12_DISPBJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,rc,rc2,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee
    ! --------------------------------------------------------------------------

    geo%sapt_ele  = 0.0d0
    geo%sapt_rep  = 0.0d0
    geo%sapt_dis = 0.0d0

    if( .not. disp_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'DISP not loaded for ffdev_energy_sapt_12_DISPBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%sapt_size
        i   = top%sapt_list(ip)%ai
        j   = top%sapt_list(ip)%aj
        nbt = top%sapt_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! dispersion coefficients
        c6  = disp_s6  * disp_pairs(agti,agtj)%c6
        c8  = disp_s8  * disp_pairs(agti,agtj)%c8
        c10 = disp_s10 * disp_pairs(agti,agtj)%c10

        select case(dispbj_mode)
            case(DISP_BJ_DRC)
                rc = damp_fa * disp_pairs(agti,agtj)%rc + damp_fb
            case(DISP_BJ_ORC)
                rc = top%nb_types(nbt)%rc
            case default
                if( .not. disp_data_loaded ) then
                    call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_gradient_nb_12_DISPBJ!')
                end if
        end select

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij = geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij = top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2  = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r   = sqrt(r2)

        rc2 = rc*rc

        r6  = r2*r2*r2
        rc6 = rc2*rc2*rc2

        r8  = r6*r2
        rc8 = rc6*rc2

        r10  = r8*r2
        rc10 = rc8*rc2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_ee =   scale2*crgij/r
        V_aa =   pa/(r6*r6)
        V_bb = - c6*r6i - c8*r8i - c10*r10i

        geo%sapt_ele = geo%sapt_ele + V_ee
        geo%sapt_rep = geo%sapt_rep + V_aa
        geo%sapt_dis = geo%sapt_dis + V_bb

    end do

end subroutine ffdev_energy_sapt_12_DISPBJ

!===============================================================================
! subroutine ffdev_gradient_nb_12_DISPBJ
!===============================================================================

subroutine ffdev_gradient_nb_12_DISPBJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3,r2a
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,rc,rc2,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee,dva,inv_scnb,inv_scee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    if( .not. disp_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'DISP not loaded for ffdev_gradient_nb_12_DISPBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i   = top%nb_list(ip)%ai
        j   = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! dispersion coefficients
        c6  = disp_s6  * disp_pairs(agti,agtj)%c6
        c8  = disp_s8  * disp_pairs(agti,agtj)%c8
        c10 = disp_s10 * disp_pairs(agti,agtj)%c10

        select case(dispbj_mode)
            case(DISP_BJ_DRC)
                rc = damp_fa * disp_pairs(agti,agtj)%rc + damp_fb
            case(DISP_BJ_ORC)
                rc = top%nb_types(nbt)%rc
            case default
                if( .not. disp_data_loaded ) then
                    call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_gradient_nb_12_DISPBJ!')
                end if
        end select

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij = geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij = top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2  = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2
        r   = sqrt(r2)

        rc2 = rc*rc

        r6  = r2*r2*r2
        rc6 = rc2*rc2*rc2

        r8  = r6*r2
        rc8 = rc6*rc2

        r10  = r8*r2
        rc10 = rc8*rc2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_ee =   scale2*crgij/r
        V_aa =   pa/(r6*r6)
        V_bb = - c6*r6i - c8*r8i - c10*r10i

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb

            dva = r2a*(V_ee + 12.0d0*V_aa &
                       - 6.0d0*c6*r6i*r6i*r6 - 8.0d0*c8*r8i*r8i*r8 - 10.0d0*c10*r10i*r10i*r10)
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb

            dva = r2a*(inv_scee*V_ee + inv_scnb*(12.0d0*V_aa &
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

end subroutine ffdev_gradient_nb_12_DISPBJ

! ------------------------------------------------------------------------------

end module ffdev_nbmode_12_DISPBJ
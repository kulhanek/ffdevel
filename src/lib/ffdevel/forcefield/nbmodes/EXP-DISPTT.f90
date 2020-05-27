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

module ffdev_nbmode_EXP_DISPTT

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_EXP_DISPTT
!===============================================================================

subroutine ffdev_energy_nb_EXP_DISPTT(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj, k
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3,upe
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,fd6,fd8,fd10,pe,arg,sump
    real(DEVDP)     :: V_aa, V_bb, V_ee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%pel14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%pel_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    if( .not. disp_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'DISP not loaded for ffdev_energy_nb_EXP_DISPBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i   = top%nb_list(ip)%ai
        j   = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa * top%nb_types(nbt)%pb)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! dispersion
        c6  = disp_s6  * disp_pairs(agti,agtj)%c6
        c8  = disp_s8  * disp_pairs(agti,agtj)%c8
        c10 = disp_s10 * disp_pairs(agti,agtj)%c10

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij =  geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij =  top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2 = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r  = sqrt(r2)

        upe = exp(-pb*r)

        arg = damp_fa * pb*r
        pe = exp(-arg)

    ! 6
        r6 = r2*r2*r2
        fd6 = 1.0d0 - pe
        sump = 1.0d0
        do k=1,6
            sump = sump * arg / real(k,DEVDP)
            fd6 = fd6 - pe*sump
        end do

    ! 8
        r8 = r6*r2
        fd8 = fd6
        sump = sump * arg / real(7,DEVDP)
        fd8 = fd8 - pe*sump
        sump = sump * arg / real(8,DEVDP)
        fd8 = fd8 - pe*sump

    ! 10
        r10 = r8*r2
        fd10 = fd8
        sump = sump * arg / real(9,DEVDP)
        fd10 = fd10 - pe*sump
        sump = sump * arg / real(10,DEVDP)
        fd10 = fd10 - pe*sump

        V_ee = scale2*crgij/r
        V_aa = pa*upe
        V_bb = - fd6*c6/r6 - fd8*c8/r8 - fd10*c10/r10

        ! FIXME
        pe = (1.0d0 - exp(-pene_fa*r/disp_pairs(agti,agtj)%rc))

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%pel_ene = geo%pel_ene + pe * V_ee
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%pel14_ene = geo%pel14_ene + inv_scee * pe * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb
        end if

    end do

end subroutine ffdev_energy_nb_EXP_DISPTT

!===============================================================================
! subroutine ffdev_energy_sapt_EXP_DISPTT
!===============================================================================

subroutine ffdev_energy_sapt_EXP_DISPTT(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj, k
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10
    real(DEVDP)     :: V_aa,V_bb,pe,V_ee,arg,upe,sump
    real(DEVDP)     :: fd6,fd8,fd10
    ! --------------------------------------------------------------------------

    geo%sapt_ele = 0.0d0
    geo%sapt_rep = 0.0d0
    geo%sapt_dis = 0.0d0

    if( .not. disp_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'DISP not loaded for ffdev_energy_sapt_EXP_DISPBJ!')
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

        ! dispersion
        c6  = disp_s6  * disp_pairs(agti,agtj)%c6
        c8  = disp_s8  * disp_pairs(agti,agtj)%c8
        c10 = disp_s10 * disp_pairs(agti,agtj)%c10

        if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
            crgij =  geo%sup_chrg(i) * geo%sup_chrg(j)
        else
            crgij =  top%atoms(i)%charge * top%atoms(j)%charge
        end if

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2 = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r  = sqrt(r2)

        upe = exp(-pb*r)

        arg = damp_fa * pb*r
        pe = exp(-arg)

    ! 6
        r6 = r2*r2*r2
        fd6 = 1.0d0 - pe
        sump = 1.0d0
        do k=1,6
            sump = sump * arg / real(k,DEVDP)
            fd6 = fd6 - pe*sump
        end do

    ! 8
        r8 = r6*r2
        fd8 = fd6
        sump = sump * arg / real(7,DEVDP)
        fd8 = fd8 - pe*sump
        sump = sump * arg / real(8,DEVDP)
        fd8 = fd8 - pe*sump

    ! 10
        r10 = r8*r2
        fd10 = fd8
        sump = sump * arg / real(9,DEVDP)
        fd10 = fd10 - pe*sump
        sump = sump * arg / real(10,DEVDP)
        fd10 = fd10 - pe*sump

        V_ee = scale2*crgij/r
        V_aa = pa*upe
        V_bb = - fd6*c6/r6 - fd8*c8/r8 - fd10*c10/r10

        ! FIXME
        pe = (1.0d0 - exp(-pene_fa*r/disp_pairs(agti,agtj)%rc))

        geo%sapt_ele = geo%sapt_ele + V_ee * (1.0d0 + pe)
        geo%sapt_rep = geo%sapt_rep + V_aa
        geo%sapt_dis = geo%sapt_dis + V_bb

    end do

end subroutine ffdev_energy_sapt_EXP_DISPTT

!===============================================================================
! subroutine ffdev_gradient_nb_EXP_DISPTT
!===============================================================================

subroutine ffdev_gradient_nb_EXP_DISPTT(top,geo)

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
    geo%pel14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%pel_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

!    scale2 = ele_qscale*ele_qscale*332.05221729d0
!
!    do ip=1,top%nb_size
!        i = top%nb_list(ip)%ai
!        j = top%nb_list(ip)%aj
!        nbt = top%nb_list(ip)%nbt
!        aLJa  = top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**12
!        bLJa  = 2.0d0*top%nb_types(nbt)%eps*top%nb_types(nbt)%r0**6
!        crgij =  geo%sup_chrg(i) * geo%sup_chrg(j)
!
!        ! calculate dx, r and r2
!        dxa1 = geo%crd(1,i) - geo%crd(1,j)
!        dxa2 = geo%crd(2,i) - geo%crd(2,j)
!        dxa3 = geo%crd(3,i) - geo%crd(3,j)
!
!        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
!        r2a = 1.0d0/r2a
!        ra  = sqrt(r2a)
!        r6a = r2a*r2a*r2a
!
!        Vela = scale2*crgij*ra
!        V_aa = aLJa*r6a*r6a
!        V_ba = bLJa*r6a
!
!        ! calculate energy
!        if( top%nb_list(ip)%dt .eq. 0 ) then
!            geo%ele_ene  = geo%ele_ene + Vela
!            geo%nb_ene   = geo%nb_ene + V_aa - V_ba
!
!            geo%nb_rep   = geo%nb_rep   + V_aa
!            geo%nb_disp  = geo%nb_disp  - V_ba
!
!            dva = r2a*(Vela + 12.0d0*V_aa - 6.0d0*V_ba)
!        else
!            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
!            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
!
!            geo%ele14_ene   = geo%ele14_ene + inv_scee*Vela
!            geo%nb14_ene    = geo%nb14_ene  + inv_scnb*(V_aa - V_ba)
!
!            geo%nb_rep      = geo%nb_rep    + inv_scnb*V_aa
!            geo%nb_disp     = geo%nb_disp   - inv_scnb*V_ba
!
!            dva   = r2a*(inv_scee*Vela + inv_scnb*(12.0d0*V_aa - 6.0d0*V_ba))
!        end if
!
!        ! calculate gradient
!        dxa1 = dva*dxa1
!        dxa2 = dva*dxa2
!        dxa3 = dva*dxa3
!        geo%grd(1,i) = geo%grd(1,i) - dxa1
!        geo%grd(2,i) = geo%grd(2,i) - dxa2
!        geo%grd(3,i) = geo%grd(3,i) - dxa3
!        geo%grd(1,j) = geo%grd(1,j) + dxa1
!        geo%grd(2,j) = geo%grd(2,j) + dxa2
!        geo%grd(3,j) = geo%grd(3,j) + dxa3
!    end do

end subroutine ffdev_gradient_nb_EXP_DISPTT

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_DISPTT

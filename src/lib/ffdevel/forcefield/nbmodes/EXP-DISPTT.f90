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

module ffdev_nbmode_EXP_XDMTT

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_EXP_XDMTT
!===============================================================================

subroutine ffdev_energy_nb_EXP_XDMTT(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj, k
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3,upe
    real(DEVDP)     :: r2,r,r6,r8,r10,scale2,c6,c8,c10,fd6,fd8,fd10,pe,arg, sump
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0
    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_EXP_XDMTT!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt
        pa  = exp(top%nb_types(nbt)%pa * top%nb_types(nbt)%pb)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! DEBUG
        ! write(*,*) agti, trim(top%atom_types(top%atoms(i)%typeid)%name), agtj, trim(top%atom_types(top%atoms(j)%typeid)%name)

        ! XDM
        c6  = disp_s6  * xdm_pairs(agti,agtj)%c6ave
        c8  = disp_s8  * xdm_pairs(agti,agtj)%c8ave
        c10 = disp_s10 * xdm_pairs(agti,agtj)%c10ave

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

        ! write(*,*) r, r6, r8, r10
        ! write(*,*) pa,pb,pe

        ! DEBUG
        ! write(*,*) fd6,fd8,fd10

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij/r

            geo%nb_ene  = geo%nb_ene + pa*upe - fd6*c6/r6 - fd8*c8/r8 - fd10*c10/r10
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene  = geo%ele14_ene + inv_scee*scale2*crgij/r

            geo%nb14_ene  = geo%nb14_ene + (pa*upe - fd6*c6/r6 - fd8*c8/r8 - fd10*c10/r10) * inv_scnb
        end if
    end do

end subroutine ffdev_energy_nb_EXP_XDMTT

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_XDMTT

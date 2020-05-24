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

module ffdev_energy

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_energy_all
! ==============================================================================

subroutine ffdev_energy_all(top,geo,skipnb)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_timers

    use ffdev_nbmode_LJ

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

    call ffdev_timers_start_timer(FFDEV_POT_ENERGY_TIMER)

    ! reset energy
    geo%bond_ene = 0.0d0
    geo%angle_ene = 0.0d0
    geo%dih_ene = 0.0d0
    geo%impropr_ene = 0.0d0
    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0
    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0
    geo%total_ene = 0.0d0
    geo%rst_energy = 0.0d0

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_energy_bonds(top,geo)
        call ffdev_energy_angles(top,geo)
        call ffdev_energy_dihedrals(top,geo)
        call ffdev_energy_impropers(top,geo)
    end if

    if( calcnb ) then
        call ffdev_timers_start_timer(FFDEV_POT_NB_ENERGY_TIMER)
        select case(nb_mode)
            case(NB_VDW_LJ)
                if( (geo%sup_chrg_loaded .eqv. .true.) .and. (ele_qsource .eq. NB_ELE_QGEO) ) then
                    call ffdev_energy_nb_LJ_qgeo(top,geo)
                else
                    call ffdev_energy_nb_LJ_qtop(top,geo)
                end if

            case(NB_VDW_12_6)
                call ffdev_energy_nb_12_6(top,geo)

            case(NB_VDW_12_XDMBJ)
                call ffdev_energy_nb_12_XDMBJ(top,geo)

            case(NB_VDW_EXP_XDMBJ)
                call ffdev_energy_nb_EXP_XDMBJ(top,geo)

            case(NB_VDW_12_D3BJ)
                call ffdev_energy_nb_12_D3BJ(top,geo)

            case(NB_VDW_EXP_TTXDM)
                call ffdev_energy_nb_EXP_TTXDM(top,geo)

            case(NB_VDW_EXP_TTD3)
                call ffdev_energy_nb_EXP_TTD3(top,geo)

            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_energy_all!')
        end select
        call ffdev_timers_stop_timer(FFDEV_POT_NB_ENERGY_TIMER)
    end if

    geo%total_ene = geo%bond_ene + geo%angle_ene + geo%dih_ene &
                  + geo%impropr_ene + geo%ele14_ene + geo%nb14_ene &
                  + geo%ele_ene + geo%nb_ene

    call ffdev_timers_stop_timer(FFDEV_POT_ENERGY_TIMER)

end subroutine ffdev_energy_all

! ==============================================================================
! subroutine ffdev_energy_sapt0
! ==============================================================================

subroutine ffdev_energy_sapt0(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_timers

    use ffdev_nbmode_LJ

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    ! reset energy
    geo%sapt0_ele = 0.0
    geo%sapt0_rep = 0.0
    geo%sapt0_disp = 0.0
    geo%sapt0_total = 0.0

    if( top%sapt0_size .le. 0 ) return ! no SAPT0 list

    call ffdev_timers_start_timer(FFDEV_POT_ENERGY_TIMER)
    call ffdev_timers_start_timer(FFDEV_POT_NB_ENERGY_TIMER)

    select case(nb_mode)
        case(NB_VDW_LJ)
            call ffdev_energy_sapt0_LJ(top,geo)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_energy_sapt0!')
    end select

    geo%sapt0_total = geo%sapt0_ele + geo%sapt0_rep + geo%sapt0_disp

    call ffdev_timers_stop_timer(FFDEV_POT_NB_ENERGY_TIMER)
    call ffdev_timers_stop_timer(FFDEV_POT_ENERGY_TIMER)

end subroutine ffdev_energy_sapt0

!===============================================================================
! subroutine ffdev_energy_bonds
!===============================================================================

subroutine ffdev_energy_bonds(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,ib,ic
    real(DEVDP)     ::  b,db
    real(DEVDP)     ::  rij(3)
    ! --------------------------------------------------------------------------

    geo%bond_ene = 0.0d0

    do ib=1,top%nbonds
        ! for each bond
        i  = top%bonds(ib)%ai
        j  = top%bonds(ib)%aj
        ic = top%bonds(ib)%bt
        ! calculate rij
        rij(:) = geo%crd(:,j) - geo%crd(:,i)

        ! calculate b and db, update energy
        b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
        db = b - top%bond_types(ic)%d0
        geo%bond_ene = geo%bond_ene + 0.5*top%bond_types(ic)%k*db**2
    end do

end subroutine ffdev_energy_bonds

!===============================================================================
! subroutine ffdev_energy_angles
!===============================================================================

subroutine ffdev_energy_angles(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,j,k,ia,ic
    real(DEVDP)     :: bjiinv, bjkinv, bji2inv, bjk2inv
    real(DEVDP)     :: scp,angv,da
    real(DEVDP)     :: rji(3),rjk(3)
    ! --------------------------------------------------------------------------

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

        ! calculate da and dv
        da = angv - top%angle_types(ic)%a0
        geo%angle_ene = geo%angle_ene + 0.5d0*top%angle_types(ic)%k*da**2
    end do

end subroutine ffdev_energy_angles

!===============================================================================
! subroutine ffdev_energy_dihedrals
!===============================================================================

subroutine ffdev_energy_dihedrals(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip,pn
    real(DEVDP)     ::  scp,phi,arg
    real(DEVDP)     ::  bjinv, bkinv, bj2inv, bk2inv
    real(DEVDP)     ::  rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
    ! -----------------------------------------------------------------------------

    geo%dih_ene = 0.0d0

    do ip=1,top%ndihedrals
        i  = top%dihedrals(ip)%ai
        j  = top%dihedrals(ip)%aj
        k  = top%dihedrals(ip)%ak
        l  = top%dihedrals(ip)%al
        ic = top%dihedrals(ip)%dt

        rji(:) = geo%crd(:,i) - geo%crd(:,j)
        rjk(:) = geo%crd(:,k) - geo%crd(:,j)
        rkl(:) = geo%crd(:,l) - geo%crd(:,k)
        rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
        rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
        rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
        rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
        rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
        rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

        bj2inv = 1.0d0/(rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
        bk2inv = 1.0d0/(rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

        ! calculate scp and phi
        scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
        if ( scp .gt.  1.0 ) then
                scp =  1.0
                phi = acos (1.0) ! const
        else if ( scp .lt. -1.0 ) then
                scp = -1.0
                phi = acos (-1.0) ! const
        else
            phi = acos ( scp )
        end if
        if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
           +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
           +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) .lt. 0) then
                    phi = -phi
        end if

        select case(top%dihedral_types(ic)%mode)
            case(DIH_COS)
                do pn=1,top%dihedral_types(ic)%n
                   ! write(*,*) i,j,k,l,pn
                    if( .not. top%dihedral_types(ic)%enabled(pn) ) cycle
                    arg = pn*phi - top%dihedral_types(ic)%g(pn)
                    if( dih_cos_only ) then
                        geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%v(pn)*cos(arg)
                    else
                        geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%v(pn)*(1.0d0+cos(arg))
                    end if
                end do
            case(DIH_GRBF)
                do pn=1,top%dihedral_types(ic)%n
                    if( .not. top%dihedral_types(ic)%enabled(pn) ) cycle
                    geo%dih_ene = geo%dih_ene + top%dihedral_types(ic)%c(pn) &
                                  * exp(-(ffdev_geometry_get_dihedral_deviation(phi,top%dihedral_types(ic)%p(pn))**2) &
                                  / top%dihedral_types(ic)%w2(pn))
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented [ffdev_energy_dihedrals]!')
        end select
    end do

end subroutine ffdev_energy_dihedrals

!===============================================================================
! subroutine ffdev_energy_impropers
!===============================================================================

subroutine ffdev_energy_impropers(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         ::  i,j,k,l,ic,ip
    real(DEVDP)     ::  scp,phi,arg
    real(DEVDP)     ::  bjinv, bkinv, bj2inv, bk2inv
    real(DEVDP)     ::  rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
    ! -----------------------------------------------------------------------------

    geo%impropr_ene = 0.0d0

    do ip=1,top%nimpropers
        i  = top%impropers(ip)%ai
        j  = top%impropers(ip)%aj
        k  = top%impropers(ip)%ak
        l  = top%impropers(ip)%al
        ic = top%impropers(ip)%dt

        rji(:) = geo%crd(:,i) - geo%crd(:,j)
        rjk(:) = geo%crd(:,k) - geo%crd(:,j)
        rkl(:) = geo%crd(:,l) - geo%crd(:,k)

        rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
        rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
        rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
        rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
        rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
        rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

        bj2inv  = 1.0d0/( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
        bk2inv  = 1.0d0/( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
        bjinv = sqrt(bj2inv)
        bkinv = sqrt(bk2inv)

        scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
        if ( scp .gt.  1.0d0 ) then
            scp =  1.0d0
        else if ( scp .lt. -1.0d0 ) then
            scp = -1.0d0
        end if
        phi = acos ( scp )
        if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
            +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
            +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) &
            .lt. 0) then
                  phi = -phi
        end if

        arg = 2*phi - top%improper_types(ic)%g
        geo%impropr_ene = geo%impropr_ene + top%improper_types(ic)%v * (1.0d0 + cos(arg))

    end do

end subroutine ffdev_energy_impropers

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
    integer         :: ip,i,j,nbt
    real(DEVDP)     :: inv_scee,inv_scnb,pa,c6,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a,scale2
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
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt
        pa  = exp(top%nb_types(nbt)%pa)
        c6  = top%nb_types(nbt)%c6 * disp_fa

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

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij*ra
            geo%nb_ene  = geo%nb_ene + pa*r6a*r6a - c6*r6a
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene  = geo%ele14_ene + inv_scee*scale2*crgij*ra
            geo%nb14_ene  = geo%nb14_ene + inv_scnb*(pa*r6a*r6a - c6*r6a)
        end if
    end do

end subroutine ffdev_energy_nb_12_6

!===============================================================================
! subroutine ffdev_energy_nb_12_XDMBJ
!===============================================================================

subroutine ffdev_energy_nb_12_XDMBJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: inv_scee,inv_scnb,pa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r2a,ra,r6,r8,r10,r12,scale2,c6,c8,c10,rc,rc2,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,r6i,r8i,r10i
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0
    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_12_XDMBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa)

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! XDM
        c6  = xdm_pairs(agti,agtj)%c6ave
        c8  = xdm_pairs(agti,agtj)%c8ave
        c10 = xdm_pairs(agti,agtj)%c10ave

        rc  = top%nb_types(nbt)%r0 ! disp_fa*xdm_pairs(agti,agtj)%rc + disp_fb

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
        r2a = 1.0d0/r2
        ra  = sqrt(r2a)

        rc2 = rc*rc

        r6 = r2*r2*r2
        rc6 = rc2*rc2*rc2

        r8 = r6*r2
        rc8 = rc6*rc2

        r10 = r8*r2
        rc10 = rc8*rc2

        r12 = r6*r6

        V_aa = pa/r12

        r6i = 1.0d0/(r6+rc6)
        r8i = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij*ra
            geo%nb_ene  = geo%nb_ene + V_aa - c6*r6i - c8*r8i - c10*r10i
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene  = geo%ele14_ene + inv_scee*scale2*crgij*ra
            geo%nb14_ene  = geo%nb14_ene + inv_scnb*(V_aa - c6*r6i - c8*r8i - c10*r10i)
        end if
    end do

end subroutine ffdev_energy_nb_12_XDMBJ

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
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0
    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. xdm_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_12_XDMBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = disp_fc * exp(top%nb_types(nbt)%pa * top%nb_types(nbt)%pb)
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

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij/r
            geo%nb_ene  = geo%nb_ene + V_aa + V_bb
            geo%nb_rep =  geo%nb_rep + V_aa
            geo%nb_disp =  geo%nb_disp + V_bb
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene  = geo%ele14_ene + inv_scee*scale2*crgij/r
            geo%nb14_ene  = geo%nb14_ene + inv_scnb*(V_aa + V_bb)

            geo%nb_rep =  geo%nb_rep + inv_scnb*V_aa
            geo%nb_disp =  geo%nb_disp + inv_scnb*V_bb
        end if
    end do

end subroutine ffdev_energy_nb_EXP_XDMBJ

!===============================================================================
! subroutine ffdev_energy_nb_12_D3BJ
!===============================================================================

subroutine ffdev_energy_nb_12_D3BJ(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_mmd3_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj
    real(DEVDP)     :: inv_scee,inv_scnb,pa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r2a,ra,r6,r8,r12,scale2,c6,c8,rc,rc2,rc6,rc8
    real(DEVDP)     :: V_aa,r6i,r8i
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0
    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. mmd3_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_12_XDMBJ!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa)

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! MMD3
        c6  = mmd3_pairs(agti,agtj)%c6ave
        c8  = disp_fc*mmd3_pairs(agti,agtj)%c8ave

        rc  = disp_fa*mmd3_pairs(agti,agtj)%rc + disp_fb

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
        r2a = 1.0d0/r2
        ra  = sqrt(r2a)

        rc2 = rc*rc

        r6 = r2*r2*r2
        rc6 = rc2*rc2*rc2

        r8 = r6*r2
        rc8 = rc6*rc2

        r12 = r6*r6

        V_aa = pa/r12

        r6i = 1.0d0/(r6+rc6)
        r8i = 1.0d0/(r8+rc8)

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij*ra
            geo%nb_ene  = geo%nb_ene + V_aa - c6*r6i - c8*r8i
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene  = geo%ele14_ene + inv_scee*scale2*crgij*ra
            geo%nb14_ene  = geo%nb14_ene + inv_scnb*(V_aa - c6*r6i - c8*r8i)
        end if
    end do

end subroutine ffdev_energy_nb_12_D3BJ

!===============================================================================
! subroutine ffdev_energy_nb_EXP_TTXDM
!===============================================================================

subroutine ffdev_energy_nb_EXP_TTXDM(top,geo)

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
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_EXP_TTXDM!')
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
        c6  = xdm_pairs(agti,agtj)%c6ave
        c8  = xdm_pairs(agti,agtj)%c8ave
        c10 = xdm_pairs(agti,agtj)%c10ave

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

        arg = disp_fa*pb*r
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

end subroutine ffdev_energy_nb_EXP_TTXDM

!===============================================================================
! subroutine ffdev_energy_nb_EXP_TTD3
!===============================================================================

subroutine ffdev_energy_nb_EXP_TTD3(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_mmd3_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, agti, agtj, k
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,scale2,c6,c8,fd6,fd8,pe,arg,sump
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%nb14_ene = 0.0d0
    geo%ele_ene = 0.0d0
    geo%nb_ene = 0.0d0
    geo%nb_rep = 0.0d0
    geo%nb_disp = 0.0d0

    if( .not. mmd3_data_loaded ) then
        call ffdev_utils_exit(DEV_ERR,1,'XDM not loaded for ffdev_energy_nb_EXP_TTD3!')
    end if

    scale2 = ele_qscale*ele_qscale*332.05221729d0

    do ip=1,top%nb_size
        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        nbt = top%nb_list(ip)%nbt

        pa  = exp(top%nb_types(nbt)%pa) ! * top%nb_types(nbt)%pb)
        pb  = top%nb_types(nbt)%pb

        agti = top%atom_types(top%atoms(i)%typeid)%glbtypeid
        agtj = top%atom_types(top%atoms(j)%typeid)%glbtypeid

        ! DEBUG
        ! write(*,*) agti, trim(top%atom_types(top%atoms(i)%typeid)%name), agtj, trim(top%atom_types(top%atoms(j)%typeid)%name)

        ! MMD3
        c6  = mmd3_pairs(agti,agtj)%c6ave
        c8  = mmd3_pairs(agti,agtj)%c8ave

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

       ! upe = exp(-pb*r)

        arg = pb*r
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

        if( top%nb_list(ip)%dt .eq. 0 ) then
            geo%ele_ene  = geo%ele_ene + scale2*crgij/r

            geo%nb_ene  = geo%nb_ene + pa/(r6*r6)- fd6*c6/r6 - fd8*c8/r8
        else
            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb

            geo%ele14_ene  = geo%ele14_ene + inv_scee*scale2*crgij/r

            geo%nb14_ene  = geo%nb14_ene + (pa/(r6*r6) - fd6*c6/r6 - fd8*c8/r8) * inv_scnb
        end if
    end do

end subroutine ffdev_energy_nb_EXP_TTD3

! ------------------------------------------------------------------------------

end module ffdev_energy

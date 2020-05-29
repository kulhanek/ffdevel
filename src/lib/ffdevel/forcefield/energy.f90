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
    use ffdev_nbmode_12_DISPBJ
    use ffdev_nbmode_EXP_DISPTT

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

    ! reset energy
    geo%bond_ene = 0.0d0
    geo%angle_ene = 0.0d0
    geo%dih_ene = 0.0d0
    geo%impropr_ene = 0.0d0

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    geo%total_ene = 0.0d0
    geo%rst_energy = 0.0d0

    call ffdev_timers_start_timer(FFDEV_POT_ENERGY_TIMER)

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_energy_bonds(top,geo)
        call ffdev_energy_angles(top,geo)
        call ffdev_energy_dihedrals(top,geo)
        call ffdev_energy_impropers(top,geo)
    end if

    if( calcnb ) then
        call ffdev_timers_start_timer(FFDEV_POT_NB_ENERGY_TIMER)

        if( top%nb_params_update ) then
            call ffdev_topology_update_nb_params(top)
        end if

        select case(nb_mode)
            case(NB_VDW_LJ)
                call ffdev_energy_nb_LJ(top,geo)

            case(NB_VDW_12_DISPBJ)
                call ffdev_energy_nb_12_DISPBJ(top,geo)

            case(NB_VDW_EXP_DISPTT)
                call ffdev_energy_nb_EXP_DISPTT(top,geo)

            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_energy_sapt! (' // &
                                ffdev_topology_nb_mode_to_string(nb_mode) // ')')
        end select

        call ffdev_timers_stop_timer(FFDEV_POT_NB_ENERGY_TIMER)
    end if

    geo%total_ene = geo%bond_ene + geo%angle_ene + geo%dih_ene + geo%impropr_ene &
                  + geo%ele14_ene + geo%rep14_ene + geo%dis14_ene  &
                  + geo%ele_ene + geo%rep_ene + geo%dis_ene

    call ffdev_timers_stop_timer(FFDEV_POT_ENERGY_TIMER)

end subroutine ffdev_energy_all

! ==============================================================================
! subroutine ffdev_energy_sapt
! ==============================================================================

subroutine ffdev_energy_sapt(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils
    use ffdev_timers

    use ffdev_nbmode_LJ
    use ffdev_nbmode_12_DISPBJ
    use ffdev_nbmode_EXP_DISPTT

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    ! reset energy
    geo%sapt_ele = 0.0
    geo%sapt_rep = 0.0
    geo%sapt_dis = 0.0
    geo%sapt_total = 0.0

    if( top%sapt_size .le. 0 ) return ! no SAPT list

    call ffdev_timers_start_timer(FFDEV_POT_ENERGY_TIMER)
    call ffdev_timers_start_timer(FFDEV_POT_NB_ENERGY_TIMER)

    if( top%nb_params_update ) then
        call ffdev_topology_update_nb_params(top)
    end if

    select case(nb_mode)
        case(NB_VDW_LJ)
            call ffdev_energy_sapt_LJ(top,geo)

        case(NB_VDW_12_DISPBJ)
            call ffdev_energy_sapt_12_DISPBJ(top,geo)

        case(NB_VDW_EXP_DISPTT)
            call ffdev_energy_sapt_EXP_DISPTT(top,geo)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_energy_sapt! (' // &
                                            ffdev_topology_nb_mode_to_string(nb_mode) // ')')
    end select

    geo%sapt_total = geo%sapt_ele + geo%sapt_rep + geo%sapt_dis

    call ffdev_timers_stop_timer(FFDEV_POT_NB_ENERGY_TIMER)
    call ffdev_timers_stop_timer(FFDEV_POT_ENERGY_TIMER)

end subroutine ffdev_energy_sapt

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

! ------------------------------------------------------------------------------

end module ffdev_energy

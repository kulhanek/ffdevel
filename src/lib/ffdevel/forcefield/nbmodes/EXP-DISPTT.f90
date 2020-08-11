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

use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_EXP_DISPTT
!===============================================================================

subroutine ffdev_energy_nb_EXP_DISPTT(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer                 :: ip,i,j,dt
    real(DEVDP)             :: r2,r,dxa1,dxa2,dxa3,inv_scee,inv_scnb
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    geo%ele_ene = 0.0d0
    geo%pen_ene = 0.0d0
    geo%ind_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    do ip=1,top%nb_size
        i    = top%nb_list(ip)%ai
        j    = top%nb_list(ip)%aj
        dt   = top%nb_list(ip)%dt

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2  = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r   = sqrt(r2)

        call ffdev_energy_nbpair_EXP_DISPTT(top%nb_list(ip),r,nbene)

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + nbene%ele_ene
            geo%pen_ene = geo%pen_ene + nbene%pen_ene
            geo%ind_ene = geo%ind_ene + nbene%ind_ene
            geo%rep_ene = geo%rep_ene + nbene%rep_ene
            geo%dis_ene = geo%dis_ene + nbene%dis_ene
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * nbene%ele_ene
            geo%rep14_ene = geo%rep14_ene + inv_scnb * nbene%rep_ene
            geo%dis14_ene = geo%dis14_ene + inv_scnb * nbene%dis_ene
        end if
    end do

end subroutine ffdev_energy_nb_EXP_DISPTT

!===============================================================================
! subroutine ffdev_energy_sapt_EXP_DISPTT
!===============================================================================

subroutine ffdev_energy_sapt_EXP_DISPTT(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer                 :: ip,i,j
    real(DEVDP)             :: r2,r,dxa1,dxa2,dxa3
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    geo%sapt_ele = 0.0d0
    geo%sapt_pen = 0.0d0
    geo%sapt_ind = 0.0d0
    geo%sapt_rep = 0.0d0
    geo%sapt_dis = 0.0d0

    do ip=1,top%sapt_size
        i    = top%sapt_list(ip)%ai
        j    = top%sapt_list(ip)%aj

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)

        call ffdev_energy_nbpair_EXP_DISPTT(top%sapt_list(ip),r,nbene)

        geo%sapt_ele = geo%sapt_ele + nbene%ele_ene
        geo%sapt_pen = geo%sapt_pen + nbene%pen_ene
        geo%sapt_ind = geo%sapt_ind + nbene%ind_ene
        geo%sapt_rep = geo%sapt_rep + nbene%rep_ene
        geo%sapt_dis = geo%sapt_dis + nbene%dis_ene
    end do

end subroutine ffdev_energy_sapt_EXP_DISPTT

!===============================================================================
! subroutine ffdev_energy_nbpair_EXP_DISPTT
!===============================================================================

subroutine ffdev_energy_nbpair_EXP_DISPTT(nbpair,r,nbene)

    use ffdev_topology_dat
    use ffdev_utils
    use ffdev_nbmode_INTEGRAL

    implicit none
    type(NB_PAIR)           :: nbpair
    real(DEVDP)             :: r
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------
    integer         :: k
    real(DEVDP)     :: r2,r2a,ra
    real(DEVDP)     :: suma,sump,pe,c6,c8,c10,r6a,r8a,r10a,fd10,fd8,fd6
    real(DEVDP)     :: z1,z2,q1,q2,pa1,pb1,pa2,pb2
    real(DEVDP)     :: V_ee,V_pe,V_in,V_aa,V_bb
    real(DEVDP)     :: exc_ij,pre_ij,Sdo,Swo,sdv,tt,dtt
    real(DEVDP)     :: pfa1,pfa2,dpfa1,dpfa2,pee,pfa,dpfa
    ! --------------------------------------------------------------------------

    nbene%ele_ene = 0.0
    nbene%pen_ene = 0.0
    nbene%ind_ene = 0.0
    nbene%rep_ene = 0.0
    nbene%dis_ene = 0.0
    nbene%tot_ene = 0.0

    V_ee = 0.0d0
    V_pe = 0.0d0
    V_aa = 0.0d0
    V_in = 0.0d0
    V_bb = 0.0d0

    ra   = 1.0d0/r
    r2   = r**2
    r2a  = ra**2

    z1  = nbpair%z1
    q1  = nbpair%q1
    pa1 = nbpair%pa1
    pb1 = nbpair%pb1

    z2  = nbpair%z2
    q2  = nbpair%q2
    pa2 = nbpair%pa2
    pb2 = nbpair%pb2

! electrostatics
    V_ee = q1*q2*ra

! exchange repulsion
    select case(exp_mode)
        case(EXP_MODE_DO)
            call ffdev_nbmode_INTEGRAL_do(r,pb1,pb2,Sdo,sdv,tt,dtt)
            exc_ij = Sdo
    ! --------------------------
        case(EXP_MODE_WO)
            call ffdev_nbmode_INTEGRAL_wo(r,pb1,pb2,Swo,sdv,tt,dtt)
            exc_ij =   (Swo**2) * ra
            tt     = 2.0d0*tt + 1.0d0
    ! --------------------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_mode in ffdev_energy_nbpair_EXP_DISPTT!')
    end select

    ! prefactor
    select case(exp_pa_mode)
        case(EXP_PA_FREEOPT)
            pre_ij = pa1*pa2
            V_aa   = pre_ij*exc_ij
    ! --------------------------
        case(EXP_PA_CHARGES)
            pre_ij = exp(k_exc)*(z1-q1)*(z2-q2)
            V_aa   = pre_ij*exc_ij
    ! --------------------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_pa_mode in ffdev_energy_nbpair_EXP_DISPTT!')
    end select

! induction part
    V_in = 0.0d0
    if( ind_enabled ) then
        select case(ind_mode)
            case(IND_MODE_K2EXC)
                V_in = - V_aa * exp(k_ind)
        ! --------------------------
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented ind_mode in ffdev_energy_nbpair_EXP_DISPTT!')
        end select
    end if

! penetration energy
    V_pe  = 0.0d0
    if( pen_enabled ) then
        select case(pen_mode)
            case(PEN_MODE_SCI)
                call ffdev_nbmode_INTEGRAL_ci_ez_damp(r,pb1*damp_pe,pfa1,dpfa1)
                call ffdev_nbmode_INTEGRAL_ci_ez_damp(r,pb2*damp_pe,pfa2,dpfa2)
                call ffdev_nbmode_INTEGRAL_ci_ee_damp(r,pb1*damp_pe,pb2*damp_pe,pfa,dpfa)
                pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 - (z1-q1)*(z2-q2)*pfa
                V_pe = pee * ra
        ! --------------------------
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_mode in ffdev_energy_nbpair_EXP_DISPTT!')
        end select
    end if

! dispersion
    ! dispersion coefficients
    c6   = nbpair%c6
    c8   = nbpair%c8
    c10  = nbpair%c10

    tt   = damp_tb * tt
    pe   = exp(-tt)

! 6
    r6a  = r2a*r2a*r2a
    sump = 1.0d0
    suma = 1.0d0
    do k=1,6
        sump = sump * tt / real(k,DEVDP)
        suma = suma + sump
    end do
    fd6  = 1.0d0 - pe*suma

! 8
    r8a  = r6a*r2a
    sump = sump * tt / real(7,DEVDP)
    suma = suma + sump

    sump = sump * tt / real(8,DEVDP)
    suma = suma + sump

    fd8  = 1.0d0 - pe*suma

! 10
    r10a = r8a*r2a
    sump = sump * tt / real(9,DEVDP)
    suma = suma + sump

    sump = sump * tt / real(10,DEVDP)
    suma = suma + sump

    fd10 = 1.0d0 - pe*suma

    V_bb = - fd6*c6*r6a - fd8*c8*r8a - fd10*c10*r10a

    nbene%ele_ene = V_ee
    nbene%pen_ene = V_pe
    nbene%ind_ene = V_in
    nbene%rep_ene = V_aa
    nbene%dis_ene = V_bb
    nbene%tot_ene = V_ee + V_pe + V_in + V_aa + V_bb

end subroutine ffdev_energy_nbpair_EXP_DISPTT

!===============================================================================
! subroutine ffdev_gradient_nb_EXP_DISPTT
!===============================================================================

subroutine ffdev_gradient_nb_EXP_DISPTT(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils
    use ffdev_nbmode_INTEGRAL

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,k,dt
    real(DEVDP)     :: inv_scee,inv_scnb,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,ra,r6a,r8a,r10a,c6,c8,c10,fd6,fd8,fd10,pe,sump
    real(DEVDP)     :: V_b6,V_b8,V_b10,r2a,dfd6,dfd8,dfd10,sumd,suma
    real(DEVDP)     :: dvee,dvpe,dvin,dvaa,dvbb,dva
    real(DEVDP)     :: z1,z2,q1,q2,pa1,pb1,pa2,pb2
    real(DEVDP)     :: V_ee,V_pe,V_in,V_aa,V_bb
    real(DEVDP)     :: exc_ij,pre_ij,Sdo,Swo,sdv,tt,dtt
    real(DEVDP)     :: pfa1,pfa2,dpfa1,dpfa2,pee,pfa,dpfa
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%ind_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    V_ee = 0.0d0
    V_pe = 0.0d0
    V_aa = 0.0d0
    V_in = 0.0d0
    V_bb = 0.0d0

    dvee  = 0.0d0
    dvpe  = 0.0d0
    dvin  = 0.0d0
    dvaa  = 0.0d0
    dvbb  = 0.0d0

    do ip=1,top%nb_size
        i    = top%nb_list(ip)%ai
        j    = top%nb_list(ip)%aj
        dt   = top%nb_list(ip)%dt

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)
        ra   = 1.0d0 / r
        r2a  = 1.0d0 / r2

        z1  = top%nb_list(ip)%z1
        q1  = top%nb_list(ip)%q1
        pa1 = top%nb_list(ip)%pa1
        pb1 = top%nb_list(ip)%pb1

        z2  = top%nb_list(ip)%z2
        q2  = top%nb_list(ip)%q2
        pa2 = top%nb_list(ip)%pa2
        pb2 = top%nb_list(ip)%pb2

    ! electrostatics
        V_ee  =   q1*q2 * ra
        dvee  = - q1*q2 * r2a

    ! exchange repulsion
        select case(exp_mode)
            case(EXP_MODE_DO)
                call ffdev_nbmode_INTEGRAL_do(r,pb1,pb2,Sdo,sdv,tt,dtt)
                exc_ij = Sdo
        ! --------------------------
            case(EXP_MODE_WO)
                call ffdev_nbmode_INTEGRAL_wo(r,pb1,pb2,Swo,sdv,tt,dtt)
                exc_ij =   (Swo**2) * ra
                sdv    = - (Swo**2) * r2a + 2.0d0*Swo*sdv * ra
                tt     = 2.0d0*tt + 1.0d0
                dtt    = 2.0d0*dtt
        ! --------------------------
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_mode in ffdev_energy_nbpair_EXP_DISPTT!')
        end select

        ! prefactor
        select case(exp_pa_mode)
            case(EXP_PA_FREEOPT)
                pre_ij = pa1*pa2
                V_aa   = pre_ij*exc_ij
                dvaa   = pre_ij*sdv
        ! --------------------------
            case(EXP_PA_CHARGES)
                pre_ij = exp(k_exc)*(z1-q1)*(z2-q2)
                V_aa   = pre_ij*exc_ij
                dvaa   = pre_ij*sdv
        ! --------------------------
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_pa_mode in ffdev_energy_nbpair_EXP_DISPTT!')
        end select

    ! induction part
        V_in = 0.0d0
        dvin = 0.0d0
        if( ind_enabled ) then
            select case(ind_mode)
                case(IND_MODE_K2EXC)
                    V_in = - V_aa * exp(k_ind)
                    dvin = - dvaa * exp(k_ind)
            ! --------------------------
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented ind_mode in ffdev_energy_nbpair_EXP_DISPTT!')
            end select
        end if

    ! penetration energy
        V_pe = 0.0d0
        dvpe = 0.0d0
        if( pen_enabled ) then
            select case(pen_mode)
                case(PEN_MODE_SCI)
                    call ffdev_nbmode_INTEGRAL_ci_ez_damp(r,pb1*damp_pe,pfa1,dpfa1)
                    call ffdev_nbmode_INTEGRAL_ci_ez_damp(r,pb2*damp_pe,pfa2,dpfa2)
                    call ffdev_nbmode_INTEGRAL_ci_ee_damp(r,pb1*damp_pe,pb2*damp_pe,pfa,dpfa)
                    pee  = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 - (z1-q1)*(z2-q2)*pfa
                    V_pe = pee * ra
                    dvpe = (- V_pe + z1*(z2-q2)*dpfa2 + (z1-q1)*z2*dpfa1 - (z1-q1)*(z2-q2)*dpfa) * ra
            ! --------------------------
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_mode in ffdev_energy_nbpair_EXP_DISPTT!')
            end select
        end if

    ! dispersion
        ! dispersion coefficients
        c6   = top%nb_list(ip)%c6
        c8   = top%nb_list(ip)%c8
        c10  = top%nb_list(ip)%c10

        tt   = damp_tb * tt
        dtt  = damp_tb * dtt

        pe   = exp(-tt)

    ! 6
        r6a  = r2a*r2a*r2a
        fd6  = 1.0d0
        dfd6 = 0.0d0
        sump = 1.0d0
        suma = 1.0d0
        sumd = 0.0d0
        do k=1,6
            sumd = sumd + sump ! sump is one item delayed in derivatives
            sump = sump * tt / real(k,DEVDP)
            suma = suma + sump
        end do
        fd6   = 1.0d0 - pe*suma
        dfd6  = (pe*suma - pe*sumd)

    ! 8
        r8a   = r6a*r2a
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * tt / real(7,DEVDP)
        suma  = suma  + sump

        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * tt / real(8,DEVDP)
        suma  = suma  + sump

        fd8   = 1.0d0 - pe*suma
        dfd8  = (pe*suma - pe*sumd)

    ! 10
        r10a  = r8a*r2a
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * tt / real(9,DEVDP)
        suma  = suma  + sump

        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * tt / real(10,DEVDP)
        suma  = suma  + sump

        fd10  = 1.0d0 - pe*suma
        dfd10 = (pe*suma - pe*sumd)

        V_b6  = - fd6*c6*r6a
        V_b8  = - fd8*c8*r8a
        V_b10 = - fd10*c10*r10a
        V_bb  = V_b6 + V_b8 + V_b10

        dvbb  = - (6.0d0*V_b6 + 8.0d0*V_b8 + 10.0d0*V_b10)*ra  &
                - (c6*r6a*dfd6 + c8*r8a*dfd8 + c10*r10a*dfd10)*dtt

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%pen_ene = geo%pen_ene + V_pe
            geo%ind_ene = geo%ind_ene + V_in
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb

            dva = ra*( dvee + dvpe + dvin + dvaa + dvbb )
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb

            ! FIXME - no penetration energy?
            dva = ra*( inv_scee*dvee + inv_scnb*(dvaa + dvbb) )
        end if

        ! calculate gradient
        dxa1 = dva*dxa1
        dxa2 = dva*dxa2
        dxa3 = dva*dxa3

        geo%grd(1,i) = geo%grd(1,i) + dxa1
        geo%grd(2,i) = geo%grd(2,i) + dxa2
        geo%grd(3,i) = geo%grd(3,i) + dxa3

        geo%grd(1,j) = geo%grd(1,j) - dxa1
        geo%grd(2,j) = geo%grd(2,j) - dxa2
        geo%grd(3,j) = geo%grd(3,j) - dxa3

    end do

end subroutine ffdev_gradient_nb_EXP_DISPTT

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_DISPTT

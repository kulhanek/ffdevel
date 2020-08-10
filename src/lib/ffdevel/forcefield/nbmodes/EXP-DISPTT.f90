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

    implicit none
    type(NB_PAIR)           :: nbpair
    real(DEVDP)             :: r
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------
    integer         :: k
    real(DEVDP)     :: pa,pb,pc,tb,er,pr,lvaa
    real(DEVDP)     :: r2,r2a,suma,sump,upe,c6,c8,c10,pe
    real(DEVDP)     :: V_aa,V_bb,r6a,r8a,r10a,arg,fd10,fd8,fd6
    real(DEVDP)     :: z1,z2,q1,q2,pepa1,pepa2,pepb1,pepb2,pee,V_pe,crgij,V_ee
    real(DEVDP)     :: pfa1,pfb1,pfa2,pfb2,pepk,pfbb
    real(DEVDP)     :: a1,a2,b2,inva2b2,V_in,pb1,pb2,Sij,t1,pbe1,pbe2,slvaa
    ! --------------------------------------------------------------------------

    nbene%ele_ene = 0.0
    nbene%pen_ene = 0.0
    nbene%ind_ene = 0.0
    nbene%rep_ene = 0.0
    nbene%dis_ene = 0.0
    nbene%tot_ene = 0.0

    r2   = r**2
    r2a  = 1.0d0/r2

! electrostatics
    crgij = nbpair%crgij
    V_ee = crgij/r

! induction - part I
    Sij = 0.0d0
    slvaa = 0.0d0 ! part for TT damping in EXP_MODE_MEDFF
    if( calc_sij ) then
        z1  = nbpair%z1
        q1  = nbpair%q1
        pb1 = nbpair%pb1
        z2  = nbpair%z2
        q2  = nbpair%q2
        pb2 = nbpair%pb2

        ! calculate overlap integral
        if( abs(pb1 - pb2) .gt. 0.1d0 ) then
            inva2b2 = 1.0d0/(pb1**2 - pb2**2)
            a1 =  4.0d0*(pb1**4)*(pb2**4)*(inva2b2**3) + (pb1**3)*(pb2**4)*(inva2b2**2)*r
            a2 = -4.0d0*(pb2**4)*(pb1**4)*(inva2b2**3) + (pb2**3)*(pb1**4)*(inva2b2**2)*r
            pbe1 = exp(-pb1*r)
            pbe2 = exp(-pb2*r)
            t1 = a1*pbe1 + a2*pbe2
            Sij = 1.0d0/(8.0d0*DEV_PI*r) * t1

            slvaa = 1.0d0/r - (1.0d0/t1)*( &
                     - a1*pbe1*pb1 + pbe1*(pb1**3)*(pb2**4)*(inva2b2**2) &
                     - a2*pbe2*pb2 + pbe2*(pb2**3)*(pb1**4)*(inva2b2**2) )
        else
            pb = 0.5d0*(pb1 + pb2)
            Sij = pb**3 / (192.0d0 * DEV_PI) * &
                  (3.0d0 + 3.0d0*pb*r + (pb*r)**2) * exp(-pb*r)
            slvaa = pb - (3.0d0*pb + 2.0d0*pb**2*r)/(3.0d0 + 3.0d0*pb*r + (pb*r)**2)
        end if

        ! complete overlap
        Sij     = Sij * (z1-q1) * (z2-q2)
    end if

! repulsion
    pa   = nbpair%pa
    pb   = nbpair%pb
    pc   = nbpair%pc

    lvaa  = 0.0d0
    select case(exp_mode)
        case(EXP_MODE_BM)
            er   = pb*r
            upe  = exp(-er)
            V_aa = pa*upe
            lvaa = pb
    ! --------------------------
        case(EXP_MODE_DO)
            er   = pb*r
            pr   = 1.0d0 + er + er**2/3.0d0
            upe  = exp(-er)
            V_aa = pa*pr*upe
            lvaa = pb - (pb + 2.0d0*pb**2*r/3.0d0)/pr
    ! --------------------------
        case(EXP_MODE_WO)
            er   = pb*r
            a1   = 1.0d0 + er/2.0d0 + er**2/12.0d0
            pr   = a1**2/r
            upe  = exp(-er)
            V_aa = pa*pr*upe
            lvaa = pb + 1.0d0/r - 2.0d0*(pb/2.0d0 + 2.0d0*(pb**2)*r/12.0d0)/a1
    ! --------------------------
        case(EXP_MODE_SC)
            er   = pb*r
            pr   = (r*DEV_A2AU)**pc ! this part must be in a.u., see pc
            upe  = exp(-er)
            V_aa = pa*pr*upe
            lvaa = pb - pc/r
    ! --------------------------
        case(EXP_MODE_MEDFF)
            V_aa = Sij * exp(k_exc)
            lvaa = slvaa
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented exp_mode in ffdev_energy_nbpair_EXP_DISPTT!')
    end select

! induction part - II
    V_in = 0.0d0
    if( ind_enabled ) then
        select case(ind_mode)
            case(IND_MODE_MEDFF)
                V_in = - Sij * exp(k_ind)
        ! --------------------------
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
            case(PEN_MODE_AMOEBA)
                z1      = nbpair%z1
                q1      = nbpair%q1
                pepa1   = nbpair%pepa1
                pepb1   = nbpair%pepb1

                z2      = nbpair%z2
                q2      = nbpair%q2
                pepa2   = nbpair%pepa2
                pepb2   = nbpair%pepb2

                pfa1 = exp(-pepa1*r)
                pfb1 = exp(-pepb1*r)
                pfa2 = exp(-pepa2*r)
                pfb2 = exp(-pepb2*r)

                pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
                      - (z1-q1)*(z2-q2)*(pfb1 + pfb2 - pfb1*pfb2)
                V_pe = pee/r
        ! --------------------------
            case(PEN_MODE_EFP_M1)
                z1      = nbpair%z1
                q1      = nbpair%q1
                pepa1   = nbpair%pepa1

                z2      = nbpair%z2
                q2      = nbpair%q2
                pepa2   = nbpair%pepa2

                pfa1 = exp(-pepa1*r)*(1.0 + 0.5d0*pepa1*r)
                pfa2 = exp(-pepa2*r)*(1.0 + 0.5d0*pepa2*r)

                if( abs(pepa1 - pepa2) .gt. 0.1d0 ) then
                    a2 = pepa1 ** 2
                    b2 = pepa2 ** 2
                    inva2b2 = 1.0d0/(b2 - a2)
                    ! order in inva2b2 is somewhere opposite, which changes the sign
                    pfbb =  + exp(-pepa1*r) * (b2*inva2b2)**2 *                 &
                                 (1.0d0 - 2.0d0*a2*inva2b2 + 0.5d0*pepa1*r)     &
                            + exp(-pepa2*r) * (a2*inva2b2)**2 *                 &
                                 (1.0d0 + 2.0d0*b2*inva2b2 + 0.5d0*pepa2*r)
                else
                    pepk = 0.5d0*(pepa1 + pepa2)
                    pfbb =       + exp(-pepk*r) *               &
                         ( 1.0d0 + 11.0d0/16.0d0*pepk*r         &
                                 + 3.0d0/16.0d0*(pepk*r)**2     &
                                 + 1.0d0/48.0d0*(pepk*r)**3 )
                end if

                pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
                      - (z1-q1)*(z2-q2)*pfbb
                V_pe = pee/r
        ! --------------------------
            case(PEN_MODE_EFP_M2)
                z1      = nbpair%z1
                q1      = nbpair%q1
                pepa1   = nbpair%pepa1

                z2      = nbpair%z2
                q2      = nbpair%q2
                pepa2   = nbpair%pepa2

                pfa1 = exp(-pepa1*r)
                pfa2 = exp(-pepa2*r)

                if( abs(pepa1 - pepa2) .gt. 0.1d0 ) then
                    a2 = pepa1 ** 2
                    b2 = pepa2 ** 2
                    inva2b2 = 1.0d0/(b2 - a2)
                    pfbb = pfa1*b2*inva2b2 - pfa2*a2*inva2b2
                else
                    pepk = 0.5d0*(pepa1 + pepa2)
                    pfbb = exp(-pepk*r)*(1.0 + 0.5d0*pepk*r)
                end if

                pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
                      - (z1-q1)*(z2-q2)*pfbb
                V_pe = pee/r
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented pen_mode in ffdev_energy_nbpair_EXP_DISPTT!')
        end select
    end if

! dispersion
    ! dispersion coefficients
    c6   = nbpair%c6
    c8   = nbpair%c8
    c10  = nbpair%c10

    ! TT damping factor
    if( disptt_exact ) then
        tb = damp_tb * lvaa
    else
        tb = nbpair%tb
    end if

    arg  = tb*r
    pe   = exp(-arg)

! 6
    r6a  = r2a*r2a*r2a
    sump = 1.0d0
    suma = 1.0d0
    do k=1,6
        sump = sump * arg / real(k,DEVDP)
        suma = suma + sump
    end do
    fd6  = 1.0d0 - pe*suma

! 8
    r8a  = r6a*r2a
    sump = sump * arg / real(7,DEVDP)
    suma = suma + sump
    sump = sump * arg / real(8,DEVDP)
    suma = suma + sump
    fd8  = 1.0d0 - pe*suma

! 10
    r10a = r8a*r2a
    sump = sump * arg / real(9,DEVDP)
    suma = suma + sump
    sump = sump * arg / real(10,DEVDP)
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

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,k,dt
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,pc,crgij,dxa1,dxa2,dxa3,upe,tb
    real(DEVDP)     :: r2,r,r6a,r8a,r10a,c6,c8,c10,fd6,fd8,fd10,pe,arg,sump
    real(DEVDP)     :: V_aa,V_b6,V_b8,V_b10,V_ee,dva,r2a,dfd6,dfd8,dfd10,sumd,suma
    real(DEVDP)     :: dvee,dvpe,dvaa,dvbb,er,pr,V_pe,lvaa,tvaa,pepk,pfbb
    real(DEVDP)     :: z1,z2,q1,q2,pfa1,pfa2,pfb1,pfb2,pepa1,pepa2,pepb1,pepb2,pee
    real(DEVDP)     :: a2,b2,inva2b2,pfe1,pfe2,pfeb,V_in,pb1,pb2,Sij,dvin,dvs
    real(DEVDP)     :: pbe,pbe1,pbe2,a1,slvaa,t1,stvaa
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
    geo%ind_ene = 0.0d0
    geo%rep_ene = 0.0d0
    geo%dis_ene = 0.0d0

    do ip=1,top%nb_size
        i    = top%nb_list(ip)%ai
        j    = top%nb_list(ip)%aj
        dt   = top%nb_list(ip)%dt

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a  = 1.0d0 / r2
        r    = sqrt(r2)

    ! electrostatics
        crgij = top%nb_list(ip)%crgij
        V_ee  = crgij/r
        dvee  = V_ee

    ! induction
        Sij = 0.0d0
        dvs = 0.0d0
        slvaa = 0.0d0
        if( calc_sij ) then
            z1  = top%nb_list(ip)%z1
            q1  = top%nb_list(ip)%q1
            pb1 = top%nb_list(ip)%pb1
            z2  = top%nb_list(ip)%z2
            q2  = top%nb_list(ip)%q2
            pb2 = top%nb_list(ip)%pb2

            ! calculate overlap integral
            if( abs(pb1 - pb2) .gt. 0.1d0 ) then
                inva2b2 = 1.0d0/(pb1**2 - pb2**2)
                a1 =  4.0d0*(pb1**4)*(pb2**4)*(inva2b2**3) + (pb1**3)*(pb2**4)*(inva2b2**2)*r
                a2 = -4.0d0*(pb2**4)*(pb1**4)*(inva2b2**3) + (pb2**3)*(pb1**4)*(inva2b2**2)*r
                pbe1 = exp(-pb1*r)
                pbe2 = exp(-pb2*r)
                t1 = a1*pbe1 + a2*pbe2
                Sij = 1.0d0/(8.0d0*DEV_PI*r) * t1

                dvs = 1.0d0/(8.0d0*DEV_PI*r) * &
                   ( -a1*pbe1*pb1 + pbe1*(pb1**3)*(pb2**4)*(inva2b2**2) &
                     -a2*pbe2*pb2 + pbe2*(pb2**3)*(pb1**4)*(inva2b2**2) ) + &
                    - 1.0d0/(8.0d0*DEV_PI*r**2) * t1

                slvaa =  1.0d0/r - (1.0d0/t1)*( &
                         - a1*pbe1*pb1 + pbe1*(pb1**3)*(pb2**4)*(inva2b2**2) &
                         - a2*pbe2*pb2 + pbe2*(pb2**3)*(pb1**4)*(inva2b2**2) )

                stvaa = -1.0d0/r  ! FIXME
            else
                pb  = 0.5d0*(pb1 + pb2)
                pbe = exp(-pb*r)
                a1  = 3.0d0 + 3.0d0*pb*r + (pb*r)**2
                Sij = pb**3 / (192.0d0 * DEV_PI) * a1 * pbe
                dvs = pb**3 / (192.0d0 * DEV_PI) * (-a1*pbe*pb +  pbe*(3.0d0*pb + 2.0d0*(pb**2)*r))
                slvaa = pb - (3.0d0*pb + 2.0d0*(pb**2)*r)/a1

                stvaa = - 2.0d0*(pb**2)/a1 + ((3.0d0*pb + 2.0d0*(pb**2)*r)/a1)**2
            end if

            ! correct for r2a
            dvs = - dvs * r

            ! complete overlap
            Sij = Sij * (z1-q1) * (z2-q2)
            dvs = dvs * (z1-q1) * (z2-q2)
        end if

    ! repulsion
        pa   = top%nb_list(ip)%pa
        pb   = top%nb_list(ip)%pb
        pc   = top%nb_list(ip)%pc

        dvaa = 0.0d0 ! derivative
        lvaa = 0.0d0 ! TT
        tvaa = 0.0d0 ! TT derivative
        select case(exp_mode)
            case(EXP_MODE_BM)
                er   = pb*r
                upe  = exp(-er)
                V_aa = pa*upe
                lvaa = pb
                tvaa = 0.0d0
                dvaa = V_aa*pb*r
         !------------
            case(EXP_MODE_DO)
                er   = pb*r
                pr   = 1.0d0 + er + er**2/3.0d0
                upe  = exp(-er)
                V_aa = pa*pr*upe
                lvaa = pb - (pb + 2.0d0*(pb**2)*r/3.0d0)/pr
                tvaa = - 2.0d0*(pb**2)/(3.0d0*pr) + (pb + 2.0d0*(pb**2)*r/3.0d0)*(pb + 2.0d0*(pb**2)*r/3.0d0)/(pr**2)
                tvaa = tvaa * r2
                dvaa = V_aa*pb*r - r*pa*upe*(pb + 2.0d0*(pb**2)*r/3.0d0)
        !------------
            case(EXP_MODE_WO)
                er   = pb*r
                a1   = 1.0d0 + er/2.0d0 + er**2/12.0d0
                pr   = a1**2/r
                upe  = exp(-er)
                V_aa = pa*pr*upe
                lvaa = pb + 1.0d0/r - 2.0d0*(pb/2.0d0 + 2.0d0*(pb**2)*r/12.0d0)/a1
!                tvaa = - 1.0d0/r**2 - (pb**2)/(3.0d0*a1) &
!                       + 2.0d0*(pb/2.0d0 + 2.0d0*(pb**2)*r/12.0d0)*(pb/2.0d0 + 2.0d0*(pb**2)*r/12.0d0)/(a1**2)
!                tvaa = tvaa * r2
                ! optimized
                tvaa = - 1.0d0 - (er**2)/(3.0d0*a1) &
                       + 2.0d0*r2*((pb/2.0d0 + 2.0d0*(pb**2)*r/12.0d0)/a1)**2
                dvaa = V_aa*pb*r &
                     - r*pa*upe*(-r2a + 5.0d0/12.0d0*pb**2 + 2.0d0/12.0d0*pb**3*r + 3.0d0/144.0d0*pb**4 * r2)
        !------------
            case(EXP_MODE_SC)
                er   = pb*r
                pr   = (r*DEV_A2AU)**pc ! this part must be in a.u., see pc
                upe  = exp(-er)
                V_aa = pa*pr*upe
                lvaa = pb - pc/r
                tvaa = pc           ! unoptimized: r*(+pc/r**2)*r
                dvaa = V_aa*pb*r  - pa*upe*pc*pr ! unoptimized: - pa*upe*pc*r**(pc-1.0d0)*r, extra r is because of r2a below
        !------------
            case(EXP_MODE_MEDFF)
                V_aa = Sij * exp(k_exc)
                lvaa = slvaa
                tvaa = stvaa * r2
                dvaa = dvs * exp(k_exc)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not defined EXP mode in ffdev_energy_nb_EXP_DISPTT!')
        end select

    ! induction part - II
        V_in = 0.0d0
        dvin = 0.0d0
        if( ind_enabled ) then
            select case(ind_mode)
                case(IND_MODE_MEDFF)
                    V_in = - Sij * exp(k_ind)
                    dvin = - dvs * exp(k_ind)
                case(IND_MODE_K2EXC)
                    V_in = - V_aa * exp(k_ind)
                    dvin = - dvaa * exp(k_ind)
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented ind_mode in ffdev_energy_nbpair_EXP_DISPTT!')
            end select
        end if

    ! penetration energy
        V_pe  = 0.0d0
        dvpe  = 0.0d0
        if( pen_enabled ) then
            select case(pen_mode)
                case(PEN_MODE_AMOEBA)
                    z1      = top%nb_list(ip)%z1
                    q1      = top%nb_list(ip)%q1
                    pepa1   = top%nb_list(ip)%pepa1
                    pepb1   = top%nb_list(ip)%pepb1

                    z2      = top%nb_list(ip)%z2
                    q2      = top%nb_list(ip)%q2
                    pepa2   = top%nb_list(ip)%pepa2
                    pepb2   = top%nb_list(ip)%pepb2

                    pfa1 = exp(-pepa1*r)
                    pfb1 = exp(-pepb1*r)
                    pfa2 = exp(-pepa2*r)
                    pfb2 = exp(-pepb2*r)

                    pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
                          - (z1-q1)*(z2-q2)*(pfb1 + pfb2 - pfb1*pfb2)
                    V_pe = pee/r

                    ! derivatives
                    dvpe = V_pe + z1*(z2-q2)*pfa2*pepa2 + (z1-q1)*z2*pfa1*pepa1 &
                         - (z1-q1)*(z2-q2)*(pfb1*pepb1 + pfb2*pepb2 - pfb1*pfb2*pepb2 - pfb1*pepb1*pfb2)
            ! --------------------------
                case(PEN_MODE_EFP_M1)
                    z1      = top%nb_list(ip)%z1
                    q1      = top%nb_list(ip)%q1
                    pepa1   = top%nb_list(ip)%pepa1

                    z2      = top%nb_list(ip)%z2
                    q2      = top%nb_list(ip)%q2
                    pepa2   = top%nb_list(ip)%pepa2

                    pfe1 = exp(-pepa1*r)
                    pfa1 = pfe1*(1.0 + 0.5d0*pepa1*r)

                    pfe2 = exp(-pepa2*r)
                    pfa2 = pfe2*(1.0 + 0.5d0*pepa2*r)

                    if( abs(pepa1 - pepa2) .gt. 0.1d0 ) then
                        a2 = pepa1 ** 2
                        b2 = pepa2 ** 2
                        inva2b2 = 1.0d0/(b2 - a2)
                        ! order in inva2b2 is somewhere opposite, which changes the sign
                        pfbb =  + pfe1 * (b2*inva2b2)**2 * (1.0d0 - 2.0d0*a2*inva2b2 + 0.5d0*pepa1*r) &
                                + pfe2 * (a2*inva2b2)**2 * (1.0d0 + 2.0d0*b2*inva2b2 + 0.5d0*pepa2*r)

                        ! derivatives - 1/r prefactor is below in r2a including the reverse sign
                        dvpe = + z1*(z2-q2)*(pfa2*pepa2 - 0.5d0*pepa2*pfe2)     &
                               + (z1-q1)*z2*(pfa1*pepa1 - 0.5d0*pepa1*pfe1)     &
                               - (z1-q1)*(z2-q2)*(                              &
                               + pfe1 * (b2*inva2b2)**2 * (1.0d0 - 2.0d0*a2*inva2b2 + 0.5d0*pepa1*r) * pepa1 &
                               - 0.5d0*pepa1*pfe1 * (b2*inva2b2)**2  &
                               + pfe2 * (a2*inva2b2)**2 * (1.0d0 + 2.0d0*b2*inva2b2 + 0.5d0*pepa2*r) * pepa2 &
                               - 0.5d0*pepa2*pfe2 * (a2*inva2b2)**2 )
                    else
                        pepk = 0.5d0*(pepa1 + pepa2)
                        pfeb = exp(-pepk*r)
                        pfbb = + pfeb * ( 1.0d0                 &
                               + 11.0d0/16.0d0*pepk*r           &
                               + 3.0d0/16.0d0*(pepk*r)**2       &
                               + 1.0d0/48.0d0*(pepk*r)**3 )

                        ! derivatives - 1/r prefactor is below in r2a including the reverse sign
                        dvpe = + z1*(z2-q2)*(pfa2*pepa2 - 0.5d0*pepa2*pfe2)     &
                               + (z1-q1)*z2*(pfa1*pepa1 - 0.5d0*pepa1*pfe1)     &
                               - (z1-q1)*(z2-q2)*(pfbb*pepk &
                               - pfeb*(11.0d0/16.0d0*pepk + 2.0d0*3.0d0/16.0d0*r*(pepk)**2 + 3.0d0/48.0d0*r**2*(pepk)**3 ) )
                    end if

                    pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
                          - (z1-q1)*(z2-q2)*pfbb
                    V_pe = pee/r

                    ! derivatives
                    dvpe = dvpe + V_pe
            ! --------------------------
                case(PEN_MODE_EFP_M2)
                    z1      = top%nb_list(ip)%z1
                    q1      = top%nb_list(ip)%q1
                    pepa1   = top%nb_list(ip)%pepa1

                    z2      = top%nb_list(ip)%z2
                    q2      = top%nb_list(ip)%q2
                    pepa2   = top%nb_list(ip)%pepa2

                    pfa1 = exp(-pepa1*r)
                    pfa2 = exp(-pepa2*r)

                    if( abs(pepa1 - pepa2) .gt. 0.1d0 ) then
                        a2 = pepa1 ** 2
                        b2 = pepa2 ** 2
                        inva2b2 = 1.0d0/(b2 - a2)
                        pfbb = pfa1*b2*inva2b2 - pfa2*a2*inva2b2

                        ! derivatives - 1/r prefactor is below in r2a including the reverse sign
                        dvpe = + z1*(z2-q2)*pfa2*pepa2 + (z1-q1)*z2*pfa1*pepa1 &
                               - (z1-q1)*(z2-q2)*(pfa1*b2*inva2b2*pepa1 - pfa2*a2*inva2b2*pepa2)
                    else
                        pepk = 0.5d0*(pepa1 + pepa2)
                        pfeb = exp(-pepk*r)
                        pfbb = pfeb*(1.0 + 0.5d0*pepk*r)

                        ! derivatives - 1/r prefactor is below in r2a including the reverse sign
                        dvpe = + z1*(z2-q2)*pfa2*pepa2 + (z1-q1)*z2*pfa1*pepa1 &
                               - (z1-q1)*(z2-q2)*(pfbb*pepk - 0.5d0*pfeb*pepk)
                    end if

                    pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
                          - (z1-q1)*(z2-q2)*pfbb
                    V_pe = pee/r

                    ! derivatives
                    dvpe = dvpe + V_pe
            end select
        end if

    ! dispersion
        ! dispersion coefficients
        c6   = top%nb_list(ip)%c6
        c8   = top%nb_list(ip)%c8
        c10  = top%nb_list(ip)%c10

        ! TT damping factor
        if( disptt_exact ) then
            tb = damp_tb * lvaa
        else
            tb = top%nb_list(ip)%tb
            tvaa = 0.0d0
        end if

        arg  = tb*r
        pe   = exp(-arg)

    ! 6
        r6a  = r2a*r2a*r2a
        fd6  = 1.0d0
        dfd6 = 0.0d0
        sump = 1.0d0
        suma = 1.0d0
        sumd = 0.0d0
        do k=1,6
            sumd = sumd + sump ! sump is one item delayed in derivatives
            sump = sump * arg / real(k,DEVDP)
            suma = suma + sump
            dfd6 = dfd6 + real(k,DEVDP)*sump
        end do
        fd6   = 1.0d0 - pe*suma
        dfd6  = (pe*suma - pe*sumd)*(arg + tvaa)  ! extra r in arg is due to r2a factor in dva

    ! 8
        r8a   = r6a*r2a
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(7,DEVDP)
        suma  = suma  + sump
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(8,DEVDP)
        suma  = suma  + sump
        fd8   = 1.0d0 - pe*suma
        dfd8  = (pe*suma - pe*sumd)*(arg + tvaa) ! extra r in arg is due to r2a factor in dva

    ! 10
        r10a  = r8a*r2a
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(9,DEVDP)
        suma  = suma  + sump
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(10,DEVDP)
        suma  = suma  + sump
        fd10  = 1.0d0 - pe*suma
        dfd10 = (pe*suma - pe*sumd)*(arg + tvaa) ! extra r in arg is due to r2a factor in dva

        V_b6  = - fd6*c6*r6a
        V_b8  = - fd8*c8*r8a
        V_b10 = - fd10*c10*r10a

        dvbb  = 6.0d0*V_b6 + 8.0d0*V_b8 + 10.0d0*V_b10  &
                  + c6*r6a*dfd6 + c8*r8a*dfd8 + c10*r10a*dfd10

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%pen_ene = geo%pen_ene + V_pe
            geo%ind_ene = geo%ind_ene + V_in
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_b6 + V_b8 + V_b10

            dva = r2a*( dvee + dvpe + dvin + dvaa + dvbb )
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * (V_b6 + V_b8 + V_b10)

            ! FIXME - no penetration energy?
            dva = r2a*( inv_scee*dvee + inv_scnb*(dvaa + dvbb) )
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

end subroutine ffdev_gradient_nb_EXP_DISPTT

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_DISPTT

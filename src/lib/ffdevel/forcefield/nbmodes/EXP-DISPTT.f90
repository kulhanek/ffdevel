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
    integer         :: ip,i,j,k,dt
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3,upe,tb
    real(DEVDP)     :: r2,r2a,r,r6a,r8a,r10a,c6,c8,c10,fd6,fd8,fd10,pe,arg,sump
    real(DEVDP)     :: V_aa,V_bb,V_ee,V_pe,er,pr
    real(DEVDP)     :: z1,z2,q1,q2,pfa1,pfa2,pfb1,pfb2,pepa1,pepa2,pepb1,pepb2,pee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
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

        r2  = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3

        r2a = 1.0d0/r2
        r   = sqrt(r2)

    ! electrostatics
        crgij = top%nb_list(ip)%crgij
        V_ee  = crgij/r

        V_pe  = 0.0d0
        if( pen_enabled ) then
            z1      = top%nb_list(ip)%z1
            q1      = top%nb_list(ip)%q1
            pepa1   = top%nb_list(ip)%pepa1
            pepb1   = top%nb_list(ip)%pepb1

            z2      = top%nb_list(ip)%z2
            q2      = top%nb_list(ip)%q2
            pepa2   = top%nb_list(ip)%pepa2
            pepb2   = top%nb_list(ip)%pepb2

            pfa1 = 1.0d0 - exp(-pepa1*r)
            pfb1 = 1.0d0 - exp(-pepb1*r)
            pfa2 = 1.0d0 - exp(-pepa2*r)
            pfb2 = 1.0d0 - exp(-pepb2*r)

            pee = z1*z2 - z1*(z2-q2)*pfa2 - (z1-q1)*z2*pfa1 &
                 + (z1-q1)*(z2-q2)*pfb1*pfb2
            V_pe = pee/r - V_ee
        end if

    ! repulsion
        pa   = top%nb_list(ip)%pa
        pb   = top%nb_list(ip)%pb

        select case(exp_mode)
            case(EXP_BM)
                er   = pb*r
                upe  = exp(-er)
                V_aa = pa*upe
            case(EXP_DO)
                er   = pb*r
                pr   = 1.0d0 + er + er**2/3.0d0
                upe  = exp(-er)
                V_aa = pa*pr*upe
            case(EXP_WO)
                er   = pb*r
                pr   = (1.0d0 + er/2.0d0 + er**2/12.0d0)**2/r
                upe  = exp(-er)
                V_aa = pa*pr*upe
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not defined EXP mode in ffdev_energy_nb_EXP_DISPTT!')
        end select

    ! dispersion
        ! dispersion coefficients
        c6   = top%nb_list(ip)%c6
        c8   = top%nb_list(ip)%c8
        c10  = top%nb_list(ip)%c10

        ! TT damping factor
        tb   = top%nb_list(ip)%tb

        arg = tb*r
        pe  = exp(-arg)

    ! 6
        r6a = r2a*r2a*r2a
        fd6 = 1.0d0 - pe
        sump = 1.0d0
        do k=1,6
            sump = sump * arg / real(k,DEVDP)
            fd6 = fd6 - pe*sump
        end do

    ! 8
        r8a   = r6a*r2a
        fd8  = fd6
        sump = sump * arg / real(7,DEVDP)
        fd8  = fd8 - pe*sump
        sump = sump * arg / real(8,DEVDP)
        fd8  = fd8 - pe*sump

    ! 10
        r10a = r8a*r2a
        fd10 = fd8
        sump = sump * arg / real(9,DEVDP)
        fd10 = fd10 - pe*sump
        sump = sump * arg / real(10,DEVDP)
        fd10 = fd10 - pe*sump

        V_bb = - fd6*c6*r6a - fd8*c8*r8a - fd10*c10*r10a

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%pen_ene = geo%pen_ene + V_pe
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            ! FIXME - what about penetration energy??
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb
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
    integer         :: ip,i,j,k
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3,tb
    real(DEVDP)     :: r2,r2a,r,r6a,r8a,r10a,c6,c8,c10
    real(DEVDP)     :: V_aa,V_bb,pe,V_ee,arg,upe,sump,suma
    real(DEVDP)     :: fd6,fd8,fd10,er,pr,V_pe
    real(DEVDP)     :: z1,z2,q1,q2,pfa1,pfa2,pfb1,pfb2,pepa1,pepa2,pepb1,pepb2,pee
    ! --------------------------------------------------------------------------

    geo%sapt_ele = 0.0d0
    geo%sapt_pen = 0.0d0
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
        r2a  = 1.0d0/r2
        r    = sqrt(r2)

    ! electrostatics
        crgij = top%sapt_list(ip)%crgij
        V_ee = crgij/r

        V_pe  = 0.0d0
        if( pen_enabled ) then
            z1      = top%sapt_list(ip)%z1
            q1      = top%sapt_list(ip)%q1
            pepa1   = top%sapt_list(ip)%pepa1
            pepb1   = top%sapt_list(ip)%pepb1
            z2      = top%sapt_list(ip)%z2
            q2      = top%sapt_list(ip)%q2
            pepa2   = top%sapt_list(ip)%pepa2
            pepb2   = top%sapt_list(ip)%pepb2

            pfa1 = 1.0d0 - exp(-pepa1*r)
            pfb1 = 1.0d0 - exp(-pepb1*r)
            pfa2 = 1.0d0 - exp(-pepa2*r)
            pfb2 = 1.0d0 - exp(-pepb2*r)

            pee = z1*z2 - z1*(z2-q2)*pfa2 - (z1-q1)*z2*pfa1 &
                 + (z1-q1)*(z2-q2)*pfb1*pfb2
            V_pe = pee/r - V_ee
        end if

    ! repulsion
        pa   = top%sapt_list(ip)%pa
        pb   = top%sapt_list(ip)%pb

        select case(exp_mode)
            case(EXP_BM)
                er   = pb*r
                upe  = exp(-er)
                V_aa = pa*upe
            case(EXP_DO)
                er   = pb*r
                pr   = 1.0d0 + er + er**2/3.0d0
                upe  = exp(-er)
                V_aa = pa*pr*upe
            case(EXP_WO)
                er   = pb*r
                pr   = (1.0d0 + er/2.0d0 + er**2/12.0d0)**2/r
                upe  = exp(-er)
                V_aa = pa*pr*upe
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not defined EXP mode in ffdev_energy_nb_EXP_DISPTT!')
        end select

    ! dispersion
        ! dispersion coefficients
        c6   = top%sapt_list(ip)%c6
        c8   = top%sapt_list(ip)%c8
        c10  = top%sapt_list(ip)%c10

        ! TT damping factor
        tb   = top%sapt_list(ip)%tb

        arg = tb*r
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

        geo%sapt_ele = geo%sapt_ele + V_ee
        geo%sapt_pen = geo%sapt_pen + V_pe
        geo%sapt_rep = geo%sapt_rep + V_aa
        geo%sapt_dis = geo%sapt_dis + V_bb

    end do

end subroutine ffdev_energy_sapt_EXP_DISPTT

!===============================================================================
! subroutine ffdev_energy_nbpair_EXP_DISPTT
!===============================================================================

real(DEVDP) function ffdev_energy_nbpair_EXP_DISPTT(nbpair,r)

    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    type(NB_PAIR)   :: nbpair
    real(DEVDP)     :: r
    ! --------------------------------------------
    integer         :: k
    real(DEVDP)     :: pa,pb,tb,er,pr
    real(DEVDP)     :: r2,r2a,suma,sump,upe,c6,c8,c10,pe
    real(DEVDP)     :: V_aa,V_bb,r6a,r8a,r10a,arg,fd10,fd8,fd6
    ! --------------------------------------------------------------------------

    r2   = r**2
    r2a  = 1.0d0/r2

! repulsion
    pa   = nbpair%pa
    pb   = nbpair%pb

    select case(exp_mode)
        case(EXP_BM)
            er   = pb*r
            upe  = exp(-er)
            V_aa = pa*upe
        case(EXP_DO)
            er   = pb*r
            pr   = 1.0d0 + er + er**2/3.0d0
            upe  = exp(-er)
            V_aa = pa*pr*upe
        case(EXP_WO)
            er   = pb*r
            pr   = (1.0d0 + er/2.0d0 + er**2/12.0d0)**2/r
            upe  = exp(-er)
            V_aa = pa*pr*upe
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not defined EXP mode in ffdev_energy_nb_EXP_DISPTT!')
    end select

! dispersion
    ! dispersion coefficients
    c6   = nbpair%c6
    c8   = nbpair%c8
    c10  = nbpair%c10

    ! TT damping factor
    tb   = nbpair%tb

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

    ! FIXME - what about penetration energy?

    ffdev_energy_nbpair_EXP_DISPTT = V_aa + V_bb

end function ffdev_energy_nbpair_EXP_DISPTT

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
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3,upe,tb
    real(DEVDP)     :: r2,r,r6a,r8a,r10a,c6,c8,c10,fd6,fd8,fd10,pe,arg,sump
    real(DEVDP)     :: V_aa,V_b6,V_b8,V_b10,V_ee,dva,r2a,dfd6,dfd8,dfd10,sumd,suma
    real(DEVDP)     :: dvee,dvpe,dvaa,dvbb,er,pr,V_pe
    real(DEVDP)     :: z1,z2,q1,q2,pfa1,pfa2,pfb1,pfb2,pepa1,pepa2,pepb1,pepb2,pee
    ! --------------------------------------------------------------------------

    geo%ele14_ene = 0.0d0
    geo%rep14_ene = 0.0d0
    geo%dis14_ene = 0.0d0

    geo%ele_ene = 0.0d0
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

        V_pe  = 0.0d0
        dvpe  = 0.0d0
        if( pen_enabled ) then
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

            ! full electrostatic with charge penetration
            pee = z1*z2 - z1*(z2-q2)*(1.0d0-pfa2) - (z1-q1)*z2*(1.0d0-pfa1) &
                 + (z1-q1)*(z2-q2)*(1.0d0-pfb1)*(1.0d0-pfb2)
            ! subtract MM energy
            V_pe = pee/r - V_ee
            ! derivatives
            dvpe = V_pe + z1*(z2-q2)*pfa2*pepa2 + (z1-q1)*z2*pfa1*pepa1 &
                 + (z1-q1)*(z2-q2)*(-pfb1*pepb1 - pfb2*pepb2 + pfb1*pfb2*pepb2 + pfb1*pepb1*pfb2)
        end if

    ! repulsion
        pa   = top%nb_list(ip)%pa
        pb   = top%nb_list(ip)%pb

        select case(exp_mode)
            case(EXP_BM)
                er   = pb*r
                upe  = exp(-er)
                V_aa = pa*upe
                dvaa = V_aa*pb*r
         !------------
            case(EXP_DO)
                er   = pb*r
                pr   = 1.0d0 + er + er**2/3.0d0
                upe  = exp(-er)
                V_aa = pa*pr*upe
                dvaa = V_aa*pb*r - r*pa*upe*(pb + 2.0d0*pb*pb*r/3.0d0)
        !------------
            case(EXP_WO)
                er   = pb*r
                pr   = (1.0d0 + er/2.0d0 + er**2/12.0d0)**2/r
                upe  = exp(-er)
                V_aa = pa*pr*upe
                dvaa = V_aa*pb*r &
                     - r*pa*upe*(-r2a + 5.0d0/12.0d0*pb**2 + 2.0d0/12.0d0*pb**3*r + 3.0d0/144.0d0*pb**4 * r2)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not defined EXP mode in ffdev_energy_nb_EXP_DISPTT!')
        end select

    ! dispersion
        ! dispersion coefficients
        c6   = top%nb_list(ip)%c6
        c8   = top%nb_list(ip)%c8
        c10  = top%nb_list(ip)%c10

        ! TT damping factor
        tb   = top%nb_list(ip)%tb

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
        dfd6  = (pe*suma - pe*sumd)*arg ! extra r in arg is due to r2a factor in dva

    ! 8
        r8a   = r6a*r2a
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(7,DEVDP)
        suma  = suma  + sump
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(8,DEVDP)
        suma  = suma  + sump
        fd8   = 1.0d0 - pe*suma
        dfd8  = (pe*suma - pe*sumd)*arg ! extra r in arg is due to r2a factor in dva

    ! 10
        r10a  = r8a*r2a
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(9,DEVDP)
        suma  = suma  + sump
        sumd  = sumd  + sump ! sump is one item delayed in derivatives
        sump  = sump  * arg / real(10,DEVDP)
        suma  = suma  + sump
        fd10  = 1.0d0 - pe*suma
        dfd10 = (pe*suma - pe*sumd)*arg ! extra r in arg is due to r2a factor in dva

        V_b6  = - fd6*c6*r6a
        V_b8  = - fd8*c8*r8a
        V_b10 = - fd10*c10*r10a

        dvbb  = 6.0d0*V_b6 + 8.0d0*V_b8 + 10.0d0*V_b10  &
                  + c6*r6a*dfd6 + c8*r8a*dfd8 + c10*r10a*dfd10

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%pen_ene = geo%pen_ene + V_pe
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_b6 + V_b8 + V_b10

            dva = r2a*( dvee + dvpe + dvaa + dvbb )
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

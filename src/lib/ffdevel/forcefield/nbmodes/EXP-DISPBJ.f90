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

module ffdev_nbmode_EXP_DISPBJ

use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_energy_nb_EXP_DISPBJ
!===============================================================================

subroutine ffdev_energy_nb_EXP_DISPBJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: inv_scee,inv_scnb,pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee,V_pe,er,pr,upe
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

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)

    ! electrostatics
        crgij = top%nb_list(ip)%crgij
        V_ee  = crgij/r

        V_pe  = 0.0d0
        if( pen_enabled ) then
            ! FIXME
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

        rc6  = top%nb_list(ip)%rc6
        rc8  = top%nb_list(ip)%rc8
        rc10 = top%nb_list(ip)%rc10

        r6   = r2*r2*r2
        r8   = r6*r2
        r10  = r8*r2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_bb = - c6*r6i - c8*r8i - c10*r10i

        ! write(*,*) 'ENE  : ',ip,r,pa,pb,c6,c8,c10,rc6,rc8,rc10,V_aa,V_bb

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

end subroutine ffdev_energy_nb_EXP_DISPBJ

!===============================================================================
! subroutine ffdev_energy_sapt_EXP_DISPBJ
!===============================================================================

subroutine ffdev_energy_sapt_EXP_DISPBJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2,r,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee,V_pe,er,pr,upe
    ! --------------------------------------------------------------------------

    geo%sapt_ele  = 0.0d0
    geo%sapt_pen  = 0.0d0
    geo%sapt_rep  = 0.0d0
    geo%sapt_dis  = 0.0d0

    do ip=1,top%sapt_size
        i     = top%sapt_list(ip)%ai
        j     = top%sapt_list(ip)%aj

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r    = sqrt(r2)

    ! electrostatics
        crgij = top%sapt_list(ip)%crgij
        V_ee =   crgij/r

        V_pe  = 0.0d0
        if( pen_enabled ) then
            ! FIXME
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

        rc6  = top%sapt_list(ip)%rc6
        rc8  = top%sapt_list(ip)%rc8
        rc10 = top%sapt_list(ip)%rc10

        r6   = r2*r2*r2
        r8   = r6*r2
        r10  = r8*r2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_bb = - c6*r6i - c8*r8i - c10*r10i

        ! write(*,*) 'SAPT: ',ip,r,pa,pb,c6,c8,c10,rc6,rc8,rc10,V_aa,V_bb

        geo%sapt_ele = geo%sapt_ele + V_ee
        geo%sapt_pen = geo%sapt_pen + V_pe
        geo%sapt_rep = geo%sapt_rep + V_aa
        geo%sapt_dis = geo%sapt_dis + V_bb

    end do

end subroutine ffdev_energy_sapt_EXP_DISPBJ

!===============================================================================
! subroutine ffdev_energy_nbpair_EXP_DISPBJ
!===============================================================================

real(DEVDP) function ffdev_energy_nbpair_EXP_DISPBJ(nbpair,r)

    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    type(NB_PAIR)   :: nbpair
    real(DEVDP)     :: r
    ! --------------------------------------------
    real(DEVDP)     :: pa,pb
    real(DEVDP)     :: r2,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,er,pr,upe
    ! --------------------------------------------------------------------------

    ! calculate distances
    r2   = r**2

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

! disersion
    ! dispersion coefficients
    c6   = nbpair%c6
    c8   = nbpair%c8
    c10  = nbpair%c10

    rc6  = nbpair%rc6
    rc8  = nbpair%rc8
    rc10 = nbpair%rc10

    r6   = r2*r2*r2
    r8   = r6*r2
    r10  = r8*r2

    r6i  = 1.0d0/(r6+rc6)
    r8i  = 1.0d0/(r8+rc8)
    r10i = 1.0d0/(r10+rc10)

    V_bb = - c6*r6i - c8*r8i - c10*r10i

    ! FIXME - what about penetration energy?

    ffdev_energy_nbpair_EXP_DISPBJ = V_aa + V_bb

end function ffdev_energy_nbpair_EXP_DISPBJ

!===============================================================================
! subroutine ffdev_gradient_nb_EXP_DISPBJ
!===============================================================================

subroutine ffdev_gradient_nb_EXP_DISPBJ(top,geo)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip,i,j,dt
    real(DEVDP)     :: pa,pb,crgij,dxa1,dxa2,dxa3,r2a
    real(DEVDP)     :: r2,r,r6,r8,r10,c6,c8,c10,rc6,rc8,rc10
    real(DEVDP)     :: V_aa,V_bb,r6i,r8i,r10i,V_ee,dva,inv_scnb,inv_scee
    real(DEVDP)     :: dvee,dvpe,dvaa,dvbb,er,pr,upe,V_pe
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

        ! calculate distances
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2   = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a  = 1.0d0/r2
        r    = sqrt(r2)

    ! electrostatics
        crgij = top%nb_list(ip)%crgij
        V_ee  = crgij/r
        dvee  = V_ee

        V_pe  = 0.0d0
        dvpe  = 0.0d0
        if( pen_enabled ) then
            ! FIXME
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

        rc6  = top%nb_list(ip)%rc6
        rc8  = top%nb_list(ip)%rc8
        rc10 = top%nb_list(ip)%rc10

        r6   = r2*r2*r2
        r8   = r6*r2
        r10  = r8*r2

        r6i  = 1.0d0/(r6+rc6)
        r8i  = 1.0d0/(r8+rc8)
        r10i = 1.0d0/(r10+rc10)

        V_bb = - c6*r6i - c8*r8i - c10*r10i
        dvbb = - 6.0d0*c6*r6i*r6i*r6 - 8.0d0*c8*r8i*r8i*r8 - 10.0d0*c10*r10i*r10i*r10

        if( dt .eq. 0 ) then
            geo%ele_ene = geo%ele_ene + V_ee
            geo%pen_ene = geo%pen_ene + V_pe
            geo%rep_ene = geo%rep_ene + V_aa
            geo%dis_ene = geo%dis_ene + V_bb

            dva = r2a*( dvee + dvpe + dvaa + dvbb )
        else
            inv_scee = glb_iscee * top%dihedral_types(dt)%inv_scee
            inv_scnb = glb_iscnb * top%dihedral_types(dt)%inv_scnb

            geo%ele14_ene = geo%ele14_ene + inv_scee * V_ee
            geo%rep14_ene = geo%rep14_ene + inv_scnb * V_aa
            geo%dis14_ene = geo%dis14_ene + inv_scnb * V_bb

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

end subroutine ffdev_gradient_nb_EXP_DISPBJ

! ------------------------------------------------------------------------------

end module ffdev_nbmode_EXP_DISPBJ

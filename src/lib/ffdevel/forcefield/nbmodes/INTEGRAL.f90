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

module ffdev_nbmode_INTEGRAL

use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_do
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_do(r,b1,b2,s,ds)

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: s,ds
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_1s1s_do(r,b1,b2,s,ds)

end subroutine ffdev_nbmode_INTEGRAL_do

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_wo
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_wo(r,b1,b2,s,ds)

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: s,ds
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_1s1s_wo(r,b1,b2,s,ds)

end subroutine ffdev_nbmode_INTEGRAL_wo

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_sci_ez
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_sci_ez(r,b1,sci,dsci)

    implicit none
    real(DEVDP)             :: r,b1
    real(DEVDP)             :: sci,dsci
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_1s1s_sci_ez(r,b1,sci,dsci)

end subroutine ffdev_nbmode_INTEGRAL_sci_ez
!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_sci_ee
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_sci_ee(r,b1,b2,sci,dsci)

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: sci,dsci
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_1s1s_sci_ee(r,b1,b2,sci,dsci)

end subroutine ffdev_nbmode_INTEGRAL_sci_ee

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_1s1s_do
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_1s1s_do(r,pb1,pb2,s,ds)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1,pb2
    real(DEVDP)             :: s,ds
    ! --------------------------------------------
    real(DEVDP)             :: inva2b2,a1,a2,pb,pbe,pbe1,pbe2,t1,pre
    ! --------------------------------------------------------------------------

    ! calculate density overlap integral
    if( abs(pb1 - pb2) .gt. 0.1d0 ) then
        inva2b2 = 1.0d0/(pb1**2 - pb2**2)
        a1      =  4.0d0*(pb1**4)*(pb2**4)*(inva2b2**3) + (pb1**3)*(pb2**4)*(inva2b2**2)*r
        a2      = -4.0d0*(pb2**4)*(pb1**4)*(inva2b2**3) + (pb2**3)*(pb1**4)*(inva2b2**2)*r
        pbe1    = exp(-pb1*r)
        pbe2    = exp(-pb2*r)
        t1      = a1*pbe1 + a2*pbe2
        pre     = 1.0d0/(8.0d0*DEV_PI*r)
        s       = pre * t1
        ds      = pre * ( &
                         -a1*pbe1*pb1 + pbe1*(pb1**3)*(pb2**4)*(inva2b2**2) &
                         -a2*pbe2*pb2 + pbe2*(pb2**3)*(pb1**4)*(inva2b2**2) &
                        ) &
                - t1 * pre / r
    else
        pb  = 0.5d0*(pb1 + pb2)
        pbe = exp(-pb*r)
        a1  = 3.0d0 + 3.0d0*pb*r + (pb*r)**2
        pre = pb**3 / (192.0d0 * DEV_PI)
        s   = pre * a1 * pbe
        ds  = pre * (-a1*pbe*pb +  pbe*(3.0d0*pb + 2.0d0*(pb**2)*r))
    end if

end subroutine ffdev_nbmode_INTEGRAL_1s1s_do

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_1s1s_wo
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_1s1s_wo(r,b1,b2,s,ds)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: s,ds
    ! --------------------------------------------

    ! --------------------------------------------------------------------------


end subroutine ffdev_nbmode_INTEGRAL_1s1s_wo

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_1s1s_sci_ez
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_1s1s_sci_ez(r,pb1,sci,dsci)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1
    real(DEVDP)             :: sci,dsci
    ! --------------------------------------------
    real(DEVDP)             :: pbe1,a1
    ! --------------------------------------------------------------------------

    pbe1 = exp(-pb1*r)
    a1   = 1.0 + 0.5d0*pb1*r
    sci  = pbe1*a1
    dsci = - a1*pbe1*pb1 + 0.5d0*pb1*pbe1

end subroutine ffdev_nbmode_INTEGRAL_1s1s_sci_ez

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_1s1s_sci_ee
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_1s1s_sci_ee(r,pb1,pb2,s,ds)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1,pb2
    real(DEVDP)             :: s,ds
    ! --------------------------------------------
    real(DEVDP)             :: inva2b2,a1,a2,pb,pbe,pbe1,pbe2,t1,pre
    ! --------------------------------------------------------------------------



end subroutine ffdev_nbmode_INTEGRAL_1s1s_sci_ee



!                    pfe1 = exp(-pepa1*r)
!                    pfa1 = pfe1*(1.0 + 0.5d0*pepa1*r)
!
!                    pfe2 = exp(-pepa2*r)
!                    pfa2 = pfe2*(1.0 + 0.5d0*pepa2*r)
!
!                    if( abs(pepa1 - pepa2) .gt. 0.1d0 ) then
!                        a2 = pepa1 ** 2
!                        b2 = pepa2 ** 2
!                        inva2b2 = 1.0d0/(b2 - a2)
!                        ! order in inva2b2 is somewhere opposite, which changes the sign
!                        pfbb =  + pfe1 * (b2*inva2b2)**2 * (1.0d0 - 2.0d0*a2*inva2b2 + 0.5d0*pepa1*r) &
!                                + pfe2 * (a2*inva2b2)**2 * (1.0d0 + 2.0d0*b2*inva2b2 + 0.5d0*pepa2*r)
!
!                        ! derivatives - 1/r prefactor is below in r2a including the reverse sign
!                        dvpe = + z1*(z2-q2)*(pfa2*pepa2 - 0.5d0*pepa2*pfe2)     &
!                               + (z1-q1)*z2*(pfa1*pepa1 - 0.5d0*pepa1*pfe1)     &
!                               - (z1-q1)*(z2-q2)*(                              &
!                               + pfe1 * (b2*inva2b2)**2 * (1.0d0 - 2.0d0*a2*inva2b2 + 0.5d0*pepa1*r) * pepa1 &
!                               - 0.5d0*pepa1*pfe1 * (b2*inva2b2)**2  &
!                               + pfe2 * (a2*inva2b2)**2 * (1.0d0 + 2.0d0*b2*inva2b2 + 0.5d0*pepa2*r) * pepa2 &
!                               - 0.5d0*pepa2*pfe2 * (a2*inva2b2)**2 )
!                    else
!                        pepk = 0.5d0*(pepa1 + pepa2)
!                        pfeb = exp(-pepk*r)
!                        pfbb = + pfeb * ( 1.0d0                 &
!                               + 11.0d0/16.0d0*pepk*r           &
!                               + 3.0d0/16.0d0*(pepk*r)**2       &
!                               + 1.0d0/48.0d0*(pepk*r)**3 )
!
!                        ! derivatives - 1/r prefactor is below in r2a including the reverse sign
!                        dvpe = + z1*(z2-q2)*(pfa2*pepa2 - 0.5d0*pepa2*pfe2)     &
!                               + (z1-q1)*z2*(pfa1*pepa1 - 0.5d0*pepa1*pfe1)     &
!                               - (z1-q1)*(z2-q2)*(pfbb*pepk &
!                               - pfeb*(11.0d0/16.0d0*pepk + 2.0d0*3.0d0/16.0d0*r*(pepk)**2 + 3.0d0/48.0d0*r**2*(pepk)**3 ) )
!                    end if
!
!                    pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
!                          - (z1-q1)*(z2-q2)*pfbb
!                    V_pe = pee/r
!
!                    ! derivatives
!                    dvpe = dvpe + V_pe


!                z1      = nbpair%z1
!                q1      = nbpair%q1
!                pepa1   = nbpair%pepa1
!
!                z2      = nbpair%z2
!                q2      = nbpair%q2
!                pepa2   = nbpair%pepa2
!
!                pfa1 = exp(-pepa1*r)*(1.0 + 0.5d0*pepa1*r)
!                pfa2 = exp(-pepa2*r)*(1.0 + 0.5d0*pepa2*r)
!
!                if( abs(pepa1 - pepa2) .gt. 0.1d0 ) then
!                    a2 = pepa1 ** 2
!                    b2 = pepa2 ** 2
!                    inva2b2 = 1.0d0/(b2 - a2)
!                    ! order in inva2b2 is somewhere opposite, which changes the sign
!                    pfbb =  + exp(-pepa1*r) * (b2*inva2b2)**2 *                 &
!                                 (1.0d0 - 2.0d0*a2*inva2b2 + 0.5d0*pepa1*r)     &
!                            + exp(-pepa2*r) * (a2*inva2b2)**2 *                 &
!                                 (1.0d0 + 2.0d0*b2*inva2b2 + 0.5d0*pepa2*r)
!                else
!                    pepk = 0.5d0*(pepa1 + pepa2)
!                    pfbb =       + exp(-pepk*r) *               &
!                         ( 1.0d0 + 11.0d0/16.0d0*pepk*r         &
!                                 + 3.0d0/16.0d0*(pepk*r)**2     &
!                                 + 1.0d0/48.0d0*(pepk*r)**3 )
!                end if
!
!                pee = + z1*(z2-q2)*pfa2 + (z1-q1)*z2*pfa1 &
!                      - (z1-q1)*(z2-q2)*pfbb
!                V_pe = pee/r

! ------------------------------------------------------------------------------

end module ffdev_nbmode_INTEGRAL

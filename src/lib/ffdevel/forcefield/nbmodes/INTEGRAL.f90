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

subroutine ffdev_nbmode_INTEGRAL_do(r,b1,b2,s,ds,tt,dtt)

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: s,ds,tt,dtt
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_do_1s1s(r,b1,b2,s,ds,tt,dtt)

end subroutine ffdev_nbmode_INTEGRAL_do

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_wo
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_wo(r,b1,b2,s,ds,tt,dtt)

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: s,ds,tt,dtt
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_wo_1s1s(r,b1,b2,s,ds,tt,dtt)

end subroutine ffdev_nbmode_INTEGRAL_wo

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_ci_ez_damp
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_ci_ez_damp(r,b1,ci,dci)

    implicit none
    real(DEVDP)             :: r,b1
    real(DEVDP)             :: ci,dci
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_ci_1s1s_ez_damp(r,b1,ci,dci)

end subroutine ffdev_nbmode_INTEGRAL_ci_ez_damp

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_ci_ee
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_ci_ee_damp(r,b1,b2,ci,dci)

    implicit none
    real(DEVDP)             :: r,b1,b2
    real(DEVDP)             :: ci,dci
    ! --------------------------------------------------------------------------

    call ffdev_nbmode_INTEGRAL_ci_1s1s_ee_damp(r,b1,b2,ci,dci)

end subroutine ffdev_nbmode_INTEGRAL_ci_ee_damp

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_do_1s1s
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_do_1s1s(r,pb1,pb2,s,ds,tt,dtt)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1,pb2
    real(DEVDP)             :: s,ds,tt,dtt
    ! --------------------------------------------
    real(DEVDP)             :: x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15
    real(DEVDP)             :: x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28
    real(DEVDP)             :: x29,x30,x31,x32
    ! --------------------------------------------------------------------------

    ! calculate density overlap integral
    if( abs(pb1 - pb2) .gt. 0.1d0 ) then
        x0  = exp(-pb1*r)
        x1  = pb1**2
        x2  = 1/x1
        x3  = 1.0*x2
        x4  = pb2**2
        x5  = 1/x4
        x6  = 1.0*x5
        x7  = x3 - x6
        x8  = x7**(-2)
        x9  = 1/pb1
        x10 = 1.0*r
        x11 = 4.0*x2*x5
        x12 = x0*(x10*x8*x9 - x11/x7**3)
        x13 = exp(-pb2*r)
        x14 = -x3 + x6
        x15 = x14**(-2)
        x16 = 1/pb2
        x17 = x13*(x10*x15*x16 - x11/x14**3)
        x18 = x12 + x17
        x19 = 0.0397887357729738/r
        x20 = r**2
        x21 = 1/x20
        x22 = x0*x8
        x23 = 1.0*x22*x9
        x24 = x13*x15
        x25 = 1.0*x16*x24
        x26 = pb1*x12
        x27 = pb2*x17
        x28 = x23 + x25 - x26 - x27
        x29 = -0.0397887357729738*x18*x21 + x19*x28
        x30 = 1/x18
        x31 = x29*x30
        x32 = 25.1327412287183*x20

        s   = x18*x19
        ds  = x29
        tt  = -x31*x32
        dtt = -50.2654824574367*r*x31 - x30*x32*(x19*(x1*x12 + x17*x4 - 2.0*x22 - 2.0*x24) &
            - 0.0795774715459477*x21*x28 + 0.0795774715459477*x18/r**3) - x29*x32*(-x23 - x25 + x26 + x27)/x18**2
    else
        x0  = pb1/2
        x1  = pb2/2
        x2  = x0 + x1
        x3  = r*x2
        x4  = x2**2
        x5  = r**2*x4
        x6  = 3.0*x3 + x5 + 3.0
        x7  = exp(-x3)
        x8  = x2**3
        x9  = 0.00165786399054058*x7*x8
        x10 = x6*x9
        x11 = 1.5*pb1
        x12 = 1.5*pb2
        x13 = 2*r*x4
        x14 = x11 + x12 + x13
        x15 = -x0 - x1
        x16 = x10*x15 + x14*x9
        x17 = 1/x8
        x18 = exp(x3)
        x19 = x16*x17*x18
        x20 = 603.18578948924/x6
        x21 = x19*x20
        x22 = r*x18*x20
        x23 = 0.00331572798108115*x7

        s   = x10
        ds  = x16
        tt  = -r*x21
        dtt = -67.0206432765822*r*x19*(-x11 - x12 - x13)/(x3 + 0.333333333333333*x5 + 1)**2 &
            - x16*x22/x4 - x17*x22*(x10*x15**2 + x14*x15*x23*x8 + x2**5*x23) - x21
    end if

end subroutine ffdev_nbmode_INTEGRAL_do_1s1s

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_wo_1s1s
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_wo_1s1s(r,pb1,pb2,s,ds,tt,dtt)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1,pb2
    real(DEVDP)             :: s,ds,tt,dtt
    ! --------------------------------------------
    real(DEVDP)             :: x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15
    real(DEVDP)             :: x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28
    real(DEVDP)             :: x29
    ! --------------------------------------------------------------------------

    if( abs(pb1 - pb2) .gt. 0.1d0 ) then
        x0  = pb1**2
        x1  = pb2**2
        x2  = x0/4.0 - x1/4.0
        x3  = r*x2
        x4  = 0.5*r
        x5  = pb1*exp(-pb2*x4)
        x6  = x5*(-2*pb2 + x3)
        x7  = pb2*exp(-pb1*x4)
        x8  = x7*(2*pb1 + x3)
        x9  = x6 + x8
        x10 = sqrt(pb1**3*pb2**3/r)
        x11 = x2**3
        x12 = x10/x11
        x13 = x12*x9
        x14 = sqrt(r)
        x15 = 0.5/x14
        x16 = r**(3.0/2.0)
        x17 = 1/x16
        x18 = x2*x5
        x19 = x2*x7
        x20 = 0.5*pb1*x8
        x21 = 0.5*pb2*x6
        x22 = x18 + x19 - x20 - x21
        x23 = x12*x15
        x24 = -0.5*x13*x17 + x22*x23
        x25 = 1.0/x10
        x26 = 1.0/x9
        x27 = x11*x24*x25*x26
        x28 = 2.0*x16
        x29 = x11*x25*x28

        s   = x13*x15
        ds  = x24
        tt  = -x27*x28
        dtt = -4.0*x14*x27 - x24*x29*(-x18 - x19 + x20 + x21)/x9**2 &
            - x26*x29*(-1.0*x12*x17*x22 + x23*(-1.0*pb1*x19 - 1.0*pb2*x18 + 0.25*x0*x8 + 0.25*x1*x6) + 1.0*x13/r**(5.0/2.0))
    else
        x0  = pb1/2 + pb2/2
        x1  = 0.5*r*x0
        x2  = exp(-x1)
        x3  = x0**2
        x4  = 0.0833333333333333*r**2*x3 + x1 + 1.0
        x5  = x2*x4
        x6  = 0.166666666666667*x3
        x7  = r*x6
        x8  = 0.25*pb1
        x9  = 0.25*pb2
        x10 = x8 + x9
        x11 = x2*(x10 + x7)
        x12 = -x8 - x9
        x13 = x11 + x12*x5
        x14 = exp(x1)
        x15 = x14/x4
        x16 = x13*x15
        x17 = r*x16

        s   = x5
        ds  = x13
        tt  = -x17
        dtt = -r*x13*x14*(x12 - x7)/x4**2 - r*x15*(2*x11*x12 + 0.0625*x2*x4*(-pb1 - pb2)**2 + x2*x6) - x10*x17 - x16
    end if

end subroutine ffdev_nbmode_INTEGRAL_wo_1s1s

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_1s1s_ci_ez_damp
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_ci_1s1s_ez_damp(r,pb1,ci,dci)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1
    real(DEVDP)             :: ci,dci
    ! --------------------------------------------
    real(DEVDP)             :: x0,x1,x2
    ! --------------------------------------------------------------------------

    ! DOI: 10.1002/jcc.20520, Model 1
    ! charge - nucleus (eq. 20a)
    ! autogenerated: share/ffdevel/forcefield/nbmodes/ci_1s1s_ez_damp.py

    x0  = pb1*r
    x1  = exp(-x0)
    x2  = x1*(x0*0.5d0 + 1.0d0)

    ci  = x2
    dci = pb1*x1*0.5d0 - pb1*x2

end subroutine ffdev_nbmode_INTEGRAL_ci_1s1s_ez_damp

!===============================================================================
! subroutine ffdev_nbmode_INTEGRAL_ci_1s1s_ee_damp
!===============================================================================

subroutine ffdev_nbmode_INTEGRAL_ci_1s1s_ee_damp(r,pb1,pb2,ci,dci)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r,pb1,pb2
    real(DEVDP)             :: ci,dci
    ! --------------------------------------------
    real(DEVDP)             :: x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
    ! --------------------------------------------------------------------------

    ! DOI: 10.1002/jcc.20520, Model 1
    ! charge - charge (eq 19a, eq 20a)

    if( abs(pb1 - pb2) .gt. 0.1d0 ) then
        ! autogenerated: share/ffdevel/forcefield/nbmodes/ci_1s1s_ee_ij_damp.py
        x0  = pb2*r
        x1  = pb2**2
        x2  = pb1**2
        x3  = x1 - x2
        x4  = 2.0d0/x3
        x5  = x3**(-2)
        x6  = x5
        x7  = (pb1**4)*exp(-x0)
        x8  = x6*x7*(0.5d0*x0 + x1*x4 + 1.0d0)
        x9  = pb1*r
        x10 = (pb2**4)*exp(-x9)
        x11 = x10*x6*(-x2*x4 + 0.5d0*x9 + 1.0d0)
        x12 = 0.5d0*x5

        ci  = x11 + x8
        dci = pb1*x10*x12 - pb1*x11 + pb2*x12*x7 - pb2*x8
    else
        ! autogenerated: share/ffdevel/forcefield/nbmodes/ci_1s1s_ee_ii_damp.py
        x0  = pb1*0.5d0
        x1  = pb2*0.5d0
        x2  = x0 + x1
        x3  = r*x2
        x4  = exp(-x3)
        x5  = r**2
        x6  = x2**2
        x7  = x2**3
        x8  = x4*(0.0208333333333333d0*r**3*x7 + 0.6875d0*x3 + 0.1875d0*x5*x6 + 1)

        ci  = x8
        dci = x4*(0.34375d0*pb1 + 0.34375d0*pb2 + 0.375d0*r*x6 + 0.0625d0*x5*x7) + x8*(-x0 - x1)
    end if

end subroutine ffdev_nbmode_INTEGRAL_ci_1s1s_ee_damp

! ------------------------------------------------------------------------------

end module ffdev_nbmode_INTEGRAL

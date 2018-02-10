! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_mmd3

use ffdev_constants
use ffdev_mmd3_dat

contains

! ==============================================================================
! subroutine ffdev_mmd3_init
! ==============================================================================

subroutine ffdev_mmd3_init

    implicit none
    ! --------------------------------------------------------------------------

    call dftd3_init(loc_dftd3_calc,loc_dftd3_input)

end subroutine ffdev_mmd3_init

! ==============================================================================
! function ffdev_mmd3_get_c6
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_c6(top,ti,tj)

    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ti,tj
    ! --------------------------------------------
    integer         :: za,zb
    real(DEVDP)     :: cna,cnb,c6
    ! --------------------------------------------------------------------------

    za  = top%atom_types(ti)%z
    cna = ffdev_mmd3_get_type_cn(top,ti)
    zb  = top%atom_types(tj)%z
    cnb = ffdev_mmd3_get_type_cn(top,tj)

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,za,zb,cna,cnb,c6)

    ffdev_mmd3_get_c6 = c6 * 627.50960803059d0 / 1.889725989d0**6

end function ffdev_mmd3_get_c6

! ==============================================================================
! function ffdev_mmd3_get_c6_by_z
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_c6_by_z(za,cna,zb,cnb)

    implicit none
    integer         :: za,zb
    real(DEVDP)     :: cna,cnb
    ! --------------------------------------------
    real(DEVDP)     :: c6
    ! --------------------------------------------------------------------------

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,za,zb,cna,cnb,c6)

    ffdev_mmd3_get_c6_by_z = c6 * 627.50960803059d0 / 1.889725989d0**6

end function ffdev_mmd3_get_c6_by_z

! ==============================================================================
! function ffdev_mmd3_get_c8
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_c8(top,ti,tj)

    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ti,tj
    ! --------------------------------------------
    integer         :: za,zb
    real(DEVDP)     :: cna,cnb,c8,c6
    ! --------------------------------------------------------------------------

    za  = top%atom_types(ti)%z
    cna = ffdev_mmd3_get_type_cn(top,ti)
    zb  = top%atom_types(tj)%z
    cnb = ffdev_mmd3_get_type_cn(top,tj)

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,za,zb,cna,cnb,c6)

    ffdev_mmd3_get_c8 = 3.0d0*c6*r2r4(za)*r2r4(zb) * 627.50960803059d0 / 1.889725989d0**8

end function ffdev_mmd3_get_c8

! ==============================================================================
! function ffdev_mmd3_get_c8_by_z
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_c8_by_z(za,cna,zb,cnb)

    implicit none
    integer         :: za,zb
    real(DEVDP)     :: cna,cnb
    ! --------------------------------------------
    real(DEVDP)     :: c8,c6
    ! --------------------------------------------------------------------------

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,za,zb,cna,cnb,c6)

    ffdev_mmd3_get_c8_by_z = 3.0d0*c6*r2r4(za)*r2r4(zb) * 627.50960803059d0 / 1.889725989d0**8

end function ffdev_mmd3_get_c8_by_z

! ==============================================================================
! function ffdev_mmd3_get_rcov
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_rcov(z1)

    implicit none
    integer         :: z1
    ! --------------------------------------------
    real(DEVDP)     :: c6
    ! --------------------------------------------------------------------------

    ! these new data (=rcov) are scaled with k2=4./3. and converted a_0 via
    ! autoang=0.52917726d0

    ffdev_mmd3_get_rcov = rcov(z1) * 1.889725989d0 ! switch back to A

end function ffdev_mmd3_get_rcov

! ==============================================================================
! function ffdev_mmd3_get_atom_cn
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_atom_cn(top,ai)

    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ai
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r,rco,rr,damp
    ! --------------------------------------------------------------------------

    ffdev_mmd3_get_atom_cn = 0.0d0

    do i=1,top%nbonds
        if( (top%bonds(i)%ai .eq. ai) .or. (top%bonds(i)%aj .eq. ai) ) then
            if( mmd3_use_frac_cn ) then
                ! bond distance
                r = top%bond_types(top%bonds(i)%bt)%d0
                ! covalent distance in au
                rco = ffdev_mmd3_get_rcov(top%atom_types(top%atoms(top%bonds(i)%ai)%typeid)%z) &
                    + ffdev_mmd3_get_rcov(top%atom_types(top%atoms(top%bonds(i)%aj)%typeid)%z)
                rr = rco/r
                ! counting function exponential has a better long-range behavior than MH
                damp=1.d0/(1.d0+exp(-mmd3_k1*(rr-1.0d0)))
                ! write(*,*) r, rr, rco, damp
                ffdev_mmd3_get_atom_cn = ffdev_mmd3_get_atom_cn + damp
            else
                ffdev_mmd3_get_atom_cn = ffdev_mmd3_get_atom_cn + 1.0d0
            end if
        end if
    end do

end function ffdev_mmd3_get_atom_cn

! ==============================================================================
! function ffdev_mmd3_get_type_cn
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_type_cn(top,ti)

    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ti
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: cn, prev_cn
    logical         :: first
    ! --------------------------------------------------------------------------

    first = .true.

    do i=1,top%natoms
        if( top%atoms(i)%typeid .eq. ti ) then
            cn = ffdev_mmd3_get_atom_cn(top,i)
            if( first ) then
                prev_cn = cn
                first = .false.
            end if
            if( prev_cn .ne. cn ) then
                call ffdev_utils_exit(DEV_OUT,1,'inconsistent CN detected in ffdev_mmd3__get_type_cn!')
            end if
        end if
    end do

    ffdev_mmd3_get_type_cn = cn

end function ffdev_mmd3_get_type_cn

! ==============================================================================
! subroutine ffdev_topology_comb_rules_to_string
! ==============================================================================

character(80) function ffdev_mmd3_damping_to_string(damping)

    use ffdev_utils

    implicit none
    integer  :: damping
    ! --------------------------------------------------------------------------

    select case(damping)
        case(MMD3_BM)
            ffdev_mmd3_damping_to_string = 'Becke–Johnson damping function (BJ)'
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdev_mmd3_damping_to_string!')
    end select

end function ffdev_mmd3_damping_to_string

! ==============================================================================
! subroutine ffdev_mmd3_print_params
! ==============================================================================

subroutine ffdev_mmd3_print_params(top)

    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: nb_mode
    ! --------------------------------------------
    integer         :: i, za, zb
    real(DEVDP)     :: cna, cnb, r0ab, c6, c8
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,510)
    write(DEV_OUT,520) ffdev_mmd3_damping_to_string(mmd3_damping)
    select case(mmd3_damping)
        case(MMD3_BM)
            write(DEV_OUT,550) mmd3_bj_a1
            write(DEV_OUT,560) mmd3_bj_a2
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdev_mmd3_print_params!')
    end select

    write(DEV_OUT,*)
    write(DEV_OUT,620)
    write(DEV_OUT,630)

    do i=1,top%nnb_types
        za  = top%atom_types(top%nb_types(i)%ti)%z
        cna = ffdev_mmd3_get_type_cn(top,top%nb_types(i)%ti)
        zb  = top%atom_types(top%nb_types(i)%tj)%z
        cnb = ffdev_mmd3_get_type_cn(top,top%nb_types(i)%tj)
        c6 = ffdev_mmd3_get_c6(top,top%nb_types(i)%ti,top%nb_types(i)%tj)
        c8 = ffdev_mmd3_get_c8(top,top%nb_types(i)%ti,top%nb_types(i)%tj)

        ! determine R0AB
        select case(mmd3_damping)
            case(MMD3_BM)
                r0ab = sqrt(c8/c6)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdev_mmd3_print_params!')
        end select

        write(DEV_OUT,640) i,top%atom_types(top%nb_types(i)%ti)%name,za,cna, &
                             top%atom_types(top%nb_types(i)%tj)%name,zb,cnb, &
                             c6, c8, r0ab

    end do

510 format('# ~~~ MMD3 parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
520 format('# Damping = ',A)
550 format('# BJ damping parameter a1 = ',F16.7)
560 format('# BJ damping parameter a2 = ',F16.7)

620 format('# ID TypA ZA  CNA  TypB ZA  CNB       C6(AB)         C8(AB)       R0(AB) ')
630 format('# -- ---- -- ----- ---- -- ----- --------------- --------------- --------')
640 format(I4,1X,A4,1X,I2,1X,F5.3,1X,A4,1X,I2,1X,F5.3,1X,E15.7,1X,E15.7,1X,F8.5)

end subroutine ffdev_mmd3_print_params

! ==============================================================================

end module ffdev_mmd3

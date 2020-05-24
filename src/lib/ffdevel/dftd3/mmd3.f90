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

    ffdev_mmd3_get_c6 = c6 * DEV_HARTREE2KCL * DEV_AU2A**6

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

    ffdev_mmd3_get_c6_by_z = c6 * DEV_HARTREE2KCL * DEV_AU2A**6

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
    real(DEVDP)     :: cna,cnb,c6
    ! --------------------------------------------------------------------------

    za  = top%atom_types(ti)%z
    cna = ffdev_mmd3_get_type_cn(top,ti)
    zb  = top%atom_types(tj)%z
    cnb = ffdev_mmd3_get_type_cn(top,tj)

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,za,zb,cna,cnb,c6)

    ffdev_mmd3_get_c8 = 3.0d0*c6*r2r4(za)*r2r4(zb) * DEV_HARTREE2KCL * DEV_AU2A**8

end function ffdev_mmd3_get_c8

! ==============================================================================
! function ffdev_mmd3_get_c8_by_z
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_c8_by_z(za,cna,zb,cnb)

    implicit none
    integer         :: za,zb
    real(DEVDP)     :: cna,cnb
    ! --------------------------------------------
    real(DEVDP)     :: c6
    ! --------------------------------------------------------------------------

    call getc6(maxc,max_elem,loc_dftd3_calc%c6ab,loc_dftd3_calc%mxc,za,zb,cna,cnb,c6)

    ffdev_mmd3_get_c8_by_z = 3.0d0*c6*r2r4(za)*r2r4(zb) * DEV_HARTREE2KCL * DEV_AU2A**8

end function ffdev_mmd3_get_c8_by_z

! ==============================================================================
! function ffdev_mmd3_get_rcov
! ==============================================================================

real(DEVDP) function ffdev_mmd3_get_rcov(z1)

    implicit none
    integer         :: z1
    ! --------------------------------------------------------------------------

    ! these new data (=rcov) are scaled with k2=4./3. and converted a_0 via
    ! autoang=0.52917726d0

    ! FIXME? what unit rcov is??
    ffdev_mmd3_get_rcov = rcov(z1) * DEV_AU2A ! switch back to A

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
                ! covalent distance in A
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
            if( abs(prev_cn-cn) .gt. 0.2d0 ) then
                write(DEV_OUT,*) prev_cn, cn
                call ffdev_utils_exit(DEV_ERR,1,'inconsistent CN detected in ffdev_mmd3_get_type_cn!')
            end if
        end if
    end do

    ffdev_mmd3_get_type_cn = cn

end function ffdev_mmd3_get_type_cn

! ==============================================================================
! subroutine ffdev_mmd3_print_params
! ==============================================================================

subroutine ffdev_mmd3_print_params(top)

    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i, za, zb
    real(DEVDP)     :: cna, cnb, r0ab, c6, c8
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,510)
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
        r0ab = sqrt(c8/c6)
        write(DEV_OUT,640) i,top%atom_types(top%nb_types(i)%ti)%name,za,cna, &
                             top%atom_types(top%nb_types(i)%tj)%name,zb,cnb, &
                             c6, c8, r0ab

    end do

510 format('# ~~~ MMD3 parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
620 format('# ID TypA ZA  CNA  TypB ZA  CNB       C6(AB)         C8(AB)       R0(AB) ')
630 format('# -- ---- -- ----- ---- -- ----- --------------- --------------- --------')
640 format(I4,1X,A4,1X,I2,1X,F5.3,1X,A4,1X,I2,1X,F5.3,1X,E15.7,1X,E15.7,1X,F8.5)

end subroutine ffdev_mmd3_print_params

! ==============================================================================
! subroutine ffdev_mmd3_run_stat
! ==============================================================================

subroutine ffdev_mmd3_run_stat()

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: ti, tj, gti, gtj, i, j, alloc_stat
    real(DEVDP) :: rms,c6,c8
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'MMD3', ':')

    ! these data are always available
    mmd3_data_loaded = .true.

    ! allocate mmd3 pairs
    allocate( mmd3_pairs(ntypes,ntypes), fflj_pairs(ntypes,ntypes), stat = alloc_stat )
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate mmd3_pairs or fflj_pairs!')
    end if

    ! clear data
    do gti=1,ntypes
        do gtj=1,ntypes
            mmd3_pairs(gti,gtj)%c6ave = 0.0d0
            mmd3_pairs(gti,gtj)%c6sig = 0.0d0
            mmd3_pairs(gti,gtj)%c8ave = 0.0d0
            mmd3_pairs(gti,gtj)%c8sig = 0.0d0
            mmd3_pairs(gti,gtj)%Rc = 0.0d0
            mmd3_pairs(gti,gtj)%num = 0

            fflj_pairs(gti,gtj)%eps = 0.0d0
            fflj_pairs(gti,gtj)%r0 = 0.0d0
            fflj_pairs(gti,gtj)%c6 = 0.0d0
            fflj_pairs(gti,gtj)%num = 0
        end do
    end do

    ! gather data
    do i=1,nsets
        if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
            call ffdev_mmd3_print_params(sets(i)%top)
        end if

        do j=1,sets(i)%top%nnb_types
            ti = sets(i)%top%nb_types(j)%ti
            tj = sets(i)%top%nb_types(j)%tj
            gti = sets(i)%top%atom_types(ti)%glbtypeid
            gtj = sets(i)%top%atom_types(tj)%glbtypeid

            c6 = ffdev_mmd3_get_c6(sets(i)%top,ti,tj)
            c8 = ffdev_mmd3_get_c8(sets(i)%top,ti,tj)

            ! accumulate data
            fflj_pairs(gti,gtj)%eps = sets(i)%top%nb_types(j)%eps
            fflj_pairs(gti,gtj)%r0  = sets(i)%top%nb_types(j)%r0
            fflj_pairs(gti,gtj)%c6  = 2.0d0 * fflj_pairs(gti,gtj)%eps * fflj_pairs(gti,gtj)%r0**6
            fflj_pairs(gti,gtj)%num = fflj_pairs(gti,gtj)%num + 1

            mmd3_pairs(gti,gtj)%c6ave = mmd3_pairs(gti,gtj)%c6ave + c6
            mmd3_pairs(gti,gtj)%c6sig = mmd3_pairs(gti,gtj)%c6sig + c6**2
            mmd3_pairs(gti,gtj)%c8ave = mmd3_pairs(gti,gtj)%c8ave + c8
            mmd3_pairs(gti,gtj)%c8sig = mmd3_pairs(gti,gtj)%c8sig + c8**2
            mmd3_pairs(gti,gtj)%num   = mmd3_pairs(gti,gtj)%num + 1

            if( ti .ne. tj ) then
                fflj_pairs(gtj,gti)%eps = fflj_pairs(gti,gtj)%eps
                fflj_pairs(gtj,gti)%r0 = fflj_pairs(gti,gtj)%r0
                fflj_pairs(gtj,gti)%c6 = fflj_pairs(gti,gtj)%c6
                fflj_pairs(gtj,gti)%num = fflj_pairs(gti,gtj)%num

                mmd3_pairs(gtj,gti)%c6ave  = mmd3_pairs(gti,gtj)%c6ave
                mmd3_pairs(gtj,gti)%c6sig  = mmd3_pairs(gti,gtj)%c6sig
                mmd3_pairs(gtj,gti)%c8ave  = mmd3_pairs(gti,gtj)%c8ave
                mmd3_pairs(gtj,gti)%c8sig  = mmd3_pairs(gti,gtj)%c8sig
                mmd3_pairs(gtj,gti)%num    = mmd3_pairs(gti,gtj)%num
            end if

        end do
    end do

    do ti=1,ntypes
        do tj=1,ntypes
            if( mmd3_pairs(ti,tj)%num .gt. 0 ) then

                rms = mmd3_pairs(ti,tj)%num*mmd3_pairs(ti,tj)%c6sig - mmd3_pairs(ti,tj)%c6ave**2;
                if( rms .gt. 0.0d0 ) then
                    rms = sqrt(rms) / real(mmd3_pairs(ti,tj)%num)
                else
                    rms = 0.0d0
                end if
                mmd3_pairs(ti,tj)%c6sig = rms
                mmd3_pairs(ti,tj)%c6ave = mmd3_pairs(ti,tj)%c6ave / real(mmd3_pairs(ti,tj)%num)

                rms = mmd3_pairs(ti,tj)%num*mmd3_pairs(ti,tj)%c8sig - mmd3_pairs(ti,tj)%c8ave**2;
                if( rms .gt. 0.0d0 ) then
                    rms = sqrt(rms) / real(mmd3_pairs(ti,tj)%num)
                else
                    rms = 0.0d0
                end if
                mmd3_pairs(ti,tj)%c8sig = rms
                mmd3_pairs(ti,tj)%c8ave = mmd3_pairs(ti,tj)%c8ave / real(mmd3_pairs(ti,tj)%num)

                if( mmd3_pairs(ti,tj)%c6ave .gt. 0 ) then
                    mmd3_pairs(ti,tj)%Rc = sqrt(mmd3_pairs(ti,tj)%c8ave/mmd3_pairs(ti,tj)%c6ave)
                end if

            end if

        end do
    end do

    call ffdev_mmd3_print

end subroutine ffdev_mmd3_run_stat

! ==============================================================================
! subroutine ffdev_mmd3_print
! ==============================================================================

subroutine ffdev_mmd3_print()

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: ti, tj
    ! --------------------------------------------------------------------------

    ! do we have MMD3 data?
    if( .not. mmd3_data_loaded ) then
        write(DEV_OUT,10)
        return
    end if

    ! dispersion data ----------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,20)
    write(DEV_OUT,*)

    write(DEV_OUT,30)
    write(DEV_OUT,40)

    do ti=1,ntypes
        tj = ti
        write(DEV_OUT,50) trim(types(ti)%name),trim(types(tj)%name),mmd3_pairs(ti,tj)%num, &
                          mmd3_pairs(ti,tj)%c6ave, mmd3_pairs(ti,tj)%c6sig, &
                          mmd3_pairs(ti,tj)%c8ave, mmd3_pairs(ti,tj)%c8sig, mmd3_pairs(ti,tj)%Rc
    end do
    write(DEV_OUT,40)
    do ti=1,ntypes
        do tj=ti+1,ntypes
        write(DEV_OUT,50) trim(types(ti)%name),trim(types(tj)%name),mmd3_pairs(ti,tj)%num, &
                          mmd3_pairs(ti,tj)%c6ave, mmd3_pairs(ti,tj)%c6sig, &
                          mmd3_pairs(ti,tj)%c8ave, mmd3_pairs(ti,tj)%c8sig, mmd3_pairs(ti,tj)%Rc
        end do
    end do
    write(DEV_OUT,40)

 10 format('>>> No MMD3 data available ....')
 20 format('# Dispersion coefficients ...')
 30 format('# TypA TypB Number         <C6>  s(C6)         <C8>  s(C8)     Rc')
 40 format('# ---- ---- ------ ------------ ------ ------------ ------ ------')
 50 format(2X,A4,1X,A4,1X,I6,1X,F12.2,1X,F6.2,1X,F12.2,1X,F6.1,1X,F6.3)

end subroutine ffdev_mmd3_print

! ==============================================================================

end module ffdev_mmd3

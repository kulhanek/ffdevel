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

module ffdev_xdm

use ffdev_xdm_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_xdm_control_nbmanip
! ==============================================================================

subroutine ffdev_xdm_control_nbmanip(string,exec)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_targetset
    use ffdev_topology_dat
    use prmfile
    use ffdev_utils

    implicit none
    character(PRMFILE_MAX_PATH) :: string
    logical                     :: exec
    ! --------------------------------------------
    integer                     :: i,mode
    ! --------------------------------------------------------------------------

    mode = APPLY_XDM_NULL
    if( trim(string) .eq. 'r0' )  mode = APPLY_XDM_R0
    if( trim(string) .eq. 'eps' ) mode = APPLY_XDM_EPS

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XDM parameters', '%')
    write(DEV_OUT,10)  trim(string)

    if( mode .eq. APPLY_XDM_NULL ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unsupported XDM mode '//trim(string)//'!')
    end if

    if( .not. exec ) return ! do not execute

    do i=1,nsets
        if( DebugFFManip ) then
            write(DEV_OUT,*)
            write(DEV_OUT,20) i
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'Original NB parameters', '*')
            call ffdev_topology_info_types(sets(i)%top,1)
        end if

        ! remix parameters
        call ffdev_xdm_apply_parameters(sets(i)%top,mode)

        if( DebugFFManip ) then
            ! new set of parameters
            write(DEV_OUT,*)
            call ffdev_utils_heading(DEV_OUT,'New NB parameters', '*')
            call ffdev_topology_info_types(sets(i)%top,2)
        end if
    end do

    ! update nb parameters
    call ffdev_targetset_reinit_nbparams

10 format('Apply XDM parameters (xdm) = ',A)
20 format('=== SET ',I2.2)

end subroutine ffdev_xdm_control_nbmanip

! ==============================================================================
! function ffdev_xdm_apply_parameters
! ==============================================================================

subroutine ffdev_xdm_apply_parameters(top,xdm_mode)

    use ffdev_utils
    use ffdev_topology
    use ffdev_parameters_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: xdm_mode
    ! --------------------------------------------
    integer         :: i, glbti, glbtj
    ! --------------------------------------------------------------------------

    do i=1,top%nnb_types
        glbti = top%atom_types(top%nb_types(i)%ti)%glbtypeid
        glbtj = top%atom_types(top%nb_types(i)%tj)%glbtypeid

        select case(xdm_mode)
            case(APPLY_XDM_EPS)
                top%nb_types(i)%eps = xdm_pairs(glbti,glbtj)%eps
            case(APPLY_XDM_R0)
                top%nb_types(i)%r0  = xdm_pairs(glbti,glbtj)%Rvdw
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported xdm_mode in ffdev_xdm_apply_parameters!')
        end select
    end do

end subroutine ffdev_xdm_apply_parameters

! ==============================================================================
! subroutine ffdev_xdm_keep_c6
! ==============================================================================

subroutine ffdev_xdm_keep_c6()

    use ffdev_parameters_dat
    use ffdev_xdm_dat

    implicit none
    integer             :: i, j
    real(DEVDP)         :: C6
    ! --------------------------------------------------------------------------

    select case(xdm_C6Mode)
        case(KEEP_XDM_C6_VIA_EPS)
            do i=1,nparams
                if( params(i)%realm .ne. REALM_VDW_R0 ) cycle
                if( params(i)%value .eq. 0.0d0 ) cycle
                ! get C6
                C6 = xdm_C6Scale * DEV_HARTREE2KCL * DEV_AU2A**6 * &
                     xdm_pairs(params(i)%ti,params(i)%tj)%c6ave
                do j=1,nparams
                    if( params(j)%realm .ne. REALM_VDW_EPS ) cycle
                    if( (params(i)%ti .eq. params(j)%ti) .and. (params(i)%tj .eq. params(j)%tj) ) then
                        params(j)%value = 0.5d0 * C6 / params(i)%value**6
                        exit
                    end if
                end do
            end do
            case(KEEP_XDM_C6_VIA_R0)
                do i=1,nparams
                    if( params(i)%realm .ne. REALM_VDW_EPS ) cycle
                    if( params(i)%value .eq. 0.0d0 ) cycle
                    ! get C6
                    C6 = xdm_C6Scale * DEV_HARTREE2KCL * DEV_AU2A**6 * &
                         xdm_pairs(params(i)%ti,params(i)%tj)%c6ave
                    do j=1,nparams
                        if( params(j)%realm .ne. REALM_VDW_R0 ) cycle
                        if( (params(i)%ti .eq. params(j)%ti) .and. (params(i)%tj .eq. params(j)%tj) ) then
                            params(j)%value = (0.5d0 * C6 / params(i)%value)**(1.0d0/6.0d0)
                            exit
                        end if
                    end do
                end do
        case(LEFT_XDM_C6)
            ! nothing to do
        case default
            ! nothing to do
    end select

end subroutine ffdev_xdm_keep_c6

! ==============================================================================
! subroutine ffdev_xdm_run_stat
! ==============================================================================

subroutine ffdev_xdm_run_stat()

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: ti, tj, ai, aj, i, j, num, alloc_stat
    real(DEVDP) :: c6sum, c8sum, c10sum, rms, pol
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XDM', ':')

    ! do we have XDM data?
    xdm_data_loaded = .false.

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_xdm_loaded ) then
                xdm_data_loaded = .true.
                exit
            end if
        end do
    end do
    if( .not. xdm_data_loaded ) then
        write(DEV_OUT,10)
        return
    end if

    ! allocate xdm pairs
    allocate( xdm_pairs(ntypes,ntypes), xdm_atoms(ntypes), stat = alloc_stat )
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,trim('Unable to allocate xdm_pairs or xdm_atoms'))
    end if

    ! clear data
    do ti=1,ntypes
        xdm_atoms(ti)%vave = 0.0d0
        xdm_atoms(ti)%vsig = 0.0d0
        xdm_atoms(ti)%v0ave = 0.0d0
        xdm_atoms(ti)%v0sig = 0.0d0
        xdm_atoms(ti)%p0ave = 0.0d0
        xdm_atoms(ti)%p0sig = 0.0d0
        xdm_atoms(ti)%Rvdw  = 0.0d0
        xdm_atoms(ti)%pol = 0.0d0
        xdm_atoms(ti)%num = 0
        do tj=1,ntypes
            xdm_pairs(ti,tj)%c6ave = 0.0d0
            xdm_pairs(ti,tj)%c6sig = 0.0d0
            xdm_pairs(ti,tj)%c8ave = 0.0d0
            xdm_pairs(ti,tj)%c8sig = 0.0d0
            xdm_pairs(ti,tj)%c10ave = 0.0d0
            xdm_pairs(ti,tj)%c10sig = 0.0d0
            xdm_pairs(ti,tj)%num = 0
        end do
    end do

    ! gather data
    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! do we have data?
            if( .not. sets(i)%geo(j)%trg_xdm_loaded ) cycle

            do ai=1,sets(i)%geo(j)%natoms
                ! get types
                ti = sets(i)%top%atom_types(sets(i)%top%atoms(ai)%typeid)%glbtypeid

                ! accumulate data
                xdm_atoms(ti)%vave  = xdm_atoms(ti)%vave + sets(i)%geo(j)%trg_xdm_vol(ai)
                xdm_atoms(ti)%vsig  = xdm_atoms(ti)%vsig + sets(i)%geo(j)%trg_xdm_vol(ai)**2
                xdm_atoms(ti)%v0ave = xdm_atoms(ti)%v0ave + sets(i)%geo(j)%trg_xdm_vol0(ai)
                xdm_atoms(ti)%v0sig = xdm_atoms(ti)%v0sig + sets(i)%geo(j)%trg_xdm_vol0(ai)**2
                xdm_atoms(ti)%p0ave = xdm_atoms(ti)%p0ave + sets(i)%geo(j)%trg_xdm_pol0(ai)
                xdm_atoms(ti)%p0sig = xdm_atoms(ti)%p0sig + sets(i)%geo(j)%trg_xdm_pol0(ai)**2
                xdm_atoms(ti)%num   = xdm_atoms(ti)%num + 1

                do aj=ai,sets(i)%geo(j)%natoms

                    ! get types
                    tj = sets(i)%top%atom_types(sets(i)%top%atoms(aj)%typeid)%glbtypeid

                    ! accumulate data
                    xdm_pairs(ti,tj)%c6ave = xdm_pairs(ti,tj)%c6ave + sets(i)%geo(j)%trg_xdm_c6(ai,aj)
                    xdm_pairs(ti,tj)%c6sig = xdm_pairs(ti,tj)%c6sig + sets(i)%geo(j)%trg_xdm_c6(ai,aj)**2
                    xdm_pairs(ti,tj)%c8ave = xdm_pairs(ti,tj)%c8ave + sets(i)%geo(j)%trg_xdm_c8(ai,aj)
                    xdm_pairs(ti,tj)%c8sig = xdm_pairs(ti,tj)%c8sig + sets(i)%geo(j)%trg_xdm_c8(ai,aj)**2
                    xdm_pairs(ti,tj)%c10ave = xdm_pairs(ti,tj)%c10ave + sets(i)%geo(j)%trg_xdm_c10(ai,aj)
                    xdm_pairs(ti,tj)%c10sig = xdm_pairs(ti,tj)%c10sig + sets(i)%geo(j)%trg_xdm_c10(ai,aj)**2
                    xdm_pairs(ti,tj)%num = xdm_pairs(ti,tj)%num + 1

                    ! complete matrix
                    if( ti .ne. tj ) then
                        xdm_pairs(tj,ti)%c6ave  = xdm_pairs(ti,tj)%c6ave
                        xdm_pairs(tj,ti)%c6sig  = xdm_pairs(ti,tj)%c6sig
                        xdm_pairs(tj,ti)%c8ave  = xdm_pairs(ti,tj)%c8ave
                        xdm_pairs(tj,ti)%c8sig  = xdm_pairs(ti,tj)%c8sig
                        xdm_pairs(tj,ti)%c10ave = xdm_pairs(ti,tj)%c10ave
                        xdm_pairs(tj,ti)%c10sig = xdm_pairs(ti,tj)%c10sig
                        xdm_pairs(tj,ti)%num    = xdm_pairs(ti,tj)%num
                    end if
                end do
            end do
        end do
    end do

    ! finish statistics for dispersion cooeficients
    do ti=1,ntypes

        if( xdm_atoms(ti)%num .gt. 0 ) then
            rms = xdm_atoms(ti)%num*xdm_atoms(ti)%vsig - xdm_atoms(ti)%vave**2;
            if( rms .gt. 0.0d0 ) then
                rms = sqrt(rms) / real(xdm_atoms(ti)%num)
            else
                rms = 0.0d0
            end if
            xdm_atoms(ti)%vsig = rms
            xdm_atoms(ti)%vave = xdm_atoms(ti)%vave / real(xdm_atoms(ti)%num)

            rms = xdm_atoms(ti)%num*xdm_atoms(ti)%v0sig - xdm_atoms(ti)%v0ave**2;
            if( rms .gt. 0.0d0 ) then
                rms = sqrt(rms) / real(xdm_atoms(ti)%num)
            else
                rms = 0.0d0
            end if
            xdm_atoms(ti)%v0sig = rms
            xdm_atoms(ti)%v0ave = xdm_atoms(ti)%v0ave / real(xdm_atoms(ti)%num)

            rms = xdm_atoms(ti)%num*xdm_atoms(ti)%p0sig - xdm_atoms(ti)%p0ave**2;
            if( rms .gt. 0.0d0 ) then
                rms = sqrt(rms) / real(xdm_atoms(ti)%num)
            else
                rms = 0.0d0
            end if
            xdm_atoms(ti)%p0sig = rms
            xdm_atoms(ti)%p0ave = xdm_atoms(ti)%p0ave / real(xdm_atoms(ti)%num)

            ! final polarizability
            xdm_atoms(ti)%pol = xdm_atoms(ti)%vave*xdm_atoms(ti)%p0ave / xdm_atoms(ti)%v0ave

            ! rvdw
            xdm_atoms(ti)%Rvdw = 2.0d0 * DEV_AU2A * xdm_rvdw_fac*xdm_atoms(ti)%pol**(1.0d0/7.0d0)
        end if

    end do

    do ti=1,ntypes
        do tj=1,ntypes
            if( xdm_pairs(ti,tj)%num .gt. 0 ) then

                rms = xdm_pairs(ti,tj)%num*xdm_pairs(ti,tj)%c6sig - xdm_pairs(ti,tj)%c6ave**2;
                if( rms .gt. 0.0d0 ) then
                    rms = sqrt(rms) / real(xdm_pairs(ti,tj)%num)
                else
                    rms = 0.0d0
                end if
                xdm_pairs(ti,tj)%c6sig = rms
                xdm_pairs(ti,tj)%c6ave = xdm_pairs(ti,tj)%c6ave / real(xdm_pairs(ti,tj)%num)

                rms = xdm_pairs(ti,tj)%num*xdm_pairs(ti,tj)%c8sig - xdm_pairs(ti,tj)%c8ave**2;
                if( rms .gt. 0.0d0 ) then
                    rms = sqrt(rms) / real(xdm_pairs(ti,tj)%num)
                else
                    rms = 0.0d0
                end if
                xdm_pairs(ti,tj)%c8sig = rms
                xdm_pairs(ti,tj)%c8ave = xdm_pairs(ti,tj)%c8ave / real(xdm_pairs(ti,tj)%num)

                rms = xdm_pairs(ti,tj)%num*xdm_pairs(ti,tj)%c10sig - xdm_pairs(ti,tj)%c10ave**2;
                if( rms .gt. 0.0d0 ) then
                    rms = sqrt(rms) / real(xdm_pairs(ti,tj)%num)
                else
                    rms = 0.0d0
                end if
                xdm_pairs(ti,tj)%c10sig = rms
                xdm_pairs(ti,tj)%c10ave = xdm_pairs(ti,tj)%c10ave / real(xdm_pairs(ti,tj)%num)
            end if

            ! rvdw
            pol = 0.5d0 * (xdm_atoms(ti)%pol + xdm_atoms(tj)%pol)
            xdm_pairs(ti,tj)%Rvdw = 2.0d0 * DEV_AU2A * xdm_rvdw_fac*pol**(1.0d0/7.0d0)

            ! eps
            if( xdm_pairs(ti,tj)%Rvdw .gt. 0 ) then
                xdm_pairs(ti,tj)%eps = xdm_C6Scale * 0.5d0 * DEV_HARTREE2KCL * DEV_AU2A**6 * &
                                       xdm_pairs(ti,tj)%c6ave / xdm_pairs(ti,tj)%Rvdw**6
            else
                xdm_pairs(ti,tj)%eps = 0.0d0
            end if
        end do
    end do

    ! dispersion data ----------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,20)
    write(DEV_OUT,*)

    write(DEV_OUT,30)
    write(DEV_OUT,40)

    do ti=1,ntypes
        tj = ti
        write(DEV_OUT,50) trim(types(ti)%name),trim(types(tj)%name),xdm_pairs(ti,tj)%num, &
                          xdm_pairs(ti,tj)%c6ave, xdm_pairs(ti,tj)%c6sig, &
                          xdm_pairs(ti,tj)%c8ave, xdm_pairs(ti,tj)%c8sig, &
                          xdm_pairs(ti,tj)%c10ave, xdm_pairs(ti,tj)%c10sig
    end do
    write(DEV_OUT,40)
    do ti=1,ntypes
        do tj=ti+1,ntypes
            write(DEV_OUT,50) trim(types(ti)%name),trim(types(tj)%name),xdm_pairs(ti,tj)%num, &
                              xdm_pairs(ti,tj)%c6ave, xdm_pairs(ti,tj)%c6sig, &
                              xdm_pairs(ti,tj)%c8ave, xdm_pairs(ti,tj)%c8sig, &
                              xdm_pairs(ti,tj)%c10ave, xdm_pairs(ti,tj)%c10sig
        end do
    end do
    write(DEV_OUT,40)

    ! atomic data ----------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,120)
    write(DEV_OUT,*)

    write(DEV_OUT,130)
    write(DEV_OUT,140)
    do ti=1,ntypes
        write(DEV_OUT,150) trim(types(ti)%name),xdm_atoms(ti)%num, &
                          xdm_atoms(ti)%p0ave, xdm_atoms(ti)%p0sig, &
                          xdm_atoms(ti)%v0ave, xdm_atoms(ti)%v0sig, &
                          xdm_atoms(ti)%vave, xdm_atoms(ti)%vsig, &
                          xdm_atoms(ti)%pol, xdm_atoms(ti)%Rvdw
    end do


    ! LJ data ----------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,220)
    write(DEV_OUT,*)

    write(DEV_OUT,230)
    write(DEV_OUT,240)
    do ti=1,ntypes
        tj = ti
        write(DEV_OUT,250) trim(types(ti)%name),trim(types(tj)%name),xdm_pairs(ti,tj)%eps, &
                           xdm_pairs(ti,tj)%Rvdw
    end do
    write(DEV_OUT,240)
    do ti=1,ntypes
        do tj=ti+1,ntypes
            write(DEV_OUT,250) trim(types(ti)%name),trim(types(tj)%name),xdm_pairs(ti,tj)%eps, &
                               xdm_pairs(ti,tj)%Rvdw
        end do
    end do
    write(DEV_OUT,240)

 10 format('>>> No XDM data available ....')

 20 format('# Dispersion coefficients ...')
 30 format('# TypA TypB Number         <C6>        s(C6)         <C8>        s(C8)        <C10>       s(C10)')
 40 format('# ---- ---- ------ ------------ ------------ ------------ ------------ ------------ ------------')
 50 format(2X,A4,1X,A4,1X,I6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6)

120 format('# Atomic data ...')
130 format('# TypA Number       <pol0>      s(pol0)         <V0>        s(V0)          <V>         s(V)          pol   Rvdw')
140 format('# ---- ------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------')
150 format(2X,A4,1X,I6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,F6.3)

220 format('# LJ parameters ...')
230 format('# TypA TypB        eps         R0')
240 format('# ---- ---- ---------- ----------')
250 format(2X,A4,1X,A4,1X,F10.5,1X,F10.5)

end subroutine ffdev_xdm_run_stat

! ------------------------------------------------------------------------------

end module ffdev_xdm

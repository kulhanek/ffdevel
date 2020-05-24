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
use ffdev_variables

contains

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
            if( sets(i)%geo(j)%sup_xdm_loaded ) then
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
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate xdm_pairs or xdm_atoms')
    end if

    ! clear data
    do ti=1,ntypes
        xdm_atoms(ti)%vave = 0.0d0
        xdm_atoms(ti)%vsig = 0.0d0
        xdm_atoms(ti)%v0ave = 0.0d0
        xdm_atoms(ti)%v0sig = 0.0d0
        xdm_atoms(ti)%p0ave = 0.0d0
        xdm_atoms(ti)%p0sig = 0.0d0
        xdm_atoms(ti)%pol = 0.0d0
        xdm_atoms(ti)%num = 0
        do tj=1,ntypes
            xdm_pairs(ti,tj)%c6ave = 0.0d0
            xdm_pairs(ti,tj)%c6sig = 0.0d0
            xdm_pairs(ti,tj)%c8ave = 0.0d0
            xdm_pairs(ti,tj)%c8sig = 0.0d0
            xdm_pairs(ti,tj)%c10ave = 0.0d0
            xdm_pairs(ti,tj)%c10sig = 0.0d0
            xdm_pairs(ti,tj)%Rc = 0.0d0
            xdm_pairs(ti,tj)%num = 0
        end do
    end do

    ! gather data
    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! do we have data?
            if( .not. sets(i)%geo(j)%sup_xdm_loaded ) cycle

            do ai=1,sets(i)%geo(j)%natoms
                ! get types
                ti = sets(i)%top%atom_types(sets(i)%top%atoms(ai)%typeid)%glbtypeid

                ! accumulate data
                xdm_atoms(ti)%vave  = xdm_atoms(ti)%vave + sets(i)%geo(j)%sup_xdm_vol(ai)
                xdm_atoms(ti)%vsig  = xdm_atoms(ti)%vsig + sets(i)%geo(j)%sup_xdm_vol(ai)**2

                xdm_atoms(ti)%v0ave = xdm_atoms(ti)%v0ave + sets(i)%geo(j)%sup_xdm_vol0(ai)
                xdm_atoms(ti)%v0sig = xdm_atoms(ti)%v0sig + sets(i)%geo(j)%sup_xdm_vol0(ai)**2

                xdm_atoms(ti)%p0ave = xdm_atoms(ti)%p0ave + sets(i)%geo(j)%sup_xdm_pol0(ai)
                xdm_atoms(ti)%p0sig = xdm_atoms(ti)%p0sig + sets(i)%geo(j)%sup_xdm_pol0(ai)**2

                xdm_atoms(ti)%num   = xdm_atoms(ti)%num + 1

                do aj=ai,sets(i)%geo(j)%natoms

                    ! get types
                    tj = sets(i)%top%atom_types(sets(i)%top%atoms(aj)%typeid)%glbtypeid

                    ! accumulate data
                    xdm_pairs(ti,tj)%c6ave = xdm_pairs(ti,tj)%c6ave +  &
                                             sets(i)%geo(j)%sup_xdm_c6(ai,aj) * DEV_HARTREE2KCL * DEV_AU2A**6
                    xdm_pairs(ti,tj)%c6sig = xdm_pairs(ti,tj)%c6sig + &
                                           ( sets(i)%geo(j)%sup_xdm_c6(ai,aj) * DEV_HARTREE2KCL * DEV_AU2A**6 )**2

                    xdm_pairs(ti,tj)%c8ave = xdm_pairs(ti,tj)%c8ave + &
                                             sets(i)%geo(j)%sup_xdm_c8(ai,aj) * DEV_HARTREE2KCL * DEV_AU2A**8
                    xdm_pairs(ti,tj)%c8sig = xdm_pairs(ti,tj)%c8sig + &
                                           ( sets(i)%geo(j)%sup_xdm_c8(ai,aj) * DEV_HARTREE2KCL * DEV_AU2A**8 )**2

                    xdm_pairs(ti,tj)%c10ave = xdm_pairs(ti,tj)%c10ave + &
                                              sets(i)%geo(j)%sup_xdm_c10(ai,aj) * DEV_HARTREE2KCL * DEV_AU2A**10
                    xdm_pairs(ti,tj)%c10sig = xdm_pairs(ti,tj)%c10sig + &
                                            ( sets(i)%geo(j)%sup_xdm_c10(ai,aj) * DEV_HARTREE2KCL * DEV_AU2A**10 )**2

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

    ! finish statistics for dispersion coefficients
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
            xdm_atoms(ti)%pol = xdm_atoms(ti)%vave * xdm_atoms(ti)%p0ave / xdm_atoms(ti)%v0ave

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

                xdm_pairs(ti,tj)%Rc = sqrt(xdm_pairs(ti,tj)%C8ave/xdm_pairs(ti,tj)%C6ave)

!                xdm_pairs(ti,tj)%Rc = ( sqrt(xdm_pairs(ti,tj)%C8ave/xdm_pairs(ti,tj)%C6ave) &
!                                      + sqrt(xdm_pairs(ti,tj)%C10ave/xdm_pairs(ti,tj)%C8ave) &
!                                      + (xdm_pairs(ti,tj)%C10ave/xdm_pairs(ti,tj)%C6ave)**(1.0/4.0) ) / 3.0d0
            end if

        end do
    end do

    call ffdev_xdm_print

 10 format('>>> No XDM data available ....')

end subroutine ffdev_xdm_run_stat

! ==============================================================================
! subroutine ffdev_xdm_print
! ==============================================================================

subroutine ffdev_xdm_print()

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_mmd3

    implicit none
    integer     :: ti, tj, ai, aj, nfsum
    real(DEVDP) :: f1,f2,f1sum,f2sum
    ! --------------------------------------------------------------------------

    ! do we have XDM data?
    if( .not. xdm_data_loaded ) then
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
        write(DEV_OUT,50) trim(types(ti)%name),trim(types(tj)%name),xdm_pairs(ti,tj)%num, &
                          xdm_pairs(ti,tj)%c6ave, xdm_pairs(ti,tj)%c6sig, &
                          xdm_pairs(ti,tj)%c8ave, xdm_pairs(ti,tj)%c8sig, &
                          xdm_pairs(ti,tj)%c10ave, xdm_pairs(ti,tj)%c10sig, xdm_pairs(ti,tj)%Rc
    end do
    write(DEV_OUT,40)
    do ti=1,ntypes
        do tj=ti+1,ntypes
            write(DEV_OUT,50) trim(types(ti)%name),trim(types(tj)%name),xdm_pairs(ti,tj)%num, &
                              xdm_pairs(ti,tj)%c6ave, xdm_pairs(ti,tj)%c6sig, &
                              xdm_pairs(ti,tj)%c8ave, xdm_pairs(ti,tj)%c8sig, &
                              xdm_pairs(ti,tj)%c10ave, xdm_pairs(ti,tj)%c10sig, xdm_pairs(ti,tj)%Rc
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
                          xdm_atoms(ti)%pol
    end do
   write(DEV_OUT,140)

    write(DEV_OUT,*)
    write(DEV_OUT,220)
    write(DEV_OUT,*)

    write(DEV_OUT,230)
    write(DEV_OUT,240)

    f1sum = 0.0
    f2sum = 0.0
    nfsum = 0

    do ti=1,ntypes
        tj = ti
        f1 = fflj_pairs(ti,tj)%c6 / xdm_pairs(ti,tj)%c6ave
        f1sum = f1sum + f1
        f2 = fflj_pairs(ti,tj)%c6 / mmd3_pairs(ti,tj)%c6ave
        f2sum = f2sum + f2
        nfsum = nfsum + 1
        write(DEV_OUT,250) trim(types(ti)%name),trim(types(tj)%name), &
                              xdm_pairs(ti,tj)%c6ave, xdm_pairs(ti,tj)%c8ave, xdm_pairs(ti,tj)%Rc, &
                              mmd3_pairs(ti,tj)%c6ave, mmd3_pairs(ti,tj)%c8ave, mmd3_pairs(ti,tj)%Rc, &
                              fflj_pairs(ti,tj)%c6,f1,f2

    end do
    write(DEV_OUT,240)
    do ti=1,ntypes
        do tj=ti+1,ntypes
            f1 = fflj_pairs(ti,tj)%c6 / xdm_pairs(ti,tj)%c6ave
            f1sum = f1sum + f1
            f2 = fflj_pairs(ti,tj)%c6 / mmd3_pairs(ti,tj)%c6ave
            f2sum = f2sum + f2
            nfsum = nfsum + 1
            write(DEV_OUT,250) trim(types(ti)%name),trim(types(tj)%name), &
                              xdm_pairs(ti,tj)%c6ave, xdm_pairs(ti,tj)%c8ave, xdm_pairs(ti,tj)%Rc,  &
                              mmd3_pairs(ti,tj)%c6ave, mmd3_pairs(ti,tj)%c8ave, mmd3_pairs(ti,tj)%Rc, &
                              fflj_pairs(ti,tj)%c6,f1,f2
        end do
    end do
    write(DEV_OUT,240)
    write(DEV_OUT,260) f1sum/nfsum, f2sum/nfsum

 10 format('>>> No XDM data available ....')

 20 format('# Dispersion coefficients ...')
 30 format('# TypA TypB Number         <C6>  s(C6)         <C8>  s(C8)        <C10> s(C10)     Rc')
 40 format('# ---- ---- ------ ------------ ------ ------------ ------ ------------ ------ ------')
 50 format(2X,A4,1X,A4,1X,I6,1X,F12.2,1X,F6.2,1X,F12.2,1X,F6.1,1X,F12.1,1X,F6.1,1X,F6.3,1X,F6.3,1X,F6.3)

120 format('# Atomic data ... (all in atomic units)')
130 format('# TypA Number       <pol0>      s(pol0)         <V0>        s(V0)          <V>         s(V)          pol')
140 format('# ---- ------ ------------ ------------ ------------ ------------ ------------ ------------ ------------')
150 format(2X,A4,1X,I6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6,1X,E12.6)

220 format('# XDM/MMD3/FFLJ comparison ...')
230 format('# TypA TypB    XDM(<C6>)    XDM(<C8>)  XDM(Rc)   MMD3(<C6>)   MMD3(<C8>) MMD3(Rc)     FFLJ(C6) FF/XDM  FF/D3')
240 format('# ---- ---- ------------ ------------ -------- ------------ ------------ -------- ------------ ------ ------')
250 format(2X,A4,1X,A4,1X,F12.2,1X,F12.2,1X,F8.3,1X,F12.2,1X,F12.2,1X,F8.3,1X,F12.2,1X,F6.2,1X,F6.2)
260 format(95X,F6.2,1X,F6.2)

end subroutine ffdev_xdm_print

! ------------------------------------------------------------------------------

end module ffdev_xdm

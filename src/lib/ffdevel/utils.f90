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

module ffdev_utils

use ffdev_sizes

implicit none

interface
    subroutine cffdev_print_errors()
    end subroutine cffdev_print_errors
end interface

contains

!===============================================================================
! Subroutine:  ffdev_utils_open
!===============================================================================

subroutine ffdev_utils_open(unitnum, filename, mystatus)

    implicit none
    integer           unitnum      ! logical unit number
    character(*)      filename     ! file name
    character(1)      mystatus     ! N - new, U - unknown, O - old, R - replace
    ! -----------------------------------------------
    integer           ierr
    character(7)      ustatus       ! status keyword
    character(7)      uposition
    ! --------------------------------------------------------------------------

    if( mystatus .eq. 'N' ) then
        ustatus = 'NEW'
        uposition = 'REWIND'
    else if( mystatus .eq. 'O' ) then
        ustatus = 'OLD'
        uposition = 'REWIND'
    else if( mystatus .eq. 'U' ) then
        ustatus = 'UNKNOWN'
        uposition = 'REWIND'
    else if( mystatus .eq. 'R' ) then
        ustatus = 'REPLACE'
        uposition = 'REWIND'
    else if( mystatus .eq. 'A' ) then
        ustatus = 'UNKNOWN'
        uposition = 'APPEND'
    else
        call ffdev_utils_exit(6, 1,'Incorrect file status!')
    endif

    open(unit = unitnum, file = filename, status = ustatus, &
    position = uposition, form = 'FORMATTED', iostat = ierr)

    if( ierr .ne. 0 ) then
        write(6, '(/,a,a)') 'Unable to open file: ',filename
        call ffdev_utils_exit(6, 1)
    endif

    return

end subroutine ffdev_utils_open

!===============================================================================
! Subroutine:  ffdev_utils_fexist
!===============================================================================

logical function ffdev_utils_fexist(filename)

    use ffdev_constants
    use ffdev_variables

    implicit none
    character(*)      filename     ! file name
    ! -----------------------------------------------
    integer           ierr
    ! --------------------------------------------------------------------------

    open(unit = DEV_TEST, file = filename, status = 'OLD', form = 'FORMATTED', iostat = ierr)

    ffdev_utils_fexist = ierr .eq. 0

    if( ffdev_utils_fexist ) close(DEV_TEST)

end function ffdev_utils_fexist

!===============================================================================
! Subroutine:   ffdev_utils_exit
!===============================================================================

subroutine ffdev_utils_exit(unitnum, errcode, message)

    implicit none
    integer                :: unitnum
    integer                :: errcode
    character(*),optional  :: message
    ! --------------------------------------------------------------------------

    if(present(message)) then
        write(unitnum,'(/,A)') '>>> ERROR: ' // message
    end if

    write(unitnum,'(/,A)') '>>> ERROR: Some fatal error occured in FFDevel!'
    write(unitnum,'(A)')   '           Look above for detailed message (if any).'
    write(unitnum,'(A,/)') '           Program execution is terminated.'

    if( errcode .eq. 0 ) then
        stop 0
    else
        stop 1
    end if

end subroutine ffdev_utils_exit

!===============================================================================
! Subroutine:   ffdev_utils_heading
!===============================================================================

subroutine ffdev_utils_heading(iunit, msg, fill)

    implicit none
    integer        :: iunit
    character(*)   :: msg
    character      :: fill
    ! -----------------------------------------------
    integer        :: n,m, i
    ! -------------------------------------------------------------------------

    if( mod(len(msg),2) .eq. 0 ) then
        n = (80 - len(msg) - 2) / 2
        m = n
    else
        n = (80 - len(msg) - 2) / 2
        m = n+1
    end if

    write(iunit,'(80a)') (fill, i=1,n), ' ', msg, ' ', (fill, i=1,m)

end subroutine ffdev_utils_heading

!===============================================================================
! helper functions
!===============================================================================

real(DEVDP) function kdelta(i,j)

    implicit none
    integer             :: i
    integer             :: j
    ! --------------------------------------------------------------------------

    if( i .eq. j ) then
        kdelta = 1.0d0
    else
        kdelta = 0.0d0
    end if

end function kdelta

!===============================================================================

real(DEVDP) function vdot(x,y)

    implicit none
    real(DEVDP)     :: x(3)
    real(DEVDP)     :: y(3)
    ! --------------------------------------------------------------------------

    vdot = x(1)*y(1) + x(2)*y(2) + x(3)*y(3)

end function vdot

!===============================================================================

function vcross(x,y)

    implicit none
    real(DEVDP)     :: x(3)
    real(DEVDP)     :: y(3)
    real(DEVDP)     :: vcross(3)
    ! --------------------------------------------------------------------------

    ! i j k i j
    ! 1 2 3 1 2
    ! 1 2 3 1 2

    ! cross-product
    vcross(1) = x(2)*y(3) - x(3)*y(2)
    vcross(2) = x(3)*y(1) - x(1)*y(3)
    vcross(3) = x(1)*y(2) - x(2)*y(1)

end function vcross

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_utils_header(progname)

    use ffdev_ver
    use ffdev_constants
    use ffdev_variables

    implicit none
    character(*)   :: progname
    ! -----------------------------------------------
    integer        :: datum(8)
    ! --------------------------------------------------------------------------

    ! get current date and time
    call date_and_time(values=datum)

    ! write header
    write (DEV_OUT,*)
    write (DEV_OUT,10)
    call ffdev_utils_heading(DEV_OUT,'*** '//trim(progname)//' ***',' ')
    write (DEV_OUT,20)
    write (DEV_OUT,50) trim(FFDEV_LIBVER)
    write (DEV_OUT,60)
    write (DEV_OUT,70) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)
    write (DEV_OUT,80)

    return

 10 format('|==============================================================================|')
 20 format('|------------------------------------------------------------------------------|')
 50 format('| Version : ',A66,' |')
 60 format('|------------------------------------------------------------------------------|')
 70 format('| Current date ',i4,'-',i2.2,'-',i2.2,' and time ',i2,':',i2.2,':',i2.2,'                                    |')
 80 format('|==============================================================================|')

end subroutine ffdev_utils_header

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_utils_footer(progname)

    use ffdev_constants
    use ffdev_variables

    implicit none
    character(*)   :: progname
    ! -----------------------------------------------
    integer        :: i
    ! --------------------------------------------------------------------------

    write (DEV_OUT,*)
    write (DEV_OUT,'(80a)') ('=',i=1,80)
    call ffdev_utils_heading(DEV_OUT,trim(progname) // ' terminated normally.',' ')
    write (DEV_OUT,'(80a)') ('=',i=1,80)
    write (DEV_OUT,*)

    return

end subroutine ffdev_utils_footer

!===============================================================================

end module ffdev_utils



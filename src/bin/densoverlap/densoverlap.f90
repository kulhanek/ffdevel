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

program ffdev_densoverlap_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_variables

    implicit none
    character(len=MAX_PATH)     :: cubefile         ! cube file
    integer                     :: sep              ! separation
    character(len=MAX_PATH)     :: string, stype
    integer                     :: read_stat
    integer                     :: NX,NY,NZ
    real(DEVDP)                 :: DX,DY,DZ
    real(DEVDP)                 :: conv_factor
    real(DEVDP)                 :: integral
    real(DEVDP),allocatable     :: density(:,:,:)
    logical                     :: gautype
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Electron Density Overlap')

    ! test number of input arguments
    if( command_argument_count() .ne. 3 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (two expected)!')
    end if

    call get_command_argument(1, cubefile)
    call get_command_argument(2, string)
    call get_command_argument(3, stype)

    gautype = .true.
    if( stype .eq. 'orca' ) then
        gautype = .false.
    end if

    ! convert to int
    sep = -1
    read(string,*,iostat = read_stat) sep
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read separation entry (integer number expected)!')
    end if

    ! process control file -----------------------------------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,10) trim(cubefile)
    write(DEV_OUT,20) sep

    ! read cube file
    call read_cube

    write(DEV_OUT,30) sep*DZ*conv_factor
    write(DEV_OUT,40) NX*NY*NZ

    if( (NZ-2*sep) .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Too large separation!')
    end if

    write(DEV_OUT,50) NX*NY*(NZ-2*sep)

    ! calculate overlap integral
    call calc_integral

    write(DEV_OUT,60) integral

    call ffdev_utils_footer('Electron Density Overlap')

 10 format('Electron density cube file    : ',A)
 20 format('Separation (voxels)           : ',I9)
 30 format('Separation (Ang)              : ',F16.6)
 40 format('Number of voxels (all)        : ',I9)
 50 format('Number of voxels (integrated) : ',I9)
 60 format('Overlap integral              : ',E20.8)

contains

!===============================================================================
! subroutine:  print_usage
!===============================================================================

subroutine print_usage()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,'(/,a,/)') '=== [usage] ===================================================================='
    write(DEV_OUT,10)
    write(DEV_OUT,*)

    return

10 format('    densoverlap <cubefile> <separation> <cubetype>')

end subroutine print_usage

!===============================================================================
! subroutine:  read_cube
!===============================================================================

subroutine read_cube()

    implicit none
    integer     :: nat,inum,alloc_stat,i,j,k
    real(DEVDP) :: rnum,rnum1,rnum2
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_CUBE,cubefile,'O')

    ! skip two comments
    read(DEV_CUBE,*,iostat = read_stat) string
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 1 from the cube file!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) string
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 2 from the cube file!')
    end if

    read(DEV_CUBE,*,iostat = read_stat) nat,rnum,rnum,rnum
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 3 from the cube file!')
    end if

    read(DEV_CUBE,*,iostat = read_stat) NX,DX,rnum1,rnum2
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 4 from the cube file!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 4 from the cube file - non-orthogonal grid!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) NY,rnum1,DY,rnum2
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 5 from the cube file!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 5 from the cube file - non-orthogonal grid!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) NZ,rnum1,rnum2,DZ
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 6 from the cube file!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 6 from the cube file - non-orthogonal grid!')
    end if

    ! read atoms
    do i=1,nat
    read(DEV_CUBE,*,iostat = read_stat) inum,rnum,rnum,rnum,rnum
        if( read_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to read line with an atom from the cube file!')
        end if
    end do

    if( (DX .ne. DY) .and. (DY .ne. DZ) ) then
        call ffdev_utils_exit(DEV_ERR,1,'The cube file has to contain a cube grid!')
    end if

    conv_factor = DEV_AU2A
    if( NX .lt. 0 ) then
        conv_factor = 1.0d0
        NX = -NX
    end if
    if( NY .lt. 0 ) then
        conv_factor = 1.0d0
        NY = -NY
    end if
    if( NZ .lt. 0 ) then
        conv_factor = 1.0d0
        NZ = -NZ
    end if

    if( (NX .ne. NY) .and. (NY .ne. NZ) ) then
        call ffdev_utils_exit(DEV_ERR,1,'The cube file has to contain a cube grid!')
    end if

    ! allocate
    allocate(density(NX,NY,NZ), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for electron density!')
    end if

    ! read density

    do i=1,NX
        do j=1,NY
            if( gautype ) then
                read(DEV_CUBE,'(6E13.5)',iostat = read_stat) (density(k,j,i),k=1,NZ)
            else
                read(DEV_CUBE,'(6E14.5)',iostat = read_stat) (density(k,j,i),k=1,NZ)
            end if
            if( read_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to read line with electron density!')
            end if
        end do
    end do

    close(DEV_CUBE)

end subroutine read_cube

!===============================================================================
! subroutine:  calc_integral
!===============================================================================

subroutine calc_integral()

    implicit none
    integer     :: i,j,k
    real(DEVDP) :: psum,lsum
    ! --------------------------------------------------------------------------

    integral = 0.d0
    do k=1,NZ-sep
        psum = 0.0d0
        do j=1,NY
            lsum = 0.0d0
            do i=1,NX
                lsum = lsum + density(i,j,k)*density(i,j,k+sep)
            end do
            psum = psum + lsum
        end do
        integral = integral + psum
    end do

    integral = integral * DX * DY * DZ * conv_factor**3

end subroutine calc_integral

!===============================================================================

end program ffdev_densoverlap_program

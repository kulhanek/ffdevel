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
    character(len=MAX_PATH)     :: cubefile1            ! cube file
    character(len=MAX_PATH)     :: cubefile2            ! cube file
    integer                     :: sep                  ! separation
    character(len=MAX_PATH)     :: string, stype
    integer                     :: read_stat,nat
    real(DEVDP)                 :: PX,PY,PZ
    integer                     :: NX,NY,NZ
    real(DEVDP)                 :: DX,DY,DZ
    real(DEVDP)                 :: conv_factor
    real(DEVDP)                 :: integral
    real(DEVDP),allocatable     :: density1(:,:,:)
    real(DEVDP),allocatable     :: density2(:,:,:)
    logical                     :: symmetrize
    logical                     :: gautype
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Electron Density Overlap')

    ! test number of input arguments
    if( (command_argument_count() .ne. 3) .and. (command_argument_count() .ne. 4)  ) then
        call print_usage()
        call ffdev_utils_exit(DEV_ERR,1,'Incorrect number of arguments was specified (three/four expected)!')
    end if

    call get_command_argument(1, cubefile1)
    if( command_argument_count() .eq. 3 ) then
        call get_command_argument(2, string)
        call get_command_argument(3, stype)
        cubefile2 = cubefile1
    else
        call get_command_argument(2, cubefile2)
        call get_command_argument(3, string)
        call get_command_argument(4, stype)
    end if

    gautype = .true.
    if( (stype .eq. 'orca') .or. (stype .eq. 'sorca') ) then
        gautype = .false.
    end if
    symmetrize = .false.
    if( (stype .eq. 'sgau') .or. (stype .eq. 'sorca') ) then
        symmetrize = .true.
    end if

    ! convert to int
    sep = -1
    read(string,*,iostat = read_stat) sep
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read separation entry (integer number expected)!')
    end if

    ! process control file -----------------------------------------------------
    write(DEV_OUT,*)
    write(DEV_OUT,10) trim(cubefile1)
    write(DEV_OUT,15) trim(cubefile2)

    ! read cube file
    call read_cube1
    call read_cube2

    write(DEV_OUT,17) nat

    write(DEV_OUT,20) sep
    write(DEV_OUT,30) sep*DZ*conv_factor
    write(DEV_OUT,40) NX*NY*NZ

    if( (NZ-2*sep) .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Too large separation!')
    end if

    write(DEV_OUT,50) NX*NY*(NZ-2*sep)

    if( symmetrize ) then
         write(DEV_OUT,55)
        call symmetrize_density(density1)
        call symmetrize_density(density2)
    end if

    ! calculate overlap integral
    call calc_integral

    write(DEV_OUT,60) integral

    call ffdev_utils_footer('Electron Density Overlap')

 10 format('Electron density cube file #1 : ',A)
 15 format('Electron density cube file #2 : ',A)
 17 format('Number of atoms               : ',I9)
 20 format('Separation (voxels)           : ',I9)
 30 format('Separation (Ang)              : ',F16.6)
 40 format('Number of voxels (all)        : ',I9)
 50 format('Number of voxels (integrated) : ',I9)
 55 format('Symmetrizing density ...')
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
    write(DEV_OUT,20)
    write(DEV_OUT,30)
    write(DEV_OUT,*)

    return

10 format('    densoverlap <cubefile1> [cubefile2] <separation> <gau/orca/sgau/sorca>')
20 format('                gau/sgau   - cube file produced by gaussian (force symmetrization)')
30 format('                orca/sorca - cube file produced by orca (force symmetrization)')

end subroutine print_usage

!===============================================================================
! subroutine:  read_cube
!===============================================================================

subroutine read_cube1()

    implicit none
    integer     :: inum,alloc_stat,i,j,k
    real(DEVDP) :: rnum,rnum1,rnum2
    ! --------------------------------------------------------------------------

    nat = 0

    call ffdev_utils_open(DEV_CUBE,cubefile1,'O')

    ! skip two comments
    read(DEV_CUBE,*,iostat = read_stat) string
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 1 from the cube file - cube1!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) string
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 2 from the cube file - cube1!')
    end if

    read(DEV_CUBE,*,iostat = read_stat) nat,PX,PY,PZ
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 3 from the cube file - cube1!')
    end if

    read(DEV_CUBE,*,iostat = read_stat) NX,DX,rnum1,rnum2
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 4 from the cube file - cube1!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 4 from the cube file - non-orthogonal grid - cube1!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) NY,rnum1,DY,rnum2
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 5 from the cube file - cube1!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 5 from the cube file - non-orthogonal grid - cube1!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) NZ,rnum1,rnum2,DZ
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 6 from the cube file - cube1!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 6 from the cube file - non-orthogonal grid - cube1!')
    end if

    ! read atoms
    do i=1,nat
    read(DEV_CUBE,*,iostat = read_stat) inum,rnum,rnum,rnum,rnum
        if( read_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to read line with an atom from the cube file - cube1!')
        end if
    end do

    if( (DX .ne. DY) .and. (DY .ne. DZ) ) then
        call ffdev_utils_exit(DEV_ERR,1,'The cube file has to contain a cube grid - cube1!')
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
        call ffdev_utils_exit(DEV_ERR,1,'The cube file has to contain a cube grid - cube1!')
    end if

    ! allocate
    allocate(density1(NX,NY,NZ), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for electron density - cube1!')
    end if

    ! read density
    do i=1,NX
        do j=1,NY
            if( gautype ) then
                read(DEV_CUBE,'(6E13.5)',iostat = read_stat) (density1(k,j,i),k=1,NZ)
            else
                read(DEV_CUBE,'(6E14.5)',iostat = read_stat) (density1(k,j,i),k=1,NZ)
            end if
            if( read_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to read line with electron density - cube1!')
            end if
        end do
    end do

    close(DEV_CUBE)

end subroutine read_cube1

!===============================================================================
! subroutine:  read_cube2
!===============================================================================

subroutine read_cube2()

    implicit none
    integer     :: inum,alloc_stat,i,j,k
    real(DEVDP) :: rnum,rnum1,rnum2,rnum3
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_CUBE,cubefile2,'O')

    ! skip two comments
    read(DEV_CUBE,*,iostat = read_stat) string
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 1 from the cube file - cube2!')
    end if
    read(DEV_CUBE,*,iostat = read_stat) string
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 2 from the cube file - cube2!')
    end if

    read(DEV_CUBE,*,iostat = read_stat) inum,rnum1,rnum2,rnum3
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 3 from the cube file - cube2!')
    end if
    if( inum .ne. nat ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - number of atoms - cube2!')
    end if

    read(DEV_CUBE,*,iostat = read_stat) inum,rnum1,rnum2,rnum3
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 4 from the cube file - cube2!')
    end if
    if( rnum2 .ne. 0 .and. rnum3 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 4 from the cube file - non-orthogonal grid - cube2!')
    end if
    if( inum .ne. NX ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - NX')
    end if
    if( rnum1 .ne. DX ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - DX')
    end if

    read(DEV_CUBE,*,iostat = read_stat) inum,rnum1,rnum2,rnum3
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 5 from the cube file - cube2!')
    end if
    if( rnum1 .ne. 0 .and. rnum3 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 5 from the cube file - non-orthogonal grid - cube2!')
    end if
    if( inum .ne. NY ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - NY')
    end if
    if( rnum2 .ne. DY ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - DY')
    end if


    read(DEV_CUBE,*,iostat = read_stat) inum,rnum1,rnum2,rnum3
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 6 from the cube file - cube2!')
    end if
    if( rnum1 .ne. 0 .and. rnum2 .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read line 6 from the cube file - non-orthogonal grid - cube2!')
    end if
    if( inum .ne. NZ ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - NZ')
    end if
    if( rnum3 .ne. DZ ) then
        call ffdev_utils_exit(DEV_ERR,1,'cube1 and cube2 inconsistent - DZ')
    end if

    ! read atoms
    do i=1,nat
    read(DEV_CUBE,*,iostat = read_stat) inum,rnum,rnum,rnum,rnum
        if( read_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to read line with an atom from the cube file - cube2!')
        end if
    end do

    ! allocate
    allocate(density2(NX,NY,NZ), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for electron density - cube2!')
    end if

    ! read density
    do i=1,NX
        do j=1,NY
            if( gautype ) then
                read(DEV_CUBE,'(6E13.5)',iostat = read_stat) (density2(k,j,i),k=1,NZ)
            else
                read(DEV_CUBE,'(6E14.5)',iostat = read_stat) (density2(k,j,i),k=1,NZ)
            end if
            if( read_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to read line with electron density - cube2!')
            end if
        end do
    end do

    close(DEV_CUBE)

end subroutine read_cube2

!===============================================================================
! subroutine:  symmetrize_density
!===============================================================================

subroutine symmetrize_density(density)

    implicit none
    real(DEVDP) :: density(:,:,:)
    ! --------------------------------------------
    integer     :: i,j,k
    real(DEVDP) :: dens
    ! --------------------------------------------------------------------------

    do k=1,NZ
        do j=1,NY
            do i=1,NX
                dens = 0.0
                dens = dens + density(i,j,k)
                dens = dens + density(k,i,j)
                dens = dens + density(j,k,i)
                dens = dens / 3.0d0
                density(i,j,k) = dens
                density(k,i,j) = dens
                density(j,k,i) = dens
            end do
        end do
    end do

end subroutine symmetrize_density

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
                lsum = lsum + density1(i,j,k)*density2(i,j,k+sep)
            end do
            psum = psum + lsum
        end do
        integral = integral + psum
    end do

    integral = integral * DX * DY * DZ * conv_factor**3

end subroutine calc_integral

!===============================================================================

end program ffdev_densoverlap_program

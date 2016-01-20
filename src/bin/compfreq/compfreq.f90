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

program ffdev_compfreq_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_topology
    use ffdev_geometry
    use ffdev_hessian
    use ffdev_hessian_utils

    implicit none
    character(len=MAX_PATH) :: crdname1     ! input coordinate name #1
    character(len=MAX_PATH) :: crdname2     ! input coordinate name #2
    type(GEOMETRY)          :: geo1
    type(GEOMETRY)          :: geo2
    double precision        :: rmsd
    integer,allocatable     :: map(:)
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Compare Frequencies')

    ! test number of input arguments
    if( command_argument_count() .ne. 2 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (two expected)!')
    end if

    call get_command_argument(1, crdname1)
    call get_command_argument(2, crdname2)

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,110) trim(crdname1)
    write(DEV_OUT,120) trim(crdname2)

    ! load coordinates #1 ------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input coordinates #1','=')

    call ffdev_geometry_init(geo1)
    call ffdev_geometry_load_point(geo1,crdname1)
    call ffdev_geometry_info_input(geo1)

    if( .not. geo1%trg_grd_loaded ) then
        call ffdev_utils_exit(DEV_OUT,1,'Point #1 does not contain Hessian!')
    end if

    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(geo1)

    ! calculate frequencies
    call ffdev_hessian_allocate(geo1)
    call ffdev_hessian_allocate_freq(geo1)
    geo1%hess = geo1%trg_hess
    call ffdev_hessian_calc_freqs(geo1)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Frequencies #1', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo1)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo1)

    ! load coordinates #2 ------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Input coordinates #2','=')

    call ffdev_geometry_init(geo2)
    call ffdev_geometry_load_point(geo2,crdname2)
    call ffdev_geometry_info_input(geo2)

    if( .not. geo2%trg_grd_loaded ) then
        call ffdev_utils_exit(DEV_OUT,1,'Point #2 does not contain Hessian!')
    end if

    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(geo2)

    ! calculate frequencies
    call ffdev_hessian_allocate(geo2)
    call ffdev_hessian_allocate_freq(geo2)
    geo2%hess = geo2%trg_hess
    call ffdev_hessian_calc_freqs(geo2)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Frequencies #2', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo2)
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs_lin(DEV_OUT,geo2)

    ! final check --------------------------------------------------------------
    if( geo1%natoms .ne. geo2%natoms ) then
        call ffdev_utils_exit(DEV_OUT,1,'Points do not have the same number of atoms!')
    end if

    ! superimpose structures
    call superimpose_pts(geo1,geo2,.true.,rmsd)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Superimposed coordinates #2', ':')
    write(DEV_OUT,130) rmsd
    write(DEV_OUT,*)
    call ffdev_geometry_print_xyz(geo2)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Superimposed vectors #2', ':')
    write(DEV_OUT,*)
    call ffdev_hessian_print_freqs(DEV_OUT,geo2)

    ! find vector mappings
    call find_mapping

    ! compare energy items
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'COMPARISON','=')

    call print_mapping

    call ffdev_utils_footer('Compare Frequencies')

110 format('Input coordinates #1 : ',A)
120 format('Input coordinates #2 : ',A)
130 format('RMSD(#2->#1) = ',F12.6)

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

10 format('    compfreq <xyzcrd#1> <xyzcrd#2>')

end subroutine print_usage

!===============================================================================
! subroutine:  superimpose_pts
!===============================================================================

subroutine superimpose_pts(ref,src,mw,rmsd)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: ref
    type(GEOMETRY)  :: src
    logical         :: mw       ! mass weighted
    real(DEVDP)     :: rmsd
    ! --------------------------------------------
    integer         :: i,info,best,j,k
    real(DEVDP)     :: x1,x2,x3,xr1,xr2,xr3,amass,totmass,itotmass
    real(DEVDP)     :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(DEVDP)     :: f(4,4),u(3,3)
    real(DEVDP)     :: eigenvalues(4),work(26*4)
    real(DEVDP)     :: x2sum,xr2sum,px,py,pz,npx,npy,npz
    ! --------------------------------------------------------------------------

    if( src%natoms .ne. ref%natoms ) then
        call ffdev_utils_exit(DEV_OUT,1,'Two geometries must have the same number of atoms!')
    end if
    if( src%natoms .lt. 2 ) then
        call ffdev_utils_exit(DEV_OUT,1,'At least two atoms must be in both structures!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    x1 = 0.0d0
    x2 = 0.0d0
    x3 = 0.0d0
    xr1 = 0.0d0
    xr2 = 0.0d0
    xr3 = 0.0d0
    totmass = 0.0d0

    do  i = 1, src%natoms
        if( src%z(i) .ne. ref%z(i) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Two geometries must have identical order of atoms!')
        end if
        if( mw ) then
            amass = pt_masses(src%z(i))
        else
            amass = 1.0d0
        end if
        ! source
        x1 = x1 + src%crd(1,i)*amass
        x2 = x2 + src%crd(2,i)*amass
        x3 = x3 + src%crd(3,i)*amass

        ! reference
        xr1 = xr1 + ref%crd(1,i)*amass
        xr2 = xr2 + ref%crd(2,i)*amass
        xr3 = xr3 + ref%crd(3,i)*amass

        totmass = totmass + amass
    end do

    itotmass = 1.0d0 / totmass
    x1 = x1 * itotmass
    x2 = x2 * itotmass
    x3 = x3 * itotmass
    xr1 = xr1 * itotmass
    xr2 = xr2 * itotmass
    xr3 = xr3 * itotmass

    ! calculate correlation matrix -------------------
    r11 = 0.0d0
    r12 = 0.0d0
    r13 = 0.0d0

    r21 = 0.0d0
    r22 = 0.0d0
    r23 = 0.0d0

    r31 = 0.0d0
    r32 = 0.0d0
    r33 = 0.0d0

    x2sum = 0.0d0
    xr2sum = 0.0d0

    do i = 1, src%natoms
        if( mw ) then
            amass = pt_masses(src%z(i))
        else
            amass = 1.0d0
        end if

        x2sum = x2sum + amass*((src%crd(1,i) - x1)**2 &
                             + (src%crd(2,i) - x2)**2 &
                             + (src%crd(3,i) - x3)**2)
        xr2sum = xr2sum + amass*((ref%crd(1,i) - xr1)**2 &
                               + (ref%crd(2,i) - xr2)**2 &
                               + (ref%crd(3,i) - xr3)**2)

        r11 = r11 + amass*(src%crd(1,i) - x1)*(ref%crd(1,i) - xr1)
        r12 = r12 + amass*(src%crd(1,i) - x1)*(ref%crd(2,i) - xr2)
        r13 = r13 + amass*(src%crd(1,i) - x1)*(ref%crd(3,i) - xr3)

        r21 = r21 + amass*(src%crd(2,i) - x2)*(ref%crd(1,i) - xr1)
        r22 = r22 + amass*(src%crd(2,i) - x2)*(ref%crd(2,i) - xr2)
        r23 = r23 + amass*(src%crd(2,i) - x2)*(ref%crd(3,i) - xr3)

        r31 = r31 + amass*(src%crd(3,i) - x3)*(ref%crd(1,i) - xr1)
        r32 = r32 + amass*(src%crd(3,i) - x3)*(ref%crd(2,i) - xr2)
        r33 = r33 + amass*(src%crd(3,i) - x3)*(ref%crd(3,i) - xr3)
    end do

    ! construct matrix for quaterion fitting
    f(1,1) =  r11 + r22 + r33
    f(1,2) =  r23 - r32
    f(1,3) =  r31 - r13
    f(1,4) =  r12 - r21

    f(2,1) =  r23 - r32
    f(2,2) =  r11 - r22 - r33
    f(2,3) =  r12 + r21
    f(2,4) =  r13 + r31

    f(3,1) =  r31 - r13
    f(3,2) =  r12 + r21
    f(3,3) = -r11 + r22 - r33
    f(3,4) =  r23 + r32

    f(4,1) =  r12 - r21
    f(4,2) =  r13 + r31
    f(4,3) =  r23 + r32
    f(4,4) = -r11 - r22 + r33

    ! calculate eignevalues and eigenvectors of matrix f
    eigenvalues(:) = 0d0

    ! now solve eigenproblem
    call dsyev('V','L', 4, f, 4, eigenvalues, work, 26*4, info)

    if( info .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to diagonalize matrix in calculate_rmsdt!')
    end if

    best = 4

    ! calculate rmsdt
    rmsd = sqrt((x2sum + xr2sum - 2.0d0*eigenvalues(best))*itotmass)

    ! rotation matrix ------------------------------
    u(1,1) = f(1,best)**2 + f(2,best)**2 - f(3,best)**2 - f(4,best)**2
    u(1,2) = 2.0d0*( f(2,best)*f(3,best) - f(1,best)*f(4,best) )
    u(1,3) = 2.0d0*( f(2,best)*f(4,best) + f(1,best)*f(3,best) )

    u(2,1) = 2.0d0*( f(2,best)*f(3,best) + f(1,best)*f(4,best) )
    u(2,2) = f(1,best)**2 - f(2,best)**2 + f(3,best)**2 - f(4,best)**2
    u(2,3) = 2.0d0*( f(3,best)*f(4,best) - f(1,best)*f(2,best) )

    u(3,1) = 2.0d0*( f(2,best)*f(4,best) - f(1,best)*f(3,best) )
    u(3,2) = 2.0d0*( f(3,best)*f(4,best) + f(1,best)*f(2,best) )
    u(3,3) = f(1,best)**2 - f(2,best)**2 - f(3,best)**2 + f(4,best)**2

    ! fit geometry --------------------------------
    do i=1,src%natoms
        ! move structure to origin
        px = src%crd(1,i) - x1
        py = src%crd(2,i) - x2
        pz = src%crd(3,i) - x3

        ! rotate structure
        npx = u(1,1)*px + u(1,2)*py + u(1,3)*pz
        npy = u(2,1)*px + u(2,2)*py + u(2,3)*pz
        npz = u(3,1)*px + u(3,2)*py + u(3,3)*pz

        ! move to reference structure COM
        src%crd(1,i) = npx + xr1
        src%crd(2,i) = npy + xr2
        src%crd(3,i) = npz + xr3

        do j=1,src%natoms
            do k=1,3
                ! it is a vector - no origin dependence
                px = src%hess(1,i,k,j)
                py = src%hess(2,i,k,j)
                pz = src%hess(3,i,k,j)

                ! rotate vector
                npx = u(1,1)*px + u(1,2)*py + u(1,3)*pz
                npy = u(2,1)*px + u(2,2)*py + u(2,3)*pz
                npz = u(3,1)*px + u(3,2)*py + u(3,3)*pz

                ! update
                src%hess(1,i,k,j) = npx
                src%hess(2,i,k,j) = npy
                src%hess(3,i,k,j) = npz
            end do
        end do
    end do

end subroutine superimpose_pts

!===============================================================================
! subroutine:  find_mapping
!===============================================================================

subroutine find_mapping()

    implicit none
    integer             :: alloc_stat, i1, j1, i2, j2, max_j, k, l
    double precision    :: max_angle, angle
    integer,allocatable :: map1(:)
    ! --------------------------------------------------------------------------

    allocate(map(3*geo1%natoms),map1(3*geo1%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate map array!')
    end if

    map(:) = 0
    map1(:) = 0

    do i1=1,geo1%natoms
        do j1=1,3
            max_j = 0
            max_angle = 0.0
            do i2=1,geo2%natoms
                do j2=1,3
                    if( map1((i2-1)*3+j2) .gt. 0 ) cycle
                    ! calculate angle
                    angle = 0.0
                    do k=1,geo2%natoms
                        do l=1,3
                            angle = angle + geo1%hess(l,k,j1,i1)*geo2%hess(l,k,j2,i2)
                        end do
                    end do
                    if( (max_angle .lt. abs(angle)) .or. (max_j .eq. 0) ) then
                        max_angle = abs(angle)
                        max_j = (i2-1)*3+j2
                    end if
                end do
            end do
            map1(max_j) = (i1-1)*3+j1
            map((i1-1)*3+j1) = max_j
        end do
    end do

    deallocate(map1)

end subroutine find_mapping

!===============================================================================
! subroutine:  print_mapping
!===============================================================================

subroutine print_mapping()

    implicit none
    integer             :: i, i1, j1, i2, j2, k, l
    double precision    :: serr, lerr, aerr, rmse, angle, diff
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,120)
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=7,3*geo1%natoms
        diff = geo2%freq(map(i))-geo1%freq(i)

        i1 = (i-1) / 3 + 1
        j1 = mod(i-1,3) + 1

        i2 = (map(i)-1) / 3 + 1
        j2 = mod(map(i)-1,3) + 1

        angle = 0.0
        do k=1,geo2%natoms
            do l=1,3
                angle = angle + geo1%hess(l,k,j1,i1)*geo2%hess(l,k,j2,i2)
            end do
        end do
        angle = acos(angle)*DEV_R2D

        write(DEV_OUT,140) i, geo1%freq(i), map(i), geo2%freq(map(i)), angle, diff
        if( serr .gt. abs(diff) ) serr = abs(diff)
        if( lerr .lt. abs(diff) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
    end do

    if( 3*geo1%natoms - 6 .gt. 0 ) then
        aerr = aerr / real(3*geo1%natoms - 6)
        rmse = sqrt(rmse / real(3*geo1%natoms - 6))
    end if

    write(DEV_OUT,130)
    write(DEV_OUT,150) serr
    write(DEV_OUT,160) lerr
    write(DEV_OUT,170) aerr
    write(DEV_OUT,180) rmse

120 format('# Indx#1  Freq#1  Indx#2  Freq#2   Angle   Diff ')
130 format('# ------ -------- ------ -------- ------- ------')
140 format(I8,1X,F8.2,1X,I6,1X,F8.2,1X,F7.1,1X,F6.1)
150 format('# Minimum unsigned difference (SUD)  = ',F9.4)
160 format('# Largest unsigned difference (MUD)  = ',F9.4)
170 format('# Average usigned difference (AD)    = ',F9.4)
180 format('# Root mean square difference (RMSD) = ',F9.4)

end subroutine

!===============================================================================

end program ffdev_compfreq_program

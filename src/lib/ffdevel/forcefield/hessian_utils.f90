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

module ffdev_hessian_utils

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_hessian_allocate
! ==============================================================================

subroutine ffdev_hessian_allocate(geo)

    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    if( geo%natoms .le. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Number of atoms has to be greater than zero for hessian!')
    end if

    allocate(geo%hess(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate array for hessian!')
    end if

    geo%hess(:,:,:,:) = 0.0d0

end subroutine ffdev_hessian_allocate

!===============================================================================
! logical function ffdev_hessian_test(geo1,geo2,alarm_treshold)
!===============================================================================

logical function ffdev_hessian_test(geo1,geo2,alarm_treshold)

    implicit none
    type(GEOMETRY)  :: geo1
    type(GEOMETRY)  :: geo2
    real(DEVDP)     :: alarm_treshold
    ! ------------------------------
    integer         :: i,j,k,l
    !---------------------------------------------------------------------------

    ffdev_hessian_test = .true.

    if( geo1%natoms .ne. geo2%natoms ) then
        ffdev_hessian_test = .false.
        return
    end if

    do i=1,geo1%natoms
        do j=1,3
            do k=1,geo1%natoms
                do l=1,3
                    if( abs(geo1%hess(j,i,l,k) - geo2%hess(j,i,l,k)) .gt. alarm_treshold ) then
                        ffdev_hessian_test = .false.
                        return
                    end if
                end do
            end do
        end do
    end do

end function ffdev_hessian_test

! ==============================================================================
! subroutine ffdev_hessian_print
! ==============================================================================

subroutine ffdev_hessian_print(fout,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer         :: fout
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,k,l
    ! --------------------------------------------------------------------------

    write(DEV_OUT,5) geo%natoms
    do i=1,geo%natoms
        write(DEV_OUT,10)
        write(DEV_OUT,20)
        do k=1,geo%natoms
            do l=1,3
                write(fout,30) k, char(119+l), i, geo%hess(l,k,1,i), geo%hess(l,k,2,i), geo%hess(l,k,3,i)
            end do
        end do
        if( i .ne. geo%natoms ) write(DEV_OUT,*)
    end do

 5 format('# FFDEVEL HESS ',I8)
10 format('# Idx  C   Idx          Hx                Hy                Hz       ')
20 format('# ---- - ------ ----------------- ----------------- -----------------')
30 format(I6,1X,A1,1X,I6,1X,F17.8,1X,F17.8,1X,F17.8)

end subroutine ffdev_hessian_print

! ==============================================================================
! subroutine ffdev_hessian_allocate_freq
! ==============================================================================

subroutine ffdev_hessian_allocate_freq(geo)

    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    if( associated(geo%nmodes) ) return     ! already allocated

    if( geo%natoms .le. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Number of atoms has to be greater than zero for frequencies!')
    end if

    allocate(geo%freq(3*geo%natoms), geo%nmodes(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arrays for frequencies!')
    end if

    geo%freq(:) = 0.0d0
    geo%nmodes(:,:,:,:) = 0.0d0

end subroutine ffdev_hessian_allocate_freq

! ==============================================================================
! subroutine ffdev_hessian_allocate_trg_hess
! ==============================================================================

subroutine ffdev_hessian_allocate_trg_hess(geo)

    use ffdev_geometry
    use ffdev_utils

    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    if( geo%natoms .le. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Number of atoms has to be greater than zero for hessian!')
    end if

    allocate(geo%trg_hess(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate array for hessian!')
    end if

    geo%trg_hess(:,:,:,:) = 0.0d0

end subroutine ffdev_hessian_allocate_trg_hess

! ==============================================================================
! subroutine ffdev_hessian_print_freqs
! ==============================================================================

subroutine ffdev_hessian_print_freqs(fout,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer         :: fout
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,j,k,l,m,r,h
    ! --------------------------------------------------------------------------

    write(fout,5) geo%natoms
    m = 0
    do while( m .lt. 3*geo%natoms )

        r = 3*geo%natoms - m
        if( r .gt. 6 ) r = 6

        if( m .ge. 6 ) write(fout,*)

        write(fout,10,ADVANCE='NO')
        do h=1,r
            write(fout,60,ADVANCE='NO') m+h
        end do
        write(fout,*)

        write(fout,20,ADVANCE='NO')
        do h=1,r
            write(fout,50,ADVANCE='NO') geo%freq(m+h)
        end do
        write(fout,*)

        write(fout,30)

        do k=1,geo%natoms
            do l=1,3
                write(fout,40,ADVANCE='NO') k, char(119+l)
                do h=1,r
                    i = (m+h-1) / 3 + 1
                    j = mod(m+h-1,3) + 1
                    write(fout,50,ADVANCE='NO') geo%nmodes(l,k,j,i)
                end do
                write(fout,*)
            end do
        end do
        m = m + r
    end do

 5 format('# FFDEVEL VIBS ',I8)
10 format('# Mode  ')
20 format('# Freq  ')
30 format('# ---- - ----------- ----------- ----------- ----------- ----------- -----------')
40 format(I6,1X,A1)
50 format(1X,F11.4)
60 format(3X,I6,3X)

end subroutine ffdev_hessian_print_freqs

! ==============================================================================
! subroutine ffdev_hessian_print_freqs_lin
! ==============================================================================

subroutine ffdev_hessian_print_freqs_lin(fout,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer         :: fout
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    write(fout,10) 3*geo%natoms
    write(fout,20)
    write(fout,30)

    do i=1, 3*geo%natoms
        write(fout,40) i,geo%freq(i)
    end do

10 format('# FFDEVEL VIBSONLY ',I8)
20 format('# Mode    Freq   ')
30 format('# ---- ----------')
40 format(I6,1X,F10.4)

end subroutine ffdev_hessian_print_freqs_lin

! ==============================================================================
! subroutine ffdev_hessian_calc_freqs
! ==============================================================================

subroutine ffdev_hessian_calc_freqs(geo)

    use ffdev_topology
    use ffdev_geometry
    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer                 :: i,j,k,l,info,lwork,alloc_stat
    real(DEVDP)             :: m1,m2,tmp(1),v,s
    real(DEVDP),allocatable :: work(:)
    ! --------------------------------------------------------------------------

    ! keep hessian as it is

    ! update hessian matrix by atom masses
    do i=1,geo%natoms
        m1 = sqrt(1.0d0 / pt_masses(geo%z(i)))
        do j=1,3
            do k=1,geo%natoms
                m2 = sqrt(1.0d0 / pt_masses(geo%z(k)))
                do l=1,3
                    geo%nmodes(j,i,l,k) =  geo%hess(j,i,l,k)*m1*m2
                end do
            end do
        end do
    end do

    ! get optimal work size
    lwork = -1
    call dsyev('V','L',geo%natoms*3, geo%nmodes, geo%natoms*3, geo%freq, tmp, lwork, info)
    lwork = tmp(1)

    ! allocate data
    allocate(work(lwork), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate array for force matrix diagonalization!')
    end if

    ! diagonalize matrix
    call dsyev('V','L',geo%natoms*3, geo%nmodes, geo%natoms*3, geo%freq, work, lwork, info)

    if( info .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Force matrix diagonalization failed!')
    end if

    deallocate(work)

    ! convert to frequencies
    do i=1,3*geo%natoms
        v = geo%freq(i)
        s = 1.0
        if( v .lt. 0.0 ) then
            s = -1.0
            v = abs(v)
        end if
        v = sqrt(v)*s
        geo%freq(i) = v*108.596538771 ! TODO - better number
    end do

end subroutine ffdev_hessian_calc_freqs

! ==============================================================================
! subroutine ffdev_hessian_calc_trg_freqs
! ==============================================================================

subroutine ffdev_hessian_calc_trg_freqs(geo)

    use ffdev_topology
    use ffdev_geometry
    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer                 :: i,j,k,l,info,lwork,alloc_stat
    real(DEVDP)             :: m1,m2,tmp(1),v,s
    real(DEVDP),allocatable :: work(:)
    ! --------------------------------------------------------------------------

    ! allocate data

    if( geo%natoms .le. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Number of atoms has to be greater than zero for target frequencies!')
    end if

    allocate(geo%trg_freq(3*geo%natoms), geo%trg_nmodes(3,geo%natoms,3,geo%natoms), &
             geo%freq_t2s_map(3*geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arrays for target frequencies!')
    end if

    geo%trg_freq(:) = 0.0d0
    geo%trg_nmodes(:,:,:,:) = 0.0d0

    ! keep hessian as it is

    ! update hessian matrix by atom masses
    do i=1,geo%natoms
        m1 = sqrt(1.0d0 / pt_masses(geo%z(i)))
        do j=1,3
            do k=1,geo%natoms
                m2 = sqrt(1.0d0 / pt_masses(geo%z(k)))
                do l=1,3
                    geo%trg_nmodes(j,i,l,k) =  geo%trg_hess(j,i,l,k)*m1*m2
                end do
            end do
        end do
    end do

    ! get optimal work size
    lwork = -1
    call dsyev('V','L',geo%natoms*3, geo%trg_nmodes, geo%natoms*3, geo%trg_freq, tmp, lwork, info)
    lwork = tmp(1)

    ! allocate data
    allocate(work(lwork), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate array for force matrix diagonalization!')
    end if

    ! diagonalize matrix
    call dsyev('V','L',geo%natoms*3, geo%trg_nmodes, geo%natoms*3, geo%trg_freq, work, lwork, info)

    if( info .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Force matrix diagonalization failed for traget hesssian!')
    end if

    deallocate(work)

    ! convert to frequencies
    do i=1,3*geo%natoms
        v = geo%trg_freq(i)
        s = 1.0
        if( v .lt. 0.0 ) then
            s = -1.0
            v = abs(v)
        end if
        v = sqrt(v)*s
        geo%trg_freq(i) = v*108.596538771 ! TODO - better number
    end do

end subroutine ffdev_hessian_calc_trg_freqs

!===============================================================================
! subroutine:  ffdev_hessian_superimpose_freqs
!===============================================================================

subroutine ffdev_hessian_superimpose_freqs(natoms,z,ref,src,srcnmodes,mw,rmsd)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    integer         :: natoms
    integer         :: z(:)
    real(DEVDP)     :: ref(:,:)
    real(DEVDP)     :: src(:,:)
    real(DEVDP)     :: srcnmodes(:,:,:,:)
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

    if( natoms .lt. 2 ) then
        call ffdev_utils_exit(DEV_OUT,1,'At least two atoms must be in the structures!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    x1 = 0.0d0
    x2 = 0.0d0
    x3 = 0.0d0
    xr1 = 0.0d0
    xr2 = 0.0d0
    xr3 = 0.0d0
    totmass = 0.0d0

    do  i = 1, natoms
        if( mw ) then
            amass = pt_masses(z(i))
        else
            amass = 1.0d0
        end if
        ! source
        x1 = x1 + src(1,i)*amass
        x2 = x2 + src(2,i)*amass
        x3 = x3 + src(3,i)*amass

        ! reference
        xr1 = xr1 + ref(1,i)*amass
        xr2 = xr2 + ref(2,i)*amass
        xr3 = xr3 + ref(3,i)*amass

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

    do i = 1, natoms
        if( mw ) then
            amass = pt_masses(z(i))
        else
            amass = 1.0d0
        end if

        x2sum = x2sum + amass*((src(1,i) - x1)**2 &
                             + (src(2,i) - x2)**2 &
                             + (src(3,i) - x3)**2)
        xr2sum = xr2sum + amass*((ref(1,i) - xr1)**2 &
                               + (ref(2,i) - xr2)**2 &
                               + (ref(3,i) - xr3)**2)

        r11 = r11 + amass*(src(1,i) - x1)*(ref(1,i) - xr1)
        r12 = r12 + amass*(src(1,i) - x1)*(ref(2,i) - xr2)
        r13 = r13 + amass*(src(1,i) - x1)*(ref(3,i) - xr3)

        r21 = r21 + amass*(src(2,i) - x2)*(ref(1,i) - xr1)
        r22 = r22 + amass*(src(2,i) - x2)*(ref(2,i) - xr2)
        r23 = r23 + amass*(src(2,i) - x2)*(ref(3,i) - xr3)

        r31 = r31 + amass*(src(3,i) - x3)*(ref(1,i) - xr1)
        r32 = r32 + amass*(src(3,i) - x3)*(ref(2,i) - xr2)
        r33 = r33 + amass*(src(3,i) - x3)*(ref(3,i) - xr3)
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
    do i=1,natoms
        ! move structure to origin
        px = src(1,i) - x1
        py = src(2,i) - x2
        pz = src(3,i) - x3

        ! rotate structure
        npx = u(1,1)*px + u(1,2)*py + u(1,3)*pz
        npy = u(2,1)*px + u(2,2)*py + u(2,3)*pz
        npz = u(3,1)*px + u(3,2)*py + u(3,3)*pz

        ! move to reference structure COM
        src(1,i) = npx + xr1
        src(2,i) = npy + xr2
        src(3,i) = npz + xr3

        do j=1,natoms
            do k=1,3
                ! it is a vector - no origin dependence
                px = srcnmodes(1,i,k,j)
                py = srcnmodes(2,i,k,j)
                pz = srcnmodes(3,i,k,j)

                ! rotate vector
                npx = u(1,1)*px + u(1,2)*py + u(1,3)*pz
                npy = u(2,1)*px + u(2,2)*py + u(2,3)*pz
                npz = u(3,1)*px + u(3,2)*py + u(3,3)*pz

                ! update
                srcnmodes(1,i,k,j) = npx
                srcnmodes(2,i,k,j) = npy
                srcnmodes(3,i,k,j) = npz
            end do
        end do
    end do

end subroutine ffdev_hessian_superimpose_freqs

!===============================================================================
! subroutine:  ffdev_hessian_find_mapping
!===============================================================================

subroutine ffdev_hessian_find_mapping(natoms,nmodes1,nmodes2,map1t2)

    use ffdev_utils

    implicit none
    integer             :: natoms
    real(DEVDP)         :: nmodes1(:,:,:,:)
    real(DEVDP)         :: nmodes2(:,:,:,:)
    integer             :: map1t2(:)
    ! -------------------------------------------
    integer             :: alloc_stat, i1, j1, i2, j2, max_j, k, l
    double precision    :: max_angle, angle
    integer,allocatable :: map1(:)
    ! --------------------------------------------------------------------------

    allocate(map1(3*natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate map array for freq mapping!')
    end if

    map1t2(:) = 0
    map1(:) = 0

    do i1=1,natoms
        do j1=1,3
            max_j = 0
            max_angle = 0.0
            do i2=1,natoms
                do j2=1,3
                    if( map1((i2-1)*3+j2) .gt. 0 ) cycle
                    ! calculate angle
                    angle = 0.0
                    do k=1,natoms
                        do l=1,3
                            angle = angle + nmodes1(l,k,j1,i1)*nmodes2(l,k,j2,i2)
                        end do
                    end do
                    if( (max_angle .lt. abs(angle)) .or. (max_j .eq. 0) ) then
                        max_angle = abs(angle)
                        max_j = (i2-1)*3+j2
                    end if
                end do
            end do
            map1(max_j) = (i1-1)*3+j1
            map1t2((i1-1)*3+j1) = max_j
        end do
    end do

    deallocate(map1)

end subroutine ffdev_hessian_find_mapping

!===============================================================================
! subroutine:  ffdev_hessian_print_mapping
!===============================================================================

subroutine ffdev_hessian_print_mapping(c12,natoms,freq1,nmodes1,freq2,nmodes2,map1t2)

    implicit none
    logical             :: c12      ! compare 1 vs 2, or TRG vs MM
    integer             :: natoms
    real(DEVDP)         :: freq1(:)
    real(DEVDP)         :: nmodes1(:,:,:,:)
    real(DEVDP)         :: freq2(:)
    real(DEVDP)         :: nmodes2(:,:,:,:)
    integer             :: map1t2(:)
    ! -------------------------------------------
    integer             :: i, i1, j1, i2, j2, k, l
    double precision    :: serr, lerr, aerr, rmse, angle, diff
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    if( c12 ) then
        write(DEV_OUT,120)
    else
        write(DEV_OUT,125)
    end if
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=7,3*natoms
        diff = freq2(map1t2(i))-freq1(i)

        i1 = (i-1) / 3 + 1
        j1 = mod(i-1,3) + 1

        i2 = (map1t2(i)-1) / 3 + 1
        j2 = mod(map1t2(i)-1,3) + 1

        angle = 0.0
        do k=1,natoms
            do l=1,3
                angle = angle + nmodes1(l,k,j1,i1)*nmodes2(l,k,j2,i2)
            end do
        end do
        angle = acos(angle)*DEV_R2D

        if( angle .gt. 90.0d0 ) then
            angle = 180.0 - angle
        end if

        write(DEV_OUT,140) i, freq1(i), map1t2(i), freq2(map1t2(i)), angle, diff
        if( serr .gt. abs(diff) ) serr = abs(diff)
        if( lerr .lt. abs(diff) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
    end do

    if( 3*natoms - 6 .gt. 0 ) then
        aerr = aerr / real(3*natoms - 6)
        rmse = sqrt(rmse / real(3*natoms - 6))
    end if

    write(DEV_OUT,130)
    write(DEV_OUT,150) serr
    write(DEV_OUT,160) lerr
    write(DEV_OUT,170) aerr
    write(DEV_OUT,180) rmse

100 format('# Individual vibrations')
120 format('# Indx#1  Freq#1  Indx#2  Freq#2    Angle   Diff')
125 format('# Id#TRG Freq#TRG  Id#MM  Freq#MM   Angle   Diff')
130 format('# ------ -------- ------ -------- ------- ------')
140 format(I8,1X,F8.2,1X,I6,1X,F8.2,1X,F7.1,1X,F6.1)
150 format('# Minimum unsigned difference (SUD)  = ',F9.4)
160 format('# Largest unsigned difference (MUD)  = ',F9.4)
170 format('# Average usigned difference (AD)    = ',F9.4)
180 format('# Root mean square difference (RMSD) = ',F9.4)

end subroutine

! ------------------------------------------------------------------------------

end module ffdev_hessian_utils

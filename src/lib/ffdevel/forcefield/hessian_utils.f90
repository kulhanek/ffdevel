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
use ffdev_variables

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
        call ffdev_utils_exit(DEV_ERR,1,'Number of atoms has to be greater than zero for hessian!')
    end if

    allocate(geo%hess(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for hessian!')
    end if

    geo%hess(:,:,:,:) = 0.0d0

end subroutine ffdev_hessian_allocate

! ==============================================================================
! subroutine ffdev_hessian_allocate_ihess
! ==============================================================================

subroutine ffdev_hessian_allocate_trg_ihess(top,geo)

    use ffdev_geometry
    use ffdev_utils
    use ffdev_topology

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    allocate(geo%trg_ihess(3,top%natoms,3,top%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array trg_ihess for hessian!')
    end if
    geo%trg_ihess = 0.0d0

    if( top%nbonds .gt. 0 ) then
        allocate(geo%trg_ihess_bonds(top%nbonds), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array trg_ihess_bonds for hessian!')
        end if
        geo%trg_ihess_bonds = 0.0d0
    end if

    if( top%nbonds .gt. 0 ) then
        allocate(geo%trg_ihess_angles(top%nangles), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array trg_ihess_angles for hessian!')
        end if
        geo%trg_ihess_angles = 0.0d0
    end if

end subroutine ffdev_hessian_allocate_trg_ihess

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
        call ffdev_utils_exit(DEV_ERR,1,'Number of atoms has to be greater than zero for frequencies!')
    end if

    allocate(geo%freq(3*geo%natoms), geo%nmodes(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate arrays for frequencies!')
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
        call ffdev_utils_exit(DEV_ERR,1,'Number of atoms has to be greater than zero for hessian!')
    end if

    allocate(geo%trg_hess(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for hessian!')
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
    lwork = INT(tmp(1))

    ! allocate data
    allocate(work(lwork), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for force matrix diagonalization!')
    end if

    ! diagonalize matrix
    call dsyev('V','L',geo%natoms*3, geo%nmodes, geo%natoms*3, geo%freq, work, lwork, info)

    if( info .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Force matrix diagonalization failed!')
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
! subroutine ffdev_hessian_calc_trg_ihess
! ==============================================================================

subroutine ffdev_hessian_calc_trg_ihess(top,geo)

    use ffdev_geometry
    use ffdev_utils
    use ffdev_topology
    use ffdev_jacobian

    implicit none
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    ! --------------------------------------------
    integer                 :: alloc_stat, n, i, m
    real(DEVDP),allocatable :: jac(:,:), ijac(:,:), ihess(:,:)
    ! --------------------------------------------------------------------------

    if( geo%natoms .le. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Number of atoms has to be greater than zero for ffdev_hessian_calc_trg_ihess!')
    end if

    m = top%natoms * 3
    n = top%nbonds + top%nangles

    if( n .le. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Number of bonds plus angles has to be greater than zero for ffdev_hessian_calc_trg_ihess!')
    end if

    ! DEBUG: write(*,*) 'm = ', m, ', n = ', n

    allocate(ihess(n,n), jac(m,n), ijac(n,m) , stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for ffdev_hessian_calc_trg_ihess!')
    end if

    if( associated(geo%trg_ihess_bonds) ) then
        geo%trg_ihess_bonds(:) = 0.0d0
    end if
    if( associated(geo%trg_ihess_angles) ) then
        geo%trg_ihess_angles(:) = 0.0d0
    end if

    ihess(:,:) = 0.0d0
    jac(:,:) = 0.0d0
    ijac(:,:) = 0.0d0

    ! calculate jacobian
    call ffdev_jacobian_calc_all(top,geo%crd,jac)
    call ffdev_jacobian_inverse(jac,ijac)
    call ffdev_jacobian_get_ihess(ijac,geo%trg_ihess,ihess)

    ! copy results to geo%trg_ihess_bonds and geo%trg_ihess_angles, only diagonal items
    do i=1,top%nbonds
        if( associated(geo%trg_ihess_bonds) ) geo%trg_ihess_bonds(i) = ihess(i,i)
    end do
    do i=1,top%nangles
        if( associated(geo%trg_ihess_angles) ) geo%trg_ihess_angles(i) = ihess(i+top%nbonds,i+top%nbonds)
    end do

    ! release working arrays
    deallocate(ijac,jac,ihess)

end subroutine ffdev_hessian_calc_trg_ihess

! ==============================================================================
! subroutine ffdev_hessian_print_trg_ihess_bonds
! ==============================================================================

subroutine ffdev_hessian_print_trg_ihess_bonds(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,ai,aj,bt
    real(DEVDP)     :: k1,k2,diff
    real(DEVDP)     :: serr,lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( top%nbonds .le. 0 ) return ! no data - exit

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    write(DEV_OUT,120)
    write(DEV_OUT,130)

    serr = 0.0d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=1,top%nbonds
        ai = top%bonds(i)%ai
        aj = top%bonds(i)%aj
        bt = top%bonds(i)%bt
        k1 = top%bond_types(bt)%k
        k2 = geo%trg_ihess_bonds(i)
        diff = k2 - k1
        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            k1,k2,diff
        if( (serr .gt. abs(diff)) .or. (i .eq. 1) ) serr = abs(diff)
        if( (lerr .lt. abs(diff)) .or. (i .eq. 1) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
    end do

    if( top%nbonds .gt. 0 ) then
        aerr = aerr / real(top%nbonds)
        rmse = sqrt(rmse / real(top%nbonds))
    end if

    write(DEV_OUT,110)
    write(DEV_OUT,150) serr
    write(DEV_OUT,160) lerr
    write(DEV_OUT,170) aerr
    write(DEV_OUT,180) rmse

100 format('# Individual bonds')
110 format('# --------------------------- = ----------------------------- -----------------------------')
120 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName K(top)    K(ihess)  diff(2-1)')
130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ ----- --------- --------- ---------')
140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.4,1X,F9.4,1X,F9.4)
150 format('# Minimum unsigned difference (SUD)  = ',F9.4)
160 format('# Largest unsigned difference (MUD)  = ',F9.4)
170 format('# Average unsigned difference (AD)   = ',F9.4)
180 format('# Root mean square difference (RMSD) = ',F9.4)


end subroutine ffdev_hessian_print_trg_ihess_bonds

!===============================================================================
! subroutine:  ffdev_hessian_print_trg_ihess_angles
!===============================================================================

subroutine ffdev_hessian_print_trg_ihess_angles(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i,ai,aj,ak,at
    real(DEVDP)     :: k1,k2,diff
    real(DEVDP)     :: serr,lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( top%nangles .le. 0 ) return ! no data - exit

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    write(DEV_OUT,120)
    write(DEV_OUT,130)

    serr = 0.0d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=1,top%nangles
        ai = top%angles(i)%ai
        aj = top%angles(i)%aj
        ak = top%angles(i)%ak
        at = top%angles(i)%at
        k1 = top%angle_types(at)%k
        k2 = geo%trg_ihess_angles(i)
        diff = k2 - k1
        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            ak, top%atoms(ak)%name, top%atom_types(top%atoms(ak)%typeid)%name, &
                            top%atoms(ak)%residx, top%atoms(ak)%resname, &
                            k1,k2,diff
        if( (serr .gt. abs(diff)) .or. (i .eq. 1) ) serr = abs(diff)
        if( (lerr .lt. abs(diff)) .or. (i .eq. 1) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
    end do

    if( top%nangles .gt. 0 ) then
        aerr = aerr / real(top%nangles)
        rmse = sqrt(rmse / real(top%nangles))
    end if

    write(DEV_OUT,110)
    write(DEV_OUT,150) serr
    write(DEV_OUT,160) lerr
    write(DEV_OUT,170) aerr
    write(DEV_OUT,180) rmse

100 format('# Individual angles')
110 format('# --------------------------- = ----------------------------- =&
           & ----------------------------- -----------------------------')
120 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName K(top)    K(ihess)   diff(2-1)')
130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ ----- =&
             & ------ ---- ---- ------ ----- --------- --------- ---------')
140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.2,1X,F9.2,1X,F9.2)
150 format('# Minimum unsigned difference (SUD)  = ',F9.2)
160 format('# Largest unsigned difference (MUD)  = ',F9.2)
170 format('# Average unsigned difference (AD)   = ',F9.2)
180 format('# Root mean square difference (RMSD) = ',F9.2)

end subroutine ffdev_hessian_print_trg_ihess_angles

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
        call ffdev_utils_exit(DEV_ERR,1,'Number of atoms has to be greater than zero for target frequencies!')
    end if

    allocate(geo%trg_freq(3*geo%natoms), geo%trg_nmodes(3,geo%natoms,3,geo%natoms), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate arrays for target frequencies!')
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
    lwork = INT(tmp(1))

    ! allocate data
    allocate(work(lwork), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate array for force matrix diagonalization!')
    end if

    ! diagonalize matrix
    call dsyev('V','L',geo%natoms*3, geo%trg_nmodes, geo%natoms*3, geo%trg_freq, work, lwork, info)

    if( info .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Force matrix diagonalization failed for traget hesssian!')
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

! ------------------------------------------------------------------------------

end module ffdev_hessian_utils

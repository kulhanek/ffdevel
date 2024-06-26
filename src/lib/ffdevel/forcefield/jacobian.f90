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

module ffdev_jacobian

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_jacobian_calc_all(top,crd,jac)
! ==============================================================================

subroutine ffdev_jacobian_calc_all(top,crd,jac)

    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd(:,:)
    real(DEVDP)     :: jac(:,:)
    ! --------------------------------------------
    integer         :: idx
    ! --------------------------------------------------------------------------

    jac(:,:) = 0.0d0

    idx = 0
    call ffdev_jacobian_bonds(top,crd,idx,jac)

    idx = idx + top%nbonds
    call ffdev_jacobian_angles(top,crd,idx,jac)

end subroutine ffdev_jacobian_calc_all

! ==============================================================================
! subroutine ffdev_jacobian_inverse(jac,ijac,fac)
! ==============================================================================

subroutine ffdev_jacobian_inverse(jac,ijac)

    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    real(DEVDP)     :: jac(:,:)
    real(DEVDP)     :: ijac(:,:)
    ! --------------------------------------------
    integer                 :: i, m, n, k, alloc_stat, info, lwork
    real(DEVDP),allocatable :: u(:,:), vt(:,:), sig(:), sig_plus(:,:)
    real(DEVDP),allocatable :: temp_mat(:,:), work(:)
    integer,allocatable     :: iwork(:)
    real(DEVDP)             :: fac
    ! --------------------------------------------------------------------------

    fac = 1d-3  ! FIXME - need to be tunable

    m = size(jac,1)
    n = size(jac,2)
    k = min(m,n)

    allocate(u(m,m),vt(n,n),sig(k),sig_plus(n,m),iwork(8*k),work(1),temp_mat(n,m), &
             stat = alloc_stat)
    if( alloc_stat .ne. 0) then
       call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate arrays I in ffdev_jacobian_inverse!')
    end if

    u(:,:)          = 0.0d0
    vt(:,:)         = 0.0d0
    sig(:)          = 0.0d0
    sig_plus(:,:)   = 0.0d0
    work(:)         = 0.0d0

    ! work size query
    lwork = -1
    call dgesdd('A', m, n, jac(1,1), m, sig(1), u(1,1), m, vt(1,1), n, work(1), &
                lwork, iwork(1), info)

    if( info .ne. 0) then
       call ffdev_utils_exit(DEV_ERR,1,'Unable to get size of working array in ffdev_jacobian_inverse!')
    end if

    ! reinit working array
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork), stat = alloc_stat)

    if( alloc_stat .ne. 0) then
       call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate arrays II in ffdev_jacobian_inverse!')
    end if

    ! do SVD
    call dgesdd('A', m, n, jac(1,1), m, sig(1), u(1,1), m, vt(1,1), n, work(1), &
                lwork, iwork(1), info)

    if( info .ne. 0) then
       call ffdev_utils_exit(DEV_ERR,1,'SVD failed in ffdev_jacobian_inverse!')
    end if

    ! set singular values that are too small to zero
    do i = 1, k
       if( sig(i) > fac*maxval(sig) ) then
          sig_plus(i,i) = 1.0d0/sig(i)
       else
          sig_plus(i,i) = 0.0d0
       end if
    end do

    ! build pseudoinverse: V*sig_plus*UT
    CALL dgemm('N', 'T', n, m, m, 1.0d0, sig_plus, n, u, m, 0.0d0, temp_mat, n)
    CALL dgemm('T', 'N', n, m, n, 1.0d0, vt, n, temp_mat, n, 0.0d0, ijac, n)

    ! clean data
    deallocate(u, vt, sig, iwork, work, sig_plus, temp_mat)

end subroutine ffdev_jacobian_inverse

! ==============================================================================
! subroutine ffdev_jacobian_get_ihess(ijac,hess,ihess)
! ==============================================================================

subroutine ffdev_jacobian_get_ihess(ijac,hess,ihess)

    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: ijac(:,:)
    real(DEVDP)             :: hess(:,:,:,:)
    real(DEVDP)             :: ihess(:,:)
    ! --------------------------------------------
    integer                 :: m,n,k,l,q,w,alloc_stat
    real(DEVDP),allocatable :: temp_mat(:,:)
    ! --------------------------------------------------------------------------

    m = size(hess,2)*3   ! 3n
    n = size(hess,4)*3   ! 3n
    k = size(ijac,1)     ! nint
    l = size(ijac,2)     ! 3n
    q = size(ihess,1)    ! nint
    w = size(ihess,2)    ! nint

    if( m .ne. n ) stop 'ffdev_jacobian_get_ihess 1'
    if( l .ne. m ) stop 'ffdev_jacobian_get_ihess 2'
    if( q .ne. w ) stop 'ffdev_jacobian_get_ihess 3'
    if( k .ne. q ) stop 'ffdev_jacobian_get_ihess 4'

    allocate(temp_mat(m,k), stat = alloc_stat)
    if( alloc_stat .ne. 0) then
       call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate arrays I in ffdev_jacobian_get_ihess!')
    end if

    CALL dgemm('N', 'T', m, k, n, 1.0d0, hess, m, ijac, k, 0.0d0, temp_mat, m)
    CALL dgemm('N', 'N', k, k, m, 1.0d0, ijac, k, temp_mat, m, 0.0d0, ihess, k)

    deallocate(temp_mat)

end subroutine ffdev_jacobian_get_ihess

! ==============================================================================
! subroutine ffdev_jacobian_bonds
! ==============================================================================

subroutine ffdev_jacobian_bonds(top,crd,idx,jac)

    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd(:,:)
    integer         :: idx
    real(DEVDP)     :: jac(:,:)
    ! --------------------------------------------
    integer         :: ib,i,j
    real(DEVDP)     :: rij(3),b,dv
    ! --------------------------------------------------------------------------

    ! bonds
    do ib=1,top%nbonds

        ! for each bond
        i  = top%bonds(ib)%ai
        j  = top%bonds(ib)%aj

        ! calculate gradient
        rij(:) = crd(:,j) - crd(:,i)
        b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
        dv = 1.0d0/b

        ! we assume no atom overlap
        jac(j*3-2,idx+ib) =  rij(1)*dv
        jac(j*3-1,idx+ib) =  rij(2)*dv
        jac(j*3-0,idx+ib) =  rij(3)*dv

        jac(i*3-2,idx+ib) = -rij(1)*dv
        jac(i*3-1,idx+ib) = -rij(2)*dv
        jac(i*3-0,idx+ib) = -rij(3)*dv
    end do

end subroutine ffdev_jacobian_bonds

! ==============================================================================
! subroutine ffdev_jacobian_angles
! ==============================================================================

subroutine ffdev_jacobian_angles(top,crd,idx,jac)

    use ffdev_topology_dat

    implicit none
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd(:,:)
    integer         :: idx
    real(DEVDP)     :: jac(:,:)
    ! --------------------------------------------
    integer         :: ia,i,j,k
    real(DEVDP)     :: scp,f1,bji2inv,bjk2inv,bjiinv,bjkinv,angv
    real(DEVDP)     :: rji(3),rjk(3),di(3),dk(3)
    ! --------------------------------------------------------------------------

    do ia=1,top%nangles

        i  = top%angles(ia)%ai
        j  = top%angles(ia)%aj
        k  = top%angles(ia)%ak

        ! calculate rji and rjk
        rji(:) = crd(:,i) - crd(:,j)
        rjk(:) = crd(:,k) - crd(:,j)

        ! calculate bjiinv and bjkinv and their squares
        bji2inv = 1.0d0/(rji(1)**2 + rji(2)**2 + rji(3)**2 )
        bjk2inv = 1.0d0/(rjk(1)**2 + rjk(2)**2 + rjk(3)**2 )
        bjiinv = sqrt(bji2inv)
        bjkinv = sqrt(bjk2inv)

        ! calculate scp and angv
        scp = ( rji(1)*rjk(1) + rji(2)*rjk(2) + rji(3)*rjk(3) )
        scp = scp * bjiinv*bjkinv
        if ( scp .gt.  1.0d0 ) then
            scp =  1.0d0
        else if ( scp .lt. -1.0d0 ) then
            scp = -1.0d0
        end if
        angv = acos(scp)

        ! gradient
        f1 = sin ( angv )
        if ( abs(f1) .lt. 1.0d-3 ) then
            ! sin(0.1 deg) = 1.7e-3
            ! this is set for angles close to 0 deg or 180 deg by 0.1 deg
            ! the aim is to avoid division be zero
            f1 = -1.0d3
        else
            f1 = -1.0d0 / f1
        end if

        di(:) = f1 * ( rjk(:)*bjiinv*bjkinv - scp*rji(:)*bji2inv )
        dk(:) = f1 * ( rji(:)*bjiinv*bjkinv - scp*rjk(:)*bjk2inv )

        jac(i*3-2,idx+ia) = + di(1)
        jac(i*3-1,idx+ia) = + di(2)
        jac(i*3-0,idx+ia) = + di(3)

        jac(j*3-2,idx+ia) = - ( di(1) + dk(1) )
        jac(j*3-1,idx+ia) = - ( di(2) + dk(2) )
        jac(j*3-0,idx+ia) = - ( di(3) + dk(3) )

        jac(k*3-2,idx+ia) = + dk(1)
        jac(k*3-1,idx+ia) = + dk(2)
        jac(k*3-0,idx+ia) = + dk(3)

    end do

end subroutine ffdev_jacobian_angles

! ------------------------------------------------------------------------------

end module ffdev_jacobian

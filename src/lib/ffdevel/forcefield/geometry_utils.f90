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

module ffdev_geometry_utils

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_geometry_utils_get_rmsd
! ==============================================================================

real(DEVDP) function ffdev_geometry_utils_get_rmsd(ref,src,mw)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: ref
    type(GEOMETRY)  :: src
    logical         :: mw       ! mass weighted
    ! --------------------------------------------
    integer        :: i,info,best
    real(DEVDP)    :: x1,x2,x3,xr1,xr2,xr3,amass,totmass,itotmass
    real(DEVDP)    :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(DEVDP)    :: f(4,4),u(3,3)
    real(DEVDP)    :: eigenvalues(4),work(26*4)
    real(DEVDP)    :: x2sum,xr2sum
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
    ffdev_geometry_utils_get_rmsd = sqrt((x2sum + xr2sum - 2.0d0*eigenvalues(best))*itotmass)

end function ffdev_geometry_utils_get_rmsd

! ==============================================================================
! subroutine ffdev_geometry_utils_get_rmsd_nofit
! ==============================================================================

real(DEVDP) function ffdev_geometry_utils_get_rmsd_nofit(ref,src,mw)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: ref
    type(GEOMETRY)  :: src
    logical         :: mw       ! mass weighted
    ! --------------------------------------------
    integer        :: i
    real(DEVDP)    :: totmass, value, dv, amass
    ! --------------------------------------------------------------------------

    if( src%natoms .ne. ref%natoms ) then
        call ffdev_utils_exit(DEV_OUT,1,'Two geometries must have the same number of atoms!')
    end if
    if( src%natoms .lt. 1 ) then
        call ffdev_utils_exit(DEV_OUT,1,'At least one atom must be in both structures!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    totmass = 0.0d0
    value = 0.0d0
    do  i = 1, src%natoms
        if( src%z(i) .ne. ref%z(i) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Two geometries must have identical order of atoms!')
        end if
        if( mw ) then
            amass = pt_masses(src%z(i))
        else
            amass = 1.0d0
        end if
        dv = ( src%crd(1,i) - ref%crd(1,i) )**2 + ( src%crd(2,i) - ref%crd(2,i) )**2 &
           + ( src%crd(3,i) - ref%crd(3,i) )**2
        value = value + amass*dv
        totmass = totmass + amass
    end do

    ffdev_geometry_utils_get_rmsd_nofit = sqrt(value / totmass)

end function ffdev_geometry_utils_get_rmsd_nofit

! ==============================================================================
! subroutine ffdev_geometry_utils_rmsdfit
! ==============================================================================

subroutine ffdev_geometry_utils_rmsdfit(ref,src,mw,rmsd)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: ref
    type(GEOMETRY)  :: src
    logical         :: mw       ! mass weighted
    real(DEVDP)     :: rmsd
    ! --------------------------------------------
    integer         :: i,info,best
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
    end do

end subroutine ffdev_geometry_utils_rmsdfit

! ------------------------------------------------------------------------------

end module ffdev_geometry_utils

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

real(DEVDP) function ffdev_geometry_utils_get_rmsd(natoms,z,ref,src,mw)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    integer        :: natoms
    integer        :: z(:)
    real(DEVDP)    :: ref(:,:)
    real(DEVDP)    :: src(:,:)
    logical        :: mw       ! mass weighted
    ! --------------------------------------------
    integer        :: i,info,best
    real(DEVDP)    :: x1,x2,x3,xr1,xr2,xr3,amass,totmass,itotmass
    real(DEVDP)    :: r11,r12,r13,r21,r22,r23,r31,r32,r33
    real(DEVDP)    :: f(4,4),u(3,3)
    real(DEVDP)    :: eigenvalues(4),work(26*4)
    real(DEVDP)    :: x2sum,xr2sum
    ! --------------------------------------------------------------------------

    if( natoms .lt. 2 ) then
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
    ffdev_geometry_utils_get_rmsd = sqrt((x2sum + xr2sum - 2.0d0*eigenvalues(best))*itotmass)

end function ffdev_geometry_utils_get_rmsd

! ==============================================================================
! subroutine ffdev_geometry_utils_get_rmsd_nofit
! ==============================================================================

real(DEVDP) function ffdev_geometry_utils_get_rmsd_nofit(natoms,z,ref,src,mw)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    integer         :: natoms
    integer         :: z(:)
    real(DEVDP)     :: ref(:,:)
    real(DEVDP)     :: src(:,:)
    logical         :: mw       ! mass weighted
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: totmass, value, dv, amass
    ! --------------------------------------------------------------------------

    if( natoms .lt. 1 ) then
        call ffdev_utils_exit(DEV_OUT,1,'At least one atom must be in both structures!')
    end if

    ! calculate geometrical centres (source and target) -------------------
    totmass = 0.0d0
    value = 0.0d0
    do  i = 1, natoms
        if( mw ) then
            amass = pt_masses(z(i))
        else
            amass = 1.0d0
        end if
        dv = ( src(1,i) - ref(1,i) )**2 + ( src(2,i) - ref(2,i) )**2 &
           + ( src(3,i) - ref(3,i) )**2
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

! ==============================================================================
! subroutine ffdev_geometry_utils_comp_bonds
! ==============================================================================

subroutine ffdev_geometry_utils_comp_bonds(c12,top,crd1,crd2)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    logical         :: c12
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd1(:,:)
    real(DEVDP)     :: crd2(:,:)
    ! --------------------------------------------
    integer         :: i, j, ai, aj, nb
    real(DEVDP)     :: d1, d2, diff
    real(DEVDP)     :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( top%nbonds .le. 0 ) return ! no data - exit

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
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

    do i=1,top%nbonds
        ai = top%bonds(i)%ai
        aj = top%bonds(i)%aj
        d1 = ffdev_geometry_get_length(crd1,ai,aj)
        d2 = ffdev_geometry_get_length(crd2,ai,aj)
        diff = d2 - d1
        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            d1,d2,diff
        if( serr .gt. abs(diff) ) serr = abs(diff)
        if( lerr .lt. abs(diff) ) lerr = abs(diff)
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
120 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName    d#1       d#2    diff(2-1)')
125 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName  d#TRG(1)   d#MM(2) diff(2-1)')
130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ ----- --------- --------- ---------')
140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.4,1X,F9.4,1X,F9.4)
150 format('# Minimum unsigned difference (SUD)  = ',F9.4)
160 format('# Largest unsigned difference (MUD)  = ',F9.4)
170 format('# Average usigned difference (AD)    = ',F9.4)
180 format('# Root mean square difference (RMSD) = ',F9.4)

    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,220)
    write(DEV_OUT,230)

    do i=1,top%nbond_types

        serr = 100d0
        lerr = 0.0d0
        aerr = 0.0d0
        rmse = 0.0d0
        nb = 0

        do j=1,top%nbonds
            if( top%bonds(j)%bt .ne. i ) cycle

            ai = top%bonds(j)%ai
            aj = top%bonds(j)%aj
            d1 = ffdev_geometry_get_length(crd1,ai,aj)
            d2 = ffdev_geometry_get_length(crd2,ai,aj)
            diff = d2 - d1

            if( serr .gt. abs(diff) ) serr = abs(diff)
            if( lerr .lt. abs(diff) ) lerr = abs(diff)
            aerr = aerr + abs(diff)
            rmse = rmse + diff**2
            nb = nb + 1
        end do

        if( nb .gt. 0 ) then
            aerr = aerr / real(nb)
            rmse = sqrt(rmse / real(nb))
        end if

        write(DEV_OUT,240) top%atom_types(top%bond_types(i)%ti)%name, &
                           top%atom_types(top%bond_types(i)%tj)%name, serr, lerr, aerr, rmse

    end do

200 format('# Bonds by types')
210 format('# ---------------------------------------------------')
220 format('# Type   Type    SUD       MUD       AD        RMSD  ')
230 format('# ---- = ---- --------- --------- --------- ---------')
240 format(2X,A4,3X,A4,1X,F9.4,1X,F9.4,1X,F9.4,1X,F9.4)

end subroutine ffdev_geometry_utils_comp_bonds

!===============================================================================
! subroutine:  ffdev_geometry_utils_comp_angles
!===============================================================================

subroutine ffdev_geometry_utils_comp_angles(c12,top,crd1,crd2)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    logical         :: c12
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd1(:,:)
    real(DEVDP)     :: crd2(:,:)
    ! --------------------------------------------
    integer         :: i, j, ai, aj, ak, nb
    real(DEVDP)     :: d1, d2, diff
    real(DEVDP)     :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( top%nangles .le. 0 ) return ! no data - exit

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
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

    do i=1,top%nangles
        ai = top%angles(i)%ai
        aj = top%angles(i)%aj
        ak = top%angles(i)%ak
        d1 = ffdev_geometry_get_angle(crd1,ai,aj,ak)
        d2 = ffdev_geometry_get_angle(crd2,ai,aj,ak)
        diff = d2 - d1
        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            ak, top%atoms(ak)%name, top%atom_types(top%atoms(ak)%typeid)%name, &
                            top%atoms(ak)%residx, top%atoms(ak)%resname, &
                            d1*DEV_R2D,d2*DEV_R2D,diff*DEV_R2D
        if( serr .gt. abs(diff) ) serr = abs(diff)
        if( lerr .lt. abs(diff) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
    end do

    if( top%nangles .gt. 0 ) then
        aerr = aerr / real(top%nangles)
        rmse = sqrt(rmse / real(top%nangles))
    end if

    write(DEV_OUT,110)
    write(DEV_OUT,150) serr*DEV_R2D
    write(DEV_OUT,160) lerr*DEV_R2D
    write(DEV_OUT,170) aerr*DEV_R2D
    write(DEV_OUT,180) rmse*DEV_R2D

100 format('# Individual angles')
110 format('# --------------------------- = ----------------------------- =&
           & ----------------------------- -----------------------------')
120 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName    a#1       a#2    diff(2-1)')
125 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName  a#TRG(1)   a#MM(2) diff(2-1)')
130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ ----- =&
             & ------ ---- ---- ------ ----- --------- --------- ---------')
140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.2,1X,F9.2,1X,F9.2)
150 format('# Minimum unsigned difference (SUD)  = ',F9.2)
160 format('# Largest unsigned difference (MUD)  = ',F9.2)
170 format('# Average usigned difference (AD)    = ',F9.2)
180 format('# Root mean square difference (RMSD) = ',F9.2)


    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,220)
    write(DEV_OUT,230)

    do i=1,top%nangle_types

        serr = 100d0
        lerr = 0.0d0
        aerr = 0.0d0
        rmse = 0.0d0
        nb = 0

        do j=1,top%nangles
            if( top%angles(j)%at .ne. i ) cycle

            ai = top%angles(j)%ai
            aj = top%angles(j)%aj
            ak = top%angles(j)%ak
            d1 = ffdev_geometry_get_angle(crd1,ai,aj,ak)
            d2 = ffdev_geometry_get_angle(crd2,ai,aj,ak)
            diff = d2 - d1
            if( serr .gt. abs(diff) ) serr = abs(diff)
            if( lerr .lt. abs(diff) ) lerr = abs(diff)
            aerr = aerr + abs(diff)
            rmse = rmse + diff**2
            nb = nb + 1
        end do

        if( nb .gt. 0 ) then
            aerr = aerr / real(nb)
            rmse = sqrt(rmse / real(nb))
        end if

        write(DEV_OUT,240) top%atom_types(top%angle_types(i)%ti)%name, &
                           top%atom_types(top%angle_types(i)%tj)%name, &
                           top%atom_types(top%angle_types(i)%tk)%name, &
                           nb, serr*DEV_R2D, lerr*DEV_R2D, aerr*DEV_R2D, rmse*DEV_R2D
    end do

200 format('# Angles by types')
210 format('# ----------------------------------------------------------------')
220 format('# Type   Type   Type Count    SUD       MUD       AD        RMSD  ')
230 format('# ---- = ---- = ---- ----- --------- --------- --------- ---------')
240 format(2X,A4,3X,A4,3X,A4,1X,I5,1X,F9.2,1X,F9.2,1X,F9.2,1X,F9.2)

end subroutine ffdev_geometry_utils_comp_angles

!===============================================================================
! subroutine:  ffdev_geometry_utils_comp_dihedrals
!===============================================================================

subroutine ffdev_geometry_utils_comp_dihedrals(c12,top,crd1,crd2,onlytyped)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_parameters_dat

    implicit none
    logical         :: c12
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd1(:,:)
    real(DEVDP)     :: crd2(:,:)
    logical         :: onlytyped
    ! --------------------------------------------
    integer         :: i, j, ai, aj, ak, al, nb
    logical         :: found
    real(DEVDP)     :: d1, d2, diff
    real(DEVDP)     :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( top%ndihedrals .le. 0 ) return ! no data - exit

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
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
    nb   = 0

    do i=1,top%ndihedrals
        ai = top%dihedrals(i)%ai
        aj = top%dihedrals(i)%aj
        ak = top%dihedrals(i)%ak
        al = top%dihedrals(i)%al

        if( onlytyped ) then
            if( .not. top%dihedral_types(top%dihedrals(i)%dt)%ffoptactive ) cycle
        end if

        d1 = ffdev_geometry_get_dihedral(crd1,ai,aj,ak,al)
        d2 = ffdev_geometry_get_dihedral(crd2,ai,aj,ak,al)
        diff = ffdev_geometry_get_dihedral_deviation(d2,d1)
        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            ak, top%atoms(ak)%name, top%atom_types(top%atoms(ak)%typeid)%name, &
                            top%atoms(ak)%residx, top%atoms(ak)%resname, &
                            al, top%atoms(al)%name, top%atom_types(top%atoms(al)%typeid)%name, &
                            top%atoms(al)%residx, top%atoms(al)%resname, &
                            d1*DEV_R2D,d2*DEV_R2D,diff*DEV_R2D
        if( serr .gt. abs(diff) ) serr = abs(diff)
        if( lerr .lt. abs(diff) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
        nb = nb + 1
    end do

    if( nb .gt. 0 ) then
        aerr = aerr / real(nb)
        rmse = sqrt(rmse / real(nb))
    end if

    write(DEV_OUT,110)
    write(DEV_OUT,150) serr*DEV_R2D
    write(DEV_OUT,160) lerr*DEV_R2D
    write(DEV_OUT,170) aerr*DEV_R2D
    write(DEV_OUT,180) rmse*DEV_R2D

100 format('# Individual dihedrals')
110 format('# --------------------------- = ----------------------------- = ----------------------------- =&
           & ----------------------------- -----------------------------')
120 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName    d#1       d#2    diff(2-1)')
125 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName  d#TRG(1)   d#MM(2) diff(2-1)')
130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ -----&
           & = ------ ---- ---- ------ ----- =&
             & ------ ---- ---- ------ ----- --------- --------- ---------')
140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,&
           A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.2,1X,F9.2,1X,F9.2)
150 format('# Minimum unsigned difference (SUD)  = ',F9.2)
160 format('# Largest unsigned difference (MUD)  = ',F9.2)
170 format('# Average usigned difference (AD)    = ',F9.2)
180 format('# Root mean square difference (RMSD) = ',F9.2)


    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,220)
    write(DEV_OUT,230)

    do i=1,top%ndihedral_types

        if( onlytyped ) then
            if( .not. top%dihedral_types(i)%ffoptactive ) cycle
        end if

        serr = 100d0
        lerr = 0.0d0
        aerr = 0.0d0
        rmse = 0.0d0
        nb = 0

        do j=1,top%ndihedrals
            if( top%dihedrals(j)%dt .ne. i ) cycle

            ai = top%dihedrals(j)%ai
            aj = top%dihedrals(j)%aj
            ak = top%dihedrals(j)%ak
            al = top%dihedrals(j)%al
            d1 = ffdev_geometry_get_dihedral(crd1,ai,aj,ak,al)
            d2 = ffdev_geometry_get_dihedral(crd2,ai,aj,ak,al)
            diff = ffdev_geometry_get_dihedral_deviation(d2,d1)
            if( serr .gt. abs(diff) ) serr = abs(diff)
            if( lerr .lt. abs(diff) ) lerr = abs(diff)
            aerr = aerr + abs(diff)
            rmse = rmse + diff**2
            nb = nb + 1
        end do

        if( nb .gt. 0 ) then
            aerr = aerr / real(nb)
            rmse = sqrt(rmse / real(nb))
        end if

        write(DEV_OUT,240) top%atom_types(top%dihedral_types(i)%ti)%name, &
                           top%atom_types(top%dihedral_types(i)%tj)%name, &
                           top%atom_types(top%dihedral_types(i)%tk)%name, &
                           top%atom_types(top%dihedral_types(i)%tl)%name, &
                           nb, serr*DEV_R2D, lerr*DEV_R2D, aerr*DEV_R2D, rmse*DEV_R2D
    end do

200 format('# Dihedrals by types')
210 format('# -----------------------------------------------------------------------')
220 format('# Type   Type   Type   Type Count    SUD       MUD       AD        RMSD  ')
230 format('# ---- = ---- = ---- = ---- ----- --------- --------- --------- ---------')
240 format(2X,A4,3X,A4,3X,A4,3X,A4,1X,I5,1X,F9.2,1X,F9.2,1X,F9.2,1X,F9.2)

end subroutine ffdev_geometry_utils_comp_dihedrals

!===============================================================================
! subroutine:  ffdev_geometry_utils_comp_impropers
!===============================================================================

subroutine ffdev_geometry_utils_comp_impropers(c12,top,crd1,crd2,lock2phase)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    logical         :: c12
    type(TOPOLOGY)  :: top
    real(DEVDP)     :: crd1(:,:)
    real(DEVDP)     :: crd2(:,:)
    logical         :: lock2phase
    ! --------------------------------------------
    integer         :: i, j, ai, aj, ak, al, nb, idt
    real(DEVDP)     :: d1, d2, diff
    real(DEVDP)     :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( top%nimpropers .le. 0 ) return ! no data - exit

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    if( c12 ) then
        write(DEV_OUT,120)
    else
        if( lock2phase ) then
            write(DEV_OUT,126)
        else
            write(DEV_OUT,125)
        end if
    end if
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=1,top%nimpropers
        ai = top%impropers(i)%ai
        aj = top%impropers(i)%aj
        ak = top%impropers(i)%ak
        al = top%impropers(i)%al
        if( lock2phase ) then
            idt = top%impropers(i)%dt
            d1 = top%improper_types(idt)%g
        else
            d1 = ffdev_geometry_get_improper(crd1,ai,aj,ak,al)
        end if
        d2 = ffdev_geometry_get_improper(crd2,ai,aj,ak,al)
        diff = ffdev_geometry_get_dihedral_deviation(d2,d1)
        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            ak, top%atoms(ak)%name, top%atom_types(top%atoms(ak)%typeid)%name, &
                            top%atoms(ak)%residx, top%atoms(ak)%resname, &
                            al, top%atoms(al)%name, top%atom_types(top%atoms(al)%typeid)%name, &
                            top%atoms(al)%residx, top%atoms(al)%resname, &
                            d1*DEV_R2D,d2*DEV_R2D,diff*DEV_R2D
        if( serr .gt. abs(diff) ) serr = abs(diff)
        if( lerr .lt. abs(diff) ) lerr = abs(diff)
        aerr = aerr + abs(diff)
        rmse = rmse + diff**2
    end do

    if( top%nimpropers .gt. 0 ) then
        aerr = aerr / real(top%nimpropers)
        rmse = sqrt(rmse / real(top%nimpropers))
    end if

    write(DEV_OUT,110)
    write(DEV_OUT,150) serr*DEV_R2D
    write(DEV_OUT,160) lerr*DEV_R2D
    write(DEV_OUT,170) aerr*DEV_R2D
    write(DEV_OUT,180) rmse*DEV_R2D

100 format('# Individual impropers')
110 format('# --------------------------- = ----------------------------- = ----------------------------- =&
           & ----------------------------- -----------------------------')
120 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName    d#1       d#2    diff(2-1)')
125 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
           &    Indx  Name Type  RIdx  RName  d#TRG(1)   d#MM(2) diff(2-1)')
126 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName    Indx  Name Type  RIdx  RName&
          &    Indx  Name Type  RIdx  RName  d#TOP(1)   d#MM(2) diff(2-1)')
130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ -----&
           & = ------ ---- ---- ------ ----- =&
             & ------ ---- ---- ------ ----- --------- --------- ---------')
140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,&
           A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.2,1X,F9.2,1X,F9.2)
150 format('# Minimum unsigned difference (SUD)  = ',F9.2)
160 format('# Largest unsigned difference (MUD)  = ',F9.2)
170 format('# Average usigned difference (AD)    = ',F9.2)
180 format('# Root mean square difference (RMSD) = ',F9.2)


    write(DEV_OUT,*)
    write(DEV_OUT,200)
    write(DEV_OUT,210)
    write(DEV_OUT,220)
    write(DEV_OUT,230)

    do i=1,top%nimproper_types

        serr = 100d0
        lerr = 0.0d0
        aerr = 0.0d0
        rmse = 0.0d0
        nb = 0

        do j=1,top%nimpropers
            if( top%impropers(j)%dt .ne. i ) cycle

            ai = top%impropers(j)%ai
            aj = top%impropers(j)%aj
            ak = top%impropers(j)%ak
            al = top%impropers(j)%al
            if( lock2phase ) then
                idt = top%impropers(j)%dt
                d1 = top%improper_types(idt)%g
            else
                d1 = ffdev_geometry_get_improper(crd1,ai,aj,ak,al)
            end if
            d2 = ffdev_geometry_get_improper(crd2,ai,aj,ak,al)
            diff = ffdev_geometry_get_dihedral_deviation(d2,d1)
            if( serr .gt. abs(diff) ) serr = abs(diff)
            if( lerr .lt. abs(diff) ) lerr = abs(diff)
            aerr = aerr + abs(diff)
            rmse = rmse + diff**2
            nb = nb + 1
        end do

        if( nb .gt. 0 ) then
            aerr = aerr / real(nb)
            rmse = sqrt(rmse / real(nb))
        end if

        write(DEV_OUT,240) top%atom_types(top%improper_types(i)%ti)%name, &
                           top%atom_types(top%improper_types(i)%tj)%name, &
                           top%atom_types(top%improper_types(i)%tk)%name, &
                           top%atom_types(top%improper_types(i)%tl)%name, &
                           nb, serr*DEV_R2D, lerr*DEV_R2D, aerr*DEV_R2D, rmse*DEV_R2D
    end do

200 format('# Impropers by types')
210 format('# -----------------------------------------------------------------------')
220 format('# Type   Type   Type   Type Count    SUD       MUD       AD        RMSD  ')
230 format('# ---- = ---- = ---- = ---- ----- --------- --------- --------- ---------')
240 format(2X,A4,3X,A4,3X,A4,3X,A4,1X,I5,1X,F9.2,1X,F9.2,1X,F9.2,1X,F9.2)

end subroutine ffdev_geometry_utils_comp_impropers

! ------------------------------------------------------------------------------

end module ffdev_geometry_utils

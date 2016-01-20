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

program ffdev_compgeo_program

    use ffdev_sizes
    use ffdev_utils
    use ffdev_constants
    use ffdev_topology
    use ffdev_geometry
    use ffdev_energy

    implicit none
    character(len=MAX_PATH) :: topname      ! input topology name
    character(len=MAX_PATH) :: crdname1     ! input coordinate name #1
    character(len=MAX_PATH) :: crdname2     ! input coordinate name #2
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo1
    type(GEOMETRY)          :: geo2
    ! --------------------------------------------------------------------------

    call ffdev_utils_header('Compare Geometry')

    ! test number of input arguments
    if( command_argument_count() .ne. 3 ) then
        call print_usage()
        call ffdev_utils_exit(DEV_OUT,1,'Incorrect number of arguments was specified (three expected)!')
    end if

    call get_command_argument(1, topname)
    call get_command_argument(2, crdname1)
    call get_command_argument(3, crdname2)

    ! paths
    write(DEV_OUT,*)
    write(DEV_OUT,100) trim(topname)
    write(DEV_OUT,110) trim(crdname1)
    write(DEV_OUT,120) trim(crdname2)

    ! load topology
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Simplified Topology','=')

    call ffdev_topology_init(top)
    call ffdev_topology_load(top,topname)
    call ffdev_topology_info(top)

    ! load coordinates #1
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates #1','=')

    call ffdev_geometry_init(geo1)
    call ffdev_geometry_load_xyz(geo1,crdname1)
    call ffdev_geometry_info_input(geo1)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo1)

    ! load coordinates #2
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'XYZ Coordinates #2','=')

    call ffdev_geometry_init(geo2)
    call ffdev_geometry_load_xyz(geo2,crdname2)
    call ffdev_geometry_info_input(geo2)

    ! check coordinates and topology
    call ffdev_geometry_check_z(top,geo2)

    ! finalize topology setup and calculate energy
    call ffdev_topology_finalize_setup(top)

    ! compare geometry
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'BONDS','=')
    call compare_bonds
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'ANGLES','=')
    call compare_angles
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'DIHEDRALS','=')
    call compare_dihedrals

    call ffdev_utils_footer('Compare Geometry')

100 format('Simplified topology : ',A)
110 format('XYZ coordinates #1  : ',A)
120 format('XYZ coordinates #2  : ',A)

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

10 format('    compgeo <stopology> <xyzcrd#1> <xyzcrd#2>')

end subroutine print_usage

!===============================================================================
! subroutine:  compare_bonds
!===============================================================================

subroutine compare_bonds()

    implicit none
    integer     :: i, j, ai, aj, nb
    real(DEVDP) :: d1, d2, diff
    real(DEVDP) :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    write(DEV_OUT,120)
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=1,top%nbonds
        ai = top%bonds(i)%ai
        aj = top%bonds(i)%aj
        d1 = ffdev_geometry_get_length(geo1,ai,aj)
        d2 = ffdev_geometry_get_length(geo2,ai,aj)
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
            d1 = ffdev_geometry_get_length(geo1,ai,aj)
            d2 = ffdev_geometry_get_length(geo2,ai,aj)
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

end subroutine compare_bonds

!===============================================================================
! subroutine:  compare_angles
!===============================================================================

subroutine compare_angles()

    implicit none
    integer     :: i, j, ai, aj, ak, nb
    real(DEVDP) :: d1, d2, diff
    real(DEVDP) :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    write(DEV_OUT,120)
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=1,top%nangles
        ai = top%angles(i)%ai
        aj = top%angles(i)%aj
        ak = top%angles(i)%ak
        d1 = ffdev_geometry_get_angle(geo1,ai,aj,ak)
        d2 = ffdev_geometry_get_angle(geo2,ai,aj,ak)
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
            d1 = ffdev_geometry_get_angle(geo1,ai,aj,ak)
            d2 = ffdev_geometry_get_angle(geo2,ai,aj,ak)
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

end subroutine compare_angles

!===============================================================================
! subroutine:  compare_dihedrals
!===============================================================================

subroutine compare_dihedrals()

    implicit none
    integer     :: i, j, ai, aj, ak, al, nb
    real(DEVDP) :: d1, d2, diff
    real(DEVDP) :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    write(DEV_OUT,120)
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0

    do i=1,top%ndihedrals
        ai = top%dihedrals(i)%ai
        aj = top%dihedrals(i)%aj
        ak = top%dihedrals(i)%ak
        al = top%dihedrals(i)%al
        d1 = ffdev_geometry_get_dihedral(geo1,ai,aj,ak,al)
        d2 = ffdev_geometry_get_dihedral(geo2,ai,aj,ak,al)
        diff = d2 - d1
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

    if( top%ndihedrals .gt. 0 ) then
        aerr = aerr / real(top%ndihedrals)
        rmse = sqrt(rmse / real(top%ndihedrals))
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
            d1 = ffdev_geometry_get_dihedral(geo1,ai,aj,ak,al)
            d2 = ffdev_geometry_get_dihedral(geo2,ai,aj,ak,al)
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

end subroutine compare_dihedrals

!===============================================================================

end program ffdev_compgeo_program

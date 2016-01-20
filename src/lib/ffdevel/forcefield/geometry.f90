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

module ffdev_geometry

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_geometry_init
! ==============================================================================

subroutine ffdev_geometry_init(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    geo%id = 0
    geo%name = ''
    geo%title = ''
    geo%natoms = 0
    geo%bond_ene = 0
    geo%angle_ene = 0
    geo%dih_ene = 0
    geo%impropr_ene = 0
    geo%ele14_ene = 0
    geo%nb14_ene = 0
    geo%ele_ene = 0
    geo%nb_ene = 0
    geo%weight = 0
    geo%trg_energy = 0
    geo%total_ene = 0
    geo%weight = 1.0

    geo%trg_ene_loaded = .false.
    geo%trg_grd_loaded = .false.
    geo%trg_hess_loaded = .false.
    geo%trg_esp_loaded = .false.
    geo%esp_npoints = 0

end subroutine ffdev_geometry_init

! ==============================================================================
! subroutine ffdev_geometry_load_xyz
! ==============================================================================

subroutine ffdev_geometry_load_xyz(geo,name)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    type(GEOMETRY)      :: geo
    character(*)        :: name
    ! --------------------------------------------
    type(XYZFILE_TYPE)  :: fin
    integer             :: alloc_stat, i
    ! --------------------------------------------------------------------------

    ! load xyz file
    call init_xyz(fin)
    call open_xyz(DEV_XYZ,name,fin,'OLD')
    call read_xyz(DEV_XYZ,fin)

    ! allocate data
    geo%natoms = fin%natoms
    geo%name   = name
    geo%title  = fin%comment

    allocate( geo%crd(3,geo%natoms), geo%z(geo%natoms), stat = alloc_stat )
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for geometry!')
    end if

    geo%crd(:,:) = fin%cvs(:,:)
    do i=1,fin%natoms
        ! is symbol?
        geo%z(i) = SearchZBySymbol(fin%symbols(i))
        if( geo%z(i) .eq. 0 ) then
            read(fin%symbols(i),*) geo%z(i)
        end if
    end do

    call close_xyz(DEV_XYZ,fin)
    call free_xyz(fin)

end subroutine ffdev_geometry_load_xyz

! ==============================================================================
! subroutine ffdev_geometry_save_xyz
! ==============================================================================

subroutine ffdev_geometry_save_xyz(geo,name)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use ffdev_topology
    use smf_periodic_table_dat

    implicit none
    type(GEOMETRY)      :: geo
    character(*)        :: name
    ! --------------------------------------------
    type(XYZFILE_TYPE)  :: fout
    integer             :: i
    ! --------------------------------------------------------------------------

    ! load xyz file
    call init_xyz(fout)
    call allocate_xyz(fout,geo%natoms)

    ! copy data
    fout%cvs(:,:) = geo%crd(:,:)
    fout%comment  = 'ffdevel xyz file'

    do i=1,geo%natoms
        fout%symbols(i) = pt_symbols(geo%z(i))
    end do

    ! write data
    call open_xyz(DEV_XYZ,name,fout,'UNKNOWN')
    call write_xyz(DEV_XYZ,fout)
    call close_xyz(DEV_XYZ,fout)

    call free_xyz(fout)

end subroutine ffdev_geometry_save_xyz

! ==============================================================================
! subroutine ffdev_geometry_load_point
! ==============================================================================

subroutine ffdev_geometry_load_point(geo,name)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    type(GEOMETRY)      :: geo
    character(*)        :: name
    ! --------------------------------------------
    integer             :: i,j,k,l,alloc_stat,read_stat
    character(len=255)  :: line,buffer
    character(len=80)   :: key,sym
    ! --------------------------------------------------------------------------

    geo%name = name

    ! open file
    call ffdev_utils_open(DEV_GEO,name,'O')

    ! load number of atoms - mandatory
    read(DEV_GEO,*,iostat = read_stat) geo%natoms
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to load number of atoms from geometry point!')
    end if

    ! load title - mandatory
    read(DEV_GEO,'(A80)',iostat = read_stat) geo%title
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to load title from geometry point!')
    end if

    ! load geometry - mandatory
    allocate( geo%z(geo%natoms), geo%crd(3,geo%natoms), stat = alloc_stat )
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for geometry!')
    end if
    do i=1,geo%natoms
        read(DEV_GEO,'(A255)',iostat = read_stat) line
        if( read_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to load line with coordinates!')
        end if
        read(line,*,iostat = read_stat) sym, geo%crd(1,i), geo%crd(2,i), geo%crd(3,i)
        if( read_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to load coordinates!')
        end if
        geo%z(i) = SearchZBySymbol(sym)
        if( geo%z(i) .eq. 0 ) then
            read(sym,*) geo%z(i)
        end if
    end do

    ! extra data - optional
    do while( .true. )
        ! read key line
        read(DEV_GEO,*,iostat = read_stat) key

        if( read_stat .lt. 0 ) exit ! end of file
        if( read_stat .gt. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to read data key token!')
        end if

        ! parse keys
        select case(trim(key))
            case('WEIGHT')
                read(DEV_GEO,*,iostat = read_stat) geo%weight
                if( read_stat .ne. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to read point weight!')
                end if
            case('ENERGY')
                read(DEV_GEO,*,iostat = read_stat) geo%trg_energy
                if( read_stat .ne. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to read point energy!')
                end if
                geo%trg_ene_loaded = .true.
            case('GRADIENT')
                allocate( geo%trg_grd(3,geo%natoms), stat = alloc_stat )
                if( alloc_stat .ne. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for trg_grd!')
                end if
                do i=1,geo%natoms
                    read(DEV_GEO,*,iostat = read_stat) geo%trg_grd(1,i), geo%trg_grd(2,i), geo%trg_grd(3,i)
                    if( read_stat .ne. 0 ) then
                        write(buffer,'(A,I3)') 'Unable to read GRADIENT entry! Gradient line = ',i
                        call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                    end if
                end do
                geo%trg_grd_loaded = .true.
            case('HESSIAN')
                allocate( geo%trg_hess(3,geo%natoms,3,geo%natoms), stat = alloc_stat )
                if( alloc_stat .ne. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for trg_hess!')
                end if
                do i=1,geo%natoms
                    do j=1,3
                        read(DEV_GEO,*,iostat = read_stat) ((geo%trg_hess(j,i,k,l), k=1,3), l=1,geo%natoms)
                        if( read_stat .ne. 0 ) then
                            write(buffer,'(A,I3)') 'Unable to read HESSIAN entry! Hessian line = ',(i-1)*3 + j
                            call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                        end if
                    end do
                end do
                geo%trg_hess_loaded = .true.
            case('ESP')
                ! read number of ESP points
                read(DEV_GEO,*,iostat = read_stat) geo%esp_npoints
                if( read_stat .ne. 0 ) then
                    write(buffer,'(A)') 'Unable to read ESP entry - number of ESP points!'
                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                end if
                if( geo%esp_npoints .le. 0 ) then
                    write(buffer,'(A)') 'Nunber of ESP points must be larger than zero!'
                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                end if
                allocate( geo%trg_esp(4,geo%esp_npoints), stat = alloc_stat )
                if( alloc_stat .ne. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for trg_esp!')
                end if
                do i=1,geo%esp_npoints
                    read(DEV_GEO,*,iostat = read_stat) geo%trg_esp(1,i), geo%trg_esp(2,i), &
                                                       geo%trg_esp(3,i), geo%trg_esp(4,i)
                    if( read_stat .ne. 0 ) then
                        write(buffer,'(A,I3)') 'Unable to read ESP entry! ESP line = ',i
                        call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                    end if
                end do
                geo%trg_esp_loaded = .true.
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unrecognized data point key ''' // trim(key) // '''!')
        end select
    end do

    ! close file
    close(DEV_GEO)

end subroutine ffdev_geometry_load_point

! ==============================================================================
! subroutine ffdev_geometry_save_point
! ==============================================================================

subroutine ffdev_geometry_save_point(geo,name)

    use ffdev_utils
    use smf_periodic_table_dat

    implicit none
    type(GEOMETRY)      :: geo
    character(*)        :: name
    ! --------------------------------------------
    integer             :: i,j,k,l,write_stat
    ! --------------------------------------------------------------------------

    geo%name = name

    ! open file
    call ffdev_utils_open(DEV_GEO,name,'U')

    ! number of atoms
    write(DEV_GEO,'(I8)',iostat = write_stat) geo%natoms
    if( write_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to save number of atoms!')
    end if

    ! title
    write(DEV_GEO,'(A80)',iostat = write_stat) geo%title
    if( write_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to save title!')
    end if

    ! save geometry - mandatory
    do i=1,geo%natoms
        write(DEV_GEO,'(A2,1X,F12.6,1X,F12.6,1X,F12.6)',iostat = write_stat) &
                       pt_symbols(geo%z(i)), geo%crd(1,i), geo%crd(2,i), geo%crd(3,i)
        if( write_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to save line with coordinates!')
        end if
    end do

    ! extra data - optional
    write(DEV_GEO,'(A)',iostat = write_stat) 'WEIGHT'
    write(DEV_GEO,'(F12.6)',iostat = write_stat) geo%weight
    if( write_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to save line with weight!')
    end if

    write(DEV_GEO,'(A)',iostat = write_stat) 'ENERGY'
    write(DEV_GEO,'(F20.6)',iostat = write_stat) geo%total_ene
    if( write_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to save line with energy!')
    end if

    write(DEV_GEO,'(A)',iostat = write_stat) 'GRADIENT'
    do i=1,geo%natoms
        write(DEV_GEO,'(F12.6,1X,F12.6,1X,F12.6)',iostat = write_stat) &
                       geo%grd(1,i), geo%grd(2,i), geo%grd(3,i)
        if( write_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to save line with gradient!')
        end if
    end do

    write(DEV_GEO,'(A)',iostat = write_stat) 'HESSIAN'
    do i=1,geo%natoms
        do j=1,3
            do k=1,geo%natoms
                do l=1,3
                    write(DEV_GEO,'(F12.6,1X)',iostat = write_stat,advance='NO') geo%hess(j,i,l,k)
                    if( write_stat .ne. 0 ) then
                        call ffdev_utils_exit(DEV_OUT,1,'Unable to save line with hessian!')
                    end if
                end do
            end do
            write(DEV_GEO,*)
        end do
    end do

    ! close file
    close(DEV_GEO)

end subroutine ffdev_geometry_save_point

! ==============================================================================
! subroutine ffdev_geometry_load_ginp
! ==============================================================================

subroutine ffdev_geometry_load_ginp(geo,name,mode)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    type(GEOMETRY)      :: geo
    character(*)        :: name
    integer             :: mode
    ! --------------------------------------------
    integer             :: alloc_stat, i, read_stat
    ! --------------------------------------------------------------------------

    mode = 0

    ! open file
    call ffdev_utils_open(DEV_GEO,name,'O')

    ! load number of atoms - mandatory
    read(DEV_GEO,*,iostat = read_stat) geo%natoms, mode
    if( read_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to load number of atoms from input file!')
    end if

    geo%name   = name
    geo%title  = 'gaussian input file'

    allocate( geo%crd(3,geo%natoms), geo%z(geo%natoms), stat = alloc_stat )
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for geometry!')
    end if

    ! read atoms
    do i=1,geo%natoms
        read(DEV_GEO,*,iostat = read_stat) geo%z(i), geo%crd(1,i), geo%crd(2,i), geo%crd(3,i)
        if( read_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to load atom coordinates!')
        end if
        ! convert from au to A
        geo%crd(:,i) = geo%crd(:,i)*DEV_AU2A
    end do

    close(DEV_GEO)


end subroutine ffdev_geometry_load_ginp

! ==============================================================================
! subroutine ffdev_geometry_save_gout
! ==============================================================================

subroutine ffdev_geometry_save_gout(geo,name,mode)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    type(GEOMETRY)      :: geo
    character(*)        :: name
    integer             :: mode
    ! --------------------------------------------
    integer             :: i,j,k,l,count,ni,nk
    double precision    :: fac
    ! --------------------------------------------------------------------------

    ! open file
    call ffdev_utils_open(DEV_GEO,name,'U')

    ! energy, dipole-moment (xyz)	    	E, Dip(I), I=1,3
    fac = DEV_KCL2HARTREE
    write(DEV_GEO,'(4D20.12)') geo%total_ene*fac, 0.0, 0.0, 0.0

    if( (mode .eq. 1) .or. (mode .eq. 2) ) then
    ! gradient on atom (xyz)	    	FX(J,I), J=1,3; I=1,NAtoms
        fac = DEV_KCL2HARTREE / DEV_A2AU
        do i = 1, geo%natoms
            write(DEV_GEO,'(3D20.12)') geo%grd(1,i)*fac, geo%grd(2,i)*fac,geo%grd(3,i)*fac
        end do
    end if

    ! polarizability	    	Polar(I), I=1,6	    	3D20.12
    write(DEV_GEO,'(3D20.12)') 0.0, 0.0, 0.0
    write(DEV_GEO,'(3D20.12)') 0.0, 0.0, 0.0

    ! dipole derivatives	    	DDip(I), I=1,9*NAtoms	    	3D20.12
    do i = 1, 3*geo%natoms
        write(DEV_GEO,'(3D20.12)') 0.0, 0.0, 0.0
    end do

    if( mode .eq. 2 ) then
        ! force constants	    	FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2	    	3D20.12
        fac = DEV_KCL2HARTREE / DEV_A2AU**2
        count = 0
        do ni=1,3*geo%natoms
            do nk=1,ni
                i = (ni-1) / 3 + 1
                j = mod(ni-1,3) + 1
                k = (nk-1) / 3 + 1
                l = mod(nk-1,3) + 1
                write(DEV_GEO,'(D20.12)',ADVANCE='NO') geo%hess(j,i,l,k)*fac
                count = count + 1
                if( mod(count,3) .eq. 0 ) write(DEV_GEO,*)
            end do
        end do
        write(DEV_GEO,*)
    end if

    close(DEV_GEO)


end subroutine ffdev_geometry_save_gout

! ==============================================================================
! subroutine ffdev_geometry_info_input
! ==============================================================================

subroutine ffdev_geometry_info_input(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,90)  trim(geo%name)
    write(DEV_OUT,95)  trim(geo%title)
    write(DEV_OUT,100) geo%natoms

 90 format('Geometry name          = ',A)
 95 format('Geometry comment       = ',A)
100 format('Number of atoms        = ',I6)

end subroutine ffdev_geometry_info_input

! ==============================================================================
! subroutine ffdev_geometry_info_point_header
! ==============================================================================

subroutine ffdev_geometry_info_point_header()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,30)
    write(DEV_OUT,40)

30 format('# ID   File                 Weight E G H P')
40 format('# ---- -------------------- ------ - - - -')

end subroutine ffdev_geometry_info_point_header

! ==============================================================================
! subroutine ffdev_geometry_info_point
! ==============================================================================

subroutine ffdev_geometry_info_point(geo)

    implicit none
    type(GEOMETRY)      :: geo
    ! --------------------------------------------
    character(len=20)   :: lname
    ! --------------------------------------------------------------------------

    lname = trim(geo%name)
    write(DEV_OUT,10) geo%id,adjustl(lname),geo%weight,geo%trg_ene_loaded, &
                      geo%trg_grd_loaded, geo%trg_hess_loaded, geo%trg_esp_loaded

! '# ---- -------------------- ------ - - - -'
  10 format(I6,1X,A20,1X,F6.3,1X,L1,1X,L1,1X,L1,1X,L1)

end subroutine ffdev_geometry_info_point

! ==============================================================================
! subroutine ffdev_geometry_info_point_header_ext
! ==============================================================================

subroutine ffdev_geometry_info_point_header_ext()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,30)
    write(DEV_OUT,40)

30 format('# ID   File                 Weight    Rel Energy    E G H P')
40 format('# ---- -------------------- ------ ---------------- - - - -')

end subroutine ffdev_geometry_info_point_header_ext

! ==============================================================================
! subroutine ffdev_geometry_info_point_ext
! ==============================================================================

subroutine ffdev_geometry_info_point_ext(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    character(len=20)   :: lname
    ! --------------------------------------------------------------------------

    lname = trim(geo%name)
    write(DEV_OUT,10) geo%id,adjustl(lname),geo%weight,geo%trg_energy,geo%trg_ene_loaded, &
                      geo%trg_grd_loaded, geo%trg_hess_loaded, geo%trg_esp_loaded

! '# ---- -------------------- ------ - - - -'
  10 format(I6,1X,A20,1X,F6.3,1X,F16.3,1X,L1,1X,L1,1X,L1,1X,L1)

end subroutine ffdev_geometry_info_point_ext

! ==============================================================================
! subroutine ffdev_geometry_info_ene
! ==============================================================================

subroutine ffdev_geometry_info_ene(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    write(DEV_OUT,100) geo%bond_ene
    write(DEV_OUT,110) geo%angle_ene
    write(DEV_OUT,120) geo%dih_ene
    write(DEV_OUT,130) geo%impropr_ene
    write(DEV_OUT,140) geo%ele14_ene
    write(DEV_OUT,150) geo%nb14_ene
    write(DEV_OUT,160) geo%ele_ene
    write(DEV_OUT,170) geo%nb_ene
    write(DEV_OUT,180) geo%total_ene

100 format('Ebonds     = ',F20.7)
110 format('Eangles    = ',F20.7)
120 format('Edihedrals = ',F20.7)
130 format('Eimpropers = ',F20.7)
140 format('E14ele     = ',F20.7)
150 format('E14vdw     = ',F20.7)
160 format('Eele       = ',F20.7)
170 format('Evdw       = ',F20.7)
180 format('Etotal     = ',F20.7)

end subroutine ffdev_geometry_info_ene

! ==============================================================================
! subroutine ffdev_geometry_info_point_header
! ==============================================================================

subroutine ffdev_geometry_print_xyz(geo)

    use smf_periodic_table_dat

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    write(DEV_OUT,5) geo%natoms
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    do i=1,geo%natoms
        write(DEV_OUT,30) i,geo%z(i),pt_symbols(geo%z(i)), &
                          geo%crd(1,i), geo%crd(2,i), geo%crd(3,i)
    end do

 5 format('# FFDEVEL XYZ ',I8)
10 format('# ID    Z Sym         X                Y               Z        ')
20 format('# ---- -- --- ---------------- ---------------- ----------------')
30 format(I6,1X,I2,1X,A3,1X,F16.8,1X,F16.8,1X,F16.8)

end subroutine ffdev_geometry_print_xyz

! ==============================================================================
! subroutine ffdev_geometry_copy
! ==============================================================================

subroutine ffdev_geometry_copy(geo1,geo2)

    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo1
    type(GEOMETRY)  :: geo2
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    geo1%natoms = geo2%natoms
    geo1%name = geo2%name
    geo1%title = geo2%title

    if( geo1%natoms .gt. 0 ) then
        allocate(geo1%crd(3,geo1%natoms), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate crd array for geometry!')
        end if
        geo1%crd(:,:) = geo2%crd(:,:)
    end if

end subroutine ffdev_geometry_copy

! ==============================================================================
! subroutine ffdev_geometry_allocate
! ==============================================================================

subroutine ffdev_geometry_allocate(geo,natoms)

    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: natoms
    ! --------------------------------------------
    integer         :: alloc_stat
    ! --------------------------------------------------------------------------

    geo%natoms = natoms
    geo%name = ''
    geo%title = ''

    if( geo%natoms .gt. 0 ) then
        allocate(geo%crd(3,geo%natoms), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate crd array for geometry!')
        end if
    end if

end subroutine ffdev_geometry_allocate

! ==============================================================================
! subroutine ffdev_geometry_check_z
! ==============================================================================

subroutine ffdev_geometry_check_z(top,geo)

    use ffdev_utils
    use ffdev_topology

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer             :: i
    character(len=128)  :: msg
    ! --------------------------------------------------------------------------

    if( top%natoms .ne. geo%natoms ) then
        call ffdev_utils_exit(DEV_OUT,1,'Topology and geometry does not have the same number of atoms!')
    end if

    do i=1,top%natoms
        if( top%atom_types(top%atoms(i)%typeid)%z .ne. geo%z(i) ) then
            write(msg,'(A,I3,A,I2,A,I2)') 'Inconsistency at position: ',i,'! Topology atom Z: ', &
                             top%atom_types(top%atoms(i)%typeid)%z,' Geometry atom Z: ',geo%z(i)
            call ffdev_utils_exit(DEV_OUT,1,trim(msg))
        end if
    end do

end subroutine ffdev_geometry_check_z

! ==============================================================================
! subroutine ffdev_geometry_get_length
! ==============================================================================

real(DEVDP) function ffdev_geometry_get_length(geo,ai,aj)

    use ffdev_utils
    use ffdev_topology

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: ai
    integer         :: aj
    ! --------------------------------------------
    real(DEVDP)     :: v
    ! --------------------------------------------------------------------------

    v = 0.0
    v = v + (geo%crd(1,ai) - geo%crd(1,aj))**2
    v = v + (geo%crd(2,ai) - geo%crd(2,aj))**2
    v = v + (geo%crd(3,ai) - geo%crd(3,aj))**2
    ffdev_geometry_get_length = sqrt(v)

end function ffdev_geometry_get_length

! ==============================================================================
! subroutine ffdev_geometry_get_angle
! ==============================================================================

real(DEVDP) function ffdev_geometry_get_angle(geo,ai,aj,ak)

    use ffdev_utils
    use ffdev_topology

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: ai
    integer         :: aj
    integer         :: ak
    ! --------------------------------------------
    real(DEVDP)     :: bjiinv, bjkinv, bji2inv, bjk2inv
    real(DEVDP)     :: scp,angv
    real(DEVDP)     :: rji(3),rjk(3)
    ! --------------------------------------------------------------------------

    ! calculate rji and rjk
    rji(:) = geo%crd(:,ai) - geo%crd(:,aj)
    rjk(:) = geo%crd(:,ak) - geo%crd(:,aj)

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
    ffdev_geometry_get_angle = angv

end function ffdev_geometry_get_angle

! ==============================================================================
! subroutine ffdev_geometry_get_dihedral
! ==============================================================================

real(DEVDP) function ffdev_geometry_get_dihedral(geo,ai,aj,ak,al)

    use ffdev_utils
    use ffdev_topology

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: ai
    integer         :: aj
    integer         :: ak
    integer         :: al
    ! --------------------------------------------
    real(DEVDP)     ::  scp,phi
    real(DEVDP)     ::  bjinv, bkinv, bj2inv, bk2inv
    real(DEVDP)     ::  rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
    ! --------------------------------------------------------------------------

    ! calculate rji and rjk
    rji(:) = geo%crd(:,ai) - geo%crd(:,aj)
    rjk(:) = geo%crd(:,ak) - geo%crd(:,aj)
    rkl(:) = geo%crd(:,al) - geo%crd(:,ak)
    rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
    rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
    rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
    rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
    rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
    rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

    bj2inv = 1.0d0/(rnj(1)**2 + rnj(2)**2 + rnj(3)**2 )
    bk2inv = 1.0d0/(rnk(1)**2 + rnk(2)**2 + rnk(3)**2 )
    bjinv = sqrt(bj2inv)
    bkinv = sqrt(bk2inv)

    ! calculate scp and phi
    scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
    if ( scp .gt.  1.0 ) then
            scp =  1.0
            phi = acos (1.0) ! const
    else if ( scp .lt. -1.0 ) then
            scp = -1.0
            phi = acos (-1.0) ! const
    else
        phi = acos ( scp )
    end if
    if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
       +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
       +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) .lt. 0) then
                phi = -phi
    end if
    ffdev_geometry_get_dihedral = phi

end function ffdev_geometry_get_dihedral

! ------------------------------------------------------------------------------

end module ffdev_geometry

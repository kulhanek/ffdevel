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
    geo%trg_crd_loaded = .false.
    geo%trg_grd_loaded = .false.
    geo%trg_hess_loaded = .false.
    geo%trg_freq_loaded = .false.
    geo%trg_esp_loaded = .false.
    geo%trg_crd_optimized = .false.
    geo%esp_npoints = 0

    geo%nrst        = 0
    geo%rst_energy  = 0.0

    geo%z           => null()
    geo%crd         => null()

    geo%grd         => null()
    geo%hess        => null()
    geo%nmodes      => null()
    geo%freq        => null()

    geo%trg_crd         => null()
    geo%trg_grd         => null()
    geo%trg_hess        => null()
    geo%trg_nmodes      => null()
    geo%trg_freq        => null()
    geo%freq_t2s_map    => null()
    geo%freq_t2s_angles => null()
    geo%freq_t2s_rmsd   => null()
    geo%trg_esp         => null()

    geo%rst         => null()

end subroutine ffdev_geometry_init

! ==============================================================================
! subroutine ffdev_geometry_init
! ==============================================================================

subroutine ffdev_geometry_destroy(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    if( associated(geo%z) ) then
        deallocate(geo%z)
    end if
    if( associated(geo%crd) ) then
        deallocate(geo%crd)
    end if

    if( associated(geo%grd) ) then
        deallocate(geo%grd)
    end if
    if( associated(geo%hess) ) then
        deallocate(geo%hess)
    end if
    if( associated(geo%freq) ) then
        deallocate(geo%freq)
    end if

    if( associated(geo%trg_crd) ) then
        deallocate(geo%trg_crd)
    end if
    if( associated(geo%trg_grd) ) then
        deallocate(geo%trg_grd)
    end if
    if( associated(geo%trg_hess) ) then
        deallocate(geo%trg_hess)
    end if
    if( associated(geo%trg_freq) ) then
        deallocate(geo%trg_freq)
    end if

    if( associated(geo%freq_t2s_map) ) then
        deallocate(geo%freq_t2s_map)
    end if
    if( associated(geo%freq_t2s_angles) ) then
        deallocate(geo%freq_t2s_angles)
    end if
    if( associated(geo%freq_t2s_rmsd) ) then
        deallocate(geo%freq_t2s_rmsd)
    end if

    if( associated(geo%trg_esp) ) then
        deallocate(geo%trg_esp)
    end if

    if( associated(geo%rst) ) then
        do i=1,geo%nrst
            deallocate(geo%rst(i)%ai)
        end do
        deallocate(geo%rst)
    end if

end subroutine ffdev_geometry_destroy

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
! subroutine ffdev_geometry_load_xyz_snapshot
! ==============================================================================

subroutine ffdev_geometry_load_xyz_snapshot(geo,fin)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    type(GEOMETRY)      :: geo
    type(XYZFILE_TYPE)  :: fin
    ! --------------------------------------------
    integer             :: alloc_stat, i
    ! --------------------------------------------------------------------------

    ! allocate data
    geo%natoms = fin%natoms
    geo%name   = 'snapshot'
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

end subroutine ffdev_geometry_load_xyz_snapshot

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
! subroutine ffdev_geometry_save_xyz_snapshot
! ==============================================================================

subroutine ffdev_geometry_save_xyz_snapshot(geo,fout)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use ffdev_topology
    use smf_periodic_table_dat

    implicit none
    type(GEOMETRY)      :: geo
    type(XYZFILE_TYPE)  :: fout
    ! --------------------------------------------
    integer             :: i
    ! --------------------------------------------------------------------------

    ! copy data
    fout%cvs(:,:) = geo%crd(:,:)
    write(fout%comment,150) geo%total_ene

    do i=1,geo%natoms
        fout%symbols(i) = pt_symbols(geo%z(i))
    end do

150 format('E=',F20.10)

end subroutine ffdev_geometry_save_xyz_snapshot

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
    ! --------------------------------------------------------------------------

    geo%name = name

    ! open file
    call ffdev_utils_open(DEV_GEO,name,'O')

    ! read 1 point
    call ffdev_geometry_load_1point(geo,.false.)

    ! close file
    close(DEV_GEO)

end subroutine ffdev_geometry_load_point

! ==============================================================================
! subroutine ffdev_geometry_num_of_points
! ==============================================================================

integer function ffdev_geometry_num_of_points(name)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    character(*)        :: name
    ! --------------------------------------------
    type(GEOMETRY)      :: geo
    ! --------------------------------------------------------------------------

    ! open file
    call ffdev_utils_open(DEV_GEO,name,'O')

    ffdev_geometry_num_of_points = 0
    do while(.true.)
        ! init empty geo
        call ffdev_geometry_init(geo)
        ! read 1 point
        call ffdev_geometry_load_1point(geo,.true.)
        ! was it read?
        if( geo%natoms .eq. 0 ) then
            call ffdev_geometry_destroy(geo)
            exit   ! no -> exit
        end if
        ffdev_geometry_num_of_points = ffdev_geometry_num_of_points + 1
        ! release data
        call ffdev_geometry_destroy(geo)
    end do

    ! close file
    close(DEV_GEO)

end function ffdev_geometry_num_of_points

! ==============================================================================
! subroutine ffdev_geometry_load_1point
! ==============================================================================

subroutine ffdev_geometry_load_1point(geo,stream)

    use smf_xyzfile
    use smf_xyzfile_type
    use ffdev_utils
    use smf_periodic_table

    implicit none
    type(GEOMETRY)      :: geo
    logical             :: stream
    ! --------------------------------------------
    integer             :: i,j,k,l,alloc_stat,read_stat
    character(len=255)  :: line,buffer
    character(len=80)   :: key,sym
    ! --------------------------------------------------------------------------

    ! load number of atoms - mandatory
    read(DEV_GEO,*,iostat = read_stat) geo%natoms
    if( read_stat .ne. 0 ) then
        if( stream ) then
            geo%natoms = 0
            return
        else
            call ffdev_utils_exit(DEV_OUT,1,'Unable to load number of atoms from geometry point!')
        end if
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
    if( geo%trg_crd_loaded ) then
        ! always load trg_geo
        allocate( geo%trg_crd(3,geo%natoms), stat = alloc_stat )
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate arays for geometry!')
        end if
        geo%trg_crd = geo%crd
    end if
    
    ! extra data - optional
    do while( .true. )
        ! read key line
        read(DEV_GEO,*,iostat = read_stat) key

        if( read_stat .lt. 0 ) exit ! end of file
        if( read_stat .gt. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to read data key token!')
        end if

        ! is it a number?
        read(key,*,iostat = read_stat) k
        if( read_stat .eq. 0 ) then
            ! backspace record and return
            backspace(DEV_GEO)
            exit
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
            case('RST')
                read(DEV_GEO,*,iostat = read_stat) geo%nrst
                if( read_stat .ne. 0 ) then
                    call ffdev_utils_exit(DEV_OUT,1,'Unable to read nrst entry!')
                end if
                if( geo%nrst .gt. 0 ) then
                    allocate( geo%rst(geo%nrst), stat = alloc_stat )
                    if( alloc_stat .ne. 0 ) then
                        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate aray for rst!')
                    end if
                    do i=1,geo%nrst
                        read(DEV_GEO,'(A255)',iostat = read_stat) line
                        if( read_stat .ne. 0 ) then
                            write(buffer,'(A,I3)') 'Unable to read RST entry! RST line = ',i
                            call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                        end if
                        geo%rst(i)%cvtype = ''
                        j = 0
                        read(line,*,iostat = read_stat) j,geo%rst(i)%cvtype
                        if( i .ne. j ) then
                            write(buffer,'(A,I3)') 'Unable to read RST entry - order mismatch! RST line = ',i
                            call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                        end if
                        select case(trim(geo%rst(i)%cvtype))
                            case('B')
                                allocate( geo%rst(i)%ai(2), stat = alloc_stat )
                                if( alloc_stat .ne. 0 ) then
                                    write(buffer,'(A,I3)') 'Unable to allocate ai for RST entry B! RST line = ',i
                                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                                end if
                                read(line,*,iostat = read_stat) j,geo%rst(i)%cvtype,geo%rst(i)%trg_value, geo%rst(i)%ai(1), &
                                                                   geo%rst(i)%ai(2)
                                if( read_stat .ne. 0 ) then
                                    write(buffer,'(A,I3)') 'Unable to read RST entry B! RST line = ',i
                                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                                end if
                            case('A')
                                allocate( geo%rst(i)%ai(3), stat = alloc_stat )
                                if( alloc_stat .ne. 0 ) then
                                    write(buffer,'(A,I3)') 'Unable to allocate ai for RST entry A! RST line = ',i
                                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                                end if
                                read(line,*,iostat = read_stat) j,geo%rst(i)%cvtype,geo%rst(i)%trg_value,geo%rst(i)%ai(1), &
                                                                   geo%rst(i)%ai(2),geo%rst(i)%ai(3)
                                if( read_stat .ne. 0 ) then
                                    write(buffer,'(A,I3)') 'Unable to read RST entry A! RST line = ',i
                                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                                end if
                                ! convert to rad
                                geo%rst(i)%trg_value = geo%rst(i)%trg_value * DEV_D2R
                            case('D')
                                allocate( geo%rst(i)%ai(4), stat = alloc_stat )
                                if( alloc_stat .ne. 0 ) then
                                    write(buffer,'(A,I3)') 'Unable to allocate ai for RST entry D! RST line = ',i
                                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                                end if
                                read(line,*,iostat = read_stat) j,geo%rst(i)%cvtype,geo%rst(i)%trg_value,geo%rst(i)%ai(1), &
                                                                   geo%rst(i)%ai(2),geo%rst(i)%ai(3),geo%rst(i)%ai(4)
                                if( read_stat .ne. 0 ) then
                                    write(buffer,'(A,I3)') 'Unable to read RST entry D! RST line = ',i
                                    call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                                end if
                                ! convert to rad
                                geo%rst(i)%trg_value = geo%rst(i)%trg_value * DEV_D2R
                            case default
                                write(buffer,'(A,A,A,I3)') 'Unsupported RST type ',trim(geo%rst(i)%cvtype),'! RST line = ',i
                                call ffdev_utils_exit(DEV_OUT,1,trim(buffer))
                        end select
                    end do
                end if
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unrecognized data point key ''' // trim(key) // '''!')
        end select
    end do

end subroutine ffdev_geometry_load_1point

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

    if( associated(geo%grd) ) then
        write(DEV_GEO,'(A)',iostat = write_stat) 'GRADIENT'
        do i=1,geo%natoms
            write(DEV_GEO,'(F12.6,1X,F12.6,1X,F12.6)',iostat = write_stat) &
                           geo%grd(1,i), geo%grd(2,i), geo%grd(3,i)
            if( write_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_OUT,1,'Unable to save line with gradient!')
            end if
        end do
    end if

    if( associated(geo%hess) ) then
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
    end if

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
    write(DEV_OUT,110) geo%nrst

 90 format('Geometry name          = ',A)
 95 format('Geometry comment       = ',A)
100 format('Number of atoms        = ',I6)
110 format('Number of CVs          = ',I6)

end subroutine ffdev_geometry_info_input

! ==============================================================================
! subroutine ffdev_geometry_info_point_header
! ==============================================================================

subroutine ffdev_geometry_info_point_header()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,30)
    write(DEV_OUT,40)

30 format('# ID   File                                     Weight E G H P')
40 format('# ---- ---------------------------------------- ------ - - - -')

end subroutine ffdev_geometry_info_point_header

! ==============================================================================
! subroutine ffdev_geometry_info_point
! ==============================================================================

subroutine ffdev_geometry_info_point(geo)

    implicit none
    type(GEOMETRY)      :: geo
    ! --------------------------------------------
    character(len=40)   :: lname
    ! --------------------------------------------------------------------------

    lname = trim(geo%name)
    write(DEV_OUT,10) geo%id,adjustl(lname),geo%weight,geo%trg_ene_loaded, &
                      geo%trg_grd_loaded, geo%trg_hess_loaded, geo%trg_esp_loaded

! '# ---- -------------------- ------ - - - -'
  10 format(I6,1X,A40,1X,F6.3,1X,L1,1X,L1,1X,L1,1X,L1)

end subroutine ffdev_geometry_info_point

! ==============================================================================
! subroutine ffdev_geometry_info_point_header_ext
! ==============================================================================

subroutine ffdev_geometry_info_point_header_ext()

    implicit none
    ! --------------------------------------------------------------------------

    write(DEV_OUT,30)
    write(DEV_OUT,40)

30 format('# ID   File                                     Weight    Rel Energy    E G H P')
40 format('# ---- ---------------------------------------- ------ ---------------- - - - -')

end subroutine ffdev_geometry_info_point_header_ext

! ==============================================================================
! subroutine ffdev_geometry_info_point_ext
! ==============================================================================

subroutine ffdev_geometry_info_point_ext(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    character(len=40)   :: lname
    ! --------------------------------------------------------------------------

    lname = trim(geo%name)
    write(DEV_OUT,10) geo%id,adjustl(lname),geo%weight,geo%trg_energy,geo%trg_ene_loaded, &
                      geo%trg_grd_loaded, geo%trg_hess_loaded, geo%trg_esp_loaded

! '# ---- -------------------- ------ - - - -'
  10 format(I6,1X,A40,1X,F6.3,1X,F16.4,1X,L1,1X,L1,1X,L1,1X,L1)

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
! geo1 <- geo2
! ==============================================================================

subroutine ffdev_geometry_copy(geo1,geo2)

    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo1
    type(GEOMETRY)  :: geo2
    ! --------------------------------------------
    integer         :: alloc_stat, i
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

    geo1%nrst = geo2%nrst
    if( geo1%natoms .gt. 0 ) then
        allocate(geo1%rst(geo1%nrst), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate rst array for geometry!')
        end if
        do i=1,geo1%nrst
            geo1%rst(i)%trg_value = geo2%rst(i)%trg_value
            geo1%rst(i)%cvtype = geo2%rst(i)%cvtype
            select case(trim(geo1%rst(i)%cvtype))
                case('B')
                    allocate(geo1%rst(i)%ai(2), stat = alloc_stat)
                    if( alloc_stat .ne. 0 ) then
                        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate ai array for geometry!')
                    end if
                case('A')
                    allocate(geo1%rst(i)%ai(3), stat = alloc_stat)
                    if( alloc_stat .ne. 0 ) then
                        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate ai array for geometry!')
                    end if
                case('D')
                    allocate(geo1%rst(i)%ai(4), stat = alloc_stat)
                    if( alloc_stat .ne. 0 ) then
                        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate ai array for geometry!')
                    end if
                case default
                    call ffdev_utils_exit(DEV_OUT,1,'Unsupported CV in ffdev_geometry_copy!')
            end select
            geo1%rst(i)%ai(:) = geo2%rst(i)%ai(:)
        end do
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

real(DEVDP) function ffdev_geometry_get_length(crd,ai,aj)

    use ffdev_utils
    use ffdev_topology

    implicit none
    real(DEVDP)     :: crd(:,:)
    integer         :: ai
    integer         :: aj
    ! --------------------------------------------
    real(DEVDP)     :: v
    ! --------------------------------------------------------------------------

    v = 0.0
    v = v + (crd(1,ai) - crd(1,aj))**2
    v = v + (crd(2,ai) - crd(2,aj))**2
    v = v + (crd(3,ai) - crd(3,aj))**2
    ffdev_geometry_get_length = sqrt(v)

end function ffdev_geometry_get_length

! ==============================================================================
! subroutine ffdev_geometry_get_angle
! ==============================================================================

real(DEVDP) function ffdev_geometry_get_angle(crd,ai,aj,ak)

    use ffdev_utils
    use ffdev_topology

    implicit none
    real(DEVDP)     :: crd(:,:)
    integer         :: ai
    integer         :: aj
    integer         :: ak
    ! --------------------------------------------
    real(DEVDP)     :: bjiinv, bjkinv, bji2inv, bjk2inv
    real(DEVDP)     :: scp,angv
    real(DEVDP)     :: rji(3),rjk(3)
    ! --------------------------------------------------------------------------

    ! calculate rji and rjk
    rji(:) = crd(:,ai) - crd(:,aj)
    rjk(:) = crd(:,ak) - crd(:,aj)

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

real(DEVDP) function ffdev_geometry_get_dihedral(crd,ai,aj,ak,al)

    use ffdev_utils
    use ffdev_topology

    implicit none
    real(DEVDP)     :: crd(:,:)
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
    rji(:) = crd(:,ai) - crd(:,aj)
    rjk(:) = crd(:,ak) - crd(:,aj)
    rkl(:) = crd(:,al) - crd(:,ak)
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

! ==============================================================================
! subroutine ffdev_geometry_get_improper
! ==============================================================================

real(DEVDP) function ffdev_geometry_get_improper(crd,ai,aj,ak,al)

    use ffdev_utils
    use ffdev_topology

    implicit none
    real(DEVDP)     :: crd(:,:)
    integer         :: ai
    integer         :: aj
    integer         :: ak
    integer         :: al
    ! --------------------------------------------
    real(DEVDP)     ::  scp,phi
    real(DEVDP)     ::  bjinv, bkinv, bj2inv, bk2inv
    real(DEVDP)     ::  rji(3),rjk(3),rkl(3),rnj(3),rnk(3)
    ! --------------------------------------------------------------------------

    rji(:) = crd(:,ai) - crd(:,aj)
    rjk(:) = crd(:,ak) - crd(:,aj)
    rkl(:) = crd(:,al) - crd(:,ak)

    rnj(1) =  rji(2)*rjk(3) - rji(3)*rjk(2)
    rnj(2) =  rji(3)*rjk(1) - rji(1)*rjk(3)
    rnj(3) =  rji(1)*rjk(2) - rji(2)*rjk(1)
    rnk(1) = -rjk(2)*rkl(3) + rjk(3)*rkl(2)
    rnk(2) = -rjk(3)*rkl(1) + rjk(1)*rkl(3)
    rnk(3) = -rjk(1)*rkl(2) + rjk(2)*rkl(1)

    bj2inv  = 1.0d0/( rnj(1)**2 + rnj(2)**2 + rnj(3)**2)
    bk2inv  = 1.0d0/( rnk(1)**2 + rnk(2)**2 + rnk(3)**2)
    bjinv = sqrt(bj2inv)
    bkinv = sqrt(bk2inv)

    scp = (rnj(1)*rnk(1)+rnj(2)*rnk(2)+rnj(3)*rnk(3))*(bjinv*bkinv)
    if ( scp .gt.  1.0d0 ) then
        scp =  1.0d0
    else if ( scp .lt. -1.0d0 ) then
        scp = -1.0d0
    end if
    phi = acos ( scp )
    if(rjk(1)*(rnj(2)*rnk(3)-rnj(3)*rnk(2)) &
        +rjk(2)*(rnj(3)*rnk(1)-rnj(1)*rnk(3)) &
        +rjk(3)*(rnj(1)*rnk(2)-rnj(2)*rnk(1)) &
        .lt. 0) then
              phi = -phi
    end if
    ffdev_geometry_get_improper = phi

end function ffdev_geometry_get_improper

!===============================================================================
! Function:  ffdev_geometry_get_dihedral_deviation
! it consider periodicity of torsion
!===============================================================================

real(DEVDP) function ffdev_geometry_get_dihedral_deviation(value1,value2)

    implicit none
    real(DEVDP)     :: value1
    real(DEVDP)     :: value2
    ! --------------------------------------------
    real(DEVDP)     :: minv,maxv,vec
    ! --------------------------------------------------------------------------

    minv = -DEV_PI
    maxv =  DEV_PI

    if( abs(value1-value2) .lt. 0.5d0*(maxv-minv) ) then
        ffdev_geometry_get_dihedral_deviation = value1 - value2
        return
    else
        ! get vector
        vec = value1 - value2
        ! shift to box center
        vec = vec + 0.5d0*(maxv+minv)
        ! image as point
        vec = vec - (maxv-minv)*floor((vec-minv)/(maxv-minv))
        ! return vector back
        ffdev_geometry_get_dihedral_deviation = vec - 0.5d0*(maxv+minv)
        
        ! debug
        ! write(PMF_DEBUG,*) 'get_deviation:',value1,value2,get_deviation
    end if
        
end function ffdev_geometry_get_dihedral_deviation

! ==============================================================================
! subroutine ffdev_geometry_get_rst_penalty
! including gardient
! ==============================================================================

subroutine ffdev_geometry_get_rst_penalty(geo)

    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    ! reset penalty
    geo%rst_energy = 0.0d0

    do i=1,geo%nrst
        select case(trim(geo%rst(i)%cvtype))
            case('B')
                call ffdev_geometry_get_bond_penalty(geo,i)
            case('A')
                call ffdev_geometry_get_angle_penalty(geo,i)
            case('D')
                call ffdev_geometry_get_torsion_penalty(geo,i)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported CV in ffdev_geometry_get_rst_penalty!')
        end select
    end do

end subroutine ffdev_geometry_get_rst_penalty

! ==============================================================================
! subroutine ffdev_geometry_get_rst_penalty_only
! without gradient
! ==============================================================================

subroutine ffdev_geometry_get_rst_penalty_only(geo)

    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: actv,tv,dev
    ! --------------------------------------------------------------------------

    ! reset penalty
    geo%rst_energy = 0.0d0

    call ffdev_geometry_get_rst_values(geo)

    do i=1,geo%nrst
        tv = geo%rst(i)%trg_value
        actv = geo%rst(i)%value
        select case(trim(geo%rst(i)%cvtype))
            case('B')
                dev = actv - tv
                geo%rst_energy = geo%rst_energy + 0.5d0*DIS_FC*dev**2
            case('A')
                dev = actv - tv
                geo%rst_energy = geo%rst_energy + 0.5d0*ANG_FC*dev**2
            case('D')
                dev = ffdev_geometry_get_dihedral_deviation(actv,tv)
                geo%rst_energy = geo%rst_energy + 0.5d0*ANG_FC*dev**2
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported CV in ffdev_geometry_get_rst_penalty_only!')
        end select
    end do

end subroutine ffdev_geometry_get_rst_penalty_only

! ==============================================================================
! subroutine ffdev_geometry_get_rst_penalty_num
! ==============================================================================

subroutine ffdev_geometry_get_rst_penalty_num(geo)

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------
    type(GEOMETRY)  :: tmp_geo
    real(DEVDP)     :: d,ene1,ene2
    integer         :: i,j
    ! --------------------------------------------------------------------------

    d = 5.0d-4  ! differentiation parameter

    ! calculate base penalty
    call ffdev_geometry_get_rst_penalty_only(geo)

    ! init temporary geometry object
    call ffdev_geometry_init(tmp_geo)
    call ffdev_geometry_copy(tmp_geo,geo)

    ! gradient by numerical differentiation
    do i=1,geo%natoms
        do j=1,3
            ! left
            tmp_geo%crd(j,i) = geo%crd(j,i) + d
            call ffdev_geometry_get_rst_penalty_only(tmp_geo)
            ene1 = tmp_geo%rst_energy

            ! right
            tmp_geo%crd(j,i) = geo%crd(j,i) - d
            call ffdev_geometry_get_rst_penalty_only(tmp_geo)
            ene2 = tmp_geo%rst_energy

            ! gradient
            geo%grd(j,i) = 0.5d0*(ene1-ene2)/d

            ! move back
            tmp_geo%crd(j,i) = geo%crd(j,i)
        end do
    end do

    ! release temporary geometry object
    call ffdev_geometry_destroy(tmp_geo)

end subroutine ffdev_geometry_get_rst_penalty_num

! ==============================================================================
! function ffdev_geometry_get_rst_values
! ==============================================================================

subroutine ffdev_geometry_get_rst_values(geo)

    use ffdev_utils

    implicit none
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,geo%nrst
        select case(trim(geo%rst(i)%cvtype))
            case('B')
                geo%rst(i)%value = &
                      ffdev_geometry_get_length(geo%crd,geo%rst(i)%ai(1),geo%rst(i)%ai(2))
            case('A')
                geo%rst(i)%value = &
                      ffdev_geometry_get_angle(geo%crd,geo%rst(i)%ai(1),geo%rst(i)%ai(2),geo%rst(i)%ai(3))
            case('D')
                geo%rst(i)%value = &
                      ffdev_geometry_get_dihedral(geo%crd,geo%rst(i)%ai(1),geo%rst(i)%ai(2), &
                                                          geo%rst(i)%ai(3),geo%rst(i)%ai(4))
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported CV in ffdev_geometry_get_rst_values!')
        end select
    end do

end subroutine ffdev_geometry_get_rst_values

! ==============================================================================
! subroutine ffdev_geometry_rstsum
! ==============================================================================

subroutine ffdev_geometry_rstsum(unid,geo)

    use ffdev_utils

    implicit none
    integer         :: unid
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: actv,tv,dev
    ! --------------------------------------------------------------------------

    if( geo%nrst .le. 0 ) then
        return
    end if

    write(unid,10)
    write(unid,20)
    write(unid,30)

    call ffdev_geometry_get_rst_values(geo)

    do i=1,geo%nrst
        select case(trim(geo%rst(i)%cvtype))
            case('B')
                actv = geo%rst(i)%value
                tv = geo%rst(i)%trg_value
                dev = actv - tv
            case('A')
                actv = geo%rst(i)%value * DEV_R2D
                tv = geo%rst(i)%trg_value * DEV_R2D
                dev = actv - tv
            case('D')
                actv = geo%rst(i)%value * DEV_R2D
                tv = geo%rst(i)%trg_value * DEV_R2D
                dev = ffdev_geometry_get_dihedral_deviation(actv,tv)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported CV in ffdev_geometry_rstsum!')
        end select
        write(unid,40) i, geo%rst(i)%cvtype, tv, actv, dev
    end do

 10 format('# Geometry Restraint Summary')
 20 format('# ID  Type  Target value Curent value  Deviation  ')
 30 format('# --- ----- ------------ ------------ ------------')
 40 format(I5,1X,A5,1X,F12.6,1X,F12.6,1X,F12.6)

end subroutine ffdev_geometry_rstsum

! ==============================================================================
! subroutine ffdev_geometry_get_bond_penalty
! ==============================================================================

subroutine ffdev_geometry_get_bond_penalty(geo,cvidx)

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: cvidx
    ! --------------------------------------------
    integer         :: i,j
    real(DEVDP)     :: b,db,dv,rij(3)
    ! --------------------------------------------------------------------------

    ! for each bond
    i  = geo%rst(cvidx)%ai(1)
    j  = geo%rst(cvidx)%ai(2)

    ! calculate rij
    rij(:) = geo%crd(:,j) - geo%crd(:,i)

    ! calculate penalty
    b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
    db = b - geo%rst(cvidx)%trg_value
    geo%rst(cvidx)%value = 0.5*DIS_FC*db**2
    geo%rst_energy = geo%rst_energy + geo%rst(cvidx)%value

    ! calculate gradient
    dv = DIS_FC*db/b
    geo%grd(:,j) = geo%grd(:,j) + rij(:)*dv
    geo%grd(:,i) = geo%grd(:,i) - rij(:)*dv

end subroutine ffdev_geometry_get_bond_penalty

! ==============================================================================
! subroutine ffdev_geometry_get_angle_penalty
! ==============================================================================

subroutine ffdev_geometry_get_angle_penalty(geo,cvidx)

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: cvidx
    ! --------------------------------------------
    integer         :: i,j,k
    real(DEVDP)     :: bji2inv,bjk2inv,bjiinv,bjkinv,scp,angv,da,dv,f1
    real(DEVDP)     :: rji(3),rjk(3),di(3),dk(3)
    ! --------------------------------------------------------------------------

    i  = geo%rst(cvidx)%ai(1)
    j  = geo%rst(cvidx)%ai(2)
    k  = geo%rst(cvidx)%ai(3)

    ! calculate rji and rjk
    rji(:) = geo%crd(:,i) - geo%crd(:,j)
    rjk(:) = geo%crd(:,k) - geo%crd(:,j)

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

    ! calculate energy
    da = angv - geo%rst(cvidx)%trg_value
    geo%rst(cvidx)%value = 0.5*ANG_FC*da**2
    geo%rst_energy = geo%rst_energy + geo%rst(cvidx)%value

    ! calculate gradient
    dv = ANG_FC*da

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

    geo%grd(:,i) = geo%grd(:,i) + dv*di(:)
    geo%grd(:,k) = geo%grd(:,k) + dv*dk(:)
    geo%grd(:,j) = geo%grd(:,j) - dv*( di(:) + dk(:) )

end subroutine ffdev_geometry_get_angle_penalty

! ==============================================================================
! subroutine ffdev_geometry_get_torsion_penalty
! ==============================================================================

subroutine ffdev_geometry_get_torsion_penalty(geo,cvidx)

    implicit none
    type(GEOMETRY)  :: geo
    integer         :: cvidx
    ! ------------------------------
    integer         :: i,j,k,l
    real(DEVDP)     :: f(3),g(3),h(3),dv,fg,hg,a2,b2,gv,scp,phi,a(3),b(3),da
    ! --------------------------------------------------------------------------

    ! source: 10.1002/(SICI)1096-987X(19960715)17:9<1132::AID-JCC5>3.0.CO;2-T

    i  = geo%rst(cvidx)%ai(1)
    j  = geo%rst(cvidx)%ai(2)
    k  = geo%rst(cvidx)%ai(3)
    l  = geo%rst(cvidx)%ai(4)

    f(:) = geo%crd(:,i) - geo%crd(:,j)
    g(:) = geo%crd(:,j) - geo%crd(:,k)
    h(:) = geo%crd(:,l) - geo%crd(:,k)

    a(1) = f(2)*g(3) - f(3)*g(2)
    a(2) = f(3)*g(1) - f(1)*g(3)
    a(3) = f(1)*g(2) - f(2)*g(1)

    b(1) = h(2)*g(3) - h(3)*g(2)
    b(2) = h(3)*g(1) - h(1)*g(3)
    b(3) = h(1)*g(2) - h(2)*g(1)

    fg = dot_product(f,g)
    hg = dot_product(h,g)
    a2 = a(1)**2 + a(2)**2 + a(3)**2
    b2 = b(1)**2 + b(2)**2 + b(3)**2
    gv = sqrt( g(1)**2 + g(2)**2 + g(3)**2 )

    ! calculate scp and phi
    scp = (a(1)*b(1)+a(2)*b(2)+a(3)*b(3))/sqrt(a2*b2)
    if ( scp .gt.  1.0 ) then
            scp =  1.0
            phi = acos(1.0)    ! const
    else if ( scp .lt. -1.0 ) then
            scp = -1.0
            phi = acos(-1.0)   ! const
    else
        phi = acos( scp )
    end if
    if( g(1)*(a(2)*b(3)-a(3)*b(2)) &
       +g(2)*(a(3)*b(1)-a(1)*b(3)) &
       +g(3)*(a(1)*b(2)-a(2)*b(1)) .gt. 0) then
                phi = -phi
    end if

    ! deviation and value
    da = ffdev_geometry_get_dihedral_deviation(phi,geo%rst(cvidx)%trg_value)
    geo%rst(cvidx)%value = 0.5*ANG_FC*da**2
    geo%rst_energy = geo%rst_energy + geo%rst(cvidx)%value

    ! contribution to grad
    dv = ANG_FC*da

    ! calculate gradient
    geo%grd(:,i) = geo%grd(:,i) + dv*( -gv/a2*a(:) )
    geo%grd(:,j) = geo%grd(:,j) + dv*(  (gv/a2 + fg/(a2*gv))*a(:) - hg/(b2*gv)*b(:) )
    geo%grd(:,k) = geo%grd(:,k) + dv*(  (hg/(b2*gv) - gv/b2)*b(:) - fg/(a2*gv)*a(:) )
    geo%grd(:,l) = geo%grd(:,l) + dv*( gv/b2*b(:) )

end subroutine ffdev_geometry_get_torsion_penalty
        
! ------------------------------------------------------------------------------

end module ffdev_geometry

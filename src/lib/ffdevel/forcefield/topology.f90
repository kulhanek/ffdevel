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

module ffdev_topology

use ffdev_topology_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_topology_init
! ==============================================================================

subroutine ffdev_topology_init(top)

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------------------------------------

    top%name = ''
    top%natoms = 0
    top%natom_types = 0
    top%nbonds = 0
    top%nbond_types = 0
    top%nangles = 0
    top%nangle_types = 0
    top%ndihedrals = 0
    top%ndihedral_types = 0
    top%ndihedral_seq_size = 0
    top%nimpropers = 0
    top%nimproper_types = 0
    top%nb_size = 0
    top%probe_size = 0
    top%nfragments = 0
    top%sapt_size = 0
    top%nsymm_classes = 0
    top%total_charge = 0
    top%nb_params_update = .true.

end subroutine ffdev_topology_init

! ==============================================================================
! subroutine ffdev_topology_load
! ==============================================================================

subroutine ffdev_topology_load(top,name)

    use prmfile
    use prmfile_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)              :: top
    character(*)                :: name
    ! --------------------------------------------
    type(PRMFILE_TYPE)          :: fin
    logical                     :: my_result, lbuff
    integer                     :: alloc_stat, io_stat, i, idx, nbuff, pn
    character(PRMFILE_MAX_LINE) :: buffer
    real(DEVDP)                 :: v,g,c,p,w,e,chrg
    ! --------------------------------------------------------------------------

    ! load topology file
    call prmfile_init(fin)
    if( .not. prmfile_read(fin,name) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to read topology "'//trim(name)//'"!')
    end if

    top%name = name

    ! read dimmensions
    if( .not. prmfile_open_section(fin,'dimensions') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [dimensions] section!')
    end if

    my_result = .true.
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'atoms',top%natoms)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'atom_types',top%natom_types)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'bonds',top%nbonds)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'bond_types',top%nbond_types)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'angles',top%nangles)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'angle_types',top%nangle_types)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'dihedrals',top%ndihedrals)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'dihedral_types',top%ndihedral_types)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'dihedral_seq_size',top%ndihedral_seq_size)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'impropers',top%nimpropers)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'improper_types',top%nimproper_types)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'nb_size',top%nb_size)
    my_result = my_result .and. prmfile_get_integer_by_key(fin,'nb_types',top%nnb_types)

    ! optional items
    lbuff = prmfile_get_integer_by_key(fin,'nb_sizeij',nbuff)
    lbuff = prmfile_get_integer_by_key(fin,'nb_size14',nbuff)

    if( .not. my_result ) then
        call ffdev_utils_exit(DEV_ERR,1,'Missing data in [dimmensions] section!')
    end if

    ! allocate arrays
    allocate( top%atoms(top%natoms), top%bonds(top%nbonds), top%bond_types(top%nbond_types), &
              top%angles(top%nangles), top%angle_types(top%nangle_types), &
              top%dihedrals(top%ndihedrals), top%dihedral_types(top%ndihedral_types), &
              top%impropers(top%nimpropers), top%improper_types(top%nimproper_types), &
              top%nb_list(top%nb_size), top%atom_types(top%natom_types), &
              top%nb_types(top%nnb_types), stat = alloc_stat )

    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate topology arrays in ffdev_topology_load!')
    end if

    ! read sections

    ! read atoms -------------------------------------
    if( .not. prmfile_open_section(fin,'atoms') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [atoms] section!')
    end if

    chrg = 0.0
    do i=1,top%natoms
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [atoms] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%atoms(i)%typeid, top%atoms(i)%name, &
                                      top%atoms(i)%residx,top%atoms(i)%resname, top%atoms(i)%charge, top%atoms(i)%symmclass
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [atoms] section!')
        end if
        if( (top%atoms(i)%typeid .le. 0) .or. (top%atoms(i)%typeid .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [atoms] section!')
        end if
        if( (top%atoms(i)%symmclass .le. 0) .or. (top%atoms(i)%symmclass .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Symmetry class out-of-legal range in [atoms] section!')
        end if
        if( top%atoms(i)%symmclass .gt. top%nsymm_classes ) then
            top%nsymm_classes = top%atoms(i)%symmclass
        end if
        top%atoms(i)%nbonds = 0
        top%atoms(i)%frgid = 0
        chrg = chrg + top%atoms(i)%charge
    end do
    top%total_charge = NINT(chrg)

    ! read types -------------------------------------
    if( .not. prmfile_open_section(fin,'atom_types') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [atom_types] section!')
    end if

    do i=1,top%natom_types
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [atom_types] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%atom_types(i)%name, top%atom_types(i)%mass, &
                                      top%atom_types(i)%z
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [atom_types] section!')
        end if
        top%atom_types(i)%probe = .false.
    end do

    ! read bonds -------------------------------------
    if( .not. prmfile_open_section(fin,'bonds') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [bonds] section!')
    end if

    do i=1,top%nbonds
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [bonds] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%bonds(i)%ai, top%bonds(i)%aj, &
                                      top%bonds(i)%bt
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [bonds] section!')
        end if
        if( (top%bonds(i)%ai .le. 0) .or. (top%bonds(i)%ai .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [bonds] section!')
        end if
        if( (top%bonds(i)%aj .le. 0) .or. (top%bonds(i)%aj .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [bonds] section!')
        end if
        if( (top%bonds(i)%bt .le. 0) .or. (top%bonds(i)%bt .gt. top%nbond_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Bond type out-of-legal range in [bonds] section!')
        end if
    end do

    ! read bond types -------------------------------------
    if( .not. prmfile_open_section(fin,'bond_types') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [bond_types] section!')
    end if

    do i=1,top%nbond_types
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [bond_types] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%bond_types(i)%ti, top%bond_types(i)%tj, &
                                      top%bond_types(i)%model, &
                                      top%bond_types(i)%d0, top%bond_types(i)%k
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [bond_types] section!')
        end if
        if( (top%bond_types(i)%ti .le. 0) .or. (top%bond_types(i)%ti .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [bond_types] section!')
        end if
        if( (top%bond_types(i)%tj .le. 0) .or. (top%bond_types(i)%tj .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [bond_types] section!')
        end if
        top%bond_types(i)%ffoptactive = .false.
    end do

    ! read angles -------------------------------------
    if( .not. prmfile_open_section(fin,'angles') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [angles] section!')
    end if

    do i=1,top%nangles
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [angles] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%angles(i)%ai, top%angles(i)%aj, &
                                      top%angles(i)%ak, top%angles(i)%at
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [angles] section!')
        end if
        if( (top%angles(i)%ai .le. 0) .or. (top%angles(i)%ai .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [angles] section!')
        end if
        if( (top%angles(i)%aj .le. 0) .or. (top%angles(i)%aj .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [angles] section!')
        end if
        if( (top%angles(i)%ak .le. 0) .or. (top%angles(i)%ak .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [angles] section!')
        end if
        if( (top%angles(i)%at .le. 0) .or. (top%angles(i)%at .gt. top%nangle_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Angle type out-of-legal range in [angles] section!')
        end if
    end do

    ! read angle types -------------------------------------
    if( .not. prmfile_open_section(fin,'angle_types') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [angle_types] section!')
    end if

    do i=1,top%nangle_types
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [angle_types] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%angle_types(i)%ti, top%angle_types(i)%tj, &
                                      top%angle_types(i)%tk, top%angle_types(i)%model, &
                                      top%angle_types(i)%a0, top%angle_types(i)%k
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [angle_types] section!')
        end if
        if( (top%angle_types(i)%ti .le. 0) .or. (top%angle_types(i)%ti .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [angle_types] section!')
        end if
        if( (top%angle_types(i)%tj .le. 0) .or. (top%angle_types(i)%tj .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [angle_types] section!')
        end if
        if( (top%angle_types(i)%tk .le. 0) .or. (top%angle_types(i)%tk .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [angle_types] section!')
        end if
        top%angle_types(i)%ffoptactive = .false.
    end do

    ! read dihedrals -------------------------------------
    if( .not. prmfile_open_section(fin,'dihedrals') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [dihedrals] section!')
    end if

    do i=1,top%ndihedrals
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [dihedrals] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%dihedrals(i)%ai, top%dihedrals(i)%aj, &
                                      top%dihedrals(i)%ak, top%dihedrals(i)%al, top%dihedrals(i)%dt
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedrals] section!')
        end if
        if( (top%dihedrals(i)%ai .le. 0) .or. (top%dihedrals(i)%ai .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [dihedrals] section!')
        end if
        if( (top%dihedrals(i)%aj .le. 0) .or. (top%dihedrals(i)%aj .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [dihedrals] section!')
        end if
        if( (top%dihedrals(i)%ak .le. 0) .or. (top%dihedrals(i)%ak .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [dihedrals] section!')
        end if
        if( (top%dihedrals(i)%al .le. 0) .or. (top%dihedrals(i)%al .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [dihedrals] section!')
        end if
        if( (top%dihedrals(i)%dt .le. 0) .or. (top%dihedrals(i)%dt .gt. top%ndihedral_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Dihedral type out-of-legal range in [dihedrals] section!')
        end if
    end do

    ! read dihedral types -------------------------------------
    if( .not. prmfile_open_section(fin,'dihedral_types') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [dihedral_types] section!')
    end if

    ! allocate
    do i=1,top%ndihedral_types
        top%dihedral_types(i)%n = top%ndihedral_seq_size
        allocate(top%dihedral_types(i)%v(top%dihedral_types(i)%n), &
                 top%dihedral_types(i)%g(top%dihedral_types(i)%n), &
                 top%dihedral_types(i)%c(top%dihedral_types(i)%n), &
                 top%dihedral_types(i)%p(top%dihedral_types(i)%n), &
                 top%dihedral_types(i)%w2(top%dihedral_types(i)%n), &
                 top%dihedral_types(i)%enabled(top%dihedral_types(i)%n), stat = alloc_stat )
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate dihedral arrays!')
        end if
        top%dihedral_types(i)%v(:) = 0.0d0
        top%dihedral_types(i)%g(:) = 0.0d0
        top%dihedral_types(i)%c(:) = 0.0d0
        top%dihedral_types(i)%p(:) = 0.0d0
        top%dihedral_types(i)%w2(:) = 0.0d0
        top%dihedral_types(i)%enabled(:) = .true.
        top%dihedral_types(i)%ffoptactive = .false.
    end do

    do i=1,top%ndihedral_types
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [dihedral_types] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, nbuff, nbuff, nbuff, nbuff
        if( io_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_types] section!')
        end if
        if( (idx .le. 0) .or. (idx .gt. top%ndihedral_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_types] section!')
        end if

        read(buffer,*,iostat=io_stat) nbuff, top%dihedral_types(idx)%ti, top%dihedral_types(idx)%tj, &
                                      top%dihedral_types(idx)%tk, top%dihedral_types(idx)%tl, &
                                      top%dihedral_types(idx)%mode, &
                                      top%dihedral_types(idx)%inv_scee, top%dihedral_types(idx)%inv_scnb

        if( top%dihedral_types(idx)%inv_scee .ne. 0.0 ) then
            top%dihedral_types(idx)%inv_scee = 1.0d0 / top%dihedral_types(idx)%inv_scee
        end if
        if( top%dihedral_types(idx)%inv_scnb .ne. 0.0 ) then
            top%dihedral_types(idx)%inv_scnb = 1.0d0 / top%dihedral_types(idx)%inv_scnb
        end if

        if( io_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_types] section!')
        end if
        if( (top%dihedral_types(idx)%ti .le. 0) .or. (top%dihedral_types(idx)%ti .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [dihedral_types] section!')
        end if
        if( (top%dihedral_types(idx)%tj .le. 0) .or. (top%dihedral_types(idx)%tj .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [dihedral_types] section!')
        end if
        if( (top%dihedral_types(idx)%tk .le. 0) .or. (top%dihedral_types(idx)%tk .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [dihedral_types] section!')
        end if
        if( (top%dihedral_types(idx)%tl .le. 0) .or. (top%dihedral_types(idx)%tl .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [dihedral_types] section!')
        end if
    end do

    if( prmfile_open_section(fin,'dihedral_seq_cos') ) then
        do while( prmfile_get_line(fin,buffer) )
            read(buffer,*,iostat=io_stat) idx, pn
            if( io_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_cos] section!')
            end if
            if( (idx .le. 0) .or. (idx .gt. top%ndihedral_types) ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_cos] section!')
            end if
            if( (pn .le. 0) .or. (pn .gt. top%ndihedral_seq_size) ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_cos] section!')
            end if
            if( top%dihedral_types(idx)%mode .ne. DIH_COS ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_cos] section!')
            end if

            e = 0
            read(buffer,*,iostat=io_stat) idx, pn, v, g, e
            if( io_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_cos] section!')
            end if

            top%dihedral_types(idx)%v(pn) = v
            top%dihedral_types(idx)%g(pn) = g
            if( e .eq. 1 ) then
                top%dihedral_types(idx)%enabled(pn) = .true.
            else
                top%dihedral_types(idx)%enabled(pn) = .false.
            end if
        end do
    end if

    if( prmfile_open_section(fin,'dihedral_seq_grbf') ) then
        do while( prmfile_get_line(fin,buffer) )
            read(buffer,*,iostat=io_stat) idx, pn
            if( io_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_grbf] section!')
            end if
            if( (idx .le. 0) .or. (idx .gt. top%ndihedral_types) ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_grbf] section!')
            end if
            if( (pn .le. 0) .or. (pn .gt. top%ndihedral_seq_size) ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_grbf] section!')
            end if
            if( top%dihedral_types(idx)%mode .ne. DIH_GRBF ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_grbf] section!')
            end if

            read(buffer,*,iostat=io_stat) idx, pn, c, p, w
            if( io_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [dihedral_seq_grbf] section!')
            end if

            top%dihedral_types(idx)%c(pn) = c
            top%dihedral_types(idx)%p(pn) = p
            top%dihedral_types(idx)%w2(pn) = w
        end do
    end if


    ! read impropers -------------------------------------
    if( .not. prmfile_open_section(fin,'impropers') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [impropers] section!')
    end if

    do i=1,top%nimpropers
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [impropers] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%impropers(i)%ai, top%impropers(i)%aj, &
                                      top%impropers(i)%ak, top%impropers(i)%al, top%impropers(i)%dt
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [impropers] section!')
        end if
        if( (top%impropers(i)%ai .le. 0) .or. (top%impropers(i)%ai .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [impropers] section!')
        end if
        if( (top%impropers(i)%aj .le. 0) .or. (top%impropers(i)%aj .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [impropers] section!')
        end if
        if( (top%impropers(i)%ak .le. 0) .or. (top%impropers(i)%ak .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [impropers] section!')
        end if
        if( (top%impropers(i)%al .le. 0) .or. (top%impropers(i)%al .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [impropers] section!')
        end if
        if( (top%impropers(i)%dt .le. 0) .or. (top%impropers(i)%dt .gt. top%nimproper_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Bond type out-of-legal range in [impropers] section!')
        end if
    end do

    ! read improper types -------------------------------------
    if( .not. prmfile_open_section(fin,'improper_types') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [improper_types] section!')
    end if

    do i=1,top%nimproper_types
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [improper_types] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%improper_types(i)%ti, top%improper_types(i)%tj, &
                                      top%improper_types(i)%tk, top%improper_types(i)%tl, &
                                      top%improper_types(i)%v, top%improper_types(i)%g
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [improper_types] section!')
        end if
        if( (top%improper_types(i)%ti .le. 0) .or. (top%improper_types(i)%ti .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [improper_types] section!')
        end if
        if( (top%improper_types(i)%tj .le. 0) .or. (top%improper_types(i)%tj .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [improper_types] section!')
        end if
        if( (top%improper_types(i)%tk .le. 0) .or. (top%improper_types(i)%tk .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [improper_types] section!')
        end if
        if( (top%improper_types(i)%tl .le. 0) .or. (top%improper_types(i)%tl .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [improper_types] section!')
        end if

        top%improper_types(i)%ffoptactive = .false.
    end do

    ! read NB list -------------------------------------
    if( .not. prmfile_open_section(fin,'nb_list') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [nb_list] section!')
    end if

    do i=1,top%nb_size
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [nb_list] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%nb_list(i)%ai, top%nb_list(i)%aj, &
                                      top%nb_list(i)%nbt, top%nb_list(i)%dt
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [nb_list] section!')
        end if
        if( (top%nb_list(i)%ai .le. 0) .or. (top%nb_list(i)%ai .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [nb_list] section!')
        end if
        if( (top%nb_list(i)%aj .le. 0) .or. (top%nb_list(i)%aj .gt. top%natoms) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom index out-of-legal range in [nb_list] section!')
        end if
        if( (top%nb_list(i)%nbt .lt. 0) .or. (top%nb_list(i)%nbt .gt. top%nnb_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'NB type out-of-legal range in [nb_list] section!')
        end if
        ! it includes zero, which is not 1,4 interaction
        if( (top%nb_list(i)%dt .lt. 0) .or. (top%nb_list(i)%dt .gt. top%ndihedral_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Dihedral type out-of-legal range in [nb_list] section!')
        end if

       top%nb_list(i)%crgij = 0.0d0
       top%nb_list(i)%pa    = 0.0d0
       top%nb_list(i)%pb    = 0.0d0
       top%nb_list(i)%c6    = 0.0d0
       top%nb_list(i)%c8    = 0.0d0
       top%nb_list(i)%c10   = 0.0d0
       top%nb_list(i)%rc6   = 0.0d0
       top%nb_list(i)%rc8   = 0.0d0
       top%nb_list(i)%rc10  = 0.0d0
       top%nb_list(i)%tb    = 0.0d0
    end do

    ! read NB types ------------------------------------
    if( .not. prmfile_open_section(fin,'nb_types') ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to open [nb_types] section!')
    end if

    do i=1,top%nnb_types
        if( .not. prmfile_get_line(fin,buffer) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Premature end of [nb_types] section!')
        end if
        read(buffer,*,iostat=io_stat) idx, top%nb_types(i)%ti, top%nb_types(i)%tj, &
                                      top%nb_types(i)%eps, top%nb_types(i)%r0
        if( (idx .ne. i) .or. (io_stat .ne. 0) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Illegal record in [nb_types] section!')
        end if

        top%nb_types(i)%pa = 0.0d0
        top%nb_types(i)%pb = 0.0d0
        top%nb_types(i)%rc = 0.0d0
        top%nb_types(i)%tb = 0.0d0

        if( (top%nb_types(i)%ti .le. 0) .or. (top%nb_types(i)%ti .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [nb_types] section!')
        end if
        if( (top%nb_types(i)%tj .le. 0) .or. (top%nb_types(i)%tj .gt. top%natom_types) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Atom type out-of-legal range in [nb_types] section!')
        end if
        top%nb_types(i)%ffoptactive = .false.
    end do

    ! check if everything was read
    if( prmfile_count_ulines(fin) .ne. 0 ) then
        write(DEV_OUT,*)
        call prmfile_dump(fin,DEV_OUT,.true.)
        call ffdev_utils_exit(DEV_ERR,1,'Unprocessed lines found in the topology file!')
    end if

    ! release the file
    call prmfile_clear(fin)

    ! create list of bonded atoms
    call ffdev_topology_gen_bonded(top)

    ! create fragments
    call ffdev_topology_gen_fragments(top)

    ! call init nbi and nbj in nb_pairs
    call ffdev_topology_init_nbij(top)

    ! we need to update nb parameters before they are used
    top%nb_params_update = .true.

    return

end subroutine ffdev_topology_load

! ==============================================================================
! subroutine ffdev_topology_save
! ==============================================================================

subroutine ffdev_topology_save(top,name)

    use prmfile
    use prmfile_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)      :: top
    character(*)        :: name
    ! --------------------------------------------
    character(len=20)   :: key
    integer             :: i,j,dm
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_TOP,name,'U')

    ! write header section
    write(DEV_TOP,10) 'dimensions'
    key = 'atoms'
    write(DEV_TOP,20) adjustl(key), top%natoms
    key = 'atom_types'
    write(DEV_TOP,20) adjustl(key), top%natom_types
    key = 'bonds'
    write(DEV_TOP,20) adjustl(key), top%nbonds
    key = 'bond_types'
    write(DEV_TOP,20) adjustl(key), top%nbond_types
    key = 'angles'
    write(DEV_TOP,20) adjustl(key), top%nangles
    key = 'angle_types'
    write(DEV_TOP,20) adjustl(key), top%nangle_types
    key = 'dihedrals'
    write(DEV_TOP,20) adjustl(key), top%ndihedrals
    key = 'dihedral_types'
    write(DEV_TOP,20) adjustl(key), top%ndihedral_types
    key = 'dihedral_seq_size'
    write(DEV_TOP,20) adjustl(key), top%ndihedral_seq_size
    key = 'impropers'
    write(DEV_TOP,20) adjustl(key), top%nimpropers
    key = 'improper_types'
    write(DEV_TOP,20) adjustl(key), top%nimproper_types
    key = 'nb_size'
    write(DEV_TOP,20) adjustl(key), top%nb_size
    key = 'nb_types'
    write(DEV_TOP,20) adjustl(key), top%nnb_types

 10 format('[',A,']')
 20 format(A20,1X,I8)

    ! atoms ----------------------------
    write(DEV_TOP,10) 'atom_types'
    do i=1,top%natom_types
        write(DEV_TOP,30)   i, adjustl(top%atom_types(i)%name), &
                               top%atom_types(i)%mass, &
                               top%atom_types(i)%z
        top%atom_types(i)%probe = .false.
    end do
 30 format(I7,1X,A4,1X,F10.4,1X,I2)

    write(DEV_TOP,10) 'atoms'
    do i=1,top%natoms
        write(DEV_TOP,40)   i, top%atoms(i)%typeid, &
                               adjustl(top%atoms(i)%name), &
                               top%atoms(i)%residx, &
                               adjustl(top%atoms(i)%resname), &
                               top%atoms(i)%charge
    end do
 40 format(I7,1X,I4,1X,A4,1X,I5,1X,A4,1X,F10.6)


    ! bonds ----------------------------

    write(DEV_TOP,10) 'bond_types'
    do i=1,top%nbond_types
        write(DEV_TOP,50)   i, top%bond_types(i)%ti, &
                               top%bond_types(i)%tj, &
                               top%bond_types(i)%model, &
                               top%bond_types(i)%d0, top%bond_types(i)%k
    end do
 50 format(I7,1X,I5,1X,I5,1X,I4,1X,F13.6,1X,F13.6)

    write(DEV_TOP,10) 'bonds'
    do i=1,top%nbonds
        write(DEV_TOP,60)   i, top%bonds(i)%ai, top%bonds(i)%aj, top%bonds(i)%bt
    end do
 60 format(I7,1X,I5,1X,I5,1X,I4)

    ! angles ----------------------------
    write(DEV_TOP,10) 'angle_types'
    do i=1,top%nangle_types
        write(DEV_TOP,70)   i, top%angle_types(i)%ti, &
                               top%angle_types(i)%tj, &
                               top%angle_types(i)%tk, &
                               top%angle_types(i)%model, &
                               top%angle_types(i)%a0, top%angle_types(i)%k
    end do
 70 format(I7,1X,I5,1X,I5,1X,I5,1X,I4,1X,F13.6,1X,F13.6)

    write(DEV_TOP,10) 'angles'
    do i=1,top%nangles
        write(DEV_TOP,80)   i, top%angles(i)%ai, top%angles(i)%aj, top%angles(i)%ak, top%angles(i)%at
    end do
 80 format(I7,1X,I5,1X,I5,1X,I5,1X,I4)

    ! dihedrals -------------------------
    write(DEV_TOP,10) 'dihedral_types'
    do i=1,top%ndihedral_types
        write(DEV_TOP,90)   i, top%dihedral_types(i)%ti, &
                               top%dihedral_types(i)%tj, &
                               top%dihedral_types(i)%tk, &
                               top%dihedral_types(i)%tl, &
                               top%dihedral_types(i)%mode,  &
                               1.0/top%dihedral_types(i)%inv_scee, 1.0/top%dihedral_types(i)%inv_scnb
    end do

 90 format(I7,1X,I5,1X,I5,1X,I5,1X,I5,1X,I4,1X,F13.6,1X,F13.6,1X,F13.6,1X,F13.6)

    ! dihedrals -------------------------
    write(DEV_TOP,10) 'dihedral_seq_cos'
    do i=1,top%ndihedral_types
        if( top%dihedral_types(i)%mode .eq. DIH_COS ) then
            do j=1,top%dihedral_types(i)%n
                dm = 0
                if( top%dihedral_types(i)%enabled(j) ) dm = 1
                write(DEV_TOP,92)   i, j, &
                                       top%dihedral_types(i)%v(j), &
                                       top%dihedral_types(i)%g(j), &
                                       dm
            end do
        end if
    end do

 92 format(I6,1X,I2,1X,F13.6,1X,F13.6,1X,I7)

    ! dihedrals -------------------------
    write(DEV_TOP,10) 'dihedral_seq_grbf'
    do i=1,top%ndihedral_types
        if( top%dihedral_types(i)%mode .eq. DIH_GRBF ) then
            do j=1,top%dihedral_types(i)%n
                write(DEV_TOP,93)   i, j, &
                                       top%dihedral_types(i)%c(j), &
                                       top%dihedral_types(i)%p(j), &
                                       top%dihedral_types(i)%w2(j)
            end do
        end if
    end do

 93 format(I6,1X,I2,1X,F13.6,1X,F13.6,1X,F13.6)

    write(DEV_TOP,10) 'dihedrals'
    do i=1,top%ndihedrals
        write(DEV_TOP,100)   i, top%dihedrals(i)%ai, top%dihedrals(i)%aj, &
                               top%dihedrals(i)%ak, top%dihedrals(i)%al, top%dihedrals(i)%dt
    end do
100 format(I7,1X,I5,1X,I5,1X,I5,1X,I5,1X,I4)

    ! impropers ------------------------
    write(DEV_TOP,10) 'improper_types'
    do i=1,top%nimproper_types
        write(DEV_TOP,110)   i, top%improper_types(i)%ti, &
                               top%improper_types(i)%tj, &
                               top%improper_types(i)%tk, &
                               top%improper_types(i)%tl, &
                               top%improper_types(i)%v, top%improper_types(i)%g
    end do

110 format(I7,1X,I5,1X,I5,1X,I5,1X,I5,1X,F13.6,1X,F13.6)

    write(DEV_TOP,10) 'impropers'
    do i=1,top%nimpropers
        write(DEV_TOP,120)   i, top%impropers(i)%ai, top%impropers(i)%aj, &
                               top%impropers(i)%ak, top%impropers(i)%al, top%impropers(i)%dt
    end do
120 format(I7,1X,I5,1X,I5,1X,I5,1X,I5,1X,I4)

    ! NB list --------------------------
    write(DEV_TOP,10) 'nb_types'
    do i=1,top%nnb_types
        write(DEV_TOP,125)   i, top%nb_types(i)%ti, top%nb_types(i)%tj, &
                               top%nb_types(i)%eps, top%nb_types(i)%r0
    end do
125 format(I7,1X,I5,1X,I5,1X,F13.7,1X,F13.7)

    write(DEV_TOP,10) 'nb_list'
    do i=1,top%nb_size
        write(DEV_TOP,130)   i, top%nb_list(i)%ai, top%nb_list(i)%aj, &
                               top%nb_list(i)%nbt, top%nb_list(i)%dt
    end do
130 format(I7,1X,I5,1X,I5,1X,I5,1X,I5)

    close(DEV_TOP)

end subroutine ffdev_topology_save

! ==============================================================================
! subroutine ffdev_topology_info
! ==============================================================================

subroutine ffdev_topology_info(top)

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------------------------------------

    write(DEV_OUT,90)  trim(top%name)
    write(DEV_OUT,100) top%natoms
    write(DEV_OUT,110) top%natom_types
    if( top%probe_size .eq. 0 ) then
    write(DEV_OUT,120) top%nbonds
    write(DEV_OUT,130) top%nbond_types
    write(DEV_OUT,140) top%nangles
    write(DEV_OUT,150) top%nangle_types
    write(DEV_OUT,160) top%ndihedrals
    write(DEV_OUT,170) top%ndihedral_types
    write(DEV_OUT,180) top%ndihedral_seq_size
    write(DEV_OUT,190) top%nimpropers
    write(DEV_OUT,200) top%nimproper_types
    end if
    write(DEV_OUT,210) top%nb_size
    write(DEV_OUT,230) top%nfragments
    if( top%probe_size .ne. 0 ) then
    write(DEV_OUT,220) top%probe_size
    end if

 90 format('Topology name                      = ',A)
100 format('Number of atoms                    = ',I6)
110 format('Number of types                    = ',I6)
120 format('Number of bonds                    = ',I6)
130 format('Number of bond types               = ',I6)
140 format('Number of angles                   = ',I6)
150 format('Number of angle types              = ',I6)
160 format('Number of dihedrals                = ',I6)
170 format('Number of dihedral types           = ',I6)
180 format('Number of dihedral seq size        = ',I6)
190 format('Number of impropers                = ',I6)
200 format('Number of improper types           = ',I6)
210 format('Number of NB size                  = ',I6)
220 format('Number of atoms in probe           = ',I6)
230 format('Number of fragments                = ',I6)

end subroutine ffdev_topology_info

! ==============================================================================
! subroutine ffdev_topology_info_types
! ==============================================================================

subroutine ffdev_topology_info_types(top,mode)

    use smf_periodic_table_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)      :: top
    integer,optional    :: mode
    ! --------------------------------------------
    integer         :: i,j,rmode
    logical         :: print
    real(DEVDP)     :: scee,scnb
    ! --------------------------------------------------------------------------

    rmode = 0 ! print all
    if( present(mode) ) then
        rmode = mode
    end if

    if( (rmode .eq. 0) .or. (rmode .eq. 1) ) then
        write(DEV_OUT,5) trim(top%name)
    end if

    if( (rmode .eq. 0) .or. (rmode .eq. 1) ) then

! atoms ----------------------------
        write(DEV_OUT,*)
        write(DEV_OUT,10)
        write(DEV_OUT,20)
        write(DEV_OUT,30)
        do i=1,top%natom_types
            write(DEV_OUT,40)   i, adjustl(top%atom_types(i)%name), &
                                   top%atom_types(i)%z, adjustl(pt_symbols(top%atom_types(i)%z)), &
                                   top%atom_types(i)%mass
        end do
    end if

    if( rmode .eq. 0 ) then
! bonds ----------------------------
        if( top%nbond_types .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,110)
            write(DEV_OUT,120)
            write(DEV_OUT,130)
            do i=1,top%nbond_types
                write(DEV_OUT,140)   i, adjustl(top%atom_types(top%bond_types(i)%ti)%name), &
                                       adjustl(top%atom_types(top%bond_types(i)%tj)%name), &
                                       top%bond_types(i)%model, &
                                       top%bond_types(i)%d0, top%bond_types(i)%k
            end do
        end if

! angles ----------------------------
        if( top%nangle_types .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,210)
            write(DEV_OUT,220)
            write(DEV_OUT,230)
            do i=1,top%nangle_types
                write(DEV_OUT,240)   i, adjustl(top%atom_types(top%angle_types(i)%ti)%name), &
                                       adjustl(top%atom_types(top%angle_types(i)%tj)%name), &
                                       adjustl(top%atom_types(top%angle_types(i)%tk)%name), &
                                       top%angle_types(i)%a0*DEV_R2D, top%angle_types(i)%k
            end do
        end if

! dihedrals -------------------------
        if( top%ndihedral_types .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,310)
            write(DEV_OUT,320)
            write(DEV_OUT,330)
            do i=1,top%ndihedral_types
                scee = 0.0
                if( top%dihedral_types(i)%inv_scee .ne. 0 ) then
                    scee = 1.0/top%dihedral_types(i)%inv_scee
                end if
                if( top%dihedral_types(i)%inv_scnb .ne. 0 ) then
                    scnb = 1.0/top%dihedral_types(i)%inv_scnb
                end if
                write(DEV_OUT,340)   i, adjustl(top%atom_types(top%dihedral_types(i)%ti)%name), &
                                       adjustl(top%atom_types(top%dihedral_types(i)%tj)%name), &
                                       adjustl(top%atom_types(top%dihedral_types(i)%tk)%name), &
                                       adjustl(top%atom_types(top%dihedral_types(i)%tl)%name), &
                                       top%dihedral_types(i)%mode, scee, scnb
            end do
        end if

! dihedrals cos mode -------------------------
        if( top%ndihedral_types .gt. 0 ) then
            print = .false.
            do i=1,top%ndihedral_types
                if( top%dihedral_types(i)%mode .eq. DIH_COS ) then
                    print = .true.
                    exit
                end if
            end do

            if( print ) then
                write(DEV_OUT,*)
                write(DEV_OUT,350)
                write(DEV_OUT,355)
                write(DEV_OUT,360)
                do i=1,top%ndihedral_types
                    if( top%dihedral_types(i)%mode .ne. DIH_COS ) cycle
                    do j=1,top%dihedral_types(i)%n
                        write(DEV_OUT,365)   i, top%dihedral_types(i)%enabled(j), j, top%dihedral_types(i)%v(j), &
                                               top%dihedral_types(i)%g(j)*DEV_R2D
                    end do
                end do
            end if
        end if

! dihedrals grbf mode -------------------------
        if( top%ndihedral_types .gt. 0 ) then
            print = .false.
            do i=1,top%ndihedral_types
                if( top%dihedral_types(i)%mode .eq. DIH_GRBF ) then
                    print = .true.
                    exit
                end if
            end do

            if( print ) then
                write(DEV_OUT,*)
                write(DEV_OUT,370)
                write(DEV_OUT,375)
                write(DEV_OUT,380)
                do i=1,top%ndihedral_types
                    if( top%dihedral_types(i)%mode .ne. DIH_GRBF ) cycle
                    do j=1,top%dihedral_types(i)%n
                        write(DEV_OUT,385)   i, top%dihedral_types(i)%enabled(j), j, top%dihedral_types(i)%c(j), &
                                               top%dihedral_types(i)%p(j)*DEV_R2D, sqrt(top%dihedral_types(i)%w2(j))*DEV_R2D
                    end do
                end do
            end if
        end if

! impropers ------------------------
        if( top%nimproper_types .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,410)
            write(DEV_OUT,420)
            write(DEV_OUT,430)
            do i=1,top%nimproper_types
                write(DEV_OUT,440)   i, adjustl(top%atom_types(top%improper_types(i)%ti)%name), &
                                       adjustl(top%atom_types(top%improper_types(i)%tj)%name), &
                                       adjustl(top%atom_types(top%improper_types(i)%tk)%name), &
                                       adjustl(top%atom_types(top%improper_types(i)%tl)%name), &
                                       top%improper_types(i)%v, top%improper_types(i)%g*DEV_R2D
            end do
        end if
    end if

    if( (rmode .eq. 0) .or. (rmode .eq. 1) .or. (rmode .eq. 2) ) then
! NB ------------------------
        if( top%nnb_types .gt. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,510)
            write(DEV_OUT,520)
            write(DEV_OUT,530)
            do i=1,top%nnb_types
                write(DEV_OUT,540)   i, adjustl(top%atom_types(top%nb_types(i)%ti)%name), &
                                       adjustl(top%atom_types(top%nb_types(i)%tj)%name), &
                                       top%nb_types(i)%eps, top%nb_types(i)%r0
            end do
        end if
    end if

  5 format('Topology name = ',A)

 10 format('# ~~~~~~~~~~~~~~~~~ atom types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
 20 format('# ID Type  Z  Sym    Mass   ')
 30 format('# -- ---- --- --- ----------')
 40 format(I4,1X,A4,1X,I3,1X,A3,1X,F10.4)

110 format('# ~~~~~~~~~~~~~~~~~ bond types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
120 format('# ID TypA TypB Mod       R0               K         ')
130 format('# -- ---- ---- --- ---------------- ----------------')
140 format(I4,1X,A4,1X,A4,1X,I3,1X,F16.6,1X,F16.6)

210 format('# ~~~~~~~~~~~~~~~~ angle types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
220 format('# ID TypA TypB TypC       A0               K         ')
230 format('# -- ---- ---- ---- ---------------- ----------------')
240 format(I4,1X,A4,1X,A4,1X,A4,1X,F16.6,1X,F16.6)

310 format('# ~~~~~~~~~~~~~~~~ dihedral types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
320 format('# ID TypA TypB TypC TypD Mod     scee       scnb    ')
330 format('# -- ---- ---- ---- ---- --- ----------- -----------')
340 format(I4,1X,A4,1X,A4,1X,A4,1X,A4,1X,I3,1X,F11.6,1X,F11.6)

350 format('# ~~~~~~~~~~~~~~~~ dihedral types - cos mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
355 format('# Type ST Pn     V           gamma    ')
360 format('# ---- -- -- ------------ ------------')
365 format(I6,1X,L2,1X,I2,1X,F12.6,1X,F12.6)

370 format('# ~~~~~~~~~~~~~~~~ dihedral types - grbf mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
375 format('# Type ST Pn       c           p             w     ')
380 format('# ---- -- -- ------------ ------------ ------------')
385 format(I6,1X,L2,1X,I2,1X,F12.6,1X,F12.6,1X,F12.6)

410 format('# ~~~~~~~~~~~~~~~~ improper types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
420 format('# ID TypA TypB TypC TypD         V             gamma      ')
430 format('# -- ---- ---- ---- ---- ---------------- ----------------')
440 format(I4,1X,A4,1X,A4,1X,A4,1X,A4,1X,F16.6,1X,F16.6)

510 format('# ~~~~~~~~~~~~~~~~~ NB types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

520 format('# ID TypA TypB        eps              R0       ')
530 format('# -- ---- ---- ---------------- ----------------')
540 format(I4,1X,A4,1X,A4,1X,F16.7,1X,F16.7,1X,F16.7)

end subroutine ffdev_topology_info_types

! ==============================================================================
! subroutine ffdev_topology_finalize_setup
! ==============================================================================

subroutine ffdev_topology_finalize_setup(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------------------------------------

    ! reserved for future usage

end subroutine ffdev_topology_finalize_setup

! ==============================================================================
! subroutine ffdev_topology_gen_bonded
! ==============================================================================

subroutine ffdev_topology_gen_bonded(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,j,k,alloc_stat
    ! --------------------------------------------------------------------------

    do i=1,top%natoms
        do j=1,top%nbonds
            if( (top%bonds(j)%ai .eq. i) .or. (top%bonds(j)%aj .eq. i) ) then
                top%atoms(i)%nbonds = top%atoms(i)%nbonds + 1
            end if
        end do
        if( top%atoms(i)%nbonds .gt. 0 ) then
            allocate(top%atoms(i)%bonded(top%atoms(i)%nbonds), stat = alloc_stat)
            if( alloc_stat .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for atom neighbours!')
            end if
            k = 0
            do j=1,top%nbonds
                if( (top%bonds(j)%ai .eq. i) .or. (top%bonds(j)%aj .eq. i) ) then
                    k = k + 1
                    if( top%bonds(j)%ai .eq. i ) then
                        top%atoms(i)%bonded(k) = top%bonds(j)%aj
                    else
                        top%atoms(i)%bonded(k) = top%bonds(j)%ai
                    end if
                end if
            end do
        end if
    end do

end subroutine ffdev_topology_gen_bonded

! ==============================================================================
! subroutine ffdev_topology_gen_fragments
! ==============================================================================

subroutine ffdev_topology_gen_fragments(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,j
    logical         :: changed,found
    ! --------------------------------------------------------------------------

    ! calculate number of fragments

    top%nfragments = 0

    do while (.true.)

        found = .false.
        ! mark first unprocessed atom
        do i=1,top%natoms
            if( top%atoms(i)%frgid .eq. 0 ) then
                top%nfragments = top%nfragments + 1
                top%atoms(i)%frgid = top%nfragments
                found = .true.
                exit
            end if
        end do

        if( .not. found ) exit

        !
        changed = .true.
        do while (changed)
            changed = .false.
            do i=1,top%natoms
                if( top%atoms(i)%frgid .eq. top%nfragments ) then
                    do j=1,top%atoms(i)%nbonds
                        if( top%atoms(top%atoms(i)%bonded(j))%frgid .eq. 0 ) then
                            top%atoms(top%atoms(i)%bonded(j))%frgid = top%nfragments
                            changed = .true.
                        end if
                    end do
                end if
            end do
        end do

    end do

end subroutine ffdev_topology_gen_fragments

! ==============================================================================
! subroutine ffdev_topology_init_nbij
! ==============================================================================

subroutine ffdev_topology_init_nbij(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,nbt,ti,tj
    ! --------------------------------------------------------------------------

    do i=1,top%nb_size
        nbt = top%nb_list(i)%nbt
        ti = top%nb_types(nbt)%ti
        tj = top%nb_types(nbt)%tj
        top%nb_list(i)%nbtii = ffdev_topology_find_nbtype_by_tindex(top,ti,ti)
        top%nb_list(i)%nbtjj = ffdev_topology_find_nbtype_by_tindex(top,tj,tj)
    end do

end subroutine ffdev_topology_init_nbij

! ==============================================================================
! function ffdev_topology_find_nbtype
! ==============================================================================

integer function ffdev_topology_find_nbtype_by_aindex(top,ai,aj)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ai,aj
    ! --------------------------------------------
    integer         :: ti,tj
    ! --------------------------------------------------------------------------

    ! convert to types
    ti = top%atoms(ai)%typeid
    tj = top%atoms(aj)%typeid

    ffdev_topology_find_nbtype_by_aindex = ffdev_topology_find_nbtype_by_tindex(top,ti,tj)

end function ffdev_topology_find_nbtype_by_aindex

! ==============================================================================
! function ffdev_topology_find_nbtype_by_tindex
! ==============================================================================

integer function ffdev_topology_find_nbtype_by_tindex(top,ti,tj)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: ti,tj
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,top%nnb_types
        if(  ( (top%nb_types(i)%ti .eq. ti) .and. (top%nb_types(i)%tj .eq. tj) ) .or. &
             ( (top%nb_types(i)%ti .eq. tj) .and. (top%nb_types(i)%tj .eq. ti) ) ) then
             ffdev_topology_find_nbtype_by_tindex = i
             return
        end if
    end do

    ! not found
    ffdev_topology_find_nbtype_by_tindex = 0
    return

end function ffdev_topology_find_nbtype_by_tindex

! ==============================================================================
! function ffdev_topology_find_nbtype_by_types
! ==============================================================================

integer function ffdev_topology_find_nbtype_by_types(top,sti,stj)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    character(*)    :: sti,stj
    ! --------------------------------------------
    integer         :: i,ti,tj
    ! --------------------------------------------------------------------------

    ti = 0
    tj = 0

    do i=1,top%natom_types
        if( top%atom_types(i)%name .eq. sti ) ti = i
        if( top%atom_types(i)%name .eq. stj ) tj = i
    end do

    ffdev_topology_find_nbtype_by_types = ffdev_topology_find_nbtype_by_tindex(top,ti,tj)

end function ffdev_topology_find_nbtype_by_types

! ==============================================================================
! function ffdev_topology_find_r_for_lj
! ==============================================================================

subroutine ffdev_topology_find_r_for_lj(eps,r0,e,r,er)

    use ffdev_utils

    implicit none
    real(DEVDP)     :: eps,r0,e,r,er
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r1,r2,a,b,r12,e12
    ! --------------------------------------------------------------------------

    ! bisection method

    a = eps*r0**12
    b = 2.0d0*eps*r0**6

    ! consider only repulsive part of function
    r1 = 1.0d0
    r2 = r0/2**(1.0d0/6.0d0)

    do i=1,1000
        r12 = (r1+r2)*0.5d0
        e12 = a/r12**12 - b/r12**6 - e
        if( e12 .gt. 0 ) then
            r1 = r12
        else
            r2 = r12
        end if
    end do

    r = (r1+r2)*0.5d0
    er = a/r**12 - b/r**6

!    write(*,*) eps,r0,e,r,er

end subroutine ffdev_topology_find_r_for_lj

! ==============================================================================
! subroutine ffdev_topology_LJ_ERA2ABC
! ==============================================================================

subroutine ffdev_topology_LJ_ERA2ABC(eps,r0,pa,c6,s6)

    implicit none
    real(DEVDP)     :: eps,r0,pa,c6,s6
    ! --------------------------------------------------------------------------

    pa = eps * r0**12
    c6 = 2.0d0 * eps * r0**6 / s6

end subroutine ffdev_topology_LJ_ERA2ABC

! ==============================================================================
! subroutine ffdev_topology_LJ_AB2ER
! ==============================================================================

subroutine ffdev_topology_LJ_AB2ER(pa,c6,eps,r0)

    implicit none
    real(DEVDP)     :: pa,c6,eps,r0
    ! --------------------------------------------
    real(DEVDP)     :: r6
    ! --------------------------------------------------------------------------

    r6  = 2.0d0*pa/c6
    r0  = r6**(1.0d0/6.0d0)
    eps = pa / (r6**2)

end subroutine ffdev_topology_LJ_AB2ER

! ==============================================================================
! function ffdev_topology_switch_to_probe_mode
! ==============================================================================

subroutine ffdev_topology_switch_to_probe_mode(top,probe_size,unique_probe_types)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: probe_size
    logical         :: unique_probe_types
    ! --------------------------------------------
    integer         :: i,j,ip,it,k,alloc_stat
    ! --------------------------------------------------------------------------

    if( probe_size .eq. 0 ) return

    top%probe_size = probe_size

    if( unique_probe_types ) then
        ! label atom types with probe
        do i=top%natoms - top%probe_size+1,top%natoms
            it = top%atoms(i)%typeid
            top%atom_types(it)%probe = .true.
        end do

        ! test probe atom type overlaps
        do i=1,top%natoms - top%probe_size
            it = top%atoms(i)%typeid
            if( top%atom_types(it)%probe ) then
                call ffdev_utils_exit(DEV_ERR,1,'Atom type of probe cannot be used in probed structure!')
            end if
        end do
    end if

    ! check covalent bonds between probe and probed structure
    do i=1,top%natoms - top%probe_size
        do j=top%natoms - top%probe_size+1,top%natoms
            do k=1,top%nbonds
                if( (top%bonds(k)%ai .eq. i) .and. (top%bonds(k)%aj .eq. j) .or. &
                    (top%bonds(k)%ai .eq. j) .and. (top%bonds(k)%aj .eq. i) ) then
                    call ffdev_utils_exit(DEV_ERR,1,'Covalent bond detected between the probe and probed structure!')
                end if
            end do
        end do
    end do

! switch to probe mode
    ! release previous NB list
    if( associated(top%nb_list) ) deallocate(top%nb_list)

    ! calculate size of new NB list
    top%nb_size = (top%natoms - top%probe_size)*top%probe_size
    allocate( top%nb_list(top%nb_size), stat = alloc_stat )

    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate topology arrays in ffdev_topology_switch_to_probe_mode!')
    end if

    ip = 1
    do i=1,top%natoms - top%probe_size
        do j=top%natoms - top%probe_size+1,top%natoms
            top%nb_list(ip)%ai    = i
            top%nb_list(ip)%aj    = j
            top%nb_list(ip)%dt    = 0
            top%nb_list(ip)%nbt   = ffdev_topology_find_nbtype_by_aindex(top,i,j)
            ip = ip + 1
        end do
    end do

    ! regenerate nbtii and nbtjj indexes
    call ffdev_topology_init_nbij(top)

! erase all NB parameters except of probe and probe/probed structure
    if( unique_probe_types ) then
        if( NBParamsMode .eq. NB_PARAMS_MODE_NORMAL ) then
            do i=1,top%nnb_types
                if( top%atom_types(top%nb_types(i)%ti)%probe .or. top%atom_types(top%nb_types(i)%tj)%probe ) cycle
                top%nb_types(i)%eps = 0.0d0
                top%nb_types(i)%r0  = 0.0d0
                top%nb_types(i)%pa  = 0.0d0
                top%nb_types(i)%pb  = 0.0d0
                top%nb_types(i)%rc  = 0.0d0
            end do
        end if
    end if

end subroutine ffdev_topology_switch_to_probe_mode

! ==============================================================================
! subroutine ffdev_topology_gen_sapt_list_for_refs
! ==============================================================================

subroutine ffdev_topology_gen_sapt_list_for_refs(top,nrefs,natomsrefs)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: nrefs
    integer         :: natomsrefs(:)
    ! --------------------------------------------
    integer         :: m,i,j,ip,alloc_status,istart,iend
    ! --------------------------------------------------------------------------

    if( top%probe_size .gt. 0 ) return        !
    if( top%sapt_size .gt. 0 ) return

    ! calculate size of the list
    istart = 1
    ip = 0 ! start from zero for sapt_size
    do m=1,nrefs
        iend = natomsrefs(m)
        do i=istart,istart+iend-1
            do j=istart+iend,top%natoms
                ip = ip + 1
            end do
        end do
        istart = istart + iend
    end do

    top%sapt_size = ip
    if( top%sapt_size .le. 0 ) return

    ! allocate list
    allocate(top%sapt_list(top%sapt_size), stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for MMSAPT in ffdev_topology_gen_sapt_list_for_refs!')
    end if

    istart = 1
    ip = 1 ! start from one for array index
    do m=1,nrefs
        iend = natomsrefs(m)
        do i=istart,istart+iend-1
            do j=istart+iend,top%natoms
                top%sapt_list(ip)%ai    = i
                top%sapt_list(ip)%aj    = j
                top%sapt_list(ip)%dt    = 0
                top%sapt_list(ip)%nbt   = ffdev_topology_find_nbtype_by_aindex(top,i,j)
                top%sapt_list(ip)%nbtii = ffdev_topology_find_nbtype_by_aindex(top,i,i)
                top%sapt_list(ip)%nbtjj = ffdev_topology_find_nbtype_by_aindex(top,j,j)
                ip = ip + 1
            end do
        end do
        istart = istart + iend
    end do

end subroutine ffdev_topology_gen_sapt_list_for_refs

! ==============================================================================
! subroutine ffdev_topology_gen_sapt_list_for_probes
! ==============================================================================

subroutine ffdev_topology_gen_sapt_list_for_probes(top)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i,j,ip,alloc_status
    ! --------------------------------------------------------------------------

    if( top%probe_size .le. 0 ) return        !
    if( top%sapt_size .gt. 0 ) return

    ! calculate size of SAPT list
    top%sapt_size = (top%natoms - top%probe_size)*top%probe_size

    allocate(top%sapt_list(top%sapt_size), stat = alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for MMSAPT in ffdev_topology_gen_sapt_list_for_probes!')
    end if

    ip = 1
    do i=1,top%natoms - top%probe_size
        do j=top%natoms - top%probe_size+1,top%natoms
            top%sapt_list(ip)%ai    = i
            top%sapt_list(ip)%aj    = j
            top%sapt_list(ip)%dt    = 0
            top%sapt_list(ip)%nbt   = ffdev_topology_find_nbtype_by_aindex(top,i,j)
            top%sapt_list(ip)%nbtii = ffdev_topology_find_nbtype_by_aindex(top,i,i)
            top%sapt_list(ip)%nbtjj = ffdev_topology_find_nbtype_by_aindex(top,j,j)
            ip = ip + 1
        end do
    end do

end subroutine ffdev_topology_gen_sapt_list_for_probes

! ==============================================================================
! subroutine ffdev_topology_damptt_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_damptt_mode_to_string(nb_mode)

    use ffdev_utils

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(DAMP_TT_COUPLED)
            ffdev_topology_damptt_mode_to_string = 'COUPLED - TB coupled by damp_pb to PB'
        case(DAMP_TT_FREEOPT)
            ffdev_topology_damptt_mode_to_string = 'FREEOPT - Use optimized TB per type'
        case(DAMP_TT_CONST)
            ffdev_topology_damptt_mode_to_string = 'CONST - constant'
        case(DAMP_TT_BFAC)
            ffdev_topology_damptt_mode_to_string = 'BFAC - Atomic b-factors'
        case(DAMP_TT_BFAC_XDM)
            ffdev_topology_damptt_mode_to_string = 'BFAC-XDM - Atomic b-factors + XDM volumes'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_damptt_mode_to_string!')
    end select

end function ffdev_topology_damptt_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_damptt_mode_from_string
! ==============================================================================

integer function ffdev_topology_damptt_mode_from_string(string)

    use ffdev_utils
    use ffdev_xdm_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('COUPLED')
            ffdev_topology_damptt_mode_from_string = DAMP_TT_COUPLED
        case('FREEOPT')
            ffdev_topology_damptt_mode_from_string = DAMP_TT_FREEOPT
        case('CONST')
            ffdev_topology_damptt_mode_from_string = DAMP_TT_CONST
        case('BFAC')
            ffdev_topology_damptt_mode_from_string = DAMP_TT_BFAC
        case('BFAC-XDM')
            ffdev_topology_damptt_mode_from_string = DAMP_TT_BFAC_XDM
            if( .not. xdm_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1,'XDM data are required for BFACXDM!')
            end if
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_damptt_mode_from_string!')
    end select

end function ffdev_topology_damptt_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_dampbj_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_dampbj_mode_to_string(nb_mode)

    use ffdev_utils

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(DAMP_BJ_DRC)
            ffdev_topology_dampbj_mode_to_string = 'DRC - Use Rc from dispersion data'
        case(DAMP_BJ_ORC)
            ffdev_topology_dampbj_mode_to_string = 'ORC - Use optimized Rc per type'
        case(DAMP_BJ_SRC)
            ffdev_topology_dampbj_mode_to_string = 'SRC - Use scaled optimized Rc per type'
        case(DAMP_BJ_DO)
            ffdev_topology_dampbj_mode_to_string = 'DO - Derived from electron density overlaps'
        case(DAMP_BJ_CONST)
            ffdev_topology_dampbj_mode_to_string = 'CONST - Constant Rc'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_dampbj_mode_to_string!')
    end select

end function ffdev_topology_dampbj_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_dampbj_mode_from_string
! ==============================================================================

integer function ffdev_topology_dampbj_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('DRC')
            ffdev_topology_dampbj_mode_from_string = DAMP_BJ_DRC
        case('ORC')
            ffdev_topology_dampbj_mode_from_string = DAMP_BJ_ORC
        case('SRC')
            ffdev_topology_dampbj_mode_from_string = DAMP_BJ_SRC
        case('DO')
            ffdev_topology_dampbj_mode_from_string = DAMP_BJ_DO
        case('CONST')
            ffdev_topology_dampbj_mode_from_string = DAMP_BJ_CONST
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_dampbj_mode_from_string!')
    end select

end function ffdev_topology_dampbj_mode_from_string

! ==============================================================================
! subroutine ffdev_topology_nb_mode_to_string
! ==============================================================================

character(80) function ffdev_topology_nb_mode_to_string(nb_mode)

    use ffdev_utils

    implicit none
    integer  :: nb_mode
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(NB_VDW_LJ)
            ffdev_topology_nb_mode_to_string = 'LJ - 12-6 Lennard-Jones potential'
        case(NB_VDW_12_DISPBJ)
            ffdev_topology_nb_mode_to_string = '12-DISPBJ - 12/Disp with BJ damping'
        case(NB_VDW_EXP_DISPBJ)
            ffdev_topology_nb_mode_to_string = 'EXP-DISPBJ - Exp/Disp with BJ damping'
        case(NB_VDW_EXP_DISPTT)
            ffdev_topology_nb_mode_to_string = 'EXP-DISPTT - Exp/Disp with TT damping'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_nb_mode_to_string!')
    end select

end function ffdev_topology_nb_mode_to_string

! ==============================================================================
! subroutine ffdev_topology_nb_mode_from_string
! ==============================================================================

integer function ffdev_topology_nb_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('LJ')
            ffdev_topology_nb_mode_from_string = NB_VDW_LJ
        case('12-DISPBJ')
            ffdev_topology_nb_mode_from_string = NB_VDW_12_DISPBJ
        case('EXP-DISPBJ')
            ffdev_topology_nb_mode_from_string = NB_VDW_EXP_DISPBJ
        case('EXP-DISPTT')
            ffdev_topology_nb_mode_from_string = NB_VDW_EXP_DISPTT
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_nb_mode_from_string!')
    end select

end function ffdev_topology_nb_mode_from_string

! FIXME

! ==============================================================================
! subroutine ffdev_topology_comb_rules_to_string
! ==============================================================================

character(80) function ffdev_topology_comb_rules_to_string(comb_rules)

    use ffdev_utils

    implicit none
    integer  :: comb_rules
    ! --------------------------------------------------------------------------

    select case(comb_rules)
        case(COMB_RULE_NONE)
            ffdev_topology_comb_rules_to_string = 'none (input data)'
        ! LJ potential
        case(COMB_RULE_LB)
            ffdev_topology_comb_rules_to_string = 'LB (Lorentz-Berthelot)'
        case(COMB_RULE_WH)
            ffdev_topology_comb_rules_to_string = 'WH (Waldman-Hagler)'
        case(COMB_RULE_KG)
            ffdev_topology_comb_rules_to_string = 'KG (Kong)'
        case(COMB_RULE_FB)
            ffdev_topology_comb_rules_to_string = 'FB (Fender-Halsey-Berthelot)'

        case(COMB_RULE_AM)
            ffdev_topology_comb_rules_to_string = 'AM (Arithmetic mean)'

        case(COMB_RULE_GS)
            ffdev_topology_comb_rules_to_string = 'GS (Gilbert-Smith)'
        case(COMB_RULE_BA)
            ffdev_topology_comb_rules_to_string = 'BA (Bohm-Ahlrichs)'
        case(COMB_RULE_VS)
            ffdev_topology_comb_rules_to_string = 'VS (Vleet-Schmidt)'

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_comb_rules_to_string!')
    end select

end function ffdev_topology_comb_rules_to_string

! ==============================================================================
! function ffdev_topology_comb_rules_from_string
! ==============================================================================

integer function ffdev_topology_comb_rules_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('IN')
            ffdev_topology_comb_rules_from_string = COMB_RULE_NONE

        case('LB')
            ffdev_topology_comb_rules_from_string = COMB_RULE_LB
        case('WH')
            ffdev_topology_comb_rules_from_string = COMB_RULE_WH
        case('KG')
            ffdev_topology_comb_rules_from_string = COMB_RULE_KG
        case('FB')
            ffdev_topology_comb_rules_from_string = COMB_RULE_FB

        case('AM')
            ffdev_topology_comb_rules_from_string = COMB_RULE_AM

        case('GS')
            ffdev_topology_comb_rules_from_string = COMB_RULE_GS
        case('BA')
            ffdev_topology_comb_rules_from_string = COMB_RULE_BA
        case('VS')
            ffdev_topology_comb_rules_from_string = COMB_RULE_VS

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_topology_comb_rules_from_string!')
    end select

end function ffdev_topology_comb_rules_from_string

! ==============================================================================
! subroutine ffdev_topology_apply_NB_comb_rules
! ==============================================================================

subroutine ffdev_topology_apply_NB_comb_rules(top,comb_rules)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: comb_rules
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(NB_VDW_LJ)
            call ffdev_topology_apply_NB_comb_rules_LJ(top,comb_rules)
        case(NB_VDW_12_DISPBJ)
            call ffdev_topology_apply_NB_comb_rules_12BJ(top,comb_rules)
        case(NB_VDW_EXP_DISPBJ)
            call ffdev_topology_apply_NB_comb_rules_EXPBJ(top,comb_rules)
        case(NB_VDW_EXP_DISPTT)
            call ffdev_topology_apply_NB_comb_rules_EXPTT(top,comb_rules)
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_topology_apply_NB_comb_rules!')
    end select

    top%nb_params_update = .true.

end subroutine ffdev_topology_apply_NB_comb_rules

! ==============================================================================
! subroutine ffdev_topology_apply_NB_comb_rules_LJ
! ==============================================================================

subroutine ffdev_topology_apply_NB_comb_rules_LJ(top,comb_rules)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: comb_rules
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: epsii,r0ii,epsjj,r0jj,epsij,r0ij,k,l
    ! --------------------------------------------------------------------------

    ! apply combining rules
    do i=1,top%nnb_types
        if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) then

            ! get type parameters
            nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
            nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)

            r0ii  = top%nb_types(nbii)%r0
            epsii = top%nb_types(nbii)%eps

            r0jj  = top%nb_types(nbjj)%r0
            epsjj = top%nb_types(nbjj)%eps

            select case(comb_rules)
                case(COMB_RULE_LB)
                    r0ij = (r0ii+r0jj)*0.5d0
                    epsij = sqrt(epsii*epsjj)
                case(COMB_RULE_WH)
                    r0ij = ((r0ii**6 + r0jj**6)*0.5d0)**(1.0d0/6.0d0)
                    epsij = sqrt( epsii*r0ii**6 * epsjj*r0jj**6 )/r0ij**6
                case(COMB_RULE_KG)
                    k = sqrt(epsii*r0ii**6 * epsjj*r0jj**6)
                    l = ( ( (epsii*r0ii**12)**(1.0d0/13.0d0) + (epsjj*r0jj**12)**(1.0d0/13.0d0) )*0.5d0 )**13
                    r0ij = (l/k)**(1.0d0/6.0d0)
                    epsij = k / (r0ij**6)
                case(COMB_RULE_FB)
                    r0ij = (r0ii+r0jj)*0.5d0
                    epsij = 2.0d0*epsii*epsjj/(epsii+epsjj)
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_apply_NB_comb_rules_LJ!')
            end select

            top%nb_types(i)%r0 = r0ij
            top%nb_types(i)%eps = epsij

        end if
    end do

end subroutine ffdev_topology_apply_NB_comb_rules_LJ

! ==============================================================================
! subroutine ffdev_topology_apply_NB_comb_rules_12
! ==============================================================================

subroutine ffdev_topology_apply_NB_comb_rules_12BJ(top,comb_rules)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: comb_rules
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: paii,paij,pajj
    real(DEVDP)     :: rcii,rcij,rcjj
    ! --------------------------------------------------------------------------

    ! apply combining rules
    do i=1,top%nnb_types
        if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) then

            ! get type parameters
            nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
            nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)

            paii = top%nb_types(nbii)%pa
            rcii = top%nb_types(nbii)%rc

            pajj = top%nb_types(nbjj)%pa
            rcjj = top%nb_types(nbjj)%rc

            select case(comb_rules)

                case(COMB_RULE_AM)
                    paij = 0.5d0 * (paii+pajj)
                    rcij = 0.5d0 * (rcii+rcjj)

                case(COMB_RULE_WH)
                    paij = 0.5d0 * (paii+pajj)
                    rcij = (0.5d0 * (rcii**6+rcjj**6))**(1.0d0/6.0d0)

                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_apply_NB_comb_rules_12BJ!')
            end select

            top%nb_types(i)%pa = paij
            top%nb_types(i)%rc = rcij

        end if
    end do

end subroutine ffdev_topology_apply_NB_comb_rules_12BJ

! ==============================================================================
! subroutine ffdev_topology_apply_NB_comb_rules_EXPBJ
! ==============================================================================

subroutine ffdev_topology_apply_NB_comb_rules_EXPBJ(top,comb_rules)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: comb_rules
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: paii,paij,pajj
    real(DEVDP)     :: pbii,pbij,pbjj
    real(DEVDP)     :: rcii,rcij,rcjj
    real(DEVDP)     :: eaii,eaij,eajj
    ! --------------------------------------------------------------------------

    ! apply combining rules
    do i=1,top%nnb_types
        if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) then

            ! get type parameters
            nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
            nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)

            paii = top%nb_types(nbii)%pa
            pbii = top%nb_types(nbii)%pb
            rcii = top%nb_types(nbii)%rc

            pajj = top%nb_types(nbjj)%pa
            pbjj = top%nb_types(nbjj)%pb
            rcjj = top%nb_types(nbjj)%rc

            select case(comb_rules)

                case(COMB_RULE_AM)
                    paij = 0.5d0 * (paii+pajj)
                    pbij = 0.5d0 * (pbii+pbjj)
                    rcij = 0.5d0 * (rcii+rcjj)

                case(COMB_RULE_VS)
                    paij = 0.5d0 * (paii+pajj)
                    pbij = sqrt(pbii*pbjj)
                    rcij = 0.5d0 * (rcii+rcjj)

                case(COMB_RULE_BA)
                    if( pbii+pbjj .gt. 0 ) then
                        pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
                    else
                        pbij = 0.5d0 * (pbii+pbjj)
                    end if

                    ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
                    eaii = exp(paii)
                    eajj = exp(pajj)
                    eaij = ( (eaii)**(1.0d0/pbii) * (eajj)**(1.0d0/pbjj) )**(pbij/2.0d0)
                    paij = log(eaij)

                    rcij = 0.5d0 * (rcii+rcjj)

                case(COMB_RULE_GS)
                    if( pbii+pbjj .gt. 0 ) then
                        pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
                    else
                        pbij = 0.5d0 * (pbii+pbjj)
                    end if

                    ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
                    eaii = exp(paii)
                    eajj = exp(pajj)
                    eaij = ( ( (eaii*pbii)**(1.0d0/pbii) * (eajj*pbjj)**(1.0d0/pbjj) )**(pbij/2.0d0) ) / pbij
                    paij = log(eaij)

                    rcij = 0.5d0 * (rcii+rcjj)

                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_apply_NB_comb_rules_EXPBJ!')
            end select

            top%nb_types(i)%pa = paij
            top%nb_types(i)%pb = pbij
            top%nb_types(i)%rc = rcij

        end if
    end do

end subroutine ffdev_topology_apply_NB_comb_rules_EXPBJ

! ==============================================================================
! subroutine ffdev_topology_apply_NB_comb_rules_EXP
! ==============================================================================

subroutine ffdev_topology_apply_NB_comb_rules_EXPTT(top,comb_rules)

    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: comb_rules
    ! --------------------------------------------
    integer         :: i,nbii,nbjj
    real(DEVDP)     :: paii,paij,pajj
    real(DEVDP)     :: pbii,pbij,pbjj
    real(DEVDP)     :: tbii,tbij,tbjj
    real(DEVDP)     :: eaii,eaij,eajj
    ! --------------------------------------------------------------------------

    ! apply combining rules
    do i=1,top%nnb_types
        if( top%nb_types(i)%ti .ne. top%nb_types(i)%tj ) then

            ! get type parameters
            nbii = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%ti,top%nb_types(i)%ti)
            nbjj = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(i)%tj,top%nb_types(i)%tj)

            paii = top%nb_types(nbii)%pa
            pbii = top%nb_types(nbii)%pb
            tbii = top%nb_types(nbii)%tb

            pajj = top%nb_types(nbjj)%pa
            pbjj = top%nb_types(nbjj)%pb
            tbjj = top%nb_types(nbjj)%tb

            select case(comb_rules)

                case(COMB_RULE_AM)
                    paij = 0.5d0 * (paii+pajj)
                    pbij = 0.5d0 * (pbii+pbjj)
                    tbij = 0.5d0 * (tbii+tbjj)

                case(COMB_RULE_GS)
                    if( pbii+pbjj .gt. 0 ) then
                        pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
                    else
                        pbij = 0.5d0 * (pbii+pbjj)
                    end if

                    ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
                    eaii = exp(paii)
                    eajj = exp(pajj)
                    eaij = ( ( (eaii*pbii)**(1.0d0/pbii) * (eajj*pbjj)**(1.0d0/pbjj) )**(pbij/2.0d0) ) / pbij
                    paij = log(eaij)

                    if( tbii+tbjj .gt. 0 ) then
                        tbij = 2.0d0 * tbii*tbjj/(tbii+tbjj)
                    else
                        tbij = 0.5d0 * (tbii+tbjj)
                    end if

                case(COMB_RULE_BA)
                    if( pbii+pbjj .gt. 0 ) then
                        pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
                    else
                        pbij = 0.5d0 * (pbii+pbjj)
                    end if

                    ! Z. Phys. D - Atoms, Molecules and Clusters 1, 91-101 (1986)
                    eaii = exp(paii)
                    eajj = exp(pajj)
                    eaij = ( (eaii)**(1.0d0/pbii) * (eajj)**(1.0d0/pbjj) )**(pbij/2.0d0)
                    paij = log(eaij)

                    if( tbii+tbjj .gt. 0 ) then
                        tbij = 2.0d0 * tbii*tbjj/(tbii+tbjj)
                    else
                        tbij = 0.5d0 * (tbii+tbjj)
                    end if

                case(COMB_RULE_VS)
                    ! DOI:10.1021/acs.jctc.6b00209J. Chem. Theory Comput.2016, 12, 3851−3870
                    paij = 0.5d0 * (paii + pajj)  ! paij is exponential of Aij, etc. ...
                    pbij = sqrt(pbii*pbjj)
                    tbij = sqrt(tbii*tbjj)

                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_topology_apply_NB_comb_rules_EXPTT!')
            end select

            top%nb_types(i)%pa = paij
            top%nb_types(i)%pb = pbij
            top%nb_types(i)%tb = tbij

        end if
    end do

end subroutine ffdev_topology_apply_NB_comb_rules_EXPTT

!===============================================================================
! subroutine ffdev_topology_apply_NB_comb_rules_PB
!===============================================================================

subroutine ffdev_topology_apply_NB_comb_rules_PB(pbii,pbjj,pbij)

    use ffdev_utils

    implicit none
    real(DEVDP)     :: pbii,pbjj,pbij
    ! --------------------------------------------------------------------------

    !write(*,*) pbii,pbjj
    select case(nb_comb_rules)
        case(COMB_RULE_AM)
            pbij = 0.5d0 * (pbii+pbjj)

        case(COMB_RULE_GS)
            if( pbii+pbjj .gt. 0 ) then
                pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
            else
                pbij = 0.5d0 * (pbii+pbjj)
            end if

        case(COMB_RULE_BA)
            if( pbii+pbjj .gt. 0 ) then
                pbij = 2.0d0 * pbii*pbjj/(pbii+pbjj)
            else
                pbij = 0.5d0 * (pbii+pbjj)
            end if

        case(COMB_RULE_VS)
            pbij = sqrt(pbii*pbjj)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented nb_comb_rules in ffdev_topology_update_nbpair_prms!')
    end select

end subroutine ffdev_topology_apply_NB_comb_rules_PB

!===============================================================================
! subroutine ffdev_topology_update_nb_params
!===============================================================================

subroutine ffdev_topology_update_nb_params(top)

    use ffdev_topology_dat
    use ffdev_utils
    use ffdev_disp_dat
    use ffdev_xdm_dat

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: ip
    ! --------------------------------------------------------------------------

    ! write(*,*) 'NB updated'

    select case(nb_mode)
        case(NB_VDW_LJ)
            ! nothing to do
        case(NB_VDW_12_DISPBJ,NB_VDW_EXP_DISPBJ)
            if( .not. disp_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1, &
                     'DISP not loaded for NB_VDW_12_DISPBJ,NB_VDW_EXP_DISPBJ in ffdev_topology_update_nb_params!')
            end if
        case(NB_VDW_EXP_DISPTT)
            if( .not. disp_data_loaded ) then
                call ffdev_utils_exit(DEV_ERR,1, &
                     'DISP not loaded for NB_VDW_EXP_DISPTT in ffdev_topology_update_nb_params!')
            end if
            select case(damptt_mode)
                case(DAMP_TT_BFAC_XDM)
                    if( .not. xdm_data_loaded ) then
                        call ffdev_utils_exit(DEV_ERR,1, &
                             'XDM not loaded for DAMP_TT_BFAC_XDM in ffdev_topology_update_nb_params!')
                     end if
                case default
                    ! nothing to do
            end select
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode in ffdev_topology_update_nb_params!')
    end select

    ! regular NB list
    do ip=1,top%nb_size
        call ffdev_topology_update_nbpair_prms(top,top%nb_list(ip))
    end do

    ! SAPT
    do ip=1,top%sapt_size
        call ffdev_topology_update_nbpair_prms(top,top%sapt_list(ip))
    end do

    ! turn off the update
    top%nb_params_update = .false.

end subroutine ffdev_topology_update_nb_params

!===============================================================================
! subroutine ffdev_topology_update_nb_params
!===============================================================================

subroutine ffdev_topology_update_nbpair_prms(top,nbpair)

    use ffdev_topology_dat
    use ffdev_utils
    use ffdev_xdm_dat
    use ffdev_disp_dat
    use ffdev_densoverlap

    implicit none
    type(TOPOLOGY)  :: top
    type(NB_PAIR)   :: nbpair
    ! --------------------------------------------
    integer         :: i,j,nbt,agti,agtj,nbtii,nbtjj,ti,tj
    real(DEVDP)     :: rc,tbii,tbjj,tbij,pbii,pbjj,pbij
    ! --------------------------------------------------------------------------

    nbpair%crgij = 0.0d0
    nbpair%pa    = 0.0d0
    nbpair%pb    = 0.0d0
    nbpair%c6    = 0.0d0
    nbpair%c8    = 0.0d0
    nbpair%c10   = 0.0d0
    nbpair%rc6   = 0.0d0
    nbpair%rc8   = 0.0d0
    nbpair%rc10  = 0.0d0
    nbpair%tb    = 0.0d0

    i       = nbpair%ai     ! not need to be defined
    j       = nbpair%aj     ! not need to be defined
    nbt     = nbpair%nbt
    nbtii   = nbpair%nbtii
    nbtjj   = nbpair%nbtjj
    ti      = top%nb_types(nbtii)%ti
    tj      = top%nb_types(nbtjj)%tj
    agti    = top%atom_types(ti)%glbtypeid
    agtj    = top%atom_types(tj)%glbtypeid

    ! electrostatics
    ! FIXME - better value
    if( (i .ne. 0 ) .and. (j .ne. 0) ) then
        nbpair%crgij = 332.05221729d0 * top%atoms(i)%charge * top%atoms(j)%charge * ele_qscale * ele_qscale
    end if

    select case(nb_mode)
        case(NB_VDW_LJ)
            ! LJ parameters
            nbpair%pa =         top%nb_types(nbt)%eps * top%nb_types(nbt)%r0**12
            nbpair%c6 = 2.0d0 * top%nb_types(nbt)%eps * top%nb_types(nbt)%r0**6

        case(NB_VDW_12_DISPBJ)
            ! ^12 repulsion
            nbpair%pa  = exp(top%nb_types(nbt)%pa)

            ! dispersion coefficients
            nbpair%c6  = disp_s6  * disp_pairs(agti,agtj)%c6
            nbpair%c8  = disp_s8  * disp_pairs(agti,agtj)%c8
            nbpair%c10 = disp_s10 * disp_pairs(agti,agtj)%c10

            ! RC for BJ damping
            rc  = 0.0d0
            select case(dampbj_mode)
                case(DAMP_BJ_DRC)
                    rc = damp_fa * disp_pairs(agti,agtj)%rc + damp_fb
                case(DAMP_BJ_ORC)
                    rc = top%nb_types(nbt)%rc
                case(DAMP_BJ_SRC)
                    rc = damp_fa * top%nb_types(nbt)%rc + damp_fb
                case(DAMP_BJ_DO)
                    rc =  ffdev_densoverlap_rc(agti,damp_fa) + ffdev_densoverlap_rc(agtj,damp_fa)
                case default
                    if( .not. disp_data_loaded ) then
                        call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_topology_update_nbpair_prms!')
                    end if
            end select

            nbpair%rc6  = rc**6
            nbpair%rc8  = rc**8
            nbpair%rc10 = rc**10

        case(NB_VDW_EXP_DISPBJ)
            ! Pauli repulsion
            if( pb_from_densoverlap ) then
                pbii = damp_pb * ffdev_densoverlap_b(agti)
                pbjj = damp_pb * ffdev_densoverlap_b(agtj)

                call ffdev_topology_apply_NB_comb_rules_PB(pbii,pbjj,pbij)
                nbpair%pb  = pbij
            else
                nbpair%pb  = damp_pb * top%nb_types(nbt)%pb
            end if

            if( couple_pa_pb_forA ) then
                nbpair%pa  = exp(top%nb_types(nbt)%pa * nbpair%pb)
            else
                nbpair%pa  = exp(top%nb_types(nbt)%pa)
            end if

            ! dispersion coefficients
            nbpair%c6  = disp_s6  * disp_pairs(agti,agtj)%c6
            nbpair%c8  = disp_s8  * disp_pairs(agti,agtj)%c8
            nbpair%c10 = disp_s10 * disp_pairs(agti,agtj)%c10

            ! RC for BJ damping
            rc  = 0.0d0
            select case(dampbj_mode)
                case(DAMP_BJ_DRC)
                    rc = damp_fa * disp_pairs(agti,agtj)%rc + damp_fb
                case(DAMP_BJ_ORC)
                    rc = top%nb_types(nbt)%rc
                case(DAMP_BJ_SRC)
                    rc = damp_fa * top%nb_types(nbt)%rc + damp_fb
                case(DAMP_BJ_CONST)
                    rc = damp_fa
                case(DAMP_BJ_DO)
                    rc =  ffdev_densoverlap_rc(agti,damp_fa) + ffdev_densoverlap_rc(agtj,damp_fa)
                case default
                    if( .not. disp_data_loaded ) then
                        call ffdev_utils_exit(DEV_ERR,1,'RC mode not implemented in ffdev_topology_update_nbpair_prms!')
                    end if
            end select

            nbpair%rc6  = rc**6
            nbpair%rc8  = rc**8
            nbpair%rc10 = rc**10

        case(NB_VDW_EXP_DISPTT)

            ! Pauli repulsion
            if( pb_from_densoverlap ) then
                pbii = damp_pb * ffdev_densoverlap_b(agti)
                pbjj = damp_pb * ffdev_densoverlap_b(agtj)

                call ffdev_topology_apply_NB_comb_rules_PB(pbii,pbjj,pbij)
                nbpair%pb  = pbij
            else
                nbpair%pb  = damp_pb * top%nb_types(nbt)%pb
            end if

            if( couple_pa_pb_forA ) then
                nbpair%pa  = exp(top%nb_types(nbt)%pa * nbpair%pb)
            else
                nbpair%pa  = exp(top%nb_types(nbt)%pa)
            end if

            ! dispersion coefficients
            nbpair%c6  = disp_s6  * disp_pairs(agti,agtj)%c6
            nbpair%c8  = disp_s8  * disp_pairs(agti,agtj)%c8
            nbpair%c10 = disp_s10 * disp_pairs(agti,agtj)%c10

            ! TT damping
            select case(damptt_mode)
                case(DAMP_TT_COUPLED)
                    nbpair%tb  = damp_fa * nbpair%pb
                case(DAMP_TT_FREEOPT)
                    nbpair%tb  = top%nb_types(nbt)%tb
                case(DAMP_TT_CONST)
                    nbpair%tb  = damp_fa
                case(DAMP_TT_BFAC)
                    tbii = damp_fa * ffdev_densoverlap_b(agti)
                    tbjj = damp_fa * ffdev_densoverlap_b(agtj)

                    call ffdev_topology_apply_NB_comb_rules_PB(tbii,tbjj,tbij)
                    nbpair%tb  = tbij

                case(DAMP_TT_BFAC_XDM)
                    tbii = damp_fa * ffdev_densoverlap_b(agti) * (xdm_atoms(agti)%vave / xdm_atoms(agti)%v0ave)**damp_fb
                    tbjj = damp_fa * ffdev_densoverlap_b(agtj) * (xdm_atoms(agtj)%vave / xdm_atoms(agtj)%v0ave)**damp_fb

                    call ffdev_topology_apply_NB_comb_rules_PB(tbii,tbjj,tbij)
                    nbpair%tb  = tbij

                case default
                    if( .not. disp_data_loaded ) then
                        call ffdev_utils_exit(DEV_ERR,1,'TT damp mode not implemented in ffdev_topology_update_nbpair_prms!')
                    end if
            end select

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode in ffdev_topology_update_nbpair_prms!')
    end select

end subroutine ffdev_topology_update_nbpair_prms

! ------------------------------------------------------------------------------

end module ffdev_topology


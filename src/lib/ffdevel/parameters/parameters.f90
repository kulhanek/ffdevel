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

module ffdev_parameters

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_parameters_init
! ==============================================================================

subroutine ffdev_parameters_init()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer     :: maxnparams, i, j, k, alloc_stat, parmid
    ! --------------------------------------------------------------------------
    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'PARAMETERS', ':')

    ! print parameters of individual topologies
    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,10) i
        call ffdev_topology_info_types(sets(i)%top)
    end do

    ! extract unique types
    call ffdev_parameters_gen_unique_types()
    call ffdev_parameters_print_types()

    ! generate parameters ------------------------------------------------------

    ! determine maximum number of parameters
    maxnparams = 0
    do i=1,nsets
        maxnparams = maxnparams + 1     ! energy offset for a set
        ! topology related data
        maxnparams = maxnparams + 2*sets(i)%top%nbond_types     ! bonds
        maxnparams = maxnparams + 2*sets(i)%top%nangle_types    ! angles
        maxnparams = maxnparams + 2*sets(i)%top%ndihedral_types*sets(i)%top%ndihedral_seq_size ! dihedrals
        maxnparams = maxnparams + 2*sets(i)%top%ndihedral_types ! dihedral scee, scnb
        maxnparams = maxnparams + 2*sets(i)%top%nimproper_types ! impropers
        maxnparams = maxnparams + 3*sets(i)%top%nnb_types       ! NB
    end do

    write(DEV_OUT,*)
    write(DEV_OUT,20) maxnparams

    allocate(params(maxnparams), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for parameter extraction!')
    end if
    do i=1,maxnparams
        allocate(params(i)%ids(nsets), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for parameter extraction!')
        end if
    end do

    nparams = 0

    ! energy offset realm ==================
    do i=1,nsets
        nparams = nparams + 1
        params(nparams)%value = sets(i)%offset
        params(nparams)%realm = REALM_EOFFSET
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ids(i) = 1
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end do

    ! bond length realm ====================
    do i=1,nsets
        do j=1,sets(i)%top%nbond_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_BOND_D0)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%bond_types(j)%d0
                params(nparams)%realm = REALM_BOND_D0
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%bond_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%bond_types(j)%tj)
                params(nparams)%tk   = 0
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! bond force realm ====================
    do i=1,nsets
        do j=1,sets(i)%top%nbond_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_BOND_K)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%bond_types(j)%k
                params(nparams)%realm = REALM_BOND_K
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%bond_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%bond_types(j)%tj)
                params(nparams)%tk   = 0
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! angle realm ==========================
    do i=1,nsets
        do j=1,sets(i)%top%nangle_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_ANGLE_A0)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%angle_types(j)%a0
                params(nparams)%realm = REALM_ANGLE_A0
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%angle_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%angle_types(j)%tj)
                params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%angle_types(j)%tk)
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! angle force realm ====================
    do i=1,nsets
        do j=1,sets(i)%top%nangle_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_ANGLE_K)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%angle_types(j)%k
                params(nparams)%realm = REALM_ANGLE_K
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%angle_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%angle_types(j)%tj)
                params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%angle_types(j)%tk)
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! dihedral v realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%ndihedral_types
            if( sets(i)%top%dihedral_types(j)%mode .ne. DIH_COS ) cycle
            do k=1,sets(i)%top%dihedral_types(j)%n
                parmid = find_parameter(sets(i)%top,j,k,REALM_DIH_V)
                if( parmid .eq. 0 ) then    ! new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%dihedral_types(j)%v(k)
                    params(nparams)%realm = REALM_DIH_V
                    params(nparams)%enabled = .false.
                    params(nparams)%identity = 0
                    params(nparams)%pn    = k
                    params(nparams)%ids(:) = 0
                    params(nparams)%ids(i) = j
                    params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%ti)
                    params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tj)
                    params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tk)
                    params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tl)
                else
                    params(parmid)%ids(i) = j ! parameter already exists, update link
                end if
            end do
        end do
    end do

    ! dihedral c realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%ndihedral_types
            if( sets(i)%top%dihedral_types(j)%mode .ne. DIH_GRBF ) cycle
            do k=1,sets(i)%top%dihedral_types(j)%n
                parmid = find_parameter(sets(i)%top,j,k,REALM_DIH_C)
                if( parmid .eq. 0 ) then    ! new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%dihedral_types(j)%c(k)
                    params(nparams)%realm = REALM_DIH_C
                    params(nparams)%enabled = .false.
                    params(nparams)%identity = 0
                    params(nparams)%pn    = k
                    params(nparams)%ids(:) = 0
                    params(nparams)%ids(i) = j
                    params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%ti)
                    params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tj)
                    params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tk)
                    params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tl)
                else
                    params(parmid)%ids(i) = j ! parameter already exists, update link
                end if
            end do
        end do
    end do

    ! dihedral g realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%ndihedral_types
            if( sets(i)%top%dihedral_types(j)%mode .ne. DIH_COS ) cycle
            do k=1,sets(i)%top%dihedral_types(j)%n
                parmid = find_parameter(sets(i)%top,j,k,REALM_DIH_G)
                if( parmid .eq. 0 ) then    ! new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%dihedral_types(j)%g(k)
                    params(nparams)%realm = REALM_DIH_G
                    params(nparams)%enabled = .false.
                    params(nparams)%identity = 0
                    params(nparams)%pn    = k
                    params(nparams)%ids(:) = 0
                    params(nparams)%ids(i) = j
                    params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%ti)
                    params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tj)
                    params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tk)
                    params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tl)
                else
                    params(parmid)%ids(i) = j ! parameter already exists, update link
                end if
            end do
        end do
    end do

    ! dihedral scee realm ==================
    do i=1,nsets
        do j=1,sets(i)%top%ndihedral_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_DIH_SCEE)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = 1.0d0 / sets(i)%top%dihedral_types(j)%inv_scee
                params(nparams)%realm = REALM_DIH_SCEE
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tj)
                params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tk)
                params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tl)
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! dihedral scnb realm ==================
    do i=1,nsets
        do j=1,sets(i)%top%ndihedral_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_DIH_SCNB)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = 1.0d0 / sets(i)%top%dihedral_types(j)%inv_scnb
                params(nparams)%realm = REALM_DIH_SCNB
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tj)
                params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tk)
                params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%dihedral_types(j)%tl)
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! imroper v realm ======================
    do i=1,nsets
        do j=1,sets(i)%top%nimproper_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_IMPR_V)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%improper_types(j)%v
                params(nparams)%realm = REALM_IMPR_V
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%tj)
                params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%tk)
                params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%tl)
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! improper g realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%nimproper_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_IMPR_G)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%improper_types(j)%g
                params(nparams)%realm = REALM_IMPR_G
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%tj)
                params(nparams)%tk   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%tk)
                params(nparams)%tl   = get_common_type_id(sets(i)%top,sets(i)%top%improper_types(j)%tl)
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! vdw A realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_A)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%nb_types(j)%A
                params(nparams)%realm = REALM_VDW_A
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%nb_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%nb_types(j)%tj)
                params(nparams)%tk   = 0
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! vdw B realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_B)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%nb_types(j)%B
                params(nparams)%realm = REALM_VDW_B
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%nb_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%nb_types(j)%tj)
                params(nparams)%tk   = 0
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    ! vdw C realm =====================
    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_C)
            if( parmid .eq. 0 ) then    ! new parameter
                nparams = nparams + 1
                params(nparams)%value = sets(i)%top%nb_types(j)%C
                params(nparams)%realm = REALM_VDW_C
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn    = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j
                params(nparams)%ti   = get_common_type_id(sets(i)%top,sets(i)%top%nb_types(j)%ti)
                params(nparams)%tj   = get_common_type_id(sets(i)%top,sets(i)%top%nb_types(j)%tj)
                params(nparams)%tk   = 0
                params(nparams)%tl   = 0
            else
                params(parmid)%ids(i) = j ! parameter already exists, update link
            end if
        end do
    end do

    write(DEV_OUT,30) nparams

    ! load parameters from file
    if( (len_trim(InpParamFileName) .gt. 0) .and. (trim(InpParamFileName) .ne. '-none-') ) then
        write(DEV_OUT,40) trim(InpParamFileName)
        call ffdev_parameters_load(InpParamFileName)
    else
        write(DEV_OUT,50)
    end if

    call ffdev_parameters_print_parameters()

 10 format('=== [topology] #',I2.2,' =============================================================')
 20 format('Estimated number of parameters (maximum) = ',I6)
 30 format('Final number of parameters               = ',I6)
 40 format('Loading input parameters from file       = ',A)
 50 format('Input parameters taken from topologies ... ')

end subroutine ffdev_parameters_init

! ------------------------------------------------------------------------------

integer function find_parameter(top,id,pn,realm)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: id
    integer         :: pn
    integer         :: realm
    ! --------------------------------------------
    integer         :: ti,tj,tk,tl
    integer         :: i
    ! --------------------------------------------------------------------------

    find_parameter = 0

    ! find global type ids for given parameter realm
    select case(realm)
        case(REALM_BOND_D0,REALM_BOND_K)
            ti = get_common_type_id(top,top%bond_types(id)%ti)
            tj = get_common_type_id(top,top%bond_types(id)%tj)
        case(REALM_ANGLE_A0,REALM_ANGLE_K)
            ti = get_common_type_id(top,top%angle_types(id)%ti)
            tj = get_common_type_id(top,top%angle_types(id)%tj)
            tk = get_common_type_id(top,top%angle_types(id)%tk)
        case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G,REALM_DIH_SCEE,REALM_DIH_SCNB)
            ti = get_common_type_id(top,top%dihedral_types(id)%ti)
            tj = get_common_type_id(top,top%dihedral_types(id)%tj)
            tk = get_common_type_id(top,top%dihedral_types(id)%tk)
            tl = get_common_type_id(top,top%dihedral_types(id)%tl)
        case(REALM_IMPR_V,REALM_IMPR_G)
            ti = get_common_type_id(top,top%improper_types(id)%ti)
            tj = get_common_type_id(top,top%improper_types(id)%tj)
            tk = get_common_type_id(top,top%improper_types(id)%tk)
            tl = get_common_type_id(top,top%improper_types(id)%tl)
        case(REALM_VDW_A,REALM_VDW_B,REALM_VDW_C)
            ti = get_common_type_id(top,top%nb_types(id)%ti)
            tj = get_common_type_id(top,top%nb_types(id)%tj)
    end select

    ! is parameter defined?
    do i=1,nparams
        if( params(i)%realm .ne. realm ) cycle
        if( params(i)%pn .ne. pn ) cycle

        select case(realm)
            case(REALM_BOND_D0,REALM_BOND_K)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj)) .or. &
                    ((params(i)%ti .eq. tj) .and. (params(i)%tj .eq. ti)) ) then
                        find_parameter = i
                end if
            case(REALM_ANGLE_A0,REALM_ANGLE_K)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj) .and. (params(i)%tk .eq. tk)) .or. &
                    ((params(i)%ti .eq. tk) .and. (params(i)%tj .eq. tj) .and. (params(i)%tk .eq. ti)) ) then
                        find_parameter = i
                end if
            case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G,REALM_DIH_SCEE,REALM_DIH_SCNB)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. tl)) .or. &
                    ((params(i)%ti .eq. tl) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. ti)) ) then
                        find_parameter = i
                end if
            case(REALM_IMPR_V,REALM_IMPR_G)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. tl)) .or. &
                    ((params(i)%ti .eq. tl) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. ti)) ) then
                        find_parameter = i
                end if
            case(REALM_VDW_A,REALM_VDW_B,REALM_VDW_C)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj)) .or. &
                    ((params(i)%ti .eq. tj) .and. (params(i)%tj .eq. ti)) ) then
                        find_parameter = i
                end if
        end select
    end do

end function find_parameter

! ------------------------------------------------------------------------------

integer function get_common_type_id(top,typeid)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: typeid
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    get_common_type_id = 0

    do i=1,ntypes
        if( top%atom_types(typeid)%name .eq. types(i)%name ) then
            get_common_type_id = i
            return
        end if
    end do

end function get_common_type_id

! ==============================================================================
! subroutine ffdev_parameters_gen_unique_types
! ==============================================================================

subroutine ffdev_parameters_gen_unique_types()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer                         :: maxnparams, i, j, k, alloc_stat
    character(len=4),allocatable    :: ltypes(:)
    character(len=4)                :: tmp
    logical                         :: changed
    ! --------------------------------------------------------------------------

    ! generate list --------------------
    maxnparams = 0
    do i=1,nsets
        maxnparams = maxnparams + sets(i)%top%natom_types
    end do

    allocate(ltypes(maxnparams), stat = alloc_stat)
    if(alloc_stat .ne. 0) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for all types!')
    end if

    k = 1
    do i=1,nsets
        do j=1,sets(i)%top%natom_types
            ltypes(k) = sets(i)%top%atom_types(j)%name
            k = k + 1
        end do
    end do

    ! sort types - brutal force
    changed = .true.
    do while( changed )
        changed = .false.
        do i=1,maxnparams-1
            if( ltypes(i) .lt. ltypes(i+1) ) then
                changed = .true.
                tmp = ltypes(i+1)
                ltypes(i+1) = ltypes(i)
                ltypes(i) = tmp
            end if
        end do
    end do

    ! count unique types
    k = 0
    tmp = ''
    do i=1,maxnparams
        if( tmp .ne. ltypes(i) ) then
            k = k + 1
            tmp = ltypes(i)
        end if
    end do
    ntypes = k

    ! allocate array for types
    allocate(types(ntypes), stat = alloc_stat)
    if(alloc_stat .ne. 0) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for common types!')
    end if

    ! copy types
    k = 0
    tmp = ''
    do i=1,maxnparams
        if( tmp .ne. ltypes(i) ) then
            k = k + 1
            types(k)%name = ltypes(i)
            tmp = ltypes(i)
        end if
    end do

    ! init set occurence
    do i=1,ntypes
        allocate(types(i)%ids(nsets), stat = alloc_stat)
        if(alloc_stat .ne. 0) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for type ids!')
        end if
        types(i)%ids(:) = 0
        do j=1,nsets
            do k=1,sets(j)%top%natom_types
                if( types(i)%name .eq. sets(j)%top%atom_types(k)%name ) then
                    types(i)%ids(j) = k
                    types(i)%z = sets(j)%top%atom_types(k)%z
                    types(i)%mass = sets(j)%top%atom_types(k)%mass
                    exit
                end if
            end do
        end do
    end do

end subroutine ffdev_parameters_gen_unique_types

! ==============================================================================
! subroutine ffdev_parameters_load
! ==============================================================================

subroutine ffdev_parameters_load(name)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    character(*)    :: name
    ! --------------------------------------------
    integer         :: i, count, j, ri, rcount, rrealm, rpn, nlp
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_PRMS,name,'O')

    read(DEV_PRMS,10) nlp
    if( nlp .ne. nparams ) then
         call ffdev_utils_exit(DEV_OUT,1,'Data mismatch between parameter file and gathered parameters from topologies!')
    end if

    do i=1,nparams
        count = 0
        do j=1,nsets
            if( params(i)%ids(j) .gt. 0 ) count = count + 1
        end do
        read(DEV_PRMS,20) ri, rcount, rrealm, rpn, params(i)%value
        if( (ri .ne. i) .or. (rcount .ne. count) .or. &
            (rrealm .ne. params(i)%realm) .or.  (rpn .ne. params(i)%pn) ) then
             call ffdev_utils_exit(DEV_OUT,1,'Data mismatch between parameter file and gathered parameters from topologies!')
        end if
    end do

    close(DEV_PRMS)

 10 format('FFDEVEL PARAMATERS',1X,I6)
 20 format(I6,1X,I3,1X,I2,1X,I2,F20.6)

end subroutine ffdev_parameters_load

! ==============================================================================
! subroutine ffdev_parameters_save
! ==============================================================================

subroutine ffdev_parameters_save(name)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    character(*)    :: name
    ! --------------------------------------------
    integer         :: i, count, j
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_PRMS,name,'U')

    write(DEV_PRMS,10) nparams

    do i=1,nparams
        count = 0
        do j=1,nsets
            if( params(i)%ids(j) .gt. 0 ) count = count + 1
        end do
        write(DEV_PRMS,20) i, count, params(i)%realm, params(i)%pn, params(i)%value
    end do

    close(DEV_PRMS)

 10 format('FFDEVEL PARAMATERS',1X,I6)
 20 format(I6,1X,I3,1X,I2,1X,I2,F20.6)

end subroutine ffdev_parameters_save

! ==============================================================================
! subroutine ffdev_parameters_save
! ==============================================================================

subroutine ffdev_parameters_save_amber(name)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    character(*)    :: name
    ! --------------------------------------------
    integer         :: datum(8)
    integer         :: i,j,it,ij,pn,max_pn
    real(DEVDP)     :: v,g, scee, scnb
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_PRMS,name,'U')

    call date_and_time(values=datum)
    write(DEV_PRMS,10) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)

    ! atom types
    write(DEV_PRMS,20) 'MASS'
    do i=1,ntypes
        write(DEV_PRMS,30) types(i)%name,types(i)%mass,0.0d0
    end do
    write(DEV_PRMS,*)

    ! bonds
    it = 0
    write(DEV_PRMS,20) 'BOND'
    do i=1,nparams
        if( params(i)%realm .ne. REALM_BOND_D0 ) cycle
        it = it + 1
        ij = 0
        do j=1,nparams
            if( params(j)%realm .ne. REALM_BOND_K ) cycle
            ij = ij + 1
            if( ij .eq. it ) exit
        end do
        write(DEV_PRMS,40) types(params(i)%ti)%name,types(params(i)%tj)%name, &
                           params(j)%value,params(i)%value
    end do
    write(DEV_PRMS,*)

    ! angles
    it = 0
    write(DEV_PRMS,20) 'ANGL'
    do i=1,nparams
        if( params(i)%realm .ne. REALM_ANGLE_A0 ) cycle
        it = it + 1
        ij = 0
        do j=1,nparams
            if( params(j)%realm .ne. REALM_ANGLE_K ) cycle
            ij = ij + 1
            if( ij .eq. it ) exit
        end do
        write(DEV_PRMS,50) types(params(i)%ti)%name,types(params(i)%tj)%name,types(params(i)%tk)%name, &
                           params(j)%value,params(i)%value*DEV_R2D
    end do
    write(DEV_PRMS,*)

    ! dihedrals
    it = 0
    write(DEV_PRMS,20) 'DIHE'
    do i=1,nparams
        if( params(i)%realm .ne. REALM_DIH_V ) cycle
        v = params(i)%value
        max_pn = 0
        do j=1,nparams
            if( params(j)%realm .ne. REALM_DIH_V ) cycle
            if( (params(i)%ti .eq. params(j)%ti) .and. &
                (params(i)%tj .eq. params(j)%tj) .and. &
                (params(i)%tk .eq. params(j)%tk) .and. &
                (params(i)%tl .eq. params(j)%tl) ) then
                    max_pn  = max_pn + 1
            end if
        end do
        pn = params(i)%pn
        if( pn .lt. max_pn ) pn = -pn
        it = it + 1
        ij = 0
        g = 0.0
        do j=1,nparams
            if( params(j)%realm .ne. REALM_DIH_G ) cycle
            ij = ij + 1
            g = params(j)%value
            if( ij .eq. it ) exit
        end do
        ij = 0
        scee = 1.2
        do j=1,nparams
            if( params(j)%realm .ne. REALM_DIH_SCEE ) cycle
            ij = ij + 1
            scee = params(j)%value
            if( ij .eq. it ) exit
        end do
        ij = 0
        scnb = 2.0
        do j=1,nparams
            if( params(j)%realm .ne. REALM_DIH_SCNB ) cycle
            ij = ij + 1
            scnb = params(j)%value
            if( ij .eq. it ) exit
        end do
        write(DEV_PRMS,60) types(params(i)%ti)%name,types(params(i)%tj)%name, &
                           types(params(i)%tk)%name, types(params(i)%tl)%name, &
                           v, g*DEV_R2D, pn, scee, scnb
    end do
    write(DEV_PRMS,*)

    ! impropers
    it = 0
    write(DEV_PRMS,20) 'IMPR'
    do i=1,nparams
        if( params(i)%realm .ne. REALM_IMPR_V ) cycle
        v = params(i)%value
        pn = 2
        it = it + 1
        ij = 0
        g = 0.0
        do j=1,nparams
            if( params(j)%realm .ne. REALM_IMPR_G ) cycle
            ij = ij + 1
            g = params(j)%value
            if( ij .eq. it ) exit
        end do
        write(DEV_PRMS,70) types(params(i)%ti)%name,types(params(i)%tj)%name, &
                           types(params(i)%tk)%name, types(params(i)%tl)%name, &
                           v, g*DEV_R2D, pn
    end do
    write(DEV_PRMS,*)

    ! hydrogen bonds
    write(DEV_PRMS,20) 'HBON'
    write(DEV_PRMS,*)

    ! nonbonded
! FIXME
!    write(DEV_PRMS,20) 'NONB'
!    do i=1,ntypes
!        write(DEV_PRMS,80) types(i)%name,types(i)%r0,types(i)%eps
!    end do
!    write(DEV_PRMS,*)

    close(DEV_PRMS)

 10 format('# parameters generated by ffdevel package ',i4,'-',i2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2)
 20 format(A)
 30 format(A4,10X,F16.6,1X,F16.6)
 40 format(A2,'-',A2,9X,F16.6,1X,F16.6)
 50 format(A2,'-',A2,'-',A2,6X,F16.6,1X,F16.6)
 60 format(A2,'-',A2,'-',A2,'-',A2,1X,'1',1X,F16.6,1X,F16.6,1X,I4,1X,'SCEE=',F4.2,1X,'SCNB=',F4.2)
 70 format(A2,'-',A2,'-',A2,'-',A2,3X,F16.6,1X,F16.6,1X,I4)
 80 format(A4,10X,F16.6,1X,F16.6)

end subroutine ffdev_parameters_save_amber

! ==============================================================================
! subroutine ffdev_parameters_print_types
! ==============================================================================

subroutine ffdev_parameters_print_types()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer     :: i,j,count
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Atom Type Summary', '=')
    write(DEV_OUT,5) ntypes

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    do i=1,ntypes
        count = 0
        do j=1,nsets
            if( types(i)%ids(j) .gt. 0 ) count = count + 1
        end do
        write(DEV_OUT,30,ADVANCE='NO') i, types(i)%name, count
        do j=1,nsets
            if( types(i)%ids(j) .gt. 0 ) then
                write(DEV_OUT,40,ADVANCE='NO') types(i)%ids(j)
            else
                write(DEV_OUT,50,ADVANCE='NO')
            end if
        end do
        write(DEV_OUT,*)
    end do

  5 format('Number of unique atom types among all target sets = ',I6)

 10 format('# ID Type Counts     Atom Type IDs in Sets')
 20 format('# -- ---- ------     -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --')
 30 format(I4,1X,A4,1X,I6,5X)
 40 format(I2,1X)
 50 format('--',1X)

end subroutine ffdev_parameters_print_types

! ==============================================================================
! subroutine ffdev_parameters_print_parameters
! ==============================================================================

subroutine ffdev_parameters_print_parameters()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer             :: i,j,count,free,act,tot
    character(len=11)   :: tmp   
    real(DEVDP)         :: scaling
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Parameter Summary', ':')

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    do i=1,nparams
        count = 0
        scaling = 1.0d0
        do j=1,nsets
            if( params(i)%ids(j) .gt. 0 ) count = count + 1
        end do
        if( params(i)%identity .gt. 0 ) then
            write(DEV_OUT,30,ADVANCE='NO') i, params(i)%enabled, params(i)%identity
        else
            write(DEV_OUT,31,ADVANCE='NO') i, params(i)%enabled
        end if

        select case(params(i)%realm)
            case(REALM_EOFFSET)
                tmp = 'E_offset'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_BOND_D0)
                tmp = 'bond_r0'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_BOND_K)
                tmp = 'bond_k'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_ANGLE_A0)
                tmp = 'angle_a0'
                scaling = DEV_R2D
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_ANGLE_K)
                tmp = 'angle_k'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_DIH_V)
                tmp = 'dih_v'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_DIH_C)
                tmp = 'dih_c'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_DIH_G)
                tmp = 'dih_gamma'
                scaling = DEV_R2D
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_DIH_SCEE)
                tmp = 'dih_scee'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_DIH_SCNB)
                tmp = 'dih_scnb'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_IMPR_V)
                tmp = 'impr_v'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_IMPR_G)
                tmp = 'impr_gamma'
                scaling = DEV_R2D
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_VDW_A)
                tmp = 'vdw_A'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_VDW_B)
                tmp = 'vdw_B'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case(REALM_VDW_C)
                tmp = 'vdw_C'
                write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdev_parameters_print_parameters!')
        end select

        write(DEV_OUT,35,ADVANCE='NO') params(i)%pn, params(i)%value*scaling, count
        do j=1,nsets
            if( params(i)%ids(j) .gt. 0 ) then
                write(DEV_OUT,40,ADVANCE='NO') params(i)%ids(j)
            else
                write(DEV_OUT,50,ADVANCE='NO')
            end if
        end do
        write(DEV_OUT,*)
    end do

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)

    do j=REALM_FIRST,REALM_LAST
        tot = 0
        free = 0
        act = 0
        do i=1,nparams

        end do
        do i=1,nparams
            if( params(i)%realm .ne. j ) cycle
            tot = tot + 1
            if( params(i)%identity .eq. 0 ) free = free + 1
            if( params(i)%enabled ) act = act + 1
        end do

        select case(j)
            case(REALM_EOFFSET)
                tmp = 'E_offset'
            case(REALM_BOND_D0)
                tmp = 'bond_r0'
            case(REALM_BOND_K)
                tmp = 'bond_k'
            case(REALM_ANGLE_A0)
                tmp = 'angle_a0'
            case(REALM_ANGLE_K)
                tmp = 'angle_k'
            case(REALM_DIH_V)
                tmp = 'dih_v'
            case(REALM_DIH_C)
                tmp = 'dih_c'
            case(REALM_DIH_G)
                tmp = 'dih_gamma'
            case(REALM_DIH_SCEE)
                tmp = 'dih_scee'
            case(REALM_DIH_SCNB)
                tmp = 'dih_scnb'
            case(REALM_IMPR_V)
                tmp = 'impr_v'
            case(REALM_IMPR_G)
                tmp = 'impr_gamma'
            case(REALM_VDW_A)
                tmp = 'vdW_A'
            case(REALM_VDW_B)
                tmp = 'vdW_B'
            case(REALM_VDW_C)
                tmp = 'vdW_C'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Not implemented in ffdev_parameters_print_parameters!')
        end select

        write(DEV_OUT,120) adjustl(tmp), tot, free, tot-free, act

    end do

    free = 0
    act = 0
    do i=1,nparams
        if( params(i)%identity .eq. 0 ) free = free + 1
        if( params(i)%enabled ) act = act + 1
    end do

    write(DEV_OUT,*)
    write(DEV_OUT,200) nparams
    write(DEV_OUT,210) free
    write(DEV_OUT,220) nparams - free
    write(DEV_OUT,230) act


 10 format('# ID ST Iden    Realm    PN   Value      Counts     IDs in Sets')
 20 format('# -- -- ---- ----------- -- ------------ ------     -- -- -- -- -- -- -- -- -- -- --')
 30 format(I4,1X,L2,1X,I4,1X)
 31 format(I4,1X,L2,1X,'----',1X)
 32 format(A11,1X)
 35 format(I2,1X,F12.3,1X,I6,5X)
 40 format(I2,1X)
 50 format('--',1X)

100 format('#    Realm    Total   Free  Ident  Active')
110 format('# ----------- ------ ------ ------ ------')
120 format(2X,A11,1X,I6,1X,I6,1X,I6,1X,I6)

200 format('Number of all parameters    = ',I6)
210 format('Number of free parameters   = ',I6)
220 format('Number of identities        = ',I6)
230 format('Number of active parameters = ',I6)

end subroutine ffdev_parameters_print_parameters

! ==============================================================================
! subroutine ffdev_targetset_disable_all_realms
! ==============================================================================

subroutine ffdev_parameters_disable_all_realms()

    use ffdev_parameters_dat

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    do i=1,nparams
        params(i)%enabled = .false.
    end do

end subroutine ffdev_parameters_disable_all_realms

! ==============================================================================
! subroutine ffdev_parameters_gather
! ==============================================================================

subroutine ffdev_parameters_gather(prms)

    use ffdev_parameters_dat

    implicit none
    real(DEVDP)     :: prms(:)
    ! --------------------------------------------
    integer         :: idx, i
    ! --------------------------------------------------------------------------

    idx = 1
    do i=1,nparams
        if( params(i)%enabled ) then
            prms(idx) = params(i)%value
            params(i)%pidx = idx
            idx = idx + 1
        else
            if( params(i)%identity .gt. 0 ) then
                params(i)%pidx = params( params(i)%identity )%pidx
            else
                params(i)%pidx = 0
            end if
        end if
    end do

end subroutine ffdev_parameters_gather

! ==============================================================================
! subroutine ffdev_parameters_scatter
! ==============================================================================

subroutine ffdev_parameters_scatter(prms)

    use ffdev_parameters_dat

    implicit none
    real(DEVDP)     :: prms(:)
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,nparams
        if( params(i)%pidx .gt. 0 ) then
            params(i)%value = prms( params(i)%pidx )
        end if
    end do

end subroutine ffdev_parameters_scatter

! ==============================================================================
! subroutine ffdev_parameters_to_tops
! ==============================================================================

subroutine ffdev_parameters_to_tops

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    integer         :: i,j
    ! --------------------------------------------------------------------------

    do i=1,nparams
        select case(params(i)%realm)
            case(REALM_EOFFSET)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%offset = params(i)%value
                    end if
                end do
            case(REALM_BOND_D0)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%bond_types(params(i)%ids(j))%d0 = params(i)%value
                    end if
                end do
            case(REALM_BOND_K)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%bond_types(params(i)%ids(j))%k = params(i)%value
                    end if
                end do
            case(REALM_ANGLE_A0)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%angle_types(params(i)%ids(j))%a0 = params(i)%value
                    end if
                end do
            case(REALM_ANGLE_K)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%angle_types(params(i)%ids(j))%k = params(i)%value
                    end if
                end do
            case(REALM_DIH_V)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%v(params(i)%pn) = params(i)%value
                    end if
                end do
            case(REALM_DIH_C)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%c(params(i)%pn) = params(i)%value
                    end if
                end do
            case(REALM_DIH_G)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%g(params(i)%pn) = params(i)%value
                    end if
                end do
            case(REALM_DIH_SCEE)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%inv_scee = 1.0d0/params(i)%value
                    end if
                end do
            case(REALM_DIH_SCNB)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%inv_scnb = 1.0d0/params(i)%value
                    end if
                end do
            case(REALM_IMPR_V)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%improper_types(params(i)%ids(j))%v = params(i)%value
                    end if
                end do
            case(REALM_IMPR_G)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%improper_types(params(i)%ids(j))%g = params(i)%value
                    end if
                end do
            case(REALM_VDW_A)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%nb_types(params(i)%ids(j))%A = params(i)%value
                    end if
                end do
            case(REALM_VDW_B)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%nb_types(params(i)%ids(j))%B = params(i)%value
                    end if
                end do
            case(REALM_VDW_C)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%nb_types(params(i)%ids(j))%C = params(i)%value
                    end if
                end do
        end select
    end do

end subroutine ffdev_parameters_to_tops

! ==============================================================================
! subroutine ffdev_params_get_lower_bounds
! ==============================================================================

subroutine ffdev_params_get_lower_bounds(tmpx)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    real(DEVDP)     :: tmpx(:)
    ! --------------------------------------------
    integer         :: i, id
    ! --------------------------------------------------------------------------

    id = 0
    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        id = id + 1
        select case(params(i)%realm)
            case(REALM_EOFFSET)
                tmpx(id) = -1000.0
            case(REALM_BOND_D0)
                tmpx(id) = params(i)%value - params(i)%value*0.3
            case(REALM_BOND_K)
                tmpx(id) = 0
            case(REALM_ANGLE_A0)
                tmpx(id) = params(i)%value - params(i)%value*0.3
            case(REALM_ANGLE_K)
                tmpx(id) = 0
            case(REALM_DIH_V)
                tmpx(id) = -100.0
            case(REALM_DIH_C)
                tmpx(id) = -100.0
            case(REALM_DIH_G)
                tmpx(id) = 0.0
            case(REALM_DIH_SCEE)
                tmpx(id) = 0.5
            case(REALM_DIH_SCNB)
                tmpx(id) = 0.5
            case(REALM_IMPR_V)
                tmpx(id) = -20.0
            case(REALM_IMPR_G)
                tmpx(id) = -DEV_PI/180
            case(REALM_VDW_A)
                tmpx(id) = 0.0d0
            case(REALM_VDW_B)
                tmpx(id) = 0.0d0
            case(REALM_VDW_C)
                tmpx(id) = 0.0d0
        end select
    end do

end subroutine ffdev_params_get_lower_bounds

! ==============================================================================
! subroutine ffdev_params_get_lower_bounds
! ==============================================================================

subroutine ffdev_params_get_upper_bounds(tmpx)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    real(DEVDP)     :: tmpx(:)
    ! --------------------------------------------
    integer         :: i, id
    ! --------------------------------------------------------------------------

    id = 0
    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        id = id + 1
        select case(params(i)%realm)
            case(REALM_EOFFSET)
                tmpx(id) = 1000.0d0
            case(REALM_BOND_D0)
                tmpx(id) = params(i)%value + params(i)%value*0.3
            case(REALM_BOND_K)
                tmpx(id) = 1500.0
            case(REALM_ANGLE_A0)
                tmpx(id) = params(i)%value + params(i)%value*0.3
            case(REALM_ANGLE_K)
                tmpx(id) = 1000.0
            case(REALM_DIH_V)
                tmpx(id) = 100.0
            case(REALM_DIH_C)
                tmpx(id) = 100.0
            case(REALM_DIH_G)
                tmpx(id) = 2.0*DEV_PI
            case(REALM_DIH_SCEE)
                tmpx(id) = 3.0
            case(REALM_DIH_SCNB)
                tmpx(id) = 3.0
            case(REALM_IMPR_V)
                tmpx(id) = 20.0
            case(REALM_IMPR_G)
                tmpx(id) = 1.05*DEV_PI
            case(REALM_VDW_A)
                tmpx(id) = 1.0d10
            case(REALM_VDW_B)
                tmpx(id) = 1.0d10
            case(REALM_VDW_C)
                tmpx(id) = 1.0d10
        end select
    end do

end subroutine ffdev_params_get_upper_bounds

! ==============================================================================
! subroutine ffdev_parameters_error
! ==============================================================================

subroutine ffdev_parameters_error(prms,error,grads)

    use ffdev_parameters_dat

    implicit none
    real(DEVDP)         :: prms(:)
    type(FFERROR_TYPE)  :: error
    real(DEVDP)         :: grads(:)
    ! --------------------------------------------------------------------------

    error%total = 0.0d0
    error%energy = 0.0d0
    error%grad = 0.0d0
    error%hess = 0.0d0

    grads(:) = 0.0d0

    call ffdev_parameters_error_num(prms,error,grads)

end subroutine ffdev_parameters_error

! ==============================================================================
! subroutine ffdev_parameters_error_num
! ==============================================================================

subroutine ffdev_parameters_error_num(prms,error,grads)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    real(DEVDP)         :: prms(:)
    type(FFERROR_TYPE)  :: error
    real(DEVDP)         :: grads(:)
    ! --------------------------------------------------------------------------
    real(DEVDP),allocatable :: tmp_prms(:)
    real(DEVDP)             :: d
    type(FFERROR_TYPE)      :: ene1,ene2
    integer                 :: i
    ! --------------------------------------------------------------------------

    d = 1.0d-4  ! differentiation parameter

    ! calculate base energy
    call ffdev_parameters_error_only(prms,error)

    ! allocate temporary geometry object
    allocate( tmp_prms(nactparms) )

    tmp_prms(:) = prms(:)

    ! gradient by numerical differentiation
    do i=1,nactparms
        ! left
        tmp_prms(i) = prms(i) + d
        call ffdev_parameters_error_only(tmp_prms,ene1)

        ! right
        tmp_prms(i) = prms(i) - d
        call ffdev_parameters_error_only(tmp_prms,ene2)

        ! gradient
        grads(i) = 0.5d0*(ene1%total-ene2%total)/d

        ! move back
        tmp_prms(i) = prms(i)
    end do

    ! release temporary geometry object
    deallocate(tmp_prms)

end subroutine ffdev_parameters_error_num

! ==============================================================================
! subroutine ffdev_parameters_error_only
! ==============================================================================

subroutine ffdev_parameters_error_only(prms,error)

    use ffdev_parameters_dat
    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_energy
    use ffdev_gradient
    use ffdev_hessian

    implicit none
    real(DEVDP)         :: prms(:)
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer         :: i,j,q,w,e,r,nene,ngrd,nhess
    real(DEVDP)     :: err,seterrene,seterrgrd,seterrhess,seterrpenalty
    ! --------------------------------------------------------------------------

    error%total = 0.0d0
    error%energy = 0.0d0
    error%grad = 0.0d0
    error%hess = 0.0d0
    error%penalty = 0.0d0

    ! scatter parameters
    call ffdev_parameters_scatter(prms)
    call ffdev_parameters_to_tops()

    ! calculate all energies, gradients, hessians
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_hess_loaded .and. EnableHessianError) then
                call ffdev_hessian_all(sets(i)%top,sets(i)%geo(j))
            else if( sets(i)%geo(j)%trg_grd_loaded .and. EnableGradientError ) then
                call ffdev_gradient_all(sets(i)%top,sets(i)%geo(j))
            else if( sets(i)%geo(j)%trg_ene_loaded .and. EnableEnergyError ) then
                call ffdev_energy_all(sets(i)%top,sets(i)%geo(j))
            end if
        end do
    end do

    ! calculate error
    do i=1,nsets
        seterrene = 0.0
        seterrgrd = 0.0
        seterrhess = 0.0
        seterrpenalty = 0.0
        nene = 0
        ngrd = 0
        nhess = 0
        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_hess_loaded .and. EnableHessianError ) then
                nhess = nhess + 1
                err = 0.0
                do q=1,sets(i)%geo(j)%natoms
                    do w=1,3
                        do e=1,sets(i)%geo(j)%natoms
                            do r=1,3
                                err = err + (sets(i)%geo(j)%hess(r,e,w,q) - sets(i)%geo(j)%trg_hess(r,e,w,q))**2
                            end do
                        end do
                    end do
                end do
                seterrhess = seterrhess + sets(i)%geo(j)%weight * sqrt(err / real(3*sets(i)%geo(j)%natoms)**2 )
            end if
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_grd_loaded .and. EnableGradientError ) then
                ngrd = ngrd + 1
                err = 0.0
                do q=1,sets(i)%geo(j)%natoms
                    do w=1,3
                       err = err + (sets(i)%geo(j)%grd(w,q) - sets(i)%geo(j)%trg_grd(w,q))**2
                    end do
                end do
                seterrgrd = seterrgrd + sets(i)%geo(j)%weight * sqrt(err / real(3*sets(i)%geo(j)%natoms))
            end if
            ! ------------------------------------------------------------------
            if( sets(i)%geo(j)%trg_ene_loaded .and. EnableEnergyError ) then
                nene = nene + 1
                err = sets(i)%geo(j)%total_ene - sets(i)%offset - sets(i)%geo(j)%trg_energy
                seterrene = seterrene + sets(i)%geo(j)%weight * err**2
            end if
            ! ------------------------------------------------------------------
            if( EnablePenaltyError ) then
                err = ffdev_parameters_geo_penalty(sets(i)%top,sets(i)%geo(j))
                seterrpenalty = seterrpenalty + err * sets(i)%geo(j)%weight
            end if
        end do
        if( nene .gt. 0 ) then
            error%energy = error%energy + sqrt(seterrene/real(nene))
        end if
        if( ngrd .gt. 0 ) then
            error%grad = error%grad + seterrgrd/real(ngrd)
        end if
        if( nhess .gt. 0 ) then
            error%hess = error%hess + seterrhess/real(nhess)
        end if
        error%penalty = error%penalty + seterrpenalty / real(sets(i)%ngeos)
    end do

    error%total = EnergyErrorWeight * error%energy &
                + GradientErrorWeight * error%grad &
                + HessianErrorWeight * error%hess &
                + PenaltyErrorWeight * error%penalty

end subroutine ffdev_parameters_error_only

! ==============================================================================
! function ffdev_parameters_geo_penalty
! ==============================================================================

real(DEVDP) function ffdev_parameters_geo_penalty(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_parameters_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ib,i,j,k,ic,ia
    real(DEVDP)     :: rij(3),rkj(3),b,db,bji2inv,bjk2inv,bjiinv,bjkinv,scp
    real(DEVDP)     :: angv,da, penalty
    ! --------------------------------------------------------------------------

    penalty = 0.0

    ! for all bonds
    do ib=1,top%nbonds
        ! for each bond
        i  = top%bonds(ib)%ai
        j  = top%bonds(ib)%aj
        ic = top%bonds(ib)%bt
        ! calculate rij
        rij(:) = geo%crd(:,j) - geo%crd(:,i)

        ! calculate b and db, update energy
        b = sqrt ( rij(1)**2 + rij(2)**2 + rij(3)**2 )
        db = b - top%bond_types(ic)%d0
        penalty = penalty + 0.5*BondD0PenaltyForceK*db**2
    end do

    ! for all angles
    do ia=1,top%nangles
        i  = top%angles(ia)%ai
        j  = top%angles(ia)%aj
        k  = top%angles(ia)%ak
        ic = top%angles(ia)%at

        ! calculate rji and rjk
        rij(:) = geo%crd(:,i) - geo%crd(:,j)
        rkj(:) = geo%crd(:,k) - geo%crd(:,j)

        ! calculate bjiinv and bjkinv and their squares
        bji2inv = 1./(rij(1)**2 + rij(2)**2 + rij(3)**2 )
        bjk2inv = 1./(rkj(1)**2 + rkj(2)**2 + rkj(3)**2 )
        bjiinv = sqrt(bji2inv)
        bjkinv = sqrt(bjk2inv)

            ! calculate scp and angv
        scp = ( rij(1)*rkj(1) + rij(2)*rkj(2) + rij(3)*rkj(3) )
        scp = scp * bjiinv*bjkinv
        if ( scp .gt.  1.0 ) then
            scp =  1.0
        else if ( scp .lt. -1.0 ) then
            scp = -1.0
        end if
        angv = acos(scp)

        ! calculate da and dv
        da = angv - top%angle_types(ic)%a0
        penalty = penalty + 0.5*AngleA0PenaltyForceK*da**2
    end do

    ffdev_parameters_geo_penalty = penalty

end function ffdev_parameters_geo_penalty

! ==============================================================================
! subroutine ffdev_parameters_bond_r0_set
! ==============================================================================

subroutine ffdev_parameters_bond_r0_init(tid,mode)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_geometry
    use prmfile
    use ffdev_utils

    implicit none
    integer         :: tid
    integer         :: mode
    ! -----------------------------------------
    integer         :: num,i,j,ai,aj,bt,k
    real(DEVDP)     :: sum,sum2,ave,lowest,v,rmsf
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    lowest = 1000000
    num = 0
    sum = 0
    sum2 = 0
    rmsf = 0.0
    ave = 0.0
    do i=1,nsets
        if( params(tid)%ids(i) .eq. 0 ) cycle
        bt = params(tid)%ids(i)
        do j=1,sets( i )%ngeos
            do k=1,sets( i )%top%nbonds
                if( sets( i )%top%bonds(k)%bt .ne. bt ) cycle
                ai = sets( i )%top%bonds(k)%ai
                aj = sets( i )%top%bonds(k)%aj
                v = ffdev_geometry_get_length( sets( i )%geo( j ), ai, aj )
                write(DEV_OUT,30) i,bt,j,ai,aj,v
                if( v .le. lowest) then
                    lowest = v
                end if
                sum = sum + v
                sum2 = sum2 + v**2
                num = num + 1
            end do
        end do
    end do
    if( num .gt. 0 ) then
        ave = sum / num
        rmsf = num*sum2 - sum*sum
        rmsf = sqrt(rmsf) / num
    end if

    write(DEV_OUT,40)
    write(DEV_OUT,50) num
    write(DEV_OUT,60) lowest
    write(DEV_OUT,70) ave
    write(DEV_OUT,80) rmsf

    ! set value
    select case(mode)
        case(0)
            params(tid)%value = lowest
        case(1)
            params(tid)%value = ave
    end select

 10 format('Set Type Pts AtomA AtomB   Value  ')
 20 format('--- ---- --- ----- ----- ---------')
 30 format(I3,1X,I4,1x,I3,1X,I5,1X,I5,1X,F9.3)
 40 format('--- ---- --- ----- ----- ---------')
 50 format('Count                    ',I9)
 60 format('Lowest                   ',F9.3)
 70 format('Average                  ',F9.3)
 80 format('RMSF                     ',F9.3)

end subroutine ffdev_parameters_bond_r0_init

! ==============================================================================
! subroutine ffdev_parameters_angle_a0_init
! ==============================================================================

subroutine ffdev_parameters_angle_a0_init(tid,mode)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_geometry
    use prmfile
    use ffdev_utils

    implicit none
    integer         :: tid
    integer         :: mode
    ! -----------------------------------------
    integer         :: num,i,j,ai,aj,ak,at,k
    real(DEVDP)     :: sum,sum2,ave,lowest,v,rmsf
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    lowest = 1000000
    num = 0
    sum = 0
    sum2 = 0
    rmsf = 0.0
    ave = 0.0
    do i=1,nsets
        if( params(tid)%ids(i) .eq. 0 ) cycle
        at = params(tid)%ids(i)
        do j=1,sets( i )%ngeos
            do k=1,sets( i )%top%nangles
                if( sets( i )%top%angles(k)%at .ne. at ) cycle
                ai = sets( i )%top%angles(k)%ai
                aj = sets( i )%top%angles(k)%aj
                ak = sets( i )%top%angles(k)%ak
                v = ffdev_geometry_get_angle( sets( i )%geo( j ), ai, aj, ak )
                write(DEV_OUT,30) i,at,j,ai,aj,ak,v*DEV_R2D
                if( v .le. lowest) then
                    lowest = v
                end if
                sum = sum + v
                sum2 = sum2 + v**2
                num = num + 1
            end do
        end do
    end do
    if( num .gt. 0 ) then
        ave = sum / num
        rmsf = num*sum2 - sum*sum
        rmsf = sqrt(rmsf) / num
    end if

    write(DEV_OUT,40)
    write(DEV_OUT,50) num
    write(DEV_OUT,60) lowest*DEV_R2D
    write(DEV_OUT,70) ave*DEV_R2D
    write(DEV_OUT,80) rmsf*DEV_R2D

    ! set value
    select case(mode)
        case(0)
            params(tid)%value = lowest
        case(1)
            params(tid)%value = ave
    end select

 10 format('Set Type Pts AtomA AtomB AtomC   Value  ')
 20 format('--- ---- --- ----- ----- ----- ---------')
 30 format(I3,1X,I4,1x,I3,1X,I5,1X,I5,1X,I5,1X,F9.3)
 40 format('--- ---- --- ----- ----- ----- ---------')
 50 format('Count                          ',I9)
 60 format('Lowest                         ',F9.3)
 70 format('Average                        ',F9.3)
 80 format('RMSF                           ',F9.3)

end subroutine ffdev_parameters_angle_a0_init

! ------------------------------------------------------------------------------

end module ffdev_parameters

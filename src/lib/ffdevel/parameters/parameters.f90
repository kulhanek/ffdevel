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
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_parameters_init
! ==============================================================================

subroutine ffdev_parameters_init()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_atomicdata
    use ffdev_nb2nb

    implicit none
    integer     :: maxnparams, i, alloc_stat
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'PARAMETERS', ':')
    write(DEV_OUT,5) nsets

    if( Verbosity .ge. DEV_VERBOSITY_MEDIUM ) then
        DEV_OUT = DEV_STD_OUTPUT
    else
        DEV_OUT = DEV_NULL
    end if

    ! print parameters of individual topologies
    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,10) i
        call ffdev_topology_info_types(sets(i)%top)
    end do

    ! extract unique types
    call ffdev_parameters_gen_unique_types()

    ! generate parameters ------------------------------------------------------

    ! determine maximum number of parameters
    maxnparams = 0
    do i=1,nsets
        if( PACAsPrms ) then
            ! sets with references are dependent on referenced topology charges
            ! use only sets with no references
            if( sets(i)%nrefs .eq. 0 )  then
                ! -1 : one class is dependent (total charge constraint))
                maxnparams = maxnparams + (sets(i)%top%nsymm_classes - 1)
            end if
        end if
        ! topology related data
        maxnparams = maxnparams + 2*sets(i)%top%nbond_types     ! bonds
        maxnparams = maxnparams + 2*sets(i)%top%nangle_types    ! angles
        maxnparams = maxnparams + 2*sets(i)%top%ndihedral_types*sets(i)%top%ndihedral_seq_size ! dihedrals
        maxnparams = maxnparams + 2*sets(i)%top%ndihedral_types ! dihedral scee, scnb
        maxnparams = maxnparams + 2*sets(i)%top%nimproper_types ! impropers
        maxnparams = maxnparams + 2*sets(i)%top%nnb_types       ! eps, r0
        maxnparams = maxnparams + 3*sets(i)%top%natom_types     ! pa, pb, rc
        maxnparams = maxnparams + 2*ntypes                      ! zeff, vdw_b0
    end do
    maxnparams = maxnparams + 13 ! ele_qscale, glb_iscee
                                 ! glb_iscnb, disp_s6, disp_s8, disp_s10, damp_fa, damp_fb, damp_pb, damp_tb, damp_pe
                                 ! k_exc, k_ind
    write(DEV_OUT,*)
    write(DEV_OUT,20) maxnparams

    allocate(params(maxnparams), stat = alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for parameter extraction!')
    end if
    do i=1,maxnparams
        allocate(params(i)%ids(nsets), stat = alloc_stat)
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for parameter extraction!')
        end if
    end do

    ! generate parameters
    call ffdev_parameters_reinit()

    ! restore output channel
    DEV_OUT = DEV_STD_OUTPUT

    write(DEV_OUT,30) nparams
    write(DEV_OUT,50)

    call ffdev_nb2nb_init_nbtypes
    call ffdev_parameters_print_types()
    call ffdev_parameters_print_charge_stat()
    call ffdev_parameters_print_parameters(PARAMS_SUMMARY_FULL)

  5 format('Number of sets (topologies)              = ',I6)
 10 format('=== [topology] #',I2.2,' =============================================================')
 20 format('Estimated number of parameters (maximum) = ',I6)
 30 format('Final number of parameters               = ',I6)
 50 format('Input parameters taken from topologies ... ')

end subroutine ffdev_parameters_init

! ==============================================================================
! subroutine ffdev_parameters_reinit
! ==============================================================================

subroutine ffdev_parameters_reinit()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_disp_dat
    use ffdev_atomicdata_dat
    use ffdev_atomicdata_db

    implicit none
    integer     :: i, j, k, parmid, ai, refid
    logical     :: use_vdw_eps, use_vdw_r0
    logical     :: use_vdw_pa, use_vdw_pb, use_vdw_rc
    logical     :: use_ele_sq, use_damp_fa, use_damp_fb, use_damp_pb, use_damp_tb, use_damp_pe
    logical     :: use_disp_s6, use_disp_s8, use_disp_s10, use_k_exc, use_k_ind, use_zeff
    logical     :: use_vdw_b0
    ! --------------------------------------------------------------------------

    nparams = 0

    ! bond length realm ====================
    do i=1,nsets
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
        do j=1,sets(i)%top%ndihedral_types
            if( sets(i)%top%dihedral_types(j)%mode .ne. DIH_COS ) cycle
            do k=1,sets(i)%top%dihedral_types(j)%n
                if( OnlyDefinedDihItems ) then
                    if( sets(i)%top%dihedral_types(j)%enabled(k) .eqv. .false. ) cycle
                else
                    ! switch back to on
                    sets(i)%top%dihedral_types(j)%enabled(k) = .true.
                end if
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
        do j=1,sets(i)%top%ndihedral_types
            if( sets(i)%top%dihedral_types(j)%mode .ne. DIH_GRBF ) cycle
            do k=1,sets(i)%top%dihedral_types(j)%n
                if( OnlyDefinedDihItems ) then
                    if( sets(i)%top%dihedral_types(j)%enabled(k) .eqv. .false. ) cycle
                else
                    ! switch to on
                    sets(i)%top%dihedral_types(j)%enabled(k) = .true.
                end if
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
        do j=1,sets(i)%top%ndihedral_types
            if( sets(i)%top%dihedral_types(j)%mode .ne. DIH_COS ) cycle
            do k=1,sets(i)%top%dihedral_types(j)%n
                if( OnlyDefinedDihItems ) then
                    if( sets(i)%top%dihedral_types(j)%enabled(k) .eqv. .false. ) cycle
                else
                    ! switch to on
                    sets(i)%top%dihedral_types(j)%enabled(k) = .true.
                end if
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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
        if( sets(i)%top%probe_size .ne. 0 ) cycle
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

! ------------------------------------------------------------------------------

    use_vdw_eps     = .false.
    use_vdw_r0      = .false.

    use_vdw_pa      = .false.
    use_vdw_pb      = .false.
    use_vdw_rc      = .false.

    use_ele_sq      = .true.

    use_damp_fa     = .false.
    use_damp_fb     = .false.
    use_damp_pb     = .false.
    use_damp_tb     = .false.
    use_damp_pe     = .false.

    use_disp_s6     = .false.
    use_disp_s8     = .false.
    use_disp_s10    = .false.

    use_k_exc       = .false.
    use_k_ind       = .false.

    use_zeff        = .false.
    use_vdw_b0      = .false.

! common setup
    select case(nb_mode)

        case(NB_VDW_LJ)
            use_vdw_eps     = .true.
            use_vdw_r0      = .true.

        case(NB_VDW_EXP_DISPBJ,NB_VDW_EXP_DISPTT)

            ! induction energy
            if( ind_enabled ) then
                use_k_ind  = .true.
            end if

            ! penetration energy
            if( pen_enabled ) then
                use_damp_pe = .true.
            end if

            use_damp_pb     = .true.

            select case(exp_pa_mode)
                case(EXP_PA_FREEOPT)
                    use_vdw_pa      = .true.
                 case(EXP_PA_CHARGES)
                    use_k_exc       = .true.
                    if( Zeff_mode .eq. AD_ZEFF_OPT ) then
                        use_zeff    = .true.
                    end if
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'exp_pb_mode not implemented in ffdev_parameters_reinit!')
            end select

            select case(exp_pb_mode)
                case(EXP_PB_FREEOPT)
                    use_vdw_pb      = .true.
                 case(EXP_PB_ADBII)
                    if( bii_source .eq. AD_BII_B0OPT ) then
                        use_vdw_b0 = .true.
                    end if
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'exp_pb_mode not implemented in ffdev_parameters_reinit!')
            end select

            use_disp_s6     = .true.
            use_disp_s8     = .true.
            use_disp_s10    = .true.

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_parameters_reinit IIa!')
    end select

! specific setup
    select case(nb_mode)

        case(NB_VDW_LJ)
            ! nothing to be here

        case(NB_VDW_EXP_DISPTT)
            use_damp_tb = .true.

        case(NB_VDW_EXP_DISPBJ)
            select case(dampbj_mode)
                case(DAMP_BJ_FREEOPT)
                    use_vdw_rc      = .true.
                case(DAMP_BJ_CONST)
                    use_damp_fa     = .true.
                case(DAMP_BJ_ADRCII)
                    use_damp_fa     = .true.
                    use_damp_fb     = .true.
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'BJ damp mode not implemented in ffdev_parameters_reinit!')
            end select

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_parameters_reinit IIc!')
    end select

    if( use_vdw_eps ) then
        ! vdw eps realm =====================
        do i=1,nsets
            do j=1,sets(i)%top%nnb_types
                if( .not. ffdev_parameters_is_nbtype_used(sets(i)%top,j) ) cycle
                parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_EPS)
                if( parmid .eq. 0 ) then    ! new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%nb_types(j)%eps
                    params(nparams)%realm = REALM_VDW_EPS
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
    end if

    if( use_vdw_r0 ) then
        ! vdw r0 realm =====================
        do i=1,nsets
            do j=1,sets(i)%top%nnb_types
                if( .not. ffdev_parameters_is_nbtype_used(sets(i)%top,j) ) cycle
                parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_R0)
                if( parmid .eq. 0 ) then    ! new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%nb_types(j)%r0
                    params(nparams)%realm = REALM_VDW_R0
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
    end if

    if( use_vdw_pa ) then
        ! rep PA realm =====================
        do i=1,nsets
            do j=1,sets(i)%top%natom_types
                parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_PA)
                if( parmid .eq. 0 ) then    ! new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%atom_types(j)%PA
                    params(nparams)%realm = REALM_VDW_PA
                    params(nparams)%enabled = .false.
                    params(nparams)%identity = 0
                    params(nparams)%pn    = 0
                    params(nparams)%ids(:) = 0
                    params(nparams)%ids(i) = j
                    params(nparams)%ti   = get_common_type_id(sets(i)%top,j)
                    params(nparams)%tj   = 0
                    params(nparams)%tk   = 0
                    params(nparams)%tl   = 0
                else
                    params(parmid)%ids(i) = j ! parameter already exists, update link
                end if
            end do
        end do
    end if

    if( use_vdw_pb ) then
        ! rep PB realm =====================
        do i=1,nsets
            do j=1,sets(i)%top%natom_types
                parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_PB)
                if( parmid .eq. 0 ) then    ! new parameter
                    ! add new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%atom_types(j)%PB
                    params(nparams)%realm = REALM_VDW_PB
                    params(nparams)%enabled = .false.
                    params(nparams)%identity = 0
                    params(nparams)%pn    = 0
                    params(nparams)%ids(:) = 0
                    params(nparams)%ids(i) = j
                    params(nparams)%ti   = get_common_type_id(sets(i)%top,j)
                    params(nparams)%tj   = 0
                    params(nparams)%tk   = 0
                    params(nparams)%tl   = 0
                else
                    params(parmid)%ids(i) = j ! parameter already exists, update link
                end if
            end do
        end do
    end if

    if( use_vdw_b0 ) then
        ! rep B0 realm =====================
        do i=1,ntypes
            parmid = find_parameter(sets(i)%top,0,types(i)%z,REALM_VDW_B0)
            if( parmid .eq. 0 ) then    ! new parameter
                ! add new parameter
                nparams = nparams + 1
                params(nparams)%value = atomicdata_b0opt(types(i)%z)
                params(nparams)%realm = REALM_VDW_B0
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn      = types(i)%z
                params(nparams)%ids(:)  = 0
                params(nparams)%ti      = 0
                params(nparams)%tj      = 0
                params(nparams)%tk      = 0
                params(nparams)%tl      = 0
            end if
        end do
    end if

    if( use_vdw_rc ) then
        ! =====================
        do i=1,nsets
            do j=1,sets(i)%top%natom_types
                parmid = find_parameter(sets(i)%top,j,0,REALM_VDW_RC)
                if( parmid .eq. 0 ) then    ! new parameter
                    ! add new parameter
                    nparams = nparams + 1
                    params(nparams)%value = sets(i)%top%atom_types(j)%RC
                    params(nparams)%realm = REALM_VDW_RC
                    params(nparams)%enabled = .false.
                    params(nparams)%identity = 0
                    params(nparams)%pn    = 0
                    params(nparams)%ids(:) = 0
                    params(nparams)%ids(i) = j
                    params(nparams)%ti   = get_common_type_id(sets(i)%top,j)
                    params(nparams)%tj   = 0
                    params(nparams)%tk   = 0
                    params(nparams)%tl   = 0
                else
                    params(parmid)%ids(i) = j ! parameter already exists, update link
                end if
            end do
        end do
    end if

    if( use_zeff ) then
        ! =====================
        do j=1,ntypes
            ! add new parameter
            nparams = nparams + 1
            params(nparams)%value = types(j)%Zeff
            params(nparams)%realm = REALM_ZEFF
            params(nparams)%enabled = .false.
            params(nparams)%identity = 0
            params(nparams)%pn    = 0
            params(nparams)%ids(:) = 0
            params(nparams)%ti   = j
            params(nparams)%tj   = 0
            params(nparams)%tk   = 0
            params(nparams)%tl   = 0
        end do
    end if

    if( use_ele_sq ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = ele_qscale
        params(nparams)%realm = REALM_ELE_SQ
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_damp_pb ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = damp_pb
        params(nparams)%realm = REALM_DAMP_PB
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_damp_tb ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = damp_tb
        params(nparams)%realm = REALM_DAMP_TB
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_damp_pe ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = damp_pe
        params(nparams)%realm = REALM_DAMP_PE
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_damp_fa ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = damp_fa
        params(nparams)%realm = REALM_DAMP_FA
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_damp_fb ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = damp_fb
        params(nparams)%realm = REALM_DAMP_FB
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_disp_s6 ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = disp_s6
        params(nparams)%realm = REALM_DISP_S6
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_disp_s8 ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = disp_s8
        params(nparams)%realm = REALM_DISP_S8
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_disp_s10 ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = disp_s10
        params(nparams)%realm = REALM_DISP_S10
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_k_exc ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = k_exc
        params(nparams)%realm = REALM_K_EXC
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    if( use_k_ind ) then
        ! =====================
        nparams = nparams + 1
        params(nparams)%value = k_ind
        params(nparams)%realm = REALM_K_IND
        params(nparams)%enabled = .false.
        params(nparams)%identity = 0
        params(nparams)%pn    = 0
        params(nparams)%ids(:) = 0
        params(nparams)%ti   = 0
        params(nparams)%tj   = 0
        params(nparams)%tk   = 0
        params(nparams)%tl   = 0
    end if

    nparams = nparams + 1
    if( glb_iscee .ne. 0.0d0 ) then
        params(nparams)%value = 1.0d0 / glb_iscee
    else
        params(nparams)%value = 0.0d0
    end if
    params(nparams)%realm = REALM_GLB_SCEE
    params(nparams)%enabled = .false.
    params(nparams)%identity = 0
    params(nparams)%pn    = 0
    params(nparams)%ids(:) = 0
    params(nparams)%ti   = 0
    params(nparams)%tj   = 0
    params(nparams)%tk   = 0
    params(nparams)%tl   = 0

    nparams = nparams + 1
    if( glb_iscnb .ne. 0.0d0 ) then
        params(nparams)%value = 1.0d0 / glb_iscnb
    else
        params(nparams)%value = 0.0d0
    end if
    params(nparams)%realm = REALM_GLB_SCNB
    params(nparams)%enabled = .false.
    params(nparams)%identity = 0
    params(nparams)%pn    = 0
    params(nparams)%ids(:) = 0
    params(nparams)%ti   = 0
    params(nparams)%tj   = 0
    params(nparams)%tk   = 0
    params(nparams)%tl   = 0

    ! charges
    if( PACAsPrms ) then
        ! reset helper index
        do i=1,nsets
            do k=1,sets(i)%top%natoms
                sets(i)%top%atoms(k)%chrg_prm_id = 0
            end do
        end do

        ! first consider non-referenced sets
        do i=1,nsets
            if( sets(i)%nrefs .gt. 0 ) cycle
            ! one symmetry class is dependent (total charge constraint)
            do j=1,sets(i)%top%nsymm_classes-1
                ! add new parameter
                nparams = nparams + 1
                params(nparams)%value = 0.0 ! will be set later
                params(nparams)%realm = REALM_PAC
                params(nparams)%enabled = .false.
                params(nparams)%identity = 0
                params(nparams)%pn     = 0
                params(nparams)%ids(:) = 0
                params(nparams)%ids(i) = j   ! use symmetry class
                params(nparams)%ti     = 0
                params(nparams)%tj     = 0
                params(nparams)%tk     = 0
                params(nparams)%tl     = 0
                do k=1,sets(i)%top%natoms
                    if( sets(i)%top%atoms(k)%symmclass .eq. j ) then
                        sets(i)%top%atoms(k)%chrg_prm_id = nparams
                        params(nparams)%value = sets(i)%top%atoms(k)%charge ! simply use the last value
                    end if
                end do
            end do
        end do

        ! then update sets with references
        do i=1,nsets
            if( sets(i)%nrefs .le. 0 ) cycle
            ai = 1
            do j=1,sets(i)%nrefs
                refid = sets(i)%refs(j)
                do k=1,sets(refid)%top%natoms
                    if( sets(refid)%top%atoms(k)%chrg_prm_id .ne. 0 ) then
                        params(sets(refid)%top%atoms(k)%chrg_prm_id)%ids(i) = sets(i)%top%atoms(ai)%symmclass
                        sets(i)%top%atoms(ai)%chrg_prm_id = sets(refid)%top%atoms(k)%chrg_prm_id
                    end if
                    ai = ai + 1
                end do
            end do
        end do
    end if

    call ffdev_parameters_update_charge_stat()

end subroutine ffdev_parameters_reinit

! ==============================================================================
! function ffdev_parameters_is_nbtype_used
! ==============================================================================

logical function ffdev_parameters_is_nbtype_used(top,nbt)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    type(TOPOLOGY)  :: top
    integer         :: nbt
    ! --------------------------------------------------------------------------

    ffdev_parameters_is_nbtype_used = .false.

    select case(NBParamsMode)
        case(NB_PARAMS_MODE_LIKE_ONLY)
            ffdev_parameters_is_nbtype_used = top%nb_types(nbt)%ti .eq. top%nb_types(nbt)%tj
            return
        case(NB_PARAMS_MODE_ALL)
            ffdev_parameters_is_nbtype_used = .true.
        case default
            call ffdev_utils_exit(DEV_ERR,1,'not implemented in ffdev_parameters_is_nbtype_used!')
    end select

    ! not found
    return

end function ffdev_parameters_is_nbtype_used

! ------------------------------------------------------------------------------

! id, pn
!  0, z   for REALM_VDW_B0


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

    if( realm .eq. REALM_VDW_B0 ) then
        do i=1,nparams
            if( params(i)%realm .ne. realm ) cycle
            ! PN is Z
            if( params(i)%pn .eq. pn ) then
                find_parameter = i
                return
            end if
        end do
        return
    end if

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
        case(REALM_VDW_EPS,REALM_VDW_R0)
            ti = get_common_type_id(top,top%nb_types(id)%ti)
            tj = get_common_type_id(top,top%nb_types(id)%tj)
        case(REALM_VDW_PA,REALM_VDW_PB,REALM_VDW_RC)
            ti = get_common_type_id(top,id)
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in find_parameter I!')
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
            case(REALM_VDW_EPS,REALM_VDW_R0)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj)) .or. &
                    ((params(i)%ti .eq. tj) .and. (params(i)%tj .eq. ti)) ) then
                        find_parameter = i
                end if
            case(REALM_VDW_PA,REALM_VDW_PB,REALM_VDW_RC)
                if( params(i)%ti .eq. ti ) then
                    find_parameter = i
                end if
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in find_parameter II!')
        end select
    end do

end function find_parameter

! ------------------------------------------------------------------------------

integer function find_parameter_by_ids(realm,pn,ti,tj,tk,tl)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer         :: realm
    integer         :: pn
    integer         :: ti,tj,tk,tl
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    find_parameter_by_ids = 0

    ! is parameter defined?
    do i=1,nparams
        if( params(i)%realm .ne. realm ) cycle
        if( params(i)%pn .ne. pn ) cycle

        select case(realm)
            case(REALM_BOND_D0,REALM_BOND_K)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj)) .or. &
                    ((params(i)%ti .eq. tj) .and. (params(i)%tj .eq. ti)) ) then
                        find_parameter_by_ids = i
                        return
                end if
            case(REALM_ANGLE_A0,REALM_ANGLE_K)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj) .and. (params(i)%tk .eq. tk)) .or. &
                    ((params(i)%ti .eq. tk) .and. (params(i)%tj .eq. tj) .and. (params(i)%tk .eq. ti)) ) then
                        find_parameter_by_ids = i
                        return
                end if
            case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G,REALM_DIH_SCEE,REALM_DIH_SCNB)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. tl)) .or. &
                    ((params(i)%ti .eq. tl) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. ti)) ) then
                        find_parameter_by_ids = i
                        return
                end if
            case(REALM_IMPR_V,REALM_IMPR_G)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. tl)) .or. &
                    ((params(i)%ti .eq. tl) .and. (params(i)%tj .eq. tj) .and. &
                     (params(i)%tk .eq. tk) .and. (params(i)%tl .eq. ti)) ) then
                        find_parameter_by_ids = i
                        return
                end if
            case(REALM_VDW_EPS,REALM_VDW_R0)
                if( ((params(i)%ti .eq. ti) .and. (params(i)%tj .eq. tj)) .or. &
                    ((params(i)%ti .eq. tj) .and. (params(i)%tj .eq. ti)) ) then
                        find_parameter_by_ids = i
                        return
                end if
            case(REALM_VDW_PA,REALM_VDW_PB,REALM_VDW_RC,REALM_ZEFF)
                if( params(i)%ti .eq. ti ) then
                    find_parameter_by_ids = i
                    return
                end if
           case(REALM_ELE_SQ,REALM_DAMP_FA,REALM_DAMP_FB, &
                REALM_DAMP_PB,REALM_DAMP_TB,REALM_DAMP_PE, &
                REALM_DISP_S6,REALM_DISP_S8,REALM_DISP_S10,REALM_K_EXC,REALM_K_IND, &
                REALM_GLB_SCEE,REALM_GLB_SCNB)
                find_parameter_by_ids = i
                return
           case(REALM_VDW_B0)
                find_parameter_by_ids = i
                return
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in find_parameter_by_ids!')
        end select
    end do

    ! not found

end function find_parameter_by_ids

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
! ------------------------------------------------------------------------------

integer function get_common_type_from_name(name)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    character(MAX_TNAME)    :: name
    ! --------------------------------------------
    integer                 :: i
    ! --------------------------------------------------------------------------

    get_common_type_from_name = 0

    do i=1,ntypes
        if(trim(types(i)%name) .eq. trim(name) ) then
            get_common_type_from_name = i
            return
        end if
    end do

end function get_common_type_from_name

! ==============================================================================
! subroutine ffdev_parameters_gen_unique_types
! ==============================================================================

subroutine ffdev_parameters_gen_unique_types()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer                             :: maxnparams, i, j, k, alloc_stat
    character(MAX_TNAME),allocatable    :: ltypes(:)
    character(MAX_TNAME)                :: tmp
    logical                             :: changed
    ! --------------------------------------------------------------------------

    ! generate list --------------------
    maxnparams = 0
    do i=1,nsets
        maxnparams = maxnparams + sets(i)%top%natom_types
    end do

    allocate(ltypes(maxnparams), stat = alloc_stat)
    if(alloc_stat .ne. 0) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for all types!')
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
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for common types!')
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

    ! init set occurrence
    do i=1,ntypes
        allocate(types(i)%ids(nsets), stat = alloc_stat)
        if(alloc_stat .ne. 0) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for type ids!')
        end if
        types(i)%ids(:) = 0
        types(i)%probe = .false.
        do j=1,nsets
            do k=1,sets(j)%top%natom_types
                if( types(i)%name .eq. sets(j)%top%atom_types(k)%name ) then
                    types(i)%ids(j) = k
                    types(i)%z = sets(j)%top%atom_types(k)%z
                    types(i)%mass = sets(j)%top%atom_types(k)%mass
                    ! this is cummulative set - if any is probe
                    types(i)%probe = types(i)%probe .or. sets(j)%top%atom_types(k)%probe
                    sets(j)%top%atom_types(k)%glbtypeid = i
                    types(i)%PA   = 0.0d0
                    types(i)%PB   = 0.0d0
                    types(i)%RC   = 0.0d0
                    types(i)%Zeff = types(i)%z
                    exit
                end if
            end do
        end do
    end do

end subroutine ffdev_parameters_gen_unique_types

! ==============================================================================
! subroutine ffdev_parameters_pac_source_to_string
! ==============================================================================

character(80) function ffdev_parameters_pac_source_to_string(pac_source)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer  :: pac_source
    ! --------------------------------------------------------------------------

    select case(pac_source)
        case(PAC_SOURCE_TOPOLOGY)
            ffdev_parameters_pac_source_to_string = 'topology'
        case(PAC_SOURCE_GEO)
            ffdev_parameters_pac_source_to_string = 'geometry'
        case(PAC_SOURCE_GEO_HIRSHFELD)
            ffdev_parameters_pac_source_to_string = 'hirshfeld'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_pac_source_to_string!')
    end select

end function ffdev_parameters_pac_source_to_string

! ==============================================================================
! subroutine ffdev_parameters_pac_source_from_string
! ==============================================================================

integer function ffdev_parameters_pac_source_from_string(string)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('topology')
            ffdev_parameters_pac_source_from_string = PAC_SOURCE_TOPOLOGY
        case('geometry')
            ffdev_parameters_pac_source_from_string = PAC_SOURCE_GEO
        case('hirshfeld')
            ffdev_parameters_pac_source_from_string = PAC_SOURCE_GEO_HIRSHFELD
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_parameters_pac_source_from_string!')
    end select

end function ffdev_parameters_pac_source_from_string

! ==============================================================================
! subroutine ffdev_parameters_update_charge_stat
! ==============================================================================

subroutine ffdev_parameters_update_charge_stat()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer     :: i,j,ti,gti,k
    real(DEVDP) :: q,sdq
    ! --------------------------------------------------------------------------

    ! reset
    types(:)%qcount = 0
    types(:)%minq   = 0
    types(:)%maxq   = 0
    types(:)%aveq   = 0
    types(:)%sdq    = 0

    select case(PACSource)
        case(PAC_SOURCE_TOPOLOGY)
            ! get data
            do i=1,nsets
                do j=1,sets(i)%top%natoms
                    ti  = sets(i)%top%atoms(j)%typeid
                    q   = sets(i)%top%atoms(j)%charge
                    gti = sets(i)%top%atom_types(ti)%glbtypeid
                    types(gti)%qcount = types(gti)%qcount + 1
                    if( types(gti)%qcount .eq. 1 ) then
                        types(gti)%minq = q
                        types(gti)%maxq = q
                    end if
                    if( q .gt. types(gti)%maxq ) types(gti)%maxq = q
                    if( q .lt. types(gti)%minq ) types(gti)%minq = q
                    types(gti)%aveq = types(gti)%aveq + q
                    types(gti)%sdq  = types(gti)%sdq + q**2
                end do
            end do
        case(PAC_SOURCE_GEO)
            ! get data
            do i=1,nsets
                do k=1,sets(i)%ngeos
                    if( .not. sets(i)%geo(k)%sup_chrg_loaded ) cycle
                    do j=1,sets(i)%top%natoms
                        ti  = sets(i)%top%atoms(j)%typeid
                        q   = sets(i)%geo(k)%sup_chrg(j)
                        gti = sets(i)%top%atom_types(ti)%glbtypeid
                        types(gti)%qcount = types(gti)%qcount + 1
                        if( types(gti)%qcount .eq. 1 ) then
                            types(gti)%minq = q
                            types(gti)%maxq = q
                        end if
                        if( q .gt. types(gti)%maxq ) types(gti)%maxq = q
                        if( q .lt. types(gti)%minq ) types(gti)%minq = q
                        types(gti)%aveq = types(gti)%aveq + q
                        types(gti)%sdq  = types(gti)%sdq + q**2
                    end do
                end do
            end do
        case(PAC_SOURCE_GEO_HIRSHFELD)
            ! get data
            do i=1,nsets
                do k=1,sets(i)%ngeos
                    if( .not. sets(i)%geo(k)%sup_hirshfeld_loaded ) cycle
                    do j=1,sets(i)%top%natoms
                        ti  = sets(i)%top%atoms(j)%typeid
                        q   = sets(i)%geo(k)%sup_hirshfeld(j)
                        gti = sets(i)%top%atom_types(ti)%glbtypeid
                        types(gti)%qcount = types(gti)%qcount + 1
                        if( types(gti)%qcount .eq. 1 ) then
                            types(gti)%minq = q
                            types(gti)%maxq = q
                        end if
                        if( q .gt. types(gti)%maxq ) types(gti)%maxq = q
                        if( q .lt. types(gti)%minq ) types(gti)%minq = q
                        types(gti)%aveq = types(gti)%aveq + q
                        types(gti)%sdq  = types(gti)%sdq + q**2
                    end do
                end do
            end do
    end select

    ! finalize
    do gti=1,ntypes
        if( types(gti)%qcount .le. 0 ) cycle
        sdq = types(gti)%qcount * types(gti)%sdq - types(gti)%aveq**2
        if( sdq .gt. 0.d0 ) then
            sdq = sqrt(sdq) / real(types(gti)%qcount)
        else
            sdq = 0.0d0
        end if
        types(gti)%sdq  = sdq
        types(gti)%aveq = types(gti)%aveq / real(types(gti)%qcount)
    end do

end subroutine ffdev_parameters_update_charge_stat

! ==============================================================================
! subroutine ffdev_parameters_print_charge_stat
! ==============================================================================

subroutine ffdev_parameters_print_charge_stat()

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Partial Atomic Charges (PAC) per Atom Types', '=')

    write(DEV_OUT,*)
    write(DEV_OUT,5) trim(ffdev_parameters_pac_source_to_string(PACSource))

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    do i=1,ntypes
        write(DEV_OUT,30) i, adjustl(types(i)%name), types(i)%qcount, &
                          types(i)%minq, types(i)%maxq, types(i)%aveq, types(i)%sdq
    end do

  5 format('PAC source: ',A)

 10 format('# ID Type #Atoms     MinQ     MaxQ      <Q>     s(Q)')
 20 format('# -- ---- ------ -------- -------- -------- --------')
 30 format(I4,1X,A4,1X,I6,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4)

end subroutine ffdev_parameters_print_charge_stat

! ==============================================================================
! subroutine ffdev_parameters_load
! ==============================================================================

subroutine ffdev_parameters_load(name,loaded,ignored)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    character(*)            :: name
    integer                 :: loaded, ignored
    ! --------------------------------------------
    integer                 :: nlp,idx,realm,ti,tj,tk,tl,parmid,pn
    character(MAX_PATH)     :: tmp
    character(MAX_TNAME)    :: sti,stj,stk,stl
    real(DEVDP)             :: value
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_PRMS,name,'O')

    read(DEV_PRMS,*) tmp ! skip the first line
    read(DEV_PRMS,*) tmp, nlp  ! nlp is read but ignored
    if( trim(tmp) .ne. 'v1' ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unsupported file format version, it should be v1!')
    end if

    ! read file by lines
    do while(.true.)
        read(DEV_PRMS,*,end=100,err=110) idx, tmp, sti, stj, stk, stl, pn, value
        realm = ffdev_parameters_get_realmid(tmp)
        ti = get_common_type_from_name(sti)
        tj = get_common_type_from_name(stj)
        tl = get_common_type_from_name(stl)
        tk = get_common_type_from_name(stk)

        parmid = find_parameter_by_ids(realm,pn,ti,tj,tk,tl)
        if( parmid .ne. 0 ) then
            params(parmid)%value = value
            loaded = loaded + 1
        else
            ignored = ignored + 1
        end if
    end do

100    close(DEV_PRMS)

    return

110 call ffdev_utils_exit(DEV_ERR,1,'Unable to read data from FFDevel parameter file!')

end subroutine ffdev_parameters_load

! ==============================================================================
! subroutine ffdev_parameters_save
! ==============================================================================

subroutine ffdev_parameters_save(name)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_ver

    implicit none
    character(*)            :: name
    ! --------------------------------------------
    integer                 :: i
    character(MAX_TNAME)    :: sti,stj,stk,stl
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_PRMS,name,'U')

    write(DEV_PRMS,5)  trim(FFDEV_LIBVER)
    write(DEV_PRMS,10) nparams

    do i=1,nparams
        sti = '--'
        stj = '--'
        stk = '--'
        stl = '--'
        if( params(i)%ti .ne. 0 ) sti = types(params(i)%ti)%name
        if( params(i)%tj .ne. 0 ) stj = types(params(i)%tj)%name
        if( params(i)%tk .ne. 0 ) stk = types(params(i)%tk)%name
        if( params(i)%tl .ne. 0 ) stl = types(params(i)%tl)%name
        write(DEV_PRMS,20) i, ffdev_parameters_get_realm_name(params(i)%realm), &
                           sti,stj,stk,stl,params(i)%pn, params(i)%value
    end do

    close(DEV_PRMS)

  5 format('# parameters created by ffdevel version : ',A66)
 10 format('v1',1X,I6)
 20 format(I3,1X,A20,1X,A4,1X,A4,1X,A4,1X,A4,1X,I2,1X,E22.16)

end subroutine ffdev_parameters_save

! ==============================================================================
! subroutine ffdev_parameters_save_amber
! FIXME - better handle dih_v/dih_g/dih_c and dih_scee, dih_scnb
! ==============================================================================

subroutine ffdev_parameters_save_amber(name)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_nb2nb
    use ffdev_nb2nb_dat

    implicit none
    character(*)    :: name
    ! --------------------------------------------
    integer         :: datum(8)
    integer         :: i,j,it,ij,pn,max_pn,idx,si
    real(DEVDP)     :: v,g, scee, scnb
    logical         :: enable_section,old_enable_section
    ! --------------------------------------------------------------------------

    call ffdev_utils_open(DEV_PRMS,name,'U')

    call date_and_time(values=datum)
    write(DEV_PRMS,10) datum(1),datum(2),datum(3),datum(5),datum(6),datum(7)

    ! atom types
    enable_section = .false.
    do i=1,ntypes
        types(i)%print_nb = .false.
    end do
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_BOND_D0) .or. ( params(i)%realm .eq. REALM_BOND_K ) ) then
            enable_section = .true.
            types(params(i)%ti)%print_nb = .true.
            types(params(i)%tj)%print_nb = .true.
            cycle
        end if
        if( (params(i)%realm .eq. REALM_ANGLE_A0) .or. ( params(i)%realm .eq. REALM_ANGLE_K ) ) then
            enable_section = .true.
            types(params(i)%ti)%print_nb = .true.
            types(params(i)%tj)%print_nb = .true.
            types(params(i)%tk)%print_nb = .true.
            cycle
        end if
        if( (params(i)%realm .eq. REALM_DIH_V) .or. ( params(i)%realm .eq. REALM_DIH_G ) .or. &
            (params(i)%realm .eq. REALM_DIH_SCEE) .or. ( params(i)%realm .eq. REALM_DIH_SCNB ) ) then
            enable_section = .true.
            types(params(i)%ti)%print_nb = .true.
            types(params(i)%tj)%print_nb = .true.
            types(params(i)%tk)%print_nb = .true.
            types(params(i)%tl)%print_nb = .true.
            cycle
        end if
        if( (params(i)%realm .eq. REALM_IMPR_V) .or. ( params(i)%realm .eq. REALM_IMPR_G ) ) then
            enable_section = .true.
            types(params(i)%ti)%print_nb = .true.
            types(params(i)%tj)%print_nb = .true.
            types(params(i)%tk)%print_nb = .true.
            types(params(i)%tl)%print_nb = .true.
            cycle
        end if
        if( (params(i)%realm .eq. REALM_VDW_EPS) .or. ( params(i)%realm .eq. REALM_VDW_R0 ) ) then
            enable_section = .true.
            types(params(i)%ti)%print_nb = .true.
            cycle
        end if
    end do
    if( enable_section ) then
        write(DEV_PRMS,20) 'MASS'
        do i=1,ntypes
            if( types(i)%print_nb ) then
                write(DEV_PRMS,30) types(i)%name,types(i)%mass,0.0d0
            end if
        end do
        write(DEV_PRMS,*)
    end if

    ! bonds
    enable_section = .false.
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_BOND_D0) .or. ( params(i)%realm .eq. REALM_BOND_K ) ) then
            enable_section = .true.
            exit
        end if
    end do
    if( enable_section ) then
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
                               0.5d0*params(j)%value,params(i)%value
        end do
        write(DEV_PRMS,*)
    end if

    ! angles
    enable_section = .false.
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_ANGLE_A0) .or. ( params(i)%realm .eq. REALM_ANGLE_K ) ) then
            enable_section = .true.
            exit
        end if
    end do
    if( enable_section ) then
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
                               0.5d0*params(j)%value,params(i)%value*DEV_R2D
        end do
        write(DEV_PRMS,*)
    end if

    ! dihedrals
    enable_section = .false.
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_DIH_V) .or. ( params(i)%realm .eq. REALM_DIH_G ) .or. &
            (params(i)%realm .eq. REALM_DIH_SCEE) .or. ( params(i)%realm .eq. REALM_DIH_SCNB ) ) then
            enable_section = .true.
            exit
        end if
    end do
    if( enable_section ) then
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
    end if

    ! dih_c realm
    old_enable_section = enable_section
    enable_section = .false.
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_DIH_C) .or. &
            (params(i)%realm .eq. REALM_DIH_SCEE) .or. ( params(i)%realm .eq. REALM_DIH_SCNB ) ) then
            enable_section = .true.
            if( .not. old_enable_section ) then
                write(DEV_PRMS,20) 'DIHE'
            end if
            exit
        end if
    end do
    if( old_enable_section .and. (.not. enable_section) ) then
        write(DEV_PRMS,*)
    end if

    if( enable_section ) then
        it = 0
        do i=1,nparams
            if( params(i)%realm .ne. REALM_DIH_C ) cycle
            if( params(i)%pn .ne. 1 ) cycle

            ! find topology for transformation
            si  = 0
            idx = 0
            do j=1,nsets
                if( params(i)%ids(j) .ne. 0 ) then
                    idx = params(i)%ids(j)
                    si = j
                    exit
                end if
            end do
            if( (idx .eq. 0) .or. (si .eq. 0) ) cycle

            ! transform the dihedral
            call ffdev_parameters_grbf2cos(sets(si)%top,idx)

            ! find max pn
            max_pn = 0
            do j=1,sets(si)%top%dihedral_types(idx)%n
                if( sets(si)%top%dihedral_types(idx)%enabled(j) ) then
                    max_pn = j
                end if
            end do

            ! write dihedral
            scee = 1.0d0 / sets(si)%top%dihedral_types(idx)%inv_scee
            scnb = 1.0d0 / sets(si)%top%dihedral_types(idx)%inv_scnb
            do j=1,sets(si)%top%dihedral_types(idx)%n
                if( .not. sets(si)%top%dihedral_types(idx)%enabled(j) ) cycle
                pn = j
                if( pn .lt. max_pn ) pn = -pn
                write(DEV_PRMS,60) sets(si)%top%atom_types(sets(si)%top%dihedral_types(idx)%ti)%name, &
                                   sets(si)%top%atom_types(sets(si)%top%dihedral_types(idx)%tj)%name, &
                                   sets(si)%top%atom_types(sets(si)%top%dihedral_types(idx)%tk)%name, &
                                   sets(si)%top%atom_types(sets(si)%top%dihedral_types(idx)%tl)%name, &
                                   sets(si)%top%dihedral_types(idx)%v(j), &
                                   sets(si)%top%dihedral_types(idx)%g(j)*DEV_R2D, pn, scee, scnb
            end do
        end do
        write(DEV_PRMS,*)
    end if

    ! impropers
    enable_section = .false.
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_IMPR_V) .or. ( params(i)%realm .eq. REALM_IMPR_G ) ) then
            enable_section = .true.
            exit
        end if
    end do
    if( enable_section ) then
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
    end if

    enable_section = .false.
    do i=1,nparams
        if( (params(i)%realm .eq. REALM_VDW_EPS) .or. ( params(i)%realm .eq. REALM_VDW_R0 ) ) then
            enable_section = .true.
            exit
        end if
    end do

    if( enable_section ) then
        ! reset data
        do i=1,ntypes
            types(i)%r0 = 0.0d0
            types(i)%eps = 0.0d0
        end do
        ! populate data
        do i=1,nnb_types
            if( nb_types(i)%gti .eq. nb_types(i)%gtj ) then
                types(i)%r0  = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%r0
                types(i)%eps = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%eps
            end if
        end do
        ! write data
        write(DEV_PRMS,20) 'NONB'
        do i=1,ntypes
            if( .not. types(i)%print_nb ) cycle
            write(DEV_PRMS,80) types(i)%name,types(i)%r0*0.5d0,types(i)%eps
        end do
        write(DEV_PRMS,*)
    end if

    close(DEV_PRMS)

 10 format('# parameters generated by ffdevel package ',i4,'-',i2.2,'-',i2.2,' ',i2,':',i2.2,':',i2.2)
 20 format(A)
 30 format(A4,10X,F16.6,1X,F16.6)
 40 format(A2,'-',A2,9X,F16.6,1X,F16.6)
 50 format(A2,'-',A2,'-',A2,6X,F16.6,1X,F16.3)
 60 format(A2,'-',A2,'-',A2,'-',A2,1X,'1',1X,F16.6,1X,F16.3,1X,I4,1X,'SCEE=',F4.2,1X,'SCNB=',F4.2)
 70 format(A2,'-',A2,'-',A2,'-',A2,3X,F16.6,1X,F16.3,1X,I4)
 80 format(A4,10X,F16.6,1X,F16.6)

end subroutine ffdev_parameters_save_amber

! ==============================================================================
! subroutine ffdev_parameters_grbf2cos
! readings:
! https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/
! ==============================================================================

subroutine ffdev_parameters_grbf2cos(top,idx)

    use ffdev_topology
    use ffdev_parameters_dat
    use ffdev_geometry
    use ffdev_utils

    implicit none
    include 'fftw3.f'
    type(TOPOLOGY)          :: top      ! input topology
    integer                 :: idx      ! dihedral index
    ! --------------------------------------------
    real(DEVDP),allocatable     :: x(:)
    complex(DEVDP),allocatable  :: y(:)
    integer                     :: i, pn, alloc_stat
    real(DEVDP)                 :: phi, ene, offset, rmse, arg
    integer*8                   :: plan
    ! --------------------------------------------------------------------------

    if( top%dihedral_types(idx)%mode .ne. DIH_GRBF ) then
        call ffdev_utils_exit(DEV_ERR,1,'Dihedral is not DIH_GRBF in ffdev_parameters_grbf2cos!')
    end if

    write(DEV_OUT,10,ADVANCE='NO')  top%atom_types(top%dihedral_types(idx)%ti)%name, &
                                    top%atom_types(top%dihedral_types(idx)%tj)%name, &
                                    top%atom_types(top%dihedral_types(idx)%tk)%name, &
                                    top%atom_types(top%dihedral_types(idx)%tl)%name

    allocate(x(GRBF2COSMaxN),y(GRBF2COSMaxN/2+1), stat = alloc_stat)
    if(alloc_stat .ne. 0) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate memory for FFTW in ffdev_parameters_grbf2cos!')
    end if

    ! calculate the dihedral potential
    do i=1,GRBF2COSMaxN
        phi = 2.0d0*DEV_PI*(i-1)/(real(GRBF2COSMaxN))
        ene = 0.0d0
        do pn=1,top%dihedral_types(idx)%n
            if( .not. top%dihedral_types(idx)%enabled(pn) ) cycle
            ene = ene + top%dihedral_types(idx)%c(pn) &
                          * exp(-(ffdev_geometry_get_dihedral_deviation(phi,top%dihedral_types(idx)%p(pn))**2) &
                          / top%dihedral_types(idx)%w2(pn))
        end do
        x(i) = ene
    end do

    ! run FFT
    call dfftw_plan_dft_r2c_1d(plan,GRBF2COSMaxN,x,y,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, x, y)
    call dfftw_destroy_plan(plan)

    ! filter them and update dihedral_type
    top%dihedral_types(idx)%mode = DIH_COS
    top%dihedral_types(idx)%enabled(:) = .false.
    do i=1,GRBF2COSMaxN
        if( 2.0d0*abs(y(i+1))/real(GRBF2COSMaxN) .gt. GRBF2COSMinV ) then
            top%dihedral_types(idx)%enabled(i) = .true.
            top%dihedral_types(idx)%g(i) = 2*DEV_PI - atan2(aimag(y(i+1)),real(y(i+1)))
            ! wrap phase into <0;360>
            top%dihedral_types(idx)%g(i) = top%dihedral_types(idx)%g(i) &
                                         - 2.0d0*DEV_PI*floor(top%dihedral_types(idx)%g(i)/(2.0d0*DEV_PI))
            top%dihedral_types(idx)%v(i) = 2.0d0*abs(y(i+1))/real(GRBF2COSMaxN)
        end if
    end do

    ! calculate offset
    offset = 0.0d0
    do i=1,GRBF2COSMaxN
        phi = 2.0d0*DEV_PI*(i-1)/(real(GRBF2COSMaxN))
        ene = 0.0d0
        do pn=1,top%dihedral_types(idx)%n
            if( .not. top%dihedral_types(idx)%enabled(pn) ) cycle
            arg = pn*phi - top%dihedral_types(idx)%g(pn)
            if( dih_cos_only ) then
                ene = ene + top%dihedral_types(idx)%v(pn)*cos(arg)
            else
                ene = ene + top%dihedral_types(idx)%v(pn)*(1.0d0+cos(arg))
            end if
        end do
        offset = offset + ene - x(i)
    end do
    offset = offset / real(GRBF2COSMaxN)

    ! write(*,*) 'offset=',offset

    ! calculate error
    rmse = 0.0d0
    do i=1,GRBF2COSMaxN
        phi = 2.0d0*DEV_PI*(i-1)/(real(GRBF2COSMaxN))
        ene = 0.0d0
        do pn=1,top%dihedral_types(idx)%n
            if( .not. top%dihedral_types(idx)%enabled(pn) ) cycle
            arg = pn*phi - top%dihedral_types(idx)%g(pn)
            if( dih_cos_only ) then
                ene = ene + top%dihedral_types(idx)%v(pn)*cos(arg)
            else
                ene = ene + top%dihedral_types(idx)%v(pn)*(1.0d0+cos(arg))
            end if
        end do
        ! write(*,*) ene- offset, x(i)
        rmse = rmse + (ene - x(i) - offset)**2
    end do
    rmse = sqrt(rmse / real(GRBF2COSMaxN))

    write(DEV_OUT,20) rmse

    deallocate(x,y)

 10 format('# converting grbf2cos for: ',A2,'-',A2,'-',A2,'-',A2)
 20 format(', RMSE= ',F10.4)

end subroutine ffdev_parameters_grbf2cos

! ==============================================================================
! subroutine ffdev_parameters_print_types
! ==============================================================================

subroutine ffdev_parameters_print_types()

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer     :: i,j,icount
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'Atom Type Summary', '=')
    write(DEV_OUT,5) ntypes

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    do i=1,ntypes
        icount = 0
        do j=1,nsets
            if( types(i)%ids(j) .gt. 0 ) icount = icount + 1
        end do
        write(DEV_OUT,30,ADVANCE='NO') i, types(i)%name, icount
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

 10 format('# ID Type  #Sets     IDs in Sets')
 20 format('# -- ---- ------     --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---')
 30 format(I4,1X,A4,1X,I6,5X)
 40 format(I3,1X)
 50 format('---',1X)

end subroutine ffdev_parameters_print_types

! ==============================================================================
! subroutine ffdev_parameters_get_realmid
! ==============================================================================

integer function ffdev_parameters_get_realmid(realm)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    character(MAX_PATH) :: realm
    ! --------------------------------------------------------------------------

    select case(trim(realm))
        case('bond_r0')
            ffdev_parameters_get_realmid = REALM_BOND_D0
        case('bond_k')
            ffdev_parameters_get_realmid = REALM_BOND_K
        case('angle_a0')
            ffdev_parameters_get_realmid = REALM_ANGLE_A0
        case('angle_k')
            ffdev_parameters_get_realmid = REALM_ANGLE_K
        case('dih_v')
            ffdev_parameters_get_realmid = REALM_DIH_V
        case('dih_c')
            ffdev_parameters_get_realmid = REALM_DIH_C
        case('dih_gamma')
            ffdev_parameters_get_realmid = REALM_DIH_G
        case('dih_scee')
            ffdev_parameters_get_realmid = REALM_DIH_SCEE
        case('dih_scnb')
            ffdev_parameters_get_realmid = REALM_DIH_SCNB
        case('impr_v')
            ffdev_parameters_get_realmid = REALM_IMPR_V
        case('impr_gamma')
            ffdev_parameters_get_realmid = REALM_IMPR_G

        case('vdw_eps')
            ffdev_parameters_get_realmid = REALM_VDW_EPS
        case('vdw_r0')
            ffdev_parameters_get_realmid = REALM_VDW_R0

        case('vdw_pa')
            ffdev_parameters_get_realmid = REALM_VDW_PA
        case('vdw_pb')
            ffdev_parameters_get_realmid = REALM_VDW_PB
        case('vdw_b0')
            ffdev_parameters_get_realmid = REALM_VDW_B0
        case('vdw_rc')
            ffdev_parameters_get_realmid = REALM_VDW_RC

        case('ele_sq')
            ffdev_parameters_get_realmid = REALM_ELE_SQ

        case('damp_fa')
            ffdev_parameters_get_realmid = REALM_DAMP_FA
        case('damp_fb')
            ffdev_parameters_get_realmid = REALM_DAMP_FB
        case('damp_pb')
            ffdev_parameters_get_realmid = REALM_DAMP_PB
        case('damp_tb')
            ffdev_parameters_get_realmid = REALM_DAMP_TB
        case('damp_pe')
            ffdev_parameters_get_realmid = REALM_DAMP_PE

        case('disp_s6')
            ffdev_parameters_get_realmid = REALM_DISP_S6
        case('disp_s8')
            ffdev_parameters_get_realmid = REALM_DISP_S8
        case('disp_s10')
            ffdev_parameters_get_realmid = REALM_DISP_S10

        case('pac')
            ffdev_parameters_get_realmid = REALM_PAC

        case('zeff')
            ffdev_parameters_get_realmid = REALM_ZEFF

        case('glb_scee')
            ffdev_parameters_get_realmid = REALM_GLB_SCEE

        case('glb_scnb')
            ffdev_parameters_get_realmid = REALM_GLB_SCNB

        case('k_exc')
            ffdev_parameters_get_realmid = REALM_K_EXC

        case('k_ind')
            ffdev_parameters_get_realmid = REALM_K_IND

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_get_realmid (' // trim(realm) // ')!')
    end select

end function ffdev_parameters_get_realmid

! ==============================================================================
! subroutine ffdev_parameters_get_realm_name
! ==============================================================================

character(MAX_PATH) function ffdev_parameters_get_realm_name(realmid)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer             :: realmid
    ! --------------------------------------------------------------------------

    select case(realmid)
        case(REALM_BOND_D0)
            ffdev_parameters_get_realm_name = 'bond_r0'
        case(REALM_BOND_K)
            ffdev_parameters_get_realm_name = 'bond_k'
        case(REALM_ANGLE_A0)
            ffdev_parameters_get_realm_name = 'angle_a0'
        case(REALM_ANGLE_K)
            ffdev_parameters_get_realm_name = 'angle_k'
        case(REALM_DIH_V)
            ffdev_parameters_get_realm_name = 'dih_v'
        case(REALM_DIH_C)
            ffdev_parameters_get_realm_name = 'dih_c'
        case(REALM_DIH_G)
            ffdev_parameters_get_realm_name = 'dih_gamma'
        case(REALM_DIH_SCEE)
            ffdev_parameters_get_realm_name = 'dih_scee'
        case(REALM_DIH_SCNB)
            ffdev_parameters_get_realm_name = 'dih_scnb'
        case(REALM_IMPR_V)
            ffdev_parameters_get_realm_name = 'impr_v'
        case(REALM_IMPR_G)
            ffdev_parameters_get_realm_name = 'impr_gamma'

        case(REALM_VDW_EPS)
            ffdev_parameters_get_realm_name = 'vdw_eps'
        case(REALM_VDW_R0)
            ffdev_parameters_get_realm_name = 'vdw_r0'

        case(REALM_VDW_PA)
            ffdev_parameters_get_realm_name = 'vdw_pa'
        case(REALM_VDW_PB)
            ffdev_parameters_get_realm_name = 'vdw_pb'
        case(REALM_VDW_B0)
            ffdev_parameters_get_realm_name = 'vdw_b0'
        case(REALM_VDW_RC)
            ffdev_parameters_get_realm_name = 'vdw_rc'

        case(REALM_ELE_SQ)
            ffdev_parameters_get_realm_name = 'ele_sq'
        case(REALM_PAC)
            ffdev_parameters_get_realm_name = 'pac'
        case(REALM_ZEFF)
            ffdev_parameters_get_realm_name = 'zeff'

        case(REALM_DAMP_FA)
            ffdev_parameters_get_realm_name = 'damp_fa'
        case(REALM_DAMP_FB)
            ffdev_parameters_get_realm_name = 'damp_fb'
        case(REALM_DAMP_PB)
            ffdev_parameters_get_realm_name = 'damp_pb'
        case(REALM_DAMP_TB)
            ffdev_parameters_get_realm_name = 'damp_tb'
        case(REALM_DAMP_PE)
            ffdev_parameters_get_realm_name = 'damp_pe'

        case(REALM_DISP_S6)
            ffdev_parameters_get_realm_name = 'disp_s6'
        case(REALM_DISP_S8)
            ffdev_parameters_get_realm_name = 'disp_s8'
        case(REALM_DISP_S10)
            ffdev_parameters_get_realm_name = 'disp_s10'

        case(REALM_GLB_SCEE)
            ffdev_parameters_get_realm_name = 'glb_scee'
        case(REALM_GLB_SCNB)
            ffdev_parameters_get_realm_name = 'glb_scnb'

        case(REALM_K_EXC)
            ffdev_parameters_get_realm_name = 'k_exc'

        case(REALM_K_IND)
            ffdev_parameters_get_realm_name = 'k_ind'

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_get_realm_name!')
    end select

end function ffdev_parameters_get_realm_name

! ==============================================================================
! subroutine ffdev_parameters_get_realm_scaling
! ==============================================================================

real(DEVDP) function ffdev_parameters_get_realm_scaling(realmid)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer             :: realmid
    ! --------------------------------------------------------------------------

    ffdev_parameters_get_realm_scaling = 1.0d0

    select case(realmid)
        case(REALM_BOND_D0,REALM_BOND_K)
            ! nothing to do
        case(REALM_ANGLE_A0)
            ffdev_parameters_get_realm_scaling = DEV_R2D
        case(REALM_ANGLE_K,REALM_DIH_V,REALM_DIH_C)
            ! nothing to do
        case(REALM_DIH_G)
            ffdev_parameters_get_realm_scaling = DEV_R2D
        case(REALM_DIH_SCEE,REALM_DIH_SCNB,REALM_IMPR_V)
            ! nothing to do
        case(REALM_IMPR_G)
            ffdev_parameters_get_realm_scaling = DEV_R2D
        case(REALM_VDW_EPS,REALM_VDW_R0)
            ! nothing to do
        case(REALM_VDW_PA,REALM_VDW_PB,REALM_VDW_B0,REALM_VDW_RC)
            ! nothing to do
        case(REALM_ELE_SQ,REALM_PAC,REALM_ZEFF)
            ! nothing to do
        case(REALM_DAMP_FA,REALM_DAMP_FB,REALM_DAMP_PB,REALM_DAMP_TB,REALM_DAMP_PE)
            ! nothing to do
        case(REALM_DISP_S6,REALM_DISP_S8,REALM_DISP_S10)
            ! nothing to do
        case(REALM_GLB_SCEE,REALM_GLB_SCNB)
            ! nothing to do
        case(REALM_K_EXC,REALM_K_IND)
            ! nothing to do
        case default
            write(DEV_ERR,*) 'realmid = ', realmid
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_get_realm_scaling!')
    end select

end function ffdev_parameters_get_realm_scaling

! ==============================================================================
! subroutine ffdev_parameters_print_parameters
! ==============================================================================

subroutine ffdev_parameters_print_parameters(mode)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer                 :: mode
    integer                 :: i,j,count,free,act,tot
    character(MAX_PATH)     :: tmp
    character(MAX_TNAME)    :: sti,stj,stk,stl
    real(DEVDP)             :: scaling
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    select case(mode)
        case(PARAMS_SUMMARY_INITIAL)
            call ffdev_utils_heading(DEV_OUT,'Input Parameter Summary', ':')
            ! print initial parameter ranges
            call ffdev_params_print_ranges

        case(PARAMS_SUMMARY_OPTIMIZED)
            call ffdev_utils_heading(DEV_OUT,'Optimized Parameter Summary', ':')
            ! print initial parameter ranges
            call ffdev_params_print_ranges

        case(PARAMS_SUMMARY_MODIFIED)
            call ffdev_utils_heading(DEV_OUT,'Modified Parameter Summary', ':')

        case(PARAMS_SUMMARY_FULL)
            call ffdev_utils_heading(DEV_OUT,'Parameter Summary', ':')

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_print_parameters!')
    end select

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    do i=1,nparams

        ! print only enabled parameters for initial/final optimization summaries
        if( .not. params(i)%enabled ) then
            if( mode .eq. PARAMS_SUMMARY_INITIAL ) cycle
            if( mode .eq. PARAMS_SUMMARY_OPTIMIZED ) cycle
        end if

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

        tmp = ffdev_parameters_get_realm_name(params(i)%realm)
        scaling = ffdev_parameters_get_realm_scaling(params(i)%realm)
        write(DEV_OUT,32,ADVANCE='NO') adjustl(tmp)

        sti = '--'
        if( params(i)%ti .ne. 0 ) sti = types(params(i)%ti)%name
        stj = '--'
        if( params(i)%tj .ne. 0 ) stj = types(params(i)%tj)%name
        stk = '--'
        if( params(i)%tk .ne. 0 ) stk = types(params(i)%tk)%name
        stl = '--'
        if( params(i)%tl .ne. 0 ) stl = types(params(i)%tl)%name

        write(DEV_OUT,35,ADVANCE='NO') trim(sti), trim(stj), trim(stk), trim(stl), &
                                        params(i)%pn, params(i)%value*scaling, count
        do j=1,nsets
            if( params(i)%ids(j) .gt. 0 ) then
                write(DEV_OUT,40,ADVANCE='NO') params(i)%ids(j)
            else
                write(DEV_OUT,50,ADVANCE='NO')
            end if
        end do
        write(DEV_OUT,*)
    end do

    free = 0
    act = 0
    do i=1,nparams
        if( params(i)%identity .eq. 0 ) free = free + 1
        if( params(i)%enabled ) act = act + 1
    end do

    if( mode .ne. PARAMS_SUMMARY_FULL ) then
        write(DEV_OUT,*)
        write(DEV_OUT,200) nparams
        if( free .ne. nparams ) then
        write(DEV_OUT,210) free
        end if
        write(DEV_OUT,230) act
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)

    do j=REALM_FIRST,REALM_LAST
        tot = 0
        free = 0
        act = 0

        do i=1,nparams
            if( params(i)%realm .ne. j ) cycle
            tot = tot + 1
            if( params(i)%identity .eq. 0 ) free = free + 1
            if( params(i)%enabled ) act = act + 1
        end do

        tmp = ffdev_parameters_get_realm_name(j)
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


 10 format('# ID ST Iden    Realm    TI TJ TK TL PN       Value      Counts     IDs in Sets')
 20 format('# -- -- ---- ----------- -- -- -- -- -- ---------------- ------     --- --- ---')
 30 format(I4,1X,L2,1X,I4,1X)
 31 format(I4,1X,L2,1X,'----',1X)
 32 format(A11,1X)


 35 format(A2,1X,A2,1X,A2,1X,A2,1X,I2,1X,F16.4,1X,I6,5X)
 40 format(I3,1X)
 50 format('---',1X)

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
    use ffdev_utils
    use ffdev_atomicdata_db

    implicit none
    integer         :: i,j,k,ai,ndep,refid
    real(DEVDP)     :: chrg
    ! --------------------------------------------------------------------------

    ! disable all types in target topologies
    do i=1,nsets
        do j=1,sets(i)%top%nbond_types
            sets(i)%top%bond_types(j)%ffoptactive = .false.
        end do
        do j=1,sets(i)%top%nangle_types
            sets(i)%top%angle_types(j)%ffoptactive = .false.
        end do
        do j=1,sets(i)%top%ndihedral_types
            sets(i)%top%dihedral_types(j)%ffoptactive = .false.
        end do
        do j=1,sets(i)%top%nimproper_types
            sets(i)%top%improper_types(j)%ffoptactive = .false.
        end do
        do j=1,sets(i)%top%nnb_types
            sets(i)%top%nb_types(j)%ffoptactive = .false.
        end do
        do j=1,sets(i)%top%natoms
            sets(i)%top%atoms(j)%ffoptactive = .false.
        end do
        do j=1,sets(i)%top%natom_types
            sets(i)%top%atom_types(j)%ffoptactive = .false.
        end do
    end do

    do i=1,nparams
        select case(params(i)%realm)

        ! bonded parameters
            case(REALM_BOND_D0)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%bond_types(params(i)%ids(j))%d0 = params(i)%value
                        sets(j)%top%bond_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_BOND_K)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%bond_types(params(i)%ids(j))%k = params(i)%value
                        sets(j)%top%bond_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_ANGLE_A0)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%angle_types(params(i)%ids(j))%a0 = params(i)%value
                        sets(j)%top%angle_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_ANGLE_K)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%angle_types(params(i)%ids(j))%k = params(i)%value
                        sets(j)%top%angle_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_DIH_V)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%v(params(i)%pn) = params(i)%value
                        ! this must be accumulatiove do to pn
                        sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive = &
                            sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive .or. params(i)%enabled
                    end if
                end do
            case(REALM_DIH_C)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%c(params(i)%pn) = params(i)%value
                        ! this must be accumulatiove do to pn
                        sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive = &
                            sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive .or. params(i)%enabled
                    end if
                end do
            case(REALM_DIH_G)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%g(params(i)%pn) = params(i)%value
                        ! this must be accumulatiove do to pn
                        sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive = &
                            sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive .or. params(i)%enabled
                    end if
                end do
            case(REALM_DIH_SCEE)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%inv_scee = 1.0d0/params(i)%value
                        ! this must be accumulatiove do to pn
                        sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive = &
                            sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive .or. params(i)%enabled
                    end if
                end do
            case(REALM_DIH_SCNB)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%dihedral_types(params(i)%ids(j))%inv_scnb = 1.0d0/params(i)%value
                        ! this must be accumulatiove do to pn
                        sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive = &
                            sets(j)%top%dihedral_types(params(i)%ids(j))%ffoptactive .or. params(i)%enabled
                    end if
                end do
            case(REALM_IMPR_V)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%improper_types(params(i)%ids(j))%v = params(i)%value
                        sets(j)%top%improper_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_IMPR_G)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%improper_types(params(i)%ids(j))%g = params(i)%value
                        sets(j)%top%improper_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do

        ! NB parameters
            case(REALM_VDW_EPS)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%nb_types(params(i)%ids(j))%eps = params(i)%value
                        sets(j)%top%nb_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_VDW_R0)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%nb_types(params(i)%ids(j))%r0 = params(i)%value
                        sets(j)%top%nb_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_VDW_PA)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%atom_types(params(i)%ids(j))%PA = params(i)%value
                        sets(j)%top%atom_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_VDW_PB)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%atom_types(params(i)%ids(j))%PB = params(i)%value
                        sets(j)%top%atom_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do
            case(REALM_VDW_B0)
                atomicdata_b0opt(params(i)%pn) = params(i)%value
            case(REALM_VDW_RC)
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        sets(j)%top%atom_types(params(i)%ids(j))%RC = params(i)%value
                        sets(j)%top%atom_types(params(i)%ids(j))%ffoptactive = params(i)%enabled
                    end if
                end do

        ! single parameters
            case(REALM_ELE_SQ)
                ele_qscale = params(i)%value

            case(REALM_DAMP_FA)
                damp_fa = params(i)%value
            case(REALM_DAMP_FB)
                damp_fb = params(i)%value
            case(REALM_DAMP_PB)
                damp_pb = params(i)%value
            case(REALM_DAMP_TB)
                damp_tb = params(i)%value
            case(REALM_DAMP_PE)
                damp_pe = params(i)%value

            case(REALM_DISP_S6)
                disp_s6 = params(i)%value
            case(REALM_DISP_S8)
                disp_s8 = params(i)%value
            case(REALM_DISP_S10)
                disp_s10 = params(i)%value

            case(REALM_GLB_SCEE)
                if( params(i)%value .ne. 0.0d0 ) then
                    glb_iscee = 1.0d0 / params(i)%value
                else
                    glb_iscee = 0.0d0
                end if
            case(REALM_GLB_SCNB)
                if( params(i)%value .ne. 0.0d0 ) then
                    glb_iscnb = 1.0d0 / params(i)%value
                else
                    glb_iscnb = 0.0d0
                end if

            case(REALM_K_EXC)
                k_exc = params(i)%value

            case(REALM_K_IND)
                k_ind = params(i)%value

            case(REALM_ZEFF)
                do j=1,ntypes
                    if( j .eq. params(i)%ti ) then
                        types(j)%Zeff = params(i)%value
                    end if
                end do

        ! PAC - partial atomic charges
            case(REALM_PAC)
                ! update independent charges in all sets
                do j=1,nsets
                    if( params(i)%ids(j) .ne. 0 ) then
                        do k=1,sets(j)%top%natoms
                            if( sets(j)%top%atoms(k)%symmclass .eq. params(i)%ids(j) ) then
                                sets(j)%top%atoms(k)%charge = params(i)%value
                                sets(j)%top%atoms(k)%ffoptactive = params(i)%enabled
                            end if
                        end do
                    end if
                end do

            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_to_tops!')
        end select
    end do

    if( PACAsPrms ) then
        ! update dependent charges in sets without references
        do i=1,nsets
            if( sets(i)%nrefs .gt. 0 ) cycle
            chrg = 0.0
            ndep = 0
            do k=1,sets(i)%top%natoms
                if( sets(i)%top%atoms(k)%symmclass .ne. sets(i)%top%nsymm_classes ) then
                    ! independent
                    chrg = chrg + sets(i)%top%atoms(k)%charge
                else
                    ! dependent
                    ndep = ndep + 1
                end if
            end do
            if( ndep .le. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1, &
                     'Consistency violation - no dependent charge was found in ffdev_parameters_to_tops!')
            end if
            chrg = (REAL(sets(i)%top%total_charge,DEVDP) - chrg)/REAL(ndep,DEVDP)
            do k=1,sets(i)%top%natoms
                if( sets(i)%top%atoms(k)%symmclass .eq. sets(i)%top%nsymm_classes ) then
                    ! dependent
                    sets(i)%top%atoms(k)%charge = chrg
                end if
            end do
        end do
        ! update dependent/independent charges in sets with references
        do i=1,nsets
            if( sets(i)%nrefs .le. 0 ) cycle
            ai = 1
            chrg = 0.0d0
            do j=1,sets(i)%nrefs
                refid = sets(i)%refs(j)
                do k=1,sets(refid)%top%natoms
                    sets(i)%top%atoms(ai)%charge = sets(refid)%top%atoms(k)%charge
                    chrg = chrg + sets(refid)%top%atoms(k)%charge
                    ai = ai + 1
                end do
            end do
            if( abs(chrg - REAL(sets(i)%top%total_charge,DEVDP)) .gt. 0.01 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Consistency violation - wrong total charge in ffdev_parameters_to_tops!')
            end if
        end do
    end if

    ! update charge stat
    call ffdev_parameters_update_charge_stat()

    ! mark NB params for update
    do i=1,nsets
        sets(i)%top%nb_params_update = .true.
    end do

end subroutine ffdev_parameters_to_tops

! ==============================================================================
! subroutine ffdev_params_reset_ranges
! ==============================================================================

subroutine ffdev_params_reset_ranges

    use ffdev_parameters_dat

    implicit none
    ! --------------------------------------------------------------------------

! === [ranges] =================================================================

     MinBondD0    =       0.5d0
     MaxBondD0    =       5.0d0
     MinBondK     =       0.0
     MaxBondK     =    1500.0d0
     MinAngleA0   =       0.0
     MaxAngleA0   =       DEV_PI
     MinAngleK    =       0.0d0
     MaxAngleK    =    1000.0d0
     MinDihV      =       0.0d0
     MaxDihV      =      50.0d0
     MinDihG      =       0.0d0
     MaxDihG      =     2*DEV_PI
     MinDihSCEE   =       0.5d0
     MaxDihSCEE   =       3.0d0
     MinDihSCNB   =       0.5d0
     MaxDihSCNB   =       3.0d0
     MinImprV     =       0.0d0
     MaxImprV     =      50.0d0
     MinImprG     =      -DEV_PI
     MaxImprG     =       DEV_PI
     MinDihC      =     -50.0d0
     MaxDihC      =      50.0d0

! non-bonded ERA
     MinVdwEps    =       0.0d0
     MaxVdwEps    =       2.0d0
     MinVdwR0     =       0.5d0
     MaxVdwR0     =       6.0d0

! non-bonded ABC
     MinVdwPA     =      0.0d0
     MaxVdwPA     =     20.0d0

     MinVdwPB     =      1.0d0
     MaxVdwPB     =      6.0d0

     MinVdwB0     =      1.0d0
     MaxVdwB0     =      6.0d0

     MinVdwRC     =      0.0d0
     MaxVdwRC     =      5.0d0

! non-bonded scaling factors
     MinEleSQ     =      0.0d0
     MaxEleSQ     =      3.0d0

! partial atomic charges
     MinPAC       =     -2.0d0
     MaxPAC       =      2.0d0

 ! zeff has special margins
     MinZeff       =     1.0d0
     MaxZeff       =    89.0d0

! vdW interactions
     MinDampFA    =      0.0d0
     MaxDampFA    =      2.0d0
     MinDampFB    =      0.0d0
     MaxDampFB    =      6.0d0

     MinDampPB    =      0.5d0
     MaxDampPB    =      1.5d0
     MinDampTB    =      0.5d0
     MaxDampTB    =      1.5d0
     MinDampPE    =      0.5d0
     MaxDampPE    =      1.5d0

! dispersion scaling
     MinDispS6    =      0.0d0
     MaxDispS6    =     10.0d0
     MinDispS8    =      0.0d0
     MaxDispS8    =     10.0d0
     MinDispS10   =      0.0d0
     MaxDispS10   =     10.0d0

     MinGlbSCEE   =      0.0d0
     MaxGlbSCEE   =      4.0d0
     MinGlbSCNB   =      0.0d0
     MaxGlbSCNB   =      4.0d0

! exchange and induction factors
     MinKExc      =     -15.0d0
     MaxKExc      =      15.0d0
     MinKInd      =     -15.0d0
     MaxKInd      =      15.0d0

end subroutine ffdev_params_reset_ranges

! ==============================================================================
! subroutine ffdev_params_print_ranges
! ==============================================================================

subroutine ffdev_params_print_ranges

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_utils

    implicit none
    integer         :: realmid, pid
    real(DEVDP)     :: sc, lmin, lmax
    logical         :: lprint
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)

    do realmid=REALM_FIRST,REALM_LAST
        lprint = .false.
        do pid=1,nparams
            if( .not. params(pid)%enabled ) cycle
            if( params(pid)%realm .ne. realmid ) cycle
            lprint = .true.
            exit
        end do
        if( lprint ) then
            sc = ffdev_parameters_get_realm_scaling(realmid)
            lmin = ffdev_params_get_lower_bound(realmid)
            lmax = ffdev_params_get_upper_bound(realmid)
            write(DEV_OUT,120) ffdev_parameters_get_realm_name(realmid), lmin*sc, lmax*sc
        end if
    end do

100 format('#   Realm      MinValue   MaxValue')
110 format('# ---------- ---------- ----------')
120 format(A12,1X,F10.2,1X,F10.2)

end subroutine ffdev_params_print_ranges

! ==============================================================================
! subroutine ffdev_params_get_lower_bounds
! ==============================================================================

subroutine ffdev_params_get_lower_bounds(tmpx)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_utils
    use ffdev_atomicdata

    implicit none
    real(DEVDP)     :: tmpx(:)
    ! --------------------------------------------
    integer         :: i, id
    ! --------------------------------------------------------------------------

    id = 0
    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        id = id + 1
        if( params(i)%realm .ne. REALM_ZEFF ) then
            tmpx(id) = ffdev_params_get_lower_bound(params(i)%realm)
        else
            tmpx(id) = ffdev_atomicdata_get_min_Zeff(params(i)%ti)
        end if
    end do
    if( id .ne. nactparms ) stop ! safety fuse

end subroutine ffdev_params_get_lower_bounds

! ==============================================================================
! function ffdev_params_get_lower_bound_for_prms
! ==============================================================================

real(DEVDP) function ffdev_params_get_lower_bound_for_prms(idx)

    use ffdev_parameters_dat
    use ffdev_atomicdata

    implicit none
    integer         :: idx
    ! --------------------------------------------------------------------------

    if( params(idx)%realm .ne. REALM_ZEFF ) then
        ffdev_params_get_lower_bound_for_prms = ffdev_params_get_lower_bound(params(idx)%realm)
    else
        ffdev_params_get_lower_bound_for_prms = ffdev_atomicdata_get_min_Zeff(params(idx)%ti)
    end if

end function ffdev_params_get_lower_bound_for_prms

! ==============================================================================
! subroutine ffdev_params_get_lower_bound
! ==============================================================================

real(DEVDP) function ffdev_params_get_lower_bound(realm)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: realm
    ! --------------------------------------------------------------------------

    select case(realm)
        case(REALM_BOND_D0)
            ffdev_params_get_lower_bound = MinBondD0
        case(REALM_BOND_K)
            ffdev_params_get_lower_bound = MinBondK
        case(REALM_ANGLE_A0)
            ffdev_params_get_lower_bound = MinAngleA0
        case(REALM_ANGLE_K)
            ffdev_params_get_lower_bound = MinAngleK
        case(REALM_DIH_V)
            ffdev_params_get_lower_bound = MinDihV
        case(REALM_DIH_C)
            ffdev_params_get_lower_bound = MinDihC
        case(REALM_DIH_G)
            ffdev_params_get_lower_bound = MinDihG
        case(REALM_DIH_SCEE)
            ffdev_params_get_lower_bound = MinDihSCEE
        case(REALM_DIH_SCNB)
            ffdev_params_get_lower_bound = MinDihSCNB
        case(REALM_IMPR_V)
            ffdev_params_get_lower_bound = MinImprV
        case(REALM_IMPR_G)
            ffdev_params_get_lower_bound = MinImprG

        case(REALM_VDW_EPS)
            ffdev_params_get_lower_bound = MinVdwEps
        case(REALM_VDW_R0)
            ffdev_params_get_lower_bound = MinVdwR0

        case(REALM_VDW_PA)
            ffdev_params_get_lower_bound = MinVdwPA
        case(REALM_VDW_PB)
            ffdev_params_get_lower_bound = MinVdwPB
        case(REALM_VDW_B0)
            ffdev_params_get_lower_bound = MinVdwB0
        case(REALM_VDW_RC)
            ffdev_params_get_lower_bound = MinVdwRC

        case(REALM_ELE_SQ)
            ffdev_params_get_lower_bound = MinEleSQ

        case(REALM_DAMP_FA)
            ffdev_params_get_lower_bound = MinDampFA
        case(REALM_DAMP_FB)
            ffdev_params_get_lower_bound = MinDampFB
        case(REALM_DAMP_PB)
            ffdev_params_get_lower_bound = MinDampPB
        case(REALM_DAMP_TB)
            ffdev_params_get_lower_bound = MinDampTB
        case(REALM_DAMP_PE)
            ffdev_params_get_lower_bound = MinDampPE

        case(REALM_DISP_S6)
            ffdev_params_get_lower_bound = MinDispS6
        case(REALM_DISP_S8)
            ffdev_params_get_lower_bound = MinDispS8
        case(REALM_DISP_S10)
            ffdev_params_get_lower_bound = MinDispS10

        case(REALM_GLB_SCEE)
            ffdev_params_get_lower_bound = MinGlbSCEE
        case(REALM_GLB_SCNB)
            ffdev_params_get_lower_bound = MinGlbSCNB

        case(REALM_PAC)
            ffdev_params_get_lower_bound = MinPAC

        case(REALM_ZEFF)
            ffdev_params_get_lower_bound = MinZeff

        case(REALM_K_EXC)
            ffdev_params_get_lower_bound = MinKExc

        case(REALM_K_IND)
            ffdev_params_get_lower_bound = MinKInd

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_params_get_lower_bounds')
    end select

end function ffdev_params_get_lower_bound

! ==============================================================================
! subroutine ffdev_params_set_lower_bound
! ==============================================================================

subroutine ffdev_params_set_lower_bound(realm,mvalue)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: realm
    real(DEVDP) :: mvalue
    ! --------------------------------------------------------------------------

    select case(realm)
        case(REALM_BOND_D0)
            MinBondD0 = mvalue
        case(REALM_BOND_K)
            MinBondK = mvalue
        case(REALM_ANGLE_A0)
            MinAngleA0 = mvalue
        case(REALM_ANGLE_K)
            MinAngleK = mvalue
        case(REALM_DIH_V)
            MinDihV = mvalue
        case(REALM_DIH_C)
            MinDihC = mvalue
        case(REALM_DIH_G)
            MinDihG = mvalue
        case(REALM_DIH_SCEE)
            MinDihSCEE = mvalue
        case(REALM_DIH_SCNB)
            MinDihSCNB = mvalue
        case(REALM_IMPR_V)
            MinImprV = mvalue
        case(REALM_IMPR_G)
            MinImprG = mvalue

        case(REALM_VDW_EPS)
            MinVdwEps = mvalue
        case(REALM_VDW_R0)
            MinVdwR0 = mvalue

        case(REALM_VDW_PA)
            MinVdwPA = mvalue
        case(REALM_VDW_PB)
            MinVdwPB = mvalue
        case(REALM_VDW_B0)
            MinVdwB0 = mvalue
        case(REALM_VDW_RC)
            MinVdwRC = mvalue

        case(REALM_ELE_SQ)
            MinEleSQ = mvalue

        case(REALM_DAMP_FA)
            MinDampFA = mvalue
        case(REALM_DAMP_FB)
            MinDampFB = mvalue
        case(REALM_DAMP_PB)
            MinDampPB = mvalue
        case(REALM_DAMP_TB)
            MinDampTB = mvalue
        case(REALM_DAMP_PE)
            MinDampPE = mvalue

        case(REALM_DISP_S6)
            MinDispS6 = mvalue
        case(REALM_DISP_S8)
            MinDispS8 = mvalue
        case(REALM_DISP_S10)
            MinDispS10 = mvalue

        case(REALM_GLB_SCEE)
            MinGlbSCEE = mvalue
        case(REALM_GLB_SCNB)
            MinGlbSCNB = mvalue

        case(REALM_PAC)
            MinPAC = mvalue
        case(REALM_ZEFF)
            MinZeff = mvalue

        case(REALM_K_EXC)
            MinKExc = mvalue

        case(REALM_K_IND)
            MinKInd = mvalue

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_params_set_lower_bound')
    end select

end subroutine ffdev_params_set_lower_bound

! ==============================================================================
! subroutine ffdev_params_get_upper_bounds
! ==============================================================================

subroutine ffdev_params_get_upper_bounds(tmpx)

    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_utils

    implicit none
    real(DEVDP)     :: tmpx(:)
    ! --------------------------------------------
    integer         :: i, id
    ! --------------------------------------------------------------------------

    id = 0
    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        id = id + 1
        if( params(i)%realm .ne. REALM_ZEFF ) then
            tmpx(id) = ffdev_params_get_upper_bound(params(i)%realm)
        else
            ! upper limit for effective nuclear charge is the number of protons
            tmpx(id) = types(params(i)%ti)%z
        end if
    end do
    if( id .ne. nactparms ) stop ! safety fuse

end subroutine ffdev_params_get_upper_bounds

! ==============================================================================
! function ffdev_params_get_upper_bound_for_prms
! ==============================================================================

real(DEVDP) function ffdev_params_get_upper_bound_for_prms(idx)

    use ffdev_parameters_dat
    use ffdev_atomicdata

    implicit none
    integer         :: idx
    ! --------------------------------------------------------------------------

    if( params(idx)%realm .ne. REALM_ZEFF ) then
        ffdev_params_get_upper_bound_for_prms = ffdev_params_get_upper_bound(params(idx)%realm)
    else
        ffdev_params_get_upper_bound_for_prms =  types(params(idx)%ti)%z
    end if

end function ffdev_params_get_upper_bound_for_prms

! ==============================================================================
! subroutine ffdev_params_get_upper_bound
! ==============================================================================

real(DEVDP) function ffdev_params_get_upper_bound(realm)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: realm
    ! --------------------------------------------------------------------------

    select case(realm)
        case(REALM_BOND_D0)
            ffdev_params_get_upper_bound = MaxBondD0
        case(REALM_BOND_K)
            ffdev_params_get_upper_bound = MaxBondK
        case(REALM_ANGLE_A0)
            ffdev_params_get_upper_bound = MaxAngleA0
        case(REALM_ANGLE_K)
            ffdev_params_get_upper_bound = MaxAngleK
        case(REALM_DIH_V)
            ffdev_params_get_upper_bound = MaxDihV
        case(REALM_DIH_C)
            ffdev_params_get_upper_bound = MaxDihC
        case(REALM_DIH_G)
            ffdev_params_get_upper_bound = MaxDihG
        case(REALM_DIH_SCEE)
            ffdev_params_get_upper_bound = MaxDihSCEE
        case(REALM_DIH_SCNB)
            ffdev_params_get_upper_bound = MaxDihSCNB
        case(REALM_IMPR_V)
            ffdev_params_get_upper_bound = MaxImprV
        case(REALM_IMPR_G)
            ffdev_params_get_upper_bound = MaxImprG

        case(REALM_VDW_EPS)
            ffdev_params_get_upper_bound = MaxVdwEps
        case(REALM_VDW_R0)
            ffdev_params_get_upper_bound = MaxVdwR0

        case(REALM_VDW_PA)
            ffdev_params_get_upper_bound = MaxVdwPA
        case(REALM_VDW_PB)
            ffdev_params_get_upper_bound = MaxVdwPB
        case(REALM_VDW_B0)
            ffdev_params_get_upper_bound = MaxVdwB0
        case(REALM_VDW_RC)
            ffdev_params_get_upper_bound = MaxVdwRC

        case(REALM_ELE_SQ)
            ffdev_params_get_upper_bound = MaxEleSQ

        case(REALM_DAMP_FA)
            ffdev_params_get_upper_bound = MaxDampFA
        case(REALM_DAMP_FB)
            ffdev_params_get_upper_bound = MaxDampFB
        case(REALM_DAMP_PB)
            ffdev_params_get_upper_bound = MaxDampPB
        case(REALM_DAMP_TB)
            ffdev_params_get_upper_bound = MaxDampTB
        case(REALM_DAMP_PE)
            ffdev_params_get_upper_bound = MaxDampPE

        case(REALM_DISP_S6)
            ffdev_params_get_upper_bound = MaxDispS6
        case(REALM_DISP_S8)
            ffdev_params_get_upper_bound = MaxDispS8
        case(REALM_DISP_S10)
            ffdev_params_get_upper_bound = MaxDispS10

        case(REALM_GLB_SCEE)
            ffdev_params_get_upper_bound = MaxGlbSCEE
        case(REALM_GLB_SCNB)
            ffdev_params_get_upper_bound = MaxGlbSCNB

        case(REALM_PAC)
            ffdev_params_get_upper_bound = MaxPAC
        case(REALM_ZEFF)
            ffdev_params_get_upper_bound = MaxZeff

        case(REALM_K_EXC)
            ffdev_params_get_upper_bound = MaxKExc
        case(REALM_K_IND)
            ffdev_params_get_upper_bound = MaxKInd

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_params_get_upper_bounds')
    end select

end function ffdev_params_get_upper_bound

! ==============================================================================
! subroutine ffdev_params_set_upper_bound
! ==============================================================================

subroutine ffdev_params_set_upper_bound(realm,mvalue)

    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: realm
    real(DEVDP) :: mvalue
    ! --------------------------------------------------------------------------

    select case(realm)
        case(REALM_BOND_D0)
            MaxBondD0 = mvalue
        case(REALM_BOND_K)
            MaxBondK = mvalue
        case(REALM_ANGLE_A0)
            MaxAngleA0 = mvalue
        case(REALM_ANGLE_K)
            MaxAngleK = mvalue
        case(REALM_DIH_V)
            MaxDihV = mvalue
        case(REALM_DIH_C)
            MaxDihC = mvalue
        case(REALM_DIH_G)
            MaxDihG = mvalue
        case(REALM_DIH_SCEE)
            MaxDihSCEE = mvalue
        case(REALM_DIH_SCNB)
            MaxDihSCNB = mvalue
        case(REALM_IMPR_V)
            MaxImprV = mvalue
        case(REALM_IMPR_G)
            MaxImprG = mvalue

        case(REALM_VDW_EPS)
            MaxVdwEps = mvalue
        case(REALM_VDW_R0)
            MaxVdwR0 = mvalue

        case(REALM_VDW_PA)
            MaxVdwPA = mvalue
        case(REALM_VDW_PB)
            MaxVdwPB = mvalue
        case(REALM_VDW_B0)
            MaxVdwB0 = mvalue
        case(REALM_VDW_RC)
            MaxVdwRC = mvalue

        case(REALM_ELE_SQ)
            MaxEleSQ = mvalue

        case(REALM_DAMP_FA)
            MaxDampFA = mvalue
        case(REALM_DAMP_FB)
            MaxDampFB = mvalue
        case(REALM_DAMP_PB)
            MaxDampPB = mvalue
        case(REALM_DAMP_TB)
            MaxDampTB = mvalue
        case(REALM_DAMP_PE)
            MaxDampPE = mvalue

        case(REALM_DISP_S6)
            MaxDispS6 = mvalue
        case(REALM_DISP_S8)
            MaxDispS8 = mvalue
        case(REALM_DISP_S10)
            MaxDispS10 = mvalue

        case(REALM_GLB_SCEE)
            MaxGlbSCEE = mvalue
        case(REALM_GLB_SCNB)
            MaxGlbSCNB = mvalue

        case(REALM_PAC)
            MaxPAC = mvalue
        case(REALM_ZEFF)
            MaxZeff = mvalue

        case(REALM_K_EXC)
            MaxKExc = mvalue
        case(REALM_K_IND)
            MaxKInd = mvalue

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_params_set_upper_bound')
    end select

end subroutine ffdev_params_set_upper_bound

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
                v = ffdev_geometry_get_length( sets(i)%geo(j)%crd, ai, aj )
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
                v = ffdev_geometry_get_angle( sets(i)%geo(j)%crd, ai, aj, ak )
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

! ==============================================================================
! subroutine ffdev_parameters_error
! ==============================================================================

subroutine ffdev_parameters_error(prms,error,grads)

    use ffdev_parameters_dat
    use ffdev_errors_dat

    implicit none
    real(DEVDP)         :: prms(:)
    type(FFERROR_TYPE)  :: error
    real(DEVDP)         :: grads(:)
    ! --------------------------------------------------------------------------

    error%total = 0.0d0
    grads(:) = 0.0d0

    call ffdev_parameters_error_num(prms,error,grads)

end subroutine ffdev_parameters_error

! ==============================================================================
! subroutine ffdev_parameters_error_num
! ==============================================================================

subroutine ffdev_parameters_error_num(prms,error,grads)

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_errors_dat

    implicit none
    real(DEVDP)         :: prms(:)
    type(FFERROR_TYPE)  :: error
    real(DEVDP)         :: grads(:)
    ! --------------------------------------------------------------------------
    real(DEVDP),allocatable :: tmp_prms(:)
    real(DEVDP)             :: d
    type(FFERROR_TYPE)      :: err1,err2
    integer                 :: i
    ! --------------------------------------------------------------------------

    d = 0.5d-5  ! differentiation parameter

    ! calculate base energy
    call ffdev_parameters_error_only(prms,error,.true.)

    ! write(*,*) 'total= ',error%total,prms

    ! allocate temporary geometry object
    allocate( tmp_prms(nactparms) )

    tmp_prms(:) = prms(:)

    ! gradient by numerical differentiation
    do i=1,nactparms
        ! left
        tmp_prms(i) = prms(i) + d
        call ffdev_parameters_error_only(tmp_prms,err1,.true.)

        ! write(*,*) ene1%total,tmp_prms

        ! right
        tmp_prms(i) = prms(i) - d
        call ffdev_parameters_error_only(tmp_prms,err2,.true.)

        ! write(*,*) ene2%total,tmp_prms
        ! write(*,*) ene1%total,ene2%total

        ! gradient
        grads(i) = 0.5d0*(err1%total-err2%total)/d

        ! move back
        tmp_prms(i) = prms(i)
    end do

    ! write(*,*) grads
    ! stop

    ! release temporary geometry object
    deallocate(tmp_prms)

end subroutine ffdev_parameters_error_num

! ==============================================================================
! subroutine ffdev_parameters_error_only
! ==============================================================================

subroutine ffdev_parameters_error_only(prms,error,opterr)

    use ffdev_errors_dat
    use ffdev_errors
    use ffdev_targetset
    use ffdev_xdm

    implicit none
    real(DEVDP)         :: prms(:)
    type(FFERROR_TYPE)  :: error
    logical             :: opterr
    ! --------------------------------------------------------------------------

    ! scatter parameters
    call ffdev_parameters_scatter(prms)

    ! distribute parameters to topologies
    call ffdev_parameters_to_tops()

    ! activate domains by enabled error or summary
    call ffdev_errors_error_setup_domains(opterr)

    ! optimize geometry if requested and then calculate energy, and optionally gradients and hessians
    ! apply combining rules for individual topologies
    call ffdev_targetset_calc_all()

    ! calculate error
    call ffdev_errors_error_only(error)

end subroutine ffdev_parameters_error_only

! ------------------------------------------------------------------------------

end module ffdev_parameters

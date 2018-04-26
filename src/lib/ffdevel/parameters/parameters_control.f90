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

module ffdev_parameters_control

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_parameters_ctrl_identities
! ==============================================================================

subroutine ffdev_parameters_ctrl_identities(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    integer                     :: niden,i
    character(PRMFILE_MAX_PATH) :: string
    logical                     :: rst
    character(50)               :: key
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    ! clear all identities
    do i=1,nparams
        params(i)%identity = 0
    end do
    write(DEV_OUT,25)

    if( .not. prmfile_open_section(fin,'identities') ) then
        write(DEV_OUT,20)
        return
    end if

    niden = 0
    rst = prmfile_first_line(fin)
    do while ( prmfile_get_line(fin,string) )
        read(string,*) key
        select case(key)
            case('dih_v')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'byrot'
                call setup_dih_identity_byrot(REALM_DIH_V,niden)
            case('dih_gamma')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'byrot'
                call setup_dih_identity_byrot(REALM_DIH_G,niden)
            case('dih_c')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'byrot'
                call setup_dih_identity_byrot(REALM_DIH_C,niden)
            case('dih_scee')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'all'
                call setup_realm_identity(REALM_DIH_SCEE,niden)
            case('dih_scnb')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'all'
                call setup_realm_identity(REALM_DIH_SCNB,niden)
            case('impr_v')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'all'
                call setup_realm_identity(REALM_IMPR_V,niden)
            case('impr_gamma')
                write(DEV_OUT,*)
                write(DEV_OUT,40) adjustl(key),'all'
                call setup_realm_identity(REALM_IMPR_G,niden)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported identity key '''//trim(key)//'''!')
        end select
    end do

    if( niden .le. 0 ) then
        write(DEV_OUT,20)
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,30) niden


10 format('=== [identities] ===============================================================')
20 format('>> INFO: No identities defined!')
25 format('clear                                           all              (initial setup)')
30 format('Number of defined identities = ', I6)
40 format(A50,1X,A)

end subroutine ffdev_parameters_ctrl_identities

! ------------------------------------------------------------------------------

subroutine setup_realm_identity(realm,niden)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer         :: realm
    integer         :: niden
    ! --------------------------------------------
    integer         :: i,j,lniden
    real(DEVDP)     :: value
    ! --------------------------------------------------------------------------

    if( (realm .eq. REALM_DIH_V) .or. (realm .eq. REALM_DIH_G) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unsupported real for setup_realm_identity!')
    end if

    lniden = 0
    do i=1,nparams
        ! skip different realms
        if( params(i)%realm .ne. realm ) cycle

        ! get root parameters
        write(DEV_OUT,10) i
        value = params(i)%value

        ! set identity for the rest
        do j=i+1,nparams
            if( params(j)%realm .ne. realm ) cycle
            params(j)%identity = i
            params(j)%enabled = .false.
            params(j)%value = value
            lniden = lniden + 1
            write(DEV_OUT,20,ADVANCE='NO') j
            if( MOD(lniden,16) .eq. 0 ) write(DEV_OUT,*)
        end do
        if( MOD(lniden,16) .ne. 0 ) write(DEV_OUT,*)
        exit
    end do

    niden = niden + lniden

 10 format('Parameter ',I4,' is root for:')
 20 format(I4,1X)

end subroutine setup_realm_identity

! ------------------------------------------------------------------------------

subroutine setup_dih_identity_byrot(realm,niden)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer         :: realm
    integer         :: niden
    ! --------------------------------------------
    integer         :: i,j,lniden,pn
    real(DEVDP)     :: value
    ! --------------------------------------------------------------------------

    if( (realm .ne. REALM_DIH_V) .and. (realm .ne. REALM_DIH_G) .and. (realm .ne. REALM_DIH_C) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unsupported realm for setup_dih_identity_byrot!')
    end if

    do i=1,nparams
        ! skip different realms
        if( params(i)%realm .ne. realm ) cycle
        ! skip already set identities
        if( params(i)%identity .gt. 0 ) cycle

        ! get root parameters
        write(DEV_OUT,10) i
        value = params(i)%value
        pn = params(i)%pn

        ! set identity for the compatible parameters
        lniden = 0
        do j=i+1,nparams
            ! skip different realms
            if( params(j)%realm .ne. realm ) cycle
            ! skip incorrect pn
            if( params(j)%pn .ne. pn ) cycle

            ! skip incompatible middle bond types
            if( .not. ( ((params(i)%tj .eq. params(j)%tj) .and. (params(i)%tk .eq. params(j)%tk)) &
                        .or. ((params(i)%tj .eq. params(j)%tk) .and. (params(i)%tk .eq. params(j)%tj)) ) ) cycle

            params(j)%identity = i
            params(j)%enabled = .false.
            params(j)%value = value

            lniden = lniden + 1
            write(DEV_OUT,20,ADVANCE='NO') j
            if( MOD(lniden,16) .eq. 0 ) write(DEV_OUT,*)
        end do
        if( (lniden .gt. 0 ) .and. (MOD(lniden,16) .ne. 0) ) write(DEV_OUT,*)
        niden = niden + lniden

    end do


 10 format('Parameter ',I4,' is root for:')
 20 format(I4,1X)

end subroutine setup_dih_identity_byrot

! ==============================================================================
! subroutine ffdev_parameters_ctrl_realms
! ==============================================================================

subroutine ffdev_parameters_ctrl_realms(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    integer                     :: i,nchanged
    character(PRMFILE_MAX_PATH) :: string,realm
    logical                     :: rst
    character(10)               :: key
    ! --------------------------------------------------------------------------

    ! disable all parameter realms
    call ffdev_parameters_disable_all_realms()

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,30) nparams,'disable all (initial setup)'

    if( .not. prmfile_open_section(fin,'parameters') ) then
        call ffdev_utils_exit(DEV_OUT,1,'No parameters to optimize! No [parameters] section found!')
    end if

    rst = prmfile_first_line(fin)
    do while ( prmfile_get_line(fin,string) )
        read(string,*) key, realm
        select case(key)
            case('enable')
                call change_realms(realm,.true.,string,nchanged)
                write(DEV_OUT,30) nchanged,trim(string)
            case('disable')
                call change_realms(realm,.false.,string,nchanged)
                write(DEV_OUT,30) nchanged,trim(string)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported action key '''//trim(key)//'''!')
        end select
    end do

    ! count active parameters
    nactparms = 0
    do i=1,nparams
        if( params(i)%enabled ) nactparms = nactparms + 1
    end do

    write(DEV_OUT,*)
    write(DEV_OUT,40) nactparms

    if( nactparms .le. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'No parameters to optimize!')
    end if

10 format('=== [parameters] ===============================================================')
30 format('Altered parameters = ',I3,' | ', A)
40 format('Number of active parameters = ',I6)

end subroutine ffdev_parameters_ctrl_realms

! ------------------------------------------------------------------------------

logical function is_realm_option(options,key)

    use prmfile

    implicit none
    character(*)    :: options
    character(*)    :: key
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: k1
    character(PRMFILE_MAX_PATH) :: k2
    character(PRMFILE_MAX_PATH) :: k3
    ! --------------------------------------------------------------------------

    is_realm_option = .false.
    read(options,*,end=100,err=100) k1,k2,k3

    is_realm_option = trim(k3) .eq. trim(key)
100 return

end function is_realm_option

! ------------------------------------------------------------------------------

subroutine change_realms(realm,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    character(*)    :: realm
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
    ! --------------------------------------------
    integer         :: realmid, i, pid, ti, tj, tk, tl
    logical         :: lenable
    character(80)   :: k1, k2, k3, sti, stj, stk, stl
    ! --------------------------------------------------------------------------

    pid = 0

    ! decode realm
    select case(realm)
        case('E_offset')
            realmid = REALM_EOFFSET
        case('bond_r0')
            realmid = REALM_BOND_D0
        case('bond_k')
            realmid = REALM_BOND_K

            if( is_realm_option(options,'zero') ) then
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%value .eq. 0.0d0 ) then
                            if( params(i)%identity .eq. 0 ) then
                                params(i)%enabled = enable
                            end if
                        end if
                    end if
                end do
            else
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%identity .eq. 0 ) then
                            params(i)%enabled = enable
                        end if
                    end if
                end do
            end if
            return

        case('angle_a0')
            realmid = REALM_ANGLE_A0
        case('angle_k')
            realmid = REALM_ANGLE_K

            if( is_realm_option(options,'zero') ) then
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%value .eq. 0.0d0 ) then
                            if( params(i)%identity .eq. 0 ) then
                                params(i)%enabled = enable
                            end if
                        end if
                    end if
                end do
            else
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%identity .eq. 0 ) then
                            params(i)%enabled = enable
                        end if
                    end if
                end do
            end if
            return

        case('dih_v')
            realmid = REALM_DIH_V

            if( is_realm_option(options,'bond') ) then
                call change_dih_realm_bond(realmid,enable,options)
                return
            end if

            if( is_realm_option(options,'pn') ) then
                call change_dih_realm_pn(realmid,enable,options)
                return
            end if

        case('dih_c')
            realmid = REALM_DIH_C

            if( is_realm_option(options,'bond') ) then
                call change_dih_realm_bond(realmid,enable,options)
                return
            end if

            if( is_realm_option(options,'pn') ) then
                call change_dih_realm_pn(realmid,enable,options)
                return
            end if

        case('dih_gamma')
            realmid = REALM_DIH_G

            if( is_realm_option(options,'bond') ) then
                call change_dih_realm_bond(realmid,enable,options)
                return
            end if

            if( is_realm_option(options,'pn') ) then
                call change_dih_realm_pn(realmid,enable,options)
                return
            end if

        case('dih_scee')
            realmid = REALM_DIH_SCEE
        case('dih_scnb')
            realmid = REALM_DIH_SCNB
        case('impr_v')
            realmid = REALM_IMPR_V
        case('impr_gamma')
            realmid = REALM_IMPR_G
        case('vdw_eps')
            realmid = REALM_VDW_EPS
        case('vdw_r0')
            realmid = REALM_VDW_R0
        case('vdw_alpha')
            realmid = REALM_VDW_ALPHA
        case('vdw_A')
            realmid = REALM_VDW_A
        case('vdw_B')
            realmid = REALM_VDW_B
        case('vdw_C6')
            realmid = REALM_VDW_C6
        case('vdw_C8')
            realmid = REALM_VDW_C8
        case('all')
            realmid = -1
        case default
            ! is it integer
            read(realm,*,end=333,err=333) pid
            goto 444

333         call ffdev_utils_exit(DEV_OUT,1,'Unsupported realm key '''//trim(realm)//'''!')
    end select

    ti = 0
    tj = 0
    tk = 0
    tl = 0

    ! this option is generic  - filter by types
    if( is_realm_option(options,'types') ) then
        ! read types
        read(options,*,end=555,err=555) k1, k2, k3, sti, stj, stk, stl
555     do i=1,ntypes
            if( types(i)%name .eq. sti ) ti = i
            if( types(i)%name .eq. stj ) tj = i
            if( types(i)%name .eq. stk ) tk = i
            if( types(i)%name .eq. stl ) tl = i
        end do
        ! we need at least two types
        if( (tj .eq. 0) .and. (tk .eq. 0) .and. (tl .eq. 0) ) return
    end if

    nchanged = 0

444 do i=1,nparams

        ! these two options have application only for NB
        ! but for simplicity we will consider them for all parameters
        if( is_realm_option(options,'like') ) then
            if( params(i)%ti .ne. params(i)%tj ) cycle
        end if
        if( is_realm_option(options,'unlike') ) then
            if( params(i)%ti .eq. params(i)%tj ) cycle
        end if
        if( is_realm_option(options,'types') ) then
            if( (tk .eq. 0) .and. (tl .eq. 0) ) then
                ! bond or NB
                if( (params(i)%tk .ne. 0) .or. (params(i)%tl .ne. 0) ) cycle   ! incompatible parameter
                if( (params(i)%ti .eq. 0) .or. (params(i)%tj .eq. 0) ) cycle   ! incompatible parameter
                if( .not. ( ( (ti .eq. params(i)%ti) .and. (tj .eq. params(i)%tj) ) .or. &
                            ( (ti .eq. params(i)%tj) .and. (tj .eq. params(i)%ti) ) ) ) cycle
            else if ( tl .eq. 0 ) then
                ! angle
                if( params(i)%tl .ne. 0 ) continue   ! incompatible parameter
                if( (params(i)%ti .eq. 0) .or. (params(i)%tj .eq. 0) .or. (params(i)%tk .eq. 0)  ) cycle   ! incompatible parameter
                if( .not. ( ( (ti .eq. params(i)%ti) .and. (tj .eq. params(i)%tj) .and. (tk .eq. params(i)%tk) ) .or. &
                            ( (ti .eq. params(i)%tk) .and. (tj .eq. params(i)%tj) .and. (tk .eq. params(i)%ti) ) ) ) cycle
            else
                ! dihedral
                if( (params(i)%ti .eq. 0) .or. (params(i)%tj .eq. 0) .or. &
                    (params(i)%tk .eq. 0) .or. (params(i)%tl .eq. 0)  ) cycle   ! incompatible parameter
                if( .not. ( ( (ti .eq. params(i)%ti) .and. (tj .eq. params(i)%tj) .and. &
                              (tk .eq. params(i)%tk) .and. (tl .eq. params(i)%tl) ) .or. &
                            ( (ti .eq. params(i)%tl) .and. (tj .eq. params(i)%tj) .and. &
                              (tk .eq. params(i)%tk) .and. (tl .eq. params(i)%ti) ) ) ) cycle

            end if
        end if

        if( (realmid .eq. -1) .or. (params(i)%realm .eq. realmid) .or. (pid .eq. i) ) then
            if( params(i)%identity .eq. 0 ) then
                lenable = enable
                select case(LastNBMode)
                    case(NB_MODE_LJ)
                        if( params(i)%realm .eq. REALM_VDW_ALPHA ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_A ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_B ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_C6 ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_C8 ) lenable = .false.
                    case(NB_MODE_EXP6)
                        if( params(i)%realm .eq. REALM_VDW_A ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_B ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_C6 ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_C8 ) lenable = .false.
                    case(NB_MODE_BP)
                        if( params(i)%realm .eq. REALM_VDW_EPS ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_R0 ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_ALPHA ) lenable = .false.
                    case(NB_MODE_EXPONLY)
                        if( params(i)%realm .eq. REALM_VDW_EPS ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_R0 ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_ALPHA ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_C6 ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_C8 ) lenable = .false.
                    case(NB_MODE_MMD3)
                        if( params(i)%realm .eq. REALM_VDW_EPS ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_R0 ) lenable = .false.
                        if( params(i)%realm .eq. REALM_VDW_ALPHA ) lenable = .false.
                end select
                params(i)%enabled = lenable
                nchanged = nchanged + 1
            end if
        end if
    end do

end subroutine change_realms

! ------------------------------------------------------------------------------

subroutine change_dih_realm_bond(realmid,enable,options)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    ! --------------------------------------------
    integer                     :: sid, aj, ak, i, j
    character(PRMFILE_MAX_PATH) :: k1
    character(PRMFILE_MAX_PATH) :: k2
    character(PRMFILE_MAX_PATH) :: k3
    logical                     :: found
    ! --------------------------------------------------------------------------

    read(options,*,end=100,err=100) k1,k2,k3,sid,aj,ak

    if( (sid .lt. 1) .or. (sid .gt. nsets) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Set ID is out-of-legal range!')
    end if

    if( (aj .lt. 1) .or. (aj .gt. sets(sid)%top%natoms) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Atom i is out-of-legal range!')
    end if

    if( (ak .lt. 1) .or. (ak .gt. sets(sid)%top%natoms) ) then
        call ffdev_utils_exit(DEV_OUT,1,'Atom j is out-of-legal range!')
    end if

    do i=1,nparams
        if( params(i)%realm .eq. realmid ) then
            ! is it in the set?
            if( params(i)%ids(sid) .le. 0 ) cycle

            ! is in torsion with given type?
            found = .false.
            do j=1,sets(sid)%top%ndihedrals
                if( ((sets(sid)%top%dihedrals(j)%aj .eq. aj) .and. (sets(sid)%top%dihedrals(j)%ak .eq. ak)) &
                   .or. ((sets(sid)%top%dihedrals(j)%ak .eq. aj) .and. (sets(sid)%top%dihedrals(j)%aj .eq. ak)) ) then
                    if( sets(sid)%top%dihedrals(j)%dt .eq. params(i)%ids(sid) ) then
                        found = .true.
                    end if
                end if
            end do
            if( .not. found ) cycle

            if( params(i)%identity .eq. 0 ) then
                params(i)%enabled = enable
            end if
        end if
    end do

    return

100 call ffdev_utils_exit(DEV_OUT,1,'Illegal option ''bond'' for dih_* realm!')

end subroutine change_dih_realm_bond

! ------------------------------------------------------------------------------

subroutine change_dih_realm_pn(realmid,enable,options)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    ! --------------------------------------------
    integer                     :: pn, i
    character(PRMFILE_MAX_PATH) :: k1
    character(PRMFILE_MAX_PATH) :: k2
    character(PRMFILE_MAX_PATH) :: k3
    ! --------------------------------------------------------------------------

    read(options,*,end=100,err=100) k1,k2,k3,pn

    do i=1,nparams
        if( params(i)%realm .eq. realmid ) then
            if( params(i)%pn .ne. pn ) cycle

            if( params(i)%identity .eq. 0 ) then
                params(i)%enabled = enable
            end if
        end if
    end do

    return

100 call ffdev_utils_exit(DEV_OUT,1,'Illegal option ''pn'' for dih_* realm!')

end subroutine change_dih_realm_pn

! ==============================================================================
! subroutine ffdev_parameters_ctrl_control
! ==============================================================================

subroutine ffdev_parameters_ctrl_control(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(PRMFILE_TYPE)          :: fin
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'control') ) then
        select case(NBParamsMode)
            case(NB_PARAMS_MODE_NORMAL)
                write(DEV_OUT,25) 'normal'
            case(NB_PARAMS_MODE_LIKE_ONLY)
                write(DEV_OUT,25) 'like-only'
            case(NB_PARAMS_MODE_LIKE_ALL)
                write(DEV_OUT,25) 'like-all'
            case(NB_PARAMS_MODE_ALL)
                write(DEV_OUT,25) 'all'
        end select
        write(DEV_OUT,35) prmfile_onoff(NBERAOnly)
        return
    end if

    if( prmfile_get_string_by_key(fin,'nb_params', string)) then
        select case(trim(string))
            case('normal')
                NBParamsMode = NB_PARAMS_MODE_NORMAL
                write(DEV_OUT,20) trim(string)
            case('like-only')
                NBParamsMode = NB_PARAMS_MODE_LIKE_ONLY
                write(DEV_OUT,20) trim(string)
            case('like-all')
                NBParamsMode = NB_PARAMS_MODE_LIKE_ALL
                write(DEV_OUT,20) trim(string)
            case('all')
                NBParamsMode = NB_PARAMS_MODE_ALL
                write(DEV_OUT,20) trim(string)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported nb_params ('//trim(string)//')')
        end select
    else
        select case(NBParamsMode)
            case(NB_PARAMS_MODE_NORMAL)
                write(DEV_OUT,25) 'normal'
            case(NB_PARAMS_MODE_LIKE_ONLY)
                write(DEV_OUT,25) 'like-only'
            case(NB_PARAMS_MODE_LIKE_ALL)
                write(DEV_OUT,25) 'like-all'
            case(NB_PARAMS_MODE_ALL)
                write(DEV_OUT,25) 'all'
        end select
    end if

    if( prmfile_get_logical_by_key(fin,'era_only', NBERAOnly)) then
        write(DEV_OUT,30) prmfile_onoff(NBERAOnly)
    else
        write(DEV_OUT,35) prmfile_onoff(NBERAOnly)
    end if

    if( prmfile_get_real8_by_key(fin,'lj2exp6_alpha', lj2exp6_alpha)) then
        write(DEV_OUT,40) lj2exp6_alpha
    else
        write(DEV_OUT,45) lj2exp6_alpha
    end if

    return

 10 format('=== [control] ==================================================================')

 20  format ('NB parameter assembly mode (nb_params) = ',a12)
 25  format ('NB parameter assembly mode (nb_params) = ',a12,'                  (default)')
 30  format ('Consider only ERA realms (era_only)    = ',a12)
 35  format ('Consider only ERA realms (era_only)    = ',a12,'                  (default)')
 40  format ('Default value of alpha (lj2exp6_alpha) = ',f12.7)
 45  format ('Default value of alpha (lj2exp6_alpha) = ',f12.7,'                  (default)')

end subroutine ffdev_parameters_ctrl_control

! ==============================================================================
! subroutine ffdev_parameters_ctrl_error
! ==============================================================================

subroutine ffdev_parameters_ctrl_error(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'error') ) then         
        write(DEV_OUT,25) prmfile_onoff(EnableEnergyError)
        write(DEV_OUT,35) EnergyErrorWeight
        write(DEV_OUT,45) prmfile_onoff(EnableGradientError)
        write(DEV_OUT,55) GradientErrorWeight
        write(DEV_OUT,65) prmfile_onoff(EnableHessianError)
        write(DEV_OUT,75) HessianErrorWeight
        
        write(DEV_OUT,115) prmfile_onoff(EnableBondError)
        write(DEV_OUT,125) BondErrorWeight
        write(DEV_OUT,215) prmfile_onoff(EnableAngleError)
        write(DEV_OUT,225) AngleErrorWeight
        write(DEV_OUT,315) prmfile_onoff(EnableTorsionError)
        write(DEV_OUT,325) TorsionErrorWeight 
        
        write(DEV_OUT,415) prmfile_onoff(EnableNBDistanceError)
        write(DEV_OUT,425) NBDistanceErrorWeight 
        return
    end if
    
    if( prmfile_get_logical_by_key(fin,'energy', EnableEnergyError)) then
        write(DEV_OUT,20) prmfile_onoff(EnableEnergyError)
    else
        write(DEV_OUT,25) prmfile_onoff(EnableEnergyError)
    end if
    if( prmfile_get_real8_by_key(fin,'energy_weight', EnergyErrorWeight)) then
        write(DEV_OUT,30) EnergyErrorWeight
    else
        write(DEV_OUT,35) EnergyErrorWeight
    end if

    if( prmfile_get_logical_by_key(fin,'gradient', EnableGradientError)) then
        write(DEV_OUT,40) prmfile_onoff(EnableGradientError)
    else
        write(DEV_OUT,45) prmfile_onoff(EnableGradientError)
    end if
    if( prmfile_get_real8_by_key(fin,'gradient_weight', GradientErrorWeight)) then
        write(DEV_OUT,50) GradientErrorWeight
    else
        write(DEV_OUT,55) GradientErrorWeight
    end if

    if( prmfile_get_logical_by_key(fin,'hessian', EnableHessianError)) then
        write(DEV_OUT,60) prmfile_onoff(EnableHessianError)
    else
        write(DEV_OUT,65) prmfile_onoff(EnableHessianError)
    end if
    if( prmfile_get_real8_by_key(fin,'hessian_weight', HessianErrorWeight)) then
        write(DEV_OUT,70) HessianErrorWeight
    else
        write(DEV_OUT,75) HessianErrorWeight
    end if
    
    if( prmfile_get_logical_by_key(fin,'bond', EnableBondError)) then
        write(DEV_OUT,110) prmfile_onoff(EnableBondError)
    else
        write(DEV_OUT,115) prmfile_onoff(EnableBondError)
    end if
    if( prmfile_get_real8_by_key(fin,'bond_weight', BondErrorWeight)) then
        write(DEV_OUT,120) BondErrorWeight
    else
        write(DEV_OUT,125) BondErrorWeight
    end if   
    
    if( prmfile_get_logical_by_key(fin,'angle', EnableAngleError)) then
        write(DEV_OUT,210) prmfile_onoff(EnableAngleError)
    else
        write(DEV_OUT,215) prmfile_onoff(EnableAngleError)
    end if
    if( prmfile_get_real8_by_key(fin,'angle_weight', AngleErrorWeight)) then
        write(DEV_OUT,220) AngleErrorWeight
    else
        write(DEV_OUT,225) AngleErrorWeight
    end if    

    if( prmfile_get_logical_by_key(fin,'torsion', EnableTorsionError)) then
        write(DEV_OUT,310) prmfile_onoff(EnableTorsionError)
    else
        write(DEV_OUT,315) prmfile_onoff(EnableTorsionError)
    end if
    if( prmfile_get_real8_by_key(fin,'torsion_weight', TorsionErrorWeight)) then
        write(DEV_OUT,320) TorsionErrorWeight
    else
        write(DEV_OUT,325) TorsionErrorWeight
    end if  
    
    
    if( prmfile_get_logical_by_key(fin,'nbdist', EnableNBDistanceError)) then
        write(DEV_OUT,410) prmfile_onoff(EnableNBDistanceError)
    else
        write(DEV_OUT,415) prmfile_onoff(EnableNBDistanceError)
    end if
    if( prmfile_get_real8_by_key(fin,'nbdist_weight', NBDistanceErrorWeight)) then
        write(DEV_OUT,420) NBDistanceErrorWeight
    else
        write(DEV_OUT,425) NBDistanceErrorWeight
    end if       
    
 10 format('=== [errors] ===================================================================')
 
 20  format ('Energy error (energy)                  = ',a12)
 25  format ('Energy error (energy)                  = ',a12,'                  (default)')
 30  format ('Energy error weight (energy_weight)    = ',f21.8)
 35  format ('Energy error weight (energy_weight)    = ',f21.8,'         (default)')

 40  format ('Gradient error (gradient)              = ',a12)
 45  format ('Gradient error (gradient)              = ',a12,'                  (default)')
 50  format ('Gradient error weight (gradient_weight)= ',f21.8)
 55  format ('Gradient error weight (gradient_weight)= ',f21.8,'         (default)')

 60  format ('Hessian error (hessian)                = ',a12)
 65  format ('Hessian error (hessian)                = ',a12,'                  (default)')
 70  format ('Hessian error weight (hessian_weight)  = ',f21.8)
 75  format ('Hessian error weight (hessian_weight)  = ',f21.8,'         (default)')

110  format ('Bond error (bond)                      = ',a12)
115  format ('Bond error (bond)                      = ',a12,'                  (default)')
120  format ('Bond error weight (bond_weight)        = ',f21.8)
125  format ('Bond error weight (bond_weight)        = ',f21.8,'         (default)')
 
210  format ('Angle error (angle)                    = ',a12)
215  format ('Angle error (angle)                    = ',a12,'                  (default)')
220  format ('Angle error weight (angle_weight)      = ',f21.8)
225  format ('Angle error weight (angle_weight)      = ',f21.8,'         (default)') 
 
310  format ('Torsion error (torsion)                = ',a12)
315  format ('Torsion error (torsion)                = ',a12,'                  (default)')
320  format ('Torsion error weight (torsion_weight)  = ',f21.8)
325  format ('Torsion error weight (torsion_weight)  = ',f21.8,'         (default)')  

410  format ('NB distance error (nbdist)             = ',a12)
415  format ('NB distance error (nbdist)             = ',a12,'                  (default)')
420  format ('NB distance error wght (nbdist_weight) = ',f21.8)
425  format ('NB distance error wght (nbdist_weight) = ',f21.8,'         (default)')  
 
 
end subroutine ffdev_parameters_ctrl_error

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_parameters_ctrl_files(fin)

    use prmfile
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! -----------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [files] ===================================================================='

    if(.not. prmfile_open_section(fin,'files')) then
        write (DEV_OUT,15) trim(InpParamFileName)
        write (DEV_OUT,25) trim(OutParamFileName)
        write (DEV_OUT,35) trim(OutAmberPrmsFileName)
        return
    end if

    if(.not. prmfile_get_string_by_key(fin,'input', InpParamFileName)) then
        write (DEV_OUT,15) trim(InpParamFileName)
    else
        write (DEV_OUT,10) trim(InpParamFileName)
    end if

    if(.not. prmfile_get_string_by_key(fin,'final', OutParamFileName)) then
        write (DEV_OUT,25) trim(OutParamFileName)
    else
        write (DEV_OUT,20) trim(OutParamFileName)
    end if

    if(.not. prmfile_get_string_by_key(fin,'amber', OutAmberPrmsFileName)) then
        write (DEV_OUT,35) trim(OutAmberPrmsFileName)
    else
        write (DEV_OUT,30) trim(OutAmberPrmsFileName)
    end if

    return

10 format('Input parameters (input)               = ',a)
15 format('Input parameters (input)               = ',a12,'                  (default)')
20 format('Final parameters file (final)          = ',a)
25 format('Final parameters file (final)          = ',a12,'                  (default)')
30 format('Final AMBER parameters (amber)         = ',a)
35 format('Final AMBER parameters (amber)         = ',a12,'                  (default)')

end subroutine ffdev_parameters_ctrl_files

! ==============================================================================
! subroutine ffdev_parameters_ctrl_ffmanip
! ==============================================================================

subroutine ffdev_parameters_ctrl_ffmanip(fin,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: noexec
    ! --------------------------------------------
    logical                     :: rst
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FFMANIP', ':')

    rst = prmfile_first_section(fin)
    do while( rst )

        ! open set section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_OUT,1,'Unable to get section name!')
        end if

        select case(string)
            case('bond_r0')
                call ffdev_parameters_ctrl_bond_r0(fin,noexec)
            case('angle_a0')
                call ffdev_parameters_ctrl_angle_a0(fin,noexec)
            case('nbmanip')
                call ffdev_parameters_ctrl_nbmanip(fin,noexec)
            case('nbload')
                call ffdev_parameters_ctrl_nbload(fin,noexec)
        end select

        rst = prmfile_next_section(fin)
    end do


end subroutine ffdev_parameters_ctrl_ffmanip

! ==============================================================================
! subroutine ffdev_parameters_ctrl_bond_r0
! ==============================================================================

subroutine ffdev_parameters_ctrl_bond_r0(fin,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: noexec
    ! --------------------------------------------
    logical                     :: rst
    integer                     :: mode, typeid, i
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [bond_r0] =================================================================='

    ! read mode if it is present
    rst = prmfile_get_string_by_key(fin,'mode',string)
    if( rst .eqv. .true. ) then
        select case(string)
            case('lowest')
                mode = 0
                write(DEV_OUT,10)
            case('average')
                mode = 1
                write(DEV_OUT,20)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported initff mode '//trim(string)//'!')
        end select
        string = 'average'
    else
        mode = 1
        write(DEV_OUT,25)
    end if

    ! do we have types specified?
    rst = prmfile_get_string_by_key(fin,'type',string)
    if( rst ) then
        rst = prmfile_first_line(fin)
        do while ( rst )
            rst = prmfile_get_integer_by_key(fin,'type',typeid)
            if( rst .eqv. .false. ) then
                call ffdev_utils_exit(DEV_OUT,1,'type must be integer number!')
            end if
            if( (typeid .le. 0) .or. (typeid .gt. nparams) ) then
                call ffdev_utils_exit(DEV_OUT,1,'type is out-of-legal range!')
            end if
            if( params(typeid)%realm .ne. REALM_BOND_D0 ) then
                call ffdev_utils_exit(DEV_OUT,1,'type is not bond_r0!')
            end if
            if( noexec .eqv. .false. ) call ffdev_parameters_bond_r0_init(typeid,mode)
            rst = prmfile_next_line(fin)
        end do
    else
        do i=1,nparams
            if( params(i)%realm .eq. REALM_BOND_D0 ) then
                if( noexec .eqv. .false. ) call ffdev_parameters_bond_r0_init(i,mode)
            end if
        end do
    end if

 10 format(/,'> Init mode = lowest')
 20 format(/,'> Init mode = average')
 25 format(/,'> Init mode = average  [default]')

end subroutine ffdev_parameters_ctrl_bond_r0

! ==============================================================================
! subroutine ffdev_parameters_ctrl_angle_a0
! ==============================================================================

subroutine ffdev_parameters_ctrl_angle_a0(fin,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: noexec
    ! --------------------------------------------
    logical                     :: rst
    integer                     :: mode, typeid, i
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [angle_a0] ================================================================='

    ! read mode if it is present
    rst = prmfile_get_string_by_key(fin,'mode',string)
    if( rst .eqv. .true. ) then
        select case(string)
            case('lowest')
                mode = 0
                write(DEV_OUT,10)
            case('average')
                mode = 1
                write(DEV_OUT,20)
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported initff mode '//trim(string)//'!')
        end select
        string = 'average'
    else
        mode = 1
        write(DEV_OUT,25)
    end if

    ! do we have types specified?
    rst = prmfile_get_string_by_key(fin,'type',string)
    if( rst ) then
        rst = prmfile_first_line(fin)
        do while ( rst )
            rst = prmfile_get_integer_by_key(fin,'type',typeid)
            if( rst .eqv. .false. ) then
                call ffdev_utils_exit(DEV_OUT,1,'type must be integer number!')
            end if
            if( (typeid .le. 0) .or. (typeid .gt. nparams) ) then
                call ffdev_utils_exit(DEV_OUT,1,'type is out-of-legal range!')
            end if
            if( params(typeid)%realm .ne. REALM_ANGLE_A0 ) then
                call ffdev_utils_exit(DEV_OUT,1,'type is not angle_a0!')
            end if
            rst = prmfile_next_line(fin)
        end do
    else
        do i=1,nparams
            if( params(i)%realm .eq. REALM_ANGLE_A0 ) then
                if( noexec .eqv. .false. ) call ffdev_parameters_angle_a0_init(i,mode)
            end if
        end do
    end if

 10 format(/,'> Init mode = lowest')
 20 format(/,'> Init mode = average')
 25 format(/,'> Init mode = average  [default]')

end subroutine ffdev_parameters_ctrl_angle_a0

! ==============================================================================
! subroutine ffdev_parameters_ctrl_nbmanip
! ==============================================================================

subroutine ffdev_parameters_ctrl_nbmanip(fin,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils
    use ffdev_topology_dat

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: noexec
    ! --------------------------------------------
    logical                     :: rst
    character(50)               :: key
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

! programatic NB change (order dependent)
    rst = prmfile_first_line(fin)
    do while( rst )
        rst = prmfile_get_field_on_line(fin,key)

        select case(trim(key))
            case('comb_rules')
                if( prmfile_get_field_on_line(fin,string) ) then
                    call ffdev_parameters_ctrl_nbmanip_comb_rules(string,noexec)
                end if
            case('nb_mode')
                if( prmfile_get_field_on_line(fin,string) ) then
                    call ffdev_parameters_ctrl_nbmanip_nb_mode(string,noexec)
                end if
            case default
                ! do nothing
        end select
        rst = prmfile_next_line(fin)
    end do

10 format('=== [nbmanip] ==================================================================')

end subroutine ffdev_parameters_ctrl_nbmanip

! ==============================================================================
! subroutine ffdev_parameters_ctrl_nbmanip_comb_rules
! ==============================================================================

subroutine ffdev_parameters_ctrl_nbmanip_comb_rules(string,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_topology_dat
    use prmfile
    use ffdev_utils

    implicit none
    character(PRMFILE_MAX_PATH) :: string
    logical                     :: noexec
    ! --------------------------------------------
    integer                     :: i,comb_rules
    ! --------------------------------------------------------------------------

    comb_rules = ffdev_topology_get_comb_rules_from_string(string)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'NB combination rules', '%')
    write(DEV_OUT,10)  ffdev_topology_comb_rules_to_string(comb_rules)

    if( noexec ) return ! do not execute

    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,20) i
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Original NB parameters', '*')
        call ffdev_topology_info_types(sets(i)%top,1)

        ! remix parameters
        call ffdev_topology_apply_NB_comb_rules(sets(i)%top,comb_rules)

        ! new set of parameters
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'New NB parameters', '*')
        call ffdev_topology_info_types(sets(i)%top,2)
    end do

10 format('Combination rules (comb_rules)   = ',A)
20 format('=== SET ',I2.2)

end subroutine ffdev_parameters_ctrl_nbmanip_comb_rules

! ==============================================================================
! subroutine ffdev_parameters_ctrl_nbmanip_nb_mode
! ==============================================================================

subroutine ffdev_parameters_ctrl_nbmanip_nb_mode(string,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_targetset
    use ffdev_topology_dat
    use prmfile
    use ffdev_utils
    use ffdev_mmd3

    implicit none
    character(PRMFILE_MAX_PATH) :: string
    logical                     :: noexec
    ! --------------------------------------------
    integer                     :: i,nb_mode
    ! --------------------------------------------------------------------------

    if( trim(string) .eq. 'ADDMMD3' ) then
        nb_mode = NB_MODE_ADDMMD3
    else
        nb_mode = ffdev_topology_nb_mode_from_string(string)
    end if

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'NB mode', '%')
    write(DEV_OUT,10)  ffdev_topology_nb_mode_to_string(nb_mode)

    if( noexec ) then
        ! do not execute but minic that we changed the mode
        if( nb_mode .eq. NB_MODE_ADDMMD3 ) then
            nb_mode = NB_MODE_MMD3
        end if
        LastNBMode = nb_mode
        return
    end if

    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,20) i

        if( (nb_mode .eq. NB_MODE_MMD3) .or. (nb_mode .eq. NB_MODE_ADDMMD3) ) then
            call ffdev_mmd3_print_params(sets(i)%top)
        end if

        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Original NB parameters', '*')
        call ffdev_topology_info_types(sets(i)%top,1)

        ! remix parameters
        call ffdev_topology_switch_nbmode(sets(i)%top,nb_mode)

        ! new set of parameters
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'New NB parameters', '*')
        call ffdev_topology_info_types(sets(i)%top,2)
    end do

    if( nb_mode .eq. NB_MODE_ADDMMD3 ) then
        nb_mode = NB_MODE_MMD3
    end if
    LastNBMode = nb_mode

    ! update nb parameters
    call ffdev_parameters_reinit_nbparams

10 format('NB mode (nb_mode)                = ',A)
20 format('=== SET ',I2.2)

end subroutine ffdev_parameters_ctrl_nbmanip_nb_mode

! ==============================================================================
! subroutine ffdev_parameters_ctrl_nbload
! ==============================================================================

subroutine ffdev_parameters_ctrl_nbload(fin,noexec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_topology_dat
    use ffdev_topology
    use prmfile
    use ffdev_utils
    use ffdev_targetset_dat
    use ffdev_targetset    

    implicit none
    type(PRMFILE_TYPE)          :: fin
    logical                     :: noexec
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: line, sti, stj, snb_mode
    real(DEVDP)                 :: a, b, c, d
    integer                     :: i,j,nbt,nb_mode
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

! print header
    write(DEV_OUT,*)
    write(DEV_OUT,510)
    write(DEV_OUT,*)
    write(DEV_OUT,520)
    write(DEV_OUT,530)

! read lines

    do while( prmfile_get_line(fin,line) )

        ! read nb_mode and types first
        read(line,*,err=100,end=100) snb_mode, sti, stj

        nb_mode = ffdev_topology_nb_mode_from_string(snb_mode)

        a = 0.0d0
        b = 0.0d0
        c = 0.0d0
        d = 0.0d0

        ! read the entire record
        select case(nb_mode)
            case(NB_MODE_LJ,NB_MODE_EXPONLY)
                read(line,*,err=100,end=100) snb_mode, sti, stj, a, b
            case(NB_MODE_EXP6,NB_MODE_BP)
                read(line,*,err=100,end=100) snb_mode, sti, stj, a, b, c
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported nb_mode in ffdev_parameters_ctrl_nbload!')
        end select

        write(DEV_OUT,540) trim(snb_mode),trim(sti),trim(stj), a,b,c,d

        if( .not. noexec ) then
            ! for each topology update given parameters
            do j=1,nsets
                if( sets(j)%top%nb_mode .ne. nb_mode ) cycle ! skip topologies with different nb_mode
                nbt = ffdev_topology_find_nbtype_by_types(sets(j)%top,sti,stj)
                if( nbt .ne. 0 ) then
                    select case(nb_mode)
                        case(NB_MODE_LJ)
                            sets(j)%top%nb_types(nbt)%eps = a
                            sets(j)%top%nb_types(nbt)%r0 = b
                        case(NB_MODE_EXP6)
                            sets(j)%top%nb_types(nbt)%eps = a
                            sets(j)%top%nb_types(nbt)%r0 = b
                            sets(j)%top%nb_types(nbt)%alpha = c
                        case(NB_MODE_BP)
                            sets(j)%top%nb_types(nbt)%a = a
                            sets(j)%top%nb_types(nbt)%b = b
                            sets(j)%top%nb_types(nbt)%c6 = c
                        case(NB_MODE_EXPONLY)
                            sets(j)%top%nb_types(nbt)%a = a
                            sets(j)%top%nb_types(nbt)%b = b
                        case default
                            call ffdev_utils_exit(DEV_OUT,1,'Unsupported nb_mode in ffdev_parameters_ctrl_nbload!')
                    end select

                end if
            end do

        end if

    end do

    if( noexec ) return

! print updated parameters
    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,20) i
        ! new set of parameters
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'New NB parameters', '*')
        call ffdev_topology_info_types(sets(i)%top,2)
    end do

! update nb parameters
    call ffdev_parameters_reinit_nbparams

    return

10 format('=== [nbload] ===================================================================')

20 format('=== SET ',I2.2)

100 call ffdev_utils_exit(DEV_OUT,1,'Unable parse line "' // trim(line) // '" in ffdev_parameters_ctrl_nbload!')

510 format('# ~~~~~~~~~~~~~~~~~ NB types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
520 format('# NBMode TypA TypB      eps/A             R0/B          alpha/C6             C8       ')
530 format('# ------ ---- ---- ---------------- ---------------- ----------------  ---------------')
540 format(2X,A6,1X,A4,1X,A4,1X,F16.7,1X,F16.7,1X,F16.7,1X,F16.7)

end subroutine ffdev_parameters_ctrl_nbload

! ------------------------------------------------------------------------------

end module ffdev_parameters_control

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
    integer                     :: i
    character(PRMFILE_MAX_PATH) :: string,realm
    logical                     :: rst
    character(10)               :: key
    ! --------------------------------------------------------------------------

    ! disable all parameter realms
    call ffdev_parameters_disable_all_realms()

    write(DEV_OUT,*)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    if( .not. prmfile_open_section(fin,'parameters') ) then
        call ffdev_utils_exit(DEV_OUT,1,'No parameters to optimize! No [parameters] section found!')
    end if

    rst = prmfile_first_line(fin)
    do while ( prmfile_get_line(fin,string) )
        read(string,*) key, realm
        select case(key)
            case('enable')
                call change_realms(realm,.true.,string)
                write(DEV_OUT,30) adjustl(key),trim(realm)
            case('disable')
                call change_realms(realm,.false.,string)
                write(DEV_OUT,30) adjustl(key),trim(realm)
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
20 format('disable                                           all            (initial setup)')
30 format(A10,40X,A)
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

subroutine change_realms(realm,enable,options)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    character(*)    :: realm
    logical         :: enable
    character(*)    :: options
    ! --------------------------------------------
    integer         :: realmid, i
    ! --------------------------------------------------------------------------

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
        case('all')
            realmid = -1
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Unsupported realm key '''//trim(realm)//'''!')
    end select

    do i=1,nparams
        if( params(i)%realm .eq. realmid ) then
            if( params(i)%identity .eq. 0 ) then
                params(i)%enabled = enable
            end if
        end if
        if( realmid .eq. -1 ) then
            if( params(i)%identity .eq. 0 ) then
                params(i)%enabled = enable
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
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string
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

    if( prmfile_get_string_by_key(fin,'comb_rules',string) ) then
        write(DEV_OUT,*)
        select case(trim(string))
            case('LB')
                FinalCombiningRule = COMB_RULE_LB
                write(DEV_OUT,40) 'LB (Lorentz-Berthelot)'
            case('WH')
                FinalCombiningRule = COMB_RULE_WH
                write(DEV_OUT,40) 'WH (Waldman-Hagler)'
            case('KG')
                FinalCombiningRule = COMB_RULE_KG
                write(DEV_OUT,40) 'KG (Kong)'
            case('FB')
                FinalCombiningRule = COMB_RULE_FB
                write(DEV_OUT,40) 'FB (Fender-Halsey-Berthelot)'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported comb_rules in ffdev_parameters_ctrl_files!')
        end select
    else
        select case(FinalCombiningRule)
            case(COMB_RULE_LB)
                write(DEV_OUT,40) 'LB (Lorentz-Berthelot)'
            case(COMB_RULE_WH)
                write(DEV_OUT,40) 'WH (Waldman-Hagler)'
            case(COMB_RULE_KG)
                write(DEV_OUT,40) 'KG (Kong)'
            case(COMB_RULE_FB)
                write(DEV_OUT,40) 'FB (Fender-Halsey-Berthelot)'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unsupported comb_rules in ffdev_parameters_ctrl_files!')
        end select
    end if

    return

10 format('Input parameters (input)               = ',a)
15 format('Input parameters (input)               = ',a12,'                  (default)')
20 format('Final parameters file (final)          = ',a)
25 format('Final parameters file (final)          = ',a12,'                  (default)')
30 format('Final AMBER parameters (amber)         = ',a)
35 format('Final AMBER parameters (amber)         = ',a12,'                  (default)')
40 format('Combining rules (comb_rules)           = ',a)
45 format('Combining rules (comb_rules)           = ',a12,'                  (default)')

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

! global parameters
    if( prmfile_get_real8_by_key(fin,'erep1',Erep1) ) then
        write(DEV_OUT,20) Erep1
    else
        write(DEV_OUT,25) Erep1
    end if
    if( prmfile_get_real8_by_key(fin,'erep2',Erep2) ) then
        write(DEV_OUT,30) Erep2
    else
        write(DEV_OUT,35) Erep2
    end if

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
20 format('Erep1 (erep1)                    = ',F10.6)
25 format('Erep1 (erep1)                    = ',F10.6,' (default)')
30 format('Erep2 (erep2)                    = ',F10.6)
35 format('Erep2 (erep2)                    = ',F10.6,' (default)')

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
    use ffdev_topology_dat
    use prmfile
    use ffdev_utils

    implicit none
    character(PRMFILE_MAX_PATH) :: string
    logical                     :: noexec
    ! --------------------------------------------
    integer                     :: i,nb_mode
    ! --------------------------------------------------------------------------

    nb_mode = ffdev_topology_nb_mode_from_string(string)

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'NB mode', '%')
    write(DEV_OUT,10)  ffdev_topology_nb_mode_to_string(nb_mode)

    if( noexec ) return ! do not execute

    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,20) i

        if( (nb_mode .eq. NB_MODE_EXPD3BJ) .or. (nb_mode .eq. NB_MODE_ADDD3BJ) ) then
            call ffdev_topology_print_dftd3_params(sets(i)%top)
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

10 format('NB mode (nb_mode)                = ',A)
20 format('=== SET ',I2.2)

end subroutine ffdev_parameters_ctrl_nbmanip_nb_mode

! ------------------------------------------------------------------------------

end module ffdev_parameters_control

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
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_parameters_ctrl_control
! ==============================================================================

subroutine ffdev_parameters_ctrl_control(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils
    use ffdev_topology_dat
    use ffdev_topology
    use ffdev_ffopt
    use ffdev_targetset_dat

    implicit none
    type(PRMFILE_TYPE)          :: fin
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    ! this is for nb_params = normal
    ApplyCombiningRules = .false.

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'control') ) then

        write(DEV_OUT,25) ffdev_parameters_prmsmode_to_string(NBParamsMode)
        write(DEV_OUT,115) prmfile_onoff(PACAsPrms)
        write(DEV_OUT,335) trim(ffdev_topology_nb_mode_to_string(PACSource))

        write(DEV_OUT,65)  prmfile_onoff(OnlyDefinedDihItems)
        write(DEV_OUT,75)  prmfile_onoff(LockDihC_PN1)
        write(DEV_OUT,85)  prmfile_onoff(ResetAllSetup)
        write(DEV_OUT,245) trim(LoadEnergy)
        write(DEV_OUT,445) trim(LoadSAPT)
        write(DEV_OUT,255) trim(LoadProbe)
        write(DEV_OUT,265) max_probe_energy
        write(DEV_OUT,235) trim(LoadCharges)
        write(DEV_OUT,95)  GlbRngSeed
        write(DEV_OUT,105) Verbosity

        if( NBParamsMode .eq. NB_PARAMS_MODE_LIKE_ONLY ) then
            ApplyCombiningRules = .true.
        end if
        return
    end if

    if( prmfile_get_string_by_key(fin,'nb_params', string)) then
        NBParamsMode = ffdev_parameters_prmsmode_from_string(string)
        write(DEV_OUT,20) ffdev_parameters_prmsmode_to_string(NBParamsMode)
    else
        write(DEV_OUT,25) ffdev_parameters_prmsmode_to_string(NBParamsMode)
    end if

    if( NBParamsMode .eq. NB_PARAMS_MODE_LIKE_ONLY ) then
        ApplyCombiningRules = .true.
    end if

    if( prmfile_get_logical_by_key(fin,'pac_as_prms', PACAsPrms)) then
        write(DEV_OUT,110) prmfile_onoff(PACAsPrms)
    else
        write(DEV_OUT,115) prmfile_onoff(PACAsPrms)
    end if

    if( prmfile_get_string_by_key(fin,'pac_source', string)) then
        PACSource = ffdev_parameters_pac_source_from_string(string)
        write(DEV_OUT,330) trim(ffdev_parameters_pac_source_to_string(PACSource))
    else
        write(DEV_OUT,335) trim(ffdev_parameters_pac_source_to_string(PACSource))
    end if

    if( prmfile_get_logical_by_key(fin,'dih_only_defined', OnlyDefinedDihItems)) then
        write(DEV_OUT,60) prmfile_onoff(OnlyDefinedDihItems)
    else
        write(DEV_OUT,65) prmfile_onoff(OnlyDefinedDihItems)
    end if

    if( prmfile_get_logical_by_key(fin,'lock_dihc_pn1', LockDihC_PN1)) then
        write(DEV_OUT,70) prmfile_onoff(LockDihC_PN1)
    else
        write(DEV_OUT,75) prmfile_onoff(LockDihC_PN1)
    end if

    if( prmfile_get_logical_by_key(fin,'resetallsetup', ResetAllSetup)) then
        write(DEV_OUT,80) prmfile_onoff(ResetAllSetup)
    else
        write(DEV_OUT,85) prmfile_onoff(ResetAllSetup)
    end if

    ! setup energy, which should be loaded
    if( prmfile_get_string_by_key(fin,'load_energy',LoadEnergy) ) then
        write(DEV_OUT,240) trim(LoadEnergy)
    else
        write(DEV_OUT,245) trim(LoadEnergy)
    end if

    ! setup sapt, which should be loaded
    if( prmfile_get_string_by_key(fin,'load_sapt',LoadSAPT) ) then
        write(DEV_OUT,440) trim(LoadSAPT)
    else
        write(DEV_OUT,445) trim(LoadSAPT)
    end if

    ! setup energy, which should be loaded for probes
    if( prmfile_get_string_by_key(fin,'load_probe',LoadProbe) ) then
        write(DEV_OUT,250) trim(LoadProbe)
    else
        write(DEV_OUT,255) trim(LoadProbe)
    end if

    if( prmfile_get_real8_by_key(fin,'max_probe_energy', max_probe_energy)) then
        write(DEV_OUT,260) max_probe_energy
    else
        write(DEV_OUT,265) max_probe_energy
    end if

    ! setup charges, which should be loaded
    if( prmfile_get_string_by_key(fin,'load_charges',LoadCharges) ) then
        write(DEV_OUT,230) trim(LoadCharges)
    else
        write(DEV_OUT,235) trim(LoadCharges)
    end if

    if( prmfile_get_integer_by_key(fin,'seed', GlbRngSeed)) then
        write(DEV_OUT,90) GlbRngSeed
    else
        write(DEV_OUT,95) GlbRngSeed
    end if

    if( prmfile_get_integer_by_key(fin,'verbosity', Verbosity)) then
        write(DEV_OUT,100) Verbosity
    else
        write(DEV_OUT,105) Verbosity
    end if

    if( GlbRngSeed .le. 0 ) then
        write(DEV_OUT,*)
        call system_clock(GlbRngSeed)
        write(DEV_OUT,96)  GlbRngSeed
    end if

    ! setup random number generator
    call random_seed(GlbRngSeed)
    call ffdev_ffopt_setup_rng(GlbRngSeed)

    return

 10 format('=== [control] ==================================================================')

 20  format ('NB parameter assembly mode (nb_params)   = ',a12)
 25  format ('NB parameter assembly mode (nb_params)   = ',a12,'                (default)')
110  format ('Partial charges as params (pac_as_prms)  = ',a12)
115  format ('Partial charges as params (pac_as_prms)  = ',a12,'                (default)')
 60  format ('Use defined dih items (dih_only_defined) = ',a12)
 65  format ('Use defined dih items (dih_only_defined) = ',a12,'                (default)')
 70  format ('Lock PN1 for dih_c (lock_dihc_pn1)       = ',a12)
 75  format ('Lock PN1 for dih_c (lock_dihc_pn1)       = ',a12,'                (default)')
 80  format ('Reset all setup (resetallsetup)          = ',a12)
 85  format ('Reset all setup (resetallsetup)          = ',a12,'                (default)')
 90  format ('Random number generator seed (seed)      = ',I12)
 95  format ('Random number generator seed (seed)      = ',I12,'                (default)')
 96  format ('Random number generator seed (seed)      = ',I12,'      (from system clock)')
100  format ('Verbosity (verbosity)                    = ',I12)
105  format ('Verbosity (verbosity)                    = ',I12,'                (default)')

240  format ('Load energy (load_energy)                = ',A12)
245  format ('Load energy (load_energy)                = ',A12,'                (default)')

440  format ('Load SAPT (load_sapt)                    = ',A12)
445  format ('Load SAPT (load_sapt)                    = ',A12,'                (default)')

250  format ('Load probe energy (load_probe)           = ',A12)
255  format ('Load probe energy (load_probe)           = ',A12,'                (default)')

260  format ('Max probe energy (max_probe_energy)      = ',F12.3)
265  format ('Max probe energy (max_probe_energy)      = ',F12.3,'                (default)')

230  format ('Load charges (load_charges)              = ',A12)
235  format ('Load charges (load_charges)              = ',A12,'                (default)')

330  format ('PAC source for stat (pac_source)         = ',A12)
335  format ('PAC source for stat (pac_source)         = ',A12,'                (default)')

end subroutine ffdev_parameters_ctrl_control

! ==============================================================================
! subroutine ffdev_parameters_prmsmode_to_string
! ==============================================================================

character(80) function ffdev_parameters_prmsmode_to_string(source)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer  :: source
    ! --------------------------------------------------------------------------

    select case(source)
        case(NB_PARAMS_MODE_LIKE_ONLY)
            ffdev_parameters_prmsmode_to_string = 'like-only - only nbii parameters are considered'
        case(NB_PARAMS_MODE_ALL)
            ffdev_parameters_prmsmode_to_string = 'all - all nbparameters are considered'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_prmsmode_to_string!')
    end select

end function ffdev_parameters_prmsmode_to_string

! ==============================================================================
! subroutine ffdev_parameters_prmsmode_from_string
! ==============================================================================

integer function ffdev_parameters_prmsmode_from_string(string)

    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('like-only')
            ffdev_parameters_prmsmode_from_string = NB_PARAMS_MODE_LIKE_ONLY
        case('all')
            ffdev_parameters_prmsmode_from_string = NB_PARAMS_MODE_ALL
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) &
                                            // '" in ffdev_parameters_prmsmode_from_string!')
    end select

end function ffdev_parameters_prmsmode_from_string

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

    ! input file cannot be setup here
    ! the correct parametr data can be loaded only after parameters are properly setup

    if(.not. prmfile_open_section(fin,'files')) then
        write (DEV_OUT,25) trim(OutParamFileName)
        write (DEV_OUT,35) trim(OutAmberPrmsFileName)
        return
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

20 format('Final parameters file (final)          = ',a)
25 format('Final parameters file (final)          = ',a12,'                  (default)')
30 format('Final AMBER parameters (amber)         = ',a)
35 format('Final AMBER parameters (amber)         = ',a12,'                  (default)')

end subroutine ffdev_parameters_ctrl_files

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_parameters_ctrl_files_exec(fin,exec)

    use prmfile
    use ffdev_parameters_dat
    use ffdev_parameters
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)      :: fin
    logical                 :: exec
    ! ---------------------------------------------
    logical                 :: rst
    character(MAX_PATH)     :: key
    character(MAX_PATH)     :: name
    integer                 :: loaded,ignored
    ! -----------------------------------------------------------------------------

    if(.not. prmfile_open_section(fin,'files')) then
        return
    end if

    write(DEV_OUT,'(/,a)') '=== [files] ===================================================================='

    loaded = 0
    ignored = 0

    rst = prmfile_first_line(fin)
    do while( rst )
        rst = prmfile_get_field_on_line(fin,key)
        select case(trim(key))
            case('load')
                rst = prmfile_get_field_on_line(fin,name)
                if( .not. rst ) then
                    call ffdev_utils_exit(DEV_ERR,1,'Name of the file with FFDevel parameters to be loaded is not provided!')
                end if
                write(DEV_OUT,10) trim(name)
                if( exec ) then
                    call ffdev_parameters_load(name,loaded,ignored)
                    call ffdev_parameters_to_tops
                end if
            case('save')
                rst = prmfile_get_field_on_line(fin,name)
                if( .not. rst ) then
                    call ffdev_utils_exit(DEV_ERR,1,'Name of the file with FFDevel parameters to be saved is not provided!')
                end if
                write(DEV_OUT,10) trim(name)
                if( exec ) then
                    call ffdev_parameters_save(name)
                end if
            case('amber')
                rst = prmfile_get_field_on_line(fin,name)
                if( .not. rst ) then
                    call ffdev_utils_exit(DEV_ERR,1,'Name of the file with AMBER parameters to be saved is not provided!')
                end if
                write(DEV_OUT,10) trim(name)
                if( exec ) then
                    call ffdev_parameters_save_amber(name)
                end if
            case default
                ! do nothing
        end select
        rst = prmfile_next_line(fin)
    end do

    write(DEV_OUT,12) loaded
    write(DEV_OUT,14) ignored

    return

10 format('Load parameters file (load)       = ',a)
12 format('     Number of loaded parameters  = ',I6)
14 format('     Number of ignored parameters = ',I6)

end subroutine ffdev_parameters_ctrl_files_exec

! ==============================================================================
! subroutine ffdev_parameters_ctrl_grbf2cos
! ==============================================================================

subroutine ffdev_parameters_ctrl_grbf2cos(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils
    use ffdev_topology_dat
    use ffdev_topology

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'grbf2cos') ) then
        write(DEV_OUT,25) GRBF2COSMaxN
        write(DEV_OUT,35) GRBF2COSMinV
        return
    end if

    if( prmfile_get_integer_by_key(fin,'max_n', GRBF2COSMaxN)) then
        write(DEV_OUT,20) GRBF2COSMaxN
    else
        write(DEV_OUT,25) GRBF2COSMaxN
    end if

    if( prmfile_get_real8_by_key(fin,'min_v', GRBF2COSMinV)) then
        write(DEV_OUT,30) GRBF2COSMinV
    else
        write(DEV_OUT,35) GRBF2COSMinV
    end if

    return

 10 format('=== [grbf2cos] =================================================================')

 20  format ('Maximum length of cos series (max_n)     = ',i12)
 25  format ('Maximum length of cos series (max_n)     = ',i12,'                (default)')
 30  format ('Minimum accepted amplitude (min_v)       = ',f12.7)
 35  format ('Minimum accepted amplitude (min_v)       = ',f12.7,'                (default)')

end subroutine ffdev_parameters_ctrl_grbf2cos

! ==============================================================================
! subroutine ffdev_parameters_ctrl_ffmanip
! ==============================================================================

subroutine ffdev_parameters_ctrl_ffmanip(fin,exec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils
    use ffdev_disp_control
    use ffdev_buried_control
    use ffdev_atomicdata_control

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: exec
    ! --------------------------------------------
    logical                     :: rst
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    rst = prmfile_first_section(fin)
    do while( rst )

        ! open set section
        if( .not. prmfile_get_section_name(fin,string) ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to get section name!')
        end if

        select case(string)
            case('bond_r0')
                call ffdev_parameters_ctrl_bond_r0(fin,exec)
            case('angle_a0')
                call ffdev_parameters_ctrl_angle_a0(fin,exec)
            case('nbsetup')
                call ffdev_parameters_ctrl_nbsetup(fin,exec)
            case('parameters')
                call ffdev_parameters_ctrl_realms(fin)
            case('setprms')
                call ffdev_parameters_ctrl_setprms(fin,exec)
            case('disp')
                call ffdev_disp_ctrl(fin,exec)
            case('buried')
                call ffdev_buried_ctrl(fin,exec)
            case('atomicdata')
                call ffdev_atomicdata_ctrl(fin,exec)
            case('files')
                call ffdev_parameters_ctrl_files_exec(fin,exec)
        end select

        rst = prmfile_next_section(fin)
    end do

end subroutine ffdev_parameters_ctrl_ffmanip

! ==============================================================================
! subroutine ffdev_parameters_ctrl_nbsetup
! ==============================================================================

subroutine ffdev_parameters_ctrl_nbsetup(fin,exec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils
    use ffdev_topology
    use ffdev_topology_tt
    use ffdev_topology_bj
    use ffdev_topology_lj
    use ffdev_topology_exp
    use ffdev_topology_pen
    use ffdev_topology_ind
    use ffdev_topology_dat
    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_nb2nb

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: exec
    ! --------------------------------------------
    logical                     :: changed
    character(PRMFILE_MAX_PATH) :: string
    integer                     :: i
    integer                     :: ldampbj_mode,lexp_pb_mode
    integer                     :: from_nb_mode,to_nb_mode,lpen_mode,lind_mode
    integer                     :: lcomb_rules,lexp_mode
    integer                     :: lexp_pa_mode
    logical                     :: lpen_enabled,lind_enabled
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'nbsetup') ) then
        write(DEV_OUT,25) ffdev_topology_nb_mode_to_string(nb_mode)

        if( nb_mode .eq. NB_VDW_EXP_DISPBJ .or. nb_mode .eq. NB_VDW_EXP_DISPTT ) then
            write(DEV_OUT,85)  ffdev_topology_EXP_exp_mode_to_string(exp_mode)
            write(DEV_OUT,65)  ffdev_topology_EXP_pa_mode_to_string(exp_pa_mode)
            write(DEV_OUT,75)  ffdev_topology_EXP_pb_mode_to_string(exp_pb_mode)
            write(DEV_OUT,95)  prmfile_onoff(pen_enabled)
            write(DEV_OUT,125) ffdev_topology_PEN_mode_to_string(pen_mode)
            write(DEV_OUT,215) prmfile_onoff(ind_enabled)
            write(DEV_OUT,225) ffdev_topology_IND_mode_to_string(ind_mode)
        end if

        if( nb_mode .eq. NB_VDW_EXP_DISPBJ ) then
            write(DEV_OUT,45) ffdev_topology_BJ_dampbj_mode_to_string(dampbj_mode)
        end if

        if( ApplyCombiningRules ) then
            if( nb_mode .eq. NB_VDW_LJ ) then
                write(DEV_OUT,35) ffdev_topology_LJ_comb_rules_to_string(lj_comb_rules)
            end if
        end if
        return
    end if

    changed = .false.

! ---------------------------
    to_nb_mode = nb_mode
    if( prmfile_get_string_by_key(fin,'nb_mode', string)) then
        to_nb_mode = ffdev_topology_nb_mode_from_string(string)
        write(DEV_OUT,20) ffdev_topology_nb_mode_to_string(to_nb_mode)
        if( exec ) then
            from_nb_mode = nb_mode
            changed = .true.
        end if
    else
        write(DEV_OUT,25) ffdev_topology_nb_mode_to_string(nb_mode)
    end if

! ---------------------------
    if( (to_nb_mode .eq. NB_VDW_EXP_DISPTT) .or. (to_nb_mode .eq. NB_VDW_EXP_DISPBJ) ) then

    ! EXP
        if( prmfile_get_string_by_key(fin,'exp_mode', string)) then
            lexp_mode = ffdev_topology_EXP_exp_mode_from_string(string)
            write(DEV_OUT,80) ffdev_topology_EXP_exp_mode_to_string(lexp_mode)
            if( exec ) then
                exp_mode = lexp_mode
                changed = .true.
            end if
        else
            write(DEV_OUT,85) ffdev_topology_EXP_exp_mode_to_string(exp_mode)
        end if
        if( prmfile_get_string_by_key(fin,'exp_pa_mode', string)) then
            lexp_pa_mode = ffdev_topology_EXP_pa_mode_from_string(string)
            write(DEV_OUT,60) ffdev_topology_EXP_pa_mode_to_string(lexp_pa_mode)
            if( exec ) then
                exp_pa_mode = lexp_pa_mode
                changed = .true.
            end if
        else
            write(DEV_OUT,65) ffdev_topology_EXP_pa_mode_to_string(exp_pa_mode)
        end if
        if( prmfile_get_string_by_key(fin,'exp_pb_mode', string)) then
            lexp_pb_mode = ffdev_topology_EXP_pb_mode_from_string(string)
            write(DEV_OUT,70) ffdev_topology_EXP_pb_mode_to_string(lexp_pb_mode)
            if( exec ) then
                exp_pb_mode = lexp_pb_mode
                changed = .true.
            end if
        else
            write(DEV_OUT,75) ffdev_topology_EXP_pb_mode_to_string(exp_pb_mode)
        end if

    ! PEN
        if( prmfile_get_logical_by_key(fin,'pen_enabled', lpen_enabled)) then
            write(DEV_OUT,90) prmfile_onoff(lpen_enabled)
            if( exec ) then
                pen_enabled = lpen_enabled
                changed = .true.
            end if
        else
            write(DEV_OUT,95) prmfile_onoff(pen_enabled)
        end if
    !----------
        if( prmfile_get_string_by_key(fin,'pen_mode', string)) then
            lpen_mode = ffdev_topology_PEN_mode_from_string(string)
            write(DEV_OUT,120) ffdev_topology_PEN_mode_to_string(lpen_mode)
            if( exec ) then
                pen_mode = lpen_mode
                changed = .true.
            end if
        else
            write(DEV_OUT,125) ffdev_topology_PEN_mode_to_string(pen_mode)
        end if

   ! IND
        if( prmfile_get_logical_by_key(fin,'ind_enabled', lind_enabled)) then
            write(DEV_OUT,210) prmfile_onoff(lind_enabled)
            if( exec ) then
                ind_enabled = lind_enabled
                changed = .true.
            end if
        else
            write(DEV_OUT,215) prmfile_onoff(ind_enabled)
        end if
    !----------
        if( prmfile_get_string_by_key(fin,'ind_mode', string)) then
            lind_mode = ffdev_topology_IND_mode_from_string(string)
            write(DEV_OUT,220) ffdev_topology_IND_mode_to_string(lind_mode)
            if( exec ) then
                ind_mode = lind_mode
                changed = .true.
            end if
        else
            write(DEV_OUT,225) ffdev_topology_IND_mode_to_string(ind_mode)
        end if

    end if

! ---------------------------
    if( to_nb_mode .eq. NB_VDW_EXP_DISPBJ ) then
        if( prmfile_get_string_by_key(fin,'dampbj_mode', string)) then
            ldampbj_mode = ffdev_topology_BJ_dampbj_mode_from_string(string)
            write(DEV_OUT,40) ffdev_topology_BJ_dampbj_mode_to_string(ldampbj_mode)
            if( exec ) then
                dampbj_mode = ldampbj_mode
                changed = .true.
            end if
        else
            write(DEV_OUT,45) ffdev_topology_BJ_dampbj_mode_to_string(dampbj_mode)
        end if
    end if


! ---------------------------
    if( ApplyCombiningRules ) then
        if( to_nb_mode .eq. NB_VDW_LJ ) then
            if( prmfile_get_string_by_key(fin,'lj_comb_rules', string)) then
                lcomb_rules = ffdev_topology_LJ_comb_rules_from_string(string)
                write(DEV_OUT,30) ffdev_topology_LJ_comb_rules_to_string(lcomb_rules)
                if( exec ) then
                    lj_comb_rules = lcomb_rules
                    changed = .true.
                end if
            else
                write(DEV_OUT,35) ffdev_topology_LJ_comb_rules_to_string(lj_comb_rules)
            end if
        end if
    end if

! ---------------------------
    if( exec .and. changed ) then
        nb_mode = to_nb_mode
        do i=1,nsets
            ! update NB parameters
            call ffdev_topology_update_nb_params(sets(i)%top)

            if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
                ! new set of parameters
                write(DEV_OUT,*)
                write(DEV_OUT,15) i
                call ffdev_utils_heading(DEV_OUT,'New NB parameters', '*')
                call ffdev_topology_info_types(sets(i)%top,2)
            end if
        end do

        call ffdev_parameters_reinit
        call ffdev_targetset_reinit_nbparams
    end if

 10 format('=== [nbsetup] ==================================================================')

 20 format('NB mode (nb_mode)                    = ',A)
 25 format('NB mode (nb_mode)                    = ',A31,' (current)')

 30 format('LJ comb. rules (lj_comb_rules)       = ',A)
 35 format('LJ combining rules (lj_comb_rules)   = ',A31,' (current)')

 40 format('BJ damping mode (dampbj_mode)        = ',A)
 45 format('BJ damping mode (dampbj_mode)        = ',A31,' (current)')

 80 format('EXP mode (exp_mode)                  = ',A)
 85 format('EXP mode (exp_mode)                  = ',A31,' (current)')

 60 format('EXP PA mode (exp_pa_mode)            = ',A)
 65 format('EXP PA mode (exp_pa_mode)            = ',A31,' (current)')

 70 format('EXP PB mode (exp_pb_mode)            = ',A)
 75 format('EXP PB mode (exp_pb_mode)            = ',A31,' (current)')

 90 format('Penetration energy (pen_enabled)     = ',A)
 95 format('Penetration energy (pen_enabled)     = ',A31,' (current)')

120 format('Pen energy mode (pen_mode)           = ',A)
125 format('Pen energy mode (pen_mode)           = ',A31,' (current)')

210 format('Induction energy (ind_enabled)       = ',A)
215 format('Induction energy (ind_enabled)       = ',A31,' (current)')

220 format('Induction energy mode (ind_mode)     = ',A)
225 format('Induction energy mode (ind_mode)     = ',A31,' (current)')

 15 format('=== SET ',I2.2)

end subroutine ffdev_parameters_ctrl_nbsetup

! ==============================================================================
! subroutine ffdev_parameters_ctrl_setprms
! ==============================================================================

subroutine ffdev_parameters_ctrl_setprms(fin,exec)

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
    logical                     :: exec
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: line, sti, stj, stk, stl, realm
    real(DEVDP)                 :: lvalue
    integer                     :: pn, realmid, ti, tj, tk, tl, parmid
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

! read lines

    write(DEV_OUT,*)
    write(DEV_OUT,15)
    write(DEV_OUT,20)

    do while( prmfile_get_line(fin,line) )

        ! read realm
        read(line,*,err=100,end=100) realm
        realmid = ffdev_parameters_get_realmid(realm)

        sti = ''
        stj = ''
        stk = ''
        stl = ''
        pn  = 0

        select case(realmid)
            case(REALM_BOND_D0,REALM_BOND_K)
                read(line,*,err=100,end=100) realm, sti, stj, lvalue
            case(REALM_ANGLE_A0,REALM_ANGLE_K)
                read(line,*,err=100,end=100) realm, sti, stj, lvalue
            case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G)
                read(line,*,err=100,end=100) realm, sti, stj, stk, stl, pn, lvalue
            case(REALM_DIH_SCEE,REALM_DIH_SCNB)
                read(line,*,err=100,end=100) realm, sti, stj, stk, stl, pn, lvalue
            case(REALM_IMPR_V,REALM_IMPR_G)
                read(line,*,err=100,end=100) realm, sti, stj, stk, stl, pn, lvalue
            case(REALM_VDW_EPS,REALM_VDW_R0)
                read(line,*,err=100,end=100) realm, sti, stj, lvalue
            case(REALM_VDW_PA,REALM_VDW_PB,REALM_VDW_RC)
                read(line,*,err=100,end=100) realm, sti, lvalue
            case(REALM_VDW_B0)
                read(line,*,err=100,end=100) realm, pn, lvalue
            case(REALM_ELE_SQ)
                read(line,*,err=100,end=100) realm, lvalue
            case(REALM_DAMP_FA,REALM_DAMP_FB,REALM_DAMP_PB,REALM_DAMP_TB,REALM_DAMP_PE)
                read(line,*,err=100,end=100) realm, lvalue
            case(REALM_DISP_S6,REALM_DISP_S8,REALM_DISP_S10)
                read(line,*,err=100,end=100) realm, lvalue
            case(REALM_GLB_SCEE,REALM_GLB_SCNB)
                read(line,*,err=100,end=100) realm, lvalue
            case(REALM_K_EXC,REALM_K_IND)
                read(line,*,err=100,end=100) realm, lvalue
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_parameters_ctrl_setprms!')
        end select

        if( exec ) then
            ti = get_common_type_from_name(sti)
            tj = get_common_type_from_name(stj)
            tk = get_common_type_from_name(stk)
            tl = get_common_type_from_name(stl)

            parmid = find_parameter_by_ids(realmid,pn,ti,tj,tk,tl)

            write(DEV_OUT,30) parmid, trim(realm), trim(sti), trim(stj), trim(stk), trim(stl), pn, lvalue

            if( parmid .le. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to find the parameter!')
            end if

            params(parmid)%value = lvalue / ffdev_parameters_get_realm_scaling(realmid)

        end if

    end do

    return

 10 format('=== [setprms] ==================================================================')

 15 format('# ID    Realm    TI TJ TK TL PN       Value     ')
 20 format('# -- ----------- -- -- -- -- -- ----------------')
 30 format(I4,1X,A11,1X,A2,1X,A2,1X,A2,1X,A2,1X,I2,1X,F16.4)

100 call ffdev_utils_exit(DEV_ERR,1,'Unable parse line "' // trim(line) // '" in ffdev_parameters_ctrl_setprms!')


end subroutine ffdev_parameters_ctrl_setprms

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
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported identity key '''//trim(key)//'''!')
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
        call ffdev_utils_exit(DEV_ERR,1,'Unsupported realm for setup_realm_identity!')
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
        call ffdev_utils_exit(DEV_ERR,1,'Unsupported realm for setup_dih_identity_byrot!')
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
    character(MAX_PATH)         :: key
    real(DEVDP)                 :: rnd, minv, maxv, rfact
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    nactparms = 0
    if( ResetAllSetup ) then
        call ffdev_parameters_disable_all_realms()
        write(DEV_OUT,35) nactparms,'all disabled by default (resetallsetup=on)'
    else
        do i=1,nparams
            if( params(i)%enabled ) nactparms = nactparms + 1
        end do
        write(DEV_OUT,35) nactparms,'(kept active from the previous setup)'
    end if

    if( .not. prmfile_open_section(fin,'parameters') ) then
        if( nactparms .eq. 0 ) then
            write(DEV_OUT,*)
            write(DEV_OUT,45)
        end if
        return
    end if

    rst = prmfile_first_line(fin)
    do while ( prmfile_get_line(fin,string) )
        read(string,*,err=555,end=555) key
        nchanged = 0
        select case(key)
            case('enable')
                read(string,*,err=556,end=556) key, realm
                call change_realms(realm,.true.,string,nchanged)
                write(DEV_OUT,30) nchanged,trim(string)
            case('disable')
                read(string,*,err=556,end=556) key, realm
                call change_realms(realm,.false.,string,nchanged)
                write(DEV_OUT,30) nchanged,trim(string)
            case('randomize')
                nchanged = 0
                do i=1,nparams
                    if( .not. params(i)%enabled ) cycle
                    call random_number(rnd)
                    minv = ffdev_params_get_lower_bound(params(i)%realm)
                    maxv = ffdev_params_get_upper_bound(params(i)%realm)
                    params(i)%value = minv + (maxv - minv)*rnd
                    nchanged = nchanged + 1
                end do
                write(DEV_OUT,30) nchanged,trim(string)
            case('randomize-small')
                rfact = 0.01
                nchanged = 0
                do i=1,nparams
                    if( .not. params(i)%enabled ) cycle
                    call random_number(rnd)
                    minv = ffdev_params_get_lower_bound(params(i)%realm)
                    maxv = ffdev_params_get_upper_bound(params(i)%realm)
                    params(i)%value = params(i)%value + rfact*rnd - rfact*0.5d0
                    if( params(i)%value .lt. minv ) params(i)%value = minv
                    if( params(i)%value .gt. maxv ) params(i)%value = maxv
                    nchanged = nchanged + 1
                end do
                write(DEV_OUT,30) nchanged,trim(string)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported action key '''//trim(key)//'''!')
        end select
    end do

    ! count active parameters
    nactparms = 0
    do i=1,nparams
        if( params(i)%enabled ) nactparms = nactparms + 1
    end do

    write(DEV_OUT,*)
    write(DEV_OUT,40) nactparms

    if( nactparms .eq. 0 ) then
        write(DEV_OUT,*)
        write(DEV_OUT,45)
    end if

    return
555 call ffdev_utils_exit(DEV_ERR,1,'Unable to read key in [prameters]!')
556 call ffdev_utils_exit(DEV_ERR,1,'Unable to read key and realm in [prameters]!')

10 format('=== [parameters] ===============================================================')
30 format('Altered parameters = ',I3,' | ', A)
35 format('Active parameters  = ',I3,' | ', A)
40 format('Number of active parameters = ',I6)
45 format('>>> WARNING: No active parameters to optimize!')

end subroutine ffdev_parameters_ctrl_realms

! ==============================================================================
! subroutine ffdev_parameters_ctrl_ranges
! ==============================================================================

subroutine ffdev_parameters_ctrl_ranges(fin)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string,realm,bound
    logical                     :: rst
    integer                     :: realmid
    real(DEVDP)                 :: lvalue, sc
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    if( .not. prmfile_open_section(fin,'ranges') ) then
        write(DEV_OUT,15)
        return
    end if

    if( ResetAllSetup ) then
        call ffdev_params_reset_ranges()
        write(DEV_OUT,'(A)') '>>> INFO: All ranges reset (resetallsetup=on)'
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,30)

    rst = prmfile_first_line(fin)
    do while ( prmfile_get_line(fin,string) )
        read(string,*,err=555,end=555) realm
        if( trim(realm) .eq. 'reset' ) then
            write(DEV_OUT,40) 'reset'
            call ffdev_params_reset_ranges()
        else
            read(string,*,err=556,end=556) realm, bound, lvalue
            realmid = ffdev_parameters_get_realmid(realm)
            write(DEV_OUT,40) ffdev_parameters_get_realm_name(realmid), trim(bound), lvalue
            sc = ffdev_parameters_get_realm_scaling(realmid)
            lvalue = lvalue / sc
            select case(trim(bound))
                case('max')
                    call ffdev_params_set_upper_bound(realmid,lvalue)
                case('min')
                    call ffdev_params_set_lower_bound(realmid,lvalue)
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'min/max required - unsupported boundary in ffdev_parameters_ctrl_ranges!')
            end select
        end if
    end do

    return

 10 format('=== [ranges] ===================================================================')
 15 format('>>> INFO: No changes in parameter ranges!')
 20 format('# Realm        Boundary     Value')
 30 format('# ------------ -------- ----------')
 40 format(A14,1X,A8,1X,F10.2)

555 call ffdev_utils_exit(DEV_ERR,1,'Unable to read key in [ranges]!')
556 call ffdev_utils_exit(DEV_ERR,1,'Unable to read realm, boundary, and/or parameter value in [ranges]!')

end subroutine ffdev_parameters_ctrl_ranges

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

logical function is_realm_nooption(options)

    use prmfile

    implicit none
    character(*)    :: options
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: k1
    character(PRMFILE_MAX_PATH) :: k2
    character(PRMFILE_MAX_PATH) :: k3
    ! --------------------------------------------------------------------------

    is_realm_nooption = .true.
    read(options,*,end=100,err=100) k1,k2,k3

    is_realm_nooption = len(k3) .eq. 0
100 return

end function is_realm_nooption

! ------------------------------------------------------------------------------

subroutine change_realms(realm,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_buried_dat

    implicit none
    character(*)    :: realm
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
    ! --------------------------------------------
    integer         :: realmid, i
    ! --------------------------------------------------------------------------

    if( realm .ne. 'all' ) then
        realmid = ffdev_parameters_get_realmid(realm)
    else
        realmid = -1
    end if

    ! no options
    if( is_realm_nooption(options) ) then
        do i=1,nparams
            if( params(i)%realm .eq. realmid ) then
                if( params(i)%identity .eq. 0 ) then
                    params(i)%enabled = enable
                    nchanged = nchanged + 1
                end if
            end if
        end do
        return
    end if

    ! supported options
    if( is_realm_option(options,'zero') ) then
        select case(realmid)
            case(REALM_BOND_K,REALM_ANGLE_K)
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%value .eq. 0.0d0 ) then
                            if( params(i)%identity .eq. 0 ) then
                                params(i)%enabled = enable
                                nchanged = nchanged + 1
                            end if
                        end if
                    end if
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''zero'' can be used only with bond_k and angle_k')
        end select
! ------------------
    else if( is_realm_option(options,'bond') ) then
        select case(realmid)
            case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G)
                call change_dih_realm_bond(realmid,enable,options,nchanged)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''bond'' can be used only with dih_v, dih_c, and dih_g')
        end select
! ------------------
    else if( is_realm_option(options,'pn') ) then
        select case(realmid)
            case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G,REALM_VDW_B0)
                call change_dih_realm_pn(realmid,enable,options,nchanged)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''pn'' can be used only with dih_v, dih_c, dih_g, wdv_b0')
        end select
! ------------------
    else if( is_realm_option(options,'like') ) then
        select case(realmid)
            case(REALM_VDW_EPS,REALM_VDW_R0)
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%ti .eq. params(i)%tj ) then
                            if( params(i)%identity .eq. 0 ) then
                                params(i)%enabled = enable
                                nchanged = nchanged + 1
                            end if
                        end if
                    end if
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''like'' can be used only with vdw_eps and wdv_r0')
        end select
! ------------------
    else if( is_realm_option(options,'unlike') ) then
        select case(realmid)
            case(REALM_VDW_EPS,REALM_VDW_R0)
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        if( params(i)%ti .ne. params(i)%tj ) then
                            if( params(i)%identity .eq. 0 ) then
                                params(i)%enabled = enable
                                nchanged = nchanged + 1
                            end if
                        end if
                    end if
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''unlike'' can be used only with vdw_eps and wdv_r0')
        end select
! ------------------
    else if( is_realm_option(options,'buried') ) then
        select case(realmid)
            case(REALM_VDW_EPS,REALM_VDW_R0)
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        ! 0.0 - buried, 1.0 - exposed
                        if( (buried_atoms(params(i)%ti)%weight .gt. 0.5d0) .and. &
                            (buried_atoms(params(i)%tj)%weight .gt. 0.5d0) ) cycle
                        if( params(i)%identity .eq. 0 ) then
                            params(i)%enabled = enable
                            nchanged = nchanged + 1
                        end if
                    end if
                end do
            case(REALM_VDW_PA,REALM_VDW_PB)
                do i=1,nparams
                    if( params(i)%realm .eq. realmid ) then
                        ! 0.0 - buried, 1.0 - exposed
                        if( buried_atoms(params(i)%ti)%weight .gt. 0.5d0 ) cycle
                        if( params(i)%identity .eq. 0 ) then
                            params(i)%enabled = enable
                            nchanged = nchanged + 1
                        end if
                    end if
                end do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''buried'' can be used only with vdw_eps, wdv_r0, vdw_pa and vdw_pb')
        end select
! ------------------
    else if( is_realm_option(options,'type') ) then
        select case(realmid)
            case(REALM_VDW_PA,REALM_VDW_PB)
                call change_dih_realm_type(realmid,enable,options,nchanged)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''type'' can be used only with vdw_pa and vdw_pb')
        end select
! ------------------
    else if( is_realm_option(options,'Z') ) then
        select case(realmid)
            case(REALM_VDW_B0)
                call change_dih_realm_Z_for_vdw_b0(realmid,enable,options,nchanged)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option ''Z'' can be used only with vdw_b0')
        end select
! ------------------
    else if( is_realm_option(options,'types') ) then
        select case(realmid)
            case(REALM_BOND_D0,REALM_BOND_K)
            case(REALM_ANGLE_A0,REALM_ANGLE_K)
            case(REALM_DIH_V,REALM_DIH_C,REALM_DIH_G)
            case(REALM_VDW_EPS,REALM_VDW_R0)
                call change_dih_realm_types(realmid,enable,options,nchanged)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Option "pn" can be used only with bond_d0, bond_k,' &
                                                // 'angle_a0, angle_k, dih_v, dih_c, dih_g, vdw_eps, and vdw_r0')
        end select
    else
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported option ('//trim(options)//')')
    end if

end subroutine change_realms

! ------------------------------------------------------------------------------

subroutine change_dih_realm_bond(realmid,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
    ! --------------------------------------------
    integer                     :: sid, aj, ak, i, j
    character(PRMFILE_MAX_PATH) :: k1
    character(PRMFILE_MAX_PATH) :: k2
    character(PRMFILE_MAX_PATH) :: k3
    logical                     :: found
    ! --------------------------------------------------------------------------

    read(options,*,end=100,err=100) k1,k2,k3,sid,aj,ak

    if( (sid .lt. 1) .or. (sid .gt. nsets) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Set ID is out-of-legal range!')
    end if

    if( (aj .lt. 1) .or. (aj .gt. sets(sid)%top%natoms) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Atom i is out-of-legal range!')
    end if

    if( (ak .lt. 1) .or. (ak .gt. sets(sid)%top%natoms) ) then
        call ffdev_utils_exit(DEV_ERR,1,'Atom j is out-of-legal range!')
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
                nchanged = nchanged + 1
            end if
        end if
    end do

    return

100 call ffdev_utils_exit(DEV_ERR,1,'Illegal option ''bond'' for dih_* realm!')

end subroutine change_dih_realm_bond

! ------------------------------------------------------------------------------

subroutine change_dih_realm_pn(realmid,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
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
                nchanged = nchanged + 1
            end if
        end if
    end do

    return

100 call ffdev_utils_exit(DEV_ERR,1,'Illegal option ''pn'' for dih_* realm!')

end subroutine change_dih_realm_pn

! ------------------------------------------------------------------------------

subroutine change_dih_realm_types(realmid,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
    ! --------------------------------------------
    integer         :: i, ti, tj, tk, tl
    character(80)   :: k1, k2, k3, sti, stj, stk, stl
    ! --------------------------------------------------------------------------

    ! read types
    sti = ''
    stj = ''
    stk = ''
    stl = ''
    ti = 0
    tj = 0
    tk = 0
    tl = 0

    read(options,*,end=555,err=555) k1, k2, k3, sti, stj, stk, stl
555     do i=1,ntypes
        if( types(i)%name .eq. sti ) ti = i
        if( types(i)%name .eq. stj ) tj = i
        if( types(i)%name .eq. stk ) tk = i
        if( types(i)%name .eq. stl ) tl = i
    end do

    ! we need at least two types
    if( (tj .eq. 0) .and. (tk .eq. 0) .and. (tl .eq. 0) ) then
        call ffdev_utils_exit(DEV_ERR,1,'At least two types are required to be defined for ''types'' option!')
    end if

    do i=1,nparams

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

        if( (realmid .eq. -1) .or. (params(i)%realm .eq. realmid) ) then
            if( params(i)%identity .eq. 0 ) then
                if( realmid .eq. REALM_DIH_C ) then
                    if( LockDihC_PN1 .and. params(i)%pn .eq. 1 ) then
                        ! do nothing
                    else
                        params(i)%enabled = enable
                        nchanged = nchanged + 1
                    end if
                else
                    params(i)%enabled = enable
                    nchanged = nchanged + 1
                end if
            end if
        end if
    end do

end subroutine change_dih_realm_types

! ------------------------------------------------------------------------------

subroutine change_dih_realm_type(realmid,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
    ! --------------------------------------------
    integer         :: i, ti
    character(80)   :: k1, k2, k3, sti
    ! --------------------------------------------------------------------------

    ! read types
    sti = ''
    ti = 0

    read(options,*,end=100,err=100) k1, k2, k3, sti
    do i=1,ntypes
        if( types(i)%name .eq. sti ) ti = i
    end do

    ! we need at least one type
    if( ti .eq. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'One type is required to be defined for ''type'' option!')
    end if

    do i=1,nparams
        if( (realmid .eq. -1) .or. (params(i)%realm .eq. realmid) ) then
            if( params(i)%ti .ne. ti ) cycle
            if( params(i)%identity .eq. 0 ) then
                params(i)%enabled = enable
                nchanged = nchanged + 1
            end if
        end if
    end do

100 call ffdev_utils_exit(DEV_ERR,1,'Inproper definition of ''type'' option!')

end subroutine change_dih_realm_type

! ------------------------------------------------------------------------------

subroutine change_dih_realm_Z_for_vdw_b0(realmid,enable,options,nchanged)

    use ffdev_parameters
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat
    use prmfile

    implicit none
    integer         :: realmid
    logical         :: enable
    character(*)    :: options
    integer         :: nchanged
    ! --------------------------------------------
    integer         :: i, pn
    character(80)   :: k1, k2, k3
    ! --------------------------------------------------------------------------

    ! read types
    pn = 0

    read(options,*,end=100,err=100) k1, k2, k3, pn

    ! we need at least one type
    if( pn .gt. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Z must be larger than zero for ''Z'' option for vdw_b0!')
    end if

    do i=1,nparams
        if( (realmid .eq. -1) .or. (params(i)%realm .eq. realmid) ) then
            if( params(i)%pn .eq. pn ) then
                if( params(i)%identity .eq. 0 ) then
                    params(i)%enabled = enable
                    nchanged = nchanged + 1
                end if
            end if
        end if
    end do

100 call ffdev_utils_exit(DEV_ERR,1,'Inproper definition of ''Z'' option for vdw_b0!')

end subroutine change_dih_realm_Z_for_vdw_b0

! ==============================================================================
! subroutine ffdev_parameters_ctrl_bond_r0
! ==============================================================================

subroutine ffdev_parameters_ctrl_bond_r0(fin,exec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: exec
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
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported initff mode '//trim(string)//'!')
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
                call ffdev_utils_exit(DEV_ERR,1,'type must be integer number!')
            end if
            if( (typeid .le. 0) .or. (typeid .gt. nparams) ) then
                call ffdev_utils_exit(DEV_ERR,1,'type is out-of-legal range!')
            end if
            if( params(typeid)%realm .ne. REALM_BOND_D0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'type is not bond_r0!')
            end if
            if( exec ) call ffdev_parameters_bond_r0_init(typeid,mode)
            rst = prmfile_next_line(fin)
        end do
    else
        do i=1,nparams
            if( params(i)%realm .eq. REALM_BOND_D0 ) then
                if( exec ) call ffdev_parameters_bond_r0_init(i,mode)
            end if
        end do
    end if

    call ffdev_parameters_to_tops

 10 format(/,'> Init mode = lowest')
 20 format(/,'> Init mode = average')
 25 format(/,'> Init mode = average  [default]')

end subroutine ffdev_parameters_ctrl_bond_r0

! ==============================================================================
! subroutine ffdev_parameters_ctrl_angle_a0
! ==============================================================================

subroutine ffdev_parameters_ctrl_angle_a0(fin,exec)

    use ffdev_parameters
    use ffdev_parameters_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    logical             :: exec
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
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported initff mode '//trim(string)//'!')
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
                call ffdev_utils_exit(DEV_ERR,1,'type must be integer number!')
            end if
            if( (typeid .le. 0) .or. (typeid .gt. nparams) ) then
                call ffdev_utils_exit(DEV_ERR,1,'type is out-of-legal range!')
            end if
            if( params(typeid)%realm .ne. REALM_ANGLE_A0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'type is not angle_a0!')
            end if
            rst = prmfile_next_line(fin)
        end do
    else
        do i=1,nparams
            if( params(i)%realm .eq. REALM_ANGLE_A0 ) then
                if( exec ) call ffdev_parameters_angle_a0_init(i,mode)
            end if
        end do
    end if

    call ffdev_parameters_to_tops

 10 format(/,'> Init mode = lowest')
 20 format(/,'> Init mode = average')
 25 format(/,'> Init mode = average  [default]')

end subroutine ffdev_parameters_ctrl_angle_a0

! ------------------------------------------------------------------------------

end module ffdev_parameters_control

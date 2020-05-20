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

module ffdev_ffopt_control

use ffdev_sizes
use ffdev_constants
use ffdev_variables

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_ffopt_ctrl_minimize(fin)

    use ffdev_ffopt_dat
    use prmfile
    use ffdev_utils
    use ffdev_ffopt

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    character(80)       :: string
    ! --------------------------------------------------------------------------

    call ffdev_ffopt_set_default

    write(DEV_OUT,'(/,a)') '=== [minimize] ================================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'minimize') ) then
        select case(OptimizationMethod)
            case(MINIMIZATION_STEEPEST_DESCENT)
                write(DEV_OUT,15) 'steepest-descent'
            case(MINIMIZATION_LBFGS)
                write(DEV_OUT,15) 'l-bfgs'
            case(MINIMIZATION_NLOPT)
                write(DEV_OUT,15) 'nlopt'
            case(MINIMIZATION_SHARK)
                write(DEV_OUT,15) 'shark'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unknown minimization method!')
        end select
        write(DEV_OUT,25) NOptSteps
        write(DEV_OUT,75) OutSamples
        write(DEV_OUT,55) MinErrorChange

        select case(OptimizationMethod)
            case(MINIMIZATION_STEEPEST_DESCENT,MINIMIZATION_LBFGS)
                write(DEV_OUT,35) MaxRMSG
                write(DEV_OUT,45) MaxG
                write(DEV_OUT,65) prmfile_onoff(PrintFinalGradient)
            case(MINIMIZATION_NLOPT,MINIMIZATION_SHARK)
                ! nothing to do
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unknown minimization method!')
        end select

        call read_opt_method(fin)
        return
    end if

    if( prmfile_get_string_by_key(fin,'method', string)) then
        select case(string)
            case('steepest-descent')
                OptimizationMethod=MINIMIZATION_STEEPEST_DESCENT
                write(DEV_OUT,10) 'steepest-descent'
            case('l-bfgs')
                OptimizationMethod=MINIMIZATION_LBFGS
                write(DEV_OUT,10) 'l-bfgs'
            case('nlopt')
                OptimizationMethod=MINIMIZATION_NLOPT
                write(DEV_OUT,10) 'nlopt'
            case('shark')
                OptimizationMethod=MINIMIZATION_SHARK
                write(DEV_OUT,10) 'shark'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unknown minimization method!')
        end select
    else
        select case(OptimizationMethod)
            case(MINIMIZATION_STEEPEST_DESCENT)
                write(DEV_OUT,15) 'steepest-descent'
            case(MINIMIZATION_LBFGS)
                write(DEV_OUT,15) 'l-bfgs'
            case(MINIMIZATION_NLOPT)
                write(DEV_OUT,15) 'nlopt'
            case(MINIMIZATION_SHARK)
                write(DEV_OUT,15) 'shark'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unknown minimization method!')
        end select
    end if

    ! minimization criteria
    if( prmfile_get_integer_by_key(fin,'nsteps', NOptSteps)) then
        write(DEV_OUT,20) NOptSteps
        if( NOptSteps .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'steps has to be grater than zero!')
        end if
    else
        write(DEV_OUT,25) NOptSteps
    end if
    if( prmfile_get_integer_by_key(fin,'outsamples', OutSamples)) then
        write(DEV_OUT,70) OutSamples
        if( OutSamples .lt. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'outsamples has to be grater than or equal zero!')
        end if
    else
        write(DEV_OUT,75) OutSamples
    end if
    if( OutSamples .eq. 0 ) OutSamples = -1

    if( prmfile_get_real8_by_key(fin,'minerrorchange', MinErrorChange)) then
        write(DEV_OUT,50) MinErrorChange
    else
        write(DEV_OUT,55) MinErrorChange
    end if

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT,MINIMIZATION_LBFGS)
            if( prmfile_get_real8_by_key(fin,'maxrmsg', MaxRMSG)) then
                write(DEV_OUT,30) MaxRMSG
                if( MaxRMSG .le. 0.0d0 ) then
                    call ffdev_utils_exit(DEV_ERR,1,'maxrmsg has to be grater than zero!')
                end if
            else
                write(DEV_OUT,35) MaxRMSG
            end if

            if( prmfile_get_real8_by_key(fin,'maxg', MaxG)) then
                write(DEV_OUT,40) MaxG
                if( MaxG .le. 0.0d0 ) then
                    call ffdev_utils_exit(DEV_ERR,1,'maxg has to be grater than zero!')
                end if
            else
                write(DEV_OUT,45) MaxG
            end if
            if( prmfile_get_logical_by_key(fin,'printfinalgrad', PrintFinalGradient)) then
                write(DEV_OUT,60) prmfile_onoff(PrintFinalGradient)
            else
                write(DEV_OUT,65) prmfile_onoff(PrintFinalGradient)
            end if
        case(MINIMIZATION_NLOPT,MINIMIZATION_SHARK)
            ! nothing to do
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unknown minimization method!')
    end select

    call read_opt_method(fin)

    return

 10  format ('Minimization method (method)           = ',a16)
 15  format ('Minimization method (method)           = ',a16,'              (default)')
 20  format ('Maximum of minimization steps (nsteps) = ',i12)
 25  format ('Maximum of minimization steps (nsteps) = ',i12,'                  (default)')
 30  format ('Maximal value of RMSDG (maxrmsg)       = ',e16.8)
 35  format ('Maximal value of RMSDG (maxrmsg)       = ',e16.8,'              (default)')
 40  format ('Max value of gradient componenr (maxg) = ',e16.8)
 45  format ('Max value of gradient componenr (maxg) = ',e16.8,'              (default)')
 50  format ('Min energy change (minerrorchange)     = ',e16.8)
 55  format ('Min energy change (minerrorchange)     = ',e16.8,'              (default)')
 60  format ('Print final gradient (printfinalgrad)  = ',a12)
 65  format ('Print final gradient (printfinalgrad)  = ',a12,'                  (default)')
 70  format ('Opt summary samples (outsamples)       = ',i12)
 75  format ('Opt summary samples (outsamples)       = ',i12,'                  (default)')

end subroutine ffdev_ffopt_ctrl_minimize

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_opt_method(fin)

    use prmfile
    use ffdev_ffopt_dat

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT)
            call read_sd_method(fin)
        case(MINIMIZATION_LBFGS)
            call read_lbfgs_method(fin)
        case(MINIMIZATION_NLOPT)
            call read_nlopt_method(fin)
        case(MINIMIZATION_SHARK)
            call read_shark_method(fin)
    end select

end subroutine read_opt_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_sd_method(fin)

    use ffdev_ffopt_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [steepest-descent] ========================================================='

    if( .not. prmfile_open_section(fin,'steepest-descent') ) then
        write(DEV_OUT,15) InitialStepSize
        write(DEV_OUT,25) MaximalStepSize
        write(DEV_OUT,35) AcceptRatio
        write(DEV_OUT,45) RejectRatio
        write(DEV_OUT,55) prmfile_onoff(AdaptiveStep)
        return
    end if

    if( prmfile_get_real8_by_key(fin,'initialstepsize', InitialStepSize)) then
        write(DEV_OUT,10) InitialStepSize
        if( InitialStepSize .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'initialstepsize has to be grater than zero!')
        end if
    else
        write(DEV_OUT,15) InitialStepSize
    end if

    if( prmfile_get_real8_by_key(fin,'maximalstepsize', MaximalStepSize)) then
        write(DEV_OUT,20) MaximalStepSize
        if( MaximalStepSize .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'maximalstepsize has to be grater than zero!')
        end if
    else
        write(DEV_OUT,25) MaximalStepSize
    end if

    if( prmfile_get_real8_by_key(fin,'acceptratio', AcceptRatio)) then
        write(DEV_OUT,30) AcceptRatio
        if( AcceptRatio .lt. 1.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'acceptratio has to be grater or equal to 1.0!')
        end if
    else
        write(DEV_OUT,35) AcceptRatio
    end if

    if( prmfile_get_real8_by_key(fin,'rejectratio', RejectRatio)) then
        write(DEV_OUT,40) RejectRatio
        if( RejectRatio .le. 0.0d0 .or. RejectRatio .ge. 1.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'rejectratio has to be in (0;1) interval!')
        end if
    else
        write(DEV_OUT,45) RejectRatio
    end if

    if( prmfile_get_logical_by_key(fin,'adaptivestep', AdaptiveStep)) then
        write(DEV_OUT,50) prmfile_onoff(AdaptiveStep)
    else
        write(DEV_OUT,55) prmfile_onoff(AdaptiveStep)
    end if

 10  format ('Initial step size (initialstepsize)    = ',f21.8)
 15  format ('Initial step size (initialstepsize)    = ',f21.8,'         (default)')
 20  format ('Maximal step size (maximalstepsize)    = ',f21.8)
 25  format ('Maximal step size (maximalstepsize)    = ',f21.8,'         (default)')
 30  format ('Accept ratio (acceptratio)             = ',f16.3)
 35  format ('Accept ratio (acceptratio)             = ',f16.3,'              (default)')
 40  format ('Reject ratio (rejectratio)             = ',f16.3)
 45  format ('Reject ratio (rejectratio)             = ',f16.3,'              (default)')
 50  format ('Adaptive step (adaptivestep)           = ',a12)
 55  format ('Adaptive step (adaptivestep)           = ',a12,'                  (default)')

 return

end subroutine read_sd_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_lbfgs_method(fin)

    use prmfile
    use ffdev_ffopt_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [l-bfgs] ==================================================================='

    if( .not. prmfile_open_section(fin,'l-bfgs') ) then
        write(DEV_OUT,15) NumberOfCorrections
        return
    end if

    if( prmfile_get_integer_by_key(fin,'numofcorrections', NumberOfCorrections)) then
        write(DEV_OUT,10) NumberOfCorrections
        if( NumberOfCorrections .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Number of corrections (numofcorrections) has to be grater than zero!')
        end if
    else
        write(DEV_OUT,15) NumberOfCorrections
    end if

    return

 10  format ('Num of corrections (numofcorrections)  = ',i12)
 15  format ('Num of corrections (numofcorrections)  = ',i12,'                  (default)')

end subroutine read_lbfgs_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_nlopt_method(fin)

    use prmfile
    use ffdev_ffopt_dat
    use ffdev_utils

    implicit none
    include 'nlopt.f'

    type(PRMFILE_TYPE)          :: fin
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [nlopt] ===================================================================='

    if( .not. prmfile_open_section(fin,'nlopt') ) then
        select case(NLOpt_Method)
            case(NLOPT_LN_COBYLA)
                write(DEV_OUT,25) 'NLOPT_LN_COBYLA'
            case(NLOPT_LN_BOBYQA)
                write(DEV_OUT,25) 'NLOPT_LN_BOBYQA'
            case(NLOPT_LN_NEWUOA)
                write(DEV_OUT,25) 'NLOPT_LN_NEWUOA'
            case(NLOPT_LN_NELDERMEAD)
                write(DEV_OUT,25) 'NLOPT_LN_NELDERMEAD'
            case(NLOPT_LN_SBPLX)
                write(DEV_OUT,25) 'NLOPT_LN_SBPLX'
            case(NLOPT_GN_DIRECT)
                write(DEV_OUT,25) 'NLOPT_GN_DIRECT'
            case(NLOPT_GN_DIRECT_L)
                write(DEV_OUT,25) 'NLOPT_GN_DIRECT_L'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported methopd in read_nlopt_method!')
        end select
        write(DEV_OUT,35) NLOpt_InitialStep
        return
    end if

    if( prmfile_get_string_by_key(fin,'algorithm', string)) then
        select case(trim(string))
            case('NLOPT_LN_COBYLA')
                NLOpt_Method = NLOPT_LN_COBYLA
                write(DEV_OUT,20) trim(string)
            case('NLOPT_LN_BOBYQA')
                NLOpt_Method = NLOPT_LN_BOBYQA
                write(DEV_OUT,20) trim(string)
            case('NLOPT_LN_NEWUOA')
                NLOpt_Method = NLOPT_LN_NEWUOA
                write(DEV_OUT,20) trim(string)
            case('NLOPT_LN_NELDERMEAD')
                NLOpt_Method = NLOPT_LN_NELDERMEAD
                write(DEV_OUT,20) trim(string)
            case('NLOPT_LN_SBPLX')
                NLOpt_Method = NLOPT_LN_SBPLX
                write(DEV_OUT,20) trim(string)
            case('NLOPT_GN_DIRECT')
                NLOpt_Method = NLOPT_GN_DIRECT
                write(DEV_OUT,20) trim(string)
            case('NLOPT_GN_DIRECT_L')
                NLOpt_Method = NLOPT_GN_DIRECT_L
                write(DEV_OUT,20) trim(string)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported methopd in read_nlopt_method!')
        end select
    else
        select case(NLOpt_Method)
            case(NLOPT_LN_COBYLA)
                write(DEV_OUT,25) 'NLOPT_LN_COBYLA'
            case(NLOPT_LN_BOBYQA)
                write(DEV_OUT,25) 'NLOPT_LN_BOBYQA'
            case(NLOPT_LN_NEWUOA)
                write(DEV_OUT,25) 'NLOPT_LN_NEWUOA'
            case(NLOPT_LN_NELDERMEAD)
                write(DEV_OUT,25) 'NLOPT_LN_NELDERMEAD'
            case(NLOPT_LN_SBPLX)
                write(DEV_OUT,25) 'NLOPT_LN_SBPLX'
            case(NLOPT_GN_DIRECT)
                write(DEV_OUT,25) 'NLOPT_GN_DIRECT'
            case(NLOPT_GN_DIRECT_L)
                write(DEV_OUT,25) 'NLOPT_GN_DIRECT_L'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported methopd in read_nlopt_method!')
        end select
    end if

    if( prmfile_get_real8_by_key(fin,'initialstep', NLOpt_InitialStep)) then
        write(DEV_OUT,30) NLOpt_InitialStep
    else
        write(DEV_OUT,35) NLOpt_InitialStep
    end if

    return
 20  format ('Optmization algorithm (algorithm)      = ',A)
 25  format ('Optmization algorithm (algorithm)      = ',A20,'          (default)')
 30  format ('Initial step (initialstep)             = ',f12.7)
 35  format ('Initial step (initialstep)             = ',f12.7,'                  (default)')

end subroutine read_nlopt_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_shark_method(fin)

    use prmfile
    use ffdev_ffopt_dat
    use ffdev_utils

    implicit none
    include 'nlopt.f'

    type(PRMFILE_TYPE)          :: fin
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [shark] ===================================================================='

    if( .not. prmfile_open_section(fin,'shark') ) then
        select case(Shark_Method)
            case(SHARK_CMA_ES)
                write(DEV_OUT,25) 'CMA-ES'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported methopd in read_shark_method!')
        end select
        write(DEV_OUT,35) Shark_InitialStep
        write(DEV_OUT,55) prmfile_onoff(Shark_EnableBoxing)
        write(DEV_OUT,65) Shark_NRuns
        select case(Shark_ParameterGuess)
            case(SHARK_GUESS_INPUT)
                write(DEV_OUT,75) 'input'
            case(SHARK_GUESS_RANDOMIZE)
                write(DEV_OUT,75) 'randomize'
            case(SHARK_GUESS_KEEP)
                write(DEV_OUT,75) 'keep'
            case(SHARK_GUESS_MIX)
                write(DEV_OUT,75) 'mix'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported guess in read_shark_method!')
        end select
        return
    end if

    if( prmfile_get_string_by_key(fin,'algorithm', string)) then
        select case(trim(string))
            case('CMA-ES')
                Shark_Method = SHARK_CMA_ES
                write(DEV_OUT,20) trim(string)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported method in read_shark_method!')
        end select
    else
        select case(Shark_Method)
            case(SHARK_CMA_ES)
                write(DEV_OUT,25) 'CMA-ES'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported method in read_shark_method!')
        end select
    end if

    if( prmfile_get_real8_by_key(fin,'initialstep', Shark_InitialStep)) then
        write(DEV_OUT,30) Shark_InitialStep
    else
        write(DEV_OUT,35) Shark_InitialStep
    end if

    if( prmfile_get_logical_by_key(fin,'boxing', Shark_EnableBoxing)) then
        write(DEV_OUT,50) prmfile_onoff(Shark_EnableBoxing)
    else
        write(DEV_OUT,55) prmfile_onoff(Shark_EnableBoxing)
    end if

    if( prmfile_get_integer_by_key(fin,'nruns', Shark_NRuns)) then
        write(DEV_OUT,60) Shark_NRuns
    else
        write(DEV_OUT,65) Shark_NRuns
    end if

    if( prmfile_get_string_by_key(fin,'guess', string)) then
        select case(trim(string))
            case('input')
                Shark_ParameterGuess = SHARK_GUESS_INPUT
                write(DEV_OUT,70) trim(string)
            case('randomize')
                Shark_ParameterGuess = SHARK_GUESS_RANDOMIZE
                write(DEV_OUT,70) trim(string)
            case('keep')
                Shark_ParameterGuess = SHARK_GUESS_KEEP
                write(DEV_OUT,70) trim(string)
            case('mix')
                Shark_ParameterGuess = SHARK_GUESS_MIX
                write(DEV_OUT,70) trim(string)
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported guess in read_shark_method!')
        end select
    else
        select case(Shark_ParameterGuess)
            case(SHARK_GUESS_INPUT)
                write(DEV_OUT,75) 'input'
            case(SHARK_GUESS_RANDOMIZE)
                write(DEV_OUT,75) 'randomize'
            case(SHARK_GUESS_KEEP)
                write(DEV_OUT,75) 'keep'
            case(SHARK_GUESS_MIX)
                write(DEV_OUT,75) 'mix'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unsupported guess in read_shark_method!')
        end select
    end if

    return

 20  format ('Optmization algorithm (algorithm)      = ',A)
 25  format ('Optmization algorithm (algorithm)      = ',A20,'          (default)')
 30  format ('Initial step (initialstep)             = ',f12.7)
 35  format ('Initial step (initialstep)             = ',f12.7,'                  (default)')
 50  format ('Enable boxing (boxing)                 = ',a12)
 55  format ('Enable boxing (boxing)                 = ',a12,'                  (default)')
 60  format ('Number of runs (nruns)                 = ',i12)
 65  format ('Number of runs (nruns)                 = ',i12,'                  (default)')
 70  format ('Initial guess parameters (guess)       = ',a12)
 75  format ('Initial guess parameters (guess)       = ',a12,'                  (default)')

end subroutine read_shark_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_ffopt_control

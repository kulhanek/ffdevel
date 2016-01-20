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
            case(MINIMIZATION_SA)
                write(DEV_OUT,15) 'simulated annealing (sa)'
        end select
        write(DEV_OUT,25) NOptSteps
        write(DEV_OUT,35) MaxRMSG
        write(DEV_OUT,45) MaxG
        write(DEV_OUT,55) MinErrorChange
        write(DEV_OUT,65) prmfile_onoff(PrintFinalGradient)
        write(DEV_OUT,75) OutSamples

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
            case('sa')
                OptimizationMethod=MINIMIZATION_SA
                write(DEV_OUT,10) 'simulated annealing (sa)'
            case default
                call ffdev_utils_exit(DEV_OUT,1,'Unknown minimization method!')
        end select
    else
        select case(OptimizationMethod)
            case(MINIMIZATION_STEEPEST_DESCENT)
                write(DEV_OUT,15) 'steepest-descent'
            case(MINIMIZATION_LBFGS)
                write(DEV_OUT,15) 'l-bfgs'
            case(MINIMIZATION_NLOPT)
                write(DEV_OUT,15) 'nlopt'
            case(MINIMIZATION_SA)
                write(DEV_OUT,15) 'simulated annealing (sa)'
        end select
    end if

    ! minimization criteria
    if( prmfile_get_integer_by_key(fin,'steps', NOptSteps)) then
        write(DEV_OUT,20) NOptSteps
        if( NOptSteps .le. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'steps has to be grater than zero!')
        end if
    else
        write(DEV_OUT,25) NOptSteps
    end if

    if( prmfile_get_real8_by_key(fin,'maxrmsg', MaxRMSG)) then
        write(DEV_OUT,30) MaxRMSG
        if( MaxRMSG .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'maxrmsg has to be grater than zero!')
        end if
    else
        write(DEV_OUT,35) MaxRMSG
    end if

    if( prmfile_get_real8_by_key(fin,'maxg', MaxG)) then
        write(DEV_OUT,40) MaxG
        if( MaxG .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'maxg has to be grater than zero!')
        end if
    else
        write(DEV_OUT,45) MaxG
    end if

    if( prmfile_get_real8_by_key(fin,'minerrorchange', MinErrorChange)) then
        write(DEV_OUT,50) MinErrorChange
    else
        write(DEV_OUT,55) MinErrorChange
    end if

    if( prmfile_get_logical_by_key(fin,'printfinalgrad', PrintFinalGradient)) then
        write(DEV_OUT,60) prmfile_onoff(PrintFinalGradient)
    else
        write(DEV_OUT,65) prmfile_onoff(PrintFinalGradient)
    end if

    if( prmfile_get_integer_by_key(fin,'outsamples', OutSamples)) then
        write(DEV_OUT,70) OutSamples
        if( OutSamples .lt. 0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'outsamples has to be grater than or equal zero!')
        end if
    else
        write(DEV_OUT,75) OutSamples
    end if
    if( OutSamples .eq. 0 ) OutSamples = -1

    call read_opt_method(fin)

    return

 10  format ('Minimization method (method)           = ',a16)
 15  format ('Minimization method (method)           = ',a16,'              (default)')
 20  format ('Maximum of minimization steps (steps)  = ',i12)
 25  format ('Maximum of minimization steps (steps)  = ',i12,'                  (default)')
 30  format ('Maximal value of RMSDG (maxrmsg)       = ',f16.3)
 35  format ('Maximal value of RMSDG (maxrmsg)       = ',f16.3,'              (default)')
 40  format ('Max value of gradient componenr (maxg) = ',f16.3)
 45  format ('Max value of gradient componenr (maxg) = ',f16.3,'              (default)')
 50  format ('Min energy change (minerrorchange)     = ',f16.3)
 55  format ('Min energy change (minerrorchange)     = ',f16.3,'              (default)')
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
        case(MINIMIZATION_SA)
            call read_sa_method(fin)
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
            call ffdev_utils_exit(DEV_OUT,1,'initialstepsize has to be grater than zero!')
        end if
    else
        write(DEV_OUT,15) InitialStepSize
    end if

    if( prmfile_get_real8_by_key(fin,'maximalstepsize', MaximalStepSize)) then
        write(DEV_OUT,20) MaximalStepSize
        if( MaximalStepSize .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'maximalstepsize has to be grater than zero!')
        end if
    else
        write(DEV_OUT,25) MaximalStepSize
    end if

    if( prmfile_get_real8_by_key(fin,'acceptratio', AcceptRatio)) then
        write(DEV_OUT,30) AcceptRatio
        if( AcceptRatio .lt. 1.0d0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'acceptratio has to be grater or equal to 1.0!')
        end if
    else
        write(DEV_OUT,35) AcceptRatio
    end if

    if( prmfile_get_real8_by_key(fin,'rejectratio', RejectRatio)) then
        write(DEV_OUT,40) RejectRatio
        if( RejectRatio .le. 0.0d0 .or. RejectRatio .ge. 1.0d0 ) then
            call ffdev_utils_exit(DEV_OUT,1,'rejectratio has to be in (0;1) interval!')
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
            call ffdev_utils_exit(DEV_OUT,1,'Number of corrections (numofcorrections) has to be grater than zero!')
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

subroutine read_sa_method(fin)

    use ffdev_utils
    use prmfile
    use ffdev_ffopt_dat

    implicit none
    type(PRMFILE_TYPE)             :: fin
    ! -----------------------------------------------
    logical                        :: read_some
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [sa] ======================================================================='

    ! is set active section?
    if( .not. prmfile_open_section(fin,'sa') ) then
        call ffdev_utils_exit(DEV_OUT,1,'>>> ERROR: [sa] section is not current section!')
    end if

    read_some = .false.

    if( prmfile_get_real8_by_key(fin,'temp', OptSA_Temp) ) then
        write(DEV_OUT,20) OptSA_Temp
        read_some = .true.
    end if

    if( prmfile_get_real8_by_key(fin,'rt', OptSA_RT) ) then
        write(DEV_OUT,30) OptSA_RT
        read_some = .true.
    end if

    if( prmfile_get_integer_by_key(fin,'ns', OptSA_NS) ) then
        write(DEV_OUT,40) OptSA_NS
        read_some = .true.
    end if

    if( prmfile_get_integer_by_key(fin,'nt', OptSA_NT) ) then
        write(DEV_OUT,50) OptSA_NT
        read_some = .true.
    end if

    if( prmfile_get_integer_by_key(fin,'maxevl', OptSA_MAXEVL) ) then
        write(DEV_OUT,60) OptSA_MAXEVL
        read_some = .true.
    end if

    if( prmfile_get_integer_by_key(fin,'iseed', OptSA_ISEED) ) then
        write(DEV_OUT,70) OptSA_ISEED
        read_some = .true.
    end if

    if( .not. read_some ) then
        write(DEV_OUT,5)
    end if

  5 format('>>> Nothing changed for simmulated annealing method ...')
 10 format('Direct chi optimization         :',A10)
 20 format('Initial temperature             :',F18.7)
 30 format('Temperature reduction factor    :',F18.7)
 40 format('Number of cycles                :',I10)
 50 format('Number of iteration to reduce T :',I10)
 60 format('Max number of function eval     :',I10)
 70 format('Random generator seed           :',I10)

end subroutine read_sa_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_ffopt_control

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

module ffdev_geoopt_control

use ffdev_sizes
use ffdev_constants
use ffdev_variables

implicit none
contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_geoopt_ctrl_minimize(fin)

    use ffdev_geoopt_dat
    use ffdev_geoopt
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    character(80)       :: string
    ! --------------------------------------------------------------------------

    call ffdev_geoopt_set_default

    write(DEV_OUT,'(/,a)') '=== [minimize] ================================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'minimize') ) then
        select case(OptimizationMethod)
            case(MINIMIZATION_STEEPEST_DESCENT)
                write(DEV_OUT,15) 'steepest-descent'
            case(MINIMIZATION_LBFGS)
                write(DEV_OUT,15) 'l-bfgs'
        end select
        write(DEV_OUT,25) NOptSteps
        write(DEV_OUT,35) MaxRMSG
        write(DEV_OUT,45) MaxG
        write(DEV_OUT,55) MinEnergyChange
        write(DEV_OUT,65) prmfile_onoff(PrintFinalGradient)
        write(DEV_OUT,75) OutSamples
        write(DEV_OUT,85) TrajSamples

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
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unknown minimization method!')
        end select
    else
        select case(OptimizationMethod)
            case(MINIMIZATION_STEEPEST_DESCENT)
                write(DEV_OUT,15) 'steepest-descent'
            case(MINIMIZATION_LBFGS)
                write(DEV_OUT,15) 'l-bfgs'
        end select
    end if

    ! minimization criteria
    if( prmfile_get_integer_by_key(fin,'steps', NOptSteps)) then
        write(DEV_OUT,20) NOptSteps
        if( NOptSteps .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'steps has to be grater than zero!')
        end if
    else
        write(DEV_OUT,25) NOptSteps
    end if

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

    if( prmfile_get_real8_by_key(fin,'minenergychange', MinEnergyChange)) then
        write(DEV_OUT,50) MinEnergyChange
    else
        write(DEV_OUT,55) MinEnergyChange
    end if

    if( prmfile_get_logical_by_key(fin,'printfinalgrad', PrintFinalGradient)) then
        write(DEV_OUT,60) prmfile_onoff(PrintFinalGradient)
    else
        write(DEV_OUT,65) prmfile_onoff(PrintFinalGradient)
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

    if( prmfile_get_integer_by_key(fin,'trajsamples', TrajSamples)) then
        write(DEV_OUT,80) TrajSamples
        if( TrajSamples .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'trajsamples has to be grater than or equal zero!')
        end if
    else
        write(DEV_OUT,85) TrajSamples
    end if
    if( TrajSamples .eq. 0 ) TrajSamples = -1

    call read_opt_method(fin)

    return

 10  format ('Energy minimization method (method)    = ',a16)
 15  format ('Energy minimization method (method)    = ',a16,'              (default)')
 20  format ('Maximum of minimization steps (steps)  = ',i12)
 25  format ('Maximum of minimization steps (steps)  = ',i12,'                  (default)')
 30  format ('Maximal value of RMSDG (maxrmsg)       = ',e16.8)
 35  format ('Maximal value of RMSDG (maxrmsg)       = ',e16.8,'              (default)')
 40  format ('Max value of gradient componenr (maxg) = ',e16.8)
 45  format ('Max value of gradient componenr (maxg) = ',e16.8,'              (default)')
 50  format ('Min energy change (minenergychange)    = ',e16.8)
 55  format ('Min energy change (minenergychange)    = ',e16.8,'              (default)')
 60  format ('Print final gradient (printfinalgrad)  = ',a12)
 65  format ('Print final gradient (printfinalgrad)  = ',a12,'                  (default)')
 70  format ('Energy summary samples (outsamples)    = ',i12)
 75  format ('Energy summary samples (outsamples)    = ',i12,'                  (default)')
 80  format ('Trajectory samples (trajsamples)       = ',i12)
 85  format ('Trajectory samples (trajsamples)       = ',i12,'                  (default)')

end subroutine ffdev_geoopt_ctrl_minimize

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_opt_method(fin)

    use prmfile
    use ffdev_geoopt_dat

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT)
            call read_sd_method(fin)
        case(MINIMIZATION_LBFGS)
            call read_lbfgs_method(fin)
    end select

end subroutine read_opt_method

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine read_sd_method(fin)

    use ffdev_geoopt_dat
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
    use ffdev_geoopt_dat
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

subroutine ffdev_geoopt_ctrl_files(fin)

    use prmfile
    use ffdev_geoopt_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! -----------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [files] ===================================================================='

    if(.not. prmfile_open_section(fin,'files')) then
        if( len(OptTopName) .eq. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Topology (topology) is not specified!')
        end if
        write (DEV_OUT,10) trim(OptTopName)
        if( len(OptCrdName) .eq. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Input coordinates (input) are required!')
        end if
        write (DEV_OUT,20) trim(OptCrdName)
        write (DEV_OUT,35) trim(OptRstName)
        write (DEV_OUT,45) trim(OptTrajName)
        return
    end if

    ! topology file
    if(.not. prmfile_get_string_by_key(fin,'topology', OptTopName)) then
        call ffdev_utils_exit(DEV_ERR,1,'Topology (topology) is not specified!')
    end if
    write (DEV_OUT,10) trim(OptTopName)

    ! input file
    if(.not. prmfile_get_string_by_key(fin,'input', OptCrdName)) then
        call ffdev_utils_exit(DEV_ERR,1,'Input coordinates (input) are required!')
    end if
    write (DEV_OUT,20) trim(OptCrdName)

    ! final coordinates
    if( .not. prmfile_get_string_by_key(fin,'final', OptRstName) ) then
        write (DEV_OUT,35) trim(OptRstName)
    else
        write (DEV_OUT,30) trim(OptRstName)
    end if

    ! trajectory file
    if( .not. prmfile_get_string_by_key(fin,'trajectory', OptTrajName)) then
        write (DEV_OUT,45) trim(OptTrajName)
    else
        write (DEV_OUT,40) trim(OptTrajName)
    end if

    return

10 format('Topology file (topology)               = ',a)
20 format('Initial coordinate file (input)        = ',a)
30 format('Final coordinate file (final)          = ',a)
35 format('Final coordinate file (final)          = ',a12,'                  (default)')
40 format('Trajectory file (trajectory)           = ',a)
45 format('Trajectory file (trajectory)           = ',a12,'                  (default)')

end subroutine ffdev_geoopt_ctrl_files

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_geoopt_control

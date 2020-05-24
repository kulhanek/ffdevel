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

module ffdev_genpoints_control

use ffdev_sizes
use ffdev_constants
use ffdev_variables

implicit none
contains

!===============================================================================
! subroutine ffdev_genpoints_ctrl_files
!===============================================================================

subroutine ffdev_genpoints_ctrl_files(fin)

    use prmfile
    use ffdev_genpoints_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [files] ===================================================================='

    if(.not. prmfile_open_section(fin,'files')) then
        call ffdev_utils_exit(DEV_ERR,1,'[files] section not found.')
    end if

    ! topology file
    if(.not. prmfile_get_string_by_key(fin,'topology', GenTopName)) then
        call ffdev_utils_exit(DEV_ERR,1,'Topology (topology) is not specified!')
    end if
    write (DEV_OUT,10) trim(GenTopName)

    ! input file
    if(.not. prmfile_get_string_by_key(fin,'input', GenCrdName)) then
        call ffdev_utils_exit(DEV_ERR,1,'Input coordinates (input) are required!')
    end if
    write (DEV_OUT,20) trim(GenCrdName)

    ! input rotors
    if( prmfile_get_string_by_key(fin,'rotors', GenRotName)) then
        write (DEV_OUT,30) trim(GenRotName)
        RotorListLoaded = .true.
    end if

    ! output points
    if(.not. prmfile_get_string_by_key(fin,'output', GenOutName)) then
        write (DEV_OUT,45) trim(GenOutName)
    else
        write (DEV_OUT,40) trim(GenOutName)
    end if

    ! global point
    if(.not. prmfile_get_string_by_key(fin,'global', GenFinName)) then
        write (DEV_OUT,55) trim(GenFinName)
    else
        write (DEV_OUT,50) trim(GenFinName)
    end if

    ! profile
    if(.not. prmfile_get_string_by_key(fin,'profile', ProfileName)) then
        write (DEV_OUT,65) trim(ProfileName)
    else
        write (DEV_OUT,60) trim(ProfileName)
    end if

    return

10 format('Topology file (topology)               = ',a)
20 format('Initial coordinate file (input)        = ',a)
30 format('Rotor bonds list (rotors)              = ',a)
40 format('Output points (output)                 = ',a)
45 format('Output points (output)                 = ',a12,'                  (default)')
50 format('Global point (global)                  = ',a)
55 format('Global point (global)                  = ',a12,'                  (default)')
60 format('Energy profile (profile)               = ',a)
65 format('Energy profile (profile)               = ',a12,'                  (default)')

end subroutine ffdev_genpoints_ctrl_files

!===============================================================================
! subroutine ffdev_genpoints_ctrl_points
!===============================================================================

subroutine ffdev_genpoints_ctrl_points(fin)

    use ffdev_genpoints_dat
    use prmfile
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------
    character(80)       :: string
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [points] ==================================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'points') ) then
        select case(GeneratorMethod)
            case(GENPOINTS_SYSTEMATIC_METHOD)
                write(DEV_OUT,15) 'systematic'
            case(GENPOINTS_STOCHASTIC_METHOD)
                write(DEV_OUT,15) 'stochastic'
            case(GENPOINTS_STOCHASTIC_OPT_METHOD)
                write(DEV_OUT,15) 'stochastic-opt'
            case(GENPOINTS_STEPBYSTEP_METHOD)
                write(DEV_OUT,15) 'stepbystep'
            case(GENPOINTS_STOCHASTIC_BY_STEPS_METHOD)
                write(DEV_OUT,15) 'stochasticbysteps'
            case(GENPOINTS_NMODES_METHOD)
                write(DEV_OUT,15) 'nmodes'
            case(GENPOINTS_SYSTEMATIC_OPT_METHOD)
                write(DEV_OUT,15) 'systematic-opt'
        end select
        write(DEV_OUT,25) MaxPoints
        write(DEV_OUT,35) MaxEnergy
        write(DEV_OUT,45) OptimizePoints
        write(DEV_OUT,55) HoldCV

        call read_genpoints_method(fin)
        return
    end if

    if( prmfile_get_string_by_key(fin,'method', string)) then
        select case(string)
            case('systematic')
                GeneratorMethod=GENPOINTS_SYSTEMATIC_METHOD
                write(DEV_OUT,10) 'systematic'
            case('stochastic')
                GeneratorMethod=GENPOINTS_STOCHASTIC_METHOD
                write(DEV_OUT,10) 'stochastic'
            case('stochastic-opt')
                GeneratorMethod=GENPOINTS_STOCHASTIC_OPT_METHOD
                write(DEV_OUT,10) 'stochastic-opt'
            case('stepbystep')
                GeneratorMethod=GENPOINTS_STEPBYSTEP_METHOD
                write(DEV_OUT,10) 'stepbystep'
            case('stochasticbysteps')
                GeneratorMethod=GENPOINTS_STOCHASTIC_BY_STEPS_METHOD
                write(DEV_OUT,10) 'stochasticbysteps'
            case('nmodes')
                GeneratorMethod=GENPOINTS_NMODES_METHOD
                write(DEV_OUT,10) 'nmodes'
            case('systematic-opt')
                GeneratorMethod=GENPOINTS_SYSTEMATIC_OPT_METHOD
                write(DEV_OUT,10) 'systematic-opt'
            case default
                call ffdev_utils_exit(DEV_ERR,1,'Unknown generator method!')
        end select
    else
        select case(GeneratorMethod)
            case(GENPOINTS_SYSTEMATIC_METHOD)
                write(DEV_OUT,15) 'systematic'
            case(GENPOINTS_STOCHASTIC_METHOD)
                write(DEV_OUT,15) 'stochastic'
            case(GENPOINTS_STOCHASTIC_OPT_METHOD)
                write(DEV_OUT,15) 'stochastic-opt'
            case(GENPOINTS_STEPBYSTEP_METHOD)
                write(DEV_OUT,15) 'stepbystep'
            case(GENPOINTS_STOCHASTIC_BY_STEPS_METHOD)
                write(DEV_OUT,15) 'stochasticbysteps'
            case(GENPOINTS_NMODES_METHOD)
                write(DEV_OUT,15) 'nmodes'
            case(GENPOINTS_SYSTEMATIC_OPT_METHOD)
                write(DEV_OUT,15) 'systematic-opt'
        end select
    end if

    if( prmfile_get_integer_by_key(fin,'max_points', MaxPoints)) then
        write(DEV_OUT,20) MaxPoints
        if( MaxPoints .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'npoints has to be grater than zero!')
        end if
    else
        write(DEV_OUT,25) MaxPoints
    end if

    if( prmfile_get_real8_by_key(fin,'max_energy', MaxEnergy)) then
        write(DEV_OUT,30) MaxEnergy
        if( MaxEnergy .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'MaxEnergy has to be grater than zero!')
        end if
    else
        write(DEV_OUT,35) MaxEnergy
    end if

    if( prmfile_get_logical_by_key(fin,'optimize', OptimizePoints)) then
        write(DEV_OUT,40) prmfile_onoff(OptimizePoints)
    else
        write(DEV_OUT,45) prmfile_onoff(OptimizePoints)
    end if

    if( prmfile_get_logical_by_key(fin,'holdcv', HoldCV)) then
        write(DEV_OUT,50) prmfile_onoff(HoldCV)
    else
        write(DEV_OUT,55) prmfile_onoff(HoldCV)
    end if

    call read_genpoints_method(fin)

    return

 10  format ('Point generator method                 = ',a16)
 15  format ('Point generator method                 = ',a16,'              (default)')
 20  format ('Maximum number of points (max_points)  = ',i12)
 25  format ('Maximum number of points (max_points)  = ',i12,'                  (default)')
 30  format ('Maximum energy (max_energy)            = ',f16.3)
 35  format ('Maximum energy (max_energy)            = ',f16.3,'              (default)')
 40  format ('Optimize points (optimize)             = ',a)
 45  format ('Optimize points (optimize)             = ',a16,'              (default)')
 50  format ('Keep rotor CV constant (holdcv)        = ',a)
 55  format ('Keep rotor CV constant (holdcv)        = ',a16,'              (default)')

end subroutine ffdev_genpoints_ctrl_points

!===============================================================================
! subroutine read_genpoints_method
!===============================================================================

subroutine read_genpoints_method(fin)

    use prmfile
    use ffdev_genpoints_dat

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    select case(GeneratorMethod)
        case(GENPOINTS_SYSTEMATIC_METHOD)
            call read_systematic_method(fin)
!        case(GENPOINTS_STOCHASTIC_METHOD)
!            call read_stochastic_method(fin)
!        case(GENPOINTS_STOCHASTIC_OPT_METHOD)
!            call read_stochastic_opt__method(fin)
        case(GENPOINTS_STEPBYSTEP_METHOD)
            call read_stepbystep_method(fin)
        case(GENPOINTS_STOCHASTIC_BY_STEPS_METHOD)
            call read_stochasticbystep_method(fin)
        case(GENPOINTS_NMODES_METHOD)
            call read_nmodes_method(fin)
    end select

end subroutine read_genpoints_method

!===============================================================================
! subroutine read_systematic_method
!===============================================================================

subroutine read_systematic_method(fin)

    use prmfile
    use ffdev_genpoints_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [systematic] ==============================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'systematic') ) then
        write(DEV_OUT,35) TickAngle*DEV_R2D
        return
    end if

    if( prmfile_get_real8_by_key(fin,'tick', TickAngle)) then
        write(DEV_OUT,30) TickAngle
        if( MaxEnergy .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'tick has to be grater than zero!')
        end if
        TickAngle = TickAngle*DEV_D2R
    else
        write(DEV_OUT,35) TickAngle*DEV_R2D
    end if

 30  format ('Dihedral angle tick (tick)             = ',f16.3)
 35  format ('Dihedral angle tick (tick)             = ',f16.3,'              (default)')

end subroutine read_systematic_method

!===============================================================================
! subroutine read_stepbystep_method
!===============================================================================

subroutine read_stepbystep_method(fin)

    use prmfile
    use ffdev_genpoints_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [stepbystep] ==============================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'stepbystep') ) then
        write(DEV_OUT,35) TickAngle*DEV_R2D
        return
    end if

    if( prmfile_get_real8_by_key(fin,'tick', TickAngle)) then
        write(DEV_OUT,30) TickAngle
        if( MaxEnergy .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'tick has to be grater than zero!')
        end if
        TickAngle = TickAngle*DEV_D2R
    else
        write(DEV_OUT,35) TickAngle*DEV_R2D
    end if

 30  format ('Dihedral angle tick (tick)             = ',f16.3)
 35  format ('Dihedral angle tick (tick)             = ',f16.3,'              (default)')

end subroutine read_stepbystep_method

!===============================================================================
! subroutine read_stochasticbystep_method
!===============================================================================

subroutine read_stochasticbystep_method(fin)

    use prmfile
    use ffdev_genpoints_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [stochasticbysteps] ========================================================'

    ! open first section
    if( .not. prmfile_open_section(fin,'stochasticbysteps') ) then
        write(DEV_OUT,35) TickAngle*DEV_R2D
        return
    end if

    if( prmfile_get_real8_by_key(fin,'tick', TickAngle)) then
        write(DEV_OUT,30) TickAngle
        if( MaxEnergy .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'tick has to be grater than zero!')
        end if
        TickAngle = TickAngle*DEV_D2R
    else
        write(DEV_OUT,35) TickAngle*DEV_R2D
    end if

 30  format ('Dihedral angle tick (tick)             = ',f16.3)
 35  format ('Dihedral angle tick (tick)             = ',f16.3,'              (default)')

end subroutine read_stochasticbystep_method

!===============================================================================
! subroutine read_nmodes_method
!===============================================================================

subroutine read_nmodes_method(fin)

    use prmfile
    use ffdev_genpoints_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [nmodes] ==================================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'nmodes') ) then
        write(DEV_OUT,35) NModeAmplitude
        write(DEV_OUT,45) NModePoints
        return
    end if

    if( prmfile_get_real8_by_key(fin,'amplitude', NModeAmplitude)) then
        write(DEV_OUT,30) NModeAmplitude
    else
        write(DEV_OUT,35) NModeAmplitude
    end if

    if( prmfile_get_integer_by_key(fin,'npoints', NModePoints)) then
        write(DEV_OUT,40) NModePoints
        if( NModePoints .le. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'npoints has to be grater than zero!')
        end if
    else
        write(DEV_OUT,45) NModePoints
    end if


 30  format ('Normal mode amplitude (amplitude)      = ',f16.3)
 35  format ('Normal mode amplitude (amplitude)      = ',f16.3,'              (default)')
 40  format ('Number of points (npoints)             = ',i12)
 45  format ('Number of points (npoints)             = ',i12,'                  (default)')

end subroutine read_nmodes_method

!===============================================================================
! subroutine read_systematic_opt_method
!===============================================================================

subroutine read_systematic_opt_method(fin)

    use prmfile
    use ffdev_genpoints_dat
    use ffdev_utils

    implicit none
    type(PRMFILE_TYPE)  :: fin
    ! --------------------------------------------------------------------------

    write(DEV_OUT,'(/,a)') '=== [systematic-opt] ==========================================================='

    ! open first section
    if( .not. prmfile_open_section(fin,'systematic-opt') ) then
        write(DEV_OUT,35) TickAngle*DEV_R2D
        return
    end if

    if( prmfile_get_real8_by_key(fin,'tick', TickAngle)) then
        write(DEV_OUT,30) TickAngle
        if( MaxEnergy .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'tick has to be grater than zero!')
        end if
        TickAngle = TickAngle*DEV_D2R
    else
        write(DEV_OUT,35) TickAngle*DEV_R2D
    end if

    if( prmfile_get_real8_by_key(fin,'force_constant', ForceConstant)) then
        write(DEV_OUT,40) ForceConstant
        if( MaxEnergy .le. 0.0d0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'force_constant has to be grater than zero!')
        end if
    else
        write(DEV_OUT,45) ForceConstant
    end if

 30  format ('Dihedral angle tick (tick)             = ',f16.3)
 35  format ('Dihedral angle tick (tick)             = ',f16.3,'              (default)')
 40  format ('Force constant (force_constant)        = ',f16.3)
 45  format ('Force constant (force_constant)        = ',f16.3,'              (default)')

end subroutine read_systematic_opt_method

!===============================================================================

end module ffdev_genpoints_control

! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_errors

use ffdev_errors_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_errors_init_all
! ==============================================================================

subroutine ffdev_errors_init_all()

    use ffdev_errors_dat
    use ffdev_utils

    use ffdev_err_bonds
    use ffdev_err_angles
    use ffdev_err_dihedrals
    use ffdev_err_impropers
    use ffdev_err_nbdists
    use ffdev_err_energy
    use ffdev_err_rmsd
    use ffdev_err_ihess
    use ffdev_err_sapt
    use ffdev_err_pacpnl
    use ffdev_err_zerograd
    use ffdev_err_probe
    use ffdev_err_pbpnl
    use ffdev_err_nbpnl
    use ffdev_err_nbr0
    use ffdev_err_nbc6
    use ffdev_err_qnb
    use ffdev_err_mue
    use ffdev_err_aimr0
    use ffdev_err_aimxdm

    implicit none
    ! --------------------------------------------------------------------------

    ! clear what should be calculated
    errors_calc_ene     = .false.
    errors_calc_sapt    = .false.
    errors_calc_grad    = .false.
    errors_calc_hess    = .false.

    ! reset all setup
    call ffdev_err_energy_init()
    call ffdev_err_sapt_init()
    call ffdev_err_probe_init()

    call ffdev_err_bonds_init()
    call ffdev_err_angles_init()
    call ffdev_err_dihedrals_init()
    call ffdev_err_impropers_init()
    call ffdev_err_nbdists_init()
    call ffdev_err_rmsd_init()
    call ffdev_err_pacpnl_init()
    call ffdev_err_zerograd_init()

    call ffdev_err_pbpnl_init()
    call ffdev_err_qnb_init()
    call ffdev_err_nbpnl_init()
    call ffdev_err_nbr0_init()
    call ffdev_err_nbc6_init()
    call ffdev_err_mue_init()

    call ffdev_err_aimr0_init()

    call ffdev_err_aimxdm_init()

end subroutine ffdev_errors_init_all

! ==============================================================================
! subroutine ffdev_errors_error_setup_domains
! ==============================================================================

subroutine ffdev_errors_error_setup_domains(opterror)

    use ffdev_errors_dat
    use ffdev_err_energy_dat
    use ffdev_err_sapt_dat
    use ffdev_err_zerograd_dat
    use ffdev_err_probe_dat

    implicit none
    logical         :: opterror
    ! --------------------------------------------------------------------------

    ! clear what should be calculated
    errors_calc_ene     = .false.
    errors_calc_sapt    = .false.
    errors_calc_grad    = .false.
    errors_calc_hess    = .false.

    errors_calc_ene  = EnableEnergyError .or. EnableProbeError
    errors_calc_sapt = EnableSAPTError
    errors_calc_grad = EnableZeroGradError

    if( .not. opterror ) then
        errors_calc_ene  = errors_calc_ene  .or. PrintEnergyErrorSummary .or. PrintProbeErrorSummary
        errors_calc_sapt = errors_calc_sapt .or. PrintSAPTErrorSummary
        errors_calc_grad = errors_calc_grad .or. PrintZeroGradErrorSummary
    end if

end subroutine ffdev_errors_error_setup_domains

! ==============================================================================
! subroutine ffdev_errors_error_only
! ==============================================================================

subroutine ffdev_errors_error_only(error)

    use ffdev_err_bonds_dat
    use ffdev_err_bonds

    use ffdev_err_angles_dat
    use ffdev_err_angles

    use ffdev_err_dihedrals_dat
    use ffdev_err_dihedrals

    use ffdev_err_nbdists_dat
    use ffdev_err_nbdists

    use ffdev_err_energy_dat
    use ffdev_err_energy

    use ffdev_err_ihess_dat
    use ffdev_err_ihess

    use ffdev_err_impropers_dat
    use ffdev_err_impropers

    use ffdev_err_rmsd_dat
    use ffdev_err_rmsd

    use ffdev_err_sapt_dat
    use ffdev_err_sapt

    use ffdev_err_pacpnl_dat
    use ffdev_err_pacpnl

    use ffdev_err_zerograd_dat
    use ffdev_err_zerograd

    use ffdev_err_probe_dat
    use ffdev_err_probe

    use ffdev_err_pbpnl_dat
    use ffdev_err_pbpnl

    use ffdev_err_qnb_dat
    use ffdev_err_qnb

    use ffdev_err_nbpnl_dat
    use ffdev_err_nbpnl

    use ffdev_err_nbr0_dat
    use ffdev_err_nbr0

    use ffdev_err_nbc6_dat
    use ffdev_err_nbc6

    use ffdev_err_mue_dat
    use ffdev_err_mue

    use ffdev_err_aimr0_dat
    use ffdev_err_aimr0

    use ffdev_err_aimxdm_dat
    use ffdev_err_aimxdm

    use ffdev_timers

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------------------------------------

    call ffdev_timers_start_timer(FFDEV_ERRORS_TIMER)

    error%total = 0.0d0
    error%energy = 0.0d0
    error%bonds = 0.0d0
    error%angles = 0.0d0
    error%dihedrals = 0.0d0
    error%impropers = 0.0d0
    error%nbdists = 0.0d0
    error%rmsd = 0.0d0
    error%ihess_bonds = 0.0d0
    error%ihess_angles = 0.0d0
    error%ihess_dihedrals = 0.0d0
    error%ihess_impropers = 0.0d0
    error%sapt_ele = 0.0d0
    error%sapt_ind = 0.0d0
    error%sapt_rep = 0.0d0
    error%sapt_dis = 0.0d0
    error%probe_ene = 0.0d0
    error%pacpnl = 0.0d0
    error%zerograd = 0.0d0
    error%pbpnl = 0.0d0
    error%nbpnl = 0.0d0
    error%nbr0 = 0.0d0
    error%nbc6 = 0.0d0
    error%qnb = 0.0d0
    error%mue = 0.0d0
    error%aimr0 = 0.0d0
    error%aimxdm = 0.0d0

! energy based errors
    if( EnableEnergyError ) then
        call ffdev_err_energy_error(error)
        error%total = error%total + error%energy*EnergyErrorWeight
    end if

    if( EnableSAPTError ) then
        call ffdev_err_sapt_error(error)
        error%total = error%total + error%sapt_ele * SAPTEleErrorWeight &
                                  + error%sapt_ind * SAPTIndErrorWeight &
                                  + error%sapt_rep * SAPTRepErrorWeight &
                                  + error%sapt_dis * SAPTDisErrorWeight
    end if

    if( EnableProbeError ) then
        call ffdev_err_probe_error(error)
        error%total = error%total + error%probe_ene * ProbeErrorWeight
    end if

! geometry based errors
    if( EnableBondsError ) then
        call ffdev_err_bonds_error(error)
        error%total = error%total + error%bonds*BondErrorsWeight
    end if

    if( EnableAnglesError ) then
        call ffdev_err_angles_error(error)
        error%total = error%total + error%angles*AngleErrorsWeight
    end if

    if( EnableDihedralsError ) then
        call ffdev_err_dihedrals_error(error)
        error%total = error%total + error%dihedrals*DihedralsErrorWeight
    end if

    if( EnableImpropersError ) then
        call ffdev_err_impropers_error(error)
        error%total = error%total + error%impropers*ImpropersErrorWeight
    end if

    if( EnableNBDistsError ) then
        call ffdev_err_nbdists_error(error)
        error%total = error%total + error%nbdists*NBDistsErrorWeight
    end if

    if( EnableRMSDError ) then
        call ffdev_err_rmsd_error(error)
        error%total = error%total + error%rmsd*RMSDErrorWeight
    end if

    if( EnablePACPnlError ) then
        call ffdev_err_pacpnl_error(error)
        error%total = error%total + error%pacpnl*PACPnlErrorWeight
    end if

    if( EnableZeroGradError ) then
        call ffdev_err_zerograd_error(error)
        error%total = error%total + error%zerograd*ZeroGradErrorWeight
    end if

    if( EnablePBPnlError ) then
        call ffdev_err_pbpnl_error(error)
        error%total = error%total + error%pbpnl*PBPnlErrorWeight
    end if

    if( EnableQNBError ) then
        call ffdev_err_qnb_error(error)
        error%total = error%total + error%qnb*QNBErrorWeight
    end if

    if( EnableNBPnlError ) then
        call ffdev_err_nbpnl_error(error)
        error%total = error%total + error%nbpnl*NBPnlErrorWeight
    end if

    if( EnableNBR0Error ) then
        call ffdev_err_nbr0_error(error)
        error%total = error%total + error%nbr0*NBR0ErrorWeight
    end if

    if( EnableNBC6Error ) then
        call ffdev_err_nbc6_error(error)
        error%total = error%total + error%nbc6*NBC6ErrorWeight
    end if

    if( EnableMUEError ) then
        call ffdev_err_mue_error(error)
        error%total = error%total + error%mue*MUEErrorWeight
    end if

    if( EnableAIMR0Error ) then
        call ffdev_err_aimr0_error(error)
        error%total = error%total + error%aimr0*AIMR0ErrorWeight
    end if

    if( EnableAIMXDMError ) then
        call ffdev_err_aimxdm_error(error)
        error%total = error%total + error%aimxdm*AIMXDMErrorWeight
    end if

    call ffdev_timers_stop_timer(FFDEV_ERRORS_TIMER)

end subroutine ffdev_errors_error_only

!===============================================================================
! subroutine ffdev_errors_ffopt_header_I
!===============================================================================

subroutine ffdev_errors_ffopt_header_I()

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_ihess_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat
    use ffdev_err_sapt_dat
    use ffdev_err_pacpnl_dat
    use ffdev_err_zerograd_dat
    use ffdev_err_probe_dat
    use ffdev_err_pbpnl_dat
    use ffdev_err_qnb_dat
    use ffdev_err_nbpnl_dat
    use ffdev_err_nbr0_dat
    use ffdev_err_nbc6_dat
    use ffdev_err_mue_dat
    use ffdev_err_aimr0_dat
    use ffdev_err_aimxdm_dat

    implicit none
    ! --------------------------------------------------------------------------

    if( EnableEnergyError ) then
        write(DEV_OUT,30,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,21,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,22,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,41,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,42,ADVANCE='NO')
    end if
    if( EnableProbeError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableBondsError ) then
        write(DEV_OUT,33,ADVANCE='NO')
    end if
    if( EnableAnglesError ) then
        write(DEV_OUT,34,ADVANCE='NO')
    end if
    if( EnableDihedralsError ) then
        write(DEV_OUT,35,ADVANCE='NO')
    end if
    if( EnableImpropersError ) then
        write(DEV_OUT,38,ADVANCE='NO')
    end if
    if( EnableNBDistsError ) then
        write(DEV_OUT,36,ADVANCE='NO')
    end if
    if( EnableRMSDError ) then
        write(DEV_OUT,39,ADVANCE='NO')
    end if
    if( EnablePACPnlError ) then
        write(DEV_OUT,43,ADVANCE='NO')
    end if
    if( EnableZeroGradError ) then
        write(DEV_OUT,44,ADVANCE='NO')
    end if
    if( EnablePBPnlError ) then
        write(DEV_OUT,60,ADVANCE='NO')
    end if
    if( EnableQNBError ) then
        write(DEV_OUT,70,ADVANCE='NO')
    end if
    if( EnableNBPnlError ) then
        write(DEV_OUT,90,ADVANCE='NO')
    end if
    if( EnableNBR0Error ) then
        write(DEV_OUT,95,ADVANCE='NO')
    end if
    if( EnableNBC6Error ) then
        write(DEV_OUT,96,ADVANCE='NO')
    end if
    if( EnableMUEError ) then
        write(DEV_OUT,80,ADVANCE='NO')
    end if
    if( EnableAIMR0Error ) then
        write(DEV_OUT,110,ADVANCE='NO')
    end if
    if( EnableAIMXDMError ) then
        write(DEV_OUT,120,ADVANCE='NO')
    end if

 30 format('       Energy')
 33 format('        Bonds')
 34 format('       Angles')
 35 format('    Dihedrals')
 36 format('       d(NBs)')
 38 format('    Impropers')
 39 format('         RMSD')
 21 format('    SAPT(Ele)')
 22 format('    SAPT(Ind)')
 41 format('    SAPT(Rep)')
 42 format('    SAPT(Dis)')
 43 format('   PACPenalty')
 44 format(' ZeroGradient')
 50 format('     ProbeEne')
 60 format('    PBPenalty')
 70 format('          QNB')
 90 format('    NBPenalty')
 95 format('    R0Penalty')
 96 format('    C6Penalty')
 80 format('          MUE')
110 format('        AIMR0')
120 format('        AIMXDM')

end subroutine ffdev_errors_ffopt_header_I

!===============================================================================
! subroutine ffdev_errors_ffopt_header_II
!===============================================================================

subroutine ffdev_errors_ffopt_header_II()

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat
    use ffdev_err_sapt_dat
    use ffdev_err_pacpnl_dat
    use ffdev_err_zerograd_dat
    use ffdev_err_probe_dat
    use ffdev_err_pbpnl_dat
    use ffdev_err_qnb_dat
    use ffdev_err_nbpnl_dat
    use ffdev_err_nbr0_dat
    use ffdev_err_nbc6_dat
    use ffdev_err_mue_dat
    use ffdev_err_aimr0_dat
    use ffdev_err_aimxdm_dat


    implicit none
    ! --------------------------------------------------------------------------

    if( EnableEnergyError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableProbeError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableBondsError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableAnglesError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableDihedralsError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableImpropersError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableNBDistsError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableRMSDError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnablePACPnlError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableZeroGradError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnablePBPnlError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableQNBError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableNBPnlError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableNBR0Error ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableNBC6Error ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableMUEError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableAIMR0Error ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if
    if( EnableAIMXDMError ) then
        write(DEV_OUT,50,ADVANCE='NO')
    end if


 50 format(' ------------')

end subroutine ffdev_errors_ffopt_header_II

!===============================================================================
! subroutine ffdev_errors_ffopt_header_II
!===============================================================================

subroutine ffdev_errors_ffopt_results(error)

    use ffdev_err_bonds_dat
    use ffdev_err_angles_dat
    use ffdev_err_dihedrals_dat
    use ffdev_err_nbdists_dat
    use ffdev_err_energy_dat
    use ffdev_err_impropers_dat
    use ffdev_err_rmsd_dat
    use ffdev_err_sapt_dat
    use ffdev_err_pacpnl_dat
    use ffdev_err_zerograd_dat
    use ffdev_err_probe_dat
    use ffdev_err_pbpnl_dat
    use ffdev_err_qnb_dat
    use ffdev_err_nbpnl_dat
    use ffdev_err_nbr0_dat
    use ffdev_err_nbc6_dat
    use ffdev_err_mue_dat
    use ffdev_err_aimr0_dat
    use ffdev_err_aimxdm_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! -----------------------------------------------------------------------------

    if( EnableEnergyError ) then
        write(DEV_OUT,15,ADVANCE='NO') EnergyErrorWeight*error%energy
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,15,ADVANCE='NO') SAPTEleErrorWeight*error%sapt_ele
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,15,ADVANCE='NO') SAPTIndErrorWeight*error%sapt_ind
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,15,ADVANCE='NO') SAPTRepErrorWeight*error%sapt_rep
    end if
    if( EnableSAPTError ) then
        write(DEV_OUT,15,ADVANCE='NO') SAPTDisErrorWeight*error%sapt_dis
    end if
    if( EnableProbeError ) then
        write(DEV_OUT,15,ADVANCE='NO') ProbeErrorWeight*error%probe_ene
    end if
    if( EnableBondsError ) then
        write(DEV_OUT,15,ADVANCE='NO') BondErrorsWeight*error%bonds
    end if
    if( EnableAnglesError ) then
        write(DEV_OUT,15,ADVANCE='NO') AngleErrorsWeight*error%angles
    end if
    if( EnableDihedralsError ) then
        write(DEV_OUT,15,ADVANCE='NO') DihedralsErrorWeight*error%dihedrals
    end if
    if( EnableImpropersError ) then
        write(DEV_OUT,15,ADVANCE='NO') ImpropersErrorWeight*error%impropers
    end if
    if( EnableNBDistsError ) then
        write(DEV_OUT,15,ADVANCE='NO') NBDistsErrorWeight*error%nbdists
    end if
    if( EnableRMSDError ) then
        write(DEV_OUT,15,ADVANCE='NO') RMSDErrorWeight*error%rmsd
    end if
    if( EnablePACPnlError ) then
        write(DEV_OUT,15,ADVANCE='NO') PACPnlErrorWeight*error%pacpnl
    end if
    if( EnableZeroGradError ) then
        write(DEV_OUT,15,ADVANCE='NO') ZeroGradErrorWeight*error%zerograd
    end if
    if( EnablePBPnlError ) then
        write(DEV_OUT,15,ADVANCE='NO') PBPnlErrorWeight*error%pbpnl
    end if
    if( EnableQNBError ) then
        write(DEV_OUT,15,ADVANCE='NO') QNBErrorWeight*error%qnb
    end if
    if( EnableNBPnlError ) then
        write(DEV_OUT,15,ADVANCE='NO') NBPnlErrorWeight*error%nbpnl
    end if
    if( EnableNBR0Error ) then
        write(DEV_OUT,15,ADVANCE='NO') NBR0ErrorWeight*error%nbr0
    end if
    if( EnableNBC6Error ) then
        write(DEV_OUT,15,ADVANCE='NO') NBC6ErrorWeight*error%nbc6
    end if
    if( EnableMUEError ) then
        write(DEV_OUT,15,ADVANCE='NO') MUEErrorWeight*error%mue
    end if
    if( EnableAIMR0Error ) then
        write(DEV_OUT,15,ADVANCE='NO') AIMR0ErrorWeight*error%aimr0
    end if
    if( EnableAIMXDMError ) then
        write(DEV_OUT,15,ADVANCE='NO') AIMXDMErrorWeight*error%aimxdm
    end if

 15 format(1X,E12.5)

end subroutine ffdev_errors_ffopt_results

! ==============================================================================
! subroutine ffdev_errors_summary
! ==============================================================================

subroutine ffdev_errors_summary(logmode)

    use ffdev_targetset_dat
    use ffdev_geometry
    use ffdev_utils

    use ffdev_err_bonds_dat
    use ffdev_err_bonds

    use ffdev_err_angles_dat
    use ffdev_err_angles

    use ffdev_err_dihedrals_dat
    use ffdev_err_dihedrals

    use ffdev_err_nbdists_dat
    use ffdev_err_nbdists

    use ffdev_err_energy_dat
    use ffdev_err_energy

    use ffdev_err_ihess_dat
    use ffdev_err_ihess

    use ffdev_err_impropers_dat
    use ffdev_err_impropers

    use ffdev_err_rmsd_dat
    use ffdev_err_rmsd

    use ffdev_err_sapt_dat
    use ffdev_err_sapt

    use ffdev_err_pacpnl_dat
    use ffdev_err_pacpnl

    use ffdev_err_zerograd_dat
    use ffdev_err_zerograd

    use ffdev_err_probe_dat
    use ffdev_err_probe

    use ffdev_err_pbpnl_dat
    use ffdev_err_pbpnl

    use ffdev_err_nbr0_dat
    use ffdev_err_nbr0

    use ffdev_err_nbc6_dat
    use ffdev_err_nbc6

    use ffdev_err_qnb_dat
    use ffdev_err_qnb

    use ffdev_err_nbpnl_dat
    use ffdev_err_nbpnl

    use ffdev_err_aimr0_dat
    use ffdev_err_aimr0

    use ffdev_err_aimxdm_dat
    use ffdev_err_aimxdm

    implicit none
    integer     :: logmode
    logical     :: printme, printsum
    integer     :: i,j
    ! --------------------------------------------------------------------------

    if( .not. (PrintEnergyErrorSummary .or. PrintSAPTErrorSummary .or. PrintZeroGradErrorSummary .or.  &
            PrintBondsErrorSummary .or. PrintAnglesErrorSummary .or. PrintDihedralsErrorSummary .or. &
            PrintImpropersErrorSummary .or. PrintProbeErrorSummary .or. &
            PrintNBDistsErrorSummary .or. PrintRMSDErrorSummary .or. PrintPACPnlErrorSummary .or. &
            PrintPBPnlErrorSummary .or. PrintQNBErrorSummary .or. PrintNBPnlErrorSummary .or. PrintNBR0ErrorSummary .or. &
            PrintNBC6ErrorSummary .or. PrintAIMR0ErrorSummary .or. PrintAIMXDMErrorSummary ) ) then
        ! no error to report
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,1)
    select case(logmode)
        case(SMMLOG_INITIAL)
            call ffdev_utils_heading(DEV_OUT,'Initial Error Statistics',':')
        case(SMMLOG_INTERMEDIATE)
            call ffdev_utils_heading(DEV_OUT,'Intermediate Error Statistics',':')
        case(SMMLOG_FINAL)
            call ffdev_utils_heading(DEV_OUT,'Final Error Statistics',':')
    end select
    write(DEV_OUT,1)

    ! individual summaries
    if( PrintPBPnlErrorSummary ) then
        call ffdev_err_pbpnl_summary
    end if

    if( PrintQNBErrorSummary ) then
        call ffdev_err_qnb_summary
    end if

    if( PrintNBPnlErrorSummary ) then
        call ffdev_err_nbpnl_summary
    end if

    if( PrintNBR0ErrorSummary ) then
        call ffdev_err_nbr0_summary
    end if

    if( PrintNBC6ErrorSummary ) then
        call ffdev_err_nbc6_summary
    end if

    if( PrintEnergyErrorSummary ) then
        call ffdev_err_energy_summary
    end if

    if( PrintSAPTErrorSummary ) then
        call ffdev_err_sapt_summary
    end if

    if( EnablePACPnlError ) then
        call ffdev_err_pacpnl_summary
    end if

    if( PrintProbeErrorSummary ) then
        call ffdev_err_probe_summary
    end if

    if( PrintZeroGradErrorSummary ) then
        call ffdev_err_zerograd_summary
    end if

    if( PrintAIMR0ErrorSummary ) then
        call ffdev_err_aimr0_summary
    end if

    if( PrintAIMXDMErrorSummary ) then
        call ffdev_err_aimxdm_summary
    end if

    ! summary per sets
    if( PrintRMSDErrorSummary ) then
        write(DEV_OUT,*)
        write(DEV_OUT,10)

        do i=1,nsets
            printme = .false.
            if( PrintRMSDErrorSummary ) then
                printsum = .false.
                call ffdev_err_rmsd_summary(sets(i),printsum)
                printme = printme .or. printsum
            end if

            if( .not. printme ) cycle

            printsum = .true.

            write(DEV_OUT,*)
            write(DEV_OUT,5) i
            if( PrintRMSDErrorSummary ) then
                call ffdev_err_rmsd_summary(sets(i),printsum)
            end if
        end do
    end if

    ! summary per points
    if( PrintBondsErrorSummary .or. PrintAnglesErrorSummary .or. PrintDihedralsErrorSummary .or. &
        PrintImpropersErrorSummary .or. &
        PrintNBDistsErrorSummary .or. PrintRMSDErrorSummary .or. PrintPACPnlErrorSummary) then

        write(DEV_OUT,*)
        write(DEV_OUT,20)

        do i=1,nsets
            do j=1,sets(i)%ngeos
                printme = .false.
                if( PrintBondsErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_bonds_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintAnglesErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_angles_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintDihedralsErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_dihedrals_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintImpropersErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_impropers_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if
                if( PrintNBDistsErrorSummary ) then
                    printsum = .false.
                    call ffdev_err_nbdists_summary(sets(i)%top,sets(i)%geo(j),printsum)
                    printme = printme .or. printsum
                end if

                if( .not. printme ) cycle

                write(DEV_OUT,*)
                write(DEV_OUT,6) i,j
                printsum = .true.
                if( PrintBondsErrorSummary ) then
                    call ffdev_err_bonds_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintAnglesErrorSummary ) then
                    call ffdev_err_angles_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintDihedralsErrorSummary ) then
                    call ffdev_err_dihedrals_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintImpropersErrorSummary ) then
                    call ffdev_err_impropers_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
                if( PrintNBDistsErrorSummary ) then
                    call ffdev_err_nbdists_summary(sets(i)%top,sets(i)%geo(j),printsum)
                end if
            end do
        end do
    end if

 1 format('# ==============================================================================')
 5 format('== [SET] #',I2.2,' ===================================================================')
 6 format('== [SET#',I5.5,']/[GEO#',I6.6,'] ====================================================')
10 format('== # SUMMARY PER SETS #')
20 format('== # SUMMARY PER POINTS #')

end subroutine ffdev_errors_summary

! ------------------------------------------------------------------------------

end module ffdev_errors

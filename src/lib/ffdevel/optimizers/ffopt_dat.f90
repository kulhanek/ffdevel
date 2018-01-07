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

module ffdev_ffopt_dat

use ffdev_sizes
use ffdev_constants
use smf_xyzfile_type

implicit none

! WARNING: default values are set in ffdev_ffopt_set_default() !!!!!!!!!!!!!!!!!
! this is because the programs can be repeated

! ------------------------------------------------------------------------------
! minimization methods
integer, parameter      :: MINIMIZATION_STEEPEST_DESCENT    = 0
integer, parameter      :: MINIMIZATION_LBFGS               = 1
integer, parameter      :: MINIMIZATION_NLOPT               = 2
integer, parameter      :: MINIMIZATION_SA                  = 3

real(DEVDP),allocatable :: FFParams(:)
real(DEVDP),allocatable :: FFParamsGrd(:)

! === [minimization] ===========================================================
integer         :: OptimizationMethod
integer         :: NOptSteps            ! max number of steps
integer         :: OutSamples           ! how often write results

! maximum number of steps is nsteps - this is becuase of change of restraints etc 
real(DEVDP)     :: MaxRMSG
real(DEVDP)     :: MaxG
real(DEVDP)     :: MinErrorChange      ! negative number - this test is switched off by default
logical         :: PrintFinalGradient

! === [steepest-descent] =======================================================
real(DEVDP)     :: InitialStepSize
real(DEVDP)     :: MaximalStepSize
real(DEVDP)     :: AcceptRatio
real(DEVDP)     :: RejectRatio
logical         :: AdaptiveStep

! === [L-BFGS] =================================================================
integer         :: NumberOfCorrections

! === [NLOPT] ==================================================================
integer(8)      :: NLoptID
integer(8)      :: NLoptDummy               ! FIXME - probleme with memory aligment?
real(DEVDP)     :: NLOpt_InitialStep        ! nlo_set_initial_step

! === [SA] =====================================================================
real(DEVDP)     :: OptSA_Temp   = 0.01d0    ! initial temperature
real(DEVDP)     :: OptSA_RT     = 0.85d0    ! the temperature reduction factor
real(DEVDP)     :: OptSA_EPS    = 0.01      ! Error tolerance for termination
integer         :: OptSA_NS     = 20        ! Number of cycles
integer         :: OptSA_NT     = 10        ! Number of iterations before temperature reduction
integer         :: OptSA_NEPS   = 4         ! Number of final function values used to decide upon termination
integer         :: OptSA_MAXEVL = 80000     ! The maximum number of function evaluations.
integer         :: OptSA_IPRINT = 1         ! controls printing inside SA
integer         :: OptSA_ISEED  = 341723    ! random generator seed

real(DEVDP),allocatable         :: sa_lb(:)
real(DEVDP),allocatable         :: sa_ub(:)
real(DEVDP),allocatable         :: sa_c(:)
real(DEVDP),allocatable         :: sa_vm(:)
real(DEVDP),allocatable         :: sa_xopt(:)
real(DEVDP),allocatable         :: sa_fstar(:)
real(DEVDP),allocatable         :: sa_xp(:)
integer,allocatable             :: sa_nacp(:)

end module ffdev_ffopt_dat

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
integer, parameter      :: MINIMIZATION_SINGLE_POINT        = -1 ! single-point (only internal)
integer, parameter      :: MINIMIZATION_STEEPEST_DESCENT    =  0
integer, parameter      :: MINIMIZATION_LBFGS               =  1
integer, parameter      :: MINIMIZATION_NLOPT               =  2
integer, parameter      :: MINIMIZATION_SHARK               =  3

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
integer         :: NLOpt_Method
real(DEVDP)     :: NLOpt_InitialStep

! === [Shark] ==================================================================
integer         :: Shark_Method
real(DEVDP)     :: Shark_InitialStep
integer         :: Shark_RngSeed
logical         :: Shark_EnableBoxing

! ------------------------------------------------------------------------------
integer, parameter      :: SHARK_CMA_ES = 1     ! CMA-ES

! ------------------------------------------------------------------------------
! working directories
integer                 :: FFOptFceEvals
real(DEVDP),allocatable :: tmp_xg(:),tmp_ub(:),tmp_lb(:)

end module ffdev_ffopt_dat

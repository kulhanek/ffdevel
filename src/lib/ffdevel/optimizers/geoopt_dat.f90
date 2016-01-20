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

module ffdev_geoopt_dat

use ffdev_sizes
use ffdev_constants
use smf_xyzfile_type

implicit none

! ------------------------------------------------------------------------------
! minimization methods
integer, parameter     :: MINIMIZATION_STEEPEST_DESCENT    = 0
integer, parameter     :: MINIMIZATION_LBFGS               = 1

! === [files] ==================================================================
character(len=MAX_PATH) :: OptTopName      ! input topology name
character(len=MAX_PATH) :: OptCrdName      ! input coordinate name
character(len=MAX_PATH) :: OptRstName  = 'final.xyz'    ! final coordinates name
character(len=MAX_PATH) :: OptTrajName = 'traj.xyz'     ! trajectory name

! === [minimization] ===========================================================
integer         :: OptimizationMethod  = MINIMIZATION_LBFGS
integer         :: NOptSteps    = 1000      ! max number of steps
integer         :: OutSamples   =   20      ! how often write results
integer         :: TrajSamples  =  100      ! how often write trajectory

! maximum number of steps is nsteps - this is becuase of change of restraints etc 
real(DEVDP)     :: MaxRMSG             = 0.01d0
real(DEVDP)     :: MaxG                = 0.01d0
real(DEVDP)     :: MinEnergyChange     = -1        ! negative number - this test is switched off by default
logical         :: PrintFinalGradient  = .false.

! === [steepest-descent] =======================================================
real(DEVDP)     :: InitialStepSize     = 0.001
real(DEVDP)     :: MaximalStepSize     = 0.010
real(DEVDP)     :: AcceptRatio         = 1.2000
real(DEVDP)     :: RejectRatio         = 0.5000
logical         :: AdaptiveStep        = .true.

! === [L-BFGS] =================================================================
integer         :: NumberOfCorrections = 50

!-------------------------------------------------------------------------------
type(XYZFILE_TYPE)  :: OptTrajFile

end module ffdev_geoopt_dat

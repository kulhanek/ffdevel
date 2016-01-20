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

module ffdev_genpoints_dat

use ffdev_sizes
use ffdev_constants

implicit none

! ------------------------------------------------------------------------------
! generator methods
integer, parameter     :: GENPOINTS_SYSTEMATIC_METHOD           = 0
integer, parameter     :: GENPOINTS_STOCHASTIC_METHOD           = 1
integer, parameter     :: GENPOINTS_STEPBYSTEP_METHOD           = 2
integer, parameter     :: GENPOINTS_STOCHASTIC_BY_STEPS_METHOD  = 3
integer, parameter     :: GENPOINTS_NMODES_METHOD               = 4
integer, parameter     :: GENPOINTS_SYSTEMATIC_OPT_METHOD       = 5
integer, parameter     :: GENPOINTS_STOCHASTIC_OPT_METHOD       = 6

! === [files] ==================================================================
character(len=MAX_PATH) :: GenTopName      ! input topology name
character(len=MAX_PATH) :: GenCrdName      ! input coordinate name
character(len=MAX_PATH) :: GenRotName      ! list of rotors
character(len=MAX_PATH) :: GenOutName = 'points.xyz'      ! output name
character(len=MAX_PATH) :: GenFinName = 'final.xyz'       ! structure with minimum energy
logical                 :: RotorListLoaded = .false.

! === [points] =================================================================
integer         :: GeneratorMethod  = GENPOINTS_SYSTEMATIC_METHOD
integer         :: MaxPoints        =  100      ! max number of points
real(DEVDP)     :: MaxEnergy        =   25.0    ! max energy to global minima
logical         :: OptimizePoints   = .false.

! === [systematic] =============================================================
! === [systematic-opt] =========================================================
real(DEVDP)     :: TickAngle        = 20.0 * DEV_PI / 180.0     ! in rad
real(DEVDP)     :: ForceConstant    = 0.05                      ! in kcal/mol/rad^2

! === [nmodes] =================================================================
real(DEVDP)     :: NModeAmplitude   = 0.4
integer         :: NModePoints      = 3

! ------------------------------------------------------------------------------

type ROTOR
    integer     :: ai
    integer     :: aj
    real(DEVDP) :: angle
end type ROTOR

integer                     :: NRotors      ! number of rotors
type(ROTOR),allocatable     :: Rotors(:)    ! list of rotor generators

! ------------------------------------------------------------------------------

type GENPOINT
    real(DEVDP),pointer :: crd(:,:)         ! geometry (3,natoms)
    real(DEVDP)         :: energy
    real(DEVDP),pointer :: angles(:)
    integer             :: active_angle
end type GENPOINT

integer                     :: NPoints      ! number of points
type(GENPOINT),allocatable  :: Points(:)    ! generated points

real(DEVDP)                 :: MinEnergy
integer                     :: MinEnergyPtsIndex   = 0 ! index to points
integer                     :: CurrPtsIndex        = 0 ! current point index
integer                     :: RejectedPoints
logical,allocatable         :: AtomMask(:)          ! what atom should be rotated
integer,allocatable         :: ProcessingStack(:)

real(DEVDP)                 :: MinOptEnergy
integer                     :: MinOptEnergyPtsIndex    ! index to points

! ------------------------------------------------------------------------------

end module ffdev_genpoints_dat

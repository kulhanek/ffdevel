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

module ffdev_parameters_dat

use ffdev_geometry_dat
use ffdev_constants

! ------------------------------------------------------------------------------

integer,parameter       :: REALM_EOFFSET    = 1
integer,parameter       :: REALM_BOND_D0    = 2
integer,parameter       :: REALM_BOND_K     = 3
integer,parameter       :: REALM_ANGLE_A0   = 4
integer,parameter       :: REALM_ANGLE_K    = 5
integer,parameter       :: REALM_DIH_V      = 6
integer,parameter       :: REALM_DIH_G      = 7
integer,parameter       :: REALM_DIH_SCEE   = 8
integer,parameter       :: REALM_DIH_SCNB   = 9
integer,parameter       :: REALM_IMPR_V     = 10
integer,parameter       :: REALM_IMPR_G     = 11
integer,parameter       :: REALM_DIH_C      = 12
integer,parameter       :: REALM_VDW_A      = 13
integer,parameter       :: REALM_VDW_B      = 14
integer,parameter       :: REALM_VDW_C      = 15

integer,parameter       :: REALM_FIRST   = REALM_EOFFSET
integer,parameter       :: REALM_LAST    = REALM_VDW_C

! ------------------------------------------------------------------------------

! parameter description
type PARAMETER_TYPE
    real(DEVDP)         :: value        ! parameter value
    integer             :: realm        ! parameter realm (bond_k,bond_r0,...)
    integer             :: pn           ! dihedral pn
    integer,pointer     :: ids(:)       ! type ids in topologies
    logical             :: enabled      ! parameter is anbled
    integer             :: identity     ! identity index
    integer             :: ti,tj,tk,tl  ! common atom types
    integer             :: pidx         ! helper index to prms() array
end type PARAMETER_TYPE

integer                             :: nparams      ! number of parameters
type(PARAMETER_TYPE),allocatable    :: params(:)    ! parameters

integer                             :: nactparms    ! number of active parameters

! ------------------------------------------------------------------------------

! parameter type
type PARM_TYPE
    character(len=4)    :: name
    integer             :: z
    real(DEVDP)         :: mass
    integer,pointer     :: ids(:)
end type PARM_TYPE

integer                     :: ntypes    ! number of types
type(PARM_TYPE),allocatable :: types(:)  ! types

! ------------------------------------------------------------------------------
! error function

type FFERROR_TYPE
    real(DEVDP)         :: total        ! total error
    real(DEVDP)         :: energy       ! energy part of error
    real(DEVDP)         :: grad         ! gradient part of error
    real(DEVDP)         :: hess         ! hessian part of error
    real(DEVDP)         :: penalty      ! geometry penalty
end type FFERROR_TYPE


type(FFERROR_TYPE)      :: FFError
real(DEVDP),allocatable :: FFParams(:)
real(DEVDP),allocatable :: FFParamsGrd(:)

! === [error] ==================================================================
logical                 :: EnableEnergyError = .true.
real(DEVDP)             :: EnergyErrorWeight = 1.0

logical                 :: EnableGradientError = .true.
real(DEVDP)             :: GradientErrorWeight = 1.0

logical                 :: EnableHessianError = .true.
real(DEVDP)             :: HessianErrorWeight = 1.0

logical                 :: EnablePenaltyError   = .false.
real(DEVDP)             :: PenaltyErrorWeight   = 1.0
real(DEVDP)             :: BondD0PenaltyForceK  = 10.0
real(DEVDP)             :: AngleA0PenaltyForceK = 10.0

! === [files] ==================================================================
character(len=MAX_PATH) :: InpParamFileName = '-none-'              ! input parameters
character(len=MAX_PATH) :: OutParamFileName = 'final.prms'          ! output parameters
character(len=MAX_PATH) :: OutAmberPrmsFileName = 'final.frcmod'    ! output Amber force field

! ------------------------------------------------------------------------------

end module ffdev_parameters_dat

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
use ffdev_variables
use ffdev_topology_dat

! ------------------------------------------------------------------------------

integer,parameter       :: REALM_BOND_D0    = 1
integer,parameter       :: REALM_BOND_K     = 2
integer,parameter       :: REALM_ANGLE_A0   = 3
integer,parameter       :: REALM_ANGLE_K    = 4
integer,parameter       :: REALM_DIH_V      = 5
integer,parameter       :: REALM_DIH_G      = 6
integer,parameter       :: REALM_DIH_SCEE   = 7
integer,parameter       :: REALM_DIH_SCNB   = 8
integer,parameter       :: REALM_IMPR_V     = 9
integer,parameter       :: REALM_IMPR_G     = 10
integer,parameter       :: REALM_DIH_C      = 11

! non-bonded - vdW setup - LJ parameters
integer,parameter       :: REALM_VDW_EPS    = 12
integer,parameter       :: REALM_VDW_R0     = 13

! exp-6 probing
integer,parameter       :: REALM_VDW_ALPHA  = 14
integer,parameter       :: REALM_VDW_ALPHA0 = 15

! non-bonded - vdW setup - repulsion
integer,parameter       :: REALM_VDW_PA     = 16
integer,parameter       :: REALM_VDW_PB     = 17
integer,parameter       :: REALM_VDW_RC     = 18

integer,parameter       :: REALM_VDW_B0     = 19

! non-bonded - vdW setup - dispersion
integer,parameter       :: REALM_DISP_S6    = 20
integer,parameter       :: REALM_DISP_S8    = 21
integer,parameter       :: REALM_DISP_S10   = 22

integer,parameter       :: REALM_DAMP_FA    = 23
integer,parameter       :: REALM_DAMP_FB    = 24

integer,parameter       :: REALM_DAMP_PB    = 25
integer,parameter       :: REALM_DAMP_TB    = 26
integer,parameter       :: REALM_DAMP_PE    = 27

! non-bonded - electrostatics
integer,parameter       :: REALM_PAC        = 28
integer,parameter       :: REALM_ELE_SQ     = 29

integer,parameter       :: REALM_ZEFF       = 30

! non-bonded - electrostatics
integer,parameter       :: REALM_GLB_SCEE   = 31
integer,parameter       :: REALM_GLB_SCNB   = 32

! Pauli repulsion K factors
integer,parameter       :: REALM_K_EXC      = 33
integer,parameter       :: REALM_K_IND      = 34

! AIM
integer,parameter       :: REALM_VDW_R0FREE = 35
integer,parameter       :: REALM_VDW_PBFREE = 36

integer,parameter       :: REALM_FIRST   = REALM_BOND_D0
integer,parameter       :: REALM_LAST    = REALM_VDW_PBFREE

! ------------------------------------------------------------------------------

! parameter description
type PARAMETER_TYPE
    real(DEVDP)         :: value        ! parameter value
    integer             :: realm        ! parameter realm (bond_k,bond_r0,...)
    integer             :: pn           ! dihedral pn
    integer,pointer     :: ids(:)       ! type ids in topologies
    logical             :: enabled      ! parameter is enabled
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
    character(MAX_TNAME)    :: name
    integer                 :: z
    real(DEVDP)             :: mass
    integer,pointer         :: ids(:)
    logical                 :: probe
    ! for amber FF output, LJ parameters
    logical                 :: print_nb
    real(DEVDP)             :: eps
    real(DEVDP)             :: r0
    ! charge statistics
    integer                 :: qcount
    real(DEVDP)             :: minq
    real(DEVDP)             :: maxq
    real(DEVDP)             :: aveq
    real(DEVDP)             :: sdq
    ! NB - parameters
    real(DEVDP)             :: Zeff
    real(DEVDP)             :: PA, PB, RC
end type PARM_TYPE

integer                     :: ntypes    = 0    ! number of types
type(PARM_TYPE),allocatable :: types(:)         ! types

! ------------------------------------------------------------------------------

! charges
integer                     :: nprmchrgs = 0    ! number of charges
real(DEVDP),allocatable     :: prmchrgs         ! values

! ------------------------------------------------------------------------------

! NBParamsMode
integer,parameter       :: NB_PARAMS_MODE_LIKE_ONLY     = 1     ! only like nb_types except of probes
integer,parameter       :: NB_PARAMS_MODE_ALL           = 2     ! all nb types

! ffdev_parameters_print_parameters modes
integer,parameter       :: PARAMS_SUMMARY_INITIAL       = 1
integer,parameter       :: PARAMS_SUMMARY_OPTIMIZED     = 2
integer,parameter       :: PARAMS_SUMMARY_MODIFIED      = 3
integer,parameter       :: PARAMS_SUMMARY_FULL          = 4

! PAC charge source
integer,parameter       :: PAC_SOURCE_TOPOLOGY          = 1
integer,parameter       :: PAC_SOURCE_GEO               = 2
integer,parameter       :: PAC_SOURCE_GEO_HIRSHFELD     = 3

! === [control] ================================================================
integer         :: NBParamsMode         = NB_PARAMS_MODE_LIKE_ONLY  ! mode for determination of NB parameters
logical         :: PACAsPrms            = .false.                   ! partial atomic charges as parameters
logical         :: OnlyDefinedDihItems  = .true.
logical         :: LockDihC_PN1         = .true.
logical         :: ResetAllSetup        = .true.
integer         :: GlbRngSeed           = 5489                      ! random number generator setup
integer         :: PACSource            = PAC_SOURCE_TOPOLOGY

! === [grbf2cos] ===============================================================
integer         :: GRBF2COSMaxN     = 4             ! max length of cos series
real(DEVDP)     :: GRBF2COSMinV     = 0.1d0         ! min amplitude of each cos item

! === [ranges] =================================================================

! values set in ffdev_params_reset_ranges

real(DEVDP)     :: MinBondD0
real(DEVDP)     :: MaxBondD0
real(DEVDP)     :: MinBondK
real(DEVDP)     :: MaxBondK
real(DEVDP)     :: MinAngleA0
real(DEVDP)     :: MaxAngleA0
real(DEVDP)     :: MinAngleK
real(DEVDP)     :: MaxAngleK
real(DEVDP)     :: MinDihV
real(DEVDP)     :: MaxDihV
real(DEVDP)     :: MinDihG
real(DEVDP)     :: MaxDihG
real(DEVDP)     :: MinDihSCEE
real(DEVDP)     :: MaxDihSCEE
real(DEVDP)     :: MinDihSCNB
real(DEVDP)     :: MaxDihSCNB
real(DEVDP)     :: MinImprV
real(DEVDP)     :: MaxImprV
real(DEVDP)     :: MinImprG
real(DEVDP)     :: MaxImprG
real(DEVDP)     :: MinDihC
real(DEVDP)     :: MaxDihC

! non-bonded LJ
real(DEVDP)     :: MinVdwEps
real(DEVDP)     :: MaxVdwEps
real(DEVDP)     :: MinVdwR0
real(DEVDP)     :: MaxVdwR0
real(DEVDP)     :: MinVdwAlpha
real(DEVDP)     :: MaxVdwAlpha
real(DEVDP)     :: MinVdwAlpha0
real(DEVDP)     :: MaxVdwAlpha0

! non-bonded ABC
real(DEVDP)     :: MinVdwPA
real(DEVDP)     :: MaxVdwPA

real(DEVDP)     :: MinVdwPB
real(DEVDP)     :: MaxVdwPB

real(DEVDP)     :: MinVdwB0
real(DEVDP)     :: MaxVdwB0

real(DEVDP)     :: MinVdwRC
real(DEVDP)     :: MaxVdwRC

! non-bonded scaling factors
real(DEVDP)     :: MinEleSQ
real(DEVDP)     :: MaxEleSQ

! partial atomic charges
real(DEVDP)     :: MinPAC
real(DEVDP)     :: MaxPAC

real(DEVDP)     :: MinZeff
real(DEVDP)     :: MaxZeff

! vdW interactions
real(DEVDP)     :: MinDampFA
real(DEVDP)     :: MaxDampFA
real(DEVDP)     :: MinDampFB
real(DEVDP)     :: MaxDampFB

real(DEVDP)     :: MinDampPB
real(DEVDP)     :: MaxDampPB
real(DEVDP)     :: MinDampTB
real(DEVDP)     :: MaxDampTB
real(DEVDP)     :: MinDampPE
real(DEVDP)     :: MaxDampPE

! dispersion scaling
real(DEVDP)     :: MinDispS6
real(DEVDP)     :: MaxDispS6
real(DEVDP)     :: MinDispS8
real(DEVDP)     :: MaxDispS8
real(DEVDP)     :: MinDispS10
real(DEVDP)     :: MaxDispS10

real(DEVDP)     :: MinGlbSCEE
real(DEVDP)     :: MaxGlbSCEE
real(DEVDP)     :: MinGlbSCNB
real(DEVDP)     :: MaxGlbSCNB

! exchange and induction factors
real(DEVDP)     :: MinKExc
real(DEVDP)     :: MaxKExc
real(DEVDP)     :: MinKInd
real(DEVDP)     :: MaxKInd

real(DEVDP)     :: MinVdwR0Free
real(DEVDP)     :: MaxVdwR0Free

real(DEVDP)     :: MinVdwPBFree
real(DEVDP)     :: MaxVdwPBFree

! === [files] ==================================================================
character(MAX_PATH) :: OutParamFileName     = 'final.prms'      ! output parameters
character(MAX_PATH) :: OutAmberPrmsFileName = 'final.frcmod'    ! output Amber force field

! ------------------------------------------------------------------------------

end module ffdev_parameters_dat

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
use ffdev_topology_dat

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
integer,parameter       :: REALM_VDW_EPS    = 13
integer,parameter       :: REALM_VDW_R0     = 14
integer,parameter       :: REALM_VDW_ALPHA  = 15
integer,parameter       :: REALM_VDW_A      = 16
integer,parameter       :: REALM_VDW_B      = 17
integer,parameter       :: REALM_VDW_C      = 18

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
    logical             :: probe
    ! derived LJ parameters
    logical             :: print_nb
    real(DEVDP)         :: eps
    real(DEVDP)         :: r0
end type PARM_TYPE

integer                     :: ntypes    ! number of types
type(PARM_TYPE),allocatable :: types(:)  ! types

! ------------------------------------------------------------------------------
! experimental/unfinished setup
integer                 :: LastNBMode                   = NB_MODE_LJ        ! determine which realms will be activated for NB

! derived setup
logical                 :: ApplyCombinationRules        = .false.           ! apply combination rules in every error evaluation

! NBParamsMode
integer,parameter       :: NB_PARAMS_MODE_NORMAL        = 0     ! normal setup depending on a probe mode of target sets
                                                                ! either all nb_types or probe/probed structure
integer,parameter       :: NB_PARAMS_MODE_LIKE_ONLY     = 1     ! only like nb_types except of probes
integer,parameter       :: NB_PARAMS_MODE_LIKE_ALL      = 2     ! only_like nb_types including probes
integer,parameter       :: NB_PARAMS_MODE_ALL           = 3     ! all nb types

! NBParamsRealms
integer,parameter       :: NB_PARAMS_REALMS_ALL         = 1     ! all realms
integer,parameter       :: NB_PARAMS_REALMS_ERA         = 2     ! ERA realms
integer,parameter       :: NB_PARAMS_REALMS_ER          = 3     ! ER realms

! === [control] ================================================================
integer                 :: NBParamsMode                 = NB_PARAMS_MODE_NORMAL     ! mode for determination of NB parameters
integer                 :: NBParamsRealms               = NB_PARAMS_REALMS_ER
integer                 :: NBCombRules                  = COMB_RULE_LB

! === [files] ==================================================================
character(len=MAX_PATH) :: InpParamFileName     = '-none-'          ! input parameters
character(len=MAX_PATH) :: OutParamFileName     = 'final.prms'      ! output parameters
character(len=MAX_PATH) :: OutAmberPrmsFileName = 'final.frcmod'    ! output Amber force field

! ------------------------------------------------------------------------------

end module ffdev_parameters_dat

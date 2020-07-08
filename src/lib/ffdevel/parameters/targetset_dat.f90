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

module ffdev_targetset_dat

use ffdev_geometry_dat
use ffdev_constants
use ffdev_variables
use ffdev_topology

! ------------------------------------------------------------------------------

! this training set associated with the topology
type TARGETSET
    type(TOPOLOGY)              :: top              ! set topology
    character(len=MAX_PATH)     :: final_stop       ! final name of topology
    character(len=MAX_PATH)     :: initial_drvene   ! driving profiles
    character(len=MAX_PATH)     :: initial_drvxyz   ! driving geometries
    character(len=MAX_PATH)     :: final_drvene     ! driving profiles
    character(len=MAX_PATH)     :: final_drvxyz     ! driving geometries
    integer                     :: ngeos            ! number of training points in the set
    type(GEOMETRY),pointer      :: geo(:)           ! training data
    integer                     :: mineneid         ! geometry with minimum of energy
    logical                     :: optgeo           ! optimize geometry
    logical                     :: keepoptgeo       ! keep optimized geometry
    logical                     :: savepts          ! save final pst
    logical                     :: savexyzr         ! save final xyzr
    logical                     :: nofreq           ! do not calculate frequencies when hessian is loaded
    character(len=MAX_PATH)     :: name             ! target name
    logical                     :: isref            ! is reference structure
    integer                     :: nrefs            ! number of references
    integer,pointer             :: refs(:)          ! references
    integer,pointer             :: natomsrefs(:)    ! number of atoms in reference sets for SAPT calculations
end type TARGETSET

integer                     :: nsets    ! number of sets
type(TARGETSET),allocatable :: sets(:)  ! all training sets

! [control] in {FFPARAMS}
real(DEVDP)                 :: max_probe_energy         = -1  ! -1 disable the filter

! [setup] ----------------------------------------------------------------------
logical                 :: DoNotCalcFreqs               = .false.   ! do not calculate frequencies when hessian is loaded
logical                 :: useOptGeometry               = .true.    ! save optimized geometry otherwise use target geometry
logical                 :: KeepPtsNames                 = .false.
logical                 :: SavePts                      = .false.
character(len=MAX_PATH) :: SavePtsPath                  = 'points'  ! storage for saved points
logical                 :: SaveXYZR                     = .false.
character(len=MAX_PATH) :: SaveXYZRPath                 = 'xyzr'
logical                 :: SaveSumLogs                  = .true.   ! for error statistics
character(len=MAX_PATH) :: SaveSumLogsPath              = '03.sumlogs'

! [optgeo] --------------------------------------------------------------------
! default values set in ffdev_targetset_ctrl_optgeo_set_default
logical                 :: GlbOptGeometryEnabled        ! force all geometry optimization in each error evaluation
logical                 :: GlbOptGeometryDisabled       ! disable all geometry optimization in each error evaluation
logical                 :: GlbShowOptProgress           ! print geometry optimization progress
logical                 :: GlbKeepOptGeometry           ! keep geometry from previous geometry optimization
integer                 :: ResetKeptOptGeoAt            ! reset initial geometry to trg for keepoptgeo = on

end module ffdev_targetset_dat

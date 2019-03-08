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
use ffdev_topology

! ------------------------------------------------------------------------------

! this training set associated with the topology
type TARGETSET
    type(TOPOLOGY)              :: top          ! set topology
    character(len=MAX_PATH)     :: final_stop   ! final name of topology
    character(len=MAX_PATH)     :: initial_drvene ! driving profiles
    character(len=MAX_PATH)     :: initial_drvxyz ! driving geometries
    character(len=MAX_PATH)     :: final_drvene ! driving profiles
    character(len=MAX_PATH)     :: final_drvxyz ! driving geometries
    integer                     :: ngeos        ! number of training points in the set
    type(GEOMETRY),pointer      :: geo(:)       ! training data
    real(DEVDP)                 :: offset       ! energy offset
    integer                     :: mineneid     ! geometry with minimum of energy
    logical                     :: optgeo       ! optimize geometry
    logical                     :: keepoptgeo   ! keep optimized geometry
    logical                     :: savegeo      ! save optimized geometry
    logical                     :: nofreq       ! do not calculate frequencies when hessian is loaded
    character(len=MAX_PATH)     :: name         ! target name
    logical                     :: isref        ! is reference structure
    integer                     :: nrefs        ! number of references
    integer,pointer             :: refs(:)      ! references
end type TARGETSET

integer                     :: nsets    ! number of sets
type(TARGETSET),allocatable :: sets(:)  ! all training sets

! [setup] ----------------------------------------------------------------------
logical                 :: OptimizeGeometry             = .false.   ! optimize geometry in each error evaluation
logical                 :: ShowOptimizationProgress     = .false.   ! print geometry optimization progress 
logical                 :: KeepOptimizedGeometry        = .true.    ! keep geometry from previous geometry optimization
logical                 :: SaveGeometry                 = .false.   
logical                 :: DoNotCalcFreqs               = .false.   ! do not calculate frequencies when hessian is loaded
character(len=MAX_PATH) :: SavePointsPath               = 'points'  ! storage for saved points

end module ffdev_targetset_dat

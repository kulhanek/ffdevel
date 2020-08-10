! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_nb2nb_dat

use ffdev_sizes
use ffdev_constants
use ffdev_variables
use ffdev_topology_dat

! ==============================================================================

! possible values for lj2exp6_alpha
real(DEVDP) :: lj2exp6_alpha    = 12.0d0    ! alpha for lj to exp-6 potential conversion
                                            ! 12.0                           - identical long-range
                                            ! 0.5d0*(19.0d0 + sqrt(73.0d0))  - identical shape in local minima

real(DEVDP) :: ljdefpb          = 3.0d0     ! default PB for r0=0
real(DEVDP) :: ljdefrc          = 2.6d0     ! default Rc


! === [nb2lj] ==================================================================
integer,parameter       :: NB2LJ_MODE_MINIMUM           = 1
integer,parameter       :: NB2LJ_MODE_OVERLAY           = 2
integer,parameter       :: NB2LJ_MODE_OVERLAY_REP       = 3
integer,parameter       :: NB2LJ_MODE_OVERLAY_DISP      = 4

integer                 :: NB2LJMode                    = NB2LJ_MODE_OVERLAY
logical                 :: NB2LJWeighted                = .false.
real(DEVDP)             :: NB2LJCutoffR                 = 10.0          ! max range for r
integer                 :: NB2LJIterGS                  = 1000          ! precision - GoldenSearch for r0, eps
integer                 :: NB2LJIterBS                  = 1000          ! precision - bisection for sigma
integer                 :: NB2LJIterOpt                 = 300           ! precision - overlay optimization via CMA-ES
real(DEVDP)             :: NB2LJSharkInitialStep        = 0.2           ! CMA-ES optimizer setup
real(DEVDP)             :: NB2LJTemp                    = 300.0         ! temp factor for weights
real(DEVDP)             :: NB2LJdr                      = 0.001         ! dr in partition function calculation, overlay calculation
real(DEVDP)             :: NB2LJdrPrint                 = 0.02          ! for printing
logical                 :: NB2LJCalcQNBIsoline          = .true.        ! add to NB pot also QNB isoline
real(DEVDP)             :: NB2LJCutoffRQNB              = 5.0           ! max r range for QNB isovalues
character(len=MAX_PATH) :: NBPotPathCore                = '04.nbpot'    ! NB potential storage
logical                 :: NB2LJIncludePen              = .true.

! working data for NB and overal calcs
integer                 :: NB2LJNParams
real(DEVDP)             :: NB2LJSigma
real(DEVDP)             :: NB2LJMinR
real(DEVDP)             :: NB2LJMaxR
type(NB_PAIR)           :: NB2LJNBPair
integer                 :: NB2LJErrFceEval
real(DEVDP),allocatable :: NB2LJprms(:),NB2LJtmp_xg(:),NB2LJtmp_ub(:),NB2LJtmp_lb(:)
! working data - QNB isoline
logical                 :: QNBModeEps = .true.
real(DEVDP)             :: QNBR0
real(DEVDP)             :: QNBEps
real(DEVDP)             :: QNBTrg
! working data
character(len=MAX_PATH) :: NBPotPathPrg

! global NB types

type GLBNB_TYPE
    integer             :: gti,gtj              ! atom types
    integer             :: setid                ! topology from set
    integer             :: nbt                  ! nb type from topology
    real(DEVDP)         :: eps, r0              ! LJ parameters
! NB parameters
    real(DEVDP)         :: SigNB                ! sigma for NB
    real(DEVDP)         :: R0NB                 ! r0 for NB
    real(DEVDP)         :: EpsNB                ! r0 for NB
    real(DEVDP)         :: QNB                  ! partition function for NB
    integer             :: num
    real(DEVDP)         :: errval
end type GLBNB_TYPE

integer                         :: nnb_types
type(GLBNB_TYPE),allocatable    :: nb_types(:)

! ------------------------------------------------------------------------------

end module ffdev_nb2nb_dat

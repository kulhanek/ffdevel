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

! === [NB2NB] ==================================================================
logical                 :: NB2NBLikeOnly                = .false.
real(DEVDP)             :: NB2NBTemp                    = 300.0         ! temp factor for weights
real(DEVDP)             :: NB2NBCutoffR                 = 10.0          ! max range for r
integer                 :: NB2NBNBins                   = 1000          ! discretion of NB potential
integer                 :: NB2NBIterGS                  = 1000          ! precision - GoldenSearch for r0, eps
integer                 :: NB2NBIterBS                  = 1000          ! precision - bisection for sigma

real(DEVDP)             :: NB2NBdrPrint                 = 0.02          ! for printing

character(len=MAX_PATH) :: NBPotPathCore                = '04.nbpot'    ! NB potential storage

logical                 :: NB2NBIncludePen              = .true.
logical                 :: NB2NBIncludeInd              = .true.

! QNBIsoline
logical                 :: NB2NBCalcQNBIsoline          = .true.        ! add to NB pot also QNB isoline

! ------------------------------------------------------------------------------

! internal setup
real(DEVDP)             :: NB2NBSharkInitialStep        = 0.2           ! CMA-ES optimizer setup
integer                 :: NB2NBIterOpt                 = 300           ! precision - isoline optimization via CMA-ES

! Gaussian quadrature for QNB calculation
logical                 :: NB2NBUseGaussQuad            = .true.
integer                 :: NB2NBGaussQuadOrder          = 60
real(DEVDP),allocatable :: NB2NBGaussQuadA(:)
real(DEVDP),allocatable :: NB2NBGaussQuadW(:)

! working data for NB and overal calcs
integer                 :: NB2NBNParams
real(DEVDP)             :: NB2NBSigma
real(DEVDP)             :: NB2NBMinR
real(DEVDP)             :: NB2NBMaxR
type(NB_PAIR)           :: NB2NBNBPair
integer                 :: NB2NBErrFceEval
real(DEVDP),allocatable :: NB2NBprms(:),NB2NBtmp_xg(:),NB2NBtmp_ub(:),NB2NBtmp_lb(:)

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
! NB parameters
    real(DEVDP)         :: SigNB                ! sigma for NB
    real(DEVDP)         :: R0NB                 ! r0 for NB
    real(DEVDP)         :: EpsNB                ! r0 for NB
    real(DEVDP)         :: QNB                  ! partition function for NB
    real(DEVDP),pointer :: NBPot(:)             ! NB potential from SigNB to SigNB+NB2NBCutoffR, with NB2NBdr step
    integer             :: num                  ! number of NB pair occurrences in all topologies
    real(DEVDP)         :: errval
end type GLBNB_TYPE

integer                         :: nnb_types
type(GLBNB_TYPE),allocatable    :: nb_types(:)

! ------------------------------------------------------------------------------

end module ffdev_nb2nb_dat

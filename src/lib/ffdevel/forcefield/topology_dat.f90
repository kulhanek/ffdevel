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

module ffdev_topology_dat

use ffdev_sizes
use ffdev_constants
use ffdev_variables

! ------------------------------------------------------------------------------

type ATOM
    real(DEVDP)             :: charge
    character(MAX_TNAME)    :: name
    integer                 :: residx
    character(MAX_RNAME)    :: resname
    integer                 :: typeid
    integer                 :: nbonds
    integer,pointer         :: bonded(:)
    integer                 :: frgid        ! fragment ID for NB error function
end type ATOM

! ------------------------------------------------------------------------------

type BOND
    integer ::  ai,aj   ! bonded atoms
    integer ::  bt      ! bond type
end type BOND

! ------------------------------------------------------------------------------

type BOND_TYPE
    integer     :: ti, tj   ! atom types
    integer     :: model    ! bond model type
    real(DEVDP) :: d0       ! equilibrium distance
    real(DEVDP) :: k        ! force constant
    logical     :: ffoptactive      ! this type is subject of ffopt
end type BOND_TYPE

! ------------------------------------------------------------------------------

type ANGLE
    integer ::  ai,aj,ak    ! bonded atoms
    integer ::  at          ! angle type
end type ANGLE

! ------------------------------------------------------------------------------

type ANGLE_TYPE
    integer     :: ti, tj, tk   ! atom types
    integer     :: model        ! angle model type
    real(DEVDP) :: a0           ! equilibrium angle
    real(DEVDP) :: k            ! force constant
    logical     :: ffoptactive      ! this type is subject of ffopt
end type ANGLE_TYPE

! ------------------------------------------------------------------------------

type DIHEDRAL
    integer ::  ai,aj,ak,al ! bonded atoms
    integer ::  dt          ! dihedral type
end type DIHEDRAL

! ------------------------------------------------------------------------------

type DIHEDRAL_TYPE
    integer             :: ti, tj, tk, tl   ! atom types
    integer             :: n                ! sequence size
    real(DEVDP)         :: inv_scee         ! 1 / 1-4 electrostatics scaling factor
    real(DEVDP)         :: inv_scnb         ! 1 / 1-4 NB scaling factor
    integer             :: mode             ! dihedral mode, 1 - cos, 2 - grbf
    real(DEVDP),pointer :: v(:)             ! cos - energies
    real(DEVDP),pointer :: g(:)             ! cos - phases
    real(DEVDP),pointer :: c(:)             ! grbf - weights
    real(DEVDP),pointer :: p(:)             ! grbf - positions
    real(DEVDP),pointer :: w2(:)            ! grbf - widths
    logical,pointer     :: enabled(:)       ! what is enabled, apllicable both to cos and grbf
    logical             :: ffoptactive      ! this type is subject of ffopt
end type DIHEDRAL_TYPE

integer,parameter       :: DIH_COS      = 1
integer,parameter       :: DIH_GRBF     = 2

! ------------------------------------------------------------------------------

type IMPROPER_TYPE
    integer     :: ti, tj, tk, tl   ! atom types
    real(DEVDP) :: v
    real(DEVDP) :: g
    logical     :: ffoptactive      ! this type is subject of ffopt
end type IMPROPER_TYPE

! ------------------------------------------------------------------------------

type NB_PAIR
    integer     ::  ai,aj       ! pair
    integer     ::  nbt         ! NB type
    integer     ::  nbtii,nbtjj ! NB types of correponding like atoms
    integer     ::  dt          ! dihedral type for scaling factors
    integer     ::  ci          ! index to grid cache
end type NB_PAIR

! ------------------------------------------------------------------------------

type ATOM_TYPE
    character(MAX_TNAME)    :: name
    real(DEVDP)             :: mass
    integer                 :: z
    logical                 :: probe            ! probe
    integer                 :: glbtypeid        ! global type id
end type ATOM_TYPE

! ------------------------------------------------------------------------------

type NB_TYPE
    integer             :: ti,tj                    ! atom types
    real(DEVDP)         :: eps, r0, alpha           ! ERA - non-bonded parameters (LJ/EXP6 vdW parameters)
    real(DEVDP)         :: PA, PB, C6               ! ABC - non-bonded parameters
    logical             :: ffoptactive              ! this type is subject of ffopt
end type NB_TYPE

! ------------------------------------------------------------------------------

type TOPOLOGY
! general info
    character(len=255)          :: name
    integer                     :: natoms
    type(ATOM),pointer          :: atoms(:)

! bonded terms
    integer                     :: nbonds
    type(BOND),pointer          :: bonds(:)
    integer                     :: nbond_types
    type(BOND_TYPE),pointer     :: bond_types(:)

    integer                     :: nangles
    type(ANGLE),pointer         :: angles(:)
    integer                     :: nangle_types
    type(ANGLE_TYPE),pointer    :: angle_types(:)

    integer                     :: ndihedrals
    type(DIHEDRAL),pointer      :: dihedrals(:)
    integer                     :: ndihedral_types
    integer                     :: ndihedral_seq_size
    type(DIHEDRAL_TYPE),pointer :: dihedral_types(:)

    integer                     :: nimpropers
    type(DIHEDRAL),pointer      :: impropers(:)
    integer                     :: nimproper_types
    type(IMPROPER_TYPE),pointer :: improper_types(:)

! non-bonded terms
    integer                     :: nb_size
    ! integer                     :: nb_mode
    type(NB_PAIR),pointer       :: nb_list(:)
    integer                     :: natom_types
    type(ATOM_TYPE),pointer     :: atom_types(:)
    integer                     :: nnb_types
    type(NB_TYPE),pointer       :: nb_types(:)
    integer                     :: probe_size           ! number of atoms in probe
    integer                     :: nfragments
end type TOPOLOGY

! ------------------------------------------------------------------------------
! global options in

logical     :: dih_cos_only     = .false.       ! .true. -> SUM Vn*cos(n*phi-gamma)
                                                ! .false. -> SUM Vn*(1+cos(n*phi-gamma))

! ==============================================================================
! ==== electrostatics
! ==============================================================================

integer,parameter   :: NB_ELE_QTOP   = 1        ! charges from topology
integer,parameter   :: NB_ELE_QGEO   = 2        ! charges from geometries

integer     :: ele_mode     = NB_ELE_QTOP
real(DEVDP) :: ele_qscale   = 1.0d0             ! scaling factor for charges

! ==============================================================================
! ==== vdW modes
! ==============================================================================

! combining rules for nb_mode
integer,parameter   :: COMB_RULE_NONE = 05    ! input data

! ####################################################################
integer,parameter   :: NB_VDW_LJ        = 1
! Lenard-Jones
! Form: Enb = eps*( (ro/r)^12 - 2*(r0/r)^6 )
! Parameters: eps, r0
! Provides: energy, gradient, Hessian

! combining rules - applicable for NB_MODE_LJ
integer,parameter   :: COMB_RULE_LB = 11    ! LB (Lorentz-Berthelot)
integer,parameter   :: COMB_RULE_WH = 12    ! WH (Waldman-Hagler)
integer,parameter   :: COMB_RULE_KG = 13    ! KG (Kong)
integer,parameter   :: COMB_RULE_FB = 14    ! FB (Fender-Halsey-Berthelot)

! ####################################################################
integer,parameter   :: NB_VDW_12_6      = 2
! Lenard-Jones
! Form: Enb = exp(PA)/r^12 - disp_fa*C6/r^6
! Parameters: PA, C6, disp_fa
! Provides: energy, gradient

! combining rules - applicable for NB_VDW_12_XDMC6
integer,parameter   :: COMB_RULE_PA1 = 15   ! geometric mean:  exp(PAIJ) = sqrt(exp(PAII)*exp(PAJJ)), PAIJ = (PAII+PAJJ)/2

! these are also applicable via backward/forward transformations
!integer,parameter   :: COMB_RULE_LB = 11    ! LB (Lorentz-Berthelot)
!integer,parameter   :: COMB_RULE_WH = 12    ! WH (Waldman-Hagler)
!integer,parameter   :: COMB_RULE_KG = 13    ! KG (Kong)
!integer,parameter   :: COMB_RULE_FB = 14    ! FB (Fender-Halsey-Berthelot)


! ####################################################################
integer,parameter   :: NB_VDW_12_XDMBJ  = 3
! Lenard-Jones
! Form: Enb = exp(PA)/r^12 - XDMC6/(r^6 + rbj^6) - XDMC8/(r^8 + rbj^8) - XDMC10/(r^10 + rbj^10)
! Parameters: PA, disp_fa, disp_fb
! Provides: energy, gradient

! rbj = disp_fa*rc + disp_fb

! combining rules - applicable for NB_VDW_12_XDMC6
! integer,parameter   :: COMB_RULE_PA1 = 15   ! geometric mean:  exp(PAIJ) = sqrt(exp(PAII)*exp(PAJJ)), PAIJ = (PAII+PAJJ)/2

! ####################################################################
integer,parameter   :: NB_VDW_TT_XDM    = 4     ! Tang–Toennies + XDM

! Tang–Toennis
! Form: Enb = exp(PA)*exp(-PB*r) - fd6*XDM_C6/r^6 - fd8*XDM_C8/r^8 - fd10*XDM_C10/r^10
! Parameters: PA, PB
! Provides: energy

! combining rules - applicable for NB_VDW_TT_XDM
integer,parameter   :: COMB_RULE_TT1 = 17   ! geometric mean:  exp(PAIJ) = sqrt(exp(PAII)*exp(PAJJ))
                                            ! arithmetic mean: PBIJ=(PBII+PBJJ)/2
integer,parameter   :: COMB_RULE_TT2 = 18

! ==============================================================================

! employed method for NB interactions
integer     :: nb_mode          = NB_VDW_LJ
integer     :: nb_comb_rules    = COMB_RULE_LB

! ------------------------------------------------------------------------------

! tuneable parameters
real(DEVDP) :: disp_fa = 1.0d0
real(DEVDP) :: disp_fb = 0.0d0

! ------------------------------------------------------------------------------

end module ffdev_topology_dat

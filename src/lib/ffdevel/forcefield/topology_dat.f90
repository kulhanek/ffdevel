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

! ------------------------------------------------------------------------------

type ATOM
    real(DEVDP)         :: charge
    character(len=4)    :: name
    integer             :: residx
    character(len=4)    :: resname
    integer             :: typeid
    integer             :: nbonds
    integer,pointer     :: bonded(:)
    integer             :: frgid        ! fragment ID for NB error function
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
    character(len=4)    :: name
    real(DEVDP)         :: mass
    integer             :: z
    logical             :: probe            ! probe
    integer             :: glbtypeid        ! global type id
end type ATOM_TYPE

! ------------------------------------------------------------------------------

type NB_TYPE
    integer             :: ti,tj                    ! atom types
    real(DEVDP)         :: eps, r0, alpha           ! LJ/EXP6 vdW parameters
    real(DEVDP)         :: PA, PB, PC, PD, PR       ! Pauli repulsion data
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
    integer                     :: nb_mode
    type(NB_PAIR),pointer       :: nb_list(:)
    integer                     :: natom_types
    type(ATOM_TYPE),pointer     :: atom_types(:)
    integer                     :: nnb_types
    type(NB_TYPE),pointer       :: nb_types(:)
    integer                     :: probe_size           ! number of atoms in probe
    integer                     :: nfragments
    ! assumed combination rules, depending on parameter optimization strategy
    ! these rules can be easily broken !!!!
    integer                     :: assumed_comb_rules
end type TOPOLOGY

! ------------------------------------------------------------------------------
! global options

logical     :: dih_cos_only = .false.   ! .true. -> SUM Vn*cos(n*phi-gamma)
                                        ! .false. -> SUM Vn*(1+cos(n*phi-gamma))

! possible values for lj2exp6_alpha
! 12.0                           - identical long-range
! 0.5d0*(19.0d0 + sqrt(73.0d0))  - identical shape in local minima
real(DEVDP) :: lj2exp6_alpha = 12.0d0   ! alpha for lj to exp-6 potential conversion

logical     :: cache_grid = .true.     ! cache grid for Eexchrep calculation

! ------------------------------------------------------------------------------

! full mode
integer,parameter   :: NB_MODE_LJ       = 1     ! Lennard-Jones potential (eps,r0)
integer,parameter   :: NB_MODE_EXP6     = 2     ! Exp-6 potential (eps,r0,alpha)

! exchange repulsion only (Pauli repulsion), gradient and hessian is not available
integer,parameter   :: NB_MODE_PAULI_DENS2  = 5     ! overlap of densities with two parameters (A,B)
integer,parameter   :: NB_MODE_PAULI_DENS3  = 6     ! overlap of densities with three parameters (A,B,C)
integer,parameter   :: NB_MODE_PAULI_DENS3P = 10    ! overlap of densities with three parameters (A,B,C - polynomial)
integer,parameter   :: NB_MODE_PAULI_DENS5  = 9     ! overlap of densities with five parameters (A,B,C,D,R)
integer,parameter   :: NB_MODE_PAULI_WAVE2  = 7     ! overlap of wavefunctions with two parameters (A,B)
integer,parameter   :: NB_MODE_PAULI_WAVE3  = 8     ! overlap of wavefunctions with three parameters (A,B,C)

! combining rules for nb_mode
integer,parameter   :: COMB_RULE_IN = 05  ! input data
! applicable for NB_MODE_LJ
integer,parameter   :: COMB_RULE_LB = 10  ! LB (Lorentz-Berthelot)
integer,parameter   :: COMB_RULE_WH = 20  ! WH (Waldman-Hagler)
integer,parameter   :: COMB_RULE_KG = 30  ! KG (Kong)
integer,parameter   :: COMB_RULE_FB = 40  ! FB (Fender-Halsey-Berthelot)

! ------------------------------------------------------------------------------

end module ffdev_topology_dat

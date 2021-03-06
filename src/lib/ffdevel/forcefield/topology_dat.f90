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
    integer                 :: frgid            ! fragment ID for NB error function
    integer                 :: symmclass
    logical                 :: ffoptactive      ! this type is subject of ffopt
    integer                 :: chrg_prm_id      ! helper variable
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
    integer     ::  nbtii,nbtjj ! NB types of corresponding like atoms
    integer     ::  dt          ! dihedral type for scaling factors
    ! acceleration data -  do not use for regular use
    ! these data are updated when nb_params_update .eqv. .true.      !
        real(DEVDP) ::  pa1,pa2
        real(DEVDP) ::  pb1,pb2
        real(DEVDP) ::  c6,c8,c10
        real(DEVDP) ::  rc6,rc8,rc10
        real(DEVDP) ::  Z1,Z2
        real(DEVDP) ::  q1,q2
end type NB_PAIR

! ------------------------------------------------------------------------------

type NB_PAIR_ENERGY
        real(DEVDP) :: tot_ene
        real(DEVDP) :: ele_ene
        real(DEVDP) :: pen_ene
        real(DEVDP) :: ind_ene
        real(DEVDP) :: rep_ene
        real(DEVDP) :: dis_ene
end type NB_PAIR_ENERGY

! ------------------------------------------------------------------------------

type ATOM_TYPE
    character(MAX_TNAME)    :: name
    real(DEVDP)             :: mass
    integer                 :: z
    logical                 :: probe            ! probe
    integer                 :: glbtypeid        ! global type id
    logical                 :: ffoptactive      ! this type is subject of ffopt
    ! NB parameters
    real(DEVDP)             :: PA
    real(DEVDP)             :: PB
    real(DEVDP)             :: RC               ! BJ damping radius
end type ATOM_TYPE

! ------------------------------------------------------------------------------

type NB_TYPE
    integer             :: ti,tj                ! atom types
    real(DEVDP)         :: eps, r0, alpha       ! LJ parameters
    logical             :: ffoptactive          ! this type is subject of ffopt
    integer             :: nbii, nbjj           ! links to like nb_types
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
    type(NB_PAIR),pointer       :: nb_list(:)
    integer                     :: natom_types
    type(ATOM_TYPE),pointer     :: atom_types(:)
    integer                     :: nnb_types
    type(NB_TYPE),pointer       :: nb_types(:)
    integer                     :: probe_size           ! number of atoms in probe
    integer                     :: nfragments
    integer                     :: nsymm_classes
    integer                     :: total_charge
    ! this update  nb parameters
    logical                     :: nb_params_update

! SAPT
    integer                     :: sapt_size
    type(NB_PAIR),pointer       :: sapt_list(:)

end type TOPOLOGY

! ------------------------------------------------------------------------------
! global options in

logical     :: dih_cos_only     = .false.       ! .true. -> SUM Vn*cos(n*phi-gamma)
                                                ! .false. -> SUM Vn*(1+cos(n*phi-gamma))

! ==============================================================================
! ==== electrostatics
! ==============================================================================

real(DEVDP) :: ele_qscale                   = 1.0d0     ! scaling factor for charges
real(DEVDP) :: glb_iscee                    = 1.0d0     ! global scaling for 1-4 electrostatics

! ==============================================================================
! ==== vdW modes
! ==============================================================================

! penetration energy
integer,parameter   :: PEN_MODE_SCI         = 2         ! Smear charge interaction

! ####################################################################

! induction energy
integer,parameter   :: IND_MODE_K2EXC       = 802       ! proportional to exchange energy

! ####################################################################

! combining rules - applicable for NB_MODE_LJ
integer,parameter   :: LJ_COMB_RULE_LB      = 11        ! LB (Lorentz-Berthelot)
integer,parameter   :: LJ_COMB_RULE_WH      = 12        ! WH (Waldman-Hagler)
integer,parameter   :: LJ_COMB_RULE_KG      = 13        ! KG (Kong)
integer,parameter   :: LJ_COMB_RULE_FH      = 14        ! FB (Fender-Halsey)

integer,parameter   :: LJ_ALPHA_CONST       = 471       ! Exp6 probe mode - constant alpha
integer,parameter   :: LJ_ALPHA_FREEOPT     = 472       ! Exp6 probe mode - free to optimize per type

! ####################################################################

! exp_mode
integer,parameter   :: EXP_MODE_DO          = 67        ! Density overlap
integer,parameter   :: EXP_MODE_WO          = 68        ! Wavefunction overlap

! PA source
integer,parameter   :: EXP_PA_FREEOPT       = 301       ! pa = free to optimize as vdw_pa
integer,parameter   :: EXP_PA_CHARGES       = 302       ! from charges with k_exc prefactor

! PB source
integer,parameter   :: EXP_PB_FREEOPT       = 301       ! pb = free to optimize as vdw_pb
integer,parameter   :: EXP_PB_ADBII         = 302       ! pb = damp_pb * bii (atomic database)

! ####################################################################

! Becke-Johnson damping
integer,parameter   :: DAMP_BJ_CONST        = 201       ! rc = damp_fa
integer,parameter   :: DAMP_BJ_FREEOPT      = 202       ! rc = free to optimize as vdw_rc
integer,parameter   :: DAMP_BJ_ADRCII       = 205       ! rc = rcii(typeA,typeB,[damp_fa,damp_fb]) (atomic database)

! ==============================================================================
! nb_mode
integer,parameter   :: NB_VDW_LJ            = 13        ! Lenard-Jones
integer,parameter   :: NB_VDW_EXP_DISPBJ    = 21        ! Exp-Becke-Johnson
integer,parameter   :: NB_VDW_EXP_DISPTT    = 30        ! Exp-Tang–Toennies

! employed method for NB interactions
integer     :: nb_mode                      = NB_VDW_LJ
! PEN
logical     :: pen_enabled                  = .false.           ! penetration energy
integer     :: pen_mode                     = PEN_MODE_SCI

! IND
logical     :: ind_enabled                  = .false.           ! induction energy
integer     :: ind_mode                     = IND_MODE_K2EXC

! LJ
integer     :: lj_comb_rules                = LJ_COMB_RULE_LB
logical     :: lj_exp6_probe                = .false.
integer     :: lj_alpha_mode                = LJ_ALPHA_CONST

! EXP
integer     :: exp_mode                     = EXP_MODE_WO
integer     :: exp_pa_mode                  = EXP_PA_FREEOPT
integer     :: exp_pb_mode                  = EXP_PB_FREEOPT

! DISP
integer     :: dampbj_mode                  = DAMP_BJ_FREEOPT

! RUNTIME setup
logical     :: ApplyCombiningRules          = .false.       ! apply combination rules in every error evaluation
! ------------------------------------------------------------------------------

real(DEVDP) :: glb_iscnb                    = 1.0d0         ! global scaling for 1-4 NB interactions

! exp6 probe mode
real(DEVDP) :: exp6_alpha = 15.5d0

! tuneable parameters
real(DEVDP) :: disp_s6  =  1.0d0
real(DEVDP) :: disp_s8  =  1.0d0
real(DEVDP) :: disp_s10 =  1.0d0

! damping
real(DEVDP) :: damp_fa  =  1.0d0    ! BJ
real(DEVDP) :: damp_fb  =  0.0d0    ! BJ
real(DEVDP) :: damp_pb  =  1.0d0    ! exchange
real(DEVDP) :: damp_tb  =  1.0d0    ! TT
real(DEVDP) :: damp_pe  =  1.0d0    ! PEN

real(DEVDP) :: k_exc    =  1.0d0
real(DEVDP) :: k_ind    =  1.0d0

! ------------------------------------------------------------------------------

end module ffdev_topology_dat

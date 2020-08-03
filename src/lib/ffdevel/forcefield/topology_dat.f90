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
    ! these data are updated when nb_params_update .eqv. .true.
        real(DEVDP) ::  crgij       !
        real(DEVDP) ::  pa
        real(DEVDP) ::  pb
        real(DEVDP) ::  c6
        real(DEVDP) ::  c8
        real(DEVDP) ::  c10
        real(DEVDP) ::  rc6
        real(DEVDP) ::  rc8
        real(DEVDP) ::  rc10
        real(DEVDP) ::  tb
    ! penetration energy
        real(DEVDP) ::  Z1,Z2
        real(DEVDP) ::  q1,q2
        real(DEVDP) ::  pepa1,pepb1
        real(DEVDP) ::  pepa2,pepb2
end type NB_PAIR

! ------------------------------------------------------------------------------

type ATOM_TYPE
    character(MAX_TNAME)    :: name
    real(DEVDP)             :: mass
    integer                 :: z
    logical                 :: probe            ! probe
    integer                 :: glbtypeid        ! global type id
    logical                 :: ffoptactive      ! this type is subject of ffopt
    ! penetration energy
    real(DEVDP)             :: pen_pa
    real(DEVDP)             :: pen_pb
end type ATOM_TYPE

! ------------------------------------------------------------------------------

type NB_TYPE
    integer             :: ti,tj                ! atom types
    real(DEVDP)         :: eps, r0              ! LJ parameters
    real(DEVDP)         :: PA, PB               ! repulsion parameters
    real(DEVDP)         :: RC                   ! BJ damping radius
    real(DEVDP)         :: TB                   ! TT damping factor
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

integer,parameter   :: PEN_PA_FREEOPT       = 784        ! test mode
integer,parameter   :: PEN_PA_CONST         = 785        ! test mode
integer,parameter   :: PEN_PA_LDO           = 786        ! test mode
integer,parameter   :: PEN_PA_DO            = 787        ! test mode

integer,parameter   :: PEN_PB_FREEOPT       = 794        ! test mode
integer,parameter   :: PEN_PB_CONST         = 795        ! test mode
integer,parameter   :: PEN_PB_COUPLED       = 796        ! test mode

! ####################################################################

! combining rules - applicable for NB_MODE_LJ
integer,parameter   :: LJ_COMB_RULE_LB      = 11        ! LB (Lorentz-Berthelot)
integer,parameter   :: LJ_COMB_RULE_WH      = 12        ! WH (Waldman-Hagler)
integer,parameter   :: LJ_COMB_RULE_KG      = 13        ! KG (Kong)
integer,parameter   :: LJ_COMB_RULE_FB      = 14        ! FB (Fender-Halsey-Berthelot)

! ####################################################################

! exp_mode
integer,parameter   :: EXP_BM               = 64        ! Born-Mayer
integer,parameter   :: EXP_DO               = 67        ! Density overlap
integer,parameter   :: EXP_WO               = 68        ! Wavefunction overlap

! PB source
integer,parameter   :: EXP_PB_FREEOPT       = 301       ! PB free to optimize
integer,parameter   :: EXP_PB_DO            = 302       ! PB from density overlap
integer,parameter   :: EXP_PB_IP            = 304       ! PB from ionization potential
integer,parameter   :: EXP_PB_IP_XDM        = 305       ! PB from ionization potential + XDM mod

! combining rules - applicable for Born-Mayer repulsion
integer,parameter   :: EXP_COMB_RULE_AM     = 391       ! arithmetic means
integer,parameter   :: EXP_COMB_RULE_GS     = 392       ! Gilbert-Smith
integer,parameter   :: EXP_COMB_RULE_BA     = 393       ! Bohm-Ahlrichs
integer,parameter   :: EXP_COMB_RULE_VS     = 394       ! Vleet-Schmidt
integer,parameter   :: EXP_COMB_RULE_D1     = 395       ! test
integer,parameter   :: EXP_COMB_RULE_D2     = 396       ! test
integer,parameter   :: EXP_COMB_RULE_W1     = 397       ! test
integer,parameter   :: EXP_COMB_RULE_W2     = 398       ! test

! ####################################################################

! Becke-Johnson damping
integer,parameter   :: DAMP_BJ_CONST        = 201       ! constant
integer,parameter   :: DAMP_BJ_FREEOPT      = 202       ! Rc free to optimize
integer,parameter   :: DAMP_BJ_DRC          = 203       ! radii from Cx
integer,parameter   :: DAMP_BJ_DO           = 204       ! derived from density overlaps

! Tang–Toennies damping
integer,parameter   :: DAMP_TT_COUPLED      = 101   ! tb = damp_tb * pb
integer,parameter   :: DAMP_TT_FREEOPT      = 102   ! tb free to optimize
integer,parameter   :: DAMP_TT_CONST        = 103   ! tb = damp_tb
integer,parameter   :: DAMP_TT_DO           = 104   ! tb = damp_tb * densoverlap_bii
integer,parameter   :: DAMP_TT_WO           = 105   ! tb = damp_tb * wfoverlap_bii
integer,parameter   :: DAMP_TT_IP           = 106   ! tb = damp_tb * f(ionization potential)
integer,parameter   :: DAMP_TT_IP_XDM       = 107   ! tb = damp_tb * f(ionization potential) * XDM mod

! ==============================================================================
! nb_mode
integer,parameter   :: NB_VDW_LJ            = 13    ! Lenard-Jones
integer,parameter   :: NB_VDW_EXP_DISPBJ    = 21    ! Exp-Becke-Johnson
integer,parameter   :: NB_VDW_EXP_DISPTT    = 30    ! Exp-Tang–Toennies

! employed method for NB interactions
integer     :: nb_mode                      = NB_VDW_LJ
! PEN
logical     :: pen_enabled                  = .false.           ! penetration energy
integer     :: pen_pa_mode                  = PEN_PA_FREEOPT
integer     :: pen_pb_mode                  = PEN_PB_COUPLED

! LJ
integer     :: lj_comb_rules                = LJ_COMB_RULE_LB
! EXP
integer     :: exp_mode                     = EXP_BM
integer     :: exp_pb_mode                  = EXP_PB_FREEOPT
integer     :: exp_comb_rules               = EXP_COMB_RULE_VS
! DISP
integer     :: dampbj_mode                  = DAMP_BJ_FREEOPT
integer     :: damptt_mode                  = DAMP_TT_FREEOPT

! derived setup
logical     :: ApplyCombiningRules          = .false.           ! apply combination rules in every error evaluation
! ------------------------------------------------------------------------------

real(DEVDP) :: glb_iscnb                    = 1.0d0         ! global scaling for 1-4 NB interactions

! tuneable parameters
real(DEVDP) :: disp_s6  =  1.0d0
real(DEVDP) :: disp_s8  =  1.0d0
real(DEVDP) :: disp_s10 =  1.0d0

! damping
real(DEVDP) :: damp_fa  =  1.0d0
real(DEVDP) :: damp_fb  =  0.0d0
real(DEVDP) :: damp_pb  =  1.0d0
real(DEVDP) :: damp_tb  =  1.0d0

! Pauli repulsion K factor for PAPNL
real(DEVDP) :: pauli_k  =  1.0d0

! penetration energy
real(DEVDP) :: pen_fa   =  1.0d0
real(DEVDP) :: pen_fb   =  1.0d0
real(DEVDP) :: pen_fc   = -5.0d0

! ------------------------------------------------------------------------------

end module ffdev_topology_dat

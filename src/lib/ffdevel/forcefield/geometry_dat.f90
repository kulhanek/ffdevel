! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_geometry_dat

use ffdev_sizes

! ------------------------------------------------------------------------------

type RESTRAINT
    real(DEVDP)                 :: trg_value
    real(DEVDP)                 :: value
    character(len=MAX_CVTYPE)   :: cvtype
    integer,pointer             :: ai(:)
end type RESTRAINT

! ------------------------------------------------------------------------------

integer,parameter   :: GEO_PROBE_ENE_HF     = 1  ! ele+pen+ind+exc
integer,parameter   :: GEO_PROBE_ENE_TOT    = 2  ! ele+pen+ind+exc+disp

! ------------------------------------------------------------------------------

type GEOMETRY
    integer                 :: id
    integer                 :: natoms           ! number of atoms
    character(len=MAX_PATH) :: name
    character(len=MAX_PATH) :: title
    real(DEVDP)             :: weight           ! contribution to all data
! geometry data
    integer,pointer         :: z(:)             ! proton numbers
    real(DEVDP),pointer     :: crd(:,:)         ! geometry (3,natoms)
! energies
    real(DEVDP)             :: bond_ene         ! bond energy
    real(DEVDP)             :: angle_ene        ! angle energy
    real(DEVDP)             :: dih_ene          ! dihedral energy
    real(DEVDP)             :: impropr_ene      ! improper energy
    real(DEVDP)             :: ele14_ene        ! 1-4 electrostatics
    real(DEVDP)             :: rep14_ene        ! 1-4 vdW - repulsion
    real(DEVDP)             :: dis14_ene        ! 1-4 vdW - dispersion
    real(DEVDP)             :: ele_ene          ! electrostatics
    real(DEVDP)             :: pen_ene          ! penetration energy (electrostatics)
    real(DEVDP)             :: ind_ene          ! induction
    real(DEVDP)             :: rep_ene          ! vdW - repulsion
    real(DEVDP)             :: dis_ene          ! vdW - dispersion
    real(DEVDP)             :: total_ene        ! total energy from MM
! sapt energies
    real(DEVDP)             :: sapt_ele         ! electrostatics
    real(DEVDP)             :: sapt_pen         ! penetration energy
    real(DEVDP)             :: sapt_ind         ! induction
    real(DEVDP)             :: sapt_rep         ! NB repulsion
    real(DEVDP)             :: sapt_dis         ! NB dispersion
    real(DEVDP)             :: sapt_total
! derivatives, etc.
    real(DEVDP),pointer     :: grd(:,:)         ! gradient (3,natoms)
    real(DEVDP),pointer     :: hess(:,:,:,:)    ! hessian (3,natoms,3,natoms) or normal modes
    real(DEVDP),pointer     :: ihess(:)         ! hessian in internal coordinates (na+nb+nd+ni)
    real(DEVDP),pointer     :: nmodes(:,:,:,:)  ! normal modes (3,natoms,3,natoms)
    real(DEVDP),pointer     :: freq(:)          ! frequencies of normal vibrations (3*natoms)
! target data
    logical                 :: trg_crd_optimized    ! indicates that the goemtry was optimized
    logical                 :: trg_ene_loaded
    logical                 :: trg_ene_generic
    logical                 :: trg_crd_loaded
    logical                 :: trg_grd_loaded
    logical                 :: trg_hess_loaded
    logical                 :: trg_ihess_loaded
    logical                 :: trg_freq_loaded
    logical                 :: trg_esp_loaded
    real(DEVDP)             :: trg_energy           ! target energy
    real(DEVDP),pointer     :: trg_crd(:,:)         ! target geometry (3,natoms)
    real(DEVDP),pointer     :: trg_grd(:,:)         ! target gradient (3,natoms)
    real(DEVDP),pointer     :: trg_hess(:,:,:,:)    ! target hessian (3,natoms,3,natoms)
    real(DEVDP),pointer     :: trg_ihess(:)         ! target hessian in internal coordinates (nb+na+nd+ni)
    real(DEVDP),pointer     :: trg_nmodes(:,:,:,:)  ! target normal modes (3,natoms,3,natoms)
    real(DEVDP),pointer     :: trg_freq(:)          ! target frequencies (3*natoms)
    integer                 :: esp_npoints          ! number of ESP points
    real(DEVDP),pointer     :: trg_esp(:,:)         ! target ESP (4,npoints)
! probe energy
    logical                 :: trg_probe_ene_loaded
    logical                 :: trg_probe_ene_generic
    real(DEVDP)             :: trg_probe_ene
    integer                 :: trg_probe_ene_mode
! target SAPT
    logical                 :: trg_sapt_loaded
    logical                 :: trg_sapt_generic
    real(DEVDP)             :: trg_sapt_ele
    real(DEVDP)             :: trg_sapt_exc
    real(DEVDP)             :: trg_sapt_ind
    real(DEVDP)             :: trg_sapt_dis
! restraints
    integer                 :: nrst                 ! number of restraints
    type(RESTRAINT),pointer :: rst(:)               ! restraints
    real(DEVDP)             :: rst_energy           ! colvar energy penalty
! supplemental data
    logical                 :: sup_xdm_loaded       ! all in a.u.
    real(DEVDP),pointer     :: sup_xdm_c6(:,:)      ! XDM C6 dispersion coefficients
    real(DEVDP),pointer     :: sup_xdm_c8(:,:)
    real(DEVDP),pointer     :: sup_xdm_c10(:,:)
    real(DEVDP),pointer     :: sup_xdm_vol(:)       ! atom volume
    real(DEVDP),pointer     :: sup_xdm_vol0(:)      ! atom free volume
    real(DEVDP),pointer     :: sup_xdm_pol0(:)      ! atom free polarizability
    logical                 :: sup_surf_loaded
    real(DEVDP),pointer     :: sup_surf_ses(:)      ! atomic radii for molecular surface
    real(DEVDP),pointer     :: sup_surf_sas(:)      ! atomic surface contribution
    real(DEVDP),pointer     :: sup_surf_atr(:)      ! atom radii
    logical                 :: sup_chrg_loaded
    logical                 :: sup_chrg_generic
    real(DEVDP),pointer     :: sup_chrg(:)          ! partial atomic charges
    logical                 :: sup_hirshfeld_loaded
    real(DEVDP),pointer     :: sup_hirshfeld(:)     ! hirshfeld atomic charges
end type GEOMETRY

! ------------------------------------------------------------------------------
! force constants for collective variables

real(DEVDP)     :: DIS_FC   = 1000.0     ! kcal/mol/A^2
real(DEVDP)     :: ANG_FC   = 1000.0     ! kcal/mol/rad^2

! which energy should be loaded
character(len=MAX_PATH) :: LoadEnergy   = ''

! which energy should be loaded for sapt
character(len=MAX_PATH) :: LoadSAPT     = ''

! which energy should be loaded for probes
character(len=MAX_PATH) :: LoadProbe    = ''

! which charges should be loaded
character(len=MAX_PATH) :: LoadCharges  = ''

! ------------------------------------------------------------------------------

integer,parameter   :: GEO_INFO_ABSENERGY = 1
integer,parameter   :: GEO_INFO_RELENERGY = 2
integer,parameter   :: GEO_INFO_PRBENERGY = 3
integer,parameter   :: GEO_INFO_NOENERGY  = 4

! ------------------------------------------------------------------------------

end module ffdev_geometry_dat

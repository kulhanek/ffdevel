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
    real(DEVDP)             :: nb14_ene         ! 1-4 vdW
    real(DEVDP)             :: ele_ene          ! electrostatics
    real(DEVDP)             :: nb_ene           ! vdW
    real(DEVDP)             :: nb_rep           ! total repulsion
    real(DEVDP)             :: nb_disp          ! total dispersion
    real(DEVDP)             :: total_ene        ! total energy from MM
! sapt0 energies
    real(DEVDP)             :: sapt0_ele         ! electrostatics
    real(DEVDP)             :: sapt0_rep         ! NB repulsion
    real(DEVDP)             :: sapt0_disp        ! NB dispersion
    real(DEVDP)             :: sapt0_total
! derivatives, etc.
    real(DEVDP),pointer     :: grd(:,:)         ! gradient (3,natoms)
    real(DEVDP),pointer     :: hess(:,:,:,:)    ! hessian (3,natoms,3,natoms) or normal modes
    real(DEVDP),pointer     :: ihess(:)         ! hessian in internal coordinates (na+nb+nd+ni)
    real(DEVDP),pointer     :: nmodes(:,:,:,:)  ! normal modes (3,natoms,3,natoms)
    real(DEVDP),pointer     :: freq(:)          ! frequencies of normal vibrations (3*natoms)
! target data
    logical                 :: trg_crd_optimized    ! indicates that the goemtry was optimized
    logical                 :: trg_ene_loaded
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
! target SAPT0
    logical                 :: trg_sapt0_loaded
    real(DEVDP)             :: trg_sapt0_ele
    real(DEVDP)             :: trg_sapt0_exch
    real(DEVDP)             :: trg_sapt0_ind
    real(DEVDP)             :: trg_sapt0_disp
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
    logical                 :: sup_chrg_loaded
    character(len=MAX_PATH) :: sup_chrg_type
    real(DEVDP),pointer     :: sup_chrg(:)          ! partial atomic charges
end type GEOMETRY

! ------------------------------------------------------------------------------
! force constants for collective variables

real(DEVDP)     :: DIS_FC   = 1000.0     ! kcal/mol/A^2
real(DEVDP)     :: ANG_FC   = 1000.0     ! kcal/mol/rad^2

! which charges should be loaded
character(len=MAX_PATH) :: LoadCharges  = 'DDEC6'

! ------------------------------------------------------------------------------

end module ffdev_geometry_dat

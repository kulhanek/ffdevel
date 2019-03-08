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

module ffdev_geometry_dat

use ffdev_sizes

! ------------------------------------------------------------------------------

type COLVAR
    real(DEVDP)                 :: trg_value
    real(DEVDP)                 :: value
    character(len=MAX_CVTYPE)   :: cvtype
    integer,pointer             :: ai(:)
end type COLVAR

! ------------------------------------------------------------------------------

type GEOMETRY
    integer                 :: id
    integer                 :: natoms           ! number of atoms
    character(len=MAX_PATH) :: name
    character(len=MAX_PATH) :: title
    integer,pointer         :: z(:)             ! proton numbers
    real(DEVDP),pointer     :: crd(:,:)         ! geometry (3,natoms)
    real(DEVDP)             :: bond_ene         ! bond energy
    real(DEVDP)             :: angle_ene        ! angle energy
    real(DEVDP)             :: dih_ene          ! dihedral energy
    real(DEVDP)             :: impropr_ene      ! improper energy
    real(DEVDP)             :: ele14_ene        ! 1-4 electrostatics
    real(DEVDP)             :: nb14_ene         ! 1-4 vdW
    real(DEVDP)             :: ele_ene          ! electrostatics
    real(DEVDP)             :: nb_ene           ! vdW
    real(DEVDP)             :: total_ene        ! total energy from MM
    real(DEVDP),pointer     :: grd(:,:)         ! gradient (3,natoms)
    real(DEVDP),pointer     :: hess(:,:,:,:)    ! hessian (3,natoms,3,natoms) or normal modes
    real(DEVDP),pointer     :: nmodes(:,:,:,:)  ! normal modes (3,natoms,3,natoms)
    real(DEVDP),pointer     :: freq(:)          ! frequencies of normal vibrations (3*natoms)
    real(DEVDP)             :: weight           ! contribution to all data
    logical                 :: trg_ene_loaded
    logical                 :: trg_crd_loaded
    logical                 :: trg_grd_loaded
    logical                 :: trg_hess_loaded
    logical                 :: trg_freq_loaded
    logical                 :: trg_esp_loaded
    logical                 :: trg_crd_optimized
    real(DEVDP)             :: trg_energy           ! target energy
    real(DEVDP),pointer     :: trg_crd(:,:)         ! target geometry (3,natoms)
    real(DEVDP),pointer     :: trg_grd(:,:)         ! target gradient (3,natoms)
    real(DEVDP),pointer     :: trg_hess(:,:,:,:)    ! target hessian (3,natoms,3,natoms)
    real(DEVDP),pointer     :: trg_nmodes(:,:,:,:)  ! target normal modes (3,natoms,3,natoms)
    real(DEVDP),pointer     :: trg_freq(:)          ! target frequencies (3*natoms)
    integer,pointer         :: freq_t2s_map(:)      ! trg_map(trg) = freq; map to freq by angle between normal vectors
    integer                 :: esp_npoints          ! number of ESP points
    real(DEVDP),pointer     :: trg_esp(:,:)         ! target ESP (4,npoints)
    integer                 :: ncvs                 ! number of collective variables
    type(COLVAR),pointer    :: cvs(:)
    real(DEVDP)             :: cvs_energy           ! colvar energy penalty
end type GEOMETRY

! ------------------------------------------------------------------------------
! force constants for collective variables

real(DEVDP)     :: DIS_FC   = 1000.0     ! kcal/mol/A^2
real(DEVDP)     :: ANG_FC   = 1000.0     ! kcal/mol/rad^2

! ------------------------------------------------------------------------------

end module ffdev_geometry_dat

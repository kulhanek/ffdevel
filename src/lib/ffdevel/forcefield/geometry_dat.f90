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

type GEOMETRY
    integer             :: id
    integer             :: natoms           ! number of atoms
    character(len=255)  :: name
    character(len=255)  :: title
    integer,pointer     :: z(:)             ! proton numbers
    real(DEVDP),pointer :: crd(:,:)         ! geometry (3,natoms)
    real(DEVDP)         :: bond_ene         ! bond energy
    real(DEVDP)         :: angle_ene        ! angle energy
    real(DEVDP)         :: dih_ene          ! dihedral energy
    real(DEVDP)         :: impropr_ene      ! improper energy
    real(DEVDP)         :: ele14_ene        ! 1-4 electrostatics
    real(DEVDP)         :: nb14_ene         ! 1-4 vdW
    real(DEVDP)         :: ele_ene          ! electrostatics
    real(DEVDP)         :: nb_ene           ! vdW
    real(DEVDP)         :: total_ene        ! total energy
    real(DEVDP),pointer :: grd(:,:)         ! gradient (3,natoms)
    real(DEVDP),pointer :: hess(:,:,:,:)    ! hessian (3,natoms,3,natoms) or normal modes
    real(DEVDP),pointer :: freq(:)          ! frequencies of normal vibrations
    real(DEVDP)         :: weight           ! contribution to all data
    logical             :: trg_ene_loaded
    logical             :: trg_grd_loaded
    logical             :: trg_hess_loaded
    logical             :: trg_esp_loaded
    real(DEVDP)         :: trg_energy           ! target energy
    real(DEVDP),pointer :: trg_grd(:,:)         ! target gradient (3,natoms)
    real(DEVDP),pointer :: trg_hess(:,:,:,:)    ! target hessian (3,natoms,3,natoms)
    integer             :: esp_npoints          ! number of ESP points
    real(DEVDP),pointer :: trg_esp(:,:)         ! target ESP (4,npoints)
end type GEOMETRY

! ------------------------------------------------------------------------------

end module ffdev_geometry_dat

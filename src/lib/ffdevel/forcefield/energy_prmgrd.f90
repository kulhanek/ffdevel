! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_energy_prmgrd

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_energy_prmgrd_all
! ==============================================================================

subroutine ffdev_energy_prmgrd_all(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    ! reset eneprmgrad
    geo%eneprmgrd(:) = 0.0d0

    ! we will simply go through all energy elemenets and

    ! bonded terms
    if( top%probe_size .eq. 0 ) then
        call ffdev_energy_prmgrd_bonds(top,geo)
        call ffdev_energy_prmgrd_angles(top,geo)
        call ffdev_energy_prmgrd_dihedrals(top,geo)
        call ffdev_energy_prmgrd_impropers(top,geo)
    end if

    ! non-bonded terms
    select case(top%nb_mode)
        case(NB_MODE_LJ)
            call ffdev_energy_prmgrd_nb_lj(top,geo)
        case(NB_MODE_BP,NB_MODE_EXP6)
            call ffdev_energy_prmgrd_nb_bp(top,geo)
        case(NB_MODE_MMD3)
            call ffdev_energy_prmgrd_nb_expd3bj(top,geo)
        case(NB_MODE_EXPONLY)
            call ffdev_energy_prmgrd_nb_exponly(top,geo)
        case default
            call ffdev_utils_exit(DEV_OUT,1,'Unsupported vdW mode in ffdev_energy_all!')
    end select


end subroutine ffdev_energy_prmgrd_all

!===============================================================================
! subroutine ffdev_energy_prmgrd_bonds
!===============================================================================

subroutine ffdev_energy_prmgrd_bonds(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    ! FIXME

end subroutine ffdev_energy_prmgrd_bonds

!===============================================================================
! subroutine ffdev_energy_prmgrd_angles
!===============================================================================

subroutine ffdev_energy_prmgrd_angles(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    ! FIXME

end subroutine ffdev_energy_prmgrd_angles

!===============================================================================
! subroutine ffdev_energy_prmgrd_dihedrals
!===============================================================================

subroutine ffdev_energy_prmgrd_dihedrals(top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_utils

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! -----------------------------------------------------------------------------

    ! FIXME

end subroutine ffdev_energy_prmgrd_dihedrals

!===============================================================================
! subroutine ffdev_energy_prmgrd_impropers
!===============================================================================

subroutine ffdev_energy_prmgrd_impropers(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! -----------------------------------------------------------------------------

    ! FIXME

end subroutine ffdev_energy_prmgrd_impropers

!===============================================================================
! subroutine ffdev_energy_prmgrd_nb_lj
!===============================================================================

subroutine ffdev_energy_prmgrd_nb_lj(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aLJa,bLJa,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a
    ! --------------------------------------------------------------------------

!    do ip=1,top%nb_size
!        i = top%nb_list(ip)%ai
!        j = top%nb_list(ip)%aj
!        nbt = top%nb_list(ip)%nbt
!        aLJa  = top%nb_types(nbt)%A
!        bLJa  = top%nb_types(nbt)%B

!        ! calculate dx, r and r2
!        dxa1 = geo%crd(1,i) - geo%crd(1,j)
!        dxa2 = geo%crd(2,i) - geo%crd(2,j)
!        dxa3 = geo%crd(3,i) - geo%crd(3,j)

!        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
!        r2a = 1.0d0/r2a
!        ra  = sqrt(r2a)
!        r6a = r2a*r2a*r2a

!        if( top%nb_list(ip)%dt .eq. 0 ) then
!            geo%nb_ene  = geo%nb_ene + aLJa*r6a*r6a - bLJa*r6a
!        else
!            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
!            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
!            geo%nb14_ene  = geo%nb14_ene + inv_scnb*(aLJa*r6a*r6a - bLJa*r6a)
!        end if
!    end do

end subroutine ffdev_energy_prmgrd_nb_lj

!===============================================================================
! subroutine ffdev_energy_prmgrd_nb_bp
!===============================================================================

subroutine ffdev_energy_prmgrd_nb_bp(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt
    real(DEVDP)     :: inv_scee,inv_scnb,aBP,bBP,cBP,crgij,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra,r6a
    ! --------------------------------------------------------------------------

!    do ip=1,top%nb_size
!        i = top%nb_list(ip)%ai
!        j = top%nb_list(ip)%aj
!        nbt = top%nb_list(ip)%nbt
!        aBP  = top%nb_types(nbt)%A
!        bBP  = top%nb_types(nbt)%B
!        cBP  = top%nb_types(nbt)%C6

!        ! calculate dx, r and r2
!        dxa1 = geo%crd(1,i) - geo%crd(1,j)
!        dxa2 = geo%crd(2,i) - geo%crd(2,j)
!        dxa3 = geo%crd(3,i) - geo%crd(3,j)

!        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
!        r2a = 1.0d0/r2a
!        ra  = sqrt(r2a)
!        r6a = r2a*r2a*r2a

!        if( top%nb_list(ip)%dt .eq. 0 ) then
!            geo%nb_ene  = geo%nb_ene + aBP*exp(-bBP/ra) - cBP*r6a
!        else
!            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
!            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
!            geo%nb14_ene  = geo%nb_ene + inv_scnb*(aBP*exp(-bBP/ra) - cBP*r6a)
!        end if
!    end do

end subroutine ffdev_energy_prmgrd_nb_bp

!===============================================================================
! subroutine ffdev_energy_nb_exponly
!===============================================================================

subroutine ffdev_energy_prmgrd_nb_exponly(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------
    integer         :: ip, i, j, nbt, pti_A, pti_B
    real(DEVDP)     :: inv_scnb,aBP,bBP,dxa1,dxa2,dxa3
    real(DEVDP)     :: r2a,ra
    ! --------------------------------------------------------------------------

    do ip=1,top%nb_size
        nbt = top%nb_list(ip)%nbt
        pti_A = top%nb_types(nbt)%pti_A
        pti_B = top%nb_types(nbt)%pti_B
        if( (pti_A .eq. 0) .and. (pti_B .eq. 0) ) cycle

        i = top%nb_list(ip)%ai
        j = top%nb_list(ip)%aj
        aBP  = top%nb_types(nbt)%A
        bBP  = top%nb_types(nbt)%B

        ! calculate dx, r and r2
        dxa1 = geo%crd(1,i) - geo%crd(1,j)
        dxa2 = geo%crd(2,i) - geo%crd(2,j)
        dxa3 = geo%crd(3,i) - geo%crd(3,j)

        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
        r2a = 1.0d0/r2a
        ra  = sqrt(r2a)

        if( top%nb_list(ip)%dt .eq. 0 ) then
            if( pti_A .ne. 0 ) then
                geo%eneprmgrd(pti_A) = geo%eneprmgrd(pti_A) + exp(-bBP/ra)
            end if
            if( pti_B .ne. 0 ) then
                geo%eneprmgrd(pti_B) = geo%eneprmgrd(pti_B) + aBP*exp(-bBP/ra)*(-1/ra)
            end if
        else
            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
            if( pti_A .ne. 0 ) then
                geo%eneprmgrd(pti_A) = geo%eneprmgrd(pti_A) + inv_scnb*exp(-bBP/ra)
            end if
            if( pti_B .ne. 0 ) then
                geo%eneprmgrd(pti_B) = geo%eneprmgrd(pti_B) + inv_scnb*aBP*exp(-bBP/ra)*(-1/ra)
            end if
        end if
    end do

end subroutine ffdev_energy_prmgrd_nb_exponly

!===============================================================================
! subroutine ffdev_energy_nb_expd3bj
!===============================================================================

subroutine ffdev_energy_prmgrd_nb_expd3bj(top,geo)

    use ffdev_topology
    use ffdev_geometry

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------


! FIXME

!    do ip=1,top%nb_size
!        i = top%nb_list(ip)%ai
!        j = top%nb_list(ip)%aj
!        nbt = top%nb_list(ip)%nbt
!        aBP  = top%nb_types(nbt)%A
!        bBP  = top%nb_types(nbt)%B
!        crgij =  top%atoms(i)%charge*top%atoms(j)%charge*332.05221729d0

!        ! calculate dx, r and r2
!        dxa1 = geo%crd(1,i) - geo%crd(1,j)
!        dxa2 = geo%crd(2,i) - geo%crd(2,j)
!        dxa3 = geo%crd(3,i) - geo%crd(3,j)

!        r2a = dxa1*dxa1 + dxa2*dxa2 + dxa3*dxa3
!        r2a = 1.0d0/r2a
!        ra  = sqrt(r2a)

!        if( top%nb_list(ip)%dt .eq. 0 ) then
!            geo%ele_ene = geo%ele_ene + crgij*ra
!            geo%nb_ene  = geo%nb_ene + aBP*exp(-bBP/ra)
!        else
!            inv_scee = top%dihedral_types(top%nb_list(ip)%dt)%inv_scee
!            inv_scnb = top%dihedral_types(top%nb_list(ip)%dt)%inv_scnb
!            geo%ele14_ene = geo%ele14_ene + inv_scee*crgij*ra
!            geo%nb14_ene  = geo%nb_ene + inv_scnb*aBP*exp(-bBP/ra)
!        end if
!    end do

end subroutine ffdev_energy_prmgrd_nb_expd3bj

! ------------------------------------------------------------------------------

end module ffdev_energy_prmgrd

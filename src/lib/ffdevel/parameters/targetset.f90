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

module ffdev_targetset

use ffdev_geometry_dat
use ffdev_constants

contains

! ==============================================================================
! subroutine ffdev_targetset_init_tops
! at this point parameters should be extracted from topologies
! ==============================================================================

subroutine ffdev_targetset_init_pts

    use ffdev_targetset_dat
    use ffdev_topology
    use ffdev_gradient_utils
    use ffdev_hessian_utils
    use ffdev_parameters_dat

    implicit none
    integer     :: i,j
    ! --------------------------------------------------------------------------

    do i=1,nsets
        call ffdev_topology_finalize_setup(sets(i)%top)
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_hess_loaded ) then
                call ffdev_hessian_allocate(sets(i)%geo(j))
                call ffdev_gradient_allocate(sets(i)%geo(j))
            else
                ! always init grd
                call ffdev_gradient_allocate(sets(i)%geo(j))
            end if
        end do
    end do

end subroutine ffdev_targetset_init_pts

! ==============================================================================
! subroutine ffdev_targetset_calc_all
! NOTE - it requires a program to be loaded
! ==============================================================================

subroutine ffdev_targetset_calc_all

    use ffdev_targetset_dat
    use ffdev_energy
    use ffdev_gradient
    use ffdev_hessian
    use ffdev_parameters_dat 
    use ffdev_parameters
    use ffdev_geoopt

    implicit none
    integer     :: i,j
    ! --------------------------------------------------------------------------

    ! apply combination rules
    if( ApplyCombinationRules ) then
        do i=1,nsets
            call ffdev_topology_apply_NB_comb_rules(sets(i)%top,sets(i)%top%assumed_comb_rules)
        end do
        call ffdev_parameters_reinit_nbparams()
    end if    
    
    ! optimize geometry
    if( OptimizeGeometry ) then
        do i=1,nsets
            do j=1,sets(i)%ngeos
                if( OptimizeOriginGeometry ) then
                    sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                end if
                if( OptimizeGeometryVerbose ) then
                    call ffdev_geoopt_run(DEV_OUT,sets(i)%top,sets(i)%geo(j))
                else
                    call ffdev_geoopt_run(DEV_NULL,sets(i)%top,sets(i)%geo(j))
                end if
                sets(i)%geo(j)%trg_crd_optimized = .true.
            end do
        end do
    else
        do i=1,nsets
            do j=1,sets(i)%ngeos
                sets(i)%geo(j)%trg_crd_optimized = .false.
            end do
        end do
    end if
    
    ! calculate all energies, gradients, hessians    
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_hess_loaded ) then
                call ffdev_hessian_all(sets(i)%top,sets(i)%geo(j))
            else if( sets(i)%geo(j)%trg_grd_loaded ) then
                call ffdev_gradient_all(sets(i)%top,sets(i)%geo(j))
            else if( sets(i)%geo(j)%trg_ene_loaded ) then
                call ffdev_energy_all(sets(i)%top,sets(i)%geo(j))
            end if
        end do
    end do

end subroutine ffdev_targetset_calc_all

! ==============================================================================
! subroutine ffdev_targetset_save_final_stops
! ==============================================================================

subroutine ffdev_targetset_save_final_stops

    use ffdev_targetset_dat
    use ffdev_topology

    implicit none
    integer     :: i
    ! --------------------------------------------------------------------------

    do i=1,nsets
        if( len_trim(sets(i)%final_name) .ne. 0 ) then
            write(DEV_OUT,10) i,trim(sets(i)%final_name)
            call ffdev_topology_save(sets(i)%top,sets(i)%final_name)
        end if
    end do

 10 format('Final topology name for set #',I2.2,' = ',A)

end subroutine ffdev_targetset_save_final_stops

! ==============================================================================
! subroutine ffdev_targetset_summary
! ==============================================================================

subroutine ffdev_targetset_summary()

    use ffdev_targetset_dat

    implicit none
    real(DEVDP)         :: toterr_ene,err_ene,toterr_grd,toterr_hess,trggrd,mmgrd,difgrd,trghess,mmhess,difhess
    integer             :: s,i,j,k,nene,ngrd,nhess,l,m
    character(len=20)   :: lname
    ! --------------------------------------------------------------------------

    write(DEV_OUT,2)

    do s=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,5) s
        write(DEV_OUT,*)
        write(DEV_OUT,10)
        write(DEV_OUT,20)

        toterr_ene = 0.0d0
        toterr_grd = 0.0d0
        toterr_hess = 0.0d0
        nene = 0
        ngrd = 0
        nhess = 0
        do i=1,sets(s)%ngeos
            err_ene = 0.0
            if( sets(s)%geo(i)%trg_ene_loaded ) then
                err_ene = sets(s)%geo(i)%total_ene - sets(s)%offset - sets(s)%geo(i)%trg_energy
                toterr_ene = toterr_ene + sets(s)%geo(i)%weight * err_ene**2
                nene = nene + 1
            end if

            trggrd = 0.0
            mmgrd = 0.0
            difgrd = 0.0
            if( sets(s)%geo(i)%trg_grd_loaded ) then
                do j=1,sets(s)%geo(i)%natoms
                    do k=1,3
                        trggrd = trggrd + sets(s)%geo(i)%trg_grd(k,j)**2
                        mmgrd = mmgrd + sets(s)%geo(i)%grd(k,j)**2
                        difgrd = difgrd + (sets(s)%geo(i)%grd(k,j) - sets(s)%geo(i)%trg_grd(k,j))**2
                    end do
                end do
                ngrd = ngrd + 1
                trggrd = sqrt(trggrd/real(3*sets(s)%geo(i)%natoms))
                mmgrd = sqrt(mmgrd/real(3*sets(s)%geo(i)%natoms))
                difgrd = sqrt(difgrd/real(3*sets(s)%geo(i)%natoms))
                toterr_grd = toterr_grd + difgrd
            end if

            trghess = 0.0
            mmhess = 0.0
            difhess = 0.0
            if( sets(s)%geo(i)%trg_hess_loaded ) then
                do j=1,sets(s)%geo(i)%natoms
                    do k=1,3
                        do l=1,sets(s)%geo(i)%natoms
                            do m=1,3
                                trghess = trghess + sets(s)%geo(i)%trg_hess(m,l,k,j)**2
                                mmhess = mmhess + sets(s)%geo(i)%hess(m,l,k,j)**2
                                difhess = difhess + (sets(s)%geo(i)%hess(m,l,k,j) - sets(s)%geo(i)%trg_hess(m,l,k,j))**2
                            end do
                        end do
                    end do
                end do
                nhess = nhess + 1
                trghess = sqrt(trghess/real(3*sets(s)%geo(i)%natoms)**2)
                mmhess = sqrt(mmhess/real(3*sets(s)%geo(i)%natoms)**2)
                difhess = sqrt(difhess/real(3*sets(s)%geo(i)%natoms)**2)
                toterr_hess = toterr_hess + difhess
            end if

            lname = trim(sets(s)%geo(i)%name)
            write(DEV_OUT,30) sets(s)%geo(i)%id,adjustl(lname),sets(s)%geo(i)%weight, &
                              sets(s)%geo(i)%total_ene-sets(s)%offset,sets(s)%geo(i)%trg_energy,err_ene, &
                              mmgrd,trggrd,difgrd,mmhess,trghess,difhess
        end do
        if( nene .gt. 0 ) then
            toterr_ene = sqrt(toterr_ene / real(nene))
        end if
        if( ngrd .gt. 0 ) then
            toterr_grd = toterr_grd / real(ngrd)
        end if
        if( nhess .gt. 0 ) then
            toterr_hess = toterr_hess / real(nhess)
        end if

        write(DEV_OUT,20)
        write(DEV_OUT,40) toterr_ene,toterr_grd,toterr_hess
    end do

 2 format('Target set summaries ...')
 5 format('=== [SET] #',I2.2,' ================================================================================&
                             &=======================================')
10 format('# ID   File                 Weight    E(MM)     E(TGR)     E(ERR)    G/c(MM)    G/c(TRG)   G/c(ERR) &
          &  H/c(MM)    H/c(TRG)   H/c(ERR) ')
20 format('# ---- -------------------- ------ ---------- ---------- ---------- ---------- ---------- ----------&
          & ---------- ---------- ----------')
30 format(I6,1X,A20,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3)
40 format(57X,F10.3,23X,F10.3,23X,F10.3)

end subroutine ffdev_targetset_summary

! ------------------------------------------------------------------------------

end module ffdev_targetset

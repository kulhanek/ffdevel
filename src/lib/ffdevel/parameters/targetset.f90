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
            else if( sets(i)%geo(j)%trg_grd_loaded .or. sets(i)%optgeo ) then
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
    use ffdev_geoopt

    implicit none
    integer     :: i,j,k
    ! --------------------------------------------------------------------------

    ! apply combination rules
    if( ApplyCombinationRules ) then
        do i=1,nsets
            call ffdev_topology_apply_NB_comb_rules(sets(i)%top,sets(i)%top%assumed_comb_rules)
        end do
        call ffdev_parameters_reinit_nbparams()
    end if    
    
    ! optimize geometry
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%optgeo .and. sets(i)%geo(j)%trg_crd_loaded ) then
                if( .not. sets(i)%keepoptgeo ) then
                    sets(i)%geo(j)%crd = sets(i)%geo(j)%trg_crd
                end if
                if( ShowOptimizationProgress ) then
                    call ffdev_geoopt_run(DEV_OUT,sets(i)%top,sets(i)%geo(j))
                else
                    call ffdev_geoopt_run(DEV_NULL,sets(i)%top,sets(i)%geo(j))
                end if
                sets(i)%geo(j)%trg_crd_optimized = .true.
            else
                sets(i)%geo(j)%trg_crd_optimized = .false.
            end if
        end do
    end do
    
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
    
   ! update energies 
    do i=1,nsets
        if( sets(i)%nrefs .gt. 0 ) then
            do j=1,sets(i)%ngeos
                do k=1,sets(i)%nrefs
                    sets(i)%geo(j)%bond_ene     = sets(i)%geo(j)%bond_ene - sets( sets(i)%refs(k) )%geo(1)%bond_ene
                    sets(i)%geo(j)%angle_ene    = sets(i)%geo(j)%angle_ene - sets( sets(i)%refs(k) )%geo(1)%angle_ene
                    sets(i)%geo(j)%dih_ene      = sets(i)%geo(j)%dih_ene - sets( sets(i)%refs(k) )%geo(1)%dih_ene
                    sets(i)%geo(j)%impropr_ene  = sets(i)%geo(j)%impropr_ene - sets( sets(i)%refs(k) )%geo(1)%impropr_ene
                    sets(i)%geo(j)%ele14_ene    = sets(i)%geo(j)%ele14_ene - sets( sets(i)%refs(k) )%geo(1)%ele14_ene
                    sets(i)%geo(j)%nb14_ene     = sets(i)%geo(j)%nb14_ene - sets( sets(i)%refs(k) )%geo(1)%nb14_ene
                    sets(i)%geo(j)%ele_ene      = sets(i)%geo(j)%ele_ene - sets( sets(i)%refs(k) )%geo(1)%ele_ene
                    sets(i)%geo(j)%nb_ene       = sets(i)%geo(j)%nb_ene - sets( sets(i)%refs(k) )%geo(1)%nb_ene
                    sets(i)%geo(j)%total_ene    = sets(i)%geo(j)%total_ene - sets( sets(i)%refs(k) )%geo(1)%total_ene  
                end do
            end do
        end if
    end do          

end subroutine ffdev_targetset_calc_all

! ==============================================================================
! subroutine ffdev_parameters_reinit_nbparams
! ==============================================================================

subroutine ffdev_parameters_reinit_nbparams()

    use ffdev_parameters_dat
    use ffdev_targetset_dat

    implicit none
    integer     :: i,j
    ! --------------------------------------------------------------------------

    ! update parameter values
    do i=1,nparams
        select case(params(i)%realm)
            case(REALM_VDW_EPS,REALM_VDW_R0,REALM_VDW_ALPHA,REALM_VDW_A,REALM_VDW_B,REALM_VDW_G, &
                 REALM_VDW_C6,REALM_VDW_C8,REALM_VDW_C10,REALM_VDW_RC)
                do j=1,nsets
                    if( params(i)%ids(j) .eq. 0 ) cycle
                    select case(params(i)%realm)
                        case(REALM_VDW_EPS)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%eps
                        case(REALM_VDW_R0)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%r0
                        case(REALM_VDW_ALPHA)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%alpha
                        case(REALM_VDW_A)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%A
                        case(REALM_VDW_B)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%B
                        case(REALM_VDW_G)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%G                           
                        case(REALM_VDW_C6)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%C6
                        case(REALM_VDW_C8)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%C8
                        case(REALM_VDW_C10)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%C10
                        case(REALM_VDW_RC)
                            params(i)%value = sets(j)%top%nb_types(params(i)%ids(j))%rc                            
                    end select
                end do
        end select
    end do

end subroutine ffdev_parameters_reinit_nbparams

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
        if( len_trim(sets(i)%final_stop) .ne. 0 ) then
            write(DEV_OUT,10) i,trim(sets(i)%final_stop)
            call ffdev_topology_save(sets(i)%top,sets(i)%final_stop)
        end if
    end do

 10 format('Final topology name for set #',I2.2,' = ',A)

end subroutine ffdev_targetset_save_final_stops

! ==============================================================================
! subroutine ffdev_targetset_summary
! ==============================================================================

subroutine ffdev_targetset_summary()

    use ffdev_targetset_dat
    use ffdev_geometry

    implicit none
    real(DEVDP)             :: toterr_ene,toterr_grd,toterr_hess,toterr_bond,toterr_angle,toterr_tors,toterr_nbs
    real(DEVDP)             :: diffene,difgrd,difhess,difbond,difangle,diftors,difnbs
    real(DEVDP)             :: d0,dt,sw,err,swsum,sswsum 
    real(DEVDP)             :: stoterr_ene,stoterr_grd,stoterr_hess,stoterr_bond,stoterr_angle,stoterr_tors,stoterr_nbs
    integer                 :: snene,sngrd,snhess,snbond,snangle,sntors,snbs
    integer                 :: i,j,k,l,m,n,nene,ngrd,nhess,ai,aj,ak,al,nbond,nangle,ntors,q,rnbds,nbs
    character(len=20)       :: lname
    character(len=MAX_PATH) :: sname
    ! --------------------------------------------------------------------------

    stoterr_ene = 0.0d0
    stoterr_grd = 0.0d0
    stoterr_hess = 0.0d0
    stoterr_bond = 0.0d0 
    stoterr_angle = 0.0d0
    stoterr_tors = 0.0d0
    stoterr_nbs = 0.0d0  
    
    snene = 0
    sngrd = 0
    snhess = 0
    snbond = 0
    snangle = 0
    sntors = 0
    snbs = 0    
    
    sswsum = 0.0
    
    ! NOTES - partial errors are not weighted
    
    write(DEV_OUT,*)       
    write(DEV_OUT,2)
    write(DEV_OUT,3)    
    
    do i=1,nsets
        write(DEV_OUT,*)
        write(DEV_OUT,5) i
        write(DEV_OUT,*)

        toterr_ene = 0.0d0
        toterr_grd = 0.0d0
        toterr_hess = 0.0d0
        toterr_bond = 0.0d0 
        toterr_angle = 0.0d0
        toterr_tors = 0.0d0
        toterr_nbs = 0.0d0
        nene = 0
        ngrd = 0
        nhess = 0
        nbond = 0
        nangle = 0
        ntors = 0
        nbs = 0
        
        do j=1,sets(i)%ngeos
        
            ! save geometry if requested
            if( sets(i)%savegeo .and. sets(i)%geo(j)%trg_crd_optimized ) then
                sname = trim(sets(i)%geo(j)%name) // 's'
                call ffdev_geometry_save_xyz(sets(i)%geo(j),sname)
            end if
        
            ! ------------------------------------------------------------------        
            diffene = 0.0
            if( sets(i)%geo(j)%trg_ene_loaded ) then
                err = sets(i)%geo(j)%total_ene - sets(i)%offset - sets(i)%geo(j)%trg_energy
                diffene = err
                toterr_ene = toterr_ene + err**2
                nene = nene + 1
                stoterr_ene = stoterr_ene + sets(i)%geo(j)%weight * err**2
                snene = snene + 1
            end if

            ! ------------------------------------------------------------------            
            difgrd = 0.0
            if( sets(i)%geo(j)%trg_grd_loaded ) then
                do k=1,sets(i)%geo(j)%natoms
                    do l=1,3
                        err = (sets(i)%geo(j)%grd(l,k) - sets(i)%geo(j)%trg_grd(l,k))**2
                        difgrd = difgrd + err**2
                        sngrd = sngrd + 1
                        stoterr_grd = stoterr_grd + sets(i)%geo(j)%weight*err**2
                    end do
                end do
                ngrd = ngrd + 1
                difgrd = sqrt(difgrd/real(3*sets(i)%geo(j)%natoms))
                toterr_grd = toterr_grd + difgrd
            end if

            ! ------------------------------------------------------------------            
            difhess = 0.0
            if( sets(i)%geo(j)%trg_hess_loaded ) then
                do k=1,sets(i)%geo(j)%natoms
                    do l=1,3
                        do m=1,sets(i)%geo(j)%natoms
                            do n=1,3
                                err = (sets(i)%geo(j)%hess(l,k,n,m) - sets(i)%geo(j)%trg_hess(l,k,n,m))**2
                                difhess = difhess + err**2
                                snhess = snhess + 1
                                stoterr_hess = stoterr_hess + sets(i)%geo(j)%weight*err**2                               
                            end do
                        end do
                    end do
                end do
                nhess = nhess + 1
                difhess = sqrt(difhess/real(3*sets(i)%geo(j)%natoms)**2)
                toterr_hess = toterr_hess + difhess
            end if
            
            ! ------------------------------------------------------------------
            difbond =  0.0
            if( sets(i)%geo(j)%trg_crd_loaded .and. sets(i)%geo(j)%trg_crd_optimized ) then
                do q=1,sets(i)%top%nbonds
                    ai = sets(i)%top%bonds(q)%ai
                    aj = sets(i)%top%bonds(q)%aj
                    d0 = ffdev_geometry_get_length(sets(i)%geo(j)%crd,ai,aj)                 
                    dt = ffdev_geometry_get_length(sets(i)%geo(j)%trg_crd,ai,aj)
                    err = d0 - dt
                    difbond = difbond + err**2
                    snbond = snbond + 1
                    stoterr_bond = stoterr_bond + sets(i)%geo(j)%weight*err**2  
                end do
                nbond = nbond + 1     
                if( sets(i)%top%nbonds .gt. 0 ) then                
                    difbond = sqrt(difbond/real(sets(i)%top%nbonds)) 
                end if
                toterr_bond = toterr_bond +  difbond              
            end if
            ! ------------------------------------------------------------------
            difangle = 0.0
            if( sets(i)%geo(j)%trg_crd_loaded .and. sets(i)%geo(j)%trg_crd_optimized ) then
                do q=1,sets(i)%top%nangles
                    ai = sets(i)%top%angles(q)%ai
                    aj = sets(i)%top%angles(q)%aj
                    ak = sets(i)%top%angles(q)%ak
                    d0 = ffdev_geometry_get_angle(sets(i)%geo(j)%crd,ai,aj,ak)                
                    dt = ffdev_geometry_get_angle(sets(i)%geo(j)%trg_crd,ai,aj,ak)
                    err = d0 - dt
                    difangle = difangle + err**2
                    snangle = snangle + 1 
                    stoterr_angle = stoterr_angle + sets(i)%geo(j)%weight*err**2
                end do
                nangle = nangle + 1  
                if( sets(i)%top%nangles .gt. 0 ) then
                    difangle = sqrt(difangle/real(sets(i)%top%nangles))
                end if
                toterr_angle = toterr_angle + difangle
            end if
            ! ------------------------------------------------------------------
            diftors = 0.0
            if( sets(i)%geo(j)%trg_crd_loaded .and. sets(i)%geo(j)%trg_crd_optimized ) then
                do q=1,sets(i)%top%ndihedrals
                    ai = sets(i)%top%dihedrals(i)%ai
                    aj = sets(i)%top%dihedrals(i)%aj
                    ak = sets(i)%top%dihedrals(i)%ak
                    al = sets(i)%top%dihedrals(i)%al
                    d0 = ffdev_geometry_get_dihedral(sets(i)%geo(j)%crd,ai,aj,ak,al)
                    dt = ffdev_geometry_get_dihedral(sets(i)%geo(j)%trg_crd,ai,aj,ak,al)
                    err = ffdev_geometry_get_dihedral_deviation(d0,dt)
                    diftors = diftors + err**2
                    sntors = sntors + 1
                    stoterr_tors = stoterr_tors + sets(i)%geo(j)%weight*err**2
                end do
                ntors = ntors + 1
                if( sets(i)%top%ndihedrals .gt. 0 ) then                
                    diftors = sqrt(diftors/real(sets(i)%top%ndihedrals))
                end if
                toterr_tors = toterr_tors + diftors
            end if
            ! ------------------------------------------------------------------
            difnbs = 0.0
            swsum = 0.0
            if( sets(i)%geo(j)%trg_crd_loaded .and. sets(i)%geo(j)%trg_crd_optimized .and. (sets(i)%top%nfragments .gt. 1) ) then
                do q=1,sets(i)%top%nb_size
                    ai = sets(i)%top%nb_list(q)%ai
                    aj = sets(i)%top%nb_list(q)%aj

                    if( sets(i)%top%atoms(ai)%frgid .eq. sets(i)%top%atoms(aj)%frgid ) cycle
                    
                    d0 = ffdev_geometry_get_length(sets(i)%geo(j)%crd,ai,aj)                 
                    dt = ffdev_geometry_get_length(sets(i)%geo(j)%trg_crd,ai,aj)
                    err = d0 - dt
                    
                    ! calculate switch function
                    sw = 1.0d0 / (1.0d0 + exp( NBDistanceSWAlpha*(dt - NBDistanceSWPosition) ) )
                    swsum = swsum + sw
                    err = err * sw
                    difnbs = difnbs + sets(i)%geo(j)%weight*err**2
                    
                    sswsum = sswsum + sw                    
                    snbs = snbs + 1
                    stoterr_nbs = stoterr_nbs + sets(i)%geo(j)%weight*err**2
                end do
                nbs = nbs + 1
                if( swsum .gt. 0 ) then
                    difnbs = sqrt(difnbs/swsum)
                    toterr_nbs = toterr_nbs + difnbs
                end if
                
            end if             
            ! ------------------------------------------------------------------
            
            if( j .eq. 1 ) then
                write(DEV_OUT,10,advance='NO')
                if( nene .gt. 0 ) then
                    write(DEV_OUT,11,advance='NO')
                end if
                if( ngrd .gt. 0 ) then
                    write(DEV_OUT,12,advance='NO')
                end if
                if( nhess .gt. 0 ) then
                    write(DEV_OUT,13,advance='NO')
                end if 
                if( nbond .gt. 0 ) then
                    write(DEV_OUT,14,advance='NO')
                end if
                if( nangle .gt. 0 ) then
                    write(DEV_OUT,15,advance='NO')
                end if  
                if( ntors .gt. 0 ) then
                    write(DEV_OUT,16,advance='NO')
                end if 
                if( nbs .gt. 0 ) then
                    write(DEV_OUT,17,advance='NO')
                end if                  
                write(DEV_OUT,*)
                
                write(DEV_OUT,20,advance='NO')
                if( nene .gt. 0 ) then
                    write(DEV_OUT,21,advance='NO')
                end if
                if( ngrd .gt. 0 ) then
                    write(DEV_OUT,22,advance='NO')
                end if
                if( ngrd .gt. 0 ) then
                    write(DEV_OUT,23,advance='NO')
                end if 
                if( nbond .gt. 0 ) then
                    write(DEV_OUT,24,advance='NO')
                end if
                if( nangle .gt. 0 ) then
                    write(DEV_OUT,25,advance='NO')
                end if  
                if( ntors .gt. 0 ) then
                    write(DEV_OUT,26,advance='NO')
                end if
                if( nbs .gt. 0 ) then
                    write(DEV_OUT,26,advance='NO')
                end if                
                write(DEV_OUT,*)
             end if

            lname = trim(sets(i)%geo(j)%name)
            write(DEV_OUT,30,advance='NO') sets(i)%geo(j)%id,adjustl(lname),sets(i)%geo(j)%weight
            if( nene .gt. 0 ) then
                write(DEV_OUT,31,advance='NO') sets(i)%geo(j)%total_ene-sets(i)%offset,sets(i)%geo(j)%trg_energy,diffene
            end if
            if( ngrd .gt. 0) then
                write(DEV_OUT,32,advance='NO') difgrd
            end if
            if( nhess .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') difhess
            end if 
            if( nbond .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') difbond
            end if 
            if( nangle .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') difangle * DEV_R2D
            end if 
            if( ntors .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') diftors * DEV_R2D
            end if  
            if( nbs .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') difnbs
            end if              
            write(DEV_OUT,*)

        end do
        
        if( sets(i)%ngeos .gt. 0 ) then
            write(DEV_OUT,20,advance='NO')
            if( nene .gt. 0 ) then
                write(DEV_OUT,21,advance='NO')
            end if
            if( ngrd .gt. 0 ) then
                write(DEV_OUT,22,advance='NO')
            end if
            if( ngrd .gt. 0 ) then
                write(DEV_OUT,23,advance='NO')
            end if 
            if( nbond .gt. 0 ) then
                write(DEV_OUT,24,advance='NO')
            end if
            if( nangle .gt. 0 ) then
                write(DEV_OUT,25,advance='NO')
            end if  
            if( ntors .gt. 0 ) then
                write(DEV_OUT,26,advance='NO')
            end if
            if( nbs .gt. 0 ) then
                write(DEV_OUT,26,advance='NO')
            end if            
            write(DEV_OUT,*)     
            
            if( nene .gt. 0 ) then
                toterr_ene = sqrt(toterr_ene / real(nene))
            end if
            if( ngrd .gt. 0 ) then
                toterr_grd = toterr_grd / real(ngrd)
            end if
            if( nhess .gt. 0 ) then
                toterr_hess = toterr_hess / real(nhess)
            end if
            if( nbond .gt. 0 ) then
                toterr_bond = toterr_bond / real(nbond)
            end if 
            if( nangle .gt. 0 ) then
                toterr_angle = toterr_angle / real(nangle)
            end if    
            if( ntors .gt. 0 ) then
                toterr_tors = toterr_tors / real(ntors)
            end if
            if( nbs .gt. 0 ) then
                toterr_nbs = toterr_nbs / real(nbs)
            end if             
            
            write(DEV_OUT,40,advance='NO')
            if( nene .gt. 0 ) then
                write(DEV_OUT,41,advance='NO') toterr_ene
            end if
            if( ngrd .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') toterr_grd
            end if
            if( ngrd .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') toterr_hess
            end if 
            if( nbond .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') toterr_bond
            end if
            if( nangle .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') toterr_angle * DEV_R2D
            end if  
            if( ntors .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') toterr_tors * DEV_R2D
            end if
            if( nbs .gt. 0 ) then
                write(DEV_OUT,32,advance='NO') toterr_nbs
            end if            
            write(DEV_OUT,*)             

        end if
    end do
    
    write(DEV_OUT,*) 

! final summary

    write(DEV_OUT,6,advance='NO')
    if( nene .gt. 0 ) then
        write(DEV_OUT,9,advance='NO')
    end if
    if( ngrd .gt. 0 ) then
        write(DEV_OUT,12,advance='NO')
    end if
    if( nhess .gt. 0 ) then
        write(DEV_OUT,13,advance='NO')
    end if 
    if( nbond .gt. 0 ) then
        write(DEV_OUT,14,advance='NO')
    end if
    if( nangle .gt. 0 ) then
        write(DEV_OUT,15,advance='NO')
    end if  
    if( ntors .gt. 0 ) then
        write(DEV_OUT,16,advance='NO')
    end if 
    if( nbs .gt. 0 ) then
        write(DEV_OUT,17,advance='NO')
    end if                  
    write(DEV_OUT,*)
                
    write(DEV_OUT,20,advance='NO')
    if( nene .gt. 0 ) then
        write(DEV_OUT,21,advance='NO')
    end if
    if( ngrd .gt. 0 ) then
        write(DEV_OUT,22,advance='NO')
    end if
    if( ngrd .gt. 0 ) then
        write(DEV_OUT,23,advance='NO')
    end if 
    if( nbond .gt. 0 ) then
        write(DEV_OUT,24,advance='NO')
    end if
    if( nangle .gt. 0 ) then
        write(DEV_OUT,25,advance='NO')
    end if  
    if( ntors .gt. 0 ) then
        write(DEV_OUT,26,advance='NO')
    end if
    if( nbs .gt. 0 ) then
        write(DEV_OUT,26,advance='NO')
    end if            
    write(DEV_OUT,*)     
    
    if( snene .gt. 0 ) then
        stoterr_ene = sqrt(stoterr_ene / real(snene))
    end if
    if( sngrd .gt. 0 ) then
        stoterr_grd = sqrt(stoterr_grd / real(sngrd))
    end if
    if( snhess .gt. 0 ) then
        stoterr_hess = sqrt(stoterr_hess / real(snhess))
    end if
    if( snbond .gt. 0 ) then
        stoterr_bond = sqrt(stoterr_bond / real(snbond))
    end if 
    if( snangle .gt. 0 ) then
        stoterr_angle = sqrt(stoterr_angle / real(snangle))
    end if    
    if( sntors .gt. 0 ) then
        stoterr_tors = sqrt(stoterr_tors / real(sntors))
    end if
    if( snbs .gt. 0 ) then
        stoterr_nbs = sqrt(stoterr_nbs / sswsum)
    end if             
    
    write(DEV_OUT,40,advance='NO')
    if( nene .gt. 0 ) then
        write(DEV_OUT,41,advance='NO') stoterr_ene
    end if
    if( ngrd .gt. 0 ) then
        write(DEV_OUT,32,advance='NO') stoterr_grd
    end if
    if( ngrd .gt. 0 ) then
        write(DEV_OUT,32,advance='NO') stoterr_hess
    end if 
    if( nbond .gt. 0 ) then
        write(DEV_OUT,32,advance='NO') stoterr_bond
    end if
    if( nangle .gt. 0 ) then
        write(DEV_OUT,32,advance='NO') stoterr_angle * DEV_R2D
    end if  
    if( ntors .gt. 0 ) then
        write(DEV_OUT,32,advance='NO') stoterr_tors * DEV_R2D
    end if
    if( nbs .gt. 0 ) then
        write(DEV_OUT,32,advance='NO') stoterr_nbs
    end if            
    write(DEV_OUT,*)
    
    write(DEV_OUT,*)    
    
 2 format('Target set summaries ...')
 3 format('NOTE: Errors reported for individual points are not weighted.')
 
 5 format('=== [SET] #',I2.2,' ==================================================================')

10 format('# ID   File                 Weight ')
 6 format('#                                  ')
20 format('# ---- -------------------- ------ ')

11 format('   E(MM)     E(TGR)     E(Err)   ')
 9 format('                        E(Err)   ')
21 format('---------- ---------- ---------- ')

12 format(' G/c(Err)  ')
22 format('---------- ')

13 format(' H/c(Err)  ')
23 format('---------- ')

14 format(' Bond [A]  ')
24 format('---------- ')

15 format('Angle[deg] ')
25 format('---------- ')

16 format('Tors [deg] ')
26 format('---------- ')

17 format('d(NBs) [A] ')
27 format('---------- ')


30 format(I6,1X,A20,1X,F6.3,1X)
31 format(F10.3,1X,F10.3,1X,F10.3,1X)
32 format(F10.3,1X)

40 format(35X)
41 format(22X,F10.3,1X)

end subroutine ffdev_targetset_summary

! ------------------------------------------------------------------------------

end module ffdev_targetset

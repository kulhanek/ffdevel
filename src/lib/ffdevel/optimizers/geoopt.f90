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

module ffdev_geoopt

use ffdev_sizes
use ffdev_constants
use ffdev_variables

implicit none
contains

!===============================================================================
! subroutine ffdev_geoopt_set_default
!===============================================================================

subroutine ffdev_geoopt_set_default()

    use ffdev_geoopt_dat

    implicit none
    ! --------------------------------------------------------------------------

! === [minimize] ===============================================================
    OptimizationMethod  = MINIMIZATION_LBFGS
    NOptSteps           =  5000      ! max number of steps
    OutSamples          =    20      ! how often write results
    TrajSamples         =   100

! maximum number of steps is nsteps - this is becuase of change of restraints etc
    MaxRMSG             = 0.001d0
    MaxG                = 0.001d0
    MinEnergyChange     = -1        ! if negative number - this test is switched off
    PrintFinalGradient  = .false.
    PrintRSTSummary     = .true.

! === [steepest-descent] =======================================================
    InitialStepSize     = 0.001
    MaximalStepSize     = 0.010
    AcceptRatio         = 1.2000
    RejectRatio         = 0.5000
    AdaptiveStep        = .true.

! === [L-BFGS] =================================================================
    NumberOfCorrections = 50

end subroutine ffdev_geoopt_set_default

!===============================================================================
! subroutine ffdev_geoopt_reset_stat_counters
!===============================================================================

subroutine ffdev_geoopt_reset_stat_counters()

    use ffdev_geoopt_dat

    implicit none
    ! --------------------------------------------------------------------------

    NumberOfRuns        = 0
    NumberOfFailedRuns  = 0
    NumberOfGrdEvals    = 0

end subroutine ffdev_geoopt_reset_stat_counters

!===============================================================================
! subroutine ffdev_geoopt_run
!===============================================================================

subroutine ffdev_geoopt_run(fout,top,geo)

    use ffdev_geoopt_dat
    use ffdev_topology
    use ffdev_geometry
    use ffdev_timers

    implicit none
    integer         :: fout
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    ! --------------------------------------------------------------------------

    call ffdev_timers_start_timer(FFDEV_GEOOPT_TIMER)

    if( PrintRSTSummary ) then
        ! print summary about RSTs
        write(fout,*)
        call ffdev_geometry_rstsum(fout,geo)
        write(fout,*)
    end if

    call write_header(fout)

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT)
            call opt_steepest_descent(fout,top,geo)
        case(MINIMIZATION_LBFGS)
            call opt_lbfgs(fout,top,geo)
    end select

    call ffdev_timers_stop_timer(FFDEV_GEOOPT_TIMER)

end subroutine ffdev_geoopt_run

!===============================================================================
! subroutine ffdev_geoopt_opentraj
!===============================================================================

subroutine ffdev_geoopt_opentraj(top)

    use ffdev_geoopt_dat
    use ffdev_topology
    use smf_periodic_table_dat
    use smf_xyzfile

    implicit none
    type(TOPOLOGY)  :: top
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    ! load xyz file
    call init_xyz(OptTrajFile)
    call allocate_xyz(OptTrajFile,top%natoms)

    ! copy data
    OptTrajFile%comment  = 'ffdevel traj xyz file'

    do i=1,top%natoms
        OptTrajFile%symbols(i) = pt_symbols(top%atom_types(top%atoms(i)%typeid)%z)
    end do

    ! write data
    call open_xyz(DEV_TRAJ,OptTrajName,OptTrajFile,'UNKNOWN')

end subroutine ffdev_geoopt_opentraj

!===============================================================================
! subroutine ffdev_geoopt_closetraj
!===============================================================================

subroutine ffdev_geoopt_closetraj()

    use ffdev_geoopt_dat
    use ffdev_topology
    use smf_xyzfile

    implicit none
    ! --------------------------------------------------------------------------

    call close_xyz(DEV_TRAJ,OptTrajFile)
    call free_xyz(OptTrajFile)

end subroutine ffdev_geoopt_closetraj

!===============================================================================
! subroutine opt_steepest_descent
!===============================================================================

subroutine opt_steepest_descent(fout,top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geoopt_dat
    use ffdev_gradient
    use ffdev_gradient_utils
    use ffdev_utils

    implicit none
    integer                 :: fout
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    ! --------------------------------------------
    integer                 :: istep,maxatom,alloc_status
    real(DEVDP)             :: rmsg, maxgrad, lastenergy, stepsize, totene
    real(DEVDP),allocatable :: tmp_xg1(:,:),tmp_xg2(:,:)
    ! --------------------------------------------------------------------------

    lastenergy  = 0.0d0
    stepsize    = InitialStepSize

    ! allocate working array
    allocate(tmp_xg1(3,top%natoms), tmp_xg2(3,top%natoms), &
             stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data for SD optimization!')
    end if

    tmp_xg1(:,:) = geo%crd(:,:)

    do istep = 1, NOptSteps

        !===============================================================================
        ! get potential energy and derivatives from FF
        call ffdev_gradient_all(top,geo)

        ! rst penalty
        call ffdev_geometry_get_rst_penalty(geo)

        totene = geo%total_ene + geo%rst_energy

        !===============================================================================
        ! check all criteria
        rmsg = ffdev_gradient_rmsg(geo,maxgrad,maxatom)

        if( istep .ne. 1 .and. abs(totene - lastenergy) .le. MinEnergyChange ) then
            write(fout,'(/,a,F16.4)') ' >>> INFO: Last energy change     : ', abs(totene - lastenergy)
            write(fout,'(a,F16.4)')   ' >>> INFO: Energy change treshold : ', MinEnergyChange
            write(fout,'(a,/)') ' >>> INFO: Energy change is below treshold! Minimization was stoped.'
            exit
        end if

        if( abs(maxgrad) .le. MaxG .and. rmsg .le. MaxRMSG ) then
            write(fout,'(/,a,F16.4)') ' >>> INFO: RMS of gradient                 : ', rmsg
            write(fout,'(a,F16.4)')   ' >>> INFO: RMS of gradient treshold        : ', MaxRMSG
            write(fout,'(a,F16.4)')   ' >>> INFO: Max gradient component          : ', abs(maxgrad)
            write(fout,'(a,F16.4)')   ' >>> INFO: Max gradient component treshold : ', MaxG
            write(fout,'(a,F16.4)')   ' >>> INFO: Last energy change              : ', abs(totene - lastenergy)
            write(fout,'(a,/)') ' >>> INFO: Gradient tresholds were satisfied! Minimization was stoped.'
            exit
        end if

        !===============================================================================
        ! print [intermediate] results (master node only)
        call write_results(fout,istep,geo,rmsg,maxgrad,maxatom,.false.)

        ! if this is last step do not update coordinates and exit cycle
        if( istep .eq. NOptSteps ) exit

        !===============================================================================
        ! correct step size and do steepest-descent minimization

        if( AdaptiveStep .and. istep .ne. 1 .and. totene .lt. lastenergy ) then
            stepsize = stepsize * AcceptRatio
            if( stepsize .gt. MaximalStepSize ) then
                stepsize = MaximalStepSize
            end if
            tmp_xg1(:,:)    = geo%crd(:,:)
            tmp_xg2(:,:)    = geo%grd(:,:)
            lastenergy      = totene
            geo%crd(:,:)       = geo%crd(:,:) - geo%grd(:,:)*stepsize/sqrt(rmsg)
            !write(DEV_OUT,'(/,a,I10,a)') '>>> INFO: Minimization step ',istep,' was accepted!'
        else if ( adaptivestep .and. istep .ne. 1 .and. totene .ge. lastenergy ) then
            ! go back and try smaller step
            stepsize = stepsize * RejectRatio
            geo%crd(:,:)       = tmp_xg1(:,:)
            geo%grd(:,:)       = tmp_xg2(:,:)
            rmsg = ffdev_gradient_rmsg(geo,maxgrad,maxatom)
            geo%crd(:,:)       = geo%crd(:,:) - geo%grd(:,:)*stepsize/sqrt(rmsg)
            !write(DEV_OUT,'(/,a,I10,a)') '>>> INFO: Minimization step ',istep,' was rejected!'
        else
            ! first step
            tmp_xg1(:,:)    = geo%crd(:,:)
            tmp_xg2(:,:)    = geo%grd(:,:)
            lastenergy      = totene
            geo%crd(:,:)    = geo%crd(:,:) - geo%grd(:,:)*stepsize/sqrt(rmsg)
        end if

    end do

    !===============================================================================
    if( fout .gt. 0 ) then
        if( istep .le. NOptSteps ) then
            write(fout,*)
            call ffdev_utils_heading(fout,'Final results', '-')
        else
            write(fout,'(/,a)') ' >>> INFO: Maximum number of minimization steps was reached!'
            write(fout,'(a,/)') ' >>> WARNING: Minimization was not completed!'
            call ffdev_utils_heading(fout,'Intermediate results', '-')
        end if
        call write_header(fout)
        ! write final results
        call write_results(fout,istep,geo,rmsg,maxgrad,maxatom,.true.) ! results to stdout

        if( PrintRSTSummary ) then
            ! print summary about RSTs
            write(fout,*)
            call ffdev_geometry_rstsum(fout,geo)
        end if

        if( PrintFinalGradient ) then
            write(fout,*)
            call ffdev_utils_heading(fout,'Final gradient', '-')
            call ffdev_gradient_print(fout,top,geo)
        end if
    end if

end subroutine opt_steepest_descent

!===============================================================================
! subroutine opt_lbfgs
!===============================================================================

subroutine opt_lbfgs(fout,top,geo)

    use ffdev_topology
    use ffdev_geometry
    use ffdev_geoopt_dat
    use ffdev_gradient
    use ffdev_gradient_utils
    use ffdev_utils
    use lbfgsmodule

    implicit none
    integer                 :: fout
    type(TOPOLOGY)          :: top
    type(GEOMETRY)          :: geo
    ! --------------------------------------------
    integer                 :: istep,maxatom,alloc_status
    real(DEVDP)             :: rmsg, maxgrad, lastenergy, eps, xtol, totene
    integer                 :: iprint(2),iflag
    real(DEVDP),allocatable :: work(:)
    real(DEVDP),allocatable :: tmp_xg(:,:)
    type(LBFGSCTX)          :: ctx
    ! --------------------------------------------------------------------------

    ! init required variables ====================
    lastenergy     = 0.0d0
    iflag          = 0
    iprint(1)      = -1
    iprint(2)      = 1
    eps            = 1.0d-16
    xtol           = 1.0d-16   ! this is unrealistic criteria - we use xbp_own based on gradient

    ! allocate working array
    allocate(work(3*top%natoms*(2*NumberOfCorrections+1)+2*NumberOfCorrections), &
             tmp_xg(3,top%natoms), &
             stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data for L-BFGS optimization!')
    end if

    ! perform minimization ========================
    do istep = 1, NOptSteps

        !===============================================================================
        ! get potential energy and derivatives from FF
        call ffdev_gradient_all(top,geo)

        ! rst penalty
        call ffdev_geometry_get_rst_penalty(geo)

        totene = geo%total_ene + geo%rst_energy

        !===============================================================================
        ! check all criteria
        rmsg = ffdev_gradient_rmsg(geo,maxgrad,maxatom)

        if( istep .ne. 1 .and. abs(totene - lastenergy) .le. MinEnergyChange ) then
            if( fout .gt. 0 ) then
                write(fout,'(/,a,/)') ' >>> INFO: Energy change is below treshold! Minimization was stoped.'
                write(fout,'(a,F16.4)') ' >>> INFO: Last energy change     : ', abs(totene - lastenergy)
                write(fout,'(a,F16.4)') ' >>> INFO: Energy change treshold : ', MinEnergyChange
            end if
            exit
        end if

        if( abs(maxgrad) .le. MaxG .and. rmsg .le. MaxRMSG ) then
            write(fout,'(/,a,/)') ' >>> INFO: Gradient tresholds were satisfied! Minimization was stoped.'
            write(fout,'(a,F16.4)') ' >>> INFO: RMS of gradient                 : ', rmsg
            write(fout,'(a,F16.4)') ' >>> INFO: RMS of gradient treshold        : ', MaxRMSG
            write(fout,'(a,F16.4)') ' >>> INFO: Max gradient component          : ', abs(maxgrad)
            write(fout,'(a,F16.4)') ' >>> INFO: Max gradient component treshold : ', MaxG
            write(fout,'(a,F16.4)') ' >>> INFO: Last energy change              : ', abs(totene - lastenergy)
            exit
        end if

        !===============================================================================
        ! print [intermediate] results (master node only)
        call write_results(fout,istep,geo,rmsg,maxgrad,maxatom,.false.)

        ! if this is last step do not update coordinates and exit cycle
        if( istep .eq. NOptSteps ) exit

        !===============================================================================
        ! do L-BFGS minimization
        call LBFGS( top%natoms*3,NumberOfCorrections, &
                    geo%crd,totene,geo%grd,&
                    .false.,tmp_xg,iprint,eps,xtol,work,iflag,ctx)

        if( iflag .eq. 0 ) exit
        if( iflag .le. 0 ) then
            write(fout,'(/,a,i2,/)') '>>> ERROR: Internal L-BFGS driver error! Code = ', iflag
            exit
        end if

        lastenergy = totene
    end do

    !===============================================================================
    if( fout .gt. 0 ) then
        if( istep .le. NOptSteps ) then
            write(fout,*)
            call ffdev_utils_heading(fout,'Final results', '-')
        else
            write(fout,'(/,a)') ' >>> INFO: Maximum number of minimization steps was reached!'
            write(fout,'(a,/)') ' >>> WARNING: Minimization was not completed!'
            call ffdev_utils_heading(fout,'Intermediate results', '-')
        end if

        call write_header(fout)
        ! write final results
        call write_results(fout,istep,geo,rmsg,maxgrad,maxatom,.true.) ! results to stdout

        if( PrintRSTSummary ) then
            ! print summary about RSTs
            write(fout,*)
            call ffdev_geometry_rstsum(fout,geo)
        end if

        if( PrintFinalGradient ) then
            write(fout,*)
            call ffdev_utils_heading(fout,'Final gradient', '-')
            call ffdev_gradient_print(fout,top,geo)
        end if
    end if

    deallocate(work,tmp_xg)

end subroutine opt_lbfgs

!===============================================================================
! subroutine ffdev_geoopt_header
!===============================================================================

subroutine write_header(fout)

    use ffdev_geoopt_dat
    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer         :: fout
    ! --------------------------------------------------------------------------

    if( fout .le. 0 ) return

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT)
            write(fout,10)
        case(MINIMIZATION_LBFGS)
            write(fout,15)
    end select

    write(fout,20)
    write(fout,30)

 10 format('# Mode = Steepest Descent')
 15 format('# Mode = L-BFGS')
 20 format('# STEP     Etot      Ebonded      Enb        Ecvs       RMSG       maxG    maxAt')
 30 format('#----- ------------ ---------- ---------- ---------- ---------- ---------- -----')
!30 format('#-------------------------------------------------------------------------------'

end subroutine write_header

!===============================================================================
! subroutine write_results
!===============================================================================

subroutine write_results(fout,istep,geo,rmsg,maxgrad,maxatom,done)

    use ffdev_geometry
    use ffdev_geoopt_dat
    use smf_xyzfile

    implicit none
    integer         :: fout
    integer         :: istep
    type(GEOMETRY)  :: geo
    real(DEVDP)     :: rmsg
    real(DEVDP)     :: maxgrad
    integer         :: maxatom
    logical         :: done
    ! --------------------------------------------
    real(DEVDP)     :: Ebn, Enb
    ! -----------------------------------------------------------------------------

    ! write energies
    if( done .or. ((OutSamples .gt. 0) .and. (mod(istep,OutSamples) .eq. 0)) .or. (istep .eq. 1) ) then
        Ebn = geo%bond_ene +geo%angle_ene+geo%dih_ene+geo%impropr_ene
        Enb = geo%ele14_ene+geo%rep14_ene+geo%dis14_ene + &
              geo%ele_ene+geo%rep_ene+geo%dis_ene
        write(fout,10) istep, geo%total_ene+geo%rst_energy,Ebn, &
                       Enb,geo%rst_energy,rmsg,maxgrad,maxatom
    end if

    ! write trajectory
    if( OptTrajFile%opened ) then
        if( (TrajSamples .gt. 0) .and. (mod(istep,TrajSamples) .eq. 0) ) then
            OptTrajFile%cvs = geo%crd
            call write_xyz(DEV_TRAJ,OptTrajFile)
        end if
    end if

 10 format(I6,1X,F12.6,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4,1X,I5)

end subroutine write_results

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_geoopt





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

module ffdev_ffopt

use ffdev_sizes
use ffdev_constants

implicit none
contains

!===============================================================================
! subroutine ffdev_ffopt_set_default
!===============================================================================

subroutine ffdev_ffopt_set_default()

    use ffdev_ffopt_dat

    implicit none
    include 'nlopt.f'
    ! --------------------------------------------------------------------------

! === [minimization] ===========================================================
    OptimizationMethod  = MINIMIZATION_NLOPT
    NOptSteps    = 10000      ! max number of steps
    OutSamples   =    20      ! how often write results

! maximum number of steps is nsteps - this is becuase of change of restraints etc
    MaxRMSG             = 0.00001d0
    MaxG                = 0.00001d0
    MinErrorChange      = 0.00001d0        ! negative number - this test is switched off by default
    PrintFinalGradient  = .false.

! === [steepest-descent] =======================================================
    InitialStepSize     = 0.001
    MaximalStepSize     = 0.010
    AcceptRatio         = 1.2000
    RejectRatio         = 0.5000
    AdaptiveStep        = .true.

! === [L-BFGS] =================================================================
    NumberOfCorrections = 20

! === [NLOPT] ==================================================================
    NLOpt_Method        = NLOPT_LN_COBYLA
    NLOpt_InitialStep   = 0.00001d0

end subroutine ffdev_ffopt_set_default

!===============================================================================
! subroutine ffdev_ffopt_run
!===============================================================================

subroutine ffdev_ffopt_run()

    use ffdev_ffopt_dat
    use ffdev_parameters_dat
    use ffdev_parameters
    use ffdev_utils
    use ffdev_targetset
    use ffdev_targetset_dat

    implicit none
    integer     :: alloc_stat, i
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'FF Parameter Optimalization', ':')

    allocate(FFParams(nactparms), FFParamsGrd(nactparms), stat=alloc_stat)
    if( alloc_stat .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate data for FF optimization!')
    end if

    ! get initial parameters
    call ffdev_parameters_gather(FFParams)

    ! initial statistics
    call ffdev_parameters_scatter(FFParams)
    call ffdev_parameters_to_tops()
    call ffdev_targetset_calc_all()
    call ffdev_targetset_summary()

    call write_header()

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT)
            call opt_steepest_descent
        case(MINIMIZATION_LBFGS)
            call opt_lbfgs
        case(MINIMIZATION_NLOPT)
            call opt_nlopt
    end select

    ! return finial statistics
    call ffdev_parameters_scatter(FFParams)
    call ffdev_parameters_to_tops()
    call ffdev_targetset_calc_all()
    call ffdev_targetset_summary()

    deallocate(FFParams,FFParamsGrd)

end subroutine ffdev_ffopt_run

!===============================================================================
! subroutine opt_steepest_descent
!===============================================================================

subroutine opt_steepest_descent()

    use ffdev_ffopt_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_parameters

    implicit none
    integer                 :: istep,alloc_status
    real(DEVDP)             :: rmsg, maxgrad, lasterror, stepsize
    real(DEVDP),allocatable :: tmp_xg1(:),tmp_xg2(:)
    ! --------------------------------------------------------------------------

    lasterror  = 0.0d0
    stepsize    = InitialStepSize

    ! allocate working array
    allocate(tmp_xg1(nactparms), tmp_xg2(nactparms), &
             stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate data for SD optimization!')
    end if

    tmp_xg1(:) = FFParams(:)

    do istep = 1, NOptSteps

        !===============================================================================
        ! get fit error
        call ffdev_parameters_error(FFParams,FFError,FFParamsGrd)

        !===============================================================================
        ! check all criteria
        rmsg = ffdev_fopt_rmsg(FFParamsGrd,maxgrad)

        if( istep .ne. 1 .and. abs(FFError%total - lasterror) .le. MinErrorChange ) then
            write(DEV_OUT,'(/,a,E16.10)') ' >>> INFO: Last error change      : ', abs(FFError%total - lasterror)
            write(DEV_OUT,'(a,E16.10)')   ' >>> INFO: Error change treshold : ', MinErrorChange
            write(DEV_OUT,'(a,/)') ' >>> INFO: Error change is below treshold! Minimization was stoped.'
            exit
        end if

        if( abs(maxgrad) .le. MaxG .and. rmsg .le. MaxRMSG ) then
            write(DEV_OUT,'(/,a,E16.10)') ' >>> INFO: RMS of gradient                 : ', rmsg
            write(DEV_OUT,'(a,E16.10)')   ' >>> INFO: RMS of gradient treshold        : ', MaxRMSG
            write(DEV_OUT,'(a,E16.10)')   ' >>> INFO: Max gradient component          : ', abs(maxgrad)
            write(DEV_OUT,'(a,E16.10)')   ' >>> INFO: Max gradient component treshold : ', MaxG
            write(DEV_OUT,'(a,E16.10)')   ' >>> INFO: Last error change               : ', abs(FFError%total - lasterror)
            write(DEV_OUT,'(a,/)') ' >>> INFO: Gradient tresholds were satisfied! Minimization was stoped.'
            exit
        end if

        !===============================================================================
        ! print [intermediate] results (master node only)
        call write_results(istep,FFError,rmsg,maxgrad,.false.)

        ! if this is last step do not update coordinates and exit cycle
        if( istep .eq. NOptSteps ) exit

        !===============================================================================
        ! correct step size and do steepest-descent minimization

        if( AdaptiveStep .and. istep .ne. 1 .and. FFError%total .lt. lasterror ) then
            stepsize = stepsize * AcceptRatio
            if( stepsize .gt. MaximalStepSize ) then
                stepsize = MaximalStepSize
            end if
            tmp_xg1(:)      = FFParams(:)
            tmp_xg2(:)      = FFParamsGrd(:)
            lasterror       = FFError%total
            FFParams(:)     = FFParams(:) - FFParamsGrd(:)*stepsize/sqrt(rmsg)
            !write(DEV_OUT,'(/,a,I10,a)') '>>> INFO: Minimization step ',istep,' was accepted!'
        else if ( adaptivestep .and. istep .ne. 1 .and. FFError%total .ge. lasterror ) then
            ! go back and try smaller step
            stepsize = stepsize * RejectRatio
            FFParams(:)       = tmp_xg1(:)
            FFParamsGrd(:)    = tmp_xg2(:)
            rmsg = ffdev_fopt_rmsg(FFParamsGrd,maxgrad)
            FFParams(:)       = FFParams(:) - FFParamsGrd(:)*stepsize/sqrt(rmsg)
            !write(DEV_OUT,'(/,a,I10,a)') '>>> INFO: Minimization step ',istep,' was rejected!'
        else
            ! first step
            tmp_xg1(:)      = FFParams(:)
            tmp_xg2(:)      = FFParamsGrd(:)
            lasterror       = FFError%total
            FFParams(:)     = FFParams(:) - FFParamsGrd(:)*stepsize/sqrt(rmsg)
        end if

    end do

    !===============================================================================
    if( istep .le. NOptSteps ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Final results', '-')
    else
        write(DEV_OUT,'(/,a)') ' >>> INFO: Maximum number of minimization steps was reached!'
        write(DEV_OUT,'(a,/)') ' >>> WARNING: Minimization was not completed!'
        call ffdev_utils_heading(DEV_OUT,'Intermediate results', '-')
    end if
    call write_header
    ! write final results
    call write_results(istep,FFError,rmsg,maxgrad,.true.) ! results to stdout

    if( PrintFinalGradient ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Final gradient', '-')
    end if

end subroutine opt_steepest_descent

!===============================================================================
! subroutine opt_lbfgs
!===============================================================================

subroutine opt_lbfgs

    use ffdev_ffopt_dat
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_parameters

    implicit none
    integer                 :: istep,alloc_status
    real(DEVDP)             :: rmsg, maxgrad, lasterror, eps, xtol
    integer                 :: iprint(2),iflag
    real(DEVDP),allocatable :: work(:)
    real(DEVDP),allocatable :: tmp_xg(:)
    ! --------------------------------------------------------------------------

    ! init required variables ====================
    lasterror     = 0.0d0
    iflag          = 0
    iprint(1)      = -1
    iprint(2)      = 1
    eps            = 1.0d-16
    xtol           = 1.0d-16   ! this is unrealistic criteria - we use xbp_own based on gradient

    ! allocate working array
    allocate(work(nactparms*(2*NumberOfCorrections+1)+2*NumberOfCorrections), &
             tmp_xg(nactparms), &
             stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate data for L-BFGS optimization!')
    end if

    ! perform minimization ========================
    do istep = 1, NOptSteps

        !===============================================================================
        ! get fit error
        call ffdev_parameters_error(FFParams,FFError,FFParamsGrd)

        !===============================================================================
        ! check all criteria
        rmsg = ffdev_fopt_rmsg(FFParamsGrd,maxgrad)

        if( istep .ne. 1 .and. abs(FFError%total - lasterror) .le. MinErrorChange ) then
            write(DEV_OUT,'(/,a,/)') ' >>> INFO: Error change is below treshold! Minimization was stoped.'
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: Last error change      : ', abs(FFError%total - lasterror)
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: Error change treshold : ', MinErrorChange
            exit
        end if

        if( abs(maxgrad) .le. MaxG .and. rmsg .le. MaxRMSG ) then
            write(DEV_OUT,'(/,a,/)') ' >>> INFO: Gradient tresholds were satisfied! Minimization was stoped.'
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: RMS of gradient                 : ', rmsg
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: RMS of gradient treshold        : ', MaxRMSG
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: Max gradient component          : ', abs(maxgrad)
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: Max gradient component treshold : ', MaxG
            write(DEV_OUT,'(a,E16.10)') ' >>> INFO: Last error change               : ', abs(FFError%total - lasterror)
            exit
        end if

        !===============================================================================
        ! print [intermediate] results (master node only)
        call write_results(istep,FFError,rmsg,maxgrad,.false.)

        ! if this is last step do not update coordinates and exit cycle
        if( istep .eq. NOptSteps ) exit

        !===============================================================================
        ! do L-BFGS minimization
        call LBFGS( nactparms, NumberOfCorrections, &
                    FFParams,FFError,FFParamsGrd,&
                    .false.,tmp_xg,iprint,eps,xtol,work,iflag)

        if( iflag .eq. 0 ) exit
        if( iflag .le. 0 ) then
            write(DEV_OUT,'(/,a,i2,/)') '>>> ERROR: Internal L-BFGS driver error! Code = ', iflag
            exit
        end if

        lasterror = FFError%total
    end do

    !===============================================================================
    if( istep .le. NOptSteps ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Final results', '-')
    else
        write(DEV_OUT,'(/,a)') ' >>> INFO: Maximum number of minimization steps was reached!'
        write(DEV_OUT,'(a,/)') ' >>> WARNING: Minimization was not completed!'
        call ffdev_utils_heading(DEV_OUT,'Intermediate results', '-')
    end if
    call write_header
    ! write final results
    call write_results(istep,FFError,rmsg,maxgrad,.true.) ! results to stdout

    if( PrintFinalGradient ) then
        write(DEV_OUT,*)
        call ffdev_utils_heading(DEV_OUT,'Final gradient', '-')
    end if

    deallocate(work,tmp_xg)

end subroutine opt_lbfgs

!===============================================================================
! subroutine opt_nlopt
!===============================================================================

subroutine opt_nlopt

    use ffdev_ffopt_dat
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_parameters

    implicit none

    include 'nlopt.f'

    integer                 :: istep,alloc_status
    real(DEVDP),allocatable :: tmp_xg(:)
    integer                 :: ires
    integer(8)              :: locoptid,oldlocoptid
    real(DEVDP)             :: final
    ! --------------------------------------------------------------------------

    NLoptID = 0
    call nlo_create(NLoptID,NLOpt_Method, nactparms)
    ! FIXME - there is a bug in passing of NLOpt_InitialStep to NLopt
    call nlo_set_initial_step(ires, NLoptID, real(NLOpt_InitialStep,DEVDP))
    call nlo_set_maxeval(ires, NLoptID, NOptSteps)
    call nlo_set_stopval(ires,NLoptID,0.0d0)
    call nlo_set_ftol_abs(ires, NLoptID, real(MinErrorChange,DEVDP))

    ! allocate working array
    allocate(tmp_xg(nactparms), stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate data for NLOPT optimization!')
    end if

    call ffdev_params_get_lower_bounds(tmp_xg)
    call nlo_set_lower_bounds(ires, NLoptID, tmp_xg)

    call ffdev_params_get_upper_bounds(tmp_xg)
    call nlo_set_upper_bounds(ires, NLoptID, tmp_xg)

    istep = 0
    call nlo_set_min_objective(ires, NLoptID, opt_nlopt_fce, istep)

    tmp_xg(:) = FFParams(:)

    call nlo_optimize(ires, NLoptID, tmp_xg, final)

    write(DEV_OUT,*)
    select case(ires)
        case(NLOPT_SUCCESS)
            write(DEV_OUT,10) ires, 'NLOPT_SUCCESS'
        case(NLOPT_STOPVAL_REACHED)
            write(DEV_OUT,10) ires, 'NLOPT_STOPVAL_REACHED'
        case(NLOPT_FTOL_REACHED)
            write(DEV_OUT,10) ires, 'NLOPT_FTOL_REACHED'
        case(NLOPT_XTOL_REACHED)
            write(DEV_OUT,10) ires, 'NLOPT_XTOL_REACHED'
        case(NLOPT_MAXEVAL_REACHED)
            write(DEV_OUT,10) ires, 'NLOPT_MAXEVAL_REACHED'
        case(NLOPT_MAXTIME_REACHED)
            write(DEV_OUT,10) ires, 'NLOPT_MAXTIME_REACHED'
        case(NLOPT_FAILURE)
            write(DEV_OUT,10) ires, 'NLOPT_FAILURE'
        case(NLOPT_INVALID_ARGS)
            write(DEV_OUT,10) ires, 'NLOPT_INVALID_ARGS'
        case(NLOPT_OUT_OF_MEMORY)
            write(DEV_OUT,10) ires, 'NLOPT_OUT_OF_MEMORY'
        case(NLOPT_ROUNDOFF_LIMITED)
            write(DEV_OUT,10) ires, 'NLOPT_ROUNDOFF_LIMITED'
        case(NLOPT_FORCED_STOP)
            write(DEV_OUT,10) ires, 'NLOPT_FORCED_STOP'
    end select

    deallocate(tmp_xg)

    call nlo_destroy(NLoptID)

10 format('NLOpt finished with return status = ',I6,' ',A)

end subroutine opt_nlopt

!===============================================================================
! subroutine ffdev_ffopt_header
!===============================================================================

subroutine opt_nlopt_fce(value, n, x, grad, need_gradient, istep)

    use ffdev_ffopt_dat
    use ffdev_parameters_dat
    use ffdev_parameters

    implicit none
    real(DEVDP)     :: value
    integer         :: n
    real(DEVDP)     :: x(n), grad(n)
    integer         :: need_gradient
    integer         :: istep
    ! --------------------------------------------
    real(DEVDP)     :: rmsg, maxgrad, lasterror, eps, xtol
    integer         :: ires
    ! --------------------------------------------------------------------------

    istep = istep + 1

    ! write(*,*) x

    FFParams(:) = x(:)
    if( need_gradient .gt. 0 ) then
        call ffdev_parameters_error(FFParams,FFError,FFParamsGrd)
        value = FFError%Total
        grad(:) = FFParamsGrd(:)
    else
        FFParamsGrd(:) = 0.0d0
        call ffdev_parameters_error_only(FFParams,FFError)
        value = FFError%Total
    end if

    rmsg = ffdev_fopt_rmsg(FFParamsGrd,maxgrad)
    call write_results(istep,FFError,rmsg,maxgrad,.false.)

    if( istep .eq. NOptSteps ) then
        call nlo_force_stop(ires, NLoptID)
    end if

end subroutine opt_nlopt_fce

!===============================================================================
! subroutine ffdev_ffopt_header
!===============================================================================

subroutine write_header()

    use ffdev_ffopt_dat
    use ffdev_topology
    use ffdev_geometry

    implicit none
    integer     :: major, minor, bugfix
    ! --------------------------------------------------------------------------

    select case(OptimizationMethod)
        case(MINIMIZATION_STEEPEST_DESCENT)
            write(DEV_OUT,10)
        case(MINIMIZATION_LBFGS)
            write(DEV_OUT,15)
        case(MINIMIZATION_NLOPT)
            call nloptv(major, minor, bugfix)
            write(DEV_OUT,17) major, minor, bugfix
    end select

    write(DEV_OUT,20)
    write(DEV_OUT,30)

 10 format('# Mode = Steepest Descent')
 15 format('# Mode = L-BFGS')
 17 format('# Mode = NLOPT v',I1,'.',I1,'.',I1)

 20 format('# STEP    Error       Err(Ene)    Err(Grad)    Err(Hess)      RMSG         maxG     ')
 30 format('#----- ------------ ------------ ------------ ------------ ------------ ------------')

end subroutine write_header

!===============================================================================
! subroutine write_results
!===============================================================================

subroutine write_results(istep,error,rmsg,maxgrad,done)

    use ffdev_ffopt_dat
    use ffdev_parameters_dat

    implicit none
    integer             :: istep
    type(FFERROR_TYPE)  :: error
    real(DEVDP)         :: rmsg
    real(DEVDP)         :: maxgrad
    logical             :: done
    ! -----------------------------------------------------------------------------

    ! write energies
    if( done .or. ((OutSamples .gt. 0) .and. (mod(istep,OutSamples) .eq. 0)) .or. (istep .eq. 1) ) then
        write(DEV_OUT,10) istep, error%total, error%energy, error%grad, error%hess, &
                          rmsg,maxgrad
    end if

 10 format(I6,1X,E12.5,1X,E12.5,1X,E12.5,1X,E12.5,1X,E12.5,1X,E12.5)

end subroutine write_results

!-------------------------------------------------------------------------------

! ==============================================================================
! function ffdev_gradient_rmsg
! ==============================================================================

real(DEVDP) function ffdev_fopt_rmsg(grd,maxgrad)

    use ffdev_parameters_dat

    implicit none
    real(DEVDP)     :: grd(:)
    real(DEVDP)     :: maxgrad
    ! ------------------------------
    integer         :: i
    real(DEVDP)     :: norm
    !------------------------------------------------------------------------------

    maxgrad = 0.0
    ffdev_fopt_rmsg = 0.0

    if( nactparms .eq. 0 ) return

    do i=1,nactparms
        norm = grd(i)**2
        if( abs(grd(i)) > abs(maxgrad) ) then
            maxgrad = grd(i)
        end if
        ffdev_fopt_rmsg = ffdev_fopt_rmsg + norm
    end do

    ffdev_fopt_rmsg = sqrt(ffdev_fopt_rmsg/real(nactparms))

    return

end function ffdev_fopt_rmsg

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_ffopt





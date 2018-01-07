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
    NLOpt_InitialStep   = 0.00001d0

! === [SA] =====================================================================
    OptSA_Temp   = 0.01d0    ! initial temperature
    OptSA_RT     = 0.85d0    ! the temperature reduction factor
    OptSA_EPS    = 0.0001    ! Error tolerance for termination
    OptSA_NS     = 20        ! Number of cycles
    OptSA_NT     = 10        ! Number of iterations before temperature reduction
    OptSA_NEPS   = 4         ! Number of final function values used to decide upon termination
    OptSA_MAXEVL = 80000     ! The maximum number of function evaluations.
    OptSA_IPRINT = 1         ! controls printing inside SA
    OptSA_ISEED  = 341723    ! random generator seed

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
        case(MINIMIZATION_SA)
            call opt_sa
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
            write(DEV_OUT,'(/,a,F16.4)') ' >>> INFO: Last error change      : ', abs(FFError%total - lasterror)
            write(DEV_OUT,'(a,F16.4)')   ' >>> INFO: Error change treshold : ', MinErrorChange
            write(DEV_OUT,'(a,/)') ' >>> INFO: Error change is below treshold! Minimization was stoped.'
            exit
        end if

        if( abs(maxgrad) .le. MaxG .and. rmsg .le. MaxRMSG ) then
            write(DEV_OUT,'(/,a,F16.4)') ' >>> INFO: RMS of gradient                 : ', rmsg
            write(DEV_OUT,'(a,F16.4)')   ' >>> INFO: RMS of gradient treshold        : ', MaxRMSG
            write(DEV_OUT,'(a,F16.4)')   ' >>> INFO: Max gradient component          : ', abs(maxgrad)
            write(DEV_OUT,'(a,F16.4)')   ' >>> INFO: Max gradient component treshold : ', MaxG
            write(DEV_OUT,'(a,F16.4)')   ' >>> INFO: Last error change               : ', abs(FFError%total - lasterror)
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
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: Last error change      : ', abs(FFError%total - lasterror)
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: Error change treshold : ', MinErrorChange
            exit
        end if

        if( abs(maxgrad) .le. MaxG .and. rmsg .le. MaxRMSG ) then
            write(DEV_OUT,'(/,a,/)') ' >>> INFO: Gradient tresholds were satisfied! Minimization was stoped.'
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: RMS of gradient                 : ', rmsg
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: RMS of gradient treshold        : ', MaxRMSG
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: Max gradient component          : ', abs(maxgrad)
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: Max gradient component treshold : ', MaxG
            write(DEV_OUT,'(a,F16.4)') ' >>> INFO: Last error change               : ', abs(FFError%total - lasterror)
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
    call nlo_create(NLoptID,NLOPT_LN_COBYLA, nactparms)
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
! subroutine opt_sa
!===============================================================================

subroutine opt_sa

    use ffdev_ffopt_dat
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_parameters

    implicit none
    integer     :: status, idx, i
    real(DEVDP) :: temp, opterror, error, r2, r2ave
    integer     :: nacc,nfcnev,nobds,ier
    ! --------------------------------------------------------------------------

    if( allocated(sa_lb) )         deallocate(sa_lb)
    if( allocated(sa_ub) )         deallocate(sa_ub)
    if( allocated(sa_c) )          deallocate(sa_c)
    if( allocated(sa_vm) )         deallocate(sa_vm)
    if( allocated(sa_xopt) )       deallocate(sa_xopt)
    if( allocated(sa_fstar) )      deallocate(sa_fstar)
    if( allocated(sa_xp) )         deallocate(sa_xp)
    if( allocated(sa_nacp) )       deallocate(sa_nacp)

    ! allocate helper arrays
    allocate(sa_lb(nactparms), &
          sa_ub(nactparms), &
          sa_c(nactparms), &
          sa_vm(nactparms), &
          sa_xopt(nactparms), &
          sa_fstar(OptSA_NEPS), &
          sa_xp(nactparms), &
          sa_nacp(nactparms), &
          stat = status)

    if( status .ne. 0 ) then
        call ffdev_utils_exit(DEV_OUT,1,'Unable to allocate memory for params in opt_sa!')
    end if

    ! init arrays
    sa_c(:)        =  2.0d0
    sa_vm(:)       = 10.0d0
    sa_xopt(:)     =  0.0d0
    sa_fstar(:)    =  0.0d0
    sa_xp(:)       =  0.0d0
    sa_nacp(:)     =  0.0d0

    ! set optimization boundary
    idx = 1
    do i=1,nparams
        if( .not. params(i)%enabled ) cycle
        select case( params(i)%realm )
            case(REALM_EOFFSET)
                sa_lb(idx)  = -1000.0
                sa_ub(idx)  =  1000.0
                idx = idx + 1
            case(REALM_DIH_V)
                sa_lb(idx)  =  -50.0
                sa_ub(idx)  =  50.0
                idx = idx + 1
            case(REALM_DIH_G)
                sa_lb(idx)  =  -2*DEV_PI
                sa_ub(idx)  =  2*DEV_PI
                idx = idx + 1
        end select
    end do


    ! initial temperature
    temp = OptSA_Temp

    write(*,*) FFParams
    write(*,*) sa_lb
    write(*,*) sa_ub

    ! start procedure
    call write_header
    call SimulatedAnnealing(nactparms,FFParams,.false.,OptSA_RT,OptSA_EPS,OptSA_NS, &
                         OptSA_NT,OptSA_NEPS,OptSA_MAXEVL,sa_lb,sa_ub,sa_c,OptSA_IPRINT, &
                         temp,sa_vm,sa_xopt,opterror, &
                         nacc,nfcnev,nobds,ier, &
                         sa_fstar,sa_xp,sa_nacp)

    write(DEV_OUT,*) 'ier=',ier

end subroutine opt_sa

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
! ------------------------------------------------------------------------------
!===============================================================================

subroutine SA_FCN(N,X,F)

    use ffdev_ffopt_dat
    use ffdev_parameters_dat
    use ffdev_parameters

    implicit none
    integer         :: N
    real(DEVDP)     :: X(:)
    real(DEVDP)     :: F
    ! ----------------------------------------
    real(DEVDP)     :: rmsg, maxgrad
    ! -----------------------------------------------------------------------------

    FFParams(:) = X(:)
    FFParamsGrd(:) = 0.0d0
    call ffdev_parameters_error_only(FFParams,FFError)
    F = FFError%Total

!    rmsg = ffdev_fopt_rmsg(FFParamsGrd,maxgrad)
!    call write_results(1,FFError,rmsg,maxgrad,.false.)

end subroutine SA_FCN

!===============================================================================
! Simulated annealing driver
!===============================================================================
! NOTE: the above GNU GPL licence is not applied for the code from this line till
! the end of this file

! ABSTRACT:
!   Simulated annealing is a global optimization method that distinguishes
!   between different local optima. Starting from an initial point, the
!   algorithm takes a step and the function is evaluated. When minimizing a
!   function, any downhill step is accepted and the process repeats from this
!   new point. An uphill step may be accepted. Thus, it can escape from local
!   optima. This uphill decision is made by the Metropolis criteria. As the
!   optimization process proceeds, the length of the steps decline and the
!   algorithm closes in on the global optimum. Since the algorithm makes very
!   few assumptions regarding the function to be optimized, it is quite
!   robust with respect to non-quadratic surfaces. The degree of robustness
!   can be adjusted by the user. In fact, simulated annealing can be used as
!   a local optimizer for difficult functions.
!
!   This implementation of simulated annealing was used in "Global Optimizatio
!   of Statistical Functions with Simulated Annealing," Goffe, Ferrier and
!   Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp.
!   65-100. Briefly, we found it competitive, if not superior, to multiple
!   restarts of conventional optimization routines for difficult optimization
!   problems.
!
!   For more information on this routine, contact its author:
!   Bill Goffe, bgoffe@whale.st.usm.edu

! NOTE: the code below was modified to be in accordance with our requirements

subroutine SimulatedAnnealing(N,X,MAX,RT,EPS,NS,NT,NEPS,MAXEVL,LB,UB,C, &
              IPRINT,T,VM,XOPT,FOPT,NACC,NFCNEV,NOBDS,IER, &
              FSTAR,XP,NACP)

!  Version: 3.2
!  Date: 1/22/94.
!  Differences compared to Version 2.0:
!     1. If a trial is out of bounds, a point is randomly selected
!        from LB(i) to UB(i). Unlike in version 2.0, this trial is
!        evaluated and is counted in acceptances and rejections.
!        All corresponding documentation was changed as well.
!  Differences compared to Version 3.0:
!     1. If VM(i) > (UB(i) - LB(i)), VM is set to UB(i) - LB(i).
!        The idea is that if T is high relative to LB & UB, most
!        points will be accepted, causing VM to rise. But, in this
!        situation, VM has little meaning; particularly if VM is
!        larger than the acceptable region. Setting VM to this size
!        still allows all parts of the allowable region to be selected.
!  Differences compared to Version 3.1:
!     1. Test made to see if the initial temperature is positive.
!     2. WRITE statements prettied up.
!     3. References to paper updated.
!
!  Synopsis:
!  This routine implements the continuous simulated annealing global
!  optimization algorithm described in Corana et al.'s article
!  "Minimizing Multimodal Functions of Continuous Variables with the
!  "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
!  no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
!  Software.
!
!  A very quick (perhaps too quick) overview of SA:
!     SA tries to find the global optimum of an N dimensional function.
!  It moves both up and downhill and as the optimization process
!  proceeds, it focuses on the most promising area.
!     To start, it randomly chooses a trial point within the step length
!  VM (a vector of length N) of the user selected starting point. The
!  function is evaluated at this trial point and its value is compared
!  to its value at the initial point.
!     In a maximization problem, all uphill moves are accepted and the
!  algorithm continues from that trial point. Downhill moves may be
!  accepted; the decision is made by the Metropolis criteria. It uses T
!  (temperature) and the size of the downhill move in a probabilistic
!  manner. The smaller T and the size of the downhill move are, the more
!  likely that move will be accepted. If the trial is accepted, the
!  algorithm moves on from that point. If it is rejected, another point
!  is chosen instead for a trial evaluation.
!     Each element of VM periodically adjusted so that half of all
!  function evaluations in that direction are accepted.
!     A fall in T is imposed upon the system with the RT variable by
!  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
!  downhill moves are less likely to be accepted and the percentage of
!  rejections rise. Given the scheme for the selection for VM, VM falls.
!  Thus, as T declines, VM falls and SA focuses upon the most promising
!  area for optimization.
!
!  The importance of the parameter T:
!     The parameter T is crucial in using SA successfully. It influences
!  VM, the step length over which the algorithm searches for optima. For
!  a small intial T, the step length may be too small; thus not enough
!  of the function might be evaluated to find the global optima. The user
!  should carefully examine VM in the intermediate output (set IPRINT =
!  1) to make sure that VM is appropriate. The relationship between the
!  initial temperature and the resulting step length is function
!  dependent.
!     To determine the starting temperature that is consistent with
!  optimizing a function, it is worthwhile to run a trial run first. Set
!  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
!  rises as well. Then select the T that produces a large enough VM.
!
!  For modifications to the algorithm and many details on its use,
!  (particularly for econometric applications) see Goffe, Ferrier
!  and Rogers, "Global Optimization of Statistical Functions with
!  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
!  Jan./Feb. 1994, pp. 65-100.
!  For more information, contact
!              Bill Goffe
!              Department of Economics and International Business
!              University of Southern Mississippi
!              Hattiesburg, MS  39506-5072
!              (601) 266-4484 (office)
!              (601) 266-4920 (fax)
!              bgoffe@whale.st.usm.edu (Internet)
!
!  As far as possible, the parameters here have the same name as in
!  the description of the algorithm on pp. 266-8 of Corana et al.
!
!  In this description, SP is single precision, DP is double precision,
!  INT is integer, L is logical and (N) denotes an array of length n.
!  Thus, DP(N) denotes a double precision array of length n.
!
!  Input Parameters:
!    Note: The suggested values generally come from Corana et al. To
!          drastically reduce runtime, see Goffe et al., pp. 90-1 for
!          suggestions on choosing the appropriate RT and NT.
!    N - Number of variables in the function to be optimized. (INT)
!    X - The starting values for the variables of the function to be
!        optimized. (DP(N))
!    MAX - Denotes whether the function should be maximized or
!          minimized. A true value denotes maximization while a false
!          value denotes minimization. Intermediate output (see IPRINT)
!          takes this into account. (L)
!    RT - The temperature reduction factor. The value suggested by
!         Corana et al. is .85. See Goffe et al. for more advice. (DP)
!    EPS - Error tolerance for termination. If the final function
!          values from the last neps temperatures differ from the
!          corresponding value at the current temperature by less than
!          EPS and the final function value at the current temperature
!          differs from the current optimal function value by less than
!          EPS, execution terminates and IER = 0 is returned. (EP)
!    NS - Number of cycles. After NS*N function evaluations, each
!         element of VM is adjusted so that approximately half of
!         all function evaluations are accepted. The suggested value
!         is 20. (INT)
!    NT - Number of iterations before temperature reduction. After
!         NT*NS*N function evaluations, temperature (T) is changed
!         by the factor RT. Value suggested by Corana et al. is
!         MAX(100, 5*N). See Goffe et al. for further advice. (INT)
!    NEPS - Number of final function values used to decide upon termi-
!           nation. See EPS. Suggested value is 4. (INT)
!    MAXEVL - The maximum number of function evaluations. If it is
!             exceeded, IER = 1. (INT)
!    LB - The lower bound for the allowable solution variables. (DP(N))
!    UB - The upper bound for the allowable solution variables. (DP(N))
!         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
!         I = 1, N, a point is from inside is randomly selected. This
!         This focuses the algorithm on the region inside UB and LB.
!         Unless the user wishes to concentrate the search to a par-
!         ticular region, UB and LB should be set to very large positive
!         and negative values, respectively. Note that the starting
!         vector X should be inside this region. Also note that LB and
!         UB are fixed in position, while VM is centered on the last
!         accepted trial set of variables that optimizes the function.
!    C - Vector that controls the step length adjustment. The suggested
!        value for all elements is 2.0. (DP(N))
!    IPRINT - controls printing inside SA. (INT)
!             Values: 0 - Nothing printed.
!                     1 - Function value for the starting value and
!                         summary results before each temperature
!                         reduction. This includes the optimal
!                         function value found so far, the total
!                         number of moves (broken up into uphill,
!                         downhill, accepted and rejected), the
!                         number of out of bounds trials, the
!                         number of new optima found at this
!                         temperature, the current optimal X and
!                         the step length VM. Note that there are
!                         N*NS*NT function evalutations before each
!                         temperature reduction. Finally, notice is
!                         is also given upon achieveing the termination
!                         criteria.
!                     2 - Each new step length (VM), the current optimal
!                         X (XOPT) and the current trial X (X). This
!                         gives the user some idea about how far X
!                         strays from XOPT as well as how VM is adapting
!                         to the function.
!                     3 - Each function evaluation, its acceptance or
!                         rejection and new optima. For many problems,
!                         this option will likely require a small tree
!                         if hard copy is used. This option is best
!                         used to learn about the algorithm. A small
!                         value for MAXEVL is thus recommended when
!                         using IPRINT = 3.
!             Suggested value: 1
!             Note: For a given value of IPRINT, the lower valued
!                   options (other than 0) are utilized.
!    ISEED1 - The first seed for the random number generator RANMAR.
!             0 .LE. ISEED1 .LE. 31328. (INT)
!    ISEED2 - The second seed for the random number generator RANMAR.
!             0 .LE. ISEED2 .LE. 30081. Different values for ISEED1
!             and ISEED2 will lead to an entirely different sequence
!             of trial points and decisions on downhill moves (when
!             maximizing). See Goffe et al. on how this can be used
!             to test the results of SA. (INT)
!
!  Input/Output Parameters:
!    T - On input, the initial temperature. See Goffe et al. for advice.
!        On output, the final temperature. (DP)
!    VM - The step length vector. On input it should encompass the
!         region of interest given the starting value X. For point
!         X(I), the next trial point is selected is from X(I) - VM(I)
!         to  X(I) + VM(I). Since VM is adjusted so that about half
!         of all points are accepted, the input value is not very
!         important (i.e. is the value is off, SA adjusts VM to the
!         correct value). (DP(N))
!
!  Output Parameters:
!    XOPT - The variables that optimize the function. (DP(N))
!    FOPT - The optimal value of the function. (DP)
!    NACC - The number of accepted function evaluations. (INT)
!    NFCNEV - The total number of function evaluations. In a minor
!             point, note that the first evaluation is not used in the
!             core of the algorithm; it simply initializes the
!             algorithm. (INT).
!    NOBDS - The total number of trial function evaluations that
!            would have been out of bounds of LB and UB. Note that
!            a trial point is randomly selected between LB and UB.
!            (INT)
!    IER - The error return number. (INT)
!          Values: 0 - Normal return; termination criteria achieved.
!                  1 - Number of function evaluations (NFCNEV) is
!                      greater than the maximum number (MAXEVL).
!                  2 - The starting value (X) is not inside the
!                      bounds (LB and UB).
!                  3 - The initial temperature is not positive.
!                  99 - Should not be seen; only used internally.
!
!  Work arrays that must be dimensioned in the calling routine:
!       RWK1 (DP(NEPS))  (FSTAR in SA)
!       RWK2 (DP(N))     (XP    "  " )
!       IWK  (INT(N))    (NACP  "  " )
!
!  Required Functions (included):
!    EXPREP - Replaces the function EXP to avoid under- and overflows.
!             It may have to be modified for non IBM-type main-
!             frames. (DP)
!    RMARIN - Initializes the random number generator RANMAR.
!    RANMAR - The actual random number generator. Note that
!             RMARIN must run first (SA does this). It produces uniform
!             random numbers on [0,1]. These routines are from
!             Usenet's comp.lang.fortran. For a reference, see
!             "Toward a Universal Random Number Generator"
!             by George Marsaglia and Arif Zaman, Florida State
!             University Report: FSU-SCRI-87-50 (1987).
!             It was later modified by F. James and published in
!             "A Review of Pseudo-random Number Generators." For
!             further information, contact stuart@ads.com. These
!             routines are designed to be portable on any machine
!             with a 24-bit or more mantissa. I have found it produces
!             identical results on a IBM 3081 and a Cray Y-MP.
!
!  Required Subroutines (included):
!    PRTVEC - Prints vectors.
!    PRT1 ... PRT10 - Prints intermediate output.
!    FCN - Function to be optimized. The form is
!            SUBROUTINE FCN(N,X,F)
!            INTEGER N
!            DOUBLE PRECISION  X(N), F
!            ...
!            function code with F = F(X)
!            ...
!            return
!            end
!          Note: This is the same form used in the multivariable
!          minimization algorithms in the IMSL edition 10 library.
!
!  Machine Specific Features:
!    1. EXPREP may have to be modified if used on non-IBM type main-
!       frames. Watch for under- and overflows in EXPREP.
!    2. Some FORMAT statements use G25.18; this may be excessive for
!       some machines.
!    3. RMARIN and RANMAR are designed to be protable; they should not
!       cause any problems.

 implicit none
 integer                    :: N
 real(DEVDP)                :: X(:)
 logical                    :: MAX
 real(DEVDP)                :: RT
 real(DEVDP)                :: EPS
 integer                    :: NS
 integer                    :: NT
 integer                    :: NEPS
 integer                    :: MAXEVL
 real(DEVDP)                :: LB(:)
 real(DEVDP)                :: UB(:)
 real(DEVDP)                :: C(:)
 integer                    :: IPRINT
 integer                    :: ISEED1
 integer                    :: ISEED2
 real(DEVDP)                :: T
 real(DEVDP)                :: VM(:)
 real(DEVDP)                :: XOPT(:)
 real(DEVDP)                :: FOPT
 integer                    :: NACC
 integer                    :: NFCNEV
 integer                    :: NOBDS
 integer                    :: IER
 real(DEVDP)                :: FSTAR(:)
 real(DEVDP)                :: XP(:)
 integer                    :: NACP(:)
 ! ------------------------------------------------
 !  Type all internal variables.
 real(DEVDP)    ::  F, FP, P, PP, RATIO
 integer        ::  NUP, NDOWN, NREJ, NNEW, LNOBDS, H, I, J, M
 logical        ::  QUIT
 !------------------------------------------------------------------------------

!  Initialize the random number generator RANMAR.
!      call RMARIN(ISEED1,ISEED2)

!  Set initial values.
      NACC = 0
      NOBDS = 0
      NFCNEV = 0
      IER = 99

      do 10, I = 1, N
         XOPT(I) = X(I)
         NACP(I) = 0
10    continue

      do 20, I = 1, NEPS
         FSTAR(I) = 1.0D+20
20    continue

!  If the initial temperature is not positive, notify the user and
!  return to the calling routine.
      if (T .LE. 0.0) then
         write(*,'(/,''  THE INITIAL TEMPERATURE IS NOT POSITIVE. RESET THE VARIABLE T. '')')
         IER = 3
         return
      end if

!  If the initial value is out of bounds, notify the user and return
!  to the calling routine.
      do 30, I = 1, N
         if ((X(I) .GT. UB(I)) .OR. (X(I) .LT. LB(I))) then
            call PRT1
            IER = 2
            return
         end if
30    continue

!  Evaluate the function with input X and return value as F.
      call SA_FCN(N,X,F)

!  If the function is to be minimized, switch the sign of the function.
!  Note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.
      if(.NOT. MAX) F = -F
      NFCNEV = NFCNEV + 1
      FOPT = F
      FSTAR(1) = F
      if(IPRINT .GE. 1) call PRT2(MAX,N,X,F,T)

!  Start the main loop. Note that it terminates if (i) the algorithm
!  succesfully optimizes the function or (ii) there are too many
!  function evaluations (more than MAXEVL).
100   NUP = 0
      NREJ = 0
      NNEW = 0
      NDOWN = 0
      LNOBDS = 0

      do 400, M = 1, NT
         do 300, J = 1, NS
            do 200, H = 1, N

!  Generate XP, the trial value of X. Note use of VM to choose XP.
               do 110, I = 1, N
                  if (I .EQ. H) then
                     XP(I) = X(I) + (RANMAR()*2.- 1.) * VM(I)
                  else
                     XP(I) = X(I)
                  end if

!  If XP is out of bounds, select a point in bounds for the trial.
                  if((XP(I) .LT. LB(I)) .OR. (XP(I) .GT. UB(I))) then
                    XP(I) = LB(I) + (UB(I) - LB(I))*RANMAR()
                    LNOBDS = LNOBDS + 1
                    NOBDS = NOBDS + 1
                    if(IPRINT .GE. 3) call PRT3(MAX,N,XP,X,FP,F)
                  end if
110            continue

!  Evaluate the function with the trial point XP and return as FP.
               call SA_FCN(N,XP,FP)
               if(.NOT. MAX) FP = -FP
               NFCNEV = NFCNEV + 1
               if(IPRINT .GE. 3) call PRT4(MAX,N,XP,X,FP,F)

!  If too many function evaluations occur, terminate the algorithm.
               if(NFCNEV .GE. MAXEVL) then
                  call PRT5
                  if (.NOT. MAX) FOPT = -FOPT
                  IER = 1
                  return
               end if

!  Accept the new point if the function value increases.
               if(FP .GE. F) then
                  if(IPRINT .GE. 3) then
                     write(*,'(''  POINT ACCEPTED'')')
                  end if
                  do 120, I = 1, N
                     X(I) = XP(I)
120               continue
                  F = FP
                  NACC = NACC + 1
                  NACP(H) = NACP(H) + 1
                  NUP = NUP + 1

!  If greater than any other point, record as new optimum.
                  if (FP .GT. FOPT) then
                  if(IPRINT .GE. 3) then
                     write(*,'(''  NEW OPTIMUM'')')
                  end if
                     do 130, I = 1, N
                        XOPT(I) = XP(I)
130                  continue
                     FOPT = FP
                     NNEW = NNEW + 1
                  end if

!  If the point is lower, use the Metropolis criteria to decide on
!  acceptance or rejection.
               else
                  P = EXPREP((FP - F)/T)
                  PP = RANMAR()
                  if (PP .LT. P) then
                     if(IPRINT .GE. 3) call PRT6(MAX)
                     do 140, I = 1, N
                        X(I) = XP(I)
140                  continue
                     F = FP
                     NACC = NACC + 1
                     NACP(H) = NACP(H) + 1
                     NDOWN = NDOWN + 1
                  else
                     NREJ = NREJ + 1
                     if(IPRINT .GE. 3) call PRT7(MAX)
                  end if
               end if

200         continue
300      continue

!  Adjust VM so that approximately half of all evaluations are accepted.
         do 310, I = 1, N
            RATIO = DFLOAT(NACP(I)) /DFLOAT(NS)
            if (RATIO .GT. .6) then
               VM(I) = VM(I)*(1. + C(I)*(RATIO - .6)/.4)
            else if (RATIO .LT. .4) then
               VM(I) = VM(I)/(1. + C(I)*((.4 - RATIO)/.4))
            end if
            if (VM(I) .GT. (UB(I)-LB(I))) then
               VM(I) = UB(I) - LB(I)
            end if
310      continue

         if(IPRINT .GE. 2) then
            call PRT8(N,VM,XOPT,X)
         end if

         do 320, I = 1, N
            NACP(I) = 0
320      continue

400   continue

      if(IPRINT .GE. 1) then
         call PRT9(MAX,NFCNEV,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
      end if

!  Check termination criteria.
      QUIT = .FALSE.
      FSTAR(1) = F
      if ((FOPT - FSTAR(1)) .LE. EPS) QUIT = .TRUE.
      do 410, I = 1, NEPS
         if (ABS(F - FSTAR(I)) .GT. EPS) QUIT = .FALSE.
410   continue

!  Terminate SA if appropriate.
      if (QUIT) then
         do 420, I = 1, N
            X(I) = XOPT(I)
420      continue
         IER = 0
         if (.NOT. MAX) FOPT = -FOPT
         if(IPRINT .GE. 1) call PRT10
         return
      end if

!  If termination criteria is not met, prepare for another loop.
      T = RT*T
      do 430, I = NEPS, 2, -1
         FSTAR(I) = FSTAR(I-1)
430   continue
      F = FOPT
      do 440, I = 1, N
         X(I) = XOPT(I)
440   continue

!  Loop again.
      GO TO 100

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

 real(DEVDP)  function  EXPREP(RDUM)
 !  This function replaces exp to avoid under- and overflows and is
 !  designed for IBM 370 type machines. It may be necessary to modify
 !  it for other machines. Note that the maximum and minimum values of
 !  EXPREP are such that they has no effect on the algorithm.

 implicit none
 real(DEVDP)    ::  RDUM
 ! -----------------------------------------------------------------------------

 if (RDUM .GT. 174.) then
    EXPREP = 3.69D+75
 else if (RDUM .LT. -180.) then
    EXPREP = 0.0
 else
    EXPREP = EXP(RDUM)
 end if

 return

end function

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

real(DEVDP) function ranmar()

 implicit none
 real(DEVDP)         :: random(1)
 !-------------------------------------------------------------------------------

 CALL RANDOM_NUMBER(random)
 ranmar = random(1)
 return

end function

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

!subroutine write_header()

! implicit none
! ! -----------------------------------------------------------------------------

! write(6,*)
! write(6,100)
! write(6,110)
! return

! 100 format('#   Step      Temp     best error   best R2    down  acc.up rej.up outbnd  new  ')
! 110 format('# -------- ---------- ------------ ---------- ------ ------ ------ ------ ------')

!end subroutine write_header

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT1
 !  This subroutine prints intermediate output, as does PRT2 through
 !  PRT10. Note that if SA is minimizing the function, the sign of the
 !  function value and the directions (up/down) are reversed in all
 !  output to correspond with the actual function optimization. This
 !  correction is because SA was written to maximize functions and
 !  it minimizes by maximizing the negative a function.
 ! -----------------------------------------------------------------------------
 implicit none
 ! -----------------------------------------------------------------------------

! to fix
!  write(*,'(/,''  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS '', &
!           /,''  (LB AND UB). EXECUTION TERMINATED WITHOUT ANY'', &
!           /,''  OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  '', &
!           /,''  LB(I) .LT. X(I) .LT. UB(I), I = 1, N. ''/)')

 return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT2(MAX,N,X,F,T)

 implicit none
 logical                    ::  MAX
 integer                    ::  N
 real(DEVDP)                ::  X(:)
 real(DEVDP)                ::  F
 real(DEVDP)                ::  T
 !------------------------------------------------
 real(DEVDP)                ::  R2
 ! -----------------------------------------------------------------------------

!  write(*,'(''  '')')
!  call PRTVEC(X,N,'INITIAL X')
!  if (MAX) then
!     write(*,'(''  INITIAL F: '',/, G25.18)') F
!  else
!     write(*,'(''  INITIAL F: '',/, G25.18)') -F
!  end if

! call eem_params_update(prms,x)

! if( OptSA_UseDirectChiOpt ) then
!    call params_solver_chi_optimize(prms,set,.true.)
! end if

! call fit_error_get_values(prms,set)
! r2 = fit_error_get_R2(set)
 R2 = 0.0
 write(6,100) 1,T,-F,R2

 return

100 format(I10,1X,E10.3,1X,E12.5,1X,F10.7)

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT3(MAX,N,XP,X,FP,F)

 implicit none
 logical     ::  MAX
 integer     ::  N
 real(DEVDP) ::  XP(:)
 real(DEVDP) ::  X(:)
 real(DEVDP) ::  FP
 real(DEVDP) ::  F
 ! -----------------------------------------------------------------------------

 write(*,'(''  '')')
 call PRTVEC(X,N,'CURRENT X')
 if (MAX) then
    write(*,'(''  CURRENT F: '',G25.18)') F
 else
    write(*,'(''  CURRENT F: '',G25.18)') -F
 end if
 call PRTVEC(XP,N,'TRIAL X')
 write(*,'(''  POINT REJECTED SINCE OUT OF BOUNDS'')')

 return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT4(MAX,N,XP,X,FP,F)

 implicit none
 logical        ::  MAX
 integer        ::  N
 real(DEVDP)    ::  XP(:)
 real(DEVDP)    ::  X(:)
 real(DEVDP)    ::  FP
 real(DEVDP)    ::  F
 ! -----------------------------------------------------------------------------

 write(*,'(''  '')')
 call PRTVEC(X,N,'CURRENT X')
 if (MAX) then
    write(*,'(''  CURRENT F: '',G25.18)') F
    call PRTVEC(XP,N,'TRIAL X')
    write(*,'(''  RESULTING F: '',G25.18)') FP
 else
    write(*,'(''  CURRENT F: '',G25.18)') -F
    call PRTVEC(XP,N,'TRIAL X')
    write(*,'(''  RESULTING F: '',G25.18)') -FP
 end if

 return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT5

 implicit none
 ! -----------------------------------------------------------------------------

! to fix
!       write(*,'(/,''  TOO MANY FUNCTION EVALUATIONS; CONSIDER ''
!      1          /,''  INCREASING MAXEVL OR EPS, OR DECREASING ''
!      2          /,''  NT OR RT. THESE RESULTS ARE LIKELY TO BE ''
!      3          /,''  POOR.'',/)')

      return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT6(MAX)

 implicit none
 logical    ::  MAX
 ! -----------------------------------------------------------------------------

 if (MAX) then
    write(*,'(''  THOUGH LOWER, POINT ACCEPTED'')')
 else
    write(*,'(''  THOUGH HIGHER, POINT ACCEPTED'')')
 end if

 return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT7(MAX)

 implicit none
 logical    ::  MAX
 ! -----------------------------------------------------------------------------

 if (MAX) then
    write(*,'(''  LOWER POINT REJECTED'')')
 else
    write(*,'(''  HIGHER POINT REJECTED'')')
 end if

 return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT8(N,VM,XOPT,X)

 implicit none
 integer        :: N
 real(DEVDP)    :: VM(:)
 real(DEVDP)    :: XOPT(:)
 real(DEVDP)    :: X(:)
 ! -----------------------------------------------------------------------------

!  write(*,'(/,'' INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT'',/)')
!  call PRTVEC(VM,N,'NEW STEP LENGTH (VM)')
!  call PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
!  call PRTVEC(X,N,'CURRENT X')
!  write(*,'('' '')')

 return

end subroutine

!===============================================================================
! ------------------------------------------------------------------------------
!===============================================================================

subroutine PRT9(MAX,NFCNEV,N,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)

 implicit none
 logical                    ::  MAX
 integer                    ::  NFCNEV
 integer                    ::  N
 real(DEVDP)                ::  T
 real(DEVDP)                ::  XOPT(:)
 real(DEVDP)                ::  VM(:)
 real(DEVDP)                ::  FOPT
 integer                    ::  NUP
 integer                    ::  NDOWN
 integer                    ::  NREJ
 integer                    ::  LNOBDS
 integer                    ::  NNEW
 !------------------------------------------------
 integer                    ::  TOTMOV
 real(DEVDP)                ::  R2
 ! -----------------------------------------------------------------------------

! TOTMOV = NUP + NDOWN + NREJ
!
! write(*,'(/,'' INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION'',/)')
! write(*,'(''  CURRENT TEMPERATURE:            '',G12.5)') T
!     if (MAX) then
!         write(*,'(''  MAX FUNCTION VALUE SO FAR:  '',G25.18)') FOPT
!         write(*,'(''  TOTAL MOVES:                '',I8)') TOTMOV
!         write(*,'(''     UPHILL:                  '',I8)') NUP
!         write(*,'(''     ACCEPTED DOWNHILL:       '',I8)') NDOWN
!         write(*,'(''     REJECTED DOWNHILL:       '',I8)') NREJ
!         write(*,'(''  OUT OF BOUNDS TRIALS:       '',I8)') LNOBDS
!         write(*,'(''  NEW MAXIMA THIS TEMPERATURE:'',I8)') NNEW
!     else
!         write(*,'(''  MIN FUNCTION VALUE SO FAR:  '',G25.18)') -FOPT
!         write(*,'(''  TOTAL MOVES:                '',I8)') TOTMOV
!         write(*,'(''     DOWNHILL:                '',I8)')  NUP
!         write(*,'(''     ACCEPTED UPHILL:         '',I8)')  NDOWN
!         write(*,'(''     REJECTED UPHILL:         '',I8)')  NREJ
!         write(*,'(''  TRIALS OUT OF BOUNDS:       '',I8)')  LNOBDS
!         write(*,'(''  NEW MINIMA THIS TEMPERATURE:'',I8)')  NNEW
!     end if
!     call PRTVEC(XOPT,N,'CURRENT OPTIMAL X')
!     call PRTVEC(VM,N,'STEP LENGTH (VM)')
!     write(*,'('' '')')

! call eem_params_update(prms,xopt)

! if( OptSA_UseDirectChiOpt ) then
!    call params_solver_chi_optimize(prms,set,.true.)
! end if

! call fit_error_get_values(prms,set)
! r2 = fit_error_get_R2(set)

 write(6,100) NFCNEV,T,-FOPT,R2,NUP,NDOWN,NREJ,LNOBDS,NNEW

 return

100 format(I10,1X,E10.3,1X,E12.5,1X,F10.7,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X)

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRT10

 implicit none
 ! -----------------------------------------------------------------------------

 write(*,'(/,''  SA ACHIEVED TERMINATION CRITERIA. IER = 0. '')')
 return

end subroutine

!===============================================================================
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!===============================================================================

subroutine PRTVEC(VECTOR,NCOLS,NAME)
 !  This subroutine prints the double precision vector named VECTOR.
 !  Elements 1 thru NCOLS will be printed. NAME is a character variable
 !  that describes VECTOR. Note that if NAME is given in the call to
 !  PRTVEC, it must be enclosed in quotes. If there are more than 10
 !  elements in VECTOR, 10 elements will be printed on each line.

 implicit none
 real(DEVDP)    :: VECTOR(NCOLS)
 integer        :: NCOLS
 CHARACTER(*)   :: NAME
 ! -----------------------------------------------
 integer        :: I, J, LL, LINES
 ! -----------------------------------------------------------------------------

      write(*,1001) NAME

      if (NCOLS .GT. 10) then
         LINES = INT(NCOLS/10.)

         do 100, I = 1, LINES
            LL = 10*(I - 1)
            write(*,1000) (VECTOR(J),J = 1+LL, 10+LL)
  100    continue

         write(*,1000) (VECTOR(J),J = 11+LL, NCOLS)
      else
         write(*,1000) (VECTOR(J),J = 1, NCOLS)
      end if


 return

 1000 FORMAT( 10(G12.5,1X))
 1001 FORMAT(/,25X,A)

end subroutine

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_ffopt





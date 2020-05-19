! ==============================================================================
! PMFLib - Library Supporting Potential of Mean Force Calculations
! ------------------------------------------------------------------------------
!    Copyright (C) 2011-2015 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2013-2015 Letif Mones, lam81@cam.ac.uk
!    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
!
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor,
!    Boston, MA  02110-1301  USA
! ==============================================================================

module ffdev_timers

use ffdev_sizes
use ffdev_constants
use ffdev_variables

implicit none

    integer     :: FFDEV_POT_TIMER                              = -30
        integer     :: FFDEV_INITIALIZATION_TIMER               = -31
        integer     :: FFDEV_FFOPT_TIMER                        = -32
            integer     :: FFDEV_ERRORS_TIMER                   = -35
            integer     :: FFDEV_POT_ENERGY_TIMER               = -40
                integer     :: FFDEV_POT_NB_ENERGY_TIMER        = -41
            integer     :: FFDEV_TARGETSET_ALL_TIMER            = -42
            integer     :: FFDEV_GEOOPT_TIMER                   = -43
                integer     :: FFDEV_POT_GRADIENT_TIMER         = -50
                    integer     :: FFDEV_POT_NB_GRADIENT_TIMER  = -51
            integer     :: FFDEV_POT_HESSIAN_TIMER              = -60

contains

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_timers_init_top

    use smf_profiling

    implicit none
    ! -------------------------------------------------------------------------

    call init_profiling(DEV_OUT,50)
    call start_timer(TOTAL_TIMER)

end subroutine ffdev_timers_init_top

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_timers_init

 use smf_profiling_dat
 use smf_profiling

 implicit none
 ! -----------------------------------------------------------------------------

 ! add standard timers --------------------------------
    FFDEV_INITIALIZATION_TIMER    = add_timer(TOTAL_TIMER,'Initialization')
    FFDEV_FFOPT_TIMER             = add_timer(TOTAL_TIMER,'FFOpt')
        FFDEV_ERRORS_TIMER            = add_timer(FFDEV_FFOPT_TIMER,'Errors')
        FFDEV_POT_ENERGY_TIMER        = add_timer(FFDEV_FFOPT_TIMER,'MM Energy')
            FFDEV_POT_NB_ENERGY_TIMER         = add_timer(FFDEV_POT_ENERGY_TIMER,'NB energy')
#ifdef _OPENMP
        FFDEV_TARGETSET_ALL_TIMER     = add_timer(FFDEV_FFOPT_TIMER,'TargetSet Evaluation')
        FFDEV_POT_HESSIAN_TIMER       = add_timer(FFDEV_FFOPT_TIMER,'MM Hessian')
#else
        FFDEV_GEOOPT_TIMER            = add_timer(FFDEV_FFOPT_TIMER,'GeoOpt')
            FFDEV_POT_GRADIENT_TIMER      = add_timer(FFDEV_GEOOPT_TIMER,'MM Gradient')
                FFDEV_POT_NB_GRADIENT_TIMER       = add_timer(FFDEV_POT_GRADIENT_TIMER,'NB gradient')
        FFDEV_POT_HESSIAN_TIMER       = add_timer(FFDEV_FFOPT_TIMER,'MM Hessian')
#endif

end subroutine ffdev_timers_init

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_timers_start_timer(id)

    use smf_profiling

    implicit none
    integer        :: id
    ! -------------------------------------------------------------------------

    ! skip timers from OpenMP regions, which were not initialized
    if( id .le. 0 ) return

    call start_timer(id)

end subroutine ffdev_timers_start_timer

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_timers_stop_timer(id)

    use smf_profiling

    implicit none
    integer        :: id
    ! -------------------------------------------------------------------------

    ! skip timers from OpenMP regions, which were not initialized
    if( id .le. 0 ) return

    call stop_timer(id)

end subroutine ffdev_timers_stop_timer

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

subroutine ffdev_timers_finalize(do_profiling)

    use smf_profiling

    implicit none
    logical :: do_profiling
    ! -------------------------------------------------------------------------

    call stop_timer(TOTAL_TIMER)
    if( do_profiling ) then
        call write_timing
        call finalize_profiling
    end if

end subroutine ffdev_timers_finalize

!===============================================================================
!-------------------------------------------------------------------------------
!===============================================================================

end module ffdev_timers

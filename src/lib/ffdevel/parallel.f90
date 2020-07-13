! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_parallel

use ffdev_constants
use ffdev_variables

contains

!===============================================================================
! subroutine ffdev_parallel_init
!===============================================================================

subroutine ffdev_parallel_init()

    use ffdev_parallel_dat
    use prmfile
    !$ use omp_lib

    implicit none
    ! --------------------------------------------------------------------------

    ParallelMode = .false.
    !$ ParallelMode = .true.

    NumOfThreads = 1
    !$ NumOfThreads = omp_get_max_threads()

    write(DEV_OUT,*)
    write(DEV_OUT,10)

    write(DEV_OUT,20) prmfile_onoff(ParallelErrorBlock)
    write(DEV_OUT,30) prmfile_onoff(ParallelGeoOptBlock)
    write(DEV_OUT,40) prmfile_onoff(ParallelMode)
    write(DEV_OUT,50) NumOfThreads

    write(DEV_OUT,10)

 10 format('#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA#PARA')
 20 format(' Parallel execution of Error block : ',A)
 30 format(' Parallel execution of GeoOpt block: ',A)
 40 format(' OpenMP enabled   : ',A)
 50 format(' Number of threads: ',I3)

end subroutine ffdev_parallel_init

!===============================================================================

end module ffdev_parallel

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

module ffdev_buried

use ffdev_buried_dat
use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_buried_run_stat
! ==============================================================================

subroutine ffdev_buried_run_stat()

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    integer     :: ti, ai, i, j, alloc_stat
    real(DEVDP) :: rms, sexp
    logical     :: buried_data_loaded
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    call ffdev_utils_heading(DEV_OUT,'BURIED ATOMS', ':')

    ! allocate buried atoms array
    if( .not. allocated(buried_atoms) ) then
        allocate( buried_atoms(ntypes), stat = alloc_stat )
        if( alloc_stat .ne. 0 ) then
            call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate buried_atoms')
        end if
    end if

    ! clear data
    do ti=1,ntypes
        buried_atoms(ti)%expave = 0
        buried_atoms(ti)%expsig = 0
        buried_atoms(ti)%weight = 0
        buried_atoms(ti)%num    = 0
    end do

    ! initialize
    buried_data_loaded = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%sup_surf_loaded ) then
                buried_data_loaded = .true.
                exit
            end if
        end do
    end do
    if( .not. buried_data_loaded ) then
        ! mark all atoms exposed if no SURF
        do ti=1,ntypes
            if( buried_atoms(ti)%num .eq. 0 ) then
                buried_atoms(ti)%expave = 1.0
            end if
        end do
        write(DEV_OUT,10)
        return
    end if

! gather data
    do i=1,nsets
        do j=1,sets(i)%ngeos
            ! do we have data?
            if( .not. sets(i)%geo(j)%sup_surf_loaded ) cycle

            do ai=1,sets(i)%geo(j)%natoms
                ! get types
                ti = sets(i)%top%atom_types(sets(i)%top%atoms(ai)%typeid)%glbtypeid

                ! by default - atom is not exposed
                sexp = 0.0
                if( sets(i)%geo(j)%sup_surf_atr(ai) .gt. 0.0d0 ) then
                    select case(surface_mode)
                        case(SURF_SESA)
                            sexp = sets(i)%geo(j)%sup_surf_ses(ai) / &
                                   (4.0 * DEV_PI * sets(i)%geo(j)%sup_surf_atr(ai) ** 2)
                        case(SURF_SASA)
                            sexp = sets(i)%geo(j)%sup_surf_sas(ai) / &
                                   (4.0 * DEV_PI * (sets(i)%geo(j)%sup_surf_atr(ai) + ProbeR) ** 2)
                        case default
                    end select
                end if
                if( sexp .lt. 0.0d0 ) sexp = 0.0d0
                if( sexp .gt. 1.0d0 ) sexp = 1.0d0

                ! accumulate data
                buried_atoms(ti)%expave  = buried_atoms(ti)%expave + sexp
                buried_atoms(ti)%expsig  = buried_atoms(ti)%expsig + sexp**2
                buried_atoms(ti)%num   = buried_atoms(ti)%num + 1
            end do
        end do
    end do

    ! mark all atoms exposed if no SURF
    do ti=1,ntypes
        if( buried_atoms(ti)%num .eq. 0 ) then
            buried_atoms(ti)%expave = 1.0
        end if
    end do

! finish data
    do ti=1,ntypes
        if( buried_atoms(ti)%num .gt. 0 ) then
            rms = buried_atoms(ti)%num*buried_atoms(ti)%expsig - buried_atoms(ti)%expave**2;
            if( rms .gt. 0.0d0 ) then
                rms = sqrt(rms) / real(buried_atoms(ti)%num)
            else
                rms = 0.0d0
            end if
            buried_atoms(ti)%expsig = rms
            buried_atoms(ti)%expave = buried_atoms(ti)%expave / real(buried_atoms(ti)%num)
        end if
        buried_atoms(ti)%weight = 1.0d0 - 1.0d0 / ( exp( (buried_atoms(ti)%expave - BuriedExp0) * BuriedBeta) + 1.0d0)
    end do

    call ffdev_buried_print

 10 format('>>> No information about buried atoms ....')

end subroutine ffdev_buried_run_stat

! ==============================================================================
! subroutine ffdev_buried_print
! ==============================================================================

subroutine ffdev_buried_print()

    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_mmd3

    implicit none
    integer     :: ti
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,20) ffdev_buried_surf_mode_to_string(surface_mode)
    write(DEV_OUT,25) ProbeR
    write(DEV_OUT,30) BuriedExp0
    write(DEV_OUT,40) BuriedBeta

    write(DEV_OUT,*)
    write(DEV_OUT,130)
    write(DEV_OUT,140)
    do ti=1,ntypes
        write(DEV_OUT,150) trim(types(ti)%name),buried_atoms(ti)%num, &
                          buried_atoms(ti)%expave, buried_atoms(ti)%expsig, &
                          buried_atoms(ti)%weight
    end do
   write(DEV_OUT,140)

 20 format("Surface mode = ",A)
 25 format("Probe radius = ",F10.4)
 30 format("Weight exp0  = ",F10.4)
 40 format("Weight beta  = ",F10.4)

130 format('# Type Number  <exp> s(exp) weight')
140 format('# ---- ------ ------ ------ ------')
150 format(2X,A4,1X,I6,1X,F6.4,1X,F6.4,1X,F6.4)

end subroutine ffdev_buried_print

! ==============================================================================
! subroutine ffdev_buried_surf_mode_from_string
! ==============================================================================

integer function ffdev_buried_surf_mode_from_string(string)

    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('SESA')
            ffdev_buried_surf_mode_from_string = SURF_SESA
        case('SASA')
            ffdev_buried_surf_mode_from_string = SURF_SASA
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_buried_surf_mode_from_string!')
    end select

end function ffdev_buried_surf_mode_from_string

! ==============================================================================
! subroutine ffdev_buried_surf_mode_to_string
! ==============================================================================

character(80) function ffdev_buried_surf_mode_to_string(mode)

    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(SURF_SESA)
            ffdev_buried_surf_mode_to_string = 'SESA'
        case(SURF_SASA)
            ffdev_buried_surf_mode_to_string = 'SASA'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_buried_surf_mode_to_string!')
    end select

end function ffdev_buried_surf_mode_to_string

! ------------------------------------------------------------------------------

end module ffdev_buried

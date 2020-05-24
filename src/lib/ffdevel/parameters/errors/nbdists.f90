! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_nbdists

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_nbdists_init
! ==============================================================================

subroutine ffdev_err_nbdists_init

    use ffdev_err_nbdists_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableNBDistsError        = .false.
    PrintNBDistsErrorSummary  = .false.
    NBDistsErrorWeight        = 1.0

    NBDistanceSWPosition         = 4.0
    NBDistanceSWAlpha            = 1.0

end subroutine ffdev_err_nbdists_init

! ==============================================================================
! subroutine ffdev_err_nbdists_error
! ==============================================================================

subroutine ffdev_err_nbdists_error(error)

    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat
    use ffdev_err_nbdists_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j,q,ai,aj,num
    real(DEVDP)         :: err,seterrnbdists
    real(DEVDP)         :: d0,dt,sw,swsum
    ! --------------------------------------------------------------------------

    error%nbdists = 0.0

    seterrnbdists = 0.0
    swsum = 0
    num = 0

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_crd_optimized .and. sets(i)%top%nfragments .gt. 1) then
                do q=1,sets(i)%top%nb_size
                    ai = sets(i)%top%nb_list(q)%ai
                    aj = sets(i)%top%nb_list(q)%aj

                    if( sets(i)%top%atoms(ai)%frgid .eq. sets(i)%top%atoms(aj)%frgid ) cycle



                    d0 = ffdev_geometry_get_length(sets(i)%geo(j)%crd,ai,aj)
                    dt = ffdev_geometry_get_length(sets(i)%geo(j)%trg_crd,ai,aj)
                    err = d0 - dt

                    ! calculate switch function
                    sw = 1.0d0 / (1.0d0 + exp( NBDistanceSWAlpha*(dt - NBDistanceSWPosition) ) )

                    seterrnbdists = seterrnbdists + sets(i)%geo(j)%weight * sw * err**2
                    swsum = swsum + sw
                    num = num + 1
                end do
            end if
        end do
    end do

    if( swsum .gt. 0 ) then
        error%nbdists = sqrt(seterrnbdists/swsum)
    end if

end subroutine ffdev_err_nbdists_error

! ==============================================================================
! subroutine ffdev_err_nbdists_summary
! ==============================================================================

subroutine ffdev_err_nbdists_summary(top,geo,printsum)

    use ffdev_topology_dat
    use ffdev_geometry_dat
    use ffdev_geometry
    use ffdev_err_nbdists_dat

    implicit none
    type(TOPOLOGY)  :: top
    type(GEOMETRY)  :: geo
    logical         :: printsum
    ! --------------------------------------------
    integer         :: q, ai, aj
    real(DEVDP)     :: d1, d2, diff, sdiff, sw, swsum
    real(DEVDP)     :: serr, lerr,aerr,rmse
    ! --------------------------------------------------------------------------

    if( .not. geo%trg_crd_optimized ) return

    if( printsum .eqv. .false. ) then
        printsum = .true.
        return
    end if

    write(DEV_OUT,*)
    write(DEV_OUT,100)
    write(DEV_OUT,110)
    write(DEV_OUT,125)
    write(DEV_OUT,130)

    serr = 100d0
    lerr = 0.0d0
    aerr = 0.0d0
    rmse = 0.0d0
    swsum= 0.0d0

    do q=1,top%nb_size
        ai = top%nb_list(q)%ai
        aj = top%nb_list(q)%aj

        if( top%atoms(ai)%frgid .eq. top%atoms(aj)%frgid ) cycle

        d1 = ffdev_geometry_get_length(geo%trg_crd,ai,aj)
        d2 = ffdev_geometry_get_length(geo%crd,ai,aj)
        diff = d2 - d1

        ! calculate switch function
        sw = 1.0d0 / (1.0d0 + exp( NBDistanceSWAlpha*(d1 - NBDistanceSWPosition) ) )
        swsum = swsum + sw
        sdiff = diff * sw

        if( serr .gt. abs(sdiff) ) serr = abs(sdiff)
        if( lerr .lt. abs(sdiff) ) lerr = abs(sdiff)
        aerr = aerr + abs(sdiff)
        rmse = rmse + sdiff**2

        write(DEV_OUT,140) ai, top%atoms(ai)%name, top%atom_types(top%atoms(ai)%typeid)%name, &
                            top%atoms(ai)%residx, top%atoms(ai)%resname, &
                            aj, top%atoms(aj)%name, top%atom_types(top%atoms(aj)%typeid)%name, &
                            top%atoms(aj)%residx, top%atoms(aj)%resname, &
                            d1,d2,diff,sw,sdiff
    end do

    if( swsum .gt. 0.0d0 ) then
        aerr = aerr / swsum
        rmse = sqrt(rmse / swsum)
    end if

    write(DEV_OUT,110)
    write(DEV_OUT,150) serr
    write(DEV_OUT,160) lerr
    write(DEV_OUT,170) aerr
    write(DEV_OUT,180) rmse

    100 format('# Individual Distances')
    110 format('# --------------------------- = ----------------------------- -------------------------------------------------')
    125 format('# Indx Name Type  RIdx  RName    Indx  Name Type  RIdx  RName  d#TRG(1)   d#MM(2) diff(2-1)  sw fce    err(2-1)')
    130 format('# ---- ---- ---- ------ ----- = ------ ---- ---- ------ ----- --------- --------- --------- --------- ---------')
    140 format(I6,1X,A4,1X,A4,1X,I6,1X,A5,3X,I6,1X,A4,1X,A4,1X,I6,1X,A5,1X,F9.4,1X,F9.4,1X,F9.4,1X,F9.4,1X,F9.4)
    150 format('# Minimum unsigned difference (SUD)  = ',F9.4)
    160 format('# Largest unsigned difference (MUD)  = ',F9.4)
    170 format('# Average usigned difference (AD)    = ',F9.4)
    180 format('# Root mean square difference (RMSD) = ',F9.4)

end subroutine ffdev_err_nbdists_summary

! ------------------------------------------------------------------------------

end module ffdev_err_nbdists

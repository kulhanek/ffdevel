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

module ffdev_err_nbpnl

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_nbpnl_init
! ==============================================================================

subroutine ffdev_err_nbpnl_init

    use ffdev_err_nbpnl_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableNBPnlError        = .false.
    PrintNBPnlErrorSummary  = .false.
    NBPnlErrorWeight        = 1.0
    NBPnlErrorTempFactor    = 0.0

end subroutine ffdev_err_nbpnl_init

! ==============================================================================
! subroutine ffdev_err_nbpnl_error
! ==============================================================================

subroutine ffdev_err_nbpnl_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_nbpnl_dat
    use ffdev_nb2nb_dat
    use ffdev_nb2nb
    use ffdev_parameters_dat
    use ffdev_targetset_dat

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i,j
    real(DEVDP)         :: dr,r,r0,eps,tote,totw,dif,w,elj,enb,paire,pairw
    logical             :: tia,tja
    ! --------------------------------------------------------------------------

    error%nbpnl = 0.0d0

    dr = NB2NBCutoffR / real(NB2NBNBins,DEVDP)

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nnb_types

        ! is it active?
        tia = .false.
        do j=1,nparams
            if( .not. params(j)%enabled ) cycle
            if( (params(j)%realm .ne. REALM_VDW_EPS) .and. (params(j)%realm .ne. REALM_VDW_R0) ) cycle
            if( ((nb_types(i)%gti .eq. params(j)%ti) .or. (nb_types(i)%gti .eq. params(j)%tj)) ) then
                tia = .true.
                exit
            end if
        end do
        tja = .false.
        do j=1,nparams
            if( .not. params(j)%enabled ) cycle
            if( (params(j)%realm .ne. REALM_VDW_EPS) .and. (params(j)%realm .ne. REALM_VDW_R0) ) cycle
            if( ((nb_types(i)%gtj .eq. params(j)%ti) .or. (nb_types(i)%gtj .eq. params(j)%tj)) ) then
                tja = .true.
                exit
            end if
        end do
        if( .not. (tia .and. tja) ) cycle

        r0  = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%r0
        eps = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%eps

        r = nb_types(i)%SigNB

        paire = 0.0d0
        pairw = 0.0d0

        do j=1,NB2NBNBins
            elj = ffdev_nb2nb_ljene_vdw_ene(r0,eps,r)
            enb = nb_types(i)%NBPot(j)
            dif = elj-enb
            w = 1.0d0
            if( NBPnlErrorTempFactor .gt. 0.0d0 ) then
                w = exp(-(enb+nb_types(i)%EpsNB)/(NBPnlErrorTempFactor*DEV_Rgas))
            end if
            paire = paire + w*dif**2
            pairw = pairw + w
            r = r + dr

            ! write(*,*) r,enb,elj
        end do

        if( pairw .gt. 0.0d0 ) then
            paire = sqrt(paire/pairw)
        end if

        nb_types(i)%errval  = paire
        tote = tote + nb_types(i)%num * paire
        totw = totw + nb_types(i)%num
    end do

    if( totw .gt. 0.0d0 ) then
        error%nbpnl = tote / totw
    end if

end subroutine ffdev_err_nbpnl_error

! ==============================================================================
! subroutine ffdev_err_nbpnl_summary
! ==============================================================================

subroutine ffdev_err_nbpnl_summary()

    use ffdev_err_nbpnl_dat
    use ffdev_utils
    use ffdev_nb2nb_dat
    use ffdev_targetset_dat
    use ffdev_parameters_dat
    use ffdev_nb2nb

    implicit none
    integer         :: i,j
    real(DEVDP)     :: tote,totw
    logical         :: tia,tja
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nnb_types

        ! is it active?
        tia = .false.
        do j=1,nparams
            if( .not. params(j)%enabled ) cycle
            if( (params(j)%realm .ne. REALM_VDW_EPS) .and. (params(j)%realm .ne. REALM_VDW_R0) ) cycle
            if( ((nb_types(i)%gti .eq. params(j)%ti) .or. (nb_types(i)%gti .eq. params(j)%tj)) ) then
                tia = .true.
                exit
            end if
        end do
        tja = .false.
        do j=1,nparams
            if( .not. params(j)%enabled ) cycle
            if( (params(j)%realm .ne. REALM_VDW_EPS) .and. (params(j)%realm .ne. REALM_VDW_R0) ) cycle
            if( ((nb_types(i)%gtj .eq. params(j)%ti) .or. (nb_types(i)%gtj .eq. params(j)%tj)) ) then
                tja = .true.
                exit
            end if
        end do
        if( tia .and. tja ) then
            tote = tote + nb_types(i)%errval * nb_types(i)%num
            totw = totw + nb_types(i)%num
            write(DEV_OUT,30) i,types(nb_types(i)%gti)%name,types(nb_types(i)%gtj)%name,nb_types(i)%num, &
                              nb_types(i)%errval
        end if
    end do

    if( totw .gt. 0.0d0 ) then
        tote = tote / totw
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40) tote
    write(DEV_OUT,45) tote*NBPnlErrorWeight

 5 format('# NB Potential Penalties')
10 format('# ID TypA TypB   Num        Error')
20 format('# -- ---- ---- ------- ----------')
30 format(I4,1X,A4,1X,A4,1X,I7,1X,F10.5)
40 format('# Final penalty      = ',F10.5)
45 format('# Final penalty w/w  = ',F10.5)

end subroutine ffdev_err_nbpnl_summary

! ------------------------------------------------------------------------------

end module ffdev_err_nbpnl



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

module ffdev_err_qnb

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_qnb_init
! ==============================================================================

subroutine ffdev_err_qnb_init

    use ffdev_err_qnb_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableQNBError        = .false.
    PrintQNBErrorSummary  = .false.
    QNBErrorWeight        = 1.0

end subroutine ffdev_err_qnb_init

! ==============================================================================
! subroutine ffdev_err_qnb_error
! ==============================================================================

subroutine ffdev_err_qnb_error(error)

    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_qnb_dat
    use ffdev_nb2nb_dat
    use ffdev_targetset_dat
    use ffdev_nb2nb

    implicit none
    type(FFERROR_TYPE)  :: error
    ! --------------------------------------------
    integer             :: i
    real(DEVDP)         :: gnb,r0,eps,glj,tote,totw
    ! --------------------------------------------------------------------------

    error%qnb = 0.0d0

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nnb_types
        gnb = NB2NBTemp*DEV_Rgas*log(nb_types(i)%QNB)

        ! take LJ parameters from topology
        r0 = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%r0
        eps = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%eps

        glj = NB2NBTemp*DEV_Rgas*log(ffdev_nb2nb_calc_QLJ(r0,eps))

        tote = tote + nb_types(i)%num*(glj-gnb)**2
        totw = totw + nb_types(i)%num
    end do

    if( totw .gt. 0.0d0 ) then
        tote = sqrt(tote/totw)
    end if

    error%qnb = tote

end subroutine ffdev_err_qnb_error

! ==============================================================================
! subroutine ffdev_err_qnb_summary
! ==============================================================================

subroutine ffdev_err_qnb_summary()

    use ffdev_err_qnb_dat
    use ffdev_utils
    use ffdev_nb2nb_dat
    use ffdev_targetset_dat
    use ffdev_parameters_dat
    use ffdev_nb2nb

    implicit none
    integer         :: i
    real(DEVDP)     :: gnb,r0,eps,glj,qnb,qlj,diff,tote,totw
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    tote = 0.0d0
    totw = 0.0d0

    do i=1,nnb_types
        qnb = nb_types(i)%QNB
        gnb = NB2NBTemp*DEV_Rgas*log(qnb)
        ! take LJ parameters from topology
        r0 = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%r0
        eps = sets(nb_types(i)%setid)%top%nb_types(nb_types(i)%nbt)%eps
        qlj = ffdev_nb2nb_calc_QLJ(r0,eps)
        glj = NB2NBTemp*DEV_Rgas*log(qlj)
        diff = glj-gnb

        tote = tote + nb_types(i)%num*diff**2
        totw = totw + nb_types(i)%num

        write(DEV_OUT,30) i,types(nb_types(i)%gti)%name,types(nb_types(i)%gtj)%name,nb_types(i)%num, &
                          qnb,qlj,gnb,glj,diff
    end do

    if( totw .gt. 0.0d0 ) then
        tote = sqrt(tote/totw)
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40) tote
    write(DEV_OUT,45) tote*QNBErrorWeight

 5 format('# QNB Penalties')
10 format('# ID TypA TypB   Num        Q(NB)      Q(LJ)      G(NB)      G(LJ)    G(diff)')
20 format('# -- ---- ---- ------- ---------- ---------- ---------- ---------- ----------')
30 format(I4,1X,A4,1X,A4,1X,I7,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4)
40 format('# Final penalty               =                                    ',F10.4)
45 format('# Final penalty (all weights) =                                    ',F10.4)

end subroutine ffdev_err_qnb_summary

! ------------------------------------------------------------------------------

end module ffdev_err_qnb



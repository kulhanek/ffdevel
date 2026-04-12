! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2026 Petr Kulhanek, kulhanek@chemi.muni.cz
!    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_err_energy

use ffdev_constants
use ffdev_variables
use ffdev_errors_dat

!===============================================================================

type, extends(ErrorFceType) :: TypeEFTEnergy

    integer         :: EnergyErrorMode
    logical         :: EnableMaxFilter
    real(DEVDP)     :: MaxTargetEnergy
    logical         :: EnableMinFilter
    real(DEVDP)     :: MinTargetEnergy

    contains
        ! executive methods
        procedure   :: init_errfce              => ffdev_err_energy_init
        procedure   :: load_errfce              => ffdev_err_energy_ctrl
        procedure   :: set_title_errfce         => ffdev_err_energy_set_title
        procedure   :: calc_errfce              => ffdev_err_energy_error
        procedure   :: print_individual_summary_errfce => ffdev_err_energy_summary
end type TypeEFTEnergy

contains

! ==============================================================================
! subroutine ffdev_err_energy_init
! ==============================================================================

subroutine ffdev_err_energy_init(err_item)

    implicit none
    class(TypeEFTEnergy)    :: err_item
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%init_errfce()

    err_item%EnergyErrorMode    = EE_ABS
    err_item%EnableMaxFilter    = .false.
    err_item%MaxTargetEnergy    = 0.0
    err_item%EnableMinFilter    = .false.
    err_item%MinTargetEnergy    = 0.0

end subroutine ffdev_err_energy_init

! ==============================================================================
! subroutine ffdev_err_energy_ctrl
! ==============================================================================

subroutine ffdev_err_energy_ctrl(err_item,fin)

    use ffdev_errors_dat
    use ffdev_utils
    use ffdev_errors_utils
    use prmfile

    implicit none
    class(TypeEFTEnergy)        :: err_item
    type(PRMFILE_TYPE)          :: fin
    ! --------------------------------------------
    character(PRMFILE_MAX_PATH) :: string
    ! --------------------------------------------------------------------------

    call err_item%ErrorFceType%load_errfce(fin)

    if( prmfile_get_string_by_key(fin,'scale', string)) then
        err_item%EnergyErrorMode = ffdev_errors_utils_scale_from_string(string)
        write(DEV_OUT,140) ffdev_errors_utils_scale_to_string(err_item%EnergyErrorMode)
    else
        write(DEV_OUT,145) ffdev_errors_utils_scale_to_string(err_item%EnergyErrorMode)
    end if

    if( prmfile_get_logical_by_key(fin,'maxfilter', err_item%EnableMaxFilter)) then
        write(DEV_OUT,150) prmfile_onoff(err_item%EnableMaxFilter)
    else
        write(DEV_OUT,155) prmfile_onoff(err_item%EnableMaxFilter)
    end if
    if( prmfile_get_real8_by_key(fin,'maxvalue', err_item%MaxTargetEnergy)) then
        write(DEV_OUT,160) err_item%MaxTargetEnergy
    else
        write(DEV_OUT,165) err_item%MaxTargetEnergy
    end if

    if( prmfile_get_logical_by_key(fin,'minfilter', err_item%EnableMinFilter)) then
        write(DEV_OUT,170) prmfile_onoff(err_item%EnableMinFilter)
    else
        write(DEV_OUT,175) prmfile_onoff(err_item%EnableMinFilter)
    end if
    if( prmfile_get_real8_by_key(fin,'minvalue', err_item%MinTargetEnergy)) then
        write(DEV_OUT,180) err_item%MinTargetEnergy
    else
        write(DEV_OUT,185) err_item%MinTargetEnergy
    end if

140  format ('Error scale (scale)                    = ',a24)
145  format ('Error scale (scale)                    = ',a24,'      (default)')
150  format ('Enable max energy filter (maxfilter)   = ',a12)
155  format ('Enable max energy filter (maxfilter)   = ',a12,'                  (default)')
160  format ('Max target rnergy (maxvalue)           = ',f21.8)
165  format ('Max target rnergy (maxvalue)           = ',f21.8,'         (default)')
170  format ('Enable min energy filter (minfilter)   = ',a12)
175  format ('Enable min energy filter (minfilter)   = ',a12,'                  (default)')
180  format ('Min target rnergy (maxvalue)           = ',f21.8)
185  format ('Min target rnergy (maxvalue)           = ',f21.8,'         (default)')

end subroutine ffdev_err_energy_ctrl

!===============================================================================
! Subroutine:  ffdev_err_energy_set_title
!===============================================================================

subroutine ffdev_err_energy_set_title(err_item)

    implicit none
    class(TypeEFTEnergy)    :: err_item
    ! --------------------------------------------------------------------------

    err_item%Title = 'Energy'

end subroutine ffdev_err_energy_set_title

! ==============================================================================
! subroutine ffdev_err_energy_error
! ==============================================================================

subroutine ffdev_err_energy_error(err_item,opterr)

    use ffdev_targetset
    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_geometry
    use ffdev_errors_dat

    implicit none
    class(TypeEFTEnergy)    :: err_item
    logical                 :: opterr
    ! --------------------------------------------
    integer                 :: i,j,nene
    real(DEVDP)             :: err,seterrene,totw
    ! --------------------------------------------------------------------------

    err_item%RepValue = 0.0d0
    if( .not. err_item%Enabled ) then
        if( opterr ) return
    end if

    seterrene = 0.0
    nene = 0
    totw = 0

    do i=1,nsets
        ! use only sets, which can provide reliable energy
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            ! ------------------------------------------------------------------
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle

            ! filters
            if( err_item%EnableMaxFilter ) then
                if( sets(i)%geo(j)%trg_energy .gt. err_item%MaxTargetEnergy ) cycle
            end if
            if( err_item%EnableMinFilter ) then
                if( sets(i)%geo(j)%trg_energy .lt. err_item%MinTargetEnergy ) cycle
            end if

            select case(err_item%EnergyErrorMode)
                case(EE_ABS)
                    nene = nene + 1
                    err = sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy
                    seterrene = seterrene + sets(i)%geo(j)%weight * err**2
                    totw = totw + sets(i)%geo(j)%weight
                case(EE_REL)
                    if( sets(i)%geo(j)%trg_energy .gt. 0 ) then
                        nene = nene + 1
                        err = (sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy) / &
                              sets(i)%geo(j)%trg_energy
                        seterrene = seterrene + sets(i)%geo(j)%weight * err**2
                        totw = totw + sets(i)%geo(j)%weight
                    end if
                case(EE_LOG)
                    if( (sets(i)%geo(j)%total_ene .gt. 0) .and. &
                        (sets(i)%geo(j)%trg_energy .gt. 0) ) then
                        nene = nene + 1
                        err = log(sets(i)%geo(j)%total_ene) - log(sets(i)%geo(j)%trg_energy)
                        seterrene = seterrene + sets(i)%geo(j)%weight * err**2
                        totw = totw + sets(i)%geo(j)%weight
                    end if
            end select

        end do
    end do

    if( totw .gt. 0 ) then
        err_item%RepValue = sqrt(seterrene/totw)
    end if

end subroutine ffdev_err_energy_error

! ==============================================================================
! subroutine ffdev_err_energy_summary
! ==============================================================================

subroutine ffdev_err_energy_summary(err_item)

    use ffdev_targetset_dat
    use ffdev_geometry

    implicit none
    class(TypeEFTEnergy) :: err_item
    ! --------------------------------------------
    real(DEVDP)         :: aerr,aserr
    real(DEVDP)         :: rerr,rserr
    real(DEVDP)         :: lerr,lserr,maxerr,atotw,rtotw,ltotw
    integer             :: i,j
    logical             :: printsum
    ! --------------------------------------------------------------------------

    if( .not. err_item%PrintSummary ) return

    printsum = .false.
    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( sets(i)%geo(j)%trg_ene_loaded .and. (sets(i)%isref .eqv. .false.) ) then
                printsum = .true.
            end if
        end do
    end do
    if( .not. printsum ) return

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10,ADVANCE='NO')
    write(DEV_OUT,50)
    write(DEV_OUT,20,ADVANCE='NO')
    write(DEV_OUT,60)

    aserr = 0.0d0
    rserr = 0.0d0
    lserr = 0.0d0
    atotw = 0.0d0
    rtotw = 0.0d0
    ltotw = 0.0d0

    maxerr = 0.0d0

    do i=1,nsets
        ! use only sets, which can provide reliable energy
        if( .not. ( (sets(i)%nrefs .ge. 1) .or. (sets(i)%top%probe_size .gt. 0) ) ) cycle

        do j=1,sets(i)%ngeos
            printsum = .false.
            if( .not. sets(i)%geo(j)%trg_ene_loaded ) cycle

            printsum = .true.

            ! absolute
            aerr  = sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy
            aserr = aserr + sets(i)%geo(j)%weight * aerr**2
            atotw  = atotw + sets(i)%geo(j)%weight

            if( abs(aerr) .gt. abs(maxerr) ) then
                maxerr = aerr
            end if

            ! relative
            rerr  = 0.0
            if( sets(i)%geo(j)%trg_energy .ne. 0 ) then
                rerr  = (sets(i)%geo(j)%total_ene - sets(i)%geo(j)%trg_energy)/sets(i)%geo(j)%trg_energy
                rserr = rserr + sets(i)%geo(j)%weight * rerr**2
                rtotw  = rtotw + sets(i)%geo(j)%weight
            end if
            ! log
            lerr  = 0.0
            if( (sets(i)%geo(j)%trg_energy .gt. 0) .and. (sets(i)%geo(j)%total_ene .gt. 0) ) then
                lerr  = log(sets(i)%geo(j)%total_ene) - log(sets(i)%geo(j)%trg_energy)
                lserr = lserr + sets(i)%geo(j)%weight * lerr**2
                ltotw  = ltotw + sets(i)%geo(j)%weight
            end if

            write(DEV_OUT,30,ADVANCE='NO') i, j, sets(i)%geo(j)%weight, &
                              sets(i)%geo(j)%total_ene, sets(i)%geo(j)%trg_energy, aerr, rerr*100.0d0, lerr

            write(DEV_OUT,70) sets(i)%geo(j)%bond_ene, sets(i)%geo(j)%angle_ene, sets(i)%geo(j)%dih_ene, &
                              sets(i)%geo(j)%impropr_ene, sets(i)%geo(j)%dih_ene + sets(i)%geo(j)%impropr_ene, &
                              sets(i)%geo(j)%ele_ene, sets(i)%geo(j)%pen_ene, sets(i)%geo(j)%ele14_ene, &
                              sets(i)%geo(j)%ele_ene + sets(i)%geo(j)%pen_ene + sets(i)%geo(j)%ele14_ene, &
                              sets(i)%geo(j)%ind_ene, &
                              sets(i)%geo(j)%rep_ene, sets(i)%geo(j)%rep14_ene, sets(i)%geo(j)%rep_ene + sets(i)%geo(j)%rep14_ene, &
                              sets(i)%geo(j)%dis_ene, sets(i)%geo(j)%dis14_ene, sets(i)%geo(j)%dis_ene + sets(i)%geo(j)%dis14_ene, &
                              sets(i)%geo(j)%bond_ene + sets(i)%geo(j)%angle_ene      &
                               + sets(i)%geo(j)%dih_ene + sets(i)%geo(j)%impropr_ene, &
                              sets(i)%geo(j)%ele_ene + sets(i)%geo(j)%pen_ene + sets(i)%geo(j)%ele14_ene  &
                               + sets(i)%geo(j)%ind_ene &
                               + sets(i)%geo(j)%rep_ene + sets(i)%geo(j)%rep14_ene  &
                               + sets(i)%geo(j)%dis_ene + sets(i)%geo(j)%dis14_ene

            if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
                call ffdev_geometry_info_ene(sets(i)%geo(j))
                write(DEV_OUT,*)
            end if

        end do
        if( printsum ) then
            write(DEV_OUT,20,ADVANCE='NO')
            write(DEV_OUT,60)
        end if
    end do

    if( atotw .gt. 0 ) then
        aserr = sqrt(aserr / atotw)
    end if
    if( rtotw .gt. 0 ) then
        rserr = sqrt(rserr / rtotw)
    end if
    if( ltotw .gt. 0 ) then
        lserr = sqrt(lserr / ltotw)
    end if

    write(DEV_OUT,35)  maxerr
    write(DEV_OUT,40)  aserr, rserr*100.0d0, lserr
    write(DEV_OUT,45)  err_item%Weight*aserr, err_item%Weight*rserr*100.0d0, err_item%Weight*lserr

 5 format('# Energy errors')
10 format('# SET GeoID Weight      E(MM)     E(TGR)     Err(E) relErr%(E)  logErr(E) | ')
20 format('# --- ----- ------ ---------- ---------- ---------- ---------- ---------- | ')
30 format(I5,1X,I5,1X,F6.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,3X)

35 format('# Maximum signed error (MSE)          =  ',F10.3)
40 format('# Root mean square error (RMSE)       =  ',F10.3,1X,F10.3,1X,F10.3)
45 format('# Final RMSE (all weights)            =  ',F10.3,1X,F10.3,1X,F10.3)

50 format('     Ebonds    Eangles      Etors      Eimps  Edih(t+i)        Eel       Epen      E14el    Etotele' &
          '       Eind       Erep     E14rep    Etotrep      Edisp    E14disp   Etotdisp        Ebn        Enb')
60 format(' ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------' &
          ' ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------')
70 format(1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3, &
          1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3,1X,F10.3)

end subroutine ffdev_err_energy_summary

! ------------------------------------------------------------------------------

end module ffdev_err_energy



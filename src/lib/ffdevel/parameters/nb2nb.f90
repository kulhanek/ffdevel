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

module ffdev_nb2nb

use ffdev_nb2nb_dat
use ffdev_constants
use ffdev_variables

implicit none

interface

! shark interfaces
    subroutine shark_create2(nactparms,initial_step,initial_params)
        integer(4)      :: nactparms
        real(8)         :: initial_step
        real(8)         :: initial_params(*)
    end subroutine shark_create2

    subroutine shark_dostep2(error)
        real(8)         :: error
    end subroutine shark_dostep2

    subroutine shark_getsol2(params)
        real(8)         :: params(*)
    end subroutine shark_getsol2

    subroutine shark_destroy2()
    end subroutine shark_destroy2

    subroutine shark_set_rngseed2(seed)
        integer(4)      :: seed
    end subroutine shark_set_rngseed2

end interface

contains

! ==============================================================================
! subroutine ffdev_nb2nb_nb2lj_mode_to_string
! ==============================================================================

character(80) function ffdev_nb2nb_nb2lj_mode_to_string(mode)

    use ffdev_nb2nb_dat
    use ffdev_utils

    implicit none
    integer  :: mode
    ! --------------------------------------------------------------------------

    select case(mode)
        case(NB2LJ_MODE_MINIMUM)
            ffdev_nb2nb_nb2lj_mode_to_string = 'MINIMUM'
        case(NB2LJ_MODE_OVERLAY)
            ffdev_nb2nb_nb2lj_mode_to_string = 'OVERLAY'
        case(NB2LJ_MODE_OVERLAY_WEIGHTED)
            ffdev_nb2nb_nb2lj_mode_to_string = 'OVERLAY-WEIGHTED'
        case(NB2LJ_MODE_OVERLAY_REP)
            ffdev_nb2nb_nb2lj_mode_to_string = 'OVERLAY-REPULSION'
        case(NB2LJ_MODE_OVERLAY_DISP)
            ffdev_nb2nb_nb2lj_mode_to_string = 'OVERLAY-DISPERSION'
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented in ffdev_nb2nb_nb2lj_mode_to_string!')
    end select

end function ffdev_nb2nb_nb2lj_mode_to_string

! ==============================================================================
! subroutine ffdev_nb2nb_nb2lj_mode_from_string
! ==============================================================================

integer function ffdev_nb2nb_nb2lj_mode_from_string(string)

    use ffdev_nb2nb_dat
    use ffdev_utils

    implicit none
    character(*)   :: string
    ! --------------------------------------------------------------------------

    select case(trim(string))
        case('MINIMUM')
            ffdev_nb2nb_nb2lj_mode_from_string = NB2LJ_MODE_MINIMUM
        case('OVERLAY')
            ffdev_nb2nb_nb2lj_mode_from_string = NB2LJ_MODE_OVERLAY
        case('OVERLAY-WEIGHTED')
            ffdev_nb2nb_nb2lj_mode_from_string = NB2LJ_MODE_OVERLAY_WEIGHTED
        case('OVERLAY-REPULSION')
            ffdev_nb2nb_nb2lj_mode_from_string = NB2LJ_MODE_OVERLAY_REP
        case('OVERLAY-DISPERSION')
            ffdev_nb2nb_nb2lj_mode_from_string = NB2LJ_MODE_OVERLAY_DISP
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Not implemented "' // trim(string) //'" in ffdev_nb2nb_nb2lj_mode_from_string!')
    end select

end function ffdev_nb2nb_nb2lj_mode_from_string

! ==============================================================================
! function ffdev_nb2nb_init_nbtypes
! ==============================================================================

subroutine ffdev_nb2nb_init_nbtypes

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat
    use ffdev_utils

    implicit none
    ! --------------------------------------------
    integer         :: i,j,idx,alloc_status,estat
    ! --------------------------------------------------------------------------

    if( ApplyCombiningRules ) then
        nnb_types = ntypes
    else
        nnb_types = ntypes * (ntypes - 1) / 2 + ntypes
    end if

    ! allocate working array
    allocate(nb_types(nnb_types),stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data for nb_types in ffdev_nb2nb_init_nbtypes!')
    end if

    idx = 1
    do i=1,ntypes
        do j=i,ntypes
            if( ApplyCombiningRules ) then
                if( i .ne. j ) cycle
            end if
            nb_types(idx)%gti   = i
            nb_types(idx)%gtj   = j
            nb_types(idx)%eps   = 0
            nb_types(idx)%r0    = 0
            nb_types(idx)%PA    = 0
            nb_types(idx)%PB    = 0
            nb_types(idx)%RC    = 0
            nb_types(idx)%TB    = 0
            nb_types(idx)%setid = 0
            nb_types(idx)%nbt   = 0
            idx = idx + 1
        end do
    end do

    call execute_command_line('mkdir -p ' // trim(NBPotPath), exitstat = estat )
    if( estat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to create NBPotPath!')
    end if

end subroutine ffdev_nb2nb_init_nbtypes

! ==============================================================================
! function ffdev_nb2nb_gather_nbtypes
! ==============================================================================

subroutine ffdev_nb2nb_gather_nbtypes

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    ! --------------------------------------------
    integer         :: i,j,k
    integer         :: gti,gtj
    ! --------------------------------------------------------------------------

    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            gti = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%ti)%glbtypeid
            gtj = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%tj)%glbtypeid

            do k=1,nnb_types
                if( nb_types(k)%setid .ne. 0 ) cycle
                if( ( (nb_types(k)%gti .eq. gti) .and. (nb_types(k)%gtj .eq. gtj) ) .or. &
                    ( (nb_types(k)%gti .eq. gtj) .and. (nb_types(k)%gtj .eq. gti) ) ) then
                    nb_types(k)%setid = i
                    nb_types(k)%nbt   = j
                    nb_types(k)%r0  = sets(i)%top%nb_types(j)%r0
                    nb_types(k)%eps = sets(i)%top%nb_types(j)%eps
                    nb_types(k)%PA  = sets(i)%top%nb_types(j)%PA
                    nb_types(k)%PB  = sets(i)%top%nb_types(j)%PB
                    nb_types(k)%TB  = sets(i)%top%nb_types(j)%TB
                    nb_types(k)%RC  = sets(i)%top%nb_types(j)%RC
                    exit
                end if
            end do
        end do
    end do

end subroutine ffdev_nb2nb_gather_nbtypes

! ==============================================================================
! function ffdev_nb2nb_scatter_nbtypes
! ==============================================================================

subroutine ffdev_nb2nb_scatter_nbtypes

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    ! --------------------------------------------
    integer         :: i,j,k
    integer         :: gti,gtj
    ! --------------------------------------------------------------------------

    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            gti = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%ti)%glbtypeid
            gtj = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%tj)%glbtypeid

            do k=1,nnb_types
                if( ( (nb_types(k)%gti .eq. gti) .and. (nb_types(k)%gtj .eq. gtj) ) .or. &
                    ( (nb_types(k)%gti .eq. gtj) .and. (nb_types(k)%gtj .eq. gti) ) ) then
                    sets(i)%top%nb_types(j)%r0  = nb_types(k)%r0
                    sets(i)%top%nb_types(j)%eps = nb_types(k)%eps
                    sets(i)%top%nb_types(j)%PA  = nb_types(k)%PA
                    sets(i)%top%nb_types(j)%PB  = nb_types(k)%PB
                    sets(i)%top%nb_types(j)%TB  = nb_types(k)%TB
                    sets(i)%top%nb_types(j)%RC  = nb_types(k)%RC
                    exit
                end if
            end do
            sets(i)%top%nb_params_update = .true.
        end do
    end do

end subroutine ffdev_nb2nb_scatter_nbtypes

! ==============================================================================
! function ffdev_nb2nb_switch_nbmode
! ==============================================================================

subroutine ffdev_nb2nb_switch_nbmode(from_nb_mode,to_nb_mode)

    use ffdev_topology_dat
    use ffdev_utils

    implicit none
    integer         :: from_nb_mode,to_nb_mode
    ! --------------------------------------------
    integer         :: nbt
    real(DEVDP)     :: pa
    ! --------------------------------------------------------------------------

    if( from_nb_mode .ne. nb_mode ) then
        call ffdev_utils_exit(DEV_ERR,1,'from_nb_mode .ne. nb_mode in ffdev_nb2nb_switch_nbmode!')
    end if

    select case(from_nb_mode)
    !---------------------------------------------
        case(NB_VDW_LJ)
            select case(to_nb_mode)
                case(NB_VDW_LJ)
                    ! nothing to do

                case(NB_VDW_EXP_DISPTT,NB_VDW_12_DISPBJ,NB_VDW_EXP_DISPBJ)
                    do nbt=1,nnb_types
                        pa = 6.0d0 * nb_types(nbt)%eps * exp(lj2exp6_alpha)/(lj2exp6_alpha - 6.0d0)
                        if( pa .gt. 0 ) then
                            pa = log(pa)
                        else
                            pa = 1.0d0
                        end if
                        nb_types(nbt)%pa = pa
                        if( nb_types(nbt)%r0 .ne. 0 ) then
                            nb_types(nbt)%pb = lj2exp6_alpha / nb_types(nbt)%r0
                        else
                            nb_types(nbt)%pb = ljdefpb
                        end if
                        nb_types(nbt)%rc = ljdefrc
                    end do
                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode(to) in ffdev_nb2nb_switch_nbmode!')
            end select
    !---------------------------------------------
        case(NB_VDW_EXP_DISPTT,NB_VDW_12_DISPBJ,NB_VDW_EXP_DISPBJ)
            select case(to_nb_mode)
                case(NB_VDW_LJ)
                    do nbt=1,nnb_types
                        call ffdev_nb2nb_nb2lj(nbt)
                    end do

                case(NB_VDW_EXP_DISPTT,NB_VDW_12_DISPBJ,NB_VDW_EXP_DISPBJ)
                    ! nothing to do

                case default
                    call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode(to) in ffdev_nb2nb_switch_nbmode!')
            end select
    !---------------------------------------------
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode(from) in ffdev_nb2nb_switch_nbmode!')
    end select

end subroutine ffdev_nb2nb_switch_nbmode

! ==============================================================================
! subroutine ffdev_nb2nb_conv_sum
! ==============================================================================

subroutine ffdev_nb2nb_conv_sum

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat

    implicit none
    integer                 :: gnbt
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,10) trim(ffdev_nb2nb_nb2lj_mode_to_string(NB2LJMode))

    write(DEV_OUT,*)
    write(DEV_OUT,20)
    write(DEV_OUT,30)

    do gnbt = 1, nnb_types
        write(DEV_OUT,40) trim(types(nb_types(gnbt)%gti)%name), trim(types(nb_types(gnbt)%gtj)%name), &
                          nb_types(gnbt)%eps,nb_types(gnbt)%r0,nb_types(gnbt)%errval
    end do

    write(DEV_OUT,*)

10 format('# NB->LJ conversion method: ',A)

20 format('# TypA TypB          Eps           R0          Err')
30 format('# ---- ---- ------------ ------------ ------------')
40 format(2X,A4,1X,A4,1X,F12.6,1X,F12.6,1X,F12.6)

end subroutine ffdev_nb2nb_conv_sum

! ==============================================================================
! subroutine ffdev_nb2nb_nb2lj
! ==============================================================================

subroutine ffdev_nb2nb_nb2lj(gnbt)

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    integer                 :: gnbt
    ! --------------------------------------------
    real(DEVDP)             :: r0,eps,sig,errval
    integer                 :: alloc_status,istep,nbt
    type(TOPOLOGY),target   :: top
    ! --------------------------------------------------------------------------

    top =  sets(nb_types(gnbt)%setid)%top
    nbt =  nb_types(gnbt)%nbt

    ! setup potential
    NB2LJNBPair%ai       = 0
    NB2LJNBPair%aj       = 0
    NB2LJNBPair%nbt      = nbt
    NB2LJNBPair%dt       = 0
    NB2LJNBPair%nbtii    = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(nbt)%ti,top%nb_types(nbt)%ti)
    NB2LJNBPair%nbtjj    = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(nbt)%tj,top%nb_types(nbt)%tj)
    call ffdev_topology_update_nbpair_prms(top,NB2LJNBPair)

    if( Verbosity .ge. DEV_VERBOSITY_MEDIUM ) then
        write(DEV_OUT,10) trim(types(nb_types(gnbt)%gti)%name), trim(types(nb_types(gnbt)%gtj)%name)
    end if

    ! find r0, eps, sigma
    call ffdev_nb2nb_find_min_for_nbpair(r0,eps)
    call ffdev_nb2nb_find_sig_for_nbpair(sig)

    NB2LJSigma = sig

    ! write the source potential
    call ffdev_nb2nb_write_source_pot(gnbt)

    ! write source potential parameters
    call ffdev_nb2nb_write_source_prms(gnbt,r0,eps,sig)

    select case(NB2LJMode)
        case(NB2LJ_MODE_MINIMUM)
            nb_types(gnbt)%r0  = r0
            nb_types(gnbt)%eps = eps

            call ffdev_nb2nb_write_LJ(gnbt,r0,eps)
            return
        case(NB2LJ_MODE_OVERLAY,NB2LJ_MODE_OVERLAY_WEIGHTED,NB2LJ_MODE_OVERLAY_REP,NB2LJ_MODE_OVERLAY_DISP)
            ! nothing to be done here
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported NB2LJMode in ffdev_nb2nb_nb2nb!')
    end select

    ! minimize overlay
    ! --------------------------------------------------------------------------

    NB2LJNParams = 2

    ! allocate working array
    allocate(NB2LJtmp_lb(NB2LJNParams),NB2LJtmp_ub(NB2LJNParams),NB2LJprms(NB2LJNParams),&
             NB2LJtmp_xg(NB2LJNParams),stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data for SHARK optimization in ffdev_nb2nb_nb2nb!')
    end if

    select case(NB2LJMode)
        case(NB2LJ_MODE_MINIMUM)
            call ffdev_utils_exit(DEV_ERR,1,'Illegal exec path in ffdev_nb2nb_nb2nb!')
        case(NB2LJ_MODE_OVERLAY,NB2LJ_MODE_OVERLAY_WEIGHTED)
            NB2LJMinR   = NB2LJSigma
            NB2LJMaxR   = NB2LJCutoffR
        case(NB2LJ_MODE_OVERLAY_REP)
            NB2LJMinR   = NB2LJSigma
            NB2LJMaxR   = r0
        case(NB2LJ_MODE_OVERLAY_DISP)
            NB2LJMinR   = r0
            NB2LJMaxR   = NB2LJCutoffR
        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported NB2LJMode in ffdev_nb2nb_nb2nb!')
    end select

    NB2LJtmp_lb(1) = sig
    NB2LJtmp_lb(2) = 0.0d0

    NB2LJtmp_ub(1) = NB2LJMaxR
    NB2LJtmp_ub(2) = 2.0d0*eps

    NB2LJprms(1) = r0
    NB2LJprms(2) = eps

    ! box data
    call ffdev_nb2nb_nb2nb_box_prms(NB2LJprms,NB2LJtmp_xg)

    ! init optimizer
    call shark_create2(NB2LJNParams,NB2LJSharkInitialStep,NB2LJtmp_xg)

    if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
        write(DEV_OUT,*)
        write(DEV_OUT,20)
        write(DEV_OUT,30)
    end if

    NB2LJErrFceEval = 0
    do istep = 1, NB2LJIterOpt
        call shark_dostep2(errval)
        if( Verbosity .ge. DEV_VERBOSITY_FULL ) then
            write(DEV_OUT,40) istep,NB2LJErrFceEval,errval
        end if
        nb_types(gnbt)%errval = errval
    end do

    ! get solution
    call shark_getsol2(NB2LJtmp_xg)

    ! unbox data
    call ffdev_nb2nb_nb2nb_unbox_prms(NB2LJtmp_xg,NB2LJprms)

    nb_types(gnbt)%r0  = NB2LJprms(1)
    nb_types(gnbt)%eps = NB2LJprms(2)

    call ffdev_nb2nb_write_LJ(gnbt,NB2LJprms(1),NB2LJprms(2))

    ! destroy optimizer
    call shark_destroy2()

    deallocate(NB2LJtmp_lb)
    deallocate(NB2LJtmp_ub)
    deallocate(NB2LJprms)
    deallocate(NB2LJtmp_xg)

10 format('Converting NB2LJ parameters for:',1X,2A,1X,2A)

20 format('#     Step FceEvals          Error')
30 format('# -------- -------- --------------')
40 format(I10,1X,I8,1X,E14.6)

end subroutine ffdev_nb2nb_nb2lj

! ==============================================================================
! subroutine ffdev_ffopt_unbox_prms
! ==============================================================================

subroutine ffdev_nb2nb_nb2nb_unbox_prms(boxed,prms)

    implicit none
    real(DEVDP)     :: boxed(:)
    real(DEVDP)     :: prms(:)
    ! --------------------------------------------
    integer         :: i
    ! --------------------------------------------------------------------------

    do i=1,NB2LJNParams
        prms(i) = NB2LJtmp_lb(i) + 0.5d0*(NB2LJtmp_ub(i)-NB2LJtmp_lb(i))*(1.0d0-cos(boxed(i)))
    end do

end subroutine ffdev_nb2nb_nb2nb_unbox_prms

! ==============================================================================
! subroutine ffdev_nb2nb_nb2nb_box_prms
! ==============================================================================

subroutine ffdev_nb2nb_nb2nb_box_prms(prms,boxed)

    implicit none
    real(DEVDP)     :: prms(:)
    real(DEVDP)     :: boxed(:)
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: sc
    ! --------------------------------------------------------------------------

    do i=1,NB2LJNParams
        sc = 1.0d0 - 2.0d0*(prms(i)-NB2LJtmp_lb(i))/(NB2LJtmp_ub(i)-NB2LJtmp_lb(i))
        if( sc .gt.  1.0d0 ) sc =  1.0d0
        if( sc .lt. -1.0d0 ) sc = -1.0d0
        boxed(i) = acos(sc)
    end do

end subroutine ffdev_nb2nb_nb2nb_box_prms

!===============================================================================
! subroutine ffdev_nb2nb_write_source_pot
!===============================================================================

subroutine ffdev_nb2nb_write_source_pot(gnbt)

    use ffdev_nb2nb_dat
    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer                 :: gnbt
    ! --------------------------------------------
    integer                 :: i
    real(DEVDP)             :: r,dr,enb
    character(len=MAX_PATH) :: name
    ! --------------------------------------------------------------------------

    name = trim(NBPotPath) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.nb'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    r  = NB2LJSigma
    dr = (NB2LJCutoffR - NB2LJSigma)/real(NB2LJIterErr,DEVDP)

    do i=1,NB2LJIterErr
        enb = ffdev_nb2nb_nbpair(NB2LJNBPair,r)
        write(DEV_NBPOT,10) r,enb
        r = r + dr
    end do

    close(DEV_NBPOT)

    10 format(F10.6,1X,E20.12)

end subroutine ffdev_nb2nb_write_source_pot

!===============================================================================
! subroutine ffdev_nb2nb_write_LJ
!===============================================================================

subroutine ffdev_nb2nb_write_LJ(gnbt,r0,eps)

    use ffdev_nb2nb_dat
    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer                 :: gnbt
    real(DEVDP)             :: r0,eps
    ! --------------------------------------------
    integer                 :: i
    real(DEVDP)             :: r,dr,enb
    character(len=MAX_PATH) :: name
    ! --------------------------------------------------------------------------

    name = trim(NBPotPath) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.lj'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    r  = NB2LJSigma
    dr = (NB2LJCutoffR - NB2LJSigma)/real(NB2LJIterErr,DEVDP)

    do i=1,NB2LJIterErr
        enb = ffdev_nb2nb_ljene(r0,eps,r)
        write(DEV_NBPOT,10) r,enb
        r = r + dr
    end do

    close(DEV_NBPOT)

    10 format(F10.6,1X,E20.12)

end subroutine ffdev_nb2nb_write_LJ

!===============================================================================
! subroutine ffdev_nb2nb_write_source_prms
!===============================================================================

subroutine ffdev_nb2nb_write_source_prms(gnbt,r0,eps,sig)

    use ffdev_nb2nb_dat
    use ffdev_utils
    use ffdev_parameters_dat

    implicit none
    integer                 :: gnbt
    real(DEVDP)             :: r0,eps,sig
    ! --------------------------------------------
    character(len=MAX_PATH) :: name
    ! --------------------------------------------------------------------------

    name = trim(NBPotPath) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.prms'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    write(DEV_NBPOT,10) r0,-eps
    write(DEV_NBPOT,10) sig,0.0d0

    close(DEV_NBPOT)

    10 format(F10.6,1X,E20.12)

end subroutine ffdev_nb2nb_write_source_prms

! ==============================================================================
! subroutine ffdev_nb2nb_find_min_for_nbpair
! ==============================================================================

subroutine ffdev_nb2nb_find_min_for_nbpair(r0,eps)

    use ffdev_nb2nb_dat

    implicit none
    real(DEVDP)     :: r0,eps
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r1,r2,r3,r4,gr,f2,f3
    ! --------------------------------------------------------------------------

    ! Golden-section search
    ! https://en.wikipedia.org/wiki/Golden-section_search

    r1 = 1.0d0
    r4 = NB2LJCutoffR
    gr = (1.0d0 + sqrt(5.0d0)) * 0.5d0

    do i=1,NB2LJIterGS
        r2 = r4 - (r4 - r1) / gr
        r3 = r1 + (r4 - r1) / gr
        f2 = ffdev_nb2nb_nbpair(NB2LJNBPair,r2)
        f3 = ffdev_nb2nb_nbpair(NB2LJNBPair,r3)
        if( f2 .lt. f3 ) then
            r4 = r3
        else
            r1 = r2
        end if
    end do

    r0  = (r1+r4)*0.5d0
    eps = - ffdev_nb2nb_nbpair(NB2LJNBPair,r0)

end subroutine ffdev_nb2nb_find_min_for_nbpair

! ==============================================================================
! subroutine ffdev_nb2nb_find_sig_for_nbpair
! ==============================================================================

subroutine ffdev_nb2nb_find_sig_for_nbpair(sig)

    use ffdev_nb2nb_dat

    implicit none
    real(DEVDP)     :: sig
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r1,r2,f1,f2,rm,fm
    ! --------------------------------------------------------------------------

    ! Bisection method
    ! https://en.wikipedia.org/wiki/Bisection_method

    r1 = 1.0d0
    r2 = NB2LJCutoffR

    f1 = ffdev_nb2nb_nbpair(NB2LJNBPair,r1)
    f2 = ffdev_nb2nb_nbpair(NB2LJNBPair,r2)

    do i=1,NB2LJIterBS
        rm = (r1 + r2)*0.5d0
        fm = ffdev_nb2nb_nbpair(NB2LJNBPair,rm)
        if( fm .gt. 0 ) then
            f1 = fm
            r1 = rm
        else
            f2 = fm
            r2 = rm
        end if
    end do

    sig = rm

end subroutine ffdev_nb2nb_find_sig_for_nbpair

! ==============================================================================
! function ffdev_nb2nb_ljene
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_ljene(r0,eps,r)

    implicit none
    real(DEVDP)     :: r0,eps
    real(DEVDP)     :: r
    ! --------------------------------------------------------------------------

    ffdev_nb2nb_ljene = eps*( (r0/r)**12 - 2.0d0*(r0/r)**6 )

end function ffdev_nb2nb_ljene

! ==============================================================================
! function ffdev_nb2nb_nbpair
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_nbpair(nbpair,r)

    use ffdev_nb2nb_dat
    use ffdev_utils

    use ffdev_nbmode_LJ
    use ffdev_nbmode_12_DISPBJ
    use ffdev_nbmode_EXP_DISPBJ
    use ffdev_nbmode_EXP_DISPTT

    implicit none
    type(NB_PAIR)   :: nbpair
    real(DEVDP)     :: r
    ! --------------------------------------------------------------------------

    ffdev_nb2nb_nbpair = 0.0d0

    select case(nb_mode)
        case(NB_VDW_LJ)
            ffdev_nb2nb_nbpair = ffdev_energy_nbpair_LJ(nbpair,r)

        case(NB_VDW_12_DISPBJ)
            ffdev_nb2nb_nbpair = ffdev_energy_nbpair_12_DISPBJ(nbpair,r)

        case(NB_VDW_EXP_DISPBJ)
            ffdev_nb2nb_nbpair = ffdev_energy_nbpair_EXP_DISPBJ(nbpair,r)

        case(NB_VDW_EXP_DISPTT)
            ffdev_nb2nb_nbpair = ffdev_energy_nbpair_EXP_DISPTT(nbpair,r)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_nb2nb_nbpair!')
    end select

end function ffdev_nb2nb_nbpair

! ------------------------------------------------------------------------------

end module ffdev_nb2nb

!===============================================================================
! subroutine opt_shark_fce2
!===============================================================================

subroutine opt_shark_fce2(n, x, errval)

    use ffdev_nb2nb_dat
    use ffdev_nb2nb

    implicit none
    integer         :: n
    real(DEVDP)     :: x(n)
    real(DEVDP)     :: errval
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r,dr,enb,elj,r0,eps,w,sumw
    ! --------------------------------------------------------------------------

    call ffdev_nb2nb_nb2nb_unbox_prms(x,NB2LJprms)

    r  = NB2LJMinR
    dr = (NB2LJMaxR - NB2LJMinR)/real(NB2LJIterErr,DEVDP)
    errval = 0.0d0

    r0  = NB2LJprms(1)
    eps = NB2LJprms(2)

    sumw = 0.0
    do i=1,NB2LJIterErr
        enb = ffdev_nb2nb_nbpair(NB2LJNBPair,r)
        elj = ffdev_nb2nb_ljene(r0,eps,r)

        if( NB2LJ_MODE_OVERLAY_WEIGHTED .eq. NB2LJ_MODE_OVERLAY_WEIGHTED ) then
            ! eps is positive, but minimum is at -eps
            w = exp(-(elj+eps)/(DEV_Rgas*NB2LJTemp))
        else
            w = 1.0
        end if

        errval = errval + w*(enb-elj)**2
        sumw = sumw + w
        r = r + dr
    end do

    errval = sqrt(errval/sumw)

    NB2LJErrFceEval = NB2LJErrFceEval + 1

end subroutine opt_shark_fce2

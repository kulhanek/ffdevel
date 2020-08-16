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
! function ffdev_nb2nb_initdirs
! ==============================================================================

subroutine ffdev_nb2nb_initdirs

    use ffdev_nb2nb_dat
    use ffdev_utils

    implicit none
    ! --------------------------------------------
    integer         :: estat
    ! --------------------------------------------------------------------------

    call execute_command_line('mkdir -p ' // trim(NBPotPathCore), exitstat = estat )
    if( estat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to create NBPotPathCore in ffdev_nb2nb_initdirs!')
    end if

end subroutine ffdev_nb2nb_initdirs

! ==============================================================================
! function ffdev_nb2nb_initdirs_for_prog
! ==============================================================================

subroutine ffdev_nb2nb_initdirs_for_prog

    use ffdev_nb2nb_dat
    use ffdev_utils

    implicit none
    integer                 :: estat
    ! --------------------------------------------------------------------------

    estat = 0

    write(NBPotPathPrg,20) trim(NBPotPathCore), CurrentProgID

    call execute_command_line('mkdir -p ' // trim(NBPotPathPrg), exitstat = estat )
    if( estat .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to create NBPotPathPrg in ffdev_nb2nb_initdirs_for_prog!')
    end if

 20 format(A,'/',I3.3,'-fin')

end subroutine ffdev_nb2nb_initdirs_for_prog

! ==============================================================================
! function ffdev_nb2nb_init_nbtypes
! ==============================================================================

subroutine ffdev_nb2nb_init_nbtypes

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat
    use ffdev_utils
    use ffdev_targetset_dat

    implicit none
    ! --------------------------------------------
    integer         :: i,j,idx,alloc_status,gti,gtj,k,l
    ! --------------------------------------------------------------------------

    if( ApplyCombiningRules .and. NB2NBLikeOnly ) then
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
            if( ApplyCombiningRules .and. NB2NBLikeOnly ) then
                if( i .ne. j ) cycle
            end if
            nb_types(idx)%gti   = i
            nb_types(idx)%gtj   = j
            nb_types(idx)%SigNB = 0
            nb_types(idx)%R0NB  = 0
            nb_types(idx)%EpsNB = 0
            nb_types(idx)%QNB   = 0
            nb_types(idx)%setid = 0
            nb_types(idx)%nbt   = 0
            nb_types(idx)%num   = 0
            allocate(nb_types(idx)%NBPot(NB2NBNBins),stat=alloc_status)
            if( alloc_status .ne. 0 ) then
                call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data for NBPot in ffdev_nb2nb_init_nbtypes!')
            end if
            idx = idx + 1
        end do
    end do

    ! count nb_types
    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            gti = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%ti)%glbtypeid
            gtj = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%tj)%glbtypeid

            do k=1,nnb_types
                if( ApplyCombiningRules ) then
                    ! any of the nbpair
                    if( ( (nb_types(k)%gti .eq. gti) .or. (nb_types(k)%gtj .eq. gtj) ) .or. &
                        ( (nb_types(k)%gti .eq. gtj) .or. (nb_types(k)%gtj .eq. gti) ) ) then
                        ! count all valid NB pairs
                        do l=1,sets(i)%top%nb_size
                            if( sets(i)%top%nb_list(l)%nbt .eq. j ) then
                                nb_types(k)%num = nb_types(k)%num + 1
                            end if
                        end do
                    end if
                else
                    ! explicit pair
                    if( ( (nb_types(k)%gti .eq. gti) .and. (nb_types(k)%gtj .eq. gtj) ) .or. &
                        ( (nb_types(k)%gti .eq. gtj) .and. (nb_types(k)%gtj .eq. gti) ) ) then
                        ! count all valid NB pairs
                        do l=1,sets(i)%top%nb_size
                            if( sets(i)%top%nb_list(l)%nbt .eq. j ) then
                                nb_types(k)%num = nb_types(k)%num + 1
                            end if
                        end do
                    end if
                end if
            end do
        end do
    end do

    call ffdev_nb2nb_setup_gaussian_quadrature

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

    ! reset setid
    do k=1,nnb_types
        nb_types(k)%setid = 0
    end do

    ! get data
    do i=1,nsets
        do j=1,sets(i)%top%nnb_types
            gti = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%ti)%glbtypeid
            gtj = sets(i)%top%atom_types(sets(i)%top%nb_types(j)%tj)%glbtypeid

            do k=1,nnb_types
                if( nb_types(k)%setid .ne. 0 ) cycle
                if( ((nb_types(k)%gti .eq. gti) .and. (nb_types(k)%gtj .eq. gtj)) .or. &
                    ((nb_types(k)%gti .eq. gtj) .and. (nb_types(k)%gtj .eq. gti)) ) then
                    nb_types(k)%setid = i
                    nb_types(k)%nbt   = j
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
                if( ((nb_types(k)%gti .eq. gti) .and. (nb_types(k)%gtj .eq. gtj)) .or. &
                    ((nb_types(k)%gti .eq. gtj) .and. (nb_types(k)%gtj .eq. gti)) ) then
                    sets(i)%top%nb_types(j)%r0  = nb_types(k)%R0NB
                    sets(i)%top%nb_types(j)%eps = nb_types(k)%EpsNB
                    exit
                end if
            end do
            sets(i)%top%nb_params_update = .true.
        end do
    end do

end subroutine ffdev_nb2nb_scatter_nbtypes

! ==============================================================================
! subroutine ffdev_nb2nb_write_all_current_pots
! ==============================================================================

subroutine ffdev_nb2nb_write_all_current_pots

    use ffdev_nb2nb_dat
    use ffdev_topology_utils
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    ! --------------------------------------------
    type(TOPOLOGY),target   :: top
    integer                 :: i,gnbt,nbt
    real(DEVDP)             :: sig,r0,eps,r
    ! --------------------------------------------------------------------------

    call ffdev_nb2nb_initdirs_for_prog

    write(DEV_OUT,*)
    write(DEV_OUT,10) trim(NBPotPathPrg)

    call ffdev_nb2nb_gather_nbtypes

    do gnbt = 1, nnb_types
        write(DEV_OUT,20) trim(types(nb_types(gnbt)%gti)%name) &
                          // '-' // trim(types(nb_types(gnbt)%gtj)%name)

        top =  sets(nb_types(gnbt)%setid)%top
        nbt =  nb_types(gnbt)%nbt

        ! setup potential
        NB2NBNBPair%ai       = 0
        NB2NBNBPair%aj       = 0
        NB2NBNBPair%nbt      = nbt
        NB2NBNBPair%dt       = 0
        NB2NBNBPair%nbtii    = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(nbt)%ti,top%nb_types(nbt)%ti)
        NB2NBNBPair%nbtjj    = ffdev_topology_find_nbtype_by_tindex(top,top%nb_types(nbt)%tj,top%nb_types(nbt)%tj)
        call ffdev_topology_update_nbpair_prms(top,NB2NBNBPair)

        ! required for print interval
        call ffdev_nb2nb_find_min_for_nbpair(r0,eps)
        call ffdev_nb2nb_find_sig_for_nbpair(sig)

        ! NB
        nb_types(gnbt)%SigNB    =  sig
        nb_types(gnbt)%R0NB     =  r0
        nb_types(gnbt)%EpsNB    =  eps
        nb_types(gnbt)%QNB      = ffdev_nb2nb_calc_QNB(sig)

        r = sig
        do i=1,NB2NBNBins
            nb_types(gnbt)%NBPot(i) = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r)
            r = r + NB2NBCutoffR / real(NB2NBNBins,DEVDP)
        end do

        ! write the source potential
        call ffdev_nb2nb_write_source_pot(gnbt)
    end do

 10 format('> Writing NB summary to: ',A)
 20 format('  > ',A)

end subroutine ffdev_nb2nb_write_all_current_pots

! ==============================================================================
! function ffdev_nb2nb_qnb_eps
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_qnb_eps(r0,qnb)

    use ffdev_nb2nb_dat
    use ffdev_parameters_dat
    use ffdev_targetset_dat
    use ffdev_utils

    implicit none
    real(DEVDP)             :: r0,qnb
    ! --------------------------------------------
    real(DEVDP)             :: errval
! DEBUG
!    real(DEVDP)             :: qlj
    integer                 :: alloc_status,istep
    ! --------------------------------------------------------------------------

    ffdev_nb2nb_qnb_eps = 0.0d0

    NB2NBNParams = 1

        ! allocate working array
    allocate(NB2NBtmp_lb(NB2NBNParams),NB2NBtmp_ub(NB2NBNParams),NB2NBprms(NB2NBNParams),&
             NB2NBtmp_xg(NB2NBNParams),stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data for SHARK optimization in ffdev_nb2nb_qnb_eps!')
    end if

    QNBModeEps  = .true.
    QNBR0       = r0
    QNBTrg      = qnb

    ! eps
    NB2NBtmp_lb(1)  = 0.0d0
    NB2NBtmp_ub(1)  = 2.0d0
    NB2NBprms(1)    = 0.10d0

    ! box data
    call ffdev_nb2nb_nb2nb_box_prms(NB2NBprms,NB2NBtmp_xg)

! DEBUG
!    write(*,*) r0, qnb

    ! init optimizer
    call shark_create3(NB2NBNParams,NB2NBSharkInitialStep,NB2NBtmp_xg)

    NB2NBErrFceEval = 0
    do istep = 1, NB2NBIterOpt
        call shark_dostep3(errval)
! DEBUG
!        write(*,*) istep,NB2NBErrFceEval,errval
    end do

    ! get solution
    call shark_getsol3(NB2NBtmp_xg)

    ! unbox data
    call ffdev_nb2nb_nb2nb_unbox_prms(NB2NBtmp_xg,NB2NBprms)
    ffdev_nb2nb_qnb_eps = NB2NBprms(1)    ! result

! DEBUG
!    qlj = ffdev_nb2nb_calc_QLJ(QNBR0,ffdev_nb2nb_qnb_eps)
!    write(*,*) QNBR0,ffdev_nb2nb_qnb_eps, QNBTrg, qlj

    ! destroy optimizer
    call shark_destroy3()

    deallocate(NB2NBtmp_lb)
    deallocate(NB2NBtmp_ub)
    deallocate(NB2NBprms)
    deallocate(NB2NBtmp_xg)

end function ffdev_nb2nb_qnb_eps

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

    do i=1,NB2NBNParams
        prms(i) = NB2NBtmp_lb(i) + 0.5d0*(NB2NBtmp_ub(i)-NB2NBtmp_lb(i))*(1.0d0-cos(boxed(i)))
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

    do i=1,NB2NBNParams
        sc = 1.0d0 - 2.0d0*(prms(i)-NB2NBtmp_lb(i))/(NB2NBtmp_ub(i)-NB2NBtmp_lb(i))
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
    real(DEVDP)             :: r,evdw
    character(len=MAX_PATH) :: name
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    name = trim(NBPotPathPrg) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.nb'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    r  = nb_types(gnbt)%SigNB

    do while(r .le. NB2NBCutoffR)
        call ffdev_nb2nb_nbpair(NB2NBNBPair,r,nbene)
        evdw = nbene%rep_ene + nbene%dis_ene
        if( NB2NBIncludePen ) then
            evdw = evdw + nbene%pen_ene
        end if
        if( NB2NBIncludeInd ) then
            evdw = evdw + nbene%ind_ene
        end if
        if( NB2NBCalcQNBIsoline ) then
            qnbeps = ffdev_nb2nb_qnb_eps(r,nb_types(gnbt)%QNB)
            ! put negative sign to compare qnbeps with NB or LJ potential
            write(DEV_NBPOT,10) r, evdw, nbene%tot_ene, nbene%ele_ene, nbene%pen_ene, nbene%ind_ene, &
                                nbene%rep_ene, nbene%dis_ene, -qnbeps
        else
            write(DEV_NBPOT,20) r, evdw, nbene%tot_ene, nbene%ele_ene, nbene%pen_ene, nbene%ind_ene, &
                                nbene%rep_ene, nbene%dis_ene
        end if
        r = r + NB2NBdrPrint
    end do

    close(DEV_NBPOT)

    ! write prms
    name = trim(NBPotPathPrg) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.nb-prms'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    write(DEV_NBPOT,20) nb_types(gnbt)%R0NB,-nb_types(gnbt)%EpsNB
    write(DEV_NBPOT,20) nb_types(gnbt)%SigNB,0.0d0

    close(DEV_NBPOT)


    10 format(F10.6,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12)
    20 format(F10.6,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12)

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
    real(DEVDP)             :: r,evdw,sig
    character(len=MAX_PATH) :: name
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    name = trim(NBPotPathPrg) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.lj'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    ! start from sigma
    sig = r0 / 2.0d0**(1.0d0/6.0d0)
    r  = sig

    do while(r .le. NB2NBCutoffR)
        call ffdev_nb2nb_ljene(r0,eps,r,nbene)
        evdw = nbene%rep_ene + nbene%dis_ene
        write(DEV_NBPOT,10) r, evdw, nbene%tot_ene, nbene%ele_ene, nbene%pen_ene, nbene%ind_ene, &
                                        nbene%rep_ene, nbene%dis_ene
        r = r + NB2NBdrPrint
    end do

    close(DEV_NBPOT)

    ! write prms
    name = trim(NBPotPathPrg) // '/' // trim(types(nb_types(gnbt)%gti)%name) &
           // '-' // trim(types(nb_types(gnbt)%gtj)%name) // '.lj-prms'

    ! open file
    call ffdev_utils_open(DEV_NBPOT,name,'U')

    write(DEV_NBPOT,20) r0,eps
    write(DEV_NBPOT,20) sig,0.0d0

    close(DEV_NBPOT)

    10 format(F10.6,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12)
    20 format(F10.6,1X,E20.12)

end subroutine ffdev_nb2nb_write_LJ

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

    r1 = 0.5d0
    r4 = NB2NBCutoffR
    gr = (1.0d0 + sqrt(5.0d0)) * 0.5d0

    do i=1,NB2NBIterGS
        r2 = r4 - (r4 - r1) / gr
        r3 = r1 + (r4 - r1) / gr
        f2 = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r2)
        f3 = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r3)
        if( f2 .lt. f3 ) then
            r4 = r3
        else
            r1 = r2
        end if
    end do

    r0  = (r1+r4)*0.5d0
    eps = - ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r0)

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

    r1 = 0.5d0
    r2 = NB2NBCutoffR

    f1 = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r1)
    f2 = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r2)

    do i=1,NB2NBIterBS
        rm = (r1 + r2)*0.5d0
        fm = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,rm)
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
! subroutine ffdev_nb2nb_setup_gaussian_quadrature
! ==============================================================================

subroutine ffdev_nb2nb_setup_gaussian_quadrature

    use ffdev_utils

    implicit none
    integer         :: i,n,m,j,alloc_status
    real(DEVDP)     :: z,p1,p2,p3,pp,z1,eps
    ! --------------------------------------------------------------------------

    n   = NB2NBGaussQuadOrder

    if( allocated(NB2NBGaussQuadA) ) then
        deallocate(NB2NBGaussQuadA,NB2NBGaussQuadW)
    end if

    allocate(NB2NBGaussQuadA(n),NB2NBGaussQuadW(n),stat=alloc_status)
    if( alloc_status .ne. 0 ) then
        call ffdev_utils_exit(DEV_ERR,1,'Unable to allocate data in ffdev_nb2nb_setup_gaussian_quadrature!')
    end if

    EPS = 3e-11

    m   = (n+1)/2

    do i=0, m-1
        z=cos(DEV_PI*(i+0.75d0)/(n+0.5d0));
        do while( abs(z-z1) > EPS)
            p1=1.0
            p2=0.0
            do j=1,n
                p3=p2
                p2=p1
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
            end do
            pp=n*(z*p1-p2)/(z*z-1.0)
            z1=z
            z=z1-p1/pp
        end do
        NB2NBGaussQuadA(i+1)  = -z
        NB2NBGaussQuadA(n-i)  = z
        NB2NBGaussQuadW(i+1)  = 2.0/((1.0-z*z)*pp*pp)
        NB2NBGaussQuadW(n-i)  = NB2NBGaussQuadW(i+1)
    end do

    ! write(*,*) NB2NBGaussQuadA, NB2NBGaussQuadW

end subroutine ffdev_nb2nb_setup_gaussian_quadrature

! ==============================================================================
! function ffdev_nb2nb_calc_QNB
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_calc_QNB(sig)

    implicit none
    real(DEVDP)     :: sig
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r,e,q,n,bi,en,dr
    ! --------------------------------------------------------------------------

    ffdev_nb2nb_calc_QNB = 0.0d0

    if( NB2NBUseGaussQuad ) then
        bi = sig
        en = NB2NBCutoffR + sig
        do i=1,NB2NBGaussQuadOrder
            r = 0.5d0*(en-bi)*NB2NBGaussQuadA(i) + 0.5d0*(en+bi)
            e = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r)
            q = 0.0d0
            if( e .lt. 0.0d0 ) then
                q = exp(-e/((DEV_Rgas*NB2NBTemp)))
            end if
            ffdev_nb2nb_calc_QNB = ffdev_nb2nb_calc_QNB + q*NB2NBGaussQuadW(i)
        end do

        ffdev_nb2nb_calc_QNB = 0.5d0*(en-bi)*ffdev_nb2nb_calc_QNB
        ffdev_nb2nb_calc_QNB = ffdev_nb2nb_calc_QNB/NB2NBCutoffR
    else
        r  = sig
        dr = NB2NBCutoffR/real(NB2NBNBins,DEVDP)
        q = 0.0d0
        n = 0.0d0
        do i=1,NB2NBNBins
            e = ffdev_nb2nb_nbpair_vdw_ene(NB2NBNBPair,r)
            if( e .lt. 0.0d0 ) then
                q = q + exp(-e/((DEV_Rgas*NB2NBTemp)))*dr
            end if
            r = r + dr
            n = n + dr
        end do

        ffdev_nb2nb_calc_QNB = 1.0d0
        if( n .gt. 0.0d0 ) then
            ffdev_nb2nb_calc_QNB = q / n
        end if
    end if

end function ffdev_nb2nb_calc_QNB

! ==============================================================================
! function ffdev_nb2nb_calc_QLJ
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_calc_QLJ(r0,eps)

    implicit none
    real(DEVDP)     :: r0,eps
    ! --------------------------------------------
    integer         :: i
    real(DEVDP)     :: r,e,q,sig,n,bi,en,dr
    ! --------------------------------------------------------------------------

    ffdev_nb2nb_calc_QLJ = 0.0d0

    sig  = r0 / 2.0d0**(1.0d0/6.0d0)

    if( NB2NBUseGaussQuad ) then
        bi = sig
        en = NB2NBCutoffR + sig
        do i=1,NB2NBGaussQuadOrder
            r = 0.5d0*(en-bi)*NB2NBGaussQuadA(i)+0.5d0*(en+bi)
            e = ffdev_nb2nb_ljene_vdw_ene(r0,eps,r)
            q = 0.0d0
            if( e .lt. 0.0d0 ) then
                q = exp(-e/((DEV_Rgas*NB2NBTemp)))
            end if
            ffdev_nb2nb_calc_QLJ = ffdev_nb2nb_calc_QLJ + q*NB2NBGaussQuadW(i)
        end do

        ffdev_nb2nb_calc_QLJ = 0.5d0*(en-bi)*ffdev_nb2nb_calc_QLJ
        ffdev_nb2nb_calc_QLJ = ffdev_nb2nb_calc_QLJ/NB2NBCutoffR
    else

        ! start from sigma derived from r0
        sig  = r0 / 2.0d0**(1.0d0/6.0d0)
        r = sig

        dr = NB2NBCutoffR/real(NB2NBNBins,DEVDP)
        q = 0.0d0
        n = 0.0d0
        do i=1,NB2NBNBins
            e = ffdev_nb2nb_ljene_vdw_ene(r0,eps,r)
            if( e .lt. 0.0d0 ) then
                q = q + exp(-e/((DEV_Rgas*NB2NBTemp)))*dr
            end if
            r = r + dr
            n = n + dr
        end do

        ffdev_nb2nb_calc_QLJ = 1.0d0
        if( n .gt. 0.0d0 ) then
            ffdev_nb2nb_calc_QLJ = q / n
        end if

    end if

end function ffdev_nb2nb_calc_QLJ

! ==============================================================================
! function ffdev_nb2nb_ljene_vdw_ene
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_ljene_vdw_ene(r0,eps,r)

    implicit none
    real(DEVDP)     :: r0,eps
    real(DEVDP)     :: r
    ! --------------------------------------------------------------------------

    ffdev_nb2nb_ljene_vdw_ene = 0.0d0
    if( r .gt. 0 ) then
        ffdev_nb2nb_ljene_vdw_ene = eps*( (r0/r)**12 - 2.0d0*(r0/r)**6 )
    end if

end function ffdev_nb2nb_ljene_vdw_ene

! ==============================================================================
! function ffdev_nb2nb_ljene
! ==============================================================================

subroutine ffdev_nb2nb_ljene(r0,eps,r,nbene)

    implicit none
    real(DEVDP)             :: r0,eps
    real(DEVDP)             :: r
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    nbene%ele_ene = 0.0
    nbene%pen_ene = 0.0
    nbene%ind_ene = 0.0
    nbene%rep_ene = 0.0
    nbene%dis_ene = 0.0
    nbene%tot_ene = 0.0

    if( r .gt. 0 ) then
        nbene%rep_ene = eps* (r0/r)**12
        nbene%dis_ene = - eps * 2.0d0 * (r0/r)**6
        nbene%tot_ene = nbene%rep_ene + nbene%dis_ene
    end if

end subroutine ffdev_nb2nb_ljene

! ==============================================================================
! function ffdev_nb2nb_nbpair
! ==============================================================================

real(DEVDP) function ffdev_nb2nb_nbpair_vdw_ene(nbpair,r)

    use ffdev_nb2nb_dat
    use ffdev_utils

    use ffdev_nbmode_LJ
    use ffdev_nbmode_EXP_DISPBJ
    use ffdev_nbmode_EXP_DISPTT

    implicit none
    type(NB_PAIR)           :: nbpair
    real(DEVDP)             :: r
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(NB_VDW_LJ)
            call ffdev_energy_nbpair_LJ(nbpair,r,nbene)

        case(NB_VDW_EXP_DISPBJ)
            call ffdev_energy_nbpair_EXP_DISPBJ(nbpair,r,nbene)

        case(NB_VDW_EXP_DISPTT)
            call ffdev_energy_nbpair_EXP_DISPTT(nbpair,r,nbene)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_nb2nb_nbpair!')
    end select

    ffdev_nb2nb_nbpair_vdw_ene = nbene%rep_ene + nbene%dis_ene

    if( NB2NBIncludePen ) then
        ffdev_nb2nb_nbpair_vdw_ene = ffdev_nb2nb_nbpair_vdw_ene + nbene%pen_ene
    end if

    if( NB2NBIncludeInd ) then
        ffdev_nb2nb_nbpair_vdw_ene = ffdev_nb2nb_nbpair_vdw_ene + nbene%ind_ene
    end if

end function ffdev_nb2nb_nbpair_vdw_ene

! ==============================================================================
! function ffdev_nb2nb_nbpair
! ==============================================================================

subroutine ffdev_nb2nb_nbpair(nbpair,r,nbene)

    use ffdev_nb2nb_dat
    use ffdev_utils

    use ffdev_nbmode_LJ
    use ffdev_nbmode_EXP_DISPBJ
    use ffdev_nbmode_EXP_DISPTT

    implicit none
    type(NB_PAIR)           :: nbpair
    real(DEVDP)             :: r
    type(NB_PAIR_ENERGY)    :: nbene
    ! --------------------------------------------------------------------------

    select case(nb_mode)
        case(NB_VDW_LJ)
            call ffdev_energy_nbpair_LJ(nbpair,r,nbene)

        case(NB_VDW_EXP_DISPBJ)
            call ffdev_energy_nbpair_EXP_DISPBJ(nbpair,r,nbene)

        case(NB_VDW_EXP_DISPTT)
            call ffdev_energy_nbpair_EXP_DISPTT(nbpair,r,nbene)

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported in ffdev_nb2nb_nbpair!')
    end select

end subroutine ffdev_nb2nb_nbpair

! ------------------------------------------------------------------------------

end module ffdev_nb2nb

!===============================================================================
! subroutine opt_shark_fce3
!===============================================================================

subroutine opt_shark_fce3(n, x, errval)

    use ffdev_nb2nb_dat
    use ffdev_nb2nb

    implicit none
    integer         :: n
    real(DEVDP)     :: x(n)
    real(DEVDP)     :: errval
    ! --------------------------------------------
    real(DEVDP)     :: eps,r0,qlj
    ! --------------------------------------------------------------------------

    call ffdev_nb2nb_nb2nb_unbox_prms(x,NB2NBprms)

    ! calculate QNB for LJ
    if( QNBModeEps ) then
        eps = NB2NBprms(1)
        r0 = QNBR0
    else
        eps = QNBEps
        r0 = NB2NBprms(1)
    end if

    qlj = ffdev_nb2nb_calc_QLJ(r0,eps)

    errval = (qlj-QNBTrg)**2

    NB2NBErrFceEval = NB2NBErrFceEval + 1

end subroutine opt_shark_fce3



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

module ffdev_err_aimxdm

use ffdev_constants
use ffdev_variables

contains

! ==============================================================================
! subroutine ffdev_err_aimxdm_init
! ==============================================================================

subroutine ffdev_err_aimxdm_init

    use ffdev_err_aimxdm_dat
    use ffdev_errors_dat

    implicit none
    ! --------------------------------------------------------------------------

    EnableAIMXDMError        = .false.
    PrintAIMXDMErrorSummary  = .false.
    AIMXDMErrorWeight        = 1.0
    AIMXDMBurriedOnly        = .false.
    TTRIJSource              = 0
    RvdWPower                = 3.0d0

end subroutine ffdev_err_aimxdm_init

! ==============================================================================
! subroutine ffdev_err_aimxdm_error
! ==============================================================================

subroutine ffdev_err_aimxdm_error(error)

    use ffdev_targetset_dat
    use ffdev_utils
    use ffdev_errors_dat
    use ffdev_err_aimxdm_dat
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat
    use ffdev_xdm_dat
    use ffdev_atomicdata_db

    implicit none
    type(FFERROR_TYPE)      :: error
    ! --------------------------------------------
    integer                 :: i,j,ip,nbt,zi,zj,ai,aj
    real(DEVDP)             :: r0,c6,c6xdm,c8xdm,c10xdm,c6eff,err_c6,err_r0
    real(DEVDP)             :: pbfree_i,pbfree_j,r0free_i,r0free_j
    real(DEVDP)             :: vaim_i,vaim_j,v0free_i,v0free_j
    real(DEVDP)             :: rii,rjj,rij
    real(DEVDP)             :: bii,bjj,bij
    real(DEVDP)             :: totc,totr,totw
    real(DEVDP)             :: w6,w8,w10,ttrij
    ! --------------------------------------------------------------------------

    error%aimxdm = 0.0d0

    totr = 0.0d0
    totc = 0.0d0
    totw = 0.0d0

    select case(nb_mode)
        case(NB_VDW_LJ)
            ! OK

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode in ffdev_err_aimxdm_summary!')
    end select

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_xdm_loaded ) cycle

            do ip=1,sets(i)%top%nb_size
                nbt = sets(i)%top%nb_list(ip)%nbt

                ai  = sets(i)%top%nb_list(ip)%ai
                aj  = sets(i)%top%nb_list(ip)%aj

                zi  = sets(i)%top%atom_types(sets(i)%top%nb_types(nbt)%ti)%z
                zj  = sets(i)%top%atom_types(sets(i)%top%nb_types(nbt)%tj)%z

                if( sets(i)%top%nb_types(nbt)%eps .eq. 0.0d0 ) cycle
                if( sets(i)%top%nb_types(nbt)%r0 .lt. 0.5d0 ) cycle

                ! r0, c6
                r0  = sets(i)%top%nb_types(nbt)%r0
                c6  = 2.0d0 * sets(i)%top%nb_types(nbt)%eps * sets(i)%top%nb_types(nbt)%r0**6

                ! XDM - pair data
                c6xdm = sets(i)%geo(j)%sup_xdm_c6(ai,aj)    * DEV_HARTREE2KCL * DEV_AU2A**6
                c8xdm = sets(i)%geo(j)%sup_xdm_c8(ai,aj)    * DEV_HARTREE2KCL * DEV_AU2A**8
                c10xdm = sets(i)%geo(j)%sup_xdm_c10(ai,aj)  * DEV_HARTREE2KCL * DEV_AU2A**10

                ! XDM - atom data
                vaim_i = sets(i)%geo(j)%sup_xdm_vol(ai)
                vaim_j = sets(i)%geo(j)%sup_xdm_vol(aj)

                v0free_i = sets(i)%geo(j)%sup_xdm_vol0(ai)
                v0free_j = sets(i)%geo(j)%sup_xdm_vol0(aj)

                ! atomic data
                r0free_i = atomicdata_vdw_r0free(zi)
                r0free_j = atomicdata_vdw_r0free(zj)

                pbfree_i = atomicdata_vdw_pbfree(zi)
                pbfree_j = atomicdata_vdw_pbfree(zj)

                ! DOI: 10.1021/ct200602x - eq 6
                bii = pbfree_i * (v0free_i / vaim_i)**(1.0d0 / 3.0d0)
                bjj = pbfree_j * (v0free_j / vaim_j)**(1.0d0 / 3.0d0)

               ! write(*,*) bii, bjj

                ! DOI: 10.1021/ct200602x - eq 5
                bij = 2.0d0 * bii * bjj / (bii + bjj)
                ! bij = sqrt(bii * bjj)


                ! DOI: 10.1021/acs.jctc.6b00027 - eq 8
                rii = r0free_i * (vaim_i / v0free_i)**(1.0d0 / RvdWPower)
                rjj = r0free_j * (vaim_j / v0free_J)**(1.0d0 / RvdWPower)

                ! Lorentz(-Berthelot) rules
                rij = 0.5d0 * (rii + rjj)

                select case(TTRIJSource)
                    case(0)
                        ttrij = r0
                    case(1)
                        ttrij = rij
                    case default
                        stop 'TTRIJSource'
                end select

                ! c6eff
                call ffdev_err_aimxdm_wx(bij*ttrij,w6,w8,w10)

                c6eff = (c6xdm*w6 + c8xdm*w8/ttrij**2 + c10xdm*w10/ttrij**4)
             !   write(*,*) w6, w8, w10, c6xdm, c6eff, c6

                write(1035,*) zi, zj, sets(i)%top%atom_types(sets(i)%top%nb_types(nbt)%ti)%name, &
                              sets(i)%top%atom_types(sets(i)%top%nb_types(nbt)%tj)%name, &
                              c6, c6eff, c6xdm

                ! errors
                err_r0 = r0 - rij
                err_c6 = c6 - c6eff

                totr = totr + err_r0**2
                totc = totc + (err_c6**2)**(1.0d0/6.0d0)

              !  write(*,*) r0,rij,err_r0,c6,c6eff,err_c6

                totw = totw + 1.0d0

            end do
        end do
    end do

    if( totw .gt. 0.0d0 ) then
        totr = sqrt( totr / totw )
        totc = sqrt( totc / totw )
    end if

   ! write(*,*) totr, totc

    error%aimxdm = totr + totc

end subroutine ffdev_err_aimxdm_error

! ==============================================================================

subroutine ffdev_err_aimxdm_wx(br,w6,w8,w10)

    real(DEVDP)             :: br
    real(DEVDP)             :: w6,w8,w10,e
    ! --------------------------------------------------------------------------

    e = exp(-br)

    w6 =  1.0d0                                                     &
        + br**1/(1.0d0)                                             &
        + br**2/(1.0d0 * 2.0d0)                                     &
        + br**3/(1.0d0 * 2.0d0 * 3.0d0)                             &
        + br**4/(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0)                     &
        + br**5/(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0 * 5.0d0)             &
        + br**6/(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0 * 5.0d0 * 6.0d0)

    w8 = w6 +                                                                       &
        + br**7/(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0 * 5.0d0 * 6.0d0 * 7.0d0)             &
        + br**8/(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0 * 5.0d0 * 6.0d0 * 7.0d0 * 8.0d0)

    w10 = w8 +                                                                                      &
        + br**9 /(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0 * 5.0d0 * 6.0d0 * 7.0d0 * 8.0d0 * 9.0d0)            &
        + br**10/(1.0d0 * 2.0d0 * 3.0d0 * 4.0d0 * 5.0d0 * 6.0d0 * 7.0d0 * 8.0d0 * 9.0d0 * 10.0d0)

    w6  = 1.0d0 - e*w6
    w8  = 1.0d0 - e*w8
    w10 = 1.0d0 - e*w10

end subroutine ffdev_err_aimxdm_wx

! ==============================================================================
! subroutine ffdev_err_aimxdm_summary
! ==============================================================================

subroutine ffdev_err_aimxdm_summary()

    use ffdev_targetset_dat
    use ffdev_err_aimxdm_dat
    use ffdev_utils
    use ffdev_parameters_dat
    use ffdev_atomicdata
    use ffdev_buried_dat
    use ffdev_xdm_dat
    use ffdev_atomicdata_db

    implicit none
    integer                 :: i,j,ip,nbt,zi,zj,ai,aj
    real(DEVDP)             :: r0,c6,c6xdm,c8xdm,c10xdm,c6eff,err_c6,err_r0
    real(DEVDP)             :: pbfree_i,pbfree_j,r0free_i,r0free_j
    real(DEVDP)             :: vaim_i,vaim_j,v0free_i,v0free_j
    real(DEVDP)             :: rii,rjj,rij
    real(DEVDP)             :: bii,bjj,bij
    real(DEVDP)             :: totc,totr,totw
    real(DEVDP)             :: w6,w8,w10,ttrij
    ! --------------------------------------------------------------------------

    write(DEV_OUT,*)
    write(DEV_OUT,5)
    write(DEV_OUT,10)
    write(DEV_OUT,20)

    totr = 0.0d0
    totc = 0.0d0
    totw = 0.0d0

    select case(nb_mode)
        case(NB_VDW_LJ)
            ! OK

        case default
            call ffdev_utils_exit(DEV_ERR,1,'Unsupported nb_mode in ffdev_err_aimxdm_summary!')
    end select

    do i=1,nsets
        do j=1,sets(i)%ngeos
            if( .not. sets(i)%geo(j)%sup_xdm_loaded ) cycle

            do ip=1,sets(i)%top%nb_size
                nbt = sets(i)%top%nb_list(ip)%nbt

                ai  = sets(i)%top%nb_list(ip)%ai
                aj  = sets(i)%top%nb_list(ip)%aj

                zi  = sets(i)%top%atom_types(sets(i)%top%nb_types(nbt)%ti)%z
                zj  = sets(i)%top%atom_types(sets(i)%top%nb_types(nbt)%tj)%z

                if( sets(i)%top%nb_types(nbt)%eps .eq. 0.0d0 ) cycle
                if( sets(i)%top%nb_types(nbt)%r0 .lt. 0.5d0 ) cycle

                ! r0, c6
                r0  = sets(i)%top%nb_types(nbt)%r0
                c6  = 2.0d0 * sets(i)%top%nb_types(nbt)%eps * sets(i)%top%nb_types(nbt)%r0**6

                ! XDM - pair data
                c6xdm = sets(i)%geo(j)%sup_xdm_c6(ai,aj)    * DEV_HARTREE2KCL * DEV_AU2A**6
                c8xdm = sets(i)%geo(j)%sup_xdm_c8(ai,aj)    * DEV_HARTREE2KCL * DEV_AU2A**8
                c10xdm = sets(i)%geo(j)%sup_xdm_c10(ai,aj)  * DEV_HARTREE2KCL * DEV_AU2A**10

                ! XDM - atom data
                vaim_i = sets(i)%geo(j)%sup_xdm_vol(ai)
                vaim_j = sets(i)%geo(j)%sup_xdm_vol(aj)

                v0free_i = sets(i)%geo(j)%sup_xdm_vol0(ai)
                v0free_j = sets(i)%geo(j)%sup_xdm_vol0(aj)

                ! atomic data
                r0free_i = atomicdata_vdw_r0free(zi)
                r0free_j = atomicdata_vdw_r0free(zj)

                pbfree_i = atomicdata_vdw_pbfree(zi)
                pbfree_j = atomicdata_vdw_pbfree(zj)

                ! DOI: 10.1021/ct200602x - eq 6
                bii = pbfree_i * (v0free_i / vaim_i)**(1.0d0 / 3.0d0)
                bjj = pbfree_j * (v0free_j / vaim_j)**(1.0d0 / 3.0d0)

               ! write(*,*) bii, bjj

                ! DOI: 10.1021/ct200602x - eq 5
                bij = 2.0d0 * bii * bjj / (bii + bjj)

                ! DOI: 10.1021/acs.jctc.6b00027 - eq 8
                rii = r0free_i * (vaim_i / v0free_i)**(1.0d0 / RvdWPower)
                rjj = r0free_j * (vaim_j / v0free_J)**(1.0d0 / RvdWPower)

                ! Lorentz(-Berthelot) rules
                rij = 0.5d0 * (rii + rjj)

                select case(TTRIJSource)
                    case(0)
                        ttrij = r0
                    case(1)
                        ttrij = rij
                    case default
                        stop 'TTRIJSource'
                end select

                ! c6eff
                call ffdev_err_aimxdm_wx(bij*ttrij,w6,w8,w10)

                c6eff = (c6xdm*w6 + c8xdm*w8/ttrij**2 + c10xdm*w10/ttrij**4)

                ! errors
                err_r0 = r0 - rij
                err_c6 = c6 - c6eff

                totr = totr + err_r0**2
                totc = totc + (err_c6**2)**(1.0d0/6.0d0)

                ! write(*,*) r0,rij,err_r0,c6,c6eff,err_c6

                totw = totw + 1.0d0

            end do
        end do
    end do

    if( totw .gt. 0.0d0 ) then
        totr = sqrt( totr / totw )
        totc = sqrt( totc / totw )
    end if

    write(DEV_OUT,20)
    write(DEV_OUT,40) totr
    write(DEV_OUT,40) totc

 5 format('# AIM XDM Penalties')
10 format('# ID TypA TypB R0Opt (DB)      V0       Vave         Raim      R0 (FF)     Diff       Error        Flag  ')
20 format('# -- ---- ---- ----------  ---------- ----------  ---------- ---------- ---------- ---------- ---------- ')
30 format(I4,1X,A4,1X,A4,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,F10.5,1X,A)
40 format('# Final penalty      =                          ',F10.5)
45 format('# Final penalty w/w  =                          ',F10.5)

end subroutine ffdev_err_aimxdm_summary

! ------------------------------------------------------------------------------

end module ffdev_err_aimxdm



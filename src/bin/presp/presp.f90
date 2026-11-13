! AMBER RESP software copyright and license notice
!                  Copyright (c) 1986, 1991, 1995,1997
!             Regents of the University of California
!                       All Rights Reserved.
!  This software provided pursuant to a license agreement containing
!  restrictions on its disclosure, duplication, and use. This software
!  contains confidential and proprietary information, and may not be
!  extracted or distributed, in whole or in part, for any purpose
!  whatsoever, without the express written permission of the authors.
!  This notice, and the associated author list, must be attached to
!  all copies, or extracts, of this software. Any additional
!  restrictions set forth in the license agreement also apply to this
!  software.

module resp_state

implicit none
integer, parameter          :: maxlgr=800, maxmol=300
integer                     :: maxq=0
integer                     :: inopt=0,ioutopt=0,iqopt=0,iunits=0,ialg=0,iwmod=0
integer                     :: nat=0,iuniq=0,nesp=0,natpl1=0,ihfree=1,irstrnt=1
character(len=8)            :: title(10),keywd(4,4)
real(kind=8), allocatable   :: apot(:,:),bpot(:),grad(:),awt(:,:),bwt(:)
real(kind=8)                :: ssvpot=0,chipot=0,vavrg=0
real(kind=8), allocatable   :: qcal(:),a(:,:),b(:),qwtval(:)
integer, allocatable        :: iwttyp(:),iqcntr(:)
real(kind=8), allocatable   :: grpchg(:)
integer, allocatable        :: lgrcnt(:,:)
integer                     :: nlgrng=0
real(kind=8), allocatable   :: q0(:),crd(:,:)
integer, allocatable        :: ivary(:),izan(:)
real(kind=8)                :: qwt=0,q0tot=0,rcond=1.0d-7
real(kind=8), allocatable   :: awork(:,:),bwork(:),scr1(:)
integer, allocatable        :: iscr1(:)
real(kind=8), allocatable   :: co(:,:)
real(kind=8)                :: cmas(3),dipol(3),dipmom,quad(6)
real(kind=8), allocatable   :: wtmol(:)
integer, allocatable        :: moleqv(:,:),ibeg(:),iend(:)
integer                     :: nmol=1
character(len=80)           :: input,output,qin,qout,punch,espot,qwts,esout
character                   :: owrite
contains

! ==============================================================================
! subroutine allocate_input_state
! ==============================================================================

subroutine allocate_input_state()

    integer                 :: u, ios, atoms, points, n_atoms, j
    character(len=256)      :: line
    ! --------------------------------------------------------------------------

    open(newunit=u, file=trim(espot), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'RESP: cannot open ESP input for sizing'
    n_atoms=0
    do
        read(u, '(2i5)', iostat=ios) atoms,points
        if (ios < 0) exit
        if (ios /= 0 .or. atoms < 1 .or. points < 0) error stop 'RESP: invalid ESP header while sizing arrays'
        n_atoms=n_atoms+atoms
        do j=1,atoms+points
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) error stop 'RESP: incomplete ESP data while sizing arrays'
        end do
    end do
    close(u)
    if (n_atoms < 1) error stop 'RESP: empty ESP input'
    maxq=n_atoms+maxlgr+1
    allocate(q0(maxq),crd(3,maxq),ivary(maxq),izan(maxq), qcal(maxq),b(maxq),qwtval(maxq),iwttyp(maxq), &
        & iqcntr(maxq),bpot(maxq),grad(maxq),bwt(maxq), bwork(maxq),scr1(maxq),iscr1(maxq),co(3,maxq))
    allocate(grpchg(maxlgr),lgrcnt(maxlgr,maxq), wtmol(maxmol),moleqv(4,maxmol),ibeg(maxmol),iend(maxmol))
    q0=0d0
    lgrcnt=0

end subroutine allocate_input_state

! ==============================================================================
! subroutine allocate_fit_matrices
! ==============================================================================

subroutine allocate_fit_matrices()

    integer                 :: n
    ! --------------------------------------------------------------------------

    n=max(iuniq+nlgrng+1,1)
    if(n>maxq) error stop 'RESP: fitted system exceeds array capacity'
    maxq=n
    allocate(apot(n,n),awt(n,n),a(n,n),awork(n,n))
    apot=0d0
    awt=0d0
    a=0d0
    awork=0d0

end subroutine allocate_fit_matrices

end module resp_state

! ==============================================================================
! program resp
! ==============================================================================

program resp

    use resp_state
    implicit none
    integer                 :: k, icycle
    ! --------------------------------------------------------------------------

    call filein
    call allocate_input_state()
    call amopen(5,input,'O','F')
    call readin
    if( nmol .gt. 1 )call mult_mol
    if( nmol .eq. 1 ) call reornt
    call allocate_fit_matrices()
    call matpot
    call data_prep
    icycle= 0
    10 continue
    call cycle(icycle)
    if( irstrnt .eq. 2 ) then
        do k= 1, iuniq
            qcal(k) = q0(k)
        end do
        qwt= 0.0d0
    else
        call charge_opt
    end if
    if (iwmod == 1) then
        call evlchi(ssvpot,bwt,awt,qcal,iuniq,maxq,chipot)
    else
        call evlchi(ssvpot,bpot,apot,qcal,iuniq,maxq,chipot)
    end if
    if(nmol .eq. 1) call elec_mom
    call pun_sum
    if( icycle .ne. 0 ) go to 10
    call pot_out
    if( ioutopt .eq. 1 ) then
        call wrt_pot
    end if

end program resp

! ==============================================================================
! subroutine readin
! ==============================================================================

subroutine readin

    use resp_state
    implicit none
    integer                 :: i, icard, ich, j
    character(len=80)       :: card
    namelist /cntrl/ ich, INOPT, ioutopt, iuniq, nmol ,IQOPT, ihfree,  qwt, irstrnt,iunits,ialg,rcond,iwmod
    do i= 1,maxlgr
        do j= 1,maxq
            lgrcnt(i,j)= 0
        end do
    end do
    nlgrng= 0
    read(5,1000) (TITLE(I),I=1,10)
    write(6, 1010) (TITLE(I),I=1,10)
    ich = 0
    iuniq = 0
    nmol = 1
    IQOPT = 0
    irstrnt = 1
    ihfree = 1
    qwt = 0.0005d0
    iunits = 0
    ialg = 0
    rcond = 1.0d-7
    iwmod = 0
    do icard=1,1000
        read(5,'(a80)',end=10) card
        if (card(3:7).eq.'cntrl') go to 9
    end do
    9 backspace (5)
    read(5,cntrl,end=10)
    go to 20
    10 continue
    write(6, 30000)
    stop
    20 continue
    if (ialg /= 0 .and. ialg /= 1) error stop 'RESP: ialg must be 0 (LU) or 1 (SVD)'
    if (iwmod /= 0 .and. iwmod /= 1) error stop 'RESP: iwmod must be 0 (legacy) or 1 (linear molecular weighting)'
    if (rcond < 0.0d0 .or. rcond >= 1.0d0) error stop 'RESP: rcond must be in [0,1)'
    write(6, 1030) INOPT, ioutopt, nmol, IQOPT, ihfree, irstrnt, iunits, qwt
    if (ialg == 1) write(6, 1031) ialg, rcond
    if (iwmod == 1) write(6, 1032) iwmod
    if( nmol .gt. 1 ) then
        return
    end if
    read(5,'(f10.0)') wtmol(1)
    write(6, 30001) wtmol(1)
    read(5,1000) (TITLE(I),I=1,10)
    write(6, 30002) (TITLE(I),I=1,10)
    read(5,'(2i5)')ich,iuniq
    write(6, 30003) ich, iuniq
    ibeg(1) = 1
    iend(1) = iuniq
    if (iwmod == 0) wtmol(1) = 1.d0
    do I=1,iuniq
        read(5,'(2i5)') IZAN(I),IVARY(I)
        write(6, 30004) i, IZAN(I), IVARY(I)
    end do
    call lagrange( ICH, ibeg(1),iend(1))
    if( IQOPT .gt. 1) then
        call amopen(3,qin,'O','F')
        write(6, 30005) iuniq
        read(3,'(8f10.6)',end=55,err=55) (Q0(I),i=1,iuniq)
        close(3)
        go to 59
        55 write(6, 30006)
        59 continue
    end if
    if( IQOPT .eq. 1 ) then
        write(6, 30007)
        do I=1,iuniq
            Q0(I)= 0.0d0
        end do
    end if
    write(6, 9088)
    write(6, 9128)
    write(6, 9168) ICH, iuniq, qwt
    write(6, 30008) nlgrng
    write(6, 30009)
    do i= 1, iuniq
        write(6, 30010) i, (lgrcnt(j,i), j= 1,nlgrng)
    end do
    if (iuniq .LE. maxq ) RETURN
    write(6, 30011)
    stop

 1000 format(10A8)
 1010 format(/,t2,'-----------------------------------------------', /,t2, &
        & '     Restrained ESP Fit 2.3  Amber 4.1', /,t2,'-----------------------------------------------', &
        & /,t2,10A8, /,t2,'-----------------------------------------------'/)
 1030 format( /t2,'inopt       = ',I5,'   ioutopt     = ',I5, /t2,'nmol        = ',I5, &
        & '   iqopt       = ',I5, /t2,'ihfree      = ',I5,'   irstrnt     = ',I5, /t2,'iunits      = ',i5, &
        & '   qwt         = ',f12.8)
 1031 format(t2,'ialg        = ',i5,'   rcond       = ',es12.5)
 1032 format(t2,'iwmod       = ',i5)
 9088 format(/,2X,38('--'),/,2X,'   ATOM   ',20X,  'COORDINATES',25X, 'CHARGE',/,23X,'X',13X, 'Y', 16X, &
        & 'Z',/,2X,38('--'))
 9128 format(2X,38('--'),/)
 9168 format(/t2,'Charge on the molecule(ich) =',I5, /t2,'Total number of atoms (iuniq)     =',I5, /t2, &
        & 'Weight factor on initial charge restraints(qwt)=' ,d16.5,/)
 30000 format('Sorry, you must use namelist input')
 30001 format(t2,'wtmol(1) = ',f12.6)
 30002 format(t2,'subtitle:',/ /,10a8)
 30003 format('ich = ',i3,'  iuniq = ',i3)
 30004 format(3i5)
 30005 format(t2,'new q0 values to be read' ,i4)
 30006 format(t2,' not enough (possibly none) q0 are given in file',/,t2,' ESP.Q0, so the remaining old ones will be used.')
 30007 format(t2,'IQOPT=1, all q0 values will be set to 0'/)
 30008 format(/t2,' there are',i3,' charge constraints:')
 30009 format(/)
 30010 format(t2,i5,5x,6(20i3,/,10x),20i3)
 30011 format(t2,'Number of atoms exceeds program dimensions')

end subroutine readin

! ==============================================================================
! subroutine mult_mol
! ==============================================================================

subroutine mult_mol

    use resp_state
    implicit none
    integer                 :: i, iat, ich, icnt, icntrs, imol, j, jcnt, mmlgr, ntmp1, ntmp2
    integer                 :: itmp(maxq),imoll(maxq)
    imol  = 0
    iuniq = 0
    mmlgr= nlgrng
    nlgrng= 0
    write(6, 30000) nmol
    do imol= 1, nmol
        read(5,'(f10.5)') wtmol(imol)
        write(6, 30001) imol, wtmol(imol)
        read(5,1000) (TITLE(I),I=1,10)
        write(6, 1010) (TITLE(I),I=1,10)
        read(5,'(2i5)') ICH,icntrs
        write(6, 1030) ICH, icntrs
        ibeg(imol)= iuniq+1
        iend(imol)= iuniq+icntrs
        if( iend(imol) .gt. maxq ) then
            write(6, 30002) maxq
            stop
        end if
        do I= ibeg(imol),iend(imol)
            read(5,'(2i5)') IZAN(I), IVARY(I)
            write(6, 30003) i, IZAN(I), IVARY(I)
            if (IVARY(i) .gt. 0) IVARY(i)= IVARY(i)+ibeg(imol)-1
        end do
        iuniq= iend(imol)
        call lagrange(ICH, ibeg(imol), iend(imol))
    end do
    call lagrange( -99, 1, iuniq)
    if( IQOPT .gt. 1 ) then
        call amopen(3,qin,'O','F')
        write(6, 30004) iuniq
        write(6, 30005)
        read(3,'(8f10.6)',end=55,err=55) (Q0(I),i=1,iuniq)
        close(3)
        goto 59
        55 write(6, 30006)
        59 continue
    end if
    if( IQOPT .eq. 1 ) then
        write(6, 30007)
        do I=1,iuniq
            Q0(I)= 0.0
        end do
    end if
    write(6, 30008)
    write(6, 30009)
    write(6, 30010)
    600 read(5,'(i5)') ntmp1
    if( ntmp1 .gt. 0 ) then
        read(5,'(16i5)') (imoll(j),itmp( j), j= 1,ntmp1)
        write(6, 30011) (imoll(j),itmp( j), j= 1,ntmp1)
        do i = 1,ntmp1
            icntrs = ibeg(imoll(i)) - 1
            itmp(i) = icntrs + itmp(i)
        end do
        do i= 2,ntmp1
            IVARY( itmp(i))= itmp(1)
        end do
        go to 600
    end if
    do i= 1,iuniq
        ntmp1= IVARY(i)
        if( ntmp1 .gt. 0 ) then
            ntmp2= IVARY(ntmp1)
            if( ntmp2 .gt. 0 ) then
                IVARY(i)= ntmp2
            end if
        end if
    end do
    write(6, 9088)
    icnt = 1
    jcnt = 1
    do IAT = 1,iuniq
        write(6, 30012) IZAN(IAT), ivary(iat)
        jcnt = jcnt + 1
        if (jcnt.gt.iend(icnt))then
            write(6, 30013)
            icnt = icnt + 1
        end if
    end do
    write(6, 9128)
    write(6, 9168) iuniq, qwt
    write(6, 30014) nlgrng
    return

 1000 format(10A8)
 1010 format(10A8)
 1030 format(/,t2,'Total charge (ich):',I3, /,t2,'Number of centers:',I3)
 9088 format(/,2X,10('--'),/,5x,'Atom   Ivary',/,2X,10('--'))
 9128 format(2X,38('--'),/)
 9168 format(/t2,'Total number of atoms      =',I5, /t2,'Weight factor on initial charge restraints=',f10.6,/)
 30000 format(/,t2,'%RESP-I-MULT_MOL,  multiple-molecule run of ',i3,' molecules')
 30001 format(/,t2,'Reading input for molecule ',i3,' weight:',f10.3)
 30002 format(t2,'ERROR: more than ',i5,' centers')
 30003 format(3i5)
 30004 format(t2,' since IQOPT=1,',i4,' new q0 values')
 30005 format(t2,' will be read in from file ESP.Q0 (unit 3)')
 30006 format(t2,'not enough (possibly none) q0 are given in file',/,t2,'ESP.Q0 (unit 3), so the remaining old ones will be used.')
 30007 format(/t2,'Iqopt =1: all q0 values will be set to 0')
 30008 format(t2,'--------------------------------')
 30009 format(t2,'reading mult_mol constraint info')
 30010 format(t2,'--------------------------------')
 30011 format(16i5)
 30012 format(2i5)
 30013 format( )
 30014 format(/,t2,'There are',i3,' charge constraints')

end subroutine mult_mol

! ==============================================================================
! subroutine lagrange
! ==============================================================================

subroutine lagrange(ncharge,ifirst,ilast)

    use resp_state
    implicit none
    integer                 :: icntrs, ifirst, ii, ilast, j, jmol, ncharge, ntmp
    real(kind=8)            :: gtemp
    integer                 :: itmp(maxq),imoll(maxq)
    ! --------------------------------------------------------------------------

    10 continue
    read(5,'(i5,f10.5)') ntmp, gtemp
    if(ntmp.eq.0) go to 20
    nlgrng= nlgrng + 1
    if( (nlgrng) .gt. maxlgr) then
        write(6, 30000) nlgrng, maxlgr
        stop
    end if
    grpchg(nlgrng) = gtemp
    read(5,'(16i5)') (imoll(j),itmp( j), j= 1,ntmp)
    write(6, 30001) (imoll(j),itmp( j), j= 1,ntmp)
    do ii = 1,ntmp
        jmol = imoll(ii)
        icntrs = ibeg(jmol) - 1
        itmp(ii) = icntrs + itmp(ii)
    end do
    do j= 1, ntmp
        if( itmp(j) .gt. 0 ) then
            itmp(j)= itmp(j) + ifirst -1
            lgrcnt(nlgrng, itmp(j))= 1
        else if( itmp(j) .lt. 0 ) then
            itmp(j)= itmp(j) - ifirst +1
            lgrcnt(nlgrng, itmp(j))= -1
        end if
    end do
    go to 10
    20 continue
    if( ncharge .gt. -99 ) then
        nlgrng= nlgrng + 1
        if( nlgrng .gt. maxlgr) then
            write(6, 30002) nlgrng, maxlgr
            stop
        end if
        grpchg(nlgrng)= float( ncharge)
        do  j= ifirst, ilast
            lgrcnt(nlgrng,j)= 1
        end do
    end if
    return

 30000 format(t2,' Too many charge-group constraints',i3,/,t2,' Maximum allowed: ',i3)
 30001 format(16i5)
 30002 format(t2,' Too many charge-group constraints',i3,/,t2,' Maximum allowed: ',i3)

end subroutine lagrange

! ==============================================================================
! subroutine matpot
! ==============================================================================

subroutine matpot

    use resp_state
    implicit none
    integer                 :: i, imol, inat, inmol, ioff, j, k
    real(kind=8)            :: espi, rij, rik, vavtmp, wt, wt2, xi, yi, zi
    ! --------------------------------------------------------------------------

    do k=1,iuniq
        bpot(k)= 0.0d0
        bwt(k) = 0.0d0
        do j=1,iuniq
            apot(j,k)= 0.0d0
            awt(j,k) = 0.0d0
        end do
    end do
    call amopen(10,espot,'O','F')
    ssvpot= 0.0d0
    vavrg = 0.0d0
    vavtmp = 0.0d0
    ioff = 1
    if( nmol .gt. 0 ) then
        inmol= nmol
    else
        inmol= 1
        ibeg(1)= 1
        iend(1)= iuniq
        wtmol(1)= 1.0d0
    end if
    do imol= 1, inmol
        read(10,'(2i5)') inat,nesp
        write(6, 30000) imol, inat, NESP
        write(6, 30001)
        do i = 1,inat
            read(10,52)crd(1,ioff),crd(2,ioff),crd(3,ioff)
            write(6, 152) i, crd(1,ioff), crd(2,ioff), crd(3,ioff)
            ioff = ioff + 1
        end do
        do i= 1,nesp
            read(10, 53, err=940, end=930) espi, xi, yi, zi
            wt= wtmol(imol)
            if (iwmod == 1) then
                wt2= wt
            else
                wt2= wt*wt
            end if
            vavtmp = vavtmp + wt*espi
            ssvpot = ssvpot + wt2*espi*espi
            vavrg  = vavrg  + vavtmp/float(NESP)
            do k= ibeg(imol),iend(imol)
                rik = sqrt( (xi - CRD(1,k))**2 + (yi - CRD(2,k))**2 + (zi - CRD(3,k))**2)
                rik = 1.d0/rik
                bpot(k)  = bpot(k) +     espi*rik
                bwt(k)   = bwt(k)  + wt2*espi*rik
                apot(k,k)= apot(k,k) +     rik*rik
                awt(k,k) = awt(k,k)  + wt2*rik*rik
                do j= k+1, iend(imol)
                    rij = sqrt((xi - CRD(1,j))**2+(yi - CRD(2,j))**2+ (zi - CRD(3,j))**2)
                    rij = 1.d0/rij
                    apot(j,k)= apot(j,k) +     rij*rik
                    awt(j,k) = awt(j,k)  + wt2*rij*rik
                end do
            end do
        end do
    end do
    do k=1,iuniq
        do j= k+1,iuniq
            awt(k,j)= awt(j,k)
            apot(k,j)= apot(j,k)
        end do
    end do
    close(10)
    write(6, 900) ssvpot
    return
    930 write(6, 30002)
    stop
    940 write(6, 30003)
    stop

 52 format(17X,3e16.7)
 152 format(1X,i4,3e16.7)
 53 format(1X,4e16.7)
 900 format(t2,'Initial ssvpot =',f10.3/)
 30000 format(/,t2,'Reading esp"s for molecule ',i3,/,t2,'total number of atoms      = ',i5,/,t2,'total number of esp points = ',i5)
 30001 format(/ /,' center     X       Y       Z ')
 30002 format(5X,'premature end of potential file')
 30003 format(5X,'Error in reading potential input file (1X,4e16.7)')

end subroutine matpot

! ==============================================================================
! subroutine data_prep
! ==============================================================================

subroutine data_prep

    use resp_state
    implicit none
    integer                 :: i, j
    real(kind=8)            :: qcntrs, tmpctr
    ! --------------------------------------------------------------------------

    NAT = 0
    do i= 1,iuniq
        if(IVARY(i) .eq. 0) then
            NAT= NAT + 1
            iqcntr( i)= NAT
        else if(IVARY(i) .gt. 0) then
            iqcntr(i)= iqcntr(IVARY(i))
            if( iqcntr( IVARY(i)) .gt. NAT ) then
                write(6, 30000)
                stop
            end if
        else
            iqcntr( i)= -1
        end if
    end do
    write(6, 30001) NAT
    if( NAT .eq. 0 ) then
        write(6, 30002)
    end if
    do i= 1,nlgrng
        iqcntr(iuniq+i)= nat + i
    end do
    natpl1= nat + nlgrng
    if( IQOPT .eq. 3 ) then
        do i= 1,iuniq-1
            qcntrs= q0(i)
            tmpctr= 1.0d0
            do j= i+1,iuniq
                if( IVARY(j) .eq. i ) then
                    qcntrs= qcntrs + q0(j)
                    tmpctr= tmpctr + 1.0d0
                end if
            end do
            if(tmpctr .gt. 0.99d0 )then
                qcntrs= qcntrs/tmpctr
                q0(i)= qcntrs
                do j= i+1,iuniq
                    if( IVARY(j) .eq. i ) q0(j)= qcntrs
                end do
            end if
        end do
    end if
    return

 30000 format(t2,'data_prep: charge vary input is screwy')
 30001 format(/t2,'Number of unique UNfrozen centers=',i5)
 30002 format(t2,'ALL charges are frozen!!!')

end subroutine data_prep

! ==============================================================================
! subroutine matbld
! ==============================================================================

subroutine matbld

    use resp_state
    implicit none
    integer                 :: i, icntr, j, jcntr, k
    ! --------------------------------------------------------------------------

    do k=1,iuniq
        b(k) = bwt(k)
        do J=1,iuniq
            a(j,k) = awt(j,k)
        end do
    end do
    do i=1,nlgrng
        b(iuniq+i)  = grpchg(i)
        do j= 1,iuniq+nlgrng
            a(iuniq+i,j)= float(lgrcnt(i,j))
            a(j,iuniq+i)= float(lgrcnt(i,j))
        end do
    end do
    continue
    call rstran
    do i= 1,natpl1
        bwork( i)   = 0.0d0
        awork( i, i)= 0.0d0
        do j= i+1,natpl1
            awork( j, i)= 0.0d0
            awork( i, j)= 0.0d0
        end do
    end do
    do i= 1,iuniq+nlgrng
        icntr= iqcntr(i)
        if(icntr .ge. 1) then
            bwork(icntr)= bwork(icntr) + b(i)
            do j = 1,iuniq+nlgrng
                jcntr = iqcntr( j)
                if( jcntr .ge. 1 ) then
                    awork(icntr,jcntr)= awork(icntr,jcntr) + a(i,j)
                else
                    bwork(icntr)= bwork(icntr)       - q0(j)*a(i,j)
                end if
            end do
        end if
    end do
    return

end subroutine matbld

! ==============================================================================
! subroutine rstran
! ==============================================================================

subroutine rstran

    use resp_state
    implicit none
    integer                 :: i, ihit, imol
    real(kind=8)            :: restraint_weight
    ! --------------------------------------------------------------------------

    imol=1
    do i = 1,iuniq
        restraint_weight=1.0d0
        if (iwmod == 1) then
            if (i > iend(imol)) imol=imol+1
            restraint_weight=wtmol(imol)
        end if
        qwtval(i)= qwt*restraint_weight
        if(ihfree.gt.0.and. IZAN(i).eq.1) then
            qwtval(i)= 0.0D0
        end if
        if(irstrnt .eq. 0) then
            a(i,i)= a(i,i) + qwtval(i)
            b(i)  = b(i)   + qwtval(i)*q0(i)
        else if (irstrnt .gt. 0 .and. qwtval(i) .gt. 0.1d-10) then
            qwtval(i) = qwtval(i)/sqrt(qcal(i)*qcal(i) + 0.01d0)
            a(i,i)= a(i,i) + qwtval(i)
        end if
    end do
    ihit = 0
    do i= 1,iuniq
        if(qwtval(i) .gt. 0.1d-10 ) ihit = 1
    end do
    if(ihit .eq. 0) irstrnt = -1
    return

end subroutine rstran

! ==============================================================================
! subroutine cycle
! ==============================================================================

subroutine cycle(icycle)

    use resp_state
    implicit none
    integer                 :: icycle, nqwt
    save nqwt
    if( INOPT .lt. 1 ) return
    if( icycle .eq. 0 ) then
        nqwt= 0
        call amopen(4,qwts,'O','F')
        read(4,'(i5)') nqwt
        if(nqwt.eq.0) then
            write(6, 30000)
            write(6, 30001)
            return
        end if
        read(4,'(f10.5)') qwt
        icycle= 1
        write(6, 30002) qwt
    else
        read(4,'(f10.5)') qwt
        icycle= icycle + 1
        write(6, 30003) icycle, qwt
        if( icycle .ge. nqwt ) icycle= 0
    end if
    return

 30000 format(' subroutine icycle: INPUT ERROR...')
 30001 format(' INOPT=1 so qwt cycling expected, but',' reading non-zero nqwt failed')
 30002 format(t2,'cycle   1: weighting factor= ',f10.4)
 30003 format(/t2,'cycle',i4,': weighting factor= ',f10.4)

end subroutine cycle

! ==============================================================================
! subroutine reornt
! ==============================================================================

subroutine reornt

    use resp_state
    implicit none
    integer                 :: i, iatom
    real(kind=8)            :: bohr, wt, xc, yc, zc, zero
    dimension IATOM(maxq),WT(35)
    data ZERO/0.0D0/,BOHR/0.52917725D0/
    data (WT(I),I=1,21) / 0.0000D0, 1.0079D0, 4.0026D0, 6.9410D0, 9.0122D0,10.8100D0,12.0110D0, 14.0067D0, &
        & 15.9994D0,18.9984D0,20.1700D0, 22.9897D0,24.3050D0,26.9815D0,28.0860D0, 30.9738D0,32.0600D0, &
        & 35.4530D0,39.9480D0, 39.0900D0,40.0800D0 /
    XC = ZERO
    YC = ZERO
    ZC = ZERO
    do I = 1,iuniq
        IATOM(I) = IZAN(I)+1
    end do
    call CMASS(XC,YC,ZC,CRD,IATOM,iuniq,WT,maxq)
    CMAS(1)= XC
    CMAS(2)= YC
    CMAS(3)= ZC
    do I = 1,3
        CMAS(I) = CMAS(I)*BOHR
    end do
    call CMOVE(XC,YC,ZC,CRD,iuniq,CO,maxq)
    call MOMIN(CO,IATOM,iuniq,WT,maxq)
    return

end subroutine reornt

! ==============================================================================
! subroutine elec_mom
! ==============================================================================

subroutine elec_mom

    use resp_state
    implicit none
    integer                 :: i
    real(kind=8)            :: bohr, debye, zero
    ! --------------------------------------------------------------------------

    ZERO=0.0D0
    BOHR=0.52917725D0
    debye=2.541765d0
    DIPOL(1) =ZERO
    DIPOL(2) =ZERO
    DIPOL(3) =ZERO
    do I =1,iuniq
        DIPOL(1) =DIPOL(1) +QCAL(I)*CRD(1,I)
        DIPOL(2) =DIPOL(2) +QCAL(I)*CRD(2,I)
        DIPOL(3) =DIPOL(3) +QCAL(I)*CRD(3,I)
    end do
    continue
    dipmom= dsqrt( DIPOL(1)*DIPOL(1) + DIPOL(2)*DIPOL(2) + DIPOL(3)*DIPOL(3))
    call QUADM(CRD,QCAL,iuniq,QUAD,maxq)
    do I = 1,3
        DIPOL(I) = DIPOL(I)*debye
        QUAD(I) = QUAD(I)*debye*BOHR
        QUAD(I+3) = QUAD(I+3)*debye*BOHR
    end do
    dipmom= dipmom*debye
    return

end subroutine elec_mom

! ==============================================================================
! subroutine cmass
! ==============================================================================

subroutine CMASS(XC,YC,ZC,C,IATOM,NATOM,WT,maxq)

    implicit none
    integer                 :: i, iatom, index, maxq, natom
    real(kind=8)            :: c, sum, sumx, sumy, sumz, wt, xc, yc, zc, zero
    dimension C(3,maxq),IATOM(maxq),WT(35)
    data ZERO/0.0D0/
    SUMX = ZERO
    SUMY = ZERO
    SUMZ = ZERO
    SUM  = ZERO
    do I =1,NATOM
        INDEX=IATOM(I)
        SUMX = SUMX + C(1,I) * WT(INDEX)
        SUMY = SUMY + C(2,I) * WT(INDEX)
        SUMZ = SUMZ + C(3,I) * WT(INDEX)
    SUM = SUM + WT(INDEX)
    end do
    XC = SUMX/SUM
    YC = SUMY/SUM
    ZC = SUMZ/SUM
    return

end subroutine cmass

! ==============================================================================
! subroutine cmove
! ==============================================================================

subroutine CMOVE(XC,YC,ZC,C,NATOM,CO,maxq)

    implicit none
    integer                 :: i, maxq, natom
    real(kind=8)            :: c, co, xc, yc, zc
    dimension C(3,maxq),CO(3,maxq)
    do I=1,NATOM
        CO(1,I) = C(1,I) - XC
        CO(2,I) = C(2,I) - YC
        CO(3,I) = C(3,I) - ZC
    end do
    return

end subroutine cmove

! ==============================================================================
! subroutine momin
! ==============================================================================

subroutine MOMIN(CO,IATOM,NATOM,WT,maxq)

    implicit none
    integer                 :: i, index, info, maxq, natom
    integer                 :: iatom(maxq)
    integer, parameter      :: lwork=128
    real(kind=8)            :: ain(3,3), co(3,maxq), eigenvalues(3), sxx, sxy, sxz, syy, syz, szz, wt(35)
    real(kind=8)            :: work(lwork), xt, xx, yt, yy, zt, zz
    ! --------------------------------------------------------------------------

    SXX=0.0d0
    SYY=0.0d0
    SZZ=0.0d0
    SXY=0.0d0
    SXZ=0.0d0
    SYZ=0.0d0
    do I=1,NATOM
        INDEX=IATOM(I)
        XX=CO(1,I)*CO(1,I)
        YY=CO(2,I)*CO(2,I)
        ZZ=CO(3,I)*CO(3,I)
        SXX=SXX + (WT(INDEX)*(YY + ZZ))
        SYY=SYY + (WT(INDEX)*(XX + ZZ))
        SZZ=SZZ + (WT(INDEX)*(XX + YY))
        SXY=SXY + (WT(INDEX)*CO(1,I)*CO(2,I))
        SXZ=SXZ + (WT(INDEX)*CO(1,I)*CO(3,I))
        SYZ=SYZ + (WT(INDEX)*CO(2,I)*CO(3,I))
    end do
    AIN(1,1)=SXX
    AIN(1,2)=-SXY
    AIN(1,3)=-SXZ
    AIN(2,1)=AIN(1,2)
    AIN(2,2)=SYY
    AIN(2,3)=-SYZ
    AIN(3,1)=AIN(1,3)
    AIN(3,2)=AIN(2,3)
    AIN(3,3)=SZZ

    call dsyev('V','U',3,ain,3,eigenvalues,work,lwork,info)
    if (info /= 0) error stop 'RESP: DSYEV failed in MOMIN'

    do I=1,NATOM
        XT=AIN(1,1)*CO(1,I) + AIN(2,1)*CO(2,I) + AIN(3,1)*CO(3,I)
        YT=AIN(1,2)*CO(1,I) + AIN(2,2)*CO(2,I) + AIN(3,2)*CO(3,I)
        ZT=AIN(1,3)*CO(1,I) + AIN(2,3)*CO(2,I) + AIN(3,3)*CO(3,I)
        CO(1,I)=XT
        CO(2,I)=YT
        CO(3,I)=ZT
    end do

end subroutine momin

! ==============================================================================
! ==============================================================================
! subroutine quadm
! ==============================================================================

subroutine QUADM(C,q,NATOM,QUAD,maxq)

    implicit none
    integer                 :: i, maxq, natom
    real(kind=8)            :: c, q, quad, qxx, qxy, qxz, qyy, qyz, qzz, rq, x2, y2, z2
    dimension C(3,maxq), q(maxq), QUAD(6)
    QXX=0.0d0
    QYY=0.0d0
    QZZ=0.0d0
    QXY=0.0d0
    QXZ=0.0d0
    QYZ=0.0d0
    do I=1,NATOM
        X2 = C(1,I)*C(1,I)
        Y2 = C(2,I)*C(2,I)
        Z2 = C(3,I)*C(3,I)
        RQ = X2 + Y2 + Z2
        QXX = QXX + Q(I)*( 3.d0* (C(1,I)*C(1,I)) - RQ )
        QYY = QYY + Q(I)*( 3.d0* (C(2,I)*C(2,I)) - RQ )
        QZZ = QZZ + Q(I)*( 3.d0* (C(3,I)*C(3,I)) - RQ )
        QXY = QXY + Q(I)*( 3.d0* (C(1,I)*C(2,I)) - RQ )
        QXZ = QXZ + Q(I)*( 3.d0* (C(1,I)*C(3,I)) - RQ )
        QYZ = QYZ + Q(I)*( 3.d0* (C(2,I)*C(3,I)) - RQ )
    end do
    QUAD(1) = QXX
    QUAD(2) = QYY
    QUAD(3) = QZZ
    QUAD(4) = QXY
    QUAD(5) = QXZ
    QUAD(6) = QYZ
    return

end subroutine quadm

! ==============================================================================
! subroutine filein
! ==============================================================================

subroutine filein

    use resp_state
    implicit none
    integer                 :: indx
    character(len=80)       :: arg
    integer                 :: iarg, narg
    ! --------------------------------------------------------------------------

    input = 'input'       ! -i
    output = 'output'     ! -o
    punch = 'punch'       ! -p
    qin= 'qin'            ! -q
    qout= 'qout'          ! -t
    espot = 'espot'       ! -e
    qwts = 'qwts'         ! -w
    esout = 'esout'       ! -s
    owrite = 'N'
    indx = command_argument_count()
    iarg = 0
    if (indx.eq.iarg) goto 20
    10 continue
    iarg = iarg + 1
    call get_command_argument(iarg,arg)
    if (arg .eq. '-O') then
        owrite = 'U'
    else if (arg .eq. '-i') then
        iarg = iarg + 1
        call get_command_argument(iarg,input)
    else if (arg .eq. '-o') then
        iarg = iarg + 1
        call get_command_argument(iarg,output)
    else if (arg .eq. '-p') then
        iarg = iarg + 1
        call get_command_argument(iarg,punch)
    else if (arg .eq. '-q') then
        iarg = iarg + 1
        call get_command_argument(iarg,qin)
    else if (arg .eq. '-t') then
        iarg = iarg + 1
        call get_command_argument(iarg,qout)
    else if (arg .eq. '-e') then
        iarg = iarg + 1
        call get_command_argument(iarg,espot)
    else if (arg .eq. '-s') then
        iarg = iarg + 1
        call get_command_argument(iarg,esout)
    else if (arg .eq. '-w') then
        iarg = iarg + 1
        call get_command_argument(iarg,qwts)
    else
        if (arg .eq. ' ') go to 20
        write(6, 30000) arg
        write(6, 9000)
        stop
    end if
    if (iarg .lt. indx) go to 10
    20 continue
    narg = iarg - 1
    if (narg .lt. 2) then
        write(6, 9000)
        stop
    end if
    if(output .ne. 'screen') then
        call amopen(6,output,owrite,'F')
    end if
    return

 9000 format(/,t2, 'usage: resp [-O] -i input -o output -p punch -q qin -t qout',' -e espot -w qwts -s esout ')
 30000 format(/,5x,'unknown flag: ',a)

end subroutine filein

! ==============================================================================
! subroutine evlchi
! ==============================================================================

subroutine evlchi(ssvpot,bpot,apot,qcal,iuniq,maxq,chipot)

    implicit none
    integer                 :: iuniq, j, k, maxq
    real(kind=8)            :: apot, bpot, chipot, cross, qcal, ssvpot, ssyclc
    dimension bpot(maxq), apot(maxq,maxq), qcal(maxq)
    cross= 0.0d0
    ssyclc= 0.0d0
    do j= 1,iuniq
        cross= cross + qcal(j)*bpot(j)
        do k= 1,iuniq
            ssyclc= ssyclc + qcal(j)*qcal(k)*apot(j,k)
        end do
    end do
    chipot = ssvpot - 2.0d0*cross + ssyclc
    return

end subroutine evlchi

! ==============================================================================
! subroutine pun_sum
! ==============================================================================

subroutine pun_sum

    use resp_state
    implicit none
    integer                 :: i, icnt, j, jcnt
    real(kind=8)            :: qcrtrn, sigma, ssvtot
    ! --------------------------------------------------------------------------

    QCRTRN= SQRT(chipot/ssvpot)
    ssvtot = ssvpot
    SIGMA = SQRT(chipot/FLOAT(NESP))
    call amopen(7,punch,owrite,'F')
    write(7, 1110) (TITLE(I),I=1,10)
    write(7, 1230) IQOPT, irstrnt, ihfree, qwt
    write(7, 1235) QCRTRN, dipmom, (QUAD(I), I=1,3)
    write(7, 1200)
    icnt = 1
    jcnt = 1
    do j=1,iuniq
        write(7, 1210) j, izan(j), q0(j), qcal(j), IVARY(j), qwtval(j)
        jcnt = jcnt + 1
        if (jcnt.gt.iend(icnt))then
            write(6, 30000)
            icnt = icnt + 1
        end if
    end do
    call amopen(19,qout,owrite,'F')
    write(19, 30001) (Qcal(I),i=1,iuniq)
    close(unit=19)
    return

 1110 format(/,10A8)
 1230 format(/,'iqopt   irstrnt  ihfree     qwt'/,3(i3,4x),f12.6)
 1235 format(/t2,'rel.rms   dipole mom       Qxx      Qyy      Qzz',/ t2,5f10.5)
 1200 format(/,10X,'Point charges before & after optimization',/,4X, 'NO',3X,'At.No.',4x,'q0',11X, &
        & 'q(opt)   IVARY  d(rstr)/dq ')
 1210 format(2x,i4,2x,i4,1x,f10.6,5x,f10.6, i5, f12.6, f12.3)
 30000 format( )
 30001 format(8f10.6)

end subroutine pun_sum

! ==============================================================================
! subroutine pot_out
! ==============================================================================

subroutine POT_OUT

    use resp_state
    implicit none
    integer                 :: i, icnt, j, jcnt
    real(kind=8)            :: au2cal, chge, qcrtrn, sigma, zero
    data ZERO/0.0D0/,AU2CAL/627.5095D0/
    write(6, 1110) (TITLE(I),I=1,10)
    write(6, 1200)
    icnt = 1
    jcnt = 1
    chge = 0.0
    do j=1,iuniq
        write(6, 1210) j, izan(j), q0(j), qcal(j), IVARY(j), qwtval(j)
        chge = chge + qcal(j)
        jcnt = jcnt + 1
        if (jcnt.gt.iend(icnt))then
            write(6, 30000)
            icnt = icnt + 1
        end if
    end do
    write(6, 30001) chge
    QCRTRN= SQRT(chipot/ssvpot)
    SIGMA = SQRT(chipot/FLOAT(NESP))
    write(6, 1040) ssvpot
    write(7, 1040) ssvpot
    write(6, 1050) chipot
    write(7, 1050) chipot
    write(6, 1080) SIGMA
    write(7, 1080) SIGMA
    write(6, 1090) QCRTRN
    write(7, 1090) QCRTRN
    if(nmol.eq.1) then
        write(6, 1340)
        write(6, 1350) ( CMAS(I), I=1,3)
        write(6, 1360)
        write(6, 1370) ( DIPOL(I), I=1,3)
        write(6, 1375) dipmom
        write(7, 1375) dipmom
        write(6, 1380)
        write(6, 1390) ( QUAD(I), I=1,3)
        write(6, 1395) ( QUAD(I), I=4,6)
    end if
    return

 1040 format(/,8X,'Statistics of the fitting:', /,2x,'The initial sum of squares (ssvpot)',t50,f15.3)
 1050 format(2x,'The residual sum of squares (chipot)',t50,f15.3)
 1080 format(2x, 'The std err of estimate (sqrt(chipot/N))',t50,f15.5)
 1090 format(2x, 'ESP relative RMS (SQRT(chipot/ssvpot))',t50,f15.5)
 1110 format(/,10A8)
 1200 format(/,10X, 'Point Charges Before & After Optimization',/ /,4X,'no.', 2x,'At.no.',4x,'q(init)', &
        & 7X,'q(opt)     ivary    d(rstr)/dq')
 1210 format(t2,2i4, 2(5x,f10.6), i7, f15.6, f12.3)
 1340 format(/,t2,'Center of Mass (Angst.):',/)
 1350 format(t2,'X  =',F10.5,5X,' Y  =',F10.5,5X,' Z  =',F10.5)
 1360 format(/,t2,'Dipole (Debye):',/)
 1370 format(t2,'X  =',F10.5,5X,' Y  =',F10.5,5X,' Z  =',F10.5)
 1375 format(/,t2,'Dipole Moment (Debye)=',F10.5)
 1380 format(/,t2,'Quadrupole (Debye*Angst.):',/)
 1390 format(t2,'Qxx =',F10.5,5X,'QYY =',F10.5,5X,'QZZ =',F10.5)
 1395 format(t2,'Qxy =',F10.5,5X,'QXZ =',F10.5,5X,'QYZ =',F10.5)
 30000 format( )
 30001 format(t2,'Sum over the calculated charges: ',f10.3)

end subroutine pot_out

! ==============================================================================
! subroutine wrt_pot
! ==============================================================================

subroutine wrt_pot

    use resp_state
    implicit none
    integer                 :: i, idum, k
    real(kind=8)            :: au2cal, bohr, chikcl, espclc, espqmi, ssvkcl, vresid, xa, xi, xik, ya, yi, yik, za, zi, zik
    data AU2CAL/627.5095d0/,BOHR/0.52917725d0/
    call amopen(10, ESPOT,'O','F')
    rewind(10)
    call amopen(20, ESOUT,owrite,'F')
    read(10,'(2i5)') idum,nesp
    ssvkcl= ssvpot*au2cal*au2cal
    chikcl= chipot*au2cal*au2cal
    write(20, 30000) nesp, izan(1), ssvkcl, chikcl
    do i = 1,idum
        read(10,100, end=930)  xi, yi, zi
    end do
    do i= 1,nesp
        read(10,101, end=930) espqmi, xi, yi, zi
        espclc= 0.0d0
        do k=1,iuniq
            Xik    = xi - CRD(1,k)
            Yik    = yi - CRD(2,k)
            Zik    = zi - CRD(3,k)
            espclc= espclc + qcal(k)/ SQRT( Xik*Xik + Yik*Yik + Zik*Zik)
        end do
        vresid= espqmi - espclc
        xa= xi*BOHR
        ya= yi*BOHR
        za= zi*BOHR
        espclc= espclc*AU2CAL
        espqmi= espqmi*AU2CAL
        vresid= vresid*AU2CAL
        write(20, 30001) xa, ya, za, espqmi, espclc, vresid
    end do
    close(10)
    close(20)
    return
    930 write(6, 30008) espot
    stop

 100 format(17x,3e16.7)
 101 format(1x,4e16.7)
 30000 format(i5,i6,4x,2f20.10)
 30001 format(3f10.5, 3f12.5)
 30008 format(1x,' unexpected eof in ',1x,a)

end subroutine wrt_pot

! ==============================================================================
! subroutine charge_opt
! ==============================================================================

subroutine charge_opt

    use resp_state
    implicit none
    integer                 :: i, icntr, irsave, jn, k, maxit, nitern
    real(kind=8)            :: gs, qchnge, qdiff, qold, qtol, xuniq
    dimension qold(maxq)
    integer                 :: istat, lwork, rank
    real(kind=8)            :: work_query(1)
    real(kind=8), allocatable :: singular_values(:), work(:)
    ! --------------------------------------------------------------------------

    gs= 0.0d0
    irsave= 0
    nitern= 0
    do i= 1,maxq
        qold(i)= 0.0d0
    end do
    qtol= 0.1d-5
    maxit=  24
    if( irstrnt .gt. 0 ) then
        irsave= irstrnt
        irstrnt= 0
        write(6, 30000)
    end if
    100 continue
    call matbld
    do jn = 1,natpl1
        if (abs(awork(jn,jn)).lt. 1.d-10) awork(jn,jn) = 1.d-10
    end do
    select case (ialg)
    case (0)
        call dgesv(natpl1,1,awork,maxq,iscr1,bwork,natpl1,istat)
    case (1)
        allocate(singular_values(natpl1))
        lwork=-1
        call dgelss(natpl1,natpl1,1,awork,maxq,bwork,natpl1,singular_values,rcond,rank,work_query,lwork,istat)
        if (istat /= 0) error stop 'RESP: DGELSS workspace query failed'
        lwork=max(1,ceiling(work_query(1)))
        allocate(work(lwork))
        call dgelss(natpl1,natpl1,1,awork,maxq,bwork,natpl1,singular_values,rcond,rank,work,lwork,istat)
        deallocate(work,singular_values)
    end select
    if (istat < 0) then
        write(6, 30001) ialg
        stop
    else if (istat > 0) then
        write(6, 30002) ialg
        stop
    end if
    do k= 1, iuniq
        icntr= iqcntr(k)
        if( icntr .ge. 1 ) then
            qcal(k) = bwork( icntr)
        else
            qcal(k) = q0(k)
        end if
    end do
    if (irstrnt .lt. 0) then
        irsave=0
        write(6, 30004)
    end if
    if( irsave .le. 0 ) then
        return
    else
        irstrnt= irsave
    end if
    qchnge= 0.0d0
    xuniq = dble(iuniq)
    do i= 1,iuniq
        qdiff = qcal(i) - qold(i)
        qchnge = qchnge + (qdiff*qdiff)
        qold(i) = qcal(i)
    end do
    qchnge = sqrt(qchnge) / xuniq
    write(6, 30005) qchnge
    if (qchnge .lt. qtol .and. nitern .gt. 1) then
        write(6, 30006) nitern
        return
    else if( nitern .ge. maxit ) then
        write(6, 30007) maxit
        return
    end if
    nitern= nitern + 1
    go to 100

 30000 format(/,t2,'Non-linear optimization requested.')
 30001 format(' chgopt: solver ialg=',i1,' has an illegal argument')
 30002 format(' chgopt: solver ialg=',i1,' failed (LU singular or SVD nonconvergence)')
 30004 format(/,t2,'WARNING: Restraints were requested, but',/,t2,'         the restraint weights were all zero')
 30005 format(t2,'qchnge =',g20.10)
 30006 format(/t2,'Convergence in',i5,' iterations')
 30007 format(t2,'after ',i5,' iterations, no convergence!')

end subroutine charge_opt

! ==============================================================================
! subroutine amopen
! ==============================================================================

subroutine amopen(lun,fname,fstat,fform)

    implicit none
    integer                 :: lun
    character(len=*)        :: fname
    character(len=1)        :: fstat
    character(len=1)        :: fform
    character(len=7)        :: stat
    character(len=11)       :: kform
    integer                 :: ios
    ! --------------------------------------------------------------------------

    if (fstat .eq. 'N') then
        stat = 'NEW'
    else if (fstat .eq. 'O') then
        stat = 'OLD'
    else if (fstat .eq. 'U') then
        stat = 'UNKNOWN'
    else
        write(6, 30000) lun
        stop
    end if
    if (fform .eq. 'U') then
        kform = 'UNFORMATTED'
    else if (fform .eq. 'F') then
        kform = 'FORMATTED'
    else
        write(6, 30001) lun
        stop
    end if
    open(unit=lun,file=fname,status=stat,form=kform,iostat=ios)
    if (ios .ne. 0) then
        if (lun .eq. 6) then
            write(0, 30002) lun, fname
            close(unit=0)
        else
            write(6, 30003) lun, fname
            close(unit=6)
            write(0, 30004) lun, fname
            close(unit=0)
        end if
        stop
    end if
    rewind(lun)
    return

 30000 format(/,2x,'amopen: bogus fstat, unit ',i4)
 30001 format(/,2x,'amopen: bogus fform, unit',i4)
 30002 format(/,2x,'Unit ',i4,' Error on OPEN: ',a)
 30003 format(/,2x,'Unit ',i4,' Error on OPEN: ',a)
 30004 format(/,2x,'Unit ',i4,' Error on OPEN: ',a)

end subroutine amopen

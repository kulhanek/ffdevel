! ==============================================================================
! This file is part of FFDevel.
!    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
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

module ffdev_constants

use ffdev_sizes

implicit none

!===============================================================================

! INTERNAL UNITS are:
! time      - fs
! length    - A
! velocity  - ?????
! energy    - kcal/mol
! force     - kcal/(mol.A)
! gradient  - kcal/(mol.A)
! mass      - g/mol

!===============================================================================

! math part ------------------------------------------------------------------
real(DEVDP),parameter   :: DEV_PI       = 3.1415926535897932384626433832795d0
real(DEVDP),parameter   :: DEV_HPI      = 0.5d0 * DEV_PI
real(DEVDP),parameter   :: DEV_R2D      = 180.0d0 / DEV_PI
real(DEVDP),parameter   :: DEV_D2R      = DEV_PI / 180.0d0

! gas constant 8.314 472(15) J mol-1 K-1
real(DEVDP), parameter  :: DEV_Rgas     = 0.0019872065d0     ! kcal mol-1 K-1 = 8.314 472 / 4184

! mass * velocity^2 -> kcal/mol
! g/mol * A^2 / fs^2 = 0.001 kg/mol * (10^-10)^2 m^2 / ((10^-15)^s^2)
! 0.001*10^-20/10^-30 * kg/mol*m^2/s^s = 10^-3*10^10 J/mol = 10^7 J/mol
! 10^7 J/mol = 10^7 / 4184 kcal/mol
! dt (fs) = vdt * sqrt( 4184 / 10^7) = 0.02045482828087295384
real(DEVDP), parameter  :: DEV_DT2VDT   = 0.02045482828087295384d0 ! sqrt(pc_cal/1d4)
real(DEVDP), parameter  :: DEV_VDT2DT   = 1.0d0 / DEV_DT2VDT

! hartree -> kcal/mol
real(DEVDP), parameter  :: DEV_HARTREE2KCL = 627.50960803059d0
real(DEVDP), parameter  :: DEV_KCL2HARTREE = 1.0d0 / DEV_HARTREE2KCL

! kJ/mol -> kcal/mol
real(DEVDP), parameter  :: DEV_KJ2KCL   = 1.0d0 / 4.184d0
real(DEVDP), parameter  :: DEV_KCL2KJ   = 1.0d0 / DEV_KJ2KCL

! eV -> kcal/mol
real(DEVDP), parameter  :: DEV_eV2KCL   = 23.06035d0
real(DEVDP), parameter  :: DEV_KCL2eV   = 1.0d0 / DEV_eV2KCL

! lambda (g A^2/(mol fs^2 CV) -> kcal/mol/CV
! 1000 * (10^10)^2 / (10^15)^2 = 1000 * 10^20 / 10^30 = 10^-7
real(DEVDP), parameter  :: DEV_L2CL     = 1e7 / 4184.d0
real(DEVDP), parameter  :: DEV_CL2L     = 1.0d0 / DEV_L2CL

! atomic unit time to fs (2.418 884 326 505(16)×10-17s)
real(DEVDP), parameter  :: DEV_AU2FS    = 2.418884326505d-2
real(DEVDP), parameter  :: DEV_FS2AU    = 1.0d0 / DEV_AU2FS

! ps to fs
real(DEVDP), parameter  :: DEV_PS2FS    = 1000.0d0
real(DEVDP), parameter  :: DEV_FS2PS    = 1.0d0 / DEV_PS2FS

! atomic unit length to A 5.291 772 108(18)×10-11m
real(DEVDP), parameter  :: DEV_AU2A     = 5.291772108d-1
real(DEVDP), parameter  :: DEV_A2AU     = 1.0d0 / DEV_AU2A

! pm to A
real(DEVDP), parameter  :: DEV_PM2A     = 0.01d0
real(DEVDP), parameter  :: DEV_A2PM     = 1.0d0 / DEV_PM2A

! nm to A
real(DEVDP), parameter  :: DEV_NM2A     = 10.0d0
real(DEVDP), parameter  :: DEV_A2NM     = 1.0d0 / DEV_NM2A

! libatoms mass to g/mol
! real(dp), parameter :: ELECTRONMASS_GPERMOL =  5.48579903e-4_dp !% grams/mol
! real(dp), parameter :: ELEM_CHARGE = 1.60217653e-19_dp !% coulombs
! real(dp), parameter :: HARTREE = 27.2113961_dp !% eV
! real(dp), parameter :: RYDBERG = 0.5_dp*HARTREE !% eV
! real(dp), parameter :: BOHR = 0.529177249_dp !% Angstrom
! real(dp), parameter :: HBAR_EVSEC = 6.5821220e-16_dp !% hbar in eV seconds
! real(dp), parameter :: HBAR_AU = 1.0_dp              !% hbar in a.u.
! real(dp), parameter :: HBAR = (HBAR_EVSEC*1e-15_dp)    !% hbar in eV fs
! real(dp), parameter :: ONESECOND = 1e15_dp           !% 1 second in fs
! real(dp), parameter :: ONESECOND_AU = (1.0_dp/(HBAR_EVSEC/(HBAR_AU*HARTREE))) !% 1 second in a.u.
! real(dp), parameter :: AU_FS = (1.0_dp/ONESECOND_AU*ONESECOND) !% a.u. time in fs
! real(dp), parameter :: MASSCONVERT = (1.0_dp/ELECTRONMASS_GPERMOL*HARTREE*AU_FS*AU_FS/(BOHR*BOHR))
! AU_FS= 1.0/1.0_dp/(6.5821220e-16/(1.0*27.2113961))*1e15
! AU_FS = (6.5821220e-16/(27.2113961))*1e15
! AU_FS = (6.5821220e-1/(27.2113961)) = 0.02418884343828283033 ----> DEV_AUT2FS
! 1.0/5.48579903e-4*27.2113961*0.02418884343828283033*0.02418884343828283033/(0.529177249*0.529177249)
! MASSCONVERT = 103.6427217590198535
real(DEVDP), parameter  :: DEV_AMU2LIBATOMSM = 103.6427217590198535d0
real(DEVDP), parameter  :: DEV_LIBATOMSM2AMU = 1.0d0 / DEV_AMU2LIBATOMSM

! a.u. to g/mol (amu)
real(DEVDP), parameter  :: DEV_AU2AMU = 1.0d0 / 1822.88842718d0
real(DEVDP), parameter  :: DEV_AMU2AU = 1.0d0 / DEV_AU2AMU

! common part ------------------------------------------------------------------
integer,parameter       :: DEV_STD_INPUT    = 5
integer,parameter       :: DEV_STD_OUTPUT   = 6
integer,parameter       :: DEV_TEST         = 1342
integer,parameter       :: DEV_DEBUG        = 1000
integer,parameter       :: DEV_XYZ          = 123
integer,parameter       :: DEV_TRAJ         = 125
integer,parameter       :: DEV_OTRAJ        = 126
integer,parameter       :: DEV_GEO          = 287
integer,parameter       :: DEV_TOP          = 289
integer,parameter       :: DEV_PRMS         = 298
integer,parameter       :: DEV_NULL         = 456
integer,parameter       :: DEV_PROFILE      = 371
integer,parameter       :: DEV_ERRSUMLOG    = 387
integer,parameter       :: DEV_NBPOT        = 149

! verbosity levels
integer,parameter       :: DEV_VERBOSITY_MINIMAL = 0
integer,parameter       :: DEV_VERBOSITY_MEDIUM  = 1
integer,parameter       :: DEV_VERBOSITY_FULL    = 2

!===============================================================================

end module ffdev_constants

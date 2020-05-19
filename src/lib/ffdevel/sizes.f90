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

module ffdev_sizes

implicit none

! ------------------------------------------------------------------------------
integer, parameter  :: MAX_PATH      = 255  ! max length of file names
integer, parameter  :: MAX_TNAME     =   4  ! max length of type name
integer, parameter  :: MAX_RNAME     =   4  ! max length of residue name
integer, parameter  :: MAX_CVTYPE    =   5  ! max length of CV type

! real numbers -----------------------------------------------------------------
integer, parameter  :: DEVDP         = 8

! ------------------------------------------------------------------------------

end module ffdev_sizes


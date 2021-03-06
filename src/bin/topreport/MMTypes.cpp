// =============================================================================
// This file is part of FFDevel.
//    Copyright (C) 2013 Petr Kulhanek, kulhanek@chemi.muni.cz
//
// FFDevel is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// FFDevel is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FFDevel. If not, see <http://www.gnu.org/licenses/>.
// =============================================================================

#include "MMTypes.hpp"
#include <math.h>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAtomType::CAtomType(void)
{
    idx = -1;
    mass = 0;
    z = -1;
    eps = 0.0;
    r0 = 0.0;
}

//------------------------------------------------------------------------------

bool CAtomType::operator != (const CAtomType& right)
{
    if( idx != right.idx ) return(true);
    if( z != right.z ) return(true);
    if( fabs(mass - right.mass) > 0.0001 ) return(true);

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CDihedralType::CDihedralType(void)
{
    idx = -1;
    at1 = -1;
    at2 = -1;
    at3 = -1;
    at4 = -1;
    scee = 0.0;
    scnb = 0.0;
}

//------------------------------------------------------------------------------

void CDihedralType::SetSeriesSize(int size)
{
    defined.resize(size);
    v0.resize(size);
    phase.resize(size);
}

//------------------------------------------------------------------------------

int CDihedralType::GetSeriesSize(void)
{
    return(v0.size());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



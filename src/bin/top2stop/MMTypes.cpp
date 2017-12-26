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

CBondType::CBondType(void)
{
    idx = -1;
    at1 = -1;
    at2 = -1;
    form = -1;
    d0 = 0;
    k = 0;
}

//------------------------------------------------------------------------------

bool CBondType::operator != (const CBondType& right)
{
    if( idx != right.idx ) return(true);

    if( ! ( ((at1 == right.at1)&&(at2 == right.at2)) ||
            ((at1 == right.at2)&&(at2 == right.at1)) )  ) return(true);

    if( form != right.form ) return(true);
    if( fabs(d0 - right.d0) > 0.0001 ) return(true);
    if( fabs(k - right.k) > 0.0001 ) return(true);

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAngleType::CAngleType(void)
{
    idx = -1;
    at1 = -1;
    at2 = -1;
    at3 = -1;
    form = -1;
    a0 = 0;
    k = 0;
}

//------------------------------------------------------------------------------

bool CAngleType::operator != (const CAngleType& right)
{
    if( idx != right.idx ) return(true);

    if( ! ( ((at1 == right.at1)&&(at2 == right.at2)&&(at3 == right.at3)) ||
            ((at1 == right.at3)&&(at2 == right.at2)&&(at3 == right.at1)) )  ) return(true);

    if( form != right.form ) return(true);
    if( fabs(a0 - right.a0) > 0.0001 ) return(true);
    if( fabs(k - right.k) > 0.0001 ) return(true);

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
    nb_processed = false;
}

//------------------------------------------------------------------------------

void CDihedralType::SetSeriesSize(int size)
{
    defined.resize(size);
    v0.resize(size);
    phase.resize(size);
    c.resize(size);
    p.resize(size);
    w2.resize(size);
    if( size > 2) {
        p[0] = -M_PI;
        p[size-1] = M_PI;
        for(int i=1; i < size-1; i++){
            p[i] = -M_PI + 2.0*M_PI*i/(size-1);
        }
        for(int i=0; i < size; i++){
            w2[i] = pow(2.0*M_PI/(size-1),2.0);
        }
    }

}

//------------------------------------------------------------------------------

int CDihedralType::GetSeriesSize(void)
{
    return(v0.size());
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



// =============================================================================
// NEMESIS - Molecular Modelling Package
// -----------------------------------------------------------------------------
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu,
//                       Jakub Stepan, xstepan3@chemi.muni.cz
//    Copyright (C) 1998-2004 Petr Kulhanek, kulhanek@chemi.muni.cz
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// =============================================================================

#include "Sphere.hpp"

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CSphere::CSphere(void)
{
    Tessellation = 0;
    SetTessellationQuality(4);
}

//---------------------------------------------------------------------------

bool CSphere::SetTessellationQuality(const unsigned int quality)
{
    if( quality == Tessellation ) return(true);
    if( (quality < 1 ) || (quality>8) )return(false);

    Tessellation = quality;

    ComputeVertices();

    return(true);
}

//---------------------------------------------------------------------------

unsigned int CSphere::GetNumberOfVertices(void) const
{
    return(Vertices.size());
}

//---------------------------------------------------------------------------

const CPoint CSphere::GetVertice(unsigned int i,const CPoint& origin,double r1)
{
    return(Vertices[i]*r1 + origin);
}

//---------------------------------------------------------------------------

void CSphere::ComputeVertices(void)
{
    Vertices.clear();

    double x = 0.525731112119133606f, z = 0.850650808352039932f;

    static double vdata[12][3] = {
        {-x, 0, z}, {x, 0, z}, {-x, 0, -z}, {x, 0, -z},
        {0, z, x}, {0, z, -x}, {0, -z, x}, {0, -z, -x},
        {z, x, 0}, {-z, x, 0}, {z, -x, 0}, {-z, -x, 0}
    };

    static int tindices[20][3] = {
        {0, 4, 1}, {0, 9, 4}, {9, 5, 4}, {4, 5, 8}, {4, 8, 1},
        {8, 10, 1}, {8, 3, 10}, {5, 3, 8}, {5, 2, 3}, {2, 7, 3},
        {7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6},
        {6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5}, {7, 2, 11}
    };

    for (unsigned int i = 0; i < 20; i++){
        ComputePartition(CPoint(vdata[tindices[i][0]][0],vdata[tindices[i][0]][1],vdata[tindices[i][0]][2]),
                         CPoint(vdata[tindices[i][1]][0],vdata[tindices[i][1]][1],vdata[tindices[i][1]][2]),
                         CPoint(vdata[tindices[i][2]][0],vdata[tindices[i][2]][1],vdata[tindices[i][2]][2]),Tessellation - 1);
    }
}

//---------------------------------------------------------------------------

void CSphere::ComputePartition(const CPoint& v1,
                               const CPoint& v2,
                               const CPoint& v3,
                               const unsigned int complexity)
{
    if (complexity == 0) {
        Vertices.push_back(v1);
        Vertices.push_back(v2);
        Vertices.push_back(v3);
    } else {
        CPoint v12,v23,v31;
        v12 = v1 + v2;
        v23 = v2 + v3;
        v31 = v3 + v1;
        v12.Normalize();
        v23.Normalize();
        v31.Normalize();
        ComputePartition(v1,v12,v31,complexity - 1);
        ComputePartition(v2,v23,v12,complexity - 1);
        ComputePartition(v3,v31,v23,complexity - 1);
        ComputePartition(v12,v23,v31,complexity - 1);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


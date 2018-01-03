#ifndef SphereH
#define SphereH
// =============================================================================
// NEMESIS - Molecular Modelling Package
// -----------------------------------------------------------------------------
//    Copyright (C) 2018 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <vector>
#include <Point.hpp>

// -----------------------------------------------------------------------------

/// sphere

class CSphere {
public:
    CSphere(void);   // default tessellation=3

    bool            SetTessellationQuality(const unsigned int quality);
    unsigned int    GetNumberOfVertices(void) const;
    const CPoint    GetVertice(unsigned int i,const CPoint& origin,double r1);

// section of private data ----------------------------------------------------
private:
    unsigned int        Tessellation;  // quality
    std::vector<CPoint> Vertices;

    void ComputeVertices(void);
    void ComputePartition(const CPoint& v1,const CPoint& v2,const CPoint& v3,
                          const unsigned int complexity);
};

// -----------------------------------------------------------------------------

#endif


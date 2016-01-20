#ifndef MMTypesH
#define MMTypesH
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

#include <string>
#include <vector>

//------------------------------------------------------------------------------

class CAtomType {
public:
    int         idx;
    std::string name;
    double      mass;
    int         z;
    double      eps;
    double      r0;
public:
    CAtomType(void);
    bool operator != (const CAtomType& right);
};

//------------------------------------------------------------------------------

class CDihedralType {
public:
    int         idx;
    int         at1;
    int         at2;
    int         at3;
    int         at4;
    std::vector<bool>   defined;
    std::vector<double> v0;
    std::vector<double> phase;
    double      scee;
    double      scnb;
public:
    CDihedralType(void);
    void SetSeriesSize(int size);
    int  GetSeriesSize(void);
};

//------------------------------------------------------------------------------

#endif

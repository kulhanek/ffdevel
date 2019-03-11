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
public:
    CAtomType(void);
    bool operator != (const CAtomType& right);
};

//------------------------------------------------------------------------------

class CBondType {
public:
    int         idx;
    int         at1;
    int         at2;
    int         form;
    double      d0;
    double      k;
public:
    CBondType(void);
    bool operator != (const CBondType& right);
};

//------------------------------------------------------------------------------

class CAngleType {
public:
    int         idx;
    int         at1;
    int         at2;
    int         at3;
    int         form;
    double      a0;
    double      k;
public:
    CAngleType(void);
    bool operator != (const CAngleType& right);
};

//------------------------------------------------------------------------------

class CDihedralType {
public:
    int                 idx;
    int                 at1;
    int                 at2;
    int                 at3;
    int                 at4;
    std::vector<bool>   defined;
    bool                grbf;
    // cos series
    std::vector<double> v0;
    std::vector<double> phase;
    // gaussian radial basis functions
    std::vector<double> c;
    std::vector<double> p;
    std::vector<double> w2;
    double              scee;
    double              scnb;
    bool                nb_processed;
public:
    CDihedralType(void);
    void    SetSeriesSize(int size);
    int     GetSeriesSize(void);
    void    Cos2GRBF(int dih_samp_freq);
    double  RMSECos2GRBF(int dih_samp_freq);
    double  GetCOSValue(double x);
    double  GetGRBFValue(double x);
    double  GetDihDeviation(double value1, double value2);
};

//------------------------------------------------------------------------------

class CDihedral{
public:
    int         at1;
    int         at2;
    int         at3;
    int         at4;
};

//------------------------------------------------------------------------------

class CDihedralTypeFilter {
public:
    std::string t1;
    std::string t2;
    std::string t3;
    std::string t4;
    bool        full;

    CDihedralTypeFilter(void);
};

//------------------------------------------------------------------------------

#endif

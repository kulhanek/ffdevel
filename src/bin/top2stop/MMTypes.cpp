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
#include <SimpleVector.hpp>
#include <SciLapack.hpp>

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
    grbf = false;
    DihCOffset = false;
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
    
    for(int i=0; i < size; i++){
        defined[i] = false;
        v0[i] = 0.0;
        phase[i] = 0.0;
        c[i] = 0.0;
        p[i] = -M_PI + 2.0*M_PI*i/size; 
        w2[i] = pow(2.0*M_PI/size,2.0);
    }
}

//------------------------------------------------------------------------------

int CDihedralType::GetSeriesSize(void)
{
    return(v0.size());
}

//------------------------------------------------------------------------------

void CDihedralType::Cos2GRBF(int dih_samp_freq)
{
    CVector         rhs;
    CFortranMatrix  A;

    A.CreateMatrix(GetSeriesSize(),GetSeriesSize());
    rhs.CreateVector(GetSeriesSize());
    
    // central matrix
    for(int i=0; i < GetSeriesSize(); i++){
        double sp1 = p[i];
        double sw1 = w2[i];
        for(int j=0; j < GetSeriesSize(); j++){
            double sp2 = p[j];
            double sw2 = w2[j];
            double a = 0;
            for(int k=0; k < (GetSeriesSize()+1)*dih_samp_freq; k++){
                double x = -M_PI + 2.0*M_PI*k/((GetSeriesSize()+1)*dih_samp_freq);
                a += exp(-(GetDihDeviation(x,sp1))*(GetDihDeviation(x,sp1))/sw1)*exp(-(GetDihDeviation(x,sp2))*(GetDihDeviation(x,sp2))/sw2);
            }
            A[i][j] = a;
        }
    }

    // rhs
    for(int i=0; i < GetSeriesSize(); i++){
        double p1 = p[i];
        double w1 = w2[i];
        double a = 0;
        for(int k=0; k < (GetSeriesSize()+1)*dih_samp_freq; k++){
            double x = -M_PI + 2.0*M_PI*k/((GetSeriesSize()+1)*dih_samp_freq);
            double value = GetCOSValue(x);
            a += exp(-(GetDihDeviation(x,p1))*(GetDihDeviation(x,p1))/w1)*value;
        }
        rhs[i] = a;
        // std::cout << p[i]*180.0/M_PI << " " << sqrt(w2[i])*180.0/M_PI << " " << a << std::endl;          
    }
    
    // solv equations
    if( CSciLapack::solvleLU(A,rhs) != 0 ){
        RUNTIME_ERROR("unable to solve transformation");
    }

    // copy final result
    for(int l=0; l < GetSeriesSize(); l++){
        c[l] = rhs[l];
    }

}
    
//------------------------------------------------------------------------------    
    
double CDihedralType::RMSECos2GRBF(int dih_samp_freq)
{    
    // calculate rmse
    double rmse = 0.0;
    for(int k=0; k < (GetSeriesSize()+1)*dih_samp_freq; k++){
        double x = -M_PI + 2.0*M_PI*k/((GetSeriesSize()+1)*dih_samp_freq);
        double value1 = GetCOSValue(x);
        double value2 = GetGRBFValue(x);
        double error = value2-value1;
        rmse += error*error;
    }
    if( (GetSeriesSize()+1)*dih_samp_freq > 0 ) {
        rmse /= (GetSeriesSize()+1)*dih_samp_freq;
    }
    rmse = sqrt(rmse);
    return(rmse);
}

//------------------------------------------------------------------------------  

double  CDihedralType::GetCOSValue(double x)
{
    double value1 = 0;
    for(int l=1; l <= GetSeriesSize(); l++){
        double arg = l*x - phase[l-1];
        value1 += v0[l-1]*(1.0+cos(arg));
    }
    value1 = DihCOffset + value1;
    return(value1);
}

//------------------------------------------------------------------------------  

double  CDihedralType::GetGRBFValue(double x)
{
    double value2 = 0;
    for(int l=0; l < GetSeriesSize(); l++){
        double c1 = c[l];
        double p1 = p[l];
        double w1 = w2[l];
        value2 += c1*exp(-(GetDihDeviation(x,p1))*(GetDihDeviation(x,p1))/w1);
    } 
    return(value2);
}

//------------------------------------------------------------------------------  

double CDihedralType::GetDihDeviation(double value1, double value2)
{
       
    double minv,maxv,vec;

    minv = -M_PI;
    maxv =  M_PI;

    if( fabs(value1-value2) <  0.5*(maxv-minv) ) {
        return(value1 - value2);
    } else {
        //! get vector
        vec = value1 - value2;
        //! shift to box center
        vec = vec + 0.5*(maxv+minv);
        //! image as point
        vec = vec - (maxv-minv)*floor((vec-minv)/(maxv-minv));
        //! return vector back
        return(vec - 0.5*(maxv+minv));
    }
}

//------------------------------------------------------------------------------

CDihedralTypeFilter::CDihedralTypeFilter(void)
{
    full = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



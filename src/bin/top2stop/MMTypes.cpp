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
    IAC = -1;
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
    mode = EDM_COS;
    DihCOffset = false;
}

//------------------------------------------------------------------------------

void CDihedralType::SetSeriesSize(int nsize,double wfac)
{
    defined.resize(nsize);
    v0.resize(nsize);
    phase.resize(nsize);
    c.resize(nsize);
    p.resize(nsize);
    w2.resize(nsize);

    for(int i=0; i < nsize; i++){
        defined[i] = false;
        v0[i] = 0.0;
        phase[i] = 0.0;
        c[i] = 0.0;
        p[i] = -M_PI + 2.0*M_PI*i/nsize;
        w2[i] = pow(wfac*2.0*M_PI/nsize,2.0);
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

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CDihedralType::Cos2GRBF(int nsamples)
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
            for(int k=0; k < nsamples; k++){
                double x = -M_PI + 2.0*M_PI*k/nsamples;
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
        for(int k=0; k < nsamples; k++){
            double x = -M_PI + 2.0*M_PI*k/nsamples;
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

double CDihedralType::RMSECos2GRBF(int nsamples)
{
    // calculate rmse
    double rmse = 0.0;
    for(int k=0; k < nsamples; k++){
        double x = -M_PI + 2.0*M_PI*k/nsamples;
        double value1 = GetCOSValue(x);
        double value2 = GetGRBFValue(x);
        double error = value2-value1;
        rmse += error*error;
    }
    if( nsamples > 0 ) {
        rmse /= nsamples;
    }
    rmse = sqrt(rmse);
    return(rmse);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

static inline int WrapIndex(int i, int n)
{
    i %= n;
    if( i < 0 ) i += n;
    return i;
}

//==============================================================================
// Solve dense linear system A x = b by Gaussian elimination with pivoting
// Used here for the periodic spline coefficient setup.
//==============================================================================

static std::vector<double> SolveLinearSystem(
                                    std::vector<std::vector<double>> A,
                                    std::vector<double> b)
{
    const int n = static_cast<int>(b.size());

    for( int k = 0; k < n; k++ ) {
        // pivot
        int piv = k;
        double amax = std::fabs(A[k][k]);
        for( int i = k+1; i < n; i++ ) {
            if( std::fabs(A[i][k]) > amax ) {
                amax = std::fabs(A[i][k]);
                piv = i;
            }
        }

        if( amax < 1e-14 ) {
            throw std::runtime_error("Singular matrix in periodic B-spline interpolation.");
        }

        if( piv != k ) {
            std::swap(A[piv], A[k]);
            std::swap(b[piv], b[k]);
        }

        // elimination
        for( int i = k+1; i < n; i++ ) {
            const double fac = A[i][k] / A[k][k];
            A[i][k] = 0.0;
            for( int j = k+1; j < n; j++ ) {
                A[i][j] -= fac * A[k][j];
            }
            b[i] -= fac * b[k];
        }
    }

    // back substitution
    std::vector<double> x(n,0.0);
    for( int i = n-1; i >= 0; i-- ) {
        double sum = b[i];
        for( int j = i+1; j < n; j++ ) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }

    return x;
}

//------------------------------------------------------------------------------

void CDihedralType::Cos2CBS(void)
{
    const int n = GetSeriesSize();

    std::vector<double> a;
    a.resize(n, 0.0);

    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::vector<double> rhs(n, 0.0);

    for( int i = 0; i < n; i++ ) {
        const int im1 = WrapIndex(i-1, n);
        const int ip1 = WrapIndex(i+1, n);

        A[i][im1] = 1.0 / 6.0;
        A[i][i  ] = 4.0 / 6.0;
        A[i][ip1] = 1.0 / 6.0;

        rhs[i] = GetCOSValue(p[i]);
    }

    a = SolveLinearSystem(A, rhs);

    for( int i = 0; i < n; i++ ) {
        c[i] = a[i];
    }
}

//------------------------------------------------------------------------------

// wrap angle to [-pi,pi)
static inline double WrapAnglePM(double x)
{
    const double twopi = 2.0 * M_PI;
    x = std::fmod(x + M_PI, twopi);
    if( x < 0.0 ) x += twopi;
    return x - M_PI;
}

//------------------------------------------------------------------------------

// shortest periodic distance in "grid units"
static inline double WrapPeriodicT(double t, int n)
{
    t -= std::round(t / n) * n;
    return t;
}

//------------------------------------------------------------------------------

static inline double CBS_B3(double t)
{
    t = std::fabs(t);

    if( t < 1.0 ) {
        return (4.0 - 6.0*t*t + 3.0*t*t*t) / 6.0;
    } else if( t < 2.0 ) {
        const double u = 2.0 - t;
        return (u*u*u) / 6.0;
    } else {
        return 0.0;
    }
}

//------------------------------------------------------------------------------

double CDihedralType::GetCBSValue(double x)
{
    const int n = GetSeriesSize();
    if( n < 3 ) return 0.0;

    const double h = 2.0 * M_PI / n;

    // periodic wrap to [-pi,pi)
    x = WrapAnglePM(x);

    double value = 0.0;

    // Only nearby 4 basis functions contribute, but for clarity
    // we keep a full loop. This is fine for small n.
    for( int i = 0; i < n; i++ ) {
        double t = (x - p[i]) / h;
        t = WrapPeriodicT(t, n);
        value += c[i] * CBS_B3(t);
    }

    return value;
}

//------------------------------------------------------------------------------

double  CDihedralType::RMSECos2CBS(int nsamples)
{
    // calculate rmse
    double rmse = 0.0;
    for(int k=0; k < nsamples; k++){
        double x = -M_PI + 2.0*M_PI*k/nsamples;
        double value1 = GetCOSValue(x);
        double value2 = GetCBSValue(x);
        double error = value2-value1;
        rmse += error*error;
    }
    if( nsamples > 0 ) {
        rmse /= nsamples;
    }
    rmse = sqrt(rmse);
    return(rmse);
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



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

#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <SimpleList.hpp>
#include <ctype.h>
#include "Cos2GRBF.hpp"
#include <map>
#include <iomanip>
#include <list>
#include <locale>
#include <FortranMatrix.hpp>
#include <Vector.hpp>
#include <Lapack.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

MAIN_ENTRY(CCos2GRBF);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCos2GRBF::CCos2GRBF(void)
{

}

//------------------------------------------------------------------------------

CCos2GRBF::~CCos2GRBF(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CCos2GRBF::Init(int argc,char* argv[])
{
    // encode program options
    int result = Options.ParseCmdLine(argc,argv);

    // should we exit or was it error?
    if( result != SO_CONTINUE ) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# cos2grbf (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    vout << "# Dihedral definition            : " << Options.GetArgInput() << endl;
    vout << "# Output file                    : " << Options.GetArgOutput() << endl;
    vout << "# Dihedral sequence size         : " << Options.GetOptDihedralSeriesSize() << endl;
    vout << "# Dihedral segment sampling      : " << Options.GetOptDihedralSamplingSize() << endl;

    return( SO_CONTINUE );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CCos2GRBF::Run(void)
{
    // read input
    if( ReadFrcMod() == false ) return(false); 

    // convert
    DihedralType.Cos2GRBF(Options.GetOptDihedralSamplingSize());
    
    // write output
     if( WriteOutput() == false ) return(false); 
    
    return(true);
}

//------------------------------------------------------------------------------

bool CCos2GRBF::ReadFrcMod(void)
{
    ifstream ifs(Options.GetArgInput());
    if( !ifs ) {
        ES_ERROR("unable to open input file");
        return(false);
    }
    
    DihedralType.SetSeriesSize(Options.GetOptDihedralSeriesSize());
    
    string line;
    while( getline(ifs,line) ){
        stringstream str(line);
        string  buf;
        int     num,n;
        double  v,p;
        str >> buf >> num >> v >> p >> n;
        if( str ){
            n = abs(n) - 1;
            if( (n < 0) || (n >= Options.GetOptDihedralSeriesSize()) ) {
                CSmallString error;
                error << "n is out of range!";
                ES_ERROR(error);
                return(false);
            }
            DihedralType.phase[n] = p*M_PI/180.0;
            DihedralType.v0[n] = v;
            DihedralType.defined[n] = true;
        }
    }
    
    return(true);
}

//------------------------------------------------------------------------------

bool CCos2GRBF::WriteOutput(void)
{
    ofstream ofs(Options.GetArgOutput());
    if( !ofs ) {
        ES_ERROR("unable to open output file");
        return(false);
    }
    
    ofs << "# original cos series ..." << endl;
    ofs << "#  N     V0         Phase      " << endl;
    ofs << "# --- ------------ ------------" << endl;
    for(int k=0; k < DihedralType.GetSeriesSize(); k++){
        if(  DihedralType.defined[k] ){
            ofs << format("# %3d %12.3f %12.3f\n")%(k+1)%DihedralType.v0[k]%(DihedralType.phase[k]*180.0/M_PI);
        }
    }
    ofs << endl;
    
    
    double rmse = DihedralType.RMSECos2GRBF(Options.GetOptDihedralSamplingSize());
    
    ofs << format("# rmse = %12.6f\n")%rmse;
    ofs << "#     phi           cos        grbf          diff    " << endl;
    ofs << "# ------------ ------------ ------------ ------------" << endl;
    
    for(int k=0; k <= (DihedralType.GetSeriesSize()+1)*Options.GetOptDihedralSamplingSize(); k++){
        double x = -M_PI + 2.0*M_PI*k/((DihedralType.GetSeriesSize()+1)*Options.GetOptDihedralSamplingSize());
        double value1 = DihedralType.GetCOSValue(x);
        double value2 = DihedralType.GetGRBFValue(x);
        ofs << format("  %12.3f %12.6f %12.6f %12.6f\n")%(x*180.0/M_PI)%value1%value2%(value2-value1);
    }
    
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCos2GRBF::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# cos2grbf terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



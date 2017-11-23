// =============================================================================
// This file is part of FFDevel.
//    Copyright (C) 2017 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include "VDWGenProbe.hpp"
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <set>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

MAIN_ENTRY(CVDWGenProbe);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CVDWGenProbe::CVDWGenProbe(void)
{
    OutputFile = NULL;
}

//------------------------------------------------------------------------------

CVDWGenProbe::~CVDWGenProbe(void)
{
    // close the output file
    if(Options.GetArgStructureName() != "-"){
        if( OutputFile != NULL ) fclose(OutputFile);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CVDWGenProbe::Init(int argc,char* argv[])
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
    vout << "# vdwgenprobe (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    return( SO_CONTINUE );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool comp_d1(const SProbeData left, const SProbeData right)
{
    return( left.D1 < right.D1 );
}

bool comp_d2(const SProbeData& left, const SProbeData& right)
{
    return( left.D2 > right.D2 );
}


//------------------------------------------------------------------------------

bool CVDWGenProbe::Run(void)
{
    // load structure
    if( LoadStructure() == false ) return(false);

    // open output file
    if(Options.GetArgOutputName() != "-") {
        if( Options.GetOptAppend() == true ){
            OutputFile = fopen(Options.GetArgOutputName(),"a");
        } else {
            OutputFile = fopen(Options.GetArgOutputName(),"w");
        }
    } else {
        OutputFile = stdout;
    }
    if( OutputFile == NULL ){
        vout << "<red>>>> ERROR: Unable to open the output file: " << Options.GetArgOutputName() << "</red>" << endl;
        return(false);
    }

    // create output structure
    StructureWithProbe.SetNumberOfAtoms(Structure.GetNumberOfAtoms()+1);
    // copy basal atoms
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        StructureWithProbe.SetPosition(i,Structure.GetPosition(i));
        StructureWithProbe.SetSymbol(i,Structure.GetSymbol(i));
    }
    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),CPoint());
    StructureWithProbe.SetSymbol(Structure.GetNumberOfAtoms(),Options.GetOptProbeSymbol());

    // generate structures
    GenAllStructures();

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GenAllStructures(void)
{
    CPoint probe;
    for(probe.x = Min.x; probe.x <= Max.x; probe.x += Options.GetOptSpacing()){
        cout << "x";
        for(probe.y = Min.y; probe.y <= Max.y; probe.y += Options.GetOptSpacing()){
            cout << "y";
            for(probe.z = Min.z; probe.z <= Max.z; probe.z += Options.GetOptSpacing()){
                if( CheckProbe(probe) ){
                    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);

                    if( StructureWithProbe.Save(OutputFile) == false ){
                        vout << "<red>>>> ERROR: Unable to save a structure to the output file: " << Options.GetArgOutputName() << "</red>" << endl;
                        return(false);
                    }
                }
            }
        }
    }
    cout << endl;
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::CheckProbe(const CPoint& probe)
{
    bool ok = false;
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        CPoint apos = Structure.GetPosition(i);
        double d = Size(apos - probe);
        if( d < Options.GetOptRMin() ) return(false);
        if( (d <= Options.GetOptRMax()) && SelectedAtoms[i] ) ok = true;
    }

    return(ok);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CVDWGenProbe::Finalize(void)
{    
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << "# ==============================================================================" << endl;
    vout << "# vdwgenprobe terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CVDWGenProbe::LoadStructure(void)
{
    if(Options.GetArgStructureName() != "-") {
        vout << "# Input XYZ structure (in)       : " << Options.GetArgStructureName() << endl;
    } else {
        vout << "# Input XYZ structure (in)       : - (standard input)" << endl;
    }

    if( Structure.Load(Options.GetArgStructureName()) == false ) {
        vout << "<red>>>> ERROR: Unable to load the XYZ structure: " << Options.GetArgStructureName() << "</red>" << endl;
        return(false);
    }
    vout << "# Number of atoms                = " << Structure.GetNumberOfAtoms() <<  endl;

    // determine min/max
    for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
        CPoint at = Structure.GetPosition(i);
        if( i == 0 ){
            Min = at;
            Max = at;
        }
        if( Min.x > at.x ) Min.x = at.x;
        if( Min.y > at.y ) Min.y = at.y;
        if( Min.z > at.z ) Min.z = at.z;

        if( Max.x < at.x ) Max.x = at.x;
        if( Max.y < at.y ) Max.y = at.y;
        if( Max.z < at.z ) Max.z = at.z;
    }

    SelectedAtoms.resize(Structure.GetNumberOfAtoms());

    int sat = 0;
    // gen list of selected atoms
    if( Options.GetOptSelectedAtoms() == NULL ){
        // select all atoms
        for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
            SelectedAtoms[i] = true;
        }
        sat = Structure.GetNumberOfAtoms();
    } else {
        string          slist(Options.GetOptSelectedAtoms());
        vector<string>  idxs;
        split(idxs,slist,is_any_of(","));
        for(int i=0; i < idxs.size(); i++){
            int indx = 0;
            stringstream str(idxs[i]);
            str >> indx;
            indx--;
            if( (indx >= 0) && (indx < Structure.GetNumberOfAtoms()) ){
                if( SelectedAtoms[indx] == false ) sat++;
                SelectedAtoms[indx] = true;
            }
        }
    }
    vout << "# Number of selected atoms       = " << sat <<  endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



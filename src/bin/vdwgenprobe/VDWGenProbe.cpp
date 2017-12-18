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
#include <openbabel/mol.h>
#include <openbabel/graphsym.h>
#include <openbabel/obconversion.h>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace OpenBabel;

MAIN_ENTRY(CVDWGenProbe);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CVDWGenProbe::CVDWGenProbe(void)
{
    OutputFile = NULL;
    NumOfProbes = 0;
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

    // load probes
    if( Options.IsOptMSMSSurfaceVerticesSet() ){
        if( LoadMSMSProbes() == false ) return(false);
    }

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


    // decode filters
    string sfilters = string(Options.GetOptFilters());
    vector<string> filters;
    split(filters,sfilters,is_any_of(","));
    vector<string>::iterator it = filters.begin();
    vector<string>::iterator ie = filters.end();
    while( it != ie ){
        string filter = *it;
        if( filter == "minmax" ){
            Filters.push_back(EF_MINMAX);
        } else if( filter == "minonly" ){
            Filters.push_back(EF_MINONLY);
        } else if( filter == "keepall" ){
            Filters.push_back(EF_KEEPALL);
        } else {
            vout << "<red>>>> ERROR: Unsupported filter: " << filter << "</red>" << endl;
            return(false);
        }
        it++;
    }

    // generate structures
    bool result;
    if( Options.GetOptGenerator() == "grid" ){
        result = GeneratorGrid();
    } else if( Options.GetOptGenerator() == "random" ){
        result = GeneratorRandom();
    } else if( Options.GetOptGenerator() == "msms-all" ){
        result = GeneratorMSMSAll();
    } else if( Options.GetOptGenerator() == "msms-random" ){
        result = GeneratorMSMSRandom();
    } else {
        vout << "<red>>>> ERROR: Unsupported generator: " << Options.GetOptGenerator() << "</red>" << endl;
        return(false);
    }

    // final stat
        vout << "# Number of generated probes     = " << NumOfProbes <<  endl;
    return(result);
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
        vout << "# Input XYZ structure (in)       = " << Options.GetArgStructureName() << endl;
    } else {
        vout << "# Input XYZ structure (in)       = - (standard input)" << endl;
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

    double buffer = Options.GetOptBuffer();
    if( Options.GetOptBuffer() == -1 ){
        buffer = Options.GetOptRMax();
    }
    Min.x -= buffer;
    Min.y -= buffer;
    Min.z -= buffer;
    Max.x += buffer;
    Max.y += buffer;
    Max.z += buffer;

// generate selected atoms
    SelectedAtoms.resize(Structure.GetNumberOfAtoms());

    vout << "# Selected atoms                 = " << Options.GetOptSelectedAtoms() <<  endl;

    int sat = 0;
    // gen list of selected atoms
    if( Options.GetOptSelectedAtoms() == "all" ){
        // select all atoms
        for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
            SelectedAtoms[i] = true;
            SelectedAtomIds.insert(i);
        }
        sat = Structure.GetNumberOfAtoms();
    } else if( Options.GetOptSelectedAtoms() == "symclasses" ){
        // read structure as OBMol
        ifstream sin;
        sin.open(Options.GetArgStructureName());
        if( !sin ) {
            vout << "<red>>>> ERROR: Unable to load (as OBMol) the XYZ structure: " << Options.GetArgStructureName() << "</red>" << endl;
            return(false);
        }

        OBConversion    conv(&sin, &cout);
        OBFormat*       obFormat = conv.FormatFromExt(Options.GetArgStructureName());
        OBMol           mol;

        if( obFormat == NULL ) {
            vout << "<red>>>> ERROR: Unsupported input file format!</red>" << endl;
            return(false);
        }

        if( ! conv.SetInFormat(obFormat) ) {
            vout << "<red>>>> ERROR: Unsupported input file format!</red>" << endl;
            return(false);
        }

        if( ! conv.Read(&mol) ) {
            vout << "<red>>>> ERROR: Unable to load obmol from the stream!</red>" << endl;
            return(false);
        }

        mol.ConnectTheDots();
        mol.PerceiveBondOrders();

        // generate symmetry classes
        OBGraphSym  graph_sym(&mol);
        std::vector<unsigned int> symclasses;
        graph_sym.GetSymmetry(symclasses);

        // select unique atoms
        sat = 0;
        std::set<unsigned int>  selected_classes;
        CSmallString sats;
        for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
            unsigned int sc = symclasses[i];
            if( selected_classes.count(sc) != 1 ){
                SelectedAtoms[i] = true;
                selected_classes.insert(sc);
                SelectedAtomIds.insert(i);
                if( sats != NULL ) sats << ",";
                sats << i+1;
                sat++;
            }
        }
        vout << "# Selected atoms                 = " << sats <<  endl;
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
                SelectedAtomIds.insert(indx);
            }
        }
    }
    vout << "# Number of selected atoms       = " << sat <<  endl;

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::LoadMSMSProbes(void)
{
        vout << "# MSMS surface vertices (in)     = " << Options.GetOptMSMSSurfaceVertices() << endl;
    ifstream ifs(Options.GetOptMSMSSurfaceVertices());
    if( !ifs ){
        vout << "<red>>>> ERROR: Unable to load the file with vertices: " << Options.GetOptMSMSSurfaceVertices() << "</red>" << endl;
        return(false);
    }

    string line;

    // skip three lines
    getline(ifs,line);
    getline(ifs,line);
    getline(ifs,line);

    while( getline(ifs,line) ){
        // get xyz
        stringstream sfs(line);
        CProbe probe;
        probe.AtomId = -1;
        probe.Selected = true;
        sfs >> probe.Pos.x >> probe.Pos.y >> probe.Pos.z;
        MSMSProbes.push_back(probe);
    }

    vout << "# Number of vertices             = " << MSMSProbes.size() <<  endl;

    // find closest atoms

    CProbe probe;
    for(size_t i=0; i < MSMSProbes.size(); i++){
        probe = MSMSProbes[i];
        double min;
        int minid = -1;
        for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
            CPoint apos = Structure.GetPosition(j);
            double d = Square(apos - probe.Pos);
            if( (j == 0) || (d < min) ){
                min = d;
                minid = j;
            }
        }
        MSMSProbes[i].AtomId = minid;
    }

    // filter by selected atoms
    int sel = 0;
    for(size_t i=0; i < MSMSProbes.size(); i++){
        if( SelectedAtomIds.count(MSMSProbes[i].AtomId) == 1 ){
            sel++;
            MSMSProbes[i].Selected = true;
        } else {
            MSMSProbes[i].Selected = false;
        }
    }

    // stat
    std::set<int>::iterator  it = SelectedAtomIds.begin();
    std::set<int>::iterator  ie = SelectedAtomIds.end();

    vout << "# Number of selected vertices    = " << sel <<  endl;
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::SaveStructure(void)
{
    NumOfProbes++;
    if( NumOfProbes == Options.GetOptMaxProbes() ){
        vout << "<red>>>> ERROR: MaxProbe limit was reached!</red>" << endl;
        return(false);
    }

    if( StructureWithProbe.Save(OutputFile) == false ){
        vout << "<red>>>> ERROR: Unable to save a structure to the output file: " << Options.GetArgOutputName() << "</red>" << endl;
        return(false);
    }
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CVDWGenProbe::GeneratorGrid(void)
{
    CPoint probe;
    for(probe.x = Min.x; probe.x <= Max.x; probe.x += Options.GetOptSpacing()){
        for(probe.y = Min.y; probe.y <= Max.y; probe.y += Options.GetOptSpacing()){
            for(probe.z = Min.z; probe.z <= Max.z; probe.z += Options.GetOptSpacing()){
                if( FilterProbe(probe) ){
                    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
                    if( SaveStructure() == false ) return(false);
                }
            }
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorRandom(void)
{
    CPoint probe;
    for(int i=0; i < Options.GetOptNumOfTrials(); i++){
        probe.x = rand()*(Max.x - Min.x)/RAND_MAX + Min.x;
        probe.y = rand()*(Max.y - Min.y)/RAND_MAX + Min.y;
        probe.z = rand()*(Max.z - Min.z)/RAND_MAX + Min.z;
        if( FilterProbe(probe) ){
            StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
            if( SaveStructure() == false ) return(false);
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorMSMSAll(void)
{
    CProbe probe;
    for(size_t i=0; i < MSMSProbes.size(); i++){
        probe = MSMSProbes[i];
        if( probe.Selected == false ) continue;
        if( FilterProbe(probe.Pos) ){
            StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
            if( SaveStructure() == false ) return(false);
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorMSMSRandom(void)
{
    std::set<int>::iterator  it = SelectedAtomIds.begin();
    std::set<int>::iterator  ie = SelectedAtomIds.end();

    vout << debug;

    while( it != ie ){
        int atomid = *it;
        // calculate number of probes
        int nprobes = 0;
        for(size_t i=0; i < MSMSProbes.size(); i++){
            CProbe probe = MSMSProbes[i];
            if( (probe.Selected == true) && (probe.AtomId == atomid) ) nprobes++;
        }

        vout << "Selected atom " << atomid + 1 << " has " << nprobes << " vertices." << endl;

        // get probes
        if( nprobes <= Options.GetOptNumOfTrials() ){
            // use all probes
            for(size_t i=0; i < MSMSProbes.size(); i++){
                CProbe probe = MSMSProbes[i];
                if( probe.AtomId == atomid ) {
                    if( FilterProbe(probe.Pos) ){
                        StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
                        if( SaveStructure() == false ) return(false);
                    }
                }
            }
        } else {
            // gen up to Options.GetOptNumOfTrials() random selections
            int ntrials = Options.GetOptNumOfTrials();
            while( ntrials > 0 ){
                int rnd = random() % nprobes;
                for(size_t i=0; i < MSMSProbes.size(); i++){
                    CProbe probe = MSMSProbes[i];
                    if( probe.AtomId == atomid ) {
                        rnd--;
                    }
                    if( rnd == 0 ){
                        if( probe.Selected == false ) break; // already selected
                        MSMSProbes[i].Selected = false;
                        if( FilterProbe(probe.Pos) ){
                            StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
                            if( SaveStructure() == false ) return(false);
                            ntrials--;
                        }
                    }
                }
            }
        }

        it++;
    }

    vout << high;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CVDWGenProbe::FilterProbe(const CPoint& probe)
{
    vector<EFilter>::iterator it = Filters.begin();
    vector<EFilter>::iterator ie = Filters.end();
    while( it != ie ){
        switch(*it){
            case EF_MINMAX:
                if( FilterMinMax(probe) == false ) return(false);
            break;
        case EF_MINONLY:
            if( FilterMinOnly(probe) == false ) return(false);
        break;
        case EF_KEEPALL:
            // nothing to be here
            break;
        }
        it++;
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::FilterMinMax(const CPoint& probe)
{
    bool ok = false;
    double dmin2 = Options.GetOptRMin();
    dmin2 = dmin2*dmin2;
    double dmax2 = Options.GetOptRMax();
    dmax2 = dmax2*dmax2;
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        CPoint apos = Structure.GetPosition(i);
        double d = Square(apos - probe);
        if( d < dmin2 ) return(false);
        if( (d <= dmax2) && SelectedAtoms[i] ) ok = true;
    }

    return(ok);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::FilterMinOnly(const CPoint& probe)
{
    double dmin2 = Options.GetOptRMin();
    dmin2 = dmin2*dmin2;
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        CPoint apos = Structure.GetPosition(i);
        double d = Square(apos - probe);
        if( d < dmin2 ) return(false);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



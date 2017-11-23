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
#include "VDWProbeFilter.hpp"
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <set>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

MAIN_ENTRY(CVDWProbeFilter);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CVDWProbeFilter::CVDWProbeFilter(void)
{
    OutputFile = NULL;
}

//------------------------------------------------------------------------------

CVDWProbeFilter::~CVDWProbeFilter(void)
{
    // close the output file
    if(Options.GetArgStructureName() != "-"){
        if( OutputFile != NULL ) fclose(OutputFile);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CVDWProbeFilter::Init(int argc,char* argv[])
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
    vout << "# vdwprobefilter (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
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

bool CVDWProbeFilter::Run(void)
{
    // load structure
    if( LoadStructure() == false ) return(false);

    // load MSMS vertices
    if( LoadProbes() == false ) return(false);

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

    vout << "# Filter mode                    : " << Options.GetOptFilterType()  << endl;

    // filter vertices
    if( Options.GetOptFilterType() == "minmax" ){
        return(MinMax());
    } else if( Options.GetOptFilterType() == "random" ){
        return(Random());
    } else {
        vout << "<red>>>> ERROR: Unsupported filter type: " << Options.GetOptFilterType() << "</red>" << endl;
    }

    return(false);
}

//------------------------------------------------------------------------------

bool CVDWProbeFilter::MinMax(void)
{
    // min/max algorithm
    vector<int>::iterator sit = AtomIDs.begin();
    vector<int>::iterator sie = AtomIDs.end();

    while( sit != sie ){
        int sat = *sit;
        CPoint apos = Structure.GetPosition(sat);

        std::vector<SProbeData> data;

        // for every probe
        int pind = 0;
        vector<CPoint>::iterator pit = Probes.begin();
        vector<CPoint>::iterator pie = Probes.end();

        while( pit != pie ){
            CPoint      ppos = *pit;
            SProbeData  sdata;
            sdata.ProbeIndex = pind;

            // get D1
            sdata.D1 = Size(ppos - apos);

            // determine D2
            bool first = true;
            for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
                if( j == sat ) continue;
                CPoint apos = Structure.GetPosition(j);
                double d2 = Size(ppos - apos);
                if( first == false ){
                    if( sdata.D2 > d2 ){
                        sdata.D2 = d2;
                    }
                } else {
                    sdata.D2 = d2;
                    first = false;
                }
            }

            data.push_back(sdata);
            pit++;
            pind++;
        }

        // sort data
        sort(data.begin(),data.end(),comp_d1);

        vector<SProbeData>::iterator rit = data.begin();
        double d1 = rit->D1;
        while( rit != data.end() ){
            if( fabs(d1-rit->D1) > 0.001 ) break;
            rit++;
        }

        sort(data.begin(),rit,comp_d2);


        // get best probe
        StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),Probes[data[0].ProbeIndex]);

        // write strcucture
        if( StructureWithProbe.Save(OutputFile) == false ){
            vout << "<red>>>> ERROR: Unable to save a structure to the output file: " << Options.GetArgOutputName() << "</red>" << endl;
            return(false);
        }

        sit++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWProbeFilter::Random(void)
{
    int nstr = Options.GetOptNumOfProbes();
    if( nstr == 0 ){
        nstr = Structure.GetNumberOfAtoms();
    }

    set<int> selected;

    // init seed
    srand (time(NULL));

    for(int i=0; i < nstr; i++){

        // https://stackoverflow.com/questions/10984974/why-do-people-say-there-is-modulo-bias-when-using-a-random-number-generator
        int probe = 0;

        do {
            probe = rand();
        } while (probe >= (RAND_MAX - RAND_MAX % Probes.size()));

        probe %= Probes.size();

        if( selected.find(probe) == selected.end() ){
            StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),Probes[probe]);

            if( StructureWithProbe.Save(OutputFile) == false ){
                vout << "<red>>>> ERROR: Unable to save a structure to the output file: " << Options.GetArgOutputName() << "</red>" << endl;
                return(false);
            }
            selected.insert(probe);
        }

    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CVDWProbeFilter::Finalize(void)
{    
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << "# ==============================================================================" << endl;
    vout << "# vdwprobefilter terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CVDWProbeFilter::LoadStructure(void)
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

    if( Options.GetOptSelectedAtoms() == NULL ){
        // select all atoms
        for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
            AtomIDs.push_back(i);
        }
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
                AtomIDs.push_back(indx);
            }
        }
    }
    vout << "# Number of selected atoms       = " << AtomIDs.size() <<  endl;

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWProbeFilter::LoadProbes(void)
{
        vout << "# MSMS surface vertices (in)     : " << Options.GetArgProbesName() << endl;
    ifstream ifs(Options.GetArgProbesName());
    if( !ifs ){
        vout << "<red>>>> ERROR: Unable to load the file with vertices: " << Options.GetArgProbesName() << "</red>" << endl;
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
        CPoint probe;
        sfs >> probe.x >> probe.y >> probe.z;
        Probes.push_back(probe);
    }

    vout << "# Number of vertices             = " << Probes.size() <<  endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



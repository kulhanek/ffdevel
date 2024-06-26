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
#include <openbabel/atom.h>
#include <openbabel/graphsym.h>
#include <openbabel/obconversion.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "Sphere.hpp"
#include <PeriodicTable.hpp>
#include <sys/time.h>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

MAIN_ENTRY(CVDWGenProbe);

//------------------------------------------------------------------------------

boost::random::mt19937 gen;

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

    if(Options.GetArgOutputName() != "-") {
        vout << "# Output file with probes        = " << Options.GetArgOutputName() << endl;
    } else {
        vout << "# Output file with probes        = - (standard output)" << endl;
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

    // seed initialization
    unsigned int seed = Options.GetOptSeed();
    if( Options.GetOptSeed() == -1 ){
        struct timeval tp;
        gettimeofday(&tp, NULL);
        seed = static_cast<unsigned int>(tp.tv_usec);
    }
        vout << "# Pseudo-random generator seed   = " << seed << endl;
        gen.seed(seed);

    // create output structure
    StructureWithProbe.SetNumberOfAtoms(Structure.GetNumberOfAtoms()+1);
    // copy basal atoms
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        StructureWithProbe.SetPosition(i,Structure.GetPosition(i));
        StructureWithProbe.SetSymbol(i,Structure.GetSymbol(i));
    }
    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),CPoint());
    StructureWithProbe.SetSymbol(Structure.GetNumberOfAtoms(),Options.GetOptProbeSymbol());


    vout << "# Probe filters                  = " << Options.GetOptFilters() << endl;

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
        } else if( filter == "maxprobes" ){
            Filters.push_back(EF_MAXPROBES);
        } else if( filter == "none" ){
            Filters.push_back(EF_NONE);
        } else {
            vout << "<red>>>> ERROR: Unsupported filter: " << filter << "</red>" << endl;
            return(false);
        }
        it++;
    }

    vout << "# Probe generator                = " << Options.GetOptGenerator() << endl;

    // generate structures
    bool result;
    if( Options.GetOptGenerator() == "grid" ){
        result = GeneratorGrid();
    } else if( Options.GetOptGenerator() == "grid-closest" ){
        result = GeneratorGridClosest();
    } else if( Options.GetOptGenerator() == "random" ){
        result = GeneratorRandom();
    } else if( Options.GetOptGenerator() == "msms-all" ){
        result = GeneratorMSMSAll();
    } else if( Options.GetOptGenerator() == "msms-random-per-atom" ){
        result = GeneratorMSMSRandomPerAtom();
    } else if( Options.GetOptGenerator() == "msms-random-per-atom-with-min-prefilter" ){
        result = GeneratorMSMSRandomPerAtomWithMinPreFilter();
    } else if( Options.GetOptGenerator() == "msms-random-per-all" ){
        result = GeneratorMSMSRandomPerAll();
    } else if( Options.GetOptGenerator() == "msms-random-per-all-with-min-prefilter" ){
        result = GeneratorMSMSRandomPerAllWithMinPreFilter();
    } else if( Options.GetOptGenerator() == "sphere-min-repulsion" ){
        result = GeneratorSphereMinRepulsion();
    } else if( Options.GetOptGenerator() == "sphere-max-distance" ){
        result = GeneratorSphereMaxDistance();
    } else if( Options.GetOptGenerator() == "vsepr" ){
        result = GeneratorVSEPR();
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
        vout << "# Input XYZ structure            = " << Options.GetArgStructureName() << endl;
    } else {
        vout << "# Input XYZ structure            = - (standard input)" << endl;
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
    Weighths.resize(Structure.GetNumberOfAtoms());
    for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
        SelectedAtoms[i] = false;
        Weighths[i] = 0.0;
    }

    vout << "# Selected atoms                 = " << Options.GetOptSelectedAtoms() <<  endl;

    // read structure as OBMol
    ifstream sin;
    sin.open(Options.GetArgStructureName());
    if( !sin ) {
        vout << "<red>>>> ERROR: Unable to load (as OBMol) the XYZ structure: " << Options.GetArgStructureName() << "</red>" << endl;
        return(false);
    }

    OBConversion    conv(&sin, &cout);
    OBFormat*       obFormat = conv.FormatFromExt(Options.GetArgStructureName());


    if( obFormat == NULL ) {
        vout << "<red>>>> ERROR: Unsupported input file format!</red>" << endl;
        return(false);
    }

    if( ! conv.SetInFormat(obFormat) ) {
        vout << "<red>>>> ERROR: Unsupported input file format!</red>" << endl;
        return(false);
    }

    if( ! conv.Read(&Mol) ) {
        vout << "<red>>>> ERROR: Unable to load obmol from the stream!</red>" << endl;
        return(false);
    }

    Mol.ConnectTheDots();
    Mol.PerceiveBondOrders();

    int sat = 0;
    // gen list of selected atoms
    if( Options.GetOptSelectedAtoms() == "all" ){
        // select all atoms
        for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
            SelectedAtoms[i] = true;
            SelectedAtomIds.insert(i);
            Weighths[i] = 1.0;
        }
        sat = Structure.GetNumberOfAtoms();
    } else if( (Options.GetOptSelectedAtoms() == "symclasses") || (Options.GetOptSelectedAtoms() == "symclasses-redundant") ){

        // generate symmetry classes
        OBGraphSym  graph_sym(&Mol);

        SymmClasses.clear();
        graph_sym.GetSymmetry(SymmClasses);

        if( SymmClasses.size() != Structure.GetNumberOfAtoms() ){
            vout << "<red>>>> ERROR: Symmetry clases: " << SymmClasses.size() << " != number of atoms: " << Structure.GetNumberOfAtoms() << "!</red>" << endl;
            return(false);
        }

        // make list of unique classes
        std::set<unsigned int>  unique_classes;
        for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
            unique_classes.insert(SymmClasses[i]);
        }

        // for each class find closest atom to previously selected atoms
        std::set<unsigned int>::iterator it = unique_classes.begin();
        std::set<unsigned int>::iterator ie = unique_classes.end();
        CPoint       com_sum;
        double       com_num = 0;
        CPoint       com;
        double       d, min_d;
        CSmallString sats;
        sat = 0;
        while( it != ie ){
            unsigned int ssc = *it;
            int      da = -1;
            for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
                unsigned int asc = SymmClasses[i];
                if( asc != ssc ) continue;
                if( com_num == 0 ){
                    da = i;
                    break;
                }
                if( da == -1 ){
                    da = i;
                    min_d = Size(Structure.GetPosition(da) - com);
                }
                d = Size(Structure.GetPosition(i) - com);
                if( d < min_d ){
                    da = i;
                    min_d = Size(Structure.GetPosition(da) - com);
                }
            }
            if( da == -1 ){
                ES_ERROR("da == -1");
            }

            double  nst = 0;
            for(int i=0; i< Structure.GetNumberOfAtoms(); i++){
                unsigned int asc = SymmClasses[i];
                if( asc != ssc ) continue;
                nst++;
            }

            com_num++;
            com_sum = com_sum + Structure.GetPosition(da);
            com = com_sum * (1.0/com_num);
            SelectedAtomIds.insert(da);
            SelectedAtoms[da] = true;
            Weighths[da] = nst;

            if( sats != NULL ) sats << ",";
            sats << da+1;
            sat++;
            it++;
        }

        // extend selection
        if( Options.GetOptSelectedAtoms() == "symclasses-redundant" ){
            std::set<int>::iterator  it = SelectedAtomIds.begin();
            std::set<int>::iterator  ie = SelectedAtomIds.end();

            while( it != ie ){
                int     atomid = *it;
                int     asc = SymmClasses[atomid];
                OBAtom* p_atom = Mol.GetAtom(atomid+1);
                it++;

                // only terminals
                if( p_atom->GetTotalValence() != 1 ) continue;

                OBBondIterator bi;

                // root atom
                OBAtom* p_na = p_atom->BeginNbrAtom(bi);

                // calc correction to weight
                double ext = 1;
                OBAtom* p_sa = p_na->BeginNbrAtom(bi);
                while( p_sa != NULL ){
                    if( p_sa != p_atom ){
                        int da = p_sa->GetIdx()-1;
                        if( SymmClasses[da] == asc ){
                            ext++;
                        }
                    }
                    p_sa = p_na->NextNbrAtom(bi);
                }

                Weighths[atomid] =  Weighths[atomid] / ext;

                // again
                p_sa = p_na->BeginNbrAtom(bi);
                while( p_sa != NULL ){
                    if( p_sa != p_atom ){
                        int da = p_sa->GetIdx()-1;
                        if( SymmClasses[da] == asc ){
                            SelectedAtomIds.insert(da);
                            SelectedAtoms[da] = true;
                            Weighths[da] = Weighths[atomid];
                        }
                    }
                    p_sa = p_na->NextNbrAtom(bi);
                }
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
                Weighths[i] = 1.0;
            }
        }
    }
    vout << "# Number of selected atoms       = " << sat <<  endl;



    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::LoadMSMSProbes(void)
{
        vout << "# MSMS surface vertices          = " << Options.GetOptMSMSSurfaceVertices() << endl;
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
        CPoint norm;
        int    num;
        probe.AtomId = -1;
        probe.Selected = true;
        probe.Used = false;
        sfs >> probe.Pos.x >> probe.Pos.y >> probe.Pos.z >> norm.x >> norm.y >> norm.z >> num >> probe.AtomId;
        probe.AtomId--;
        MSMSProbes.push_back(probe);
    }

    vout << "# Number of vertices             = " << MSMSProbes.size() <<  endl;

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

    vout << "# Number of selected vertices    = " << sel <<  endl;
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::SaveStructure(double w)
{
    NumOfProbes++;

    CSmallString comment;
    CSmallString sw;
    sw.DoubleToStr(w,"%10.6f");
    comment << "W= " << sw;
    StructureWithProbe.SetComment(comment);

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
                EFilterResult result = FilterProbe(probe);
                if( result == EFR_STOP ) return(false);
                if( result == EFR_OK ){
                    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
                    if( SaveStructure() == false ) return(false);
                }
            }
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorGridClosest(void)
{
    CPoint probe;
    for(probe.x = Min.x; probe.x <= Max.x; probe.x += Options.GetOptSpacing()){
        for(probe.y = Min.y; probe.y <= Max.y; probe.y += Options.GetOptSpacing()){
            for(probe.z = Min.z; probe.z <= Max.z; probe.z += Options.GetOptSpacing()){
                // determine closest atom
                double min;
                int minid = -1;
                for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
                    CPoint apos = Structure.GetPosition(j);
                    double d = Square(apos - probe);
                    if( (j == 0) || (d < min) ){
                        min = d;
                        minid = j;
                    }
                }
                // is atom selected?
                if( SelectedAtomIds.count(minid) != 1 ) continue;

                // filter probe
                EFilterResult result = FilterProbe(probe);
                if( result == EFR_STOP ) return(false);
                if( result == EFR_OK ){
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
    double spacing = Options.GetOptSpacing();
    int nx = (Max.x - Min.x) / spacing;
    int ny = (Max.y - Min.y) / spacing;
    int nz = (Max.z - Min.z) / spacing;
    boost::random::uniform_int_distribution<> dist1(0, nx-1);
    boost::random::uniform_int_distribution<> dist2(0, ny-1);
    boost::random::uniform_int_distribution<> dist3(0, nz-1);

    CPoint probe;
    for(int i=0; i < Options.GetOptNumOfTrials(); i++){
        probe.x = dist1(gen)*spacing + Min.x;
        probe.y = dist2(gen)*spacing + Min.y;
        probe.z = dist3(gen)*spacing + Min.z;
        EFilterResult result = FilterProbe(probe);
        if( result == EFR_STOP ) return(false);
        if( result == EFR_OK ){
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
        EFilterResult result = FilterProbe(probe.Pos);
        if( result == EFR_STOP ) return(false);
        if( result == EFR_OK ){
            StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
            if( SaveStructure() == false ) return(false);
        }
    }
    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorMSMSRandomPerAtomWithMinPreFilter(void)
{
    // run min filter
    for(size_t i=0; i < MSMSProbes.size(); i++){
        CProbe probe = MSMSProbes[i];
        if( FilterMinOnly(probe.Pos) != EFR_OK ) MSMSProbes[i].Selected = false;
    }

    // generate probes
    return( GeneratorMSMSRandomPerAtom() );
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorMSMSRandomPerAtom(void)
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
        it++;

        vout << "Selected atom " << atomid + 1 << " has " << nprobes << " vertices." << endl;

        if( nprobes <= Options.GetOptMaxProbes() ) continue;

        // get probes
        if( nprobes <= Options.GetOptNumOfTrials() ){
            // use all probes
            for(size_t i=0; i < MSMSProbes.size(); i++){
                CProbe probe = MSMSProbes[i];
                if( (probe.Selected == true) && (probe.AtomId == atomid) ) {
                    EFilterResult result = FilterProbe(probe.Pos);
                    if( result == EFR_STOP ) return(false);
                    if( result == EFR_OK ){
                        StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
                        if( SaveStructure() == false ) return(false);
                    }
                }
            }
        } else {
            // gen up to Options.GetOptNumOfTrials() random selections
            int ntrials = Options.GetOptNumOfTrials();
            while( ntrials > 0 ){
                boost::random::uniform_int_distribution<> dist(0, nprobes-1);
                int rnd = dist(gen);
                for(size_t i=0; i < MSMSProbes.size(); i++){
                    CProbe probe = MSMSProbes[i];
                    if( (probe.Selected == true) && (probe.AtomId == atomid) ) {
                        rnd--;
                        if( rnd == 0 ){
                            if( probe.Used == true ) break;
                            MSMSProbes[i].Used = true;
                            EFilterResult result = FilterProbe(probe.Pos);
                            if( result == EFR_STOP ) return(false);
                            if( result == EFR_OK ){
                                StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
                                if( SaveStructure() == false ) return(false);
                                ntrials--;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    vout << high;

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorMSMSRandomPerAllWithMinPreFilter(void)
{
    // run min filter
    for(size_t i=0; i < MSMSProbes.size(); i++){
        CProbe probe = MSMSProbes[i];
        if( FilterMinOnly(probe.Pos) != EFR_OK ) MSMSProbes[i].Selected = false;
    }

    // generate probes
    return( GeneratorMSMSRandomPerAll() );
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorMSMSRandomPerAll(void)
{

    // calculate number of probes
    int nprobes = 0;
    for(size_t i=0; i < MSMSProbes.size(); i++){
        CProbe probe = MSMSProbes[i];
        if( probe.Selected == true ) nprobes++;
    }

    if( nprobes <= Options.GetOptMaxProbes() ) return(true);

    // get probes
    if( nprobes <= Options.GetOptNumOfTrials() ){
        // use all probes
        for(size_t i=0; i < MSMSProbes.size(); i++){
            CProbe probe = MSMSProbes[i];
            if( probe.Selected == true ) {
                EFilterResult result = FilterProbe(probe.Pos);
                if( result == EFR_STOP ) return(false);
                if( result == EFR_OK ){
                    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
                    if( SaveStructure() == false ) return(false);
                }
            }
        }
    } else {
        // gen up to Options.GetOptNumOfTrials() random selections
        int ntrials = Options.GetOptNumOfTrials();
        while( ntrials > 0 ){
            boost::random::uniform_int_distribution<> dist(0, nprobes-1);
            int rnd = dist(gen);
            for(size_t i=0; i < MSMSProbes.size(); i++){
                CProbe probe = MSMSProbes[i];
                if( probe.Selected == true ) {
                    rnd--;
                    if( rnd == 0 ){
                        if( probe.Used == true ) break;
                        MSMSProbes[i].Used = true;
                        EFilterResult result = FilterProbe(probe.Pos);
                        if( result == EFR_STOP ) return(false);
                        if( result == EFR_OK ){
                            StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe.Pos);
                            if( SaveStructure() == false ) return(false);
                            ntrials--;
                        }
                        break;
                    }
                }
            }
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorSphereMinRepulsion(void)
{
    std::set<int>::iterator  it = SelectedAtomIds.begin();
    std::set<int>::iterator  ie = SelectedAtomIds.end();

    CSphere sphere;
    if( sphere.SetTessellationQuality(Options.GetOptSphereQuality()) == false ){
        ES_ERROR("unable to initialize sphere");
    }

    vector<int> z;
    for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
        z.push_back(PeriodicTable.SearchZBySymbol(Structure.GetSymbol(j)));
    }

    double rprobe = PeriodicTable.GetVdWRadius(PeriodicTable.SearchZBySymbol(Options.GetOptProbeSymbol()));

    while( it != ie ){
        int    vertex = -1;
        CPoint origin;

        if( Options.GetOptAllWithinSymmClass() == true ) {
            double min = std::numeric_limits<double>::max();
            for(int atomid=0; atomid <= Structure.GetNumberOfAtoms(); atomid++ ) {
                if( SymmClasses[atomid] != SymmClasses[*it] ) continue; // consider atoms from the same symmetry class
                CPoint lorigin = Structure.GetPosition(atomid);
                for(unsigned int i=0; i < sphere.GetNumberOfVertices(); i++){
                    double interaction = 0.0;
                    for(double r = Options.GetOptRMin(); r <= Options.GetOptRMax(); r += Options.GetOptSpacing()){
                        CPoint v1 = sphere.GetVertice(i,lorigin,r);
                        for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
                            double rvdw = PeriodicTable.GetVdWRadius(z[j]);
                            CPoint apos = Structure.GetPosition(j);
                            double d = Size(apos - v1);
                            double r0 = rprobe + rvdw;
                            interaction = interaction + pow(r0/d,12);
                        }
                    }
                    if( min > interaction ){
                        vertex = i;
                        min = interaction;
                        origin = lorigin;
                    }
                }
            }
        } else {
            int atomid = *it;
            origin = Structure.GetPosition(atomid);
            double min = std::numeric_limits<double>::max();
            for(unsigned int i=0; i < sphere.GetNumberOfVertices(); i++){
                double interaction = 0.0;
                for(double r = Options.GetOptRMin(); r <= Options.GetOptRMax(); r += Options.GetOptSpacing()){
                    CPoint v1 = sphere.GetVertice(i,origin,r);
                    for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
                        double rvdw = PeriodicTable.GetVdWRadius(z[j]);
                        CPoint apos = Structure.GetPosition(j);
                        double d = Size(apos - v1);
                        double r0 = rprobe + rvdw;
                        interaction = interaction + pow(r0/d,12);
                    }
                }
                if( min > interaction ){
                    vertex = i;
                    min = interaction;
                }
            }
        }

        if( vertex >= 0 ) {
            for(double r = Options.GetOptRMin(); r <= Options.GetOptRMax(); r += Options.GetOptSpacing()){
                CPoint probe = sphere.GetVertice(vertex,origin,r);
                EFilterResult result = FilterProbe(probe);
                if( result == EFR_STOP ) return(false);
                if( result == EFR_OK ){
                    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
                    if( SaveStructure() == false ) return(false);
                }
            }
        }
        it++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorSphereMaxDistance(void)
{
    std::set<int>::iterator  it = SelectedAtomIds.begin();
    std::set<int>::iterator  ie = SelectedAtomIds.end();

    CSphere sphere;
    if( sphere.SetTessellationQuality(Options.GetOptSphereQuality()) == false ){
        ES_ERROR("unable to initialize sphere");
    }

    vector<int> z;
    for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
        z.push_back(PeriodicTable.SearchZBySymbol(Structure.GetSymbol(j)));
    }

    double rprobe = PeriodicTable.GetVdWRadius(PeriodicTable.SearchZBySymbol(Options.GetOptProbeSymbol()));

    while( it != ie ){
        int atomid = *it;
        CPoint origin = Structure.GetPosition(atomid);

        double dmax = 0.0;
        int    vertex;
        for(unsigned int i=0; i < sphere.GetNumberOfVertices(); i++){
            double dmintot = 0.0;
            for(double r = Options.GetOptRMin(); r <= Options.GetOptRMax(); r += Options.GetOptSpacing()){
                CPoint v1 = sphere.GetVertice(i,origin,r);
                double dmin = std::numeric_limits<double>::max();
                for(int j=0; j < Structure.GetNumberOfAtoms(); j++){
                    double rvdw = PeriodicTable.GetVdWRadius(z[j]);
                    CPoint apos = Structure.GetPosition(j);
                    double d = Size(v1 - apos);
                    if( d < dmin ){
                        dmin = d;
                    }
                }
                dmintot = dmintot + dmin;
            }
            if( dmintot > dmax ){
                vertex = i;
                dmax = dmintot;
            }
        }
        for(double r = Options.GetOptRMin(); r <= Options.GetOptRMax(); r += Options.GetOptSpacing()){
            CPoint probe = sphere.GetVertice(vertex,origin,r);
            EFilterResult result = FilterProbe(probe);
            if( result == EFR_STOP ) return(false);
            if( result == EFR_OK ){
                StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
                if( SaveStructure() == false ) return(false);
            }
        }
    it++;
    }

    return(true);

}

//------------------------------------------------------------------------------

bool CVDWGenProbe::GeneratorVSEPR(void)
{
    std::set<int>::iterator  it = SelectedAtomIds.begin();
    std::set<int>::iterator  ie = SelectedAtomIds.end();

    while( it != ie ){
        int     atomid = *it;
        it++;
        OBAtom* p_atom = Mol.GetAtom(atomid+1);
        int     valence = p_atom->GetTotalValence();
        CPoint  origin = Structure.GetPosition(atomid);
        CPoint  dir;
        bool    bothsides = false;
        switch(valence){
            case 1:{
                OBBondIterator bi;
                OBAtom* p_na = p_atom->BeginNbrAtom(bi);
                if( p_na == NULL ) continue;
                dir = origin - Structure.GetPosition(p_na->GetIdx()-1);
                dir.Normalize();
            }
            break;
            case 2:{
                OBBondIterator bi;
                OBAtom* p_na = p_atom->BeginNbrAtom(bi);
                CPoint first = Structure.GetPosition(p_na->GetIdx()-1) - origin;
                first.Normalize();
                p_na = p_atom->NextNbrAtom(bi);
                CPoint second = Structure.GetPosition(p_na->GetIdx()-1) - origin;
                second.Normalize();

                if( Angle(first,second) > 160.0*M_PI/180.0 ){
                    // linear
                    ES_ERROR("linear AB2 - not implemented");
                } else {
                    // like H2O
                    dir = - first - second;
                    dir.Normalize();
                    bothsides = true;
                }
            }
            break;
            case 3:{
                CPoint sum;
                OBBondIterator bi;
                OBAtom* p_na = p_atom->BeginNbrAtom(bi);
                while( p_na != NULL ){
                    CPoint diff = Structure.GetPosition(p_na->GetIdx()-1) - origin;
                    diff.Normalize();
                    sum += diff;
                    p_na = p_atom->NextNbrAtom(bi);
                }
                if( Size(sum) > 0.5 ){
                    // like NH3
                    dir = -sum;
                    dir.Normalize();
                } else {
                    // planar
                    OBAtom* p_na = p_atom->BeginNbrAtom(bi);
                    CPoint first = Structure.GetPosition(p_na->GetIdx()-1) - origin;
                    first.Normalize();
                    p_na = p_atom->NextNbrAtom(bi);
                    CPoint second = Structure.GetPosition(p_na->GetIdx()-1) - origin;
                    second.Normalize();
                    p_na = p_atom->NextNbrAtom(bi);
                    CPoint third = Structure.GetPosition(p_na->GetIdx()-1) - origin;
                    third.Normalize();

                    CPoint v1 = CrossDot(first,second);
                    CPoint v2 = CrossDot(first,third);
                    CPoint v3 = CrossDot(second,third);
                    v2 = v2*VectDot(v1,v2);
                    v3 = v3*VectDot(v1,v3);

                    dir = v1 + v2 + v3;
                    dir.Normalize();
                    bothsides = true;
                }
            }
            break;
            default:
                ES_WARNING("atom ignored");
                continue;
        }

        for(double r = Options.GetOptRMin(); r <= Options.GetOptRMax(); r += Options.GetOptSpacing()){
            CPoint probe;
            probe = origin + dir*r;
            double w = Weighths[atomid];
            if( bothsides ) w /= 2.0;

            EFilterResult result = FilterProbe(probe);
            if( result == EFR_STOP ) return(false);

            if( result == EFR_OK ){
                StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
                if( SaveStructure(w) == false ) return(false);
            }
            if( bothsides ){
                probe = origin - dir*r;
                EFilterResult result = FilterProbe(probe);
                if( result == EFR_STOP ) return(false);
                if( result == EFR_OK ){
                    StructureWithProbe.SetPosition(Structure.GetNumberOfAtoms(),probe);
                    if( SaveStructure(w) == false ) return(false);
                }
            }
        }
    }

    return(true);

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

EFilterResult CVDWGenProbe::FilterProbe(const CPoint& probe)
{
    vector<EFilter>::iterator it = Filters.begin();
    vector<EFilter>::iterator ie = Filters.end();
    while( it != ie ){
        switch(*it){
            case EF_MINMAX:{
                EFilterResult result = FilterMinMax(probe);
                if( result != EFR_OK ) return(result);
                }
                break;
            case EF_MINONLY:{
                EFilterResult result = FilterMinOnly(probe);
                if( result != EFR_OK ) return(result);
                }
                break;
            case EF_MAXPROBES:{
                EFilterResult result = FilterMaxProbes(probe);
                if( result != EFR_OK ) return(result);
                }
                break;
            case EF_NONE:
                // nothing to be here
                break;
        }
        it++;
    }
    return(EFR_OK);
}

//------------------------------------------------------------------------------

EFilterResult CVDWGenProbe::FilterMinMax(const CPoint& probe)
{
    bool ok = false;
    double dmin2 = Options.GetOptRMin();
    dmin2 = dmin2*dmin2;
    double dmax2 = Options.GetOptRMax();
    dmax2 = dmax2*dmax2;
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        CPoint apos = Structure.GetPosition(i);
        double d = Square(apos - probe);
        if( d < dmin2 ) return(EFR_NEXT);
        if( (d <= dmax2) && SelectedAtoms[i] ) ok = true;
    }

    return(EFR_OK);
}

//------------------------------------------------------------------------------

EFilterResult CVDWGenProbe::FilterMinOnly(const CPoint& probe)
{
    double dmin2 = Options.GetOptRMin();
    dmin2 = dmin2*dmin2;
    for(int i=0; i < Structure.GetNumberOfAtoms(); i++){
        CPoint apos = Structure.GetPosition(i);
        double d = Square(apos - probe);
        if( d < dmin2 ) return(EFR_NEXT);
    }

    return(EFR_OK);
}

//------------------------------------------------------------------------------

EFilterResult CVDWGenProbe::FilterMaxProbes(const CPoint& probe)
{
    if( NumOfProbes == Options.GetOptMaxProbes() ){
        vout << "<red>>>> ERROR: MaxProbes limit was reached!</red>" << endl;
        return(EFR_STOP);
    }
    return(EFR_OK);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



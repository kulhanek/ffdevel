#ifndef VDWGenProbeH
#define VDWGenProbeH
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

#include <stdio.h>
#include "VDWGenProbeOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <vector>
#include <set>
#include <Point.hpp>
#include <XYZStructure.hpp>
#include <openbabel/mol.h>

using namespace OpenBabel;

//------------------------------------------------------------------------------

class SProbeData {
public:
    int     ProbeIndex;
    double  D1;
    double  D2;
};

//------------------------------------------------------------------------------

enum EFilter{
    EF_MINMAX = 0,
    EF_MINONLY =1,
    EF_MAXPROBES = 2,
    EF_NONE = 3
};

//------------------------------------------------------------------------------

enum EFilterResult{
    EFR_OK = 0,
    EFR_STOP =1,
    EFR_NEXT = 2
};

//------------------------------------------------------------------------------

class CProbe {
public:
    CPoint  Pos;
    int     AtomId;
    bool    Selected;
    bool    Used;
};

//------------------------------------------------------------------------------

class CVDWGenProbe {
public:
    // constructor
    CVDWGenProbe(void);
    ~CVDWGenProbe(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CVDWGenProbeOptions     Options;            // program options
    CXYZStructure           Structure;          // input structure
    std::vector<CProbe>     MSMSProbes;         // probes from MSMS surface
    std::vector<bool>       SelectedAtoms;      // consider only these atoms
    std::vector<double>     Weighths;           // point weights
    std::set<int>           SelectedAtomIds;
    FILE*                   OutputFile;         // output file
    CXYZStructure           StructureWithProbe; // output structure
    CPoint                  Min,Max;            // min/max conners
    std::vector<EFilter>    Filters;
    int                     NumOfProbes;
    OBMol                   Mol;


    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    /// load/save structure
    bool LoadStructure(void);
    bool SaveStructure(double w=1.0);
    bool LoadMSMSProbes(void);

    /// generators
    bool GeneratorGrid(void);
    bool GeneratorGridClosest(void);
    bool GeneratorRandom(void);
    bool GeneratorMSMSAll(void);
    bool GeneratorMSMSRandomPerAtom(void);
    bool GeneratorMSMSRandomPerAtomWithMinPreFilter(void);
    bool GeneratorMSMSRandomPerAll(void);
    bool GeneratorMSMSRandomPerAllWithMinPreFilter(void);
    bool GeneratorSphereMinRepulsion(void);
    bool GeneratorVSEPR(void);

    /// filters
    EFilterResult FilterProbe(const CPoint& probe);
    EFilterResult FilterMinMax(const CPoint& probe);
    EFilterResult FilterMinOnly(const CPoint& probe);
    EFilterResult FilterMaxProbes(const CPoint& probe);
};

//------------------------------------------------------------------------------

#endif

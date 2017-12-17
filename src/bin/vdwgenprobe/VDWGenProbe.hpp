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
#include <Point.hpp>
#include <XYZStructure.hpp>

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
    EF_MINONLY =1
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
    CVDWGenProbeOptions     Options;    // program options
    CXYZStructure           Structure;  // input structure
    std::vector<bool>       SelectedAtoms;    // consider only these atoms
    FILE*                   OutputFile; // output file
    CXYZStructure           StructureWithProbe;  // output structure
    CPoint                  Min,Max;    // min/max conners
    std::vector<EFilter>    Filters;
    int                     NumOfProbes;

    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    /// load/save structure
    bool LoadStructure(void);
    bool SaveStructure(void);

    /// generators
    bool GeneratorGrid(void);
    bool GeneratorRandom(void);

    /// filters
    bool FilterProbe(const CPoint& probe);
    bool FilterMinMax(const CPoint& probe);
    bool FilterMinOnly(const CPoint& probe);
};

//------------------------------------------------------------------------------

#endif

#ifndef VDWProbeFilterH
#define VDWProbeFilterH
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
#include "VDWProbeFilterOptions.hpp"
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

class CVDWProbeFilter {
public:
    // constructor
    CVDWProbeFilter(void);
    ~CVDWProbeFilter(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CVDWProbeFilterOptions  Options;    // program options
    CXYZStructure           Structure;  // input structure
    std::vector<CPoint>     Probes;     // probes from MSMS surface
    std::vector<int>        AtomIDs;    // consider only these atoms
    FILE*                   OutputFile; // output file
    CXYZStructure           StructureWithProbe;  // output structure

    // output ------------------------------------
    CTerminalStr            Console;
    CVerboseStr             vout;

    /// load structure
    bool LoadStructure(void);

    /// load probes
    bool LoadProbes(void);

    /// filters
    bool MinMax(void);
    bool Random(void);
};

//------------------------------------------------------------------------------

#endif

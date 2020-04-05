#ifndef TopFilterH
#define TopFilterH
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

#include <stdio.h>
#include <AmberTopology.hpp>
#include <AmberRestart.hpp>
#include "TopFilterOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>

//------------------------------------------------------------------------------

class CTopFilter {
public:
    // constructor
    CTopFilter(void);
    ~CTopFilter(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CTopFilterOptions   Options;            // program options
    CAmberTopology      Topology;           // input topology
    CAmberRestart       Restart;            // input coordinates

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// topology manipulation
    bool LoadTopologyAndCoordinates(void);
    bool SaveTopology(void);
    bool FilterTopology(void);

    /// filters
    bool FilterLAngles(void);
    bool FilterZAngles(void);
    bool FilterLDihedrals(void);
    bool FilterZDihedrals(void);
};

//------------------------------------------------------------------------------

#endif

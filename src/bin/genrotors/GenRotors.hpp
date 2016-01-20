#ifndef GenRotorsH
#define GenRotorsH
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
#include "GenRotorsOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include <ostream>
#include "openbabel/mol.h"

//------------------------------------------------------------------------------

class CGenRotors {
public:
    // constructor
    CGenRotors(void);
    ~CGenRotors(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CGenRotorsOptions   Options;            // program options
    CAmberTopology      Topology;           // input topology
    CAmberRestart       Restart;            // input restart
    OpenBabel::OBMol    OBMolecule;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    /// load topology and restart
    bool LoadTopologyAndRestart(void);

    /// generate babel molecule
    void GenBabelMolecule(void);

    /// write rotatable bonds
    void WriteRotatableBonds(std::ostream& sout);
};

//------------------------------------------------------------------------------

#endif

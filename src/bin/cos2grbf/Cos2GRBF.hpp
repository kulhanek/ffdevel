#ifndef Cos2GRBFH
#define Cos2GRBFH
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
#include "Cos2GRBFOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include "../top2stop/MMTypes.hpp"

//------------------------------------------------------------------------------

class CCos2GRBF {
public:
    // constructor
    CCos2GRBF(void);
    ~CCos2GRBF(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CCos2GRBFOptions    Options;            // program options
    CAmberTopology      Topology;           // input topology

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;
    
    CDihedralType       DihedralType;
    
    bool ReadFrcMod(void);
    bool WriteOutput(void);

};

//------------------------------------------------------------------------------

#endif

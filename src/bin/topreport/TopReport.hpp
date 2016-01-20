#ifndef TopReportH
#define TopReportH
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
#include "TopReportOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include "MMTypes.hpp"
#include <map>

//------------------------------------------------------------------------------

class CTopReport {
public:
    // constructor
    CTopReport(void);
    ~CTopReport(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CTopReportOptions   Options;            // program options
    CAmberTopology      Topology;           // input topology

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    // types -------------------------------------
    std::map<int,CAtomType>     AtomTypes;
    std::map<int,CDihedralType> DihedralTypes;

    /// load topology
    bool LoadTopology(void);

    /// write all data
    void WriteAll(std::ostream& sout);

    /// output methods
    void WriteAtoms(std::ostream& sout);
    void WriteBonds(std::ostream& sout);
    void WriteAngles(std::ostream& sout);
    void WriteDihedrals(std::ostream& sout);
    void WriteNBList(std::ostream& sout);

    /// helper methods
    int           FindAtomTypeIdx(int atidx);
    CDihedralType FindDihedralByTypes(int ia,int ib,int ic,int id);
};

//------------------------------------------------------------------------------

#endif

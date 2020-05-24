#ifndef Top2STopH
#define Top2STopH
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
#include "Top2STopOptions.hpp"
#include <VerboseStr.hpp>
#include <TerminalStr.hpp>
#include "MMTypes.hpp"
#include <map>
#include <list>

#include <openbabel/mol.h>

//------------------------------------------------------------------------------

class CTop2STop {
public:
    // constructor
    CTop2STop(void);
    ~CTop2STop(void);

// main methods ----------------------------------------------------------------
    /// init options
    int Init(int argc,char* argv[]);

    /// main part of program
    bool Run(void);

    /// finalize program
    void Finalize(void);

// section of private data -----------------------------------------------------
private:
    CTop2STopOptions    Options;            // program options
    CAmberTopology      Topology;           // input topology
    CAmberRestart       Coords;

    // output ------------------------------------
    CTerminalStr        Console;
    CVerboseStr         vout;

    // types -------------------------------------
    std::map<int,CAtomType>             AtomTypes;
    std::map<int,CBondType>             BondTypes;
    std::map<int,CAngleType>            AngleTypes;
    std::map<int,CDihedralType>         DihedralTypes;
    std::map<int,CDihedralType>         ImproperTypes;
    std::map<int,int>                   NBTypes;
    std::list<CDihedral>                UniqueDihedrals;
    std::list<CDihedral>                UniqueImpropers;
    std::vector<CDihedralTypeFilter>    DihFilters;
    std::vector<unsigned int>           SymmClasses;

    int dih_mode;       // dihedral mode: 1 - cos; 2 - grbf
    int dih_samp_freq;  // dihedral sampling frequency

    // dimmensions
    int natoms;
    int natom_types;
    int nbonds;
    int nbond_types;
    int nangles;
    int nangle_types;
    int ndihedrals;
    int ndihedral_types;
    int ndihedral_seq_size;
    int nimpropers;
    int nimproper_types;
    int nb_size;
    int nb_sizeij;
    int nb_size14;
    int nnb_types;
    int nsymm_classes;

    /// load topology
    bool LoadTopology(void);
    bool LoadCoords(void);
    bool LoadDihFilters(void);

    /// write all data
    void WriteAll(std::ostream& sout);

    /// output methods
    void WriteAtomTypes(std::ostream& sout);
    void WriteAtoms(std::ostream& sout);

    void WriteBondTypes(std::ostream& sout);
    void WriteBonds(std::ostream& sout);

    void WriteAngleTypes(std::ostream& sout);
    void WriteAngles(std::ostream& sout);

    void WriteDihedralTypes(std::ostream& sout);
    void WriteDihedralSeqCosMode(std::ostream& sout);
    void TransformCosToGRBF(void);
    void SolveTransformation(int type);
    void WriteDihedralSeqGRBFMode(std::ostream& sout);
    void WriteDihedrals(std::ostream& sout);

    void WriteImproperTypes(std::ostream& sout);
    void WriteImpropers(std::ostream& sout);

    void WriteNBTypes(std::ostream& sout);
    void WriteNBListRebuild(std::ostream& sout);
    void WriteNBListKeep(std::ostream& sout);

    void WriteDimensions(std::ostream& sout);

    void TransformTypes(void);

    /// helper methods
    int           FindAtomTypeIdx(int atidx);
    CDihedralType FindDihedralByTypes(int ia,int ib,int ic,int id);
    CDihedralType FindImproperByTypes(int ia,int ib,int ic,int id);
};

//------------------------------------------------------------------------------

#endif

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

#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <SimpleList.hpp>
#include "GenRotors.hpp"
#include <openbabel/rotor.h>
#include <iomanip>
#include <set>

//------------------------------------------------------------------------------

using namespace std;

MAIN_ENTRY(CGenRotors);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGenRotors::CGenRotors(void)
{

}

//------------------------------------------------------------------------------

CGenRotors::~CGenRotors(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CGenRotors::Init(int argc,char* argv[])
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
    vout << "# genrotors (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgTopologyName() != "-") {
        vout << "# AMBER topology file (in)       : " << Options.GetArgTopologyName() << endl;
    } else {
        vout << "# AMBER topology file (in)       : - (standard input)" << endl;
    }
    if(Options.GetArgCoordinateName() != "-") {
        vout << "# AMBER coordinates file (in)    : " << Options.GetArgCoordinateName() << endl;
    } else {
        vout << "# AMBER coordinates file (in)    : - (standard input)" << endl;
    }
    if(Options.GetArgRotorName() != "-") {
        vout << "# Freely rotatable bonds (out)   : " << Options.GetArgRotorName() << endl;
    } else {
        vout << "# Freely rotatable bonds (out)   : - (standard output)" << endl;
    }
        vout << "# Include terminal rotors        : " << bool_to_str(Options.GetOptIncludeTerminals()) << endl;

    return( SO_CONTINUE );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CGenRotors::Run(void)
{
    // load topology and coordinates
    if( LoadTopologyAndRestart() == false ) return(false);

    // generate babel structure
    GenBabelMolecule();

    // write results
    if( Options.GetArgRotorName() == "-" ){
        WriteRotatableBonds(cout);
    } else {
        ofstream ofs;
        ofs.open(Options.GetArgRotorName());
        if( ! ofs ){
            vout << "<red>>>> ERROR: Unable to open file with rotatable bonds: " << Options.GetArgRotorName() << "</red>" << endl;
            return(false);
        }
        WriteRotatableBonds(ofs);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGenRotors::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# genrotors terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CGenRotors::LoadTopologyAndRestart(void)
{
    if( Topology.Load(Options.GetArgTopologyName(),true) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER topology: " << Options.GetArgTopologyName() << "</red>" << endl;
        return(false);
    }

    Restart.AssignTopology(&Topology);

    if( Restart.Load(Options.GetArgCoordinateName(),true) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER coordinates: " << Options.GetArgCoordinateName() << "</red>" << endl;
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CGenRotors::GenBabelMolecule(void)
{
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        // create new atom
        OpenBabel::OBAtom* p_obatom = OBMolecule.NewAtom();
        if( p_atom->GetAtomicNumber() > 0 ){
            p_obatom->SetAtomicNum(p_atom->GetAtomicNumber());
        } else {
            p_obatom->SetAtomicNum(p_atom->GuessZ());
        }
        CPoint pos = Restart.GetPosition(i);
        p_obatom->SetVector(pos.x,pos.y,pos.z);
    }

    OBMolecule.ConnectTheDots();
    OBMolecule.PerceiveBondOrders();
}

//------------------------------------------------------------------------------

void CGenRotors::WriteRotatableBonds(std::ostream& sout)
{
    sout << "[rotors]" << endl;

    OpenBabel::OBRotorList rot_list;
    rot_list.FindRotors(OBMolecule);

    OpenBabel::OBRotorIterator it = rot_list.BeginRotors();
    OpenBabel::OBRotorIterator ie = rot_list.EndRotors();

    std::set<OpenBabel::OBBond*>    rot_bonds;

    while( it != ie ){
        OpenBabel::OBRotor* p_rot = *it;
        OpenBabel::OBBond*  p_bond = p_rot->GetBond();
        if( p_bond ){
            rot_bonds.insert(p_bond);
            sout << right << setw(6) << p_bond->GetBeginAtomIdx() << " ";
            sout << right << setw(6) << p_bond->GetEndAtomIdx() << endl;
        }
        it++;
    }

    if( ! Options.GetOptIncludeTerminals() ) return;

    OpenBabel::OBBondIterator ibt = OBMolecule.BeginBonds();
    OpenBabel::OBBondIterator ibe = OBMolecule.EndBonds();

    while( ibt != ibe ){
        OpenBabel::OBBond*  p_bond = *ibt;
        ibt++;
        if( rot_bonds.count(p_bond) > 0 ) continue; // already detected
        if( p_bond->IsInRing() ) continue;  // part of ring
        if( ! p_bond->IsSingle() ) continue;    // not single bond

        bool terminal1 = false;
        bool terminal2 = false;

        // is terminal bond?
        OpenBabel::OBAtom* p_fatom = p_bond->GetBeginAtom();
        OpenBabel::OBAtom* p_satom = p_bond->GetEndAtom();

        OpenBabel::OBBondIterator   atit;
        OpenBabel::OBAtom* p_nbatom;
        p_nbatom = p_fatom->BeginNbrAtom(atit);
        while( p_nbatom != NULL ){
            if( p_nbatom == p_satom ){
                p_nbatom = p_fatom->NextNbrAtom(atit);
                continue;
            }
            terminal1 = true;
            if( p_nbatom->CountBondsOfOrder(1) != 1 ){
                terminal1 = false;
                break;
            }
            p_nbatom = p_fatom->NextNbrAtom(atit);
        }

        p_nbatom = p_satom->BeginNbrAtom(atit);
        while( p_nbatom != NULL ){
            if( p_nbatom == p_fatom ){
                p_nbatom = p_satom->NextNbrAtom(atit);
                continue;
            }
            terminal2 = true;
            if( p_nbatom->CountBondsOfOrder(1) != 1 ){
                terminal2 = false;
                break;
            }
            p_nbatom = p_satom->NextNbrAtom(atit);
        }
        if( (! terminal1) && (! terminal2) )  continue;
        rot_bonds.insert(p_bond);

        int count1 = 0;
        p_nbatom = p_fatom->BeginNbrAtom(atit);
        while( p_nbatom != NULL ){
            count1++;
            p_nbatom = p_fatom->NextNbrAtom(atit);
        }

        int count2 = 0;
        p_nbatom = p_satom->BeginNbrAtom(atit);
        while( p_nbatom != NULL ){
            count2++;
            p_nbatom = p_satom->NextNbrAtom(atit);
        }

        if( count1 > count2 ){
            sout << right << setw(6) << p_bond->GetBeginAtomIdx() << " ";
            sout << right << setw(6) << p_bond->GetEndAtomIdx() << endl;
        } else {
            sout << right << setw(6) << p_bond->GetEndAtomIdx() << " ";
            sout << right << setw(6) << p_bond->GetBeginAtomIdx() << endl;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



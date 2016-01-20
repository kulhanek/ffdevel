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
#include <ctype.h>
#include "TopReport.hpp"
#include <map>
#include <iomanip>
#include <list>
#include <locale>

//------------------------------------------------------------------------------

using namespace std;

MAIN_ENTRY(CTopReport);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CTopReport::CTopReport(void)
{
}

//------------------------------------------------------------------------------

CTopReport::~CTopReport(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CTopReport::Init(int argc,char* argv[])
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
    vout << "# topreport (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgTopologyName() != "-") {
        vout << "# AMBER topology file (in) : " << Options.GetArgTopologyName() << endl;
    } else {
        vout << "# AMBER topology file (in) : - (standard input)" << endl;
    }
    if(Options.GetArgTopReportName() != "-") {
        vout << "# Topology report (out)    : " << Options.GetArgTopReportName() << endl;
    } else {
        vout << "# Topology report (out)    : - (standard output)" << endl;
    }

    return( SO_CONTINUE );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CTopReport::Run(void)
{
    // load topology
    if( LoadTopology() == false ) return(false);

    // save topology
    if( Options.GetArgTopReportName() == "-" ){
        WriteAll(cout);
    } else {
        ofstream ifs;
        ifs.open(Options.GetArgTopReportName());
        if( ! ifs ){
            vout << "<red>>>> ERROR: Unable to open simplified topology: " << Options.GetArgTopReportName() << "</red>" << endl;
            return(false);
        }
        WriteAll(ifs);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CTopReport::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# topreport terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CTopReport::LoadTopology(void)
{
    if( Topology.Load(Options.GetArgTopologyName(),true) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER topology: " << Options.GetArgTopologyName() << "</red>" << endl;
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CTopReport::WriteAll(ostream& sout)
{
    WriteAtoms(sout);
    WriteBonds(sout);
    WriteAngles(sout);
    WriteDihedrals(sout);
    WriteNBList(sout);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CTopReport::WriteAtoms(ostream& sout)
{
    std::map<std::string,CAtomType> unique_types;

    // extract atom types
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        CAtomType mmtype;
        mmtype.name = p_atom->GetType();
        mmtype.mass = p_atom->GetMass();
        if( p_atom->GetAtomicNumber() > 0 ){
            mmtype.z = p_atom->GetAtomicNumber();
        } else {
            mmtype.z = p_atom->GuessZ();
        }
        if( unique_types[mmtype.name].idx == -1 ){

            if( Topology.NonBondedList.GetNonBondedType(p_atom,p_atom) < 0 ) {
                RUNTIME_ERROR("10-12 interaction is not supported");
            }

            double a,b,eps,rii;

            a = Topology.NonBondedList.GetAParam(p_atom,p_atom);
            b = Topology.NonBondedList.GetBParam(p_atom,p_atom);

            eps = 0.0;
            rii = 0.0;
            if( a != 0.0 ) {
                eps = b*b / (4.0 * a);
                rii = pow(2*a/b,1.0/6.0) * 0.5;
            }
            mmtype.eps = eps;
            mmtype.r0 = rii;

            unique_types[mmtype.name] = mmtype;
        } else {
            if( unique_types[mmtype.name] != mmtype ){
                RUNTIME_ERROR("atom type mismatch");
            }
        }
    }

    // list types and set their indexes
    std::map<std::string,CAtomType>::iterator uit = unique_types.begin();
    std::map<std::string,CAtomType>::iterator uie = unique_types.end();

    int idx = 1;
    while( uit != uie ){
        CAtomType mmtype = uit->second;
        mmtype.idx = idx;
        AtomTypes[mmtype.idx] = mmtype;
        uit++;
        idx++;
    }

    sout << "[atoms]" << endl;
    sout << "! Index Name Type ResID ResN     Charge       Mass  Z           eps            r*" << endl;

    // write atoms
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        sout << right << setw(7) << i+1 << " ";
        sout << left << setw(4) << p_atom->GetName() << " ";
        sout << left << setw(4) << p_atom->GetType() << " ";
        sout << right << setw(5) << p_atom->GetResidue()->GetIndex()+1 << " ";
        sout << left << setw(4) << p_atom->GetResidue()->GetName() << " ";
        sout << right << fixed << setw(10) << setprecision(6) << p_atom->GetStandardCharge() << " ";
        int tidx = FindAtomTypeIdx(i);
        sout << right << fixed << setw(10) << setprecision(4) << AtomTypes[tidx].mass << " ";
        sout << right << setw(2) << AtomTypes[tidx].z << " ";
        sout << right << fixed << setw(13) << setprecision(7) << AtomTypes[tidx].eps << " ";
        sout << right << fixed << setw(13) << setprecision(7) << AtomTypes[tidx].r0 << endl;
    }
}

//------------------------------------------------------------------------------

int CTopReport::FindAtomTypeIdx(int atidx)
{
    CAmberAtom* p_atom = Topology.AtomList.GetAtom(atidx);

    std::map<int,CAtomType>::iterator tit = AtomTypes.begin();
    std::map<int,CAtomType>::iterator tie = AtomTypes.end();

    while( tit != tie ){
        CAtomType mmtype = tit->second;
        if( mmtype.name == p_atom->GetType() ) return(mmtype.idx);
        tit++;
    }
    return(-1);
}

//------------------------------------------------------------------------------

void CTopReport::WriteBonds(ostream& sout)
{
    sout << "[bonds]" << endl;
    sout << "! 0.5*k(d-d0)^2" << endl;
    sout << "! Index AtomA NameA TypeA AtomB NameB TypeB            d0             K" << endl;

    // print bonds
    for(int i=0; i < Topology.BondList.GetNumberOfBonds(); i++) {
        CAmberBond* p_bond = Topology.BondList.GetBond(i);
        sout << right << setw(7) << i+1 << " ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(p_bond->GetIB());
        sout << right << setw(5) << p_bond->GetIB()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";

        p_atom = Topology.AtomList.GetAtom(p_bond->GetJB());
        sout << right << setw(5) << p_bond->GetJB()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";

        CAmberBondType* p_btype = Topology.BondList.GetBondType(p_bond->GetICB());
        sout << right << fixed << setw(13) << setprecision(6) << p_btype->GetREQ() << " ";
        sout << right << fixed << setw(13) << setprecision(6) << 2.0*p_btype->GetRK() << endl;
    }
}

//------------------------------------------------------------------------------

void CTopReport::WriteAngles(ostream& sout)
{
    sout << "[angles]" << endl;
    sout << "! 0.5*k(a-a0)^2" << endl;
    sout << "! Index AtomA NameA TypeA AtomB NameB TypeB AtomC NameC TypeC            a0             K" << endl;

    // print angles
    for(int i=0; i < Topology.AngleList.GetNumberOfAngles(); i++) {
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(i);
        sout << right << setw(7) << i+1 << " ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(p_angle->GetIT());
        sout << right << setw(5) << p_angle->GetIT()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_angle->GetJT());
        sout << right << setw(5) << p_angle->GetJT()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_angle->GetKT());
        sout << right << setw(5) << p_angle->GetKT()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        CAmberAngleType* p_atype = Topology.AngleList.GetAngleType(p_angle->GetICT());
        sout << right << fixed << setw(13) << setprecision(6) << p_atype->GetTEQ() << " ";
        sout << right << fixed << setw(13) << setprecision(6) << 2.0*p_atype->GetTK() << endl;
    }
}

//------------------------------------------------------------------------------

CDihedralType CTopReport::FindDihedralByTypes(int ia,int ib,int ic,int id)
{
    int iat = FindAtomTypeIdx(ia);
    int ibt = FindAtomTypeIdx(ib);
    int ict = FindAtomTypeIdx(ic);
    int idt = FindAtomTypeIdx(id);

    std::map<int,CDihedralType>::iterator tit = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator tie = DihedralTypes.end();

    while( tit != tie ){
        CDihedralType mmtype = tit->second;
        if( (iat == mmtype.at1)&&(ibt == mmtype.at2)&&(ict == mmtype.at3)&&(idt == mmtype.at4) ) return(mmtype);
        if( (iat == mmtype.at4)&&(ibt == mmtype.at3)&&(ict == mmtype.at2)&&(idt == mmtype.at1) ) return(mmtype);
        tit++;
    }
    return(CDihedralType());
}

//------------------------------------------------------------------------------

void CTopReport::WriteDihedrals(ostream& sout)
{
    // generate list of dihedral types
    int idx = 1;
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++) {
        CAmberDihedral* p_dihedral = Topology.DihedralList.GetDihedral(i);
        int ip = p_dihedral->GetIP();
        int jp = p_dihedral->GetJP();
        int kp = p_dihedral->GetKP();
        int lp = p_dihedral->GetLP();

        // only normal and terminal dihedrals are processed
        if( (p_dihedral->GetType() != 0) && (p_dihedral->GetType() != -1) ) continue;

        CDihedralType dtype;

        dtype = FindDihedralByTypes(ip,jp,kp,lp);

        int icp = p_dihedral->GetICP();
        CAmberDihedralType* p_dihedral_type = Topology.DihedralList.GetDihedralType(icp);
        int pn = abs((int)(p_dihedral_type->GetPN()-1));

        if( dtype.idx == -1 ){
            // new type
            dtype.SetSeriesSize(20);
            dtype.idx = idx;
            idx++;
            dtype.at1 = FindAtomTypeIdx(ip);
            dtype.at2 = FindAtomTypeIdx(jp);
            dtype.at3 = FindAtomTypeIdx(kp);
            dtype.at4 = FindAtomTypeIdx(lp);
        } else {
            // previous type
            if( pn >= dtype.GetSeriesSize() ){
                RUNTIME_ERROR("too short dihedral series");
            }
        }

        // add new data
        dtype.defined[pn] = true;
        dtype.v0[pn]     = p_dihedral_type->GetPK();
        dtype.phase[pn] = p_dihedral_type->GetPHASE();
        dtype.scee  = p_dihedral_type->GetSCEE();
        dtype.scnb  = p_dihedral_type->GetSCNB();

        DihedralTypes[dtype.idx] = dtype;
    }

    sout << "[dihedrals]" << endl;
    sout << "! Index AtomA NameA TypeA AtomB NameB TypeB AtomC NameC TypeC AtomD NameD TypeD pn T            v0         phase" << endl;

    // list unique dihedrals
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++) {
        CAmberDihedral* p_dihedral = Topology.DihedralList.GetDihedral(i);

        sout << right << setw(7) << i+1 << " ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(p_dihedral->GetIP());
        sout << right << setw(5) << p_dihedral->GetIP()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_dihedral->GetJP());
        sout << right << setw(5) << p_dihedral->GetJP()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_dihedral->GetKP());
        sout << right << setw(5) << p_dihedral->GetKP()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_dihedral->GetLP());
        sout << right << setw(5) << p_dihedral->GetLP()+1 << " ";
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        CAmberDihedralType* p_dtype = Topology.DihedralList.GetDihedralType(p_dihedral->GetICP());
        sout << right << setw(2) << (int)p_dtype->GetPN() << " ";
        if( (p_dihedral->GetType() != 0) && (p_dihedral->GetType() != -1) ) {
            sout << "i ";
        } else {
            sout << "  ";
        }
        sout << right << fixed << setw(13) << setprecision(6) << p_dtype->GetPK() << " ";
        sout << right << fixed << setw(13) << setprecision(6) << p_dtype->GetPHASE() << endl;
//        if( (p_dihedral->GetType() != 0) && (p_dihedral->GetType() != -1) ) {
//            sout << endl;
//        } else {
//            sout << right << fixed << setw(13) << setprecision(6) << p_dtype->GetSCEE() << " ";
//            sout << right << fixed << setw(13) << setprecision(6) << p_dtype->GetSCNB() << endl;
//        }
    }
}

//------------------------------------------------------------------------------

void CTopReport::WriteNBList(ostream& sout)
{
    sout << "[nb_list]" << endl;
    sout << "! Index AtomA NameA TypeA AtomB NameB TypeB          scee          scnb" << endl;

    int idx = 1;
    int excidx = 0;
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms() - 1; i++) {
        CAmberAtom* p_atom1 = Topology.AtomList.GetAtom(i);
        int nexcluded = p_atom1->GetNUMEX();
        for(int j=i+1; j < Topology.AtomList.GetNumberOfAtoms(); j++) {
            CAmberAtom* p_atom2 = Topology.AtomList.GetAtom(j);
            // is excluded
            if( nexcluded > 0 ){
                if( Topology.NonBondedList.GetNATEX(excidx) == j ){
                    excidx++;
                    nexcluded--;
                    continue;
                }
            }

            bool found = false;
            // is 1-4?
            CAmberDihedral* p_dih;
            for(int k=0; k < Topology.DihedralList.GetNumberOfDihedrals(); k++) {
                p_dih = Topology.DihedralList.GetDihedral(k);
                if( ((p_dih->GetIP() == i)&&(p_dih->GetLP() == j)) ||
                    ((p_dih->GetIP() == j)&&(p_dih->GetLP() == i)) ){
                        found = true;
                        break;
                }
            }
            int type = 0;
            if( found == true ){
                CDihedralType dtype;
                int ip = p_dih->GetIP();
                int jp = p_dih->GetJP();
                int kp = p_dih->GetKP();
                int lp = p_dih->GetLP();
                dtype = FindDihedralByTypes(ip,jp,kp,lp);
                if( dtype.idx == -1 ){
                    RUNTIME_ERROR("1-4 nb interaction dihedral not found");
                }
                type = dtype.idx;
            }

            sout << right << setw(7) << idx << " ";
            sout << right << setw(5) << i+1 << " ";
            sout << left << setw(5) << p_atom1->GetName() << " ";
            sout << left << setw(5) << p_atom1->GetType() << " ";
            sout << right << setw(5) << j+1 << " ";
            sout << left << setw(5) << p_atom2->GetName() << " ";
            sout << left << setw(5) << p_atom2->GetType() << " ";
            if( type > 0 ){
                sout << right << fixed << setw(13) << setprecision(6) << DihedralTypes[type].scee << " ";
                sout << right << fixed << setw(13) << setprecision(6) << DihedralTypes[type].scnb;
            }
            sout << endl;
            idx++;
        }
        if( nexcluded != 0 ){
            CSmallString error;
            error << "incorrectly processed NB exclusion list, atom " << i << ", remaining " << nexcluded;
            RUNTIME_ERROR(error);
        }

    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



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
#include "Top2STop.hpp"
#include <map>
#include <iomanip>
#include <list>
#include <locale>
#include <FortranMatrix.hpp>
#include <Vector.hpp>
#include <SciLapack.hpp>
#include <boost/algorithm/string.hpp>

#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/graphsym.h>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace OpenBabel;

MAIN_ENTRY(CTop2STop);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CTop2STop::CTop2STop(void)
{
    natoms = 0;
    natom_types = 0;
    nbonds = 0;
    nbond_types = 0;
    nangles = 0;
    nangle_types = 0;
    ndihedrals = 0;
    ndihedral_types = 0;
    ndihedral_seq_size = 0;
    nimpropers = 0;
    nimproper_types = 0;
    nb_size = 0;
    nb_size14 = 0;
    nb_sizeij = 0;
    nnb_types = 0;
    dih_mode = 0;
    dih_samp_freq = 5;
    nsymm_classes = 0;

    Coords.AssignTopology(&Topology);
}

//------------------------------------------------------------------------------

CTop2STop::~CTop2STop(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CTop2STop::Init(int argc,char* argv[])
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
    vout << "# top2stop (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgTopologyName() != "-") {
        vout << "# AMBER topology file (in)       : " << Options.GetArgTopologyName() << endl;
    } else {
        vout << "# AMBER topology file (in)       : - (standard input)" << endl;
    }
    if(Options.GetArgSTopologyName() != "-") {
        vout << "# Simplified topology file (out) : " << Options.GetArgSTopologyName() << endl;
    } else {
        vout << "# Simplified topology file (out) : - (standard output)" << endl;
    }
        vout << "# Dihedral sequence size         : " << Options.GetOptDihedralSeriesSize() << endl;
        vout << "# Dihedral sequence mode         : " << Options.GetOptDihedralMode() << endl;

    return( SO_CONTINUE );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CTop2STop::Run(void)
{
    // load topology
    if( LoadTopology() == false ) return(false);

    if( Options.IsOptCrdNameSet() ){
        if( LoadCoords() == false ) return(false);
    }

    // load filters
    if( Options.IsOptDihedralTypesSet() ){
        if( LoadDihFilters()  == false ) return(false);
    }

    // transform atom types
    TransformTypes();

    // save topology
    if( Options.GetArgSTopologyName() == "-" ){
        WriteAll(cout);
    } else {
        ofstream ifs;
        ifs.open(Options.GetArgSTopologyName());
        if( ! ifs ){
            vout << "<red>>>> ERROR: Unable to open simplified topology: " << Options.GetArgSTopologyName() << "</red>" << endl;
            return(false);
        }
        WriteAll(ifs);
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CTop2STop::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# top2stop terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CTop2STop::LoadTopology(void)
{
    if( Topology.Load(Options.GetArgTopologyName(),true) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER topology: " << Options.GetArgTopologyName() << "</red>" << endl;
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CTop2STop::LoadCoords(void)
{
    if( Coords.Load(Options.GetOptCrdName(),true) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER coords: " << Options.GetOptCrdName() << "</red>" << endl;
        return(false);
    }

    // convert to OBMol
    OBMol  obmol;

    int at_lid = 1;
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++){
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        OBAtom* p_ob_atom = obmol.NewAtom();
        p_ob_atom->SetAtomicNum(p_atom->GetAtomicNumber());
        CPoint pos = Coords.GetPosition(i);
        p_ob_atom->SetVector(pos.x, pos.y, pos.z);
        at_lid++;
    }

    // add bonds
    for(int i=0; i < Topology.BondList.GetNumberOfBonds(); i++){
        CAmberBond*  p_bond = Topology.BondList.GetBond(i);
        int ob_a1_id = p_bond->GetIB() + 1;
        int ob_a2_id = p_bond->GetJB() + 1;
        int order = 1;
        // create bond
        obmol.AddBond(ob_a1_id, ob_a2_id, order);
    }

    // generate symmetry classes
    OBGraphSym  graph_sym(&obmol);
    graph_sym.GetSymmetry(SymmClasses);

    return(true);
}

//------------------------------------------------------------------------------

bool CTop2STop::LoadDihFilters(void)
{
    vout << endl;
    vout << "# Loading dihedral type filters from: " << Options.GetOptDihedralTypes() << endl;
    ifstream ifs(Options.GetOptDihedralTypes());
    if( !ifs ){
        vout << "<red>>>> ERROR: Unable to open dihedral type filters: " << Options.GetOptDihedralTypes() << "</red>" << endl;
        return(false);
    }

    vout << endl;
    vout << "# TA   TB   TC   TD" << endl;
    vout << "# -- ---- ---- ----" << endl;

    string line;
    while(getline(ifs,line)){
        stringstream sfs(line);
        CDihedralTypeFilter filter;
        sfs >> filter.t1 >> filter.t2 >> filter.t3 >> filter.t4;
        if( filter.t3.empty() &&  filter.t3.empty() ){
            vout << setw(4) << right << "X" << " " << setw(4) << right << filter.t1 << " ";
            vout << setw(4) << right << filter.t2 << " " << setw(4) << right << "X";
            filter.full = false;
        } else {
            vout        << setw(4) << right << filter.t1 << " " << setw(4) << right << filter.t2;
            vout << " " << setw(4) << right << filter.t3 << " " << setw(4) << right << filter.t4;
            filter.full = true;
        }
        vout << endl;
        DihFilters.push_back(filter);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CTop2STop::WriteAll(ostream& sout)
{
    WriteAtomTypes(sout);
    WriteAtoms(sout);

    WriteBondTypes(sout);
    WriteBonds(sout);

    WriteAngleTypes(sout);
    WriteAngles(sout);

    WriteDihedralTypes(sout);
    WriteDihedrals(sout);

    WriteImproperTypes(sout);
    WriteImpropers(sout);

    WriteNBTypes(sout);
    if( Options.GetOptRebuildNBList() == true ){
        WriteNBListRebuild(sout);
    } else {
        WriteNBListKeep(sout);
    }

    WriteDimensions(sout);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CTop2STop::TransformTypes(void)
{
    if( Options.GetOptTransform() == "none" ) return;

    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);

        std:string type = p_atom->GetType();

        if( Options.GetOptTransform() == "both" ){
            if( type.size() > 0 ){
                if( isupper(type[0]) == true ){
                    type[0] = tolower(type[0]);
                } else {
                    type[0] = toupper(type[0]);
                }
            }
            if( type.size() > 1 ){
                if( isupper(type[1]) == true ){
                    type[1] = tolower(type[1]);
                } else {
                    type[1] = toupper(type[1]);
                }
            }
        }

        if( Options.GetOptTransform() == "first" ){
            if( type.size() > 0 ){
                if( isupper(type[0]) == true ){
                    type[0] = tolower(type[0]);
                } else {
                    type[0] = toupper(type[0]);
                }
            }
        }

        if( Options.GetOptTransform() == "second" ){
            if( type.size() > 1 ){
                if( isupper(type[1]) == true ){
                    type[1] = tolower(type[1]);
                } else {
                    type[1] = toupper(type[1]);
                }
            }
        }

        p_atom->SetType(type.c_str());
    }
}

//------------------------------------------------------------------------------

void CTop2STop::WriteAtomTypes(ostream& sout)
{
    std::map<std::string,CAtomType> unique_types;

    // extract atom types
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        CAtomType mmtype;
        mmtype.name = p_atom->GetType();
        mmtype.mass = p_atom->GetMass();
        mmtype.IAC = p_atom->GetIAC();
        if( p_atom->GetAtomicNumber() > 0 ){
            mmtype.z = p_atom->GetAtomicNumber();
        } else {
            mmtype.z = p_atom->GuessZ();
        }
        if( unique_types[mmtype.name].idx == -1 ){

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

    // write types
    sout << "[atom_types]" << endl;
    sout << "! Index Type       Mass  Z  " << endl;

    std::map<int,CAtomType>::iterator tit = AtomTypes.begin();
    std::map<int,CAtomType>::iterator tie = AtomTypes.end();

    while( tit != tie ){
        CAtomType mmtype = tit->second;
        sout << right << setw(7) << mmtype.idx << " ";
        sout << left << setw(4) << mmtype.name << " ";
        sout << right << fixed << setw(10) << setprecision(4) << mmtype.mass << " ";
        sout << right << setw(2) << mmtype.z <<  endl;
        tit++;
    }

    natom_types = AtomTypes.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteAtoms(ostream& sout)
{
    sout << "[atoms]" << endl;
    sout << "! Index Type Name ResID ResN     Charge SymmClass ! Type" << endl;

    // write atoms
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        sout << right << setw(7) << i+1 << " ";
        sout << right << setw(4) << FindAtomTypeIdx(i) << " ";
        sout << left  << setw(4) << p_atom->GetName() << " ";
        sout << right << setw(5) << p_atom->GetResidue()->GetIndex()+1 << " ";
        sout << left  << setw(4) << p_atom->GetResidue()->GetName() << " ";
        sout << right << fixed   << setw(10) << setprecision(6) << p_atom->GetStandardCharge() << " ";
        int symm_class = i + 1;
        if( SymmClasses.size() != 0 ){
            symm_class = SymmClasses[i];
        }
        sout << right << setw(9) << symm_class << " ! ";
        sout << left << setw(4) << p_atom->GetType() << endl;
    }

    natoms = Topology.AtomList.GetNumberOfAtoms();
}

//------------------------------------------------------------------------------

int CTop2STop::FindAtomTypeIdx(int atidx)
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

void CTop2STop::WriteBondTypes(ostream& sout)
{
    // generate list of bond types
    for(int i=0; i < Topology.BondList.GetNumberOfBonds(); i++) {
        CAmberBond* p_bond = Topology.BondList.GetBond(i);
        CAmberBondType* p_btype = Topology.BondList.GetBondType(p_bond->GetICB());
        CBondType btype;
        btype.idx = p_bond->GetICB()+1;
        btype.at1 = FindAtomTypeIdx(p_bond->GetIB());
        btype.at2 = FindAtomTypeIdx(p_bond->GetJB());
        btype.form = 1; // always harmonic
        btype.d0 = p_btype->GetREQ();
        btype.k = p_btype->GetRK();
        if( BondTypes[btype.idx].idx == -1 ){
            BondTypes[btype.idx] = btype;
        } else {
            if( BondTypes[btype.idx] != btype ){
                RUNTIME_ERROR("bond type mismatch");
            }
        }
    }

    sout << "[bond_types]" << endl;
    sout << "! 0.5*k(d-d0)^2" << endl;
    sout << "! Index TypeA TypeB Form            d0             K ! TypeA TypeB" << endl;

    std::map<int,CBondType>::iterator it = BondTypes.begin();
    std::map<int,CBondType>::iterator ie = BondTypes.end();

    while( it != ie ){
        CBondType btype = it->second;
        sout << right << setw(7) << btype.idx << " ";
        sout << right << setw(5) << btype.at1 << " ";
        sout << right << setw(5) << btype.at2 << " ";
        sout << right << setw(4) << btype.form << " ";
        sout << right << fixed << setw(13) << setprecision(6) << btype.d0 << " ";
        sout << right << fixed << setw(13) << setprecision(6) << 2.0*btype.k << " ! ";
        sout << left << setw(5) << AtomTypes[btype.at1].name << " ";
        sout << left << setw(5) << AtomTypes[btype.at2].name << endl;
        it++;
    }

    nbond_types = BondTypes.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteBonds(ostream& sout)
{
    sout << "[bonds]" << endl;
    sout << "! Index AtomA AtomB Type ! AtomA TypeA AtomB TypeB" << endl;

    // print bonds
    for(int i=0; i < Topology.BondList.GetNumberOfBonds(); i++) {
        CAmberBond* p_bond = Topology.BondList.GetBond(i);
        sout << right << setw(7) << i+1 << " ";
        sout << right << setw(5) << p_bond->GetIB()+1 << " ";
        sout << right << setw(5) << p_bond->GetJB()+1 << " ";
        sout << right << setw(4) << p_bond->GetICB()+1 << " ! ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(p_bond->GetIB());
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_bond->GetJB());
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << endl;
    }

    nbonds = Topology.BondList.GetNumberOfBonds();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteAngleTypes(ostream& sout)
{
    // generate list of angle types
    for(int i=0; i < Topology.AngleList.GetNumberOfAngles(); i++) {
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(i);
        CAmberAngleType* p_atype = Topology.AngleList.GetAngleType(p_angle->GetICT());
        CAngleType atype;
        atype.idx = p_angle->GetICT()+1;
        atype.at1 = FindAtomTypeIdx(p_angle->GetIT());
        atype.at2 = FindAtomTypeIdx(p_angle->GetJT());
        atype.at3 = FindAtomTypeIdx(p_angle->GetKT());
        atype.form = 1; // always harmonic
        atype.a0 = p_atype->GetTEQ();
        atype.k = p_atype->GetTK();
        if( AngleTypes[atype.idx].idx == -1 ){
            AngleTypes[atype.idx] = atype;
        } else {
            if( AngleTypes[atype.idx] != atype ){
                RUNTIME_ERROR("angle type mismatch");
            }
        }
    }

    sout << "[angle_types]" << endl;
    sout << "! 0.5*k(a-a0)^2" << endl;
    sout << "! Index TypeA TypeB TypeC Form            a0             K ! TypeA TypeB TypeC" << endl;

    std::map<int,CAngleType>::iterator it = AngleTypes.begin();
    std::map<int,CAngleType>::iterator ie = AngleTypes.end();

    while( it != ie ){
        CAngleType atype = it->second;
        sout << right << setw(7) << atype.idx << " ";
        sout << right << setw(5) << atype.at1 << " ";
        sout << right << setw(5) << atype.at2 << " ";
        sout << right << setw(5) << atype.at3 << " ";
        sout << right << setw(4) << atype.form << " ";
        sout << right << fixed << setw(13) << setprecision(6) << atype.a0 << " ";
        sout << right << fixed << setw(13) << setprecision(6) << 2.0*atype.k << " ! ";
        sout << left << setw(5) << AtomTypes[atype.at1].name << " ";
        sout << left << setw(5) << AtomTypes[atype.at2].name << " ";
        sout << left << setw(5) << AtomTypes[atype.at3].name << endl;
        it++;
    }

    nangle_types = AngleTypes.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteAngles(ostream& sout)
{
    sout << "[angles]" << endl;
    sout << "! Index AtomA AtomB AtomC Type ! AtomA TypeA AtomB TypeB AtomC TypeC" << endl;

    // print angles
    for(int i=0; i < Topology.AngleList.GetNumberOfAngles(); i++) {
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(i);
        sout << right << setw(7) << i+1 << " ";
        sout << right << setw(5) << p_angle->GetIT()+1 << " ";
        sout << right << setw(5) << p_angle->GetJT()+1 << " ";
        sout << right << setw(5) << p_angle->GetKT()+1 << " ";
        sout << right << setw(4) << p_angle->GetICT()+1 << " ! ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(p_angle->GetIT());
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_angle->GetJT());
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(p_angle->GetKT());
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << endl;
    }

    nangles = Topology.AngleList.GetNumberOfAngles();
}

//------------------------------------------------------------------------------

CDihedralType CTop2STop::FindDihedralByTypes(int ia,int ib,int ic,int id)
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

CDihedralType CTop2STop::FindImproperByTypes(int ia,int ib,int ic,int id)
{
    int iat = FindAtomTypeIdx(ia);
    int ibt = FindAtomTypeIdx(ib);
    int ict = FindAtomTypeIdx(ic);
    int idt = FindAtomTypeIdx(id);

    std::map<int,CDihedralType>::iterator tit = ImproperTypes.begin();
    std::map<int,CDihedralType>::iterator tie = ImproperTypes.end();

    while( tit != tie ){
        CDihedralType mmtype = tit->second;
        if( (iat == mmtype.at1)&&(ibt == mmtype.at2)&&(ict == mmtype.at3)&&(idt == mmtype.at4) ) return(mmtype);
        if( (iat == mmtype.at4)&&(ibt == mmtype.at3)&&(ict == mmtype.at2)&&(idt == mmtype.at1) ) return(mmtype);
        tit++;
    }
    return(CDihedralType());
}

//------------------------------------------------------------------------------

void CTop2STop::WriteDihedralTypes(ostream& sout)
{

    ndihedral_seq_size = Options.GetOptDihedralSeriesSize();

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
            dtype.SetSeriesSize(ndihedral_seq_size);
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

    dih_mode = 0;
    if( Options.GetOptDihedralMode() == "cos" ){
        dih_mode = 1;
    } else if( Options.GetOptDihedralMode() == "grbf" ){
        dih_mode = 2;
    }
    if( dih_mode == 0 ){
        RUNTIME_ERROR("unsupported dihedral mode - p1");
    }

    ndihedral_types = DihedralTypes.size();

    if( dih_mode == 2 ){
        TransformCosToGRBF();
    }

    sout << "[dihedral_types]" << endl;
    sout << "! Index TypeA TypeB TypeC TypeD Form          scee          scnb ! TypeA TypeB TypeC TypeD" << endl;

    std::map<int,CDihedralType>::iterator it = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator ie = DihedralTypes.end();

    while( it != ie ){
        CDihedralType dtype = it->second;
        sout << right << setw(7) << dtype.idx << " ";
        sout << right << setw(5) << dtype.at1 << " ";
        sout << right << setw(5) << dtype.at2 << " ";
        sout << right << setw(5) << dtype.at3 << " ";
        sout << right << setw(5) << dtype.at4 << " ";
        if( dtype.grbf ){
            sout << right << setw(4) << 2 << " ";
        } else {
            sout << right << setw(4) << 1 << " ";
        }
        sout << right << fixed << setw(13) << setprecision(6) << dtype.scee << " ";
        sout << right << fixed << setw(13) << setprecision(6) << dtype.scnb << " ! ";
        sout << left << setw(5) << AtomTypes[dtype.at1].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at2].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at3].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at4].name << endl;
        it++;
    }

    int costypes = 0;
    int grbftypes = 0;
    std::map<int,CDihedralType>::iterator tit = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator tie = DihedralTypes.end();

    while( tit != tie ){
        CDihedralType dtype = tit->second;
        if( dtype.grbf ){
            grbftypes++;
        } else {
            costypes++;
        }
        tit++;
    }

    if( costypes > 0 ) WriteDihedralSeqCosMode(sout);
    if( grbftypes > 0 ) WriteDihedralSeqGRBFMode(sout);
}

//------------------------------------------------------------------------------

void CTop2STop::WriteDihedralSeqCosMode(ostream& sout)
{
    sout << "[dihedral_seq_cos]" << endl;
    sout << "! Type pn            v0         phase defined" << endl;

    std::map<int,CDihedralType>::iterator it = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator ie = DihedralTypes.end();

    while( it != ie ){
        CDihedralType dtype = it->second;
        if( dtype.grbf == false ){
            for(int i=0; i < dtype.GetSeriesSize(); i++ ){
                sout << right << setw(6) << dtype.idx << " ";
                sout << right << setw(2) << i+1 << " ";
                if( Options.GetOptZeroDihPhase() ) {
                    if( fabs(dtype.phase[i] - M_PI) < 0.1 ) {
                        dtype.v0[i] = - dtype.v0[i];
                        dtype.phase[i] = dtype.phase[i] - M_PI;
                    }

                }
                sout << right << fixed << setw(13) << setprecision(6) << dtype.v0[i] << " ";
                sout << right << fixed << setw(13) << setprecision(6) << dtype.phase[i];
                if( dtype.defined[i] ){
                    sout << "       1";
                } else {
                    sout << "       0";
                }
                sout << endl;
            }
        }
        it++;
    }
}

//------------------------------------------------------------------------------

void CTop2STop::TransformCosToGRBF(void)
{
    vout << endl;
    vout << "Transforming dihedral cos series to rgbf series ..." << endl;
    vout << "   Number of series (dihedral types) = " << ndihedral_types << endl;
    vout << "   Series size                       = " << ndihedral_seq_size << endl;
    vout << "   Number of training points         = " << (ndihedral_seq_size+1)*dih_samp_freq << endl;

    std::map<int,CDihedralType>::iterator it = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator ie = DihedralTypes.end();

    vout << endl;
    while(it != ie){
        CDihedralType type =  it->second;

        bool ok = false;
        if( DihFilters.size() > 0 ){
            for(size_t i=0; i < DihFilters.size(); i++ ){
                if( DihFilters[i].full ){
                    if( ((trim_copy(AtomTypes[DihedralTypes[type.idx].at1].name) == DihFilters[i].t1) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at2].name) == DihFilters[i].t2) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at3].name) == DihFilters[i].t3) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at4].name) == DihFilters[i].t4)) ||
                        ((trim_copy(AtomTypes[DihedralTypes[type.idx].at1].name) == DihFilters[i].t4) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at2].name) == DihFilters[i].t3) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at3].name) == DihFilters[i].t2) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at4].name) == DihFilters[i].t1)) ){
                        ok = true;
                        break;
                    }
                } else {
                    if( ((trim_copy(AtomTypes[DihedralTypes[type.idx].at2].name) == DihFilters[i].t1) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at3].name) == DihFilters[i].t2)) ||
                        ((trim_copy(AtomTypes[DihedralTypes[type.idx].at2].name) == DihFilters[i].t2) &&
                         (trim_copy(AtomTypes[DihedralTypes[type.idx].at3].name) == DihFilters[i].t1)) ){
                        ok = true;
                        break;
                    }
                }
            }
        } else {
            ok = true;
        }

        if( ok ){
            SolveTransformation(type.idx);
            DihedralTypes[type.idx].grbf = true;
        }
        it++;
    }
}

//------------------------------------------------------------------------------

void CTop2STop::SolveTransformation(int type)
{
    vout << "   fitting dihedral type " << setw(4) << type << " ... final error = ";

    // transform
    DihedralTypes[type].DihCOffset = Options.GetOptDihCOffset();
    DihedralTypes[type].Cos2GRBF(dih_samp_freq);

    // calculate rmse
    double rmse = DihedralTypes[type].RMSECos2GRBF(dih_samp_freq);

    vout << fixed << setw(13) << setprecision(8) << rmse << endl;
}

//------------------------------------------------------------------------------

void CTop2STop::WriteDihedralSeqGRBFMode(ostream& sout)
{
    sout << "[dihedral_seq_grbf]" << endl;
    sout << "! Type pn             c             p            w2" << endl;

    std::map<int,CDihedralType>::iterator it = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator ie = DihedralTypes.end();

    while( it != ie ){
        CDihedralType dtype = it->second;
        if( dtype.grbf == true ){
            for(int i=0; i < dtype.GetSeriesSize(); i++ ){
                sout << right << setw(6) << dtype.idx << " ";
                sout << right << setw(2) << i+1 << " ";
                sout << right << fixed << setw(13) << setprecision(6) << dtype.c[i] << " ";
                sout << right << fixed << setw(13) << setprecision(6) << dtype.p[i] << " ";
                sout << right << fixed << setw(13) << setprecision(6) << dtype.w2[i];
                sout << endl;
            }
        }
        it++;
    }
}

//------------------------------------------------------------------------------

void CTop2STop::WriteDihedrals(ostream& sout)
{
    // list unique dihedrals
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++) {
        CAmberDihedral* p_dihedral = Topology.DihedralList.GetDihedral(i);
        int ip = p_dihedral->GetIP();
        int jp = p_dihedral->GetJP();
        int kp = p_dihedral->GetKP();
        int lp = p_dihedral->GetLP();

        // only normal and terminal dihedrals are processed
        if( (p_dihedral->GetType() != 0) && (p_dihedral->GetType() != -1) ) continue;

        std::list<CDihedral>::iterator  it = UniqueDihedrals.begin();
        std::list<CDihedral>::iterator  ie = UniqueDihedrals.end();

        bool found = false;
        while( it != ie ){
            CDihedral dih = *it;
            if( ((dih.at1 == ip)&&(dih.at2 == jp)&&(dih.at3 == kp)&&(dih.at4 == lp)) ||
                ((dih.at4 == ip)&&(dih.at3 == jp)&&(dih.at2 == kp)&&(dih.at1 == lp)) ){
                    found = true;
                    break;
                }
            it++;
        }
        if( found == true ) continue;
        CDihedral dih;
        dih.at1 = ip;
        dih.at2 = jp;
        dih.at3 = kp;
        dih.at4 = lp;
        UniqueDihedrals.push_back(dih);
    }

    sout << "[dihedrals]" << endl;
    sout << "! Index AtomA AtomB AtomC AtomD Type ! AtomA TypeA AtomB TypeB AtomC TypeC AtomD TypeD" << endl;

    std::list<CDihedral>::iterator  it = UniqueDihedrals.begin();
    std::list<CDihedral>::iterator  ie = UniqueDihedrals.end();

    int idx = 1;
    while( it != ie ){
        CDihedral dih = *it;
        int ip = dih.at1;
        int jp = dih.at2;
        int kp = dih.at3;
        int lp = dih.at4;
        CDihedralType dtype;
        dtype = FindDihedralByTypes(ip,jp,kp,lp);
        if( dtype.idx == -1 ){
            RUNTIME_ERROR("dihedral not defined");
        }

        sout << right << setw(7) << idx << " ";
        sout << right << setw(5) << ip+1 << " ";
        sout << right << setw(5) << jp+1 << " ";
        sout << right << setw(5) << kp+1 << " ";
        sout << right << setw(5) << lp+1 << " ";
        sout << right << setw(4) << dtype.idx << " ! ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(ip);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(jp);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(kp);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(lp);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << endl;
        it++;
        idx++;
    }

    ndihedrals = UniqueDihedrals.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteImproperTypes(ostream& sout)
{
    // generate list of dihedral types
    int idx = 1;
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++) {
        CAmberDihedral* p_dihedral = Topology.DihedralList.GetDihedral(i);
        int ip = p_dihedral->GetIP();
        int jp = p_dihedral->GetJP();
        int kp = p_dihedral->GetKP();
        int lp = p_dihedral->GetLP();

        // only improper and terminal improper dihedrals are processed
        if( (p_dihedral->GetType() != 1)&&(p_dihedral->GetType() != -2) ) continue;

        CDihedralType dtype;

        dtype = FindImproperByTypes(ip,jp,kp,lp);

        int icp = p_dihedral->GetICP();
        CAmberDihedralType* p_dihedral_type = Topology.DihedralList.GetDihedralType(icp);

        if( dtype.idx == -1 ){
            // new type
            dtype.SetSeriesSize(1);
            dtype.idx = idx;
            idx++;
            dtype.at1 = FindAtomTypeIdx(ip);
            dtype.at2 = FindAtomTypeIdx(jp);
            dtype.at3 = FindAtomTypeIdx(kp);
            dtype.at4 = FindAtomTypeIdx(lp);
        }

        // add new data
        dtype.v0[0]     = p_dihedral_type->GetPK();
        dtype.phase[0]  = p_dihedral_type->GetPHASE();
        dtype.scee      = p_dihedral_type->GetSCEE();
        dtype.scnb      = p_dihedral_type->GetSCNB();

        ImproperTypes[dtype.idx] = dtype;
    }

    sout << "[improper_types]" << endl;
    sout << "! Index TypeA TypeB TypeC TypeD            v0         phase ! TypeA TypeB TypeC TypeD" << endl;

    std::map<int,CDihedralType>::iterator it = ImproperTypes.begin();
    std::map<int,CDihedralType>::iterator ie = ImproperTypes.end();

    while( it != ie ){
        CDihedralType dtype = it->second;
        sout << right << setw(7) << dtype.idx << " ";
        sout << right << setw(5) << dtype.at1 << " ";
        sout << right << setw(5) << dtype.at2 << " ";
        sout << right << setw(5) << dtype.at3 << " ";
        sout << right << setw(5) << dtype.at4 << " ";
        sout << right << fixed << setw(13) << setprecision(6) << dtype.v0[0] << " ";
        sout << right << fixed << setw(13) << setprecision(6) << dtype.phase[0] << " ! ";
        sout << left << setw(5) << AtomTypes[dtype.at1].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at2].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at3].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at4].name << endl;
        it++;
    }

    nimproper_types = ImproperTypes.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteImpropers(ostream& sout)
{

    // list unique dihedrals
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++) {
        CAmberDihedral* p_dihedral = Topology.DihedralList.GetDihedral(i);
        int ip = p_dihedral->GetIP();
        int jp = p_dihedral->GetJP();
        int kp = p_dihedral->GetKP();
        int lp = p_dihedral->GetLP();

        // only improper and terminal improper dihedrals are processed
        if( (p_dihedral->GetType() != 1)&&(p_dihedral->GetType() != -2) ) continue;

        std::list<CDihedral>::iterator  it = UniqueImpropers.begin();
        std::list<CDihedral>::iterator  ie = UniqueImpropers.end();

        bool found = false;
        while( it != ie ){
            CDihedral dih = *it;
            if( ((dih.at1 == ip)&&(dih.at2 == jp)&&(dih.at3 == kp)&&(dih.at4 == lp)) ||
                ((dih.at4 == ip)&&(dih.at3 == jp)&&(dih.at2 == kp)&&(dih.at1 == lp)) ){
                    found = true;
                    break;
                }
            it++;
        }
        if( found == true ) continue;
        CDihedral dih;
        dih.at1 = ip;
        dih.at2 = jp;
        dih.at3 = kp;
        dih.at4 = lp;
        UniqueImpropers.push_back(dih);
    }

    sout << "[impropers]" << endl;
    sout << "! Index AtomA AtomB AtomC AtomD Type ! AtomA TypeA AtomB TypeB AtomC TypeC AtomD TypeD" << endl;

    std::list<CDihedral>::iterator  it = UniqueImpropers.begin();
    std::list<CDihedral>::iterator  ie = UniqueImpropers.end();

    int idx = 1;
    while( it != ie ){
        CDihedral dih = *it;
        int ip = dih.at1;
        int jp = dih.at2;
        int kp = dih.at3;
        int lp = dih.at4;
        CDihedralType dtype;
        dtype = FindImproperByTypes(ip,jp,kp,lp);
        if( dtype.idx == -1 ){
            RUNTIME_ERROR("improper not defined");
        }

        sout << right << setw(7) << idx << " ";
        sout << right << setw(5) << ip+1 << " ";
        sout << right << setw(5) << jp+1 << " ";
        sout << right << setw(5) << kp+1 << " ";
        sout << right << setw(5) << lp+1 << " ";
        sout << right << setw(4) << dtype.idx << " ! ";
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(ip);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(jp);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(kp);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << " ";
        p_atom = Topology.AtomList.GetAtom(lp);
        sout << left << setw(5) << p_atom->GetName() << " ";
        sout << left << setw(5) << p_atom->GetType() << endl;
        it++;
        idx++;
    }

    nimpropers = UniqueImpropers.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteNBTypes(ostream& sout)
{
    // Topology.NonBondedList.GetNumberOfTypes() - contains only unique NB types not all!!!
    // thus we must rebuild the list ...

    NBTypes.clear();

    sout << "[nb_types]" << endl;
    sout << "! Index TypA TypB           eps            r0 ! TypeA TypeB" << endl;

    int idx = 1;
    for(int i=1; i <= natom_types; i++){
        for(int j=i; j <= natom_types; j++){

            int iac_a = AtomTypes[i].IAC-1;
            int iac_b = AtomTypes[j].IAC-1;

            // ij parameters
            int icoij = Topology.NonBondedList.GetICOIndex(Topology.NonBondedList.GetNumberOfTypes()*iac_a + iac_b);

            NBTypes[(i-1)*natom_types+j] = idx;
            NBTypes[(j-1)*natom_types+i] = idx;

            double aij,bij,epsij,rij;

            aij = Topology.NonBondedList.GetAParam(icoij);
            bij = Topology.NonBondedList.GetBParam(icoij);

            epsij = 0.0;
            rij = 0.0;
            if( aij != 0.0 ) {
                epsij = bij*bij / (4.0 * aij);
                rij = pow(2*aij/bij,1.0/6.0) * 0.5;
            }

            int it = AtomTypes[i].idx;
            int jt = AtomTypes[j].idx;

            sout << right << setw(7) << idx << " ";
            sout << left << setw(4) << it << " ";
            sout << left << setw(4) << jt << " ";
            sout << right << fixed << setw(13) << setprecision(7) << epsij << " ";
            sout << right << fixed << setw(13) << setprecision(7) << rij*2.0 << " ! ";
            sout << left << setw(5) << AtomTypes[it].name << " ";
            sout << left << setw(5) << AtomTypes[jt].name << endl;
            idx++;
        }
    }

    nnb_types = idx-1;
}

//------------------------------------------------------------------------------

// keep existing exclusion list from AMBER topology

void CTop2STop::WriteNBListKeep(ostream& sout)
{
    sout << "[nb_list]" << endl;
    sout << "! Index AtomA AtomB  Type Dihed ! AtomA TypeA AtomB TypeB" << endl;

    nb_size = 0;
    nb_size14 = 0;
    nb_sizeij = 0;
    int idx = 1;
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom1 = Topology.AtomList.GetAtom(i);
        for(int j=i+1; j < Topology.AtomList.GetNumberOfAtoms(); j++) {
            CAmberAtom* p_atom2 = Topology.AtomList.GetAtom(j);

            // is exluded
            if( Topology.IsNBPairExcluded(i,j) == true ) continue;

            int ati = FindAtomTypeIdx(i);
            int atj = FindAtomTypeIdx(j);
            int nbidx = NBTypes[(ati-1)*natom_types+atj];

            sout << right << setw(7) << idx << " ";
            sout << right << setw(5) << i+1 << " ";
            sout << right << setw(5) << j+1 << " ";
            sout << right << setw(5) << nbidx << " ";
            sout << right << setw(5) << 0 << " ! ";
            sout << left << setw(5) << p_atom1->GetName() << " ";
            sout << left << setw(5) << p_atom1->GetType() << " ";
            sout << left << setw(5) << p_atom2->GetName() << " ";
            sout << left << setw(5) << p_atom2->GetType() << endl;
            nb_size++;
            nb_sizeij++;
            idx++;
        }
    }

    // go through unique list of dihedrals to generate 1-4 NB interactions

    std::list<CDihedral>::iterator  it = UniqueDihedrals.begin();
    std::list<CDihedral>::iterator  ie = UniqueDihedrals.end();

    list<CNBPair>   NBPairs14;

    while( it != ie ){
        CDihedral dih = *it;
        it++;
        int ip = dih.at1;
        int jp = dih.at2;
        int kp = dih.at3;
        int lp = dih.at4;

        bool found = false;

        // is it already in the list (6-membered rings)
        std::list<CNBPair>::iterator  nit = NBPairs14.begin();
        std::list<CNBPair>::iterator  nie = NBPairs14.end();

        while( nit != nie ){
            CNBPair pair = *nit;
            if( ( (pair.at1 == ip) && (pair.at2 == lp)) || ((pair.at1 == lp) && (pair.at2 == ip)) ) {
                found = true;
                break;
            }
            nit++;
        }
        if( found == true ) continue;

        // exclude 1-2 and 1-3 bonded (in 4-and 5-membered rings)

        // is bonded?
        for(int k=0; k < Topology.BondList.GetNumberOfBonds(); k++) {
            CAmberBond* p_bond = Topology.BondList.GetBond(k);
            if( ((p_bond->GetIB() == ip)&&(p_bond->GetJB() == lp)) ||
                ((p_bond->GetIB() == lp)&&(p_bond->GetJB() == ip)) ){
                    found = true;
                    break;
            }
        }
        if( found == true ) continue;

        // is 1-3?
        for(int k=0; k < Topology.AngleList.GetNumberOfAngles(); k++) {
            CAmberAngle* p_angle = Topology.AngleList.GetAngle(k);
            if( ((p_angle->GetIT() == ip)&&(p_angle->GetKT() == lp)) ||
                ((p_angle->GetIT() == lp)&&(p_angle->GetKT() == ip)) ){
                    found = true;
                    break;
            }
        }
        if( found == true ) continue;

        // get dihedral type
        CDihedralType dtype;
        dtype = FindDihedralByTypes(ip,jp,kp,lp);
        if( dtype.idx == -1 ){
            RUNTIME_ERROR("1-4 nb interaction dihedral not found");
        }
        int didx = dtype.idx;

        CAmberAtom* p_atom1 = Topology.AtomList.GetAtom(ip);
        CAmberAtom* p_atom2 = Topology.AtomList.GetAtom(lp);

        int ati = FindAtomTypeIdx(ip);
        int atj = FindAtomTypeIdx(lp);
        int nbidx = NBTypes[(ati-1)*natom_types+atj];

        sout << right << setw(7) << idx << " ";
        sout << right << setw(5) << ip+1 << " ";
        sout << right << setw(5) << lp+1 << " ";
        sout << right << setw(5) << nbidx << " ";
        sout << right << setw(5) << didx << " ! ";
        sout << left << setw(5) << p_atom1->GetName() << " ";
        sout << left << setw(5) << p_atom1->GetType() << " ";
        sout << left << setw(5) << p_atom2->GetName() << " ";
        sout << left << setw(5) << p_atom2->GetType() << endl;

        CNBPair pair;
        pair.at1 = ip;
        pair.at2 = lp;
        NBPairs14.push_back(pair);

        nb_size++;
        nb_size14++;
        idx++;
    }
}

//------------------------------------------------------------------------------

// rebuild exclusion list

void CTop2STop::WriteNBListRebuild(ostream& sout)
{
    sout << "[nb_list]" << endl;
    sout << "! Index AtomA AtomB  Type Dihed ! AtomA TypeA AtomB TypeB" << endl;

    nb_size = 0;
    nb_size14 = 0;
    nb_sizeij = 0;
    int idx = 1;
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom1 = Topology.AtomList.GetAtom(i);
        for(int j=i+1; j < Topology.AtomList.GetNumberOfAtoms(); j++) {
            CAmberAtom* p_atom2 = Topology.AtomList.GetAtom(j);

            bool found = false;
            // is bonded?
            for(int k=0; k < Topology.BondList.GetNumberOfBonds(); k++) {
                CAmberBond* p_bond = Topology.BondList.GetBond(k);
                if( ((p_bond->GetIB() == i)&&(p_bond->GetJB() == j)) ||
                    ((p_bond->GetIB() == j)&&(p_bond->GetJB() == i)) ){
                        found = true;
                        break;
                }
            }
            if( found == true ) continue;

            // is 1-3?
            for(int k=0; k < Topology.AngleList.GetNumberOfAngles(); k++) {
                CAmberAngle* p_angle = Topology.AngleList.GetAngle(k);
                if( ((p_angle->GetIT() == i)&&(p_angle->GetKT() == j)) ||
                    ((p_angle->GetIT() == j)&&(p_angle->GetKT() == i)) ){
                        found = true;
                        break;
                }
            }
            if( found == true ) continue;

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
            int didx = 0;
            if( (found == true) && ( (p_dih->GetType() == 0) || (p_dih->GetType() == -1) ) ){
                CDihedralType dtype;
                int ip = p_dih->GetIP();
                int jp = p_dih->GetJP();
                int kp = p_dih->GetKP();
                int lp = p_dih->GetLP();
                dtype = FindDihedralByTypes(ip,jp,kp,lp);
                if( dtype.idx == -1 ){
                    RUNTIME_ERROR("1-4 nb interaction dihedral not found");
                }
                didx = dtype.idx;
                nb_size14++;
            } else {
                nb_sizeij++;
            }

            int ati = FindAtomTypeIdx(i);
            int atj = FindAtomTypeIdx(j);
            int nbidx = NBTypes[(ati-1)*natom_types+atj];

            sout << right << setw(7) << idx << " ";
            sout << right << setw(5) << i+1 << " ";
            sout << right << setw(5) << j+1 << " ";
            sout << right << setw(5) << nbidx << " ";
            sout << right << setw(5) << didx << " ! ";
            sout << left << setw(5) << p_atom1->GetName() << " ";
            sout << left << setw(5) << p_atom1->GetType() << " ";
            sout << left << setw(5) << p_atom2->GetName() << " ";
            sout << left << setw(5) << p_atom2->GetType() << endl;
            nb_size++;
            idx++;
        }
    }
}

//------------------------------------------------------------------------------

void CTop2STop::WriteDimensions(std::ostream& sout)
{
    sout << "[dimensions]" << endl;
    sout << "! Scope           Size" << endl;
    sout << "atoms             " << natoms << endl;
    sout << "atom_types        " << natom_types << endl;
    sout << "bonds             " << nbonds << endl;
    sout << "bond_types        " << nbond_types << endl;
    sout << "angles            " << nangles << endl;
    sout << "angle_types       " << nangle_types << endl;
    sout << "dihedrals         " << ndihedrals << endl;
    sout << "dihedral_types    " << ndihedral_types << endl;
    sout << "dihedral_seq_size " << ndihedral_seq_size << endl;
    sout << "impropers         " << nimpropers << endl;
    sout << "improper_types    " << nimproper_types << endl;
    sout << "nb_sizeij         " << nb_sizeij << endl;
    sout << "nb_size14         " << nb_size14 << endl;
    sout << "nb_size           " << nb_size << endl;
    sout << "nb_types          " << nnb_types << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



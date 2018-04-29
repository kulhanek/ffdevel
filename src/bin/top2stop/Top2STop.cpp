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
#include <Lapack.hpp>

//------------------------------------------------------------------------------

using namespace std;

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
    nnb_types = 0;
    dih_mode = 0;
    dih_samp_freq = 5;
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
    WriteNBList(sout);
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
    sout << "! Index Type Name ResID ResN     Charge   Type" << endl;

    // write atoms
    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        sout << right << setw(7) << i+1 << " ";
        sout << right << setw(4) << FindAtomTypeIdx(i) << " ";
        sout << left << setw(4) << p_atom->GetName() << " ";
        sout << right << setw(5) << p_atom->GetResidue()->GetIndex()+1 << " ";
        sout << left << setw(4) << p_atom->GetResidue()->GetName() << " ";
        sout << right << fixed << setw(10) << setprecision(6) << p_atom->GetStandardCharge() << " ! ";
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
    sout << "! Index TypeA TypeB Form            d0             K   TypeA TypeB" << endl;

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
    sout << "! Index AtomA AtomB Type   AtomA TypeA AtomB TypeB" << endl;

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
    sout << "! Index TypeA TypeB TypeC Form            a0             K   TypeA TypeB TypeC" << endl;

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
    sout << "! Index AtomA AtomB AtomC Type   AtomA TypeA AtomB TypeB AtomC TypeC" << endl;

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

    sout << "[dihedral_types]" << endl;
    sout << "! Index TypeA TypeB TypeC TypeD Form          scee          scnb   TypeA TypeB TypeC TypeD" << endl;

    std::map<int,CDihedralType>::iterator it = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator ie = DihedralTypes.end();

    dih_mode = 0;
    if( Options.GetOptDihedralMode() == "cos" ){
        dih_mode = 1;
    } else if( Options.GetOptDihedralMode() == "grbf" ){
        dih_mode = 2;
    }
    if( dih_mode == 0 ){
        RUNTIME_ERROR("unsupported dihedral mode - p1");
    }

    while( it != ie ){
        CDihedralType dtype = it->second;
        sout << right << setw(7) << dtype.idx << " ";
        sout << right << setw(5) << dtype.at1 << " ";
        sout << right << setw(5) << dtype.at2 << " ";
        sout << right << setw(5) << dtype.at3 << " ";
        sout << right << setw(5) << dtype.at4 << " ";
        sout << right << setw(4) << dih_mode << " ";
        sout << right << fixed << setw(13) << setprecision(6) << dtype.scee << " ";
        sout << right << fixed << setw(13) << setprecision(6) << dtype.scnb << " ! ";
        sout << left << setw(5) << AtomTypes[dtype.at1].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at2].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at3].name << " ";
        sout << left << setw(5) << AtomTypes[dtype.at4].name << endl;
        it++;
    }

    ndihedral_types = DihedralTypes.size();

    switch(dih_mode){
        case 1:
            WriteDihedralSeqCosMode(sout);
            break;
        case 2:
            TransformCosToGRBF();
            WriteDihedralSeqGRBFMode(sout);
            break;
        default:
            RUNTIME_ERROR("unsupported dihedral mode - p2");
    }
}

//------------------------------------------------------------------------------

void CTop2STop::WriteDihedralSeqCosMode(ostream& sout)
{
    sout << "[dihedral_seq_cos]" << endl;
    sout << "! Type pn            v0         phase" << endl;

    std::map<int,CDihedralType>::iterator it = DihedralTypes.begin();
    std::map<int,CDihedralType>::iterator ie = DihedralTypes.end();

    while( it != ie ){
        CDihedralType dtype = it->second;
        for(int i=0; i < dtype.GetSeriesSize(); i++ ){
            sout << right << setw(6) << dtype.idx << " ";
            sout << right << setw(2) << i+1 << " ";
            if( Options.GetOptZeroDihPhase() ) {
                if( fabs(dtype.phase[i] - M_PI) < 0.1 ) {
                    dtype.v0[i] = - dtype.v0[i];
                    dtype.phase[i] = dtype.phase[i] - M_PI;
                }
                sout << right << fixed << setw(13) << setprecision(6) << dtype.v0[i] << " ";
                sout << right << fixed << setw(13) << setprecision(6) << dtype.phase[i];
            } else {
                sout << right << fixed << setw(13) << setprecision(6) << dtype.v0[i] << " ";
                sout << right << fixed << setw(13) << setprecision(6) << dtype.phase[i];
            }
            sout << endl;
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
        SolveTransformation(type.idx);
        it++;
    }
}

//------------------------------------------------------------------------------

void CTop2STop::SolveTransformation(int type)
{
    vout << "   fitting dihedral type " << setw(4) << type << " ... final error = ";

    // transform
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
        for(int i=0; i < dtype.GetSeriesSize(); i++ ){
            sout << right << setw(6) << dtype.idx << " ";
            sout << right << setw(2) << i+1 << " ";
            sout << right << fixed << setw(13) << setprecision(6) << dtype.c[i] << " ";
            sout << right << fixed << setw(13) << setprecision(6) << dtype.p[i] << " ";
            sout << right << fixed << setw(13) << setprecision(6) << dtype.w2[i];
            sout << endl;
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
    sout << "! Index AtomA AtomB AtomC AtomD Type   AtomA TypeA AtomB TypeB AtomC TypeC AtomD TypeD" << endl;

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
    sout << "! Index TypeA TypeB TypeC TypeD            v0         phase   TypeA TypeB TypeC TypeD" << endl;

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
    std::list<CDihedral> UniqueDihedrals;

    // list unique dihedrals
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++) {
        CAmberDihedral* p_dihedral = Topology.DihedralList.GetDihedral(i);
        int ip = p_dihedral->GetIP();
        int jp = p_dihedral->GetJP();
        int kp = p_dihedral->GetKP();
        int lp = p_dihedral->GetLP();

        // only improper and terminal improper dihedrals are processed
        if( (p_dihedral->GetType() != 1)&&(p_dihedral->GetType() != -2) ) continue;

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

    sout << "[impropers]" << endl;
    sout << "! Index AtomA AtomB AtomC AtomD Type   AtomA TypeA AtomB TypeB AtomC TypeC AtomD TypeD" << endl;

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

    nimpropers = UniqueDihedrals.size();
}

//------------------------------------------------------------------------------

void CTop2STop::WriteNBTypes(ostream& sout)
{
    map<int,int> Types;

    for(int i=0; i < Topology.AtomList.GetNumberOfAtoms(); i++) {
        CAmberAtom* p_atom = Topology.AtomList.GetAtom(i);
        Types[p_atom->GetIAC()-1] = FindAtomTypeIdx(i);
    }

    sout << "[nb_types]" << endl;
    sout << "! Index TypA TypB           eps            r0         alpha ! TypeA TypeB" << endl;

    int idx = 1;
    for(int i=0; i < Topology.NonBondedList.GetNumberOfTypes(); i++){
        for(int j=i; j < Topology.NonBondedList.GetNumberOfTypes(); j++){

            // ij parameters
            int icoij = Topology.NonBondedList.GetICOIndex(Topology.NonBondedList.GetNumberOfTypes()*i + j);

            NBTypes[icoij] = idx;

            double aij,bij,epsij,rij;

            aij = Topology.NonBondedList.GetAParam(icoij);
            bij = Topology.NonBondedList.GetBParam(icoij);

            epsij = 0.0;
            rij = 0.0;
            if( aij != 0.0 ) {
                epsij = bij*bij / (4.0 * aij);
                rij = pow(2*aij/bij,1.0/6.0) * 0.5;
            }

            int it = Types[i];
            int jt = Types[j];

            sout << right << setw(7) << idx << " ";
            sout << left << setw(4) << it << " ";
            sout << left << setw(4) << jt << " ";
            sout << right << fixed << setw(13) << setprecision(7) << epsij << " ";
            sout << right << fixed << setw(13) << setprecision(7) << rij*2.0 << " ";
            sout << right << fixed << setw(13) << setprecision(7) << 0.0 << " ! ";
            sout << left << setw(5) << AtomTypes[it].name << " ";
            sout << left << setw(5) << AtomTypes[jt].name << endl;
            idx++;
        }
    }

    nnb_types = idx-1;
}

//------------------------------------------------------------------------------

void CTop2STop::WriteNBList(ostream& sout)
{
    sout << "[nb_list]" << endl;
    sout << "! Index AtomA AtomB  Type Dihed   AtomA TypeA AtomB TypeB" << endl;

    nb_size = 0;
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
                didx = dtype.idx;
            }

            int ico = Topology.NonBondedList.GetICOIndex(p_atom1,p_atom2);
            int nbidx = NBTypes[ico];

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
    sout << "nb_size           " << nb_size << endl;
    sout << "nb_types          " << nnb_types << endl;

//    integer,parameter               :: NB_MODE_LJ = 1   ! Lennard-Jones potential
//    integer,parameter               :: NB_MODE_BP = 2   ! Buckingham potential

    sout << "nb_mode           " << 1 << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



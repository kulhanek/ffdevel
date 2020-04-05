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
#include "TopFilter.hpp"
#include <vector>
#include <iomanip>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//------------------------------------------------------------------------------

MAIN_ENTRY(CTopFilter);

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CTopFilter::CTopFilter(void)
{
}

//------------------------------------------------------------------------------

CTopFilter::~CTopFilter(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CTopFilter::Init(int argc,char* argv[])
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
    vout << "# topfilter (FFDevel utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgInTopologyName() != "-") {
        vout << "# AMBER topology file (in)   : " << Options.GetArgInTopologyName() << endl;
    } else {
        vout << "# AMBER topology file (in)   : - (standard input)" << endl;
    }
        vout << "# AMBER coordinate file (in) : " << Options.GetArgInCoordinatesName() << endl;

    if(Options.GetArgOutTopologyName() != "-") {
        vout << "# AMBER topology file (out)  : " << Options.GetArgOutTopologyName() << endl;
    } else {
        vout << "# AMBER topology file (out)  : - (standard output)" << endl;
    }
        vout << "# Applied filters            : " << Options.GetArgFilter() << endl;

    return( SO_CONTINUE );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CTopFilter::Run(void)
{
    // load topology
    if( LoadTopologyAndCoordinates() == false ) return(false);

    // filter
    if( FilterTopology() == false ) return(false);

    // save topology
    if( SaveTopology() == false ) return(false);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CTopFilter::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# topfilter terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CTopFilter::LoadTopologyAndCoordinates(void)
{
    if( Topology.Load(Options.GetArgInTopologyName(),true) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER topology: " << Options.GetArgInTopologyName() << "</red>" << endl;
        return(false);
    }
    Restart.AssignTopology(&Topology);
    if( Restart.Load(Options.GetArgInCoordinatesName()) == false ) {
        vout << "<red>>>> ERROR: Unable to load AMBER coordinates: " << Options.GetArgInCoordinatesName() << "</red>" << endl;
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CTopFilter::SaveTopology(void)
{
    if( Topology.Save(Options.GetArgOutTopologyName(),AMBER_VERSION_7, true) == false ) {
        vout << "<red>>>> ERROR: Unable to save AMBER topology: " << Options.GetArgOutTopologyName() << "</red>" << endl;
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CTopFilter::FilterTopology(void)
{
    string sfilters(Options.GetArgFilter());
    vector<string> filters;
    split(filters,sfilters,is_any_of("/"));

    vector<string>::iterator    it = filters.begin();
    vector<string>::iterator    ie = filters.end();

    while( it != ie ){
        string filter = *it;

        bool result = false;
        if( filter == "langles" ){
            result = FilterLAngles();
        } else if( filter == "zangles" ){
            result = FilterZAngles();
        } else if( filter == "ldihedrals" ){
            result =FilterLDihedrals();
        } else if( filter == "zdihedrals" ){
            result =FilterZDihedrals();
        }  else {
            vout << "<red>>>> ERROR: Unsupported filter: " << filter << "</red>" << endl;
        }
        if( result == false ) return(false);

        it++;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

class CAngle {
public:
    int         tidx; // type index
    int         ai,aj,ak;
};

//------------------------------------------------------------------------------

class CAngleType {
public:
    int         nidx; // new type index
    double      a0;
    double      k;
};

//------------------------------------------------------------------------------

bool CTopFilter::FilterLAngles(void)
{
    vout << endl;
    vout << ">>> FILTER: langles - remove straight angles from topology" << endl;

    vector<CAngle>  angles;
    vout << endl;
    vout << "    Number of angles          : " << Topology.AngleList.GetNumberOfAngles() << endl;
    int nsa = 0;
    for(int i=0; i < Topology.AngleList.GetNumberOfAngles(); i++){
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(i);
        CAngle       angle;
        angle.ai = p_angle->GetIT();
        angle.aj = p_angle->GetJT();
        angle.ak = p_angle->GetKT();
        angle.tidx = p_angle->GetICT();

        // get angle value
        CPoint pi = Restart.GetPosition(p_angle->GetIT());
        CPoint pj = Restart.GetPosition(p_angle->GetJT());
        CPoint pk = Restart.GetPosition(p_angle->GetKT());

        double value = Angle(pi-pj,pk-pj);

        if( value > (180.0-10.0)*M_PI/180.0 ){
            angle.tidx = -1;
            vout << "       ... found straight angle id " << i+1 << " value "
                 << fixed << setprecision(1) << value*180.0/M_PI
                 << " atoms " << p_angle->GetIT()+1 << "-" << p_angle->GetJT()+1
                 << "-" << p_angle->GetKT()+1 << " type " << p_angle->GetICT()+1 << endl;
            nsa++;
        }
        angles.push_back(angle);
    }
    vout << "    Number of straight angles : " << nsa << endl;

    vector<CAngleType>  atypes;
    vout << endl;
    vout << "    Number of angle types          : " << Topology.AngleList.GetNumberOfAngleTypes() << endl;

    // find all angles types that are involved only in the straight angles
    int nsat = 0;
    for(int i=0; i < Topology.AngleList.GetNumberOfAngleTypes(); i++){
        CAmberAngleType* p_atype = Topology.AngleList.GetAngleType(i);
        CAngleType atype;
        atype.a0 = p_atype->GetTEQ();
        atype.k = p_atype->GetTK();
        atype.nidx = i;
        bool found = false;
        for(unsigned int j=0; j < angles.size(); j++){
            if( angles[j].tidx == i ){
                found = true;
                break;
            }
        }
        if( found == false ){
            atype.nidx = -1; // disable type
            vout << "       ... found angle type id " << i+1 << " participated only in straight angles" << endl;
            nsat++;
        }
        atypes.push_back(atype);
    }
    vout << "    Number of straight angle types : " << nsat << endl;

    // renumber angle types
    int nat = 0;
    for(unsigned int i=0; i < atypes.size(); i++ ){
        if( atypes[i].nidx >= 0 ){
            atypes[i].nidx = nat;
            nat++;
        }
    }

    vout << endl;
    vout << "    Rebuilding topology ..." << endl;
    vout << "    Removing " << nsa << " angles ..." << endl;

    // count angles with and without hydrogens
    int nah = 0;
    int na = 0;
    for(unsigned int i=0; i < angles.size(); i++){
        if( i < Topology.AngleList.GetNumberOfAnglesWithoutHydrogen() ){
            if( angles[i].tidx >= 0 ) na++;
        } else {
            if( angles[i].tidx >= 0 ) nah++;
        }
    }
    vout << "       and keeping angles without hydrogens " << na << endl;
    vout << "       and keeping angles with hydrogens " << nah << endl;
    vout << "    Removing " << nsat << " types ..." << endl;

//    int NTHETH;   //NTHETH : number of angles containing hydrogen
//    int MTHETA;   //MTHETA : number of angles not containing hydrogen
//    int NUMANG;   //NUMANG : number of unique angle types
//    int NTHETA;   //NTHETA : MTHETA + number of constraint angles
//    int NGPER;    //NGPER  : number of angles to be perturbed
//    int MGPER;    //MGPER  : number of angles with atoms completely in perturbed group

//   void CAmberAngleList::InitFields(int iNTHETH, int iMTHETA,
//        int iNUMANG, int iNTHETA, int iNGPER, int iMGPER)

    Topology.AngleList.InitFields(nah,na,nat,na,0,0);

    // set angles
    int idx = 0;
    for(unsigned int i=0; i < angles.size(); i++){
        if( angles[i].tidx < 0 ) continue;
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(idx);
        p_angle->SetIT(angles[i].ai);
        p_angle->SetJT(angles[i].aj);
        p_angle->SetKT(angles[i].ak);
        p_angle->SetICT(atypes[angles[i].tidx].nidx);
        idx++;
    }

    // set types
    idx = 0;
    for(unsigned int i=0; i < atypes.size(); i++){
        if( atypes[i].nidx < 0 ) continue;
        CAmberAngleType* p_atype = Topology.AngleList.GetAngleType(idx);
        p_atype->SetTEQ(atypes[i].a0);
        p_atype->SetTK(atypes[i].k);
        idx++;
    }

    vout << endl;
    vout << "<red><b>    NOTE: Exclussion list is not modified!" << endl;
    vout <<         "          Atoms from removed angles should be excluded from NB interactions!" << endl;
    vout <<         "          To achieve that please setup corectly control file for sander/pmemd!</b></red>" << endl;
    return(true);
}

//------------------------------------------------------------------------------

bool CTopFilter::FilterZAngles(void)
{
    vout << endl;
    vout << ">>> FILTER: zangles - remove zero angles from topology" << endl;

    vector<CAngle>  angles;
    vout << endl;
    vout << "    Number of angles          : " << Topology.AngleList.GetNumberOfAngles() << endl;
    int nsa = 0;
    for(int i=0; i < Topology.AngleList.GetNumberOfAngles(); i++){
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(i);
        CAngle       angle;
        angle.ai = p_angle->GetIT();
        angle.aj = p_angle->GetJT();
        angle.ak = p_angle->GetKT();
        angle.tidx = p_angle->GetICT();

        // get angle value
        CPoint pi = Restart.GetPosition(p_angle->GetIT());
        CPoint pj = Restart.GetPosition(p_angle->GetJT());
        CPoint pk = Restart.GetPosition(p_angle->GetKT());

        double value = Angle(pi-pj,pk-pj);

        if( value < 10.0*M_PI/180.0 ){
            angle.tidx = -1;
            vout << "       ... found zero angle id " << i+1 << " value "
                 << fixed << setprecision(1) << value*180.0/M_PI
                 << " atoms " << p_angle->GetIT()+1 << "-" << p_angle->GetJT()+1
                 << "-" << p_angle->GetKT()+1 << " type " << p_angle->GetICT()+1 << endl;
            nsa++;
        }
        angles.push_back(angle);
    }
    vout << "    Number of zero angles : " << nsa << endl;

    vector<CAngleType>  atypes;
    vout << endl;
    vout << "    Number of angle types          : " << Topology.AngleList.GetNumberOfAngleTypes() << endl;

    // find all angles types that are involved only in the zero angles
    int nsat = 0;
    for(int i=0; i < Topology.AngleList.GetNumberOfAngleTypes(); i++){
        CAmberAngleType* p_atype = Topology.AngleList.GetAngleType(i);
        CAngleType atype;
        atype.a0 = p_atype->GetTEQ();
        atype.k = p_atype->GetTK();
        atype.nidx = i;
        bool found = false;
        for(unsigned int j=0; j < angles.size(); j++){
            if( angles[j].tidx == i ){
                found = true;
                break;
            }
        }
        if( found == false ){
            atype.nidx = -1; // disable type
            vout << "       ... found angle type id " << i+1 << " participated only in zero angles" << endl;
            nsat++;
        }
        atypes.push_back(atype);
    }
    vout << "    Number of zero angle types : " << nsat << endl;

    // renumber angle types
    int nat = 0;
    for(unsigned int i=0; i < atypes.size(); i++ ){
        if( atypes[i].nidx >= 0 ){
            atypes[i].nidx = nat;
            nat++;
        }
    }

    vout << endl;
    vout << "    Rebuilding topology ..." << endl;
    vout << "    Removing " << nsa << " angles ..." << endl;

    // count angles with and without hydrogens
    int nah = 0;
    int na = 0;
    for(unsigned int i=0; i < angles.size(); i++){
        if( i < Topology.AngleList.GetNumberOfAnglesWithoutHydrogen() ){
            if( angles[i].tidx >= 0 ) na++;
        } else {
            if( angles[i].tidx >= 0 ) nah++;
        }
    }
    vout << "       and keeping angles without hydrogens " << na << endl;
    vout << "       and keeping angles with hydrogens " << nah << endl;
    vout << "    Removing " << nsat << " types ..." << endl;

//    int NTHETH;   //NTHETH : number of angles containing hydrogen
//    int MTHETA;   //MTHETA : number of angles not containing hydrogen
//    int NUMANG;   //NUMANG : number of unique angle types
//    int NTHETA;   //NTHETA : MTHETA + number of constraint angles
//    int NGPER;    //NGPER  : number of angles to be perturbed
//    int MGPER;    //MGPER  : number of angles with atoms completely in perturbed group

//   void CAmberAngleList::InitFields(int iNTHETH, int iMTHETA,
//        int iNUMANG, int iNTHETA, int iNGPER, int iMGPER)

    Topology.AngleList.InitFields(nah,na,nat,na,0,0);

    // set angles
    int idx = 0;
    for(unsigned int i=0; i < angles.size(); i++){
        if( angles[i].tidx < 0 ) continue;
        CAmberAngle* p_angle = Topology.AngleList.GetAngle(idx);
        p_angle->SetIT(angles[i].ai);
        p_angle->SetJT(angles[i].aj);
        p_angle->SetKT(angles[i].ak);
        p_angle->SetICT(atypes[angles[i].tidx].nidx);
        idx++;
    }

    // set types
    idx = 0;
    for(unsigned int i=0; i < atypes.size(); i++){
        if( atypes[i].nidx < 0 ) continue;
        CAmberAngleType* p_atype = Topology.AngleList.GetAngleType(idx);
        p_atype->SetTEQ(atypes[i].a0);
        p_atype->SetTK(atypes[i].k);
        idx++;
    }

    vout << endl;
    vout << "<red><b>    NOTE: Exclussion list is not modified!" << endl;
    vout <<         "          Atoms from removed angles should be excluded from NB interactions!" << endl;
    vout <<         "          To achieve that please setup corectly control file for sander/pmemd!</b></red>" << endl;
    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

class CDihedral {
public:
    int         tidx; // type index
    int         ai,aj,ak,al;
    int         type;
};

//------------------------------------------------------------------------------

class CDihedralType {
public:
    int         nidx; // new type index
    double      v0;
    double      gamma;
    int         pn;
    double      scee;
    double      scnb;
};

//------------------------------------------------------------------------------

bool CTopFilter::FilterLDihedrals(void)
{
    vout << endl;
    vout << ">>> FILTER: ldihedrals - remove straight angle dihedrals from topology" << endl;

    vector<CDihedral>  dihedrals;
    vout << endl;
    vout << "    Number of dihedrals          : " << Topology.DihedralList.GetNumberOfDihedrals() << endl;
    int nsa = 0;
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++){
        CAmberDihedral* p_dih = Topology.DihedralList.GetDihedral(i);
        CDihedral  dih;
        dih.ai = p_dih->GetIP();
        dih.aj = p_dih->GetJP();
        dih.ak = p_dih->GetKP();
        dih.al = p_dih->GetLP();
        dih.tidx = p_dih->GetICP();
        dih.type = p_dih->GetType();

        // get angle value
        CPoint pi = Restart.GetPosition(p_dih->GetIP());
        CPoint pj = Restart.GetPosition(p_dih->GetJP());
        CPoint pk = Restart.GetPosition(p_dih->GetKP());
        CPoint pl = Restart.GetPosition(p_dih->GetLP());

        double value1 = Angle(pi-pj,pk-pj);
        double value2 = Angle(pj-pk,pl-pk);

        if( (value1 > (180.0-10.0)*M_PI/180.0) || (value2 > (180.0-10.0)*M_PI/180.0) ){
            dih.tidx = -1;
            vout << "       ... found straight angle dihedral id " << i+1 << " values ";
            vout << fixed << setprecision(1) << value1*180.0/M_PI << ";";
            vout << fixed << setprecision(1) << value2*180.0/M_PI;

            vout  << " atoms " << p_dih->GetIP()+1 << "-" << p_dih->GetJP()+1
                 << "-" << p_dih->GetKP()+1 << "-" << p_dih->GetLP()+1 << " type " << p_dih->GetICP()+1 << endl;
            nsa++;
        }
        dihedrals.push_back(dih);
    }
    vout << "    Number of straight angle dihedrals : " << nsa << endl;

    vector<CDihedralType>  dtypes;
    vout << endl;
    vout << "    Number of dihedral types          : " << Topology.DihedralList.GetNumberOfDihedralTypes() << endl;

    // find all angles types that are involved only in the straight angles
    int nsat = 0;
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedralTypes(); i++){
        CAmberDihedralType* p_dtype = Topology.DihedralList.GetDihedralType(i);
        CDihedralType dtype;
        dtype.v0 = p_dtype->GetPK();
        dtype.gamma = p_dtype->GetPHASE();
        dtype.scee = p_dtype->GetSCEE();
        dtype.scnb = p_dtype->GetSCNB();
        dtype.pn = p_dtype->GetPN();
        dtype.nidx = i;
        bool found = false;
        for(unsigned int j=0; j < dihedrals.size(); j++){
            if( dihedrals[j].tidx == i ){
                found = true;
                break;
            }
        }
        if( found == false ){
            dtype.nidx = -1; // disable type
            vout << "       ... found dihedral type id " << i+1 << " participated only in straight angle dihedrals" << endl;
            nsat++;
        }
        dtypes.push_back(dtype);
    }
    vout << "    Number of straight angle dihedral types : " << nsat << endl;

    // renumber dihedral types
    int nat = 0;
    for(unsigned int i=0; i < dtypes.size(); i++ ){
        if( dtypes[i].nidx >= 0 ){
            dtypes[i].nidx = nat;
            nat++;
        }
    }

    vout << endl;
    vout << "    Rebuilding topology ..." << endl;
    vout << "    Removing " << nsa << " dihedrals ..." << endl;

    // count dihedrals with and without hydrogens
    int nah = 0;
    int na = 0;
    for(unsigned int i=0; i < dihedrals.size(); i++){
        if( i < Topology.DihedralList.GetNumberOfDihedralsWithoutHydrogen() ){
            if( dihedrals[i].tidx >= 0 ) na++;
        } else {
            if( dihedrals[i].tidx >= 0 ) nah++;
        }
    }
    vout << "       and keeping dihedrals without hydrogens " << na << endl;
    vout << "       and keeping dihedrals with hydrogens " << nah << endl;
    vout << "    Removing " << nsat << " types ..." << endl;

//    int NPHIH;    //NPHIH  : number of dihedrals containing hydrogen
//    int MPHIA;    //MPHIA  : number of dihedrals not containing hydrogen
//    int NPTRA;    //NPTRA  : number of unique dihedral types
//    int NPHIA;    //NPHIA  : MPHIA + number of constraint dihedrals
//    int NDPER;    //NDPER  : number of dihedrals to be perturbed
//    int MDPER;    //MDPER  : number of dihedrals with atoms completely in perturbed groups

//    /// init all fields - old data are destroyed
//    void InitFields(int iNPHIH, int iMPHIA, int iNPTRA,
//                    int iNPHIA, int iNDPER, int iMDPER);

    Topology.DihedralList.InitFields(nah,na,nat,na,0,0);

    // set dihedrals
    int idx = 0;
    for(unsigned int i=0; i < dihedrals.size(); i++){
        if( dihedrals[i].tidx < 0 ) continue;
        CAmberDihedral* p_dih = Topology.DihedralList.GetDihedral(idx);
        p_dih->SetIP(dihedrals[i].ai);
        p_dih->SetJP(dihedrals[i].aj);
        p_dih->SetKP(dihedrals[i].ak);
        p_dih->SetLP(dihedrals[i].al);
        p_dih->SetICP(dtypes[dihedrals[i].tidx].nidx);
        p_dih->SetType(dihedrals[i].type);
        idx++;
    }

    // set types
    idx = 0;
    for(unsigned int i=0; i < dtypes.size(); i++){
        if( dtypes[i].nidx < 0 ) continue;
        CAmberDihedralType* p_dtype = Topology.DihedralList.GetDihedralType(idx);
        p_dtype->SetPK(dtypes[i].v0);
        p_dtype->SetPHASE(dtypes[i].gamma);
        p_dtype->SetPN(dtypes[i].pn);
        p_dtype->SetSCEE(dtypes[i].scee);
        p_dtype->SetSCNB(dtypes[i].scnb);
        idx++;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CTopFilter::FilterZDihedrals(void)
{
    vout << endl;
    vout << ">>> FILTER: zdihedrals - remove zero angle dihedrals from topology" << endl;

    vector<CDihedral>  dihedrals;
    vout << endl;
    vout << "    Number of dihedrals          : " << Topology.DihedralList.GetNumberOfDihedrals() << endl;
    int nsa = 0;
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedrals(); i++){
        CAmberDihedral* p_dih = Topology.DihedralList.GetDihedral(i);
        CDihedral  dih;
        dih.ai = p_dih->GetIP();
        dih.aj = p_dih->GetJP();
        dih.ak = p_dih->GetKP();
        dih.al = p_dih->GetLP();
        dih.tidx = p_dih->GetICP();
        dih.type = p_dih->GetType();

        // DEBUG
        //   vout << format("%10d-%10d-%10d-%10d %d")%dih.ai%dih.aj%dih.ak%dih.al%dih.tidx << endl;

        // get angle value
        CPoint pi = Restart.GetPosition(p_dih->GetIP());
        CPoint pj = Restart.GetPosition(p_dih->GetJP());
        CPoint pk = Restart.GetPosition(p_dih->GetKP());
        CPoint pl = Restart.GetPosition(p_dih->GetLP());

        double value1 = Angle(pi-pj,pk-pj);
        double value2 = Angle(pj-pk,pl-pk);

        if( (value1 < 10.0*M_PI/180.0) || (value2 < 10.0*M_PI/180.0) ){
            dih.tidx = -1;
            vout << "       ... found zero angle dihedral id " << i+1 << " values ";
            vout << fixed << setprecision(1) << value1*180.0/M_PI << ";";
            vout << fixed << setprecision(1) << value2*180.0/M_PI;

            vout  << " atoms " << p_dih->GetIP()+1 << "-" << p_dih->GetJP()+1
                 << "-" << p_dih->GetKP()+1 << "-" << p_dih->GetLP()+1 << " type " << p_dih->GetICP()+1 << endl;
            nsa++;
        }
        dihedrals.push_back(dih);
    }
    vout << "    Number of zero angle dihedrals : " << nsa << endl;

    vector<CDihedralType>  dtypes;
    vout << endl;
    vout << "    Number of dihedral types          : " << Topology.DihedralList.GetNumberOfDihedralTypes() << endl;

    // find all angles types that are involved only in the zero angles
    int nsat = 0;
    for(int i=0; i < Topology.DihedralList.GetNumberOfDihedralTypes(); i++){
        CAmberDihedralType* p_dtype = Topology.DihedralList.GetDihedralType(i);
        CDihedralType dtype;
        dtype.v0 = p_dtype->GetPK();
        dtype.gamma = p_dtype->GetPHASE();
        dtype.scee = p_dtype->GetSCEE();
        dtype.scnb = p_dtype->GetSCNB();
        dtype.pn = p_dtype->GetPN();
        dtype.nidx = i;
        bool found = false;
        for(unsigned int j=0; j < dihedrals.size(); j++){
            if( dihedrals[j].tidx == i ){
                found = true;
                break;
            }
        }
        if( found == false ){
            dtype.nidx = -1; // disable type
            vout << "       ... found dihedral type id " << i+1 << " participated only in zero angle dihedrals" << endl;
            nsat++;
        }
        dtypes.push_back(dtype);
    }
    vout << "    Number of zero angle dihedral types : " << nsat << endl;

    // renumber dihedral types
    int nat = 0;
    for(unsigned int i=0; i < dtypes.size(); i++ ){
        if( dtypes[i].nidx >= 0 ){
            dtypes[i].nidx = nat;
            nat++;
        }
    }

    vout << endl;
    vout << "    Rebuilding topology ..." << endl;
    vout << "    Removing " << nsa << " dihedrals ..." << endl;

    // count dihedrals with and without hydrogens
    int nah = 0;
    int na = 0;
    for(unsigned int i=0; i < dihedrals.size(); i++){
        if( i < Topology.DihedralList.GetNumberOfDihedralsWithoutHydrogen() ){
            if( dihedrals[i].tidx >= 0 ) na++;
        } else {
            if( dihedrals[i].tidx >= 0 ) nah++;
        }
    }
    vout << "       and keeping dihedrals without hydrogens " << na << endl;
    vout << "       and keeping dihedrals with hydrogens " << nah << endl;
    vout << "    Removing " << nsat << " types ..." << endl;

//    int NPHIH;    //NPHIH  : number of dihedrals containing hydrogen
//    int MPHIA;    //MPHIA  : number of dihedrals not containing hydrogen
//    int NPTRA;    //NPTRA  : number of unique dihedral types
//    int NPHIA;    //NPHIA  : MPHIA + number of constraint dihedrals
//    int NDPER;    //NDPER  : number of dihedrals to be perturbed
//    int MDPER;    //MDPER  : number of dihedrals with atoms completely in perturbed groups

//    /// init all fields - old data are destroyed
//    void InitFields(int iNPHIH, int iMPHIA, int iNPTRA,
//                    int iNPHIA, int iNDPER, int iMDPER);

    Topology.DihedralList.InitFields(nah,na,nat,na,0,0);

    // set dihedrals
    int idx = 0;
    for(unsigned int i=0; i < dihedrals.size(); i++){
        if( dihedrals[i].tidx < 0 ) continue;
        CAmberDihedral* p_dih = Topology.DihedralList.GetDihedral(idx);
        p_dih->SetIP(dihedrals[i].ai);
        p_dih->SetJP(dihedrals[i].aj);
        p_dih->SetKP(dihedrals[i].ak);
        p_dih->SetLP(dihedrals[i].al);
        p_dih->SetICP(dtypes[dihedrals[i].tidx].nidx);
        p_dih->SetType(dihedrals[i].type);
        idx++;
    }

    // set types
    idx = 0;
    for(unsigned int i=0; i < dtypes.size(); i++){
        if( dtypes[i].nidx < 0 ) continue;
        CAmberDihedralType* p_dtype = Topology.DihedralList.GetDihedralType(idx);
        p_dtype->SetPK(dtypes[i].v0);
        p_dtype->SetPHASE(dtypes[i].gamma);
        p_dtype->SetPN(dtypes[i].pn);
        p_dtype->SetSCEE(dtypes[i].scee);
        p_dtype->SetSCNB(dtypes[i].scnb);
        idx++;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



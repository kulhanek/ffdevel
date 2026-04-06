#ifndef Top2STopOptionsH
#define Top2STopOptionsH
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

#include <SimpleOptions.hpp>

//------------------------------------------------------------------------------

class CTop2STopOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CTop2STopOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "top2stop"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "Convert the AMBER topology to the simplified topology used by the ffoptimize program and other ffdevel utilities."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    "1.0"
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,TopologyName)
    CSO_ARG(CSmallString,STopologyName)
    // options ------------------------------
    CSO_OPT(CSmallString,CrdName)
    CSO_OPT(int,DihedralSeriesSize)
    CSO_OPT(CSmallString,DihedralMode)
    CSO_OPT(int,NDihSamples)
    CSO_OPT(double,GWidthFactor)
    CSO_OPT(CSmallString,DihedralTypes)
    CSO_OPT(bool,ZeroDihPhase)
    CSO_OPT(CSmallString,Transform)
    CSO_OPT(bool,RebuildNBList)
    CSO_OPT(double,DihCOffset)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, TopologyName, NULL, true, "TOPOLOGY",
                "AMBER topology file name. If the name is '-', the AMBER topology is read from the standard input.")
    // -------------------------------------------
        CSO_MAP_ARG(CSmallString, STopologyName, NULL, true, "STOPOLOGY",
                "Simplified topology file name. If the name is '-', the simplified topology is written to the standard output.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, CrdName, NULL, false, 'c', "coords", "NAME",
                "Input file with coordinates.")
    // -------------------------------------------
        CSO_MAP_OPT(int, DihedralSeriesSize, 4, false, 'd', "dihsize", "SIZE",
                "Size of the dihedral cosine series.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, DihedralMode, "cos", false, 'm', "dihmode", "MODE",
                "Dihedral representation mode: cos - cosine series; grbf - Gaussian radial basis functions; cbs - cubic B-spline.")
    // -------------------------------------------
        CSO_MAP_OPT(int, NDihSamples, 180, false, 's', "ndihsamples", NULL,
                "Number of samples for the 2pi rotation.")
    // -------------------------------------------
        CSO_MAP_OPT(double, GWidthFactor, 1.0, false, 'w', "wfactor", NULL,
                "Gaussian width modulation factor for GRBF.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, DihedralTypes, NULL, false, 'f', "dihfilters", "NAME",
                "File with atom types defining dihedral angles for GRBF transformation. Multiple filters can be specified, "
                "each on a separate line, either as two atom types for the central bond or four atom types for the exact dihedral type.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, ZeroDihPhase, false, false, 'z', "zerophase", NULL,
                "Transform all dihedral phases to zero, if applicable.")
    // -------------------------------------------
        CSO_MAP_OPT(CSmallString, Transform, "none", false, 't', "transform", "MODE",
                "Transform atom type names by letter capitalization. Allowed modes: none, first, second, both.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, RebuildNBList, false, false, 'r', "rebuild", NULL,
                "Rebuild the non-bonded list from scratch.")
    // -------------------------------------------
        CSO_MAP_OPT(double, DihCOffset, 0.0, false, 'o', "offset", "NUM",
                "Offset applied to dih_c.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Verbose, false, false, 'v', "verbose", NULL,
                "Increase output verbosity.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Version, false, false, '\0', "version", NULL,
                "Output version information and exit.")
    // -------------------------------------------
        CSO_MAP_OPT(bool, Help, false, false, 'h', "help", NULL,
                "Display this help and exit.")
    // -------------------------------------------
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
    virtual int CheckArguments(void);
};

//------------------------------------------------------------------------------

#endif

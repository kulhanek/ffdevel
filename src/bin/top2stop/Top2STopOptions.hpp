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
    CSO_OPT(int,DihedralSeriesSize)
    CSO_OPT(CSmallString,DihedralMode)
    CSO_OPT(CSmallString,DihedralTypes)
    CSO_OPT(bool,ZeroDihPhase)
    CSO_OPT(CSmallString,Transform)
    CSO_OPT(bool,RebuildNBList)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                TopologyName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "TOPOLOGY",                           /* parametr name */
                "AMBER topology name. If the name is '-' then the AMBER topology is read from the standard input.")   /* argument description */
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                STopologyName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "STOPOLOGY",                           /* parametr name */
                "Simplified topology name. If the name is '-' then the simplified topology is written to the standard output.")   /* argument description */
// description of options -----------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                DihedralSeriesSize,                        /* option name */
                4,                          /* default value */
                false,                          /* is option mandatory */
                'd',                           /* short option name */
                "dihsize",                      /* long option name */
                "SIZE",                           /* parametr name */
                "dihedral series size")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                DihedralMode,                        /* option name */
                "cos",                          /* default value */
                false,                          /* is option mandatory */
                'm',                           /* short option name */
                "dihmode",                      /* long option name */
                NULL,                           /* parametr name */
                "dihedral mode: cos - cosine series, grbf - gaussian radial basis functions")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                DihedralTypes,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                'f',                           /* short option name */
                "dihfilters",                      /* long option name */
                NULL,                           /* parametr name */
                "name of file with atoms types defining dihedrals for grbf transformation, multiple filters can be specified each on a line, "
                "either two atom types for central bond or four atom types for exact dihedaral type")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                ZeroDihPhase,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'z',                           /* short option name */
                "zerophase",                      /* long option name */
                NULL,                           /* parametr name */
                "transform all dihedral phases to zero if applicable")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(CSmallString,                           /* option type */
                Transform,                        /* option name */
                "none",                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "transform",                      /* long option name */
                "MODE",                           /* parametr name */
                "transform atom type names by letter capitalization. Allowed modes are: none, first, second, both.")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                RebuildNBList,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'r',                           /* short option name */
                "rebuild",                      /* long option name */
                NULL,                           /* parametr name */
                "rebuild NB list from scratch")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Verbose,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'v',                           /* short option name */
                "verbose",                      /* long option name */
                NULL,                           /* parametr name */
                "increase output verbosity")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Version,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "version",                      /* long option name */
                NULL,                           /* parametr name */
                "output version information and exit")   /* option description */
    //----------------------------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                Help,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'h',                           /* short option name */
                "help",                      /* long option name */
                NULL,                           /* parametr name */
                "display this help and exit")   /* option description */
    CSO_MAP_END

// final operation with options ------------------------------------------------
private:
    virtual int CheckOptions(void);
    virtual int FinalizeOptions(void);
    virtual int CheckArguments(void);
};

//------------------------------------------------------------------------------

#endif

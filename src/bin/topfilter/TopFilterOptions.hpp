#ifndef TopFilterOptionsH
#define TopFilterOptionsH
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

class CTopFilterOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CTopFilterOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "topfilter"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "Manipulate with some FF interactions present in AMBER topology."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    "1.0"
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,InTopologyName)
    CSO_ARG(CSmallString,InCoordinatesName)
    CSO_ARG(CSmallString,Filter)
    CSO_ARG(CSmallString,OutTopologyName)
    // options ------------------------------
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                InTopologyName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "INTOP",                           /* parametr name */
                "Input AMBER topology name. If the name is '-' then the AMBER topology is read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                InCoordinatesName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "INCRD",                           /* parametr name */
                "Input AMBER coordinate file name.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Filter,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "FILTERS",                           /* parametr name */
                "Applied filters (slash separated): langles - remove all straight angles, znagles - remove all zero angles, ldihedrals - remove all straight angle dihedrals, zdihedrals - remove all zero angle dihedrals")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                OutTopologyName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "OUTPUT",                           /* parametr name */
                "Output AMBER topology name. If the name is '-' then the AMBER topology is written to the standard output.")   /* argument description */
// description of options -----------------------------------------------------
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

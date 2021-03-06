#ifndef GenRotorsOptionsH
#define GenRotorsOptionsH
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

class CGenRotorsOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CGenRotorsOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "top2stop"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "Convert the AMBER topology to the simplified topology used by the ffoptimize program."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    "1.0"
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,TopologyName)
    CSO_ARG(CSmallString,CoordinateName)
    CSO_ARG(CSmallString,RotorName)
    // options ------------------------------
    CSO_OPT(bool,IncludeTerminals)
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
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                CoordinateName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "COORDINATES",                           /* parametr name */
                "AMBER coordinates name. If the name is '-' then the AMBER coordinates are read from the standard input.")   /* argument description */
    //----------------------------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                RotorName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "ROTORS",                           /* parametr name */
                "Output file name contaning freely rotatable bonds. If the name is '-' then the file is written to the standard output.")   /* argument description */

// description of options -----------------------------------------------------
    CSO_MAP_OPT(bool,                           /* option type */
                IncludeTerminals,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                't',                           /* short option name */
                "includeterminals",                      /* long option name */
                NULL,                           /* parametr name */
                "include terminal rotatable bonds")   /* option description */
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

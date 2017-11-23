#ifndef VDWGenProbeOptionsH
#define VDWGenProbeOptionsH
// =============================================================================
// This file is part of FFDevel.
//    Copyright (C) 2017 Petr Kulhanek, kulhanek@chemi.muni.cz
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

class CVDWGenProbeOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CVDWGenProbeOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "vdwprobefilter"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "Create multiple-xyz file containing the structure and probe generated on a rectangulated grid around the structure."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    "1.0"
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,StructureName)
    CSO_ARG(CSmallString,OutputName)
    // options ------------------------------
    CSO_OPT(double,RMin)
    CSO_OPT(double,RMax)
    CSO_OPT(double,Spacing)
    CSO_OPT(int,NumOfProbes)
    CSO_OPT(CSmallString,ProbeSymbol)
    CSO_OPT(CSmallString,SelectedAtoms)
    CSO_OPT(bool,Append)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                StructureName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "STRUCTURE",                           /* parametr name */
                "Structure in the xyz format. If the name is '-' then the structure is read from the standard input.")   /* argument description */
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                OutputName,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "OUTPUT",                           /* parametr name */
                "Multiple-xyz file containing the structure and probe. If the name is '-' then the structure is written to the standard output.")   /* argument description */
// description of options -----------------------------------------------------

    CSO_MAP_OPT(double,                           /* option type */
                RMin,                        /* option name */
                2.4,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "min",                      /* long option name */
                "NUM",                           /* parametr name */
                "minimal distance from any atom")   /* option description */
    CSO_MAP_OPT(double,                           /* option type */
                RMax,                        /* option name */
                5.0,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "max",                      /* long option name */
                "NUM",                           /* parametr name */
                "maximal distance from any atom")   /* option description */
    CSO_MAP_OPT(double,                           /* option type */
                Spacing,                        /* option name */
                0.1,                          /* default value */
                false,                          /* is option mandatory */
                '\0',                           /* short option name */
                "spacing",                      /* long option name */
                "NUM",                           /* parametr name */
                "grid point spacing")   /* option description */
    CSO_MAP_OPT(int,                           /* option type */
                NumOfProbes,                        /* option name */
                0,                          /* default value */
                false,                          /* is option mandatory */
                'n',                           /* short option name */
                "number",                      /* long option name */
                "NUM",                           /* parametr name */
                "number of generated probes, which are selected randomly; zero means all possible probes")   /* option description */
    CSO_MAP_OPT(CSmallString,                           /* option type */
                ProbeSymbol,                        /* option name */
                "Ar",                          /* default value */
                false,                          /* is option mandatory */
                'p',                           /* short option name */
                "probe",                      /* long option name */
                "NAME",                           /* parametr name */
                "probe symbol")   /* option description */
    CSO_MAP_OPT(CSmallString,                           /* option type */
                SelectedAtoms,                        /* option name */
                NULL,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "selected",                      /* long option name */
                "LIST",                           /* parametr name */
                "comma separated list of atom indexes (counted from 1)")   /* option description */
    CSO_MAP_OPT(bool,                           /* option type */
                Append,                        /* option name */
                false,                          /* default value */
                false,                          /* is option mandatory */
                'a',                           /* short option name */
                "append",                      /* long option name */
                NULL,                           /* parametr name */
                "append structures to the output file if it exists")   /* option description */
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

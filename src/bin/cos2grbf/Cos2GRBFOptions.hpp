#ifndef Cos2GRBFOptionsH
#define Cos2GRBFOptionsH
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

class CCos2GRBFOptions : public CSimpleOptions {
public:
    // constructor - tune option setup
    CCos2GRBFOptions(void);

// program name and description -----------------------------------------------
    CSO_PROG_NAME_BEGIN
    "cos2grbf"
    CSO_PROG_NAME_END

    CSO_PROG_DESC_BEGIN
    "Convert the cosinus dihedral series to rational gaussian basis functions."
    CSO_PROG_DESC_END

    CSO_PROG_VERS_BEGIN
    "1.0"
    CSO_PROG_VERS_END

// list of all options and arguments ------------------------------------------
    CSO_LIST_BEGIN
    // arguments ----------------------------
    CSO_ARG(CSmallString,Input)
    CSO_ARG(CSmallString,Output)
    // options ------------------------------
    CSO_OPT(int,DihedralSeriesSize)
    CSO_OPT(int,DihedralSamplingSize)
    CSO_OPT(bool,Help)
    CSO_OPT(bool,Version)
    CSO_OPT(bool,Verbose)
    CSO_LIST_END

    CSO_MAP_BEGIN
// description of arguments ---------------------------------------------------
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Input,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "COS",                           /* parametr name */
                "Name of file with cosinus series in the AMBER frcmod format.")  /* argument description */
    CSO_MAP_ARG(CSmallString,                   /* argument type */
                Output,                          /* argument name */
                NULL,                           /* default value */
                true,                           /* is argument mandatory */
                "OUT",                           /* parametr name */
                "Output file name.")   /* argument description */
// description of options -----------------------------------------------------
    CSO_MAP_OPT(int,                           /* option type */
                DihedralSeriesSize,                        /* option name */
                4,                          /* default value */
                false,                          /* is option mandatory */
                'd',                           /* short option name */
                "dihsize",                      /* long option name */
                "SIZE",                           /* parametr name */
                "dihedral series size")   /* option description */
    CSO_MAP_OPT(int,                           /* option type */
                DihedralSamplingSize,                        /* option name */
                10,                          /* default value */
                false,                          /* is option mandatory */
                's',                           /* short option name */
                "sampling",                      /* long option name */
                NULL,                           /* parametr name */
                "number of samples per segment")   /* option description */
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

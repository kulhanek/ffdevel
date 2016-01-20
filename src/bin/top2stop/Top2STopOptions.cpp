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

#include "Top2STopOptions.hpp"
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CTop2STopOptions::CTop2STopOptions(void)
{
    SetShowMiniUsage(true);
}

//------------------------------------------------------------------------------

int CTop2STopOptions::CheckOptions(void)
{
    if( GetOptDihedralMode() == "grbf" ){
        if( GetOptDihedralSeriesSize() < 6 ){
            if( IsError == false ) fprintf(stderr,"\n");
            fprintf(stderr,"%s: dihedral series size must have at least six members in grbf mode\n",
                    (const char*)GetProgramName() );
            IsError = true;
        }
    }
    if( GetOptDihedralMode() == "cos" ){
        if( GetOptDihedralSeriesSize() < 1 ){
            if( IsError == false ) fprintf(stderr,"\n");
            fprintf(stderr,"%s: dihedral series size must have at least one member in cos mode\n",
                    (const char*)GetProgramName() );
            IsError = true;
        }
    }

    if( IsError ) return(SO_OPTS_ERROR);
    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CTop2STopOptions::FinalizeOptions(void)
{
    bool ret_opt = false;

    if( GetOptHelp() == true ) {
        PrintUsage();
        ret_opt = true;
    }

    if( GetOptVersion() == true ) {
        PrintVersion();
        ret_opt = true;
    }

    if( ret_opt == true ) {
        printf("\n");
        return(SO_EXIT);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

int CTop2STopOptions::CheckArguments(void)
{
    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

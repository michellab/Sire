/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#ifndef SIRESEARCH_HELPER_FUNCS
#define SIRESEARCH_HELPER_FUNCS

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireSearch
{
    SIRESEARCH_EXPORT void install_search_parser();

    SIRESEARCH_EXPORT int get_min_protein_residues();
    SIRESEARCH_EXPORT void set_min_protein_residues(int nres);

    SIRESEARCH_EXPORT QSet<QString> get_protein_residue_names();
    SIRESEARCH_EXPORT void set_protein_residue_names(const QSet<QString> &names);

    SIRESEARCH_EXPORT void set_token(const QString &token, const QString &search);
    SIRESEARCH_EXPORT QString get_token(const QString &token);
    SIRESEARCH_EXPORT void delete_token(const QString &token);
    SIRESEARCH_EXPORT void delete_all_tokens();
}

SIRE_EXPOSE_FUNCTION( SireSearch::get_min_protein_residues )
SIRE_EXPOSE_FUNCTION( SireSearch::set_min_protein_residues )
SIRE_EXPOSE_FUNCTION( SireSearch::get_protein_residue_names )
SIRE_EXPOSE_FUNCTION( SireSearch::set_protein_residue_names )
SIRE_EXPOSE_FUNCTION( SireSearch::set_token )
SIRE_EXPOSE_FUNCTION( SireSearch::get_token )
SIRE_EXPOSE_FUNCTION( SireSearch::delete_token )
SIRE_EXPOSE_FUNCTION( SireSearch::delete_all_tokens )

SIRE_END_HEADER

#endif

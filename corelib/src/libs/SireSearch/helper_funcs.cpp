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

#include "helper_funcs.h"

#include "SireSearch/parser.h"

#include "SireMol/core.h"

#include <QReadWriteLock>

struct ProtRes
{
    int min_res = 5;
    QSet<QString> res_names = {"gly", "ala", "val", "leu", "ile", "ser",
                               "thr", "asp", "asn", "lys", "glu", "gln",
                               "arg", "his", "phe", "cys", "trp", "tyr",
                               "met", "pro", "ash", "glh", "cyx", "hid",
                               "hie", "hip"};
    QReadWriteLock lock;
};

Q_GLOBAL_STATIC(ProtRes, protres);


namespace SireSearch
{
    SIRESEARCH_EXPORT void install_search_parser()
    {
        SireSearch::parser::SearchParser::install();
    }

    SIRESEARCH_EXPORT int get_min_protein_residues()
    {
        auto p = protres();
        QReadLocker lkr(&(p->lock));
        return p->min_res;
    }

    SIRESEARCH_EXPORT void set_min_protein_residues(int nres)
    {
        if (nres < 0)
            nres = 0;

        auto p = protres();
        QWriteLocker lkr(&(p->lock));
        p->min_res = nres;
    }

    SIRESEARCH_EXPORT QSet<QString> get_protein_residue_names()
    {
        auto p = protres();
        QReadLocker lkr(&(p->lock));
        return p->res_names;
    }

    SIRESEARCH_EXPORT void set_protein_residue_names(const QSet<QString> &names)
    {
        // have to lower-case the names
        QSet<QString> r;

        for (const auto &name : names)
        {
            r.insert(name.toLower());
        }

        auto p = protres();
        QWriteLocker lkr(&(p->lock));
        p->res_names = r;
    }

    SIRESEARCH_EXPORT void set_token(const QString &token, const QString &search)
    {
        parser::SearchParser::globalParser().set_token(token, search);
    }

    SIRESEARCH_EXPORT QString get_token(const QString &token)
    {
        return parser::SearchParser::globalParser().get_token(token);
    }

    SIRESEARCH_EXPORT void delete_token(const QString &token)
    {
        parser::SearchParser::globalParser().delete_token(token);
    }

    SIRESEARCH_EXPORT void delete_all_tokens()
    {
        parser::SearchParser::globalParser().delete_all_tokens();
    }
}

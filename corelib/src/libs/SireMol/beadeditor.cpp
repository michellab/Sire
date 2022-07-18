/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2010  Christopher Woods
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

#include "beadeditor.h"
#include "mover.hpp"
#include "selector.hpp"
#include "atom.h"
#include "residue.h"
#include "chain.h"
#include "cutgroup.h"
#include "segment.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<BeadEditor> r_beadeditor;

QDataStream &operator<<(QDataStream &ds, const BeadEditor &beadeditor)
{
    writeHeader(ds, r_beadeditor, 1);

    ds << static_cast<const Bead&>(beadeditor);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, BeadEditor &beadeditor)
{
    VersionID v = readHeader(ds, r_beadeditor);

    if (v == 1)
    {
        ds >> static_cast<Bead&>(beadeditor);
    }
    else
        throw version_error(v, "1", r_beadeditor, CODELOC);

    return ds;
}

/** Null constructor */
BeadEditor::BeadEditor() : ConcreteProperty<BeadEditor,BeadEditorBase>()
{}

/** Constructor an editor for the passed bead */
BeadEditor::BeadEditor(const Bead &bead)
           : ConcreteProperty<BeadEditor,BeadEditorBase>(bead)
{}

/** Copy constructor */
BeadEditor::BeadEditor(const BeadEditor &other)
           : ConcreteProperty<BeadEditor,BeadEditorBase>(other)
{}

/** Destructor */
BeadEditor::~BeadEditor()
{}

/** Copy assignment operator */
BeadEditor& BeadEditor::operator=(const Bead &bead)
{
    Bead::operator=(bead);
    return *this;
}

/** Copy assignment operator */
BeadEditor& BeadEditor::operator=(const BeadEditor &other)
{
    Bead::operator=(other);
    return *this;
}

const char* BeadEditor::typeName()
{
    return QMetaType::typeName( qMetaTypeId<BeadEditor>() );
}

/** Return a string representation of this editor */
QString BeadEditor::toString() const
{
    return QObject::tr("Editor{ %1 }").arg(Bead::toString());
}

/** Commit the changes */
Bead BeadEditor::commit() const
{
    return Bead(*this);
}

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

#ifndef SIREMOL_BEADEDITOR_H
#define SIREMOL_BEADEDITOR_H

#include "bead.h"
#include "editor.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{
class BeadEditor;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::BeadEditor&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::BeadEditor&);

namespace SireMol
{

class BeadEditor;
typedef Editor<BeadEditor, Bead> BeadEditorBase;

/** This is an editor for a single bead in a molecule

    @author Christopher Woods
*/
class SIREMOL_EXPORT BeadEditor 
        : public SireBase::ConcreteProperty< BeadEditor,Editor<BeadEditor,Bead> >
{

friend QDataStream& ::operator<<(QDataStream&, const BeadEditor&);
friend QDataStream& ::operator>>(QDataStream&, BeadEditor&);

public:
    BeadEditor();
    
    BeadEditor(const Bead &bead);
    
    BeadEditor(const BeadEditor &other);
    
    ~BeadEditor();
    
    BeadEditor& operator=(const Bead &bead);
    BeadEditor& operator=(const BeadEditor &other);
    
    static const char* typeName();

    QString toString() const;
    
    Bead commit() const;
};

}

Q_DECLARE_METATYPE( SireMol::BeadEditor );

SIRE_EXPOSE_CLASS( SireMol::BeadEditor )

SIRE_EXPOSE_ALIAS( (SireMol::Editor<SireMol::BeadEditor, SireMol::Bead>),
                    SireMol::BeadEditorBase )

SIRE_END_HEADER

#endif

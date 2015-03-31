/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#include "combinespaces.h"

#include "space.h"
#include "combinedspace.h"

#include "SireBase/properties.h"

#include "SireStream/datastream.h"

using namespace SireVol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<CombineSpaces> r_combinespaces;

/** Serialise to a binary datastream */
QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, 
                                       const CombineSpaces &combinespaces)
{
    writeHeader(ds, r_combinespaces, 1);
    
    ds << static_cast<const CombineProperties&>(combinespaces);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, CombineSpaces &combinespaces)
{
    VersionID v = readHeader(ds, r_combinespaces);
    
    if (v == 1)
    {
        ds >> static_cast<CombineProperties&>(combinespaces);
    }
    else
        throw version_error( v, "1", r_combinespaces, CODELOC );
        
    return ds;
}

/** Constructor */
CombineSpaces::CombineSpaces() : ConcreteProperty<CombineSpaces,CombineProperties>()
{}

/** Construct to use just as single space, from the supplied source */
CombineSpaces::CombineSpaces(const PropertyName &source)
              : ConcreteProperty<CombineSpaces,CombineProperties>(source)
{}

/** Construct to combine together the two spaces specified by the 
    two supplied sources */
CombineSpaces::CombineSpaces(const PropertyName &source0, const PropertyName &source1)
              : ConcreteProperty<CombineSpaces,CombineProperties>(source0, source1)
{}

/** Construct to combine together the spaces from the passed sources */
CombineSpaces::CombineSpaces(const QList<PropertyName> &sources)
              : ConcreteProperty<CombineSpaces,CombineProperties>(sources)
{}

/** Construct to combine together the spaces from the passed sources */
CombineSpaces::CombineSpaces(const QVector<PropertyName> &sources)
              : ConcreteProperty<CombineSpaces,CombineProperties>(sources)
{}

/** Construct to combine together the spaces from the passed sources */
CombineSpaces::CombineSpaces(const QList<QString> &sources)
              : ConcreteProperty<CombineSpaces,CombineProperties>(sources)
{}

/** Construct to combine together the spaces from the passed sources */
CombineSpaces::CombineSpaces(const QVector<QString> &sources)
              : ConcreteProperty<CombineSpaces,CombineProperties>(sources)
{}

/** Copy constructor */
CombineSpaces::CombineSpaces(const CombineSpaces &other)
              : ConcreteProperty<CombineSpaces,CombineProperties>(other)
{}

/** Destructor */
CombineSpaces::~CombineSpaces()
{}

/** Copy assignment operator */
CombineSpaces& CombineSpaces::operator=(const CombineSpaces &other)
{
    CombineProperties::operator=(other);
    return *this;
}

/** Comparison operator */
bool CombineSpaces::operator==(const CombineSpaces &other) const
{
    return CombineProperties::operator==(other);
}

/** Comparison operator */
bool CombineSpaces::operator!=(const CombineSpaces &other) const
{
    return CombineProperties::operator!=(other);
}

/** Update this combined space by extracting the required space 
    properties from 'properties'
    
    \throw SireBase::missing_property
    \throw SireError::invalid_cast
    \throw SireError::incompatible_error
*/
void CombineSpaces::updateFrom(const Properties &properties)
{
    if (this->isEmpty())
        return;
        
    else if (this->count() == 1)
    {
        this->setCombinedProperty( properties.property(this->at(0))
                                             .asA<Space>() );
    }
    else
    {
        QList<SpacePtr> spaces;
        
        for (CombineProperties::const_iterator it = this->constBegin();
             it != this->constEnd();
             ++it)
        {
            spaces.append( properties.property(*it).asA<Space>() );
        }
        
        this->setCombinedProperty( CombinedSpace(spaces) );
    }
}

const char* CombineSpaces::typeName()
{
    return QMetaType::typeName( qMetaTypeId<CombineSpaces>() );
}

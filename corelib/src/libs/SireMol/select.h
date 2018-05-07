/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#ifndef SIREMOL_SELECT_H
#define SIREMOL_SELECT_H

#include "SireBase/property.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireMol
{
class Select;
}

QDataStream& operator<<(QDataStream&, const SireMol::Select&);
QDataStream& operator>>(QDataStream&, SireMol::Select&);

namespace SireMol
{

namespace parser
{
/** This is the base class of all of the select objects. It is a private
    object that should only be used by Select

    @author Christopher Woods
*/
class SelectEngine
{

public:
    SelectEngine();
    virtual ~SelectEngine();
};
} //end of namespace parser

/** This is the only publicly visible selector class. This provides a 
    front-end interface to selecting atoms and molecules
    
    @author Christopher Woods
*/
class SIREMOL_EXPORT Select : public SireBase::ConcreteProperty<Select,SireBase::Property>
{

friend QDataStream& ::operator<<(QDataStream&, const Select&);
friend QDataStream& ::operator>>(QDataStream&, Select&);

public:
    Select();
    Select(const QString &str);
    
    Select(const Select &other);
    
    ~Select();
    
    Select& operator=(const Select &other);
    
    bool operator==(const Select &other) const;
    bool operator!=(const Select &other) const;
    
    Select* clone() const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;

private:
    /** The actual search string */
    QString search_string;

    /** The underlying engine used to perform the selection */
    boost::shared_ptr<SireMol::parser::SelectEngine> e;
};

}

Q_DECLARE_METATYPE( SireMol::Select )

SIRE_EXPOSE_CLASS( SireMol::Select )

SIRE_END_HEADER

#endif

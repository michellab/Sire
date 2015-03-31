/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#ifndef SIREMOL_ELEMENT_H
#define SIREMOL_ELEMENT_H

#include "SireUnits/dimensions.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class Element;
}

QDataStream& operator<<(QDataStream&, const SireMol::Element&);
QDataStream& operator>>(QDataStream&, SireMol::Element&);

namespace SireMol
{

class ElementData;

/**
  * This class is used to represent a chemical element. The implementation
  * of this class is such that it is quick and easy to pass and copy
  * Elements, while the storage requirements are very low (just a single
  * pointer). The actual element data is held in a private ElementData
  * class that is declared in elementdb.cpp
  *
  * @author Christopher Woods
  */
class SIREMOL_EXPORT Element
{

public:
    Element();
    
    Element(QString element);
    Element(const char *element);
    Element(unsigned int nprotons);
    Element(int nprotons);
    
    Element(const Element &element);
    
    ~Element();
    
    static const char* typeName();
    
    const char* what() const
    {
        return Element::typeName();
    }
    
    bool operator==(const Element &other) const;
    bool operator!=(const Element &other) const;
    
    bool operator>(const Element &other) const;
    bool operator<(const Element &other) const;
    bool operator>=(const Element &other) const;
    bool operator<=(const Element &other) const;
    
    const Element& operator=(const Element &element);

    int nProtons() const;
    
    QString symbol() const;
    QString name() const;
    
    QString toString() const;
    
    SireUnits::Dimension::Length covalentRadius() const;
    SireUnits::Dimension::Length bondOrderRadius() const;
    SireUnits::Dimension::Length vdwRadius() const;
    
    int maxBonds() const;
    
    SireUnits::Dimension::MolarMass mass() const;
    
    float red() const;
    float green() const;
    float blue() const;
    
    int group() const;
    int period() const;
    
    bool nobleGas() const;
    bool halogen() const;
    bool biological() const;
    bool alkaliMetal() const;
    bool alkaliEarthMetal() const;
    bool transitionMetal() const;
    bool lanthanide() const;
    bool actinide() const;
    bool rareEarth() const;
    
    static Element elementWithMass(const SireUnits::Dimension::MolarMass &mass);
    
    static Element biologicalElement(const QString &name);
    
private:
    /** Pointer to the object containing the data for this element */
    ElementData *eldata;
};
  
#ifndef SIRE_SKIP_INLINE_FUNCTIONS

inline uint qHash(const Element &element)
{
    return element.nProtons();
}
        
/** Comparison operators. These are very quick */
inline bool Element::operator==(const Element &other) const
{
    return eldata == other.eldata;
}
    
/** Comparison operators. These are very quick */
inline bool Element::operator!=(const Element &other) const
{
    return eldata != other.eldata;
}

#endif //SIRE_SKIP_INLINE_FUNCTIONS

}

Q_DECLARE_METATYPE(SireMol::Element);
Q_DECLARE_TYPEINFO(SireMol::Element, Q_MOVABLE_TYPE);

SIRE_EXPOSE_CLASS( SireMol::Element )

SIRE_END_HEADER

#endif

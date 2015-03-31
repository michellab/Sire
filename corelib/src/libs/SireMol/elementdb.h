/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_ELEMENTDB_H
#define SIREMOL_ELEMENTDB_H

#include <QString>
#include <QHash>

#include "element.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

/**
This singleton class holds a database of element properties and information indexed by atomic number or element symbol. You must initialise this database before you can create any Element objects.
 
@author Christopher Woods
*/
class ElementDB
{
friend class Element;

public:
    /** You must call this function before you can access the database */
    static void initialise();
protected:
    ElementDB();
    ~ElementDB();
    
    /** Populate this database with data. This function is defined
        in element-data.h, which is automatically generated from the
        openbabel data. */
    void populate();

    /** Import an individual element to the database */
    void import(ElementData *element);
      
    /** Return the ElementData for the string 'element'. See the QString
        constructor for Element for a description of how 'element' is 
        interpreted. This will never return a null pointer. */
    ElementData* element(const QString &element) const;
    
    /** Return the element with proton number 'nprotons'. This will
        never return a null pointer. */
    ElementData* element(int nprotons) const;
    
    /** Pointer to the single instance of this class */
    static ElementDB *db;

    /** A hash indexing the elements by proton number */
    QHash<int, ElementData*> protonindex;
    
    /** A hash indexing the elements by symbol */
    QHash<QString, ElementData*> symbolindex;
    /** A hash indexing the elements by their name */
    QHash<QString, ElementData*> nameindex;
};

}

SIRE_END_HEADER

#endif

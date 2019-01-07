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

#include <QRegExp>
#include <QDataStream>
#include <cmath>

#include <limits>

#include "elementdb.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireStream;
using namespace SireMol;
using namespace SireUnits;
using namespace SireUnits::Dimension;

static const RegisterMetaType<Element> r_element(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Element &element)
{
    //this seems slightly over the top, adding magic and version to element,
    //as element has only a single quint32 data member. However, this may
    //change (e.g. perhaps to qint64 to allow negative proton numbers) and
    //the best way to help maintain backwards compatibility is if the
    //magic and version numbers are present

    writeHeader(ds, r_element, 1) << quint32( element.nProtons() );
    return ds;
}

/** Deserialise from a binary data stream */
QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Element &element)
{
    VersionID v = readHeader(ds, r_element);

    if (v == 1)
    {
        quint32 nprots;
        ds >> nprots;
        element = SireMol::Element(nprots);
    }
    else
        throw version_error(v, "1", r_element, CODELOC);

    return ds;
}

/** This is a basic (private) class that holds the actual element data */
namespace SireMol
{

class ElementData
{
public:
    ElementData(int protnum, QString name, QString symbol, int group,
                int period, double covalent_radii, double bond_order_radii,
                double vdw_radii, int maxbonds,
                double mass, double elec_neg, float red, float green, float blue);
    ~ElementData();

    Length cov_rad, bond_rad, vdw_rad;
    MolarMass mss;
    double electro;

    float r, g, b;
    QString symb;
    QString name;

    uchar protnum;  //none of these will exceed 255! (not in my lifetime, anyway!)
    uchar maxbonds;
    uchar group;
    uchar period;
};
}

ElementData::ElementData(int pnum, QString nam, QString sym, int grp, int per,
                         double crad, double brad, double vdw, int mxb,
                         double m, double elec, float rd, float grn, float blu)
        : cov_rad(crad * angstrom), 
          bond_rad(brad * angstrom), 
          vdw_rad(vdw * angstrom), 
          mss(m * g_per_mol), 
          electro(elec),
          r(rd), g(grn), b(blu), symb(sym), name(nam),
          protnum(static_cast<uchar>(pnum)), maxbonds(static_cast<uchar>(mxb)),
          group(static_cast<uchar>(grp)), period(static_cast<uchar>(per))
{}

ElementData::~ElementData()
{}

ElementDB* ElementDB::db = 0;

#include <QObject>
#include "element-data.h"

///////////
/////////// Implementation of Element
///////////

/** Construct a dummy element */
Element::Element()
{
    eldata = ElementDB::db->element(0);
}

/** Construct an element from the string 'element'. If the string
    contains 1 or 2 characters, then it is interpreted as an IUPAC
    chemical symbol (e.g. 'C', or 'Al'), else it is interpreted as
    the name of the element in the local language of the application.
    If the string cannot be interpreted then the element is set to
    the dummy element. Note 'ca' is Calcium, not C-alpha! */
Element::Element(QString element)
{
    eldata = ElementDB::db->element(element);
}

/** Overload so that 'const char*' is interpreted as a QString, and not
    as an unsigned int */
Element::Element(const char *element)
{
    eldata = ElementDB::db->element( QString(element) );
}

/** Construct an element with proton number 'nprot'. If there is no
    element with this number of protons, or if the proton number is 0,
    then a dummy elements is constructed */
Element::Element(unsigned int nprots)
{
    eldata = ElementDB::db->element(nprots);
}

/** Overload to disambiguate the call to Element(int) */
Element::Element(int nprots)
{
    eldata = ElementDB::db->element(nprots);
}

/** Copy constructor. This is very quick as it involves copying
    only a single pointer. */
Element::Element(const Element &element) : eldata(element.eldata)
{}

/** Destructor */
Element::~Element()
{}

/** Assignment operator */
const Element& Element::operator=(const Element &element)
{
    eldata = element.eldata;
    return *this;
}

/** Return the number of protons in the element */    
int Element::nProtons() const
{
    return eldata->protnum;
}

/** Return the IUPAC symbol for the element */
QString Element::symbol() const
{
    return eldata->symb;
}

/** Return the name of the element in the local language */
QString Element::name() const
{
    return eldata->name;
}

/** Return the group number of this element (IUPAC group, from 1-18)
    (lanthanides and actinides have a group number of 0 - this should
    not be too big a problem as I would be surprised to hear of anyone
    using this code to simulate them...) */
int Element::group() const
{
    return eldata->group;
}

/** Return the period (the row number) of the element (IUPAC period, from 1-7) */
int Element::period() const
{
    return eldata->period;
}

/** Return the element's covalent radius */
Length Element::covalentRadius() const
{
    return eldata->cov_rad;
}

/** Return the bond order radius */
Length Element::bondOrderRadius() const
{
    return eldata->bond_rad;
}

/** Return the van der waals radius */
Length Element::vdwRadius() const
{
    return eldata->vdw_rad;
}

/** Return the maximum number of simultaneous bonds that this
    element can form */
int Element::maxBonds() const
{
    return eldata->maxbonds;
}

/** Return the average mass of this element */
MolarMass Element::mass() const
{
    return eldata->mss;
}

/** Return the element's electronegativity */
/*double Element::electroNegativity() const
{
    return eldata->electro;
}*/

/** Return the red colour components (0.0->1.0) for
    the colour of this element */
float Element::red() const
{
    return eldata->r;
}

/** Return the green colour components (0.0->1.0) for
    the colour of this element */
float Element::green() const
{
    return eldata->g;
}

/** Return the blue colour components (0.0->1.0) for
    the colour of this element */
float Element::blue() const
{
    return eldata->b;
}

/**
  * Now the implementation of the ElementDB class
  *
  */

//use a static global variable to initialise the ElementDB at the point
//the SireMol library is loaded
class EDB_loader
{
public:
    EDB_loader()
    {
        ElementDB::initialise();
    }

    ~EDB_loader()
    {}
};

static EDB_loader loader;

void ElementDB::initialise()
{
    if (db)
        return;

    db = new ElementDB();
}

ElementDB::ElementDB()
{
    this->populate();
}

ElementDB::~ElementDB()
{
    //delete the data
    for (QHash<QString,ElementData*>::iterator it = symbolindex.begin();
            it != symbolindex.end();
            ++it)
    {
        delete it.value();
    }

    //clear the hashes
    protonindex.clear();
    symbolindex.clear();
    nameindex.clear();
}

void ElementDB::import(ElementData *element)
{
    if (not element)
        return;

    symbolindex.insert(element->symb.toLower(),element);
    protonindex.insert(element->protnum,element);
    nameindex.insert(element->name.toLower(),element);
}

ElementData* ElementDB::element(int z) const
{
    if (protonindex.contains(z))
        return protonindex.value(z);
    else
        return protonindex.value(0);
}

ElementData* ElementDB::element(const QString &s) const
{
    const QRegExp numregexp("[0-9]");
    //lowercase the name, remove spaces and numbers
    QString el = s.toLower().remove(numregexp).trimmed();

    if (symbolindex.contains(el))
        return symbolindex.value(el);
    else
    {
        //try and interpret this string as the element name
        if (nameindex.contains(el))
            return nameindex.value(el);
        else
        {
            //try to guess the name...
            if (symbolindex.contains(el.left(2)))
                return symbolindex.value(el.left(2));
            else if (symbolindex.contains(el.left(1)))
                return symbolindex.value(el.left(1));
            else
                //no luck, return a dummy atom
                return protonindex.value(0);
        }
    }
}

/** Return a biological element that has been guessed from the passed name.
    Note that if no biological element was guessed, then the nearest
    non-biological element match is used. A biological element is one that
    is in the first couple of rows (proton number < 18) and is not a noble gas. */
Element Element::biologicalElement(const QString &name)
{
    //guess an element with this name...
    Element elmnt(name);

    //is this a biological element? - if so, return it!
    if (elmnt.biological())
        return elmnt;

    //try to guess the atom from just the first two letters...
    Element elmnt2(name.left(2));

    if (elmnt2.biological())
        return elmnt2;

    //try to guess the atom from just the first letter...
    Element elmnt3(name.left(1));

    if (elmnt3.biological())
        return elmnt3;

    //we couldn't find anything - return the original, non-biological guess
    return elmnt;
}

/** Return whether or not this is a noble gas */
bool Element::nobleGas() const
{
    return group() == 18;
}

/** Return whether or not this is a halogen */
bool Element::halogen() const
{
    return group() == 17;
}

/** Return whether or not this is biological
    (in first three periods and not a noble gas, or a halogen)
    (this does preclude iron, potassium and calcium, which are
    rather biological... :-) */
bool Element::biological() const
{
    return (period() <= 3 and not nobleGas()) or halogen();
}

/** Return whether or not this is an alkali metal (group 1 or 2) */
bool Element::alkaliMetal() const
{
    return nProtons() != 1 and group() == 1;
}

/** Return whether or not this is an alkali earth metal (group 2) */
bool Element::alkaliEarthMetal() const
{
    return group() == 2;
}

/** Return whether or not this is a transition metal */
bool Element::transitionMetal() const
{
    return group() >= 3 and group() <= 12;
}

/** Return whether or not this is a lanthanide */
bool Element::lanthanide() const
{
    return nProtons() >= 58 and nProtons() <= 71;
}

/** Return whether or not this is an actinide */
bool Element::actinide() const
{
    return nProtons() >= 90 and nProtons() <= 103;
}

/** Return whether or not this is a rare earth element (e.g. a lanthanide or actinide) */
bool Element::rareEarth() const
{
    return lanthanide() or actinide();
}

/** Sorting operators. Elements are compared based on their
    proton numbers, with elements with greater numbers being higher.
    The functions are also very quick. */
bool Element::operator>(const Element &other) const
{
    return eldata->protnum > other.eldata->protnum;
}

/** Sorting operators. Elements are compared based on their
    proton numbers, with elements with greater numbers being higher.
    The functions are also very quick. */
bool Element::operator<(const Element &other) const
{
    return eldata->protnum < other.eldata->protnum;
}

/** Sorting operators. Elements are compared based on their
    proton numbers, with elements with greater numbers being higher.
    The functions are also very quick. */
bool Element::operator>=(const Element &other) const
{
    return eldata->protnum >= other.eldata->protnum;
}

/** Sorting operators. Elements are compared based on their
    proton numbers, with elements with greater numbers being higher.
    The functions are also very quick. */
bool Element::operator<=(const Element &other) const
{
    return eldata->protnum <= other.eldata->protnum;
}

/** Return a string representation of the Element */
QString Element::toString() const
{
    return QObject::tr("%1 (%2, %3)").arg(name(),symbol()).arg(nProtons());
}

/** Return an element which has the closest mass to 'mass' (in atomic 
    mass units, g mol-1) */
Element Element::elementWithMass(const MolarMass &molar_mass)
{
    double mass = molar_mass.to( g_per_mol );

    //round up the mass to the nearest integer to see if we can match to
    //a core, biological type
    int i_mass = int( mass + 0.5 );    //this rounds to the nearest int

    switch(i_mass)
    {
        case 1:  //hydrogen
            return Element(1);
        case 12:  //carbon
            return Element(6);
        case 16:  //oxygen
            return Element(8);
        case 14:  //nitrogen
            return Element(7);
        case 32:  //sulfur
            return Element(16);
    }

    //the quick test failed, so now we will do a long lookup
    double diff = std::numeric_limits<double>::max();
    ElementData *bestmatch = 0;

    //loop over all of the elements...
    foreach (ElementData *element, ElementDB::db->protonindex.values())
    {
        double testdiff = std::abs(mass - element->mss);
        if (testdiff < diff)
        {
            bestmatch = element;
            diff = testdiff;

            if (diff == 0.0)
                return Element(bestmatch->protnum);
        }
    }

    //return the closest match
    if (bestmatch)
        return Element(bestmatch->protnum);
    else
        return Element(0);
}

const char* Element::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Element>() );
}

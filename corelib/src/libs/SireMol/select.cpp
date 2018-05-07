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

#include "SireMol/select.h"
#include "SireMol/parser.h"

#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;
using namespace SireMol;

///////////
/////////// Implementation of SelectBase
///////////

static const RegisterMetaType<SelectBase> r_base( MAGIC_ONLY, SelectBase::typeName() );

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const SelectBase &select)
{
    writeHeader(ds, r_base, 1);
    ds << static_cast<const Property&>(select);
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, SelectBase &select)
{
    VersionID v = readHeader(ds, r_base);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(select);
    }
    else
        throw version_error( v, "1", r_base, CODELOC );
    
    return ds;
}

/** Constructor */
SelectBase::SelectBase() : Property()
{}

/** Copy constructor */
SelectBase::SelectBase(const SelectBase &other) : Property(other)
{}

/** Destructor */
SelectBase::~SelectBase()
{}

const char* SelectBase::typeName()
{
    return "SireMol::SelectBase";
}

SelectBase& SelectBase::operator=(const SelectBase &other)
{
    Property::operator=(other);
    return *this;
}

bool SelectBase::operator==(const SelectBase &other) const
{
    return Property::operator==(other);
}

bool SelectBase::operator!=(const SelectBase &other) const
{
    return not SelectBase::operator==(other);
}

const Select& SelectBase::null()
{
    return *(create_shared_null<Select>());
}

///////////
/////////// Implementation of Select
///////////

static const RegisterMetaType<Select> r_select;

QDataStream SIREMOL_EXPORT &operator<<(QDataStream &ds, const Select &select)
{
    writeHeader(ds, r_select, 1);
    
    SharedDataStream sds(ds);
    
    sds << select.search_string << select.p
        << static_cast<const SelectBase&>(select);
    
    return ds;
}

QDataStream SIREMOL_EXPORT &operator>>(QDataStream &ds, Select &select)
{
    VersionID v = readHeader(ds, r_select);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> select.search_string >> select.p
            >> static_cast<SelectBase&>(select);
    }
    else
        throw version_error(v, "1", r_select, CODELOC);
    
    return ds;
}

/** Construct an empty selection (will select nothing) */
Select::Select() : ConcreteProperty<Select,SelectBase>()
{}

/** Construct a selection based on the passed string */
Select::Select(const QString &str) : ConcreteProperty<Select,SelectBase>()
{
    p = parse(str);
    search_string = str;
}

/** Copy constructor */
Select::Select(const Select &other)
       : ConcreteProperty<Select,SelectBase>(other),
         search_string(other.search_string), p(other.p)
{}

/** Destructor */
Select::~Select()
{}

Select& Select::operator=(const Select &other)
{
    search_string = other.search_string;
    p = other.p;
    SelectBase::operator=(other);
    return *this;
}

bool Select::operator==(const Select &other) const
{
    return p == other.p and SelectBase::operator==(other);
}

bool Select::operator!=(const Select &other) const
{
    return not Select::operator==(other);
}

Select* Select::clone() const
{
    return new Select(*this);
}

const char* Select::what() const
{
    return Select::typeName();
}

const char* Select::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Select>() );
}

QString Select::toString() const
{
    if (p.isEmpty())
    {
        return QObject::tr("Select::null");
    }
    else
        return search_string;
}

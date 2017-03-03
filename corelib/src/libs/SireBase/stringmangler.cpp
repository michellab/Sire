/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2008  Christopher Woods
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

#include <QMutex>

#include "stringmangler.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireBase;
using namespace SireStream;

////////
//////// Implementation of StringMangler
////////

static const RegisterMetaType<StringMangler> r_stringmangler( MAGIC_ONLY,
                                                  "SireBase::StringMangler" );
                                                  
/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds, 
                                        const StringMangler &stringmangler)
{
    writeHeader(ds, r_stringmangler, 1);
    ds << static_cast<const Property&>(stringmangler);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, StringMangler &stringmangler)
{
    VersionID v = readHeader(ds, r_stringmangler);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(stringmangler);
    }
    else
        throw version_error( v, "1", r_stringmangler, CODELOC );
    
    return ds;
}

/** Constructor */
StringMangler::StringMangler() : Property()
{}

/** Copy constructor */
StringMangler::StringMangler(const StringMangler &other)
              : Property(other)
{}

/** Destructor */
StringMangler::~StringMangler()
{}

/** Mangle the input string */
QString StringMangler::operator()(const QString &input) const
{
    return this->mangle(input);
}

////////
//////// Implementation of NoMangling
////////

static const RegisterMetaType<NoMangling> r_nomangle;

/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds,
                                        const NoMangling &nomangle)
{
    writeHeader(ds, r_nomangle, 1);
    
    ds << static_cast<const StringMangler&>(nomangle);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, NoMangling &nomangle)
{
    VersionID v = readHeader(ds, r_nomangle);
    
    if (v == 1)
    {
        ds >> static_cast<StringMangler&>(nomangle);
    }
    else
        throw version_error(v, "1", r_nomangle, CODELOC);
        
    return ds;
}

/** Constructor */
NoMangling::NoMangling() : ConcreteProperty<NoMangling,StringMangler>()
{}

/** Copy constructor */
NoMangling::NoMangling(const NoMangling &other)
           : ConcreteProperty<NoMangling,StringMangler>(other)
{}

/** Destructor */
NoMangling::~NoMangling()
{}

/** Copy assignment operator */
NoMangling& NoMangling::operator=(const NoMangling&)
{
    return *this;
}

/** Comparison operator */
bool NoMangling::operator==(const NoMangling&) const
{
    return true;
}

/** Comparison operator */
bool NoMangling::operator!=(const NoMangling&) const
{
    return false;
}

const char* NoMangling::typeName()
{
    return QMetaType::typeName( qMetaTypeId<NoMangling>() );
}

/** Mangle the string - remove all initial and trailing spaces */
QString NoMangling::mangle(const QString &input) const
{
    return input;
}

const NoMangling& StringMangler::null()
{
    return *(create_shared_null<NoMangling>());
}

////////
//////// Implementation of TrimString
////////

static const RegisterMetaType<TrimString> r_trimstring;

/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds,
                                        const TrimString &trimstring)
{
    writeHeader(ds, r_trimstring, 1);
    
    ds << static_cast<const StringMangler&>(trimstring);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, TrimString &trimstring)
{
    VersionID v = readHeader(ds, r_trimstring);
    
    if (v == 1)
    {
        ds >> static_cast<StringMangler&>(trimstring);
    }
    else
        throw version_error(v, "1", r_trimstring, CODELOC);
        
    return ds;
}

/** Constructor */
TrimString::TrimString() : ConcreteProperty<TrimString,StringMangler>()
{}

/** Copy constructor */
TrimString::TrimString(const TrimString &other)
           : ConcreteProperty<TrimString,StringMangler>(other)
{}

/** Destructor */
TrimString::~TrimString()
{}

/** Copy assignment operator */
TrimString& TrimString::operator=(const TrimString&)
{
    return *this;
}

/** Comparison operator */
bool TrimString::operator==(const TrimString&) const
{
    return true;
}

/** Comparison operator */
bool TrimString::operator!=(const TrimString&) const
{
    return false;
}

const char* TrimString::typeName()
{
    return QMetaType::typeName( qMetaTypeId<TrimString>() );
}

/** Mangle the string - remove all initial and trailing spaces */
QString TrimString::mangle(const QString &input) const
{
    return input.trimmed();
}

////////
//////// Implementation of UpperCaseString
////////

static const RegisterMetaType<UpperCaseString> r_upperstring;

/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds,
                                        const UpperCaseString &upperstring)
{
    writeHeader(ds, r_upperstring, 1);
    
    ds << static_cast<const StringMangler&>(upperstring);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, UpperCaseString &upperstring)
{
    VersionID v = readHeader(ds, r_upperstring);
    
    if (v == 1)
    {
        ds >> static_cast<StringMangler&>(upperstring);
    }
    else
        throw version_error(v, "1", r_upperstring, CODELOC);
        
    return ds;
}

/** Constructor */
UpperCaseString::UpperCaseString() : ConcreteProperty<UpperCaseString,StringMangler>()
{}

/** Copy constructor */
UpperCaseString::UpperCaseString(const UpperCaseString &other)
           : ConcreteProperty<UpperCaseString,StringMangler>(other)
{}

/** Destructor */
UpperCaseString::~UpperCaseString()
{}

/** Copy assignment operator */
UpperCaseString& UpperCaseString::operator=(const UpperCaseString&)
{
    return *this;
}

/** Comparison operator */
bool UpperCaseString::operator==(const UpperCaseString&) const
{
    return true;
}

/** Comparison operator */
bool UpperCaseString::operator!=(const UpperCaseString&) const
{
    return false;
}

const char* UpperCaseString::typeName()
{
    return QMetaType::typeName( qMetaTypeId<UpperCaseString>() );
}

/** Mangle the string - remove all initial and trailing spaces */
QString UpperCaseString::mangle(const QString &input) const
{
    return input.trimmed().toUpper();
}

////////
//////// Implementation of LowerCaseString
////////

static const RegisterMetaType<LowerCaseString> r_lowerstring;

/** Serialise to a binary datastream */
QDataStream SIREBASE_EXPORT &operator<<(QDataStream &ds,
                                        const LowerCaseString &lowerstring)
{
    writeHeader(ds, r_lowerstring, 1);
    
    ds << static_cast<const StringMangler&>(lowerstring);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREBASE_EXPORT &operator>>(QDataStream &ds, LowerCaseString &lowerstring)
{
    VersionID v = readHeader(ds, r_lowerstring);
    
    if (v == 1)
    {
        ds >> static_cast<StringMangler&>(lowerstring);
    }
    else
        throw version_error(v, "1", r_lowerstring, CODELOC);
        
    return ds;
}

/** Constructor */
LowerCaseString::LowerCaseString() : ConcreteProperty<LowerCaseString,StringMangler>()
{}

/** Copy constructor */
LowerCaseString::LowerCaseString(const LowerCaseString &other)
           : ConcreteProperty<LowerCaseString,StringMangler>(other)
{}

/** Destructor */
LowerCaseString::~LowerCaseString()
{}

/** Copy assignment operator */
LowerCaseString& LowerCaseString::operator=(const LowerCaseString&)
{
    return *this;
}

/** Comparison operator */
bool LowerCaseString::operator==(const LowerCaseString&) const
{
    return true;
}

/** Comparison operator */
bool LowerCaseString::operator!=(const LowerCaseString&) const
{
    return false;
}

const char* LowerCaseString::typeName()
{
    return QMetaType::typeName( qMetaTypeId<LowerCaseString>() );
}

/** Mangle the string - remove all initial and trailing spaces */
QString LowerCaseString::mangle(const QString &input) const
{
    return input.trimmed().toLower();
}

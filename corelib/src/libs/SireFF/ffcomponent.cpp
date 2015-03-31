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

#include "ffcomponent.h"
#include "ff.h"

#include "SireStream/datastream.h"

#include "SireError/errors.h"

#include <QRegExp>

using namespace SireFF;
using namespace SireCAS;
using namespace SireStream;

static const RegisterMetaType<FFComponent> r_ffcomp( MAGIC_ONLY,
                                                     "SireFF::FFComponent" );
                                                     
/** Serialise to a binary datastream */
QDataStream SIREFF_EXPORT &operator<<(QDataStream &ds, const FFComponent &ffcomp)
{
    writeHeader(ds, r_ffcomp, 1);
    ds << static_cast<const Symbol&>(ffcomp);

    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREFF_EXPORT &operator>>(QDataStream &ds, FFComponent &ffcomp)
{
    VersionID v = readHeader(ds, r_ffcomp);
    
    if (v == 1)
    {
        ds >> static_cast<Symbol&>(ffcomp);
    }
    else
        throw version_error(v, "1", r_ffcomp, CODELOC);
        
    return ds;
}

static QRegExp name_regexp( "E\\_\\{(.+)\\}\\^\\{(.+)\\}" );

/** Constructor */
FFComponent::FFComponent(const FFName &ffname, const QString &name)
            : Symbol( FFComponent::symbolName(ffname, name) )
{}

/** Construct from a passed symbol 

    \throw SireError::incompatible_error
*/
FFComponent::FFComponent(const Symbol &symbol, const QString &name)
            : Symbol(symbol)
{
    QRegExp local_copy = name_regexp;
    
    if (local_copy.indexIn( Symbol::toString() ) == -1)
        //we could not interpret the symbol
        throw SireError::incompatible_error( QObject::tr(
            "The symbol \"%1\" is not a valid FFComponent symbol. "
            "FFComponent symbols have the form E^{<FF_NAME>}^{<COMPONENT_NAME>} "
            "where <FF_NAME> is the name of the forcefield that contains this "
            "component and <COMPONENT_NAME> is the name of the component itself.")
                .arg(symbol.toString()), CODELOC );
                                
    if (name != local_copy.cap(2))
        throw SireError::incompatible_error( QObject::tr(
            "The symbol \"%1\" is for the wrong component! (%2). "
            "The correct component is called \"%3\".")
                .arg(symbol.toString(), local_copy.cap(2), name),
                    CODELOC );
}

/** Copy constructor */
FFComponent::FFComponent(const FFComponent &other)
            : Symbol(other)
{}

/** Destructor */
FFComponent::~FFComponent()
{}

/** Return the name of the forcefield that owns this component */
FFName FFComponent::forceFieldName() const
{
    QRegExp local_copy = name_regexp;
    
    if (local_copy.indexIn( Symbol::toString() ) == -1)
        //we could not interpret the name!
        throw SireError::program_bug( QObject::tr(
            "Could not interpret this symbol (%1) as an FFComponent.")
                .arg(Symbol::toString()), CODELOC );
                
    return FFName( local_copy.cap(1) );
}

/** Return the name of the component of the potential energy
    surface that this symbol represents */
QString FFComponent::componentName() const
{
    QRegExp local_copy = name_regexp;
    
    if (name_regexp.indexIn( Symbol::toString() ) == -1)
        //we could not interpret the name!
        throw SireError::program_bug( QObject::tr(
            "Could not interpret this symbol (%1) as an FFComponent.")
                .arg(Symbol::toString()), CODELOC );

    return local_copy.cap(1);
}

/** Protected function used to set the energy of the component with symbol
    'symbol' in the forcefield 'ff' to the value 'value' */
void FFComponent::setEnergy(FF &ff, const Symbol &symbol, double value) const
{
    ff.setComponent(symbol, value);
}

/** Protected function used to change the energy of the component with symbol
    'symbol' in the forcefield 'ff' by 'delta' */
void FFComponent::changeEnergy(FF &ff, const Symbol &symbol, double delta) const
{
    ff.changeComponent(symbol, delta);
}

//////
////// Implementation of SingleComponent
//////

static const RegisterMetaType<SingleComponent> r_single;

/** Serialise to a binary datastream */
QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const SingleComponent &single)
{
    writeHeader(ds, r_single, 1);
    ds << static_cast<const FFComponent&>(single);
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, SingleComponent &single)
{
    VersionID v = readHeader(ds, r_single);
    
    if (v == 1)
    {
        ds >> static_cast<FFComponent&>(single);
    }
    else
        throw version_error(v, "1", r_single, CODELOC);
        
    return ds;
}

/** Constructor */
SingleComponent::SingleComponent(const FFName &ffname)
                : FFComponent(ffname, QLatin1String("total"))
{}

/** Construct using the passed forcefield name and suffix */
SingleComponent::SingleComponent(const FFName &ffname, const QString &suffix)
                 : FFComponent(ffname, QString("total_{%1}").arg(suffix))
{}

/** Construct from a symbol

    \throw SireError::incompatible_error
*/
SingleComponent::SingleComponent(const SireCAS::Symbol &symbol)
                 : FFComponent(symbol, QLatin1String("total"))
{}

/** Copy constructor */  
SingleComponent::SingleComponent(const SingleComponent &other)
                 : FFComponent(other)
{}
  
/** Destructor */  
SingleComponent::~SingleComponent()
{}

/** Set the energy in the forcefield 'ff' to equal to the passed SingleEnergy */
void SingleComponent::setEnergy(FF &ff, const SingleEnergy &nrg) const
{
    FFComponent::setEnergy(ff, this->total(), nrg);
}

/** Change the energy in the forcefield 'ff' by 'delta' */
void SingleComponent::changeEnergy(FF &ff, const SingleEnergy &delta) const
{
    FFComponent::changeEnergy(ff, this->total(), delta);
}

const char* SingleComponent::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SingleComponent>() );
}

/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIREBASE_CPUID_H
#define SIREBASE_CPUID_H

#include "property.h"

#include <QHash>

SIRE_BEGIN_HEADER

namespace SireBase
{
class CPUID;
}

SIREBASE_EXPORT QDataStream& operator<<(QDataStream&, const SireBase::CPUID&);
SIREBASE_EXPORT QDataStream& operator>>(QDataStream&, SireBase::CPUID&);

namespace SireBase
{

/** This class obtains and displays the capabilities and ID of
    the CPU at runtime
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT CPUID : public ConcreteProperty<CPUID,Property>
{

friend SIREBASE_EXPORT QDataStream& ::operator<<(QDataStream&, const CPUID&);
friend SIREBASE_EXPORT QDataStream& ::operator>>(QDataStream&, CPUID&);

public:
    CPUID();
    CPUID(const CPUID &other);
    
    ~CPUID();
    
    CPUID& operator=(const CPUID &other);
    
    bool operator==(const CPUID &other) const;
    bool operator!=(const CPUID &other) const;
    
    CPUID* clone() const;
    
    const char* what() const;
    static const char* typeName();
    
    QString toString() const;
    
    bool supports(const QString &feature) const;

    QStringList supportableFeatures() const;
    
    QStringList supportedFeatures() const;
    
    QString vendor() const;
    QString brand() const;
    
    int clockSpeed() const;
    int numCores() const;
    
    bool supportsSSE2() const;
    bool supportsAVX() const;
    
private:
    QHash<QString,QString>* getCPUID();

    /** A simple dictionary of key-value pairs for the CPU */
    QHash<QString,QString> props;
    
    static QHash<QString,QString> *global_props;
};

}

Q_DECLARE_METATYPE( SireBase::CPUID )

SIRE_EXPOSE_CLASS( SireBase::CPUID )

SIRE_END_HEADER

#endif

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

#ifndef SIREFF_FFDETAIL_H
#define SIREFF_FFDETAIL_H

#include "SireBase/property.h"
#include "SireBase/properties.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class FFDetail;
}

QDataStream& operator<<(QDataStream&, const SireFF::FFDetail&);
QDataStream& operator>>(QDataStream&, SireFF::FFDetail&);

namespace SireFF
{

using SireBase::Property;
using SireBase::Properties;

/** This is the base class of the classes that provide details about
    the forcefield of a molecule or system. The details include the 
    canonical name of the forcefield. When specialised into
    SireMM::MMDetail it describes whether it uses LJ parameters,
    combining rules, default 1-4 scale factors, functional forms
    of internal parameters etc. etc.
    
    @author Christopher Woods
*/
class SIREFF_EXPORT FFDetail : public SireBase::Property
{

friend QDataStream& ::operator<<(QDataStream&, const FFDetail&);
friend QDataStream& ::operator>>(QDataStream&, FFDetail&);

public:
    FFDetail();
    FFDetail(const QString &name);
    
    FFDetail(const FFDetail &other);
    
    virtual ~FFDetail();
    
    virtual FFDetail* clone() const=0;
    
    FFDetail& operator=(const FFDetail &other);
    
    bool operator==(const FFDetail &other) const;
    bool operator!=(const FFDetail &other) const;
    
    static const char* typeName();
    
    QString name() const;

    Properties properties() const;

    bool isNull() const;

    static SireBase::PropertyPtr get(QString forcefield);
    static QStringList forcefields();

protected:
    void setProperty(const QString &key, const Property &value);
    const Property& property(const QString &key) const;

    static SireBase::PropertyPtr registerForceField(const FFDetail &ff);

private:
    /** The store of all of the forcefield properties */
    Properties props;
};

}

SIRE_EXPOSE_CLASS( SireFF::FFDetail )

SIRE_END_HEADER

#endif

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

#ifndef SIRECAS_CONSTANT_H
#define SIRECAS_CONSTANT_H

#include "exbase.h"

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Constant;
}

QDataStream& operator<<(QDataStream&, const SireCAS::Constant&);
QDataStream& operator>>(QDataStream&, SireCAS::Constant&);

namespace SireCAS
{

/**
This class represents a constant value (e.g. a number).

@author Christopher Woods
*/
class SIRECAS_EXPORT Constant : public ExBase
{
public:
    Constant();
    Constant(const Constant &other);

    ~Constant();

    ///////
    /////// Virtual functions - you may wish to override these
    /////// in your derived class
    ///////

    Expression differentiate(const Symbol &symbol) const;
    Expression integrate(const Symbol &symbol) const;

    ///////
    /////// Pure-virtual functions - these must be overridden
    /////// in your derived class
    ///////

    bool operator==(const ExBase &other) const;

    uint hash() const;

    static const char* typeName();

    const char* what() const
    {
        return Constant::typeName();
    }

    Constant* clone() const;

    QString toString() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

    Expression substitute(const Identities &identities) const;

    Symbols symbols() const;
    Functions functions() const;

    Expressions children() const;

    QList<Factor> expand(const Symbol &symbol) const;
};

}

Q_DECLARE_METATYPE(SireCAS::Constant)

SIRE_EXPOSE_CLASS( SireCAS::Constant )

SIRE_END_HEADER

#endif

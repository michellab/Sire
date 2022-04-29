/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
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

#ifndef SIREMOL_EVALUATORM_H
#define SIREMOL_EVALUATORM_H

#include "evaluator.h"
#include "partialmolecule.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
class EvaluatorM;
class SelectorMol;

template<class T>
class SelectorM;
}

SIREMOL_EXPORT QDataStream& operator<<(QDataStream&, const SireMol::EvaluatorM&);
SIREMOL_EXPORT QDataStream& operator>>(QDataStream&, SireMol::EvaluatorM&);

namespace SireMol
{

/** This is a multi-molecule version of the evaluator class */
class SIREMOL_EXPORT EvaluatorM
    : public SireBase::ConcreteProperty<EvaluatorM, SireBase::Property>
{

friend SIREMOL_EXPORT QDataStream& ::operator<<(QDataStream&, const EvaluatorM&);
friend SIREMOL_EXPORT QDataStream& ::operator>>(QDataStream&, EvaluatorM&);

public:
    EvaluatorM();
    EvaluatorM(const SelectorMol &mols);

    template<class T>
    EvaluatorM(const SelectorM<T> &views);

    EvaluatorM(const EvaluatorM &other);

    virtual ~EvaluatorM();

    static const char* typeName();

    virtual const char* what() const
    {
        return EvaluatorM::typeName();
    }

    virtual EvaluatorM* clone() const
    {
        return new EvaluatorM(*this);
    }

    EvaluatorM& operator=(const EvaluatorM &other);

    bool operator==(const EvaluatorM &other) const;
    bool operator!=(const EvaluatorM &other) const;

    QString toString() const;

    bool isEmpty() const;

    int nAtoms() const;
    int nMolecules() const;

    SireUnits::Dimension::MolarMass mass() const;
    SireUnits::Dimension::MolarMass mass(const SireBase::PropertyMap &map) const;

    SireUnits::Dimension::Charge charge() const;
    SireUnits::Dimension::Charge charge(const PropertyMap &map) const;

    SireMaths::Vector center() const;
    SireMaths::Vector center(const PropertyMap &map) const;

    SireVol::AABox aaBox() const;
    SireVol::AABox aaBox(const PropertyMap &map) const;

    SireMaths::Sphere boundingSphere() const;
    SireMaths::Sphere boundingSphere(const PropertyMap &map) const;

    SireMaths::Vector centroid() const;
    SireMaths::Vector centroid(const PropertyMap &map) const;

    SireMaths::Vector centerOfGeometry() const;
    SireMaths::Vector centerOfGeometry(const PropertyMap &map) const;

    SireMaths::Vector centerOfMass() const;
    SireMaths::Vector centerOfMass(const PropertyMap &map) const;

protected:
    /** The actual views */
    QList<SireMol::PartialMolecule> vws;
};

} // end of namespace SireMol

Q_DECLARE_METATYPE( SireMol::EvaluatorM )

SIRE_EXPOSE_CLASS( SireMol::EvaluatorM )

SIRE_END_HEADER

#endif

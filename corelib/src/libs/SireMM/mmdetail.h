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

#ifndef SIREMM_MMDETAIL_H
#define SIREMM_MMDETAIL_H

#include "SireFF/ffdetail.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
class MMDetail;
}

QDataStream& operator<<(QDataStream&, const SireMM::MMDetail&);
QDataStream& operator>>(QDataStream&, SireMM::MMDetail&);

namespace SireMM
{

/** This class holds most of the data the describes different types
    of molecular mechanics forcefields, e.g. Amber FF99, OPLS etc.
    
    @author Christopher Woods
*/
class SIREMM_EXPORT MMDetail : public SireBase::ConcreteProperty<MMDetail,SireFF::FFDetail>
{

friend QDataStream& ::operator<<(QDataStream&, const MMDetail&);
friend QDataStream& ::operator>>(QDataStream&, MMDetail&);

public:
    MMDetail();
    MMDetail(QString name, QString combining_rules,
             double scale14elec, double scale14vdw,
             QString elecstyle, QString vdwstyle,
             QString bondstyle, QString anglestyle,
             QString dihedralstyle);
    
    MMDetail(const MMDetail &other);
    
    ~MMDetail();
    
    MMDetail& operator=(const MMDetail &other);
    
    bool operator==(const MMDetail &other) const;
    bool operator!=(const MMDetail &other) const;
    
    MMDetail* clone() const;
    
    static const char* typeName();
    const char* what() const;
    
    QString toString() const;
    
    QString combiningRules() const;
    bool usesArithmeticCombiningRules() const;
    bool usesGeometricCombiningRules() const;
    
    double electrostatic14ScaleFactor() const;
    double vdw14ScaleFactor() const;
    
    bool usesCoulombCharges() const;
    QString electrostaticStyle() const;
    
    bool usesLJTerm() const;
    QString vdwStyle() const;
    
    bool usesHarmonicBonds() const;
    QString bondStyle() const;
    
    bool usesHarmonicAngles() const;
    QString angleStyle() const;
    
    bool usesCosineDihedrals() const;
    QString dihedralStyle() const;
    
    bool isAmberStyle() const;
    
    bool isCompatibleWith(const FFDetail &other) const;
    
    static MMDetail guessFrom(QString combrule, QString elecstyle, QString vdwstyle,
                              double elec14, double vdw14, QString bondstyle,
                              QString anglestyle, QString dihedralstyle);
};

}

Q_DECLARE_METATYPE( SireMM::MMDetail )

SIRE_EXPOSE_CLASS( SireMM::MMDetail )

SIRE_END_HEADER

#endif

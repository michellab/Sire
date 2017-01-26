/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include <QFile>
#include <QTextStream>
#include <QHash>

#include "amber2.h"

#include "SireBase/parallel.h"
#include "SireSystem/system.h"
#include "SireError/errors.h"

using namespace SireIO;
using namespace SireMol;
using namespace SireSystem;
using namespace SireStream;
using namespace SireBase;

// The partial charges in the top file are not in electrons
static const double AMBERCHARGECONV = 18.2223;

static const double AMBER14COUL = 1.0 / 1.2 ;
static const double AMBER14LJ = 0.50 ;

static const RegisterMetaType<Amber2> r_amber2(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const Amber2 &amber2)
{
    writeHeader(ds, r_amber2, 1);
    
    ds << amber2.coul_14scl << amber2.lj_14scl;
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, Amber2 &amber2)
{
    VersionID v = readHeader(ds, r_amber2);

    if (v == 1)
    {
        ds >> amber2.coul_14scl >> amber2.lj_14scl;
    }
    else
        throw version_error( v, "1", r_amber2, CODELOC );

    return ds;
}

/** Constructor */
Amber2::Amber2() : coul_14scl(AMBER14COUL), lj_14scl(AMBER14LJ)
{}

/** Copy constructor */
Amber2::Amber2(const Amber2 &other)
       : coul_14scl(other.coul_14scl), lj_14scl(other.lj_14scl)
{}

/** Destructor */
Amber2::~Amber2()
{}

/** Copy assignment operator */
Amber2& Amber2::operator=(const Amber2 &other)
{
    coul_14scl = other.coul_14scl;
    lj_14scl = other.lj_14scl;
    return *this;
}

/** Comparison operator */
bool Amber2::operator==(const Amber2 &other) const
{
    return coul_14scl == other.coul_14scl and lj_14scl == other.lj_14scl;
}

/** Comparison operator */
bool Amber2::operator!=(const Amber2 &other) const
{
    return not operator==(other);
}

const char* Amber2::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Amber2>() );
}

void Amber2::set14Factors(double coul_14, double lj_14)
{
    coul_14scl = coul_14;
    lj_14scl = lj_14;
}

double Amber2::coulomb14Factor() const
{
    return coul_14scl;
}

double Amber2::lj14Factor() const
{
    return lj_14scl;
}

const char* Amber2::what() const
{
    return Amber2::typeName();
}

/** Read in the molecules from the passed Amber 7 format restart and 
    topology/parameter contents from the files, using the passed CuttingFunction to break
    molecules into parts, and the passed PropertyMap to assign data to
    molecular properties. The molecules and associated data are retruned
    in the passed SireSystem::System.
*/
SireSystem::System Amber2::readRst7Parm7(const QStringList &rstlines,
                                         const QStringList &prmlines,
                                         const SireMol::CuttingFunction &cutting_function,
                                         const PropertyMap &map) const
{


    return System();
}

/** Function that reads the entire contents of the file 'filename', returning each
    line as a QStringList */
QStringList readLines(const QString &filename)
{
    QFile f(filename);
    
    if (not f.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        throw SireError::file_error(f, CODELOC);
    }
    
    QStringList lines;
    
    QTextStream ts(&f);
    
    while (not ts.atEnd())
    {
        lines.append( ts.readLine() );
    }
    
    return lines;
}

/** Read in the molecules from the passed Amber 7 format restart and 
    topology/parameter files, using the passed CuttingFunction to break
    molecules into parts, and the passed PropertyMap to assign data to
    molecular properties. The molecules and associated data are retruned
    in the passed SireSystem::System.
*/
SireSystem::System Amber2::readRst7Parm7(const QString &rstfile,
                                         const QString &prmfile,
                                         const SireMol::CuttingFunction &cutting_function,
                                         const PropertyMap &map) const
{
    QStringList rstlines;
    QStringList prmlines;

    qDebug() << "READ";
    tbb::parallel_invoke( [&](){ rstlines = readLines(rstfile); },
                          [&](){ prmlines = readLines(prmfile); } );
    qDebug() << "READ COMPLETE" << rstlines.count() << prmlines.count();

    return this->readRst7Parm7(rstlines, prmlines, cutting_function, map);
}

/** Read in the molecules from the passed Amber 7 format restart and
    topology/parameter lines, using the passed CuttingFunction to break
    molecules into parts, and the passed PropertyMap to assign data to
    molecular properties. The molecules and associated data are retruned
    in the passed SireSystem::System.
*/
SireSystem::System Amber2::readRst7Parm7(const QStringList &rstlines,
                                         const QStringList &prmlines,
                                         const PropertyMap &map,
                                         const SireMol::CuttingFunction &cutting_function) const
{
    return this->readRst7Parm7(rstlines, prmlines, cutting_function, map);
}

/** Read in the molecules from the passed Amber 7 format restart and 
    topology/parameter files, using the passed CuttingFunction to break
    molecules into parts, and the passed PropertyMap to assign data to
    molecular properties. The molecules and associated data are retruned
    in the passed SireSystem::System.
*/
SireSystem::System Amber2::readRst7Parm7(const QString &rstfile,
                                         const QString &prmfile,
                                         const PropertyMap &map,
                                         const SireMol::CuttingFunction &cutting_function) const
{
    return this->readRst7Parm7(rstfile, prmfile, cutting_function, map);
}

/** Write the molecules in the passed system to the Amber 7 format
    restart and topology/parameter files called rstfile and prmfile,
    using the passed PropertyMap to specify which molecular properties
    should be used */
void Amber2::writeRst7Parm7(const SireSystem::System &system,
                            const QString &rstfile, const QString &prmfile,
                            const PropertyMap &map) const
{
}

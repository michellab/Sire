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
#include <QElapsedTimer>

#include "amber2.h"

#include "SireBase/parallel.h"
#include "SireBase/unittest.h"

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

//////////////
////////////// Implementation of AmberParm7
//////////////

static const RegisterMetaType<AmberParm7> r_amberparm7(NO_ROOT);

/** Serialise to a binary datastream */
QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const AmberParm7 &amberparm7)
{
    writeHeader(ds, r_amberparm7, 1);
    
    SharedDataStream sds(ds);
    
    sds << amberparm7.lnes << amberparm7.flag_to_line;
    
    return ds;
}

/** Read from a binary datastream */
QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, AmberParm7 &amberparm7)
{
    VersionID v = readHeader(ds, r_amberparm7);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> amberparm7.lnes >> amberparm7.flag_to_line;
    }
    else
        throw version_error(v, "1", r_amberparm7, CODELOC);
    
    return ds;
}

/** Constructor */
AmberParm7::AmberParm7()
{}

AmberParm7::FLAG_TYPE flagType(const QStringList &lines, const QPair<qint64,qint64> &index)
{
    //TODO TOMORROW
    #WARNING TODO TOMORROW
}

/** Return the flag type for the data associated with the passed flag.
    This returns UNKNOWN if this is not known */
AmberParm7::FLAG_TYPE AmberParm7::flagType(const QString &flag) const
{
    return ::flagType(lines, flag_to_line.value(flag));
}

/** Return the integer data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QList<qint64> AmberParm7::intData(const QString &flag) const
{
    auto it = int_data.constFind(flag);
    
    if (it != int_data.constEnd())
    {
        return it.value();
    }
    
    if (flag_to_line.contains(flag))
    {
        if (float_data.contains(flag))
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the float data for flag '%1' to integer data!")
                    .arg(flag), CODELOC );
        else if (string_data.contains(flag))
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the string data for flag '%1' to integer data!")
                    .arg(flag), CODELOC );
        else
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the data for flag '%1' to integer data!")
                    .arg(flag), CODELOC );
    }

    return QList<qint64>();
}

/** Return the float data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QList<double> AmberParm7::floatData(const QString &flag) const
{
    auto it = float_data.constFind(flag);
    
    if (it != float_data.constEnd())
    {
        return it.value();
    }
    
    if (flag_to_line.contains(flag))
    {
        if (int_data.contains(flag))
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the integer data for flag '%1' to float data!")
                    .arg(flag), CODELOC );
        else if (string_data.contains(flag))
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the string data for flag '%1' to float data!")
                    .arg(flag), CODELOC );
        else
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the data for flag '%1' to float data!")
                    .arg(flag), CODELOC );
    }

    return QList<double>();
}

/** Return the string data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QStringList AmberParm7::stringData(const QString &flag) const
{
    auto it = string_data.constFind(flag);
    
    if (it != string_data.constEnd())
    {
        return it.value();
    }
    
    if (flag_to_line.contains(flag))
    {
        if (float_data.contains(flag))
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the float data for flag '%1' to string data!")
                    .arg(flag), CODELOC );
        else if (int_data.contains(flag))
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the integer data for flag '%1' to string data!")
                    .arg(flag), CODELOC );
        else
            throw SireError::invalid_cast( QObject::tr(
                "Cannot convert the data for flag '%1' to string data!")
                    .arg(flag), CODELOC );
    }

    return QStringList();
}

QList<qint64> readIntData(const QStringList &lines, const QPair<qint64,qint64> &index)
{
    return QList<qint64>();
}

QList<double> readFloatData(const QStringList &lines, const QPair<qint64,qint64> &index)
{
    return QList<double>();
}

QStringList readStringData(const QStringList &lines, const QPair<qint64,qint64> &index)
{
    return QStringList();
}

/** Process all of the flags */
void AmberParm7::processAllFlags()
{
    QMutex int_mutex, float_mutex, string_mutex;
    
    const QStringList flags = flag_to_line.keys();
    
    tbb::parallel_for( tbb::blocked_range<int>(0,flags.count()),
                       [&](tbb::blocked_range<int> r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            const QString &flag = flags[i];
            const QPair<qint64,qint64> index = flag_to_line.value(flag);
            FLAG_TYPE flag_type = ::flagType(lines, index);
            
            switch(flag_type)
            {
                case INT:
                {
                    QList<qint64> data = readIntData(lines, index);
                    QMutexLocker lkr(&int_mutex);
                    int_data.insert(flag, data);
                    break;
                }
                case FLOAT:
                {
                    QList<double> data = readFloatData(lines, index);
                    QMutexLocker lkr(&float_mutex);
                    float_data.insert(flag, data);
                    break;
                }
                case STRING:
                {
                    QStringList data = readStringData(lines, index);
                    QMutexLocker lkr(&string_mutex);
                    string_data.insert(flag, data);
                    break;
                }
                default:
                    break;
            }
        }
    });
}

/** Construct by reading from the file called 'filename' */
AmberParm7::AmberParm7(const QString &filename)
{
    //first, open the file and read the lines
    QFile f(filename);
    
    if (not f.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        throw SireError::file_error(f, CODELOC);
    }
    
    // as we are reading, look out for any FLAGs, so that we
    // can record their locations
    QTextStream ts(&f);
    int i = 0;
    QString last_flag = QString::null;
    
    while (not ts.atEnd())
    {
        QString line = ts.readLine();
    
        lnes.append(line);
        
        if (line[0] == '%')
        {
            //this is a control line
            if (line.startsWith("%FLAG"))
            {
                //this is a new flag - close any open old flag
                if (not last_flag.isNull())
                {
                    if (flag_to_line.contains(last_flag))
                    {
                        flag_to_line[last_flag].second = i - flag_to_line[last_flag].first;
                    }
                    
                    last_flag = QString::null;
                }
                
                //find the new flag
                QStringList words = line.split(" ");
                
                QString flag = words[1];
                
                if (flag_to_line.contains(flag))
                    throw SireError::file_error( QObject::tr(
                        "The file '%1' does not look like a valid Amber Parm7 file, "
                        "as the FLAG '%2' is duplicated! (on lines %3 and %4)")
                            .arg(filename).arg(flag)
                            .arg(flag_to_line[flag].first)
                            .arg(i), CODELOC );
                
                flag_to_line.insert( flag, QPair<qint64,qint64>(i,-1) );
                last_flag = flag;
            }
        }
        
        ++i;
    }
    
    if (not last_flag.isNull())
    {
        flag_to_line[last_flag].second = i - flag_to_line[last_flag].first;
        last_flag = QString::null;
    }

    this->processAllFlags();
}

/** Construct by converting from the passed system, using the passed property
    map to find the right properties */
AmberParm7::AmberParm7(const System &system, const PropertyMap &map)
{
}

/** Copy constructor */
AmberParm7::AmberParm7(const AmberParm7 &other)
           : lnes(other.lnes), flag_to_line(other.flag_to_line)
{}

/** Destructor */
AmberParm7::~AmberParm7()
{}

/** Copy assignment operator */
AmberParm7& AmberParm7::operator=(const AmberParm7 &other)
{
    if (this != &other)
    {
        lnes = other.lnes;
        flag_to_line = other.flag_to_line;
    }
    
    return *this;
}

/** Comparison operator */
bool AmberParm7::operator==(const AmberParm7 &other) const
{
    return lnes == other.lnes;
}

/** Comparison operator */
bool AmberParm7::operator!=(const AmberParm7 &other) const
{
    return not operator==(other);
}

const char* AmberParm7::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberParm7>() );
}

const char* AmberParm7::what() const
{
    return AmberParm7::typeName();
}

/** Return an AmberParm7 object read from the passed file */
AmberParm7 AmberParm7::read(const QString &filename)
{
    return AmberParm7(filename);
}

/** Return an AmberParm7 object created from the passed System */
AmberParm7 AmberParm7::write(const System &system, const PropertyMap &map)
{
    return AmberParm7(system, map);
}

/** Return the raw lines of the Parm7 file */
QStringList AmberParm7::lines() const
{
    return lnes;
}

/** Return the lines that correspond to the passed flag. This returns an
    empty list of there are no lines associated with the passed flag */
QStringList AmberParm7::lines(const QString &flag) const
{
    auto it = flag_to_line.constFind(flag);
    
    if (it != flag_to_line.constEnd())
    {
        const int start = it->first;
        const int count = it->second;
        
        qDebug() << start << count;
        
        SireBase::assert_true( start >= 0 and start < lnes.count(), CODELOC );
        SireBase::assert_true( count > 0 and start+count < lnes.count(), CODELOC );

        return lnes.mid(start,count);
    }
    else
        return QStringList();
}

/** Return all of the flags that are held in this file */
QStringList AmberParm7::flags() const
{
    return flag_to_line.keys();
}

//////////////
////////////// Implementation of Amber2
//////////////

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
System Amber2::readRst7Parm7(const QString &rstfile,
                             const QString &prmfile,
                             const CuttingFunction &cutting_function,
                             const PropertyMap &map) const
{
    QStringList rstlines;
    AmberParm7 prm7;

    qDebug() << "READ";

    QElapsedTimer t;
    t.start();

    tbb::parallel_invoke( [&](){ rstlines = readLines(rstfile); },
                          [&](){ prm7 = AmberParm7::read(prmfile); } );

    qint64 ns = t.nsecsElapsed();

    qDebug() << "READ COMPLETE" << rstlines.count() << prm7.lines().count();
    qDebug() << "TOOK" << (0.000001*ns) << "ms";

    return System();
}

/** Read in the molecules from the passed Amber 7 format restart and 
    topology/parameter files, using the passed CuttingFunction to break
    molecules into parts, and the passed PropertyMap to assign data to
    molecular properties. The molecules and associated data are retruned
    in the passed SireSystem::System.
*/
System Amber2::readRst7Parm7(const QString &rstfile,
                             const QString &prmfile,
                             const PropertyMap &map,
                             const CuttingFunction &cutting_function) const
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

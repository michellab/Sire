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
#include <QRegularExpression>

#include "amber2.h"

#include "SireMol/molecule.h"

#include "SireBase/parallel.h"
#include "SireBase/unittest.h"

#include "SireSystem/system.h"
#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

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
////////////// Implementation of MoleculeParser
//////////////

static const RegisterMetaType<MoleculeParser> r_parser( MAGIC_ONLY, MoleculeParser::typeName() );

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const MoleculeParser &parser)
{
    writeHeader(ds, r_parser, 1);
    ds << static_cast<const Property&>(parser);
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, MoleculeParser &parser)
{
    VersionID v = readHeader(ds, r_parser);
    
    if (v == 1)
    {
        ds >> static_cast<Property&>(parser);
    }
    else
        throw version_error(v,"1", r_parser, CODELOC);
    
    return ds;
}

/** Constructor */
MoleculeParser::MoleculeParser() : Property()
{}

/** Copy constructor */
MoleculeParser::MoleculeParser(const MoleculeParser &other)
               : Property(other)
{}

/** Destructor */
MoleculeParser::~MoleculeParser()
{}

const char* MoleculeParser::typeName()
{
    return "SireIO::MoleculeParser";
}

MoleculeParser& MoleculeParser::operator=(const MoleculeParser &other)
{
    Property::operator=(other);
    return *this;
}

bool MoleculeParser::operator==(const MoleculeParser &other) const
{
    return Property::operator==(other);
}

bool MoleculeParser::operator!=(const MoleculeParser &other) const
{
    return Property::operator!=(other);
}

//////////////
////////////// Implementation of AmberRst
//////////////

static const RegisterMetaType<AmberRst> r_rst;

QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const AmberRst &rst)
{
    writeHeader(ds, r_rst, 1);
    
    SharedDataStream sds(ds);
    
    sds << static_cast<const MoleculeParser&>(rst);
    
    return ds;
}

QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, AmberRst &rst)
{
    VersionID v = readHeader(ds, r_rst);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> static_cast<MoleculeParser&>(rst);
    }
    else
        throw version_error(v, "1", r_rst, CODELOC);
    
    return ds;
}

/** Constructor */
AmberRst::AmberRst() : ConcreteProperty<AmberRst,MoleculeParser>()
{}

/** Copy constructor */
AmberRst::AmberRst(const AmberRst &other)
         : ConcreteProperty<AmberRst,MoleculeParser>(other)
{}

/** Destructor */
AmberRst::~AmberRst()
{}

AmberRst& AmberRst::operator=(const AmberRst &other)
{
    if (this != &other)
    {
        MoleculeParser::operator=(other);
    }

    return *this;
}

bool AmberRst::operator==(const AmberRst &other) const
{
    return MoleculeParser::operator==(other);
}

bool AmberRst::operator!=(const AmberRst &other) const
{
    return MoleculeParser::operator!=(other);
}

const char* AmberRst::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberRst>() );
}

const char* AmberRst::what() const
{
    return AmberRst::typeName();
}

QString AmberRst::toString() const
{
    return QObject::tr("AmberRst::null");
}

//////////////
////////////// Implementation of AmberParm
//////////////

static const RegisterMetaType<AmberParm> r_parm;

/** Serialise to a binary datastream */
QDataStream SIREIO_EXPORT &operator<<(QDataStream &ds, const AmberParm &parm)
{
    writeHeader(ds, r_parm, 1);
    
    SharedDataStream sds(ds);
    
    sds << parm.lnes << parm.flag_to_line
        << parm.int_data << parm.float_data << parm.string_data
        << static_cast<const MoleculeParser&>(parm);
    
    return ds;
}

/** Function called after loading the AmberParm from a binary stream
    to populate all of the calculated member data */
void AmberParm::rebuildAfterReload()
{
    pointers = this->intData("POINTERS");
    
    if (pointers.isEmpty())
    {
        this->operator=( AmberParm() );
    }
}

/** Read from a binary datastream */
QDataStream SIREIO_EXPORT &operator>>(QDataStream &ds, AmberParm &parm)
{
    VersionID v = readHeader(ds, r_parm);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        parm = AmberParm();
        
        sds >> parm.lnes >> parm.flag_to_line
            >> parm.int_data >> parm.float_data >> parm.string_data
            >> static_cast<MoleculeParser&>(parm);
        
        parm.rebuildAfterReload();
    }
    else
        throw version_error(v, "1", r_parm, CODELOC);
    
    return ds;
}

/** Constructor */
AmberParm::AmberParm() : ConcreteProperty<AmberParm,MoleculeParser>()
{
    for (int i=0; i<32; ++i)
    {
        pointers.append(0);
    }
}

/** Internal class used by AmberParm to hold format description */
class AmberFormat
{
public:
    AmberParm::FLAG_TYPE flag_type;
    int num_values;
    int field_width;
    int point_width;

    AmberFormat() : flag_type(AmberParm::UNKNOWN),
                    num_values(0), field_width(0), point_width(0)
    {}
    
    AmberFormat(const QString &line)
    {
        QRegularExpression re("%FORMAT\\((\\d+)(\\w)(\\d+)\\.?(\\d+)?\\)");
        
        auto m = re.match(line);
        
        if (not m.hasMatch())
        {
            throw SireError::io_error( QObject::tr(
                    "Could not extract the format from the line '%1'. "
                    "Expected to read a line similar to '%FORMAT(5E16.8)'.")
                        .arg(line), CODELOC );
        }
        
        num_values = m.captured(1).toInt();
        field_width = m.captured(3).toInt();
        
        QString typ = m.captured(2).toLower();
        
        if (typ == "a")
            flag_type = AmberParm::STRING;
        else if (typ == "i")
            flag_type = AmberParm::INTEGER;
        else if (typ == "e")
            flag_type = AmberParm::FLOAT;
        else
            flag_type = AmberParm::UNKNOWN;
        
        if (not m.captured(4).isNull())
        {
            point_width = m.captured(4).toInt();
        }
    }
    
    ~AmberFormat()
    {}
    
    QString toString() const
    {
        switch(flag_type)
        {
            case AmberParm::STRING:
                return QString("AmberFormat( %1 x string[width = %2] )")
                            .arg(num_values).arg(field_width);
            case AmberParm::INTEGER:
                return QString("AmberFormat( %1 x integer[width = %2] )")
                            .arg(num_values).arg(field_width);
            case AmberParm::FLOAT:
                return QString("AmberFormat( %1 x float[width = %2, precision = %3] )")
                            .arg(num_values).arg(field_width).arg(point_width);
            default:
                return QString("AmberFormat( UNKNOWN )");
        }
    }
    
    /** The width of each field (number of columns) */
    int width() const
    {
        return field_width;
    }
    
    /** The maximum number of items per line */
    int numValues() const
    {
        return num_values;
    }
    
    /** The number of values after the decimal point for float values */
    int pointWidth() const
    {
        return point_width;
    }
    
    /** Return the number of values to read from the passed line */
    int numValues( const QString &line ) const
    {
        return qMin( num_values, line.length() / field_width );
    }
};

AmberParm::FLAG_TYPE flagType(const QStringList &lines, const QPair<qint64,qint64> &index)
{
    AmberFormat f(lines[index.first+1]);
    return f.flag_type;
}

/** Return the flag type for the data associated with the passed flag.
    This returns UNKNOWN if this is not known */
AmberParm::FLAG_TYPE AmberParm::flagType(const QString &flag) const
{
    if (flag_to_line.contains(flag))
    {
        auto index = flag_to_line.value(flag);
    
        if (index.first > 0)
            return AmberFormat(lnes[index.first-1]).flag_type;
    }

    return AmberParm::UNKNOWN;
}

/** Return the integer data for the passed flag. This returns an empty
    list if there is no data associated with this flag. This raises
    an invalid_cast error if data exists, but it is the wrong type */
QList<qint64> AmberParm::intData(const QString &flag) const
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
QList<double> AmberParm::floatData(const QString &flag) const
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
QStringList AmberParm::stringData(const QString &flag) const
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

QList<qint64> readIntData(const QStringList &lines, AmberFormat format,
                          const QPair<qint64,qint64> &index,
                          QStringList *errors=0)
{
    QList<qint64> data;

    //read in all of the lines...
    const int strt = index.first;
    const int end = index.first + index.second;
    
    for (int i=strt; i<end; ++i)
    {
        const QString &line = lines[i];
        
        int nvalues = format.numValues(line);
        
        if (errors and nvalues < format.numValues() and i != end-1)
        {
            //one of the data lines has too little data
            errors->append( QObject::tr("Too few data values on line %1: "
                "Expected %2 values but only saw %3.")
                    .arg(i+1).arg(format.numValues()).arg(nvalues) );
        }
        
        for (int j=0; j<nvalues; ++j)
        {
            auto ref = line.midRef( j*format.width(), format.width() );
            
            if (not ref.isNull())
            {
                bool ok = true;
                qint64 value = ref.toLong(&ok);
                
                if (ok)
                    data.append(value);
                else if (errors)
                    errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 (%3) into an integer!")
                            .arg(i+1).arg(j+1).arg(ref.toString()) );
            }
            else if (errors)
            {
                errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 as the string is null!")
                            .arg(i+1).arg(j+1) );
            }
        }
    }

    return data;
}

QList<double> readFloatData(const QStringList &lines, AmberFormat format,
                            const QPair<qint64,qint64> &index,
                            QStringList *errors)
{
    QList<double> data;

    //read in all of the lines...
    const int strt = index.first;
    const int end = index.first + index.second;
    
    for (int i=strt; i<end; ++i)
    {
        const QString &line = lines[i];
        
        int nvalues = format.numValues(line);
        
        if (errors and nvalues < format.numValues() and i != end-1)
        {
            //one of the data lines has too little data
            errors->append( QObject::tr("Too few data values on line %1: "
                "Expected %2 values but only saw %3.")
                    .arg(i+1).arg(format.numValues()).arg(nvalues) );
        }
        
        for (int j=0; j<nvalues; ++j)
        {
            auto ref = line.midRef( j*format.width(), format.width() );
            
            if (not ref.isNull())
            {
                bool ok = true;
                double value = ref.toDouble(&ok);
                
                if (ok)
                    data.append(value);
                else if (errors)
                    errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 (%3) into an integer!")
                            .arg(i+1).arg(j+1).arg(ref.toString()) );
            }
            else if (errors)
            {
                errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 as the string is null!")
                            .arg(i+1).arg(j+1) );
            }
        }
    }

    return data;
}

QStringList readStringData(const QStringList &lines, AmberFormat format,
                           const QPair<qint64,qint64> &index,
                           QStringList *errors)
{
    QStringList data;

    //read in all of the lines...
    const int strt = index.first;
    const int end = index.first + index.second;
    
    for (int i=strt; i<end; ++i)
    {
        const QString &line = lines[i];
        
        int nvalues = format.numValues(line);
        
        if (errors and nvalues < format.numValues() and i != end-1)
        {
            //one of the data lines has too little data
            errors->append( QObject::tr("Too few data values on line %1: "
                "Expected %2 values but only saw %3.")
                    .arg(i+1).arg(format.numValues()).arg(nvalues) );
        }
        
        for (int j=0; j<nvalues; ++j)
        {
            auto ref = line.midRef( j*format.width(), format.width() );
            
            if (not ref.isNull())
            {
                data.append( ref.toString() );
            }
            else if (errors)
            {
                errors->append( QObject::tr("Failed to convert the data "
                       "on line %1, column %2 as the string is null!")
                            .arg(i+1).arg(j+1) );
            }
        }
    }

    return data;
}

/** Process all of the flags */
void AmberParm::processAllFlags()
{
    QMutex int_mutex, float_mutex, string_mutex;
    
    const QStringList flags = flag_to_line.keys();
    
    QStringList global_errors;
    
    tbb::parallel_for( tbb::blocked_range<int>(0,flags.count()),
                       [&](tbb::blocked_range<int> r)
    {
        QStringList local_errors;
    
        for (int i=r.begin(); i<r.end(); ++i)
        {
            const QString &flag = flags[i];
            const QPair<qint64,qint64> index = flag_to_line.value(flag);
            
            //the format for the data is on the preceeding line
            const AmberFormat format(lnes[index.first-1]);
            
            switch(format.flag_type)
            {
                case INTEGER:
                {
                    QList<qint64> data = readIntData(lnes, format, index, &local_errors);
                    QMutexLocker lkr(&int_mutex);
                    int_data.insert(flag, data);
                    break;
                }
                case FLOAT:
                {
                    QList<double> data = readFloatData(lnes, format, index, &local_errors);
                    QMutexLocker lkr(&float_mutex);
                    float_data.insert(flag, data);
                    break;
                }
                case STRING:
                {
                    QStringList data = readStringData(lnes, format, index, &local_errors);
                    QMutexLocker lkr(&string_mutex);
                    string_data.insert(flag, data);
                    break;
                }
                default:
                    break;
            }
        }
        
        if (not local_errors.isEmpty())
        {
            QMutexLocker lkr(&int_mutex);
            global_errors += local_errors;
        }
    });
    
    if (not global_errors.isEmpty())
    {
        throw SireError::io_error( QObject::tr(
            "There were errors parsing the file:%1")
                .arg(global_errors.join("\n")), CODELOC );
    }
    
    pointers = int_data.value("POINTERS");
    
    if (pointers.count() < 31)
    {
        throw SireError::io_error( QObject::tr(
            "There was no, or an insufficient 'POINTERS' section in the file! (%1)")
                .arg(pointers.count()), CODELOC );
    }
}

/** Construct by reading from the file called 'filename' */
AmberParm::AmberParm(const QString &filename) : ConcreteProperty<AmberParm,MoleculeParser>()
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
                
                //skip the FLAG line, and the FORMAT line that must come immediately after
                flag_to_line.insert( flag, QPair<qint64,qint64>(i+2,-1) );
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

    //now process all of the flag data
    this->processAllFlags();
    
    //finally, make sure that we have been constructed sane
    this->assertSane();
}

/** Construct by converting from the passed system, using the passed property
    map to find the right properties */
AmberParm::AmberParm(const System &system, const PropertyMap &map)
{
    this->operator=(AmberParm());
}

/** Copy constructor */
AmberParm::AmberParm(const AmberParm &other)
           : ConcreteProperty<AmberParm,MoleculeParser>(other),
             lnes(other.lnes), flag_to_line(other.flag_to_line),
             int_data(other.int_data), float_data(other.float_data),
             string_data(other.string_data), pointers(other.pointers)
{}

/** Destructor */
AmberParm::~AmberParm()
{}

/** Copy assignment operator */
AmberParm& AmberParm::operator=(const AmberParm &other)
{
    if (this != &other)
    {
        lnes = other.lnes;
        flag_to_line = other.flag_to_line;
        int_data = other.int_data;
        float_data = other.float_data;
        string_data = other.string_data;
        pointers = other.pointers;
        MoleculeParser::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool AmberParm::operator==(const AmberParm &other) const
{
    return lnes == other.lnes and MoleculeParser::operator==(other);
}

/** Comparison operator */
bool AmberParm::operator!=(const AmberParm &other) const
{
    return not operator==(other);
}

const char* AmberParm::typeName()
{
    return QMetaType::typeName( qMetaTypeId<AmberParm>() );
}

const char* AmberParm::what() const
{
    return AmberParm::typeName();
}

/** Return a string representation of this object */
QString AmberParm::toString() const
{
    if (lnes.isEmpty())
        return QObject::tr("AmberParm::null");
    else
    {
        return QObject::tr("AmberParm( title() = '%1', nMolecules() = %2, "
                           "nResidues() = %3, nAtoms() = %4 )")
                    .arg(title()).arg(nMolecules()).arg(nResidues()).arg(nAtoms());
    }
}

/** Return the title of the parameter file */
QString AmberParm::title() const
{
    return string_data.value("TITLE").join("");
}

/** Return the total number of atoms in the file */
int AmberParm::nAtoms() const
{
    return pointers[0];
}

/** Return the total number of atoms in the ith molecule in the file */
int AmberParm::nAtoms(int i) const
{
    return 0;
}

/** Return the number of distinct atom types */
int AmberParm::nTypes() const
{
    return pointers[1];
}

/** Return the total number of bonds */
int AmberParm::nBonds() const
{
    return nBondsWithHydrogen() + nBondsNoHydrogen();
}

/** Return the number of bonds including hydrogen */
int AmberParm::nBondsWithHydrogen() const
{
    return pointers[2];
}

/** Return the number of bonds no containing hydrogen */
int AmberParm::nBondsNoHydrogen() const
{
    return pointers[3];
}

/** Return the number of angles */
int AmberParm::nAngles() const
{
    return nAnglesWithHydrogen() + nAnglesNoHydrogen();
}

/** Return the number of angles containing hydrogen */
int AmberParm::nAnglesWithHydrogen() const
{
    return pointers[4];
}

/** Return the number of angles without hydrogen */
int AmberParm::nAnglesNoHydrogen() const
{
    return pointers[5];
}

/** Return the number of dihedrals */
int AmberParm::nDihedrals() const
{
    return nDihedralsWithHydrogen() + nDihedralsNoHydrogen();
}

/** Return the number of dihedrals containing hydrogen */
int AmberParm::nDihedralsWithHydrogen() const
{
    return pointers[6];
}

/** Return the number of dihedrals without hydrogen */
int AmberParm::nDihedralsNoHydrogen() const
{
    return pointers[7];
}

/** Return the number of excluded atoms */
int AmberParm::nExcluded() const
{
    return pointers[9];
}

/** Return the number of residues */
int AmberParm::nResidues() const
{
    return pointers[10];
}

/** Return the number of molecules in the file */
int AmberParm::nMolecules() const
{
    QList<qint64> atoms_per_mol = int_data.value("ATOMS_PER_MOLECULE");
    
    if (atoms_per_mol.isEmpty())
    {
        if (nAtoms() == 0)
            return 0;
        else
            return 1;
    }
    else
        return atoms_per_mol.count();
}

/** Return the first index of the atom in each molecule, together with
    the number of atoms in that molecule */
QList< QPair<int,int> > AmberParm::moleculeIndicies() const
{
    QList< QPair<int,int> > idxs;
    
    QList<qint64> atoms_per_mol = int_data.value("ATOMS_PER_MOLECULE");
    
    if (atoms_per_mol.isEmpty())
    {
        if (nAtoms() > 0)
        {
            idxs.append( QPair<int,int>(0,nAtoms()) );
        }
    }
    else
    {
        qint64 nats = 0;
    
        for (auto a : atoms_per_mol)
        {
            idxs.append( QPair<int,int>(nats, a) );
            nats += a;
        }
    }
    
    return idxs;
}

/** Return an AmberParm object read from the passed file */
AmberParm AmberParm::read(const QString &filename)
{
    return AmberParm(filename);
}

/** Return an AmberParm object created from the passed System */
AmberParm AmberParm::write(const System &system, const PropertyMap &map)
{
    return AmberParm(system, map);
}

/** Return the raw lines of the Parm7 file */
QStringList AmberParm::lines() const
{
    return lnes;
}

/** Return the lines that correspond to the passed flag. This returns an
    empty list of there are no lines associated with the passed flag */
QStringList AmberParm::lines(const QString &flag) const
{
    auto it = flag_to_line.constFind(flag);
    
    if (it != flag_to_line.constEnd())
    {
        const int start = it->first;
        const int count = it->second;
        
        SireBase::assert_true( start >= 0 and start < lnes.count(), CODELOC );
        SireBase::assert_true( count > 0 and start+count < lnes.count(), CODELOC );

        return lnes.mid(start,count);
    }
    else
        return QStringList();
}

/** Return all of the flags that are held in this file */
QStringList AmberParm::flags() const
{
    return flag_to_line.keys();
}

template<class FUNC, class Exception>
void live_test(FUNC function, QList<boost::shared_ptr<Exception>> &errors)
{
    try
    {
        return function();
    }
    catch(const Exception &e)
    {
        errors.append( boost::shared_ptr<Exception>(e.clone()) );
    }
}

/** Run through all of the data that has been read and perform a series
    of tests that will see if the prm7 data is sane. If any test fails,
    then an exception will be thrown */
void AmberParm::assertSane() const
{
    QList<boost::shared_ptr<SireError::exception>> errors;
    
    int natoms = this->nAtoms();
    
    live_test( [&](){ SireBase::assert_equal(this->stringData("ATOM_NAME").count(),
                                             natoms, CODELOC);}, errors );
    
    live_test( [&](){ SireBase::assert_equal(this->floatData("CHARGE").count(),
                                             natoms, CODELOC);}, errors );
    
    live_test( [&](){ SireBase::assert_equal(this->intData("ATOMIC_NUMBER").count(),
                                             natoms, CODELOC);}, errors );
    
    live_test( [&](){ SireBase::assert_equal(this->floatData("MASS").count(),
                                             natoms, CODELOC);}, errors );
    
    live_test( [&](){ SireBase::assert_equal(this->intData("ATOM_TYPE_INDEX").count(),
                                             natoms, CODELOC);}, errors );
    
    if (not errors.isEmpty())
    {
        for (auto error : errors)
        {
            qDebug() << error->toString();
        }
    
        throw SireError::io_error( QObject::tr(
            "Sanity tests failed for the loaded Amber prm7 format file"), CODELOC );
    }
}

/** Internal function used to get the molecule that starts at index 'start_idx'
    in the file, and that has 'natoms' atoms */
Molecule AmberParm::getMolecule(int start_idx, int natoms) const
{
    return Molecule();
}

/** Return the ith molecule that is described by this AmberParm file. Note
    that this molecule won't have any coordinate data, as this is not
    provided in this file */
Molecule AmberParm::getMolecule(int idx) const
{
    const QList< QPair<int,int> > mol_idxs = this->moleculeIndicies();
    idx = Index(idx).map(mol_idxs.count());
    return this->getMolecule( mol_idxs[idx].first, mol_idxs[idx].second );
}

/** Return the System that is described by this AmberParm file. Note that
    the molecules in this system don't have any coordinates (as these aren't
    provided by the file */
System AmberParm::toSystem() const
{
    const QList< QPair<int,int> > mol_idxs = this->moleculeIndicies();

    const int nmols = mol_idxs.count();
    
    if (nmols == 0)
        return System();
    
    QVector<Molecule> global_mols(nmols);
    QMutex global_mutex;
    
    tbb::parallel_for( tbb::blocked_range<int>(0,nmols),
                       [&](tbb::blocked_range<int> r)
    {
        QVector<Molecule> local_mols( r.end()-r.begin() );
        
        //create and populate all of the molecules
        for (int i=r.begin(); i<r.end(); ++i)
        {
            local_mols[i-r.begin()] = this->getMolecule( mol_idxs[i].first,
                                                         mol_idxs[i].second );
        }
        
        //copy them over to global memory
        QMutexLocker lkr(&global_mutex);
        for (int i=r.begin(); i<r.end(); ++i)
        {
            global_mols[i] = local_mols[i-r.begin()];
        }
    });
    
    MoleculeGroup mols("all");
    
    System system( this->title() );
    system.add(mols);
    
    return system;
}

//////////////
////////////// Implementation of Amber2
//////////////

static const RegisterMetaType<Amber2> r_amber2;

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
Amber2::Amber2()
       : ConcreteProperty<Amber2,MoleculeParser>(),
         coul_14scl(AMBER14COUL), lj_14scl(AMBER14LJ)
{}

/** Copy constructor */
Amber2::Amber2(const Amber2 &other)
       : ConcreteProperty<Amber2,MoleculeParser>(other),
         coul_14scl(other.coul_14scl), lj_14scl(other.lj_14scl)
{}

/** Destructor */
Amber2::~Amber2()
{}

/** Copy assignment operator */
Amber2& Amber2::operator=(const Amber2 &other)
{
    if (this != &other)
    {
        coul_14scl = other.coul_14scl;
        lj_14scl = other.lj_14scl;
        
        MoleculeParser::operator=(other);
    }
        
    return *this;
}

/** Comparison operator */
bool Amber2::operator==(const Amber2 &other) const
{
    return coul_14scl == other.coul_14scl and
           lj_14scl == other.lj_14scl and
           MoleculeParser::operator==(other);
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

/** Return a string representation of this object */
QString Amber2::toString() const
{
    return QObject::tr("Amber2( coulomb14Factor() = %1, lj14Factor() = %2 )")
                .arg(coulomb14Factor()).arg(lj14Factor());
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
System Amber2::readRstParm(const QString &rstfile,
                           const QString &prmfile,
                           const CuttingFunction &cutting_function,
                           const PropertyMap &map) const
{
    QStringList rstlines;
    AmberParm prm;

    qDebug() << "READ";

    QElapsedTimer t;
    t.start();

    tbb::parallel_invoke( [&](){ rstlines = readLines(rstfile); },
                          [&](){ prm = AmberParm::read(prmfile); } );

    qint64 ns = t.nsecsElapsed();

    qDebug() << "READ COMPLETE" << rstlines.count() << prm.lines().count();
    qDebug() << "TOOK" << (0.000001*ns) << "ms";

    return System();
}

/** Read in the molecules from the passed Amber 7 format restart and 
    topology/parameter files, using the passed CuttingFunction to break
    molecules into parts, and the passed PropertyMap to assign data to
    molecular properties. The molecules and associated data are retruned
    in the passed SireSystem::System.
*/
System Amber2::readRstParm(const QString &rstfile,
                           const QString &prmfile,
                           const PropertyMap &map,
                           const CuttingFunction &cutting_function) const
{
    return this->readRstParm(rstfile, prmfile, cutting_function, map);
}

/** Write the molecules in the passed system to the Amber 7 format
    restart and topology/parameter files called rstfile and prmfile,
    using the passed PropertyMap to specify which molecular properties
    should be used */
void Amber2::writeRstParm(const SireSystem::System &system,
                          const QString &rstfile, const QString &prmfile,
                          const PropertyMap &map) const
{
}

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

#ifndef SIREMATHS_HISTOGRAM_H
#define SIREMATHS_HISTOGRAM_H

#include <QVector>

#include "SireUnits/dimensions.h"
#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{
class HistogramBin;
class HistogramValue;
class Histogram;
}

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::HistogramBin&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::HistogramBin&);

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::HistogramValue&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::HistogramValue&);

SIREMATHS_EXPORT QDataStream& operator<<(QDataStream&, const SireMaths::Histogram&);
SIREMATHS_EXPORT QDataStream& operator>>(QDataStream&, SireMaths::Histogram&);

namespace SireMaths
{

/** This class represents a single histogram bin */
class SIREMATHS_EXPORT HistogramBin
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const HistogramBin&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, HistogramBin&);

public:
    HistogramBin();
    HistogramBin(double minval, double maxval);
    
    HistogramBin(const HistogramBin &other);
    
    HistogramBin& operator=(const HistogramBin &other);
    
    bool operator==(const HistogramBin &other) const;
    bool operator!=(const HistogramBin &other) const;
    
    double minimum() const;
    double middle() const;
    double maximum() const;

    QString toString() const;

private:
    /** The minimum value of the bin */
    double minval;
    
    /** The maximum value of the bin */
    double maxval;
};

/** This class represents a single histogram bin with its associated value */
class SIREMATHS_EXPORT HistogramValue : public HistogramBin
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const HistogramValue&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, HistogramValue&);

public:
    HistogramValue();
    HistogramValue(const HistogramBin &bin, double value);
    
    HistogramValue(const HistogramValue &other);
    
    ~HistogramValue();
    
    HistogramValue& operator=(const HistogramValue &other);
    
    bool operator==(const HistogramValue &other) const;
    bool operator!=(const HistogramValue &other) const;
    
    double value() const;

    QString toString() const;

private:
    /** The actual value in the bin */
    double val;
};

/** This class holds a simple one-dimensional (x,y) histogram of values.

    @author Christopher Woods
*/
class SIREMATHS_EXPORT Histogram
        : public SireBase::ConcreteProperty<Histogram,SireBase::Property>
{

friend SIREMATHS_EXPORT QDataStream& ::operator<<(QDataStream&, const Histogram&);
friend SIREMATHS_EXPORT QDataStream& ::operator>>(QDataStream&, Histogram&);

public:
    Histogram();
    
    Histogram(double binwidth);
    Histogram(double binwidth, const QVector<double> &values);

    Histogram(const Histogram &other);
    
    ~Histogram();
    
    Histogram& operator=(const Histogram &other);
    
    static const char* typeName();
    
    const char* what() const;
    
    Histogram* clone() const;

    bool operator==(const Histogram &other) const;
    bool operator!=(const Histogram &other) const;
    
    Histogram& operator+=(const Histogram &other);
    Histogram& operator+=(double value);
    Histogram& operator+=(const QVector<double> &values);
    
    Histogram operator+(const Histogram &other) const;
    Histogram operator+(double value) const;
    Histogram operator+(const QVector<double> &values) const;
    
    QString toString() const;
    
    HistogramValue operator[](int i) const;

    HistogramValue at(int i) const;
    
    int count() const;
    int size() const;
    
    QVector<HistogramValue> values() const;
    QVector<HistogramValue> normalDistribution() const;
    
    void accumulate(double value);
    void accumulate(double value, double weight);
    
    void accumulate(const QVector<double> &values);
    
    void accumulate(const Histogram &other);

    void add(double value);
    void add(double value, double weight);
    
    void add(const QVector<double> &values);
    
    void add(const Histogram &other);

    double sumOfBins() const;
    
    double binWidth() const;
    
    double mean() const;
    double median() const;
    double mode() const;

    double skew() const;
    double kirtosis() const;

    double meanOfSquares() const;

    double range() const;

    double standardDeviation() const;

    static double tValue(int nsamples, double level);
    double tValue(double level) const;

    double standardError() const;
    double standardError(double level) const;
    
    double maximumValue() const;
    double minimumValue() const;

    Histogram normalise() const;

    Histogram resize(double binwidth) const;

private:
    /** The values in each of the bins */
    QHash<qint64,double> binvals;
    
    /** The bin width */
    double binwidth;
    
    /** The mean value */
    double avgval;
    
    /** The mean square value */
    double avgval2;
    
    /** The sum of weights */
    double sum_of_bins;
};

}

Q_DECLARE_METATYPE( SireMaths::Histogram )

SIRE_EXPOSE_CLASS( SireMaths::HistogramBin )
SIRE_EXPOSE_CLASS( SireMaths::HistogramValue )
SIRE_EXPOSE_CLASS( SireMaths::Histogram )

SIRE_END_HEADER

#endif

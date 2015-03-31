/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
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

#include <limits>
#include <cmath>

#ifdef SIRE_USE_SSE
   #ifdef __SSE__
       #include <emmintrin.h>   // CONDITIONAL_INCLUDE
   #else
       #undef SIRE_USE_SSE
   #endif
#endif

#include "cartesian.h"
#include "coordgroup.h"

#include "SireBase/countflops.h"

#include "SireMaths/rangenerator.h"

#include "SireError/errors.h"
#include "SireVol/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>

using namespace SireVol;
using namespace SireBase;
using namespace SireUnits::Dimension;
using namespace SireStream;

using boost::tuple;

static const RegisterMetaType<Cartesian> r_cartesian;

/** Serialise to a binary datastream */
QDataStream SIREVOL_EXPORT &operator<<(QDataStream &ds, const Cartesian &cart)
{
    writeHeader(ds, r_cartesian, 1)
                 << static_cast<const Space&>(cart);

    return ds;
}

/** Deserialise from a binary datastream */
QDataStream SIREVOL_EXPORT &operator>>(QDataStream &ds, Cartesian &cart)
{
    VersionID v = readHeader(ds, r_cartesian);

    if (v == 1)
    {
        ds >> static_cast<Space&>(cart);
    }
    else
        throw version_error(v, "1", r_cartesian, CODELOC);

    return ds;
}

/** Construct a default Cartesian volume */
Cartesian::Cartesian() : ConcreteProperty<Cartesian,Space>()
{}

/** Copy constructor */
Cartesian::Cartesian(const Cartesian &other) 
          : ConcreteProperty<Cartesian,Space>(other)
{}

/** Destructor */
Cartesian::~Cartesian()
{}

/** Copy assignment operator */
Cartesian& Cartesian::operator=(const Cartesian &other)
{
    Space::operator=(other);
    return *this;
}

/** Comparison operator */
bool Cartesian::operator==(const Cartesian &other) const
{
    return other.what() == Cartesian::typeName();
}

/** Comparison operator */
bool Cartesian::operator!=(const Cartesian &other) const
{
    return other.what() != Cartesian::typeName();
}

/** Return a string representation of this space */
QString Cartesian::toString() const
{
    return QObject::tr("Infinite cartesian space");
}

/** A Cartesian space is not periodic */
bool Cartesian::isPeriodic() const
{
    return false;
}

/** A Cartesian space is cartesian */
bool Cartesian::isCartesian() const
{
    return true;
}

/** Throw an exception as an infinite space doesn't have a volume! */
SireUnits::Dimension::Volume Cartesian::volume() const
{
    throw SireVol::incompatible_space( QObject::tr(
        "It is not possible to calculate the volume of an infinite space!"),
            CODELOC );

    return SireUnits::Dimension::Volume(0);
}

/** Throw an exception as an infinite space doesn't have a volume! */
SpacePtr Cartesian::setVolume(SireUnits::Dimension::Volume) const
{
    throw SireVol::incompatible_space( QObject::tr(
        "It is not possible to change the volume of an infinite space!"),
            CODELOC );

    return SpacePtr();
}

/** Calculate the distance between two points */
double Cartesian::calcDist(const Vector &point0, const Vector &point1) const
{
    return Vector::distance(point0, point1);
}

/** Calculate the distance squared between two points */
double Cartesian::calcDist2(const Vector &point0, const Vector &point1) const
{
    return Vector::distance2(point0, point1);
}

/** Populate the matrix 'mat' with the distances between all points in
    the group 'group'. Return the shortest distance between points. */
double Cartesian::calcDist(const CoordGroup &group, DistMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());
    double tmpdist;

    int n = group.count();
    const Vector *array = group.constData();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n,n);

    for (int i=0; i<n; ++i)
    {
        const Vector &point = array[i];
        mat.setOuterIndex(i);

        //zero diagonal
        mat[i] = 0;

        for (int j=i+1; j<n; ++j)
        {
            //calculate the distance between the two points
            tmpdist = Vector::distance(point,array[j]);

            //store the minimum distance
            mindist = qMin(tmpdist,mindist);

            //place this distance into the matrix
            mat[j] = tmpdist;
            mat(j,i) = tmpdist;
        }
    }

    return mindist;
}

/** Populate the matrix 'mat' with the distances^2 between all points in
    the group 'group'. Return the shortest distance between points. */
double Cartesian::calcDist2(const CoordGroup &group, DistMatrix &mat) const
{
    double mindist2(std::numeric_limits<double>::max());
    double tmpdist2;

    int n = group.count();
    const Vector *array = group.constData();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n,n);

    for (int i=0; i<n; ++i)
    {
        const Vector &point = array[i];
        mat.setOuterIndex(i);

        //zero diagonal
        mat[i] = 0;

        for (int j=i+1; j<n; ++j)
        {
            //calculate the distance between the two points
            tmpdist2 = Vector::distance2(point,array[j]);

            //store the minimum distance
            mindist2 = qMin(tmpdist2,mindist2);

            //place this distance into the matrix
            mat[j] = tmpdist2;
            mat(j,i) = tmpdist2;
        }
    }

    return sqrt(mindist2);
}

/** Populate the matrix 'mat' with the inverse distances between all points in
    the group 'group'. Return the smallest distance between points. */
double Cartesian::calcInvDist(const CoordGroup &group, DistMatrix &mat) const
{
    double mindist(0);
    double tmpdist;

    int n = group.count();
    const Vector *array = group.constData();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n,n);

    for (int i=0; i<n; ++i)
    {
        const Vector &point = array[i];
        mat.setOuterIndex(i);

        //zero diagonal
        mat[i] = 0;

        for (int j=i+1; j<n; ++j)
        {
            //calculate the distance between the two points
            tmpdist = Vector::invDistance(point,array[j]);

            //store the minimum distance
            mindist = qMax(tmpdist,mindist);

            //place this distance into the matrix
            mat[j] = tmpdist;
            mat(j,i) = tmpdist;
        }
    }

    return 1.0 / mindist;
}

/** Populate the matrix 'mat' with the inverse distances^2 between all points in
    the group 'group'. Return the smallest distance between points. */
double Cartesian::calcInvDist2(const CoordGroup &group, DistMatrix &mat) const
{
    double mindist2(0);
    double tmpdist2;

    int n = group.count();
    const Vector *array = group.constData();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n,n);

    for (int i=0; i<n; ++i)
    {
        const Vector &point = array[i];
        mat.setOuterIndex(i);

        //zero diagonal
        mat[i] = 0;

        for (int j=i+1; j<n; ++j)
        {
            //calculate the distance between the two points
            tmpdist2 = Vector::invDistance2(point,array[j]);

            //store the minimum distance
            mindist2 = qMax(tmpdist2,mindist2);

            //place this distance into the matrix
            mat[j] = tmpdist2;
            mat(j,i) = tmpdist2;
        }
    }

    return sqrt( 1.0 / mindist2 );
}

/** Populate the matrix 'mat' with the distances between all of the
    points of the two CoordGroups. Return the shortest distance between the two
    CoordGroups. */
double Cartesian::calcDist(const CoordGroup &group0, const CoordGroup &group1,
                           DistMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());

    const int n0 = group0.count();
    const int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    #ifdef SIRE_TIME_ROUTINES
    int nflops = 0;
    #endif

    #ifdef SIRE_USE_SSE
    {
        //version of the code for processors with SSE2 or above
        const int remainder = n1 % 2;

        for (int i=0; i<n0; ++i)
        {
            const Vector& point0 = array0[i];
            mat.setOuterIndex(i);
            
            const double *p0 = point0.constData();
            
            __m128d sse_x0 = _mm_load1_pd( p0 );
            __m128d sse_y0 = _mm_load1_pd( p0 + 1 );
            __m128d sse_z0 = _mm_load1_pd( p0 + 2 );
            
            __m128d sse_mindist = _mm_load1_pd(&mindist);

            // Process points in pairs
            for (int j=0; j < n1-remainder; j+=2)
            {
                const Vector &point1 = array1[j];
                const Vector &point2 = array1[j+1];
                
                __m128d sse_x1 = _mm_setr_pd( point1.x(), point2.x() );
                __m128d sse_y1 = _mm_setr_pd( point1.y(), point2.y() );
                __m128d sse_z1 = _mm_setr_pd( point1.z(), point2.z() );
                
                __m128d delta = _mm_sub_pd(sse_x0, sse_x1);      // 2 flops
                __m128d tmpdist = _mm_mul_pd(delta, delta);      // 2 flops

                delta = _mm_sub_pd(sse_y0, sse_y1);              // 2 flops
                delta = _mm_mul_pd(delta, delta);                // 2 flops
                tmpdist = _mm_add_pd(tmpdist, delta);            // 2 flops
                
                delta = _mm_sub_pd(sse_z0, sse_z1);              // 2 flops
                delta = _mm_mul_pd(delta, delta);                // 2 flops
                tmpdist = _mm_add_pd(tmpdist, delta);            // 2 flops
                
                tmpdist = _mm_sqrt_pd(tmpdist);       // 2 flops

                #ifdef SIRE_TIME_ROUTINES
                nflops += 18;
                #endif

                sse_mindist = _mm_min_pd( sse_mindist, tmpdist );

                //place this distance into the matrix
                mat[j]   = *((const double*)&tmpdist);
                mat[j+1] = *( ((const double*)&tmpdist) + 1 );
            }
            
            mindist = qMin( *((const double*)&sse_mindist),
                            *( ((const double*)&sse_mindist) + 1 ) );
            
            if (remainder == 1)
            {
                const double tmpdist = Vector::distance(point0, array1[n1-1]);
                
                #ifdef SIRE_TIME_ROUTINES
                nflops += 9;
                #endif
                
                mindist = qMin(tmpdist,mindist);
                mat[n1-1] = tmpdist;
            }
        }
    }
    #else
    {
        //version suitable for all processors
        for (int i=0; i<n0; ++i)
        {
            const Vector &point0 = array0[i];
            mat.setOuterIndex(i);
            
            for (int j=0; j<n1; ++j)
            {
                const double tmpdist = Vector::distance(point0, array1[j]);
                
                #ifdef SIRE_TIME_ROUTINES
                nflops += 9;
                #endif
                
                mindist = qMin(mindist, tmpdist);
                
                mat[j] = tmpdist;
            }
        }
    }
    #endif

    #ifdef SIRE_TIME_ROUTINES
    ADD_FLOPS(nflops);
    #endif

    //return the minimum distance
    return mindist;
}

/** Populate the matrix 'mat' with the distances between all of the
    points of the passed CoordGroup and 'point'. Returns the shortest distance. */
double Cartesian::calcDist(const CoordGroup &group, const Vector &point,
                           DistMatrix &mat) const
{
    double mindist(std::numeric_limits<double>::max());

    const int n = group.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(1, n);

    //get raw pointers to the array - this provides more efficient access
    const Vector *array = group.constData();

    #ifdef SIRE_TIME_ROUTINES
    int nflops = 0;
    #endif

    #ifdef SIRE_USE_SSE
    {
        //version of the code for processors with SSE2 or above
        const int remainder = n % 2;

        mat.setOuterIndex(0);
        
        const double *p = point.constData();
        
        __m128d sse_x0 = _mm_load1_pd( p );
        __m128d sse_y0 = _mm_load1_pd( p + 1 );
        __m128d sse_z0 = _mm_load1_pd( p + 2 );
            
        __m128d sse_mindist = _mm_load1_pd(&mindist);

        // Process points in pairs
        for (int j=0; j < n-remainder; j+=2)
        {
            const Vector &point1 = array[j];
            const Vector &point2 = array[j+1];
                
            __m128d sse_x1 = _mm_setr_pd( point1.x(), point2.x() );
            __m128d sse_y1 = _mm_setr_pd( point1.y(), point2.y() );
            __m128d sse_z1 = _mm_setr_pd( point1.z(), point2.z() );
                
            __m128d delta = _mm_sub_pd(sse_x0, sse_x1);      // 2 flops
            __m128d tmpdist = _mm_mul_pd(delta, delta);      // 2 flops

            delta = _mm_sub_pd(sse_y0, sse_y1);              // 2 flops
            delta = _mm_mul_pd(delta, delta);                // 2 flops
            tmpdist = _mm_add_pd(tmpdist, delta);            // 2 flops
                
            delta = _mm_sub_pd(sse_z0, sse_z1);              // 2 flops
            delta = _mm_mul_pd(delta, delta);                // 2 flops
            tmpdist = _mm_add_pd(tmpdist, delta);            // 2 flops
                
            tmpdist = _mm_sqrt_pd(tmpdist);       // 2 flops

            #ifdef SIRE_TIME_ROUTINES
            nflops += 18;
            #endif

            sse_mindist = _mm_min_pd( sse_mindist, tmpdist );

            //place this distance into the matrix
            mat[j]   = *((const double*)&tmpdist);
            mat[j+1] = *( ((const double*)&tmpdist) + 1 );
        }
            
        mindist = qMin( *((const double*)&sse_mindist),
                        *( ((const double*)&sse_mindist) + 1 ) );
            
        if (remainder == 1)
        {
            const double tmpdist = Vector::distance(point, array[n-1]);
                
            #ifdef SIRE_TIME_ROUTINES
            nflops += 9;
            #endif
                
            mindist = qMin(tmpdist,mindist);
            mat[n-1] = tmpdist;
        }
    }
    #else
    {
        //version suitable for all processors
        mat.setOuterIndex(0);
            
        for (int j=0; j<n; ++j)
        {
            const double tmpdist = Vector::distance(point, array[j]);
                
            #ifdef SIRE_TIME_ROUTINES
            nflops += 9;
            #endif
                
            mindist = qMin(mindist, tmpdist);
                
            mat[j] = tmpdist;
        }
    }
    #endif

    #ifdef SIRE_TIME_ROUTINES
    ADD_FLOPS(nflops);
    #endif

    //return the minimum distance
    return mindist;
}

/** Populate the matrix 'mat' with the distances squared between all of the
    points of the passed CoordGroup and 'point'. Returns the shortest distance. */
double Cartesian::calcDist2(const CoordGroup &group, const Vector &point,
                            DistMatrix &mat) const
{
    double mindist2(std::numeric_limits<double>::max());

    const int n = group.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(1, n);

    //get raw pointers to the array - this provides more efficient access
    const Vector *array = group.constData();

    #ifdef SIRE_TIME_ROUTINES
    int nflops = 0;
    #endif

    #ifdef SIRE_USE_SSE
    {
        //version of the code for processors with SSE2 or above
        const int remainder = n % 2;

        mat.setOuterIndex(0);
        
        const double *p = point.constData();
        
        __m128d sse_x0 = _mm_load1_pd( p );
        __m128d sse_y0 = _mm_load1_pd( p + 1 );
        __m128d sse_z0 = _mm_load1_pd( p + 2 );
            
        __m128d sse_mindist2 = _mm_load1_pd(&mindist2);

        // Process points in pairs
        for (int j=0; j < n-remainder; j+=2)
        {
            const Vector &point1 = array[j];
            const Vector &point2 = array[j+1];
                
            __m128d sse_x1 = _mm_setr_pd( point1.x(), point2.x() );
            __m128d sse_y1 = _mm_setr_pd( point1.y(), point2.y() );
            __m128d sse_z1 = _mm_setr_pd( point1.z(), point2.z() );
                
            __m128d delta = _mm_sub_pd(sse_x0, sse_x1);      // 2 flops
            __m128d tmpdist2 = _mm_mul_pd(delta, delta);     // 2 flops

            delta = _mm_sub_pd(sse_y0, sse_y1);              // 2 flops
            delta = _mm_mul_pd(delta, delta);                // 2 flops
            tmpdist2 = _mm_add_pd(tmpdist2, delta);          // 2 flops
                
            delta = _mm_sub_pd(sse_z0, sse_z1);              // 2 flops
            delta = _mm_mul_pd(delta, delta);                // 2 flops
            tmpdist2 = _mm_add_pd(tmpdist2, delta);          // 2 flops

            #ifdef SIRE_TIME_ROUTINES
            nflops += 16;
            #endif

            sse_mindist2 = _mm_min_pd( sse_mindist2, tmpdist2 );

            //place this distance into the matrix
            mat[j]   = *((const double*)&tmpdist2);
            mat[j+1] = *( ((const double*)&tmpdist2) + 1 );
        }
            
        mindist2 = qMin( *((const double*)&sse_mindist2),
                         *( ((const double*)&sse_mindist2) + 1 ) );
            
        if (remainder == 1)
        {
            const double tmpdist2 = Vector::distance2(point, array[n-1]);
                
            #ifdef SIRE_TIME_ROUTINES
            nflops += 9;
            #endif
                
            mindist2 = qMin(tmpdist2,mindist2);
            mat[n-1] = tmpdist2;
        }
    }
    #else
    {
        //version suitable for all processors
        mat.setOuterIndex(0);
            
        for (int j=0; j<n; ++j)
        {
            const double tmpdist2 = Vector::distance2(point, array[j]);
                
            #ifdef SIRE_TIME_ROUTINES
            nflops += 8;
            #endif
                
            mindist2 = qMin(mindist2, tmpdist2);
                
            mat[j] = tmpdist2;
        }
    }
    #endif

    #ifdef SIRE_TIME_ROUTINES
    ADD_FLOPS(nflops);
    #endif

    //return the minimum distance
    return sqrt(mindist2);
}

/** Populate the matrix 'mat' with the distances^2 between all of the
    points of the two CoordGroups. Return the shortest distance between the
    two CoordGroups. */
double Cartesian::calcDist2(const CoordGroup &group0, const CoordGroup &group1,
                            DistMatrix &mat) const
{
    double mindist2(std::numeric_limits<double>::max());

    const int n0 = group0.count();
    const int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        const Vector& point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two points
            const double tmpdist = Vector::distance2(point0,array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist2 = qMin(tmpdist,mindist2);

            //place this distance into the matrix
            mat[j] = tmpdist;
        }
    }

    //return the minimum distance
    return sqrt(mindist2);
}

/** Populate the matrix 'mat' with the inverse distances between all of the
    points of the two CoordGroups. Return the shortest distance between
    the two CoordGroups. */
double Cartesian::calcInvDist(const CoordGroup &group0, const CoordGroup &group1,
                              DistMatrix &mat) const
{
    double maxinvdist(0);
    double tmpdist;

    int n0 = group0.count();
    int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        const Vector &point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two points
            tmpdist = Vector::invDistance(point0,array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            maxinvdist = qMax(tmpdist,maxinvdist);

            //place this distance into the matrix
            mat[j] = tmpdist;
        }
    }

    return 1.0 / maxinvdist;
}

/** Populate the matrix 'mat' with the inverse distances^2 between all of the
    points of the two CoordGroups. Return the shortest distance between
    the two CoordGroups. */
double Cartesian::calcInvDist2(const CoordGroup &group0, const CoordGroup &group1,
                               DistMatrix &mat) const
{
    double maxinvdist2(0);
    double tmpdist;

    int n0 = group0.count();
    int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        const Vector& point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two points
            tmpdist = Vector::invDistance2(point0, array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            maxinvdist2 = qMax(tmpdist,maxinvdist2);

            //place this distance into the matrix
            mat[j] = tmpdist;
        }
    }

    return sqrt( 1.0 / maxinvdist2 );
}

/** Calculate the distance vector between two points */
DistVector Cartesian::calcDistVector(const Vector &point0, 
                                     const Vector &point1) const
{
    return point1 - point0;
}
    
/** Populate the matrix 'distmat' with all of the interpoint distance vectors
    between all points within the CoordGroup. This is *not* a symmetrical matrix,
    as the direction from point A to point B is the negative of the 
    direction from point B to point A. This returns the shortest distance
    between two points in the group (that is not the self-self distance) */
double Cartesian::calcDistVectors(const CoordGroup &group, 
                                  DistVectorMatrix &mat) const
{
    double mindist = std::numeric_limits<double>::max();

    const int n = group.count();
    const Vector *array = group.constData();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n,n);

    for (int i=0; i<n; ++i)
    {
        const Vector &point = array[i];
        mat.setOuterIndex(i);

        //zero diagonal
        mat[i] = DistVector();

        for (int j=i+1; j<n; ++j)
        {
            mat[j] = (array[j] - point);

            //store the minimum distance
            mindist = qMin(mat[j].length(),mindist);

            //remember to store the negative of the distance
            mat(j,i) = -(mat[j]);
        }
    }

    return mindist;
}

/** Populate the matrix 'distmat' between all the points of the two CoordGroups
    'group1' and 'group2' - the returned matrix has the vectors pointing
    from each point in 'group1' to each point in 'group2'. This returns
    the shortest distance between two points in the group */
double Cartesian::calcDistVectors(const CoordGroup &group0, const CoordGroup &group1,
                                  DistVectorMatrix &mat) const
{
    double mindist = std::numeric_limits<double>::max();

    const int n0 = group0.count();
    const int n1 = group1.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(n0, n1);

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        const Vector& point0 = array0[i];
        mat.setOuterIndex(i);

        for (int j=0; j < n1; ++j)
        {
            mat[j] = (array1[j] - point0);
            mindist = qMin(mat[j].length(),mindist);
        }
    }

    //return the minimum distance
    return mindist;
}

/** Populate the matrix 'distmat' between all the points of the passed
    CoordGroup with 'point' - the returned matrix has the vectors pointing
    from the point, to each point in 'group'. This returns the shortest distance. */
double Cartesian::calcDistVectors(const CoordGroup &group, const Vector &point,
                                  DistVectorMatrix &mat) const
{
    double mindist = std::numeric_limits<double>::max();

    const int n = group.count();

    //redimension the matrix to hold all of the pairs
    mat.redimension(1, n);

    //get raw pointer to the array - this provides more efficient access
    const Vector *array = group.constData();

    mat.setOuterIndex(0);

    for (int j=0; j < n; ++j)
    {
        mat[j] = (array[j] - point);
        mindist = qMin(mat[j].length(),mindist);
    }

    //return the minimum distance
    return mindist;
}

/** Calculate the angle between the passed three points. This should return
    the acute angle between the points, which should lie between 0 and 180 degrees */
Angle Cartesian::calcAngle(const Vector &point0, const Vector &point1,
                           const Vector &point2) const
{
    return Vector::angle(point0, point1, point2);
}

/** Calculate the torsion angle between the passed four points. This should
    return the torsion angle measured clockwise when looking down the 
    torsion from point0-point1-point2-point3. This will lie between 0 and 360 
    degrees */
Angle Cartesian::calcDihedral(const Vector &point0, const Vector &point1,
                              const Vector &point2, const Vector &point3) const
{
    return Vector::dihedral(point0, point1, point2, point3);
}

/** Return whether or not two groups enclosed by the AABoxes 'aabox0'
    and 'aabox1' are definitely beyond the cutoff distance */
bool Cartesian::beyond(double dist, const AABox &aabox0, const AABox &aabox1) const
{
    #ifdef SIRE_TIME_ROUTINES
    ADD_FLOPS( 8 );
    #endif

    return Vector::distance2(aabox0.center(), aabox1.center()) >
                SireMaths::pow_2(dist + aabox0.radius() + aabox1.radius());
}

/** Return whether or not these two groups are definitely beyond the cutoff distance. */
bool Cartesian::beyond(double dist, const CoordGroup &group0,
                       const CoordGroup &group1) const
{
    return Cartesian::beyond(dist, group0.aaBox(), group1.aaBox());
}

/** Return the minimum distance between the points in 'group0' and 'group1'. */
double Cartesian::minimumDistance(const CoordGroup &group0,
                                  const CoordGroup &group1) const
{
    double mindist2(std::numeric_limits<double>::max());

    int n0 = group0.count();
    int n1 = group1.count();

    //get raw pointers to the arrays - this provides more efficient access
    const Vector *array0 = group0.constData();
    const Vector *array1 = group1.constData();

    for (int i=0; i<n0; ++i)
    {
        const Vector& point0 = array0[i];

        for (int j=0; j<n1; ++j)
        {
            //calculate the distance between the two points
            double tmpdist = Vector::distance2(point0,array1[j]);

            //store the minimum distance, the value expected to be the minimum
            //value is most efficiently placed as the second argument
            mindist2 = qMin(tmpdist,mindist2);
        }
    }

    //return the minimum distance
    return sqrt(mindist2);
}

/** Return the minimum distance between the two passed boxes */
double Cartesian::minimumDistance(const AABox &box0, const AABox &box1) const
{
    Vector delta = (box0.center() - box1.center());
    delta = Vector( std::abs(delta.x()), std::abs(delta.y()), std::abs(delta.z()) );
    
    delta -= box0.halfExtents();
    delta -= box1.halfExtents();
    
    delta = delta.max( Vector(0) );
    
    return delta.length();
}

/** Return the minimum distance between a point and a box */
double Cartesian::minimumDistance(const Vector &point, const AABox &box) const
{
    Vector delta = point - box.center();
    delta = Vector( std::abs(delta.x()), std::abs(delta.y()), std::abs(delta.z()) );
    
    delta -= box.halfExtents();
    
    delta = delta.max( Vector(0) );
    
    return delta.length();
}

/** Return the minimum distance between points within the group 'group'. */
double Cartesian::minimumDistance(const CoordGroup &group) const
{
    double mindist2(std::numeric_limits<double>::max());

    int n = group.count();
    const Vector *array = group.constData();

    for (int i=0; i<n; ++i)
    {
        const Vector &point = array[i];

        for (int j=i+1; j<n; ++j)
        {
            //calculate the distance between the two points
            double tmpdist2 = Vector::distance2(point,array[j]);

            //store the minimum distance
            mindist2 = qMin(tmpdist2,mindist2);
        }
    }

    return sqrt(mindist2);
}

/** Return the minimum image copy of 'group' with respect to 'center'.
    In this case, as this is not a periodic space, this just returns
    'group' */
CoordGroup Cartesian::getMinimumImage(const CoordGroup &group, const Vector&) const
{
    return group;
}

/** Return the minimum image copy of 'groups' with respect to 'center'.
    In this case, as this is not a periodic space, this just returns
    'groups' */
CoordGroupArray Cartesian::getMinimumImage(const CoordGroupArray &groups,
                                           const Vector&, bool) const
{
    return groups;
}

/** A cartesian space is not periodic, so this just returns the input 'aabox' */
AABox Cartesian::getMinimumImage(const AABox &aabox, const Vector&) const
{
    return aabox;
}

/** A cartesian space is not periodic, so this just returns the input 'point' */
Vector Cartesian::getMinimumImage(const Vector &point, const Vector&) const
{
    return point;
}

/** Return all periodic images of 'point' with respect to 'center' within
    'dist' distance of 'center' */
QVector<Vector> Cartesian::getImagesWithin(const Vector &point, const Vector &center,
                                           double dist) const
{
    QVector<Vector> points;

    if ( Vector::distance(point,center) < dist )
    {
        points.append(point);
    }
    
    return points;
}

/** Return a list of copies of CoordGroup 'group' that are within
    'distance' of the CoordGroup 'center', translating 'group' so that
    it has the right coordinates to be around 'center'. As this is not
    a periodic space, this will merely return a copy of 'group' if
    it is within the specified distance. */
QList< tuple<double,CoordGroup> >
Cartesian::getCopiesWithin(const CoordGroup &group, const CoordGroup &center,
                           double dist) const
{
    QList< tuple<double,CoordGroup> > closegroups;

    if (not this->beyond(dist, group, center))
    {
        //calculate the minimum distance - do this via the
        //minimumDistance function
        double mindist = this->minimumDistance(group, center);

        if (mindist <= dist)
        {
            closegroups.append( tuple<double,CoordGroup>(mindist,group) );
        }
    }

    return closegroups;
}

/** Return a random point in this space - this can be truly anywhere!!!

    (well, it is limited to within -10^20 and 10^20 angstroms)
*/
Vector Cartesian::getRandomPoint(const Vector&, const RanGenerator &generator) const
{
    static const double max_rand = 1.0e20;
    static const double min_rand = -max_rand;

    return Vector( generator.rand(min_rand, max_rand),
                   generator.rand(min_rand, max_rand),
                   generator.rand(min_rand, max_rand) );
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at the origin */
Vector Cartesian::getBoxCenter(const Vector&) const
{
	return Vector(0,0,0);
}

/** Return the center of the box that contains the point 'p' assuming
    that the center for the central box is located at 'center' */
Vector Cartesian::getBoxCenter(const Vector&, const Vector &center) const
{
	return center;
}

const char* Cartesian::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Cartesian>() );
}

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

#include "sphere.h"
#include "rangenerator.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include <QDebug>
#include <QElapsedTimer>

using namespace SireStream;
using namespace SireMaths;

static const RegisterMetaType<Sphere> r_sphere(NO_ROOT);

/** Serialise to a binary data stream */
QDataStream &operator<<(QDataStream &ds, const Sphere &sphere)
{
    writeHeader(ds, r_sphere, 1) << sphere._center << sphere._radius;

    return ds;
}

/** Deserialise from a binary data stream */
QDataStream &operator>>(QDataStream &ds, Sphere &sphere)
{
    VersionID v = readHeader(ds, r_sphere);

    if (v == 1)
    {
        ds >> sphere._center >> sphere._radius;
    }
    else
        throw version_error(v, "1", r_sphere, CODELOC);

    return ds;
}

/** Create a default sphere, centered on the origin with zero radius */
Sphere::Sphere() : _center(0.0), _radius(0.0)
{}

/** Create a sphere with radius 'radius', centered on the origin */
Sphere::Sphere(const double &radius) : _center(0.0), _radius(radius)
{
    if (radius < 0.0)
        _radius = -radius;
}

/** Create a sphere centered at 'position' and with radius 'radius' */
Sphere::Sphere(const Vector &position, const double &radius)
       : _center(position), _radius(radius)
{
    if (radius < 0.0)
        _radius = -radius;
}

/** Copy constructor */
Sphere::Sphere(const Sphere &other)
       : _center(other._center), _radius(other._radius)
{}

/** Destructor */
Sphere::~Sphere()
{}

/** Comparison operator */
bool Sphere::operator==(const Sphere &other) const
{
    return _radius == other._radius and _center == other._center;
}

/** Comparison operator */
bool Sphere::operator!=(const Sphere &other) const
{
    return _radius != other._radius or _center != other._center;
}

QString Sphere::toString() const
{
    return QObject::tr("Sphere( center() == %1, radius == %2 )")
                .arg(center().toString()).arg(radius());
}

/** Return a copy of this sphere translated by 'delta' */
Sphere Sphere::translate(const Vector &delta) const
{
    return Sphere( _center + delta, _radius );
}

/** Set the position of the center of this sphere */
void Sphere::setPosition(const Vector &position)
{
    _center = position;
}

/** Set the position of the center of this sphere */
void Sphere::setCenter(const Vector &center)
{
    _center = center;
}

/** Set the radius of this sphere */
void Sphere::setRadius(double radius)
{
    _radius = radius;

    if (radius < 0.0)
        _radius = 0.0;
}

/** Return the volume of this sphere */
double Sphere::volume() const
{
    return (4.0/3.0) * SireMaths::pi * SireMaths::pow_3(_radius);
}

/** Return the surface area of this sphere */
double Sphere::surfaceArea() const
{
    return 4.0 * SireMaths::pi * SireMaths::pow_2(_radius);
}

/** Return whether or not this sphere intersects with 'other' */
bool Sphere::intersects(const Sphere &other) const
{
    return Vector::distance(_center, other._center) < (_radius + other._radius);
}

/** Return whether or not this sphere contains the point 'point'
    (returns true even if the point is just on the surface of the sphere) */
bool Sphere::contains(const Vector &point) const
{
    return Vector::distance(_center, point) <= _radius;
}

/** Return whether or not this sphere contains 'other' */
bool Sphere::contains(const Sphere &other) const
{
    return Vector::distance(_center, other._center) + other._radius <= _radius;
}

/** Return the volume of space formed at the intersection of this sphere with 'other' */
double Sphere::intersectionVolume(const Sphere &other) const
{
    if (this->radius() < other.radius())
        return other.intersectionVolume(*this);

    const double dist = Vector::distance(_center, other._center);
    const double sum_radii = _radius + other._radius;
    const double big_radius = this->radius();
    const double small_radius = other.radius();
    
    if (dist >= sum_radii)
    {
        //no intersection
        return 0;
    }
    else if (dist + small_radius <= big_radius)
    {
        //the small sphere is entirely enclosed in the big sphere
        return (4.0/3.0) * SireMaths::pi * SireMaths::pow_3(small_radius);
    }
    else
    {
        //the intersection volume is the two spherical caps cut off at
        //the top of two spheres
        
        // we can work out the height of the two caps by imagining
        // two triangles formed by the intersection of the spheres
        // (x axis is along the sphere centers, r0 is radius 0, r1 is radius 1,
        //  h0 is the cap distance for sphere 0, h1 is the cap distance for sphere 1,
        //  and d is the distance between the center of the spheres
        
        //         O         O
        //   r0  / |  x    x | \  r1
        //     /   |         |   \
        //    O----O         O----O
        //     r0-h0         r1-h1
        
        //  this gives us three equations with three unknowns (x, h0 and h1)
        //
        //  (1) d = r0 - h0 + r1 - h1
        //  (2) (r0-h0)**2 + x**2 = r0**2
        //  (3) (r1-h1)**2 + x**2 = r1**2
        //
        //  Rearranging these to solve for h0 and h1 gives us;
        //
        //  (i)  h0 = [ r1**2 - (d - r0)**2 ] / 2d
        //  (ii) h1 = r0 - h0 + r1 - d
        
        const double d = dist;
        const double r0 = _radius;
        const double r1 = other._radius;
        
        const double h0 = (r1*r1 - SireMaths::pow_2(d-r0)) / (2*d);
        const double h1 = r0 - h0 + r1 - d;
        
        //we can check this is true by calculating x...
        //const double x0 = r0*r0 - SireMaths::pow_2(r0-h0);
        //const double x1 = r1*r1 - SireMaths::pow_2(r1-h1);
        //qDebug() << h0 << h1 << x0 << x1;
        
        // now the volume will be
        // volume(sphere 0) + volume(sphere 1) - cap(sphere 0) - cap(sphere 1)
        //
        // (where the cap has volume (1/3) pi h**2 (3r - h)
        const double cap0 = (1.0/3.0) * SireMaths::pi * h0 * h0 * (3*r0 - h0);
        const double cap1 = (1.0/3.0) * SireMaths::pi * h1 * h1 * (3*r1 - h1);
        
        return cap0 + cap1;
    }
}

SIRE_ALWAYS_INLINE double my_arctan(double val)
{
    const double result = std::atan(val);
    
    //this outputs the angle in the range -pi/2 to pi/2
    //We need this to be in the range 0 < pi
    
    if (result < 0)
        return SireMaths::pi + result;
    else
        return result;
}

/** Return the volume of intersection of this sphere with the two other spheres.
    This returns the volume of space covered by all three spheres. */
double Sphere::intersectionVolume(const Sphere &sphere1, const Sphere &sphere2) const
{
    //process in descending order of radii
    if (this->radius() < sphere1.radius())
    {
        if (this->radius() < sphere2.radius())
        {
            if (sphere1.radius() < sphere2.radius())
            {
                return sphere2.intersectionVolume(sphere1, *this);
            }
            else
            {
                return sphere1.intersectionVolume(sphere2, *this);
            }
        }
        else
        {
            return sphere1.intersectionVolume(*this, sphere2);
        }
    }
    else if (this->radius() < sphere2.radius())
    {
        return sphere2.intersectionVolume(*this, sphere1);
    }
    else if (sphere1.radius() < sphere2.radius())
    {
        return this->intersectionVolume(sphere2, sphere1);
    }
    
    const Sphere &sphere0 = *this;
    
    if (sphere0.radius() < sphere1.radius() or sphere0.radius() < sphere2.radius()
           or sphere1.radius() < sphere2.radius())
    {
        throw SireError::program_bug( QObject::tr(
                "Spheres should be arranged in size order? %1 %2 %3")
                    .arg(sphere0.toString()).arg(sphere1.toString())
                    .arg(sphere2.toString()), CODELOC );
    }
    
    //all three spheres must be intersecting...
    if ( not (sphere0.intersects(sphere1) and sphere0.intersects(sphere2) and
              sphere1.intersects(sphere2)) )
    {
        //they are not intersecting
        return 0;
    }
    
    //if one sphere contains another, then the intersection must be the volume
    //of intersection of the smaller spheres
    if (sphere0.contains(sphere1) or sphere0.contains(sphere2))
    {
        return sphere1.intersectionVolume(sphere2);
    }
    else if (sphere1.contains(sphere2))
    {
        return this->intersectionVolume(sphere1);
    }
    
    //use the method derived by Gibson and Scheraga in J. Phys. Chem. 1987, 91, 4121-4122
    const Sphere &A = sphere0;
    const Sphere &B = sphere1;
    const Sphere &C = sphere2;
    
    //quote from paper; From figure 1, a, b, and c are the distances between the centers
    //of the spheres B and C, C and A, and A and B respectively
    const double a = Vector::distance(B.center(), C.center());
    const double b = Vector::distance(C.center(), A.center());
    const double c = Vector::distance(A.center(), B.center());
    
    //alpha, beta and gamma are the radii of the spheres centered at A, B and C
    const double alpha = A.radius();
    const double beta = B.radius();
    const double gamma = C.radius();
    
    const double a2 = a*a;
    const double b2 = b*b;
    const double c2 = c*c;
    
    const double alpha2 = alpha*alpha;
    const double beta2 = beta*beta;
    const double gamma2 = gamma*gamma;

    //to simplify the expression of the volume of the triple intersection we define
    
    //epsilon1 = (beta^2 - gamma^2) / a^2
    const double epsilon1 = (beta2 - gamma2) / a2;
    
    //epsilon2 = (gamma^2 - alpha^2) / b^2
    const double epsilon2 = (gamma2 - alpha2) / b2;
    
    //epsilon3 = (alpha^2 - beta^2) / c^2
    const double epsilon3 = (alpha2 - beta2) / c2;
    
    //w^2 = (alpha^2 a^2 + beta^2 b^2 + gamma^2 c^2)(a^2 + b^2 + c^2) -
    //         2(alpha^2 a^4 + beta^2 b^4 + gamma^2 c^4) +
    //             a^2 b^2 c^2 (epsilon1 epsilon2 + epsilon2 epsilon3 + epsilon3 epsilon1 - 1)
    
    const double w2 = ((alpha2*a2 + beta2*b2 + gamma2*c2)*(a2 + b2 + c2)) -
                      (2.0*(alpha2*a2*a2 + beta2*b2*b2 + gamma2*c2*c2)) +
                      (a2*b2*c2*(epsilon1*epsilon2 + epsilon2*epsilon3 + epsilon3*epsilon1 - 1.0));
    
    //we can use w2 to work out how the spheres intersect
    qDebug() << "w2" << w2;
    
    // if w2 == 0 then the spheres intersect at a single point, and the volume is zero
    if (w2 == 0)
    {
        return 0;
    }
    // if w2 > 0 then the spheres intersect in two points, and we can use
    // the method in the paper to find the volume
    else if (w2 > 0)
    {
        const double w = std::sqrt(w2);
        const double two_w = 2.0 * w;
        
        //q1 = a [ b^2 + c^2 - a^2 + beta^2 + gamma^2 - 2 alpha^2 + epsilon1 (b^2 - c^2) ]
        const double q1 = a * (b2 + c2 - a2 + beta2 + gamma2 - 2.0*alpha2 +
                                 epsilon1*(b2-c2));

        //q2 = b [ c^2 + a^2 - b^2 + gamma^2 + alpha^2 - 2 beta^2 + epsilon2 (c^2 - a^2) ]
        const double q2 = b * (c2 + a2 - b2 + gamma2 + alpha2 - 2.0*beta2 +
                                 epsilon2*(c2-a2));
        
        //q3 = c [ a^2 + b^2 - c^2 + alpha^2 + beta^2 - 2 gamma^2 + epsilon3 (a^2 - b^2) ]
        const double q3 = c * (a2 + b2 - c2 + alpha2 + beta2 - 2.0*gamma2 +
                                 epsilon3*(a2-b2));

        //we can now calculate the volume of the intersection using the formula in the paper
        // V_ABC =
        //    (w/6) -
        //    (a/2)[ beta^2 + gamma^2 - a^2( 1/6 - epsilon1^2/2) ] tan-1 (2w/q1) -
        //    (b/2)[ gamma^2 + alpha^2 - b^2( 1/6 - epsilon2^2/2) ] tan-1 (2w/q2) -
        //    (c/2)[ alpha^2 + beta^2 - c^2( 1/6 - epsilon3^2/2) ] tan-1 (2w/q3) +
        //    (2/3 alpha^3)[ tan-1{(bw/alpha q2)(1 - epsilon2)} +
        //                   tan-1{(cw/alpha q3)(1 + epsilon3)} ] +
        //    (2/3 beta^3)[ tan-1{(cw/beta q3)(1 - epsilon3)} +
        //                  tan-1{(aw/beta q1)(1 + epsilon1)} ] +
        //    (2/3 gamma^3)[ tan-1{(aw/gamma q1)(1 - epsilon1)} +
        //                   tan-1{(bw/gamma q2)(1 + epsilon2)} ]
        //
        //    where tan-1 returns values in the range  0 <= tan-1 <= pi
        
        const double alpha3 = alpha2 * alpha;
        const double beta3 = beta2 * beta;
        const double gamma3 = gamma2 * gamma;
        const double two_thirds = 2.0 / 3.0;
        const double one_sixth = 1.0 / 6.0;

        const double v_abc = (w / 6.0) -
                    ((a / 2.0) * (beta2 + gamma2 - a2*(one_sixth - 0.5*epsilon1*epsilon1))
                               * my_arctan( two_w / q1 )) -
                    ((b / 2.0) * (gamma2 + alpha2 - b2*(one_sixth - 0.5*epsilon2*epsilon2))
                               * my_arctan( two_w / q2 )) -
                    ((c / 2.0) * (alpha2 + beta2 - c2*(one_sixth - 0.5*epsilon3*epsilon3))
                               * my_arctan( two_w / q3 )) +
                 ((two_thirds * alpha3) * (my_arctan( ((b*w)/(alpha*q2)) * (1.0-epsilon2) ) +
                                           my_arctan( ((c*w)/(alpha*q3)) * (1.0+epsilon3) )) ) +
                 ((two_thirds * beta3) * (my_arctan( ((c*w)/(beta*q3)) * (1.0-epsilon3) ) +
                                          my_arctan( ((a*w)/(beta*q1)) * (1.0+epsilon1) )) ) +
                 ((two_thirds * gamma3) * (my_arctan( ((a*w)/(gamma*q1)) * (1.0-epsilon1) ) +
                                           my_arctan( ((b*w)/(gamma*q2)) * (1.0+epsilon2) )) );

        //this should have worked...
        return v_abc;
    }
    else
    {
        qDebug() << "SPECIAL CASE";

        //special cases... there are three possibilities
        //define;
        //t^2 = (a+b+c)(-a+b+c)(a-b+c)(a+b-c)
        const double t2 = (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c);
        
        //t_abg^2 = (a+beta+gamma)(-a+beta+gamma)(a-beta+gamma)(a+beta-gamma)
        const double t_abg2 = (a+beta+gamma)*(-a+beta+gamma)*(a-beta+gamma)*(a+beta-gamma);
        
        //cyclically permute a,b,c alpha,beta,gamma in above to get t_bga^2 and t_cab
        const double t_bga2 = (b+gamma+alpha)*(-b+gamma+alpha)*(b-gamma+alpha)*(b+gamma-alpha);
        const double t_cab2 = (c+alpha+beta)*(-c+alpha+beta)*(c-alpha+beta)*(c+alpha-beta);
        
        //now set p1 and p2 via
        //p1 = [(b^2 - c^2 + beta^2 - gamma^2)^2 + (t - t_abg)^2] / 4 a^2 - alpha^2
        //p2 = [(b^2 - c^2 + beta^2 - gamma^2)^2 + (t + t_abg)^2] / 4 a^2 - alpha^2
        
        const double p1 = (SireMaths::pow_2(b2 - c2 + beta2 - gamma2) +
                            SireMaths::pow_2(std::sqrt(t2) - std::sqrt(t_abg2))) / (4.0*a2)
                                - alpha2;
        const double p2 = (SireMaths::pow_2(b2 - c2 + beta2 - gamma2) +
                            SireMaths::pow_2(std::sqrt(t2) + std::sqrt(t_abg2))) / (4.0*a2)
                                - alpha2;
    
        //now calculate p3 and p4 by cyclically permuting a,b,c and alpha,beta,gamma
        const double p3 = (SireMaths::pow_2(c2 - a2 + gamma2 - alpha2) +
                            SireMaths::pow_2(std::sqrt(t2) - std::sqrt(t_bga2))) / (4.0*b2)
                                - beta2;
        const double p4 = (SireMaths::pow_2(c2 - a2 + gamma2 - alpha2) +
                            SireMaths::pow_2(std::sqrt(t2) + std::sqrt(t_bga2))) / (4.0*b2)
                                - beta2;

        //repeat again to get p5 and p6
        const double p5 = (SireMaths::pow_2(a2 - b2 + alpha2 - beta2) +
                            SireMaths::pow_2(std::sqrt(t2) - std::sqrt(t_cab2))) / (4.0*c2)
                                - gamma2;
        const double p6 = (SireMaths::pow_2(a2 - b2 + alpha2 - beta2) +
                            SireMaths::pow_2(std::sqrt(t2) + std::sqrt(t_cab2))) / (4.0*c2)
                                - gamma2;
    
        qDebug() << p1 << p2 << p3 << p4 << p5 << p6;
    
        // if all pi are > 0 then we have case 1 - each pairwise intersection lies outside
        // the third sphere, and so there is no triple intersection
        if (p1 >= 0 and p2 >= 0 and p3 >= 0 and p4 >= 0 and p5 >= 0 and p6 >= 0)
        {
            return 0;
        }
        //if p1 and p2 < 0, and all others and > 0 we have case 2 - the sphere centered
        // on C contains the pairwise intersection of A and B
        else if (p1 < 0 and p2 < 0 and p3 >= 0 and p4 >= 0 and p5 >= 0 and p6 >= 0)
        {
            return A.intersectionVolume(B);
        }
        //similarly, if p3 and p4 are < 0 but the others are greater, then have case 2
        // where sphere A contains the intersection of B and C
        else if (p3 < 0 and p4 < 0 and p1 >= 0 and p2 >= 0 and p5 >= 0 and p6 >= 0)
        {
            return B.intersectionVolume(C);
        }
        //similarly, if p5 and p6 are < 0 but the others are greater, then have case 2
        // where sphere B contains the intersection of C and A
        else if (p5 < 0 and p6 < 0 and p1 >= 0 and p2 >= 0 and p3 >= 0 and p4 >= 0)
        {
            return C.intersectionVolume(A);
        }
        //else if p1 and p2 > 0 but the others are negative, then we have a triple
        //intersection even though the surfaces of the spheres have no common point
        else if (p1 > 0 and p2 > 0 and p3 <= 0 and p4 <= 0 and p5 <= 0 and p6 <= 0)
        {
            return A.intersectionVolume(B) + A.intersectionVolume(C) - A.volume();
        }
        //else if p3 and p4 > 0 then we have the same case, but cyclically permuted
        else if (p3 > 0 and p4 > 0 and p1 <= 0 and p2 <= 0 and p5 <= 0 and p6 <= 0)
        {
            return B.intersectionVolume(C) + B.intersectionVolume(A) - B.volume();
        }
        //else if p5 and p6 > 0 then we have the same case, but cyclically permuted
        else if (p5 > 0 and p6 > 0 and p1 <= 0 and p2 <= 0 and p3 <= 0 and p4 <= 0)
        {
            return C.intersectionVolume(A) + C.intersectionVolume(B) - C.volume();
        }
        else
        {
            //don't know how we got here
            qDebug() << "HOW DID WE GET HERE?"
                     << p1 << p2 << p3 << p4 << p5 << p6;
            
            return 0;
        }
    }
}

/** Return the combined volume of the passed array of spheres. This calculates the volume
    analytically using the inclusion/exclusion principle only up to third order
    (intersections of three spheres). The volume is calculated as;
    
     sum{ volume of all spheres - 
          volume of all pair intersections + 
          volume of all triple intersections } 
          
    This performs quite well in most cases, but will have increasing error as
    the number of quadruple or higher intersections increases. If this would be a problem,
    (e.g. exact vdw volume of benzene rings) then use combinedVolumeMC that 
    uses Monte Carlo integration to find the approximate volume of the spheres.
*/
double Sphere::combinedVolume(const QVector<SireMaths::Sphere> &spheres)
{
    if (spheres.isEmpty())
    {
        return 0;
    }
    else if (spheres.count() == 1)
    {
        return spheres.at(0).volume();
    }
    else if (spheres.count() == 2)
    {
        //volume is sum of volumes minus the volume of overlap
        return spheres.at(0).volume() + spheres.at(1).volume()
                 - spheres.at(0).intersectionVolume(spheres.at(1));
    }
    else
    {
        QElapsedTimer timer;
        timer.start();
    
        //calculate the volume of all spheres
        double total_volume = 0;
        
        for (int i=0; i<spheres.count(); ++i)
        {
            total_volume += spheres.constData()[i].volume();
        }
        
        //now find all intersecting pairs of spheres and subtract their volume
        //from the total
        QVector< QSet<int> > intersecting_pairs( spheres.count() );
        
        for (int i=0; i<spheres.count()-1; ++i)
        {
            const Sphere &sphere0 = spheres.constData()[i];
        
            for (int j=i+1; j<spheres.count(); ++j)
            {
                const Sphere &sphere1 = spheres.constData()[j];
            
                if (sphere0.intersects(sphere1))
                {
                    total_volume -= sphere0.intersectionVolume(sphere1);
                    intersecting_pairs[i].insert(j);
                    intersecting_pairs[j].insert(i);
                }
            }
        }
        
        //now find all of the sets of triple intersections
        for (int i=0; i<spheres.count(); ++i)
        {
            if (intersecting_pairs[i].count() > 1)
            {
                //this sphere is intersecting with more than one other sphere - are
                //these spheres intersecting with each other too? If they are, then
                //there is a good chance that there is a combined volume between these
                //spheres
                QList<int> intersects = intersecting_pairs[i].toList();
                qSort(intersects);
                
                const int sphere0 = i;
                
                for (int j=0; j<intersects.count()-1; ++j)
                {
                    const int sphere1 = intersects[j];
                    
                    if (sphere1 > sphere0)
                    {
                        for (int k=j+1; k<intersects.count(); ++k)
                        {
                            //sphere2 must be larger than sphere1 as we sorted the
                            //intersects list
                            const int sphere2 = intersects[k];
                        
                            if (intersecting_pairs[sphere1].contains(sphere2))
                            {
                                //sphere0 intersects with sphere1 and sphere2
                                //sphere1 intersects with sphere2
                                //so all three are intersecting together
                                total_volume += spheres[sphere0]
                                                    .intersectionVolume(spheres[sphere1],
                                                                        spheres[sphere2]);
                            }
                        }
                    }
                }
            }
        }
        
        qint64 nsecs = timer.nsecsElapsed();
        
        qDebug() << "Analytic volume calculation took" << (0.000001*nsecs) << "ms";
        
        return total_volume;
    }
}

/** Return an approximation of the combined volume of the set of passed spheres using
    Monte Carlo sampling. If 'nsamples' is greater than zero, then the set number of
    samples will be used. Otherwise, enough samples will be used to converge the volume
    to within a good approximation */
double Sphere::combinedVolumeMC(const QVector<Sphere> &spheres, qint64 nsamples)
{
    if (spheres.isEmpty())
        return 0;

    //create a box that encompasses all of the spheres
    Vector minbox = spheres.at(0).center();
    Vector maxbox = spheres.at(1).center();
    
    for (int i=0; i<spheres.count(); ++i)
    {
        minbox = minbox.min( spheres.at(i).center() - Vector(spheres.at(i).radius()) );
        maxbox = maxbox.max( spheres.at(i).center() + Vector(spheres.at(i).radius()) );
    }

    //now generate points randomly in this box...
    RanGenerator rand;
    
    quint64 n_in_sphere = 0;
    quint64 n_tests = 0;
    
    if (nsamples > 0)
    {
        for (qint64 i=0; i<nsamples; ++i)
        {
            Vector test( rand.rand( minbox.x(), maxbox.x() ),
                         rand.rand( minbox.y(), maxbox.y() ),
                         rand.rand( minbox.z(), maxbox.z() ) );
            
            //is this point inside any of the spheres?
            for (int j=0; j<spheres.count(); ++j)
            {
                if (spheres.constData()[j].contains(test))
                {
                    n_in_sphere += 1;
                    break;
                }
            }
        }
        
        n_tests = nsamples;
    }
    else
    {
        double ratio = 0;
        double old_ratio = 0;
        
        QElapsedTimer t;
        t.start();
        qint64 last_timeout = 2000;
        
        while (true)
        {
            n_tests += 1;
        
            Vector test( rand.rand( minbox.x(), maxbox.x() ),
                         rand.rand( minbox.y(), maxbox.y() ),
                         rand.rand( minbox.z(), maxbox.z() ) );
            
            //is this point inside any of the spheres?
            for (int j=0; j<spheres.count(); ++j)
            {
                if (spheres.constData()[j].contains(test))
                {
                    n_in_sphere += 1;
                    break;
                }
            }
            
            if (n_tests > 1000 and (n_tests % 1000 == 0))
            {
                ratio = double(n_in_sphere) / double(n_tests);
                
                const double conv_level = 1e-8;
                
                if (std::abs(ratio - old_ratio) < conv_level)
                {
                    //converged
                    break;
                }
                
                if (t.elapsed() > last_timeout)
                {
                    last_timeout += 2000;
                    qDebug() << "Solving... SAMPLES:" << n_tests
                             << "CONVERGENCE:" << std::abs(ratio-old_ratio)
                             << "DESIRED LEVEL:" << conv_level;
                }
                
                old_ratio = ratio;
            }
        }
    }

    const double total_volume = double(n_in_sphere) * (maxbox.x()-minbox.x()) *
                                                      (maxbox.y()-minbox.y()) *
                                                      (maxbox.z()-minbox.z()) / double(n_tests);

    return total_volume;
}

const char* Sphere::typeName()
{
    return QMetaType::typeName( qMetaTypeId<Sphere>() );
}

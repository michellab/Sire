/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2013  Christopher Woods
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

#ifndef SIREMATHS_SPHEREPROPERTY_H
#define SIREMATHS_SPHEREPROPERTY_H

#include "SireBase/podproperty.hpp"

#include "SireMaths/sphere.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

typedef SireBase::PODProperty<Sphere> SphereProperty;
typedef SireBase::PODArrayProperty<Sphere> SphereArrayProperty;

SIREMATHS_EXPORT SireBase::PropertyPtr wrap(const Sphere &sphere);
SIREMATHS_EXPORT SireBase::PropertyPtr wrap(const QVector<Sphere> &spheres);
SIREMATHS_EXPORT SireBase::PropertyPtr wrap(const QList<Sphere> &spheres);

}

Q_DECLARE_METATYPE(SireMaths::SphereProperty);
Q_DECLARE_METATYPE(SireMaths::SphereArrayProperty);

SIRE_EXPOSE_FUNCTION(SireMaths::wrap)

SIRE_EXPOSE_ALIAS(SireBase::PODProperty<SireMaths::Sphere>,
                  SireMaths::SphereProperty);
SIRE_EXPOSE_ALIAS(SireBase::PODArrayProperty<SireMaths::Sphere>,
                  SireMaths::SphereArrayProperty);

#ifdef SIRE_INSTANTIATE_TEMPLATES
template class SireBase::PODProperty<SireMaths::Sphere>;
template class SireBase::PODArrayProperty<SireMaths::Sphere>;
#endif

SIRE_END_HEADER

#endif

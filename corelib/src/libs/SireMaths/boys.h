/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
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

#ifndef SIREMATHS_BOYS_H
#define SIREMATHS_BOYS_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMaths
{

double boys_f0(double x);
double boys_f1(double x);
double boys_f2(double x);
double boys_f3(double x);

double boys(double m, double x);
double boys(int m, double x);

void multi_boys_2(double x, double boys[2]);
void multi_boys_3(double x, double boys[3]);
void multi_boys_4(double x, double boys[4]);
void multi_boys_n(double x, double boys[], int n);

void multi_boys_2(double x, double boys[2], int start);
void multi_boys_3(double x, double boys[3], int start);
void multi_boys_4(double x, double boys[4], int start);
void multi_boys_n(double x, double boys[], int n, int start);

QVector<double> multi_boys(double x, int n);
QVector<double> multi_boys(double x, int n, int start);

}

SIRE_EXPOSE_FUNCTION( SireMaths::boys )
SIRE_EXPOSE_FUNCTION( SireMaths::multi_boys )

SIRE_END_HEADER

#endif

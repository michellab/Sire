/*
 * regress.h
 *
 * Polynomial regression for free energy estimates
 * Copyright (C) 2008   Conrad Shyu (shyu4751@yahoo.com)
 * Department of Physics, University of Idaho
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author's comments
 * -----------------
 * This program fits a polynomial regression model to powers of a single
 * predictor by the method of linear least squares. Interpolation and
 * calculation of areas under the curve are also given.
 *
 * written by Conrad Shyu (shyu4751@yahoo.com)
 * Department of Physics
 * University of Idaho, Moscow, ID 83844
 *
 * first created on December 27, 2007
 * revised on December 30, 2007
*/

#ifndef _REGRESS_H
#define _REGRESS_H

#include "sireglobal.h"
#include "SireStream/datastream.h"
#include <list>
#include <cmath>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>

using namespace std;

const double TOLERANCE_LEVEL = 1e-30;
const unsigned int POLYNOMIAL_DEGREE = 10;

typedef struct
{
    double x;   // positions on the x-axis
    double y;   // values on the y-axis, y=f(x)
} stREGRESS;

class Regress
{
public:
    Regress(
        const list<stREGRESS>&,                 // structure for ti data
        unsigned int = POLYNOMIAL_DEGREE );     // maximum degree of polynomial
    Regress(
        const vector<double>&,                  // lambda values
        const vector<double>&,                  // dG/dl
        unsigned int = POLYNOMIAL_DEGREE );     // maximum degree of polynomial
    ~Regress() {};

    void PrintMatrix() const;                   // print out the matrix
    bool GetEstimate( const string&, unsigned int ) const;
    const vector<double>& GetPolynomial( bool = false ) const;
    double DoIntegral( bool = false ) const;
    double DoQuadrature( bool = false ) const;

    const list<stREGRESS>& LoadData(
        const list<stREGRESS>&,                 // structure for ti data
        unsigned int = POLYNOMIAL_DEGREE );     // maximum degree of polynomial
    const list<stREGRESS>& LoadData(
        const vector<double>&,                  // lambda values
        const vector<double>&,                  // dG/dl
        unsigned int = POLYNOMIAL_DEGREE );     // maximum degree of polynomial

private:
    list<stREGRESS> sample;
    vector<double> factor;      // coefficients of the polynomial
    vector<double> matrix;      // matrix for solving linear equations
    unsigned int degree;

    void ClearData();
    void DoPolynomial( unsigned int );  // construct the polynomial using regression
    bool SetPivot( unsigned int );      // apply partial pivoting to the matrix
    void Exchange( double&, double& );  // swap the contents of two variables
    unsigned int SetMatrix();           // construct the matrix
    unsigned int Translate( unsigned int, unsigned int ) const;
};  // class definition for nonlinear regression

#endif  // _REGRESS_H

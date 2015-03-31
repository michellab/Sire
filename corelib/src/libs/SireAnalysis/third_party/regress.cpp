/*
 * regress.cpp
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

#include "regress.h"

#include "SireMaths/errors.h"

/*
 * class constructor
*/
Regress::Regress(
    const list<stREGRESS>& _sample,
    unsigned int _degree )
{
    LoadData( _sample, _degree );
}   // end of class constructor

/*
 * class constructor
*/
Regress::Regress(
    const vector<double>& _x,
    const vector<double>& _y,
    unsigned int _degree )
{
    LoadData( _x, _y, _degree );
}   // end of class constructor

/*
 * reset and initialize essential variables
*/
const list<stREGRESS>& Regress::LoadData(
    const list<stREGRESS>& _sample,
    unsigned int _degree )
{
    stREGRESS unit; ClearData();

    for ( list<stREGRESS>::const_iterator i = _sample.begin(); !( i == _sample.end() ); i++ )
    {
        unit.x = ( *i ).x; unit.y = ( *i ).y; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation with regression and polynomials
    DoPolynomial( _degree ); return( sample );
}   // end of LoadData()

/*
 * reset and initialize essential variables
*/
const list<stREGRESS>& Regress::LoadData(
    const vector<double>& _x,       // lambda values
    const vector<double>& _y,       // dg/dl; free energy difference
    unsigned int _degree )          // degree of polynomial
{
    stREGRESS unit; ClearData();

    for ( unsigned int i = 0; i < _x.size(); ++i )
    {
        unit.x = _x[ i ]; unit.y = _y[ i ]; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation with regression and polynomials
    DoPolynomial( _degree ); return( sample );
}   // end of LoadData()

/*
 * clear all contents
*/
void Regress::ClearData()
{
    sample.clear(); factor.clear(); matrix.clear();
}   // end of ClearData()

/*
 * translate two dimensional coordinate into one dimension
*/
unsigned int Regress::Translate(
    unsigned int _r, unsigned int _c ) const
{
    return( ( _r * factor.size() + _c ) );
}   // end of Translate()

/*
 * swap the contents of two variables
*/
void Regress::Exchange(
    double& _a, double& _b )
{
    double swap = _a; _a = _b; _b = swap;
}   // end of Exchange()

/*
 * set the matrix pivoting elements
 * note: matrix pivoting has been verified to work correctly on december 29, 2007
*/
bool Regress::SetPivot(
    unsigned int _r )       // current diagonal position
{
    for ( unsigned i = ( _r + 1 ); i < factor.size(); ++i )
    {
        if ( !( matrix[ Translate( i, _r ) ] > matrix[ Translate( _r, _r ) ] ) )
        {
            continue;
        }   // search the largest value for pivot

        for ( unsigned int j = 0; j < factor.size(); ++j )
        {
            Exchange( matrix[ Translate( i, j ) ], matrix[ Translate( _r, j ) ] );
        }   // swap the elements in the matrix

        Exchange( factor[ i ], factor[ _r ] );
    }   // perform partial pivoting on the matrix

    return( ( fabs( matrix[ Translate( _r, _r ) ] ) > TOLERANCE_LEVEL ) ? true : false );
}   // end of SetPivot()

/*
 * construct the regression matrix
 *
 * note: the construction of the linear equation matrix has been verified to produce
 * correct results on december 27, 2007
*/
unsigned int Regress::SetMatrix()
{
    vector<double> a; a.clear(); a.resize( ( 2 * degree + 1 ), 0.0 );
    vector<double> b; b.clear(); b.resize( ( degree + 1 ), 0.0 );

    for ( list<stREGRESS>::iterator i = sample.begin(); !( i == sample.end() ); i++ )
    {
        for ( unsigned int j = 0; j < a.size(); ++j )
        {
            a[ j ] += pow( ( *i ).x, static_cast<double>( j ) );
        }   // calculate the coefficients for the x elements

        for ( unsigned int k = 0; k < b.size(); ++k )
        {
            b[ k ] += pow( ( *i ).x, static_cast<double>( k ) ) * ( *i ).y;
        }   // calculate the coefficients for the y elements
    }   // iterate through all elements in the list to construct the matrix

    for ( unsigned int u = 0; u < b.size(); ++u )
    {
        factor.push_back( b[ u ] );

        for ( unsigned int v = 0; v < b.size(); ++v )
        {
            matrix.push_back( a[ ( u + v ) ] );
        }   // copy the row elements
    }   // assign the elements into the matrix data structure

    return( factor.size() );    // return the size of the matrix
}   // end of SetMatrix()

/*
 * print out the matrix
*/
void Regress::PrintMatrix() const
{
    for ( unsigned int i = 0; i < factor.size(); ++i )
    {
        for ( unsigned int j = 0; j < factor.size(); ++j )
        {
            printf( "%8.6f ", matrix[ Translate( i, j ) ] );
        }   // print out each element in the array

        printf( "| %8.6f\n", factor[ i ] );
    }   // print out the matrix
}   // end of PrintMatrix()

/*
 * Gaussian elimination is an algorithm that can be used to determine the solutions of a
 * system of linear equations, to find the rank of a matrix, and to calculate the inverse
 * of an invertible square matrix. Gaussian elimination is named after German
 * mathematician and scientist Carl Friedrich Gauss.
 *
 * Elementary row operations are used throughout the algorithm. The algorithm has two
 * parts, each of which considers the rows of the matrix in order. The first part reduces
 * the matrix to row echelon form while the second reduces the matrix further to reduced
 * row echelon form. The first part alone is sufficient for many applications.
*/
void Regress::DoPolynomial(
    unsigned int _degree )
{
    degree = ( _degree > POLYNOMIAL_DEGREE ) ? POLYNOMIAL_DEGREE : _degree;
    SetMatrix();        // construct the regression matrix
    double ratio; unsigned int u;

    for ( unsigned int s = 0; s < factor.size(); ++s )
    {
        if ( !SetPivot( s ) )
        {
            throw SireMaths::math_error( QObject::tr(
                   "Cannot fit the polynomial as the matrix is singular."), CODELOC );
            // cout << "Error: matrix is singular" << endl; // exit( 1 );
        }   // apply partial pivoting to the matrix and check for singularity

        for ( unsigned int i = ( s + 1 ); i < factor.size(); ++i )
        {
            ratio = matrix[ Translate( i, s ) ] / matrix[ Translate( s, s ) ];
            factor[ i ] -= factor[ s ] * ratio;

            for ( unsigned int j = s; j < factor.size(); ++j )
            {
                matrix[ Translate( i, j ) ] -= ( matrix[ Translate( s, j ) ] * ratio );
            }   // successively remove previous terms
        }   // perform forward elimination on the matrix
    }   // perform gauss elimination with partial pivoting

    for ( unsigned int offset = 0; offset < factor.size(); ++offset )
    {
        u = factor.size() - offset - 1;
        ratio = factor[ u ] / matrix[ Translate( u, u ) ];
        matrix[ Translate( u, u ) ] = 1.0; factor[ u ] = ratio;

        for ( unsigned int v = 0; v < u; ++v )
        {
            factor[ v ] -= ( matrix[ Translate( v, u ) ] * ratio );
            matrix[ Translate( v, u ) ] = 0.0;
        }   // update the solution to the linear equations
    }   // perform backward substitution
}   // end of DoPolynomial()

/*
 * return the newton interpolating polynomial
*/
const vector<double>& Regress::GetPolynomial(
    bool _print ) const
{
    if ( _print )
    {
        printf( "Degree, Coefficients\n" );

        for ( unsigned int i = 0; i < factor.size(); ++i )
        {
            printf( "%6d, %.8f\n", i, factor[ i ] );
        }   // print out the polynomial model
    }   // print out the coefficients of the polynomial, if necessary

    return( factor );
}   // end of GetPolynomial()

/*
 * perform integration on the newton polynomial
*/
double Regress::DoIntegral(
    bool _print ) const
{
    double lower = ( sample.front() ).x;
    double upper = ( sample.back() ).x;
    double area = 0.0;
    double power;

    for ( unsigned int i = 0; i < factor.size(); ++i )
    {
        power = static_cast<double>( i + 1.0 );
        area += ( ( pow( upper, power ) / power ) * factor[ i ] -
            ( pow( lower, power ) / power ) * factor[ i ] );
    }   // perform the integration with given upper and lower bounds

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoIntegral()

/*
 * calculate the area under the curve using quadrature
*/
double Regress::DoQuadrature(
    bool _print ) const
{
    list<stREGRESS>::const_iterator a = sample.begin();
    list<stREGRESS>::const_iterator b = sample.begin(); b++;
    double area = 0.0;

    while ( !( b == sample.end() ) )
    {
        area += ( ( *b ).y + ( *a ).y ) * 0.5 * ( ( *b ).x - ( *a ).x );
        a = b; b++;
    }   // iterate through the entire list

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoQuadrature()

/*
 * get the estimate of f(x) with a given value
*/
bool Regress::GetEstimate(
    const string& _file, unsigned int _step ) const
{
    ofstream ofs( _file.c_str(), ios::trunc );

    if ( ofs.bad() )
    {
        cout << "file " << _file << "cannot be opened" << endl; return( false );
    }   // make sure the file stream has been opened successfully

    char buffer[ 80 ]; double y = 0.0; double x = 0.0;
    double step = 1.0 / static_cast<double>( _step );

    for ( unsigned int s = 0; !( s > _step ); ++s )
    {
        for ( unsigned int i = 0; i < factor.size(); ++i )
        {
            y += ( pow( x, static_cast<double>( i ) ) * factor[ i ] );
        }   // substitute the given value and calculate the estimate

        sprintf( buffer, "%.4f, %.8f", x, y ); ofs << buffer << endl;
        y = 0.0; x += step;
    }   // iterate through the entire interval

    ofs.close(); return( true );
}   // end of GetEstimate()

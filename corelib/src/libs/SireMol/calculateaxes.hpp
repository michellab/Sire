/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIREMOL_CALCULATEAXES_HPP
#define SIREMOL_CALCULATEAXES_HPP

#include "SireMaths/axisset.h"
#include "SireMaths/maths.h"
#include "SireMaths/vector.h"
#include "SireMaths/matrix.h"
#include "SireMaths/quaternion.h"

#include "atom.h"

SIRE_BEGIN_HEADER

namespace SireMol
{

using SireMaths::AxisSet;
using SireMaths::Matrix;
using SireMaths::Vector;
using SireMaths::Matrix4;
using SireMaths::Quaternion;

/** 
This template function is used to calculate the primary, aligned axes of a collection of atoms. The template type 'T' must possess functions that allow index based axes to the atoms in the collection. 

The reason that this is a template function is so that I don't have to write separate code to do this for both an AtomSet and an AtomArray (and perhaps any other atom collection that may be written!)

@author Christopher Woods
*/

template<class T>
SIRE_OUTOFLINE_TEMPLATE
AxisSet SIREMOL_EXPORT calculateAxes(const T &atoms)
{
    if (atoms.count() == 0)
        return AxisSet();
    else if (atoms.count() == 1)
        return AxisSet(Matrix(),atoms[0]);
    
    //First calculate the center of geometry of the atoms - we do this by calculating
    //the average coordinates of the atoms...
    Vector cog;
    
    int sz = atoms.count();
    for (int i=0; i<sz; i++)
    {
        cog += atoms[i].vector();
    }

    cog /= sz;
        
    //Now loop over each atom in the set and calculate the metric tensor for each atom.
    //Sum all of these together, and calculate the principal axes of the resulting matrix
    Matrix tensor_sum = Matrix::zero();
    
    for (int i=0; i<sz; i++)
    {
        Vector ivec = atoms[i].vector() - cog;
        tensor_sum += ivec.metricTensor();
    }
        
    //Get the principal axes of this matrix...
    tensor_sum = tensor_sum.getPrincipalAxes();
    
    return AxisSet(tensor_sum,cog);
}

/** This function calculates the center of mass and the principal inertia tensor for 'atoms'. The template class 'T' must have the same characteristics as that used by the 'calculateAxes' function. */
template<class T>
SIRE_OUTOFLINE_TEMPLATE
AxisSet SIREMOL_EXPORT calculateInertiaAxes(const T &atoms)
{
    if (atoms.count() == 0)
        return AxisSet();
    else if (atoms.count() == 1)
        return AxisSet(Matrix(),atoms[0]);
    
    //First calculate the center of mass of the atoms - we do this by calculating
    //the mass-weighted average coordinates of the atoms...
    Vector com(0.0);
    
    int sz = atoms.count();
    double totalmass = 0.0;
    for (int i=0; i<sz; i++)
    {
        const Atom &atom = atoms[i];
        com += atom.mass() * atom.vector();
        totalmass += atom.mass();
    }

    com /= totalmass;
        
    //Now loop over each atom in the set and calculate the metric tensor for each atom.
    //Sum all of these together, and calculate the principal axes of the resulting matrix
    Matrix inertia_tensor = Matrix::zero();
    
    for (int i=0; i<sz; i++)
    {
        const Atom &atom = atoms[i];
        Vector ivec = atom.vector() - com;
        inertia_tensor += atom.mass() * ivec.metricTensor();
    }
        
    //Get the principal axes of this matrix...
    inertia_tensor = inertia_tensor.getPrincipalAxes();
    
    return AxisSet(inertia_tensor,com);
}    

/** Return the AxisSet needed to align moveatoms against refatoms (minimise the RMSD). 
    This uses the quaternion alignment algorithm presented by Simon K. Kearsley 
    in Acta Cryst. A 45, pp208-210.

    This function aligns atom 'i' in moveatom against atom 'i' in refatom. If there are
    more atoms in moveatoms than in refatoms, then the extra atoms are ignored. If 
    there are less atoms in moveatoms compared to refatoms then the extra atoms in 
    refatoms are ignored. 
 */
template<class S, class T>
SIRE_OUTOFLINE_TEMPLATE
AxisSet SIREMOL_EXPORT calculateAlignmentAxes(const S &moveatoms, const T &refatoms)
{
    int nats0 = moveatoms.count();
    int nats1 = refatoms.count();
    
    //get the maximum number of shared atoms between these two sets.
    int nats = SIRE_MIN(nats0,nats1);
    
    if (nats == 0)
        return AxisSet();
    else if (nats == 1)
        return AxisSet(Matrix(), refatoms[0]);
        
    //calculate the center of geometry of both molecules
    Vector refcog(0.0), movcog(0.0);
    for (int i=0; i<nats; i++)
    {
        refcog += refatoms[i];
        movcog += moveatoms[i];
    }
    
    refcog /= nats;
    movcog /= nats;
    
    //now calculate the 4*4 quaternion matrix...
    /*Matrix4 qmat(0.0);
    
    for (int i=0; i<nats; i++)
    {
        const Atom &refatom = refatoms[i];
        const Atom &movatom = moveatoms[i];
        
        Vector m = refatom.vector() - movatom.vector();
        Vector p = refatom.vector() + movatom.vector();
        
        qmat += Matrix4( m.length2(), p.y()*m.z()-m.y()*p.z(), 
                         m.x()*p.z()-p.x()*m.z(), p.x()*m.y()-m.x()*p.y(),
                         
                         p.y()*m.z()-m.y()*p.z(), p.y()*p.y()+p.z()*p.z()+m.x()*m.x(),
                         m.x()*m.y()-p.x()*p.y(), m.x()*m.z()-p.x()*p.z(),
                         
                         m.x()*p.z()-p.x()*m.z(), m.x()*m.y()-p.x()*p.y(),
                         p.x()*p.x()+p.z()*p.z()*m.y()*m.y(), m.y()*m.z()-p.y()*p.z(),
                         
                         p.x()*m.y()-m.x()*p.y(), m.x()*m.z()-p.x()*p.z(),
                         m.y()*m.z()-p.y()*p.z(), p.x()*p.x()+p.y()*p.y()+m.z()*m.z() );
    }
    
    //diagonalise this matrix
    qmat = qmat.getPrincipalAxes();
    
    //the rotation matrix that maximises the RMSD is the smallest (4th) eigenvector
    Quaternion q = qmat.column0();
    
    //convert this into the rotation matrix
    Matrix mat = q.toMatrix();
    
    return AxisSet(mat, refcog - movcog);
    */
    
    return AxisSet();
}

}

SIRE_END_HEADER

#endif

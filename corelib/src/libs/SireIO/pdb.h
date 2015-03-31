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

#ifndef SIREIO_PDB_H
#define SIREIO_PDB_H

#include "iobase.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class PDB;
}

QDataStream& operator<<(QDataStream&, const SireIO::PDB&);
QDataStream& operator>>(QDataStream&, SireIO::PDB&);

class QTextStream;

namespace SireIO
{

using SireMol::MoleculeView;

/** This class holds all of the sources and default values of the
    properties and parameters used by the PDB reader/writer
    
    @author Christopher Woods
*/
class SIREIO_EXPORT PDBParameters : public IOParametersBase
{
public:
    PDBParameters();
    ~PDBParameters();

    //// Locations of the properties into which to place data

    /** Return the name of the property that will contain all of
        the animation frames of the molecule (if there are any!) 
        
        default == "animation-frames"
    */
    const PropertyName& animationFrames() const
    {
        return animation_property;
    }

    /** Return the name of the property that will contain all of the 
        alternative atom positions as read from the PDB file 
        
        default == "alternative"
    */
    const PropertyName& alternatives() const
    {
        return alternatives_property;
    }
    
    /** Return the name of the property that will contain all of the
        residue insertion codes
        
        default == "icode"
    */
    const PropertyName& iCode() const
    {
        return icode_property;
    }
    
    /** Return the name of the property that will contain all of the 
        (temperature) b-factors
        
        default == "b-factor"
    */
    const PropertyName& bFactor() const
    {
        return bfactor_property;
    }
    
    /** Return the name of the property that will contain all of the
        formal charges of the atoms as read from the PDB file
        
        default == "formal-charge"
    */
    const PropertyName& formalCharge() const
    {
        return formalcharge_property;
    }

    /** Return the name of the property that will contain all of the
        names of the atoms as exactly as they were read from the 
        PDB file
        
        default == "PDB-atom-name"
    */
    const PropertyName& pdbAtomName() const
    {
        return pdbatomnames_property;
    }

    /** Return the name of the property that will contain all of 
        the names of the residues in exactly the format as read from the 
        PDB file
        
        default == "PDB-residue-name"
    */
    const PropertyName& pdbResidueName() const
    {
        return pdbresnames_property;
    }

    /** Return the name of the property that will contain all of 
        the names of the chains in exactly the format as read from the 
        PDB file
        
        default == "PDB-chain-name"
    */
    const PropertyName& pdbChainName() const
    {
        return pdbchainnames_property;
    }

    /** Return the name of the property that will contain all of 
        the names of the segments in exactly the format as read from the 
        PDB file
        
        default == "PDB-segment-name"
    */
    const PropertyName& pdbSegmentName() const
    {
        return pdbsegnames_property;
    }

    ////// Parameters that control the reader or writer
    
    /** The function used to select or skip animation frames from each molecule
    
        source  == "animation-frame-selector"
        default == PropertyBase::none
    */
    const PropertyName& animationFrameSelector() const
    {
        return frame_selector;
    }
    
    /** The function used to process atom names. Must be a StringMangler()
    
        source  == "atom-name-mangler"
        default == TrimString()
    */
    const PropertyName& atomNameMangler() const
    {
        return atomname_mangler;
    }
    
    /** The function used to process residue names. Must be a StringMangler()
    
        source  == "residue-name-mangler"
        default == TrimString()
    */
    const PropertyName& residueNameMangler() const
    {
        return resname_mangler;
    }
    
    /** The function used to process chain names. Must be a StringMangler()
    
        source  == "chain-name-mangler"
        default == TrimString()
    */
    const PropertyName& chainNameMangler() const
    {
        return chainname_mangler;
    }
    
    /** The function used to process segment names. Must be a StringMangler()
    
        source  == "segment-name-mangler"
        default == TrimString()
    */
    const PropertyName& segmentNameMangler() const
    {
        return segname_mangler;
    }

private:
    ///////
    /////// Properties that hold the data of the molecule
    ///////
    
    /** The default name of the animation property */
    static PropertyName animation_property;
    
    /** The default name in which to store the alternative atoms.
        This is also used to store the occupancy */
    static PropertyName alternatives_property;
    
    /** The default name of the residue insertion code property */
    static PropertyName icode_property;
    
    /** The default name of the property containing the temperature
        factor */
    static PropertyName bfactor_property;
    
    /** The property containing the formal charges on the atoms */
    static PropertyName formalcharge_property;
    
    /** The property containing the original PDB atom names */
    static PropertyName pdbatomnames_property;

    /** The property containing the original PDB residue names */
    static PropertyName pdbresnames_property;
    
    /** The property containing the original PDB chain names */
    static PropertyName pdbchainnames_property;
    
    /** The property containing the original PDB segment names */
    static PropertyName pdbsegnames_property;
    
    ///////
    /////// Parameters that control how the reader works
    ///////
    
    /** The function used to decide which animation frames
        to read (or skip) */
    static PropertyName frame_selector;
    
    /** The mangler used for atom names */
    static PropertyName atomname_mangler;
    
    /** The mangler used for residue names */
    static PropertyName resname_mangler;
    
    /** The mangler used for chain names */
    static PropertyName chainname_mangler;
    
    /** The mangler used for segment names */
    static PropertyName segname_mangler;
};

/** This is a IOBase object that has been specialised to read and 
    write PDB format files.

    @author Christopher Woods
*/
class SIREIO_EXPORT PDB : public SireBase::ConcreteProperty<PDB,IOBase>
{
public:
    PDB();
    PDB(const PDB &other);
    
    ~PDB();

    static const char* typeName();

    const char* what() const
    {
        return PDB::typeName();
    }

    static const PDBParameters& parameters()
    {
        return pdbparams;
    }

    PDB& operator=(const PDB &other);
    
    bool operator==(const PDB &other) const;
    bool operator!=(const PDB &other) const;

protected:
    MoleculeGroup readMols(const QByteArray &data,
                           const PropertyMap &map) const;

    QByteArray writeMols(const MoleculeGroup &molgroup,
                         const PropertyMap &map) const;

    QByteArray writeMols(const Molecules &molecules,
                         const PropertyMap &map) const;

    int writeMolecule(QTextStream &ts, const MoleculeView &molview,
                      int atomnum, const PropertyMap &map) const;

private:
    /** All of the default sources and parameters used to 
        control the reading and writing of PDB molecules */
    static PDBParameters pdbparams;
};

}

Q_DECLARE_METATYPE( SireIO::PDB );

SIRE_EXPOSE_CLASS( SireIO::PDBParameters )
SIRE_EXPOSE_CLASS( SireIO::PDB )

SIRE_END_HEADER

#endif

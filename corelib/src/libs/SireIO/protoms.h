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

#ifndef SIREIO_PROTOMS_H
#define SIREIO_PROTOMS_H

#include <QString>
#include <QStringList>

#include "SireBase/properties.h"
#include "SireBase/propertymap.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
class ProtoMS;
}

SIREIO_EXPORT QDataStream& operator<<(QDataStream&, const SireIO::ProtoMS&);
SIREIO_EXPORT QDataStream& operator>>(QDataStream&, SireIO::ProtoMS&);

class QTextStream;

namespace SireMM
{
class TwoAtomFunctions;
class ThreeAtomFunctions;
class FourAtomFunctions;
class CLJNBPairs;
}

namespace SireMol
{
class Molecule;
class Molecules;
class MolEditor;
class Perturbation;
class GeometryPerturbation;
class Connectivity;
class ConnectivityEditor;

typedef SireBase::PropPtr<Perturbation> PerturbationPtr;
typedef SireBase::PropPtr<GeometryPerturbation> GeomPertPtr;
}

namespace SireMove
{
class ZMatrix;
}

namespace SireBase
{
class TempDir;
}

namespace SireIO
{

using SireBase::PropertyName;
using SireBase::PropertyMap;

using SireMol::Molecule;
using SireMol::Molecules;
using SireMol::MolEditor;

using SireMove::ZMatrix;

namespace detail
{
class ProtoMSWorkspace;
}

/** This class holds all of the source and default values of the 
    properties used by the ProtoMS parameter reader
    
    @author Christopher Woods
*/
class SIREIO_EXPORT ProtoMSParameters
{
public:
    ProtoMSParameters();
    ~ProtoMSParameters();
    
    //// Locations of the properties into which to place data

    /** Return the name of the property that will contain the 
        partial atomic charges
        
        default == "charge"
    */
    const PropertyName& charge() const
    {
        return charge_property;
    }

    /** Return the name of the property that will contain the 
        atomic Lennard Jones parameters
        
        default == "LJ"
    */
    const PropertyName& lj() const
    {
        return lj_property;
    }

    /** Return the name of the property that will contain the 
        initial partial atomic charges
        
        default == "initial_charge"
    */
    const PropertyName& initialCharge() const
    {
        return initial_charge_property;
    }

    /** Return the name of the property that will contain the 
        initial atomic Lennard Jones parameters
        
        default == "initial_LJ"
    */
    const PropertyName& initialLJ() const
    {
        return initial_lj_property;
    }

    /** Return the name of the property that will contain the 
        final_partial atomic charges
        
        default == "final_charge"
    */
    const PropertyName& finalCharge() const
    {
        return final_charge_property;
    }

    /** Return the name of the property that will contain the 
        final atomic Lennard Jones parameters
        
        default == "final_LJ"
    */
    const PropertyName& finalLJ() const
    {
        return final_lj_property;
    }

    /** Return the name of the property that will contain the 
        molecular connectivity
        
        default == "connectivity"
    */
    const PropertyName& connectivity() const
    {
        return connectivity_property;
    }

    /** Return the name of the property that will contain the 
        coordinates
        
        default == "coordinates"
    */
    const PropertyName& coordinates() const
    {
        return coords_property;
    }

    /** Return the name of the property that will contain the 
        bond parameters
        
        default == "bond"
    */
    const PropertyName& bond() const
    {
        return bond_property;
    }

    /** Return the name of the property that will contain the 
        angle parameters
        
        default == "angle"
    */
    const PropertyName& angle() const
    {
        return angle_property;
    }

    /** Return the name of the property that will contain the 
        dihedral parameters
        
        default == "dihedral"
    */
    const PropertyName& dihedral() const
    {
        return dihedral_property;
    }

    /** Return the name of the property that will contain the 
        Urey-Bradley parameters
        
        default == "Urey-Bradley"
    */
    const PropertyName& ureyBradley() const
    {
        return ub_property;
    }

    /** Return the name of the property that will contain the 
        z-matrix
        
        default == "zmatrix"
    */
    const PropertyName& zmatrix() const
    {
        return zmatrix_property;
    }

    /** Return the name of the property that will contain the
        non-bonded pairs
        
        default == "intrascale"
    */
    const PropertyName& nonBonded() const
    {
        return nb_property;
    }

    /** Return the name of the property that will contain
        the perturbations property
        
        default == "perturbations"
    */
    const PropertyName& perturbations() const
    {
        return perts_property;
    }

private:
    ///////
    /////// Properties that hold the data of the molecule
    ///////
    
    /** The default name of the partial charge property */
    static PropertyName charge_property;

    /** The default name of the LJ property */
    static PropertyName lj_property;
    
    /** The default name of the initial partial charge property */
    static PropertyName initial_charge_property;

    /** The default name of the initial LJ property */
    static PropertyName initial_lj_property;
    
    /** The default name of the final partial charge property */
    static PropertyName final_charge_property;

    /** The default name of the final LJ property */
    static PropertyName final_lj_property;

    /** The default name of the connectivity property */
    static PropertyName connectivity_property;

    /** The default name of the coordinates property */
    static PropertyName coords_property;

    /** The default name of the bond property */
    static PropertyName bond_property;

    /** The default name of the angle property */
    static PropertyName angle_property;

    /** The default name of the dihedral property */
    static PropertyName dihedral_property;

    /** The default name of the Urey-Bradley property */
    static PropertyName ub_property;

    /** The default name of the zmatrix property */
    static PropertyName zmatrix_property;
    
    /** The default name of the non-bonded property */
    static PropertyName nb_property;
    
    /** The default name of the perturbations property */
    static PropertyName perts_property;
};

/** This class is used to read in ProtoMS parameter files and
    parameterise passed molecules.
 
    @author Christopher Woods
*/
class SIREIO_EXPORT ProtoMS
{

friend SIREIO_EXPORT QDataStream& ::operator<<(QDataStream&, const SireIO::ProtoMS&);
friend SIREIO_EXPORT QDataStream& ::operator>>(QDataStream&, SireIO::ProtoMS&);

public:
    enum { PROTEIN = 1,     // a ProtoMS protein molecule
           SOLUTE  = 2,     // a ProtoMS solute molecule
           SOLVENT = 3  };  // a ProtoMS solvent molecule

    ProtoMS();
    ProtoMS(const QString &protoms);
    
    ~ProtoMS();

    static const char* typeName();

    const char* what() const
    {
        return ProtoMS::typeName();
    }
    
    ProtoMS* clone() const;

    static const ProtoMSParameters& parameters()
    {
        return protoms_parameters;
    }

    void setExecutable(const QString &protoms);

    void addParameterFile(const QString &paramfile);
    
    QStringList parameterFiles() const;

    QString parameterisationCommandFile(const Molecule &molecule,
                                        int type) const;
    
    Molecule parameterise(const Molecule &molecule, int type,
                          const PropertyMap &map = PropertyMap());

    Molecules parameterise(const Molecules &molecules, int type,
                           const PropertyMap &map = PropertyMap());

private:
    QString writeShellFile(const SireBase::TempDir &tempdir,
                           const QString &cmdfile) const;
    QString writeCommandFile(const SireBase::TempDir &tempdir, 
                             const Molecule &molecule, int type) const;
    
    void processZMatrixLine(const QStringList &words, 
                            const Molecule &mol, int type,
                            ZMatrix &zmatrix,
                            detail::ProtoMSWorkspace &workspace) const;
                     
    void processZMatrixPertLine(const QStringList &words, const Molecule &mol, int type,
                                QList<SireMol::GeomPertPtr> &geom_perturbations,
                                const ZMatrix &zmatrix,
                                const PropertyMap &pert_map,
                                detail::ProtoMSWorkspace &workspace) const;
                                          
    void processAtomLine(const QStringList &words,
                         MolEditor &editmol, int type,
                         const QString &charge_property,
                         const QString &lj_property,
                         detail::ProtoMSWorkspace &workspace) const;

    void processAtomPertLine(const QStringList &words, MolEditor &mol, int type,
                             const QString &initial_charge_property,
                             const QString &final_charge_property,
                             const QString &initial_lj_property, 
                             const QString &final_lj_property,
                             detail::ProtoMSWorkspace &workspace) const;
    
    void processBondLine(const QStringList &words,
                         const Molecule &molecule, int type,
                         SireMM::TwoAtomFunctions &bondfuncs,
                         detail::ProtoMSWorkspace &workspace) const;
    
    SireMol::PerturbationPtr 
    processBondPertLine(const QStringList &words,
                        const Molecule &molecule, int type,
                        const PropertyName &bond_property,
                        detail::ProtoMSWorkspace &workspace) const;

    void processConnectLine(const QStringList &words,
                            const Molecule &molecule, int type,
                            SireMol::ConnectivityEditor &connectivity,
                            detail::ProtoMSWorkspace &workspace) const;
    
    void processAngleLine(const QStringList &words,
                          const Molecule &molecule, int type,
                          SireMM::ThreeAtomFunctions &anglefuncs,
                          detail::ProtoMSWorkspace &workspace) const;
    
    SireMol::PerturbationPtr
    processAnglePertLine(const QStringList &words,
                         const Molecule &molecule, int type,
                         const PropertyName &angle_property,
                         detail::ProtoMSWorkspace &workspace) const;
    
    QString processDihedralLine(QTextStream &ts,
                                const QStringList &words,
                                const Molecule &molecule, int type,
                                SireMM::FourAtomFunctions &dihedralfuncs,
                                const PropertyName &dihedral_property,
                                QList<SireMol::PerturbationPtr> &perturbations,
                                detail::ProtoMSWorkspace &workspace) const;

    void processBondDeltaLine(const QStringList &words,
                              const Molecule &molecule, int type,
                              ZMatrix &zmatrix,
                              detail::ProtoMSWorkspace &workspace) const;

    void processAngleDeltaLine(const QStringList &words,
                               const Molecule &molecule, int type,
                               ZMatrix &zmatrix,
                               detail::ProtoMSWorkspace &workspace) const;

    void processDihedralDeltaLine(const QStringList &words,
                                  const Molecule &molecule, int type,
                                  ZMatrix &zmatrix,
                                  detail::ProtoMSWorkspace &workspace) const;
    
    void processUBLine(const QStringList &words,
                       const Molecule &molecule, int type,
                       SireMM::TwoAtomFunctions &ubfuncs,
                       detail::ProtoMSWorkspace &workspace) const;
    
    SireMol::PerturbationPtr 
    processUBPertLine(const QStringList &words,
                      const Molecule &molecule, int type,
                      const PropertyName &ub_property,
                      detail::ProtoMSWorkspace &workspace) const;
    
    void processNBLine(const QStringList &words,
                       const Molecule &molecule, int type,
                       SireMM::CLJNBPairs &nbpairs,
                       detail::ProtoMSWorkspace &workspace) const;
    
    Molecule runProtoMS(const Molecule &molecule, int type,
                        const PropertyMap &map) const;

    /** The default properties used to store the parameters */
    static ProtoMSParameters protoms_parameters;

    /** The list of parameter files that will be used to 
        parameterise the molecules */
    QStringList paramfiles;
    
    /** The full path to the ProtoMS executable that will
        be used to perform the parameterisation */
    QString protoms_exe;
};

}

Q_DECLARE_METATYPE( SireIO::ProtoMS )

SIRE_EXPOSE_CLASS( SireIO::ProtoMS )
SIRE_EXPOSE_CLASS( SireIO::ProtoMSParameters )

SIRE_END_HEADER

#endif

/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2017  Christopher Woods
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

#include "pdbparams.h"

#include "SireMol/mover.hpp"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace SireMM;
using namespace SireMol;
using namespace SireBase;
using namespace SireStream;

static const RegisterMetaType<PDBParams> r_params;

QDataStream SIREMM_EXPORT &operator<<(QDataStream &ds, const PDBParams &params)
{
    writeHeader(ds, r_params, 1);
    
    SharedDataStream sds(ds);
    
    sds << params.molinfo << params.propmap;
    
    return ds;
}

QDataStream SIREMM_EXPORT &operator>>(QDataStream &ds, PDBParams &params)
{
    VersionID v = readHeader(ds, r_params);
    
    if (v == 1)
    {
        SharedDataStream sds(ds);
        
        sds >> params.molinfo >> params.propmap;
    }
    else
        throw version_error(v, "1", r_params, CODELOC);

    return ds;
}

/** Construct an empty set of parameters */
PDBParams::PDBParams() : ConcreteProperty<PDBParams,MoleculeProperty>()
{}

/** Construct so that the parameters for the passed molecule 'molecule' can be created */
PDBParams::PDBParams(const MoleculeView &mol, const PropertyMap &map)
           : ConcreteProperty<PDBParams,MoleculeProperty>()
{
    //get the underlying molecule data
    const auto moldata = mol.data();
    
    //get the parameters property name
    const auto param_name = map["parameters"];
    
    //if possible, start from the existing parameters and update from there
    if (moldata.hasProperty(param_name))
    {
        const Property &param_prop = moldata.property(param_name);
        
        if (param_prop.isA<PDBParams>())
        {
            this->operator=(param_prop.asA<PDBParams>());
            
            if (propmap == map and this->isCompatibleWith(moldata.info()))
            {
                this->_pvt_updateFrom(moldata);
                return;
            }
        }
    }
    
    //otherwise construct this parameter from scratch
    this->operator=(PDBParams());
    
    molinfo = MoleculeInfo(moldata.info());
    propmap = map;
    this->_pvt_createFrom(moldata);
}

/** Constructor when we only have the molecule info */
PDBParams::PDBParams(const MoleculeInfo &info)
           : ConcreteProperty<PDBParams,MoleculeProperty>(), molinfo(info)
{}

/** Constructor when we only have the molecule info */
PDBParams::PDBParams(const MoleculeInfoData &info)
           : ConcreteProperty<PDBParams,MoleculeProperty>(), molinfo(info)
{}

/** Copy constructor */
PDBParams::PDBParams(const PDBParams &other)
           : ConcreteProperty<PDBParams,MoleculeProperty>(other),
             molinfo(other.molinfo), propmap(other.propmap)
{}

/** Destructor */
PDBParams::~PDBParams()
{}

const char* PDBParams::typeName()
{
    return QMetaType::typeName( qMetaTypeId<PDBParams>() );
}

/** Assignment operator */
PDBParams& PDBParams::operator=(const PDBParams &other)
{
    if (this != &other)
    {
        MoleculeProperty::operator=(other);
        molinfo = other.molinfo;
        propmap = other.propmap;
    }
    
    return *this;
}

/** Comparison operator */
bool PDBParams::operator==(const PDBParams &other) const
{
    return molinfo == other.molinfo and
           propmap == other.propmap and
           MoleculeProperty::operator==(other);
}

/** Comparison operator */
bool PDBParams::operator!=(const PDBParams &other) const
{
    return not operator==(other);
}

/** Return the MoleculeInfo that describes the molecule whose parameters
    are stored in this object */
MoleculeInfo PDBParams::info() const
{
    return molinfo;
}

/** Return whether or not these parameters are compatible with the molecule
    whose info is in 'molinfo' */
bool PDBParams::isCompatibleWith(const SireMol::MoleculeInfoData &info) const
{
    //they are only compatible if they are the same molinfo, so have the same UID
    return molinfo.UID() == info.UID();
}

/** Return a string representation of these parameters */
QString PDBParams::toString() const
{
    return QObject::tr("PDBParams::null");
}

/** Set the property map that is used to map data to molecular properties */
void PDBParams::setPropertyMap(const PropertyMap &map)
{
    propmap = map;
}

/** Return the property map used to map data to molecular properties */
const PropertyMap& PDBParams::propertyMap() const
{
    return propmap;
}

/** Update these parameters from the contents of the passed molecule. This
    will only work if these parameters are compatible with this molecule */
void PDBParams::updateFrom(const MoleculeView &molview)
{
    this->assertCompatibleWith(molview);
    this->_pvt_updateFrom(molview.data());
}

/** Validate that the parameters in this object are correct and consistent.
    This returns a list of errors if there are any problems. An empty list
    means that everything is ok. */
QStringList PDBParams::validate() const
{
    return QStringList();
}

/** Construct the data in this object by extracting it from the passed 
    MoleculeData. */
void PDBParams::_pvt_createFrom(const MoleculeData &moldata)
{}

/** Update this set of parameters from the passed object */
void PDBParams::_pvt_updateFrom(const MoleculeData &moldata)
{
    //for the moment we will just create everything from scratch.
    //However, one day we will optimise this and take existing
    //data that doesn't need to be regenerated.
    PropertyMap oldmap = propmap;
    const auto info = molinfo;
    
    this->operator=(PDBParams());

    propmap = oldmap;
    molinfo = info;

    this->_pvt_createFrom(moldata);
}

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

#ifndef SIREFF_FORCEFIELD_H
#define SIREFF_FORCEFIELD_H

#include "ff.h"

#include "SireBase/property.h"

SIRE_BEGIN_HEADER

namespace SireFF
{
class NullFF;
}

SIREFF_EXPORT QDataStream& operator<<(QDataStream&, const SireFF::NullFF&);
SIREFF_EXPORT QDataStream& operator>>(QDataStream&, SireFF::NullFF&);

namespace SireFF
{

/** This is a completely null forcefield. It has zero energy!

    @author Christopher Woods
*/
class SIREFF_EXPORT NullFF : public SireBase::ConcreteProperty<NullFF,FF>
{

friend const NullFF& FF::null();

friend SIREFF_EXPORT QDataStream& ::operator<<(QDataStream&, const NullFF&);
friend SIREFF_EXPORT QDataStream& ::operator>>(QDataStream&, NullFF&);

public:
    NullFF();
    NullFF(const NullFF &other);
    
    ~NullFF();
    
    static const char* typeName();

    const FFComponent& components() const;
    
    NullFF& operator=(const NullFF &other);
    
    bool operator==(const NullFF &other) const;
    bool operator!=(const NullFF &other) const;
    
    const MoleculeGroup& at(MGNum mgnum) const;
    
    QString toString() const;

    bool setProperty(const QString &name, const Property &value);
    const Property& property(const QString &name) const;
    bool containsProperty(const QString &name) const;
    const Properties& properties() const;
    void mustNowRecalculateFromScratch();

    bool needsAccepting() const;
    void accept();

protected:
    void recalculateEnergy();
    void group_add(quint32 i, const MoleculeView &molview,
                   const PropertyMap &map);
    void group_add(quint32 i, const ViewsOfMol &molviews, 
                   const PropertyMap &map);
    void group_add(quint32 i, const Molecules &molecules, 
                   const PropertyMap &map);
    void group_add(quint32 i, const MoleculeGroup &molgroup, 
                   const PropertyMap &map);
    bool group_addIfUnique(quint32 i, const MoleculeView &molview, 
                           const PropertyMap &map);
    ViewsOfMol group_addIfUnique(quint32 i, const ViewsOfMol &molviews, 
                                 const PropertyMap &map);
    QList<ViewsOfMol> group_addIfUnique(quint32 i, const Molecules &molecules, 
                                        const PropertyMap &map);
    QList<ViewsOfMol> group_addIfUnique(quint32 i, const MoleculeGroup &molgroup, 
                                        const PropertyMap &map);

    bool group_remove(quint32 i, const MoleculeView &molview);
    ViewsOfMol group_remove(quint32 i, const ViewsOfMol &molviews);
    QList<ViewsOfMol> group_remove(quint32 i, const Molecules &molecules);
    QList<ViewsOfMol> group_remove(quint32 i, const MoleculeGroup &molgroup);
    
    bool group_removeAll(quint32 i, const MoleculeView &molview);
    ViewsOfMol group_removeAll(quint32 i, const ViewsOfMol &molviews);
    QList<ViewsOfMol> group_removeAll(quint32 i, const Molecules &molecules);
    QList<ViewsOfMol> group_removeAll(quint32 i, const MoleculeGroup &molgroup);

    ViewsOfMol group_remove(quint32 i, MolNum molnum);
    QList<ViewsOfMol> group_remove(quint32 i, const QSet<MolNum> &molnums);

    void group_removeAll(quint32 i);

    bool group_update(quint32 i, const MoleculeData &moldata, bool auto_commit);

    QList<Molecule> group_update(quint32 i, const Molecules &molecules, bool auto_commit);
    QList<Molecule> group_update(quint32 i, const MoleculeGroup &molgroup, bool auto_commit);
    
    bool group_setContents(quint32 i, const MoleculeView &molview, 
                           const PropertyMap &map);
    bool group_setContents(quint32 i, const ViewsOfMol &molviews, 
                           const PropertyMap &map);
    bool group_setContents(quint32 i, const Molecules &molecules, 
                           const PropertyMap &map);
    bool group_setContents(quint32 i, const MoleculeGroup &molgroup, 
                           const PropertyMap &map);

    void _pvt_updateName();

    const MoleculeGroup& getGroup(MGNum mgnum) const;
    
    void getGroups(const QList<MGNum> &mgnums,
                   QVarLengthArray<const MoleculeGroup*,10> &groups) const;

    QHash<MGNum,const MoleculeGroup*> getGroups() const;
    
    void group_setName(quint32 i, const QString &new_name);

    void reindex();

private:
    NullFF(bool);
    
    static NullFF& getSharedNullFF();

    Properties props;
};

typedef SireBase::PropPtr<FF> FFPtr;
typedef FFPtr ForceField;

}

Q_DECLARE_METATYPE(SireFF::NullFF);

SIRE_EXPOSE_CLASS( SireFF::NullFF )

SIRE_EXPOSE_PROPERTY( SireFF::FFPtr, SireFF::FF )

SIRE_END_HEADER

#endif

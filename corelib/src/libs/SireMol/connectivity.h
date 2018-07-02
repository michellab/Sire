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

#ifndef SIREMOL_CONNECTIVITY_H
#define SIREMOL_CONNECTIVITY_H

#include <QHash>
#include <QVector>
#include <QSet>

#include <boost/tuple/tuple.hpp>

#include "bondhunter.h"
#include "molviewproperty.h"
#include "moleculeinfo.h"

#include "SireBase/property.h"
#include "SireBase/shareddatapointer.hpp"

SIRE_BEGIN_HEADER

namespace SireMol
{
class ConnectivityBase;
class Connectivity;
class ConnectivityEditor;
class MoleculeInfo;
}

QDataStream& operator<<(QDataStream&, const SireMol::ConnectivityBase&);
QDataStream& operator>>(QDataStream&, SireMol::ConnectivityBase&);

QDataStream& operator<<(QDataStream&, const SireMol::Connectivity&);
QDataStream& operator>>(QDataStream&, SireMol::Connectivity&);

QDataStream& operator<<(QDataStream&, const SireMol::ConnectivityEditor&);
QDataStream& operator>>(QDataStream&, SireMol::ConnectivityEditor&);

namespace SireMol
{

using boost::tuple;

class AtomIdx;
class ResIdx;

class AtomSelection;

class AtomID;
class ResID;

class BondID;
class AngleID;
class DihedralID;
class ImproperID;

class MoleculeData;
class MoleculeInfoData;

class ConnectivityEditor;

/** The base class of Connectivity and ConnectivityEditor

    @author Christopher Woods
*/
class SIREMOL_EXPORT ConnectivityBase : public MolViewProperty
{

friend QDataStream& ::operator<<(QDataStream&, const SireMol::ConnectivityBase&);
friend QDataStream& ::operator>>(QDataStream&, SireMol::ConnectivityBase&);

public:
    ~ConnectivityBase();

    static const char* typeName()
    {
        return "SireMol::ConnectivityBase";
    }

    QString toString() const;

    MoleculeInfo info() const;

    bool isCompatibleWith(const MoleculeInfoData &molinfo) const;

    bool areConnected(AtomIdx atom0, AtomIdx atom1) const;
    bool areConnected(const AtomID &atom0, const AtomID &atom1) const;

    bool areConnected(ResIdx res0, ResIdx res1) const;
    bool areConnected(const ResID &res0, const ResID &res1) const;

    bool areConnected(CGIdx cg0, CGIdx cg1) const;
    bool areConnected(const CGID &cg0, const CGID &cg1) const;

    bool areBonded(AtomIdx atom0, AtomIdx atom1) const;
    bool areAngled(AtomIdx atom0, AtomIdx atom2) const;
    bool areDihedraled(AtomIdx atom0, AtomIdx atom3) const;

    bool areBonded(const AtomID &atom0, const AtomID &atom1) const;
    bool areAngled(const AtomID &atom0, const AtomID &atom2) const;
    bool areDihedraled(const AtomID &atom0, const AtomID &atom3) const;

    int connectionType(AtomIdx atom0, AtomIdx atom1) const;
    int connectionType(const AtomID &atom0, const AtomID &atom1) const;

    QList<AtomIdx> findPath(AtomIdx atom0, AtomIdx atom1) const;
    QList< QList<AtomIdx> > findPaths(AtomIdx atom0, AtomIdx atom1) const;

    QList<AtomIdx> findPath(const AtomID &atom0, const AtomID &atom1) const;
    QList< QList<AtomIdx> > findPaths(const AtomID &atom0, const AtomID &atom1) const;

    bool inRing(AtomIdx atom0, AtomIdx atom1) const;
    bool inRing(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2) const;
    bool inRing(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3) const;

    bool inRing(const AtomID &atom0, const AtomID &atom1) const;
    bool inRing(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2) const;
    bool inRing(const AtomID &atom0, const AtomID &atom1,
                const AtomID &atom2, const AtomID &atom3) const;

    bool inRing(const BondID &bond) const;
    bool inRing(const AngleID &angle) const;
    bool inRing(const DihedralID &dihedral) const;

    int nConnections() const;
    int nConnections(AtomIdx atomidx) const;
    int nConnections(const AtomID &atomid) const;

    int nConnections(ResIdx residx) const;
    int nConnections(const ResID &resid) const;

    int nConnections(ResIdx res0, ResIdx res1) const;
    int nConnections(const ResID &res0, const ResID &res1) const;

    const QSet<AtomIdx>& connectionsTo(AtomIdx atomidx) const;
    const QSet<AtomIdx>& connectionsTo(const AtomID &atomid) const;

    const QSet<ResIdx>& connectionsTo(ResIdx residx) const;
    const QSet<ResIdx>& connectionsTo(const ResID &resid) const;

    tuple<AtomSelection,AtomSelection> split(AtomIdx atom0, AtomIdx atom1) const;
    tuple<AtomSelection,AtomSelection> split(const AtomID &atom0,
                                             const AtomID &atom1) const;

    tuple<AtomSelection,AtomSelection> split(const BondID &bond) const;

    tuple<AtomSelection,AtomSelection>
    split(AtomIdx atom0, AtomIdx atom1, const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const AtomID &atom0, const AtomID &atom1,
          const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const BondID &bond, const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection> split(AtomIdx atom0, AtomIdx atom1,
                                             AtomIdx atom2) const;
    tuple<AtomSelection,AtomSelection> split(const AtomID &atom0,
                                             const AtomID &atom1,
                                             const AtomID &atom2) const;

    tuple<AtomSelection,AtomSelection> split(const AngleID &angle) const;

    tuple<AtomSelection,AtomSelection>
    split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2,
          const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2,
          const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const AngleID &angle, const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3) const;

    tuple<AtomSelection,AtomSelection>
    split(const AtomID &atom0, const AtomID &atom1,
          const AtomID &atom2, const AtomID &atom3) const;

    tuple<AtomSelection,AtomSelection>
    split(const DihedralID &dihedral) const;

    tuple<AtomSelection,AtomSelection>
    split(AtomIdx atom0, AtomIdx atom1, AtomIdx atom2, AtomIdx atom3,
          const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const AtomID &atom0, const AtomID &atom1,
          const AtomID &atom2, const AtomID &atom3,
          const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const DihedralID &dihedral,
          const AtomSelection &selected_atoms) const;

    tuple<AtomSelection,AtomSelection>
    split(const ImproperID &improper) const;

    tuple<AtomSelection,AtomSelection>
    split(const ImproperID &improper,
          const AtomSelection &selected_atoms) const;

    QList<BondID> getBonds() const;
    QList<BondID> getBonds(const AtomID &atom) const;
    QList<AngleID> getAngles() const;
    QList<AngleID> getAngles(const AtomID &atom0) const;
    QList<AngleID> getAngles(const AtomID &atom0, const AtomID &atom1) const;
    QList<DihedralID> getDihedrals() const;
    QList<DihedralID> getDihedrals(const AtomID &atom0) const;
    QList<DihedralID> getDihedrals(const AtomID &atom0, const AtomID &atom1) const;
    QList<DihedralID> getDihedrals(const AtomID &atom0, const AtomID &atom1, const AtomID &atom2) const;

    QVector< QVector<bool> > getBondMatrix(int order) const;
    QVector< QVector<bool> > getBondMatrix(int start, int end) const;

protected:
    ConnectivityBase();
    ConnectivityBase(const MoleculeInfo &molinfo);
    ConnectivityBase(const MoleculeData &moldata);

    ConnectivityBase(const ConnectivityBase &other);

    ConnectivityBase& operator=(const ConnectivityBase &other);

    bool operator==(const ConnectivityBase &other) const;
    bool operator!=(const ConnectivityBase &other) const;

    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const AtomMatcher &atommatcher) const;
    SireBase::PropertyPtr _pvt_makeCompatibleWith(const MoleculeInfoData &molinfo,
                                                  const QHash<AtomIdx,AtomIdx> &map) const;


    /** The which atoms are connected to which other atoms
        in this molecule */
    QVector< QSet<AtomIdx> > connected_atoms;

    /** Which residues are connected to which other residues */
    QVector< QSet<ResIdx> > connected_res;

    /** The info object that describes the molecule */
    MoleculeInfo minfo;

private:
    const QSet<AtomIdx>& _pvt_connectedTo(AtomIdx atomidx) const;

    QList< QList<AtomIdx> > _pvt_findPaths(AtomIdx cursor, const AtomIdx end_atom,
                                           QSet<AtomIdx> &done) const;

    void traceRoute(AtomIdx start, QSet<AtomIdx> &root,
                    QSet<AtomIdx> &group) const;

    void traceRoute(const AtomSelection &selected_atoms,
                    AtomIdx start, QSet<AtomIdx> &root,
                    QSet<AtomIdx> &group) const;

    tuple<AtomSelection,AtomSelection>
    selectGroups(const QSet<AtomIdx> &group0,
                 const QSet<AtomIdx> &group1) const;
};

/** This class contains the connectivity of the molecule, namely which
atoms are connected to which other atoms. The connectivity is used
to move parts of the molecule (e.g. moving an atom also moves all
of the atoms that it is connected to), and to automatically generate
the internal geometry of the molecule (e.g. to auto-generate
all of the bonds, angles and dihedrals). Note that the connectivity
is not the same as the bonding - the connectivity is used to move
parts of the molecule (e.g. moving an atom should move all of the
atoms it is connected to) and also to auto-generate internal angles
(e.g. auto-generation of bonds, angles and dihedrals)

    @author Christopher Woods

*/
class SIREMOL_EXPORT Connectivity
            : public SireBase::ConcreteProperty<Connectivity,
                                                ConnectivityBase>
{

friend QDataStream& ::operator<<(QDataStream&, const SireMol::Connectivity&);
friend QDataStream& ::operator>>(QDataStream&, SireMol::Connectivity&);

public:
    Connectivity();

    Connectivity(const MoleculeInfo &molinfo);
    Connectivity(const MoleculeData &moldata);

    Connectivity(const MoleculeView &molview,
                 const BondHunter &bondhunter = CovalentBondHunter(),
                 const PropertyMap &map = PropertyMap());

    Connectivity(const ConnectivityEditor &editor);
    Connectivity(const Connectivity &other);

    ~Connectivity();

    static const char* typeName();

    Connectivity& operator=(const Connectivity &other);
    Connectivity& operator=(const ConnectivityEditor &editor);

    bool operator==(const Connectivity &other) const;
    bool operator!=(const Connectivity &other) const;

    ConnectivityEditor edit() const;

protected:
    friend class ConnectivityBase;
    Connectivity(const ConnectivityBase &other);

private:
    void squeeze();
};

/** An editor that can be used to edit a Connectivity object

    @author Christopher Woods
*/
class SIREMOL_EXPORT ConnectivityEditor
        : public SireBase::ConcreteProperty<ConnectivityEditor,
                                            ConnectivityBase>
{

friend QDataStream& ::operator<<(QDataStream&, const SireMol::ConnectivityEditor&);
friend QDataStream& ::operator>>(QDataStream&, SireMol::ConnectivityEditor&);

public:
    ConnectivityEditor();

    ConnectivityEditor(const Connectivity &connectivity);

    ConnectivityEditor(const ConnectivityEditor &other);

    ~ConnectivityEditor();

    static const char* typeName();

    ConnectivityEditor& operator=(const ConnectivityBase &other);

    bool operator==(const ConnectivityEditor &other) const;
    bool operator!=(const ConnectivityEditor &other) const;

    ConnectivityEditor& connect(AtomIdx atom0, AtomIdx atom1);
    ConnectivityEditor& disconnect(AtomIdx atom0, AtomIdx atom1);

    ConnectivityEditor& connect(const AtomID &atom0, const AtomID &atom1);
    ConnectivityEditor& disconnect(const AtomID &atom0, const AtomID &atom1);

    ConnectivityEditor& disconnectAll(AtomIdx atomidx);
    ConnectivityEditor& disconnectAll(ResIdx residx);

    ConnectivityEditor& disconnectAll(const AtomID &atomid);
    ConnectivityEditor& disconnectAll(const ResID &resid);

    ConnectivityEditor& disconnectAll();

    Connectivity commit() const;
};

}

Q_DECLARE_METATYPE(SireMol::Connectivity);
Q_DECLARE_METATYPE(SireMol::ConnectivityEditor);

SIRE_EXPOSE_CLASS( SireMol::ConnectivityBase )
SIRE_EXPOSE_CLASS( SireMol::Connectivity )
SIRE_EXPOSE_CLASS( SireMol::ConnectivityEditor )

SIRE_END_HEADER

#endif

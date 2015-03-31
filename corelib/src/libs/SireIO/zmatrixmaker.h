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

#ifndef SIREIO_ZMATRIXMAKER_H
#define SIREIO_ZMATRIXMAKER_H

#include "iobase.h"

SIRE_BEGIN_HEADER

namespace SireIO
{
  class ZmatrixMaker;
}

QDataStream& operator<<(QDataStream&, const SireIO::ZmatrixMaker&);
QDataStream& operator>>(QDataStream&, SireIO::ZmatrixMaker&);

namespace SireMol
{
  class Molecule;
}


namespace SireIO
{
  /** 
      This internal class is used to store information describing a ZmatrixLineTemplate.
      
      @author Julien Michel
   */
  using SireMol::Molecule;

  class ZmatrixLineTemplate
  {
  public:
    ZmatrixLineTemplate( const QString &atom = " ", const QString &bond = " ", const QString &angle = " ", 
		 const QString &dihedral = " ", const double &bondDelta = 0.0, const double &angleDelta = 0.0, 
		 const double &dihedralDelta = 0.0);
    
    bool has(const QString &atom );
    const QString getAtom();
    const QString getBond();
    const QString getAngle();
    const QString getDihedral();
    const double getBondDelta();
    const double getAngleDelta();
    const double getDihedralDelta();
    void setBondDelta( double bondDelta);
    void setAngleDelta( double angleDelta);
    void setDihedralDelta( double dihedralDelta);
    QString toString();
  private:
    QString atom;
    QString bond;
    QString angle;
    QString dihedral;
    double bondDelta;
    double angleDelta;
    double dihedralDelta;
  };
  /** 
      This internal class is used to store information describing a ZmatrixTemplate, which is 
      an hash of ZmatrixLineTemplate with a name
      
      @author Julien Michel
   */
  class ZmatrixTemplate
  {
  public:
    ZmatrixTemplate( const QString &name= " ");

    QList<ZmatrixLineTemplate> getZmatrix();
    ZmatrixLineTemplate getZmatrixLineTemplate(const QString &atom);
    const QString getName();

    void addZmatrixLineTemplate( const QString &atom, const ZmatrixLineTemplate &zmatline);
    void setBondDelta( const QString &atom, const QString &bond, double bondDelta  );
    void setAngleDelta( const QString &atom, const QString &bond, const QString &angle, double angleDelta);
    void setDihedralDelta( const QString &atom, const QString &bond, 
			   const QString &angle, const QString &dihedral, double dihedralDelta);

  private:
    QString name;
    QHash< QString, ZmatrixLineTemplate > zmatrix;
  };

  /** 
      This internal class is used to store information describing a ZmatrixResidue template. It inherits 
      from ZmatrixTemplate.

      @author Julien Michel
  */
  class ZmatrixResidue : public ZmatrixTemplate
  {
  public:
    ZmatrixResidue( const QString &name= " " );

    void setRotation ( const double rot);
    void setTranslation( const double trans);
    void addChain( const QString &name, const ZmatrixTemplate &chain);
    QStringList getBBatoms();
    void addBBatom( const QString &bbAtom);

    const double getRotation();
    const double getTranslation();
    QList<ZmatrixTemplate> getChains();
    const ZmatrixTemplate getChain(const QString &name);
    
  private:
    double rotate;
    double translate;
    QStringList bbatoms;
    // We use 'first' , 'middle' and 'last'
    QHash< QString, ZmatrixTemplate> backbone;
  };


  /** This class is used to read templates describing how a residue can be moved using zmatrix moves 
   and to create a zmatrix property for residues whose template is available.

   @author Julien Michel
*/
  class SIREIO_EXPORT ZmatrixMaker
  {
    friend QDataStream& ::operator<<(QDataStream&, const SireIO::ZmatrixMaker&);
    friend QDataStream& ::operator>>(QDataStream&, SireIO::ZmatrixMaker&);

  public:
    ZmatrixMaker();
    ~ZmatrixMaker();
    
    /** Read the contents of an input file to create a set of ZmatrixResidues */
    void loadTemplates( const QString &templatesfile);
    /** Add the property 'z-matrix' to molecule */
    Molecule applyTemplates( Molecule &molecule);
    
    static const char* typeName();
    
    const char* what() const
    {
      return ZmatrixMaker::typeName();
    }

  private:
    /** The hash of residues for which a template is available*/
    QHash< QString, ZmatrixResidue > residues;
  };
  
}

Q_DECLARE_METATYPE( SireIO::ZmatrixMaker )

SIRE_EXPOSE_CLASS( SireIO::ZmatrixMaker )

SIRE_END_HEADER

#endif

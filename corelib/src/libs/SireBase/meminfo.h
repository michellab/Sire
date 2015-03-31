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

#ifndef SIREBASE_MEMINFO_H
#define SIREBASE_MEMINFO_H

#include "sireglobal.h"

#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireBase
{

namespace detail
{
class MemInfoPvt;
}

/** This class holds information about the current memory usage
    of the process 
    
    @author Christopher Woods
*/
class SIREBASE_EXPORT MemInfo
{
public:
    MemInfo();
    
    MemInfo(const MemInfo &other);
    
    ~MemInfo();
    
    MemInfo& operator=(const MemInfo &other);
    
    QString toString() const;
    
    quint64 allocatedBytes() const;
    quint64 mMappedBytes() const;
    
    quint64 usedBytes() const;
    
    quint64 totalSystemMemory() const;
    quint64 totalVirtualMemory() const;
    
    quint64 usedSystemMemory() const;
    quint64 usedVirtualMemory() const;
    
    static MemInfo takeMeasurement();
    
    static void startMonitoring(int ms=1000);
    static void startMonitoring(const QString &filename, int ms=1000);
    
    static void stopMonitoring();
    
private:  
    /** PIMPL pointer */
    boost::shared_ptr<detail::MemInfoPvt> d;
};

}

SIRE_EXPOSE_CLASS(SireBase::MemInfo)

SIRE_END_HEADER

#endif

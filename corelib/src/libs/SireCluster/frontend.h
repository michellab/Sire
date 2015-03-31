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

#ifndef SIRECLUSTER_FRONTEND_H
#define SIRECLUSTER_FRONTEND_H

#include "sireglobal.h"

#include <QUuid>
#include <QMutex>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Frontend;
class Backend;

class WorkPacket;

/** This is the base class of all Frontends - a Frontend is an object
    that you can use locally that can control a Backend that is either
    local or remote
    
    @author Christopher Woods
*/
class FrontendBase : public boost::noncopyable
{

friend class Frontend;

public:
    FrontendBase();
    
    virtual ~FrontendBase();
    
    virtual bool isLocal() const=0;
    
    virtual QUuid UID()=0;
    
    virtual void startJob(const WorkPacket &workpacket)=0;
    
    virtual void stopJob()=0;
    virtual void abortJob()=0;
    
    virtual void wait()=0;
    virtual bool wait(int timeout)=0;
    
    virtual float progress()=0;
    virtual WorkPacket interimResult()=0;
    
    virtual WorkPacket result()=0;

private:
    /** Mutex to protect access to this Frontend */
    QMutex datamutex;
};

/** This is the generic holder of a Frontend - a Frontend is an object
    that allows us to communicate with Backend, which may be local or remote
    
    @author Christopher Woods
*/
class Frontend
{
public:
    Frontend();
    Frontend(const boost::shared_ptr<FrontendBase> &ptr);
    
    Frontend(const Backend &backend);
    
    Frontend(const Frontend &other);
    
    ~Frontend();

    Frontend& operator=(const Frontend &other);
    
    bool operator==(const Frontend &other) const;
    bool operator!=(const Frontend &other) const;
    
    static Frontend tryAcquire(const Backend &backend);
    
    bool isLocal() const;
    
    bool isNull() const;
    
    QUuid UID();
    
    void startJob(const WorkPacket &workpacket);
    
    void stopJob();
    void abortJob();
    
    void wait();
    bool wait(int timeout);
    
    float progress();
    WorkPacket interimResult();
    
    WorkPacket result();

private:
    /** Pointer to the private implementation of this class */
    boost::shared_ptr<FrontendBase> d;
};

}

SIRE_END_HEADER

#endif

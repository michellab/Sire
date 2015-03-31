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

#ifndef SIRECLUSTER_PROMISE_H
#define SIRECLUSTER_PROMISE_H

#include <QByteArray>

#include "SireStream/streamdata.hpp"

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Node;
class WorkPacket;

namespace detail
{
class PromisePvt;
}

/** This class provides a handle to the (future) result of
    a piece of work that is being conducted on a node
    in the cluster
    
    @author Christopher Woods
*/
class SIRECLUSTER_EXPORT Promise
{

friend class Node;

public:
    Promise();
    
    Promise(const Promise &other);
    
    ~Promise();
    
    Promise& operator=(const Promise &other);
    
    bool operator==(const Promise &other) const;
    bool operator!=(const Promise &other) const;
    
    void abort();
    void stop();
    
    void wait();
    bool wait(int timeout);
    
    bool isNull() const;
    
    bool isRunning();
    
    bool isError();
    void throwError();
    
    bool wasStopped();
    
    bool wasAborted();
    
    float progress();

    WorkPacket input();
    WorkPacket interimResult();
    WorkPacket result();
    
protected:
    Promise(const Node &node, const WorkPacket &initial_workpacket);
    
    /** Pointer to the private implementation */
    boost::shared_ptr<detail::PromisePvt> d;
};

}

SIRE_EXPOSE_CLASS( SireCluster::Promise )

SIRE_END_HEADER

#endif

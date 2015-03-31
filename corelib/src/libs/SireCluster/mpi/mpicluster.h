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

#ifndef SIRECLUSTER_MPI_MPICLUSTER_H
#define SIRECLUSTER_MPI_MPICLUSTER_H

#ifdef SIRE_USE_MPI

#include "sireglobal.h"

#include <QUuid>
#include <QList>

SIRE_BEGIN_HEADER

namespace SireCluster
{

class Frontend;
class Backend;

namespace MPI
{

class Message;
class Reply;
class P2PComm;

/** This class provides the global interface to all of the
    MPI nodes in the cluster (and, on the root node, the 
    global registry of all nodes available via MPI)
    
    This is a private class which is only available internally
    to SireCluster if MPI is available

    @author Christopher Woods
*/
class MPICluster
{
public:
    static void start();
    static void shutdown();

    static void sync();

    static bool isRunning();

    static void registerBackend(const Backend &backend);

    static P2PComm createP2P(int master_rank, int slave_rank);
    
    static Frontend getFrontend();
    static Frontend getFrontend(const QUuid &uid);
    
    static QList<Frontend> getFrontends(int n);
    
    static QList<QUuid> UIDs();

    static bool hasBackend(const QUuid &uid);

    static int master();
    static int getRank();
    static int getCount();
    static bool isMaster();

    static void send(const Message &message);
    static void received(const Message &message);

    //functions called by MPI messages
    static void registerBackend(int rank, const QUuid &uid);
    static void informedShutdown();

    static Reply getReply(const Message &message);
    
    static void postResult(const QUuid &subject_uid, int sender,
                           const QByteArray &result_data);
                           
    static void postError(const QUuid &subject_uid, int sender,
                          const QByteArray &message_data,
                          const QByteArray &error_data);
};

} // end of namespace MPI

} // end of namespace SireCluster

SIRE_END_HEADER

#endif // SIRE_USE_MPI

#endif


from Sire.Cluster import *

uids = Cluster.UIDs()

print("The number of nodes in your cluster is %d." % len(uids))

print("The UIDs are;")

for uid in uids:
    print(uid)

#get a Nodes object that contains just the current thread
nodes = Nodes()
this_thread = nodes.borrowThisThread()

print(nodes)

node = nodes.getNode()

if (node.isNull()):
    print("Strange - I couldn't get a node!")
    assert( not node.isNull() )

if not node.isLocal():
    print("I'm running on a non-local node!")

#start the job, but don't autodelete the node
promise = node.startJob( WorkTest(0, 3), False )

promise.wait()

result = promise.result()

#start the job, and autodelete the node
promise = node.startJob( WorkTest(0,3) )

result = promise.result()

#add up to 150 more nodes
nodes.addNodes(150)

print(nodes)

promises = []

print("Starting lots of jobs...")
for i in range(0, nodes.count()):
    node = nodes.getNode()
    promises.append( node.startJob(WorkTest(0,10)) )

print("Waiting for them to finish...")
for promise in promises:
    result = promise.result()

print("All done!")

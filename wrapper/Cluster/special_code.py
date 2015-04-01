###############################################
#
# This file contains special code to help
# with the wrapping of SireVol classes
#
#

from pyplusplus.module_builder import call_policies

implicitly_convertible = [ ("SireCluster::WorkPacketBase","SireCluster::WorkPacket") ]

def fixMB(mb):
   mb.add_declaration_code("#include \"SireCluster/workpacket.h\"")


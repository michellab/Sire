

#include "SireMol/molecule.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atommasses.h"
#include "SireMol/cgname.h"
#include "SireMol/resname.h"
#include "SireMol/chainname.h"
#include "SireMol/segname.h"
#include "SireMol/atomname.h"
#include "SireMol/mover.hpp"
#include "SireMol/moleditor.h"

#include "create_test_molecule.h"

using namespace SireMol;

namespace SireSystem
{
    Molecule create_test_molecule()
    {
        auto m = Molecule();
        m = m.edit().add(CGName("1")).add(AtomName("H")).renumber(1).molecule()
                    .add(ChainName("A"))
                    .add(ResName("R")).renumber(1).molecule()
                    .add(SegName("S")).molecule()
                    .atom(AtomName("H"))
                    .reparent(ResName("R")).reparent(SegName("S")).molecule()
                    .commit();

        return m;
    }
}

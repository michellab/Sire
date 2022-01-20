

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
        m = m.edit().add(CGName("0")).add(AtomName("H")).molecule()
                    .add(ChainName("A"))
                    .add(ResName("H")).molecule()
                    .add(SegName("A")).molecule()
                    .atom(AtomName("H"))
                    .reparent(ResName("H")).reparent(SegName("A")).molecule()
                    .commit();

        return m;
    }
}

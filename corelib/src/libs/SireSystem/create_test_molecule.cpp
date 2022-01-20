

#include "SireMol/molecule.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atommasses.h"
#include "SireMol/cgname.h"
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
        m = m.edit().add(CGName("0")).add(AtomName("H")).molecule().commit();

        return m;
    }
}

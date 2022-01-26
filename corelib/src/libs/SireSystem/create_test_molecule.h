#ifndef SIRESYSTEM_CREATE_TEST_MOLECULE_H
#define SIRESYSTEM_CREATE_TEST_MOLECULE_H

#include "sireglobal.h"

SIRE_BEGIN_HEADER

namespace SireMol
{
    class Molecule;
}

namespace SireSystem
{
    SIRESYSTEM_EXPORT SireMol::Molecule create_test_molecule();
}

SIRE_EXPOSE_FUNCTION(SireSystem::create_test_molecule)

SIRE_END_HEADER

#endif

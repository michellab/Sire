#ifndef SIREMOL_MOVER_METAID_H
#define SIREMOL_MOVER_METAID_H

#include "atom.h"
#include "chain.h"
#include "cutgroup.h"
#include "residue.h"
#include "segment.h"

#include "mover.hpp"
#include "selector.hpp"

Q_DECLARE_METATYPE( SireMol::Mover< SireMol::Selector<SireMol::Atom> > )
Q_DECLARE_METATYPE( SireMol::Mover< SireMol::Selector<SireMol::Chain> > )
Q_DECLARE_METATYPE( SireMol::Mover< SireMol::Selector<SireMol::CutGroup> > )
Q_DECLARE_METATYPE( SireMol::Mover< SireMol::Selector<SireMol::Residue> > )
Q_DECLARE_METATYPE( SireMol::Mover< SireMol::Selector<SireMol::Segment> > )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Selector<SireMol::Segment> >, 
                   SireMol::Mover_Selector_Segment_ )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Selector<SireMol::Residue> >, 
                   SireMol::Mover_Selector_Residue_ )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Selector<SireMol::CutGroup> >, 
                   SireMol::Mover_Selector_CutGroup_ )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Selector<SireMol::Chain> >, 
                   SireMol::Mover_Selector_Chain_ )

SIRE_EXPOSE_ALIAS( SireMol::Mover<SireMol::Selector<SireMol::Atom> >, 
                   SireMol::Mover_Selector_Atom_ )

#ifdef SIRE_INSTANTIATE_TEMPLATES

#include "mover.hpp"

template class SireMol::Selector<SireMol::Segment>;
template class SireMol::Mover< SireMol::Selector<SireMol::Segment> >;
template class SireMol::Selector<SireMol::Residue>;
template class SireMol::Mover< SireMol::Selector<SireMol::Residue> >;
template class SireMol::Selector<SireMol::CutGroup>;
template class SireMol::Mover< SireMol::Selector<SireMol::CutGroup> >;
template class SireMol::Selector<SireMol::Chain>;
template class SireMol::Mover< SireMol::Selector<SireMol::Chain> >;
template class SireMol::Selector<SireMol::Atom>;
template class SireMol::Mover< SireMol::Selector<SireMol::Atom> >;

#endif



#endif


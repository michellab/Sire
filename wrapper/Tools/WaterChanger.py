#!/bin/env python
# -*- coding: utf-8 -*-

import Sire.Stream

from Sire.Mol import *
from Sire.MM import *
from Sire.Units import *
from Sire.IO import *
from Sire.Maths import *

global _t4p_template
_t4p_template = None

def getT4PTemplate():
    global _t4p_template

    if _t4p_template is None:
        t = Molecule("T4P")
        t = t.edit() \
             .add( CGName("1") ).molecule() \
             .add( ResName("T4P") ).molecule() \
             .add( AtomName("O00") ) \
               .reparent( CGName("1") ).reparent( ResName("T4P") ).molecule() \
             .add( AtomName("H01") ) \
               .reparent( CGName("1") ).reparent( ResName("T4P") ).molecule() \
             .add( AtomName("H02") ) \
               .reparent( CGName("1") ).reparent( ResName("T4P") ).molecule() \
             .add( AtomName("M03") ) \
               .reparent( CGName("1") ).reparent( ResName("T4P") ).molecule()

        t = t.commit()
        t = t.edit().atom( AtomName("O00") ) \
                    .setProperty("element", Element("O")) \
                    .setProperty("coordinates", Vector(4.293, 30.743, -23.281)) \
                    .setProperty("charge", 0*mod_electron) \
                    .setProperty("LJ", LJParameter(3.15363*angstrom, 0.1550*kcal_per_mol)) \
                    .molecule() \
                    .atom( AtomName("H01") ) \
                    .setProperty("element", Element("H")) \
                    .setProperty("coordinates", Vector(3.678, 30.673, -24.011)) \
                    .setProperty("LJ", LJParameter.dummy()) \
                    .setProperty("charge", 0.520*mod_electron).molecule() \
                    .atom( AtomName("H02") ) \
                    .setProperty("element", Element("H")) \
                    .setProperty("coordinates", Vector(5.061, 30.247, -23.565)) \
                    .setProperty("LJ", LJParameter.dummy()) \
                    .setProperty("charge", 0.520*mod_electron).molecule() \
                    .atom( AtomName("M03") ) \
                    .setProperty("element", Element(0)) \
                    .setProperty("coordinates", Vector(4.312, 30.671, -23.410)) \
                    .setProperty("LJ", LJParameter.dummy()) \
                    .setProperty("charge", -1.040*mod_electron).molecule() \
                    .commit()

        c = Connectivity(t)
        c = c.edit().disconnectAll() \
                    .connect( AtomName("O00"), AtomName("H01") ) \
                    .connect( AtomName("O00"), AtomName("H02") ) \
                    .connect( AtomName("O00"), AtomName("M03") ).commit()

        t = t.edit().setProperty("connectivity", c).commit()

        _t4p_template = t

    return _t4p_template   

def convertTip3PtoTip4P(tip3p):
    """Function creates a copy of the passed TIP3P water molecule
       and returns it converted to a TIP4P water molecule"""

    # first, get the coordinates of the oxygen and two hydrogens
    o_coords = None
    h1_coords = None
    h2_coords = None

    if tip3p.nAtoms() != 3:
        print("Why does this TIP3P water not have three atoms? %d" % tip3p.nAtoms())
        return None

    for atom in tip3p.atoms():
        if atom.property("element").nProtons() == 1:
            # this is one of the hydrogens
            if h1_coords:
                h2_coords = atom.property("coordinates")
            else:
                h1_coords = atom.property("coordinates")
        elif atom.property("element").nProtons() == 8:
            # this is the oxygen
            o_coords = atom.property("coordinates")
        else:
            print("Why is a non oxygen or hydrogen atom in TIP3P? %s" % atom)
            return None

    if o_coords is None or h1_coords is None or h2_coords is None:
        print("Why is this TIP3P molecule missing two hydrogens and an oxygen?")
        return None

    # The M03 atom is 0.15 A along the bond that lies in the plane
    # of the water, in the center of the HOH angle
    m_coords = Vector.generate(0.15, o_coords, 52.26*degrees, h1_coords,  
                               0*degrees, h2_coords)

    tip4p = getT4PTemplate()
    tip4p = tip4p.edit().renumber( tip3p.number() ) \
                 .atom( AtomName("O00") ).setProperty("coordinates", o_coords).molecule() \
                 .atom( AtomName("H01") ).setProperty("coordinates", h1_coords).molecule() \
                 .atom( AtomName("H02") ).setProperty("coordinates", h2_coords).molecule() \
                 .atom( AtomName("M03") ).setProperty("coordinates", m_coords).molecule().commit()

    return tip4p

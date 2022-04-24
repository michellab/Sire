
#include "Python.h"
#include "boost/python.hpp"

#include "SireSystem/system.h"

#include "SireBase/propertymap.h"
#include "SireBase/propertylist.h"

#include "SireIO/moleculeparser.h"

#include "SireMol/moleculegroup.h"
#include "SireMol/bondhunter.h"
#include "SireMol/connectivity.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/element.h"
#include "SireMol/atomelements.h"

#include "Helpers/release_gil_policy.hpp"

#include <QDebug>

using namespace SireBase;
using namespace SireIO;
using namespace SireSystem;
using namespace SireMol;

System load_molecules(const QStringList &files,
                      const PropertyMap &map=PropertyMap())
{
    boost::python::release_gil_policy::release_gil_no_raii();

    try
    {
        auto mols = MoleculeParser::load(files, map);

        // This is an opinionated loader - we must have atom elements
        // and a connectivity defined
        auto grp = MoleculeGroup("all");

        for (const auto &molnum : mols.molNums())
        {
            auto mol = mols[molnum].molecule();

            if (not mol.hasProperty("element"))
            {
                auto editor = mol.edit();

                for (int i=0; i<mol.nAtoms(); ++i)
                {
                    auto atom = editor.atom(AtomIdx(i));
                    atom.setProperty("element",
                                     Element::biologicalElement(atom.name()));
                    editor = atom.molecule();
                }

                mol = editor.commit();
            }

            if (not mol.hasProperty("connectivity"))
            {
                try
                {
                    auto hunter = CovalentBondHunter();
                    mol = mol.edit().setProperty("connectivity",
                                                 hunter(mol) ).commit();
                }
                catch(...)
                {
                    qDebug() << "Failed to auto-generate the connectivity";
                }
            }

            // we now want to break the molecule up into sub-molecules,
            // based on the connectivity
            // NOT IMPLEMENTED YET

            grp.add(mol);
        }

        auto s = System();
        s.setName(mols.name());
        s.add(grp);

        for (const auto &key : mols.propertyKeys())
        {
            s.setProperty(key, mols.property(key));
        }

        s.setProperty("filenames", StringArrayProperty(files));

        boost::python::release_gil_policy::acquire_gil_no_raii();
        return s;
    }
    catch(SireError::exception &e)
    {
        boost::python::release_gil_policy::acquire_gil_no_raii();
        throw SireError::io_error(
            QObject::tr("Cannot load the molecules: %1")
                .arg(e.why()), CODELOC);
    }
}

void register_SireIO_load_function()
{
    boost::python::def("load_molecules",
                       &load_molecules,
                       ( boost::python::arg("filenames"),
                         boost::python::arg("map")=PropertyMap() ),
                       "Load molecules from the passed files.");
}

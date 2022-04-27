
__all__ = ["MultiMolView", "Selector_Molecule_"]


class MultiMolView:
    def __init__(self):
        """Base class of all multi-molecular views"""
        from Sire.Mol import MoleculeGroup
        self._views = MoleculeGroup()

    def _set_from(self, views):
        """Internal function used to set the list of views. All of
           the views should be of the same type!
        """
        from Sire.Mol import MoleculeGroup
        if views.what() != MoleculeGroup.typeName():
            raise TypeError(f"Setting from the wrong type! {views.what()}")

        self._views = views

    def __getitem__(self, idx):
        """Return the item at the specified index. By default,
           we index by molecule
        """
        result = self.molecules(idx)

        if len(result) == 1:
            return self._views.viewAt(0).toUnit()
        else:
            return result

    def molecules(self, idx=None):
        """Return the molecules from this view that match the
           passed index (or all molecules if the index is None)
        """
        m = Selector_Molecule_()

        molnums = []

        if idx is None:
            molnums = self._views.molNums()

        elif type(idx) is int:
            molnums = [self._views.molNumAt(idx)]

        elif type(idx) is str:
            try:
                from Sire.Mol import MolName
                molnums = self._views.map(MolName(idx))
            except KeyError:
                molnums = self._views.search(idx).molNums()

        elif type(idx) is slice:
            raise TypeError("No slicing yet!")

        else:
            for i in idx:
                molnums.append(self._views.molNumAt(i))

        for molnum in molnums:
            m._views.add(self._views.at(molnum).molecule())

        return m

    def molecule(self, idx):
        """Return the molecule that matches the passed index"""
        mols = self.molecules(idx)

        if mols._views.nMolecules() > 1:
            raise KeyError(f"More than one molecule matches {idx}. "
                           f"Number of matches is {len(mols)}.")

        if mols._views.nMolecules() == 0:
            raise KeyError(f"There are no matches for {idx}.")

        return mols._views.moleculeAt(0).molecule()

    def nMolecules(self):
        """Return the number of molecules in this selector"""
        return self._views.nMolecules()

    def nAtoms(self):
        """Return the number of atoms in this selector"""
        return self._views.nAtoms()

    def nResidues(self):
        """Return the number of residues in this selector"""
        return self._views.nResidues()

    def nChains(self):
        """Return the number of chains in this selector"""
        return self._views.nChains();

    def nSegments(self):
        """Return the number of segments in this selector"""
        return self._views.nSegments()


class Selector_Molecule_(MultiMolView):
    def __init__(self, molgroup=None):
        """This is a Selector class that contains multiple molecules"""
        super().__init__()

        if molgroup is not None:
            self._construct_from_group(molgroup)

    def _construct_from_group(self, molgroup):
        """Internal function used to construct this MultiMolView
           from the passed molecule group
        """
        from Sire.Mol import MoleculeGroup
        mols = MoleculeGroup()

        for molnum in molgroup.molNums():
            mols.add(molgroup[molnum].molecule())

        self._set_from(mols)

    def __str__(self):
        n = len(self)

        if n == 0:
            return "Selector<SireMol::Molecule>::empty"

        parts = []

        if n < 10:
            for i in range(0, n):
                parts.append(f"{i}:  {self.__getitem__(i)}")
        else:
            for i in range(0, 5):
                parts.append(f"{i}:  {self.__getitem__(i)}")

            parts.append("...")

            for i in range(n-5, n):
                parts.append(f"{i}:  {self.__getitem__(i)}")

        return "Selector<SireMol::Molecule>( size=%d\n%s\n)" % \
                                                (n, "\n".join(parts))

    def names(self):
        """Return the names of all of the molecules"""
        n = []

        for molnum in self._views.molNums():
            n.append(self._views[molnum].molecule().name())

        return n

    def numbers(self):
        """Return the numbers of all of the molecules"""
        n = []

        for molnum in self._views.molNums():
            n.append(self._views[molnum].molecule().number())

        return n

    def __len__(self):
        return self.nMolecules()

    def count(self):
        return len(self)


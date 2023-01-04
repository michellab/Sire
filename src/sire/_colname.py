
__all__ = ["colname", "colnames"]

_col_funcs = None


def colname(obj, component=None):
    global _col_funcs

    if _col_funcs is None:
        from .mol import Atom, Residue, Chain, Segment, CutGroup, Molecule
        from .mm import Bond, Angle, Dihedral, Improper

        def _colname_atom(atom):
            return f"{atom.name().value()}:{atom.number().value()}"

        def _colname_residue(res):
            return f"{res.name().value()}:{res.number().value()}"

        def _colname_bond(bond):
            return f"{_colname_atom(bond.atom0())}=>{_colname_atom(bond.atom1())}"

        def _colname_angle(ang):
            return f"{_colname_atom(ang.atom0())}<={_colname_atom(ang.atom1())}=>{_colname_atom(ang.atom2())}"

        def _colname_dihedral(dih):
            return f"{_colname_atom(dih.atom0())}<={_colname_atom(dih.atom1())}={_colname_atom(dih.atom2())}=>{_colname_atom(dih.atom3())}"

        def _colname_improper(imp):
            return f"{_colname_atom(imp.atom0())}<={_colname_atom(imp.atom1())}=>{_colname_atom(imp.atom2())}--{_colname_atom(imp.atom3())}"

        def _colname_molecule(mol):
            return f"{mol.name().value()}:{mol.number().value()}"

        _col_funcs = {}
        _col_funcs[Atom] = _colname_atom
        _col_funcs[Residue] = _colname_residue
        _col_funcs[Bond] = _colname_bond
        _col_funcs[Angle] = _colname_angle
        _col_funcs[Dihedral] = _colname_dihedral
        _col_funcs[Improper] = _colname_improper
        _col_funcs[Molecule] = _colname_molecule

    try:
        n = _col_funcs[type(obj)](obj)
    except KeyError:
        n = str(obj)

    if component is None:
        return n
    else:
        return f"{component}({n})"


def colnames(views):
    return [colname(view) for view in views]

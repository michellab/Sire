
from typing import Union as _Union
from typing import List as _List

__all__ = ["load", "save", "create", "smiles", "expand",
           "tutorial_url", "supported_formats"]


tutorial_url = "https://siremol.org/m"

_range = range


def supported_formats():
    """Return a string that describes all of the molecular file formats
       that are supported by Sire
    """
    from Sire.IO import MoleculeParser
    return MoleculeParser.supportedFormats()


def _create_dir(directory):
    import os

    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.isdir(directory):
        raise IOError(f"{directory} is not a directory!")


def _get_gromacs_dir():
    import os

    if "GROMACS_HOME" in os.environ:
        gromacs_dir = os.environ["GROMACS_HOME"]
        if os.path.exists(gromacs_dir) and os.path.isdir(gromacs_dir):
            return gromacs_dir

    import Sire.Config

    gromacs_dir = os.path.join(Sire.Config.share_directory, "gromacs")

    if os.path.exists(gromacs_dir):
        return gromacs_dir

    # it doesn't exist, so we need to download it
    gromacs_tbz2 = os.path.join(Sire.Config.share_directory, "gromacs.tar.bz2")

    if not os.path.exists(gromacs_tbz2):
        try:
            import urllib.request
            urllib.request.urlretrieve(f"{tutorial_url}/gromacs.tar.bz2",
                                       gromacs_tbz2)
        except Exception:
            # we cannot download - just give up
            return None

    if not os.path.exists(gromacs_tbz2):
        return None

    try:
        import tarfile
        t = tarfile.open(gromacs_tbz2, "r|bz2")
        t.extractall(path=Sire.Config.share_directory)
    except Exception:
        return None

    if os.path.exists(gromacs_dir):
        return gromacs_dir
    else:
        return None


def _resolve_path(path, directory):
    import os

    if os.path.exists(path) and os.path.isfile(path):
        if path.endswith(".gz"):
            # unzip the file first
            unzipped = path[0:-3]

            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            _create_dir(directory)
            unzipped = os.path.join(directory, os.path.basename(path))
            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            print(f"Unzipping '{path}'...")

            import gzip
            import shutil
            with gzip.open(path, 'rb') as f_in:
                with open(unzipped, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            return [os.path.abspath(unzipped)]

        elif path.endswith(".bz2"):
            # unzip the file first
            unzipped = path[0:-4]

            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            _create_dir(directory)
            unzipped = os.path.join(directory, os.path.basename(path))
            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            print(f"Unzipping '{path}'...")
            import bz2
            import shutil
            with bz2.open(path, 'rb') as f_in:
                with open(unzipped, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            return [os.path.abspath(unzipped)]

        else:
            return [os.path.abspath(path)]

    if path.startswith("http"):
        # try to download this from the internet
        _create_dir(directory)
        filename = os.path.join(directory, path.split("/")[-1])

        if os.path.exists(filename):
            if os.path.isfile(filename):
                print(f"Using cached download of '{path}'...")
                return _resolve_path(filename, directory=directory)
            else:
                raise IOError(
                    f"Cannot overwrite {filename} as it is an "
                    "existing directory!")

        print(f"Downloading from '{path}'...")

        try:
            import urllib.request
            urllib.request.urlretrieve(path, filename)
        except Exception as e:
            raise IOError(f"Unable to download '{path}': {e}")

        if os.path.exists(filename) and os.path.isfile(filename):
            return _resolve_path(filename, directory=directory)
        else:
            raise IOError(f"Downloaded file does not exist? {filename}")
    else:
        if len(path) == 4:
            # the first character should be a number
            try:
                int(path[0])
                is_code = True
            except Exception:
                is_code = False

            if is_code:
                code = path.lower()
                # https://files.rcsb.org/download/4hhb.pdb.gz
                return _resolve_path(f"https://files.rcsb.org/download/{path}.pdb.gz",
                                     directory=directory)

    # this may be a globbed path
    import glob

    matches = glob.glob(path)

    if len(matches) > 0:
        paths = []
        for match in matches:
            paths += _resolve_path(match, directory=directory)

        return paths

    raise IOError(f"Cannot find file '{path}'")


def expand(base: str, path: _Union[str, _List[str]], *args):
    """Expand the set of paths with the supplied base.

       Args:
        base (str):
            The base to be prepended to all paths

        path (str or list[str]):
            The filename (or names) that will be prepended
            with the base.

        Returns:
            list[str]:
            The list of expanded filenames or URLs

        Examples:
            >>> expand("https://siremol.org/m", "urea.gro", "urea.top")
            ["https://siremol.org/m/urea.gro", "https://siremol.org/n/urea.top"]

            >>> expand("input", ["ala.top", "ala.crd"])
            ["input/ala.top", "input/ala.crd"]
    """
    import os

    if type(path) is not list:
        paths = [path]
    else:
        paths = path

    for arg in args:
        paths.append(arg)

    expanded = []

    for path in paths:
        expanded.append(os.path.join(base, path))

    return expanded


def load(path: _Union[str, _List[str]], *args, **kwargs):
    """Load the molecular system at 'path'. This can be a filename
       of a URL. If it is a URL, then the file will be downloaded
       to the current directory and loaded from there.

       Args:
        path (str or list[str]):
            The filename (or names) or the URL or URLS of the molecular
            system to load. This allows multiple paths to be input
            as some molecular file formats split molecular information
            across multiple files. Multiple paths can also be passed
            as multiple arguments to this function.

        log (dict):
            Optional dictionary that you can pass in that will be populated
            with any error messages or warnings from the parsers as they
            attempt to load in the molecular data. This can be helpful
            in diagnosing why your file wasn't loaded.

        directory (str):
            Optional directory which will be used when creating any
            files (e.g. as a download from a URL or which unzipping files)

       Returns:
            Sire.System.System:
            The molecules that have been loaded are returned as
            a Sire.System.System

       Examples:
            >>> mols = load("caffeine.pdb")

            >>> mols = load(["ala.crd", "ala.top"])

            >>> mols = load("ala.crd", "ala.top")

            >>> mols = load("https://something")

            >>> log = []
            >>> mols = load("caffeine.pdb", log=log)
            Exception
            (look at 'log' to find out what went wrong in detail)
    """
    if type(path) is not list:
        paths = [path]
    else:
        paths = path

    for arg in args:
        paths.append(arg)

    if "log" in kwargs:
        log = kwargs["log"]
    else:
        log = {}

    if "directory" in kwargs:
        directory = kwargs["directory"]
    else:
        directory = "."

    p = []

    for i in range(0, len(paths)):
        # resolve the paths, downloading as needed
        p += _resolve_path(paths[i], directory=directory)

    paths = p

    import Sire.IO
    from Sire.Base import PropertyMap, StringProperty

    map = PropertyMap()
    map.set("GROMACS_PATH", StringProperty(_get_gromacs_dir()))

    mols = Sire.IO.MoleculeParser.load(paths, map=map)

    # This is an opinionated loader - we must have atom elements
    # and a connectivity defined
    from Sire.System import System
    from Sire.Mol import MoleculeGroup

    grp = MoleculeGroup("all")

    for mol in mols:
        c = None

        if not mol.hasProperty("element"):
            from Sire.Mol import Element
            c = mol.cursor()

            for atom in c.atoms():
                atom["element"] = Element.biologicalElement(atom.name().value())
                print(f"guess {atom.name()} is a {atom['element']}")

            mol = c.commit()

        if not mol.hasProperty("connectivity"):
            from Sire.Mol import CovalentBondHunter
            hunter = CovalentBondHunter()

            try:
                connectivity = hunter(mol)
                if c is None:
                    c = mol.cursor()

                c["connectivity"] = connectivity
            except Exception:
                print("Failed to auto-generate the connectivity")
                pass

        if c is not None:
            mol = c.commit()

        # we now want to break the molecule up into sub-molecules,
        # based on the connectivity
        grp.add(mol)

    s = System()
    s.setName(mols.name())
    s.add(grp)

    for key in mols.propertyKeys():
        s.setProperty(key, mols.property(key))

    s.setProperty("filenames", paths)

    return s


def save(molecules, filename: str, format: _Union[str, _List[str]]=None,
         log={}) -> _List[str]:
    """Save the passed molecules to a file called 'filename'. If the format
       is not specified, then the format will be guessed from the
       filename. If the format is specified, and is a list, then multiple
       files will be written, one for each specified format.

       Args:
        molecules (:class:`Sire.System.System`, :class:`Sire.Mol.Molecule`, List[:class:`Sire.Mol.Molecule`] etc.)
            The molecule (or molecules) that should be written to the file.
            This can be anything that can be converted to a :class:`Sire.System.System`,
            i.e. a single :class:`~Sire.Mol.Molecule` (or :class:`~Sire.Mol.MoleculeView`), or a list of
            Molecules (or MoleculeViews)

        filename (str):
            The name of the file to which to write the file. Extensions
            will be automatically added if they are needed to match
            the formats of the file (or files) that are written.

        format (str or list(str)):
            The format (or formats) that should be used to write the
            file (or files). If the format isn't specified, then it
            will be guessed from the extension used for `filename`.
            If this doesn't have an extension, then it will be guessed
            based on the formats used to load the molecule originally.
            If it still isn't available, then PDB will be used.

        log (dict):
            Optional dictionary that you can pass in that will be populated
            with any error messages or warnings from the parsers as they
            attempt to write the molecular data. This can be helpful
            in diagnosing why your file wasn't saved.

       Returns:
            list[str]:
            The absolute paths/name(s) of the files that have been written.

       Examples:
            >>> save(molecules, "molecules.pdb")
            ["/path/to/molecules.pdb"]

            >>> save([mol1, mol2, mol3], "molecules.sdf")
            ["/path/to/molecules.sdf"]

            >>> save(mols, "ala", format=["top", "crd"])
            ["/path/to/ala.top", "/path/to/ala.crd"]

            >>> log = {}
            >>> save(mols, "broken.top", log=log)
            Exception
            (look at `log` to find in detail what went wrong)
    """
    from Sire.IO import MoleculeParser
    from Sire.Base import PropertyMap, StringProperty

    p = PropertyMap()

    if format is not None:
        if type(format) is str:
            format = [format]

        p.set("fileformat", StringProperty(",".join(format)))

    if molecules.what() != "SireSystem::System":
        from Sire.System import System
        from Sire.Mol import MoleculeGroup
        s = System()
        m = MoleculeGroup("all")
        m.add(molecules)
        s.add(m)
        molecules = s

    return MoleculeParser.save(molecules, filename, map=p)


def create(args):
    """Create a Molecule from the passed arguments
    """
    pass


def smiles(args):
    """Create a Molecule from the passed smiles string"""
    pass

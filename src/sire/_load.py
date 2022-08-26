
from typing import Union as _Union
from typing import List as _List

__all__ = ["load", "save", "create", "smiles", "expand",
           "tutorial_url", "load_test_files", "supported_formats"]


tutorial_url = "https://siremol.org/m"

_range = range


def supported_formats():
    """Return a string that describes all of the molecular file formats
       that are supported by Sire
    """
    from .legacy.IO import MoleculeParser

    try:
        return MoleculeParser.supportedFormats()
    except AttributeError:
        return MoleculeParser.supported_formats()


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

    from .config import share_directory

    gromacs_dir = os.path.join(share_directory, "gromacs")

    if os.path.exists(gromacs_dir):
        return gromacs_dir

    # it doesn't exist, so we need to download it
    gromacs_tbz2 = os.path.join(share_directory, "gromacs.tar.bz2")

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
        t.extractall(path=share_directory)
    except Exception:
        return None

    if os.path.exists(gromacs_dir):
        return gromacs_dir
    else:
        return None


def _resolve_path(path, directory, silent=False):
    import os

    if os.path.exists(path) and os.path.isfile(path):
        if path.endswith(".gz"):
            # unzip the file first
            unzipped = path[0:-3]

            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            _create_dir(directory)
            unzipped = os.path.join(directory, os.path.basename(path)[0:-3])
            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            if not silent:
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
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            _create_dir(directory)
            unzipped = os.path.join(directory, os.path.basename(path)[0:-4])
            if os.path.exists(unzipped) and os.path.isfile(unzipped):
                if not silent:
                    print(f"Using cached unzipped file '{unzipped}'...")
                return [os.path.abspath(unzipped)]

            if not silent:
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
                if not silent:
                    print(f"Using cached download of '{path}'...")
                return _resolve_path(filename, directory=directory, silent=silent)
            else:
                raise IOError(
                    f"Cannot overwrite {filename} as it is an "
                    "existing directory!")

        if not silent:
            print(f"Downloading from '{path}'...")

        if not filename.endswith(".bz2"):
            # try the bz2 file first
            try:
                import urllib.request
                urllib.request.urlretrieve(f"{path}.bz2", f"{filename}.bz2")
                have_downloaded_file = True
                filename = f"{filename}.bz2"
            except Exception as e:
                have_downloaded_file = False
        else:
            have_downloaded_file = False

        if not have_downloaded_file:
            try:
                import urllib.request
                urllib.request.urlretrieve(path, filename)
            except Exception as e:
                raise IOError(f"Unable to download '{path}': {e}")

        if os.path.exists(filename) and os.path.isfile(filename):
            return _resolve_path(filename, directory=directory, silent=silent)
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
                                     directory=directory, silent=silent)

    # this may be a globbed path
    import glob

    matches = glob.glob(path)

    if len(matches) > 0:
        paths = []
        for match in matches:
            paths += _resolve_path(match, directory=directory, silent=silent)

        return paths

    raise IOError(f"Cannot find file '{path}'")


def expand(base: str, path: _Union[str, _List[str]], *args, **kwargs):
    """Expand the set of paths with the supplied base.

       Args:
        base (str):
            The base to be prepended to all paths

        path (str or list[str]):
            The filename (or names) that will be prepended
            with the base.

        suffix (str):
            A suffix to attach to all files, e.g. ".bz2"

        Returns:
            list[str]:
            The list of expanded filenames or URLs

        Examples:
            >>> expand("https://siremol.org/m", "urea.gro", "urea.top")
            ["https://siremol.org/m/urea.gro", "https://siremol.org/n/urea.top"]

            >>> expand("input", ["ala.top", "ala.crd"])
            ["input/ala.top", "input/ala.crd"]
    """
    if "suffix" in kwargs:
        suffix = kwargs["suffix"]
    else:
        suffix = None

    if type(path) is not list:
        paths = [path]
    else:
        paths = path

    for arg in args:
        paths.append(arg)

    expanded = []

    if base.startswith("http"):
        join = lambda x, y: f"{x}/{y}"
    else:
        import os
        join = os.path.join

    for path in paths:
        if suffix is None:
            expanded.append(join(base, path))
        else:
            expanded.append(join(base, f"{path}{suffix}"))

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
            sire.system.System:
            The molecules that have been loaded are returned as
            a sire.system.System

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

    if "silent" in kwargs:
        silent = kwargs["silent"]
    else:
        silent = False

    p = []

    for i in range(0, len(paths)):
        # resolve the paths, downloading as needed
        p += _resolve_path(paths[i], directory=directory, silent=silent)

    paths = p

    if len(paths) == 0:
        raise IOError("No valid files specified. Nothing to load?")

    from .io import load_molecules
    from .base import wrap
    return load_molecules(paths, map={"GROMACS_PATH":_get_gromacs_dir()})


def save(molecules, filename: str, format: _Union[str, _List[str]]=None,
         log={}) -> _List[str]:
    """Save the passed molecules to a file called 'filename'. If the format
       is not specified, then the format will be guessed from the
       filename. If the format is specified, and is a list, then multiple
       files will be written, one for each specified format.

       Args:
        molecules (:class:`sire.system.System`, :class:`sire.mol.Molecule`, List[:class:`sire.mol.Molecule`] etc.)
            The molecule (or molecules) that should be written to the file.
            This can be anything that can be converted to a :class:`sire.system.System`,
            i.e. a single :class:`~sire.mol.Molecule` (or :class:`~sire.mol.MoleculeView`), or a list of
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
    from .legacy.IO import MoleculeParser
    from .legacy.Base import PropertyMap, StringProperty

    p = PropertyMap()

    if format is not None:
        if type(format) is str:
            format = [format]

        p.set("fileformat", StringProperty(",".join(format)))

    if molecules.what() != "SireSystem::System":
        from .legacy.System import System
        from .legacy.Mol import MoleculeGroup
        s = System()
        m = MoleculeGroup("all")
        m.add(molecules)
        s.add(m)
        molecules = s

    return MoleculeParser.save(molecules, filename, map=p)


def load_test_files(files: _Union[_List[str], str], *args):
    """Load the passed files that are part of the unit testing
       and return the resulting molecules. This will cache the files
       into a directory called "../cache" so that downloads can be shared
       between tests. You should only need this function if you
       are writing unit tests.

       Args:
        files (str or list[str])
            The list of files to load from the tutorial website. This
            will automatically add on the tutorial URL and compression suffix

       Returns:
        sire.system.System
            The loaded molecules
    """
    if not type(files) is list:
        files = [files]

    for arg in args:
        files.append(arg)

    import os
    d = os.path.abspath(os.path.curdir)

    if (d.endswith("tests")):
        # we are running in the tests directory, so cache downloads here
        cache_dir = os.path.join(d, "cache")
    else:
        d2 = os.path.split(d)[0]
        if d2.endswith("tests"):
            # we are a subdirectory of the parent directory
            cache_dir = os.path.join(d2, "cache")
        else:
            cache_dir = os.path.join(d, "cache")

    files = expand(tutorial_url, files, suffix=".bz2")
    return load(files, directory=cache_dir, silent=True)


def create(args):
    """Create a Molecule from the passed arguments
    """
    pass


def smiles(args):
    """Create a Molecule from the passed smiles string"""
    pass

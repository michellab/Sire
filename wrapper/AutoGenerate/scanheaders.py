##########################################
#
#  This script scans all of the header
#  files in a specified directory to
#  get the list of classes that should
#  be exposed in Python
#

import sys
import os
import re
import pickle
import string

from glob import glob

def getFiles(dir, pattern):
    files = glob("%s/%s" % (dir,pattern))

    trimmed_files = []

    for file in files:
        trimmed_files.append( file[len(dir)+1:] )

    return trimmed_files

def getIncludes(file, lines):

    includes = { "\"%s\"" % file : 1 }

    for line in lines:
        if line.find("CONDITIONAL_INCLUDE") == -1:
            m = re.search(r"#include\s+([\"|<].*[\"|>])", line)

            if m:
                includes[m.groups()[0]] = 1

    ret = list(includes.keys())
    ret.sort()
    return ret

def getDependencies(dir, file):
    """Return the list of header files included by the .cpp or .c file
       that corresponds to 'file'"""

    try:
        #is there a corresponding .cpp file?
        file_cpp = re.sub(r"h(p*)$", "cpp", file)

        lines = open("%s/%s" % (dir,file_cpp), "r").readlines()
        return getIncludes(file, lines)

    except:
        pass

    try:
        #is there a corresponding .c file?
        file_c = re.sub(r"h(p*)$", "c", file)

        lines = open("%s/%s" % (dir,file_c), "r").readlines()
        return getIncludes(file, lines)

    except:
        pass

    return [ "\"%s\"" % file ]

class Properties:
    def __init__(self):
        self._dependencies = {}
        self._properties = []

    def addProperty(self, property, alias):
        self._properties.append( (property,alias) )

    def properties(self):
        return self._properties

    def addDependency(self, headerfile, dir, module_dir):
        deps = getDependencies(dir, headerfile) + getDependencies(module_dir, headerfile)

        for dep in deps:
            self._dependencies[dep] = 1

    def dependencies(self):
        return list(self._dependencies.keys())

skip_metatypes = [ "QVariant",
                   "SireCAS::ExpressionBase",
                   "SireMaths::Rational",
                   "SireCluster::Node",
                   "SireCluster::Nodes" ]

class HeaderInfo:
    def __init__(self, filename, dir, module_dir):
        self._filename = filename
        self._dependencies = getDependencies(dir, filename) + getDependencies(module_dir, filename)
        self._classes = []
        self._functions = []
        self._aliases = {}
        self._properties = []
        self._metatypes = []

    def addClass(self, classname):
        self._classes.append(classname)

    def addFunction(self, func):
        self._functions.append(func)

    def addMetaType(self, classname):
        #don't register some types
        if classname in skip_metatypes:
            return

        self._metatypes.append(classname)

    def addAlias(self, classname, alias):
        self._aliases[classname] = alias

    def addProperty(self, prop, propbase):
        self._properties.append( (prop, propbase) )

    def dependencies(self):
        return self._dependencies

    def classes(self):
        return self._classes

    def functions(self):
        return self._functions

    def metaTypes(self):
        return self._metatypes

    def aliases(self):
        return self._aliases

    def properties(self):
        return self._properties

    def hasProperties(self):
        return len(self._properties) > 0

match_class = r"SIREN*_EXPOSE_CLASS\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_alias = r"SIREN*_EXPOSE_ALIAS\(\s*\n*\s*\(?([<>,\-\s\w\d:]+)\)?\s*\n*,\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_function = r"SIREN*_EXPOSE_FUNCTION\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_property = r"SIREN*_EXPOSE_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_atom_property = r"SIRE_EXPOSE_ATOM_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_res_property = r"SIRE_EXPOSE_RESIDUE_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_cg_property = r"SIRE_EXPOSE_CUTGROUP_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_chain_property = r"SIRE_EXPOSE_CHAIN_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_seg_property = r"SIRE_EXPOSE_SEGMENT_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_bead_property = r"SIRE_EXPOSE_BEAD_PROPERTY\(\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*,\s*\n*\s*\(?([<>,\s\w\d:]+)\)?\s*\n*\s*\)"
match_metatype = r"Q_DECLARE_METATYPE\(\s*\n*\s*\(?([<>.\s\w\d:]+)\)?\s*\n*\s*\)"

db = {}

def add_doc(function, docs, db):
    # now try to extract the function signature
    # - first, remove any derived bases from constructor lines
    function = function.split(" : ")[0]

    match = re.search(r"([\d_\w]+)::([\d_\w]+)\((.*)\)", function)

    if match:
        cls = match.groups()[0].lstrip().rstrip()
        nam = match.groups()[1].lstrip().rstrip()
        targs = match.groups()[2].split(",")

        args = []

        i = 0
        while i < len(targs):
            arg = targs[i].lstrip().rstrip()

            if arg.find("<") != -1:
                nopen = arg.count("<") - arg.count(">")

                #template - see if there are multiple template arguments
                while nopen > 0 and (i+1) < len(targs):
                    next_arg = targs[i+1].lstrip().rstrip()
                    arg += ", " + next_arg
                    i += 1
                    nopen += next_arg.count("<") - next_arg.count(">")

            if len(arg) > 0:
                args.append(arg)

            i += 1

        if not cls in db:
            db[cls] = {}

        if not nam in db[cls]:
            db[cls][nam] = {}

        nargs = len(args)

        if not nargs in db[cls][nam]:
            db[cls][nam][nargs] = []

        db[cls][nam][nargs].append( (args, docs) )

def extract_docs(filename, db):
    read_from = open(filename, "r")

    # concat to one line
    file_str = ''
    for line in read_from.readlines():
        file_str += line

    # break off all code
    file_str = re.sub('{', '\n{', file_str)

    # remove '//' comments
    file_str = re.sub('//.*', '', file_str)

    # split on '\n'
    file_as_list = file_str.splitlines(True)

    # strip everything
    for index in range(len(file_as_list)):
        file_as_list[index] = file_as_list[index].strip()

    doc = ""
    function = ""
    mode = 0   # 0 is nothing, 1 is comment, 2 is function

    # add in newlines where appropriate
    for line in file_as_list:
        line = line.lstrip().rstrip()

        added = False

        if line.startswith('/*'):
            mode = 1
            doc = line
            added = True

        elif line.startswith('{'):
            #print("code... (%s)" % mode)
            if mode == 2:
                # hopefully we have seen a function
                add_doc(function,doc,db)
                doc = ""
                function = ""

            mode = 0

        if line.endswith('*/'):
            # end of comment - look for a function
            mode = 2
            if not added:
                doc += "\n" + line
                added = True
            #print("completed comment\n%s" % doc)

        elif line.endswith(';') or line.endswith('}'):
            #print("line... (%s)" % mode)
            mode = 0

        else:
            if not added:
                if mode == 1:
                    doc += "\n" + line
                elif mode == 2:
                    if len(function) == 0:
                        function = line
                    else:
                        function += " " + line

                    #print("Function | %s" % function)

def scanFiles(dir, module_dir, atom_properties, cg_properties,
                               res_properties, chain_properties, seg_properties,
                               bead_properties):
    """Scan the header files in the passed directory to get information
       about all of the exposed classes, returning a list of all of
       the classes that are being exposed, and placing meta information
       into the directory 'module_dir'"""

    h_files = getFiles(dir, "*.h") + getFiles(module_dir, "*.h")
    hpp_files = getFiles(dir, "*.hpp") + getFiles(module_dir, "*.hpp")
    cpp_files = getFiles(dir, "*.cpp")

    #dictionary mapping files to exposed classes
    active_files = {}

    #the list of exposed classes
    exposed_classes = []

    #the list of classes that have been registered with QMetaType
    meta_classes = []

    #database of all documentation
    doc_db = {}

    #read through each .cpp file, looking for documentation
    for file in cpp_files:
        try:
            extract_docs("%s/%s" % (dir,file), doc_db)
        except Exception as e:
            print("Problem parsing %s | %s" % (file,e))
            pass

    #read each file, looking for SIRE_EXPOSE_FUNCTION or SIRE_EXPOSE_CLASS
    for file in h_files + hpp_files:
        if file.find("sirenglobal.h") != -1:
            continue

        try:
            lines = open("%s/%s" % (dir,file), "r").readlines()
        except:
            lines = open("%s/%s" % (module_dir,file), "r").readlines()

        text = " ".join(lines)

        for m in re.finditer(match_class, text):
            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            active_files[file].addClass(m.groups()[0].strip())
            exposed_classes.append(m.groups()[0].strip())

        for m in re.finditer(match_alias, text):
            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            active_files[file].addClass(m.groups()[0].strip())
            active_files[file].addAlias(m.groups()[0].strip(), m.groups()[1].strip())
            exposed_classes.append(m.groups()[0].strip())

        for m in re.finditer(match_function, text):
            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            active_files[file].addFunction(m.groups()[0].strip())

        for m in re.finditer(match_metatype, text):
            #don't match the 'errors.h' files, as these are wrapped separately
            if file == "errors.h":
                continue

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            active_files[file].addMetaType(m.groups()[0].strip())

        for m in re.finditer(match_property, text):
            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            active_files[file].addProperty(m.groups()[0].strip(), m.groups()[1].strip())

        for m in re.finditer(match_atom_property, text):
            atom_properties.addDependency(file, dir, module_dir)
            atom_properties.addProperty(m.groups()[0].strip(), m.groups()[1].strip())

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            classname = m.groups()[1].split("::")[-1].strip()

            active_files[file].addClass(classname)
            active_files[file].addAlias(classname, classname)
            exposed_classes.append(classname)

        for m in re.finditer(match_cg_property, text):
            cg_properties.addDependency(file, dir, module_dir)
            cg_properties.addProperty(m.groups()[0].strip(), m.groups()[1].strip())

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            classname = m.groups()[1].split("::")[-1].strip()

            active_files[file].addClass(classname)
            active_files[file].addAlias(classname, classname)
            exposed_classes.append(classname)

        for m in re.finditer(match_res_property, text):
            res_properties.addDependency(file, dir, module_dir)
            res_properties.addProperty(m.groups()[0].strip(), m.groups()[1].strip())

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            classname = m.groups()[1].split("::")[-1].strip()

            active_files[file].addClass(classname)
            active_files[file].addAlias(classname, classname)
            exposed_classes.append(classname)

        for m in re.finditer(match_chain_property, text):
            chain_properties.addDependency(file, dir, module_dir)
            chain_properties.addProperty(m.groups()[0].strip(), m.groups()[1].strip())

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            classname = m.groups()[1].split("::")[-1].strip()

            active_files[file].addClass(classname)
            active_files[file].addAlias(classname, classname)
            exposed_classes.append(classname)

        for m in re.finditer(match_seg_property, text):
            seg_properties.addDependency(file, dir, module_dir)
            seg_properties.addProperty(m.groups()[0].strip(), m.groups()[1].strip())

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            classname = m.groups()[1].split("::")[-1].strip()

            active_files[file].addClass(classname)
            active_files[file].addAlias(classname, classname)
            exposed_classes.append(classname)

        for m in re.finditer(match_bead_property, text):
            bead_properties.addDependency(file, dir, module_dir)
            bead_properties.addProperty(m.groups()[0].strip(), m.groups()[1].strip())

            if file not in active_files:
                active_files[file] = HeaderInfo(file, dir, module_dir)

            classname = m.groups()[1].split("::")[-1].strip()

            active_files[file].addClass(classname)
            active_files[file].addAlias(classname, classname)
            exposed_classes.append(classname)

    #now add each active file to a single header file that can be parsed by Py++
    FILE = open("%s/active_headers.h" % module_dir, "w")

    print("#ifndef ACTIVE_HEADERS_H\n" + \
                 "#define ACTIVE_HEADERS_H\n\n" + \
                 "#ifdef GCCXML_PARSE\n", file=FILE)

    files = list(active_files.keys())
    files.sort()

    for file in files:
        print("#include \"%s\"" % file, file=FILE)

    print("\n#endif\n\n#endif", file=FILE)

    FILE.close()

    #now write out the active_files data structure so it can be
    #used by other scripts
    FILE = open("%s/active_headers.data" % module_dir,"wb")
    pickle.dump(active_files, FILE)
    FILE.close()

    #now write out the documentation data so it can be used by other scripts
    FILE = open("%s/docs.data" % module_dir,"wb")
    pickle.dump(doc_db, FILE)
    FILE.close()

    return exposed_classes

if __name__ == "__main__":

    modules = { "Analysis" : "SireAnalysis",
                "Base" : "SireBase",
                "CAS" : "SireCAS",
                "Cluster" : "SireCluster",
                "FF" : "SireFF",
                "ID" : "SireID",
                "IO" : "SireIO",
                "Maths" : "SireMaths",
                "MM" : "SireMM",
                "Mol" : "SireMol",
                "Move" : "SireMove",
                "Search" : "SireSearch",
                "Squire" : "Squire",
                "Stream" : "SireStream",
                "System" : "SireSystem",
                "Units" : "SireUnits",
                "Vol" : "SireVol" }

    if len(sys.argv) < 3:
        print("USAGE: python scanheaders.py input_path output_path")
        sys.exit(-1)

    siredir = sys.argv[1]
    outdir = sys.argv[2]

    exposed_classes = {}

    atom_properties = Properties()
    cg_properties = Properties()
    res_properties = Properties()
    chain_properties = Properties()
    seg_properties = Properties()
    bead_properties = Properties()

    for module in modules:

        try:
            os.makedirs( "%s/%s" % (outdir,module) )
        except:
            pass

        FILE = open("%s/%s/module_info" % (outdir,module), "w")

        print("Module %s" % module, file=FILE)
        print("Source %s" % modules[module], file=FILE)
        print("Root %s" % siredir, file=FILE)

        FILE.close()

        module_classes = scanFiles( "%s/%s" % (siredir,modules[module]),
                                    "%s/%s" % (outdir,module),
                                    atom_properties, cg_properties, res_properties,
                                    chain_properties, seg_properties, bead_properties )

        for clas in module_classes:
            exposed_classes[clas] = 1


    #write the set of exposed classes to a data file to be used
    #by other scripts
    pickle.dump( exposed_classes, open("%s/classdb.data" % outdir, "wb") )

    pickle.dump( atom_properties, open("%s/Mol/atomprops.data" % outdir, "wb") )
    pickle.dump( cg_properties, open("%s/Mol/cgprops.data" % outdir, "wb") )
    pickle.dump( res_properties, open("%s/Mol/resprops.data" % outdir, "wb") )
    pickle.dump( chain_properties, open("%s/Mol/chainprops.data" % outdir, "wb") )
    pickle.dump( seg_properties, open("%s/Mol/segprops.data" % outdir, "wb") )
    pickle.dump( bead_properties, open("%s/Mol/beadprops.data" % outdir, "wb") )

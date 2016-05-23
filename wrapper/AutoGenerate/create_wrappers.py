####################################################
#
# This script uses Py++ to create the Python 
# wrappers for Sire. This script should be run
# in a directory that contains the results
# of scanheaders.py
#
#

import os
import sys
import pickle
import re
import logging

from glob import glob

sys.path.append("../AutoGenerate")
from scanheaders import *
from doxygen import *

from pyplusplus.module_builder import module_builder_t
from pyplusplus.decl_wrappers import calldef_wrapper
from pyplusplus.code_creators import class_t
from pyplusplus.code_creators import algorithm
from pyplusplus.code_creators import free_function_t
from pyplusplus.code_creators import mem_fun_t
from pyplusplus.decl_wrappers import call_policies
import pyplusplus

from pygccxml.declarations.matchers import access_type_matcher_t
from pygccxml import declarations
import pygccxml

all_exposed_classes = {}

def _generate_bases(self, base_creators):
    """This is a new version of the Py++ generate_bases function that only
       adds in bases that are known to be exposed (known via the global list
       of exposed classes)"""

    bases = []
    assert isinstance( self.declaration, declarations.class_t )
    
    for base_desc in self.declaration.bases:
        assert isinstance( base_desc, declarations.hierarchy_info_t )
        if base_desc.access != declarations.ACCESS_TYPES.PUBLIC:
            continue

        try:
            #only include bases that are in the global list
            if (base_desc.related_class.demangled in all_exposed_classes):
                bases.append( algorithm.create_identifier( self, base_desc.related_class.decl_string ) )
        except:
            # doesn't work with CastXML
            demangled = "%s::%s" % (base_desc.related_class.parent.name,
                                    base_desc.related_class.name)

            if (demangled in all_exposed_classes):
                bases.append( algorithm.create_identifier( self, base_desc.related_class.decl_string ) )

    if not bases:
        return None
    
    bases_identifier = algorithm.create_identifier( self, '::boost::python::bases' )
    
    return declarations.templates.join( bases_identifier, bases )

class_t._generate_bases = _generate_bases    

####
#### Override the free_function functions so that we fix a compile bug using xlC on AIX
#### Overloaded function signatures output by Py++ look like this;
####
####  typedef void (*my_function_type)( args );
####  def( "my_function", my_function_type( &my_function ) );
####
####  This breaks when there are multiple overload of "my_function" as the xlC compiler
####  fails with "The call does not match any parameter list for "bp::def"" errors
####
####  The solution is for Py++ to create a variable of type my_function_type and pass this to def, e.g.
####
####  typedef void (*my_function_type)( args );
####  my_function_type my_function_value( &my_function );
####
####  def( "my_function", my_function_value );
####
####  This compiles property using xlC. The below code changes free_function_t and
####  mem_function_t to create the xlC compatible code, rather than the original Py++ code
####
def _create_function_type_alias_code( self, exported_class_alias=None  ):
    f_type = self.declaration.function_type()
    falias = self.function_type_alias
    fname = declarations.full_name( self.declaration, with_defaults=False )
    fvalue = re.sub("_type$", "_value", falias )

    return "typedef %s;\n%s %s( &%s );" % (f_type.create_typedef( falias, with_defaults=False ),
                                           falias, fvalue, fname)

free_function_t.create_function_type_alias_code = _create_function_type_alias_code
mem_fun_t.create_function_type_alias_code = _create_function_type_alias_code

def _create_function_ref_code(self, use_function_alias=False):
    fname = declarations.full_name( self.declaration, with_defaults=False )
    if use_function_alias:
        falias = self.function_type_alias
        fvalue = re.sub("_type$", "_value", falias)
        return fvalue
    elif self.declaration.create_with_signature:
        return '(%s)( &%s )' % ( self.declaration.function_type().partial_decl_string, fname )
    else:
        return '&%s' % fname

free_function_t.create_function_ref_code = _create_function_ref_code
mem_fun_t.create_function_ref_code = _create_function_ref_code

#fix broken "operators" function
def operators( self, name=None, symbol=None, function=None, return_type=None, arg_types=None, decl_type=None, header_dir=None, header_file=None, recursive=None ):
    """Please see L{decl_wrappers.scopedef_t} class documentation"""
    return self.global_ns.operators( name=name
                                     , symbol=symbol
                                     , function=function
                                     , return_type=return_type
                                     , arg_types=arg_types
                                     , header_dir=header_dir
                                     , header_file=header_file
                                     , recursive=recursive )

module_builder_t.operators = operators

def has_datastream_operators(mb, c):
   """Return whether or not the class has QDataStream streaming operators"""

   try:
       d = mb.operators(arg_types=["::QDataStream &","%s &" % c.decl_string])
       return len(d) > 0

   except:
       return False

def has_function(c, funcname):
   """Recursively move through this class and its bases to find
      if it has a function called 'funcname'"""
   
   try:
       c.decl(funcname)
       return True
   except:
       
       for base in c.bases:
           if has_function(base.related_class, funcname):
               return True
       
       return False

def find_class(mb, classname):
   for clas in mb.classes():
       if str(clas).find("%s [class]" % classname) != -1 or \
          str(clas).find("%s [struct]" % classname) != -1 or \
          str(clas).find("%s [typedef]" % classname) != -1:
           return clas
       else:
           for alias in clas.aliases:
               if str(alias).find("%s [class]" % classname) != -1 or \
                  str(alias).find("%s [struct]" % classname) != -1 or \
                  str(alias).find("%s [typedef]" % classname) != -1:
                   print("Found %s as alias %s of %s" % (classname, clas, alias))
                   return clas

   print("Cannot find the class %s" % classname)
   raise "Cannot find the class %s" % classname

def export_function(mb, function, includes):
   """Do all the work necessary to allow the function 'function'
      to be exported, adding the header files in 'includes'
      to the generated C++"""

   name = function.split("::")[-1]
   root = "::".join(function.split("::")[0:-1]) + "::"

   try:
       for f in mb.free_functions(name):
           demangled = f.demangled
       
           if demangled:
               if demangled.find(root) != -1:
                   f.include()

                   for include in includes:
                       f.add_declaration_code("#include %s" % include)
   except:
       for f in mb.free_functions(name):
           # demangled doesn't work with CastXML      
           if root.find(str(f.parent.name)) != -1:
               f.include()

               for include in includes:
                   f.add_declaration_code("#include %s" % include)

def has_clone_function(t):
    c = None

    try:
        fullname = " ".join(str(t.base).split(" ")[0:-1])
        c = find_class(mb, fullname)
    except:
        print("WARNING!!! Couldn't find the class for %s" % (t))
        return False

    try:
        c.mem_funs("clone")
        return True
    except:
        return False


def is_Index_T_(c):
    for base in c.bases:
        if str(base.related_class).find("SireID::Index_T_") != -1:
            return True

    return False


def fix_Index_T_(c):
   """The Index_T_<T> classes need extra work to wrap... This function does that work"""
   try:
       c.operators("==").exclude()
   except:
       pass

has_copy_function = {}


def export_class(mb, classname, aliases, includes, special_code, auto_str_function=True):
   """Do all the work necessary to allow the class called 'classname'
      to be exported, using the supplied aliases, and using the 
      supplied special code, and adding the header files in 'includes'
      to the generated C++"""

   #find the class in the declarations   
   c = find_class(mb, classname)
   
   #include the class in the wrapper
   c.include()

   #write out all function signatures
   c.calldefs().create_with_signature = True
   c.always_expose_using_scope = True
   
   #add the extra include files for this class
   for include_file in includes:
       c.add_declaration_code("#include %s" % include_file)

   #ensure that the list of bases includes *all* bases,
   # - this is to fix problems with typeerror being
   #   thrown for derived types
   c.bases = c.recursive_bases

   #exclude any "clone" functions
   try:
       c.decls( "clone" ).exclude()
   except:
       pass

   #now replace copy_const_reference with clone_const_reference for suitable classes
   funs = []

   try:
       #all copy_const_reference call policies with clone_const_reference
       #funs = c.mem_funs( lambda f: declarations.is_reference( f.return_type ) )
       funs = c.mem_funs( lambda f: f.return_type.decl_string.endswith("&") )
   except:
       pass

   print("ALL FUNCTIONS: %s" % funs)

   for f in funs:
       print("%s : %s : %s" % (f, f.return_type, has_clone_function(f.return_type)))
       if has_clone_function(f.return_type):
           f.call_policies = call_policies.custom_call_policies( \
                 "bp::return_value_policy<bp::clone_const_reference>", \
                 "Helpers/clone_const_reference.hpp" )

   #also add any operator[] or operator() functions
   try:
       #funs = c.operators( lambda f: declarations.is_reference( f.return_type ) )
       funs = c.operators( lambda f: f.return_type.decl_string.endswith("&") )
   except:
       pass

   for f in funs:
       if (str(f).find("[]") != -1) or (str(f).find("()") != -1):
           if has_clone_function(f.return_type):
               f.call_policies = call_policies.custom_call_policies( \
                   "bp::return_value_policy<bp::clone_const_reference>", \
                   "Helpers/clone_const_reference.hpp" )

   #remove any declarations that return a pointer to something
   #(special code is needed in these cases!)
   for decl in c.decls():
       try:
           if str(decl.return_type) != "char const *":
               rt = str(decl.return_type)
               if rt.endswith("*") or \
                  rt.endswith("::iterator") or rt.endswith("::const_iterator"):
                   decl.exclude()
       except:
           pass

   try:
       for o in c.operators():
           if o.call_policies is None:
               o.exclude()
   except:
       print("Class %s has no operators" % classname)

   #run any class specific code
   if (classname in special_code):
     print("Running special code for %s" % classname)
     special_code[classname](c)
   else:
     print("No special code needed for %s" % classname)

   #if this is a noncopyable class then remove all constructors!
   if c.noncopyable:
      c.constructors().exclude()
   else:
      #if there is a copy-constructor then ensure that
      #it is exposed!
      decls = c.decls()
      
      made_copy_function = False

      for decl in decls:
          if made_copy_function:
              break

          try:
              if decl.is_copy_constructor:
                  #create a __copy__ function
                  class_name = re.sub(r"\s\[class\]","",str(c))
                  class_name = re.sub(r"\s\[struct\]","",class_name)
                  
                  if not (class_name in has_copy_function):
                      has_copy_function[class_name] = True

                      print("Creating a copy function for class %s" % class_name)
                      made_copy_function = True

                      c.add_declaration_code( \
                           "%s __copy__(const %s &other){ return %s(other); }" \
                              % (class_name, class_name, class_name) )
                       
                      c.add_registration_code( "def( \"__copy__\", &__copy__)" )
                  
                      c.add_registration_code( "def( \"__deepcopy__\", &__copy__)" )
                      c.add_registration_code( "def( \"clone\", &__copy__)" )

                      #only do this once for the class
                      break
                  
          except AttributeError:
              pass

   #If this is an Index_T_ class then fix the operators
   if is_Index_T_(c):
      fix_Index_T_(c)

   #if this class can be streamed to a QDataStream then add
   #streaming operators
   if has_datastream_operators(mb,c):
       c.add_declaration_code( "#include \"Qt/qdatastream.hpp\"" )

       c.add_registration_code(
            """def( \"__rlshift__\", &__rlshift__QDataStream< %s >,
                    bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )""" % c.decl_string )
       c.add_registration_code(
            """def( \"__rrshift__\", &__rrshift__QDataStream< %s >,
                    bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() )""" % c.decl_string )
           

   #is there a "toString" function for this class?
   if auto_str_function:
       if has_function(c, "toString"):
           #there is a .toString() member function - we can thus use the 
           #templer __str__ function provided in the Helpers directory
           c.add_declaration_code( "#include \"Helpers/str.hpp\"" )
       
           c.add_registration_code( "def( \"__str__\", &__str__< %s > )" % c.decl_string )
           c.add_registration_code( "def( \"__repr__\", &__str__< %s > )" % c.decl_string )   

       else:
           #there is no .toString() function
           # - instead create a new __str__ that just returns a pretty form
           #   of the name of the class
           name = str(c.decl_string)
       
           if name.startswith("::"):
               name = name[2:]
       
           c.add_declaration_code( "const char* pvt_get_name(const %s&){ return \"%s\";}" % (name,name) )
       
           c.add_registration_code("def( \"__str__\", &pvt_get_name)")
           c.add_registration_code("def( \"__repr__\", &pvt_get_name)")
           
   #is there a "count" or "size" function for this class?
   if has_function(c, "size"):
       c.add_declaration_code( "#include \"Helpers/len.hpp\"" )
       c.add_registration_code("def( \"__len__\", &__len_size< %s > )" % c.decl_string )
   elif has_function(c, "count"):
       c.add_declaration_code( "#include \"Helpers/len.hpp\"" )
       c.add_registration_code("def( \"__len__\", &__len_count< %s > )" % c.decl_string )

   #is there a python-style getitem function?
   if has_function(c, "getitem"):
       c.add_registration_code("def( \"__getitem__\", &%s::getitem )" % c.decl_string )

   #is there a .hash() function?
   if has_function(c, "hash"):
       c.add_registration_code("def( \"__hash__\", &%s::hash )" % c.decl_string )

   #provide an alias for this class
   if (classname in aliases):
      c.alias = " ".join( aliases[classname].split("::")[1:] )

def register_implicit_conversions(mb, implicitly_convertible):
    """This function sets the wrapper generator to use only the implicit conversions
       that have been specifically specified in 'implicitly_convertible'"""

    #remove all existing implicit conversions
    try:
        mb.constructors().allow_implicit_conversion = False
    except:
        pass

    try:
        mb.casting_operators().exclude()
    except:
        pass    

    #add our manual implicit conversions to the declaration section
    for conversion in implicitly_convertible:
       mb.add_registration_code("bp::implicitly_convertible< %s, %s >();" % conversion)

def write_wrappers(mb, module, huge_classes, header_files = [] ):
   """This function performs the actual work of writing the wrappers to disk"""
                   
   #make sure that the protected and private member functions and 
   #data aren't wrapped
   try:
       mb.calldefs( access_type_matcher_t( 'protected' ) ).exclude()
   except:
       pass

   try:
       mb.calldefs( access_type_matcher_t( 'private' ) ).exclude()
   except:
       pass
 
   #build a code creator - this must be done after the above, as
   #otherwise our modifications above won't take effect
   mb.build_code_creator( module_name="_%s" % module,
                          doc_extractor=doxygen_doc_extractor() )

   #get rid of the standard headers
   mb.code_creator.replace_included_headers( header_files )

   #give each piece of code the GPL license header
   mb.code_creator.license = "// (C) Christopher Woods, GPL >= 2 License\n"

   #use local directory paths
   mb.code_creator.user_defined_directories.append(".")

   mb.split_module( ".", huge_classes )

def needPropertyWrappers(active_headers):
   for header in active_headers:
       if active_headers[header].hasProperties():
           return True

   return False

def writePropertyWrappers(mb, sourcedir, active_headers):
   """This function writes the property wrappers that are required for this module"""
   
   #are there any wrappers?
   if not needPropertyWrappers(active_headers):
       return

   #create the files
   FILE = open("%s_properties.h" % sourcedir, "w")

   print("#ifndef %s_PROPERTIES_H" % sourcedir, file=FILE)
   print("#define %s_PROPERTIES_H\n" % sourcedir, file=FILE)
   print("void register_%s_properties();\n" % sourcedir, file=FILE)
   print("#endif", file=FILE)

   FILE.close()

   FILE = open("%s_properties.cpp" % sourcedir, "w")

   print("#include <Python.h>", file=FILE)
   print("#include <boost/python.hpp>\n", file=FILE)
   print("#include \"Base/convertproperty.hpp\"", file=FILE)
   print("#include \"%s_properties.h\"\n" % sourcedir, file=FILE)

   for header in active_headers:
       active_header = active_headers[header]
       if active_header.hasProperties():
           for dependency in active_header.dependencies():
               print("#include %s" % dependency, file=FILE)

   print("void register_%s_properties()" % sourcedir, file=FILE)
   print("{", file=FILE)
   
   for header in active_headers:
       active_header = active_headers[header]
       if active_header.hasProperties():
           for property in active_header.properties():
               print("    register_property_container< %s, %s >();" % (property[0], property[1]), file=FILE)

   print("}", file=FILE)

   FILE.close()

   mb.add_declaration_code("#include \"%s_properties.h\"" % sourcedir)
   mb.add_registration_code("register_%s_properties();" % sourcedir)

def fixMB(mb):
   pass

if __name__ == "__main__":

    #read in the information about this module
    lines = open("module_info", "r").readlines()

    module = lines[0].split()[1]
    sourcedir = lines[1].split()[1]
    rootdir = lines[2].split()[1]

    #load up the dictionary of all exposed classes
    all_exposed_classes = pickle.load( open("../classdb.data", "rb") )
    
    #load up the active headers object
    active_headers = pickle.load( open("active_headers.data", "rb") )

    #get the special code, big classes and implicit conversions
    implicitly_convertible = []
    special_code = {}
    huge_classes = []

    if os.path.exists("special_code.py"):
        sys.path.append(".")
        from special_code import *
        
    sire_include_dirs = [ rootdir, "%s/%s" % (rootdir,sourcedir) ]

    # All of the headers must be installed in the pkgs/sire-*** directory
    # - lets locate this directory
    dir = glob( "%s/../pkgs/sire-*" % os.path.dirname(sys.executable) )

    if len(dir) == 0:
        print("Cannot find the Sire directory. Please use the python that comes with "
              "the anaconda/miniconda python, and that you have already installed the "
              "Sire corelib into the anaconda/miniconda distribution.")
        sys.exit(-1)

    qtdir = "%s/bundled/include" % os.path.abspath(dir[0])
    boostdir = qtdir
    gsldir = qtdir
    openmm_include_dir = "%s/../include" % os.path.dirname(sys.executable)

    need_input = False 

    if (qtdir is None):
        print("You must set the environmental variable QTDIR to the location " + \
              "of the Qt4 header files")
        need_input = True

    if (boostdir is None):
        print("You must set the environmental variable BOOSTDIR to the location " + \
              "of the boost header files")
        need_input = True

    if (gsldir is None):
        print("You must set the environmental variable GSLDIR to the location " + \
              "of the GSL header files")
        need_input = True

    if (need_input):
        print("Cannot continue as I don't know where the header files are")
        sys.exit(-1)

    qt_include_dirs = []

    qt_include_dirs = [ qtdir, "%s/QtCore" % qtdir ]
    boost_include_dirs = [ boostdir ]
    gsl_include_dirs = [ gsldir ]

    generator_path, generator_name = pygccxml.utils.find_xml_generator()

    print("%s | %s" % (generator_path, generator_name))

    if openmm_include_dir is not None:
        if os.path.exists("%s/OpenMM.h" % openmm_include_dir):
            print("Generating wrappers including OpenMM from %s" % openmm_include_dir)
            openmm_include_dirs = [openmm_include_dir]
        else:
            print("Cannot find %s/OpenMM.h - disabling generation of OpenMM wrappers." % openmm_include_dir)
            openmm_include_dirs = None

    if os.getenv("VERBOSE"):
        pygccxml.utils.loggers.cxx_parser.setLevel(logging.DEBUG)

    if openmm_include_dirs is None:
        #construct a module builder that will build all of the wrappers for this module
        xml_generator_config = pygccxml.parser.xml_generator_configuration_t( 
                                xml_generator_path=generator_path,
                                xml_generator=generator_name,
                                compiler="gcc",
                                cflags = "-m64 -fPIC",
                                include_paths = sire_include_dirs + qt_include_dirs +
                                           boost_include_dirs + gsl_include_dirs,
                                define_symbols = ["GCCXML_PARSE", "__PIE__",
                                                  "SIRE_SKIP_INLINE_FUNCTIONS",
                                                  "SIREN_SKIP_INLINE_FUNCTIONS",
                                                  "SIRE_INSTANTIATE_TEMPLATES",
                                                  "SIREN_INSTANTIATE_TEMPLATES"]
                         )

        mb = module_builder_t( files = [ "active_headers.h" ],
                               gccxml_config=xml_generator_config )
    else:
        #construct a module builder that will build all of the wrappers for this module
        xml_generator_config = pygccxml.parser.xml_generator_configuration_t(
                                xml_generator_path=generator_path,
                                xml_generator=generator_name,
                                compiler="gcc",
                                cflags = "-m64 -fPIC",
                                include_paths = sire_include_dirs + qt_include_dirs +
                                           boost_include_dirs + gsl_include_dirs +
                                           openmm_include_dirs,
                                define_symbols = ["GCCXML_PARSE", "__PIE__",
                                                  "SIRE_USE_OPENMM",
                                                  "SIRE_SKIP_INLINE_FUNCTIONS",
                                                  "SIREN_SKIP_INLINE_FUNCTIONS",
                                                  "SIRE_INSTANTIATE_TEMPLATES",
                                                  "SIREN_INSTANTIATE_TEMPLATES"]
                         )

        mb = module_builder_t( files = [ "active_headers.h" ],
                               gccxml_config=xml_generator_config )


    #get rid of all virtual python functions - this is to stop slow wrapper code
    #from being generated for C++ virtual objects
    for calldef in mb.calldefs():
        try:
            calldef.virtuality = declarations.VIRTUALITY_TYPES.NOT_VIRTUAL
        except:
            pass

    #add calls to additional hand-written code
    if os.path.exists("%s_containers.h" % sourcedir):
        mb.add_declaration_code( "#include \"%s_containers.h\"" % sourcedir )
        mb.add_registration_code( "register_%s_containers();" % sourcedir, tail=False )

    mb.calldefs().create_with_signature = True

    metaheaders = []

    #export each of the classes in this module in turn
    for header in active_headers:
        classes = active_headers[header].classes()
        includes = active_headers[header].dependencies()
        aliases = active_headers[header].aliases()
        functions = active_headers[header].functions()

        for clas in classes:
            print("Trying to export the class %s" % clas)
            export_class(mb, clas, aliases, includes, special_code)

        for func in functions:
            export_function(mb, func, includes)

        if len(active_headers[header].metaTypes()) > 0:
           metaheaders.append(header)

    if len(metaheaders) > 0:
        mb.add_declaration_code( "#include \"%s_registrars.h\"" % sourcedir )
        mb.add_registration_code( "register_%s_objects();" % sourcedir, tail=False )

        FILE = open("%s_registrars.h" % sourcedir, "w")

        print(r"//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!", file=FILE)
        print("#ifndef PYWRAP_%s_REGISTRARS_H" % sourcedir, file=FILE)
        print("#define PYWRAP_%s_REGISTRARS_H" % sourcedir, file=FILE)
        print("void register_%s_objects();" % sourcedir, file=FILE)
        print("#endif\n", file=FILE)
        FILE.close()

        FILE = open("%s_registrars.cpp" % sourcedir, "w")

        print(r"//WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!", file=FILE)
        print("#include <Python.h>\n", file=FILE)
        print("#include \"%s_registrars.h\"\n" % sourcedir, file=FILE)
        
        for header in metaheaders:
            print("#include \"%s\"" % header, file=FILE)

        print("\n#include \"Helpers/objectregistry.hpp\"\n", file=FILE)

        print("void register_%s_objects()\n{\n" % sourcedir, file=FILE)

        for header in metaheaders:
            metatypes = active_headers[header].metaTypes()

            for metatype in metatypes:
                print("    ObjectRegistry::registerConverterFor< %s >();" % metatype, file=FILE)

        print("\n}\n", file=FILE)

        FILE.close()

    #now export all of the namespace-level operators
    for operator in mb.operators():
        p = str(operator.parent)
        if p.find(module) != -1 and p.find("[namespace]") != -1:
            operator.include()        

    #write the code that wraps up the Property classes
    writePropertyWrappers(mb, sourcedir, active_headers)

    #remove all implicit implicit conversions and add the explicit implicit conversions (!)
    register_implicit_conversions(mb, implicitly_convertible)

    #now perform any last-minute fixes
    fixMB(mb)

    write_wrappers(mb, module, huge_classes)  

    #now write a CMakeFile that contains all of the autogenerated files
    FILE = open("CMakeAutogenFile.txt", "w")

    print("# WARNING - AUTOGENERATED FILE - CONTENTS WILL BE OVERWRITTEN!", file=FILE)
    print("set ( PYPP_SOURCES", file=FILE)

    pyppfiles = glob("*.pypp.cpp")

    for pyppfile in pyppfiles:
        print("       %s" % pyppfile, file=FILE)

    if os.path.exists("%s_containers.cpp" % sourcedir):
        print("       %s_containers.cpp" % sourcedir, file=FILE)

    if os.path.exists("%s_properties.cpp" % sourcedir):
        print("       %s_properties.cpp" % sourcedir, file=FILE)

    if os.path.exists("%s_registrars.cpp" % sourcedir):
        print("       %s_registrars.cpp" % sourcedir, file=FILE)

    print("    )", file=FILE)

    FILE.close()


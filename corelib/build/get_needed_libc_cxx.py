
from Sire.Base import *

import os
import sys

# Get a list of all of the Sire libraries (installed and bundled)

def isLibrary(filename):
    return filename.endswith(".so") or filename.endswith(".dylib")

libs = []

for filename in os.listdir(getLibDir()):
    filename = "%s/%s" % (getLibDir(),filename)
    if isLibrary(filename):
        libs.append(filename)

for filename in os.listdir(getBundledLibDir()):
    filename = "%s/%s" % (getBundledLibDir(),filename)
    if isLibrary(filename):
        libs.append(filename)


try:
    objdump = findExe("objdump")
    #print(objdump.absoluteFilePath())
except:
    sys.stderr.write("Could not find 'objdump'\n")
    sys.stderr.flush()
    objdump = None

libc_words = {}
libcxx_words = {}
glibc_words = {}
glibcxx_words = {}

def getLibVersion(lib):

    sys.stderr.write("Checking library %s...\n" % lib)
    sys.stderr.flush()

    # Do we have objdump?
    if objdump:
        PIPE = os.popen("%s -T %s" % (objdump.absoluteFilePath(),lib), "r")
        lines = PIPE.readlines()
        for line in lines:
           for word in line.split():
               if word.find("GLIBCXX") != -1:
                   glibcxx_words[word] = 1
               elif word.find("GLIBC") != -1:
                   glibc_words[word] = 1
    
for lib in libs:
    getLibVersion(lib)

libc_versions = list(libc_words.keys())
libcxx_versions = list(libcxx_words.keys())
glibc_versions = list(glibc_words.keys())
glibcxx_versions = list(glibcxx_words.keys())

libc_versions.sort()
libcxx_versions.sort()
glibc_versions.sort()
glibcxx_versions.sort()

print( "LIBC    VERSION STRINGS: %s" % libc_versions )
print("LIBC++   VERSION STRINGS: %s" % libcxx_versions )
print("GLIBC    VERSION STRINGS: %s" % glibc_versions )
print("GLIBC++  VERSION STRINGS: %s" % glibcxx_versions )




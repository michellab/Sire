/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2007  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "getbacktrace.h"

#include <QObject>
#include <QString>
#include <QRegExp>
#include <QDebug>

#ifdef __INTEL_COMPILER
  // This causes a segfault when using intel's compiler. Probably
  // just an incompatibility or my bug
  #undef _HAVE_EXECINFO_H_
#endif

#ifdef _HAVE_EXECINFO_H_
  #include <execinfo.h>
  #include <cxxabi.h>
#endif

namespace SireError
{
        
//need to extract the symbol from the output of 'backtrace_symbols'

//a typical output from backtrace_symbols will look like;
//a.out(_ZNK1A12getBackTraceEv+0x12) [0x804ad36]

// This needs to be split into;
//  (1) The program or library containing the symbol (a.out)
//  (2) The symbol itself (_ZNK1A12getBackTraceEv)
//  (3) The offset? +0x12
//  (4) The symbol address ([0x804ad36])

// This is achieved by the following regexp
//              (unit )  (symbol)   (offset)          (address)
QRegExp regexp("([^(]+)\\(([^)^+]+)(\\+[^)]+)\\)\\s(\\[[^]]+\\])");

//However, on OS X the output looks something like this;
//2 libSireBase.0.dylib 0x00da01a5 _ZNK8SireBase10PropertiesixERKNS_12PropertyNameE + 595
//
// Word 2 is the library, word 3 is the symbol address, word 4 is the symbol
// itself and word 6 is the offset(?)

/** Obtain a backtrace and return as a QStringList.
    This is not well-optimised, requires compilation with "-rdynamic" on linux
    and doesn't do a great job of demangling the symbol names. It is sufficient
    though to work out call trace. */
QStringList SIREERROR_EXPORT getBackTrace()
{
    //now get the backtrace of the code at this point
    //(we can only do this if we have 'execinfo.h'
#ifdef _HAVE_EXECINFO_H_
    
    //create a void* array to hold the function addresses. We will only go at most 128 deep
    void *func_addresses[128];
    int nfuncs = backtrace(func_addresses, 128);

    //now get the function names associated with these symbols. This should work for elf
    //binaries, though additional linker options may need to have been called 
    //(e.g. -rdynamic for GNU ld. See the glibc documentation for 'backtrace')
    char **symbols = backtrace_symbols(func_addresses, nfuncs);
    
    //save all of the function names onto the QStringList....
    //(note that none of this will work if we have run out of memory)
    QStringList ret;

    if (nfuncs == 1)
    {
        //we have probably been compiled with -fomit-frame-pointer, so this
        //has only been able to get the backtrace back to the current function
        ret.append( QObject::tr("This is an incomplete backtrace as it looks "
                    "like this code was compiled without a frame pointer\n"
                    "(e.g. using -fomit-frame-pointer)") );
    }

    for (int i=0; i<nfuncs; i++)
    {
        if (regexp.indexIn(symbols[i]) != -1)
        {
            //get the library or app that contains this symbol
            QString unit = regexp.cap(1);
            //get the symbol
            QString symbol = regexp.cap(2);
            //get the offset
            QString offset = regexp.cap(3);
            //get the address
            QString address = regexp.cap(4);
        
            //now try and demangle the symbol
            int stat;
            char *demangled = 
                    abi::__cxa_demangle(qPrintable(symbol),0,0,&stat);
        
            if (demangled)
            {
                symbol = demangled;
                delete demangled;
            }


            //put this all together
            ret.append( QString("(%1) %2 (%3 +%4)\n  -- %5\n")
                                .arg(QString::number(i), 3)
                                .arg(unit).arg(address,offset)
                                .arg(symbol) );
        }
        else
        {
            //split line into words
            QStringList words = QString(symbols[i])
                                    .split(" ", QString::SkipEmptyParts);
            
            if (words.count() == 6 and words[4] == "+")
            {
                //this is probably an OS X line...

                //get the library or app that contains this symbol
                QString unit = words[1];
                //get the symbol
                QString symbol = words[3];
                //get the offset
                QString offset = words[5];
                //get the address
                QString address = words[2];
        
                //now try and demangle the symbol
                int stat;
                char *demangled = 
                        abi::__cxa_demangle(qPrintable(symbol),0,0,&stat);
        
                if (demangled)
                {
                    symbol = demangled;
                    delete demangled;
                }

                //put this all together
                ret.append( QString("(%1) %2 (%3 +%4)\n  -- %5\n")
                                    .arg(QString::number(i), 3)
                                    .arg(unit).arg(address,offset)
                                    .arg(symbol) );
            }
            else
                //I don't recognise this string - just add the raw
                //string to the backtrace
                ret.append(symbols[i]);
        }
    }
    
    //we now need to release the memory of the symbols array. Since it was allocated using
    //malloc, we must release it using 'free'
    free(symbols);

    return ret;

#else
    return QStringList( QObject::tr(
                "Backtrace is not available on this system. Backtrace is "
                "available on Linux, Mac OS X (>=10.5) and any other system "
                "that provides the backtrace_symbols() API found in "
                "execinfo.h")
                      );
#endif

}

}

#ifndef SIREGLOBAL_H
#define SIREGLOBAL_H

/** This file contains some global macros that are used
    to improve the portability of the code. */

//include qglobal.h - this includes all of the definitions
//that are needed for these macros
#include <qconfig.h>
#include <qglobal.h>

#ifdef __cplusplus

#include <qmetatype.h>
#include <QString>
#include <QList>
#include <QSet>
#include <QVector>

#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
    // Newer Qt moves QString::SkipEmptyParts to the Qt namespace,
    // and deprecates the functions that use the QString value.
    // However, Qt < 5.14 doesn't have SkipEmptyParts defined in
    // the Qt namespace, so we need to do this here
    namespace Qt
    {
        const auto SkipEmptyParts = QString::SkipEmptyParts;
    }

    // Older Qt also uses a now deprecated way to initialise QSets
    template<class T>
    inline QSet<T> convert_to_qset(const QList<T> &l)
    {
        return l.toSet();
    }

    template<class T>
    inline QSet<T> convert_to_qset(const QVector<T> &v)
    {
        return v.toList().toSet();
    }

    template<class T>
    inline QList<T> convert_to_qlist(const QVector<T> &v)
    {
        return v.toList();
    }
#else
    //Newer Qt also has a different way of initialising QSet!
    template<class T>
    inline QSet<T> convert_to_qset(const QList<T> &l)
    {
        return QSet<T>(l.begin(), l.end());
    }

    template<class T>
    inline QSet<T> convert_to_qset(const QVector<T> &v)
    {
        return QSet<T>(v.begin(), v.end());
    }

    template<class T>
    inline QList<T> convert_to_qlist(const QVector<T> &v)
    {
        return QList<T>(v.begin(), v.end());
    }
#endif

#include <boost/current_function.hpp>
//macro used to return the current file and line as a QString
#ifdef CODELOC
#undef CODELOC
#endif

#define CODELOC QString("FILE: %1, LINE: %2, FUNCTION: %3").arg(__FILE__).arg(__LINE__).arg(BOOST_CURRENT_FUNCTION)

//copy the QT_BEGIN_HEADER and QT_END_HEADER
//to SIRE_BEGIN_HEADER and SIRE_END_HEADER. This
//will allow me to change their definition at some
//future point in time without having to mess with
//Qt
#define SIRE_BEGIN_HEADER QT_BEGIN_HEADER
#define SIRE_END_HEADER QT_END_HEADER

//create keywords that control whether classes or functions are
//exposed to a scripting language
#define SIRE_EXPOSE_FUNCTION(f)  /* Exposing function #1 */
#define SIRE_EXPOSE_CLASS(c)     /* Exposing class #1 */
#define SIRE_EXPOSE_ALIAS(c,a)   /* Exposing class #1 as alias #2 */
#define SIRE_EXPOSE_PROPERTY(c,a)  /* Exposing property #1 of base class #2 */
#define SIRE_EXPOSE_ATOM_PROPERTY(c,a) /* Exposing atom property #1 with alias #2 */
#define SIRE_EXPOSE_CUTGROUP_PROPERTY(c,a) /* Exposing CutGroup property #1 with alias #2 */
#define SIRE_EXPOSE_RESIDUE_PROPERTY(c,a) /* Exposing residue property #1 with alias #2 */
#define SIRE_EXPOSE_CHAIN_PROPERTY(c,a) /* Exposing chain property #1 with alias #2 */
#define SIRE_EXPOSE_SEGMENT_PROPERTY(c,a) /* Exposing segment property #1 with alias #2 */
#define SIRE_EXPOSE_BEAD_PROPERTY(c,a)  /* Exposing bead property #1 with alias #2 */

#ifdef _WIN32
#define SIRE_EXPORT __declspec(dllexport)
#define SIRE_IMPORT __declspec(dllimport)
#else
//create the keyword used to export a symbol - this
//is a copy of Q_DECL_EXPORT, which will definitely
//be set to the correct value
#ifndef SIRE_NO_VISIBILITY_AVAILABLE
#define SIRE_EXPORT Q_DECL_EXPORT

//create the keyword to import a symbol - this is copied
//from Q_DECL_IMPORT, which will definitely be set to
//the correct value
#define SIRE_IMPORT Q_DECL_IMPORT
#else
#define SIRE_EXPORT
#define SIRE_IMPORT
#endif
#endif

//create the keyword to fix symbol visibility problems for out-of-line
//template functions
#define SIRE_OUTOFLINE_TEMPLATE Q_OUTOFLINE_TEMPLATE

//do the same for inline template functions
#define SIRE_INLINE_TEMPLATE Q_INLINE_TEMPLATE

//Sire is much more picky about what the compiler can support
//than Qt - use the results of the Qt portability tests
//to check that the C++ compiler is up to scratch

#ifdef Q_BROKEN_TEMPLATE_SPECIALIZATION
#error Sire requires a C++ compiler that can handle templates in \
       a standards compliant manner. Please use a more recent version \
       of your compiler, or try gcc >= 3.0
#endif

#ifdef QT_NO_PARTIAL_TEMPLATE_SPECIALIZATION
#error Sire requires that the C++ compiler supports \
       partial template specialization. Please use a more \
     recent version of your compiler, or try gcc >= 3.0
#endif

#ifdef Q_NO_USING_KEYWORD
#error Sire requires that the C++ compiler supports the 'using' \
       keyword. Please use a more recent version of your compiler \
     or try gcc >= 3.0
#endif

#ifdef QT_NO_MEMBER_TEMPLATES
#error Sire requires that the C++ compiler supports the use of \
       member templates. Please use a more recent version of your \
     compiler or try gcc >= 3.0
#endif

#ifdef QT_NO_THREAD
#error Sire requires that Qt was built with threading enabled. \
       Please recompile Qt with threading enabled.
#endif

#ifdef QT_NO_TEXTSTREAM
#error Sire requires that Qt was built with QTextStream enabled \
       Please recompile Qt with QTextStream enabled.
#endif

/** These functions convert a pointer into an integer.

    toInt(void*) : Returns an integer from the pointer. Returns
                   an integer of the same size as the pointer.

    toUInt(void*) : Returns an unsigned integer of the same size as
                    the pointer.

    toInt32(void*)   : Returns a 32bit integer from a pointer. This
                       will truncate if the pointer is larger than 32bits.

    toUInt32(void*)  : Same as toInt32, but returns the pointer as an
                       unsigned 32bit integer.

    toInt64(void*) : Returns a 64bit integer

    toUInt64(void*) : Returns a 64bit unsigned integer.

    @author Christopher Woods
*/
#ifdef QT_POINTER_SIZE
  #if QT_POINTER_SIZE == 4

    // Functions used if we have 32bit pointers
    SIRE_ALWAYS_INLINE qint32 toInt(const void *ptr)
    {
        return qint32(ptr);
    }

    SIRE_ALWAYS_INLINE quint32 toUInt(const void *ptr)
    {
        return quint32(ptr);
    }

    SIRE_ALWAYS_INLINE qint32 toInt32(const void *ptr)
    {
        return toInt(ptr);
    }

    SIRE_ALWAYS_INLINE quint32 toUInt32(const void *ptr)
    {
        return toUInt(ptr);
    }

    SIRE_ALWAYS_INLINE qint64 toInt64(const void *ptr)
    {
        return qint64( toInt(ptr) );
    }

    SIRE_ALWAYS_INLINE quint64 toUInt64(const void *ptr)
    {
        return quint64( toUInt(ptr) );
    }

  #elif QT_POINTER_SIZE == 8

    // Functions used if we have 64bit pointers
    SIRE_ALWAYS_INLINE qint64 toInt(const void *ptr)
    {
        return qint64(ptr);
    }

    SIRE_ALWAYS_INLINE quint64 toUInt(const void *ptr)
    {
        return quint64(ptr);
    }

    SIRE_ALWAYS_INLINE qint32 toInt32(const void *ptr)
    {
        return qint32( toInt(ptr) );
    }

    SIRE_ALWAYS_INLINE quint32 toUInt32(const void *ptr)
    {
        return qint32( toUInt(ptr) );
    }

    SIRE_ALWAYS_INLINE qint64 toInt64(const void *ptr)
    {
        return toInt(ptr);
    }

    SIRE_ALWAYS_INLINE quint64 toUInt64(const void *ptr)
    {
        return toUInt(ptr);
    }

  #else
    #fatal Invalid QT_POINTER_SIZE value stored!
  #endif
#else
  #fatal No QT_POINTER_SIZE macro defined!
#endif

// MSVC does not support the #warning preprocessor directive
#ifdef _MSC_VER
#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)
#endif

//I now define seperate SIRE_EXPORT macros for each of the different Sire libraries.
// SIREANALYSIS_EXPORT definitions
#ifdef _WIN32
#ifdef SIREANALYSIS_BUILD
#define SIREANALYSIS_EXPORT SIRE_EXPORT
#else
#define SIREANALYSIS_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREANALYSIS_EXPORT
#define SIREANALYSIS_EXPORT SIRE_EXPORT
#endif
// SIREANALYSIS_EXPORT end definitions
// SIREBASE_EXPORT definitions
#ifdef _WIN32
#ifdef SIREBASE_BUILD
#define SIREBASE_EXPORT SIRE_EXPORT
#else
#define SIREBASE_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREBASE_EXPORT
#define SIREBASE_EXPORT SIRE_EXPORT
#endif
// SIREBASE_EXPORT end definitions
// SIRECAS_EXPORT definitions
#ifdef _WIN32
#ifdef SIRECAS_BUILD
#define SIRECAS_EXPORT SIRE_EXPORT
#else
#define SIRECAS_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIRECAS_EXPORT
#define SIRECAS_EXPORT SIRE_EXPORT
#endif
// SIRECAS_EXPORT end definitions
// SIRECLUSTER_EXPORT definitions
#ifdef _WIN32
#ifdef SIRECLUSTER_BUILD
#define SIRECLUSTER_EXPORT SIRE_EXPORT
#else
#define SIRECLUSTER_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIRECLUSTER_EXPORT
#define SIRECLUSTER_EXPORT SIRE_EXPORT
#endif
// SIRECLUSTER_EXPORT end definitions
// SIREERROR_EXPORT definitions
#ifdef _WIN32
#ifdef SIREERROR_BUILD
#define SIREERROR_EXPORT SIRE_EXPORT
#else
#define SIREERROR_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREERROR_EXPORT
#define SIREERROR_EXPORT SIRE_EXPORT
#endif
// SIREERROR_EXPORT end definitions
// SIREFF_EXPORT definitions
#ifdef _WIN32
#ifdef SIREFF_BUILD
#define SIREFF_EXPORT SIRE_EXPORT
#else
#define SIREFF_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREFF_EXPORT
#define SIREFF_EXPORT SIRE_EXPORT
#endif
// SIREFF_EXPORT end definitions
// SIREID_EXPORT definitions
#ifdef _WIN32
#ifdef SIREID_BUILD
#define SIREID_EXPORT SIRE_EXPORT
#else
#define SIREID_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREID_EXPORT
#define SIREID_EXPORT SIRE_EXPORT
#endif
// SIREID_EXPORT end definitions
// SIREIO_EXPORT definitions
#ifdef _WIN32
#ifdef SIREIO_BUILD
#define SIREIO_EXPORT SIRE_EXPORT
#else
#define SIREIO_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREIO_EXPORT
#define SIREIO_EXPORT SIRE_EXPORT
#endif
// SIREIO_EXPORT end definitions
// SIREMATHS_EXPORT definitions
#ifdef _WIN32
#ifdef SIREMATHS_BUILD
#define SIREMATHS_EXPORT SIRE_EXPORT
#else
#define SIREMATHS_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREMATHS_EXPORT
#define SIREMATHS_EXPORT SIRE_EXPORT
#endif
// SIREMATHS_EXPORT end definitions
// SIREMM_EXPORT definitions
#ifdef _WIN32
#ifdef SIREMM_BUILD
#define SIREMM_EXPORT SIRE_EXPORT
#else
#define SIREMM_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREMM_EXPORT
#define SIREMM_EXPORT SIRE_EXPORT
#endif
// SIREMM_EXPORT end definitions
// SIREMOL_EXPORT definitions
#ifdef _WIN32
#ifdef SIREMOL_BUILD
#define SIREMOL_EXPORT SIRE_EXPORT
#else
#define SIREMOL_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREMOL_EXPORT
#define SIREMOL_EXPORT SIRE_EXPORT
#endif
// SIREMOL_EXPORT end definitions
// SIREMOVE_EXPORT definitions
#ifdef _WIN32
#ifdef SIREMOVE_BUILD
#define SIREMOVE_EXPORT SIRE_EXPORT
#else
#define SIREMOVE_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREMOVE_EXPORT
#define SIREMOVE_EXPORT SIRE_EXPORT
#endif
// SIREMOVE_EXPORT end definitions
// SIRESEARCH_EXPORT definitions
#ifdef _WIN32
#ifdef SIRESEARCH_BUILD
#define SIRESEARCH_EXPORT SIRE_EXPORT
#else
#define SIRESEARCH_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIRESEARCH_EXPORT
#define SIRESEARCH_EXPORT SIRE_EXPORT
#endif
// SIRESEARCH_EXPORT end definitions
// SIRESTREAM_EXPORT definitions
#ifdef _WIN32
#ifdef SIRESTREAM_BUILD
#define SIRESTREAM_EXPORT SIRE_EXPORT
#else
#define SIRESTREAM_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIRESTREAM_EXPORT
#define SIRESTREAM_EXPORT SIRE_EXPORT
#endif
// SIRESTREAM_EXPORT end definitions
// SIRESYSTEM_EXPORT definitions
#ifdef _WIN32
#ifdef SIRESYSTEM_BUILD
#define SIRESYSTEM_EXPORT SIRE_EXPORT
#else
#define SIRESYSTEM_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIRESYSTEM_EXPORT
#define SIRESYSTEM_EXPORT SIRE_EXPORT
#endif
// SIRESYSTEM_EXPORT end definitions
// SIREUNITS_EXPORT definitions
#ifdef _WIN32
#ifdef SIREUNITS_BUILD
#define SIREUNITS_EXPORT SIRE_EXPORT
#else
#define SIREUNITS_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREUNITS_EXPORT
#define SIREUNITS_EXPORT SIRE_EXPORT
#endif
// SIREUNITS_EXPORT end definitions
// SIREVOL_EXPORT definitions
#ifdef _WIN32
#ifdef SIREVOL_BUILD
#define SIREVOL_EXPORT SIRE_EXPORT
#else
#define SIREVOL_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SIREVOL_EXPORT
#define SIREVOL_EXPORT SIRE_EXPORT
#endif
// SIREVOL_EXPORT end definitions
// SQUIRE_EXPORT definitions
#ifdef _WIN32
#ifdef SQUIRE_BUILD
#define SQUIRE_EXPORT SIRE_EXPORT
#else
#define SQUIRE_EXPORT SIRE_IMPORT
#endif
#endif
#ifndef SQUIRE_EXPORT
#define SQUIRE_EXPORT SIRE_EXPORT
#endif
// SQUIRE_EXPORT end definitions

//Add functions that are used to register all public
//types (this allows them to be streamed to a binary file and
//created dynamically)
#include <QMetaType>

namespace Sire
{
typedef quint32 MagicID;

SIRESTREAM_EXPORT MagicID getMagic(const char *type_name);

/** Enum used to register with only the magic ID */
enum MagicOnlyEnum{ MAGIC_ONLY = 1 };
/** Enum used to register without using the streaming operators */
enum NoStreamEnum{ NO_STREAM = 2 };
/** Enum used to register a class without a root class */
enum NoRootEnum{ NO_ROOT = 3 };

class RegisterMetaTypeBase
{
public:
    RegisterMetaTypeBase(MagicID magic, const char *type_name, int ID)
          : magicid(magic), typnam(type_name), id(ID)
    {}

    ~RegisterMetaTypeBase()
    {}

    MagicID magicID() const
    {
        return magicid;
    }

    const char* typeName() const
    {
        return typnam;
    }

    QString typeNameString() const
    {
        return QLatin1String(typnam);
    }

    int ID() const
    {
        return id;
    }

private:
    MagicID magicid;
    const char *typnam;
    int id;
};

namespace detail
{
    SIREERROR_EXPORT const QHash< QString, QSet<QString> > branchClasses();
    SIREERROR_EXPORT const QHash< QString, QSet<QString> > leafClasses();
    SIREERROR_EXPORT const QSet<QString> rootlessClasses();
    SIREERROR_EXPORT void registerLeaf(const QString &type_name, const char *root_class);
    SIREERROR_EXPORT void registerBranch(const QString &type_name, const char *root_class);
    SIREERROR_EXPORT void registerRootless(const QString &type_name);
}

/** This is used to register the type 'T' - this
    needs to be called once for each public type.

    @author Christopher Woods
*/
template<class T>
class RegisterMetaType : public RegisterMetaTypeBase
{
public:

    /** Use this constructor to register a concrete class */
    RegisterMetaType()
        : RegisterMetaTypeBase( getMagic( QMetaType::typeName( qMetaTypeId<T>() ) ),
                                QMetaType::typeName( qMetaTypeId<T>() ),
                                qMetaTypeId<T>() )
    {
        qRegisterMetaType<T>(this->typeName());
        qRegisterMetaTypeStreamOperators<T>(this->typeName());
        detail::registerLeaf(this->typeNameString(), T::ROOT::typeName());
    }

    /** Use this constructor to register a concrete class with no root */
    RegisterMetaType(NoRootEnum )
        : RegisterMetaTypeBase( getMagic( QMetaType::typeName( qMetaTypeId<T>() ) ),
                                QMetaType::typeName( qMetaTypeId<T>() ),
                                qMetaTypeId<T>() )
    {
        qRegisterMetaType<T>(this->typeName());
        qRegisterMetaTypeStreamOperators<T>(this->typeName());
        detail::registerRootless(this->typeNameString());
    }

    /** Use this construct to register a pure-virtual class */
    RegisterMetaType(MagicOnlyEnum, const char *type_name)
        : RegisterMetaTypeBase( getMagic(type_name), type_name, 0 )
    {
        detail::registerBranch(this->typeNameString(), T::ROOT::typeName());
    }

    /** Use this construct to register a pure-virtual class with no registered root */
    RegisterMetaType(MagicOnlyEnum, NoRootEnum, const char *type_name)
        : RegisterMetaTypeBase( getMagic(type_name), type_name, 0 )
    {
        detail::registerRootless(this->typeNameString());
    }

    /** Use this constructor to register a class that cannot be streamed */
    RegisterMetaType(NoStreamEnum )
        : RegisterMetaTypeBase( getMagic( QMetaType::typeName( qMetaTypeId<T>() ) ),
                                QMetaType::typeName( qMetaTypeId<T>() ),
                                qMetaTypeId<T>() )
    {
        qRegisterMetaType<T>(this->typeName());
        detail::registerLeaf(this->typeNameString(), T::ROOT::typeName());
    }

    /** Use this constructor to register a class that cannot be streamed
        and that has no registered root */
    RegisterMetaType(NoStreamEnum, NoRootEnum )
        : RegisterMetaTypeBase( getMagic( QMetaType::typeName( qMetaTypeId<T>() ) ),
                                QMetaType::typeName( qMetaTypeId<T>() ),
                                qMetaTypeId<T>() )
    {
        qRegisterMetaType<T>(this->typeName());
        detail::registerRootless(this->typeNameString());
    }

    /** Destructor */
    ~RegisterMetaType()
    {}
};

} // namespace Sire

namespace SireStream
{
class XMLStream; // This class is used to stream objects to and from XML files

} // namespace SireStream

using Sire::RegisterMetaTypeBase;
using Sire::RegisterMetaType;
using Sire::MAGIC_ONLY;
using Sire::NO_STREAM;
using Sire::NO_ROOT;
using Sire::MagicID;
using Sire::getMagic;
using SireStream::XMLStream;

class QDataStream;
class QTextStream;

#else  // else #ifdef __cplusplus

//copied directly from qglobal.h /////////////////
#ifndef Q_DECL_EXPORT
#  ifdef Q_OS_WIN
#    define Q_DECL_EXPORT __declspec(dllexport)
#  elif defined(QT_VISIBILITY_AVAILABLE)
#    define Q_DECL_EXPORT __attribute__((visibility("default")))
#  endif
#  ifndef Q_DECL_EXPORT
#    define Q_DECL_EXPORT
#  endif
#endif
#ifndef Q_DECL_IMPORT
#  ifdef Q_OS_WIN
#    define Q_DECL_IMPORT __declspec(dllimport)
#  else
#    define Q_DECL_IMPORT
#  endif
#endif
///////////////////////////////////////////////////

#ifndef SIRE_NO_VISIBILITY_AVAILABLE
#define SIRE_EXPORT Q_DECL_EXPORT
#define SIRE_IMPORT Q_DECL_IMPORT
#else
#define SIRE_EXPORT
#define SIRE_IMPORT
#endif

#endif // #ifdef __cplusplus

#endif // SIREGLOBAL_H

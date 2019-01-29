/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2006  Christopher Woods
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

#ifndef SIRECAS_TRIGFUNCS_H
#define SIRECAS_TRIGFUNCS_H

#include "singlefunc.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class Cos;
class Sin;
class Tan;
class Csc;
class Sec;
class Cot;
}

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Cos&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Cos&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Sin&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Sin&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Tan&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Tan&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Csc&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Csc&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Sec&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Sec&);

SIRECAS_EXPORT QDataStream& operator<<(QDataStream&, const SireCAS::Cot&);
SIRECAS_EXPORT QDataStream& operator>>(QDataStream&, SireCAS::Cot&);

namespace SireCAS
{

/** Cosine */
class SIRECAS_EXPORT Cos : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Cos&);
friend QDataStream& ::operator>>(QDataStream&, Cos&);

public:
    Cos();
    Cos(const Expression &ex);

    Cos(const Cos &other);

    ~Cos();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Cos::typeName();
    }

    Cos* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Cos(arg));
    }

    QString stringRep() const
    {
        return "cos";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;

};

/** Sine */
class SIRECAS_EXPORT Sin : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Sin&);
friend QDataStream& ::operator>>(QDataStream&, Sin&);

public:
    Sin();
    Sin(const Expression &ex);

    Sin(const Sin &other);

    ~Sin();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Sin::typeName();
    }

    Sin* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Sin(arg));
    }

    QString stringRep() const
    {
        return "sin";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Tangent */
class SIRECAS_EXPORT Tan : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Tan&);
friend QDataStream& ::operator>>(QDataStream&, Tan&);

public:
    Tan();
    Tan(const Expression &ex);

    Tan(const Tan &other);

    ~Tan();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Tan::typeName();
    }

    Tan* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Tan(arg));
    }

    QString stringRep() const
    {
        return "tan";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Secant */
class SIRECAS_EXPORT Sec : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Sec&);
friend QDataStream& ::operator>>(QDataStream&, Sec&);

public:
    Sec();
    Sec(const Expression &ex);

    Sec(const Sec &other);

    ~Sec();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Sec::typeName();
    }

    Sec* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Sec(arg));
    }

    QString stringRep() const
    {
        return "sec";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Cosecant */
class SIRECAS_EXPORT Csc : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Csc&);
friend QDataStream& ::operator>>(QDataStream&, Csc&);

public:
    Csc();
    Csc(const Expression &ex);

    Csc(const Csc &other);

    ~Csc();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return Csc::typeName();
    }

    Csc* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Csc(arg));
    }

    QString stringRep() const
    {
        return "csc";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Cotangent */
class SIRECAS_EXPORT Cot : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const Cot&);
friend QDataStream& ::operator>>(QDataStream&, Cot&);

public:
    Cot();
    Cot(const Expression &ex);

    Cot(const Cot &other);

    ~Cot();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char *what() const
    {
        return Cot::typeName();
    }

    Cot* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(Cot(arg));
    }

    QString stringRep() const
    {
        return "cot";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

}

Q_DECLARE_METATYPE(SireCAS::Cos)
Q_DECLARE_METATYPE(SireCAS::Sin)
Q_DECLARE_METATYPE(SireCAS::Tan)
Q_DECLARE_METATYPE(SireCAS::Csc)
Q_DECLARE_METATYPE(SireCAS::Sec)
Q_DECLARE_METATYPE(SireCAS::Cot)

SIRE_EXPOSE_CLASS( SireCAS::Cos )
SIRE_EXPOSE_CLASS( SireCAS::Sin )
SIRE_EXPOSE_CLASS( SireCAS::Tan )
SIRE_EXPOSE_CLASS( SireCAS::Csc )
SIRE_EXPOSE_CLASS( SireCAS::Sec )
SIRE_EXPOSE_CLASS( SireCAS::Cot )

SIRE_END_HEADER

#endif

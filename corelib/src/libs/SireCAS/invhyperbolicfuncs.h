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

#ifndef SIRECAS_INVHYPERBOLICFUNCS_H
#define SIRECAS_INVHYPERBOLICFUNCS_H

#include "singlefunc.h"

SIRE_BEGIN_HEADER

namespace SireCAS
{
class ArcCosh;
class ArcSinh;
class ArcTanh;
class ArcCsch;
class ArcSech;
class ArcCoth;
}

QDataStream& operator<<(QDataStream&, const SireCAS::ArcCosh&);
QDataStream& operator>>(QDataStream&, SireCAS::ArcCosh&);

QDataStream& operator<<(QDataStream&, const SireCAS::ArcSinh&);
QDataStream& operator>>(QDataStream&, SireCAS::ArcSinh&);

QDataStream& operator<<(QDataStream&, const SireCAS::ArcTanh&);
QDataStream& operator>>(QDataStream&, SireCAS::ArcTanh&);

QDataStream& operator<<(QDataStream&, const SireCAS::ArcCsch&);
QDataStream& operator>>(QDataStream&, SireCAS::ArcCsch&);

QDataStream& operator<<(QDataStream&, const SireCAS::ArcSech&);
QDataStream& operator>>(QDataStream&, SireCAS::ArcSech&);

QDataStream& operator<<(QDataStream&, const SireCAS::ArcCoth&);
QDataStream& operator>>(QDataStream&, SireCAS::ArcCoth&);

namespace SireCAS
{

/** Inverse-hyperbolic-cosine */
class SIRECAS_EXPORT ArcCosh : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const ArcCosh&);
friend QDataStream& ::operator>>(QDataStream&, ArcCosh&);

public:
    ArcCosh();
    ArcCosh(const Expression &ex);

    ArcCosh(const ArcCosh &other);

    ~ArcCosh();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ArcCosh::typeName();
    }

    ArcCosh* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(ArcCosh(arg));
    }

    QString stringRep() const
    {
        return "arccosh";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;

};

/** Inverse-hyperbolic-sine */
class SIRECAS_EXPORT ArcSinh : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const ArcSinh&);
friend QDataStream& ::operator>>(QDataStream&, ArcSinh&);

public:
    ArcSinh();
    ArcSinh(const Expression &ex);

    ArcSinh(const ArcSinh &other);

    ~ArcSinh();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ArcSinh::typeName();
    }

    ArcSinh* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(ArcSinh(arg));
    }

    QString stringRep() const
    {
        return "arcsinh";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Inverse-hyperbolic-tangent */
class SIRECAS_EXPORT ArcTanh : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const ArcTanh&);
friend QDataStream& ::operator>>(QDataStream&, ArcTanh&);

public:
    ArcTanh();
    ArcTanh(const Expression &ex);

    ArcTanh(const ArcTanh &other);

    ~ArcTanh();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ArcTanh::typeName();
    }

    ArcTanh* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(ArcTanh(arg));
    }

    QString stringRep() const
    {
        return "arctanh";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Inverse-hyperbolic-secant */
class SIRECAS_EXPORT ArcSech : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const ArcSech&);
friend QDataStream& ::operator>>(QDataStream&, ArcSech&);

public:
    ArcSech();
    ArcSech(const Expression &ex);

    ArcSech(const ArcSech &other);

    ~ArcSech();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ArcSech::typeName();
    }

    ArcSech* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(ArcSech(arg));
    }

    QString stringRep() const
    {
        return "arcsech";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Inverse-hyperbolic-cosecant */
class SIRECAS_EXPORT ArcCsch : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const ArcCsch&);
friend QDataStream& ::operator>>(QDataStream&, ArcCsch&);

public:
    ArcCsch();
    ArcCsch(const Expression &ex);

    ArcCsch(const ArcCsch &other);

    ~ArcCsch();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ArcCsch::typeName();
    }

    ArcCsch* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(ArcCsch(arg));
    }

    QString stringRep() const
    {
        return "arccsch";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

/** Inverse-hyperbolic-cotangent */
class SIRECAS_EXPORT ArcCoth : public SingleFunc
{

friend QDataStream& ::operator<<(QDataStream&, const ArcCoth&);
friend QDataStream& ::operator>>(QDataStream&, ArcCoth&);

public:
    ArcCoth();
    ArcCoth(const Expression &ex);

    ArcCoth(const ArcCoth &other);

    ~ArcCoth();

    bool operator==(const ExBase &other) const;

    static const char* typeName();

    const char* what() const
    {
        return ArcCoth::typeName();
    }

    ArcCoth* clone() const;

    double evaluate(const Values &values) const;
    Complex evaluate(const ComplexValues &values) const;

protected:
    Expression functionOf(const Expression &arg) const
    {
        if (arg == argument())
            return Expression(*this);
        else
            return Expression(ArcCoth(arg));
    }

    QString stringRep() const
    {
        return "arccoth";
    }

    uint magic() const;

    Expression diff() const;
    Expression integ() const;
};

}

Q_DECLARE_METATYPE(SireCAS::ArcCosh)
Q_DECLARE_METATYPE(SireCAS::ArcSinh)
Q_DECLARE_METATYPE(SireCAS::ArcTanh)
Q_DECLARE_METATYPE(SireCAS::ArcCsch)
Q_DECLARE_METATYPE(SireCAS::ArcSech)
Q_DECLARE_METATYPE(SireCAS::ArcCoth)

SIRE_EXPOSE_CLASS( SireCAS::ArcCosh )
SIRE_EXPOSE_CLASS( SireCAS::ArcSinh )
SIRE_EXPOSE_CLASS( SireCAS::ArcTanh )
SIRE_EXPOSE_CLASS( SireCAS::ArcCsch )
SIRE_EXPOSE_CLASS( SireCAS::ArcSech )
SIRE_EXPOSE_CLASS( SireCAS::ArcCoth )

SIRE_END_HEADER

#endif

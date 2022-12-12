/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2018  Christopher Woods
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

#include "ast.h"
#include "idengine.h"

#include "SireMol/core.h"

#include "SireError/errors.h"

using namespace parser_idengine;

namespace AST
{
    QString expression_to_string(const ExpressionVariant &expression)
    {
        return boost::apply_visitor( qstring_visitor(), expression );
    }

    QString idobject_to_string(IDObject obj)
    {
        switch(obj)
        {
        case ID_UNKNOWN:
            return "unknown";
        case ATOM:
            return "atom";
        case CUTGROUP:
            return "cutgroup";
        case RESIDUE:
            return "residue";
        case CHAIN:
            return "chain";
        case SEGMENT:
            return "segment";
        case MOLECULE:
            return "molecule";
        case BOND:
            return "bond";
        default:
            return "unknown";
        }
    }

    QString idnumtype_to_string(IDNumType typ)
    {
        switch(typ)
        {
        case ID_NUMBER:
            return "number";
        case ID_INDEX:
            return "index";
        default:
            return "unknown";
        }
    }

    QString idoperation_to_string(IDOperation op)
    {
        switch(op)
        {
        case ID_AND:
            return "and";
        case ID_OR:
            return "or";
        default:
            return "unknown";
        }
    }

    QString idcomparison_to_string(IDComparison cmp)
    {
        switch(cmp)
        {
        case ID_CMP_LT:
            return "<";
        case ID_CMP_LE:
            return "<=";
        case ID_CMP_AE:
            return "=~";
        case ID_CMP_EQ:
            return "==";
        case ID_CMP_NE:
            return "!=";
        case ID_CMP_GT:
            return ">";
        case ID_CMP_GE:
            return ">=";
        default:
            return "unknown";
        }
    }

    QString idbondtoken_to_string(IDBondToken token)
    {
        switch(token)
        {
        case ID_BOND_FROM:
            return "from";
        case ID_BOND_TO:
            return "to";
        default:
            return "unknown";
        }
    }

    QString idtoken_to_string(IDToken token)
    {
        switch(token)
        {
        case ID_WHERE:
            return "where";
        case ID_WITH:
            return "with";
        case ID_IN:
            return "in";
        default:
            return "unknown";
        }
    }

    QString idcoordtype_to_string(IDCoordType typ)
    {
        switch(typ)
        {
        case ID_COORD_CENTER:
            return "coords.center";
        case ID_COORD_CENTER_X:
            return "coords.center.x";
        case ID_COORD_CENTER_Y:
            return "coords.center.y";
        case ID_COORD_CENTER_Z:
            return "coords.center.z";
        case ID_COORD_MAX:
            return "coords.max";
        case ID_COORD_MAX_X:
            return "coords.max.x";
        case ID_COORD_MAX_Y:
            return "coords.max.y";
        case ID_COORD_MAX_Z:
            return "coords.max.z";
        case ID_COORD_MIN:
            return "coords.min";
        case ID_COORD_MIN_X:
            return "coords.min.x";
        case ID_COORD_MIN_Y:
            return "coords.min.y";
        case ID_COORD_MIN_Z:
            return "coords.min.z";
        case ID_COORD_X:
            return "coords.x";
        case ID_COORD_Y:
            return "coords.y";
        case ID_COORD_Z:
            return "coords.z";
        default:
            return "unknown";
        }
    }

    QString LengthValue::toString() const
    {
        return (value * unit).toString();
    }

    QString VectorValue::toString() const
    {
        if (_c == 0)
            return "unset";
        else if (_c == 1)
            return x.toString();
        else if (_c == 2)
            return QString("( %1, %2 )").arg(x.toString()).arg(y.toString());
        else
            return QString("( %1, %2, %3 )")
                        .arg(x.toString()).arg(y.toString()).arg(z.toString());
    }

    QString RegExpValue::toString() const
    {
        QString qstr = QString::fromStdString(value);

        if (is_case_sensitive)
            return QObject::tr("/%1/").arg(qstr);
        else
            return QObject::tr("/%1/i").arg(qstr);
    }

    QString CompareValue::toString() const
    {
        return QString("%1 %2").arg(idcomparison_to_string(compare))
                               .arg(value);
    }

    SireBase::Slice RangeValue::toSlice() const
    {
        int _start = 0;

        if (start.get() != 0)
            _start = *start;

        if (stop.get() == 0)
        {
            if (_c == 0)
                // single value
                return SireBase::Slice::fromStartStop(_start, _start);
            else if (step.get() == 0)
                return SireBase::Slice::fromStart(_start);
            else
                return SireBase::Slice::fromStart(_start, *step);
        }
        else
        {
            if (step.get() == 0)
                return SireBase::Slice::fromStartStop(_start, *stop);
            else
                return SireBase::Slice::fromStartStop(_start, *stop, *step);
        }
    }

    QString RangeValue::toString() const
    {
        if (_c == 0)
        {
            if (start.get() != 0)
                return QString::number(*start);
            else
                return "0";
        }

        QStringList parts;

        if (start.get() != 0)
            parts.append(QString::number(*start));
        else
            parts.append("");

        if (stop.get() != 0)
            parts.append(QString::number(*stop));
        else
            parts.append("");

        if (step.get() != 0)
            parts.append(QString::number(*step));
        else
            parts.append("");

        return parts.join(":");
    }

    QString Expression::toString() const
    {
        return boost::apply_visitor( qstring_visitor(), value );
    }

    QString ExpressionPart::toString() const
    {
        return boost::apply_visitor( qstring_visitor(), value );
    }

    QString IDUser::toString() const
    {
        return QString("{ %1 => %2 }").arg( QString::fromStdString(token) )
                                      .arg( boost::apply_visitor( qstring_visitor(), value ) );
    }

    QString Node::toString() const
    {
        return values.toString();
    }

    QString IDNull::toString() const
    {
        return QObject::tr("null");
    }

    QString IDName::toString() const
    {
        QStringList lines;
        for (const auto &value : values)
        {
            lines.append( value.toString() );
        }

        return QObject::tr("%1name %2")
                        .arg( idobject_to_string(name) )
                        .arg( lines.join(",") );
    }

    QString IDAll::toString() const
    {
        return QString("all %1").arg( idobject_to_string(name) );
    }

    QString IDWater::toString() const
    {
        return QString("water");
    }

    QString IDPerturbable::toString() const
    {
        return QString("perturbable");
    }

    QString IDNumber::toString() const
    {
        QStringList lines;
        for (const auto &value : values)
        {
            lines.append( boost::apply_visitor( qstring_visitor(), value ) );
        }

        return QObject::tr("%1%2 %3")
                    .arg( idobject_to_string(name) )
                    .arg( idnumtype_to_string(numtype) )
                    .arg( lines.join(",") );
    }

    QString IDElement::toString() const
    {
        QStringList lines;

        for (const auto &element : values)
        {
            if (element == "biological")
                lines.append(element);
            else
                lines.append( SireMol::Element(element).symbol() );
        }

        return QObject::tr("element %1").arg(lines.join(","));
    }

    QString IDBinary::toString() const
    {
        return QObject::tr("(%1 %2 %3)")
                    .arg(part0.toString())
                    .arg(idoperation_to_string(operation))
                    .arg(part1.toString());
    }

    QString IDProperty::toString() const
    {
        auto p = QString::fromStdString(this->property);
        auto v = QString("True");

        if (not this->value.empty())
        {
            v = QString::fromStdString(this->value);
        }

        return QObject::tr("%1 property %2 %3 %4")
                    .arg(idobject_to_string(name))
                    .arg(p)
                    .arg(idcomparison_to_string(compare))
                    .arg(v);
    }

    QString IDBond::toString() const
    {
        if (to_token != ID_BOND_UNKNOWN)
            return QObject::tr("bonds %1 %2 %3 %4")
                .arg(idbondtoken_to_string(from_token))
                .arg(from_value.toString())
                .arg(idbondtoken_to_string(to_token))
                .arg(to_value.toString());
        else
            return QObject::tr("bonds %1 %2")
                .arg(idbondtoken_to_string(from_token))
                .arg(from_value.toString());
    }

    QString IDWith::toString() const
    {
        return QObject::tr("(%1) %2 (%3)")
                    .arg(value0.toString())
                    .arg(idtoken_to_string(token))
                    .arg(value1.toString());
    }

    QString IDWhereWithin::toString() const
    {
        return QString("is within %1 of %2")
                    .arg(distance.toString())
                    .arg(value.toString());
    }

    QString IDWhereCompare::toString() const
    {
        return QString("%1 %2")
                .arg(idcomparison_to_string(compare))
                .arg(value.toString());
    }

    QString IDWhere::toString() const
    {
        return QString("%1s where %2 %3")
                    .arg(idobject_to_string(name))
                    .arg(idcoordtype_to_string(typ))
                    .arg(boost::apply_visitor( qstring_visitor(), value ));
    }

    QString IDJoin::toString() const
    {
        return QString("join (%1)").arg(value.toString());
    }

    QString IDNot::toString() const
    {
        return QString("not (%1)").arg(value.toString());
    }

    QString IDSubscript::toString() const
    {
        return QString("(%1)[%2]").arg(value.toString()).arg(range.toString());
    }

    QString IDWithin::toString() const
    {
        return QObject::tr("%1s within %2 of %3")
                    .arg(idobject_to_string(name))
                    .arg(distance.toString())
                    .arg(value.toString());
    }

    QString IDWithinVector::toString() const
    {
        return QObject::tr("%1s within %2 of %3")
                    .arg(idobject_to_string(name))
                    .arg(distance.toString())
                    .arg(value.toString());
    }

    SelectEnginePtr Node::toEngine() const
    {
        return values.toEngine();
    }

    SelectEnginePtr Expression::toEngine() const
    {
        return boost::apply_visitor( engine_visitor(), value );
    }

    SelectEnginePtr ExpressionPart::toEngine() const
    {
        return boost::apply_visitor( engine_visitor(), value );
    }

    SelectEnginePtr IDUser::toEngine() const
    {
        return boost::apply_visitor( engine_visitor(), value );
    }

    SelectEnginePtr IDName::toEngine() const
    {
        return IDNameEngine::construct(name,values);
    }

    SelectEnginePtr IDNull::toEngine() const
    {
        return SelectEnginePtr();
    }

    SelectEnginePtr IDElement::toEngine() const
    {
        return IDElementEngine::construct(values);
    }

    SelectEnginePtr IDNumber::toEngine() const
    {
        switch(numtype)
        {
        case ID_NUMBER:
            return IDNumberEngine::construct(name,values);
        case ID_INDEX:
            return IDIndexEngine::construct(name,values);
        default:
            return SelectEnginePtr();
        }
    }

    SelectEnginePtr IDAll::toEngine() const
    {
        return IDAllEngine::construct(name);
    }

    SelectEnginePtr IDWater::toEngine() const
    {
        return IDWaterEngine::construct();
    }

    SelectEnginePtr IDPerturbable::toEngine() const
    {
        return IDPerturbableEngine::construct();
    }

    SelectEnginePtr IDBinary::toEngine() const
    {
        switch(operation)
        {
        case ID_AND:
            return IDAndEngine::construct(part0.toEngine(),part1.toEngine());
        case ID_OR:
            return IDOrEngine::construct(part0.toEngine(),part1.toEngine());
        default:
            return SelectEnginePtr();
        }
    }

    SelectEnginePtr IDWith::toEngine() const
    {
        return IDWithEngine::construct(value0.toEngine(),
                                       token, value1.toEngine());
    }

    SelectEnginePtr IDProperty::toEngine() const
    {
        QString v;

        if (not this->value.empty())
        {
            v = QString::fromStdString(this->value);
        }
        else
        {
            v = "True";
        }

        return IDPropertyEngine::construct(this->name,
                                           QString::fromStdString(this->property),
                                           this->compare,
                                           v);
    }

    SelectEnginePtr IDBond::toEngine() const
    {
        return IDBondEngine::construct(from_token, from_value.toEngine(),
                                       to_token, to_value.toEngine());
    }

    SelectEnginePtr IDWhereWithin::toEngine(IDObject name, IDCoordType typ) const
    {
        return IDDistanceEngine::construct(name, typ, distance.value * distance.unit,
                                           value.toEngine());
    }

    SelectEnginePtr IDWhereCompare::toEngine(IDObject name, IDCoordType typ) const
    {
        qDebug() << "NOT YET IMPLEMENTED IDWhereCompare!";
        return SelectEnginePtr();
    }

    class where_engine_visitor : public boost::static_visitor<SelectEnginePtr>
    {
    public:
        where_engine_visitor(IDObject o, IDCoordType t)
            : boost::static_visitor<SelectEnginePtr>(), obj(o), typ(t)
        {}

        IDObject obj;
        IDCoordType typ;

        template<class T>
        SelectEnginePtr operator()(const T &value) const
        {
            return value.toEngine(obj,typ);
        }
    };

    SelectEnginePtr IDWhere::toEngine() const
    {
        return boost::apply_visitor( where_engine_visitor(name,typ), value );
    }

    SelectEnginePtr IDJoin::toEngine() const
    {
        return IDJoinEngine::construct(value.toEngine());
    }

    SelectEnginePtr IDNot::toEngine() const
    {
        return IDNotEngine::construct(value.toEngine());
    }

    SelectEnginePtr IDSubscript::toEngine() const
    {
        return IDSubScriptEngine::construct(value.toEngine(),range);
    }

    SelectEnginePtr IDWithin::toEngine() const
    {
        return IDDistanceEngine::construct(name, distance.value * distance.unit, value.toEngine());
    }

    SelectEnginePtr IDWithinVector::toEngine() const
    {
        return IDDistanceVectorEngine::construct(name, distance.value * distance.unit, value);
    }
}

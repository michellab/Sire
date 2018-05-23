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

#include "SireError/errors.h"

using namespace parser_idengine;

namespace AST
{
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

    QString RangeValue::toString() const
    {
        if (start == end)
            return QString::number(start);
        else if (step == 1)
            return QString("%1:%2").arg(start).arg(end);
        else
            return QString("%1:%2:%3").arg(start).arg(end).arg(step);
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
        QStringList lines;
        
        for (const auto value : values)
        {
            lines.append( value.toString() );
        }
        
        return lines.join("; ");
    }

    QString IDName::toString() const
    {
        QStringList lines;
        for (const auto value : values)
        {
            lines.append( value.toString() );
        }
    
        return QObject::tr("%1name %2")
                        .arg( idobject_to_string(name) )
                        .arg( lines.join(",") );
    }

    QString IDNumber::toString() const
    {
        QStringList lines;
        for (const auto value : values)
        {
            lines.append( boost::apply_visitor( qstring_visitor(), value ) );
        }
        
        return QObject::tr("%1%2 %3")
                    .arg( idobject_to_string(name) )
                    .arg( idnumtype_to_string(numtype) )
                    .arg( lines.join(",") );
    }

    QString IDBinary::toString() const
    {
        return QObject::tr("(%1 %2 %3)")
                    .arg(part0.toString())
                    .arg(idoperation_to_string(operation))
                    .arg(part1.toString());
    }

    QString IDWith::toString() const
    {
        return QObject::tr("%1s %2 %3")
                    .arg(idobject_to_string(name))
                    .arg(idtoken_to_string(token))
                    .arg(value.toString());
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
        return QString("{%1}[%2]").arg(value.toString()).arg(range.toString());
    }

    QString IDWithin::toString() const
    {
        return QObject::tr("%1s within %2 of %3")
                    .arg(idobject_to_string(name))
                    .arg(distance.toString())
                    .arg(value.toString());
    }

    SelectEnginePtr Node::toEngine() const
    {
        if (values.size() == 0)
            return SelectEnginePtr();
        else if (values.size() == 1)
            return values.at(0).toEngine();
        else
        {
            QList<SelectEnginePtr> engines;
        
            for (const auto value : values)
            {
                engines.append( value.toEngine() );
            }

            return SelectEnginePtr( new IDOrEngine(engines) );
        }
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
        return SelectEnginePtr( new IDNameEngine(name,values) );
    }
    
    SelectEnginePtr IDNumber::toEngine() const
    {
        switch(numtype)
        {
        case ID_NUMBER:
            return SelectEnginePtr( new IDNumberEngine(name,values) );
        case ID_INDEX:
            return SelectEnginePtr( new IDIndexEngine(name,values) );
        default:
            return SelectEnginePtr();
        }
    }
    
    SelectEnginePtr IDBinary::toEngine() const
    {
        switch(operation)
        {
        case ID_AND:
            return SelectEnginePtr( new IDAndEngine(part0.toEngine(),part1.toEngine()) );
        case ID_OR:
            return SelectEnginePtr( new IDOrEngine(part0.toEngine(),part1.toEngine()) );
        default:
            return SelectEnginePtr();
        }
    }
    
    SelectEnginePtr IDWith::toEngine() const
    {
        return SelectEnginePtr( new IDWithEngine(name, token, value.toEngine()) );
    }
    
    SelectEnginePtr IDWhereWithin::toEngine() const
    {
        return SelectEnginePtr();
    }
    
    SelectEnginePtr IDWhereCompare::toEngine() const
    {
        return SelectEnginePtr();
    }
    
    SelectEnginePtr IDWhere::toEngine() const
    {
        return SelectEnginePtr();
    }

    SelectEnginePtr IDJoin::toEngine() const
    {
        return SelectEnginePtr( new IDJoinEngine(value.toEngine()) );
    }

    SelectEnginePtr IDNot::toEngine() const
    {
        return SelectEnginePtr( new IDNotEngine(value.toEngine()) );
    }
    
    SelectEnginePtr IDSubscript::toEngine() const
    {
        return SelectEnginePtr( new IDSubScriptEngine(value.toEngine(),range) );
    }
    
    SelectEnginePtr IDWithin::toEngine() const
    {
        return SelectEnginePtr();
    }
}

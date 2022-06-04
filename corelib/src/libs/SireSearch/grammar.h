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

#ifndef SIRESEARCH_GRAMMAR_H

using qi::lit;
using qi::lexeme;
using qi::eps;
using qi::_1;
using qi::int_;
using qi::double_;
using qi::on_error;
using qi::fail;
using namespace qi::labels;
using qi::as_string;

using phoenix::construct;
using phoenix::val;

using boost::spirit::ascii::char_;

UserTokens getUserTokens();

/** This is the grammar that enables skipping of spaces, newlines and comments */
template<typename IteratorT>
class SkipperGrammar : public qi::grammar<IteratorT>
{
public:
    SkipperGrammar() : SkipperGrammar::base_type( rule )
    {
        lineCommentRule  = qi::lit( "//" ) >>
                           *(qi::char_ -qi::eol) >>
                           qi::eol;
        blockCommentRule = qi::lit( "/*" ) >>
                           *(qi::char_ -qi::lit( "*/" ) ) >>
                           qi::lit( "*/" );
        spaceRule        = qi::space;
        rule             = spaceRule | lineCommentRule | blockCommentRule;
    }

    qi::rule<IteratorT> lineCommentRule;
    qi::rule<IteratorT> blockCommentRule;
    qi::rule<IteratorT> spaceRule;
    qi::rule<IteratorT> rule;
};

/** This is a quoted string grammar that will parse quoted strings and also
    auto-escape characters */
template<typename IteratorT, typename SkipperT>
class ValueGrammar : public qi::grammar<IteratorT, std::string(), SkipperT>
{
public:
    ValueGrammar() : ValueGrammar::base_type( rule, "String" )
    {
        escapedStringRule %= qi::lexeme[
             qi::lit( "'" ) >>
             *( escapeCharSymbols | ( qi::char_ - qi::char_( "'" ) ) ) >>
             qi::lit( "'" ) ];

        rawStringRule %= qi::lexeme[
                    +( qi::alnum |
                       qi::char_( '.' ) |
                       qi::char_( '/' ) |
                       qi::char_( '_' ) |
                       qi::char_( '-' )
                      ) ];

        rule %= rawStringRule | escapedStringRule;

        escapeCharSymbols.add( "\\a", '\a' )
                             ( "\\b", '\b' )
                             ( "\\f", '\f' )
                             ( "\\n", '\n' )
                             ( "\\r", '\r' )
                             ( "\\t", '\t' )
                             ( "\\v", '\v' )
                             ( "\\\\", '\\' )
                             ( "\\\'", '\'' )
                             ( "\\\"", '\"' );

        escapedStringRule.name( "Escaped String" );
        rawStringRule.name( "Escaped String" );

        escapeCharSymbols.name( "Escaped Chars" );
    }

    qi::rule<IteratorT, std::string(), SkipperT>   escapedStringRule;
    qi::rule<IteratorT, std::string(), SkipperT>   rawStringRule;
    qi::rule<IteratorT, std::string(), SkipperT>   rule;
    qi::symbols<const char, const char>            escapeCharSymbols;
};

/** This the main grammar for the selection statements */
template<typename IteratorT, typename SkipperT>
class Grammar : public qi::grammar<IteratorT, AST::Node(), SkipperT>
{
public:
    Grammar() : Grammar::base_type( nodeRule, "Node" )
    {
        /////
        ///// first define all of the tokens recognised by the grammar
        /////

        //all of the different words to match "all"
        all_token.add( "*", AST::IDAll() )
                     ( "atoms", AST::IDAll(AST::ATOM) )
                     ( "bonds", AST::IDAll(AST::BOND) )
                     ( "residues", AST::IDAll(AST::RESIDUE) )
                     ( "chains", AST::IDAll(AST::CHAIN) )
                     ( "segments", AST::IDAll(AST::SEGMENT) )
                     ( "molecules", AST::IDAll(AST::MOLECULE) )
                     ( "cutgroups", AST::IDAll(AST::CUTGROUP) )
                     ;

        // all of the different tokens to match "water"
        water_token.add( "water", AST::IDWater() )
                       ( "WATER", AST::IDWater() )
                       ( "wat", AST::IDWater() )
                       ( "WAT", AST::IDWater() )
                       ( "waters", AST::IDWater() )
                       ( "WATERS", AST::IDWater() );

        // all of the different tokens to match "protein"
        protein_token.add( "protein", AST::IDProtein() )
                         ( "proteins", AST::IDProtein() )
                         ( "PROTEIN", AST::IDProtein() )
                         ( "PROTEINS", AST::IDProtein() );

        // all of the different tokens to match "perturbable"
        pert_token.add( "perturbable", AST::IDPerturbable() )
                      ( "PERTURBABLE", AST::IDPerturbable() )
                      ( "pert", AST::IDPerturbable() )
                      ( "PERT", AST::IDPerturbable() );

        //all of the different object names
        name_token.add( "atomnam", AST::ATOM )
                      ( "atomname", AST::ATOM )
                      ( "cgname", AST::CUTGROUP )
                      ( "cgnam", AST::CUTGROUP )
                      ( "resnam", AST::RESIDUE )
                      ( "resname", AST::RESIDUE )
                      ( "chainnam", AST::CHAIN )
                      ( "chainname", AST::CHAIN )
                      ( "segnam", AST::SEGMENT )
                      ( "segname", AST::SEGMENT )
                      ( "molnam", AST::MOLECULE )
                      ( "molname", AST::MOLECULE );

        //all of the different object numbers
        number_token
          .add( "atomnum", QPair<AST::IDObject,AST::IDNumType>(AST::ATOM,AST::ID_NUMBER) )
              ( "atomidx", QPair<AST::IDObject,AST::IDNumType>(AST::ATOM,AST::ID_INDEX) )
              ( "cgnum", QPair<AST::IDObject,AST::IDNumType>(AST::CUTGROUP,AST::ID_NUMBER) )
              ( "cgidx", QPair<AST::IDObject,AST::IDNumType>(AST::CUTGROUP,AST::ID_INDEX) )
              ( "resnum", QPair<AST::IDObject,AST::IDNumType>(AST::RESIDUE,AST::ID_NUMBER) )
              ( "residx", QPair<AST::IDObject,AST::IDNumType>(AST::RESIDUE,AST::ID_INDEX) )
              ( "chainnum", QPair<AST::IDObject,AST::IDNumType>(AST::CHAIN,AST::ID_NUMBER) )
              ( "chainidx", QPair<AST::IDObject,AST::IDNumType>(AST::CHAIN,AST::ID_INDEX) )
              ( "segnum", QPair<AST::IDObject,AST::IDNumType>(AST::SEGMENT,AST::ID_NUMBER) )
              ( "segidx", QPair<AST::IDObject,AST::IDNumType>(AST::SEGMENT,AST::ID_INDEX) )
              ( "molnum", QPair<AST::IDObject,AST::IDNumType>(AST::MOLECULE,AST::ID_NUMBER) )
              ( "molidx", QPair<AST::IDObject,AST::IDNumType>(AST::MOLECULE,AST::ID_INDEX) )
              ;

        //all of the different types of logical operation
        op_token.add( "and", AST::ID_AND )
                    ( "AND", AST::ID_AND )
                    ( "or", AST::ID_OR )
                    ( "OR", AST::ID_OR );

        //all of the different value comparison tokens
        cmp_token.add( "<=", AST::ID_CMP_LE )
                     ( "<", AST::ID_CMP_LT )
                     ( "==", AST::ID_CMP_EQ )
                     ( "!=", AST::ID_CMP_NE )
                     ( ">=", AST::ID_CMP_GE )
                     ( ">", AST::ID_CMP_GT )
                     ( "=~", AST::ID_CMP_AE );

        //all of the different object identification tokens
        obj_token.add( "atoms",  AST::ATOM )
                     ( "atom", AST::ATOM )
                     ( "cutgroups", AST::CUTGROUP )
                     ( "cutgroup", AST::CUTGROUP )
                     ( "group", AST::CUTGROUP )
                     ( "groups", AST::CUTGROUP )
                     ( "residues", AST::RESIDUE )
                     ( "residue", AST::RESIDUE )
                     ( "res", AST::RESIDUE )
                     ( "chains", AST::CHAIN )
                     ( "chain", AST::CHAIN )
                     ( "segments", AST::SEGMENT )
                     ( "segment", AST::SEGMENT )
                     ( "segs", AST::SEGMENT )
                     ( "seg", AST::SEGMENT )
                     ( "molecules", AST::MOLECULE )
                     ( "molecule", AST::MOLECULE )
                     ( "mol", AST::MOLECULE )
                     ( "mols", AST::MOLECULE )
                     ( "bond", AST::BOND )
                     ( "bonds", AST::BOND )
                    ;

        //all of the different length unit tokens
        length_token.add( "Angstroms", SireUnits::angstrom )
                        ( "Angstrom", SireUnits::angstrom )
                        ( "angstroms", SireUnits::angstrom )
                        ( "angstrom", SireUnits::angstrom )
                        ( "A", SireUnits::angstrom )
                        ( "picometers", SireUnits::picometer )
                        ( "picometer", SireUnits::picometer )
                        ( "pm", SireUnits::picometer )
                        ( "nanometers", SireUnits::nanometer )
                        ( "nanometer", SireUnits::nanometer )
                        ( "nm", SireUnits::nanometer )
                        ;

        //all of the different mass tokens
        mass_token.add( "g_per_mol", SireUnits::g_per_mol )
                      ( "mg_per_mol", SireUnits::mg_per_mol )
                      ( "kg_per_mol", SireUnits::kg_per_mol )
                      ;

        //all of the different charge tokens
        charge_token.add( "e", SireUnits::mod_electron )
                        ( "coulomb", SireUnits::coulomb )
                        ;

        //all of the different "with" and "in" expression tokens
        with_token.add( "with", AST::ID_WITH )
                      ( "in", AST::ID_IN )
                    ;

        //all of the different bond tokens
        bond_token.add( "to", AST::ID_BOND_TO )
                      ( "from", AST::ID_BOND_FROM )
                    ;

        //all of the different types of coordinates tokens
        coord_token.add( "center", AST::ID_COORD_CENTER )
                       ( "coords.center", AST::ID_COORD_CENTER )
                       ( "center.x", AST::ID_COORD_CENTER_X )
                       ( "coords.center.x", AST::ID_COORD_CENTER_X )
                       ( "center.y", AST::ID_COORD_CENTER_Y )
                       ( "coords.center.y", AST::ID_COORD_CENTER_Y )
                       ( "center.z", AST::ID_COORD_CENTER_Z )
                       ( "coords.center.z", AST::ID_COORD_CENTER_Z )
                       ( "max", AST::ID_COORD_MAX )
                       ( "coords.max", AST::ID_COORD_MAX )
                       ( "max.x", AST::ID_COORD_MAX_X )
                       ( "coords.max.x", AST::ID_COORD_MAX_X )
                       ( "max.y", AST::ID_COORD_MAX_Y )
                       ( "coords.max.y", AST::ID_COORD_MAX_Y )
                       ( "max.z", AST::ID_COORD_MAX_Z )
                       ( "coords.max.z", AST::ID_COORD_MAX_Z )
                       ( "min", AST::ID_COORD_MIN )
                       ( "coords.min", AST::ID_COORD_MIN )
                       ( "min.x", AST::ID_COORD_MIN_X )
                       ( "coords.min.x", AST::ID_COORD_MIN_X )
                       ( "min.y", AST::ID_COORD_MIN_Y )
                       ( "coords.min.y", AST::ID_COORD_MIN_Y )
                       ( "min.z", AST::ID_COORD_MIN_Z )
                       ( "coords.min.z", AST::ID_COORD_MIN_Z )
                       ( "x", AST::ID_COORD_X )
                       ( "coords.x", AST::ID_COORD_X )
                       ( "y", AST::ID_COORD_Y )
                       ( "coords.y", AST::ID_COORD_Y )
                       ( "z", AST::ID_COORD_Z )
                       ( "coords.z", AST::ID_COORD_Z )
                    ;

        //now add in all of the element tokens
        for (int i=0; i<=111; ++i)  //loop through all known elements
        {
            Element e(i);

            //add tokens for the capitalised symbol, and lowercase symbol and name
            element_token.add( e.symbol().toLatin1().constData(), e );
            element_token.add( e.symbol().toLower().toLatin1().constData(), e );
            element_token.add( e.name().toLower().toLatin1().constData(), e );
        }

        //now get all of the user tokens (user-identified sub-expressions)
        user_token = getUserTokens();

        /////
        ///// Now define all of the grammar rules
        /////

        //root rule to read a node as a set of expressions
        nodeRule %= expressionsRule;

        //a set of expressions is a list of expression rules separated by semicolons
        expressionsRule %= ( expressionRule % qi::lit( ';' ) );

        //an expression is either a binary or a expression
        expressionRule %= binaryRule2 | binaryRule | withRule2 | withRule | expressionPartRule;

        //a binary is two expressions separated by an op_token (and/or)
        binaryRule %= (expressionPartRule >> op_token >> expressionPartRule) |
                      ( qi::lit('(') >> binaryRule >> qi::lit(')') );

        //allow multiple op_tokens, e.g. a and b and c
        binaryRule2 %= binaryRule >> op_token >> binaryRule |
                       binaryRule >> op_token >> expressionPartRule |
                       expressionPartRule >> op_token >> binaryRule |
                       (qi::lit('(') >> binaryRule2 >> qi::lit(')') );

        //a withRule is two expressions separated by a "with" or "in"
        withRule %= (expressionPartRule >> with_token >> expressionPartRule) |
                    ( qi::lit('(') >> withRule >> qi::lit(')') );

        //allow multiple with_tokens, e.g. atoms in molecules with resname ALA
        withRule2 %= withRule >> with_token >> withRule |
                     expressionPartRule >> with_token >> withRule |
                     withRule >> with_token >> expressionPartRule |
                     (qi::lit('(') >> withRule2 >> qi::lit(')') );

        //an expression is either a subscript, name, number, within, where, not
        //or user-identified expression, optionally surrounded by parenthesis '( )'
        expressionPartRule %= subscriptRule | idNameRule | idNumberRule | idElementRule |
                              propertyRule | bondRule | all_token | water_token | pert_token | protein_token | withinRule |
                              withinVectorRule | whereRule | notRule | joinRule |
                              massRule | massCmpRule | chargeRule | chargeCmpRule |
                              massObjRule | massObjCmpRule | chargeObjRule | chargeObjCmpRule |
                              countRule | user_token |
                              ( qi::lit('(') >> expressionPartRule >> qi::lit(')') );

        //grammar that specifies a list of names (comma-separated)
        nameValuesRule %= ( nameValueRule % qi::lit( ',' ) );

        //grammar for a single name (string or regular expression)
        nameValueRule %= regExpRule | stringRule;

        //grammar for a regular expression (identified using '/')
        regExpRule = eps [ _val = AST::RegExpValue() ] >>
                     (
                        lexeme[ "/" >> as_string[+(char_ - "/")][ _val += _1 ] >> "/" ]
                        >> -qi::lit("i")[ _val *= 1 ]
                     )
                     ;

        //grammar for a set of integers (either as ranges or comparisons)
        rangeValuesRule %= ( (compareValueRule | rangeValueRule) % qi::lit( ',' ) );

        //grammar for an integer or range (e.g. 0:10, or 5)
        rangeValueRule = eps [ _val = AST::RangeValue() ] >>
                            (
                                -(int_[ _val += _1 ]) >>
                                -(qi::lit(":")[ _val *= 1 ]) >>
                                -(int_[ _val += _1 ]) >>
                                -(qi::lit(":")[ _val *= 1 ]) >>
                                -(int_[ _val += _1 ])
                            )
                            ;

        massValueRule = eps [ _val = AST::IDMass() ] >>
                            (
                                double_[ _val += _1 ] >>
                                mass_token[ _val += _1 ]
                            )
                            |
                            (
                                double_[ _val += _1 ]
                            )
                            ;

        chargeValueRule = eps [ _val = AST::IDCharge() ] >>
                              (
                                  double_[ _val += _1 ] >>
                                  charge_token[ _val += _1 ]
                              )
                              |
                              (
                                  double_[ _val += _1 ]
                              )
                              ;

        //allow looking for mass
        massRule %= qi::lit("mass") >> massValueRule;
        massCmpRule %= qi::lit("mass") >> cmp_token >> massValueRule;
        massObjRule %= obj_token >> qi::lit("mass") >> massValueRule;
        massObjCmpRule %= obj_token >> qi::lit("mass") >> cmp_token >> massValueRule;

        //allow looking for charge
        chargeRule %= qi::lit("charge") >> chargeValueRule;
        chargeCmpRule %= qi::lit("charge") >> cmp_token >> chargeValueRule;
        chargeObjRule %= obj_token >> qi::lit("charge") >> chargeValueRule;
        chargeObjCmpRule %= obj_token >> qi::lit("charge") >> cmp_token >> chargeValueRule;

        //grammar for a comparison (e.g. x > 5)
        compareValueRule %= cmp_token >> int_;

        //grammar for a length/distance (with optional unit)
        lengthValueRule = eps [ _val = AST::LengthValue() ] >>
                            (
                                double_[ _val += _1 ] >>
                                length_token[ _val += _1 ]
                            )
                            |
                            (
                                double_[ _val += _1 ]
                            )
                            ;

        //grammar for a vector/point (with optional units, and optionally in brackets '(')
        vectorValueRule = eps[ _val = AST::VectorValue() ] >>
                            (
                                lengthValueRule[ _val += _1 ] >>
                                qi::repeat(0,2)[( ',' >> lengthValueRule[ _val += _1 ] )]
                            )
                            |
                            (
                                qi::lit('(') >>
                                lengthValueRule[ _val += _1 ] >>
                                qi::repeat(0,2)[( ',' >> lengthValueRule[ _val += _1 ] )] >>
                                qi::lit(')')
                            )
                            ;

        //grammar for an individual name assigned to name values
        idNameRule  %= name_token >> nameValuesRule;

        //grammer for an individual numbers assigned to number values
        idNumberRule = eps [ _val = AST::IDNumber() ] >>
                            (
                                number_token[ _val += _1 ] >>
                                rangeValuesRule[ _val += _1 ]
                            )
                            ;

        //grammar for selecting by chemical element
        idElementRule %= qi::lit("element") >> ( element_token % qi::lit(",") );

        //allow searching by molecular property
        propertyRule = eps [ _val = AST::IDProperty() ] >>
                            (
                                qi::lit("property")[ _val /= 1 ] >>
                                stringRule[ _val += _1 ] >>
                                cmp_token[ _val += _1 ] >>
                                stringRule[ _val *= _1 ]
                            )
                            |
                            (
                                qi::lit("property")[ _val /= 1 ] >>
                                stringRule[ _val += _1 ]
                            )
                            |
                            (
                                obj_token[ _val += _1 ] >>
                                qi::lit("property") >>
                                stringRule[ _val += _1 ] >>
                                cmp_token[ _val += _1 ] >>
                                stringRule[ _val *= _1 ]
                            )
                            |
                            (
                                obj_token[ _val += _1 ] >>
                                qi::lit("property") >>
                                stringRule[ _val += _1 ]
                            )
                            ;

        //allow looking for bonds
        bondRule %= (qi::lit("bonds") >> bond_token >> expressionRule
                                      >> bond_token >> expressionRule) |
                    (qi::lit("bonds") >> bond_token >> expressionRule);

        //grammar for a "not" expression
        notRule %= qi::lit("not") >> expressionRule;

        //grammar for a "join" expression
        joinRule %= qi::lit("join") >> expressionRule;

        //grammar for a "within" expression
        withinRule %= obj_token >> qi::lit("within") >> lengthValueRule
                                >> qi::lit("of") >> expressionRule;

        //grammar for a "within" expression comparing with a vector position.
        withinVectorRule %= obj_token >> qi::lit("within") >> lengthValueRule
                                      >> qi::lit("of") >> vectorValueRule;

        //grammar to enable subscripting of an expression
        subscriptRule %= qi::lit("{") >> expressionRule >> qi::lit("}") >>
                         qi::lit("[") >> rangeValueRule >> qi::lit("]");

        //grammar for a "where" expression
        whereRule %= obj_token >> qi::lit("where") >> coord_token >>
                        (whereWithinRule | whereCompareRule);

        //sub-grammar for a "where within" expression
        whereWithinRule %= qi::lit("within") >> lengthValueRule >> qi::lit("of")
                              >> expressionRule;

        //sub-grammar for a "where comparison" expression
        whereCompareRule %= cmp_token >> vectorValueRule;

        //grammar for a count expression e.g. count(atoms) < 5
        countRule = eps [ _val = AST::IDCount() ] >>
                            (
                                qi::lit("count(") >>
                                expressionRule[ _val += _1 ] >>
                                qi::lit(")") >>
                                cmp_token[ _val += _1 ] >>
                                qi::int_[ _val += _1 ]
                            )
                            ;

        /////
        ///// name all of the rules to simplify error messages
        /////
        nodeRule.name( "Node" );
        idNameRule.name( "Name" );
        idNumberRule.name( "Number" );
        idElementRule.name( "Element" );
        binaryRule.name( "Binary" );
        binaryRule2.name( "Binary2" );
        withRule.name( "With" );
        withRule2.name( "With2" );
        withinRule.name( "Within" );
        withinVectorRule.name( "Within Vector" );
        notRule.name( "Not" );
        joinRule.name( "Join" );
        subscriptRule.name( "Subscript" );
        whereRule.name( "Where" );
        whereWithinRule.name( "Where Within" );
        whereCompareRule.name( "Where Compare" );
        countRule.name( "Count Rule" );
        expressionsRule.name( "Expressions" );
        expressionRule.name( "Expression" );
        expressionPartRule.name( "Expression Part" );
        nameValuesRule.name( "Name Values" );
        nameValueRule.name( "Name Value" );
        rangeValuesRule.name( "Range Values" );
        compareValueRule.name( "Compare Value" );
        rangeValueRule.name( "Range Value" );
        lengthValueRule.name( "Length Value" );
        vectorValueRule.name( "Vector Value" );
        massValueRule.name( "Mass Value" );
        chargeValueRule.name( "Charge Value" );
        massRule.name( "Mass Rule" );
        massCmpRule.name( "Mass Compare Rule" );
        massObjRule.name( "Mass Object Rule" );
        massObjCmpRule.name( "Mass Object Compare Rule" );
        chargeRule.name( "Charge Rule" );
        chargeCmpRule.name( "Charge Compare Rule" );
        chargeObjRule.name( "Charge Object Rule" );
        chargeObjCmpRule.name( "Charge Object Compare Rule" );
        stringRule.name( "String" );
        regExpRule.name( "RegExp" );
        bondRule.name( "Bond" );
        propertyRule.name( "Property" );

        //action on failure to parse the string using the grammar
        on_error<fail>
        (
            nodeRule
          , std::cout
                << val("Error! Expecting ")
                << _4                               // what failed?
                << val(" here: \"")
                << construct<std::string>(_3, _2)   // iterators to error-pos, end
                << val("\"")
                << std::endl
        );
    }

    qi::rule<IteratorT, AST::Node(), SkipperT> nodeRule;
    qi::rule<IteratorT, AST::IDName(), SkipperT> idNameRule;
    qi::rule<IteratorT, AST::IDNumber(), SkipperT> idNumberRule;
    qi::rule<IteratorT, AST::IDElement(), SkipperT> idElementRule;
    qi::rule<IteratorT, AST::IDBinary(), SkipperT> binaryRule;
    qi::rule<IteratorT, AST::IDBinary(), SkipperT> binaryRule2;
    qi::rule<IteratorT, AST::IDBond(), SkipperT> bondRule;
    qi::rule<IteratorT, AST::IDProperty(), SkipperT> propertyRule;
    qi::rule<IteratorT, AST::IDWith(), SkipperT> withRule;
    qi::rule<IteratorT, AST::IDWith(), SkipperT> withRule2;
    qi::rule<IteratorT, AST::IDWithin(), SkipperT> withinRule;
    qi::rule<IteratorT, AST::IDWithinVector(), SkipperT> withinVectorRule;
    qi::rule<IteratorT, AST::IDNot(), SkipperT> notRule;
    qi::rule<IteratorT, AST::IDJoin(), SkipperT> joinRule;
    qi::rule<IteratorT, AST::IDSubscript(), SkipperT> subscriptRule;
    qi::rule<IteratorT, AST::IDMass(), SkipperT> massRule;
    qi::rule<IteratorT, AST::IDCmpMass(), SkipperT> massCmpRule;
    qi::rule<IteratorT, AST::IDObjMass(), SkipperT> massObjRule;
    qi::rule<IteratorT, AST::IDObjCmpMass(), SkipperT> massObjCmpRule;
    qi::rule<IteratorT, AST::IDCharge(), SkipperT> chargeRule;
    qi::rule<IteratorT, AST::IDCmpCharge(), SkipperT> chargeCmpRule;
    qi::rule<IteratorT, AST::IDObjCharge(), SkipperT> chargeObjRule;
    qi::rule<IteratorT, AST::IDObjCmpCharge(), SkipperT> chargeObjCmpRule;

    qi::rule<IteratorT, AST::IDWhere(), SkipperT> whereRule;
    qi::rule<IteratorT, AST::IDWhereWithin(), SkipperT> whereWithinRule;
    qi::rule<IteratorT, AST::IDWhereCompare(), SkipperT> whereCompareRule;
    qi::rule<IteratorT, AST::IDCount(), SkipperT> countRule;

    qi::rule<IteratorT, AST::Expressions(), SkipperT> expressionsRule;
    qi::rule<IteratorT, AST::Expression(), SkipperT> expressionRule;

    qi::rule<IteratorT, AST::ExpressionPart(), SkipperT> expressionPartRule;

    qi::rule<IteratorT, AST::NameValues(), SkipperT> nameValuesRule;
    qi::rule<IteratorT, AST::NameValue(), SkipperT> nameValueRule;

    qi::rule<IteratorT, AST::RangeValues(), SkipperT> rangeValuesRule;
    qi::rule<IteratorT, AST::CompareValue(), SkipperT> compareValueRule;
    qi::rule<IteratorT, AST::RangeValue(), SkipperT> rangeValueRule;

    qi::rule<IteratorT, AST::IDMass(), SkipperT> massValueRule;
    qi::rule<IteratorT, AST::IDCharge(), SkipperT> chargeValueRule;

    qi::rule<IteratorT, AST::LengthValue(), SkipperT> lengthValueRule;
    qi::rule<IteratorT, AST::VectorValue(), SkipperT> vectorValueRule;

    qi::symbols<char,AST::IDObject> name_token;
    qi::symbols<char,QPair<AST::IDObject,AST::IDNumType> > number_token;
    qi::symbols<char,AST::IDOperation> op_token;
    qi::symbols<char,AST::IDObject> obj_token;
    qi::symbols<char,AST::IDToken> with_token;
    qi::symbols<char,AST::IDBondToken> bond_token;
    qi::symbols<char,SireUnits::Dimension::Length> length_token;
    qi::symbols<char,SireUnits::Dimension::MolarMass> mass_token;
    qi::symbols<char,SireUnits::Dimension::Charge> charge_token;
    qi::symbols<char,AST::IDComparison> cmp_token;
    qi::symbols<char,AST::IDCoordType> coord_token;
    qi::symbols<char,SireMol::Element> element_token;
    qi::symbols<char,AST::IDAll> all_token;
    qi::symbols<char,AST::IDWater> water_token;
    qi::symbols<char,AST::IDPerturbable> pert_token;
    qi::symbols<char,AST::IDProtein> protein_token;
    UserTokens user_token;

    ValueGrammar<IteratorT, SkipperT> stringRule;
    qi::rule<IteratorT, AST::RegExpValue(), SkipperT> regExpRule;
};

#endif

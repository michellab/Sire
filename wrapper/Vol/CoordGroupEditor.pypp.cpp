// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CoordGroupEditor.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/quickcopy.hpp"

#include "SireError/errors.h"

#include "SireMaths/align.h"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "coordgroup.h"

#include <QDebug>

#include "coordgroup.h"

SireVol::CoordGroupEditor __copy__(const SireVol::CoordGroupEditor &other){ return SireVol::CoordGroupEditor(other); }

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_CoordGroupEditor_class(){

    { //::SireVol::CoordGroupEditor
        typedef bp::class_< SireVol::CoordGroupEditor, bp::bases< SireVol::CoordGroupBase > > CoordGroupEditor_exposer_t;
        CoordGroupEditor_exposer_t CoordGroupEditor_exposer = CoordGroupEditor_exposer_t( "CoordGroupEditor", "This class is used to edit a CoordGroup. This class is used when you want to\nmake several small changes to a CoordGroup, but do not want the CoordGroup to\nupdate its internal state after each change (e.g. you are moving each point in\nturn, and do not want the AABox to be updated for every step)\n\nYou use a CoordGroupEditor like this;\n\node\n\ncreate a CoordGroup with space for 100 coordinates\nCoordGroup coordgroup(100);\n\ncreate an editor for this group\nCoordGroupEditor editor = coordgroup.edit();\n\nmanipulate each coordinate in turn\nfor (int i=0; i<100; ++i)\neditor[i] = Vector(i,i+1,i+2);\n\ncommit the changes\ncoordgroup = editor.commit();\n\nE:ndcode\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope CoordGroupEditor_scope( CoordGroupEditor_exposer );
        CoordGroupEditor_exposer.def( bp::init< SireVol::CoordGroup const & >(( bp::arg("other") ), "Construct from a CoordGroup") );
        CoordGroupEditor_exposer.def( bp::init< SireVol::CoordGroupEditor const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireVol::CoordGroupEditor::changeFrame
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*changeFrame_function_type)( ::SireMaths::AxisSet const &,::SireMaths::AxisSet const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireVol::CoordGroupEditor::changeFrame );
            
            CoordGroupEditor_exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("from_frame"), bp::arg("to_frame") )
                , bp::return_self< >()
                , "Change the coordinate frame of all of the coordinates from\nfrom_frame to to_frame. This has the same effect\nas mapInto if from_frame is the unit cartesian set" );
        
        }
        { //::SireVol::CoordGroupEditor::changeFrame
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*changeFrame_function_type)( ::quint32,::SireMaths::AxisSet const &,::SireMaths::AxisSet const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireVol::CoordGroupEditor::changeFrame );
            
            CoordGroupEditor_exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("i"), bp::arg("from_frame"), bp::arg("to_frame") )
                , bp::return_self< >()
                , "Map the coordinates of the ith point from the frame\nfrom_frame to the frame to_frame\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::commit
        
            typedef ::SireVol::CoordGroup ( ::SireVol::CoordGroupEditor::*commit_function_type)(  ) ;
            commit_function_type commit_function_value( &::SireVol::CoordGroupEditor::commit );
            
            CoordGroupEditor_exposer.def( 
                "commit"
                , commit_function_value
                , "Return a CoordGroup which is a copy of this group." );
        
        }
        { //::SireVol::CoordGroupEditor::mapInto
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*mapInto_function_type)( ::SireMaths::AxisSet const & ) ;
            mapInto_function_type mapInto_function_value( &::SireVol::CoordGroupEditor::mapInto );
            
            CoordGroupEditor_exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("axes") )
                , bp::return_self< >()
                , "Map all of these coordinates into the axis set axes" );
        
        }
        { //::SireVol::CoordGroupEditor::mapInto
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*mapInto_function_type)( ::quint32,::SireMaths::AxisSet const & ) ;
            mapInto_function_type mapInto_function_value( &::SireVol::CoordGroupEditor::mapInto );
            
            CoordGroupEditor_exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("i"), bp::arg("axes") )
                , bp::return_self< >()
                , "Map the coordinates of the ith point into the axis set axes\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::operator=
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*assign_function_type)( ::SireVol::CoordGroup const & ) ;
            assign_function_type assign_function_value( &::SireVol::CoordGroupEditor::operator= );
            
            CoordGroupEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("cgroup") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireVol::CoordGroupEditor::operator=
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*assign_function_type)( ::SireVol::CoordGroupEditor const & ) ;
            assign_function_type assign_function_value( &::SireVol::CoordGroupEditor::operator= );
            
            CoordGroupEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireVol::CoordGroupEditor::operator[]
        
            typedef ::SireMaths::Vector & ( ::SireVol::CoordGroupEditor::*__getitem___function_type)( ::quint32 ) ;
            __getitem___function_type __getitem___function_value( &::SireVol::CoordGroupEditor::operator[] );
            
            CoordGroupEditor_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_internal_reference< >()
                , "" );
        
        }
        { //::SireVol::CoordGroupEditor::rotate
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*rotate_function_type)( ::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupEditor::rotate );
            
            CoordGroupEditor_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("point") )
                , bp::return_self< >()
                , "Rotate this group by the Quaternion quat about the point point" );
        
        }
        { //::SireVol::CoordGroupEditor::rotate
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*rotate_function_type)( ::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupEditor::rotate );
            
            CoordGroupEditor_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("point") )
                , bp::return_self< >()
                , "Rotate (and perhaps scale and shear) this group by the matrix rotmat\nabout the point point" );
        
        }
        { //::SireVol::CoordGroupEditor::rotate
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*rotate_function_type)( ::quint32,::SireMaths::Quaternion const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupEditor::rotate );
            
            CoordGroupEditor_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("i"), bp::arg("quat"), bp::arg("point") )
                , bp::return_self< >()
                , "Rotate the ith point in the group using the quaternion quat about the\npoint point\nThrow: SireError::index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::rotate
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*rotate_function_type)( ::quint32,::SireMaths::Matrix const &,::SireMaths::Vector const & ) ;
            rotate_function_type rotate_function_value( &::SireVol::CoordGroupEditor::rotate );
            
            CoordGroupEditor_exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("i"), bp::arg("rotmat"), bp::arg("point") )
                , bp::return_self< >()
                , "Rotate the ith point in the group using the matrix rotmat about the\npoint point\nThrow: SireError::index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::setCoordinates
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*setCoordinates_function_type)( ::QVector< SireMaths::Vector > const & ) ;
            setCoordinates_function_type setCoordinates_function_value( &::SireVol::CoordGroupEditor::setCoordinates );
            
            CoordGroupEditor_exposer.def( 
                "setCoordinates"
                , setCoordinates_function_value
                , ( bp::arg("newcoords") )
                , bp::return_self< >()
                , "Set the coordinates of the CoordGroup to newcoords - this must\nhave the same number of points as this CoordGroup or an exception\nwill be thrown\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupEditor::setCoordinates
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*setCoordinates_function_type)( ::SireVol::CoordGroupBase const & ) ;
            setCoordinates_function_type setCoordinates_function_value( &::SireVol::CoordGroupEditor::setCoordinates );
            
            CoordGroupEditor_exposer.def( 
                "setCoordinates"
                , setCoordinates_function_value
                , ( bp::arg("newcoords") )
                , bp::return_self< >()
                , "Set the coordinates of the CoordGroup to newcoords - this\nmust have the same number of points as this CoordGroup or an\nexception will be thrown\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireVol::CoordGroupEditor::setCoordinates
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*setCoordinates_function_type)( ::quint32,::SireMaths::Vector const & ) ;
            setCoordinates_function_type setCoordinates_function_value( &::SireVol::CoordGroupEditor::setCoordinates );
            
            CoordGroupEditor_exposer.def( 
                "setCoordinates"
                , setCoordinates_function_value
                , ( bp::arg("i"), bp::arg("newcoords") )
                , bp::return_self< >()
                , "Set the coordinates of the ith atom to newcoords\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::transform
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*transform_function_type)( ::SireMaths::Transform const & ) ;
            transform_function_type transform_function_value( &::SireVol::CoordGroupEditor::transform );
            
            CoordGroupEditor_exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("t") )
                , bp::return_self< >()
                , "Transform this group by the transformation t" );
        
        }
        { //::SireVol::CoordGroupEditor::transform
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*transform_function_type)( ::quint32,::SireMaths::Transform const & ) ;
            transform_function_type transform_function_value( &::SireVol::CoordGroupEditor::transform );
            
            CoordGroupEditor_exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("i"), bp::arg("t") )
                , bp::return_self< >()
                , "Rotate the ith point in the group using the matrix rotmat about the\npoint point\nThrow: SireError::index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::translate
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*translate_function_type)( ::SireMaths::Vector const & ) ;
            translate_function_type translate_function_value( &::SireVol::CoordGroupEditor::translate );
            
            CoordGroupEditor_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta") )
                , bp::return_self< >()
                , "Translate this CoordGroup by delta" );
        
        }
        { //::SireVol::CoordGroupEditor::translate
        
            typedef ::SireVol::CoordGroupEditor & ( ::SireVol::CoordGroupEditor::*translate_function_type)( ::quint32,::SireMaths::Vector const & ) ;
            translate_function_type translate_function_value( &::SireVol::CoordGroupEditor::translate );
            
            CoordGroupEditor_exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("i"), bp::arg("delta") )
                , bp::return_self< >()
                , "Translate the ith point in the group by delta\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireVol::CoordGroupEditor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireVol::CoordGroupEditor::typeName );
            
            CoordGroupEditor_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireVol::CoordGroupEditor::what
        
            typedef char const * ( ::SireVol::CoordGroupEditor::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireVol::CoordGroupEditor::what );
            
            CoordGroupEditor_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        CoordGroupEditor_exposer.staticmethod( "typeName" );
        CoordGroupEditor_exposer.def( "__copy__", &__copy__);
        CoordGroupEditor_exposer.def( "__deepcopy__", &__copy__);
        CoordGroupEditor_exposer.def( "clone", &__copy__);
        CoordGroupEditor_exposer.def( "__str__", &__str__< ::SireVol::CoordGroupEditor > );
        CoordGroupEditor_exposer.def( "__repr__", &__str__< ::SireVol::CoordGroupEditor > );
        CoordGroupEditor_exposer.def( "__len__", &__len_size< ::SireVol::CoordGroupEditor > );
    }

}

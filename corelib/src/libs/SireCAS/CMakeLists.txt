########################################
#
# CMake file for library:SireCAS
#
########################################

# Other Sire libraries
include_directories(${CMAKE_SOURCE_DIR}/src/libs)

# Define the headers in SireCAS
set ( SIRECAS_HEADERS
      abs.h
      complexvalues.h
      conditional.h
      constant.h
      doublefunc.h
      errors.h
      exbase.h
      exp.h
      expressionbase.h
      expression.h
      expressionproperty.h
      expressions.h
      function.h
      functions.h
      functionsignature.h
      hyperbolicfuncs.h
      i.h
      identities.h
      integrationconstant.h
      invhyperbolicfuncs.h
      invtrigfuncs.h
      minmax.h
      powerconstant.h
      power.h
      product.h
      singlefunc.h
      sum.h
      symbol.h
      symbols.h
      symbolcomplex.h
      symbolexpression.h
      symbolvalue.h
      trigfuncs.h
      values.h
    )

# Define the sources in SireCAS
set ( SIRECAS_SOURCES

      register_sirecas.cpp

      abs.cpp
      complexvalues.cpp           
      conditional.cpp
      constant.cpp
      doublefunc.cpp                
      errors.cpp                  
      exbase.cpp                  
      exp.cpp                     
      expressionbase.cpp          
      expression.cpp              
      expressionproperty.cpp
      expressions.cpp             
      function.cpp                
      functions.cpp               
      functionsignature.cpp       
      hyperbolicfuncs.cpp         
      i.cpp                       
      identities.cpp              
      integrationconstant.cpp     
      invhyperbolicfuncs.cpp      
      invtrigfuncs.cpp 
      minmax.cpp           
      powerconstant.cpp           
      power.cpp                   
      product.cpp                 
      singlefunc.cpp              
      sum.cpp                     
      symbol.cpp                  
      symbolexpression.cpp        
      trigfuncs.cpp               
      values.cpp
      
      ${SIRECAS_HEADERS}
                        
    )

add_library (SireCAS ${SIRECAS_SOURCES})

set_target_properties (SireCAS PROPERTIES
                       VERSION ${SIRE_VERSION}
                       SOVERSION ${S_VERSION_MAJOR}
                      )

# Link to other Sire libraries, and GSL
target_link_libraries (SireCAS
                       SireMaths
                       SireBase
                       SireStream
                       )

# installation
install( TARGETS SireCAS EXPORT SireLibraries
         RUNTIME DESTINATION ${SIRE_BIN}
         LIBRARY DESTINATION ${SIRE_LIBS}
         ARCHIVE DESTINATION ${SIRE_ARCHIVES}
       )

install( FILES ${SIRECAS_HEADERS}
         DESTINATION ${SIRE_INCLUDES}/SireCAS )


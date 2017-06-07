// ==============================================================================
//      ___   __   ____             
//     / _ | / /  / __/             
//    / __ |/ /__/ _/               
//   /_/ |_/____/___/  application  
//
// License:		     BSD License
//					 license: ALEApplication/license.txt
//==============================================================================

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: April 2017 $
 * ***************************************************************************/

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "ale_application_variables.h"
#include "ale_application.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;



BOOST_PYTHON_MODULE(KratosALEApplication)
{

    class_<KratosALEApplication,
           KratosALEApplication::Pointer,
           bases<KratosApplication>, boost::noncopyable >("KratosALEApplication")
           ;

    AddCustomStrategiesToPython();
    AddCustomUtilitiesToPython();

    //registering variables in python

    //MESH_VELOCITY currently put to the core since used in other applications
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(MESH_VELOCITY)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);

}

}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined

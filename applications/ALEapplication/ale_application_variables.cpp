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

 #include "ale_application_variables.h"

 namespace Kratos
 {

   //KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_VELOCITY) //currently put to the core since used in other applications
   KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
   KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
   KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);

 }  // namespace Kratos.

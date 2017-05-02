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


 #if !defined(KRATOS_ALE_APPLICATION_VARIABLES_H_INCLUDED)
 #define  KRATOS_ALE_APPLICATION_VARIABLES_H_INCLUDED

 // System includes

 // External includes


 // Project includes
 #include "includes/define.h"
 #include "includes/variables.h"

 namespace Kratos
 {
   //MESH_VELOCITY currently put to the core since used in other applications
   //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_VELOCITY)
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
   KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);


 } // namespace Kratos.

 #endif // KRATOS_ALE_APPLICATION_VARIABLES_H_INCLUDED

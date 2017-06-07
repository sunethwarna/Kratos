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


// External includes


// Project includes
//#include "includes/define.h"
#include "ale_application.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/prism_3d_15.h"
#include "geometries/prism_3d_6.h"


namespace Kratos
{
//

KratosALEApplication::KratosALEApplication():
    mStructuralMeshMovingElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mStructuralMeshMovingElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mStructuralMeshMovingElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mStructuralMeshMovingElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mStructuralMeshMovingElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mStructuralMeshMovingElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) )
{}


void KratosALEApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "KRATOS    ___   __   ____             " << std::endl;
    std::cout << "         / _ | / /  / __/             " << std::endl;
    std::cout << "        / __ |/ /__/ _/               " << std::endl;
    std::cout << "       /_/ |_/____/___/  application  " << std::endl;
    std::cout << "Initializing KratosALEApplication...  " << std::endl;


    // register local flags
    

    // register elements
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D3N", mStructuralMeshMovingElement2D3N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement2D4N", mStructuralMeshMovingElement2D4N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D4N", mStructuralMeshMovingElement3D4N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D8N", mStructuralMeshMovingElement3D8N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D6N", mStructuralMeshMovingElement3D6N);
    KRATOS_REGISTER_ELEMENT("StructuralMeshMovingElement3D15N", mStructuralMeshMovingElement3D15N);


    // register variables
    //MESH_VELOCITY currently put to the core since used in other applications
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_VELOCITY);
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_DISPLACEMENT);
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_REACTION);
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_RHS);

}

}  // namespace Kratos.

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
#include <boost/python.hpp>


// Project includes
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/ball_vertex_meshmoving.h"
#include "custom_utilities/ball_vertex_meshmoving3D.h"

#include "custom_utilities/ale_flags.h"


namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;


    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_< BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("BallVertexMeshMoving2D", init<	>() )
    .def("ConstructSystem",&BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("BuildAndSolveSystem",&BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >::BuildAndSolveSystem)
    .def("ClearSystem",&BallVertexMeshMoving< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;


    class_< BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >,  boost::noncopyable >	("BallVertexMeshMoving3D", init<	>() )
    .def("ConstructSystem",&BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("BuildAndSolveSystem",&BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >::BuildAndSolveSystem)
    .def("ClearSystem",&BallVertexMeshMoving3D< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;

}





}  // namespace Python.

} // Namespace Kratos

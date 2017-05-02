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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/structural_meshmoving_strategy.h"

//linear solvers
#include "linear_solvers/linear_solver.h"



namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef SolvingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > BaseSolvingStrategyType;

    class_< StructuralMeshMovingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >,
            bases< BaseSolvingStrategyType >,  boost::noncopyable >
            ("StructuralMeshMovingStrategy",
             init<ModelPart&, LinearSolverType::Pointer, int, bool, bool >() )
            .def("MoveNodes",&StructuralMeshMovingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::MoveNodes)
            .def("UpdateReferenceMesh",&StructuralMeshMovingStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::UpdateReferenceMesh)
            ;
}

}  // namespace Python.

} // Namespace Kratos

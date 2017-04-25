//   
//   Project Name:        
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//

// External includes 
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp> 

// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//strategies
#include "solving_strategies/strategies/solving_strategy.h"
#include "custom_strategies/strategies/dam_newton_raphson_global_damage_strategy.hpp"

//builders and solvers

//schemes
#include "custom_strategies/schemes/incrementalupdate_static_smoothing_scheme.hpp"
#include "custom_strategies/schemes/incrementalupdate_static_damped_smoothing_scheme.hpp"
#include "custom_strategies/schemes/bossak_displacement_smoothing_scheme.hpp"
#include "custom_strategies/schemes/dam_UP_scheme.hpp"
#include "custom_strategies/schemes/dam_P_scheme.hpp"


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
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;

    //custom scheme types
    typedef IncrementalUpdateStaticSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticSmoothingSchemeType; 
    typedef IncrementalUpdateStaticDampedSmoothingScheme< SparseSpaceType, LocalSpaceType >  IncrementalUpdateStaticDampedSmoothingSchemeType;
    typedef BossakDisplacementSmoothingScheme< SparseSpaceType, LocalSpaceType >  BossakDisplacementSmoothingSchemeType;
    typedef DamUPScheme< SparseSpaceType, LocalSpaceType >  DamUPSchemeType;
    typedef DamPScheme< SparseSpaceType, LocalSpaceType >  DamPSchemeType;

    //custom strategy types
    typedef DamNewtonRaphsonGlobalDamageStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType > DamNewtonRaphsonGlobalDamageStrategyType;


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	// Schemes
    class_< IncrementalUpdateStaticSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("IncrementalUpdateStaticSmoothingScheme",
        init< >());

    class_< IncrementalUpdateStaticDampedSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("IncrementalUpdateStaticDampedSmoothingScheme",
        init< double, double >());

    class_< BossakDisplacementSmoothingSchemeType, bases< BaseSchemeType >, boost::noncopyable >("BossakDisplacementSmoothingScheme",
        init< double, double, double >());

	class_< DamUPSchemeType, bases< BaseSchemeType >,  boost::noncopyable >("DamUPScheme",
        init< double, double, double, double >());
    
    class_< DamPSchemeType, bases< BaseSchemeType >,  boost::noncopyable >("DamPScheme",
        init< double, double >());

	// Strategies
    class_< DamNewtonRaphsonGlobalDamageStrategyType, bases< BaseSolvingStrategyType >, boost::noncopyable >("DamNewtonRaphsonGlobalDamageStrategy", 
        init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer,
                BuilderAndSolverType::Pointer, int, bool, bool, bool >());
}

}  // namespace Python.
} // Namespace Kratos


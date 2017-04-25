//   
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:        Lorenzo Gracia $
//   Date:                $Date:              April 2017 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_DAM_NEWTON_RAPHSON_GLOBAL_DAMAGE_STRATEGY)
#define KRATOS_DAM_NEWTON_RAPHSON_GLOBAL_DAMAGE_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class DamNewtonRaphsonGlobalDamageStrategy : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(DamNewtonRaphsonGlobalDamageStrategy);
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    using MotherType::mpScheme;
    using MotherType::mpBuilderAndSolver;
    using MotherType::mpA; //Tangent matrix
    using MotherType::mpb; //Residual vector of iteration i
    using MotherType::mpDx; //Delta x of iteration i
    using MotherType::mMaxIterationNumber;
    using MotherType::mInitializeWasPerformed;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    DamNewtonRaphsonGlobalDamageStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {

        }

    //------------------------------------------------------------------------------------

    ///Destructor
    virtual ~DamNewtonRaphsonGlobalDamageStrategy() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void FinalizeSolutionStep()
	{
		KRATOS_TRY

        MotherType::FinalizeSolutionStep();
        
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;
        KRATOS_WATCH(mDx)

        mpBuilderAndSolver->BuildRHS(mpScheme, BaseType::GetModelPart(), mb);
        KRATOS_WATCH(mb)

        double NonLinearPotencialEnergy;
        double LinearPotencialEnergy;

        // Computing the potential energy using the damage stress

        NonLinearPotencialEnergy = 1.0;
        KRATOS_WATCH(NonLinearPotencialEnergy)

        // Computing the potential energy using effective stress (linear)
        LinearPotencialEnergy = 1.0;
        KRATOS_WATCH(LinearPotencialEnergy)

        // Compùting the global damage of the structure according to potential Energy (Hanganu, Onate...)
        double global_damage = 1.0 - (NonLinearPotencialEnergy/LinearPotencialEnergy);
        KRATOS_WATCH(global_damage)

		KRATOS_CATCH("")
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


protected:
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class DamNewtonRaphsonGlobalDamageStrategy

} // namespace Kratos

#endif // KRATOS_DAM_NEWTON_RAPHSON_GLOBAL_DAMAGE_STRATEGY  defined

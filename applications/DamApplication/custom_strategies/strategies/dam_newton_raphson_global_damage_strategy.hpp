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
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag) , mr_model_part(model_part)
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
        TSystemVectorType& mb1 = *mpb;
        TSystemVectorType& mb2 = *mpb;

        TSparseSpace::SetToZero(mb1);
        double NonLinearPotencialEnergy = 0.0;
        double LinearPotencialEnergy = 0.0 ;

        BaseType::GetModelPart().GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] = 1;

        // Computing the potential energy using the damage stress
        if (BaseType::GetModelPart().GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] == 1) 
        {

            const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
            ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();

            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
                it->Set(ACTIVE,false);
            }

            mr_model_part.GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] = 1;
            mpBuilderAndSolver->BuildRHS(mpScheme, mr_model_part, mb1);
            NonLinearPotencialEnergy = TSparseSpace::Dot(mDx, mb1);

            KRATOS_WATCH(NonLinearPotencialEnergy)
  
            BaseType::GetModelPart().GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] = 2;
        } 

        // Computing the potential energy using effective stress (linear)
        TSparseSpace::SetToZero(mb2);
        if (BaseType::GetModelPart().GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] == 2) 
        {
            mr_model_part.GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] = 2;
            mpBuilderAndSolver->BuildRHS(mpScheme, mr_model_part, mb2);
            LinearPotencialEnergy = TSparseSpace::Dot(mDx, mb2);

            KRATOS_WATCH(LinearPotencialEnergy)
        } 

        // Compùting the global damage of the structure according to potential Energy (Hanganu, Onate...)
        double global_damage = 1.0 - (NonLinearPotencialEnergy/LinearPotencialEnergy);
        KRATOS_WATCH(global_damage)

        // Reactivation of all conditions for the next step and setting the flag to 0
        const int nconditions = static_cast<int>(mr_model_part.Conditions().size());
        ModelPart::ConditionsContainerType::iterator cond_begin = mr_model_part.ConditionsBegin();
        for (int k = 0; k < nconditions; k++)
        {
            ModelPart::ConditionsContainerType::iterator it = cond_begin + k;
            it->Set(ACTIVE,true);
        }

        BaseType::GetModelPart().GetProcessInfo()[COMPUTE_GLOBAL_DAMAGE] = 0;

		KRATOS_CATCH("")
	}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    ModelPart& mr_model_part;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class DamNewtonRaphsonGlobalDamageStrategy

} // namespace Kratos

#endif // KRATOS_DAM_NEWTON_RAPHSON_GLOBAL_DAMAGE_STRATEGY  defined

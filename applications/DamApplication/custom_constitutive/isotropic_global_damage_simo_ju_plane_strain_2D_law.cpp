//
//   Project Name:        KratosDamApplication            $
//   Created by:          $Author:         Lorenzo Gracia $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// Project includes
#include "custom_constitutive/isotropic_global_damage_simo_ju_plane_strain_2D_law.hpp"

namespace Kratos
{

//Default Constructor
IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::IsotropicGlobalDamageSimoJuPlaneStrain2DLaw()
    : SimoJuLocalDamagePlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new IsotropicDamageFlowRule(mpYieldCriterion) );

}

//----------------------------------------------------------------------------------------

//Second Constructor
IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : SimoJuLocalDamagePlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw) {}

//----------------------------------------------------------------------------------------

//Copy Constructor
IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(const IsotropicGlobalDamageSimoJuPlaneStrain2DLaw& rOther)
    : SimoJuLocalDamagePlaneStrain2DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::~IsotropicGlobalDamageSimoJuPlaneStrain2DLaw() {}

//----------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::Clone() const
{
    IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::Pointer p_clone(new IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    Vector& rStrainVector = rValues.GetStrainVector();

    // Initialize main variables

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    // Computing the linear branch
    if(CurrentProcessInfo[COMPUTE_GLOBAL_DAMAGE]==2)
    {   
        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                //COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                this->CalculateLinearElasticMatrix(ConstitutiveMatrix,YoungModulus,PoissonCoefficient);
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();
                this->CalculateLinearElasticMatrix(ConstitutiveMatrix,YoungModulus,PoissonCoefficient);

                noalias(rStressVector) = prod(ConstitutiveMatrix,rStrainVector);

            }
        }
        else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();
            Matrix ConstitutiveMatrix(VoigtSize,VoigtSize);
	        noalias(ConstitutiveMatrix) = ZeroMatrix(VoigtSize,VoigtSize);
            this->CalculateLinearElasticMatrix(ConstitutiveMatrix,YoungModulus,PoissonCoefficient);

            noalias(rStressVector) = prod(ConstitutiveMatrix,rStrainVector);
            
        }
    }
    else
    {
        // ReturnMappingVariables
        FlowRule::RadialReturnVariables ReturnMappingVariables;
        ReturnMappingVariables.initialize();
        // Strain and Stress matrices
        const unsigned int Dim = this->WorkingSpaceDimension();
        Matrix AuxMatrix(Dim,Dim);
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
        ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
        ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
        // CharacteristicSize
        double CharacteristicSize = 1.0;
        this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
        ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if(Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                // COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector EffectiveStressVector(VoigtSize);

                this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,EffectiveStressVector,LinearElasticMatrix,rStrainVector);

                this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();
                
                this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
                
                this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);

            }
        }
        else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();
            
            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);

        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{    
    // Check
    rValues.CheckAllParameters();
    
    // Get values for the constitutive law
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    Vector& rStrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = rStrainVector.size();
    Vector EffectiveStressVector(VoigtSize);
    
    // Initialize main variables

    // LinearElasticMatrix
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);
    
    if(CurrentProcessInfo[COMPUTE_GLOBAL_DAMAGE] !=2)
    {
        // ReturnMappingVariables
        FlowRule::RadialReturnVariables ReturnMappingVariables;
        ReturnMappingVariables.initialize();
        // Strain and Stress matrices
        const unsigned int Dim = this->WorkingSpaceDimension();
        Matrix AuxMatrix(Dim,Dim);
        noalias(AuxMatrix) = MathUtils<double>::StrainVectorToTensor(rStrainVector);
        ReturnMappingVariables.StrainMatrix.resize(Dim,Dim,false);
        noalias(ReturnMappingVariables.StrainMatrix) = AuxMatrix;
        ReturnMappingVariables.TrialIsoStressMatrix.resize(Dim,Dim,false);
        // CharacteristicSize
        double CharacteristicSize = 1.0;
        this->CalculateCharacteristicSize(CharacteristicSize,rValues.GetElementGeometry());
        ReturnMappingVariables.CharacteristicSize = CharacteristicSize;

        if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
        {
            ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,false); // Restore sate variable = false
            
            this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,rStrainVector);
        }
        else // No convergence is achieved. Restore state variable to equilibrium
        {
            ReturnMappingVariables.Options.Set(FlowRule::RETURN_MAPPING_COMPUTED,true); // Restore sate variable = true
            
            this->UpdateInternalStateVariables(ReturnMappingVariables,EffectiveStressVector,LinearElasticMatrix,rStrainVector);
        }

        if(rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();
            
            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
        }
    }
}


} // Namespace Kratos
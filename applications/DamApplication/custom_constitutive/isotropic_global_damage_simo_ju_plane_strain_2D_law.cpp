//
//   Project Name:        KratosDamApplication            $
//   Created by:          $Author:         Lorenzo Gracia $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>
#include <cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/isotropic_global_damage_simo_ju_plane_strain_2D_law.hpp"

#include "dam_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::IsotropicGlobalDamageSimoJuPlaneStrain2DLaw()
    : IsotropicDamageSimoJuPlaneStrain2DLaw()
{
  mpHardeningLaw   = HardeningLaw::Pointer( new ExponentialDamageHardeningLaw() );
  mpYieldCriterion = YieldCriterion::Pointer( new SimoJuYieldCriterion(mpHardeningLaw) );
  mpFlowRule       = FlowRule::Pointer( new IsotropicDamageFlowRule(mpYieldCriterion) );

}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw)
    : IsotropicDamageSimoJuPlaneStrain2DLaw(pFlowRule, pYieldCriterion, pHardeningLaw)
{

}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(const IsotropicGlobalDamageSimoJuPlaneStrain2DLaw& rOther)
    : IsotropicDamageSimoJuPlaneStrain2DLaw(rOther)
{

}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::Clone() const
{
    IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::Pointer p_clone(new IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::~IsotropicGlobalDamageSimoJuPlaneStrain2DLaw()
{
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//*****************************MATERIAL RESPONSES*************************************
//************************************************************************************

void IsotropicGlobalDamageSimoJuPlaneStrain2DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{

    //a.-Check if the constitutive parameters are passed correctly to the law calculation
    CheckParameters(rValues);

    //b.- Get Values to compute the constitutive law:
    Flags& Options                        = rValues.GetOptions();
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();   
    Vector& rStrainVector                 = rValues.GetStrainVector();


    //0.- Initialize parameters
    FlowRule::RadialReturnVariables ReturnMappingVariables;
    ReturnMappingVariables.initialize(); //it has to be called at the start

    // Initialize variables from the process information
    ReturnMappingVariables.DeltaTime = CurrentProcessInfo[DELTA_TIME];
    
    if(CurrentProcessInfo[IMPLEX] == 1)	
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,true);
    else
        ReturnMappingVariables.Options.Set(FlowRule::IMPLEX_ACTIVE,false);

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

    //1.- Lame constants
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix LinearElasticMatrix (VoigtSize,VoigtSize);
    this->CalculateLinearElasticMatrix(LinearElasticMatrix,YoungModulus,PoissonCoefficient);


    if(Options.Is( ConstitutiveLaw::COMPUTE_STRAIN ))
    {
        //1.-Compute total deformation gradient
        const Matrix&   DeformationGradientF  = rValues.GetDeformationGradientF();

        //2.-Push-Forward Left Cauchy-Green tensor b to the new configuration
        Matrix LeftCauchyGreenMatrix = prod(DeformationGradientF,trans(DeformationGradientF));

        //3.-Almansi Strain: e= 0.5*(1-invFT*invF)
        this->CalculateAlmansiStrain(LeftCauchyGreenMatrix,rStrainVector);
    }

    //2.-Calculate Total Kirchhoff stress
    
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
            
            // This is for computing the linear response
            if ( CurrentProcessInfo[COMPUTE_GLOBAL_DAMAGE] == 2)
            {
                noalias(rStressVector) = prod(LinearElasticMatrix, rStrainVector);
            }
            else
            {
                this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
            }

            this->CalculateConstitutiveTensor(rConstitutiveMatrix, ReturnMappingVariables, LinearElasticMatrix);
        }
    }
    else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS) && Options.IsNot( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ))
    {
        // COMPUTE_STRESS
        Vector& rStressVector = rValues.GetStressVector();

        // This is for computing the linear response
        if ( CurrentProcessInfo[COMPUTE_GLOBAL_DAMAGE] == 2)
        {
            noalias(rStressVector) = prod(LinearElasticMatrix, rStrainVector);
        }
        else
        {
            this->CalculateReturnMapping(ReturnMappingVariables,AuxMatrix,rStressVector,LinearElasticMatrix,rStrainVector);
        }
    
    }

    if( Options.Is( ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE ) )
    {
        Vector& rStressVector = rValues.GetStressVector();
        
        this->UpdateInternalStateVariables(ReturnMappingVariables,rStressVector,LinearElasticMatrix,rStrainVector);
    }
}

} // Namespace Kratos

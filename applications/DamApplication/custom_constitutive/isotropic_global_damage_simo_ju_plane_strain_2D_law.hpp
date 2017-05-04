//
//   Project Name:        KratosDamApplication            $
//   Created by:          $Author:       Lorenzo Gracia   $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRAIN_2D_LAW_H_INCLUDED


// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/simo_ju_local_damage_plane_strain_2D_law.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/custom_flow_rules/local_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"
#include "dam_application_variables.h"


namespace Kratos
{

class IsotropicGlobalDamageSimoJuPlaneStrain2DLaw : public SimoJuLocalDamagePlaneStrain2DLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IsotropicGlobalDamageSimoJuPlaneStrain2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    IsotropicGlobalDamageSimoJuPlaneStrain2DLaw();
    
    /// Second Constructor
    IsotropicGlobalDamageSimoJuPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    IsotropicGlobalDamageSimoJuPlaneStrain2DLaw (const IsotropicGlobalDamageSimoJuPlaneStrain2DLaw& rOther);

    /// Destructor
    virtual ~IsotropicGlobalDamageSimoJuPlaneStrain2DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues);
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SimoJuLocalDamagePlaneStrain2DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SimoJuLocalDamagePlaneStrain2DLaw )
    }

}; // Class IsotropicGlobalDamageSimoJuPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRAIN_2D_LAW_H_INCLUDED defined

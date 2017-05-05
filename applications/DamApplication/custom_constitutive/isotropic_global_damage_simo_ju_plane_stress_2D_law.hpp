//
//   Project Name:        KratosDamApplication            $
//   Created by:          $Author:       Lorenzo Gracia   $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRESS_2D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRESS_2D_LAW_H_INCLUDED


// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/simo_ju_local_damage_plane_stress_2D_law.hpp"

#include "dam_application_variables.h"


namespace Kratos
{

class IsotropicGlobalDamageSimoJuPlaneStress2DLaw : public SimoJuLocalDamagePlaneStress2DLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IsotropicGlobalDamageSimoJuPlaneStress2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    IsotropicGlobalDamageSimoJuPlaneStress2DLaw();
    
    /// Second Constructor
    IsotropicGlobalDamageSimoJuPlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    IsotropicGlobalDamageSimoJuPlaneStress2DLaw (const IsotropicGlobalDamageSimoJuPlaneStress2DLaw& rOther);

    /// Destructor
    virtual ~IsotropicGlobalDamageSimoJuPlaneStress2DLaw();

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Serialization

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SimoJuLocalDamagePlaneStress2DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SimoJuLocalDamagePlaneStress2DLaw )
    }

}; // Class IsotropicGlobalDamageSimoJuPlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRESS_2D_LAW_H_INCLUDED defined

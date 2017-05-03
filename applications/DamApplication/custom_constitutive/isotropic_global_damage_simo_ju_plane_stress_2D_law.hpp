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

// External includes

// Project includes
#include "custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"


namespace Kratos
{


class KRATOS_API(SOLID_MECHANICS_APPLICATION) IsotropicGlobalDamageSimoJuPlaneStress2DLaw : public IsotropicDamageSimoJuPlaneStress2DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef FlowRule::Pointer                FlowRulePointer;
    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of IsotropicGlobalDamageSimoJuPlaneStress2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( IsotropicGlobalDamageSimoJuPlaneStress2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    IsotropicGlobalDamageSimoJuPlaneStress2DLaw();


    IsotropicGlobalDamageSimoJuPlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    IsotropicGlobalDamageSimoJuPlaneStress2DLaw (const IsotropicGlobalDamageSimoJuPlaneStress2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //IsotropicGlobalDamageSimoJuPlaneStress2DLaw& operator=(const IsotropicGlobalDamageSimoJuPlaneStress2DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~IsotropicGlobalDamageSimoJuPlaneStress2DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

       /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);


    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{
    ///@}


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsotropicDamageSimoJuPlaneStress2DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsotropicDamageSimoJuPlaneStress2DLaw )
    }



}; // Class IsotropicGlobalDamageSimoJuPlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_PLANE_STRESS_2D_LAW_H_INCLUDED defined

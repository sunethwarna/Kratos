//
//   Project Name:        KratosDamApplication            $
//   Created by:          $Author:       Lorenzo Gracia   $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                 May 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_3D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/isotropic_damage_simo_ju_3D_law.hpp"


namespace Kratos
{


class KRATOS_API(SOLID_MECHANICS_APPLICATION) IsotropicGlobalDamageSimoJu3DLaw : public IsotropicDamageSimoJu3DLaw
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
     * Counted pointer of IsotropicGlobalDamageSimoJu3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( IsotropicGlobalDamageSimoJu3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    IsotropicGlobalDamageSimoJu3DLaw();


    IsotropicGlobalDamageSimoJu3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    IsotropicGlobalDamageSimoJu3DLaw (const IsotropicGlobalDamageSimoJu3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //IsotropicGlobalDamageSimoJu3DLaw& operator=(const IsotropicGlobalDamageSimoJu3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~IsotropicGlobalDamageSimoJu3DLaw();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, IsotropicDamageSimoJu3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, IsotropicDamageSimoJu3DLaw )
    }



}; // Class IsotropicGlobalDamageSimoJu3DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_GLOBAL_DAMAGE_SIMO_JU_3D_LAW_H_INCLUDED defined

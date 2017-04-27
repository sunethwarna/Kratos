//
//   Project Name:        
//   Last modified by:    $Author:      
//   Date:                $Date:          
//   Revision:            $Revision:        
//

#if !defined(KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_GLOBAL_DAMAGE_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_GLOBAL_DAMAGE_ELEMENT_H_INCLUDED

/* Project includes */
#include "includes/serializer.h"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"

#include "custom_utilities/element_utilities.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

class SmallDisplacementThermoMechanicGlobalDamageElement : public SmallDisplacementThermoMechanicElement
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementThermoMechanicGlobalDamageElement );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallDisplacementThermoMechanicGlobalDamageElement();
    
    // Constructor 1
    SmallDisplacementThermoMechanicGlobalDamageElement(IndexType NewId, GeometryType::Pointer pGeometry);
    
    // Constructor 2
    SmallDisplacementThermoMechanicGlobalDamageElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);
    
    // Destructor
    virtual ~SmallDisplacementThermoMechanicGlobalDamageElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    void CalculateElementalSystem(LocalSystemComponents& rLocalSystem, ProcessInfo& rCurrentProcessInfo);

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementThermoMechanicElement )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementThermoMechanicElement )
    }
    
    
}; // Class SmallDisplacementThermoMechanicGlobalDamageElement

} // namespace Kratos

#endif // KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_GLOBAL_DAMAGE_ELEMENT_H_INCLUDED  defined 

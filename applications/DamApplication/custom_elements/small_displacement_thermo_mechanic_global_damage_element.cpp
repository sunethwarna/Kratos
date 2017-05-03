//
//   Project Name:        
//   Last modified by:    $Author:          
//   Date:                $Date:            
//   Revision:            $Revision:     
//

/* Project includes */
#include "custom_elements/small_displacement_thermo_mechanic_global_damage_element.hpp"


namespace Kratos
{

// Default Constructor
SmallDisplacementThermoMechanicGlobalDamageElement::SmallDisplacementThermoMechanicGlobalDamageElement() : SmallDisplacementThermoMechanicElement() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallDisplacementThermoMechanicGlobalDamageElement::SmallDisplacementThermoMechanicGlobalDamageElement( IndexType NewId, GeometryType::Pointer pGeometry ) : SmallDisplacementThermoMechanicElement( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallDisplacementThermoMechanicGlobalDamageElement::SmallDisplacementThermoMechanicGlobalDamageElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : SmallDisplacementThermoMechanicElement( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallDisplacementThermoMechanicGlobalDamageElement::~SmallDisplacementThermoMechanicGlobalDamageElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallDisplacementThermoMechanicGlobalDamageElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallDisplacementThermoMechanicGlobalDamageElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallDisplacementThermoMechanicGlobalDamageElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
							 ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Vector VolumeForce(dimension);
    noalias(VolumeForce) = ZeroVector(dimension);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values,PointNumber);

        //compute stresses and constitutive parameters
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

        //calculating weights for integration on the "reference configuration"
        double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


        //if ( dimension == 2 ) IntegrationWeight *= GetProperties()[THICKNESS];

        if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
            //contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );

        }

        if ( rLocalSystem.CalculationFlags.Is(SmallDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
            // For computing the global damage the contribution of volume forces are not needed but it is necessary for the normal calculus 
            if (rCurrentProcessInfo[COMPUTE_GLOBAL_DAMAGE] == 0)
            {
                //contribution to external forces
                VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );
            }

            this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, IntegrationWeight );

        }

    }

    KRATOS_CATCH( "" )
}


} // Namespace Kratos
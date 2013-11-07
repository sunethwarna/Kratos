//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYM_SPATIAL_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_SPATIAL_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/large_displacement_element.hpp"


namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

/// Axisym Spatial Lagrangian Element 2D geometries.

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D
 */

class AxisymSpatialLagrangianElement
    : public LargeDisplacementElement
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of AxisymSpatialLagrangianElement
    KRATOS_CLASS_POINTER_DEFINITION(AxisymSpatialLagrangianElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    AxisymSpatialLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymSpatialLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymSpatialLagrangianElement(AxisymSpatialLagrangianElement const& rOther);

    /// Destructor.
    virtual ~AxisymSpatialLagrangianElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymSpatialLagrangianElement& operator=(AxisymSpatialLagrangianElement const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

    // //************* GETTING METHODS

    //SET

    /**
     * Set a double  Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);


    //GET:

    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);


    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize();

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //int Check(const ProcessInfo& rCurrentProcessInfo);


    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    /**
     * Container for historical total elastic deformation measure
     */
    std::vector< Matrix > mDeformationGradientF0;

    /**
     * Container for the total deformation gradient determinants
     */
    Vector mDeterminantF0;

    ///@}
    ///@name Protected Operators
    ///@{
    AxisymSpatialLagrangianElement() : LargeDisplacementElement()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    GeneralVariables& rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight);

    /**
     * Calculation of the Total Mass of the Element
     */
    double& CalculateTotalMass(double& rTotalMass);


    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rK,
                                     GeneralVariables & rVariables,
                                     double& rIntegrationWeight
                                    );



    /**
     * Initialize Element General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    virtual void SetGeneralVariables(GeneralVariables& rVariables,
                                     ConstitutiveLaw::Parameters& rValues,
                                     const int & rPointNumber);

    /**
     * Finalize Element Internal Variables
     */
    virtual void FinalizeStepVariables(GeneralVariables & rVariables, const double& rPointNumber );


    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables,
                                     const double& rPointNumber);


    /**
     * Calculate Radius in the current and deformed geometry
     */
    void CalculateRadius(double & rCurrentRadius,
                         double & rReferenceRadius,
                         const Vector& rN);

    /**
     * Calculation of the Deformation Gradient F
     */
    void CalculateDeformationGradient(const Matrix& rDN_DX,
                                      Matrix& rF,
                                      Matrix& rDeltaPosition,
                                      double & rCurrentRadius,
                                      double & rReferenceRadius);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
                                            Matrix& rDN_DX,
                                            Vector& rN,
                                            double & rCurrentRadius);

    /**
     * Calculation of the Green Lagrange Strain Vector
     */
    void CalculateGreenLagrangeStrain(const Matrix& rF,
                                      Vector& rStrainVector);

    /**
     * Calculation of the Almansi Strain Vector
     */
    void CalculateAlmansiStrain(const Matrix& rF,
                                Vector& rStrainVector);



    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
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

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AxisymSpatialLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYM_SPATIAL_LAGRANGIAN_ELEMENT_H_INCLUDED  defined 

// ==============================================================================
//      ___   __   ____             
//     / _ | / /  / __/             
//    / __ |/ /__/ _/               
//   /_/ |_/____/___/  application  
//
// License:		     BSD License
//					 license: ALEApplication/license.txt
//==============================================================================

/* ****************************************************************************
 *  Projectname:         $KratosALEApplication
 *  Last Modified by:    $Author: A.Winterstein@tum.de $
 *  Date:                $Date: April 2017 $
 * ***************************************************************************/

#if !defined(KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED)
#define KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

/// This class implements a structural similarity mesh-updating scheme
/**
 * This mesh-updating scheme treats the mesh as a solid and therefore
 * solves the equations of solid mechanics using a corrotational formulation
 * and a linear elastic consitutive law. The stiffness of the elements
 * depends on their size and can be controlled by the Jacobian Determinant
 * weightened by an exponent. Therefore the current model part is copied and
 * thus it is not necessary to create the mesh moving element explicity in the
 * input file for the geometry. It is implemented for all types of elements.
 */

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

// template<unsigned int TDim>
class StructuralMeshMovingElement : public Element
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of StructuralMeshMovingElement
    KRATOS_CLASS_POINTER_DEFINITION(StructuralMeshMovingElement);

    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::DofsVectorType DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{
    StructuralMeshMovingElement(IndexType NewId, GeometryType::Pointer pGeometry);

    StructuralMeshMovingElement(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties);

    virtual ~StructuralMeshMovingElement()
    {
    }

    ///@}
    ///@name Operators
    ///@{
    /// Assignment operator.
    ///@}

    ///@name Operations
    ///@{
    /**
    * Returns the currently selected integration method
    * @return current integration method selected
    */
    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    BaseType::Pointer Create(IndexType NewId,
                             NodesArrayType const& rThisNodes,
                             PropertiesType::Pointer pProperties) const;

    BaseType::Pointer Create(IndexType NewId,
                             GeometryType::Pointer pGeom,
                             PropertiesType::Pointer pProperties) const;

    void CalculateLocalSystem(Matrix& rLeftHandSideMatrix,
                              Vector& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo);

    /**
    * Sets on rResult the ID's of the element degrees of freedom
    */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
    * Sets on rElementalDofList the degrees of freedom of the considered element
    * geometry
    */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    /// Initialize initial values of the element (must be called before
    /// calculation is done)
    void Initialize();

    void CalculateLocalSize(SizeType& LocalSize);

    Matrix SetAndModifyConstitutiveLaw(const int& rdimension, const double& rPointNumber);

    void CalculateBMatrix(Matrix& rB, const int& rdimension, const double& rPointNumber);

    void CheckElementMatrixDimension(Matrix& rLeftHandSideMatrix,
                                     Vector& rRightHandSideVector);

    void CalculateRightHandSide(Vector& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo);

    void CalculateAndAddRHS(Vector& rRightHandSideVector);

    void GetDisplacementValues(Vector& rValues, const int Step = 0);

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
    ///@}

    ///@name Protected Operators
    ///@{
    ///@}

    ///@name Protected Operations
    ///@{


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
    GeometryType::JacobiansType mInvJ0;
    Vector mDetJ0;
    double mTotalDomainInitialSize;
    ///@}
    ///@name Member Variables
    ///@{
    ///@}

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

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Un accessible methods
    ///@{
    ///@}

    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    StructuralMeshMovingElement();

    ///@}

}; // Class StructuralMeshMovingElement

///@}

///@name Type Definitions
///@{
///@}
}
// namespace Kratos.

#endif // KRATOS_STRUCTURAL_MESHMOVING_ELEMENT_INCLUDED

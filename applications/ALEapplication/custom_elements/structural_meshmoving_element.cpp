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

// System includes
// External includes
// Project includes
#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element.h"
#include "ale_application.h"

namespace Kratos
{

/******************************************************************************/
/******************************************************************************/
StructuralMeshMovingElement::StructuralMeshMovingElement( )
        : Element( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}


/******************************************************************************/
/******************************************************************************/
StructuralMeshMovingElement::StructuralMeshMovingElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

/******************************************************************************/
/******************************************************************************/
StructuralMeshMovingElement::StructuralMeshMovingElement(IndexType NewId,
                                                         GeometryType::Pointer pGeometry,
                                                         PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}


/******************************************************************************/
/******************************************************************************/
Element::Pointer StructuralMeshMovingElement::Create(IndexType NewId,
                                                     NodesArrayType const& rThisNodes,
                                                     PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();

    return BaseType::Pointer(new StructuralMeshMovingElement(
        NewId, rGeom.Create(rThisNodes), pProperties));
}

/******************************************************************************/
/******************************************************************************/
Element::Pointer StructuralMeshMovingElement::Create(IndexType NewId,
                                                     GeometryType::Pointer pGeom,
                                                     PropertiesType::Pointer pProperties) const
{
    return BaseType::Pointer(new StructuralMeshMovingElement(NewId, pGeom, pProperties));
}


/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::GetDisplacementValues(Vector& rValues, const int Step)
{
    KRATOS_TRY;

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType LocalSize;
    CalculateLocalSize(LocalSize);

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize, false);

    if (dimension == 2)
    {
        SizeType Index = 0;

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
        }
    }
    else if (dimension == 3)
    {
        SizeType Index = 0;

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_X, Step);
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_Y, Step);
            rValues[Index++] =
                rGeom[iNode].FastGetSolutionStepValue(MESH_DISPLACEMENT_Z, Step);
        }
    }

    KRATOS_CATCH("");
}

/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::Initialize()
{
    KRATOS_TRY;

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();

    GeometryType::JacobiansType J0;
    J0 = GetGeometry().Jacobian(J0, GetGeometry().GetDefaultIntegrationMethod());

    mInvJ0.resize(integration_points.size());
    mDetJ0.resize(integration_points.size(), false);
    mTotalDomainInitialSize = 0.00;

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        MathUtils<double>::InvertMatrix(
            J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber]);

        mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
    }

    KRATOS_CATCH("");
}

/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::CalculateLocalSize(SizeType& LocalSize)
{
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const int dimension = this->GetGeometry().WorkingSpaceDimension();

    LocalSize = NumNodes * dimension;

}

/******************************************************************************/
/******************************************************************************/
Matrix StructuralMeshMovingElement::SetAndModifyConstitutiveLaw(
    const int& rdimension, const double& rPointNumber)
{
    KRATOS_TRY;

    double YoungsModulus = 1.0f;
    const double PoissonCoefficient = 0.3f;

    Vector DetJ;
    this->GetGeometry().DeterminantOfJacobian(DetJ);

    // Stiffening of elements using Jacobian determinants and exponent between
    // 0.0 and 2.0
    const double J0 = 10; // Factor influences how far the displacement spreads
                           // into the fluid mesh
    const double Xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                           // elements; 0 = no stiffening
    double DetJMag = mTotalDomainInitialSize; 
    const double Quotient = J0 / DetJMag;
    YoungsModulus *= (DetJMag * pow(Quotient, Xi));

    double Lambda = (YoungsModulus * PoissonCoefficient) /
                    ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));
    double Mue = YoungsModulus / (2 * (1 - PoissonCoefficient));

    Matrix ConstitutiveMatrix;

    if (rdimension == 2)
    {
        ConstitutiveMatrix = ZeroMatrix(3, 3);

        ConstitutiveMatrix(0, 0) = Lambda + 2 * Mue;
        ConstitutiveMatrix(1, 1) = ConstitutiveMatrix(0, 0);
        ConstitutiveMatrix(2, 2) = Mue;
        ConstitutiveMatrix(0, 1) = Lambda;
        ConstitutiveMatrix(1, 0) = Lambda;
    }

    else if (rdimension == 3)
    {
        ConstitutiveMatrix = ZeroMatrix(6, 6);

        double ConstitutiveMatrixComponent01 =
            (YoungsModulus * (1 - PoissonCoefficient) /
             ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient)));

        double ConstitutiveMatrixComponent02 =
            YoungsModulus * (1 - PoissonCoefficient) / (1 + PoissonCoefficient);

        double ConstitutiveMatrixComponent03 =
            YoungsModulus * PoissonCoefficient /
            ((1 + PoissonCoefficient) * (1 - 2 * PoissonCoefficient));

        ConstitutiveMatrix(0, 0) = ConstitutiveMatrixComponent01;
        ConstitutiveMatrix(1, 1) = ConstitutiveMatrixComponent01;
        ConstitutiveMatrix(2, 2) = ConstitutiveMatrixComponent01;
        ConstitutiveMatrix(3, 3) = ConstitutiveMatrixComponent02;
        ConstitutiveMatrix(4, 4) = ConstitutiveMatrixComponent02;
        ConstitutiveMatrix(5, 5) = ConstitutiveMatrixComponent02;

        ConstitutiveMatrix(0, 1) = ConstitutiveMatrixComponent03;
        ConstitutiveMatrix(1, 0) = ConstitutiveMatrixComponent03;
        ConstitutiveMatrix(0, 2) = ConstitutiveMatrixComponent03;
        ConstitutiveMatrix(2, 0) = ConstitutiveMatrixComponent03;
        ConstitutiveMatrix(1, 2) = ConstitutiveMatrixComponent03;
        ConstitutiveMatrix(2, 1) = ConstitutiveMatrixComponent03;
    }

    return ConstitutiveMatrix;

    KRATOS_CATCH("");
}

/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::CalculateBMatrix(Matrix& rB,
    const int& rdimension, const double& rPointNumber)
{
    KRATOS_TRY;

    GeometryType::JacobiansType j;

    GetGeometry().Jacobian(j, GetGeometry().GetDefaultIntegrationMethod());

    Matrix F = prod(j[rPointNumber], mInvJ0[rPointNumber]);

    GeometryType::ShapeFunctionsGradientsType DN_De =
        this->GetGeometry().ShapeFunctionsLocalGradients();

    Matrix DN_DX = prod(DN_De[rPointNumber], mInvJ0[rPointNumber]);

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    if (rdimension == 2)
    {
        rB = ZeroMatrix(3, NumNodes * 2);

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = 2 * iNode;

            rB(0, Index + 0) = DN_DX(iNode, 0);
            rB(0, Index + 1) = 0.0;
            rB(1, Index + 0) = 0.0;
            rB(1, Index + 1) = DN_DX(iNode, 1);
            rB(2, Index + 0) = DN_DX(iNode, 1);
            rB(2, Index + 1) = DN_DX(iNode, 0);
        }
    }

    else if (rdimension == 3)
    {
        rB = ZeroMatrix(6, NumNodes * 3);

        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType index = 3 * iNode;

            rB(0, index + 0) = DN_DX(iNode, 0);
            rB(1, index + 1) = DN_DX(iNode, 1);
            rB(2, index + 2) = DN_DX(iNode, 2);

            rB(3, index + 0) = DN_DX(iNode, 1);
            rB(3, index + 1) = DN_DX(iNode, 0);

            rB(4, index + 1) = DN_DX(iNode, 2);
            rB(4, index + 2) = DN_DX(iNode, 1);

            rB(5, index + 0) = DN_DX(iNode, 2);
            rB(5, index + 2) = DN_DX(iNode, 0);
        }
    }

    KRATOS_CATCH("");
}


/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::CheckElementMatrixDimension(Matrix& rLeftHandSideMatrix,
                                                              Vector& rRightHandSideVector)
{
    KRATOS_TRY;

    SizeType LocalSize;
    CalculateLocalSize(LocalSize);

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    KRATOS_CATCH("");
}

/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::CalculateLocalSystem(Matrix& rLeftHandSideMatrix,
                                                       Vector& rRightHandSideVector,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const int dimension = this->GetGeometry().WorkingSpaceDimension();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();

    CheckElementMatrixDimension(rLeftHandSideMatrix, rRightHandSideVector);

    //Integration loop over element
    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        Matrix B;

        CalculateBMatrix(B, dimension, PointNumber);

        Matrix ConstitutiveMatrix = SetAndModifyConstitutiveLaw(dimension, PointNumber);

        // Compute LHS
        Matrix matrix_tmp;
        matrix_tmp = IntegrationWeight * prod(ConstitutiveMatrix, B);
        noalias(rLeftHandSideMatrix) += prod(trans(B), matrix_tmp);

        // Compute RHS
        SizeType LocalSize;
        CalculateLocalSize(LocalSize);
        Vector LastValues;
        this->GetDisplacementValues(LastValues, 0);
        noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, LastValues);
    }

    KRATOS_CATCH("");
}


/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::EquationIdVector(EquationIdVectorType& rResult,
                                                   ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.size();
    int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType LocalSize;
    CalculateLocalSize(LocalSize);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    int pos = this->GetGeometry()[0].GetDofPosition(MESH_DISPLACEMENT_X);
    if (dimension == 2)
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rResult[Index] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
            rResult[Index + 1] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
        }
    else
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rResult[Index] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_X,pos).EquationId();
            rResult[Index + 1] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Y,pos+1).EquationId();
            rResult[Index + 2] = rGeom[iNode].GetDof(MESH_DISPLACEMENT_Z,pos+2).EquationId();
        }

    KRATOS_CATCH("");
}


/******************************************************************************/
/******************************************************************************/
void StructuralMeshMovingElement::GetDofList(DofsVectorType& rElementalDofList,
                                             ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.size();
    int dimension = GetGeometry().WorkingSpaceDimension();
    SizeType LocalSize;
    CalculateLocalSize(LocalSize);

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    if (dimension == 2)
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rElementalDofList[Index] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_X);
            rElementalDofList[Index + 1] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Y);
        }
    else
        for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
        {
            SizeType Index = iNode * dimension;
            rElementalDofList[Index] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_X);
            rElementalDofList[Index + 1] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Y);
            rElementalDofList[Index + 2] = rGeom[iNode].pGetDof(MESH_DISPLACEMENT_Z);
        }

    KRATOS_CATCH("");
}


/******************************************************************************/
/******************************************************************************/
// Called in function "CalculateReactions" within the block builder and solver
void StructuralMeshMovingElement::CalculateRightHandSide(Vector& rRightHandSideVector,
                                                         ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const int dimension = this->GetGeometry().WorkingSpaceDimension();

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();
    SizeType LocalSize;
    CalculateLocalSize(LocalSize);

    Matrix LeftHandSideMatrix = ZeroMatrix(LocalSize, LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);
    rRightHandSideVector.clear();

    // Get nodal solutions
    GeometryType& rGeom = this->GetGeometry();
    Vector MeshDisplacement;
    MeshDisplacement.resize(LocalSize, false);

    for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
    {
        MeshDisplacement[dimension * iNode + 0] =
            rGeom[iNode].GetSolutionStepValue(MESH_DISPLACEMENT_X);
        MeshDisplacement[dimension * iNode + 1] =
            rGeom[iNode].GetSolutionStepValue(MESH_DISPLACEMENT_Y);
        MeshDisplacement[dimension * iNode + 2] =
            rGeom[iNode].GetSolutionStepValue(MESH_DISPLACEMENT_Z);
    }

    for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber)
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        Matrix B;
        CalculateBMatrix(B, dimension, PointNumber);

        Matrix ConstitutiveMatrix = SetAndModifyConstitutiveLaw(dimension, PointNumber);

        Vector StressVector =
            prod(ConstitutiveMatrix, Vector(prod(B, MeshDisplacement)));
        noalias(rRightHandSideVector) += IntegrationWeight * prod(trans(B), StressVector);
    }

    KRATOS_CATCH("");
}


/******************************************************************************/
/******************************************************************************/
// Only called in the CalculateLocalSystem Function above
void StructuralMeshMovingElement::CalculateAndAddRHS(Vector& rRightHandSideVector)
{
    KRATOS_TRY;

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const int dimension = this->GetGeometry().WorkingSpaceDimension();
        SizeType LocalSize;
    CalculateLocalSize(LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    GeometryType& rGeom = this->GetGeometry();
    for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
    {
        // Note that we need to divide by the neighbours since the final RHS
        // needs
        // to be the desired value
        //(the addition below is done as many times as elements we have where
        // this node appears)

        int number_neighbours = rGeom[iNode].GetValue(NEIGHBOUR_ELEMENTS).size();

        rRightHandSideVector[dimension * iNode + 0] +=
            rGeom[iNode].GetSolutionStepValue(MESH_RHS_X) / number_neighbours;
        rRightHandSideVector[dimension * iNode + 1] +=
            rGeom[iNode].GetSolutionStepValue(MESH_RHS_Y) / number_neighbours;
        rRightHandSideVector[dimension * iNode + 2] +=
            rGeom[iNode].GetSolutionStepValue(MESH_RHS_Z) / number_neighbours;
    }

    KRATOS_CATCH("");
}
} // Namespace Kratos

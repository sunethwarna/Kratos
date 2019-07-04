// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes


// External includes

// Project includes
#include "includes/checks.h"
#include "includes/define.h"
#include "custom_elements/laplacian_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
LaplacianElement::LaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
LaplacianElement::LaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer LaplacianElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer LaplacianElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<LaplacianElement>(NewId, pGeom, pProperties);
}

LaplacianElement::~LaplacianElement()
{
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS


    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients();
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues();

    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_points,dim);
    Matrix InvJ0(dim,dim);
    Vector temp(number_of_points);

    Vector heat_flux_local(number_of_points);
    Vector nodal_conductivity(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        heat_flux_local[node_element] = r_geometry[node_element].FastGetSolutionStepValue(HEAT_FLUX);
        nodal_conductivity[node_element] = r_geometry[node_element].FastGetSolutionStepValue(CONDUCTIVITY);
    }

    r_geometry.Jacobian(J0);
    double DetJ0;

    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point],InvJ0,DetJ0);

        //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point); //these are the N which correspond to the gauss point "i_point"
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double conductivity_gauss = inner_prod(N, nodal_conductivity);
        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * conductivity_gauss * prod(DN_DX, trans(DN_DX)); //

        // Calculating the local RHS
        const double qgauss = inner_prod(N, heat_flux_local);

        noalias(rRightHandSideVector) += IntToReferenceWeight*qgauss*N;
    }


    // RHS = ExtForces - K*temp;
    for (unsigned int i=0; i<number_of_points; i++)
        temp[i] = r_geometry[i].GetSolutionStepValue(TEMPERATURE) ; //this includes the - sign

    //axpy_prod(rLeftHandSideMatrix, temp, rRightHandSideVector, false);  //RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
    }
}

//************************************************************************************
//************************************************************************************
int LaplacianElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geom = GetGeometry();
    for (unsigned int i = 0; i < r_geom.PointsNumber(); i++)
    {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEAT_FLUX, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(CONDUCTIVITY, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TEMPERATURE, r_node);
    }

    return Element::Check(rCurrentProcessInfo);
}

} // Namespace Kratos



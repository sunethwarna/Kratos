//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_EVM_K_EPSILON_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_element_data.h"

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

namespace EvmKepsilonModelUtilities
{
/// Node type
using NodeType = Node<3>;
using GeometryType = Geometry<NodeType>;


double CalculateTurbulentViscosity(const double C_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double f_mu);

double CalculateFmu(const double y_plus);

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate);

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy);

double CalculateGamma(const double C_mu,
                      const double f_mu,
                      const double turbulent_kinetic_energy,
                      const double turbulent_kinematic_viscosity);

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos

#endif
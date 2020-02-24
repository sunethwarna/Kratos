//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/variables.h"
#include "includes/variables_derivatives.h"


namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(VariableTimeDerivativeHasMethods, KratosCoreFastSuite)
{
    const bool checktrue = VariablesDerivatives<Variable<array_1d<double, 3>>>::HasTimeDerivative(DISPLACEMENT);
    KRATOS_CHECK(checktrue);
    const bool checkfalse = VariablesDerivatives<Variable<array_1d<double, 3>>>::HasTimeDerivative(VECTOR_LAGRANGE_MULTIPLIER);
    KRATOS_CHECK_IS_FALSE(checkfalse);
}

KRATOS_TEST_CASE_IN_SUITE(VariableTimeDerivativeGetMethods, KratosCoreFastSuite)
{
    const auto& r_velocity = VariablesDerivatives<Variable<array_1d<double, 3>>>::GetFirstTimeDerivative(DISPLACEMENT);
    KRATOS_CHECK(r_velocity.Name() == "VELOCITY");
    const auto& r_acceleration = VariablesDerivatives<Variable<array_1d<double, 3>>>::GetSecondTimeDerivative(DISPLACEMENT);
    KRATOS_CHECK(r_acceleration.Name() == "ACCELERATION");
}

KRATOS_TEST_CASE_IN_SUITE(VariableResidualHasMethods, KratosCoreFastSuite)
{
    const bool checktrue = VariablesDerivatives<Variable<array_1d<double, 3>>>::HasResidualVariable(DISPLACEMENT);
    KRATOS_CHECK(checktrue);
    const bool checkfalse = VariablesDerivatives<Variable<array_1d<double, 3>>>::HasResidualVariable(VECTOR_LAGRANGE_MULTIPLIER);
    KRATOS_CHECK_IS_FALSE(checkfalse);
}

KRATOS_TEST_CASE_IN_SUITE(VariableResidualGetMethods, KratosCoreFastSuite)
{
    const auto& r_force_residual = VariablesDerivatives<Variable<array_1d<double, 3>>>::GetResidualVariable(DISPLACEMENT);
    KRATOS_CHECK(r_force_residual.Name() == "FORCE_RESIDUAL");
    const auto& r_moment_residual = VariablesDerivatives<Variable<array_1d<double, 3>>>::GetResidualVariable(ROTATION);
    KRATOS_CHECK(r_moment_residual.Name() == "MOMENT_RESIDUAL");
}

}
}  // namespace Kratos.

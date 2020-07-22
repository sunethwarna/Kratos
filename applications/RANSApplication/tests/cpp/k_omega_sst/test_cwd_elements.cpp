//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"

// Application includes
#include "../stabilization_method_test_utilities.h"
#include "custom_utilities/test_utilities.h"
#include "rans_application_variables.h"
#include "test_utilities.h"

namespace Kratos
{
namespace Testing
{
namespace
{
ModelPart& RansKOmegaSSTKCWD2D3N_SetUp(
    Model& rModel)
{
    ModelPart& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTK2D3N_SetUp(
        rModel, "RansKOmegaSSTKCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    return r_model_part;
}

ModelPart& RansKOmegaSSTOmegaCWD2D3N_SetUp(
    Model& rModel)
{
    ModelPart& r_model_part = KOmegaSSTTestUtilities::RansKOmegaSSTOmega2D3N_SetUp(
        rModel, "RansKOmegaSSTOmegaCWD2D3N");

    StabilizationMethodTestUtilities::InitializeCrossWindStabilizationConstants(
        r_model_part.GetProcessInfo());

    return r_model_part;
}

} // namespace

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_KINETIC_ENERGY);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 5.37345251410611801646e+00;
    ref_RHS[1] = 4.41760934465202481647e+00;
    ref_RHS[2] = 3.70832202716588232860e+00;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = 5.37345251410611801646e+00;
    ref_RHS[1] = 4.41760934465202481647e+00;
    ref_RHS[2] = 3.70832202716588232860e+00;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 4.28929893595158500830e+03;
    ref_LHS(0, 1) = -3.88524890504886707276e+03;
    ref_LHS(0, 2) = 1.62018598368927030151e+02;
    ref_LHS(1, 0) = -3.88408619762882608484e+03;
    ref_LHS(1, 1) = 8.46659492799086729065e+03;
    ref_LHS(1, 2) = -3.85622796801704043901e+03;
    ref_LHS(2, 0) = 1.62289207973309089539e+02;
    ref_LHS(2, 1) = -3.85636322449375893484e+03;
    ref_LHS(2, 2) = 4.44109764046969576157e+03;

    ref_RHS[0] = 5.69857108792817816720e+04;
    ref_RHS[1] = -2.70857677327735931613e+05;
    ref_RHS[2] = 1.10122239470008848002e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.50532337313251862732e-01;
    ref_M(0, 1) = 4.19634783864842250689e-02;
    ref_M(0, 2) = 4.17743588711753238707e-02;
    ref_M(1, 0) = 4.12513970335363494568e-02;
    ref_M(1, 1) = 2.49710814830827881883e-01;
    ref_M(1, 2) = 4.15718439184399671249e-02;
    ref_M(2, 0) = 4.15488302453199884190e-02;
    ref_M(2, 1) = 4.16586587559343163312e-02;
    ref_M(2, 2) = 2.49986887089775289272e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTKCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTKCWD2D3N_SetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 4.28929893595158500830e+03;
    ref_D(0, 1) = -3.88524890504886707276e+03;
    ref_D(0, 2) = 1.62018598368927030151e+02;
    ref_D(1, 0) = -3.88408619762882608484e+03;
    ref_D(1, 1) = 8.46659492799086729065e+03;
    ref_D(1, 2) = -3.85622796801704043901e+03;
    ref_D(2, 0) = 1.62289207973309089539e+02;
    ref_D(2, 1) = -3.85636322449375893484e+03;
    ref_D(2, 2) = 4.44109764046969576157e+03;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_EquationIdVector, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestEquationIdVector<ModelPart::ElementsContainerType>(r_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_GetDofList, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    RansApplicationTestUtilities::TestGetDofList<ModelPart::ElementsContainerType>(
        r_model_part, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateLocalSystem, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS;
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalSystem(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -6.35234246208723379823e+03;
    ref_RHS[1] = -6.32799836853508077184e+03;
    ref_RHS[2] = -6.33783278963220982405e+03;
    ref_LHS = ZeroMatrix(3, 3);

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateRightHandSide, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    Vector RHS, ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateRightHandSide(
        RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_RHS[0] = -6.35234246208723379823e+03;
    ref_RHS[1] = -6.32799836853508077184e+03;
    ref_RHS[2] = -6.33783278963220982405e+03;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateLocalVelocityContribution, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    Matrix LHS, ref_LHS(3, 3);
    Vector RHS(3, 0.0), ref_RHS(3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateLocalVelocityContribution(
        LHS, RHS, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_LHS(0, 0) = 1.07790365761357202246e+03;
    ref_LHS(0, 1) = -6.93686918600685430647e+02;
    ref_LHS(0, 2) = 1.01110635938086389274e+02;
    ref_LHS(1, 0) = -6.92524211180644670094e+02;
    ref_LHS(1, 1) = 2.03122752185442345763e+03;
    ref_LHS(1, 2) = -7.12934423988043818099e+02;
    ref_LHS(2, 0) = 1.01381245542468462872e+02;
    ref_LHS(2, 1) = -7.13069680464762313932e+02;
    ref_LHS(2, 2) = 9.81813831919713720708e+02;

    ref_RHS[0] = 1.33866736437548097456e+05;
    ref_RHS[1] = -1.29646252541568363085e+06;
    ref_RHS[2] = 4.70260314797010854818e+05;

    KRATOS_CHECK_VECTOR_NEAR(RHS, ref_RHS, 1e-12);
    KRATOS_CHECK_MATRIX_NEAR(LHS, ref_LHS, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateMassMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    Matrix M, ref_M(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateMassMatrix(
        M, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_M(0, 0) = 2.50470129386400064408e-01;
    ref_M(0, 1) = 4.19235197346189275569e-02;
    ref_M(0, 2) = 4.16201571536654948980e-02;
    ref_M(1, 0) = 4.12993034971783079534e-02;
    ref_M(1, 1) = 2.49741025773473374061e-01;
    ref_M(1, 2) = 4.16934349516055735574e-02;
    ref_M(2, 0) = 4.15631627117620111589e-02;
    ref_M(2, 1) = 4.16683781937928071626e-02;
    ref_M(2, 2) = 2.50019365400278581468e-01;

    KRATOS_CHECK_MATRIX_NEAR(M, ref_M, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(RansKOmegaSSTOmegaCWD2D3N_CalculateDampingMatrix, KratosRansFastSuite)
{
    // Setup:
    Model model;
    auto& r_model_part = RansKOmegaSSTOmegaCWD2D3N_SetUp(model);

    // Test:
    Matrix D, ref_D(3, 3);
    auto& r_element = r_model_part.Elements().front();
    r_element.CalculateDampingMatrix(
        D, static_cast<const ProcessInfo&>(r_model_part.GetProcessInfo()));

    // setting reference values
    ref_D(0, 0) = 1.07790365761357202246e+03;
    ref_D(0, 1) = -6.93686918600685430647e+02;
    ref_D(0, 2) = 1.01110635938086389274e+02;
    ref_D(1, 0) = -6.92524211180644670094e+02;
    ref_D(1, 1) = 2.03122752185442345763e+03;
    ref_D(1, 2) = -7.12934423988043818099e+02;
    ref_D(2, 0) = 1.01381245542468462872e+02;
    ref_D(2, 1) = -7.13069680464762313932e+02;
    ref_D(2, 2) = 9.81813831919713720708e+02;

    KRATOS_CHECK_MATRIX_NEAR(D, ref_D, 1e-12);
}

} // namespace Testing
} // namespace Kratos.

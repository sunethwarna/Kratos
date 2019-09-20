// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
typedef array_1d<double, 3> Vector3;
typedef std::size_t IndexType;

// General pourpose
KRATOS_CREATE_VARIABLE(int, INTEGRATION_ORDER); // The integration order considered on the element
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_MATERIAL_AXIS_1)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_MATERIAL_AXIS_2)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_MATERIAL_AXIS_3)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CENTER_OF_GRAVITY)
KRATOS_CREATE_VARIABLE(double, MASS_MOMENT_OF_INERTIA)
KRATOS_CREATE_VARIABLE(Matrix, ELASTICITY_TENSOR)

// Generalized eigenvalue problem
KRATOS_CREATE_VARIABLE(int, BUILD_LEVEL)
KRATOS_CREATE_VARIABLE(Vector, EIGENVALUE_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, EIGENVECTOR_MATRIX)

KRATOS_CREATE_VARIABLE(Matrix, MODAL_MASS_MATRIX)
KRATOS_CREATE_VARIABLE(Matrix, MODAL_STIFFNESS_MATRIX)

// Geometrical
KRATOS_CREATE_VARIABLE(double, AREA)
KRATOS_CREATE_VARIABLE(double, IT)
KRATOS_CREATE_VARIABLE(double, IY)
KRATOS_CREATE_VARIABLE(double, IZ)
KRATOS_CREATE_VARIABLE(double, CROSS_AREA)
KRATOS_CREATE_VARIABLE(double, MEAN_RADIUS)
KRATOS_CREATE_VARIABLE(int, SECTION_SIDES)
KRATOS_CREATE_VARIABLE(Matrix, GEOMETRIC_STIFFNESS)
KRATOS_CREATE_VARIABLE(Matrix, LOCAL_ELEMENT_ORIENTATION)
KRATOS_CREATE_VARIABLE(double, MATERIAL_ORIENTATION_ANGLE)
KRATOS_CREATE_VARIABLE(bool, USE_CONSISTENT_MASS_MATRIX)
KRATOS_CREATE_VARIABLE(Vector, CONDENSED_DOF_LIST)

// Truss generalized variables
KRATOS_CREATE_VARIABLE(double, TRUSS_PRESTRESS_PK2)
KRATOS_CREATE_VARIABLE(double, HARDENING_MODULUS_1D)
KRATOS_CREATE_VARIABLE(double, TANGENT_MODULUS)
KRATOS_CREATE_VARIABLE(double, PLASTIC_ALPHA)

// Beam generalized variables
KRATOS_CREATE_VARIABLE(double, AREA_EFFECTIVE_Y)
KRATOS_CREATE_VARIABLE(double, AREA_EFFECTIVE_Z)
KRATOS_CREATE_VARIABLE(double, INERTIA_ROT_Y)
KRATOS_CREATE_VARIABLE(double, INERTIA_ROT_Z)
KRATOS_CREATE_VARIABLE(Vector, LOCAL_AXES_VECTOR)
KRATOS_CREATE_VARIABLE(double, TORSIONAL_INERTIA)
KRATOS_CREATE_VARIABLE(double, I22)
KRATOS_CREATE_VARIABLE(double, I33)
KRATOS_CREATE_VARIABLE(double, LUMPED_MASS_ROTATION_COEFFICIENT)

// Shell generalized variables
KRATOS_CREATE_VARIABLE(bool, STENBERG_SHEAR_STABILIZATION_SUITABLE)
KRATOS_CREATE_VARIABLE(double, SHELL_OFFSET)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRAIN)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRAIN_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_CURVATURE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_CURVATURE_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_FORCE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_FORCE_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_MOMENT)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_MOMENT_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRESS_TOP_SURFACE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRESS_TOP_SURFACE_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRESS_MIDDLE_SURFACE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRESS_MIDDLE_SURFACE_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRESS_BOTTOM_SURFACE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_STRESS_BOTTOM_SURFACE_GLOBAL)
KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS_TOP_SURFACE)
KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS_MIDDLE_SURFACE)
KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS_BOTTOM_SURFACE)
KRATOS_CREATE_VARIABLE(double, SHEAR_ANGLE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_4PLY_THROUGH_THICKNESS)
KRATOS_CREATE_VARIABLE(double, TSAI_WU_RESERVE_FACTOR)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_LAMINA_STRENGTHS)

// Shell energies
KRATOS_CREATE_VARIABLE(double, SHELL_ELEMENT_MEMBRANE_ENERGY)
KRATOS_CREATE_VARIABLE(double, SHELL_ELEMENT_BENDING_ENERGY)
KRATOS_CREATE_VARIABLE(double, SHELL_ELEMENT_SHEAR_ENERGY)
KRATOS_CREATE_VARIABLE(double, SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION)
KRATOS_CREATE_VARIABLE(double, SHELL_ELEMENT_BENDING_ENERGY_FRACTION)
KRATOS_CREATE_VARIABLE(double, SHELL_ELEMENT_SHEAR_ENERGY_FRACTION)

// Membrane variables
KRATOS_CREATE_VARIABLE(Matrix, MEMBRANE_PRESTRESS)
KRATOS_CREATE_VARIABLE(Vector, PRESTRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, PRESTRESS_AXIS_1_GLOBAL)
KRATOS_CREATE_VARIABLE(Vector, PRESTRESS_AXIS_2_GLOBAL)
KRATOS_CREATE_VARIABLE(Matrix, PRESTRESS_AXIS_1)
KRATOS_CREATE_VARIABLE(Matrix, PRESTRESS_AXIS_2)
KRATOS_CREATE_VARIABLE(std::string, PROJECTION_TYPE_COMBO)


// Formfinding
KRATOS_CREATE_VARIABLE(double, LAMBDA_MAX)
KRATOS_CREATE_VARIABLE(bool, IS_FORMFINDING)
KRATOS_CREATE_VARIABLE(Matrix, BASE_REF_1)
KRATOS_CREATE_VARIABLE(Matrix, BASE_REF_2)

// Cross section
KRATOS_CREATE_VARIABLE(ShellCrossSection::Pointer, SHELL_CROSS_SECTION)
KRATOS_CREATE_VARIABLE(int, SHELL_CROSS_SECTION_OUTPUT_PLY_ID)
KRATOS_CREATE_VARIABLE(double, SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION)
KRATOS_CREATE_VARIABLE(Matrix, SHELL_ORTHOTROPIC_LAYERS)

// Nodal stiffness for the nodal concentrated element
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_INITIAL_DISPLACEMENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_DISPLACEMENT_STIFFNESS)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_INITIAL_ROTATION)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_STIFFNESS)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_DAMPING_RATIO)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_DAMPING_RATIO)

// For explicit central difference scheme
KRATOS_CREATE_VARIABLE(double, MASS_FACTOR)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MIDDLE_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACTIONAL_ACCELERATION)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACTIONAL_ANGULAR_ACCELERATION)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MIDDLE_ANGULAR_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_INERTIA)
KRATOS_CREATE_VARIABLE(double, NODAL_DISPLACEMENT_DAMPING)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATION_DAMPING)

// CONDITIONS
/* Beam conditions */
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(POINT_MOMENT)

// Adding the SPRISM EAS variables
KRATOS_CREATE_VARIABLE(double, ALPHA_EAS);
KRATOS_CREATE_VARIABLE(bool, CONSIDER_IMPLICIT_EAS_SPRISM_ELEMENT);
KRATOS_CREATE_VARIABLE(bool, CONSIDER_TOTAL_LAGRANGIAN_SPRISM_ELEMENT);
KRATOS_CREATE_VARIABLE(bool, PURE_EXPLICIT_RHS_COMPUTATION);

// Reset equations ids "flag"
KRATOS_CREATE_VARIABLE(bool, RESET_EQUATION_IDS);

// Adding the SPRISM additional variables
KRATOS_CREATE_VARIABLE(double, ANG_ROT);

// Adding the SPRISM variable to deactivate the quadratic interpolation
KRATOS_CREATE_VARIABLE(bool, CONSIDER_QUADRATIC_SPRISM_ELEMENT);

// Additional strain measures
KRATOS_CREATE_VARIABLE(Vector, HENCKY_STRAIN_VECTOR);
KRATOS_CREATE_VARIABLE(Matrix, HENCKY_STRAIN_TENSOR);

KRATOS_CREATE_VARIABLE(double, VON_MISES_STRESS)

KRATOS_CREATE_VARIABLE(Matrix, REFERENCE_DEFORMATION_GRADIENT);
KRATOS_CREATE_VARIABLE(double, REFERENCE_DEFORMATION_GRADIENT_DETERMINANT);

// Rayleigh variables
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_ALPHA)
KRATOS_CREATE_VARIABLE(double, RAYLEIGH_BETA)

// System damping
KRATOS_CREATE_VARIABLE(double,  SYSTEM_DAMPING_RATIO )
KRATOS_CREATE_VARIABLE(double,  SECOND_SYSTEM_DAMPING_RATIO )

// Nodal load variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

// Condition load variables
KRATOS_CREATE_VARIABLE(Vector, POINT_LOADS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, LINE_LOADS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, SURFACE_LOADS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, POSITIVE_FACE_PRESSURES_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, NEGATIVE_FACE_PRESSURES_VECTOR)

// Response function variables
KRATOS_CREATE_VARIABLE(double, RESPONSE_VALUE)

// Constitutive laws variables
KRATOS_CREATE_VARIABLE(bool, INELASTIC_FLAG)
KRATOS_CREATE_VARIABLE(double, INFINITY_YIELD_STRESS)
KRATOS_CREATE_VARIABLE(double, YIELD_STRESS_TENSION)
KRATOS_CREATE_VARIABLE(Vector, PLASTIC_STRAIN_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, BACK_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Matrix, PLASTIC_DEFORMATION_GRADIENT)
KRATOS_CREATE_VARIABLE(double, YIELD_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, DILATANCY_ANGLE)
KRATOS_CREATE_VARIABLE(int, SOFTENING_TYPE)
KRATOS_CREATE_VARIABLE(int, SOFTENING_TYPE_COMPRESSION)
KRATOS_CREATE_VARIABLE(int, HARDENING_CURVE)
KRATOS_CREATE_VARIABLE(int, MAX_NUMBER_NL_CL_ITERATIONS)
KRATOS_CREATE_VARIABLE(double, VISCOUS_PARAMETER)
KRATOS_CREATE_VARIABLE(double, DELAY_TIME)
KRATOS_CREATE_VARIABLE(double, MAXIMUM_STRESS)
KRATOS_CREATE_VARIABLE(double, MAXIMUM_STRESS_POSITION)
KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS)
KRATOS_CREATE_VARIABLE(double, FRICTION_ANGLE)
KRATOS_CREATE_VARIABLE(double, COHESION)
KRATOS_CREATE_VARIABLE(double, DAMAGE)
KRATOS_CREATE_VARIABLE(double, THRESHOLD)
KRATOS_CREATE_VARIABLE(Matrix, INTEGRATED_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Matrix, PLASTIC_STRAIN_TENSOR)
KRATOS_CREATE_VARIABLE(Matrix, BACK_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Vector, CURVE_FITTING_PARAMETERS)
KRATOS_CREATE_VARIABLE(Vector, PLASTIC_STRAIN_INDICATORS)
KRATOS_CREATE_VARIABLE(double, EQUIVALENT_PLASTIC_STRAIN)
KRATOS_CREATE_VARIABLE(Vector, KINEMATIC_PLASTICITY_PARAMETERS)
KRATOS_CREATE_VARIABLE(int, KINEMATIC_HARDENING_TYPE)
KRATOS_CREATE_VARIABLE(bool, CONSIDER_PERTURBATION_THRESHOLD)
KRATOS_CREATE_VARIABLE(int, TANGENT_OPERATOR_ESTIMATION)
KRATOS_CREATE_VARIABLE(Matrix, TENSION_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Matrix, COMPRESSION_STRESS_TENSOR)
KRATOS_CREATE_VARIABLE(Vector, TENSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, COMPRESSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, EFFECTIVE_TENSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(Vector, EFFECTIVE_COMPRESSION_STRESS_VECTOR)
KRATOS_CREATE_VARIABLE(double, EXPONENTIAL_SATURATION_YIELD_STRESS)
KRATOS_CREATE_VARIABLE(double, ACCUMULATED_PLASTIC_STRAIN)
KRATOS_CREATE_VARIABLE(Vector, HIGH_CYCLE_FATIGUE_COEFFICIENTS)
KRATOS_CREATE_VARIABLE(double, FATIGUE_REDUCTION_FACTOR)
KRATOS_CREATE_VARIABLE(int, NUMBER_OF_CYCLES)
KRATOS_CREATE_VARIABLE(double, WOHLER_STRESS)
KRATOS_CREATE_VARIABLE(Vector, HARDENING_MODULI_VECTOR)

// D+D- Damage Constitutive laws variables
KRATOS_CREATE_VARIABLE(double, DAMAGE_TENSION)
KRATOS_CREATE_VARIABLE(double, DAMAGE_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, THRESHOLD_TENSION)
KRATOS_CREATE_VARIABLE(double, THRESHOLD_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS_TENSION)
KRATOS_CREATE_VARIABLE(double, UNIAXIAL_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_DAMAGE_PROCESS)

// D+D- Damage Constitutive laws variables, additional Masonry 2D & 3D
KRATOS_CREATE_VARIABLE(double, DAMAGE_ONSET_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, BIAXIAL_COMPRESSION_MULTIPLIER)
KRATOS_CREATE_VARIABLE(double, FRACTURE_ENERGY_TENSION)
KRATOS_CREATE_VARIABLE(double, RESIDUAL_STRESS_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_C1)
KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_C2)
KRATOS_CREATE_VARIABLE(double, BEZIER_CONTROLLER_C3)
KRATOS_CREATE_VARIABLE(double, YIELD_STRAIN_COMPRESSION)
KRATOS_CREATE_VARIABLE(double, SHEAR_COMPRESSION_REDUCTOR)
KRATOS_CREATE_VARIABLE(double, TRIAXIAL_COMPRESSION_COEFFICIENT)

// Adjoint Variables
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_DISPLACEMENT)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_ROTATION)
KRATOS_CREATE_VARIABLE(double, PERTURBATION_SIZE)
KRATOS_CREATE_VARIABLE(bool, ADAPT_PERTURBATION_SIZE)

// Variables for output of sensitivities
KRATOS_CREATE_VARIABLE( double, CROSS_AREA_SENSITIVITY );
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD_SENSITIVITY );
KRATOS_CREATE_VARIABLE(double, I22_SENSITIVITY );
KRATOS_CREATE_VARIABLE(double, I33_SENSITIVITY );
KRATOS_CREATE_VARIABLE(double, THICKNESS_SENSITIVITY );
KRATOS_CREATE_VARIABLE(double, YOUNG_MODULUS_SENSITIVITY );
KRATOS_CREATE_VARIABLE(double, AREA_EFFECTIVE_Y_SENSITIVITY );
KRATOS_CREATE_VARIABLE(double, AREA_EFFECTIVE_Z_SENSITIVITY );
KRATOS_CREATE_VARIABLE(bool, IS_ADJOINT );

// Variables to for computing parts of sensitivity analysis
KRATOS_CREATE_VARIABLE( int, TRACED_STRESS_TYPE );
KRATOS_CREATE_VARIABLE( Matrix, STRESS_DISP_DERIV_ON_NODE);
KRATOS_CREATE_VARIABLE( Matrix, STRESS_DISP_DERIV_ON_GP );
KRATOS_CREATE_VARIABLE( Matrix, STRESS_DESIGN_DERIVATIVE_ON_NODE);
KRATOS_CREATE_VARIABLE( Matrix, STRESS_DESIGN_DERIVATIVE_ON_GP );
KRATOS_CREATE_VARIABLE( Vector, STRESS_ON_GP );
KRATOS_CREATE_VARIABLE( Vector, STRESS_ON_NODE );
KRATOS_CREATE_VARIABLE( std::string, DESIGN_VARIABLE_NAME );
}

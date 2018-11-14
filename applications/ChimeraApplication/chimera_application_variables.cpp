#include "chimera_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CHIM_NEUMANN_COND);

KRATOS_CREATE_VARIABLE(MpcDataPointerVectorType, MPC_DATA_CONTAINER)

//KRATOS_CREATE_VARIABLE(bool, IS_WEAK);

KRATOS_CREATE_VARIABLE(double, BOUNDARY_NODE);
KRATOS_CREATE_VARIABLE(double, FLUX);
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(TRACTION);
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHEAR_FORCE);
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE);


KRATOS_CREATE_VARIABLE(int,PATCH_INDEX)
KRATOS_CREATE_VARIABLE(double,TAUONE)
KRATOS_CREATE_VARIABLE(double,TAUTWO)
KRATOS_CREATE_VARIABLE(double,PRESSURE_MASSMATRIX_COEFFICIENT)

//KRATOS_CREATE_VARIABLE(double,Y_WALL)
KRATOS_CREATE_VARIABLE(double,SUBSCALE_PRESSURE)
KRATOS_CREATE_VARIABLE(double, C_DES)
//    KRATOS_CREATE_VARIABLE(double, C_SMAGORINSKY)
KRATOS_CREATE_VARIABLE(double, CHARACTERISTIC_VELOCITY)


KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY)

// Non-Newtonian constitutive relations
KRATOS_CREATE_VARIABLE(double, REGULARIZATION_COEFFICIENT)

KRATOS_CREATE_VARIABLE(double, BINGHAM_SMOOTHER)
KRATOS_CREATE_VARIABLE(double, GEL_STRENGTH )

// Q-Criterion (for vortex visualization)
KRATOS_CREATE_VARIABLE(double, Q_VALUE)

// Vorticity
KRATOS_CREATE_VARIABLE(double, VORTICITY_MAGNITUDE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT)

// For swimming DEM
KRATOS_CREATE_VARIABLE(Vector, NODAL_WEIGHTS)

// Embedded fluid variables
KRATOS_CREATE_VARIABLE(int, EMBEDDED_IS_ACTIVE)
KRATOS_CREATE_VARIABLE(double, EMBEDDED_WET_PRESSURE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(EMBEDDED_WET_VELOCITY)

}

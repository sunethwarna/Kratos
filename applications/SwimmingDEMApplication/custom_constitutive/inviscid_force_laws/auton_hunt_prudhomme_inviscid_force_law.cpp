// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "auton_hunt_prudhomme_inviscid_force_law.h"

namespace Kratos {

    InviscidForceLaw::Pointer AutonHuntPrudhommeInviscidForceLaw::Clone() const {
        AutonHuntPrudhommeInviscidForceLaw::Pointer p_clone(new AutonHuntPrudhommeInviscidForceLaw(*this));
        return p_clone;
    }

    AutonHuntPrudhommeInviscidForceLaw::AutonHuntPrudhommeInviscidForceLaw(Parameters& r_parameters)
        :mDoApplyFaxenCorrections(false)
    {

    }

    void AutonHuntPrudhommeInviscidForceLaw::Initialize(const ProcessInfo& r_process_info){}

    std::string AutonHuntPrudhommeInviscidForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Standard inviscid force law";
        return type_of_law;
    }

    void AutonHuntPrudhommeInviscidForceLaw::ComputeForce(Geometry<Node<3> >& r_geometry,
                                                const double fluid_density,
                                                const double displaced_volume,
                                                array_1d<double, 3>& virtual_mass_plus_undisturbed_flow_force,
                                                const ProcessInfo& r_current_process_info)
    {
        const array_1d<double, 3>& fluid_acc = r_geometry[0].FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
        const double radius = r_geometry[0].FastGetSolutionStepValue(RADIUS);
        const double virtual_mass_coeff = 0.5;
        const double fluid_mass = displaced_volume * fluid_density;
        mLastVirtualMassAddedMass = virtual_mass_coeff * fluid_mass;

        array_1d<double, 3> slip_acc = fluid_acc; // the particle acceleration is assumed to be treated implicitly through the added_mass

        if (mDoApplyFaxenCorrections) {
            const array_1d<double, 3>& fluid_vel_laplacian_rate = r_geometry[0].FastGetSolutionStepValue(FLUID_VEL_LAPL_RATE_PROJECTED);
            noalias(slip_acc) -= 0.1 * radius * radius * fluid_vel_laplacian_rate; // add Faxen term
        }

        noalias(virtual_mass_plus_undisturbed_flow_force) = fluid_mass * (virtual_mass_coeff * slip_acc + fluid_acc);
    }

} // namespace Kratos

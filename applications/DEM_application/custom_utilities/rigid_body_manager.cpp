// Created by: $Author: Salva Latorre, latorre@cimne.upc.edu

// System includes
#include <iostream>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "rigid_body_manager.h"

namespace Kratos
{
    /// Turn back information as a string.
    std::string RigidBodyManager::Info() const { return "";}

    /// Print information about this object.
    void RigidBodyManager::PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    void RigidBodyManager::PrintData(std::ostream& rOStream) const {}

} // namespace Kratos

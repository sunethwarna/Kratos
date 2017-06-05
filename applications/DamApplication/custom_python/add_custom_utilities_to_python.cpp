//   
//   Project Name:                  KratosDamApplication $
//   Last Modified by:    $Author:        Lorenzo Gracia $
//   Date:                $Date:           November 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/streamlines_output_3D_utilities.hpp"
#include "custom_utilities/global_joint_stress_utility.hpp"

namespace Kratos
{
	
namespace Python
{

void  AddCustomUtilitiesToPython() 
{
    using namespace boost::python;
    
    class_< StreamlinesOutput3DUtilities > ("StreamlinesOutput3DUtilities", init<>())
    .def("ComputeOutputStep",&StreamlinesOutput3DUtilities::ComputeOutputStep)
    ;

    class_< GlobalJointStressUtility > ("GlobalJointStressUtility", init<>())
    .def("ComputingGlobalStress",&GlobalJointStressUtility::ComputingGlobalStress)
    ;
  
}

}  // namespace Python.
} // Namespace Kratos

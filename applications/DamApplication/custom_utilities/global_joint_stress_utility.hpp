//
//   Project Name:        KratosDamApplication    $
//   Last modified by:    $Author: Lorenzo Gracia $
//   Date:                $Date:        June 2017 $
//   Revision:            $Revision:          0.0 $
//

#if !defined(KRATOS_GLOBAL_JOINT_STRESS_UTILITIES )
#define  KRATOS_GLOBAL_JOINT_STRESS_UTILITIES

#include <cmath>

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

#include "dam_application_variables.h"

namespace Kratos
{
    
class GlobalJointStressUtility
{
    
public:
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    GlobalJointStressUtility() {}

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    ~GlobalJointStressUtility() {}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Transforming local stress vector in global coordinates 
    void ComputingGlobalStress( ModelPart& r_model_part)
    {
        const int nelements = r_model_part.GetMesh().Elements().size();
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
        GeometryData::IntegrationMethod MyIntegrationMethod;
        const ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        boost::numeric::ublas::bounded_matrix<double,3,3> RotationPlane;

        //Definition of the plane
        array_1d<double, 3> pmid0;
        array_1d<double, 3> pmid1;
        array_1d<double, 3> pmid2;

        // Coordinates of interest plane
        pmid0[0]=247.67;
        pmid0[1]=-368.92;
        pmid0[2]=242.6;

        pmid1[0]=247.15;
        pmid1[1]=-385.15;
        pmid1[2]=242.6;

        pmid2[0]=313.88;
        pmid2[1]=-373.26;
        pmid2[2]=242.6;

        this->CalculateRotationMatrix(RotationPlane,pmid0,pmid1,pmid2);
            
        for(int k = 0; k<nelements; k++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;
            Element::GeometryType& Geom = it->GetGeometry();
            boost::numeric::ublas::bounded_matrix<double,3,3> RotationMatrix;

            //Define mid-plane points for prism_interface_3d_6
            noalias(pmid0) = 0.5 * (Geom.GetPoint( 0 ) + Geom.GetPoint( 3 ));
            noalias(pmid1) = 0.5 * (Geom.GetPoint( 1 ) + Geom.GetPoint( 4 ));
            noalias(pmid2) = 0.5 * (Geom.GetPoint( 2 ) + Geom.GetPoint( 5 ));

            this->CalculateRotationMatrix(RotationMatrix,pmid0,pmid1,pmid2);
            MyIntegrationMethod = it->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = Geom.IntegrationPoints(MyIntegrationMethod);
            unsigned int NumGPoints = IntegrationPoints.size();
            std::vector<array_1d<double,3>> LocalStressVector;
            array_1d<double,3> GlobalStressVector;
            it->GetValueOnIntegrationPoints(LOCAL_STRESS_VECTOR,LocalStressVector,CurrentProcessInfo);
            array_1d<double,3> TotalStressElement;

            for(unsigned int GPoint=0; GPoint<NumGPoints; GPoint++)
            {
                GlobalStressVector= prod(RotationMatrix,LocalStressVector[GPoint]);
                noalias(TotalStressElement) += GlobalStressVector;
            }

            GlobalStressVector = prod(RotationPlane,GlobalStressVector);

        }

    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Currently is just working for prism element

   void CalculateRotationMatrix(boost::numeric::ublas::bounded_matrix<double,3,3>& rRotationMatrix, array_1d<double, 3>& pmid0, array_1d<double, 3>& pmid1, array_1d<double, 3>& pmid2 )
    {
        KRATOS_TRY
        
        //Unitary vector in local x direction
        array_1d<double, 3> Vx;
        noalias(Vx) = pmid1 - pmid0;
        double inv_norm_x = 1.0/norm_2(Vx);
        Vx[0] *= inv_norm_x;
        Vx[1] *= inv_norm_x;
        Vx[2] *= inv_norm_x;
        
        //Unitary vector in local z direction
        array_1d<double, 3> Vy;
        noalias(Vy) = pmid2 - pmid0;
        array_1d<double, 3> Vz;
        MathUtils<double>::CrossProduct(Vz, Vx, Vy);
        double inv_norm_z = 1.0/norm_2(Vz);
        Vz[0] *= inv_norm_z;
        Vz[1] *= inv_norm_z;
        Vz[2] *= inv_norm_z;
            
        //Unitary vector in local y direction
        MathUtils<double>::CrossProduct( Vy, Vz, Vx);
    
        //Rotation Matrix
        rRotationMatrix(0,0) = Vx[0];
        rRotationMatrix(0,1) = Vx[1];
        rRotationMatrix(0,2) = Vx[2];
    
        rRotationMatrix(1,0) = Vy[0];
        rRotationMatrix(1,1) = Vy[1];
        rRotationMatrix(1,2) = Vy[2];
    
        rRotationMatrix(2,0) = Vz[0];
        rRotationMatrix(2,1) = Vz[1];
        rRotationMatrix(2,2) = Vz[2];
    
        KRATOS_CATCH( "" )
        
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_GLOBAL_JOINT_STRESS_UTILITIES defined */
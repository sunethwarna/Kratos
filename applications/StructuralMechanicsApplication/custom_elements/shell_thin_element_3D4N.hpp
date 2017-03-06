// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//

#if !defined(SHELL_THIN_ELEMENT_3D4N_H_INCLUDED )
#define  SHELL_THIN_ELEMENT_3D4N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "utilities/quaternion.h"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"



/*
// Project includes
#include "includes/element.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "utilities/quaternion.h"


#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"
#include "custom_utilities/shellq4_local_coordinate_system.hpp" mine!
//#include "geometries/quadrilateral_3d_4.h"

*/

namespace Kratos
{


	///@name Kratos Globals
	///@{
	///@}

	///@name Type Definitions
	///@{
	///@}
	//class Kratos::ShellQ4_CorotationalCoordinateTransformation;

	class ShellQ4_CoordinateTransformation;

	//class ShellQ4_LocalCoordinateSystem;

	///@name  Enum's
	///@{
	///@}

	///@name  Functions
	///@{
	///@}

	///@name Kratos Classes
	///@{

	/** \brief ShellThinElement3D4N
	*
	* This element represents a 4-node Shell element
	* based on the Assumed Natural DEviatoric Strain (ANDES) by Felippa.
	* This element is formulated for small strains,
	* but can be used in Geometrically nonlinear problems
	* involving large displacements and rotations
	* using a Corotational Coordinate Transformation.
	* Material nonlinearity is handled by means of the cross section object.
	*/
	class ShellThinElement3D4N : public Element
	{
	public:

		///@name Type Definitions
		///@{

		KRATOS_CLASS_POINTER_DEFINITION(ShellThinElement3D4N);

		typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

		typedef ShellQ4_CoordinateTransformation CoordinateTransformationBaseType;

		typedef boost::shared_ptr<CoordinateTransformationBaseType> CoordinateTransformationBasePointerType;

		typedef array_1d<double, 3> Vector3Type;

		//typedef Quaternion<double> QuaternionType;

		///@}

		///@name Classes
		///@{
		/** \brief JacobianOperator
		*
		* This class is a utility to compute at a given integration point,
		* the Jacobian, its inverse, its determinant
		* and the derivatives of the shape functions in the local
		* cartesian coordinate system.
		*/


		class JacobianOperator
		{
		public:

			JacobianOperator();

			void Calculate(const ShellQ4_LocalCoordinateSystem & CS, const Matrix & dN);

			inline const Matrix & Jacobian()const
			{
				return mJac;
			}

			inline const Matrix & Inverse()const
			{
				return mInv;
			}

			inline const Matrix & XYDerivatives()const
			{
				return mXYDeriv;
			}

			inline const double Determinant()const
			{
				return mDet;
			}

		private:

			Matrix mJac;     //!< Jacobian matrix
			Matrix mInv;    //!< Inverse of the Jacobian matrix 
			Matrix mXYDeriv; //*!< Shape function derivatives in cartesian coordinates
			double mDet;     //*!< Determinant of the Jacobian matrix 
		};



		// TODO: Add Calculation Data

		///@}

		///@name Life Cycle
		///@{

		ShellThinElement3D4N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			bool NLGeom = false);

		ShellThinElement3D4N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties,
			bool NLGeom = false);

		ShellThinElement3D4N(IndexType NewId,
			GeometryType::Pointer pGeometry,
			PropertiesType::Pointer pProperties,
			CoordinateTransformationBasePointerType pCoordinateTransformation);

		virtual ~ShellThinElement3D4N();

		///@}

		///@name Operations
		///@{

		// Basic

		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

		IntegrationMethod GetIntegrationMethod() const;

		void Initialize();

		//void ResetConstitutiveLaw();

		void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

		void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo);

		//int Check(const ProcessInfo& rCurrentProcessInfo);

		//void CleanMemory();

		//void GetValuesVector(Vector& values, int Step = 0);

		//void GetFirstDerivativesVector(Vector& values, int Step = 0);

		//void GetSecondDerivativesVector(Vector& values, int Step = 0);

		//void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		//void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		//void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

		//void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

		//void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

		//void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);

		void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo);

		//void CalculateRightHandSide(VectorType& rRightHandSideVector,
		//	ProcessInfo& rCurrentProcessInfo);

		// Results calculation on integration points
		/*
		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable, std::vector<array_1d<double, 3> >& rValues, const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6> >& rVariable, std::vector<array_1d<double, 6> >& rValues, const ProcessInfo& rCurrentProcessInfo);
		*/
		///@}

		///@name Public specialized Access - Temporary
		///@{

		//void SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections);

		///@}

	protected:

		///@name Protected Lyfe Cycle
		///@{

		/**
		* Protected empty constructor
		*/
		ShellThinElement3D4N() : Element()
		{
		}

		///@}

	private:

		///@name Private Classes
		///@{

		class CalculationData
		{

		public:

			// ---------------------------------------
			// calculation-constant data
			// ----------------------------------------
			// these data are allocated and constructed
			// at the beginning of the calculation

			ShellQ4_LocalCoordinateSystem LCS0; /*!< reference coordinate system */
			ShellQ4_LocalCoordinateSystem LCS;  /*!< current coordinate system */

			Vector s_xi, s_eta, s_24, s_13;
			array_1d<double, 4> N_pw;
			array_1d<Vector, 4> r_cartesian;
			array_1d<double, 4> dA;			/** \brief dA
											*
											* data.dA[PointNumber] = mDetJ0[PointNumber] * IntegrationWeight;
											*
											*/


			MatrixType L;
			MatrixType L_bend;
			MatrixType L_bend_mod;

			MatrixType Q1;
			MatrixType Q2;
			MatrixType Q3;
			MatrixType Q4;
			MatrixType QF;

			MatrixType B_h_1;
			MatrixType B_h_2;
			MatrixType B_h_3;
			MatrixType B_h_4;
			MatrixType B_h_bar;

			array_1d<double, 4> DKQ_a;
			array_1d<double, 4> DKQ_b;
			array_1d<double, 4> DKQ_c;
			array_1d<double, 4> DKQ_d;
			array_1d<double, 4> DKQ_e;
			MatrixType DKQ_indices;
			array_1d<Matrix, 4> DKQ_invJac;
			array_1d<Matrix, 4> DKQ_jac;
			array_1d<double, 4> DKQ_jac_det;

			MatrixType K_bend_temp;

			array_1d<Matrix, 4> B_i_bend;
			MatrixType B_h_bend_bar;
			MatrixType B_higher_test;
			MatrixType B_bend_DKQ;

			MatrixType H_mem;
			MatrixType H_mem_mod;

			MatrixType T_13 = Matrix(3, 3, 0.0);
			MatrixType T_24 = Matrix(3, 3, 0.0);

			double hMean;
			double TotalArea;
			double TotalVolume;
			double alpha = 1.5;
			std::vector< array_1d<double, 4> > gpLocations;

			MatrixType dNxy; /*!< shape function cartesian derivatives */ // Maybe DELETE????

			VectorType globalDisplacements; /*!< global displacement vector */
			VectorType localDisplacements;  /*!< local displacement vector */

			bool CalculateRHS; /*!< flag for the calculation of the right-hand-side vector */
			bool CalculateLHS; /*!< flag for the calculation of the left-hand-side vector */

							   // ---------------------------------------
							   // calculation-variable data
							   // ---------------------------------------
							   // these data are updated during the
							   // calculations

			double beta0;
			size_t gpIndex;

			// ---------------------------------------
			// calculation-variable data
			// ---------------------------------------
			// these data are updated during the
			// calculations, but they are allocated
			// only once(the first time they are used)
			// to avoid useless re-allocations

			MatrixType B;   /*!< total strain-displacement matrix at the current integration point */
			MatrixType D;   /*!< section constitutive matrix at the current integration point */
			MatrixType BTD; /*!< auxiliary matrix to store the product B'*D */

			VectorType generalizedStrains;  /*!< generalized strain vector at the current integration point */
			VectorType generalizedStresses; /*!< generalized stress vector at the current integration point */

			VectorType N; /*!< shape function vector at the current integration point */
			JacobianOperator jacOp;

			ShellCrossSection::Parameters SectionParameters; /*!< parameters for cross section calculations */

			array_1d< Vector3Type, 3 > Sig;

		public:

			const ProcessInfo& CurrentProcessInfo;

		public:

			CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
				const ProcessInfo& rCurrentProcessInfo);

		};

		///@}

		///@name Private Operations
		///@{

		//void DecimalCorrection(Vector& a);

		void SetupOrientationAngles();

		void InitializeCalculationData(CalculationData& data);

		void CalculateQMatrices(CalculationData& data);

		void CalculateB_h_mats(CalculationData& data);

		void CalculateBMatrix(CalculationData& data);

		void CalculateBeta0(CalculationData& data);

		void CalculateSectionResponse(CalculationData& data);

		void CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS);

		//void ApplyCorrectionToRHS(CalculationData& data, VectorType& RHS);

		//void AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector);

		void CalculateAll(MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo,
			const bool LHSrequired,
			const bool RHSrequired);

		/*
		bool TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double, 3> >& rVariable,
		std::vector<array_1d<double, 3> >& rValues,
		const ProcessInfo& rCurrentProcessInfo);

		bool TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo);
		*/

		void printMatrix(Matrix& matrixIn, std::string stringIn);

		void printVector(Vector& matrixIn, std::string stringIn);

		void printMatrix(const Matrix& matrixIn, std::string stringIn);

		void printDouble(std::string stringIn, double doubleIn);

		///@}

		///@name Static Member Variables
		///@{
		///@}

		///@name Member Variables
		///@{

		CoordinateTransformationBasePointerType mpCoordinateTransformation; /*!< The Coordinate Transformation */

		CrossSectionContainerType mSections; /*!< Container for cross section associated to each integration point */

		IntegrationMethod mThisIntegrationMethod; /*!< Currently selected integration method */

												  ///@}

												  ///@name Serialization
												  ///@{

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const;

		virtual void load(Serializer& rSerializer);

		///@}

		///@name Private  Access
		///@{
		///@}

		///@name Private Inquiry
		///@{
		///@}

		///@name Un accessible methods
		///@{
		///@}

	};

}
#endif // SHELL_THIN_ELEMENT_3D4N_H_INCLUDED

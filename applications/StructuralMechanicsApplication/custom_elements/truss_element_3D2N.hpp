// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:   
//                   
//                   
//

#if !defined(KRATOS_TRUSS_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_TRUSS_ELEMENT_3D2N_H_INCLUDED


#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{

	class TrussElement3D2N : public Element
	{
	public:
		KRATOS_CLASS_POINTER_DEFINITION(TrussElement3D2N);


		typedef Element BaseType;
		typedef BaseType::GeometryType GeometryType;
		typedef BaseType::NodesArrayType NodesArrayType;
		typedef BaseType::PropertiesType PropertiesType;
		typedef BaseType::IndexType IndexType;
		typedef BaseType::SizeType SizeType;
		typedef BaseType::MatrixType MatrixType;
		typedef BaseType::VectorType VectorType;
		typedef BaseType::EquationIdVectorType EquationIdVectorType;
		typedef BaseType::DofsVectorType DofsVectorType;


		TrussElement3D2N(IndexType NewId, 
						GeometryType::Pointer pGeometry);
		TrussElement3D2N(IndexType NewId,
						GeometryType::Pointer pGeometry,
						PropertiesType::Pointer pProperties);


		virtual ~TrussElement3D2N();


		BaseType::Pointer Create(
			IndexType NewId,
			NodesArrayType const& rThisNodes,
			PropertiesType::Pointer pProperties) const;

		void EquationIdVector(
			EquationIdVectorType& rResult,
			ProcessInfo& rCurrentProcessInfo) override;

		void GetDofList(
			DofsVectorType& rElementalDofList,
			ProcessInfo& rCurrentProcessInfo) override;

		void Initialize();

		MatrixType CreateElementStiffnessMatrix();

		void CalculateOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<double>& rVariable,
			std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void UpdateInternalForces(
			VectorType& rinternalForces);
		void CreateTransformationMatrix(
			Matrix& rRotationMatrix);
		double CalculateCurrentLength();

		void CalculateOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rOutput,
			const ProcessInfo& rCurrentProcessInfo) override;

		void GetValueOnIntegrationPoints(
			const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo) override;

		void CalculateLocalSystem(
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;


		void CalculateRightHandSide(
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateLeftHandSide(
			MatrixType& rLeftHandSideMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateMassMatrix(
			MatrixType& rMassMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void CalculateDampingMatrix(
			MatrixType& rDampingMatrix,
			ProcessInfo& rCurrentProcessInfo) override;

		void GetValuesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetSecondDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		void GetFirstDerivativesVector(
			Vector& rValues,
			int Step = 0) override;

		int  Check(
			const ProcessInfo& rCurrentProcessInfo);


		double CalculateGreenLagrangeStrain();
		double CalculateReferenceLength();

		VectorType CalculateBodyForces();  


	private:
		double mPreStress, mArea, mYoungsModulus, mLength, mDensity;
		double mCurrentLength;
		MatrixType mLHS;
		int mIterCount = 0; 
		bool mIsCable;
		bool mIsCompressed;

		TrussElement3D2N() {};

		friend class Serializer;
	};


}


#endif

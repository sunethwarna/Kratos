// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Peter Wilson
//					 (inspired by Massimo Petracca's shells)

#include "shell_thick_element_3D3N.hpp"
#include "custom_utilities/shellt3_corotational_coordinate_transformation.hpp"
#include "structural_mechanics_application_variables.h"

#include "geometries/triangle_3d_3.h"

#include <string>
#include <iomanip>

#define OPT_NUM_NODES 3
#define OPT_STRAIN_SIZE 6
#define OPT_NUM_DOFS 18

//----------------------------------------
// preprocessors for the integration
// method used by this element.

#define OPT_1_POINT_INTEGRATION

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_1
#define OPT_NUM_GP 1
#else
#define OPT_INTEGRATION_METHOD Kratos::GeometryData::GI_GAUSS_2
#define OPT_NUM_GP 3
#endif // OPT_1_POINT_INTEGRATION

//----------------------------------------
// preprocessors to handle the output
// in case of 3 integration points

//#define OPT_USES_INTERIOR_GAUSS_POINTS

#ifdef OPT_1_POINT_INTEGRATION
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X)
#else
#define OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(X) Utilities::InterpToStandardGaussPoints(X)
#endif // OPT_USES_INTERIOR_GAUSS_POINTS
#endif // OPT_1_POINT_INTEGRATION

//#define OPT_AVARAGE_RESULTS

namespace Kratos
{

	namespace Utilities
	{
		inline void InterpToStandardGaussPoints(double& v1, double& v2, double& v3)
		{
			double vg1 = v1;
			double vg2 = v2;
			double vg3 = v3;
#ifdef OPT_AVARAGE_RESULTS
			v1 = (vg1 + vg2 + vg3) / 3.0;
			v2 = (vg1 + vg2 + vg3) / 3.0;
			v3 = (vg1 + vg2 + vg3) / 3.0;
#else
			v1 = (2.0*vg1) / 3.0 - vg2 / 3.0 + (2.0*vg3) / 3.0;
			v2 = (2.0*vg1) / 3.0 + (2.0*vg2) / 3.0 - vg3 / 3.0;
			v3 = (2.0*vg2) / 3.0 - vg1 / 3.0 + (2.0*vg3) / 3.0;
#endif // OPT_AVARAGE_RESULTS
		}

		inline void InterpToStandardGaussPoints(std::vector< double >& v)
		{
			if (v.size() != 3) return;
			InterpToStandardGaussPoints(v[0], v[1], v[2]);
		}

		inline void InterpToStandardGaussPoints(std::vector< array_1d<double, 3> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 3; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< array_1d<double, 6> >& v)
		{
			if (v.size() != 3) return;
			for (size_t i = 0; i < 6; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< Vector >& v)
		{
			if (v.size() != 3) return;
			size_t ncomp = v[0].size();
			for (int i = 1; i < 3; i++)
				if (v[i].size() != ncomp)
					return;
			for (size_t i = 0; i < ncomp; i++)
				InterpToStandardGaussPoints(v[0][i], v[1][i], v[2][i]);
		}

		inline void InterpToStandardGaussPoints(std::vector< Matrix >& v)
		{
			if (v.size() != 3) return;
			size_t nrows = v[0].size1();
			size_t ncols = v[0].size2();
			for (int i = 1; i < 3; i++)
				if (v[i].size1() != nrows || v[i].size2() != ncols)
					return;
			for (size_t i = 0; i < nrows; i++)
				for (size_t j = 0; j < ncols; j++)
					InterpToStandardGaussPoints(v[0](i, j), v[1](i, j), v[2](i, j));
		}

	}

	// =====================================================================================
	//
	// CalculationData
	//
	// =====================================================================================

	ShellThickElement3D3N::CalculationData::CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
		const ProcessInfo& rCurrentProcessInfo)
		: LCS0(pCoordinateTransformation->CreateReferenceCoordinateSystem())
		, LCS(pCoordinateTransformation->CreateLocalCoordinateSystem())
		, CurrentProcessInfo(rCurrentProcessInfo)

	{
	}

	// =====================================================================================
	//
	// Class ShellThickElement3D3N
	//
	// =====================================================================================

	ShellThickElement3D3N::ShellThickElement3D3N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		bool NLGeom)
		: Element(NewId, pGeometry)
		, mpCoordinateTransformation(NLGeom ?
			new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
			new ShellT3_CoordinateTransformation(pGeometry))
	{
		mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
	}

	ShellThickElement3D3N::ShellThickElement3D3N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties,
		bool NLGeom)
		: Element(NewId, pGeometry, pProperties)
		, mpCoordinateTransformation(NLGeom ?
			new ShellT3_CorotationalCoordinateTransformation(pGeometry) :
			new ShellT3_CoordinateTransformation(pGeometry))
	{
		mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
	}

	ShellThickElement3D3N::ShellThickElement3D3N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties,
		CoordinateTransformationBasePointerType pCoordinateTransformation)
		: Element(NewId, pGeometry, pProperties)
		, mpCoordinateTransformation(pCoordinateTransformation)
	{
		mThisIntegrationMethod = OPT_INTEGRATION_METHOD;
	}

	ShellThickElement3D3N::~ShellThickElement3D3N()
	{
	}

	Element::Pointer ShellThickElement3D3N::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
		return boost::make_shared< ShellThickElement3D3N >(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom));
		//     return Element::Pointer( new ShellThickElement3D3N(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom)) );
	}

	ShellThickElement3D3N::IntegrationMethod ShellThickElement3D3N::GetIntegrationMethod() const
	{
		return mThisIntegrationMethod;
	}

	void ShellThickElement3D3N::Initialize()
	{
		KRATOS_TRY

			const GeometryType & geom = GetGeometry();
		const PropertiesType & props = GetProperties();

		if (geom.PointsNumber() != OPT_NUM_NODES)
			KRATOS_THROW_ERROR(std::logic_error, "ShellThickElement3D3N Element - Wrong number of nodes", geom.PointsNumber());

		const GeometryType::IntegrationPointsArrayType & integrationPoints = geom.IntegrationPoints(GetIntegrationMethod());
		if (integrationPoints.size() != OPT_NUM_GP)
			KRATOS_THROW_ERROR(std::logic_error, "ShellThickElement3D3N Element - Wrong integration scheme", integrationPoints.size());

		if (mSections.size() != OPT_NUM_GP)
		{
			const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

			ShellCrossSection::Pointer theSection;
			if (props.Has(SHELL_CROSS_SECTION))
			{
				theSection = props[SHELL_CROSS_SECTION];
			}
			else
			{
				theSection = ShellCrossSection::Pointer(new ShellCrossSection());
				theSection->BeginStack();
				theSection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
				theSection->EndStack();
			}

			mSections.clear();
			for (int i = 0; i < OPT_NUM_GP; i++)
			{
				ShellCrossSection::Pointer sectionClone = theSection->Clone();
				sectionClone->SetSectionBehavior(ShellCrossSection::Thick);
				sectionClone->InitializeCrossSection(props, geom, row(shapeFunctionsValues, i));
				mSections.push_back(sectionClone);
			}
		}

		mpCoordinateTransformation->Initialize();

		this->SetupOrientationAngles();

		KRATOS_CATCH("")
	}

	void ShellThickElement3D3N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		if (rResult.size() != OPT_NUM_DOFS)
			rResult.resize(OPT_NUM_DOFS, false);

		GeometryType & geom = this->GetGeometry();

		for (SizeType i = 0; i < geom.size(); i++)
		{
			int index = i * 6;
			NodeType & iNode = geom[i];

			rResult[index] = iNode.GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = iNode.GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = iNode.GetDof(DISPLACEMENT_Z).EquationId();

			rResult[index + 3] = iNode.GetDof(ROTATION_X).EquationId();
			rResult[index + 4] = iNode.GetDof(ROTATION_Y).EquationId();
			rResult[index + 5] = iNode.GetDof(ROTATION_Z).EquationId();
		}
	}

	void ShellThickElement3D3N::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
	{
		ElementalDofList.resize(0);
		ElementalDofList.reserve(OPT_NUM_DOFS);

		GeometryType & geom = this->GetGeometry();

		for (SizeType i = 0; i < geom.size(); i++)
		{
			NodeType & iNode = geom[i];

			ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(iNode.pGetDof(DISPLACEMENT_Z));

			ElementalDofList.push_back(iNode.pGetDof(ROTATION_X));
			ElementalDofList.push_back(iNode.pGetDof(ROTATION_Y));
			ElementalDofList.push_back(iNode.pGetDof(ROTATION_Z));
		}
	}

	int ShellThickElement3D3N::Check(const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			GeometryType& geom = GetGeometry();

		// verify that the variables are correctly initialized
		if (DISPLACEMENT.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "");
		if (ROTATION.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "ROTATION has Key zero! (check if the application is correctly registered", "");
		if (VELOCITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "");
		if (ACCELERATION.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "");
		if (DENSITY.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "");
		if (SHELL_CROSS_SECTION.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "SHELL_CROSS_SECTION has Key zero! (check if the application is correctly registered", "");
		if (THICKNESS.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "");
		if (CONSTITUTIVE_LAW.Key() == 0)
			KRATOS_THROW_ERROR(std::invalid_argument, "CONSTITUTIVE_LAW has Key zero! (check if the application is correctly registered", "");

		// verify that the dofs exist
		for (unsigned int i = 0; i<geom.size(); i++)
		{
			if (geom[i].SolutionStepsDataHas(DISPLACEMENT) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "missing variable DISPLACEMENT on node ", geom[i].Id());
			if (geom[i].HasDofFor(DISPLACEMENT_X) == false || geom[i].HasDofFor(DISPLACEMENT_Y) == false || geom[i].HasDofFor(DISPLACEMENT_Z) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id());
			if (geom[i].SolutionStepsDataHas(ROTATION) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "missing variable ROTATION on node ", geom[i].Id());
			if (geom[i].HasDofFor(ROTATION_X) == false || geom[i].HasDofFor(ROTATION_Y) == false || geom[i].HasDofFor(ROTATION_Z) == false)
				KRATOS_THROW_ERROR(std::invalid_argument, "missing one of the dofs for the variable ROTATION on node ", geom[i].Id());

			if (geom[i].GetBufferSize() < 2)
				KRATOS_THROW_ERROR(std::logic_error, "This Element needs at least a buffer size = 2", "");
		}

		// check properties
		if (this->pGetProperties() == NULL)
			KRATOS_THROW_ERROR(std::logic_error, "Properties not provided for element ", this->Id());

		const PropertiesType & props = this->GetProperties();

		if (props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
		{
			const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
			if (section == NULL)
				KRATOS_THROW_ERROR(std::logic_error, "SHELL_CROSS_SECTION not provided for element ", this->Id());

			section->Check(props, geom, rCurrentProcessInfo);
		}
		else // ... allow the automatic creation of a homogeneous section from a material and a thickness
		{
			if (!props.Has(CONSTITUTIVE_LAW))
				KRATOS_THROW_ERROR(std::logic_error, "CONSTITUTIVE_LAW not provided for element ", this->Id());
			const ConstitutiveLaw::Pointer& claw = props[CONSTITUTIVE_LAW];
			if (claw == NULL)
				KRATOS_THROW_ERROR(std::logic_error, "CONSTITUTIVE_LAW not provided for element ", this->Id());

			if (!props.Has(THICKNESS))
				KRATOS_THROW_ERROR(std::logic_error, "THICKNESS not provided for element ", this->Id());
			if (props[THICKNESS] <= 0.0)
				KRATOS_THROW_ERROR(std::logic_error, "wrong THICKNESS value provided for element ", this->Id());

			ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
			dummySection->BeginStack();
			dummySection->AddPly(props[THICKNESS], 0.0, 5, this->pGetProperties());
			dummySection->EndStack();
			dummySection->SetSectionBehavior(ShellCrossSection::Thick);
			dummySection->Check(props, geom, rCurrentProcessInfo);
		}

		return 0;

		KRATOS_CATCH("")
	}

	void ShellThickElement3D3N::CleanMemory()
	{
	}

	void ShellThickElement3D3N::GetValuesVector(Vector& values, int Step)
	{
		if (values.size() != OPT_NUM_DOFS)
			values.resize(OPT_NUM_DOFS, false);

		const GeometryType & geom = GetGeometry();

		for (SizeType i = 0; i < geom.size(); i++)
		{
			const NodeType & iNode = geom[i];
			const array_1d<double, 3>& disp = iNode.FastGetSolutionStepValue(DISPLACEMENT, Step);
			const array_1d<double, 3>& rot = iNode.FastGetSolutionStepValue(ROTATION, Step);

			int index = i * 6;
			values[index] = disp[0];
			values[index + 1] = disp[1];
			values[index + 2] = disp[2];

			values[index + 3] = rot[0];
			values[index + 4] = rot[1];
			values[index + 5] = rot[2];
		}
	}

	void ShellThickElement3D3N::GetFirstDerivativesVector(Vector& values, int Step)
	{
		if (values.size() != OPT_NUM_DOFS)
			values.resize(OPT_NUM_DOFS, false);

		const GeometryType & geom = GetGeometry();

		for (SizeType i = 0; i < geom.size(); i++)
		{
			const NodeType & iNode = geom[i];
			const array_1d<double, 3>& vel = iNode.FastGetSolutionStepValue(VELOCITY, Step);

			int index = i * 6;
			values[index] = vel[0];
			values[index + 1] = vel[1];
			values[index + 2] = vel[2];
			values[index + 3] = 0.0;
			values[index + 4] = 0.0;
			values[index + 5] = 0.0;
		}
	}

	void ShellThickElement3D3N::GetSecondDerivativesVector(Vector& values, int Step)
	{
		if (values.size() != OPT_NUM_DOFS)
			values.resize(OPT_NUM_DOFS, false);

		const GeometryType & geom = GetGeometry();

		for (SizeType i = 0; i < geom.size(); i++)
		{
			const NodeType & iNode = geom[i];
			const array_1d<double, 3>& acc = iNode.FastGetSolutionStepValue(ACCELERATION, Step);

			int index = i * 6;
			values[index] = acc[0];
			values[index + 1] = acc[1];
			values[index + 2] = acc[2];
			values[index + 3] = 0.0;
			values[index + 4] = 0.0;
			values[index + 5] = 0.0;
		}
	}

	void ShellThickElement3D3N::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
		mpCoordinateTransformation->InitializeNonLinearIteration(CurrentProcessInfo);

		const GeometryType & geom = this->GetGeometry();
		const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
		for (SizeType i = 0; i < mSections.size(); i++)
			mSections[i]->InitializeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), CurrentProcessInfo);
	}

	void ShellThickElement3D3N::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
		mpCoordinateTransformation->FinalizeNonLinearIteration(CurrentProcessInfo);

		const GeometryType & geom = this->GetGeometry();
		const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
		for (SizeType i = 0; i < mSections.size(); i++)
			mSections[i]->FinalizeNonLinearIteration(GetProperties(), geom, row(shapeFunctionsValues, i), CurrentProcessInfo);
	}

	void ShellThickElement3D3N::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		const PropertiesType& props = GetProperties();
		const GeometryType & geom = GetGeometry();
		const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

		for (SizeType i = 0; i < mSections.size(); i++)
			mSections[i]->InitializeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

		mpCoordinateTransformation->InitializeSolutionStep(CurrentProcessInfo);
	}

	void ShellThickElement3D3N::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		const PropertiesType& props = GetProperties();
		const GeometryType& geom = GetGeometry();
		const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());

		for (SizeType i = 0; i < mSections.size(); i++)
			mSections[i]->FinalizeSolutionStep(props, geom, row(shapeFunctionsValues, i), CurrentProcessInfo);

		mpCoordinateTransformation->FinalizeSolutionStep(CurrentProcessInfo);
	}

	void ShellThickElement3D3N::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		if ((rMassMatrix.size1() != OPT_NUM_DOFS) || (rMassMatrix.size2() != OPT_NUM_DOFS))
			rMassMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);
		noalias(rMassMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);

		// Compute the local coordinate system.

		ShellT3_LocalCoordinateSystem referenceCoordinateSystem(
			mpCoordinateTransformation->CreateReferenceCoordinateSystem());

		// lumped area

		double lump_area = referenceCoordinateSystem.Area() / 3.0;

		// Calculate avarage mass per unit area
		double av_mass_per_unit_area = 0.0;
		for (size_t i = 0; i < OPT_NUM_GP; i++)	//this might be just 1
			av_mass_per_unit_area += mSections[i]->CalculateMassPerUnitArea();
		av_mass_per_unit_area /= double(OPT_NUM_GP);

		// loop on nodes
		for (size_t i = 0; i < 3; i++)
		{
			size_t index = i * 6;

			double nodal_mass = av_mass_per_unit_area * lump_area;

			// translational mass
			rMassMatrix(index, index) = nodal_mass;
			rMassMatrix(index + 1, index + 1) = nodal_mass;
			rMassMatrix(index + 2, index + 2) = nodal_mass;

			// rotational mass - neglected for the moment...
		}
	}

	void ShellThickElement3D3N::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		if ((rDampingMatrix.size1() != OPT_NUM_DOFS) || (rDampingMatrix.size2() != OPT_NUM_DOFS))
			rDampingMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);

		noalias(rDampingMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);
	}


	void ShellThickElement3D3N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	void ShellThickElement3D3N::CalculateRightHandSide(VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		Matrix dummy;
		CalculateAll(dummy, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	// =====================================================================================
	//
	// Class ShellThickElement3D3N - Results on Gauss Points
	//
	// =====================================================================================

	void ShellThickElement3D3N::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
		std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rValues.size() != OPT_NUM_GP)
			rValues.resize(OPT_NUM_GP);

		for (int i = 0; i < OPT_NUM_GP; i++)
			mSections[i]->GetValue(rVariable, rValues[i]);

		OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
	}

	void ShellThickElement3D3N::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ShellThickElement3D3N::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(rVariable, rValues, rCurrentProcessInfo)) return;
	}

	void ShellThickElement3D3N::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3> >& rVariable,
		std::vector<array_1d<double, 3> >& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (TryGetValueOnIntegrationPoints_MaterialOrientation(rVariable, rValues, rCurrentProcessInfo)) return;
	}

	void ShellThickElement3D3N::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6> >& rVariable,
		std::vector<array_1d<double, 6> >& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
	}

	void ShellThickElement3D3N::SetCrossSectionsOnIntegrationPoints(std::vector< ShellCrossSection::Pointer >& crossSections)
	{
		KRATOS_TRY
			if (crossSections.size() != OPT_NUM_GP)
				KRATOS_THROW_ERROR(std::logic_error, "The number of cross section is wrong", crossSections.size());
		mSections.clear();
		for (SizeType i = 0; i <crossSections.size(); i++)
			mSections.push_back(crossSections[i]);
		this->SetupOrientationAngles();
		KRATOS_CATCH("")
	}



	// =====================================================================================
	//
	// Class ShellThickElement3D3N - Private methods
	//
	// =====================================================================================

	void ShellThickElement3D3N::DecimalCorrection(Vector& a)
	{
		double norm = norm_2(a);
		double tolerance = std::max(norm * 1.0E-12, 1.0E-12);
		for (SizeType i = 0; i < a.size(); i++)
			if (std::abs(a(i)) < tolerance)
				a(i) = 0.0;
	}

	void ShellThickElement3D3N::SetupOrientationAngles()
	{
		ShellT3_LocalCoordinateSystem lcs(mpCoordinateTransformation->CreateReferenceCoordinateSystem());

		Vector3Type normal;
		noalias(normal) = lcs.Vz();

		Vector3Type dZ;
		dZ(0) = 0.0;
		dZ(1) = 0.0;
		dZ(2) = 1.0; // for the moment let's take this. But the user can specify its own triad! TODO

		Vector3Type dirX;
		MathUtils<double>::CrossProduct(dirX, dZ, normal);

		// try to normalize the x vector. if it is near zero it means that we need
		// to choose a default one.
		double dirX_norm = dirX(0)*dirX(0) + dirX(1)*dirX(1) + dirX(2)*dirX(2);
		if (dirX_norm < 1.0E-12)
		{
			dirX(0) = 1.0;
			dirX(1) = 0.0;
			dirX(2) = 0.0;
		}
		else if (dirX_norm != 1.0)
		{
			dirX_norm = std::sqrt(dirX_norm);
			dirX /= dirX_norm;
		}

		Vector3Type elem_dirX = lcs.Vx();

		// now calculate the angle between the element x direction and the material x direction.
		Vector3Type& a = elem_dirX;
		Vector3Type& b = dirX;
		double a_dot_b = a(0)*b(0) + a(1)*b(1) + a(2)*b(2);
		if (a_dot_b < -1.0) a_dot_b = -1.0;
		if (a_dot_b >  1.0) a_dot_b = 1.0;
		double angle = std::acos(a_dot_b);

		// if they are not counter-clock-wise, let's change the sign of the angle
		if (angle != 0.0)
		{
			const MatrixType& R = lcs.Orientation();
			if (dirX(0)*R(1, 0) + dirX(1)*R(1, 1) + dirX(2)*R(1, 2) < 0.0)
				angle = -angle;
		}

		for (CrossSectionContainerType::iterator it = mSections.begin(); it != mSections.end(); ++it)
			(*it)->SetOrientationAngle(angle);
	}
	
	void ShellThickElement3D3N::CalculateSectionResponse(CalculationData& data)
	{
		const array_1d<double, 3>& loc = data.gpLocations[0];
		data.N(0) = 1.0 - loc[1] - loc[2];
		data.N(1) = loc[1];
		data.N(2) = loc[2];

		ShellCrossSection::Pointer& section = mSections[0];
		data.SectionParameters.SetShapeFunctionsValues(data.N);
		section->CalculateSectionResponse(data.SectionParameters, ConstitutiveLaw::StressMeasure_PK2);
		
		if (data.basicTriCST == false)
		{
			//add in shear stabilization
			double h2 = data.hMean*data.hMean;
			double shearStabilisation = (h2) / (h2 + data.alpha*data.h_e*data.h_e);
			data.D(6, 6) *= shearStabilisation;
			data.D(6, 7) *= shearStabilisation;
			data.D(7, 6) *= shearStabilisation;
			data.D(7, 7) *= shearStabilisation;
		}
	}

	void ShellThickElement3D3N::InitializeCalculationData(CalculationData& data)
	{
		//-------------------------------------
		// Computation of all stuff that remain
		// constant throughout the calculations

		//-------------------------------------
		// geometry data

		const double x12 = data.LCS0.X1() - data.LCS0.X2();
		const double x23 = data.LCS0.X2() - data.LCS0.X3();
		const double x31 = data.LCS0.X3() - data.LCS0.X1();
		const double x21 = -x12;
		const double x32 = -x23;
		const double x13 = -x31;

		const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
		const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
		const double y31 = data.LCS0.Y3() - data.LCS0.Y1();
		const double y21 = -y12;
		const double y32 = -y23;
		const double y13 = -y31;

		const double A = 0.5*(y21*x13 - x21*y13);
		const double A2 = 2.0*A;

		// Note: here we compute the avarage thickness,
		// since L is constant over the element.
		// Now it is not necessary to compute the avarage
		// because the current implementation of the cross section
		// doesn't have a variable thickness
		// (for example as a function of the spatial coordinates...).
		// This is just a place-holder for future
		// implementation of a variable thickness

		double h = 0.0;
		for (unsigned int i = 0; i < mSections.size(); i++)
			h += mSections[i]->GetThickness();
		h /= (double)mSections.size();

		data.hMean = h;
		data.TotalArea = A;
		data.TotalVolume = A * h;
		data.dA = A;

		// this is the integration weight
		// used during the gauss loop.
		// it is dArea because it will
		// multiply section stress resultants
		// and section constitutive matrices
		// that already take into accout the
		// thickness

		// crete the integration point locations
		if (data.gpLocations.size() != 0) data.gpLocations.clear();
		data.gpLocations.resize(OPT_NUM_GP);
#ifdef OPT_1_POINT_INTEGRATION
		array_1d<double, 3>& gp0 = data.gpLocations[0];
		gp0[0] = 1.0 / 3.0;
		gp0[1] = 1.0 / 3.0;
		gp0[2] = 1.0 / 3.0;
#else
		array_1d<double, 3>& gp0 = data.gpLocations[0];
		array_1d<double, 3>& gp1 = data.gpLocations[1];
		array_1d<double, 3>& gp2 = data.gpLocations[2];
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
		gp0[0] = 1.0 / 6.0;
		gp0[1] = 1.0 / 6.0;
		gp0[2] = 2.0 / 3.0;
		gp1[0] = 2.0 / 3.0;
		gp1[1] = 1.0 / 6.0;
		gp1[2] = 1.0 / 6.0;
		gp2[0] = 1.0 / 6.0;
		gp2[1] = 2.0 / 3.0;
		gp2[2] = 1.0 / 6.0;
#else
		gp0[0] = 0.5;
		gp0[1] = 0.5;
		gp0[2] = 0.0;
		gp1[0] = 0.0;
		gp1[1] = 0.5;
		gp1[2] = 0.5;
		gp2[0] = 0.5;
		gp2[1] = 0.0;
		gp2[2] = 0.5;
#endif // OPT_USES_INTERIOR_GAUSS_POINTS

#endif // OPT_1_POINT_INTEGRATION

		// cartesian derivatives
		data.dNxy.resize(3, 2, false);
		data.dNxy(0, 0) = (y13 - y12) / A2;
		data.dNxy(0, 1) = (x12 - x13) / A2;
		data.dNxy(1, 0) = -y13 / A2;
		data.dNxy(1, 1) = x13 / A2;
		data.dNxy(2, 0) = y12 / A2;
		data.dNxy(2, 1) = -x12 / A2;

		data.N.resize(3, false);

		

	
		// --------------------------------------------------------------------------
		// Total Formulation - as per Efficient Co-Rotational 3-Node Shell Element paper (2016)
		// --------------------------------------------------------------------------

		data.B.resize(8, 18, false);
		data.B.clear();

		//Membrane components
			//node 1
			data.B(0, 0) = y23;
			data.B(1, 1) = x32;
			data.B(2, 0) = x32;
			data.B(2, 1) = y23;

			//node 2
			data.B(0, 6) = y31;
			data.B(1, 7) = x13;
			data.B(2, 6) = x13;
			data.B(2, 7) = y31;

			//node 3
			data.B(0, 12) = y12;
			data.B(1, 13) = x21;
			data.B(2, 12) = x21;
			data.B(2, 13) = y12;

		//Bending components
			//node 1
			data.B(3, 4) = y23;
			data.B(4, 3) = -1.0 * x32;
			data.B(5, 3) = -1.0 * y23;
			data.B(5, 4) = x32;

			//node 2
			data.B(3, 10) = y31;
			data.B(4, 9) = -1.0 * x13;
			data.B(5, 9) = -1.0 * y31;
			data.B(5, 10) = x13;

			//node 3
			data.B(3, 16) = y12;
			data.B(4, 15) = -1.0 * x21;
			data.B(5, 15) = -1.0 * y12;
			data.B(5, 16) = x21;

		//Shear components
			if (data.basicTriCST == false)
			{
				// Use DSG method

				const double a = x21;
				const double b = y21;
				const double c = y31;
				const double d = x31;
				//node 1
				data.B(6, 2) = b - c;
				data.B(6, 4) = A;

				data.B(7, 2) = d - a;
				data.B(7, 3) = -1.0 * A;

				//node 2
				data.B(6, 8) = c;
				data.B(6, 9) = -1.0 * b*c / 2.0;
				data.B(6, 10) = a*c / 2.0;

				data.B(7, 8) = -1.0 * d;
				data.B(7, 9) = b*d / 2.0;
				data.B(7, 10) = -1.0 * a*d / 2.0;

				//node 3
				data.B(6, 14) = -1.0 * b;
				data.B(6, 15) = b*c / 2.0;
				data.B(6, 16) = b*d / 2.0;

				data.B(7, 14) = a;
				data.B(7, 15) = -1.0 * a*c / 2.0;
				data.B(7, 16) = a*d / 2.0;
			}
			else
			{
				// Basic CST displacement derived shear
				// strain displacement matrix.
				// Only for testing!

				const Matrix & shapeFunctions = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

				//node 1
				data.B(6, 2) = y23;
				data.B(6, 3) = shapeFunctions(0, 0);

				data.B(7, 2) = x32;
				data.B(7, 4) = shapeFunctions(0, 0);

				//node 2
				data.B(6, 8) = y31;
				data.B(6, 9) = shapeFunctions(0, 1);

				data.B(7, 8) = x13;
				data.B(7, 10) = shapeFunctions(0, 1);

				//node 3
				data.B(6, 14) = y12;
				data.B(6, 15) = shapeFunctions(0, 2);

				data.B(7, 14) = x21;
				data.B(7, 16) = shapeFunctions(0, 2);
			}
			

		//Final multiplication
			data.B /= (A2);

		//determine longest side length (eqn 22) - for shear correction
		Vector P12 = Vector(data.LCS0.P1() - data.LCS0.P2());
		Vector P13 = Vector(data.LCS0.P1() - data.LCS0.P3());
		Vector P23 = Vector(data.LCS0.P2() - data.LCS0.P3());
		data.h_e = std::sqrt(inner_prod(P12, P12));
		double edge_length = std::sqrt(inner_prod(P13, P13));
		if (edge_length > data.h_e) { data.h_e = edge_length; }
		edge_length = std::sqrt(inner_prod(P23, P23));
		if (edge_length > data.h_e) { data.h_e = edge_length; }


		//--------------------------------------
		// Calculate material matrices
		// 
		
		//allocate and setup
		data.D.resize(8, 8, false);
		data.generalizedStrains.resize(8, false);
		data.generalizedStresses.resize(8, false);
		data.SectionParameters.SetElementGeometry(GetGeometry());
		data.SectionParameters.SetMaterialProperties(GetProperties());
		data.SectionParameters.SetProcessInfo(data.CurrentProcessInfo);
		data.SectionParameters.SetGeneralizedStrainVector(data.generalizedStrains);
		data.SectionParameters.SetGeneralizedStressVector(data.generalizedStresses);
		data.SectionParameters.SetConstitutiveMatrix(data.D);
		data.SectionParameters.SetShapeFunctionsDerivatives(data.dNxy);
		Flags& options = data.SectionParameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, data.CalculateLHS);

		//--------------------------------------
		// calculate the displacement vector
		// in global and local coordinate systems

		data.globalDisplacements.resize(OPT_NUM_DOFS, false);
		GetValuesVector(data.globalDisplacements);

		data.localDisplacements =
			mpCoordinateTransformation->CalculateLocalDisplacements(
				data.LCS, data.globalDisplacements);
	}

	void ShellThickElement3D3N::AddBodyForces(CalculationData& data, VectorType& rRightHandSideVector)
	{
		const GeometryType& geom = GetGeometry();

		// Get shape functions
#ifdef OPT_USES_INTERIOR_GAUSS_POINTS
		const Matrix & N = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
#else
		Matrix N(3, 3);
		for (unsigned int igauss = 0; igauss < OPT_NUM_GP; igauss++)
		{
			const array_1d<double, 3>& loc = data.gpLocations[igauss];
			N(igauss, 0) = 1.0 - loc[1] - loc[2];
			N(igauss, 1) = loc[1];
			N(igauss, 2) = loc[2];
		}
#endif // !OPT_USES_INTERIOR_GAUSS_POINTS

		// auxiliary
		array_1d<double, 3> bf;

		// gauss loop to integrate the external force vector
		for (unsigned int igauss = 0; igauss < OPT_NUM_GP; igauss++)
		{
			// get mass per unit area
			double mass_per_unit_area = mSections[igauss]->CalculateMassPerUnitArea();
			

			// interpolate nodal volume accelerations to this gauss point
			// and obtain the body force vector
			bf.clear();
			for (unsigned int inode = 0; inode < 3; inode++)
			{
				if (geom[inode].SolutionStepsDataHas(VOLUME_ACCELERATION)) //temporary, will be checked once at the beginning only
				{
					bf += N(igauss, inode) * geom[inode].FastGetSolutionStepValue(VOLUME_ACCELERATION);
					
				}
			}
			bf *= (mass_per_unit_area * data.dA);
			

			// add it to the RHS vector
			for (unsigned int inode = 0; inode < 3; inode++)
			{
				unsigned int index = inode * 6;
				double iN = N(igauss, inode);
				rRightHandSideVector[index + 0] += iN * bf[0];
				rRightHandSideVector[index + 1] += iN * bf[1];
				rRightHandSideVector[index + 2] += iN * bf[2];
				
			}
			
		}
	}

	void ShellThickElement3D3N::CalculateAll(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool LHSrequired,
		const bool RHSrequired)
	{
		KRATOS_TRY
		// Resize the Left Hand Side if necessary,
		// and initialize it to Zero

			if ((rLeftHandSideMatrix.size1() != OPT_NUM_DOFS) || (rLeftHandSideMatrix.size2() != OPT_NUM_DOFS))
				rLeftHandSideMatrix.resize(OPT_NUM_DOFS, OPT_NUM_DOFS, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(OPT_NUM_DOFS, OPT_NUM_DOFS);

		// Resize the Right Hand Side if necessary,
		// and initialize it to Zero

		if (rRightHandSideVector.size() != OPT_NUM_DOFS)
			rRightHandSideVector.resize(OPT_NUM_DOFS, false);
		noalias(rRightHandSideVector) = ZeroVector(OPT_NUM_DOFS);

		// Initialize common calculation variables
		//std::cout << "before initialize" << std::endl;
		CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
		data.CalculateLHS = LHSrequired;
		data.CalculateRHS = RHSrequired;
		InitializeCalculationData(data);
		CalculateSectionResponse(data);

		// Calulate element stiffness
		Matrix BTD = Matrix(18, 8, 0.0);
		data.D *= data.TotalArea;
		BTD = prod(trans(data.B), data.D);
		noalias(rLeftHandSideMatrix) += prod(BTD, data.B);

		//add in z_rot artificial stiffness
		double z_rot_multiplier = 0.001;
		double max_stiff = 0.0; //max diagonal stiffness
		for (int i = 0; i < 18; i++)
		{
			if (rLeftHandSideMatrix(i, i) > max_stiff) { max_stiff = rLeftHandSideMatrix(i, i); }
		}
		for (int i = 0; i < 3; i++)
		{
			rLeftHandSideMatrix(6 * i + 5, 6 * i + 5) = z_rot_multiplier*max_stiff;
		}

		// Add RHS term
		rRightHandSideVector -= prod(rLeftHandSideMatrix, data.localDisplacements);

		// Let the CoordinateTransformation finalize the calculation.
		// This will handle the transformation of the local matrices/vectors to
		// the global coordinate system.
		mpCoordinateTransformation->FinalizeCalculations(data.LCS,
			data.globalDisplacements,
			data.localDisplacements,
			rLeftHandSideMatrix,
			rRightHandSideVector,
			RHSrequired,
			LHSrequired);

		// Add body forces contributions. This doesn't depend on the coordinate system
		AddBodyForces(data, rRightHandSideVector);
		KRATOS_CATCH("")
	}

	bool ShellThickElement3D3N::TryGetValueOnIntegrationPoints_MaterialOrientation(const Variable<array_1d<double, 3> >& rVariable,
		std::vector<array_1d<double, 3> >& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		// Check the required output

		int ijob = 0;
		if (rVariable == MATERIAL_ORIENTATION_DX)
			ijob = 1;
		else if (rVariable == MATERIAL_ORIENTATION_DY)
			ijob = 2;
		else if (rVariable == MATERIAL_ORIENTATION_DZ)
			ijob = 3;

		// quick return

		if (ijob == 0) return false;

		// resize output

		if (rValues.size() != OPT_NUM_GP)
			rValues.resize(OPT_NUM_GP);

		// Compute the local coordinate system.

		ShellT3_LocalCoordinateSystem localCoordinateSystem(
			mpCoordinateTransformation->CreateLocalCoordinateSystem());

		Vector3Type eZ = localCoordinateSystem.Vz();

		// Gauss Loop

		if (ijob == 1)
		{
			Vector3Type eX = localCoordinateSystem.Vx();
			for (int i = 0; i < OPT_NUM_GP; i++)
			{
				QuaternionType q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
				q.RotateVector3(eX, rValues[i]);
			}
		}
		else if (ijob == 2)
		{
			Vector3Type eY = localCoordinateSystem.Vy();
			for (int i = 0; i < OPT_NUM_GP; i++)
			{
				QuaternionType q = QuaternionType::FromAxisAngle(eZ(0), eZ(1), eZ(2), mSections[i]->GetOrientationAngle());
				q.RotateVector3(eY, rValues[i]);
			}
		}
		else if (ijob == 3)
		{
			for (int i = 0; i < OPT_NUM_GP; i++)
			{
				noalias(rValues[i]) = eZ;
			}
		}

		return true;
	}

	bool ShellThickElement3D3N::TryGetValueOnIntegrationPoints_GeneralizedStrainsOrStresses(const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		// Check the required output
		
		int ijob = 0;
		bool bGlobal = false;
		if (rVariable == SHELL_STRAIN)
		{
			ijob = 1;
		}
		else if (rVariable == SHELL_STRAIN_GLOBAL)
		{
			ijob = 1;
			bGlobal = true;
		}
		else if (rVariable == SHELL_CURVATURE)
		{
			ijob = 2;
		}
		else if (rVariable == SHELL_CURVATURE_GLOBAL)
		{
			ijob = 2;
			bGlobal = true;
		}
		else if (rVariable == SHELL_FORCE)
		{
			ijob = 3;
		}
		else if (rVariable == SHELL_FORCE_GLOBAL)
		{
			ijob = 3;
			bGlobal = true;
		}
		else if (rVariable == SHELL_MOMENT)
		{
			ijob = 4;
		}
		else if (rVariable == SHELL_MOMENT_GLOBAL)
		{
			ijob = 4;
			bGlobal = true;
		}

		// quick return

		if (ijob == 0) return false;

		// resize output

		if (rValues.size() != OPT_NUM_GP)
			rValues.resize(OPT_NUM_GP);

		// Just to store the rotation matrix for visualization purposes

		Matrix R(8, 8);
		Matrix aux33(3, 3);

		// Initialize common calculation variables

		CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
		data.CalculateLHS = true;
		data.CalculateRHS = true;
		InitializeCalculationData(data);
		CalculateSectionResponse(data);


		// set the current integration point index
		data.gpIndex = 0;
		ShellCrossSection::Pointer& section = mSections[0];

		// compute strains
		noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

		//compute stresses
		if (ijob >2)
		{
			noalias(data.generalizedStresses) = prod(data.D, data.generalizedStrains);
			DecimalCorrection(data.generalizedStresses);
		}


		// adjust output
		DecimalCorrection(data.generalizedStrains);
		//std::cout << data.generalizedStrains << std::endl;

		// store the results, but first rotate them back to the section coordinate system.
		// we want to visualize the results in that system not in the element one!
		if (section->GetOrientationAngle() != 0.0 && !bGlobal)
		{
			if (ijob > 2)
			{

				section->GetRotationMatrixForGeneralizedStresses(-(section->GetOrientationAngle()), R);
				data.generalizedStresses = prod(R, data.generalizedStresses);
			}
			else
			{
				section->GetRotationMatrixForGeneralizedStrains(-(section->GetOrientationAngle()), R);
				data.generalizedStrains = prod(R, data.generalizedStrains);
			}
		}


		// Gauss Loop.

		for (size_t i = 0; i < OPT_NUM_GP; i++)
		{
			

			// save results
			Matrix & iValue = rValues[i];
			if (iValue.size1() != 3 || iValue.size2() != 3)
				iValue.resize(3, 3, false);

			if (ijob == 1) // strains
			{
				iValue(0, 0) = data.generalizedStrains(0);
				iValue(1, 1) = data.generalizedStrains(1);
				iValue(2, 2) = 0.0;
				iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(2);
				iValue(0, 2) = iValue(2, 0) = 0.5 * data.generalizedStrains(7);
				iValue(1, 2) = iValue(2, 1) = 0.5 * data.generalizedStrains(6);
			}
			else if (ijob == 2) // curvatures
			{
				iValue(0, 0) = data.generalizedStrains(3);
				iValue(1, 1) = data.generalizedStrains(4);
				iValue(2, 2) = 0.0;
				iValue(0, 1) = iValue(1, 0) = 0.5 * data.generalizedStrains(5);
				iValue(0, 2) = iValue(2, 0) = 0.0;
				iValue(1, 2) = iValue(2, 1) = 0.0;
			}
			else if (ijob == 3) // forces
			{
				iValue(0, 0) = data.generalizedStresses(0);
				iValue(1, 1) = data.generalizedStresses(1);
				iValue(2, 2) = 0.0;
				iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(2);
				iValue(0, 2) = iValue(2, 0) = data.generalizedStresses(7);
				iValue(1, 2) = iValue(2, 1) = data.generalizedStresses(6);
			}
			else if (ijob == 4) // moments
			{
				iValue(0, 0) = data.generalizedStresses(3);
				iValue(1, 1) = data.generalizedStresses(4);
				iValue(2, 2) = 0.0;
				iValue(0, 1) = iValue(1, 0) = data.generalizedStresses(5);
				iValue(0, 2) = iValue(2, 0) = 0.0;
				iValue(1, 2) = iValue(2, 1) = 0.0;
			}

			// if requested, rotate the results in the global coordinate system
			if (bGlobal)
			{
				const Matrix& RG = data.LCS.Orientation();
				noalias(aux33) = prod(trans(RG), iValue);
				noalias(iValue) = prod(aux33, RG);
			}
		}

		OPT_INTERPOLATE_RESULTS_TO_STANDARD_GAUSS_POINTS(rValues);
		//Utilities::InterpToStandardGaussPoints(rValues);

		return true;
	}

	void ShellThickElement3D3N::printMatrix(Matrix& matrixIn, std::string stringIn)
	{
		std::cout << "\n" << stringIn << std::endl;
		for (unsigned i = 0; i<matrixIn.size1(); ++i)
		{
			std::cout << "| ";
			for (unsigned j = 0; j<matrixIn.size2(); ++j)
			{
				std::cout << std::fixed << std::setprecision(1) << std::setw(8) << matrixIn(i, j) << " | ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// =====================================================================================
	//
	// Class ShellThickElement3D3N - Serialization
	//
	// =====================================================================================

	void ShellThickElement3D3N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.save("CTr", mpCoordinateTransformation);
		rSerializer.save("Sec", mSections);
		rSerializer.save("IntM", (int)mThisIntegrationMethod);
	}

	void ShellThickElement3D3N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.load("CTr", mpCoordinateTransformation);
		rSerializer.load("Sec", mSections);
		int temp;
		rSerializer.load("IntM", temp);
		mThisIntegrationMethod = (IntegrationMethod)temp;
	}

}

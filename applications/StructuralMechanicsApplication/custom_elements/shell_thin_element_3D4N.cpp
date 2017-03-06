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

#include "shell_thin_element_3D4N.hpp"
#include "custom_utilities/shellq4_corotational_coordinate_transformation.hpp"
#include "structural_mechanics_application_variables.h"

#include "geometries/quadrilateral_3d_4.h"

#include <string>
#include <iomanip>

#define OPT_NUM_NODES 4
#define OPT_STRAIN_SIZE 6
#define OPT_NUM_DOFS 24
#define OPT_NUM_GP 4

namespace Kratos
{

	namespace Utilities
	{
		template<class TVec>
		inline void ShapeFunc(double xi, double eta, TVec & N)
		{
			N(0) = 0.25 * (1.0 - xi) * (1.0 - eta); // node 1
			N(1) = 0.25 * (1.0 + xi) * (1.0 - eta); // node 2
			N(2) = 0.25 * (1.0 + xi) * (1.0 + eta); // node 3
			N(3) = 0.25 * (1.0 - xi) * (1.0 + eta); // node 4
		}

		template<class TMat>
		inline void ShapeFunc_NaturalDerivatives(double xi, double eta, TMat & dN)
		{
			dN(0, 0) = -(1.0 - eta) * 0.25;
			dN(1, 0) = (1.0 - eta) * 0.25;
			dN(2, 0) = (1.0 + eta) * 0.25;
			dN(3, 0) = -(1.0 + eta) * 0.25;

			dN(0, 1) = -(1.0 - xi)  * 0.25;
			dN(1, 1) = -(1.0 + xi)  * 0.25;
			dN(2, 1) = (1.0 + xi)  * 0.25;
			dN(3, 1) = (1.0 - xi)  * 0.25;
		}

		inline double dN_seren_dxi(int actualNodeNumber, double xi, double eta)
		{
			double returnValue;
			switch (actualNodeNumber)
			{
			case 1:
				returnValue = -(-eta + 1)*(-0.25*xi + 0.25) - 0.25*(-eta + 1)*(-eta - xi - 1);
				break;
			case 2:
				returnValue = (-eta + 1)*(0.25*xi + 0.25) + 0.25*(-eta + 1)*(-eta + xi - 1);
				break;
			case 3:
				returnValue = (eta + 1)*(0.25*xi + 0.25) + 0.25*(eta + 1)*(eta + xi - 1);
				break;
			case 4:
				returnValue = -(eta + 1)*(-0.25*xi + 0.25) - 0.25*(eta + 1)*(eta - xi - 1);
				break;
			case 5:
				returnValue = -1.0*xi*(-eta + 1);
				break;
			case 6:
				returnValue = -0.5*eta*eta + 0.5;
				break;
			case 7:
				returnValue = -1.0*xi*(eta + 1);
				break;
			case 8:
				returnValue = 0.5*eta*eta - 0.5;
				break;
			default:
				std::cout << "FALLING INTO DEFAULT LOOP" << std::endl;
				returnValue = 1000;	//for error checking - REMOVE!!!!!
			}

			return returnValue;
		}

		inline double dN_seren_deta(int actualNodeNumber, double xi, double eta)
		{
			double returnValue;
			switch (actualNodeNumber)
			{
			case 1:
				returnValue = -(-eta + 1)*(-0.25*xi + 0.25) - (-0.25*xi + 0.25)*(-eta - xi - 1);
				break;
			case 2:
				returnValue = -(-eta + 1)*(0.25*xi + 0.25) - (0.25*xi + 0.25)*(-eta + xi - 1);
				break;
			case 3:
				returnValue = (eta + 1)*(0.25*xi + 0.25) + (0.25*xi + 0.25)*(eta + xi - 1);
				break;
			case 4:
				returnValue = (eta + 1)*(-0.25*xi + 0.25) + (-0.25*xi + 0.25)*(eta - xi - 1);
				break;
			case 5:
				returnValue = 0.5*xi*xi - 0.5;
				break;
			case 6:
				returnValue = -1.0*eta*(xi + 1);
				break;
			case 7:
				returnValue = -0.5*xi*xi + 0.5;
				break;
			case 8:
				returnValue = -1.0*eta*(-xi + 1);
				break;
			default:
				std::cout << "FALLING INTO DEFAULT LOOP" << std::endl;
				returnValue = 1000;	//for error checking - REMOVE!!!!!	
			}

			return returnValue;
		}

	}


	// =====================================================================================
	//
	// Class JacobianOperator
	//
	// =====================================================================================

	ShellThinElement3D4N::JacobianOperator::JacobianOperator()
		: mJac(2, 2, 0.0)
		, mInv(2, 2, 0.0)
		, mXYDeriv(4, 2, 0.0)
		, mDet(0.0)
	{
	}

	void ShellThinElement3D4N::JacobianOperator::Calculate(const ShellQ4_LocalCoordinateSystem & CS, const Matrix & dN)
	{
		mJac(0, 0) = dN(0, 0) * CS.X1() + dN(1, 0) * CS.X2() + dN(2, 0) * CS.X3() + dN(3, 0) * CS.X4();
		mJac(0, 1) = dN(0, 0) * CS.Y1() + dN(1, 0) * CS.Y2() + dN(2, 0) * CS.Y3() + dN(3, 0) * CS.Y4();
		mJac(1, 0) = dN(0, 1) * CS.X1() + dN(1, 1) * CS.X2() + dN(2, 1) * CS.X3() + dN(3, 1) * CS.X4();
		mJac(1, 1) = dN(0, 1) * CS.Y1() + dN(1, 1) * CS.Y2() + dN(2, 1) * CS.Y3() + dN(3, 1) * CS.Y4();

		mDet = mJac(0, 0) * mJac(1, 1) - mJac(1, 0) * mJac(0, 1);
		double mult = 1.0 / mDet;

		mInv(0, 0) = mJac(1, 1) * mult;
		mInv(0, 1) = -mJac(0, 1) * mult;
		mInv(1, 0) = -mJac(1, 0) * mult;
		mInv(1, 1) = mJac(0, 0) * mult;

		noalias(mXYDeriv) = prod(dN, trans(mInv));
	}


	// =====================================================================================
	//
	// Class ShellThinElement3D4N
	//
	// =====================================================================================

	ShellThinElement3D4N::ShellThinElement3D4N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		bool NLGeom)
		: Element(NewId, pGeometry)
		, mpCoordinateTransformation(NLGeom ?
			new ShellQ4_CorotationalCoordinateTransformation(pGeometry) :
			new ShellQ4_CoordinateTransformation(pGeometry))
	{
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	ShellThinElement3D4N::ShellThinElement3D4N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties,
		bool NLGeom)
		: Element(NewId, pGeometry, pProperties)
		, mpCoordinateTransformation(NLGeom ?
			new ShellQ4_CorotationalCoordinateTransformation(pGeometry) :
			new ShellQ4_CoordinateTransformation(pGeometry))
	{
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	ShellThinElement3D4N::ShellThinElement3D4N(IndexType NewId,
		GeometryType::Pointer pGeometry,
		PropertiesType::Pointer pProperties,
		CoordinateTransformationBasePointerType pCoordinateTransformation)
		: Element(NewId, pGeometry, pProperties)
		, mpCoordinateTransformation(pCoordinateTransformation)
	{
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	ShellThinElement3D4N::~ShellThinElement3D4N()
	{
	}

	//Basic methods

	Element::Pointer ShellThinElement3D4N::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		GeometryType::Pointer newGeom(GetGeometry().Create(ThisNodes));
		return boost::make_shared< ShellThinElement3D4N >(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom));
		//     return Element::Pointer( new ShellThinElement3D4N(NewId, newGeom, pProperties, mpCoordinateTransformation->Create(newGeom)) );
	}

	ShellThinElement3D4N::IntegrationMethod ShellThinElement3D4N::GetIntegrationMethod() const
	{
		return mThisIntegrationMethod;
	}

	void ShellThinElement3D4N::Initialize()
	{
		KRATOS_TRY

			const GeometryType & geom = GetGeometry();
		const PropertiesType & props = GetProperties();

		if (geom.PointsNumber() != OPT_NUM_NODES)
			KRATOS_THROW_ERROR(std::logic_error, "ShellThinElement3D4N Element - Wrong number of nodes", geom.PointsNumber());

		const GeometryType::IntegrationPointsArrayType & integrationPoints = geom.IntegrationPoints(GetIntegrationMethod());
		if (integrationPoints.size() != OPT_NUM_GP)
			KRATOS_THROW_ERROR(std::logic_error, "ShellThinElement3D4N Element - Wrong integration scheme", integrationPoints.size());

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
				//std::cout << "TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
				ShellCrossSection::Pointer sectionClone = theSection->Clone();
				sectionClone->SetSectionBehavior(ShellCrossSection::Thin);
				sectionClone->InitializeCrossSection(props, geom, row(shapeFunctionsValues, i));
				mSections.push_back(sectionClone);
			}
		}

		mpCoordinateTransformation->Initialize();

		this->SetupOrientationAngles();

		KRATOS_CATCH("")
	}

	void ShellThinElement3D4N::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
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

	void ShellThinElement3D4N::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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

	void ShellThinElement3D4N::SetupOrientationAngles()
	{
		ShellQ4_LocalCoordinateSystem lcs(mpCoordinateTransformation->CreateReferenceCoordinateSystem());

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


	void ShellThinElement3D4N::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)
	{
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
	}

	void ShellThinElement3D4N::InitializeCalculationData(CalculationData& data)
	{
		KRATOS_TRY
			//-------------------------------------
			// Computation of all stuff that remain
			// constant throughout the calculations

			//-------------------------------------
			// geometry data
			//std::cout << "Printing x's \n " << data.LCS0.X1() << " , " << data.LCS0.X2() << " , " << data.LCS0.X3() << " , " << data.LCS0.X4() << std::endl;
			//std::cout << "Printing y's \n " << data.LCS0.Y1() << " , " << data.LCS0.Y2() << " , " << data.LCS0.Y3() << " , " << data.LCS0.Y4() << std::endl;
			//std::cout << "Printing z's \n " << data.LCS0.Z1() << " , " << data.LCS0.Z2() << " , " << data.LCS0.Z3() << " , " << data.LCS0.Z4() << std::endl;

		const double x12 = data.LCS0.X1() - data.LCS0.X2();
		const double x13 = data.LCS0.X1() - data.LCS0.X3();
		const double x23 = data.LCS0.X2() - data.LCS0.X3();
		const double x24 = data.LCS0.X2() - data.LCS0.X4();
		const double x34 = data.LCS0.X3() - data.LCS0.X4();
		const double x41 = data.LCS0.X4() - data.LCS0.X1();

		const double x21 = -x12;
		const double x31 = -x13;
		const double x32 = -x23;
		const double x42 = -x24;
		const double x43 = -x34;
		const double x14 = -x41;

		const double y12 = data.LCS0.Y1() - data.LCS0.Y2();
		const double y13 = data.LCS0.Y1() - data.LCS0.Y3();
		const double y23 = data.LCS0.Y2() - data.LCS0.Y3();
		const double y24 = data.LCS0.Y2() - data.LCS0.Y4();
		const double y34 = data.LCS0.Y3() - data.LCS0.Y4();
		const double y41 = data.LCS0.Y4() - data.LCS0.Y1();

		const double y21 = -y12;
		const double y31 = -y13;
		const double y32 = -y23;
		const double y42 = -y24;
		const double y43 = -y34;
		const double y14 = -y41;

		for (int i = 0; i < 4; i++)
		{
			data.r_cartesian[i] = Vector(3, 0.0);
		}
		data.r_cartesian[0] = data.LCS0.P1();
		data.r_cartesian[1] = data.LCS0.P2();
		data.r_cartesian[2] = data.LCS0.P3();
		data.r_cartesian[3] = data.LCS0.P4();

		for (int i = 0; i < 4; i++)
		{
			//std::cout << "Printing Point " << i << ": " << data.r_cartesian[i] << std::endl;
		}

		//Precalculate dA to be multiplied with material matrix
		GeometryType & geom = GetGeometry();
		const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(mThisIntegrationMethod);
		data.dA.clear();
		for (int gp = 0; gp < 4; gp++)
		{
			//getting informations for integration
			double IntegrationWeight = integration_points[gp].Weight();

			// Compute Jacobian, Inverse of Jacobian, Determinant of Jacobian
			// and Shape functions derivatives in the local coordinate system
			data.jacOp.Calculate(data.LCS0, geom.ShapeFunctionLocalGradient(gp));

			// compute the 'area' of the current integration point
			data.dA[gp] = IntegrationWeight * data.jacOp.Determinant();

		}






		// Note: here we compute the avarage thickness,
		// since L is constant over the element.
		// Now it is not necessary to compute the avarage
		// because the current implementation of the cross section
		// doesn't have a variable thickness
		// (for example as a function of the spatial coordinates...).
		// This is just a place-holder for future
		// implementation of a variable thickness

		const double A = data.LCS0.Area();
		//printDouble("Area = ", A);
		double h = 0.0;
		for (unsigned int i = 0; i < mSections.size(); i++)
			h += mSections[i]->GetThickness();
		h /= (double)mSections.size();
		//data.hMean = h;
		//data.TotalArea = A;
		//data.TotalVolume = A * h;


		// Unit vectors s_xi and s_eta (eqn 5.2.25)
		// eqns 5.2.29 -> 5.2.31
		data.s_xi.resize(3, false);
		data.s_xi.clear();
		data.s_eta.resize(3, false);
		data.s_eta.clear();

		//set values of SFs to xi = 1 and eta = 0
		Utilities::ShapeFunc(1, 0, data.N_pw);
		for (int i = 0; i < 4; i++)
		{
			data.s_xi(0) += data.r_cartesian[i][0] * data.N_pw(i);
			data.s_xi(1) += data.r_cartesian[i][1] * data.N_pw(i);
		}
		double s_xi_mag = std::sqrt(inner_prod(data.s_xi, data.s_xi));
		data.s_xi = data.s_xi / s_xi_mag;
		//set values of SFs to xi = 0 and eta = 1
		Utilities::ShapeFunc(0, 1, data.N_pw);
		for (int i = 0; i < 4; i++)
		{
			data.s_eta(0) += data.r_cartesian[i][0] * data.N_pw(i);
			data.s_eta(1) += data.r_cartesian[i][1] * data.N_pw(i);
		}
		double s_eta_mag = std::sqrt(inner_prod(data.s_eta, data.s_eta));
		data.s_eta = data.s_eta / s_eta_mag;


		//Template constants
		const double alpha_6 = data.alpha / 6.0;
		const double alpha_3 = data.alpha / 3.0;

		// calculate L - Lumping matrix (ref eqn 5.2.4)
		// for the construction of the basic
		// stiffness. Transpose from the presented version to allow combination of B matrices later

		double L_mult = 0.5 / A;
		data.L.resize(3, 12, false);
		data.L.clear();

		//j = 1, i=4, k=2
		//ki = 24, ij = 41 
		data.L(0, 0) = L_mult * y24;
		data.L(1, 0) = 0.00;
		data.L(2, 0) = L_mult * x42;
		data.L(0, 1) = 0.00;
		data.L(1, 1) = L_mult * x42;
		data.L(2, 1) = L_mult * y24;
		data.L(0, 2) = L_mult * alpha_6*(y41*y41 - y21*y21);
		data.L(1, 2) = L_mult *  alpha_6*(x41*x41 - x21*x21);
		data.L(2, 2) = L_mult * alpha_3*(x21*y21 - x41*y41);

		//j = 2, i=1, k=3
		//ki = 31, ij = 12 
		data.L(0, 3) = L_mult * y31;
		data.L(1, 3) = 0.00;
		data.L(2, 3) = L_mult * x13;
		data.L(0, 4) = 0.00;
		data.L(1, 4) = L_mult * x13;
		data.L(2, 4) = L_mult * y31;
		data.L(0, 5) = L_mult * alpha_6*(y12*y12 - y32*y32);
		data.L(1, 5) = L_mult * alpha_6*(x12*x12 - x32*x32);
		data.L(2, 5) = L_mult * alpha_3*(x32*y32 - x12*y12);

		//j = 3, i=2, k=4
		//ki = 42, ij = 23 
		data.L(0, 6) = L_mult * y42;
		data.L(1, 6) = 0.00;
		data.L(2, 6) = L_mult * x24;
		data.L(0, 7) = 0.00;
		data.L(1, 7) = L_mult * x24;
		data.L(2, 7) = L_mult * y42;
		data.L(0, 8) = L_mult * alpha_6*(y23*y23 - y43*y43);
		data.L(1, 8) = L_mult * alpha_6*(x23*x23 - x43*x43);
		data.L(2, 8) = L_mult * alpha_3*(x43*y43 - x23*y23);

		//j = 4, i=3, k=1
		//ki = 13, ij = 34 
		data.L(0, 9) = L_mult * y13;
		data.L(1, 9) = 0.00;
		data.L(2, 9) = L_mult * x31;
		data.L(0, 10) = 0.00;
		data.L(1, 10) = L_mult * x31;
		data.L(2, 10) = L_mult * y13;
		data.L(0, 11) = L_mult * alpha_6*(y34*y34 - y14*y14);
		data.L(1, 11) = L_mult * alpha_6*(x34*x34 - x14*x14);
		data.L(2, 11) = L_mult * alpha_3*(x14*y14 - x34*y34);


		//--------------------------------------
		// calculate H - matrix
		// for the construction of the
		// higher order stiffness
		// modified for rearranged DOF ordering!

		//H_mem transformation matrix 'H' (ref eqn 5.2.26)
		data.H_mem.resize(7, 12, false);
		data.H_mem.clear();

		//Section 3.1.1 - H_theta_v transformation matrix (ref eqn 5.2.14)
		Matrix H_theta_v = Matrix(5, 12, 0.0);

		//ref eqn 5.2.12
		double J_mag = 0;
		J_mag += (data.LCS0.X1()*data.LCS0.Y2() - data.LCS0.X2()*data.LCS0.Y1());
		J_mag += (data.LCS0.X2()*data.LCS0.Y3() - data.LCS0.X3()*data.LCS0.Y2());
		J_mag += (data.LCS0.X3()*data.LCS0.Y4() - data.LCS0.X4()*data.LCS0.Y3());
		J_mag += (data.LCS0.X4()*data.LCS0.Y1() - data.LCS0.X1()*data.LCS0.Y4());
		J_mag /= 8;
		double f_coefficient = J_mag * 16;

		//ref eqn 5.2.14 - modified for rearranged DOF ordering!
		Matrix H_theta_v_mod = Matrix(5, 12, 0.0);
		H_theta_v_mod(0, 2) = 0.75;
		H_theta_v_mod(0, 5) = -0.25;
		H_theta_v_mod(0, 8) = -0.25;
		H_theta_v_mod(0, 11) = -0.25;

		H_theta_v_mod(1, 2) = -0.25;
		H_theta_v_mod(1, 5) = 0.75;
		H_theta_v_mod(1, 8) = -0.25;
		H_theta_v_mod(1, 11) = -0.25;

		H_theta_v_mod(2, 2) = -0.25;
		H_theta_v_mod(2, 5) = -0.25;
		H_theta_v_mod(2, 8) = 0.75;
		H_theta_v_mod(2, 11) = -0.25;

		H_theta_v_mod(3, 2) = -0.25;
		H_theta_v_mod(3, 5) = -0.25;
		H_theta_v_mod(3, 8) = -0.25;
		H_theta_v_mod(3, 11) = 0.75;

		H_theta_v_mod(4, 0) = x42 / f_coefficient;
		H_theta_v_mod(4, 1) = y42 / f_coefficient;
		H_theta_v_mod(4, 2) = 0.25;
		H_theta_v_mod(4, 3) = x13 / f_coefficient;
		H_theta_v_mod(4, 4) = y13 / f_coefficient;
		H_theta_v_mod(4, 5) = 0.25;
		H_theta_v_mod(4, 6) = x24 / f_coefficient;
		H_theta_v_mod(4, 7) = y24 / f_coefficient;
		H_theta_v_mod(4, 8) = 0.25;
		H_theta_v_mod(4, 9) = x31 / f_coefficient;
		H_theta_v_mod(4, 10) = y31 / f_coefficient;
		H_theta_v_mod(4, 11) = 0.25;

		//vector v_h - ref eqn 5.2.21
		Vector v_h = Vector(4, 0.0);
		double variable_a = 0.0;
		v_h[0] = 1 * 1 - variable_a;
		v_h[1] = -1 * 1 - variable_a;
		v_h[2] = -1 * -1 - variable_a;
		v_h[3] = 1 * -1 - variable_a;

		//H_v_t transformation matrix (ref eqn 5.2.22, 5.2.23) - modified for rearranged DOF ordering!
		Matrix H_v_t = Matrix(2, 12, 0.0);

		Matrix H_v_t_mod = Matrix(2, 12, 0.0);
		H_v_t_mod(0, 0) = v_h[0];
		H_v_t_mod(1, 1) = v_h[0];

		H_v_t_mod(0, 3) = v_h[1];
		H_v_t_mod(1, 4) = v_h[1];

		H_v_t_mod(0, 6) = v_h[2];
		H_v_t_mod(1, 7) = v_h[2];

		H_v_t_mod(0, 9) = v_h[3];
		H_v_t_mod(1, 10) = v_h[3];

		//Combine into H mat - eqn 5.2.26
		data.H_mem_mod.resize(7, 12, false);
		data.H_mem_mod.clear();

		for (int col_index = 0; col_index < 12; col_index++)
		{
			for (int row_index = 0; row_index < 5; row_index++)
			{
				data.H_mem_mod(row_index, col_index) = H_theta_v_mod(row_index, col_index);
			}

			for (int row_index = 0; row_index < 2; row_index++)
			{
				data.H_mem_mod(row_index + 5, col_index) = H_v_t_mod(row_index, col_index);
			}
		}

		//--------------------------------------
		// calculate Q1,Q2,Q3,Q4 - matrices
		// for the construction of the
		// higher order stiffness

		CalculateQMatrices(data);


		//--------------------------------------
		// calculate T_13, T_24 -
		// transformation matrices
		// for the construction of the
		// higher order stiffness
		// ref eqn 5.2.36 and 5.2.37
		Matrix T_13_inv = Matrix(3, 3, 0.0);
		double T_13_inv_det = 0;
		array_1d<Vector, 3> vecList;
		vecList[0] = (data.s_xi);
		vecList[1] = (data.s_eta);
		vecList[2] = (data.s_24);

		for (int i = 0; i < 3; i++)
		{
			T_13_inv(i, 0) += vecList[i][0] * vecList[i][0];
			T_13_inv(i, 1) += vecList[i][1] * vecList[i][1];
			T_13_inv(i, 2) += vecList[i][0] * vecList[i][1];
		}
		T_13_inv_det = MathUtils<double>::Det(T_13_inv);
		MathUtils<double>::InvertMatrix(T_13_inv, data.T_13, T_13_inv_det);

		Matrix T_24_inv = Matrix(T_13_inv);
		double T_24_inv_det = 0;
		vecList[2] = data.s_13;
		for (int i = 2; i < 3; i++)
		{
			T_24_inv(i, 0) = vecList[i][0] * vecList[i][0];
			T_24_inv(i, 1) = vecList[i][1] * vecList[i][1];
			T_24_inv(i, 2) = vecList[i][0] * vecList[i][1];
		}
		T_24_inv_det = MathUtils<double>::Det(T_24_inv);
		MathUtils<double>::InvertMatrix(T_24_inv, data.T_24, T_24_inv_det);


		//--------------------------------------
		// calculate B_h_bar
		// for the construction of the
		// higher order MEMBRANE stiffness

		CalculateB_h_mats(data);













		//--------------------------------------
		//
		//     BENDING PART OF THE ELEMENT
		//
		//--------------------------------------

		//Calculate edge normal vectors for eqn 5.3.15
		//ref eqn 5.1.9 for calc of normal vector from edge vec
		//std::cout << "\nCurrently the z component of edge normals are set to 0!!" << std::endl;
		Vector s_12 = Vector(data.LCS0.P1() - data.LCS0.P2());
		const double l_12 = std::sqrt(inner_prod(s_12, s_12));

		Vector s_23 = Vector(data.LCS0.P2() - data.LCS0.P3());
		const double l_23 = std::sqrt(inner_prod(s_23, s_23));

		Vector s_34 = Vector(data.LCS0.P3() - data.LCS0.P4());
		const double l_34 = std::sqrt(inner_prod(s_34, s_34));

		Vector s_41 = Vector(data.LCS0.P4() - data.LCS0.P1());
		const double l_41 = std::sqrt(inner_prod(s_41, s_41));

		Vector s_13 = Vector(data.LCS0.P1() - data.LCS0.P3());
		const double l_13 = std::sqrt(inner_prod(s_13, s_13));

		Vector s_24 = Vector(data.LCS0.P2() - data.LCS0.P4());
		const double l_24 = std::sqrt(inner_prod(s_24, s_24));





		//--------------------------------------
		// calculate DKQ bending stiffness

		bool DKQ = true;
		if (DKQ == true)
		{
			data.DKQ_a.clear();
			data.DKQ_b.clear();
			data.DKQ_c.clear();
			data.DKQ_d.clear();
			data.DKQ_e.clear();

			//assemble a_k - eqn 3.86a : a[0] = a_5
			data.DKQ_a[0] = -1 * x12 / l_12 / l_12;
			data.DKQ_a[1] = -1 * x23 / l_23 / l_23;
			data.DKQ_a[2] = -1 * x34 / l_34 / l_34;
			data.DKQ_a[3] = -1 * x41 / l_41 / l_41;

			//assemble b_k - eqn 3.86b : b[0] = b_5
			data.DKQ_b[0] = 3 / 4 * x12 * y12 / l_12 / l_12;
			data.DKQ_b[1] = 3 / 4 * x23 * y23 / l_23 / l_23;
			data.DKQ_b[2] = 3 / 4 * x34 * y34 / l_34 / l_34;
			data.DKQ_b[3] = 3 / 4 * x41 * y41 / l_41 / l_41;

			//assemble c_k - eqn 3.86c : c[0] = c_5
			data.DKQ_c[0] = (x12 * x12 / 4 - y12 * y12 / 2) / l_12 / l_12;
			data.DKQ_c[1] = (x23 * x23 / 4 - y23 * y23 / 2) / l_23 / l_23;
			data.DKQ_c[2] = (x34 * x34 / 4 - y34 * y34 / 2) / l_34 / l_34;
			data.DKQ_c[3] = (x41 * x41 / 4 - y41 * y41 / 2) / l_41 / l_41;

			//assemble d_k - eqn 3.86d : d[0] = d_5
			data.DKQ_d[0] = -1 * y12 / l_12 / l_12;
			data.DKQ_d[1] = -1 * y23 / l_23 / l_23;
			data.DKQ_d[2] = -1 * y34 / l_34 / l_34;
			data.DKQ_d[3] = -1 * y41 / l_41 / l_41;

			//assemble e_k - eqn 3.86e : e[0] = e_5
			data.DKQ_e[0] = (-1 * x12 * x12 / 2 + y12 * y12 / 4) / l_12 / l_12;
			data.DKQ_e[1] = (-1 * x23 * x23 / 2 + y23 * y23 / 4) / l_23 / l_23;
			data.DKQ_e[2] = (-1 * x34 * x34 / 2 + y34 * y34 / 4) / l_34 / l_34;
			data.DKQ_e[3] = (-1 * x41 * x41 / 2 + y41 * y41 / 4) / l_41 / l_41;

			//prepare DKT assembly indices - for eqn 3.92 a->f
			data.DKQ_indices.resize(4, 2);
			data.DKQ_indices.clear();
			// 1st col = r, 2nd col = s
			data.DKQ_indices(0, 0) = 5;	//actual node number, not index!
			data.DKQ_indices(0, 1) = 8;
			data.DKQ_indices(1, 0) = 6;
			data.DKQ_indices(1, 1) = 5;
			data.DKQ_indices(2, 0) = 7;
			data.DKQ_indices(2, 1) = 6;
			data.DKQ_indices(3, 0) = 8;
			data.DKQ_indices(3, 1) = 7;

			//custom jacobian as per eqn 14
			GeometryType & geom = GetGeometry();
			const GeometryType::IntegrationPointsArrayType & integrationPoints = geom.IntegrationPoints(GetIntegrationMethod());
			for (int i = 0; i < 4; i++)
			{
				double xi = integrationPoints[i].Coordinate(1);		//set to current parametric integration location
				double eta = integrationPoints[i].Coordinate(2);		//set to current parametric integration location
				Matrix DKQ_temp = Matrix(2, 2, 0.0);
				DKQ_temp(0, 0) = x21 + x34 + eta*(x12 + x34);
				DKQ_temp(0, 1) = y21 + y34 + eta*(y12 + y34);
				DKQ_temp(1, 0) = x32 + x41 + xi*(x12 + x34);
				DKQ_temp(1, 1) = y32 + y41 + xi*(y12 + y34);
				DKQ_temp = DKQ_temp / 4;
				//printMatrix(DKQ_temp, "Printing DKQ_temp");
				double det = MathUtils<double>::Det(DKQ_temp);
				Matrix DKQ_temp_inv = Matrix(2, 2, 0.0);
				DKQ_temp_inv(0, 0) = DKQ_temp(1, 1);
				DKQ_temp_inv(1, 1) = DKQ_temp(0, 0);
				DKQ_temp_inv(0, 1) = -1 * DKQ_temp(0, 1);
				DKQ_temp_inv(1, 0) = -1 * DKQ_temp(1, 0);
				DKQ_temp_inv = DKQ_temp_inv / det;
				data.DKQ_invJac[i] = Matrix(DKQ_temp_inv);
			}

		}






		//--------------------------------------
		// calculate the displacement vector
		// in global and local coordinate systems

		data.globalDisplacements.resize(OPT_NUM_DOFS, false);
		GetValuesVector(data.globalDisplacements);

		data.localDisplacements =
			mpCoordinateTransformation->CalculateLocalDisplacements(
				data.LCS, data.globalDisplacements);


		//--------------------------------------
		// Finally allocate all auxiliary
		// matrices to be used later on
		// during the element integration.
		// Just to avoid re-allocations

		data.B.resize(OPT_STRAIN_SIZE, OPT_NUM_DOFS, false);
		data.D.resize(OPT_STRAIN_SIZE, OPT_STRAIN_SIZE, false);
		data.BTD.resize(OPT_NUM_DOFS, OPT_STRAIN_SIZE, false);

		data.generalizedStrains.resize(OPT_STRAIN_SIZE, false);
		data.generalizedStresses.resize(OPT_STRAIN_SIZE, false);

		//--------------------------------------
		// Initialize the section parameters

		data.SectionParameters.SetElementGeometry(GetGeometry());
		data.SectionParameters.SetMaterialProperties(GetProperties());
		data.SectionParameters.SetProcessInfo(data.CurrentProcessInfo);

		data.SectionParameters.SetGeneralizedStrainVector(data.generalizedStrains);
		data.SectionParameters.SetGeneralizedStressVector(data.generalizedStresses);
		data.SectionParameters.SetConstitutiveMatrix(data.D);

		//this isn't constant for a quad!!!!
		//data.SectionParameters.SetShapeFunctionsDerivatives(data.dNxy);

		Flags& options = data.SectionParameters.GetOptions();
		options.Set(ConstitutiveLaw::COMPUTE_STRESS, data.CalculateRHS);
		options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, data.CalculateLHS);

		KRATOS_CATCH("")
	}

	void ShellThinElement3D4N::CalculateQMatrices(CalculationData& data)
	{
		//Calculate Q Matrices according to eqns 5.2.32 -> 5.2.35

		// Template constants defined in eqn 5.2.41
		double rho_1 = 0.1;
		double rho_2 = -0.1;
		double rho_3 = -0.1;
		double rho_4 = 0.1;
		double rho_5 = 0;
		double rho_6 = 0.5;
		double rho_7 = 0;
		double rho_8 = -0.5;
		double beta_1 = 0.6;
		double beta_2 = 0;

		// Preparing entries for Q Matrices
		// eqns 5.2.29 -> 5.2.31

		//5.2.29
		Vector r_xi = Vector(0.5*(data.r_cartesian[1] + data.r_cartesian[2] - data.r_cartesian[0] - data.r_cartesian[3]));
		Vector r_eta = Vector(0.5*(data.r_cartesian[2] + data.r_cartesian[3] - data.r_cartesian[0] - data.r_cartesian[1]));
		double l_xi = std::sqrt(inner_prod(r_xi, r_xi));
		double l_eta = std::sqrt(inner_prod(r_eta, r_eta));

		array_1d<double, 4> d_xi_i;
		for (int i = 0; i < 4; i++)
		{
			d_xi_i[i] = std::sqrt(inner_prod((MathUtils<double>::CrossProduct(data.r_cartesian[i], data.s_xi)), (MathUtils<double>::CrossProduct(data.r_cartesian[i], data.s_xi))));
		}

		array_1d<double, 4> d_eta_i;
		for (int i = 0; i < 4; i++)
		{
			d_eta_i[i] = std::sqrt(inner_prod((MathUtils<double>::CrossProduct(data.r_cartesian[i], data.s_eta)), (MathUtils<double>::CrossProduct(data.r_cartesian[i], data.s_eta))));
		}


		Vector chi_xi_i = Vector(4, 0.0);
		for (int i = 0; i < 4; i++)
		{
			chi_xi_i[i] = d_xi_i[i] / l_xi;
		}

		Vector chi_eta_i = Vector(4, 0.0);
		for (int i = 0; i < 4; i++)
		{
			chi_eta_i[i] = d_eta_i[i] / l_eta;
		}


		//5.2.30
		Vector r_24 = Vector(data.LCS0.P2() - data.LCS0.P4());
		double l_24 = std::sqrt(inner_prod(r_24, r_24));
		Vector e_24 = Vector(r_24 / l_24);

		Vector r_13 = Vector(data.LCS0.P1() - data.LCS0.P3());
		double l_13 = std::sqrt(inner_prod(r_13, r_13));
		Vector e_13 = Vector(r_13 / l_13);

		double d_24 = std::sqrt(inner_prod(MathUtils<double>::CrossProduct(r_13*-1, e_24), MathUtils<double>::CrossProduct(r_13, e_24)));
		double d_13 = std::sqrt(inner_prod(MathUtils<double>::CrossProduct(r_13*-1, e_24), MathUtils<double>::CrossProduct(r_13, e_24)));

		double chi_24 = d_24 / 2 / l_24;
		double chi_13 = d_13 / 2 / l_13;

		//5.2.31
		double chi_xi_t = l_eta / l_xi;
		double chi_eta_t = l_xi / l_eta;

		//calc chi_bars, in 5.2.32 -> 5.2.35
		// average of all chi_xi_i or chi_eta_i
		double chi_xi_bar = 0;
		for (int i = 0; i < 4; i++)
		{
			chi_xi_bar += chi_xi_i[i] / 4;
		}

		double chi_eta_bar = 0;
		for (int i = 0; i < 4; i++)
		{
			chi_eta_bar += chi_eta_i[i] / 4;
		}

		data.s_13 = Vector(r_13 / l_13);
		data.s_24 = Vector(r_24 / l_24);


		double c_13_xi = MathUtils<double>::Dot(data.s_13, data.s_xi);
		double c_13_eta = MathUtils<double>::Dot(data.s_13, data.s_eta);
		double c_24_xi = MathUtils<double>::Dot(data.s_24, data.s_xi);
		double c_24_eta = MathUtils<double>::Dot(data.s_24, data.s_eta);



		data.Q1.resize(3, 7, false);
		data.Q1.clear();
		data.Q2.resize(3, 7, false);
		data.Q2.clear();
		data.Q3.resize(3, 7, false);
		data.Q3.clear();
		data.Q4.resize(3, 7, false);
		data.Q4.clear();

		//Preparing Q1
		int Q_index = 1;

		//row 1
		data.Q1(0, 0) = rho_1 * chi_xi_i[Q_index - 1];
		data.Q1(0, 1) = rho_2 * chi_xi_i[Q_index - 1];
		data.Q1(0, 2) = rho_3 * chi_xi_i[Q_index - 1];
		data.Q1(0, 3) = rho_4 * chi_xi_i[Q_index - 1];

		data.Q1(0, 4) = data.alpha * chi_xi_t;

		data.Q1(0, 5) = -1 * beta_1 * chi_xi_i[Q_index - 1] / chi_xi_bar / l_xi;

		data.Q1(0, 6) = 0;

		//row 2
		data.Q1(1, 0) = -1 * rho_1 * chi_eta_i[Q_index - 1];
		data.Q1(1, 1) = -1 * rho_4 * chi_eta_i[Q_index - 1];
		data.Q1(1, 2) = -1 * rho_3 * chi_eta_i[Q_index - 1];
		data.Q1(1, 3) = -1 * rho_2 * chi_eta_i[Q_index - 1];

		data.Q1(1, 4) = -1 * data.alpha * chi_eta_t;

		data.Q1(1, 5) = 0;

		data.Q1(1, 6) = -1 * beta_1 * chi_eta_i[Q_index - 1] / chi_eta_bar / l_eta;

		//row 3
		data.Q1(2, 0) = rho_5 * chi_24;
		data.Q1(2, 1) = rho_6 * chi_24;
		data.Q1(2, 2) = rho_7 * chi_24;
		data.Q1(2, 3) = rho_8 * chi_24;

		data.Q1(2, 4) = 0;

		data.Q1(2, 5) = beta_2 * c_24_xi / l_24;

		data.Q1(2, 6) = -1 * beta_2 * c_24_eta / l_24;


		//Preparing Q2
		Q_index = 2;

		//row 1
		data.Q2(0, 0) = -1 * rho_2 * chi_xi_i[Q_index - 1];
		data.Q2(0, 1) = -1 * rho_1 * chi_xi_i[Q_index - 1];
		data.Q2(0, 2) = -1 * rho_4 * chi_xi_i[Q_index - 1];
		data.Q2(0, 3) = -1 * rho_3 * chi_xi_i[Q_index - 1];

		data.Q2(0, 4) = -1 * data.alpha * chi_xi_t;

		data.Q2(0, 5) = -1 * beta_1 * chi_xi_i[Q_index - 1] / chi_xi_bar / l_xi;

		data.Q2(0, 6) = 0;

		//row 2
		data.Q2(1, 0) = rho_4 * chi_eta_i[Q_index - 1];
		data.Q2(1, 1) = rho_1 * chi_eta_i[Q_index - 1];
		data.Q2(1, 2) = rho_2 * chi_eta_i[Q_index - 1];
		data.Q2(1, 3) = rho_3 * chi_eta_i[Q_index - 1];

		data.Q2(1, 4) = data.alpha * chi_eta_t;

		data.Q2(1, 5) = 0;

		data.Q2(1, 6) = beta_1 * chi_eta_i[Q_index - 1] / chi_eta_bar / l_eta;

		//row 3
		data.Q2(2, 0) = rho_8 * chi_13;
		data.Q2(2, 1) = rho_5 * chi_13;
		data.Q2(2, 2) = rho_6 * chi_13;
		data.Q2(2, 3) = rho_7 * chi_13;

		data.Q2(2, 4) = 0;

		data.Q2(2, 5) = -1 * beta_2 * c_13_xi / l_13;

		data.Q2(2, 6) = beta_2 * c_13_eta / l_13;

		//printMatrix(Q[1], "Q[1]");

		//TEMPLATE 3
		Q_index = 3;

		//row 1
		data.Q3(0, 0) = rho_3 * chi_xi_i[Q_index - 1];
		data.Q3(0, 1) = rho_4 * chi_xi_i[Q_index - 1];
		data.Q3(0, 2) = rho_1 * chi_xi_i[Q_index - 1];
		data.Q3(0, 3) = rho_2 * chi_xi_i[Q_index - 1];

		data.Q3(0, 4) = data.alpha * chi_xi_t;

		data.Q3(0, 5) = beta_1 * chi_xi_i[Q_index - 1] / chi_xi_bar / l_xi;

		data.Q3(0, 6) = 0;

		//row 2
		data.Q3(1, 0) = -1 * rho_3 * chi_eta_i[Q_index - 1];
		data.Q3(1, 1) = -1 * rho_2 * chi_eta_i[Q_index - 1];
		data.Q3(1, 2) = -1 * rho_1 * chi_eta_i[Q_index - 1];
		data.Q3(1, 3) = -1 * rho_4 * chi_eta_i[Q_index - 1];

		data.Q3(1, 4) = -1 * data.alpha * chi_eta_t;

		data.Q3(1, 5) = 0;

		data.Q3(1, 6) = beta_1 * chi_eta_i[Q_index - 1] / chi_eta_bar / l_eta;

		//row 3
		data.Q3(2, 0) = rho_7 * chi_13;
		data.Q3(2, 1) = rho_8 * chi_13;
		data.Q3(2, 2) = rho_5 * chi_13;
		data.Q3(2, 3) = rho_6 * chi_13;

		data.Q3(2, 4) = 0;

		data.Q3(2, 5) = -1 * beta_2 * c_13_xi / l_13;

		data.Q3(2, 6) = beta_2 * c_13_eta / l_13;

		//TEMPLATE 4
		Q_index = 4;

		//row 1
		data.Q4(0, 0) = -1 * rho_4 * chi_xi_i[Q_index - 1];
		data.Q4(0, 1) = -1 * rho_3 * chi_xi_i[Q_index - 1];
		data.Q4(0, 2) = -1 * rho_2 * chi_xi_i[Q_index - 1];
		data.Q4(0, 3) = -1 * rho_1 * chi_xi_i[Q_index - 1];

		data.Q4(0, 4) = -1 * data.alpha * chi_xi_t;

		data.Q4(0, 5) = beta_1 * chi_xi_i[Q_index - 1] / chi_xi_bar / l_xi;

		data.Q4(0, 6) = 0;

		//row 2
		data.Q4(1, 0) = rho_2 * chi_eta_i[Q_index - 1];
		data.Q4(1, 1) = rho_3 * chi_eta_i[Q_index - 1];
		data.Q4(1, 2) = rho_4 * chi_eta_i[Q_index - 1];
		data.Q4(1, 3) = rho_1 * chi_eta_i[Q_index - 1];

		data.Q4(1, 4) = data.alpha * chi_eta_t;

		data.Q4(1, 5) = 0;

		data.Q4(1, 6) = -1 * beta_1 * chi_eta_i[Q_index - 1] / chi_eta_bar / l_eta;

		//row 3
		data.Q4(2, 0) = rho_6 * chi_13;
		data.Q4(2, 1) = rho_7 * chi_13;
		data.Q4(2, 2) = rho_8 * chi_13;
		data.Q4(2, 3) = rho_5 * chi_13;

		data.Q4(2, 4) = 0;

		data.Q4(2, 5) = beta_2 * c_13_xi / l_13;

		data.Q4(2, 6) = -1 * beta_2 * c_13_eta / l_13;
	}

	void ShellThinElement3D4N::CalculateB_h_mats(CalculationData& data) // FOR MEMBRANE PART
	{
		data.B_h_1.resize(3, 7, false);
		data.B_h_1.clear();
		data.B_h_2.resize(3, 7, false);
		data.B_h_2.clear();
		data.B_h_3.resize(3, 7, false);
		data.B_h_3.clear();
		data.B_h_4.resize(3, 7, false);
		data.B_h_4.clear();
		noalias(data.B_h_1) = prod(data.T_13, data.Q1);
		noalias(data.B_h_2) = prod(data.T_24, data.Q2);
		noalias(data.B_h_3) = prod(data.T_13, data.Q3);
		noalias(data.B_h_4) = prod(data.T_24, data.Q4);
	}

	void ShellThinElement3D4N::CalculateBMatrix(CalculationData& data)
	{
		//---------------------------------------------
		// geom data
		GeometryType & geom = GetGeometry();

		//---------------------------------------------
		// set to zero the total B matrix
		data.B.clear();


		//---------------------------------------------
		// membrane basic part L
		// already computed.it is constant over the
		// element

		//---------------------------------------------
		// membrane higher order part B_h
		// already calculated in calculate B_h_mats func
		Matrix B_h = Matrix(3, 7, 0.0);
		const Matrix & shapeFunctionsValues = geom.ShapeFunctionsValues(GetIntegrationMethod());
		B_h += shapeFunctionsValues(data.gpIndex, 0) * data.B_h_1;
		B_h += shapeFunctionsValues(data.gpIndex, 1) * data.B_h_2;
		B_h += shapeFunctionsValues(data.gpIndex, 2) * data.B_h_3;
		B_h += shapeFunctionsValues(data.gpIndex, 3) * data.B_h_4;


		//---------------------------------------------
		// combine membrane entries by transforming 
		// B_h
		Matrix B_mem = Matrix(3, 12, 0.0);
		Matrix B_hH_mem = Matrix(prod(B_h, data.H_mem_mod));
		B_mem += data.L + B_hH_mem;








		// ---------------------- DKQ BENDING FORMULATION HERE -----------------------

		//---------------------------------------------
		Matrix B_bend_total = Matrix(3, 12, 0.0);


		Vector dpsiX_dxi = Vector(12, 0.0);
		Vector dpsiX_deta = Vector(12, 0.0);
		Vector dpsiY_dxi = Vector(12, 0.0);
		Vector dpsiY_deta = Vector(12, 0.0);
		const GeometryType::IntegrationPointsArrayType & integrationPoints = geom.IntegrationPoints(GetIntegrationMethod());
		double xi = integrationPoints[data.gpIndex].Coordinate(1);		//set to current parametric integration location
		double eta = integrationPoints[data.gpIndex].Coordinate(2);		//set to current parametric integration location

		double ar, br, cr, dr, er;
		double as, bs, cs, ds, es;
		int r, s;
		for (int node = 0; node < 4; node++)
		{
			r = data.DKQ_indices(node, 0);	//actual node number r retrieved
			s = data.DKQ_indices(node, 1);	//actual node number s retrieved

			ar = data.DKQ_a[r - 5];	//converted to index, where node 5 = index 0
			br = data.DKQ_b[r - 5];
			cr = data.DKQ_c[r - 5];
			dr = data.DKQ_d[r - 5];
			er = data.DKQ_e[r - 5];

			as = data.DKQ_a[s - 5];
			bs = data.DKQ_b[s - 5];
			cs = data.DKQ_c[s - 5];
			ds = data.DKQ_d[s - 5];
			es = data.DKQ_e[s - 5];

			//eqn 3.92 a->f
			// d( ) / dxi
			// Compute vector of dPsi_x/dxi
			dpsiX_dxi[3 * node] = 1.5 * (ar*Utilities::dN_seren_dxi(r, xi, eta) - as*Utilities::dN_seren_dxi(s, xi, eta));

			dpsiX_dxi[3 * node + 1] = br*Utilities::dN_seren_dxi(r, xi, eta) + bs*Utilities::dN_seren_dxi(s, xi, eta);

			dpsiX_dxi[3 * node + 2] = Utilities::dN_seren_dxi(node + 1, xi, eta) - cr*Utilities::dN_seren_dxi(r, xi, eta) - cs*Utilities::dN_seren_dxi(s, xi, eta);
			// Compute vector of dPsi_y/dxi
			dpsiY_dxi[3 * node] = 1.5 * (dr*Utilities::dN_seren_dxi(r, xi, eta) - ds*Utilities::dN_seren_dxi(s, xi, eta));

			dpsiY_dxi[3 * node + 1] = -1 * Utilities::dN_seren_dxi(node + 1, xi, eta) + er*Utilities::dN_seren_dxi(r, xi, eta) + es*Utilities::dN_seren_dxi(s, xi, eta);

			dpsiY_dxi[3 * node + 2] = -1 * br*Utilities::dN_seren_dxi(r, xi, eta) - bs*Utilities::dN_seren_dxi(s, xi, eta);

			// d( ) / deta
			// Compute vector of dPsi_x/deta
			dpsiX_deta[3 * node] = 1.5 * (ar*Utilities::dN_seren_deta(r, xi, eta) - as*Utilities::dN_seren_deta(s, xi, eta));
			dpsiX_deta[3 * node + 1] = br*Utilities::dN_seren_deta(r, xi, eta) + bs*Utilities::dN_seren_deta(s, xi, eta);
			dpsiX_deta[3 * node + 2] = Utilities::dN_seren_deta(node + 1, xi, eta) - cr*Utilities::dN_seren_deta(r, xi, eta) - cs*Utilities::dN_seren_deta(s, xi, eta);
			// Compute vector of dPsi_y/deta
			dpsiY_deta[3 * node] = 1.5 * (dr*Utilities::dN_seren_deta(r, xi, eta) - ds*Utilities::dN_seren_deta(s, xi, eta));
			dpsiY_deta[3 * node + 1] = -1 * Utilities::dN_seren_deta(node + 1, xi, eta) + er*Utilities::dN_seren_deta(r, xi, eta) + es*Utilities::dN_seren_deta(s, xi, eta);
			dpsiY_deta[3 * node + 2] = -1 * br*Utilities::dN_seren_deta(r, xi, eta) - bs*Utilities::dN_seren_deta(s, xi, eta);
		}

		double j11, j12, j21, j22;
		j11 = data.DKQ_invJac[data.gpIndex](0, 0);
		j12 = data.DKQ_invJac[data.gpIndex](0, 1);
		j21 = data.DKQ_invJac[data.gpIndex](1, 0);
		j22 = data.DKQ_invJac[data.gpIndex](1, 1);

		//try assembling like in original formulation
		Matrix B_bend_DKQ = Matrix(3, 12, 0.0);
		//B_bend_DKQ.clear();
		for (int col = 0; col < 12; col++)
		{
			B_bend_DKQ(0, col) += j11*dpsiX_dxi(col) + j12*dpsiX_deta(col);

			B_bend_DKQ(1, col) += j21*dpsiY_dxi(col) + j22*dpsiY_deta(col);

			B_bend_DKQ(2, col) += j11*dpsiY_dxi(col) + j12*dpsiY_deta(col) + j21*dpsiX_dxi(col) + j22*dpsiX_deta(col);
		}

		B_bend_total.clear();
		B_bend_total += B_bend_DKQ;



		//---------------------------------------------
		// assemble the membrane contribution
		// and the bending contribution
		// into the combined B matrix


		for (int nodeid = 0; nodeid < 4; nodeid++)
		{
			int i = nodeid * 3;
			int j = nodeid * 6;

			// membrane map: [0,1,5] <- [0,1,2]

			data.B(0, j) = B_mem(0, i);
			data.B(0, j + 1) = B_mem(0, i + 1);
			data.B(0, j + 5) = B_mem(0, i + 2);

			data.B(1, j) = B_mem(1, i);
			data.B(1, j + 1) = B_mem(1, i + 1);
			data.B(1, j + 5) = B_mem(1, i + 2);

			data.B(2, j) = B_mem(2, i);
			data.B(2, j + 1) = B_mem(2, i + 1);
			data.B(2, j + 5) = B_mem(2, i + 2);

			// bending map: [2,3,4] <- [0,1,2]

			data.B(3, j + 2) = B_bend_total(0, i);
			data.B(3, j + 3) = B_bend_total(0, i + 1);
			data.B(3, j + 4) = B_bend_total(0, i + 2);

			data.B(4, j + 2) = B_bend_total(1, i);
			data.B(4, j + 3) = B_bend_total(1, i + 1);
			data.B(4, j + 4) = B_bend_total(1, i + 2);

			data.B(5, j + 2) = B_bend_total(2, i);
			data.B(5, j + 3) = B_bend_total(2, i + 1);
			data.B(5, j + 4) = B_bend_total(2, i + 2);
		}

	}

	void ShellThinElement3D4N::CalculateBeta0(CalculationData& data)
	{
		data.beta0 = 1.0; // to be changed!
	}

	void ShellThinElement3D4N::CalculateSectionResponse(CalculationData& data)
	{

		// DO I NEED THIS??? CHECK!
		GeometryType & geom = GetGeometry();
		const Matrix & shapeFunctions = geom.ShapeFunctionsValues();
		Vector iN(shapeFunctions.size2());
		noalias(iN) = row(shapeFunctions, data.gpIndex);

		data.jacOp.Calculate(data.LCS0, geom.ShapeFunctionLocalGradient(data.gpIndex));
		data.SectionParameters.SetShapeFunctionsDerivatives(data.jacOp.XYDerivatives());

		ShellCrossSection::Pointer& section = mSections[data.gpIndex];
		data.SectionParameters.SetShapeFunctionsValues(iN);
		section->CalculateSectionResponse(data.SectionParameters, ConstitutiveLaw::StressMeasure_PK2);
	}

	void ShellThinElement3D4N::CalculateGaussPointContribution(CalculationData& data, MatrixType& LHS, VectorType& RHS)
	{
		//std::cout << "----------- Gauss Point " << data.gpIndex << std::endl;
		// calculate beta0
		CalculateBeta0(data);

		// calculate the total strain displ. matrix
		CalculateBMatrix(data);

		// compute generalized strains
		noalias(data.generalizedStrains) = prod(data.B, data.localDisplacements);

		// calculate section response
		CalculateSectionResponse(data);

		//printMatrix(data.D, "printing D matrix EARLY");

		// multiply the section tangent matrices and stress resultants by 'dA'
		data.D *= data.dA[data.gpIndex];
		//printMatrix(data.D, "printing D matrix LATE");

		// Add all contributions to the Stiffness Matrix
		data.BTD.clear();
		noalias(data.BTD) = prod(trans(data.B), data.D);
		noalias(LHS) += prod(data.BTD, data.B);

		// Add all contributions to the residual vector
		//noalias(RHS) -= prod(trans(data.B), data.generalizedStresses);
	}

	void ShellThinElement3D4N::CalculateAll(MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo,
		const bool LHSrequired,
		const bool RHSrequired)
	{
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
		CalculationData data(mpCoordinateTransformation, rCurrentProcessInfo);
		data.CalculateLHS = LHSrequired;
		data.CalculateRHS = RHSrequired;
		InitializeCalculationData(data);



		// Gauss Loop.
		for (size_t i = 0; i < OPT_NUM_GP; i++)
		{
			data.gpIndex = i;
			CalculateGaussPointContribution(data, rLeftHandSideMatrix, rRightHandSideVector);
		}

		//ApplyCorrectionToRHS(data, rRightHandSideVector);

		// Let the CoordinateTransformation finalize the calculation.
		// This will handle the transformation of the local matrices/vectors to
		// the global coordinate system.

		//printMatrix(rLeftHandSideMatrix, "Printing element K mat BEFORE finialize calc");

		mpCoordinateTransformation->FinalizeCalculations(data.LCS,
			data.globalDisplacements,
			data.localDisplacements,
			rLeftHandSideMatrix,
			rRightHandSideVector,
			RHSrequired,
			LHSrequired);

		// Add body forces contributions. This doesn't depend on the coordinate system

		//AddBodyForces(data, rRightHandSideVector);

		//printMatrix(rLeftHandSideMatrix, "Printing element K mat after finialize calc");
	}

	void ShellThinElement3D4N::printMatrix(Matrix& matrixIn, std::string stringIn)
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
			//std::cout << "|" << std::endl;
		}
		std::cout << std::endl;
	}

	void ShellThinElement3D4N::printVector(Vector& matrixIn, std::string stringIn)
	{
		std::cout << "\n" << stringIn << std::endl;
		for (unsigned i = 0; i<matrixIn.size(); ++i)
		{
			std::cout << "| ";

			std::cout << std::fixed << std::setprecision(6) << std::setw(12) << matrixIn(i) << " | ";

			std::cout << std::endl;
			//std::cout << "|" << std::endl;
		}
		std::cout << std::endl;
	}

	void ShellThinElement3D4N::printMatrix(const Matrix& matrixIn, std::string stringIn)
	{
		std::cout << "\n" << stringIn << std::endl;
		for (unsigned i = 0; i<matrixIn.size1(); ++i)
		{
			std::cout << "| ";
			for (unsigned j = 0; j<matrixIn.size2(); ++j)
			{
				std::cout << std::setprecision(2) << std::setw(10) << matrixIn(i, j) << " | ";
			}
			std::cout << std::endl;
			//std::cout << "|" << std::endl;
		}
		std::cout << std::endl;
	}

	void ShellThinElement3D4N::printDouble(std::string stringIn, double doubleIn)
	{
		std::cout << stringIn << doubleIn << std::endl;
	}

	// =====================================================================================
	//
	// CalculationData
	//
	// =====================================================================================

	ShellThinElement3D4N::CalculationData::CalculationData(const CoordinateTransformationBasePointerType& pCoordinateTransformation,
		const ProcessInfo& rCurrentProcessInfo)
		: LCS0(pCoordinateTransformation->CreateReferenceCoordinateSystem())
		, LCS(pCoordinateTransformation->CreateLocalCoordinateSystem())
		, CurrentProcessInfo(rCurrentProcessInfo)
	{
	}

	// =====================================================================================
	//
	// Class ShellThinElement3D4N - Serialization
	//
	// =====================================================================================

	void ShellThinElement3D4N::save(Serializer& rSerializer) const
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.save("CTr", mpCoordinateTransformation);
		rSerializer.save("Sec", mSections);
		rSerializer.save("IntM", (int)mThisIntegrationMethod);
	}

	void ShellThinElement3D4N::load(Serializer& rSerializer)
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
		rSerializer.load("CTr", mpCoordinateTransformation);
		rSerializer.load("Sec", mSections);
		int temp;
		rSerializer.load("IntM", temp);
		mThisIntegrationMethod = (IntegrationMethod)temp;
	}
}
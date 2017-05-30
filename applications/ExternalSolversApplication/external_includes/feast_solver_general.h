//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:		 BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $ExternalSolversApplication   $
//   Last modified by:    $Author:	peter.wilson@tum.de $
//   Date:                $Date:			  June 2017 $
//   Revision:            $Revision:                0.0 $
//
//

// System includes
#include <iostream>
#include <complex>
#include <vector>
#include <unordered_set>
#include <algorithm>

// External includes
#include <boost/smart_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/value_type/complex.hpp>
#include <amgcl/solver/skyline_lu.hpp>
extern "C" {
#include "feast.h"
}

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"
#include "includes/ublas_interface.h"
#include "spaces/ublas_space.h"

#if !defined(KRATOS_FEAST_SOLVER_GENERAL)
#define  KRATOS_FEAST_SOLVER_GENERAL

namespace Kratos {

///@name Kratos Classes
///@{

/*

template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class SkylineLUSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SkylineLUSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TSparseSpaceType::DataType DataType;
    
    typedef typename amgcl::backend::builtin<DataType>::matrix BuiltinMatrixType;
    
    typedef amgcl::solver::skyline_lu<DataType> SolverType;

    ~SkylineLUSolver()
    {
        Clear();
    }

    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Clear();

        pBuiltinMatrix = amgcl::adapter::zero_copy(
                rA.size1(),
                rA.index1_data().begin(),
                rA.index2_data().begin(),
                rA.value_data().begin());

        pSolver = boost::make_shared<SolverType>(*pBuiltinMatrix);
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        std::vector<DataType> x(rX.size());
        std::vector<DataType> b(rB.size());

        std::copy(std::begin(rB), std::end(rB), std::begin(b));

        (*pSolver)(b, x);

        std::copy(std::begin(x), std::end(x), std::begin(rX));

        return true;
    }

    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        Clear();
    }

    void Clear() override
    {
        pSolver.reset();
        pBuiltinMatrix.reset();
    }

private:

    boost::shared_ptr<BuiltinMatrixType> pBuiltinMatrix;

    boost::shared_ptr<SolverType> pSolver;
};

*/

/// Adapter to FEAST eigenvalue problem solver.
template<class TSparseSpaceType, class TDenseSpaceType,
        class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class FEASTSolverGeneral: public LinearSolver<TSparseSpaceType, TDenseSpaceType,
        TReordererType> {

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( FEASTSolverGeneral);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType SparseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef std::complex<double> ComplexType;

    typedef boost::numeric::ublas::compressed_matrix<ComplexType> ComplexSparseMatrixType;

    typedef boost::numeric::ublas::matrix<ComplexType> ComplexDenseMatrixType;

    typedef boost::numeric::ublas::vector<ComplexType> ComplexVectorType;

    typedef UblasSpace<ComplexType, ComplexSparseMatrixType, ComplexVectorType> ComplexSparseSpaceType;

    typedef UblasSpace<ComplexType, ComplexDenseMatrixType, ComplexVectorType> ComplexDenseSpaceType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;


    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor for built-in linear solver.
    /**
     * Parameters let the user control the settings of the FEAST library.
     */
	FEASTSolverGeneral(Parameters::Pointer pParam) : mpParam(pParam)
    {

		std::cout << "\n\n========= USING GENERAL FEAST SOLVER ==========\n\n" << std::endl;

        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {
                "solver_type": "skyline_lu"
            }
        })");

        mpParam->RecursivelyValidateAndAssignDefaults(default_params);

        if (mpParam->GetValue("linear_solver_settings")["solver_type"].GetString() != "skyline_lu")
            KRATOS_ERROR << "built-in solver type must be used with this constructor" << std::endl;

        mpLinearSolver = boost::make_shared<SkylineLUSolver<ComplexSparseSpaceType, ComplexDenseSpaceType>>();
    }

    /// Constructor for externally provided linear solver.
    /**
     * Parameters let the user control the settings of the FEAST library.
     * Warning: For iterative solvers, very small tolerances (~1e-15)
     *          may be needed for FEAST to work properly. Common iterative 
     *          solvers normally don't perform efficiently with FEAST 
     *          (M. Galgon et al., Parallel Computing (49) 2015 153-163).
     */
	FEASTSolverGeneral(Parameters::Pointer pParam, ComplexLinearSolverType::Pointer pLinearSolver)
        : mpParam(pParam), mpLinearSolver(pLinearSolver)
    {
		std::cout << "\n\n========= USING GENERAL FEAST SOLVER ==========\n\n" << std::endl;

        Parameters default_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "number_of_eigenvalues": 0,
            "search_dimension": 10,
            "linear_solver_settings": {}
        })");

        // don't validate linear_solver_settings here
        mpParam->ValidateAndAssignDefaults(default_params);
    }

    /// Deleted copy constructor.
	FEASTSolverGeneral(const FEASTSolverGeneral& Other) = delete;

    /// Destructor.
    virtual ~FEASTSolverGeneral() {}

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
	FEASTSolverGeneral& operator=(const FEASTSolverGeneral& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Solve the generalized eigenvalue problem.
    /**
	 * det(K - lamda M)z = 0
     * K are M are both real square matrices.
     */
    virtual void Solve(
			SparseMatrixType& K,
			SparseMatrixType& M,
			DenseVectorType& rEigenvalues,				// changed
			DenseMatrixType& rEigenvectors)		// changed
    {
		std::cout << "In solve function" << std::endl;

        const auto SystemSize = K.size1();

        Parameters& FEAST_Settings = *mpParam;
        const double EigenvalueRangeMin = FEAST_Settings["lambda_min"].GetDouble();
        const double EigenvalueRangeMax = FEAST_Settings["lambda_max"].GetDouble();
		//const double Emid = (EigenvalueRangeMin + EigenvalueRangeMax) / 2.0;
		const ComplexType Emid = ComplexType((EigenvalueRangeMin + EigenvalueRangeMax) / 2.0, 0.0);
		const double Eradius = EigenvalueRangeMax - Emid.real();

        int SearchDimension = FEAST_Settings["search_dimension"].GetInt();
        int NumEigenvalues = FEAST_Settings["number_of_eigenvalues"].GetInt();

        rEigenvalues.resize(SearchDimension,false);
        rEigenvectors.resize(SearchDimension,2*SystemSize,false);

		ComplexVectorType Eigenvalues = ComplexVectorType(SearchDimension);
		Eigenvalues.clear();
		ComplexDenseMatrixType Eigenvectors = ComplexDenseMatrixType(SearchDimension, SystemSize);
		Eigenvectors.clear();

        if (FEAST_Settings["perform_stochastic_estimate"].GetBool())
        {
            // this estimates the number of eigenvalues in the area of Emid and Eradius
            Calculate(M,K,Emid,Eradius,SearchDimension,
                    NumEigenvalues,Eigenvalues,Eigenvectors,true);

            std::cout << "Estimated number of eigenvalues = " << NumEigenvalues << std::endl;

            // recommended estimate of search dimension from FEAST documentation
            SearchDimension = NumEigenvalues + NumEigenvalues/2 + 1;
            FEAST_Settings["search_dimension"].SetInt(SearchDimension);
        }
        if (FEAST_Settings["solve_eigenvalue_problem"].GetBool())
        {
            // this attempts to solve the generalized eigenvalue problem
            Calculate(M,K,Emid,Eradius,SearchDimension,
                    NumEigenvalues,Eigenvalues,Eigenvectors,false);

            Eigenvalues.resize(NumEigenvalues,true);
            Eigenvectors.resize(NumEigenvalues,SystemSize,true);

        }
        FEAST_Settings["number_of_eigenvalues"].SetInt(NumEigenvalues);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "FEAST solver (general).";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Parameters::Pointer mpParam;

    ComplexLinearSolverType::Pointer mpLinearSolver;

    ///@}
    ///@name Private Operations
    ///@{

    /// Wrapper for FEAST library.
    void Calculate(
            SparseMatrixType& rMassMatrix,
            SparseMatrixType& rStiffnessMatrix,
			ComplexType Emid,
            double Eradius,
            int SearchDimension,
            int& rNumEigenvalues,
			ComplexVectorType& rEigenvalues,
			ComplexDenseMatrixType& rEigenvectors,
            const bool PerformStochasticEstimate)
    {
		KRATOS_TRY

		std::cout << "In calculate function" << std::endl;

        int FEAST_Params[64] = {};
        int NumIter, Info, SystemSize;
        double Epsout;
        DenseVectorType Residual(2*SearchDimension);										// good
        std::vector<std::complex<double> > IntegrationNodes, IntegrationWeights;
        SystemSize = static_cast<int>(rMassMatrix.size1());
		matrix<std::complex<double>, column_major> work(SystemSize,2*SearchDimension);	// good?
        matrix<std::complex<double>, column_major> zwork(SystemSize,SearchDimension);	// good?
        matrix<std::complex<double>, column_major> Aq(SearchDimension,SearchDimension);				// good?
        matrix<std::complex<double>, column_major> Bq(SearchDimension,SearchDimension);				// good?
        std::complex<double> Ze;
        ComplexSparseMatrixType Az;
        ComplexVectorType b(SystemSize);
        ComplexVectorType x(SystemSize);

        this->InitializeFEASTSystemMatrix(rMassMatrix, rStiffnessMatrix, Az);

        Parameters& FEAST_Settings = *mpParam;

        // initialize FEAST eigenvalue solver (see FEAST documentation for details)
        feastinit(FEAST_Params);
        if (FEAST_Settings["print_feast_output"].GetBool())
            FEAST_Params[0] = 1;
        FEAST_Params[2] = 8; // stopping convergence criteria 10^-FEAST_Params[2]
        FEAST_Params[28] = 1;// not sure if this is needed
        if (PerformStochasticEstimate)
        {
            FEAST_Params[1] = 4; // number of quadrature points (default: 8)
            FEAST_Params[13] = 2;
        }

        IntegrationNodes.resize(FEAST_Params[1]);
        IntegrationWeights.resize(FEAST_Params[1]);

        // get quadrature nodes and weights
		/*
        zfeast_contour(&EigenvalueRangeMin,
                &EigenvalueRangeMax,
                &FEAST_Params[1],
                &FEAST_Params[15],
                &FEAST_Params[17],
                (double *)IntegrationNodes.data(),
                (double *)IntegrationWeights.data());
				*/

		// get quadrature nodes and weights
		zfeast_gcontour((double *)&Emid,
				&Eradius,
				&FEAST_Params[7],
				&FEAST_Params[16],
				&FEAST_Params[17],
				&FEAST_Params[18],
				(double *)IntegrationNodes.data(),
				(double *)IntegrationWeights.data());

        int ijob = -1;
        // solve the eigenvalue problem
        while (ijob != 0)
        {
            // FEAST's reverse communication interface

			/*
            dfeast_srcix(&ijob,&SystemSize,(double *)&Ze,(double *)work.data().begin(),
                    (double *)zwork.data().begin(),(double *)Aq.data().begin(),
                    (double *)Bq.data().begin(),FEAST_Params,&Epsout,&NumIter,
                    &EigenvalueRangeMin,&EigenvalueRangeMax,&SearchDimension,
                    (double *)rEigenvalues.data().begin(),
                    (double *)rEigenvectors.data().begin(),
                    &rNumEigenvalues,(double *)Residual.data().begin(),&Info,
                    (double *)IntegrationNodes.data(),
                    (double *)IntegrationWeights.data());*/

			std::cout << "Before grcix call" << std::endl;

			// FEAST's reverse communication interface
			dfeast_grcix(&ijob,								// 
				&SystemSize,								// 
				(double *)&Ze,								// 
				(double *)work.data().begin(),				// workr - check?
				(double *)zwork.data().begin(),				// workc - check?
				(double *)Aq.data().begin(),				// Aq - check?
				(double *)Bq.data().begin(),				// Bq - check?
				FEAST_Params,								//
				&Epsout,									// 
				&NumIter,									// 
				(double *)&Emid,							// 
				&Eradius,									// 
				&SearchDimension,							// 
				(double *)rEigenvalues.data().begin(),		// lambda - check
				(double *)rEigenvectors.data().begin(),		// Q - check
				&rNumEigenvalues,							// 
				(double *)Residual.data().begin(),			// res - check
				&Info,										// 
				(double *)IntegrationNodes.data(),			// 
				(double *)IntegrationWeights.data());		// 
			
			std::cout << "After grcix call" << std::endl;

            switch (ijob)
            {
                case 10:
                {
                    // set up quadrature matrix (ZeM-K) and solver
                    this->CalculateFEASTSystemMatrix(Ze, rMassMatrix, rStiffnessMatrix, Az);
                    mpLinearSolver->Clear();
                    mpLinearSolver->Initialize(Az,x,b);
                    mpLinearSolver->InitializeSolutionStep(Az, x, b);
                } break;
                case 11:
                {
                    // solve the linear system for one quadrature point:
					// Az * Qz = zwork
                    for (int j=0; j < FEAST_Params[22]; j++)
                    {
                        for (int i=0; i < SystemSize; i++)
                            b[i] = zwork(i,j);
                        mpLinearSolver->Solve(Az,x,b);
                        for (int i=0; i < SystemSize; i++)
                            zwork(i,j) = x[i];
                    }
                } break;
				case 20:
					std::cout << "\n\n ----------- I'VE ENDED UP IN CASE 20! ------------- \n\n" << std::endl;
					// Nothing programmed here so far.
					break;
				case 21:
					// solve the linear system for one quadrature point
					// Az^H * Qz = zwork. SAME AS CASE 11 ATM ????????????????????
					for (int j = 0; j < FEAST_Params[22]; j++)
					{
						for (int i = 0; i < SystemSize; i++)
							b[i] = zwork(i, j);
						mpLinearSolver->Solve(Az, x, b);
						for (int i = 0; i < SystemSize; i++)
							zwork(i, j) = x[i];
					}
					break; 
                case 30:
                {
                    // multiply Kx
                    for (int i=0; i < FEAST_Params[24]; i++)
                    {
                        int k = FEAST_Params[23]-1+i;
                        noalias(column(work,k)) = prod(rStiffnessMatrix,row(rEigenvectors,k));
                    }
                } break;
				case 31:
					// SAME AS 30 BUT WITH INDICES CHANGED????????????????
					for (int i = 0; i < FEAST_Params[34]; i++)
					{
						int k = FEAST_Params[33] - 1 + i;
						noalias(column(work, k)) = prod(rStiffnessMatrix, row(rEigenvectors, k));
					}
					break;

                case 40:
                {
                    // multiply Mx
                    for (int i=0; i < FEAST_Params[24]; i++)
                    {
                        int k = FEAST_Params[23]-1+i;
                        noalias(column(work,k)) = prod(rMassMatrix,row(rEigenvectors,k));
                    }
				} break;
				case 41:
					// multiply Mx
					// SAME AS 40 BUT WITH INDICES CHANGED??????
					for (int i = 0; i < FEAST_Params[34]; i++)
					{
						int k = FEAST_Params[33] - 1 + i;
						noalias(column(work, k)) = prod(rMassMatrix, row(rEigenvectors, k));
					}
					break;
            } // switch
        } // while

        KRATOS_CATCH("")
    }

    /**
     * Initialize CSR matrix structure for FEAST system matrix: C = z * B - A.
     */
    void InitializeFEASTSystemMatrix(const SparseMatrixType& B,
                                     const SparseMatrixType& A,
                                     ComplexSparseMatrixType& C)
    {
        C.resize(B.size1(), B.size2(), false);

        std::vector<std::unordered_set<std::size_t> > indices(C.size1());

        // indices for row begin / end
        C.index1_data()[0] = 0;
        for (std::size_t i = 0; i < C.size1(); ++i)
        {
            std::size_t row_begin, row_end;
            indices[i].reserve(40); // initialize C's indices

            row_begin = B.index1_data()[i];
            row_end = B.index1_data()[i + 1];
            indices[i].insert(B.index2_data().begin() + row_begin,
                    B.index2_data().begin() + row_end); // insert B's column indices for row i

            row_begin = A.index1_data()[i];
            row_end = A.index1_data()[i + 1];
            indices[i].insert(A.index2_data().begin() + row_begin,
                    A.index2_data().begin() + row_end); // insert A's column indices for row i

            // C.index1_data()[i+1] = number of non-zeros in rows <= i
            C.index1_data()[i + 1] = C.index1_data()[i] + indices[i].size();
        }

        // C.index1_data()[C.size1()] = number of non-zeros
        C.reserve(C.index1_data()[C.size1()]);

        // column indices
        std::size_t k = 0;
        for (std::size_t i = 0; i < C.size1(); ++i)
        {
            for (std::size_t j : indices[i])
                C.index2_data()[k++] = j; // fill C's column indices

            indices[i].clear();

            std::sort(C.index2_data().begin() + C.index1_data()[i],
                    C.index2_data().begin() + C.index1_data()[i + 1]);
        }

        C.set_filled(C.size1() + 1, C.index1_data()[C.size1()]);
    }

    /**
     * Calculate FEAST system matrix: C = z * B - A. Similar to FEAST's zdaddcsr subroutine.
     */
    void CalculateFEASTSystemMatrix(std::complex<double> z,
                                    SparseMatrixType& B,
                                    SparseMatrixType& A,
                                    ComplexSparseMatrixType& C)
    {
        std::size_t jb, ja;
        const std::size_t dimension = B.size1();

        std::size_t ptr = 0;
        for (std::size_t i = 0; i < dimension; ++i)
        {
            std::size_t b_ptr = B.index1_data()[i];
            std::size_t a_ptr = A.index1_data()[i];
            while (b_ptr < B.index1_data()[i + 1] || a_ptr < A.index1_data()[i + 1])
            {
                jb = (b_ptr < B.index1_data()[i + 1]) ?
                        B.index2_data()[b_ptr] : dimension;
                ja = (a_ptr < A.index1_data()[i + 1]) ?
                        A.index2_data()[a_ptr] : dimension;

                if (jb < ja)
                {
                    C.value_data()[ptr] = z * B(i, jb).ref();
                    b_ptr++;
                }
                else if (jb > ja)
                {
                    C.value_data()[ptr] = -A(i, ja).ref();
                    a_ptr++;
                }
                else
                { // jb == ja
                    C.value_data()[ptr] = z * B(i, jb).ref() - A(i, ja).ref();
                    b_ptr++;
                    a_ptr++;
                }
                ptr++;
            }
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class FEASTSolverGeneral

///@}

///@name Input and output
///@{

/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        FEASTSolverGeneral<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const FEASTSolverGeneral<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}// namespace Kratos.

#endif // KRATOS_FEAST_SOLVER_GENERAL  defined

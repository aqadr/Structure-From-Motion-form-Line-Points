#include "MathHelper.h"


	template<typename PType>
	void sfm::lineNormalizedCoordinate(const Matrix<PType, 3, 3> &K, Matrix<PType, 4, 1> &EndPtsIm)
	{
		Matrix<PType, 4, 1> EndPtsTmp = EndPtsIm;
		Matrix<PType, 2, 1> PtStart = EndPtsTmp.block<2, 1>(0, 0);
		Matrix<PType, 2, 1> PtEnd = EndPtsTmp.block<2, 1>(2, 0);

		Matrix<PType, 3, 1> PtStart_H; PtStart_H.block<2, 1>(0, 0) = PtStart; PtStart_H[2] = 1.0;
		Matrix<PType, 3, 1> PtEnd_H; PtEnd_H.block<2, 1>(0, 0) = PtEnd; PtEnd_H[2] = 1.0;

		Matrix<PType, 3, 3> invK = K.inverse();
		Matrix<PType, 3, 1> PtStart_HN = invK*PtStart_H;
		Matrix<PType, 3, 1> PtEnd_HN = invK*PtEnd_H;

		EndPtsIm.block<2, 1>(0, 0) = PtStart_HN.block<2, 1>(0, 0);
		EndPtsIm.block<2, 1>(2, 0) = PtEnd_HN.block<2, 1>(0, 0);
	}

	//////////////////////////////////////////////////////
	bool sfm::normalize2DPoints(std::vector<Eigen::Vector2d> &vPts, Eigen::Matrix3d &T)
	{
		// shift origins to centroids
		double cu = 0, cv = 0;
		double m = static_cast<double>(vPts.size());
		for (std::vector<Eigen::Vector2d>::iterator it = vPts.begin(); it != vPts.end(); ++it)
		{
			cu += (*it)[0];
			cv += (*it)[1];
		}

		cu /= m;
		cv /= m;

		for (size_t i = 0; i<vPts.size(); ++i)
		{
			vPts[i][0] -= cu;
			vPts[i][1] -= cv;
		}

		//! scale features such that mean distance from origin is sqrt(2)
		double sp = 0;
		for (std::vector<Eigen::Vector2d>::iterator it = vPts.begin(); it != vPts.end(); it++)
		{
			Eigen::Vector2d pt = *it;
			sp += sqrt(pt[0] * pt[0] + pt[1] * pt[1]);
		}

		if (fabs(sp)<1e-10)
			return false;

		sp = sqrt(2.0)*m / sp;
		for (size_t i = 0; i<vPts.size(); ++i)
			vPts[i] *= sp;

		//! compute corresponding transformation matrices
		T << sp, 0, -sp*cu, 0, sp, -sp*cv, 0, 0, 1;

		//! return true on success
		return true;
	}


	////////////////////////////////////////////
	
	void sfm::constraintMinimization(MatrixXd A, MatrixXd G, VectorXd &X)
	{
		if (A.rows() >= A.cols())
		{
			Eigen::JacobiSVD<MatrixXd> svd(G, Eigen::ComputeThinU | Eigen::ComputeFullV);
			MatrixXd UMat = svd.matrixU();
			const int row = UMat.rows();

			const int r = svd.rank();
			MatrixXd U_hat(row, r);

			for (size_t i = 0; i<row; ++i)
				for (size_t j = 0; j<r; ++j)
					U_hat(i, j) = UMat(i, j);

			MatrixXd A_hat = A*U_hat;
			Eigen::JacobiSVD<MatrixXd> svd2(A_hat, Eigen::ComputeThinU | Eigen::ComputeThinV);
			MatrixXd V2 = svd2.matrixV();
			MatrixXd X2 = V2.col(V2.cols() - 1);
			X = U_hat*X2;

			return;
		}
		//! Extend A with rows of zeros to make it square. It's a hack, but is
		//! necessary until Eigen supports SVD with more columns than rows.
		MatrixXd A_extended(A.cols(), A.cols());
		A_extended.block(A.rows(), 0, A.cols() - A.rows(), A.cols()).setZero();
		A_extended.block(0, 0, A.rows(), A.cols()) = (A);

		MatrixXd G_extended(G.cols(), G.cols());
		G_extended.block(G.rows(), 0, G.cols() - G.rows(), G.cols()).setZero();
		G_extended.block(0, 0, G.rows(), G.cols()) = (G);
		constraintMinimization(A_extended, G_extended, X);
	}
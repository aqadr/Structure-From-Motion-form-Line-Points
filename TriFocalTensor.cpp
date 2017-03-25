#include "TriFocalTensor.h"
#include "MathHelper.h"
namespace sfm
{
	
	////////////////////////////////////////////////////////////////////
	TrifocalTensor::TrifocalTensor(const std::vector<Vector2d> &ptsFrame1,
		const std::vector<Vector2d> &ptsFrame2,
		const std::vector<Vector2d> &ptsFrame3)
	{
		ptsCorrespondences.resize(0);
		ptsCorrespondences.reserve(ptsFrame1.size());

		for (size_t i = 0; i < ptsFrame1.size(); ++i)
		{
			Vector2d pt1 = ptsFrame1[i];
			Vector2d pt2 = ptsFrame2[i];
			Vector2d pt3 = ptsFrame3[i];

			corrThree pt(pt1, pt2, pt3);
			ptsCorrespondences.push_back(pt);
		}
	}


	//////////////////////////////////////////////////////
	void TrifocalTensor::epipolsFromTrifocalTensor(Eigen::Matrix<double, 3, 9> T, Vector3d &e1, Vector3d &e2)
	{
		Matrix3d T1 = T.block<3, 3>(0, 0);
		Matrix3d T2 = T.block<3, 3>(0, 3);
		Matrix3d T3 = T.block<3, 3>(0, 6);
		Matrix3d T1_t = T1.transpose();
		Matrix3d T2_t = T2.transpose();
		Matrix3d T3_t = T3.transpose();

		Eigen::JacobiSVD<Matrix3d> svd1(T1, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::JacobiSVD<Matrix3d> svd2(T2, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::JacobiSVD<Matrix3d> svd3(T3, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::JacobiSVD<Matrix3d> svd1_t(T1_t, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::JacobiSVD<Matrix3d> svd2_t(T2_t, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::JacobiSVD<Matrix3d> svd3_t(T3_t, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Vector3d u1 = svd1_t.matrixV().col(2);
		Vector3d v1 = svd1.matrixV().col(2);

		Vector3d u2 = svd2_t.matrixV().col(2);  //need to double check a couple of things
		Vector3d v2 = svd2.matrixV().col(2);

		Vector3d u3 = svd3_t.matrixV().col(2);
		Vector3d v3 = svd3.matrixV().col(2);

		Matrix3d U_Mat = Matrix3d::Zero();
		U_Mat.block<3, 1>(0, 0) = u1;
		U_Mat.block<3, 1>(0, 1) = u2;
		U_Mat.block<3, 1>(0, 2) = u3;

		Matrix3d V_Mat = Matrix3d::Zero();
		V_Mat.block<3, 1>(0, 0) = v1;
		V_Mat.block<3, 1>(0, 1) = v2;
		V_Mat.block<3, 1>(0, 2) = v3;

		Eigen::JacobiSVD<Matrix3d> svdV(V_Mat.transpose(), Eigen::ComputeFullV);
		Eigen::JacobiSVD<Matrix3d> svdU(U_Mat.transpose(), Eigen::ComputeFullV);
		e2 = svdV.matrixV().col(2);
		e1 = svdU.matrixV().col(2);
		e2[0] /= e2[2];
		e2[1] /= e2[2];
		e2[2] = 1.0;

		e1[0] /= e1[2];
		e1[1] /= e1[2];
		e1[2] = 1.0;
	}


	/////////////////////////////////////
	void TrifocalTensor::EFromEpipoles(Eigen::Matrix<double, 27, 18> &E,
		const Eigen::Vector3d &e2,
		const Eigen::Vector3d &e3)
	{
		E = Eigen::Matrix<double, 27, 18>::Zero();
		Eigen::Matrix<double, 9, 3> e3Block = Eigen::Matrix<double, 9, 3>::Zero();
		e3Block.block<3, 1>(0, 0) = e3;
		e3Block.block<3, 1>(3, 1) = e3;
		e3Block.block<3, 1>(6, 2) = e3;

		Eigen::Matrix<double, 9, 3> e2Block = Eigen::Matrix<double, 9, 3>::Zero();

		e2Block(0, 0) = -e2[0];
		e2Block(1, 1) = -e2[0];
		e2Block(2, 2) = -e2[0];

		e2Block(3, 0) = -e2[1];
		e2Block(4, 1) = -e2[1];
		e2Block(5, 2) = -e2[1];

		e2Block(6, 0) = -e2[2];
		e2Block(7, 1) = -e2[2];
		e2Block(8, 2) = -e2[2];

		E.block<9, 3>(0, 0) = e3Block;
		E.block<9, 3>(9, 3) = e3Block;
		E.block<9, 3>(18, 6) = e3Block;

		E.block<9, 3>(0, 9) = e2Block;
		E.block<9, 3>(9, 12) = e2Block;
		E.block<9, 3>(18, 15) = e2Block;
	}


	/////////////////////////////////////////////
	void TrifocalTensor::computeEFromEpipoles(Eigen::MatrixXd &E, const Eigen::Vector3d &e2, const Eigen::Vector3d &e3)
	{
		E = Eigen::MatrixXd::Zero(27, 18);
		Eigen::Matrix<double, 9, 3> e3Block = Eigen::Matrix<double, 9, 3>::Zero();
		e3Block.block<3, 1>(0, 0) = e3;
		e3Block.block<3, 1>(3, 1) = e3;
		e3Block.block<3, 1>(6, 2) = e3;

		Eigen::Matrix<double, 9, 3> e2Block = Eigen::Matrix<double, 9, 3>::Zero();

		e2Block(0, 0) = -e2[0];
		e2Block(1, 1) = -e2[0];
		e2Block(2, 2) = -e2[0];

		e2Block(3, 0) = -e2[1];
		e2Block(4, 1) = -e2[1];
		e2Block(5, 2) = -e2[1];

		e2Block(6, 0) = -e2[2];
		e2Block(7, 1) = -e2[2];
		e2Block(8, 2) = -e2[2];

		E.block<9, 3>(0, 0) = e3Block;
		E.block<9, 3>(9, 3) = e3Block;
		E.block<9, 3>(18, 6) = e3Block;

		E.block<9, 3>(0, 9) = e2Block;
		E.block<9, 3>(9, 12) = e2Block;
		E.block<9, 3>(18, 15) = e2Block;
	}


	/////////////////////////////////////////////////
	void TrifocalTensor::Tfrom3Ps(Eigen::Matrix<double, 3, 9>&T, const Eigen::Matrix<double, 3, 4> &P1, const Eigen::Matrix<double, 3, 4> &P2,
		const Eigen::Matrix<double, 3, 4> &P3)
	{
		Eigen::Matrix<double, 3, 1> b4 = P3.block<3, 1>(0, 3);
		Eigen::Matrix<double, 3, 1> a4 = P2.block<3, 1>(0, 3);
		Eigen::Matrix3d T1 = P2.block<3, 1>(0, 0)*b4.transpose() - (a4*P3.block<3, 1>(0, 0).transpose());
		Eigen::Matrix3d T2 = P2.block<3, 1>(0, 1)*b4.transpose() - (a4*P3.block<3, 1>(0, 1).transpose());
		Eigen::Matrix3d T3 = P2.block<3, 1>(0, 2)*b4.transpose() - (a4*P3.block<3, 1>(0, 2).transpose());
		T.block<3, 3>(0, 0) = T1;
		T.block<3, 3>(0, 3) = T2;
		T.block<3, 3>(0, 6) = T3;
	}


	


	////////////////////////////////////////////////////////////
	bool TrifocalTensor::computeTrifocalTensorPoints(const std::vector<Eigen::Vector2d> &vImPts1_,
		const std::vector<Eigen::Vector2d> &vImPts2_,
		const std::vector<Eigen::Vector2d> &vImPts3_)
	{
		if (vImPts1_.size() != vImPts2_.size() || vImPts1_.size() != vImPts3_.size())
		return false;

		if (vImPts1_.size() <= 7)
		return false;

		std::vector<Eigen::Vector2d> vImPts1 = vImPts1_;
		std::vector<Eigen::Vector2d> vImPts2 = vImPts2_;
		std::vector<Eigen::Vector2d> vImPts3 = vImPts3_;

		Matrix3d H1, H2, H3;
		sfm::normalize2DPoints(vImPts1, H1);
		sfm::normalize2DPoints(vImPts2, H2);
		sfm::normalize2DPoints(vImPts3, H3);

		//now create the n*27 matrix
		//each set of points correspondences provide 4 equations
		int n = 4 * vImPts1_.size();
		Eigen::MatrixXd A = MatrixXd::Zero(n, 27);
		//populate the matrix with the points correspondences
		for (size_t ii = 0; ii<vImPts1.size(); ii++)
		{
			Vector3d pt1 = Vector3d(vImPts1[ii][0], vImPts1[ii][1], 1.0);
			Vector3d pt2 = Vector3d(vImPts2[ii][0], vImPts2[ii][1], 1.0);
			Vector3d pt3 = Vector3d(vImPts3[ii][0], vImPts3[ii][1], 1.0);

			size_t i = 4 * ii;
			//equation number 1
			A(i, 0) = pt1[0];
			A(i, 2) = -pt1[0] * pt3[0];
			A(i, 6) = -pt1[0] * pt2[0];
			A(i, 8) = (pt1[0] * pt2[0] * pt3[0]);

			A(i, 9) = pt1[1];
			A(i, 11) = -pt1[1] * pt3[0];
			A(i, 15) = -pt1[1] * pt2[0];
			A(i, 17) = (pt1[1] * pt2[0] * pt3[0]);

			A(i, 18) = pt1[2];
			A(i, 20) = -pt1[2] * pt3[0];;
			A(i, 24) = -pt1[2] * pt2[0];
			A(i, 26) = (pt1[2] * pt2[0] * pt3[0]);

			//equation number 2
			A(i + 1, 3) = pt1[0];
			A(i + 1, 5) = -pt1[0] * pt3[0];
			A(i + 1, 6) = -pt1[0] * pt2[1];
			A(i + 1, 8) = pt1[0] * pt2[1] * pt3[0];

			A(i + 1, 12) = pt1[1];
			A(i + 1, 14) = -pt1[1] * pt3[0];
			A(i + 1, 15) = -pt1[1] * pt2[1];
			A(i + 1, 17) = pt1[1] * pt2[1] * pt3[0];

			A(i + 1, 21) = pt1[2];
			A(i + 1, 23) = -pt1[2] * pt3[0];
			A(i + 1, 24) = -pt1[2] * pt2[1];
			A(i + 1, 26) = pt1[2] * pt2[1] * pt3[0];

			//equation number 3
			A(i + 2, 1) = pt1[0];
			A(i + 2, 2) = -pt1[0] * pt3[1];
			A(i + 2, 7) = -pt1[0] * pt2[0];
			A(i + 2, 8) = pt1[0] * pt2[0] * pt3[1];

			A(i + 2, 10) = pt1[1];
			A(i + 2, 11) = -pt1[1] * pt3[1];
			A(i + 2, 16) = -pt1[1] * pt2[0];
			A(i + 2, 17) = pt1[1] * pt2[0] * pt3[1];

			A(i + 2, 19) = pt1[2];
			A(i + 2, 20) = -pt1[2] * pt3[1];
			A(i + 2, 25) = -pt1[2] * pt2[0];
			A(i + 2, 26) = pt1[2] * pt2[0] * pt3[1];

			//equation number 4
			A(i + 3, 4) = pt1[0];
			A(i + 3, 5) = -pt1[0] * pt3[1];
			A(i + 3, 7) = -pt1[0] * pt2[1];
			A(i + 3, 8) = pt1[0] * pt2[1] * pt3[1];

			A(i + 3, 13) = pt1[1];
			A(i + 3, 14) = -pt1[1] * pt3[1];
			A(i + 3, 16) = -pt1[1] * pt2[1];
			A(i + 3, 17) = pt1[1] * pt2[1] * pt3[1];

			A(i + 3, 22) = pt1[2];
			A(i + 3, 23) = -pt1[2] * pt3[1];
			A(i + 3, 25) = -pt1[2] * pt2[1];
			A(i + 3, 26) = pt1[2] * pt2[1] * pt3[1];

		}

		/*QString filename = "C:/SFM_Simulation/output/MatA.txt";
		QFile file(filename);
		if (file.open(QIODevice::ReadWrite))
		{
			QTextStream stream(&file);
			for (size_t ii = 0; ii < A.rows(); ++ii)
			{
				for (size_t jj = 0; jj < A.cols(); ++jj)
				{
					stream << QString::number(A(ii, jj));
					stream << "\t";
				}
					stream << "\n";
			}
		}
		file.close();
		*/
		Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

		//Extract trifocal tensor from the column of V corresponding to
		//smallest singular value.
		Eigen::Matrix<double, 27, 1> t = svd.matrixV().col(26);

		//crate the trifocal tensor matrix from the 27 vector
		Eigen::Matrix<double, 3, 9>T0;
		Eigen::Matrix<double, 3, 3> T1, T2, T3;
		T1 << t(0, 0), t(1, 0), t(2, 0), t(3, 0), t(4, 0), t(5, 0), t(6, 0), t(7, 0), t(8, 0);
		T2 << t(9, 0), t(10, 0), t(11, 0), t(12, 0), t(13, 0), t(14, 0), t(15, 0), t(16, 0), t(17, 0);
		T3 << t(18, 0), t(19, 0), t(20, 0), t(21, 0), t(22, 0), t(23, 0), t(24, 0), t(25, 0), t(26, 0);


		T0.block<3, 3>(0, 0) = T1;
		T0.block<3, 3>(0, 3) = T2;
		T0.block<3, 3>(0, 6) = T3;


		//Ensure that that trifocal tensor is geomatrically valid by retrieving
		//its epipoles and performing algebraic minimization

		Eigen::Vector3d e1, e2;
		epipolsFromTrifocalTensor(T0, e1, e2);
		//std::cout << "check epipoles from trifocaltensor :" << e1 << "\n" << e2 << std::endl;
		Eigen::MatrixXd E;
		computeEFromEpipoles(E, e1, e2);
		Eigen::VectorXd tFinal;
		sfm::constraintMinimization(A, E, tFinal);


		Eigen::Matrix<double, 3, 9>TFinal;
		T1 << tFinal(0, 0), tFinal(1, 0), tFinal(2, 0), tFinal(3, 0), tFinal(4, 0), tFinal(5, 0), tFinal(6, 0), tFinal(7, 0), tFinal(8, 0);
		T2 << tFinal(9, 0), tFinal(10, 0), tFinal(11, 0), tFinal(12, 0), tFinal(13, 0), tFinal(14, 0), tFinal(15, 0), tFinal(16, 0), tFinal(17, 0);
		T3 << tFinal(18, 0), tFinal(19, 0), tFinal(20, 0), tFinal(21, 0), tFinal(22, 0), tFinal(23, 0), tFinal(24, 0), tFinal(25, 0), tFinal(26, 0);
		TFinal.block<3, 3>(0, 0) = T1;
		TFinal.block<3, 3>(0, 3) = T2;
		TFinal.block<3, 3>(0, 6) = T3;


		Matrix3d inv_H2 = H2.inverse();
		Matrix3d inv_H3 = H3.inverse();
		Matrix3d Y1 = inv_H2*T1*inv_H3.transpose();
		Matrix3d Y2 = inv_H2*T2*inv_H3.transpose();
		Matrix3d Y3 = inv_H2*T3*inv_H3.transpose();

		Matrix3d T1_F = H1(0, 0)*Y1 + H1(1, 0)*Y2 + H1(2, 0)*Y3;
		Matrix3d T2_F = H1(0, 1)*Y1 + H1(1, 1)*Y2 + H1(2, 1)*Y3;
		Matrix3d T3_F = H1(0, 2)*Y1 + H1(1, 2)*Y2 + H1(2, 2)*Y3;

		std::cout << "check T: " << T1_F <<" \n\n"<<T2_F<<"  \n\n"<<T3_F<< std::endl;
	}

	
	////////////////////////////////////////////////////////////
	/*bool TrifocalTensor::computeTrifocalTensorFromLinesLinear(const std::vector<Eigen::Vector4d> &vLines1,
								const std::vector<Eigen::Vector4d> &vLines2,
								const std::vector<Eigen::Vector4d> &vLines3)
	{
		if (vLines1.size() != vLines2.size() || vLines2.size() != vLines2.size())
			std::cout << "we have big issues: " << std::endl;

		if (vLines1.size() < 13)
			std::cout << "can not generate the tri focal tensor :" << std::endl;

		int N = vLines1.size();
		MatrixXd A = Eigen::MatrixXd::Zero(N * 2, 27);

		Matrix3d H1, H2, H3;
		std::vector<Eigen::Vector4d> vLines1_ = vLines1;
		std::vector<Eigen::Vector4d> vLines2_ = vLines2;
		std::vector<Eigen::Vector4d> vLines3_ = vLines3;

		PoseFromLines::normalize2DLinesv2(vLines1_, H1);
		PoseFromLines::normalize2DLinesv2(vLines2_, H2);
		PoseFromLines::normalize2DLinesv2(vLines3_, H3);

		//std::cout << "check H1 H2 and H3: " << H1 << "\n\n" << H2 << "\n\n" << H3 << std::endl;

		//generate the lines first from the vectors
		std::vector<Eigen::Vector3d> vLines2D1;
		std::vector<Eigen::Vector3d> vLines2D2;
		std::vector<Eigen::Vector3d> vLines2D3;
		for (size_t i = 0; i < N; ++i)
		{
			Eigen::Vector4d EndPts1 = vLines1_[i];
			Vector3d l1;
			Line2D::computeLineCoeffv2(EndPts1, l1);
			
			vLines2D1.push_back(l1);

			Eigen::Vector4d EndPts2 = vLines2_[i];
			Vector3d l2;
			Line2D::computeLineCoeffv2(EndPts2, l2);
			vLines2D2.push_back(l2);

			Eigen::Vector4d EndPts3 = vLines2_[i];
			Vector3d l3;
			Line2D::computeLineCoeffv2(EndPts2, l3);
			vLines2D3.push_back(l3);
		}

		for (std::vector<Eigen::Vector3d>::iterator it = vLines2D1.begin(); it != vLines2D1.end(); ++it)
		{
			Vector3d &l = *it;
			std::cout << "line before: " << l << std::endl;
			*it = H1.inverse()*l;
			std::cout << "line after: " << l << std::endl;

		}

		for (std::vector<Eigen::Vector3d>::iterator it = vLines2D2.begin(); it != vLines2D2.end(); ++it)
		{
			Vector3d &l = *it;
			*it = H2.inverse()*l;
		}

		for (std::vector<Eigen::Vector3d>::iterator it = vLines2D3.begin(); it != vLines2D3.end(); ++it)
		{
			Vector3d &l = *it;
			*it = H3.inverse()*l;
		}

		
		std::cout << "check H1 H2 and H3: " << H1 << "\n\n" << H2 << "\n\n" << H3 << std::endl;

		for (size_t i = 0; i < N; ++i)
		{
			Vector3d l1 = vLines2D1[i];
			Vector3d l2 = vLines2D2[i];
			Vector3d l3 = vLines2D3[i];

			//explicitly declare variables to make it easier
			double l11 = l2[0], l12 = l2[1], l13 = l2[2], l_1 = l1[0], l_2 = l1[1], l_3 = l1[2],l21=l3[0];
			double l22 = l3[1], l23 = l3[2];

			size_t ii = 2 * i;

			//equation number 1
			//(-l11*l21*l_3)*T1_11
			A(ii, 0) = (-l11*l21*l_3);

			//(-l11*l22*l_3)*T1_12 
			A(ii, 1) = (-l11*l22*l_3);

			//(-l11*l23*l_3)*T1_13 
			A(ii, 2) = (-l11*l23*l_3);

			//(-l12*l21*l_3)*T1_21 
			A(ii, 3) = (-l12*l21*l_3);

			//(-l12*l22*l_3)*T1_22 
			A(ii, 4) = (-l12*l22*l_3);

			//(-l12*l23*l_3)*T1_23
			A(ii, 5) = (-l12*l23*l_3);

			A(ii, 6) = (-l13*l21*l_3);

			//(-l13*l22*l_3)*T1_32 
			A(ii, 7) = (-l13*l22*l_3);

			//(-l13*l23*l_3)*T1_33 
			A(ii, 8) = (-l13*l23*l_3);

			//(l11*l21*l_1)*T3_11 
			A(ii, 18) = (l11*l21*l_1);

			//(l11*l22*l_1)*T3_12 
			A(ii, 19) = (l11*l22*l_1);

			//(l11*l23*l_1)*T3_13
			A(ii, 20) = (l11*l23*l_1);

			//(l12*l21*l_1)*T3_21
			A(ii, 21) = (l12*l21*l_1);

			//(l12*l22*l_1)*T3_22 
			A(ii, 22) = (l12*l22*l_1);

			//(l12*l23*l_1)*T3_23
			A(ii, 23) = (l12*l23*l_1);

			//(l13*l21*l_1)*T3_31
			A(ii, 24) = (l13*l21*l_1);

			//(l13*l22*l_1)*T3_32
			A(ii, 25) = (l13*l22*l_1);

			//(l13*l23*l_1)*T3_33
			A(ii, 26) = (l13*l23*l_1);

			//equation number 2
			A(ii + 1, 0) = (l11*l21*l_2);   //T1_11
			A(ii + 1, 1) = (l11*l22*l_2);   //T1_12
			A(ii + 1, 2) = (l11*l23*l_2);   // T1_13
			A(ii + 1, 3) = (l12*l21*l_2);      // T1_21
			//(l12*l22*l_2)*T1_22 
			A(ii + 1, 4) = (l12*l22*l_2);

			//(l12*l23*l_2)*T1_23 
			A(ii + 1, 5) = (l12*l23*l_2);

			//(l13*l21*l_2)*T1_31 
			A(ii + 1, 6) = (l13*l21*l_2);

			//(l13*l22*l_2)*T1_32 
			A(ii + 1, 7) = (l13*l22*l_2);

			//(l13*l23*l_2)*T1_33
			A(ii + 1, 8) = (l13*l23*l_2);   // T1_33

			//(-l11*l21*l_1)*T2_11
			A(ii + 1, 9) = (-l11*l21*l_1);

			//(-l11*l22*l_1)*T2_12 
			A(ii + 1, 10) = (-l11*l22*l_1);

			//(-l11*l23*l_1)*T2_13 
			A(ii + 1, 11) = (-l11*l23*l_1);

			//(-l12*l21*l_1)*T2_21
			A(ii + 1, 12) = (-l12*l21*l_1);

			//(-l12*l22*l_1)*T2_22 
			A(ii + 1, 13) = (-l12*l22*l_1);

			//(-l12*l23*l_1)*T2_23
			A(ii + 1, 14) = (-l12*l23*l_1);

			//(-l13*l21*l_1) *T2_31
			A(ii + 1, 15) = (-l13*l21*l_1);// *T2_31

			//(-l13*l22*l_1)*T2_32
			A(ii + 1, 16) = (-l13*l22*l_1);

			//(-l13*l23*l_1)*T2_33
			A(ii + 1, 17) = (-l13*l23*l_1);

		}

		Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

		//Extract trifocal tensor from the column of V corresponding to
		//smallest singular value.
		Eigen::Matrix<double, 27, 1> t = svd.matrixV().col(26);
		
		// T0 = permute( reshape(t0,3,3,3) , [2 1 3] );
		//crate the trifocal tensor matrix from the 27 vector
		Eigen::Matrix<double, 3, 9>T0;
		Eigen::Matrix<double, 3, 3> T1, T2, T3;
		T1 << t(0, 0), t(1, 0), t(2, 0), t(3, 0), t(4, 0), t(5, 0), t(6, 0), t(7, 0), t(8, 0);
		T2 << t(9, 0), t(10, 0), t(11, 0), t(12, 0), t(13, 0), t(14, 0), t(15, 0), t(16, 0), t(17, 0);
		T3 << t(18, 0), t(19, 0), t(20, 0), t(21, 0), t(22, 0), t(23, 0), t(24, 0), t(25, 0), t(26, 0);


		T0.block<3, 3>(0, 0) = T1;
		T0.block<3, 3>(0, 3) = T2;
		T0.block<3, 3>(0, 6) = T3;


		//Ensure that that trifocal tensor is geomatrically valid by retrieving
		//its epipoles and performing algebraic minimization

		Eigen::Vector3d e1, e2;
		epipolsFromTrifocalTensor(T0, e1, e2);
		//std::cout << "check epipoles from trifocaltensor :" << e1 << "\n" << e2 << std::endl;
		Eigen::MatrixXd E;
		computeEFromEpipoles(E, e1, e2);
		//std::cout<<"check E size: "<<E.rows()<<"  "<<E.cols()<<std::endl;
		//std::cout << "check E: \n" << E << std::endl;
		Eigen::VectorXd tFinal;
		sfm::constraintMinimization(A, E, tFinal);

		std::cout << "check tf: " << tFinal << std::endl;
		Eigen::Matrix<double, 3, 9>TFinal;
		T1 << tFinal(0, 0), tFinal(1, 0), tFinal(2, 0), tFinal(3, 0), tFinal(4, 0), tFinal(5, 0), tFinal(6, 0), tFinal(7, 0), tFinal(8, 0);
		T2 << tFinal(9, 0), tFinal(10, 0), tFinal(11, 0), tFinal(12, 0), tFinal(13, 0), tFinal(14, 0), tFinal(15, 0), tFinal(16, 0), tFinal(17, 0);
		T3 << tFinal(18, 0), tFinal(19, 0), tFinal(20, 0), tFinal(21, 0), tFinal(22, 0), tFinal(23, 0), tFinal(24, 0), tFinal(25, 0), tFinal(26, 0);
		TFinal.block<3, 3>(0, 0) = T1;
		TFinal.block<3, 3>(0, 3) = T2;
		TFinal.block<3, 3>(0, 6) = T3;


		Matrix3d inv_H2 = H2.inverse();
		Matrix3d inv_H3 = H3.inverse();
		Matrix3d Y1 = inv_H2*T1*inv_H3.transpose();
		Matrix3d Y2 = inv_H2*T2*inv_H3.transpose();
		Matrix3d Y3 = inv_H2*T3*inv_H3.transpose();

		Matrix3d T1_F = H1(0, 0)*Y1 + H1(1, 0)*Y2 + H1(2, 0)*Y3;
		Matrix3d T2_F = H1(0, 1)*Y1 + H1(1, 1)*Y2 + H1(2, 1)*Y3;
		Matrix3d T3_F = H1(0, 2)*Y1 + H1(1, 2)*Y2 + H1(2, 2)*Y3;

		std::cout << "check T: " << T1_F << " \n\n" << T2_F << "  \n\n" << T3_F << std::endl;

		return true;

	}
	*/


}





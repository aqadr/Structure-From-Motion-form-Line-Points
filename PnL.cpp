/****************************************************************************
*
*   Copyright (c) 2016 Ashraf Qadir All rights reserved.
*   Author: Ashraf Qadir, <ashraf.qadr@gmail.com>
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
* OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
****************************************************************************/

/**
* @file PnL.cpp
* @brief Pose from 3D-2D line correspondences from the paper
* @brief 
* @author Ashraf Qadir, <ashraf.qadr@gmail.com>
*/

#include "PnL.h"
#include <set>

PnL::PnL()
{
	coeffVector = Eigen::Matrix<double, 9, 1>::Zero();
}


//////////////////////////////////////////////////////////////////
bool PnL::computeRotationP3L(const std::vector<Matrix<double, 6, 1> > &v3Lines,
	const std::vector<Matrix<double, 4, 1> > &v2Lines,
	const Matrix<double, 3, 3> &K, Matrix3d &Rcm)
{
	if (v3Lines.size() != 3 || v2Lines.size() != 3)
		return false;

	Matrix3d Rmw;
	bool success=rotationWorldtoModel(v3Lines, Rmw);

	//step 1: get the direction unit vector m1, m2, m3 for the lines
	Vector3d v0 = v3Lines[0].block<3, 1>(3, 0) - v3Lines[0].block<3, 1>(0, 0);
	Vector3d v1 = v3Lines[1].block<3, 1>(3, 0) - v3Lines[1].block<3, 1>(0, 0);
	Vector3d v2 = v3Lines[2].block<3, 1>(3, 0) - v3Lines[2].block<3, 1>(0, 0);
	v0.normalize();
	v1.normalize();
	v2.normalize();

	//step 2: get 3 points on 3 lines
	Vector3d P0 = v3Lines[0].block<3, 1>(0, 0);
	Vector3d P1 = v3Lines[1].block<3, 1>(0, 0);
	Vector3d P2 = v3Lines[2].block<3, 1>(0, 0);

	//step 3: generate the normals with the 2D lines
	Matrix<double, 4, 1> lEndPtsL0 = v2Lines[0];
	Matrix<double, 4, 1> lEndPtsL1 = v2Lines[1];
	Matrix<double, 4, 1> lEndPtsL2 = v2Lines[2];
	Line2D::lineNormalizedCoordinate(K, lEndPtsL0);
	Line2D::lineNormalizedCoordinate(K, lEndPtsL1);
	Line2D::lineNormalizedCoordinate(K, lEndPtsL2);

	Vector3d p1(lEndPtsL0(0, 0), lEndPtsL0(1, 0), 1);
	Vector3d p2(lEndPtsL0(2, 0), lEndPtsL0(3, 0), 1);


	Vector3d n0 = p1.cross(p2);
	n0.normalize();

	p1 = Vector3d(lEndPtsL1(0, 0), lEndPtsL1(1, 0), 1);
	p2 = Vector3d(lEndPtsL1(2, 0), lEndPtsL1(3, 0), 1);

	Vector3d n1 = p1.cross(p2);
	n1.normalize();

	p1 = Vector3d(lEndPtsL2(0, 0), lEndPtsL2(1, 0), 1);
	p2 = Vector3d(lEndPtsL2(2, 0), lEndPtsL2(3, 0), 1);

	Vector3d n2 = p1.cross(p2);
	n2.normalize();


	//step 4: then get the rotation matrix from world to model
	//v0 need to be [0 0 1]'
	//we also need to transform v1 and v2 
	Vector3d v0_ = Rmw*v0;
	Vector3d v1_ = Rmw*v1;
	Vector3d v2_ = Rmw*v2;

	Matrix3d vMat;
	vMat.block<3, 1>(0, 0) = v0_;
	vMat.block<3, 1>(0, 1) = v1_;
	vMat.block<3, 1>(0, 2) = v2_;
	Eigen::FullPivLU<Eigen::Matrix<double, 3, 3> > lu_decomp(vMat);
	int rank = lu_decomp.rank();


	//step 5: generate the coefficients
	Matrix3d R_;
	success = generateRotMat(n0, v0_, R_);
	if (!success)
	{
		//TODO use a proper error output
		std::cout << "could not generate the initial rotation matrix: " << std::endl;
		return false;
	}
	
	//step 6: solve the polynomial and get the real roots
	Matrix<double, 1, 3> n1_ = n1.transpose() *R_;
	Matrix<double, 1, 3> n2_ = n2.transpose() *R_;

	//std::cout << "check n1_ and n2_: " << n1_ << " \n" << n2_ << std::endl;


	constructPolynomCoeffs(R_, n1, n2, v1_, v2_);
	success = computeRoots();
	if (realRoots.size() <= 0)
	{
		std::cout << "could not compute real roots for the system" << std::endl;
		return false;
	}

	//step 7: get the rotation from the cos alpha 
	//in order to get rotation we need to 
	Matrix3d Rcw;
	for (size_t i = 0; i < realRoots.size(); ++i)
	{
		double calpha = realRoots.at(i);
		double salpha = sqrt(1 - calpha*calpha);
		//std::cout << "check salpha: " << salpha << std::endl;
		double alpha = atan(salpha / calpha);

		double sigma1 = v1_[1] * n1_(0, 1)*calpha + v1_[1] * n1_(0, 2)*salpha + (v1_[0] * n1_(0, 0));
		double sigma2 = v1_[0] * n1_(0, 1)*calpha + v1_[0] * n1_(0, 2)*salpha - (v1_[1] * n1_(0, 0));
		double sigma3 = v1_[2] * n1_(0, 2)*calpha - v1_[2] * n1_(0, 1)*salpha;

		double sigma4 = v2_[1] * n2_(0, 1)*calpha + v2_[1] * n2_(0, 2)*salpha + v2_[0] * n2_(0, 0);
		double sigma5 = v2_[0] * n2_(0, 1)*calpha + v2_[0] * n2_(0, 2)*salpha - v2_[1] * n2_(0, 0);
		double sigma6 = v2_[2] * n2_(0, 2)*calpha - v2_[2] * n2(0, 1)*salpha;

		double beta;
		if (sigma3 == 0 && (sigma1 != 0 && sigma2 != 0))
			beta = atan(-sigma1 / sigma2);
		else if (sigma6 == 0 && (sigma4 != 0 && sigma5 != 0))
			beta = atan(-sigma4 / sigma5);
		else if ((sigma1 != 0 && sigma5 != 0) || (sigma2 != 0 && sigma4 != 0))
		{
			double cBeta = (sigma2*sigma6 - sigma3*sigma5) / (sigma1*sigma5 - sigma2*sigma4);
			beta = acos(cBeta);
		}
		
		//std::cout << "check sigmas :" << sigma1 << " " << sigma2 << "  " << sigma3 << " " << sigma4 << " " << sigma5 << " " << sigma6 << std::endl;
		//std::cout << "check alpha and beta: " << alpha << "\n" << beta << std::endl;

		//now create the rotation matrix to transform from world to camera coordinate to 
		double cb = cos(beta);
		double sb = sin(beta);
		Matrix3d Rz = Matrix3d::Identity();
		Rz(1, 1) = Rz(0, 0) = cb;
		Rz(0, 1) = -sb;
		Rz(1, 0) = sb;

		Matrix3d Rx = Matrix3d::Identity();
		Rx(1, 1) = Rx(2, 2) = calpha;
		Rx(1, 2) = salpha;
		Rx(2, 1) = -salpha;

		Matrix3d R_Final = R_*Rx*Rz*Rmw;
		std::cout << "check R_Final: " << R_Final << std::endl;

	}

	return success;
	//step 9: solve for translation
	//Vector3d t;
	//success=computeTranslation(Rcw, n0, n1, n2, P0, P1, P2,t);

}

////////////////////////////////////////////////////////////////////////////
bool PnL::rotationWorldtoModel(const std::vector<Matrix<double, 6, 1> > &v3Lines, Matrix3d &Rmw)
{
	if (v3Lines.size() != 3)
	{
		std::cout << "3 lines are required " << std::endl;
		return false;
	}
	Matrix<double, 6, 1> l3d0 = v3Lines[0];
	Matrix<double, 6, 1> l3d1 = v3Lines[1];
	Matrix<double, 6, 1> l3d2 = v3Lines[2];

	Vector3d l3p1 = l3d0.block<3, 1>(0, 0);
	Vector3d l3p2 = l3d0.block<3, 1>(3, 0);

	Vector3d m0 = l3p2 - l3p1;
	m0.normalize();

	l3p1 = Vector3d(l3d1.block<3, 1>(0, 0));
	l3p2 = Vector3d(l3d1.block<3, 1>(3, 0));
	Vector3d m1 = l3p2 - l3p1;
	m1.normalize();

	Vector3d e1 = m0.cross(m1);
	Vector3d n1 = e1.cross(m0);
	Rmw.block<1, 3>(0, 0) = n1.transpose();
	Rmw.block<1, 3>(1, 0) = e1.transpose();
	Rmw.block<1, 3>(2, 0) = m0.transpose();

	return true;
}

////////////////////////////////////////////////////////////////////////////
bool PnL::generateRotMat(const Vector3d &n1, const Vector3d &m1, Matrix3d &R_)
{
	//get a vector e1 orthogonal to both n1 and m1
	Vector3d e1 = n1.cross(m1);
	e1.normalize();
	//std::cout << "check e1: " << e1 << std::endl;

	//generate a second vector m1_ as the n1xe1
	Vector3d m1_ = n1.cross(e1);
	m1_.normalize();

	//construct the rotation matrix
	Matrix3d R;
	R.block<3, 1>(0, 0) = n1;
	R.block<3, 1>(0, 1) = e1;
	R.block<3, 1>(0, 2) = m1_;

	//std::cout << "check R inside: " << R << std::endl;

	//check if determinant of R is +1
	double det = R.determinant();
	if (det < 0)
	{
		R_.block<3, 1>(0, 0) = R.block<3, 1>(0, 0);
		R_.block<3, 1>(0, 1) = R.block<3, 1>(0, 2);
		R_.block<3, 1>(0, 2) = R.block<3, 1>(0, 1);
		return true;
	}
	
	R_ = R;
	return true;
}


///////////////////////////////////////////////////////////////////
void PnL::constructPolynomCoeffs(const Matrix3d &R_, const Vector3d &n1_,
	const Vector3d &n2_, const Vector3d &v1_,
	const Vector3d &v2_)
{
	Matrix<double, 1, 3> n1 = n1_.transpose()*R_;
	Matrix<double, 1, 3> n2 = n2_.transpose()*R_;
	
	//! simplify the elements of the vectors to create the coefficients. will not need later
	double nx1 = n1(0,0); double ny1 = n1(0,1); double nz1 = n1(0,2);
	double nx2 = n2(0,0); double ny2 = n2(0,1); double nz2 = n2(0,2);
	double vx1 = v1_[0]; double vy1 = v1_[1]; double vz1 = v1_[2];
	double vx2 = v2_[0]; double vy2 = v2_[1]; double vz2 = v2_[2];

	double u4 = (2 * pow(ny1, 2) * pow(ny2, 2) * vx1*vx2*vy1*vy2 - pow(ny1, 2) * pow(ny2, 2) * pow(vx1, 2) * pow(vz2, 2) - pow(ny1, 2) * pow(ny2, 2) *
		pow(vx1, 2) * pow(vy2, 2) + 2 * pow(ny1, 2) * pow(ny2, 2) * vx1*vx2*vz1*vz2 - pow(ny1, 2) * pow(ny2, 2) * pow(vx2, 2) * pow(vy1, 2) - pow(ny1, 2) *
		pow(ny2, 2) * pow(vx2, 2) * pow(vz1, 2) - pow(ny1, 2) * pow(ny2, 2) * pow(vy1, 2) * pow(vz2, 2) +
		2 * pow(ny1, 2) * pow(ny2, 2) * vy1*vy2*vz1*vz2 - pow(ny1, 2) * pow(ny2, 2) * pow(vy2, 2) * pow(vz1, 2) +
		pow(ny1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2) + pow(ny1, 2) * pow(nz2, 2) * pow(vx1, 2) *
		pow(vz2, 2) - 2 * pow(ny1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2 - 2 * pow(ny1, 2) * pow(nz2, 2) *
		vx1*vx2*vz1*vz2 + pow(ny1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2) + pow(ny1, 2) * pow(nz2, 2) *
		pow(vx2, 2) * pow(vz1, 2) + pow(ny1, 2) * pow(nz2, 2) * pow(vy1, 2) * pow(vz2, 2) - 2 * pow(ny1, 2) *
		pow(nz2, 2) * vy1*vy2*vz1*vz2 + pow(ny1, 2) * pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2) +
		4 * ny1*ny2*nz1*nz2*pow(vx1, 2) * pow(vy2, 2) + 4 * ny1*ny2*nz1*nz2*pow(vx1, 2) * pow(vz2, 2) -
		8 * ny1*ny2*nz1*nz2*vx1*vx2*vy1*vy2 - 8 * ny1*ny2*nz1*nz2*vx1*vx2*vz1*vz2 +
		4 * ny1*ny2*nz1*nz2*pow(vx2, 2) * pow(vy1, 2) + 4 * ny1*ny2*nz1*nz2*pow(vx2, 2) * pow(vz1, 2)
		+ 4 * ny1*ny2*nz1*nz2*pow(vy1, 2) * pow(vz2, 2) - 8 * ny1*ny2*nz1*nz2*vy1*vy2*vz1*vz2
		+ 4 * ny1*ny2*nz1*nz2*pow(vy2, 2) * pow(vz1, 2) + pow(ny2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vy2, 2)
		+ pow(ny2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vz2, 2) - 2 * pow(ny2, 2) * pow(nz1, 2) * vx1*vx2*vy1*vy2
		- 2 * pow(ny2, 2) * pow(nz1, 2) * vx1*vx2*vz1*vz2 + pow(ny2, 2) * pow(nz1, 2) * pow(vx2, 2) * pow(vy1, 2)
		+ pow(ny2, 2) * pow(nz1, 2) * pow(vx2, 2) * pow(vz1, 2) + pow(ny2, 2) * pow(nz1, 2) * pow(vy1, 2)
		* pow(vz2, 2) - 2 * pow(ny2, 2) * pow(nz1, 2) * vy1*vy2*vz1*vz2 + pow(ny2, 2) * pow(nz1, 2)
		* pow(vy2, 2) * pow(vz1, 2) - pow(nz1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2) - pow(nz1, 2)
		* pow(nz2, 2) * pow(vx1, 2) * pow(vz2, 2) + 2 * pow(nz1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2 +
		2 * pow(nz1, 2) * pow(nz2, 2) * vx1*vx2*vz1*vz2 - pow(nz1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2)
		- pow(nz1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vz1, 2) - pow(nz1, 2) * pow(nz2, 2) * pow(vy1, 2)
		* pow(vz2, 2) + 2 * pow(nz1, 2) * pow(nz2, 2) * vy1*vy2*vz1*vz2 - pow(nz1, 2) * pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2));

	double v3 = (4 * pow(ny1, 2) * ny2*nz2*vx1*vx2*vy1*vy2 - 2 * pow(ny1, 2) * ny2*nz2*pow(vx1, 2) * pow(vz2, 2)
		- 2 * pow(ny1, 2) * ny2*nz2*pow(vx1, 2) * pow(vy2, 2) + 4 * pow(ny1, 2) * ny2*nz2*vx1*vx2*vz1*vz2
		- 2 * pow(ny1, 2) * ny2*nz2*pow(vx2, 2) * pow(vy1, 2) - 2 * pow(ny1, 2) * ny2*nz2*pow(vx2, 2) * pow(vz1, 2)
		- 2 * pow(ny1, 2) * ny2*nz2*pow(vy1, 2) * pow(vz2, 2) + 4 * pow(ny1, 2) * ny2*nz2*vy1*vy2*vz1*vz2
		- 2 * pow(ny1, 2) * ny2*nz2*pow(vy2, 2) * pow(vz1, 2) - 2 * ny1*pow(ny2, 2) * nz1*pow(vx1, 2) * pow(vy2, 2)
		- 2 * ny1*pow(ny2, 2) * nz1*pow(vx1, 2) * pow(vz2, 2) + 4 * ny1*pow(ny2, 2) * nz1*vx1*vx2*vy1*vy2
		+ 4 * ny1*pow(ny2, 2) * nz1*vx1*vx2*vz1*vz2 - 2 * ny1*pow(ny2, 2) * nz1*pow(vx2, 2) * pow(vy1, 2)
		- 2 * ny1*pow(ny2, 2) * nz1*pow(vx2, 2) * pow(vz1, 2) - 2 * ny1*pow(ny2, 2) * nz1*pow(vy1, 2) * pow(vz2, 2)
		+ 4 * ny1*pow(ny2, 2) * nz1*vy1*vy2*vz1*vz2 - 2 * ny1*pow(ny2, 2) * nz1*pow(vy2, 2) * pow(vz1, 2)
		+ 2 * ny1*nz1*pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2) + 2 * ny1*nz1*pow(nz2, 2) * pow(vx1, 2) * pow(vz2, 2)
		- 4 * ny1*nz1*pow(nz2, 2) * vx1*vx2*vy1*vy2 - 4 * ny1*nz1*pow(nz2, 2) * vx1*vx2*vz1*vz2
		+ 2 * ny1*nz1*pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2) + 2 * ny1*nz1*pow(nz2, 2) * pow(vx2, 2) * pow(vz1, 2)
		+ 2 * ny1*nz1*pow(nz2, 2) * pow(vy1, 2) * pow(vz2, 2) - 4 * ny1*nz1*pow(nz2, 2) * vy1*vy2*vz1*vz2
		+ 2 * ny1*nz1*pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2) + 2 * ny2*pow(nz1, 2) * nz2*pow(vx1, 2) * pow(vy2, 2)
		+ 2 * ny2*pow(nz1, 2) * nz2*pow(vx1, 2) * pow(vz2, 2) - 4 * ny2*pow(nz1, 2) * nz2*vx1*vx2*vy1*vy2
		- 4 * ny2*pow(nz1, 2) * nz2*vx1*vx2*vz1*vz2 + 2 * ny2*pow(nz1, 2) * nz2*pow(vx2, 2) * pow(vy1, 2)
		+ 2 * ny2*pow(nz1, 2) * nz2*pow(vx2, 2) * pow(vz1, 2) + 2 * ny2*pow(nz1, 2) * nz2*pow(vy1, 2) * pow(vz2, 2)
		- 4 * ny2*pow(nz1, 2) * nz2*vy1*vy2*vz1*vz2 + 2 * ny2*pow(nz1, 2) * nz2*pow(vy2, 2) * pow(vz1, 2));


	double u3 = (2 * nx2*pow(ny1, 2) * ny2*vx1*pow(vx2, 2) * vy1 - 2 * nx2*pow(ny1, 2) * ny2*pow(vx1, 2) * vx2*vy2
		- 2 * nx2*pow(ny1, 2) * ny2*vx1*vy1*pow(vy2, 2) - 2 * nx2*vz1*vz2*pow(ny1, 2) * ny2*vx1*vy2
		+ 2 * nx2*pow(ny1, 2) * ny2*vx2*pow(vy1, 2) * vy2 + 2 * nx2*vz1*vz2*pow(ny1, 2) * ny2*vx2*vy1
		+ 2 * nx1*ny1*pow(ny2, 2) * pow(vx1, 2) * vx2*vy2 - 2 * nx1*ny1*pow(ny2, 2) * vx1*pow(vx2, 2) * vy1
		+ 2 * nx1*ny1*pow(ny2, 2) * vx1*vy1*pow(vy2, 2) + 2 * nx1*vz1*vz2*ny1*pow(ny2, 2) * vx1*vy2
		- 2 * nx1*ny1*pow(ny2, 2) * vx2*pow(vy1, 2) * vy2 - 2 * nx1*vz1*vz2*ny1*pow(ny2, 2) * vx2*vy1
		+ 4 * nx2*ny1*nz1*nz2*pow(vx1, 2) * vx2*vy2 - 4 * nx2*ny1*nz1*nz2*vx1*pow(vx2, 2) * vy1
		+ 4 * nx2*ny1*nz1*nz2*vx1*vy1*pow(vy2, 2) + 4 * nx2*vz1*vz2*ny1*nz1*nz2*vx1*vy2
		- 4 * nx2*ny1*nz1*nz2*vx2*pow(vy1, 2) * vy2 - 4 * nx2*vz1*vz2*ny1*nz1*nz2*vx2*vy1
		- 2 * nx1*ny1*pow(nz2, 2) * pow(vx1, 2) * vx2*vy2 + 2 * nx1*ny1*pow(nz2, 2) * vx1*pow(vx2, 2) * vy1
		- 2 * nx1*ny1*pow(nz2, 2) * vx1*vy1*pow(vy2, 2) - 2 * nx1*vz1*vz2*ny1*pow(nz2, 2) * vx1*vy2
		+ 2 * nx1*ny1*pow(nz2, 2) * vx2*pow(vy1, 2) * vy2 + 2 * nx1*vz1*vz2*ny1*pow(nz2, 2) * vx2*vy1
		+ 2 * nx2*ny2*pow(nz1, 2) * pow(vx1, 2) * vx2*vy2 - 2 * nx2*ny2*pow(nz1, 2) * vx1*pow(vx2, 2) * vy1
		+ 2 * nx2*ny2*pow(nz1, 2) * vx1*vy1*pow(vy2, 2) + 2 * nx2*vz1*vz2*ny2*pow(nz1, 2) * vx1*vy2
		- 2 * nx2*ny2*pow(nz1, 2) * vx2*pow(vy1, 2) * vy2 - 2 * nx2*vz1*vz2*ny2*pow(nz1, 2) * vx2*vy1
		- 4 * nx1*ny2*nz1*nz2*pow(vx1, 2) * vx2*vy2 + 4 * nx1*ny2*nz1*nz2*vx1*pow(vx2, 2) * vy1
		- 4 * nx1*ny2*nz1*nz2*vx1*vy1*pow(vy2, 2) - 4 * nx1*vz1*vz2*ny2*nz1*nz2*vx1*vy2
		+ 4 * nx1*ny2*nz1*nz2*vx2*pow(vy1, 2) * vy2 + 4 * nx1*vz1*vz2*ny2*nz1*nz2*vx2*vy1);

	double v2 = (2 * nx2*pow(ny1, 2) * nz2*vx1*pow(vx2, 2) * vy1
		- 2 * nx2*pow(ny1, 2) * nz2*pow(vx1, 2) * vx2*vy2 - 2 * nx2*pow(ny1, 2) * nz2*vx1*vy1*pow(vy2, 2)
		- 2 * nx2*vz1*vz2*pow(ny1, 2) * nz2*vx1*vy2 + 2 * nx2*pow(ny1, 2) * nz2*vx2*pow(vy1, 2) * vy2
		+ 2 * nx2*vz1*vz2*pow(ny1, 2) * nz2*vx2*vy1 - 4 * nx2*ny1*ny2*nz1*pow(vx1, 2) * vx2*vy2
		+ 4 * nx2*ny1*ny2*nz1*vx1*pow(vx2, 2) * vy1 - 4 * nx2*ny1*ny2*nz1*vx1*vy1*pow(vy2, 2)
		- 4 * nx2*vz1*vz2*ny1*ny2*nz1*vx1*vy2 + 4 * nx2*ny1*ny2*nz1*vx2*pow(vy1, 2) * vy2
		+ 4 * nx2*vz1*vz2*ny1*ny2*nz1*vx2*vy1 + 4 * nx1*ny1*ny2*nz2*pow(vx1, 2) * vx2*vy2
		- 4 * nx1*ny1*ny2*nz2*vx1*pow(vx2, 2) * vy1 + 4 * nx1*ny1*ny2*nz2*vx1*vy1*pow(vy2, 2)
		+ 4 * nx1*vz1*vz2*ny1*ny2*nz2*vx1*vy2 - 4 * nx1*ny1*ny2*nz2*vx2*pow(vy1, 2) * vy2
		- 4 * nx1*vz1*vz2*ny1*ny2*nz2*vx2*vy1 + 2 * nx1*pow(ny2, 2) * nz1*pow(vx1, 2) * vx2*vy2
		- 2 * nx1*pow(ny2, 2) * nz1*vx1*pow(vx2, 2) * vy1 + 2 * nx1*pow(ny2, 2) * nz1*vx1*vy1*pow(vy2, 2)
		+ 2 * nx1*vz1*vz2*pow(ny2, 2) * nz1*vx1*vy2 - 2 * nx1*pow(ny2, 2) * nz1*vx2*pow(vy1, 2) * vy2
		- 2 * nx1*vz1*vz2*pow(ny2, 2) * nz1*vx2*vy1 + 2 * nx2*pow(nz1, 2) * nz2*pow(vx1, 2) * vx2*vy2
		- 2 * nx2*pow(nz1, 2) * nz2*vx1*pow(vx2, 2) * vy1 + 2 * nx2*pow(nz1, 2) * nz2*vx1*vy1*pow(vy2, 2)
		+ 2 * nx2*vz1*vz2*pow(nz1, 2) * nz2*vx1*vy2 - 2 * nx2*pow(nz1, 2) * nz2*vx2*pow(vy1, 2) * vy2
		- 2 * nx2*vz1*vz2*pow(nz1, 2) * nz2*vx2*vy1 - 2 * nx1*nz1*pow(nz2, 2) * pow(vx1, 2) * vx2*vy2
		+ 2 * nx1*nz1*pow(nz2, 2) * vx1*pow(vx2, 2) * vy1 - 2 * nx1*nz1*pow(nz2, 2) * vx1*vy1*pow(vy2, 2)
		- 2 * nx1*vz1*vz2*nz1*pow(nz2, 2) * vx1*vy2 + 2 * nx1*nz1*pow(nz2, 2) * vx2*pow(vy1, 2) * vy2
		+ 2 * nx1*vz1*vz2*nz1*pow(nz2, 2) * vx2*vy1);

	double u2 = (pow(nx1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vx2, 2) - pow(nx1, 2) * pow(ny2, 2) * pow(vx1, 2) * pow(vz2, 2)
		- 2 * pow(nx1, 2) * pow(ny2, 2) * vx1*vx2*vy1*vy2 - pow(nx1, 2) * pow(ny2, 2) * pow(vy1, 2) * pow(vy2, 2)
		- pow(nx1, 2) * pow(ny2, 2) * pow(vy1, 2) * pow(vz2, 2) - pow(nx1, 2) * pow(ny2, 2) * pow(vx1, 2) * pow(vx2, 2)
		+ pow(nx1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vz2, 2) + 2 * pow(nx1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2
		+ pow(nx1, 2) * pow(nz2, 2) * pow(vy1, 2) * pow(vy2, 2) + pow(nx1, 2) * pow(nz2, 2) * pow(vy1, 2) * pow(vz2, 2)
		+ 2 * nx1*nx2*ny1*ny2*pow(vx1, 2) * pow(vx2, 2) - 2 * nx1*nx2*ny1*ny2*pow(vx1, 2) * pow(vy2, 2)
		+ 8 * nx1*nx2*ny1*ny2*vx1*vx2*vy1*vy2 + 2 * nx1*nx2*ny1*ny2*vx1*vx2*vz1*vz2
		- 2 * nx1*nx2*ny1*ny2*pow(vx2, 2) * pow(vy1, 2) + 2 * nx1*nx2*ny1*ny2*pow(vy1, 2) * pow(vy2, 2)
		+ 2 * nx1*nx2*ny1*ny2*vy1*vy2*vz1*vz2 - 2 * nx1*nx2*nz1*nz2*pow(vx1, 2) * pow(vx2, 2)
		+ 2 * nx1*nx2*nz1*nz2*pow(vx1, 2) * pow(vy2, 2) - 8 * nx1*nx2*nz1*nz2*vx1*vx2*vy1*vy2
		- 2 * nx1*nx2*nz1*nz2*vx1*vx2*vz1*vz2 + 2 * nx1*nx2*nz1*nz2*pow(vx2, 2) * pow(vy1, 2)
		- 2 * nx1*nx2*nz1*nz2*pow(vy1, 2) * pow(vy2, 2) - 2 * nx1*nx2*nz1*nz2*vy1*vy2*vz1*vz2
		- pow(nx2, 2) * pow(ny1, 2) * pow(vx1, 2) * pow(vx2, 2) - 2 * pow(nx2, 2) * pow(ny1, 2) * vx1*vx2*vy1*vy2
		- pow(nx2, 2) * pow(ny1, 2) * pow(vx2, 2) * pow(vz1, 2) - pow(nx2, 2) * pow(ny1, 2) * pow(vy1, 2) * pow(vy2, 2)
		- pow(nx2, 2) * pow(ny1, 2) * pow(vy2, 2) * pow(vz1, 2) + pow(nx2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vx2, 2)
		+ 2 * pow(nx2, 2) * pow(nz1, 2) * vx1*vx2*vy1*vy2 + pow(nx2, 2) * pow(nz1, 2) * pow(vx2, 2) * pow(vz1, 2)
		+ pow(nx2, 2) * pow(nz1, 2) * pow(vy1, 2) * pow(vy2, 2) + pow(nx2, 2) * pow(nz1, 2) * pow(vy2, 2) * pow(vz1, 2)
		+ pow(ny1, 2) * pow(ny2, 2) * pow(vx1, 2) * pow(vz2, 2) - 2 * pow(ny1, 2) * pow(ny2, 2) * vx1*vx2*vz1*vz2
		+ pow(ny1, 2) * pow(ny2, 2) * pow(vx2, 2) * pow(vz1, 2) + pow(ny1, 2) * pow(ny2, 2) * pow(vy1, 2) * pow(vz2, 2)
		- 2 * pow(ny1, 2) * pow(ny2, 2) * vy1*vy2*vz1*vz2 + pow(ny1, 2) * pow(ny2, 2) * pow(vy2, 2) * pow(vz1, 2)
		- pow(ny1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2) + 2 * pow(ny1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2
		+ 2 * pow(ny1, 2) * pow(nz2, 2) * vx1*vx2*vz1*vz2 - pow(ny1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2)
		- 2 * pow(ny1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vz1, 2) + 2 * pow(ny1, 2) * pow(nz2, 2) * vy1*vy2*vz1*vz2
		- 2 * pow(ny1, 2) * pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2) - 4 * ny1*ny2*nz1*nz2*pow(vx1, 2) * pow(vy2, 2)
		- 4 * ny1*ny2*nz1*nz2*pow(vx1, 2) * pow(vz2, 2) + 8 * ny1*ny2*nz1*nz2*vx1*vx2*vy1*vy2 + 8 * ny1*ny2*nz1*nz2*vx1*vx2*vz1*vz2
		- 4 * ny1*ny2*nz1*nz2*pow(vx2, 2) * pow(vy1, 2) - 4 * ny1*ny2*nz1*nz2*pow(vx2, 2) * pow(vz1, 2)
		- 4 * ny1*ny2*nz1*nz2*pow(vy1, 2) * pow(vz2, 2) + 8 * ny1*ny2*nz1*nz2*vy1*vy2*vz1*vz2
		- 4 * ny1*ny2*nz1*nz2*pow(vy2, 2) * pow(vz1, 2) - pow(ny2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vy2, 2)
		- 2 * pow(ny2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vz2, 2) + 2 * pow(ny2, 2) * pow(nz1, 2) * vx1*vx2*vy1*vy2
		+ 2 * pow(ny2, 2) * pow(nz1, 2) * vx1*vx2*vz1*vz2 - pow(ny2, 2) * pow(nz1, 2) * pow(vx2, 2) * pow(vy1, 2)
		- 2 * pow(ny2, 2) * pow(nz1, 2) * pow(vy1, 2) * pow(vz2, 2) + 2 * pow(ny2, 2) * pow(nz1, 2) * vy1*vy2*vz1*vz2
		+ 2 * pow(nz1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2) + pow(nz1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vz2, 2)
		- 4 * pow(nz1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2 - 2 * pow(nz1, 2) * pow(nz2, 2) * vx1*vx2*vz1*vz2
		+ 2 * pow(nz1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2) + pow(nz1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vz1, 2)
		+ pow(nz1, 2) * pow(nz2, 2) * pow(vy1, 2) * pow(vz2, 2) - 2 * pow(nz1, 2) * pow(nz2, 2) * vy1*vy2*vz1*vz2
		+ pow(nz1, 2) * pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2));

	double v1 = (2 * nx1*nx2*ny1*nz2*pow(vx1, 2) * pow(vx2, 2) - 2 * pow(nx1, 2) * ny2*nz2*pow(vx1, 2) * pow(vz2, 2)
		- 4 * pow(nx1, 2) * ny2*nz2*vx1*vx2*vy1*vy2 - 2 * pow(nx1, 2) * ny2*nz2*pow(vy1, 2) * pow(vy2, 2)
		- 2 * pow(nx1, 2) * ny2*nz2*pow(vy1, 2) * pow(vz2, 2) - 2 * pow(nx1, 2) * ny2*nz2*pow(vx1, 2) * pow(vx2, 2)
		- 2 * nx1*nx2*ny1*nz2*pow(vx1, 2) * pow(vy2, 2) + 8 * nx1*nx2*ny1*nz2*vx1*vx2*vy1*vy2 + 2 * nx1*nx2*ny1*nz2*vx1*vx2*vz1*vz2
		- 2 * nx1*nx2*ny1*nz2*pow(vx2, 2) * pow(vy1, 2) + 2 * nx1*nx2*ny1*nz2*pow(vy1, 2) * pow(vy2, 2)
		+ 2 * nx1*nx2*ny1*nz2*vy1*vy2*vz1*vz2 + 2 * nx1*nx2*ny2*nz1*pow(vx1, 2) * pow(vx2, 2)
		- 2 * nx1*nx2*ny2*nz1*pow(vx1, 2) * pow(vy2, 2) + 8 * nx1*nx2*ny2*nz1*vx1*vx2*vy1*vy2
		+ 2 * nx1*nx2*ny2*nz1*vx1*vx2*vz1*vz2 - 2 * nx1*nx2*ny2*nz1*pow(vx2, 2) * pow(vy1, 2)
		+ 2 * nx1*nx2*ny2*nz1*pow(vy1, 2) * pow(vy2, 2) + 2 * nx1*nx2*ny2*nz1*vy1*vy2*vz1*vz2
		- 2 * pow(nx2, 2) * ny1*nz1*pow(vx1, 2) * pow(vx2, 2) - 4 * pow(nx2, 2) * ny1*nz1*vx1*vx2*vy1*vy2
		- 2 * pow(nx2, 2) * ny1*nz1*pow(vx2, 2) * pow(vz1, 2) - 2 * pow(nx2, 2) * ny1*nz1*pow(vy1, 2) * pow(vy2, 2)
		- 2 * pow(nx2, 2) * ny1*nz1*pow(vy2, 2) * pow(vz1, 2) - 2 * pow(ny1, 2) * ny2*nz2*vx1*vx2*vz1*vz2
		+ 2 * pow(ny1, 2) * ny2*nz2*pow(vx2, 2) * pow(vz1, 2) - 2 * pow(ny1, 2) * ny2*nz2*vy1*vy2*vz1*vz2
		+ 2 * pow(ny1, 2) * ny2*nz2*pow(vy2, 2) * pow(vz1, 2) + 2 * ny1*pow(ny2, 2) * nz1*pow(vx1, 2) * pow(vz2, 2)
		- 2 * ny1*pow(ny2, 2) * nz1*vx1*vx2*vz1*vz2 + 2 * ny1*pow(ny2, 2) * nz1*pow(vy1, 2) * pow(vz2, 2)
		- 2 * ny1*pow(ny2, 2) * nz1*vy1*vy2*vz1*vz2 - 2 * ny1*nz1*pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2)
		+ 4 * ny1*nz1*pow(nz2, 2) * vx1*vx2*vy1*vy2 + 2 * ny1*nz1*pow(nz2, 2) * vx1*vx2*vz1*vz2
		- 2 * ny1*nz1*pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2) - 2 * ny1*nz1*pow(nz2, 2) * pow(vx2, 2) * pow(vz1, 2)
		+ 2 * ny1*nz1*pow(nz2, 2) * vy1*vy2*vz1*vz2 - 2 * ny1*nz1*pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2)
		- 2 * ny2*pow(nz1, 2) * nz2*pow(vx1, 2) * pow(vy2, 2) - 2 * ny2*pow(nz1, 2) * nz2*pow(vx1, 2) * pow(vz2, 2)
		+ 4 * ny2*pow(nz1, 2) * nz2*vx1*vx2*vy1*vy2 + 2 * ny2*pow(nz1, 2) * nz2*vx1*vx2*vz1*vz2
		- 2 * ny2*pow(nz1, 2) * nz2*pow(vx2, 2) * pow(vy1, 2) - 2 * ny2*pow(nz1, 2) * nz2*pow(vy1, 2) * pow(vz2, 2)
		+ 2 * ny2*pow(nz1, 2) * nz2*vy1*vy2*vz1*vz2);

	double u1 = (2 * pow(nx1, 2) * nx2*ny2*pow(vx1, 2) * vx2*vy2 - 2 * pow(nx1, 2) * nx2*ny2*vx1*pow(vx2, 2) * vy1
		+ 2 * pow(nx1, 2) * nx2*ny2*vx1*vy1*pow(vy2, 2) - 2 * pow(nx1, 2) * nx2*ny2*vx2*pow(vy1, 2) * vy2
		- 2 * nx1*pow(nx2, 2) * ny1*pow(vx1, 2) * vx2*vy2 + 2 * nx1*pow(nx2, 2) * ny1*vx1*pow(vx2, 2) * vy1
		- 2 * nx1*pow(nx2, 2) * ny1*vx1*vy1*pow(vy2, 2) + 2 * nx1*pow(nx2, 2) * ny1*vx2*pow(vy1, 2) * vy2
		- 2 * vz1*vz2*nx1*ny1*pow(ny2, 2) * vx1*vy2 + 2 * vz1*vz2*nx1*ny1*pow(ny2, 2) * vx2*vy1
		+ 2 * nx1*ny1*pow(nz2, 2) * pow(vx1, 2) * vx2*vy2 - 2 * nx1*ny1*pow(nz2, 2) * vx1*pow(vx2, 2) * vy1
		+ 2 * nx1*ny1*pow(nz2, 2) * vx1*vy1*pow(vy2, 2) + 2 * vz1*vz2*nx1*ny1*pow(nz2, 2) * vx1*vy2
		- 2 * nx1*ny1*pow(nz2, 2) * vx2*pow(vy1, 2) * vy2 - 2 * vz1*vz2*nx1*ny1*pow(nz2, 2) * vx2*vy1
		+ 4 * nx1*ny2*nz1*nz2*pow(vx1, 2) * vx2*vy2 - 4 * nx1*ny2*nz1*nz2*vx1*pow(vx2, 2) * vy1
		+ 4 * nx1*ny2*nz1*nz2*vx1*vy1*pow(vy2, 2) + 2 * vz1*vz2*nx1*ny2*nz1*nz2*vx1*vy2
		- 4 * nx1*ny2*nz1*nz2*vx2*pow(vy1, 2) * vy2 - 2 * vz1*vz2*nx1*ny2*nz1*nz2*vx2*vy1
		+ 2 * vz1*vz2*nx2*pow(ny1, 2) * ny2*vx1*vy2 - 2 * vz1*vz2*nx2*pow(ny1, 2) * ny2*vx2*vy1
		- 4 * nx2*ny1*nz1*nz2*pow(vx1, 2) * vx2*vy2 + 4 * nx2*ny1*nz1*nz2*vx1*pow(vx2, 2) * vy1
		- 4 * nx2*ny1*nz1*nz2*vx1*vy1*pow(vy2, 2) - 2 * vz1*vz2*nx2*ny1*nz1*nz2*vx1*vy2
		+ 4 * nx2*ny1*nz1*nz2*vx2*pow(vy1, 2) * vy2 + 2 * vz1*vz2*nx2*ny1*nz1*nz2*vx2*vy1
		- 2 * nx2*ny2*pow(nz1, 2) * pow(vx1, 2) * vx2*vy2 + 2 * nx2*ny2*pow(nz1, 2) * vx1*pow(vx2, 2) * vy1
		- 2 * nx2*ny2*pow(nz1, 2) * vx1*vy1*pow(vy2, 2) - 2 * vz1*vz2*nx2*ny2*pow(nz1, 2) * vx1*vy2
		+ 2 * nx2*ny2*pow(nz1, 2) * vx2*pow(vy1, 2) * vy2 + 2 * vz1*vz2*nx2*ny2*pow(nz1, 2) * vx2*vy1);

	double v0 = (2 * pow(nx1, 2) * nx2*nz2*pow(vx1, 2) * vx2*vy2 - 2 * pow(nx1, 2) * nx2*nz2*vx1*pow(vx2, 2) * vy1
		+ 2 * pow(nx1, 2) * nx2*nz2*vx1*vy1*pow(vy2, 2) - 2 * pow(nx1, 2) * nx2*nz2*vx2*pow(vy1, 2) * vy2
		- 2 * nx1*pow(nx2, 2) * nz1*pow(vx1, 2) * vx2*vy2 + 2 * nx1*pow(nx2, 2) * nz1*vx1*pow(vx2, 2) * vy1
		- 2 * nx1*pow(nx2, 2) * nz1*vx1*vy1*pow(vy2, 2) + 2 * nx1*pow(nx2, 2) * nz1*vx2*pow(vy1, 2) * vy2
		+ 2 * nx1*nz1*pow(nz2, 2) * pow(vx1, 2) * vx2*vy2 - 2 * nx1*nz1*pow(nz2, 2) * vx1*pow(vx2, 2) * vy1
		+ 2 * nx1*nz1*pow(nz2, 2) * vx1*vy1*pow(vy2, 2) - 2 * nx1*nz1*pow(nz2, 2) * vx2*pow(vy1, 2) * vy2
		- 2 * ny1*ny2*vz1*vz2*nx1*nz2*vx1*vy2 + 2 * ny1*ny2*vz1*vz2*nx1*nz2*vx2*vy1
		- 2 * nx2*pow(nz1, 2) * nz2*pow(vx1, 2) * vx2*vy2 + 2 * nx2*pow(nz1, 2) * nz2*vx1*pow(vx2, 2) * vy1
		- 2 * nx2*pow(nz1, 2) * nz2*vx1*vy1*pow(vy2, 2) + 2 * nx2*pow(nz1, 2) * nz2*vx2*pow(vy1, 2) * vy2
		+ 2 * ny1*ny2*vz1*vz2*nx2*nz1*vx1*vy2 - 2 * ny1*ny2*vz1*vz2*nx2*nz1*vx2*vy1);

	double u0 = -pow(nx1, 2) * pow(nx2, 2) * pow(vx1, 2) * pow(vy2, 2) + 2 * pow(nx1, 2) * pow(nx2, 2) * vx1*vx2*vy1*vy2 - pow(nx1, 2) * pow(nx2, 2) * pow(vx2, 2) * pow(vy1, 2) + pow(nx1, 2) * pow(ny2, 2) * pow(vx1, 2) * pow(vz2, 2) + pow(nx1, 2) * pow(ny2, 2) * pow(vy1, 2) * pow(vz2, 2) - pow(nx1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vx2, 2) - 2 * pow(nx1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2 - pow(nx1, 2) * pow(nz2, 2) * pow(vy1, 2) * pow(vy2, 2) - 2 * nx1*nx2*ny1*ny2*vx1*vx2*vz1*vz2 - 2 * nx1*nx2*ny1*ny2*vy1*vy2*vz1*vz2 + 2 * nx1*nx2*nz1*nz2*pow(vx1, 2) * pow(vx2, 2) - 2 * nx1*nx2*nz1*nz2*pow(vx1, 2) * pow(vy2, 2) + 8 * nx1*nx2*nz1*nz2*vx1*vx2*vy1*vy2 - 2 * nx1*nx2*nz1*nz2*pow(vx2, 2) * pow(vy1, 2) + 2 * nx1*nx2*nz1*nz2*pow(vy1, 2) * pow(vy2, 2) + pow(nx2, 2) * pow(ny1, 2) * pow(vx2, 2) * pow(vz1, 2) + pow(nx2, 2) * pow(ny1, 2) * pow(vy2, 2) * pow(vz1, 2) - pow(nx2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vx2, 2) - 2 * pow(nx2, 2) * pow(nz1, 2) * vx1*vx2*vy1*vy2 - pow(nx2, 2) * pow(nz1, 2) * pow(vy1, 2) * pow(vy2, 2) + pow(ny1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vz1, 2) + pow(ny1, 2) * pow(nz2, 2) * pow(vy2, 2) * pow(vz1, 2) - 2 * ny1*ny2*nz1*nz2*vx1*vx2*vz1*vz2 - 2 * ny1*ny2*nz1*nz2*vy1*vy2*vz1*vz2 + pow(ny2, 2) * pow(nz1, 2) * pow(vx1, 2) * pow(vz2, 2) + pow(ny2, 2) * pow(nz1, 2) * pow(vy1, 2) * pow(vz2, 2) - pow(nz1, 2) * pow(nz2, 2) * pow(vx1, 2) * pow(vy2, 2) + 2 * pow(nz1, 2) * pow(nz2, 2) * vx1*vx2*vy1*vy2 - pow(nz1, 2) * pow(nz2, 2) * pow(vx2, 2) * pow(vy1, 2);

	v0 = -v0;
	v1 = -v1;
	v2 = -v2;
	v3 = -v3;


	coeffVector[8] = u4*u4 + v3*v3;
	coeffVector[7] = 2 * u3*u4 + 2 * v2*v3;
	coeffVector[6] = 2 * u2*u4 + u3*u3 + 2 * v1*v3 + v2*v2 - v3*v3;
	coeffVector[5] = 2 * u1*u4 + 2 * u2*u3 + 2 * v0*v3 + 2 * v1*v2 - 2 * v2*v3;
	coeffVector[4] = 2 * u0*u4 + 2 * u1*u3 + u2*u2 + 2 * v0*v2 + v1*v1 - 2 * v1*v3 - v2*v2;
	coeffVector[3] = 2 * u0*u3 + 2 * u1*u2 + 2 * v0*v1 - 2 * v0*v3 - 2 * v1*v2;
	coeffVector[2] = 2 * u0*u2 + u1*u1 + v0*v0 - 2 * v0*v2 - v1*v1;
	coeffVector[1] = 2 * u0*u1 - 2 * v0*v1;
	coeffVector[0] = u0*u0 - v0*v0;

	std::cout << "Polynomial: ";
	for (int i = 0; i<8; ++i) { std::cout << coeffVector[i] << ".x^" << i << "+ "; }
	std::cout << coeffVector[8] << ".x^8" << std::endl;

}


//////////////////////////////////////////////////////////
bool PnL::computeRoots()
{
	if (coeffVector.norm() == 0)
		return false;
	pSolver.compute(coeffVector);
	Eigen::Matrix<std::complex<double>, 8, 1>roots = pSolver.roots();
	
	if (realRoots.size() > 0)
		realRoots.resize(0);

	pSolver.realRoots(realRoots, 0.00001);

	for (size_t i = 0; i < realRoots.size(); ++i)
	{
		std::cout << "check real roots: " << realRoots[i] << std::endl;
	}
	return true;
}


///////////////////////////////////////////////////////////////////
bool PnL::computeTranslation(const Matrix3d &R, const Vector3d &n0,
	const Vector3d &n1, const Vector3d &n2,
	const Vector3d &P0, const Vector3d &P1,
	const Vector3d &P2, Vector3d &t)
{
	//construct the matrix N from the three normals
	Matrix3d N;
	N.block<1, 3>(0, 0) = n0.transpose();
	N.block<1, 3>(1, 0) = n1.transpose();
	N.block<1, 3>(2, 0) = n2.transpose();

	Eigen::FullPivLU<Eigen::Matrix<double, 3, 3> > lu_decomp(N);
	int rank = lu_decomp.rank();
	if (rank<3)   //three lines are not linearly independent
		return false;

	double t1 = -(P0[0] * n0[0] * n1[1] * n2[2] * R(0, 0) - P0[0] * n0[0] * n1[2] * n2[1] * R(0, 0) + P0[1] * n0[0] * n1[1] * n2[2] * R(0, 1)
		- P0[1] * n0[0] * n1[2] * n2[1] * R(0, 1) + P0[2] * n0[0] * n1[1] * n2[2] * R(0, 2) - P0[2] * n0[0] * n1[2] * n2[1] * R(0, 2)
		- P1[0] * n0[1] * n1[0] * n2[2] * R(0, 0) + P1[0] * n0[2] * n1[0] * n2[1] * R(0, 0) + P0[0] * n0[1] * n1[1] * n2[2] * R(1, 0)
		- P0[0] * n0[1] * n1[2] * n2[1] * R(1, 0) - P1[1] * n0[1] * n1[0] * n2[2] * R(0, 1) + P1[1] * n0[2] * n1[0] * n2[1] * R(0, 1)
		+ P0[1] * n0[1] * n1[1] * n2[2] * R(1, 1) - P0[1] * n0[1] * n1[2] * n2[1] * R(1, 1) - P1[2] * n0[1] * n1[0] * n2[2] * R(0, 2)
		+ P1[2] * n0[2] * n1[0] * n2[1] * R(0, 2) + P0[2] * n0[1] * n1[1] * n2[2] * R(1, 2) - P0[2] * n0[1] * n1[2] * n2[1] * R(1, 2)
		+ P2[0] * n0[1] * n1[2] * n2[0] * R(0, 0) - P2[0] * n0[2] * n1[1] * n2[0] * R(0, 0) - P1[0] * n0[1] * n1[1] * n2[2] * R(1, 0)
		+ P1[0] * n0[2] * n1[1] * n2[1] * R(1, 0) + P0[0] * n0[2] * n1[1] * n2[2] * R(2, 0) - P0[0] * n0[2] * n1[2] * n2[1] * R(2, 0)
		+ P2[1] * n0[1] * n1[2] * n2[0] * R(0, 1) - P2[1] * n0[2] * n1[1] * n2[0] * R(0, 1) - P1[1] * n0[1] * n1[1] * n2[2] * R(1, 1)
		+ P1[1] * n0[2] * n1[1] * n2[1] * R(1, 1) + P0[1] * n0[2] * n1[1] * n2[2] * R(2, 1) - P0[1] * n0[2] * n1[2] * n2[1] * R(2, 1)
		+ P2[2] * n0[1] * n1[2] * n2[0] * R(0, 2) - P2[2] * n0[2] * n1[1] * n2[0] * R(0, 2) - P1[2] * n0[1] * n1[1] * n2[2] * R(1, 2)
		+ P1[2] * n0[2] * n1[1] * n2[1] * R(1, 2) + P0[2] * n0[2] * n1[1] * n2[2] * R(2, 2) - P0[2] * n0[2] * n1[2] * n2[1] * R(2, 2)
		+ P2[0] * n0[1] * n1[2] * n2[1] * R(1, 0) - P2[0] * n0[2] * n1[1] * n2[1] * R(1, 0) - P1[0] * n0[1] * n1[2] * n2[2] * R(2, 0)
		+ P1[0] * n0[2] * n1[2] * n2[1] * R(2, 0) + P2[1] * n0[1] * n1[2] * n2[1] * R(1, 1) - P2[1] * n0[2] * n1[1] * n2[1] * R(1, 1)
		- P1[1] * n0[1] * n1[2] * n2[2] * R(2, 1) + P1[1] * n0[2] * n1[2] * n2[1] * R(2, 1) + P2[2] * n0[1] * n1[2] * n2[1] * R(1, 2)
		- P2[2] * n0[2] * n1[1] * n2[1] * R(1, 2) - P1[2] * n0[1] * n1[2] * n2[2] * R(2, 2) + P1[2] * n0[2] * n1[2] * n2[1] * R(2, 2)
		+ P2[0] * n0[1] * n1[2] * n2[2] * R(2, 0) - P2[0] * n0[2] * n1[1] * n2[2] * R(2, 0) + P2[1] * n0[1] * n1[2] * n2[2] * R(2, 1)
		- P2[1] * n0[2] * n1[1] * n2[2] * R(2, 1) + P2[2] * n0[1] * n1[2] * n2[2] * R(2, 2) - P2[2] * n0[2] * n1[1] * n2[2] * R(2, 2)) /
		(n0[0] * n1[1] * n2[2] - n0[0] * n1[2] * n2[1] - n0[1] * n1[0] * n2[2] + n0[1] * n1[2] * n2[0] + n0[2] * n1[0] * n2[1]
			- n0[2] * n1[1] * n2[0]);

	double t2 = (P0[0] * n0[0] * n1[0] * n2[2] * R(0, 0) - P0[0] * n0[0] * n1[2] * n2[0] * R(0, 0) + P0[1] * n0[0] * n1[0] * n2[2] * R(0, 1)
		- P0[1] * n0[0] * n1[2] * n2[0] * R(0, 1) + P0[2] * n0[0] * n1[0] * n2[2] * R(0, 2) - P0[2] * n0[0] * n1[2] * n2[0] * R(0, 2)
		- P1[0] * n0[0] * n1[0] * n2[2] * R(0, 0) + P1[0] * n0[2] * n1[0] * n2[0] * R(0, 0) + P0[0] * n0[1] * n1[0] * n2[2] * R(1, 0)
		- P0[0] * n0[1] * n1[2] * n2[0] * R(1, 0) - P1[1] * n0[0] * n1[0] * n2[2] * R(0, 1) + P1[1] * n0[2] * n1[0] * n2[0] * R(0, 1)
		+ P0[1] * n0[1] * n1[0] * n2[2] * R(1, 1) - P0[1] * n0[1] * n1[2] * n2[0] * R(1, 1) - P1[2] * n0[0] * n1[0] * n2[2] * R(0, 2)
		+ P1[2] * n0[2] * n1[0] * n2[0] * R(0, 2) + P0[2] * n0[1] * n1[0] * n2[2] * R(1, 2) - P0[2] * n0[1] * n1[2] * n2[0] * R(1, 2)
		+ P2[0] * n0[0] * n1[2] * n2[0] * R(0, 0) - P2[0] * n0[2] * n1[0] * n2[0] * R(0, 0) - P1[0] * n0[0] * n1[1] * n2[2] * R(1, 0)
		+ P1[0] * n0[2] * n1[1] * n2[0] * R(1, 0) + P0[0] * n0[2] * n1[0] * n2[2] * R(2, 0) - P0[0] * n0[2] * n1[2] * n2[0] * R(2, 0)
		+ P2[1] * n0[0] * n1[2] * n2[0] * R(0, 1) - P2[1] * n0[2] * n1[0] * n2[0] * R(0, 1) - P1[1] * n0[0] * n1[1] * n2[2] * R(1, 1)
		+ P1[1] * n0[2] * n1[1] * n2[0] * R(1, 1) + P0[1] * n0[2] * n1[0] * n2[2] * R(2, 1) - P0[1] * n0[2] * n1[2] * n2[0] * R(2, 1)
		+ P2[2] * n0[0] * n1[2] * n2[0] * R(0, 2) - P2[2] * n0[2] * n1[0] * n2[0] * R(0, 2) - P1[2] * n0[0] * n1[1] * n2[2] * R(1, 2)
		+ P1[2] * n0[2] * n1[1] * n2[0] * R(1, 2) + P0[2] * n0[2] * n1[0] * n2[2] * R(2, 2) - P0[2] * n0[2] * n1[2] * n2[0] * R(2, 2)
		+ P2[0] * n0[0] * n1[2] * n2[1] * R(1, 0) - P2[0] * n0[2] * n1[0] * n2[1] * R(1, 0) - P1[0] * n0[0] * n1[2] * n2[2] * R(2, 0)
		+ P1[0] * n0[2] * n1[2] * n2[0] * R(2, 0) + P2[1] * n0[0] * n1[2] * n2[1] * R(1, 1) - P2[1] * n0[2] * n1[0] * n2[1] * R(1, 1)
		- P1[1] * n0[0] * n1[2] * n2[2] * R(2, 1) + P1[1] * n0[2] * n1[2] * n2[0] * R(2, 1) + P2[2] * n0[0] * n1[2] * n2[1] * R(1, 2)
		- P2[2] * n0[2] * n1[0] * n2[1] * R(1, 2) - P1[2] * n0[0] * n1[2] * n2[2] * R(2, 2) + P1[2] * n0[2] * n1[2] * n2[0] * R(2, 2)
		+ P2[0] * n0[0] * n1[2] * n2[2] * R(2, 0) - P2[0] * n0[2] * n1[0] * n2[2] * R(2, 0) + P2[1] * n0[0] * n1[2] * n2[2] * R(2, 1)
		- P2[1] * n0[2] * n1[0] * n2[2] * R(2, 1) + P2[2] * n0[0] * n1[2] * n2[2] * R(2, 2) - P2[2] * n0[2] * n1[0] * n2[2] * R(2, 2))
		/ (n0[0] * n1[1] * n2[2] - n0[0] * n1[2] * n2[1] - n0[1] * n1[0] * n2[2] + n0[1] * n1[2] * n2[0] + n0[2] * n1[0] * n2[1] - n0[2] * n1[1] * n2[0]);

	double t3 = -(P0[0] * n0[0] * n1[0] * n2[1] * R(0, 0) - P0[0] * n0[0] * n1[1] * n2[0] * R(0, 0) + P0[1] * n0[0] * n1[0] * n2[1] * R(0, 1)
		- P0[1] * n0[0] * n1[1] * n2[0] * R(0, 1) + P0[2] * n0[0] * n1[0] * n2[1] * R(0, 2) - P0[2] * n0[0] * n1[1] * n2[0] * R(0, 2)
		- P1[0] * n0[0] * n1[0] * n2[1] * R(0, 0) + P1[0] * n0[1] * n1[0] * n2[0] * R(0, 0) + P0[0] * n0[1] * n1[0] * n2[1] * R(1, 0)
		- P0[0] * n0[1] * n1[1] * n2[0] * R(1, 0) - P1[1] * n0[0] * n1[0] * n2[1] * R(0, 1) + P1[1] * n0[1] * n1[0] * n2[0] * R(0, 1)
		+ P0[1] * n0[1] * n1[0] * n2[1] * R(1, 1) - P0[1] * n0[1] * n1[1] * n2[0] * R(1, 1) - P1[2] * n0[0] * n1[0] * n2[1] * R(0, 2)
		+ P1[2] * n0[1] * n1[0] * n2[0] * R(0, 2) + P0[2] * n0[1] * n1[0] * n2[1] * R(1, 2) - P0[2] * n0[1] * n1[1] * n2[0] * R(1, 2)
		+ P2[0] * n0[0] * n1[1] * n2[0] * R(0, 0) - P2[0] * n0[1] * n1[0] * n2[0] * R(0, 0) - P1[0] * n0[0] * n1[1] * n2[1] * R(1, 0)
		+ P1[0] * n0[1] * n1[1] * n2[0] * R(1, 0) + P0[0] * n0[2] * n1[0] * n2[1] * R(2, 0) - P0[0] * n0[2] * n1[1] * n2[0] * R(2, 0)
		+ P2[1] * n0[0] * n1[1] * n2[0] * R(0, 1) - P2[1] * n0[1] * n1[0] * n2[0] * R(0, 1) - P1[1] * n0[0] * n1[1] * n2[1] * R(1, 1)
		+ P1[1] * n0[1] * n1[1] * n2[0] * R(1, 1) + P0[1] * n0[2] * n1[0] * n2[1] * R(2, 1) - P0[1] * n0[2] * n1[1] * n2[0] * R(2, 1)
		+ P2[2] * n0[0] * n1[1] * n2[0] * R(0, 2) - P2[2] * n0[1] * n1[0] * n2[0] * R(0, 2) - P1[2] * n0[0] * n1[1] * n2[1] * R(1, 2)
		+ P1[2] * n0[1] * n1[1] * n2[0] * R(1, 2) + P0[2] * n0[2] * n1[0] * n2[1] * R(2, 2) - P0[2] * n0[2] * n1[1] * n2[0] * R(2, 2)
		+ P2[0] * n0[0] * n1[1] * n2[1] * R(1, 0) - P2[0] * n0[1] * n1[0] * n2[1] * R(1, 0) - P1[0] * n0[0] * n1[2] * n2[1] * R(2, 0)
		+ P1[0] * n0[1] * n1[2] * n2[0] * R(2, 0) + P2[1] * n0[0] * n1[1] * n2[1] * R(1, 1) - P2[1] * n0[1] * n1[0] * n2[1] * R(1, 1)
		- P1[1] * n0[0] * n1[2] * n2[1] * R(2, 1) + P1[1] * n0[1] * n1[2] * n2[0] * R(2, 1) + P2[2] * n0[0] * n1[1] * n2[1] * R(1, 2)
		- P2[2] * n0[1] * n1[0] * n2[1] * R(1, 2) - P1[2] * n0[0] * n1[2] * n2[1] * R(2, 2) + P1[2] * n0[1] * n1[2] * n2[0] * R(2, 2)
		+ P2[0] * n0[0] * n1[1] * n2[2] * R(2, 0) - P2[0] * n0[1] * n1[0] * n2[2] * R(2, 0) + P2[1] * n0[0] * n1[1] * n2[2] * R(2, 1)
		- P2[1] * n0[1] * n1[0] * n2[2] * R(2, 1) + P2[2] * n0[0] * n1[1] * n2[2] * R(2, 2) - P2[2] * n0[1] * n1[0] * n2[2] * R(2, 2))
		/ (n0[0] * n1[1] * n2[2] - n0[0] * n1[2] * n2[1] - n0[1] * n1[0] * n2[2] + n0[1] * n1[2] * n2[0] + n0[2] * n1[0] * n2[1]
			- n0[2] * n1[1] * n2[0]);

	Matrix<double, 3, 1>B;

	B(0, 0) = -n0[0] * (P0[0] * R(0, 0) + P0[1] * R(0, 1) + P0[2] * R(0, 2)) - n0[1] * (P0[0] * R(1, 0) + P0[1] * R(1, 1) + P0[2] * R(1, 2)) - n0[2] * (P0[0] * R(2, 0) + P0[1] * R(2, 1) + P0[2] * R(2, 2));
	B(1, 0) = -n1[0] * (P1[0] * R(0, 0) + P1[1] * R(0, 1) + P1[2] * R(0, 2)) - n1[1] * (P1[0] * R(1, 0) + P1[1] * R(1, 1) + P1[2] * R(1, 2)) - n1[2] * (P1[0] * R(2, 0) + P1[1] * R(2, 1) + P1[2] * R(2, 2));
	B(2, 0) = -n2[0] * (P2[0] * R(0, 0) + P2[1] * R(0, 1) + P2[2] * R(0, 2)) - n2[1] * (P2[0] * R(1, 0) + P2[1] * R(1, 1) + P2[2] * R(1, 2)) - n2[2] * (P2[0] * R(2, 0) + P2[1] * R(2, 1) + P2[2] * R(2, 2));

	Matrix3d A;
	A << n0[0], n0[1], n0[2], n1[0], n1[1], n1[2], n2[0], n2[1], n2[2];

	Vector3d x = A.ldlt().solve(B);
}

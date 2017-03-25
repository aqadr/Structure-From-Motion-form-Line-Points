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
 * @file TriFocalTensor.h
 * @brief Trifocal tensor estimation from a set of corresponding points in three views
 * @author Ashraf Qadir, <ashraf.qadr@gmail.com>
 */

#ifndef TRIFOCALTENSOR_H
#define TRIFOCALTENSOR_H

//std libraries
#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <cmath>
#include <numeric>
#include <QString>
#include <QFile.h>
#include <QTextStream.h>

//Eigen
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SparseCore>
#include <Eigen/SVD>
#include <Eigen/StdVector>
#include <Eigen/Dense>


#include <boost/concept_check.hpp>
#include "PoseFromLine.h"
#include "MathHelper.h"

using namespace std;
using namespace Eigen;


namespace sfm 
{

template<typename T>
inline T Square( T x )
{
    return x * x;
}

class TrifocalTensor
{
public:
	TrifocalTensor()
	{}

	~TrifocalTensor()
	{
	}

	TrifocalTensor(const std::vector<Vector2d> &ptsFrame1,
		const std::vector<Vector2d> &ptsFrame2,
		const std::vector<Vector2d> &ptsFrame3);

	
	//!function for computing the epipole from trifocal tensor T. Refer to Hartley and Zisserman p373
	//! @param T   = 3x3x3 trifocal tensor is a 3x9 matrix
	//! @output e2 = 3x1 the epipole in image 2
	//! @output e3 = 3x1 the epipole in image 3
	static void epipolsFromTrifocalTensor(Eigen::Matrix<double, 3, 9> T, Vector3d &e1, Vector3d &e2);

	
	//! Function to build the relationship matrix which represent
	//! T_i^jk = a_i^j * b_4^k  -  a_4^j * b_i^k
	//! as t = E * aa, where aa = [a'(:) ; b'(:)], (note: for representation only)
	static void EFromEpipoles(Eigen::Matrix<double, 27, 18> &E,
		const Eigen::Vector3d &e2,
		const Eigen::Vector3d &e3);


	//! compute E matrix from the epipoles
	static void computeEFromEpipoles(Eigen::MatrixXd &E, const Eigen::Vector3d &e2, const Eigen::Vector3d &e3);


	//! trifocal tensor from three camera matrices
	//! P=[I|0], P1=[R1 | t1], P2=[R2 | t2]
	//! Ti = ai*b4^T-a4*bi^T
	static void Tfrom3Ps(Eigen::Matrix<double, 3, 9>&T, const Eigen::Matrix<double, 3, 4> &P1, const Eigen::Matrix<double, 3, 4> &P2,
		const Eigen::Matrix<double, 3, 4> &P3);


	

	//! compute trifocal tensor from atleast 7 points in three views
	bool computeTrifocalTensorPoints(const std::vector<Eigen::Vector2d> &vImPts1,
		  const std::vector<Eigen::Vector2d> &vImPts2,
		  const std::vector<Eigen::Vector2d> &vImPts3);

	//! compute trifocal tensor from lines. Atleast 13 lines are needed
	//bool computeTrifocalTensorFromLinesLinear(const std::vector<Eigen::Vector4d> &vImLines1,
	//	const std::vector<Eigen::Vector4d> &vImLines2,
	//	const std::vector<Eigen::Vector4d> &vImLines3);


		

	struct corrThree
	{
		Vector2d pt1;
		Vector2d pt2;
		Vector2d pt3;

		corrThree()
		{
		}

		corrThree(const Vector2d &_pt1, const Vector2d &_pt2, const Vector2d &_pt3)
			:pt1(_pt1), pt2(_pt2), pt3(_pt3)
		{
		}
	};
protected:
	std::vector<corrThree> ptsCorrespondences;

};

}


#endif

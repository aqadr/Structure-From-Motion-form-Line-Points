#pragma once
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
* @file MathHelper.h
* @brief helper class for performing normalization of points and lines
* @author Ashraf Qadir, <ashraf.qadr@gmail.com>
*/

//Eigen
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Dense>

using namespace Eigen;
namespace sfm
{
	template<typename PType>
	void lineNormalizedCoordinate(const Matrix<PType, 3, 3> &K, Matrix<PType, 4, 1> &EndPtsIm);

	bool normalize2DPoints(std::vector<Eigen::Vector2d> &vPts, Eigen::Matrix3d &T);

	/**Function computes the vector x that minimizes ||Ax|| subject to the
	* condition ||x||=1 and x = G*y, where G has rank r.
	* Implementation of Hartley and Zisserman A5.6 on p595 (2nd Ed)
	* @param: A - The constraint matrix, ||Ax|| to be minimized.
	* @param: G - The condition matrix, x = G*y.
	* @output: x - The vector that minimizes ||Ax|| subject to the condition ||x||=1 and x = G*v, where G has rank r.
	* @output: v - The vector that makes up x=G*v
	*/
	void constraintMinimization(MatrixXd A, MatrixXd G, VectorXd &X);

}

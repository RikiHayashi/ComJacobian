/*                                                                                                                                                                         
																																										   The BSD License (BSD)

																																										   Copyright (c) 2016 RDC lab and Chiba Institute of Technology.

																																										   Redistribution and use in source and binary forms, with or without
																																										   modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Jacobian.h"

// 全リンクがあワールド座標系の原点回りに作り出すモーメント
Matrix<double,3,1> calcMC(Link *ulink, int rootlink)
{
	Matrix<double,3,1> mc;

	if(rootlink == -1){
		mc = Matrix<double,3,1>::Zero();
	}else{
		mc = ulink[rootlink].m * (ulink[rootlink].p + ulink[rootlink].R * ulink[rootlink].c);
		mc = mc + calcMC(ulink, ulink[rootlink].sister) + calcMC(ulink, ulink[rootlink].child);
	}

	return mc;
}

// ロボットの質量を計算
double calcTotalMass(Link *ulink, int rootlink)
{
	double M = 0;

	if(rootlink == -1){
		M = 0.0f;
	}else{
		M = ulink[rootlink].m + calcTotalMass(ulink, ulink[rootlink].sister) + calcTotalMass(ulink, ulink[rootlink].child);
	}

	return M;
}

// ロボットの重心位置計算
Matrix<double,3,1> calcCoM(Link *ulink)
{
	double M = calcTotalMass(ulink, WAIST);
	Matrix<double,3,1> MC = calcMC(ulink, WAIST); 
	
	Matrix<double,3,1> com = MC / M;

	return com;
}

MatrixXd calcJacobian(Link *ulink, std::vector<int> idx)
{
	std::size_t jsize = idx.size();
	Matrix<double,3,1> target = ulink[idx.back()].p;
	MatrixXd J;

	J.resize(6,jsize);

	for(std::size_t i=0;i<jsize;i++)
	{
		int j = idx[i];
		Matrix<double,3,1> a = ulink[j].R * ulink[j].a;
		Matrix<double,3,1> b = a.cross(target - ulink[j].p);
		J(0,i) = b(0); J(1,i) = b(1); J(2,i) = b(2);
		J(3,i) = a(0); J(4,i) = a(1); J(5,i) = a(2);
	}

	return J;
}

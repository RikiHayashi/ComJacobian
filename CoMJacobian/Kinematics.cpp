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
#include "Kinematics.h"
#include <Eigen/SVD>

using namespace MotionControl;

Matrix<double,3,3> Kinematics::Rodrigues(Matrix<double,3,1> a, double q)
{
    return AngleAxisd(q,a).toRotationMatrix();
}

Matrix<double,3,1> Kinematics::rot2omega(Matrix<double,3,3> R)
{
	double alpha = (R(0,0)+R(1,1)+R(2,2)-1)/2;
	double th;
	Matrix<double,3,1> vector_R(Matrix<double,3,1>::Zero());

	if(fabs(alpha-1) < eps)
		return Matrix<double,3,1>::Zero();

	th = acos(alpha);
	vector_R << R(2,1)-R(1,2), R(0,2)-R(2,0), R(1,0)-R(0,1);
	return 0.5*th/sin(th)*vector_R;
}

vector<int> Kinematics::FindRoute(int to)
{
	vector<int> idx;
	int link_num = to;

	while(link_num != 0)
	{
		idx.push_back(link_num);
		link_num = ulink[link_num].parent;
	}
	reverse(idx.begin(), idx.end());
	return idx;
}

vector<int> Kinematics::CustomLinkPath(int to, int from)
{
	vector<int> idx;
	vector<int> idy;
	int link_num = from;

	while(link_num != 0)
	{
		idx.push_back(link_num);
		link_num = ulink[link_num].parent;
	}

	link_num = to;
	while(link_num != 0)
	{
		idy.push_back(link_num);
		link_num = ulink[link_num].parent;
	}
	reverse(idy.begin(), idy.end());

	std::copy(idy.begin(), idy.end(), std::back_inserter(idx));
	
	return idx;
}

Matrix<double,6,1> Kinematics::calcVWerr(Link Cref, Link Cnow)
{
	Matrix<double,3,1> perr = Cref.p - Cnow.p;
	Matrix<double,3,3> Rerr = Cref.R - Cnow.R;
	Matrix<double,3,1> werr = Cnow.R * rot2omega(Rerr);
	Matrix<double,6,1> err;

	err << perr,werr;
	return err;
}

void Kinematics::calcForwardKinematics(int rootlink)
{
	if(rootlink == -1)
		return ;
	if(rootlink != 0)
	{
		int parent = ulink[rootlink].parent;
		ulink[rootlink].p = ulink[parent].R * ulink[rootlink].b + ulink[parent].p;
		ulink[rootlink].R = ulink[parent].R * Rodrigues(ulink[rootlink].a, ulink[rootlink].q);
		ulink[rootlink].p_com = ulink[rootlink].p +  ulink[rootlink].R * ulink[rootlink].compos;
	}
	if(rootlink == 0)
		ulink[rootlink].p_com = ulink[rootlink].compos;

	calcForwardKinematics(ulink[rootlink].sister);
	calcForwardKinematics(ulink[rootlink].child);
}

MatrixXd Kinematics::calcJacobian(vector<int> idx)
{
	size_t jsize = idx.size();
	Matrix<double,3,1> target = ulink[idx.back()].p;
	Matrix<double,6,12> J = MatrixXd::Zero(6,12);

	for(size_t i=0;i<jsize;i++)
	{
		int j = idx[i];
		Matrix<double,3,1> a = ulink[j].R * ulink[j].a;
		Matrix<double,3,1> b = a.cross(target - ulink[j].p);
		J(0,i) = b(0); J(1,i) = b(1); J(2,i) = b(2);
		J(3,i) = a(0); J(4,i) = a(1); J(5,i) = a(2);
	}

	return J;
}

template <typename t_matrix>
t_matrix Kinematics::PseudoInverse(const t_matrix& m, const double &tolerance)
{
	typedef JacobiSVD<t_matrix> TSVD;
	unsigned int svd_opt(ComputeThinU | ComputeThinV);
	if(m.RowsAtCompileTime!=Dynamic || m.ColsAtCompileTime!=Dynamic)
		svd_opt= ComputeFullU | ComputeFullV;
	TSVD svd(m, svd_opt);
	const typename TSVD::SingularValuesType &sigma(svd.singularValues());
	typename TSVD::SingularValuesType sigma_inv(sigma.size());
	for(long i=0; i<sigma.size(); ++i)
	{
		if(sigma(i) > tolerance)
			sigma_inv(i)= 1.0/sigma(i);
		else
			sigma_inv(i)= 0.0;
	}
	return svd.matrixV()*sigma_inv.asDiagonal()*svd.matrixU().transpose();
}

bool Kinematics::calcInverseKinematics(int to, Link target)
{
	MatrixXd J, dq;
	Matrix<double,6,1> err;

	ColPivHouseholderQR<MatrixXd> QR; //QR分解?
	const double dampingConstantSqr = 1.0e-12;
	const double lambda = 0.5;
	const int iteration = 100;
	
	calcForwardKinematics(WAIST);
	
	vector<int> idx = FindRoute(to);
	const int jsize = idx.size();
	
	J.resize(6,jsize); dq.resize(jsize,1);

	for(int n=0;n<iteration;n++){
		J = calcJacobian(idx);
		err = calcVWerr(target, ulink[to]);
		if(err.norm() < eps) return true;

		MatrixXd JJ = J*J.transpose()+dampingConstantSqr*MatrixXd::Identity(J.rows(),J.rows());
		dq = J.transpose() * QR.compute(JJ).solve(err) * lambda;
		
		for(size_t nn=0;nn<jsize;nn++){
			int j = idx[nn];
			ulink[j].q += dq(nn);
		}
		calcForwardKinematics(WAIST);
	}
	return false;
}

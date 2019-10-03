#include "Kinematics.h"
#include "ComJacobian.h"
#include "func.h"

using namespace std;
using namespace MotionControl;

	template <typename t_matrix>
t_matrix PseudoInverse(const t_matrix& m, const double &tolerance=1.e-6)
{
	using namespace Eigen;
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

int main(int argc, char* argv[])
{
	//Initialize
	Link ulink[JOINT_NUM];
	Kinematics kine(ulink);	
	SetJointInfo(ulink);
	CoMJacobian CJ(ulink);

	//Set Angle to Joint
	for(int i=0;i<JOINT_NUM;i++)
		ulink[i].q = 0.0;

	ulink[RLEG_JOINT0].q = deg2rad(0);
	ulink[RLEG_JOINT1].q = deg2rad(0);
	ulink[RLEG_JOINT2].q = deg2rad(15);
	ulink[RLEG_JOINT3].q = deg2rad(-30);
	ulink[RLEG_JOINT4].q = deg2rad(15);
	ulink[RLEG_JOINT5].q = deg2rad(0);
	ulink[LLEG_JOINT2].q = deg2rad(-15);
	ulink[LLEG_JOINT3].q = deg2rad(30);
	ulink[LLEG_JOINT4].q = deg2rad(-15);

	//Calculation Forward Kinematics
	kine.calcForwardKinematics(WAIST);
	Link RLEG_LINK = ulink[RLEG_JOINT1];
	Link LLEG_LINK = ulink[LLEG_JOINT5];

	MatrixXd comJ_r(3,6);
	MatrixXd comJ_l;
	Matrix<double,6,6> comJ;

	comJ_r = CJ.resultCJ(RLEG_JOINT5);
	comJ_l = CJ.resultCJ(LLEG_JOINT5);

	cout << "\nPseudo:" << PseudoInverse(comJ_r) << endl;

	comJ << comJ_r(0), comJ_r(3), comJ_r(6), comJ_r(9), comJ_r(12),comJ_r(15),
			comJ_r(1), comJ_r(4), comJ_r(7), comJ_r(10), comJ_r(13),comJ_r(16),
			comJ_r(2), comJ_r(5), comJ_r(8), comJ_r(11), comJ_r(14),comJ_r(17),
			comJ_r(2), comJ_r(5), comJ_r(8), comJ_r(11), comJ_r(14),comJ_r(15),
			comJ_r(1), comJ_r(4), comJ_r(7), comJ_r(10), comJ_r(13),comJ_r(16),
			comJ_r(0), comJ_r(3), comJ_r(6), comJ_r(9), comJ_r(12),comJ_r(17);
	comJ = comJ_r.transpose() * comJ_r;

	cout << "comJ\n" << comJ << endl;
	comJ = comJ.inverse();

	Matrix<double,6,1> dp_r = Matrix<double,6,1>::Zero();
	Matrix<double,6,1> dp_l = Matrix<double,6,1>::Zero();

	Matrix<double,6,3> ref_comJ_r;
	Matrix<double,6,3> ref_comJ_l;
	
	ref_comJ_r << PseudoInverse(comJ_r);
	
	ref_comJ_l <<	comJ(18), comJ(24), comJ(30),
					comJ(19), comJ(25), comJ(31),
					comJ(20), comJ(26), comJ(32),
					comJ(21), comJ(27), comJ(33),
					comJ(22), comJ(28), comJ(34),
					comJ(23), comJ(29), comJ(35);
	ref_comJ_l << PseudoInverse(comJ_l);

	cout << "ref_comJ_r\n" << ref_comJ_r << endl;
	cout << "ref_comJ_l\n" << ref_comJ_l << endl;

	Matrix<double,3,1> t;
	Matrix<double,1,3> t_i;

	t << 0.000143,0.001037,0;
	t_i << 0.11,0.17,0;

	cout << "ref_r" << ref_comJ_r * t << endl;
	cout << "ref_l" << ref_comJ_l * t << endl;
	
	//cout << "ref_r" << t_i * comJ_r << endl;

	return 0;
}

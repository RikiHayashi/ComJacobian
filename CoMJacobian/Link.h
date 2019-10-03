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
/*
 * @file		Link.h
 * @brief		Link Information
 * @author	Ryu Yamamoto
 * @date		2016/02/26
 */

#ifndef _LINK_
#define _LINK_

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "LinkParameter.h"

using namespace std;
using namespace Eigen;

namespace MotionControl
{
	struct Link
	{
		string joint_name;					
		int sister;						
		int child;						
		int parent;						
		Matrix<double,3,1> p;	
		Matrix<double,3,1> p_com;	
		Matrix<double,3,3> R;	
		Matrix<double,3,1> v;	
		Matrix<double,3,1> w;	
		Matrix<double,3,1> a;	
		Matrix<double,3,1> b;	
		double q;							
		double dq;						
		double ddq;						
		double vartex;
		double face;
		double m;

		Matrix<double,3,1> c;
		Matrix<double,3,3> I;

		double mass;
		Matrix<double,3,1> compos;

		Link() : p(Matrix<double,3,1>::Zero()), R(Matrix<double,3,3>::Identity()),
						v(Matrix<double,3,1>::Zero()), w(Matrix<double,3,1>::Zero()),
						a(Matrix<double,3,1>::Zero()), b(Matrix<double,3,1>::Zero()),
						compos(Matrix<double,3,1>::Zero()),
						q(0.0), dq(0.0), ddq(0.0),
						c(Matrix<double,3,1>::Zero()), I(Matrix<double,3,3>::Identity()) 
		{
		}
	};

	extern "C" void SetJointInfo(struct Link *link);
};

#endif

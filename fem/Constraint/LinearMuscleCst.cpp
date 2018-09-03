#include "LinearMuscleCst.h"
using namespace FEM;

LinearMuscleCst::
LinearMuscleCst(const std::string& name,double k,int i0,int i1,int i2,int i3,double volume,const Eigen::Matrix3d& invDm,const Eigen::Vector3d& fiber_direction)
	:Cst(name,k),mi0(i0),mi1(i1),mi2(i2),mi3(i3),mVolume(volume),mInvDm(invDm),mFiberDirection(fiber_direction),mActivationLevel(0.0)
{

}
std::shared_ptr<Cst>
LinearMuscleCst::
Clone()
{
	return Create(mName,mStiffness,mi0,mi1,mi2,mi3,mVolume,mInvDm,mFiberDirection);
}
std::shared_ptr<LinearMuscleCst>
LinearMuscleCst::
Create(const std::string& name,double k,int i0,int i1,int i2,int i3,double volume,const Eigen::Matrix3d& invDm,const Eigen::Vector3d& fiber_direction)
{
	auto c = new LinearMuscleCst(name,k,i0,i1,i2,i3,volume,invDm,fiber_direction);
	return std::shared_ptr<LinearMuscleCst>(c);
}
void
LinearMuscleCst::
EvaluatePotentialEnergy(const Eigen::VectorXd& x)
{
	ComputeF(x);
	Computep0();

	mE = 0.5*mVolume*mStiffness*((mF*mFiberDirection - mp0).squaredNorm());
}
void
LinearMuscleCst::
EvaluateGradient(const Eigen::VectorXd& x)
{
	ComputeF(x);
	Computep0();

	Eigen::Matrix3d P;
	ComputeP(P);

	P = mVolume*P*(mInvDm.transpose());
	mg.block<3,1>(0*3,0) = -(P.block<3,1>(0,0) + P.block<3,1>(0,1) + P.block<3,1>(0,2));
	mg.block<3,1>(1*3,0) = P.block<3,1>(0,0);
	mg.block<3,1>(2*3,0) = P.block<3,1>(0,1);
	mg.block<3,1>(3*3,0) = P.block<3,1>(0,2);
}
void
LinearMuscleCst::
EvaluateHessian(const Eigen::VectorXd& x)
{
	ComputeF(x);
	Computep0();

	Tensor3333 dPdF;

	ComputedPdF(dPdF);

	dPdF = mVolume*dPdF*mInvDm.transpose();

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			//H_ij
			Eigen::Vector3d df0,df1,df2;
			df0 = dPdF(0,j)*mInvDm.transpose().block<3,1>(0,i);
			df1 = dPdF(1,j)*mInvDm.transpose().block<3,1>(0,i);
			df2 = dPdF(2,j)*mInvDm.transpose().block<3,1>(0,i);

			mH.block<3,1>(i*3,j*3+0) = df0;
			mH.block<3,1>(i*3,j*3+1) = df1;
			mH.block<3,1>(i*3,j*3+2) = df2;
		}

		mH.block<3,3>(i*3,9) = -mH.block<3,3>(i*3,0)-mH.block<3,3>(i*3,3)-mH.block<3,3>(i*3,6);
	}
	for(int i=0;i<4;i++)
	{
		mH.block<3,3>(9,i*3) = -mH.block<3,3>(0,i*3)-mH.block<3,3>(3,i*3)-mH.block<3,3>(6,i*3);
	}
}
void
LinearMuscleCst::
GetPotentialEnergy(double& e)
{
	e += mE;
}
void
LinearMuscleCst::
GetGradient(Eigen::VectorXd& g)
{
	g.block<3,1>(mi0*3,0) += mg.block<3,1>(0*3,0);
	g.block<3,1>(mi1*3,0) += mg.block<3,1>(1*3,0);
	g.block<3,1>(mi2*3,0) += mg.block<3,1>(2*3,0);
	g.block<3,1>(mi3*3,0) += mg.block<3,1>(3*3,0);
}
void
LinearMuscleCst::
GetHessian(std::vector<Eigen::Triplet<double>>& h_triplets)
{
	int idx[4] = {mi0,mi1,mi2,mi3};

	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			//H.block [i,j] --- 3x3 matrix
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*idx[j]+0, mH(3*i+0, 3*j+0)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*idx[j]+1, mH(3*i+0, 3*j+1)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*idx[j]+2, mH(3*i+0, 3*j+2)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*idx[j]+0, mH(3*i+1, 3*j+0)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*idx[j]+1, mH(3*i+1, 3*j+1)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*idx[j]+2, mH(3*i+1, 3*j+2)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*idx[j]+0, mH(3*i+2, 3*j+0)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*idx[j]+1, mH(3*i+2, 3*j+1)));
			h_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*idx[j]+2, mH(3*i+2, 3*j+2)));
		}
	}
}
void
LinearMuscleCst::
Evaluatedgda(const Eigen::VectorXd& x)
{
	ComputeF(x);
	Computep0();

	Eigen::Matrix3d dP_da;
	Eigen::Vector3d dp0;
	Computedp0(dp0);
	dP_da = - mVolume*mStiffness*dp0*(mFiberDirection.transpose())*mInvDm.transpose();;

	mdgda.block<3,1>(0*3,0) = -(dP_da.block<3,1>(0,0) + dP_da.block<3,1>(0,1) + dP_da.block<3,1>(0,2));
	mdgda.block<3,1>(1*3,0) = dP_da.block<3,1>(0,0);
	mdgda.block<3,1>(2*3,0) = dP_da.block<3,1>(0,1);
	mdgda.block<3,1>(3*3,0) = dP_da.block<3,1>(0,2);
}	
void
LinearMuscleCst::
Getdgda(Eigen::VectorXd& dgda)
{
	dgda.block<3,1>(mi0*3,0) += mdgda.block<3,1>(0*3,0);
	dgda.block<3,1>(mi1*3,0) += mdgda.block<3,1>(1*3,0);
	dgda.block<3,1>(mi2*3,0) += mdgda.block<3,1>(2*3,0);
	dgda.block<3,1>(mi3*3,0) += mdgda.block<3,1>(3*3,0);
}


void
LinearMuscleCst::
EvaluateDVector(const Eigen::VectorXd& x)
{
	ComputeF(x);
	Computep0();
}
void
LinearMuscleCst::
GetDVector(int& index,Eigen::VectorXd& d)
{
	d.block<3,1>(3*index,0) = mp0;

	index++;
}
void
LinearMuscleCst::
EvaluateJMatrix(int& index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	Eigen::MatrixXd Ai(3,12);

	Eigen::Vector3d v = mInvDm*mFiberDirection;

	double a,b,c;
	a = v[0];
	b = v[1];
	c = v[2];

	Ai<<
		-(a+b+c),0,0,a,0,0,b,0,0,c,0,0,
		0,-(a+b+c),0,0,a,0,0,b,0,0,c,0,
		0,0,-(a+b+c),0,0,a,0,0,b,0,0,c;

	Eigen::MatrixXd MuAiT = mVolume*mStiffness*Ai.transpose();

	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+0, 3*index+0, MuAiT(3*0+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+0, 3*index+1, MuAiT(3*0+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+0, 3*index+2, MuAiT(3*0+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+1, 3*index+0, MuAiT(3*0+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+1, 3*index+1, MuAiT(3*0+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+1, 3*index+2, MuAiT(3*0+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+2, 3*index+0, MuAiT(3*0+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+2, 3*index+1, MuAiT(3*0+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+2, 3*index+2, MuAiT(3*0+2,3*0+2)));

	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+0, 3*index+0, MuAiT(3*1+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+0, 3*index+1, MuAiT(3*1+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+0, 3*index+2, MuAiT(3*1+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+1, 3*index+0, MuAiT(3*1+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+1, 3*index+1, MuAiT(3*1+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+1, 3*index+2, MuAiT(3*1+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+2, 3*index+0, MuAiT(3*1+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+2, 3*index+1, MuAiT(3*1+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi1+2, 3*index+2, MuAiT(3*1+2,3*0+2)));

	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+0, 3*index+0, MuAiT(3*2+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+0, 3*index+1, MuAiT(3*2+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+0, 3*index+2, MuAiT(3*2+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+1, 3*index+0, MuAiT(3*2+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+1, 3*index+1, MuAiT(3*2+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+1, 3*index+2, MuAiT(3*2+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+2, 3*index+0, MuAiT(3*2+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+2, 3*index+1, MuAiT(3*2+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi2+2, 3*index+2, MuAiT(3*2+2,3*0+2)));

	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+0, 3*index+0, MuAiT(3*3+0,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+0, 3*index+1, MuAiT(3*3+0,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+0, 3*index+2, MuAiT(3*3+0,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+1, 3*index+0, MuAiT(3*3+1,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+1, 3*index+1, MuAiT(3*3+1,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+1, 3*index+2, MuAiT(3*3+1,3*0+2)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+2, 3*index+0, MuAiT(3*3+2,3*0+0)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+2, 3*index+1, MuAiT(3*3+2,3*0+1)));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi3+2, 3*index+2, MuAiT(3*3+2,3*0+2)));

	index++;
}
void
LinearMuscleCst::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	Eigen::MatrixXd Ai(3,12);

	Eigen::Vector3d v = mInvDm*mFiberDirection;

	double a,b,c;
	a = v[0];
	b = v[1];
	c = v[2];

	Ai<<
		-(a+b+c),0,0,a,0,0,b,0,0,c,0,0,
		0,-(a+b+c),0,0,a,0,0,b,0,0,c,0,
		0,0,-(a+b+c),0,0,a,0,0,b,0,0,c;

	auto MuAiTAi = mVolume*mStiffness*((Ai.transpose())*Ai);

	int idx[4] = {mi0,mi1,mi2,mi3};

	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			//H.block [i,j] --- 3x3 matrix
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*idx[j]+0, MuAiTAi(3*i+0, 3*j+0)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*idx[j]+1, MuAiTAi(3*i+0, 3*j+1)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*idx[j]+2, MuAiTAi(3*i+0, 3*j+2)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*idx[j]+0, MuAiTAi(3*i+1, 3*j+0)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*idx[j]+1, MuAiTAi(3*i+1, 3*j+1)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*idx[j]+2, MuAiTAi(3*i+1, 3*j+2)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*idx[j]+0, MuAiTAi(3*i+2, 3*j+0)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*idx[j]+1, MuAiTAi(3*i+2, 3*j+1)));
			L_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*idx[j]+2, MuAiTAi(3*i+2, 3*j+2)));
		}
	}
}
int
LinearMuscleCst::
GetNumHessianTriplets()
{
	return 144;
}
void
LinearMuscleCst::
AddOffset(int offset)
{
	mi0 +=offset;
	mi1 +=offset;
	mi2 +=offset;
	mi3 +=offset;
}


void
LinearMuscleCst::
ComputeF(const Eigen::VectorXd& x)
{
	mX.block<3,1>(0*3,0) = x.block<3,1>(mi0*3,0);
	mX.block<3,1>(1*3,0) = x.block<3,1>(mi1*3,0);
	mX.block<3,1>(2*3,0) = x.block<3,1>(mi2*3,0);
	mX.block<3,1>(3*3,0) = x.block<3,1>(mi3*3,0);

	Eigen::Matrix3d Ds;
	Ds.block<3,1>(0,0) = mX.block<3,1>(1*3,0) - mX.block<3,1>(0*3,0);
	Ds.block<3,1>(0,1) = mX.block<3,1>(2*3,0) - mX.block<3,1>(0*3,0);
	Ds.block<3,1>(0,2) = mX.block<3,1>(3*3,0) - mX.block<3,1>(0*3,0);

	mF = Ds*mInvDm;
}
void
LinearMuscleCst::
ComputeP(Eigen::Matrix3d& P)
{
	P = mStiffness *(mF*mFiberDirection*mFiberDirection.transpose() - mp0*mFiberDirection.transpose());
}
void
LinearMuscleCst::
ComputedPdF(Tensor3333& dPdF)
{
	Tensor3333 dFdF;
	dFdF.SetIdentity();

	for(int i =0;i<3;i++)
		for(int j=0;j<3;j++)
			dPdF(i,j) = mStiffness*dFdF(i,j)*mFiberDirection*mFiberDirection.transpose();
}
void
LinearMuscleCst::
Computep0()
{
	mp0 = (1.0-mActivationLevel)*mF*mFiberDirection;
}
void
LinearMuscleCst::
Computedp0(Eigen::Vector3d& dp0)
{
	dp0 = -mF*mFiberDirection;
}
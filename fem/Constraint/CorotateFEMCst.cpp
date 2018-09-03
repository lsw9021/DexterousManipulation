#include "CorotateFEMCst.h"
#include <Eigen/SVD>
#include <Eigen/LU>

using namespace FEM;
CorotateFEMCst::
CorotateFEMCst(const std::string& name,double k,double poisson_ratio,int i0,int i1,int i2,int i3,double volume,const Eigen::Matrix3d& invDm)
	:Cst(name,k),mPoissonRatio(poisson_ratio),mi0(i0),mi1(i1),mi2(i2),mi3(i3),mVolume(volume),mInvDm(invDm),
	mMu(k/((1.0+poisson_ratio))),mLambda(k*poisson_ratio/((1.0+poisson_ratio)*(1-2.0*poisson_ratio)))
{
	mE = 0.0;
	mg.setZero();
	mH.setZero();
	md.setZero();
	mX.setZero();
	mF.setZero();
	mR.setZero();
	mU.setZero();
	mV.setZero();
	mD.setZero();
}
std::shared_ptr<Cst>
CorotateFEMCst::
Clone()
{
	return Create(mName,mStiffness,mPoissonRatio,mi0,mi1,mi2,mi3,mVolume,mInvDm);
}
std::shared_ptr<CorotateFEMCst>
CorotateFEMCst::
Create(const std::string& name,double k,double poisson_ratio,int i0,int i1,int i2,int i3,double volume,const Eigen::Matrix3d& invDm)
{
	auto c = new CorotateFEMCst(name,k,poisson_ratio,i0,i1,i2,i3,volume,invDm);
	return std::shared_ptr<CorotateFEMCst>(c);
}
void
CorotateFEMCst::
EvaluatePotentialEnergy(const Eigen::VectorXd& x)
{
	ComputeF(x);

	double vol_preserve = (mD-Eigen::Matrix3d::Identity()).trace();
	mE = mVolume*(0.5*mMu*((mF - mR).norm())+0.5*mLambda*vol_preserve*vol_preserve);
}
void
CorotateFEMCst::
EvaluateGradient(const Eigen::VectorXd& x)
{
	ComputeF(x);

	Eigen::Matrix3d P;
	ComputeP(P);

	P = mVolume*P*(mInvDm.transpose());
	mg.block<3,1>(0*3,0) = -(P.block<3,1>(0,0) + P.block<3,1>(0,1) + P.block<3,1>(0,2));
	mg.block<3,1>(1*3,0) = P.block<3,1>(0,0);
	mg.block<3,1>(2*3,0) = P.block<3,1>(0,1);
	mg.block<3,1>(3*3,0) = P.block<3,1>(0,2);
}
void
CorotateFEMCst::
EvaluateHessian(const Eigen::VectorXd& x)
{
	ComputeF(x);

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
CorotateFEMCst::
GetPotentialEnergy(double& e)
{
	e += mE;
}
void
CorotateFEMCst::
GetGradient(Eigen::VectorXd& g)
{
	g.block<3,1>(mi0*3,0) += mg.block<3,1>(0*3,0);
	g.block<3,1>(mi1*3,0) += mg.block<3,1>(1*3,0);
	g.block<3,1>(mi2*3,0) += mg.block<3,1>(2*3,0);
	g.block<3,1>(mi3*3,0) += mg.block<3,1>(3*3,0);
}
void
CorotateFEMCst::
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
CorotateFEMCst::
EvaluateDVector(const Eigen::VectorXd& x)
{
	ComputeF(x);

	md = mR;
	if(mD(2,2)<0)
		md.block<3,1>(0,2) = -mR.block<3,1>(0,2);

	Eigen::Vector3d S = mD.diagonal();
	Eigen::Vector3d D;
	D.setZero();
	for(int i=0;i<10;i++)
	{
		double CD = (S[0]+D[0])*(S[1]+D[1])*(S[2]+D[2])-1;
		Eigen::Vector3d gradCD( (S[1]+D[1])*(S[2]+D[2]),
								(S[0]+D[0])*(S[2]+D[2]),
								(S[0]+D[0])*(S[1]+D[1]));

		D = (gradCD.dot(D) -CD)/(gradCD.squaredNorm())*gradCD;

	}

	md_volume = mU*((S+D).asDiagonal())*mV.transpose();
}
void
CorotateFEMCst::
GetDVector(int& index,Eigen::VectorXd& d)
{
	d.block<3,1>(3*(index+0),0) = md.block<3,1>(0,0);
	d.block<3,1>(3*(index+1),0) = md.block<3,1>(0,1);
	d.block<3,1>(3*(index+2),0) = md.block<3,1>(0,2);
	index+=3;

	d.block<3,1>(3*(index+0),0) = md_volume.block<3,1>(0,0);
	d.block<3,1>(3*(index+1),0) = md_volume.block<3,1>(0,1);
	d.block<3,1>(3*(index+2),0) = md_volume.block<3,1>(0,2);
	index+=3;
}
void
CorotateFEMCst::
EvaluateJMatrix(int& index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	Eigen::MatrixXd Ai(3*3,3*4);
	double d11 = mInvDm(0,0);
	double d12 = mInvDm(0,1);
	double d13 = mInvDm(0,2);
	double d21 = mInvDm(1,0);
	double d22 = mInvDm(1,1);
	double d23 = mInvDm(1,2);
	double d31 = mInvDm(2,0);
	double d32 = mInvDm(2,1);
	double d33 = mInvDm(2,2);

	Ai<<
		-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,0,
		0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,
		0,0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,
		-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,0,
		0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,
		0,0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,
		-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,0,
		0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,
		0,0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33;

	Eigen::MatrixXd MuAiT = mMu*mVolume*Ai.transpose();
	int idx[4] = {mi0,mi1,mi2,mi3};
	//MuAiT --- 12x9 matrix
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<3;j++)
		{
			//MuAiT.block [i,j] -- 3x3 matrix
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*(index+j)+0, MuAiT(3*i+0, 3*j+0)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*(index+j)+1, MuAiT(3*i+0, 3*j+1)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*(index+j)+2, MuAiT(3*i+0, 3*j+2)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*(index+j)+0, MuAiT(3*i+1, 3*j+0)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*(index+j)+1, MuAiT(3*i+1, 3*j+1)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*(index+j)+2, MuAiT(3*i+1, 3*j+2)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*(index+j)+0, MuAiT(3*i+2, 3*j+0)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*(index+j)+1, MuAiT(3*i+2, 3*j+1)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*(index+j)+2, MuAiT(3*i+2, 3*j+2)));
		}
	}
	index+=3;

	MuAiT = (MuAiT*mPoissonRatio).eval();
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<3;j++)
		{
			//MuAiT.block [i,j] -- 3x3 matrix
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*(index+j)+0, MuAiT(3*i+0, 3*j+0)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*(index+j)+1, MuAiT(3*i+0, 3*j+1)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+0, 3*(index+j)+2, MuAiT(3*i+0, 3*j+2)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*(index+j)+0, MuAiT(3*i+1, 3*j+0)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*(index+j)+1, MuAiT(3*i+1, 3*j+1)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+1, 3*(index+j)+2, MuAiT(3*i+1, 3*j+2)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*(index+j)+0, MuAiT(3*i+2, 3*j+0)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*(index+j)+1, MuAiT(3*i+2, 3*j+1)));
			J_triplets.push_back(Eigen::Triplet<double>(3*idx[i]+2, 3*(index+j)+2, MuAiT(3*i+2, 3*j+2)));
		}
	}
	index+=3;
}
void
CorotateFEMCst::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	Eigen::MatrixXd Ai(3*3,3*4);
	double d11 = mInvDm(0,0);
	double d12 = mInvDm(0,1);
	double d13 = mInvDm(0,2);
	double d21 = mInvDm(1,0);
	double d22 = mInvDm(1,1);
	double d23 = mInvDm(1,2);
	double d31 = mInvDm(2,0);
	double d32 = mInvDm(2,1);
	double d33 = mInvDm(2,2);

	Ai<<
		-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,0,
		0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,0,
		0,0,-d11-d21-d31,0,0,d11,0,0,d21,0,0,d31,
		-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,0,
		0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,0,
		0,0,-d12-d22-d32,0,0,d12,0,0,d22,0,0,d32,
		-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,0,
		0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33,0,
		0,0,-d13-d23-d33,0,0,d13,0,0,d23,0,0,d33;

	Eigen::MatrixXd MuAiTAi = mMu*mVolume*((Ai.transpose())*Ai);
	int idx[4] = {mi0,mi1,mi2,mi3};
	//MuAiT --- 12x12 matrix
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			//MuAiTAi.block [i,j] -- 3x3 matrix
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

	MuAiTAi = (MuAiTAi*mPoissonRatio).eval();
	for(int i =0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			//MuAiTAi.block [i,j] -- 3x3 matrix
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
CorotateFEMCst::
GetNumHessianTriplets()
{
	return 288;
}

void
CorotateFEMCst::
AddOffset(int offset)
{
	mi0 +=offset;
	mi1 +=offset;
	mi2 +=offset;
	mi3 +=offset;
}
void
CorotateFEMCst::
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

	Eigen::JacobiSVD<Eigen::Matrix3d> svd(mF, Eigen::ComputeFullU | Eigen::ComputeFullV);

	mD(0,0) = svd.singularValues()[0];
	mD(1,1) = svd.singularValues()[1];
	mD(2,2) = svd.singularValues()[2];

	mU = svd.matrixU();
	mV = svd.matrixV();
	mR = mU*mV.transpose();
}

void
CorotateFEMCst::
ComputeP(Eigen::Matrix3d& P)
{
	P = mMu*(mF - mR) +
		mLambda*((mR.transpose()*mF-Eigen::Matrix3d::Identity()).trace())*mR;
}
void
CorotateFEMCst::
ComputedPdF(Tensor3333& dPdF)
{
	Tensor3333 dFdF,dRdF;
	dFdF.SetIdentity();

	//Compute dRdF
	for(int i =0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			Eigen::Matrix3d M = mU.transpose()*dFdF(i,j)*mV;
			if(fabs(mD(0,0)-mD(1,1))<1E-4 && fabs(mD(0,0)-mD(2,2))<1E-4)
			{
				Eigen::Matrix3d off_diag_M;
				off_diag_M.setZero();
				for(int a=0;a<3;a++)
					for(int b=0;b<3;b++)		
					{
						if(a==b)
							continue;
						else
							off_diag_M(a,b) = M(a,b) / mD(0,0);
					}

				dRdF(i,j) = mU*off_diag_M*mV.transpose();
			}
			else
			{
				Eigen::Vector2d unknown_side, known_side;
				Eigen::Matrix2d known_matrix;
				Eigen::Matrix3d U_tilde, V_tilde;
				U_tilde.setZero(); V_tilde.setZero();
				Eigen::Matrix2d reg;
				reg.setZero();
				reg(0, 0) = reg(1, 1) = 1E-4;
				for (unsigned int row = 0; row < 3; row++)
				{
					for (unsigned int col = 0; col < row; col++)
					{
						known_side = Eigen::Vector2d(M(col, row), M(row, col));
						known_matrix.block<2, 1>(0, 0) = Eigen::Vector2d(-mD(row,row), mD(col,col));
						known_matrix.block<2, 1>(0, 1) = Eigen::Vector2d(-mD(col,col), mD(row,row));

						if (fabs(mD(row,row) - mD(col,col) < 1E-4))
							known_matrix += reg;
						else
							assert(fabs(known_matrix.determinant()) > 1E-6);

						unknown_side = known_matrix.inverse() * known_side;
						U_tilde(row, col) = unknown_side[0];
						U_tilde(col, row) = -U_tilde(row, col);
						V_tilde(row, col) = unknown_side[1];
						V_tilde(col, row) = -V_tilde(row, col);
					}
				}
				Eigen::Matrix3d deltaU = mU*U_tilde;
				Eigen::Matrix3d deltaV = V_tilde*mV.transpose();

				dRdF(i, j) = deltaU*mV.transpose() + mU*deltaV;
			}
		}
	
	Tensor3333 lambda_term;
	for(int i =0;i<3;i++)
		for(int j=0;j<3;j++)
		{
			lambda_term(i,j) =
				(dRdF(i,j).transpose()*mF+mR.transpose()*dFdF(i,j)).trace()*mR +
				(mR.transpose()*mF-Eigen::Matrix3d::Identity()).trace()*dRdF(i,j);
		}

	dPdF = mMu*(dFdF - dRdF) + mLambda*lambda_term;
}

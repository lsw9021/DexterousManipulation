#include "BarycentricAttachmentCst.h"
#include <iostream>
using namespace FEM;

void 
BarycentricAttachmentCst::
EvaluatePotentialEnergy(const Eigen::VectorXd& x)
{

}
void 
BarycentricAttachmentCst::
EvaluateGradient(const Eigen::VectorXd& x)
{

}
void 
BarycentricAttachmentCst::
EvaluateHessian(const Eigen::VectorXd& x)
{

}
void 
BarycentricAttachmentCst::
GetPotentialEnergy(double& e)
{

}
void 
BarycentricAttachmentCst::
GetGradient(Eigen::VectorXd& g)
{

}
void 
BarycentricAttachmentCst::
GetHessian(std::vector<Eigen::Triplet<double>>& h_triplets)
{

}
void 
BarycentricAttachmentCst::
EvaluateDVector(const Eigen::VectorXd& x)
{
}
void 
BarycentricAttachmentCst::
GetDVector(int& index,Eigen::VectorXd& d)
{
	d.block<3,1>(3*index,0) = mP;
	index++;
}
void 
BarycentricAttachmentCst::
EvaluateJMatrix(int& index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	Eigen::MatrixXd Ai(3,12);

	for(int i =0;i<4;i++)
		Ai.block<3,3>(0,i*3) = mC[i]*Eigen::Matrix3d::Identity();

	Eigen::MatrixXd MuAiT = mStiffness*Ai.transpose();

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
BarycentricAttachmentCst::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	Eigen::MatrixXd Ai(3,12);

	for(int i =0;i<4;i++)
		Ai.block<3,3>(0,i*3) = mC[i]*Eigen::Matrix3d::Identity();
	Eigen::MatrixXd MuAiTAi = mStiffness*(Ai.transpose()*Ai);

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
BarycentricAttachmentCst::
GetNumHessianTriplets()
{
	return 144;
}
void 
BarycentricAttachmentCst::
AddOffset(int offset)
{
	mi0 +=offset;
	mi1 +=offset;
	mi2 +=offset;
	mi3 +=offset;
}

std::shared_ptr<Cst>
BarycentricAttachmentCst::
Clone()
{
	return nullptr;
}
std::shared_ptr<BarycentricAttachmentCst>
BarycentricAttachmentCst::
Create(const std::string& name,double k,const Eigen::Vector3d& p,const Eigen::Vector4d& c,int i0,int i1,int i2,int i3)
{
	auto bc = new BarycentricAttachmentCst(name,k,p,c,i0,i1,i2,i3);
	return std::shared_ptr<BarycentricAttachmentCst>(bc);
}

BarycentricAttachmentCst::
BarycentricAttachmentCst(const std::string& name,double k,const Eigen::Vector3d& p,const Eigen::Vector4d& c,int i0,int i1,int i2,int i3)
	:Cst(name,k),mC(c),mi0(i0),mi1(i1),mi2(i2),mi3(i3),mP(p)
{
	mE = 0.0;
	mg.setZero();
	mH.setZero();
}
#include "AttachmentCst.h"

using namespace FEM;

AttachmentCst::
AttachmentCst(const std::string& name,double k,int i0,const Eigen::Vector3d& p)
	:Cst(name,k),mi0(i0),mp(p)
{
	mE = 0.0;
	mg.setZero();
	mH.setZero();
	md.setZero();
}
std::shared_ptr<Cst>
AttachmentCst::
Clone()
{
	return Create(mName,mStiffness,mi0,mp);
}
std::shared_ptr<AttachmentCst>
AttachmentCst::
Create(const std::string& name,double k,int i0,const Eigen::Vector3d& p)
{
	auto c = new AttachmentCst(name,k,i0,p);
	return std::shared_ptr<AttachmentCst>(c);
}
void
AttachmentCst::
EvaluatePotentialEnergy(const Eigen::VectorXd& x)
{
	mE = 0.5 * mStiffness * ((x.block<3,1>(mi0*3,0) - mp).squaredNorm());
}
void
AttachmentCst::
EvaluateGradient(const Eigen::VectorXd& x)
{
	mg = mStiffness*(x.block<3,1>(mi0*3,0) - mp);
}
void
AttachmentCst::
EvaluateHessian(const Eigen::VectorXd& x)
{
	Eigen::Vector3d kv;

	kv[0] = mStiffness;
	kv[1] = mStiffness;
	kv[2] = mStiffness;
	mH = kv.asDiagonal();
}
void
AttachmentCst::
GetPotentialEnergy(double& e)
{
	e += mE;
}
void
AttachmentCst::
GetGradient(Eigen::VectorXd& g)
{
	g.block<3,1>(mi0*3,0) += mg;
}
void
AttachmentCst::
GetHessian(std::vector<Eigen::Triplet<double>>& h_triplets)
{
	h_triplets.push_back(Eigen::Triplet<double>(3*mi0+0, 3*mi0+0, mH(0, 0)));
	h_triplets.push_back(Eigen::Triplet<double>(3*mi0+1, 3*mi0+1, mH(1, 1)));
	h_triplets.push_back(Eigen::Triplet<double>(3*mi0+2, 3*mi0+2, mH(2, 2)));
}
void
AttachmentCst::
EvaluateDVector(const Eigen::VectorXd& x)
{
	md = mp;
}
void
AttachmentCst::
GetDVector(int& index,Eigen::VectorXd& d)
{
	d.block<3,1>(index*3,0) = md;
	index++;
}
void
AttachmentCst::
EvaluateJMatrix(int& index, std::vector<Eigen::Triplet<double>>& J_triplets)
{
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+0, 3*index+0, mStiffness));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+1, 3*index+1, mStiffness));
	J_triplets.push_back(Eigen::Triplet<double>(3*mi0+2, 3*index+2, mStiffness));
	index++;
}
void
AttachmentCst::
EvaluateLMatrix(std::vector<Eigen::Triplet<double>>& L_triplets)
{
	L_triplets.push_back(Eigen::Triplet<double>(3*mi0+0, 3*mi0+0, mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(3*mi0+1, 3*mi0+1, mStiffness));
	L_triplets.push_back(Eigen::Triplet<double>(3*mi0+2, 3*mi0+2, mStiffness));
}
int
AttachmentCst::
GetNumHessianTriplets()
{
	return 3;
}

void
AttachmentCst::
AddOffset(int offset)
{
	mi0 += offset;
}
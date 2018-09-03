#include "IKOptimization.h"
#include <Eigen/Geometry>
#include <iostream>
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;
void
IKOptimization::
ChangeSkeleton(const dart::dynamics::SkeletonPtr& skeleton)
{
	mSkeleton = skeleton;
	mSolution.setZero();
	mTargets.clear();
}
void 
IKOptimization::
AddTargetPositions(AnchorPoint ap,const Eigen::Vector3d& target)
{
	bool isAlreadyExist = false;
	for(int i =0;i<mTargets.size();i++)
	{
		if(mTargets[i].first.first ==ap.first){
			isAlreadyExist = true;
			mTargets[i].first.second = ap.second;
			mTargets[i].second = target;
		}
	}
	if(!isAlreadyExist)
		mTargets.push_back(std::make_pair(ap,target));
	// for(int i =0;i<mTargets.size();i++)
	// {
	// 	std::cout<<"Target "<<i<<" : "<<mTargets[i].first.first->getName()<<" -> "<<mTargets[i].second.transpose()<<std::endl;
	// }
}
const std::vector<std::pair<AnchorPoint,Eigen::Vector3d>>&
IKOptimization::
GetTargets()
{
	return mTargets;
}
IKOptimization::
IKOptimization(const SkeletonPtr& skeleton)
	:mSkeleton(skeleton),mSolution(skeleton->getPositions())
{	
	mInitialPositions = mSolution;
	mSkeleton->computeForwardKinematics(true,false,false);
	w_reg = 1E-2;
}
void
IKOptimization::
ClearTarget()
{
	mTargets.resize(0);
}
const Eigen::VectorXd&
IKOptimization::
GetSolution()
{
	return mSolution;
}
void
IKOptimization::
SetSolution(const Eigen::VectorXd& sol)
{
	mSolution = sol;
}
IKOptimization::
~IKOptimization()
{

}

bool					
IKOptimization::
get_nlp_info(	Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
												Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	n = mSkeleton->getNumDofs();
	m = 0;
	nnz_jac_g = 0;
	nnz_h_lag = n;
	index_style = TNLP::C_STYLE;
	return true;
}


bool					
IKOptimization::
get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
												Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) 	
{
	for(int i =0;i<n;i++)
	{
		x_l[i] = mSkeleton->getDof(i)->getPositionLowerLimit();
		x_u[i] = mSkeleton->getDof(i)->getPositionUpperLimit();
	}
	
	return true;
}

bool					
IKOptimization::
get_starting_point(	Ipopt::Index n, bool init_x, Ipopt::Number* x,
													bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
													Ipopt::Index m, bool init_lambda,
													Ipopt::Number* lambda) 
{
	mSavePositions = mSkeleton->getPositions();
	
	if(init_x)
	for(int i =0;i<n;i++){
		x[i] = mSkeleton->getDof(i)->getPosition();
		//x[i] = mInitialPositions[i];
	}

	mSavePositions = mSkeleton->getPositions();
	// std::cout<<"SAVE : "<<mSavePositions.transpose()<<std::endl;
	return true;
}
Eigen::AngleAxisd GetDiff(const Eigen::Quaterniond& diff)
{
	Eigen::AngleAxisd diff1,diff2;
	diff1 = Eigen::AngleAxisd(diff);

	if(diff1.angle()>3.141592)
	{
		diff2.axis() = -diff1.axis();
		diff2.angle() = 3.141592*2 - diff1.angle();	
	}
	else
		diff2 = diff1;
	return diff2;
}
bool					
IKOptimization::
eval_f(	Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) 
{
	Eigen::VectorXd q(n);
	for(int i =0;i<n;i++)
		q[i] = x[i];
	// std::cout<<q.transpose()<<std::endl;
	mSkeleton->setPositions(q);
	mSkeleton->computeForwardKinematics(true,false,false);
	obj_value = 0;
	for(auto& target : mTargets)
	{
		// std::cout<<(target.first.first->getTransform()*target.first.second).transpose()<<std::endl;
		// std::cout<<target.second.transpose()<<std::endl<<std::endl;
		obj_value += 0.5*(target.first.first->getTransform()*target.first.second - target.second).squaredNorm();
	}
	obj_value+= 0.5*w_reg*((q-mInitialPositions).block<6,1>(5,0)).squaredNorm();
	// std::cout<<"OBJ : "<<obj_value<<std::endl;
	return true;
}

bool					
IKOptimization::
eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) 
{
	Eigen::VectorXd q(n),g(n);
	for(int i =0;i<n;i++)
		q[i] = x[i];

	g.setZero();
	// std::cout<<q.transpose()<<std::endl;
	mSkeleton->setPositions(q);
	mSkeleton->computeForwardKinematics(true,false,false);
	for(auto& target: mTargets)
	{
		dart::math::LinearJacobian J = mSkeleton->getLinearJacobian(target.first.first,target.first.second);
		Eigen::MatrixXd J_inv = J.transpose()*(J*J.transpose()).inverse();

		Eigen::Vector3d x_minus_x_target = target.first.first->getTransform()*target.first.second - target.second;

		g += J_inv*x_minus_x_target;

	}
	g.block<6,1>(5,0)+= w_reg*(q-mInitialPositions).block<6,1>(5,0);
	for(int i =0;i<n;i++)
		grad_f[i] = g[i];
	// std::cout<<"grad : "<<g.transpose()<<std::endl;
	return true;
}

bool					
IKOptimization::
eval_g(	Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) 
{
	return true;
}

bool					
IKOptimization::
eval_jac_g( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
											Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
											Ipopt::Number* values) 
{
	return true;
}

bool					
IKOptimization::
eval_h( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
										Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
										bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
										Ipopt::Index* jCol, Ipopt::Number* values) 
{
	int nnz = 0;

	if(values == NULL)
	{
		for(int i=0;i<n;i++)
		{
			iRow[nnz] = i;
			jCol[nnz++] = i;
		}
	}
	else
	{
		for(int i=0;i<n;i++)
		{
			values[nnz++] = 1.0;
		}
	}

	return true;
}

void 					
IKOptimization::
finalize_solution(	Ipopt::SolverReturn status,
													Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
													Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
													Ipopt::Number obj_value,
													const Ipopt::IpoptData* ip_data,
													Ipopt::IpoptCalculatedQuantities* ip_cq) 
{
	for(int i=0;i<n;i++)
		mSolution[i] = x[i];
	mSkeleton->setPositions(mSavePositions);
	mSkeleton->computeForwardKinematics(true,false,false);
}

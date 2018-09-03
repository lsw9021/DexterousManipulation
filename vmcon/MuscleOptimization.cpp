#include "MuscleOptimization.h"
#include "IntegratedWorld.h"
#include "MusculoSkeletalSystem.h"
#include "Ball.h"
#include <IpIpoptData.hpp>
using namespace Ipopt;
MuscleOptimization::
MuscleOptimization(const std::shared_ptr<IntegratedWorld>& iw)
	:mIntegratedWorld(iw),mWeightTracking(1.0),mWeightEffort(0.1),mSparseUpdateCount(0)
{
	int num_muscles =  mIntegratedWorld->GetMusculoSkeletalSystem()->GetNumMuscles();
	int num_muscle_forces =  mIntegratedWorld->GetMusculoSkeletalSystem()->GetNumMuscleForces();
	int dofs 		=  mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();  
	mSolution.resize(dofs +num_muscles);
	mQddDesired.resize(dofs);

	mJt.resize(dofs,(num_muscle_forces)*3);
	mA.resize(num_muscle_forces*3,num_muscles);
	mP.resize(num_muscle_forces*3);

	mM_minus_JtA.resize(dofs,dofs+num_muscles);
	mJtp_minus_c.resize(dofs);	
	mSolution.setZero();
	mQddDesired.setZero();

	mJt.setZero();
	mP.setZero();

	mM_minus_JtA.setZero();
	mJtp_minus_c.setZero();
}
MuscleOptimization::
~MuscleOptimization()
{

}
void
MuscleOptimization::
Update(const Eigen::VectorXd& qdd_desired)
{
	mQddDesired = qdd_desired;

	//Update Jt
	auto& skel = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();
	int dofs = skel->getNumDofs();
	int index = 0;
	mJt.setZero();
	for(auto& muscle : mIntegratedWorld->GetMusculoSkeletalSystem()->GetMuscles())
	{
		
		auto& origin_way_points = muscle->origin_way_points;
		auto& insertion_way_points = muscle->insertion_way_points;
		int n = muscle->GetNumForces();
		int no = origin_way_points.size();
		mJt.block(0,index*3,dofs,muscle->GetNumForces()*3) = muscle->GetJacobianTranspose();
		// std::cout<<mJt<<std::endl<<std::endl;
		// mJt.block(0,(index+0)*3,dofs,3) = skel->getLinearJacobian(muscle->origin_way_points.back().first,muscle->origin_way_points.back().second).transpose();
		// mJt.block(0,(index+1)*3,dofs,3) = skel->getLinearJacobian(muscle->insertion_way_points.back().first,muscle->insertion_way_points.back().second).transpose();
		// for(int j=0;j<no;j++)
		// 	mJt.block(0,(index+j)*3,dofs,3) = skel->getLinearJacobian(muscle->origin_way_points[j].first,muscle->origin_way_points[j].second).transpose();

		// for(int j=no;j<n;j++)
		// 	mJt.block(0,(index+j)*3,dofs,3) = skel->getLinearJacobian(muscle->insertion_way_points[j-no].first,muscle->insertion_way_points[j-no].second).transpose();

		index +=n;
	}
	// std::cout<<std::endl;
	int num_muscles =  mIntegratedWorld->GetMusculoSkeletalSystem()->GetNumMuscles();
	mIntegratedWorld->GetMusculoSkeletalSystem()->ComputeForceDerivative(mIntegratedWorld->GetSoftWorld(),mA);
	mP = mIntegratedWorld->GetMusculoSkeletalSystem()->ComputeForce(mIntegratedWorld->GetSoftWorld());

	mP = mP-mA*mIntegratedWorld->GetMusculoSkeletalSystem()->GetActivationLevels();

	mM_minus_JtA.setZero();
	mM_minus_JtA.block(0,0,dofs,dofs)= mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getMassMatrix();
	mM_minus_JtA.block(0,dofs,dofs,num_muscles)= -mJt*mA;
	mJtp_minus_c = mJt*mP - mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getCoriolisAndGravityForces();
}

const Eigen::VectorXd&
MuscleOptimization::
GetSolution()
{
	return mSolution;
}

void
MuscleOptimization::
UpdateConstraints(const Eigen::VectorXd& act)
{
	// return;
	auto prev_act = mIntegratedWorld->GetMusculoSkeletalSystem()->GetActivationLevels();
	int num_muscles =  mIntegratedWorld->GetMusculoSkeletalSystem()->GetNumMuscles();
	int dofs 		=  mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs(); 
	// if(mSparseUpdateCount%10==0)
	// {
	// 	if( (act-prev_act).norm()>1E-5 || mSparseUpdateCount ==0)
	// 	{
			mIntegratedWorld->GetMusculoSkeletalSystem()->SetActivationLevels(act);
			mIntegratedWorld->GetSoftWorld()->TimeStepping(false);
			//Update A,p
			mIntegratedWorld->GetMusculoSkeletalSystem()->ComputeForceDerivative(mIntegratedWorld->GetSoftWorld(),mA);
			mP = mIntegratedWorld->GetMusculoSkeletalSystem()->ComputeForce(mIntegratedWorld->GetSoftWorld());

			mP = mP-mA*act;
			//Update Cache
			mM_minus_JtA.block(0,0,dofs,dofs)= mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getMassMatrix();
			mM_minus_JtA.block(0,dofs,dofs,num_muscles)= -mJt*mA;
			mJtp_minus_c = mJt*mP - mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getCoriolisAndGravityForces();
	// 	}
	// }
	mSparseUpdateCount++;
}
bool
MuscleOptimization::
get_nlp_info(	Index& n, Index& m, Index& nnz_jac_g,Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	m = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();
	n = m + mIntegratedWorld->GetMusculoSkeletalSystem()->GetNumMuscles();

	nnz_jac_g = n*m;			//g : full matrix
	nnz_h_lag = n;				//H : full matrix

	index_style = TNLP::C_STYLE;
	return true;
}
bool
MuscleOptimization::
get_bounds_info(Index n, Number* x_l, Number* x_u,Index m, Number* g_l, Number* g_u)
{
	for(int i=0;i<m;i++)
	{
		x_l[i] =-1E5;
		x_u[i] =1E5;
	}
	for(int i=m;i<n;i++)
	{
		x_l[i] = 0.0;
		x_u[i] = 1.0;
	}
	for(int i =0;i<m;i++)
		g_l[i]=g_u[i] =0.0;

	return true;
}
bool
MuscleOptimization::
get_starting_point(	Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
	for(int i =0;i<n;i++)
		x[i] = mSolution[i];

	return true;
}
bool
MuscleOptimization::
eval_f(	Index n, const Number* x, bool new_x, Number& obj_value)
{
	double track = 0.0;
	double effort = 0.0;

	int m = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();
	Eigen::VectorXd qdd(m);

	for(int i=0;i<m;i++)
		qdd[i] = x[i];

	track = mWeightTracking*((qdd-mQddDesired).squaredNorm());

	for(int i=m;i<n;i++)
		effort += x[i]*x[i];
	effort *= mWeightEffort;

	obj_value = track + effort;

	return true;
}
bool
MuscleOptimization::
eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	int m = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();

	for(int i=0;i<m;i++)
		grad_f[i] = 2.0*mWeightTracking*(x[i]-mQddDesired[i]);
	for(int i=m;i<n;i++)
		grad_f[i] = 2.0*mWeightEffort*(x[i]);

	return true;
}
bool
MuscleOptimization::
eval_g(	Index n, const Number* x, bool new_x, Index m, Number* g)
{
	Eigen::VectorXd eigen_x(n), activation(n-m),eigen_g(m);
	for(int i=0;i<n;i++)
		eigen_x[i] = x[i];

	activation = eigen_x.tail(n-m);
	UpdateConstraints(activation);
	eigen_g = mM_minus_JtA*eigen_x - mJtp_minus_c;
	for(int i = 0;i<m;i++)
		g[i] = eigen_g[i];

	return true;
}
bool
MuscleOptimization::
eval_jac_g( Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
	int nnz = 0;

	if(values == NULL)
	{
		for(int i =0;i<m;i++)
		{
			for(int j =0;j<n;j++)
			{
				iRow[nnz] = i;
				jCol[nnz++] = j;
			}
		}
	}
	else
	{
		Eigen::VectorXd eigen_x(n), activation(n-m),eigen_g(m);
		for(int i=0;i<n;i++)
			eigen_x[i] = x[i];

		activation = eigen_x.tail(n-m);
		UpdateConstraints(activation);
		for(int i =0;i<m;i++)
		{
			for(int j =0;j<n;j++)
			{
				values[nnz++] = mM_minus_JtA(i,j);
			}
		}
	}

	return true;

}
bool
MuscleOptimization::
eval_h( Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	int nnz = 0;

	if(values == NULL)
	{
		for(int i =0;i<n;i++)
		{
			iRow[nnz] = i;
			jCol[nnz++] = i;
		}

	}
	else
	{
		for(int i =0;i<m;i++)
			values[nnz++] =  2.0*obj_factor*mWeightTracking;
		for(int i=m;i<n;i++)
			values[nnz++] = 2.0*obj_factor*mWeightEffort;
	}

	return true;
}
Ipopt::Index num_iter=0;
void
MuscleOptimization::
finalize_solution(	SolverReturn status,Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda,Number obj_value,const IpoptData* ip_data,IpoptCalculatedQuantities* ip_cq)
{
	num_iter+=ip_data->iter_count();
	// std::cout<<ip_data->iter_count()<<std::endl;
	// std::cout<<num_iter<<std::endl;
	for(int i=0;i<n;i++)
		mSolution[i] = x[i];
	mSparseUpdateCount = 0;
}

#include "HLC.h"
#include "LLC.h"
#include "Record.h"
#include "MusculoSkeletalSystem.h"
#include "Ball.h"
#include "IntegratedWorld.h"
#include "IKOptimization.h"
#include "BezierCurve.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <fstream>

#include <tinyxml.h>

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;


HighLevelController::
HighLevelController(const std::string& weight_path,int maxiteration)
	:mMaxIteration(maxiteration)
{
	std::ifstream param(weight_path);
	param>>w_pos_track>>w_vel_track;
	param.close();
}

void
HighLevelController::
Initialze(
	const Eigen::Vector3d& pos_desired,
		const Eigen::Vector3d& vel_desired,
		const std::shared_ptr<IntegratedWorld>& iw,
		const std::vector<std::pair<AnchorPoint,Eigen::Vector3d>>& ik_targets,
		const std::map<std::string,Eigen::Vector3d>& catch_point,
		int index,const std::string& body_name,
		int next_index,const std::string& next_body_name,bool next_ball_initially_attached,double t_hold)
{
	mReferenceIntegratedWorld = iw;
	mIntegratedWorld = mReferenceIntegratedWorld->Clone();

	mCatchPoint = catch_point;
	mBallIndex = index;
	mBody = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getBodyNode(body_name);
	mNextBallIndex = next_index;
	mNextBody = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getBodyNode(next_body_name);
	mNextBallInitiallyAttached = next_ball_initially_attached;
	mT_Hold = t_hold;

	std::cout<<"p_t : "<<pos_desired.transpose()<<std::endl;
	std::cout<<"v_t : "<<vel_desired.transpose()<<std::endl;

	mBallTargetPosition = pos_desired;
	mBallTargetVelocity = vel_desired;


	if(GetRawPtr(mIKOptimization)==nullptr)
	{	
		mIKOptimization = new IKOptimization(mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton());
		mIKSolver = new IpoptApplication();
		mIKSolver->Options()->SetStringValue("jac_c_constant", "yes");
		mIKSolver->Options()->SetStringValue("hessian_constant", "yes");
		mIKSolver->Options()->SetIntegerValue("print_level", 2);
		mIKSolver->Options()->SetIntegerValue("max_iter", 100);
		mIKSolver->Options()->SetNumericValue("tol", 1e-4);
	}
	else
		static_cast<IKOptimization*>(GetRawPtr(mIKOptimization))->ChangeSkeleton(mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton());
	

	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));

	for(auto target : ik_targets)
	{
		ik->AddTargetPositions(std::make_pair(mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getBodyNode(target.first.first->getName()),target.first.second),target.second);
	}

	mIKSolver->Initialize();
	mIKSolver->OptimizeTNLP(mIKOptimization);

	mLLC = std::make_shared<LowLevelController>(mIntegratedWorld);
}

std::vector<Eigen::VectorXd>
HighLevelController::
Solve(const Eigen::VectorXd& x0)
{
	Eigen::VectorXd x_star = x0;
	std::vector<Eigen::VectorXd> ret;
	
	
	for(int i =0;i<mMaxIteration;i++)
	{
		double f = EvalEnergy(x_star,true);

		ret = mReferenceMotions;
		if(f<1E-4)
			break;
	
	
		Eigen::VectorXd fx = EvalGradient(x_star);
	
		Eigen::VectorXd dir = -fx;
		double alpha = 0.2+mBallTargetVelocity[1]*0.05;
		Eigen::VectorXd x_next;
		
		bool is_updated = false;
	
		for(int j=0;j<10;j++)
		{
			x_next = x_star + alpha*dir;
			double f_next = EvalEnergy(x_next);
	
			if(f_next<f)
			{
				f = f_next;
				is_updated =true;
				break;
			}
			alpha*=0.5;
		}
	
		if(is_updated)
			x_star = x_next;

		if(!is_updated)// || f<0.9E-4)
			break;
	}
	return ret;
}

double
HighLevelController::
EvalEnergy(const Eigen::VectorXd& x,bool bake)
{
	std::vector<Eigen::VectorXd> motions;
	GenerateMotions(x,motions);

	Simulate(motions,bake);
	mReferenceMotions = motions;
	Eigen::Vector3d ball_pos = mIntegratedWorld->GetBalls()[mBallIndex]->GetPosition();
	Eigen::Vector3d ball_vel = mIntegratedWorld->GetBalls()[mBallIndex]->GetVelocity();
	if(bake)
		std::cout<<"velocity : "<<mIntegratedWorld->GetBalls()[mBallIndex]->GetVelocity().transpose()<<std::endl;
	double val = 0.0;
	val += 0.5*w_pos_track*(mBallTargetPosition - ball_pos).squaredNorm();
	val += 0.5*w_vel_track*(mBallTargetVelocity - ball_vel).squaredNorm();
	return val;
}

Eigen::VectorXd
HighLevelController::
EvalGradient(const Eigen::VectorXd& x)
{
	Eigen::VectorXd fx(x.rows());
	fx.setZero();

	Eigen::VectorXd x_i;
	double delta = 0.01;
	double fx_minus,fx_plus;
	double x_i_plus,x_i_minus;
	for(int i = 0;i<x.rows();i++)
	{
		if(i>2)
		{
			x_i = x;
			fx_minus =0;
			fx_plus =0;
			x_i_minus = x[i] - delta;
			x_i_plus = x[i] + delta;
			
			x_i[i] = x_i_minus;
			fx_minus = EvalEnergy(x_i);
			x_i[i] = x_i_plus;
			fx_plus = EvalEnergy(x_i);

			fx[i] = (fx_plus - fx_minus)/(x_i_plus - x_i_minus);
		}
		else
		{
			fx[i] = 0.0;
		}
	}
	fx.block<3,1>(3,0) *= 2.0;
	return fx;
}

void
HighLevelController::
GenerateMotions(const Eigen::VectorXd& x,std::vector<Eigen::VectorXd>& motions)
{
	auto& skel = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();
	auto& ball = mIntegratedWorld->GetBalls()[mBallIndex];

	BezierCurve bc;
	bc.Initialize(x.block<3,1>(0,0),x.block<3,1>(3,0),x.block<3,1>(6,0),1);

	std::vector<std::pair<Eigen::VectorXd,double>> coarse_motions;

	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	mIntegratedWorld->SetFrom(mReferenceIntegratedWorld.get());
	ik->SetSolution(mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getPositions());
	auto save_target = ik->GetTargets();
	Eigen::Vector3d p_hb = ball->GetPosition() - mBody->getCOM();
	Eigen::Isometry3d T_hand = mBody->getTransform();
	p_hb = T_hand.linear().inverse()*p_hb;
	AnchorPoint ap = std::make_pair(mBody,p_hb);
	for(int i =0;i<11;i++)
	{
		double tt = ((double)i)/((double)10);
		
		Eigen::Vector3d p_ee = bc.GetPosition(tt);
		ik->AddTargetPositions(ap,p_ee);
		mIKSolver->Initialize();
		mIKSolver->OptimizeTNLP(mIKOptimization);
		Eigen::VectorXd sol = ik->GetSolution();

		coarse_motions.push_back(std::make_pair(sol,tt*mT_Hold));
	}
	
	for(auto& target : save_target){
		ik->AddTargetPositions(target.first,target.second);
	}

	Eigen::VectorXd p;
	double time_elapsed = 0;
	while(true)
	{
		int k =0,k1 =0;
		for(int i =0;i<coarse_motions.size();i++)
		{
			if(coarse_motions[i].second<time_elapsed)
				k=i;
		}
		if(k ==coarse_motions.size()-1)
			break;
		k1 = k+1;
		double t = time_elapsed-coarse_motions[k].second;
		double dt = coarse_motions[k1].second-coarse_motions[k].second;

		p = (1.0-t/dt)*(coarse_motions[k].first) + (t/dt)*(coarse_motions[k1].first);

		motions.push_back(p);
		time_elapsed+= mIntegratedWorld->GetSoftWorld()->GetTimeStep();
	}

}

void
HighLevelController::
Simulate(const std::vector<Eigen::VectorXd>& motions,bool bake)
{
	mIntegratedWorld->SetFrom(mReferenceIntegratedWorld.get());
	auto& balls = mIntegratedWorld->GetBalls();

	Eigen::VectorXd target_positions,target_velocities;
	
	if(mNextBallInitiallyAttached)
		balls[mNextBallIndex]->Attach(mIntegratedWorld->GetRigidWorld(),mNextBody);	
	else
		balls[mNextBallIndex]->Release(mIntegratedWorld->GetRigidWorld());

	for(int t = 0;t<motions.size();t++)
	{
		target_positions = motions[t];
		if(t==0)
			target_velocities = (motions[t+1]-motions[t])/mIntegratedWorld->GetSoftWorld()->GetTimeStep();
		else
			target_velocities = (motions[t]-motions[t-1])/mIntegratedWorld->GetSoftWorld()->GetTimeStep();

		Eigen::VectorXd u = mLLC->Solve(target_positions,target_velocities);
		mIntegratedWorld->TimeStepping(u,bake);

		if(balls[mNextBallIndex]->IsReleased()) // check if already attached.
		{
			Eigen::Vector3d body_position = mNextBody->getTransform()*mCatchPoint[mNextBody->getName()];
			Eigen::Vector3d ball_position = balls[mNextBallIndex]->GetPosition();
			if((body_position-ball_position).norm()<5E-2){
				balls[mNextBallIndex]->Attach(mIntegratedWorld->GetRigidWorld(),mNextBody);
			}
		}
	}

}

#include "FSM.h"
#include "JugglingInfo.h"
#include "MusculoSkeletalSystem.h"
#include "IntegratedWorld.h"
#include "IKOptimization.h"
#include "Ball.h"
#include "BezierCurve.h"

#include "HLC.h"
#include "LLC.h"

using namespace dart::dynamics;
using namespace dart::simulation;

using namespace Ipopt;
FSM::
FSM(const std::shared_ptr<IntegratedWorld>& iw,const std::vector<int>& sequences,int ball_size)
	:mIntegratedWorld(iw),mJugglingInfo(new JugglingInfo(sequences,ball_size))
{
	const auto& skel = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();
	auto& balls = mIntegratedWorld->GetBalls();
	Eigen::Isometry3d T_hand = skel->getBodyNode("HandR")->getTransform();
	Eigen::Isometry3d T_ball = balls[0]->GetSkeleton()->getBodyNode(0)->getTransform();
	Eigen::Vector3d p_hand = T_hand.translation();
	Eigen::Vector3d p_ball = T_ball.translation();

	Eigen::Vector3d left,right;
	right = T_hand.linear().inverse()*(p_ball-p_hand);
	right[2] -= 0.03;
	left = right;
	left[0] = -left[0];

	mCatchPoint.insert(std::make_pair("HandR",right));
	mCatchPoint.insert(std::make_pair("HandL",left));

	mIKOptimization = new IKOptimization(skel);

	mIKSolver = new IpoptApplication();
	mIKSolver->Options()->SetStringValue("jac_c_constant", "yes");
	mIKSolver->Options()->SetStringValue("hessian_constant", "yes");
	mIKSolver->Options()->SetIntegerValue("print_level", 2);
	mIKSolver->Options()->SetIntegerValue("max_iter", 100);
	mIKSolver->Options()->SetNumericValue("tol", 1e-4);

	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));
	ik->SetSolution(skel->getPositions());
	Eigen::Vector3d l_loc = skel->getBodyNode("HandL")->getTransform()*Eigen::Vector3d::Zero();
	Eigen::Vector3d r_loc = skel->getBodyNode("HandR")->getTransform()*Eigen::Vector3d::Zero();

	ik->AddTargetPositions(std::make_pair(skel->getBodyNode("HandL"),Eigen::Vector3d::Zero()),l_loc);
	ik->AddTargetPositions(std::make_pair(skel->getBodyNode("HandR"),Eigen::Vector3d::Zero()),r_loc);

	mHandPoint.insert(std::make_pair("HandR",r_loc));
	mHandPoint.insert(std::make_pair("HandL",l_loc));

	mPhase = 0;
	mHLC = std::make_shared<HighLevelController>("../vmcon/export/HLC_Param.txt",15);
	
	for(int i =0;i<100;i++)
		mMotions.push_back(ik->GetSolution());
}

void
FSM::
GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v)
{
	auto& skel = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();
	auto ball = mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()];
	auto bn_from = skel->getBodyNode(mJugglingInfo->From());
	auto bn_to = skel->getBodyNode(mJugglingInfo->To());
	//Check need to update.
	bool need_update = false;
	//Check Catch Phase Finished.
	if(mPhase == 0)
	{
		Eigen::Vector3d body_position = bn_from->getTransform()*mCatchPoint[bn_from->getName()];
		Eigen::Vector3d ball_position = ball->GetPosition();

		if(ball->IsReleased())
			GenerateCatchMotions();

		if((body_position-ball_position).norm()<5E-2 || !ball->IsReleased())
		{
			ball->Attach(mIntegratedWorld->GetRigidWorld(),bn_from);
			need_update = true;
			mPhase = 1;
		}
	}

	//Look opposite Hand
	mJugglingInfo->CountPlusPlus();
	ball = mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()];
	bn_from = skel->getBodyNode(mJugglingInfo->From());

	if(ball->IsReleased()) // check if already attached.
	{
		Eigen::Vector3d body_position = bn_from->getTransform()*mCatchPoint[bn_from->getName()];
		Eigen::Vector3d ball_position = ball->GetPosition();

		if((body_position-ball_position).norm()<5E-2){
			ball->Attach(mIntegratedWorld->GetRigidWorld(),bn_from);
		}
	}
	mJugglingInfo->CountMinusMinus();
		
	//Check Swing Phase Finished.
	ball = mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()];
	bn_from = skel->getBodyNode(mJugglingInfo->From());
	if(mPhase ==1 &&!need_update)
	{
		//Swing Phase End.
		if(mCount == mMotions.size())
		{
			ball->Release(mIntegratedWorld->GetRigidWorld());
			std::cout<<"Released Velocity : "<<mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()]->releasedVelocity.transpose()<<std::endl;

			//Look Ahead for targeting next ball. (count+2 -> next target)
			mJugglingInfo->CountPlusPlus();
			mJugglingInfo->CountPlusPlus();
			if(mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()]->IsReleased())
				GenerateCatchMotions();
			//to next state.
			mJugglingInfo->CountMinusMinus();
			
			
			need_update = true;
			mPhase = 0;
		}
	}

	if(need_update)
		if(mPhase == 0)
			GenerateCatchMotions();
		else
			GenerateSwingMotions();

	p = mMotions[mCount];
	if(mCount == 0)
		v = (mMotions[mCount+1] - mMotions[mCount])/mIntegratedWorld->GetSoftWorld()->GetTimeStep();
	else
		v = (p-mMotions[mCount-1])/mIntegratedWorld->GetSoftWorld()->GetTimeStep();

	mCount++;
}

void
FSM::
GenerateCatchMotions()
{
	auto& skel = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();
	auto& ball = mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()];
	auto bn_from = skel->getBodyNode(mJugglingInfo->From());
	auto bn_to = skel->getBodyNode(mJugglingInfo->To());

	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));

	Eigen::Vector3d target;

	ball->ComputeFallingPosition(mHandPoint[bn_from->getName()][1],target);

	//No Solution
	if(target.norm()<1E-6){
		mMotions.push_back(mMotions.back());
		return;
	}

	AnchorPoint ap = std::make_pair(bn_from,mCatchPoint[bn_from->getName()]);

	ik->AddTargetPositions(ap,target);
	Eigen::VectorXd save_position = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getPositions();

	mIKSolver->Initialize();
	mIKSolver->OptimizeTNLP(mIKOptimization);	
	
	mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->setPositions(save_position);
	mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->computeForwardKinematics(true,false,false);
	
	mMotions.clear();
	for(int i =0;i<100;i++)
		mMotions.push_back(ik->GetSolution());
	mCount = 0;
}
void
FSM::
GenerateSwingMotions()
{
	auto& skel = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();
	auto& ball = mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()];
	auto bn_from = skel->getBodyNode(mJugglingInfo->From());
	auto bn_to = skel->getBodyNode(mJugglingInfo->To());
	double t_hold = mJugglingInfo->GetT_hold();
	double t_free = mJugglingInfo->GetT_free();

	Eigen::Vector3d dir = bn_to->getCOM() - ball->GetPosition();	

	Eigen::Vector3d p_des = ball->GetPosition();
	Eigen::Vector3d target = mHandPoint[bn_to->getName()]; //bn_to->getTransform()*mCatchPoint[bn_to->getName()];
	p_des[0] *=0.6;
	target[0] *= 1.1;
	Eigen::Vector3d v_des = mJugglingInfo->GetTargetVelocity(p_des,target);

	//Initial Guess
	Eigen::VectorXd x0(9);
	x0.block<3,1>(0,0) = p_des;
	x0.block<3,1>(3,0) = p_des;
	x0.block<3,1>(6,0) = p_des;

	//Get Next ball info
	mJugglingInfo->CountPlusPlus();
	int next_index = mJugglingInfo->GetBallIndex();
	std::string next_body = mJugglingInfo->From();	
	bool next_ball_initially_attached = !mIntegratedWorld->GetBalls()[mJugglingInfo->GetBallIndex()]->IsReleased();
	mJugglingInfo->CountMinusMinus();

	std::string body = mJugglingInfo->From();
	IKOptimization* ik = static_cast<IKOptimization*>(GetRawPtr(mIKOptimization));

	mHLC->Initialze(
		p_des,v_des,
		mIntegratedWorld,
		ik->GetTargets(),
		mCatchPoint,
		mJugglingInfo->GetBallIndex(),body,
		next_index,next_body,next_ball_initially_attached,
		mJugglingInfo->GetT_hold());

	mMotions = mHLC->Solve(x0);

	mCount = 0;
}

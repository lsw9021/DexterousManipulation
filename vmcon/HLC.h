#ifndef __HIGH_LEVEL_CONTROLLER_H__
#define __HIGH_LEVEL_CONTROLLER_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

#include "fem/fem.h"

#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
class IntegratedWorld;
class Controller;
class LowLevelController;
class HighLevelController
{
public:
	HighLevelController(const std::string& weight_path,int max_iteration = 10);

	void Initialze(
		const Eigen::Vector3d& pos_desired,
		const Eigen::Vector3d& vel_desired,
		const std::shared_ptr<IntegratedWorld>& iw,
		const std::vector<std::pair<AnchorPoint,Eigen::Vector3d>>& ik_targets,
		const std::map<std::string,Eigen::Vector3d>& catch_point,
		int index,const std::string& body_name,
		int next_index,const std::string& next_body_name,bool next_ball_initially_attached,double t_hold);

	std::vector<Eigen::VectorXd> Solve(const Eigen::VectorXd& x0);
	
private:
	double EvalEnergy(const Eigen::VectorXd& x,bool bake = false);
	Eigen::VectorXd EvalGradient(const Eigen::VectorXd& x);
	void GenerateMotions(const Eigen::VectorXd& x,std::vector<Eigen::VectorXd>& motions);
	void Simulate(const std::vector<Eigen::VectorXd>& motions,bool bake = false);

	std::shared_ptr<LowLevelController>			mLLC;
	std::shared_ptr<IntegratedWorld>			mIntegratedWorld,mReferenceIntegratedWorld;

	Ipopt::SmartPtr<Ipopt::TNLP>			 	mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mIKSolver;
	std::vector<Eigen::VectorXd>				mReferenceMotions;

	double 										w_pos_track,w_vel_track;
	
	double										mT_Hold;
	std::map<std::string,Eigen::Vector3d>		mCatchPoint;
	int 										mBallIndex;
	dart::dynamics::BodyNode* 					mBody;
	int 										mNextBallIndex;
	dart::dynamics::BodyNode* 					mNextBody;
	bool 										mNextBallInitiallyAttached;

	Eigen::Vector3d								mBallTargetPosition;
	Eigen::Vector3d								mBallTargetVelocity;

	int 										mMaxIteration;
};

#endif
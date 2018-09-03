#ifndef __FSM_H__
#define __FSM_H__

#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

#include "fem/fem.h"

#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>

typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
class JugglingInfo;
class IntegratedWorld;
class HighLevelController;
class FSM
{
public:
	FSM(const std::shared_ptr<IntegratedWorld>& iw,const std::vector<int>& sequences,int ball_size);

	void GetMotion(Eigen::VectorXd& p,Eigen::VectorXd& v);
private:
	void GenerateCatchMotions();
	void GenerateSwingMotions();

	std::shared_ptr<IntegratedWorld>			mIntegratedWorld;

//For Juggling
	std::shared_ptr<JugglingInfo> 				mJugglingInfo;
	int 										mCount;
	int 										mPhase; // 0: catch, 1:swing
	std::vector<Eigen::VectorXd> 				mMotions;
//For Swing Phase
	std::shared_ptr<HighLevelController>		mHLC;

//For Catch Phase
	Ipopt::SmartPtr<Ipopt::TNLP>			 	mIKOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 	mIKSolver;
	std::map<std::string,Eigen::Vector3d>		mCatchPoint;
	std::map<std::string,Eigen::Vector3d>		mHandPoint;
};

#endif
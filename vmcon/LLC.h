#ifndef __LowLevelCONTROLLER__H__
#define __LowLevelCONTROLLER__H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

#include "fem/fem.h"

#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>	

class MuscleOptimization;
class IntegratedWorld;
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;

class LowLevelController
{
public:
	LowLevelController(const std::shared_ptr<IntegratedWorld>& iw);

	void SetKp(double kp);
	Eigen::VectorXd Solve(const Eigen::VectorXd& p_d,const Eigen::VectorXd& v_d);
private:
	std::shared_ptr<IntegratedWorld> 					mIntegratedWorld;
	
	Eigen::VectorXd 									mKp,mKv,mKjt;
	Eigen::VectorXd										mTargetPositions;
	Eigen::VectorXd										mTargetVelocities;

#ifdef USE_MUSCLE
	Ipopt::SmartPtr<Ipopt::TNLP> 			 			mMuscleOptimization;
	Ipopt::SmartPtr<Ipopt::IpoptApplication> 			mMuscleOptimizationSolver;

	Eigen::VectorXd ComputeActivationLevels(const Eigen::VectorXd& qdd_desired);
#endif
	Eigen::VectorXd ComputePDForces();
};


#endif
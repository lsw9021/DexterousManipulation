#ifndef __BALL_H__
#define __BALL_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
class Ball
{
public:
	dart::dynamics::SkeletonPtr 				skeleton;
	dart::constraint::WeldJointConstraintPtr 	constraint;

	bool 										isReleased;

	Eigen::Vector3d								releasedPoint;
	Eigen::Vector3d								releasedVelocity;
public:
	Ball(const dart::dynamics::SkeletonPtr& skel);

	Eigen::Vector3d GetPosition();
	Eigen::Vector3d GetVelocity();
	void Release(const dart::simulation::WorldPtr& world);
	void Attach(const dart::simulation::WorldPtr& world,dart::dynamics::BodyNode* bn);
	void ComputeFallingPosition(double h,Eigen::Vector3d& fp);

	dart::dynamics::SkeletonPtr GetSkeleton(){return skeleton;};
	dart::constraint::WeldJointConstraintPtr GetConstraint(){return constraint;};
	bool IsReleased(){return isReleased;};
};

#endif
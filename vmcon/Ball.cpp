#include "Ball.h"
using namespace dart::dynamics;
using namespace dart::simulation;
Ball::
Ball(const SkeletonPtr& skel)
	:isReleased(true),constraint(nullptr),skeleton(skel),releasedPoint(Eigen::Vector3d::Zero()),releasedVelocity(Eigen::Vector3d::Zero())
{
	
}

void
Ball::
ComputeFallingPosition(double h,Eigen::Vector3d& fp)
{
	Eigen::Vector3d p = skeleton->getBodyNode(0)->getCOM();
	Eigen::Vector3d v = skeleton->getBodyNode(0)->getCOMLinearVelocity();

	double dx = h-p[1];
	double g = -9.8;
	double v2_plus_2gdx = v[1]*v[1] + 2.0*g*dx;
	if(v2_plus_2gdx<0)
	{
		// std::cout<<"no solution"<<std::endl;
		fp.setZero();
		return;
	}

	double t1 = (-v[1] - sqrt(v2_plus_2gdx))/g;
	double t2 = (-v[1] + sqrt(v2_plus_2gdx))/g;
	double t = (t1>t2?t1:t2);
	if(t<0.0)
		t=0.05;

	fp = p+t*v;
	if(fp[1]>h)
	fp[1] = h;
}
Eigen::Vector3d
Ball::
GetPosition()
{
	return skeleton->getBodyNode(0)->getCOM();
}
Eigen::Vector3d
Ball::
GetVelocity()
{
	return skeleton->getBodyNode(0)->getCOMLinearVelocity();
}
void
Ball::
Release(const dart::simulation::WorldPtr& world)
{
	if(!isReleased){
		world->getConstraintSolver()->removeConstraint(constraint);
		constraint = nullptr;
		isReleased = true;
		releasedPoint = skeleton->getBodyNode(0)->getCOM();
		releasedVelocity = skeleton->getBodyNode(0)->getCOMLinearVelocity();
	}
}

void
Ball::
Attach(const dart::simulation::WorldPtr& world,dart::dynamics::BodyNode* bn)
{
	if(isReleased)
	{
		constraint.reset();
		constraint = std::make_shared<dart::constraint::WeldJointConstraint>(skeleton->getBodyNode(0),bn);
		isReleased = false;
		world->getConstraintSolver()->addConstraint(constraint);
	}
}
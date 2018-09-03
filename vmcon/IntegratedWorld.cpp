#include "IntegratedWorld.h"
#include "DART_helper.h"
#include "MusculoSkeletalSystem.h"
#include "Record.h"
#include "Ball.h"
#include <fstream>
#include <sstream>
#include <tinyxml.h>
#include <boost/filesystem.hpp>
 
using namespace FEM;
using namespace dart::dynamics;
using namespace dart::simulation;

IntegratedWorld::
IntegratedWorld()
{

}
void
IntegratedWorld::
TimeStepping(const Eigen::VectorXd& u,bool bake)
{
	// mMusculoSkeletalSystem->TransformAttachmentPoints();
	// mSoftWorld->TimeStepping(false);
#ifdef USE_MUSCLE
	if(u.rows()!=0)
		mMusculoSkeletalSystem->SetActivationLevels(u);
	mMusculoSkeletalSystem->TransformAttachmentPoints();
	mSoftWorld->TimeStepping();
#endif
	
	double nn = mSoftWorld->GetTimeStep() / mRigidWorld->getTimeStep();
	for(int i =0; i<nn;i++)
	{
#ifdef USE_MUSCLE
		mMusculoSkeletalSystem->ApplyForcesToSkeletons(mSoftWorld);
#else
		if(u.rows()!=0)
		mMusculoSkeletalSystem->GetSkeleton()->setForces(
			mMusculoSkeletalSystem->GetSkeleton()->getMassMatrix()*u +
			mMusculoSkeletalSystem->GetSkeleton()->getCoriolisAndGravityForces());
#endif
		mRigidWorld->step();
	}
	
	if(bake)
	{
		mRecords.push_back(std::make_shared<Record>());
		auto rec = mRecords.back();
		rec->Set(mRigidWorld,mSoftWorld,mMusculoSkeletalSystem);
	}
	
}

void
IntegratedWorld::
Initialize()
{
	//Init Soft world
	mSoftWorld = FEM::World::Create(
		FEM::IntegrationMethod::PROJECTIVE_QUASI_STATIC,	//Integration Method
		// FEM::IntegrationMethod::QUASI_STATIC,	//Integration Method
		// FEM::IntegrationMethod::PROJECTIVE_DYNAMICS,	//Integration Method
		1.0/500.0,							//Time Step
		10,								//Max Iteration
		Eigen::Vector3d(0,-9.81,0),					//Gravity
		0.999								//Damping
		);
	//Init Rigid world
	mRigidWorld = std::make_shared<dart::simulation::World>();
	mRigidWorld->setGravity(Eigen::Vector3d(0,-9.81,0));
	// mRigidWorld->setGravity(Eigen::Vector3d(0,0,0));
	mRigidWorld->setTimeStep(1.0/1000.0);
	mRigidWorld->checkCollision();

	//Init MusculoSkeletalSystem
	mMusculoSkeletalSystem = std::make_shared<MusculoSkeletalSystem>();
	MakeSkeleton(mMusculoSkeletalSystem);
	MakeMuscles("../vmcon/export/muscle_params.xml",mMusculoSkeletalSystem);
	mMusculoSkeletalSystem->Initialize(mSoftWorld,mRigidWorld);

	MakeBalls(mRigidWorld,mMusculoSkeletalSystem,mBalls,5);
	dart::constraint::WeldJointConstraintPtr cons1 = 
		std::make_shared<dart::constraint::WeldJointConstraint>(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandR"));
	dart::constraint::WeldJointConstraintPtr cons2 = 
		std::make_shared<dart::constraint::WeldJointConstraint>(mMusculoSkeletalSystem->GetSkeleton()->getBodyNode("HandL"));
	// mRigidWorld->getConstraintSolver()->addConstraint(cons1);
	// mRigidWorld->getConstraintSolver()->addConstraint(cons2);
	mSoftWorld->Initialize();
}

std::shared_ptr<IntegratedWorld>
IntegratedWorld::
Clone()
{
	auto new_iw = std::make_shared<IntegratedWorld>();
	new_iw->Initialize();
	new_iw->SetFrom(this);
	
	return new_iw;
}

void
IntegratedWorld::
SetFrom(IntegratedWorld* other)
{
	this->mSoftWorld->SetPositions(other->GetSoftWorld()->GetPositions());

	for(int i =0;i<mBalls.size();i++)
		this->mBalls[i]->Release(mRigidWorld);

	for(int i =0;i<mRigidWorld->getNumSkeletons();i++)
	{
		this->mRigidWorld->getSkeleton(i)->setPositions(other->GetRigidWorld()->getSkeleton(i)->getPositions());
		this->mRigidWorld->getSkeleton(i)->setVelocities(other->GetRigidWorld()->getSkeleton(i)->getVelocities());
		this->mRigidWorld->getSkeleton(i)->computeForwardKinematics(true,false,false);
	}
	this->mMusculoSkeletalSystem->TransformAttachmentPoints();
	this->mMusculoSkeletalSystem->SetActivationLevels(other->GetMusculoSkeletalSystem()->GetActivationLevels());

	for(int i =0;i<mBalls.size();i++)
	{
		if(!other->GetBalls()[i]->IsReleased())
		{
			this->mBalls[i]->Attach(
				this->mRigidWorld,
				this->mMusculoSkeletalSystem->GetSkeleton()->getBodyNode(other->mBalls[i]->GetConstraint()->getBodyNode2()->getName()));
		}
	}
}
#ifndef __INTEGRATED_WORLD_H__
#define __INTEGRATED_WORLD_H__

#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

#include "fem/fem.h"

class Ball;
class MusculoSkeletalSystem;
class Record;

class IntegratedWorld
{
public:
	IntegratedWorld();


	void Initialize();
	void TimeStepping(const Eigen::VectorXd& u = Eigen::VectorXd::Zero(0),bool bake = true);

	std::shared_ptr<IntegratedWorld> 	Clone();
	void								SetFrom(IntegratedWorld* other);

	const FEM::WorldPtr& GetSoftWorld(){return mSoftWorld;};	
	const dart::simulation::WorldPtr& GetRigidWorld() {return mRigidWorld;};
	const std::shared_ptr<MusculoSkeletalSystem>& GetMusculoSkeletalSystem(){return mMusculoSkeletalSystem;};
	const std::vector<std::shared_ptr<Ball>>& GetBalls(){return mBalls;};
	const std::vector<std::shared_ptr<Record>>& GetRecords(){return mRecords;};
	void AddRecord(const std::shared_ptr<Record>& rec) {mRecords.push_back(rec);};
private:
	std::vector<std::shared_ptr<Record>> 				mRecords;
	FEM::WorldPtr 										mSoftWorld;
	dart::simulation::WorldPtr 							mRigidWorld;
	std::shared_ptr<MusculoSkeletalSystem> 				mMusculoSkeletalSystem;
	std::vector<std::shared_ptr<Ball>>					mBalls;
};



#endif
#ifndef __MUSCULO_SKELETAL_SYSTEM_H__
#define __MUSCULO_SKELETAL_SYSTEM_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include "fem/fem.h"

typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
Eigen::Vector3d GetPoint(const AnchorPoint& ap);
class Ball;

struct Muscle
{
	Muscle(){};
	int GetNumForces();
	Eigen::MatrixXd GetJacobianTranspose();
	std::vector<int>									origin_force_indices,insertion_force_indices;
	void TransferForce(const Eigen::Vector3d& f_origin_tilda,const Eigen::Vector3d& f_insertion_tilda,Eigen::VectorXd& f);
	void SetActivationLevel(double a);
	void Initialize();
	std::string 										name;
	FEM::MeshPtr										mesh;
	std::vector<AnchorPoint>							origin_way_points,insertion_way_points;
	double												activation_level;
	

	FEM::AttachmentCstPtr								origin,insertion;
	std::vector<FEM::LinearMuscleCstPtr>				muscle_csts;
	std::vector<FEM::CstPtr>							csts;
	
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
};

class MusculoSkeletalSystem
{
public:
	MusculoSkeletalSystem();
	void AddMuscle(
		const std::string& name,
		const std::vector<AnchorPoint>& origin,
		const std::vector<AnchorPoint>& insertion,
		int origin_index,int insertion_index,
		const Eigen::Vector3d& fiber_direction,
		const FEM::MeshPtr& mesh);

	void Initialize(const FEM::WorldPtr& soft_world,const dart::simulation::WorldPtr& rigid_world);

	void SetActivationLevels(const Eigen::VectorXd& a);
	void TransformAttachmentPoints();
	void ApplyForcesToSkeletons(const std::shared_ptr<FEM::World>& soft_world);

	void ComputeForceDerivative(const FEM::WorldPtr& world,Eigen::SparseMatrix<double>& J);
	Eigen::VectorXd ComputeForce(const FEM::WorldPtr& world);
	
	int 									GetNumMuscles() 		{return mMuscles.size();}
	std::vector<std::shared_ptr<Muscle>>&	GetMuscles()			{return mMuscles;}
	dart::dynamics::SkeletonPtr&			GetSkeleton()			{return mSkeleton;}
	Eigen::VectorXd 						GetActivationLevels()	{return mActivationLevels;}
	int	 									GetNumMuscleForces();
private:
	//Muscles and Skeleton
	std::vector<std::shared_ptr<Muscle>>	mMuscles;
	dart::dynamics::SkeletonPtr 			mSkeleton;

	//Material Properties
	double	mTendonStiffness;
	double	mMuscleStiffness;
	double	mYoungsModulus;
	double	mPoissonRatio;

	Eigen::VectorXd							mActivationLevels;

	std::vector<FEM::LinearMuscleCstPtr>	mAllMuscleConstraints;
};

void MakeMuscles(const std::string& path,std::shared_ptr<MusculoSkeletalSystem>& ms);
void MakeSkeleton(std::shared_ptr<MusculoSkeletalSystem>& ms);
void MakeBalls(dart::simulation::WorldPtr& world,const std::shared_ptr<MusculoSkeletalSystem>& ms,std::vector<std::shared_ptr<Ball>>& ball,int num);
#endif
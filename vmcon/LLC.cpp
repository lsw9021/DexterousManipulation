#include "LLC.h"
#include "IKOptimization.h"
#include "MuscleOptimization.h"
#include "IntegratedWorld.h"
#include "MusculoSkeletalSystem.h"
#include "Ball.h"
// #include "FSM/Machine.h"

using namespace FEM;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace Ipopt;

LowLevelController::
LowLevelController(const std::shared_ptr<IntegratedWorld>& iw)
	:mIntegratedWorld(iw)
{
	dart::math::seedRand();
	int dof = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();
	double k = 500;

	mKp = Eigen::VectorXd::Constant(dof,k);
	mKv = Eigen::VectorXd::Constant(dof,2*sqrt(k));
	mKjt = Eigen::VectorXd::Constant(dof,0.1);
	for(int i =0;i<6;i++)
	{
		mKp[dof-1-i] = 2.0*mKp[dof-1-i];	
		mKv[dof-1-i] = sqrt(2.0)*mKv[dof-1-i];
	}
	
	mTargetPositions = Eigen::VectorXd::Constant(dof,0.0);
	mTargetVelocities = Eigen::VectorXd::Constant(dof,0.0);

#ifdef USE_MUSCLE
	mMuscleOptimization = new MuscleOptimization(mIntegratedWorld);
	mMuscleOptimizationSolver = new IpoptApplication();
	
	// mMuscleOptimizationSolver->Options()->SetStringValue("mu_strategy", "adaptive");
	mMuscleOptimizationSolver->Options()->SetStringValue("jac_c_constant", "yes");
	mMuscleOptimizationSolver->Options()->SetStringValue("hessian_constant", "yes");
	// mMuscleOptimizationSolver->Options()->SetStringValue("mehrotra_algorithm", "yes");
	mMuscleOptimizationSolver->Options()->SetIntegerValue("print_level", 2);
	mMuscleOptimizationSolver->Options()->SetIntegerValue("max_iter", 100);
	mMuscleOptimizationSolver->Options()->SetNumericValue("tol", 1e-2);

	mMuscleOptimizationSolver->Initialize();
	mMuscleOptimizationSolver->OptimizeTNLP(mMuscleOptimization);
#endif
}
void
LowLevelController::
SetKp(double kp)
{
	int dof = mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();
	mKp = Eigen::VectorXd::Constant(dof,kp);
	mKv = Eigen::VectorXd::Constant(dof,2*sqrt(kp));
	mKjt = Eigen::VectorXd::Constant(dof,0.1);
	for(int i =0;i<6;i++)
	{
		mKp[dof-1-i] = 2.0*mKp[dof-1-i];	
		mKv[dof-1-i] = sqrt(2.0)*mKv[dof-1-i];
	}
	
}
Eigen::VectorXd
LowLevelController::
Solve(const Eigen::VectorXd& p_d,const Eigen::VectorXd& v_d)
{
	Eigen::VectorXd u;
	mTargetPositions = p_d;
	mTargetVelocities = v_d;
	
#ifdef USE_MUSCLE
	mIntegratedWorld->GetMusculoSkeletalSystem()->TransformAttachmentPoints();
	mIntegratedWorld->GetSoftWorld()->TimeStepping(false);
	u = ComputeActivationLevels(ComputePDForces());
#else
	u = ComputePDForces();
#endif
	return u;
}
Eigen::VectorXd
LowLevelController::
ComputePDForces()
{
	auto& skel =mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();

	Eigen::VectorXd pos_m = mTargetPositions;
	Eigen::VectorXd vel_m = mTargetVelocities;

	Eigen::VectorXd pos = skel->getPositions();
	Eigen::VectorXd vel = skel->getVelocities();

	Eigen::VectorXd pos_diff(pos.rows());

	pos_diff = skel->getPositionDifferences(pos_m,pos);
	for(int i = 0;i<pos_diff.rows();i++)
			pos_diff[i] = dart::math::wrapToPi(pos_diff[i]);
	double time_step = mIntegratedWorld->GetRigidWorld()->getTimeStep();
	Eigen::Vector6d F_b_r = skel->getBodyNode("HandR")->getConstraintImpulse()*(1.0/time_step);
	Eigen::Vector6d F_b_l = skel->getBodyNode("HandL")->getConstraintImpulse()*(1.0/time_step);
	Eigen::VectorXd jaco_t_r = skel->getInvMassMatrix()*(skel->getJacobian(skel->getBodyNode("HandR"),Eigen::Vector3d(0,0,0)).transpose()*F_b_r);
	Eigen::VectorXd jaco_t_l = skel->getInvMassMatrix()*(skel->getJacobian(skel->getBodyNode("HandL"),Eigen::Vector3d(0,0,0)).transpose()*F_b_l);   
	
	Eigen::VectorXd qdd_desired = pos_diff.cwiseProduct(mKp) + (mTargetVelocities - skel->getVelocities()).cwiseProduct(mKv);//- jaco_t_l.cwiseProduct(mKjt) - jaco_t_r.cwiseProduct(mKjt);
				
	return qdd_desired;
}
#ifdef USE_MUSCLE
Eigen::VectorXd
LowLevelController::
ComputeActivationLevels(const Eigen::VectorXd& qdd_desired)
{
	auto& skel =mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton();

	static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->Update(qdd_desired);

	mMuscleOptimizationSolver->ReOptimizeTNLP(mMuscleOptimization);	

	Eigen::VectorXd solution =  static_cast<MuscleOptimization*>(GetRawPtr(mMuscleOptimization))->GetSolution();
	Eigen::VectorXd qdd = solution.head(skel->getNumDofs());
	Eigen::VectorXd activation = solution.tail(mIntegratedWorld->GetMusculoSkeletalSystem()->GetNumMuscles());
	// std::cout<<(qdd-qdd_desired).norm()<<std::endl;
	
	// mIntegratedWorld->GetMusculoSkeletalSystem()->SetActivationLevels(activation);
	// mIntegratedWorld->GetMusculoSkeletalSystem()->TransformAttachmentPoints();
	// mIntegratedWorld->GetSoftWorld()->TimeStepping(false);
	// mIntegratedWorld->GetMusculoSkeletalSystem()->ApplyForcesToSkeletons(mIntegratedWorld->GetSoftWorld());
	// mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->computeForwardDynamics();
	// std::cout<<mIntegratedWorld->GetMusculoSkeletalSystem()->GetSkeleton()->getAccelerations().transpose()<<std::endl;
	// std::cout<<qdd.transpose()<<std::endl;
	// std::cout<<std::endl;
	return activation;
}
#endif

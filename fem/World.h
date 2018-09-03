#ifndef __FEM_WORLD_H__
#define __FEM_WORLD_H__
#include "Constraint/Constraint.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace FEM
{
enum IntegrationMethod
{
	NEWTON_METHOD,
	QUASI_STATIC,
	PROJECTIVE_DYNAMICS,
	PROJECTIVE_QUASI_STATIC
};
class World
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	void Initialize();
	void TimeStepping(bool integrate = true);

	void AddBody(const Eigen::VectorXd& x0, const std::vector<std::shared_ptr<Cst>>& constraints, double m = 1.0);
	void AddConstraint(const std::shared_ptr<Cst>& c);
	void RemoveConstraint(const std::shared_ptr<Cst>& c);

	void Computedxda(Eigen::VectorXd& dx_da,const Eigen::VectorXd& dg_da);

	IntegrationMethod GetIntegrationMethod(){return mIntegrationMethod;}
	int GetMaxIteration() {return mMaxIteration;}
	int GetNumVertices(){return mNumVertices;}
	double GetTime() {return mTime;}
	double GetTimeStep() {return mTimeStep;}
	const Eigen::VectorXd& GetPositions() {return mPositions;}
	std::vector<std::shared_ptr<Cst>>& GetConstraints() {return mConstraints;}

	void SetTime(double t) {mTime = t;};
	void SetPositions(const Eigen::VectorXd& X) {mPositions = X;};
	World(const World& other) = delete;
	World& operator=(const World& other) = delete;
	std::shared_ptr<World> Clone();
	static std::shared_ptr<World> Create(
			IntegrationMethod im = NEWTON_METHOD,
			double time_step = 1.0/120.0,
			int max_iteration = 100,
			const Eigen::Vector3d& gravity = Eigen::Vector3d(0,-9.81,0),
			double damping_coeff = 0.999);
private:
	World(	IntegrationMethod im = NEWTON_METHOD,
			double time_step = 1.0/120.0,
			int max_iteration = 100,
			const Eigen::Vector3d& gravity = Eigen::Vector3d(0,-9.81,0),
			double damping_coeff = 0.999);

	Eigen::VectorXd IntegrateNewtonMethod();
	Eigen::VectorXd IntegrateQuasiStatic();
	Eigen::VectorXd IntegrateProjectiveDynamics();
	Eigen::VectorXd IntegrateProjectiveQuasiStatic();

	void FactorizeLLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& llt_solver);
	void FactorizeLDLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldlt_solver);

	//For Newton method, Quasi-static
	double EvaluateEnergy(const Eigen::VectorXd& x);
	void EvaluateGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g);
	void EvaluateHessian(const Eigen::VectorXd& x,Eigen::SparseMatrix<double>& H);
	
	double EvaluateConstraintsEnergy(const Eigen::VectorXd& x);
	void EvaluateConstraintsGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g);
	void EvaluateConstraintsHessian(const Eigen::VectorXd& x,Eigen::SparseMatrix<double>& H);
	
	double ComputeStepSize(const Eigen::VectorXd& x, const Eigen::VectorXd& g,const Eigen::VectorXd& d);
	//For Projective Dynamics, Projective Quasi-static
	void EvaluateDVector(const Eigen::VectorXd& x,Eigen::VectorXd& d);
	void EvaluateJMatrix(Eigen::SparseMatrix<double>& J);
	void EvaluateLMatrix(Eigen::SparseMatrix<double>& L);

	void Precompute();
	//For detailing
	void InversionFree(Eigen::VectorXd& x);

	void IntegratePositionsAndVelocities(const Eigen::VectorXd& x_next);
private:
	bool mIsInitialized;
	int mNumVertices;
	int mConstraintDofs;
	int mMaxIteration;

	double mTimeStep,mTime;
	double mDampingCoefficient;
	Eigen::Vector3d mGravity;

	IntegrationMethod mIntegrationMethod;

	std::vector<double> mUnitMass;
	std::vector<std::shared_ptr<Cst>> mConstraints;

	Eigen::VectorXd mPositions,mVelocities;
	Eigen::VectorXd mExternalForces;

	Eigen::SparseMatrix<double> mMassMatrix,mInvMassMatrix,mIdentityMatrix;
	
	//For Projective Dynamics, Projective Quasi-static
	Eigen::VectorXd mq;
	Eigen::SparseMatrix<double> mJ,mL;

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> mSolver;
};
};

#endif